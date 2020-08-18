#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by Paul Minner <minner.paul@gmail.com>
#            Charles Reynolds <reynolds12@llnl.gov>             
# LLNL-CODE-704098
# All rights reserved.
# This file is part of Curvallis. 
# For details, see https://github.com/llnl/Curvallis.
# Please also Curvallis/LICENSE.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import math
import numpy
import bisect
import scipy
import scipy.interpolate
import scipy.misc
import scipy.integrate
from curvallis.version import version as VERSION_STRING


class TriLocalSmoother(object):
    """TriLocalSmoother creates and object for smoothing 1D data via local averaging.
    It tries to avoid shifting the graph by biasing closer points.
    This version strictly stops at teh edges, so the 2 edge points don't get averaged at all.
    """

    def __init__(self, in_numPts = 5, in_repeat = 1):
        self.maxNumPts = in_numPts
        self.repeat = in_repeat
        if(self.maxNumPts % 2 != 1):
            raise ValueError("TriLocalSmoother: numpoints must be odd, %d is not odd" % self.maxNumPts)

    def makeTriLocalCoefficents(self, xdata):
        midpoint = int(len(xdata) / 2)

        distance_xdata = list(map(lambda x: math.fabs(x - xdata[midpoint]), xdata))
        sumtot = sum(distance_xdata)
        avg_xdata = list(map(lambda x: 1 - (x / sumtot), distance_xdata))  # We want farther away points to have less importance, so avg - 1

        sumtot = sum(avg_xdata)
        coeff = map(lambda x: x / sumtot, avg_xdata)  # Normalize it again because the flip around threw it off

        return coeff

    def triLocalSmoothing(self, xdata, ydata, xmin, xmax):
        avgYdata = []

        ylen = len(ydata)


        xMinIndex = 0
        xMaxIndex = len(xdata)
        if (xmin > 0):
          xMinIndex = bisect.bisect_left(xdata, xmin)
        if(xmax >= 0):
          xMaxIndex = bisect.bisect_right(xdata, xmax)

        orig_lpts = int(self.maxNumPts / 2)  # Half the eval points left of center
        orig_rpts = self.maxNumPts - orig_lpts  # One at center, and the remainder at right (right includes center pt).

        for ii in range(xMinIndex, xMaxIndex):
            lpts = orig_lpts
            rpts = orig_rpts
            if (ii < lpts):
                lpts = ii
                rpts = lpts + 1  # rpts includes the center point, so it's always at least 1
            elif (ii + (rpts - 1) > ylen - 1):
                rpts = ylen - ii
                lpts = rpts - 1

            if (lpts == 0 or rpts == 1):
                xpts = [xdata[ii]]
                ypts = [ydata[ii]]
                coeff = [1]
            else:
                xpts = xdata[ii - lpts: ii + rpts]
                ypts = ydata[ii - lpts: ii + rpts]
                coeff = list(self.makeTriLocalCoefficents(xpts))

            sumtot = 0

            for xx in range(0, len(xpts)):
                sumtot += coeff[xx] * ypts[xx]

            avgYdata.append(sumtot)

        ydata[xMinIndex:xMaxIndex] = avgYdata
        return ydata

    def applySmooth(self, xdata, ydata, xmin, xmax):

      yTemp = ydata
      xTemp = xdata

      for ii in range(0,self.repeat):
          yTemp = self.triLocalSmoothing(xTemp, yTemp, xmin, xmax)


      return (xTemp, yTemp)


class AcuteAngleRepair(object):
    """AcuteAndleRepair creates and object for for repairing outliers.
       Some functions are smooth except for one point that's out of line.
       AcuteAndleRepair finds that point, deletes it, and replaces it with
       interpolation.
       If multiple acute angles occur in sequence, an optional averaging
       algorithm can be used to repair.
    """

    def __init__(self, in_minAngle = 90, in_repeat = 1):
        self.minAngle = math.radians(in_minAngle) #The minimum angle to allow.  Comes in as degrees, internally radians are used
        self.repeat = in_repeat

    #Thanks pythagorus
    def lineLength(self, x1, y1, x2, y2):
        xdiff = x1 -x2
        ydiff = y1 - y2
        tmp = (xdiff * xdiff) + (ydiff * ydiff)
        tmp = math.fabs(tmp)
        return math.sqrt(tmp)

    #calculates angle on point 2:
    def calcAngle(self, x1, y1, x2, y2, x3, y3):
        sumy = y1 + y2 + y3  #Normalize the angles.  This avoids huge mismatches between x & y
        y1 = y1 / sumy
        y2 = y2 / sumy
        y3 = y3 / sumy
        sumx = x1 + x2 + x3
        x1 = x1 / sumx
        x2 = x2 / sumx
        x3 = x3 / sumx

        #print "1:", x1, y1
        #print "2:", x2, y2
        #print "3:", x3, y3
        z = self.lineLength(x2, y2, x3, y3) #x and z are the adjacent lines
        x = self.lineLength(x2, y2, x1, y1)
        y = self.lineLength(x1, y1, x3, y3) #y is the opposite line
        #print x,y,z
        top = x*x + z*z - y*y
        bottom = 2 * x * z
        whole = top/bottom
        
        #print "top %g bottom %g whole %g" % ( top, bottom, whole)
        if(whole > 1 or whole < -1):
            return math.pi

        angle = math.acos( whole)
        return angle

    #Just use two cartesian points to find the y value of an x between them
    def linInterp(self, x0, y0, x1, y1, newx):
        #print "linInterp:", x0, y0, x1, y1
        newy = y0 + ((y1 - y0) * ((newx - x0) / (x1 - x0))) #From wikipedia
        return newy

    #avgPoints replaces a block of bad points by averaging them, both in X and Y.  So you
    #End up with n-1 points at the end.  
    def avgPts(self, xdata, ydata, startIndex, endIndex):
        resXdata = []
        resYdata = []
        for ii in range(startIndex, endIndex):
            resXdata.append((xdata[ii] + xdata[ii+1])/2)
            resYdata.append((ydata[ii] + ydata[ii+1])/2)

        return (resXdata, resYdata)

    def acuteAngleRepair(self, xdata, ydata, xmin, xmax):
        avgYdata = []

        ylen = len(ydata)


        xMinIndex = 0
        xMaxIndex = len(xdata)
        if (xmin > 0):
          xMinIndex = bisect.bisect_left(xdata, xmin)
        if(xmax >= 0):
          xMaxIndex = bisect.bisect_right(xdata, xmax)

        badAngles = []  #Bad angles is an array of bool.  Each point gets a 0 (good) or 1 (bad)  Bad points must be replaced.
        badAngles.append(0) #Can't check the endpoint.
        #Step through all the data, you can't check the endpoints because they don't form a triangle.
        for ii in range(xMinIndex+1, xMaxIndex-1):
            angle = self.calcAngle(xdata[ii-1], ydata[ii-1], xdata[ii], ydata[ii], xdata[ii+1], ydata[ii+1])
            #print "x: %g angle: %g" % ( xdata[ii], angle)
            if(angle < self.minAngle):
                badAngles.append(1)
            else:
                badAngles.append(0)

        #for ii in range(0, len(badAngles)):
            #print xdata[xMinIndex + ii], badAngles[ii]

        resXdata = []
        resYdata = []
        ii = 0
        while ii < len(badAngles):
            if(badAngles[ii] == 0):
                resXdata.append(xdata[xMinIndex+ii])
                resYdata.append(ydata[xMinIndex+ii])
            else:
#                if(badAngles[ii+1] == 0):  #If this is a one off bad angle, just linearly interpolate
                resXdata.append(xdata[xMinIndex+ii])
                resYdata.append(self.linInterp( xdata[xMinIndex+ii-1], ydata[xMinIndex+ii-1],
                                                xdata[xMinIndex+ii+1], ydata[xMinIndex+ii+1],
                                                xdata[xMinIndex+ii]))
                """
                else:  #Otherwise we have a block of bad points, send them off to be dealt with
                    startIndex = ii
                    endIndex = ii+1;
                    while(badAngles[endIndex+1] == 1):  #Find the end of the block of bad points
                        endIndex += 1
#                    print "Last bad index, should be 1", badAngles[endIndex], startIndex, endIndex
                    (newXPts, newYPts) = self.avgPts(xdata, ydata, startIndex + xMinIndex, endIndex + xMinIndex)
                    resXdata.extend(newXPts)
                    resYdata.extend(newYPts)
#                    print "ii", ii
                    ii = endIndex+1  #Jump past the points we just filled in
#                    print "after ii", ii
                    continue
                    """
            ii += 1

        xdata[xMinIndex:xMaxIndex-1] = resXdata
        ydata[xMinIndex:xMaxIndex-1] = resYdata
        return (xdata, ydata)


    def applySmooth(self, xdata, ydata, xmin, xmax):

      yTemp = ydata
      xTemp = xdata

      for ii in range(0,self.repeat):
          (xTemp, yTemp) = self.acuteAngleRepair(xTemp, yTemp, xmin, xmax)


      return (xTemp, yTemp)



class IntegralSmoother(object):
    """ Integral Smoother is a smoothing algorithm that works by taking the derivative of the function
    then taking it's integral to get back to the original function.  The integral smooths things out a bit
    so you end up at a smoother function (usually).

    You can also pass in another smoother as "derivRepair" this will be used to smooth the derivative.
    """


    def __init__(self, in_xMatchPoint = -1, in_derivRepair = None, in_interpString = 'cubic'):
        self.xMatchPoint = in_xMatchPoint;
        self.derivRepair = in_derivRepair
        self.interpString = in_interpString

    #Since we can't get derivatives right at the end point, this gives a small, but valid shift
    #That we should be able to get the endpoint at
    def endNumber(self, aa,bb, diff):
        tmp = (bb + aa )/ 2
        tmp -= aa
#        tmp *= diff * 10
        return tmp

    #xmin and xmax are the range to do smoothing on.
    def applySmooth(self, xdata, ydata, xmin, xmax):
        diff = 1e-8

        xMinIndex = 0
        xMaxIndex = len(xdata)
        if (xmin > 0):
          xMinIndex = bisect.bisect_left(xdata, xmin)
        if(xmax >= 0):
          xMaxIndex = bisect.bisect_right(xdata, xmax)

        interpolator = scipy.interpolate.interp1d(xdata, ydata, kind=self.interpString, bounds_error=False)

        if(self.xMatchPoint >= 0 and xmin > 0 and self.xMatchPoint < xmin):
             raise ValueError("xMatchPoint %g must be inside xmin, xmax range [%g, %g]" % (self.xMatchPoint, xmin, xmax))
        if(self.xMatchPoint >= 0 and xmax > 0 and self.xMatchPoint > xmax):
             raise ValueError("xMatchPoint %g must be inside xmin, xmax range [%g, %g]" % (self.xMatchPoint, xmin, xmax))
        localXMatch = self.xMatchPoint
        if(localXMatch < 0):  #Default to xmin
             localXMatch = xdata[xMinIndex]


        yMatchPoint = interpolator(localXMatch)

        dydxData = []
        locXData = xdata[xMinIndex : xMaxIndex]
        if(xMinIndex == 0): #Python can't take derivitives on the boundry point, so use this case
             datapt1 = self.endNumber(xdata[0], xdata[1], diff) + xdata[0]
             dydxData.append(scipy.misc.derivative(interpolator,
                                           datapt1,
                                           dx=(diff * xdata[0]) ))
        else:
             dydxData.append(scipy.misc.derivative(interpolator,xdata[xMinIndex], dx=(1e-8 * xdata[xMinIndex]) ))


        for xx in range(xMinIndex+1, xMaxIndex -1):  #End points are done seperately
             dydxData.append(scipy.misc.derivative(interpolator,xdata[xx], dx=(1e-8 * xdata[xx]) ))

        if(xMaxIndex == len(xdata)):
            datapt2 = xdata[-1] - self.endNumber(xdata[-2], xdata[-1], diff)
            dydxData.append(scipy.misc.derivative(interpolator,
                                                           datapt2,
                                                           dx=(diff * xdata[-1]) ))
        else:
            dydxData.append(scipy.misc.derivative(interpolator,xdata[xMaxIndex], dx=(1e-8 * xdata[xMaxIndex]) ))

        if(self.derivRepair):
            dydxData = self.derivRepair.applySmooth(locXData, dydxData, -1, -1)

        dydxInterp = scipy.interpolate.interp1d(locXData, dydxData, kind=self.interpString, bounds_error=False)

        matchResult = scipy.integrate.quad(dydxInterp, locXData[0], localXMatch)
        constYAddr = yMatchPoint - matchResult[0]; #constant of integration


        outYData = []
        for xx in range(0, len(dydxData)):
             result = scipy.integrate.quad(dydxInterp, locXData[0], locXData[xx])
             yval = result[0] + constYAddr;
             outYData.append(yval)

        ydata[xMinIndex:xMaxIndex] = outYData
        return (xdata, ydata)



class BSplineSmoother(object):
    """BSplineSmoother creates an object for smoothing 1D data via B-Spline smoothing"""
    def __init__(self, detail=0):
        # The detail of the curve. Higher values give greater accuracy but slower operation. The reverse is also true.
        # It is recommended to keep detail
        self.detail = detail
    def applySmooth(self, datax, datay, xmin, xmax):
        # l= number of input data points
        l = len(datax)
        if(l < 4):
            print("Error: Must have at least 4 points selected")
            return (datax,datay)
        t = numpy.linspace(0,1,l-2,endpoint=True)
        t = numpy.append([0,0,0],t)
        t = numpy.append(t, [1,1,1])
        # t= knots, c= coefficients, k= degree of spline
        tck = [t,[datax,datay],3]
        # The linspace that will be used to compute the spline
        u3 = numpy.linspace(0,1,l*self.detail,endpoint=True)
        out = scipy.interpolate.splev(u3,tck)
        # Moves the input points to their place on the curve.
        newy = []
        for i in range(l):
            x = datax[i]
            for i in range(len(out[0])):
                if(out[0][i] == x):
                    newy.append(out[1][i])
                    break
                elif(i < len(out[0])-1 and out[0][i] > x):
                    # y = m*x+b
                    x1, y1 = out[0][i], out[1][i]
                    x2, y2 = out[0][i+1], out[1][i+1]
                    m = (y2-y1)/(x2-x1)
                    b = -1*((m*x1)-y1)
                    y = (m*x)+b
                    newy.append(y)
                    break
                elif(out[0][i] > x):
                    newy.append(out[1][i])
                    break
        return (datax,newy)

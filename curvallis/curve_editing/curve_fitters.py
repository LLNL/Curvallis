# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by Paul Minner <minner.paul@gmail.com>
#            Charles Reynolds <reynolds12@llnl.gov>             
# LLNL-CODE-704098
# All rights reserved.
# This file is part of Curvallis. 
# For details, see https://github.com/llnl/Curvallis.
# Please also Curvallis/LICENSE.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""This module holds the various curve-fitting methods used by curve_editor.py
"""

# TODO: energy functions for AP1, AP2, AP3
# TODO: refactor convenience functions
# TODO: refactor fit functions
# TODO: Add generic initial parameters for johnson, kumari_dass
# TODO: fit curves for energy vs. density datasets


from __future__ import print_function

from abc import ABCMeta, abstractmethod
import numpy as np
import scipy.optimize as opt
import scipy.special as spec
import scipy.interpolate as interp
from pylab import polyfit
import math
import bisect
import re
from curvallis.window import fitter_info_window
from curvallis.window import update_fitter_info_window
from curvallis.version import version as VERSION_STRING

_min_polynomial_degree = 1
_max_polynomial_degree = 12
rho = '\u03c1'
naught = '\u2080'

def define_args(parser):
    fitter_args = parser.add_argument_group(
        title='Curve fitter arguments',
        description='These arguments select and parameterize the curve '
                    'fitter formulas.')

    fitter_args.add_argument(
        '--fit_type',
        choices=factory.get_sorted_fit_names(),
        metavar='<fit>',
        nargs='+',
        help=str(factory.get_sorted_fit_names()) +
             ' Select the type of fit to apply to each region in sequential order[default: %(default)s]')
    fitter_args.add_argument(
        '--refine_fit',
        choices=factory.get_sorted_refine_fit_names(),
        metavar='<fit>',
        nargs='+',
        help=str(factory.get_sorted_refine_fit_names()) +
             ' Select the type of fit to apply to the error of each fit type[default: %(default)s]')
    fitter_args.add_argument(
        '--delta_p_guess', dest='k0',
        type=float, metavar='<float>',
        help='Initial guess for bulk modulus [default: %(default)s]')
    fitter_args.add_argument(
        '--k0_guess', dest='k0',
        type=float, metavar='<float>',
        help='Initial guess for bulk modulus [default: %(default)s]')
    fitter_args.add_argument(
        '--k0_prime_guess', dest='k0_prime',
        type=float, metavar='<float>',
        help='Initial guess for derivative of bulk modulus [default: %(default)s]')
    fitter_args.add_argument(
        '--k0_prime_prime_guess', dest='k0_prime_prime',
        type=float, metavar='<float>',
        help='Initial guess for second derivative of bulk modulus [default: %(default)s]')
    fitter_args.add_argument(
        '--lam_guess', dest='lam',
        type=float, metavar='<float>',
        help='Initial guess for lambda [default: %(default)s]')
    fitter_args.add_argument(
        '--rho0_guess', dest='rho0',
        type=float, metavar='<positive float>',
        help='Initial guess for equilibrium density [default: %(default)s]')
    fitter_args.add_argument(
        '--scale_derivative_by', dest='derivative_scale',
        type=float, metavar='<float>',
        help='Scaling factor for derivative plot [default: %(default)s]')
    fitter_args.add_argument(
        '--scale_second_derivative_by', dest='second_derivative_scale',
        type=float, metavar='<float>',
        help='Scaling factor for second derivative plot [default: %(default)s]')
    fitter_args.add_argument(
        '--scale_integral_by', dest='integral_scale',
        type=float, metavar='<float>',
        help='Scaling factor for integral plot [default: %(default)s]')
    fitter_args.add_argument(
        '--xref', dest='x_integral_ref',
        type=float, metavar='<float>',
        help='Reference x value to determine integration constant [default: %(default)s]')
    fitter_args.add_argument(
        '--y_axis',
        choices=['P', 'E'], metavar='<P or E>',
        help='Choose Pressure or Energy for the integral y axis [default: %(default)s]')
    fitter_args.add_argument(
        '--yref', dest='y_integral_ref',
        type=float, metavar='<float>',
        help='Reference y value to determine integration constant [default: %(default)s]')
    fitter_args.add_argument(
        '--e0_guess', dest='e0',
        type=float, metavar='<float>',
        help='Initial guess for E0 [default: %(default)s]')
    fitter_args.add_argument(
        '--theta0_guess', dest='theta0',
        type=float, metavar='<float>',
        help='Initial guess for theta0 [default: %(default)s]')
    fitter_args.add_argument(
        '--gamma0_guess', dest='gamma0',
        type=float, metavar='<float>',
        help='Initial guess for gamma0 [default: %(default)s]')
    # fitter_args.add_argument(
    #     '--z_guess', dest='z',
    #     type=float, metavar='<float>',
    #     help='Initial guess for z parameter [default: %(default)s]')

    parser.set_defaults(
        fit_type=['poly5'],
        refine_fit=['none'],
        derivative_scale=1000,
        second_derivative_scale=1000,
        integral_scale=.001,
        polynomial_order=5,
        x_integral_ref=0,
        y_axis='P',
        y_integral_ref=0,

        a=1.0,
        b=1.0,
        c=1.0,
        d=1.0,
        e=1.0,
        f=1.0,
        delta_p=1.0,
        k0=None,  # Calculated in guess_coefficients
        k0_prime=None,  # May be calcualted in guess_coefficients
        k0_prime_prime=1.0,
        e0=None,  # Energy Curve, calcualted in guess_coefficients
        kneg1=1.0e11,
        kneg0=1.0e11,
        k1=1.0e11,
        k2=1.0e11,
        k3=1.0e11,
        k4=1.0e11,
        k5=1.0e11,
        lam=1.0,
        rho0=None,  # Same
        q=None,  # For Theta and Gamma fitters
        theta0=None,  # For Theta fitters
        gamma0=None,  # For gamma fitters
        z=1.0,
        c1=1.0,
        c2=1.0,
        c3=1.0,
        c4=1.0,
        cn=[1.0] * 12,  # 12 is maximum power of eseries
    )

# ------------------------------------------------------------------------------
# Supporting functions

def _rno_norm(x, rho0):
    return np.asarray(x) / rho0


def _eta(x, rho0):
    return spec.cbrt(_rno_norm(x, rho0))


def _eta2(x, rho0):
    return np.square(_eta(x, rho0))


def _eta3(x, rho0):
    return _rno_norm(x, rho0)


def _eta4(x, rho0):
    return np.power(_eta(x, rho0), 4.0)


def _eta5(x, rho0):
    return np.power(_eta(x, rho0), 5.0)


def _eta6(x, rho0):
    return np.square(_rno_norm(x, rho0))


def _eta7(x, rho0):
    return np.power(_eta(x, rho0), 7.0)


def _eta8(x, rho0):
    return np.power(_eta(x, rho0), 8.0)


def _mu(x, rho0):
    return np.asarray(x) / rho0 - 1.0


def _mu2(x, rho0):
    return np.square(_mu(x, rho0))


def _mu3(x, rho0):
    # noinspection PyTypeChecker
    return np.power(_mu(x, rho0), 3.0)


def _mu4(x, rho0):
    # noinspection PyTypeChecker
    return np.power(_mu(x, rho0), 4.0)


def _pFG0(rho0, z):
    #     return (.02337) * np.power(z * np.asarray(x), (5.0 / 3.0))
    return .02337 * (rho0 * z) ** (5.0 / 3.0)


def _c1(k0, rho0, z):
    return -np.log(3.0 * k0 / _pFG0(rho0, z))


def _c2(k0, k0_prime, rho0, z):
    return 1.5 * (k0_prime - 3) - _c1(k0, rho0, z)


def _c0x(x, k0, rho0, z):
    return _c1(k0, rho0, z) * _eta(x, rho0)


class Factory(object):
    """ Collects all the classes and can create an object of any of them.
    """

    def __init__(self):
        self._classes = dict()
        self._poly_classes = dict()  # classes that expect users to specify polynomial order in the fit_type name

    def register(self, name, classs):
        assert name not in self._classes
        self._classes[name] = classs
        if issubclass(classs, PolyBase):
            self._poly_classes[classs.name_prefix] = classs

    def get_sorted_fit_names(self):
        names = list(self._classes)

        # Only allows polynomials up to 12 to avoid crashing
        for pre in self._poly_classes.keys():
            for i in range(_min_polynomial_degree, _max_polynomial_degree + 1):
                names.append(pre + str(i))
            names.remove(pre)

        names.sort()
        return names

    def get_sorted_refine_fit_names(self):
        names = []
        for i in range(0, 12):
            names.append(Poly_Original.name_prefix + str(i + 1))
        for i in range(3, 12):
            names.append('eseries' + str(i + 1))
        names.append('eseries')
        names.append('none')
        names.append('highp')
        return names

    def make_object_of_class(self, name, args):
        r = re.compile('({0})\d+'.format('|'.join(self._poly_classes.keys())))
        polymatch = r.match(name)
        if polymatch:
            return self._classes[polymatch.group(1)](args, name)
        elif name[:7] == 'eseries':
            return self._classes[name[:7]](args, name)
        else:
            return self._classes[name](args, name)


factory = Factory()

class MetaPoly(ABCMeta):
    @property
    def name_prefix(cls):
        return cls._name_prefix


class PolyBase(object, metaclass=MetaPoly):
    """
    Base class for fitters that specify polynomial order in the fit_type name
    Must be defined before any call to factory.register()
    """

    @property
    def name_prefix(self):
        return type(self).name_prefix


class Base_Fit_Class(object):
    """ Call fit to (re-)calculate the parameters for the function, which are
    stored in the object.  Call f, derivative, and integral to get the y value
    for a point to be plotted.
    """
    __metaclass__ = ABCMeta

    def __init__(self, args, name):
        # The list of parameters parameters for _f.  This is a list so that it
        # can be passed into opt.curve_fit:
        super(Base_Fit_Class, self).__init__()
        self.derivative_scale = args.derivative_scale
        self.second_derivative_scale = args.second_derivative_scale
        self.integral_scale = args.integral_scale
        self.maxfev = 50000
        self.use_energy_integral = args.y_axis == 'E'
        self.x_integral_ref = args.x_integral_ref
        self.y_integral_ref = args.y_integral_ref
        self.name = name
        self._is_first_fit = True  # Is first regression fit

    @staticmethod
    @abstractmethod
    def _f(x, *coeffs):
        """ This is the internal, parametrized function, called by fit() and
        f()
        """

    @abstractmethod
    def _set_coefficients(self, coeffs):
        """ Sets this fitter's coefficients from a list.

        :param coeffs: list
        """

    @abstractmethod
    def _get_coefficients(self):
        """ Gets this fitters coefficients as a list, for use by fit().

        :return: list
        """

    @abstractmethod
    def _print_coefficients(self):
        """ Prints coefficients used by specific equation
        """

    def bounds(self):
        """ Gives bounds, if any, to how parameters for the
            'fit_to_points()' function should be optimized. Default
            is no bounds, though if values are given by the user that
            shouldn't be changed, their bounds should be set here
            to prevent change
        """
        return (-np.inf, np.inf)

    def fit_to_points(self, points):
        """ Derive the coefficients for this fit function that best fit the
        data.  Call this whenever the data changes.  Args to _f are provided via
        the p0 parameter.
        """
        x_values = [point[0] for point in points]
        y_values = [point[1] for point in points]
        new_coefficients, unused = opt.curve_fit(
            f=self._f,
            xdata=x_values,
            ydata=y_values,
            p0=self._get_coefficients(),
            bounds=self.bounds(),
            maxfev=self.maxfev)
        self._set_coefficients(new_coefficients)
        #        print("curve_fit parameters out = ", new_coefficients)
        self._print_coefficients()
        print("")

    @staticmethod
    @abstractmethod
    def guess_coefficients(self, points):
        """ Does the first guess of coefficients for this fitter
        """

    def func(self, x):
        """ Return f(x) using the coefficients for f derived by fit().  Call this
        to plot each point on the curve.
        """
        return self._f(x, *self._get_coefficients())

    @abstractmethod
    def _derivative(self, x):
        """ Return the derivative of f(x).
        """

    def derivative(self, x):
        """ Return the scaled derivative of f(x). Call this to plot each point
            on the curve.
        """
        return self.derivative_scale * self._derivative(x)

    @abstractmethod
    def _second_derivative(self, x):
        """ Return the second derivative of f(x).
        """

    def second_derivative(self, x):
        """ Return the scaled second derivative of f(x). Call this to plot each point
            on the curve.
        """
        return self.second_derivative_scale * self._second_derivative(x)

    @abstractmethod
    def _energy_integral(self, x):
        """ Return the energy integral of f(x).
        """

    @abstractmethod
    def _pressure_integral(self, x):
        """ Return the pressure integral of f(x).
        """

    def integral(self, x):
        """ Return the integral of f(x).  Call this to plot each point on the
        curve.
        """
        if self.use_energy_integral:
            integral_func = self._energy_integral
        else:
            integral_func = self._pressure_integral
        return self.integral_scale * \
               (integral_func(x) +
                (self.y_integral_ref - integral_func(self.x_integral_ref)))


class Pressure_Fit_Class(Base_Fit_Class):
    """ An abstract superclass for all the fitters for pressure curves.
    It defines some shared functionality.
    """

    def __init__(self, args, name):
        super(Pressure_Fit_Class, self).__init__(args, name)

    def guess_coefficients(self, points):

        # Check if selected plotter uses coefficients
        if ((not hasattr(self, 'k0')) or
                (not hasattr(self, 'k0_prime')) or
                (not hasattr(self, 'rho0'))):

            # If not, set defaults for required coefficients
            if (not hasattr(self, 'k0')):
                self.k0 = 1.0e11
            if (not hasattr(self, 'k0_prime')):
                self.k0_prime = 20.0
            if (not hasattr(self, 'rho0')):
                self.rho0 = 5.0

            return

        # Find point closest to zero
        zero_point = points[0]
        zero_point_index = 0
        for idx in range(0, len(points)):
            if (abs(points[idx][1]) < abs(zero_point[1])):
                zero_point = points[idx]
                zero_point_index = idx

        parabola_points = []
        # Create list of points around zero point
        if zero_point_index > 25 and len(points) > zero_point_index + 25:
            startidx = zero_point_index - 25
            for i in range(startidx, startidx + 50):
                parabola_points.append(points[i])
        elif zero_point_index > 0 and len(points) > zero_point_index + 25:  # Or just take what points we have
            for i in range(0, zero_point_index * 2):
                parabola_points.append(points[i])
        #        else:  zero_point_index is 0, nothing we can do

        # Use the points around the zero point to fit a parabola to guess
        # coefficients. This only makes sense for a pressure curve
        if len(parabola_points) > 0:
            xcoords = [point[0] for point in parabola_points]
            ycoords = [point[1] for point in parabola_points]
            a, b, c = polyfit(xcoords, ycoords, 2)
            rho0 = zero_point[0]
            k0 = ((2.0 * a * rho0) + b) * rho0
            k0_prime = 1 + ((2 * a * rho0 * rho0) / k0)

            if self.rho0 is None:
                self.rho0 = rho0
            if self.k0 is None:
                self.k0 = k0
            if self.k0_prime is None:
                self.k0_prime = k0_prime

        else:
            print("Unable to guess coefficients, using defaults.")

            if self.rho0 is None:
                self.rho0 = 5.0
            if self.k0 is None:
                self.k0 = 1.0e11
            if self.k0_prime is None:
                self.k0_prime = 20.0

        print("Coefficient Guesses: ", self.k0,
              self.k0_prime, self.rho0)

        self._is_first_fit = False


class Energy_Fit_Class(Base_Fit_Class):
    """ An abstract superclass for all the fitters for pressure curves.
    It defines some shared functionality.
    """

    def __init__(self, args, name):
        super(Energy_Fit_Class, self).__init__(args, name)

    def guess_coefficients(self, points):

        if (len(points) < 10):
            print("The data does not have many points in it.  The fit will probably be poor.");

        if (self.rho0 == None):
            # Find the minimum point, that's rho0. ( unless rho0 isn't in the region)
            min_point = points[0]
            min_point_index = 0
            for idx in range(0, len(points)):
                if (points[idx][1] < min_point[1]):
                    min_point = points[idx]
                    min_point_index = idx
        else:  # a rho0 guess exists, so try to find it in the data
            xcoords = [point[0] for point in points]
            ycoords = [point[1] for point in points]
            min_point_index = bisect.bisect(xcoords,
                                            self.rho0) - 1  # bisect seems to find the index after the target for some reason.
            min_point = points[min_point_index]

        parabola_points = []
        # Create list of points around min point to make a parabola if possible
        if min_point_index <= 0 or min_point_index >= len(points) - 1:
            # rho0 is at or off the end of the data range, so we can't make a parabola
            if (self.rho0 == None):
                print("rho0 doesn't appear to be in the data, please provide a guess for rho0.")
                self.rho0 = 5.0
            else:
                print(
                    "Your guess for rho0 doesn't appear to be in the data range.\n Cannot use it to guess other coefficients. The fit may be poor");

        else:  # Get points as many points as possible around the min_point, up to 10.
            self.rho0 = min_point[0]
            start_idx = min_point_index - 10
            if (start_idx < 0):
                start_idx = 0
            end_idx = min_point_index + 10
            if (end_idx > len(points)):
                end_idx = len_points
            parabola_points.extend(points[start_idx:end_idx])

        # Use the points around rho0 to fit a parabola to guess
        # coefficients. This only makes sense for an energy curve
        if len(parabola_points) > 0:
            xcoords_tup = [point[0] for point in parabola_points]
            ycoords = [point[1] for point in parabola_points]
            xcoords = list(map(lambda rho: rho - self.rho0,
                               xcoords_tup))  # Shift the rho axis so rho0 is 0 (otherwise the parabolo fit doesn't work)
            a, b, c = polyfit(xcoords, ycoords, 2)

            if self.k0 is None:
                self.k0 = a  # a is the coefficent of x^2
            if self.k0_prime is None:
                self.k0_prime = 5.0
            if self.e0 is None:
                self.e0 = c  # c should be the parabola minimum point on the y axis

        else:
            print("Unable to guess coefficients, using defaults.")

            if self.rho0 is None:
                self.rho0 = 5.0
            if self.k0 is None:
                self.k0 = 1.0e11
            if self.k0_prime is None:
                self.k0_prime = 5.0
            if self.e0 is None:
                self.e0 = 1.0

        print("Coefficient Guesses: ", self.k0,
              self.k0_prime, self.rho0)

        self._is_first_fit = False


# ------------------------------------------------------------------------------
# EOS models tested against MEOS Equations
# These models have all been tested against the MEOS data to make sure they
# produce the same curve.

class Birch_Murnaghan3(Pressure_Fit_Class):
    """ Fit function for third-order Birch-Murnaghan EOS (BE2?)
    Jeanloz, Raymond. "Finite-Strain Equation of State for High-Pressure Phases"
        Geophysical Research Letters 8(12) 1219-1222, December 1981
    J.-P. Poirier and A. Tarantola, "A logarithmic equation of state"
        Physics of the Earth and Planetary Interiors 109 (1998) 1-8
    """

    def __init__(self, args, name):
        super(Birch_Murnaghan3, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        xi = (.75 * (k0_prime - 4.0))
        return (1.5 * k0) * \
               (_eta7(x, rho0) - _eta5(x, rho0)) * \
               (1.0 + xi * (_eta2(x, rho0) - 1.0))

    def _derivative(self, x):
        xi = .75 * (4.0 - self.k0_prime)
        return ((self.k0 / (2.0 * self.rho0)) *
                ((7.0 - 14.0 * xi) * _eta4(x, self.rho0) +
                 5.0 * (xi - 1.0) +
                 9.0 * xi * _eta6(x, self.rho0)))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        xi = .75 * (4.0 - self.k0_prime)
        return (3.0 * self.k0 * self.rho0 / 80.0) * \
               _eta8(x, self.rho0) * \
               (2.0 * _eta2(x, self.rho0) *
                ((5.0 * _eta2(x, self.rho0) - 12.0) * xi + 6.0) +
                15.0 * (xi - 1.0))

    def _pressure_integral(self, x):
        xi = .75 * (4.0 - self.k0_prime)
        return (3.0 * self.k0 * _eta2(x, self.rho0) / (8.0 * self.rho0)) * \
               (6.0 * (xi - 1.0) +
                _eta2(x, self.rho0) *
                (3.0 + 2.0 * (_eta2(x, self.rho0) - 3.0) * xi))


factory.register('birch3', Birch_Murnaghan3)


class Birch_Murnaghan4(Pressure_Fit_Class):
    """ Fit function for fourth-order Birch-Murnaghan EOS
        Jeanloz, Raymond. "Finite-Strain Equation of State for High-Pressure Phases"
            Geophysical Research Letters 8(12) 1219-1222, December 1981
        J.-P. Poirier and A. Tarantola, "A logarithmic equation of state"
            Physics of the Earth and Planetary Interiors 109 (1998) 1-8
    """

    def __init__(self, args, name):
        super(Birch_Murnaghan4, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.k0_prime_prime = args.k0_prime_prime
        self.rho0 = args.rho0
        self.BMurn3 = Birch_Murnaghan3(args, name)

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, k0_prime_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.k0_prime_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("Bpp = {};".format(self.k0_prime_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("Bpp = {};\n".format(self.k0_prime_prime)) + ("rho0 = {};".format(self.rho0)))

    # Birch_Murnaghan4 doesn't seem to converge well, use BMurn3 as inital guess.
    def guess_coefficients(self, points):
        self.BMurn3.guess_coefficients(points)
        self.BMurn3.fit_to_points(points)
        (self.k0, self.k0_prime, self.rho0) = self.BMurn3._get_coefficients()

    @staticmethod
    def _f(rho, *coeffs):
        (k0, k0_prime, k0_prime_prime, rho0) = coeffs
        strain = (1.0 / 2.0) * pow((rho / rho0), (2.0 / 3.0)) - 1.0

        term1 = (3.0 / 2.0) * (k0 * k0_prime_prime + k0_prime * (k0_prime - 7) + (143.0 / 9.0)) * pow(strain, 2.0)
        term2 = 1 + (3.0 / 2.0) * (k0_prime - 4) * strain
        bracket_term = term1 + term2
        outer_term = 3 * k0 * strain * pow(1 + 2 * strain, (5.0 / 2.0))
        return outer_term * bracket_term

    def _derivative(self, x):
        xi = .75 * (4.0 - self.k0_prime)
        zeta = (3.0 / 8.0) * \
               (self.k0 * self.k0_prime_prime +
                self.k0_prime * (self.k0_prime - 7.0) +
                (143.0 / 9.0))
        return ((self.k0 / (2.0 * self.rho0)) *
                (9.0 * (xi - 3 * zeta) * _eta6(x, self.rho0) +
                 5.0 * (xi - zeta - 1.0) +
                 11.0 * zeta * _eta8(x, self.rho0) +
                 7.0 * (1.0 - 2.0 * xi + 3.0 * zeta) * _eta4(x, self.rho0)))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        xi = .75 * (4.0 - self.k0_prime)
        zeta = (3.0 / 8.0) * \
               (self.k0 * self.k0_prime_prime +
                self.k0_prime * (self.k0_prime - 7.0) +
                (143.0 / 9.0))
        return (3.0 * self.k0 * _eta8(x, self.rho0) * self.rho0 / 560.0) * \
               (70.0 * _eta4(x, self.rho0) * (xi - 3.0 * zeta) +
                84.0 * _eta2(x, self.rho0) * (3.0 * zeta - 2.0 * xi + 1.0) +
                60.0 * _eta6(x, self.rho0) * zeta +
                105.0 * (xi - zeta - 1.0))

    def _pressure_integral(self, x):
        xi = .75 * (4.0 - self.k0_prime)
        zeta = (3.0 / 8.0) * \
               (self.k0 * self.k0_prime_prime +
                self.k0_prime * (self.k0_prime - 7.0) +
                (143.0 / 9.0))
        return 3.0 * self.k0 * _eta2(x, self.rho0) * \
               (4.0 * _eta4(x, self.rho0) * (xi - 3.0 * zeta) +
                12.0 * (-1.0 + xi - zeta) +
                3.0 * _eta6(x, self.rho0) * zeta +
                _eta2(x, self.rho0) * (6.0 - 12.0 * xi + 18.0 * zeta)) / \
               (16.0 * self.rho0)


factory.register('birch4', Birch_Murnaghan4)


class Vinet(Pressure_Fit_Class):
    """ Fit function for third-order Vinet EOS
        R. Jeanloz, "Universal equation of state"
            Physical Review B 38(1) 805
    """

    def __init__(self, args, name):
        super(Vinet, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        veta = 1.5 * (k0_prime - 1.0)
        return 3.0 * \
               k0 * \
               ((1.0 - _eta(rho0, x)) / _eta2(rho0, x)) * \
               np.exp(veta * (1.0 - _eta(rho0, x)))

    def _derivative(self, x):
        return (np.exp(-1.5 * (-1.0 + self.k0_prime) * (-1.0 + _eta(x, self.rho0))) *
                self.k0 *
                (-2.0 +
                 (2.5 - 1.5 * self.k0_prime) * _eta(x, self.rho0) +
                 (-1.5 + 1.5 * self.k0_prime) * _eta2(x, self.rho0))) / \
               (_eta5(x, self.rho0) * self.rho0)

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return 2.0 * \
               np.exp(-1.5 * (-1.0 + self.k0_prime) * (-1.0 + _eta(x, self.rho0))) * \
               self.k0 * \
               (5.0 +
                3.0 * self.k0_prime(-1.0 + _eta(x, self.rho0)) -
                3.0 * _eta(x, self.rho0)) * self.rho0 / (self.k0_prime - 1.0) ** 2.0

    def _pressure_integral(self, x):
        return -(1.0 / (1280.0 * _eta5(x, self.rho0) * self.rho0)) * \
               9.0 * \
               np.exp(-1.5 * (-1.0 + _eta(x, self.rho0)) * (-1.0 + self.k0_prime)) * \
               self.k0 * \
               (2.0 *
                (128.0 -
                 16.0 * _eta(x, self.rho0) * (7.0 + 3.0 * self.k0_prime) -
                 6.0 * _eta3(x, self.rho0) * (-1.0 + self.k0_prime) ** 2.0 * (7.0 + 3.0 * self.k0_prime) +
                 9.0 * _eta4(x, self.rho0) * (-1.0 + self.k0_prime) ** 3.0 * (7.0 + 3.0 * self.k0_prime) +
                 8.0 * _eta2(x, self.rho0) * (-7.0 + 4 * self.k0_prime + 3 * self.k0_prime ** 2.0)) +
                27.0 *
                np.exp(1.5 * _eta(x, self.rho0) * (-1.0 + self.k0_prime)) *
                _eta5(x, self.rho0) *
                (-1.0 + self.k0_prime) ** 4.0 *
                (7.0 + 3.0 * self.k0_prime) *
                spec.expi(-1.5 * _eta(x, self.rho0) * (-1.0 + self.k0_prime)))


factory.register('vinet', Vinet)


class Murnaghan(Pressure_Fit_Class):
    """ Fit function for Murnaghan EOS (MU2)
        FD Stacey, et al., "Finite Strain Theories and Comparisons with Seismological Data
            Geophysical Surveys 4 (1981) 189-232
    """

    def __init__(self, args, name):
        super(Murnaghan, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        return (k0 / k0_prime) * \
               (np.power(_rno_norm(x, rho0), k0_prime) - 1.0)

    def _derivative(self, x):
        # noinspection PyTypeChecker
        return ((self.k0 / self.rho0) *
                np.power(_rno_norm(x, self.rho0), (self.k0_prime - 1.0)))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return (self.k0 *
                self.rho0 *
                _rno_norm(x, self.rho0)) * \
               (np.power(_rno_norm(x, self.rho0), self.k0_prime) - self.k0_prime - 1.0) / \
               (self.k0_prime * (1.0 + self.k0_prime))

    def _pressure_integral(self, x):
        return self.k0 * \
               (np.power(_rno_norm(x, self.rho0), self.k0_prime) - self.k0_prime - 1.0) / \
               (self.k0_prime * self.rho0 * _rno_norm(x, self.rho0) * (self.k0_prime - 1.0))


factory.register('murnaghan', Murnaghan)


# ------------------------------------------------------------------------------
# Pressure Fitter from Sandia wanted by Christine; hasn't been
# tested yet in MEOS


class SandiaPC(Pressure_Fit_Class):
    def __init__(self, args, name):
        super(SandiaPC, self).__init__(args, name)
        self.kneg1 = args.kneg1
        self.kneg0 = args.kneg0
        self.k1 = args.k1
        self.k2 = args.k2
        self.k3 = args.k3
        self.k4 = args.k4
        self.k5 = args.k5
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.kneg1, self.kneg0, self.k1, self.k2, self.k3, self.k4, self.k5) = coeffs[:-1]

    def _get_coefficients(self):
        return self.kneg1, self.kneg0, self.k1, self.k2, self.k3, self.k4, self.k5, self.rho0

    def _print_coefficients(self):
        print("kneg1 = {};".format(self.kneg1))
        print("k0 = {};".format(self.kneg0))
        print("k1 = {};".format(self.k1))
        print("k2 = {};".format(self.k2))
        print("k3 = {};".format(self.k3))
        print("k4 = {};".format(self.k4))
        print("k5 = {};".format(self.k5))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("kneg1 = {};\n".format(self.kneg1)) + ("k0 = {};\n".format(self.kneg0)) +
                                  ("k1 = {};\n".format(self.k1)) + ("k2 = {};\n".format(self.k2)) +
                                  ("k3 = {};\n".format(self.k3)) + ("k4 = {};\n".format(self.k4)) +
                                  ("k5 = {};\n".format(self.k5)) + ("rho0 = {};".format(self.rho0)))

    def bounds(self):
        """ Since rho0 is required to be given by the user, rho0 should
            not be changed during the parameter optimizing process.
            This is done by making the bounds of rho0 to change so
            small that the user won't see it and won't affect
            calculations
        """
        return ([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, self.rho0],
                [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, self.rho0 + 1e10])

    @staticmethod
    def _f(x, *coeffs):
        (kneg1, kneg0, k1, k2, k3, k4, k5, rho0) = coeffs
        y = _eta(x, rho0)
        return kneg1 * y ** (-1) + \
               kneg0 * 1 + \
               k1 * y + \
               k2 * y ** 2 + \
               k3 * y ** 3 + \
               k4 * y ** 4 + \
               k5 * y ** 5

    def _derivative(self, x):
        y = _eta(x, rho0)
        return -(self.kneg1 * y ** (-2)) + \
               self.k1 / self.rho0 + \
               2 * self.k2 * y + \
               3 * self.k3 * y ** 2 + \
               4 * self.k4 * y ** 3 + \
               5 * self.k5 * y ** 4

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x): raise RuntimeError(
        "NOT YET IMPLEMENTED")  # Needed to implement the abstract method, even though it's not of use for this fitter

    def _pressure_integral(self, x): raise RuntimeError(
        "NOT YET IMPLEMENTED")  # Needed to implement the abstract method, even though it's not of use for this fitter


factory.register('sandiapc', SandiaPC)


# ------------------------------------------------------------------------------
# Working EOS models:

# Energy curve is integral of fit function divided by rho^2 for pressure
# graphs


class Anton_Schmidt(Pressure_Fit_Class):
    """ Fit function for Anton-Schmidt EOS
        H. Anton and P.C. Schmidt, "Theoretical investigations of the elastic
        constants in Laves phases"
            Intermetallics 5 (1997) 449-465
    """

    def __init__(self, args, name):
        super(Anton_Schmidt, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        n = -.5 * k0_prime
        # noinspection PyTypeChecker
        return k0 * \
               np.power(_rno_norm(x, rho0), -n) * \
               np.log(_rno_norm(x, rho0))

    #        return k0 * \
    #               np.power(_rno_norm(x, rho0), -n) * \
    #               np.log(_rno_norm(x, rho0))

    def _derivative(self, x):
        n = -.5 * self.k0_prime
        # noinspection PyTypeChecker
        return (self.k0 / self.rho0) * \
               (np.power(_rno_norm(x, self.rho0), -n - 1.0) *
                (1.0 - n * np.log(_rno_norm(x, self.rho0))))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        n = -.5 * self.k0_prime
        # noinspection PyTypeChecker
        return -(self.k0 *
                 np.power(_rno_norm(x, self.rho0), 1 - n) *
                 self.rho0 *
                 (1.0 + (-1.0 + n) * np.log(_rno_norm(x, self.rho0))) /
                 (-1.0 + n) ** 2.0)

    def _pressure_integral(self, x):
        n = -.5 * self.k0_prime
        # noinspection PyTypeChecker
        return -(self.k0 *
                 np.power(_rno_norm(x, self.rho0), -1 - n) *
                 (-1.0 + (1.0 + n) * np.log(_rno_norm(x, self.rho0))) /
                 (self.rho0 * (1.0 + n) ** 2.0))


factory.register('anton', Anton_Schmidt)


class Bardeen(Pressure_Fit_Class):
    """ Fit function for Bardeen EOS
        I. Alkammash, "Evaluation of pressure and bulk modulus for alkali
        halides under high pressure and temperature using different EOS"
            Journal of the Association of Arab Universities for Basic and
            Applied Sciences 14 (2013) 38-45
    """

    def __init__(self, args, name):
        super(Bardeen, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        return 3.0 * \
               k0 * \
               (_eta5(x, rho0) - _eta4(x, rho0)) * \
               (1.0 + 1.5 * (k0_prime - 3.0) * (_eta(x, rho0) - 1.0))

    def _derivative(self, x):
        return ((-self.k0 / self.rho0) *
                ((27.0 - 9.0 * self.k0_prime) * _eta3(x, self.rho0) +
                 (15.0 * self.k0_prime - 50.0) * _eta2(x, self.rho0) +
                 (22.0 - 6.0 * self.k0_prime) * _eta(x, self.rho0)))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return (3.0 * self.k0 * self.rho0 * _eta7(x, self.rho0) / 56.0) * \
               (-132.0 +
                36.0 * self.k0_prime +
                21.0 * (10.0 - 3.0 * self.k0_prime) * _eta(x, self.rho0) +
                28.0 * (-3.0 + self.k0_prime) * _eta2(x, self.rho0))

    def _pressure_integral(self, x):
        return (9.0 * self.k0 * ((3.0 * self.k0_prime - 11.0) * _eta(x, self.rho0) +
                                 (10.0 - 3.0 * self.k0_prime) * _eta2(x, self.rho0) +
                                 (self.k0_prime - 3.0) * _eta3(x, self.rho0))) / \
               (2 * self.rho0)


factory.register('bardeen', Bardeen)


class Birch_Murnaghan2(Pressure_Fit_Class):
    """ Fit function for second-order Murnaghan EOS (BE1?)
    R. Jeanloz, "Finite-Strian Equation of State for High-Pressure Phases"
        Geophysical Research Letters 8 (12) 1219-1222 (1981)
    J.-P. Poirier and A. Tarantola, "A logarithmic equation of state"
        Physics of the Earth and Planetary Interiors 109 (1998) 1-8
    """

    def __init__(self, args, name):
        super(Birch_Murnaghan2, self).__init__(args, name)
        self.k0 = args.k0
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, rho0) = coeffs
        return (1.5 * k0) * \
               (_eta7(x, rho0) - _eta5(x, rho0))

    def _derivative(self, x):
        return ((self.k0 * _eta2(x, self.rho0) / (2.0 * self.rho0)) *
                (7.0 * _eta2(x, self.rho0) - 5.0))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return (9.0 * self.k0 * self.rho0 / 80.0) * \
               (4.0 * _eta2(x, self.rho0) - 5.0) * \
               _eta8(x, self.rho0)

    def _pressure_integral(self, x):
        return (9.0 * self.k0 / (self.rho0 * 8.0)) * \
               (_eta2(x, self.rho0) - 2.0) * \
               _eta2(x, self.rho0)


factory.register('birch2', Birch_Murnaghan2)


class Johnson_Holmquist(Pressure_Fit_Class):
    """ Fit function for Johnson-Holmquist EOS, compression mode only
        G. McIntosh, "The Johnson-Holmquist Ceramic Model as used in LS-DYNA2D"
            Report DREV-TM-9822 Department of National Defence, Canada
    """

    def __init__(self, args, name):
        super(Johnson_Holmquist, self).__init__(args, name)
        self.k1 = args.k1
        self.k2 = args.k2
        self.k3 = args.k3
        self.delta_p = args.delta_p
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k1, self.k2, self.k3, self.delta_p, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k1, self.k2, self.k3, self.delta_p, self.rho0

    def _print_coefficients(self):
        print("k1 = {};".format(self.k1))
        print("k2 = {};".format(self.k2))
        print("k3 = {};".format(self.k3))
        print("delta_p = {};".format(self.delta_p))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("k1 = {};\n".format(self.k1)) + ("k2 = {};\n".format(self.k2)) +
                                  ("k3 = {};\n".format(self.k3)) + ("delta_p = {};\n".format(self.delta_p)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k1, k2, k3, delta_p, rho0) = coeffs
        return k1 * _mu(x, rho0) + \
               k2 * _mu2(x, rho0) + \
               k3 * _mu3(x, rho0) + \
               delta_p

    def _derivative(self, x):
        return ((self.k1 +
                 2 * self.k2 * _mu(x, self.rho0) +
                 3 * self.k3 * _mu2(x, self.rho0)) /
                self.rho0)

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return (self.rho0 / 12.0) * \
               (4 * self.k2 -
                3 * self.k3 +
                4 * self.k2 * _mu3(x, self.rho0) +
                3 * self.k3 * _mu4(x, self.rho0) +
                12 * self.delta_p * (1 + _mu(x, self.rho0)) +
                6 * self.k1 * (_mu2(x, self.rho0) - 1))

    def _pressure_integral(self, x):
        return (-2 * self.delta_p +
                2 * self.k1 * (1 + np.log(x) * (1 + _mu(x, self.rho0))) +
                2 * self.k2 * (-2 * np.log(x) * (1 + _mu(x, self.rho0)) +
                               _mu(x, self.rho0) * (2 + _mu(x, self.rho0))) +
                self.k3 * (-3 + 6 * np.log(x) * (1 + _mu(x, self.rho0)) +
                           _mu(x, self.rho0) * (-9 + (-3 + _mu(x, self.rho0)) * _mu(x, self.rho0)))) / \
               (2 * (1 + _mu(x, self.rho0)) * self.rho0)


factory.register('johnson', Johnson_Holmquist)


class Kumari_Dass(Pressure_Fit_Class):
    """ Fit function for Kumari-Dass Equation of State
        UCD PhysWiki: http://physwiki.ucdavis.edu/Condensed_Matter/
        Equations_of_state/Kumari-Dass_equation_of_state
    """

    def __init__(self, args, name):
        super(Kumari_Dass, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.lam = args.lam
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.lam, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.lam, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("lam = {};".format(self.lam))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("lam = {};\n".format(self.lam)) + ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, lam, rho0) = coeffs
        return (1.0 / lam) * \
               (lam *
                k0 *
                np.power(_rno_norm(x, rho0), (lam * k0 - k0_prime) + k0_prime) /
                (lam * k0 + k0_prime))

    def _derivative(self, x):
        # noinspection PyTypeChecker
        return ((self.k0 / self.rho0) *
                (-self.k0_prime + self.k0 * self.lam) *
                np.power(_rno_norm(x, self.rho0), (self.lam * self.k0 - self.k0_prime - 1.0)) /
                (self.k0_prime + self.k0 * self.lam))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        # noinspection PyTypeChecker
        return ((self.k0 *
                 self.lam *
                 np.power(_rno_norm(x, self.rho0), (1.0 - self.k0_prime + self.k0 * self.lam)) *
                 self.rho0) /
                (1.0 - self.k0_prime + self.k0 * self.lam) +
                self.k0_prime * _rno_norm(x, self.rho0) * self.rho0) / \
               (self.k0_prime * self.lam + self.k0 * self.lam ** 2.0)

    def _pressure_integral(self, x):
        return ((-self.k0_prime -
                 self.k0_prime ** 2.0 +
                 self.k0 * self.k0_prime * self.lam -
                 self.k0 * self.lam * np.power(_rno_norm(x, self.rho0), (self.k0 * self.lam - self.k0_prime))) /
                (self.lam * _rno_norm(x, self.rho0) *
                 self.rho0 * (self.k0_prime + self.k0 * self.lam) *
                 (1.0 + self.k0_prime - self.k0 * self.lam)))


factory.register('kumari', Kumari_Dass)


class Logarithmic2(Pressure_Fit_Class):
    """ Fit function for second-order Logarithmic EOS
        J.-P. Poirier and A. Tarantola, "A logarithmic equation of state"
            Physics of the Earth and Planetary Interiors 109 (1998) 1-8
    """

    def __init__(self, args, name):
        super(Logarithmic2, self).__init__(args, name)
        self.k0 = args.k0
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, rho0) = coeffs
        return k0 * \
               _rno_norm(x, rho0) * \
               np.log(_rno_norm(x, rho0))

    def _derivative(self, x):
        return (self.k0 / self.rho0) * \
               (1.0 + np.log(_rno_norm(x, self.rho0)))

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return .25 * self.k0 * \
               np.square(_rno_norm(x, self.rho0)) * \
               self.rho0 * (2.0 * np.log(_rno_norm(x, self.rho0)) - 1.0)

    def _pressure_integral(self, x):
        return self.k0 * \
               np.square(np.log(_rno_norm(x, self.rho0))) / \
               (2 * self.rho0)


factory.register('log2', Logarithmic2)


class Logarithmic3(Pressure_Fit_Class):
    """ Fit function for third-order Logarithmic EOS
        J.-P. Poirier and A. Tarantola, "A logarithmic equation of state"
            Physics of the Earth and Planetary Interiors 109 (1998) 1-8
    """

    def __init__(self, args, name):
        super(Logarithmic3, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        return k0 * \
               _rno_norm(x, rho0) * \
               (np.log(_rno_norm(x, rho0)) +
                (k0_prime / 2.0 - 1.0) * np.square(np.log(_rno_norm(x, rho0))))

    def _derivative(self, x):
        return self.k0 * \
               (2.0 +
                2.0 * (-1.0 + self.k0_prime) * np.log(_rno_norm(x, self.rho0)) +
                (-2.0 + self.k0_prime) * np.square(np.log(_rno_norm(x, self.rho0)))) / \
               (2.0 * self.rho0)

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return (self.k0 * np.square(_rno_norm(x, self.rho0)) * self.rho0 / 8.0) * \
               (-4.0 +
                self.k0_prime -
                2.0 * (self.k0_prime - 4.0) * np.log(_rno_norm(x, self.rho0)) +
                2.0 * (self.k0_prime - 2.0) * np.square(np.log(_rno_norm(x, self.rho0))))

    def _pressure_integral(self, x):
        return self.k0 * \
               np.square(np.log(_rno_norm(x, self.rho0))) * \
               (3.0 + (-2.0 + self.k0_prime) * np.log(_rno_norm(x, self.rho0))) / \
               (6.0 * self.rho0)


factory.register('log', Logarithmic3)


class Shankar(Pressure_Fit_Class):
    """ Fit function for Shankar EOS
        JC Bhatt, et al., "High Pressure Equation of State for Nanomaterials"
            ISRN Nanotechnology 2013 404920 1-6
    """

    def __init__(self, args, name):
        super(Shankar, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0) = coeffs
        return k0 * \
               _rno_norm(x, rho0) * \
               ((1.0 - 1.0 / _rno_norm(x, rho0)) +
                ((k0_prime + 1.0) / 2.0) * (1.0 - 1.0 / np.square(_rno_norm(x, rho0))))

    def _derivative(self, x):
        return (self.k0 * (1.0 +
                           self.k0_prime +
                           (3.0 + self.k0_prime) * np.square(_rno_norm(x, self.rho0)))) / \
               (2.0 * np.square(_rno_norm(x, self.rho0)) * self.rho0)

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        return .25 * \
               self.k0 * \
               self.rho0 * \
               (_rno_norm(x, self.rho0) * (-4.0 + (3.0 + self.k0_prime) * _rno_norm(x, self.rho0)) -
                2.0 * (1.0 + self.k0_prime) * np.log(_rno_norm(x, self.rho0) * self.rho0))

    def _pressure_integral(self, x):
        return (self.k0 *
                (1.0 +
                 self.k0_prime +
                 4.0 * _rno_norm(x, self.rho0) +
                 2.0 *
                 (3.0 + self.k0_prime) *
                 np.square(_rno_norm(x, self.rho0)) *
                 np.log(_rno_norm(x, self.rho0) * self.rho0))) / \
               (4.0 * np.square(_rno_norm(x, self.rho0)) * self.rho0)


factory.register('shank', Shankar)


# ------------------------------------------------------------------------------
# Incomplete EOS models


class Ap1(Pressure_Fit_Class):
    """ Fit function for first-order Adapted Polynomial EOS (AP1)
        W. Holzapfel, "Equations of state for regular solids"
        High Pressure Research 22(1) 209-216 (2002)
    """

    def __init__(self, args, name):
        super(Ap1, self).__init__(args, name)
        self.k0 = args.k0
        self.rho0 = args.rho0
        self.z = args.z

    def _set_coefficients(self, coeffs):
        (self.k0, self.rho0, self.z) = coeffs

    def _get_coefficients(self):
        return self.k0, self.rho0, self.z

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("rho0 = {};".format(self.rho0))
        print("z = {};".format(self.z))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("rho0 = {};\n".format(self.rho0)) +
                                              ("z = {};".format(self.z)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, rho0, z) = coeffs
        return 3.0 * \
               k0 * \
               ((1.0 - _eta(x, rho0)) / _eta5(x, rho0)) * \
               np.exp(_c1(k0, rho0, z) * (1.0 - _eta(x, rho0)))

    # noinspection PyUnusedLocal
    def _derivative(self, x):
        return np.exp(_c1(self.k0, self.rho0, self.z) * (1.0 - _eta(x, self.rho0))) * \
               self.k0 * \
               (-5.0 +
                _c1(self.k0, self.rho0, self.z) * (5.0 -
                                                   11.0 * _eta(x, self.rho0) +
                                                   6.0 * _eta2(x, self.rho0)) +
                4.0 * _eta(x, self.rho0)) / \
               (_eta8(x, self.rho0) * self.rho0)

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("Not implemented")

    def _pressure_integral(self, x):
        return (4.5 *
                self.k0 /
                (self.rho0 * _eta2(x, self.rho0)) *
                np.exp(_c1(self.k0, self.rho0, self.z) * (1.0 - _eta(x, self.rho0))) *
                (1.0 -
                 (_c1(self.k0, self.rho0, self.z) + 2.0) *
                 _eta(x, self.rho0) *
                 (1.0 -
                  (_c0x(x, self.k0, self.rho0, self.z) *
                   np.exp(_c0x(x, self.k0, self.rho0, self.z)) *
                   spec.expi(_c0x(x, self.k0, self.rho0, self.z))))))


factory.register('ap1', Ap1)


class Ap2(Pressure_Fit_Class):
    """ Fit function for second-order Adapted Polynomial EOS
        Jun Jiuxun, et al. "Equation of state for solids with high accuracy
        and satisfying the limitation condition at high pressure"
        Physica B 371 (2006) 257-271
    """

    def __init__(self, args, name):
        super(Ap2, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0
        self.z = args.z

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0, self.z) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0, self.z

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        print("z = {};".format(self.z))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};\n".format(self.rho0)) + ("z = {};".format(self.z)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0, z) = coeffs
        return 3 * k0 * \
               ((1 - _eta(x, rho0)) / _eta5(x, rho0)) * \
               np.exp(_c1(x, k0, z) * (1 - _eta(x, rho0))) * \
               (1 + _c2(x, k0, k0_prime, z) * _eta(x, rho0) * (1 - _eta(x, rho0)))

    def _derivative(self, x):
        # very complicated
        raise RuntimeError("Not implemented")

    def _second_derivative(self, x):
        raise RuntimeError("Not implemented")

    def _energy_integral(self, x):
        # no closed form?
        raise RuntimeError("Not implemented")

    def _pressure_integral(self, x):
        raise RuntimeError("Not implemented")


factory.register('ap2', Ap2)


# def ap3(self, x, k0, k0_prime, k0_prime_prime, z, rho0):
#     """ Fit function for third-order Adapted Polynomial EOS
#         Jun Jiuxun, et al. "Equation of state for solids with high accuracy
#         and satisfying the limitation condition at high pressure"
#         Physica B 371 (2006) 257-271
#     """
#     a0 = .02337  # GPa nm^5
#     _pFG0 = a0 * np.power(z * x, (5 / 3))
#     _c1 = -np.log(3 * k0 / _pFG0)
#     _c2 = 1.5 * (k0_prime - 3) - _c1
#     c3 = (1 / 6) * (-61 + 9 * np.square(k0_prime) + k0 * k0_prime_prime -
#                     24 * _c1 - 3 * np.square(_c1) - 6 * _c1 * _c2 - 18 * _c2)
#     result = 3 * k0 * ((1 - _eta(x, rho0)) / _eta5(x, rho0)) * np.exp(_c1 * (1 - _eta(x, rho0))) * \
#         (1 + _c2 * _eta(x, rho0) * (1 - _eta(x, rho0)) +
#          c3 * _eta(x, rho0) * np.square(1 - _eta(x, rho0)))
#     return result

# ------------------------------------------------------------------------------


class old_CurveInteractor(object):
    def __init__(self):
        self.args = None
        self.regions = None
        self.lines = None

    # Private methods: =========================================================


def _print_polynomial(coeffs):
    print('Calculated polynomial is:')
    length = len(coeffs)
    calculated_polynomial = ""
    for index in range(0, length):
        x_order = length - index - 1
        if index == 0:
            print('%.3E * x**%d ' % (coeffs[index], x_order), end='')
            calculated_polynomial += '%.3E * x**%d ' % (coeffs[index], x_order)
        elif x_order > 0:
            print('%+.3E * x**%d ' % (coeffs[index], x_order), end='')
            calculated_polynomial += '%+.3E * x**%d ' % (coeffs[index], x_order)
        else:
            print('%+.3E' % coeffs[index])
            calculated_polynomial += '%+.3E' % coeffs[index]
    update_fitter_info_window(-2, False, calculated_polynomial)
    print()


# Old polynomial fit calculator:


class Poly_Original(PolyBase):
    _name_prefix = 'poly'

    """ This is substantially different from the rest of the fit classes above,
    but presents the same public methods as Base_Fit_Class:
    """

    def __init__(self, args, name):
        #        self._order = args.polynomial_order
        self._order = int(name[len(self.name_prefix):])
        self._f = None
        self._is_first_fit = True
        self.derivative_scale = args.derivative_scale
        self.second_derivative_scale = args.second_derivative_scale
        self.integral_scale = args.integral_scale
        self.x_integral_ref = args.x_integral_ref
        self.y_integral_ref = args.y_integral_ref

    def _set_poly(self, coeffs):
        # Create a polynomial function using the coefficients, to be used to
        # calculate a new curve:
        self._f = np.poly1d(coeffs)
        self._der = np.polyder(self._f)
        self._int = np.polyint(self._f)
        self._scnd_der = np.polyder(self._der)

    def fit_to_points(self, points):
        """ Calculate the coefficients and create the polynomial function.
        Doesn't seem to mind if the points are not in order by x.

             Do a least-squares polynomial fit using Ax = y
                A: Vandermonde matrix, derived from x values
                x: coefficients vector, to be found
                y: y values vector

        :param points:   list of x, y tuples
        :return:         polynomial function that returns a "good" y for an x
                         within the given x_values range
        """

        x_values = [point[0] for point in points]
        y_values = [point[1] for point in points]
        # Form the Vandermonde matrix:
        # (x[0]^0, x[0]^1, ... x[0]^order)
        # (x[1]^0, x[1]^1, ... x[1]^order)
        # ...
        # (x[n]^0, x[n]^1, ... x[n]^order)
        A = np.vander(x_values, self._order + 1)

        # Find the x (coeffs) that minimizes the norm of Ax-y
        # _vars are unused
        coeffs, _residuals, _rank, _sing_vals = np.linalg.lstsq(A, y_values, rcond=-1)

        _print_polynomial(coeffs)

        self._set_poly(coeffs)

    def guess_coefficients(self, points):
        pass

    def func(self, x):
        return self._f(x)

    def derivative(self, x):
        return self.derivative_scale * self._der(x)

    def second_derivative(self, x):
        return self.second_derivative_scale * self._scnd_der(x)

    def integral(self, x):
        return self.integral_scale * \
               (self._int(x) +
                (self.y_integral_ref - self._int(self.x_integral_ref)))


factory.register(Poly_Original.name_prefix, Poly_Original)


# --------------------------------------------------------------------------------
# Energy EOS models
# These modules are for energy instead of pressure. All models were obtained from
# The main MEOS Code. No integrals or derivatives are currently calculated. All 
# energy model names begin with an 'E'.

class EBirch_Murnaghan3(Energy_Fit_Class):
    def __init__(self, args, name):
        super(EBirch_Murnaghan3, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0
        self.e0 = args.e0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0, self.e0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0, self.e0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        print("E0 = {};".format(self.e0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};\n".format(self.rho0)) + ("E0 = {};".format(self.e0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0, e0) = coeffs
        multiplier = ((9 * (1 / rho0) * k0) / 16)
        term1_0 = pow(x / rho0, (2.0 / 3.0))
        term1_1 = term1_0 - 1
        term1 = pow(term1_1, 3) * k0_prime
        term2_1 = pow(term1_1, 2)
        term2_2 = 6 - (4 * term1_0)
        term2 = term2_1 * term2_2
        block = multiplier * (term1 + term2)
        fcold = e0 + block
        return fcold

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('ebirch3', EBirch_Murnaghan3)


class EBirch_Murnaghan4(Energy_Fit_Class):
    def __init__(self, args, name):
        super(EBirch_Murnaghan4, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.k0_prime_prime = args.k0_prime_prime
        self.rho0 = args.rho0
        self.e0 = args.e0
        self.EBMurn3 = EBirch_Murnaghan3(args, name)

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.k0_prime_prime, self.rho0, self.e0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.k0_prime_prime, self.rho0, self.e0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("Bpp = {};".format(self.k0_prime_prime))
        print("rho0 = {};".format(self.rho0))
        print("E0 = {};".format(self.e0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("Bpp = {};\n".format(self.k0_prime_prime)) + ("rho0 = {};\n".format(self.rho0)) +
                                  ("E0 = {};".format(self.e0)))

    # Birch_Murnaghan4 doesn't seem to converge well, use BMurn3 as inital guess.
    # New guessing algorithm works well now, so taking out BMurn3 stuff
    #    def guess_coefficients(self, points):
    #        self.EBMurn3.guess_coefficients(points)
    #        self.EBMurn3.fit_to_points(points)
    #        (self.k0, self.k0_prime, self.rho0, self.e0) = self.EBMurn3._get_coefficients()

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, k0_prime_prime, rho0, e0) = coeffs
        strain = (1.0 / 2.0) * (pow((x / rho0), (2.0 / 3.0)) - 1.0)
        strain_2 = pow(strain, 2)
        inner2_term = k0 * k0_prime_prime + k0_prime * (k0_prime - 7.0) + 143.0 / 9.0
        inner_term = 1.0 + (k0_prime - 4.0) * strain + (3.0 / 4.0) * inner2_term * strain_2
        fcold = e0 + (9.0 / 2.0) * k0 * (1 / rho0) * strain_2 * inner_term
        return fcold

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('ebirch4', EBirch_Murnaghan4)


class EMurnaghan(Energy_Fit_Class):
    def __init__(self, args, name):
        super(EMurnaghan, self).__init__(args, name)
        self.k0 = args.k0
        self.k0_prime = args.k0_prime
        self.rho0 = args.rho0
        self.e0 = args.e0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0, self.e0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0, self.e0

    def _print_coefficients(self):
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        print("E0 = {};".format(self.e0))
        update_fitter_info_window(-1, False, ("B0 = {};\n".format(self.k0)) + ("Bp = {};\n".format(self.k0_prime)) +
                                  ("rho0 = {};\n".format(self.rho0)) + ("E0 = {};".format(self.e0)))

    @staticmethod
    def _f(x, *coeffs):
        (k0, k0_prime, rho0, e0) = coeffs
        term3 = pow(x / rho0, k0_prime) / (k0_prime - 1.0) + 1
        fcold = e0 + k0 * x / k0_prime * term3 - k0 * (1 / rho0) / (k0_prime - 1.0)
        return fcold

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('emurnaghan', EMurnaghan)


# ----------------------------------------------------------------------------
# eSeries
# This model is specifically a refined fit type.
# It can be any order between 4 and 12.

class ESeries(Energy_Fit_Class):
    def __init__(self, args, name):
        super(ESeries, self).__init__(args, name)
        self.rho0 = args.rho0
        self.cn = args.cn
        if len(name) > 7:
            self.order = int(name[7:])
        else:
            # Set default order to 7
            self.order = 7

    def _set_coefficients(self, coeffs):
        (self.rho0, self.cn[0], self.cn[1], self.cn[2],
         self.cn[3], self.cn[4], self.cn[5], self.cn[6], self.cn[7],
         self.cn[8], self.cn[9], self.cn[10], self.cn[11],
         self.order) = coeffs

    def _get_coefficients(self):
        if self.rho0 == None:
            self.rho0 = 5.0
        return (self.rho0, self.cn[0], self.cn[1],
                self.cn[2], self.cn[3], self.cn[4], self.cn[5],
                self.cn[6], self.cn[7], self.cn[8], self.cn[9],
                self.cn[10], self.cn[11], self.order)

    def _print_coefficients(self):
        print("rho0 = {};".format(self.rho0))
        fitter_info_text = ("rho0 = {};\n".format(self.rho0))
        for i in range(0, int(self.order) - 3):
            print("C{} = {};".format(i + 4, self.cn[i]))
            fitter_info_text += ("C{} = {};\n".format(i + 4, self.cn[i]))
        update_fitter_info_window(-1, False, fitter_info_window[0:len(fitter_info_text)-1])

    @staticmethod
    def _f(x, *coeffs):
        cn = [None] * 12  # 12 is the maximum order of eseries
        (rho0, cn[0], cn[1], cn[2], cn[3],
         cn[4], cn[5], cn[6], cn[7], cn[8],
         cn[9], cn[10], cn[11], order) = coeffs

        answer = 0
        for i in range(0, int(order) - 3):
            answer += ((float(cn[i]) / float(factorial(i + 4))) *
                       pow(0.5 * (pow(np.asarray(x) / rho0,
                                      2.0 / 3.0) - 1.0), i + 4))

        return answer

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('eseries', ESeries)

# Allow option for none, doesn't acually use Poly_Original
factory.register('none', Poly_Original)


class EVinet(Energy_Fit_Class):
    """ Fit function for third-order Vinet Energy EOS
        R. Jeanloz, "Universal equation of state"
            Physical Review B 38(1) 805
    """

    def __init__(self, args, name):
        super(EVinet, self).__init__(args, name)
        self.k0 = args.k0  # B0 (Bulk Modulus)
        self.k0_prime = args.k0_prime  # B'
        self.rho0 = args.rho0  # rho0
        self.e0 = args.e0  # E0

    def _set_coefficients(self, coeffs):
        (self.k0, self.k0_prime, self.rho0, self.e0) = coeffs

    def _get_coefficients(self):
        return self.k0, self.k0_prime, self.rho0, self.e0

    def _print_coefficients(self):
        print("E0 = {};".format(self.e0))
        print("B0 = {};".format(self.k0))
        print("Bp = {};".format(self.k0_prime))
        print("rho0 = {};".format(self.rho0))
        update_fitter_info_window(-1, False, ("E0 = {};\n".format(self.e0)) + ("B0 = {};\n".format(self.k0)) +
                                  ("Bp = {};\n".format(self.k0_prime)) + ("rho0 = {};".format(self.rho0)))

    @staticmethod
    def _f(rho, *coeffs):
        (k0, k0_prime, rho0, e0) = coeffs
        #        x = pow(rho0/rho, 1.0/3.0 )
        #        term1 = np.exp(-1.5 * (k0_prime - 1.0) *(x  - 1.0));
        #        term2 = 5.0+3.0* k0_prime * (x-1) - 3.0*x;
        #        term3 = (2.0*k0) / (rho0 * pow( (k0_prime - 1.0), 2.0));
        #        term4 = (4.0*k0) / (rho0 * pow( (k0_prime - 1.0), 2.0));
        #        fcold = e0 + term4 -term3 * term2 * term1;
        #        return fcold;

        # Christine Version:
        #        for thisRho in rho:
        #            print(thisRho, rho0, rho0/thisRho, pow(rho0/thisRho, 1.0/3.0 ))
        x = 1.5 * (k0_prime - 1) * (pow(rho0 / rho, 1.0 / 3.0) - 1);
        term1 = (4 * k0) / (rho0 * pow(k0_prime - 1, 2.0));
        term2 = 1 - (1 + x) * np.exp(-x);
        return term1 * term2 + e0

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('evinet', EVinet)


# -----------------------------------------------------------------------------------------
# High Pressure Correction Equation
# This model doesn't fit any other descriptions because it is specifically a refined fit
# type for the energy version of birch murnaghan. 

class EHighP(Energy_Fit_Class):
    def __init__(self, args, name):
        super(EHighP, self).__init__(args, name)
        self.rho0 = args.rho0
        self.c1 = args.c1
        self.c2 = args.c2
        self.c3 = args.c3
        self.c4 = args.c4

    def _set_coefficients(self, coeffs):
        (self.rho0, self.c1, self.c2, self.c3, self.c4) = coeffs

    def _get_coefficients(self):
        if self.rho0 == None:
            self.rho0 = 1.0
        return self.rho0, self.c1, self.c2, self.c3, self.c4

    def _print_coefficients(self):
        print("rho0 = {};".format(self.rho0))
        print("c1 = {};".format(self.c1))
        print("c2 = {};".format(self.c2))
        print("c3 = {};".format(self.c3))
        print("c4 = {};".format(self.c4))
        update_fitter_info_window(-1, False, ("rho0 = {};\n".format(self.rho0)) + ("c1 = {};\n".format(self.c1)) +
                                  ("c2 = {};\n".format(self.c2)) + ("c3 = {};\n".format(self.c3)) +
                                  ("c4 = {};".format(self.c4)))

    @staticmethod
    def _f(x, *coeffs):
        (rho0, c1, c2, c3, c4) = coeffs
        xx = .5 * (((x / rho0) ** (2.0 / 3.0)) - 1.0)
        answer = (c1 * xx + c2 * ((xx ** 2) / 2) +
                  c3 * ((xx ** 3) / (3 * 2)) + c4 * ((xx ** 4) / (4 * 3 * 2)))
        return answer

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('highp', EHighP)


class ThetaBP(Base_Fit_Class):
    """ Fit function for theta > rho0  
        Leonid Burakovsky and Dean L Preston, "Analytic model of the Gruneisen parameter all densities"
             J. Phys. C 65, 1581 2004
    """

    def __init__(self, args, name):
        super(ThetaBP, self).__init__(args, name)
        self.theta0 = args.theta0  # theta at rho0
        self.rho0 = args.rho0  # rho0
        self.q = args.q  # q exponent
        self.c1 = args.c1  # coefficent 1
        self.c2 = args.c2  # coefficent 2

    def _set_coefficients(self, coeffs):
        (self.theta0, self.rho0, self.q, self.c1, self.c2) = coeffs

    def _get_coefficients(self):
        return self.theta0, self.rho0, self.q, self.c1, self.c2

    def _print_coefficients(self):
        print("theta0 = {};".format(self.theta0))
        print("rho0 = {};".format(self.rho0))
        print("q = {};".format(self.q))
        print("c1 = {};".format(self.c1))
        print("c2 = {};".format(self.c2))
        update_fitter_info_window(-1, False, ("theta0 = {};\n".format(self.theta0)) +
                                  ("rho0 = {};\n".format(self.rho0)) + ("q = {};\n".format(self.q)) +
                                  ("c1 = {};\n".format(self.c1)) + ("c2 = {};".format(self.c2)))

    def guess_coefficients(self,
                           points):  # We don't have a guess method for ThetaBP yet.  It only handles the compression side
        if (self.theta0 == None):
            self.theta0 = 1000
        if (self.rho0 == None):
            self.rho0 = 3.0
        self.q = 1.5
        self.c1 = 1.0
        self.c2 = 1.0

    @staticmethod
    def _f(rho, *coeffs):
        (theta0, rho0, q, c1, c2) = coeffs
        # theta = theta0 * sqrt(rho/rho0) * np.exp( -3 * c1 * ( pow(rho, -1.0/3.0) - pow(rho0, -1.0/3.0) ) - (c2 / q) * (pow(rho, -q) - pow(rho0, -q) ))
        inner_term1 = -3 * c1 * (pow(rho, -1.0 / 3.0) - pow(rho0, -1.0 / 3.0))
        inner_term2 = (c2 / q) * (pow(rho, -q) - pow(rho0, -q))
        exp_term = np.exp(inner_term1 - inner_term2)
        theta = theta0 * np.sqrt(rho / rho0) * exp_term

        return theta

    def _derivative(self, x):
        term1 = self.c1 / pow(x, 1.0 / 3.0)
        term2 = self.c2 / pow(x, self.q)
        return 0.5 + term1 + term2

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('thetabp', ThetaBP)


class GammaRho(Base_Fit_Class):
    """ Fit function for gamma on density g/cc  
        Reference?
    """

    def __init__(self, args, name):
        super(GammaRho, self).__init__(args, name)
        self.gamma0 = args.gamma0  # gamma at rho0
        self.rho0 = args.rho0  # rho0
        self.q = args.q  # q exponent
        self.c1 = args.c1  # coefficent 1
        self.c2 = args.c2  # coefficent 2

    def _set_coefficients(self, coeffs):
        (self.gamma0, self.rho0, self.q, self.c1, self.c2) = coeffs

    def _get_coefficients(self):
        return self.gamma0, self.rho0, self.q, self.c1, self.c2

    def _print_coefficients(self):
        print("gamma0 = {};".format(self.gamma0))
        print("rho0 = {};".format(self.rho0))
        print("q = {};".format(self.q))
        print("c1 = {};".format(self.c1))
        print("c2 = {};".format(self.c2))
        update_fitter_info_window(-1, False, ("gamma0 = {};\n".format(self.gamma0)) +
                                  ("rho0 = {};\n".format(self.rho0)) + ("q = {};\n".format(self.q)) +
                                  ("c1 = {};\n".format(self.c1)) + ("c2 = {};".format(self.c2)))

    def guess_coefficients(self, points):
        # If we don't have either gamma0 or rho0, just up the first point for both.
        if (self.gamma0 == None and self.rho0 == None):
            self.rho0 = points[0][0]
            self.gamma0 = points[0][1]
        elif (self.gamma0 == None):  # If we have rho0 but no gamma, find gamma0 at that point
            x_values = [point[0] for point in points]
            y_values = [point[1] for point in points]
            interpFunc = interp.interp1d(x_values, y_values);
            self.gamma0 = interpFunc(rh0)
        elif (self.rho0 == None):
            raise RuntimeError("Please define rho0 if gamma0 is defined")
        if (self.q == None):  # Calculate a guess for q
            yy = (self.gamma0 / self.c2) - (1 / (2 * self.c2)) - ((self.c1 / self.c2) * pow(self.rho0, -1.0 / 3.0))
            self.q = -math.log(yy) / math.log(self.rho0)

    @staticmethod
    def _f(rho, *coeffs):
        (gamma0, rho0, q, c1, c2) = coeffs
        term1 = c1 / pow(rho, 1.0 / 3.0)
        term2 = c2 / pow(rho, q)
        return 0.5 + term1 + term2

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('gammaRho', GammaRho)


class GammaV(Base_Fit_Class):
    """ Fit function for gamma on volume cc/g
        We just treat rho0 and V0 here, rather than using a different name.  Pass in V0 under the name rho0
        Reference?
    """

    def __init__(self, args, name):
        super(GammaV, self).__init__(args, name)
        self.gamma0 = args.gamma0  # theta at v0
        self.rho0 = args.rho0  # v0
        self.q = args.q  # q exponent
        self.c1 = args.c1  # coefficent 1
        self.c2 = args.c2  # coefficent 2

    def _set_coefficients(self, coeffs):
        (self.gamma0, self.rho0, self.q, self.c1, self.c2) = coeffs

    def _get_coefficients(self):
        return self.gamma0, self.rho0, self.q, self.c1, self.c2

    def _print_coefficients(self):
        print("gamma0 = {};".format(self.gamma0))
        print("V0 = {};".format(self.rho0))
        print("q = {};".format(self.q))
        print("c1 = {};".format(self.c1))
        print("c2 = {};".format(self.c2))
        update_fitter_info_window(-1, False, ("gamma0 = {};\n".format(self.gamma0)) +
                                  ("V0 = {};\n".format(self.rho0)) + ("q = {};\n".format(self.q)) +
                                  ("c1 = {};\n".format(self.c1)) + ("c2 = {};".format(self.c2)))

    def guess_coefficients(self, points):
        # If we don't have either gamma0 or rho0, just up the first point for both.
        if (self.gamma0 == None and self.rho0 == None):
            self.rho0 = points[0][0]
            self.gamma0 = points[0][1]
        elif (self.gamma0 == None):  # If we have rho0 but no gamma, find gamma0 at that point
            x_values = [point[0] for point in points]
            y_values = [point[1] for point in points]
            interpFunc = interp.interp1d(x_values, y_values);
            self.gamma0 = interpFunc(rh0)
        elif (self.rho0 == None):
            raise RuntimeError("Please define rho0 if gamma0 is defined")

        if (self.q == None):  # Calculate a guess for q
            yy = (self.gamma0 / self.c2) - (1 / (2 * self.c2)) - (self.c1 * pow(self.rho0, 1.0 / 3.0)) / self.c2
            self.q = - math.log(yy) / math.log(self.rho0)

    @staticmethod
    def _f(x, *coeffs):
        (gamma0, rho0, q, c1, c2) = coeffs
        term1 = c1 * pow(x, 1.0 / 3.0)
        term2 = c2 * pow(x, q)
        return 0.5 + term1 + term2

    def _derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _second_derivative(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _energy_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")

    def _pressure_integral(self, x):
        raise RuntimeError("NOT YET IMPLEMENTED")


factory.register('gammaV', GammaV)


class GammaPoly(PolyBase):
    """
    Fits gamma as a function of density (rho) or volume
    Separate polynomials are fitted for high pressure versus low pressure data (relative to some reference rho0)
        High pressure means high density or low volume
        Low pressure means low density or high volume
    """
    _name_prefix = 'gammapoly'

    # TODO: wire in rho_is_density so external app can use.  Maybe just subclass GammaV and pass False to __init__?
    def __init__(self, args, name, rho_is_density=True):
        """

        :param args:
        :param name:
        :param rho_is_density: If True, x values and rho0 are densities.  If False, x values and rho0 are unit volumes.
        """
        self._order = int(name[len(self.name_prefix):])
        self._is_first_fit = True
        self.rho0 = args.rho0
        self._overlap = args.overlap
        # create two fitters: one for high P data and one for low P data
        fname = '{0}{1}'.format(Poly_Original.name_prefix, self._order)
        self._hiP_fitter = Poly_Original(args, fname)
        self._loP_fitter = Poly_Original(args, fname)
        self._highP_f = None
        self._lowP_f = None
        # rho_is_density indicates whether x values and rho0 are density (True) or volume (False)
        self._rho_is_density = rho_is_density

    def _get_highP_lowP_x_indices(self, x):
        """
        Separates independent variable into high pressure and low pressure regimes
        If self.rho0 is present, it will be in both
        """
        # determine the correct comparison function for (x, rho0) depending on whether using density or volume
        hiP_comp = np.greater_equal if self._rho_is_density else np.less_equal
        loP_comp = np.less_equal if self._rho_is_density else np.greater_equal
        hiP_indices = np.flatnonzero(hiP_comp(x, self.rho0))
        loP_indices = np.flatnonzero(loP_comp(x, self.rho0))
        return hiP_indices, loP_indices

    def _get_highP_fit_x(self, x):
        """ for high pressure points, fit a polynomial in rho/rho0 if rho_is_density, or V0/V (the reciprocal) otherwise
        """
        x1_x2 = (x, self.rho0) if self._rho_is_density else (self.rho0, x)
        return np.divide(*x1_x2)  # use numpy to correctly handle divide by 0

    def _get_lowP_fit_x(self, x):
        """ for low pressure points, fit a polynomial in rho0/rho if rho_is_density, or V/V0 (the reciprocal) otherwise
        """
        # use np.divide to handle divide by zero
        x1_x2 = (self.rho0, x) if self._rho_is_density else (x, self.rho0)
        return np.divide(*x1_x2)  # use numpy to correctly handle divide by 0

    def _get_highP_lowP_fit_points(self, points):
        """ Gets high pressure and low pressure points for fitting
        x values in points are converted to the appropriate expression depending on whether x is density or volume
        Includes some overlap with the other domain as dictated by the overlap parameter

        :param points: numpy array of (x,y) points to fit, sorted by x
                    and with x in terms of density or volume as indicated by the rho_is_density arg
        :return: Two arrays of (fit_x, y) 2-tuples, where fit_x is in terms of some ratio with the reference density
                    the first array gives points for the high pressure fit
                    the second array gives points for the low pressure fit
        """

        def get_overlap_point_count(pt_count):
            return self._overlap if pt_count > self._overlap else pt_count

        def add_points_to_end(orig_points_idx, points_to_add_idx, point_count):
            # append the first point_count points from points_to_add to the end of original_points
            return np.append(orig_points_idx, points_to_add_idx[0:point_count]) if point_count else orig_points_idx

        def add_points_to_beginning(orig_points_idx, points_to_add_idx, point_count):
            # insert the first point_count points from points_to_add at the beginning of original_points
            return np.insert(orig_points_idx, 0,
                             points_to_add_idx[-1 * point_count:]) if point_count else orig_points_idx

        points = np.array(points)  # wrap for easier indexing; if already a numpy.ndarray, this won't do anything
        hiP_idx, loP_idx = self._get_highP_lowP_x_indices([x for (x, y) in points])
        # add overlap points to non-empty domains
        add_hi2lo_ct = get_overlap_point_count(len(hiP_idx)) if len(loP_idx) else 0
        add_lo2hi_ct = get_overlap_point_count(len(loP_idx)) if len(hiP_idx) else 0
        if self._rho_is_density:  # x values are in density (high pressure means high density)
            upd_hiP_idx = add_points_to_beginning(hiP_idx, loP_idx, add_lo2hi_ct)
            upd_loP_idx = add_points_to_end(loP_idx, hiP_idx,
                                            add_hi2lo_ct)  # must make second object since hiP_idx changed
        else:  # x values are in unit volume (high pressure means low unit volume)
            upd_hiP_idx = add_points_to_end(hiP_idx, loP_idx, add_lo2hi_ct)
            upd_loP_idx = add_points_to_beginning(loP_idx, hiP_idx, add_hi2lo_ct)
        hiP_fit_pts = np.array([(self._get_highP_fit_x(x), y) for (x, y) in points[upd_hiP_idx]])
        loP_fit_pts = np.array([(self._get_lowP_fit_x(x), y) for (x, y) in points[upd_loP_idx]])
        return hiP_fit_pts, loP_fit_pts

    def guess_coefficients(self, points):
        pass

    def fit_to_points(self, points):
        def _get_fit_description(is_hiP):
            out = '{0} pressure fit, where x = {1}'
            if self._rho_is_density:
                x = rho + '/' + rho + naught if is_hiP else rho + naught + '/' + rho
            else:
                x = 'V' + naught + '/V' if is_hiP else 'V/V' + naught
            return out.format('High' if is_hiP else 'Low', x)

        highP_points, lowP_points = self._get_highP_lowP_fit_points(points)
        print(_get_fit_description(True))
        self._hiP_fitter.fit_to_points(highP_points)
        print(_get_fit_description(False))
        self._loP_fitter.fit_to_points(lowP_points)
        self._is_first_fit = False

    def _eval(self, x, hiP_fn, loP_fn):
        x = np.asarray([x] if np.isscalar(x) else x)
        hiP_idx, loP_idx = self._get_highP_lowP_x_indices(x)
        y = np.empty(x.shape, np.float)
        y[hiP_idx] = hiP_fn(self._get_highP_fit_x(x[hiP_idx]))
        y[loP_idx] = loP_fn(self._get_lowP_fit_x(x[loP_idx]))
        return y

    def func(self, x):
        # TODO: if rho0 is among values of x, will return lowP fit.  Make this configurable?
        return self._eval(x, self._hiP_fitter.func, self._loP_fitter.func)

    def derivative(self, x):
        return self._eval(x, self._hiP_fitter.derivative, self._loP_fitter.derivative)

    def second_derivative(self, x):
        return self._eval(x, self._hiP_fitter.second_derivative, self._loP_fitter.second_derivative)

    def integral(self, x):
        return self._eval(x, self._hiP_fitter.integral, self._loP_fitter.integral)


factory.register(GammaPoly.name_prefix, GammaPoly)


class GammaPolyV(GammaPoly):
    _name_prefix = 'gammapolyv'

    def __init__(self, args, name):
        super().__init__(args, name, rho_is_density=False)


factory.register(GammaPolyV.name_prefix, GammaPolyV)

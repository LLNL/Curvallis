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

from __future__ import print_function
import copy

import curve_fitters
import io
import lines
import smoothers 
from pylab import polyfit
import numpy as np
from operator import itemgetter
from math import log

_INFINITY = float ('inf')

# TODO: recalculate the equation for a region when a point is moved in that region and then plot the region

def define_args(parser):
    parser.add_argument(
        '--do_derivative',
        action='store_true',
        help='Calculate and plot the derivative of the fit function '
             '[default: %(default)s]')
    parser.add_argument(
        '--do_integral',
        action='store_true',
        help='Calculate and plot the integral of the fit function '
             '[default: %(default)s]')
    parser.add_argument(
        '--points_in_fit_curve', action='store', type=int,
        help='Calculate this many points in the fit curve when writing the '
             'curve to a file [default: %(default)s]', metavar='<count>')
    parser.add_argument(
        '--points_in_user_curve', action='store', type=int,
        help='Calculate this many points in the user inputted curves when writing the '
             'curves to a file [default: %(default)s]', metavar='<count>')
    parser.add_argument(
        '--region_bound',
        action='append',
        nargs='+', type=float,
        help='Select where region boundaries are supposed to go in sequential order'
             '[default: Evenly Spaced]',
        metavar='<bound>')
    parser.add_argument(
        '--overlap',
        type=int,
        help='Select number of points to overlap for region fitted curves. '
             '[default: %(default)s]',
        metavar='<int>')
    parser.add_argument(
        '--numpoints',
        type=int,
        help='Select the number of points for the local average trilocal smoothing.'
             '[default: %(default)s]',
        metavar='<int>')
    parser.add_argument(
        '--repeat',
        type=int,
        help='Select the number of times to repeat trilocal smoothing.'
             '[default: %(default)s]',
        metavar='<int>')
    parser.add_argument(
        '--matchpt',
        type=float,
        help='Select the matchpoint for integral and tri-integral smoothing.'
             '[default: %(default)s]',
        metavar='<float>')
    parser.add_argument(
        '--interp',
        help='Select the interpolator for integral and tri-integral smoothing.'
             '[default: %(default)s]',
        metavar='<arg>')
    parser.add_argument(
        '--angle',
        type=int,
        help='Select the angle for acute angle repair.'
             '[default: %(default)s]',
        metavar='<int>')

    parser.set_defaults(
        do_derivative=False,
        do_integral=False,
        points_in_fit_curve=100,
        points_in_user_curve=100,
        polynomial_order=5,
        overlap=2,
        numpoints=5,
        repeat=10,
        matchpt=-1.0,
        angle=90,
        interp="cubic",
    )


class _Line_Set_With_Fit(lines.Line_Set):
    """ Has own fitter object, plus various fit curve lines
    """
    def __init__(self, name, ax, fitter,fitter2,  ghost_set, x_low_limit, 
                 x_high_limit, allow_xy_move, args, is_eos_data):
        super(_Line_Set_With_Fit, self).__init__(name, ax, allow_xy_move, is_eos_data)
        self._args = args
        self._fitter = copy.deepcopy(fitter)
        self._fitter2 = copy.deepcopy(fitter2)
        self._ghost_set = ghost_set
        self._update_points = False
        self._x_high_limit = x_high_limit
        self._x_low_limit = x_low_limit
        self._x_view_high_limit = x_high_limit
        self._x_view_low_limit = x_low_limit
        self._logscale = False
        self._ax = ax
        if not is_eos_data:
            self.fit_curve = lines.Line(ax, lines.line_attributes['fit_curve'])
            self.derivative_curve = lines.Line(ax, lines.line_attributes['derivative'])
            self.integral_curve = lines.Line(ax, lines.line_attributes['integral'])

    def _calc_x_values(self, x_first, x_last, x_count):
        """ For a given x range, interpolate x_count x values (without any
          cumulative rounding errors) and return them.

        :param x_first:
        :param x_last:
        :param x_count:
        :return: list
        """
        if self._logscale == True:
            x_first = log(x_first,10)
            x_last = log(x_last,10)

        result = []
        x_range = x_last - x_first
        for i in range(x_count):
            # Calculate each x without any cumulative errors:
            portion = float(i) / float(x_count-1)
            x = x_first + (portion * x_range)

            if self._logscale == True:
                result.append(pow(10,x))
            else:
                result.append(x)
        return result

    def calc_fit_points_for_range(self, x_first, x_last, point_count):
        """ For a given x range, interpolate x_count x values and calculate
        a y value for each, returning the calculated x and y values.

        """
        x_values = self._calc_x_values(x_first, x_last, point_count)
        y_values = []

        if (self._fitter2 == 'none'):
            for x in x_values:
                y_values.append(self._fitter.func(x))
        else:
            for x in x_values:
                y_values.append(self._fitter.func(x) - self._fitter2.func(x))

        return zip(x_values, y_values)

    def _calc_derivative_points_for_range(self, x_first, x_last, point_count):
        """ For a given x range, interpolate x_count x values and calculate
        a y value for each, returning the calculated x and y values.

        """
        x_values = self._calc_x_values(x_first, x_last, point_count)
        y_values = []

        if (self._fitter2 == 'none'):
            for x in x_values:
                y_values.append(self._fitter.derivative(x))
        else:
            for x in x_values:
                y_values.append(self._fitter.derivative(x) - self._fitter2.derivative(x))

        return zip(x_values, y_values)

    def _calc_integral_points_for_range(self, x_first, x_last, point_count):
        """ For a given x range, interpolate x_count x values and calculate
        a y value for each, returning the calculated x and y values.

        """
        x_values = self._calc_x_values(x_first, x_last, point_count)
        y_values = []

        if (self._fitter2 == 'none'):
            for x in x_values:
                y_values.append(self._fitter.integral(x))
        else:
            for x in x_values:
                y_values.append(self._fitter.integral(x) - self._fitter2.integral(x))

        return zip(x_values, y_values)


    def get_ghost_points(self):
        return self._ghost_set

    def calculate_fit(self):
        """ Derive the coefficients for the function for this region's subset of
        points

        Call at startup AFTER creating self.lines.movable.id, and whenever a
        point is moved.
        """
        if not self._is_eos_data and self._fitter != 'none':
            points = self.get_ghost_points()

            # guess better coefficients than defaults for first fit
            if self._fitter._is_first_fit:
                self._fitter.guess_coefficients(points)

            # print ('points:\n%s' % points)
            point_count = self.get_movable_point_count()
            if point_count >= 2:
                self._fitter.fit_to_points(points)

                # Run if combined fit type
                if (self._fitter2 != 'none'):
                    # Calculate error from first fit to use for refined fit
                    orig_points, y_vals = zip(*points)
                    err_points = []
                    for i in range(0, point_count):
                        err_points.append(self._fitter.func(orig_points[i]) - y_vals[i])
                        points2 = zip(orig_points, err_points)
                    self._fitter2.fit_to_points(points2)
            else:
                print('!!! Not doing fit for line set %s because num points = %s' %
                      (self._name, point_count))

    def get_movable_point_count(self):
        """
        :return: int
        """
        return len(self.get_movable_points())

    def plot_curves(self):
        """ Plot a smooth fit curve for the current zoom x range.

        Call after calling calculate_fit, and whenever xlim changes.
        """
        if not self._is_eos_data and self.get_movable_point_count() >=2 and self._fitter != 'none':
            # calculate new limits of view
            x_limit = self._ax.get_xlim()

            # Keep line drawn within region boundaries
            if self._x_low_limit < x_limit[0]:
                self._x_view_low_limit = x_limit[0]
            else:
                self._x_view_low_limit = self._x_low_limit
            if self._x_high_limit > x_limit[1]:
                self._x_view_high_limit = x_limit[1]
            else:
                self._x_view_high_limit = self._x_high_limit

            fit_points = self.calc_fit_points_for_range(
                x_first=self._x_view_low_limit,
                x_last=self._x_view_high_limit,
                point_count=self._args.points_in_fit_curve)

            # Need animated = True for curve plots to get the last curve to go away:
            if self.fit_curve._id == None:
                self.fit_curve.plot_xy_data(fit_points,
                                            animated=True)
            else:
                self.fit_curve.set_xy_data(fit_points)
            if self._args.do_derivative:
                derivative_points = self._calc_derivative_points_for_range(
                    x_first=self._x_low_limit,
                    x_last=self._x_high_limit,
                    point_count=self._args.points_in_fit_curve)
                self.derivative_curve.plot_xy_data(derivative_points,
                                                   animated=True)
            if self._args.do_integral:
                integral_points = self._calc_integral_points_for_range(
                    x_first=self._x_low_limit,
                    x_last=self._x_high_limit,
                    point_count=self._args.points_in_fit_curve)
                self.integral_curve.plot_xy_data(integral_points,
                                                 animated=True)

    def undo(self):
        self.movable.undo()
        # Undo highlighting if in block select mode
        # Must do here to check if block select mode is active
        if (self._in_set == True and self.movable._last_highlight != None 
            and self.movable._last_highlight != []):
            self.movable._highlight.set_data(zip(*self.movable._last_highlight))
        self.movable.draw()

        # Undo information for any fit curves as well
        if not self._is_eos_data:
            self.fit_curve.undo()
            self.derivative_curve.undo()
            self.integral_curve.undo()
        

# Point movement ###############################################################

    def finish_move_point(self):
        """ Recalculate and redraw the fit curve.
        """
        super(_Line_Set_With_Fit, self).finish_move_point()
        self.calculate_fit()
        self.plot_curves()

    def check_move_point(self):
        if (self._update_points == True):
            self.calculate_fit()
            self.plot_curves()
            self._update_points = False

# END Point movement ###########################################################


class _Line_Sets(object):
    """
    Similar to io.Data_Sets.  One Line_Set object per named Data_Set.
    """
    def __init__(self, names, ax, fitter, fitter2, ghost_sets, x_low_limit, 
                 x_high_limit, allow_xy_move, args, is_eos_data):
        self._args = args
        self._ax = ax
        self._is_eos_data = is_eos_data
        self._fitter = copy.deepcopy(fitter)
        self._fitter2 = copy.deepcopy(fitter2)
        self._ghost_sets = ghost_sets
        self._moving_point_set_name = None
        self._sets = {}
        self._region_boundary = lines.Line(ax, lines.line_attributes['region_boundary'])
        self._x_high_limit = x_high_limit
        self._x_low_limit = x_low_limit

        for name in names:
            self._sets[name] = _Line_Set_With_Fit(
                ax=self._ax,
                fitter=self._fitter,
                fitter2=self._fitter2,
                ghost_set = self._ghost_sets._sets[name],
                name=name,
                x_low_limit=self._x_low_limit,
                x_high_limit=self._x_high_limit,
                allow_xy_move=allow_xy_move,
                args=self._args,
                is_eos_data=self._is_eos_data)

    def calculate_fit(self):
        for line_set in self._sets.itervalues():
            line_set.calculate_fit()

    def draw(self):
        self._region_boundary.draw()
        for lines in self._sets.itervalues():
            lines.draw()

    def get_data_sets(self):
        """ Get the data set from each line set and return them.  Each name
        should appear only once.

        :return: Data_Sets
        """
        result = io.Data_Sets()
        for line_set in self._sets.itervalues():
            result.add_set(name=line_set.get_name(),
                           points=line_set.get_movable_points())
        return result

    def get_a_line_set(self):
        """ For when just any line set will do:

        :return: Line_Set
        """
        return self._sets.values()[0]

    def get_fit_curve_points(self):
        """ Return the line set's ONLY fit curve's points.
        """
        line_set = self._get_only_line_set()
        # The current fit curve data is custom-generated for the current zoom's .
        # x range. Calculate the curve for the the movable line's full x range:
        return line_set.calc_fit_points_for_range(
            x_first=self._x_low_limit,
            x_last=self._x_high_limit,
            point_count=self._args.points_in_fit_curve)

    def get_info(self, indent=''):
        result  = indent + 'Line sets count: %s\n' % len(self._sets)
        indent += '    '
        for name, line_set in self._sets.iteritems():
            result += indent + 'Line set "%s":' % name
            result += line_set.get_info(indent + '    ')
        return result

    def _get_only_line_set(self):
        """ For when there should be only one line set (non-EOS data):
        :return: list of (x,y) points
        """
        assert len(self._sets) == 1
        return self._sets.values()[0]

    def len(self):
        return len(self._sets)

    def plot_boundary_line(self):
        min_y, max_y = self._ax.get_ylim()
        x = self._x_high_limit
        self._region_boundary.plot_xy_data(((x, min_y), (x, max_y)))

    def plot_curves(self):
        for line_set in self._sets.itervalues():
            line_set.plot_curves()

    def plot_movable_xy_data(self, data_sets, visible=True, animated=True):
        for name, data_set in data_sets.iteritems():
            self._sets[name].movable.plot_xy_data(data_set, visible, animated)

    def plot_original_xy_data(self, data_sets, visible=False, animated=True):
        for name, data_set in data_sets.iteritems():
            self._sets[name].original.plot_xy_data(data_set, visible, animated)

    def set_allow_xy_move(self, allow):
        for line_set in self._sets.itervalues():
            line_set.set_allow_xy_move(allow)

    def toggle_original_line_visibility(self):
        for line_set in self._sets.itervalues():
            line_set.toggle_original_line_visibility()

    def replot(self):
        for line_set in self._sets.itervalues():
            line_set.replot()

    def update_ghost_points(self, newlist):
        for line_set in self._sets.itervalues():
            oldlist = line_set._ghost_set
            if (not np.array_equal(oldlist, newlist)):
                line_set._ghost_set = newlist
                line_set._update_points = True
            else:
                line_set._update_points = False

    def find_set_span(self):
        xmin = None
        xmax = None
        for line_set in self._sets.itervalues():
            if line_set._in_set == True:
                if line_set._set_points[0][0] < xmin or xmin == None:
                    xmin = line_set._set_points[0][0]
                if line_set._set_points[-1][0] > xmax or xmax == None:
                    xmax = line_set._set_points[-1][0]
        return xmin, xmax

    def undo(self):
        for line_set in self._sets.itervalues():
            line_set.undo()

# Point movement ###############################################################

    def attempt_begin_move_point(self, event):
        """ Select the first line set where the cursor is within EPSILON of one
          of its points.
        """
        for name, line_set in self._sets.iteritems():
            line_set.attempt_begin_move_point(event)
            if line_set.move_point_in_progress():
                self._moving_point_set_name = name
                # Reset undo information for all other lines
                for name, line_set in self._sets.iteritems():
                    if name != self._moving_point_set_name:
                        line_set.movable._last_data = None
                break

    def attempt_get_set(self, xmin, xmax, ymin, ymax):
        for name, line_set in self._sets.iteritems():
            line_set.attempt_get_set(xmin, xmax, ymin, ymax)

    def begin_move_set(self, event):
        for name, line_set in self._sets.iteritems():
            if line_set._in_set == True:
                line_set.begin_move_set(event)

    def finish_move_set(self):
        for name, line_set in self._sets.iteritems():
            if line_set._in_set == True:
                line_set.finish_move_set()

    def cancel_move_set(self):
        for name, line_set in self._sets.iteritems():
            line_set.cancel_move_set()

    def move_set(self, event):
        for name, line_set in self._sets.iteritems():
            if line_set._in_set == True:
                line_set.move_set(event)

    def rotate_set(self, event, xmin, xmax, ymin, ymax):
        for name, line_set in self._sets.iteritems():
            if line_set._in_set == True:
                line_set.rotate_set(event, xmin, xmax, ymin, ymax)

    def cancel_any_move_points(self):
        self._moving_point_set_name = None
        for line_set in self._sets.itervalues():
            line_set.cancel_any_move_points()

    def finish_move_point(self):
        if self.move_point_in_progress():
            self._sets[self._moving_point_set_name].finish_move_point()

    def check_move_point(self):
        for name, line_set in self._sets.iteritems():
            line_set.check_move_point()

    def move_point(self, event):
        if self.move_point_in_progress():
            self._sets[self._moving_point_set_name].move_point(event)

    def moving_point_info(self, xy_only=False):
        if self.move_point_in_progress():
            return self._sets[self._moving_point_set_name].moving_point_info(xy_only)
        else:
            return '(no moving point)'

    def move_point_in_progress(self):
        return self._moving_point_set_name is not None

    def add_point(self, event):
        # FIX: Will add point to every line if more than one
        # Currently disabled for 2d plots
        for name, line_set in self._sets.iteritems():
            line_set.add_point(event)

    def remove_point(self, event):
        # FIX: Will try to remove point from every line if more than one
        # Currently disabled for 2d plots
        for name, line_set in self._sets.iteritems():
            line_set.remove_point(event)

# END Point movement ###########################################################

class _Region(object):
    """ A region is a subset in the X axis (a vertical slice) of the plotting
    area, with a separately calculated fit curve.
    """
    def __init__(self, id_in, ax, data_sets, ghost_sets, x_low_limit, x_high_limit, fitter,
                 fitter2, is_eos_data, allow_xy_move, args):
        """
        :param id_in:   Region number
        :param ax:
        :param data_sets:  Data_Sets, with points between x_low_limit and x_high_limit
        :param x_low_limit:
        :param x_high_limit:
        :param fitter:
        :param is_eos_data: Boolean

        """
        self._args=args
        self._ax = ax
        self._is_eos_data = is_eos_data
        self._ghost_sets = ghost_sets

        if fitter != 'none':
            self._fitter = curve_fitters.factory.make_object_of_class(
            name=fitter, args=args)
        else:
            self._fitter = 'none'
        if fitter2 != 'none':
            self._fitter2 = curve_fitters.factory.make_object_of_class(
            name=fitter2, args=args)
        else:
            self._fitter2 = 'none'

        self._id = id_in
        self._moving_point_data_set_name = None
        self._moving_set_min = None
        self._moving_set_max = None
        self._x_high_limit = x_high_limit
        self._x_low_limit = x_low_limit

        self._line_sets = _Line_Sets(
                names=data_sets.get_names(),
                ax=ax,
                fitter=self._fitter,
                fitter2=self._fitter2,
                ghost_sets=self._ghost_sets,
                x_low_limit=x_low_limit,
                x_high_limit=x_high_limit,
                allow_xy_move=allow_xy_move,
                args=self._args,
                is_eos_data=self._is_eos_data)
        self._line_sets.plot_original_xy_data(data_sets)
        self._line_sets.plot_movable_xy_data(data_sets)
        print (self.get_info())
        self._line_sets.plot_boundary_line()
        print("Approximate Fit (No ghost points yet):")
        self.calculate_fit()
        self._line_sets.plot_curves()

    def calculate_fit(self):
        self._line_sets.calculate_fit()

    def display_point_is_below_region(self, x, y):
        data_x = self._line_sets.get_a_line_set().movable.point_to_data_space(x, y)[0]
        return data_x < self._x_low_limit

    def display_point_is_in_region(self, x, y):
        data_x = self._line_sets.get_a_line_set().movable.point_to_data_space(x, y)[0]
        return self._x_low_limit <= data_x <= self._x_high_limit

    def get_data_sets(self):
        """ Get the data set from each line set and return them.

        :return: Data_Sets
        """
        return self._line_sets.get_data_sets()

    def get_info(self, indent=''):
        result  = indent + 'REGION "%s" INFO:\n' % self._id
        indent += '    '
        result += indent + 'X Range: %.15E .. %.15E\n' % (self._x_low_limit, self._x_high_limit)
        result += self._line_sets.get_info(indent)
        # result += 'Points:\n%s' % self._sets.movable.get_xy_data()
        return result

    def get_x_low_limit(self):
        return self._x_low_limit

    def get_x_high_limit(self):
        return self._x_high_limit

    def draw(self):
        self._line_sets.draw()

    def get_fit_curve_points(self):
        """ Return the region's ONLY fit curve's points.
        """
        return self._line_sets.get_fit_curve_points()

    def plot_curves(self):
        self._line_sets.plot_curves()

    def set_allow_xy_move(self, allow):
        self._line_sets.set_allow_xy_move(allow)

    def toggle_original_line_visibility(self):
        self._line_sets.toggle_original_line_visibility()

    def replot(self):
        self._line_sets.replot()

    def undo(self):
        self._line_sets.undo()

    # Point movement ###########################################################

    def attempt_begin_move_point(self, event):
        self._line_sets.attempt_begin_move_point(event)

    def attempt_get_set(self, xmin, xmax, ymin, ymax):
        self._line_sets.attempt_get_set(xmin, xmax, ymin, ymax)

    def begin_move_set(self, event):
        #Set edges to keep points in region
        xmin, xmax = self._line_sets.find_set_span()
        if (xmin != None and xmax != None):
            line_set = self._line_sets._sets.itervalues().next()
            xmin_disp = line_set.movable.points_to_display_space([[xmin, 0]])[0][0]
            xmax_disp = line_set.movable.points_to_display_space([[xmax, 0]])[0][0]

            self._moving_set_min = event.x - xmin_disp
            self._moving_set_max = xmax_disp - event.x

        self._line_sets.begin_move_set(event)

    def finish_move_set(self):
        self._line_sets.finish_move_set()

    def cancel_move_set(self):
        self._line_sets.cancel_move_set()

    def move_set(self, event):
        """ Move set of points in the data and replot
        """
        def keep_x_in_region():
            line_set = self._line_sets._sets.itervalues().next()

            lowest = event.x - self._moving_set_min
            highest = event.x + self._moving_set_max
            
            low_limit_disp = line_set.movable.points_to_display_space([[self._x_low_limit, 0]])[0][0]
            high_limit_disp = line_set.movable.points_to_display_space([[self._x_high_limit, 0]])[0][0]

            if lowest < low_limit_disp:
                event.x = low_limit_disp + self._moving_set_min
            if highest > high_limit_disp:
                event.x = high_limit_disp - self._moving_set_max

        keep_x_in_region()
        self._line_sets.move_set(event)

    def rotate_set(self, event, xmin, xmax, ymin, ymax):
        #TODO: MAKE SURE SET STAYS IN REGION
        self._line_sets.rotate_set(event, xmin, xmax, ymin, ymax)

    def cancel_any_move_points(self):
        self._line_sets.cancel_any_move_points()

    def finish_move_point(self):
        self._line_sets.finish_move_point()

    def check_move_point(self):
        self._line_sets.check_move_point()

    def move_point(self, event):
        """ Move the point in the data and replot the line.
        """
        def keep_x_in_region():
            event.xdata = max(event.xdata, self._x_low_limit)
            event.xdata = min(self._x_high_limit, event.xdata)

        if self.move_point_in_progress():
            keep_x_in_region()
            self._line_sets.move_point(event)

    def moving_point_info(self, xy_only=False):
        if self.move_point_in_progress():
            point_info = self._line_sets.moving_point_info(xy_only)
            if xy_only:
                return point_info
            else:
                return 'region %s, %s' % (self._id, point_info)
        else:
            return '(no moving point)'

    def move_point_in_progress(self):
        return  self._line_sets.move_point_in_progress()

    def add_point(self, event):
        self._line_sets.add_point(event)

    def remove_point(self, event):
        self._line_sets.remove_point(event)

    # END Point movement #######################################################


class Regions(object):
    """ Manages all the regions in the plotting area.
    """
    def __init__(self, ax, args, input_data_sets, xy_limits, io_manager):
        """ Create the regions, their boundaries, their fit curve funcs,
        their fit curve plots, etc.
        The regions are ordered.  [0] is the first region, with the lowest x values.

        :param ax:
        :param args::
        :param input_data_sets: io.Data_Sets
        :param xy_limits.:      io.XY_Limits
        :param io_manager:      io.Manager
        """
        self._args = args
        self._ax = ax
        self._io_manager = io_manager
        self._is_eos_data = self._io_manager.is_eos_data()
        self._allow_xy_move = not self._is_eos_data
        self._moving_point_region_index = None
        self._moving_set_region_index = None
        self._data_sets = input_data_sets.get_copy()
        self._lowest_boundary_line = lines.Line(ax, lines.line_attributes['region_boundary'])
        self._regions=list()
        self._xlim_callback_active = False
        self._x_max = xy_limits.x_max
        self._x_min = xy_limits.x_min
        self._fitter = args.fit_type
        self._fitter2 = args.refine_fit
        if self._args.region_bound != None:
            region_count = len(self._args.region_bound[0]) + 1
        else:
            region_count = 1
        self._create_regions(region_count)
        print ("regions:")
        self.calculate_fits()
#        print ("regions:\n%s" % pprint.pformat(self.regions))
        self.plot_curves()

    def calculate_fits(self):
        """ Derives a fit function for each region's subset of points
        """
        for region in self._regions:
            region.calculate_fit()

    def _create_regions(self, region_count):
        """ Split self._data_sets up into region_count X-slices and create a
        region for each slice in self._regions.  Divide up the regions up evenly
        over the X range.
        """

        def ascending(mylist):
            for i in range (0, len(mylist)-1):
                if float(mylist[i]) > float(mylist[i+1]):
                    return False
            return True

        def calculate_x_boundaries():
            # noinspection PyUnusedLocal
            result = [_INFINITY for unused in range(region_count + 1)]
            result[0] = self._x_min
            result[-1] = self._x_max
            x_range = self._x_max - self._x_min
            for boundary_index in range(1, region_count):
                # Calculate the fraction each time to prevent rounding error:
                x_range_fraction = float(boundary_index) / float(region_count)
                result[boundary_index] = self._x_min + x_range * x_range_fraction
            return result

        def input_x_boundaries():
            # Check boundaries were entered in ascending order
            assert ascending(self._args.region_bound[0]), "Region boundaries must be entered in ascending order."

            result = [_INFINITY for unused in range(region_count + 1)]
            result[0] = self._x_min
            result[-1] = self._x_max
            for i in range(1, region_count):
                result[i] = float(self._args.region_bound[0][i-1])
                assert self._x_min < result[i] < self._x_max, "Region boundaries are not in range of the data. (%E to %E)" % (self._x_min, self._x_max)

            return result

        if self._args.region_bound != None:
            x_boundaries = input_x_boundaries()
        else:
            x_boundaries = calculate_x_boundaries()
        self._plot_lowest_boundary_line(x_boundaries[0])
        for region_index in range(region_count):
            x_low = x_boundaries[region_index]
            x_high = x_boundaries[region_index + 1]
            print ('Creating region %s for X range %.15E .. %.15E' %
                   (region_index, x_low, x_high))
            if len(self._fitter) > region_index:
                myfitter = copy.copy(self._fitter[region_index])
            else:
                #Default value for fitter
                myfitter = 'poly5'
            if len(self._fitter2) > region_index:
                myfitter2 = copy.copy(self._fitter2[region_index])
            else:
                #Default value for refine_fit (fitter2)
                myfitter2 = 'none'
            self._regions.append(
                _Region(
                    id_in=region_index,
                    ax=self._ax,
                    data_sets=self._data_sets.get_x_slice(x_low, x_high),
                    ghost_sets=self._data_sets.get_x_slice(x_low, x_high),
                    x_low_limit=x_low,
                    x_high_limit=x_high,
                    fitter=myfitter,
                    fitter2=myfitter2,
                    is_eos_data=self._is_eos_data,
                    allow_xy_move = self._allow_xy_move,
                    args=self._args))
        self._update_ghost_points()

    def draw(self):
        self._lowest_boundary_line.draw()
        for region in self._regions:
            region.draw()

    def _get_fit_curve_points(self):
        """ Return a list of all the regions' ONLY fit curve points, concatenated.
        """
        assert not self._is_eos_data
        points=[]
        for region in self._regions:
            points.extend(region.get_fit_curve_points())

        return points

    def _get_data_sets(self):
        """ Get the data sets from every region, glue them back together as if
        there were only one region, and return that.

        :return: Data_Sets
        """
        result = io.Data_Sets()
        for region in self._regions:
            region_data_sets = region.get_data_sets()
            for name, region_data_set in region_data_sets.iteritems():
                result.add_or_append_to_set(name, region_data_set)
        return result


    def get_region_index(self, display_x, display_y):
        """ Given a point in display space, return the index of the region it
         is in.

        Instead of raising an exception if x is out of range, just return the
        highest or lowest region.

        :param display_x: (int) x value in display space
        :param display_y: (int) y value in display space
        :return:          (int) index into self.regions
        """
        for index in range(len(self._regions)):
            if self._regions[index].display_point_is_in_region(display_x, display_y):
                return index
        if self._regions[0].display_point_is_below_region(display_x, display_y):
            return 0
        else:
            return len(self._regions) - 1

    def plot_curves(self):
        """ Plots a fit curve for each region. Called by xlim_callback
        """
        # The draw op in plot_fit_curve may call matplotlib.Axes.set_xscale ->
        # autoscale_view -> set_xbound -> set_xlim -> self.xlim_changed_callback
        # -> self.plot_fit_curves
        # SO: disable the set_xlim callback temporarily so we don't have
        # infinite recursion:
        if not self._xlim_callback_active:
            self._xlim_callback_active = True
            for region in self._regions:
                region.plot_curves()
            self._xlim_callback_active = False

    def _plot_lowest_boundary_line(self, x):
        min_y, max_y = self._ax.get_ylim()
        self._lowest_boundary_line.plot_xy_data(((x, min_y), (x, max_y)))

    def toggle_allow_xy_move(self):
        self._allow_xy_move = not self._allow_xy_move
        for region in self._regions:
            region.set_allow_xy_move(self._allow_xy_move)
        print ("Allow xy move: %s" % self._allow_xy_move)

    def toggle_original_line_visibility(self):
        for region in self._regions:
            region.toggle_original_line_visibility()

    def toggle_points(self):
        """
        Disable point markers on movable lines.
        """
        if (lines.line_attributes['movable']['marker'] != None or
            lines.line_attributes['movable_no_curves']['marker'] != None):
            lines.line_attributes['movable']['marker'] = None
            lines.line_attributes['movable_no_curves']['marker'] = None
        else:
            lines.line_attributes['movable']['marker'] = 'o'
            lines.line_attributes['movable_no_curves']['marker'] = 'o'
            
        for region in self._regions:
            region.replot()

    def write_output_files(self):
        data_sets = self._get_data_sets()
        self._io_manager.write_movable_data_sets(data_sets)
        if not self._is_eos_data:
            fit_points = self._get_fit_curve_points()
            filtered_points = [fit_points[0]] #don't miss the 1st point

            # Remove any duplicate fit points at region boundaries
            for i in range(1, len(fit_points)-2):
                if fit_points[i][0] != fit_points[i-1][0]:
                    filtered_points.append(fit_points[i])
                    
            io.write_point_file(self._args.curve_output_file_name,
                                filtered_points)

    def smooth_data(self, smooth_type, xmin, xmax, ymin, ymax):
        """
        Apply smoothing algorithm to all movable data lines.
        """
        print ("Smoothing...")
        for region in self._regions:
            for line_set in region._line_sets._sets.itervalues():
                orig_line = line_set.movable.get_x_data_y_data()
                line = []
                line.append([])
                line.append([])
                #Only alter points in region
                for i in range(0, len(orig_line[0])):
                    #Find start of region
                    if ((xmin < orig_line[0][i] < xmax) and (ymin < orig_line[1][i] < ymax)):
                        start = i
                        break
                for i in range(0, len(orig_line[0])):
                    #Find end of region
                    if ((xmin < orig_line[0][i] < xmax) and (ymin < orig_line[1][i] < ymax)):
                        end = i
                for i in range(0, len(orig_line[0])):
                    #Store region points in line
                    if ((xmin < orig_line[0][i] < xmax) and (ymin < orig_line[1][i] < ymax)):
                        line[0].append(orig_line[0][i])
                        line[1].append(orig_line[1][i])
                if (smooth_type == "trilocal"):
                    smoother = smoothers.TriLocalSmoother(self._args.numpoints, self._args.repeat)
                elif (smooth_type == "integral"):
                    smoother = smoothers.IntegralSmoother(self._args.matchpt, None, self._args.interp)
                # Tri-integral smoothing disabled because it never worked correctly.
#                elif (smooth_type == "triintegral"):
#                    trismoother = smoothers.TriLocalSmoother(self._args.numpoints, self._args.repeat)
#                    smoother = smoothers.IntegralSmoother(self._args.matchpt, trismoother, self._args.interp)
                elif (smooth_type == "acute"):
                    smoother = smoothers.AcuteAngleRepair(self._args.angle, 1)
                else:
                    print ("Invalid smooth type")
                    return
                smooth_line = smoother.applySmooth(list(line[0]), list(line[1]), line[0][0], line[0][-1])
                #Add points outside of smoothed region
                for i in range(len(orig_line[0])):
                    if i < start or i > end:
                        smooth_line[0].insert(i, orig_line[0][i])
                        smooth_line[1].insert(i, orig_line[1][i])
                # Replot the line with changed points highlighted
                line_set.movable.set_x_data_y_data_moved_points(smooth_line)

                # Rehighlight selected points
                datapoints = line_set.movable.get_xy_data()
                new_set_points = []
                for i in range(0, len(line_set._set_index)):
                    new_set_points.append(datapoints[line_set._set_index[i]])
                line_set.movable.set_block_select(new_set_points)

        self._update_ghost_points()
        self.check_move_point()
        print ("Done")

    def undo(self):
        """ Undo the last point manipulation of any kind
        """
        for region in self._regions:
            region.undo()
            

    # Point movement ###########################################################

    def attempt_begin_move_point(self, event):
        """ Find (region index, point index) for the event point if it is
        within epsilon tolerance. If it is not, return None for the point index.
        """
        region_index = self.get_region_index(event.x, event.y)
        region = self._regions[region_index]
        region.attempt_begin_move_point(event)
        if region.move_point_in_progress():
            self._moving_point_region_index = region_index
            self._print_move_point_begin()

    def attempt_get_set(self, xmin, xmax, ymin, ymax, xdisp, ydisp):
        """ Find points within bounds and store them.
        """
        region_index = self.get_region_index(xdisp, ydisp)
        self._moving_set_region_index = region_index
        region = self._regions[region_index]
        region.attempt_get_set(xmin, xmax, ymin, ymax)

    def begin_move_set(self, event):
        """ Pass Lines the first event xy data
        """
        self._regions[self._moving_set_region_index].begin_move_set(event)

    def finish_move_set(self):
        """ Recalculate and redraw fit curve
        """
        self._update_ghost_points()
        self._regions[self._moving_set_region_index].check_move_point()
        
        # Update adjacent regions if 1d data with fitted line
        if (self._args.in_eos_file_base == None):
            if (self._moving_set_region_index != 0):
                self._regions[self._moving_set_region_index-1].check_move_point()
            if (self._moving_set_region_index != len(self._regions)-1):
                self._regions[self._moving_set_region_index+1].check_move_point()

        # Re-order points
        self._regions[self._moving_set_region_index].finish_move_set()

    def cancel_move_set(self):
        """ Reset set indexes to []
        """
        if (self._moving_set_region_index != None):
            region = self._regions[self._moving_set_region_index]
            region.cancel_move_set()
            self._moving_set_region_index = None
            print ("Block Select Disabled")    

    def move_set(self, event):
        self._regions[self._moving_set_region_index].move_set(event)

    def rotate_set(self, event, xmin, xmax, ymin, ymax):
        self._regions[self._moving_set_region_index].rotate_set(event, xmin, xmax, ymin, ymax)

    def cancel_any_move_points(self):
        """ Reset any moving point indexes.
        """
        self._moving_point_region_index = None
        for region in self._regions:
            region.cancel_any_move_points()

    def finish_move_point(self):
        """ Recalculate and redraw the fit curve.
        """
        if self.move_point_in_progress():
            self._print_move_point_end()
            self._update_ghost_points()
            self._regions[self._moving_point_region_index].finish_move_point()

            # Update adjacent regions if 1d data with fitted line
            if (self._args.in_eos_file_base == None):
                if (self._moving_point_region_index != 0):
                    self._regions[self._moving_point_region_index-1].check_move_point()
                if (self._moving_point_region_index != len(self._regions)-1):
                    self._regions[self._moving_point_region_index+1].check_move_point()
            self._moving_point_region_index = None

    def move_point(self, event):
        """ Move the point in the data and redraw the line.
        """
        if self.move_point_in_progress():
            self._regions[self._moving_point_region_index].move_point(event)

    def check_move_point(self):
        for region in self._regions:
            region.check_move_point()

    def moving_point_info(self, xy_only=False):
        if self.move_point_in_progress():
            return self._regions[self._moving_point_region_index].\
                moving_point_info(xy_only)
        else:
            return '(no moving point)'

    def move_point_in_progress(self):
        return self._moving_point_region_index is not None

    def _print_move_point_begin(self):
        print ('Moving point: %s ' % (self.moving_point_info()))

    def _print_move_point_end(self):
        print ('to: %s ' % (self.moving_point_info(xy_only=True)))

    def _update_ghost_points(self):
        '''Update points used to calculate fitted curve. 
           Uses a few points on either side of region
        '''
        # Return if 2d data
        if (self._args.in_eos_file_base != None):
            return

        # Get changed points
        pointlist = []
        for i in range (0, len (self._regions)):
            pointlist.append(
                sorted(self._regions[i]._line_sets.get_data_sets().get_only_set()))

        # Pass list with N extra points on either side
        for i in range(0, len (self._regions)):
            newlist = copy.copy(pointlist[i])
            for j in range(0, self._args.overlap):
                if ((i != 0) and
                    (len (pointlist[i-1]) > j)):
                    newlist.insert(0, pointlist[i-1][-1-j])
                if ((i != len (self._regions)-1) and
                    (len (pointlist[i+1]) > j)):
                    newlist.append(pointlist[i+1][j])

            self._regions[i]._line_sets.update_ghost_points(newlist)

    def _add_point(self, event):
        """
        Adds a point when left clicking
        """
        region_index = self.get_region_index(event.x, event.y)
        region = self._regions[region_index]
        region.add_point(event)

    def _remove_point(self, event):
        """
        removes a point when right clicking
        """
        region_index = self.get_region_index(event.x, event.y)
        region = self._regions[region_index]
        region.remove_point(event)
        # Redraw fit curve
        self._moving_set_region_index = region_index
        self.finish_move_set()

    # END Point movement #######################################################



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

from abc import ABCMeta, abstractmethod
import copy
import os
import re
from operator import itemgetter
import numpy as np
import pprint
import curvallis.curve_editing.eos_data_io as eos_data_io
import curvallis.curve_editing.configargparse as configargparse
import curvallis.curve_editing.curve_fitters as cf
from curvallis.version import version as VERSION_STRING

""" Supports getting point data into and out of the curve_editor.  Can read from
a file, generate, or provide predefined points.  Can write to file(s).

    Input points:
        Initial points to be plotted and moved.  May consist of multiple curves,
        or sets of points, to be managed separately.
    Background points:
        Unchanging points to display as reference data
"""

_INFINITY = float ('inf')
_NEG_INFINITY = float('-inf')
_NO_DECIMATE = 0
_POINT_FILE_HEADER_LINE = '#rho                  Data'

# Command line arguments =======================================================

def define_args(parser):
    """ Define the input, output, scaling, and decimation command-line arguments.
    """
    # Would like this to be a mutually-exclusive group with one arg required,
    # but those don't have title and description parms and don't group their
    # help separately. Mutual exclusion is tested for in _process_args.
    input_group = parser.add_argument_group(
        title='Inputs',
        description='These arguments each specify a different input data '
                    'source.  Use exactly one.')
    input_group.add_argument(
        '--in_eos_file_base', metavar='<base path>. ',
        help='Read source data points from EOS data <base path>'
             '[default: %(default)s]. '
             'If this is set, --out_eos_file_base must also be set.')
    input_group.add_argument(
        '--input_file', metavar='<path>',
        type=configargparse.FileType('r'),
        help='Read this input file for source data points')
    input_group.add_argument(
        '--parabola_in', action='store_true',
        help='(For developer testing) Use the internal noisy input parabola '
             'as source data [default: %(default)s]')
    input_group.add_argument(
        '--predefined_in', action='store_true',
        help='(For developer testing) Use the internal predefined data as '
             'source data [default: %(default)s]')

    # Would like this to be a mutually-exclusive group, but those don't
    # have title and description parms and don't group their help separately.
    # Mutual exclusion is tested for in _process_args.
    output_group = parser.add_argument_group(
        title='Outputs',
        description='These arguments each specify where the changed data '
                    'is written to.  Use no more than one.')
    output_group.add_argument(
        '--out_eos_file_base', metavar='<base path>',
        help='Use the EOS file <base path> for '
             'the adjusted/moved data points [default: %(default)s]. '
             'If this is set, --in_eos_file_base must also be set.')
    output_group.add_argument(
        '--output_file_name', metavar='<path>',
        # Not setting "type=argparse.FileType('w')" because that opens the
        # file for writing in the argument parser, and that conflicts with
        # opening it for writing later:
        help='Use this output file for the adjusted/moved data points '
             '[default: %(default)s]')
    output_group.add_argument(
        '--pressure_file_name', metavar='<path>',
        # See '--output_file_name' comment
        help='Use this output file for the pressure value at each point '
             '[default: %(default)s]')
    output_group.add_argument(
        '--bulk_mod_file_name', metavar='<path>',
        # See '--output_file_name' comment
        help='Use this output file for the bulk modulus value at each point '
             '[default: %(default)s]')
    output_group.add_argument(
        '--bulk_mod_prime_file_name', metavar='<path>',
        # See '--output_file_name' comment
        help='Use this output file for the bulk modulus prime value at each point '
             '[default: %(default)s]')
    output_group.add_argument(
        '--gamma_file_name', metavar='<path>',
        # See '--output_file_name' comment
        help='Use this output file for the Gruneisen gamma value at each point '
             '[default: %(default)s]')
    shifts = parser.add_argument_group(
        title='Shifts, Limits, and Point Exclusion',
        description='These arguments limit, exclude or alter the input data.')
    shifts.add_argument(
        '--decimate', metavar='<count>',
        type=int,
        help='Only use <count> points from the input data, minimum 2 '
             '[default: all points]')
    shifts.add_argument(
        '--step', metavar='<step>',
        type=int,
        help='Only use every <step>th input point.  e.g. --step 2 only' +
             ' uses every other point, minimum 1 [default: %(default)s]')
    shifts.add_argument(
        '--t_step', metavar='<step>',
        type=int,
        help='Only use every <step>th isotherm.  e.g. --step 2 only' +
             ' uses every other isotherm, minimum 1 [default: %(default)s]')
    shifts.add_argument(
        '--x_include','--xinclude', nargs=2, metavar=('<lo>', '<hi>'),
        type=float,
        help='Only consume input points when <lo> > x > <hi> [default: %(default)s]')
    shifts.add_argument(
        '--x_scale', '--xscale', nargs=3, metavar=('<scale>', '<lo>', '<hi>'),
        type=float,
        help='Multiply the x value of each input point by <scale> '
             'when <lo> > x > <hi> [default: %(default)s] '
             '>>NOTE<< applied AFTER --x_include')
    shifts.add_argument(
        '--x_shift', '--xshift', nargs=3,
        type=float,
        help='Add <shift> to the x value of each input point '
             'when <lo> > x > <hi> [default: %(default)s] '
             '>>NOTE<< applied AFTER --x_scale',
        metavar=('<shift>', '<lo>', '<hi>'))
    shifts.add_argument(
        '--y_include','--yinclude', nargs=2, metavar=('<lo>', '<hi>'),
        type=float,
        help='Only consume input points when <lo> > y > <hi> [default: '
             '%(default)s]')
    shifts.add_argument(
        '--y_scale', '--yscale', nargs=3, metavar=('<scale>', '<lo>', '<hi>'),
        type=float,
        help='Multiply the y value of each input point by <scale> '
             'when <lo> > x > <hi> [default: %(default)s].  '
             '>>NOTE<< uses x to trigger, not y.  '
             '>>NOTE<< applied BEFORE x adjustments, and AFTER --y_include.  ')
    shifts.add_argument(
        '--y_shift', '--yshift', nargs=3, metavar=('<shift>', '<lo>', '<hi>'),
        type=float,
        help='Add <shift> to the y value of each input point '
             'when <lo> > x > <hi> [default: %(default)s] '
             '>>NOTE<< uses x to trigger, not y.  ')
    shifts.add_argument(
        '--t_include','--tinclude', nargs=2, metavar=('<lo>', '<hi>'),
        type=float,
        help='Only consume input points when <lo> > t > <hi> [default: '
             '%(default)s]')
    shifts.add_argument(
        '--v_axis', action='store_true',
        help='Use v as x axis instead or rho. Note curve fits don\'t work like this.')

    # TODO: Add the rest of the arguments in alphabetical order, so they come
    # out that way in help. (Even though the arguments are defined in different
    # modules)
    parser.add_argument(
        '--background_file', metavar='<path>', nargs='+',
        type=configargparse.FileType('r'),
        help='Use this input file for background (unchanging) data points '
             '[default: %(default)s]')
    parser.add_argument(
        '--curve_output_file_name', metavar='<path>',
        # Not setting "type=argparse.FileType('w')" because that opens the file
        # for writing in the argument parser, and that conflicts with opening it
        # for writing later.
        help='Use this output file for data points generated from the '
             'current fit curve [default: %(default)s]')
    parser.add_argument(
        '--eos_function',
        help='Only read in the data for this function'
             '[default: %(default)s]')
    # parser.add_argument(
    #     '--yaml', action='store_true',
    #     dest='use_yaml_io',
    #     help='Use YAML format for input and output files [default: %(default)s]')

    output_group = parser.add_argument_group(
        title='View',
        description='Change the view but don\'t change the actual data')
    output_group.add_argument(
        '--x_max', metavar='<num>',
        help='Set the maximum x value for the window')
    output_group.add_argument(
        '--x_min', metavar='<num>',
        help='Set the minimum x value for the window')
    output_group.add_argument(
        '--y_max', metavar='<num>',
        help='Set the maximum y value for the window')
    output_group.add_argument(
        '--y_min', metavar='<num>',
        help='Set the minimum y value for the window')


    parser.set_defaults(
        # Input group:
        in_eos_file_base=None,
        input_file=None,
        parabola_in=False,
        predefined_in=False,

        # Output group:
        out_eos_file_base=None,
        output_file_name='moved_points_out.dat',
        pressure_file_name='E2P.dat',
        bulk_mod_file_name='P2B.dat',
        bulk_mod_prime_file_name='P2Bprime.dat',
        gamma_file_name='Theta2Gamma.dat',

        # Shifts group:
        decimate=_NO_DECIMATE,
        step=1,
        t_step=1,
        x_include=[_NEG_INFINITY, _INFINITY],
        x_scale=[1.0, _NEG_INFINITY, _INFINITY],
        x_shift=[0.0, _NEG_INFINITY, _INFINITY],
        y_include=[_NEG_INFINITY, _INFINITY],
        y_scale=[1.0, _NEG_INFINITY, _INFINITY],
        y_shift=[0.0, _NEG_INFINITY, _INFINITY],
        t_include=[_NEG_INFINITY, _INFINITY],

        # other:
        background_file=[],
        curve_output_file_name='fit_curve_out.dat',
        use_yaml_io=False,
        eos_function='all',
    )

def process_args (parser, args):
    def check_no_config_file():
        if args.no_config_file:
            if args.config_file:
                parser.error("cannot use both config_file and no_config_file")

    def check_decimate_and_step():
        if args.decimate != parser.get_default('decimate') and \
                        args.decimate < 2:
            parser.error("decimate value was %s, must be at least 2." %
                              args.decimate)
        if args.step < 1:
            parser.error("step value was %s, must be at least 1." %
                              args.step)

    def check_input_arg():
        """ Make sure there is exactly one of these.
        """
        arg_names = ['in_eos_file_base', 'input_file', 'parabola_in',
                     'predefined_in']
        num_set = 0
        for arg in arg_names:
            if getattr(args, arg) != parser.get_default(arg):
                num_set += 1
        if num_set != 1:
            parser.error('must set exactly one of: %s' % arg_names)

    def check_output_arg():
        """ Make sure there is no more than one of these.
        """
        arg_names = ['out_eos_file_base', 'output_file_name']
        num_set = 0
        for arg in arg_names:
            if getattr(args, arg) != parser.get_default(arg):
                num_set += 1
        if num_set > 1:
            parser.error('must set no more than one of: %s' % arg_names)

    def check_eos_args():
        """ If there is an EOS input arg, there must be an EOS output arg,
        and vice-versa.
        """
        if (args.in_eos_file_base is None) != (args.out_eos_file_base is None):
            parser.error('Both --in_eos_file_base and --out_eos_file_base '
                         'must be set, or neither must be set.')

    def check_required_rho0():
        """ If any fitter is sandiapc or gammapoly, the user must
            be forced to enter the constant "rho0" instead of
            Curvallis gussing it.
        """
        r = re.compile('('+'|'.join([cf.GammaPoly.name_prefix, cf.GammaPolyV.name_prefix])+')\d+')
        if args.rho0 is None and ("sandiapc" in args.fit_type or list(filter(r.match, args.fit_type))):
            parser.error('If using fitter "sandiapc", "{0}", or "{1}", you must give a value for "rho0_guess".'.format(
                cf.GammaPoly.name_prefix, cf.GammaPolyV.name_prefix))

    def check_region_divisions():
        """ Prevent the user from trying to define regions with both the
            'region_bound' and the 'region_data_points' command line arguments
            at the same time.
        """
        if args.region_bound != None and args.region_data_points != None:
            parser.error('Either "--region_bound" or "--region_data_points" '
                         'is used to define regions, not both at the same time.')

    def check_numpoints():
        """ Prevent the user from giving the numpoint argument an even number,
            which causes problems later in the program.
        """
        if args.numpoints % 2 == 0:
            parser.error('"--numpoints" must be an odd number, not "%s".' %
                          args.numpoints)

    check_no_config_file()
    check_decimate_and_step()
    check_input_arg()
    check_output_arg()
    check_eos_args()
    check_required_rho0()
    check_region_divisions()
    check_numpoints()

class XY_Limits(object):
    """ Contains an x range and a y range.
    """
    def __init__(self):
        self.x_max = _NEG_INFINITY
        self.x_min = _INFINITY
        self.y_max = _NEG_INFINITY
        self.y_min = _INFINITY

        self.view_x_max = self.x_max
        self.view_x_min = self.x_min
        self.view_y_max = self.y_max
        self.view_y_min = self.y_min

    def update_using_x_values_y_values(self, x_values, y_values):
        """ Incorporate a new set of x_values' and y_values' mins and maxes.  May
        or may not change self.
        """
        self.x_max = max(self.x_max, *x_values)
        self.x_min = min(self.x_min, *x_values)
        self.y_max = max(self.y_max, *y_values)
        self.y_min = min(self.y_min, *y_values)

    def update_using_xy_values(self, xy_values):
        """
        :param xy_values: [[x,y], [x,y], [x,y]...]
        """
        if len(xy_values) > 0:
            x_values = [point[0] for point in xy_values]
            y_values = [point[1] for point in xy_values]
            self.update_using_x_values_y_values(x_values, y_values)

    def update_using_limits(self, other_limits):
        """
        :param other_limits: XY_Limits
        """
        self.x_max = max(self.x_max, other_limits.x_max)
        self.x_min = min(self.x_min, other_limits.x_min)
        self.y_max = max(self.y_max, other_limits.y_max)
        self.y_min = min(self.y_min, other_limits.y_min)

class Data_Sets(object):
    """ One or more lists of points.  If more than one, each has a name.  Each
    is a list of x,y tuples.
    """
    # TODO: Is a data set [[x1,y1], [x2,y2]] or [[x1,x2], [y1,y2]]?
    def __init__(self):
        self._sets = {}
        self._all_sets_xy_limits = XY_Limits()

    def add_set(self, points, name):
        """ Adds a new data set and updates extent
        """
        assert name not in self._sets
        self._sets[name] = points
        self._expand_limits_using(points)

    def add_or_append_to_set(self, name, points):
        """ Adds a new data set or adds data to an existing one, and updates
        extent.
        """
        if name in self._sets:
            for i in range(0, len(points)):
                self._sets[name].append(points[i])
        else:
            self._sets[name] = points
        self._expand_limits_using(points)

    def decimate(self, args):
        for name, data_set in self._sets.items():
            self._sets[name] = _decimate(data_set, name, args)

    def do_include_scale_and_shift(self, args):
        for name, data_set in self._sets.items():
            self._sets[name] = _do_include_scale_and_shift(data_set, name, args)

    def _expand_limits_using(self, points):
        self._all_sets_xy_limits.update_using_xy_values(points)

    def get_copy(self):
        return copy.deepcopy(self)

    def get_name_set_items(self):
        return self._sets.items()

    def get_names(self):
        return self._sets.keys()

    def get_num_points_in_fullest_set(self):
        if self.num_sets() == 0:
            result = 0
        else:
            result = _NEG_INFINITY
            for data_set in self._sets:
                result = max(result, len(data_set))
        return int(result)

    def get_num_points_in_sparsest_set(self):
        if self.num_sets() == 0:
            result = 0
        else:
            result = _INFINITY
            for data_set in self._sets:
                result = min(result, len(data_set))
        return int(result)

    def get_only_set(self):
        assert self.num_sets() == 1
        return list(self._sets.values())[0]

    def get_set(self, name):
        return self._sets[name]

    def get_set_values(self):
        return self._sets.values()

    def get_x_slice(self, x_low, x_high):
        """
        :param x_low:
        :param x_high:
        :return: Data_Sets with all points having x_low <= x <= x_high
        """
        result = Data_Sets()
        for name, points in self._sets.items():
            sliced_points = []
            for point in points:
                if x_low <= point[0] <= x_high:
                    sliced_points.append(point)
            result.add_set(sliced_points, name)
        return result

    def get_y_slice(self, y_low, y_high):
        """
        :param y_low:
        :param y_high:
        :return: Data_Sets with all points having y_low <= y <= y_high
        """
        result = Data_Sets()
        for name, points in self._sets.items():
            sliced_points = []
            for point in points:
                if y_low <= point[1] <= y_high:
                    sliced_points.append(point)
            result.add_set(sliced_points, name)
        return result

    def get_xy_limits(self):
        return self._all_sets_xy_limits

    def get_specific_xy_limits (self, x_max, x_min, y_max, y_min):
        
        resultx = self.get_x_slice(x_min, x_max)
        result = resultx.get_y_slice(y_min, y_max)

        return result._all_sets_xy_limits
#        return self._all_sets_xy_limits

    def num_sets(self):
        return len(self._sets)

    def print_stats(self, prefix):
        for name, data_set in self._sets.items():
            print_data_stats(data_set, '%s: %s' % (prefix, name))

    def sort_by_x(self):
        for name, data_set in self._sets.items():
            data_set.sort(key=itemgetter(0))
            self._sets[name] = data_set

class Manager(object):
    """ Manages the input, background, and output objects.
    """
    def __init__(self, args):
        self._io_adapter = None
        self._args=args
        self._do_curves = self._args.in_eos_file_base is None
        check_output_file(self._args.output_file_name)
        if self._do_curves:
            check_output_file(self._args.curve_output_file_name)
        if self._args.print_E2P:
            check_output_file(self._args.pressure_file_name)
        if self._args.print_P2B:
            check_output_file(self._args.bulk_mod_file_name)
        if self._args.print_P2Bprime:
            check_output_file(self._args.bulk_mod_prime_file_name)
        if self._args.print_theta2gamma:
            check_output_file(self._args.gamma_file_name)
        self._input = None
        self._is_eos_data = self._args.in_eos_file_base is not None
        self._background = None
        self._output = None
        # TODO: move to process_args?
        if self._args.predefined_in:
            self._io_adapter = Predefined_Adapter(self._args)
        elif self._args.parabola_in:
            self._io_adapter = Noisy_Parabola_Adapter(self._args)
        elif self._args.input_file is not None:
            self._io_adapter = File_Adapter(self._args)
        elif self._is_eos_data:
            self._io_adapter = EOS_Files_Adapter(self._args)
        else:
            raise RuntimeError('No input data parameter found.')

    def get_movable_data_sets(self):
        """ Gets data either from a file, a noisy parabola, a few predefined
        points, or an EOS data file.
    
        :return: Data_Sets
        """
        result = self._io_adapter.get_data_sets()
        
        if (self._args.v_axis == True):
            keylist = list(result._sets.keys())
            for i in range(0, len(keylist)):
                newlist = []
                for j in range(0, len(result._sets[keylist[i]])):
                    xval = 1 / result._sets[keylist[i]][j][0]
                    yval = result._sets[keylist[i]][j][1]
                    newlist.append((xval, yval))
                result._sets[keylist[i]] = newlist
                result._all_sets_xy_limits.update_using_xy_values(
                    result._sets[keylist[i]])

        result.print_stats('input points')
        result.do_include_scale_and_shift(self._args)
        result.sort_by_x()
        result.decimate(self._args)
        return result

    def get_background_data_sets(self):
        """
        :return: Data_Sets
        """
        result = Data_Sets()
        i = 0
        for back_file in self._args.background_file:
            result.add_set(get_points_from_file(back_file),
                           name='background points'+str(i))
            result.print_stats('background_points'+str(i))
            i = i + 1
        return result

    def is_eos_data(self):
        """
        :return: bool
        """
        return self._is_eos_data

    def write_movable_data_sets(self, data_sets):
        """
        Writes data sets(s) to file(s)
        """
        self._io_adapter.write_changed_data_sets(data_sets)

class IO_Adapter(object):
    __metaclass__ = ABCMeta
    def __init__(self, args):
        self._args = args

    @abstractmethod
    def get_data_sets(self):
        """
        :return: Data_Sets
        """

    def write_changed_data_sets(self, data_sets):
        """
        Default implementation writes the only data set out to a file.
        :param data_sets: Data_Sets containing one data set
        :return:
        """
        write_point_file(self._args.output_file_name, data_sets.get_only_set())


class Predefined_Adapter(IO_Adapter):
    """ Predefined data in one set.
    """
    def get_data_sets(self):
        """ Return a few predefined points

        :return:  Data_Sets with one set
        """
        result = Data_Sets()
        points = [
            (-3.75, -2.0),
            # (0.35, -1.1),
            (0.375, 2.0),
            (0.85, 1.15),
            (1.58, -2.57),
            (2.0, -0.5),
            # This point is out of order, but points is sorted in __init__:
            (0.35, -1.1),
            (2.2, 3.2),
            (3, -4.0)]
        result.add_set(points, name='predefined_points')
        return result


class Noisy_Parabola_Adapter(IO_Adapter):
    """ A noisy parabola of points in one set.
    """
    def get_data_sets(self):
        """ Generate initial data as a noisy parabola

        :return:  Data_Sets with one set
        """
        result = Data_Sets()
        num_pairs = 50

        # generate a noisy parabola
        x = np.linspace(0,100,num_pairs)
        # noinspection PyTypeChecker
        parabola = x**2
        noise = np.random.normal(0, 300, num_pairs)
        y = parabola + noise
        result.add_set(points=[[x[i],y[i]] for i in range(len(x))], name='noisy_parabola_points')
        return result


class File_Adapter(IO_Adapter):
    """ A file with one point per line, in one set.
    """
    def get_data_sets(self):
        """ Get points from a file

        :return:  Data_Sets with one set
        """
        input_file = self._args.input_file
        # Check if input file is a valid data file
        ########################################
        opened_input_file = open(input_file.name,"r")
        for line in opened_input_file:
            if(line[0] == '['):
                print("Error: --input_file\nA data file was expected but a non-data file was given.")
                exit()
        opened_input_file.close()
        ########################################
        result = Data_Sets()
        result.add_set(points=get_points_from_file(input_file),
                       name=os.path.basename(input_file.name))
        return result

class EOS_Files_Adapter(IO_Adapter):
    """ A pair of .info, .dat EOS-formatted files, with one or more data sets.

    _eos_data: A name-indexed dictionary of sections. Each section has a
    temperature-indexed dictionary of isotherms. Each isotherm is a list of x,y
    points.
    """

    def __init__(self, args):
        super(EOS_Files_Adapter, self).__init__(args)
        self._eos_data = None

    def get_data_sets(self):
        """ Get one or more data sets from a pair of EOS files.  Collapse the
        list of sections, each containing a list of isotherms, into one list of
        data sets.

        :return: Data_Sets
        """
        result = Data_Sets()
        eos_file_base = self._args.in_eos_file_base
        self._eos_data = eos_data_io.Data(
#                use_function=self._args.eos_function)
                use_function='all')
        self._eos_data.read(eos_file_base)
        for section in self._eos_data.sections.values():
            # Skip over sections if eos_function defined
            if self._args.eos_function != 'all' and section.name != self._args.eos_function:
                continue
            # Skip over N isotherms defined by t_step
            newlist = list(section.isotherms.values())[::self._args.t_step]
            for isotherm in newlist:
                # Check if within t_include bounds
                if self._args.t_include[0] < isotherm.temp < self._args.t_include[1]:
                    result.add_set(
                        points=isotherm.points,
                        name=self._to_data_set_name(section.name, isotherm.temp))
        return result

    @staticmethod
    def _to_data_set_name(section_name, isotherm_temp):
        """
        :param section_name:  str
        :param isotherm_temp: float
        :return:              str
        """
        return '%s: %.15E' % (section_name, isotherm_temp)

    def write_changed_data_sets(self, data_sets):
        """ Write data sets back out to data files.  First, expand the flat list
        of data sets back into a list of sections, each containing a list of
        isotherms.

        :param data_sets: Data_Sets
        """
        assert self._eos_data is not None
        for name, points in data_sets.get_name_set_items():
            section_name, isotherm_temp = self._to_section_name_isotherm_temp(name)
            isotherm = eos_data_io._Isotherm(temperature=isotherm_temp,
                                          num_points=len(points))
            isotherm.points=points
            self._eos_data.sections[section_name].isotherms[isotherm_temp] = \
                    isotherm


        self._eos_data.write(self._args.out_eos_file_base)

    @staticmethod
    def _to_section_name_isotherm_temp(name):
        """ Inverse of _to_data_set_name
        :param name: str
        :return: str, float
        """
        words = name.split(':')
        name = words[0]
        temp = float(words[1])
        return name, temp


# Reading input data ===========================================================

def get_points_from_file(input_file):
    """ Get points data from a file if the file is not None

    :file:   input file object
    :return: list of x,y tuples
    """
    if input_file is None:
        result = []
    else:
        result = []
        print ('Current directory is "%s"' % os.getcwd())
        print ('Reading data from "%s"' % input_file.name)
        line_number=0
        with open(input_file.name) as in_file:
            for line in in_file:
                line_number += 1
                try:
                    strings = line.split()
                    values = list(range(2))
                    values_found = 0
                    for substr in strings:
                        # Split can give us strings of length 0:
                        if len(substr) == 0:
                            pass
                        # This and all remaining substrings are part of a comment:
                        elif substr[0] == '#':
                            break
                        else:
                            # Skip anything but a good float:
                            try:
                                value = float(substr)
                                # Count all the values found so we can check
                                # for more or less than two:
                                values_found +=1
                                # Only want the first two values:
                                if values_found <= 2:
                                    values[values_found - 1] = value
                                else:
                                    print ('!!! Warning: File "%s", line %s, "%s"' %
                                           (input_file.name, line_number, line))
                                    print ('!!! Found more than two float values.'
                                           '  Extra float value %.15E ignored.' % value)
                            except ValueError:
                                print ('!!! Warning: File "%s", line %s, "%s"' %
                                       (input_file.name, line_number, line))
                                print ('!!! Found non-float string "%s" - ignored.' % substr)
                    # Only use the values in this line if there were at least two of them:
                    if values_found == 1:
                        print ('!!! Warning: File "%s", line %s, "%s"' %
                               (input_file.name, line_number, line))
                        print ('!!! Found only one float value.  Line ignored.')
                    elif values_found >= 2:
                        result.append(tuple(values))
                except Exception:
                    print ('*** In file %s at line %s ("%s"), got exception:' % (input_file.name, line_number, line))
                    raise
        print('Lines read: %s' % line_number)
    return result

def _do_include_scale_and_shift (points_in, name, args):
    """ Apply x_include, y_include, y_scale, y_shift, x_scale, then x_shift.
    Do y scale and shift first because their logic depends on the value of x,
    which is changed by x scale and shift.

    Could do this in six individual loops to make it faster, but the
    requirements for this code change a lot, and clarity,
    maintainability, and robustness take priority.

    :points_in: list of [x,y] points
    :name     ; name of this isotherm.  Like #
    :returns    list of [x,y] points, possibly modified
    """
    class Changes(object):
        def __init__(self):
            self.excluded=0
            self.shifted=0
            self.scaled=0

    class Point_changes(object):
        def __init__(self):
            self.x=Changes()
            self.y=Changes()

    point_changes = Point_changes()
    points_changed = False

    def apply_x_scale (point):
        scale, low_lim, hi_lim = args.x_scale
        if scale != 1.0:
            if low_lim < point[0] < hi_lim:
                point[0] = point[0] * scale
                point_changes.x.scaled += 1
        return point

    def apply_x_shift (point):
        shift, low_lim, hi_lim = args.x_shift
        if shift != 0.0:
            if low_lim < point[0] < hi_lim:
                point[0] = point[0] + shift
                point_changes.x.shifted += 1
        return point

    def apply_y_scale (point):
        """Note that we are checking x for the limit, not y.
        """
        scale, low_lim, hi_lim = args.y_scale
        if scale != 1.0:
            if low_lim < point[0] < hi_lim:
                point[1] = point[1] * scale
                point_changes.y.scaled += 1
        return point

    def apply_y_shift (point):
        """Note that we are checking x for the limit, not y.
        """
        shift, low_lim, hi_lim = args.y_shift
        if shift != 0.0:
            if low_lim < point[0] < hi_lim:
                point[1] = point[1] + shift
                point_changes.y.shifted += 1
        return point

    def prefix_print(msg):
        print('--- %s: %s' % (name, msg))

    points_out = list()
    for point in points_in:
        if args.x_include[0] < point[0] < args.x_include[1]:
            if args.y_include[0] < point[1] < args.y_include[1]:
                lpoint = list(point)
                # Do y scale and shift first because their logic  depends on the
                # value of x, which is changed by x scale and shift.
                point = apply_y_scale(lpoint)
                point = apply_y_shift(lpoint)
                point = apply_x_scale(lpoint)
                point = apply_x_shift(lpoint)
                points_out.append(tuple(lpoint))
            else:
                point_changes.y.excluded += 1
        else:
            point_changes.x.excluded += 1

    if point_changes.x.excluded > 0:
        points_changed = True
        prefix_print ('Excluded %s points with X not between %.15E and %.15E' %
                      (tuple([point_changes.x.excluded] + args.x_include)))
    if point_changes.y.excluded > 0:
        points_changed = True
        prefix_print ('Excluded %s points with Y not between %.15E and %.15E' %
                      (tuple([point_changes.y.excluded] + args.y_include)))
    if point_changes.x.scaled > 0:
        points_changed = True
        prefix_print ('Scaled %s points\' X values (with %.15E < X < %.15E) by %.15E' %
                      (tuple([point_changes.x.scaled] + args.x_scale)))
    if point_changes.y.scaled > 0:
        points_changed = True
        prefix_print ('Scaled %s points\' Y values (with %.15E < X < %.15E) by %.15E' %
                      (tuple([point_changes.y.scaled] + args.y_scale)))
    if point_changes.x.shifted > 0:
        points_changed = True
        prefix_print ('Shifted %s points\' X values (with %.15E < X < %.15Es) by %.15E' %
                      (tuple([point_changes.x.shifted] + args.x_shift)))
    if point_changes.y.shifted > 0:
        points_changed = True
        prefix_print ('Shifted %s points\' Y values (with %.15E < X < %.15E) by %.15E' %
                      (tuple([point_changes.y.shifted] + args.y_shift)))
    if points_changed:
        print_data_stats(points_out,
                         'Points after excluding, scaling, and/or shifting')
    return points_out

# Data Stats ===================================================================

def print_data_stats(xy_values, name):
    print('---------------------------------------------------------------')
    print ('%s:' % name)
    if len(xy_values) == 0:
        print ('None')
    else:
        print ('Number of points: %s' % len(xy_values))
        if len(xy_values) > 0:
            x_vals = [point[0] for point in xy_values]
            y_vals = [point[1] for point in xy_values]
            print ('X range: %.15E .. %.15E' % (min(x_vals), max(x_vals)))
            print ('Y range: %.15E .. %.15E' % (min(y_vals), max(y_vals)))
    print('---------------------------------------------------------------')

# Decimation ===================================================================

def _decimate(points_in, name, args):
    """ Call _decimate_reduce and _decimate_step

    :param points_in: sequence of (x, y) values
    :return:          sequence of possibly fewer(x, y) values
    """
    result = points_in
    if args.step > 1:
        result = _decimate_by_stepping(result, args.step, name)
    if args.decimate != _NO_DECIMATE:
        result = _decimate_to_limit(result, args.decimate, name)
    return result

def _decimate_by_stepping (points, step, name):
    """ Return every 2nd, 3rd... point as specified in step.  Be sure
    to preserve the first and last points

    :param points: sequence of (x, y) values
    :return:       sequence of possibly fewer(x, y) values
    """
    # Preserve the first and last input points (need to preserve both if the
    # step is bigger than len(points):
    result = [points[0]] + points[step:-2:step] + [points[-1]]
    print_data_stats(result, '%s: Points, decimated stepwise' % name)
    return result

def _decimate_to_limit(points, max_count, name):
    """ Return at most max_count points, evenly distributed, including the
    first and last points.  If max_count > len(points, return points unchaged.

    :param points:    sequence of (x, y) values
    :param max_count: maximum number of points to return
    :return:          sequence of possibly fewer(x, y) values
    """
    if len(points) <= max_count:
        result = points
    else:
        in_size = len(points)
        out_size = max_count
        result = []
        # count goes from 0 to out_size - 1, so portion goes from from 0.0 to 1.0:
        for count in range (out_size):
            # Calculate each index without any cumulative errors:
            portion = float(count) / float(out_size - 1)
            index = int(portion * float(in_size - 1))
            result.append(points [index])
    print_data_stats(result, '%s: Points, decimated to limit' % name)
    return result

# Output =======================================================================

def check_output_file(file_name):
    """ If the output file exists, ask if it's ok to overwrite it.  If so,
    delete it. If not, exit cleanly.

    :param file_name: string
    """
    if os.path.exists (file_name):
        want_to = input('Do you want to overwrite the file "%s"? ' % file_name)
        if 'Y' in want_to or 'y' in want_to:
            os.remove(file_name)
        else:
            print ('Exiting')
            exit(0)

def write_point_file(file_name, points):
    print('Writing points to %s' % file_name)
    try:
        with open(file_name, 'w') as out_file:
            out_file.write (_POINT_FILE_HEADER_LINE + '\n')
            for point in points:
                x = point[0]
                y = point[1]
                out_file.write ('% .15E % .15E\n' % (x, y))

        print('DONE Writing points to %s' % file_name)
    except IOError as e:
        print('*** EXCEPTION %s %s' % (type(e), e))

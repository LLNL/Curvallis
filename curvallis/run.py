#!/usr/bin/env python
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


"""
Created on Fri Mar 13 15:34:13 2015

@author: reynolds12
@author: Paul Minner (minner2)
"""

from __future__ import print_function

import sys
import math

from matplotlib import pyplot, rcParams
from matplotlib.backend_bases import NavigationToolbar2, FigureManagerBase
from matplotlib.widgets import RectangleSelector
from curve_editing import curve_fitters, io, lines, regions, configargparse
from Tkinter import Tk, Label, Button, Entry

VERSION_STRING = '2015-12-01 11:02AM'

# Overwrite Panning and Zooming Functions
PAN_ENABLED = False
ZOOM_ENABLED = False

old_pan = NavigationToolbar2.pan
def new_pan(self, *args, **kwargs):
    global PAN_ENABLED
    global ZOOM_ENABLED
    ZOOM_ENABLED = False
    PAN_ENABLED = not PAN_ENABLED
    old_pan(self, *args, **kwargs)
NavigationToolbar2.pan = new_pan

old_zoom = NavigationToolbar2.zoom
def new_zoom(self, *args, **kwargs):
    global ZOOM_ENABLED
    global PAN_ENABLED
    PAN_ENABLED = False
    ZOOM_ENABLED = not ZOOM_ENABLED
    old_zoom(self, *args, **kwargs)
NavigationToolbar2.zoom = new_zoom

old_home = NavigationToolbar2.home
def new_home(self, *args, **kwargs):
    self.push_current()
    old_home(self, *args, **kwargs)
    self.push_current()
NavigationToolbar2.home = new_home

class CurveInteractor(object):
    """ Calculate a curve to fit the data, let user move points and recalculate.
    """
    RANGE_MARGIN = 0.1  # relative plot margin around data on initial plot

    # Init support =============================================================

    def __init__(self):
        self._args = None
        self._ax = None
        # Used in pyplot callbacks.  Not related to background_*:
        self.background = None
        self._background_data_sets = io.Data_Sets()
        self._background_line = None
        self._canvas = None
        self._figure = None
        self._input_data_sets = io.Data_Sets()
        self._io_manager = None
        self._parser = None
        self._quit_pending = False
        self._regions = None
        self._xlim_callback_id = None
        self._xy_limits = io.XY_Limits()
        self._define_args()
        self._add_points = False
        # interactive equation plotting
        self._iplot = {}
        # Disable to create own scaling funcion
        rcParams['keymap.xscale'] = ''
        rcParams['keymap.yscale'] = ''
        # Variables for block selection
        self._selector = None
        self._find_set = False
        self._move_set = False
        self._xmin = None
        self._xmax = None
        self._ymin = None
        self._ymax = None
        self._last_click = None

    def _define_args(self):
        parser = configargparse.ArgumentParser(
            description='Curve Editor.  Fits a curve to the input data, then '
                        'lets user change points and recalculate the curve. '
                        '>>Config File<< ',
            epilog='You only need to provide enough of each argument name to '
                   'make it unique. e.g. --i and --input_file are equivalent.  '
                   'You must provide exactly one source data argument.',
            default_config_files=['curve_editor.ini'],
            # args_for_setting_config_path=['--config_file'],
            # config_arg_help_message='Configuration file path(s) [default: %(default)s].'
            #                         '  Command line parms override config file parms.'
        )
        parser.add_argument('-v', '--version',
                            action='version', version='%(prog)s ' + VERSION_STRING)
        # noinspection PyProtectedMember
        parser.add_argument(
            '--config_file', is_config_file=True,
            help='Configuration file path(s) [default: %s].  Command line parms '
                 'override config file parms.' %
            parser._default_config_files, metavar='<path>')
        parser.add_argument( '--no_config_file', 
            help='Do not use any config file, even if a default is provided',
            action='store_true')
        # config_file default is implicitly 'config_file.ini'
        io.define_args(parser)
        regions.define_args(parser)
        curve_fitters.define_args(parser)
        self._parser = parser

    # END Init support =========================================================

    # Run support: =============================================================

    def run(self):
        """ Create original line, movable line, background "line" (just points),
        and fitted curve lines. Create figure and subplot, register the callbacks,
        and show the plot.

        The movable line and the curve lines are in the regions.
        """
        self._process_args()
        self._io_manager = io.Manager(self._args)
        # Read the data sets outside of Regions so we can discover the x, y
        # range of the data to use for the initial GUI size:
        self._input_data_sets = self._io_manager.get_movable_data_sets()
        self._background_data_sets = self._io_manager.get_background_data_sets()
        self._setup_GUI()
        self._regions = regions.Regions(
            ax=self._ax,
            args=self._args,
            input_data_sets=self._input_data_sets,
            xy_limits=self._xy_limits,
            io_manager=self._io_manager)
        self._register_callbacks()
        self._print_keymap()
        # Create rectangle selector for selecting multiple points
        self._selector = RectangleSelector(self._ax, self.line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle button
                                           spancoords='pixels')
        self._selector.set_active(False)

        # Main routine: display all figures and block until all figures have
        # been closed:
        pyplot.show()

    def _setup_GUI(self):
        self._figure = pyplot.figure('Curve Editor', figsize=(12, 8))
        self._canvas = self._figure.canvas
        # One row, one column, first subplot:
        self._ax = self._figure.add_subplot(1, 1, 1)
        # Minimize margins:
        self._figure.tight_layout()
        self._ax.set_xlabel('X')
        self._ax.set_ylabel('Y')
#        self._background_line = lines.Line(
#            self._ax, lines.line_attributes['background_points'])
        self._background_line = []
        self._plot_background_data()
        # Default logscale if 2d eos data
        if self._args.in_eos_file_base:
            self._ax.set_xscale('log')
            self._ax.set_yscale('log')
        # Needs self.ax set first:
        self._set_xlim_ylim()

    def _process_args(self):
        def fix_negative_scientific_notation_parms():
            """ argparse does not handle an argument with multiple float values
            properly if a value is negative and in scientific notation.  One
            workaround is to quote the number, with a leading space, e.g.
            " -1.0E15".  Here, we do that automatically.  We also handle -inf.

            Must be called before parser.parse_args()
            """
            for i, arg in enumerate(sys.argv):
                if (arg[0] == '-') and (arg[1].isdigit() or arg[1:4] == 'inf'):
                    sys.argv[i] = ' ' + arg

        fix_negative_scientific_notation_parms()
        self._args = self._parser.parse_args()
        io.process_args(self._parser, self._args)

    def _set_xlim_ylim(self):
        """ Find the max and min x and y of the input values.
            Set the canvas plotting limits to those, plus some padding.
        """
        def expanded_range(low, high, factor):
            """ Return low and high moved apart by abs(max-min) * factor
            """
            rangee = high - low
            move = rangee * factor / 2
            # print ('expanded_range: input range is %.15E .. %.15E' % (low, high))
            low -= move
            high += move
            # print ('expanded_range: output range is %.15E .. %.15E' % (low, high))
            return low, high

        self._xy_limits.update_using_limits(self._input_data_sets.get_xy_limits())
        self._xy_limits.update_using_limits(self._background_data_sets.get_xy_limits())

        # Set user inputted view limits
        if self._args.x_max:
            self._xy_limits.view_x_max = float (self._args.x_max)
        else:
            self._xy_limits.view_x_max = self._xy_limits.x_max
        if self._args.x_min:
            self._xy_limits.view_x_min = float (self._args.x_min)
        else:
            self._xy_limits.view_x_min = self._xy_limits.x_min
        if self._args.y_max:
            self._xy_limits.view_y_max = float (self._args.y_max)
        else:
            self._xy_limits.view_y_max = self._xy_limits.y_max
        if self._args.y_min:
            self._xy_limits.view_y_min = float (self._args.y_min)
        else:
            self._xy_limits.view_y_min = self._xy_limits.y_min

        # determine new fit based on view limits
        view_limits = self._input_data_sets.get_specific_xy_limits (
            self._xy_limits.view_x_max, self._xy_limits.view_x_min, 
            self._xy_limits.view_y_max, self._xy_limits.view_y_min)

        self._xy_limits.view_x_max = view_limits.x_max
        self._xy_limits.view_x_min = view_limits.x_min
        self._xy_limits.view_y_max = view_limits.y_max
        self._xy_limits.view_y_min = view_limits.y_min

        # Transform limits to display coordinates
        xmin_disp = self._ax.transData.transform((self._xy_limits.view_x_min,0))[0]
        xmax_disp = self._ax.transData.transform((self._xy_limits.view_x_max,0))[0]
        ymin_disp = self._ax.transData.transform((0,self._xy_limits.view_y_min))[1]
        ymax_disp = self._ax.transData.transform((0,self._xy_limits.view_y_max))[1]

        # Calculate expanded range in display coordinates
        xmin_exp, xmax_exp = expanded_range(xmin_disp, xmax_disp, self.RANGE_MARGIN)
        ymin_exp, ymax_exp = expanded_range(ymin_disp, ymax_disp, self.RANGE_MARGIN)

        # Convert expanded range data to data coordinates
        xmin_data = self._ax.transData.inverted().transform((xmin_exp,0))[0]
        xmax_data = self._ax.transData.inverted().transform((xmax_exp,0))[0]
        ymin_data = self._ax.transData.inverted().transform((0,ymin_exp))[1]
        ymax_data = self._ax.transData.inverted().transform((0,ymax_exp))[1]

        self._ax.set_xlim(xmin_data, xmax_data)
        self._ax.set_ylim(ymin_data, ymax_data)

    def _plot_background_data(self):
        if self._background_data_sets.num_sets > 0:
            # Plot each background data line
            for back_set in self._background_data_sets.get_set_iterator():
                if len(back_set) > 0:
                    self._background_line.append(
                        lines.Line(self._ax, lines.line_attributes['background_points']))
                    self._background_line[-1].plot_xy_data(back_set)

    # END Run support ==========================================================

    # Callbacks ================================================================

    # noinspection PyUnusedLocal
    def draw_callback(self, event):
        """ Event value is not used here.
        Called by plot system directly e.g. when resizing window.
        """
        #print('Got draw callback')
        self._background = self._canvas.copy_from_bbox(self._ax.bbox)
        self._draw()

    def button_press_callback(self, event):
        """ Whenever a mouse button is pressed
        """
        # Quit if panning or zooming
        global PAN_ENABLED
        global ZOOM_ENABLED
        if (PAN_ENABLED or ZOOM_ENABLED):
            # Disable block selection while zooming or panning
            self._move_set = False
            self._find_set = False
            self._selector.set_active(False)
            self._regions.cancel_move_set()
            self._canvas.draw()
            return

        if self._find_set == True:
            return
        elif self._move_set == True:
            if event.button == 1:
                # Move a region of points
                self._check_move_set(event)
            else:
                # Rotate a region of points
                self._check_rotate_set(event)
        elif self._add_points == True:
            # Add / Remove points
            if event.button == 1:
                self._regions._add_point(event)
            else:
                self._regions._remove_point(event)
            self._draw()
            self._attempt_begin_move_point(event)
        else:
            # Move Individual Points
            if event.button == 1:
                if event.inaxes is not None:
                    self._attempt_begin_move_point(event)

    def motion_notify_callback(self, event):
        """ On mouse movement
        """
        if (self._move_set == True and event.button == 1):
            # Move set of points
            if event.inaxes is not None:
                self._attempt_move_set(event)
        elif (self._move_set == True and event.button == 3):
            # Rotate set of points
            self._rotate_set(event)
        elif event.button == 1:
            # Move individual point
            if event.inaxes is not None:
                self._move_point(event)

    def button_release_callback(self, event):
        """ Whenever a mouse button is released
        """
        if self._move_set == True and event.button == 1:
            # Finish moving set
            self._finish_move_set(event)
        elif self._move_set == True and event.button == 3:
            # Finish rotating set
            self._finish_rotate_set(event)
        elif event.button == 1:
            # Don't care if the cursor is on the plot (event.inaxes):
            self._finish_move_point()

    def key_press_callback(self, event):
        """ Whenever a key is pressed
        """
        print ('Pressed %s' % event.key)
        if self._quit_pending:
            if event.key == 'q':
                print ('Quitting.')
                exit(0)
            else:
                self._quit_pending = False
                print ('Quit cancelled.')
        elif event.key == 'q':
            print ('Quit requested.  Press q again to quit, any other key to cancel.')
            self._quit_pending = True
        elif event.key == 'm':
            self._print_keymap()
        elif event.inaxes:
            if event.key == 'c':
                self._regions.plot_curves()
                self._set_xlim_ylim()
                self._canvas.draw()
            elif event.key == 't':
                self._regions.toggle_original_line_visibility()
                self._canvas.draw()
            elif event.key == 'b':
                self._regions.toggle_points()
                self._canvas.draw()
            elif event.key == 'w':
                self._regions.write_output_files()
            elif event.key == 'y':
                self._regions.toggle_allow_xy_move()
            elif event.key == 'a':
                # Toggle adding/removing points
                if self._args.in_eos_file_base is not None:
                    # Don't allow for multi line data
                    print ("Option Disabled: Not a 1d Plot")
                else:
                    self._add_points = not self._add_points
                    if self._add_points == True:
                        print ("Add/Remove Mode Enabled")
                    else:
                        print ("Add/Remove Mode Disabled")
            elif event.key == 'e':
                #Toggle moving block of points
                if self._move_set == True:
                    self._regions.cancel_move_set()
                    self._move_set = False
                    self._find_set = False
                    self._last_click = None
                    self._canvas.draw()
                    self._selector.set_active(False)
                else:
                    self._find_set = not self._find_set
                    if self._find_set == True:
                        self._selector.set_active(True)
                        print ("Block Select Enabled")
                    else:
                        self._selector.set_active(False)
                        print ("Block Select Disabled")
            elif event.key == 'Q':
                self._get_equation()
            elif event.key == 'u':
                #Undo the last point manipulation
                self._regions.undo()
                self._canvas.draw()
            elif event.key == 'Z':
                if self._move_set == True:
                    self._regions.smooth_data("trilocal", self._xmin, self._xmax, 
                                              self._ymin, self._ymax)
                    self._canvas.draw()
                else:
                    print ("Select a region to smooth by pressing 'e'.")
            elif event.key == 'X':
                if self._move_set == True:
                    self._regions.smooth_data("integral", self._xmin, self._xmax, 
                                              self._ymin, self._ymax)
                    self._canvas.draw()
                else:
                    print ("Select a region to smooth by pressing 'e'.")
            elif event.key == 'C':
                if self._move_set == True:
                    self._regions.smooth_data("triintegral", self._xmin, self._xmax, 
                                              self._ymin, self._ymax)
                    self._canvas.draw()
                else:
                    print ("Select a region to smooth by pressing 'e'.")
            elif event.key == 'V':
                if self._move_set == True:
                    self._regions.smooth_data("acute", self._xmin, self._xmax, 
                                              self._ymin, self._ymax)
                    self._canvas.draw()
                else:
                    print ("Select a region to smooth by pressing 'e'.")
            elif event.key == 'k' or event.key == 'L':
                # Set graph x scale
                if self._ax.get_xscale() == 'linear':
                    self._ax.set_xscale('log')
                else:
                    self._ax.set_xscale('linear')                    
                self._set_xlim_ylim()
                self._set_logscale()
                self._canvas.draw()
            elif event.key == 'l':
                # Set graph y scale
                if self._ax.get_yscale() == 'linear':
                    self._ax.set_yscale('log')
                else:
                    self._ax.set_yscale('linear')
                self._set_xlim_ylim()
                self._canvas.draw()
            elif event.key == 'H':
                #Increase background point marker size
                for i in range(len(self._background_line)):
                    self._background_line[i].set_marker_size(self._background_line[i].get_marker_size() * 1.25)
                self._canvas.draw()
            elif event.key == 'J':
                #Decrease background point marker size
                for i in range(len(self._background_line)):
                    self._background_line[i].set_marker_size(self._background_line[i].get_marker_size() * 0.8)
                self._canvas.draw()

    def xlim_changed_callback(self, event):
        """ xlim is changed by a zoom or a pan

        :param event: matplotlib.axes.Axes object
        """
        #print('xlim changed.  New limit is: %.15E, %.15E' % event.get_xlim())
        self._regions.plot_curves()
        # Plot interactive equation
        self._plot_icurves()
        self._canvas.draw()

    def line_select_callback(self, eclick, erelease):
        """Press and release events for block selecting.
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        self._selector.set_active(False)
        self._find_set = False
        self._move_set = True

        self._xmin = min(x1, x2)
        self._xmax = max(x1, x2)
        self._ymin = min(y1, y2)
        self._ymax = max(y1, y2)

        self._regions.attempt_get_set(self._xmin, self._xmax, self._ymin, 
                                      self._ymax, eclick.x, eclick.y)  
        self._canvas.draw()

    def _register_callbacks(self):
        self._canvas.mpl_connect('draw_event', self.draw_callback)
        self._canvas.mpl_connect('button_press_event', self.button_press_callback)
        self._canvas.mpl_connect('button_release_event', self.button_release_callback)
        self._canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self._canvas.mpl_connect('key_press_event', self.key_press_callback)
        self._canvas.mpl_connect('line_select_event', self.line_select_callback)
        self._xlim_callback_id = self._ax.callbacks.connect(
                                'xlim_changed',         self.xlim_changed_callback)

    # END Callbacks ============================================================

    # Callback support: ========================================================

    def _draw(self):
        """ Called by self.motion_notify_callback and self.draw_callback
        """
        # Erases all but the latest line:
        self._canvas.restore_region(self._background)
        for line in self._background_line:
            line.draw()
        self._regions.draw()

        # TODO: Find out why using this gets a runtime error:
        # self._ax.draw()
        # THIS does infinite recursion:
        # self.canvas.draw()
        # Leaves old point there until button is released:
        self._canvas.blit(self._ax.bbox)

    def _attempt_begin_move_point(self, event):
        """
        Start moving the point nearest the event point if it is
        within EPSILON tolerance.
        """
        # Two button presses CAN happen without an intervening release.
        # Be sure we start clean:
        self._regions.cancel_any_move_points()
        # event.x and y are in display coords.  Get x in data coords to find
        # the right region:
        # data_x = self.original_line.point_to_data_space(event.x, event.y)[0]
        self._regions.attempt_begin_move_point(event)

    def _finish_move_set(self, event):
        """
        Redraw lines and setup click boundaries for next move
        """
        # Check that mouse is within canvas
        self._check_in_canvas(event)

        if event.xdata != None and event.ydata != None:
            xdiff = event.xdata - self._last_click[0]
            ydiff = event.ydata - self._last_click[1]
            
            self._xmin = self._xmin + xdiff
            self._xmax = self._xmax + xdiff
            self._ymin = self._ymin + ydiff
            self._ymax = self._ymax + ydiff

        self._regions.finish_move_set()
        self._canvas.draw()

    def _finish_rotate_set(self, event):
        """
        Same as finish_move_set but don't set new boundaries
        """
        self._regions.finish_move_set()
        self._canvas.draw()

    def _check_in_canvas(self, event):
        """
        Check if event happened in canvas. 
        If not, set event coordinates in canvas
        """
        bbox = self._ax.get_window_extent().transformed(self._figure.dpi_scale_trans.inverted())
        width, height = bbox.width * self._figure.dpi, bbox.height * self._figure.dpi

        if event.x > width:
            event.x = width
        if event.x < 0:
            event.x = 0
        if event.y > height:
            event.y = height
        if event.y < 0:
            event.y = 0

    def _check_move_set(self, event):
        """
        move set if user clicked in range
        """
        if (self._xmin < event.xdata < self._xmax and 
            self._ymin < event.ydata < self._ymax):
            self._last_click = [event.xdata, event.ydata]
            self._regions.begin_move_set(event)
        else:
            self._regions.cancel_move_set()
            self._move_set = False
            self._last_click = None
            self._canvas.draw()

    def _check_rotate_set(self, event):
        """
        Prepare to rotate set
        """
        self._last_click = [event.xdata, event.ydata]
        self._regions.begin_move_set(event)

    def _rotate_set(self, event):
        """
        rotate set
        """
        if self._regions._moving_set_region_index != None:
            self._regions.rotate_set(event, self._xmin, self._xmax, self._ymin, self._ymax)
            self._draw()

    def _move_point(self, event):
        """ Move the point in the data and redraw the line.
        """
        if self._regions.move_point_in_progress():
            self._regions.move_point(event)
            self._draw()

    def _attempt_move_set(self, event):
        """
        Prepare to move set.
        """
        if self._regions._moving_set_region_index != None:
            self._regions.move_set(event)
            self._draw()

    def _finish_move_point(self):
        """ Recalculate and redraw the fit curve.
        """
        # Two button releases CAN happen without an intervening push.  Be sure
        # there is a move to finish:
        if self._regions.move_point_in_progress():
            self._regions.finish_move_point()
            self._draw()

    def _set_logscale(self):
        """
        Set flag to determine if graph is logscale x
        """
        for region in self._regions._regions:
            for lineset in region._line_sets._sets.itervalues():
                lineset._logscale = not lineset._logscale

    def _get_equation(self):
        """
        Set up interface for plotting equations interactively
        """
        def plot_callback(e):
            # Button to plot equation
            if '_' in e.get():
                # Underscores are a possible security risk
                print("Underscores not allowed in equation.")
            elif e.get() not in self._iplot:
                self._iplot.update({e.get():self._ax.plot([],[],'--')[0]})
                # Disallow changing equation until it's been deleted
                e.configure(state='disabled')
                self._plot_icurves()
                self._canvas.draw()
                
        def delete_callback(e):
            # Button to delete equation
            if e.get() in self._iplot:
                self._iplot[e.get()].remove()
                self._iplot.pop(e.get(), None)
                # Allow changing equation once the old one has been deleted
                e.configure(state='normal')
                self._canvas.draw()

        def write_callback(e, e2):
            # Button to write file
            if e2.get() == '':
                print("Cannot Write: Please specify a filename.")
            else:
                try:
                    io.write_point_file(e2.get(), self._iplot[e.get()].get_xydata().tolist())
                except:
                    print("Cannot Write: Please make sure the equation is plotted.")

        def close_callback(e):
            # Delete line when closing window
            delete_callback(e)
            textbox.destroy()

        # Set up Equation Plotting window
        textbox = Tk()
        textbox.title("Interactive Plotter")
        l = Label(textbox, text="Enter an equation to plot in terms of x.")
        l2 = Label(textbox, text="F(x) = ")
        l3 = Label(textbox, text="Filename: ")
        e = Entry(textbox, width=35)
        e2 = Entry(textbox, width=35)
        b1 = Button(textbox, text="Plot", width=10, 
                   command=lambda:plot_callback(e))
        b2 = Button(textbox, text="Delete", width=10, 
                   command=lambda:delete_callback(e))
        b3 = Button(textbox, text="Write", width=10,
                    command=lambda:write_callback(e, e2))
        textbox.protocol("WM_DELETE_WINDOW", lambda:close_callback(e))
        l.grid(columnspan=4, pady=10, padx=100)
        l2.grid(column=0, row=1)
        l3.grid(column=0, row=2)
        e.grid(column=1, columnspan=2, row=1)
        e2.grid(column=1, columnspan=2, row=2, pady=5)
        b1.grid(column=0, pady=10)
        b2.grid(column=1, row=3, pady=10, padx=34)
        b3.grid(column=2, row=3)
        
    def _plot_icurves(self):
        """
        Plot interactive user inputted equations.
        Not in region file because this curve
        is plotted wherever the user's screen
        is.
        """
        xlimit = self._ax.get_xlim()
        xdata = []
        for i in range(0, self._args.points_in_user_curve):
            # Calculate N evenly spaced x values between screen limits
            xdata.append(xlimit[0]+(((abs(xlimit[1] - xlimit[0]))/
                                     (self._args.points_in_user_curve-1)) * i))
        for eq, plot in self._iplot.items():
            # Plot every user inputted equation
            ydata=[]
            # Make eval relatively safe by restricting namespace
            ns = vars(math).copy()
            ns['__builtins__'] = None
            for i in range(0, self._args.points_in_user_curve):
                # Evaluate the function for each x value
                x = xdata[i]
                ns['x'] = x
                ydata.append(eval(eq, ns))
            # Plot the data
            plot.set_data([[xdata],[ydata]])

    # END Callback support======================================================

    @staticmethod
    def _print_keymap():
        print('===============================================================')
        print('matplotlib keys:')
        print('===============================================================')
        def keys_for(action):
            return str(rcParams['keymap.%s' % action])
        for action in ('fullscreen', 'home', 'back', 'forward', 'pan', 'zoom',
                      'save', 'quit',):
            print(action + ': %s' % keys_for(action))
        print('show_all_axes: %s' % keys_for('all_axes'))
        print()
        print('Drag points to update line')
        print('Press k to toggle x log scale')
        print('Press l to toggle y log scale') 
        print('Press q then q again to quit')
        print('Press t to toggle original line on and off [default: off]')
        print('Press b to toggle points on and off [default: on]')
        print('Press w to write the the current points to a file')
        print('Press y to toggle xy move capability on and off')
        print('Press a to toggle adding and removing points with left and right click \n[default: off]')
        print('Press e to toggle selecting a block of points')
        print('Press u to undo the last point manipulation')
        print('Press <shift> H to increase size of background markers')
        print('Press <shift> J to decrease size of background markers')
        print('Press <shift> Q to enter an equation to plot')
        print('Press <shift> Z for trilocal smoothing')
        print('Press <shift> X for integral smoothing')
#        print('Press <shift> C for trintegral smoothing')
        print('Press  <shift> V for acute repair smoothing')
        print()
        print('===============================================================')
        print('Make sure focus is on the plotting window and the cursor is')
        print('also in the plotting window when using these keys.')
        print()
        print('Press "m" to show these keys again')
        print('===============================================================')
        print()

if __name__ == '__main__':

    CurveInteractor().run()

def main():
    CurveInteractor().run()

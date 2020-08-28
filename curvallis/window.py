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

"""
Windowing system for curvallis
@author: Eric Heinke (sudo-Eric)
"""

from curvallis.version import version as VERSION_STRING
from tkinter import ttk
import tkinter

show_window_errors = False


def define_args(parser):
    parser.add_argument("--show_window_errors", help="shows windowing errors when they occur", action="store_true")


def process_args(args):
    global show_window_errors
    show_window_errors = args.show_window_errors


class WindowedDisplay(object):
    """WindowedDisplay generates a display window that shows blocks of information
    separated by a horizontal line and optionally a header."""

    def __init__(self, info_blocks, block_headers, window_name="WindowedDisplay", resizeable=None,
                 max_window_width=0.4, max_window_height=0.8, window_size_errors=None, scrollbar_width=20):
        if resizeable is None:
            resizeable = [False, False]
        if window_size_errors is None:
            window_size_errors = [False, False]
        self._main_window_open = False
        self._info_blocks = info_blocks
        self._block_headers = block_headers
        self._resizeable = resizeable
        self._window_name = window_name
        self.max_window_width = max_window_width
        self.max_window_height = max_window_height
        self.window_size_errors = window_size_errors
        self.scrollbar_width = scrollbar_width  # not fully working
        self._main_window = None

    def display_main_window(self):
        def on_closing():
            self._main_window_open = False
            self._main_window.destroy()

        if self._main_window_open:
            self.bring_forward()
            return

        if not (len(self._block_headers) == len(self._info_blocks)):
            self.window_error("Text length mismatch error",
                              "The lengths of the text headers and the text blocks do not match", True)

        # Calculate size of horizontal dividing line
        horizontal_line = '=' * self._get_longest_line_width(2)
        # Font and font size (MUST use fixed width font / constant width font)
        text_font = ("Courier", 12)
        # Create Textbox
        self._main_window = tkinter.Tk()
        # resize window to fit content
        self.resize_to_content()
        self.set_resizability(self._resizeable)
        # Set function for when window is closed
        self._main_window.protocol("WM_DELETE_WINDOW", on_closing)
        # Set title of window
        self.set_window_name(self._window_name)
        # Generate scrollbar
        ttk.Style().configure("Vertical.TScrollbar", arrowsize=self.scrollbar_width)
        container = ttk.Frame(self._main_window)
        canvas = tkinter.Canvas(container,
                                width=self._main_window.winfo_reqwidth(), height=self._main_window.winfo_reqheight())
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview, cursor="sb_v_double_arrow",
                                  style="Vertical.TScrollbar")
        scrollable_frame = ttk.Frame(canvas)
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        # Generate textboxes and place them in window
        for i in range(len(self._info_blocks)):
            block = self._info_blocks[i]
            header = self._block_headers[i]
            if (len(header) > 0):
                tkinter.Label(scrollable_frame, text=horizontal_line + "\n" + header + "\n" + horizontal_line,
                              font=text_font).pack()
            else:
                tkinter.Label(scrollable_frame, text=horizontal_line, font=text_font).pack()
            self._place_lines(scrollable_frame, block, text_font)
        tkinter.Label(scrollable_frame, text=horizontal_line, font=text_font).pack()
        # Place items in window
        container.pack(fill="both", expand=True)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        # Declare that the help window has been created
        self._main_window_open = True

    ##################################################
    def _place_lines(self, window, lines, text_font):
        for i in range(len(lines)):
            tkinter.Label(window, text=lines[i], font=text_font).pack(anchor='w')

    def _number_of_lines(self):
        total_lines = 0
        for block in self._info_blocks:
            total_lines += len(block)
        return total_lines

    def _get_header_space(self):
        header_space = len(self._block_headers)
        for i in self._block_headers:
            if (len(i) > 0):
                header_space += 2
        return header_space + 1

    def _get_longest_line_width(self, extra_spacing):
        max_length = 0
        for block in self._info_blocks:
            for line in block:
                if (len(line) > max_length):
                    max_length = len(line)
        for line in self._block_headers:
            if (len(line) > max_length):
                max_length = len(line)
        max_length += extra_spacing
        return max_length

    ##################################################

    def get_window_name(self):
        return self._window_name

    def set_window_name(self, new_name):
        # Set title of window
        self._window_name = new_name
        self._main_window.title(self._window_name)

    def set_resizability(self, resizability):
        # Enable / Disable window resizing (Horizontal, Vertical)
        self._resizeable = resizability
        self._main_window.resizable(self._resizeable[0], self._resizeable[1])

    def get_block_headers(self):
        return self._block_headers

    def set_block_headers(self, new_block_headers):
        self._block_headers = new_block_headers

    def get_info_blocks(self):
        return self._info_blocks

    def set_info_blocks(self, new_info_blocks):
        self._info_blocks = new_info_blocks

    def set_window_size(self, width, height, x_offset=30, y_offset=30):  # Might also need to resize canvases
        # width x height + x_offset + y_offset
        self._main_window.geometry("%dx%d+%d+%d" % (width, height, x_offset, y_offset))

    def resize_to_content(self):
        window_width = 10 * self._get_longest_line_width(2)
        window_height = (self._get_header_space() + self._number_of_lines()) * 22.5
        # Host system screen resolution
        host_screen = (self._main_window.winfo_screenwidth(), self._main_window.winfo_screenheight())
        # Check if window can fit on screen (40% of width, 80% of height)
        if window_width > (host_screen[0] * self.max_window_width):
            window_width = int(host_screen[0] * self.max_window_width)
            self.window_error("Window too wide", self._window_name + " window exceeds " +
                              str(self.max_window_width * 100) + "% of screen.", self.window_size_errors[0])
        if window_height > (host_screen[1] * self.max_window_height):
            window_height = int(host_screen[1] * self.max_window_height)
            self.window_error("Window too tall", self._window_name + " window exceeds " +
                              str(self.max_window_height * 100) + "% of screen.", self.window_size_errors[1])
        self.set_window_size(window_width + self.scrollbar_width, window_height)

    def redraw(self):
        if self._main_window_open:
            self._main_window.destroy()
            self._main_window_open = False
            self.display_main_window()
        else:
            self.window_error("Redraw error", "Unable to redraw window that does not exist.")

    def bring_forward(self):
        if self._main_window_open:
            self._main_window.lift()
            self._main_window.focus_force()
        else:
            self.window_error("Window action error", "Unable to bring nonexistent window to front.")

    def is_open(self):
        return self._main_window_open

    def window_error(self, error, error_message, fatal=False):
        if show_window_errors:
            print("Error: " + error)
            print(error_message)
        if fatal:
            self._main_window.destroy()


########################################################################################################################
def key_mappings_window_header():
    return ["Matplotlib Keys:",
            "Curvallis Keys:",
            ""]


def key_mappings_window_text():
    return [
        [  # MatPlotLib Keys
            'Press h or r to reset to home view',
            'Press p to toggle pan/zoom mode [left mouse: pan, right mouse: zoom]',
            'Press o to toggle Zoom-to-rectangle',
            'Press <ctrl> s to save a screenshot of the MarPlotLib display',
            'Press f or <ctrl> f to toggle fullscreen',
            'Press <ctrl> w to immediately quit',
            # 'Constrain pan/zoom to x axis	hold x when panning/zooming with mouse',
            # 'Constrain pan/zoom to y axis	hold y when panning/zooming with mouse',
            # 'Preserve aspect ratio	hold CONTROL when panning/zooming with mouse',
            'Press g to toggle major grids when mouse is over plot',
            'Press <shift> G to toggle minor grids when mouse is over plot',
            'Press <shift> L or k to toggle x axis scale (log/linear)',
            '\twhen mouse is over plot',
            # 'Press l to toggle y axis scale (log/linear) when mouse is over plot'
        ], [  # Curvallis Keys
            'Press q twice to quit',
            'Press w to write the the current points to a file',
            'Press u to undo the last point manipulation',
            'Drag points to update line',
            'Press t to toggle original line on and off [default: off]',
            'Press b to toggle points on and off [default: on]',
            'Press y to toggle xy move capability on and off',
            'Press a to toggle adding and removing points with left ',
            '\tand right click respectively [default: off]',
            'Press e to toggle selecting a block of points',
            'Press i to enlarge figure margins',
            'Press <shift> I to decrease figure margins',
            'Press <shift> H to increase size of background markers',
            'Press <shift> J to decrease size of background markers',
            'Press <shift> Q to enter an equation to plot',
            'Press <shift> Z for trilocal smoothing',
            'Press <shift> X for integral smoothing',
            'Press <shift> B for B-spline smoothing',
            'Press <shift> V for acute repair smoothing'
        ], [  # Post Keymapping Lines
            'Make sure focus is on the plotting window and the cursor is',
            'also in the plotting window when using these keys.',
            '',
            'Press "F1" to show these keys again',
            '',
            'More key mappings can be found at:',
            'https://github.com/LLNL/Curvallis#interactive-commands',
        ]]


key_mappings_window = WindowedDisplay(key_mappings_window_text(), key_mappings_window_header(),
                                      "Key Mappings", [False, True], 0.40, 0.80, [True, False])

fitter_info_window = WindowedDisplay([], [],
                                     "Fitter Information", [False, True], 0.40, 0.80, [False, False])

########################################################################################################################
fitter_info_window_working_index = -1

# This function is used to update the fitter info window in a more streamline manner.
def update_fitter_info_window(index, reset_text, new_text):
    """ Updates the contents of the fitter info window. Updated contents
        are not shown until the window is refreshed.
    """
    # index = -1: Use the same index as was last used
    # index = -2: Set index to the largest valid index and clear block entirely
    global fitter_info_window_working_index
    if index == -1:
        index = fitter_info_window_working_index
    elif index < 0:
        pass
    else:
        fitter_info_window_working_index = index
    info_blocks = fitter_info_window.get_info_blocks()
    if index == -1 or index < -2 or index >= len(info_blocks):
        fitter_info_window.window_error("Display Variable Update Error",
                                        "An attempt was made to access a display variable out of bounds.")
        return
    if (index >= 0):
        if reset_text:
            info_blocks[index] = [info_blocks[index][0]]
    else:
        index = -1
        if reset_text:
            info_blocks[index] = []
    lines = new_text.split('\n')
    for line in lines:
        info_blocks[index].append(line)

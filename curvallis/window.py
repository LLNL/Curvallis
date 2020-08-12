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


class WindowedDisplay(object):
    """WindowedDisplay generates a display window that shows blocks of information
    separated by a horizontal line and optionally a header."""

    def __init__(self, info_blocks, block_headers, window_name="WindowedDisplay", resizeable=None,
                 max_window_width=0.4, max_window_height=0.8, window_size_errors=None):
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
        self._main_window = None

    def display_main_window(self):
        def on_closing():
            self._main_window_open = False
            self._main_window.destroy()

        if (self._main_window_open):
            self._main_window.lift()
            self._main_window.focus_force()
            return

        ##################################################
        def place_lines(window, lines, text_font):
            for i in range(len(lines)):
                tkinter.Label(window, text=lines[i], font=text_font).pack(anchor='w')

        def longest_number_of_lines():
            total_lines = 0
            for block in self._info_blocks:
                total_lines += len(block)
            total_lines -= len(self._info_blocks) - 1
            return total_lines

        def get_header_space():
            header_space = 0
            for i in self._block_headers:
                if (len(i) > 0):
                    header_space += 1
            return header_space + 1

        def get_longest_line_width(extra_spacing):
            max_length = 0
            for block in self._info_blocks:
                for line in block:
                    if (len(line) > max_length):
                        max_length = len(line)
            max_length += extra_spacing
            return max_length

        ##################################################
        # Calculate size of horizontal dividing line
        horizontal_line = '=' * get_longest_line_width(2)
        # Font and font size (MUST use fixed width font / constant width font)
        text_font = ("Courier", 12)
        # Scrollbar width
        scrollbar_width = 20  # not fully working
        # Window Size (width, height, x_offset, y_offset)
        window_size = (
            (10 * get_longest_line_width(2)), ((get_header_space() + longest_number_of_lines()) * 22.5), 30, 30)
        # Create Textbox
        self._main_window = tkinter.Tk()
        # Host system screen resolution
        host_screen = (self._main_window.winfo_screenwidth(), self._main_window.winfo_screenheight())
        # Check if window can fit on screen (40% of width, 80% of height)
        if (window_size[0] > (host_screen[0] * self.max_window_width)):
            window_size = ((host_screen[0] * self.max_window_width), host_screen[1], window_size[2], window_size[3])
            if(self.window_size_errors[0]):
                print("Error: window too wide")
                print(self._window_name + " window exceeds " + str(self.max_window_width*100) + "% of screen.")
                self._main_window.destroy()
                return
        if (window_size[1] > (host_screen[1] * self.max_window_height)):
            window_size = (window_size[0], int(host_screen[1] * self.max_window_height), window_size[2], window_size[3])
            if(self.window_size_errors[1]):
                print("Error: window too tall")
                print(self._window_name + " window exceeds " + str(self.max_window_height*100) + "% of screen.")
                self._main_window.destroy()
                return
        # width x height + x_offset + y_offset:
        self._main_window.geometry(
            "%dx%d+%d+%d" % (window_size[0] + scrollbar_width, window_size[1], window_size[2], window_size[3]))
        # Enable / Disable window resizing (Horizontal, Vertical)
        self._main_window.resizable(self._resizeable[0], self._resizeable[1])
        # Set function for when window is closed
        self._main_window.protocol("WM_DELETE_WINDOW", on_closing)
        # Set title of window
        self._main_window.title(self._window_name)
        # Generate scrollbar
        ttk.Style().configure("Vertical.TScrollbar", arrowsize=scrollbar_width)
        container = ttk.Frame(self._main_window)
        canvas = tkinter.Canvas(container, width=window_size[0], height=window_size[1])
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
            place_lines(scrollable_frame, block, text_font)
        tkinter.Label(scrollable_frame, text=horizontal_line, font=text_font).pack()
        # Place items in window
        container.pack(fill="both", expand=True)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        # Declare that the help window has been created
        self._main_window_open = True

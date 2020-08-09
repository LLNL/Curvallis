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
Windowing system for curvallis
@author: Eric Heinke (sudo-Eric)
"""

from curvallis.version import version as VERSION_STRING
from tkinter import ttk
import tkinter

fitting_window_open = False

def _generate_window():
    """
    Set up interface for displaying fitting information
    """
    global fitting_window_open

    def on_closing():
        global fitting_window_open
        fitting_window_open = False
        main_window.destroy()

    if (fitting_window_open):
        global main_window
        main_window.lift()
        main_window.focus_force()
        return

    info_blocks = [
        [
            'Fitter type: none',
            'Moving point: region 0, \nat: 1.506464681207422E+01, -8.986979761105912E+12',
            'to: 1.458308225514779E+01, 7.536867438430522E+14'
        ], [
            'Fitter type: poly3',
            'Moving point: region 1, \nat: 6.870196867832244E+01, 4.190320426764570E+13',
            'to: 6.705440132823236E+01, 6.569662922040228E+14',
            'Calculated polynomial is:',
            '-1.010E+09 * x**3 +1.678E+11 * x**2 -5.598E+12 * x**1 +4.349E+13'
        ], [
            'Fitter type: ebirch3',
            'Moving point: region 2, \nat: 1.189499107923780E+02, 1.186430432031626E+14',
            'to: 1.622801285349414E+02, 9.084394664654992E+14',
            'B0 = 479956838652108.9;',
            'Bp = 4.012249080480707;',
            'rho0 = 19.848030178310996;',
            'E0 = 1952.2421183931094;'
        ]]

    ##################################################
    def placelines(window, lines, text_font):
        for i in range(len(lines)):
            l = tkinter.Label(window, text=lines[i], font=text_font).pack(anchor='w')

    def longest_number_of_lines():
        total_lines = 0
        for block in info_blocks:
            total_lines += len(block)
        return total_lines

    def get_longest_line_width(extra_spacing):
        # Need to add check for lines containing '\n'
        max_length = 0
        for block in info_blocks:
            for line in block:
                if (len(line) > max_length):
                    max_length = len(line)
        max_length += extra_spacing
        return max_length

    ##################################################
    # Calculate size of vertical deviding line
    vertical_line = '|'
    vertical_line += '\n|' * (longest_number_of_lines()-1)
    # Font and font size (MUST use fixed width font / constant width font)
    text_font = ("Courier", 12)
    # Scrollbar width
    scrollbar_width = 20  # not fully working
    # Calculate the width of the window
    window_width = 10 * get_longest_line_width(2)
    # Window Size (width, height, x_offset, y_offset)
    window_size = (window_width, ((1 + longest_number_of_lines()) * 22.5), 30, 30)
    # Create Textbox
    main_window = tkinter.Tk()
    # Host system screen resolution
    host_screen = (main_window.winfo_screenwidth(), main_window.winfo_screenheight())
    # Check if window can fit on screen (40% of width, 80% of height)
    if (window_size[0] > (host_screen[0] * 0.8)):
        window_size = ((host_screen[0] * 0.8), host_screen[1], window_size[2], window_size[3])
        print("error: window too wide\nkey mappings window excedes 40% of screen")
        return
    if (window_size[1] > (host_screen[1] * 0.5)):
        window_size = (window_size[0], int(host_screen[1] * 0.5), window_size[2], window_size[3])
    # width x height + x_offset + y_offset:
    main_window.geometry(
        "%dx%d+%d+%d" % (window_size[0] + scrollbar_width, window_size[1], window_size[2], window_size[3]))
    # Enable / Disable window resizing (Horizontal, Vertical)
    main_window.resizable(True, True)
    # Set function for when window is closed
    main_window.protocol("WM_DELETE_WINDOW", on_closing)
    # Set title of window
    main_window.title("Fitter Info")
    # Generate scrollbar
    style = ttk.Style().configure("Horizontal.TScrollbar", arrowsize=scrollbar_width)
    container = ttk.Frame(main_window)
    canvas = tkinter.Canvas(container, width=window_size[0], height=window_size[1])
    scrollbar = ttk.Scrollbar(container, orient="horizontal", command=canvas.xview, cursor="sb_h_double_arrow",
                              style="Horizontal.TScrollbar")
    scrollable_frame = ttk.Frame(canvas)
    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(
            scrollregion=canvas.bbox("all")
        )
    )
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(xscrollcommand=scrollbar.set)
    # Generate textboxes and place them in window
    lines_from_top = 0
    for i in range(len(info_blocks)):
        block = info_blocks[i]
        t = tkinter.Label(scrollable_frame, text=vertical_line + "\n" + block[0] + "\n" + vertical_line,
                          font=text_font).pack()
        placelines(scrollable_frame, block[1:len(block)], text_font)
        lines_from_top += len(block) + 2
    '''block = info_blocks[-1]
    t = tkinter.Label(scrollable_frame, text=vertical_line, font=text_font).pack()
    placelines(scrollable_frame, block, text_font)
    t = tkinter.Label(scrollable_frame, text=vertical_line, font=text_font).pack()'''
    # Place items in window
    container.pack(fill="both", expand=True)
    canvas.pack(side="top", fill="both", expand=True)
    scrollbar.pack(side="bottom", fill="x")
    # Declare that the help window has been created
    fitting_window_open = True

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

import numpy as np
from operator import itemgetter
import math
from curvallis.version import version as VERSION_STRING

# These line attributes are declared here at the top of the file so they are
# easy to find and change:
line_attributes = dict(
    background_points=dict(
        linestyle='',                 marker='+', markerfacecolor='green',
        markersize=20),
    derivative=dict(
        linestyle='-', color='orange',),
    fit_curve=dict(
        linestyle='-', color='red',),
    integral=dict(
        linestyle='-', color='purple',),
    movable=dict(
        linestyle='-', color='blue',  marker='o', markerfacecolor='yellow'),
    movable_no_curves=dict(
        linestyle='-',                marker='o', markerfacecolor='yellow'),
    original=dict(
        linestyle='-', color='black', marker='o', markerfacecolor='white'),
    original_no_curves=dict(
        linestyle='',                 marker='o', markerfacecolor='orange'),
    region_boundary=dict(
        linestyle=':', color='black',))


class Line(object):
    """ Contains info about one line.
    """
    def __init__(self, ax, attributes):
        """
        :param ax:         matplotlib axes in a figure
        :param attributes: line attributes dictionary
        """
        # ax is the figure to draw this on.
        self._ax = ax
        self._attributes = attributes
        self._id = None
        # Plots moved points in a different color.
        self._moved_points = None
        # Plots highlighted points from block select in a different color.
        self._highlight = None
        # Stores last data information for undo command
        self._last_data = None
        self._last_highlight = None

    def draw(self):
        """ Actually draws the line on the figure.
        """
        if self._id is not None:
            self._ax.draw_artist(self._id)
        if self._highlight is not None:
            self._ax.draw_artist(self._highlight)
        if self._moved_points is not None:
            self._ax.draw_artist(self._moved_points)

    def get_x_data_y_data(self):
        """
        :return:  [[x1, x2, x3...],[y1, y2, y3...]]
        """
        x_data_y_data = [[data[0] for data in self.get_xy_data()],
                         [data[1] for data in self.get_xy_data()]]
        return x_data_y_data

    def get_xy_data(self):
        """
        :return: [[x1, y1],[x2, y2],[x3, y3]...]
        """
        xy_data = self._id.get_xydata().tolist()
        return xy_data

    def get_point_count(self):
        return len(self.get_xy_data())

    def get_marker_size(self):
        return self._id.get_markersize()

    def point_to_data_space(self, x, y):
        """ Convert a point in display coordinates, as in an callback event, to
        a point in data coordinates.
        :return: x, y
        """
        invert = self._id.get_transform().inverted()
        return invert.transform_point((x, y))

    def points_to_display_space(self, points):
        """ Convert a list of points in data coordinates to display coordinates.

        :param points: list of x,y pairs
        :return:
        """
        return self._id.get_transform().transform(points)

    def plot_x_data_y_data(self, x_data_y_data, visible=True, animated=True):
        """ Creates a line, plots it, and stores its id

        Calling this from _plot_fit_curve repeatedly leaves the old curve
        lying around (line doesn't get finalized?) with BOTH the ax.plot
        and the ax.add_line approach.

        :param x_data_y_data:  [[x1, x2, x3...],[y1, y2, y3...]]
        :param visible:  boolean
        """
        x_data_y_data = list(x_data_y_data)
        assert len(x_data_y_data) == 2, "No points in region."

        #Actual Line
        self._id = self._ax.plot(
            x_data_y_data[0], x_data_y_data[1],
            visible=visible,
            animated=animated,
            **self._attributes)[0]
        #Set highlighted points (Starts blank)
        self._highlight = self._ax.plot(
            [],[],
            visible=visible,
            animated=animated,
            linestyle='',
            marker='o', markerfacecolor='green')[0]
        #Moved points (Starts blank)
        self._moved_points = self._ax.plot(
            [],[],
            visible=visible,
            animated=animated,
            linestyle='',
            marker='o', markerfacecolor='red')[0]

    def plot_xy_data(self, xy_data, visible=True, animated=True):
        """ See plot_x_data_y_data()

        :param xy_data: [[x1, y1],[x2, y2],[x3, y3]...]
        :param visible:  boolean
        :param animated:  boolean
        """
        x_data_y_data = [[data[0] for data in xy_data],
                         [data[1] for data in xy_data]]
        self.plot_x_data_y_data(x_data_y_data, visible, animated)

    def set_block_select(self, xy_data):
        """
        Highlight block select points in a different color
        """
        # Set last highlight for undo information
        self._last_highlight = self._highlight.get_xydata().tolist()
        if xy_data != []:
            self._highlight.set_data([data[0] for data in xy_data],
                                     [data[1] for data in xy_data])
        else:
            self._highlight.set_data([],[])            

    def set_block_select_no_undo(self, xy_data):
        """
        Same as set_block_select, but don't set undo information.
        """
        if xy_data != []:
            self._highlight.set_data([data[0] for data in xy_data],
                                     [data[1] for data in xy_data])
        else:
            self._highlight.set_data([],[])            

    def set_marker_size(self, pts):
        self._id.set_markersize(pts)

    def set_x_data_y_data(self, x_data_y_data):
        """
        :param x_data_y_data: [[x1, x2, x3...],[y1, y2, y3...]]
        """
        #Set undo information first
        self._last_data = self.get_x_data_y_data()
        #Set new line data
        self._id.set_data(x_data_y_data)
        #Disable moved point highlighting
        self._moved_points.set_data([[],[]])

    def set_x_data_y_data_no_undo(self, x_data_y_data):
        """
        Same as set_x_data_y_data except doesn't store undo information
        Used for moving point by hand while dragging the mouse.
        """
        #Set new line data
        self._id.set_data(x_data_y_data)
        #Disable moved points highlights
        self._moved_points.set_data([[],[]])

    def set_x_data_y_data_moved_points(self, x_data_y_data):
        """ Highlights points which have been moved since
            points have last been set.
            Used after automatic smoothing algorithms
        """
        # Get the old data
        old_points = self.get_x_data_y_data()
        # Store undo information
        self._last_data = old_points

        newpoints = []
        newpoints.append([])
        newpoints.append([])

        #Check for differences in new and old data
        old_points = list(old_points)
        x_data_y_data = list(x_data_y_data)
        for i in range(0, min(len(old_points[0]), len(x_data_y_data[0]))):
            if (old_points[0][i] != x_data_y_data[0][i] or
                old_points[1][i] != x_data_y_data[1][i]):
                newpoints[0].append(x_data_y_data[0][i])
                newpoints[1].append(x_data_y_data[1][i])

        self._id.set_data(x_data_y_data)
        self._moved_points.set_data(newpoints)

    def set_xy_data(self, xy_data):
        """
        :param xy_data: [[x1, y1],[x2, y2],[x3, y3]...]
        """
        self.set_x_data_y_data([[data[0] for data in xy_data],
                                [data[1] for data in xy_data]])

    def set_xy_data_no_undo(self, xy_data):
        self.set_x_data_y_data_no_undo([[data[0] for data in xy_data],
                                        [data[1] for data in xy_data]])

    def undo(self):
        """ Undo the last movement
        """
        if self._last_data != None:
            orig = self._last_data
            self.set_x_data_y_data(self._last_data)
            self._last_data = orig
            self.draw()

    def toggle_visible(self):
        self._id.set_visible(not self._id.get_visible())

    def replot(self):
        self.plot_x_data_y_data(self.get_x_data_y_data(), visible=True, animated=True)

class Line_Set(object):
    """Contains the original and movable lines, and point movement.
    """
    # Max pixel distance to count as "close" to a point when selecting to move:
    EPSILON = 5
    def __init__(self, name, ax, allow_xy_move, is_eos_data):
        """ Create the lines and give each some line attributes:
        See matplotlib.lines or matplotlib.axes for more attributes:
        """
        self._allow_xy_move = allow_xy_move
        self._is_eos_data = is_eos_data
        self._moving_point_index = None
        self._name = name
        self._in_set = False
        self._set_points = []
        self._set_index = []
        self._last_event = None
        if is_eos_data:
            self.movable  = Line(ax, line_attributes['movable_no_curves'])
            self.original = Line(ax, line_attributes['original_no_curves'])
        else:
            self.movable = Line(ax, line_attributes['movable'])
            self.original = Line(ax, line_attributes['original'])

    def draw(self):
        """ Redraw all the lines
        """
        for item in self.__dict__.values():
            if isinstance(item, Line):
                item.draw()

    def get_info(self, indent=''):
        result  = indent + 'Name: "%s"\n' % self._name
        result += indent + 'Point Count: %s\n' % self.movable.get_point_count()
        return result

    def get_movable_points(self):
        """ Return the movable line's points as a list of [x,y] pairs
        """
        return self.movable.get_xy_data()

    def get_name(self):
        return self._name

    def set_allow_xy_move(self, allow):
        self._allow_xy_move = allow

    def toggle_original_line_visibility(self):
        self.original.toggle_visible()

    def replot(self):
        self.movable.replot()

    # Point movement ###########################################################

    def attempt_begin_move_point(self, event):
        """ Record the index of the closest point to the event point, if it is
        also within EPSILON of the event point.  After success,
        move_point_in_progress returns True.
        """
        data_points = self.movable.get_xy_data()
        display_points = self.movable.points_to_display_space(data_points)
        # Could below also be x_values, y_values = zip(*display_points)? NO!
        # Below keeps x_values and x_values as an ndarray.  Zip makes tuples.
        x_values, y_values = display_points[:, 0], display_points[:, 1]
        distances = np.sqrt((x_values - event.x) ** 2 +
                            (y_values - event.y) ** 2)
        closest_index = distances.argmin()
        if distances[closest_index] <= self.EPSILON:
            self._moving_point_index = closest_index
            # Set undo information
            self.movable.set_xy_data(data_points)

    def attempt_get_set(self, xmin, xmax, ymin, ymax):
        data_points = self.movable.get_xy_data()
        self._set_index = []
        self._set_points = []
        i = 0
        for point in data_points:
            if xmin < point[0] < xmax and ymin < point[1] < ymax:
                self._in_set = True
                self._set_points.append(point)
                self._set_index.append(i)
            i += 1
        # Set undo information
        self.movable.set_xy_data(data_points)
        self.movable.set_block_select(self._set_points)

    def begin_move_set(self, event):
        # Set undo information
        self.movable.set_xy_data(self.movable.get_xy_data())
        self.movable.set_block_select(self._set_points)
        self._last_event = [event.x, event.y]
    def move_set(self, event):
        """ Move set of points in data and replot line
        """
        datapoints = self.movable.get_xy_data()
        if self._allow_xy_move == True:
            xchange = event.x - self._last_event[0]
        else:
            xchange = 0
        ychange = event.y - self._last_event[1]

        for i in range(0, len(self._set_index)):
            point = self.movable.points_to_display_space([datapoints[self._set_index[i]]])[0]
            datapoints[self._set_index[i]] = self.movable.point_to_data_space(
                (point[0] + xchange), (point[1] + ychange))
            self._set_points[i] = datapoints[self._set_index[i]]

        self.movable.set_xy_data_no_undo(datapoints)
        self.movable.set_block_select_no_undo(self._set_points)
        self._last_event = [event.x, event.y]

    def rotate_set(self, event, xmin, xmax, ymin, ymax):
        """ Rotate set of points in data and replot line
        """
        # Impossible to rotate and keep same x values
        if self._allow_xy_move == False:
            print ("Toggle xy move capability to rotate points")
            return

        datapoints = self.movable.get_xy_data()

        # Pivot point for rotation. 
        point1 = self.movable.points_to_display_space([[xmin + (xmax - xmin) / 2, ymin + (ymax - ymin) / 2]])[0]

        # Previous and current mouse drag positions
        point2 = self._last_event
        point3 = [event.x, event.y]

        vec1 = [(point2[0] - point1[0]), (point2[1] - point1[1])]
        vec2 = [(point3[0] - point1[0]), (point3[1] - point1[1])]

        # Calculate angle from previous and current mouse position
        angle = math.acos((vec1[0]*vec2[0] + vec1[1]*vec2[1]) / 
                          (math.sqrt(vec1[0]**2 + vec1[1]**2) * math.sqrt(vec2[0]**2 + vec2[1]**2)))
        if ((point3[1] < point2[1] and point3[0] > point1[0]) or
            (point3[1] > point2[1] and point3[0] < point1[0])):
            angle = angle * -1

        # Convert to display space, multipy point by rotation matrix, and convert back to data space
        for i in range(0, len(self._set_index)):
            point = self.movable.points_to_display_space([datapoints[self._set_index[i]]])[0]
            point[0] = point[0] - point1[0]
            point[1] = point[1] - point1[1]
            newx = point[0] * math.cos(angle) - point[1] * math.sin(angle)
            newy = point[0] * math.sin(angle) + point[1] * math.cos(angle)
            newx += point1[0]
            newy += point1[1]
            datapoints[self._set_index[i]] = self.movable.point_to_data_space(
                newx, newy)
            self._set_points[i] = datapoints[self._set_index[i]]

        self.movable.set_xy_data_no_undo(datapoints)
        self.movable.set_block_select_no_undo(self._set_points)
        self._last_event = [event.x, event.y]            

    def finish_move_set(self):
        self.sort_points_by_x()
        self._set_index = []
        data_points = self.movable.get_xy_data()
        i = 0
        j = 0
        # Recalculate set indexes if points are in different order
        for point in data_points:
            if np.array_equal(point, self._set_points[i]):
                self._set_index.append(j)
                if np.array_equal(self._set_points[i], self._set_points[-1]):
                    break
                i += 1
            j += 1

    def cancel_move_set(self):
        self._set_points = []
        self._set_index = []
        self._in_set = False
        self._last_event = None
        self.sort_points_by_x()
        self.movable.set_block_select([])

    def cancel_any_move_points(self):
        """ Reset any moving point indexes.
        """
        self._moving_point_index = None

    def finish_move_point(self):
        """ Recalculate and redraw the fit curve.
        """
        if self.move_point_in_progress():
            if self._allow_xy_move:
                self.sort_points_by_x()
            self._moving_point_index = None

    def sort_points_by_x(self):
        """ Order points by x so that moved line draws nicely and doesn't
        cross itself.
        """
        points = self.movable.get_xy_data()
        sorted_points = sorted(points, key=itemgetter(0))
        self.movable.set_xy_data_no_undo(sorted_points)

    def move_point(self, event):
        """ Move the point in the data and replot the line.
        """
        if self.move_point_in_progress():
            points = self.movable.get_xy_data()
            points[self._moving_point_index][1] = event.ydata
            if self._allow_xy_move:
                points[self._moving_point_index][0] = event.xdata
            self.movable.set_xy_data_no_undo(points)

    def move_point_in_progress(self):
        return self._moving_point_index is not None

    def moving_point_info(self, xy_only=False):
        if self.move_point_in_progress():
            x, y = self.movable.get_xy_data()[self._moving_point_index]
            if xy_only:
                if len(self._name.split(" ")) > 1:
                    #Print if 2d data
                    return '%s at (%.15E,%.15E) = %.15E' % \
                        (self._name.split(":")[0], x,
                         float(self._name.split(" ")[1]), y)
                else:
                    #Print if 1d data
                    return '%.15E, %.15E' % (x,y)
            else:
                if len(self._name.split(" ")) > 1:
                    #Print if 2d data
                    return '\nat: %s at (%.15E,%.15E) = %.15E' % \
                        (self._name.split(":")[0], x,
                         float(self._name.split(" ")[1]), y)
                else:
                    #Print if 1d data
                    return '\nat: %.15E, %.15E' % (x,y)
                
        else:
            return '(no point move in progress)'

    def add_point(self, event):
        done = False
        points = self.movable.get_xy_data()
        for i in range (0, len(points)):
            if points[i][0] > event.xdata:
                points.insert(i, [event.xdata,event.ydata])
                done = True
                break
        if done == False:
            points.append([event.xdata,event.ydata])
        self.movable.set_xy_data_no_undo(points)

    def remove_point(self, event):
        data_points = self.movable.get_xy_data()
        display_points = self.movable.points_to_display_space(data_points)
        x_values, y_values = display_points[:, 0], display_points[:, 1]
        distances = np.sqrt((x_values - event.x) ** 2 +
                            (y_values - event.y) ** 2)
        closest_index = distances.argmin()
        if distances[closest_index] <= self.EPSILON:
            remove_point_index = closest_index
            data_points.remove(data_points[remove_point_index])
            self.movable.set_xy_data(data_points)
        else:
            return
        


    # END Point movement #######################################################

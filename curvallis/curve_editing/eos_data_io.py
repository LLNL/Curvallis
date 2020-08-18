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

# Using OrderedDict to make output sections, and isotherms come out in the same
# order they were read in, to support comparison of outputs diring unit  testing.
from collections import OrderedDict
from curvallis.version import version as VERSION_STRING

"""This module provides input and output and a data structure for MEOS EOS
 output data.
"""

class _Reader_Writer_Base(object):
    """Provides operations and state used when reading or writing MEOS EOS
    output data files.  The data is stored in a .dat file. The .dat file
    contains each function's data.

    The .dat file contains a data section for each function in the .info file.
    The first line of each section is the function name, the number of isotherms,
    and the number of rho,func values for every isotherm.
    Each data line contains the current temperature, rho, and function value.
    The temperature value is the same in all the data lines for a given isotherm.
    """

    #These are strings in the files:
    _INFO_HEADER_LINE      = 'info'
    _INFO_SECTION_NAME     = 'material_info'
    _RHO_COUNT_FIELD_NAME  = 'numrho'
    _TEMP_COUNT_FIELD_NAME = 'numtemp'

    def __init__(self, base_name):
        self._base_name = base_name
        self._data_file = None
        self._data_file_name = '%s' % self._base_name

    ## Error handling support ##
    @staticmethod
    def _assert_eq(l, r):
        assert l == r, '%s%s /= %s%s' % (l, type(l), r, type(r))

    @staticmethod
    def _assert_in(l, r):
        assert l in r, "'%s' not in %s" % (l, r)

    @staticmethod
    def _assert_not_in(l, r):
        assert l not in r, "'%s' not not in %s" % (l, r)


class _EqIfFieldsEq(object):
    """ Objects of this class are equal if their components are equal.
    """
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __ne__(self, other):
        return not self.__eq__(other)


class _Isotherm(_EqIfFieldsEq):
    """ Contains the data for one isotherm: temperature and the rho-func
    values.
    """
    def __init__(self, temperature, num_points):
        super(_Isotherm, self).__init__()
        self.temp = temperature
        self.num_points = num_points
        # self.points later gets filled by the parser:
        self.points = list(range(self.num_points))

    def write(self, file_out):
        """ Write out the isotherm curve points for one temperature. No header,
        No blank lines between.

        :param file_out: INPUT parameter. The output file.
        """
        for point in self.points:
            # PyCharm mistakenly thinks point is an integer, so:
            # noinspection PyUnresolvedReferences
            file_out.writelines(' %.15E %.15E %.15E\n' %
                               (self.temp, point[0], point[1]))


class _Section(_EqIfFieldsEq):
    """ Contains the data for one section: the section name, fields, and isotherms.
    The number of isotherms is not known at create time.
    """
    def __init__(self, name):
        super(_Section, self).__init__()
        self.name = name
        self.num_isotherms = 0
        self.isotherm_num_points = 0
        self.fields = OrderedDict()
        self.isotherms = OrderedDict()

    def write(self, data_file):
        self._write_data(data_file)

    def _write_data(self, file_out):
        """ Write out one data section, with a trailing blank line.
        """
        file_out.writelines('%s    %s   %s\n' %
            (self.name, self.num_isotherms, self.isotherm_num_points))
        for isotherm in self.isotherms.values():
            isotherm.write(file_out)
        file_out.writelines('\n')

    def _values_to_string(self, values):
        """ Returns stringified value(s) separated by spaces.  Does not return
         quoted strings.
        """
        if isinstance(values, (list, tuple)):
            result = ''
            for value in values:
                result = '%s  %s' % (result, self._value_to_string(value))
        else:
            result = self._value_to_string(values)
        return result

    @staticmethod
    def _value_to_string(value):
        if isinstance(value, float):
            return '%1.15E' % value
        else:
            #Don't return any quotes for strings:
            return '%s' % value


class _Parser(_Reader_Writer_Base):
    """Provides file reading support for MEOS EOS output data.
    """
    def __init__(self, base_name, use_function):
        super(_Parser, self).__init__(base_name)
        self._sections = OrderedDict()
        self._read_called = False
        self._current_file = None
        self._current_line = None
        self._current_line_num = 0
        self._current_words = None
        self._current_field_name = None
        self._current_section = None
        self._current_section_name = None
        self._current_isotherm_temp = 0.0
        self._current_rho_num = 0
        self._use_function = use_function

    def read(self):
        """May only be called once.

        :return: dict of _Sections
        """
        self._assert_eq(self._read_called, False)
        self._read_called = True
        self._parse_data_file()
        return self._sections

    def _section_wanted(self, name):
        return self._use_function == name or self._use_function == 'all' or \
                self._use_function == self._INFO_SECTION_NAME

    def _current_section_wanted(self):
        return self._section_wanted(self._current_section_name)

    def _start_new_section(self):
        """ Reuse an existing section or create a new one.
        """
        name = self._current_field_name
        if self._section_wanted(name):
            if name not in self._sections:
                self._sections[name] = _Section(name)
            self._current_section = self._sections[name]
        self._current_section_name = name

    def _parse_data_file(self):
        message = 'Reading EOS data from "%s"' % self._data_file_name
        print (message)
        self._at_eof = False
        with open(self._data_file_name, 'r') as data_file:
            self._current_file = data_file
            self._current_line_num = 0
            self._parse_data_lines()
        print ('DONE ' + message)

    def _parse_data_lines(self):
        self._get_next_line_words()
        while not self._at_eof:
            self._parse_data_line()
            self._get_next_line_words()

    def _parse_data_line(self):
        if self._current_words[0].isalpha():
            self._start_section_data()
        else:
            if self._current_section_wanted():
                self._parse_data_values()

    def _start_section_data(self):
        """ Make an existing section the current section and create a new
        "isotherms" data component. Raise an exception if:
        - dimensions don't match
        - called more than once for a given section
        - the section does not exist
        - starting the same sections data more than once.
        - the section name is 'isotherms'
        """
        self._assert_word_count_eq(3)
        name = self._current_words[0]
        self._current_section_name = name
        self._current_isotherm_temp = 0.0
        self._current_rho_num = 0
        if self._section_wanted(name):
            data_num_temps = int(self._current_words[1])
            data_num_rhos = int(self._current_words[2])
            if name not in self._sections:
                self._sections[name] = _Section(name)
            self._current_section = self._sections[name]
            self._current_section.num_isotherms = data_num_temps
            self._current_section.isotherm_num_points = data_num_rhos

    def _parse_data_values(self):
        if len(self._current_words) == 3:        
            #Multiple Isotherms in this section
            temp = float(self._current_words[0])
            rho  = float(self._current_words[1])
            func = float(self._current_words[2])
            self._start_new_isotherm_if_needed(temp)
        else:
            #Single Isotherm in this section
            self._assert_word_count_eq(2)
            if (len(self._current_section.isotherms) != 0):
                assert self._current_isotherm_temp == 0.0,\
                    "Line doesn't contain 3 coordinates."
            rho  = float(self._current_words[0])
            func = float(self._current_words[1])
            self._start_new_isotherm_if_needed(float(0.0))
        self._current_section.isotherms[self._current_isotherm_temp].\
            points[self._current_rho_num] = [rho, func]
        self._current_rho_num += 1

    def _start_new_isotherm_if_needed(self, temp):
        """ Start a new isotherm if the temperature changed.  Be sure we have
        the right number of points, and that this is not a duplicate isotherm.
        """
        needed = False
        # Create the first one without checking the previous one:
        if len(self._current_section.isotherms) == 0:
            needed = True
        else:
            if temp != self._current_isotherm_temp:
                # Zero based (-1) rho_num post-increment (+1) is num_rhos:
                self._assert_eq(self._current_rho_num,
                                self._current_section.isotherm_num_points)
                self._assert_not_in(temp, self._current_section.isotherms)
                needed = True
        if needed:
            self._current_section.isotherms[temp] =_Isotherm(
                    temperature=temp,
                    num_points=self._current_section.isotherm_num_points)
            self._current_isotherm_temp = temp
            self._current_rho_num = 0

    def _get_next_line_words(self):
        """ Skip blank lines and get a list of words from the next non-blank line.
        Gets a max of 3 words. The last word includes all remaining chars on the line.

        len(line) == 0 at EOF.
        len(line) > 0 and len(line.split()) == 0 on a blank line.
        """
        self._current_line = None
        self._current_words = None
        done = False
        while not done:
            self._current_line = self._current_file.readline()
            self._current_line_num += 1
            if len(self._current_line) == 0:
                self._current_words = ()
                self._at_eof = True
                done = True
            else:
                # 2 means number of splits here, not number of words!
                self._current_words = self._current_line.split(None, 2)
                if len(self._current_words) > 0:
                    self._current_field_name = self._current_words[0]
                    done = True

    ## Error handling support ##

    def _assert_field_not_set(self):
        assert self._current_field_name not in self._current_section.fields,\
            '%s\nField %s.%s set again' %\
            (self._line_info(), self._current_section_name, self._current_field_name)

    def _assert_word_count_ge(self, count):
        assert len(self._current_words) >= count,\
            '%s\nExpected at least %s words, found %s' %\
            (self._line_info(), count, len(self._current_words))

    def _assert_word_count_eq(self, count):
        assert len(self._current_words) == count,\
            '%s\nExpected %s words, found %s' %\
            (self._line_info(), count, len(self._current_words))

    def _assert_field_name_eq(self, name):
        assert self._current_field_name == name,\
            '%s\nExpected field name "%s", found "%s"' %\
            (self._line_info(), name, self._current_words[0])

    def _line_info(self):
        return 'File %s\nLine %s: "%s"' %\
               (self._current_file.name, self._current_line_num, self._current_line)


class _Writer(_Reader_Writer_Base):
    """Provides file writing support for MEOS EOS output data.
    """
    def __init__(self, base_name):
        super(_Writer, self).__init__(base_name)

    def write(self, sections):
        message = 'Writing EOS data to "%s"' % self._data_file_name
        print (message)
        self._data_file = self._create_file(self._data_file_name)
        self._write_func_sections(sections)
        print ('DONE ' + message)

    @staticmethod
    def _create_file(name):
        try:
            result = open(name, 'w')
        except IOError as e:
            print('*** EXCEPTION %s %s while opening "%s" for writing' %
                  (type(e), e, name))
            result = None
        return result

    def _write_func_sections(self, sections):
        """ Write each function section, in the same order they were added
        """
        for name in sections.keys():
            if name != self._INFO_SECTION_NAME:
                sections[name].write(data_file=self._data_file)


class Data(_EqIfFieldsEq):
    """ Contains the actual data, and encapsulates I/O
    """

    def __init__(self, use_function):
        """

        :param use_function: str - which function to read in, or all.
        :return:
        """
        super(Data, self).__init__()
        self._use_function = use_function
        self.sections = {}

    def read(self, base_name):
        self.sections = _Parser(base_name,
                                self._use_function).read()

    def write(self, base_name):
        _Writer(base_name).write(self.sections)


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
Setup file for curvallis
"""

from setuptools import setup

setup(name='curvallis',
      version='1.0.0',
      packages=['curvallis', 'curvallis.curve_editing'],
      description='A MatPlotLib based curve editing tool',
      long_description = open("README.rst").read(),
      license='BSD',
      url='https://github.com/llnl/curvallis',
      maintainer='Paul Minner',
      maintainer_email='minner.paul@gmail.com',
      install_requires=['matplotlib', 'numpy', 'scipy', 'argparse'],
      entry_points = {
        'console_scripts': [
            'curvallis = curvallis.run:main'
            ]
        }
      )

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
Setup file for curvallis
@author: Eric Heinke (sudo-Eric)
"""

from curvallis.version import version as VERSION_STRING
from sys import version_info as version
from sys import argv
import subprocess
import platform
import argparse
import shutil
import sys

# https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/
# https://docs.python-guide.org/dev/virtualenvs/

# Begin Functions
##################################################
def vprint(message): # For verbose output
    if(args.verbose):
        print(message)
def display_help(): # Help message
    print(" [dir]\t\tThe directory you want to create the virtual python")
    print("\t\tenvironment in (leave blank for current directory)")
    exit()
def version_info(): # Version information
    return ("Virtual python environment tools version: " + VERSION_STRING)
def check_pip():
    try:
        subprocess.check_output([sys.executable,"-m" "pip","--version"])
        vprint("pip is installed.")
        return True
    except Exception as e:
        vprint(e)
        return False
def check_virtualenv():
    try:
        #subprocess.check_output(["virtualenv","--version"])
        subprocess.check_output([sys.executable,"-m" "virtualenv","--version"])
        vprint("virtualenv is installed.")
        return True
    except Exception as e:
        vprint(e)
        return False
def check_venv():
    try:
        subprocess.check_output([sys.executable,"-m","venv","-h"])
        vprint("venv is installed.")
        return True
    except Exception as e:
        vprint(e)
        return False
def create_virtualenv(directory):
    try:
        print("Creating virtual environment at " + directory)
        subprocess.check_output([sys.executable,"-m","virtualenv",directory])
        print("Virtual environment successfully created.")
    except Exception as e:
        vprint(e)
def create_venv(directory):
    try:
        print("Creating virtual environment at " + directory)
        subprocess.check_output([sys.executable,"-m","venv",directory])
        print("Virtual environment successfully created.")
    except Exception as e:
        vprint(e)
def install_virtualenv():
    try:
        print("Installing virtualenv...")
        subprocess.check_output([sys.executable,"-m","pip","install","--user","virtualenv"])
        print("Installed.")
    except Exception as e:
        vprint(e)
def install_venv():
    try:
        print("Installing venv...")
        subprocess.check_output([sys.executable,"-m","pip","install","--user","venv"])
        print("Installed.")
    except Exception as e:
        vprint(e)
##################################################
# Begin functions

# Start initial variables
##################################################
py_ver = [version.major,version.minor,version.micro,version.releaselevel]
##################################################
# End initial variables

# Begin Argument decoder
##################################################
# https://docs.python.org/3/library/argparse.html
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
# Visplay version information
parser.add_argument("--version", action="version",version = version_info())
# Enable verbose output
parser.add_argument("-v", "--verbose", help="enable verbose output", action="store_true")
# Create virtual python environment
group.add_argument("-c", help="create virtual python environment", metavar="DIR", type=str, action="store")
group.add_argument("-r", help="remove virtual python environment", metavar="DIR", type=str, action="store")
group.add_argument("-a", help="get instructions for activating virtual python environment", action="store_true")
group.add_argument("-d", help="get instructions for deactivating virtual python environment", action="store_true")
args = parser.parse_args()
##################################################
# End Argument decoder

# Begin virtual environment tools
##################################################
print("python " + str(py_ver[0]) + "." + str(py_ver[1]) + "." + str(py_ver[2]))
if(not(check_pip())):                               # Check if pip is installed
    print("pip is not installed.")
    print("Please install pip.")
    exit()
vprint("System: " + platform.system() + ", " + platform.release())
##################################################
if(not(args.c == None)):                            # Creating virtual python environment
    if((py_ver[0] == 2 and py_ver[1] == 7) \
    or (py_ver[0] == 3 and py_ver[1] < 3)):         # If python version is 2.7 or 3.0-2
        if(platform.system() == "Linux"):               # For Linux:
            if(not(check_virtualenv())):                    # Check if virtualenv is installed
                install_virtualenv()                        # Install virtualenv if needed
            create_virtualenv(args.c)                       # Create virtualenv
        elif(platform.system() == "Windows"):           # For Windows:
            if(not(check_virtualenv())):                    # Check if virtualenv is installed
                install_virtualenv()                        # Install virtualenv if needed
            create_virtualenv(args.c)                       # Create virtualenv
        elif(platform.system() == "Darwin"):            # For MacOS / OSX
            if(not(check_virtualenv())):                    # Check if virtualenv is installed
                install_virtualenv()                        # Install virtualenv if needed
            create_virtualenv(args.c)                       # Create virtualenv
        else:                                           # For unknown OS
            print("Unknown OS. No environment instructions available.")
            exit()
    elif(py_ver[0] == 3 and py_ver[1] >= 3):        # If python version is 3.3 or above
        if(platform.system() == "Linux"):               # For Linux:
            if(not(check_venv())):                          # Check if venv is installed
                install_venv()                              # Install venv if needed
            create_venv(args.c)                             # Create venv
        elif(platform.system() == "Windows"):           # For Windows:
            if(not(check_venv())):                          # Check if venv is installed
                install_venv()                              # Install venv if needed
            create_venv(args.c)                             # Create venv
        elif(platform.system() == "Darwin"):            # For MacOS / OSX
            if(not(check_venv())):                          # Check if venv is installed
                install_venv()                              # Install venv if needed
            create_venv(args.c)                             # Create venv
        else:                                           # For unknown OS
            print("Unknown OS. No environment instructions available.")
            exit()
    else:                                           # Unsupported version of python
        print("Unsupported version of python (" + str(py_ver[0]) + "." + str(py_ver[1]) + "." + str(py_ver[2]))
        exit()
##################################################
elif(not(args.r == None)):                          # Removing virtual python environment
    if(len(argv) < 3):  # Check if directory specified
        print("A directory must be specified.")
        display_help()
    print("Removing virtual python environment...")
    if(platform.system() == "Linux"):               # For Linux:
        try:
            shutil.rmtree(args.r)
            print("Virtual environment successfully removed.")
        except Exception as e:
            vprint(e)
    elif(platform.system() == "Windows"):           # For Windows:
        try:
            shutil.rmtree(args.r)
            print("Virtual environment successfully removed.")
        except Exception as e:
            vprint(e)
    elif(platform.system() == "Darwin"):            # For MacOS / OSX
        try:
            shutil.rmtree(args.r)
            print("Virtual environment successfully removed.")
        except Exception as e:
            vprint(e)
    else:                                           # For unknown OS
        print("Unknown OS. No environment instructions available.")
        exit()
##################################################
elif(not(args.a == None)):                           # Activating virtual python environment
    if(platform.system() == "Linux" \
    or platform.system() == "Darwin"):               # For Linux and MacOS:
        print("To activate the virtual environment, run the following command:")
        print("\"source venv/bin/activate\" where venv is the directory containing the virtual environment")
    elif(platform.system() == "Windows"):           # For Windows:
        print("To activate the virtual environemtn, run the following command:")
        print("\"venv\\Scripts\\activate.bat\" where venv is the directory containing the virtual environment")
    else:                                           # For unknown OS
        print("Unknown OS. No environment instructions available.")
        exit()
##################################################
elif(not(args.r == None)):                           # Deactivating virtual python environment
    if(platform.system() == "Linux" \
    or platform.system() == "Darwin"):               # For Linux and MacOS:
        print("To deactivate the virtual environment, run the following command:")
        print("\"deactivate\"")
    elif(platform.system() == "Windows"):           # For Windows:
        print("To deactivate the virtual environment, run the following command:")
        print("\"venv\\Scripts\\deactivate.bat\" where venv is the directory containing the virtual environment")
    else:                                           # For unknown OS
        print("Unknown OS. No environment instructions available.")
        exit()
##################################################
# End virtual environment tools


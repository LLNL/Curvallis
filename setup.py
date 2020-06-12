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
Written by: Eric Heinke
"""

from sys import version_info as version
from sys import argv
import subprocess
import platform
import sys
import pip
import os

# Start initial variables
##################################################
version_number = "1.1"
version_string = "6/12/2020"
version_message = "Tested with Linux (Ubuntu 18) and Windows 10 using python 2.7.17 and python 3.8."
package_installer_level = 1 # 1,2,3
verbose = False
modules_needed = []
generic = False
py_ver = [version.major,version.minor,version.micro,version.releaselevel]
py_path = ""
# Format for modules: ["module_name","module_version"] where "0" means latest version.
py27_modules = [["Tkinter","0"],["scipy","0"],["numpy","0"],["matplotlib","0"],["argparse","0"]]#,["pylab","0"]
py3_modules = [["tkinter","0"],["scipy","0"],["numpy","0"],["matplotlib","0"],["argparse","0"]]#,["pylab","0"]
##################################################
# End initial variables

# Start Functions
##################################################
def vprint(message): # For verbose output
    if(verbose):
        print(message)
def display_help(): # Help message
    print("usage: setup.py [option] ... [-env] ...")
    print("Options and arguments:")
    print(" -h,   --help\t\tDisplay this help message and exit")
    print("       --version\tOutput version information and exit")
    print(" -v,   --verbose\tEnable verbose output")
    print(" -py,  --python\t\tManually set python version")
    print("       --path\t\tManually set python executable's location (path)")
    print(" -p,   --PACKAGE\tSpecify aditional packages to install")
    print("       --GENERIC\tRun installer in generic mode")
    print(" -env, --ENVIRONMENT\tRun virtual python environment tools")
    print(" -l,   --level\t\tChanges the level of integration of the installer (only change is the default of 1 fails)")
    exit()
def version_info(): # Version information
    print("Curvallis readdyness tool version " + version_number + " (" + version_string + ")")
    print(version_message)
def decode_python_version(ver): # Decode python version
    tmp = ""
    new_ver = []
    for i in ver:
        if(i.isdigit()):
            tmp += i
        elif(i == '.'):
            new_ver.append(int(tmp))
            tmp = ""
        else:
            print("Python version format error. Python version must be entered in the " + \
            "format of x, x.x, or x.x.x.")
            exit()
    new_ver.append(int(tmp))
    if(len(new_ver) > 3):
        print("Python version format error. Python version must be entered in the " + \
        "format of x, x.x, or x.x.x.")
        exit()
    for i in range(len(new_ver)):
        py_ver[i] = new_ver[i]
    for i in range(len(new_ver), 3):
        py_ver[i] = 0
    py_ver[3] = "custom"
def module_checker(module_name): # Checks if module is installed
    try:
        exec ("import " + module_name)
        return True
    except:
        return False
def group_module_checker(modules): # Checks if group of module are installed
    needed = []
    for module in modules:
        if(not(module_checker(module[0]))):
            needed.append(module)
    return needed
def py_ver_as_string(): # Prints python version as a string
    return str(py_ver[0]) + "." + str(py_ver[1]) + "." + str(py_ver[2]) + " (" + str(py_ver[3]) + ")"
def install_package(package,version="0",method=3): # Package installer
    if(method == 1):
        if(hasattr(pip,"main")):
            pip.main(["install",package])
        else:
            pip._internal.main(['install', package])
    elif(method == 2):
        if(py_path == "python" or py_path == "python3"):
            subprocess.check_call([sys.executable,"-m","pip","install",package])
        else:
            subprocess.check_call([py_path,"-m","pip","install",package])
    elif(method == 3):
        if(version == "0"):
            os.system(py_path + " -m pip install " + package)
        else:
            os.system(py_path + " -m pip install " + package + "==" + version)
    # python -m pip install module==version
def add_package(name,version):
    if(version[0] == "-"):
        py27_modules.append([name,"0"])
        py3_modules.append([name,"0"])
        return 1
    else:
        py27_modules.append([name,version])
        py3_modules.append([name,version])
        return 2
##################################################
# End Functions

# Begin Argument decoder
##################################################
if(len(argv) == 1):
    pass
else:
    i = 1
    while(i < len(argv)):
        arg = argv[i].upper()
        if(arg == "--HELP" or arg == "-H"):         # Displays help message
            display_help()
        elif(arg == "--VERSION"):                   # Display version information
            version_info()
            exit()
        elif(arg == "--VERBOSE" or arg == "-V"):    # Enable verbose output
            verbose = True
        elif(arg == "--PYTHON" or arg == "-PY"):    # Set python version
            decode_python_version(argv[i+1])
            i += 1
        elif(arg == "--PATH"):                      # Set python executable location / path
            py_path = argv[i+1]
            i += 1
        elif(arg == "--GENERIC"):                   # Runs the installer in generic mode
                                                    # (attempt to install all modules regardless)
            generic = True
        elif(arg == "--PACKAGE" or arg == "-P"):    # Adds a module to the list of modules that need to be installed
            if(len(argv) >= i+3):
                i += add_package(argv[i+1],argv[i+2])
            else:
                i += add_package(argv[i+1],"0")
        elif(arg == "--ENVIRONMENT" or arg == "-ENV"):# Run virtual python environment tools
            vprint(platform.system() + ", " + platform.release())
            pass_args = [sys.executable,"venv_tools.py"]
            for a in range(i+1,len(argv)):
                pass_args.append(argv[a])
            subprocess.check_call(pass_args)
            exit()
        elif(arg == "--LEVEL" or arg == "-L"):      # Change install level
            tmp = int(argv[i+1])
            if(tmp < 1 or tmp > 3):
                print("Install level must be a number between 1 and 3.")
                exit()
            package_installer_level = tmp
            i += 1
        else:                                       # Bad argument given
            print("Invalid argument: " + argv[i])
            display_help()
        i += 1
##################################################
# End Argument decoder

# Start normal and verbose information about python, the os, and manually set variables
##################################################
version_info()
vprint("You are currently running Python " + py_ver_as_string() + " on " + \
        platform.system() + " " + platform.release() + ".")
if(not(py_path == "")):
    vprint("Manual override to python version " + py_ver_as_string() + ", located at " + py_path + ".")
elif(not(py_ver[0] == version.major and py_ver[1] == version.minor and py_ver[2] == version.micro)):
    vprint("Manual override to python version " + py_ver_as_string() + ".")
##################################################
# End normal and verbose information about python, the os, and manually set variables

# Start module checker
##################################################
if(py_ver[0] == 2 and py_ver[1] == 7):
    if(py_path == ""):
        py_path = "python"
    modules_needed = group_module_checker(py27_modules)
    if(generic):
        modules_needed = py27_modules
elif(py_ver[0] == 3):
    if(py_path == ""):
        py_path = "python3"
    modules_needed = group_module_checker(py3_modules)
    if(generic):
        modules_needed = py3_modules
else:
    print("Unsupported python version " + py_ver_as_string())
    exit()
if(len(modules_needed) == 0):
    print("All modules are already installed.")
    exit()
vprint("The following modules will be installed:")
for module in modules_needed:
    if(module[1] == "0"):
        vprint(" * " + module[0])
    else:
        vprint(" * " + module[0] + " version " + module[1])
##################################################
# End module checker

# Begin module installer
##################################################
print("Beginning instalation...")
for i in modules_needed:
    install_package(i[0],i[1],package_installer_level)
print("Done.")
##################################################
# End module installer

# Begin module verifier
##################################################
print("Verifying installed modules...")
modules_not_installed = []
if(py_ver[0] == 2 and py_ver[1] == 7):
    modules_not_installed = group_module_checker(py27_modules)
elif(py_ver[0] == 3):
    modules_not_installed = group_module_checker(py3_modules)
if(len(modules_not_installed) > 0):
    print("The following modules were not / could not be installed:")
    for module in modules_not_installed:
        if(module[1] == "0"):
            print(" * " + module[0])
        else:
            print(" * " + module[0] + " version " + module[1])
else:
    print("All modules were installed successfully.")
##################################################
# End module verifier


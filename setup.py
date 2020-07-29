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
@author: Eric Heinke (sudo-Eric)
"""

from curvallis.version import version as VERSION_STRING
from sys import version_info as version
from sys import argv
import subprocess
import platform
import argparse
import sys
import pip
import os

# Start Functions
##################################################
def vprint(message): # For verbose output
    if(args.verbose):
        print(message)
def version_info(): # Version information
    return ("Curvallis readiness tool version: " + VERSION_STRING)
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
def py_ver_as_string(): # Returns python version as a string
    return str(py_ver[0]) + "." + str(py_ver[1]) + "." + str(py_ver[2]) + " (" + str(py_ver[3]) + ")"
def install_package(package,version="0"): # Package installer
    try:
        if(args.level == 1):
            if(hasattr(pip,"main")):
                pip.main(["install",package])
            else:
                pip._internal.main(['install', package])
        elif(args.level == 2):
            if(version == "0"):
                subprocess.check_call([args.path,"-m","pip","install",package])
            else:
                subprocess.check_call([args.path,"-m","pip","install",package + "==" + version])
        elif(args.level == 3):
            if(version == "0"):
                os.system(args.path + " -m pip install " + package)
            else:
                os.system(args.path + " -m pip install " + package + "==" + version)
    except Exception:
        print("An error occured while installing " + package)
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
def decode_package(package):
    if("==" in package):
        return [package[0:package.find("==")], package[package.find("==")+2:len(package)]]
    else:
        return [package,"0"]
def check_pip():
    try:
        subprocess.check_output([sys.executable,"-m" "pip","--version"])
        vprint("pip is installed.")
        return True
    except Exception as e:
        vprint(e)
        return False
def write_hook(hook_code):
    hook = open(".git/hooks/pre-commit", 'w')
    for line in hook_code:
        hook.write(line + "\n")
    hook.close()
    subprocess.check_call(["chmod","+x","./.git/hooks/pre-commit"])
def update_hook(hook_code):
    pre_commit_file = open(".git/hooks/pre-commit",'r')
    old_hook = []
    for line in pre_commit_file:
        old_hook.append(line)
    pre_commit_file.close()
    new_hook = []
    i = 0
    while(i < len(old_hook)):
        if(old_hook[i] != "# versioning begin\n"):
            new_hook.append(old_hook[i])
        else:
            while(old_hook[i] != "# versioning end\n"):
                i += 1
        i += 1
    for line in hook_code:
        new_hook.append(line + "\n")
    pre_commit_file = open(".git/hooks/pre-commit", 'w')
    for line in new_hook:
        pre_commit_file.write(line)
    pre_commit_file.close()
##################################################
# End Functions

# Start initial variables
##################################################
modules_needed = []
py_ver = [0,0,0,"()"]
# Format for modules: ["module_name","module_version"] where "0" means latest version.
py27_modules = [["Tkinter","0"],["scipy","0"],["numpy","0"],["matplotlib","0"],["argparse","0"]]
py3_modules = [["tkinter","0"],["scipy","0"],["numpy","0"],["matplotlib","0"],["argparse","0"]]
##################################################
# End initial variables

# Begin Argument decoder
##################################################
# https://docs.python.org/3/library/argparse.html
parser = argparse.ArgumentParser()
# Visplay version information
parser.add_argument("--version", action="version",version = version_info())
# Enable verbose output
parser.add_argument("-v", "--verbose", help="enable verbose output", action="store_true")
# Set python version
parser.add_argument("-py", "--python", help="manually set python version", type=str)
# Set python executable location / path
parser.add_argument("--path", help="manually set pythons executable's location (path)", type=str, default=sys.executable)
# Adds modules to the list of modules that need to be installed
parser.add_argument("-p", "--package", help="specify aditional packages to install", type=str, default=[], nargs='+')
# Change installer level
parser.add_argument("-l", "--level", help="Changes the level of integration of the installer (only change if the default of 1 fails)", type=int, default=1)
# Runs the installer in generic mode (attempt to install all modules regardless)
parser.add_argument("--generic", help="run installer in generic mode", action="store_true")
# Run virtual python environment tools
parser.add_argument("-env", "--environment", help="run virtual python environment tools", nargs=argparse.REMAINDER)
# Allow for pre-commit hook for versioning to be added / updated
parser.add_argument('--versioning', help="add or update the pre-commit hook that allows for versioning", action="store_true")
args = parser.parse_args()
if(not(args.environment == None)):
    subprocess.check_call([sys.executable, "venv_tools.py"] + args.environment)
    exit()
if(args.level < 1 or args.level > 3):
    print("Installation level must be between 1 and 3.")
    exit()
if(args.python == None):
    py_ver = [version.major,version.minor,version.micro,version.releaselevel]
else:
    decode_python_version(args.python)
for i in args.package:
    modules_needed.append(decode_package(i))
##################################################
# End Argument decoder

# Start normal and verbose information about python, the os, and manually set variables
##################################################
version_info()
vprint("You are currently running Python " + py_ver_as_string() + " on " + \
        platform.system() + " " + platform.release() + ".")
vprint("Installer configured for python version " + py_ver_as_string() + ", located at " + args.path + ".")
##################################################
# End normal and verbose information about python, the os, and manually set variables

# Start versioning pre-commit hook installer / updater
##################################################
if(args.versioning):
    hook_code = [
        "# versioning begin",
        "dt=\"$(date -u)\"",
        "{ echo -n \"# This file is auto-generated. Any changes made to this file will be overwritten.\nversion = '\"; echo -n $dt; echo \"'\n\"; } > curvallis/version.py",
        "git add curvallis/version.py",
        "# versioning end"]
    if(os.path.exists(".git")): # Check if git has been initialized in current directory
        print("Git directory found.")
        if(os.path.exists(".git/hooks/pre-commit")): # Check if pre-commit hook already exists
            print("Pre-commit hook already exists. Updating relevant code.")
            update_hook(hook_code)
        else:   # Create pre-commit hook
            print("No existing pre-commit hook found. Creating new hook.")
            hook_code.insert(0,"#!/bin/sh\n")
            write_hook(hook_code)
    else:
        print("Git has not been initialized yet. Please initialize git with \"git init\".")
        exit()
##################################################
# End versioning pre-commit hook installer / updater

# Start module checker
##################################################
if(not(check_pip())): # Check if pip is installed
    print("pip is not installed.")
    print("Please install pip.")
    exit()
append_modules = []
if(py_ver[0] == 2 and py_ver[1] == 7):
    if(args.generic):
        append_modules = py27_modules
    else:
        append_modules = group_module_checker(py27_modules)
elif(py_ver[0] == 3):
    if(args.generic):
        append_modules = py3_modules
    else:
        append_modules = group_module_checker(py3_modules)
else:
    print("Unsupported python version " + py_ver_as_string())
    exit()
for i in append_modules:
    modules_needed.append(i)
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
print("Beginning installation...")
for i in modules_needed:
    install_package(i[0],i[1])
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
    print("You may need to install these modules manually or install them through your OS's package manager.")
    print("You could also try to upgrade pip to the latest version.") # python -m pip install --upgrade pip
else:
    print("All modules were installed successfully.")
##################################################
# End module verifier


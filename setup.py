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
import urllib.request
from sys import argv
import subprocess
import platform
import argparse
import datetime
import sys
import pip
import os


# Start Functions
##################################################
def vprint(message):  # For verbose output
    if args.verbose:
        print(message)


def version_info():  # Version information
    return "Curvallis readiness tool version: " + VERSION_STRING


def decode_python_version(ver):  # Decode python version
    tmp = ""
    new_ver = []
    for i in ver:
        if i.isdigit():
            tmp += i
        elif i == '.':
            new_ver.append(int(tmp))
            tmp = ""
        else:
            print("Python version format error. Python version must be entered in the " + \
                  "format of x, x.x, or x.x.x.")
            exit()
    new_ver.append(int(tmp))
    if len(new_ver) > 3:
        print("Python version format error. Python version must be entered in the " + \
              "format of x, x.x, or x.x.x.")
        exit()
    for i in range(len(new_ver)):
        py_ver[i] = new_ver[i]
    for i in range(len(new_ver), 3):
        py_ver[i] = 0
    py_ver[3] = "custom"


def module_checker(module_name):  # Checks if module is installed
    try:
        exec("import " + module_name)
        return True
    except:
        return False


def group_module_checker(modules):  # Checks if group of module are installed
    needed = []
    for module in modules:
        if not (module_checker(module[0])):
            needed.append(module)
    return needed


def py_ver_as_string():  # Returns python version as a string
    return str(py_ver[0]) + "." + str(py_ver[1]) + "." + str(py_ver[2]) + " (" + str(py_ver[3]) + ")"


def install_package(package, version="0"):  # Package installer
    try:
        if args.level == 1:
            if hasattr(pip, "main"):
                pip.main(["install", package])
            else:
                pip._internal.main(['install', package])
        elif args.level == 2:
            if version == "0":
                subprocess.check_call([args.path, "-m", "pip", "install", package])
            else:
                subprocess.check_call([args.path, "-m", "pip", "install", package + "==" + version])
        elif args.level == 3:
            if version == "0":
                os.system(args.path + " -m pip install " + package)
            else:
                os.system(args.path + " -m pip install " + package + "==" + version)
    except Exception:
        print("An error occurred while installing " + package)
    # python -m pip install module==version


def add_package(name, version):
    if version[0] == "-":
        py27_modules.append([name, "0"])
        py3_modules.append([name, "0"])
        return 1
    else:
        py27_modules.append([name, version])
        py3_modules.append([name, version])
        return 2


def decode_package(package):
    if "==" in package:
        return [package[0:package.find("==")], package[package.find("==") + 2:len(package)]]
    else:
        return [package, "0"]


def check_pip():
    try:
        subprocess.check_output([sys.executable, "-m" "pip", "--version"])
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
    if args.os == "posix":
        subprocess.check_call(["chmod", "+x", "./.git/hooks/pre-commit"])


def update_hook(hook_code):
    pre_commit_file = open(".git/hooks/pre-commit", 'r')
    old_hook = []
    for line in pre_commit_file:
        old_hook.append(line)
    pre_commit_file.close()
    new_hook = []
    i = 0
    while i < len(old_hook):
        if old_hook[i] != "# versioning begin\n":
            new_hook.append(old_hook[i])
        else:
            while old_hook[i] != "# versioning end\n":
                i += 1
        i += 1
    for line in hook_code:
        new_hook.append(line + "\n")
    pre_commit_file = open(".git/hooks/pre-commit", 'w')
    for line in new_hook:
        pre_commit_file.write(line)
    pre_commit_file.close()


# Check for updates at the given location. Return True if update available, False otherwise.
def check_for_updates(github_version_url):
    try:
        print("Checking for new version...")
        latest_version = urllib.request.urlopen(github_version_url)
        latest_version = str(latest_version.read())
        latest_version = latest_version[latest_version.find("version") + 11:latest_version.find("version") + 36]
    except Exception as e:
        print("An error occurred while collecting latest version info:")
        print(e)
        exit()
    # except urllib.error.HTTPError as e:
    # except urllib.error.URLError as e:
    # except urllib.error.ContentTooShortError as e:
    # from version import version as latest_version
    from curvallis.version import version as current_version

    # Convert date-time string to datetime object
    try:
        current_version = ISO8601_decode(current_version)
        latest_version = ISO8601_decode(latest_version)
    except ValueError:
        print("There was an error decoding version information.")
        exit()
    if current_version == latest_version:
        print("The latest version of Curvallis is already installed.")
        return False
    elif current_version < latest_version:
        print("A new version of Curvallis is available.")  # A new version of Curvallis was found
        return True
    elif current_version > latest_version:
        print("The current version of Curvallis is newer than the latest version.")
        if os.path.exists(".git"):
            versioning_setup = False
            pre_commit_file = open(".git/hooks/pre-commit", 'r')
            for line in pre_commit_file:
                if line == "# versioning begin\n":
                    versioning_setup = True
            pre_commit_file.close()
            if versioning_setup:
                print("Setup is unable to automatically update Curvallis.")
            else:
                print("Automatic versioning is nto set up. If you are developing Curvallis, please run this setup "
                      "file with \"--versioning\" flag to set up automatic versioning.")
        else:
            print("Git has not been initialized yet. Please initialize git with \"git init\".")
        return False
    else:
        print("An error occurred while comparing versions.")
        return False


# Converts an ISO8601 fate-time (YYY-MM-DDThh:mm:ss+hh:mm) to a datetime object
# https://en.wikipedia.org/wiki/ISO_8601
def ISO8601_decode(date_time):
    date_time = date_time[:22] + date_time[22 + 1:]
    return datetime.datetime.strptime(date_time, "%Y-%m-%dT%H:%M:%S%z")


##################################################
# End Functions

# Start initial variables
##################################################
modules_needed = []
py_ver = [0, 0, 0, "()"]
# Format for modules: ["module_name","module_version"] where "0" means latest version.
py27_modules = [["Tkinter", "0"], ["scipy", "0"], ["numpy", "0"], ["matplotlib", "0"], ["argparse", "0"]]
py3_modules = [["tkinter", "0"], ["scipy", "0"], ["numpy", "0"], ["matplotlib", "0"], ["argparse", "0"]]
curvallis_github_url = "https://github.com/%s/Curvallis.git"
curvallis_github_version_url = "https://raw.githubusercontent.com/%s/Curvallis/master/curvallis/version.py"
curvallis_github_download_page_url = "https://github.com/%s/Curvallis/archive/refs/heads/master.zip"
##################################################
# End initial variables

# Begin Argument decoder
##################################################
# https://docs.python.org/3/library/argparse.html
parser = argparse.ArgumentParser()
# Visplay version information
parser.add_argument("--version", action="version", version=version_info())
# Display debug info (for helping debug user problems)
parser.add_argument('--get-info', help="Display some useful information about current setup",
                    action="store", nargs="?", default="NONE", dest="get_info")
# Enable verbose output
parser.add_argument("-v", "--verbose", help="enable verbose output", action="store_true")
# Set python version
parser.add_argument("-py", "--python", help="manually set python version", type=str)
# Set python executable location / path
parser.add_argument("--path", help="manually set pythons executable's location (path)", type=str,
                    default=sys.executable)
# Set operating system type as either posix or nt
parser.add_argument("--os", help="manually set operating system (\"nt\" for Windows and \"posix\" for Linux/MacOS)",
                    type=str, default=os.name)
# Adds modules to the list of modules that need to be installed
parser.add_argument("-p", "--package", help="specify aditional packages to install", type=str, default=[], nargs='+')
# Change installer level
parser.add_argument("-l", "--level",
                    help="Changes the level of integration of the installer (only change if the default of 1 fails)",
                    type=int, default=1)
# Runs the installer in generic mode (attempt to install all modules regardless)
parser.add_argument("--generic", help="run installer in generic mode", action="store_true")
# Run virtual python environment tools
parser.add_argument("-env", "--environment", help="run virtual python environment tools", nargs=argparse.REMAINDER)
# Allow for pre-commit hook for versioning to be added / updated
parser.add_argument("--versioning", help="add or update the pre-commit hook that allows for versioning",
                    action="store_true")
# Check to sss if a new version of Curvallis is available
parser.add_argument("--check-update", help="check for new version of Curvallis", action="store_true")
# Update Curvallis to the latest version if an update is available
parser.add_argument("--update", help="update to new version of Curvallis if available", action="store_true")
# Set GitHub repo for version checking and updating
parser.add_argument("--repo", help="manually set which repo is used for updates (enter only username, eg. LLNL)",
                    type=str, default="LLNL")
args = parser.parse_args()

# Display useful info about configuration
if args.get_info is None or args.get_info.upper() != "NONE":
    if args.get_info is None:
        args.get_info = "STANDARD"
    if args.get_info.upper() == "HELP":
        print("%s --get-info [NONE, STANDARD, FULL, HELP]" % parser.prog)
        exit()
    print("========== Curvallis Debug Info ==========")
    print("Python version: %d.%d.%d %s" % (version.major, version.minor, version.micro, version.releaselevel))
    print("Curvallis version: %s" % VERSION_STRING)
    try: from matplotlib import __version__ as matplotlib_version
    except: matplotlib_version = "Not installed"
    print("Matplotlib version: %s" % matplotlib_version)
    try: from tkinter import TkVersion as tkinter_version
    except: tkinter_version = "Not installed"
    print("Tkinter version version: %s" % tkinter_version)
    try: from scipy import __version__ as scipy_version
    except: scipy_version = "Not installed"
    print("Scipy version: %s" % scipy_version)
    try: from numpy import __version__ as numpy_version
    except: numpy_version = "Not installed"
    print("Numpy version: %s" % numpy_version)
    if args.get_info.upper() == "FULL":
        print("System: %s %s" % (platform.system(), platform.release()))
        print("System type: %s" % args.os)
        print("Git initialized: %s" % os.path.exists(".git"))
        if os.path.exists(".git" + os.path.sep + "config"):
            tmp_file = open(".git" + os.path.sep + "config", 'r')
            lines = tmp_file.readlines()
            tmp_file.close()
            for i in range(len(lines)):
                if lines[i].strip() == "[remote \"origin\"]":
                    line = lines[i+1].strip()
                    line = line[line.find('=')+1:].strip()
                    print("Git repository URL: %s" % line)
                    break
    exit()

if not (args.environment is None):
    try:
        subprocess.check_call([sys.executable, "venv_tools.py"] + args.environment)
    finally:
        exit()
if args.level < 1 or args.level > 3:
    print("Installation level must be between 1 and 3.")
    exit()
if args.python is None:
    py_ver = [version.major, version.minor, version.micro, version.releaselevel]
else:
    decode_python_version(args.python)
for i in args.package:
    modules_needed.append(decode_package(i))
if args.os != 'posix' and args.os != 'nt':
    print("Unknown or unsupported operating system detected.")
    exit()
curvallis_github_url = curvallis_github_url % args.repo
curvallis_github_version_url = curvallis_github_version_url % args.repo
curvallis_github_download_page_url = curvallis_github_download_page_url % args.repo
##################################################
# End Argument decoder

# Start normal and verbose information about python, the os, and manually set variables
##################################################
version_info()
vprint("You are currently running Python " + py_ver_as_string() + " on "
       + platform.system() + " " + platform.release() + ".")
vprint("Installer configured for python version " + py_ver_as_string() + ", located at " + args.path + ".")
##################################################
# End normal and verbose information about python, the os, and manually set variables

# Start update checker
##################################################
if args.check_update:
    print(curvallis_github_version_url)
    update_available = check_for_updates(curvallis_github_version_url)
    if update_available:
        prompt = input("Would you like to update Curvallis? ").upper()
        if not (prompt == "Y" or prompt == "YES"):
            exit()
        # If 'Y' or 'YES', continue to the update code below
        args.update = True
    else:
        exit()
##################################################
# End update checker

# Start update installer
##################################################
if args.update:
    # If update has not already been checked for, check for it.
    if not args.check_update and False:
        update_available = check_for_updates(curvallis_github_version_url)
        if not update_available:
            print("There are no updates available for Curvallis.")
            exit()
    # Check if git has been set up for Curvallis
    if os.path.exists(".git"):
        print("Git directory found.")
        if args.os == 'posix' or args.os == 'nt':
            subprocess.check_call(["git", "pull", curvallis_github_url, "master"])
    else:
        # If git has not been set up, ask if the user wants to do a manual update
        print("Git directory not found.")
        if args.os == 'nt':
            print("Alternate update method is not supported on windows.")
            exit()
        prompt = input("Would you like to switch to the alternate update method?\n"
                       "Warning: This will overwrite your current installation (y/n) ").upper()
        if not (prompt == "Y" or prompt == "YES"):
            exit()
        if args.os == 'posix':
            try:
                if os.path.exists("master.zip"):
                    os.remove("master.zip")
                # Download Curvallis archive
                print("Downloading Curvallis...")
                subprocess.check_call(["wget", "-q", curvallis_github_download_page_url])
                # Unzip archive
                print("Unzipping files...")
                subprocess.check_call(["unzip", "-o", "-q", "master.zip"])
                # Copy uncompressed archive files to correct place, overwriting existing files
                print("Processing files...")
                for file in os.listdir("./Curvallis-master"):
                    subprocess.check_call(["cp", "-rf", "./Curvallis-master/" + file, "."])
                # Remove temporary/unneeded files
                print("Cleaning up...")
                subprocess.check_call(["rm", "master.zip"])
                subprocess.check_call(["rm", "-rf", "./Curvallis-master"])
                print("Done.\nUpdate successful!")
            except Exception as e:
                print("An error has occurred while trying to update.")
                vprint(e)
        elif args.os == 'nt':
            pass  # Not supported yet
    exit()
    # check if git is initialized or not
    # subprocess.check_call([args.path, "git pull https://github.com/LLNL/Curvallis.git master"])
##################################################
# End update installer

# Start versioning pre-commit hook installer / updater
##################################################
if args.versioning:
    if args.os == "posix" or args.os == "nt":
        hook_code = [
            "# versioning begin",
            "dt=\"$(date -u --iso-8601=seconds)\"",
            "{ echo -n \"# This file is auto-generated. Any changes made to this file will be overwritten.\nversion = "
            "'\"; echo -n $dt; echo \"'\n\"; } > curvallis/version.py",
            "git add curvallis/version.py",
            "# versioning end"]
        if os.path.exists(".git"):  # Check if git has been initialized in current directory
            print("Git directory found.")
            if os.path.exists(".git/hooks/pre-commit"):  # Check if pre-commit hook already exists
                print("Pre-commit hook already exists. Updating relevant code.")
                update_hook(hook_code)
            else:  # Create pre-commit hook
                print("No existing pre-commit hook found. Creating new hook.")
                hook_code.insert(0, "#!/bin/sh\n")
                write_hook(hook_code)
        else:
            print("Git has not been initialized yet. Please initialize git.")
    exit()
##################################################
# End versioning pre-commit hook installer / updater

# Start module checker
##################################################
if not check_pip():  # Check if pip is installed
    print("pip is not installed.")
    print("Please install pip.")
    exit()
append_modules = []
# If python version 2.7.X is detected
if py_ver[0] == 2 and py_ver[1] == 7:
    if args.generic:
        append_modules = py27_modules
    else:
        append_modules = group_module_checker(py27_modules)
    # Python 2.7.X is not longer supported. Display message and exit.
    print("Unsupported python version " + py_ver_as_string())
    print("Python 2.7 is no longer supported by Curvallis. Please try again with Python 3.")
    exit()
# If python version 3.X.X is detected
elif py_ver[0] == 3:
    if args.generic:
        append_modules = py3_modules
    else:
        append_modules = group_module_checker(py3_modules)
else:
    print("Unsupported python version " + py_ver_as_string())
    exit()
for i in append_modules:
    modules_needed.append(i)
if len(modules_needed) == 0:
    print("All modules are already installed.")
    exit()
vprint("The following modules will be installed:")
for module in modules_needed:
    if module[1] == "0":
        vprint(" * " + module[0])
    else:
        vprint(" * " + module[0] + " version " + module[1])
##################################################
# End module checker

# Begin module installer
##################################################
print("Beginning installation...")
for i in modules_needed:
    install_package(i[0], i[1])
print("Done.")
##################################################
# End module installer

# Begin module verifier
##################################################
print("Verifying installed modules...")
modules_not_installed = []
if py_ver[0] == 2 and py_ver[1] == 7:
    modules_not_installed = group_module_checker(py27_modules)
elif py_ver[0] == 3:
    modules_not_installed = group_module_checker(py3_modules)
if len(modules_not_installed) > 0:
    print("The following modules were not / could not be installed:")
    for module in modules_not_installed:
        if module[1] == "0":
            print(" * " + module[0])
        else:
            print(" * " + module[0] + " version " + module[1])
    print("You may need to install these modules manually or install them through your OS's package manager.")
    print("You could also try to upgrade pip to the latest version.")  # python -m pip install --upgrade pip
else:
    print("All modules were installed successfully.")
##################################################
# End module verifier

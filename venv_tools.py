from sys import version_info as version
from sys import argv
import subprocess
import platform
import shutil
import sys
import os

# https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/
# https://docs.python-guide.org/dev/virtualenvs/

# Start initial variables
##################################################
version_number = "1.1"
version_string = "6/12/2020"
version_message = "Tested with Linux (Ubuntu 18) and Windows 10 using python 2.7.17 and python 3.8."
verbose = False
py_ver = [version.major,version.minor,version.micro,version.releaselevel]
venv_dir = ""
##################################################
# End initial variables

# Begin General Functions
##################################################
def vprint(message): # For verbose output
    if(verbose):
        print(message)
def display_help(): # Help message
    print("usage: venv_tools.py [option] ... [dir]")
    print("Options and arguments:")
    print(" -h,   --help\t\tDisplay this help message and exit")
    print("       --version\tOutput version information and exit")
    print(" -v,   --verbose\tEnable verbose output")
    print(" -c,   --create\t\tCreate virtual python environment")
    print(" -r,   --remove\t\tRemove virtual python environment")
    print(" -a,   --activate\tGet instructions for activating virtual python environment")
    print(" -d,   --deactivate\tGet instructions for deactivating virtual python environment")
    print("")
    print(" [dir]\t\tThe directory you want to create the virtual python")
    print("\t\tenvironment in (leave blank for current directory)")
    exit()
def version_info(): # Version information
    print("Curvallis readdyness tool version " + version_number + " (" + version_string + ")")
    print(version_message)
##################################################
# End General Functions

# Begin virtual python functions
##################################################
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
# Begin virtual python functions

# Begin Argument decoder
##################################################
function = ''
if(len(argv) == 1):
    display_help()
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
        elif(arg == "--CREATE" or arg == "-C"):     # Create virtual python environment
            function = 'c'
        elif(arg == "--REMOVE" or arg == "-R"):     # Remove virtual python environment
            function = 'r'
        elif(arg == "--ACTIVATE" or arg == "-A"):   # Activate virtual python environment
            function = 'a'
        elif(arg == "--DEACTIVATE" or arg == "-D"): # Deactivate virtual python environment
            function = 'd'
        elif(i == len(argv)-1 and not(arg[0] == '-')):# directory given
            venv_dir = argv[i]
        else:                                       # Bad argument given
            print("Invalid argument: " + argv[i])
            display_help()
        i += 1
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
if(function == 'c'):                            # Creating virtual python environment
    if((py_ver[0] == 2 and py_ver[1] == 7) \
    or (py_ver[0] == 3 and py_ver[1] < 3)):         # If python version is 2.7 or 3.0-2
        if(platform.system() == "Linux"):               # For Linux:
            if(not(check_virtualenv())):                    # Check if virtualenv is installed
                install_virtualenv()                        # Install virtualenv if needed
            create_virtualenv(venv_dir)                     # Create virtualenv
        elif(platform.system() == "Windows"):           # For Windows:
            if(not(check_virtualenv())):                    # Check if virtualenv is installed
                install_virtualenv()                        # Install virtualenv if needed
            create_virtualenv(venv_dir)                     # Create virtualenv
        elif(platform.system() == "Darwin"):            # For MacOS / OSX
            if(not(check_virtualenv())):                    # Check if virtualenv is installed
                install_virtualenv()                        # Install virtualenv if needed
            create_virtualenv(venv_dir)                     # Create virtualenv
        else:                                           # For unknown OS
            print("Unknown OS. No environment instructions available.")
            exit()
    elif(py_ver[0] == 3 and py_ver[1] >= 3):        # If python version is 3.3 or above
        if(platform.system() == "Linux"):               # For Linux:
            if(not(check_venv())):                          # Check if venv is installed
                install_venv()                              # Install venv if needed
            create_venv(venv_dir)                           # Create venv
        elif(platform.system() == "Windows"):           # For Windows:
            if(not(check_venv())):                          # Check if venv is installed
                install_venv()                              # Install venv if needed
            create_venv(venv_dir)                           # Create venv
        elif(platform.system() == "Darwin"):            # For MacOS / OSX
            if(not(check_venv())):                          # Check if venv is installed
                install_venv()                              # Install venv if needed
            create_venv(venv_dir)                           # Create venv
        else:                                           # For unknown OS
            print("Unknown OS. No environment instructions available.")
            exit()
    else:                                           # Unsupported version of python
        print("Unsupported version of python (" + str(py_ver[0]) + "." + str(py_ver[1]) + "." + str(py_ver[2]))
        exit()
##################################################
elif(function == 'r'):                          # Removing virtual python environment
    if(len(argv) < 3):  # Check if directory specified
        print("A directory must be specified.")
        display_help()
    print("Removing virtual python environment...")
    if(platform.system() == "Linux"):               # For Linux:
        try:
            shutil.rmtree(venv_dir)
            #subprocess.check_output(["rm","-rf",venv_dir])# Delete virtual python's directory
            print("Virtual environment successfully removed.")
        except Exception as e:
            vprint(e)
    elif(platform.system() == "Windows"):           # For Windows:
        try:
            shutil.rmtree(venv_dir)                     # Delete virtual python's directory
            print("Virtual environment successfully removed.")
        except Exception as e:
            vprint(e)
    elif(platform.system() == "Darwin"):            # For MacOS / OSX
        try:
            shutil.rmtree(venv_dir)
            #subprocess.check_output(["rm","-rf",venv_dir])# Delete virtual python's directory
            print("Virtual environment successfully removed.")
        except Exception as e:
            vprint(e)
    else:                                           # For unknown OS
        print("Unknown OS. No environment instructions available.")
        exit()
##################################################
elif(function == 'a'):                          # Activating virtual python environment
    if(platform.system() == "Linux" \
    or platform.system() == "Darwin"):               # For Linux and MacOS:
        print("To activate the virtual environment, run the following command:")
        print("\"source venv/bin/activate\" where venv is the directory containing the virtual environemtn")
    elif(platform.system() == "Windows"):           # For Windows:
        print("To activate the virtual environemtn, run the following command:")
        print("\"venv\\Scripts\\activate.bat\" where venv is the directory containing the virtual environemtn")
    else:                                           # For unknown OS
        print("Unknown OS. No environment instructions available.")
        exit()
##################################################
elif(function == 'd'):                          # Deactivating virtual python environment
    if(platform.system() == "Linux" \
    or platform.system() == "Darwin"):               # For Linux and MacOS:
        print("To deactivate the virtual environment, run the following command:")
        print("\"deactivate\"")
    elif(platform.system() == "Windows"):           # For Windows:
        print("To deactivate the virtual environment, run the following command:")
        print("\"venv\\Scripts\\deactivate.bat\" where venv is the directory containing the virtual environemtn")
    else:                                           # For unknown OS
        print("Unknown OS. No environment instructions available.")
        exit()
##################################################
# End virtual environment tools


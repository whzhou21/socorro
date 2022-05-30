#!/usr/bin/env python

"""
Install.py is a helper tool to download, unpack, and build this library
This is used to automate steps described in the README file in this dir
"""

# *** IWYU Modules *************************************************** #

from __future__ import print_function
from argparse import ArgumentParser
import os, shutil, subprocess, sys

# *** Library version ************************************************ #

version = "master"

# *** Help message *************************************************** #

HELP_MESSAGE = """
Examples:

Syntax from lib dir: python Install.py -d
                 or: python Install.py -b
                 or: python Install.py -c

Syntax from src dir: make lib-fftmpi args="-d"
                 or: make lib-fftmpi args="-b"
                 or: make lib-fftmpi args="-c"
"""

# *** Input arguments ************************************************ #

parser = ArgumentParser(prog = 'Install.py', description = "Helper script to download and build the FFTMPI library")

parser.add_argument("-d", action = "store_true", help = "download the FFTMPI library to ../lib/fftmpi")
parser.add_argument("-b", action = "store_true", help = "build the FFTMPI library")
parser.add_argument("-c", action = "store_true", help = "clean the FFTMPI build space")

args = parser.parse_args()

cloneflag = args.d
buildflag = args.b
cleanflag = args.c

clonepath = os.path.join(os.path.abspath(os.path.expanduser(".")),"fftmpi")
buildpath = os.path.join(clonepath,"src")
libprefix = os.path.join(buildpath,"bin")

# *** Print the help message ***************************************** #

if not cloneflag and not buildflag and not cleanflag:
   parser.print_help()
   sys.exit(HELP_MESSAGE)

# *** Download the library ******************************************* #

cmd = 'git clone --depth 1 --branch %s https://github.com/lammps/fftmpi.git %s' % (version,clonepath)

if cloneflag:
   print("Downloading the FFTMPI library ...")
   if os.path.isdir(clonepath):
      print("ERROR: The destination path '" + clonepath + "' already exists.")
      sys.exit(1)
   try:
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
   except subprocess.CalledProcessError as e:
      print("Download failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Build the library ********************************************** #

cmd = 'cd %s && make lib CC=mpicxx CCFLAGS="-O3" fft=FFTW3 p=DOUBLE' % (buildpath)

if buildflag:
   print("Building the FFTMPI library ...")
   if not os.path.isdir(buildpath):
      print("ERROR: The destination path '" + buildpath + "' does not exist.")
      sys.exit(1)
   try:
      txt = subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
      print(txt)
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Clean the build space ****************************************** #

cmd = 'cd %s && make clean-all' % (buildpath)

if cleanflag:
   print("Cleaning the FFTMPI build space ...")
   if not os.path.isdir(buildpath):
      print("ERROR: The destination path '" + buildpath + "' does not exist.")
      sys.exit(1)
   try:
      txt = subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
      print(txt)
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** End of the file ************************************************ #

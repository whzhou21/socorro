#!/usr/bin/env python

"""
Install.py is a helper tool to download, unpack, and build this library
This is used to automate steps described in the README file in this dir
"""

# *** IWYU Modules *************************************************** #

from __future__ import print_function
from argparse import ArgumentParser
import os, shutil, subprocess, sys

# *** Library settings *********************************************** #

libname = "FFTMPI"
version = "master"

# *** Help message *************************************************** #

HELP_MESSAGE = """
Examples:

Syntax from lib dir: python Install.py -d
                 or: python Install.py -c
                 or: python Install.py -b
"""

# *** Input arguments ************************************************ #

parser = ArgumentParser(prog = 'Install.py', description = "Helper script to download and build the " + libname + " library")

parser.add_argument("-b", action = "store_true", help = "build the " + libname + " library")
parser.add_argument("-c", action = "store_true", help = "clean the " + libname + " build space")
parser.add_argument("-d", action = "store_true", help = "download the " + libname + " library to lib/" + libname.lower())

args = parser.parse_args()

buildflag = args.b
cleanflag = args.c
cloneflag = args.d
buildpath = os.path.join(os.path.abspath(os.path.expanduser(".")), libname.lower())

# *** Print the help message ***************************************** #

if not buildflag and not cleanflag and not cloneflag:
   parser.print_help()
   sys.exit(HELP_MESSAGE)

# *** Download the library ******************************************* #

url = 'https://github.com/lammps/fftmpi.git'
cmd = 'git clone --depth 1 --branch %s %s %s' % (version, url, buildpath)

if cloneflag:
   print("Downloading the " + libname + " library ...")
   try:
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
   except subprocess.CalledProcessError as e:
      print("Download failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Clean the build space ****************************************** #

cmd = 'cd %s/src/ && make clean-all' % (buildpath)

if cleanflag:
   print("Cleaning the " + libname + " build space ...")
   try:
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Build the library ********************************************** #

bin = os.path.join(buildpath, "bin")
cmd = 'cd %s/src/ && make lib CC=mpicxx CCFLAGS="-O3" fft=FFTW3 p=DOUBLE' % (buildpath)

if buildflag:
   print("Building the " + libname + " library ...")
   try:
      txt = subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
      print(txt)
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** End of the file ************************************************ #

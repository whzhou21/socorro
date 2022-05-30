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

version = "v2.2.1"

# *** Help message *************************************************** #

HELP_MESSAGE = """
Examples:

Syntax from lib dir: python Install.py -d
                 or: python Install.py -b
                 or: python Install.py -c

Syntax from src dir: make lib-scalapack args="-d"
                 or: make lib-scalapack args="-b"
                 or: make lib-scalapack args="-c"
"""

# *** Input arguments ************************************************ #

parser = ArgumentParser(prog = 'Install.py', description = "Helper script to download and build the SCALAPACK library")

parser.add_argument("-d", action = "store_true", help = "download the SCALAPACK library to ../lib/scalapack") 
parser.add_argument("-b", action = "store_true", help = "build the SCALAPACK library")
parser.add_argument("-c", action = "store_true", help = "clean the SCALAPACK build space")

args = parser.parse_args()

cloneflag = args.d
buildflag = args.b
cleanflag = args.c

clonepath = os.path.join(os.path.abspath(os.path.expanduser(".")),"scalapack")
buildpath = clonepath
libprefix = os.path.join(buildpath,"bin")

# *** Print the help message ***************************************** #

if not cloneflag and not buildflag and not cleanflag:
   parser.print_help()
   sys.exit(HELP_MESSAGE)

# *** Download the library ******************************************* #

cmd = 'git clone --depth 1 --branch %s https://github.com/Reference-ScaLAPACK/scalapack.git %s' % (version,clonepath)

if cloneflag:
   print("Downloading the SCALAPACK library ...")
   if os.path.isdir(clonepath):
      print("ERROR: The destination path '" + clonepath + "' already exists.")
      sys.exit(1)
   try:
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
   except subprocess.CalledProcessError as e:
      print("Download failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Build the library ********************************************** #

cmd = 'cd %s && cp SLmake.inc.example SLmake.inc && make lib FC=mpifort FCFLAGS="-O3 -fallow-argument-mismatch" CC=mpicc CCFLAGS="-O3 -Wno-error=implicit-function-declaration"' % (buildpath)

if buildflag:
   print("Building the SCALAPACK library ...")
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

cmd = 'cd %s && cp SLmake.inc.example SLmake.inc && make clean' % (buildpath)

if cleanflag:
   print("Cleaning the SCALAPACK build space ...")
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

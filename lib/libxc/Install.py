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

libname = "LIBXC"
version = "2.1.2"

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

url = 'https://gitlab.com/libxc/libxc.git'
cmd = 'git clone --depth 1 --branch %s %s %s' % (version, url, buildpath)

if cloneflag:
   print("Downloading the " + libname + " library ...")
   try:
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
   except subprocess.CalledProcessError as e:
      print("Download failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Clean the build space ****************************************** #

cmd = 'cd %s && make clean' % (buildpath)

if cleanflag:
   print("Cleaning the " + libname + " build space ...")
   try:
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** Build the library ********************************************** #

bin = os.path.join(buildpath, "bin")
cmd = 'cd %s && autoreconf -i && ./configure FC=mpifort FCFLAGS="-O3" CC=mpicc CFLAGS="-O3" --prefix=%s && make && make install' % (buildpath, bin)

if buildflag:
   print("Building the " + libname + " library ...")
   try:
      txt = subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
      print(txt)
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** End of the file ************************************************ #

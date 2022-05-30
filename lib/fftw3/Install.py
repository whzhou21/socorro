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

libname = "FFTW3"
version = "3.3.10"

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

src = 'fftw-' + version
out = src + '.tar.gz'
url = 'https://fftw.org/pub/fftw/' + out

if cloneflag:
   print("Downloading the " + libname + " library ...")
   success = False
   if not success:
      cmd = 'curl -L -o "%s" %s' % (out,url)
      try:
         subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
         success = True
      except subprocess.CalledProcessError as e:
         print("Download using curl failed with: %s" % e.output.decode('UTF-8'))
   if not success:
      cmd = 'wget -O "%s" %s' % (out,url)
      try:
         subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
         success = True
      except subprocess.CalledProcessError as e:
         print("Download using wget failed with: %s" % e.output.decode('UTF-8'))
   if not success:
      print('Failed to download source code with "curl" or "wget"')
      sys.exit(1)
   else:
      cmd = 'tar -zxf %s && mv %s %s && rm %s' % (out,src,buildpath,out)
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')

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
cmd = 'cd %s && ./configure CC=mpicc CFLAGS="-O3" --enable-threads --prefix=%s && make && make install' % (buildpath, bin)

if buildflag:
   print("Building the " + libname + " library ...")
   try:
      txt = subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
      print(txt)
   except subprocess.CalledProcessError as e:
      print("Make failed with:\n %s" % e.output.decode('UTF-8'))
      sys.exit(1)

# *** End of the file ************************************************ #

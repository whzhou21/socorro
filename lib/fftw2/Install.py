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

version = "2.1.5"

# *** Help message *************************************************** #

HELP_MESSAGE = """
Examples:

Syntax from lib dir: python Install.py -d
                 or: python Install.py -b
                 or: python Install.py -c

Syntax from src dir: make lib-fftw2 args="-d"
                 or: make lib-fftw2 args="-b"
                 or: make lib-fftw2 args="-c"
"""

# *** Input arguments ************************************************ #

parser = ArgumentParser(prog = 'Install.py', description = "Helper script to download and build the FFTW2 library")

parser.add_argument("-d", action = "store_true", help = "download the FFTW2 library to ../lib/fftw2")
parser.add_argument("-b", action = "store_true", help = "build the FFTW2 library")
parser.add_argument("-c", action = "store_true", help = "clean the FFTW2 build space")

args = parser.parse_args()

cloneflag = args.d
buildflag = args.b
cleanflag = args.c

clonepath = os.path.join(os.path.abspath(os.path.expanduser(".")),"fftw2")
buildpath = clonepath
libprefix = os.path.join(buildpath,"bin")

# *** Print the help message ***************************************** #

if not buildflag and not cleanflag and not cloneflag:
   parser.print_help()
   sys.exit(HELP_MESSAGE)

# *** Download the library ******************************************* #

src = 'fftw-' + version
out = src + '.tar.gz'
url = 'https://fftw.org/pub/fftw/' + out

if cloneflag:
   print("Downloading the FFTW2 library ...")
   if os.path.isdir(clonepath):
      print("ERROR: The destination path '" + clonepath + "' already exists.")
      sys.exit(1)
   success = False
   if not success:
      cmd = 'curl -L -o %s %s' % (out,url)
      try:
         subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
         success = True
      except subprocess.CalledProcessError as e:
         print("Download using curl failed with: %s" % e.output.decode('UTF-8'))
   if not success:
      cmd = 'wget -O %s %s' % (out,url)
      try:
         subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
         success = True
      except subprocess.CalledProcessError as e:
         print("Download using wget failed with: %s" % e.output.decode('UTF-8'))
   if not success:
      print('Failed to download source code with "curl" or "wget"')
      sys.exit(1)
   else:
      cmd = 'tar -zxf %s && mv %s %s && rm %s' % (out,src,clonepath,out)
      subprocess.check_output(cmd, shell = True, stderr = subprocess.STDOUT).decode('UTF-8')
      try:
         f = open(buildpath+"/fftw/fftw.h.in", "r")
         contents = f.readlines()
         f.close()
         contents.insert(23, '#define fftw_destroy_plan fftw_destroy_planK\n')
         f = open(buildpath+"/fftw/fftw.h.in", "w")
         contents = "".join(contents)
         f.write(contents)
         f.close()
      except subprocess.CalledProcessError as e:
         print("Kludge failed with:\n %s" % e.output.decode('UTF-8'))
         sys.exit(1)

# *** Build the library ********************************************** #

cmd = 'cd %s && ./configure --prefix=%s CC=mpicc CFLAGS="-O3" && make && make install' % (buildpath,libprefix)

if buildflag:
   print("Building the FFTW2 library ...")
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

cmd = 'cd %s && make clean' % (buildpath)

if cleanflag:
   print("Cleaning the FFTW2 build space ...")
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

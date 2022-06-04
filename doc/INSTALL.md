To build Socorro, execute the following:

0) Required Packages/Compilers
-----------------------------------------------------------
Socorro requires both a Fortran 90/95 compiler and a C compiler.
There are also several other packages required:
    -FFTW 2.1.5 (http://www.fftw.org).  Newer versions
          won't work.  We're in the process of switching
          the interface to support 3.0 and above.
    -BLAS library
    -LAPACK library
    -Gnu make or equivalent
    -Perl
    -Some kind of MPI


1) Create a make.conf
-----------------------------------------------------------
You will need to manually copy a make.conf from
the makefiles directory to the socorro root directory
(here).  After that edit it to your hearts content.

Most of the variables are pretty obvious what they stand for
but two critical variables MPLIBS and UPLIBS aren't.  The
first MPLIBS (Multi-Processor Libraries) should contain all
the libraries needed for socorro, including MPI.  The other,
UPLIBS (UniProcessor Libraries), should NOT contain anything
related to MPI.

Once you've made a stab at configuring make.conf
you are ready to continue.

2) Run configure
----------------------------------------------------------
Currently configure creates the files cpointer_mod.f90 and
ctof_io.h.  These file take into account the endian-ness 
of your machine and provide the C <-> Fortran interface
passing data. It also creates all the local Makefiles.

To skip making the C <-> F90 interface add the option

   --no-interface

when invoking configure.  You will need to manually create
or copy cpointer_mod.f90 and ctof_io.h in the src directory
before running configure.

Configure is run by typing:

  ./configure

with the optional option mentioned above.

If this works fine you can skip to step 3.

Cross-compiling or C<->F90 problems
--------------------------------------------------
If you are cross-compiling or mci doesn't work you will 
need to manually run the "mci" program to generate the 
2 files mentioned above.  You can do this by first 
changing to the tools/config directory and typing

  cc -o mci make_c_interface.c

on the compile machine. Replace "cc" with your C 
compiler and any other misc options you need.

Next run the executable "mci" on the destination
machine.  The machine you are going to be running 
socorro on. This will create the 2 files
cpointer_mod.f90 and ctof_io.h.  These files
should then be copied to the socorro source
directory "src".

After you've done this re-run configure skipping the
interface part:

   ./configure --no-interface

This will create the appropriate file links in the
different directories.

3) Build socrroro
-------------------------------------------------

Just type

   make

and hope for the best.  This will build socorro and 
also the routine taginfo located in tools/taginfo.




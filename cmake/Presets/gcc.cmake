# *** Presets that will explicitly request GCC compilers

set( CMAKE_C_COMPILER       "gcc-11"      CACHE STRING "" FORCE )
set( CMAKE_CXX_COMPILER     "g++-11"      CACHE STRING "" FORCE )
set( CMAKE_Fortran_COMPILER "gfortran-11" CACHE STRING "" FORCE )

set( MPI_C_COMPILER         "/usr/local/Cellar/open-mpi/4.1.3/bin/mpicc"   CACHE STRING "" FORCE )
set( MPI_CXX_COMPILER       "/usr/local/Cellar/open-mpi/4.1.3/bin/mpicxx"  CACHE STRING "" FORCE )
set( MPI_Fortran_COMPILER   "/usr/local/Cellar/open-mpi/4.1.3/bin/mpifort" CACHE STRING "" FORCE )

# *** 

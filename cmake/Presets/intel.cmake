# *** Presets that will explicitly request Intel compilers

set( CMAKE_C_COMPILER       "icc"     CACHE STRING "" FORCE )
set( CMAKE_CXX_COMPILER     "icpc"    CACHE STRING "" FORCE )
set( CMAKE_Fortran_COMPILER "ifort"   CACHE STRING "" FORCE )

set( MPI_C_COMPILER         "mpicc"   CACHE STRING "" FORCE )
set( MPI_CXX_COMPILER       "mpicxx"  CACHE STRING "" FORCE )
set( MPI_Fortran_COMPILER   "mpifort" CACHE STRING "" FORCE )

# *** 

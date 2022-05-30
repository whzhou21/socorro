# *** CMake version requirements ************************************* #

cmake_minimum_required( VERSION 3.20 )

# *** Project name and language requirements ************************* #

project( socorro Fortran C CXX )

set( CMAKE_C_STANDARD 11 )
set( CMAKE_C_STANDARD_REQUIRED ON )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# *** Project directories ******************************************** #

get_filename_component( SOCORRO_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE )

set( SOCORRO_BLD_DIR ${CMAKE_CURRENT_BINARY_DIR} )
set( SOCORRO_DLC_DIR ${SOCORRO_BLD_DIR}/Downloads )
set( SOCORRO_DOC_DIR ${SOCORRO_DIR}/doc )
set( SOCORRO_LIB_DIR ${SOCORRO_DIR}/lib )
set( SOCORRO_SRC_DIR ${SOCORRO_DIR}/src )

set( CMAKE_Fortran_MODULE_DIRECTORY ${SOCORRO_BLD_DIR}/Modules )
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules )

include( PreventInSourceBuilds )

# *** Project source files and target ******************************** #

set( BUILD_NAME "socorro.x" CACHE STRING "Name of the executable file (Does not enable any features)" )
if ( BUILD_NAME )
   set( SOCORRO_BINARY "${BUILD_NAME}" )
else ()
   set( SOCORRO_BINARY "socorro.x" )
endif()

file( GLOB SOURCE_FILES ${SOCORRO_SRC_DIR}/*.F90 )
add_executable( ${SOCORRO_BINARY} ${SOURCE_FILES} )

# *** Project build options ****************************************** #

set( CMAKE_Fortran_FLAGS_DEBUG          "-O0 -g -Wall" )
set( CMAKE_Fortran_FLAGS_RELEASE        "-O3" )
set( CMAKE_Fortran_FLAGS_MINSIZEREL     "-Os" )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g" )

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
   set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE )
endif()

if ( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS )
   set( CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the build type, options are: Debug Release RelWithDebInfo MinSizeRel" FORCE )
   set( CMAKE_BUILD_TYPE_VALUES None Debug Release RelWithDebInfo MinSizeRel )
   set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_BUILD_TYPE_VALUES} )
endif()

if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" )
   set( CMAKE_TUNE_DEFAULT "-xHost" )
elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" )
   set( CMAKE_TUNE_DEFAULT "" )
   target_compile_options( ${SOCORRO_BINARY} PRIVATE -ffree-line-length-none )
else ()
   set( CMAKE_TUNE_DEFAULT "" )
endif()

set( CMAKE_TUNE_FLAGS "${CMAKE_TUNE_DEFAULT}" CACHE STRING "Compiler and machine specific optimization flags (compilation only)" )

separate_arguments( CMAKE_TUNE_FLAGS )
foreach( TUNE_FLAG ${CMAKE_TUNE_FLAGS} )
   target_compile_options( ${SOCORRO_BINARY} PRIVATE ${TUNE_FLAG} )
endforeach()

# *** Build the project with MPI support ***************************** #

find_package( MPI REQUIRED )
target_link_libraries( ${SOCORRO_BINARY} MPI::MPI_Fortran )

# *** Build the project with OMP support ***************************** #

find_package( OpenMP QUIET )
option( BUILD_WITH_OMP "Build with OpenMP support" ${OpenMP_FOUND} )

if ( BUILD_WITH_OMP )
   find_package( OpenMP REQUIRED )
   target_link_libraries( ${SOCORRO_BINARY} OpenMP::OpenMP_Fortran )
   target_compile_definitions( ${SOCORRO_BINARY} PRIVATE -D_USE_OMP_ )
endif()

# *** Build the project with FFT support ***************************** #

include( Packages/FFTW3 )
include( Packages/FFTMPI )

# *** Build the project with LAPACK support ************************** #

include( Packages/LAPACK )

# *** Build the project with LIBXC support *************************** #

include( Packages/LIBXC )

# *** Build the project with SCALAPACK support *********************** #

option( BUILD_WITH_SCALAPACK "Build with SCALAPACK support" OFF )

if ( BUILD_WITH_SCALAPACK )
   include( Packages/SCALAPACK )
else ()
   unset( PKG_SCALAPACK CACHE )
   unset( PKG_SCALAPACK_PATH CACHE )
   unset( PKG_SCALAPACK_VERSION CACHE )
endif()

# *** Print a summary of the project build configuration ************* #

include( FeatureSummary )

if ( CMAKE_BUILD_TYPE )
   string( TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE )
endif()

get_target_property( OPTIONS ${SOCORRO_BINARY} COMPILE_OPTIONS )
if ( NOT OPTIONS )
   set( OPTIONS "" )
endif()

get_target_property( DEFINES ${SOCORRO_BINARY} COMPILE_DEFINITIONS )
if ( NOT DEFINES )
   set( DEFINES "" )
endif()

feature_summary( DESCRIPTION "\nThe following packages have been found and configured:\n" WHAT PACKAGES_FOUND )

message(">> Build Configuration:
   Core OS:       ${CMAKE_SYSTEM_NAME}
   Generator:     ${CMAKE_GENERATOR} using ${CMAKE_MAKE_PROGRAM}
   Build type:    ${CMAKE_BUILD_TYPE}
   Executable:    ${SOCORRO_BINARY}
\n")

message(">> Fortran Summary:
   Path:          ${CMAKE_Fortran_COMPILER}
   Version:       ${CMAKE_Fortran_COMPILER_VERSION}
   Vendor:        ${CMAKE_Fortran_COMPILER_ID}
")
if ( CMAKE_BUILD_TYPE AND BUILD_WITH_OMP )
   message("   Options:       ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}} ${OPTIONS} ${OpenMP_Fortran_FLAGS}\n")
else ()
   message("   Options:       ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}} ${OPTIONS} \n")
endif()
message("   Definitions:   ${DEFINES}\n\n")

message(">> C Summary:
   Path:          ${CMAKE_C_COMPILER}
   Version:       ${CMAKE_C_COMPILER_VERSION}
   Vendor:        ${CMAKE_C_COMPILER_ID}
\n")

message(">> CXX Summary:
   Path:          ${CMAKE_CXX_COMPILER}
   Version:       ${CMAKE_CXX_COMPILER_VERSION}
   Vendor:        ${CMAKE_CXX_COMPILER_ID}
\n")

# *** End of the CMakeLists.txt build file *************************** #

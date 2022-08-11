# *** CMake version requirements ************************************* #

cmake_minimum_required( VERSION 3.20 )

# *** Project name and language requirements ************************* #

project( socorro Fortran C CXX )

set( CMAKE_C_FLAGS "-O3" )
set( CMAKE_C_STANDARD 11 )
set( CMAKE_C_STANDARD_REQUIRED ON )

set( CMAKE_CXX_FLAGS "-O3" )
set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# *** Project directories ******************************************** #

get_filename_component( SOCORRO_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE )

set( SOCORRO_DOC_DIR ${SOCORRO_DIR}/doc )
set( SOCORRO_LIB_DIR ${SOCORRO_DIR}/lib )
set( SOCORRO_SRC_DIR ${SOCORRO_DIR}/src )

set( SOCORRO_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} )
set( SOCORRO_DOWNLOADS_DIR ${SOCORRO_BUILD_DIR}/Downloads )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules )
set( CMAKE_Fortran_MODULE_DIRECTORY ${SOCORRO_BUILD_DIR}/Modules )

# *** Project source files and executable **************************** #

set( BUILD_NAME "socorro" CACHE STRING "Name of the executable file (Does not enable any features)" )
set( SOCORRO_EXE "${BUILD_NAME}" )

file( GLOB SOURCE_FILES ${SOCORRO_SRC_DIR}/*.F90 )
add_executable( ${SOCORRO_EXE} ${SOURCE_FILES} )

file( GLOB SOURCE_CLIBS ${SOCORRO_SRC_DIR}/*.c )
add_library( clib ${SOURCE_CLIBS} )

add_dependencies( ${SOCORRO_EXE} clib )
target_link_libraries( ${SOCORRO_EXE} clib )

# *** Project build options ****************************************** #

include( PreventInSourceBuilds )

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
   set( CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE )
endif()

set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall" )
set( CMAKE_Fortran_FLAGS_MINSIZEREL "-Os" )
set( CMAKE_Fortran_FLAGS_NONE "" )
set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
set( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g" )

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the build type, options are: Debug Release RelWithDebInfo MinSizeRel" FORCE )
   set( CMAKE_BUILD_TYPE_VALUES None Debug Release RelWithDebInfo MinSizeRel )
   set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_BUILD_TYPE_VALUES} )
endif()

if ( (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel") OR (${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM") )
   set( CMAKE_TUNE_DEFAULT "-xHost" )
elseif ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" )
   set( CMAKE_TUNE_DEFAULT "" )
   target_compile_options( ${SOCORRO_EXE} PRIVATE -ffree-line-length-none -fallow-argument-mismatch )
else ()
   set( CMAKE_TUNE_DEFAULT "" )
endif()

set( CMAKE_TUNE_FLAGS "${CMAKE_TUNE_DEFAULT}" CACHE STRING "Compiler and machine specific optimization flags (compilation only)" )

separate_arguments( CMAKE_TUNE_FLAGS )
foreach( TUNE_FLAG ${CMAKE_TUNE_FLAGS} )
   target_compile_options( ${SOCORRO_EXE} PRIVATE ${TUNE_FLAG} )
endforeach()

# *** Build the project with required MPI support ******************** #

find_package( MPI REQUIRED )
target_link_libraries( clib MPI::MPI_C )
target_link_libraries( ${SOCORRO_EXE} MPI::MPI_Fortran )

# *** Build the project with optional OMP support ******************** #

find_package( OpenMP QUIET )
option( BUILD_WITH_OMP "Build with OpenMP support" ${OpenMP_FOUND} )
if ( BUILD_WITH_OMP )
   find_package( OpenMP REQUIRED )
   target_link_libraries( ${SOCORRO_EXE} OpenMP::OpenMP_Fortran )
endif()

# *** Build the project with the required libraries ****************** #

message( "\n" )

include( Packages/FFTW3 )
include( Packages/LIBXC )
include( Packages/LAPACK )
include( Packages/SCALAPACK )

# *** Print a summary of the project build configuration ************* #

include( FeatureSummary )

if ( CMAKE_BUILD_TYPE )
   string( TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE )
endif()

get_target_property( OPTIONS ${SOCORRO_EXE} COMPILE_OPTIONS )
if ( NOT OPTIONS )
   set( OPTIONS "" )
endif()

get_target_property( DEFINES ${SOCORRO_EXE} COMPILE_DEFINITIONS )
if ( NOT DEFINES )
   set( DEFINES "" )
endif()

feature_summary( DESCRIPTION "\nThe following packages have been found and configured:\n" WHAT PACKAGES_FOUND )

message(">> Build Configuration:
   Core OS:       ${CMAKE_SYSTEM_NAME}
   Generator:     ${CMAKE_GENERATOR} using ${CMAKE_MAKE_PROGRAM}
   Build type:    ${CMAKE_BUILD_TYPE}
   Executable:    ${SOCORRO_EXE}
\n")

message(">> Fortran Summary:
   Path:          ${CMAKE_Fortran_COMPILER}
   Version:       ${CMAKE_Fortran_COMPILER_VERSION}
   Vendor:        ${CMAKE_Fortran_COMPILER_ID}
   Options:       ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}} ${OPTIONS} ${DEFINES}
\n")

message(">> C Summary:
   Path:          ${CMAKE_C_COMPILER}
   Version:       ${CMAKE_C_COMPILER_VERSION}
   Vendor:        ${CMAKE_C_COMPILER_ID}
   Options:       ${CMAKE_C_FLAGS}
\n")

message(">> CXX Summary:
   Path:          ${CMAKE_CXX_COMPILER}
   Version:       ${CMAKE_CXX_COMPILER_VERSION}
   Vendor:        ${CMAKE_CXX_COMPILER_ID}
   Options:       ${CMAKE_CXX_FLAGS}
\n")

# *** End of the file ************************************************ #

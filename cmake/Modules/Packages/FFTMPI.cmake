# *** 

set( PKG_FFTMPI "DownloadSource" CACHE STRING "Choose the library retrieval type, options are: CMakeFindPackage DownloadSource ExistingSource" )
set( PKG_OPTIONS CMakeFindPackage DownloadSource ExistingSource )
set_property( CACHE PKG_FFTMPI PROPERTY STRINGS ${PKG_OPTIONS} )

# *** 

if ( PKG_FFTMPI STREQUAL "CMakeFindPackage" )

   message( STATUS "Existing package requested - we will link to a pre-compiled FFTMPI" )

elseif ( PKG_FFTMPI STREQUAL "DownloadSource" )

   message( STATUS "Download package requested - we will build our own FFTMPI" )
   unset( PKG_FFTMPI_PATH CACHE )
   set( PKG_FFTMPI_VERSION "master" CACHE STRING "Version of the library to be downloaded (Specify a git tag)" )

   include( ExternalProject )
   ExternalProject_Add( fftmpi
      PREFIX            ${SOCORRO_DLC_DIR}/fftmpi-${PKG_FFTMPI_VERSION}
      GIT_REPOSITORY    "https://github.com/lammps/fftmpi.git"
      GIT_TAG           ${PKG_FFTMPI_VERSION}
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      SOURCE_SUBDIR     src/
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     make lib fft=FFTW3 CC=${MPI_CXX_COMPILER}
      BUILD_IN_SOURCE   YES
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( fftmpi BINARY_DIR )

   file( GLOB FFTMPI_SOURCES ${BINARY_DIR}/fft3d_wrap.f90 )
   target_sources( ${SOCORRO_BINARY} PRIVATE ${FFTMPI_SOURCES} )

   add_library( SOCORRO::FFTMPI UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::FFTMPI PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libfft3dmpi.a" )
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::FFTMPI )
   add_dependencies( SOCORRO::FFTMPI fftmpi )

elseif ( PKG_FFTMPI STREQUAL "ExistingSource" )

   message( STATUS "Existing package requested - we will build our own FFTMPI" )
   unset( PKG_FFTMPI_VERSION CACHE )
   set( PKG_FFTMPI_PATH "${SOCORRO_LIB_DIR}/fftmpi/" CACHE STRING "Path to the top-level FFTMPI folder to be built" )

   include( ExternalProject )
   ExternalProject_Add( fftmpi
      DOWNLOAD_COMMAND  ""
      SOURCE_DIR        ${PKG_FFTMPI_PATH}
      PREFIX            ${SOCORRO_BLD_DIR}/build_fftmpi
#      CMAKE_ARGS        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
#                        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
#                        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      INSTALL_COMMAND   ""
   )
#   ExternalProject_Get_Property( libxc BINARY_DIR )

#   add_library( SOCORRO::LIBXCF03 UNKNOWN IMPORTED )
#   set_target_properties( SOCORRO::LIBXCF03 PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libxcf03.a" )
#   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LIBXCF03 )
#   add_dependencies( SOCORRO::LIBXCF03 libxc )
#   target_compile_definitions( ${SOCORRO_BINARY} PRIVATE -DUSE_LIBXC )

endif()

# *** 

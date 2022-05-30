# *** 

set( PKG_FFTW3 "DownloadSource" CACHE STRING "Choose the library retrieval type, options are: CMakeFindPackage DownloadSource ExistingSource" )
set( PKG_OPTIONS CMakeFindPackage DownloadSource ExistingSource )
set_property( CACHE PKG_FFTW3 PROPERTY STRINGS ${PKG_OPTIONS} )

# *** 

if ( PKG_FFTW3 STREQUAL "CMakeFindPackage" )

   message( STATUS "Existing package requested - we will link to a pre-compiled FFTW3" )

elseif ( PKG_FFTW3 STREQUAL "DownloadSource" )

   message( STATUS "Download package requested - we will build our own FFTW3" )
   unset( PKG_FFTW3_PATH CACHE )
   set( PKG_FFTW3_VERSION "3.3.10" CACHE STRING "Version of the library to be downloaded (Specify a version)" )

   include( ExternalProject )
   ExternalProject_Add( fftw3
      PREFIX            ${SOCORRO_DLC_DIR}/fftw-${PKG_FFTW3_VERSION}
      URL               https://fftw.org/pub/fftw/fftw-${PKG_FFTW3_VERSION}.tar.gz
      URL_MD5           "8ccbf6a5ea78a16dbc3e1306e234cc5c"
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CMAKE_ARGS        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                        -DBUILD_SHARED_LIBS=OFF
                        -DBUILD_TESTS=OFF
      BUILD_IN_SOURCE   NO
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( fftw3 BINARY_DIR )

   add_library( SOCORRO::FFTW3 UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::FFTW3 PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libfftw3.a" )
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::FFTW3 )
   add_dependencies( SOCORRO::FFTW3 fftw3 )

   target_include_directories( ${SOCORRO_BINARY} PUBLIC ${BINARY_DIR} )

elseif ( PKG_FFTW3 STREQUAL "ExistingSource" )

   message( STATUS "Existing package requested - we will build our own FFTW3" )
   unset( PKG_FFTW3_VERSION CACHE )
   set( PKG_FFTW3_PATH "${SOCORRO_LIB_DIR}/fftw/" CACHE STRING "Path to the top-level FFTMPI folder to be built" )

   include( ExternalProject )
   ExternalProject_Add( fftw3
      DOWNLOAD_COMMAND  ""
      SOURCE_DIR        ${PKG_FFTW3_PATH}
#      PREFIX            ${SOCORRO_BLD_DIR}/build_fftmpi
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

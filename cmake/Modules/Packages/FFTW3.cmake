# *** Initialize package options ************************************* #

set( PKG_FFTW3 "Download" CACHE STRING "Choose the library retrieval method, options are: Download Existing Find_PKG" )
set( PKG_OPTIONS Download Existing Find_PKG )
set_property( CACHE PKG_FFTW3 PROPERTY STRINGS ${PKG_OPTIONS} )

# *** Package options ************************************************ #

if ( PKG_FFTW3 STREQUAL "Download" )

   message( STATUS "Download package requested - we will download, compile, and link to FFTW3 (Version: 3.3.10)" )

   unset( PKG_FFTW3_PREFIX CACHE )

   list( APPEND FFTW3_BUILD_ARGS "-DBUILD_SHARED_LIBS=OFF" )
   list( APPEND FFTW3_BUILD_ARGS "-DBUILD_TESTS=OFF" )
   list( APPEND FFTW3_BUILD_ARGS "-DCMAKE_BUILD_TYPE=Release" )
   list( APPEND FFTW3_BUILD_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}" )
   list( APPEND FFTW3_BUILD_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}" )
   list( APPEND FFTW3_BUILD_ARGS "-DCMAKE_INSTALL_LIBDIR=lib" )
   list( APPEND FFTW3_BUILD_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>" )

   if ( BUILD_WITH_OMP )
      list( APPEND FFTW3_BUILD_ARGS "-DENABLE_THREADS=ON")
   else ()
      list( APPEND FFTW3_BUILD_ARGS "-DENABLE_THREADS=OFF")
   endif()

   include( ExternalProject )
   ExternalProject_Add( fftw3
      PREFIX            ${SOCORRO_DOWNLOADS_DIR}/fftw-3.3.10
      URL               https://fftw.org/pub/fftw/fftw-3.3.10.tar.gz
      URL_MD5           "8ccbf6a5ea78a16dbc3e1306e234cc5c"
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CMAKE_ARGS        ${FFTW3_BUILD_ARGS}
      BUILD_IN_SOURCE   NO
   )
   ExternalProject_Get_Property( fftw3 INSTALL_DIR )
   file( MAKE_DIRECTORY ${INSTALL_DIR}/include )

   add_library( SOCORRO::FFTW3 UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::FFTW3 PROPERTIES
      IMPORTED_LOCATION "${INSTALL_DIR}/lib/libfftw3.a"
      INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
   )
   add_dependencies( SOCORRO::FFTW3 fftw3 )
   target_link_libraries( clib SOCORRO::FFTW3 )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::FFTW3 )

   if ( BUILD_WITH_OMP )
      add_library( SOCORRO::FFTW3_THREADS UNKNOWN IMPORTED )
      set_target_properties( SOCORRO::FFTW3_THREADS PROPERTIES
         IMPORTED_LOCATION "${INSTALL_DIR}/lib/libfftw3_threads.a"
      )
      add_dependencies( SOCORRO::FFTW3_THREADS fftw3 )
      target_link_libraries( ${SOCORRO_EXE} SOCORRO::FFTW3_THREADS )
   endif()

elseif ( PKG_FFTW3 STREQUAL "Existing" )

   message( STATUS "Existing package requested - we will link to a user-specified and pre-compiled FFTW3" )

   set( PKG_FFTW3_PREFIX "${SOCORRO_LIB_DIR}/fftw3/fftw3/bin" CACHE STRING "Absolute path to the FFTW3 installation directory" )

   add_library( SOCORRO::FFTW3 UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::FFTW3 PROPERTIES
      IMPORTED_LOCATION "${PKG_FFTW3_PREFIX}/lib/libfftw3.a"
      INTERFACE_INCLUDE_DIRECTORIES "${PKG_FFTW3_PREFIX}/include"
   )
   target_link_libraries( clib SOCORRO::FFTW3 )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::FFTW3 )

   if ( BUILD_WITH_OMP )
      add_library( SOCORRO::FFTW3_THREADS UNKNOWN IMPORTED )
      set_target_properties( SOCORRO::FFTW3_THREADS PROPERTIES
         IMPORTED_LOCATION "${PKG_FFTW3_PREFIX}/lib/libfftw3_threads.a"
      )
      target_link_libraries( ${SOCORRO_EXE} SOCORRO::FFTW3_THREADS )
   endif()

elseif ( PKG_FFTW3 STREQUAL "Find_PKG" )

   message( STATUS "Existing package requested - we will search for and link to a pre-compiled FFTW3" )

   unset( PKG_FFTW3_PREFIX CACHE )

endif()

# *** End of the file ************************************************ #

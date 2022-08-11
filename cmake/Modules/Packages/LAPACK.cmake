# *** Initialize package options ************************************* #

set( PKG_LAPACK "Download" CACHE STRING "Choose the library retrieval method, options are: Download Existing Find_PKG" )
set( PKG_OPTIONS Download Existing Find_PKG )
set_property( CACHE PKG_LAPACK PROPERTY STRINGS ${PKG_OPTIONS} )

# *** Package options ************************************************ #

if ( PKG_LAPACK STREQUAL "Download" )

   message( STATUS "Download package requested - we will download, compile, and link to LAPACK (Version: 3.10.1)" )

   unset( PKG_LAPACK_PREFIX CACHE )

   list( APPEND LAPACK_BUILD_ARGS "-DCMAKE_BUILD_TYPE=Release" )
   list( APPEND LAPACK_BUILD_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}" )
   list( APPEND LAPACK_BUILD_ARGS "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}" )
   list( APPEND LAPACK_BUILD_ARGS "-DCMAKE_INSTALL_LIBDIR=lib" )
   list( APPEND LAPACK_BUILD_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>" )

   include( ExternalProject )
   ExternalProject_Add( lapack
      PREFIX            ${SOCORRO_DOWNLOADS_DIR}/lapack-v3.10.1
      GIT_REPOSITORY    "https://github.com/Reference-LAPACK/lapack.git"
      GIT_TAG           v3.10.1
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CMAKE_ARGS        ${LAPACK_BUILD_ARGS}
      BUILD_IN_SOURCE   NO
   )
   ExternalProject_Get_Property( lapack INSTALL_DIR )
   file( MAKE_DIRECTORY ${INSTALL_DIR}/include )

   add_library( SOCORRO::BLAS UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::BLAS PROPERTIES
      IMPORTED_LOCATION "${INSTALL_DIR}/lib/libblas.a"
   )
   add_dependencies( SOCORRO::BLAS lapack )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::BLAS )

   add_library( SOCORRO::LAPACK UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LAPACK PROPERTIES
      IMPORTED_LOCATION "${INSTALL_DIR}/lib/liblapack.a"
   )
   add_dependencies( SOCORRO::LAPACK lapack )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::LAPACK )

elseif ( PKG_LAPACK STREQUAL "Existing" )

   message( STATUS "Existing package requested - we will link to a user-specified and pre-compiled LAPACK" )

   set( PKG_LAPACK_PREFIX "${SOCORRO_LIB_DIR}/lapack/lapack" CACHE STRING "Absolute path to the LAPACK installation directory" )

   add_library( SOCORRO::BLAS UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::BLAS PROPERTIES
      IMPORTED_LOCATION "${PKG_LAPACK_PREFIX}/librefblas.a"
   )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::BLAS )

   add_library( SOCORRO::LAPACK UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LAPACK PROPERTIES
      IMPORTED_LOCATION "${PKG_LAPACK_PREFIX}/liblapack.a"
   )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::LAPACK )

elseif ( PKG_LAPACK STREQUAL "Find_PKG" )

   message( STATUS "Existing package requested - we will search for and link to a pre-compiled LAPACK" )

   unset( PKG_LAPACK_PREFIX CACHE )

endif()

# *** End of the file ************************************************ #

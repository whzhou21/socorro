# *** Initialize package options ************************************* #

set( PKG_LIBXC "Download" CACHE STRING "Choose the library retrieval type, options are: Download Existing Find_PKG" )
set( PKG_OPTIONS Download Existing Find_PKG )
set_property( CACHE PKG_LIBXC PROPERTY STRINGS ${PKG_OPTIONS} )

# *** Package options ************************************************ #

if ( PKG_LIBXC STREQUAL "Download" )

   message( STATUS "Download package requested - we will download, compile, and link to LIBXC (Version: 4.3.4)" )

   unset( PKG_LIBXC_PREFIX CACHE )

   list( APPEND LIBXC_BUILD_ARGS "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}" )
   list( APPEND LIBXC_BUILD_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}" )
   list( APPEND LIBXC_BUILD_ARGS "-DBUILD_TESTING=OFF" )
   list( APPEND LIBXC_BUILD_ARGS "-DCMAKE_BUILD_TYPE=Release" )
   list( APPEND LIBXC_BUILD_ARGS "-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>" )
   list( APPEND LIBXC_BUILD_ARGS "-DENABLE_FORTRAN=ON" )
   list( APPEND LIBXC_BUILD_ARGS "-DENABLE_XHOST=OFF" )

   include( ExternalProject )
   ExternalProject_Add( libxc
      PREFIX            ${SOCORRO_DOWNLOADS_DIR}/libxc-4.3.4
      GIT_REPOSITORY    "https://gitlab.com/libxc/libxc.git"
      GIT_TAG           4.3.4
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CMAKE_ARGS        ${LIBXC_BUILD_ARGS}
      BUILD_IN_SOURCE   NO
   )
   ExternalProject_Get_Property( libxc INSTALL_DIR )
   file( MAKE_DIRECTORY ${INSTALL_DIR}/include )

   add_library( SOCORRO::LIBXC UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LIBXC PROPERTIES
      IMPORTED_LOCATION "${INSTALL_DIR}/lib/libxc.a"
      INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
   )
   add_dependencies( SOCORRO::LIBXC libxc )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::LIBXC )

   add_library( SOCORRO::LIBXCF90 UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LIBXCF90 PROPERTIES
      IMPORTED_LOCATION "${INSTALL_DIR}/lib/libxcf90.a"
   )
   add_dependencies( SOCORRO::LIBXCF90 libxc )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::LIBXCF90 )

elseif ( PKG_LIBXC STREQUAL "Existing" )

   message( STATUS "Existing package requested - we will link to a user-specified and pre-compiled LIBXC" )

   set( PKG_LIBXC_PREFIX "${SOCORRO_LIB_DIR}/libxc/libxc/bin" CACHE STRING "Absolute path to the LIBXC installation directory" )

   add_library( SOCORRO::LIBXC UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LIBXC PROPERTIES
      IMPORTED_LOCATION "${PKG_LIBXC_PREFIX}/lib/libxc.a"
      INTERFACE_INCLUDE_DIRECTORIES "${PKG_LIBXC_PREFIX}/include"
   )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::LIBXC )

   add_library( SOCORRO::LIBXCF90 UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LIBXCF90 PROPERTIES
      IMPORTED_LOCATION "${PKG_LIBXC_PREFIX}/lib/libxcf90.a"
   )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::LIBXC )

elseif ( PKG_LIBXC STREQUAL "Find_PKG" )

   message( STATUS "Existing package requested - we will search for and link to a pre-compiled LIBXC" )

   unset( PKG_LIBXC_PREFIX CACHE )

endif()

# *** End of the file ************************************************ #

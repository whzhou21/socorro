# *** 

set( PKG_LIBXC "DownloadSource" CACHE STRING "Choose the library retrieval type, options are: CMakeFindPackage DownloadSource ExistingSource" )
set( PKG_OPTIONS CMakeFindPackage DownloadSource ExistingSource )
set_property( CACHE PKG_LIBXC PROPERTY STRINGS ${PKG_OPTIONS} )

# *** 

if ( PKG_LIBXC STREQUAL "CMakeFindPackage" )

   message( STATUS "Existing package requested - we will link to a pre-compiled LIBXC" )

elseif ( PKG_LIBXC STREQUAL "DownloadSource" )

   message( STATUS "Download package requested - we will build our own LIBXC" )
   unset( PKG_LIBXC_PATH CACHE )
   set( PKG_LIBXC_VERSION "5.1.7" CACHE STRING "Version of the library to be downloaded (Specify a git tag)" )

   include( ExternalProject )
   ExternalProject_Add( libxc
      PREFIX            ${SOCORRO_DLC_DIR}/libxc-${PKG_LIBXC_VERSION}
      GIT_REPOSITORY    "https://gitlab.com/libxc/libxc.git"
      GIT_TAG           ${PKG_LIBXC_VERSION}
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CMAKE_ARGS        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                        -DBUILD_TESTING=OFF
                        -DENABLE_FORTRAN=ON
                        -DENABLE_XHOST=OFF
      BUILD_IN_SOURCE   NO
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( libxc BINARY_DIR )

   add_library( SOCORRO::LIBXC UNKNOWN IMPORTED )
   add_library( SOCORRO::LIBXCF03 UNKNOWN IMPORTED )

   set_target_properties( SOCORRO::LIBXC PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libxc.a" )
   set_target_properties( SOCORRO::LIBXCF03 PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libxcf03.a" )

   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LIBXC )
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LIBXCF03 )

   add_dependencies( SOCORRO::LIBXC libxc )
   add_dependencies( SOCORRO::LIBXCF03 libxc )

   target_include_directories( ${SOCORRO_BINARY} PUBLIC ${BINARY_DIR} )

elseif ( PKG_LIBXC STREQUAL "ExistingSource" )

   message( STATUS "Existing package requested - we will build our own LIBXC" )
   unset( PKG_LIBXC_VERSION CACHE )
   set( PKG_LIBXC_PATH "${SOCORRO_LIB_DIR}/libxc/" CACHE STRING "Path to the top-level LIBXC folder to be built" )

   include( ExternalProject )
   ExternalProject_Add( libxc
      DOWNLOAD_COMMAND  ""
      SOURCE_DIR        ${PKG_LIBXC_PATH}
      PREFIX            ${SOCORRO_BLD_DIR}/build_libxc
      CMAKE_ARGS        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DBUILD_FPIC=ON
                        -DBUILD_SHARED_LIBS=OFF
                        -DBUILD_TESTING=OFF
                        -DDISABLE_FXC=OFF
                        -DDISABLE_KXC=ON
                        -DDISABLE_LXC=ON
                        -DDISABLE_VXC=OFF
                        -DENABLE_FORTRAN=ON
                        -DENABLE_GENERIC=OFF
                        -DENABLE_PYTHON=OFF
                        -DENABLE_XHOST=OFF
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( libxc BINARY_DIR )

   add_library( SOCORRO::LIBXC UNKNOWN IMPORTED )
   add_library( SOCORRO::LIBXCF03 UNKNOWN IMPORTED )

   set_target_properties( SOCORRO::LIBXC PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libxc.a" )
   set_target_properties( SOCORRO::LIBXCF03 PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libxcf03.a" )

   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LIBXC )
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LIBXCF03 )

   add_dependencies( SOCORRO::LIBXC libxc )
   add_dependencies( SOCORRO::LIBXCF03 libxc )

   target_include_directories( ${SOCORRO_BINARY} PUBLIC ${BINARY_DIR} )
   target_compile_definitions( ${SOCORRO_BINARY} PRIVATE -DUSE_LIBXC )

endif()

# *** 

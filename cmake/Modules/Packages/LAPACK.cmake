# ***

set( PKG_LAPACK "DownloadSource" CACHE STRING "Choose the library retrieval type, options are: CMakeFindPackage DownloadSource ExistingSource" )
set( PKG_OPTIONS CMakeFindPackage DownloadSource ExistingSource )
set_property( CACHE PKG_LAPACK PROPERTY STRINGS ${PKG_OPTIONS} )

# *** 

if ( PKG_LAPACK STREQUAL "CMakeFindPackage" )

   message( STATUS "Existing package requested - we will link to a pre-compiled LAPACK" )

elseif ( PKG_LAPACK STREQUAL "DownloadSource" )

   message( STATUS "Download package requested - we will build our own LAPACK" )
   unset( PKG_LAPACK_PATH CACHE )
   set( PKG_LAPACK_VERSION "v3.10.1" CACHE STRING "Version of the library to be downloaded (Specify a git tag)" )

   include( ExternalProject )
   ExternalProject_Add( lapack
      PREFIX            ${SOCORRO_DLC_DIR}/lapack-${PKG_LAPACK_VERSION}
      GIT_REPOSITORY    "https://github.com/Reference-LAPACK/lapack.git"
      GIT_TAG           ${PKG_LAPACK_VERSION}
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CMAKE_ARGS        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      BUILD_IN_SOURCE   NO
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( lapack BINARY_DIR )

   add_library( SOCORRO::BLAS UNKNOWN IMPORTED )
   add_library( SOCORRO::LAPACK UNKNOWN IMPORTED )

   set_target_properties( SOCORRO::BLAS PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/lib/libblas.a" )
   set_target_properties( SOCORRO::LAPACK PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/lib/liblapack.a" )

   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::BLAS )
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LAPACK )

   add_dependencies( SOCORRO::BLAS lapack )
   add_dependencies( SOCORRO::LAPACK lapack )

elseif ( PKG_LAPACK STREQUAL "ExistingSource" )

   message( STATUS "Existing package requested - we will build our own LAPACK" )
   unset( PKG_LAPACK_VERSION CACHE )
   set( PKG_LAPACK_PATH "${SOCORRO_LIB_DIR}/lapack/" CACHE STRING "Path to the top-level LAPACK folder to be built" )

   include( ExternalProject )
   ExternalProject_Add( lapack
      DOWNLOAD_COMMAND  ""
      SOURCE_DIR        ${PKG_LAPACK_PATH}
      PREFIX            ${SOCORRO_BLD_DIR}/build_lapack
      CMAKE_ARGS        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
                        -DBUILD_SINGLE=OFF
                        -DBUILD_COMPLEX=OFF
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( lapack BINARY_DIR )

   add_library( SOCORRO::LAPACK UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::LAPACK PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/lib/liblapack.a" )
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LAPACK )
   add_dependencies( SOCORRO::LAPACK lapack )
   target_compile_definitions( ${SOCORRO_BINARY} PRIVATE -DUSE_LAPACK )

endif()

# *** 

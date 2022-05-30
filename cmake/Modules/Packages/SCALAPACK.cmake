# *** 

set( PKG_SCALAPACK "DownloadSource" CACHE STRING "Choose the library retrieval type, options are: CMakeFindPackage DownloadSource ExistingSource" )
set( PKG_OPTIONS CMakeFindPackage DownloadSource ExistingSource )
set_property( CACHE PKG_SCALAPACK PROPERTY STRINGS ${PKG_OPTIONS} )

# *** 

if ( PKG_SCALAPACK STREQUAL "CMakeFindPackage" )

   message( STATUS "Existing package requested - we will link to a pre-compiled SCALAPACK" )

elseif ( PKG_SCALAPACK STREQUAL "DownloadSource" )

   message( STATUS "Download package requested - we will build our own SCALAPACK" )
   unset( PKG_SCALAPACK_PATH CACHE )
   set( PKG_SCALAPACK_VERSION "v2.2.1" CACHE STRING "Version of the library to be downloaded (Specify a git tag)" )

   include( ExternalProject )
   ExternalProject_Add( scalapack
      PREFIX            ${SOCORRO_DLC_DIR}/scalapack-${PKG_SCALAPACK_VERSION}
      GIT_REPOSITORY    "https://github.com/Reference-ScaLAPACK/scalapack.git"
      GIT_TAG           ${PKG_SCALAPACK_VERSION}
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/SLmake.inc.example <SOURCE_DIR>/SLmake.inc
      COMMAND           make lib CC=${MPI_C_COMPILER} FC=${MPI_Fortran_COMPILER}
                        CCFLAGS+="-Wno-error=implicit-function-declaration" FCFLAGS+="-fallow-argument-mismatch"
      BUILD_IN_SOURCE   YES
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( scalapack BINARY_DIR )

   add_library( SOCORRO::SCALAPACK UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::SCALAPACK PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libscalapack.a")
   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::SCALAPACK )
   add_dependencies( SOCORRO::SCALAPACK scalapack )

   target_compile_definitions( ${SOCORRO_BINARY} PRIVATE -D_USE_SCALAPACK_ )

elseif ( PKG_SCALAPACK STREQUAL "ExistingSource" )

   message( STATUS "Existing package requested - we will build our own SCALAPACK" )
   unset( PKG_SCALAPACK_VERSION CACHE )
   set( PKG_SCALAPACK_PATH "${SOCORRO_LIB_DIR}/scalapack/" CACHE STRING "Path to the top-level SCALAPACK folder to be built" )

   include( ExternalProject )
   ExternalProject_Add( scalapack
      DOWNLOAD_COMMAND  ""
      SOURCE_DIR        ${PKG_SCALAPACK_PATH}
      PREFIX            ${SOCORRO_BLD_DIR}/build_scalapack
#      CMAKE_ARGS        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
#                        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
#                        -DBUILD_TESTING=OFF
#                        -DENABLE_FORTRAN=ON
#                        -DENABLE_XHOST=OFF
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( scalapack BINARY_DIR )

#   add_library(SOCORRO::LIBXC UNKNOWN IMPORTED)
#   set_target_properties(SOCORRO::LIBXC PROPERTIES IMPORTED_LOCATION "${BINARY_DIR}/libxc.a")
#   target_link_libraries( ${SOCORRO_BINARY} SOCORRO::LIBXC )
#   add_dependencies(SOCORRO::LIBXC libxc)
   target_compile_definitions( ${SOCORRO_BINARY} PRIVATE -DUSE_SCALAPACK )

endif()

# *** 

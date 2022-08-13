# *** Initialize package options ************************************* #

set( PKG_SCALAPACK "Download" CACHE STRING "Choose the library retrieval type, options are: Download Existing Find_PKG" )
set( PKG_OPTIONS Download Existing Find_PKG )
set_property( CACHE PKG_SCALAPACK PROPERTY STRINGS ${PKG_OPTIONS} )

# *** Package options ************************************************ #

if ( PKG_SCALAPACK STREQUAL "Download" )

   message( STATUS "Download package requested - we will download, compile, and link to SCALAPACK (Version: 2.2.1)" )

   unset( PKG_SCALAPACK_PREFIX CACHE )

   if ( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0 )
      set ( FC_EXTRA "-fallow-argument-mismatch" )
   endif()

   ExternalProject_Add( scalapack
      PREFIX            ${SOCORRO_DOWNLOADS_DIR}/scalapack-v2.2.1
      GIT_REPOSITORY    "https://github.com/Reference-ScaLAPACK/scalapack.git"
      GIT_TAG           v2.2.1
      GIT_SHALLOW       YES
      GIT_PROGRESS      YES
      UPDATE_COMMAND    ""
      PATCH_COMMAND     ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     ${CMAKE_COMMAND} -E copy_if_different <SOURCE_DIR>/SLmake.inc.example <SOURCE_DIR>/SLmake.inc
      COMMAND           make lib CC=${MPI_C_COMPILER} FC=${MPI_Fortran_COMPILER}
                        CCFLAGS+="-Wno-error=implicit-function-declaration" FCFLAGS+="${FC_EXTRA}"
      BUILD_IN_SOURCE   YES
      INSTALL_COMMAND   ""
   )
   ExternalProject_Get_Property( scalapack BINARY_DIR )

   add_library( SOCORRO::SCALAPACK UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${BINARY_DIR}/libscalapack.a"
   )
   add_dependencies( SOCORRO::SCALAPACK scalapack )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::SCALAPACK )

elseif ( PKG_SCALAPACK STREQUAL "Existing" )

   message( STATUS "Existing package requested - we will link to a user-specified and pre-compiled SCALAPACK" )

   set( PKG_SCALAPACK_PREFIX "${SOCORRO_LIB_DIR}/scalapack/scalapack" CACHE STRING "Absolute path to the SCALAPACK installation directory" )

   add_library( SOCORRO::SCALAPACK UNKNOWN IMPORTED )
   set_target_properties( SOCORRO::SCALAPACK PROPERTIES
      IMPORTED_LOCATION "${PKG_SCALAPACK_PREFIX}/libscalapack.a"
   )
   target_link_libraries( ${SOCORRO_EXE} SOCORRO::SCALAPACK )

elseif ( PKG_SCALAPACK STREQUAL "Find_PKG" )

   message( STATUS "Existing package requested - we will search for and link to a pre-compiled SCALAPACK" )

   unset( PKG_SCALAPACK_PREFIX CACHE )

endif()

# *** End of the file ************************************************ #

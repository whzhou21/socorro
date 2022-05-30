!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!                                                                                                                                  !
!  Copyright (2003). See the README file in the top-level directory.                                                               !
!  This software is distributed with the GNU General Public License.                                                               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

module version_mod

   use kind_mod
   use utils_mod

   implicit none ; private

   !* Tag for the current stable version of Socorro

   character(line_len), parameter :: socorro_release = "(13 May 2022)"

   !* Publicly available parameters and procedures

   public :: x_version

   !* Interfaces for the publicly available procedures

   interface x_version
      module procedure x_version_
   end interface x_version

contains


   !* Method to retrieve the current stable version of Socorro

   function x_version_() result( v )

      character(len_trim(socorro_release)) :: v

      v = trimstr(socorro_release)

   end function x_version_


end module version_mod

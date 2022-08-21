!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module version_mod
!doc$ module version_mod

!     Comments ...

      use kind_mod
      use utils_mod

!cod$
      implicit none ; private

      character(line_len), parameter :: socorro_release = "(13 May 2022)"

!doc$
      public :: x_version

!cod$
      interface x_version
         module procedure x_version_
      end interface

      contains

! *** Public routines

      function x_version_() result( v )
!doc$ function x_version_()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         character(:), allocatable :: v

         v = trimstr(socorro_release)

      end function x_version_

      end module version_mod

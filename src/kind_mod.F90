!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module kind_mod
!doc$ module kind_mod

!     The module which holds kind info.

!cod$
      implicit none ; public

!doc$

      ! Data types

      integer, parameter :: tag_sz = 8
      integer, parameter :: line_len = 132

      integer, parameter :: single = kind(1.0e0)
      integer, parameter :: double = kind(1.0d0)

      integer, parameter :: long = selected_int_kind(9)
      integer, parameter :: longlong = selected_int_kind(16)

      integer, parameter :: sizeof_single = 4
      integer, parameter :: sizeof_double = 8

      integer, parameter :: sizeof_long = 4
      integer, parameter :: sizeof_longlong = 8

      integer, parameter :: CPTR = selected_int_kind(16)

!cod$
      end module

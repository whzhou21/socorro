!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module ghost_mod
!doc$ module ghost_mod

!     This module defines the unique identifier system that
!     is used to determine identity among various components in
!     Socorro.

!cod$
      implicit none
      private

      integer :: gid_current = 0

      type, public :: ghost
        integer :: gid
      end type

!doc$
      public :: operator(.eq.), operator(.ne.)
      public :: x_ghost

!cod$

      interface operator(.eq.)  ! also recognizes ==
        module procedure eq_g
      end interface
      interface operator(.ne.)  ! also recognizes /=
        module procedure ne_g
      end interface
      interface x_ghost
        module procedure get_ghost
      end interface

      contains

      function get_ghost() result(g)
!doc$ function x_ghost() result(g)
        type(ghost) :: g
!       effects: Returns a new unique ghost.

!cod$
        gid_current = gid_current + 1
        g%gid = gid_current
      end function

      function eq_g(g1,g2)
!doc$ operator(.eq.)(g1,g2)
        type(ghost), intent(in) :: g1, g2
        logical :: eq_g
!       effects: Tests ghosts g1 and g2 for equality.

!cod$
        eq_g = (g1%gid == g2%gid)
      end function

      function ne_g(g1,g2)
!doc$ operator(.ne.)(g1,g2)
        type(ghost), intent(in) :: g1, g2
        logical :: ne_g
!       effects: Tests ghosts g1 and g2 for inequality.

!cod$
        ne_g = (g1%gid /= g2%gid)
      end function

    end module

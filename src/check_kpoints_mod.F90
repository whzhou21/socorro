!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module check_kpoints_mod
!doc$ module check_kpoints_mod

!     Some info ...
!     ...

      use kind_mod
      use mpi_mod
      use error_mod
      use diary_mod
      use crystal_mod
      use symmetry_mod
      use kpoints_mod
      use timing_mod

!cod$
      implicit none ; private

      type(crystal_obj) :: cr
      type(point_group_obj) :: dg, lg
      type(space_group_obj) :: sg
      type(kpoints_obj) :: kp

!doc$
      public :: check_kpoints

!cod$
      interface check_kpoints
         module procedure check_kpoints_
      end interface

      contains

! public routines

      subroutine check_kpoints_()
!doc$ subroutine check_kpoints()
!        effects:
!        errors:
!        requires:

!cod$
         call start_timer("check_kpoints: total time")

         call my(crystal(),cr) ; if (error()) goto 100

         call my(point_group(x_lattice(cr)),lg) ; if (error()) goto 100
         call my(space_group(lg,x_atoms(cr),x_lattice(cr)),sg) ; if (error()) goto 100

         call my(point_group(sg,parity=.true.),dg) ; if (error()) goto 100
         call my(kpoints(x_lattice(cr),lg,dg),kp) ; if (error()) goto 100

         call diary(cr)
         call diary(lg)
         call diary(sg)
         call diary(kp)

         call glean(thy(cr))
         call glean(thy(lg))
         call glean(thy(sg))
         call glean(thy(dg))
         call glean(thy(kp))

100      if ( error("Exit check_kpoints") ) continue
         if ( .not.error() ) call stop_timer("check_kpoints: total time")

      end subroutine

      end module check_kpoints_mod

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module check_kpoints_mod
!doc$ module check_kpoints_mod

!     This check_kpoints module determines the symmetry for a set of
!     atomic coordinates in a parallelpiped and the Monkhorst-Pack
!     special k-points for a set of Monkhorst-Pack parameters.

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
      type(kpoints_obj) :: kp
      type(point_group_obj) :: dg, lg
      type(space_group_obj) :: sg

!doc$
      public :: check_kpoints

!cod$
      interface check_kpoints
         module procedure check_kpoints_
      end interface

      contains

! Public routines

      subroutine check_kpoints_()
!doc$ subroutine check_kpoints()
!        effects:
!        errors:
!        requires:

!cod$
         call start_timer("check_kpoints: total time")

         call my(crystal(),cr) ; if ( error() ) goto 900

         call my(point_group(x_lattice(cr)),lg) ; if ( error() ) goto 900
         call my(space_group(lg,x_atoms(cr),x_lattice(cr)),sg) ; if ( error() ) goto 900

         call my(point_group(sg,parity=.true.),dg) ; if ( error() ) goto 900
         call my(kpoints(x_lattice(cr),lg,dg),kp) ; if ( error() ) goto 900

         call diary(cr)
         call diary(lg)
         call diary(sg)
         call diary(kp)

         call glean(thy(cr))
         call glean(thy(lg))
         call glean(thy(sg))
         call glean(thy(dg))
         call glean(thy(kp))

900      if ( error("Exiting check_kpoints") ) continue
         if ( .not.error() ) call stop_timer("check_kpoints: total time")

      end subroutine

      end module check_kpoints_mod

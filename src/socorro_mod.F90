!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module socorro_mod
!doc$ module socorro_mod

!     This socorro module determines if the pre-processing checks
!     or the main electronic structure calculations are performed
!     and executes them accordingly.

      use arg_mod
      use born_oppenheimer_mod
      use check_kpoints_mod
      use check_symmetry_mod
      use config_fh_mod
      use config_sc_mod
      use config_td_mod
      use ehrenfest_mod
      use error_mod
      use kind_mod
      use many_body_theory_mod
      use relax_mod
      use system_mod
      use transition_state_mod

!cod$
      implicit none ; private

      character(line_len) :: mode
      logical :: found
      type(config_fh_obj) :: cfg_fh
      type(config_sc_obj) :: cfg_sc
      type(config_td_obj) :: cfg_td

!doc$
      public :: socorro

!cod$
      interface socorro
         module procedure socorro_
      end interface

      contains

! Public routines

      subroutine socorro_()
!doc$ subroutine socorro()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         ! Create the runtime environment

         call system_start() ; if ( error() ) goto 900

         ! Determine check k-points status

         call arglc("check_kpoints",mode,found) ; if ( .not.found ) mode = "no"
         select case ( trim( mode ) )
         case ( "y" , "yes" )
            call check_kpoints() ; goto 900
         end select

         ! Determine check symmetry status

         call arglc("check_symmetry",mode,found) ; if ( .not.found ) mode = "no"
         select case ( trim( mode ) )
         case ( "y" , "yes" )
            call check_symmetry() ; goto 900
         end select

         ! Determine main calculation type

         call arglc("config_type",mode,found) ; if ( .not.found ) mode = "none"
         select case ( trim( mode ) )
         case ( "fh" , "fixed-hamiltonian" )
            call my(config_fh(),cfg_fh) ; if ( error() ) goto 900
            call diary(cfg_fh)
            call decompose(cfg_fh)
            call glean(thy(cfg_fh))
         case ( "sc" , "self-consistent" )
            call my(config_sc(),cfg_sc)       ; if ( error() ) goto 900
            call diary(cfg_sc)
            call forces(cfg_sc)               ; if ( error() ) goto 100
            call diary_forces(cfg_sc)
            if ( born_oppenheimer_dynamics(cfg_sc) ) then
               if ( error() ) goto 100
               call diary(cfg_sc)
            end if
            if ( many_body_theory(cfg_sc) ) then
               if ( error() ) goto 100
               call diary(cfg_sc)
            end if
            if ( optimize_structure(cfg_sc) ) then
               if ( error() ) goto 100
               call diary(cfg_sc)
            end if
            if ( transition_state(cfg_sc) ) then
               if ( error() ) goto 100
               call diary(cfg_sc)
            end if
            call pressure(cfg_sc)             ; if ( error() ) goto 100
            call diary_pressure(cfg_sc)
            call stress_tensor(cfg_sc)        ; if ( error() ) goto 100
            call diary_stress_tensor(cfg_sc)
            call decompose(cfg_sc)            ; if ( error() ) goto 100
            call write_els_potential(cfg_sc)  ; if ( error() ) goto 100
            call write_restart(cfg_sc)
100         call glean(thy(cfg_sc))
         case ( "td" , "time-dependent" )
            call ehrenfest_dynamics()
         case ( "none" )
            if ( error(.true.,"A config_type was not provided") ) continue
         case default
            if ( error(.true.,"The config_type was not recognized") ) continue
         end select

         ! Destroy the runtime environment

900      if ( error("Exiting socorro") ) continue
         call system_stop()

      end subroutine

      end module socorro_mod

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module socorro_mod
!doc$ module socorro_mod

!     Some info ...
!     ...

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

      logical :: found
      character(line_len) :: mode
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

! public routines

      subroutine socorro_()
!doc$ subroutine socorro()
!       effects:
!       errors:
!       requires:

!cod$
         call system_start() ; if ( error() ) goto 200

         call arglc("check_kpoints",mode,found) ; if ( .not.found ) mode = "no"
         select case ( trim( mode ) )
         case ( "y" , "yes" )
            call check_kpoints() ; goto 200
         end select

         call arglc("check_symmetry",mode,found) ; if ( .not.found ) mode = "no"
         select case ( trim( mode ) )
         case ( "y" , "yes" )
            call check_symmetry() ; goto 200
         end select

         call arglc("config_type",mode,found) ; if ( .not.found ) mode = "self-consistent"
         select case ( trim( mode ) )
         case ( "sc" , "self-consistent" )
            call my(config_sc(),cfg_sc)       ; if ( error() ) goto 200
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
         case ( "fh" , "fixed-hamiltonian" )
            call my(config_fh(),cfg_fh)       ; if ( error() ) goto 200
            call diary(cfg_fh)
            call decompose(cfg_fh)
            call glean(thy(cfg_fh))
         case default
            if ( error(.true.,"ERROR: config_type is not recognized") ) continue
         end select

200      if ( error("Exiting socorro") ) continue
         call system_stop()

      end subroutine

      end module socorro_mod

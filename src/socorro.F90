!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!                                                                                                                                  !
!  Copyright (2020). See the README file in the top-level directory.                                                               !
!  This software is distributed with the GNU General Public License.                                                               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      program socorro

      use arg_mod
      use kind_mod
      use system_mod
      use error_mod
      use config_sc_mod
      use config_td_mod
      use config_fh_mod
      use relax_mod
      use born_oppenheimer_mod
      use many_body_theory_mod
      use transition_state_mod
      use ehrenfest_mod

      implicit none

      logical :: found
      character(line_len) :: mode
      type(config_fh_obj) :: cfg_fh
      type(config_sc_obj) :: cfg_sc
      type(config_td_obj) :: cfg_td

      call system_start() ; if ( error() ) goto 200

      call arglc("config_type",mode,found) ; if ( .not.found ) mode = "self-consistent"

      select case ( trim(mode) )
      case ( "self-consistent" , "sc" )
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
         if ( optimize_crystal(cfg_sc) ) then
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
100      call glean(thy(cfg_sc))
      case ( "time-dependent" , "td" )
         call ehrenfest_dynamics()
      case ( "fixed-hamiltonian" , "fh" )
         call my(config_fh(),cfg_fh)       ; if ( error() ) goto 200
         call diary(cfg_fh)
         call decompose(cfg_fh)
         call glean(thy(cfg_fh))
      case default
         if ( error(.true.,"ERROR: config_type is not recognized") ) continue
      end select

200   if ( error("Exiting socorro") ) continue
      call system_stop() ; stop

      end program

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module atomic_density_mod
!doc$ module atomic_density_mod

      use kind_mod
      use error_mod
      use tagio_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use symmetry_mod
      use atomic_operators_mod
      use atomic_density_ncp_mod
      use atomic_density_paw_mod

!     One datatype is defined here: type(atomic_density_obj).

!cod$
      implicit none
      private

      type, public :: atomic_density_obj
        private
        integer :: type
        type(atomic_density_ncp_obj) :: ad_ncp
        type(atomic_density_paw_obj) :: ad_paw
      end type

!doc$
      public :: atomic_density
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_type
      public :: x_atomic_density_ncp
      public :: x_atomic_density_paw
      public :: trivial
      public :: distance
      public :: symmetrize
      public :: merge_atomic_density
      public :: add_atomic_density
      public :: atomic_hartree_density
      public :: guess_density
      public :: atomic_energy
      public :: atomic_forces
      public :: atomic_pressure
      public :: atomic_stress_tensor
      public :: extract_density
      public :: insert_density
      public :: get_normalization
      public :: write_restart

!cod$
      interface atomic_density
        module procedure constructor_ad
      end interface
      interface update
        module procedure update_ad
      end interface
      interface my
        module procedure my_ad, my_new_ad
      end interface
      interface thy
        module procedure thy_ad
      end interface
      interface glean
        module procedure glean_ad
      end interface
      interface bequeath
        module procedure bequeath_ad
      end interface
      interface assignment(=)
        module procedure assign_ad
      end interface
      interface x_ref
        module procedure ad_ref
      end interface
      interface x_ghost
        module procedure ad_ghost
      end interface
      interface x_type
        module procedure ad_type
      end interface
      interface x_atomic_density_ncp
        module procedure ad_atomic_density_ncp
      end interface
      interface x_atomic_density_paw
        module procedure ad_atomic_density_paw
      end interface
      interface trivial
        module procedure trivial_ad
      end interface
      interface distance
        module procedure distance_ad
      end interface
      interface symmetrize
        module procedure symmetrize_ad
      end interface
      interface merge_atomic_density
        module procedure merge_atomic_density_ad
      end interface
      interface add_atomic_density
        module procedure add_atomic_density_ad
      end interface
      interface atomic_hartree_density
        module procedure atomic_hartree_density_ad
      end interface
      interface guess_density
        module procedure guess_density_ad
      end interface
      interface atomic_energy
        module procedure atomic_energy_ad
      end interface
      interface atomic_forces
        module procedure atomic_forces_ad
      end interface
      interface atomic_pressure
        module procedure atomic_pressure_ad
      end interface
      interface atomic_stress_tensor
        module procedure atomic_stress_tensor_ad
      end interface
      interface extract_density
        module procedure extract_density_ad
      end interface
      interface insert_density
        module procedure insert_density_ad
      end interface
      interface get_normalization
        module procedure get_normalization_ad
      end interface
      interface write_restart
        module procedure write_restart_ad
      end interface

      contains

! public routines

      function constructor_ad(ao,restf,empty) result(ad)
!doc$ function atomic_density(ao,restf,empty) result(ad)
        type(atomic_operators_obj) :: ao
        type(tagio_obj), optional :: restf
        logical, optional :: empty
        type(atomic_density_obj) :: ad
!       effects: Constructs a new ad.

!cod$
        ad%type = x_type(ao)

        select case (ad%type)
        case (NCP)
          if (present(restf)) then
            call my(atomic_density_ncp(x_atomic_operators_ncp(ao),restf),ad%ad_ncp) ; if (error()) goto 100
          else
            call my(atomic_density_ncp(x_atomic_operators_ncp(ao)),ad%ad_ncp) ; if (error()) goto 100
          end if
          call bequeath(thy(ad%ad_ncp))
        case (PAW)
          if (present(restf)) then
            call my(atomic_density_paw(x_atomic_operators_paw(ao),restf),ad%ad_paw) ; if (error()) goto 100
          else
            if (present(empty)) then
              call my(atomic_density_paw(x_atomic_operators_paw(ao),empty=empty),ad%ad_paw) ; if (error()) goto 100
            else
              call my(atomic_density_paw(x_atomic_operators_paw(ao)),ad%ad_paw) ; if (error()) goto 100
            end if
          end if
          call bequeath(thy(ad%ad_paw))
        end select

100     if (error("Exit atomic_density_mod::constructor_ad")) continue

      end function

      subroutine update_ad(ad,ao)
!doc$ subroutine update(ad,ao)
        type(atomic_density_obj) :: ad
        type(atomic_operators_obj) :: ao
!       modifies: ad
!       requires: ad%type = ao%type.
!       effects: Updates ad.

!cod$
        select case (ad%type)
        case (NCP)
          call update(ad%ad_ncp,x_atomic_operators_ncp(ao))
        case (PAW)
          call update(ad%ad_paw,x_atomic_operators_paw(ao))
        end select

        if (error("Exit atomic_density_mod::update_ad")) continue

      end subroutine

      subroutine my_ad(ad)
!doc$ subroutine my(ad)
        type(atomic_density_obj) :: ad

!cod$
       select case (ad%type)
        case (NCP)
          call my(ad%ad_ncp)
        case (PAW)
          call my(ad%ad_paw)
        end select
      end subroutine

      subroutine my_new_ad(adi,ad)
!doc$ subroutine my(adi,ad)
        type(atomic_density_obj) :: adi, ad

!cod$
        ad%type = adi%type
        select case (ad%type)
        case (NCP)
          call my(adi%ad_ncp,ad%ad_ncp)
        case (PAW)
          call my(adi%ad_paw,ad%ad_paw)
        end select
      end subroutine

      function thy_ad(ad) result(ado)
!doc$ function thy(ad) result(ado)
        type(atomic_density_obj) :: ad, ado

!cod$
        ado%type = ad%type
        select case (ad%type)
        case (NCP)
          call my(thy(ad%ad_ncp),ado%ad_ncp)
          call bequeath(thy(ado%ad_ncp))
        case (PAW)
          call my(thy(ad%ad_paw),ado%ad_paw)
          call bequeath(thy(ado%ad_paw))
        end select
      end function

      subroutine glean_ad(ad)
!doc$ subroutine glean(ad)
        type(atomic_density_obj) :: ad

!cod$
        select case (ad%type)
        case (NCP)
          call glean(ad%ad_ncp)
        case (PAW)
          call glean(ad%ad_paw)
        end select
      end subroutine

      subroutine bequeath_ad(ad)
!doc$ subroutine bequeath(ad)
        type(atomic_density_obj) :: ad

!cod$
        select case (ad%type)
        case (NCP)
          call bequeath(ad%ad_ncp)
        case (PAW)
          call bequeath(ad%ad_paw)
        end select
      end subroutine

      subroutine assign_ad(ad,ad2)
!doc$ subroutine assignment(=)(ad,ad2)
        type(atomic_density_obj), intent(inout) :: ad
        type(atomic_density_obj), intent(in) :: ad2
!       requires: ad and ad2 have the same type.

!cod$
        select case (ad%type)
        case (NCP)
          ad%ad_ncp = ad2%ad_ncp
        case (PAW)
          ad%ad_paw = ad2%ad_paw
        end select
      end subroutine

      function ad_ref(ad) result(r)
!doc$ function x_ref(ad) result(r)
        type(atomic_density_obj) :: ad
        integer, dimension(2) :: r
!       effects: Returns ad%ref and ad%o%ref.

!cod$
        select case (ad%type)
        case (NCP)
          r = x_ref(ad%ad_ncp)
        case (PAW)
          r = x_ref(ad%ad_paw)
        end select
      end function

      function ad_ghost(ad) result(g)
!doc$ function x_ghost(ad) result(g)
        type(atomic_density_obj) :: ad
        type(ghost) :: g
!       effects: Returns the ghost of ad.

!cod$
        select case (ad%type)
        case (NCP)
          g = x_ghost(ad%ad_ncp)
        case (PAW)
          g = x_ghost(ad%ad_paw)
        end select
      end function

      function ad_type(ad) result(t)
!doc$ function x_type(ad) result(t)
        type(atomic_density_obj) :: ad
        integer :: t
!       effects: Returns ad%type.

!cod$
        t = ad%type
        call glean(ad)
      end function

      function ad_atomic_density_ncp(ad) result(ad_ncp)
!doc$ function x_atomic_density_ncp(ad) result(ad_ncp)
        type(atomic_density_obj) :: ad
        type(atomic_density_ncp_obj) :: ad_ncp
!       requires: x_type(ad) = NCP
!       effects: Returns ad%ad_ncp.

!cod$
        call my(ad%ad_ncp,ad_ncp)
        call bequeath(thy(ad_ncp))
      end function

      function ad_atomic_density_paw(ad) result(ad_paw)
!doc$ function x_atomic_density_paw(ad) result(ad_paw)
        type(atomic_density_obj) :: ad
        type(atomic_density_paw_obj) :: ad_paw
!       requires: x_type(ad) = PAW
!       effects: Returns ad%ad_paw.

!cod$
        call my(ad%ad_paw,ad_paw)
        call bequeath(thy(ad_paw))
      end function

      function trivial_ad(ad) result(t)
!doc$ function trivial(ad) result(t)
        type(atomic_density_obj) :: ad
        logical :: t
!       effects: Returns the triviality of ad.

!cod$
        select case (ad%type)
        case (NCP)
          t = .true.
        case (PAW)
          t = .false.
        end select
        call glean(ad)
      end function

      function distance_ad(ad1,ad2) result(d)
!doc$ function distance(ad1,ad2) result(d)
        type(atomic_density_obj) :: ad1, ad2
        real(double) :: d
!       requires: ad1 and ad2 have the same types.
!       effects: Returns the norm of the difference between ad1 and ad2.

!cod$
        select case (ad1%type)
        case (NCP)
          d = 0.0_double
          call glean(ad1)
          call glean(ad2)
        case (PAW)
          d = distance(ad1%ad_paw,ad2%ad_paw)
        end select
        if (error("Exit atomic_density_mod::distance_ad")) continue
      end function

      subroutine symmetrize_ad(ad,sg)
!doc$ subroutine symmetrize(ad,sg)
        type(atomic_density_obj) :: ad
        type(space_group_obj) :: sg
!       effects: Symmetrizes ad over the operations in sg.

!cod$
        select case (ad%type)
        case (NCP)
          call glean(ad)
          call glean(sg)
        case (PAW)
          call symmetrize(ad%ad_paw,sg)
        end select
      end subroutine

      subroutine merge_atomic_density_ad(ad)
!doc$ subroutine merge_atomic_density(ad)
        type(atomic_density_obj) :: ad
!       effects: Merges ad contributions from different kgroups.

!cod$
        select case (ad%type)
        case (NCP)
          call merge_atomic_density(ad%ad_ncp)
        case (PAW)
          call merge_atomic_density(ad%ad_paw)
        end select
      end subroutine

      subroutine add_atomic_density_ad(ad,pdots,weights)
!doc$ subroutine add_atomic_density(ad,pdots,weights)
        type(atomic_density_obj) :: ad
        complex(double), dimension(:,:), intent(in) :: pdots
        real(double), dimension(:), intent(in) :: weights
!       modifies: ad
!       effects: Extracts an atomic density from pdots and adds it to ad.

!cod$
        select case (ad%type)
        case (NCP)
          call add_atomic_density(ad%ad_ncp,pdots,weights)
        case (PAW)
          call add_atomic_density(ad%ad_paw,pdots,weights)
        end select
      end subroutine

      function atomic_hartree_density_ad(ad,lay) result(ahd)
!doc$ function atomic_hartree_density(ad,lay) result(ahd)
        type(atomic_density_obj) :: ad
        type(layout_obj) :: lay
        type(grid_obj) :: ahd
!       effects: Returns the atomic contribution to the hartree density.

!cod$
        select case (ad%type)
        case (NCP)
          call my(atomic_hartree_density(ad%ad_ncp,lay),ahd)
        case (PAW)
          call my(atomic_hartree_density(ad%ad_paw,lay),ahd)
        end select
        call bequeath(thy(ahd))
      end function

      function guess_density_ad(ad,lay) result(gd)
!doc$ function guess_density(ad,lay) result(gd)
        type(atomic_density_obj) :: ad
        type(layout_obj) :: lay
        type(grid_obj) :: gd
!       effects: Returns a guess for the valence density constructed from ad.

!cod$
        select case (ad%type)
        case (NCP)
          call my(guess_density(ad%ad_ncp,lay),gd)
        case (PAW)
          call my(guess_density(ad%ad_paw,lay),gd)
        end select
        call bequeath(thy(gd))
      end function

      subroutine atomic_energy_ad(ad,e)
!doc$ subroutine atomic_energy(ad,e)
        type(atomic_density_obj) :: ad
        real(double), intent(out) :: e
!       modifies: e
!       effects: Returns atomic contributions to the energy.
!       errors: Passes errors.

!cod$
        select case (ad%type)
        case (NCP)
          call atomic_energy(ad%ad_ncp,e)
        case (PAW)
          call atomic_energy(ad%ad_paw,e)
        end select
      end subroutine

      subroutine atomic_forces_ad(ad,den,xcp,ahd,ccd,f)
!doc$ subroutine atomic_forces(ad,den,xcp,ahd,ccd,f)
        type(atomic_density_obj) :: ad
        type(grid_obj) :: den, xcp, ahd, ccd
        real(double), dimension(:,:), intent(out) :: f
!       requires: Consistent layouts. f be dimension(3,number-of-atoms).
!       modifies: f
!       effects: Returns unsymmetrized atomic contributions to the forces.

!cod$
        select case (ad%type)
        case (NCP)
          call atomic_forces(ad%ad_ncp,den,xcp,ccd,f)
        case (PAW)
          call atomic_forces(ad%ad_paw,den,xcp,ahd,f)
        end select
      end subroutine

      subroutine atomic_pressure_ad(ad,den,xcp,p)
!doc$ subroutine atomic_pressure(ad,den,xcp,p)
        type(atomic_density_obj) :: ad
        type(grid_obj) :: den, xcp
        real(double), intent(out) :: p
!       requires: Consistent layouts.
!       modifies: p
!       effects: Returns atomic contributions to the pressure.
!       errors: If ad type is PAW.

!cod$
        select case (ad%type)
        case (NCP)
          call atomic_pressure(ad%ad_ncp,den,xcp,p)
        case (PAW)
          call glean(den)
          call glean(xcp)
          if (error(.true.,"ERROR: atomic_pressure_ad is not yet implemented for type = PAW")) continue
        end select
        if (error("Exit atomic_density_mod::atomic_pressure_ad")) continue
      end subroutine

      subroutine atomic_stress_tensor_ad(ad,den,xcp,s)
!doc$ subroutine atomic_stress_tensor(ad,den,xcp,s)
        type(atomic_density_obj) :: ad
        type(grid_obj) :: den, xcp
        real(double), dimension(:,:), intent(out) :: s
!       requires: Consistent layouts. s be dimension(3,3).
!       modifies: s
!       effects: Returns atomic contributions to the stress_tensor.
!       errors: If ad type is PAW.

!cod$
        select case (ad%type)
        case (NCP)
          call atomic_stress_tensor(ad%ad_ncp,den,xcp,s)
        case (PAW)
          call glean(den)
          call glean(xcp)
          if (error(.true.,"ERROR: atomic_stress_tensor_ad is not yet implemented for type = PAW")) continue
        end select
        if (error("Exit atomic_density_mod::atomic_stress_tensor_ad")) continue
      end subroutine

      subroutine extract_density_ad(ad,c1d)
!doc$ subroutine extract_density(ad,c1d)
        type(atomic_density_obj) :: ad
        complex(double), dimension(:), pointer :: c1d
!       requires: c1d be nullified or associated.
!       modifies: c1d
!       effects: Returns ad%o%wij in c1d.

!cod$
        select case (ad%type)
        case (NCP)
          if (associated( c1d )) deallocate( c1d ) ; allocate( c1d(0) )
          call glean(ad)
        case (PAW)
          call extract_density(ad%ad_paw,c1d)
        end select
        if (error("Exit atomic_density_mod::extract_density_ad")) continue
      end subroutine

      subroutine insert_density_ad(c1d,ad)
!doc$ subroutine insert_density(c1d,ad)
        complex(double), dimension(:), pointer :: c1d
        type(atomic_density_obj) :: ad
!       modifies: ad
!       effects: Inserts c1d into ad%o%wij.

!cod$
        select case (ad%type)
        case (NCP)
          call glean(ad)
        case (PAW)
          call insert_density(c1d,ad%ad_paw)
        end select
        if (error("Exit atomic_density_mod::insert_density_ad")) continue
      end subroutine

      subroutine get_normalization_ad(ad,ne)
!doc$ subroutine get_normalization(ad,ne)
        type(atomic_density_obj) :: ad
        real(double) :: ne
!       effects: Returns the number of electrons in ad.

!cod$
        select case (ad%type)
        case (NCP)
          ne = 0.0_double
          call glean(ad)
        case (PAW)
          call get_normalization(ad%ad_paw,ne)
        end select
        if (error("Exit atomic_density_mod::get_normalization_ad")) continue
      end subroutine

      subroutine write_restart_ad(ad,nrestf)
!doc$ subroutine write_restart(ad,nrestf)
        type(atomic_density_obj) :: ad
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ad restart information to nrestf.

!cod$
        select case (ad%type)
        case (NCP)
          call write_restart(ad%ad_ncp,nrestf)
        case (PAW)
          call write_restart(ad%ad_paw,nrestf)
        end select
      end subroutine

      end module

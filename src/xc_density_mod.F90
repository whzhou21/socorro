!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module xc_density_mod
!doc$ module xc_density_mod

      use kind_mod
      use error_mod
      use ghost_mod
      use layout_mod
      use symmetry_mod
      use grid_mod
      use xc_type_mod
      use xc_density_native_mod
      use xc_density_libxc_mod

!     One datatype is defined here: type(xc_density_obj).

!     xc_density_mod is a wrapper for density-dependent native and libxc exchange-correlation
!     routines that operate on grid data.

!cod$
      implicit none
      private

      type, public :: xc_density_obj
        private
        integer :: source                                    ! Source of the (density-dependent) functionals.
        type(xc_density_native_obj) :: xcd_native
        type(xc_density_libxc_obj)  :: xcd_libxc
      end type

!doc$
      public :: xc_density
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_xc_type
      public :: x_layout
      public :: x_space_group
      public :: uses_gradient
      public :: uses_laplacian
      public :: xc_energy
      public :: xc_potential
      public :: xc_energy_and_potential
      public :: xc_grid_pressure
      public :: xc_grid_stress_tensor

!cod$
      interface xc_density
        module procedure constructor_xcd
      end interface
      interface update
        module procedure update_xcd
      end interface
      interface my
        module procedure my_xcd, my_new_xcd
      end interface
      interface thy
        module procedure thy_xcd
      end interface
      interface glean
        module procedure glean_xcd
      end interface
      interface bequeath
        module procedure bequeath_xcd
      end interface
      interface assignment(=)
        module procedure assign_xcd
      end interface
      interface x_ref
        module procedure xcd_ref
      end interface
      interface x_ghost
        module procedure xcd_ghost
      end interface
      interface x_xc_type
        module procedure xcd_xc_type
      end interface
      interface x_layout
        module procedure xcd_layout
      end interface
      interface x_space_group
        module procedure xcd_space_group
      end interface
      interface uses_gradient
        module procedure q_uses_gradient
      end interface
      interface uses_laplacian
        module procedure q_uses_laplacian
      end interface
      interface xc_energy
        module procedure xc_energy_xcd
      end interface
      interface xc_potential
        module procedure xc_potential_xcd
      end interface
      interface xc_energy_and_potential
        module procedure xc_energy_and_potential_xcd
      end interface
      interface xc_grid_pressure
        module procedure xc_grid_pressure_xcd
      end interface
      interface xc_grid_stress_tensor
        module procedure xc_grid_stress_tensor_xcd
      end interface

      contains

! public routines

      function constructor_xcd(xct,lay,sg) result(xcd)
!doc$ function xc_density(xct,lay,sg) result(xcd)
        type(xc_type_obj) :: xct
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_density_obj) :: xcd
!       effects: Constructs a new xcd.
!       errors: Passes errors

!cod$ 
        call my(xct)
        call my(lay)
        call my(sg)

        xcd%source = x_ddf_source(xct)
        select case(xcd%source)
        case (DDFS_NATIVE)
          call my(xc_density_native(xct,lay,sg),xcd%xcd_native)
          call bequeath(thy(xcd%xcd_native))
        case (DDFS_LIBXC)
          call my(xc_density_libxc(xct,lay,sg),xcd%xcd_libxc)
          call bequeath(thy(xcd%xcd_libxc))
        end select

100     call glean(thy(xct))
        call glean(thy(lay))
        call glean(thy(sg))

        if (error("Exit xc_density_mod::constructor_xcd")) continue

      end function 

      subroutine update_xcd(xcd,lay,sg)
!doc$ subroutine update(xcd,lay,sg)
        type(xc_density_obj) :: xcd
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
!       effects: Updates xcd.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call update(xcd%xcd_native,lay,sg)
        case (DDFS_LIBXC)
          call update(xcd%xcd_libxc,lay,sg)
        end select
        if (error("Exit xc_density_mod::update_xcd")) continue
      end subroutine

      subroutine my_xcd(xcd)
!doc$ subroutine my(xcd)
        type(xc_density_obj) :: xcd

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(xcd%xcd_native)
        case (DDFS_LIBXC)
          call my(xcd%xcd_libxc)
        end select
      end subroutine

      subroutine my_new_xcd(xcdi,xcd)
!doc$ subroutine my(xcdi,xcd)
        type(xc_density_obj) :: xcdi, xcd

!cod$
        xcd%source = xcdi%source
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(xcdi%xcd_native,xcd%xcd_native)
        case (DDFS_LIBXC)
          call my(xcdi%xcd_libxc,xcd%xcd_libxc)
        end select
      end subroutine

      function thy_xcd(xcd) result(xcdo)
!doc$ function thy(xcd) result(xcdo)
        type(xc_density_obj) :: xcd, xcdo

!cod$
        xcdo%source = xcd%source
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(thy(xcd%xcd_native),xcdo%xcd_native)
          call bequeath(thy(xcdo%xcd_native))
        case (DDFS_LIBXC)
          call my(thy(xcd%xcd_libxc),xcdo%xcd_libxc)
          call bequeath(thy(xcdo%xcd_libxc))
        end select
      end function

      subroutine glean_xcd(xcd)
!doc$ subroutine glean(xcd)
        type(xc_density_obj) :: xcd

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call glean(xcd%xcd_native)
        case (DDFS_LIBXC)
          call glean(xcd%xcd_libxc)
        end select
      end subroutine

      subroutine bequeath_xcd(xcd)
!doc$ subroutine bequeath(xcd)
        type(xc_density_obj) :: xcd

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call bequeath(xcd%xcd_native)
        case (DDFS_LIBXC)
          call bequeath(xcd%xcd_libxc)
        end select
      end subroutine
 
      subroutine assign_xcd(xcd,xcd2)
!doc$ subroutine assignment(=)(xcd,xcd2)
        type(xc_density_obj), intent(inout) :: xcd
        type(xc_density_obj), intent(in) :: xcd2

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          xcd%xcd_native = xcd2%xcd_native
        case (DDFS_LIBXC)
          xcd%xcd_libxc = xcd2%xcd_libxc
        end select
      end subroutine

      function xcd_ref(xcd) result(r)
!doc$ function x_ref(xcd) result(r)
        type(xc_density_obj) :: xcd
        integer, dimension(2) :: r
!       effects: Returns %ref and %o%ref of the underlying data structure.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          r = x_ref(xcd%xcd_native)
        case (DDFS_LIBXC)
          r = x_ref(xcd%xcd_libxc)
        end select
      end function
 
      function xcd_ghost(xcd) result(g)
!doc$ function x_ghost(xcd) result(g)
        type(xc_density_obj) :: xcd
        type(ghost) :: g
!       effects: Returns the ghost of the underlying data structure.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          g = x_ghost(xcd%xcd_native)
        case (DDFS_LIBXC)
          g = x_ghost(xcd%xcd_libxc)
        end select
      end function

      function xcd_xc_type(xcd) result(xct)
!doc$ function x_xc_type(xcd) result(xct)
        type(xc_density_obj) :: xcd
        type(xc_type_obj) :: xct
!       effects: Returns the xct of xcd.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(x_xc_type(xcd%xcd_native),xct)
        case (DDFS_LIBXC)
          call my(x_xc_type(xcd%xcd_libxc),xct)
        end select
        call bequeath(thy(xct))
      end function

      function xcd_layout(xcd) result(lay)
!doc$ function x_layout(xcd) result(lay)
        type(xc_density_obj) :: xcd
        type(layout_obj) :: lay
!       effects: Returns the lay of xcd.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(x_layout(xcd%xcd_native),lay)
        case (DDFS_LIBXC)
          call my(x_layout(xcd%xcd_libxc),lay)
        end select
        call bequeath(thy(lay))
      end function

      function xcd_space_group(xcd) result(sg)
!doc$ function x_space_group(xcd) result(sg)
        type(xc_density_obj) :: xcd
        type(space_group_obj) :: sg
!       effects: Returns the sg of xcd.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(x_space_group(xcd%xcd_native),sg)
        case (DDFS_LIBXC)
          call my(x_space_group(xcd%xcd_libxc),sg)
        end select
        call bequeath(thy(sg))
      end function

      function q_uses_gradient(xcd) result(ug)
!doc$ function uses_gradient(xcd) result(ug)
        type(xc_density_obj) :: xcd
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        select case (xcd%source)
        case (DDFS_NATIVE)
          ug = uses_gradient(xcd%xcd_native)
        case (DDFS_LIBXC)
          ug = uses_gradient(xcd%xcd_libxc)
        end select
      end function

      function q_uses_laplacian(xcd) result(ul)
!doc$ function uses_laplacian(xcd) result(ul)
        type(xc_density_obj) :: xcd
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        select case (xcd%source)
        case (DDFS_NATIVE)
          ul = uses_laplacian(xcd%xcd_native)
        case (DDFS_LIBXC)
          ul = uses_laplacian(xcd%xcd_libxc)
        end select
      end function

      function xc_energy_xcd(xcd,n_total,n_valence) result(fxc)
!doc$ function xc_energy(xcd,n_total,n_valence) result(fxc)
        type(xc_density_obj) :: xcd
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        real(double) :: fxc
!       requires: n_total and n_valence data be symmetrized.
!       effects: Returns the exchange-correlation energy.
!       errors: Passes errors.

!cod$
        call my(n_total)
        call my(n_valence)

        select case (xcd%source)
        case (DDFS_NATIVE)
          fxc = xc_energy(xcd%xcd_native,n_total)
        case (DDFS_LIBXC)
          fxc = xc_energy(xcd%xcd_libxc,n_total,n_valence)
        end select 
        call glean(thy(n_total))
        call glean(thy(n_valence))
        if (error("Exit xc_density_mod::xc_energy_xcd")) continue
      end function 

      function xc_potential_xcd(xcd,n_total,n_valence) result(vxc_g)
!doc$ function xc_potential(xcd,n_total,n_valence) result(vxc_g)
        type(xc_density_obj) :: xcd
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        type(grid_obj) :: vxc_g
!       requires: n_total and n_valence data be symmetrized with respect to xcd%o%sg.
!       effects: Returns the symmetrized/filtered exchange-correlation potential.
!       errors: Passes errors.

!cod$ 
        call my(n_total) 
        call my(n_valence) 
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(xc_potential(xcd%xcd_native,n_total),vxc_g) ; if (error()) goto 100
        case (DDFS_LIBXC)
          call my(xc_potential(xcd%xcd_libxc,n_total,n_valence),vxc_g) ; if (error()) goto 100
        end select 
        call bequeath(thy(vxc_g))
        call glean(thy(n_valence))
        call glean(thy(n_total))
100     if (error("Exit xc_density_mod::cd_potential_xcd")) continue
      end function

      function xc_energy_and_potential_xcd(xcd,n_total,n_valence,fxc) result(vxc_g)
!doc$ function xc_energy_and_potential(xcd,n_total,n_valence,fxc) result(vxc_g)
        type(xc_density_obj) :: xcd
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        real(double), intent(out) :: fxc
        type(grid_obj) :: vxc_g
!       requires: n_total and n_valence data be symmetrized. Functional be semilocal.
!       effects: Returns the symmetrized/filtered exchange-correlation potential and energy.
!       errors: Passes errors.

!cod$
        call my(n_total)
        call my(n_valence)
        select case (xcd%source)
        case (DDFS_NATIVE)
          call my(xc_energy_and_potential(xcd%xcd_native,n_total,fxc),vxc_g) ; if (error()) goto 100
        case (DDFS_LIBXC)
          call my(xc_energy_and_potential(xcd%xcd_libxc,n_total,n_valence,fxc),vxc_g) ; if (error()) goto 100
        end select 
        call glean(thy(n_valence))
        call glean(thy(n_total))
        call bequeath(thy(vxc_g))
100     if (error("Exit xc_density_mod::xc_energy_and_potential_xcd")) continue
      end function

      subroutine xc_grid_pressure_xcd(xcd,n,p)
!doc$ subroutine xc_grid_pressure(xcd,n,p)
        type(xc_density_obj) :: xcd
        type(grid_obj) :: n
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the exchange-correlation grid potential.

!cod$
        call my(n)

        select case (xcd%source)
        case (DDFS_NATIVE)
          call xc_grid_pressure(xcd%xcd_native,n,p)
        case (DDFS_LIBXC)
          if (error(.true.,"ERROR: grid_pressure routine is not currently available")) continue
!          call xc_grid_pressure(xcd%xcd_libxc,n,p)
        end select 
        call glean(thy(n))
        if (error("Exit xc_density_mod::xc_grid_pressure_xcd")) continue
      end subroutine

      subroutine xc_grid_stress_tensor_xcd(xcd,n_g,s)
!doc$ subroutine xc_grid_stress_tensor(xcd,n_g,s)
        type(xc_density_obj) :: xcd
        type(grid_obj) :: n_g
        real(double), dimension(:,:), intent(out) :: s
!       modifies: s
!       effects: Returns stress tensor contributions due to the exchange-correlation grid potential.

!cod$
        call my(n_g)
        select case (xcd%source)
        case (DDFS_NATIVE)
          call xc_grid_stress_tensor(xcd%xcd_native,n_g,s)
        case (DDFS_LIBXC)
          if (error(.true.,"ERROR: grid_stress_tensor routine is not currently available")) continue
!          call xc_grid_stress_tensor(xcd%xcd_libxc,n_g,s)
        end select 
        call glean(thy(n_g))
        if (error("Exit xc_density_mod::xc_grid_stress_tensor_xcd")) continue
      end subroutine

      end module

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module xc_mod
!doc$ module xc_mod

      use kind_mod
      use error_mod
      use ghost_mod
      use layout_mod
      use symmetry_mod
      use grid_mod
      use xc_type_mod
      use xc_density_mod
      use xc_orbital_mod
      use electrons_sc_mod
      use multivector_mod

!     One datatype is defined here: type(xc_obj).  

!     xc_mod is a wrapper for exchange-correlation routines that operate on grid data.

!cod$
      implicit none
      private

      type, public :: xc_obj
        private
        integer :: dependence                                ! functional dependence
        type(xc_density_obj) :: xc_density
        type(xc_orbital_obj) :: xc_orbital
      end type

!doc$
      public :: xc
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
      public :: x_functional_dependence
      public :: x_hybrid_mixing
      public :: uses_gradient
      public :: uses_laplacian
      public :: xc_energy
      public :: xc_potential
      public :: xc_energy_and_potential
      public :: xc_grid_pressure
      public :: xc_grid_stress_tensor
      public :: xc_derivative
      public :: xc_energy_and_derivative

!cod$
      interface xc
        module procedure constructor_xc
      end interface
      interface update
        module procedure update_xc
      end interface
      interface my
        module procedure my_xc, my_new_xc
      end interface
      interface thy
        module procedure thy_xc
      end interface
      interface glean
        module procedure glean_xc
      end interface
      interface bequeath
        module procedure bequeath_xc
      end interface
      interface assignment(=)
        module procedure assign_xc
      end interface
      interface x_ref
        module procedure xc_ref
      end interface
      interface x_ghost
        module procedure xc_ghost
      end interface
      interface x_xc_type
        module procedure xc_xc_type
      end interface
      interface x_layout
        module procedure xc_layout
      end interface
      interface x_space_group
        module procedure xc_space_group
      end interface
      interface x_functional_dependence
        module procedure xc_functional_dependence
      end interface
      interface x_hybrid_mixing
        module procedure xc_hybrid_mixing
      end interface
      interface uses_gradient
        module procedure q_uses_gradient
      end interface
      interface uses_laplacian
        module procedure q_uses_laplacian
      end interface
      interface xc_energy
        module procedure xcd_energy_xc, xco_energy_xc
      end interface
      interface xc_potential
        module procedure xcd_potential_xc
      end interface
      interface xc_energy_and_potential
        module procedure xcd_energy_and_potential_xc
      end interface
      interface xc_grid_pressure
        module procedure xcd_grid_pressure_xc
      end interface
      interface xc_grid_stress_tensor
        module procedure xcd_grid_stress_tensor_xc
      end interface
      interface xc_derivative
        module procedure xco_derivative_xc
      end interface
      interface xc_energy_and_derivative
        module procedure xco_energy_and_derivative_xc
      end interface

      contains

! public routines

      function constructor_xc(xc_type,lay,sg) result(xc)
!doc$ function xc(xc_type,lay,sg) result(xc)
        type(xc_type_obj) :: xc_type
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_obj) :: xc
!       effects: Constructs a new xc.
!       errors: Passes errors

!cod$ 
        call my(xc_type)
        call my(lay)
        call my(sg)

        xc%dependence = x_functional_dependence(xc_type)
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(xc_density(xc_type,lay,sg),xc%xc_density)
          call bequeath(thy(xc%xc_density))
        case (FD_ORBITAL)
          call my(xc_orbital(xc_type,lay,sg),xc%xc_orbital)
          call bequeath(thy(xc%xc_orbital))
        case (FD_HYBRID)
          call my(xc_density(xc_type,lay,sg),xc%xc_density)
          call bequeath(thy(xc%xc_density))
          call my(xc_orbital(xc_type,lay,sg),xc%xc_orbital)
          call bequeath(thy(xc%xc_orbital))
        end select

100     call glean(thy(xc_type))
        call glean(thy(lay))
        call glean(thy(sg))

        if (error("Exit xc_mod::constructor_xc")) continue

      end function 

      subroutine update_xc(xc,lay,sg)
!doc$ subroutine update(xc,lay,sg)
        type(xc_obj) :: xc
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
!       effects: Updates xc.

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call update(xc%xc_density,lay,sg)
        case (FD_ORBITAL)
          call update(xc%xc_orbital,lay,sg)
        case (FD_HYBRID)
          call update(xc%xc_density,lay,sg)
          call update(xc%xc_orbital,lay,sg)
        end select
        if (error("Exit xc_mod::update_xc")) continue
      end subroutine

      subroutine my_xc(xc)
!doc$ subroutine my(xc)
        type(xc_obj) :: xc

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(xc%xc_density)
        case (FD_ORBITAL)
          call my(xc%xc_orbital)
        case (FD_HYBRID)
          call my(xc%xc_density)
          call my(xc%xc_orbital)
        end select
      end subroutine

      subroutine my_new_xc(xci,xc)
!doc$ subroutine my(xci,xc)
        type(xc_obj) :: xci, xc

!cod$
        xc%dependence = xci%dependence
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(xci%xc_density,xc%xc_density)
        case (FD_ORBITAL)
          call my(xci%xc_orbital,xc%xc_orbital)
        case (FD_HYBRID)
          call my(xci%xc_density,xc%xc_density)
          call my(xci%xc_orbital,xc%xc_orbital)
        end select
      end subroutine

      function thy_xc(xc) result(xco)
!doc$ function thy(xc) result(xco)
        type(xc_obj) :: xc, xco

!cod$
        xco%dependence = xc%dependence
        select case (xc%dependence)
        case (FD_DENSITY)
          call my(thy(xc%xc_density),xco%xc_density)
          call bequeath(thy(xco%xc_density))
        case (FD_ORBITAL)
          call my(thy(xc%xc_orbital),xco%xc_orbital)
          call bequeath(thy(xco%xc_orbital))
        case (FD_HYBRID)
          call my(thy(xc%xc_density),xco%xc_density)
          call bequeath(thy(xco%xc_density))
          call my(thy(xc%xc_orbital),xco%xc_orbital)
          call bequeath(thy(xco%xc_orbital))
        end select
      end function

      subroutine glean_xc(xc)
!doc$ subroutine glean(xc)
        type(xc_obj) :: xc

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call glean(xc%xc_density)
        case (FD_ORBITAL)
          call glean(xc%xc_orbital)
        case (FD_HYBRID)
          call glean(xc%xc_density)
          call glean(xc%xc_orbital)
        end select
      end subroutine

      subroutine bequeath_xc(xc)
!doc$ subroutine bequeath(xc)
        type(xc_obj) :: xc

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call bequeath(xc%xc_density)
        case (FD_ORBITAL)
          call bequeath(xc%xc_orbital)
        case (FD_HYBRID)
          call bequeath(xc%xc_density)
          call bequeath(xc%xc_orbital)
        end select
      end subroutine
 
      subroutine assign_xc(xc,xc2)
!doc$ subroutine assignment(=)(xc,xc2)
        type(xc_obj), intent(inout) :: xc
        type(xc_obj), intent(in) :: xc2

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          xc%xc_density = xc2%xc_density
        case (FD_ORBITAL)
          xc%xc_orbital = xc2%xc_orbital
        case (FD_HYBRID)
          xc%xc_density = xc2%xc_density
          xc%xc_orbital = xc2%xc_orbital
        end select
      end subroutine

      function xc_ref(xc) result(r)
!doc$ function x_ref(xc) result(r)
        type(xc_obj) :: xc
        integer, dimension(2) :: r
!       effects: Returns %ref and %o%ref for the underlying data structure
!                (the density-dependent data structure when dependence = FD_HYBRID).

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          r = x_ref(xc%xc_density)
        case (FD_ORBITAL)
          r = x_ref(xc%xc_orbital)
        case (FD_HYBRID)
          r = x_ref(xc%xc_density)
        end select
      end function
 
      function xc_ghost(xc) result(g)
!doc$ function x_ghost(xc) result(g)
        type(xc_obj) :: xc
        type(ghost) :: g
!       effects: Returns the ghost of the underlying data structure (the 
!                density-dependent data structure when dependence = FD_HYBRID).

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          g = x_ghost(xc%xc_density)
        case (FD_ORBITAL)
          g = x_ghost(xc%xc_orbital)
        case (FD_HYBRID)
          g = x_ghost(xc%xc_density)
        end select
      end function

      function xc_xc_type(xc) result(xc_type)
!doc$ function x_xc_type(xc) result(xc_type)
        type(xc_obj) :: xc
        type(xc_type_obj) :: xc_type
!       effects: Returns the xc_type of xc.

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(x_xc_type(xc%xc_density),xc_type)
        case (FD_ORBITAL)
          call my(x_xc_type(xc%xc_orbital),xc_type)
        case (FD_HYBRID)
          call my(x_xc_type(xc%xc_density),xc_type)
        end select
        call bequeath(thy(xc_type))
      end function

      function xc_layout(xc) result(lay)
!doc$ function x_layout(xc) result(lay)
        type(xc_obj) :: xc
        type(layout_obj) :: lay
!       effects: Returns the lay of xc.

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(x_layout(xc%xc_density),lay)
        case (FD_ORBITAL)
          call my(x_layout(xc%xc_orbital),lay)
        case (FD_HYBRID)
          call my(x_layout(xc%xc_density),lay)
        end select
        call bequeath(thy(lay))
      end function

      function xc_space_group(xc) result(sg)
!doc$ function x_space_group(xc) result(sg)
        type(xc_obj) :: xc
        type(space_group_obj) :: sg
!       effects: Returns the sg of xc.

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(x_space_group(xc%xc_density),sg)
        case (FD_ORBITAL)
          call my(x_space_group(xc%xc_orbital),sg)
        case (FD_HYBRID)
          call my(x_space_group(xc%xc_density),sg)
        end select
        call bequeath(thy(sg))
      end function

      function xc_functional_dependence(xc) result(fd)
!doc$ function x_functional_dependence(xc) result(fd)
        type(xc_obj) :: xc
        integer :: fd
!       effects: Returns the functional dependence of xc.

!cod$
        fd = xc%dependence
      end function

      function xc_hybrid_mixing(xc) result(hm)
!doc$ function x_hybrid_mixing(xc) result(hm)
        type(xc_obj) :: xc
        real(double) :: hm
!       effects: Returns the hybrid mixing parameter of xc.

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          if (error(.true.,"ERROR: accessor is only valid with hybrid functionals")) continue
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: accessor is only valid with hybrid functionals")) continue
        case (FD_HYBRID)
          hm = x_hybrid_mixing(x_xc_type(xc))
        end select
      end function

      function q_uses_gradient(xc) result(ug)
!doc$ function uses_gradient(xc) result(ug)
        type(xc_obj) :: xc
        logical :: ug
!       effects: Returns .true. iff the density-dependent functional uses gradients of the density.

!cod$ 
        select case(xc%dependence)
        case (FD_DENSITY)
          ug = uses_gradient(xc%xc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: query is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          ug = uses_gradient(xc%xc_density)
        end select
        if (error("Exit xc_mod::q_uses_gradient")) continue
      end function

      function q_uses_laplacian(xc) result(ul)
!doc$ function uses_laplacian(xc) result(ul)
        type(xc_obj) :: xc
        logical :: ul
!       effects: Returns .true. iff the density-dependent functional uses laplacians of the density.

!cod$
        select case(xc%dependence)
        case (FD_DENSITY)
          ul = uses_laplacian(xc%xc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: query is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          ul = uses_laplacian(xc%xc_density)
        end select
        if (error("Exit xc_mod::q_uses_laplacian")) continue
      end function

      function xcd_energy_xc(xc,n_total,n_valence) result(fxc)
!doc$ function xc_energy(xc,n_total,n_valence) result(fxc)
        type(xc_obj) :: xc
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        real(double) :: fxc
!       requires: n_total and n_valence data be symmetrized.
!       effects: Returns the exchange-correlation energy.
!       errors: Passes errors.

!cod$
        call my(n_total)
        call my(n_valence)
        select case(xc%dependence)
        case (FD_DENSITY)
          fxc = xc_energy(xc%xc_density,n_total,n_valence)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          fxc = xc_energy(xc%xc_density,n_total,n_valence)
        end select
        call glean(thy(n_total))
        call glean(thy(n_valence))
        if (error("Exit xc_mod::xcd_energy_xc")) continue
      end function 

      function xcd_potential_xc(xc,n_total,n_valence) result(vxc_g)
!doc$ function xc_potential(xc,n_total,n_valence) result(vxc_g)
        type(xc_obj) :: xc
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        type(grid_obj) :: vxc_g
!       requires: n_total and n_valence data be symmetrized with respect to xcd%o%sg.
!       effects: Returns the symmetrized/filtered exchange-correlation potential.
!       errors: Passes errors.

!cod$ 
        call my(n_total) 
        call my(n_valence) 
        select case(xc%dependence)
        case (FD_DENSITY)
          call my(xc_potential(xc%xc_density,n_total,n_valence),vxc_g) ; if (error()) goto 100
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) goto 100
        case (FD_HYBRID)
          call my(xc_potential(xc%xc_density,n_total,n_valence),vxc_g) ; if (error()) goto 100
        end select
        call bequeath(thy(vxc_g))
100     call glean(thy(n_total))
        call glean(thy(n_valence))
        if (error("Exit xc_mod::xcd_potential_xc")) continue
      end function

      function xcd_energy_and_potential_xc(xc,n_total,n_valence,fxc) result(vxc_g)
!doc$ function xc_energy_and_potential(xc,n_total,n_valence,fxc) result(vxc_g)
        type(xc_obj) :: xc
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

        select case(xc%dependence)
        case (FD_DENSITY)
          call my(xc_energy_and_potential(xc%xc_density,n_total,n_valence,fxc),vxc_g) ; if (error()) goto 100
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) goto 100
        case (FD_HYBRID)
          call my(xc_energy_and_potential(xc%xc_density,n_total,n_valence,fxc),vxc_g) ; if (error()) goto 100
        end select
        call bequeath(thy(vxc_g))
100     call glean(thy(n_total))
        call glean(thy(n_valence))
        if (error("Exit xc_mod::xcd_energy_and_potential_xc")) continue
      end function

      subroutine xcd_grid_pressure_xc(xc,n_g,p)
!doc$ subroutine xc_grid_pressure(xc,n_g,p)
        type(xc_obj) :: xc
        type(grid_obj) :: n_g
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the exchange-correlation grid potential.

!cod$
        call my(n_g)
        select case(xc%dependence)
        case (FD_DENSITY)
          call xc_grid_pressure(xc%xc_density,n_g,p)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: routine is not valid with hybrid functionals")) continue
        end select
        call glean(thy(n_g))
        if (error("Exit xc_mod::xcd_grid_pressure_xc")) continue
      end subroutine

      subroutine xcd_grid_stress_tensor_xc(xc,n_g,s)
!doc$ subroutine xc_grid_stress_tensor(xc,n_g,s)
        type(xc_obj) :: xc
        type(grid_obj) :: n_g
        real(double), dimension(:,:), intent(out) :: s
!       modifies: s
!       effects: Returns stress tensor contributions due to the exchange-correlation grid potential.

!cod$
        call my(n_g)
        select case(xc%dependence)
        case (FD_DENSITY)
          call xc_grid_stress_tensor(xc%xc_density,n_g,s)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: routine is not valid with hybrid functionals")) continue
        end select
        call glean(thy(n_g))
        if (error("Exit xc_mod::xcd_grid_stress_tensor_xc")) continue
      end subroutine

      function xco_energy_xc(xc,el) result (e)
!doc$ function xc_energy(xc,el) result(e)
        type(xc_obj) :: xc
        type(electrons_sc_obj) :: el
        real(double) :: e
!       effects: Returns the exx energy of el.

!cod$
        call my(el)
        select case(xc%dependence)
        case (FD_DENSITY)
          if (error(.true.,"ERROR: routine is only valid with orbital-dependent functionals")) continue
        case (FD_ORBITAL)
          call exx_energy(xc%xc_orbital,el,e)
        case (FD_HYBRID)
          call exx_energy(xc%xc_orbital,el,e)
        end select
        call glean(thy(el))
        if (error("Exit xc_mod::xco_energy_xc")) continue
      end function

      subroutine xco_derivative_xc(xc,el,ik1,mvo)
!doc$ subroutine xc_derivative(xc,el,ik1,mvo)
        type(xc_obj) :: xc
        type(electrons_sc_obj) :: el
        integer :: ik1
        type(multivector_obj) :: mvo
!       effects: Augments mvo with ik1 contributions to the exx functional derivative.

!cod$
        call my(el)
        call my(mvo)
        select case(xc%dependence)
        case (FD_DENSITY)
          if (error(.true.,"ERROR: routine is only valid with orbital-dependent functionals")) continue
        case (FD_ORBITAL)
          call exx_derivative(xc%xc_orbital,el,ik1,mvo)
        case (FD_HYBRID)
          call exx_derivative(xc%xc_orbital,el,ik1,mvo)
        end select
        call glean(thy(el))
        call glean(thy(mvo))
        if (error("Exit xc_mod::xco_derivative_xc")) continue
      end subroutine

      subroutine xco_energy_and_derivative_xc(xc,el,ik1,e,mvo)
!doc$ subroutine xc_energy_and_derivative(xc,el,ik1,e,mvo)
        type(xc_obj) :: xc
        type(electrons_sc_obj) :: el
        integer :: ik1
        real(double) :: e
        type(multivector_obj) :: mvo
!       effects: Augments e and mvo with ik1 contributions to the exx functional and its derivative.

!cod$
        call my(el)
        call my(mvo)
        select case(xc%dependence)
        case (FD_DENSITY)
          if (error(.true.,"ERROR: routine is only valid with orbital-dependent functionals")) continue
        case (FD_ORBITAL)
          call exx_energy_and_derivative(xc%xc_orbital,el,ik1,e,mvo)
        case (FD_HYBRID)
          call exx_energy_and_derivative(xc%xc_orbital,el,ik1,e,mvo)
        end select
        call glean(thy(el))
        call glean(thy(mvo))
        if (error("Exit xc_mod::xco_energy_and_derivative_xc")) continue
      end subroutine

      end module

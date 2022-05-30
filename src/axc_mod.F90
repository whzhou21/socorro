! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module axc_mod
!doc$ module axc_mod

      use kind_mod
      use error_mod
      use ghost_mod
      use xc_type_mod
      use axc_density_mod

!     One datatype is defined here: type(axc_obj).  

!     axc_mod is a wrapper for exchange-correlation routines that operate on atomic data.

!cod$
      implicit none
      private

      type, public :: axc_obj
        private
        integer :: dependence                                ! functional dependence
        type(axc_density_obj) :: axc_density
      end type

!doc$
      public :: axc
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_xc_type
      public :: x_functional_dependence
      public :: uses_gradient
      public :: uses_laplacian
      public :: xc_energy_density
      public :: xc_derivatives

!cod$
      interface axc
        module procedure constructor_axc
      end interface
      interface update
        module procedure update_axc
      end interface
      interface my
        module procedure my_axc, my_new_axc
      end interface
      interface thy
        module procedure thy_axc
      end interface
      interface glean
        module procedure glean_axc
      end interface
      interface bequeath
        module procedure bequeath_axc
      end interface
      interface assignment(=)
        module procedure assign_axc
      end interface
      interface x_ref
        module procedure axc_ref
      end interface
      interface x_ghost
        module procedure axc_ghost
      end interface
      interface x_xc_type
        module procedure axc_xc_type
      end interface
      interface x_functional_dependence
        module procedure axc_functional_dependence
      end interface
      interface uses_gradient
        module procedure q_uses_gradient
      end interface
      interface uses_laplacian
        module procedure q_uses_laplacian
      end interface
      interface xc_energy_density
        module procedure xcd_energy_density_axc
      end interface
      interface xc_derivatives
        module procedure xcd_derivatives_axc
      end interface

      contains

! public routines

      function constructor_axc(xc_type) result(axc)
!doc$ function axc(xc_type) result(axc)
        type(xc_type_obj) :: xc_type
        type(axc_obj) :: axc
!       effects: Constructs a new axc.
!       errors: Passes errors

!cod$ 
        call my(xc_type)

        axc%dependence = x_functional_dependence(xc_type)
        select case(axc%dependence)
        case (FD_DENSITY)
          call my(axc_density(xc_type),axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select

100     call glean(thy(xc_type))

        if (error("Exit axc_mod::constructor_axc")) continue

      end function 

      subroutine update_axc(axc)
!doc$ subroutine update(axc)
        type(axc_obj) :: axc
!       effects: Updates axc.

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call update(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::update_axc")) continue
      end subroutine

      subroutine my_axc(axc)
!doc$ subroutine my(axc)
        type(axc_obj) :: axc

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call my(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::my_axc")) continue
      end subroutine

      subroutine my_new_axc(axci,axc)
!doc$ subroutine my(axci,axc)
        type(axc_obj) :: axci, axc

!cod$
        axc%dependence = axci%dependence
        select case(axc%dependence)
        case (FD_DENSITY)
          call my(axci%axc_density,axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::my_new_axc")) continue
      end subroutine

      function thy_axc(axc) result(axco)
!doc$ function thy(axc) result(axco)
        type(axc_obj) :: axc, axco

!cod$
        axco%dependence = axc%dependence
        select case (axc%dependence)
        case (FD_DENSITY)
          call my(thy(axc%axc_density),axco%axc_density)
          call bequeath(thy(axco%axc_density))
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::thy_axc")) continue
      end function

      subroutine glean_axc(axc)
!doc$ subroutine glean(axc)
        type(axc_obj) :: axc

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call glean(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::glean_axc")) continue
      end subroutine

      subroutine bequeath_axc(axc)
!doc$ subroutine bequeath(axc)
        type(axc_obj) :: axc

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call bequeath(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::bequeath_axc")) continue
      end subroutine
 
      subroutine assign_axc(axc,axc2)
!doc$ subroutine assignment(=)(axc,axc2)
        type(axc_obj), intent(inout) :: axc
        type(axc_obj), intent(in) :: axc2

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          axc%axc_density = axc2%axc_density
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::assign_axc")) continue
      end subroutine

      function axc_ref(axc) result(r)
!doc$ function x_ref(axc) result(r)
        type(axc_obj) :: axc
        integer, dimension(2) :: r
!       effects: Returns the reference counts of the underlying data structure
!                (the density-dependent data structure when dependence = FD_HYBRID).

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          r = x_ref(axc%axc_density)
        case (FD_ORBITAL)
          r = 0
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          r = 0
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::axc_ref")) continue
      end function
 
      function axc_ghost(axc) result(g)
!doc$ function x_ghost(axc) result(g)
        type(axc_obj) :: axc
        type(ghost) :: g
!       effects: Returns the ghost of the underlying data structure
!                (the density-dependent data structure when dependence = FD_HYBRID).

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          g = x_ghost(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::axc_ghost")) continue
      end function

      function axc_xc_type(axc) result(xc_type)
!doc$ function x_xc_type(axc) result(xc_type)
        type(axc_obj) :: axc
        type(xc_type_obj) :: xc_type
!       effects: Returns the xc_type of axc.

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call my(x_xc_type(axc%axc_density),xc_type)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: orbital-dependent functionals are not yet implemented")) goto 100
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) goto 100
        end select
        call bequeath(thy(xc_type))
100     if (error("Exit axc_mod::axc_xc_type")) continue
      end function

      function axc_functional_dependence(axc) result(fd)
!doc$ function x_functional_dependence(axc) result(fd)
        type(axc_obj) :: axc
        integer :: fd
!       effects: Returns the functional dependence of axc.

!cod$
        fd = axc%dependence
      end function

      function q_uses_gradient(axc) result(ug)
!doc$ function uses_gradient(axc) result(ug)
        type(axc_obj) :: axc
        logical :: ug
!       effects: Returns .true. iff the density-dependent functional uses gradients of the density.

!cod$ 
        select case(axc%dependence)
        case (FD_DENSITY)
          ug = uses_gradient(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: query is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::q_uses_gradient")) continue
      end function

      function q_uses_laplacian(axc) result(ul)
!doc$ function uses_laplacian(axc) result(ul)
        type(axc_obj) :: axc
        logical :: ul
!       effects: Returns .true. iff the density-dependent functional uses laplacians of the density.

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          ul = uses_laplacian(axc%axc_density)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: query is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::q_uses_laplacian")) continue
      end function

      subroutine xcd_energy_density_axc(axc,n,dn,lapl,exc)
!doc$ subroutine xc_energy_density(axc,n,dn,lapl,exc)
        type(axc_obj) :: axc
        real(double), dimension(:), pointer :: n, dn, lapl, exc
!       requires: Pointers be nullified or associated.
!       effects: Returns exc with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call xc_energy_density(axc%axc_density,n,dn,lapl,exc)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::xcd_energy_density_axc")) continue
      end subroutine

      subroutine xcd_derivatives_axc(axc,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xc_derivatives(axc,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(axc_obj) :: axc
        real(double), dimension(:), pointer :: n, dn, lapl, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Pointers be nullified or associated.
!       effects: Returns derivatives with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        select case(axc%dependence)
        case (FD_DENSITY)
          call xc_derivatives(axc%axc_density,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        case (FD_ORBITAL)
          if (error(.true.,"ERROR: routine is only valid with density-dependent functionals")) continue
        case (FD_HYBRID)
          if (error(.true.,"ERROR: hybrid functionals are not yet implemented")) continue
        end select
        if (error("Exit axc_mod::xcd_derivatives_axc")) continue
      end subroutine

      end module

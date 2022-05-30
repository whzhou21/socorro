! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module axc_density_mod
!doc$ module axc_density_mod

      use kind_mod
      use error_mod
      use ghost_mod
      use xc_type_mod
      use axc_density_native_mod
      use axc_density_libxc_mod

!     One datatype is defined here: type(axc_density_obj).  

!     axc_density_mod is a wrapper for density-dependent native and libxc exchange-correlation
!     routines that operate on atomic data.

!cod$
      implicit none
      private

      type, public :: axc_density_obj
        private
        integer :: source                                    ! Source of the (density-dependent) functionals.
        type(axc_density_native_obj) :: axcd_native
        type(axc_density_libxc_obj)  :: axcd_libxc
      end type

!doc$
      public :: axc_density
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_xc_type
      public :: uses_gradient
      public :: uses_laplacian
      public :: xc_energy_density
      public :: xc_derivatives

!cod$
      interface axc_density
        module procedure constructor_axcd
      end interface
      interface update
        module procedure update_axcd
      end interface
      interface my
        module procedure my_axcd, my_new_axcd
      end interface
      interface thy
        module procedure thy_axcd
      end interface
      interface glean
        module procedure glean_axcd
      end interface
      interface bequeath
        module procedure bequeath_axcd
      end interface
      interface assignment(=)
        module procedure assign_axcd
      end interface
      interface x_ref
        module procedure axcd_ref
      end interface
      interface x_ghost
        module procedure axcd_ghost
      end interface
      interface x_xc_type
        module procedure axcd_xc_type
      end interface
      interface uses_gradient
        module procedure q_uses_gradient
      end interface
      interface uses_laplacian
        module procedure q_uses_laplacian
      end interface
      interface xc_energy_density
        module procedure xc_energy_density_axcd
      end interface
      interface xc_derivatives
        module procedure xc_derivatives_axcd
      end interface

      contains

! public routines

      function constructor_axcd(xct) result(axcd)
!doc$ function axc_density(xct) result(axcd)
        type(xc_type_obj) :: xct
        type(axc_density_obj) :: axcd
!       effects: Constructs a new axcd.
!       errors: Passes errors

!cod$ 
        call my(xct)

        axcd%source = x_ddf_source(xct)
        select case(axcd%source)
        case (DDFS_NATIVE)
          call my(axc_density_native(xct),axcd%axcd_native)
        case (DDFS_LIBXC)
          call my(axc_density_libxc(xct),axcd%axcd_libxc)
        end select

100     call glean(thy(xct))

        if (error("Exit axc_density_mod::constructor_axcd")) continue

      end function 

      subroutine update_axcd(axcd)
!doc$ subroutine update(axcd)
        type(axc_density_obj) :: axcd
!       effects: Updates axcd.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          call update(axcd%axcd_native)
        case (DDFS_LIBXC)
          call update(axcd%axcd_libxc)
        end select
        if (error("Exit axc_density_mod::update_axcd")) continue
      end subroutine

      subroutine my_axcd(axcd)
!doc$ subroutine my(axcd)
        type(axc_density_obj) :: axcd

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          call my(axcd%axcd_native)
        case (DDFS_LIBXC)
          call my(axcd%axcd_libxc)
        end select
      end subroutine

      subroutine my_new_axcd(axcdi,axcd)
!doc$ subroutine my(axcdi,axcd)
        type(axc_density_obj) :: axcdi, axcd

!cod$
        axcd%source = axcdi%source
        select case (axcd%source)
        case (DDFS_NATIVE)
          call my(axcdi%axcd_native,axcd%axcd_native)
        case (DDFS_LIBXC)
          call my(axcdi%axcd_libxc,axcd%axcd_libxc)
        end select
      end subroutine

      function thy_axcd(axcd) result(axcdo)
!doc$ function thy(axcd) result(axcdo)
        type(axc_density_obj) :: axcd, axcdo

!cod$
        axcdo%source = axcd%source
        select case (axcd%source)
        case (DDFS_NATIVE)
          call my(thy(axcd%axcd_native),axcdo%axcd_native)
          call bequeath(thy(axcdo%axcd_native))
        case (DDFS_LIBXC)
          call my(thy(axcd%axcd_libxc),axcdo%axcd_libxc)
          call bequeath(thy(axcdo%axcd_libxc))
        end select
      end function

      subroutine glean_axcd(axcd)
!doc$ subroutine glean(axcd)
        type(axc_density_obj) :: axcd

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          call glean(axcd%axcd_native)
        case (DDFS_LIBXC)
          call glean(axcd%axcd_libxc)
        end select
      end subroutine

      subroutine bequeath_axcd(axcd)
!doc$ subroutine bequeath(axcd)
        type(axc_density_obj) :: axcd

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          call bequeath(axcd%axcd_native)
        case (DDFS_LIBXC)
          call bequeath(axcd%axcd_libxc)
        end select
      end subroutine
 
      subroutine assign_axcd(axcd,axcd2)
!doc$ subroutine assignment(=)(axcd,axcd2)
        type(axc_density_obj), intent(inout) :: axcd
        type(axc_density_obj), intent(in) :: axcd2

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          axcd%axcd_native = axcd2%axcd_native
        case (DDFS_LIBXC)
          axcd%axcd_libxc = axcd2%axcd_libxc
        end select
      end subroutine

      function axcd_ref(axcd) result(r)
!doc$ function x_ref(axcd) result(r)
        type(axc_density_obj) :: axcd
        integer, dimension(2) :: r
!       effects: Returns %ref and %o%ref of the underlying data structure.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          r = x_ref(axcd%axcd_native)
        case (DDFS_LIBXC)
          r = x_ref(axcd%axcd_libxc)
        end select
      end function
 
      function axcd_ghost(axcd) result(g)
!doc$ function x_ghost(axcd) result(g)
        type(axc_density_obj) :: axcd
        type(ghost) :: g
!       effects: Returns the ghost of the underlying data structure.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          g = x_ghost(axcd%axcd_native)
        case (DDFS_LIBXC)
          g = x_ghost(axcd%axcd_libxc)
        end select
      end function

      function axcd_xc_type(axcd) result(xct)
!doc$ function x_xc_type(axcd) result(xct)
        type(axc_density_obj) :: axcd
        type(xc_type_obj) :: xct
!       effects: Returns the xct of the underlying data structure.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          call my(x_xc_type(axcd%axcd_native),xct)
        case (DDFS_LIBXC)
          call my(x_xc_type(axcd%axcd_libxc),xct)
        end select
        call bequeath(thy(xct))
      end function

      function q_uses_gradient(axcd) result(ug)
!doc$ function uses_gradient(axcd) result(ug)
        type(axc_density_obj) :: axcd
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        select case (axcd%source)
        case (DDFS_NATIVE)
          ug = uses_gradient(axcd%axcd_native)
        case (DDFS_LIBXC)
          ug = uses_gradient(axcd%axcd_libxc)
        end select

      end function

      function q_uses_laplacian(axcd) result(ul)
!doc$ function uses_laplacian(axcd) result(ul)
        type(axc_density_obj) :: axcd
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          ul = uses_laplacian(axcd%axcd_native)
        case (DDFS_LIBXC)
          ul = uses_laplacian(axcd%axcd_libxc)
        end select

      end function

      subroutine xc_energy_density_axcd(axcd,n,dn,lapl,exc)
!doc$ subroutine xc_energy_density(axcd,n,dn,lapl,exc)
        type(axc_density_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lapl, exc
!       requires: Pointers be nullified or associated.
!       effects: Returns the energy density with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
           call xc_energy_density(axcd%axcd_native,n,dn,lapl,exc)
        case (DDFS_LIBXC)
           call xc_energy_density(axcd%axcd_libxc,n,dn,lapl,exc)
        end select 
        if (error("Exit axc_density_mod::xc_energy_density_axcd")) continue
      end subroutine

      subroutine xc_derivatives_axcd(axcd,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xc_derivatives(axcd,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(axc_density_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lapl, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Pointers be nullified or associated.
!       effects: Returns functional derivatives with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        select case (axcd%source)
        case (DDFS_NATIVE)
          call xc_derivatives(axcd%axcd_native,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        case (DDFS_LIBXC)
          call xc_derivatives(axcd%axcd_libxc,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select 
        if (error("Exit axc_density_mod::xc_derivatives_axcd")) continue
      end subroutine

      end module

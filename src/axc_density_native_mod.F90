!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module axc_density_native_mod
!doc$ module axc_density_native_mod

!     One datatype is available here: type(axc_density_native_obj).

!     axc_density_native_mod encapsulates native routines used to compute density-dependent
!     exchange-correlation quantities for atomic data.

      use kind_mod
      use error_mod
      use ghost_mod
      use xc_type_mod
      use xc_functional_mod

!cod$
      implicit none
      private 

      real(double), parameter :: mindens = 1.0e-20_double

      type :: axc_density_native_rep
        integer :: ref                                        ! reference count
        type(ghost) :: g                                      ! ghost
        integer :: etype                                      ! exchange type
        integer :: ctype                                      ! correlation type
        integer :: derivative_method                          ! method for computing derivitives
        logical :: uses_gradient                              ! indicates whether or not the gradient is used
        logical :: uses_laplacian                             ! indicates whether or not the laplacian is used
        type(xc_type_obj) :: xct                              ! xc_type object
      end type

      type, public :: axc_density_native_obj
        private
        integer :: ref
        type(axc_density_native_rep), pointer :: o
      end type

      type :: dens_der
        real(double), dimension(:), pointer :: n
        real(double), dimension(:), pointer :: dn
        real(double), dimension(:), pointer :: lapl
      end type

      public :: axc_density_native
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

      interface axc_density_native
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
        module procedure axcd_energy_density_native
      end interface
      interface xc_derivatives
        module procedure axcd_derivatives_native
      end interface

      contains

! public routines

      function constructor_axcd(xct) result(axcd)
!doc$ function axc_density_native(xct) result(axcd)
        type(xc_type_obj) :: xct
        type(axc_density_native_obj) :: axcd
!       effects: Constructs a new xcd.
!       errors: Passes errors

!cod$ 
        call my(xct)

        axcd%ref = 0
        allocate( axcd%o )
        axcd%o%ref = 0
        axcd%o%g = x_ghost()

        call my(xct,axcd%o%xct)

        axcd%o%etype = x_exchange_type(axcd%o%xct)
        axcd%o%ctype = x_correlation_type(axcd%o%xct)
        axcd%o%derivative_method = x_derivative_method(axcd%o%xct)
        axcd%o%uses_gradient = uses_gradient(axcd%o%xct)
        axcd%o%uses_laplacian = uses_laplacian(axcd%o%xct)

100     call glean(thy(xct))

        if (error("Exit axc_density_native_mod::constructor_axcd")) continue

      end function 

      subroutine update_axcd(axcd)
!doc$ subroutine update(axcd)
        type(axc_density_native_obj) :: axcd
!       modifies: axcd
!       effects: Updates axcd.
!       notes: axcd%o%xct is not allowed to change.

!cod$
        call my(axcd)

        continue

100     call glean(thy(axcd))

      end subroutine

      subroutine my_axcd(axcd)
!doc$ subroutine my(axcd)
        type(axc_density_native_obj) :: axcd

!cod$
        axcd%ref = axcd%ref + 1
        axcd%o%ref = axcd%o%ref + 1
      end subroutine

      subroutine my_new_axcd(axcdi,axcd)
!doc$ subroutine my(axcdi,axcd)
        type(axc_density_native_obj) :: axcdi, axcd

!cod$
        axcd%ref = 1
        axcd%o => axcdi%o
        axcd%o%ref = axcd%o%ref + 1
      end subroutine
      
      function thy_axcd(axcd) result(axcdo)
!doc$ function thy(axcd) result(axcdo)
        type(axc_density_native_obj) :: axcd, axcdo

!cod$
        axcd%ref = axcd%ref - 1
        axcd%o%ref = axcd%o%ref - 1
        axcdo%ref = axcd%ref
        axcdo%o => axcd%o
      end function

      subroutine glean_axcd(axcd)
!doc$ subroutine glean(axcd)
        type(axc_density_native_obj) :: axcd

!cod$
        if (axcd%o%ref < 1) then
          call glean(thy(axcd%o%xct))
          deallocate( axcd%o )
        end if
      end subroutine

      subroutine bequeath_axcd(axcd)
!doc$ subroutine bequeath(axcd)
        type(axc_density_native_obj) :: axcd

!cod$
        continue
      end subroutine

      subroutine assign_axcd(axcd,axcd2)
!doc$ subroutine assign(axcd,axcd2)
        type(axc_density_native_obj), intent(inout) :: axcd
        type(axc_density_native_obj), intent(in) :: axcd2

!cod$
        type(axc_density_native_obj) :: axcdt
        call my(axcd2)
        axcdt%o => axcd%o
        axcd%o%ref = axcd%o%ref - axcd%ref
        axcd%o => axcd2%o
        axcd%o%ref = axcd%o%ref + axcd%ref
        call glean(axcdt)
        call glean(thy(axcd2))
      end subroutine

      function axcd_ref(axcd) result(r)
!doc$ function x_ref(axcd) result(r)
        type(axc_density_native_obj) :: axcd
        integer, dimension(2) :: r
!       effects: Returns axcd%ref and axcd%o%ref.

!cod$
        r(1) = axcd%ref
        r(2) = axcd%o%ref
        call glean(axcd)
      end function

      function axcd_ghost(axcd) result(g)
!doc$ function x_ghost(axcd) result(g)
        type(axc_density_native_obj) :: axcd
        type(ghost) :: g
!       effects: Returns the ghost of axcd.

!cod$
        call my(axcd)
        g = axcd%o%g
        call glean(thy(axcd))
      end function

      function axcd_xc_type(axcd) result(xct)
!doc$ function x_xc_type(axcd) result(xct)
        type(axc_density_native_obj) :: axcd
        type(xc_type_obj) :: xct
!       effects: Returns the xc_type of axcd.

!cod$
        call my(axcd)
        call my(axcd%o%xct,xct)
        call glean(thy(axcd))
        call bequeath(thy(xct))
      end function

      function q_uses_gradient(axcd) result(ug)
!doc$ function uses_gradient(axcd) result(ug)
        type(axc_density_native_obj) :: axcd
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        call my(axcd)
        ug = axcd%o%uses_gradient 
        call glean(thy(axcd))
      end function

      function q_uses_laplacian(axcd) result(ul)
!doc$ function uses_laplacian(axcd) result(ul)
        type(axc_density_native_obj) :: axcd
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        call my(axcd)
        ul = axcd%o%uses_laplacian
        call glean(thy(axcd))
      end function

      subroutine axcd_energy_density_native(axcd,n,dn,lapl,exc)
!doc$ subroutine xc_energy_density(axcd,n,dn,lapl,exc)
        type(axc_density_native_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lapl, exc
!       requires: Pointers be nullified or associated.
!       effects: Returns exc with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        type(dens_der) :: nder

        call my(axcd)

        call form_nder_i(n,dn,lapl,nder)
        call xce_analytic(axcd%o%etype,axcd%o%ctype,nder%n,nder%dn,nder%lapl,exc) ; if (error()) goto 100
        call kill_nder_i(nder)

        call glean(thy(axcd))

100     if (error("Exit axc_density_native_mod::axcd_energy_density_native")) continue

      end subroutine

      subroutine axcd_derivatives_native(axcd,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xc_derivatives(axcd,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(axc_density_native_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lapl, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Pointers be nullified or associated.
!       effects: Returns derivatives with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        type(dens_der) :: nder

        call my(axcd)

        call form_nder_i(n,dn,lapl,nder)
        select case (axcd%o%derivative_method)
        case (XCD_ANALYTICAL)
          call xcd_analytic(axcd%o%etype,axcd%o%ctype,nder%n,nder%dn,nder%lapl,dfxcdn,dfxcddnodn,dfxcdlapl) ; if (error()) goto 100
        case (XCD_NUMERICAL)
          call xcd_numerical_i(axcd,nder%n,nder%dn,nder%lapl,dfxcdn,dfxcddnodn,dfxcdlapl) ; if (error()) goto 100
        end select
        call kill_nder_i(nder)

        call glean(thy(axcd))

100     if (error("Exit axc_density_native_mod::axcd_derivatives_native")) continue

      end subroutine

! private routines

      subroutine xcd_numerical_i(axcd,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(axc_density_native_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lp, dfxcdn, dfxcddnodn, dfxcdlapl

        real(double), parameter :: maxdelta = 1.0e-4_double
        real(double), dimension(:), allocatable :: delta
        real(double), dimension(:), pointer :: exc_md, exc_pd, md, pd

        call my(axcd)

        nullify( exc_md, exc_pd, md, pd )

        allocate( delta(size(n)) )
        allocate( exc_md(size(n)), exc_pd(size(n)), md(size(n)), pd(size(n)) )

        ! derivatives wrt density
        delta = n/100.0_double
        where (delta >= maxdelta) delta = maxdelta
        pd = n + delta
        call xce_analytic(axcd%o%etype,axcd%o%ctype,pd,dn,lp,exc_pd) ; if (error()) goto 100
        md = n - delta
        call xce_analytic(axcd%o%etype,axcd%o%ctype,md,dn,lp,exc_md) ; if (error()) goto 100
        dfxcdn = 0.5_double*(exc_pd*pd - exc_md*md)/delta

        ! derivatives wrt |gradient density|
        if (axcd%o%uses_gradient) then
          delta = dn/100.0_double
          where (delta >= maxdelta) delta = maxdelta
          pd = dn + delta
          call xce_analytic(axcd%o%etype,axcd%o%ctype,n,pd,lp,exc_pd) ; if (error()) goto 100
          md = dn - delta
          call xce_analytic(axcd%o%etype,axcd%o%ctype,n,md,lp,exc_md) ; if (error()) goto 100
          dfxcddnodn = 0.5_double*(exc_pd*n - exc_md*n)/delta/dn
        end if

        ! derivatives wrt laplacian
        if (axcd%o%uses_laplacian) then
          delta = abs(lp/100.0_double)
          where (delta >= maxdelta) delta = maxdelta
          pd = lp + delta
          call xce_analytic(axcd%o%etype,axcd%o%ctype,n,dn,pd,exc_pd) ; if (error()) goto 100
          md = lp - delta
          call xce_analytic(axcd%o%etype,axcd%o%ctype,n,dn,md,exc_md) ; if (error()) goto 100
          dfxcdlapl = 0.5_double*(exc_pd*n - exc_md*n)/delta
        end if

100     if (allocated( delta )) deallocate( delta )
        if (associated( exc_md )) deallocate( exc_md )
        if (associated( exc_pd )) deallocate( exc_pd )
        if (associated( md )) deallocate( md )
        if (associated( pd )) deallocate( pd )

        call glean(thy(axcd))

        if (error("Exit axc_density_native_mod::xcd_numerical_i")) continue

      end subroutine

      subroutine form_nder_i(n,dn,lapl,nder)
        real(double), dimension(:), pointer :: n, dn, lapl
        type(dens_der) :: nder

        nullify( nder%n, nder%dn, nder%lapl )

        allocate( nder%n(size(n)) )
        nder%n = n
        where (nder%n < mindens) nder%n = mindens
        if (associated( dn )) then
          allocate( nder%dn(size(dn)) )
          nder%dn = dn
          where (nder%n <= mindens) nder%dn = mindens
        end if
        if (associated( lapl )) then
          allocate( nder%lapl(size(lapl)) )
          nder%lapl = lapl
          where (nder%n <= mindens) nder%lapl = mindens
        end if

      end subroutine

      subroutine kill_nder_i(nder)
        type(dens_der) :: nder

        if (associated( nder%n )) deallocate(nder%n)
        if (associated( nder%dn )) deallocate(nder%dn)
        if (associated( nder%lapl )) deallocate(nder%lapl)

      end subroutine

      subroutine own_i(xcd)
        type(axc_density_native_obj) :: xcd
        type(axc_density_native_obj) :: xcd_t
        if (xcd%ref < xcd%o%ref) then
          allocate( xcd_t%o )
          xcd_t%o%ref               = 0
          xcd_t%o%g                 = xcd%o%g
          xcd_t%o%etype             = xcd%o%etype
          xcd_t%o%ctype             = xcd%o%ctype
          xcd_t%o%derivative_method = xcd%o%derivative_method
          xcd_t%o%uses_gradient     = xcd%o%uses_gradient
          xcd_t%o%uses_laplacian    = xcd%o%uses_laplacian
          call my(xcd%o%xct,xcd_t%o%xct)
          xcd%o%ref = xcd%o%ref - xcd%ref
          xcd%o => xcd_t%o
          xcd%o%ref = xcd%o%ref + xcd%ref
        end if
      end subroutine

      end module

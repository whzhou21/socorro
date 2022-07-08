!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module xc_density_native_mod
!doc$ module xc_density_native_mod

!     One datatype is available here: type(xc_density_native_obj).

!     xc_density_native_mod encapsulates native routines used to compute density-dependent
!     exchange-correlation quantities for grid data.

      use kind_mod
      use mpi_mod
      use error_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use lattice_mod
      use symmetry_mod
      use xc_type_mod
      use xc_functional_mod

!cod$
      implicit none
      private 

      real(double), parameter :: mindens = 1.0e-20_double

      type :: xc_density_native_rep
        integer :: ref                                        ! reference count
        type(ghost) :: g                                      ! ghost
        integer :: etype                                      ! exchange type
        integer :: ctype                                      ! correlation type
        integer :: derivative_method                          ! method for computing derivitives
        integer :: potential_method                           ! method for computing the potential
        logical :: uses_gradient                              ! indicates whether or not the gradient is used
        logical :: uses_laplacian                             ! indicates whether or not the laplacian is used
        real(double) :: volume                                ! cell volume
        type(xc_type_obj) :: xct                              ! xc_type object
        type(layout_obj) :: lay                               ! layout object
        type(space_group_obj) :: sg                           ! space group object
        real(double), dimension(:,:,:), pointer :: gx         ! x component of reciprocal-space mesh values
        real(double), dimension(:,:,:), pointer :: gy         ! y component of reciprocal-space mesh values
        real(double), dimension(:,:,:), pointer :: gz         ! z component of reciprocal-space mesh values
      end type

      type, public :: xc_density_native_obj
        private
        integer :: ref
        type(xc_density_native_rep), pointer :: o
      end type

      type :: dens_der
        real(double), dimension(:,:,:), pointer :: n
        real(double), dimension(:,:,:), pointer :: dn
        real(double), dimension(:,:,:), pointer :: dnx
        real(double), dimension(:,:,:), pointer :: dny
        real(double), dimension(:,:,:), pointer :: dnz
        real(double), dimension(:,:,:), pointer :: lapl
        real(double), dimension(:,:,:), pointer :: d2nxx
        real(double), dimension(:,:,:), pointer :: d2nxy
        real(double), dimension(:,:,:), pointer :: d2nxz
        real(double), dimension(:,:,:), pointer :: d2nyy
        real(double), dimension(:,:,:), pointer :: d2nyz
        real(double), dimension(:,:,:), pointer :: d2nzz
      end type

      public :: xc_density_native
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

      interface xc_density_native
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
        module procedure xcd_energy_native
      end interface
      interface xc_potential
        module procedure xcd_potential_native
      end interface
      interface xc_energy_and_potential
        module procedure xcd_energy_and_potential_native
      end interface
      interface xc_grid_pressure
        module procedure xcd_grid_pressure_native
      end interface
      interface xc_grid_stress_tensor
        module procedure xcd_grid_stress_tensor_native
      end interface

      contains

! public routines

      function constructor_xcd(xct,lay,sg) result(xcd)
!doc$ function xc_density_native(xct,lay,sg) result(xcd)
        type(xc_type_obj) :: xct
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_density_native_obj) :: xcd
!       effects: Constructs a new xcd.
!       errors: Passes errors

!cod$ 
        call my(xct)
        call my(lay)
        call my(sg)

        xcd%ref = 0
        allocate( xcd%o )
        xcd%o%ref = 0
        xcd%o%g = x_ghost()

        call my(xct,xcd%o%xct)
        call my(lay,xcd%o%lay)
        call my(sg,xcd%o%sg)

        xcd%o%etype = x_exchange_type(xcd%o%xct)
        xcd%o%ctype = x_correlation_type(xcd%o%xct)
        xcd%o%derivative_method = x_derivative_method(xcd%o%xct)
        xcd%o%potential_method = x_potential_method(xcd%o%xct)
        xcd%o%uses_gradient = uses_gradient(xcd%o%xct)
        xcd%o%uses_laplacian = uses_laplacian(xcd%o%xct)

        xcd%o%volume = x_cell_volume(x_lattice(xcd%o%lay))

        nullify( xcd%o%gx )
        nullify( xcd%o%gy )
        nullify( xcd%o%gz )
        if (xcd%o%uses_gradient .or. xcd%o%uses_laplacian) then
          call fmesh(xcd%o%gx,xcd%o%gy,xcd%o%gz,xcd%o%lay,D_TYPE,SGROUP)
        end if

100     call glean(thy(xct))
        call glean(thy(lay))
        call glean(thy(sg))

        if (error("Exit xc_density_native_mod::constructor_xcd")) continue

      end function 

      subroutine update_xcd(xcd,lay,sg)
!doc$ subroutine update(xcd,lay,sg)
        type(xc_density_native_obj) :: xcd
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
!       modifies: xcd
!       effects: Updates xcd.
!       note: xcd%o%xct is not allowed to change.

!cod$
        logical :: layout_change, space_group_change

        call my(xcd)
        call my(lay)
        call my(sg)

        layout_change = (x_ghost(lay) /= x_ghost(xcd%o%lay))
        space_group_change = (x_ghost(sg) /= x_ghost(xcd%o%sg))

        if (layout_change .or. space_group_change) then

          call own_i(xcd)
          xcd%o%g = x_ghost()

          if (layout_change) then
            xcd%o%lay = lay
            if (xcd%o%uses_gradient .or. xcd%o%uses_laplacian) then
              deallocate( xcd%o%gx )
              deallocate( xcd%o%gy )
              deallocate( xcd%o%gz )
              call fmesh(xcd%o%gx,xcd%o%gy,xcd%o%gz,xcd%o%lay,D_TYPE,SGROUP)
            else
              nullify( xcd%o%gx )
              nullify( xcd%o%gy )
              nullify( xcd%o%gz )
            end if
          end if

          if (space_group_change) xcd%o%sg = sg

        end if

100     call glean(thy(xcd))
        call glean(thy(lay))
        call glean(thy(sg))

      end subroutine

      subroutine my_xcd(xcd)
!doc$ subroutine my(xcd)
        type(xc_density_native_obj) :: xcd

!cod$
        xcd%ref = xcd%ref + 1
        xcd%o%ref = xcd%o%ref + 1
      end subroutine

      subroutine my_new_xcd(xcdi,xcd)
!doc$ subroutine my(xcdi,xcd)
        type(xc_density_native_obj) :: xcdi, xcd

!cod$
        xcd%ref = 1
        xcd%o => xcdi%o
        xcd%o%ref = xcd%o%ref + 1
      end subroutine
      
      function thy_xcd(xcd) result(xcdo)
!doc$ function thy(xcd) result(xcdo)
        type(xc_density_native_obj) :: xcd, xcdo

!cod$
        xcd%ref = xcd%ref - 1
        xcd%o%ref = xcd%o%ref - 1
        xcdo%ref = xcd%ref
        xcdo%o => xcd%o
      end function

      subroutine glean_xcd(xcd)
!doc$ subroutine glean(xcd)
        type(xc_density_native_obj) :: xcd

!cod$
        if (xcd%o%ref < 1) then
          call glean(thy(xcd%o%xct))
          call glean(thy(xcd%o%lay))
          call glean(thy(xcd%o%sg))
          if (associated( xcd%o%gx )) deallocate( xcd%o%gx )
          if (associated( xcd%o%gy )) deallocate( xcd%o%gy )
          if (associated( xcd%o%gz )) deallocate( xcd%o%gz )
          deallocate( xcd%o )
        end if
      end subroutine

      subroutine bequeath_xcd(xcd)
!doc$ subroutine bequeath(xcd)
        type(xc_density_native_obj) :: xcd

!cod$
        continue
      end subroutine

      subroutine assign_xcd(xcd,xcd2)
!doc$ subroutine assign(xcd,xcd2)
        type(xc_density_native_obj), intent(inout) :: xcd
        type(xc_density_native_obj), intent(in) :: xcd2

!cod$
        type(xc_density_native_obj) :: xcdt
        call my(xcd2)
        xcdt%o => xcd%o
        xcd%o%ref = xcd%o%ref - xcd%ref
        xcd%o => xcd2%o
        xcd%o%ref = xcd%o%ref + xcd%ref
        call glean(xcdt)
        call glean(thy(xcd2))
      end subroutine

      function xcd_ref(xcd) result(r)
!doc$ function x_ref(xcd) result(r)
        type(xc_density_native_obj) :: xcd
        integer, dimension(2) :: r
!       effects: Returns xcd%ref and xcd%o%ref.

!cod$
        r(1) = xcd%ref
        r(2) = xcd%o%ref
        call glean(xcd)
      end function

      function xcd_ghost(xcd) result(g)
!doc$ function x_ghost(xcd) result(g)
        type(xc_density_native_obj) :: xcd
        type(ghost) :: g
!       effects: Returns the ghost of xcd.

!cod$
        call my(xcd)
        g = xcd%o%g
        call glean(thy(xcd))
      end function

      function xcd_xc_type(xcd) result(xct)
!doc$ function x_xc_type(xcd) result(xct)
        type(xc_density_native_obj) :: xcd
        type(xc_type_obj) :: xct
!       effects: Returns the xc_type of xcd.

!cod$
        call my(xcd)
        call my(xcd%o%xct,xct)
        call glean(thy(xcd))
        call bequeath(thy(xct))
      end function

      function xcd_layout(xcd) result(lay)
!doc$ function x_layout(xcd) result(lay)
        type(xc_density_native_obj) :: xcd
        type(layout_obj) :: lay
!       effects: Returns the lay of xcd.

!cod$
        call my(xcd)
        call my(xcd%o%lay,lay)
        call glean(thy(xcd))
        call bequeath(thy(lay))
      end function

      function xcd_space_group(xcd) result(sg)
!doc$ function x_space_group(xcd) result(sg)
        type(xc_density_native_obj) :: xcd
        type(space_group_obj) :: sg
!       effects: Returns the sg of xcd.

!cod$
        call my(xcd)
        call my(xcd%o%sg,sg)
        call glean(thy(xcd))
        call bequeath(thy(sg))
      end function

      function q_uses_gradient(xcd) result(ug)
!doc$ function uses_gradient(xcd) result(ug)
        type(xc_density_native_obj) :: xcd
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        call my(xcd)
        ug = xcd%o%uses_gradient 
        call glean(thy(xcd))
      end function

      function q_uses_laplacian(xcd) result(ul)
!doc$ function uses_laplacian(xcd) result(ul)
        type(xc_density_native_obj) :: xcd
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        call my(xcd)
        ul = xcd%o%uses_laplacian
        call glean(thy(xcd))
      end function

      function xcd_energy_native(xcd,n_g) result(fxc)
!doc$ function xc_energy(xcd,n_g) result(fxc)
        type(xc_density_native_obj) :: xcd
        type(grid_obj) :: n_g
        real(double) :: fxc
!       requires: n_g data be symmetrized.
!       effects: Returns the exchange-correlation energy.
!       errors: Passes errors.

!cod$
        real(double) :: e_local, e_global
        real(double), dimension(:,:,:), pointer :: n, exc
        type(layout_obj) :: lay
        type(dens_der) :: nder

        call my(xcd)
        call my(n_g)

        call my(x_layout(n_g),lay)

        nullify( n, exc )

        call take(n,n_g,RD_KIND)

        call alloc(exc,lay,D_TYPE,SGROUP)

        call form_nder_i(xcd,lay,n,nder)
        call xce_analytic(xcd%o%etype,xcd%o%ctype,nder%n,nder%dn,nder%lapl,exc) ; if (error()) goto 100
        call kill_nder_i(nder)

        e_local = sum(n*exc)
        call allreduce(CONFIG,MPI_SUM,e_local,e_global)
        fxc = e_global*(xcd%o%volume/product(x_dims(lay)))

        call put(n,n_g,RD_KIND)

100     if (associated( n )) deallocate( n )
        if (associated( exc )) deallocate( exc )

        call glean(thy(lay))

        call glean(thy(xcd))
        call glean(thy(n_g))

        if (error("Exit xc_density_native_mod::xcd_energy_native")) continue

      end function

      function xcd_potential_native(xcd,n_g) result(vxc_g)
!doc$ function xc_potential(xcd,n_g) result(vxc_g)
        type(xc_density_native_obj) :: xcd
        type(grid_obj) :: n_g
        type(grid_obj) :: vxc_g
!       requires: n_g data be symmetrized with respect to xcd%o%sg. Functional be density dependent.
!       effects: Returns the symmetrized/filtered exchange-correlation potential.
!       errors: Passes errors.

!cod$ 
        real(double), dimension(:,:,:), pointer :: n, vxc, dfxcdn, dfxcddnodn, dfxcdlapl
        type(layout_obj) :: lay
        type(dens_der) :: nder

        call my(xcd)
        call my(n_g) 

        nullify( n, vxc, dfxcdn, dfxcddnodn, dfxcdlapl )

        call my(x_layout(n_g),lay)
        call my(grid(lay,SGROUP),vxc_g)

        call take(n,n_g,RD_KIND)

        call alloc(vxc,lay,D_TYPE,SGROUP)
        call alloc(dfxcdn,lay,D_TYPE,SGROUP)
        if (xcd%o%uses_gradient) call alloc(dfxcddnodn,lay,D_TYPE,SGROUP)
        if (xcd%o%uses_laplacian) call alloc(dfxcdlapl,lay,D_TYPE,SGROUP)

        call form_nder_i(xcd,lay,n,nder)
        select case (xcd%o%derivative_method)
        case (XCD_ANALYTICAL)
          call xcd_analytic(xcd%o%etype,xcd%o%ctype,nder%n,nder%dn,nder%lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        case (XCD_NUMERICAL)
          call xcd_numerical_i(xcd,nder%n,nder%dn,nder%lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select
        if (error()) goto 100
        select case (xcd%o%potential_method)
        case (XCD_SIMPLE)
          vxc = dfxcdn
        case (XCD_WHITE_BIRD)
          call white_bird_i(xcd,lay,nder,dfxcdn,dfxcddnodn,dfxcdlapl,vxc)
        end select
        call kill_nder_i(nder)

        call put(vxc,vxc_g,RD_KIND)
        call symmetrize_grid(xcd%o%sg,vxc_g)
        call filter(vxc_g)

        call put(n,n_g,RD_KIND)

100     if (associated( n )) deallocate( n )
        if (associated( vxc )) deallocate( vxc )
        if (associated( dfxcdn )) deallocate( dfxcdn )
        if (associated( dfxcddnodn )) deallocate( dfxcddnodn )
        if (associated( dfxcdlapl )) deallocate( dfxcdlapl )

        call glean(thy(lay))
        call bequeath(thy(vxc_g))
      
        call glean(thy(xcd))
        call glean(thy(n_g))

        if (error("Exit xc_density_native_mod::xcd_potential_native")) continue

      end function

      function xcd_energy_and_potential_native(xcd,n_g,fxc) result(vxc_g)
!doc$ function xc_energy_and_potential(xcd,n_g,fxc) result(vxc_g)
        type(xc_density_native_obj) :: xcd
        type(grid_obj) :: n_g
        real(double), intent(out) :: fxc
        type(grid_obj) :: vxc_g
!       requires: n_g data be symmetrized.
!       effects: Returns the symmetrized/filtered exchange-correlation potential and energy.
!       errors: Passes errors.

!cod$
        real(double) :: e_local, e_global
        real(double), dimension(:,:,:), pointer :: n, exc, vxc, dfxcdn, dfxcddnodn, dfxcdlapl
        type(layout_obj) :: lay
        type(dens_der) :: nder

        call my(xcd)
        call my(n_g)

        call my(x_layout(n_g),lay)
        call my(grid(lay,SGROUP),vxc_g)

        nullify( n, exc, vxc, dfxcdn, dfxcddnodn, dfxcdlapl )

        call take(n,n_g,RD_KIND)

        call alloc(exc,lay,D_TYPE,SGROUP)
        call alloc(vxc,lay,D_TYPE,SGROUP)
        call alloc(dfxcdn,lay,D_TYPE,SGROUP)
        if (xcd%o%uses_gradient) call alloc(dfxcddnodn,lay,D_TYPE,SGROUP)
        if (xcd%o%uses_laplacian) call alloc(dfxcdlapl,lay,D_TYPE,SGROUP)

        call form_nder_i(xcd,lay,n,nder)
        select case (xcd%o%derivative_method)
        case (XCD_ANALYTICAL)
          call xcq_analytic(xcd%o%etype,xcd%o%ctype,nder%n,nder%dn,nder%lapl,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        case (XCD_NUMERICAL)
          call xcq_numerical_i(xcd,nder%n,nder%dn,nder%lapl,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select
        if (error()) goto 100
        select case (xcd%o%potential_method)
        case (XCD_SIMPLE)
          vxc = dfxcdn
        case (XCD_WHITE_BIRD)
          call white_bird_i(xcd,lay,nder,dfxcdn,dfxcddnodn,dfxcdlapl,vxc)
        end select
        call kill_nder_i(nder)

        call put(vxc,vxc_g,RD_KIND)
        call symmetrize_grid(xcd%o%sg,vxc_g)
        call filter(vxc_g)

        e_local = sum(n*exc)
        call allreduce(CONFIG,MPI_SUM,e_local,e_global)
        fxc = e_global*(xcd%o%volume/product(x_dims(lay)))

        call put(n,n_g,RD_KIND)

100     if (associated( n )) deallocate( n )
        if (associated( exc )) deallocate( exc )
        if (associated( vxc )) deallocate( vxc )
        if (associated( dfxcdn )) deallocate( dfxcdn )
        if (associated( dfxcddnodn )) deallocate( dfxcddnodn )
        if (associated( dfxcdlapl )) deallocate( dfxcdlapl )

        call glean(thy(lay))
        call bequeath(thy(vxc_g))

        call glean(thy(xcd))
        call glean(thy(n_g))

        if (error("Exit xc_density_native_mod::xcd_energy_and_potential_native")) continue

      end function

      subroutine xcd_grid_pressure_native(xcd,n,p)
!doc$ subroutine xc_grid_pressure(xcd,n,p)
        type(xc_density_native_obj) :: xcd
        type(grid_obj) :: n
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the exchange-correlation grid potential.

!cod$
        real(double), dimension(:,:), pointer :: s

        call my(xcd)
        call my(n)

        allocate( s(3,3) )
        call xcd_grid_stress_tensor_native(xcd,n,s) ; if (error()) goto 100
        p = -(s(1,1) + s(2,2) + s(3,3))/3.0_double

100     if (associated( s )) deallocate( s )

        call glean(thy(xcd))
        call glean(thy(n))

        if (error("Exit xc_density_native_mod::xcd_grid_pressure_native")) continue

      end subroutine

      subroutine xcd_grid_stress_tensor_native(xcd,n_g,s)
!doc$ subroutine xc_grid_stress_tensor(xcd,n_g,s)
        type(xc_density_native_obj) :: xcd
        type(grid_obj) :: n_g
        real(double), dimension(:,:), intent(out) :: s
!       modifies: s
!       effects: Returns stress tensor contributions due to the exchange-correlation grid potential.

!cod$
        logical :: stress
        real(double) :: sum_r
        real(double), dimension(:,:), allocatable :: s_local
        real(double), dimension(:,:,:), pointer :: n, exc, vxc, dfxcdn, dfxcddnodn, dfxcdlapl
        type(layout_obj) :: lay
        type(dens_der) :: nder

        call my(xcd)
        call my(n_g)

        nullify( n, exc, vxc, dfxcdn, dfxcddnodn, dfxcdlapl )

        call my(x_layout(n_g),lay)

        call take(n,n_g,RD_KIND)

        allocate( s_local(3,3) )
        call alloc(exc,lay,D_TYPE,SGROUP)
        call alloc(vxc,lay,D_TYPE,SGROUP)
        call alloc(dfxcdn,lay,D_TYPE,SGROUP)
        if (xcd%o%uses_gradient) call alloc(dfxcddnodn,lay,D_TYPE,SGROUP)
        if (xcd%o%uses_laplacian) call alloc(dfxcdlapl,lay,D_TYPE,SGROUP)

        stress = .true.
        call form_nder_i(xcd,lay,n,nder,stress)
        select case (xcd%o%derivative_method)
        case (XCD_ANALYTICAL)
          call xcq_analytic(xcd%o%etype,xcd%o%ctype,nder%n,nder%dn,nder%lapl,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        case (XCD_NUMERICAL)
          call xcq_numerical_i(xcd,nder%n,nder%dn,nder%lapl,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select
        if (error()) goto 100
        select case (xcd%o%potential_method)
        case (XCD_SIMPLE)
          vxc = dfxcdn
        case (XCD_WHITE_BIRD)
          call white_bird_i(xcd,lay,nder,dfxcdn,dfxcddnodn,dfxcdlapl,vxc)
        end select
        sum_r = sum(n*(exc - vxc))                                           ! contributions from the density
        s_local(1,1) = sum_r
        s_local(1,2) = 0.0_double
        s_local(1,3) = 0.0_double
        s_local(2,2) = sum_r
        s_local(2,3) = 0.0_double
        s_local(3,3) = sum_r
        if (xcd%o%uses_gradient) then                                        ! contributions from the gradient
          s_local(1,1) = s_local(1,1) - sum(dfxcddnodn*nder%dnx**2)
          s_local(1,2) = s_local(1,2) - sum(dfxcddnodn*nder%dnx*nder%dny)
          s_local(1,3) = s_local(1,3) - sum(dfxcddnodn*nder%dnx*nder%dnz)
          s_local(2,2) = s_local(2,2) - sum(dfxcddnodn*nder%dny**2)
          s_local(2,3) = s_local(2,3) - sum(dfxcddnodn*nder%dny*nder%dnz)
          s_local(3,3) = s_local(3,3) - sum(dfxcddnodn*nder%dnz**2) 
        end if
        if ( xcd%o%uses_laplacian ) then                                      ! contributions from the laplacian
          s_local(1,1) = s_local(1,1) - sum(dfxcdlapl*nder%d2nxx)
          s_local(1,2) = s_local(1,2) - sum(dfxcdlapl*nder%d2nxy)
          s_local(1,3) = s_local(1,3) - sum(dfxcdlapl*nder%d2nxz)
          s_local(2,2) = s_local(2,2) - sum(dfxcdlapl*nder%d2nyy)
          s_local(2,3) = s_local(2,3) - sum(dfxcdlapl*nder%d2nyz)
          s_local(3,3) = s_local(3,3) - sum(dfxcdlapl*nder%d2nzz)   
        end if               
        call kill_nder_i(nder)

        s_local(2,1) = s_local(1,2)
        s_local(3,1) = s_local(1,3)
        s_local(3,2) = s_local(2,3)
        call allreduce(CONFIG,MPI_SUM,s_local,s)
        s = s/product(x_dims(lay))

100     if (allocated( s_local )) deallocate( s_local )
        if (associated( n )) deallocate( n )
        if (associated( exc )) deallocate( exc )
        if (associated( vxc )) deallocate( vxc )
        if (associated( dfxcdn )) deallocate( dfxcdn )
        if (associated( dfxcddnodn )) deallocate( dfxcddnodn )
        if (associated( dfxcdlapl )) deallocate( dfxcdlapl )

        call glean(thy(lay))

        call glean(thy(xcd))
        call glean(thy(n_g))

        if (error("Exit xc_density_native_mod::xcd_grid_stress_tensor_native")) continue

      end subroutine

! private routines

      subroutine xcd_numerical_i(xcd,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(xc_density_native_obj) :: xcd
        real(double), dimension(:,:,:), pointer :: n, dn, lp, dfxcdn, dfxcddnodn, dfxcdlapl

        integer :: n1, n2, n3
        real(double), parameter :: maxdelta = 1.0e-4_double
        real(double), dimension(:,:,:), allocatable :: delta
        real(double), dimension(:,:,:), pointer :: exc_md, exc_pd, md, pd

        call my(xcd)

        nullify( exc_md, exc_pd, md, pd )

        n1 = size(n,1)
        n2 = size(n,2)
        n3 = size(n,3)

        allocate( delta(n1,n2,n3) )
        allocate( exc_md(n1,n2,n3), exc_pd(n1,n2,n3), md(n1,n2,n3), pd(n1,n2,n3) )

        ! derivatives wrt density
        delta = n/100.0_double
        where (delta >= maxdelta) delta = maxdelta
        pd = n + delta
        call xce_analytic(xcd%o%etype,xcd%o%ctype,pd,dn,lp,exc_pd) ; if (error()) goto 100
        md = n - delta
        call xce_analytic(xcd%o%etype,xcd%o%ctype,md,dn,lp,exc_md) ; if (error()) goto 100
        dfxcdn = 0.5_double*(exc_pd*pd - exc_md*md)/delta

        ! derivatives wrt |gradient density|
        if (xcd%o%uses_gradient) then
          delta = dn/100.0_double
          where (delta >= maxdelta) delta = maxdelta
          pd = dn + delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,pd,lp,exc_pd) ; if (error()) goto 100
          md = dn - delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,md,lp,exc_md) ; if (error()) goto 100
          dfxcddnodn = 0.5_double*(exc_pd*n - exc_md*n)/delta/dn
        end if

        ! derivatives wrt laplacian
        if (xcd%o%uses_laplacian) then
          delta = abs(lp/100.0_double)
          where (delta >= maxdelta) delta = maxdelta
          pd = lp + delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,dn,pd,exc_pd) ; if (error()) goto 100
          md = lp - delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,dn,md,exc_md) ; if (error()) goto 100
          dfxcdlapl = 0.5_double*(exc_pd*n - exc_md*n)/delta
        end if

100     if (allocated( delta )) deallocate( delta )
        if (associated( exc_md )) deallocate( exc_md )
        if (associated( exc_pd )) deallocate( exc_pd )
        if (associated( md )) deallocate( md )
        if (associated( pd )) deallocate( pd )

        call glean(thy(xcd))

        if (error("Exit xc_density_native_mod::xcd_numerical_i")) continue

      end subroutine

      subroutine xcq_numerical_i(xcd,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(xc_density_native_obj) :: xcd
        real(double), dimension(:,:,:), pointer :: n, dn, lp, exc, dfxcdn, dfxcddnodn, dfxcdlapl

        integer :: n1, n2, n3
        real(double), parameter :: maxdelta = 1.0e-4_double
        real(double), dimension(:,:,:), allocatable :: delta
        real(double), dimension(:,:,:), pointer :: exc_md, exc_pd, md, pd

        call my(xcd)

        nullify( exc_md, exc_pd, md, pd )

        n1 = size(n,1)
        n2 = size(n,2)
        n3 = size(n,3)

        allocate( delta(n1,n2,n3) )
        allocate( exc_md(n1,n2,n3), exc_pd(n1,n2,n3), md(n1,n2,n3), pd(n1,n2,n3) )

        ! exc
        call xce_analytic(xcd%o%etype,xcd%o%ctype,n,dn,lp,exc) ; if (error()) goto 100

        ! derivatives wrt density
        delta = n/100.0_double
        where (delta >= maxdelta) delta = maxdelta
        pd = n + delta
        call xce_analytic(xcd%o%etype,xcd%o%ctype,pd,dn,lp,exc_pd) ; if (error()) goto 100
        md = n - delta
        call xce_analytic(xcd%o%etype,xcd%o%ctype,md,dn,lp,exc_md) ; if (error()) goto 100
        dfxcdn = 0.5_double*(exc_pd*pd - exc_md*md)/delta

        ! derivatives wrt |gradient density|
        if (xcd%o%uses_gradient) then
          delta = dn/100.0_double
          where (delta >= maxdelta) delta = maxdelta
          pd = dn + delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,pd,lp,exc_pd) ; if (error()) goto 100
          md = dn - delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,md,lp,exc_md) ; if (error()) goto 100
          dfxcddnodn = 0.5_double*(exc_pd*n - exc_md*n)/delta/dn
        end if

        ! derivatives wrt laplacian
        if (xcd%o%uses_laplacian) then
          delta = abs(lp/100.0_double)
          where (delta >= maxdelta) delta = maxdelta
          pd = lp + delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,dn,pd,exc_pd) ; if (error()) goto 100
          md = lp - delta
          call xce_analytic(xcd%o%etype,xcd%o%ctype,n,dn,md,exc_md) ; if (error()) goto 100
          dfxcdlapl = 0.5_double*(exc_pd*n - exc_md*n)/delta
        end if

100     if (allocated( delta )) deallocate( delta )
        if (associated( exc_md )) deallocate( exc_md )
        if (associated( exc_pd )) deallocate( exc_pd )
        if (associated( md )) deallocate( md )
        if (associated( pd )) deallocate( pd )

        call glean(thy(xcd))

        if (error("Exit xc_density_native_mod::xcq_numerical_i")) continue

      end subroutine

      subroutine white_bird_i(xcd,lay,nder,dfxcdn,dfxcddnodn,dfxcdlapl,vxc)
        type(xc_density_native_obj) :: xcd
        type(layout_obj) :: lay
        type(dens_der) :: nder
        real(double), dimension(:,:,:), pointer :: dfxcdn, dfxcddnodn, dfxcdlapl, vxc

        real(double), dimension(:,:,:), pointer :: r1
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: gr1

        call my(xcd)
        call my(lay)

        nullify( r1, c1, c2 )

        call my(grid(lay,SGROUP),gr1)

        call alloc(r1,lay,D_TYPE,SGROUP)

        ! contribution due to the density
          vxc = dfxcdn

        ! contribution due to the gradient
        if (xcd%o%uses_gradient) then
          call alloc(c2,lay,D_TYPE,SGROUP)
          r1 = nder%dnx*dfxcddnodn
          call put(r1,gr1,RD_KIND)
          call take(c1,gr1,CDF_KIND)
          c2 = c1*xcd%o%gx
          deallocate( c1 )
          call alloc(r1,lay,D_TYPE,SGROUP)
          r1 = nder%dny*dfxcddnodn
          call put(r1,gr1,RD_KIND)
          call take(c1,gr1,CDF_KIND)
          c2 = c2 + c1*xcd%o%gy
          deallocate( c1 )
          call alloc(r1,lay,D_TYPE,SGROUP)
          r1 = nder%dnz*dfxcddnodn
          call put(r1,gr1,RD_KIND)
          call take(c1,gr1,CDF_KIND)
          c2 = c2 + c1*xcd%o%gz
          deallocate( c1 )
          c2 = (0.0_double,1.0_double)*c2
          call put(c2,gr1,CDF_KIND)
          call take(r1,gr1,RD_KIND)

          vxc = vxc - r1
        end if

        ! contribution due to the laplacian
        if (xcd%o%uses_laplacian) then
          r1 = dfxcdlapl
          call put(r1,gr1,RD_KIND)
          call take(c1,gr1,CDF_KIND)
          c1 = c1*(xcd%o%gx**2 + xcd%o%gy**2 + xcd%o%gz**2)
          call put(c1,gr1,CDF_KIND)
          call take(r1,gr1,RD_KIND)
          vxc = vxc - r1
        end if

100     if (associated( r1 )) deallocate( r1 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(gr1))

        call glean(thy(xcd))
        call glean(thy(lay))

        if (error("Exit xc_density_native_mod::white_bird_i")) continue

      end subroutine

      subroutine form_nder_i(xcd,lay,n,nder,stress)
        type(xc_density_native_obj) :: xcd
        type(layout_obj) :: lay
        real(double), dimension(:,:,:) :: n
        type(dens_der) :: nder
        logical, intent(in), optional :: stress

        logical :: stress_on
        real(double), dimension(:,:,:), pointer :: r1
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: gr1

        call my(xcd)
        call my(lay)

        nullify( nder%n )
        nullify( nder%dn )
        nullify( nder%lapl )
        nullify( nder%dnx, nder%dny, nder%dnz )
        nullify( nder%d2nxx, nder%d2nxy, nder%d2nxz, nder%d2nyy, nder%d2nyz, nder%d2nzz )

        nullify( r1, c1, c2 )

        stress_on = .false.
        if (present(stress)) stress_on = stress

        call my(grid(lay,SGROUP),gr1)

        call alloc(nder%n,lay,D_TYPE,SGROUP)
        nder%n = n
        where (nder%n < mindens) nder%n = mindens

        call alloc(r1,lay,D_TYPE,SGROUP)
        r1 = n
        call put(r1,gr1,RD_KIND)
        call take(c1,gr1,CDF_KIND)

        if (xcd%o%uses_gradient) then
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = (0.0_double,1.0_double)*c1*xcd%o%gx
          call put(c2,gr1,CDF_KIND)
          call take(nder%dnx,gr1,RD_KIND)
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = (0.0_double,1.0_double)*c1*xcd%o%gy
          call put(c2,gr1,CDF_KIND)
          call take(nder%dny,gr1,RD_KIND)
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = (0.0_double,1.0_double)*c1*xcd%o%gz
          call put(c2,gr1,CDF_KIND)
          call take(nder%dnz,gr1,RD_KIND)
          call alloc(nder%dn,lay,D_TYPE,SGROUP)
          where (nder%n <= mindens)
            nder%dnx = mindens
            nder%dny = mindens
            nder%dnz = mindens
          end where
          nder%dn = sqrt(nder%dnx**2 + nder%dny**2 + nder%dnz**2)
        end if

        if (xcd%o%uses_laplacian) then
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = -c1*(xcd%o%gx**2 + xcd%o%gy**2 + xcd%o%gz**2)
          call put(c2,gr1,CDF_KIND)
          call take(nder%lapl,gr1,RD_KIND)
          where (nder%n <= mindens) nder%lapl = mindens
          if (stress_on) then
            call alloc(c2,lay,D_TYPE,SGROUP)
            c2 = -c1*xcd%o%gx**2
            call put(c2,gr1,CDF_KIND)
            call take(nder%d2nxx,gr1,RD_KIND)
            call alloc(c2,lay,D_TYPE,SGROUP)
            c2 = -c1*xcd%o%gx*xcd%o%gy
            call put(c2,gr1,CDF_KIND)
            call take(nder%d2nxy,gr1,RD_KIND)
            call alloc(c2,lay,D_TYPE,SGROUP)
            c2 = -c1*xcd%o%gx*xcd%o%gz
            call put(c2,gr1,CDF_KIND)
            call take(nder%d2nxz,gr1,RD_KIND)
            call alloc(c2,lay,D_TYPE,SGROUP)
            c2 = -c1*xcd%o%gy**2
            call put(c2,gr1,CDF_KIND)
            call take(nder%d2nyy,gr1,RD_KIND)
            call alloc(c2,lay,D_TYPE,SGROUP)
            c2 = -c1*xcd%o%gy*xcd%o%gz
            call put(c2,gr1,CDF_KIND)
            call take(nder%d2nyz,gr1,RD_KIND)
            call alloc(c2,lay,D_TYPE,SGROUP)
            c2 = -c1*xcd%o%gz**2
            call put(c2,gr1,CDF_KIND)
            call take(nder%d2nzz,gr1,RD_KIND)
            where (nder%n <= mindens) 
              nder%d2nxx = mindens
              nder%d2nxy = mindens
              nder%d2nxz = mindens
              nder%d2nyy = mindens
              nder%d2nyz = mindens
              nder%d2nzz = mindens
            end where
          end if
        end if

        if (associated( r1 )) deallocate( r1 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(gr1))

        call glean(thy(xcd))
        call glean(thy(lay))

      end subroutine

      subroutine kill_nder_i(nder)
        type(dens_der) :: nder

        if (associated( nder%n)) deallocate( nder%n )
        if (associated( nder%dn)) deallocate( nder%dn )
        if (associated( nder%dnx)) deallocate( nder%dnx )
        if (associated( nder%dny)) deallocate( nder%dny )
        if (associated( nder%dnz)) deallocate( nder%dnz )
        if (associated( nder%d2nxx)) deallocate( nder%d2nxx )
        if (associated( nder%d2nxy)) deallocate( nder%d2nxy )
        if (associated( nder%d2nxz)) deallocate( nder%d2nxz )
        if (associated( nder%d2nyy)) deallocate( nder%d2nyy )
        if (associated( nder%d2nyz)) deallocate( nder%d2nyz )
        if (associated( nder%d2nzz)) deallocate( nder%d2nzz )

      end subroutine

      subroutine own_i(xcd)
        type(xc_density_native_obj) :: xcd
        type(xc_density_native_obj) :: xcd_t
        if (xcd%ref < xcd%o%ref) then
          allocate( xcd_t%o )
          xcd_t%o%ref               = 0
          xcd_t%o%g                 = xcd%o%g
          xcd_t%o%etype             = xcd%o%etype
          xcd_t%o%ctype             = xcd%o%ctype
          xcd_t%o%derivative_method = xcd%o%derivative_method
          xcd_t%o%potential_method  = xcd%o%potential_method
          xcd_t%o%uses_gradient     = xcd%o%uses_gradient
          xcd_t%o%uses_laplacian    = xcd%o%uses_laplacian
          xcd_t%o%volume            = xcd%o%volume
          call my(xcd%o%xct,xcd_t%o%xct)
          call my(xcd%o%lay,xcd_t%o%lay)
          call my(xcd%o%sg,xcd_t%o%sg)
          if (xcd%o%uses_gradient .or. xcd%o%uses_laplacian) then
            allocate( xcd_t%o%gx(size(xcd%o%gx,1),size(xcd%o%gx,2),size(xcd%o%gx,3)) )
              xcd_t%o%gx = xcd%o%gx
            allocate( xcd_t%o%gy(size(xcd%o%gy,1),size(xcd%o%gy,2),size(xcd%o%gy,3)) )
              xcd_t%o%gy = xcd%o%gy
            allocate( xcd_t%o%gz(size(xcd%o%gz,1),size(xcd%o%gz,2),size(xcd%o%gz,3)) )
              xcd_t%o%gz = xcd%o%gz
          else
            nullify( xcd_t%o%gx )
            nullify( xcd_t%o%gy )
            nullify( xcd_t%o%gz )
          end if
          xcd%o%ref = xcd%o%ref - xcd%ref
          xcd%o => xcd_t%o
          xcd%o%ref = xcd%o%ref + xcd%ref
        end if
      end subroutine

      end module

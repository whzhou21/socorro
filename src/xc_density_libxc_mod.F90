!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module xc_density_libxc_mod
!doc$ module xc_density_libxc_mod

!     One datatype is available here: type(xc_density_libxc_obj).

!     xc_density_libxc_mod encapsulates libxc routines used to compute density-dependent
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
      use xc_f90_types_m
      use libxc_funcs_m
      use xc_f90_lib_m

!cod$
      implicit none
      private 

      real(double), parameter :: mindens = 1.0e-20_double
      real(double), parameter :: Hartrees2Rydbergs = 2.0_double

      type :: xc_density_libxc_rep
        integer :: ref                                        ! reference count
        type(ghost) :: g                                      ! ghost
        integer :: xc_id                                      ! exchange-correlation type
        integer :: x_id                                       ! exchange type
        integer :: c_id                                       ! correlation type
        integer :: potential_method                           ! method for computing the potential
        logical :: uses_gradient                              ! indicates whether or not the gradient is used
        logical :: uses_laplacian                             ! indicates whether or not the laplacian is used
        integer :: spin_status                                ! spin status indicator
        real(double) :: volume                                ! cell volume
        type(xc_type_obj) :: xct                              ! xc_type object
        type(layout_obj) :: lay                               ! layout object
        type(space_group_obj) :: sg                           ! space group object
        real(double), dimension(:,:,:), pointer :: gx         ! x component of reciprocal-space mesh values
        real(double), dimension(:,:,:), pointer :: gy         ! y component of reciprocal-space mesh values
        real(double), dimension(:,:,:), pointer :: gz         ! z component of reciprocal-space mesh values
      end type

      type, public :: xc_density_libxc_obj
        private
        integer :: ref
        type(xc_density_libxc_rep), pointer :: o
      end type

      type :: dens_der
        real(double), dimension(:,:,:,:), pointer :: n
        real(double), dimension(:,:,:,:), pointer :: dnx
        real(double), dimension(:,:,:,:), pointer :: dny
        real(double), dimension(:,:,:,:), pointer :: dnz
        real(double), dimension(:,:,:,:), pointer :: sigma
      end type

      public :: xc_density_libxc
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

      interface xc_density_libxc
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
        module procedure xcd_energy_libxc
      end interface
      interface xc_potential
        module procedure xcd_potential_libxc
      end interface
      interface xc_energy_and_potential
        module procedure xcd_energy_and_potential_libxc
      end interface
      interface xc_grid_pressure
        module procedure xcd_grid_pressure_libxc
      end interface
      interface xc_grid_stress_tensor
        module procedure xcd_grid_stress_tensor_libxc
      end interface

      contains

! public routines

      function constructor_xcd(xct,lay,sg) result(xcd)
!doc$ function xc_density_libxc(xct,lay,sg) result(xcd)
        type(xc_type_obj) :: xct
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_density_libxc_obj) :: xcd
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

        xcd%o%xc_id = x_functional_type(xcd%o%xct)
        xcd%o%x_id = x_exchange_type(xcd%o%xct)
        xcd%o%c_id = x_correlation_type(xcd%o%xct)
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

        select case (mpi_nsgroups())
        case (1)
          xcd%o%spin_status = XC_UNPOLARIZED
        case (2)
          xcd%o%spin_status = XC_POLARIZED
        end select

100     call glean(thy(xct))
        call glean(thy(lay))
        call glean(thy(sg))

        if (error("Exit xc_density_libxc_mod::constructor_xcd")) continue

      end function 

      subroutine update_xcd(xcd,lay,sg)
!doc$ subroutine update(xcd,lay,sg)
        type(xc_density_libxc_obj) :: xcd
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
        type(xc_density_libxc_obj) :: xcd

!cod$
        xcd%ref = xcd%ref + 1
        xcd%o%ref = xcd%o%ref + 1
      end subroutine

      subroutine my_new_xcd(xcdi,xcd)
!doc$ subroutine my(xcdi,xcd)
        type(xc_density_libxc_obj) :: xcdi, xcd

!cod$
        xcd%ref = 1
        xcd%o => xcdi%o
        xcd%o%ref = xcd%o%ref + 1

      end subroutine
      
      function thy_xcd(xcd) result(xcdo)
!doc$ function thy(xcd) result(xcdo)
        type(xc_density_libxc_obj) :: xcd, xcdo

!cod$
        xcd%ref = xcd%ref - 1
        xcd%o%ref = xcd%o%ref - 1
        xcdo%ref = xcd%ref
        xcdo%o => xcd%o
      end function

      subroutine glean_xcd(xcd)
!doc$ subroutine glean(xcd)
        type(xc_density_libxc_obj) :: xcd

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
        type(xc_density_libxc_obj) :: xcd

!cod$
        continue
      end subroutine

      subroutine assign_xcd(xcd,xcd2)
!doc$ subroutine assign(xcd,xcd2)
        type(xc_density_libxc_obj), intent(inout) :: xcd
        type(xc_density_libxc_obj), intent(in) :: xcd2

!cod$
        type(xc_density_libxc_obj) :: xcdt
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
        type(xc_density_libxc_obj) :: xcd
        integer, dimension(2) :: r
!       effects: Returns xcd%ref and xcd%o%ref.

!cod$
        r(1) = xcd%ref
        r(2) = xcd%o%ref
        call glean(xcd)
      end function

      function xcd_ghost(xcd) result(g)
!doc$ function x_ghost(xcd) result(g)
        type(xc_density_libxc_obj) :: xcd
        type(ghost) :: g
!       effects: Returns the ghost of xcd.

!cod$
        call my(xcd)
        g = xcd%o%g
        call glean(thy(xcd))
      end function

      function xcd_xc_type(xcd) result(xct)
!doc$ function x_xc_type(xcd) result(xct)
        type(xc_density_libxc_obj) :: xcd
        type(xc_type_obj)  :: xct 
!       effects: Returns the xc_type of xcd.

!cod$
        call my(xcd)
        call my(xcd%o%xct,xct)
        call glean(thy(xcd))
        call bequeath(thy(xct))
      end function

      function xcd_layout(xcd) result(lay)
!doc$ function x_layout(xcd) result(lay)
        type(xc_density_libxc_obj) :: xcd
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
        type(xc_density_libxc_obj) :: xcd
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
        type(xc_density_libxc_obj) :: xcd
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        call my(xcd)
        ug = xcd%o%uses_gradient 
        call glean(thy(xcd))
      end function

      function q_uses_laplacian(xcd) result(ul)
!doc$ function uses_laplacian(xcd) result(ul)
        type(xc_density_libxc_obj) :: xcd
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        call my(xcd)
        ul = xcd%o%uses_laplacian
        call glean(thy(xcd))
      end function

      function xcd_energy_libxc(xcd,n_total,n_valence) result(fxc)
!doc$ function xc_energy(xcd,n_total,n_valence) result(fxc)
        type(xc_density_libxc_obj) :: xcd
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        real(double) :: fxc
!       requires: n_total and n_valence data be symmetrized.
!       effects: Returns the exchange-correlation energy
!       errors: Passes errors.

!cod$
        integer :: np
        real(double) :: e_local, e_global, semilocal_fraction, exx_fraction
        real(double), dimension(:,:,:), pointer :: ntot, exc, ex, ec
        real(double), dimension(:,:,:), pointer :: nval, ex_val
        type(layout_obj) :: lay
        type(dens_der) :: nder_total
        type(dens_der) :: nder_valence
        type(xc_f90_pointer_t) :: xc_func, xc_info

        call my(xcd)
        call my(n_total)
        call my(n_valence)

        call my(x_layout(n_total),lay)

        nullify( ntot )
        nullify( nval )
        nullify( exc )
        nullify( ex )
        nullify( ec )
        nullify( ex_val )

        call take(ntot,n_total,RD_KIND)
        call take(nval,n_valence,RD_KIND)

        call alloc(exc,lay,D_TYPE,SGROUP)
        exc = 0.0_double
        
        if (uses_nlcc(xcd%o%xct)) then
          call alloc(ex_val,lay,D_TYPE,SGROUP)
          ex_val = 0.0_double
        endif

        call form_nder_i(xcd,lay,ntot,nder_total)
        if (uses_nlcc(xcd%o%xct)) then
          call form_nder_i(xcd,lay,nval,nder_valence)
        endif

        np = size(ntot,1)*size(ntot,2)*size(ntot,3)

        exx_fraction = x_hybrid_mixing(xcd%o%xct)
        semilocal_fraction = 1.0_double - exx_fraction

        ! contribution from xc
        select case(xc_f90_family_from_id(xcd%o%xc_id))
        case (XC_FAMILY_LDA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_lda_exc(xc_func,np,nder_total%n(1,1,1,1),exc(1,1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_GGA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),exc(1,1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_HYB_GGA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),exc(1,1,1))
          call xc_f90_func_end(xc_func)
        end select

        ! contribution from x
        select case(xc_f90_family_from_id(xcd%o%x_id))
        case (XC_FAMILY_LDA)
          call alloc(ex,lay,D_TYPE,SGROUP)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_lda_exc(xc_func,np,nder_total%n(1,1,1,1),ex(1,1,1))
          call xc_f90_func_end(xc_func)
          select case (x_functional_dependence(xcd%o%xct))
          case (FD_DENSITY)
            exc = exc + ex
          case (FD_HYBRID)
            ! In the case of a hybrid, there is an option to include a compensation 
            ! term to account for the nonlinear core correction of the exx
            ! fraction of exchange (alpha)
            !   
            !   Ex = alpha*(Ex_exx - Ex_semiloc(valence)) + Ex_semilocal(total)
            !
            if (uses_nlcc(xcd%o%xct)) then
              exc = exc + ex
              call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
              call xc_f90_lda_exc(xc_func,np,nder_valence%n(1,1,1,1),ex_val(1,1,1))
              call xc_f90_func_end(xc_func)
            else
              exc = exc + ex*semilocal_fraction
            endif
          end select
          deallocate( ex )
        case (XC_FAMILY_GGA)
          call alloc(ex,lay,D_TYPE,SGROUP)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),ex(1,1,1))
          call xc_f90_func_end(xc_func)
          select case (x_functional_dependence(xcd%o%xct))
          case (FD_DENSITY)
            exc = exc + ex
          case (FD_HYBRID)
            if (uses_nlcc(xcd%o%xct)) then
              exc = exc + ex
              call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
              call xc_f90_gga_exc(xc_func,np,nder_valence%n(1,1,1,1),nder_valence%sigma(1,1,1,1),ex_val(1,1,1))
              call xc_f90_func_end(xc_func)
            else
              exc = exc + ex*semilocal_fraction
            endif
          end select
          deallocate( ex )
        case (XC_FAMILY_HYB_GGA)
          call alloc(ex,lay,D_TYPE,SGROUP)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),ex(1,1,1))
          call xc_f90_func_end(xc_func)
          if (uses_nlcc(xcd%o%xct)) then
            call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
            call xc_f90_gga_exc(xc_func,np,nder_valence%n(1,1,1,1),nder_valence%sigma(1,1,1,1),ex_val(1,1,1))
            call xc_f90_func_end(xc_func)
            exc = exc + ex
          else
            exc = exc + ex
          endif
          deallocate( ex )
        end select
           
        ! contribution from c
        select case(xc_f90_family_from_id(xcd%o%c_id))
        case (XC_FAMILY_LDA)
          call alloc(ec,lay,D_TYPE,SGROUP)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
          call xc_f90_lda_exc(xc_func,np,nder_total%n(1,1,1,1),ec(1,1,1))
          call xc_f90_func_end(xc_func)
          exc = exc + ec
          deallocate( ec )
        case (XC_FAMILY_GGA)
          call alloc(ec,lay,D_TYPE,SGROUP)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),ec(1,1,1))
          call xc_f90_func_end(xc_func)
          exc = exc + ec
          deallocate( ec )
        end select

        exc = exc*Hartrees2Rydbergs
        if (uses_nlcc(xcd%o%xct)) then
          ex_val = ex_val*Hartrees2Rydbergs
        endif

        call kill_nder_i(nder_total)
        if (uses_nlcc(xcd%o%xct)) then
          call kill_nder_i(nder_valence)
        endif

        e_local = sum(ntot*exc)
        if (uses_nlcc(xcd%o%xct)) then
          e_local = e_local - exx_fraction*sum(nval*ex_val)
        endif

        call allreduce(CONFIG,MPI_SUM,e_local,e_global)
        fxc = e_global*(xcd%o%volume/product(x_dims(lay)))

        call put(ntot,n_total,RD_KIND)
        call put(nval,n_valence,RD_KIND)

100     if (associated( ntot )) deallocate( ntot )
        if (associated( nval )) deallocate( nval )
        if (associated( exc )) deallocate( exc )
        if (associated( ex ))  deallocate( ex )
        if (associated( ec ))  deallocate( ec )
        if (associated( ex_val )) deallocate( ex_val )

        call glean(thy(lay))

        call glean(thy(xcd))
        call glean(thy(n_total))
        call glean(thy(n_valence))

        if (error("Exit xc_density_libxc_mod::xcd_energy_libxc")) continue

      end function

      function xcd_potential_libxc(xcd,n_total,n_valence) result(vxc_g)
!doc$ function xc_potential(xcd,n_total,n_valence) result(vxc_g)
        type(xc_density_libxc_obj) :: xcd
        type(grid_obj) :: n_total
        type(grid_obj) :: n_valence
        type(grid_obj) :: vxc_g
!       requires: n_total and n_valence data be symmetrized with respect to xcd%o%sg.
!       effects: Returns the symmetrized/filtered exchange-correlation potential.
!       errors: Passes errors.

!cod$ 
        integer :: nx, ny, nz, np
        real(double) :: semilocal_fraction, exx_fraction
        real(double), dimension(:,:,:),   pointer :: ntot, nval, vxc_sg
        real(double), dimension(:,:,:,:), pointer :: vxc
        real(double), dimension(:,:,:,:), pointer :: dfxcdn, dfxdn, dfcdn
        real(double), dimension(:,:,:,:), pointer :: dfxcdsigma, dfxdsigma, dfcdsigma
        real(double), dimension(:,:,:,:), pointer :: vxc_val, dfxdn_val, dfxdsigma_val 
        type(layout_obj) :: lay
        type(dens_der) :: nder_total, nder_valence
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        call my(xcd)
        call my(n_total) 
        call my(n_valence) 

        nullify( ntot )
        nullify( nval )
        nullify( vxc )
        nullify( vxc_sg )
        nullify( dfxcdn )
        nullify( dfxdn )
        nullify( dfcdn )
        nullify( dfxcdsigma )
        nullify( dfxdsigma )
        nullify( dfcdsigma )
        nullify( vxc_val )
        nullify( dfxdn_val )
        nullify( dfxdsigma_val )

        call my(x_layout(n_total),lay)
        call my(grid(lay,SGROUP),vxc_g)

        call take(ntot,n_total,RD_KIND)
        call take(nval,n_valence,RD_KIND)

        nx = size(ntot,1)
        ny = size(ntot,2)
        nz = size(ntot,3)

        np = nx*ny*nz

        select case (xcd%o%spin_status)
        case (XC_UNPOLARIZED)
          allocate( dfxcdn(1,nx,ny,nz) )
          if (xcd%o%uses_gradient) allocate( dfxcdsigma(1,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) then
             allocate( dfxdn_val(1,nx,ny,nz) )
             if (xcd%o%uses_gradient) allocate( dfxdsigma_val(1,nx,ny,nz) )
          endif
        case (XC_POLARIZED)
          allocate( dfxcdn(2,nx,ny,nz) )
          if (xcd%o%uses_gradient) allocate( dfxcdsigma(3,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) then
            allocate( dfxdn_val(2,nx,ny,nz) )
            if (xcd%o%uses_gradient) allocate( dfxdsigma_val(3,nx,ny,nz) )
          endif
        end select
        dfxcdn = 0.0_double
        if (xcd%o%uses_gradient) then
          dfxcdsigma = 0.0_double
        end if

        if (uses_nlcc(xcd%o%xct)) then
          dfxdn_val = 0.0_double
          if (xcd%o%uses_gradient) dfxdsigma_val = 0.0_double
        endif

        call form_nder_i(xcd,lay,ntot,nder_total)
        if (uses_nlcc(xcd%o%xct)) call form_nder_i(xcd,lay,nval,nder_valence)

        exx_fraction = x_hybrid_mixing(xcd%o%xct)
        semilocal_fraction = 1.0_double - exx_fraction

        ! contribution from xc
        select case(xc_f90_family_from_id(xcd%o%xc_id))
        case (XC_FAMILY_LDA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_lda_vxc(xc_func,np,nder_total%n(1,1,1,1),dfxcdn(1,1,1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_GGA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),dfxcdn(1,1,1,1),dfxcdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_HYB_GGA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),dfxcdn(1,1,1,1),dfxcdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
        end select
 
        ! contribution from x
        select case(xc_f90_family_from_id(xcd%o%x_id))
        case (XC_FAMILY_LDA)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfxdn(2,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_lda_vxc(xc_func,np,nder_total%n(1,1,1,1),dfxdn(1,1,1,1))
          call xc_f90_func_end(xc_func)
          select case (x_functional_dependence(xcd%o%xct))
          case (FD_DENSITY)
            dfxcdn = dfxcdn + dfxdn
          case (FD_HYBRID)
            if (uses_nlcc(xcd%o%xct)) then
              dfxcdn = dfxcdn + dfxdn
              call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
              call xc_f90_lda_vxc(xc_func,np,nder_valence%n(1,1,1,1),dfxdn_val(1,1,1,1))
              call xc_f90_func_end(xc_func)
            else
              dfxcdn = dfxcdn + dfxdn*semilocal_fraction
            endif
          end select
          deallocate( dfxdn )
        case (XC_FAMILY_GGA)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn(1,nx,ny,nz) )
            allocate( dfxdsigma(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfxdn(2,nx,ny,nz) )
            allocate( dfxdsigma(3,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),dfxdn(1,1,1,1),dfxdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
          select case (x_functional_dependence(xcd%o%xct))
          case (FD_DENSITY)
            dfxcdn = dfxcdn + dfxdn
            dfxcdsigma = dfxcdsigma + dfxdsigma
          case (FD_HYBRID)
            if (uses_nlcc(xcd%o%xct)) then
              dfxcdn = dfxcdn + dfxdn
              dfxcdsigma = dfxcdsigma + dfxdsigma
              call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
              call xc_f90_gga_vxc(xc_func,np,nder_valence%n(1,1,1,1),nder_valence%sigma(1,1,1,1),&
                                  dfxdn_val(1,1,1,1),dfxdsigma_val(1,1,1,1))
              call xc_f90_func_end(xc_func)
            else
              dfxcdn = dfxcdn + dfxdn*semilocal_fraction
              dfxcdsigma = dfxcdsigma + dfxdsigma*semilocal_fraction
            endif
          end select
          deallocate( dfxdn )
          deallocate( dfxdsigma )
        case (XC_FAMILY_HYB_GGA)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn(1,nx,ny,nz) )
            allocate( dfxdsigma(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfxdn(2,nx,ny,nz) )
            allocate( dfxdsigma(3,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),dfxdn(1,1,1,1),dfxdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn = dfxcdn + dfxdn
          dfxcdsigma = dfxcdsigma + dfxdsigma
          deallocate( dfxdn )
          deallocate( dfxdsigma )
        end select
        
        ! contribution from c
        select case(xc_f90_family_from_id(xcd%o%c_id))
        case (XC_FAMILY_LDA)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfcdn(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfcdn(2,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
          call xc_f90_lda_vxc(xc_func,np,nder_total%n(1,1,1,1),dfcdn(1,1,1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn = dfxcdn + dfcdn
          deallocate( dfcdn )
        case (XC_FAMILY_GGA)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfcdn(1,nx,ny,nz) )
            allocate( dfcdsigma(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfcdn(2,nx,ny,nz) )
            allocate( dfcdsigma(3,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder_total%n(1,1,1,1),nder_total%sigma(1,1,1,1),dfcdn(1,1,1,1),dfcdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn = dfxcdn + dfcdn
          dfxcdsigma = dfxcdsigma + dfcdsigma
          deallocate( dfcdn )
          deallocate( dfcdsigma )
        end select

        dfxcdn = dfxcdn*Hartrees2Rydbergs
        if (xcd%o%uses_gradient) dfxcdsigma = dfxcdsigma*Hartrees2Rydbergs

        if (uses_nlcc(xcd%o%xct)) then
          dfxdn_val = dfxdn_val*Hartrees2Rydbergs
          if (xcd%o%uses_gradient) dfxdsigma_val = dfxdsigma_val*Hartrees2Rydbergs
        endif

        select case (xcd%o%spin_status)
        case (XC_UNPOLARIZED)
          allocate( vxc(1,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) allocate( vxc_val(1,nx,ny,nz) )
        case (XC_POLARIZED)
          allocate( vxc(2,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) allocate( vxc_val(2,nx,ny,nz) )
        end select

        select case (xcd%o%potential_method)
        case (XCD_SIMPLE)
          vxc = dfxcdn 
          if (uses_nlcc(xcd%o%xct)) vxc = vxc - exx_fraction*dfxdn_val
        case (XCD_WHITE_BIRD)
          call white_bird_i(xcd,lay,nder_total,dfxcdn,dfxcdsigma,vxc)
          if (uses_nlcc(xcd%o%xct)) then
             call white_bird_i(xcd,lay,nder_valence,dfxdn_val,dfxdsigma_val,vxc_val)
             vxc = vxc - exx_fraction*vxc_val    
          endif
        end select
        call kill_nder_i(nder_total)
        if (uses_nlcc(xcd%o%xct)) call kill_nder_i(nder_valence)

        call alloc(vxc_sg,lay,D_TYPE,SGROUP)
        select case (xcd%o%spin_status)
        case (XC_UNPOLARIZED)
          vxc_sg(:,:,:) = vxc(1,:,:,:)
        case (XC_POLARIZED)
          select case (mpi_mysgroup())
          case (1)
            vxc_sg(:,:,:) = vxc(1,:,:,:)
          case (2)
            vxc_sg(:,:,:) = vxc(2,:,:,:)
          end select
        end select

        call put(vxc_sg,vxc_g,RD_KIND)
        call symmetrize_grid(xcd%o%sg,vxc_g)
        call filter(vxc_g)

        call put(ntot,n_total,RD_KIND)
        call put(nval,n_valence,RD_KIND)

100     if (associated( ntot ))       deallocate( ntot )
        if (associated( nval ))       deallocate( nval )
        if (associated( vxc ))        deallocate( vxc )
        if (associated( vxc_sg ))     deallocate( vxc_sg )
        if (associated( dfxcdn ))     deallocate( dfxcdn )
        if (associated( dfxdn ))      deallocate( dfxdn )
        if (associated( dfcdn ))      deallocate( dfcdn )
        if (associated( dfxcdsigma )) deallocate( dfxcdsigma )
        if (associated( dfxdsigma ))  deallocate( dfxdsigma )
        if (associated( dfcdsigma ))  deallocate( dfcdsigma )
        if (associated( vxc_val ))    deallocate( vxc_val )
        if (associated( dfxdn_val ))  deallocate( dfxdn_val )
        if (associated( dfxdsigma_val )) deallocate( dfxdsigma_val )

        call glean(thy(lay))
        call bequeath(thy(vxc_g))
      
        call glean(thy(xcd))
        call glean(thy(n_total))
        call glean(thy(n_valence))

        if (error("Exit xc_density_libxc_mod::xcd_potential_libxc")) continue

      end function

      function xcd_energy_and_potential_libxc(xcd,n_total,n_valence,fxc) result(vxc_g)
!doc$ function xc_energy_and_potential(xcd,n_total,n_valence,fxc) result(vxc_g)
        type(xc_density_libxc_obj) :: xcd
        type(grid_obj)             :: n_total
        type(grid_obj)             :: n_valence
        real(double), intent(out)  :: fxc
        type(grid_obj)             :: vxc_g
!       requires: n_total and n_valence data be symmetrized.
!       effects: Returns the symmetrized/filtered exchange-correlation potential and energy.
!       errors: Passes errors.

!cod$
        integer :: nx, ny, nz, np
        real(double) :: e_local, e_global, semilocal_fraction, exx_fraction
        real(double), dimension(:,:,:), pointer   :: ntot, vxc_sg, exc, ex, ec
        real(double), dimension(:,:,:), pointer   :: nval, ex_val
        real(double), dimension(:,:,:,:), pointer :: vxc
        real(double), dimension(:,:,:,:), pointer :: dfxcdn, dfxdn, dfcdn
        real(double), dimension(:,:,:,:), pointer :: dfxcdsigma, dfxdsigma, dfcdsigma
        real(double), dimension(:,:,:,:), pointer :: vxc_val, dfxdn_val, dfxdsigma_val
        type(layout_obj) :: lay
        type(dens_der) :: nder_total
        type(dens_der) :: nder_valence
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        call my(xcd)
        call my(n_total)
        call my(n_valence)

        call my(x_layout(n_total),lay)
        call my(grid(lay,SGROUP),vxc_g)

        nullify( ntot )
        nullify( nval )
        nullify( exc )
        nullify( ex )
        nullify( ec )
        nullify( vxc )
        nullify( vxc_sg )
        nullify( dfxcdn )
        nullify( dfxdn )
        nullify( dfcdn )
        nullify( dfxcdsigma )
        nullify( dfxdsigma )
        nullify( dfcdsigma )
        nullify( ex_val )
        nullify( vxc_val )
        nullify( dfxdn_val )
        nullify( dfxdsigma_val )

        call take(ntot,n_total,RD_KIND)
        call take(nval,n_valence,RD_KIND)

        nx = size(ntot,1)
        ny = size(ntot,2)
        nz = size(ntot,3)

        np = nx*ny*nz

        call alloc(exc,lay,D_TYPE,SGROUP)
        exc = 0.0_double

        if (uses_nlcc(xcd%o%xct)) then
          call alloc(ex_val,lay,D_TYPE,SGROUP)
          ex_val = 0.0_double
        endif 

        select case (xcd%o%spin_status)
        case (XC_UNPOLARIZED)
          allocate( dfxcdn(1,nx,ny,nz) )
          if (xcd%o%uses_gradient) allocate( dfxcdsigma(1,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) then
             allocate( dfxdn_val(1,nx,ny,nz) )
             if (xcd%o%uses_gradient) allocate( dfxdsigma_val(1,nx,ny,nz) )
          endif
        case (XC_POLARIZED)
          allocate( dfxcdn(2,nx,ny,nz) )
          if (xcd%o%uses_gradient) allocate( dfxcdsigma(3,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) then
            allocate( dfxdn_val(2,nx,ny,nz) )
            if (xcd%o%uses_gradient) allocate( dfxdsigma_val(3,nx,ny,nz) )
          endif
        end select
        dfxcdn = 0.0_double
        if (xcd%o%uses_gradient) then
          dfxcdsigma = 0.0_double
        end if

        if (uses_nlcc(xcd%o%xct)) then
          dfxdn_val = 0.0_double
          if (xcd%o%uses_gradient) dfxdsigma_val = 0.0_double
        endif

        call form_nder_i(xcd,lay,ntot,nder_total)
        if (uses_nlcc(xcd%o%xct)) call form_nder_i(xcd,lay,nval,nder_valence)

        exx_fraction = x_hybrid_mixing(xcd%o%xct)
        semilocal_fraction = 1.0_double - exx_fraction

        ! contribution from xc
        select case(xc_f90_family_from_id(xcd%o%xc_id))
        case (XC_FAMILY_LDA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_lda_exc_vxc(xc_func,np,nder_total%n(1,1,1,1),exc(1,1,1),dfxcdn(1,1,1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_GGA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_gga_exc_vxc(xc_func,np,nder_total%n(1,1,1,1), &
                                  nder_total%sigma(1,1,1,1),exc(1,1,1),dfxcdn(1,1,1,1),dfxcdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_HYB_GGA)
          call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
          call xc_f90_gga_exc_vxc(xc_func,np,nder_total%n(1,1,1,1), &
                                  nder_total%sigma(1,1,1,1),exc(1,1,1),dfxcdn(1,1,1,1),dfxcdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
        end select

        ! contribution from x
        select case(xc_f90_family_from_id(xcd%o%x_id))
        case (XC_FAMILY_LDA)
          call alloc(ex,lay,D_TYPE,SGROUP)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfxdn(2,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_lda_exc_vxc(xc_func,np,nder_total%n(1,1,1,1),ex(1,1,1),dfxdn(1,1,1,1))
          call xc_f90_func_end(xc_func)
          select case (x_functional_dependence(xcd%o%xct))
          case (FD_DENSITY)
            exc = exc + ex
            dfxcdn = dfxcdn + dfxdn
          case (FD_HYBRID)
            if (uses_nlcc(xcd%o%xct)) then
              exc = exc + ex
              dfxcdn = dfxcdn + dfxdn
              call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
              call xc_f90_lda_exc_vxc(xc_func,np,nder_valence%n(1,1,1,1),ex_val(1,1,1),dfxdn_val(1,1,1,1))
              call xc_f90_func_end(xc_func)
            else
              exc = exc + ex*semilocal_fraction
              dfxcdn = dfxcdn + dfxdn*semilocal_fraction
            endif
          end select
          deallocate( ex )
          deallocate( dfxdn )
        case (XC_FAMILY_GGA)
          call alloc(ex,lay,D_TYPE,SGROUP)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn(1,nx,ny,nz) )
            allocate( dfxdsigma(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfxdn(2,nx,ny,nz) )
            allocate( dfxdsigma(3,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_gga_exc_vxc(xc_func,np,nder_total%n(1,1,1,1), &
                                  nder_total%sigma(1,1,1,1),ex(1,1,1),dfxdn(1,1,1,1),dfxdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
          select case (x_functional_dependence(xcd%o%xct))
          case (FD_DENSITY)
            exc = exc + ex
            dfxcdn = dfxcdn + dfxdn
            dfxcdsigma = dfxcdsigma + dfxdsigma
          case (FD_HYBRID)
            if (uses_nlcc(xcd%o%xct)) then
              exc = exc + ex
              dfxcdn = dfxcdn + dfxdn
              dfxcdsigma = dfxcdsigma + dfxdsigma
              call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
              call xc_f90_gga_exc_vxc(xc_func,np,nder_valence%n(1,1,1,1), &
                                      nder_valence%sigma(1,1,1,1),ex_val(1,1,1),&
                                      dfxdn_val(1,1,1,1),dfxdsigma_val(1,1,1,1))
              call xc_f90_func_end(xc_func)
            else
              exc = exc + ex*semilocal_fraction
              dfxcdn = dfxcdn + dfxdn*semilocal_fraction
              dfxcdsigma = dfxcdsigma + dfxdsigma*semilocal_fraction
            endif
          end select
          deallocate( ex )
          deallocate( dfxdn )
          deallocate( dfxdsigma )
        case (XC_FAMILY_HYB_GGA)
          call alloc(ex,lay,D_TYPE,SGROUP)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn(1,nx,ny,nz) )
            allocate( dfxdsigma(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfxdn(2,nx,ny,nz) )
            allocate( dfxdsigma(3,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
          call xc_f90_gga_exc_vxc(xc_func,np,nder_total%n(1,1,1,1), &
                                  nder_total%sigma(1,1,1,1),ex(1,1,1),dfxdn(1,1,1,1),dfxdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
          exc = exc + ex
          dfxcdn = dfxcdn + dfxdn
          dfxcdsigma = dfxcdsigma + dfxdsigma
          deallocate( ex )
          deallocate( dfxdn )
          deallocate( dfxdsigma )
        end select
 
        ! contribution from c
        select case(xc_f90_family_from_id(xcd%o%c_id))
        case (XC_FAMILY_LDA)
          call alloc(ec,lay,D_TYPE,SGROUP)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfcdn(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfcdn(2,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
          call xc_f90_lda_exc_vxc(xc_func,np,nder_total%n(1,1,1,1),ec(1,1,1),dfcdn(1,1,1,1))
          call xc_f90_func_end(xc_func)
          exc = exc + ec
          dfxcdn = dfxcdn + dfcdn
          deallocate( ec )
          deallocate( dfcdn )
        case (XC_FAMILY_GGA)
          call alloc(ec,lay,D_TYPE,SGROUP)
          select case (xcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfcdn(1,nx,ny,nz) )
            allocate( dfcdsigma(1,nx,ny,nz) )
          case (XC_POLARIZED)
            allocate( dfcdn(2,nx,ny,nz) )
            allocate( dfcdsigma(3,nx,ny,nz) )
          end select
          call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
          call xc_f90_gga_exc_vxc(xc_func,np,nder_total%n(1,1,1,1), &
                                  nder_total%sigma(1,1,1,1),ec(1,1,1),dfcdn(1,1,1,1),dfcdsigma(1,1,1,1))
          call xc_f90_func_end(xc_func)
          exc = exc + ec
          dfxcdn = dfxcdn + dfcdn
          dfxcdsigma = dfxcdsigma + dfcdsigma
          deallocate( ec )
          deallocate( dfcdn )
          deallocate( dfcdsigma )
        end select

        exc = exc*Hartrees2Rydbergs
        dfxcdn = dfxcdn*Hartrees2Rydbergs
        if (xcd%o%uses_gradient) dfxcdsigma = dfxcdsigma*Hartrees2Rydbergs

        if (uses_nlcc(xcd%o%xct)) then
          ex_val = ex_val*Hartrees2Rydbergs
          dfxdn_val = dfxdn_val*Hartrees2Rydbergs
          if (xcd%o%uses_gradient) dfxdsigma_val = dfxdsigma_val*Hartrees2Rydbergs
        endif

        select case (xcd%o%spin_status)
        case (XC_UNPOLARIZED)
          allocate( vxc(1,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) allocate( vxc_val(1,nx,ny,nz) )
        case (XC_POLARIZED)
          allocate( vxc(2,nx,ny,nz) )
          if (uses_nlcc(xcd%o%xct)) allocate( vxc_val(2,nx,ny,nz) )
        end select

        select case (xcd%o%potential_method)
        case (XCD_SIMPLE)
          vxc = dfxcdn 
          if (uses_nlcc(xcd%o%xct)) vxc = vxc - exx_fraction*dfxdn_val
        case (XCD_WHITE_BIRD)
          call white_bird_i(xcd,lay,nder_total,dfxcdn,dfxcdsigma,vxc)
          if (uses_nlcc(xcd%o%xct)) then
             call white_bird_i(xcd,lay,nder_valence,dfxdn_val,dfxdsigma_val,vxc_val)
             vxc = vxc - exx_fraction*vxc_val
          endif
        end select
        call kill_nder_i(nder_total)
        if (uses_nlcc(xcd%o%xct)) call kill_nder_i(nder_valence)


        call alloc(vxc_sg,lay,D_TYPE,SGROUP)
        select case (xcd%o%spin_status)
        case (XC_UNPOLARIZED)
          vxc_sg(:,:,:) = vxc(1,:,:,:)
        case (XC_POLARIZED)
          select case (mpi_mysgroup())
          case (1)
            vxc_sg(:,:,:) = vxc(1,:,:,:)
          case (2)
            vxc_sg(:,:,:) = vxc(2,:,:,:)
          end select
        end select

        call put(vxc_sg,vxc_g,RD_KIND)
        call symmetrize_grid(xcd%o%sg,vxc_g)
        call filter(vxc_g)

        e_local = sum(ntot*exc)
        if (uses_nlcc(xcd%o%xct)) then
           e_local = e_local - exx_fraction*sum(nval*ex_val)
        endif
        call allreduce(CONFIG,MPI_SUM,e_local,e_global)
        fxc = e_global*(xcd%o%volume/product(x_dims(lay)))

        call put(ntot,n_total,RD_KIND)
        call put(nval,n_valence,RD_KIND)

100     if (associated( ntot )) deallocate( ntot )
        if (associated( nval )) deallocate( nval )
        if (associated( exc )) deallocate( exc )
        if (associated( ex )) deallocate( ex )
        if (associated( ec )) deallocate( ec )
        if (associated( vxc )) deallocate( vxc )
        if (associated( vxc_sg )) deallocate( vxc_sg )
        if (associated( dfxcdn )) deallocate( dfxcdn )
        if (associated( dfxdn )) deallocate( dfxdn )
        if (associated( dfcdn )) deallocate( dfcdn )
        if (associated( dfxcdsigma )) deallocate( dfxcdsigma )
        if (associated( dfxdsigma )) deallocate( dfxdsigma )
        if (associated( dfcdsigma )) deallocate( dfcdsigma )
        if (associated( ex_val )) deallocate( ex_val )
        if (associated( vxc_val ))    deallocate( vxc_val )
        if (associated( dfxdn_val ))  deallocate( dfxdn_val )
        if (associated( dfxdsigma_val )) deallocate( dfxdsigma_val )

        call glean(thy(lay))
        call bequeath(thy(vxc_g))

        call glean(thy(xcd))
        call glean(thy(n_total))
        call glean(thy(n_valence))

        if (error("Exit xc_density_libxc_mod::xcd_energy_and_potential_libxc")) continue

      end function

      subroutine xcd_grid_pressure_libxc(xcd,n,p)
!doc$ subroutine xc_grid_pressure(xcd,n,p)
        type(xc_density_libxc_obj) :: xcd
        type(grid_obj) :: n
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the exchange-correlation grid potential.

!cod$
        real(double), dimension(:,:), pointer :: s

        call my(xcd)
        call my(n)

        allocate( s(3,3) )
        call xcd_grid_stress_tensor_libxc(xcd,n,s) ; if (error()) goto 100
        p = -(s(1,1) + s(2,2) + s(3,3))/3.0_double

100     if (associated( s )) deallocate( s )

        call glean(thy(xcd))
        call glean(thy(n))

        if (error("Exit xc_density_libxc_mod::xcd_grid_pressure_libxc")) continue

      end subroutine

      subroutine xcd_grid_stress_tensor_libxc(xcd,n_g,s)
!doc$ subroutine xc_grid_stress_tensor(xcd,n_g,s)
        type(xc_density_libxc_obj) :: xcd
        type(grid_obj) :: n_g
        real(double), dimension(:,:), intent(out) :: s
!       modifies: s
!       effects: Returns stress tensor contributions due to the exchange-correlation grid potential.

!cod$
        logical :: stress
        real(double) :: sum_r
        real(double), dimension(:,:), allocatable :: s_local
        real(double), dimension(:,:,:), pointer :: n, exc, vxc, dfxcdn, dfxcddnodn, dfxcdlapl
        real(double), dimension(:,:,:), pointer :: ex, dfxdn, dfxddnodn, dfxdlapl
        real(double), dimension(:,:,:), pointer :: ec, dfcdn, dfcddnodn, dfcdlapl
        type(layout_obj) :: lay
        type(dens_der) :: nder
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info
        integer :: arraysize

        if (error(.true.,'Error: stress tensor not supported with libxc functionals (yet)')) goto 200

!        call my(xcd)
!        call my(n_g)

!        nullify( n, exc, vxc, dfxcdn, dfxcddnodn, dfxcdlapl )
!        nullify( ex, dfxdn, dfxddnodn, dfxdlapl )
!        nullify( ec, dfcdn, dfcddnodn, dfcdlapl )


!        call my(x_layout(n_g),lay)

!        call take(n,n_g,RD_KIND)

!        allocate( s_local(3,3) )

!        call alloc(exc,lay,D_TYPE,SGROUP)
!        call alloc(ex,lay,D_TYPE,SGROUP)
!        call alloc(ec,lay,D_TYPE,SGROUP)
!        call alloc(vxc,lay,D_TYPE,SGROUP)
!        call alloc(dfxcdn,lay,D_TYPE,SGROUP)
!        call alloc(dfxdn,lay,D_TYPE,SGROUP)
!        call alloc(dfcdn,lay,D_TYPE,SGROUP)
!        if (xcd%o%uses_gradient) call alloc(dfxcddnodn,lay,D_TYPE,SGROUP)
!        if (xcd%o%uses_gradient) call alloc(dfxddnodn,lay,D_TYPE,SGROUP)
!        if (xcd%o%uses_gradient) call alloc(dfcddnodn,lay,D_TYPE,SGROUP)
!        if (xcd%o%uses_laplacian) call alloc(dfxcdlapl,lay,D_TYPE,SGROUP)
!        if (xcd%o%uses_laplacian) call alloc(dfxdlapl,lay,D_TYPE,SGROUP)
!        if (xcd%o%uses_laplacian) call alloc(dfcdlapl,lay,D_TYPE,SGROUP)

!        vxc = 0.0
!        exc = 0.0
!        ex = 0.0
!        ec = 0.0
!        dfxcdn = 0.0
!        dfxdn = 0.0
!        dfcdn = 0.0
!        if (xcd%o%uses_gradient) then
!           dfxcddnodn = 0.0
!           dfxddnodn = 0.0
!           dfcddnodn = 0.0
!        end if
!        if (xcd%o%uses_laplacian) then
!           dfxcdlapl = 0.0
!           dfxdlapl = 0.0
!           dfcdlapl = 0.0
!        end if


!        stress = .true.
!        call form_nder_i(xcd,lay,n,nder,stress)

!        arraysize = size(n,1)*size(n,2)*size(n,3)

        !** Calc contribution from xc

!        call xc_f90_func_init(xc_func,xc_info,xcd%o%xc_id,xcd%o%spin_status)
        
!        select case(xc_f90_family_from_id(xcd%o%xc_id))
!        case (XC_FAMILY_LDA)
!           call xc_f90_lda_exc_vxc(xc_func, arraysize, n(1,1,1), exc(1,1,1), dfxcdn(1,1,1))
!        case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
!           call xc_f90_gga_exc_vxc(xc_func, arraysize, n(1,1,1), nder%dn(1,1,1),exc(1,1,1),dfxcdn(1,1,1),dfxcddnodn(1,1,1))
!        end select
        
!        call xc_f90_func_end(xc_func)
        
        !** Calc contribution from x
!        call xc_f90_func_init(xc_func,xc_info,xcd%o%x_id,xcd%o%spin_status)
        
!        select case(xc_f90_family_from_id(xcd%o%x_id))
!        case (XC_FAMILY_LDA)
!           call xc_f90_lda_exc_vxc(xc_func, arraysize, n(1,1,1), ex(1,1,1), dfxdn(1,1,1))
!           exc = exc + ex
!           dfxcdn = dfxcdn + dfxdn
!        case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
!           call xc_f90_gga_exc_vxc(xc_func, arraysize, n(1,1,1), nder%dn(1,1,1),ex(1,1,1),dfxdn(1,1,1),dfxddnodn(1,1,1))
!           exc = exc + ex
!           dfxcdn = dfxcdn + dfxdn
!           dfxcddnodn = dfxcddnodn + dfxddnodn
!        end select
        
!        call xc_f90_func_end(xc_func)
        
        !** Calc contribution from c
!        call xc_f90_func_init(xc_func,xc_info,xcd%o%c_id,xcd%o%spin_status)
        
!        select case(xc_f90_family_from_id(xcd%o%c_id))
!        case (XC_FAMILY_LDA)
!           call xc_f90_lda_exc_vxc(xc_func, arraysize, n(1,1,1), ec(1,1,1), dfcdn(1,1,1))
!           exc = exc + ec
!           dfxcdn = dfxcdn + dfcdn
!        case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA)
!           call xc_f90_gga_exc_vxc(xc_func, arraysize, n(1,1,1), nder%dn(1,1,1),ec(1,1,1),dfcdn(1,1,1),dfcddnodn(1,1,1))
!           exc = exc + ec
!           dfxcdn = dfxcdn + dfcdn
!           dfxcddnodn = dfxcddnodn + dfcddnodn
!        end select
 
!        call xc_f90_func_end(xc_func)

        !** Convert from Atomic Mass Units to Atomic Rydberg Units (Hartrees to Rydbergs)
!        exc = exc*Hartrees2Rydbergs
!        dfxcdn = dfxcdn*Hartrees2Rydbergs
!        if (xcd%o%uses_gradient) dfxcddnodn = dfxcddnodn*Hartrees2Rydbergs

!        select case (xcd%o%potential_method)
!        case (XCD_SIMPLE)
!          vxc = dfxcdn
!        case (XCD_WHITE_BIRD)
!          call white_bird_i(xcd,lay,nder,dfxcdn,dfxcddnodn,dfxcdlapl,vxc)
!        end select
!        sum_r = sum(n*(exc - vxc))                                           ! contributions from the density
!        s_local(1,1) = sum_r
!        s_local(1,2) = 0.0_double
!        s_local(1,3) = 0.0_double
!        s_local(2,2) = sum_r
!        s_local(2,3) = 0.0_double
!        s_local(3,3) = sum_r
!        if (xcd%o%uses_gradient) then                                         ! contributions from the gradient
!          s_local(1,1) = s_local(1,1) - sum(dfxcddnodn*nder%dnx**2)
!          s_local(1,2) = s_local(1,2) - sum(dfxcddnodn*nder%dnx*nder%dny)
!          s_local(1,3) = s_local(1,3) - sum(dfxcddnodn*nder%dnx*nder%dnz)
!          s_local(2,2) = s_local(2,2) - sum(dfxcddnodn*nder%dny**2)
!          s_local(2,3) = s_local(2,3) - sum(dfxcddnodn*nder%dny*nder%dnz)
!          s_local(3,3) = s_local(3,3) - sum(dfxcddnodn*nder%dnz**2) 
!        end if
!        if ( xcd%o%uses_laplacian ) then                                      ! contributions from the laplacian
!          s_local(1,1) = s_local(1,1) - sum(dfxcdlapl*nder%d2nxx)
!          s_local(1,2) = s_local(1,2) - sum(dfxcdlapl*nder%d2nxy)
!          s_local(1,3) = s_local(1,3) - sum(dfxcdlapl*nder%d2nxz)
!          s_local(2,2) = s_local(2,2) - sum(dfxcdlapl*nder%d2nyy)
!          s_local(2,3) = s_local(2,3) - sum(dfxcdlapl*nder%d2nyz)
!          s_local(3,3) = s_local(3,3) - sum(dfxcdlapl*nder%d2nzz)   
!        end if               
!        call kill_nder_i(nder)
!
!        s_local(2,1) = s_local(1,2)
!        s_local(3,1) = s_local(1,3)
!        s_local(3,2) = s_local(2,3)
!        call allreduce(CONFIG,MPI_SUM,s_local,s)
!        s = s/product(x_dims(lay))
         s = 0.0_double  ! Remove when commenting is gone

!100     if (allocated( s_local )) deallocate( s_local )
!        if (associated( n )) deallocate( n )
!        if (associated( exc )) deallocate( exc )
!        if (associated( vxc )) deallocate( vxc )
!        if (associated( dfxcdn )) deallocate( dfxcdn )
!        if (associated( dfxcddnodn )) deallocate( dfxcddnodn )
!        if (associated( dfxcdlapl )) deallocate( dfxcdlapl )

!        call glean(thy(lay))

!        call glean(thy(xcd))
!        call glean(thy(n_g))

200     if (error("Exit xc_density_libxc_mod::xcd_grid_stress_tensor_libxc")) continue

      end subroutine

! private routines

      subroutine white_bird_i(xcd,lay,nder,dfxcdn,dfxcdsigma,vxc)
        type(xc_density_libxc_obj) :: xcd
        type(layout_obj) :: lay
        type(dens_der) :: nder
        real(double), dimension(:,:,:,:), pointer :: dfxcdn, dfxcdsigma, vxc

        real(double), dimension(:,:,:), pointer :: r1
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: gr1

        call my(xcd)
        call my(lay)

        nullify( r1 )
        nullify( c1 )
        nullify( c2 )

        call my(grid(lay,SGROUP),gr1)

        ! contribution due to the density
        vxc = dfxcdn

        ! contribution due to the gradient
        !** Because LibXC works with contracted gradients, the GGA potential must be calculated 
        !   with a different expression than for the native Socorro implementation:
        !  
        !      vxc = dfxc/dn - 2*del.(dfxc/dsigma * del.n)
        !
        !
        !   For the spin-polarized version:
        !      
        !      vxc_up = dfxc/dn_up - del.(dfxc/dsigma_1 * del.n_up  +  dfxc/dsigma_2 * del.n_dn)
        !
        !      vxc_dn = dfxc/dn_dn - del.(dfxc/dsigma_2 * del.n_up  +  dfxc/dsigma_3 * del.n_dn)
        !
        !      where: 
        !
        !        sigma_1 = (del n_up).(del n_up)
        !        sigma_2 = (del n_up).(del n_dn)
        !        sigma_3 = (del n_dn).(del n_dn)
        !
        ! end

        if (xcd%o%uses_gradient) then
           select case (mpi_nsgroups())
           case (1) ! Spin-independent case.   
              ! form x component of (dfxc/dsigma * dn/dx)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = dfxcdsigma(1,:,:,:)*nder%dnx(1,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the x component of the divergence d/dx (dfxc/dsigma * dn/dx)
              call alloc(c2,lay,D_TYPE,SGROUP)
              c2 = c1*xcd%o%gx
              deallocate( c1 )


              ! form y component of (dfxc/dsigma * dn/dy)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = dfxcdsigma(1,:,:,:)*nder%dny(1,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the y component of the divergence: d/dy(dfxc/dsigma * dn/dy)
              c2 = c2 + c1*xcd%o%gy
              deallocate( c1 )

              ! form z component of (dfxc/dsigma * dn/dz)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = dfxcdsigma(1,:,:,:)*nder%dnz(1,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the z component of the divergence d/dz(dfxc/dsigma * dn/dz)
              c2 = c2 + c1*xcd%o%gz
              deallocate( c1 )

              ! Move del.(dfxc/dsigma * del.n) back to real space. 
              c2 = (0.0_double,1.0_double)*c2

              call put(c2,gr1,CDF_KIND)
              call take(r1,gr1,RD_KIND)

              !** Note factor of two to account for the contracted gradient.
              vxc(1,:,:,:) = vxc(1,:,:,:) - 2.0_double*r1

           case (2) ! Spin dependent.
              ! First calculate the contribution to vxc_up:

              ! Form x component of dfxc/dsigma_1 * dn_up/dx  +  dfxc/dsigma_2 * dn_dn/dx)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = 2.0_double*dfxcdsigma(1,:,:,:)*nder%dnx(1,:,:,:) + dfxcdsigma(2,:,:,:)*nder%dnx(2,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the x component of the divergence d/dx (2*dfxc/dsigma_1 * dn_up/dx  +  dfxc/dsigma_2 * dn_dn/dx))
              call alloc(c2,lay,D_TYPE,SGROUP)
              c2 = c1*xcd%o%gx
              deallocate( c1 )

              ! Form y component of dfxc/dsigma_1 * dn_up/dy  +  dfxc/dsigma_2 * dn_dn/dy)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = 2.0_double*dfxcdsigma(1,:,:,:)*nder%dny(1,:,:,:) + dfxcdsigma(2,:,:,:)*nder%dny(2,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the y component of the divergence d/dy (2*dfxc/dsigma_1 * dn_up/dy  +  dfxc/dsigma_2 * dn_dn/dy))
              c2 = c2 + c1*xcd%o%gy
              deallocate( c1 )

              ! Form z component of dfxc/dsigma_1 * dn_up/dz  +  dfxc/dsigma_2 * dn_dn/dz)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = 2.0_double*dfxcdsigma(1,:,:,:)*nder%dnz(1,:,:,:) + dfxcdsigma(2,:,:,:)*nder%dnz(2,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the z component of the divergence d/dz (2*dfxc/dsigma_1 * dn_up/dz  +  dfxc/dsigma_2 * dn_dn/dz))
              c2 = c2 + c1*xcd%o%gz
              deallocate( c1 )

              ! Move del.(dfxc/dsigma * del.n) back to real space. 
              c2 = (0.0_double,1.0_double)*c2
              call put(c2,gr1,CDF_KIND)
              call take(r1,gr1,RD_KIND)

              ! Spin up component of vxc
              vxc(1,:,:,:) = vxc(1,:,:,:) - r1

              !** Next calculate the contribution to the down spin component of vxc:
              ! Form x component of dfxc/dsigma_2 * dn_up/dx  +  dfxc/dsigma_3 * dn_dn/dx)
              r1 = dfxcdsigma(2,:,:,:)*nder%dnx(1,:,:,:) + 2.0_double*dfxcdsigma(3,:,:,:)*nder%dnx(2,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the x component of the divergence d/dx (dfxc/dsigma_2 * dn_up/dx  +  2*dfxc/dsigma_3 * dn_dn/dx))
              call alloc(c2,lay,D_TYPE,SGROUP)
              c2 = c1*xcd%o%gx
              deallocate( c1 )

              ! Form y component of dfxc/dsigma_2 * dn_up/dy  +  dfxc/dsigma_3 * dn_dn/dy)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = dfxcdsigma(2,:,:,:)*nder%dny(1,:,:,:) + 2.0_double*dfxcdsigma(3,:,:,:)*nder%dny(2,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the y component of the divergence d/dy (dfxc/dsigma_2 * dn_up/dy  +  2*dfxc/dsigma_3 * dn_dn/dy))
              c2 = c2 + c1*xcd%o%gy
              deallocate( c1 )

              ! Form z component of dfxc/dsigma_2 * dn_up/dz  +  dfxc/dsigma_3 * dn_dn/dz)
              call alloc(r1,lay,D_TYPE,SGROUP)
              r1 = dfxcdsigma(2,:,:,:)*nder%dnz(1,:,:,:) + 2.0_double*dfxcdsigma(3,:,:,:)*nder%dnz(2,:,:,:)
              call put(r1,gr1,RD_KIND)
              call take(c1,gr1,CDF_KIND)
              ! calculate the z component of the divergence d/dz (dfxc/dsigma_2 * dn_up/dz  +  2*dfxc/dsigma_3 * dn_dn/dz))
              c2 = c2 + c1*xcd%o%gz
              deallocate( c1 )

              ! Move del.(dfxc/dsigma * del.n) back to real space. 
              c2 = (0.0_double,1.0_double)*c2
              call put(c2,gr1,CDF_KIND)
              call take(r1,gr1,RD_KIND)

              ! Spin down component of vxc
              vxc(2,:,:,:) = vxc(2,:,:,:) - r1

           end select

        end if

        ! contribution due to the laplacian
        if (xcd%o%uses_laplacian) then
           if (error(.true.,'laplacian not supported for LIBXC')) goto 100
          !r1 = dfxcdlapl
          !call put(r1,gr1,RD_KIND)
          !call take(c1,gr1,CDF_KIND)
          !c1 = c1*(xcd%o%gx**2 + xcd%o%gy**2 + xcd%o%gz**2)
          !call put(c1,gr1,CDF_KIND)
          !call take(r1,gr1,RD_KIND)
          !vxc = vxc - r1
        end if

100     if (associated( r1 )) deallocate( r1 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(gr1))

        call glean(thy(xcd))
        call glean(thy(lay))

        if (error("Exit xc_density_libxc_mod::white_bird_i")) continue

      end subroutine

      subroutine form_nder_i(xcd,lay,n,nder,stress)
        type(xc_density_libxc_obj) :: xcd
        type(layout_obj) :: lay
        real(double), dimension(:,:,:) :: n
        type(dens_der) :: nder
        logical, intent(in), optional :: stress

        logical :: stress_on
        integer :: nx, ny, nz
        real(double), dimension(:,:,:), pointer :: r1
        complex(double), dimension(:,:,:), pointer :: c1, c2
        real(double), dimension(:,:,:,:), pointer :: n_sg, dnx_sg, dny_sg, dnz_sg
        real(double), dimension(:,:,:), pointer :: dnx, dny, dnz
        type(grid_obj) :: gr1

        call my(xcd)
        call my(lay)

        nullify( nder%n )
        nullify( nder%dnx )
        nullify( nder%dny )
        nullify( nder%dnz )
        nullify( nder%sigma )

        nullify( r1 )
        nullify( c1 )
        nullify( c2 )

        nullify( n_sg )
        nullify( dnx_sg )
        nullify( dny_sg )
        nullify( dnz_sg )

        nullify( dnx ) 
        nullify( dny ) 
        nullify( dnz ) 

        nx = size(n,1)
        ny = size(n,2)
        nz = size(n,3)

        stress_on = .false.
        if (present(stress)) stress_on = stress

        call my(grid(lay,SGROUP),gr1)

        ! initialize the density array
        select case (mpi_nsgroups())
        case (1)

          allocate( nder%n(1,nx,ny,nz) )
          nder%n(1,:,:,:) = n

        case (2)

          allocate( nder%n(2,nx,ny,nz) )

          allocate( n_sg(2,nx,ny,nz) )
          n_sg = 0.0_double
          if (mpi_mysgroup() == 1) then
            n_sg(1,:,:,:) = n(:,:,:)
          elseif (mpi_mysgroup() == 2) then
            n_sg(2,:,:,:) = n(:,:,:)
          end if
          call xcomm_rank_allreduce(XSGROUP,MPI_SUM,n_sg,nder%n)

          deallocate( n_sg )

        end select

        ! set a minimum density for numerical stability
        where (nder%n < mindens) nder%n = mindens

        call alloc(r1,lay,D_TYPE,SGROUP)
        r1 = n
        call put(r1,gr1,RD_KIND)
        call take(c1,gr1,CDF_KIND)

        ! initialize the density-gradient arrays
        if (xcd%o%uses_gradient) then

          ! calculate the gradient and store in dnx, dny, dnz
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = (0.0_double,1.0_double)*c1*xcd%o%gx
          call put(c2,gr1,CDF_KIND)
          call take(dnx,gr1,RD_KIND)
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = (0.0_double,1.0_double)*c1*xcd%o%gy
          call put(c2,gr1,CDF_KIND)
          call take(dny,gr1,RD_KIND)
          call alloc(c2,lay,D_TYPE,SGROUP)
          c2 = (0.0_double,1.0_double)*c1*xcd%o%gz
          call put(c2,gr1,CDF_KIND)
          call take(dnz,gr1,RD_KIND)

          ! load the gradient arrays
          select case (mpi_nsgroups())
          case (1)

            allocate( nder%dnx(1,nx,ny,nz) )
            allocate( nder%dny(1,nx,ny,nz) )
            allocate( nder%dnz(1,nx,ny,nz) )
            allocate( nder%sigma(1,nx,ny,nz) )

            nder%dnx(1,:,:,:) = dnx(:,:,:)
            nder%dny(1,:,:,:) = dny(:,:,:)
            nder%dnz(1,:,:,:) = dnz(:,:,:)

            ! set a minimum gradient for numerical stability
            where (nder%n <= mindens)
              nder%dnx = mindens
              nder%dny = mindens
              nder%dnz = mindens
            end where

            nder%sigma(1,:,:,:) = nder%dnx(1,:,:,:)*nder%dnx(1,:,:,:) + &
                                  nder%dny(1,:,:,:)*nder%dny(1,:,:,:) + &
                                  nder%dnz(1,:,:,:)*nder%dnz(1,:,:,:)

          case (2)

            allocate( nder%dnx(2,nx,ny,nz) )
            allocate( nder%dny(2,nx,ny,nz) )
            allocate( nder%dnz(2,nx,ny,nz) )
            allocate( nder%sigma(3,nx,ny,nz) )

            allocate( dnx_sg(2,nx,ny,nz) )
            allocate( dny_sg(2,nx,ny,nz) )
            allocate( dnz_sg(2,nx,ny,nz) )
            dnx_sg = 0.0_double
            dny_sg = 0.0_double
            dnz_sg = 0.0_double
            if (mpi_mysgroup() == 1) then
              dnx_sg(1,:,:,:) = dnx(:,:,:) 
              dny_sg(1,:,:,:) = dny(:,:,:)  
              dnz_sg(1,:,:,:) = dnz(:,:,:)  
            elseif (mpi_mysgroup() == 2) then
              dnx_sg(2,:,:,:) = dnx(:,:,:) 
              dny_sg(2,:,:,:) = dny(:,:,:) 
              dnz_sg(2,:,:,:) = dnz(:,:,:) 
            end if
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,dnx_sg,nder%dnx)
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,dny_sg,nder%dny)
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,dnz_sg,nder%dnz)

            ! set a minimum gradient for numerical stability
            where (nder%n <= mindens)
              nder%dnx = mindens
              nder%dny = mindens
              nder%dnz = mindens
            end where
             
            ! create the sigma array needed in LibXC
            ! sigma(1) = (del n up).(del n up)
            nder%sigma(1,:,:,:) = nder%dnx(1,:,:,:)*nder%dnx(1,:,:,:) + &
                                  nder%dny(1,:,:,:)*nder%dny(1,:,:,:) + &
                                  nder%dnz(1,:,:,:)*nder%dnz(1,:,:,:)

            ! sigma(2) = (del n up).(del n down)
            nder%sigma(2,:,:,:) = nder%dnx(1,:,:,:)*nder%dnx(2,:,:,:) + &
                                  nder%dny(1,:,:,:)*nder%dny(2,:,:,:) + &
                                  nder%dnz(1,:,:,:)*nder%dnz(2,:,:,:)

            ! sigma(3) = (del n down).(del n down)
            nder%sigma(3,:,:,:) = nder%dnx(2,:,:,:)*nder%dnx(2,:,:,:) + &
                                  nder%dny(2,:,:,:)*nder%dny(2,:,:,:) + &
                                  nder%dnz(2,:,:,:)*nder%dnz(2,:,:,:)

          end select
          
        end if
        
        if (xcd%o%uses_laplacian) then
          if (error(.true.,"Error: laplacian not supported with libxc")) goto 100
!          call alloc(c2,lay,D_TYPE,SGROUP)
!          c2 = -c1*(xcd%o%gx**2 + xcd%o%gy**2 + xcd%o%gz**2)
!          call put(c2,gr1,CDF_KIND)
!          call take(nder%lapl,gr1,RD_KIND)
!          where (nder%n <= mindens) nder%lapl = mindens
!          if (stress_on) then
!            call alloc(c2,lay,D_TYPE,SGROUP)
!            c2 = -c1*xcd%o%gx**2
!            call put(c2,gr1,CDF_KIND)
!            call take(nder%d2nxx,gr1,RD_KIND)
!            call alloc(c2,lay,D_TYPE,SGROUP)
!            c2 = -c1*xcd%o%gx*xcd%o%gy
!            call put(c2,gr1,CDF_KIND)
!            call take(nder%d2nxy,gr1,RD_KIND)
!            call alloc(c2,lay,D_TYPE,SGROUP)
!            c2 = -c1*xcd%o%gx*xcd%o%gz
!            call put(c2,gr1,CDF_KIND)
!            call take(nder%d2nxz,gr1,RD_KIND)
!            call alloc(c2,lay,D_TYPE,SGROUP)
!            c2 = -c1*xcd%o%gy**2
!            call put(c2,gr1,CDF_KIND)
!            call take(nder%d2nyy,gr1,RD_KIND)
!            call alloc(c2,lay,D_TYPE,SGROUP)
!            c2 = -c1*xcd%o%gy*xcd%o%gz
!            call put(c2,gr1,CDF_KIND)
!            call take(nder%d2nyz,gr1,RD_KIND)
!            call alloc(c2,lay,D_TYPE,SGROUP)
!            c2 = -c1*xcd%o%gz**2
!            call put(c2,gr1,CDF_KIND)
!            call take(nder%d2nzz,gr1,RD_KIND)
!            where (nder%n <= mindens)
!              nder%d2nxx = mindens
!              nder%d2nxy = mindens
!              nder%d2nxz = mindens
!              nder%d2nyy = mindens
!              nder%d2nyz = mindens
!              nder%d2nzz = mindens
!            end where
!          end if
        end if

        if (associated( r1 ))     deallocate( r1 )
        if (associated( c1 ))     deallocate( c1 )
        if (associated( c2 ))     deallocate( c2 )
        if (associated( n_sg ))   deallocate( n_sg )
        if (associated( dnx_sg )) deallocate( dnx_sg )
        if (associated( dny_sg )) deallocate( dny_sg )
        if (associated( dnz_sg )) deallocate( dnz_sg )        
        if (associated( dnx ))    deallocate( dnx )        
        if (associated( dny ))    deallocate( dny )        
        if (associated( dnz ))    deallocate( dnz )        

        call glean(thy(gr1))

        call glean(thy(xcd))
        call glean(thy(lay))

100     if (error("xcd_libxc::form_nder_i")) continue

      end subroutine

      subroutine kill_nder_i(nder)
        type(dens_der) :: nder

        if (associated( nder%n))     deallocate( nder%n )
        if (associated( nder%dnx))   deallocate( nder%dnx )
        if (associated( nder%dny))   deallocate( nder%dny )
        if (associated( nder%dnz))   deallocate( nder%dnz )
        if (associated( nder%sigma)) deallocate( nder%sigma )

      end subroutine

      subroutine own_i(xcd)
        type(xc_density_libxc_obj) :: xcd
        type(xc_density_libxc_obj) :: xcd_t
        if (xcd%ref < xcd%o%ref) then
          allocate( xcd_t%o )
          xcd_t%o%ref               = 0
          xcd_t%o%g                 = xcd%o%g
          xcd_t%o%xc_id             = xcd%o%xc_id
          xcd_t%o%x_id              = xcd%o%x_id
          xcd_t%o%c_id              = xcd%o%c_id
          xcd_t%o%potential_method  = xcd%o%potential_method
          xcd_t%o%uses_gradient     = xcd%o%uses_gradient
          xcd_t%o%uses_laplacian    = xcd%o%uses_laplacian
          xcd_t%o%spin_status       = xcd%o%spin_status
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

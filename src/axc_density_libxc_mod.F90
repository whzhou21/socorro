! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module axc_density_libxc_mod
!doc$ module axc_density_libxc_mod

!     One datatype is available here: type(axc_density_libxc_obj).

!     axc_density_libxc_mod encapsulates libxc routines used to compute density-dependent
!     exchange-correlation quantities for atomic data.

      use kind_mod
      use mpi_mod
      use error_mod
      use ghost_mod
      use xc_type_mod
      use xc_f90_types_m
      use libxc_funcs_m
      use xc_f90_lib_m

!cod$
      implicit none
      private 

      real(double), parameter :: mindens = 1.0e-20_double
      real(double), parameter :: Hartrees2Rydbergs = 2.0_double

      type :: axc_density_libxc_rep
        integer :: ref                                        ! reference count
        type(ghost) :: g                                      ! ghost
        integer :: xc_id                                      ! exchange-correlation type
        integer :: x_id                                       ! exchange type
        integer :: c_id                                       ! correlation type
        logical :: uses_gradient                              ! indicates whether or not the gradient is used
        logical :: uses_laplacian                             ! indicates whether or not the laplacian is used
        integer :: spin_status                                ! spin status indicator
        type(xc_type_obj) :: xct                              ! xc_type object
      end type

      type, public :: axc_density_libxc_obj
        private
        integer :: ref
        type(axc_density_libxc_rep), pointer :: o
      end type

      type :: dens_der
        real(double), dimension(:,:), pointer :: n
        real(double), dimension(:,:), pointer :: dn
        real(double), dimension(:,:), pointer :: sigma
      end type

      public :: axc_density_libxc
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

      interface axc_density_libxc
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
        module procedure axcd_energy_density_libxc
      end interface
      interface xc_derivatives
        module procedure axcd_derivatives_libxc
      end interface

      contains

! public routines

      function constructor_axcd(xct) result(axcd)
!doc$ function axc_density_libxc(xct) result(axcd)
        type(xc_type_obj) :: xct
        type(axc_density_libxc_obj) :: axcd
!       effects: Constructs a new axcd.
!       errors: Passes errors

!cod$ 
        call my(xct)

        axcd%ref = 0
        allocate( axcd%o )
        axcd%o%ref = 0
        axcd%o%g = x_ghost()

        call my(xct,axcd%o%xct)

        axcd%o%xc_id = x_functional_type(axcd%o%xct)
        axcd%o%x_id = x_exchange_type(axcd%o%xct)
        axcd%o%c_id = x_correlation_type(axcd%o%xct)
        axcd%o%uses_gradient = uses_gradient(axcd%o%xct)
        axcd%o%uses_laplacian = uses_laplacian(axcd%o%xct)

        select case (mpi_nsgroups())
        case (1)
          axcd%o%spin_status = XC_UNPOLARIZED
        case (2)
          axcd%o%spin_status = XC_POLARIZED
        end select

100     call glean(thy(xct))

        if (error("Exit axc_density_libxc_mod::constructor_axcd")) continue

      end function 

      subroutine update_axcd(axcd)
!doc$ subroutine update(axcd)
        type(axc_density_libxc_obj) :: axcd
!       modifies: axcd
!       effects: Updates axcd.
!       note: axcd%o%xct is not allowed to change.

!cod$
        call my(axcd)

        continue

100     call glean(thy(axcd))

      end subroutine

      subroutine my_axcd(axcd)
!doc$ subroutine my(axcd)
        type(axc_density_libxc_obj) :: axcd

!cod$
        axcd%ref = axcd%ref + 1
        axcd%o%ref = axcd%o%ref + 1
      end subroutine

      subroutine my_new_axcd(axcdi,axcd)
!doc$ subroutine my(axcdi,axcd)
        type(axc_density_libxc_obj) :: axcdi, axcd

!cod$
        axcd%ref = 1
        axcd%o => axcdi%o
        axcd%o%ref = axcd%o%ref + 1
      end subroutine
      
      function thy_axcd(axcd) result(axcdo)
!doc$ function thy(axcd) result(axcdo)
        type(axc_density_libxc_obj) :: axcd, axcdo

!cod$
        axcd%ref = axcd%ref - 1
        axcd%o%ref = axcd%o%ref - 1
        axcdo%ref = axcd%ref
        axcdo%o => axcd%o
      end function

      subroutine glean_axcd(axcd)
!doc$ subroutine glean(axcd)
        type(axc_density_libxc_obj) :: axcd

!cod$
        if (axcd%o%ref < 1) then
          call glean(thy(axcd%o%xct))
          deallocate( axcd%o )
        end if
      end subroutine

      subroutine bequeath_axcd(axcd)
!doc$ subroutine bequeath(axcd)
        type(axc_density_libxc_obj) :: axcd

!cod$
        continue
      end subroutine

      subroutine assign_axcd(axcd,axcd2)
!doc$ subroutine assign(axcd,axcd2)
        type(axc_density_libxc_obj), intent(inout) :: axcd
        type(axc_density_libxc_obj), intent(in) :: axcd2

!cod$
        type(axc_density_libxc_obj) :: axcdt
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
        type(axc_density_libxc_obj) :: axcd
        integer, dimension(2) :: r
!       effects: Returns axcd%ref and axcd%o%ref.

!cod$
        r(1) = axcd%ref
        r(2) = axcd%o%ref
        call glean(axcd)
      end function

      function axcd_ghost(axcd) result(g)
!doc$ function x_ghost(axcd) result(g)
        type(axc_density_libxc_obj) :: axcd
        type(ghost) :: g
!       effects: Returns the ghost of axcd.

!cod$
        call my(axcd)
        g = axcd%o%g
        call glean(thy(axcd))
      end function

      function axcd_xc_type(axcd) result(xct)
!doc$ function x_xc_type(axcd) result(xct)
        type(axc_density_libxc_obj) :: axcd
        type(xc_type_obj)  :: xct 
!       effects: Returns the xc_type of axcd.

!cod$
        call my(axcd)
        call my(axcd%o%xct,xct)
        call glean(thy(axcd))
        call bequeath(thy(xct))
      end function

      function q_uses_gradient(axcd) result(ug)
!doc$ function uses_gradient(axcd) result(ug)
        type(axc_density_libxc_obj) :: axcd
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        call my(axcd)
        ug = axcd%o%uses_gradient 
        call glean(thy(axcd))
      end function

      function q_uses_laplacian(axcd) result(ul)
!doc$ function uses_laplacian(axcd) result(ul)
        type(axc_density_libxc_obj) :: axcd
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        call my(axcd)
        ul = axcd%o%uses_laplacian
        call glean(thy(axcd))
      end function

      subroutine axcd_energy_density_libxc(axcd,n,dn,lapl,exc)
!doc$ subroutine xc_energy_density(axcd,n,dn,lapl,exc)
        type(axc_density_libxc_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lapl, exc
!       requires: Pointers be nullified or associated.
!       effects: Returns exc with respect to n and dn (lapl not used at this time).
!       errors: Passes errors.

!cod$
        integer :: np
        real(double), dimension(:), pointer :: ex, ec
        type(dens_der) :: nder
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        call my(axcd)

        nullify( ex )
        nullify( ec )

        call form_nder_i(n,dn,nder)

        np = size(n)

        exc = 0.0_double

        ! contribution from xc
        select case(xc_f90_family_from_id(axcd%o%xc_id))
        case (XC_FAMILY_LDA)
          call xc_f90_func_init(xc_func,xc_info,axcd%o%xc_id,axcd%o%spin_status)
          call xc_f90_lda_exc(xc_func,np,nder%n(1,1),exc(1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_GGA)
          call xc_f90_func_init(xc_func,xc_info,axcd%o%xc_id,axcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder%n(1,1),nder%sigma(1,1),exc(1))
          call xc_f90_func_end(xc_func)
        end select

        ! contribution from x
        select case(xc_f90_family_from_id(axcd%o%x_id))
        case (XC_FAMILY_LDA)
          allocate( ex(np) )
          call xc_f90_func_init(xc_func,xc_info,axcd%o%x_id,axcd%o%spin_status)
          call xc_f90_lda_exc(xc_func,np,nder%n(1,1),ex(1))
          call xc_f90_func_end(xc_func)
          exc = exc + ex
          deallocate( ex )
        case (XC_FAMILY_GGA)
          allocate( ex(np) )
          call xc_f90_func_init(xc_func,xc_info,axcd%o%x_id,axcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder%n(1,1),nder%sigma(1,1),ex(1))
          call xc_f90_func_end(xc_func)
          exc = exc + ex
          deallocate( ex )
        end select

        ! contribution from c
        select case(xc_f90_family_from_id(axcd%o%c_id))
        case (XC_FAMILY_LDA)
          allocate( ec(np) )
          call xc_f90_func_init(xc_func,xc_info,axcd%o%c_id,axcd%o%spin_status)
          call xc_f90_lda_exc(xc_func,np,nder%n(1,1),ec(1))
          call xc_f90_func_end(xc_func)
          exc = exc + ec
          deallocate( ec )
        case (XC_FAMILY_GGA)
          allocate( ec(np) )
          call xc_f90_func_init(xc_func,xc_info,axcd%o%c_id,axcd%o%spin_status)
          call xc_f90_gga_exc(xc_func,np,nder%n(1,1),nder%sigma(1,1),ec(1))
          call xc_f90_func_end(xc_func)
          exc = exc + ec
          deallocate( ec )
        end select

        exc = exc*Hartrees2Rydbergs

        call kill_nder_i(nder)

100     if (associated( ex )) deallocate( ex )
        if (associated( ec )) deallocate( ec )

        call glean(thy(axcd))

        if (error("Exit axc_density_libxc_mod::axcd_energy_density_libxc")) continue

      end subroutine

      subroutine axcd_derivatives_libxc(axcd,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xc_derivatives(axcd,n,dn,lapl,dfxcdn,dfxcddnodn,dfxcdlapl)
        type(axc_density_libxc_obj) :: axcd
        real(double), dimension(:), pointer :: n, dn, lapl, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Pointers be nullified or associated.
!       effects: Returns derivatives with respect to n, dn, and lapl.
!       errors: Passes errors.

!cod$
        integer :: np
        real(double), dimension(:,:), pointer :: dfxcdn_c, dfxdn_c, dfcdn_c
        real(double), dimension(:,:), pointer :: dfxcdsigma, dfxdsigma, dfcdsigma
        type(dens_der) :: nder
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        call my(axcd)

        nullify( dfxcdn_c )
        nullify( dfxdn_c )
        nullify( dfcdn_c )
        nullify( dfxcdsigma )
        nullify( dfxdsigma )
        nullify( dfcdsigma )

        call form_nder_i(n,dn,nder)

        np = size(n)

        select case (axcd%o%spin_status)
        case (XC_UNPOLARIZED)
          allocate( dfxcdn_c(1,np) )
          allocate( dfxcdsigma(1,np) )
        case (XC_POLARIZED)
          allocate( dfxcdn_c(2,np) )
          allocate( dfxcdsigma(3,np) )
        end select
        dfxcdn_c = 0.0_double
        dfxcdsigma = 0.0_double

        if (associated(dfxcddnodn)) dfxcddnodn = 0.0_double

        ! contribution from xc
        select case(xc_f90_family_from_id(axcd%o%xc_id))
        case (XC_FAMILY_LDA)
          call xc_f90_func_init(xc_func,xc_info,axcd%o%xc_id,axcd%o%spin_status)
          call xc_f90_lda_vxc(xc_func,np,nder%n(1,1),dfxcdn_c(1,1))
          call xc_f90_func_end(xc_func)
        case (XC_FAMILY_GGA)
          call xc_f90_func_init(xc_func,xc_info,axcd%o%xc_id,axcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder%n(1,1),nder%sigma(1,1),dfxcdn_c(1,1),dfxcdsigma(1,1))
          call xc_f90_func_end(xc_func)
        end select

        ! contribution from x
        select case(xc_f90_family_from_id(axcd%o%x_id))
        case (XC_FAMILY_LDA)
          select case (axcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn_c(1,np) )
          case (XC_POLARIZED)
            allocate( dfxdn_c(2,np) )
          end select
          call xc_f90_func_init(xc_func,xc_info,axcd%o%x_id,axcd%o%spin_status)
          call xc_f90_lda_vxc(xc_func,np,nder%n(1,1),dfxdn_c(1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn_c = dfxcdn_c + dfxdn_c
          deallocate( dfxdn_c )
        case (XC_FAMILY_GGA)
          select case (axcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfxdn_c(1,np) )
            allocate( dfxdsigma(1,np) )
          case (XC_POLARIZED)
            allocate( dfxdn_c(2,np) )
            allocate( dfxdsigma(3,np) )
          end select
          call xc_f90_func_init(xc_func,xc_info,axcd%o%x_id,axcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder%n(1,1),nder%sigma(1,1),dfxdn_c(1,1),dfxdsigma(1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn_c = dfxcdn_c + dfxdn_c
          dfxcdsigma = dfxcdsigma + dfxdsigma
          deallocate( dfxdn_c )
          deallocate( dfxdsigma )
        end select

        ! contribution from c
        select case(xc_f90_family_from_id(axcd%o%c_id))
        case (XC_FAMILY_LDA)
          select case (axcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfcdn_c(1,np) )
          case (XC_POLARIZED)
            allocate( dfcdn_c(2,np) )
          end select
          call xc_f90_func_init(xc_func,xc_info,axcd%o%c_id,axcd%o%spin_status)
          call xc_f90_lda_vxc(xc_func,np,nder%n(1,1),dfcdn_c(1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn_c = dfxcdn_c + dfcdn_c
          deallocate( dfcdn_c )
        case (XC_FAMILY_GGA)
          select case (axcd%o%spin_status)
          case (XC_UNPOLARIZED)
            allocate( dfcdn_c(1,np) )
            allocate( dfcdsigma(1,np) )
          case (XC_POLARIZED)
            allocate( dfcdn_c(2,np) )
            allocate( dfcdsigma(3,np) )
          end select
          call xc_f90_func_init(xc_func,xc_info,axcd%o%c_id,axcd%o%spin_status)
          call xc_f90_gga_vxc(xc_func,np,nder%n(1,1),nder%sigma(1,1),dfcdn_c(1,1),dfcdsigma(1,1))
          call xc_f90_func_end(xc_func)
          dfxcdn_c = dfxcdn_c + dfcdn_c
          dfxcdsigma = dfxcdsigma + dfcdsigma
          deallocate( dfcdn_c )
          deallocate( dfcdsigma )
        end select

        dfxcdn_c = dfxcdn_c*Hartrees2Rydbergs
        dfxcdsigma = dfxcdsigma*Hartrees2Rydbergs

        select case (axcd%o%spin_status)
        case (XC_UNPOLARIZED)
          dfxcdn(:) = dfxcdn_c(1,:)
        case (XC_POLARIZED)
          if (mpi_mysgroup() == 1) then
            dfxcdn(:) = dfxcdn_c(1,:)
          else
            dfxcdn(:) = dfxcdn_c(2,:)
          end if
        end select

        !** Socorro expects the derivative of exc wrt the |grad n| but libxc returns the
        !    derivative wrt sigma (the contracted gradients):
        !
        !        For the spin independent case:
        !        sigma = (del n).(del n)  
        !     
        !        For the spin dependent case:
        !        sigma_1 = (del n_up).(del n_up)
        !        sigma_2 = (del n_up).(del n_dn)
        !        sigma_3 = (del n_dn).(del n_dn)
        !
        !    So must transform  
        !
        !        dExc/dsigma -->  dexc/d|del n|
        !
        !        dExc/d|del n| = 2*|del n| dExc/dsigma   For the spin independent case
        !
        !        dExc/d|del n_up| = 2*|del n_up| dExc/dsigma_1    For spin dependent, spin up
        !        dExc/d|del n_dn| = 2*|del n_dn| dExc/dsigma_3    For spin dependent, spin down
        !

        if (associated(dfxcddnodn)) then
          select case (axcd%o%spin_status)
          case (XC_UNPOLARIZED)
            dfxcddnodn(:) = 2.0_double*nder%dn(1,:)*dfxcdsigma(1,:)
          case (XC_POLARIZED)
            if (mpi_mysgroup() == 1) then
              dfxcddnodn(:) = 2.0_double*nder%dn(1,:)*dfxcdsigma(1,:)
            else
              dfxcddnodn(:) = 2.0_double*nder%dn(2,:)*dfxcdsigma(3,:)
            end if
          end select
        end if

        call kill_nder_i(nder)

100     if (associated( dfxcdn_c ))   deallocate( dfxcdn_c )
        if (associated( dfxdn_c ))    deallocate( dfxdn_c )
        if (associated( dfcdn_c ))    deallocate( dfcdn_c )
        if (associated( dfxcdsigma )) deallocate( dfxcdsigma )
        if (associated( dfxdsigma ))  deallocate( dfxdsigma )
        if (associated( dfcdsigma ))  deallocate( dfcdsigma )

        call glean(thy(axcd))

        if (error("Exit axc_density_libxc_mod::axcd_derivatives")) continue

      end subroutine

! private routines

      subroutine form_nder_i(n,dn,nder)
        real(double), dimension(:), pointer :: n, dn  !, lapl
        type(dens_der) :: nder

        real(double), dimension(:,:), pointer :: n_sg, dn_sg  !, lapl_sg
        integer :: numr
        
        nullify( nder%n ) 
        nullify( nder%dn ) 
        nullify( nder%sigma ) 

        nullify( n_sg )
        nullify( dn_sg )

        numr = size(n)

        select case (mpi_nsgroups())
        case (1)
          allocate( nder%n(1,numr) )
          nder%n(1,:) = n
        case (2)
          allocate( n_sg(2,numr) )
          n_sg = 0.0_double
          if (mpi_mysgroup() == 1) then
            n_sg(1,:) = n(:)
          elseif (mpi_mysgroup() == 2) then
            n_sg(2,:) = n(:)
          end if
          allocate( nder%n(2,numr) )
          call xcomm_pair_allreduce(XSGROUP,MPI_SUM,n_sg,nder%n)
        end select
        where (nder%n < mindens) nder%n = mindens

        if (associated( dn )) then

          select case (mpi_nsgroups())
          case (1)
            allocate( nder%dn(1,numr) )
            nder%dn(1,:) = dn
            where (nder%n <= mindens) nder%dn = mindens
            allocate( nder%sigma(1,numr) )
            nder%sigma(1,:) = dn(:)*dn(:)
          case (2)
            allocate( dn_sg(2,numr) )
            dn_sg = 0.0_double
            nder%dn = 0.0_double
            nder%sigma = 0.0_double
            if (mpi_mysgroup() == 1) then
              dn_sg(1,:) = dn(:)
            elseif (mpi_mysgroup() == 2) then
              dn_sg(2,:) = dn(:)
            end if
            allocate( nder%dn(2,numr) )
            call xcomm_pair_allreduce(XSGROUP,MPI_SUM,dn_sg,nder%dn)
            where (nder%n(1,:) <= mindens) nder%dn(1,:) = mindens
            where (nder%n(2,:) <= mindens) nder%dn(2,:) = mindens
            allocate( nder%sigma(3,numr) )
            nder%sigma(1,:) = nder%dn(1,:)*nder%dn(1,:)
            nder%sigma(2,:) = nder%dn(1,:)*nder%dn(2,:)
            nder%sigma(3,:) = nder%dn(2,:)*nder%dn(2,:)
          end select

        end if

        if (associated( n_sg ))  deallocate( n_sg )
        if (associated( dn_sg )) deallocate( dn_sg )

      end subroutine

      subroutine kill_nder_i(nder)
        type(dens_der) :: nder

        if (associated( nder%n ))     deallocate( nder%n )
        if (associated( nder%dn ))    deallocate( nder%dn )
        if (associated( nder%sigma )) deallocate( nder%sigma )

      end subroutine

      subroutine own_i(axcd)
        type(axc_density_libxc_obj) :: axcd
        type(axc_density_libxc_obj) :: axcd_t
        if (axcd%ref < axcd%o%ref) then
          allocate( axcd_t%o )
          axcd_t%o%ref               = 0
          axcd_t%o%g                 = axcd%o%g
          axcd_t%o%xc_id             = axcd%o%xc_id
          axcd_t%o%x_id              = axcd%o%x_id
          axcd_t%o%c_id              = axcd%o%c_id
          axcd_t%o%uses_gradient     = axcd%o%uses_gradient
          axcd_t%o%uses_laplacian    = axcd%o%uses_laplacian
          axcd_t%o%spin_status       = axcd%o%spin_status
          call my(axcd%o%xct,axcd_t%o%xct)
          axcd%o%ref = axcd%o%ref - axcd%ref
          axcd%o => axcd_t%o
          axcd%o%ref = axcd%o%ref + axcd%ref
        end if
      end subroutine

      end module

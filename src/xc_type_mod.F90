!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module xc_type_mod
!doc$ module xc_type_mod

!     One datatype is available here: type(xc_type_obj).

!     xc_type_mod encapsulates parameters that define the formulation and methods
!     of exchange and correlation

      use arg_mod
      use diary_mod
      use error_mod
      use ghost_mod
      use io_mod
      use kind_mod
      use libxc_funcs_m
      use mpi_mod
      use tagio_mod
      use xc_f90_lib_m
      use xc_f90_types_m

!cod$
      implicit none
      private 

      logical, parameter      :: L_NULL = .false.     ! Logical null value
      integer, parameter      :: I_NULL = 0           ! Integer null value
      real(double), parameter :: R_NULL = 0.0_double  ! Real null value

      ! functional dependence
      integer, parameter, public :: FD_DENSITY = 1                    ! density
      integer, parameter, public :: FD_ORBITAL = 2                    ! (OEP) orbital
      integer, parameter, public :: FD_HYBRID  = 3                    ! density and orbital

      ! density-dependent functional source
      integer, parameter, public :: DDFS_NATIVE = 1                    ! Native functionals
      integer, parameter, public :: DDFS_LIBXC  = 2                    ! Libxc functionals

      ! native: density-dependent exchange functionals
      integer, parameter, public :: E_LDA   =  1
      integer, parameter, public :: E_PW91  =  2
      integer, parameter, public :: E_PBE   =  3
      integer, parameter, public :: E_BLYP  =  4
      integer, parameter, public :: E_AM05  =  5                      ! R. Armiento, A. E. Mattsson, PRB 72, 085108 (2005)
      integer, parameter, public :: E_LSDA  = 21

      ! native: density-dependent correlation functionals
      integer, parameter, public :: C_PZ   = 1                        ! PRB 23, 5048 (1981)
      integer, parameter, public :: C_PW   = 2                        ! PRB 45, 13244 (1992)
      integer, parameter, public :: C_LYP  = 3                        ! PRB 37, 785 (1988)
      integer, parameter, public :: C_NG   = 4                        ! not given

      ! native/libxc: method used to compute functional derivatives
      integer, parameter, public :: XCD_ANALYTICAL = 1                ! analytical derivatives
      integer, parameter, public :: XCD_NUMERICAL  = 2                ! numerical derivatives

      ! native/libxc: method used to compute the potential
      integer, parameter, public :: XCD_SIMPLE     = 1                ! simple method (LDA functionals)
      integer, parameter, public :: XCD_WHITE_BIRD = 2                ! White and Bird method (GGA functionals)

      ! EXX: Coulomb kernel
      integer, parameter, public :: CK_NORMAL     = 1                 ! standard Coulomb kernel
      integer, parameter, public :: CK_ATTENUATED = 2                 ! attenuated Coulomb kernel
      integer, parameter, public :: CK_SCREENED   = 3                 ! screened (HSE) Coulomb kernel

      ! EXX: auxiliary type
      integer, parameter, public :: AT_LEGACY                = 1      ! no treatment of integrable divergence (gamma point only)
      integer, parameter, public :: AT_STRUCTURE_DEPENDENT   = 2      ! structure-dependent treatment of integrable divergence
      integer, parameter, public :: AT_STRUCTURE_INDEPENDENT = 3      ! structure-independent treatment of integrable divergence

      ! EXX: structure-dependent form (lattice type)
      integer, parameter, public :: SDF_SC  = 1                       ! simple cubic lattice
      integer, parameter, public :: SDF_BCC = 2                       ! body-centered cubic lattice
      integer, parameter, public :: SDF_FCC = 3                       ! face-centered cubic lattice

      ! EXX: structure-independent formulation
      integer, parameter, public :: SIF_SFH = 1                       ! Sorouri, Foulkes, Hine
      integer, parameter, public :: SIF_CRG = 2                       ! Carrier, Rohra, Gorling

      type :: xc_type_rep
        integer :: ref                    ! reference count
        type(ghost) :: g                  ! ghost
        integer :: dependence             ! dependence of the exchange-correlation functionals
        integer :: ddf_source             ! source of the density-dependent functionals
        integer :: functional             ! native/libxc: exchange-correlation form
        integer :: exchange               ! native/libxc: exchange form
        integer :: correlation            ! native/libxc: correlation form
        logical :: uses_gradient          ! native/libxc: indicates whether or not the density gradient is used
        logical :: uses_laplacian         ! native/libxc: indicates whether or not the density laplacian is used
        logical :: uses_nlcc              ! hybrid: applies non linear core correction to hybrid NCP calculation
        integer :: derivative_method      ! native/libxc: method used to compute density derivitives
        integer :: potential_method       ! native/libxc: method used to compute the (density-dependent) potential
        integer :: coulomb_kernel         ! EXX: Coulomb kernel type
        integer :: auxiliary_type         ! EXX: type of the auxiliary equation for k-point integration
        integer :: sd_aux_form            ! EXX: structure-dependent auxiliary form
        integer :: si_aux_form            ! EXX: structure-independent auxiliary form
        real(double) :: hybrid_mixing     ! hybrid mixing parameter (the fraction of EXX)
        real(double) :: omega_orb         ! range separation parameter for screened EXX (Hartree Fock)
        real(double) :: omega_den         ! range separation parameter for semilocal exchange
      end type

      type, public :: xc_type_obj
        private
        integer :: ref
        type(xc_type_rep), pointer :: o
      end type

      public :: xc_type
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_functional_dependence
      public :: x_ddf_source 
      public :: x_functional_type
      public :: x_exchange_type
      public :: x_correlation_type  
      public :: uses_gradient
      public :: uses_laplacian
      public :: uses_nlcc
      public :: x_derivative_method  
      public :: x_potential_method
      public :: x_coulomb_kernel
      public :: x_auxiliary_type
      public :: x_sd_aux_form
      public :: x_si_aux_form
      public :: x_hybrid_mixing
      public :: x_omega_orb
      public :: x_omega_den
      public :: diary
      public :: write_restart

      interface xc_type
        module procedure constructor_xct
      end interface
      interface my
        module procedure my_xct, my_new_xct
      end interface
      interface thy
        module procedure thy_xct
      end interface
      interface glean
        module procedure glean_xct
      end interface
      interface bequeath
        module procedure bequeath_xct
      end interface
      interface assignment(=)
        module procedure assign_xct
      end interface
      interface x_ref
        module procedure xct_ref
      end interface
      interface x_ghost
        module procedure xct_ghost
      end interface
      interface x_functional_dependence
        module procedure xct_functional_dependence
      end interface
      interface x_ddf_source
        module procedure xct_ddf_source
      end interface
      interface x_functional_type
        module procedure xct_functional_type
      end interface
      interface x_exchange_type
        module procedure xct_exchange_type
      end interface
      interface x_correlation_type
        module procedure xct_correlation_type
      end interface
      interface uses_gradient
        module procedure uses_gradient_xct
      end interface
      interface uses_laplacian
        module procedure uses_laplacian_xct
      end interface
      interface uses_nlcc
        module procedure uses_nlcc_xct
      end interface
      interface x_derivative_method
        module procedure xct_derivative_method
      end interface
      interface x_potential_method
        module procedure xct_potential_method
      end interface
      interface x_coulomb_kernel
        module procedure xct_coulomb_kernel
      end interface
      interface x_auxiliary_type
        module procedure xct_auxiliary_type
      end interface
      interface x_sd_aux_form
        module procedure xct_sd_aux_form
      end interface
      interface x_si_aux_form
        module procedure xct_si_aux_form
      end interface
      interface x_hybrid_mixing
        module procedure xct_hybrid_mixing
      end interface
      interface x_omega_orb
        module procedure xct_omega_orb
      end interface
      interface x_omega_den
        module procedure xct_omega_den
      end interface
      interface diary
        module procedure diary_xct
      end interface
      interface write_restart
        module procedure write_restart_xct
      end interface

      contains

! public routines

      function constructor_xct(restf) result(xct)
!doc$ function xc_type(restf) result(xct)
        type(tagio_obj), optional :: restf
        type(xc_type_obj) :: xct
!       effects: Constructs a new xct.

!cod$ 
        logical :: found, found_omega_orb, found_omega_den, have_defaults
        character(1) :: tios
        character(line_len) :: tag
        integer :: spin_status
        integer(long) :: dsize, iosl, ndata, s4
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info
        real(double) :: omega, alpha, beta

        if (present(restf)) call my(restf)

        xct%ref = 0
        allocate( xct%o )
        xct%o%ref = 0
        xct%o%g = x_ghost()

        ! open the XC_TYPE block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"XC_TYPE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: XC_TYPE block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)
        end if

        ! functional dependence
        call arglc("xctype_dependence",tag,found)
        if (found) then
          select case(trim(tag))
          case ("density","d")
            xct%o%dependence = FD_DENSITY
          case ("orbital","o")
            xct%o%dependence = FD_ORBITAL
          case ("hybrid","h")
            xct%o%dependence = FD_HYBRID
          case default
            if (error(.true.,"ERROR: xctype_dependence was not recognized")) goto 100
          end select
        else
          if (present(restf)) then
            if (i_access(restf)) tios = findfirsttag(restf,"DEPENDENCE")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (tios == TAG_NORMAL) then
              if (i_access(restf)) then
                dsize = sizeof_long
                ndata = 1
                call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                xct%o%dependence = s4
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%dependence)
            else
              call warn("WARNING: Legacy restart file - taking functional dependence to be density")
              xct%o%dependence = FD_DENSITY
            end if
          else
            xct%o%dependence = FD_DENSITY
          end if
        end if

        ! density-dependent-functional source
        select case (xct%o%dependence)
        case (FD_DENSITY)
          call arglc("xctype_ddf_source",tag,found)
          if (found) then
            select case(trim(tag))
            case ("native","n")
              xct%o%ddf_source = DDFS_NATIVE
            case ("libxc","l")
              xct%o%ddf_source = DDFS_LIBXC
            case default
              if (error(.true.,"ERROR: xctype_ddf_source was not recognized")) goto 100
            end select
          else
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"DDF_SOURCE")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (tios == TAG_NORMAL) then
                if (i_access(restf)) then
                  dsize = sizeof_long
                  ndata = 1
                  call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  xct%o%ddf_source = s4
                end if
                if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%ddf_source)
              else
                call warn("WARNING: Legacy restart file - taking density-dependent functional source to be native")
                xct%o%ddf_source = DDFS_NATIVE
              end if
            else
              xct%o%ddf_source = DDFS_NATIVE
            end if
          end if
        case (FD_ORBITAL)
          xct%o%ddf_source = I_NULL
        case (FD_HYBRID)
          call arglc("xctype_ddf_source",tag,found)
          if (found) then
            select case(trim(tag))
            case ("native","n")
              if (error(.true.,"ERROR: native xctype_ddf_source is not valid with hybrid functionals")) goto 100
            case ("libxc","l")
              xct%o%ddf_source = DDFS_LIBXC
            case default
              if (error(.true.,"ERROR: xctype_ddf_source was not recognized")) goto 100
            end select
          else
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"DDF_SOURCE")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR: DDF_SOURCE tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_long
                ndata = 1
                call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                xct%o%ddf_source = s4
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%ddf_source)
            else
              xct%o%ddf_source = DDFS_LIBXC
            end if
          end if
        end select

        ! spin status
        select case (xct%o%ddf_source)
        case (DDFS_LIBXC)
          select case (mpi_nsgroups())
          case (1)
            spin_status = XC_UNPOLARIZED
          case (2)
            spin_status = XC_POLARIZED
          end select
        end select

        ! null values
        xct%o%functional        = I_NULL
        xct%o%exchange          = I_NULL
        xct%o%correlation       = I_NULL
        xct%o%uses_gradient     = L_NULL
        xct%o%uses_laplacian    = L_NULL
        xct%o%uses_nlcc         = L_NULL
        xct%o%derivative_method = I_NULL
        xct%o%potential_method  = I_NULL
        xct%o%coulomb_kernel    = I_NULL
        xct%o%auxiliary_type    = I_NULL
        xct%o%sd_aux_form       = I_NULL
        xct%o%si_aux_form       = I_NULL
        xct%o%hybrid_mixing     = R_NULL
        xct%o%omega_orb         = R_NULL
        xct%o%omega_den         = R_NULL

        ! determine parameters for the density-dependent functionals
        select case (xct%o%ddf_source)
        case (DDFS_NATIVE)

          ! functional, exchange and correlation
          call arglc("xctype_functional",tag,found)
          if (found) then
            xct%o%exchange = exchange_i(tag) ; if (error()) goto 100
            xct%o%correlation = correlation_i(xct%o%exchange)
          else
            call arglc("xctype_exchange",tag,found)
            if (found) then
              xct%o%exchange = exchange_i(tag) ; if (error()) goto 100
              call arglc("xctype_correlation",tag,found)
              if (found) then
                xct%o%correlation = correlation_i(xct%o%exchange,tag) ; if (error()) goto 100
              else
                xct%o%correlation = C_NG
              end if
            else
              call warn("WARNING: Looking for legacy functional/correlation tags")
              call arglc("functional",tag,found)
              if (found) then
                xct%o%exchange = exchange_i(tag) ; if (error()) goto 100
                call arglc("correlation",tag,found)
                if (found) then
                  xct%o%correlation = correlation_i(xct%o%exchange,tag) ; if (error()) goto 100
                else
                  xct%o%correlation = correlation_i(xct%o%exchange)
                end if
              else
                if (present(restf)) then
                  if (i_access(restf)) tios = findfirsttag(restf,"EXCHANGE")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (error(tios == TAG_NOT_FOUND,"ERROR: EXCHANGE tag was not found")) goto 100
                  if (i_access(restf)) then
                    dsize = sizeof_long
                    ndata = 1
                    call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                    xct%o%exchange = s4
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%exchange)
                  if (i_access(restf)) tios = findfirsttag(restf,"CORRELATION")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (error(tios == TAG_NOT_FOUND,"ERROR: CORRELATION tag was not found")) goto 100
                  if (i_access(restf)) then
                    dsize = sizeof_long
                    ndata = 1
                    call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                    xct%o%correlation = s4
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%correlation)
                else
                  if (error(.true.,"ERROR: No exchange-correlation information was found")) goto 100
                end if
              end if
            end if
          end if

          ! dependence on gradient and laplacian
          xct%o%uses_gradient = gradient_dependent_i(xct%o%exchange)
          xct%o%uses_laplacian = laplacian_dependent_i(xct%o%exchange)
           
          ! default method of computing functional derivatives
          xct%o%derivative_method = default_derivative_method_i(xct%o%exchange,xct%o%correlation)

        case (DDFS_LIBXC)

          ! functional, exchange and correlation
          call arg("xctype_functional",xct%o%functional,found)
          if (found) then
            if (error(.not.libxc_xc_is_valid_i(xct%o%functional),'ERROR: xctype_functional is not valid')) goto 100
            select case (xc_f90_family_from_id(xct%o%functional))
            case (XC_FAMILY_LDA,XC_FAMILY_GGA)
              if (error(xct%o%dependence == FD_HYBRID,"ERROR: xctype_functional is not valid for a hybrid calculation")) goto 100
            case (XC_FAMILY_HYB_GGA)
              if (error(xct%o%dependence == FD_DENSITY,"ERROR: xctype_functional is only valid for a hybrid calculation")) goto 100
            end select
          else
            call arg("xctype_exchange",xct%o%exchange,found)
            if (found) then
              if (error(.not.libxc_x_is_valid_i(xct%o%exchange),'ERROR: xctype_exchange is not valid')) goto 100
              select case (xc_f90_family_from_id(xct%o%exchange))
              case (XC_FAMILY_HYB_GGA)
                if (error(xct%o%dependence == FD_DENSITY,"ERROR: xctype_exchange is only valid for a hybrid calculation")) goto 100
              end select
              call arg("xctype_correlation",xct%o%correlation,found)
              if (found) then
                if (error(.not.libxc_c_is_valid_i(xct%o%correlation),'ERROR: xctype_correlation is not valid')) goto 100
              end if
              select case (xc_f90_family_from_id(xct%o%correlation))
              case (XC_FAMILY_HYB_GGA)
                if (error(.true.,"ERROR: not prepared for a xctype_correlation in XC_FAMILY_HYB_GGA")) goto 100
              end select
            else
              if (present(restf)) then
                if (i_access(restf)) tios = findfirsttag(restf,"FUNCTIONAL")
                if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                if (error(tios == TAG_NOT_FOUND,"ERROR: FUNCTIONAL tag was not found")) goto 100
                if (i_access(restf)) then
                  dsize = sizeof_long
                  ndata = 1
                  call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  xct%o%functional = s4
                end if
                if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%functional)
                if (i_access(restf)) tios = findfirsttag(restf,"EXCHANGE")
                if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                if (error(tios == TAG_NOT_FOUND,"ERROR: EXCHANGE tag was not found")) goto 100
                if (i_access(restf)) then
                  dsize = sizeof_long
                  ndata = 1
                  call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  xct%o%exchange = s4
                end if
                if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%exchange)
                if (i_access(restf)) tios = findfirsttag(restf,"CORRELATION")
                if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                if (error(tios == TAG_NOT_FOUND,"ERROR: CORRELATION tag was not found")) goto 100
                if (i_access(restf)) then
                  dsize = sizeof_long
                  ndata = 1
                  call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  xct%o%correlation = s4
                end if
                if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%correlation)
              else
                if (error(.true.,"ERROR: exchange-correlation information was not found")) goto 100
              end if
            end if
          end if

          ! dependence on gradient and laplacian
          xct%o%uses_gradient = libxc_gradient_dependent_i(xct%o%functional,xct%o%exchange,xct%o%correlation)
          xct%o%uses_laplacian = libxc_laplacian_dependent_i(xct%o%functional,xct%o%exchange,xct%o%correlation)

          ! default method of computing functional derivatives
          xct%o%derivative_method = libxc_derivative_method_i(xct%o%functional,xct%o%exchange,xct%o%correlation)
 
        end select

        ! method of computing (density-dependent) functional derivatives 
        select case (xct%o%ddf_source)
        case (DDFS_NATIVE)
          call arglc("xctype_derivative_method",tag,found)
          if (found) then
            select case (trim(tag))
            case ("analytical","a")
              if (xct%o%derivative_method == XCD_NUMERICAL) then
                call warn("WARNING: analytical derivatives are not available - using numerical derivatives")
              end if
            case ("numerical","n")
              if (xct%o%derivative_method == XCD_ANALYTICAL) then
                call warn("WARNING: numerical derivatives are being used when analytical derivatives are available")
                xct%o%derivative_method = XCD_NUMERICAL
              end if
            case default
               call warn("WARNING: xctype_derivative_method was not recognized - using default method")
            end select
          else
            call warn("WARNING: Looking for legacy derivatives tag")
            call arglc("derivatives",tag,found)
            if (found) then
              select case (trim(tag))
              case ("analytical","a")
                if (xct%o%derivative_method == XCD_NUMERICAL) then
                  call warn("WARNING: analytical derivatives are not available - using numerical derivatives")
                end if
              case ("numerical","n")
                if (xct%o%derivative_method == XCD_ANALYTICAL) then
                  call warn("WARNING: numerical derivatives are being used when analytical derivatives are available")
                  xct%o%derivative_method = XCD_NUMERICAL
                end if
              case default
                call warn("WARNING: derivative_method was not recognized - using default method")
              end select
            else
              if (present(restf)) then
                if (i_access(restf)) tios = findfirsttag(restf,"DERIVATIVE_METHOD")
                if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                if (tios == TAG_NORMAL) then
                  if (i_access(restf)) then
                    dsize = sizeof_long
                    ndata = 1
                    call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                    xct%o%derivative_method = s4
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%derivative_method)
                else
                  call warn("WARNING: Legacy restart file - using the default derivative_method")
                end if
              else
                call warn("WARNING: derivative_method was not found - using default method")
              end if
            end if
          end if
        case (DDFS_LIBXC)
          call arglc("xctype_derivative_method",tag,found)
          if (found) then
            select case (trim(tag))
            case ("analytical","a")
              if (xct%o%derivative_method == XCD_NUMERICAL) then
                call warn("WARNING: analytical derivatives are not available - using numerical derivatives")
              end if
            case ("numerical","n")
              if (xct%o%derivative_method == XCD_ANALYTICAL) then
                call warn("WARNING: numerical derivatives are being used when analytical derivatives are available")
                xct%o%derivative_method = XCD_NUMERICAL
              end if
            case default
               call warn("WARNING: xctype_derivative_method was not recognized - using default method")
            end select
          else
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"DERIVATIVE_METHOD")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR: DERIVATIVE_METHOD tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_long
                ndata = 1
                call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                xct%o%derivative_method = s4
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%derivative_method)
            else
              call warn("WARNING: derivative_method was not found - using default method")
            end if
          end if
        end select

        ! method of computing the potential when using density-dependent functionals
        select case (xct%o%ddf_source)
        case (DDFS_NATIVE,DDFS_LIBXC)
          xct%o%potential_method = potential_method_i(xct%o%dependence,xct%o%uses_gradient,xct%o%uses_laplacian)
        end select

        ! parameters pertaining to EXX calculations
        select case (xct%o%dependence)
        case (FD_ORBITAL,FD_HYBRID)

          ! Coulomb kernel type
          call arglc("exx_coulomb_kernel",tag,found)
          if (found) then
            select case (trim(tag))
            case ("normal")
              xct%o%coulomb_kernel = CK_NORMAL
            case ("attenuated")
              xct%o%coulomb_kernel = CK_ATTENUATED
            case ("screened","hse","sx")
              xct%o%coulomb_kernel = CK_SCREENED
            case default
              if (error(.true.,"ERROR: exx_coulomb_kernel was not recognized")) goto 100
            end select
          else
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"COULOMB_KERNEL")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR:COULOMB_KERNEL tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_long
                ndata = 1
                call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                xct%o%coulomb_kernel = s4
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%coulomb_kernel)
            else
              xct%o%coulomb_kernel = CK_NORMAL
            end if
          end if

          select case (xct%o%coulomb_kernel)
          case (CK_NORMAL)

            ! Determine the auxiliary type
            call arglc("exx_auxiliary_type",tag,found)
            if (found) then
              select case (trim(tag))
              case ("legacy")
                xct%o%auxiliary_type = AT_LEGACY
              case ("structure_dependent","sd")
                xct%o%auxiliary_type = AT_STRUCTURE_DEPENDENT
              case ("structure_independent","si")
                if (error(.true.,"ERROR: structure independent auxiliary_type is not currently supported")) goto 100
                xct%o%auxiliary_type = AT_STRUCTURE_INDEPENDENT
              case default
                if (error(.true.,"ERROR: exx_auxiliary_type was not recognized")) goto 100
              end select
            else
              if (present(restf)) then
                if (i_access(restf)) tios = findfirsttag(restf,"AUXILIARY_TYPE")
                if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                if (tios == TAG_NORMAL) then
                  if (i_access(restf)) then
                    dsize = sizeof_long
                    ndata = 1
                    call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                    xct%o%auxiliary_type = s4
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%auxiliary_type)
                end if
              else
                xct%o%auxiliary_type = AT_LEGACY
              end if
            end if

            select case (xct%o%auxiliary_type)
            case (AT_STRUCTURE_DEPENDENT)

              ! Determine the structure-dependent auxiliary form
              call arglc("exx_sd_aux_form",tag,found)
              if (found) then
                select case (trim(tag))
                case ("sc")
                  xct%o%sd_aux_form = SDF_SC
                case ("bcc")
                  xct%o%sd_aux_form = SDF_BCC
                case ("fcc")
                  xct%o%sd_aux_form = SDF_FCC
                case default
                  if (error(.true.,"ERROR: exx_sd_aux_form was not recognized")) goto 100
                end select
              else
                if (present(restf)) then
                  if (i_access(restf)) tios = findfirsttag(restf,"SD_AUX_FORM")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (tios == TAG_NORMAL) then
                    if (i_access(restf)) then
                      dsize = sizeof_long
                      ndata = 1
                      call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                      xct%o%sd_aux_form = s4
                    end if
                    if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%sd_aux_form)
                  end if
                else
                  if (error(.true.,"ERROR: exx_sd_aux_form was not found")) goto 100
                end if
              end if

            case (AT_STRUCTURE_INDEPENDENT)

              if (error(.true.,"ERROR: structure_independent auxiliary forms are not currently available")) goto 100

              ! Determine the structure-independent form
              call arglc("exx_si_aux_form",tag,found)
              if (found) then
                select case (trim(tag))
                case ("sfh")
                  xct%o%si_aux_form = SIF_SFH
                case ("crg")
                  xct%o%si_aux_form = SIF_CRG
                case default
                  if (error(.true.,"ERROR: exx_si_aux_form was not recognized")) goto 100
                end select
              else
                if (present(restf)) then
                  if (i_access(restf)) tios = findfirsttag(restf,"SI_AUX_FORM")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (tios == TAG_NORMAL) then
                    if (i_access(restf)) then
                      dsize = sizeof_long
                      ndata = 1
                      call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                      xct%o%si_aux_form = s4
                    end if
                    if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%si_aux_form)
                  end if
                else
                  xct%o%si_aux_form = SIF_SFH
                end if
              end if

            end select

          case (CK_SCREENED)
            
            ! Set default parameters for functionals with screened exchange for known functionals
            ! These parameters cannot be read from libxc 1.x.x. and earlier versions of libxc 2.x.x
            ! For later versions of libxc 2.x.x, these parameters will replaced by values read from libxc
            select case (xct%o%ddf_source)
            case (DDFS_LIBXC)

               have_defaults = .false.
               
               ! Functional definitions in libxc_funcs.f90
               select case(xct%o%functional)

               case(XC_HYB_GGA_XC_HSE06)
                  !** According to Krukau et al [JCP 125, 224106, (2006)], 0.11 is the best value for
                  !   the screening length for a variety of properties.

                  have_defaults = .true.
                  xct%o%omega_orb = 0.11_double
                  xct%o%omega_den = 0.11_double
                  xct%o%hybrid_mixing = 0.25_double

               case(XC_HYB_GGA_XC_HSE03)
                  !  For HSE03, there is confusion in the literature due to an error in the
                  !   original paper [Heyd et al JCP 118 8207 (2003)] where the authors 
                  !   (erroneously) reported the screening length parameters as:
                  !
                  !     omega^PBE = omega^HF = 0.15 
                  !
                  !   but as reported in their 2006 erratum [Heyd et al JCP 124 219906 (2006)],
                  !   the parameters actually used in that 2003 paper were:
                  !
                  !     omega^PBE = 0.15*cbrt(2)  (cbrt is the cube root)
                  !     omega^HF = 0.15/sqrt(2)  
                  
                  have_defaults = .true.
                  xct%o%omega_orb = 0.15_double/sqrt(2.0_double)
                  xct%o%omega_den = 0.15_double*(2.0_double**(1.0_double/3.0_double))
                  xct%o%hybrid_mixing = 0.25_double

               case(XC_HYB_GGA_XC_CAM_B3LYP)

                  have_defaults = .true.
                  xct%o%omega_orb = 0.33_double
                  xct%o%omega_den = 0.33_double
                  xct%o%hybrid_mixing = 0.25_double

               case(XC_HYB_GGA_XC_TUNED_CAM_B3LYP)

                  have_defaults = .true.
                  xct%o%omega_orb = 0.150
                  xct%o%omega_den = 0.150
                  xct%o%hybrid_mixing = 0.25_double

               end select


            end select

          end select

        end select 

        ! Determine parameters for hybrids and screened hybrids
        select case(xct%o%dependence)
        case (FD_HYBRID)
          select case (xct%o%ddf_source)
          case (DDFS_LIBXC)
            select case(xc_f90_family_from_id(xct%o%functional))
            case (XC_FAMILY_HYB_GGA)
              call xc_f90_func_init(xc_func,xc_info,xct%o%functional,spin_status)

              !  Beginning of code for libxc 1.x.x
              !  if (xct%o%coulomb_kernel == CK_SCREENED) then
              !    if (error(.not. have_defaults ,"ERROR: Screened exchange parameters are not known for this functional")) goto 100
              !  else
              !    call xc_f90_hyb_gga_exx_coef(xc_func,xct%o%hybrid_mixing)
              !  end if
              !  End of code for libxc 1.x.x

              !  Code for earlier versions of libxc 2.x.x
              !  if (xct%o%coulomb_kernel == CK_SCREENED) then
              !    if (error(.not. have_defaults ,"ERROR: Screened exchange parameters are not known for this functional")) goto 100
              !  else
              !    call xc_f90_hyb_exx_coef(xc_func,xct%o%hybrid_mixing)
              !  end if
              !  End of code for earlier versions of libxc 2.x.x

              !  Beginning of code for later versions of libxc 2.x.x
              call xc_f90_hyb_cam_coef(xc_func, omega, alpha, beta)
              xct%o%omega_orb = omega
              xct%o%omega_den = omega
              if ((alpha > 0.0_double) .and. .not. (beta > 0.0_double)) then
                if (error(xct%o%coulomb_kernel == CK_SCREENED,"ERROR: functional is inconsistent with screeened exchange")) goto 100
                xct%o%hybrid_mixing = alpha
              elseif ((beta > 0.0_double) .and. .not. (alpha > 0.0_double)) then 
                if (error(xct%o%coulomb_kernel == CK_NORMAL,"ERROR: functional requires screened exchange")) goto 100         
                if (error(xct%o%coulomb_kernel == CK_ATTENUATED,"ERROR: functional requires screened exchange")) goto 100
                xct%o%hybrid_mixing = beta
              elseif ((beta > 0.0_double) .and. (alpha > 0.0_double)) then
                if (error(.true.,"ERROR: functionals with both HF exchange and screened exchange are not implemented")) goto 100
              else
                if (error(.true.,"ERROR: mixing coefficents indicate functional is not a hybrid")) goto 100
              end if
              !  End of code for later versions of libxc 2.x.x

              call xc_f90_func_end(xc_func)
            end select

            select case(xc_f90_family_from_id(xct%o%exchange))
            case (XC_FAMILY_HYB_GGA)
              call xc_f90_func_init(xc_func,xc_info,xct%o%exchange,spin_status)

              !  Beginning of code for libxc 1.x.x
              !  if (xct%o%coulomb_kernel == CK_SCREENED) then
              !    if (error(.not. have_defaults ,"ERROR: Screened exchange parameters are not known for this functional")) goto 100
              !  else
              !    call xc_f90_hyb_gga_exx_coef(xc_func,xct%o%hybrid_mixing)
              !  end if
              !  End of code for libxc 1.x.x

              !  Code for earlier versions of libxc 2.x.x
              !  if (xct%o%coulomb_kernel == CK_SCREENED) then
              !    if (error(.not. have_defaults ,"ERROR: Screened exchange parameters are not known for this functional")) goto 100
              !  else
              !    call xc_f90_hyb_exx_coef(xc_func,xct%o%hybrid_mixing)
              !  end if
              !  End of code for earlier versions of libxc 2.x.x

              !  Beginning of code for later versions of libxc 2.x.x
              call xc_f90_hyb_cam_coef(xc_func, omega, alpha, beta)
              xct%o%omega_orb = omega
              xct%o%omega_den = omega
              if ((alpha > 0.0_double) .and. .not. (beta > 0.0_double)) then
                if (error(xct%o%coulomb_kernel == CK_SCREENED,"ERROR: functional is inconsistent with screeened exchange")) goto 100
                xct%o%hybrid_mixing = alpha
              elseif ((beta > 0.0_double) .and. .not. (alpha > 0.0_double)) then
                if (error(xct%o%coulomb_kernel == CK_NORMAL,"ERROR: functional requires screened exchange")) goto 100
                if (error(xct%o%coulomb_kernel == CK_ATTENUATED,"ERROR: functional requires screened exchange")) goto 100
                xct%o%hybrid_mixing = beta
              elseif ((beta > 0.0_double) .and. (alpha > 0.0_double)) then
                if (error(.true.,"ERROR: functionals with both HF exchange and screened exchange are not implemented")) goto 100
              else
                if (error(.true.,"ERROR: mixing coefficents indicate functional is not a hybrid")) goto 100
              end if
              !  End of code for later versions of libxc 2.x.x

              call xc_f90_func_end(xc_func)
            case (XC_FAMILY_LDA,XC_FAMILY_GGA)
              call arg("xctype_hybrid_mixing",xct%o%hybrid_mixing,found)
              if (found) then
                if ((xct%o%hybrid_mixing < 0.0_double) .or. (xct%o%hybrid_mixing > 1.0_double)) then
                  if (error(.true.,"ERROR: xctype_hybrid_mixing is not valid")) goto 100
                end if
              else
                if (present(restf)) then
                  if (i_access(restf)) tios = findfirsttag(restf,"HYBRID_MIXING")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (error(tios == TAG_NOT_FOUND,"ERROR: HYBRID_MIXING tag was not found")) goto 100
                  if (i_access(restf)) then
                    dsize = sizeof_double
                    ndata = 1
                    call readf(xct%o%hybrid_mixing,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,xct%o%hybrid_mixing)
                else
                  if (error(.not.found,"ERROR: hybrid-mixing information was not found")) goto 100
                end if
              end if
            end select

            !  For screened exchange, allow the the user to override the default values for the screening lengths
            if (xct%o%coulomb_kernel == CK_SCREENED) then
               call arg("xctype_omega_orb", xct%o%omega_orb, found_omega_orb)
               if (found_omega_orb) then
                  call warn("Warning: overriding default orbital screening length")
                  if (xct%o%omega_orb < 0.0_double) then
                     if (error(.true.,"ERROR: omega_orb, screening length < 0. Not valid")) goto 100
                  end if
               end if
               call arg("xctype_omega_den", xct%o%omega_den, found_omega_den)
               if (found_omega_den) then
                  call warn("Warning: overriding default density screening length")
                  if (xct%o%omega_den < 0.0_double) then
                     if (error(.true.,"ERROR: omega_den, screening length < 0. Not valid")) goto 100
                  end if
               end if
            end if

          end select

          ! Determine if one should apply non linear core corrections to the exchange 
          call arglc("xctype_nlcc",tag,found)
          if (found) then
              select case(trim(tag))
              case("t","true",".true.","y","yes","on")
                 xct%o%uses_nlcc = .true.
              case("f","false",".false.","n","no","off")
                 xct%o%uses_nlcc = .false.
              case default
                 if (error(.true.,"Unrecognized argument to xctype_nlcc in argvf")) goto 200
              end select
          else
              xct%o%uses_nlcc = .false.
          endif

        end select

        ! close the XC_TYPE block
100     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

200     if (present(restf)) call glean(thy(restf))

        if (error("Exit xc_type_mod::constructor_xc")) continue

      end function 

      subroutine my_xct(xct)
!doc$ subroutine my(xct)
        type(xc_type_obj) :: xct

!cod$
        xct%ref = xct%ref + 1
        xct%o%ref = xct%o%ref + 1
      end subroutine

      subroutine my_new_xct(xcti,xct)
!doc$ subroutine my(xcti,xct)
        type(xc_type_obj) :: xcti, xct

!cod$
        xct%ref = 1
        xct%o => xcti%o
        xct%o%ref = xct%o%ref + 1
      end subroutine
      
      function thy_xct(xct) result(xcto)
!doc$ function thy(xct) result(xcto)
        type(xc_type_obj) :: xct, xcto

!cod$
        xct%ref = xct%ref - 1
        xct%o%ref = xct%o%ref - 1
        xcto%ref = xct%ref
        xcto%o => xct%o
      end function

      subroutine glean_xct(xct)
!doc$ subroutine glean(xct)
        type(xc_type_obj) :: xct

!cod$
        if (xct%o%ref < 1) then
          deallocate( xct%o )
        end if
      end subroutine

      subroutine bequeath_xct(xct)
!doc$ subroutine bequeath(xct)
        type(xc_type_obj) :: xct

!cod$
        continue
      end subroutine

      subroutine assign_xct(xct,xct2)
!doc$ subroutine assign(xct,xct2)
        type(xc_type_obj), intent(inout) :: xct
        type(xc_type_obj), intent(in) :: xct2

!cod$
        type(xc_type_obj) :: xctt
        call my(xct2)
        xctt%o => xct%o
        xct%o%ref = xct%o%ref - xct%ref
        xct%o => xct2%o
        xct%o%ref = xct%o%ref + xct%ref
        call glean(xctt)
        call glean(thy(xct2))
      end subroutine

      function xct_ref(xct) result(r)
!doc$ function x_ref(xct) result(r)
        type(xc_type_obj) :: xct
        integer, dimension(2) :: r
!       effects: Returns xct%ref and xct%o%ref.

!cod$
        r(1) = xct%ref
        r(2) = xct%o%ref
        call glean(xct)
      end function

      function xct_ghost(xct) result(g)
!doc$ function x_ghost(xct) result(g)
        type(xc_type_obj) :: xct
        type(ghost) :: g
!       effects: Returns the ghost of xct.

!cod$
        call my(xct)
        g = xct%o%g
        call glean(thy(xct))
      end function

      function xct_functional_dependence(xct) result(fd)
!doc$ function x_functional_dependence(xct) result(fd)
        type(xc_type_obj) :: xct
        integer :: fd
!       effects: Returns the functional dependence.
!cod$
        call my(xct)
        fd = xct%o%dependence
        call glean(thy(xct))
      end function

      function xct_ddf_source(xct) result(s)
!doc$ function x_ddf_source(xct) result(s)
        type(xc_type_obj) :: xct
        integer :: s
!       effects: Returns the source of the density-dependent functional.
!cod$
        call my(xct)
        s = xct%o%ddf_source
        call glean(thy(xct))
      end function

      function xct_functional_type(xct) result(ft)
!doc$ function x_functional_type(xct) result(ft)
        type(xc_type_obj) :: xct
        integer :: ft
!       effects: Returns the functional type.
!cod$
        call my(xct)
        ft = xct%o%functional
        call glean(thy(xct))
      end function

      function xct_exchange_type(xct) result(et)
!doc$ function x_exchange_type(xct) result(et)
        type(xc_type_obj) :: xct
        integer :: et
!       effects: Returns the exchange type.
!cod$
        call my(xct)
        et = xct%o%exchange
        call glean(thy(xct))
      end function

      function xct_correlation_type(xct) result(ct)
!doc$ function x_correlation_type(xct) result(ct)
        type(xc_type_obj) :: xct
        integer :: ct
!       effects: Returns the correlation type.
!cod$
        call my(xct)
        ct = xct%o%correlation
        call glean(thy(xct))
      end function

      function uses_gradient_xct(xct) result(ug)
!doc$ function uses_gradient(xct) result(ug)
        type(xc_type_obj) :: xct
        logical :: ug
!       effects: Returns .true. iff the functional uses gradients of the density.

!cod$ 
        call my(xct)
        ug = xct%o%uses_gradient 
        call glean(thy(xct))
      end function

      function uses_laplacian_xct(xct) result(ul)
!doc$ function uses_laplacian(xct) result(ul)
        type(xc_type_obj) :: xct
        logical :: ul
!       effects: Returns .true. iff the functional uses laplacians of the density.

!cod$
        call my(xct)
        ul = xct%o%uses_laplacian
        call glean(thy(xct))
      end function

      function uses_nlcc_xct(xct) result(nlcc)
!doc$ function uses_nlcc(xct) result(ul)
        type(xc_type_obj) :: xct
        logical :: nlcc
!       effects: Returns .true. iff a non linear core correction is applied to the exchange term

!cod$
        call my(xct)
        nlcc = xct%o%uses_nlcc
        call glean(thy(xct))
      end function

      function xct_coulomb_kernel(xct) result(ck)
!doc$ function x_coulomb_kernel(xct) result(ck)
        type(xc_type_obj) :: xct
        integer :: ck
!       effects: Returns the coulomb_kernel.

!cod$
        call my(xct)
        ck = xct%o%coulomb_kernel
        call glean(thy(xct))
      end function

      function xct_auxiliary_type(xct) result(at)
!doc$ function x_auxiliary_type(xct) result(at)
        type(xc_type_obj) :: xct
        integer :: at
!       effects: Returns the auxiliary_type.

!cod$
        call my(xct)
        at = xct%o%auxiliary_type
        call glean(thy(xct))
      end function

      function xct_sd_aux_form(xct) result(saf)
!doc$ function x_sd_aux_form(xct) result(saf)
        type(xc_type_obj) :: xct
        integer :: saf
!       effects: Returns the sd_aux_form

!cod$
        call my(xct)
        saf = xct%o%sd_aux_form
        call glean(thy(xct))
      end function

      function xct_si_aux_form(xct) result(saf)
!doc$ function x_si_aux_form(xct) result(saf)
        type(xc_type_obj) :: xct
        integer :: saf
!       effects: Returns the si_aux_form.

!cod$
        call my(xct)
        saf = xct%o%si_aux_form
        call glean(thy(xct))
      end function

      function xct_derivative_method(xct) result(dm)
!doc$ function x_derivative_method(xct) result(dm)
        type(xc_type_obj) :: xct
        integer :: dm
!       effects: Returns the derivative method.

!cod$
        call my(xct)
        dm = xct%o%derivative_method
        call glean(thy(xct))
      end function

      function xct_potential_method(xct) result(pm)
!doc$ function x_potential_method(xct) result(pm)
        type(xc_type_obj) :: xct
        integer :: pm
!       effects: Returns the potential method.

!cod$
        call my(xct)
        pm = xct%o%potential_method
        call glean(thy(xct))
      end function

      function xct_hybrid_mixing(xct) result(hm)
!doc$ function x_hybrid_mixing(xct) result(hm)
        type(xc_type_obj) :: xct
        real(double) :: hm
!       effects: Returns the hybrid_mixing.

!cod$
        call my(xct)
        hm = xct%o%hybrid_mixing
        call glean(thy(xct))
      end function

      function xct_omega_orb(xct) result(omega)
!doc$ function x_omega_orb(xct) result(omega)
        type(xc_type_obj) :: xct
        real(double) :: omega
!       effects: Returns the screening length for screened EXX

!cod$
        call my(xct)
        omega = xct%o%omega_orb
        call glean(thy(xct))
      end function

      function xct_omega_den(xct) result(omega)
!doc$ function x_omega_den(xct) result(omega)
        type(xc_type_obj) :: xct
        real(double) :: omega
!       effects: Returns the screening length for screened EXX

!cod$
        call my(xct)
        omega = xct%o%omega_den
        call glean(thy(xct))
      end function

      subroutine diary_xct(xct)
!doc$ subroutine diary(xct)
        type(xc_type_obj) :: xct
!       modifies: Output stream
!       effects: Writes exchange-correlation information to the diary file.

!cod$ 
        character(line_len)     :: cs, es, str
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info
        integer                :: spin_status

        call my(xct)

        if (i_access( diaryfile() )) then

          select case(xct%o%dependence)
          case (FD_DENSITY)
            write(x_unit(diaryfile()),'(/,t4,"Density-dependent functional:")')
            select case(xct%o%ddf_source)
            case (DDFS_NATIVE)
              es = estring_i(xct%o)
              cs = cstring_i(xct%o)
              write(x_unit(diaryfile()),'(/,t6,"Exchange    (native): ",a)') trim(es)
              write(x_unit(diaryfile()),'(  t6,"Correlation (native): ",a)') trim(cs)
            case (DDFS_LIBXC)
              select case (mpi_nsgroups())
              case (1)
                spin_status = XC_UNPOLARIZED
              case (2)
                spin_status = XC_POLARIZED
              end select
              if (xct%o%functional > 0) then
                call xc_f90_func_init(xc_func,xc_info,xct%o%functional,spin_status)
                call xc_f90_info_name(xc_info,str)
                write(x_unit(diaryfile()),'(/,t6,"Exchange-Correlation (libxc): ",a)') trim(str)
                call xc_f90_func_end(xc_func)
              end if
              if (xct%o%exchange > 0) then
                call xc_f90_func_init(xc_func,xc_info,xct%o%exchange,spin_status)
                call xc_f90_info_name(xc_info, str)
                write(x_unit(diaryfile()),'(/,t6,"Exchange    (libxc): ",a)') trim(str)
                call xc_f90_func_end(xc_func)
              end if
              if (xct%o%correlation > 0) then
                call xc_f90_func_init(xc_func,xc_info,xct%o%correlation,spin_status)
                call xc_f90_info_name(xc_info, str)
                write(x_unit(diaryfile()),'(t6,"Correlation (libxc): ",a)') trim(str)
                call xc_f90_func_end(xc_func)
              end if
            end select
          case (FD_ORBITAL)
            write(x_unit(diaryfile()),'(/,t4,"Orbital-dependent functional:")')
            write(x_unit(diaryfile()),'(/,t6,"Exchange    (native): EXX")')
            write(x_unit(diaryfile()),'(  t6,"Correlation (native): None")')
          case (FD_HYBRID)
            write(x_unit(diaryfile()),'(/,t4,"Hybrid functional:")')
            write(x_unit(diaryfile()),'(/,t6,"Exchange (native): EXX fraction = ",f0.3)') xct%o%hybrid_mixing
            select case(xct%o%ddf_source)
            case (DDFS_LIBXC)
              select case (mpi_nsgroups())
              case (1)
                spin_status = XC_UNPOLARIZED
              case (2)
                spin_status = XC_POLARIZED
              end select
              if (xct%o%functional > 0) then
                call xc_f90_func_init(xc_func,xc_info,xct%o%functional,spin_status)
                call xc_f90_info_name(xc_info,str)
                write(x_unit(diaryfile()),'(t6,"Exchange-Correlation (libxc): ",a," fraction = ",f0.3)') &
                     & trim(str), (1.0_double - xct%o%hybrid_mixing)
                call xc_f90_func_end(xc_func)
              end if
              if (xct%o%exchange > 0) then
                call xc_f90_func_init(xc_func,xc_info,xct%o%exchange,spin_status)
                call xc_f90_info_name(xc_info, str)
                write(x_unit(diaryfile()),'(t6,"Exchange    (libxc): ",a," fraction = ",f0.3)') &
                     & trim(str), (1.0_double - xct%o%hybrid_mixing)
                call xc_f90_func_end(xc_func)
              end if
              if (xct%o%correlation > 0) then
                call xc_f90_func_init(xc_func,xc_info,xct%o%correlation,spin_status)
                call xc_f90_info_name(xc_info, str)
                write(x_unit(diaryfile()),'(t6,"Correlation (libxc): ",a)') trim(str)
                call xc_f90_func_end(xc_func)
              end if
            end select
          end select
           
          select case (xct%o%dependence)
          case (FD_DENSITY,FD_HYBRID)
            select case (xct%o%derivative_method)
            case (XCD_ANALYTICAL)
              write(x_unit(diaryfile()),'(/,t6,"Analytical expressions for the functional derivatives")')
            case (XCD_NUMERICAL)
              write(x_unit(diaryfile()),'(/,t6,"Numerical expressions for the functional derivatives")')
            end select
          end select

          select case (xct%o%coulomb_kernel)
          case (CK_NORMAL)
            write(x_unit(diaryfile()),'(/,t4,"Using a standard Coulomb kernel in EXX calculations:")')
            select case (xct%o%auxiliary_type)
            case (AT_LEGACY)
              write(x_unit(diaryfile()),'(t6,"- integrable divergences are ignored")')
            case (AT_STRUCTURE_DEPENDENT)
              write(x_unit(diaryfile()),'(t6,"- integrable divergences are treated with a WCB formula")')
            end select
          case (CK_ATTENUATED)
            write(x_unit(diaryfile()),'(/,t4,"Using an attenuated Coulomb kernel in EXX calculations")')
          case (CK_SCREENED)
            write(x_unit(diaryfile()),'(/,t4,"Using a screened Coulomb kernel in EXX calculations")')
            write(x_unit(diaryfile()),'(t6,"Orbital screening parameter, aka omega^HF  (a0^-1):   ",f0.3)') xct%o%omega_orb
            write(x_unit(diaryfile()),'(t6,"Semilocal screening parameter, aka omega^PBE (a0^-1): ",f0.3)') xct%o%omega_den
          end select

          select case (xct%o%dependence)
          case (FD_HYBRID)
            if (xct%o%uses_nlcc) then
              write(x_unit(diaryfile()),'(t4,"Applying nonlinear correction to EXX portion of exchange")')
            endif
          end select

        end if

        call glean(thy(xct))

      end subroutine

      subroutine write_restart_xct(xct,nrestf)
!doc$ subroutine write_restart(xct,nrestf)
        type(xc_type_obj) :: xct
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes xct information to nrestf.

!cod$      
        integer(long) :: dsize, iosl, ndata, s4

        call my(xct)
        call my(nrestf)

        if (i_access(nrestf)) then

          call startblock(nrestf,"XC_TYPE")

          ! (functional) dependence
          call writetag(nrestf,"DEPENDENCE")
          dsize = sizeof_long
          ndata = 1
          s4 = xct%o%dependence
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! density-dependent functional source
          call writetag(nrestf,"DDF_SOURCE")
          dsize = sizeof_long
          ndata = 1
          s4 = xct%o%ddf_source
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! functional (exchange-correlation) type
          call writetag(nrestf,"FUNCTIONAL")
          dsize = sizeof_long
          ndata = 1
          s4 = xct%o%functional
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! exchange type
          call writetag(nrestf,"EXCHANGE")
          dsize = sizeof_long
          ndata = 1
          s4 = xct%o%exchange
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! correlation type
          call writetag(nrestf,"CORRELATION")
          dsize = sizeof_long
          ndata = 1
          s4 = xct%o%correlation
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! derivative method
          call writetag(nrestf,"DERIVATIVE_METHOD")
          dsize = sizeof_long
          ndata = 1
          s4 = xct%o%derivative_method
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! Coulomb-kernel type
          call writetag(nrestf,"COULOMB_KERNEL")
          s4 =  xct%o%coulomb_kernel
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! auxiliary type
          call writetag(nrestf,"AUXILIARY_TYPE")
          s4 =  xct%o%auxiliary_type
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! structure-dependent auxiliary form
          call writetag(nrestf,"SD_AUX_FORM")
          s4 =  xct%o%sd_aux_form
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! structure-independent auxiliary form
          call writetag(nrestf,"SI_AUX_FORM")
          s4 =  xct%o%si_aux_form
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! hybrid mixing
          call writetag(nrestf,"HYBRID_MIXING")
          dsize = sizeof_double
          ndata = 1
          call writef(xct%o%hybrid_mixing,dsize,ndata,x_tagfd(nrestf),iosl)

          ! screening length for screened EXX
          call writetag(nrestf,"OMEGA_ORB")
          dsize = sizeof_double
          ndata = 1
          call writef(xct%o%omega_orb,dsize,ndata,x_tagfd(nrestf),iosl)

          ! screening length for semilocal exchange
          call writetag(nrestf,"OMEGA_DEN")
          dsize = sizeof_double
          ndata = 1
          call writef(xct%o%omega_den,dsize,ndata,x_tagfd(nrestf),iosl)

          call endblock(nrestf)

        end if

        call glean(thy(xct))
        call glean(thy(nrestf))

        if (error("Exit xc_type_mod::write_restart_xct")) continue

      end subroutine

! private routines

      subroutine own_i(xct)
        type(xc_type_obj) :: xct
        type(xc_type_obj) :: xct_t
        if (xct%ref < xct%o%ref) then
          allocate( xct_t%o )
          xct_t%o%ref               = 0
          xct_t%o%g                 = xct%o%g
          xct_t%o%dependence        = xct%o%dependence
          xct_t%o%ddf_source        = xct%o%ddf_source
          xct_t%o%functional        = xct%o%functional
          xct_t%o%exchange          = xct%o%exchange
          xct_t%o%correlation       = xct%o%correlation
          xct_t%o%uses_gradient     = xct%o%uses_gradient
          xct_t%o%uses_laplacian    = xct%o%uses_laplacian
          xct_t%o%uses_nlcc         = xct%o%uses_nlcc
          xct_t%o%derivative_method = xct%o%derivative_method
          xct_t%o%potential_method  = xct%o%potential_method
          xct_t%o%coulomb_kernel    = xct%o%coulomb_kernel
          xct_t%o%auxiliary_type    = xct%o%auxiliary_type
          xct_t%o%sd_aux_form       = xct%o%sd_aux_form
          xct_t%o%si_aux_form       = xct%o%si_aux_form
          xct_t%o%hybrid_mixing     = xct%o%hybrid_mixing
          xct_t%o%omega_orb         = xct%o%omega_orb
          xct_t%o%omega_den         = xct%o%omega_den
          xct%o%ref = xct%o%ref - xct%ref
          xct%o => xct_t%o
          xct%o%ref = xct%o%ref + xct%ref
        end if
      end subroutine

      function estring_i(xctr) result(es)
        type(xc_type_rep) :: xctr
        character(line_len) :: es

        select case (xctr%exchange)
        case (E_LDA)
          es = "LDA"
        case (E_PW91)
          es = "PW91"
        case (E_PBE)
          es = "PBE"
        case (E_BLYP)
          es = "BLYP"
        case (E_AM05)
          es = "AM05"
        case (E_LSDA)
          es = "LSDA"
        end select

      end function

      function cstring_i(xctr) result(cs)
        type(xc_type_rep) :: xctr
        character(line_len) :: cs

        select case (xctr%correlation)
        case (C_PZ)
          cs = "PZ"
        case (C_PW)
          cs = "PW"
        case (C_LYP)
          cs = "LYP"
        case (C_NG)
          cs = "none"
        end select

      end function

      function exchange_i(tag) result(f)
        character(line_len), intent(in) :: tag
        integer :: f

        select case (mpi_nsgroups())
        case (1)
          select case (trim(tag))
          case ("lda")
            f = E_LDA
          case ("pw91")
            f = E_PW91
          case ("pbe")
            f = E_PBE
          case ("blyp")
            f = E_BLYP
          case ("am05")
            f = E_AM05
          case default
            if (error(.true.,"ERROR: unknown exchange type")) continue
          end select
        case (2)
          select case (trim(tag))
          case ("lsda")
            f = E_LSDA
          case default
            if (error(.true.,"ERROR: unknown exchange type")) continue
          end select
        end select

        if (error("Exit xc_type_mod::exchange_i")) continue

      end function

      function correlation_i(f,tag) result(c)
        integer, intent(in) :: f
        character(line_len), intent(in), optional :: tag
        integer :: c

        if (present(tag)) then
          select case (f)
          case (E_LDA)
            c = C_PZ
            select case (trim(tag))
            case ("pz")
              continue
            case ("pw")
              c = C_PW
            case ("lyp")
              call warn("WARNING: improper correlation - using default")
            case default
              call warn("WARNING: unknown correlation - using default")
            end select
          case (E_PW91)
            c = C_PW
            select case (trim(tag))
            case ("pw")
              continue
            case ("pz")
              call warn("WARNING: non-standard correlation being used")
              c = C_PZ
            case ("lyp")
              call warn("WARNING: improper correlation - using default")
            case default
              call warn("WARNING: unknown correlation - using default")
            end select
          case (E_PBE)
            c = C_PW
            select case (trim(tag))
            case ("pw")
              continue
            case ("pz")
              call warn("WARNING: non-standard correlation being used")
              c = C_PZ
            case ("lyp")
              call warn("WARNING: improper correlation - using default")
            case default
              call warn("WARNING: unknown correlation - using default")
            end select
          case (E_BLYP)
            c = C_LYP
            select case (trim(tag))
            case ("lyp")
              continue
            case ("pz","pw")
              call warn("WARNING: improper correlation - using default")
            case default
              call warn("WARNING: unknown correlation - using default")
            end select
          case (E_AM05)
            c = C_PW
            select case (trim(tag))
            case ("pw")
              continue
            case ("pz")
              call warn("WARNING: non-standard correlation being used")
              c = C_PZ
            case ("lyp")
              call warn("WARNING: improper correlation - using default")
            case default
              call warn("WARNING: unknown correlation - using default")
            end select
          case (E_LSDA)
            c = C_PW
            select case (trim(tag))
            case ("pw")
              continue
            case ("pz","lyp")
              call warn("WARNING: improper correlation - using default")
            case default
              call warn("WARNING: unknown correlation - using default")
            end select
          end select
        else
          select case (f)
          case (E_LDA)
            c = C_PZ
          case (E_PW91)
            c = C_PW
          case (E_PBE)
            c = C_PW
          case (E_BLYP)
            c = C_LYP
          case (E_AM05)
            c = C_PW
          case (E_LSDA)
            c = C_PW
          end select
        end if

      end function

      function libxc_xc_is_valid_i(xc) result(valid)
        integer :: xc
        logical :: valid

        integer :: xc_kind, xc_family, spin_status
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        select case (mpi_nsgroups())
        case (1)
          spin_status = XC_UNPOLARIZED
        case (2)
          spin_status = XC_POLARIZED
        end select
              
        call xc_f90_func_init(xc_func,xc_info,xc,spin_status)
        
        xc_kind = xc_f90_info_kind(xc_info)
        if (xc_kind == XC_EXCHANGE_CORRELATION) then
          valid = .true.
        else
          valid = .false.
        end if

        xc_family = xc_f90_info_family(xc_info)
        select case(xc_family)
        case (XC_FAMILY_UNKNOWN)
          valid = .false.
          call warn("WARNING: xc - XC_FAMILY UNKNOWN")
        case (XC_FAMILY_NONE)
          valid = .false.
          call warn("WARNING: xc - XC_FAMILY_NONE")
        end select

        call xc_f90_func_end(xc_func)

      end function

      function libxc_x_is_valid_i(x) result(valid)
        integer :: x
        logical :: valid

        integer :: x_kind, x_family, spin_status
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        select case (mpi_nsgroups())
        case (1)
          spin_status = XC_UNPOLARIZED
        case (2)
          spin_status = XC_POLARIZED
        end select
              
        call xc_f90_func_init(xc_func,xc_info,x,spin_status)
        
        x_kind = xc_f90_info_kind(xc_info)
        if (x_kind == XC_EXCHANGE) then
          valid = .true.
        else
          valid = .false.
        end if

        x_family = xc_f90_info_family(xc_info)
        select case(x_family)
        case (XC_FAMILY_UNKNOWN)
          valid = .false.
          call warn("WARNING: x - XC_FAMILY UNKNOWN")
        case (XC_FAMILY_NONE)
          valid = .false.
          call warn("WARNING: x - XC_FAMILY NONE")
        end select

        call xc_f90_func_end(xc_func)

      end function

      function libxc_c_is_valid_i(c) result(valid)
        integer :: c
        logical :: valid

        integer :: c_kind, c_family, spin_status
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        select case (mpi_nsgroups())
        case (1)
          spin_status = XC_UNPOLARIZED
        case (2)
          spin_status = XC_POLARIZED
        end select
              
        call xc_f90_func_init(xc_func,xc_info,c,spin_status)

        c_kind = xc_f90_info_kind(xc_info)
        if (c_kind == XC_CORRELATION) then
          valid = .true.
        else
          valid = .false.
        end if

        c_family = xc_f90_info_family(xc_info)
        select case(c_family)
        case (XC_FAMILY_UNKNOWN)
          valid = .false.
          call warn("WARNING: c - XC_FAMILY UNKNOWN")
        case (XC_FAMILY_NONE)
          valid = .false.
          call warn("WARNING: c - XC_FAMILY NONE")
        end select

        call xc_f90_func_end(xc_func)

      end function       

      function gradient_dependent_i(f) result(gd)
        integer, intent(in) :: f
        logical :: gd

        select case (f)
        case (E_LDA)
          gd = .false.
        case (E_PW91)
          gd = .true.
        case (E_PBE)
          gd = .true.
        case (E_BLYP)
          gd = .true.
        case (E_AM05)
          gd = .true.
        case (E_LSDA)
          gd = .false.
        end select

      end function

      function laplacian_dependent_i(f) result(ld)
        integer, intent(in) :: f
        logical :: ld

        select case (f)
        case (E_LDA)
          ld = .false.
        case (E_PW91)
          ld = .false.
        case (E_PBE)
          ld = .false.
        case (E_BLYP)
          ld = .true.
        case (E_AM05)
          ld = .false.
        case (E_LSDA)
          ld = .false.
        end select

      end function

      function potential_method_i(fd,ug,ul) result(pm)
        integer, intent(in) :: fd
        logical :: ug, ul
        integer :: pm

        select case (fd)
        case (FD_DENSITY,FD_HYBRID)
          pm = XCD_SIMPLE
          if (ug .or. ul) pm = XCD_WHITE_BIRD
        case (FD_ORBITAL)
          pm = I_NULL
        end select

      end function

      function default_derivative_method_i(e,c) result(dm)
        integer, intent(in) :: e, c
        integer :: dm

        dm = XCD_NUMERICAL
        select case (e)
        case (E_LDA)
          select case (c)
          case (C_PZ,C_PW)
            dm = XCD_ANALYTICAL
          end select
        case (E_PW91)
          select case (c)
          case (C_PZ,C_PW)
            dm = XCD_ANALYTICAL
          end select
        case (E_PBE)
          select case (c)
          case (C_PZ,C_PW)
            dm = XCD_ANALYTICAL
          end select
        case (E_BLYP)
          dm = XCD_ANALYTICAL
        case (E_AM05)
          select case (c)
          case (C_PZ,C_PW)
            dm = XCD_ANALYTICAL
          end select
        case (E_LSDA)
          select case (c)
          case (C_PW)
            dm = XCD_ANALYTICAL
          end select
        end select

      end function

      function libxc_gradient_dependent_i(xc,x,c) result(gd)
        integer, intent(in) :: xc, x, c
        logical :: gd

        integer :: spin_status
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        gd = .false.

        select case (mpi_nsgroups())
        case (1)
          spin_status = XC_UNPOLARIZED
        case (2)
          spin_status = XC_POLARIZED
        end select

        ! gradient dependence of xc
        if (0 < xc) then
          call xc_f90_func_init(xc_func,xc_info,xc,spin_status)
          select case(xc_f90_info_family(xc_info))
          case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA,XC_FAMILY_MGGA)
            gd = .true.
          end select
          call xc_f90_func_end(xc_func)
        end if

        ! gradient dependence of x
        if (0 < x) then
          call xc_f90_func_init(xc_func,xc_info,x,spin_status)
          select case(xc_f90_info_family(xc_info))
          case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA,XC_FAMILY_MGGA)
            gd = .true.
          end select
          call xc_f90_func_end(xc_func)
        end  if

        ! gradient dependence of c
        if (0 < c) then
          call xc_f90_func_init(xc_func,xc_info,c,spin_status)
          select case(xc_f90_info_family(xc_info))
          case (XC_FAMILY_GGA,XC_FAMILY_HYB_GGA,XC_FAMILY_MGGA)
            gd = .true.
          end select
          call xc_f90_func_end(xc_func)
        end if

      end function

      function libxc_laplacian_dependent_i(xc,x,c) result(ld)
        integer, intent(in) :: xc, x, c
        logical :: ld

        integer :: spin_status
        type(xc_f90_pointer_t) :: xc_func
        type(xc_f90_pointer_t) :: xc_info

        ld = .false.

        select case (mpi_nsgroups())
        case (1)
          spin_status = XC_UNPOLARIZED
        case (2)
          spin_status = XC_POLARIZED
        end select

        ! laplacian dependence of xc
        if (0 < xc) then
          call xc_f90_func_init(xc_func,xc_info,xc,spin_status)
          select case(xc_f90_info_family(xc_info))
          case (XC_FAMILY_MGGA)
            ld = .true.
          end select
          call xc_f90_func_end(xc_func)
        end if

        ! laplacian dependence of x
        if (0 < x) then
          call xc_f90_func_init(xc_func,xc_info,x,spin_status)
          select case(xc_f90_info_family(xc_info))
          case (XC_FAMILY_MGGA)
            ld = .true.
          end select
          call xc_f90_func_end(xc_func)
        end if

        ! laplacian dependence of c
        if (0 < c) then
          call xc_f90_func_init(xc_func,xc_info,c,spin_status)
          select case(xc_f90_info_family(xc_info))
          case (XC_FAMILY_MGGA)
            ld = .true.
          end select
          call xc_f90_func_end(xc_func)
        end if

      end function

      function libxc_derivative_method_i(xc,x,c) result(dm)
        integer, intent(in) :: xc, x, c
        integer :: dm

        dm = XCD_ANALYTICAL

      end function

      end module

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module time_propagator_mod
!doc$ module time_propagator_mod

!     One datatype is available here: type(time_propagator_obj)

!     Time_propagator_mod encapsulates the time_propagator object
!     and is roughly independent of basis.  Tags are used to read
!     in various quantities
!
! tddft_method {msd2 | msd4 | msd6 | ch | sop | sop2 | sop4 | lanczos | sil }
!   Indicates the numerical scheme used to propagate the wavefunctions
!
!     msdX    Multistep differencing schemes of order X
!     ch      Chebychev expansion
!     sop     Split Operator Method (2nd Order)
!     sop2    Split Operator Method (2nd Order)
!     sop4    Split Operator Method (4th Order)
!     sil     Short Iterative Lanczos method
!     lanczos Short Iterative Lanczos method
!
! tddft_step_size dt  
!   dt should be entered in units of fs (10E-15 seconds)
!
! tddft_tol tol
!   tol is the tolerance for the error in the propagation of the wavefunctions 
!   during each time step.
!
! tddft_calc_j  {yes|y|true|on|no|n|false|off|
!   calculate the current density
!
! tddft_hartree_vecpot yes|on|true|no|off|false
!   Determines whether the full current density (j) should be calculated at each
!      time step. Default is "false". The current density is required in order to
!      calculate the hartree vector potential.

!
      use arg_mod
      use diary_mod
      use error_mod
      use ghost_mod
      use operators_mod
!      use io_mod
      use kind_mod
      use math_mod
      use mpi_mod
      use grid_mod
      use multibasis_mod
      use multivector_mod
      use tagio_mod
      use timing_mod

!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = KGROUP

      !** Available Time-Evolving Schemes **
      integer, parameter :: NONE  = 0
      integer, parameter :: MSD2  = 1
      integer, parameter :: MSD4  = 2
      integer, parameter :: MSD6  = 3
      integer, parameter :: CH    = 10
      integer, parameter :: SOP2  = 11
      integer, parameter :: SOP4  = 12
      integer, parameter :: SIL   = 13

      !** Conversion Constants
      real(double), parameter :: FS_2_ARU = 1.0_double/0.048377687
      real(double), parameter :: ARU_2_FS = 0.048377687

      !** Default values
      integer                 :: default_method = CH
      real(double), parameter :: default_dt = 1.0d-3    ! In units of fs
      real(double), parameter :: default_tol = 1.0d-15  ! tolerance for propgation methods
      
      !** Define time_propagator object
      type :: time_propagator_rep
         integer      :: ref
         type(ghost)  :: g

         ! State variables
         integer      :: num_hamiltonian_ops       ! Number of H|Psi> operations per td step
         integer      :: num_mvs                   ! Number of electron objects needed by the time evolution operator
         logical      :: is_mvs_allocated          ! flag that indicates whether the mvs objects have been allocated
         logical      :: is_tmpv_init              ! Flag that indicates if tmp_v mvec has been initialized
         type(multivector_obj)          :: tmp_v   ! cached multivector
         type(multivector_obj), pointer :: mvs(:)  ! Array of multivector objects

         ! The following parameters can be read in from argvf
         integer      :: method       ! Method used to expand the propagation operator (chebychev, split operator, etc..)
         real(double) :: dt           ! Time step (\delta t) of simulation
         real(double) :: tol          ! Error tolerance - only used for some of the methods (ch, sil)
         logical      :: need_current_density ! True if this calculation requires the current density.

      end type

      type, public :: time_propagator_obj
         private
         integer :: ref
         type(time_propagator_rep), pointer :: o
      end type

!doc$
      public :: time_propagator
      public :: update
      public :: propagate
      public :: x_num_hamiltonian_ops
      public :: my, thy, glean, bequeath, assignment(=)
      public :: write_restart
      public :: need_current_density
!cod$

      ! Define interface to dbesjn in case user wants to use an external call to a function that 
      !   returns a bessel function of the first kind for the chebychev propagator.
      interface
         real(8) function dbesjn(n,x)
           real(8) x
           integer(4) n
         end function
      end interface

      interface time_propagator
        module procedure constructor_tp
      end interface
      interface update
        module procedure update_tp
      end interface
      interface propagate
        module procedure propagate_tp
      end interface
      interface my
        module procedure my_tp, my_new_tp
      end interface
      interface thy
        module procedure thy_tp
      end interface
      interface glean
        module procedure glean_tp
      end interface
      interface bequeath
        module procedure bequeath_tp
      end interface
      interface assignment(=)
        module procedure assign_tp
      end interface
      interface x_num_hamiltonian_ops
         module procedure get_num_hamiltonian_ops
      end interface
      interface write_restart
         module procedure write_restart_tp
      end interface
      interface need_current_density
         module procedure need_current_density_tp
      end interface


!*********************************************************************************************************

      contains

      function constructor_tp(mb,restf) result(tp)
!doc$ function time_propagator() result(tp)
        type(tagio_obj), intent(inout), optional :: restf
        type(multibasis_obj), optional           :: mb
        type(time_propagator_obj)                :: tp
!       effects: Constructs a new time_propagator_obj object
!       errors: Unrecognized time propagation method. Parameters out of bounds.
!cod$
        ! local vars
!        integer :: num_remaining_steps
        !------------------------------------------------------------

        tp%ref = 0
        allocate( tp%o )
        tp%o%ref = 0
        tp%o%g = x_ghost()

        tp%o%num_hamiltonian_ops = 1

        !** Read in all parameters from argvf file
        call read_parameters_i(tp%o); if (error()) goto 100

        !** Check whether restart info is being read in from a file
        if (present(restf) .and. present(mb)) then
           call my(restf)
           call my(mb)
           call read_restart_i(tp%o,mb,restf); if (error()) goto 100
           call glean(thy(mb))
           call glean(thy(restf))
        else
           !** Initialize all the bookeeping vars
           tp%o%is_mvs_allocated = .false.   !** set the allocated flag to false because we can't allocate
                                             !   multivectors until the step routine is called with access
                                             !   to some sort of multibasis
        end if

        tp%o%is_tmpv_init = .false.  ! Always initialize to false - even if restarting

100     if (error("Exit time_propagator_mod::constructor_tp")) continue

      end function

!*********************************************************************************************************

      subroutine read_restart_i(tpr,mb,restf)
        type(time_propagator_rep)      :: tpr
        type(multibasis_obj)           :: mb
        type(tagio_obj), intent(inout) :: restf

        character(1)       :: tios
        character(line_len) :: usage
        integer            :: num_mvs, imvs
        integer(long)      :: s4, dsize, ios, ndata

        call my(mb)
        call my(restf)

        if (i_access(restf)) tios = findnexttag(restf,"TIME_PROPAGATOR")
        if (i_comm(restf)) call broadcast(MOD_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: TIME_PROPAGATOR block was not found")) goto 200

        if (i_access(restf)) call openblock(restf)

        if (i_access(restf)) tios = findnexttag(restf,"PARAMETERS")
        if (i_comm(restf)) call broadcast(MOD_SCOPE,tios)
        if(error(tios == tag_not_found,"ERROR: PARAMETERS tag not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_long ; ndata = 1
          call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
          num_mvs = s4
        end if
        if (i_comm(restf)) call broadcast(MOD_SCOPE,num_mvs)
        if (error(num_mvs /= tpr%num_mvs,"ERROR: num_mvs in restart file not consistent with params")) goto 100

        if (i_access(restf)) tios = findnexttag(restf,"MVS_HISTORY")
        if (i_comm(restf)) call broadcast(MOD_SCOPE,tios)
        if(error(tios == tag_not_found,"ERROR:MVS_HISTORY tag not found")) goto 100

        allocate( tpr%mvs(num_mvs) )
        tpr%is_mvs_allocated = .true.

        usage = "normal"
        do imvs = 1,num_mvs
          call my(multivector(mb,usage,restf=restf),tpr%mvs(imvs)) ; if (error()) exit
        end do

100     if (i_access(restf)) call closeblock(restf)

200     call glean(thy(mb))
        call glean(thy(restf))

        if (error("Exit tddft_mod::read_restart_i")) continue

      end subroutine

!*********************************************************************************************************

      subroutine read_parameters_i(tpr)
        type(time_propagator_rep)      :: tpr

        ! local vars
        character(line_len) :: tag
        logical            :: fnd, fnd_alt
        integer            :: j
        integer            :: status, i

        !** Read the parameter defining the method **
        call arglc("tddft_method",tag,fnd)
        if (.not.fnd) tag = "none"
        select case (tag(1:len_trim(tag)))
        case ("msd2")
           tpr%method = MSD2
           tpr%num_mvs = 2
        case ("msd4")
           tpr%method = MSD4
           tpr%num_mvs = 5
        case ("msd6")
           tpr%method = MSD6
           tpr%num_mvs = 7
        case ("ch")
           tpr%method = CH
           tpr%num_mvs = 4
        case ("sop","sop2")
           tpr%method = SOP2
           tpr%num_mvs = 1
        case ("sop4")
           tpr%method = SOP4
           tpr%num_mvs = 1
        case ("sil","lanczos")
           tpr%method = SIL
           tpr%num_mvs = 2
        case default
           if (error(.true.,"constructor_tp: unrecognized tddft_method")) goto 100
        end select

        allocate(tpr%mvs(tpr%num_mvs),STAT=status)
        if (error(status /= 0,"constructor_tp: could not allocate tpr%mvs()")) goto 100
        
        ! Size of each time step - N.B. the step size is read in assuming units of fs and
        !    then converted to atomic rydberg time units (4.84E-17s)
        call arg("tddft_step_size",tpr%dt,fnd)
        if (.not.fnd) tpr%dt = default_dt*FS_2_ARU
        if (error(tpr%dt<0,"ERROR! step_size < 0")) goto 100
        tpr%dt = tpr%dt*FS_2_ARU

!        !** Find the total number of time steps 
!        tpr%num_steps = int(tpr%runtime/tpr%dt)

        !** Get the tolerance for the error
        call arg("tddft_tol",tpr%tol,fnd)
        if (.not.fnd) tpr%tol = default_tol
        if (error(tpr%tol<0,"ERROR: tol < 0")) goto 100


        !** Read in whether the current density should be calculated
        call arglc("tddft_calc_j",tag,fnd)
        if (.not.fnd) tag = "false"
        select case (tag)
        case ("yes","y","on","true",".true.")
           tpr%need_current_density = .true.
        end select

        !** Read in whether the hartree vector potential is to be calculated.
        !     if so, the current density must be calculated. Note that if this 
        !     is set to true it overrides the tddft_calc_j option if it has 
        !     been set to false.
        call arglc("tddft_hartree_vecpot",tag,fnd)
        if (.not.fnd) tag = "false"
        select case (tag)
        case ("yes","y","on","true",".true.")
           tpr%need_current_density = .true.
        end select
        
100     if (error("Exit tddft_mod::read_parameters_i")) continue
        
      end subroutine


!*********************************************************************************************************

      subroutine update_tp(tp,max_dir,res_tol)
!doc$ subroutine update(tp,max_dir,res_tol)
        type(time_propagator_obj), intent(inout) :: tp
        integer, intent(in), optional :: max_dir
        real(double), intent(in), optional :: res_tol
!       modifies: tp
!       effects: nothing

!cod$
        continue

      end subroutine

!*********************************************************************************************************

      subroutine propagate_tp(tp,h,v,evals)
!doc$ subroutine propagate(tp,h,v,evals)
        type(time_propagator_obj), intent(inout)  :: tp
        type(hamiltonian_obj)                     :: h
        type(multivector_obj), intent(inout)      :: v
        real(double), dimension(:), intent(inout) :: evals
!       modifies: v
!       effects: Evolves the wave functions in v by a single time step.

!cod$
        ! Local vars
        logical :: debug
        integer :: i
        !-----------------------------------

        debug = .false.

        if (debug) call warn( "propagate starting ")


        if (debug) call warn( "propagate my... ")
        call my(tp)
        call my(h)
        call my(v)

        !** Initialize the cached mvec if necessary
        if (.not. tp%o%is_tmpv_init) then
           call my(v,tp%o%tmp_v)
           tp%o%is_tmpv_init = .true.
        end if

        if (tp%o%method == NONE) GOTO 100

        if (debug) call warn( "propagate own_i ")
        call own_i(tp)

        !** Allocate the temporary multivectors and fill them with the current wavefuncs in v
        if (.not.(tp%o%is_mvs_allocated)) then
           do i=1,tp%o%num_mvs
              call my(v,tp%o%mvs(i))
           end do
           tp%o%is_mvs_allocated = .true.
        end if

        !** Calculate the vector potential at this time
!        tp%o%A_prev = tp%o%A
!        tp%o%A = get_external_vecpot(tp, .true.)

        select case (tp%o%method)
        case (MSD2)
          call msd2_i(tp,h,v) ; if (error()) goto 100

        case (MSD4)
           call msd4_i(tp,h,v) ; if (error()) goto 100

        case (MSD6)
           call msd6_i(tp,h,v) ; if (error()) goto 100

        case (CH)
           call ch_i(tp,h,v) ; if (error()) goto 100

        case (SOP2)
           call sop2_i(tp,h,v) ; if (error()) goto 100

        case (SOP4)
           call sop4_i(tp,h,v) ; if (error()) goto 100

        case (SIL)
           call sil_i(tp,h,v) ; if (error()) goto 100

        case default
           call ch_i(tp,h,v) ; if (error()) goto 100
        end select

        if (debug) call warn( "propagate elapsed_time ")

        !** Update the expectation values of the Hamiltonian (necessary to update energy)
        if (x_ghost(x_multibasis(v)) /= x_ghost(x_multibasis(tp%o%tmp_v))) tp%o%tmp_v = v
        call apply_hamiltonian(v,tp%o%tmp_v,h) ; if (error()) goto 100  ! H*v --> hv
        call multiply(v,tp%o%tmp_v,evals)

100     call glean(thy(tp))
        call glean(thy(h))
        call glean(thy(v))

        !deallocate(s)

        if (error("Exit time_propagator_mod::propagate_tp")) continue

        if (debug) call warn( "propagate ending" )

      end subroutine


!*********************************************************************************************************

      !** Propagates using the 2nd order multistep differencing scheme.
      !    v(t+dt) = v(t-dt) - 2*i*dt*H*v(t)
      subroutine msd2_i(tp, h, v)
        type(time_propagator_obj)            :: tp
        type(hamiltonian_obj)                :: h
        type(multivector_obj), intent(inout) :: v

        ! local vars
        type(multivector_rep), pointer :: worm_v, worm_tmp, worm_hv, worm_prev1
        complex(double), pointer       :: tmp(:,:)
        complex(double)                :: i

        logical :: debug
        !-------------------------------------------------------------
        !-------------------------------------------------------------
        !-------------------------------------------------------------

        i = (0.0_double, 1.0_double)

        call my(tp)
        call my(h)
        call my(v)

        debug = .false.

        call apply_hamiltonian(v,tp%o%mvs(2),h) ; if (error()) GOTO 100

        worm_v => wormhole(v)
        worm_prev1 => wormhole(tp%o%mvs(1))
        worm_hv => wormhole(tp%o%mvs(2))

        worm_prev1%mat = worm_prev1%mat - 2*i*tp%o%dt*worm_hv%mat

        tmp => worm_prev1%mat
        worm_prev1%mat => worm_v%mat
        worm_v%mat => tmp

        nullify(worm_v, worm_hv, worm_prev1, tmp)

100     call glean(thy(tp))
        call glean(thy(h))
        call glean(thy(v))

      end subroutine




!*********************************************************************************************************

      !** Propagates using the 4th order multistep differencing scheme
      !    v(t+dt) = v(t-3dt) - 4*i*dt H {- (1/3)*v(t-dt)
      !                                   + (2/3)*[v(t) + v(t-2dt)]}
      subroutine msd4_i(tp, h, v)
        type(time_propagator_obj)            :: tp
        type(hamiltonian_obj)                :: h
        type(multivector_obj), intent(inout) :: v

        ! local vars
        type(multivector_rep), pointer :: worm_v, worm_hv, worm_tmp
        type(multivector_rep), pointer :: worm_prev1, worm_prev2, worm_prev3
        complex(double), pointer       :: tmp(:,:)
        complex(double)                :: i

        !-------------------------------------------------------------

        call my(tp)
        call my(h)
        call my(v)

        worm_v     => wormhole(v)
        worm_prev1 => wormhole(tp%o%mvs(1))
        worm_prev2 => wormhole(tp%o%mvs(2))
        worm_prev3 => wormhole(tp%o%mvs(3))
        worm_tmp   => wormhole(tp%o%mvs(4))
        worm_hv    => wormhole(tp%o%mvs(5))

        i = (0.0_double,1.0_double)

        ! Construct the sum in curly braces.
        worm_tmp%mat = -(1.0_double/3.0_double)*worm_prev1%mat &
                       +(2.0_double/3.0_double)*(worm_v%mat + worm_prev2%mat)

        call apply_hamiltonian(tp%o%mvs(4),tp%o%mvs(5),h) ; if (error()) GOTO 100

        worm_tmp%mat = worm_prev3%mat - 4.0_double*i*tp%o%dt*worm_hv%mat

        tmp => worm_prev3%mat
        worm_prev3%mat => worm_prev2%mat
        worm_prev2%mat => worm_prev1%mat
        worm_prev1%mat => worm_v%mat
        worm_v%mat => worm_tmp%mat
        worm_tmp%mat => tmp    !** store v(t-3dt) in worm_tmp

        nullify(worm_v,worm_hv,worm_prev1,worm_prev2,worm_prev3,worm_tmp,tmp)

100     call glean(thy(tp))
        call glean(thy(h))
        call glean(thy(v))

      end subroutine


!*********************************************************************************************************

      !** Propagates using the 6th order multistep differencing scheme
      !    v(t+dt) = v(t-5dt) - 6*i*dt H {  (13/10)* v(t-2dt) 
      !                                   -  (7/10)*[v(t-dt) + v(t-3dt)] 
      !                                   + (11/20)*[v(t)    + v(t-4dt)]}
      subroutine msd6_i(tp, h, v)
        type(time_propagator_obj)            :: tp
        type(hamiltonian_obj)                :: h
        type(multivector_obj), intent(inout) :: v

        ! local vars
        type(multivector_rep), pointer :: worm_v, worm_hv, worm_tmp
        type(multivector_rep), pointer :: worm_prev1, worm_prev2, worm_prev3, worm_prev4, worm_prev5
        complex(double), pointer       :: tmp(:,:)
        complex(double)                :: i

        !-------------------------------------------------------------

        call my(tp)
        call my(h)
        call my(v)

        worm_v     => wormhole(v)
        worm_prev1 => wormhole(tp%o%mvs(1))
        worm_prev2 => wormhole(tp%o%mvs(2))
        worm_prev3 => wormhole(tp%o%mvs(3))
        worm_prev4 => wormhole(tp%o%mvs(4))
        worm_prev5 => wormhole(tp%o%mvs(5))
        worm_tmp   => wormhole(tp%o%mvs(6))
        worm_hv    => wormhole(tp%o%mvs(7))

        i = (0.0_double,1.0_double)

        ! Construct the sum in curly braces.
        worm_tmp%mat =   (13.0_double/10.0_double)* worm_prev2%mat &
                       -  (7.0_double/10.0_double)*(worm_prev1%mat + worm_prev3%mat) &
                       + (11.0_double/20.0_double)*(worm_v%mat     + worm_prev4%mat) 

        call apply_hamiltonian(tp%o%mvs(6),tp%o%mvs(7),h) ; if (error()) GOTO 100

        worm_tmp%mat = worm_prev5%mat - 6.d0*(0.d0, 1.d0)*tp%o%dt*worm_hv%mat

        tmp => worm_prev5%mat
        worm_prev5%mat => worm_prev4%mat
        worm_prev4%mat => worm_prev3%mat
        worm_prev3%mat => worm_prev2%mat
        worm_prev2%mat => worm_prev1%mat
        worm_prev1%mat => worm_v%mat
        worm_v%mat => worm_tmp%mat
        worm_tmp%mat => tmp

        nullify(worm_prev1,worm_prev2,worm_prev3,worm_prev4,worm_prev5)
        nullify(worm_v,worm_hv,worm_tmp,tmp)

100     call glean(thy(tp))
        call glean(thy(h))
        call glean(thy(v))

      end subroutine

!*********************************************************************************************************
!*********************************************************************************************************

      !** ch_i
      !
      ! Propagate using the Chebyshev expansion of the time evolution operator
      !
      !                      __
      !                      \
      !    v(dt) = exp(-i*z) /_ a(x) phi(-i*Hnorm)v(0)
      !                       k  k      k
      !
      !           z = (Emax + Emin)*dt/2
      !
      !           x = (Emax - Emin)*dt/2
      !
      !           a(x) = J(x)   ,    a(x) = 2*J(x)   for k > 0
      !            0      0           k        k
      !
      !                       2*H
      !           H    =  ----------- - (z/x)*I
      !            norm   Emax - Emin
      !
      !
      !          phi   = -2*i*H * phi  +  phi
      !            k+1         norm  k      k-1
      !
      !                phi  = v(0)        phi  = -i*H * v(0)
      !                   0                  1       norm
      !
      subroutine ch_i(tp,h,v)
        type(time_propagator_obj) :: tp
        type(hamiltonian_obj)     :: h
        type(multivector_obj)     :: v

        ! local vars
        type(multivector_rep), pointer :: worm_v, worm_phi, worm_phi_minus, worm_phi_plus
        type(multivector_obj), pointer :: phi, phi_minus, phi_plus
        complex(double), pointer       :: tmp(:,:)
        complex(double), parameter     :: i = (0.0_double, 1.0_double)
        complex(double)                :: phase
        integer                        :: k,n
        real(double)                   :: dA, max_dt, emin, emax, egrid, esum, x, z
        real(double), pointer          :: evals(:)
        !-------------------------------------------------------------------------------------------

        call my(tp)
        call my(v)
        call my(h)

        !** Find the time step based on how much the potential changes
        dA = 0.5_double

        allocate(evals( x_n_bands(v) ))

        emin = lp_minimum(x_h_common(h))
        emax = lp_maximum(x_h_common(h))

!write(*,*) 'time_propagator::ch_i  emin, emax = ', emin, emax

        egrid = emax - emin
        esum  = emax + emin
        z =  esum*tp%o%dt/2.0_double
        x = egrid*tp%o%dt/2.0_double

        n = ceiling(x)

        !** upper limit on the number of iterations
        n = 40
!       N.M.
!       n = 24

        worm_v         => wormhole(v)
        worm_phi_minus => wormhole(tp%o%mvs(1))
        worm_phi       => wormhole(tp%o%mvs(2))
        worm_phi_plus  => wormhole(tp%o%mvs(3))

        !** phi = v(0)
        !      0
        worm_phi_minus%mat = worm_v%mat

        !** calculate the zeroth term
        k = 0
        worm_v%mat = besjn(k,x)*worm_phi_minus%mat
!        worm_v%mat = dbesjn(k,x)*worm_phi_minus%mat

        !** phi = -i*H   * v(0)
        !      1      norm
        call hnorm_i(tp%o%mvs(1),tp%o%mvs(2),h,emax,emin) ; if (error()) goto 100
        worm_phi%mat = -i*worm_phi%mat

        !** calculate the first term
        k = 1
        worm_v%mat = worm_v%mat + 2*besjn(k,x)*worm_phi%mat
!        worm_v%mat = worm_v%mat + 2*dbesjn(k,x)*worm_phi%mat

!        write(*,*) "   k, besj ", k, besjn(k,x)

        !** calculate the remaining terms
        do while (k <= n)
           k = k + 1

           call hnorm_i(tp%o%mvs(2),tp%o%mvs(3),h,emax,emin)

           worm_phi_plus%mat = -2.0_double*i*worm_phi_plus%mat + worm_phi_minus%mat

           worm_v%mat = worm_v%mat + 2.0_double*besjn(k,x)*worm_phi_plus%mat
!           worm_v%mat = worm_v%mat + 2.0_double*dbesjn(k,x)*worm_phi_plus%mat
!           write(*,*) "   k, besj ", k, besjn(k,x)

           tmp => worm_phi_minus%mat
           worm_phi_minus%mat => worm_phi%mat
           worm_phi%mat => worm_phi_plus%mat
           worm_phi_plus%mat => tmp

           !** set k so that the loop exits if the bessel function is small enough
           !    the tolerance is read in or set to some default at construction
!           if (dbesjn(k,x) < tp%o%tol) n = k-1
!          N.M.
           if (besjn(k,x) < tp%o%tol) n = k-1
        end do

        tp%o%num_hamiltonian_ops = k
!        write(*,*) "  J(x) = ", besjn(k,x)

        worm_v%mat = exp(-i*z)*worm_v%mat

100     call glean(thy(tp))
        call glean(thy(v))
        call glean(thy(h))

        deallocate(evals)

        if (error("Exit time_propagator_mod::ch_i")) continue

      end subroutine

!*******************************************************************************************************

      !** hnorm_i
      !
      !                    2*H    (emax + x)
      !           H    =  ----- - ---------- * I
      !            norm   Egrid   (emax - x)
      !
      !           z = (Emax + Emin)*dt/2
      !
      !           x = (Emax - Emin)*dt/2
      !
      !           Egrid = Emax - Emin
      !
      subroutine hnorm_i(v,hv,h,emax,emin)
        type(multivector_obj)     :: v
        type(multivector_obj)     :: hv
        type(hamiltonian_obj)     :: h
        real(double)              :: emax
        real(double)              :: emin

        !local vars
        real(double)              :: egrid, esum
        type(multivector_rep), pointer :: worm_v, worm_hv

        real(double), pointer     :: evals(:)

        !------------------------------------------

        call my(v)
        call my(hv)
        call my(h)

        allocate(evals( x_n_bands(v) ))

        egrid = emax - emin
        esum  = emax + emin

        worm_v => wormhole(v)

!        call apply_hamiltonian(v,hv,h,A0,A1) ; if (error()) GOTO 100
        call apply_hamiltonian(v,hv,h) ; if (error()) GOTO 100

        call multiply(v,hv,evals)
!        write(*,*) "   hnorm ==> <psi|H|psi> = ", evals(1), evals(2), evals(3), evals(4)

        worm_hv => wormhole(hv)

!        write(*,*) "  hnorm: emax, emin ", emax, emin
!        write(*,*) "  hnorm: 2/egrid, esum/egrid ", 2.0_double/egrid, esum/egrid

        worm_hv%mat = 2.0_double*worm_hv%mat/egrid - worm_v%mat*esum/egrid

        call multiply(v,hv,evals)
!        write(*,*) "   hnorm ==> <psi|Hnorm|psi> = ", evals(1), evals(2), evals(3), evals(4)

!        worm_hv%mat = 2.0_double*worm_hv%mat/egrid
!        worm_hv%mat = worm_hv%mat - esum/egrid*worm_v%mat

        deallocate(evals)

100     call glean(thy(v))
        call glean(thy(hv))
        call glean(thy(h))

        if (error("Exit time_propagator_mod::h_norm_i")) continue

      end subroutine



!*******************************************************************************



!*********************************************************************************************************
!*********************************************************************************************************

      !** sop2_i
      !
      ! Propagate using the Suzuki-Trotter 2nd Order Split Operator method described by Sugino and 
      !  Miyamoto in their 1999 paper: PRB 59, p2579
      !
      !     
      !     
      !    v(dt) = S(dt,t) v(0)
      !             2           
      !
      !                              ASCENDING                              DESCENDING
      !                          __           ps                        __           ps
      ! S(x,t) = exp(-i*x*T/2) { || exp(-i*x*V(t)/2) } exp( -i*V(t) ) { || exp(-i*x*V(t)/2) } exp(-i*x*T/2)
      !  2                       p            p                 Hxc                  p
      !
      !           ps 
      !   Where  V    are the pseudopotentials (there are p pseudopotentials)
      !           p
      !
      subroutine sop2_i(tp,h,v)
        type(time_propagator_obj) :: tp
        type(hamiltonian_obj)     :: h
        type(multivector_obj)     :: v

        ! local vars
        type(multivector_rep), pointer :: worm_v, worm_new_v
        real(double)                   :: dt
        complex(double), pointer       :: tmp(:,:)

        complex(double), parameter     :: i = (0.0_double, 1.0_double)
        complex(double)                :: check_tmp, check_tmp_local

        !-------------------------------------------------------------------------------------------

        call my(tp)
        call my(v)
        call my(h)

        ! size of the time step
        dt = tp%o%dt

        worm_v         => wormhole(v)
        worm_new_v     => wormhole(tp%o%mvs(1))

        ! call the apply_split_operator routine which returns mvs(1) as the time evolved multivector
!        call apply_split_operator(v,tp%o%mvs(1),h,dt) ; if (error()) GOTO 100

        ! perform a swap operation so that the time evolved multivector is returned in v.
        tmp        => worm_v%mat
        worm_v%mat => worm_new_v%mat
        worm_new_v%mat => tmp

!write(*,*) 'exact <psi|U|psi> = ', exp(-i*dt*eigenvals(1,2))
check_tmp_local = dot_product(worm_new_v%mat(:,2),worm_v%mat(:,2))
call allreduce(MOD_SCOPE,mpi_sum,check_tmp_local,check_tmp)


!write(*,*) ' sopd <psi|U|psi> = ', check_tmp

!check_tmp = check_tmp - exp(-i*dt*eigenvals(1,2))
!write(*,*) '             diff = ', check_tmp, abs(check_tmp)

        nullify(worm_v, worm_new_v, tmp)

100     call glean(thy(tp))
        call glean(thy(v))
        call glean(thy(h))

      end subroutine


!*********************************************************************************************************
!*********************************************************************************************************

      !** sop4_i
      !
      ! Propagate using the Suzuki-Trotter 4th Order Split Operator method described by Sugino and 
      !  Miyamoto in their 1999 paper: PRB 59, p2579
      !
      !     
      !     
      !    v(dt) = S(p1*dt, t1) S(p1*dt, t1) S(p1*dt, t1) S(p1*dt, t1) S(p1*dt, t1) v(0)
      !             2            2            2            2            2
      !
      !
      !                              ASCENDING                              DESCENDING
      !                          __           ps                        __           ps
      ! S(x,t) = exp(-i*x*T/2) { || exp(-i*x*V(t)/2) } exp( -i*V(t) ) { || exp(-i*x*V(t)/2) } exp(-i*x*T/2)
      !  2                       p            p                 Hxc                  p
      !
      !           ps 
      !   Where  V    are the pseudopotentials (there are p pseudopotentials)
      !           p
      !
      subroutine sop4_i(tp,h,v)
        type(time_propagator_obj) :: tp
        type(hamiltonian_obj)     :: h
        type(multivector_obj)     :: v

        ! local vars
        type(multivector_rep), pointer :: worm_v, worm_new_v
        real(double)                   :: dt, p1,p2,p3
        complex(double), pointer       :: tmp(:,:)


        complex(double), parameter     :: i = (0.0_double, 1.0_double)
        complex(double)                :: check_tmp, check_tmp_local

        !-------------------------------------------------------------------------------------------

        call my(tp)
        call my(v)
        call my(h)


        dt = tp%o%dt

        worm_v         => wormhole(v)
        worm_new_v     => wormhole(tp%o%mvs(1))


        p1 = 1.d0/(4.d0 - 4.d0**(1.d0/3.d0))
        p2 = p1
        p3 = 1 - 4*p1
       
!        call apply_split_operator(v,tp%o%mvs(1),h,p1*dt) ; if (error()) GOTO 100
        tmp        => worm_v%mat
        worm_v%mat => worm_new_v%mat
        worm_new_v%mat => tmp

!        call apply_split_operator(v,tp%o%mvs(1),h,p2*dt) ; if (error()) GOTO 100
        tmp        => worm_v%mat
        worm_v%mat => worm_new_v%mat
        worm_new_v%mat => tmp

!        call apply_split_operator(v,tp%o%mvs(1),h,p3*dt) ; if (error()) GOTO 100
        tmp        => worm_v%mat
        worm_v%mat => worm_new_v%mat
        worm_new_v%mat => tmp

!        call apply_split_operator(v,tp%o%mvs(1),h,p2*dt) ; if (error()) GOTO 100
        tmp        => worm_v%mat
        worm_v%mat => worm_new_v%mat
        worm_new_v%mat => tmp

!        call apply_split_operator(v,tp%o%mvs(1),h,p1*dt) ; if (error()) GOTO 100
        tmp        => worm_v%mat
        worm_v%mat => worm_new_v%mat
        worm_new_v%mat => tmp

        
!write(*,*) 'exact <psi|U|psi> = ', exp(-i*dt*eigenvals(1,2))
!check_tmp_local = dot_product(worm_new_v%mat(:,2),worm_v%mat(:,2))
!call allreduce(mpi_sum,check_tmp_local,check_tmp)
!write(*,*) mpi_myproc(), 'sop4d <psi|U|psi> = ', check_tmp

!check_tmp = check_tmp - exp(-i*dt*eigenvals(1,2))
!write(*,*) '             diff = ', check_tmp, abs(check_tmp)

        nullify(worm_v, worm_new_v, tmp)

100     call glean(thy(tp))
        call glean(thy(v))
        call glean(thy(h))

      end subroutine


      !** sil_i
      !
      ! Propagate using the Short Iterative Lanczos Method
      !                                   ->
      !  1) Construct an Arnoldi basis, {w_j}m from the m'th order Krylov Subspace
      !      of the Hamiltonian operator.
      !   
      !     Let w_0 = v(0)
      !         b_0 = 0
      !     
      !     do j=1,m   
      !       |w_j> = H|w_(j-1)> - b_(j-1)*|w_(j-2)>
      !       a_(j-1) = <w_j|w_(j-1)>
      !       |w_j> = |w_j> - a_(j-1)*|w_(j-1)>
      !       b_(j) = sqrt(||w_j||^2)
      !       if (a_(j) < tol) exit
      !       |w_(j)> = |w_j>/b_(j)
      !     end do
      !
      !  2) The projection of the Hamiltonian on this reduced arnoldi basis, {w_j}m, 
      !      is a symmetric tridiagonal matrix. a_j and b_j are the alpha and beta 
      !      matrix elements.
      !            _                                    _
      !           |a_1 b_1                               |
      !           |b_1 a_2 b_2                           |
      !           |    b_2 a_3 b_3                       |
      !           |        b_3 a_4                       |
      !      Tm = |                 .                    | 
      !           |                    .                 | 
      !           |                       .              |
      !           |                                b_m-1 |
      !           |_                         b_m-1  a_m _| 
      !
      !
      !  3) Diagonalize the reduced Hamiltonian, Tm, in the reduced basis of length m.
      !   
      !        Tm|d_j> = e_j*|d_j>
      !     
      !     where |d_j> and e_j are the jth eigenvector and eigenvalue of Tm
      !
      !  4) Propagate the wavefunctions in this basis:       
      !         
      !        |d_j(t)> = exp(-i*e_j*dt)|d_j(0)>
      !        
      !  5) Project the propagated wavefunction onto the full basis:
      !                 
      !                 m
      !                 __
      !                 \  
      !        |v(t)> = /_  ( d_j(t) )*|w_j>
      !                 j=1
      !        
      !        
      subroutine sil_i(tp,h,v)
        type(time_propagator_obj) :: tp
        type(hamiltonian_obj)     :: h
        type(multivector_obj)     :: v

        ! local vars
        integer                 :: ib, nb, j, k,mysz
        integer                 :: num_cached_vecs ! number of cached multivectors in tp%o%mvs(:)
        integer                 :: num_basis_vecs  ! number of arnoldi basis vectors
        integer,parameter       :: maxiters = 25   ! maximum number of iterations for a given calculation.
        real(double)            :: maxerr
        real(double), pointer   :: errors(:)
        real(double), pointer   :: alphas(:,:)  ! Holds all the alpha matrix elements
        real(double), pointer   :: betas(:,:)   ! Holds all the beta matrix elemetns
        real(double), pointer   :: diag(:)      ! Holds the diagonal elements
        real(double), pointer   :: subdiag(:)   ! Holds the sub diagonal elements
        real(double), pointer   :: beta(:)      ! matrix elements of Tm, corresponds to b_j above
        real(double), pointer   :: invbeta(:)   ! 1/beta
        real(double), pointer   :: evecs(:,:) ! matrix that holds the evecs of the tridiagonal matrix
        real(double), pointer   :: evals(:)   ! vectors that holds the evals of the tridiagonal matrix
        real(double), pointer   :: one_r(:)     ! vector with each element set to 1.0
        complex(double),pointer :: one_c(:)     ! vector with each element set to (1.0,0.0) 
        complex(double),pointer :: alpha(:)     ! matrix elements of Tm, corresponds to c_j above
        complex(double),pointer :: new_av(:,:)  ! propagated wavefunction in the arnoldi basis
        complex(double),pointer :: new_dv(:)  ! propagated wavefunction in the diagonalized basis
        complex(double), parameter :: i = (0.0_double, 1.0_double)
        integer, pointer        :: basis_size(:) ! number of arnoldi basis vectors for each band
        logical                 :: isfinished
        !-------------------------------------------------------------------------------------------

        call my(tp)
        call my(v)
        call my(h)

        !** allocate the various vectors
        nb = x_n_bands(v)

        num_cached_vecs = size(tp%o%mvs) 
        allocate( alpha(nb), beta(nb), invbeta(nb),one_r(nb), one_c(nb) )
        allocate( errors(nb), basis_size(nb) )
!        num_basis_vecs = num_cached_vecs + 1
!        allocate( alphas(num_basis_vecs,nb), betas(num_basis_vecs,nb) )
        allocate( alphas(num_cached_vecs,nb), betas(num_cached_vecs,nb) )
        
        one_c = (1.0_double,0.0_double)
        one_r = 1.0_double
        basis_size = -1   ! Initialize the basis size to -1. Will be reset once the number of basis
                          !  functions has been found.
        
        !** Calculate the 1st arnoldi basis vector and the first set of betas and alphas

        !** Form H|v>
        call apply_hamiltonian(v,tp%o%mvs(1),h) ; if (error()) goto 100
        
        !** For the initialization step, wj = Hv, so don't bother setting wj
        !** Calculate the alpha matrix element for each evec - alpha_j = <v_j|H|v_j>
        call multiply(v,tp%o%mvs(1),alpha)
       
        !** set wj = wj - alpha*v
        call combine(one_c,tp%o%mvs(1),-alpha,v)

        !** Calculate the beta matrix element.
        call multiply(tp%o%mvs(1),beta)
        beta = sqrt(beta)

        where (beta < tp%o%tol) basis_size = 1

        isfinished = .true.
        do ib = 1,nb
           if (basis_size(ib) < 0) isfinished = .false.
        end do
        
        if (.not. isfinished) then
           where (tp%o%tol < beta) 
              invbeta = 1.0d0/beta
           elsewhere
              invbeta = 1.0d0
           end where
           ! set the first arnoldi basis vector: w1 = w1/beta1
           call portion(invbeta,tp%o%mvs(1))
        end if

        !** Load the first elements of the alphas and betas
        alphas(1,:) = real(alpha)
        !alphas(1,:) = abs(alpha)
        betas(1,:)  = beta

        j = 1
        !** Calculate the error associated with each band.
        errors(:) = tp%o%dt*beta
        maxerr = max(errors)
        
        if (maxerr < tp%o%tol) then
           isfinished = .true.
        end if

        !************************************************************************
        !** Start the loop that fills out the arnoldi basis.  *******************
        do while ( .not.isfinished )  

           j = j + 1

           if (num_cached_vecs < j) then
              call reallocate_memory_i(j,tp%o,alphas,betas)
              num_cached_vecs = j
           end if
           
           ! Apply the Hamiltonian 
           call apply_hamiltonian(tp%o%mvs(j-1),tp%o%mvs(j),h) ; if (error()) GOTO 100

           ! Set |w_j> = H|w_(j-1)> - beta_j*|w_(j-2)>
           if (j < 3) then
              call combine(one_r,tp%o%mvs(j),-beta,v)
           else
              call combine(one_r,tp%o%mvs(j),-beta,tp%o%mvs(j-2))
           end if

           ! Calculate the new alpha: a_(j-1) = <w_(j-1)|H|w_(j-1)> = <w_j|w_(j-1)> 
           call multiply(tp%o%mvs(j),tp%o%mvs(j-1),alpha)

           ! set |w_j> = |w_j> - alpha_(j-1)*|w_(j-1)>
           call combine(one_c,tp%o%mvs(j),-alpha,tp%o%mvs(j-1))

           !** Calculate the beta matrix element.  
           !     beta = sqrt( || H|w_(j-1)> - alpha_(j-1)|w_(j-1)> - beta_(j-1)|w_(j-2)> ||)
           call multiply(tp%o%mvs(j),beta)
           beta = sqrt(beta)

           !** Normalize the new vector.
           do ib=1,nb 
              if ( tp%o%tol < beta(ib) ) then
                 invbeta(ib) = 1.0d0/beta(ib)
              else
                 invbeta(ib) = 0.0_double
                 if (basis_size(ib) < 0) basis_size(ib) = j
              end if

              !** zero out the wavefunction so it doesn't contribute to the solution
              if ( 0 < basis_size(ib) ) invbeta(ib) = 0.0_double
           end do
           call portion(invbeta,tp%o%mvs(j))

           !** Load the elements of the diagonal and subdiagonal vectors
           alphas(j,:)  = real(alpha)
           !alphas(j,:)  = abs(alpha)
           betas(j,:) = beta

           !** Check for Stop conditions:

           !** 1) If no vectors are left then the the calculation is finished
           isfinished = .true.
           do ib = 1,nb
              if (basis_size(ib) < 0) isfinished = .false.
           end do

           !** 2) If the largest error is small enough, then the calculation is finished
           errors(:) = (tp%o%dt)**j/factorial(j)*product(betas(1:j,1))
           maxerr = max(errors)
           if (maxerr < tp%o%tol) isfinished = .true.

           !** 3) Total number of iterations is larger than some fixed number
           if (maxiters < j) isfinished = .true.



        end do !** End big loop *************************************************
        !************************************************************************

        tp%o%num_hamiltonian_ops = j
        num_basis_vecs = j

        allocate( new_av(num_basis_vecs,nb) )

        !** Iterate through all the bands
        do ib=1,nb

           if (basis_size(ib) < 0) basis_size(ib) = num_basis_vecs
           mysz = basis_size(ib)

           !** Allocate the various structures that will be needed
           allocate(evecs(mysz,mysz),evals(mysz))
           allocate(diag(mysz),subdiag(mysz))
           allocate( new_dv(mysz) )
           
           !** Load the diagonal and subdiagonal elements
           diag = alphas(1:mysz,ib)
           subdiag = betas(1:mysz,ib)
           
           !** Diagonalize the Arnoldi basis
           call diagonalize_i(diag,subdiag,evals,evecs)
        
           !** Propagate the wavefunctions in the diagonalized reduced basis
           new_dv(:) = evecs(1,:)*exp(-i*tp%o%dt*evals(:))
           !new_dv(k) = evecs(k,1)*exp(-i*tp%o%dt*evals(k))

           !** Transform the solution back to the arnoldi basis
           new_av(:,ib) = 0.0_double
           do k=1,mysz
              new_av(k,ib) = sum(evecs(k,1:mysz)*new_dv(1:mysz))
           end do

           deallocate(evecs,evals,diag,subdiag,new_dv)

        end do

        !** Project the Solution back out to the full planewave basis:
        ! The first term is just the original wavefunction times the first coefficient
        call portion(new_av(1,:),v)

        ! Iterate through the remaining basis vectors
        do j=2, num_basis_vecs
           ! set |v> = |v> + new_av_(j-1)*|w_(j-1)>
           call combine(one_c,v,new_av(j,:),tp%o%mvs(j-1))
        end do

        !** Normalize v  !!! This probably shouldn't be necessary !!!
!        call multiply(v,beta)
!        invbeta = 1.0/sqrt(beta)
!        call portion(invbeta,v)

100     if (associated(errors)) deallocate(errors)
        if (associated(alphas)) deallocate(alphas)
        if (associated(betas)) deallocate(betas)
        if (associated(diag)) deallocate(diag)
        if (associated(subdiag)) deallocate(subdiag)
        if (associated(beta)) deallocate(beta)
        if (associated(invbeta)) deallocate(invbeta)
        if (associated(evecs)) deallocate(evecs)
        if (associated(evals)) deallocate(evals)
        if (associated(one_r)) deallocate(one_r)
        if (associated(one_c)) deallocate(one_c)
        if (associated(alpha)) deallocate(alpha)
        if (associated(new_av)) deallocate(new_av)
        if (associated(new_dv)) deallocate(new_dv)
        if (associated(basis_size)) deallocate(basis_size)

        call glean(thy(tp))
        call glean(thy(v))
        call glean(thy(h))

        if (error("Exit time_propagator_mod::sil_i")) continue

      end subroutine


      subroutine reallocate_memory_i(new_num_vecs,tpr,alphas,betas)
        integer                   :: new_num_vecs
        type(time_propagator_rep) :: tpr
        real(double), pointer     :: alphas(:,:)
        real(double), pointer     :: betas(:,:)
        ! local vars -----------------------------------------------------------
        integer                        :: nb, iv, nv
        integer                        :: cur_num_vecs
        real(double), pointer          :: tmp_mat(:,:)  ! Holds all the alpha matrix elements
        type(multivector_obj), pointer :: tmp_mvs(:)
        !-----------------------------------------------------------------------

        !** Get the number of basis vectors and check for errors
        cur_num_vecs = size(tpr%mvs)
        if (error(cur_num_vecs /= size(alphas,1),'size alphas /= size mvs')) goto 100
        if (error(cur_num_vecs /= size(betas,1),'size betas /= size mvs')) goto 100
        
        !** Get the number of bands and check for errors
        nb = size(alphas,2)
        if (error(nb /= size(betas,2),'size(alphas,2) /= size(betas,2)')) goto 100

        !** Allocate temporary array to move alphas and betas info
        allocate( tmp_mat(new_num_vecs,nb) )
        
        !** Reload the alphas array
        nv = min(new_num_vecs,cur_num_vecs)
        
        do iv=1,nv
           tmp_mat(iv,:) = alphas(iv,:)
        end do

        deallocate (alphas)
        allocate(alphas(new_num_vecs,nb))
        alphas = 0.0_double

        do iv=1,nv
           alphas(iv,:) = tmp_mat(iv,:)
        end do
        
        !** Reload the betas array
        do iv=1,nv
           tmp_mat(iv,:) = betas(iv,:)
        end do

        deallocate (betas)
        allocate(betas(new_num_vecs,nb))
        betas = 0.0_double

        do iv=1,nv
           betas(iv,:) = tmp_mat(iv,:)
        end do
        
        !** Allocate the temporary array for the multivector
        allocate(tmp_mvs(nv))

        !** Load the temporary multivector array
        do iv=1,nv
           call my(tpr%mvs(iv),tmp_mvs(iv))
        end do

        do iv=1,cur_num_vecs
           call glean(thy(tpr%mvs(iv)))
        end do
        
        deallocate(tpr%mvs)
        allocate(tpr%mvs(new_num_vecs))

        do iv=1,nv
           call my(tmp_mvs(iv),tpr%mvs(iv))
        end do

        do iv=nv+1,new_num_vecs
           call my(tmp_mvs(1),tpr%mvs(iv))
        end do

        tpr%num_mvs = new_num_vecs


100     deallocate(tmp_mat)
        do iv=1,size(tmp_mvs)
           call glean(thy(tmp_mvs(iv)))
        end do
        deallocate(tmp_mvs)
        if (error('tprop::reallocate_memory_i - Exiting')) continue

      end subroutine


      function max(v) result(maxval)
        real(double), pointer :: v(:)
        real(double)          :: maxval

        real(double) :: tmp
        integer      :: i

        if (0 < size(v)) then
           
           tmp = v(1)
           do i=2,size(v)
              if (tmp < v(i)) tmp = v(i)
           end do
        else
           tmp = -1.0_double
        end if
        
        maxval = tmp

      end function


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

      subroutine my_tp(tp)
!doc$ subroutine my<type(time_propagator)>(tp)
        type(time_propagator_obj) :: tp

!cod$
        tp%ref = tp%ref + 1
        tp%o%ref = tp%o%ref + 1

      end subroutine

!*******************************************************************************

      subroutine my_new_tp(tpi,tp)
!doc$ subroutine my<type(time_propagator)>(tpi,tp)
        type(time_propagator_obj) :: tpi, tp

!cod$
        tp%ref = 1
        tp%o => tpi%o
        tp%o%ref = tp%o%ref + 1
      end subroutine

!*******************************************************************************

      function thy_tp(tp)
!doc$  function thy<type(time_propagator)>(tp)
        type(time_propagator_obj) :: tp, thy_tp

!cod$
        tp%ref = tp%ref - 1
        tp%o%ref = tp%o%ref - 1
        thy_tp%ref = tp%ref
        thy_tp%o => tp%o

      end function
      
!*********************************************************************************************************

      subroutine glean_tp(tp)
!doc$ subroutine glean<type(time_propagator)>(tp)
        type(time_propagator_obj) :: tp

!cod$
        integer :: i

        if (tp%o%ref < 1) then
           if (associated(tp%o%mvs)) then
              do i=1, tp%o%num_mvs
                 if (tp%o%is_mvs_allocated) call glean(thy(tp%o%mvs(i)))
              end do
              deallocate(tp%o%mvs)
           end if
           if (tp%o%is_tmpv_init) call glean(thy(tp%o%tmp_v))
           deallocate( tp%o )
        end if

      end subroutine

!*********************************************************************************************************

      subroutine bequeath_tp(tp)
!doc$ subroutine bequeath<type(time_propagator)>(tp)
        type(time_propagator_obj) :: tp

!cod$
        continue
      end subroutine
    
!*********************************************************************************************************

      subroutine assign_tp(tp,tp2) 
!doc$ subroutine assign<type(time_propagator)>(tp,tp2) 
        type(time_propagator_obj), intent(inout) :: tp
        type(time_propagator_obj), intent(in) :: tp2

        type(time_propagator_obj) :: tpt


        call my(tp2)
        tpt%o => tp%o
        tp%o%ref = tp%o%ref - tp%ref
        tp%o => tp2%o
        tp%o%ref = tp%o%ref + tp%ref
        call glean(tpt)
        call glean(thy(tp2))

      end subroutine

!*********************************************************************************************************

!      function is_finished_tp(tp) result(isfinished)
!doc$ function is_finished(tp) result(isfinished)
!        type(time_propagator_obj) :: tp
!        logical :: isfinished
!       effects: Returns whether this time_propagation object has finished iterating.  

!cod$
!        call my(tp)
!        isfinished = .true.
!        if (tp%o%elapsed_time < tp%o%runtime) isfinished = .false.
!        call glean(thy(tp))
        
!      end function is_finished_tp


!*********************************************************************************************************

      function get_num_hamiltonian_ops(tp) result(num)
!doc$ function   x_num_hamiltonian_ops(tp) result(num)
        type(time_propagator_obj) :: tp
        integer :: num
!       effects: Returns the number of H|Psi> operations per current time step
!         
!cod$
        call my(tp)
        num = tp%o%num_hamiltonian_ops
        call glean(thy(tp))

      end function 


!*********************************************************************************************************

!      function get_istep(tp) result(istep)
!doc$ function   x_istep(tp) result(time)
!        type(time_propagator_obj) :: tp
!        real(double)              :: istep
!       effects: Returns the current step number
!         
!cod$        

!        call my(tp)
!        istep = tp%o%istep
!        call glean(thy(tp))

!      end function 

!*********************************************************************************************************

      function need_current_density_tp(tp) result(need_currdens)
!doc$ function need_current_density(tp) result(need_currdens)
        type(time_propagator_obj) :: tp
        logical                   :: need_currdens
!       effects: Returns a logical that is true iff this calculation
!                requires the full current density to be calculated.
!cod$        

        call my(tp)
        need_currdens = tp%o%need_current_density
        call glean(thy(tp))

      end function 

!*********************************************************************************************************

      subroutine write_restart_tp(tp,nrestf)
!doc$ subroutine write_restart(tp,nrestf)
        type(time_propagator_obj)      :: tp
        type(tagio_obj), intent(inout) :: nrestf
!       modifies: nrestf
!       effects: Writes tp restart information to nrestf.
!       errors: Passes errors.

!cod$
        integer :: imvs
        integer(long) :: s4, dsize, ios, ndata 

        call my(tp)
        call my(nrestf)

        if (i_access(nrestf)) then
          call startblock(nrestf,"TIME_PROPAGATOR")
          call writetag(nrestf,"PARAMETERS")
          s4 = tp%o%num_mvs ; dsize = sizeof_long ; ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)           
        end if

        if (i_access(nrestf)) call writetag(nrestf,"MVS_HISTORY")
        do imvs = 1,tp%o%num_mvs
          call write_restart(tp%o%mvs(imvs),nrestf) ; if (error()) exit
        end do

        if (i_access(nrestf)) call endblock(nrestf)


        call glean(thy(tp))
        call glean(thy(nrestf))

        if (error("Exit time_propagator_mod::write_restart_tp")) continue

      end subroutine


      !********************************************************************************

      subroutine own_i(tp)
        type(time_propagator_obj), intent(inout) :: tp
        type(time_propagator_obj) :: tpt
        integer :: i
        if (tp%ref < tp%o%ref) then
          allocate( tpt%o )
          tpt%o%ref = 0
          tpt%o%g = tp%o%g

          ! params read in from argvf
          tpt%o%method = tp%o%method
          tpt%o%dt = tp%o%dt
          tpt%o%tol = tp%o%tol
          tpt%o%need_current_density = tp%o%need_current_density

          ! State variables
          tpt%o%num_hamiltonian_ops = tp%o%num_hamiltonian_ops
          tpt%o%num_mvs = tp%o%num_mvs

          tpt%o%is_tmpv_init = tp%o%is_tmpv_init
          call my(tp%o%tmp_v,tpt%o%tmp_v)

          tpt%o%is_mvs_allocated = tp%o%is_mvs_allocated
          allocate(tpt%o%mvs(size(tp%o%mvs)))
          if (tp%o%is_mvs_allocated) then
            allocate(tpt%o%mvs(size(tp%o%mvs)))
            do i = 1,size(tp%o%mvs)
              call my(tp%o%mvs(i),tpt%o%mvs(i))
            end do
          end if

          tp%o%ref = tp%o%ref - tp%ref
          tp%o => tpt%o
          tp%o%ref = tp%o%ref + tp%ref
        end if
      end subroutine

      !********************************************************************************

      subroutine diagonalize_i(diag,subdiag,evals,evecs)
        real(double), intent(inout) :: diag(:)
        real(double), intent(inout) :: subdiag(:)
        real(double), intent(inout) :: evals(:)
        real(double), intent(inout) :: evecs(:,:)

        character(1) :: jobz
        integer      :: n
        integer      :: ierr
        integer      :: ldz
        real(double), allocatable :: work(:)

        Call Start_Timer("tprop: diagonalize_i")

        n = size(diag)
        if (error(size(subdiag) /= n, "ERROR: dimensions of diag and subdiag don't agree")) goto 200
        if (error(size(evals) /= n, "ERROR: dimensions of diag and evals don't agree")) goto 200
        if (error(size(evecs,1) /= n, "ERROR: dimensions of diag and evecs columns don't agree")) goto 200
        if (error(size(evecs,1) /= size(evecs,2), "ERROR: evecs matrix is not square")) goto 200

        if (n < 1) then
           ldz = 1
        else
           ldz = n
        end if

        if (n < 2) then
           allocate(work(1))
        else
           allocate(work(2*n-2))
        end if

        jobz = 'v'
        call dstev(jobz,n,diag,subdiag,evecs,ldz,work,ierr)
!        write(iobuf,*) "ERROR: dstev ierr = ",ierr
!        if (error(ierr /= 0,iobuf)) goto 100
        if (error(ierr /= 0,"ERROR: dstev failed")) goto 100

        evals = diag
        
100     deallocate( work )

        Call Stop_Timer("tprop: diagonalize_i")

200     if (error("Exit time_propagator_mod::diagonalize_i")) continue


      end subroutine

!*********************************************************************************************************

      subroutine write_params_out_i(tpr)
        type(time_propagator_rep) :: tpr
        
        write(*,*) '*******************************************'
        write(*,*) '*********  Time Prop Params  **************'
        write(*,*) 'method = ', tpr%method
        write(*,*) 'dt = ', tpr%dt
        write(*,*) 'tol = ', tpr%tol
        write(*,*) 'need_current_density = ', tpr%need_current_density
        write(*,*) 'num_hamiltonian_ops = ', tpr%num_hamiltonian_ops
        write(*,*) 'num_mvs = ', tpr%num_mvs
        write(*,*) 'is_mvs_alloc = ', tpr%is_mvs_allocated
        write(*,*) 'is_tmpv_init = ', tpr%is_tmpv_init
        write(*,*) '*******************************************'

      end subroutine write_params_out_i

!*********************************************************************************************************

    end module

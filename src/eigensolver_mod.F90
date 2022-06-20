!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module eigensolver_mod
!doc$ module eigensolver_mod

!     One datatype is available here: type(eigensolver_obj).

!     eigensolver_mod encapsulates eigensolver methods and associated parameters.

      use kind_mod
      use arg_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use ghost_mod
      use math_mod
      use multivector_mod
      use operators_mod
      use timing_mod

!cod$
      implicit none
      private

      ! solver_method
      integer, parameter :: CG  = 1   ! Conjugate Gradient solver
      integer, parameter :: GCG = 2   ! Grassmann Conjugate Gradient solver
      integer, parameter :: BD  = 3   ! Block-Davidson solver

      ! diagonalization_method
      integer, parameter :: LAPACK     = 1   ! LAPACK diagonalization
      integer, parameter :: SCALAPACK  = 2   ! ScaLAPACK diagonalization

      ! profile types
      integer, parameter :: FIXED   = 1
      integer, parameter :: STEPPED = 2

      ! rf modes
      integer, parameter :: NOT_ENABLED = 0
      integer, parameter :: ENABLED  = 1

      ! rf mode status
      integer, parameter :: DONE    = 0
      integer, parameter :: PENDING = 1

      integer :: blocksize                  ! blocksize for ScaLAPACK diagonalization method

      type :: eigensolver_rep
        integer :: ref
        type(ghost) :: g
        type(file_obj) :: f                            ! solver_report
        logical :: report                              ! indicates whether or not to print a report
        integer :: solver_method                       ! solver method
        integer :: diagonalization_method              ! diagonalization method
        integer :: dir_profile                         ! directions profile type
        integer :: tol_profile                         ! residual tolerance profile type
        integer :: max_steps                           ! maximum number of steps
        integer :: step                                ! step
        integer :: rf_mode                             ! restart f mode
        integer :: rfm_status                          ! restart f mode status
        integer :: rfm_dir                             ! restart f mode directions
        real(double) :: rfm_tol                        ! restart f mode tolerance
        integer, dimension(:), pointer :: dir          ! step-dependent number of solver directions
        real(double), dimension(:), pointer :: tol     ! step-dependent residual tolerance for solver
      end type

      type, public :: eigensolver_obj
        private
        integer :: ref
        type(eigensolver_rep), pointer :: o
      end type

!doc$
      public :: eigensolver
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: eigensolve
      public :: diary

!cod$
      interface eigensolver
        module procedure constructor_es
      end interface
      interface update
        module procedure update_es
      end interface
      interface my
        module procedure my_es, my_new_es
      end interface
      interface thy
        module procedure thy_es
      end interface
      interface glean
        module procedure glean_es
      end interface
      interface bequeath
        module procedure bequeath_es
      end interface
      interface assignment(=)
        module procedure assign_es
      end interface
      interface eigensolve
        module procedure eigensolve_es
      end interface
      interface diary
        module procedure diary_es
      end interface

      contains

! public routines

      function constructor_es() result(es)
!doc$ function eigensolver() result(es)
        type(eigensolver_obj) :: es
!       effects: Constructs a new es.
!       errors: Input keywords or parameters.

!cod$
        logical :: found
        character(line_len) :: mode, tag
        integer :: d1, d2, d3, is, ios, nl, s12, s23
        real(double) :: t1, t2, t3

        es%ref = 0
        allocate( es%o )
        es%o%ref = 0
        es%o%g = x_ghost()

        ! determine the reporting status
        call arglc("solver_report",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("off")
          es%o%report = .false.
        case ("on")
          es%o%report = .true.
        case default
          if (error(.true.,"WARNING: solver_report not recognized")) continue
        end select

        ! determine the rf_mode
        es%o%rf_mode = NOT_ENABLED
        call arglc("config_type",tag,found)
        if (.not.found) tag = "self-consistent"
        select case (trim(tag))
        case ("self-consistent","sc")
          call arglc("restart",mode,found)
          if (.not.found) mode = "off"
          select case (trim(mode))
          case ("f","ef")
            call arglc("solver_rf_mode",mode,found)
            if (.not.found) mode = "on"
            select case (trim(mode))
            case ("on")
              es%o%rf_mode = ENABLED
              es%o%rfm_status = PENDING
            end select
          end select
        end select

        ! get the maximum number of solver steps
        call arg("config_steps",es%o%max_steps,found)
        if (.not.found) es%o%max_steps = 40 ! Note: this default must be consistent with the one in config_sc_mod.f90

        ! get solver method
        call arglc("solver_method",tag,found)
        if (.not.found) tag = "bd"
        select case (trim(tag))
        case ("conjugate-gradients","cg")
          es%o%solver_method = CG
        case ("grassman-conjugate-gradients","gcg")
          es%o%solver_method = GCG
        case ("blocked-davidson","bd")
          es%o%solver_method = BD
        case default
          if (error(.true.,"ERROR: unrecognized solver_method")) goto 100
        end select

        ! get the solver directions profile type
        call arglc("solver_dir_profile",tag,found)
        if (.not.found) tag = "fixed"
        select case (trim(tag))
        case ("fixed")
          es%o%dir_profile = FIXED
        case ("stepped")
          es%o%dir_profile = STEPPED
        case default
          if (error(.true.,"ERROR: solver_dir_profile not recognized")) goto 100
        end select

        ! set solver directions profile values
        allocate( es%o%dir(es%o%max_steps) )
        call arg("solver_directions",d1,found)
        if (.not.found) then
          call arglc("config_type",tag,found)
          if (.not.found) tag = "self-consistent"
          select case (trim(tag))
          case ("self-consistent","sc")
            select case (es%o%solver_method)
            case (CG)
              d1 = 3
            case (GCG)
              d1 = 10
            case (BD)
              d1 = 3
            end select
          case ("fixed-hamiltonian","fh")
            select case (es%o%solver_method)
            case (CG)
              d1 = 5
            case (GCG)
              d1 = 12
            case (BD)
              d1 = 5
            end select
          end select
        end if
        if (error(d1 < 1,"ERROR: solver_directions < 1")) goto 100
        select case (es%o%dir_profile)
        case (FIXED)
          es%o%dir = d1
        case (STEPPED)
          call arg("solver_dir_levels",nl,found)
            if (error(.not.found,"ERROR: solver_dir_levels not found")) goto 100
            if (error(nl < 2,"ERROR: solver_dir_levels < 2")) goto 100
            if (error(nl > 3,"ERROR: solver_dir_levels > 3")) goto 100
          call arg("solver_dir_step_1-2",s12,found)
            if (error(.not.found,"ERROR: solver_dir_step_1-2 not found")) goto 100
            if (error(s12 < 2,"ERROR: solver_dir_step_1-2 < 2")) goto 100
            if (error(s12 > (es%o%max_steps - 1),"ERROR: solver_dir_step_1-2 too large")) goto 100
          call arg("solver_dir_2",d2,found)
            if (error(.not.found,"ERROR: solver_dir_2 not found")) goto 100
            if (error(d2 < 1,"ERROR: solver_dir_2 < 1")) goto 100
            if (error(d2 == d1,"ERROR: solver_dir_2 = solver_directions")) goto 100
          do is = 1,(s12-1)
            es%o%dir(is) = d1
          end do
          if (nl == 2) then
            do is = s12,es%o%max_steps
              es%o%dir(is) = d2
            end do
          else
            call arg("solver_dir_step_2-3",s23,found)
              if (error(.not.found,"ERROR: solver_dir_step_2-3 not found")) goto 100
              if (error(s23 <= s12,"ERROR: solver_dir_step_2-3 <= solver_dir_step_1-2")) goto 100
              if (error(s23 >= (es%o%max_steps - 1),"ERROR: solver_dir_step_2-3 too large")) goto 100
            call arg("solver_dir_3",d3,found)
              if (error(.not.found,"ERROR: solver_dir_3 not found")) goto 100
              if (error(d3 < 1,"ERROR: solver_dir_3 < 1")) goto 100
              if (error(d3 == d2,"ERROR: solver_dir_3 = solver_dir_2")) goto 100
            do is = s12,(s23-1)
              es%o%dir(is) = d2
            end do
            do is = s23,es%o%max_steps
              es%o%dir(is) = d3
            end do
          end if
        end select
        
        ! get the restart f mode directions
        select case (es%o%rf_mode)
        case (ENABLED)
          call arg("solver_rfm_directions",es%o%rfm_dir,found)
          if (.not.found) then
            select case (es%o%solver_method)
            case (CG)
              es%o%rfm_dir = 15
            case (GCG)
              es%o%rfm_dir = 15
            case (BD)
              es%o%rfm_dir = 15
            end select
          end if
          if (error(es%o%rfm_dir <= 1,"ERROR: solver_rfm_directions <= 1")) goto 100
        end select

        ! get the solver tolerance profile type
        call arglc("solver_tol_profile",tag,found)
        if (.not.found) tag = "fixed"
        select case (trim(tag))
        case ("fixed")
          es%o%tol_profile = FIXED
        case ("stepped")
          es%o%tol_profile = STEPPED
        case default
          if (error(.true.,"ERROR: solver_tol_profile not recognized")) goto 100
        end select

        ! set solver tolerance profile values
        allocate( es%o%tol(es%o%max_steps) )
        call arg("solver_tolerance",t1,found)
        if (.not.found) then
          call arglc("config_type",tag,found)
          if (.not.found) tag = "self-consistent"
          select case (trim(tag))
          case ("self-consistent","sc")
            select case (es%o%solver_method)
            case (CG)
              t1 = 1.0e-4_double
            case (GCG)
              t1 = 1.0e-5_double
            case (BD)
              t1 = 1.0e-4_double
            end select
          case ("fixed-hamiltonian","fh")
            select case (es%o%solver_method)
            case (CG)
              t1 = 5.0e-6_double
            case (GCG)
              t1 = 1.0e-8_double
            case (BD)
              t1 = 5.0e-6_double
            end select
          end select
        end if
        if (error(t1 < 0.0_double,"ERROR: solver_tolerance < 0")) goto 100
        select case (es%o%tol_profile)
        case (FIXED)
          es%o%tol = t1
        case (STEPPED)
          call arg("solver_tol_levels",nl,found)
            if (error(.not.found,"ERROR: solver_tol_levels not found")) goto 100
            if (error(nl > 3,"ERROR: solver_tol_levels > 3")) goto 100
          call arg("solver_tol_step_1-2",s12,found)
            if (error(.not.found,"ERROR: solver_tol_step_1-2 not found")) goto 100
            if (error(s12 < 2,"ERROR: solver_tol_step_1-2 < 2")) goto 100
            if (error(s12 > (es%o%max_steps - 1),"ERROR: solver_tol_step_1-2 out of bounds")) goto 100
          call arg("solver_tol_2",t2,found)
            if (error(.not.found,"ERROR: solver_tol_2 not found")) goto 100
            if (error(t2 < 0.0_double,"ERROR: solver_tol_2 < 0")) goto 100
            if (error(t2 == t1,"ERROR: solver_tol_2 = solver_tolerance")) goto 100
          do is = 1,(s12-1)
            es%o%tol(is) = t1
          end do
          if (nl == 2) then
            do is = s12,es%o%max_steps
              es%o%tol(is) = t2
            end do
          else
            call arg("solver_tol_step_2-3",s23,found)
              if (error(.not.found,"ERROR: solver_tol_step_2-3 not found")) goto 100
              if (error(s23 <= s12,"ERROR: solver_tol_step_2-3 <= solver_tol_step_1-2")) goto 100
              if (error(s23 >= (es%o%max_steps - 1),"ERROR: solver_tol_step_2-3 out of bounds")) goto 100
            call arg("solver_tol_3",t3,found)
              if (error(.not.found,"ERROR: solver_tol_3 not found")) goto 100
              if (error(t3 < 0.0_double,"ERROR: solver_tol_3 < 0")) goto 100
              if (error(t3 == t2,"ERROR: solver_tol_3 = solver_tol_2")) goto 100
            do is = s12,(s23-1)
              es%o%tol(is) = t2
            end do
            do is = s23,es%o%max_steps
              es%o%tol(is) = t3
            end do
          end if
        end select

        ! get the restart f mode tolerance
        select case (es%o%rf_mode)
        case (ENABLED)
          call arg("solver_rfm_tolerance",es%o%rfm_tol,found)
          if (.not.found) then
            select case (es%o%solver_method)
            case (CG)
              es%o%rfm_tol = 1.0e-5_double
            case (GCG)
              es%o%rfm_tol = 1.0e-8_double
            case (BD)
              es%o%rfm_tol = 1.0e-5_double
            end select
          end if
          if (error(es%o%rfm_tol < 0.0_double,"ERROR: solver_tolerance < 0")) goto 100
        end select

        ! get the diagonalization method
        call arglc("diagonalization_method",tag,found)
        if (.not.found) tag = "lapack"
        select case (trim(tag))
        case ("lapack")
          es%o%diagonalization_method = LAPACK
        case ("scalapack")
          es%o%diagonalization_method = SCALAPACK
        case default
          if (error(.true.,"ERROR: unrecognized diagonalization_method")) goto 100
        end select

        ! get scalapack specific values
        select case (es%o%diagonalization_method)
        case (SCALAPACK)
          call arg("blocksize",blocksize,found)
          if (.not.found) blocksize = 30
        case default
          blocksize = 0
        end select

        ! initialize the step counter
        es%o%step = 0

        ! open the report file
        if (es%o%report) then
          call my(file(trim(solver_report_path)),es%o%f)
          if (i_access(es%o%f)) open(x_unit(es%o%f),file=x_name(es%o%f),status='unknown',iostat=ios)
        end if

100     if (error("Exit eigensolver_mod::constructor_es")) continue

      end function

      subroutine update_es(es)
!doc$ subroutine update(es)
        type(eigensolver_obj) :: es
!       modifies: es
!       effects: Resets es%o%step to 0.

!cod$
        call my(es)

        call own_i(es)
        es%o%g = x_ghost()

        es%o%step = 0

        call glean(thy(es))

        if (error("Exit eigensolver_mod::update_es")) continue

      end subroutine

      subroutine my_es(es)
!doc$ subroutine my(es)
        type(eigensolver_obj) :: es

!cod$
        es%ref = es%ref + 1
        es%o%ref = es%o%ref + 1
      end subroutine

      subroutine my_new_es(esi,es)
!doc$ subroutine my(esi,es)
        type(eigensolver_obj) :: esi, es

!cod$
        es%ref = 1
        es%o => esi%o
        es%o%ref = es%o%ref + 1
      end subroutine

      function thy_es(es) result(eso)
!doc$  function thy(es) result(eso)
        type(eigensolver_obj) :: es, eso

!cod$
        es%ref = es%ref - 1
        es%o%ref = es%o%ref - 1
        eso%ref = es%ref
        eso%o => es%o
      end function
      
      subroutine glean_es(es)
!doc$ subroutine glean(es)
        type(eigensolver_obj) :: es

!cod$
        if (es%o%ref < 1) then
          if (es%o%report) then
            if (i_access(es%o%f)) close(x_unit(es%o%f))
            call glean(thy(es%o%f))
          end if
          deallocate( es%o%dir )
          deallocate( es%o%tol )
          deallocate( es%o )
        end if
      end subroutine

      subroutine bequeath_es(es)
!doc$ subroutine bequeath(es)
        type(eigensolver_obj) :: es

!cod$
        continue
      end subroutine
    
      subroutine assign_es(es,es2) 
!doc$ subroutine assign(es,es2) 
        type(eigensolver_obj), intent(inout) :: es
        type(eigensolver_obj), intent(in) :: es2

!cod$
        type(eigensolver_obj) :: est
        call my(es2)
        est%o => es%o
        es%o%ref = es%o%ref - es%ref
        es%o => es2%o
        es%o%ref = es%o%ref + es%ref
        call glean(est)
        call glean(thy(es2))
      end subroutine

      subroutine eigensolve_es(es,h,v,evals,rnorm)
!doc$ subroutine eigensolve(es,h,v,evals,rnorm)
        type(eigensolver_obj) :: es
        type(hamiltonian_obj) :: h
        type(multivector_obj) :: v
        real(double), dimension(:), intent(inout) :: evals
        real(double), intent(out) :: rnorm
!       requires: CG: Consistent dimensions. v contains orthonormal vectors.
!                 GCG: Consistent dimensions.
!       modifies: v, evals, rnorm
!       effects: Solves the eigenproblem given by h with initial guess v.
!       errors: Passes errors.

!cod$
        call my(es)
        call my(h)
        call my(v)

        es%o%step = es%o%step + 1
        if (error(es%o%step > size(es%o%dir),"ERROR: step is out of bounds")) goto 100

        select case (es%o%solver_method)
        case (CG)
          if (overlap_is_identity(h)) then
            call solve_cg_i(es,h,v,evals,rnorm) ; if (error()) goto 100
          else
            if (error(.true.,"ERROR: CG solver is not supported for non-identity overlap")) goto 100
          end if
        case (GCG)
          if (overlap_is_identity(h)) then
            call solve_gcg_i(es,h,v,evals,rnorm) ; if (error()) goto 100
          else
            if (error(.true.,"ERROR: GCG solver is not supported for non-identity overlap")) goto 100
          end if
        case (BD)
          if (overlap_is_identity(h)) then  
            call solve_bd_i(es,h,v,evals,rnorm) ; if (error()) goto 100
          else
            call solve_bdgen_i(es,h,v,evals,rnorm) ; if(error()) goto 100
          end if
        end select

        call glean(thy(es))
        call glean(thy(h))
        call glean(thy(v))

100     if (error("Exit eigensolver_mod::eigensolve_es")) continue

      end subroutine

      subroutine diary_es(es)
!doc$ subroutine diary(es)
        type(eigensolver_obj) :: es
!       effects: Writes es information to the diary.

!cod$
        integer :: is
        call my(es)
        if (i_access( diaryfile() )) then
          select case (es%o%solver_method)
          case (CG)
            write(x_unit(diaryfile()),'(/,t4,"Conjugate gradient solver:")')
          case (GCG)
            write(x_unit(diaryfile()),'(/,t4,"Grassman conjugate gradient solver:")')
          case (BD)
            write(x_unit(diaryfile()),'(/,t4,"Block Davidson solver:")')
          end select
          select case (es%o%rf_mode)
          case (ENABLED)
            select case (es%o%dir_profile)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t6,"first iteration:")')
              write(x_unit(diaryfile()),'(t8,"maximum directions = ",i0)') es%o%rfm_dir
              write(x_unit(diaryfile()),'(t6,"subsequent iterations:")')
              write(x_unit(diaryfile()),'(t8,"maximum directions = ",i0)') es%o%dir(1)
            case (STEPPED)
              write(x_unit(diaryfile()),'(/,t6,"first iteration:")')
              write(x_unit(diaryfile()),'(t8,"maximum directions = ",i0)') es%o%rfm_dir
              write(x_unit(diaryfile()),'(t6,"subsequent iterations:")')
              write(x_unit(diaryfile()),'(t8,"initial maximum directions = ",i0)') es%o%dir(2)
              do is = 3,size(es%o%dir)
                if (es%o%dir(is) /= es%o%dir(is-1)) then
                  write(x_unit(diaryfile()),'(t8,"new maximum directions     = ",i0," at step ",i0)') es%o%dir(is), is
                end if
              end do
            end select
            select case (es%o%tol_profile)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t6,"first iteration:")')
              write(x_unit(diaryfile()),'(t8,"residual tolerance = ",es8.2)') es%o%rfm_tol
              write(x_unit(diaryfile()),'(t6,"subsequent iterations:")')
              write(x_unit(diaryfile()),'(t8,"residual tolerance = ",es8.2)') es%o%tol(1)
            case (STEPPED)
              write(x_unit(diaryfile()),'(/,t6,"first iteration:")')
              write(x_unit(diaryfile()),'(t8,"residual tolerance = ",es8.2)') es%o%rfm_tol
              write(x_unit(diaryfile()),'(t6,"subsequent iterations:")')
              write(x_unit(diaryfile()),'(t8,"initial residual tolerance = ",es8.2)') es%o%tol(2)
              do is = 3,size(es%o%tol)
                if (es%o%tol(is) /= es%o%tol(is-1)) then
                  write(x_unit(diaryfile()),'(t8,"new residual tolerance     = ",es8.2," at step ",i0)') es%o%tol(is), is
                end if
              end do
            end select
          case (NOT_ENABLED)
            select case (es%o%dir_profile)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t6,"maximum directions = ",i0)') es%o%dir(1)
            case (STEPPED)
              write(x_unit(diaryfile()),'(/,t6,"initial maximum directions = ",i0)') es%o%dir(1)
              do is = 2,size(es%o%dir)
                if (es%o%dir(is) /= es%o%dir(is-1)) then
                  write(x_unit(diaryfile()),'(t6,"new maximum directions     = ",i0," at step ",i0)') es%o%dir(is), is
                end if
              end do
            end select
            select case (es%o%tol_profile)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t6,"residual tolerance = ",es8.2)') es%o%tol(1)
            case (STEPPED)
              write(x_unit(diaryfile()),'(/,t6,"initial residual tolerance = ",es8.2)') es%o%tol(1)
              do is = 2,size(es%o%tol)
                if (es%o%tol(is) /= es%o%tol(is-1)) then
                  write(x_unit(diaryfile()),'(t6,"new residual tolerance     = ",es8.2," at step ",i0)') es%o%tol(is), is
                end if
              end do
            end select
          end select
          select case (es%o%diagonalization_method)
          case (LAPACK)
            write(x_unit(diaryfile()),'(/,t4,"LAPACK diagonalization")')
          case (SCALAPACK)
            write(x_unit(diaryfile()),'(/,t4,"ScaLAPACK diagonalization")')
          end select
        end if
        call glean(thy(es))
      end subroutine

! local routines

      subroutine own_i(es)
        type(eigensolver_obj) :: es
        type(eigensolver_obj) :: est
        if (es%ref < es%o%ref) then
          allocate( est%o )
          est%o%ref = 0
          est%o%g = es%o%g
          est%o%report = es%o%report
          if (es%o%report) call my(es%o%f,est%o%f)
          est%o%solver_method = es%o%solver_method
          est%o%diagonalization_method = es%o%diagonalization_method
          est%o%dir_profile = es%o%dir_profile
          est%o%tol_profile = es%o%tol_profile
          est%o%max_steps = es%o%max_steps
          est%o%step = es%o%step
          est%o%rf_mode = es%o%rf_mode
          select case (es%o%rf_mode)
          case (ENABLED)
            est%o%rfm_status = es%o%rfm_status
            est%o%rfm_dir = es%o%rfm_dir
            est%o%rfm_tol = es%o%rfm_tol
          end select
          allocate( est%o%dir(size(es%o%dir)) ) ; est%o%dir = es%o%dir
          allocate( est%o%tol(size(es%o%tol)) ) ; est%o%tol = es%o%tol
          es%o%ref = es%o%ref - es%ref
          es%o => est%o
          es%o%ref = es%o%ref + es%ref
        end if
      end subroutine

      subroutine solve_cg_i(es,h,v,evals,rnorm)
        type(eigensolver_obj) :: es
        type(hamiltonian_obj) :: h
        type(multivector_obj) :: v
        real(double), dimension(:), intent(inout) :: evals
        real(double), intent(out) :: rnorm

        logical :: converged
        integer :: ib, id, max_dir, nb, nd
        real(double) :: res_tol, theta_max
        complex(double) :: cm1, cp1, cz
        real(double), dimension(:), allocatable :: chc, dhd, re_dhc, dl, gl, rl, gamma, gdotg
        real(double), dimension(:), allocatable :: de_t, de_tt, e_s1, e_c1, d1e_0, d2e_0, theta
        real(double), dimension(:), allocatable :: rv1, rv2, rvm1, rvz, rvp1
        complex(double), dimension(:), allocatable :: cv, cvm1, cvz, cvp1
        complex(double), dimension(:,:), allocatable :: cmat
        type(multivector_obj) :: c, hc, d, hd, f

        call my(es)
        call my(h)
        call my(v)

        nb = x_n_bands(v)

        allocate( chc(nb), dhd(nb), re_dhc(nb), dl(nb), gl(nb), rl(nb), gamma(nb), gdotg(nb) )
        allocate( de_t(nb), de_tt(nb), e_s1(nb), e_c1(nb), d1e_0(nb), d2e_0(nb), theta(nb) )
        allocate( rv1(nb), rv2(nb), rvm1(nb), rvz(nb), rvp1(nb), cv(nb), cvm1(nb), cvz(nb), cvp1(nb) )
        allocate( cmat(nb,nb) )

        cm1 = cmplx(-1,0,double)
        cz  = cmplx( 0,0,double)
        cp1 = cmplx(+1,0,double)
        rvm1 = -1.0_double
        rvz  =  0.0_double
        rvp1 = +1.0_double
        cvm1 = cmplx(-1,0,double)
        cvz  = cmplx( 0,0,double)
        cvp1 = cmplx(+1,0,double)
        theta_max = 1.0_double

        call my(v,c)  ! Note: real copies are made here (with memory allocated)
        call my(v,hc)
        call my(v,d)
        call my(v,hd)
        call my(v,f)

        call apply_hamiltonian(c,hc,h) ! H*c --> hc
        call multiply(c,hc,chc)  ! <c|hc> --> chc
        call combine(chc,c,rvm1,hc,d)  ! chc*c - hc --> d

        call get_parameters_i(es%o,max_dir,res_tol)

        nd = 0
        do id = 1,max_dir

          nd = nd + 1

          call precondition(d,c,h)  ! precondition(d,{c,H}) --> d
          call overlap(v,d,cmat)  ! (v' x d) --> cmat
          call transform(cp1,d,cm1,v,cmat)  ! 1.0*d - 1.0*(v x cmat) --> d
          call multiply(d,gl) ; gl = sqrt(gl)  ! sqrt(,d|d>) --> gl

          select case (id)
          case (1)
            f = d
            dl = gl
            gdotg = gl*gl
          case default
            gamma = gl*gl/gdotg
            call combine(gamma,f,rvp1,d)  ! gamma*f + 1.0*d --> f
            gdotg = gl*gl
            call multiply(c,f,cv)  ! <c|f> --> cv
            cv = -cv
            call combine(cvp1,f,cv,c,d)  ! (1.0,0.0)*f + cv*c --> d
            call multiply(d,dl) ; dl = sqrt(dl)  ! sqrt(<d|d>) --> dl
          end select

          rv1 = 1.0_double/dl
          call portion(rv1,d)  ! d/dl --> d

          call apply_hamiltonian(d,hd,h)  ! H*d --> hd
          call multiply(d,hd,dhd)  ! <d|hd> --> dhd
          call multiply(d,hc,re_dhc)  ! Real(<d|hc>) --> re_dhc

          do ib = 1,nb
            d1e_0(ib) = 2.0_double*re_dhc(ib)
            d2e_0(ib) = 2.0_double*(dhd(ib) - chc(ib))
            e_s1(ib) =  d1e_0(ib)/2.0_double
            e_c1(ib) = -d2e_0(ib)/4.0_double
            if (e_c1(ib) == 0.0_double) then
              if (e_s1(ib) == 0.0_double) then
                theta(ib) = 0.0_double
              else
                theta(ib) = pi/2.0_double
              end if
            else
              theta(ib) = (atan(e_s1(ib)/e_c1(ib)))/2.0_double
              de_t(ib)  = -e_c1(ib) + e_c1(ib)*cos(2.0_double*theta(ib)) + e_s1(ib)*sin(2.0_double*theta(ib))
              de_tt(ib) = -e_c1(ib) + e_c1(ib)*cos(2.0_double*theta(ib) + pi) + e_s1(ib)*sin(2.0_double*theta(ib) + pi)
              if (de_tt(ib) < de_t(ib)) then
                theta(ib) = theta(ib) + ( pi/2.0_double )
                de_t(ib) = de_tt(ib)
              end if
            end if
          end do

          rv1 = cos(theta) ; rv2 = sin(theta)
          call combine(rv1,c,rv2,d)  ! rv1*c + rv2*d --> c
          call combine(rv1,hc,rv2,hd)  ! rv1*hc + rv2*hd --> hc
          chc = chc*rv1**2 + dhd*rv2**2 + 2.0_double*re_dhc*rv1*rv2

          call combine(chc,c,rvm1,hc,d)  ! chc*c - hc --> d
          call multiply(d,rl) ; rl = sqrt(rl)  ! sqrt(<d|d>) --> rl

          converged = .true.
          do ib = 1,nb
            converged = (converged .and. (rl(ib) < res_tol) )
          end do
          if (converged) exit

        end do

        rnorm = sum(rl)/real(nb,double)

        call overlap(c,cmat)
        call inverse_cholesky(cmat) ; if (error()) goto 100
        call transform(cz,d,cp1,c,cmat)
        call transform(cz,hd,cp1,hc,cmat)

        call overlap(d,hd,cmat)
        select case (es%o%diagonalization_method)
        case (LAPACK)
          call diagonalize_lapack_i(cmat,evals)
        case (SCALAPACK)
          call diagonalize_scalapack_i(cmat,evals)
        end select
        if (error()) goto 100
        call transform(cz,v,cp1,d,cmat)

        if (es%o%report) call report_i(es%o,max_dir,nd,res_tol,rnorm)

100     deallocate( chc, dhd, re_dhc, dl, rl, gl, gamma, gdotg )
        deallocate( de_t, de_tt, e_s1, e_c1, d1e_0, d2e_0, theta )
        deallocate( rv1, rv2, rvm1, rvz, rvp1, cv, cvm1, cvz, cvp1 )
        deallocate( cmat )

        call glean(thy(c))
        call glean(thy(hc))
        call glean(thy(d))
        call glean(thy(hd))
        call glean(thy(f))

        call glean(thy(es))
        call glean(thy(h))
        call glean(thy(v))

        if (error("Exit eigensolver_mod::solve_cg_i")) continue

      end subroutine

      subroutine solve_gcg_i(es,h,v,evals,rnorm)
        type(eigensolver_obj) :: es
        type(hamiltonian_obj) :: h
        type(multivector_obj) :: v
        real(double), dimension(:), intent(inout) :: evals
        real(double), intent(out) :: rnorm

        integer :: id, istop, max_dir, nb
        real(double) :: a, b, c, dmag, rho, rho_old, t, tol, x
        complex(double) :: c1, c2
        complex(double), dimension(:,:), allocatable :: vhv, tmp
        complex(double), dimension(:,:), allocatable :: d, dhd, dd, dg
        complex(double), dimension(:,:), allocatable :: miniH
        type(multivector_obj) :: tmpv, dir

        call my(es)
        call my(h)
        call my(v)

        nb = x_n_bands(v)

        call my(v,dir)  ! Note: real copies are made here (with memory allocated)
        call my(v,tmpv)

        allocate( d(nb,nb) )     ;     d = cmplx(0,0,double)
        allocate( vhv(nb,nb) )   ;   vhv = cmplx(0,0,double)
        allocate( miniH(nb,nb) ) ; miniH = cmplx(0,0,double)
        allocate( tmp(nb,nb) )   ;   tmp = cmplx(0,0,double)
        allocate( dhd(nb,nb) )   ;   dhd = cmplx(0,0,double)
        allocate( dd(nb,nb) )    ;    dd = cmplx(0,0,double)
        allocate( dg(nb,nb) )    ;    dg = cmplx(0,0,double)

        ! get first orthotransformation
        call overlap(v,d)
        call inverse_cholesky(d) ; if (error()) goto 100

        ! apply h and project (finding vhv) to get new direction
        call apply_hamiltonian(v,dir,h) ! dir = H*v
        call overlap(v,dir,vhv)   ! vhv = v*dir
        c1 = cmplx(1,0,double) ; c2 = cmplx(-1,0,double)
        call transform(c1,dir,c2,v,matmul(d,matmul(transpose(conjg(d)),vhv))) ! dir = dir - v*d*d'*vhv
         
        ! get new residual
        rho_old = 1.0_double
        call overlap(dir,dd) ! dd = dir'*dir
        rho = trace(real(matmul(transpose(conjg(d)),matmul(dd,d)),double))
        dmag = sqrt(rho)

        call get_parameters_i(es%o,max_dir,tol)

        istop = 0
        do id = 1,max_dir

         if (sqrt(rho) <= tol) exit

          istop = istop + 1
         
          if (mod(istop,10) == 1) then
!            write(iobuf,*) "Grassmann iteration ",istop; call diary
!            write(iobuf,*) 'wave residual :',sqrt(rho) ; call diary

            ! make Ritz approximations to report
            miniH = matmul(transpose(conjg(d)),matmul(vhv,d))
            select case (es%o%diagonalization_method)
            case (LAPACK)
              call diagonalize_lapack_i(miniH,evals)
            case (SCALAPACK)
              call diagonalize_scalapack_i(miniH,evals)
            end select
            if (error()) goto 100
!            write(iobuf,*) 'Ritz values :', evals ; call diary
          end if

          ! compute a,b,c coefficients for linesearch approx
          call apply_hamiltonian(dir,tmpv,h) ! tmpv = H*dir
          call overlap(dir,tmpv,dhd) ! dhd = dir'*tmpv
          a = trace(real(matmul(transpose(conjg(d)),matmul(dhd,d)),double))
          call overlap(v,tmpv,dg) ! dg = v'*tmpv
          b = trace(real(matmul(transpose(conjg(d)),matmul(dg,d)),double))
          miniH = matmul(transpose(conjg(d)),matmul(vhv,d)) ! (because miniH gets clobbered by the diagonalization)
          tmp = matmul(transpose(conjg(d)),matmul(dd,d))
          c = trace(real(matmul(tmp,miniH),double))
          t = b/abs(a - c) ! line search approximation.
          t = t*dmag; ! a hack to keep it compact
          t = t/sqrt(1+t**2)
          t = t/dmag;

          ! do updates to v in place and to dir (stored in tmpv), muck about in order to do
          c1 = cmplx(0,0,double) ; c2 = cmplx(1,0,double)
          call combine(c1,tmpv,c2,dir) ! tmpv = dir + t*v*d'*d*dir'*dir
          c1 = cmplx(1,0,double) ; c2 = cmplx(t,0,double)
          call transform(c1,tmpv,c2,v,matmul(d,matmul(transpose(conjg(d)),dd)) )
          c1 = cmplx(1,0,double) ; c2 = cmplx(-t,0,double)
          call combine(c1,v,c2,dir) ! v = v + t * dir

          ! get new orthotransformation
          call overlap(v,d)
          call inverse_cholesky(d) ; if (error()) goto 100

          ! apply h and project (finding vhv) to get new direction
          call apply_hamiltonian(v,dir,h) ! dir = H*v
          call overlap(v,dir,vhv) ! vhv = v*dir
          c1 = cmplx(1,0,double) ; c2 = cmplx(-1,0,double)
          call transform(c1,dir,c2,v,matmul(d,matmul(transpose(conjg(d)),vhv)) ) ! dir = dir - v*d*d'*vhv
         
          ! get new residual
          rho_old = rho 
          call overlap(dir,dd) ! dd = dir'*dir
          rho = trace(real(matmul(transpose(conjg(d)),matmul(dd,d)),double))
          dmag = sqrt(rho)

          ! do check to make sure new grad is approx perp to old search direction
          ! this determines whether we are in the linear regime of the operator
          call overlap(dir,tmpv,dg) ! dg = dir'*tmpv
          x = trace(real(matmul(transpose(conjg(d)),matmul(dg,d)),double))

          ! if we are in the linear regime then conjugate the gradient
          if (abs(x) < 1.0e-2_double ) then
            c1 = cmplx(1,0,double) ; c2 = cmplx(rho/rho_old,0,double)
            call combine(c1,dir,c2,tmpv) ! dir = dir + rho/rho_old * tmpv
            call overlap(dir,dd) ! dd = dir'*dir
            dmag = sqrt(trace(real(matmul(transpose(conjg(d)),matmul(dd,d)),double)))
          end if

        end do

        ! perform re-orthogonalization
        call transform(v,d)

        ! get final stats
        miniH = matmul(transpose(conjg(d)),matmul(vhv,d))
        select case (es%o%diagonalization_method)
        case (LAPACK)
          call diagonalize_lapack_i(miniH,evals)
        case (SCALAPACK)
          call diagonalize_scalapack_i(miniH,evals)
        end select
        if (error()) goto 100
        call transform(v,miniH)
!        write(iobuf,*) '% of iterations :',istop ; call diary
!        call diary("Final values:")
!        write(iobuf,*) 'subspace eigenvalues :',evals ; call diary
!        write(iobuf,*) 'wave residual :', sqrt(rho) ; call diary

        rnorm = sqrt(rho)

100     deallocate( tmp, dhd, dd, dg, vhv, d, miniH )

        call glean(thy(tmpv))
        call glean(thy(dir))

        call glean(thy(es))
        call glean(thy(h))
        call glean(thy(v))

        if (error("Exit eigensolver_mod::solve_gcg_i")) continue

      end subroutine

      subroutine solve_bd_i(es,h,v,evals,rnorm)
        type(eigensolver_obj) :: es
        type(hamiltonian_obj) :: h
        type(multivector_obj) :: v
        real(double), dimension(:), intent(inout) :: evals
        real(double), intent(out) :: rnorm

        logical :: converged
        integer :: ib, id, max_dir, nb, nd
        real(double) :: res_tol
        real(double), dimension(:), allocatable :: res, res_norm, evals_red 
        complex(double) :: cz, cp1 
        complex(double), dimension(:), allocatable :: evals_cmplx
        complex(double), dimension(:,:), allocatable :: cmat, h_mat, o_mat
        type(multivector_obj) :: hv, r, hr  

        call my(es)
        call my(h)
        call my(v)

        cz = cmplx(0,0,double)
        cp1 = cmplx(1,0,double)

        nb = x_n_bands(v)

        allocate( res(nb) )
        allocate( res_norm(nb) )
        allocate( evals_red(2*nb) )
        allocate( evals_cmplx(nb) )
        allocate( cmat(nb,nb) )
        allocate( h_mat(2*nb,2*nb) )
        allocate( o_mat(2*nb,2*nb) )

        call my(v,hv)  ! Note: real copies are made here (with memory allocated)
        call my(v,r)
        call my(v,hr)

        call apply_hamiltonian(v,hv,h) ; if (error()) goto 100  ! H*v --> hv
        call multiply(v,hv,evals)
        evals_cmplx = cmplx(evals,0,double)
        call residual(v,hv,evals_cmplx,r) ; if (error()) goto 100
        call multiply(r,res) ; if (error()) goto 100
        res = sqrt(res)

        call get_parameters_i(es%o,max_dir,res_tol)

        nd = 0
        do id = 1,max_dir

          nd = nd + 1

          res_norm = 1.0_double
          do ib = 1,nb
            if (machine_zero < res(ib)) res_norm(ib) = 1.0_double/res(ib)
          end do
          call portion(res_norm,r) ; if (error()) goto 100

          call precondition(r,v,h) ; if (error()) goto 100  ! precondition(r,{v,H}) --> r

          ! Construct the reduced hamiltonian. The reduced hamiltonian has dimensions
          !  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
          !             _                  _ 
          !            |                    |
          !            | <v|H|v>   <v|H|r>  |
          !    h_mat = |                    |
          !            | *******   <r|H|r>  | 
          !            |_                  _|

          call apply_hamiltonian(r,hr,h) ; if (error()) goto 100  ! H*r --> hr
          if (id == 1) then
            call overlap(v,hv,cmat) ; if (error()) goto 100  ! <v|H|v> --> cmat
            h_mat(1:nb,1:nb) = cmat
          else
            h_mat(1:nb,1:nb) = cmplx(0,0,double)
            do ib = 1,nb
              h_mat(ib,ib) = cmplx(evals(ib),0,double)
            end do
          end if
          call overlap(v,hr,cmat) ; if (error()) goto 100  ! <v|H|r> --> cmat
          h_mat(1:nb,nb+1:2*nb) = cmat
          call overlap(r,hr,cmat) ; if (error()) goto 100  ! <r|H|r> --> cmat
          h_mat(nb+1:2*nb,nb+1:2*nb) = cmat

          ! Construct the reduced overlap matrix which has dimenstions 2nb x 2nb
          !   and is constructed by filling in four nb x nb blocks one at a time:
          !             _               _ 
          !            |                 |
          !            |  <v|v>   <v|r>  |
          !    o_mat = |                 |
          !            |  *****   <r|r>  | 
          !            |_               _|

          o_mat(1:nb,1:nb) = cmplx(0,0,double)
          do ib = 1,nb
            o_mat(ib,ib) = cmplx(1,0,double)
          end do
          call overlap(v,r,cmat) ; if (error()) goto 100  ! <v|r> --> cmat
          o_mat(1:nb,nb+1:2*nb) = cmat
          call overlap(r,cmat) ; if (error()) goto 100  ! <r|r> --> cmat
          o_mat(nb+1:2*nb,nb+1:2*nb) = cmat

          select case (es%o%diagonalization_method)
          case (LAPACK)
            call diagonalize_gen_lapack_i(h_mat,o_mat,evals_red)
          case (SCALAPACK)
            call diagonalize_gen_scalapack_i(h_mat,o_mat,evals_red)
          end select
          if (error()) goto 100

          evals = evals_red(1:nb)
          
          cmat = h_mat(1:nb,1:nb)
          call transform(v,cmat) ; if (error()) goto 100
          call transform(hv,cmat) ; if (error()) goto 100
          cmat = h_mat(nb+1:2*nb,1:nb)
          call transform(cp1,v,cp1,r,cmat) ; if (error()) goto 100
          call transform(cp1,hv,cp1,hr,cmat) ; if (error()) goto 100

          evals_cmplx = cmplx(evals,0,double)
          call residual(v,hv,evals_cmplx,r)
          
          call multiply(r,res) ; if (error()) goto 100
          res = sqrt(res)

          converged = .true.
          do ib = 1,nb
            converged = (converged .and. (res(ib) < res_tol) )
          end do
          if (converged) exit

        end do

        rnorm = sum(res)/real(nb,double)

        if (es%o%report) call report_i(es%o,max_dir,nd,res_tol,rnorm)

100     if (allocated( res )) deallocate( res )
        if (allocated( res_norm )) deallocate( res_norm )
        if (allocated( evals_red )) deallocate( evals_red )
        if (allocated( evals_cmplx )) deallocate( evals_cmplx )
        if (allocated( cmat )) deallocate( cmat )
        if (allocated( h_mat )) deallocate( h_mat )
        if (allocated( o_mat )) deallocate( o_mat )

        call glean(thy(hv))
        call glean(thy(r))
        call glean(thy(hr))

        call glean(thy(es))
        call glean(thy(h))
        call glean(thy(v))

        if (error("Exit eigensolver_mod::solve_bd_i")) continue

      end subroutine

      subroutine solve_bdgen_i(es,h,v,evals,rnorm)
        type(eigensolver_obj) :: es
        type(hamiltonian_obj) :: h
        type(multivector_obj) :: v
        real(double), dimension(:), intent(inout) :: evals
        real(double), intent(out) :: rnorm

        logical :: converged
        integer :: ib, id, max_dir, nb, nd
        real(double) :: res_tol
        real(double), dimension(:), allocatable :: res, res_norm, evals_red
        complex(double) :: cz, cp1
        complex(double), dimension(:), allocatable :: evals_cmplx
        complex(double), dimension(:,:), allocatable :: cmat, h_mat, s_mat
        type(multivector_obj) :: hv, ov, r, hr, or

        call my(es)
        call my(h)
        call my(v)

        cz  = ( 0.0_double,0.0_double)
        cp1 = (+1.0_double,0.0_double)

        nb = x_n_bands(v)

        allocate( res(nb) )
        allocate( res_norm(nb) )
        allocate( evals_red(2*nb) )
        allocate( evals_cmplx(nb) )
        allocate( cmat(nb,nb) )
        allocate( h_mat(2*nb,2*nb) )
        allocate( s_mat(2*nb,2*nb) )

        call my(v,hv)  ! Note: real copies are made here (with memory allocated)
        call my(v,ov)  
        call my(v,r)
        call my(v,hr)
        call my(v,or)

        call apply_hamiltonian(v,hv,h) ! H*v --> hv

        call apply_overlap(v,ov,h) ! O*v --> ov

        call multiply(v,hv,evals)
        evals_cmplx = evals
        call multiply(v,ov,evals)
        evals_cmplx = evals_cmplx/evals
        call residual(ov,hv,evals_cmplx,r)  ! eval*ov - hv --> r

        call multiply(r,res)
        res = sqrt(res)

        call get_parameters_i(es%o,max_dir,res_tol)

        nd = 0
        do id = 1,max_dir

          nd = nd + 1

          res_norm = 1.0_double
          do ib = 1,nb
            if (machine_zero < res(ib)) res_norm(ib) = 1.0_double/res(ib)
          end do
          call portion(res_norm,r) ; if (error()) goto 100

          call precondition(r,v,h)  ! precondition(r,{v,H}) --> r

          ! Construct the reduced hamiltonian. The reduced hamiltonian has dimensions of 
          !  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
          !
          !             _                  _ 
          !            |                    |
          !            | <v|H|v>   <v|H|r>  |
          !    h_mat = |                    |
          !            | <r|H|v>   <r|H|r>  | 
          !            |_                  _|
          !

          call apply_hamiltonian(r,hr,h) ! H*r --> hr

          if (id == 1) then
             call overlap(v,hv,cmat) ! <v|H|v> --> cmat
             h_mat(1:nb,1:nb) = cmat(1:nb,1:nb)
          else
             h_mat(1:nb,1:nb) = cmplx(0,0,double)
             do ib=1,nb
                h_mat(ib,ib) = cmplx(evals(ib),0,double)
             end do
          end if

          call overlap(v,hr,cmat) ! <v|H|r> --> cmat
          h_mat(1:nb,1+nb:2*nb) = cmat(1:nb,1:nb)

          call overlap(r,hr,cmat) ! <r|H|r> --> cmat
          h_mat(1+nb:2*nb,1+nb:2*nb) = cmat(1:nb,1:nb)

          ! Construct the reduced overlap matrix which has dimenstions of
          !  2nb x 2nb and is constructed by filling in four nb x nb blocks one at a time:
          !
          !             _                   _ 
          !            |                     |
          !            |  <v|O|v>   <v|O|r>  |
          !    s_mat = |                     |
          !            |  <r|O|v>   <r|O|r>  | 
          !            |_                   _|
          !

          call apply_overlap(r,or,h) ! O*r --> or

          s_mat(1:nb,1:nb) = cmplx(0,0,double)
          do ib=1,nb
             s_mat(ib,ib) = cmplx(1,0,double)
          end do

          call overlap(v,or,cmat) ! <v|or> --> cmat
          s_mat(1:nb,1+nb:2*nb) = cmat(1:nb,1:nb)

          call overlap(r,or,cmat) ! <r|or> --> cmat
          s_mat(1+nb:2*nb,1+nb:2*nb) = cmat(1:nb,1:nb)

          select case (es%o%diagonalization_method)
          case (LAPACK)
            call diagonalize_gen_lapack_i(h_mat,s_mat,evals_red)
          case (SCALAPACK)
            call diagonalize_gen_scalapack_i(h_mat,s_mat,evals_red)
          end select
          if (error()) goto 100

          evals(1:nb) = evals_red(1:nb)

          cmat = h_mat(1:nb,1:nb)
          call transform(v,  cmat)
          call transform(hv, cmat)
          call transform(ov, cmat)
          cmat = h_mat(1+nb:2*nb,1:nb)
          call transform(cp1, v,  cp1, r, cmat)
          call transform(cp1, hv, cp1, hr, cmat)
          call transform(cp1, ov, cp1, or, cmat)

          evals_cmplx = evals
          call residual(ov,hv,evals_cmplx,r) ! eval*ov - hv --> r   A.K.A. the residual

          call multiply(r,res) ; if (error()) goto 100
          res = sqrt(res)

          converged = .true.
          do ib = 1,nb
            converged = (converged .and. (res(ib) < res_tol) )
          end do
          if (converged) exit

        end do

        rnorm = sum(res)/real(nb,double)

        if (es%o%report) call report_i(es%o,max_dir,nd,res_tol,rnorm)

100     if (allocated( res )) deallocate( res )
        if (allocated( res_norm )) deallocate( res_norm )
        if (allocated( evals_red )) deallocate( evals_red )
        if (allocated( evals_cmplx )) deallocate( evals_cmplx )
        if (allocated( cmat )) deallocate( cmat )
        if (allocated( h_mat )) deallocate( h_mat )
        if (allocated( s_mat )) deallocate( s_mat )

        call glean(thy(hv))
        call glean(thy(ov))
        call glean(thy(r))
        call glean(thy(hr))
        call glean(thy(or))

        call glean(thy(es))
        call glean(thy(h))
        call glean(thy(v))

        if (error("Exit eigensolver_mod::solve_bdgen_i")) continue

      end subroutine

      subroutine diagonalize_lapack_i(a,w)
        complex(double), dimension(:,:), intent(inout) :: a
        real(double), dimension(:), intent(inout) :: w

        logical :: i_participate
        character(1) :: jobz, uplo
        integer, save :: n0 = 0
        integer, save :: lwork
        integer :: info, n
        real(double), dimension(:), allocatable :: rwork
        complex(double), dimension(:), allocatable :: work

        call start_timer("eigensolver: diagonalize")

        i_participate = mpi_isroot(KGROUP)

        if (i_participate) then

          n = size(a,1)
          if (n /= n0) lwork = 65*n
          allocate( work(lwork), rwork(3*n-2) )

          jobz = 'v'
          uplo = 'u'
          call zheev(jobz,uplo,n,a,n,w,work,lwork,rwork,info)
          if (error(info /= 0,"ERROR: zheev returned error number ",info)) goto 100

          if (n /= n0) then
            lwork = int(real(work(1),double))
            n0 = n
          end if

        end if

100     call sync_kgroup_process_errors() ; if (error()) goto 200

        call broadcast(KGROUP,a)
        call broadcast(KGROUP,w)

        if (allocated( work )) deallocate( work )
        if (allocated( rwork )) deallocate( rwork )

200     if (error("Exit eigensolver_mod::diagonalize_lapack_i")) continue

        if (.not.error()) call stop_timer("eigensolver: diagonalize")

      end subroutine

      subroutine diagonalize_scalapack_i(a,w)
        complex(double), dimension(:,:), intent(inout) :: a
        real(double), dimension(:), intent(inout) :: w

        logical :: i_participate
        character(1) :: jobz, order, uplo
        integer :: comm, d1, d2, i0, i1, info, lwork, lrwork, na, nb, np, npc, npr, numroc, pc, pr
        integer, dimension(9) :: desc_a, desc_as
        complex(double), dimension(:), allocatable :: work, rwork
        complex(double), dimension(:,:), allocatable :: as, zs

        call start_timer("eigensolver: diagonalize")

        ! Determine processor grid
        np = mpi_nprocs(KGROUP)
        npr = int(sqrt(real(np)))
        npc = np/npr

        ! Create blacs context
        comm = mpi_comm(KGROUP)
        order = 'r'
        call blacs_gridinit(comm,order,npr,npc)
        call blacs_gridinfo(comm,npr,npc,pr,pc)

        i_participate = ((pr >= 0) .and. (pr < npr) .and. (pc >= 0) .and. (pc < npc))

        if (i_participate) then

          i0 = 0
          i1 = 1

          jobz = 'v'
          uplo = 'u'
          na = size(a,1)
          nb = blocksize

          ! Set up descriptors and matrices
          d1 = numroc(na,nb,pr,i0,npr)
          d2 = numroc(na,nb,pc,i0,npc)
          allocate( as(d1,d2), zs(d1,d2) )
          call descinit(desc_as,na,na,nb,nb,i0,i0,comm,max(1,d1),info)
          if (error(info /= 0,"ERROR: first call to descinit returned error number ",info)) goto 100
          call descinit(desc_a,na,na,na,na,i0,i0,comm,na,info)
          if (error(info /= 0,"ERROR: second call to descinit returned error number ",info)) goto 100
          call pzgemr2d(na,na,a,i1,i1,desc_a,as,i1,i1,desc_as,comm)

          ! Allocate work space
          lwork = -1 ; lrwork = -1
          allocate( work(1), rwork(1) )
          call pzheev(jobz,uplo,na,as,i1,i1,desc_as,w,zs,i1,i1,desc_as,work,lwork,rwork,lrwork,info)
          if (error(info /= 0,"ERROR: first call to pzheev returned error number ",info)) goto 100
          lwork = int(real(work(1),double))
          lrwork = int(real(rwork(1),double))
          deallocate( work, rwork )
          allocate( work(lwork), rwork(lrwork) )

          ! Diagonalize matrix
          call pzheev(jobz,uplo,na,as,i1,i1,desc_as,w,zs,i1,i1,desc_as,work,lwork,rwork,lrwork,info)
          if (error(info /= 0,"ERROR: second call to pzheev returned error number ",info)) goto 100
          call pzgemr2d(na,na,zs,i1,i1,desc_as,a,i1,i1,desc_a,comm)

          call blacs_gridexit(comm)

        end if

100     call sync_kgroup_process_errors() ; if (error()) goto 200

        call broadcast(KGROUP,a)
        call broadcast(KGROUP,w)

200     if (allocated( as )) deallocate( as )
        if (allocated( zs )) deallocate( zs )
        if (allocated( work )) deallocate( work )
        if (allocated( rwork )) deallocate( rwork )

        if (error("Exit eigensolver_mod::diagonalize_scalapack_i")) continue

        if (.not.error()) call stop_timer("eigensolver: diagonalize")

      end subroutine

      subroutine diagonalize_gen_lapack_i(a,b,w)
        complex(double), dimension(:,:), intent(inout) :: a
        complex(double), dimension(:,:), intent(inout) :: b
        real(double), dimension(:), intent(inout) :: w

        logical :: i_participate
        character(1) :: jobz, uplo
        integer, save :: n0 = 0
        integer, save :: lwork
        integer :: info, n, type
        real(double), dimension(:), allocatable :: rwork
        complex(double), dimension(:), allocatable :: work

        call start_timer("eigensolver: diagonalize_gen")

        i_participate = mpi_isroot(KGROUP)

        if (i_participate) then

          n = size(a,1)
          if (n /= n0) lwork = 65*n
          allocate( work(lwork), rwork(3*n-2) )

          type = 1
          jobz = 'v'
          uplo = 'u'
          call zhegv(type,jobz,uplo,n,a,n,b,n,w,work,lwork,rwork,info)
          if (error(info /= 0,"ERROR: zhegv returned error number ",info)) goto 100

          if (n /= n0) then
            lwork = int(real(work(1),double))
            n0 = n
          end if

        end if

100     call sync_kgroup_process_errors() ; if (error()) goto 200

        call broadcast(KGROUP,a)
        call broadcast(KGROUP,w)

        if (allocated( work )) deallocate( work )
        if (allocated( rwork )) deallocate( rwork )

200     if (error("Exit eigensolver_mod::diagonalize_gen_lapack_i")) continue

        if (.not.error()) call stop_timer("eigensolver: diagonalize_gen")

      end subroutine

      subroutine diagonalize_gen_scalapack_i(a,b,w)
        complex(double), dimension(:,:), intent(inout) :: a, b
        real(double), dimension(:), intent(inout) :: w

        logical :: found, i_participate
        character(1) :: cmach, jobz, order, range, uplo
        integer :: i0, i1, comm, d1, d2, info, lwork, lrwork, liwork, mf, na, nb, np, npc, npr, pc, pr, type, nec, nef, numroc
        integer, dimension(9) :: desc_a, desc_as
        integer, dimension(:), allocatable :: clustr, fail, iwork
        real(double) :: abstol, orfac, pdlamch, r0
        real(double), dimension(:), allocatable :: gap
        complex(double), dimension(:), allocatable :: work, rwork
        complex(double), dimension(:,:), allocatable :: as, bs, zs

        call start_timer("eigensolver: diagonalize_gen")

        ! Create processor grid
        np = mpi_nprocs(KGROUP)
        npr = int(sqrt(real(np)))
        npc = np/npr

        ! Create blacs context
        order = 'r'
        comm = mpi_comm(KGROUP)
        call blacs_gridinit(comm,order,npr,npc)
        call blacs_gridinfo(comm,npr,npc,pr,pc)
        i_participate = ((pr >= 0) .and. (pr < npr) .and. (pc >= 0) .and. (pc < npc))

        if (i_participate) then

          i0 = 0
          i1 = 1
          r0 = 0.0_double

          type = 1
          jobz = 'v'
          range = 'a'
          uplo = 'u'
          na = size(a,1)
          nb = blocksize

          ! Set up descriptors and matrices
          d1 = numroc(na,nb,pr,i0,npr)
          d2 = numroc(na,nb,pc,i0,npc)
          allocate( as(d1,d2), bs(d1,d2), zs(d1,d2) )
          call descinit(desc_as,na,na,nb,nb,i0,i0,comm,max(1,d1),info)
          if (error(info /= 0,"ERROR: first call to descinit returned error number ",info)) goto 100
          call descinit(desc_a,na,na,na,na,i0,i0,comm,na,info)
          if (error(info /= 0,"ERROR: second call to descinit returned error number ",info)) goto 100
          call pzgemr2d(na,na,a,i1,i1,desc_a,as,i1,i1,desc_as,comm)
          call pzgemr2d(na,na,b,i1,i1,desc_a,bs,i1,i1,desc_as,comm)

          allocate( fail(na), clustr(2*npr*npc), gap(npr*npc) )

          ! Allocate work arrays
          lwork = -1 ; lrwork = -1 ; liwork = -1
          allocate( work(1), rwork(1), iwork(1))
          call pzhegvx(type,jobz,range,uplo,na,as,i1,i1,desc_as,bs,i1,i1,desc_as,r0,r0,i0,i0,abstol,nef,nec,w,orfac, &
                       & zs,i1,i1,desc_as,work,lwork,rwork,lrwork,iwork,liwork,fail,clustr,gap,info)
          if (error(info /= 0,"ERROR: first call to pzhegvx returned error number ",info)) goto 100
!          call notify("work(1) return value",int(real(work(1),double)))
!          call notify("rwork(1) return value",int(real(rwork(1),double)))
!          call notify("iwork(1) return value",iwork(1))
          call arg("memory_factor",mf,found)
          if (.not.found) mf = 100
          lwork = mf*int(real(work(1),double))
          lrwork = mf*int(real(rwork(1),double))
          liwork = mf*iwork(1)
          deallocate( work, rwork, iwork )
          allocate( work(lwork), rwork(lrwork), iwork(liwork))

          ! Diagonalize matrix
          cmach = 's'
          abstol = 2*pdlamch(comm,cmach)
          orfac = -1.0_double
          call pzhegvx(type,jobz,range,uplo,na,as,i1,i1,desc_as,bs,i1,i1,desc_as,r0,r0,i0,i0,abstol,nef,nec,w,orfac, &
                       & zs,i1,i1,desc_as,work,lwork,rwork,lrwork,iwork,liwork,fail,clustr,gap,info)
          if (error(info /= 0,"ERROR: second call to pzhegvx returned error number ",info)) goto 100
          call pzgemr2d(na,na,zs,i1,i1,desc_as,a,i1,i1,desc_a,comm)

          call blacs_gridexit(comm)

        end if

100     call sync_kgroup_process_errors() ; if (error()) goto 200

        call broadcast(KGROUP,a)
        call broadcast(KGROUP,w)

200     if (allocated( as )) deallocate( as )
        if (allocated( bs )) deallocate( bs )
        if (allocated( zs )) deallocate( zs )
        if (allocated( work )) deallocate( work )
        if (allocated( rwork )) deallocate( rwork )
        if (allocated( iwork )) deallocate( iwork )
        if (allocated( fail )) deallocate( fail )
        if (allocated( clustr )) deallocate( clustr )
        if (allocated( gap )) deallocate( gap )

        if (error("Exit eigensolver_mod::diagonalize_gen_scalapack_i")) continue

        if (.not.error()) call stop_timer("eigensolver: diagonalize_gen")

      end subroutine

      subroutine get_parameters_i(esr,dir,tol)
        type(eigensolver_rep) :: esr
        integer, intent(out) :: dir
        real(double), intent(out) :: tol

        select case (esr%rf_mode)
        case (ENABLED)
          select case(esr%rfm_status)
          case (PENDING)
            dir = esr%rfm_dir
            tol = esr%rfm_tol
            esr%rfm_status = DONE
          case (DONE)
            dir = esr%dir(esr%step)
            tol = esr%tol(esr%step)
          end select
        case (NOT_ENABLED)
          dir = esr%dir(esr%step)
          tol = esr%tol(esr%step)
        end select

      end subroutine

      subroutine report_i(esr,max_dir,num_dir,res_tol,res_norm)
        type(eigensolver_rep) :: esr
        integer, intent(in) :: max_dir, num_dir
        real(double), intent(in) :: res_tol, res_norm

        if (i_access(esr%f)) then
          if (esr%step == 1) then
            write(x_unit(esr%f),'("step ",i0,":")') esr%step
          else
            write(x_unit(esr%f),'(/,"step ",i0,":")') esr%step
          end if
          write(x_unit(esr%f),'(t4,"maximum directions = ",i0)') max_dir
          write(x_unit(esr%f),'(t4,"actual directions  = ",i0)') num_dir
          write(x_unit(esr%f),'(t4,"residual tolerance = ",es8.2)') res_tol
          write(x_unit(esr%f),'(t4,"residual norm      = ",es8.2,/)') res_norm
        end if

      end subroutine

      end module

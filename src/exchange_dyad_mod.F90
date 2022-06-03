!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module exchange_dyad_mod
!doc$ module exchange_dyad_mod

!     One datatype is available here: type(exchange_dyad_obj).
!     exchange_dyad_mod encapsulates an iterative method for finding the Optimized Effective Potential.

      use kind_mod
      use mpi_mod
      use io_mod
      use tagio_mod
      use error_mod
      use diary_mod
      use arg_mod
      use ghost_mod
      use xc_type_mod
      use xc_mod
      use electrons_sc_mod
      use wavefunctions_es_mod
      use multibasis_mod
      use multivector_mod
      use kpoints_mod
      use dyad_mod
      use dyad_kpoint_mod

!cod$
      implicit none
      private

!doc$
      public :: exchange_dyad

!cod$
      interface exchange_dyad
        module procedure exchange_dyad_xc
      end interface

      contains

! public routines

      function exchange_dyad_xc(xc,el,e) result(xd)
!doc$ function exchange_dyad(xc,el,e) result(xd)
        type(xc_obj) :: xc
        type(electrons_sc_obj) :: el
        real(double) :: e
        type(dyad_obj) :: xd
!       requires: 
!       effects: Builds a dyad_obj encoding a compact representation of the exchange operator
!       modifies: e
!       errors: Passes errors.

!cod$

        character(line_len), parameter :: init = "zeros"
        complex(double), parameter :: cp1 = (+1.0_double,0.0_double)
        complex(double), parameter :: cmh = (-0.5_double,0.0_double)
        integer :: ik, nb, nk
        real(double) :: ek, kg_ek, sg_ek
        complex(double), dimension(:,:), allocatable :: e_mat
        type(kpoints_obj) :: kp
        type(multibasis_obj) :: mb
        type(multivector_obj) :: e_mv, t_mv

        call my(xc)
        call my(el)

        call my(x_kpoints(el),kp)
        nk = x_n_kpoints(kp)
        nb = x_n_bands(el)

        call my(dyad(nk),xd)

        allocate(e_mat(nb,nb))

        kg_ek = 0.0_double

        select case(x_functional_dependence(xc))
        case (FD_ORBITAL,FD_HYBRID)
          if (mpi_nkgroups() > 1) call distribute(el)
        end select

        do ik = 1,nk

          if (mpi_mykgroup() /= x_kgroup_index(el,ik)) cycle

          call my(x_multivector(x_wavefunctions(el,ik)),e_mv)
          call my(x_multibasis(e_mv),mb)
          call my(multivector(mb,init),t_mv)

          select case(x_functional_dependence(xc))
          case (FD_ORBITAL)
            call xc_energy_and_derivative(xc,el,ik,kg_ek,t_mv)
          case(FD_HYBRID)
            call xc_energy_and_derivative(xc,el,ik,kg_ek,t_mv)
            call portion(x_hybrid_mixing(xc),t_mv)
          end select

          call overlap(e_mv,t_mv,e_mat)

!         Consistency Check
!         if (i_access(diaryfile())) write(x_unit(diaryfile()),'("check sum 1 =",f15.10,f15.10)') sum(e_mat)

          call transform(cp1,t_mv,cmh,e_mv,e_mat)

          call update(xd,ik,dyad_kpoint(t_mv,e_mv))

!         Consistency Check
!         t_mv = apply(x_dyad_kpoint(xd,ik),e_mv)
!         call overlap(e_mv,t_mv,e_mat)
!         if (i_access(diaryfile())) write(x_unit(diaryfile()),'("check sum 2 =",f15.10,f15.10)') sum(e_mat)

          call glean(thy(e_mv))
          call glean(thy(mb))
          call glean(thy(t_mv))

        end do

        select case(x_functional_dependence(xc))
        case (FD_ORBITAL,FD_HYBRID)
          if (mpi_nkgroups() > 1) call release(el)
        end select

        call xcomm_allreduce(XKGROUP,MPI_SUM,kg_ek,sg_ek)
        call xcomm_allreduce(XSGROUP,MPI_SUM,sg_ek,ek)

        select case(x_functional_dependence(xc))
        case (FD_ORBITAL)
           e = e + ek
        case (FD_HYBRID)
           e = e + x_hybrid_mixing(xc)*ek
        end select

100     if (allocated( e_mat )) deallocate( e_mat )
        call glean(thy(kp))
        call glean(thy(el))
        call glean(thy(xc))
        call bequeath(thy(xd))

        if (error("Exit exchange_dyad_mod::exchange_dyad_xc")) continue

      end function

      end module

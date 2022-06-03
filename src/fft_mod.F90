!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module fft_mod
!doc$ module fft_mod

      use kind_mod
      use error_mod
      use diary_mod
      use timing_mod
      use point_blas_mod
      use omp_lib
      use mpi_mod
      use fft3d_wrap
      use,intrinsic :: iso_c_binding

!     This module encapsulates fast-Fourier transform routines that utilize fftw3
!     interfaces for transforms on serial data and fftw2 interfaces for transforms
!     on distributed data.
!
!     FFT directions:  -1: r -> q  Forward transform
!                      +1: q -> r  Backward transform

!cod$
      implicit none
      private

      include 'fftw3.f'

!doc$
      ! non-indexed serial FFT plan
      type, public :: fft_serial_plan
        private
        integer(sizeof_double) :: fplan, bplan
        real(double) :: norm
        complex(double), dimension(:,:,:), pointer :: data
      end type

      ! indexed (pruned) serial FFT plan
      type, public :: fft_serial_plan_i
        private
        integer(sizeof_double), dimension(:), pointer :: fplan, bplan
        real(double) :: norm
        integer, dimension(:,:), pointer :: fstart, bstart
        complex(double), dimension(:,:,:), pointer :: data
      end type

      ! distributed FFT plan
      type, public :: fft_distributed_plan
        private
        real(double) :: plan
        real(double) :: norm
        type(c_ptr) :: pl
      end type

!cod$
      integer, parameter :: R_TO_Q = FFTW_FORWARD
      integer, parameter :: Q_TO_R = FFTW_BACKWARD

      integer :: fft_nthreads

!doc$
      public :: R_TO_Q
      public :: Q_TO_R

      public :: fft_start

      public :: assignment(=)
      public :: datalink

      public :: get_nf

      public :: fft_create_serial_plan
      public :: fft_destroy_serial_plan
      public :: fft_serial

      public :: fft_create_distributed_plan
      public :: fft_destroy_distributed_plan
      public :: fft_distributed

!cod$
      interface assignment(=)
        module procedure assign_fft_ni, assign_fft_i
      end interface
      interface datalink
        module procedure datalink_fft_ni, datalink_fft_i
      end interface

      interface fft_create_serial_plan
        module procedure fft_create_serial_plan_ni, fft_create_serial_plan_i
      end interface
      interface fft_destroy_serial_plan
        module procedure fft_destroy_serial_plan_ni, fft_destroy_serial_plan_i
      end interface
      interface fft_serial
        module procedure fft_serial_ni, fft_serial_i
      end interface

      contains

      subroutine fft_start()
!doc$ subroutine fft_start()
        integer :: fft_rthreads

        fft_nthreads = omp_get_max_threads()

        call dfftw_init_threads(fft_rthreads)
        if (error(fft_rthreads == 0,"FATAL ERROR: fft_rthreads = 0")) goto 100

        call dfftw_plan_with_nthreads(fft_nthreads)

100     if (error("Exit fft_mod::fft_start")) continue

      end subroutine

      subroutine assign_fft_ni(plan,plan2)
!doc$ subroutine assign(plan,plan2)
        type(fft_serial_plan), intent(inout) :: plan
        type(fft_serial_plan), intent(in) :: plan2

!cod$
        plan%fplan = plan2%fplan
        plan%bplan = plan2%bplan
        plan%norm = plan2%norm
        plan%data => plan2%data
      end subroutine

      subroutine assign_fft_i(plan,plan2)
!doc$ subroutine assign(plan,plan2)
        type(fft_serial_plan_i), intent(inout) :: plan
        type(fft_serial_plan_i), intent(in) :: plan2

!cod$
        integer :: np
        np = size(plan2%fplan)
        allocate( plan%fplan(np) )
        allocate( plan%bplan(np) )
        plan%fplan = plan2%fplan
        plan%bplan = plan2%bplan
        plan%norm = plan2%norm
        allocate( plan%fstart(3,np) )
        allocate( plan%bstart(3,np) )
        plan%fstart = plan2%fstart
        plan%bstart = plan2%bstart
        plan%data => plan2%data
      end subroutine

      function datalink_fft_ni(plan) result(data)
!doc$ function datalink(plan) result(data)
        type(fft_serial_plan) :: plan
        complex(double), dimension(:,:,:), pointer :: data
!       advice: On entry, data should be nullified
!       effects: Points data at plan%data

!cod$
        data => plan%data
      end function

      function datalink_fft_i(plan) result(data)
!doc$ function datalink(plan) result(data)
        type(fft_serial_plan_i) :: plan
        complex(double), dimension(:,:,:), pointer :: data
!       advice: On entry, data should be nullified
!       effects: Points data at plan%data

!cod$
        data => plan%data
      end function

      subroutine get_nf(nf,n2,n3,n5,n7)
!doc$ subroutine get_nf(nf,n2,n3,n5,n7)
        integer, intent (inout) :: nf
        integer, intent (out) :: n2, n3, n5, n7
!       requires: nf > 0
!       modifies: nf, n2, n3, n5, n7
!       effects : Overwrites nf with the smallest integer >= nf such that
!                            nf = (2**n2)*(3**n3)*(5**n5)*(7**n7).

!cod$
        integer :: nfd
        do
          nfd = nf
          n2 = 0
          do while ( mod(nfd,2) == 0 )
            nfd = nfd/2
            n2 = n2 + 1
          end do
          n3 = 0
          do while ( mod(nfd,3) == 0 )
            nfd = nfd/3
            n3 = n3 + 1
          end do
          n5 = 0
          do while ( mod(nfd,5) == 0 )
            nfd = nfd/5
            n5 = n5 + 1
          end do
          n7 = 0
          do while ( mod(nfd,7) == 0 )
            nfd = nfd/7
            n7 = n7 + 1
          end do
          if (nfd == 1) exit
          nf = nf + 1
        end do
      end subroutine

      subroutine fft_create_serial_plan_ni(dims,plan)
!doc$ subroutine fft_create_serial_plan(dims,plan)
        integer, dimension(3), intent(in) :: dims
        type(fft_serial_plan) :: plan

!cod$
        integer :: rank
        integer, dimension(3) :: n

        allocate( plan%data(dims(1),dims(2),dims(3)) )

        plan%norm = product(dims)

        rank = 3
        n = dims

        ! plan for a non-indexed R_TO_Q (forward) FFT
        call dfftw_plan_dft(plan%fplan,rank,n,plan%data,plan%data,R_TO_Q,FFTW_ESTIMATE)

        ! plan for a non-indexed Q_TO_R (backward) FFT
        call dfftw_plan_dft(plan%bplan,rank,n,plan%data,plan%data,Q_TO_R,FFTW_ESTIMATE)

100     if (error("Exit fft_mod::fft_create_serial_plan_ni")) continue

      end subroutine

      subroutine fft_create_serial_plan_i(dims,index,plan)
!doc$ subroutine fft_create_serial_plan(dims,index,plan)
        integer, dimension(3), intent(in) :: dims
        integer, dimension(:,:), intent(in) :: index
        type(fft_serial_plan_i) :: plan
        
!cod$
        integer :: i1, i2, i3, ip, np
        integer :: rank, howmany, stride, distance
        integer, dimension(1) :: n, nembed

        allocate( plan%data(dims(1),dims(2),dims(3)) )

        plan%norm = product(dims)

        ! determine the number of plans
        np = 0
        do i2 = 1,dims(2)
          np = np + 1
        end do
        do i3 = 1,dims(3)
          select case (index(1,i3))
          case (1)
            np = np + 2
          case (2)
            np = np + 3
          end select
        end do

        ! allocate space for the plans
        allocate( plan%fplan(np) )
        allocate( plan%bplan(np) )

        ! allocate space for the starting indices
        allocate( plan%fstart(3,np) )
        allocate( plan%bstart(3,np) )

        rank = 1

        ! plans for an indexed R_TO_Q (forward) FFT

        i1 = 1
        ip = 0

        do i2 = 1,dims(2)
          ip = ip + 1
          n = dims(3)
          howmany = dims(1)
          i3 = 1
          plan%fstart(:,ip) = (/i1,i2,i3/)
          nembed = n
          stride = dims(1)*dims(2)
          distance = 1
          call dfftw_plan_many_dft(plan%fplan(ip),rank,n,howmany, &
                                                & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                & R_TO_Q,FFTW_ESTIMATE)
        end do
        do i3 = 1,dims(3)
          select case (index(1,i3))
          case (1)
            ip = ip + 1
            n = dims(2)
            howmany = dims(1)
            i2 = 1
            plan%fstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = dims(1)
            distance = 1
            call dfftw_plan_many_dft(plan%fplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & R_TO_Q,FFTW_ESTIMATE)
            ip = ip + 1
            n = dims(1)
            howmany = index(3,i3)
            i2 = index(2,i3)
            plan%fstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = 1
            distance = dims(1)
            call dfftw_plan_many_dft(plan%fplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & R_TO_Q,FFTW_ESTIMATE)
          case (2)
            ip = ip + 1
            n = dims(2)
            howmany = dims(1)
            i2 = 1
            plan%fstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = dims(1)
            distance = 1
            call dfftw_plan_many_dft(plan%fplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & R_TO_Q,FFTW_ESTIMATE)
            ip = ip + 1
            n = dims(1)
            howmany = index(3,i3)
            i2 = index(2,i3)
            plan%fstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = 1
            distance = dims(1)
            call dfftw_plan_many_dft(plan%fplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & R_TO_Q,FFTW_ESTIMATE)
            ip = ip + 1
            n = dims(1)
            howmany = index(5,i3)
            i2 = index(4,i3)
            plan%fstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = 1
            distance = dims(1)
            call dfftw_plan_many_dft(plan%fplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & R_TO_Q,FFTW_ESTIMATE)
          end select
        end do

        ! create plans for an indexed Q_TO_R (backward) FFT

        i1 = 1
        ip = 0

        do i3 = 1,dims(3)
          select case (index(1,i3))
          case (1)
            ip = ip + 1
            n = dims(1)
            howmany = index(3,i3)
            i2 = index(2,i3)
            plan%bstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = 1
            distance = dims(1)
            call dfftw_plan_many_dft(plan%bplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & Q_TO_R,FFTW_ESTIMATE)
            ip = ip + 1
            n = dims(2)
            howmany = dims(1)
            i2 = 1
            plan%bstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = dims(1)
            distance = 1
            call dfftw_plan_many_dft(plan%bplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & Q_TO_R,FFTW_ESTIMATE)
          case (2)
            ip = ip + 1
            n = dims(1)
            howmany = index(3,i3)
            i2 = index(2,i3)
            plan%bstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = 1
            distance = dims(1)
            call dfftw_plan_many_dft(plan%bplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & Q_TO_R,FFTW_ESTIMATE)
            ip = ip + 1
            n = dims(1)
            howmany = index(5,i3)
            i2 = index(4,i3)
            plan%bstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = 1
            distance = dims(1)
            call dfftw_plan_many_dft(plan%bplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & Q_TO_R,FFTW_ESTIMATE)
            ip = ip + 1
            n = dims(2)
            howmany = dims(1)
            i2 = 1
            plan%bstart(:,ip) = (/i1,i2,i3/)
            nembed = n
            stride = dims(1)
            distance = 1
            call dfftw_plan_many_dft(plan%bplan(ip),rank,n,howmany, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                  & Q_TO_R,FFTW_ESTIMATE)
          end select
        end do
        do i2 = 1,dims(2)
          ip = ip + 1
          n = dims(3)
          howmany = dims(1)
          i3 = 1
          plan%bstart(:,ip) = (/i1,i2,i3/)
          nembed = n
          stride = dims(1)*dims(2)
          distance = 1
          call dfftw_plan_many_dft(plan%bplan(ip),rank,n,howmany, &
                                                & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                & plan%data(i1,i2,i3),nembed,stride,distance, &
                                                & Q_TO_R,FFTW_ESTIMATE)
        end do

100     if (error("Exit fft_mod::fft_create_serial_plan_i")) continue

      end subroutine

      subroutine fft_destroy_serial_plan_ni(plan)
!doc$ subroutine fft_destroy_serial_plan(plan)
        type(fft_serial_plan) :: plan

!cod$
         deallocate( plan%data )

         call dfftw_destroy_plan(plan%fplan)
         call dfftw_destroy_plan(plan%bplan)

      end subroutine

      subroutine fft_destroy_serial_plan_i(plan)
!doc$ subroutine fft_destroy_serial_plan(plan)
        type(fft_serial_plan_i) :: plan

!cod$
         integer :: ip

         deallocate( plan%data )

         do ip = 1,size(plan%fplan)
           call dfftw_destroy_plan(plan%fplan(ip))
           call dfftw_destroy_plan(plan%bplan(ip))
         end do
         deallocate( plan%fplan )
         deallocate( plan%bplan )

         deallocate( plan%fstart )
         deallocate( plan%bstart )

      end subroutine

      subroutine fft_serial_ni(dir,plan)
!doc$ subroutine fft_serial(dir,plan)
        integer, intent(in) :: dir
        type(fft_serial_plan) :: plan
!       requires: plan be initialized and dir = -1 (R_TO_Q) or +1 (Q_TO_R).
!       modifies: data
!       effects : Performs a serial FFT using dir as the sign of the exponent.

!cod$
        real(double) :: scale

        call start_timer("fft: serial")

        select case (dir)
        case (R_TO_Q)
          call dfftw_execute_dft(plan%fplan,plan%data,plan%data)
          scale = 1.0_double/plan%norm
          call point_mxs(plan%data,scale)
        case (Q_TO_R)
          call dfftw_execute_dft(plan%bplan,plan%data,plan%data)
        end select

        if (.not.error()) call stop_timer("fft: serial")

      end subroutine

      subroutine fft_serial_i(dir,plan)
!doc$ subroutine fft_serial(dir,plan)
        integer, intent(in) :: dir
        type(fft_serial_plan_i) :: plan
!       requires: Shape(data) > 0, dir = -1 (R_TO_Q) or +1 (Q_TO_R), plan be initialized.
!       modifies: plan%data
!       effects : Performs (pruned) orthogonal FFT with dir used as the sign of the exponent.

!cod$
        integer :: i1, i2, i3, ip
        real(double) :: scale

        call start_timer("fft: serial")

        select case (dir)
        case (R_TO_Q)
          do ip = 1,size(plan%fplan)
            i1 = plan%fstart(1,ip)
            i2 = plan%fstart(2,ip)
            i3 = plan%fstart(3,ip)
            call dfftw_execute_dft(plan%fplan(ip),plan%data(i1,i2,i3),plan%data(i1,i2,i3))
          end do
          scale = 1.0_double/plan%norm
          call point_mxs(plan%data,scale)
        case (Q_TO_R)
          do ip = 1,size(plan%bplan)
            i1 = plan%bstart(1,ip)
            i2 = plan%bstart(2,ip)
            i3 = plan%bstart(3,ip)
            call dfftw_execute_dft(plan%bplan(ip),plan%data(i1,i2,i3),plan%data(i1,i2,i3))
          end do
        end select

        if (.not.error()) call stop_timer("fft: serial")

      end subroutine

      subroutine fft_create_distributed_plan(comm,dims,locdims,base,plan)
!doc$ subroutine fft_create_distributed_plan(comm,dims,locdims,base,plan)
        integer :: comm
        integer, dimension(3) :: dims, base, locdims
        type(fft_distributed_plan) :: plan
        
!cod$
        integer :: howmany

        howmany = 1
        call create_nfplan_i(comm,dims,base,base+locdims-1,howmany,plan)
        plan%norm = product(dims)

      end subroutine

      subroutine fft_destroy_distributed_plan(plan)
!doc$ subroutine fft_destroy_distributed_plan(plan)
        type(fft_distributed_plan) :: plan
        
!cod$
        !!call fft3d_destroy(plan%pl)
        call fft_3d_destroy_plan(plan%plan)

      end subroutine

      subroutine fft_distributed(data,dir,plan)
!doc$ subroutine fft_distribured(data,dir,plan)
        complex(double), dimension(:,:,:), intent(inout), target :: data
        integer, intent(in) :: dir
        type(fft_distributed_plan) :: plan
!       requires: dir = -1 or +1 and plan be compatible with the data.
!       modifies: data
!       effects : Performs a distributed FFT using dir as the sign of the exponent.

!cod$
        real(double) :: scale

        call start_timer("fft: distributed")

        select case (dir)
        case (R_TO_Q)
      !!    call fft3d_compute(plan%pl,c_loc(data),c_loc(data),Q_TO_R)
          call fft_3d(data(1,1,1),data(1,1,1),dir,plan%plan)
          scale = 1.0_double/plan%norm
          call point_mxs(data,scale)
        case (Q_TO_R)
      !!    call fft3d_compute(plan%pl,c_loc(data),c_loc(data),R_TO_Q)
          call fft_3d(data(1,1,1),data(1,1,1),dir,plan%plan)
        end select

        if (.not.error()) call stop_timer("fft: distributed")

      end subroutine

      subroutine create_nfplan_i(comm,nf,fnf,lnf,howmany,plan)
        integer :: comm
        integer, dimension(3) :: nf, fnf, lnf
        integer :: howmany
        type(fft_distributed_plan) :: plan

        integer, parameter :: precision = 2, permute = 0, scaled = 0
        integer :: nf1, nf2, nf3, fnf1, fnf2, fnf3, lnf1, lnf2, lnf3
        integer :: nloc, fftsize, sendsize, recvsize

        nf1 = nf(1)
        nf2 = nf(2)
        nf3 = nf(3)

        fnf1 = fnf(1)
        fnf2 = fnf(2)
        fnf3 = fnf(3)

        lnf1 = lnf(1)
        lnf2 = lnf(2)
        lnf3 = lnf(3)

        !call fft3d_create(comm,precision,plan%pl)
        !call fft3d_set(plan%pl,"scale",scaled)
        !call fft3d_set(plan%pl,"pack",0)
        !call fft3d_setup(plan%pl,nf1,nf2,nf3,fnf1,lnf1,fnf2,lnf2,fnf3,lnf3,&
        !                         fnf1,lnf1,fnf2,lnf2,fnf3,lnf3,permute,fftsize,sendsize,recvsize)

        call fft_3d_create_plan(comm,nf1,nf2,nf3,fnf1,lnf1,fnf2,lnf2,fnf3,lnf3,fnf1, &
                              & lnf1,fnf2,lnf2,fnf3,lnf3,howmany,scaled,permute,nloc,plan%plan)
      end subroutine

      end module

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module mpi_mod
!doc$ module mpi_mod

!     One data type is available here: type(mpi_obj)

!     This module provides simplified interfaces to MPI calls and a private data structure (the_mpi) with the mpi environment.
!
!     For runs having multiple configurations:
!                  - Each process belongs to at least two distinct communicators: the full communicator
!                      (world_comm) and a communicator for one of the configurations (config_comm).
!                  - Each configuration is currently required to have the same number of processes. That is,
!                      mod(number of world processes,number of configurations) must equal 0.
!                  - Each process is identified by a configuration id and a process rank in that configuration -
!                      ( configuration id , rank ). The processes with configuration id 3 are designated (3,*).
!                      The mapping to the world rank is world rank = (configuration id - 1) * np + configuration
!                      rank, where np is the number of processes per configuration.
!                  - Each process from a configuration also belongs to the cross-configuration communicator
!                      (xconfig_comm) with the corresponding processes from other configurations, (*,3) for example.

!     For runs having multiple sgroups (spins):
!                  - Each process belongs to at least two distinct communicators: the full communicator,
!                      a communicator for one of the configurations, and a communicator for one of the
!                      sgroups (sgroup_comm).
!                  - Each sgroup is currently required to have the same number of processes. That is,
!                      mod(number of configuration processes,number of sgroups) must equal 0.
!                  - Each process is identified by a sgroup id and a process rank in that sgroup -
!                      ( sgroup id , rank ). The processes with sgroup id 3 are designated (3,*).
!                  - Each process from an sgroup also belongs to the cross-sgroup communicator
!                      (xsgroup_comm) with the corresponding processes from other sgroups, (*,3) for example.

!
!     For runs having multiple kgroups (sets of k-points):
!                  - Each process belongs to at least two distinct communicators: the full communicator,
!                      a communicator for one of the configurations, a communicator for one of the sgroups,
!                      and a communicator for one of the kgroups (kgroup_comm).
!                  - Each kgroup is currently required to have the same number of processes. That is,
!                      mod(number of configuration processes,number of kgroups) must equal 0.
!                  - Each process is identified by a kgroup id and a process rank in that kgroup -
!                      ( kgroup id , rank ). The processes with kgroup id 3 are designated (3,*).
!                  - Each process from a kgroup also belongs to the cross-kgroup communicator
!                      (xkgroup_comm) with the corresponding processes from other kgroups, (*,3) for example.

      use error_mod
      use kind_mod
      use mpi
      use utils_mod
#if defined(_OPENMP)
      use omp_lib
#endif

!cod$
      implicit none ; private

      integer, parameter :: WORLD   = 1
      integer, parameter :: CONFIG  = 2
      integer, parameter :: XCONFIG = 3
      integer, parameter :: SGROUP  = 4
      integer, parameter :: XSGROUP = 5
      integer, parameter :: KGROUP  = 6
      integer, parameter :: XKGROUP = 7

      type :: mpi_obj
                                 ! WORLD COMMUNICATOR VARIABLES (*,*):
        integer :: world_comm       ! communicator
        logical :: world_first      ! rank 0 process
        integer :: world_nprocs     ! number of processes
        integer :: world_myproc     ! rank
                                 ! CONFIG COMMUNICATOR VARIABLES (configuration id,*):
        integer :: nconfigs         ! number of configurations
        integer :: myconfig         ! configuration id
        integer :: config_comm      ! communicator
        logical :: config_first     ! rank 0 process
        integer :: config_nprocs    ! number of processes
        integer :: config_myproc    ! rank
                                 ! XCONFIG COMMUNICATOR VARIABLES (*,rank):
        integer :: xconfig_comm     ! communicator
        logical :: xconfig_first    ! rank 0 process
        integer :: xconfig_nprocs   ! number of processes (same as nconfigs)
        integer :: xconfig_myproc   ! rank
                                 ! SGROUP COMMUNICATOR VARIABLES (sgroup id,*):
        integer :: nsgroups         ! number of sgroups
        integer :: mysgroup         ! sgroup id
        integer :: sgroup_comm      ! communicator
        logical :: sgroup_first     ! rank 0 process
        integer :: sgroup_nprocs    ! number of processes
        integer :: sgroup_myproc    ! rank
                                 ! XSGROUP COMMUNICATOR VARIABLES (*,rank):
        integer :: xsgroup_comm     ! communicator
        logical :: xsgroup_first    ! rank 0 process
        integer :: xsgroup_nprocs   ! number of processes (same as nsgroups)
        integer :: xsgroup_myproc   ! rank
                                 ! KGROUP COMMUNICATOR VARIABLES (kgroup id,*):
        integer :: nkgroups         ! number of kgroups
        integer :: mykgroup         ! kgroup id
        integer :: kgroup_comm      ! communicator
        logical :: kgroup_first     ! rank 0 process
        integer :: kgroup_nprocs    ! number of processes
        integer :: kgroup_myproc    ! rank
                                 ! XKGROUP COMMUNICATOR VARIABLES (*,rank):
        integer :: xkgroup_comm     ! communicator
        logical :: xkgroup_first    ! rank 0 process
        integer :: xkgroup_nprocs   ! number of processes (same as nkgroups)
        integer :: xkgroup_myproc   ! rank
                                 ! PROCESS VARIABLES
        integer :: error            ! error status
        integer :: omp_nthreads     ! number of OMP threads

      end type

!doc$
      public :: mpi__start
      public :: mpi_stop
      public :: mpi_config_split
      public :: mpi_sgroup_split
      public :: mpi_kgroup_split
      public :: mpi_comm
      public :: mpi_isroot
      public :: mpi_nprocs
      public :: mpi_myproc
      public :: mpi_nconfigs
      public :: mpi_myconfig
      public :: mpi_nsgroups
      public :: mpi_mysgroup
      public :: mpi_nkgroups
      public :: mpi_mykgroup
      public :: x_nthreads
      public :: barrier
      public :: vote
      public :: broadcast_seh
      public :: broadcast
      public :: reduce
      public :: allreduce
      public :: allgather
      public :: allgatherv
      public :: alltoallv
      public :: blocking_send
      public :: nonblocking_send
      public :: blocking_recv
      public :: xcomm_broadcast
      public :: xcomm_reduce
      public :: xcomm_allreduce
      public :: xcomm_allgather
      public :: xcomm_rank_allreduce
      public :: xcomm_pair_allreduce
      public :: MPI_SUM
      public :: MPI_PROD
      public :: MPI_MAX
      public :: MPI_MIN
      public :: MPI_LAND
      public :: MPI_LOR
      public :: MPI_LXOR
      public :: WORLD
      public :: CONFIG
      public :: XCONFIG
      public :: SGROUP
      public :: XSGROUP
      public :: KGROUP
      public :: XKGROUP

!cod$
      type(mpi_obj) :: the_mpi

      interface broadcast_seh
        module procedure broadcast_ch_seh, &
                         broadcast_log0_seh, &
                         broadcast_int0_seh
      end interface
      interface broadcast
        module procedure broadcast_ch, &
                         broadcast_log0, broadcast_log1, broadcast_log2, broadcast_log3,  &
                         broadcast_int0, broadcast_int1, broadcast_int2, broadcast_int3, &
                         broadcast_dpr0, broadcast_dpr1, broadcast_dpr2, broadcast_dpr3, broadcast_dpr4, &
                         broadcast_dpc0, broadcast_dpc1, broadcast_dpc2, broadcast_dpc3, broadcast_dpc4
      end interface
      interface reduce
        module procedure reduce_log0, &
                         reduce_int0, reduce_int1, reduce_int2, reduce_int3, &
                         reduce_dpr0, reduce_dpr1, reduce_dpr2, reduce_dpr3, &
                         reduce_dpc0, reduce_dpc1, reduce_dpc2, reduce_dpc3
      end interface
      interface allreduce
        module procedure allreduce_log0, allreduce_log1, &
                         allreduce_int0, allreduce_int1, allreduce_int2, allreduce_int3, &
                         allreduce_dpr0, allreduce_dpr1, allreduce_dpr2, allreduce_dpr3, &
                         allreduce_dpc0, allreduce_dpc1, allreduce_dpc2, allreduce_dpc3, allreduce_dpc4, allreduce_dpc2_3
      end interface
      interface allgather
        module procedure allgather_int0
      end interface
      interface allgatherv
        module procedure allgatherv_int1, allgatherv_int2, &
                         allgatherv_dpr1, allgatherv_dpr2, allgatherv_dpr3, &
                         allgatherv_dpc1
      end interface
      interface alltoallv
        module procedure alltoallv_dpc1_2, alltoallv_dpc2_1
      end interface
      interface blocking_send
        module procedure blocking_send_int0, blocking_send_dpc1
      end interface
      interface nonblocking_send
        module procedure nonblocking_send_int0, nonblocking_send_dpc1
      end interface
      interface blocking_recv
        module procedure blocking_recv_int0, blocking_recv_dpc1
      end interface
      interface xcomm_broadcast
        module procedure xcomm_broadcast_dpc2
      end interface
      interface xcomm_reduce
        module procedure xcomm_reduce_int0, xcomm_reduce_int1, &
                         xcomm_reduce_dpr1,                    xcomm_reduce_dpr3, xcomm_reduce_dpr4, &
                         xcomm_reduce_dpc1, xcomm_reduce_dpc2, xcomm_reduce_dpc3, xcomm_reduce_dpc4
      end interface
      interface xcomm_allreduce
        module procedure xcomm_allreduce_log0, &
                         xcomm_allreduce_int0, xcomm_allreduce_int1, &
                         xcomm_allreduce_dpr0, xcomm_allreduce_dpr1, xcomm_allreduce_dpr2, xcomm_allreduce_dpr3, &
                                                                                           xcomm_allreduce_dpr4, &
                         xcomm_allreduce_dpc1,                                             xcomm_allreduce_dpc3
      end interface
      interface xcomm_allgather
        module procedure xcomm_allgather_int0, &
                         xcomm_allgather_dpr0, xcomm_allgather_dpr1, xcomm_allgather_dpr2
      end interface
      interface xcomm_rank_allreduce
        module procedure xcomm_rank_allreduce_dpr1, xcomm_rank_allreduce_dpr2, xcomm_rank_allreduce_dpr3, &
                                                                               xcomm_rank_allreduce_dpr4, &
                                                                               xcomm_rank_allreduce_dpc3
      end interface
      interface xcomm_pair_allreduce
        module procedure xcomm_pair_allreduce_dpr0, xcomm_pair_allreduce_dpr1, xcomm_pair_allreduce_dpr2
      end interface

      contains

      subroutine mpi__start()
!doc$ subroutine mpi__start()
!        effects: Starts the MPI system.
!        errors: Problems starting MPI.
!        notes: Error conditions trigger a program halt.

!cod$
         integer :: provided
         character(line_len) :: val

         ! Initialize the WORLD communicator

         the_mpi%error = 0

#if defined(_OPENMP)
         call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,the_mpi%error)
#else
         call MPI_Init(the_mpi%error)
#endif
         if ( the_mpi%error /= 0 ) then
            write(6,'("ERROR: There was a problem with MPI_Init()")')
            stop
         end if

         the_mpi%world_comm = MPI_COMM_WORLD

         call MPI_Comm_rank(the_mpi%world_comm,the_mpi%world_myproc,the_mpi%error)
         if ( the_mpi%error /= 0 ) then
            write(6,'("ERROR: There was a problem with MPI_Comm_rank()")')
            call mpi_stop() ; stop
         end if

         the_mpi%world_first = ( the_mpi%world_myproc == 0 )

         call MPI_Comm_size(the_mpi%world_comm,the_mpi%world_nprocs,the_mpi%error)
         if ( the_mpi%error /= 0 ) then
            if (the_mpi%world_first) write(6,'("ERROR: There was a problem with MPI_Comm_size()")')
            call mpi_stop() ; stop
         end if

         ! Initialize the CONFIG and XCONFIG variables

         the_mpi%nconfigs = 1
         the_mpi%myconfig = 1

         the_mpi%config_comm   = the_mpi%world_comm
         the_mpi%config_first  = the_mpi%world_first
         the_mpi%config_nprocs = the_mpi%world_nprocs
         the_mpi%config_myproc = the_mpi%world_myproc

         the_mpi%xconfig_comm   = the_mpi%world_comm
         the_mpi%xconfig_first  = the_mpi%world_first
         the_mpi%xconfig_nprocs = 1
         the_mpi%xconfig_myproc = 0

         ! Initialize the SGROUP and XSGROUP variables

         the_mpi%nsgroups = 1
         the_mpi%mysgroup = 1

         the_mpi%sgroup_comm   = the_mpi%config_comm
         the_mpi%sgroup_first  = the_mpi%config_first
         the_mpi%sgroup_nprocs = the_mpi%config_nprocs
         the_mpi%sgroup_myproc = the_mpi%config_myproc

         the_mpi%xsgroup_comm   = the_mpi%config_comm
         the_mpi%xsgroup_first  = the_mpi%config_first
         the_mpi%xsgroup_nprocs = 1
         the_mpi%xsgroup_myproc = 0

         ! Initialize the KGROUP and XKGROUP variables

         the_mpi%nkgroups = 1
         the_mpi%mykgroup = 1

         the_mpi%kgroup_comm   = the_mpi%sgroup_comm
         the_mpi%kgroup_first  = the_mpi%sgroup_first
         the_mpi%kgroup_nprocs = the_mpi%sgroup_nprocs
         the_mpi%kgroup_myproc = the_mpi%sgroup_myproc

         the_mpi%xkgroup_comm   = the_mpi%sgroup_comm
         the_mpi%xkgroup_first  = the_mpi%sgroup_first
         the_mpi%xkgroup_nprocs = 1
         the_mpi%xkgroup_myproc = 0

         ! Initialize the use of OpenMP threads

#if defined(_OPENMP)
         call get_environment_variable("OMP_NUM_THREADS",value=val)
         select case ( trimstr( val ) )
         case ( "" )
            the_mpi%omp_nthreads = 1
         case default
            the_mpi%omp_nthreads = omp_get_max_threads()
         end select
         call MPI_Bcast(the_mpi%omp_nthreads,1,MPI_INTEGER,0,the_mpi%world_comm,the_mpi%error)
         call omp_set_num_threads(the_mpi%omp_nthreads)
#else
         the_mpi%omp_nthreads = 0
#endif

      end subroutine

      subroutine mpi_stop()
!doc$ subroutine mpi_stop()
!        effects: Stops the MPI system.

!cod$
         call MPI_Finalize(the_mpi%error)
      end subroutine

      subroutine mpi_config_split(nc)
!doc$ subroutine mpi_config_split(nc)
        integer, intent(in) :: nc
!       requires: nc > 1. mod(the_mpi%world_nprocs,nc) == 0.
!       modifies: the_mpi
!       effects: Sets up the_mpi to allow simultaneous calculations for multiple configurations.
!       errors: Problems spawning the new communicators.
!       notes: Error conditions trigger a program halt.
!              This routine should be called only from sys_mod.

!cod$
        integer :: ip, color, key, myc, nppc, config_error, world_error

        the_mpi%nconfigs = nc

        nppc = (the_mpi%world_nprocs + the_mpi%nconfigs - 1)/the_mpi%nconfigs
        myc = the_mpi%world_myproc/nppc + 1
        the_mpi%myconfig = myc

        call MPI_COMM_SPLIT(the_mpi%world_comm,myc,0,the_mpi%config_comm,world_error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SPLIT in mpi_config_split: state = ",i0)') world_error
          call mpi_stop() ; stop
        end if

        call MPI_COMM_RANK(the_mpi%config_comm,the_mpi%config_myproc,config_error)
        call MPI_ALLREDUCE(config_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_RANK in mpi_config_split: state = ",i0)') config_error
          call mpi_stop() ; stop
        end if
        the_mpi%config_first = ( the_mpi%config_myproc == 0 )

        call MPI_COMM_SIZE(the_mpi%config_comm,the_mpi%config_nprocs,config_error)
        call MPI_ALLREDUCE(config_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SIZE in mpi_config_split: state = ",i0)') config_error
          call mpi_stop() ; stop
        end if

        do ip = 0,(the_mpi%config_nprocs - 1)
          if (the_mpi%config_myproc == ip) then
            color = the_mpi%config_myproc
            key = the_mpi%myconfig
            call MPI_COMM_SPLIT(the_mpi%world_comm,color,key,the_mpi%xconfig_comm,config_error)
          end if
        end do
        call MPI_ALLREDUCE(config_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SPLIT for xconfig: state = ",i0)') config_error
          call mpi_stop() ; stop
        end if

        do ip = 0,(the_mpi%config_nprocs - 1)
          if (the_mpi%config_myproc == ip) then
            call MPI_COMM_RANK(the_mpi%xconfig_comm,the_mpi%xconfig_myproc,config_error)
          end if
        end do
        call MPI_ALLREDUCE(config_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_RANK for xconfig: state = ",i0)') config_error
          call mpi_stop() ; stop
        end if
        the_mpi%xconfig_first = ( the_mpi%xconfig_myproc == 0 )

        the_mpi%xconfig_nprocs = the_mpi%nconfigs

        the_mpi%nsgroups = 1
        the_mpi%mysgroup = 1

        the_mpi%sgroup_comm   = the_mpi%config_comm
        the_mpi%sgroup_first  = the_mpi%config_first
        the_mpi%sgroup_nprocs = the_mpi%config_nprocs
        the_mpi%sgroup_myproc = the_mpi%config_myproc

        the_mpi%xsgroup_comm   = the_mpi%config_comm
        the_mpi%xsgroup_first  = the_mpi%config_first
        the_mpi%xsgroup_nprocs = 1
        the_mpi%xsgroup_myproc = 0

        the_mpi%nkgroups = 1
        the_mpi%mykgroup = 1

        the_mpi%kgroup_comm   = the_mpi%sgroup_comm
        the_mpi%kgroup_first  = the_mpi%sgroup_first
        the_mpi%kgroup_nprocs = the_mpi%sgroup_nprocs
        the_mpi%kgroup_myproc = the_mpi%sgroup_myproc

        the_mpi%xkgroup_comm   = the_mpi%sgroup_comm
        the_mpi%xkgroup_first  = the_mpi%sgroup_first
        the_mpi%xkgroup_nprocs = 1
        the_mpi%xkgroup_myproc = 0

      end subroutine

      subroutine mpi_sgroup_split(ng)
!doc$ subroutine mpi_sgroup_split(ng)
        integer, intent(in) :: ng
!       requires: ng = 2. mod(the_mpi%config_nprocs,ng) == 0.
!       modifies: the_mpi
!       effects: Sets up the_mpi to allow calculations for two sgroups.
!       errors: Problems spawning the new communicators.
!       notes: Error conditions trigger a program halt.

!cod$
        integer :: ip, color, key, myg, nppg, config_error, sgroup_error, world_error

        the_mpi%nsgroups = ng

        nppg = the_mpi%config_nprocs/the_mpi%nsgroups
        myg = the_mpi%config_myproc/nppg + 1
        the_mpi%mysgroup = myg

        call MPI_COMM_SPLIT(the_mpi%config_comm,myg,0,the_mpi%sgroup_comm,config_error)
        call MPI_ALLREDUCE(config_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SPLIT in mpi_sgroup_split: state = ",i0)') config_error
          call mpi_stop() ; stop
        end if

        call MPI_COMM_RANK(the_mpi%sgroup_comm,the_mpi%sgroup_myproc,sgroup_error)
        call MPI_ALLREDUCE(sgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_RANK in mpi_sgroup_split: state = ",i0)') sgroup_error
          call mpi_stop() ; stop
        end if
        the_mpi%sgroup_first = ( the_mpi%sgroup_myproc == 0 )

        call MPI_COMM_SIZE(the_mpi%sgroup_comm,the_mpi%sgroup_nprocs,sgroup_error)
        call MPI_ALLREDUCE(sgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SIZE in mpi_sgroup_split: state = ",i0)') sgroup_error
          call mpi_stop() ; stop
        end if

        do ip = 0,(the_mpi%sgroup_nprocs - 1)
          if (the_mpi%sgroup_myproc == ip) then
            color = the_mpi%sgroup_myproc
            key = the_mpi%mysgroup
            call MPI_COMM_SPLIT(the_mpi%config_comm,color,key,the_mpi%xsgroup_comm,sgroup_error)
          end if
        end do
        call MPI_ALLREDUCE(sgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SPLIT for xsgroup: state = ",i0)') sgroup_error
          call mpi_stop() ; stop
        end if

        do ip = 0,(the_mpi%sgroup_nprocs - 1)
          if (the_mpi%sgroup_myproc == ip) then
            call MPI_COMM_RANK(the_mpi%xsgroup_comm,the_mpi%xsgroup_myproc,sgroup_error)
          end if
        end do
        call MPI_ALLREDUCE(sgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_RANK for xsgroup: state = ",i0)') sgroup_error
          call mpi_stop() ; stop
        end if
        the_mpi%xsgroup_first = ( the_mpi%xsgroup_myproc == 0 )

        the_mpi%xsgroup_nprocs = the_mpi%nsgroups

        the_mpi%nkgroups = 1
        the_mpi%mykgroup = 1

        the_mpi%kgroup_comm   = the_mpi%sgroup_comm
        the_mpi%kgroup_first  = the_mpi%sgroup_first
        the_mpi%kgroup_nprocs = the_mpi%sgroup_nprocs
        the_mpi%kgroup_myproc = the_mpi%sgroup_myproc

        the_mpi%xkgroup_comm   = the_mpi%sgroup_comm
        the_mpi%xkgroup_first  = the_mpi%sgroup_first
        the_mpi%xkgroup_nprocs = 1
        the_mpi%xkgroup_myproc = 0

      end subroutine

      subroutine mpi_kgroup_split(ng)
!doc$ subroutine mpi_kgroup_split(ng)
        integer, intent(in) :: ng
!       requires: ng > 1. mod(the_mpi%sgroup_nprocs,ng) == 0.
!       modifies: the_mpi
!       effects: Sets up the_mpi to allow calculations for multiple kgroups.
!       errors: Problems spawning the new communicators.
!       notes: Error conditions trigger a program halt.

!cod$
        integer :: ip, color, key, myg, nppg, kgroup_error, sgroup_error, world_error

        the_mpi%nkgroups = ng

        nppg = (the_mpi%sgroup_nprocs + the_mpi%nkgroups - 1)/the_mpi%nkgroups
        myg = the_mpi%sgroup_myproc/nppg + 1
        the_mpi%mykgroup = myg

        call MPI_COMM_SPLIT(the_mpi%sgroup_comm,myg,0,the_mpi%kgroup_comm,sgroup_error)
        call MPI_ALLREDUCE(sgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SPLIT in mpi_kgroup_split: state = ",i0)') sgroup_error
          call mpi_stop() ; stop
        end if

        call MPI_COMM_RANK(the_mpi%kgroup_comm,the_mpi%kgroup_myproc,kgroup_error)
        call MPI_ALLREDUCE(kgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_RANK in mpi_kgroup_split: state = ",i0)') kgroup_error
          call mpi_stop() ; stop
        end if
        the_mpi%kgroup_first = ( the_mpi%kgroup_myproc == 0 )

        call MPI_COMM_SIZE(the_mpi%kgroup_comm,the_mpi%kgroup_nprocs,kgroup_error)
        call MPI_ALLREDUCE(kgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SIZE in mpi_kgroup_split: state = ",i0)') kgroup_error
          call mpi_stop() ; stop
        end if

        do ip = 0,(the_mpi%kgroup_nprocs - 1)
          if (the_mpi%kgroup_myproc == ip) then
            color = the_mpi%kgroup_myproc
            key = the_mpi%mykgroup
            call MPI_COMM_SPLIT(the_mpi%sgroup_comm,color,key,the_mpi%xkgroup_comm,kgroup_error)
          end if
        end do
        call MPI_ALLREDUCE(kgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_SPLIT for xkgroup: state = ",i0)') kgroup_error
          call mpi_stop() ; stop
        end if

        do ip = 0,(the_mpi%kgroup_nprocs - 1)
          if (the_mpi%kgroup_myproc == ip) then
            call MPI_COMM_RANK(the_mpi%xkgroup_comm,the_mpi%xkgroup_myproc,kgroup_error)
          end if
        end do
        call MPI_ALLREDUCE(kgroup_error,world_error,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        if (world_error /= 0) then
          write(6,'("ERROR: Problem with MPI_COMM_RANK for xkgroup: state = ",i0)') kgroup_error
          call mpi_stop() ; stop
        end if
        the_mpi%xkgroup_first = ( the_mpi%xkgroup_myproc == 0 )

        the_mpi%xkgroup_nprocs = the_mpi%nkgroups

      end subroutine

      function mpi_comm(sc) result(c)
!doc$ function mpi_comm(sc) result(c)
        integer, intent(in) :: sc
        integer :: c
!       effects: Returns the_mpi%sc_comm.

!cod$
        select case (sc)
        case (WORLD)
          c = the_mpi%world_comm
        case (CONFIG)
          c = the_mpi%config_comm
        case (XCONFIG)
          c = the_mpi%xconfig_comm
        case (SGROUP)
          c = the_mpi%sgroup_comm
        case (XSGROUP)
          c = the_mpi%xsgroup_comm
        case (KGROUP)
          c = the_mpi%kgroup_comm
        case (XKGROUP)
          c = the_mpi%xkgroup_comm
        end select
      end function 

      function mpi_isroot(sc) result(f)
!doc$ function mpi_isroot(sc) result(f)
        integer, intent(in) :: sc
        logical :: f
!       effects: Returns the_mpi%sc_first.

!cod$
        select case (sc)
        case (WORLD)
          f = the_mpi%world_first
        case (CONFIG)
          f = the_mpi%config_first
        case (XCONFIG)
          f = the_mpi%xconfig_first
        case (SGROUP)
          f = the_mpi%sgroup_first
        case (XSGROUP)
          f = the_mpi%xsgroup_first
        case (KGROUP)
          f = the_mpi%kgroup_first
        case (XKGROUP)
          f = the_mpi%xkgroup_first
        end select
      end function

      function mpi_nprocs(sc) result(n)
!doc$ function mpi_nprocs(sc) result(n)
        integer, intent(in) :: sc
        integer :: n
!       effects: Returns the_mpi%sc_nprocs.

!cod$
        select case (sc)
        case (WORLD)
          n = the_mpi%world_nprocs
        case (CONFIG)
          n = the_mpi%config_nprocs
        case (XCONFIG)
          n = the_mpi%xconfig_nprocs
        case (SGROUP)
          n = the_mpi%sgroup_nprocs
        case (XSGROUP)
          n = the_mpi%xsgroup_nprocs
        case (KGROUP)
          n = the_mpi%kgroup_nprocs
        case (XKGROUP)
          n = the_mpi%xkgroup_nprocs
        end select
      end function 

      function mpi_myproc(sc) result(m)
!doc$ function mpi_myproc(sc) result(m)
        integer, intent(in) :: sc
        integer :: m
!       effects: Returns the_mpi%sc_myproc.

!cod$
        select case (sc)
        case (WORLD)
          m = the_mpi%world_myproc
        case (CONFIG)
          m = the_mpi%config_myproc
        case (XCONFIG)
          m = the_mpi%xconfig_myproc
        case (SGROUP)
          m = the_mpi%sgroup_myproc
        case (XSGROUP)
          m = the_mpi%xsgroup_myproc
        case (KGROUP)
          m = the_mpi%kgroup_myproc
        case (XKGROUP)
          m = the_mpi%xkgroup_myproc
        end select
      end function 

      function mpi_nconfigs() result(n)
!doc$ function mpi_nconfigs() result(n)
        integer :: n
!       effects: Returns the_mpi%nconfigs

!cod$
        n = the_mpi%nconfigs
      end function 

      function mpi_myconfig() result(n)
!doc$ function mpi_myconfig() result(n)
        integer :: n
!       effects: Returns the_mpi%myconfig

!cod$
        n = the_mpi%myconfig
      end function 

      function mpi_nsgroups() result(n)
!doc$ function mpi_nsgroups() result(n)
        integer :: n
!       effects: Returns the_mpi%nsgroups

!cod$
        n = the_mpi%nsgroups
      end function 

      function mpi_mysgroup() result(n)
!doc$ function mpi_mysgroup() result(n)
        integer :: n
!       effects: Returns the_mpi%mysgroup.

!cod$
        n = the_mpi%mysgroup
      end function 

      function mpi_nkgroups() result(n)
!doc$ function mpi_nkgroups() result(n)
        integer :: n
!       effects: Returns the_mpi%nkgroups

!cod$
        n = the_mpi%nkgroups
      end function 

      function mpi_mykgroup() result(n)
!doc$ function mpi_mykgroup() result(n)
        integer :: n
!       effects: Returns the_mpi%mykgroup.

!cod$
        n = the_mpi%mykgroup
      end function 

      function x_nthreads() result( nt )
!doc$ function x_nthreads() result(n)
        integer :: nt
!       effects: Returns the_mpi%nthreads.

!cod$
         nt = the_mpi%omp_nthreads
      end function x_nthreads

      subroutine barrier(sc)
!doc$ subroutine barrier(sc)
        integer, intent(in) :: sc
!       effects: Creates an sc synchronization barrier.

!cod$
        select case (sc)
        case (WORLD)
          call MPI_BARRIER(the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BARRIER(the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BARRIER(the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BARRIER(the_mpi%kgroup_comm,the_mpi%error)
        end select
      end subroutine

      function vote(sc,l) result(v)
!doc$ function vote(sc,l) result(v)
        integer, intent(in) :: sc
        logical :: l
        integer :: v
!       effects: Returns the number of process in sc asserting l.

!cod$
        integer :: sig, tally
        sig = 0
        if (l) sig = 1
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(sig,tally,1,MPI_INTEGER,MPI_SUM,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(sig,tally,1,MPI_INTEGER,MPI_SUM,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(sig,tally,1,MPI_INTEGER,MPI_SUM,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(sig,tally,1,MPI_INTEGER,MPI_SUM,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
        v = tally
100     if (error("Exit mpi_mod::vote")) continue
      end function

!doc$ subroutine broadcast_seh(data)
!       character, logical, integer :: data
!          rank 1 for character
!          rank 0 for logical, integer
!       requires: data has the same dimension on all processes.
!       modifies: data
!       effects : data is overwritten with values from the first world process.
!       notes: Special error handling (seh) is done.

!cod$ 

      subroutine broadcast_ch_seh(data)
        character(*) :: data
        call MPI_BCAST(data,len(data),MPI_CHARACTER,0,the_mpi%world_comm,the_mpi%error)
        if (the_mpi%error /= 0) then
          write(6,'("ERROR: problem with mpi_mod::broadcast_ch_seh")')
          call mpi_stop() ; stop
        end if
      end subroutine

      subroutine broadcast_log0_seh(data)
        logical :: data
        call MPI_BCAST(data,1,MPI_LOGICAL,0,the_mpi%world_comm,the_mpi%error)
        if (the_mpi%error /= 0) then
          write(6,'("ERROR: problem with mpi_mod::broadcast_log0_seh")')
          call mpi_stop() ; stop
        end if
      end subroutine

      subroutine broadcast_int0_seh(data)
        integer :: data
        call MPI_BCAST(data,1,MPI_INTEGER,0,the_mpi%world_comm,the_mpi%error)
        if (the_mpi%error /= 0) then
          write(6,'("ERROR: problem with mpi_mod::broadcast_int0_seh")')
          call mpi_stop() ; stop
        end if
      end subroutine

!doc$ subroutine broadcast(sc,data)
!       integer, intent(in) (WORLD, CONFIG, SGROUP, or KGROUP) :: sc
!       character, logical, integer, real(double), or complex(double) :: data
!          rank 1 for character
!          rank 0,1,2,3 for logical, integer, real(double), and complex(double)
!       requires: data has the same dimension on all sc processes.
!       modifies: data
!       effects : data is overwritten with values from the first sc process.

!cod$ 

      subroutine broadcast_ch(sc,data)
        integer, intent(in) :: sc
        character(*) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data,len(data),MPI_CHARACTER,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data,len(data),MPI_CHARACTER,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data,len(data),MPI_CHARACTER,0,the_mpi%Sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data,len(data),MPI_CHARACTER,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_ch")) continue
      end subroutine

      subroutine broadcast_log0(sc,data)
        integer, intent(in) :: sc
        logical :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data,1,MPI_LOGICAL,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data,1,MPI_LOGICAL,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data,1,MPI_LOGICAL,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data,1,MPI_LOGICAL,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_log0")) continue
      end subroutine

      subroutine broadcast_log1(sc,data)
        integer, intent(in) :: sc
        logical, dimension(1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1),size(data),MPI_LOGICAL,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1),size(data),MPI_LOGICAL,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1),size(data),MPI_LOGICAL,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1),size(data),MPI_LOGICAL,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_log1")) continue
      end subroutine

      subroutine broadcast_log2(sc,data)
        integer, intent(in) :: sc
        logical, dimension(1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1),size(data),MPI_LOGICAL,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1),size(data),MPI_LOGICAL,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_LOGICAL,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_LOGICAL,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_log2")) continue
      end subroutine

      subroutine broadcast_log3(sc,data)
        integer, intent(in) :: sc
        logical, dimension(1:,1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1,1),size(data),MPI_LOGICAL,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1,1),size(data),MPI_LOGICAL,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_LOGICAL,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_LOGICAL,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_log3")) continue
      end subroutine

      subroutine broadcast_int0(sc,data)
        integer, intent(in) :: sc
        integer :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data,1,MPI_INTEGER,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data,1,MPI_INTEGER,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data,1,MPI_INTEGER,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data,1,MPI_INTEGER,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_int0")) continue
      end subroutine

      subroutine broadcast_int1(sc,data)
        integer, intent(in) :: sc
        integer, dimension(1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1),size(data),MPI_INTEGER,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1),size(data),MPI_INTEGER,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1),size(data),MPI_INTEGER,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1),size(data),MPI_INTEGER,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_int1")) continue
      end subroutine

      subroutine broadcast_int2(sc,data)
        integer, intent(in) :: sc
        integer, dimension(1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1),size(data),MPI_INTEGER,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1),size(data),MPI_INTEGER,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_INTEGER,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_INTEGER,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_int2")) continue
      end subroutine

      subroutine broadcast_int3(sc,data)
        integer, intent(in) :: sc
        integer, dimension(1:,1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1,1),size(data),MPI_INTEGER,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1,1),size(data),MPI_INTEGER,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_INTEGER,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_INTEGER,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_int3")) continue
      end subroutine

      subroutine broadcast_dpr0(sc,data)
        integer, intent(in) :: sc
        real(double) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data,1,MPI_DOUBLE_PRECISION,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpr0")) continue
      end subroutine

      subroutine broadcast_dpr1(sc,data)
        integer, intent(in) :: sc
        real(double), dimension(1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpr1")) continue
      end subroutine

      subroutine broadcast_dpr2(sc,data)
        integer, intent(in) :: sc
        real(double), dimension(1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpr2")) continue
      end subroutine

      subroutine broadcast_dpr3(sc,data)
        integer, intent(in) :: sc
        real(double), dimension(1:,1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpr3")) continue
      end subroutine

      subroutine broadcast_dpr4(sc,data)
        integer, intent(in) :: sc
        real(double), dimension(1:,1:,1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_PRECISION,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpr4")) continue
      end subroutine

      subroutine broadcast_dpc0(sc,data)
        integer, intent(in) :: sc
        complex(double) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data,1,MPI_DOUBLE_COMPLEX,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpc0")) continue
      end subroutine

      subroutine broadcast_dpc1(sc,data)
        integer, intent(in) :: sc
        complex(double), dimension(1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpc1")) continue
      end subroutine

      subroutine broadcast_dpc2(sc,data)
        integer, intent(in) :: sc
        complex(double), dimension(1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpc2")) continue
      end subroutine

      subroutine broadcast_dpc3(sc,data)
        integer, intent(in) :: sc
        complex(double), dimension(1:,1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpc3")) continue
      end subroutine

      subroutine broadcast_dpc4(sc,data)
        integer, intent(in) :: sc
        complex(double), dimension(1:,1:,1:,1:) :: data
        select case (sc)
        case (WORLD)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_BCAST(data(1,1,1,1),size(data),MPI_DOUBLE_COMPLEX,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::broadcast_dpc4")) continue
      end subroutine

!doc$ subroutine reduce(sc,op,data,tot)
!       integer, intent(in) (WORLD, CONFIG, SGROUP, or KGROUP) :: sc
!       integer, intent(in) (MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN, or MPI_LAND) :: op
!       logical, integer, real(double), or complex(double) :: data
!          rank 0 for logical
!          rank 0,1,2,3 for integer, real(double), and complex(double)
!       same type and rank as data :: tot
!       requires: data has the same rank and dimensions on all sc processes.
!       modifies: tot
!       effects: tot on the first sc proccess is overwritten with data values from all sc processes reduced by op.

!cod$
 
      subroutine reduce_log0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        logical :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data,tot,1,MPI_LOGICAL,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data,tot,1,MPI_LOGICAL,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data,tot,1,MPI_LOGICAL,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data,tot,1,MPI_LOGICAL,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_log0")) continue
      end subroutine

      subroutine reduce_int0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_int0")) continue
      end subroutine

      subroutine reduce_int1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_int1")) continue
      end subroutine

      subroutine reduce_int2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_int2")) continue
      end subroutine

      subroutine reduce_int3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_int3")) continue
      end subroutine

      subroutine reduce_dpr0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double) :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpr0")) continue
      end subroutine

      subroutine reduce_dpr1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpr1")) continue
      end subroutine

      subroutine reduce_dpr2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpr2")) continue
      end subroutine

      subroutine reduce_dpr3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpr3")) continue
      end subroutine

      subroutine reduce_dpc0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double) :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpc0")) continue
      end subroutine

      subroutine reduce_dpc1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpc1")) continue
      end subroutine

      subroutine reduce_dpc2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpc2")) continue
      end subroutine

      subroutine reduce_dpc3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::reduce_dpc3")) continue
      end subroutine

!doc$ subroutine allreduce(sc,op,data,tot,i)
!       integer, intent(in) (WORLD, CONFIG, SGROUP, or KGROUP) :: sc
!       integer, intent(in) (MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN, MPI_LAND, MPI_LOR, or MPI_LXOR) :: op
!       logical, integer, real(double), or complex(double) :: data
!          rank 0,1 for logical
!          rank 0,1,2,3 for integer, real(double), and complex(double)
!       same type as data :: tot
!          same rank as data with one exception
!          exception: data has rank 2, tot has rank 3 and both are complex(double)
!       integer :: i
!          only passed if data has rank 2 and tot has rank 3
!       requires: data has the same rank and dimensions on all processes.
!       modifies: tot
!       effects: tot on each process is overwritten by data from all sc processes, reduced
!                by op. If data has rank 2 complex(double) and tot has rank 3 complex(double),
!                the i'th rank 2 instance of tot is overwritten.

!cod$
 
      subroutine allreduce_log0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        logical :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_log0")) continue
      end subroutine

      subroutine allreduce_log1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        logical, dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_LOGICAL,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_LOGICAL,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_LOGICAL,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_LOGICAL,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_log1")) continue
      end subroutine

      subroutine allreduce_int0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_int0")) continue
      end subroutine

      subroutine allreduce_int1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_int1")) continue
      end subroutine

      subroutine allreduce_int2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_INTEGER,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_int2")) continue
      end subroutine

      subroutine allreduce_int3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_INTEGER,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_int3")) continue
      end subroutine

      subroutine allreduce_dpr0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double) :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpr0")) continue
      end subroutine

      subroutine allreduce_dpr1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpr1")) continue
      end subroutine

      subroutine allreduce_dpr2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpr2")) continue
      end subroutine

      subroutine allreduce_dpr3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpr3")) continue
      end subroutine

      subroutine allreduce_dpc0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double) :: data, tot
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_COMPLEX,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpc0")) continue
      end subroutine

      subroutine allreduce_dpc1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpc1")) continue
      end subroutine

      subroutine allreduce_dpc2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpc2")) continue
      end subroutine

      subroutine allreduce_dpc3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpc3")) continue
      end subroutine

      subroutine allreduce_dpc4(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpc4")) continue
      end subroutine

      subroutine allreduce_dpc2_3(sc,op,data,tot,i)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:) :: data
        complex(double), dimension(1:,1:,1:) :: tot
        integer :: i
        integer :: n
        n = size(data)
        select case (sc)
        case (WORLD)
          call MPI_ALLREDUCE(data(1,1),tot(1,1,i),n,MPI_DOUBLE_COMPLEX,op,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLREDUCE(data(1,1),tot(1,1,i),n,MPI_DOUBLE_COMPLEX,op,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1,i),n,MPI_DOUBLE_COMPLEX,op,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLREDUCE(data(1,1),tot(1,1,i),n,MPI_DOUBLE_COMPLEX,op,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allreduce_dpc2_3")) continue
      end subroutine

!doc$ subroutine allgather(sc,data,tot)
!       integer, intent(in) (WORLD, CONFIG, SGROUP, or KGROUP) :: sc
!       integer :: data
!       integer, dimension(:) :: tot
!       requires: size(tot) be the same on all sc processes.
!       modifies: tot
!       effects: Gathers data elements into tot across sc.

!cod$

      subroutine allgather_int0(sc,data,tot)
        integer, intent(in) :: sc
        integer :: data
        integer, dimension(:) :: tot
        select case (sc)
        case (WORLD)
          call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgather_int0")) continue
      end subroutine

!doc$ subroutine allgatherv(sc,data,tot,count,disp)
!       integer, intent(in) (WORLD, CONFIG, SGROUP, or KGROUP) :: sc
!       integer, real(double) rank 1,2,3; complex(double) rank 1 :: data
!       same type and rank as data :: tot
!       integer, dimension(:), intent(in) :: count, disp
!       requires: size(tot) >= sum(size(data))
!       modifies: tot
!       effects: Gathers data elements into tot across sc.

!cod$

      subroutine allgatherv_int1(sc,data,tot,count,disp)
        integer, intent(in) :: sc
        integer, dimension(:) :: data, tot
        integer, dimension(:), intent(in) :: count, disp
        integer :: n
        integer :: dummy(1)
        n = size(data)
        if (n < 1) then
          select case (sc)
          case (WORLD)
            call MPI_ALLGATHERV(dummy,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%world_comm,the_mpi%error)
          case (CONFIG)
            call MPI_ALLGATHERV(dummy,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%config_comm,the_mpi%error)
          case (SGROUP)
            call MPI_ALLGATHERV(dummy,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%sgroup_comm,the_mpi%error)
          case (KGROUP)
            call MPI_ALLGATHERV(dummy,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%kgroup_comm,the_mpi%error)
          end select
        else
          select case (sc)
          case (WORLD)
            call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%world_comm,the_mpi%error)
          case (CONFIG)
            call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%config_comm,the_mpi%error)
          case (SGROUP)
            call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%sgroup_comm,the_mpi%error)
          case (KGROUP)
            call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%kgroup_comm,the_mpi%error)
          end select
        end if
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgatherv_int1")) continue
      end subroutine

      subroutine allgatherv_int2(sc,data,tot,count,disp)
        integer, intent(in) :: sc
        integer, dimension(:,:) :: data, tot
        integer, dimension(:), intent(in) :: count, disp
        integer :: n
        n = size(data)
        select case (sc)
        case (WORLD)
          call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLGATHERV(data,n,MPI_INTEGER,tot,count,disp,MPI_INTEGER,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgatherv_int2")) continue
      end subroutine

      subroutine allgatherv_dpr1(sc,data,tot,count,disp)
        integer, intent(in) :: sc
        real(double), dimension(:) :: data, tot
        integer, dimension(:), intent(in) :: count, disp
        integer :: n
        real(double) :: dummy(1)
        n = size(data)
        if (n < 1) then
          select case (sc)
          case (WORLD)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%world_comm,the_mpi%error)
          case (CONFIG)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%config_comm,the_mpi%error)
          case (SGROUP)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%sgroup_comm,the_mpi%error)
          case (KGROUP)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%kgroup_comm,the_mpi%error)
          end select
        else
          select case (sc)
          case (WORLD)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%world_comm,the_mpi%error)
          case (CONFIG)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%config_comm,the_mpi%error)
          case (SGROUP)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%sgroup_comm,the_mpi%error)
          case (KGROUP)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%kgroup_comm,the_mpi%error)
          end select
        end if
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgatherv_dpr1")) continue
      end subroutine

      subroutine allgatherv_dpr2(sc,data,tot,count,disp)
        integer, intent(in) :: sc
        real(double), dimension(:,:) :: data, tot
        integer, dimension(:), intent(in) :: count, disp
        integer :: n
        n = size(data)
        select case (sc)
        case (WORLD)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgatherv_dpr2")) continue
      end subroutine

      subroutine allgatherv_dpr3(sc,data,tot,count,disp)
        integer, intent(in) :: sc
        real(double), dimension(:,:,:) :: data, tot
        integer, dimension(:), intent(in) :: count, disp
        integer :: n
        n = size(data)
        select case (sc)
        case (WORLD)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLGATHERV(data,n,MPI_DOUBLE_PRECISION,tot,count,disp,MPI_DOUBLE_PRECISION,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgatherv_dpr3")) continue
      end subroutine

      subroutine allgatherv_dpc1(sc,data,tot,count,disp)
        integer, intent(in) :: sc
        complex(double), dimension(:) :: data, tot
        integer, dimension(:), intent(in) :: count, disp
        integer :: n
        complex(double) :: dummy(1)
        n = size(data)
        if (n < 1) then
          select case (sc)
          case (WORLD)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%world_comm,the_mpi%error)
          case (CONFIG)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%config_comm,the_mpi%error)
          case (SGROUP)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%sgroup_comm,the_mpi%error)
          case (KGROUP)
            call MPI_ALLGATHERV(dummy,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%kgroup_comm,the_mpi%error)
          end select
        else
          select case (sc)
          case (WORLD)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%world_comm,the_mpi%error)
          case (CONFIG)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%config_comm,the_mpi%error)
          case (SGROUP)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%sgroup_comm,the_mpi%error)
          case (KGROUP)
            call MPI_ALLGATHERV(data,n,MPI_DOUBLE_COMPLEX,tot,count,disp,MPI_DOUBLE_COMPLEX,the_mpi%kgroup_comm,the_mpi%error)
          end select
        end if
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) continue
        if (error("Exit mpi_mod::allgatherv_dpc1")) continue
      end subroutine

!doc$ subroutine alltoallv(sc,data_in,count_in,disp_in,data_out,count_out,disp_out)
!       integer, intent(in) (WORLD, CONFIG, SGROUP, or KGROUP) :: sc
!       complex(double) type with rank 1 or 2 :: data_in
!       same type as data with rank 2 or 1 :: data_out
!       integer, dimension(:), intent(in) :: count_in, disp_in, count_out, disp_out
!       modifies: data_out
!       effects: Disperses data_in to data_out according to mpi_alltoallv routine.

!cod$

      subroutine alltoallv_dpc1_2(sc,data_in,count_in,disp_in,data_out,count_out,disp_out)
        integer, intent(in) :: sc
        complex(double), dimension(:) :: data_in
        complex(double), dimension(:,:) :: data_out
        integer, dimension(:), intent(in) :: count_in, disp_in, count_out, disp_out
        select case (sc)
        case (WORLD)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_ALLTOALLV,the_mpi%error")) continue
        if (error("Exit mpi_mod::alltoallv_dpc1_2")) continue
      end subroutine

      subroutine alltoallv_dpc2_1(sc,data_in,count_in,disp_in,data_out,count_out,disp_out)
        integer, intent(in) :: sc
        complex(double), dimension(:,:) :: data_in
        complex(double), dimension(:) :: data_out
        integer, dimension(:), intent(in) :: count_in, disp_in, count_out, disp_out
        select case (sc)
        case (WORLD)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%world_comm,the_mpi%error)
        case (CONFIG)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%config_comm,the_mpi%error)
        case (SGROUP)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%sgroup_comm,the_mpi%error)
        case (KGROUP)
          call MPI_ALLTOALLV(data_in,count_in,disp_in,MPI_DOUBLE_COMPLEX, &
                             data_out,count_out,disp_out,MPI_DOUBLE_COMPLEX,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_ALLTOALLV,the_mpi%error")) continue
        if (error("Exit mpi_mod::alltoallv_dpc2_1")) continue
      end subroutine

!doc$ subroutine blocking_send(sc,data,dest,tag)
!       integer, intent(in) (KGROUP) :: sc
!       integer, intent(in) :: data
!       complex(double), dimension(:), intent(in) :: data
!       integer, intent(in) :: dest, tag
!       effects: Sends data to dest in sc with ID tag, waiting for completion.

!cod$

      subroutine blocking_send_int0(sc,data,dest,tag)
        integer, intent(in) :: sc
        integer, intent(in) :: data
        integer, intent(in) :: dest, tag
        integer :: count
        count = 1
        select case (sc)
        case (KGROUP)
          call MPI_SEND(data,count,MPI_INTEGER,dest,tag,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_SEND,the_mpi%error")) continue
        if (error("Exit mpi_mod::blocking_send_int0")) continue
      end subroutine

      subroutine blocking_send_dpc1(sc,data,dest,tag)
        integer, intent(in) :: sc
        complex(double), dimension(:), intent(in) :: data
        integer, intent(in) :: dest, tag
        integer :: count
        count = size(data)
        select case (sc)
        case (KGROUP)
          call MPI_SEND(data,count,MPI_DOUBLE_COMPLEX,dest,tag,the_mpi%kgroup_comm,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_SEND,the_mpi%error")) continue
        if (error("Exit mpi_mod::blocking_send_dpc1")) continue
      end subroutine

!doc$ subroutine nonblocking_send(sc,data,dest,tag)
!       integer, intent(in) (KGROUP) :: sc
!       integer, intent(in) :: data
!       complex(double), dimension(:), intent(in) :: data
!       integer, intent(in) :: dest, tag
!       effects: Sends data to dest in sc with ID tag, not waiting for completion.

!cod$

      subroutine nonblocking_send_int0(sc,data,dest,tag)
        integer, intent(in) :: sc
        integer, intent(in) :: data
        integer, intent(in) :: dest, tag
        integer :: count, request
        count = 1
        select case (sc)
        case (KGROUP)
          call MPI_ISEND(data,count,MPI_INTEGER,dest,tag,the_mpi%kgroup_comm,request,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_ISEND,the_mpi%error")) continue
        if (error("Exit mpi_mod::nonblocking_send_int0")) continue
      end subroutine

      subroutine nonblocking_send_dpc1(sc,data,dest,tag)
        integer, intent(in) :: sc
        complex(double), dimension(:), intent(in) :: data
        integer, intent(in) :: dest, tag
        integer :: count, request
        count = size(data)
        select case (sc)
        case (KGROUP)
          call MPI_ISEND(data,count,MPI_DOUBLE_COMPLEX,dest,tag,the_mpi%kgroup_comm,request,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_ISEND,the_mpi%error")) continue
        if (error("Exit mpi_mod::nonblocking_send_dpc1")) continue
      end subroutine

!doc$ subroutine blocking_recv(sc,data,origin,tag)
!       integer, intent(in) (KGROUP) :: sc
!       integer, intent(in) :: data
!       complex(double), dimension(:), intent(in) :: data
!       integer, intent(in) :: origin, tag
!       effects: Receives data from origin in sc with ID tag, waiting for completion.

!cod$

      subroutine blocking_recv_int0(sc,data,origin,tag)
        integer, intent(in) :: sc
        integer, intent(in) :: data
        integer, intent(in) :: origin, tag
        integer :: count
        integer status(MPI_STATUS_SIZE)
        count = 1
        select case (sc)
        case (KGROUP)
          call MPI_RECV(data,count,MPI_INTEGER,origin,tag,the_mpi%kgroup_comm,status,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_RECV,the_mpi%error")) continue
        if (error("Exit mpi_mod::blocking_recv_int0")) continue
      end subroutine

      subroutine blocking_recv_dpc1(sc,data,origin,tag)
        integer, intent(in) :: sc
        complex(double), dimension(:), intent(in) :: data
        integer, intent(in) :: origin, tag
        integer :: count
        integer status(MPI_STATUS_SIZE)
        count = size(data)
        select case (sc)
        case (KGROUP)
          call MPI_RECV(data,count,MPI_DOUBLE_COMPLEX,origin,tag,the_mpi%kgroup_comm,status,the_mpi%error)
        end select
        if (error(the_mpi%error /= 0,"ERROR: Non-zero error status from MPI_RECV,the_mpi%error")) continue
        if (error("Exit mpi_mod::blocking_recv_dpc1")) continue
      end subroutine

!doc$ subroutine xcomm_broadcast(sc,root,data)
!       integer, intent(in) (XCONFIG, XSGROUP, or XKGROUP) :: sc
!       integer, intent(in) :: root
!       complex(double), rank 2 :: data
!       modifies: data
!       effects: data on (sc,root) is broadcast to all xsc processes.

!cod$
 
      subroutine xcomm_broadcast_dpc2(sc,root,data)
        integer, intent(in) :: sc, root
        complex(double), dimension(1:,1:) :: data
        integer :: n
        n = size(data)
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() > 1) then
            call MPI_BCAST(data(1,1),n,MPI_DOUBLE_COMPLEX,root,the_mpi%xconfig_comm,the_mpi%error)
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (2)
            call MPI_BCAST(data(1,1),n,MPI_DOUBLE_COMPLEX,root,the_mpi%xsgroup_comm,the_mpi%error)
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() > 1) then
            call MPI_BCAST(data(1,1),n,MPI_DOUBLE_COMPLEX,root,the_mpi%xkgroup_comm,the_mpi%error)
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_broadcast_dpc2")) continue
      end subroutine

!doc$ subroutine xcomm_reduce(sc,op,data,tot)
!       integer, intent(in) (XCONFIG, XSGROUP, or XKGROUP) :: sc
!       integer (MPI_SUM) :: op
!       integer, rank 0, 1; real(double), rank 1, 3; complex(double), rank 1 :: data, tot
!       modifies: tot
!       effects: tot on (sc,0) is overwritten by data from (sc-1,0) reduced by op.

!cod$
 
      subroutine xcomm_reduce_int0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer :: data, tot
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = 0
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = 0
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = 0
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = 0
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = 0
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data,tot,1,MPI_INTEGER,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = 0
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_int0")) continue
      end subroutine

      subroutine xcomm_reduce_int1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = 0
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = 0
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = 0
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = 0
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = 0
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_INTEGER,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = 0
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_int1")) continue
      end subroutine

      subroutine xcomm_reduce_dpr1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpr1")) continue
      end subroutine

      subroutine xcomm_reduce_dpr3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpr3")) continue
      end subroutine

      subroutine xcomm_reduce_dpr4(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = 0.0_double
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = 0.0_double
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpr4")) continue
      end subroutine

      subroutine xcomm_reduce_dpc1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpc1")) continue
      end subroutine

      subroutine xcomm_reduce_dpc2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpc2")) continue
      end subroutine

      subroutine xcomm_reduce_dpc3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpc3")) continue
      end subroutine

      subroutine xcomm_reduce_dpc4(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            if (the_mpi%world_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%config_first) then
              call MPI_REDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xconfig_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            if (the_mpi%config_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_REDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xsgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            if (the_mpi%sgroup_first) then
              tot = data
            else
              tot = (0.0_double,0.0_double)
            end if
          else
            if (the_mpi%kgroup_first) then
              call MPI_REDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_COMPLEX,op,0,the_mpi%xkgroup_comm,the_mpi%error)
            else
              tot = (0.0_double,0.0_double)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_reduce_dpc4")) continue
      end subroutine

!doc$ subroutine xcomm_allreduce(sc,op,data,tot)
!       integer, intent(in) (XCONFIG, XSGROUP, or XKGROUP) :: sc
!       integer (MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN, MPI_LAND, MPI_LOR, or MPI_LXOR) :: op
!       logical, rank 0; integer, rank 1; real(double), rank 0, 1, 2, or 3; complex(double), rank 1 or 3 :: data, tot
!       modifies: tot
!       effects: tot on (sc,*) is overwritten by data from (sc-1,0) reduced by op.

!cod$
 
      subroutine xcomm_allreduce_log0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        logical :: data, tot
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_LOGICAL,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_log0")) continue
      end subroutine

      subroutine xcomm_allreduce_int0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer :: data, tot
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_INTEGER,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_int0")) continue
      end subroutine

      subroutine xcomm_allreduce_int1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        integer, dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_INTEGER,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_int1")) continue
      end subroutine

      subroutine xcomm_allreduce_dpr0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double) :: data, tot
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpr0")) continue
      end subroutine

      subroutine xcomm_allreduce_dpr1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpr1")) continue
      end subroutine

      subroutine xcomm_allreduce_dpr2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpr2")) continue
      end subroutine

      subroutine xcomm_allreduce_dpr3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpr3")) continue
      end subroutine

      subroutine xcomm_allreduce_dpr4(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpr4")) continue
      end subroutine




      subroutine xcomm_allreduce_dpc1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpc1")) continue
      end subroutine

      subroutine xcomm_allreduce_dpc3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allreduce_dpc3")) continue
      end subroutine

!doc$ subroutine xcomm_allgather(sc,data,tot)
!       integer, intent(in) (XCONFIG, XSGROUP, or XKGROUP) :: sc
!       integer, rank 0; real(double), rank 0, 1, 2 :: data
!       integer, rank 1; real(double), rank 1, 2, 3 :: tot
!       modifies: tot
!       effects: Gathers data from processes (sc,0) into tot on processes (sc,*).

!cod$
 
      subroutine xcomm_allgather_int0(sc,data,tot)
        integer, intent(in) :: sc
        integer :: data
        integer, dimension(:) :: tot
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot(1) = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot(1) = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot(1) = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLGATHER(data,1,MPI_INTEGER,tot,1,MPI_INTEGER,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allgather_int0")) continue
      end subroutine

      subroutine xcomm_allgather_dpr0(sc,data,tot)
        integer, intent(in) :: sc
        real(double) :: data
        real(double), dimension(:) :: tot
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot(1) = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLGATHER(data,1,MPI_DOUBLE_PRECISION,tot,1,MPI_DOUBLE_PRECISION,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot(1) = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLGATHER(data,1,MPI_DOUBLE_PRECISION,tot,1,MPI_DOUBLE_PRECISION,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot(1) = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLGATHER(data,1,MPI_DOUBLE_PRECISION,tot,1,MPI_DOUBLE_PRECISION,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allgather_dpr0")) continue
      end subroutine

      subroutine xcomm_allgather_dpr1(sc,data,tot)
        integer, intent(in) :: sc
        real(double), dimension(:) :: data
        real(double), dimension(:,:) :: tot
        integer :: n
        n = size(data)
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot(:,1) = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLGATHER(data,n,MPI_DOUBLE_PRECISION,tot,n,MPI_DOUBLE_PRECISION,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot(:,1) = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLGATHER(data,n,MPI_DOUBLE_PRECISION,tot,n,MPI_DOUBLE_PRECISION,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot(:,1) = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLGATHER(data,n,MPI_DOUBLE_PRECISION,tot,n,MPI_DOUBLE_PRECISION,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xconfig_allgather_dpr1")) continue
      end subroutine

      subroutine xcomm_allgather_dpr2(sc,data,tot)
        integer, intent(in) :: sc
        real(double), dimension(:,:) :: data
        real(double), dimension(:,:,:) :: tot
        integer :: n
        n = size(data)
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot(:,:,1) = data
          else
            if (the_mpi%config_first) then
              call MPI_ALLGATHER(data,n,MPI_DOUBLE_PRECISION,tot,n,MPI_DOUBLE_PRECISION,the_mpi%xconfig_comm,the_mpi%error)
            end if
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(CONFIG,tot)
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot(:,:,1) = data
          case (2)
            if (the_mpi%sgroup_first) then
              call MPI_ALLGATHER(data,n,MPI_DOUBLE_PRECISION,tot,n,MPI_DOUBLE_PRECISION,the_mpi%xsgroup_comm,the_mpi%error)
            end if
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(SGROUP,tot)
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot(:,:,1) = data
          else
            if (the_mpi%kgroup_first) then
              call MPI_ALLGATHER(data,n,MPI_DOUBLE_PRECISION,tot,n,MPI_DOUBLE_PRECISION,the_mpi%xkgroup_comm,the_mpi%error)
            end if
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
            call broadcast(KGROUP,tot)
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_allgather_dpr2")) continue
      end subroutine

!doc$ subroutine xcomm_rank_allreduce(sc,op,data,tot)
!       integer, intent(in) (XCONFIG, XSGROUP, or XKGROUP) :: sc
!       integer (MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN, MPI_LAND, MPI_LOR, or MPI_LXOR) :: op
!       real(double), rank 1, 2, 3 or 4; complex(double), rank 3 :: data, tot
!       modifies: tot
!       effects: allreduce operation across sc for each rank followed by a broadcast of the_mpi%error.

!cod$

      subroutine xcomm_rank_allreduce_dpr1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_rank_allreduce_dpr1")) continue
      end subroutine

      subroutine xcomm_rank_allreduce_dpr2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_rank_allreduce_dpr2")) continue
      end subroutine

      subroutine xcomm_rank_allreduce_dpr3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_rank_allreduce_dpr3")) continue
      end subroutine

      subroutine xcomm_rank_allreduce_dpr4(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xconfig_comm,the_mpi%error)
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1,1,1),tot(1,1,1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xkgroup_comm,the_mpi%error)
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_rank_allreduce_dpr4")) continue
      end subroutine

      subroutine xcomm_rank_allreduce_dpc3(sc,op,data,tot)
        integer, intent(in) :: sc, op
        complex(double), dimension(1:,1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XCONFIG)
          if (mpi_nconfigs() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xconfig_comm,the_mpi%error)
            call broadcast(CONFIG,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xsgroup_comm,the_mpi%error)
            call broadcast(SGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case (XKGROUP)
          if (mpi_nkgroups() == 1) then
            tot = data
          else
            call MPI_ALLREDUCE(data(1,1,1),tot(1,1,1),n,MPI_DOUBLE_COMPLEX,op,the_mpi%xkgroup_comm,the_mpi%error)
            call broadcast(KGROUP,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end if
        end select
100     if (error("Exit mpi_mod::xcomm_rank_allreduce_dpc3")) continue
      end subroutine

!doc$ subroutine xcomm_pair_allreduce(sc,op,data,tot)
!       integer, intent(in) (XSGROUP) :: sc
!       integer (MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN, MPI_LAND, MPI_LOR, or MPI_LXOR) :: op
!       real(double), rank 0, 1 or 2 :: data, tot
!       modifies: tot
!       effects: allreduce operation across sc for each xsc pair.
!       notes: This routine performs the same function as xcomm_rank_allreduce, but does not broadcast
!              the_mpi%error. Thus, it can be used in cases where not all processes participate.

!cod$

      subroutine xcomm_pair_allreduce_dpr0(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double) :: data, tot
        select case (sc)
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data,tot,1,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case default
          if (error(.true.,"ERROR: sc is not equal to XSGROUP")) continue
        end select
100     if (error("Exit mpi_mod::xcomm_pair_allreduce_dpr0")) continue
      end subroutine

      subroutine xcomm_pair_allreduce_dpr1(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1),tot(1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case default
          if (error(.true.,"ERROR: sc is not equal to XSGROUP")) continue
        end select
100     if (error("Exit mpi_mod::xcomm_pair_allreduce_dpr1")) continue
      end subroutine

      subroutine xcomm_pair_allreduce_dpr2(sc,op,data,tot)
        integer, intent(in) :: sc, op
        real(double), dimension(1:,1:) :: data, tot
        integer :: n
        n = min(size(data),size(tot))
        select case (sc)
        case (XSGROUP)
          select case (mpi_nsgroups())
          case (1)
            tot = data
          case (2)
            call MPI_ALLREDUCE(data(1,1),tot(1,1),n,MPI_DOUBLE_PRECISION,op,the_mpi%xsgroup_comm,the_mpi%error)
            if (error(the_mpi%error /= 0,"ERROR: Non-zero error status",the_mpi%error)) goto 100
          end select
        case default
          if (error(.true.,"ERROR: sc is not equal to XSGROUP")) continue
        end select
100     if (error("Exit mpi_mod::xcomm_pair_allreduce_dpr2")) continue
      end subroutine

      end module

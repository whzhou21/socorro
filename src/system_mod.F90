!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module system_mod
!doc$ module system_mod

!     sys_mod provides routines for starting and stopping the runtime environment.

      use kind_mod
      use mpi_mod
      use path_mod
      use arg_mod
      use error_mod
      use interrupt_mod
      use io_mod
      use diary_mod
      use fft_mod
      use math_mod
      use timing_mod
      use point_blas_mod

!cod$
      implicit none ; private

!doc$
      public :: system_start
      public :: system_stop

!cod$

      interface system_start
         module procedure system_start_
      end interface

      interface system_stop
         module procedure system_stop_
      end interface

      contains

      subroutine system_start_()
!doc$ subroutine system_start()
!       effects: Calls routines that start the runtime environment.
!       errors: configurations < 1 or > 98. Passes errors.
!       notes: Errors here and in called routines cause program to halt.

!cod$
        logical :: found
        character(line_len) :: ef_status, tag
        integer :: nc, ng

        ! start the mpi system
        call mpi__start()      

        ! read the arguments file
        call arg_start()

        ! set the number of configurations
        call arg('configurations',nc,found)
        if (.not.found) nc = 1
        if (nc < 1) then
          if (mpi_first(WORLD)) write(6,'("ERROR: configurations < 1")')
          call mpi_stop() ; stop
        elseif (nc == 1) then
          continue
        elseif ((nc > 1) .and. (nc < 99)) then
          if (mod(mpi_nprocs(WORLD),nc) /= 0) then
            if (mpi_first(WORLD)) write(6,'("ERROR: non-equal division of world processes among configurations")')
            call mpi_stop() ; stop
          end if
          call mpi_config_split(nc)
        else
          if (mpi_first(WORLD)) write(6,'("ERROR: configurations > 98")')
          call mpi_stop() ; stop
        end if

        ! determine the transition-state method
        call arglc('ts_method',tag,found)
        if (.not.found) tag = "none"
        select case (trim(tag))
        case ("dimer","dmr")
          if (nc /= 2) then
            if (mpi_first(WORLD)) write(6,'("ERROR: configurations /= 2 for a dimer calculation")')
            call mpi_stop() ; stop
          end if
        case ("neb")
          if (nc == 1) then
            if (mpi_first(WORLD)) write(6,'("ERROR: configurations == 1 for an neb calculation")')
            call mpi_stop() ; stop
          end if
        end select

        ! set the number of spin groups
        call arg('sgroups',ng,found)
        if (.not.found) ng = 1
        if (ng < 1) then
          if (mpi_first(WORLD)) write(6,'("ERROR: sgroups < 1")')
          call mpi_stop() ; stop
        elseif (ng == 1) then
          continue
        elseif (ng == 2) then
          if (mod(mpi_nprocs(CONFIG),ng) /= 0) then
            if (mpi_first(WORLD)) write(6,'("ERROR: non-equal division of config processes among sgroups")')
            call mpi_stop() ; stop
          end if
          call mpi_sgroup_split(ng)
        else
          if (mpi_first(WORLD)) write(6,'("ERROR: sgroups > 2")')
          call mpi_stop() ; stop
        end if

        ! set the number of k-point groups
        call arg('kgroups',ng,found)
        if (.not.found) ng = 1
        if (ng < 1) then
          if (mpi_first(WORLD)) write(6,'("ERROR: kgroups < 1")')
          call mpi_stop() ; stop
        elseif (ng == 1) then
          continue
        else
          if (mod(mpi_nprocs(SGROUP),ng) /= 0) then
            if (mpi_first(WORLD)) write(6,'("ERROR: non-equal division of sgroup processes among kgroups")')
            call mpi_stop() ; stop
          end if
          call mpi_kgroup_split(ng)
        end if

        ! set the file paths
        call set_paths(mpi_nconfigs(),mpi_myconfig(),mpi_myproc(CONFIG))

        ! start the error-handling system
        call arglc('error_file_mode',tag,found)
        if (.not.found) tag = "first"
        select case (trim(tag))
        case ("all")
          ef_status = "on"
        case ("first")
          ef_status = "off"
          if (mpi_first(CONFIG)) ef_status = "on"
        case ("none")
          ef_status = "off"
        case default
          if (mpi_first(WORLD)) write(6,'("ERROR: error_file_mode tag must be all, first, or none")')
          call mpi_stop() ; stop
        end select
        call error_start(mpi_comm(CONFIG),mpi_comm(SGROUP),mpi_comm(KGROUP),mpi_myproc(CONFIG),mpi_first(WORLD),tag,ef_status)

        ! write MPI information to the error files
        select case (trim(tag))
        case ("all")
          call warn("MPI information:")
          call warn(" ")
          call warn("WORLD")
          call notify("communicator         ",mpi_comm(WORLD))
          call notify("processes            ",mpi_nprocs(WORLD))
          call notify("  rank               ",mpi_myproc(WORLD))
          call warn(" ")
          call warn("CONFIG")
          call notify("number of configs  ",mpi_nconfigs())
          call notify("communicator       ",mpi_comm(CONFIG))
          call notify("processes          ",mpi_nprocs(CONFIG))
          call notify("  config number    ",mpi_myconfig())
          call notify("  rank             ",mpi_myproc(CONFIG))
          call warn(" ")
          call warn("  XCONFIG")
          call notify("  communicator           ",mpi_comm(XCONFIG))
          call notify("  processes              ",mpi_nprocs(XCONFIG))
          call notify("    rank                 ",mpi_myproc(XCONFIG))
          call warn(" ")
          call warn("SGROUP")
          call notify("number of sgroups  ",mpi_nsgroups())
          call notify("communicator       ",mpi_comm(SGROUP))
          call notify("processes          ",mpi_nprocs(SGROUP))
          call notify("  sgroup number    ",mpi_mysgroup())
          call notify("  rank             ",mpi_myproc(SGROUP))
          call warn(" ")
          call warn("  XSGROUP")
          call notify("  communicator           ",mpi_comm(XSGROUP))
          call notify("  processes              ",mpi_nprocs(XSGROUP))
          call notify("    rank                 ",mpi_myproc(XSGROUP))
          call warn(" ")
          call warn("KGROUP")
          call notify("number of kgroups  ",mpi_nkgroups())
          call notify("communicator       ",mpi_comm(KGROUP))
          call notify("processes          ",mpi_nprocs(KGROUP))
          call notify("  kgroup number    ",mpi_mykgroup())
          call notify("  rank             ",mpi_myproc(KGROUP))
          call warn(" ")
          call warn("  XKGROUP")
          call notify("  communicator           ",mpi_comm(XKGROUP))
          call notify("  processes              ",mpi_nprocs(XKGROUP))
          call notify("    rank                 ",mpi_myproc(XKGROUP))
          call warn(" ")
        end select

        ! start the input/output system
        call io_start() ; if(error()) goto 100

        ! initialize the diary file
        call diary_start()
        call diary_socorro_env()
        call sync_configuration_errors() ; if (error()) goto 100

        ! initialize the fft threads
        call fft_start()

        call init_machine_constants()
        call inittimers()

        ! zero timers (or risk getting incorrect timing)

!        call start_timer("config_sc: constructor")
!        call stop_timer("config_sc: constructor")

!        call start_timer("config_sc: update")
!        call stop_timer("config_sc: update")

!        call start_timer("eigensolver: diagonalize")
!        call stop_timer("eigensolver: diagonalize")

!        call start_timer("eigensolver: diagonalize_gen")
!        call stop_timer("eigensolver: diagonalize_gen")

!        call start_timer("electrons_sc: constructor")
!        call stop_timer("electrons_sc: constructor")

!        call start_timer("electrons_sc: update")
!        call stop_timer("electrons_sc: update")

        call start_timer("fft: distributed")
        call stop_timer("fft: distributed")

        call start_timer("fft: serial")
        call stop_timer("fft: serial")

!        call start_timer("fields_sc: constructor")
!        call stop_timer("fields_sc: constructor")

!        call start_timer("fields_sc: update")
!        call stop_timer("fields_sc: update")

        call start_timer("multivector: zgemm")
        call stop_timer("multivector: zgemm")

        call start_timer("multivector: gather_wf")
        call stop_timer("multivector: gather_wf")

        call start_timer("multivector: scatter_wf")
        call stop_timer("multivector: scatter_wf")

!        call start_timer("multivector: exchange_operator_mv")
!        call stop_timer("multivector: exchange_operator_mv")

!        call start_timer("multivector: exchange_operator_2mv")
!        call stop_timer("multivector: exchange_operator_2mv")

!        call start_timer("multivector: p2p calculation")
!        call stop_timer("multivector: p2p calculation")

!        call start_timer("multivector: p2p calculation barrier")
!        call stop_timer("multivector: p2p calculation barrier")

!        call start_timer("multivector: p2p communication")
!        call stop_timer("multivector: p2p communication")

!        call start_timer("multivector: kernel_mm_cr")
!        call stop_timer("multivector: kernel_mm_cr")

!        call start_timer("multivector: kernel_mmm_ccc")
!        call stop_timer("multivector: kernel_mmm_ccc")

!        call start_timer("multivector: kernel_rsum_mmm_ccc")
!        call stop_timer("multivector: kernel_rsum_mmm_ccc")

!        call start_timer("multivector: kernel_csum_portion_mmm_ccc_rs")
!        call stop_timer("multivector: kernel_csum_portion_mmm_ccc_rs")

        call start_timer("multibasis: band_remap")
        call stop_timer("multibasis: band_remap")

!        call start_timer("multibasis: spair_remap")
!        call stop_timer("multibasis: spair_remap")

!        call start_timer("multibasis: lpair_remap")
!        call stop_timer("multibasis: lpair_remap")

        call start_timer("operators: zgemm")
        call stop_timer("operators: zgemm")

        call start_timer("operators: apply_hamiltonian")
        call stop_timer("operators: apply_hamiltonian")

        ! The tag point_mode is currently hardwired to the value f90 so that point
        ! blas operations are performed with (normal) Fortran 90 constructs. It is
        ! noted that, consistent with this hardwiring, the point_mode tag is not
        ! described in the README file.
        call arglc('point_mode',tag,found)
        if (.not.found) tag = "f90"
        select case (trim(tag))
        case ("f77")
          call point_mode(POINT_F77)
        case ("f90")
          call point_mode(POINT_F90)
        end select

        call random_seed()

        call start_timer("Socorro: total time")

100      if ( error("Exiting system_mod::system_start_()") ) continue

      end subroutine system_start_

      subroutine system_stop_()
!doc$ subroutine system_stop()
!        effects: Calls routines to stop the runtime system.

!cod$
         if ( .not.error() ) then
            call stop_timer("Socorro: total time")
            call print_all_timers()
         end if

         call interrupt()
         call diary_stop()
         call io_stop()
         call error_stop()
         call arg_stop()
	 call kill_all_timers()
         call mpi_stop()

      end subroutine system_stop_

      end module system_mod

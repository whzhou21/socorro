!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module error_mod
!doc$ module error_mod

!     error_mod: - maintains an error file for each process.
!                - maintains variables representing the configuration error state and the error state for each process.
!                - provides a routine for modifying the process and configuration error states. The routine can be used
!                    to evaluate a test condition that may or may not be synchronized among the processes.
!                - provides a routine for checking the configuration error state and optionally posting a message to the
!                     error file if the configuration error state is set.
!                - provides routines for posting messages to the error file independent of the configuration error state.
!                - provides routines for synchronizing the process error states and the configuration error states.
!                - currently does not provide routines for resetting the error states (once true, always true).

!     Note: mpif.h is included here to gain access to MPI routines and thereby allow mpi_mod to use error_mod.

!     Note: This module should not be changed without first consulting AFW.

      use kind_mod
      use mpi
      use path_mod

!cod$
      implicit none ; private

      ! error file status
      integer, parameter :: EF_OFF = 0
      integer, parameter :: EF_ON  = 1

      ! error file mode
      integer, parameter :: EH_NONE  = 0
      integer, parameter :: EH_ALL   = 1
      integer, parameter :: EH_FIRST = 2

      integer :: file_mode           ! mode of error file handling
      integer :: file_status         ! status of process error files
      integer :: error_unit          ! unit number of the process error file
      logical :: kp_error_state      ! kgroup process error state
      logical :: k_error_state       ! whole kgroup error state
      integer :: kgroup_comm         ! MPI kgroup communicator
      logical :: sp_error_state      ! sgroup process error state
      logical :: s_error_state       ! whole sgroup error state
      integer :: sgroup_comm         ! MPI sgroup communicator
      logical :: cp_error_state      ! config process error state
      logical :: c_error_state       ! whole config error state
      integer :: config_comm         ! MPI config communicator

!doc$
      public :: error_start
      public :: error_stop
      public :: error
      public :: warn
      public :: notify
      public :: sync_kgroup_process_errors
      public :: sync_sgroup_process_errors
      public :: sync_config_process_errors
      public :: sync_configuration_errors

!cod$
      interface error
        module procedure error_test_msg, error_test_msg_n, error_msg
      end interface
      interface notify
        module procedure notify_msg_integer, notify_msg_real, notify_msg_complex
      end interface

      contains

      subroutine error_start(c_comm,s_comm,k_comm,c_myp,w_first,mode,status)
!doc$ subroutine error_start(c_comm,s_comm,k_comm,c_myp,w_first,mode,status)
        integer, intent(in) :: c_comm, s_comm, k_comm, c_myp
        logical :: w_first
        character(line_len) :: mode, status
!       requires: MPI be initialized. c_comm, s_comm, and k_comm be the MPI configuration, sgroup, and kgroup communicator.
!       effects: Starts the error handling system.
!       errors: Problems opening an error file.

!cod$
        logical :: p_ex, w_ex, p_op, w_op, p_switch, w_switch
        integer :: ios, mpi_error

        ! determine the error handling mode
        select case (trim(mode))
        case ("none")
          file_mode = EH_NONE
        case ("all")
          file_mode = EH_ALL
        case ("first")
          file_mode = EH_FIRST
        end select

        ! determine the process file status
        select case (trim(status))
        case ("off")
          file_status = EF_OFF
        case ("on")
          file_status = EF_ON
        end select

        select case (file_mode)
        case (EH_ALL)

          ! Find an error_unit common to all processes
          error_unit = first_unit
          do
            if (error_unit > last_unit) then
              if (w_first) write(6,'("ERROR: unable to find a unit for the error files")')
              call MPI_FINALIZE(mpi_error) ; stop
            end if
            inquire(error_unit,exist=p_ex)
            call MPI_ALLREDUCE(p_ex,w_ex,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_error)
            if (.not.w_ex) then
              error_unit = error_unit + 1
              cycle
            end if
            inquire(error_unit,opened=p_op)
            call MPI_ALLREDUCE(p_op,w_op,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_error)
            if (w_op) then
              error_unit = error_unit + 1
              cycle
            end if
            exit
          end do

          ! Open the error files.
          open(error_unit,file=trim(p_error_path),status='unknown',iostat=ios)
          if (ios /= 0) p_op = .true.
          call MPI_ALLREDUCE(p_op,w_op,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_error)
          if (w_op) then
            if (w_first) write(6,'("ERROR: unable to open an error file")')
            call MPI_FINALIZE(mpi_error) ; stop
          end if

          write(error_unit,'("PROCESS ",i4.4,/)') c_myp
          call flush_i()

        case (EH_FIRST)

          ! Find an error_unit common to all first processes
          error_unit = first_unit
          do
            if (error_unit > last_unit) then
              if (w_first) write(6,'("ERROR: unable to find a unit for the error files")')
              call MPI_FINALIZE(mpi_error) ; stop
            end if
            select case (file_status)
            case (EF_ON)
              inquire(error_unit,exist=p_ex)
            case (EF_OFF)
              p_ex = .false.
            end select
            call MPI_ALLREDUCE(p_ex,w_ex,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_error)
            if (.not.w_ex) then
              error_unit = error_unit + 1
              cycle
            end if
            select case (file_status)
            case (EF_OFF)
              p_op = .false.
            case (EF_ON)
              inquire(error_unit,opened=p_op)
            end select
            call MPI_ALLREDUCE(p_op,w_op,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_error)
            if (w_op) then
              error_unit = error_unit + 1
              cycle
            end if
            exit
          end do

          ! Open the error files.
          select case (file_status)
          case (EF_ON)
            open(error_unit,file=trim(f_error_path),status='unknown',iostat=ios)
            if (ios /= 0) p_op = .true.
          case (EF_OFF)
            p_op = .false.
          end select
          call MPI_ALLREDUCE(p_op,w_op,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_error)
          if (w_op) then
            if (w_first) write(6,'("ERROR: unable to open an error file")')
            call MPI_FINALIZE(mpi_error) ; stop
          end if

        end select

        config_comm = c_comm
        cp_error_state = .false.
        c_error_state = .false.

        sgroup_comm = s_comm
        sp_error_state = .false.
        s_error_state = .false.

        kgroup_comm = k_comm
        kp_error_state = .false.
        k_error_state = .false.

      end subroutine

      subroutine error_stop()
!doc$ subroutine error_stop()
!       effects: Stops the error handling system.

!cod$
        logical :: op

        select case (file_mode)
        case (EH_ALL)
          inquire(error_unit,opened=op)
          if (op) then
            if (cp_error_state) then
              write(error_unit,'(/,"ERRORS WERE DETECTED BY THIS PROCESS")')
            else
              write(error_unit,'(/,"NO ERRORS WERE DETECTED BY THIS PROCESS")')
            end if
            close(error_unit)
          end if
        case (EH_FIRST)
          select case (file_status)
          case (EF_ON)
            inquire(error_unit,opened=op)
            if (op) then
              if (cp_error_state) then
                write(error_unit,'(/,"ERRORS WERE DETECTED")')
              else
                write(error_unit,'(/,"NO ERRORS WERE DETECTED")')
              end if
              close(error_unit)
            end if
          end select
        end select

      end subroutine

      function error_test_msg(test,msg) result(ces)
!doc$ function error(test,msg) result(ces)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        logical :: ces
!       requires: MPI be initialized. test has the same value on all processes.
!       modifies: Process error state. Configuration error state.
!       effects: If test = .true., msg is written to the error file.
!                The process error state is set equal to (process error state .or. test).
!                ces is set equal to the configuration error state.
!       errors: Configuration error state = .true. on entry. MPI problem.

!cod$
        if (c_error_state) then
          select case (file_status)
          case (EF_ON)
            write(error_unit,'("ERROR: c_error_state = .true. on entry into error_mod::error_test_msg")')
            call flush_i()
          end select
          goto 100
        end if

        if (test) then
          select case (file_status)
          case (EF_ON)
            write(error_unit,'(a)') msg
            call flush_i()
          end select
        end if

        kp_error_state = (kp_error_state .or. test)
        k_error_state = kp_error_state

        sp_error_state = (sp_error_state .or. test)
        s_error_state = sp_error_state

        cp_error_state = (cp_error_state .or. test)
        c_error_state = cp_error_state

100     ces = c_error_state

      end function

      function error_test_msg_n(test,msg,n) result(ces)
!doc$ function error(test,msg,n) result(ces)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        integer, intent(in) :: n
        logical :: ces
!       requires: MPI be initialized.
!       modifies: Process error state. Configuration error state.
!       effects: If test = .true., msg and n are written to the error file.
!                The process error state is set equal to (process error state .or. test).
!                ces is set equal to the configuration error state.
!       errors: Configuration error state = .true. on entry. MPI problem.
!       notes: It is assumed that test has the same value on all processes calling the routine.

!cod$
        if (c_error_state) then
          select case (file_status)
          case (EF_ON)
            write(error_unit,'("ERROR: c_error_state = .true. on entry into error_mod::error_test_msg_n")')
            call flush_i()
          end select
          goto 100
        end if

        if (test) then
          select case (file_status)
          case (EF_ON)
            write(error_unit,'(a," ",i0)') msg, n
            call flush_i()
          end select
        end if

        kp_error_state = (kp_error_state .or. test)
        k_error_state = kp_error_state

        sp_error_state = (sp_error_state .or. test)
        s_error_state = sp_error_state

        cp_error_state = (cp_error_state .or. test)
        c_error_state = cp_error_state

100     ces = c_error_state

      end function

      function error_msg(msg) result(ces)
!doc$ function error(msg) result(ces)
        character(*), intent(in), optional :: msg
        logical :: ces
!       effects: If the configuration error state = .true., msg is written to the error file.
!                ces is set equal to the configuration error state.

!cod$
        if (present(msg)) then
          if (c_error_state) then
            select case (file_status)
            case (EF_ON)
              write(error_unit,'(a)') msg
              call flush_i()
            end select
          end if
        end if

        ces = c_error_state

      end function

      subroutine warn(msg)
!doc$ subroutine warn(msg)
        character(*), intent(in) :: msg
!       effects: Writes msg to the error file.

!cod$
        select case (file_status)
        case (EF_ON)
          write(error_unit,'(a)') trim(msg)
          call flush_i()
        end select

      end subroutine

      subroutine notify_msg_integer(msg,n)
!doc$ subroutine notify(msg,n)
        character(*), intent(in) :: msg
        integer, intent(in) :: n
!       modifies: Error file.
!       effects: Writes msg and n to error file.

!cod$
        select case (file_status)
        case (EF_ON)
          write(error_unit,'(a," = ",i0)') msg, n
          call flush_i()
        end select

      end subroutine

      subroutine notify_msg_real(msg,r)
!doc$ subroutine notify(msg,r)
        character(*), intent(in) :: msg
        real(double), intent(in) :: r
!       modifies: Error file.
!       effects: Writes msg and r to error file.

!cod$
        select case (file_status)
        case (EF_ON)
          write(error_unit,'(a," = ",f0.5)') msg, r
          call flush_i()
        end select

      end subroutine

      subroutine notify_msg_complex(msg,c)
!doc$ subroutine notify(msg,c)
        character(*), intent(in) :: msg
        complex(double), intent(in) :: c
!       modifies: Error file.
!       effects: Writes msg and c to error file.

!cod$
        select case (file_status)
        case (EF_ON)
          write(error_unit,'(a," = ",2f0.5)') msg, c
          call flush_i()
        end select

      end subroutine

      subroutine sync_kgroup_process_errors()
!doc$ subroutine sync_kgroup_process_errors()
!       requires: Call from a location where all kgroup processes are present.
!       modifies: k_error_state
!       effects: Sets the kgroup state = (.or. kgroup process error states).
!       errors: MPI problem.

!cod$
        integer :: mpi_stat

        call MPI_ALLREDUCE(kp_error_state,k_error_state,1,MPI_LOGICAL,MPI_LOR,kgroup_comm,mpi_stat)

        if (mpi_stat == 0) then
          if (k_error_state) then
            if (.not.kp_error_state) then
              select case (file_status)
              case (EF_ON)
                write(error_unit,'("ERROR: Process error communicated during kgroup synchronization")')
                call flush_i()
              end select
            end if
          end if
        else
          select case (file_status)
          case (EF_ON)
            write(error_unit,'("ERROR: mpi_stat /= 0 in error_mod::sync_kgroup_process_errors")')
            call flush_i()
          end select
          k_error_state = .true.
        end if

      end subroutine

      subroutine sync_sgroup_process_errors()
!doc$ subroutine sync_sgroup_process_errors()
!       requires: Call from a location where all sgroup processes are present.
!       modifies: s_error_state
!       effects: Sets the sgroup state = (.or. sgroup process error states).
!       errors: MPI problem.

!cod$
        integer :: mpi_stat

        call MPI_ALLREDUCE(sp_error_state,s_error_state,1,MPI_LOGICAL,MPI_LOR,sgroup_comm,mpi_stat)

        if (mpi_stat == 0) then
          if (s_error_state) then
            if (.not.sp_error_state) then
              select case (file_status)
              case (EF_ON)
                write(error_unit,'("ERROR: Process error communicated during sgroup synchronization")')
                call flush_i()
              end select
            end if
          end if
        else
          select case (file_status)
          case (EF_ON)
            write(error_unit,'("ERROR: mpi_stat /= 0 in error_mod::sync_sgroup_process_errors")')
            call flush_i()
          end select
          s_error_state = .true.
        end if

      end subroutine

      subroutine sync_config_process_errors()
!doc$ subroutine sync_config_process_errors()
!       requires: Call from a location where all config processes are present.
!       modifies: c_error_state
!       effects: Sets the config error state = (.or. config process error states).
!       errors: MPI problem.

!cod$
        integer :: mpi_stat

        call MPI_ALLREDUCE(cp_error_state,c_error_state,1,MPI_LOGICAL,MPI_LOR,config_comm,mpi_stat)
        if (mpi_stat == 0) then
          if (c_error_state) then
            if (.not.cp_error_state) then
              select case (file_status)
              case (EF_ON)
                write(error_unit,'("ERROR: Process error communicated during config synchronization")')
                call flush_i()
              end select
            end if
          end if
        else
          select case (file_status)
          case (EF_ON)
            write(error_unit,'("ERROR: mpi_stat /= 0 in error_mod::sync_config_process_errors")')
            call flush_i()
          end select
          c_error_state = .true.
        end if

      end subroutine

      subroutine sync_configuration_errors()
!doc$ subroutine sync_configuration_errors()
!       requires: Call from a location where all configurations are present.
!       modifies: c_error_state
!       effects: Sets the configuration error state = (.or. configuration error states).
!       errors: MPI problem.

!cod$
        logical :: w_error_state
        integer :: mpi_stat

        call MPI_ALLREDUCE(c_error_state,w_error_state,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_stat)
        if (mpi_stat == 0) then
          if (w_error_state) then
            if (.not.c_error_state) then
              select case (file_status)
              case (EF_ON)
                write(error_unit,'("ERROR: Configuration error communicated during world synchronization")')
                call flush_i()
              end select
              c_error_state = .true.
            end if
          end if
        else
          select case (file_status)
          case (EF_ON)
            write(error_unit,'("ERROR: mpi_stat /= 0 in error_mod::sync_configuration_errors")')
            call flush_i()
          end select
          c_error_state = .true.
        end if

      end subroutine

      subroutine flush_i()

        endfile(error_unit)
        backspace(error_unit)

      end subroutine

      end module

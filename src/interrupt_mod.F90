!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module interrupt_mod
!doc$ module interrupt_mod

!     Interrupt_mod encapsulates procedures for interrupting a run.

      use kind_mod
      use mpi_mod
      use path_mod
      use error_mod
      use io_mod
      use utils_mod

!cod$
      implicit none ; private

!doc$
      public :: user_abort
      public :: user_stop
      public :: interrupt
      public :: interrupt_stop

!cod$

      interface interrupt_stop
         module procedure interrupt_stop_
      end interface interrupt_stop

      contains

      function user_abort() result(ua)
!doc$ function user_abort() result(ua)
        logical :: ua
!       effects: Returns .true. iff file "stop_name" exists and contains the command "ABORT".

!cod$
        logical :: ex
        character(line_len) :: cmd
        integer :: ios
        type(file_obj) :: f

        call my(file(trim(stop_name)),f)

        ua = .false.

        if (i_access(f)) inquire(file=x_name(f),exist=ex)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ex)
        if (ex) then
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
          if (i_access(f)) read(x_unit(f),'(a)',iostat=ios) cmd
          if (i_comm(f)) call broadcast(FILE_SCOPE,cmd)
          select case (trim(adjustl(cmd)))
          case ("ABORT","Abort","abort","A","a")
            ua = .true.
          end select
          if (i_access(f)) close(unit=x_unit(f))
        end if

        call glean(thy(f))

      end function

      function user_stop() result(us)
!doc$ function user_stop() result(us)
        logical :: us
!       effects: Returns .true. iff file "stop_name" exists and contains the command "STOP".

!cod$
        logical :: ex
        character(line_len) :: cmd
        integer :: ios
        type(file_obj) :: f

        call my(file(trim(stop_name)),f)

        us = .false.

        if (i_access(f)) inquire(file=x_name(f),exist=ex)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ex)
        if (ex) then
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
          if (i_access(f)) read(x_unit(f),'(a)',iostat=ios) cmd
          if (i_comm(f)) call broadcast(FILE_SCOPE,cmd)
          select case (trim(adjustl(cmd)))
          case ("STOP","Stop","stop","S","s")
            us = .true.
          end select
          if (i_access(f)) close(unit=x_unit(f))
        end if

        call glean(thy(f))

      end function

      subroutine interrupt()
!doc$ subroutine interrupt()
!       effects: If present, removes the interrupt file.

!cod$
        logical :: ex
        integer :: ios

        if (mpi_isroot(WORLD)) then
          inquire(file=trim(stop_name),exist=ex)
          if (ex) then
            open(unit=INTERRUPT_UNIT,file=trim(stop_name),status='old',iostat=ios)
            if (ios /= 0) then
              call warn("WARNING: The interrupt file could not be opened")
              goto 100
            end if
            close(unit=INTERRUPT_UNIT,status='delete')
          end if
        end if

100     call barrier(WORLD)

      end subroutine

      !* Method to terminate execution due to the provided fatal error message

      subroutine interrupt_stop_( file , line , mesg )

         character(*) :: file, mesg
         integer :: line

         if ( mpi_isroot( world ) ) then
            write(*,'(/,"ERROR: ",a," (src/",a,":",a,")")') trimstr(mesg),basename(file),num2str(line)
         end if

         call mpi_stop()
         stop

      end subroutine interrupt_stop_

      end module

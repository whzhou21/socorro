!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module diary_mod
!doc$ module diary_mod

!     diary_mod encapsulates a special file_obj (diary_file) and provides routines for writing
!     character strings to the diary_file.

      use error_mod
      use io_mod
      use kind_mod
      use path_mod
      use mpi_mod
      use utils_mod
      use version_mod

!cod$
      implicit none ; private

! Version numbers

      character(18), parameter :: CHECK_SYMMETRY_RELEASE = "check_symmetry 1.0"
      character(17), parameter :: CHECK_KPOINTS_RELEASE = "check_kpoints 1.0"

      type(file_obj) :: diary_file

!doc$
      public :: diary_start
      public :: diary_stop
      public :: diaryfile
      public :: diary
      public :: diary_socorro_env
      public :: diary_check_symmetry_env
      public :: diary_check_kpoints_env

!cod$
      interface diary
        module procedure diary_message
      end interface

      contains

      subroutine diary_start()
!doc$ subroutine diary_start()
!       effects: Starts the diary system.
!       errors: Problems opening a diary file. Passes errors.

!cod$
        integer :: ios

        call my(file(trim(diary_path),first_only=.true.),diary_file) ; if (error()) goto 100

        if (i_access(diary_file)) open(x_unit(diary_file),file=x_name(diary_file),status='unknown',iostat=ios)
        if (i_comm(diary_file)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: Problems opening the diary file")) goto 100

100     if (error("Exit diary_mod::diary_start")) continue

      end subroutine

      subroutine diary_stop()
!doc$ subroutine diary_stop()
!       effects: Stops the diary system.

!cod$
        logical :: op

        if (i_access(diary_file)) then
          inquire(x_unit(diary_file),opened=op)
          if (op) close(x_unit(diary_file))
        end if
        call glean(thy(diary_file))

      end subroutine

      function diaryfile() result(f)
!doc$ function diaryfile() result(f)
        type(file_obj) :: f
!       effects: Returns the diary file.

!cod$
        call my(diary_file,f)
        call bequeath(thy(f))
      end function

      subroutine diary_message(msg)
!doc$ subroutine diary(msg)
        character(*), intent(in), optional :: msg
!       requires: diary_file be open.
!       effects: Writes trim(msg) or trim(iobuf) to the diary file.

!cod$
        if (present(msg)) then
          if (i_access(diary_file)) write(x_unit(diary_file),'(a)') trim(msg)
        else
          if (i_access(diary_file)) write(x_unit(diary_file),'(a)') trim(iobuf)
        end if

        if (i_access(diary_file)) call flushbuf(diary_file)

      end subroutine

      subroutine diary_socorro_env()
!doc$ subroutine diary_socorro_env()
!        requires: diary_file be open.
!        effects: Writes the socorro environment to the diary_file and error files.

!cod$
         integer :: nc, ncp, nkg, nkgp, np, nsg, nsgp, nt, un
         character(line_len) :: date, time
         character(:), allocatable :: dd, mo, yy, hh, mm, ss

         un = x_unit(diary_file)

         call date_and_time(date,time)

         dd = trimstr(date(7:8))
         mo = trimstr(date(5:6))
         yy = trimstr(date(1:4))
         hh = trimstr(time(1:2))
         mm = trimstr(time(3:4))
         ss = trimstr(time(5:6))

         np = mpi_nprocs(WORLD)
         nt = x_nthreads()

         nc = mpi_nconfigs()
         ncp = mpi_nprocs(CONFIG)

         nsg = mpi_nsgroups()
         nsgp = mpi_nprocs(SGROUP)

         nkg = mpi_nkgroups()
         nkgp = mpi_nprocs(KGROUP)

         if ( i_access(diary_file) ) then

            write(un,'("Socorro ",a)') trimstr(x_version())
            write(un,'(/,t4,"Started on ",a," ",a," ",a," at ",a,":",a,":",a)') dd,mo,yy,hh,mm,ss
            write(un,'(/,"Runtime environment:")')

            if ( nc == 1 ) then
               write(un,'(/,t4,i0," configuration is running")') nc
               if (ncp == 1) then
                  write(un,'(/,t4,i0," MPI process is working on this configuration")') ncp
               else
                  write(un,'(/,t4,i0," MPI processes are working on this configuration")') ncp
               end if
            else
               write(un,'(/,t4,i0," configurations are running")') nc
               write(un,'(t4,i0," MPI processes are working on the configurations")') np
               if ( ncp == 1 ) then
                  write(un,'(t4,i0," MPI process is working on this configuration")') ncp
               else
                  write(un,'(t4,i0," MPI processes are working on this configuration")') ncp
               end if
            end if

            if ( nt == 0 ) then
               continue
            else if ( nt == 1 ) then
               write(un,'(t4,a," OMP thread is running per MPI process")') num2str(nt)
            else
               write(un,'(t4,a, " OMP threads are running per MPI process")') num2str(nt)
            end if

            if ( nsg == 1 ) then
               write(un,'(/,t4,a," spin group is running")') num2str(nsg)
               if ( nsgp == 1 ) then
                  write(un,'(t4,a," MPI process is working on this spin group")') num2str(nsgp)
               else
                  write(un,'(t4,a," MPI processes are working on this spin group")') num2str(nsgp)
               end if
            else
               write(un,'(/,t4,a," spin groups are running")') num2str(nsg)
               if ( nsgp == 1 ) then
                  write(un,'(t4,a," MPI process is working on each spin group")') num2str(nsgp)
               else
                  write(un,'(t4,a," MPI processes are working on each spin group")') num2str(nsgp)
               end if
            end if

            if ( nkg == 1 ) then
               write(un,'(/,t4,a," k-point group is running")') num2str(nkg)
               if ( nkgp == 1 ) then
                  write(un,'(t4,a," MPI process is working on this k-point group")') num2str(nkgp)
               else
                  write(un,'(t4,a," MPI processes are working on this k-point group")') num2str(nkgp)
               end if
            else
               write(un,'(/,t4,a," k-point groups are running")') num2str(nkg)
               if ( nkgp == 1 ) then
                  write(un,'(t4,a," MPI process is working on each k-point group")') num2str(nkgp)
               else
                  write(un,'(t4,a," MPI processes are working on each k-point group")') num2str(nkgp)
               end if
            end if

            call flush(un)

         end if

      end subroutine

      subroutine diary_check_symmetry_env()
!doc$ subroutine diary_check_symmetry_env()
!       requires: diary_file be open.
!       effects: Writes the check_symmetry environment to diaryfile.

!cod$
        integer :: np

        np = mpi_nprocs(WORLD)

        call diary(CHECK_SYMMETRY_RELEASE)
        if (i_access(diary_file)) then
          write(x_unit(diary_file),'(/,"Run-time environment:")')
          if (np == 1) then
            write(x_unit(diary_file),'(t4,i0," processor is working on this configuration")') np
          else
            write(x_unit(diary_file),'(t4,i0," processors are working on this configuration")') np
          end if
          call flushbuf(diary_file)
        end if

      end subroutine

      subroutine diary_check_kpoints_env()
!doc$ subroutine diary_check_kpoints_env()
!       requires: diary_file be open.
!       effects: Writes the check_kpoints environment to diaryfile.

!cod$
        integer :: np

        np = mpi_nprocs(WORLD)

        call diary(CHECK_KPOINTS_RELEASE)
        if (i_access(diary_file)) then
          write(x_unit(diary_file),'(/,"Run-time environment:")')
          if (np == 1) then
            write(x_unit(diary_file),'(t4,i0," processor is working on this configuration")') np
          else
            write(x_unit(diary_file),'(t4,i0," processors are working on this configuration")') np
          end if
          call flushbuf(diary_file)
        end if

      end subroutine

      end module

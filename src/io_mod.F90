!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module io_mod
!doc$ module io_mod

!     One datatype is available here: type(file_obj)

!     io_mod: - provides routines for constructing and accessing a file_obj.
!             - provides access to the standard input (input) and output (output) files.
!             - provides a buffer (iobuf) that can be used to pass strings to the error, diary, and output files.

!     Note: This module should not be changed without consulting AFW.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use,intrinsic :: iso_fortran_env,only : input_unit, output_unit

!cod$
      implicit none
      private

      logical, parameter :: default_mangle = .true.
      logical, parameter :: default_first_only = .true.

      logical, dimension(:), allocatable :: unit_table  ! .false. means the unit is available for use

      character(10*line_len) :: iobuf

      integer, parameter, public :: FILE_SCOPE = CONFIG
      integer, public :: INTERRUPT_UNIT
      integer, public :: IOSTAT_OK
      integer, public :: IOSTAT_EOL
      integer, public :: IOSTAT_EOF

      type :: file_rep
        integer :: ref
        integer :: unit
        logical :: first_only
        character(line_len) :: name
        character(line_len) :: manglename
      end type

      type, public :: file_obj
        private
        integer :: ref
        type(file_rep), pointer :: o
      end type

      type(file_obj) :: input, output

!doc$
      public :: io_start
      public :: io_stop
      public :: file
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: min_unit
      public :: max_unit
      public :: unit_in_use
      public :: x_unit
      public :: x_name
      public :: i_access
      public :: i_comm
      public :: display
      public :: input
      public :: output
      public :: iobuf
      public :: flushbuf

!cod$
      interface file
         module procedure constructor_file
      end interface
      interface my
         module procedure my_file, my_new_file
      end interface
      interface thy
         module procedure thy_file
      end interface
      interface glean
         module procedure glean_file
      end interface
      interface bequeath
         module procedure bequeath_file
      end interface
      interface assignment(=)
        module procedure assign_file
      end interface
      interface x_unit
         module procedure f_unit
      end interface
      interface x_name
         module procedure f_name
      end interface
      interface i_access
        module procedure i_access_file
      end interface
      interface i_comm
        module procedure i_comm_file
      end interface
      interface flushbuf
         module procedure f_flushbuf
      end interface

      contains

      subroutine io_start()
!doc$ subroutine io_start()
!       effects: Initializes the I/O system.

!cod$
        logical :: c_ex, p_ex, c_op, p_op, in_use
        character :: c
        integer :: iu, ios

        ! Construct the unit table
        allocate( unit_table(first_unit:last_unit) )
        do iu = first_unit,last_unit
          unit_table(iu) = .false.
          inquire(iu,exist=p_ex)
          call allreduce(WORLD,MPI_LOR,p_ex,c_ex)
          if (.not.c_ex) unit_table(iu) = .true.
          inquire(iu,opened=p_op)
          call allreduce(WORLD,MPI_LOR,p_op,c_op)
          if (c_op) unit_table(iu) = .true.
        end do

        ! Construct input and output file objects
        call my(file(name="input",first_only=.true.,unit=input_unit),input) ; if (error()) goto 200
        call my(file(name="output",first_only=.true.,unit=output_unit),output)

        ! Reserve an interrupt unit
        in_use = .true.
        do iu = first_unit,last_unit
          if (.not.unit_table(iu)) then
            unit_table(iu) = .true.
            in_use = .false.
            INTERRUPT_UNIT = iu
            exit
          end if
        end do
        if (error(in_use,"ERROR: no unit was found")) goto 200

        ! Determine iostat codes
        in_use = .true.
        do iu = first_unit,last_unit
          if (.not.unit_table(iu)) then
            unit_table(iu) = .true.
            in_use = .false.
            exit
          end if
        end do
        if (error(in_use,"ERROR: no unit was found")) goto 200
        if (mpi_isroot(WORLD)) then
          open(unit=iu,file="io_scratch",status="new",iostat=ios)
          if (error(ios /= 0,"ERROR: unable to open file")) goto 100
          write(iu,*) "this is a test"
          write(iu,*) "second line"
          rewind(iu)
          read(iu,'(a1)',advance='no',iostat=ios) c
          IOSTAT_OK = ios
          do while (ios == IOSTAT_OK)
            read(iu,'(a1)',advance='no',iostat=ios) c
          end do
          IOSTAT_EOL = ios
          do while ((ios == IOSTAT_OK) .or. (ios == IOSTAT_EOL))
            read(iu,'(a1)',advance='no',iostat=ios) c
          end do
          IOSTAT_EOF = ios
          close(iu,iostat=ios,status="delete")
        end if
        unit_table(iu) = .false.

100     call sync_configuration_errors() ; if (error()) goto 200

        call broadcast(WORLD,IOSTAT_OK)
        call broadcast(WORLD,IOSTAT_EOL)
        call broadcast(WORLD,IOSTAT_EOF)

200     if (error("Exit io_mod::io_start")) continue

      end subroutine

      subroutine io_stop()
!doc$ subroutine io_stop()
!       effects: Stops the I/O system.

!cod$
        call glean(thy(input))
        call glean(thy(output))
        if (allocated( unit_table )) deallocate( unit_table )
      end subroutine

      function constructor_file(name,first_only,unit) result(f)
!doc$ function file(name,first_only,unit) result(f)
        character(*), intent(in) :: name
        logical, intent(in), optional :: first_only
        integer, intent(in), optional :: unit
        type(file_obj) :: f
!       effects: Constructs a new f.
!       errors: If unit is present: Unit out of range. Unit not available. Open file associated with unit.
!               If unit is not present: Open file associated with the selected unit. No available units.

!cod$
        logical :: c_op, p_op
        integer :: un

        if (present(unit)) then
          un = unit
          if (error((un < first_unit) .or. (un > last_unit),"ERROR: requested unit is out of range")) goto 100
          if ((un /= input_unit) .and. (un /= output_unit))  then
            if (error(unit_table(un),"ERROR: requested unit is not available")) goto 100
            inquire(un,opened=p_op)
            call allreduce(FILE_SCOPE,MPI_LOR,p_op,c_op)
            if (error(c_op,"ERROR: an open file is associated with the requested unit")) goto 100
          end if
        else
          un = first_unit
          do
            if (.not.unit_table(un)) then
              inquire(un,opened=p_op)
              call allreduce(FILE_SCOPE,MPI_LOR,p_op,c_op)
              if (error(c_op,"ERROR: an open file is associated with an available unit")) goto 100
              exit
            end if
            un = un + 1
            if (error(un > last_unit,"ERROR: no units are available")) goto 100
          end do
        end if

        unit_table(un) = .true.

        f%ref = 0
        allocate( f%o )
        f%o%ref = 0
        f%o%unit = un
        if (present(first_only)) then
          f%o%first_only = first_only
        else
          f%o%first_only = default_first_only
        end if
        f%o%name = name
        write(f%o%manglename,'(a,"_",I3.3)') trim(name), mpi_myproc(FILE_SCOPE)

100     if (error("Exit io_mod::constructor_file")) continue

      end function

      subroutine my_file(f)
!doc$ subroutine my(f)
        type(file_obj) :: f

!cod$
        f%ref = f%ref + 1
        f%o%ref = f%o%ref + 1
      end subroutine

      subroutine my_new_file(fi,f)
!doc$ subroutine my(fi,f)
        type(file_obj) :: fi, f

!cod$
        f%ref = 1
        f%o => fi%o
        f%o%ref = f%o%ref + 1
      end subroutine

      function thy_file(f) result(fo)
!doc$ function thy(f) result(fo)
        type(file_obj) :: f, fo

!cod$
        f%ref = f%ref - 1
        f%o%ref = f%o%ref - 1
        fo%ref = f%ref
        fo%o => f%o
      end function

      subroutine glean_file(f)
!doc$ subroutine glean(f)
        type(file_obj) :: f
!       requires: f be closed.

!cod$
        if (f%o%ref < 1) then
          unit_table(f%o%unit) = .false.
          deallocate( f%o )
        end if
      end subroutine

      subroutine bequeath_file(f)
!doc$ subroutine bequeath(f)
        type(file_obj) :: f

!cod$
        continue
      end subroutine

      subroutine assign_file(f,f2)
!doc$ subroutine assignment(=)(f,f2)
        type(file_obj), intent(inout) :: f
        type(file_obj), intent(in) :: f2

!cod$
        type(file_obj) :: ft
        call my(f2)
        ft%o => f%o
        f%o%ref = f%o%ref - f%ref
        f%o => f2%o
        f%o%ref = f%o%ref + f%ref
        call glean(ft)
        call glean(thy(f2))
      end subroutine

      function min_unit() result(u)
!doc$ function min_unit() result(u)
        integer :: u
!       effects: Returns the lowest allowed file unit.

!cod$
        u = first_unit
      end function

      function max_unit() result(u)
!doc$ function max_unit() result(u)
        integer :: u
!       effects: Returns the largest allowed file unit.

!cod$
        u = last_unit
      end function

      function unit_in_use(u) result(l)
!doc$ function unit_in_use(u) result(l)
        integer, intent(in) :: u
        logical :: l
!       effects: Returns the use status of unit u.
!       errors: u out of range.

!cod$
        if (error((u < first_unit) .or. (u > last_unit),"ERROR: u is out of range")) goto 100
        l = unit_table(u)
100     if (error("Exit io_mod::unit_in_use")) continue
      end function

      function f_unit(f) result(u)
!doc$ function x_unit(f) result(u)
        type(file_obj) :: f
        integer :: u
!       effects: Returns the unit number of f.

!cod$
        call my(f)
        u = f%o%unit
        call glean(thy(f))
      end function

      function f_name(f,mangle) result(n)
!doc$ function x_name(f,mangle) result(n)
        type(file_obj) :: f
        logical, intent(in), optional :: mangle
        character(line_len) :: n
!       effects: Returns the name of f.

!cod$
        logical :: local_mangle
        local_mangle = default_mangle
        call my(f)
        if (present(mangle)) local_mangle = mangle
        if (f%o%first_only .or. (.not.local_mangle)) then 
          n = f%o%name
        else
          n = f%o%manglename
        end if
        call glean(thy(f))
      end function 

      function i_access_file(f) result(l)
!doc$ function i_access(f) result(l)
        type(file_obj) :: f
        logical :: l
!       effects: Returns .true. if access to f is required.

!cod$
        call my(f)
        l = (.not.f%o%first_only .or. mpi_isroot(FILE_SCOPE))
        call glean(thy(f))
      end function

      function i_comm_file(f) result(l)
!doc$ function i_comm(f) result(l)
        type(file_obj) :: f
        logical :: l
!       effects: Returns .true. if communication with f is required.

!cod$
        call my(f)
        l = f%o%first_only
        call glean(thy(f))
      end function

      subroutine display(msg)
!doc$ subroutine display(msg)
        character(*), intent(in), optional :: msg
!       effects: Writes trim(msg) or trim(iobuf) to standard output.

!cod$
        if (present(msg)) then
          if (i_access(output)) write(x_unit(output),'(a)') trim(msg)
        else
          if (i_access(output)) write(x_unit(output),'(a)') trim(iobuf)
        end if
      end subroutine

      subroutine f_flushbuf(f)
!doc$ subroutine flushbuf(f)
        type(file_obj) :: f
!       effects: flushes the buffer by closing and reopening the file

!cod$
        call my(f)
        if (mpi_isroot(FILE_SCOPE)) then
          close(unit=f%o%unit)
          open(unit=f%o%unit,file=trim(f%o%name),status='old',position='append')
        end if
        call glean(thy(f))
      end subroutine

      end module

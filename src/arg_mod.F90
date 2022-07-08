!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module arg_mod
!doc$ module arg_mod

!     arg module provides routines for reading and distributing control parameters. The control file is
!     read once in routine arg_start and the results are cached for later retrieval. Note that unit 4 is
!     associated with the control file in arg_start and then dissociated once the contents have been read
!     and cached.

      use error_mod
      use interrupt_mod
      use kind_mod
      use mpi_mod
      use path_mod
      use utils_mod

!cod$
      implicit none ; private

      type :: arg_param
        private
        character(line_len) :: var        ! keyword
        character(line_len) :: val        ! value associated with the keyword
      end type
      type(arg_param), dimension(:), allocatable :: arg_params

!doc$
      public :: arg_start
      public :: arg_stop
      public :: arg
      public :: arglc

!cod$
      interface arg
        module procedure read_arg_log, &
                         read_arg_ch, &
                         read_arg_int, read_arg_int_1d, &
                         read_arg_dpr, read_arg_dpr_1d, &
                         read_arg_dpc
      end interface
      interface arglc
        module procedure read_arglc
      end interface

      contains

! public routines

      subroutine arg_start()
!doc$ subroutine arg_start()
!       effects: Starts the argument retrieval system.
!       errors: Problems locating or opening the arguments file.

!cod$
        logical :: found
        character(line_len) :: buffer, value, infile
        integer :: ii, nlines, iostatus, inflag, iarg, argc, argunit

        ! Process the command-line arguments

        inflag = 0

        if ( mpi_isroot( world ) ) argc = command_argument_count()
        call broadcast_seh(argc)
        if ( argc == 0 ) call interrupt_stop(FLERR,"No command-line arguments given")

        iarg = 1
        do while ( iarg < argc + 1 )
           call get_command_argument(iarg,value)
           if ( trimstr(value) == "-in" .or. trimstr(value) == "-i" ) then
              if ( iarg + 1 > argc ) call interrupt_stop(FLERR,"Invalid command-line argument")
              inflag = iarg + 1
              iarg = iarg + 2
           else
              iarg = iarg + 1
           end if
        end do

        ! Check the status of the -in command-line switch

        if ( inflag == 0 ) call interrupt_stop(FLERR,"The '-in' command-line switch was not found")

        ! Process the input-file arguments

        if ( mpi_isroot( world ) ) then
           call get_command_argument(inflag,infile)
           inquire(file=trim(infile),exist=found)
        end if
        call broadcast_seh(found)
        if ( .not.found ) call interrupt_stop(FLERR,"The input file '"//trimstr(infile)//"' was not found")

        if ( mpi_isroot( world ) ) open(newunit=argunit,file=trim(infile),status="unknown",iostat=iostatus)
        call broadcast_seh(iostatus)
        if ( iostatus /= 0 ) call interrupt_stop(FLERR,"The input file '"//trimstr(infile)//"' could not be opened")

        ! Determine the number of lines in the control file

        if ( mpi_isroot( world ) ) then
           nlines = -1
           iostatus = 0
           do while ( iostatus == 0 )
              nlines = nlines + 1
              read(argunit,'(a)',iostat=iostatus) buffer
           end do
           rewind( argunit )
        end if

        call broadcast_seh(nlines)
        allocate( arg_params( nlines ) )

        ! Read and broadcast the contents of the control file

        do ii = 1,nlines
           if ( mpi_isroot( world ) ) read(argunit,'(a)') buffer
           call broadcast_seh(buffer)
           call arg_split_i(buffer,arg_params(ii))
        end do

        if ( mpi_isroot( world ) ) close( argunit )

      end subroutine

      subroutine arg_stop()
!doc$ subroutine arg_stop()
!       effects: Stops the argument retrieval system.

!cod$
        if (allocated( arg_params )) deallocate( arg_params )
      end subroutine

      subroutine read_arg_log(name,value,found)
!doc$ subroutine arg(name,value,found)
        character(*), intent(in) :: name
        logical, intent(out):: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        integer :: ios = 0
        var = name
        call read_arg_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        read(buf,*,iostat=ios) value
        if (error(ios /= 0,"ERROR: "//buf(1:len_trim(buf))//" to logical" )) goto 100
100     if (error("Exit arg_mod::read_arg_log")) continue
      end subroutine

      subroutine read_arg_ch(name,value,found)
!doc$ subroutine arg(name,value,found)
        character(*), intent(in) :: name
        character(*), intent(out) :: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        var = name
        call read_arg_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        value = buf
100     if (error("Exit arg_mod::read_arg_ch")) continue
      end subroutine

      subroutine read_arg_int(name,value,found)
!doc$ subroutine arg(name,value,found)
        character(*), intent(in) :: name
        integer, intent(out) :: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        integer :: ios = 0
        var = name
        call read_arg_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        read(buf,*,iostat=ios) value
        if (error(ios /= 0,"ERROR: "//buf(1:len_trim(buf))//" to integer" )) goto 100
100     if (error("Exit arg_mod::read_arg_int")) continue
      end subroutine

      subroutine read_arg_int_1d(name,value,found)
!doc$ subroutine arg(name,value,found)
        character(*), intent(in) :: name
        integer, dimension(:), intent(out) :: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        integer :: ios = 0
        var = name
        call read_arg_i(var,buf,fnd) 
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        read(buf,*,iostat=ios) value
        if (error(ios /= 0,"ERROR: "//buf(1:len_trim(buf))//" to integer(:)" )) goto 100
100     if (error("Exit arg_mod::read_arg_int_1d")) continue
      end subroutine

      subroutine read_arg_dpr(name,value,found)
!doc$ subroutine arg(name,value,found)
        character(*), intent(in) :: name
        real(double), intent(out):: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        integer :: ios = 0
        var = name
        call read_arg_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        read(buf,*,iostat=ios) value
        if (error(ios /= 0,"ERROR: "//buf(1:len_trim(buf))//" to real(double)" )) goto 100
100     if (error("Exit arg_mod::read_arg_dpr")) continue
      end subroutine

      subroutine read_arg_dpr_1d(name,value,found)
!doc$ subroutine arg(name,value,found)
      character(*), intent(in) :: name
      real(double), dimension(:), intent(out) :: value
      logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        integer :: ios = 0
        var = name
        call read_arg_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        read(buf,*,iostat=ios) value
        if (error(ios /= 0,"ERROR: "//buf(1:len_trim(buf))//" to real(double)(:)" )) goto 100
100     if (error("Exit arg_mod::read_arg_dpr_1d")) continue
      end subroutine

      subroutine read_arg_dpc(name,value,found)
!doc$ subroutine arg(name,value,found)
        character(*), intent(in) :: name
        complex(double), intent(out) :: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        integer :: ios = 0
        var = name
        call read_arg_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        read(buf,*,iostat=ios) value
        if (error(ios /= 0,"ERROR: "//buf(1:len_trim(buf))//" to complex" )) goto 100
100     if (error("Exit arg_mod::read_arg_dpc")) continue
      end subroutine

      subroutine read_arglc(name,value,found)
!doc$ subroutine arglc(name,value,found)
        character(*), intent(in) :: name
        character(*), intent(out) :: value
        logical, intent(out), optional :: found

!cod$
        character(line_len) :: var, buf
        logical :: fnd
        var = name
        call read_arglc_i(var,buf,fnd)
        if (present(found)) then
          found = fnd
          if (.not.fnd) goto 100
        else
          if (error(.not.fnd,"ERROR: "//name//" not found")) goto 100
        end if
        value = buf
100     if (error("Exit arg_mod::read_arglc")) continue
      end subroutine

! private routines

      subroutine arg_split_i(s,p)
        character(line_len), intent(inout) :: s
        type(arg_param), intent(out) :: p
        integer :: k
        call trim_comments_i(s)
        call trim_whitespace_i(s)
        do k = 1,line_len
          if (s(k:k) == ' ') exit
        end do
        p%var = s(1:k-1); p%val=s(k+1:line_len)
        call upcase_arg_i(p%var)
        call trim_whitespace_i(p%val)
      end subroutine

      subroutine trim_comments_i(s)
        character(line_len) :: s
        integer :: k
        do k = 1,line_len
          if (s(k:k) == '!') exit
        end do
        s = s(1:k-1)
      end subroutine

      subroutine trim_whitespace_i(s)
        character(line_len) :: s
        integer :: k
        do k = 1,line_len
          if (s(k:k) == ' ') cycle
          if (s(k:k) == '	') cycle
          exit
        end do
        s = s(k:line_len)
      end subroutine

      subroutine read_arg_i(var,val,found)
        character(line_len), intent(inout)  :: var
        character(line_len), intent(out)  :: val
        logical, intent(out) :: found
        integer :: i, j
        found = .false.
        call trim_whitespace_i(var)
        call upcase_arg_i(var)
        j = len_trim(var)
        do i = 1,size(arg_params)
          if (var(1:j)==arg_params(i)%var(1:j)) then
            found=.true.
            val = arg_params(i)%val
            exit
          end if
        end do
      end subroutine

      subroutine read_arglc_i(var,val,found)
        character(line_len), intent(inout)  :: var
        character(line_len), intent(out)  :: val
        logical, intent(out) :: found
        integer :: i, j
        found = .false.
        val = ''
        call trim_whitespace_i(var)
        call upcase_arg_i(var)
        j = len_trim(var)
        do i = 1,size(arg_params)
          if (var(1:j)==arg_params(i)%var(1:j)) then
            found=.true.
            val = arg_params(i)%val
            exit
          end if
        end do
        call locase_arg_i(val)
      end subroutine

      subroutine upcase_arg_i(s)
        character(line_len) :: s
        integer i, lstr, c
        integer, parameter :: a = iachar('a'), z = iachar('z'), shift = (iachar('A') - iachar('a'))
        lstr = len_trim(s)
        do i = 1,lstr
          c = iachar(s(i:i))
          if ((c >= a) .and. (c <= z)) c = c + shift
          s(i:i) = achar(c)
        end do
      end subroutine

      subroutine locase_arg_i(s)
        character(line_len) :: s
        integer i, lstr, c
        integer, parameter :: a = iachar('A'), z = iachar('Z'), shift = (iachar('A') - iachar('a'))
        lstr = len_trim(s)
        do i = 1,lstr
          c = iachar(s(i:i))
          if ((c >= a) .and. (c <= z)) c = c - shift
          s(i:i) = achar(c)
        end do
      end subroutine

      end module

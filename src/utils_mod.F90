!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module utils_mod
!doc$ module utils_mod

!     Comments...

      use kind_mod

!cod$
      implicit none ; private

!doc$
      public :: basename
      public :: num2str
      public :: time_and_date
      public :: trimstr

!cod$
      interface basename
         module procedure basename_
      end interface

      interface num2str
         module procedure int2str_ , real2str_
      end interface

      interface time_and_date
         module procedure time_and_date_
      end interface

      interface trimstr
         module procedure trimstr_
      end interface

      contains

! *** Public routines

      function basename_( path ) result( fname )
!doc$ function socorro()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         character(*) :: path
         character(:), allocatable :: fname

         character(line_len) :: str
         integer :: i, j
         logical :: back = .true.

         i = index(path,"/",back) + 1
         j = len_trim(path)

         write(str,'(a)') path(i:j)
         fname = trimstr(str)

      end function basename_

      function int2str_( val ) result( str )
!doc$ function socorro()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         integer :: val
         character(:), allocatable :: str

         character(line_len) :: tmp

         write(tmp,'(i0)') val
         str = trimstr(tmp)

      end function int2str_

      function real2str_( val , nd ) result( str )
!doc$ function socorro()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         real(double) :: val
         integer, optional :: nd
         character(:), allocatable :: str

         character(line_len) :: tmp, fmt

         if (present(nd)) then
            write(fmt,'("f128.",i0)') nd
         else
            write(fmt,'("f128.",i0)') 5
         end if

         write(tmp,'('//trimstr(fmt)//')') val
         str = trimstr(tmp)

      end function real2str_

      function time_and_date_() result( r )
!doc$ function time_and_date()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         character(line_len) :: date, time, tmp
         character(:), allocatable :: d, n, y, h, m, s, month, r

         call date_and_time(date,time)

         d = trimstr(date(7:8))
         n = trimstr(date(5:6))
         y = trimstr(date(1:4))

         h = trimstr(time(1:2))
         m = trimstr(time(3:4))
         s = trimstr(time(5:6))

         select case (trimstr(n))
         case ( "01" )
            month = "Jan"
         case ( "02" )
            month = "Feb"
         case ( "03" )
            month = "Mar"
         case ( "04" )
            month = "Apr"
         case ( "05" )
            month = "May"
         case ( "06" )
            month = "June"
         case ( "07" )
            month = "July"
         case ( "08" )
            month = "Aug"
         case ( "09" )
            month = "Sep"
         case ( "10" )
            month = "Oct"
         case ( "11" )
            month = "Nov"
         case ( "12" )
            month = "Dec"
         end select

         write(tmp,'(a," ",a," ",a," at ",a,":",a,":",a)') d,month,y,h,m,s
         r = trimstr(tmp)

      end function time_and_date_

      function trimstr_( strin ) result( strout )
!doc$ function socorro()
!        effects:
!        errors:
!        requires:
!        notes:

!cod$
         character(*) :: strin
         character(:), allocatable :: strout

         strout = trim(adjustl(strin))

      end function trimstr_

      end module utils_mod

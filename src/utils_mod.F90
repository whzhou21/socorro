!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!                                                                                                                                  !
!  Copyright (2020). See the README file in the top-level directory.                                                               !
!  This software is distributed with the GNU General Public License.                                                               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

module utils_mod

   use kind_mod

   implicit none ; private

   !* Publicly available parameters and procedures

   public :: basename
   public :: num2str
   public :: trimstr

   !* Interfaces for the publicly available procedures

   interface basename
      module procedure basename_
   end interface basename

   interface num2str
      module procedure int2str_ , real2str_
   end interface num2str

   interface trimstr
      module procedure trimstr_
   end interface trimstr

contains


   !* Method to retrieve the base filename from a path

   function basename_( path ) result( fname )

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


   !* Method to convert an integer to a string

   function int2str_( val ) result( str )

      integer :: val
      character(:), allocatable :: str

      character(line_len) :: tmp

      write(tmp,'(i0)') val
      str = trimstr(tmp)

   end function int2str_


   !* Method to convert a real to a string

   function real2str_( val , nd ) result( str )

      real(double) :: val
      integer, optional :: nd
      character(:), allocatable :: str

      character(line_len) :: tmp, fmt

      if ( present( nd ) ) then
         write(fmt,'("f128.",i0)') nd
      else
         write(fmt,'("f128.",i0)') 5
      end if

      write(tmp,'('//trimstr(fmt)//')') val
      str = trimstr(tmp)

   end function real2str_


   !* Method to remove any leading and trailing whitespace from a string

   function trimstr_( str1 ) result ( str2 )

      character(*) :: str1
      character(:), allocatable :: str2

      str2 = trim(adjustl(str1))

   end function trimstr_


end module utils_mod

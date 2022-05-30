! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module intstack_mod
!doc$ module INTSTACK_MOD

      use kind_mod
      implicit none
      private

      type, public :: intstack_obj
         private
         integer(LONGLONG) :: n
         type (intstack_obj), pointer :: next
      end type intstack_obj

      public :: createstack, deletestack, push, pop

      interface createstack
         module procedure createintstack
      end interface

      interface push
         module procedure push_int
      end interface

      interface pop
         module procedure pop_int
      end interface

      interface deletestack
         module procedure deleteintstack
      end interface

!*********************************************************
      contains
!*********************************************************



      subroutine createintstack(head)
        type (intstack_obj), pointer :: head

        nullify(head)

        return
      end subroutine

!*********************************************************

      subroutine deleteintstack(head)
        type (intstack_obj), pointer :: head

        logical :: Ok
        integer(LONGLONG) :: i

        Ok = .TRUE.
        do while (ok)
           i = pop(head, Ok)
!write(*,*) 'deleteintstack: i=',i, ' Ok=',Ok
        end do 

        return
      end subroutine

!*********************************************************

      integer(LONGLONG) function pop_int(head, Ok)
        type (intstack_obj), pointer :: head
        logical,         intent(out) :: Ok
 
        type (intstack_obj), pointer :: ptr

        if (associated(head)) then
           ok = .TRUE.
           pop_int = head%n
           ptr => head
           head => head%next
           deallocate(ptr)
        else
           pop_int = -1
           ok = .FALSE.
        end if

        return
      end function pop_int

!*********************************************************

      subroutine push_int(head, value)
        type (intstack_obj),  pointer :: head
        integer(LONGLONG), intent(in) :: value
 
        type (intstack_obj), pointer :: ptr

        allocate(ptr)

        ptr%n = value
        ptr%next => head

        head => ptr;

        return
      end subroutine

end module

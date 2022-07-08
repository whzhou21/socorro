!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module dyad_kpoint_mod
!doc$ module dyad_kpoint_mod

!     Defines a type(dyad_kpoint_obj), which encapsulates a k-point specific "dyad": an operator parameterized by a pair of
!     multivectors, V and W, whose action on an arbitrary multivector X is defined by multiplication by the symmetrized outer
!     product of V and W.  If multivectors are considered to be matrices with x_n_bands(X) columns and x_n_gvectors(X) rows,
!     then the action of the dyad_kpoint_obj defined by V and W on X is given by (VW' + WV')X, where "'" indicates the
!     transpose of the preceding matrix.  Dyads may either contain data or be empty.  Calling the apply() function on an
!     empty dyad_kpoint is an error.
!

      use kind_mod
      use error_mod
      use multivector_mod
      use ghost_mod

!cod$
      implicit none
      private

      type :: dyad_kpoint_rep
        integer :: ref
        type(ghost) :: g
        type(multivector_obj), pointer :: dv
        type(multivector_obj), pointer :: dw
      end type

      type, public :: dyad_kpoint_obj
        private
        integer :: ref
        type(dyad_kpoint_rep), pointer :: o
      end type

!doc$
      public :: dyad_kpoint, update
      public :: my, thy, glean, bequeath, assignment(=)
      public :: is_empty, x_v, x_w, apply, x_ghost

!cod$

      interface dyad_kpoint
        module procedure constructor_dk
      end interface
      interface update
        module procedure update_dk
      end interface

      interface my
        module procedure my_dk, my_new_dk
      end interface
      interface thy
        module procedure thy_dk
      end interface
      interface glean
        module procedure glean_dk
      end interface
      interface bequeath
        module procedure bequeath_dk
      end interface
      interface assignment(=)
        module procedure assign_dk
      end interface

      interface is_empty
        module procedure is_empty_dk
      end interface
      interface x_v
        module procedure dk_v
      end interface
      interface x_w
        module procedure dk_w
      end interface
      interface apply
        module procedure apply_dk
      end interface

      interface x_ghost
        module procedure dk_ghost
      end interface

      contains

      function constructor_dk(dv,dw) result (dk)
!doc$ function dyad_kpoint(dv,dw) result(dk)
        type(dyad_kpoint_obj) :: dk
        type(multivector_obj), optional :: dv
        type(multivector_obj), optional :: dw
!       effects: Creates a new dk.
!       errors: dv and dw are not both present or both absent

!cod$
        if (present(dv)) call my(dv)
        if (present(dw)) call my(dw)

        dk%ref = 0
        allocate( dk%o )
        dk%o%ref = 0
        dk%o%g = x_ghost()
        nullify(dk%o%dv)
        nullify(dk%o%dw)

        if (error((present(dv) .neqv. present(dw)),"ERROR: dv and dw must become defined together")) goto 100

        if (present(dv)) then
           allocate(dk%o%dv)
           call my(dv,dk%o%dv)
           allocate(dk%o%dw)
           call my(dw,dk%o%dw)
        end if

100     if (present(dv)) call glean(thy(dv))
        if (present(dw)) call glean(thy(dw))

        if (error("Exit dyad_kpoint_mod::constructor_dk")) continue

      end function

      subroutine update_dk(dk,dv,dw)
!doc$ subroutine update(dk, dv, dw)
        type(dyad_kpoint_obj) :: dk
        type(multivector_obj), optional :: dv
        type(multivector_obj), optional :: dw
!       modifies: dk
!       effects: Updates dk.
!       errors: The dyad_kpoint is empty and either dv or dw, but not both, are present

!cod$
        call my(dk)
        if (present(dv)) call my(dv)
        if (present(dw)) call my(dw)

        if (present(dv)) then
          if (associated(dk%o%dv)) then
             if ( x_ghost(dk%o%dv) /= x_ghost(dv) ) then
                call own_i(dk)
                dk%o%g = x_ghost()
                dk%o%dv = dv
             end if
          else
             if (error((.not. present(dw)), "ERROR: dv and dw must become defined together")) goto 100
             call own_i(dk)
             dk%o%g = x_ghost()
             allocate(dk%o%dv)
             call my(dv,dk%o%dv)
          end if
        end if

        if (present(dw)) then
          if (associated(dk%o%dw)) then
             if ( x_ghost(dk%o%dw) /= x_ghost(dw) ) then
                call own_i(dk)
                dk%o%g = x_ghost()
                dk%o%dw = dw
             end if
          else
             if (error((.not. present(dv)), "ERROR: dv and dw must become defined together")) goto 100
             call own_i(dk)
             dk%o%g = x_ghost()
             allocate(dk%o%dw)
             call my(dw,dk%o%dw)
          end if
        end if

100     if (present(dv)) call glean(thy(dv))
        if (present(dw)) call glean(thy(dw))

        call glean(thy(dk))

        if (error("Exit dyad_kpoint_mod::update_dk")) continue

      end subroutine

      subroutine my_dk(dk)
!doc$ subroutine my(dk)
        type(dyad_kpoint_obj) :: dk

!cod$
        dk%ref = dk%ref + 1
        dk%o%ref = dk%o%ref + 1
      end subroutine

      subroutine my_new_dk(dki,dk)
!doc$ subroutine my(dki,dk)
        type(dyad_kpoint_obj) :: dki, dk

!cod$
        dk%ref = 1
        dk%o => dki%o
        dk%o%ref = dk%o%ref + 1
      end subroutine

      function thy_dk(dk) result(dko)
!doc$ function thy(dk) result(dko)
        type(dyad_kpoint_obj) :: dk
        type(dyad_kpoint_obj) :: dko

!cod$
        dk%ref = dk%ref - 1
        dk%o%ref = dk%o%ref - 1
        dko%ref = dk%ref
        dko%o => dk%o
      end function

      subroutine glean_dk(dk)
!doc$ subroutine glean(dk)
        type(dyad_kpoint_obj) :: dk

!cod$
        if (dk%o%ref < 1) then
          if (associated(dk%o%dv)) then
            call glean(thy(dk%o%dv))
            deallocate(dk%o%dv)
          end if 
          if (associated(dk%o%dw)) then
            call glean(thy(dk%o%dw))
            deallocate(dk%o%dw)
          end if
          deallocate( dk%o )
        end if
      end subroutine

      subroutine bequeath_dk(dk)
!doc$ subroutine bequeath(dk)
        type(dyad_kpoint_obj) :: dk

!cod$
        continue
      end subroutine

      subroutine assign_dk(dk,dk2)
!doc$ subroutine assignment(=)(dk,dk2)
        type(dyad_kpoint_obj), intent(inout) :: dk
        type(dyad_kpoint_obj), intent(in) :: dk2

!cod$
        type(dyad_kpoint_obj) :: dkt

        call my(dk2)
        dkt%o => dk%o
        dk%o%ref = dk%o%ref - dk%ref
        dk%o => dk2%o
        dk%o%ref = dk%o%ref + dk%ref
        call glean(dkt)
        call glean(thy(dk2))
      end subroutine

      subroutine own_i(dk)
        type(dyad_kpoint_obj) :: dk

        type(dyad_kpoint_obj) :: dkt
        if (dk%ref < dk%o%ref) then
          allocate( dkt%o )
          dkt%o%ref = 0
          if (associated(dk%o%dv)) then
             allocate(dkt%o%dv)
             call my(dk%o%dv,dkt%o%dv)
          else
             nullify(dkt%o%dv)
          end if
          if (associated(dk%o%dw)) then
             allocate(dkt%o%dw)
             call my(dk%o%dw,dkt%o%dw)
          else
             nullify(dkt%o%dw)
          end if
          dkt%o%g = dk%o%g
          dk%o%ref = dk%o%ref - dk%ref
          dk%o => dkt%o
          dk%o%ref = dk%o%ref + dk%ref
        end if
      end subroutine

      function is_empty_dk(dk) result(ie)
!doc$ function is_empty(dk) result(ie)
        type(dyad_kpoint_obj) :: dk
        logical :: ie
!       effects: Returns .true. iff dk is empty

!cod$
        call my(dk)
        ie = .not. associated(dk%o%dv)
        call glean(thy(dk))
      end function

      function dk_v(dk) result(dv)
!doc$ function x_v(dk) result(dv)
        type(dyad_kpoint_obj) :: dk
        type(multivector_obj) :: dv
!       effects: Returns the multivector dv contained in dk.
!       errors: dk is empty

!cod$
        call my(dk)
        if (error(is_empty(dk),"ERROR: dk is empty")) goto 100
        call my(dk%o%dv,dv)
        call bequeath(thy(dv))

100     call glean(thy(dk))
        if (error("Exit dyad_kpoint_mod::x_v")) continue

      end function

      function dk_w(dk) result(dw)
!doc$ function x_w(dk) result(dw)
        type(dyad_kpoint_obj) :: dk
        type(multivector_obj) :: dw
!       effects: Returns the multivector dw contained in dk.
!       errors: dk is empty

!cod$
        call my(dk)
        if (error(is_empty(dk),"ERROR: dk is empty")) goto 100
        call my(dk%o%dw,dw)
        call bequeath(thy(dw))

100     call glean(thy(dk))
        if (error("Exit dyad_kpoint_mod::x_w")) continue

      end function

      function dk_ghost(dk) result(g)
!doc$ function x_ghost(dk) result(g)
        type(dyad_kpoint_obj) :: dk
        type(ghost) :: g
!       effects: Returns the ghost of dk.
        
!cod$
        call my(dk)
        g = dk%o%g
        call glean(thy(dk))
      end function

      function apply_dk(dk,mv) result(rmv)
!doc$ function apply(dk,mv) result(rmv)
        type(dyad_kpoint_obj) :: dk
        type(multivector_obj) :: mv
        type(multivector_obj) :: rmv
!       requires: mv%o%usage = NORMAL
!       effects: Returns the result of applying the dyad_kpoint operator to mv
!       errors: dk is empty

!cod$
        complex(double), parameter :: c0  = ( 0.0_double,0.0_double)
        complex(double), parameter :: cp1 = (+1.0_double,0.0_double)
        character(line_len), parameter :: init = "zeros"
        complex(double), dimension(:,:), allocatable :: cmat

        call my(dk)
        call my(mv)
        allocate( cmat(x_n_bands(mv),x_n_bands(mv)))

        if (error(is_empty(dk),"ERROR: dk is empty")) goto 100

        call my(multivector(x_multibasis(mv),init),rmv) ; if (error()) goto 100
        call overlap(dk%o%dw,mv,cmat)
        call transform(c0,rmv,cp1,dk%o%dv,cmat)
        call overlap(dk%o%dv,mv,cmat)
        call transform(cp1,rmv,cp1,dk%o%dw,cmat)
        call bequeath(thy(rmv))

100     deallocate(cmat)
        call glean(thy(mv))
        call glean(thy(dk))
        if (error("Exit dyad_kpoint_mod::apply")) continue

      end function

      end module

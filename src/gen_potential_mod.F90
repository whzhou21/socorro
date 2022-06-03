!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module gen_potential_mod
!doc$ module gen_potential_mod

!     Defines a type(gen_potential_obj), which encapsulates the potential data involved in a calculation.  Data is "potential data"
!     if it captures the part of the electronic operators that change from one self-consistency iteration to the next.  It is
!     calculated from the density data in fields_mod and passed back down to operators.  Currently, the gen_potential_obj is a
!     container for a type(grid_obj), a type(atomic_potential_obj), and optionally a type(dyad_obj)

      use kind_mod
      use error_mod
      use grid_mod
      use ghost_mod
      use atomic_potential_mod
      use dyad_mod

!cod$
      implicit none
      private

      type :: gen_potential_rep
        integer :: ref
        type(ghost) :: g
        type(grid_obj) :: gpot
        type(atomic_potential_obj) :: apot
        type(dyad_obj), pointer :: dpot
      end type

      type, public :: gen_potential_obj
        private
        integer :: ref
        type(gen_potential_rep), pointer :: o
      end type

!doc$
      public :: gen_potential, update
      public :: my, thy, glean, bequeath, assignment(=)
      public :: x_grid_potential, x_atomic_potential, x_ghost
      public :: has_dyadic_potential, x_dyadic_potential

!cod$

      interface gen_potential
        module procedure constructor_gp
      end interface
      interface update
        module procedure update_gp
      end interface

      interface my
        module procedure my_gp, my_new_gp
      end interface
      interface thy
        module procedure thy_gp
      end interface
      interface glean
        module procedure glean_gp
      end interface
      interface bequeath
        module procedure bequeath_gp
      end interface
      interface assignment(=)
        module procedure assign_gp
      end interface

      interface x_grid_potential
        module procedure gp_grid_potential
      end interface
      interface x_atomic_potential
        module procedure gp_atomic_potential
      end interface
      interface x_ghost
        module procedure gp_ghost
      end interface

      interface has_dyadic_potential
        module procedure has_dyadic_potential_gp
      end interface
      interface x_dyadic_potential
        module procedure gp_dyadic_potential
      end interface

      contains

      function constructor_gp(gpot,apot,dpot) result (gp)
!doc$ function gen_potential(gpot,apot,dpot) result(gp)
        type(grid_obj) :: gpot
        type(gen_potential_obj) :: gp
        type(atomic_potential_obj) :: apot
        type(dyad_obj), optional :: dpot
!       effects: Creates a new gp.

!cod$
        call my(gpot)
        call my(apot)
        if (present(dpot)) call my(dpot)

        gp%ref = 0
        allocate( gp%o )
        gp%o%ref = 0
        gp%o%g = x_ghost()

        call my(gpot,gp%o%gpot)
        call my(apot,gp%o%apot)
        if (present(dpot)) then
           allocate(gp%o%dpot)
           call my(dpot,gp%o%dpot)
        else
           nullify(gp%o%dpot)
        end if

100     call glean(thy(gpot))
        call glean(thy(apot))
        if (present(dpot)) call glean(thy(dpot))

        if (error("Exit gen_potential_mod::constructor_gp")) continue

      end function

      subroutine update_gp(gp,gpot,apot,dpot)
!doc$ subroutine update(gp,gpot,apot)
        type(gen_potential_obj) :: gp
        type(grid_obj), optional :: gpot
        type(atomic_potential_obj), optional :: apot
        type(dyad_obj), optional :: dpot
!       modifies: gp
!       effects: Updates gp.

!cod$
        call my(gp)

        if (present(gpot)) call my(gpot)
        if (present(apot)) call my(apot)
        if (present(dpot)) call my(dpot)

        if (present(gpot)) then
          if ( x_ghost(gp%o%gpot) /= x_ghost(gpot) ) then
             call own_i(gp)
             gp%o%g = x_ghost()
             gp%o%gpot = gpot
          end if
        end if

        if (present(apot)) then
          if ( x_ghost(gp%o%apot) /= x_ghost(apot) ) then
             call own_i(gp)
             gp%o%g = x_ghost()
             gp%o%apot = apot
          end if
        end if

        if (present(dpot)) then
          if (associated(gp%o%dpot)) then
             if ( x_ghost(gp%o%dpot) /= x_ghost(dpot) ) then
                call own_i(gp)
                gp%o%g = x_ghost()
                gp%o%dpot = dpot
             end if
          else
             call own_i(gp)
             gp%o%g = x_ghost()
             allocate(gp%o%dpot)
             call my(dpot,gp%o%dpot)
          end if
        end if

100     if (present(gpot)) call glean(thy(gpot))
        if (present(apot)) call glean(thy(apot))
        if (present(dpot)) call glean(thy(dpot))

        call glean(thy(gp))

        if (error("Exit gen_potential_mod::update_gp")) continue

      end subroutine

      subroutine my_gp(gp)
!doc$ subroutine my(gp)
        type(gen_potential_obj) :: gp

!cod$
        gp%ref = gp%ref + 1
        gp%o%ref = gp%o%ref + 1
      end subroutine

      subroutine my_new_gp(gpi,gp)
!doc$ subroutine my(gpi,gp)
        type(gen_potential_obj) :: gpi, gp

!cod$
        gp%ref = 1
        gp%o => gpi%o
        gp%o%ref = gp%o%ref + 1
      end subroutine

      function thy_gp(gp) result(gpo)
!doc$ function thy(gp) result(gpo)
        type(gen_potential_obj) :: gp
        type(gen_potential_obj) :: gpo

!cod$
        gp%ref = gp%ref - 1
        gp%o%ref = gp%o%ref - 1
        gpo%ref = gp%ref
        gpo%o => gp%o
      end function

      subroutine glean_gp(gp)
!doc$ subroutine glean(gp)
        type(gen_potential_obj) :: gp

!cod$
        integer :: i
        if (gp%o%ref < 1) then
          call glean(thy(gp%o%gpot))
          call glean(thy(gp%o%apot))
          if (associated(gp%o%dpot)) then
            call glean(thy(gp%o%dpot))
            deallocate(gp%o%dpot)
          end if 
          deallocate( gp%o )
        end if
      end subroutine

      subroutine bequeath_gp(gp)
!doc$ subroutine bequeath(gp)
        type(gen_potential_obj) :: gp

!cod$
        continue
      end subroutine

      subroutine assign_gp(gp,gp2)
!doc$ subroutine assignment(=)(gp,gp2)
        type(gen_potential_obj), intent(inout) :: gp
        type(gen_potential_obj), intent(in) :: gp2

!cod$
        type(gen_potential_obj) :: gpt

        call my(gp2)
        gpt%o => gp%o
        gp%o%ref = gp%o%ref - gp%ref
        gp%o => gp2%o
        gp%o%ref = gp%o%ref + gp%ref
        call glean(gpt)
        call glean(thy(gp2))
      end subroutine

      subroutine own_i(gp)
        type(gen_potential_obj) :: gp

        integer :: i
        type(gen_potential_obj) :: gpt
        if (gp%ref < gp%o%ref) then
          allocate( gpt%o )
          gpt%o%ref = 0
          call my(gp%o%gpot,gpt%o%gpot)
          call my(gp%o%apot,gpt%o%apot)
          if (associated(gp%o%dpot)) then
             allocate(gpt%o%dpot)
             call my(gp%o%dpot,gpt%o%dpot)
          else
             nullify(gpt%o%dpot)
          end if
          gpt%o%g = gp%o%g
          gp%o%ref = gp%o%ref - gp%ref
          gp%o => gpt%o
          gp%o%ref = gp%o%ref + gp%ref
        end if
      end subroutine

      function gp_grid_potential(gp) result(gpot)
!doc$ function x_grid_potential(gp) result(gpot)
        type(gen_potential_obj) :: gp
        type(grid_obj) :: gpot
!       effects: Returns the gpot contained in gp.

!cod$
        call my(gp)
        call my(gp%o%gpot,gpot)
        call glean(thy(gp))
        call bequeath(thy(gpot))
      end function


      function gp_atomic_potential(gp) result(apot)
!doc$ function x_atomic_potential(gp) result(apot)
        type(gen_potential_obj) :: gp
        type(atomic_potential_obj) :: apot
!       effects: Returns the apot contained in gp.

!cod$
        call my(gp)
        call my(gp%o%apot,apot)
        call glean(thy(gp))
        call bequeath(thy(apot))
      end function

      function has_dyadic_potential_gp(gp) result(hd)
!doc$ function has_dyadic_potential(gp) result(hd)
        type(gen_potential_obj) :: gp
        logical :: hd
!       effects: Returns .true. iff dpot has been initialized

!cod$
        hd = associated(gp%o%dpot)
      end function

      function gp_dyadic_potential(gp) result(dpot)
!doc$ function x_dyadic_potential(gp) result(dpot)
        type(gen_potential_obj) :: gp
        type(dyad_obj) :: dpot
!       effects: Returns the dpot contained in gp.
!       errors: The dpot has not been initialized

!cod$
        call my(gp)
        if (error((.not. associated(gp%o%dpot)),"ERROR: dpot has not been initialized")) goto 100
        call my(gp%o%dpot,dpot)
        call bequeath(thy(dpot))

100     call glean(thy(gp))
        if (error("Exit gen_potential_mod::x_dyadic_potential")) continue

      end function

      function gp_ghost(gp) result(g)
!doc$ function x_ghost(gp) result(g)
        type(gen_potential_obj) :: gp
        type(ghost) :: g
!       effects: Returns the ghost of gp.
        
!cod$
        call my(gp)
        g = gp%o%g
        call glean(thy(gp))
      end function

      end module

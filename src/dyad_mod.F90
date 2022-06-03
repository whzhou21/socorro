!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module dyad_mod
!doc$ module dyad_mod

!     Defines a type(dyad_obj), which encapsulates a dyad_kpoint_obj for each k-point.  In normal usage, each processor
!     will only have data for the k-points assigned to its k-group, and the dyad_obj's for the other k-points will be "empty".

      use kind_mod
      use error_mod
      use dyad_kpoint_mod
      use ghost_mod

!cod$
      implicit none
      private

      type :: dyad_rep
        integer :: ref
        type(ghost) :: g
        type(dyad_kpoint_obj), dimension(:), pointer :: dyads
      end type

      type, public :: dyad_obj
        private
        integer :: ref
        type(dyad_rep), pointer :: o
      end type

!doc$
      public :: dyad, update
      public :: my, thy, glean, bequeath, assignment(=)
      public :: x_dyad_kpoint, x_n_kpoints, x_ghost

!cod$

      interface dyad
        module procedure constructor_do
      end interface
      interface update
        module procedure update_do
      end interface

      interface my
        module procedure my_do, my_new_do
      end interface
      interface thy
        module procedure thy_do
      end interface
      interface glean
        module procedure glean_do
      end interface
      interface bequeath
        module procedure bequeath_do
      end interface
      interface assignment(=)
        module procedure assign_do
      end interface

      interface x_dyad_kpoint
        module procedure do_dyad_kpoint
      end interface
      interface x_n_kpoints
        module procedure do_n_kpoints
      end interface
      interface x_ghost
        module procedure do_ghost
      end interface

      contains

      function constructor_do(nk) result (do)
!doc$ function dyad(nk) result(do)
        type(dyad_obj) :: do
        integer :: nk
!       effects: Creates a new do.

!cod$
        integer :: i

        do%ref = 0
        allocate( do%o )
        do%o%ref = 0
        do%o%g = x_ghost()

        allocate(do%o%dyads(nk))
        do i = 1, size(do%o%dyads)
           call my(dyad_kpoint(),do%o%dyads(i))
        end do

      end function

      subroutine update_do(do,ik,dk)
!doc$ subroutine update(do,ik,dk)
        type(dyad_obj) :: do
        integer :: ik
        type(dyad_kpoint_obj) :: dk
!       modifies: do
!       effects: Updates do.
!       errors: ik is out of range

!cod$
        call my(do)
        call my(dk)

        if (error((ik .lt. 1) .or. (ik .gt. size(do%o%dyads)),"ERROR: Index ik is out of range")) goto 100

        call own_i(do)
        do%o%g = x_ghost()
        do%o%dyads(ik) = dk

 100    call glean(thy(dk))
        call glean(thy(do))

        if (error("Exit dyad_mod::update_do")) continue

      end subroutine

      subroutine my_do(do)
!doc$ subroutine my(do)
        type(dyad_obj) :: do

!cod$
        do%ref = do%ref + 1
        do%o%ref = do%o%ref + 1
      end subroutine

      subroutine my_new_do(doi,do)
!doc$ subroutine my(doi,do)
        type(dyad_obj) :: doi, do

!cod$
        do%ref = 1
        do%o => doi%o
        do%o%ref = do%o%ref + 1
      end subroutine

      function thy_do(do) result(doo)
!doc$ function thy(do) result(doo)
        type(dyad_obj) :: do
        type(dyad_obj) :: doo

!cod$
        do%ref = do%ref - 1
        do%o%ref = do%o%ref - 1
        doo%ref = do%ref
        doo%o => do%o
      end function

      subroutine glean_do(do)
!doc$ subroutine glean(do)
        type(dyad_obj) :: do

!cod$
        integer :: i
        if (do%o%ref < 1) then
          do i = 1, size(do%o%dyads)
            call glean(thy(do%o%dyads(i)))
          end do
          deallocate(do%o%dyads)
          deallocate( do%o )
        end if
      end subroutine

      subroutine bequeath_do(do)
!doc$ subroutine bequeath(do)
        type(dyad_obj) :: do

!cod$
        continue
      end subroutine

      subroutine assign_do(do,do2)
!doc$ subroutine assignment(=)(do,do2)
        type(dyad_obj), intent(inout) :: do
        type(dyad_obj), intent(in) :: do2

!cod$
        type(dyad_obj) :: dot

        call my(do2)
        dot%o => do%o
        do%o%ref = do%o%ref - do%ref
        do%o => do2%o
        do%o%ref = do%o%ref + do%ref
        call glean(dot)
        call glean(thy(do2))
      end subroutine

      subroutine own_i(do)
        type(dyad_obj) :: do

        integer :: i
        type(dyad_obj) :: dot
        if (do%ref < do%o%ref) then
          allocate( dot%o )
          dot%o%ref = 0
          allocate(dot%o%dyads(size(do%o%dyads)))
          do i = 1, size(do%o%dyads)
             call my(do%o%dyads(i),dot%o%dyads(i))
          end do
          dot%o%g = do%o%g
          do%o%ref = do%o%ref - do%ref
          do%o => dot%o
          do%o%ref = do%o%ref + do%ref
        end if
      end subroutine

      function do_dyad_kpoint(do,ik) result(dk)
!doc$ function x_dyad_kpoint(do,ik) result(dk)
        type(dyad_obj) :: do
        integer :: ik
        type(dyad_kpoint_obj) :: dk
!       effects: Returns the ik'th dyad_kpoint contained in do.
!       errors: ik is out of range

!cod$
        call my(do)
        if (error(((ik .lt. 1) .or. (ik .gt. size(do%o%dyads))),"ERROR: Index ik is out of range")) goto 100
        call my(do%o%dyads(ik),dk)
        call bequeath(thy(dk))

100     call glean(thy(do))
        if (error("Exit dyad_mod::x_dyad_kpoint")) continue

      end function

      function do_n_kpoints(do) result(nk)
!doc$ function x_n_kpoints(do) result(nk)
        type(dyad_obj) :: do
        integer :: nk
!       effects: Returns the number of dyad_kpoint_obj's contained in do.

!cod$
        call my(do)
        nk = size(do%o%dyads)
100     call glean(thy(do))
        if (error("Exit dyad_mod::x_n_kpoints")) continue

      end function

      function do_ghost(do) result(g)
!doc$ function x_ghost(do) result(g)
        type(dyad_obj) :: do
        type(ghost) :: g
!       effects: Returns the ghost of do.
        
!cod$
        call my(do)
        g = do%o%g
        call glean(thy(do))
      end function

      end module

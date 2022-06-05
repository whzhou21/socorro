!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module check_symmetry_mod
!doc$ module check_symmetry_mod

!     Some info ...
!     ...

      use kind_mod
      use mpi_mod
      use error_mod
      use arg_mod
      use diary_mod
      use atoms_mod
      use crystal_mod
      use symmetry_mod
      use timing_mod

!cod$
      implicit none ; private

      logical :: found, sym_atoms
      character(line_len), parameter :: sym_cr_path = "sym_crystal"
      integer :: ia, na
      real(double), dimension(:,:), allocatable :: pos
      type(atoms_obj) :: sym_at
      type(crystal_obj) :: cr, sym_cr
      type(point_group_obj) :: lg
      type(space_group_obj) :: sg

!doc$
      public :: check_symmetry

!cod$
      interface check_symmetry
         module procedure check_symmetry_
      end interface

      contains

! public routines

      subroutine check_symmetry_()
!doc$ subroutine check_kpoints()
!        effects:
!        errors:
!        requires:

!cod$
         call start_timer("check_symmetry: total time")

         call my(crystal(),cr) ; if (error()) goto 100

         call my(point_group(x_lattice(cr)),lg) ; if (error()) goto 100
         call my(space_group(lg,x_atoms(cr),x_lattice(cr)),sg) ; if (error()) goto 100

         call diary(cr)
         call diary(lg)
         call diary(sg)

         call arg("symmetrize_atoms",sym_atoms,found)
         if (.not.found) sym_atoms = .false.
         if (sym_atoms) then
            call my(x_atoms(cr),sym_at)
            na = x_n_atoms(sym_at)
            allocate ( pos(3,na) )
            do ia = 1,na
               pos(:,ia) = x_position(sym_at,ia)
            end do
            call symmetrize_coordinates(sg,pos) ; if (error()) goto 100
            call move(sym_at,pos) ; if (error()) goto 100
            call my(crystal(sym_at,x_lattice(cr),x_name(cr)),sym_cr) ; if (error()) goto 100
            call save(sym_cr,path=sym_cr_path)
            call diary(sym_cr)
            call glean(thy(sym_at))
            call glean(thy(sym_cr))
         end if

         if (allocated( pos )) deallocate( pos )
         call glean(thy(cr))
         call glean(thy(lg))
         call glean(thy(sg))

100      if ( error("Exit check_symmetry") ) continue
         if ( .not.error() ) call stop_timer("check_symmetry: total time")

      end subroutine

      end module check_symmetry_mod

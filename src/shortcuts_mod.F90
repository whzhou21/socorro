!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module shortcuts_mod
!doc$ module shortcuts_mod

!     shortcuts_mod encapsulates simple procedures used by other modules. It provides
!     a way of skipping intervening levels, for example.

      use kind_mod
      use mpi_mod
      use error_mod
      use atoms_mod
      use lattice_mod
      use crystal_mod
      use external_mod
      use config_sc_mod
      use config_td_mod

!cod$
      implicit none
      private

!doc$
      public :: x_n_atoms
      public :: x_cart_position
      public :: x_cart_positions
      public :: x_lattice_vector
      public :: x_type
      public :: update_config

!cod$
      interface x_n_atoms
        module procedure n_atoms_cfg, n_atoms_cfg_td, n_atoms_cry
      end interface
      interface x_cart_position
        module procedure cart_position_cfg, cart_position_cfg_td, cart_position_cry
      end interface
      interface x_cart_positions
        module procedure cart_positions_cfg, cart_positions_cfg_td, cart_positions_cry
      end interface
      interface x_lattice_vector
        module procedure lattice_vector_cfg, lattice_vector_cfg_td
      end interface
      interface x_type
        module procedure atom_type_cfg, atom_type_cfg_td
      end interface
      interface update_config
        module procedure update_config_pos_cart, update_config_td_pos_cart
      end interface

      contains

      function n_atoms_cfg(cfg) result(na)
!doc$ function x_n_atoms(cfg) result(na)
        type(config_sc_obj) :: cfg
        integer :: na
!       effects: Returns the number of atoms in cfg.

!cod$
        call my(cfg)
        na = x_n_atoms(x_atoms(x_crystal(x_external(cfg))))
        call glean(thy(cfg))
      end function

      function n_atoms_cfg_td(cfg) result(na)
!doc$ function x_n_atoms(cfg) result(na)
        type(config_td_obj) :: cfg
        integer :: na
!       effects: Returns the number of atoms in cfg.

!cod$
        call my(cfg)
        na = x_n_atoms(x_atoms(x_crystal(x_external(cfg))))
        call glean(thy(cfg))
      end function

      function n_atoms_cry(cry) result(na)
!doc$ function x_n_atoms(cry) result(na)
        type(crystal_obj) :: cry
        integer :: na
!       effects: Returns the number of atoms in cry.

!cod$
        call my(cry)
        na = x_n_atoms(x_atoms(cry))
        call glean(thy(cry))
      end function

      function cart_position_cfg(cfg,ia) result(pos)
!doc$ function x_cart_position(cfg,ia) result(pos)
        type(config_sc_obj) :: cfg
        integer, intent(in) :: ia
        real(double), dimension(3) :: pos
!       requires: ia be between 1 and x_n_atoms(cfg)
!       effects: Returns the cartesian position of the ia'th atom in cfg.

!cod$
        call my(cfg)
        pos = x_position(x_atoms(x_crystal(x_external(cfg))),ia)
        pos = lat2r(x_lattice(x_crystal(x_external(cfg))),pos)
        call glean(thy(cfg))
      end function

      function cart_position_cfg_td(cfg,ia) result(pos)
!doc$ function x_cart_position(cfg,ia) result(pos)
        type(config_td_obj) :: cfg
        integer, intent(in) :: ia
        real(double), dimension(3) :: pos
!       requires: ia be between 1 and x_n_atoms(cfg)
!       effects: Returns the cartesian position of the ia'th atom in cfg.

!cod$
        call my(cfg)
        pos = x_position(x_atoms(x_crystal(x_external(cfg))),ia)
        pos = lat2r(x_lattice(x_crystal(x_external(cfg))),pos)
        call glean(thy(cfg))
      end function

      function cart_position_cry(cry,ia) result(pos)
!doc$ function x_cart_position(cry,ia) result(pos)
        type(crystal_obj) :: cry
        integer, intent(in) :: ia
        real(double), dimension(3) :: pos
!       requires: ia be between 1 and x_n_atoms(cry)
!       effects: Returns the cartesian position of the ia atom in cry.

!cod$
        call my(cry)
        pos = x_position(x_atoms(cry),ia)
        pos = lat2r(x_lattice(cry),pos)
        call glean(thy(cry))
      end function

      subroutine cart_positions_cfg(cfg,pos)
!doc$ subroutine x_cart_positions(cfg,pos)
        type(config_sc_obj) :: cfg
        real(double), dimension(:,:), intent(inout) :: pos
!       effects: Returns the cartesian positions of the atoms in cfg.

!cod$
        integer :: ia, na
        call my(cfg)
        na = x_n_atoms(x_atoms(x_crystal(x_external(cfg))))
        do ia = 1,na
          pos(:,ia) = x_cart_position(cfg,ia)
        end do
        call glean(thy(cfg))
      end subroutine

      subroutine cart_positions_cfg_td(cfg,pos)
!doc$ subroutine x_cart_positions(cfg,pos)
        type(config_td_obj) :: cfg
        real(double), dimension(:,:), intent(inout) :: pos
!       effects: Returns the cartesian positions of the atoms in cfg.

!cod$
        integer :: ia, na
        call my(cfg)
        na = x_n_atoms(x_atoms(x_crystal(x_external(cfg))))
        do ia = 1,na
          pos(:,ia) = x_cart_position(cfg,ia)
        end do
        call glean(thy(cfg))
      end subroutine

      subroutine cart_positions_cry(cry,pos)
!doc$ subroutine x_cart_positions(cry,pos)
        type(crystal_obj) :: cry
        real(double), dimension(:,:), intent(inout) :: pos
!       effects: Returns the cartesian positions of the atoms in cry.

!cod$
        integer :: ia, na
        call my(cry)
        na = x_n_atoms(x_atoms(cry))
        do ia = 1,na
          pos(:,ia) = x_cart_position(cry,ia)
        end do
        call glean(thy(cry))
      end subroutine

      function lattice_vector_cfg(cfg,i) result(vector)
!doc$ function x_lattice_vector(cfg) result(vector)
        type(config_sc_obj) :: cfg
        integer, intent(in) :: i
        real(double), dimension(3) :: vector
!       effects: Returns the ith lattice vector

!cod$
        call my(cfg)
        vector = x_lattice_vector(x_lattice(x_crystal(x_external(cfg))),i)
        call glean(thy(cfg))
      end function

      function lattice_vector_cfg_td(cfg,i) result(vector)
!doc$ function x_lattice_vector(cfg) result(vector)
        type(config_td_obj) :: cfg
        integer, intent(in) :: i
        real(double), dimension(3) :: vector
!       effects: Returns the ith lattice vector

!cod$
        call my(cfg)
        vector = x_lattice_vector(x_lattice(x_crystal(x_external(cfg))),i)
        call glean(thy(cfg))
      end function

      function atom_type_cfg(cfg,ia) result(atg)
!doc$ function x_type(cfg,ia) result(atg)
        type(config_sc_obj) :: cfg
        integer, intent(in) :: ia
        character(tag_sz) :: atg
!       requires: ia be between 1 and x_n_atoms(cfg)
!       effects: Returns the identifying atom tag of the ia'th atom in cfg.

!cod$
        call my(cfg)
        atg = x_type(x_atoms(x_crystal(x_external(cfg))),ia)
        call glean(thy(cfg))
      end function

      function atom_type_cfg_td(cfg,ia) result(atg)
!doc$ function x_type(cfg,ia) result(atg)
        type(config_td_obj) :: cfg
        integer, intent(in) :: ia
        character(tag_sz) :: atg
!       requires: ia be between 1 and x_n_atoms(cfg)
!       effects: Returns the identifying atom tag of the ia'th atom in cfg.

!cod$
        call my(cfg)
        atg = x_type(x_atoms(x_crystal(x_external(cfg))),ia)
        call glean(thy(cfg))
      end function

      subroutine update_config_pos_cart(cfg,pos_cart)
!doc$ subroutine update_config(cfg,pos_cart)
        type(config_sc_obj) :: cfg
        real(double), dimension(:,:), intent(in) :: pos_cart
!       modifies: cfg
!       requires: Positions be in cartesian coordinates.
!       effects: Updates atom positions in cfg.
!       errors: Passes errors.

!cod$
        integer :: ia, na
        real(double), dimension(size(pos_cart,1),size(pos_cart,2)) :: pos_lat
        type(lattice_obj) :: lat
        type(atoms_obj) :: at_tmp
        type(crystal_obj) :: cr_tmp
        type(external_obj) :: ex_tmp

        call my(cfg)

        call my(x_lattice(x_crystal(x_external(cfg))),lat)
        na = size(pos_cart,2)
        do ia = 1,na
          pos_lat(:,ia) = r2lat(lat,pos_cart(:,ia))
        end do
        call glean(thy(lat))

        call my(x_external(cfg),ex_tmp)
        call my(x_crystal(ex_tmp),cr_tmp)
        call my(x_atoms(cr_tmp),at_tmp)

        call move(at_tmp,pos_lat) ; if (error()) goto 100
        call update(cr_tmp,at_tmp) ; if (error()) goto 100
        call update(ex_tmp,cr_tmp) ; if (error()) goto 100
        call update(cfg,ex_tmp) ; if (error()) goto 100

100     call glean(thy(at_tmp))
        call glean(thy(cr_tmp))
        call glean(thy(ex_tmp))

        call glean(thy(cfg))

        if (error("Exit shortcuts_mod::update_config_pos_cart")) continue

      end subroutine

      subroutine update_config_td_pos_cart(cfg,pos_cart)
!doc$ subroutine update_config(cfg,pos_cart)
        type(config_td_obj) :: cfg
        real(double), dimension(:,:), intent(in) :: pos_cart
!       modifies: cfg
!       requires: Positions be in cartesian coordinates.
!       effects: Updates atom positions in cfg.
!       errors: Passes errors.

!cod$
        integer :: ia, na
        real(double), dimension(size(pos_cart,1),size(pos_cart,2)) :: pos_lat
        type(lattice_obj) :: lat
        type(atoms_obj) :: at_tmp
        type(crystal_obj) :: cr_tmp
        type(external_obj) :: ex_tmp

        call my(cfg)

        call my(x_lattice(x_crystal(x_external(cfg))),lat)
        na = size(pos_cart,2)
        do ia = 1,na
          pos_lat(:,ia) = r2lat(lat,pos_cart(:,ia))
        end do
        call glean(thy(lat))

        call my(x_external(cfg),ex_tmp)
        call my(x_crystal(ex_tmp),cr_tmp)
        call my(x_atoms(cr_tmp),at_tmp)

        call move(at_tmp,pos_lat) ; if (error()) goto 100
        call update(cr_tmp,at_tmp) ; if (error()) goto 100
        call update(ex_tmp,cr_tmp) ; if (error()) goto 100
        call update(cfg,ex_tmp) ; if (error()) goto 100

100     call glean(thy(at_tmp))
        call glean(thy(cr_tmp))
        call glean(thy(ex_tmp))

        call glean(thy(cfg))

        if (error("Exit shortcuts_mod::update_config_pos_cart")) continue

      end subroutine

      end module

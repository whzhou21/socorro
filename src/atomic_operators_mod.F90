!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module atomic_operators_mod
!doc$ module atomic_operators_mod

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use diary_mod
      use arg_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use crystal_mod
      use symmetry_mod
      use xc_type_mod
      use atomic_operators_ncp_mod
      use atomic_operators_paw_mod

!     One datatype is defined here: type(atomic_operators_obj).
     
!cod$
      implicit none
      private

      integer, parameter, public :: NCP = 1
      integer, parameter, public :: PAW = 2

      type, public :: atomic_operators_obj
        private
        integer :: type
        type(atomic_operators_ncp_obj) :: ao_ncp
        type(atomic_operators_paw_obj) :: ao_paw
      end type

!doc$
      public :: atomic_operators
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_type
      public :: x_crystal
      public :: x_atomic_operators
      public :: x_atomic_operators_ncp
      public :: x_atomic_operators_paw
      public :: x_n_types
      public :: x_n_atoms
      public :: x_type_name
      public :: x_type_valence
      public :: x_valence_electrons
      public :: x_type_matching_radius
      public :: x_atom_matching_radius
      public :: x_n_projectors
      public :: x_n_type_projectors
      public :: x_n_atom_projectors
      public :: x_projector_l
      public :: x_projector_m
      public :: x_projector_type
      public :: x_projector_atom
      public :: x_projector_index_in_type
      public :: projector_f_values
      public :: projector_stress_f_values
      public :: projector_r_value
      public :: projector_r_gradients
      public :: x_projector_radius
      public :: x_type_projector_radius
      public :: atomic_grid_potential
      public :: atomic_xc_density
      public :: diary_representation
      public :: diary_angular_mesh
      public :: diary_occupation
      public :: diary_rs_projectors
      public :: write_restart

!cod$
      interface atomic_operators
        module procedure constructor_ao, constructor_ao_ncp, constructor_ao_paw
      end interface
      interface update
        module procedure update_ao
      end interface
      interface my
        module procedure my_ao, my_new_ao
      end interface
      interface thy
        module procedure thy_ao
      end interface
      interface glean
        module procedure glean_ao
      end interface
      interface bequeath
        module procedure bequeath_ao
      end interface
      interface assignment(=)
        module procedure assign_ao
      end interface
      interface x_ref
        module procedure ao_ref
      end interface
      interface x_ghost
        module procedure ao_ghost
      end interface
      interface x_type
        module procedure ao_type
      end interface
      interface x_crystal
        module procedure ao_crystal
      end interface
      interface x_atomic_operators
        module procedure ao_atomic_operators
      end interface
      interface x_atomic_operators_ncp
        module procedure ao_atomic_operators_ncp
      end interface
      interface x_atomic_operators_paw
        module procedure ao_atomic_operators_paw
      end interface
      interface x_n_types
        module procedure ao_n_types
      end interface
      interface x_n_atoms
        module procedure ao_n_atoms
      end interface
      interface x_type_name
        module procedure ao_type_name
      end interface
      interface x_type_valence
        module procedure ao_type_valence
      end interface
      interface x_valence_electrons
        module procedure ao_valence_electrons
      end interface
      interface x_type_matching_radius
        module procedure ao_type_matching_radius
      end interface
      interface x_atom_matching_radius
        module procedure ao_atom_matching_radius
      end interface
      interface x_n_projectors
        module procedure ao_n_projectors
      end interface
      interface x_n_type_projectors
        module procedure ao_n_type_projectors
      end interface
      interface x_n_atom_projectors
        module procedure ao_n_atom_projectors
      end interface
      interface x_projector_l
        module procedure ao_projector_l, ao_type_projector_l
      end interface
      interface x_projector_m
        module procedure ao_projector_m, ao_type_projector_m
      end interface
      interface x_projector_type
        module procedure ao_projector_type
      end interface
      interface x_projector_index_in_type
        module procedure ao_projector_index_in_type
      end interface
      interface x_projector_atom
        module procedure ao_projector_atom
      end interface
      interface projector_f_values
        module procedure ao_projector_f_values, ao_type_projector_f_values
      end interface
      interface projector_stress_f_values
        module procedure ao_projector_stress_f_values
      end interface
      interface projector_r_value
        module procedure ao_projector_r_value
      end interface
      interface projector_r_gradients
        module procedure ao_projector_r_gradients
      end interface
      interface x_projector_radius
        module procedure ao_projector_radius
      end interface
      interface x_type_projector_radius
        module procedure ao_type_projector_radius
      end interface
      interface atomic_grid_potential
        module procedure atomic_grid_potential_ao
      end interface
      interface atomic_xc_density
        module procedure atomic_xc_density_ao
      end interface
      interface diary_representation
        module procedure diary_representation_ao
      end interface
      interface diary_angular_mesh
        module procedure diary_angular_mesh_ao
      end interface
      interface diary_occupation
        module procedure diary_occupation_ao
      end interface
      interface diary_rs_projectors
        module procedure diary_rs_projectors_ao
      end interface
      interface write_restart
        module procedure write_restart_ao
      end interface

      contains

! public routines

      function constructor_ao(cr,lay,sg,xct,restf) result(ao)
!doc$ function atomic_operators(cr,lay,sg,xct,restf) result(ao)
        type(crystal_obj) :: cr
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_type_obj) :: xct
        type(tagio_obj), optional :: restf
        type(atomic_operators_obj) :: ao
!       effects: Creates a new ao.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: tag

        if (present(restf)) call my(restf)

        if (present(restf)) then

          ! open the ATOMIC_TYPE block
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_TYPE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_TYPE block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)

          ! find the atomic representation tag
          if (i_access(restf)) tios = findfirsttag(restf,"PAW")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (tios == TAG_NOT_FOUND) then
            if (i_access(restf)) tios = findfirsttag(restf,"NCP")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: atomic representation  was not found")) goto 100
            ao%type = NCP
          else
            ao%type = PAW
          end if

          ! close the ATOMIC_TYPE block
100       if (i_access(restf)) call closeblock(restf)
          if (error()) goto 200

        else

          ! find the atomic representation tag
          call arglc("atomic_representation",tag,found)
          if (.not.found) tag = "ncp"
          select case (trim(tag))
          case ("ncp")
            ao%type = NCP
          case ("paw")
            ao%type = PAW
          case default
            if (error(.true.,"ERROR: atomic_representation tag was not recognized")) goto 200
          end select

        end if

        ! construct the atomic_operators
        select case (ao%type)
        case (NCP)
          if (present(restf)) then
            call my(atomic_operators_ncp(cr,lay,restf),ao%ao_ncp) ; if (error()) goto 200
          else
            call my(atomic_operators_ncp(cr,lay),ao%ao_ncp)       ; if (error()) goto 200
          end if
          call glean(sg)
          call bequeath(thy(ao%ao_ncp))
        case (PAW)
          if (present(restf)) then
            call my(atomic_operators_paw(cr,lay,sg,xct,restf),ao%ao_paw) ; if (error()) goto 200
          else
            call my(atomic_operators_paw(cr,lay,sg,xct),ao%ao_paw)       ; if (error()) goto 200
          end if
          call bequeath(thy(ao%ao_paw))
        end select

200     if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_operators_mod::constructor_ao")) continue

      end function

      function constructor_ao_ncp(ao_ncp) result(ao)
!doc$ function atomic_operators(ao_ncp) result(ao)
        type(atomic_operators_ncp_obj) :: ao_ncp
        type(atomic_operators_obj) :: ao
!       effects: Creates a new ao.

!cod$
        ao%type = NCP
        call my(ao_ncp,ao%ao_ncp)
        call bequeath(thy(ao%ao_ncp))
      end function

      function constructor_ao_paw(ao_paw) result(ao)
!doc$ function atomic_operators(ao_paw) result(ao)
        type(atomic_operators_paw_obj) :: ao_paw
        type(atomic_operators_obj) :: ao
!       effects: Creates a new ao.

!cod$
        ao%type = PAW
        call my(ao_paw,ao%ao_paw)
        call bequeath(thy(ao%ao_paw))
      end function

      subroutine update_ao(ao,cr,lay,sg)
!doc$ subroutine update(ao,cr,lay,sg)
        type(atomic_operators_obj) :: ao
        type(crystal_obj) :: cr
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
!       effects: Updates ao.

!cod$
        select case (ao%type)
        case (NCP)
          call update(ao%ao_ncp,cr,lay)
        case (PAW)
          call update(ao%ao_paw,cr,lay,sg)
        end select
        if (error("Exit atomic_operators_mod::update_ao")) continue
      end subroutine

      subroutine my_ao(ao)
!doc$ subroutine my(ao)
        type(atomic_operators_obj) :: ao

!cod$
        select case (ao%type)
        case (NCP)
          call my(ao%ao_ncp)
        case (PAW)
          call my(ao%ao_paw)
        end select
      end subroutine

      subroutine my_new_ao(aoi,ao)
!doc$ subroutine my(aoi,ao)
        type(atomic_operators_obj) :: aoi, ao

!cod$
        ao%type = aoi%type
        select case (ao%type)
        case (NCP)
          call my(aoi%ao_ncp,ao%ao_ncp)
        case (PAW)
          call my(aoi%ao_paw,ao%ao_paw)
        end select
      end subroutine

      function thy_ao(ao) result(aoo)
!doc$ function thy(ao) result(aoo)
        type(atomic_operators_obj) :: ao, aoo

!cod$
        aoo%type = ao%type
        select case (ao%type)
        case (NCP)
          call my(thy(ao%ao_ncp),aoo%ao_ncp)
          call bequeath(thy(aoo%ao_ncp))
        case (PAW)
          call my(thy(ao%ao_paw),aoo%ao_paw)
          call bequeath(thy(aoo%ao_paw))
        end select
      end function

      subroutine glean_ao(ao)
!doc$ subroutine glean(ao)
        type(atomic_operators_obj) :: ao

!cod$
        select case (ao%type)
        case (NCP)
          call glean(ao%ao_ncp)
        case (PAW)
          call glean(ao%ao_paw)
        end select
      end subroutine

      subroutine bequeath_ao(ao)
!doc$ subroutine bequeath(ao)
        type(atomic_operators_obj) :: ao

!cod$
        select case (ao%type)
        case (NCP)
          call bequeath(ao%ao_ncp)
        case (PAW)
          call bequeath(ao%ao_paw)
        end select
      end subroutine
 
      subroutine assign_ao(ao,ao2)
!doc$ subroutine assignment(=)(ao,ao2)
        type(atomic_operators_obj), intent(inout) :: ao
        type(atomic_operators_obj), intent(in) :: ao2
!       requires: ao and ao2 have the same type.

!cod$
        select case (ao%type)
        case (NCP)
          ao%ao_ncp = ao2%ao_ncp
        case (PAW)
          ao%ao_paw = ao2%ao_paw
        end select
      end subroutine

      function ao_ref(ao) result(r)
!doc$ function x_ref(ao) result(r)
        type(atomic_operators_obj) :: ao
        integer, dimension(2) :: r
!       effects: Returns ao%ref and ao%o%ref.

!cod$
        select case (ao%type)
        case (NCP)
          r = x_ref(ao%ao_ncp)
        case (PAW)
          r = x_ref(ao%ao_paw)
        end select
      end function
 
      function ao_ghost(ao) result(g)
!doc$ function x_ghost(ao) result(g)
        type(atomic_operators_obj) :: ao
        type(ghost) :: g
!       effects: Returns the ghost of ao.

!cod$
        select case (ao%type)
        case (NCP)
          g = x_ghost(ao%ao_ncp)
        case (PAW)
          g = x_ghost(ao%ao_paw)
        end select
      end function

      function ao_type(ao) result(t)
!doc$ function x_type(ao) result(t)
        type(atomic_operators_obj) :: ao
        integer :: t
!       effects: Returns ao%type.

!cod$
        t = ao%type
        call glean(ao)
      end function
    
      function ao_crystal(ao) result(cr)
!doc$ function x_crystal(ao) result(cr)
        type(atomic_operators_obj) :: ao
        type(crystal_obj) :: cr
!       effects: Returns the crystal of ao.

!cod$
        select case (ao%type)
        case (NCP)
          call my(x_crystal(ao%ao_ncp),cr)
        case (PAW)
          call my(x_crystal(ao%ao_paw),cr)
        end select
        call bequeath(thy(cr))
      end function

      function ao_atomic_operators(ao) result(aoc)
!doc$ function x_atomic_operators(ao) result(aoc)
        type(atomic_operators_obj) :: ao, aoc
!       effects: Returns a copy of ao.

!cod$
        call my(ao,aoc)
        call bequeath(thy(aoc))
      end function

       function ao_atomic_operators_ncp(ao) result(ao_ncp)
!doc$  function x_atomic_operators_ncp(ao) result(ao_ncp)
        type(atomic_operators_obj) :: ao
        type(atomic_operators_ncp_obj) :: ao_ncp
!       requires: x_type(ao) = NCP
!       effects: Returns ao%ao_ncp.

!cod$
        call my(ao%ao_ncp,ao_ncp)
        call bequeath(thy(ao_ncp))
      end function

       function ao_atomic_operators_paw(ao) result(ao_paw)
!doc$  function x_atomic_operators_paw(ao) result(ao_paw)
        type(atomic_operators_obj) :: ao
        type(atomic_operators_paw_obj) :: ao_paw
!       requires: x_type(ao) = PAW
!       effects: Returns ao%ao_paw.

!cod$
        call my(ao%ao_paw,ao_paw)
        call bequeath(thy(ao_paw))
      end function

      function ao_n_types(ao) result(nt)
!doc$ function x_n_types(ao) result(nt)
        type(atomic_operators_obj) :: ao
        integer :: nt
!       effects: Returns the number of atom types.

!cod$
        select case (ao%type)
        case (NCP)
          nt = x_n_types(ao%ao_ncp)
        case (PAW)
          nt = x_n_types(ao%ao_paw)
        end select
      end function

      function ao_n_atoms(ao) result(na)
!doc$ function x_n_atoms(ao) result(na)
        type(atomic_operators_obj) :: ao
        integer :: na
!       effects: Returns the number of atoms.

!cod$
        select case (ao%type)
        case (NCP)
          na = x_n_atoms(ao%ao_ncp)
        case (PAW)
          na = x_n_atoms(ao%ao_paw)
        end select
      end function

      function ao_type_name(ao,it) result(tag)
!doc$ function x_type_name(ao,it) result(tag)
        type(atomic_operators_obj) :: ao
        integer :: it
        character(tag_sz) :: tag
!       effects: Returns the name of type it.

!cod$
        select case (ao%type)
        case (NCP)
          tag = x_type_name(ao%ao_ncp,it)
        case (PAW)
          tag = x_type_name(ao%ao_paw,it)
        end select
      end function

      function ao_type_valence(ao,it) result(v)
!doc$ function x_type_valence(ao,it) result(v)
        type(atomic_operators_obj) :: ao
        integer :: it
        real(double) :: v
!       effects: Returns the valence of type it.

!cod$
        select case (ao%type)
        case (NCP)
          v = x_type_valence(ao%ao_ncp,it)
        case (PAW)
          v = x_type_valence(ao%ao_paw,it)
        end select
      end function

      function ao_valence_electrons(ao) result(ve)
!doc$ function x_valence_electrons(ao) result(ve)
        type(atomic_operators_obj) :: ao
        real(double) :: ve
!       effects: Returns the sum of the atomic valences.

!cod$
        select case (ao%type)
        case (NCP)
          ve = x_valence_electrons(ao%ao_ncp)
        case (PAW)
          ve = x_valence_electrons(ao%ao_paw)
        end select
      end function
    
      function ao_type_matching_radius(ao,it) result(mr)
!doc$ function x_type_matching_radius(ao,it) result(mr)
        type(atomic_operators_obj) :: ao
        integer :: it
        real(double) :: mr
!       effects: Returns the matching radius of atom type it.

!cod$
        select case (ao%type)
        case (NCP)
          mr = 0.0_double
        case (PAW)
          mr = x_type_matching_radius(ao%ao_paw,it)
        end select
      end function
    
      function ao_atom_matching_radius(ao,ia) result(mr)
!doc$ function x_atom_matching_radius(ao,ia) result(mr)
        type(atomic_operators_obj) :: ao
        integer :: ia
        real(double) :: mr
!       effects: Returns the matching radius of atom ia.

!cod$
        select case (ao%type)
        case (NCP)
          mr = 0.0_double
        case (PAW)
          mr = x_atom_matching_radius(ao%ao_paw,ia)
        end select
      end function
    
      function ao_n_projectors(ao) result(np)
!doc$ function x_n_projectors(ao) result(np)
        type(atomic_operators_obj) :: ao
        integer :: np
!       effects: Returns the total number of projectors.

!cod$
        select case (ao%type)
        case (NCP)
          np = x_n_projectors(ao%ao_ncp)
        case (PAW)
          np = x_n_projectors(ao%ao_paw)
        end select
      end function

      function ao_n_type_projectors(ao,it) result(np)
!doc$ function x_n_type_projectors(ao,it) result(np)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: it
        integer :: np
!       effects: Returns the number of projectors for atom type it.
!       errors: Passes errrors.

!cod$
        select case (ao%type)
        case (NCP)
          np = x_n_type_projectors(ao%ao_ncp,it)
        case (PAW)
          np = x_n_type_projectors(ao%ao_paw,it)
        end select
      end function

      function ao_n_atom_projectors(ao,ia) result(np)
!doc$ function x_n_atom_projectors(ao,ia) result(np)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ia
        integer :: np
!       effects: Returns the number of projectors for atom ia.
!       errors: Passes errrors.

!cod$
        select case (ao%type)
        case (NCP)
          np = x_n_atom_projectors(ao%ao_ncp,ia)
        case (PAW)
          np = x_n_atom_projectors(ao%ao_paw,ia)
        end select
      end function

      function ao_projector_l(ao,ip) result(l)
!doc$ function x_projector_l(ao,ip) result(l)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        integer :: l
!       effects: Returns the l value of projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          l = x_projector_l(ao%ao_ncp,ip)
        case (PAW)
          l = x_projector_l(ao%ao_paw,ip)
        end select
      end function

      function ao_type_projector_l(ao,it,ip) result(l)
!doc$ function x_projector_l(ao,it,ip) result(l)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: it, ip
        integer :: l
!       effects: Returns the l value of projector ip of type it.

!cod$
        select case (ao%type)
        case (NCP)
          l = x_projector_l(ao%ao_ncp,it,ip)
        case (PAW)
          l = x_projector_l(ao%ao_paw,it,ip)
        end select
      end function

      function ao_projector_m(ao,ip) result(m)
!doc$ function x_projector_m(ao,ip) result(m)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        integer :: m
!       effects: Returns the m value of projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          m = x_projector_m(ao%ao_ncp,ip)
        case (PAW)
          m = x_projector_m(ao%ao_paw,ip)
        end select
      end function

      function ao_type_projector_m(ao,it,ip) result(m)
!doc$ function x_projector_m(ao,it,ip) result(m)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: it, ip
        integer :: m
!       effects: Returns the m value of projector ip of type it.

!cod$
        select case (ao%type)
        case (NCP)
          m = x_projector_m(ao%ao_ncp,it, ip)
        case (PAW)
          m = x_projector_m(ao%ao_paw,it,ip)
        end select
      end function

      function ao_projector_type(ao,ip) result(it)
!doc$ function x_projector_type(ao,ip) result(it)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        integer :: it
!       effects: Returns the type of projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          it = x_projector_type(ao%ao_ncp,ip)
        case (PAW)
          it = x_projector_type(ao%ao_paw,ip)
        end select
      end function

      function ao_projector_atom(ao,ip) result(ia)
!doc$ function x_projector_atom(ao,ip) result(ia)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        integer :: ia
!       effects: Returns the indes of the atom associated with projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          ia = x_projector_atom(ao%ao_ncp,ip)
        case (PAW)
          ia = x_projector_atom(ao%ao_paw,ip)
        end select
      end function

      function ao_projector_index_in_type(ao,ip) result(i)
!doc$ function x_projector_index_in_type(ao,ip) result(i)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        integer :: i
!       effects: Returns the in-type index of projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          i = x_projector_index_in_type(ao%ao_ncp,ip)
        case (PAW)
          i = x_projector_index_in_type(ao%ao_paw,ip)
        end select
      end function

      function ao_projector_f_values(ao,q,ip) result(pfv)
!doc$ function projector_f_values(ao,q,ip) result(pfv)
        type(atomic_operators_obj) :: ao
        real(double), dimension(:), intent(in)  :: q
        integer, intent(in) :: ip
        real(double), dimension(size(q)) :: pfv
!       effects: Returns the q Fourier coefficients of projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          pfv = projector_f_values(ao%ao_ncp,q,ip)
        case (PAW)
          pfv = projector_f_values(ao%ao_paw,q,ip)
        end select
      end function

      function ao_type_projector_f_values(ao,q,it,ip) result(pfv)
!doc$ function projector_f_values(ao,q,it,ip) result(pfv)
        type(atomic_operators_obj) :: ao
        real(double), dimension(:), intent(in)  :: q
        integer, intent(in) :: it, ip
        real(double), dimension(size(q)) :: pfv
!       effects: Returns the q Fourier coefficients of projector ip of type it.

!cod$
        select case (ao%type)
        case (NCP)
          pfv = projector_f_values(ao%ao_ncp,q,it,ip)
        case (PAW)
          pfv = projector_f_values(ao%ao_paw,q,it,ip)
        end select
      end function

      function ao_projector_stress_f_values(ao,q,it,ip) result(psfv)
!doc$ function projector_stress_f_values(ao,q,it,ip) result(psfv)
        type(atomic_operators_obj) :: ao
        real(double), dimension(:), intent(in)  :: q
        integer, intent(in) :: it, ip
        real(double), dimension(2,size(q)) :: psfv
!       effects: Returns the sress tensor contributions from the q Fourier coefficients of projector ip of type it.

!cod$
        select case (ao%type)
        case (NCP)
          psfv = projector_stress_f_values(ao%ao_ncp,q,it,ip)
        case (PAW)
         if (error(.true.,"ERROR: stress tensor calculation is not yet implemented for PAW")) continue
        end select
        if (error("Exit atomic_operators_mod::ao_projector_stress_f_values")) continue
      end function

      function ao_projector_r_value(ao,ip,r,gi,go) result(prv)
!doc$ function projector_r_value(ao,ip,r,gi,go) result(prv)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        real(double), intent(in) :: r
        real(double), intent(in) :: gi, go
        real(double) :: prv
!       effects: Returns the value of projector ip at r.

!cod$
        select case (ao%type)
        case (NCP)
          prv = projector_r_value(ao%ao_ncp,ip,r,gi,go)
        case (PAW)
          prv = projector_r_value(ao%ao_paw,ip,r,gi,go)
        end select
      end function

      function ao_projector_r_gradients(ao,ip,r) result(prg)
!doc$ function projector_r_gradients(ao,ip,r) result(prg)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        real(double), intent(in) :: r
        real(double), dimension(2) :: prg
!       effects: Returns the gradients of projector ip at r.

!cod$
        select case (ao%type)
        case (NCP)
          prg = projector_r_gradients(ao%ao_ncp,ip,r)
        case (PAW)
          prg = projector_r_gradients(ao%ao_paw,ip,r)
        end select
      end function

      function ao_projector_radius(ao,ip) result(r)
!doc$ function x_projector_radius(ao,ip) result(r)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: ip
        real(double) :: r
!       effects: Returns the optimization radius of projector ip.

!cod$
        select case (ao%type)
        case (NCP)
          r = x_projector_radius(ao%ao_ncp,ip)
        case (PAW)
          r = x_projector_radius(ao%ao_paw,ip)
        end select
      end function

      function ao_type_projector_radius(ao,it) result(r)
!doc$ function x_type_projector_radius(ao,it) result(r)
        type(atomic_operators_obj) :: ao
        integer, intent(in) :: it
        real(double) :: r
!       effects: Returns the optimization radius of type it projectors.

!cod$
        select case (ao%type)
        case (NCP)
          r = x_type_projector_radius(ao%ao_ncp,it)
        case (PAW)
          r = x_type_projector_radius(ao%ao_paw,it)
        end select
      end function

      function atomic_grid_potential_ao(ao,lay) result(gp)
!doc$ function atomic_grid_potential(ao,lay) result(gp)
        type(atomic_operators_obj) :: ao
        type(layout_obj) :: lay
        type(grid_obj) :: gp
!       effects: Returns the atomic contribution to the grid potential.

!cod$
        select case (ao%type)
        case (NCP)
          call my(atomic_grid_potential(ao%ao_ncp,lay),gp)
        case (PAW)
          call my(atomic_grid_potential(ao%ao_paw,lay),gp)
        end select
        call bequeath(thy(gp))
      end function

      function atomic_xc_density_ao(ao,lay) result(xcd)
!doc$ function atomic_xc_density(ao,lay) result(xcd)
        type(atomic_operators_obj) :: ao
        type(layout_obj) :: lay
        type(grid_obj) :: xcd
!       effects: Returns the atomic contribution to the exchange-correlation density.

!cod$
        select case (ao%type)
        case (NCP)
          call my(atomic_xc_density(ao%ao_ncp,lay),xcd)
        case (PAW)
          call my(atomic_xc_density(ao%ao_paw,lay),xcd)
        end select
        call bequeath(thy(xcd))
      end function

      subroutine diary_representation_ao(ao)
!doc$ subroutine diary_representation(ao)
        type(atomic_operators_obj) :: ao
!       effects: Writes information about the atomic representation to the diary.

!cod$
        select case (ao%type)
        case (NCP)
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Using the Norm-Conserving Pseudopotential method")')
        case (PAW)
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Using the Projector-Augmented Wave method")')
        end select
        call glean(ao)
      end subroutine

      subroutine diary_angular_mesh_ao(ao)
!doc$ subroutine diary_angular_mesh(ao)
        type(atomic_operators_obj) :: ao
!       effects: Writes information about the angular mesh to the diary.

!cod$
        select case (ao%type)
        case (NCP)
          call glean(ao)
        case (PAW)
          call diary_angular_mesh(ao%ao_paw)
        end select
      end subroutine

      subroutine diary_occupation_ao(ao)
!doc$ subroutine diary_occupation(ao)
        type(atomic_operators_obj) :: ao
!       effects: Writes information about the initial atomic occupation to the diary.

!cod$
        select case (ao%type)
        case (NCP)
          call glean(ao)
        case (PAW)
          call diary_occupation(ao%ao_paw)
        end select
      end subroutine

      subroutine diary_rs_projectors_ao(ao)
!doc$ subroutine diary_rs_projectors(ao)
        type(atomic_operators_obj) :: ao
!       requires: Real-space projectors be in use.
!       effects: Writes real-space projector information to the diary.

!cod$
        select case (ao%type)
        case (NCP)
          call diary_rs_projectors(ao%ao_ncp)
        case (PAW)
          call diary_rs_projectors(ao%ao_paw)
        end select
      end subroutine

      subroutine write_restart_ao(ao,nrestf)
!doc$ subroutine write_restart(ao,nrestf)
        type(atomic_operators_obj) :: ao
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ao restart information to nrestf.

!cod$
        if (i_access(nrestf)) then

          ! start the ATOMIC_TYPE block
          call startblock(nrestf,"ATOMIC_TYPE")

          select case (ao%type)
          case (NCP)
            call writetag(nrestf,"NCP")
          case (PAW)
            call writetag(nrestf,"PAW")
          end select

          ! end the ATOMIC_TYPE block
          if (i_access(nrestf)) call endblock(nrestf)

        end if

        select case (ao%type)
        case (NCP)
          call write_restart(ao%ao_ncp,nrestf)
        case (PAW)
          call write_restart(ao%ao_paw,nrestf)
        end select

        if (error("Exit atomic_operators_mod::write_restart_ao")) continue

      end subroutine

      end module

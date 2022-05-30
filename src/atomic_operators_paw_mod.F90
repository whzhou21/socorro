! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module atomic_operators_paw_mod
!doc$ module atomic_operators_paw_mod

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use arg_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use lattice_mod
      use atoms_mod
      use crystal_mod
      use math_mod
      use symmetry_mod
      use paw_data_mod
      use xc_type_mod
      use axc_mod
      use diary_mod

!     One datatype is available here: type(atomic_operators_paw_obj)

!cod$
      implicit none
      private

      type :: gp_mat
        complex(double), dimension(:,:), pointer :: mat
      end type

      type :: type_data
        type(paw_data_obj) :: pawd                                ! paw data object
        character(tag_sz) :: tag                                  ! tag
        real(double) :: valence                                   ! number of valence electrons
        real(double) :: matching_radius                           ! matching radius
        integer, dimension(:), pointer :: pn                      ! projector basis index
        integer, dimension(:), pointer :: pl                      ! projector l values
        integer, dimension(:), pointer :: pm                      ! projector m values
        integer, dimension(:), pointer :: bs                      ! basis start indexes
        integer, dimension(:), pointer :: bl                      ! basis l values
        real(double), dimension(:), pointer :: bo                 ! basis initial occupations
        real(double), dimension(:,:,:), pointer :: gpff           ! grid potential form factors
        real(double), dimension(:,:,:), pointer :: cdff           ! coretail density form factors
        real(double), dimension(:,:,:,:), pointer :: hdff         ! hat density form factors
        real(double), dimension(:,:,:,:), pointer :: dmff         ! density matrix form factors
        complex(double), dimension(:,:), pointer :: r2c           ! real_to_complex Ylm transformation
        complex(double), dimension(:,:), pointer :: c2r           ! complex_to_real Ylm transformation
        complex(double), dimension(:,:), pointer :: tr2c          ! transpose of real_to_complex Ylm transformation
        complex(double), dimension(:,:), pointer :: tc2r          ! transpose of complex_to_real Ylm transformation
        integer :: theta_points                                   ! number of angular points in the theta dimension
        integer :: phi_points                                     ! number of angular points in the phi dimension
        real(double), dimension(:), pointer :: v_weight           ! weights for angular integration mesh
        complex(double), dimension(:,:,:), pointer :: ylmij       ! products of pairs of spherical harmonics
        complex(double), dimension(:,:,:), pointer :: dylmij_dt   ! derivative of ylmij with respect to angle theta
        complex(double), dimension(:,:,:), pointer :: dylmij_dp   ! derivative of ylmij with respect to angle phi
      end type

      type :: atom_data
        integer :: ti                                       ! type index
        integer :: bi                                       ! base index in projector list
        integer :: valence                                  ! valence
        complex(double), dimension(:,:), pointer :: c_oij   ! overlap matrix for complex-valued spherical harmonics
        complex(double), dimension(:,:), pointer :: r_oij   ! overlap matrix for real-valued spherical harmonics
      end type

      type :: projector_data
        integer :: ai    ! atom index
        integer :: ti    ! type index
        integer :: tpi   ! type projector index
        integer :: n     ! basis index
        integer :: l     ! l value
        integer :: m     ! m value
      end type

      type :: unique_atom_data
        integer :: ai    ! atom index
        integer :: ti    ! type index
      end type

      type :: vector_data
        integer :: ai    ! atom index
        integer :: ui    ! unique atom index
        integer :: ti    ! type index
        integer :: tvi   ! type vector index
      end type

      type :: atomic_operators_paw_rep
        integer :: ref
        type(ghost) :: g
        type(ghost) :: g_atoms                                   ! atoms ghost
        type(ghost) :: g_layout                                  ! layout ghost
        type(crystal_obj) :: cr                                  ! crystal object
        type(space_group_obj) :: sg                              ! space group object
        type(axc_obj) :: axc                                     ! atomic exchange-correlation object
        type(type_data), dimension(:), pointer :: tdata          ! type data
        type(atom_data), dimension(:), pointer :: adata          ! atom data
        type(projector_data), dimension(:), pointer :: pdata     ! projector data
        type(unique_atom_data), dimension(:), pointer :: udata   ! vector data
        type(vector_data), dimension(:), pointer :: vdata        ! vector data
        type(grid_obj) :: ctd                                    ! coretail density
      end type

      type, public :: atomic_operators_paw_obj
        private
        integer :: ref
        type(atomic_operators_paw_rep), pointer :: o
      end type

!doc$
      public :: atomic_operators_paw
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_layout_ghost
      public :: x_crystal
      public :: x_space_group
      public :: x_axc
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
      public :: x_atom_type
      public :: x_atom_base
      public :: x_atom_valence
      public :: x_n_basis
      public :: x_basis_start
      public :: x_basis_offset
      public :: x_basis_l
      public :: x_basis_occupation
      public :: x_type_l_max
      public :: x_atom_l_max
      public :: projector_f_values
      public :: projector_r_value
      public :: projector_r_gradients
      public :: x_projector_radius
      public :: x_type_projector_radius
      public :: x_n_unique_atoms
      public :: x_unique_atom
      public :: x_unique_type
      public :: x_n_vectors
      public :: x_vector_atom
      public :: x_vector_unique_atom
      public :: x_vector_type
      public :: x_vector_weight
      public :: atomic_grid_potential
      public :: atomic_xc_density
      public :: type_projector_data
      public :: type_r2c
      public :: type_c2r
      public :: type_tr2c
      public :: type_tc2r
      public :: atom_c_oij
      public :: atom_r_oij
      public :: vector_ylmij
      public :: vector_dylmij_dt
      public :: vector_dylmij_dp
      public :: diary_rs_projectors
      public :: diary_angular_mesh
      public :: diary_occupation
      public :: add_coretail_density
      public :: grid_potential_ff
      public :: coretail_density_ff
      public :: hat_density_ff
      public :: density_matrix_ff
      public :: write_restart

!cod$
      interface atomic_operators_paw
        module procedure constructor_ao
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
      interface x_layout_ghost
        module procedure ao_layout_ghost
      end interface
      interface x_crystal
        module procedure ao_crystal
      end interface
      interface x_space_group
        module procedure ao_space_group
      end interface
      interface x_axc
        module procedure ao_axc
      end interface
      interface x_n_types
        module procedure ao_n_types
      end interface
      interface x_type_name
        module procedure ao_type_name
      end interface
      interface x_type_valence
        module procedure ao_type_valence
      end interface
      interface x_n_atoms
        module procedure ao_n_atoms
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
      interface x_projector_atom
        module procedure ao_projector_atom
      end interface
      interface x_projector_index_in_type
        module procedure ao_projector_index_in_type
      end interface
      interface x_atom_type
        module procedure ao_atom_type
      end interface
      interface x_atom_base
        module procedure ao_atom_base
      end interface
      interface x_atom_valence
        module procedure ao_atom_valence
      end interface
      interface x_n_basis
        module procedure ao_n_basis
      end interface
      interface x_basis_start
        module procedure ao_basis_start
      end interface
      interface x_basis_offset
        module procedure ao_basis_offset
      end interface
      interface x_basis_l
        module procedure ao_basis_l
      end interface
      interface x_basis_occupation
        module procedure ao_basis_occupation
      end interface
      interface x_type_l_max
        module procedure ao_type_l_max
      end interface
      interface x_atom_l_max
        module procedure ao_atom_l_max
      end interface
      interface projector_f_values
        module procedure ao_projector_f_values, ao_type_projector_f_values
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
      interface x_n_unique_atoms
        module procedure ao_n_unique_atoms
      end interface
      interface x_unique_atom
        module procedure ao_unique_atom
      end interface
      interface x_unique_type
        module procedure ao_unique_type
      end interface
      interface x_n_vectors
        module procedure ao_n_vectors
      end interface
      interface x_vector_atom
        module procedure ao_vector_atom
      end interface
      interface x_vector_unique_atom
        module procedure ao_vector_unique_atom
      end interface
      interface x_vector_type
        module procedure ao_vector_type
      end interface
      interface x_vector_weight
        module procedure ao_vector_weight
      end interface
      interface atomic_grid_potential
        module procedure atomic_grid_potential_ao
      end interface
      interface atomic_xc_density
        module procedure atomic_xc_density_ao
      end interface
      interface diary_rs_projectors
        module procedure diary_rs_projectors_ao
      end interface
      interface diary_angular_mesh
        module procedure diary_angular_mesh_ao
      end interface
      interface diary_occupation
        module procedure diary_occupation_ao
      end interface
      interface add_coretail_density
        module procedure add_coretail_density_ao
      end interface
      interface write_restart
        module procedure write_restart_ao
      end interface

      contains

! public routines

      function constructor_ao(cr,lay,sg,xct,restf) result(ao)
!doc$ function atomic_operators_paw(cr,lay,sg,xct,restf) result(ao)
        type(crystal_obj) :: cr
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_type_obj) :: xct
        type(tagio_obj), optional :: restf
        type(atomic_operators_paw_obj) :: ao 
!       effects: Constructs a new ao.
!       errors: Passes errors.

!cod$
        character(1) :: tios

        call my(cr)
        call my(lay)
        call my(sg)
        call my(xct)
        if (present(restf)) call my(restf)

        ao%ref = 0
        allocate( ao%o )
        ao%o%ref = 0
        ao%o%g = x_ghost()
        ao%o%g_atoms = x_ghost(x_atoms(cr))
        ao%o%g_layout = x_ghost(lay)

        call my(cr,ao%o%cr)
        call my(sg,ao%o%sg)
        call my(axc(xct),ao%o%axc)

        ! open the ATOMIC_OPERATORS block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_OPERATORS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_OPERATORS block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)
        end if

        if (present(restf)) then
          call form_tdata_i(ao%o,lay,restf) ; if (error()) goto 100
        else
          call form_tdata_i(ao%o,lay)       ; if (error()) goto 100
        end if
        call form_adata_i(ao%o)
        call form_pdata_i(ao%o)
        call form_udata_i(ao%o) ; if (error()) goto 100
        call form_vdata_i(ao%o)

        call my(grid(lay,SGROUP),ao%o%ctd)
        call form_coretail_density_i(ao%o)

        ! close the ATOMIC_OPERATORS block
100     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

200     call glean(thy(cr))
        call glean(thy(lay))
        call glean(thy(sg))
        call glean(thy(xct))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_operators_paw_mod::constructor_ao")) continue

      end function 

      subroutine update_ao(ao,cr,lay,sg)
!doc$ subroutine update(ao,cr,lay,sg)
        type(atomic_operators_paw_obj) :: ao 
        type(crystal_obj) :: cr
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
!       requires: Atom types in cr not change.
!       modifies: ao
!       effects: Updates ao.
!       errors: Passes errors.

!cod$
        logical :: atoms_change, crystal_change, layout_change
 
        call my(ao)
        call my(cr)
        call my(lay)
        call my(sg)
 
        atoms_change = ( ao%o%g_atoms /= x_ghost(x_atoms(cr)) )
        crystal_change = ( x_ghost(ao%o%cr) /= x_ghost(cr) )
        layout_change = ( ao%o%g_layout /= x_ghost(lay) )

        if (layout_change .or. crystal_change) then
           call own_i(ao)
           ao%o%g = x_ghost()
           ao%o%g_atoms = x_ghost(x_atoms(cr))
           ao%o%g_layout = x_ghost(lay)
           ao%o%sg = sg
        end if

        if (crystal_change) ao%o%cr = cr

        if (layout_change .or. (crystal_change .and. .not.atoms_change)) then
          call form_factors_i(ao%o,lay) ; if (error()) goto 100
        end if

        if (atoms_change) then
          ao%o%ctd = grid(lay,SGROUP)
          call form_coretail_density_i(ao%o)
        end if

100     call glean(thy(ao))
        call glean(thy(cr))
        call glean(thy(lay))
        call glean(thy(sg))
 
        if (error("Exit atomic_operators_paw_mod::update_ao")) continue
 
      end subroutine

      subroutine my_ao(ao)
!doc$ subroutine my(ao)
        type(atomic_operators_paw_obj) :: ao

!cod$
        ao%ref = ao%ref + 1
        ao%o%ref = ao%o%ref + 1
      end subroutine
 
      subroutine my_new_ao(aoi,ao)
!doc$ subroutine my(aoi,ao)
        type(atomic_operators_paw_obj) :: aoi, ao

!cod$
        ao%ref = 1
        ao%o => aoi%o
        ao%o%ref = ao%o%ref + 1
      end subroutine

      function thy_ao(ao) result(aoo)
!doc$ function thy(ao) result(aoo)
        type(atomic_operators_paw_obj) :: ao, aoo

!cod$
        ao%ref = ao%ref - 1 
        ao%o%ref = ao%o%ref - 1
        aoo%ref = ao%ref
        aoo%o => ao%o
      end function
 
      subroutine glean_ao(ao)
!doc$ subroutine glean(ao)
        type(atomic_operators_paw_obj) :: ao

!cod$
        integer :: i
        if (ao%o%ref < 1) then
          call glean(thy(ao%o%cr))
          call glean(thy(ao%o%sg))
          call glean(thy(ao%o%axc))
          if (associated( ao%o%tdata )) then
            do i = 1,size(ao%o%tdata)
              call glean(thy(ao%o%tdata(i)%pawd))
              if (associated( ao%o%tdata(i)%pn )) deallocate( ao%o%tdata(i)%pn )
              if (associated( ao%o%tdata(i)%pl )) deallocate( ao%o%tdata(i)%pl )
              if (associated( ao%o%tdata(i)%pm )) deallocate( ao%o%tdata(i)%pm )
              if (associated( ao%o%tdata(i)%bs )) deallocate( ao%o%tdata(i)%bs )
              if (associated( ao%o%tdata(i)%bl )) deallocate( ao%o%tdata(i)%bl )
              if (associated( ao%o%tdata(i)%bo )) deallocate( ao%o%tdata(i)%bo )
              if (associated( ao%o%tdata(i)%gpff )) deallocate( ao%o%tdata(i)%gpff )
              if (associated( ao%o%tdata(i)%cdff )) deallocate( ao%o%tdata(i)%cdff )
              if (associated( ao%o%tdata(i)%hdff )) deallocate( ao%o%tdata(i)%hdff )
              if (associated( ao%o%tdata(i)%dmff )) deallocate( ao%o%tdata(i)%dmff )
              if (associated( ao%o%tdata(i)%r2c )) deallocate( ao%o%tdata(i)%r2c )
              if (associated( ao%o%tdata(i)%c2r )) deallocate( ao%o%tdata(i)%c2r )
              if (associated( ao%o%tdata(i)%tr2c )) deallocate( ao%o%tdata(i)%tr2c )
              if (associated( ao%o%tdata(i)%tc2r )) deallocate( ao%o%tdata(i)%tc2r )
              if (associated( ao%o%tdata(i)%v_weight )) deallocate( ao%o%tdata(i)%v_weight )
              if (associated( ao%o%tdata(i)%ylmij )) deallocate( ao%o%tdata(i)%ylmij )
              if (associated( ao%o%tdata(i)%dylmij_dt )) deallocate( ao%o%tdata(i)%dylmij_dt )
              if (associated( ao%o%tdata(i)%dylmij_dp )) deallocate( ao%o%tdata(i)%dylmij_dp )
            end do
            deallocate( ao%o%tdata )
          end if
          if (associated( ao%o%adata )) then
            do i = 1,size(ao%o%adata)
              if (associated( ao%o%adata(i)%c_oij )) deallocate( ao%o%adata(i)%c_oij )
              if (associated( ao%o%adata(i)%r_oij )) deallocate( ao%o%adata(i)%r_oij )
            end do
          end if
          deallocate( ao%o%adata )
          if (associated( ao%o%pdata )) deallocate( ao%o%pdata )
          if (associated( ao%o%udata )) deallocate( ao%o%udata )
          if (associated( ao%o%vdata )) deallocate( ao%o%vdata )
          call glean(thy(ao%o%ctd))
          deallocate( ao%o )
        end if

      end subroutine

      subroutine bequeath_ao(ao)
!doc$ subroutine bequeath(ao)
        type(atomic_operators_paw_obj) :: ao

!cod$
        continue
      end subroutine
 
      subroutine assign_ao(ao,ao2)
!doc$ subroutine assignment(=)(ao,ao2)
        type(atomic_operators_paw_obj), intent(inout) :: ao
        type(atomic_operators_paw_obj), intent(in) :: ao2

!cod$
        type(atomic_operators_paw_obj) :: aot
        call my(ao2)
        aot%o => ao%o
        ao%o%ref = ao%o%ref - ao%ref
        ao%o => ao2%o
        ao%o%ref = ao%o%ref + ao%ref
        call glean(aot)
        call glean(thy(ao2))
      end subroutine
 
      function ao_ref(ao) result(r)
!doc$ function x_ref(ao) result(r)
        type(atomic_operators_paw_obj) :: ao
        integer, dimension(2) :: r
!       effects: Returns ao%ref and ao%o%ref.

!cod$
        r(1) = ao%ref
        r(2) = ao%o%ref
        call glean(ao)
      end function 

      function ao_ghost(ao) result(g)
!doc$ function x_ghost(ao) result(g)
        type(atomic_operators_paw_obj) :: ao
        type(ghost) :: g
!       effects: Returns the ghost of ao.

!cod$ 
        call my(ao)
        g = ao%o%g
        call glean(thy(ao))
      end function 

      function ao_layout_ghost(ao) result(g)
!doc$ function x_layout_ghost(ao) result(g)
        type(atomic_operators_paw_obj) :: ao
        type(ghost) :: g
!       effects: Returns the layout ghost of ao.

!cod$ 
        call my(ao)
        g = ao%o%g_layout
        call glean(thy(ao))
      end function 

      function ao_crystal(ao) result(cr)
!doc$ function x_crystal(ao) result(cr)
        type(atomic_operators_paw_obj) :: ao
        type(crystal_obj) :: cr
!       effects: Returns ao%o%cr.

!cod$ 
        call my(ao)
        call my(ao%o%cr,cr)
        call bequeath(thy(cr))
        call glean(thy(ao))
      end function 

      function ao_space_group(ao) result(sg)
!doc$ function x_space_group(ao) result(sg)
        type(atomic_operators_paw_obj) :: ao
        type(space_group_obj) :: sg
!       effects: Returns ao%o%sg.

!cod$ 
        call my(ao)
        call my(ao%o%sg,sg)
        call bequeath(thy(sg))
        call glean(thy(ao))
      end function 

      function ao_axc(ao) result(axc)
!doc$ function x_axc(ao) result(axc)
        type(atomic_operators_paw_obj) :: ao
        type(axc_obj) :: axc
!       effects: Returns ao%o%axc.

!cod$ 
        call my(ao)
        call my(ao%o%axc,axc)
        call bequeath(thy(axc))
        call glean(thy(ao))
      end function 

      function ao_n_types(ao) result(nt)
!doc$ function x_n_types(ao) result(nt)
        type(atomic_operators_paw_obj) :: ao
        integer :: nt
!       effects: Returns the number of atom types.

!cod$
        call my(ao)
        nt = size(ao%o%tdata)
        call glean(thy(ao))
      end function

      function ao_type_name(ao,it) result(tag)
!doc$ function x_type_name(ao,it) result(tag)
        type(atomic_operators_paw_obj) :: ao
        integer :: it
        character(tag_sz) :: tag
!       effects: Returns the name of type it.

!cod$
        call my(ao)
        tag = ao%o%tdata(it)%tag
        call glean(thy(ao))
      end function

      function ao_type_valence(ao,it) result(v)
!doc$ function x_type_valence(ao,it) result(v)
        type(atomic_operators_paw_obj) :: ao
        integer :: it
        real(double) :: v
!       effects: Returns the valence of type it.

!cod$
        call my(ao)
        v = ao%o%tdata(it)%valence
        call glean(thy(ao))
      end function

      function ao_n_atoms(ao) result(na)
!doc$ function x_n_atoms(ao) result(na)
        type(atomic_operators_paw_obj) :: ao
        integer :: na
!       effects: Returns the number of atoms in ao.
!       errors: Passes errors.

!cod$
        call my(ao)
        na = size(ao%o%adata)
        call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_n_atoms")) continue
      end function

      function ao_valence_electrons(ao) result(ve)
!doc$ function x_valence_electrons(ao) result(ve)
        type(atomic_operators_paw_obj) :: ao
        real(double) :: ve
!       effects: Returns the sum of atom valences.

!cod$
        integer :: ia
        call my(ao)
        ve = 0.0_double
        do ia = 1,size(ao%o%adata)
          ve = ve + ao%o%adata(ia)%valence
        end do
        call glean(thy(ao))
      end function

      function ao_type_matching_radius(ao,it) result(mr)
!doc$ function x_type_matching_radius(ao,it) result(mr)
        type(atomic_operators_paw_obj) :: ao
        integer :: it
        real(double) :: mr
!       effects: Returns the matching radius of atom type it.

!cod$
        call my(ao)
        mr = ao%o%tdata(it)%matching_radius
        call glean(thy(ao))
      end function

      function ao_atom_matching_radius(ao,ia) result(mr)
!doc$ function x_atom_matching_radius(ao,ia) result(mr)
        type(atomic_operators_paw_obj) :: ao
        integer :: ia
        real(double) :: mr
!       effects: Returns the matching radius of atom ia.

!cod$
        call my(ao)
        mr = ao%o%tdata(ao%o%adata(ia)%ti)%matching_radius
        call glean(thy(ao))
      end function

      function ao_n_projectors(ao) result(np)
!doc$ function x_n_projectors(ao) result(np)
        type(atomic_operators_paw_obj) :: ao 
        integer :: np
!       effects: Returns the number of projectors.

!cod$ 
        call my(ao)
        if (associated( ao%o%pdata )) then
          np = size(ao%o%pdata)
        else
          np = 0
        end if
        call glean(thy(ao))
      end function

      function ao_n_type_projectors(ao,it) result(np)
!doc$ function x_n_type_projectors(ao,it) result(np)
        type(atomic_operators_paw_obj) :: ao 
        integer, intent(in) :: it
        integer :: np
!       effects: Returns the number of projectors for the type with index it.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error( ((it < 1) .or. (it > size(ao%o%tdata))),"ERROR: it is out of range")) goto 100
        if (associated( ao%o%tdata(it)%pn )) then
          np = size(ao%o%tdata(it)%pn)
        else
          np = 0
        end if
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_n_type_projectors")) continue
      end function

      function ao_n_atom_projectors(ao,ia) result(np)
!doc$ function x_n_atom_projectors(ao,ia) result(np)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        integer :: np
!       effects: Returns the number of projectors for the atom with index ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error( ((ia < 1) .or. (ia > size(ao%o%adata))),"ERROR: ia is out of range")) goto 100
        if (associated( ao%o%tdata(ao%o%adata(ia)%ti)%pn )) then
          np = size(ao%o%tdata(ao%o%adata(ia)%ti)%pn)
        else
          np = 0
        end if
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_n_atom_projectors")) continue
      end function

      function ao_projector_l(ao,ip) result(l)
!doc$ function x_projector_l(ao,ip) result(l)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        integer :: l
!       requires: At least one projector exist.
!       effects: Returns the l angular momentum quantum number of projector ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        l = ao%o%pdata(ip)%l
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_l")) continue
      end function

      function ao_type_projector_l(ao,it,ip) result(l)
!doc$ function x_projector_l(ao,it,ip) result(l)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it, ip
        integer :: l
!       requires: At least one projector exist.
!       effects: Returns the l angular momentum quantum number of projector ip of type it.
!       errors: it or ip out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pl)),"ERROR: ip is out of range")) goto 100
        l = ao%o%tdata(it)%pl(ip)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_type_projector_l")) continue
      end function

      function ao_projector_m(ao,ip) result(m)
!doc$ function x_projector_m(ao,ip) result(m)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        integer :: m
!       requires: At least one projector exist.
!       effects: Returns the m angular momentum quantum number of projector ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        m = ao%o%pdata(ip)%m
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_m")) continue
      end function

      function ao_type_projector_m(ao,it,ip) result(m)
!doc$ function x_projector_m(ao,it,ip) result(m)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it, ip
        integer :: m
!       requires: At least one projector exist.
!       effects: Returns the l angular momentum quantum number of projector ip of type it.
!       errors: it or ip out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pm)),"ERROR: ip is out of range")) goto 100
        m = ao%o%tdata(it)%pm(ip)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_type_projector_m")) continue
      end function

      function ao_projector_type(ao,ip) result(ti)
!doc$ function x_projector_type(ao,ip) result(ti)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        integer :: ti
!       requires: At least one projector exist.
!       effects:  Returns the type index of the projector with global index ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        ti = ao%o%pdata(ip)%ti
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_type")) continue
      end function

      function ao_projector_atom(ao,ip) result(ai)
!doc$ function x_projector_atom(ao,ip) result(ai)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        integer :: ai
!       requires: At least one projector exist.
!       effects: Returns the atom index of the projector with global index ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        ai = ao%o%pdata(ip)%ai
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_atom")) continue
      end function

      function ao_projector_index_in_type(ao,ip) result(tpi)
!doc$ function x_projector_index_in_type(ao,ip) result(tpi)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        integer :: tpi
!       requires: At least one projector exist.
!       effects:  Returns the type index of the projector with global index ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        tpi = ao%o%pdata(ip)%tpi
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_index_in_type")) continue
      end function

      function ao_atom_type(ao,ia) result(ti)
!doc$ function x_atom_type(ao,ia) result(ti)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        integer :: ti
!       effects: Returns the type index of atom ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        ti = ao%o%adata(ia)%ti
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_atom_type")) continue
      end function

      function ao_atom_base(ao,ia) result(bi)
!doc$ function x_atom_base(ao,ia) result(bi)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        integer :: bi
!       requires: At least one projector exist for atom ia.
!       effects: Returns the base index of atom ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        bi = ao%o%adata(ia)%bi
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_atom_base")) continue
      end function

      function ao_atom_valence(ao,ia) result(v)
!doc$ function x_atom_valence(ao,ia) result(v)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        real(double) :: v
!       effects: Returns the valence of atom ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        v = ao%o%adata(ia)%valence
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_atom_valence")) continue
      end function

      function ao_n_basis(ao,it) result(n)
!doc$ function x_n_basis(ao,it) result(n)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        integer :: n
!       effects: Returns the number of basis functions of type it.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        n = size(ao%o%tdata(it)%bs)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_n_basis")) continue
      end function

      function ao_basis_start(ao,it,ib) result(s)
!doc$ function x_basis_start(ao,it,ib) result(s)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it, ib
        integer :: s
!       effects: Returns the starting index of the ib basis function of type it.
!       errors: it or ib out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(ao%o%tdata(it)%bs)),"ERROR: ib is out of range")) goto 100
        s = ao%o%tdata(it)%bs(ib)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_basis_start")) continue
      end function

      function ao_basis_offset(ao,it,ib) result(o)
!doc$ function x_basis_offset(ao,it,ib) result(o)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it, ib
        integer :: o
!       effects: Returns the offset of the ib basis function of type it.
!       errors: it or ib out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(ao%o%tdata(it)%bs)),"ERROR: ib is out of range")) goto 100
        o = ao%o%tdata(it)%bs(ib) + ao%o%tdata(it)%bl(ib)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_basis_start")) continue
      end function

      function ao_basis_l(ao,it,ib) result(l)
!doc$ function x_basis_l(ao,it,ib) result(l)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it, ib
        integer :: l
!       effects: Returns the l value of the ib basis function of type it.
!       errors: it or ib out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(ao%o%tdata(it)%bs)),"ERROR: ib is out of range")) goto 100
        l = ao%o%tdata(it)%bl(ib)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_basis_l")) continue
      end function

      function ao_basis_occupation(ao,it,ib) result(o)
!doc$ function x_basis_occupation(ao,it,ib) result(o)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it, ib
        real(double) :: o
!       effects: Returns the occupation of the ib basis function of type it.
!       errors: it or ib out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(ao%o%tdata(it)%bs)),"ERROR: ib is out of range")) goto 100
        o = ao%o%tdata(it)%bo(ib)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_basis_occupation")) continue
      end function

      function ao_type_l_max(ao,it) result(l)
!doc$ function x_type_l_max(ao,it) result(l)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        integer :: l
!       requires: At least one type projector exist.
!       effects: Returns the maximum l value among the type it projectors.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        l = maxval(ao%o%tdata(it)%pl)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_type_l_max")) continue
      end function

      function ao_atom_l_max(ao,ia) result(l)
!doc$ function x_atom_l_max(ao,ia) result(l)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        integer :: l
!       requires: At least one atom projector exist.
!       effects: Returns the maximum l value among the atom ia projectors.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        l = maxval(ao%o%tdata(ao%o%adata(ia)%ti)%pl)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_atom_l_max")) continue
      end function

      function ao_projector_f_values(ao,q,ip) result(pfv)
!doc$ function projector_f_values(ao,q,ip) result(pfv)
        type(atomic_operators_paw_obj) :: ao
        real(double), dimension(:), intent(in) :: q
        integer, intent(in) :: ip
        real(double), dimension(size(q)) :: pfv
!       requires: At least one projector exist.
!       effects: Returns the q Fourier coefficients of the projector with global index ip.
!       errors: ip out of range. Passes errors.

!cod$
        integer :: n, ti

        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        ti = ao%o%pdata(ip)%ti
        n = ao%o%pdata(ip)%n
        pfv = projector_f_values(ao%o%tdata(ti)%pawd,q,n)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_f_values")) continue

      end function 

      function ao_type_projector_f_values(ao,q,it,ip) result(pfv)
!doc$ function projector_f_values(ao,q,it,ip) result(pfv)
        type(atomic_operators_paw_obj) :: ao
        real(double), dimension(:), intent(in) :: q
        integer, intent(in) :: it, ip
        real(double), dimension(size(q)) :: pfv
!       requires: At least one projector exist.
!       effects: Returns the q Fourier coefficients of the projector with with type it and type index ip.
!       errors: it or ip out of range. Passes errors.

!cod$
        integer :: ni

        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pn)),"ERROR: ip is out of range")) goto 100
        ni = ao%o%tdata(it)%pn(ip)
        pfv = projector_f_values(ao%o%tdata(it)%pawd,q,ni)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_type_projector_f_values")) continue

      end function

      function ao_projector_r_value(ao,ip,r,gi,go) result(prv)
!doc$ function projector_r_value(ao,ip,r,gi,go) result(prv)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        real(double), intent(in) :: r
        real(double), intent(in) :: gi, go
        real(double) :: prv
!       requires: At least one projector exist.
!       effects: Returns the value of the ip optimized projector at r.
!       errors: ip out of range. Passes errors.

!cod$
        integer :: n, ti

        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        ti = ao%o%pdata(ip)%ti
        n = ao%o%pdata(ip)%n
        prv = projector_r_value(ao%o%tdata(ti)%pawd,r,n,gi,go) ; if (error()) goto 100
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_r_value")) continue

      end function

      function ao_projector_r_gradients(ao,ip,r) result(prg)
!doc$ function projector_r_gradients(ao,ip,r) result(prg)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        real(double), intent(in) :: r
        real(double), dimension(2) :: prg
!       requires: At least one projector exist.
!       effects: Returns the components of the gradient of the ip optimized projector at r.
!       errors: ip out of range. Passes errors.

!cod$
        integer :: n, ti

        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        ti = ao%o%pdata(ip)%ti
        n = ao%o%pdata(ip)%n
        prg = projector_r_gradients(ao%o%tdata(ti)%pawd,r,n) ; if (error()) goto 100
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_r_gradients")) continue

      end function

      function ao_projector_radius(ao,ip) result(r)
!doc$ function x_projector_radius(ao,ip) result(r)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ip
        real(double) :: r
!       requires: At least one projector exist.
!       effects: Returns the optimization radius of the projector with global index ip.
!       errors: ip out of range. Passes errors.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        r = x_radius(ao%o%tdata(ao%o%pdata(ip)%ti)%pawd)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_projector_radius")) continue
      end function

      function ao_type_projector_radius(ao,it) result(r)
!doc$ function x_type_projector_radius(ao,it) result(r)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        real(double) :: r
!       effects: Returns the optimization radius of type it projectors.
!       errors: it out of range. Passes errors.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        r = x_radius(ao%o%tdata(it)%pawd)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_type_projector_radius")) continue
      end function

      function ao_n_unique_atoms(ao) result(n)
!doc$ function x_n_unique_atoms(ao) result(n)
        type(atomic_operators_paw_obj) :: ao
        integer :: n
!       effects: Returns the number of unique atoms.

!cod$
        call my(ao)
        n = size(ao%o%udata)
        call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_n_unique_atoms")) continue
      end function

      function ao_unique_atom(ao,iu) result(i)
!doc$ function x_unique_atom(ao,iu) result(i)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iu
        integer :: i
!       effects: Returns the atom index for unique atom iu.
!       errors: iu out of range.

!cod$
        call my(ao)
        if (error((iu < 1) .or. (iu > size(ao%o%udata)),"ERROR: iu is out of range")) goto 100
        i = ao%o%udata(iu)%ai
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_unique_atom")) continue
      end function

      function ao_unique_type(ao,iu) result(i)
!doc$ function x_unique_type(ao,iu) result(i)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iu
        integer :: i
!       effects: Returns the type index for unique atom iu.
!       errors: iu out of range.

!cod$
        call my(ao)
        if (error((iu < 1) .or. (iu > size(ao%o%udata)),"ERROR: iu is out of range")) goto 100
        i = ao%o%udata(iu)%ti
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_unique_type")) continue
      end function

      function ao_n_vectors(ao) result(n)
!doc$ function x_n_vectors(ao) result(n)
        type(atomic_operators_paw_obj) :: ao
        integer :: n
!       effects: Returns the number of angular integration vectors.

!cod$
        call my(ao)
        n = size(ao%o%vdata)
        call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_n_vectors")) continue
      end function

      function ao_vector_atom(ao,iv) result(i)
!doc$ function x_vector_atom(ao,iv) result(i)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        integer :: i
!       effects: Returns the atom index for angular integration vector iv.
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        i = ao%o%vdata(iv)%ai
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_vector_atom")) continue
      end function

      function ao_vector_unique_atom(ao,iv) result(i)
!doc$ function x_vector_unique_atom(ao,iv) result(i)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        integer :: i
!       effects: Returns the unique atom index for angular integration vector iv.
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        i = ao%o%vdata(iv)%ui
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_vector_unique_atom")) continue
      end function

      function ao_vector_type(ao,iv) result(i)
!doc$ function x_vector_type(ao,iv) result(i)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        integer :: i
!       effects: Returns the type index for angular integration vector iv.
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        i = ao%o%vdata(iv)%ti
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_vector_type")) continue
      end function

      function ao_vector_weight(ao,iv) result(r)
!doc$ function x_vector_weight(ao,iv) result(r)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        real(double) :: r
!       effects: Returns ao%o%tdata(ao%o%vdata(iv)%ti)%v_weight(ao%o%vdata(iv)%tvi).
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        r = ao%o%tdata(ao%o%vdata(iv)%ti)%v_weight(ao%o%vdata(iv)%tvi)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::ao_vector_weight")) continue
      end function

      function atomic_grid_potential_ao(ao,lay) result(pot)
!doc$ function atomic_grid_potential(ao,lay) result(pot)
        type(atomic_operators_paw_obj) :: ao
        type(layout_obj) :: lay
        type(grid_obj) :: pot
!       effects: Returns the filtered atomic grid potential.
!       errors: Inconsistent layouts.

!cod$
        integer :: i1, i2, i3, ia, it
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        complex(double) :: igr
        complex(double), dimension(:,:,:), pointer :: c_pot
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(ao)
        call my(lay)

        call my(x_lattice(ao%o%cr),lat)
        call my(x_atoms(ao%o%cr),ats)
        call my(grid(lay,SGROUP),pot)

        nullify( gx, gy, gz, c_pot )

        if (error(x_ghost(lay) /= ao%o%g_layout,"ERROR: inconsistent layouts")) goto 100

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call alloc(c_pot,lay,D_TYPE,SGROUP)
        c_pot = (0.0_double,0.0_double)
        do ia = 1,x_n_atoms(ao)
          it = x_atom_type(ao,ia)
          pos = lat2r(lat,x_position(ats,ia))
          do i3 = 1,size(gx,3)
          do i2 = 1,size(gx,2)
          do i1 = 1,size(gx,1)
            if (ao%o%tdata(it)%gpff(i1,i2,i3) == 0.0_double) cycle
            igr = (0.0_double,1.0_double)*(pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3))
            c_pot(i1,i2,i3) = c_pot(i1,i2,i3) + ao%o%tdata(it)%gpff(i1,i2,i3)*exp(-igr)
          end do
          end do
          end do
        end do
        call put(c_pot,pot,CDF_KIND)

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_pot )) deallocate( c_pot )

        call glean(thy(lat))
        call glean(thy(ats))
        call bequeath(thy(pot))

        call glean(thy(ao))
        call glean(thy(lay))

        if (error("Exit atomic_operators_paw_mod::atomic_grid_potential_ao")) continue

      end function

      function atomic_xc_density_ao(ao,lay) result(den)
!doc$ function atomic_xc_density(ao,lay) result(den)
        type(atomic_operators_paw_obj) :: ao
        type(layout_obj) :: lay
        type(grid_obj) :: den
!       effects: Returns the filtered coretail density.

!cod$
        real(double) :: spin_factor
        real(double), dimension(:,:,:), pointer :: r1

        call my(ao)
        call my(lay)

        call my(ao%o%ctd,den)
        select case (mpi_nsgroups())
        case (2)
          call take(r1,den,RD_KIND)
          spin_factor = 1.0_double/real(mpi_nsgroups(),double)
          r1 = spin_factor*r1
          call put(r1,den,RD_KIND)
        end select

        call glean(thy(ao))
        call glean(thy(lay))

        call bequeath(thy(den))

      end function

      subroutine type_projector_data(ao,it,pawd)
!doc$ subroutine type_projector_data(ao,it,pawd)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        type(paw_data_obj), pointer :: pawd
!       modifies: pawd
!       effects: Points pawd to ao%o%tdata(it)%pawd.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        pawd => ao%o%tdata(it)%pawd
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::type_projector_data")) continue

      end subroutine

      subroutine type_r2c(ao,it,p)
!doc$ subroutine type_r2c(ao,it,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(it)%r2c.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        p => ao%o%tdata(it)%r2c
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::type_r2c")) continue

      end subroutine

      subroutine type_c2r(ao,it,p)
!doc$ subroutine type_c2r(ao,it,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(it)%c2r.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        p => ao%o%tdata(it)%c2r
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::type_c2r")) continue

      end subroutine

      subroutine type_tr2c(ao,it,p)
!doc$ subroutine type_tr2c(ao,it,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(it)%tr2c.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        p => ao%o%tdata(it)%tr2c
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::type_tr2c")) continue

      end subroutine

      subroutine type_tc2r(ao,it,p)
!doc$ subroutine type_tc2r(ao,it,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(it)%tc2r.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        p => ao%o%tdata(it)%tc2r
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::type_tc2r")) continue

      end subroutine

      subroutine atom_c_oij(ao,ia,p)
!doc$ subroutine atom_c_oij(ao,ia,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%adata(ia)%c_oij.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        p => ao%o%adata(ia)%c_oij
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::atom_c_oij")) continue

      end subroutine

      subroutine atom_r_oij(ao,ia,p)
!doc$ subroutine atom_r_oij(ao,ia,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: ia
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%adata(ia)%r_oij.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        p => ao%o%adata(ia)%r_oij
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::atom_r_oij")) continue

      end subroutine

      subroutine vector_ylmij(ao,iv,p)
!doc$ subroutine vector_ylmij(ao,iv,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(ao%o%vdata(iv)%ti)%ylmij(:,:,ao%o%vdata(iv)%tvi).
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        p => ao%o%tdata(ao%o%vdata(iv)%ti)%ylmij(:,:,ao%o%vdata(iv)%tvi)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::vector_ylmij")) continue
      end subroutine

      subroutine vector_dylmij_dt(ao,iv,p)
!doc$ subroutine vector_dylmij_dt(ao,iv,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(ao%o%vdata(iv)%ti)%dylmij_dt(:,:,ao%o%vdata(iv)%tvi).
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        p => ao%o%tdata(ao%o%vdata(iv)%ti)%dylmij_dt(:,:,ao%o%vdata(iv)%tvi)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::vector_dylmij_dt")) continue
      end subroutine

      subroutine vector_dylmij_dp(ao,iv,p)
!doc$ subroutine vector_dylmij_dp(ao,iv,p)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: iv
        complex(double), dimension(:,:), pointer :: p
!       modifies: p
!       effects: Points p to ao%o%tdata(ao%o%vdata(iv)%ti)%dylmij_dp(:,:,ao%o%vdata(iv)%tvi).
!       errors: iv out of range.

!cod$
        call my(ao)
        if (error((iv < 1) .or. (iv > size(ao%o%vdata)),"ERROR: iv is out of range")) goto 100
        p => ao%o%tdata(ao%o%vdata(iv)%ti)%dylmij_dp(:,:,ao%o%vdata(iv)%tvi)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_paw_mod::vector_dylmij_dp")) continue
      end subroutine

      subroutine diary_rs_projectors_ao(ao)
!doc$ subroutine diary_rs_projectors(ao)
        type(atomic_operators_paw_obj) :: ao
!       requires: Real-space projectors be in use.
!       effects: Writes real-space projector information to the diary.

!cod$
        integer :: it
        call my(ao)
        do it = 1,size(ao%o%tdata)
          call diary_rs_projectors(ao%o%tdata(it)%pawd)
        end do
        call glean(thy(ao))
      end subroutine

      subroutine diary_angular_mesh_ao(ao)
!doc$ subroutine diary_angular_mesh(ao)
        type(atomic_operators_paw_obj) :: ao
!       effects: Writes angular mesh information to the diary.

!cod$
        character(tag_sz) :: tag
        integer :: it

        call my(ao)
        if (i_access(diaryfile())) then
          if (size(ao%o%tdata) == 1) then
            write(x_unit(diaryfile()),'(/,t4,"Angular mesh for atomic exchange and correlation:")')
          else
            write(x_unit(diaryfile()),'(/,t4,"Angular meshes for atomic exchange and correlation:")')
          end if
          do it = 1,size(ao%o%tdata)
            tag = ao%o%tdata(it)%tag
            select case (len_trim(tag))
            case (1)
              write(x_unit(diaryfile()),'(/,t6,a1,":")') trim(tag)
            case (2)
              write(x_unit(diaryfile()),'(/,t6,a2,":")') trim(tag)
            case (3)
              write(x_unit(diaryfile()),'(/,t6,a3,":")') trim(tag)
            case (4)
              write(x_unit(diaryfile()),'(/,t6,a4,":")') trim(tag)
            case (5)
              write(x_unit(diaryfile()),'(/,t6,a5,":")') trim(tag)
            case default
              write(x_unit(diaryfile()),'(/,t6,a6,":")') trim(tag)
            end select
            if ((ao%o%tdata(it)%theta_points < 10) .and. (ao%o%tdata(it)%phi_points < 10)) then 
              write(x_unit(diaryfile()),'(t8,i1," points along theta")') ao%o%tdata(it)%theta_points
              write(x_unit(diaryfile()),'(t8,i1," points along phi")') ao%o%tdata(it)%phi_points
            else
              write(x_unit(diaryfile()),'(t8,i2," points along theta")') ao%o%tdata(it)%theta_points
              write(x_unit(diaryfile()),'(t8,i2," points along phi")') ao%o%tdata(it)%phi_points
            end if
          end do
        end if
        call glean(thy(ao))
      end subroutine

      subroutine diary_occupation_ao(ao)
!doc$ subroutine diary_occupation(ao)
        type(atomic_operators_paw_obj) :: ao
!       effects: Writes initial atomic occupation information to the diary.

!cod$
        character(tag_sz) :: tag
        integer :: ib, it, mb, msg, nsg, nt
        real(double), dimension(:,:,:), allocatable :: o1, o2

        call my(ao)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        nt = size(ao%o%tdata)

        ! find the maximum number of basis functions
        mb = 0
        do it = 1,nt
          mb = max(mb,size(ao%o%tdata(it)%bo))
        end do

        ! Gather all occupations to the rank 0 CONFIG process
        allocate( o1(nsg,nt,mb), o2(nsg,nt,mb) )
        o1 = 0.0_double
        do it = 1,nt
          do ib = 1,size(ao%o%tdata(it)%bo)
            o1(msg,it,ib) = ao%o%tdata(it)%bo(ib)
          end do
        end do
        call xcomm_reduce(XSGROUP,MPI_SUM,o1,o2) ; if (error()) goto 100

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,t4,"Initial atomic occupations:")')
          do it = 1,nt
            tag = ao%o%tdata(it)%tag
            select case (nsg)
            case (1)
              select case (len_trim(tag))
              case (1)
                write(x_unit(diaryfile()),'(/,t6,a1,":")') trim(tag)
              case (2)
                write(x_unit(diaryfile()),'(/,t6,a2,":")') trim(tag)
              case (3)
                write(x_unit(diaryfile()),'(/,t6,a3,":")') trim(tag)
              case (4)
                write(x_unit(diaryfile()),'(/,t6,a4,":")') trim(tag)
              case (5)
                write(x_unit(diaryfile()),'(/,t6,a5,":")') trim(tag)
              case default
                write(x_unit(diaryfile()),'(/,t6,a6,":")') trim(tag)
              end select
              do ib = 1,size(ao%o%tdata(it)%bl)
                write(x_unit(diaryfile()),'(t8,"l = ",i1,": ",f9.6)') ao%o%tdata(it)%bl(ib), o2(1,it,ib)
              end do
            case (2)
              select case (len_trim(tag))
              case (1)
                write(x_unit(diaryfile()),'(/,t6,a1,":",7x,"spin group 1",4x,"spin group 2")') trim(tag)
              case (2)
                write(x_unit(diaryfile()),'(/,t6,a2,":",7x,"spin group 1",4x,"spin group 2")') trim(tag)
              case (3)
                write(x_unit(diaryfile()),'(/,t6,a3,":",7x,"spin group 1",4x,"spin group 2")') trim(tag)
              case (4)
                write(x_unit(diaryfile()),'(/,t6,a4,":",7x,"spin group 1",4x,"spin group 2")') trim(tag)
              case (5)
                write(x_unit(diaryfile()),'(/,t6,a5,":",7x,"spin group 1",4x,"spin group 2")') trim(tag)
              case default
                write(x_unit(diaryfile()),'(/,t6,a6,":",7x,"spin group 1",4x,"spin group 2")') trim(tag)
              end select
              do ib = 1,size(ao%o%tdata(it)%bl)
                write(x_unit(diaryfile()),'(t8,"l = ",i1,":   ",f9.6,9x,f9.6)') ao%o%tdata(it)%bl(ib), o2(1,it,ib), o2(2,it,ib)
              end do
            end select
          end do
        end if

100     if (allocated( o1 )) deallocate( o1 )
        if (allocated( o2 )) deallocate( o2 )

        call glean(thy(ao))

      end subroutine

      subroutine add_coretail_density_ao(ao,den)
!doc$ subroutine add_coretail_density_ao(ao,den)
        type(atomic_operators_paw_obj) :: ao
        type(grid_obj) :: den
!       modifies: den
!       effect: Adds the coretail density to den.
!       errors: Passes errors.

!cod$
        call my(ao)
        call my(den)
        call saxpby(1.0_double,den,1.0_double,ao%o%ctd)
        call glean(thy(ao))
        call glean(thy(den))
        if (error("Exit atomic_operators_paw_mod::add_coretail_density_ao")) continue
      end subroutine

      subroutine grid_potential_ff(ao,it,gpff)
!doc$ subroutine grid_potential_ff(ao,it,gpff)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:), pointer :: gpff
!       requires: it in range.
!       modifies: gpff
!       effects: Points gpff to the grid potential form factors of type it.

!cod$
        call my(ao)
        gpff => ao%o%tdata(it)%gpff
        call glean(thy(ao))
      end subroutine

      subroutine coretail_density_ff(ao,it,cdff)
!doc$ subroutine coretail_density_ff(ao,it,cdff)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:), pointer :: cdff
!       requires: it in range.
!       modifies: cdff
!       effects: Points cdff to the coretail density form factors of type it.

!cod$
        call my(ao)
        cdff => ao%o%tdata(it)%cdff
        call glean(thy(ao))
      end subroutine

      subroutine hat_density_ff(ao,it,hdff)
!doc$ subroutine hat_density_ff(ao,it,hdff)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:,:), pointer :: hdff
!       requires: it in range.
!       modifies: hdff
!       effects: Points hdff to the hat density form factors of type it.

!cod$
        call my(ao)
        hdff => ao%o%tdata(it)%hdff
        call glean(thy(ao))
      end subroutine

      subroutine density_matrix_ff(ao,it,dmff)
!doc$ subroutine density_matrix_ff(ao,it,dmff)
        type(atomic_operators_paw_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:,:), pointer :: dmff
!       requires: it in range.
!       modifies: dmff
!       effects: Points dmff to the density matrix form factors of type it.

!cod$
        call my(ao)
        dmff => ao%o%tdata(it)%dmff
        call glean(thy(ao))
      end subroutine

      subroutine write_restart_ao(ao,nrestf)
!doc$ subroutine write_restart(ao,nrestf)
        type(atomic_operators_paw_obj) :: ao
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ao restart information to nrestf.

!cod$
        integer :: it, nt
        integer(long) :: dsize, iosl, ndata, s4
        integer(long), dimension(:), allocatable :: v4

        call my(ao)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the ATOMIC_OPERATORS block
          call startblock(nrestf,"ATOMIC_OPERATORS")

          ! write the number of types
          call writetag(nrestf,"NUMBER_OF_TYPES")
          nt = size(ao%o%tdata)
          s4 = nt ; dsize = sizeof_long ; ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          allocate( v4(nt) )

          ! write the theta_points
          call writetag(nrestf,"THETA_POINTS")
          do it = 1,nt
            v4(it) = ao%o%tdata(it)%theta_points
          end do
          dsize = sizeof_long ; ndata = nt
          call writef(v4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! write the phi_points
          call writetag(nrestf,"PHI_POINTS")
          do it = 1,nt
            v4(it) = ao%o%tdata(it)%phi_points
          end do
          dsize = sizeof_long ; ndata = nt
          call writef(v4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! end the ATOMIC_OPERATORS block
          call endblock(nrestf)

        end if

        if (allocated( v4 )) deallocate( v4 )

        call glean(thy(ao))
        call glean(thy(nrestf))

        if (error("Exit atomic_operators_paw_mod::write_restart_ao")) continue

      end subroutine


      function atomic_exchange_density_ao(ao,pdots1, pdots2, lay) result(den)
!doc$ function atomic_exchange_density(ao,lay) result(den)
        type(atomic_operators_paw_obj) :: ao
        complex(double), dimension(:), intent(in) :: pdots1, pdots2
        type(layout_obj) :: lay
        type(grid_obj) :: den
!       effects: Returns the "compensation density" correcting the product of two smooth wavefunctions
!       errors: If size(pdots1) and/or size(pdots2) is not equal to x_n_projectors(ao)

!cod$
        integer :: np, na, ia, it, ml, n1, n2, nm, ab
        complex(double), dimension(:), pointer :: tv1, tv2
        type(gp_mat), dimension(:), pointer :: qlm

        integer :: qb, i, bi, bj, mi, mj, l, m, omi, omj
        integer :: nt, i1, i2, i3
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        real(double), dimension(:,:,:,:), pointer :: hdff
        complex(double) :: c0, c1, igr
        complex(double), dimension(13) :: ylm
        complex(double), dimension(:,:,:), pointer :: c_den
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(ao)
        call my(lay)

        np = x_n_projectors(ao)
        if (error(np /= size(pdots1),"ERROR: size(pdots1) is inconsistent with this ao")) goto 200
        if (error(np /= size(pdots2),"ERROR: size(pdots2) is inconsistent with this ao")) goto 200

        na = x_n_atoms(ao)

        allocate( qlm(na) )
        do ia = 1,na
           it = x_atom_type(ao,ia)
           ml = x_type_l_max(ao,it)
           n1 = 2*ml + 1
           n2 = 2*(2*ml) + 1
           allocate( qlm(ia)%mat(n1,n2) )
           qlm(ia)%mat = (0.0_double,0.0_double)
        end do

        do ia = 1, na
          nm = x_n_atom_projectors(ao,ia)
          allocate(tv1(nm),tv2(nm))

          ab = x_atom_base(ao,ia)
          tv1 = pdots1(ab:ab+nm-1)  !  This is a prime suspect for an errot
          tv2 = conjg(pdots2(ab:ab+nm-1))

          it = x_atom_type(ao,ia)
          tv1 = matmul(ao%o%tdata(it)%tc2r,tv1)
          tv2 = matmul(tv2,ao%o%tdata(it)%tr2c)

          ml = x_type_l_max(ao,it)
          qb = 2*ml + 1

          do i = 1,x_qvlm_size(ao%o%tdata(it)%pawd)
            call orbital_decode(ao%o%tdata(it)%pawd,i,bi,bj,mi,mj,l)
            m = mj - mi
            omi = x_basis_offset(ao,it,bi) + mi
            omj = x_basis_offset(ao,it,bj) + mj
            qlm(ia)%mat(l+1,qb+m) = qlm(ia)%mat(l+1,qb+m) + x_aqlm(ao%o%tdata(it)%pawd,i)*tv1(omi)*tv2(omj)
          end do

          deallocate(tv1,tv2)
        end do

        call my(x_lattice(x_crystal(ao)),lat)
        call my(x_atoms(x_crystal(ao)),ats)
        call my(grid(lay,SGROUP),den)

        nullify( gx, gy, gz, c_den )

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call alloc(c_den,lay,D_TYPE,SGROUP)
        c_den = (0.0_double,0.0_double)
        nt = x_n_types(ao)
        do it = 1, nt
          call hat_density_ff(ao,it,hdff)
          qb = 2*x_type_l_max(ao,it) + 1
          do ia = 1, na
            if (x_atom_type(ao,ia) /= it) cycle
            pos = lat2r(lat,x_position(ats,ia))
            do i3 = 1,size(gx,3)
            do i2 = 1,size(gx,2)
            do i1 = 1,size(gx,1)
              if (hdff(i1,i2,i3,1) == 0.0_double) cycle
              c0 = (0.0_double,0.0_double)
              do l = 0,(size(hdff,4) - 1)
                c1 = hdff(i1,i2,i3,l+1)*(0.0_double,-1.0_double)**l
                ylm = spharm(gx(i1,i2,i3),gy(i1,i2,i3),gz(i1,i2,i3),l,.true.)
                do m = -l,+l
                  c0 = c0 + c1*ylm(l+m+1)*qlm(ia)%mat(l+1,qb+m)
                end do
              end do
              igr = (0.0_double,1.0_double)*(pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3))
              c_den(i1,i2,i3) = c_den(i1,i2,i3) + c0*exp(-igr)
            end do
            end do
            end do
          end do
        end do
        call put(c_den,den,CDF_KIND)

        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )
        nullify( hdff )

        do ia = 1,na
           deallocate(qlm(ia)%mat)
        end do
        deallocate(qlm)

        call glean(thy(lat))
        call glean(thy(ats))
        call bequeath(thy(den))
 
200     call glean(thy(ao))
        call glean(thy(lay))

        if (error("Exit atomic_density_paw_mod::atomic_exchange_density_ao")) continue

      end function


! private routines

      subroutine form_tdata_i(aor,lay,restf)
        type(atomic_operators_paw_rep) :: aor
        type(layout_obj) :: lay
        type(tagio_obj), optional :: restf

        logical :: found, match
        character(1) :: tios
        character(tag_sz) :: tag
        character(11+tag_sz) :: pp_tag
        character(13+tag_sz) :: tp_tag
        character(19+tag_sz) :: bo_tag
        character(tag_sz), dimension(:), allocatable :: tags
        integer :: ia, ib, ios, ip, iph, it, ith, iv, lm_size, maxl, na, nb, np, nt, nv
        integer :: li, lj, mi, mj, oi, oj, osi, osj
        integer :: msg, nsg
        integer :: r_nt
        integer(long) :: dsize, iosl, ndata, s4
        integer, dimension(3) :: nlm
        integer, dimension(:), allocatable :: r_theta_points, r_phi_points
        integer(long), dimension(:), allocatable :: v4
        real(double), parameter :: ll = -1.0_double, ul = +1.0_double
        real(double) :: angle, sin_theta, cos_phi, sin_phi, spin_factor
        real(double), dimension(3) :: va
        real(double), dimension(:), allocatable :: cos_theta, weight
        real(double), dimension(:,:), allocatable :: v
        complex(double), dimension(13) :: ylmi, ylmj, dylmi_dt, dylmj_dt, dylmi_dp, dylmj_dp
        type(file_obj) :: f

        call my(lay)
        if (present(restf)) call my(restf)

        nsg = mpi_nsgroups()
        spin_factor = 1.0_double/real(nsg,double)
        msg = mpi_mysgroup()

        na = x_n_atoms(x_atoms(aor%cr))

        ! Determine the number of atom types and their tags.
        allocate( tags(na) )
        nt = 0
        do ia = 1,na
          tag = x_type(x_atoms(aor%cr),ia)
          match = .false.
          do it = 1,nt
            match = ( match .or. (tag == tags(it)) )
          end do
          if (match) cycle
          nt = nt + 1
          tags(nt) = tag
        end do

        ! Allocate aor%tdata and set the tag data.
        allocate( aor%tdata(nt) )
        do it = 1,nt
          aor%tdata(it)%tag = tags(it)
        end do
        deallocate( tags )

        ! Construct the paw data objects.
        do it = 1,nt
          call my(file(trim(paw_path)//aor%tdata(it)%tag),f)
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: problem opening file = "//x_name(f))) goto 100
          call my(paw_data(f),aor%tdata(it)%pawd)
          if (i_access(f)) close(x_unit(f))
100       call glean(thy(f)) ; if (error()) goto 200
        end do

        ! Set the valence data.
        do it = 1,nt
          aor%tdata(it)%valence = x_valence_electrons(aor%tdata(it)%pawd)
        end do

        ! Set the matching radius data.
        do it = 1,nt
          aor%tdata(it)%matching_radius = x_matching_radius(aor%tdata(it)%pawd)
        end do

        ! Set the projector data.
        do it = 1,nt
          np = x_nlm_size(aor%tdata(it)%pawd)
          if (np == 0) then
            nullify( aor%tdata(it)%pn )
            nullify( aor%tdata(it)%pl )
            nullify( aor%tdata(it)%pm )
          else
            allocate( aor%tdata(it)%pn(np) )
            allocate( aor%tdata(it)%pl(np) )
            allocate( aor%tdata(it)%pm(np) )
            do ip = 1,np
              nlm = x_nlm(aor%tdata(it)%pawd,ip)
              aor%tdata(it)%pn(ip) = nlm(1)
              aor%tdata(it)%pl(ip) = nlm(2)
              aor%tdata(it)%pm(ip) = nlm(3)
            end do
          end if
        end do

        ! Set the basis data.
        do it = 1,nt
          nb = x_basis_size(aor%tdata(it)%pawd)
          allocate( aor%tdata(it)%bs(nb) )
          allocate( aor%tdata(it)%bl(nb) )
          allocate( aor%tdata(it)%bo(nb) )
          do ib = 1,nb
            aor%tdata(it)%bs(ib) = x_nl_base(aor%tdata(it)%pawd,ib)
            aor%tdata(it)%bl(ib) = x_l_value(aor%tdata(it)%pawd,ib)
          end do
          select case (nsg)
          case (1)
            bo_tag = "atomic_occupations_"//aor%tdata(it)%tag
            call arg(trim(bo_tag),aor%tdata(it)%bo,found)
          case (2)
            select case (msg)
            case (1)
              bo_tag = "atomic_occupations_sg1_"//aor%tdata(it)%tag
              call arg(trim(bo_tag),aor%tdata(it)%bo,found)
            case (2)
              bo_tag = "atomic_occupations_sg2_"//aor%tdata(it)%tag
              call arg(trim(bo_tag),aor%tdata(it)%bo,found)
            end select
          end select
          if (.not.found) then
            do ib = 1,nb
              aor%tdata(it)%bo(ib) = spin_factor*x_occupation(aor%tdata(it)%pawd,ib)
            end do
          end if
        end do

        ! Construct the real-to-complex and complex-to-real Ylm transformation matrices.
        call form_r2c_and_c2r_i(aor) ; if (error()) goto 200

        ! Construct the grid-potential, core-density, hat-density, and density-matrix form factors.
        do it = 1,nt
          nullify( aor%tdata(it)%gpff )
          nullify( aor%tdata(it)%cdff )
          nullify( aor%tdata(it)%hdff )
          nullify( aor%tdata(it)%dmff )
        end do
        call form_factors_i(aor,lay) ; if (error()) goto 200

        ! read the number of types from the restart file
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_TYPES")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_TYPES tag was not found")) goto 200
          if (i_access(restf)) then
            dsize = sizeof_long ; ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            r_nt = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_nt)
          if (error(r_nt /= nt,"ERROR: different numbers of types")) goto 200
        end if

        ! read the theta and phi points from the restart file
        if (present(restf)) then

          allocate( v4(nt), r_theta_points(nt), r_phi_points(nt) )
          dsize = sizeof_long ; ndata = nt

          if (i_access(restf)) tios = findfirsttag(restf,"THETA_POINTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: THETA_POINTS tag was not found")) goto 200
          if (i_access(restf)) then
            call readf(v4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            r_theta_points = v4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_theta_points)

          if (i_access(restf)) tios = findfirsttag(restf,"PHI_POINTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: PHI_POINTS tag was not found")) goto 200
          if (i_access(restf)) then
            call readf(v4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            r_phi_points = v4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_phi_points)

        end if

        ! Construct quantities needed for angular integration.
        do it = 1,nt

          tp_tag = "theta_points_"//aor%tdata(it)%tag
          call arg(trim(tp_tag),aor%tdata(it)%theta_points,found)
          if (found) then
            if (error(aor%tdata(it)%theta_points <= 0,"ERROR: theta_points <= 0")) goto 200
          else
            if (present(restf)) then
              aor%tdata(it)%theta_points = r_theta_points(it)
            else
              aor%tdata(it)%theta_points = 12
            end if
          end if

          pp_tag = "phi_points_"//aor%tdata(it)%tag
          call arg(trim(pp_tag),aor%tdata(it)%phi_points,found)
          if (found) then
            if (error(aor%tdata(it)%phi_points <= 0,"ERROR: phi_points <= 0")) goto 200
          else
            if (present(restf)) then
              aor%tdata(it)%phi_points = r_phi_points(it)
            else
              aor%tdata(it)%phi_points = 12
            end if
          end if

          if (allocated( cos_theta )) deallocate( cos_theta ) ; allocate( cos_theta(aor%tdata(it)%theta_points) )
          if (allocated( weight )) deallocate( weight ) ; allocate( weight(aor%tdata(it)%theta_points) )
          call gauss_legendre(ll,ul,cos_theta,weight)
          if (allocated( v )) deallocate( v ) ; allocate( v(3,aor%tdata(it)%theta_points*aor%tdata(it)%phi_points) )
          allocate( aor%tdata(it)%v_weight(aor%tdata(it)%theta_points*aor%tdata(it)%phi_points) )
          nv = 0 
          do ith = 1,aor%tdata(it)%theta_points
            sin_theta = sqrt(1.0_double - cos_theta(ith)*cos_theta(ith))
            do iph = 1,aor%tdata(it)%phi_points
              angle = two_pi*real(iph-1,double)/real(aor%tdata(it)%phi_points,double)
              cos_phi = cos(angle)
              sin_phi = sin(angle)
              nv = nv + 1
              v(:,nv) = (/sin_theta*cos_phi,sin_theta*sin_phi,cos_theta(ith)/)
              aor%tdata(it)%v_weight(nv) = four_pi*weight(ith)/real(2*aor%tdata(it)%phi_points,double)
            end do
          end do

          maxl = maxval(aor%tdata(it)%pl)
          lm_size = maxl**2 + 2*maxl + 1

          if (error(uses_laplacian(aor%axc),"ERROR: Functionals with laplacians are not supported with PAW")) goto 200

          if (uses_gradient(aor%axc)) then
            allocate( aor%tdata(it)%ylmij(lm_size,lm_size,nv) )
            allocate( aor%tdata(it)%dylmij_dt(lm_size,lm_size,nv) )
            allocate( aor%tdata(it)%dylmij_dp(lm_size,lm_size,nv) )
            do li = 0,maxl
              oi = li**2 + li + 1
              osi = li + 1
              do lj = 0,maxl
                oj = lj**2 + lj + 1
                osj = lj + 1
                do iv = 1,nv
                  va = v(:,iv)
                  ylmi = spharm(va(1),va(2),va(3),li,.true.)
                  ylmj = spharm(va(1),va(2),va(3),lj,.true.)
                  call gradylm(va(1),va(2),va(3),li,dylmi_dt,dylmi_dp)
                  call gradylm(va(1),va(2),va(3),lj,dylmj_dt,dylmj_dp)
                  do mi = -li,+li
                  do mj = -lj,+lj
                    aor%tdata(it)%ylmij(oi+mi,oj+mj,iv) = conjg(ylmi(osi+mi))*ylmj(osj+mj)
                    aor%tdata(it)%dylmij_dt(oi+mi,oj+mj,iv) = conjg(dylmi_dt(osi+mi))*ylmj(osj+mj) +  &
                                                              conjg(ylmi(osi+mi))*dylmj_dt(osj+mj)
                    aor%tdata(it)%dylmij_dp(oi+mi,oj+mj,iv) = conjg(dylmi_dp(osi+mi))*ylmj(osj+mj) +  &
                                                              conjg(ylmi(osi+mi))*dylmj_dp(osj+mj)
                  end do
                  end do
                end do
              end do
            end do
          else
            allocate( aor%tdata(it)%ylmij(lm_size,lm_size,nv) )
            nullify( aor%tdata(it)%dylmij_dt, aor%tdata(it)%dylmij_dp )
            do li = 0,maxl
              oi = li**2 + li + 1
              osi = li + 1
              do lj = 0,maxl
                oj = lj**2 + lj + 1
                osj = lj + 1
                do iv = 1,nv
                  va = v(:,iv)
                  ylmi = spharm(va(1),va(2),va(3),li,.true.)
                  ylmj = spharm(va(1),va(2),va(3),lj,.true.)
                  do mi = -li,+li
                  do mj = -lj,+lj
                    aor%tdata(it)%ylmij(oi+mi,oj+mj,iv) = conjg(ylmi(osi+mi))*ylmj(osj+mj)
                  end do
                  end do
                end do
              end do
            end do
          end if

        end do

200     if (allocated( tags )) deallocate( tags )
        if (allocated( r_theta_points )) deallocate( r_theta_points )
        if (allocated( r_phi_points )) deallocate( r_phi_points )
        if (allocated( v4 )) deallocate( v4 )
        if (allocated( cos_theta )) deallocate( cos_theta )
        if (allocated( weight )) deallocate( weight )
        if (allocated( v )) deallocate( v )

        call glean(thy(lay))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_operators_paw_mod::form_tdata_i")) continue

      end subroutine

      subroutine form_adata_i(aor)
        type(atomic_operators_paw_rep) :: aor

        character(tag_sz) :: tag
        integer :: ia, it, na, np
        real(double), parameter :: tol_val_sum = 1.0e-8_double
        real(double) :: occ_sum, occ_sum_sg, val_sum

        ! Allocate aor%adata.
        na = x_n_atoms(x_atoms(aor%cr))
        allocate( aor%adata(na) )

        ! Set the data.
        np = 0
        do ia = 1,na
          tag = x_type(x_atoms(aor%cr),ia)
          do it = 1,size(aor%tdata)
            if (tag == aor%tdata(it)%tag) then
              aor%adata(ia)%ti = it
              exit
            end if
          end do
          aor%adata(ia)%bi = np + 1
          aor%adata(ia)%valence = aor%tdata(aor%adata(ia)%ti)%valence
          if (associated( aor%tdata(aor%adata(ia)%ti)%pn )) np = np + size(aor%tdata(aor%adata(ia)%ti)%pn)
        end do

        ! Check the initial occupations.
        val_sum = 0.0_double
        occ_sum_sg = 0.0_double
        do ia = 1,na
          val_sum = val_sum + aor%adata(ia)%valence
          it = aor%adata(ia)%ti
          occ_sum_sg = occ_sum_sg + sum(aor%tdata(it)%bo)
        end do
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,occ_sum_sg,occ_sum) ; if (error()) goto 100
        if (error(occ_sum .out. nbhd(val_sum,tol_val_sum),"ERROR: initial occupations sum is incorrect")) goto 100

        call form_oij_i(aor)

100     if (error("Exit atomic_operators_paw_mod::form_adata_i")) continue

      end subroutine

      subroutine form_pdata_i(aor)
        type(atomic_operators_paw_rep) :: aor

        integer :: ia, ip, it, itp, np

        ! Determine the total number of projectors.
        np = 0
        do ia = 1,size(aor%adata)
          np = np + size(aor%tdata(aor%adata(ia)%ti)%pn)
        end do

        ! Allocate aor%pdata and set the data.
        if (np == 0) then
          nullify( aor%pdata )
        else
          allocate( aor%pdata(np) )
          ip = 0
          do ia = 1,size(aor%adata)
            it = aor%adata(ia)%ti
            if (associated( aor%tdata(it)%pn )) then
              do itp = 1,size(aor%tdata(it)%pn)
                ip = ip + 1
                aor%pdata(ip)%ai = ia
                aor%pdata(ip)%ti = it
                aor%pdata(ip)%tpi = itp
                aor%pdata(ip)%n = aor%tdata(it)%pn(itp)
                aor%pdata(ip)%l = aor%tdata(it)%pl(itp)
                aor%pdata(ip)%m = aor%tdata(it)%pm(itp)
              end do
            end if
          end do
        end if

        if (error("Exit atomic_operators_paw_mod::form_pdata_i")) continue

      end subroutine

      subroutine form_udata_i(aor)
        type(atomic_operators_paw_rep) :: aor

        logical :: found
        character(line_len) :: tag
        integer :: ia, iu, na, nu
        integer, dimension(:), pointer :: sua

        nullify( sua )

        na = size(aor%adata)

        call arglc("atomic_symmetry",tag,found)
        if (.not.found) tag = "on"
        select case (trim(tag))
        case ("off")
          allocate( aor%udata(na) )
          do ia = 1,na
            aor%udata(ia)%ai = ia
            aor%udata(ia)%ti = aor%adata(ia)%ti
          end do
        case ("on")
          allocate( sua(na) )
          call symmetry_unique_atoms(aor%sg,sua)
          nu = size(sua)
          allocate( aor%udata(nu) )
          do iu = 1,nu
            ia = sua(iu)
            aor%udata(iu)%ai = ia
            aor%udata(iu)%ti = aor%adata(ia)%ti
          end do
        case default
          if (error(.true.,"ERROR: atomic_symmetry tag is not recognized")) goto 100
        end select

100     if (associated( sua )) deallocate( sua )

        if (error("Exit atomic_operators_paw_mod::form_udata_i")) continue

      end subroutine

      subroutine form_vdata_i(aor)
        type(atomic_operators_paw_rep) :: aor

        integer :: i, ia, it, iu, nu
        integer :: first_v, last_v, iv, ivt, nv, nv_proc

        nu = size(aor%udata)

        ! Determine the total number of vectors.
        nv = 0
        do iu = 1,nu
          ia = aor%udata(iu)%ai
          it = aor%adata(ia)%ti
          nv = nv + size(aor%tdata(it)%ylmij,3)
        end do

        ! Divide up the vectors and set the data.
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,nv,first_v,last_v,nv_proc)
        allocate( aor%vdata(nv_proc) )
        iv = 0
        do iu = 1,nu
          ia = aor%udata(iu)%ai
          it = aor%adata(ia)%ti
          do ivt = 1,size(aor%tdata(it)%ylmij,3)
            iv = iv + 1
            if ((iv < first_v) .or. (iv > last_v)) cycle
            i = iv - first_v + 1
            aor%vdata(i)%ai = ia
            aor%vdata(i)%ui = iu
            aor%vdata(i)%ti = it
            aor%vdata(i)%tvi = ivt
          end do
        end do

        if (error("Exit atomic_operators_paw_mod::form_vdata_i")) continue

      end subroutine

      subroutine own_i(ao)
        type(atomic_operators_paw_obj) :: ao

        type(atomic_operators_paw_obj) :: aot
        integer :: i, n1, n2, n3

        if (ao%ref < ao%o%ref) then
          allocate ( aot%o )
          aot%o%ref = 0
          aot%o%g = ao%o%g
          aot%o%g_atoms = ao%o%g_atoms
          aot%o%g_layout = ao%o%g_layout
          call my(ao%o%cr,aot%o%cr)
          call my(ao%o%sg,aot%o%sg)
          call my(ao%o%axc,aot%o%axc)
          allocate( aot%o%tdata(size(ao%o%tdata)) )
          do i = 1,size(aot%o%tdata)
            call my(ao%o%tdata(i)%pawd,aot%o%tdata(i)%pawd)
            aot%o%tdata(i)%tag = ao%o%tdata(i)%tag
            aot%o%tdata(i)%valence = ao%o%tdata(i)%valence
            aot%o%tdata(i)%matching_radius = ao%o%tdata(i)%matching_radius
            if (associated( ao%o%tdata(i)%pn )) then
              allocate( aot%o%tdata(i)%pn(size(ao%o%tdata(i)%pn)) ) ; aot%o%tdata(i)%pn = ao%o%tdata(i)%pn
              allocate( aot%o%tdata(i)%pl(size(ao%o%tdata(i)%pl)) ) ; aot%o%tdata(i)%pl = ao%o%tdata(i)%pl
              allocate( aot%o%tdata(i)%pm(size(ao%o%tdata(i)%pm)) ) ; aot%o%tdata(i)%pm = ao%o%tdata(i)%pm
            else
              nullify( aot%o%tdata(i)%pn )
              nullify( aot%o%tdata(i)%pl )
              nullify( aot%o%tdata(i)%pm )
            end if
            allocate( aot%o%tdata(i)%bs(size(ao%o%tdata(i)%bs)) ) ; aot%o%tdata(i)%bs = ao%o%tdata(i)%bs
            allocate( aot%o%tdata(i)%bl(size(ao%o%tdata(i)%bl)) ) ; aot%o%tdata(i)%bl = ao%o%tdata(i)%bl
            allocate( aot%o%tdata(i)%bo(size(ao%o%tdata(i)%bo)) ) ; aot%o%tdata(i)%bo = ao%o%tdata(i)%bo
            allocate( aot%o%tdata(i)%gpff(size(ao%o%tdata(i)%gpff,1),size(ao%o%tdata(i)%gpff,2),size(ao%o%tdata(i)%gpff,3)) )
            aot%o%tdata(i)%gpff = ao%o%tdata(i)%gpff
            allocate( aot%o%tdata(i)%cdff(size(ao%o%tdata(i)%cdff,1),size(ao%o%tdata(i)%cdff,2),size(ao%o%tdata(i)%cdff,3)) )
            aot%o%tdata(i)%cdff = ao%o%tdata(i)%cdff
            allocate( aot%o%tdata(i)%hdff(size(ao%o%tdata(i)%hdff,1),size(ao%o%tdata(i)%hdff,2),size(ao%o%tdata(i)%hdff,3), &
                                          size(ao%o%tdata(i)%hdff,4)) )
            aot%o%tdata(i)%hdff = ao%o%tdata(i)%hdff
            allocate( aot%o%tdata(i)%dmff(size(ao%o%tdata(i)%dmff,1),size(ao%o%tdata(i)%dmff,2),size(ao%o%tdata(i)%dmff,3), &
                                          size(ao%o%tdata(i)%dmff,4)) )
            aot%o%tdata(i)%dmff = ao%o%tdata(i)%dmff
            allocate( aot%o%tdata(i)%r2c(size(ao%o%tdata(i)%r2c,1),size(ao%o%tdata(i)%r2c,2)) )
            aot%o%tdata(i)%r2c = ao%o%tdata(i)%r2c
            allocate( aot%o%tdata(i)%c2r(size(ao%o%tdata(i)%c2r,1),size(ao%o%tdata(i)%c2r,2)) )
            aot%o%tdata(i)%c2r = ao%o%tdata(i)%c2r
            allocate( aot%o%tdata(i)%tr2c(size(ao%o%tdata(i)%tr2c,1),size(ao%o%tdata(i)%tr2c,2)) )
            aot%o%tdata(i)%tr2c = ao%o%tdata(i)%tr2c
            allocate( aot%o%tdata(i)%tc2r(size(ao%o%tdata(i)%tc2r,1),size(ao%o%tdata(i)%tc2r,2)) )
            aot%o%tdata(i)%tc2r = ao%o%tdata(i)%tc2r
            aot%o%tdata(i)%theta_points = ao%o%tdata(i)%theta_points
            aot%o%tdata(i)%phi_points = ao%o%tdata(i)%phi_points
            allocate( aot%o%tdata(i)%v_weight(size(ao%o%tdata(i)%v_weight)) )
            aot%o%tdata(i)%v_weight = ao%o%tdata(i)%v_weight
            n1 = size(ao%o%tdata(i)%ylmij,1) ; n2 = size(ao%o%tdata(i)%ylmij,2) ; n3 = size(ao%o%tdata(i)%ylmij,3)
            allocate( aot%o%tdata(i)%ylmij(n1,n2,n3) )
            aot%o%tdata(i)%ylmij = ao%o%tdata(i)%ylmij
            if (associated(ao%o%tdata(i)%dylmij_dt)) then
              allocate( aot%o%tdata(i)%dylmij_dt(n1,n2,n3) )
              aot%o%tdata(i)%dylmij_dt = ao%o%tdata(i)%dylmij_dt
            else
              nullify( aot%o%tdata(i)%dylmij_dt )
            end if
            if (associated(ao%o%tdata(i)%dylmij_dp)) then
              allocate( aot%o%tdata(i)%dylmij_dp(n1,n2,n3) )
              aot%o%tdata(i)%dylmij_dp = ao%o%tdata(i)%dylmij_dp
            else
              nullify( aot%o%tdata(i)%dylmij_dp )
            end if
          end do
          allocate( aot%o%adata(size(ao%o%adata)) )
          do i = 1,size(aot%o%adata)
            aot%o%adata(i)%ti = ao%o%adata(i)%ti
            aot%o%adata(i)%bi = ao%o%adata(i)%bi
            aot%o%adata(i)%valence = ao%o%adata(i)%valence
            allocate( aot%o%adata(i)%c_oij(size(ao%o%adata(i)%c_oij,1),size(ao%o%adata(i)%c_oij,2)) )
            aot%o%adata(i)%c_oij = ao%o%adata(i)%c_oij
            allocate( aot%o%adata(i)%r_oij(size(ao%o%adata(i)%r_oij,1),size(ao%o%adata(i)%r_oij,2)) )
            aot%o%adata(i)%r_oij = ao%o%adata(i)%r_oij
          end do
          if (associated( ao%o%pdata )) then
            allocate( aot%o%pdata(size(ao%o%pdata)) )
            do i = 1,size(aot%o%pdata)
              aot%o%pdata(i)%ai = ao%o%pdata(i)%ai
              aot%o%pdata(i)%ti = ao%o%pdata(i)%ti
              aot%o%pdata(i)%tpi = ao%o%pdata(i)%tpi
              aot%o%pdata(i)%n = ao%o%pdata(i)%n
              aot%o%pdata(i)%l = ao%o%pdata(i)%l
              aot%o%pdata(i)%m = ao%o%pdata(i)%m
            end do
          else
            nullify( aot%o%pdata )
          end if
          if (associated( ao%o%udata )) then
            allocate( aot%o%udata(size(ao%o%udata)) )
            do i = 1,size(aot%o%udata)
              aot%o%udata(i)%ai = ao%o%udata(i)%ai
              aot%o%udata(i)%ti = ao%o%udata(i)%ti
            end do
          else
            nullify( aot%o%udata )
          end if
          if (associated( ao%o%vdata )) then
            allocate( aot%o%vdata(size(ao%o%vdata)) )
            do i = 1,size(aot%o%vdata)
              aot%o%vdata(i)%ai = ao%o%vdata(i)%ai
              aot%o%vdata(i)%ui = ao%o%vdata(i)%ui
              aot%o%vdata(i)%ti = ao%o%vdata(i)%ti
              aot%o%vdata(i)%tvi = ao%o%vdata(i)%tvi
            end do
          else
            nullify( aot%o%vdata )
          end if
          call my(ao%o%ctd,aot%o%ctd)
          ao%o%ref = ao%o%ref - ao%ref
          ao%o => aot%o
          ao%o%ref = ao%o%ref + ao%ref
        end if
      end subroutine
 
      subroutine form_r2c_and_c2r_i(aor)
        type(atomic_operators_paw_rep) :: aor

        integer :: ic, ir, it, mc, mr, nc, np, nr
        real(double), parameter :: minus_one = -1.0_double
        real(double) :: norm
        complex(double), parameter :: i = (0.0_double,1.0_double)

        norm = 1.0_double/sqrt(2.0_double)

        do it = 1,size(aor%tdata)
          np = size(aor%tdata(it)%pn)
          allocate( aor%tdata(it)%r2c(np,np), aor%tdata(it)%c2r(np,np), aor%tdata(it)%tr2c(np,np), aor%tdata(it)%tc2r(np,np) )
          do ic = 1,np
            do ir = 1,np
              nr = aor%tdata(it)%pn(ir)
              nc = aor%tdata(it)%pn(ic)
              mr = aor%tdata(it)%pm(ir)
              mc = aor%tdata(it)%pm(ic)
              if ((nr == nc) .and. (abs(mr) == abs(mc))) then
                if ((mr < 0) .and. (mc < 0)) then
                  aor%tdata(it)%r2c(ir,ic) = i*norm*minus_one**abs(mr)
                elseif ((mr < 0) .and. (mc > 0)) then
                  aor%tdata(it)%r2c(ir,ic) = norm*minus_one**abs(mr)
                elseif ((mr == 0) .and. (mc == 0)) then
                  aor%tdata(it)%r2c(ir,ic) = 1.0_double
                elseif ((mr > 0) .and. (mc < 0)) then
                  aor%tdata(it)%r2c(ir,ic) = -i*norm
                elseif ((mr > 0) .and. (mc > 0)) then
                  aor%tdata(it)%r2c(ir,ic) = norm
                end if
              else
                aor%tdata(it)%r2c(ir,ic) = (0.0_double,0.0_double)
              end if
            end do
          end do
          aor%tdata(it)%tr2c = transpose(aor%tdata(it)%r2c)
          aor%tdata(it)%c2r = conjg(aor%tdata(it)%tr2c)
          aor%tdata(it)%tc2r = conjg(aor%tdata(it)%r2c)
        end do

        if (error("Exit atomic_operators_paw_mod::form_r2c_and_c2r_i")) continue

      end subroutine

      subroutine form_factors_i(aor,lay)
        type(atomic_operators_paw_rep) :: aor
        type(layout_obj) :: lay
!       requires: Relevant pointers be nullified or associated.

        integer :: i1, i2, i3, it, n1, n2, n3
        real(double) :: cutoff, gm, volume
        real(double), dimension(:,:,:), pointer :: g2

        call my(lay)

        nullify( g2 )

        cutoff = x_cutoff(lay)
        volume = x_cell_volume(x_lattice(aor%cr))

        call fdel(g2,lay,D_TYPE,SGROUP)

        n1 = size(g2,1)
        n2 = size(g2,2)
        n3 = size(g2,3)

        do it = 1,size(aor%tdata)

          if (associated( aor%tdata(it)%gpff )) deallocate( aor%tdata(it)%gpff )
          if (associated( aor%tdata(it)%cdff )) deallocate( aor%tdata(it)%cdff )
          if (associated( aor%tdata(it)%hdff )) deallocate( aor%tdata(it)%hdff )
          if (associated( aor%tdata(it)%dmff )) deallocate( aor%tdata(it)%dmff )

          allocate( aor%tdata(it)%gpff(n1,n2,n3) )
          allocate( aor%tdata(it)%cdff(n1,n2,n3) )
          allocate( aor%tdata(it)%hdff(n1,n2,n3,x_hatden_size(aor%tdata(it)%pawd)) )
          allocate( aor%tdata(it)%dmff(n1,n2,n3,x_denmat_size(aor%tdata(it)%pawd)) )

          ! construct the form factors
          do i3 = 1,n3
          do i2 = 1,n2
          do i1 = 1,n1
            if (g2(i1,i2,i3) > cutoff) then
              aor%tdata(it)%gpff(i1,i2,i3) = 0.0_double
              aor%tdata(it)%cdff(i1,i2,i3) = 0.0_double
              aor%tdata(it)%hdff(i1,i2,i3,:) = 0.0_double
              aor%tdata(it)%dmff(i1,i2,i3,:) = 0.0_double
            else
              gm = sqrt(g2(i1,i2,i3))
              aor%tdata(it)%gpff(i1,i2,i3) = local_f_value(aor%tdata(it)%pawd,gm)/volume
              aor%tdata(it)%cdff(i1,i2,i3) = coretail_f_value(aor%tdata(it)%pawd,gm)/volume
              aor%tdata(it)%hdff(i1,i2,i3,:) = hatden_f_value(aor%tdata(it)%pawd,gm)/volume
              aor%tdata(it)%dmff(i1,i2,i3,:) = denmat_f_value(aor%tdata(it)%pawd,gm)/volume
            end if
          end do
          end do
          end do

        end do

100     if (associated( g2 )) deallocate( g2 )

        call glean(thy(lay))

        if (error("Exit atomic_operators_paw_mod::form_factors_i")) continue

      end subroutine

      subroutine form_oij_i(aor)
        type(atomic_operators_paw_rep) :: aor

        integer :: ia, it, np

        do ia = 1,size(aor%adata)
          it = aor%adata(ia)%ti
          np = size(aor%tdata(it)%pn)
          allocate( aor%adata(ia)%c_oij(np,np), aor%adata(ia)%r_oij(np,np) )
          aor%adata(ia)%c_oij = cmplx(x_oij(aor%tdata(it)%pawd),0,double)
          aor%adata(ia)%r_oij = matmul(aor%tdata(it)%c2r,matmul(aor%adata(ia)%c_oij,aor%tdata(it)%r2c))
        end do

        if (error("Exit atomic_operators_paw_mod::form_oij_i")) continue

      end subroutine

      subroutine form_coretail_density_i(aor)
        type(atomic_operators_paw_rep) :: aor

        integer :: i1, i2, i3, ia, it
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        complex(double) :: igr
        complex(double), dimension(:,:,:), pointer :: c_ctd
        type(layout_obj) :: lay
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        nullify( gx, gy, gz, c_ctd )

        call my(x_layout(aor%ctd),lay)
        call my(x_lattice(aor%cr),lat)
        call my(x_atoms(aor%cr),ats)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call alloc(c_ctd,lay,D_TYPE,SGROUP)
        c_ctd = (0.0_double,0.0_double)
        do ia = 1,size(aor%adata)
          it = aor%adata(ia)%ti
          pos = lat2r(lat,x_position(ats,ia))
          do i3 = 1,size(gx,3)
          do i2 = 1,size(gx,2)
          do i1 = 1,size(gx,1)
            if (aor%tdata(it)%cdff(i1,i2,i3) == 0.0_double) cycle
            igr = (0.0_double,1.0_double)*(pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3))
            c_ctd(i1,i2,i3) = c_ctd(i1,i2,i3) + aor%tdata(it)%cdff(i1,i2,i3)*exp(-igr)
          end do
          end do
          end do
        end do
        call put(c_ctd,aor%ctd,CDF_KIND)

        call transform(aor%ctd,RD_KIND)

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_ctd )) deallocate( c_ctd )

        call glean(thy(lay))
        call glean(thy(lat))
        call glean(thy(ats))

        if (error("Exit atomic_operators_paw_mod::form_coretail_density_i")) continue

      end subroutine

      end module

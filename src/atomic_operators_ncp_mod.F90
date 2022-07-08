!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module atomic_operators_ncp_mod
!doc$ module atomic_operators_ncp_mod

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use ghost_mod
      use arg_mod
      use layout_mod
      use grid_mod
      use lattice_mod
      use atoms_mod
      use crystal_mod
      use ewald_mod
      use math_mod
      use symmetry_mod
      use ncp_data_mod
      use diary_mod

!     One datatype is available here: type(atomic_operators_ncp_obj).

!cod$
      implicit none
      private

      type :: type_data
        type(ncp_data_obj) :: ncpd                          ! ncp data object
        character(tag_sz) :: tag                            ! tag
        real(double) :: valence                             ! number of valence electrons
        integer, dimension(:), pointer :: pl                ! projector l values
        integer, dimension(:), pointer :: pm                ! projector m values
        integer, dimension(:), pointer :: pc                ! projector channnel,c, index
        complex(double), dimension(:,:), pointer :: r2c           ! real_to_complex Ylm transformation
        complex(double), dimension(:,:), pointer :: c2r           ! complex_to_real Ylm transformation
        complex(double), dimension(:,:), pointer :: tr2c          ! transpose of real_to_complex Ylm transformation
        complex(double), dimension(:,:), pointer :: tc2r          ! transpose of complex_to_real Ylm transformation
        real(double), dimension(:), pointer :: pkbf         ! projector Kleinman-Bylander factors
        real(double), dimension(:,:,:), pointer :: vdff     ! valence density form factors
        real(double), dimension(:,:,:), pointer :: gpff     ! grid potential form factors
        real(double), dimension(:,:,:), pointer :: cdff     ! core density form factors (optional)
        real(double), dimension(:,:,:), pointer :: sgpff    ! stress grid potential form factors
        real(double), dimension(:,:,:), pointer :: scdff    ! stress core density form factors (optional)
      end type

      type :: atom_data
        integer :: ti          ! type index
        integer :: valence     ! valence
        integer :: bi          ! base projector index
      end type

      type :: projector_data
        integer :: ai           ! atom index
        integer :: ti           ! type index
        integer :: tpi          ! type projector index
        integer :: l            ! l value
        integer :: m            ! m value
        integer :: c            ! channel index
        real(double) :: kbf     ! Kleinman-Bylander factor
      end type

      type :: atomic_operators_ncp_rep
        integer :: ref
        type(ghost) :: g
        type(ghost) :: g_atoms                                 ! atoms ghost
        type(ghost) :: g_layout                                ! layout ghost
        real(double) :: energy                                 ! Ewald energy
        type(crystal_obj) :: cr                                ! crystal
        type(type_data), dimension(:), pointer :: tdata        ! type data
        type(atom_data), dimension(:), pointer :: adata        ! atom data
        type(projector_data), dimension(:), pointer :: pdata   ! projector data
      end type

      type, public :: atomic_operators_ncp_obj
        private
        integer :: ref
        type(atomic_operators_ncp_rep), pointer :: o
      end type

!doc$
      public :: atomic_operators_ncp
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
      public :: x_n_types
      public :: x_n_atoms
      public :: x_type_name
      public :: x_type_valence
      public :: x_valence_electrons
      public :: x_n_projectors
      public :: x_n_type_projectors
      public :: x_n_atom_projectors
      public :: x_projector_l
      public :: x_projector_m
      public :: x_projector_kbf
      public :: x_projector_type
      public :: x_projector_atom
      public :: x_projector_index_in_type
      public :: x_atom_type
      public :: x_atom_base
      public :: x_atom_valence
      public :: apply_projector_kbf
      public :: x_atom_energy
      public :: projector_f_values
      public :: projector_stress_f_values
      public :: projector_r_value
      public :: projector_r_gradients
      public :: x_projector_radius
      public :: x_type_projector_radius
      public :: x_ewald_energy
      public :: atomic_grid_potential
      public :: atomic_xc_density
      public :: diary_rs_projectors
      public :: core_density
      public :: valence_density_ff
      public :: grid_potential_ff
      public :: core_density_ff
      public :: stress_grid_potential_ff
      public :: stress_core_density_ff
      public :: write_restart
!      public :: type_r2c
!      public :: type_c2r
      public :: type_tr2c
      public :: type_tc2r


!cod$
      interface atomic_operators_ncp
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
      interface x_projector_kbf
        module procedure ao_projector_kbf, ao_type_projector_kbf
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
      interface apply_projector_kbf
         module procedure ao_apply_projector_kbf
      end interface
      interface x_atom_energy
         module procedure ao_atom_energy
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
      interface x_ewald_energy
        module procedure ao_ewald_energy
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
      interface core_density
        module procedure core_density_ao
      end interface
      interface write_restart
        module procedure write_restart_ao
      end interface

      contains

! public routines

      function constructor_ao(cr,lay,restf) result(ao)
!doc$ function atomic_operators_ncp(cr,lay,restf) result(ao)
        type(crystal_obj) :: cr
        type(layout_obj) :: lay
        type(tagio_obj), optional :: restf
        type(atomic_operators_ncp_obj) :: ao
!       effects: Initializes ao using cr, lay, and (optionally) restf.
!       errors: Passes errors.

!cod$
        logical :: match
        character(1) :: tios
        character(tag_sz) :: tag
        character(tag_sz), dimension(:), allocatable :: tags
        integer :: ia, il, ios, ip, it, itp, l, m, na, np, nt, c, nl
        integer :: r_nt
        integer(long) :: dsize, iosl, ndata, s4
        real(double) :: kbf
        real(double), dimension(:), allocatable :: q
        type(file_obj) :: f

        call my(cr)
        call my(lay)
        if (present(restf)) call my(restf)

        ao%ref = 0
        allocate( ao%o )
        ao%o%ref = 0
        ao%o%g = x_ghost()
        ao%o%g_atoms = x_ghost(x_atoms(cr))
        ao%o%g_layout = x_ghost(lay)

        call my(cr,ao%o%cr)

        ! open the ATOMIC_OPERATORS block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_OPERATORS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_OPERATORS block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)
        end if

        ! get the number of atoms
        na = x_n_atoms(x_atoms(ao%o%cr))

        ! extract the atom types
        allocate( tags(na) )
        nt = 0
        do ia = 1,na
          tag = x_type(x_atoms(ao%o%cr),ia)
          match = .false.
          do it = 1,nt
            match = ( match .or. (tag == tags(it)) )
          end do
          if (match) cycle
          nt = nt + 1
          tags(nt) = tag
        end do
        allocate( ao%o%tdata(nt) )
        do it = 1,size(ao%o%tdata)
          ao%o%tdata(it)%tag = tags(it)
        end do
        deallocate( tags )

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

        ! read the ncp files
        do it = 1,size(ao%o%tdata)
          call my(ncp_data(ao%o%tdata(it)%tag),ao%o%tdata(it)%ncpd)
          ao%o%tdata(it)%valence = x_valence_electrons(ao%o%tdata(it)%ncpd)
          if (i_access( diaryfile() )) then
!          write(x_unit(diaryfile()), '(/,"New Valence of atom is =", f5.2)') ao%o%tdata(it)%valence
          end if
          if (error()) goto 200
        end do

        ! collect type data
        do it = 1,size(ao%o%tdata)
          np = 0
          do il = 1,x_nlcomps(ao%o%tdata(it)%ncpd)
            l = x_lcomp(ao%o%tdata(it)%ncpd,il) ; if (error()) goto 200
            c = x_l_channels(ao%o%tdata(it)%ncpd,il); if (error()) goto 200
            np = np + c*(2*l + 1)
          end do
          if (np == 0) then
            nullify( ao%o%tdata(it)%pl )
            nullify( ao%o%tdata(it)%pm )
            nullify( ao%o%tdata(it)%pc )
            nullify( ao%o%tdata(it)%pkbf )
          else
            allocate( ao%o%tdata(it)%pl(np) )
            allocate( ao%o%tdata(it)%pm(np) )
            allocate( ao%o%tdata(it)%pc(np) )
            allocate( ao%o%tdata(it)%pkbf(np) )
            ip = 0
            do il = 1,x_nlcomps(ao%o%tdata(it)%ncpd)
              l = x_lcomp(ao%o%tdata(it)%ncpd,il) ; if (error()) goto 200
              kbf = x_kbparam(ao%o%tdata(it)%ncpd,il); if (error()) goto 200
              nl = x_l_channels(ao%o%tdata(it)%ncpd,il)
              do m = -l,l
                 do c = 1,nl
                   ip = ip + 1
                   ao%o%tdata(it)%pl(ip) = l
                   ao%o%tdata(it)%pm(ip) = m
                   ao%o%tdata(it)%pc(ip) = c
                   ao%o%tdata(it)%pkbf(ip) = kbf
                 end do  
              end do
            end do
          end if
        end do

        call form_r2c_and_c2r_i(ao%o)

        ! collect atom data
        allocate( ao%o%adata(na) )
        np = 0
        do ia = 1,na
          tag = x_type(x_atoms(ao%o%cr),ia)
          do it = 1,size(ao%o%tdata)
            if (tag == ao%o%tdata(it)%tag) then
              ao%o%adata(ia)%ti = it
              exit
            end if
          end do
          ao%o%adata(ia)%valence = ao%o%tdata(ao%o%adata(ia)%ti)%valence
          ao%o%adata(ia)%bi = np +1
          if (associated( ao%o%tdata(it)%pl )) np = np + size(ao%o%tdata(it)%pl)
        end do
        ! collect projector data
        if (np == 0) then
          nullify( ao%o%pdata )
        else
          allocate( ao%o%pdata(np) )
          ip = 0
          do ia = 1,na
            it = ao%o%adata(ia)%ti
            if (associated( ao%o%tdata(it)%pl )) then
              do itp = 1,size(ao%o%tdata(it)%pl)
                ip = ip + 1
                ao%o%pdata(ip)%ai = ia
                ao%o%pdata(ip)%ti = it
                ao%o%pdata(ip)%tpi = itp
                ao%o%pdata(ip)%l = ao%o%tdata(it)%pl(itp)
                ao%o%pdata(ip)%m = ao%o%tdata(it)%pm(itp)
                ao%o%pdata(ip)%c = ao%o%tdata(it)%pc(itp)
                ao%o%pdata(ip)%kbf = ao%o%tdata(it)%pkbf(itp)
              end do
            end if
          end do
        end if
        ! compute the form factors
        do it = 1,size(ao%o%tdata)
          nullify( ao%o%tdata(it)%vdff )
          nullify( ao%o%tdata(it)%gpff )
          nullify( ao%o%tdata(it)%cdff )
          nullify( ao%o%tdata(it)%sgpff )
          nullify( ao%o%tdata(it)%scdff )
        end do
        call form_factors_i(ao%o,lay) ; if (error()) goto 200
        ! compute the Ewald energy
        allocate( q(na) )
        do ia = 1,na
          q(ia) = ao%o%adata(ia)%valence
        end do
        call ewald_energy(ao%o%cr,q,ao%o%energy)
        ! close the ATOMIC_OPERATORS block
200     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

300      if (allocated( tags )) deallocate( tags )
        if (allocated( q )) deallocate( q )

!        do it = 1,size(ao%o%tdata)
!        call glean(thy(ao%o%tdata(it)%ncpd))
!        end do

        call glean(thy(cr))
        call glean(thy(lay))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_operators_ncp_mod::constructor_ao")) continue
      end function

      subroutine update_ao(ao,cr,lay)
!doc$ subroutine update(ao,cr,lay)
        type(atomic_operators_ncp_obj) :: ao
        type(crystal_obj) :: cr
        type(layout_obj) :: lay
!       modifies: ao
!       requires: Atom types in cr do not change.
!       effects: Updates ao.
!       errors: Passes errors.

!cod$
        logical :: atoms_change, crystal_change, layout_change
        integer :: ia, na
        real(double), dimension(:), allocatable :: q

        call my(ao)
        call my(cr)
        call my(lay)

        atoms_change = ( ao%o%g_atoms /= x_ghost(x_atoms(cr)) )
        crystal_change = ( x_ghost(ao%o%cr) /= x_ghost(cr) )
        layout_change = ( ao%o%g_layout /= x_ghost(lay) )

        if (layout_change .or. crystal_change) then
           call own_i(ao)
           ao%o%g = x_ghost()
           ao%o%g_atoms = x_ghost(x_atoms(cr))
           ao%o%g_layout = x_ghost(lay)
        end if

        if (crystal_change) then
          ao%o%cr = cr
          na = size(ao%o%adata)
          allocate( q(na) )
          do ia = 1,na
            q(ia) = ao%o%adata(ia)%valence
          end do
          call ewald_energy(ao%o%cr,q,ao%o%energy)
        end if

        if (layout_change .or. (crystal_change .and. .not.atoms_change)) then
          call form_factors_i(ao%o,lay) ; if (error()) goto 100
        end if

100     if (allocated( q )) deallocate( q )

        call glean(thy(ao))
        call glean(thy(cr))
        call glean(thy(lay))

        if (error("Exit atomic_operators_ncp_mod::update_ao")) continue

      end subroutine

      subroutine my_ao(ao)
!doc$ subroutine my(ao)
        type(atomic_operators_ncp_obj) :: ao

!cod$
        ao%ref = ao%ref + 1
        ao%o%ref = ao%o%ref + 1
      end subroutine

      subroutine my_new_ao(aoi,ao)
!doc$ subroutine my(aoi,ao)
        type(atomic_operators_ncp_obj) :: aoi
        type(atomic_operators_ncp_obj) :: ao

!cod$
        ao%ref = 1
        ao%o => aoi%o
        ao%o%ref = ao%o%ref + 1
      end subroutine

      function thy_ao(ao) result(aoo)
!doc$ function thy(ao) result(aoo)
        type(atomic_operators_ncp_obj) :: ao, aoo

!cod$
        ao%ref = ao%ref - 1
        ao%o%ref = ao%o%ref - 1
        aoo%ref = ao%ref
        aoo%o => ao%o
      end function

      subroutine glean_ao(ao)
!doc$ subroutine glean(ao)
        type(atomic_operators_ncp_obj) :: ao

!cod$
        integer :: i
        if (ao%o%ref < 1) then
          call glean(thy(ao%o%cr))
          do i = 1,size(ao%o%tdata)
            call glean(thy(ao%o%tdata(i)%ncpd))
            if (associated( ao%o%tdata(i)%pl )) deallocate( ao%o%tdata(i)%pl )
            if (associated( ao%o%tdata(i)%pm )) deallocate( ao%o%tdata(i)%pm )
            if (associated( ao%o%tdata(i)%pc )) deallocate( ao%o%tdata(i)%pc )
            if (associated( ao%o%tdata(i)%pkbf )) deallocate( ao%o%tdata(i)%pkbf )
            if (associated( ao%o%tdata(i)%vdff )) deallocate( ao%o%tdata(i)%vdff )
            if (associated( ao%o%tdata(i)%gpff )) deallocate( ao%o%tdata(i)%gpff )
            if (associated( ao%o%tdata(i)%cdff )) deallocate( ao%o%tdata(i)%cdff )
            if (associated( ao%o%tdata(i)%sgpff )) deallocate( ao%o%tdata(i)%sgpff )
            if (associated( ao%o%tdata(i)%scdff )) deallocate( ao%o%tdata(i)%scdff )
            if (associated( ao%o%tdata(i)%r2c )) deallocate( ao%o%tdata(i)%r2c )
            if (associated( ao%o%tdata(i)%c2r )) deallocate( ao%o%tdata(i)%c2r )
            if (associated( ao%o%tdata(i)%tr2c )) deallocate( ao%o%tdata(i)%tr2c )
            if (associated( ao%o%tdata(i)%tc2r )) deallocate( ao%o%tdata(i)%tc2r )
         end do
          deallocate( ao%o%tdata )
          if (associated( ao%o%adata )) deallocate( ao%o%adata )
          if (associated( ao%o%pdata )) deallocate( ao%o%pdata )
          deallocate( ao%o )
        end if
      end subroutine

      subroutine bequeath_ao(ao)
!doc$ subroutine bequeath(ao)
        type(atomic_operators_ncp_obj) :: ao

!cod$
        continue
      end subroutine
    
      subroutine assign_ao(ao,ao2)
!doc$ subroutine assignment(=)(ao,ao2)
        type(atomic_operators_ncp_obj), intent(inout) :: ao
        type(atomic_operators_ncp_obj), intent(in) :: ao2

!cod$
        type(atomic_operators_ncp_obj) :: aot
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
        type(atomic_operators_ncp_obj) :: ao
        integer, dimension(2) :: r
!       effects: Returns ao%ref and ao%o%ref.

!cod$
        r(1) = ao%ref
        r(2) = ao%o%ref
        call glean(ao)
      end function

      function ao_ghost(ao) result(g)
!doc$ function x_ghost(ao) result(g)
        type(atomic_operators_ncp_obj) :: ao
        type(ghost) :: g
!       effects: Returns the ghost of ao.

!cod$
        call my(ao)
        g = ao%o%g
        call glean(thy(ao))
      end function
 
      function ao_layout_ghost(ao) result(g)
!doc$ function x_layout_ghost(ao) result(g)
        type(atomic_operators_ncp_obj) :: ao
        type(ghost) :: g
!       effects: Returns the layout ghost of ao.

!cod$
        call my(ao)
        g = ao%o%g_layout
        call glean(thy(ao))
      end function

      function ao_crystal(ao) result(cr)
!doc$ function x_crystal(ao) result(cr)
        type(atomic_operators_ncp_obj) :: ao
        type(crystal_obj) :: cr
!       effects: Returns ao%o%cr.

!cod$
        call my(ao)
        call my(ao%o%cr,cr)
        call bequeath(thy(cr))
        call glean(thy(ao))
      end function

      function ao_n_types(ao) result(nt)
!doc$ function x_n_types(ao) result(nt)
        type(atomic_operators_ncp_obj) :: ao
        integer :: nt
!       effects: Returns the number of atom types.

!cod$
        call my(ao)
        nt = size(ao%o%tdata)
        call glean(thy(ao))
      end function

      function ao_n_atoms(ao) result(n)
!doc$ function x_n_atoms(ao) result(n)
        type(atomic_operators_ncp_obj) :: ao
        integer :: n
!       effects: Returns the number of atoms.

!cod$
        call my(ao)
        n = size(ao%o%adata)
        call glean(thy(ao))
      end function

      function ao_type_name(ao,it) result(tag)
!doc$ function x_type_name(ao,it) result(tag)
        type(atomic_operators_ncp_obj) :: ao
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
        type(atomic_operators_ncp_obj) :: ao
        integer :: it
        real(double) :: v
!       effects: Returns the valence of type it.

!cod$
        call my(ao)
        v = ao%o%tdata(it)%valence
        call glean(thy(ao))
      end function

      function ao_valence_electrons(ao) result(ve)
!doc$ function x_valence_electrons(ao) result(ve)
        type(atomic_operators_ncp_obj) :: ao
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

      function ao_n_projectors(ao) result(n)
!doc$ function x_n_projectors(ao) result(n)
        type(atomic_operators_ncp_obj) :: ao
        integer :: n
!       effects: Returns the total number of projectors.

!cod$
        call my(ao)
        if (associated( ao%o%pdata )) then
          n = size(ao%o%pdata)
        else
          n = 0
        end if
        call glean(thy(ao))
      end function

      function ao_n_type_projectors(ao,it) result(n)
!doc$ function x_n_type_projectors(ao,it) result(n)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it
        integer :: n
!       effects: Returns the number of projectors for atom type it.
!       errors: it out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100 
        if (associated( ao%o%tdata(it)%pl )) then
          n = size(ao%o%tdata(it)%pl)
        else
          n = 0
        end if
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_n_type_projectors")) continue
      end function

      function ao_n_atom_projectors(ao,ia) result(n)
!doc$ function x_n_atom_projectors(ao,ia) result(n)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ia
        integer :: n
!       effects: Returns the number of projectors for atom ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100 
        if (associated( ao%o%tdata(ao%o%adata(ia)%ti)%pl )) then
          n = size(ao%o%tdata(ao%o%adata(ia)%ti)%pl)
        else
          n = 0
        end if
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_n_atom_projectors")) continue
      end function

      function ao_projector_l(ao,ip) result(l)
!doc$ function x_projector_l(ao,ip) result(l)
        type(atomic_operators_ncp_obj) :: ao
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
        if (error("Exit atomic_operators_ncp_mod::ao_projector_l")) continue
      end function

      function ao_type_projector_l(ao,it,ip) result(l)
!doc$ function x_projector_l(ao,it,ip) result(l)
        type(atomic_operators_ncp_obj) :: ao
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
        if (error("Exit atomic_operators_ncp_mod::ao_type_projector_l")) continue
      end function

      function ao_projector_m(ao,ip) result(m)
!doc$ function x_projector_m(ao,ip) result(m)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        integer :: m
!       requires: At least one projector exist.
!       effects: Returns the m angular momentum quantum number of projector ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: projector index is out of range")) goto 100
        m = ao%o%pdata(ip)%m
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_m")) continue
      end function

      function ao_type_projector_m(ao,it,ip) result(m)
!doc$ function x_projector_m(ao,it,ip) result(m)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it, ip
        integer :: m
!       requires: At least one projector exist.
!       effects: Returns the l angular momentum quantum number of projector ip of type it.
!       errors: it or ip out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pl)),"ERROR: ip is out of range")) goto 100
        m = ao%o%tdata(it)%pm(ip)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_type_projector_m")) continue
      end function

      function ao_projector_kbf(ao,ip) result(kbf)
!doc$ function x_projector_kbf(ao,ip) result(kbf)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        real(double) :: kbf
!       requires: At least one projector exist.
!       effects:  Returns the Kleinman-Bylander factor of projector ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        kbf = ao%o%pdata(ip)%kbf
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_kbf")) continue
      end function

      function ao_type_projector_kbf(ao,it,ip) result(kbf)
!doc$ function x_projector_kbf(ao,it,ip) result(kbf)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it, ip
        real(double) :: kbf
!       requires: At least one projector exist.
!       effects: Returns the Kleinman-Bylander factor of projector ip of type it.
!       errors: it or ip out of range.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pl)),"ERROR: ip is out of range")) goto 100
        kbf = ao%o%tdata(it)%pkbf(ip)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_type_projector_kbf")) continue
      end function

      function ao_projector_type(ao,ip) result(it)
!doc$ function x_projector_type(ao,ip) result(it)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        integer :: it
!       requires: At least one projector exist.
!       effects: Returns the type index of projector ip.
!       errors: ip out of range.
!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        it = ao%o%pdata(ip)%ti
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_type")) continue
      end function

      function ao_projector_atom(ao,ip) result(ia)
!doc$ function x_projector_atom(ao,ip) result(ia)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        integer :: ia
!       requires: At least one projector exist.
!       effects:  Returns the atom index of projector ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        ia = ao%o%pdata(ip)%ai
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_atom")) continue
      end function


      function ao_projector_index_in_type(ao,ip) result(tpi)
!doc$ function x_projector_index_in_type(ao,ip) result(tpi)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        integer :: tpi
!       requires: At least one projector exist.
!       effects:  Returns the index of projector ip.
!       errors: ip out of range.

!cod$
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        tpi = ao%o%pdata(ip)%tpi
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_index_in_type")) continue
      end function

      function ao_atom_type(ao,ia) result(ti)
!doc$ function x_atom_type(ao,ia) result(ti)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ia
        integer :: ti
!       effects: Returns the type index of atom ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        ti = ao%o%adata(ia)%ti
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_atom_type")) continue
      end function

      function ao_atom_base(ao,ia) result(bi)
!doc$ function x_atom_base(ao,ia) result(bi)
        type(atomic_operators_ncp_obj) :: ao
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
        if (error("Exit atomic_operators_ncp_mod::ao_atom_base")) continue
      end function

      function ao_atom_valence(ao,ia) result(v)
!doc$ function x_atom_valence(ao,ia) result(v)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ia
        real(double) :: v
!       effects: Returns the valence of atom ia.
!       errors: ia out of range.

!cod$
        call my(ao)
        if (error((ia < 1) .or. (ia > size(ao%o%adata)),"ERROR: ia is out of range")) goto 100
        v = ao%o%adata(ia)%valence
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_atom_valence")) continue
      end function

      subroutine ao_apply_projector_kbf(ao,pdots)
!doc$ subroutine apply_projector_kbf(ao,pdots) result(kbf_pdots)
        type(atomic_operators_ncp_obj) :: ao
        complex(double), dimension(:), intent(inout) :: pdots
!       
!       effects: Returns new pdots altered by the respected kbfactors
!       

!cod$
        integer :: ia, il, im, ip, na, it, l, c, start_index, end_index, last_index, i, p_index, p_index_end, ntp
        complex(double), dimension(:), allocatable :: pdots_kbf, pdots_r2c
        call my(ao)
        last_index = 0
        !do i = 1,size(pdots)
        !   if (i_access(diaryfile())) write(x_unit(diaryfile()), '(/,"before kbf pdots", i3, "is", f15.10, f15.10)') i, real(pdots(i)), aimag(pdots(i))
        !end do
        p_index = 1
        do ia = 1,size(ao%o%adata)
           it = ao%o%adata(ia)%ti
           ntp = size(ao%o%tdata(it)%pl)
           p_index_end = p_index+ntp-1
           allocate( pdots_r2c(ntp) )
           pdots_r2c = pdots(p_index:p_index_end)
           pdots(p_index:p_index_end) = matmul(ao%o%tdata(it)%r2c, pdots_r2c)
           do il = 1,x_nlcomps(ao%o%tdata(it)%ncpd)
               l = x_lcomp(ao%o%tdata(it)%ncpd,il)
               c = x_l_channels(ao%o%tdata(it)%ncpd,il)
               do im = -l,l
                  start_index = 1 + last_index
                  end_index = start_index + (c-1)
                  allocate( pdots_kbf(c) )
                  pdots_kbf = pdots(start_index:end_index)
                  call apply_projector_kbf(ao%o%tdata(it)%ncpd, pdots_kbf, il)
                  pdots(start_index:end_index) = pdots_kbf
                  deallocate(pdots_kbf)
                  last_index = end_index
               end do
           end do
           pdots_r2c = pdots(p_index:p_index_end)
           pdots(p_index:p_index_end)= matmul(ao%o%tdata(it)%c2r, pdots_r2c)
           p_index = p_index + ntp
           deallocate(pdots_r2c)
        end do    
        !do i = 1,size(pdots)
        !   if (i_access(diaryfile())) write(x_unit(diaryfile()), '(/,"after kbf pdots", i3, "is", f15.10, f15.10)') i, real(pdots(i)), aimag(pdots(i))
        !end do
        call glean(thy(ao))
      end subroutine

      function ao_atom_energy(ao,wij, ia) result(energy)
!doc$ subroutine x_atom_energy(ao,wij, ia) result(kbf_pdots)
        type(atomic_operators_ncp_obj) :: ao
        complex(double), dimension(:, :), intent(in) :: wij
        integer, intent(in) :: ia
        real(double) :: energy 
!
!       effects: Returns energy associated with atomic density
!

!cod$
        integer :: il, im, ip, na, it, l, c, start_index, end_index, last_index
        complex(double), dimension(:, :), allocatable :: wij_kbf
        complex(double), dimension(:,:), allocatable :: wij_new
        call my(ao)
        last_index = 0
           it = ao%o%adata(ia)%ti
           allocate ( wij_new(size(wij,1), size(wij,2)))
           do il = 1,x_nlcomps(ao%o%tdata(it)%ncpd)
               l = x_lcomp(ao%o%tdata(it)%ncpd,il)
               c = x_l_channels(ao%o%tdata(it)%ncpd,il)
               !wij_new = matmul(ao%o%tdata(it)%r2c, wij)
               do im = -l,l
                  start_index = 1 + last_index
                  end_index = start_index + (c-1)
                  allocate( wij_kbf(c,size(wij,2)) )
                  wij_kbf = wij(start_index:end_index, :)
                  call apply_projector_kbf(ao%o%tdata(it)%ncpd, wij_kbf, il)
                  wij_new(start_index:end_index, :) = wij_kbf
                  deallocate(wij_kbf)
                  last_index = end_index
               end do
           end do
           !wij_new = matmul(ao%o%tdata(it)%c2r, wij_new)
        energy = trace(wij_new)
        deallocate ( wij_new )
        call glean(thy(ao))
      end function


      function ao_projector_f_values(ao,q,ip) result(pfv)
!doc$ function projector_f_values(ao,q,ip) result(pfv)
        type(atomic_operators_ncp_obj) :: ao
        real(double), dimension(:), intent(in)  :: q
        integer, intent(in) :: ip
        real(double), dimension(size(q)) :: pfv
!       requires: At least one projector exist.
!       effects: Returns the Fourier coefficients of projector ip.
!       errors: ip out of range. Passes errors.

!cod$
        integer :: it, l, c
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        it = ao%o%pdata(ip)%ti
        l = ao%o%tdata(it)%pl(ao%o%pdata(ip)%tpi)
        c = ao%o%tdata(it)%pc(ao%o%pdata(ip)%tpi)
        pfv = projector_f_values(ao%o%tdata(it)%ncpd,q,l,c) ; if (error()) goto 100
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_f_values")) continue
      end function

      function ao_type_projector_f_values(ao,q,it,ip) result(pfv)
!doc$ function projector_f_values(ao,q,it,ip,nl) result(pfv)
        type(atomic_operators_ncp_obj) :: ao
        real(double), dimension(:), intent(in)  :: q
        integer, intent(in) :: it, ip
        real(double), dimension(size(q)) :: pfv
!       requires: At least one projector exist.
!       effects: Returns the Fourier coefficients of projector ip of type it.
!       errors: it out of range. ip out of range. Passes errors.

!cod$
        integer :: l, c
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pl)),"ERROR: ip is out of range")) goto 100
        l = ao%o%tdata(it)%pl(ip)
        c = ao%o%tdata(it)%pc(ip)
        pfv = projector_f_values(ao%o%tdata(it)%ncpd,q,l,c) ; if (error()) goto 100
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_type_projector_f_values")) continue
      end function

      function ao_projector_stress_f_values(ao,q,it,ip) result(psfv)
!doc$ function projector_stress_f_values(ao,q,it,ip) result(psfv)
        type(atomic_operators_ncp_obj) :: ao
        real(double), dimension(:), intent(in) :: q
        integer, intent(in) :: it, ip
        real(double), dimension(2,size(q)) :: psfv
!       requires: At least one projector exist.
!       effects: Returns the stress tensor contributions from the Fourier coefficients of projector ip of type it.
!       errors: it out of range. ip out of range. Passes errors.

!cod$
        integer :: l, c
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        if (error((ip < 1) .or. (ip > size(ao%o%tdata(it)%pl)),"ERROR: ip is out of range")) goto 100
        l = ao%o%tdata(it)%pl(ip)
        c = ao%o%tdata(it)%pc(ip)
        psfv = projector_stress_f_values(ao%o%tdata(it)%ncpd,q,l,c) ; if (error()) goto 100
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_stress_f_values")) continue
      end function

      function ao_projector_r_value(ao,ip,r,gi,go) result(prv)
!doc$ function projector_r_value(ao,ip,r,gi,go) result(prv)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        real(double), intent(in) :: r
        real(double), intent(in) :: gi, go
        real(double) :: prv
!       requires: At least one projector exist.
!       effects: Returns the value of projector ip at r.
!       errors: ip out of range. Passes errors.

!cod$
        integer :: l, type
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        type = ao%o%pdata(ip)%ti
        l = ao%o%tdata(type)%pl(ao%o%pdata(ip)%tpi)
        prv = projector_r_value(ao%o%tdata(type)%ncpd,r,l,gi,go) ; if (error()) goto 100
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_r_value")) continue
      end function

      function ao_projector_r_gradients(ao,ip,r) result(prg)
!doc$ function projector_r_gradients(ao,ip,r) result(prg)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ip
        real(double), intent(in) :: r
        real(double), dimension(2) :: prg
!       requires: At least one projector exist.
!       effects: Returns the gradient of projector ip at r.
!       errors: ip out of range. Passes errors.

!cod$
        integer :: l, type
        call my(ao)
        if (error((ip < 1) .or. (ip > size(ao%o%pdata)),"ERROR: ip is out of range")) goto 100
        type = ao%o%pdata(ip)%ti
        l = ao%o%tdata(type)%pl(ao%o%pdata(ip)%tpi)
        prg = projector_r_gradients(ao%o%tdata(type)%ncpd,r,l)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_r_gradients")) continue
      end function

      function ao_projector_radius(ao,ipg) result(radius)
!doc$ function x_projector_radius(ao,ipg) result(radius)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: ipg
        real(double) :: radius
!       requires: At least one projector exist.
!       effects:  Returns the optimization radius of global projector ipg.
!       errors: ipg out of range. Passes errors.

!cod$
        call my(ao)
        if (error((ipg < 1) .or. (ipg > size(ao%o%pdata)),"ERROR: ipg is out of range")) goto 100
        radius = x_radius(ao%o%tdata(ao%o%pdata(ipg)%ti)%ncpd)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_projector_radius")) continue
      end function

      function ao_type_projector_radius(ao,it) result(radius)
!doc$ function x_type_projector_radius(ao,it) result(radius)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it
        real(double) :: radius
!       effects: Returns the optimization radius of type it projectors.
!       errors: it out of range. Passes errors.

!cod$
        call my(ao)
        if (error((it < 1) .or. (it > size(ao%o%tdata)),"ERROR: it is out of range")) goto 100
        radius = x_radius(ao%o%tdata(it)%ncpd)
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::ao_type_projector_radius")) continue
      end function

      function ao_ewald_energy(ao) result(e)
!doc$ function x_ewald_energy(ao) result(e)
        type(atomic_operators_ncp_obj) :: ao
        real(double) :: e
!       effects: Returns the Ewald energy.

!cod$
        call my(ao)
        e = ao%o%energy
        call glean(thy(ao))
      end function



      subroutine type_tr2c(ao,it,p)
!doc$ subroutine type_tr2c(ao,it,p)
        type(atomic_operators_ncp_obj) :: ao
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
        type(atomic_operators_ncp_obj) :: ao
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


      function atomic_grid_potential_ao(ao,lay) result(pot)
!doc$ function atomic_grid_potential(ao,lay) result(pot)
        type(atomic_operators_ncp_obj) :: ao
        type(layout_obj) :: lay
        type(grid_obj) :: pot
!       effects: Returns the filtered atomic grid potential.
!       errors: Passes errors.

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
            igr = (0.0_double,1.0_double)*( pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3) )
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

        if (error("Exit atomic_operators_ncp_mod::atomic_grid_potential_ao")) continue

      end function

      function atomic_xc_density_ao(ao,lay) result(den)
!doc$ function atomic_xc_density(ao,lay) result(den)
        type(grid_obj) :: den
        type(atomic_operators_ncp_obj) :: ao
        type(layout_obj) :: lay
!       effects: Returns the filtered core density or an empty grid if there is no core density.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, it
        real(double) :: spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        complex(double) :: igr
        complex(double), dimension(:,:,:), pointer :: c_den
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(ao)
        call my(lay)

        call my(x_lattice(ao%o%cr),lat)
        call my(x_atoms(ao%o%cr),ats)
        call my(grid(lay,SGROUP),den)

        nullify( gx, gy, gz, c_den )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        if (.not.core_density(ao)) goto 100

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call alloc(c_den,lay,D_TYPE,SGROUP)
        c_den = (0.0_double,0.0_double)
        do ia = 1,x_n_atoms(ao)
          it = x_atom_type(ao,ia)
          pos = lat2r(lat,x_position(ats,ia))
          if (associated(ao%o%tdata(it)%cdff)) then
            do i3 = 1,size(gx,3)
            do i2 = 1,size(gx,2)
            do i1 = 1,size(gx,1)
              igr = (0.0_double,1.0_double)*( pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3) )
              c_den(i1,i2,i3) = c_den(i1,i2,i3) + spin_factor*ao%o%tdata(it)%cdff(i1,i2,i3)*exp(-igr)
            end do
            end do
            end do
          end if
        end do
        call put(c_den,den,CDF_KIND)

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )

        call glean(thy(lat))
        call glean(thy(ats))
        call bequeath(thy(den))

        call glean(thy(ao))
        call glean(thy(lay))

        if (error("Exit atomic_operators_ncp_mod::atomic_xc_density_ao")) continue

      end function

      subroutine diary_rs_projectors_ao(ao)
!doc$ subroutine diary_rs_projectors(ao)
        type(atomic_operators_ncp_obj) :: ao
!       requires: Real-space projectors be in use.
!       effects: Writes real-space projector information to the diary.

!cod$
        integer :: it
        call my(ao)
        do it = 1,size(ao%o%tdata)
          call diary_rs_projectors(ao%o%tdata(it)%ncpd)
        end do
        call glean(thy(ao))
      end subroutine

      function core_density_ao(ao) result(cd)
!doc$ function core_density(ao,it) result(cd)
        type(atomic_operators_ncp_obj) :: ao
        logical :: cd
!       effects: Returns .true. iff any type has core density.

!cod$
        integer :: it
        call my(ao)
        cd = .false.
        do it = 1,size(ao%o%tdata)
          if (associated(ao%o%tdata(it)%cdff)) then
            cd = .true.
            exit
          end if
        end do
100     call glean(thy(ao))
        if (error("Exit atomic_operators_ncp_mod::core_density_ao")) continue
      end function

      subroutine valence_density_ff(ao,it,vdff)
!doc$ subroutine valence_density_ff(ao,it,vdff)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:), pointer :: vdff
!       requires: it in range.
!       modifies: vdff
!       effects: Points vdff to the valence density form factors of type it.

!cod$
        call my(ao)
        vdff => ao%o%tdata(it)%vdff
        call glean(thy(ao))
      end subroutine

      subroutine grid_potential_ff(ao,it,gpff)
!doc$ subroutine grid_potential_ff(ao,it,gpff)
        type(atomic_operators_ncp_obj) :: ao
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

      subroutine core_density_ff(ao,it,cdff)
!doc$ subroutine core_density_ff(ao,it,cdff)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:), pointer :: cdff
!       requires: it in range.
!       modifies: cdff
!       effects: Points cdff to the core density form factors of type it or nullifies cdff if they do not exist.

!cod$
        call my(ao)
        if (associated( ao%o%tdata(it)%cdff )) then
          cdff => ao%o%tdata(it)%cdff
        else
          nullify( cdff )
        end if
        call glean(thy(ao))
      end subroutine

      subroutine stress_grid_potential_ff(ao,it,sgpff)
!doc$ subroutine stress_grid_potential_ff(ao,it,sgpff)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:), pointer :: sgpff
!       requires: it in range. Stress grid potential form factors exist.
!       modifies: sgpff
!       effects: Points sgpff to the stress grid potential form factors of type it.

!cod$
        call my(ao)
        sgpff => ao%o%tdata(it)%sgpff
        call glean(thy(ao))
      end subroutine

      subroutine stress_core_density_ff(ao,it,scdff)
!doc$ subroutine stress_core_density_ff(ao,it,scdff)
        type(atomic_operators_ncp_obj) :: ao
        integer, intent(in) :: it
        real(double), dimension(:,:,:), pointer :: scdff
!       requires: it in range. Stress core density form factors exist.
!       modifies: den
!       effects: Points scdff to the stress core density form factors of type it or nullifies scdff if they do not exist.

!cod$
        call my(ao)
        if (associated( ao%o%tdata(it)%scdff )) then
          scdff => ao%o%tdata(it)%scdff
        else
          nullify( scdff )
        end if
        call glean(thy(ao))
      end subroutine

      subroutine write_restart_ao(ao,nrestf)
!doc$ subroutine write_restart(ao,nrestf)
        type(atomic_operators_ncp_obj) :: ao
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ao restart information to nrestf.

!cod$
        integer :: nt
        integer(long) :: dsize, iosl, ndata, s4

        call my(ao)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the ATOMIC_OPERATORS block
          call startblock(nrestf,"ATOMIC_OPERATORS")

          ! write the number of types
          call writetag(nrestf,"NUMBER_OF_TYPES")
          nt = size(ao%o%tdata)
          s4 = nt
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! end the ATOMIC_OPERATORS block
          call endblock(nrestf)

        end if

        call glean(thy(ao))
        call glean(thy(nrestf))

        if (error("Exit atomic_operators_ncp_mod::write_restart_ao")) continue

      end subroutine

! private routines

      subroutine own_i(ao)
        type(atomic_operators_ncp_obj) :: ao
        type(atomic_operators_ncp_obj) :: aot
        integer :: i
        if (ao%ref < ao%o%ref) then
          allocate( aot%o )
          aot%o%ref = 0
          aot%o%g_atoms = ao%o%g_atoms
          aot%o%g_layout = ao%o%g_layout
          aot%o%energy = ao%o%energy
          call my(ao%o%cr,aot%o%cr)
          allocate( aot%o%tdata(size(ao%o%tdata)) )
          do i = 1,size(aot%o%tdata)
            call my(ao%o%tdata(i)%ncpd,aot%o%tdata(i)%ncpd)
            aot%o%tdata(i)%tag = ao%o%tdata(i)%tag
            aot%o%tdata(i)%valence = ao%o%tdata(i)%valence
            if (associated( ao%o%tdata(i)%pl )) then
              allocate( aot%o%tdata(i)%pl(size(ao%o%tdata(i)%pl)) ) ; aot%o%tdata(i)%pl = ao%o%tdata(i)%pl
              allocate( aot%o%tdata(i)%pm(size(ao%o%tdata(i)%pm)) ) ; aot%o%tdata(i)%pm = ao%o%tdata(i)%pm
              allocate( aot%o%tdata(i)%pc(size(ao%o%tdata(i)%pc)) ) ; aot%o%tdata(i)%pc = ao%o%tdata(i)%pc
              allocate( aot%o%tdata(i)%pkbf(size(ao%o%tdata(i)%pkbf)) ) ; aot%o%tdata(i)%pkbf = ao%o%tdata(i)%pkbf
            else
              nullify( aot%o%tdata(i)%pl )
              nullify( aot%o%tdata(i)%pm )
              nullify( aot%o%tdata(i)%pc )
              nullify( aot%o%tdata(i)%pkbf )
            end if
            allocate( aot%o%tdata(i)%vdff(size(ao%o%tdata(i)%vdff,1),size(ao%o%tdata(i)%vdff,2),size(ao%o%tdata(i)%vdff,3)) )
            aot%o%tdata(i)%vdff = ao%o%tdata(i)%vdff
            allocate( aot%o%tdata(i)%gpff(size(ao%o%tdata(i)%gpff,1),size(ao%o%tdata(i)%gpff,2),size(ao%o%tdata(i)%gpff,3)) )
            aot%o%tdata(i)%gpff = ao%o%tdata(i)%gpff
            if (associated( ao%o%tdata(i)%cdff )) then
              allocate( aot%o%tdata(i)%cdff(size(ao%o%tdata(i)%cdff,1),size(ao%o%tdata(i)%cdff,2),size(ao%o%tdata(i)%cdff,3)) )
              aot%o%tdata(i)%cdff = ao%o%tdata(i)%cdff
            else
              nullify( aot%o%tdata(i)%cdff )
            end if
            if (associated( ao%o%tdata(i)%sgpff )) then
              allocate( aot%o%tdata(i)%sgpff(size(ao%o%tdata(i)%sgpff,1),size(ao%o%tdata(i)%sgpff,2),size(ao%o%tdata(i)%sgpff,3)) )
              aot%o%tdata(i)%sgpff = ao%o%tdata(i)%sgpff
            else
              nullify( aot%o%tdata(i)%sgpff )
            end if
            if (associated( ao%o%tdata(i)%scdff )) then
              allocate( aot%o%tdata(i)%scdff(size(ao%o%tdata(i)%scdff,1),size(ao%o%tdata(i)%scdff,2),size(ao%o%tdata(i)%scdff,3)) )
              aot%o%tdata(i)%scdff = ao%o%tdata(i)%scdff
            else
              nullify( aot%o%tdata(i)%scdff )
            end if
            if (associated( ao%o%tdata(i)%r2c )) then
               allocate( aot%o%tdata(i)%r2c(size(ao%o%tdata(i)%r2c,1),size(ao%o%tdata(i)%r2c,2)))
               allocate( aot%o%tdata(i)%tr2c(size(ao%o%tdata(i)%r2c,1),size(ao%o%tdata(i)%r2c,2)))
               allocate( aot%o%tdata(i)%c2r(size(ao%o%tdata(i)%r2c,1),size(ao%o%tdata(i)%r2c,2)))
               allocate( aot%o%tdata(i)%tc2r(size(ao%o%tdata(i)%r2c,1),size(ao%o%tdata(i)%r2c,2)))
               aot%o%tdata(i)%r2c = ao%o%tdata(i)%r2c
               aot%o%tdata(i)%tr2c = ao%o%tdata(i)%tr2c
               aot%o%tdata(i)%c2r = ao%o%tdata(i)%c2r
               aot%o%tdata(i)%tc2r = ao%o%tdata(i)%tc2r
             else
                nullify( aot%o%tdata(i)%r2c )
                nullify( aot%o%tdata(i)%c2r )
                nullify( aot%o%tdata(i)%tr2c )
                nullify( aot%o%tdata(i)%tc2r )
             end if
          end do
          allocate( aot%o%adata(size(ao%o%adata)) )
          do i = 1,size(aot%o%adata)
            aot%o%adata(i)%ti = ao%o%adata(i)%ti
            aot%o%adata(i)%bi = ao%o%adata(i)%bi
            aot%o%adata(i)%valence = ao%o%adata(i)%valence
          end do
          if (associated( ao%o%pdata )) then
            allocate( aot%o%pdata(size(ao%o%pdata)) )
            do i = 1,size(aot%o%pdata)
              aot%o%pdata(i)%ai = ao%o%pdata(i)%ai
              aot%o%pdata(i)%ti = ao%o%pdata(i)%ti
              aot%o%pdata(i)%tpi = ao%o%pdata(i)%tpi
              aot%o%pdata(i)%l = ao%o%pdata(i)%l
              aot%o%pdata(i)%m = ao%o%pdata(i)%m
              aot%o%pdata(i)%c = ao%o%pdata(i)%c
              aot%o%pdata(i)%kbf = ao%o%pdata(i)%kbf
            end do
          else
            nullify( aot%o%pdata )
          end if
          aot%o%g = ao%o%g
          ao%o%ref = ao%o%ref - ao%ref
          ao%o => aot%o
          ao%o%ref = ao%o%ref + ao%ref
        end if
      end subroutine

      subroutine form_r2c_and_c2r_i(aor)
        type(atomic_operators_ncp_rep) :: aor

        integer :: ic, ir, it, mc, mr, lc, np, lr, cr, cc
        real(double), parameter :: minus_one = -1.0_double
        real(double) :: norm
        complex(double), parameter :: i = (0.0_double,1.0_double)

        norm = 1.0_double/sqrt(2.0_double)

        do it = 1,size(aor%tdata)
          np = size(aor%tdata(it)%pl)
          allocate( aor%tdata(it)%r2c(np,np), aor%tdata(it)%c2r(np,np), aor%tdata(it)%tr2c(np,np), aor%tdata(it)%tc2r(np,np) )
          do ic = 1,np
            do ir = 1,np
              lr = aor%tdata(it)%pl(ir)
              lc = aor%tdata(it)%pl(ic)
              mr = aor%tdata(it)%pm(ir)
              mc = aor%tdata(it)%pm(ic)
              cr = aor%tdata(it)%pc(ir)
              cc = aor%tdata(it)%pc(ic)
              if ((lr == lc) .and. (abs(mr) == abs(mc)) .and. (cr == cc)) then
                if ((mr < 0) .and. (mc < 0)) then
                  aor%tdata(it)%r2c(ir,ic) = -i*norm
                elseif ((mr < 0) .and. (mc > 0)) then
                  aor%tdata(it)%r2c(ir,ic) = norm
                elseif ((mr == 0) .and. (mc == 0)) then
                  aor%tdata(it)%r2c(ir,ic) = 1.0_double
                elseif ((mr > 0) .and. (mc < 0)) then
                  aor%tdata(it)%r2c(ir,ic) = i*minus_one**abs(mr)*norm
                elseif ((mr > 0) .and. (mc > 0)) then
                  aor%tdata(it)%r2c(ir,ic) = norm*minus_one**abs(mr)
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
        type(atomic_operators_ncp_rep) :: aor
        type(layout_obj) :: lay
!       requires: Relevant pointers be nullified or associated.

        real(double), parameter :: gauss_width = 1.0_double

        logical :: found, pressure, stress_tensor, lattice_relaxation
        character(line_len) :: tag
        integer :: i1, i2, i3, it, n1, n2, n3
        real(double) :: cutoff, gm, volume
        real(double), dimension(:,:,:), pointer :: g2

        call my(lay)

        nullify( g2 )

        ! Determine if pressure and stress-tensor related coefficients are needed.
        pressure = .false.
        call arg("pressure",tag,found)
        if (.not.found) tag = "OFF"
        select case (trim(tag))
        case ("ON","On","on",".TRUE.",".true.")
          pressure = .true.
        end select
        stress_tensor = .false.
        call arg("stress_tensor",tag,found)
        if (.not.found) tag = "OFF"
        select case (trim(tag))
        case ("ON","On","on",".TRUE.",".true.")
          stress_tensor = .true.
        end select
        call arg("lattice_relaxation",tag,found)
         if (.not.found) tag = "OFF"
         select case (trim(tag))
         case ( "none", "NONE", "None", "off", "OFF", "Off", "FALSE", "false", "False" )
            lattice_relaxation = .false.
         case ( "shape", "SHAPE", "Shape" )
            lattice_relaxation = .true.
         case ( "shape_quench", "SHAPE_QUENCH", "Shape_quench" )
            lattice_relaxation = .true.
         case ( "shape_quench_z", "SHAPE_QUENCH_Z", "Shape_quench_z")
            lattice_relaxation = .true.
         case default
           if (error(.true.,"ERROR: unrecognized lattice_relaxation")) goto 100
         end select


        cutoff = x_cutoff(lay)
        volume = x_cell_volume(x_lattice(aor%cr))

        call fdel(g2,lay,D_TYPE,SGROUP)

        n1 = size(g2,1)
        n2 = size(g2,2)
        n3 = size(g2,3)

        do it = 1,size(aor%tdata)

          if (associated( aor%tdata(it)%vdff )) deallocate( aor%tdata(it)%vdff )
          if (associated( aor%tdata(it)%gpff )) deallocate( aor%tdata(it)%gpff )
          if (associated( aor%tdata(it)%cdff )) deallocate( aor%tdata(it)%cdff )
          if (associated( aor%tdata(it)%sgpff )) deallocate( aor%tdata(it)%sgpff )
          if (associated( aor%tdata(it)%scdff )) deallocate( aor%tdata(it)%scdff )

          allocate( aor%tdata(it)%vdff(n1,n2,n3) )
          allocate( aor%tdata(it)%gpff(n1,n2,n3) )
          if (core_density(aor%tdata(it)%ncpd)) then
            allocate( aor%tdata(it)%cdff(n1,n2,n3) )
          end if
          if (pressure .or. stress_tensor .or. lattice_relaxation) then
            allocate( aor%tdata(it)%sgpff(n1,n2,n3) )
            if (core_density(aor%tdata(it)%ncpd)) then
              allocate( aor%tdata(it)%scdff(n1,n2,n3) )
            end if
          end if

          ! construct the form factors
          do i3 = 1,n3
          do i2 = 1,n2
          do i1 = 1,n1
            if (g2(i1,i2,i3) > cutoff) then
              aor%tdata(it)%vdff(i1,i2,i3) = 0.0_double
              aor%tdata(it)%gpff(i1,i2,i3) = 0.0_double
              if (core_density(aor%tdata(it)%ncpd)) then
                aor%tdata(it)%cdff(i1,i2,i3) = 0.0_double
              end if
              if (pressure .or. stress_tensor .or. lattice_relaxation) then
                aor%tdata(it)%sgpff(i1,i2,i3) = 0.0_double
                if (core_density(aor%tdata(it)%ncpd)) then
                  aor%tdata(it)%scdff(i1,i2,i3) = 0.0_double
                end if
              end if
            else
              gm = sqrt(g2(i1,i2,i3))
              if (x_valence(aor%tdata(it)%ncpd)) then
                aor%tdata(it)%vdff(i1,i2,i3) = valence_f_value(aor%tdata(it)%ncpd,gm)/volume
              else
                aor%tdata(it)%vdff(i1,i2,i3) = exp(-gauss_width*g2(i1,i2,i3))/volume
              end if
              aor%tdata(it)%gpff(i1,i2,i3) = local_f_value(aor%tdata(it)%ncpd,gm,aor%tdata(it)%valence)/volume
              if (core_density(aor%tdata(it)%ncpd)) then
                aor%tdata(it)%cdff(i1,i2,i3) = core_f_value(aor%tdata(it)%ncpd,gm)/volume
              end if
              if (pressure .or. stress_tensor .or. lattice_relaxation) then
                aor%tdata(it)%sgpff(i1,i2,i3) = local_stress_f_value(aor%tdata(it)%ncpd,gm,aor%tdata(it)%valence)/volume
                if (core_density(aor%tdata(it)%ncpd)) then
                  aor%tdata(it)%scdff(i1,i2,i3) = core_stress_f_value(aor%tdata(it)%ncpd,gm)/volume
                end if
              end if
            end if
          end do
          end do
          end do

        end do

100     if (associated( g2 )) deallocate( g2 )

        call glean(thy(lay))

        if (error("Exit atomic_operators_ncp_mod::form_factors_i")) continue

      end subroutine

      end module

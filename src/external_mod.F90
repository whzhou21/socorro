! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module external_mod
!doc$ module external_mod

!     One data type is available here: type(external_obj).

!     external_mod contains the "exterior" components of a configuration.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use diary_mod
      use ghost_mod
      use arg_mod
      use layout_mod
      use grid_mod
      use atoms_mod
      use crystal_mod
      use symmetry_mod
      use xc_type_mod
      use atomic_operators_mod

!cod$
      implicit none
      private

      integer, parameter :: INPUT_CRYSTAL = 1
      integer, parameter :: RESTART_FILE  = 2

      type external_rep
        integer :: ref
        type(ghost) :: g
        integer :: init_method
        type(crystal_obj) :: crystal
        type(layout_obj) :: layout
        type(xc_type_obj) :: xc_type
        type(atomic_operators_obj) :: ao
        type(point_group_obj) :: lattice_group
        type(space_group_obj) :: space_group
        type(point_group_obj) :: double_group
      end type

      type, public :: external_obj
        integer :: ref
        type(external_rep), pointer :: o
      end type
 
!doc$
      public :: external
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_crystal
      public :: x_layout
      public :: x_xc_type
      public :: x_atomic_operators
      public :: x_lattice_group
      public :: x_space_group
      public :: x_double_group
      public :: valence_electrons
      public :: symmetrize_vectors
      public :: symmetrize_tensor
      public :: symmetrize_grid
      public :: diary
      public :: write_restart

!cod$
      interface external
        module procedure constructor_ext
      end interface
      interface update
        module procedure update_ext
      end interface
      interface my
        module procedure my_ext, my_new_ext
      end interface
      interface thy
        module procedure thy_ext
      end interface
      interface glean
        module procedure glean_ext
      end interface
      interface bequeath
        module procedure bequeath_ext
      end interface
      interface assignment(=)
        module procedure assign_ext
      end interface
      interface x_ref
        module procedure ext_ref
      end interface
      interface x_ghost
        module procedure ext_ghost
      end interface
      interface x_crystal
        module procedure ext_crystal
      end interface
      interface x_layout
        module procedure ext_layout
      end interface
      interface x_xc_type
        module procedure ext_xc_type
      end interface
      interface x_atomic_operators
        module procedure ext_atomic_operators
      end interface
      interface x_lattice_group
        module procedure ext_lattice_group
      end interface
      interface x_space_group
        module procedure ext_space_group
      end interface
      interface x_double_group
        module procedure ext_double_group
      end interface
      interface valence_electrons
        module procedure q_valence_electrons
      end interface
      interface symmetrize_vectors
        module procedure symmetrize_vectors_ext
      end interface
      interface symmetrize_tensor
        module procedure symmetrize_tensor_ext
      end interface
      interface symmetrize_grid
        module procedure symmetrize_grid_ext
      end interface
      interface diary
        module procedure diary_ext, diary_crystal_step
      end interface
      interface write_restart
        module procedure write_restart_ext
      end interface

      contains

! public routines

      function constructor_ext(restf,cr) result(ext)
!doc$ function external(restf,cr) result(ext)
        type(tagio_obj), optional :: restf
        type(crystal_obj), optional :: cr
        type(external_obj) :: ext
!       effects: Constructs a new ext.
!       errors: Passes errors.

!cod$
        logical :: found, sym_atoms
        character(1) :: tios
        integer :: ia
        real(double), dimension(:,:), allocatable :: pos
        type(atoms_obj) :: at_tmp

        if (error("  Error on entry")) then
          ext%ref = 0
          allocate( ext%o )
          ext%o%ref = 0
          goto 999
        end if

        if (present(restf)) call my(restf)
        if (present(cr)) call my(cr)

        ext%ref = 0
        allocate( ext%o )
        ext%o%ref = 0
        ext%o%g = x_ghost()

        ! open the EXTERNAL block
        if (present(restf)) then
          if (error(present(cr),"ERROR: optional arguments restf and cr are incompatible")) goto 300
          if (i_access(restf)) tios = findfirsttag(restf,"EXTERNAL")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: EXTERNAL block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)
        end if

        ! form the crystal, layout, and xc_type
        if (present(restf)) then
          ext%o%init_method = RESTART_FILE
          call my(crystal(restf),ext%o%crystal) ; if (error()) goto 200
          call my(layout(x_lattice(ext%o%crystal),restf=restf),ext%o%layout) ; if (error()) goto 200
          call my(xc_type(restf),ext%o%xc_type) ; if (error()) goto 200
        else
          ext%o%init_method = INPUT_CRYSTAL
          if (present(cr)) then
             call my(cr,ext%o%crystal)
          else
             call my(crystal(),ext%o%crystal)
          end if
          call my(layout(x_lattice(ext%o%crystal)),ext%o%layout) ; if (error()) goto 200
          call my(xc_type(),ext%o%xc_type) ; if (error()) goto 300
        end if

        ! form the lattice group
        call my(point_group(x_lattice(ext%o%crystal)),ext%o%lattice_group) ; if (error()) goto 200

        ! form the space group
        call my(space_group(ext%o%lattice_group,x_atoms(ext%o%crystal), &
                & x_lattice(ext%o%crystal),ext%o%layout),ext%o%space_group) ; if (error()) goto 200

        ! symmetrize the atom positions
        call arg("symmetrize_atoms",sym_atoms,found)
        if (.not.found) sym_atoms = .false.
        if (sym_atoms) then
          call my(x_atoms(ext%o%crystal),at_tmp)
          allocate( pos(3,x_n_atoms(x_atoms(ext%o%crystal))) )
          do ia = 1,x_n_atoms(x_atoms(ext%o%crystal))
            pos(:,ia) = x_position(x_atoms(ext%o%crystal),ia)
          end do
          call symmetrize_coordinates(ext%o%space_group,pos) ; if (error()) goto 200
          call move(at_tmp,pos) ; if (error()) goto 100
          call update(ext%o%crystal,at_tmp) ; if (error()) goto 100
          call update(ext%o%space_group,ext%o%lattice_group,x_atoms(ext%o%crystal),x_lattice(ext%o%crystal),ext%o%layout)
100       call glean(thy(at_tmp))
          if (error()) goto 200
        end if

        ! form the double group
        call my(point_group(ext%o%space_group,parity=.true.),ext%o%double_group) ; if (error()) goto 200

        ! form atomic_operators
        if (present(restf)) then
          call my(atomic_operators(ext%o%crystal,ext%o%layout,ext%o%space_group, &
                  & ext%o%xc_type,restf),ext%o%ao) ; if (error()) goto 200
        else
           call my(atomic_operators(ext%o%crystal,ext%o%layout,ext%o%space_group, &
                  & ext%o%xc_type),ext%o%ao) ; if (error()) goto 200
        end if

        ! flag unsupported PAW capability
        select case (x_type(ext%o%ao))
        case (PAW)
          select case (x_functional_dependence(ext%o%xc_type))
          case (FD_ORBITAL,FD_HYBRID)
            if (error(.true.,"ERROR: hybrid and orbital-dependent functionals are not supported in the PAW method")) goto 200
          end select
        end select

        call diary_construction_i(ext%o,sym_atoms)

        ! close the EXTERNAL block
200     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)           
        end if

300     if (allocated( pos )) deallocate( pos )

        if (present(cr)) call glean(thy(cr))
        if (present(restf)) call glean(thy(restf))

999     if (error("Exit external_mod::constructor_ext")) continue

      end function

      subroutine update_ext(ext,cr)
!doc$ subroutine update(ext,cr)
        type(external_obj) :: ext
        type(crystal_obj), optional :: cr
!       modifies: ext
!       effects: Updates ext.
!       errors: Inconsistent lattices. Passes errors.

!cod$ 
        logical :: crystal_change
        type(ghost) :: g_old
        call my(ext)
        if (present(cr)) then
          call my(cr)
          if (error(x_ghost(x_lattice(cr)) /= x_ghost(x_lattice(ext%o%layout)),"ERROR: inconsistent lattices")) goto 100
          crystal_change = ( x_ghost(ext%o%crystal) /= x_ghost(cr) )
        else
          crystal_change = .false.
        end if
        if (crystal_change) then
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Updating external:")')
          call own_i(ext)
          ext%o%g = x_ghost()
          g_old = x_ghost(ext%o%lattice_group)
          call update(ext%o%lattice_group,x_lattice(cr)) ; if (error()) goto 100
          if (x_ghost(ext%o%lattice_group) /= g_old) call diary(ext%o%lattice_group)
          g_old = x_ghost(ext%o%space_group)
          call update(ext%o%space_group,ext%o%lattice_group,x_atoms(cr),x_lattice(cr),ext%o%layout) ; if (error()) goto 100
          if (x_ghost(ext%o%space_group) /= g_old) then
            call diary(ext%o%space_group)
            call update(ext%o%double_group,ext%o%space_group,parity=.true.) ; if (error()) goto 100
          end if
          ext%o%crystal = cr
          call update(ext%o%ao,ext%o%crystal,ext%o%layout,ext%o%space_group)
        end if
100     if (present(cr)) call glean(thy(cr))
        call glean(thy(ext))
        if (error("Exit external_mod::update_ext")) continue
      end subroutine

      subroutine my_ext(ext)
!doc$ subroutine my(ext)
        type(external_obj) :: ext

!cod$
        ext%ref = ext%ref + 1
        ext%o%ref = ext%o%ref + 1
      end subroutine

      subroutine my_new_ext(exti,ext)
!doc$ subroutine my(exti,ext)
        type(external_obj) :: exti, ext

!cod$
        ext%ref = 1
        ext%o => exti%o
        ext%o%ref = ext%o%ref + 1
      end subroutine

      function thy_ext(ext) result (exto)
!doc$ function thy(ext) result(exto)
        type(external_obj) :: ext, exto

!cod$
        ext%ref = ext%ref - 1
        ext%o%ref = ext%o%ref - 1
        exto%ref = ext%ref
        exto%o => ext%o
      end function
      
      subroutine glean_ext(ext)
!doc$ subroutine glean(ext)
        type(external_obj) :: ext

!cod$
        if (ext%o%ref < 1) then
          call glean(thy(ext%o%crystal))
          call glean(thy(ext%o%layout))
          call glean(thy(ext%o%xc_type))
          call glean(thy(ext%o%ao))
          call glean(thy(ext%o%lattice_group))
          call glean(thy(ext%o%space_group))
          call glean(thy(ext%o%double_group))
          deallocate( ext%o )
        end if
      end subroutine

      subroutine bequeath_ext(ext)
!doc$ subroutine bequeath(ext)

!cod$
        type(external_obj) :: ext
        continue
      end subroutine

      subroutine assign_ext(ext,ext2)
!doc$ subroutine assignment(=)(ext,ext2)
        type(external_obj), intent(inout) :: ext
        type(external_obj), intent(in) :: ext2

!cod$
        type(external_obj) extt
        call my(ext2)
        extt%o => ext%o
        ext%o%ref = ext%o%ref - ext%ref
        ext%o => ext2%o
        ext%o%ref = ext%o%ref + ext%ref
        call glean(extt)
        call glean(thy(ext2))
      end subroutine

      function ext_ref(ext) result(r)
!doc$ function x_ref(ext) result(r)
        type(external_obj) :: ext
        integer, dimension(2) :: r
!       effects: Returns ext%ref and ext%o%ref.

!cod$
        r(1) = ext%ref
        r(2) = ext%o%ref
        call glean(ext)
      end function

      function ext_ghost(ext) result(g)
!doc$ function x_ghost(ext) result(g)
        type(external_obj) :: ext
        type(ghost) :: g
!       effects: Returns ext%o%g.

!cod$
        call my(ext)
        g = ext%o%g
        call glean(thy(ext))
      end function

      function ext_crystal(ext) result(cr)
!doc$ function x_crystal(ext) result(cr)
        type(external_obj) :: ext
        type(crystal_obj) :: cr
!       effects: Returns ext%o%crystal.

!cod$
        call my(ext)
        call my(ext%o%crystal,cr)
        call glean(thy(ext))
        call bequeath(thy(cr))
      end function

      function ext_layout(ext) result(lay)
!doc$ function x_layout(ext) result(lay)
        type(external_obj) :: ext
        type(layout_obj) :: lay
!       effects: Returns ext%o%layout.
!doc$
        call my(ext)
        call my(ext%o%layout,lay)
        call glean(thy(ext))
        call bequeath(thy(lay))
      end function

      function ext_xc_type(ext) result(xct)
!doc$ function x_xc_type(ext) result(xct)
        type(external_obj) :: ext
        type(xc_type_obj) :: xct
!       effects: Returns ext%o%xc_type.
!doc$
        call my(ext)
        call my(ext%o%xc_type,xct)
        call glean(thy(ext))
        call bequeath(thy(xct))
      end function

      function ext_atomic_operators(ext) result(ao)
!doc$ function x_atomic_operators(ext) result(ao)
        type(external_obj) :: ext
        type(atomic_operators_obj) :: ao
!       effects: Returns ext%o%ao.
!       errors: Passes errors.

!cod$
        call my(ext)
        call my(ext%o%ao,ao)
        call glean(thy(ext))
        call bequeath(thy(ao))
        if (error("Exit external_mod::ext_atomic_operators")) continue
      end function

      function ext_lattice_group(ext) result(lg)
!doc$ function x_lattice_group(ext) result(lg)
        type(external_obj) :: ext
        type(point_group_obj) :: lg
!       effects: Returns ext%o%lattice_group.

!cod$
        call my(ext)
        call my(ext%o%lattice_group,lg)
        call bequeath(thy(lg))
        call glean(thy(ext))
      end function

      function ext_space_group(ext) result(sg)
!doc$ function x_space_group(ext) result(sg)
        type(external_obj) :: ext
        type(space_group_obj) :: sg
!       effects: Returns a pointer to ext%o%space_group.

!cod$
        call my(ext)
        call my(ext%o%space_group,sg)
        call glean(thy(ext))
        call bequeath(thy(sg))
      end function

      function ext_double_group(ext) result(dg)
!doc$ function x_double_group(ext) result(dg)
        type(external_obj) :: ext
        type(point_group_obj) :: dg
!       effects: Returns a pointer to ext%o%double_group.

!cod$
        call my(ext)
        call my(ext%o%double_group,dg)
        call glean(thy(ext))
        call bequeath(thy(dg))
      end function

      function q_valence_electrons(ext) result(ve)
!doc$ function valence_electrons(ext) result(ve)
        type(external_obj) :: ext
        real(double) :: ve
!       effects: Returns the number of valence electrons.

!cod$
        call my(ext)
        ve = x_valence_electrons(ext%o%ao)
        call glean(thy(ext))
      end function

      subroutine symmetrize_vectors_ext(ext,v)
!doc$ subroutine symmetrize_vectors(ext,v)
        type(external_obj) :: ext
        real(double), dimension(:,:), intent(inout) :: v
!       modifies : v
!       effects: Overwrites vectors in v with symmetrized versions.
!       errors: Passes errors.

!cod$
        call my(ext)
        call symmetrize_vectors(ext%o%space_group,v) ; if (error()) goto 100
100     call glean(thy(ext))
        if (error("Exit external_mod::symmetrize_vectors_ext")) continue
      end subroutine
     
      subroutine symmetrize_tensor_ext(ext,t)
!doc$ subroutine symmetrize_tensor(ext,t)
        type(external_obj) :: ext
        real(double), dimension(3,3), intent(inout) :: t
!       modifies : t
!       effects: Overwrites t with a symmetrized version.
!       errors: Passes errors.

!cod$
        call my(ext)
        call symmetrize_tensor(ext%o%space_group,t) ; if (error()) goto 100
100     call glean(thy(ext))
        if (error("Exit external_mod::symmetrize_tensor_ext")) continue
      end subroutine
     
      subroutine symmetrize_grid_ext(ext,g)
!doc$ subroutine symmetrize_grid(ext,g)
        type(external_obj) :: ext
        type(grid_obj) :: g
!       modifies : g
!       effects: Overwrites g with a symmetrized version.
!       errors: Passes errors.
!doc$
        call my(ext)
        call symmetrize_grid(ext%o%space_group,g) ; if (error()) goto 100
100     call glean(thy(ext))
        if (error("Exit external_mod::symmetrize_grid_ext")) continue
      end subroutine

      subroutine write_restart_ext(ext,nrestf)
!doc$ subroutine write_restart(ext,nrestf)
        type(external_obj) :: ext
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ext information to nrestf.
!       errors: Passes errors.

!cod$
        call my(ext)
        call my(nrestf)

        ! start the EXTERNAL block
        if (i_access(nrestf)) call startblock(nrestf,"EXTERNAL")

        ! write crystal information
        call write_restart(ext%o%crystal,nrestf) ; if (error()) goto 100

        ! write layout information
        call write_restart(ext%o%layout,nrestf) ; if (error()) goto 100

        ! write xc_type information
        call write_restart(ext%o%xc_type,nrestf) ; if (error()) goto 100

        ! write atomic_operators information
        call write_restart(ext%o%ao,nrestf) ; if (error()) goto 100

        ! end the EXTERNAL block
        if (i_access(nrestf)) call endblock(nrestf)

100     call glean(thy(ext))
        call glean(thy(nrestf))

        if (error("Exit external_mod::write_restart_ext")) continue

      end subroutine

      subroutine diary_ext(ext)
!doc$ subroutine diary(ext)
        type(external_obj) :: ext
!       modifies : diary stream
!       effects: Writes ext information to the diary.
!       errors: Passes errors.

!cod$
        call my(ext)
        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"External:")')
        call diary(ext%o%crystal)
        call diary(ext%o%lattice_group,lat=x_lattice(ext%o%crystal))
        call diary(ext%o%space_group,lat=x_lattice(ext%o%crystal))
        call diary(ext%o%layout)
100     call glean(thy(ext))
        if (error("Exit external_mod::diary_ext")) continue
      end subroutine

      subroutine diary_crystal_step(ext,step)
!doc$ subroutine diary(ext)
        type(external_obj) :: ext
        integer, intent(in) :: step
!       modifies : diary stream
!       effects: Prints atom-related data to diaryf for step in an atom relaxation or MD sequence.
!       errors: Passes errors.

!cod$
        call my(ext)
        call diary(ext%o%crystal,step)
        call diary_atom_neighbors_i(ext%o)
100     call glean(thy(ext))
        if (error("Exit external_mod::diary_crystal_step")) continue
      end subroutine

! private routines

      subroutine own_i(ext)
        type(external_obj) :: ext
        type(external_obj) :: ext_t
        if (ext%ref < ext%o%ref) then
          allocate( ext_t%o )
          ext_t%o%ref = 0
          ext_t%o%init_method = ext%o%init_method
          call my(ext%o%crystal,ext_t%o%crystal)
          call my(ext%o%layout,ext_t%o%layout)
          call my(ext%o%xc_type,ext_t%o%xc_type)
          call my(ext%o%ao,ext_t%o%ao)
          call my(ext%o%lattice_group,ext_t%o%lattice_group)
          call my(ext%o%space_group,ext_t%o%space_group)
          call my(ext%o%double_group,ext_t%o%double_group)
          ext_t%o%g = ext%o%g
          ext%o%ref = ext%o%ref - ext%ref
          ext%o => ext_t%o
          ext%o%ref = ext%o%ref + ext%ref
        end if
      end subroutine

      subroutine diary_construction_i(extr,sym_atoms)
        type(external_rep) :: extr
        logical, intent(in) :: sym_atoms
!       effects: Writes extr information to the diary.

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"External object construction:")')
        call diary_representation(extr%ao)
        if (i_access( diaryfile() )) then
          select case (extr%init_method)
          case (RESTART_FILE)
            write(x_unit(diaryfile()),'(/,t4,"Restart-file initialization")')
          end select
        end if
        call diary(extr%crystal)
        call diary_matching_radii_i(extr)
        call diary_atom_neighbors_i(extr)
        call diary(extr%lattice_group,lat=x_lattice(extr%crystal))
        call diary(extr%space_group,lat=x_lattice(extr%crystal))
        if (sym_atoms) then
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Symmetrizing the atom coordinates")')
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(t4,"Updating the space group")')
        end if
        call diary(extr%layout)
        call diary(extr%xc_type)
        if (error("Exit external_mod::diary_construction_i")) continue
      end subroutine

      subroutine diary_matching_radii_i(extr)
        type(external_rep) :: extr
!       effects: Prints atom matching radii to the diary.

        integer :: it

        if (x_type(extr%ao) == PAW) then
          if (i_access( diaryfile() )) then
            write(x_unit(diaryfile()),'(/,t4,"Atom matching radii:",/)')
            do it = 1,x_n_types(extr%ao)
              write(x_unit(diaryfile()),'(t7,a6,":",2x,f6.3)') trim(x_type_name(extr%ao,it)), x_type_matching_radius(extr%ao,it)
            end do
          end if
        end if
        if (error("Exit external_mod::diary_matching_radii_i")) continue
      end subroutine

      subroutine diary_atom_neighbors_i(extr)
        type(external_rep) :: extr
!       effects: Prints neighbor list to the diary and gives a warning if atoms overlap.

        logical, dimension(:), allocatable :: atom_overlap
        integer :: ia, ia2, in, na, nn
        real(double) :: mr, mr2

        na = x_n_atoms(x_atoms(extr%crystal))
        allocate( atom_overlap(na) )
        atom_overlap = .false.
        do ia = 1,na
          mr = x_atom_matching_radius(extr%ao,ia)
          do in = 1,x_n_neighbors(extr%crystal,ia)
            ia2 = x_neighbor_index(extr%crystal,ia,in)
            mr2 = x_atom_matching_radius(extr%ao,ia2)
            if (x_neighbor_distance(extr%crystal,ia,in) < (mr + mr2)) atom_overlap(ia) = .true.
            if (atom_overlap(ia)) exit
          end do
        end do
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Nearest-neighbor atoms (range = ",f0.3,"):")') x_neighbor_range(extr%crystal)
          write(x_unit(diaryfile()),'(/,t7,"atom",5x,"neighbor:separation",/)')
          do ia = 1,na
            write(x_unit(diaryfile()),'(t7,i4,2x)',advance="no") ia
            nn = x_n_neighbors(extr%crystal,ia)
            do in = 1,nn
              write(x_unit(diaryfile()),'(2x,i4,":",f5.2)',advance="no") x_neighbor_index(extr%crystal,ia,in), &
                                                                         x_neighbor_distance(extr%crystal,ia,in)
              if ((nn > in) .and. (mod(in,8) == 0)) write(x_unit(diaryfile()),'(/,t13)',advance="no")
            end do
            if (atom_overlap(ia)) then
              write(x_unit(diaryfile()),'(2x,"Warning: atom overlap")')
            else
              write(x_unit(diaryfile()),*)
            end if
          end do
        end if
        if (allocated( atom_overlap )) deallocate( atom_overlap )
        if (error("Exit external_mod::diary_atom_neighbors_i")) continue
      end subroutine

      end module

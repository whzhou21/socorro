! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module atoms_mod
!doc$ module atoms_mod

!     One datatype is available here: type(atoms_obj).

!     atoms_mod contains information about the atoms in the simulation cell.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use ghost_mod
      use math_mod
      use lattice_mod
      use tagio_mod

!cod$
      implicit none
      private

      type :: atom
        character(tag_sz) :: tag                 ! string indicating the associated atomic_operator_mod data
        real(double), dimension(3) :: position   ! position in the lattice representation
        real(double), dimension(3) :: velocity   ! velocity
      end type

      type :: atoms_rep
        integer :: ref
        type(ghost) :: g
        integer :: natoms                              ! number of atoms
        type(atom), dimension(:), pointer :: atom_list ! atom data: note that the size may be > natoms
      end type

      type, public :: atoms_obj
        private
        integer :: ref
        type(atoms_rep), pointer :: o
      end type

!doc$
      public :: atoms
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_type
      public :: x_n_atoms
      public :: x_position
      public :: add
      public :: remove
      public :: move
      public :: save
      public :: diary
      public :: write_restart

!cod$
      interface atoms
        module procedure constructor_at
      end interface
      interface my
        module procedure my_at, my_new_at
      end interface
      interface thy
        module procedure thy_at
      end interface
      interface glean
        module procedure glean_at
      end interface
      interface bequeath
        module procedure bequeath_at
      end interface
      interface assignment(=)
        module procedure assign_at
      end interface
      interface x_ref
        module procedure at_ref
      end interface
      interface x_ghost
        module procedure at_ghost
      end interface
      interface x_type
        module procedure at_type, at_type_i
      end interface
      interface x_n_atoms
        module procedure at_n_atoms
      end interface
      interface x_position
        module procedure at_position, at_position_i
      end interface
      interface add
        module procedure add_at
      end interface
      interface remove
        module procedure remove_at
      end interface
      interface move
        module procedure move_at, move_at_i
      end interface
      interface save
        module procedure save_at
      end interface
      interface diary
        module procedure diary_at
      end interface
      interface write_restart
        module procedure write_restart_at
      end interface

      contains

! public routines

      function constructor_at(restf) result(at)
!doc$ function atoms(restf) result(at)
        type(tagio_obj), optional :: restf
        type(atoms_obj) :: at
!       requires: restf be positioned inside the CRYSTAL block.
!       effects: Creates a new at.
!       errors: restart-file tags not found.

!cod$
        character(1) :: tios
        integer :: ia
        integer(long) :: dsize, iosl, ndata, s4

        if (present(restf)) call my(restf)

        at%ref = 0
        allocate( at%o )
        at%o%ref = 0
        at%o%g = x_ghost()

        ! open the ATOMS block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMS block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)
        end if

        if (present(restf)) then

          ! read the number of atoms
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_ATOMS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_ATOMS tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long ; ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            at%o%natoms = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,at%o%natoms)

          ! read the atom information
          if (i_access(restf)) tios = findfirsttag(restf,"ATOM_INFORMATION")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: ATOM_INFORMATION tag was not found")) goto 100
          allocate( at%o%atom_list(at%o%natoms) )
          if (i_access(restf)) then
            do ia = 1,at%o%natoms
              dsize = 1 ; ndata = tag_sz
              call readf(at%o%atom_list(ia)%tag,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              dsize = sizeof_double ; ndata = 3
              call readf(at%o%atom_list(ia)%position,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              call readf(at%o%atom_list(ia)%velocity,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end do
          end if
          if (i_comm(restf)) then
            do ia = 1,at%o%natoms
              call broadcast(FILE_SCOPE,at%o%atom_list(ia)%tag)
              call broadcast(FILE_SCOPE,at%o%atom_list(ia)%position)
              call broadcast(FILE_SCOPE,at%o%atom_list(ia)%velocity)
            end do
          end if

        else

          ! create an empty data structure
          at%o%natoms = 0
          allocate( at%o%atom_list(1) )

        end if

        ! close the ATOMS block
100     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

200     if (present(restf)) call glean(thy(restf))

        if (error("Exit: constructor_at")) continue

      end function

      subroutine my_at(at)
!doc$ subroutine my(at)
        type(atoms_obj) :: at

!cod$
        at%ref = at%ref + 1
        at%o%ref = at%o%ref + 1
      end subroutine

      subroutine my_new_at(ati,at)
!doc$ subroutine my(ati,at)
      type(atoms_obj) :: ati, at

!cod$
        at%ref = 1
        at%o => ati%o
        at%o%ref = at%o%ref + 1
      end subroutine

      function thy_at(at) result(ato)
!doc$ function thy(at) result(ato)
        type(atoms_obj) :: at, ato

!cod$
        at%ref = at%ref - 1
        at%o%ref = at%o%ref - 1
        ato%ref = at%ref
        ato%o => at%o
      end function

      subroutine glean_at(at)
!doc$ subroutine glean(at)
        type(atoms_obj) :: at

!cod$
        if (at%o%ref < 1) then
          deallocate( at%o%atom_list )
          deallocate( at%o )
        end if
      end subroutine

      subroutine bequeath_at(at)
!doc$ subroutine bequeath(at)
        type(atoms_obj) :: at

!cod$
        continue
      end subroutine

      subroutine assign_at(at,at2)
!doc$ subroutine assignment(=)(at,at2)
        type(atoms_obj), intent(inout) :: at
        type(atoms_obj), intent(in) :: at2

!cod$
        type(atoms_obj) :: att
        call my(at2)
        att%o => at%o
        at%o%ref = at%o%ref - at%ref
        at%o => at2%o
        at%o%ref = at%o%ref + at%ref
        call glean(att)
        call glean(thy(at2))
      end subroutine

      function at_ref(at) result(r)
!doc$ function x_ref(at) result(r)
        type(atoms_obj) :: at
        integer, dimension(2) :: r
!       effects: Returns at%ref and at%o%ref.

!cod$
        r(1) = at%ref
        r(2) = at%o%ref
        call glean(at)
      end function

      function at_ghost(at) result(g)
!doc$ function x_ghost(at) result(g)
        type(atoms_obj) :: at
        type(ghost) :: g
!       effects: Returns the ghost of at.

!cod$
        call my(at)
        g = at%o%g
        call glean(thy(at))
      end function

      function at_type(at) result(atg)
!doc$ function x_type(at) result(atg)
        type(atoms_obj) :: at
        character(tag_sz), dimension(at%o%natoms) :: atg
!       effects: Returns strings identifying the atom types.

!cod$
        integer :: ia
        call my(at)
        do ia = 1,at%o%natoms
          atg(ia) = at%o%atom_list(ia)%tag
        end do
        call glean(thy(at))
      end function

      function at_type_i(at,ia) result(atg)
!doc$ function x_type(at,ia) result(atg)
        type(atoms_obj) :: at
        integer, intent(in) :: ia
        character(tag_sz) :: atg
!       effects: Returns a string identifying the type of atom ia.
!       errors: ia out of range.

!cod$
        call my(at)
        if (error((ia < 1) .or. (ia > at%o%natoms),"ERROR: ia is out of range")) goto 100
        atg = at%o%atom_list(ia)%tag
100     call glean(thy(at))
        if (error("Exit atoms_mod::at_type_i")) continue
      end function

      function at_n_atoms(at) result(na)
!doc$ function x_n_atoms(at) result(na)
        type(atoms_obj) :: at
        integer :: na
!       effects: Returns the number of atoms in at.

!cod$
        call my(at)
        na = at%o%natoms
        call glean(thy(at))
      end function

      function at_position(at) result(pos)
!doc$ function x_position(at) result(pos)
        type(atoms_obj) :: at
        real(double), dimension(3,at%o%natoms) :: pos
!       effects: Returns the atom positions in the lattice representation.

!cod$
        integer :: ia
        call my(at)
        do ia = 1,at%o%natoms
          pos(:,ia) = at%o%atom_list(ia)%position
        end do
        call glean(thy(at))
      end function

      function at_position_i(at,ia) result(pos)
!doc$ function x_position(at,ai) result(pos)
        type(atoms_obj) :: at
        real(double), dimension(3) :: pos
        integer, intent(in) :: ia
!       effects: Returns the position of atom ia in the lattice representation.
!       errors: ia out of range

!cod$
        call my(at)
        pos = real((/0,0,0/),double)
        if (error((ia < 1) .or. (ia > at%o%natoms),"ERROR: ia is out of range")) goto 100
        pos = at%o%atom_list(ia)%position
100     call glean(thy(at))
        if (error("Exit atoms_mod::at_position_i")) continue
      end function

      subroutine add_at(at,tag,pos,velocity)
!doc$ subroutine add(at,tag,pos,velocity)
        type(atoms_obj) :: at
        character(tag_sz), intent(in) :: tag
        real(double), dimension(3), intent(in) :: pos
        real(double), dimension(3), intent(in), optional :: velocity
!       modifies: at
!       requires: pos be in the lattice representation.
!       effects: Adds new atom at position pos. Constrains the atom position to lie in the range [0,1).

!cod$
        type(atom), dimension(:), pointer :: tmp
        call my(at)
        call own_i(at)
        at%o%g = x_ghost()
        if (at%o%natoms == size(at%o%atom_list)) then
          allocate( tmp(at%o%natoms*2+1) )
          tmp(1:at%o%natoms) = at%o%atom_list(1:at%o%natoms)
          deallocate( at%o%atom_list )
          at%o%atom_list => tmp
        end if
        at%o%natoms = at%o%natoms + 1
        at%o%atom_list(at%o%natoms)%tag = tag
        at%o%atom_list(at%o%natoms)%position = pos
        call centralize_position_i(at%o%atom_list(at%o%natoms)%position)
        if (any(at%o%atom_list(at%o%natoms)%position /= pos)) call warn("WARNING: an atom was translated into the central cell")
        if (present(velocity)) then
          at%o%atom_list(at%o%natoms)%velocity = velocity
        else
          at%o%atom_list(at%o%natoms)%velocity = real((/0,0,0/),double)
        end if
        call glean(thy(at))
      end subroutine

      subroutine remove_at(at,ia)
!doc$ subroutine remove(at,ia)
        type(atoms_obj) :: at
        integer, intent(in) :: ia
!       modifies: at
!       effects: Removes atom ia.
!       errors: ia out of range.

!cod$
        call my(at)
        call own_i(at)
        at%o%g = x_ghost()
        if (error((ia < 1) .or. (ia > at%o%natoms),"ERROR: ia is out of range")) goto 100
        at%o%atom_list(ia:at%o%natoms-1) = at%o%atom_list(ia+1:at%o%natoms)
        at%o%natoms = at%o%natoms - 1
100     call glean(thy(at))
        if (error("Exit atoms_mod::remove_at")) continue
      end subroutine

      subroutine move_at(at,pos)
!doc$ subroutine move(at,pos)
        type(atoms_obj) :: at
        real(double), dimension(1:,1:), intent(in) :: pos
!       modifies: at
!       requires: pos be in the lattice representation.
!       effects: Changes the positions of the atoms to those in pos.
!       errors: If shape(pos)/=(/3,x_n_atoms(atoms)/).

!cod$
        integer :: ia
        call my(at)
        call own_i(at)
        at%o%g = x_ghost()
        if (error((size(pos,1) /= 3) .or. (size(pos,2) /= at%o%natoms),"ERROR: incorrect dimensions for pos")) goto 100
        do ia = 1,at%o%natoms
          at%o%atom_list(ia)%position = pos(:,ia)
          call centralize_position_i(at%o%atom_list(ia)%position)
        end do
100     call glean(thy(at))
        if (error("Exit atoms_mod::move_at")) continue
      end subroutine

      subroutine move_at_i(at,ia,pos)
!doc$ subroutine move(at,ia,pos)
        type(atoms_obj) :: at
        integer, intent(in) :: ia
        real(double), dimension(3), intent(in) :: pos
!       modifies: at
!       requires: pos be in the lattice representation.
!       effects: Changes the position of atom ia to pos.
!       errors: ia out of range.

!cod$
        call my(at)
        call own_i(at)
        at%o%g = x_ghost()
        if (error((ia < 1) .or. (ia > at%o%natoms),"ERROR: ia is out of range")) goto 100
        at%o%atom_list(ia)%position = pos
100     call glean(thy(at))
        if (error("Exit atoms_mod::move_at_i")) continue
      end subroutine

      subroutine save_at(at,prefix)
!doc$ subroutine save(at,prefix)
        type(atoms_obj) :: at
        character(*), intent(in) :: prefix
!       effects: Writes atoms information to the file prefix//"atoms".
!       errors: If problem opening the file.

!cod$
        integer :: ia, ios
        type(file_obj) :: f
        call my(at)
        call my(file(prefix//"atoms"),f)
        if (i_access(f)) open(x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100
        if (i_access(f)) then
          write(x_unit(f),*) at%o%natoms
          do ia = 1,at%o%natoms
            write(x_unit(f),*) at%o%atom_list(ia)%tag, at%o%atom_list(ia)%position
          end do
        end if
        if (i_access(f)) close(x_unit(f))
100     call glean(thy(f))
        call glean(thy(at))
999     if (error("Exit atoms_mod::save_at")) continue
      end subroutine

      subroutine diary_at(at,lat)
!doc$ subroutine diary(at)
        type(atoms_obj) :: at
        type(lattice_obj) :: lat
!       effects: Writes atom information to diary.

!cod$
        integer :: ia
        real(double), dimension(3) :: pos_cart, pos_lat
        call my(at)
        call my(lat)
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Atom Positions:")')
          write(x_unit(diaryfile()),'(/,t8,"atom",4x," type ",11x,"x",10x,"y",10x,"z",15x,"a1",9x,"a2",9x,"a3")')
          write(x_unit(diaryfile()),'(t7,"----------------    ---------------------------------     -----------", &
                                       &  "----------------------")')
          do ia = 1,at%o%natoms
            pos_lat = at%o%atom_list(ia)%position
            pos_cart = lat2r(lat,pos_lat)
            write(x_unit(diaryfile()),'(t7,i4,6x,a6,2(5x,3(f9.5,2x)))') ia, at%o%atom_list(ia)%tag, pos_cart, pos_lat
          end do
        end if
        call glean(thy(at))
        call glean(thy(lat))
      end subroutine

      subroutine write_restart_at(at,nrestf)
!doc$ subroutine write_restart(at,nrestf)
        type(atoms_obj) :: at
        type(tagio_obj) :: nrestf
!       effects: Writes at information to nrestf.
!       errors: Passes errors.

!cod$
        integer :: ia
        integer(long) :: dsize, iosl, ndata, s4

        call my(at)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the ATOMS block
          call startblock(nrestf,"ATOMS")

          ! write the number of atoms
          call writetag(nrestf,"NUMBER_OF_ATOMS")
          s4 = at%o%natoms ; dsize = sizeof_long ; ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! write the atom tags, positions, and velocities
          call writetag(nrestf,"ATOM_INFORMATION")
          do ia = 1,at%o%natoms
            dsize = 1 ; ndata = tag_sz
            call writef(at%o%atom_list(ia)%tag,dsize,ndata,x_tagfd(nrestf),iosl)
            dsize = sizeof_double ; ndata = 3
            call writef(at%o%atom_list(ia)%position,dsize,ndata,x_tagfd(nrestf),iosl)
            call writef(at%o%atom_list(ia)%velocity,dsize,ndata,x_tagfd(nrestf),iosl)
          end do

          ! end the ATOMS block
          call endblock(nrestf)

        end if

100     call glean(thy(at))
        call glean(thy(nrestf))

        if (error("Exit atoms_mod::write_restart_at")) continue

      end subroutine

! private routines

      subroutine own_i(at)
        type(atoms_obj) :: at
        type(atoms_obj) :: att
        integer :: ia
        if (at%ref < at%o%ref) then
          allocate( att%o )
          att%o%ref = 0
          att%o%natoms = at%o%natoms
          allocate( att%o%atom_list(at%o%natoms) )
          do ia = 1,at%o%natoms
            att%o%atom_list(ia)%tag = at%o%atom_list(ia)%tag
            att%o%atom_list(ia)%position = at%o%atom_list(ia)%position
            att%o%atom_list(ia)%velocity = at%o%atom_list(ia)%velocity
          end do
          att%o%g = at%o%g
          at%o%ref = at%o%ref - at%ref
          at%o => att%o
          at%o%ref = at%o%ref + at%ref
        end if
      end subroutine

      subroutine centralize_position_i(pos)
        real(double), dimension(:) :: pos

        logical :: inside
        integer :: ic
        real(double), parameter :: tol_zero = 1.0e-9_double
        inside = .false.
        do while (.not.inside)
          do ic = 1,3
            if ((pos(ic) > 0.0_double) .or. (pos(ic) .in. nbhd(0.0_double,tol_zero))) then
              inside = .true.
            else
              pos(ic) = pos(ic) + 1.0_double
              inside = .false.
              exit
            end if
            if ((pos(ic) < 1.0_double) .and. (pos(ic) .out. nbhd(1.0_double,tol_zero))) then
              inside = .true.
            else
              pos(ic) = pos(ic) - 1.0_double
              inside = .false.
              exit
            end if
          end do
        end do

      end subroutine

      end module

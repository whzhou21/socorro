!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module crystal_mod
!doc$ module crystal_mod

!     One datatype is available here: type(crystal_obj).

!     Crystal_mod contains components needed to describe the crystal
!     structure including the lattice, atom positions, and atom types.

      use kind_mod
      use arg_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use ghost_mod
      use math_mod
      use lattice_mod
      use atoms_mod
      use tagio_mod
      use timing_mod

!cod$
      implicit none
      private

      real(double), parameter :: tol_equal = 1.0e-8_double
      real(double), parameter :: neighbor_range_default = 5.0_double
      character(len("cartesian")), parameter :: CARTESIAN_STR = "cartesian"
      character(len("lattice")), parameter :: LATTICE_STR = "lattice"

      type :: crystal_rep
        integer :: ref
        type(ghost) :: g
        character(line_len) :: name
        type(atoms_obj) :: atoms
        type(lattice_obj) :: lattice
        real(double) :: neighbor_range
        type(ghost) :: neighbor_list_ghost
        integer, dimension(:), pointer :: n_neighbors
        integer, dimension(:,:), pointer :: neighbor_index
        real(double), dimension(:,:), pointer :: neighbor_distance
        real(double), dimension(:,:,:), pointer :: neighbor_vector
      end type

      type, public :: crystal_obj
        private
        integer :: ref
        type(crystal_rep), pointer :: o
      end type

!doc$
      public :: crystal
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_name
      public :: x_atoms
      public :: x_lattice
      public :: x_n_neighbors
      public :: x_neighbor_index
      public :: x_neighbor_range
      public :: x_neighbor_distance
      public :: x_neighbor_vector
      public :: save
      public :: diary
      public :: write_restart

!cod$
      interface crystal
        module procedure constructor_cr_1, constructor_cr_2, constructor_cr_3
      end interface
      interface update
        module procedure update_cr
      end interface
      interface my
        module procedure my_cr, my_new_cr
      end interface
      interface thy
        module procedure thy_cr
      end interface
      interface glean
        module procedure glean_cr
      end interface
      interface bequeath
        module procedure bequeath_cr
      end interface
      interface assignment(=)
        module procedure assign_cr
      end interface
      interface x_ref
        module procedure cr_ref
      end interface
      interface x_ghost
        module procedure cr_ghost
      end interface
      interface x_name
        module procedure cr_name
      end interface
      interface x_atoms
        module procedure cr_atoms
      end interface
      interface x_lattice
        module procedure cr_lattice
      end interface
      interface x_n_neighbors
        module procedure cr_n_neighbors
      end interface
      interface x_neighbor_range
        module procedure cr_neighbor_range
      end interface
      interface x_neighbor_index
        module procedure cr_neighbor_index
      end interface
      interface x_neighbor_distance
        module procedure cr_neighbor_distance
      end interface
      interface x_neighbor_vector
        module procedure cr_neighbor_vector
      end interface
      interface save
        module procedure save_crys
      end interface
      interface diary
        module procedure diary_cr, diary_crystal_step
      end interface
      interface write_restart
        module procedure write_restart_cr
      end interface

      contains

      function constructor_cr_1(path) result(cr)
!doc$ function crystal(path) result(cr)
        character(*), optional :: path
        type(crystal_obj) :: cr
!       effects: Constructs cr from file trim(crystal_path) or trim(path) if present.
!       errors: Trouble opening the file. File data missing or incorrect. Lattice vectors linearly dependent. Passes errors.

!cod$
        logical :: cartesian_mode, exist_file, linear_dependence
        character(line_len) :: tag_ls
        character(tag_sz) :: tag_ts
        integer :: ia, ios, na
        real(double) :: latc
        real(double), dimension(3) :: pos, v1, v2, v3
        real(double), dimension(3,3) :: mat
        type(file_obj) :: f
        logical :: found

        if (error("  Error on entry")) then
          cr%ref = 0
          allocate( cr%o )
          cr%o%ref = 0
          goto 999
        end if

        cr%ref = 0
        allocate( cr%o )
        cr%o%ref = 0
        cr%o%g = x_ghost()

        if (present(path)) then
          call my(file(trim(path)),f)
        else
          call my(file(trim(crystal_path)),f)
        end if
        if (i_access(f)) inquire(file=x_name(f),exist=exist_file)
        if (i_comm(f)) call broadcast(FILE_SCOPE,exist_file)
        if (error(.not.exist_file,"ERROR: file does not exist")) goto 200

        if (i_access(f)) open(x_unit(f),file=x_name(f),status='old',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 200

        if (i_access(f)) read(x_unit(f),'(a)',iostat=ios) tag_ls
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read name")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,tag_ls)
        cr%o%name = tag_ls(1:len_trim(tag_ls))

        if (i_access(f)) read(x_unit(f),*,iostat=ios) latc
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read lattice constant")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,latc)
        if (i_access(f)) read(x_unit(f),*,iostat=ios) v1
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read lattice vector")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,v1)
        if (i_access(f)) read(x_unit(f),*,iostat=ios) v2
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read lattice vector")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,v2)
        if (i_access(f)) read(x_unit(f),*,iostat=ios) v3
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read lattice vector")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,v3)
        mat(:,1) = v1
        mat(:,2) = v2
        mat(:,3) = v3
        linear_dependence = ( determinant(mat) .in. nbhd(0.0_double,tol_equal) )
        if (error(linear_dependence,"ERROR: lattice vectors are linearly dependent")) goto 100
        call my(lattice(v1,v2,v3,latc),cr%o%lattice)

        if (i_access(f)) read(x_unit(f),'(a)',iostat=ios) tag_ls
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read designator")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,tag_ls)
        select case (tag_ls(1:len_trim(tag_ls)))
        case ( "cartesian", "Cartesian", "CARTESIAN", "cart", "Cart", "CART" )
          cartesian_mode = .true.
        case ( "lattice", "Lattice", "LATTICE", "lat", "Lat", "LAT" )
          cartesian_mode = .false.
        case default
          if (error(.true.,"ERROR: unrecognized designator")) goto 100
        end select

        if (i_access(f)) read(x_unit(f),*,iostat=ios) na
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read the number of atoms")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,na)
        if (error(na <= 0,"ERROR: incorrect number of atoms")) goto 100
        call my(atoms(),cr%o%atoms)
        do ia = 1,na
          if (i_access(f)) read(x_unit(f),*,iostat=ios) tag_ts, pos
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to read atom information")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,tag_ts)
          if (i_comm(f)) call broadcast(FILE_SCOPE,pos)
          if (cartesian_mode) pos = r2lat(cr%o%lattice,x_lattice_constant(cr%o%lattice)*pos)
          call add(cr%o%atoms,tag_ts,pos)
        end do

        call arg("neighbor_range",cr%o%neighbor_range,found)
        if (.not.found) cr%o%neighbor_range = neighbor_range_default
        if (error(cr%o%neighbor_range <= 0,"ERROR: neighbor_range <= 0")) goto 100

        cr%o%neighbor_list_ghost = x_ghost()
        allocate( cr%o%n_neighbors(na) )
        allocate( cr%o%neighbor_index(1,na) )
        allocate( cr%o%neighbor_distance(1,na) )
        allocate( cr%o%neighbor_vector(1,1,na) )

100     if (i_access(f)) close(x_unit(f))
200     call glean(thy(f))

999     if (error("Exit crystal_mod::constructor_cr_1")) continue

      end function
      
      function constructor_cr_2(restf) result(cr)
!doc$ function crystal(restf) result(cr)
        type(tagio_obj) :: restf
        type(crystal_obj) :: cr
!       requires: restf be positioned inside the EXTERNAL block.
!       effects: Constructs a new cr from information in restf.
!       errors: CRYSTAL block or TAGS not found. Passes errors.

!cod$
        logical :: found
        character(1) :: tios
        integer :: na
        integer(long) :: dsize, iosl, ndata

        call my(restf)

        cr%ref = 0
        allocate( cr%o )
        cr%o%ref = 0
        cr%o%g = x_ghost()

        ! open the CRYSTAL block
        if (i_access(restf)) tios = findfirsttag(restf,"CRYSTAL")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: CRYSTAL block was not found")) goto 200
        if (i_access(restf)) call openblock(restf)

        ! determine the name
        if (i_access(restf)) tios = findfirsttag(restf,"NAME")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: NAME tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = 1 ; ndata = line_len
          call readf(cr%o%name,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,cr%o%name)

        ! form the lattice
        call my(lattice(restf),cr%o%lattice) ; if (error()) goto 100

        ! form the atoms
        call my(atoms(restf),cr%o%atoms)

        ! close the CRYSTAL block
100     if (i_access(restf)) call closeblock(restf)
        if (error()) goto 200

        ! determine the neighbor range
        call arg("neighbor_range",cr%o%neighbor_range,found)
        if (.not.found) cr%o%neighbor_range = neighbor_range_default
        if (error(cr%o%neighbor_range <= 0,"ERROR: neighbor_range <= 0")) goto 200

        ! allocate arrays for the neighbor information
        cr%o%neighbor_list_ghost = x_ghost()
        na = x_n_atoms(cr%o%atoms)
        allocate( cr%o%n_neighbors(na) )
        allocate( cr%o%neighbor_index(1,na) )
        allocate( cr%o%neighbor_distance(1,na) )
        allocate( cr%o%neighbor_vector(1,1,na) )

200     call glean(thy(restf))

        if (error("Exit crystal_mod::constructor_cr_2")) continue

      end function
      
      function constructor_cr_3(at,lat,name) result(cr)
!doc$ function crystal(at,lat,name) result(cr)
        type(atoms_obj) :: at
        type(lattice_obj) :: lat
        character(*), intent(in), optional :: name
        type(crystal_obj) :: cr
!       effects: Constructs a new cr from at, lat, and (optionally) name.
!       errors: Passes errors.

!cod$
        logical :: found
        integer :: na

        if (error("  Error on entry")) then
          cr%ref = 0
          allocate( cr%o )
          cr%o%ref = 0
          goto 999
        end if

        call my(at)
        call my(lat)

        cr%ref = 0
        allocate( cr%o )
        cr%o%ref = 0

        cr%o%g = x_ghost()

        if (present(name)) then
          cr%o%name = name
        else
          cr%o%name = ""
        end if

        call my(at,cr%o%atoms)
        call my(lat,cr%o%lattice)

        call arg("neighbor_range",cr%o%neighbor_range,found)
        if (.not.found) cr%o%neighbor_range = neighbor_range_default
        if (error(cr%o%neighbor_range <= 0,"ERROR: neighbor_range <= 0")) goto 100

        cr%o%neighbor_list_ghost = x_ghost()
        na = x_n_atoms(cr%o%atoms)
        allocate( cr%o%n_neighbors(na) )
        allocate( cr%o%neighbor_index(1,na) )
        allocate( cr%o%neighbor_distance(1,na) )
        allocate( cr%o%neighbor_vector(1,1,na) )

100     call glean(thy(at))
        call glean(thy(lat))

999     if (error("Exit crystal_mod::constructor_cr_3")) continue

      end function
      
      subroutine update_cr(cr,atoms,lattice,neighbor_range)
!doc$ subroutine update(cr,atoms,lattice,neighbor_range)
        type(crystal_obj) :: cr
        type(atoms_obj), optional :: atoms
        type(lattice_obj), optional :: lattice
        real(double), intent(in), optional :: neighbor_range
!       modifies: cr
!       effects: Updates cr with respect to its dependencies:
!                 - atoms_change: atom positions
!                 - lattice_change: lattice vectors (not yet implemented)
!                 - range_change: range in neighbor list
!       errors: For changes to cr%o%lattice.

!cod$ 
        logical :: any_change, atoms_change, lattice_change, range_change

        call my(cr)

        if (present(atoms)) then
          call my(atoms)
          atoms_change = ( x_ghost(cr%o%atoms) /= x_ghost(atoms) )
        else
          atoms_change = .false.
        end if
        if (present(lattice)) then
          call my(lattice)
          lattice_change = ( x_ghost(cr%o%lattice) /= x_ghost(lattice) )
        else
          lattice_change = .false.
        end if
        if (present(neighbor_range)) then
          range_change = ( neighbor_range /= cr%o%neighbor_range )
        else
          range_change = .false.
        end if
        any_change = ( atoms_change .or. lattice_change .or. range_change )
        if (any_change) then
          call own_i(cr)
          cr%o%g = x_ghost()
          if (atoms_change) cr%o%atoms = atoms
          if (lattice_change) then
            if (error(.true.,"ERROR: lattice changes are not currently allowed")) goto 100
          end if
          if (range_change) then
             cr%o%neighbor_range = neighbor_range
             if (error(cr%o%neighbor_range <= 0,"ERROR: neighbor_range <= 0")) goto 100
!            Do not calculate the actual neighbor list until it is needed
             cr%o%neighbor_list_ghost = x_ghost()  ! Make sure neighbor_list_ghost will not match cr%o%g
          end if
          call save(cr)
        end if

100     call glean(thy(cr))
        if (present(atoms)) call glean(thy(atoms))
        if (present(lattice)) call glean(thy(lattice))

        if (error("Exit crystal_mod::update_cr")) continue

      end subroutine

      subroutine my_cr(cr)
!doc$ subroutine my(cr)
        type(crystal_obj) :: cr

!cod$
        cr%ref = cr%ref + 1
        cr%o%ref = cr%o%ref + 1
      end subroutine

      subroutine my_new_cr(cri,cr)
!doc$ subroutine my(cri,cr)
        type(crystal_obj) :: cri, cr

!cod$
        cr%ref = 1
        cr%o => cri%o
        cr%o%ref = cr%o%ref + 1
      end subroutine

      function thy_cr(cr) result(cro)
!doc$ function thy(cr) result(cro)
        type(crystal_obj) :: cr, cro

!cod$
        cr%ref = cr%ref - 1
        cr%o%ref = cr%o%ref - 1
        cro%ref = cr%ref
        cro%o => cr%o
      end function
      
      subroutine glean_cr(cr)
!doc$ subroutine glean(cr)
        type(crystal_obj) :: cr

!cod$
        if (cr%o%ref < 1) then
          call glean(thy(cr%o%lattice))
          call glean(thy(cr%o%atoms))
          deallocate ( cr%o%n_neighbors)
          deallocate ( cr%o%neighbor_index)
          deallocate ( cr%o%neighbor_distance )
          deallocate ( cr%o%neighbor_vector )
          deallocate( cr%o )
        end if
      end subroutine

      subroutine bequeath_cr(cr)
!doc$ subroutine bequeath(cr)
        type(crystal_obj) :: cr

!cod$
        continue
      end subroutine

      subroutine assign_cr(cr,cr2)
!doc$ subroutine assignment(=)(cr,cr2)
        type(crystal_obj), intent(inout) :: cr
        type(crystal_obj), intent(in) :: cr2

!cod$
        type(crystal_obj) :: crt
        call my(cr2)
        crt%o => cr%o
        cr%o%ref = cr%o%ref - cr%ref
        cr%o => cr2%o
        cr%o%ref = cr%o%ref + cr%ref
        call glean(crt)
        call glean(thy(cr2))
      end subroutine

      function cr_ref(cr) result(r)
!doc$ function x_ref(cr) result(r)
        type(crystal_obj) :: cr
        integer, dimension(2) :: r
!       effects: Returns cr%ref and cr%o%ref.

!cod$
        r(1) = cr%ref
        r(2) = cr%o%ref
        call glean(cr)
      end function
      
      function cr_ghost(cr) result(g)
!doc$ function x_ghost(cr) result(g)
        type(crystal_obj) :: cr
        type(ghost) :: g
!       effects: Returns the ghost of cr.

!cod$
        call my(cr)
        g = cr%o%g
        call glean(thy(cr))
      end function

      function cr_lattice(cr) result(lat)
!doc$ function x_lattice(cr) result(lat)
        type(crystal_obj) :: cr
        type(lattice_obj) :: lat
!       effects: Returns the lattice in cr.

!cod$
        call my(cr)
        call my(cr%o%lattice,lat)
        call bequeath(thy(lat))
        call glean(thy(cr))
      end function
      
      function cr_name(cr) result(n)
!doc$ function x_name(cr) result(n)
        type(crystal_obj) :: cr
        character(line_len) :: n
!       effects: Returns the name of cr.

!cod$
        call my(cr)
        n = cr%o%name
        call glean(thy(cr))
      end function

      function cr_atoms(cr) result(at)
!doc$ function x_atoms(cr) result(at)
        type(crystal_obj) :: cr
        type(atoms_obj) :: at
!       effects: Returns the atoms_obj of cr.

!cod$
        call my(cr)
        call my(cr%o%atoms,at)
        call bequeath(thy(at))
        call glean(thy(cr))
      end function

      function cr_n_neighbors(cr,ia) result(nn)
!doc$ function x_n_neighbors(cr,ia) result(nn)
        type(crystal_obj) :: cr
        integer, intent(in) :: ia
        integer :: nn
!       effects: Returns the number of neighbors of atom ia.
!       errors: ia out of range.

!cod$
        call my(cr)
        if (error((ia < 1) .or. (ia > x_n_atoms(cr%o%atoms)),"ERROR: atom index out of range")) goto 100
        if (cr%o%neighbor_list_ghost /= cr%o%g) call find_neighbors_i(cr)
        nn = cr%o%n_neighbors(ia)
100     call glean(thy(cr))
        if (error("Exit crystal_mod::cr_n_neighbors")) continue

      end function

      function cr_neighbor_range(cr) result(nc)
!doc$ function x_neighbor_range(cr) result(nc)
        type(crystal_obj) :: cr
        real(double) :: nc
!       effects: Returns the range used to compute the neighbor list.
!       errors: none

!cod$
        call my(cr)
        nc = cr%o%neighbor_range
100     call glean(thy(cr))
        if (error("Exit crystal_mod::cr_neighbor_range")) continue
      end function

      function cr_neighbor_index(cr,ia,in) result(ni)
!doc$ function x_neighbor_index(cr,ia,in) result(ni)
        type(crystal_obj) :: cr
        integer, intent(in) :: ia, in
        integer :: ni
!       effects: Returns the atom index of the in neighbor of the atom ia.
!       errors: ia or in out of range.

!cod$
        call my(cr)
        if (error((ia < 1) .or. (ia > x_n_atoms(cr%o%atoms)),"ERROR: atom index out of range")) goto 100
        if (error((in < 1) .or. (in > cr%o%n_neighbors(ia)),"ERROR: neighbor index out of range")) goto 100
        ni = cr%o%neighbor_index(in,ia)
100     call glean(thy(cr))
        if (error("Exit crystal_mod::cr_neighbor_index")) continue
      end function

      function cr_neighbor_distance(cr,ia,in) result(nd)
!doc$ function x_neighbor_distance(cr,ia,in) result(nd)
        type(crystal_obj) :: cr
        integer, intent(in) :: ia, in
        real(double) :: nd
!       effects: Returns the distance from the ia atom to its in neighbor.
!       errors: ia or in out of range.

!cod$
        call my(cr)
        if (error((ia < 1) .or. (ia > x_n_atoms(cr%o%atoms)),"ERROR: atom index out of range")) goto 100
        if (error((in < 1) .or. (in > cr%o%n_neighbors(ia)),"ERROR: neighbor index out of range")) goto 100
        nd = cr%o%neighbor_distance(in,ia)
100     call glean(thy(cr))
        if (error("Exit crystal_mod::cr_neighbor_distance")) continue
      end function

      function cr_neighbor_vector(cr,ia,in) result(nv)
!doc$ function x_neighbor_vector(cr,ia,in) result(nv)
        type(crystal_obj) :: cr
        integer, intent(in) :: ia, in
        real(double), dimension(3) :: nv
!       effects: Returns the vector from the ia atom to its in neighbor in cartesian coordinates.
!       errors: ia or in out of bounds.

!cod$
        call my(cr)
        if (error((ia < 1) .or. (ia > x_n_atoms(cr%o%atoms)),"ERROR: atom index out of range")) goto 100
        if (error((in < 1) .or. (in > cr%o%n_neighbors(ia)),"ERROR: neighbor index out of range")) goto 100
        nv = cr%o%neighbor_vector(:,in,ia)
100     call glean(thy(cr))
        if (error("Exit crystal_mod::cr_neighbor_vector")) continue
      end function

      subroutine save_crys(cr,cp_mode,path)
!doc$ subroutine save(cr,cp_mode,path)
        type(crystal_obj) :: cr
        logical, intent(in), optional :: cp_mode
        character(line_len), optional :: path
!       modifies: file system
!       effects: Saves crystal information to a file.
!       errors: File I/O problems.

!cod$      
        logical :: cartesian_mode, exist_file
        character(tag_sz) :: tag
        integer :: ia, ios, na
        real(double) :: latc
        real(double), dimension(3) :: v
        type(file_obj) :: f

        call my(cr)

        if (present(cp_mode)) then
          cartesian_mode = cp_mode
        else
          cartesian_mode = .false.
        end if
        if (present(path)) then
          call my(file(trim(path)),f)
          if (i_access(f)) open(x_unit(f),file=x_name(f),status='replace',iostat=ios)
        else
          call my(file(trim(new_crystal_path)),f)
          if (i_access(f)) then
            inquire(file=x_name(f),exist=exist_file)
            if (exist_file) then
              open(x_unit(f),file=x_name(f),position='append',status='old',iostat=ios)
            else
              open(x_unit(f),file=x_name(f),status='new',iostat=ios)
            end if
          end if
        end if
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100
        if (i_access(f)) write(x_unit(f),'(a)') trim(cr%o%name)
        latc = x_lattice_constant(cr%o%lattice)
        if (i_access(f)) write(x_unit(f),'(f15.10)') latc
        v = lat2r(cr%o%lattice,real((/1,0,0/),double))/latc
        if (i_access(f)) write(x_unit(f),'(3f17.10)') v
        v = lat2r(cr%o%lattice,real((/0,1,0/),double))/latc
        if (i_access(f)) write(x_unit(f),'(3f17.10)') v
        v = lat2r(cr%o%lattice,real((/0,0,1/),double))/latc
        if (i_access(f)) write(x_unit(f),'(3f17.10)') v
        na = x_n_atoms(cr%o%atoms)
        if (i_access(f)) then
          if (cartesian_mode) then
            write(x_unit(f),'(a)') CARTESIAN_STR
            write(x_unit(f),'(i4)') na
            do ia = 1,na
              tag = x_type(cr%o%atoms,ia)
              write(x_unit(f),'(a8,3f17.10)') trim(tag), lat2r(cr%o%lattice,x_position(cr%o%atoms,ia))
            end do
          else
            write(x_unit(f),'(a)') LATTICE_STR
            write(x_unit(f),'(i4)') na
            do ia = 1,na
              tag = x_type(cr%o%atoms,ia)
              write(x_unit(f),'(a8,3f17.10)') trim(tag), x_position(cr%o%atoms,ia)
            end do
          end if
        end if
        if (i_access(f)) close(x_unit(f))

100     call glean(thy(cr))
        call glean(thy(f))

        if (error("Exit crystal_mod::save_crys")) continue

      end subroutine

      subroutine diary_cr(cr)
!doc$ subroutine diary(cr)
        type(crystal_obj) :: cr
!       modifies: output stream
!       effects: Prints crystal data to diary.
!       errors: Passes errors.

!cod$
        call my(cr)
        call diary(cr%o%lattice)
        call diary(cr%o%atoms,cr%o%lattice) ; if (error()) goto 100
100     call glean(thy(cr))
        if (error("Exit crystal_mod::diary_cr")) continue
      end subroutine

      subroutine diary_crystal_step(cr,step)
!doc$ subroutine diary(cr,step)
        type(crystal_obj) :: cr
        integer, intent(in) :: step
!       modifies: output stream
!       effects: Prints atom-related crystal data to diary for step in an atom relaxation sequence.
!       errors: Passes errors.

!cod$
        call my(cr)
        if (i_access( diaryfile() )) then
          if (step < 10) then
            write(x_unit(diaryfile()),'(/,t2,"Atom information for step #",i1,":")') step
          elseif (step < 100) then
            write(x_unit(diaryfile()),'(/,t2,"Atom information for step #",i2,":")') step
          elseif (step < 1000) then
            write(x_unit(diaryfile()),'(/,t2,"Atom information for step #",i3,":")') step
          else
            write(x_unit(diaryfile()),'(/,t2,"Atom information for step #",i4,":")') step
          end if
        end if
        call diary(cr%o%atoms,cr%o%lattice) ; if (error()) goto 100
100     call glean(thy(cr))
        if (error("Exit crystal_mod::diary_crystal_step")) continue
      end subroutine

      subroutine write_restart_cr(cr,nrestf)
!doc$ subroutine write_restart(cr,nrestf)
        type(crystal_obj) :: cr
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes cr information to nrestf.

!cod$      
        integer(long) :: dsize, iosl, ndata

        call my(cr)
        call my(nrestf)

        ! start the CRYSTAL block
        if (i_access(nrestf)) call startblock(nrestf,"CRYSTAL")

        ! write the NAME
        if (i_access(nrestf)) then
          call writetag(nrestf,"NAME")
          dsize = 1 ; ndata = line_len
          call writef(cr%o%name,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! write the lattice information
        call write_restart(cr%o%lattice,nrestf)

        ! write the atoms information
        call write_restart(cr%o%atoms,nrestf)

        ! end the CRYSTAL block
        if (i_access(nrestf)) call endblock(nrestf)

        call glean(thy(cr))
        call glean(thy(nrestf))

        if (error("Exit crystal_mod::write_restart_cr")) continue

      end subroutine

! private routines

      subroutine own_i(cr)
        type(crystal_obj) :: cr
        type(crystal_obj) :: crt
        if (cr%ref < cr%o%ref) then
          allocate( crt%o )
          crt%o%ref = 0
          crt%o%g = cr%o%g
          crt%o%name = cr%o%name
          call my(cr%o%atoms,crt%o%atoms)
          call my(cr%o%lattice,crt%o%lattice)
          crt%o%neighbor_range = cr%o%neighbor_range
          crt%o%neighbor_list_ghost = cr%o%neighbor_list_ghost
          allocate (crt%o%n_neighbors(size(cr%o%n_neighbors)) )
          crt%o%n_neighbors = cr%o%n_neighbors
          allocate (crt%o%neighbor_index(size(cr%o%neighbor_index,1),size(cr%o%neighbor_index,2)))
          crt%o%neighbor_index = cr%o%neighbor_index
          allocate (crt%o%neighbor_distance(size(cr%o%neighbor_distance,1),size(cr%o%neighbor_distance,2)) )
          crt%o%neighbor_distance = cr%o%neighbor_distance
          allocate (crt%o%neighbor_vector(size(cr%o%neighbor_vector,1),size(cr%o%neighbor_vector,2), size(cr%o%neighbor_vector,3)) )
          crt%o%neighbor_vector = cr%o%neighbor_vector
          cr%o%ref = cr%o%ref - cr%ref
          cr%o => crt%o
          cr%o%ref = cr%o%ref + cr%ref
        end if
      end subroutine

      subroutine find_neighbors_i(cr)
        type(crystal_obj) :: cr
!       effects : generates a nearest-neighbor list and writes it to the diary

        integer :: i, i1, i2, i3, ia, idx_t, inc, j, m1, m2, m3, na, nat
        real(double) :: sep_t
        real(double), dimension(3) :: vec_t
        real(double), dimension(3) :: a1, a2, a3
        integer, allocatable, dimension(:) :: idx
        real(double), allocatable, dimension(:) :: sep
        real(double), allocatable, dimension(:,:) :: vec
        real(double), allocatable, dimension(:,:) :: atom_x
        integer :: max_neighbors

        call my(cr)

        a1 = lat2r(cr%o%lattice,real((/1,0,0/),double))
        a2 = lat2r(cr%o%lattice,real((/0,1,0/),double))
        a3 = lat2r(cr%o%lattice,real((/0,0,1/),double))
        m1 = ceiling( cr%o%neighbor_range*norm(cross_product(a2,a3))/x_cell_volume(cr%o%lattice) )
        m2 = ceiling( cr%o%neighbor_range*norm(cross_product(a3,a1))/x_cell_volume(cr%o%lattice) )
        m3 = ceiling( cr%o%neighbor_range*norm(cross_product(a1,a2))/x_cell_volume(cr%o%lattice) )

        na = x_n_atoms(cr%o%atoms)
        allocate( atom_x(3,na) )
        do ia = 1,na
          atom_x(:,ia) = lat2r(cr%o%lattice,x_position(cr%o%atoms,ia))
        end do
        do ia = 1,na
          nat = 0
          do idx_t = 1,na
            do i1 = -m1,+m1
            do i2 = -m2,+m2
            do i3 = -m3,+m3
              if ((idx_t == ia) .and. (i1 == 0) .and. (i2 == 0) .and. (i3 == 0)) cycle
              vec_t = atom_x(:,idx_t) - atom_x(:,ia) + a1*real(i1,double) + a2*real(i2,double) + a3*real(i3,double)
              sep_t = norm(vec_t)
              if (sep_t <= cr%o%neighbor_range) nat = nat + 1
            end do
            end do
            end do
          end do
          cr%o%n_neighbors(ia) = nat
        end do
        max_neighbors = maxval(cr%o%n_neighbors)

        allocate( idx(max_neighbors), sep(max_neighbors), vec(3,max_neighbors) )
        if (max_neighbors /= size(cr%o%neighbor_index, 1)) then
           deallocate( cr%o%neighbor_index )    ; allocate( cr%o%neighbor_index(max_neighbors,na) )
           deallocate( cr%o%neighbor_distance ) ; allocate( cr%o%neighbor_distance(max_neighbors,na) )
           deallocate( cr%o%neighbor_vector )   ; allocate( cr%o%neighbor_vector(3,max_neighbors,na) )
        end if

        idx = 0
        sep = 0.0_double
        vec = 0.0_double
        do ia = 1,na
          nat = 0
          do idx_t = 1,na
            do i1 = -m1,+m1
            do i2 = -m2,+m2
            do i3 = -m3,+m3
              if ((idx_t == ia) .and. (i1 == 0) .and. (i2 == 0) .and. (i3 == 0)) cycle
              vec_t = atom_x(:,idx_t) - atom_x(:,ia) + a1*real(i1,double) + a2*real(i2,double) + a3*real(i3,double)
              sep_t = norm(vec_t)
              if (sep_t <= cr%o%neighbor_range) then
                nat = nat + 1
                idx(nat) = idx_t
                sep(nat) = sep_t
                vec(:,nat) = vec_t
              end if
            end do
            end do
            end do
          end do
          inc = 1
          do
            inc = 3*inc + 1
            if (inc > nat) exit
          end do
          do
            inc = inc/3
            do i = inc+1,nat
              sep_t = sep(i)
              vec_t = vec(:,i)
              idx_t = idx(i)
              j = i
              do
                if (sep(j-inc) <= sep_t) exit
                sep(j) = sep(j-inc)
                vec(:,j) = vec(:,j-inc)
                idx(j) = idx(j-inc)
                j = j - inc
                if (j <= inc) exit
              end do
              sep(j) = sep_t
              vec(:,j) = vec_t
              idx(j) = idx_t
            end do
            if (inc <= 1) exit
          end do
          cr%o%neighbor_index(:,ia) = idx(:)
          cr%o%neighbor_distance(:,ia) = sep(:)
          cr%o%neighbor_vector(:,:,ia) = vec(:,:)
        end do

        cr%o%neighbor_list_ghost = cr%o%g

        if (allocated( idx )) deallocate( idx)
        if (allocated( vec )) deallocate( vec)
        if (allocated( sep )) deallocate( sep)
        if (allocated( atom_x )) deallocate( atom_x)

        call glean(thy(cr))

      end subroutine

      end module

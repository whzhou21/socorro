!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module symmetry_mod
!doc$ module symmetry_mod

!     Two datatypes are available here: type(point_group_obj) and type(space_group_obj).

!     The space group is derived from the positions of the atoms and the point group of the
!     lattice. It is a finite subgroup of the Poincare subgroup given by point_group+translations.

!     The point group of a lattice is the set of point transformations (excluding translations)
!     that map lattice points (defined by all integer combinations of the lattice vectors) into
!     themselves. It is a finite subgroup of the orthogonal group in real space.  The point group
!     of a space group is the set of point transformations. It may or may not map the set of atoms
!     into themselves depending on whether the corresponding translations are or are not zero.

!     The purpose of point and space groups is to provide facilities for symmetrizing grids, atoms
!     vectors (vectors at atom positions such as forces) and space tensors (such as the stress tensor).

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use ghost_mod
      use diary_mod
      use arg_mod
      use lattice_mod
      use math_mod
      use atoms_mod
      use grid_mod
      use layout_mod
      use timing_mod

!cod$
      implicit none
      private

      integer, parameter :: l_max = 6

      integer, parameter :: OFF  = 0
      integer, parameter :: USER = 1
      integer, parameter :: AUTO = 2
      integer, parameter :: FULL = 3

      integer, parameter :: LATT = 10
      integer, parameter :: SGRP = 11

      real(double), parameter :: tol_identity = 1.0e-8_double
      real(double), parameter :: tol_equal = 1.0e-8_double

      type :: point_trans
        real(double), dimension(3,3) :: point_op
      end type

      type :: point_group_rep
        integer :: ref
        type(ghost) :: g
        logical :: save                                  ! save switch
        integer :: mode                                  ! mode (AUTO or USER)
        integer :: source                                ! source (LATT or SGRP)
        type(ghost) :: g_source                          ! source ghost
        type(point_trans), dimension(:), pointer :: sym  ! point group operations in the lattice representation
      end type

      type, public :: point_group_obj
        private
        integer :: ref
        type(point_group_rep), pointer :: o
      end type

      type :: space_trans
        real(double), dimension(3,3) :: point_op
        real(double), dimension(:,:), pointer :: translation
      end type

      type :: atoms_trans
        integer, dimension(:,:), pointer :: perm
      end type

      type :: field_trans
        integer :: phase_class
        integer, dimension(:), pointer :: perm
      end type

      type :: field_phase
        complex(double), dimension(:), pointer :: phase
      end type

      type :: angular_momentum
        complex(double), dimension(:,:), pointer :: mat
      end type

      type :: local_trans
        type(angular_momentum), dimension(:), pointer :: l
      end type

      type :: space_group_rep
        integer :: ref
        type(ghost) :: g
        logical :: save                                    ! save switch
        integer :: mode                                    ! mode (FULL, AUTO, or OFF)
        type(ghost) :: g_point_group                       ! point-group ghost (used only with FULL and AUTO modes)
        type(ghost) :: g_atoms                             ! atoms ghost                     "
        type(ghost) :: g_layout                            ! layout ghost                    "
        real(double) :: tol_aperm                          ! tolerance for atom permutations
        type(space_trans), dimension(:), pointer :: sym    ! space-group operations in the lattice representation
        integer, dimension(:,:), pointer :: gsmap          ! 3-d <-> 1-d G-vector map (used with non-trivial space groups)
        integer, dimension(:), pointer :: gscount          ! mesh points                            "
        integer, dimension(:), pointer :: gsdisp           ! mesh point offsets                     "
        type(atoms_trans), dimension(:), pointer :: atom   ! atom transformations                   "
        type(field_trans), dimension(:), pointer :: field  ! field transformations                  "
        type(field_phase), dimension(:), pointer :: class  ! field phases                           "
        type(local_trans), dimension(:), pointer :: local  ! local transformations                  "
      end type

      type, public :: space_group_obj
        private
        integer :: ref
        type(space_group_rep), pointer :: o
      end type

!doc$
      public :: point_group
      public :: space_group
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_n_operations
      public :: trivial
      public :: symmetrize_coordinates
      public :: symmetrize_vectors
      public :: symmetrize_tensor
      public :: symmetrize_spherical_tensor
      public :: symmetrize_grid
      public :: closure
      public :: reduction
      public :: kpoint_map
      public :: check_monkhorst_pack
      public :: invariant_site
      public :: symmetry_unique_atoms
      public :: distribute_spherical_tensor
      public :: distribute_energy
      public :: diary
      public :: diary_atom_stars
      public :: save

!cod$
      interface point_group
        module procedure constructor_pg_lat, constructor_pg_sg
      end interface
      interface space_group
        module procedure constructor_sg
      end interface
      interface update
        module procedure update_pg_lat, update_pg_sg, update_sg
      end interface
      interface my
        module procedure my_pg, my_sg, my_new_pg, my_new_sg
      end interface
      interface thy
        module procedure thy_pg, thy_sg
      end interface
      interface glean
        module procedure glean_pg, glean_sg
      end interface
      interface bequeath
        module procedure bequeath_pg, bequeath_sg
      end interface
      interface assignment(=)
        module procedure assign_pg, assign_sg
      end interface
      interface x_ref
        module procedure pg_ref, sg_ref
      end interface
      interface x_ghost
        module procedure pg_ghost, sg_ghost
      end interface
      interface x_n_operations
        module procedure pg_n_operations, sg_n_operations
      end interface
      interface trivial
        module procedure trivial_pg, trivial_sg
      end interface
      interface symmetrize_coordinates
        module procedure symmetrize_coordinates_sg
      end interface
      interface symmetrize_vectors
        module procedure symmetrize_vectors_sg
      end interface
      interface symmetrize_tensor
        module procedure symmetrize_tensor_sg
      end interface
      interface symmetrize_spherical_tensor
        module procedure symmetrize_spherical_tensor_sg
      end interface
      interface symmetrize_grid
        module procedure symmetrize_grid_sg
      end interface
      interface closure
        module procedure closure_pg
      end interface
      interface reduction
        module procedure reduction_pg
      end interface
      interface kpoint_map
        module procedure kpoint_map_pg
      end interface
      interface check_monkhorst_pack
        module procedure check_monkhorst_pack_pg
      end interface
      interface invariant_site
        module procedure invariant_site_sg
      end interface
      interface symmetry_unique_atoms
        module procedure symmetry_unique_atoms_sg
      end interface
      interface distribute_spherical_tensor
        module procedure distribute_spherical_tensor_sg
      end interface
      interface distribute_energy
        module procedure distribute_energy_sg
      end interface
      interface diary
        module procedure diary_pg, diary_sg
      end interface
      interface diary_atom_stars
        module procedure diary_atom_stars_sg
      end interface
      interface save
        module procedure save_pg, save_sg
      end interface

      contains

! point group public routines

      function constructor_pg_lat(lat) result(pg)
!doc$ function point_group(lat) result(pg)
        type(lattice_obj) :: lat
        type(point_group_obj) :: pg
!       requires: Magnitudes of lat vectors related by a symmetry operator can differ by no more than tol_vector.
!       effects: Constructs a new pg from lat or reads pg from a file.
!       errors: Invalid lattice group. Passes errors.

!cod$
        logical :: fnd
        character(tag_sz) :: tag

        call my(lat)

        pg%ref = 0
        allocate( pg%o )
        pg%o%ref = 0
        pg%o%g = x_ghost()
        pg%o%source = LATT

        call arglc("lattice_symmetry",tag,fnd)
        if (.not.fnd) tag = "auto"
        select case(trim(tag))
        case("auto")
          pg%o%mode = AUTO
        case("user")
          pg%o%mode = USER
        case default
          if (error(.true.,"ERROR: lattice_symmetry tag is not recognized")) goto 100
        end select

        call arg("save_lattice_group",pg%o%save,fnd)
        if (.not.fnd) pg%o%save = .false.

        nullify( pg%o%sym )
        select case (pg%o%mode)
        case (AUTO)
          call generate_lattice_group_i(pg%o,lat) ; if (error()) goto 100
        case (USER)
          call read_lattice_group_i(pg%o,lat) ; if (error()) goto 100
          if (error(.not.valid_lattice_group_i(pg%o,lat),"ERROR: User-supplied lattice group is not valid")) goto 100
        end select

        if (pg%o%save) call save_pg_i(pg%o)

100     call glean(thy(lat))

        if (error("Exit symmetry_mod::constructor_pg_lat")) continue

      end function

      subroutine update_pg_lat(pg,lat)
!doc$ subroutine update(pg,lat)
        type(point_group_obj) :: pg
        type(lattice_obj) :: lat
!       modifies: pg
!       effects: Updates pg with respect to lat.
!                  For mode = AUTO, updates pg if lat symmetry is lower than pg symmetry.
!                  For mode = USER, confirms that lat symmetry is lower than pg symmetry.
!       errors: Incorrect source. Invalid pg. Passes errors.

!cod$
        logical :: remake, source_change

        call my(pg)
        call my(lat)

        select case (pg%o%source)
        case (LATT)
          source_change = (pg%o%g_source /= x_ghost(lat))
          select case (pg%o%mode)
          case (AUTO)
            remake = .false.
            if (source_change) remake = (remake .or. .not.valid_lattice_group_i(pg%o,lat))
            if (remake) then
              call own_pg_i(pg)
              pg%o%g = x_ghost()
              call generate_lattice_group_i(pg%o,lat) ; if (error()) goto 100
              if (pg%o%save) call save_pg_i(pg%o)
            end if
          case (USER)
            if (source_change) then
              if (error(.not.valid_lattice_group_i(pg%o,lat),"ERROR: User-supplied lattice group is not valid")) goto 100
            end if
          end select
        case default
          if (error(.true.,"ERROR: incorrect source")) goto 100
        end select

100     call glean(thy(pg))
        call glean(thy(lat))

        if (error("Exit symmetry_mod::update_pg_lat")) continue

      end subroutine

      function constructor_pg_sg(sg,parity) result(pg)
!doc$ function point_group(sg,parity) result(pg)
        type(space_group_obj) :: sg
        logical, intent(in), optional :: parity
        type(point_group_obj) :: pg
!       effects: Constructs a new pg from sg. If parity is present and true, inversion is added if it is not already present.

!cod$
        logical :: par

        call my(sg)

        pg%ref = 0
        allocate( pg%o )
        pg%o%ref = 0
        pg%o%g = x_ghost()
        pg%o%source = SGRP

        pg%o%mode = FULL
        pg%o%save = .false.

        par = .false.
        if (present(parity)) par = parity

        nullify( pg%o%sym )
        call extract_point_group_i(pg%o,sg,par)

        call glean(thy(sg))

        if (error("Exit symmetry_mod::constructor_pg_sg")) continue

      end function

      subroutine update_pg_sg(pg,sg,parity)
!doc$ subroutine update(pg,sg,parity)
        type(point_group_obj) :: pg
        type(space_group_obj) :: sg
        logical, intent(in), optional :: parity
!       modifies: pg
!       effects: Updates pg with respect to sg. If parity is present and true, inversion is added if it is not already present.
!       errors: Incorrect source.

!cod$
        logical :: par, parity_change, source_change
        type(point_group_obj) :: pgn

        call my(pg)
        call my(sg)

        select case (pg%o%source)
        case (SGRP)
          select case (pg%o%mode)
          case (FULL)
            source_change = (pg%o%g_source /= x_ghost(sg))
            parity_change = .false.
            if (present(parity)) parity_change = (parity .neqv. parity_space_group_i(sg))
            if (source_change) then
              pgn%ref = 0
              allocate( pgn%o )
              pgn%o%ref = 0
              pgn%o%g = x_ghost()
              pgn%o%source = pg%o%source
              pgn%o%mode = pg%o%mode
              pgn%o%save = pg%o%save
              if (present(parity)) then
                par = parity
              else
                par = parity_space_group_i(sg)
              end if
              nullify( pgn%o%sym )
              call extract_point_group_i(pgn%o,sg,par)
              if (same_point_group_i(pg%o,pgn%o)) then
                pg%o%g_source = x_ghost(sg)
              else
                pg = pgn
              end if
              call glean(pgn)
            else
              if (parity_change) then
                call own_pg_i(pg)
                pg%o%g = x_ghost()
                call extract_point_group_i(pg%o,sg,parity)
              end if
            end if
          end select
        case default
          if (error(.true.,"ERROR: incorrect source")) goto 100
        end select

100     call glean(thy(pg))
        call glean(thy(sg))

        if (error("Exit symmetry_mod::update_pg_sg")) continue

      end subroutine

      subroutine my_pg(pg)
!doc$ subroutine my(pg)
        type(point_group_obj) :: pg

!cod$
        pg%ref = pg%ref + 1
        pg%o%ref = pg%o%ref + 1
      end subroutine

      subroutine my_new_pg(pgi,pg)
!doc$ subroutine my(pgi,pg)
        type(point_group_obj) :: pgi, pg

!cod$
        pg%ref = 1
        pg%o => pgi%o
        pg%o%ref = pg%o%ref + 1
      end subroutine

      function thy_pg(pg) result(pgo)
!doc$ function thy(pg) result(pgo)
        type(point_group_obj) :: pg, pgo

!cod$
        pg%ref = pg%ref - 1
        pg%o%ref = pg%o%ref - 1
        pgo%ref = pg%ref
        pgo%o => pg%o
      end function

      subroutine glean_pg(pg)
!doc$ subroutine glean(pg)
        type(point_group_obj) :: pg

!cod$
        if (pg%o%ref < 1) then
          if (associated( pg%o%sym )) deallocate( pg%o%sym )
          deallocate( pg%o )
        end if
      end subroutine

      subroutine bequeath_pg(pg)
!doc$ subroutine bequeath(pg)
        type(point_group_obj) :: pg

!cod$
        continue
      end subroutine

      subroutine assign_pg(pg,pgi)
!doc$ subroutine assignment(=)(pg,pgi)
        type(point_group_obj), intent(inout) :: pg
        type(point_group_obj), intent(in) :: pgi

!cod$
        type(point_group_obj) :: pgt
        call my(pgi)
        pgt%o => pg%o
        pg%o%ref = pg%o%ref - pg%ref
        pg%o => pgi%o
        pg%o%ref = pg%o%ref + pg%ref
        call glean(pgt)
        call glean(thy(pgi))
      end subroutine

      function pg_ref(pg) result(r)
!doc$ function x_ref(pg) result(r)
        type(point_group_obj) :: pg
        integer, dimension(2) :: r
!       effects: Returns pg%ref and pg%o%ref.

!cod$
        r(1) = pg%ref
        r(2) = pg%o%ref
        call glean(pg)
      end function

      function pg_ghost(pg) result(g)
!doc$ function x_ghost(pg) result(g)
        type(point_group_obj) :: pg
        type(ghost) :: g
!       effects: Returns the ghost of pg.

!cod$
        call my(pg)
        g = pg%o%g
        call glean(thy(pg))
      end function

      function pg_n_operations(pg) result(n)
!doc$ function x_n_operations(pg) result(n)
        type(point_group_obj) :: pg
        integer :: n
!       effects: Returns the number of symmetry operations in pg.

!cod$
        call my(pg)
        n = size(pg%o%sym)
        call glean(thy(pg))
      end function

      function trivial_pg(pg) result(t)
!doc$ function trivial(pg) result(t)
        type(point_group_obj) :: pg
        logical :: t
!       effects: Returns .true. if pg is trivial.

!cod$
        call my(pg)
        select case (pg%o%source)
        case (LATT)
          t = (size(pg%o%sym) == 2)
        case (SGRP)
          t = (size(pg%o%sym) == 1)
        end select
        call glean(thy(pg))
      end function

      subroutine closure_pg(pg,lat,v)
!doc$ subroutine closure(pg,lat,v)
        type(point_group_obj) :: pg
        type(lattice_obj) :: lat
        real(double), dimension(:,:), pointer :: v
!       modifies: v
!       requires: v be in reciprocal-lattice representation and confined within the the first Brillouin zone.
!       effects: v is overwritten with the superset (no duplications) of vectors obtained by acting on v with pg.

!cod$
        logical :: mapped
        integer :: is, iv, ivt, ns, nv, nvt
        integer, dimension(:), allocatable :: mag_num, mag_num_t
        real(double), dimension(3) :: vr
        real(double), dimension(:,:), allocatable :: vt

        call my(pg)
        call my(lat)

        nv = size(v,2)
        ns = size(pg%o%sym)

        ! Compute the magnitudes of the vectors, v.
        allocate( vt(3,nv), mag_num(nv) )
        do iv = 1,nv
          vt(:,iv) = lat2f(lat,v(:,iv))
        end do
        call get_mag_list_i(vt,mag_num)
        deallocate( vt )

        allocate( vt(3,ns*nv), mag_num_t(ns*nv) )
        vt(:,1:nv) = v
        mag_num_t(1:nv) = mag_num

        nvt = nv
        do ivt = 1,nvt
          do is = 1,ns
            vr = matmul(vt(:,ivt),pg%o%sym(is)%point_op)
            mapped = .false.
            do iv = 1,nv
              if (mag_num_t(iv) /= mag_num_t(ivt)) cycle
              if (vr .in. zone(vt(:,iv),tol_equal)) then
                mapped = .true.
                exit
              end if
            end do
            if (.not.mapped) then
              nv = nv + 1
              mag_num_t(nv) = mag_num_t(ivt)
              vt(:,nv) = vr
            end if
          end do
        end do

        deallocate( v )
        allocate( v(3,nv) )
        v = vt(:,1:nv)

        if (allocated( mag_num )) deallocate( mag_num )
        if (allocated( mag_num_t )) deallocate( mag_num_t )
        if (allocated( vt )) deallocate( vt )

        call glean(thy(pg))
        call glean(thy(lat))

      end subroutine
      
      subroutine reduction_pg(pg,lat,v,wt)
!doc$ subroutine reduction(pg,lat,v,wt)
        type(point_group_obj) :: pg
        type(lattice_obj) :: lat
        real(double), dimension(:,:), pointer :: v
        integer, dimension(:), pointer :: wt
!       modifies: v and wt
!       requires: v be in reciprocal-lattice representation and confined within the first Brillouin zone.
!                 wt be nullified or allocated.
!       effects: v is overwritten with the subset (no duplications) of vectors with the property that every vector of an
!                input v is obtained by acting on one of the output v with an operator in pg. wt contains the number of
!                input v that map into the output v.
!       error: none

!cod$
        logical :: mapped
        integer :: i, j, inc, is, ns, iv, nv, iv_t, nv_t, wtr
        integer, dimension(:), allocatable :: mag, mag_t, wt_t
        real(double), dimension(3) :: vr
        real(double), dimension(:,:), allocatable :: vc, v_t

        call my(pg)
        call my(lat)

        if (associated( wt )) deallocate( wt )

        ns = size(pg%o%sym)
        nv = size(v,2)

        allocate( vc(3,nv), mag(nv) )
        do iv = 1,nv
          vc(:,iv) = lat2f(lat,v(:,iv))
        end do
        call get_mag_list_i(vc,mag)
        deallocate( vc )

        allocate( v_t(3,nv) )
        allocate( wt_t(nv) )
        allocate( mag_t(nv) )
        nv_t = 0
        do iv = 1,nv                                                             ! REDUCE V UNDER PG
          mapped = .false.
          SYM_OPS: do is = 1,ns
            vr = matmul(v(:,iv),pg%o%sym(is)%point_op)
            do iv_t = 1,nv_t
              if (mag_t(iv_t) /= mag(iv)) cycle
              if (vr .in. zone(v_t(:,iv_t),tol_equal)) then
                mapped = .true.
                wt_t(iv_t) = wt_t(iv_t) + 1
                exit SYM_OPS
              end if
            end do
          end do SYM_OPS
          if (.not.mapped) then
            nv_t = nv_t + 1
            v_t(:,nv_t) = v(:,iv)
            wt_t(nv_t) = 1
            mag_t(nv_t) = mag(iv)
          end if
        end do

        deallocate( v )
        nv = nv_t
        allocate( v(3,nv) )
        allocate( wt(nv) )
        v = v_t(:,1:nv)                                                          ! STORE THE RESULTS
        wt = wt_t(1:nv)

        inc = 1                                                                  ! SORT ACCORTING TO WT
        do
          inc = 3*inc + 1
          if (inc > nv) exit
        end do
        do
          inc = inc/3
          do i = inc+1,nv
            vr = v(:,i)
            wtr = wt(i)
            j = i
            do
              if (wt(j-inc) <= wtr) exit
              v(:,j) = v(:,j-inc)
              wt(j) = wt(j-inc)
              j = j - inc
              if (j <= inc) exit
            end do
            v(:,j) = vr
            wt(j) = wtr
          end do
          if (inc <= 1) exit
        end do

100     if (allocated( mag )) deallocate( mag )
        if (allocated( mag_t )) deallocate( mag_t )
        if (allocated( wt_t )) deallocate( wt_t )
        if (allocated( vc )) deallocate( vc )
        if (allocated( v_t )) deallocate( v_t )

        call glean(thy(pg))
        call glean(thy(lat))

        if (error("Exit symmetry_mod::reduction_pg")) continue

      end subroutine

      subroutine kpoint_map_pg(pg,lat,v,vm)
!doc$ subroutine kpoint_map(pg,lat,v,vm)
        type(point_group_obj) :: pg
        type(lattice_obj) :: lat
        real(double), dimension(:,:), pointer :: v
        logical, dimension(:,:), pointer :: vm
!       modifies: vm
!       requires: v be a set of Generalized Monkhorst-Pack (gmp) vectors in lattice coordinates.
!                 vm be nullified or allocated.
!       effects: vm gives the mapping of the vectors (within the set) under pg.
!       error: none

!cod$
        integer :: is, ns, iv, nv, ivt
        integer, dimension(:), allocatable :: mag
        real(double), dimension(3) :: vr
        real(double), dimension(:,:), allocatable :: vc

        call my(pg)
        call my(lat)

        if (associated( vm )) deallocate( vm )

        ns = size(pg%o%sym)
        nv = size(v,2)

        allocate( vc(3,nv) )
        allocate( mag(nv) )
        do iv = 1,nv
          vc(:,iv) = lat2f(lat,v(:,iv))
        end do
        call get_mag_list_i(vc,mag)
        deallocate( vc )

        allocate( vm(nv,nv) )
        vm = .false.
        do iv = 1,nv
          do is = 1,ns
            vr = matmul(v(:,iv),pg%o%sym(is)%point_op)
            do ivt = 1,nv
              if (mag(ivt) /= mag(iv)) cycle
              if (vr .in. zone(v(:,ivt),tol_equal)) then
                vm(iv,ivt) = .true.
              end if
            end do
          end do
        end do

100     if (allocated( mag )) deallocate( mag )
        if (allocated( vc )) deallocate( vc )

        call glean(thy(pg))
        call glean(thy(lat))

        if (error("Exit symmetry_mod::kpoint_map_pg")) continue

      end subroutine

      subroutine check_monkhorst_pack_pg(pg,q1,q2,q3)
!doc$ subroutine check_monkhorst_pack(pg,q1,q2,q3)
        type(point_group_obj) :: pg
        integer, intent(inout) :: q1, q2, q3
!       requires: pg be generated from a lattice.
!       effects: Issues a warning if Monkhorst-Pack parameters are different for two reciprocal lattice vectors related
!                by an operation in pg.

!cod$
        logical :: related_12, related_13, related_23
        integer :: is, ns
        real(double), dimension(3) :: b1, b2, b3

        call my(pg)

        ns = size(pg%o%sym)

        b1 = real((/1,0,0/),double)
        b2 = real((/0,1,0/),double)
        b3 = real((/0,0,1/),double)

        related_12 = .false.
        do is = 1,ns
          related_12 = related_12 .or. (matmul(b1,pg%o%sym(is)%point_op) .in. nbhd(b2,tol_equal))
          if (related_12) exit
        end do
        related_13 = .false.
        do is = 1,ns
          related_13 = related_13 .or. (matmul(b1,pg%o%sym(is)%point_op) .in. nbhd(b3,tol_equal))
          if (related_13) exit
        end do
        related_23 = .false.
        do is = 1,ns
          related_23 = related_23 .or. (matmul(b2,pg%o%sym(is)%point_op) .in. nbhd(b3,tol_equal))
          if (related_23) exit
        end do

        if ( related_12 .and. (q1 /= q2) ) call warn("WARNING: q1 should equal q2")
        if ( related_13 .and. (q1 /= q3) ) call warn("WARNING: q1 should equal q3")
        if ( related_23 .and. (q2 /= q3) ) call warn("WARNING: q2 should equal q3")

100     call glean(thy(pg))

        if (error("Exit symmetry_mod::check_monkhorst_pack_pg")) continue

      end subroutine

      subroutine diary_pg(pg,list,lat)
!doc$ subroutine diary(pg,list,lat)
        type(point_group_obj) :: pg
        logical, intent(in), optional :: list
        type(lattice_obj), optional :: lat
!       modifies: diary stream
!       effects: Diaries pg. If list is present and true, pg operations are listed.
!       errors: Passes errors.

!cod$
        character(4) :: pg_type
        character(60) :: lattice_type
        logical :: fnd, list_local
        integer :: is, ns
        integer, dimension(10) :: op_type
        real(double), dimension(3,3) :: tm

        call my(pg)
        if (present(lat)) call my(lat)

        if (present(list)) then
          list_local = list
        else
          select case (pg%o%source)
          case (LATT)
            call arg("list_lattice_group",list_local,fnd)
            if (.not.fnd) list_local = .false.
          case (SGRP)
            list_local = .false.
          end select
        end if

        ns = size(pg%o%sym)
        op_type = 0
        do is = 1,ns
          op_type = op_type + classify_op_i(pg%o%sym(is)%point_op)
        end do

        select case (pg%o%source)
        case (LATT)
          call get_lattice_type_i(op_type,lattice_type) ; if (error()) goto 100
          if (i_access(diaryfile())) then
            select case (pg%o%mode)
            case (AUTO)
              write(x_unit(diaryfile()),'(/,t4,"Lattice symmetry (automatic):")')
            case (USER)
              write(x_unit(diaryfile()),'(/,t4,"Lattice symmetry (user supplied):")')
            end select
            write(x_unit(diaryfile()),'(/,t6,a)') lattice_type
          end if
        case (SGRP)
          call get_point_group_type_i(op_type,pg_type) ; if (error()) goto 100
          if (i_access(diaryfile())) then
            if (ns == 1) then
              write(x_unit(diaryfile()),'(/,t6,"The point group is ",a," with ",i0," operation")') trim(pg_type), ns
            else
              write(x_unit(diaryfile()),'(/,t6,"The point group is ",a," with ",i0," operations ")') trim(pg_type), ns
            end if
          end if
        end select
        if (list_local) then
          if (i_access(diaryfile())) then
            if (present(lat)) then
              write(x_unit(diaryfile()),'(/,t6,"Point-group operations in the cartesian representation:")')
              do is = 1,ns
                tm = lat2r(lat,pg%o%sym(is)%point_op)
                write(x_unit(diaryfile()),'(/,t16,3f6.1)') tm(1,:)
                write(x_unit(diaryfile()),'(t6,i2,".",7x,3f6.1)') is, tm(2,:)
                write(x_unit(diaryfile()),'(t16,3f6.1)') tm(3,:)
              end do
            else
              write(x_unit(diaryfile()),'(/,t6,"Point-group operations in the lattice representation:")')
              do is = 1,ns
                write(x_unit(diaryfile()),'(/,t16,3f6.1)') pg%o%sym(is)%point_op(1,:)
                write(x_unit(diaryfile()),'(t6,i2,".",7x,3f6.1)') is, pg%o%sym(is)%point_op(2,:)
                write(x_unit(diaryfile()),'(t16,3f6.1)') pg%o%sym(is)%point_op(3,:)
              end do
            end if
          end if
        end if

100     call glean(thy(pg))
        if (present(lat)) call glean(thy(lat))

        if (error("Exit symmetry_mod::diary_pg")) continue

      end subroutine

      subroutine save_pg(pg)
!doc$ subroutine save(pg)
        type(point_group_obj) :: pg
!       modifies: file system
!       effects: Saves the pg symmetry operations.
!       errors: File I/O problems.

!cod$
        integer :: ios, is, ns
        type(file_obj) :: f

        call my(pg)

        select case (pg%o%source)
        case (LATT)
          call my(file(trim(new_lattice_group_path)),f)
        case (SGRP)
          call my(file(trim(new_point_group_path)),f)
        end select

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100

        ns = size(pg%o%sym)
        if (i_access(f)) then
          write(x_unit(f),'(i0)') ns
          do is = 1,ns
            write(x_unit(f),'(a)') " "
            write(x_unit(f),'(3i5)') nint(pg%o%sym(is)%point_op(1,:))
            write(x_unit(f),'(3i5)') nint(pg%o%sym(is)%point_op(2,:))
            write(x_unit(f),'(3i5)') nint(pg%o%sym(is)%point_op(3,:))
          end do
        end if

        if (i_access(f)) close(x_unit(f))

100     call glean(thy(pg))
        call glean(thy(f))

        if (error("Exit symmetry_mod::save_pg")) continue

      end subroutine

! space group public routines

      function constructor_sg(pg,ats,lat,lay) result(sg)
!doc$ function space_group(pg,ats,lat,lay) result(sg)
        type(point_group_obj) :: pg
        type(atoms_obj) :: ats
        type(lattice_obj) :: lat
        type(layout_obj), optional :: lay
        type(space_group_obj) :: sg
!       requires: pg be a subgroup of the group associated with lat. lay be present unless the call is from check_symmetry.
!       effects: Constructs a new sg. Symmetrization structures are formed only if sg is not trivial.
!       errors: Unrecognized symmetry tag. Passes errors.

!cod$
        logical :: fnd
        character(tag_sz) :: tag
        real(double) :: tol_in
        real(double), dimension(3) :: mag

        call my(pg)
        call my(ats)
        call my(lat)
        if (present(lay)) call my(lay)

        sg%ref = 0
        allocate( sg%o )
        sg%o%ref = 0
        sg%o%g = x_ghost()

        call arglc("symmetry",tag,fnd)
        if (.not.fnd) tag = "auto"
        select case(trim(tag))
        case("full")
          sg%o%mode = FULL
        case("auto")
          sg%o%mode = AUTO
        case("off")
          sg%o%mode = OFF
        case default
          if (error(.true.,"ERROR: symmetry tag was not recognized")) goto 100
        end select

        call arg("sg_tolerance",tol_in,fnd)
        if (fnd) then
          mag(1) = norm(lat2r(lat,real((/1,0,0/),double)))
          mag(2) = norm(lat2r(lat,real((/0,1,0/),double)))
          mag(3) = norm(lat2r(lat,real((/0,0,1/),double)))
          sg%o%tol_aperm = tol_in/minval(mag)
          if (sg%o%tol_aperm > 1.0e-4_double) call warn("WARNING: the space-group tolerance is very loose")
        else
          sg%o%tol_aperm = 1.0e-8_double
        end if

        call arg("save_space_group",sg%o%save,fnd)
        if (.not.fnd) sg%o%save = .false.

        nullify( sg%o%sym )
        nullify( sg%o%atom )
        nullify( sg%o%field )
        nullify( sg%o%class )
        nullify( sg%o%gsmap )
        nullify( sg%o%gscount )
        nullify( sg%o%gsdisp )
        nullify( sg%o%local )

        select case (sg%o%mode)
        case (FULL,AUTO)
          if (present(lay)) then
            call generate_space_group_i(sg%o,pg,ats,lat,lay) ; if (error()) goto 100
          else
            call generate_space_group_i(sg%o,pg,ats,lat) ; if (error()) goto 100
          end if
        case (OFF)
          call form_trivial_space_group_i(sg%o)
        end select

        if (sg%o%save) call save_sg_i(sg%o)

100     call glean(thy(pg))
        call glean(thy(ats))
        call glean(thy(lat))
        if (present(lay)) call glean(thy(lay))

        if (error("Exit symmetry_mod::constructor_sg")) continue

      end function

      subroutine update_sg(sg,pg,ats,lat,lay)
!doc$ subroutine update(sg,pg,ats,lat,lay)
        type(space_group_obj) :: sg
        type(point_group_obj) :: pg
        type(atoms_obj) :: ats
        type(lattice_obj) :: lat
        type(layout_obj) :: lay
!       modifies: sg
!       requires: pg be a subgroup of the group associated with lat.
!       effects: Updates sg with respect to its dependencies.
!       errors: Passes errors.

!cod$
        logical :: atoms_change, layout_change, point_group_change, remake
        type(space_group_obj) :: sgn

        call my(sg)
        call my(pg)
        call my(ats)
        call my(lat)
        call my(lay)

        select case (sg%o%mode)
        case (FULL)
          point_group_change = (sg%o%g_point_group /= x_ghost(pg))
          atoms_change = (sg%o%g_atoms /= x_ghost(ats))
          layout_change = (sg%o%g_layout /= x_ghost(lay))
          if (point_group_change .or. layout_change) then
            call own_sg_i(sg)
            sg%o%g = x_ghost()
            call generate_space_group_i(sg%o,pg,ats,lat,lay) ; if (error()) goto 100
            if (sg%o%save) call save_sg_i(sg%o)
          else
            if (atoms_change) then
              sgn%ref = 0
              allocate( sgn%o )
              sgn%o%ref = 0
              sgn%o%g = x_ghost()
              sgn%o%mode = sg%o%mode
              sgn%o%tol_aperm = sg%o%tol_aperm
              sgn%o%save = sg%o%save
              nullify( sgn%o%sym )
              nullify( sgn%o%atom )
              nullify( sgn%o%field )
              nullify( sgn%o%class )
              nullify( sgn%o%gsmap )
              nullify( sgn%o%gscount )
              nullify( sgn%o%gsdisp )
              nullify( sgn%o%local )
              call generate_space_group_i(sgn%o,pg,ats,lat,lay) ; if (error()) goto 100
              if (same_space_group_i(sg%o,sgn%o)) then
                sg%o%g_atoms = x_ghost(ats)
              else
                sg = sgn
                if (sg%o%save) call save_sg_i(sg%o)
              end if
              call glean(sgn)
            end if
          end if
        case (AUTO)
          point_group_change = (sg%o%g_point_group /= x_ghost(pg))
          atoms_change = (sg%o%g_atoms /= x_ghost(ats))
          layout_change = (sg%o%g_layout /= x_ghost(lay))
          remake = .false.
          if (point_group_change .or. layout_change) remake = .true.
          if (atoms_change) remake = (remake .or. .not.valid_space_group_i(sg%o,ats))
          if (remake) then
            call own_sg_i(sg)
            sg%o%g = x_ghost()
            call generate_space_group_i(sg%o,pg,ats,lat,lay) ; if (error()) goto 100
            if (sg%o%save) call save_sg_i(sg%o)
          else
            if (atoms_change) sg%o%g_atoms = x_ghost(ats)
          end if
        end select

100     call glean(thy(sg))
        call glean(thy(pg))
        call glean(thy(ats))
        call glean(thy(lat))
        call glean(thy(lay))

        if (error("Exit symmetry_mod::update_sg")) continue

      end subroutine

      subroutine my_sg(sg)
!doc$ subroutine my(sg)
        type(space_group_obj) :: sg

!cod$
        sg%ref = sg%ref + 1
        sg%o%ref = sg%o%ref + 1
      end subroutine

      subroutine my_new_sg(sgi,sg)
!doc$ subroutine my(sgi,sg)
        type(space_group_obj) :: sgi, sg

!cod$
        sg%ref = 1
        sg%o => sgi%o
        sg%o%ref = sg%o%ref + 1
      end subroutine

      function thy_sg(sg) result(sgo)
!doc$ function thy(sg) result(sgo)
        type(space_group_obj) :: sg, sgo

!cod$
        sg%ref = sg%ref - 1
        sg%o%ref = sg%o%ref - 1
        sgo%ref = sg%ref
        sgo%o => sg%o
      end function

      subroutine glean_sg(sg)
!doc$ subroutine glean(sg)
        type(space_group_obj) :: sg

!cod$
        integer :: ic, is, l
        if (sg%o%ref < 1) then
          if (associated( sg%o%sym )) then
            do is = 1,size(sg%o%sym)
              if (associated( sg%o%sym(is)%translation )) deallocate( sg%o%sym(is)%translation )
            end do
            deallocate( sg%o%sym )
          end if
          if (associated( sg%o%gsmap )) deallocate( sg%o%gsmap )
          if (associated( sg%o%gscount )) deallocate( sg%o%gscount )
          if (associated( sg%o%gsdisp )) deallocate( sg%o%gsdisp )
          if (associated( sg%o%atom )) then
            do is = 1,size(sg%o%atom)
              if (associated( sg%o%atom(is)%perm )) deallocate( sg%o%atom(is)%perm )
            end do
            deallocate( sg%o%atom )
          end if
          if (associated( sg%o%field )) then
            do is = 1,size(sg%o%field)
              if (associated( sg%o%field(is)%perm )) deallocate( sg%o%field(is)%perm )
            end do
            deallocate( sg%o%field )
          end if
          if (associated( sg%o%class )) then
            do ic = 1,size(sg%o%class)
              if (associated( sg%o%class(ic)%phase )) deallocate( sg%o%class(ic)%phase )
            end do
            deallocate( sg%o%class )
          end if
          if (associated( sg%o%local )) then
            do is = 1,size(sg%o%local)
              if (associated( sg%o%local(is)%l )) then
                do l = 0,size(sg%o%local(is)%l)-1
                  if (associated( sg%o%local(is)%l(l)%mat )) deallocate( sg%o%local(is)%l(l)%mat )
                end do
                deallocate( sg%o%local(is)%l )
              end if
            end do
            deallocate( sg%o%local )
          end if
          deallocate( sg%o )
        end if
      end subroutine

      subroutine bequeath_sg(sg)
!doc$ subroutine bequeath(sg)
         type(space_group_obj) :: sg

!cod$
         continue
      end subroutine

      subroutine assign_sg(sg,sgi)
!doc$ subroutine assignment(=)(sg,sgi)
        type(space_group_obj), intent(inout) :: sg
        type(space_group_obj), intent(in) :: sgi
!       effects: sg = sgi

!cod$
        type(space_group_obj) :: sgt
        call my(sgi)
        sgt%o => sg%o
        sg%o%ref = sg%o%ref - sg%ref
        sg%o => sgi%o
        sg%o%ref = sg%o%ref + sg%ref
        call glean(sgt)
        call glean(thy(sgi))
      end subroutine

      function sg_ref(sg) result(r)
!doc$ function x_ref(sg) result(r)
        type(space_group_obj) :: sg
        integer, dimension(2) :: r
!       effects: Returns sg%ref and sg%o%ref.

!cod$
        r(1) = sg%ref
        r(2) = sg%o%ref
        call glean(sg)
      end function

      function sg_ghost(sg) result(g)
!doc$ function x_ghost(sg) result(g)
        type(space_group_obj) :: sg
        type(ghost) :: g
!       effects: Returns the ghost of sg.

!cod$
        call my(sg)
        g = sg%o%g
        call glean(thy(sg))
      end function

      function sg_n_operations(sg) result(n)
!doc$ function x_n_operations(sg) result(n)
        type(space_group_obj) :: sg
        integer :: n
!       effects: Returns the number of symmetry operations in sg.

!cod$
        call my(sg)
        n = size(sg%o%sym)
        call glean(thy(sg))
      end function

      function trivial_sg(sg) result(t)
!doc$ function trivial(sg) result(t)
        type(space_group_obj) :: sg
        logical :: t
!       effects: Returns .true. if sg consists of the identity operation and the zero translation.

!cod$
        call my(sg)
        t = ((size(sg%o%sym) == 1) .and. (size(sg%o%sym(1)%translation,2) == 1))
        call glean(thy(sg))
      end function

      subroutine symmetrize_coordinates_sg(sg,c)
!doc$ subroutine symmetrize_coordinates(sg,c)
        type(space_group_obj) :: sg
        real(double), dimension(:,:), intent(inout) :: c
!       modifies: c
!       requires: c be in the lattice representation.
!       effects: Overwrites c with a symmetrized version.
!       errors: Improper dimensions for c.

!cod$
        integer :: ia, iap, is, it, na, ns, nt
        real(double), dimension(3) :: ct
        real(double), dimension(:,:), allocatable :: c_sum

        call my(sg)

        if (.not.trivial(sg)) then
          ns = size(sg%o%atom)
          na = size(sg%o%atom(1)%perm,1)
          nt = size(sg%o%atom(1)%perm,2)
          if (error(size(c,1) /= 3,"ERROR: improper first dimension for c")) goto 100
          if (error(size(c,2) /= na,"ERROR: improper second dimension for c")) goto 100
          allocate( c_sum(3,na) )
          c_sum = 0.0_double
          do is = 1,ns
            do it = 1,nt
              do ia = 1,na
                iap = sg%o%atom(is)%perm(ia,it)
                ct = matmul(sg%o%sym(is)%point_op,c(:,ia)) + sg%o%sym(is)%translation(:,it)
                c_sum(:,iap) = c_sum(:,iap) + ( ct - real(nint(ct - c(:,iap)),double) )
              end do
            end do
          end do
          c = c_sum/real(ns*nt,double)
        end if

100     if (allocated( c_sum )) deallocate( c_sum )

        call glean(thy(sg))

        if (error("Exit symmetry_mod::symmetrize_coordinates_sg")) continue

      end subroutine

      subroutine symmetrize_vectors_sg(sg,v)
!doc$ subroutine symmetrize_vectors(sg,v)
        type(space_group_obj) :: sg
        real(double), dimension(:,:), intent(inout) :: v
!       modifies: v
!       requires: v be in the lattice representation.
!       effects: Overwrites v with a symmetrized version.
!       errors: Improper dimensions for v.

!cod$
        integer :: ia, iap, is, it, na, ns, nt
        real(double), dimension(3) :: vt
        real(double), dimension(:,:), allocatable :: v_sum

        call my(sg)

        if (.not.trivial(sg)) then
          ns = size(sg%o%atom)
          na = size(sg%o%atom(1)%perm,1)
          nt = size(sg%o%atom(1)%perm,2)
          if (error(size(v,1) /= 3,"ERROR: improper first dimension for v")) goto 100
          if (error(size(v,2) /= na,"ERROR: improper second dimension for v")) goto 100
          allocate( v_sum(3,na) )
          v_sum = 0.0_double
          do is = 1,ns
            do ia = 1,na
              vt = matmul(sg%o%sym(is)%point_op,v(:,ia))
              do it = 1,nt
                iap = sg%o%atom(is)%perm(ia,it)
                v_sum(:,iap) = v_sum(:,iap) + vt
              end do
            end do
          end do
          v = v_sum/real(ns*nt,double)
        end if

100     if (allocated( v_sum )) deallocate( v_sum )

        call glean(thy(sg))

        if (error("Exit symmetry_mod::symmetrize_vectors_sg")) continue

      end subroutine

      subroutine symmetrize_tensor_sg(sg,t)
!doc$ subroutine symmetrize_tensor(sg,t)
        type(space_group_obj) :: sg
        real(double), dimension(:,:), intent(inout) :: t
!       modifies: t
!       requires: t be dimension(3,3) and be in the lattice representation.
!       effects: Overwrites t with a symmetrized version.
!       errors: Passes errors.

!cod$
        integer :: is, ns
        real(double), dimension(3,3) :: t_sum

        call my(sg)

        if (.not.trivial(sg)) then
          ns = size(sg%o%sym)
          t_sum = 0.0_double
          do is = 1,ns
            t_sum = t_sum + matmul(matmul(sg%o%sym(is)%point_op,t),transpose(sg%o%sym(is)%point_op))
          end do
          t = t_sum/real(ns,double)
        end if

100     call glean(thy(sg))

        if (error("Exit symmetry_mod::symmetrize_tensor_sg")) continue

      end subroutine

      subroutine symmetrize_spherical_tensor_sg(sg,st,mask)
!doc$ subroutine symmetrize_spherical_tensor(sg,st,mask)
        type(space_group_obj) :: sg
        complex(double), dimension(:,:,:) :: st
        logical, dimension(:), intent(in) :: mask
!       requires: size(mask) = size(st,3). l1 and l2 < l_max.
!       requires: st dimensions 1 and 2 be odd and dimension 3 be the number of atoms.
!       modifies: st
!       effects: Overwrites st with a symmetrized version.
!       errors: None.

!cod$
        integer :: ia, iap, is, it, l1, l2, na, ns, nt
        complex(double), dimension(:,:), allocatable :: c_mat, t_mat
        complex(double), dimension(:,:,:), allocatable :: st_sum

        call my(sg)

        if (.not.trivial(sg)) then

          ns = size(sg%o%atom)
          na = size(sg%o%atom(1)%perm,1)
          nt = size(sg%o%atom(1)%perm,2)

          allocate( c_mat(size(st,2),size(st,2)) )
          allocate( t_mat(size(st,1),size(st,1)) )
          allocate( st_sum(size(st,1),size(st,2),size(st,3)) )

          l1 = (size(st,1) - 1)/2
          l2 = (size(st,2) - 1)/2

          st_sum = (0.0_double,0.0_double)
          do is = 1,ns
            c_mat = conjg(sg%o%local(is)%l(l2)%mat)
            t_mat = transpose(sg%o%local(is)%l(l1)%mat)
            do ia = 1,na
              if (.not.mask(ia)) cycle
              do it = 1,nt
                iap = sg%o%atom(is)%perm(ia,it)
                st_sum(:,:,ia) = st_sum(:,:,ia) + matmul(t_mat,matmul(st(:,:,iap),c_mat))
              end do
            end do
          end do
          do ia = 1,na
            if (mask(ia)) st(:,:,ia) = st_sum(:,:,ia)/real(ns*nt,double)
          end do
        end if

100     if (allocated( c_mat)) deallocate( c_mat )
        if (allocated( t_mat)) deallocate( t_mat )
        if (allocated( st_sum)) deallocate( st_sum )

        call glean(thy(sg))

        if (error("Exit symmetry_mod::symmetrize_spherical_tensor_sg")) continue

      end subroutine

      subroutine symmetrize_grid_sg(sg,g)
!doc$ subroutine symmetrize_grid(sg,g)
        type(space_group_obj) :: sg
        type(grid_obj) :: g
!       requires: g be filtered
!       modifies: g
!       effects: Replaces g data with a symmetrized version.
!       errors: Passes errors.

!cod$
        integer :: i, i1, i2, i3, ic, is, n
        real(double) :: scale
        complex(double), dimension(:), allocatable :: c1, c1p
        complex(double), dimension(:,:,:), pointer :: c3

        call my(sg)
        call my(g)

        nullify( c3 )

        if (.not.trivial(sg)) then

          n = size(sg%o%gsmap,2)
          allocate( c1(n) )
          n = size(sg%o%field(1)%perm)
          allocate( c1p(n) )

          n = size(sg%o%field)
          scale = 1.0_double/real(n,double)

          ! Extract g data in Complex Serial Fourier representation
          call take(c3,g,CSF_KIND) ; if (error()) goto 100

          ! Convert non-zero 3-d array data to 1-d array data and scale
          do i = 1,size(c1)
            i1 = sg%o%gsmap(1,i)
            i2 = sg%o%gsmap(2,i)
            i3 = sg%o%gsmap(3,i)
            c1(i) = scale*c3(i1,i2,i3)
          end do

          ! Symmetrize the data
          c1p = (0.0_double,0.0_double)
          do is = 1,size(sg%o%field)
            ic = sg%o%field(is)%phase_class
            do i = 1,size(c1p)
              c1p(i) = c1p(i) + c1(sg%o%field(is)%perm(i))*sg%o%class(ic)%phase(i)
            end do
          end do

          ! Gather the symmetrized data
          call allgatherv(SGROUP,c1p,c1,sg%o%gscount,sg%o%gsdisp)

          ! Convert 1-d array data to 3-d array data
          do i = 1,size(c1)
            i1 = sg%o%gsmap(1,i)
            i2 = sg%o%gsmap(2,i)
            i3 = sg%o%gsmap(3,i)
            c3(i1,i2,i3) = c1(i)
          end do

          ! Insert the 3-d data and convert to Distributed Fourier representation
          call put(c3,g,CSF_KIND)
          call transform(g,CDF_KIND)

        end if

100     if (allocated( c1 )) deallocate( c1 )
        if (allocated( c1p )) deallocate( c1p )
        if (associated( c3 )) deallocate( c3 )

        call glean(thy(sg))
        call glean(thy(g))

        if (error("Exit symmetry_mod::symmetrize_grid_sg")) continue

      end subroutine

      function invariant_site_sg(sg,s) result(l)
!doc$ function invariant_site(sg,s) result(l)
        type(space_group_obj) :: sg
        real(double), dimension(3), intent(in) :: s
        logical :: l
!       effects: Determines if s is invariant wrt sg.

!cod$
        integer :: is, it
        real(double), dimension(3) :: rs, rst

        call my(sg)

        l = .true.

        do is = 1,size(sg%o%sym)
          rs = matmul(sg%o%sym(is)%point_op,s)
          do it = 1,size(sg%o%sym(is)%translation,2)
            rst = rs + sg%o%sym(is)%translation(:,it)
            if (.not.(rst .in. zone(s,tol_equal))) then
              l = .false.
              exit
            end if
          end do
          if (.not.l) exit
        end do

100     call glean(thy(sg))

        if (error("Exit symmetry_mod::invariant_site_sg")) continue

      end function

      subroutine symmetry_unique_atoms_sg(sg,su)
!doc$ subroutine symmetry_unique_atoms(sg,su)
        type(space_group_obj) :: sg
        integer, dimension(:), pointer :: su
!       requires: su dimension be the number of atoms.
!       modifies: su
!       effects: Re-allocates su to the number of symmetry-unique atoms and returns their indices.
!       errors: None.

!cod$
        logical, dimension(:), allocatable :: unique
        integer :: i, ia, iap, is, it, na, ns, nsu, nt

        call my(sg)

        if (trivial(sg)) then

          do ia = 1,size(su)
            su(ia) = ia
          end do

        else

          na = size(sg%o%atom(1)%perm,1)
          ns = size(sg%o%atom)
          nt = size(sg%o%atom(1)%perm,2)

          allocate( unique(na) )
          unique = .true.
          do ia = 1,na
            if (.not.unique(ia)) cycle
            do is = 1,ns
              do it = 1,nt
                iap = sg%o%atom(is)%perm(ia,it)
                if ((iap == ia) .or. .not.unique(iap)) cycle
                unique(iap) = .false.
              end do
            end do
          end do

          nsu = 0
          do ia = 1,na
            if (.not.unique(ia)) cycle
            nsu = nsu + 1
          end do
          deallocate( su )
          allocate( su(nsu) )

          i = 0
          do ia = 1,na
            if (.not.unique(ia)) cycle
            i = i + 1
            su(i) = ia
          end do

        end if

        if (allocated( unique )) deallocate( unique )

        call glean(thy(sg))

        if (error("Exit symmetry_mod::symmetry_unique_atoms_sg")) continue

      end subroutine

      subroutine distribute_spherical_tensor_sg(sg,ia,st,filled)
!doc$ subroutine distribute_spherical_tensor(sg,ia,st,filled)
        type(space_group_obj) :: sg
        integer :: ia
        complex(double), dimension(:,:,:) :: st
        logical, dimension(:) :: filled
!       requires: 1 < ia < number of atoms. filled dimension be the number of atoms. l1 and l2 <= l_max.
!       requires: st dimensions 1 and 2 be odd and dimension 3 be the number of atoms.
!       modifies: st
!       effects: Add contributions to st from non-symmetry-unique atoms.
!       errors: None.

!cod$
        integer :: iap, is, it, l1, l2
        complex(double), dimension(:,:), allocatable :: c_mat, t_mat

        call my(sg)

        if (.not.trivial(sg)) then

          allocate( c_mat(size(st,1),size(st,1)) )
          allocate( t_mat(size(st,2),size(st,2)) )

          l1 = (size(st,1) - 1)/2
          l2 = (size(st,2) - 1)/2

          filled = .false.
          filled(ia) = .true.

          do is = 1,size(sg%o%atom)
            c_mat = conjg(sg%o%local(is)%l(l1)%mat)
            t_mat = transpose(sg%o%local(is)%l(l2)%mat)
            do it = 1,size(sg%o%atom(is)%perm,2)
              iap = sg%o%atom(is)%perm(ia,it)
              if (filled(iap)) cycle
              st(:,:,iap) = matmul(c_mat,matmul(st(:,:,ia),t_mat))
              filled(iap) = .true.
            end do
          end do

        end if

        if (allocated( c_mat )) deallocate( c_mat )
        if (allocated( t_mat )) deallocate( t_mat )

        call glean(thy(sg))

        if (error("Exit symmetry_mod::distribute_spherical_tensor_sg")) continue

      end subroutine

      subroutine distribute_energy_sg(sg,ia,e)
!doc$ subroutine distribute_energy(sg,ia,e)
        type(space_group_obj) :: sg
        integer :: ia
        real(double), dimension(:) :: e
!       requires: 1 < ia < number of atoms. e dimension be the number of atoms.
!       modifies: e
!       effects: Adds contributions to e from non-symmetry-unique atoms.
!       errors: None.

!cod$
        logical, dimension(:), allocatable :: filled
        integer :: iap, is, it

        call my(sg)

        if (.not.trivial(sg)) then

          allocate( filled(size(e)) )
          filled = .false.
          filled(ia) = .true.

          do is = 1,size(sg%o%atom)
            do it = 1,size(sg%o%atom(is)%perm,2)
              iap = sg%o%atom(is)%perm(ia,it)
              if (filled(iap)) cycle
              e(iap) = e(ia)
              filled(iap) = .true.
            end do
          end do

        end if

        if (allocated( filled )) deallocate( filled )

        call glean(thy(sg))

        if (error("Exit symmetry_mod::distribute_energy_sg")) continue

      end subroutine

      subroutine diary_sg(sg,list,lat)
!doc$ subroutine diary(sg,list,lat)
        type(space_group_obj) :: sg
        logical, intent(in), optional :: list
        type(lattice_obj), optional :: lat
!       modifies: diary stream
!       effects: Diaries sg. If list is present and true, sg operations are listed.
!       errors: Passes errors.

!cod$
        character(4) :: pg_type
        logical :: fnd, list_local
        integer :: is, it, ns, nt
        integer, dimension(10) :: op_type
        real(double), dimension(3) :: tv
        real(double), dimension(3,3) :: tm

        call my(sg)
        if (present(lat)) call my(lat)

        if (present(list)) then
          list_local = list
        else
          call arg("list_space_group",list_local,fnd)
          if (.not.fnd) list_local = .false.
        end if

        ns = size(sg%o%sym)
        nt = size(sg%o%sym(1)%translation,2)

        op_type = 0
        do is = 1,ns
           op_type = op_type + classify_op_i(sg%o%sym(is)%point_op)
        end do
        call get_point_group_type_i(op_type,pg_type) ; if (error()) goto 100
        if (i_access(diaryfile())) then
          select case (sg%o%mode)
          case (FULL)
            write(x_unit(diaryfile()),'(/t4,"Crystal symmetry (full updates with tolerance =",es8.1,"):")') sg%o%tol_aperm
          case (AUTO)
            write(x_unit(diaryfile()),'(/t4,"Crystal symmetry (automatic with tolerance =",es8.1,"):")') sg%o%tol_aperm
          case (OFF)
            write(x_unit(diaryfile()),'(/t4,"Crystal symmetry (off):")')
          end select
          if (ns*nt == 1) then
            write(x_unit(diaryfile()),'(/,t6,"The space group has ",i0," operation with point group ",a)') ns*nt, trim(pg_type)
          else
            write(x_unit(diaryfile()),'(/,t6,"The space group has ",i0," operations with point group ",a)') ns*nt, trim(pg_type)
          end if
          if (list_local) then
            if (present(lat)) then
              write(x_unit(diaryfile()),'(/,t6,"Space-group operations in the cartesian representation:")')
              do is = 1,ns
                tm = lat2r(lat,sg%o%sym(is)%point_op)
                write(x_unit(diaryfile()),'(/,t16,3f6.1)') tm(1,:)
                write(x_unit(diaryfile()),'(t6,i2,".",7x,3f6.1,7x,sp,3f15.10)') is, tm(2,:)
                write(x_unit(diaryfile()),'(  t16,3f6.1)') tm(3,:)
                do it = 1,nt
                  tv = lat2r(lat,sg%o%sym(is)%translation(:,it))
                  write(x_unit(diaryfile()),'(t34,3f15.10)') tv
                end do
              end do
            else
              write(x_unit(diaryfile()),'(/,t6,"Space-group operations in the lattice representation:")')
              do is = 1,ns
                write(x_unit(diaryfile()),'(/,t16,3f6.1)') sg%o%sym(is)%point_op(1,:)
                write(x_unit(diaryfile()),'(t6,i2,".",7x,3f6.1,7x,sp,3f15.10)') is, sg%o%sym(is)%point_op(2,:)
                write(x_unit(diaryfile()),'(  t16,3f6.1)') sg%o%sym(is)%point_op(3,:)
                do it = 1,nt
                  write(x_unit(diaryfile()),'(t34,3f15.10)') sg%o%sym(is)%translation(:,it)
                end do
              end do
            end if
          end if
        end if

100     call glean(thy(sg))
        if (present(lat)) call glean(thy(lat))

        if (error("Exit symmetry_mod::diary_sg")) continue

      end subroutine

      subroutine diary_atom_stars_sg(sg,at)
        type(space_group_obj) :: sg
        type(atoms_obj) :: at

        integer :: ia, is, it, na, ns, nt
        real(double), dimension(3) :: p, rp

        call my(sg)
        call my(at)

        na = x_n_atoms(at)
        ns = size(sg%o%sym)

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,t6,"Atom coordinates:",/)')
          do ia = 1,na
            p = x_position(at,ia)
            write(x_unit(diaryfile()),'(t8,i0,".",4x,3f20.15)') ia, p
          end do
          write(x_unit(diaryfile()),'(/,t6,"Rotated/shifted atom coordinates:")')
          do ia = 1,na
            write(x_unit(diaryfile()),'(/,t8,"Atom #",i0,":")') ia
            p = x_position(at,ia)
            if (.not.trivial(sg)) then
              do is = 1,ns
                rp = matmul(sg%o%sym(is)%point_op,p)
                nt = size(sg%o%sym(is)%translation,2)
                do it = 1,nt
                  write(x_unit(diaryfile()),'(t10,3f20.15,";",2i4)') (rp + sg%o%sym(is)%translation(:,it)), is, it
                end do
              end do
            else
              write(x_unit(diaryfile()),'(t10,3f20.15)') p
            end if
          end do
        end if

        call barrier(FILE_SCOPE)

        call glean(thy(sg))
        call glean(thy(at))

      end subroutine

      subroutine save_sg(sg)
!doc$ subroutine save(sg)
        type(space_group_obj) :: sg
!       modifies: file system
!       effects: Saves the sg symmetry operations.
!       errors: File I/O problems.

!cod$
        integer :: ios, is, it, ns, nt
        type(file_obj) :: f

        call my(sg)

        call my(file(trim(new_space_group_path)),f)

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100

        ns = size(sg%o%sym)
        nt = size(sg%o%sym(1)%translation,2)
        if (i_access(f)) then
          write(x_unit(f),'(i0,2x,i0)') ns, nt
          do is = 1,ns
            write(x_unit(f),'(a)') " "
            write(x_unit(f),'(3i5)') nint(sg%o%sym(is)%point_op(1,:))
            write(x_unit(f),'(3i5)') nint(sg%o%sym(is)%point_op(2,:))
            write(x_unit(f),'(3i5)') nint(sg%o%sym(is)%point_op(3,:))
            do it = 1,nt
              write(x_unit(f),'(3f16.10)') sg%o%sym(is)%translation(:,it)
            end do
          end do
        end if

        if (i_access(f)) close(x_unit(f))

100     call glean(thy(sg))
        call glean(thy(f))

        if (error("Exit symmetry_mod::save_sg")) continue

      end subroutine

! point group local routines

      subroutine own_pg_i(pg)
        type(point_group_obj) :: pg
        type(point_group_obj) :: pgt
        integer :: is
        if (pg%ref < pg%o%ref) then
          allocate( pgt%o )
          pgt%o%ref = 0
          pgt%o%g = pg%o%g
          pgt%o%save = pg%o%save
          pgt%o%mode = pg%o%mode
          pgt%o%source = pg%o%source
          pgt%o%g_source = pg%o%g_source
          allocate( pgt%o%sym(size(pg%o%sym)) )
          do is = 1,size(pg%o%sym)
            pgt%o%sym(is)%point_op = pg%o%sym(is)%point_op
          end do
          pg%o%ref = pg%o%ref - pg%ref
          pg%o => pgt%o
          pg%o%ref = pg%o%ref + pg%ref
        end if
      end subroutine

      subroutine generate_lattice_group_i(pgr,lat)
        type(point_group_rep) :: pgr
        type(lattice_obj) :: lat

        real(double), parameter :: tol_vector = 1.0e-8_double
        real(double), parameter :: tol_matrix = 2.0e-8_double

        integer :: i1, i2, i3, is, na, ns
        real(double) :: vmag
        real(double), dimension(3) :: mag
        real(double), dimension(3,3) :: mat, test_mat, unitary_form
        real(double), dimension(:,:), allocatable :: a
        real(double), dimension(:,:,:), allocatable :: tpo

        call my(lat)

        if (associated( pgr%sym )) deallocate( pgr%sym )

        pgr%g_source = x_ghost(lat)

        allocate( a(3,27), tpo(3,3,48) )

        mag(1) = norm(lat2r(lat,real((/1,0,0/),double)))  ! MAKE A LIST OF VECTORS WHOSE MAGNITUDES
        mag(2) = norm(lat2r(lat,real((/0,1,0/),double)))  ! EQUAL ONE OF THE THREE LATTICE VECTORS.
        mag(3) = norm(lat2r(lat,real((/0,0,1/),double)))
        na = 0
        do i3 = -1,+1
        do i2 = -1,+1
        do i1 = -1,+1
          vmag = norm(lat2r(lat,real((/i1,i2,i3/),double)))
          if ( any(abs(mag - vmag) < tol_vector) ) then
            na = na + 1
            a(:,na) = real((/i1,i2,i3/),double)
          end if
        end do
        end do
        end do

        mat(1,:) = lat2r(lat,real((/1,0,0/),double))  ! FORM CANDIDATE TRANSFORMATION MATRICES FROM TRIPLETS OF
        mat(2,:) = lat2r(lat,real((/0,1,0/),double))  ! VECTORS FROM THE LIST. TEST THAT THE TRANSFORMATIONS ARE
        mat(3,:) = lat2r(lat,real((/0,0,1/),double))  ! UNITARY AND SAVE.
        unitary_form = matmul(mat,transpose(mat))
        ns = 0
        do i1 = 1,na
        do i2 = 1,na
        do i3 = 1,na
          mat = reshape((/a(:,i1),a(:,i2),a(:,i3)/),(/3,3/))
          test_mat = matmul(matmul(transpose(mat),unitary_form),mat)
          if (test_mat .in. nbhd(unitary_form,tol_matrix)) then
            ns = ns + 1
            if (error(ns > 48,"ERROR: more than 48 point group operations")) goto 100
            tpo(:,:,ns) = mat
          end if
        end do
        end do
        end do

        allocate( pgr%sym(ns) )
        do is = 1,ns
          pgr%sym(is)%point_op = tpo(:,:,is)
        end do

100     if (allocated( a )) deallocate( a )
        if (allocated( tpo )) deallocate( tpo )

        call glean(thy(lat))

        if (error("Exit symmetry_mod::generate_lattice_group_i")) continue

      end subroutine

      subroutine read_lattice_group_i(pgr,lat)
        type(point_group_rep) :: pgr
        type(lattice_obj) :: lat

        logical :: exist_file
        character(line_len) :: dummy
        integer :: ios, is, ns
        integer, dimension(:,:,:), allocatable :: tpo
        type(file_obj) :: f

        call my(lat)
        call my(file(trim(lattice_group_path)),f)

        if (associated( pgr%sym )) deallocate( pgr%sym )

        pgr%g_source = x_ghost(lat)

        if (i_access(f)) inquire(file=x_name(f),exist=exist_file)
        if (i_comm(f)) call broadcast(FILE_SCOPE,exist_file)
        if (error(.not.exist_file,"ERROR: lattice group file was not found")) goto 200

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 200

        if (i_access(f)) read(x_unit(f),*,iostat=ios) ns
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to read the number of symmetry operations")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,ns)

        allocate( tpo(3,3,ns) )
        if (i_access(f)) then
          do is = 1,ns
            read(x_unit(f),'(a)') dummy
            read(x_unit(f),*) tpo(1,:,is)
            read(x_unit(f),*) tpo(2,:,is)
            read(x_unit(f),*) tpo(3,:,is)
          end do
        end if
        if (i_comm(f)) call broadcast(FILE_SCOPE,tpo)

        allocate( pgr%sym(ns) )
        do is = 1,ns
          pgr%sym(is)%point_op = real(tpo(:,:,is),double)
        end do

        if (error(.not.valid_lattice_group_i(pgr,lat),"ERROR: Lattice group is not valid")) goto 100

100     if (i_access(f)) close(x_unit(f))

200     if (allocated( tpo )) deallocate( tpo )

        call glean(thy(lat))
        call glean(thy(f))

        if (error("Exit symmetry_mod::read_lattice_group_i")) continue

      end subroutine

      function valid_lattice_group_i(pgr,lat) result(vlg)
        type(point_group_rep) :: pgr
        type(lattice_obj) :: lat
        logical :: vlg

        integer :: is, ns
        real(double), dimension(3,3) :: mat, test_mat, unitary_form

        call my(lat)

        mat(1,:) = lat2r(lat,real((/1,0,0/),double))
        mat(2,:) = lat2r(lat,real((/0,1,0/),double))
        mat(3,:) = lat2r(lat,real((/0,0,1/),double))
        unitary_form = matmul(mat,transpose(mat))

        ns = size(pgr%sym)

        vlg = .true.
        do is = 1,ns
          mat = pgr%sym(is)%point_op
          test_mat = matmul(matmul(transpose(mat),unitary_form),mat)
          if (.not.(test_mat .in. nbhd(unitary_form,tol_equal))) then
            vlg = .false.
            exit
          end if
        end do

        call glean(thy(lat))

      end function

      function same_point_group_i(pgr1,pgr2) result(spg)
        type(point_group_rep) :: pgr1, pgr2
        logical :: spg

        integer :: is, ns

        spg = .true.

        if (size(pgr1%sym) /= size(pgr2%sym)) then
          spg = .false.
          goto 100
        end if

        if (size(pgr1%sym) == 1) goto 100

        ns = size(pgr1%sym)

        do is = 1,ns
          spg = spg .and. (pgr1%sym(is)%point_op .in. nbhd(pgr2%sym(is)%point_op,tol_equal))
          if (.not.spg) exit
        end do

100     continue

      end function

      function compare_point_ops_i(pg1,pg2) result(spo)
        type(point_group_obj) :: pg1, pg2
        logical :: spo
        integer :: is1, is2, ns1, ns2
        call my(pg1)
        call my(pg2)
        ns1 = size(pg1%o%sym)
        ns2 = size(pg2%o%sym)
        if (ns1 /= ns2) then
          spo = .false.
        else
          do is1 = 1,ns1
            do is2 = 1,ns2
              spo = ( pg1%o%sym(is1)%point_op .in. nbhd(pg2%o%sym(is2)%point_op,tol_equal) )
              if (spo) exit
            end do
            if (.not.spo) exit
          end do
        end if
        call glean(thy(pg1))
        call glean(thy(pg2))
      end function

      subroutine save_pg_i(pgr)
        type(point_group_rep) :: pgr

        integer :: ios, is, ns
        type(file_obj) :: f

        select case (pgr%source)
        case (LATT)
          call my(file(trim(new_lattice_group_path)),f)
        case (SGRP)
          call my(file(trim(new_point_group_path)),f)
        end select

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100

        ns = size(pgr%sym)
        if (i_access(f)) then
          write(x_unit(f),'(i0)') ns
          do is = 1,ns
            write(x_unit(f),'(a)') " "
            write(x_unit(f),'(3i5)') nint(pgr%sym(is)%point_op(1,:))
            write(x_unit(f),'(3i5)') nint(pgr%sym(is)%point_op(2,:))
            write(x_unit(f),'(3i5)') nint(pgr%sym(is)%point_op(3,:))
          end do
        end if

        if (i_access(f)) close(x_unit(f))

100     call glean(thy(f))

        if (error("Exit symmetry_mod::save_pg_i")) continue

      end subroutine

      subroutine extract_point_group_i(pgr,sg,parity)
        type(point_group_rep) :: pgr
        type(space_group_obj) :: sg
        logical, intent(in) :: parity

        logical :: par
        integer :: is, ns

        call my(sg)

        if (associated( pgr%sym )) deallocate( pgr%sym )

        pgr%g_source = x_ghost(sg)

        ns = size(sg%o%sym)
        if (parity) then
          par = parity_space_group_i(sg)
          if (par) then
            allocate( pgr%sym(ns) )
            do is = 1,ns
              pgr%sym(is)%point_op = sg%o%sym(is)%point_op
            end do
          else
            allocate( pgr%sym(2*ns) )
            do is = 1,ns
              pgr%sym(is)%point_op = sg%o%sym(is)%point_op
              pgr%sym(is+ns)%point_op = -sg%o%sym(is)%point_op
            end do
          end if
        else
          allocate( pgr%sym(ns) )
          do is = 1,ns
            pgr%sym(is)%point_op = sg%o%sym(is)%point_op
          end do
        end if

        call glean(thy(sg))

      end subroutine

! space group local routines

      function parity_space_group_i(sg) result(psg)
        type(space_group_obj) :: sg
        logical :: psg

        integer :: is, ns

        call my(sg)

        ns = size(sg%o%sym)
        psg = .false.
        do is = 1,ns
          psg = psg .or. (trace(sg%o%sym(is)%point_op) .in. nbhd(-3.0_double,tol_equal))
          if (psg) exit
        end do

        call glean(thy(sg))

      end function

      subroutine generate_space_group_i(sgr,pg,ats,lat,lay)
        type(space_group_rep) :: sgr
        type(point_group_obj) :: pg
        type(atoms_obj) :: ats
        type(lattice_obj) :: lat
        type(layout_obj), optional :: lay

        logical :: identity, mapped, mapped_local, new_type
        character(4) :: pg_type
        character(tag_sz) :: type_min
        character(tag_sz), dimension(:), allocatable :: temp_type
        integer :: i1, i2, i3, ia, ip, ipt, is, it, l, na, ns, ns_pg, nt, slot
        integer :: base_num, first_a, last_a, na_proc
        integer, dimension(1) :: min_slot
        integer, dimension(10) :: op_type
        integer, dimension(:), allocatable :: map, map_proc, count, disp, type_count
        integer, dimension(:,:), allocatable :: temp_map, trans_map
        real(double), dimension(3) :: v, trans
        real(double), dimension(:,:), allocatable :: atom_rot_base, pure_trans, temp_tau
        real(double), dimension(:,:,:), allocatable :: atom_rot, temp_sym

        call my(pg)
        call my(ats)
        call my(lat)
        if (present(lay)) call my(lay)

        sgr%g_point_group = x_ghost(pg)
        sgr%g_atoms = x_ghost(ats)
        if (present(lay)) then
          sgr%g_layout = x_ghost(lay)
        else
          sgr%g_layout = x_ghost()
        end if

        if (associated( sgr%sym )) then
          do is = 1,size(sgr%sym)
            if (associated( sgr%sym(is)%translation )) deallocate( sgr%sym(is)%translation )
          end do
          deallocate( sgr%sym )
        end if
        if (associated( sgr%atom )) then
          do is = 1,size(sgr%atom)
            if (associated( sgr%atom(is)%perm )) deallocate( sgr%atom(is)%perm )
          end do
          deallocate( sgr%atom )
        end if
        if (associated( sgr%field )) then
          do is = 1,size(sgr%field)
            if (associated( sgr%field(is)%perm )) deallocate( sgr%field(is)%perm )
          end do
          deallocate( sgr%field )
        end if
        if (associated( sgr%class )) then
          do is = 1,size(sgr%class)
            if (associated( sgr%class(is)%phase )) deallocate( sgr%class(is)%phase )
          end do
          deallocate( sgr%class )
        end if
        if (associated( sgr%local )) then
          do is = 1,size(sgr%local)
            if (associated( sgr%local(is)%l )) then
              do l = 0,size(sgr%local(is)%l)-1
                if (associated( sgr%local(is)%l(l)%mat )) deallocate( sgr%local(is)%l(l)%mat )
              end do
              deallocate( sgr%local(is)%l )
            end if
          end do
          deallocate( sgr%local )
        end if
        if (associated( sgr%gscount )) deallocate( sgr%gscount )
        if (associated( sgr%gsdisp )) deallocate( sgr%gsdisp )

        na = x_n_atoms(ats)
        ns_pg = size(pg%o%sym)

        if (error(na < 1,"ERROR: atoms must have at least one atom in it")) goto 100

        allocate( type_count(na), temp_type(na) )                              ! SELECT A BASE ATOM FOR DETERMINING TRANSLATIONS
        type_count = na
        nt = 0
        do ia = 1,na
          new_type = .true.
          do it = 1,nt
            if (x_type(ats,ia) == temp_type(it)) then
              new_type = .false.
              slot = it
              exit
            end if
          end do
          if (new_type) then
            nt = nt + 1
            type_count(nt) = 1
            temp_type(nt) = x_type(ats,ia)
          else
            type_count(slot) = type_count(slot) + 1
          end if
        end do
        min_slot = minloc(type_count)
        type_min = temp_type(min_slot(1))
        deallocate( temp_type, type_count )
        do ia = 1,na
          if (x_type(ats,ia) == type_min) then
            base_num = ia
            exit
          end if
        end do

        call subdivide(mpi_myproc(CONFIG),mpi_nprocs(CONFIG),1,na,first_a,last_a,na_proc) ! DIVIDE THE ATOMS AMONG PROCESSORS

        allocate( temp_sym(3,3,ns_pg), temp_tau(3,ns_pg) )
        allocate( temp_map(na,ns_pg), pure_trans(3,na), trans_map(na,na) )
        allocate( atom_rot(3,na_proc,ns_pg), atom_rot_base(3,ns_pg) )
        allocate( map(na), map_proc(na_proc) )
        allocate( count(0:mpi_nprocs(CONFIG)-1), disp(0:mpi_nprocs(CONFIG)-1) )

        call allgather(CONFIG,na_proc,count)                                              ! SET UP THE COMMUNICATION ARRAYS
        disp = 0
        do ip = 0,(mpi_nprocs(CONFIG)-1)
          do ipt = 0,(ip - 1)
            disp(ip) = disp(ip) + count(ipt)
          end do
        end do

        do is = 1,ns_pg                                                                   ! PRECOMPUTE TRANSFORMED ATOM POSITIONS
          atom_rot_base(:,is) = matmul(pg%o%sym(is)%point_op,x_position(ats,base_num))
          do ia = 1,na_proc
            atom_rot(:,ia,is) = matmul(pg%o%sym(is)%point_op,x_position(ats,ia+first_a-1))
          end do
        end do

        ns = 0
        nt = 0
        mapped_local = .true.
        do is = 1,ns_pg                                                                     ! IDENTIFY SPACE-GROUP OPERATIONS
          identity = (trace(pg%o%sym(is)%point_op) .in. nbhd(3.0_double,tol_identity))
          do i1 = 1,na
            if (x_type(ats,i1) /= x_type(ats,base_num)) cycle
            trans = x_position(ats,i1) - atom_rot_base(:,is)                                ! CANDIDATE TRANSLATION
            do i2 = 1,na_proc
              v = atom_rot(:,i2,is) + trans                                                 ! CANDIDATE ROTATION + TRANSLATION
              mapped_local = .false.
              do i3 = 1,na
                if (x_type(ats,i3) == x_type(ats,i2+first_a-1)) then                        ! CHECK FOR MAPPING
                  if (v .in. zone(x_position(ats,i3),sgr%tol_aperm)) then
                    mapped_local = .true.
                    map_proc(i2) = i3
                    exit
                  end if
                end if
              end do
              if (.not. mapped_local) exit
            end do
            call allreduce(CONFIG,MPI_LAND,mapped_local,mapped)
            call allgatherv(CONFIG,map_proc,map,count,disp)
            if (mapped) then                                                                ! SAVE THE VALIDATED OPERATION
              if (identity) then
                nt = nt + 1
                pure_trans(:,nt) = modulo(trans,1.0_double)
                trans_map(:,nt) = map
                if (nt == 1) then
                  ns = ns + 1
                  temp_sym(:,:,ns) = pg%o%sym(is)%point_op
                  temp_tau(:,ns) = modulo(trans,1.0_double)
                  temp_map(:,ns) = map
                end if
              else
                ns = ns + 1
                temp_sym(:,:,ns) = pg%o%sym(is)%point_op
                temp_tau(:,ns) = modulo(trans,1.0_double)
                temp_map(:,ns) = map
                exit
              end if
            end if
          end do
        end do

        if (error(ns == 0,"ERROR: no point operations were found")) goto 100
        if (error(nt == 0,"ERROR: no translations were found")) goto 100

        call sort_by_magnitude_i(lat,nt,pure_trans,trans_map)                              ! SORT THE PURE TRANSLATIONS
        call group_by_translation_i(lat,ns,temp_sym,temp_tau,temp_map)                     ! GROUP THE SPACE-GROUP OPERATIONS

        allocate( sgr%sym(ns) )                                                            ! STORE THE SPACE-GROUP OPERATIONS
        do is = 1,ns
          sgr%sym(is)%point_op = temp_sym(:,:,is)
          allocate( sgr%sym(is)%translation(3,nt) )
          do it = 1,nt
            sgr%sym(is)%translation(:,it) = temp_tau(:,is) + pure_trans(:,it)
          end do
        end do

        op_type = 0                                                                        ! CHECK THE POINT GROUP
        do is = 1,ns
           op_type = op_type + classify_op_i(sgr%sym(is)%point_op)
        end do
        call get_point_group_type_i(op_type,pg_type)
        if (error()) then
          if (i_access(diaryfile())) then
            write(x_unit(diaryfile()),'(/,t6,"Space-group point operations in the lattice representation:")')
            do is = 1,ns
              write(x_unit(diaryfile()),'(/,t16,3f6.1)') sgr%sym(is)%point_op(1,:)
              write(x_unit(diaryfile()),'(t6,i2,".",7x,3f6.1)') is, sgr%sym(is)%point_op(2,:)
              write(x_unit(diaryfile()),'(t16,3f6.1)') sgr%sym(is)%point_op(3,:)
            end do
          end if
          goto 100
        end if

        if (present(lay)) then
          if (ns*nt > 1) then
            call form_symmetrizing_structures_i(sgr,temp_tau,temp_map,trans_map,lat,lay) ; if (error()) goto 100
          end if
          if (ns > 1) then
            call check_mesh_dimensions_i(sgr,lay) ; if (error()) goto 100
          end if
        else
          allocate( sgr%atom(ns) )
          call form_atom_permutations_i(sgr,temp_map,trans_map) ; if (error()) goto 100
        end if

100     if (allocated( pure_trans )) deallocate( pure_trans )
        if (allocated( trans_map )) deallocate( trans_map )
        if (allocated( temp_sym )) deallocate( temp_sym )
        if (allocated( temp_tau )) deallocate( temp_tau )
        if (allocated( temp_map )) deallocate( temp_map )
        if (allocated( temp_type )) deallocate( temp_type )
        if (allocated( type_count )) deallocate( type_count )
        if (allocated( count )) deallocate( count )
        if (allocated( disp )) deallocate( disp )
        if (allocated( map )) deallocate( map )
        if (allocated( map_proc )) deallocate( map_proc )
        if (allocated( atom_rot )) deallocate( atom_rot )
        if (allocated( atom_rot_base )) deallocate( atom_rot_base )

        call glean(thy(pg))
        call glean(thy(ats))
        call glean(thy(lat))
        if (present(lay)) call glean(thy(lay))

        if (error("Exit symmetry_mod::generate_space_group_i")) continue

      end subroutine

      subroutine form_trivial_space_group_i(sgr)
        type(space_group_rep) :: sgr

        integer :: is, l

        if (associated( sgr%sym )) then
          do is = 1,size(sgr%sym)
            if (associated( sgr%sym(is)%translation )) deallocate( sgr%sym(is)%translation )
          end do
          deallocate( sgr%sym )
        end if
        if (associated( sgr%atom )) then
          do is = 1,size(sgr%atom)
            if (associated( sgr%atom(is)%perm )) deallocate( sgr%atom(is)%perm )
          end do
          deallocate( sgr%atom )
        end if
        if (associated( sgr%field )) then
          do is = 1,size(sgr%field)
            if (associated( sgr%field(is)%perm )) deallocate( sgr%field(is)%perm )
          end do
          deallocate( sgr%field )
        end if
        if (associated( sgr%class )) then
          do is = 1,size(sgr%class)
            if (associated( sgr%class(is)%phase )) deallocate( sgr%class(is)%phase )
          end do
          deallocate( sgr%class )
        end if
        if (associated( sgr%local )) then
          do is = 1,size(sgr%local)
            if (associated( sgr%local(is)%l )) then
              do l = 0,size(sgr%local(is)%l)-1
                if (associated( sgr%local(is)%l(l)%mat )) deallocate( sgr%local(is)%l(l)%mat )
              end do
              deallocate( sgr%local(is)%l )
            end if
          end do
          deallocate( sgr%local )
        end if
        if (associated( sgr%gsmap )) deallocate( sgr%gsmap )
        if (associated( sgr%gscount )) deallocate( sgr%gscount )
        if (associated( sgr%gsdisp )) deallocate( sgr%gsdisp )

        allocate( sgr%sym(1) )
        sgr%sym(1)%point_op = reshape(real((/1,0,0,0,1,0,0,0,1/),double),(/3,3/))
        allocate( sgr%sym(1)%translation(3,1) )
        sgr%sym(1)%translation(:,1) = (/0.0_double,0.0_double,0.0_double/)

      end subroutine

      function valid_space_group_i(sgr,at) result(vsg)
        type(space_group_rep) :: sgr
        type(atoms_obj) :: at
        logical :: vsg

        integer :: ia, iap, ii, is, it, na, ns, nt
        real(double), dimension(3) :: p, pt

        call my(at)

        vsg = .true.

        if ((size(sgr%sym) == 1) .and. (size(sgr%sym(1)%translation,2) == 1)) goto 100

        ns = size(sgr%atom)
        na = size(sgr%atom(1)%perm,1)
        nt = size(sgr%atom(1)%perm,2)

        if (x_n_atoms(at) /= na) then
          vsg = .false.
          goto 100
        end if

        if (nt > 1) then                                                             ! CHECK THE PURE TRANSLATIONS
          do is = 1,ns
            if (trace(sgr%sym(is)%point_op) .in. nbhd(3.0_double,tol_identity)) then
              ii = is
              exit
            end if
          end do
          pure_trans: do ia = 1,na
            p = x_position(at,ia)
            do it = 1,nt
              iap = sgr%atom(ii)%perm(ia,it)
              if (x_type(at,iap) /= x_type(at,ia)) then
                vsg = .false.
                exit pure_trans
              end if
              pt = x_position(at,ia) + sgr%sym(ii)%translation(:,it)
              if (.not.(pt .in. zone(x_position(at,iap),sgr%tol_aperm))) then
                vsg = .false.
                exit pure_trans
              end if
            end do
          end do pure_trans
        end if
        if (.not.vsg) goto 100

        prim_ops: do ia = 1,na                                                       ! CHECK THE PRIMITIVE OPERATIONS
          p = x_position(at,ia)
          do is = 1,ns
            iap = sgr%atom(is)%perm(ia,1)
            if (x_type(at,iap) /= x_type(at,ia)) then
              vsg = .false.
              exit prim_ops
            end if
            pt = matmul(sgr%sym(is)%point_op,p) + sgr%sym(is)%translation(:,1)
            if (.not.(pt .in. zone(x_position(at,iap),sgr%tol_aperm))) then
              vsg = .false.
              exit prim_ops
            end if
          end do
        end do prim_ops

100     call glean(thy(at))

      end function

      function same_space_group_i(sgr1,sgr2) result(ssg)
        type(space_group_rep) :: sgr1, sgr2
        logical :: ssg

        integer :: is, it, ns, nt

        ssg = .true.

        if (size(sgr1%sym) /= size(sgr2%sym)) then
          ssg = .false.
          goto 100
        end if

        if (size(sgr1%sym(1)%translation,2) /= size(sgr2%sym(1)%translation,2)) then
          ssg = .false.
          goto 100
        end if

        if ((size(sgr1%sym) == 1) .and. (size(sgr1%sym(1)%translation,2) == 1)) goto 100

        if (size(sgr1%atom(1)%perm,1) /= size(sgr2%atom(1)%perm,1)) then
          ssg = .false.
          goto 100
        end if

        ns = size(sgr1%sym)
        nt = size(sgr1%sym(1)%translation,2)

        outer_loop: do is = 1,ns
          ssg = ssg .and. (sgr1%sym(is)%point_op .in. nbhd(sgr2%sym(is)%point_op,tol_equal))
          if (.not.ssg) exit
          do it = 1,nt
            ssg = ssg .and. (sgr1%sym(is)%translation(:,it) .in. zone(sgr2%sym(is)%translation(:,it),tol_equal))
            if (.not.ssg) exit outer_loop
          end do
          ssg = ssg .and. all(sgr1%atom(is)%perm == sgr2%atom(is)%perm)
          if (.not.ssg) exit
        end do outer_loop

100     continue

      end function

      function compare_space_ops_i(sg1,sg2) result(sso)
        type(space_group_obj) :: sg1, sg2
        logical :: sso
        logical :: spo, st
        integer :: is1, is2, it1, it2, ns1, ns2, nt1, nt2
        call my(sg1)
        call my(sg2)
        ns1 = size(sg1%o%sym)
        ns2 = size(sg2%o%sym)
        nt1 = size(sg1%o%sym(1)%translation,2)
        nt2 = size(sg2%o%sym(1)%translation,2)
        if ((ns1 /= ns2) .or. (nt1 /= nt2)) then
          sso = .false.
        else
          do is1 = 1,ns1
            do is2 = 1,ns2
              spo = ( sg1%o%sym(is1)%point_op .in. nbhd(sg2%o%sym(is2)%point_op,tol_equal) )
              if (spo) exit
            end do
            if (spo) then
              do it1 = 1,nt1
                do it2 = 1,nt2
                  st = ( sg1%o%sym(is1)%translation(:,it1) .in. zone(sg2%o%sym(is2)%translation(:,it2),tol_equal) )
                  if (st) exit
                end do
                if (.not.st) exit
              end do
              sso = (spo .and. st)
            else
              sso = .false.
            end if
            if (.not.sso) exit
          end do
        end if
        call glean(thy(sg1))
        call glean(thy(sg2))
      end function

      function compare_atom_perms_i(sg1,sg2) result(sap)
        type(space_group_obj) :: sg1, sg2
        logical :: sap
!       requires: Same numbers of point operations, atoms and translations. Same point operations.
        logical :: spo
        integer :: is1, is2, na1, na2, ns, ns1, ns2, nt1, nt2
        call my(sg1)
        call my(sg2)
        ns1 = size(sg1%o%atom) ; ns2 = size(sg2%o%atom)
        if (error(ns1 /= ns2,"ERROR: different numbers of point operations")) goto 100
        ns = size(sg1%o%atom)
        na1 = size(sg1%o%atom(1)%perm,1) ; na2 = size(sg2%o%atom(1)%perm,1)
        if (error(na1 /= na2,"ERROR: different numbers of atoms")) goto 100
        nt1 = size(sg1%o%atom(1)%perm,2) ; nt2 = size(sg2%o%atom(1)%perm,2)
        if (error(nt1 /= nt2,"ERROR: different numbers of translations")) goto 100
        do is1 = 1,ns
          do is2 = 1,ns
            spo = ( sg1%o%sym(is1)%point_op .in. nbhd(sg2%o%sym(is2)%point_op,tol_equal) )
            if (spo) then
              sap = all(sg1%o%atom(is1)%perm == sg2%o%atom(is2)%perm)
              exit
            end if
          end do
          if (error(.not.spo,"ERROR: cannot match point operations")) goto 100
          if (.not.sap) exit
        end do
100     call glean(thy(sg1))
        call glean(thy(sg2))
        if (error("Exit symmetry_mod::compare_atom_perms_i")) continue
      end function

      function compare_field_ops_i(sg1,sg2) result(sfo)
        type(space_group_obj) :: sg1, sg2
        logical :: sfo
!       requires: Same point operations.
        logical :: spo
        integer :: is1, i_ps, is2, n_ps, n_ps1, n_ps2, ns, ns1, ns2, pc1, pc2
        real(double), parameter :: tol_equal_phase = 1.0e-10_double
        call my(sg1)
        call my(sg2)
        ns1 = size(sg1%o%field) ; ns2 = size(sg2%o%field)
        if (error(ns1 /= ns2,"ERROR: different numbers of point operations")) goto 100
        ns = ns1
        n_ps1 = size(sg1%o%field(1)%perm) ; n_ps2 = size(sg2%o%field(1)%perm)
        if (n_ps1 /= n_ps2) then
          sfo = .false.
          goto 100
        end if
        n_ps = n_ps1
        do is1 = 1,ns
          do is2 = 1,ns
            spo = ( sg1%o%sym(is1)%point_op .in. nbhd(sg2%o%sym(is2)%point_op,tol_equal) )
            if (spo) then
              sfo = all(sg1%o%field(is1)%perm == sg2%o%field(is2)%perm)
              if (.not.sfo) exit
              pc1 = sg1%o%field(is1)%phase_class
              pc2 = sg2%o%field(is2)%phase_class
              do i_ps = 1,n_ps
                sfo = ( sg1%o%class(pc1)%phase(i_ps) .in. nbhd(sg2%o%class(pc2)%phase(i_ps),tol_equal_phase) )
                if (.not.sfo) exit
              end do
              exit
            end if
          end do
          if (error(.not.spo,"ERROR: cannot match point operations")) goto 100
          if (.not.sfo) exit
        end do
        sfo = sfo .and. ( size(sg1%o%class) /= size(sg2%o%class) )
100     call glean(thy(sg1))
        call glean(thy(sg2))
        if (error("Exit symmetry_mod::compare_field_ops_i")) continue
      end function

      subroutine own_sg_i(sg)
        type(space_group_obj) :: sg
        type(space_group_obj) :: sgt
        integer :: ic, is, l, na, nc, nl, np, ns, nt
        if (sg%ref < sg%o%ref) then
          allocate( sgt%o )
          sgt%o%ref = 0
          sgt%o%g = sg%o%g
          sgt%o%save = sg%o%save
          sgt%o%mode = sg%o%mode
          select case (sg%o%mode)
          case (FULL,AUTO)
            sgt%o%g_point_group = sg%o%g_point_group
            sgt%o%g_atoms = sg%o%g_atoms
            sgt%o%g_layout = sg%o%g_layout
          end select
          sgt%o%tol_aperm = sg%o%tol_aperm
          ns = size(sg%o%sym)
          nt = size(sg%o%sym(1)%translation,2)
          allocate( sgt%o%sym(ns) )
          do is = 1,ns
            sgt%o%sym(is)%point_op = sg%o%sym(is)%point_op
            allocate( sgt%o%sym(is)%translation(3,nt) )
            sgt%o%sym(is)%translation = sg%o%sym(is)%translation
          end do
          if (.not.trivial(sg)) then
            np = size(sg%o%gsmap,2)
            allocate( sgt%o%gsmap(3,np) )
            sgt%o%gsmap = sg%o%gsmap
            np = size(sg%o%gscount)
            allocate( sgt%o%gscount(0:np-1) )
            sgt%o%gscount = sg%o%gscount
            np = size(sg%o%gsdisp)
            allocate( sgt%o%gsdisp(0:np-1) )
            sgt%o%gsdisp = sg%o%gsdisp
            allocate( sgt%o%atom(ns) )
            na = size(sg%o%atom(1)%perm,1)
            do is = 1,ns
              allocate( sgt%o%atom(is)%perm(na,nt) )
              sgt%o%atom(is)%perm = sg%o%atom(is)%perm
            end do
            allocate( sgt%o%field(ns) )
            do is = 1,ns
              sgt%o%field(is)%phase_class = sg%o%field(is)%phase_class
              np = size(sg%o%field(is)%perm)
              allocate( sgt%o%field(is)%perm(np) )
              sgt%o%field(is)%perm = sg%o%field(is)%perm
            end do
            nc = size(sg%o%class)
            allocate( sgt%o%class(nc) )
            do ic = 1,nc
              np = size(sg%o%class(ic)%phase)
              allocate( sgt%o%class(ic)%phase(np) )
              sgt%o%class(ic)%phase = sg%o%class(ic)%phase
            end do
            allocate( sgt%o%local(ns) )
            nl = size(sg%o%local(1)%l)
            do is = 1,ns
              allocate( sgt%o%local(is)%l(0:nl-1) )
              do l = 0,(nl-1)
                allocate( sgt%o%local(is)%l(l)%mat(-l:l,-l:l) )
                sgt%o%local(is)%l(l)%mat = sg%o%local(is)%l(l)%mat
              end do
            end do
          end if
          sg%o%ref = sg%o%ref - sg%ref
          sg%o => sgt%o
          sg%o%ref = sg%o%ref + sg%ref
        end if
      end subroutine

      subroutine sort_by_magnitude_i(lat,nt,pure_trans,trans_map)
        type(lattice_obj) :: lat
        integer, dimension(:,:), intent(inout) :: trans_map
        real(double), dimension(:,:), intent(inout) :: pure_trans

        integer :: i, idx, inc, it, j, na, nt
        integer, dimension(:), allocatable :: idxs
        integer, dimension(:,:), allocatable :: map_old
        real(double) :: mag
        real(double), dimension(:), allocatable :: mags
        real(double), dimension(:,:), allocatable :: trans_old

        call my(lat)

        na = size(trans_map,1)

        allocate( mags(nt), idxs(nt) )
        do it = 1,nt
          mags(it) = norm(lat2r(lat,pure_trans(:,it)))
          idxs(it) = it
        end do

        inc = 1
        do
          inc = 3*inc + 1
          if (inc > nt) exit
        end do
        do
          inc = inc/3
          do i = inc+1,nt
            mag = mags(i)
            idx = idxs(i)
            j = i
            do
              if (mags(j-inc) <= mag) exit
              mags(j) = mags(j-inc)
              idxs(j) = idxs(j-inc)
              j = j - inc
              if (j <= inc) exit
            end do
            mags(j) = mag
            idxs(j) = idx
          end do
          if (inc <= 1) exit
        end do

        allocate( trans_old(3,nt), map_old(na,nt) )
        trans_old = pure_trans(:,1:nt)
        map_old   = trans_map(:,1:nt)
        do it = 1,nt
          pure_trans(:,it) = trans_old(:,idxs(it))
          trans_map(:,it) = map_old(:,idxs(it))
        end do
        deallocate( map_old, trans_old )

        deallocate( idxs, mags )

        call glean(thy(lat))

      end subroutine

      subroutine group_by_translation_i(lat,ns,temp_sym,temp_tau,temp_map)
        type(lattice_obj) :: lat
        integer, intent(in) :: ns
        integer, dimension(:,:), intent(inout) :: temp_map
        real(double), dimension(:,:), intent(inout) :: temp_tau
        real(double), dimension(:,:,:), intent(inout) :: temp_sym

        logical :: duplicate
        integer :: i, inc, is, it, j, na, ntd, ntu, tau_c
        integer, dimension(:), allocatable :: index, tau_count
        integer, dimension(:,:), allocatable :: map_old
        real(double) :: mag
        real(double), dimension(3) :: tau_u
        real(double), dimension(:), allocatable :: mags
        real(double), dimension(:,:), allocatable :: tau_old, tau_unique
        real(double), dimension(:,:,:), allocatable :: sym_old

        call my(lat)

        na = size(temp_map,1)

        allocate( tau_unique(3,ns), tau_count(ns) )

        ntu = 0
        tau_count = 0
        do is = 1,ns
          duplicate = .false.
          do it = 1,ntu
            if ( temp_tau(:,is) .in. nbhd(tau_unique(:,it),tol_equal) ) then
              duplicate = .true.
              ntd = it
              exit
            end if
          end do
          if (duplicate) then
            tau_count(ntd) = tau_count(ntd) + 1
          else
            ntu = ntu + 1
            tau_count(ntu) = tau_count(ntu) + 1
            tau_unique(:,ntu) = temp_tau(:,is)
          end if
        end do

        allocate( mags(ntu) )

        do it = 1,ntu
          mags(it) = norm(lat2r(lat,tau_unique(:,it)))
        end do

        inc = 1
        do
          inc = 3*inc + 1
          if (inc > ntu) exit
        end do
        do
          inc = inc/3
          do i = inc+1,ntu
            mag = mags(i)
            tau_u = tau_unique(:,i)
            tau_c = tau_count(i)
            j = i
            do
              if (mags(j-inc) <= mag) exit
              mags(j) = mags(j-inc)
              tau_unique(:,j) = tau_unique(:,j-inc)
              tau_count(j) = tau_count(j-inc)
              j = j - inc
              if (j <= inc) exit
            end do
            mags(j) = mag
            tau_unique(:,j) = tau_u
            tau_count(j) = tau_c
          end do
          if (inc <= 1) exit
        end do

        deallocate( mags )

        allocate( sym_old(3,3,ns), tau_old(3,ns), map_old(na,ns) )
        sym_old = temp_sym(:,:,1:ns)
        tau_old = temp_tau(:,1:ns)
        map_old = temp_map(:,1:ns)

        allocate( index(ntu) )
        do it = 1,ntu
          index(it) = 0
          do i = 1,it-1
            index(it) = index(it) + tau_count(i)
          end do
        end do

        do is = 1,ns
          do it = 1,ntu
            if ( tau_old(:,is) .in. nbhd(tau_unique(:,it),tol_equal) ) then
              index(it) = index(it) + 1
              temp_sym(:,:,index(it)) = sym_old(:,:,is)
              temp_tau(:,index(it))   = tau_old(:,is)
              temp_map(:,index(it))   = map_old(:,is)
              exit
            end if
          end do
        end do

        deallocate( index )
        deallocate( map_old, tau_old, sym_old )
        deallocate( tau_count, tau_unique )

        call glean(thy(lat))

      end subroutine

      subroutine check_mesh_dimensions_i(sgr,lay)
        type(space_group_rep) :: sgr
        type(layout_obj) :: lay
!       effects: Checks that mesh dimensions are equal along lattice vectors related by symmetry.

        logical :: related_12, related_13, related_23
        integer :: is, ns
        integer, dimension(3) :: dims
        real(double), dimension(3) :: a1, a2, a3

        call my(lay)

        ns = size(sgr%sym)
        dims = x_dims(lay)

        a1 = real((/1,0,0/),double)
        a2 = real((/0,1,0/),double)
        a3 = real((/0,0,1/),double)

        related_12 = .false.
        do is = 1,ns
          related_12 = related_12 .or. (matmul(sgr%sym(is)%point_op,a1) .in. nbhd(a2,tol_equal))
          if (related_12) exit
        end do
        related_13 = .false.
        do is = 1,ns
          related_13 = related_13 .or. (matmul(sgr%sym(is)%point_op,a1) .in. nbhd(a3,tol_equal))
          if (related_13) exit
        end do
        related_23 = .false.
        do is = 1,ns
          related_23 = related_23 .or. (matmul(sgr%sym(is)%point_op,a2) .in. nbhd(a3,tol_equal))
          if (related_23) exit
        end do

        if (error(related_12 .and. (dims(1) /= dims(2)),"ERROR: mesh dimensions 1 and 2 are not equal")) continue
        if (error(related_13 .and. (dims(1) /= dims(3)),"ERROR: mesh dimensions 1 and 3 are not equal")) continue
        if (error(related_23 .and. (dims(2) /= dims(3)),"ERROR: mesh dimensions 2 and 3 are not equal")) goto 100

100     call glean(thy(lay))

        if (error("Exit symmetry_mod::check_mesh_dimensions_i")) continue

      end subroutine

      subroutine form_symmetrizing_structures_i(sgr,tau,pg_map,trans_map,lat,lay)
        type(space_group_rep) :: sgr
        real(double), dimension(:,:), intent(in) :: tau
        integer, dimension(:,:), intent(in) :: pg_map
        integer, dimension(:,:), intent(in) :: trans_map
        type(lattice_obj) :: lat
        type(layout_obj) :: lay
!       effects: Forms structures used to symmetrize objects.

        integer :: ns

        call my(lat)
        call my(lay)

        ns = size(sgr%sym)

        allocate( sgr%atom(ns) )
        call form_atom_permutations_i(sgr,pg_map,trans_map)

        call distribute_field_points_i(sgr,lay)
        allocate( sgr%field(ns) )
        call form_field_permutations_i(sgr,lay) ; if (error()) goto 100
        call form_class_phases_i(sgr,tau,lay)   ; if (error()) goto 100

        allocate( sgr%local(ns) )
        call form_local_transforms_i(sgr,lat)

100     call glean(thy(lat))
        call glean(thy(lay))

        if (error("Exit symmetry_mod::form_symmetrizing_structures_i")) continue

      end subroutine

      subroutine form_atom_permutations_i(sgr,pg_map,trans_map)
        type(space_group_rep) :: sgr
        integer, dimension(:,:), intent(in) :: pg_map
        integer, dimension(:,:), intent(in) :: trans_map

        integer :: ia, is, it, na, ns, nt

        na = size(pg_map,1)
        ns = size(sgr%sym)
        nt = size(sgr%sym(1)%translation,2)

        do is = 1,ns
          allocate( sgr%atom(is)%perm(na,nt) )
          do it = 1,nt
            do ia = 1,na
              sgr%atom(is)%perm(ia,it) = trans_map(pg_map(ia,is),it)
            end do
          end do
        end do

      end subroutine

      subroutine distribute_field_points_i(sgr,lay)
        type(space_group_rep) :: sgr
        type(layout_obj) :: lay

        integer :: i, j, i1, i2, i3, f_ps, l_ps, n_ps, n_ss
        integer, dimension(3) :: a, m
        real(double) :: cutoff, g2

        call my(lay)

        m = x_dims(lay)

        cutoff = x_cutoff(lay)

        ! Count the number of G-vectors in the sub-set (ss) having energies within the cutoff
        n_ss = 0
        do i3 = 1,m(3)
        do i2 = 1,m(2)
        do i1 = 1,m(1)
          a = (/i1,i2,i3/)
          g2 = a2f2(a,lay,S_TYPE)
          if (g2 > cutoff) cycle
          n_ss = n_ss + 1
        end do
        end do
        end do

        ! Form a 3-d <-> 1d mapping for use in form_field_permutations_i and symmetrize_grid
        allocate( sgr%gsmap(3,n_ss) )
        n_ss = 0
        do i3 = 1,m(3)
        do i2 = 1,m(2)
        do i1 = 1,m(1)
          a = (/i1,i2,i3/)
          g2 = a2f2(a,lay,S_TYPE)
          if (g2 > cutoff) cycle
          n_ss = n_ss + 1
          sgr%gsmap(:,n_ss) = a
        end do
        end do
        end do

        ! Divide the G-vectors into SGROUP process-sets (ps)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,n_ss,f_ps,l_ps,n_ps)

        ! Form indexing arrays for use in symmetrize_grid
        allocate( sgr%gscount(0:mpi_nprocs(SGROUP)-1), sgr%gsdisp(0:mpi_nprocs(SGROUP)-1) )
        call allgather(SGROUP,n_ps,sgr%gscount)
        sgr%gsdisp = 0
        do i = 0,(mpi_nprocs(SGROUP)-1)
          do j = 0,(i-1)
            sgr%gsdisp(i) = sgr%gsdisp(i) + sgr%gscount(j)
          end do
        end do

        call glean(thy(lay))

      end subroutine

      subroutine form_field_permutations_i(sgr,lay)
        type(space_group_rep) :: sgr
        type(layout_obj) :: lay

        integer :: i1, i2, i3, i_ps, f_ps, l_ps, n_ps, is, ns, i_ss, n_ss 
        integer, dimension(3) :: a, m, o, ta, tgl
        integer, dimension(3,3) :: mat
        integer, dimension(:,:), allocatable :: gl
        integer, dimension(:,:,:), allocatable :: map
        real(double) :: cutoff, g2

        call my(lay)

        o = (/1,1,1/)
        m = x_dims(lay)

        ns = size(sgr%field)
        n_ss = size(sgr%gsmap,2)

        cutoff = x_cutoff(lay)

        ! Form a 3-d to 1-d mapping for the sub-set (ss) of G-vectors having energies within the cutoff
        allocate( map(m(1),m(2),m(3)) )
        i_ss = 0
        do i3 = 1,m(3)
        do i2 = 1,m(2)
        do i1 = 1,m(1)
          a = (/i1,i2,i3/)
          g2 = a2f2(a,lay,S_TYPE)
          if (g2 > cutoff) then
            map(i1,i2,i3) = 0
          else
            i_ss = i_ss + 1
            map(i1,i2,i3) = i_ss
          end if
        end do
        end do
        end do

        ! Divide the ss of G-vectors into SGROUP process-sets (ps)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,n_ss,f_ps,l_ps,n_ps)

        ! Form the G-vectors (lattice representation) in the ps
        allocate( gl(3,n_ps) )
        do i_ps = 1,n_ps
          i_ss = i_ps + f_ps - 1
          a = sgr%gsmap(:,i_ss)
          gl(:,i_ps) = a2fl(a,lay,S_TYPE)
        end do

        ! Form the G-vector permutations
        do is = 1,ns
          allocate( sgr%field(is)%perm(n_ps) )
          ! inverse transform: (TAT^-1)*tgl = gl where (TA^-1T^-1) = (TA^~T^-1) = (SAS^-1)^~
          mat = nint(transpose(sgr%sym(is)%point_op))
          do i_ps = 1,n_ps
            tgl(1) = mat(1,1)*gl(1,i_ps) + mat(1,2)*gl(2,i_ps) + mat(1,3)*gl(3,i_ps)
            tgl(2) = mat(2,1)*gl(1,i_ps) + mat(2,2)*gl(2,i_ps) + mat(2,3)*gl(3,i_ps)
            tgl(3) = mat(3,1)*gl(1,i_ps) + mat(3,2)*gl(2,i_ps) + mat(3,3)*gl(3,i_ps)
            ta = fl2a(tgl,lay,S_TYPE)
            if (error(any(ta < o) .or. any(ta > m),"ERROR: transformed 3-d array indices are out-of-bounds")) goto 100
            i_ss = map(ta(1),ta(2),ta(3))
            if (error((i_ss < 1) .or. (i_ss > n_ss),"ERROR: transformed 1-d array indices are out-of-bounds")) goto 100
            sgr%field(is)%perm(i_ps) = i_ss
          end do
        end do

100     if (allocated( gl )) deallocate( gl )
        if (allocated( map )) deallocate( map )

        call glean(thy(lay))

        if (error("Exit symmetry_mod::form_field_permutations_i")) continue

      end subroutine

      subroutine form_class_phases_i(sgr,tau,lay)
        type(space_group_rep) :: sgr
        real(double), dimension(:,:), intent(in) :: tau
        type(layout_obj) :: lay

        integer :: ic, nc, i_ps, f_ps, l_ps, n_ps, is, ns, i_ss, n_ss, it, nt
        integer, dimension(3) :: a
        integer, dimension(:,:), allocatable :: gl
        real(double) :: dp, scale
        real(double), dimension(3) :: rgl
        real(double), dimension(:,:,:), allocatable :: tl
        complex(double) :: c0
        complex(double), dimension(:), allocatable :: pt

        call my(lay)

        ns = size(sgr%field)
        nt = size(sgr%sym(1)%translation,2)

        call form_class_structure_i(sgr,tau)
        nc = 0
        do is = 1,ns
          nc = max(nc,sgr%field(is)%phase_class)
        end do

        ! Divide the sub-set (ss) of G-vectors with energy within the cutoff into SGROUP process-sets (ps)
        n_ss = size(sgr%gsmap,2)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,n_ss,f_ps,l_ps,n_ps)

        ! Form the G-vectors (lattice representation) in the ps
        allocate( gl(3,n_ps) )
        do i_ps = 1,n_ps
          i_ss = i_ps + f_ps - 1
          a = sgr%gsmap(:,i_ss)
          gl(:,i_ps) = a2fl(a,lay,S_TYPE)
        end do

        ! Form the translations (lattice representation)
        allocate( tl(3,nt,nc) )
        ic = 1
        tl(:,:,ic) = sgr%sym(1)%translation
        do is = 2,ns
          if (sgr%field(is)%phase_class > sgr%field(is-1)%phase_class) then
            ic = ic + 1
            tl(:,:,ic) = sgr%sym(is)%translation
          end if
        end do

        ! Form the class phases
        scale = 1.0_double/real(nt,double)
        c0 = 2.0_double*pi*cmplx(0,-1,double)
        allocate( sgr%class(nc) )
        allocate( pt(n_ps) )
        do ic = 1,nc
          allocate( sgr%class(ic)%phase(n_ps) )
          pt = cmplx(0,0,double)
          do it = 1,nt
            do i_ps = 1,n_ps
              rgl = real(gl(:,i_ps),double)
              dp = tl(1,it,ic)*rgl(1) + tl(2,it,ic)*rgl(2) + tl(3,it,ic)*rgl(3)
              pt(i_ps) = pt(i_ps) + exp(c0*dp)
            end do
          end do
          sgr%class(ic)%phase = scale*pt
        end do

        if (allocated( gl )) deallocate( gl )
        if (allocated( tl )) deallocate( tl )
        if (allocated( pt )) deallocate( pt )

        call glean(thy(lay))

      end subroutine

      subroutine form_class_structure_i(sgr,tau)
        type(space_group_rep) :: sgr
        real(double), dimension(:,:), intent(in) :: tau

        logical :: duplicate
        integer :: is, it, ns, ntd, ntu
        real(double), allocatable, dimension(:,:) :: tau_unique

        ns = size(sgr%field)
        allocate( tau_unique(3,ns) )

        ntu = 0
        do is = 1,ns
          duplicate = .false.
          do it = 1,ntu
            if ( tau(:,is) .in. nbhd(tau_unique(:,it),tol_equal) ) then
              duplicate = .true.
              ntd = it
              exit
            end if
          end do
          if (duplicate) then
            sgr%field(is)%phase_class = ntd
          else
            ntu = ntu + 1
            sgr%field(is)%phase_class = ntu
            tau_unique(:,ntu) = tau(:,is)
          end if
        end do

        if (allocated( tau_unique )) deallocate( tau_unique )

      end subroutine

      subroutine form_local_transforms_i(sgr,lat)
        type(space_group_rep) :: sgr
        type(lattice_obj) :: lat

        integer :: is, l, mc, mr, ns, sgn
        real(double) :: cos_b
        real(double), dimension(3,3) :: mat
        complex(double) :: exp_mia, exp_mig

        call my(lat)

        ns = size(sgr%sym)

        do is = 1,ns
          allocate( sgr%local(is)%l(0:l_max) )
          mat = lat2r(lat,sgr%sym(is)%point_op)
          sgn = nint(determinant(mat))
          mat = real(sgn,double)*mat
          call euler_angles_i(mat,cos_b,exp_mia,exp_mig)
          do l = 0,l_max
            allocate( sgr%local(is)%l(l)%mat(-l:l,-l:l) )
            do mc = -l,l
              do mr = -l,l
                sgr%local(is)%l(l)%mat(mr,mc) = real(sgn**l,double)*(exp_mia**mr)*d_beta_i(cos_b,l,mr,mc)*(exp_mig**mc)
              end do
            end do
          end do
        end do

        call glean(thy(lat))

      end subroutine

      subroutine get_point_group_type_i(op_type,st)
        character(*) :: st
        integer, dimension(10) :: op_type

        if (     all(op_type == (/1,0,0,0,0,0,0,0,0,0/)) ) then
           st = 'C_1'
        elseif ( all(op_type == (/1,0,0,0,0,1,0,0,0,0/)) ) then
           st = 'S_2'
        elseif ( all(op_type == (/1,1,0,0,0,0,0,0,0,0/)) ) then
           st = 'C_2'
        elseif ( all(op_type == (/1,0,0,0,0,0,1,0,0,0/)) ) then
           st = 'C_1h'
        elseif ( all(op_type == (/1,1,0,0,0,1,1,0,0,0/)) ) then
           st = 'C_2h'
        elseif ( all(op_type == (/1,3,0,0,0,0,0,0,0,0/)) ) then
           st = 'D_2'
        elseif ( all(op_type == (/1,1,0,0,0,0,2,0,0,0/)) ) then
           st = 'C_2v'
        elseif ( all(op_type == (/1,3,0,0,0,1,3,0,0,0/)) ) then
           st = 'D_2h'
        elseif ( all(op_type == (/1,1,0,2,0,0,0,0,0,0/)) ) then
           st = 'C_4'
        elseif ( all(op_type == (/1,1,0,0,0,0,0,0,2,0/)) ) then
           st = 'S_4'
        elseif ( all(op_type == (/1,1,0,2,0,1,1,0,2,0/)) ) then
           st = 'C_4h'
        elseif ( all(op_type == (/1,5,0,2,0,0,0,0,0,0/)) ) then
           st = 'D_4'
        elseif ( all(op_type == (/1,1,0,2,0,0,4,0,0,0/)) ) then
           st = 'C_4v'
        elseif ( all(op_type == (/1,3,0,0,0,0,2,0,2,0/)) ) then
           st = 'D_2d'
        elseif ( all(op_type == (/1,5,0,2,0,1,5,0,2,0/)) ) then
           st = 'D_4h'
        elseif ( all(op_type == (/1,0,2,0,0,0,0,0,0,0/)) ) then
           st = 'C_3'
        elseif ( all(op_type == (/1,0,2,0,0,1,0,2,0,0/)) ) then
           st = 'S_6'
        elseif ( all(op_type == (/1,3,2,0,0,0,0,0,0,0/)) ) then
           st = 'D_3'
        elseif ( all(op_type == (/1,0,2,0,0,0,3,0,0,0/)) ) then
           st = 'C_3v'
        elseif ( all(op_type == (/1,3,2,0,0,1,3,2,0,0/)) ) then
           st = 'D_3d'
        elseif ( all(op_type == (/1,1,2,0,2,0,0,0,0,0/)) ) then
           st = 'C_6'
        elseif ( all(op_type == (/1,0,2,0,0,0,1,0,0,2/)) ) then
           st = 'C_3h'
        elseif ( all(op_type == (/1,1,2,0,2,1,1,2,0,2/)) ) then
           st = 'C_6h'
        elseif ( all(op_type == (/1,7,2,0,2,0,0,0,0,0/)) ) then
           st = 'D_6'
        elseif ( all(op_type == (/1,1,2,0,2,0,6,0,0,0/)) ) then
           st = 'C_6v'
        elseif ( all(op_type == (/1,3,2,0,0,0,4,0,0,2/)) ) then
           st = 'D_3h'
        elseif ( all(op_type == (/1,7,2,0,2,1,7,2,0,2/)) ) then
           st = 'D_6h'
        elseif ( all(op_type == (/1,3,8,0,0,0,0,0,0,0/)) ) then
           st = 'T'
        elseif ( all(op_type == (/1,3,8,0,0,1,3,8,0,0/)) ) then
           st = 'T_h'
        elseif ( all(op_type == (/1,9,8,6,0,0,0,0,0,0/)) ) then
           st = 'O'
        elseif ( all(op_type == (/1,3,8,0,0,0,6,0,6,0/)) ) then
           st = 'T_d'
        elseif ( all(op_type == (/1,9,8,6,0,1,9,8,6,0/)) ) then
           st = 'O_h'
        else
           if (error(.true.,"ERROR: unknown point-group type")) goto 100
        end if

100     if (error("Exit symmetry_mod::get_point_group_type_i")) continue

      end subroutine

      subroutine get_lattice_type_i(op_type,st)
        character(*) :: st
        integer, dimension(10) :: op_type

        if (     all(op_type == (/1,0,0,0,0,1,0,0,0,0/)) ) then
           st = 'The lattice is triclinic with point group S_2'
        elseif ( all(op_type == (/1,1,0,0,0,1,1,0,0,0/)) ) then
           st = 'The lattice is monoclinic with point group C_2h'
        elseif ( all(op_type == (/1,3,0,0,0,1,3,0,0,0/)) ) then
           st = 'The lattice is orthorhombic with point group D_2h'
        elseif ( all(op_type == (/1,5,0,2,0,1,5,0,2,0/)) ) then
           st = 'The lattice is tetragonal with point group D_4h'
        elseif ( all(op_type == (/1,3,2,0,0,1,3,2,0,0/)) ) then
           st = 'The lattice is trigonal (rhombohedral) with point group D_3d'
        elseif ( all(op_type == (/1,7,2,0,2,1,7,2,0,2/)) ) then
           st = 'The lattice is hexagonal with point group D_6h'
        elseif ( all(op_type == (/1,9,8,6,0,1,9,8,6,0/)) ) then
           st = 'The lattice is cubic with point group O_h'
        else
           if (error(.true.,"ERROR: unknown lattice type")) goto 100
        end if

100     if (error("Exit symmetry_mod::get_lattice_type_i")) continue

      end subroutine

      function classify_op_i(mat) result(cl)
        integer, dimension(10) :: cl
        real(double), dimension(3,3), intent(in) :: mat  ! represented in lattice coordinates
!       classification scheme is (E,C2,C3,C4,C6,I,M,S6,S4,S3)

        cl = 0
        select case(nint(determinant(mat)))
        case(1)
           select case(nint(trace(mat)))
           case(3)
              cl(1) = 1
           case(-1)
              cl(2) = 1
           case(0)
              cl(3) = 1
           case(+1)
              cl(4) = 1
           case(+2)
              cl(5) = 1
           case default
              if (error(.true.,"ERROR: incorrect trace")) goto 100
           end select
        case(-1)
           select case(nint(trace(mat)))
           case(-3)
              cl(6) = 1
           case(+1)
              cl(7) = 1
           case(0)
              cl(8) = 1
           case(-1)
              cl(9) = 1
           case(-2)
              cl(10) = 1
           case default
              if (error(.true.,"ERROR: incorrect trace")) goto 100
           end select
        case default
           if (error(.true.,"ERROR: incorrect determinant")) goto 100
        end select

100     if (error("Exit symmetry_mod::classify_op_i")) return

      end function

      subroutine get_mag_list_i(v,mag_num)
        integer, dimension(:), intent(out) :: mag_num
        real(double), dimension(:,:), intent(in) :: v

        logical :: mapped
        integer :: iv, ivt, nv, nvt
        real(double), dimension(:), allocatable :: mag, mag_list

        nv = size(v,2)

        allocate( mag(nv), mag_list(nv) )
        nvt = 0
        do iv = 1,nv
          mag(iv) = norm(v(:,iv))
          mapped = .false.
          do ivt = 1,nvt
            if (mag(iv) .in. nbhd(mag_list(ivt),tol_equal)) then
              mapped = .true.
              mag_num(iv) = ivt
              exit
            end if
          end do
          if (.not.mapped) then
            nvt = nvt + 1
            mag_num(iv) = nvt
            mag_list(nvt) = mag(iv)
          end if
        end do
        deallocate( mag_list, mag )

      end subroutine

      subroutine euler_angles_i(mat,cos_b,exp_mia,exp_mig)
        real(double), dimension(3,3), intent(in) :: mat
        real(double), intent(out) :: cos_b
        complex(double), intent(out) :: exp_mia, exp_mig
!       requires: mat be in the {x,y,z} representation and must represent a proper rotation (determinant(mat) = +1).

        real(double), parameter :: tol_cos_b = 1.0e-6_double
        real(double) :: mag_sin_b

        cos_b = mat(3,3)
        if (abs(cos_b) .in. nbhd(1.0_double,tol_cos_b)) then
          exp_mia = cmplx(mat(1,1),-mat(2,1),double)/cos_b
          exp_mig = cmplx(1,0,double)
        else
          mag_sin_b = sqrt(1.0_double - cos_b**2)
          exp_mia = cmplx(mat(1,3),-mat(2,3),double)/mag_sin_b
          exp_mig = cmplx(-mat(3,1),-mat(3,2),double)/mag_sin_b
        end if

      end subroutine

      function d_beta_i(cos_b,l,mp,m) result(db)
        real(double), intent(in) :: cos_b
        integer, intent(in) :: l, mp, m
        real(double) :: db
!       requires: -1 <= cos_b <= +1, l >= 0, -l <= mp <= +l, -l <= m <= +l.

        real (double), parameter :: tol_cos_b = 1.0e-6_double
        integer :: a, b, c, i
        real (double) :: m_tan_hb_2, mag_cos_hb, mag_tan_hb, s, sum, term

        if (cos_b .in. nbhd(1.0_double,tol_cos_b)) then

          db = 0.0_double
          if (mp == m) db = 1.0_double

        elseif (cos_b .in. nbhd(-1.0_double,tol_cos_b)) then

          db = 0.0_double
          if (mp == -m) db = real((-1)**(l+m),double)

        else

          mag_cos_hb = sqrt( 0.5_double*(1.0_double + cos_b) )
          mag_tan_hb = sqrt(1.0_double - cos_b**2)/(1.0_double + cos_b)
          m_tan_hb_2 = -(mag_tan_hb**2)
          if (mp >= m) then
            a =  (mp-l)
            b = -(m+l)
            c =  (mp-m+1)
            s = 1.0_double
          else
            a =  (m-l)
            b = -(mp+l)
            c =  (m-mp+1)
            s = real((-1)**(c-1),double)
          end if
          sum = 1.0_double
          term = 1.0_double
          i = 0
          do while( ((a+i) < 0) .and. ((b+i) < 0) )
            term = term*real((a+i)*(b+i),double)*m_tan_hb_2/real((1+i)*(c+i),double)
            sum = sum + term
            i = i + 1
          end do
          db = real(factorial_i(c-a-1),double)/real(factorial_i(-b),double)
          db = db*(real(factorial_i(c-b-1),double)/real(factorial_i(-a),double))
          db = sqrt(db)
          db = db*real((-1)**(c-1),double)/real(factorial_i(c-1),double)
          db = db*s*sum
          db = db*(mag_cos_hb**(c-a-b-1))
          db = db*(mag_tan_hb**(c-1))
          
        end if

      end function

      function factorial_i(n) result(f)
        integer, intent(in) :: n
        integer :: f
!       requires: n >= 0

        integer :: i
        f = 1
        do i = 2,n
          f = f*i
        end do
      end function

      subroutine save_sg_i(sgr)
        type(space_group_rep) :: sgr

        integer :: ios, is, it, ns, nt
        type(file_obj) :: f

        call my(file(trim(new_space_group_path)),f)

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100

        ns = size(sgr%sym)
        nt = size(sgr%sym(1)%translation,2)
        if (i_access(f)) then
          write(x_unit(f),'(i0,2x,i0)') ns, nt
          do is = 1,ns
            write(x_unit(f),'(a)') " "
            write(x_unit(f),'(3i5)') nint(sgr%sym(is)%point_op(1,:))
            write(x_unit(f),'(3i5)') nint(sgr%sym(is)%point_op(2,:))
            write(x_unit(f),'(3i5)') nint(sgr%sym(is)%point_op(3,:))
            do it = 1,nt
              write(x_unit(f),'(3f16.10)') sgr%sym(is)%translation(:,it)
            end do
          end do
        end if

        if (i_access(f)) close(x_unit(f))

100     call glean(thy(f))

        if (error("Exit symmetry_mod::save_sg_i")) continue

      end subroutine

      end module

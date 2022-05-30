! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module lattice_mod
!doc$ module lattice_mod

!     One datatype is available here: type(lattice_obj).

!     Lattice_mod contains information about the lattice vectors {a1,a2,a3} and the
!     reciprocal lattice vectors {b1,b2,b3} defined such that ai*bj = 2*pi*delta_ij.

      use kind_mod
      use mpi_mod
      use arg_mod
      use error_mod
      use io_mod
      use diary_mod
      use math_mod
      use ghost_mod
      use tagio_mod

!cod$
      implicit none
      private

      integer, parameter :: PS = 1   ! Position Space
      integer, parameter :: RS = 2   ! Reciprocal Space

      type :: lattice_rep
        integer :: ref
        type(ghost) :: g
        real(double) :: lconstant
        real(double), dimension(3,3) :: vectors
        real(double), dimension(3,3) :: ivectors
      end type

      type, public :: lattice_obj
        private
        integer :: ref
        type(lattice_rep), pointer :: o
      end type

!doc$
      public :: lattice
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_cell_volume
      public :: x_lattice_constant
      public :: x_lattice_vector
      public :: r2lat
      public :: r2lat_tensor
      public :: lat2r
      public :: lat2r_tensor
      public :: f2lat
      public :: lat2f
      public :: save
      public :: cubic_cell
      public :: wigner_seitz_vectors
      public :: maximum_sphere_radius
      public :: diary
      public :: write_restart

!cod$
      interface lattice
        module procedure constructor_lat_1, constructor_lat_2, constructor_lat_3
      end interface
      interface my
        module procedure my_lat, my_new_lat
      end interface
      interface thy
        module procedure thy_lat
      end interface
      interface glean
        module procedure glean_lat
      end interface
      interface bequeath
        module procedure bequeath_lat
      end interface
      interface assignment(=)
        module procedure assign_lat
      end interface
      interface x_ref
        module procedure lat_ref
      end interface
      interface x_ghost
        module procedure lat_ghost
      end interface
      interface x_cell_volume
        module procedure lat_cell_volume
      end interface
      interface x_lattice_constant
        module procedure lat_lattice_constant
      end interface
      interface x_lattice_vector
        module procedure lat_lattice_vector
      end interface
      interface r2lat
        module procedure r2lat_vector, r2lat_matrix
      end interface
      interface lat2r
        module procedure lat2r_vector, lat2r_matrix
      end interface
      interface f2lat
        module procedure f2lat_vector, f2lat_matrix
      end interface
      interface lat2f
        module procedure lat2f_vector, lat2f_matrix
      end interface
      interface save
        module procedure save_lat
      end interface
      interface diary
        module procedure diary_lat
      end interface
      interface write_restart
        module procedure write_restart_lat
      end interface

      contains

! public routines

      function constructor_lat_1(v1,v2,v3,latc) result(lat)
!doc$ function lattice(v1,v2,v3,latc) result(lat)
        type(lattice_obj) :: lat
        real(double), intent(in) :: latc
        real(double), dimension(3), intent(in) :: v1, v2, v3
!       requires: v1, v2, v3 be linearly independent.
!       effects: Creates a new lat based on v1, v2, v3 and latc.

!cod$
        real(double), dimension(3,3) :: iv

        lat%ref = 0
        allocate( lat%o )
        lat%o%ref = 0
        lat%o%g = x_ghost()
        lat%o%lconstant = latc
        lat%o%vectors(:,1) = latc*v1
        lat%o%vectors(:,2) = latc*v2
        lat%o%vectors(:,3) = latc*v3
        iv = inverse(lat%o%vectors)
        lat%o%ivectors = iv
      end function

      function constructor_lat_2(prefix) result(lat)
!doc$ function lattice(prefix) result(lat)
        type(lattice_obj) :: lat
        character(*), intent(in) :: prefix
!       requires: v1, v2, v3 be linearly independent.
!       effects: Creates a new lat based on v1, v2, v3, and latc read from file ./prefix/lattice.
!       errors: Format problems and I/O problems.

!cod$
        logical :: exist_file
        integer :: ios
        real(double) :: latc
        real(double), dimension(3) :: v1, v2, v3
        real(double), dimension(3,3) :: iv
        type(file_obj) :: f

        call my(file(prefix//"lattice"),f)
        if (i_access(f)) inquire(file=x_name(f),exist=exist_file)
        if (i_comm(f)) call broadcast(FILE_SCOPE,exist_file)
        if (error(.not.exist_file,"ERROR: file does not exist")) goto 200

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 200

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

        lat%ref = 0
        allocate( lat%o )
        lat%o%ref = 0
        lat%o%g = x_ghost()
        lat%o%lconstant = latc
        lat%o%vectors(:,1) = latc*v1
        lat%o%vectors(:,2) = latc*v2
        lat%o%vectors(:,3) = latc*v3
        iv = inverse(lat%o%vectors)
        lat%o%ivectors = iv

100     if (i_access(f)) close(x_unit(f))
200     call glean(thy(f))

        if (error("Exit lattice_mod::constructor_lat_2")) continue

      end function

      function constructor_lat_3(restf) result(lat)
!doc$ function lattice(restf) result(lat)
        type(tagio_obj) :: restf
        type(lattice_obj) :: lat
!       requires: restf be positioned inside the CRYSTAL block.
!       effects: Creates a new lat from information in restf.
!       errors: LATTICE tag not found.

!cod$
        character(1) :: tios
        integer(long) :: dsize, iosl, ndata

        call my(restf)

        lat%ref = 0
        allocate( lat%o )
        lat%o%ref = 0
        lat%o%g = x_ghost()

        ! open the LATTICE block
        if (i_access(restf)) tios = findfirsttag(restf,"LATTICE")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: LATTICE block was not found")) goto 200
        if (i_access(restf)) call openblock(restf)

        ! read the lattice constant       
        if (i_access(restf)) tios = findfirsttag(restf,"LATTICE_CONSTANT")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: LATTICE_CONSTANT tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_double ; ndata = 1
          call readf(lat%o%lconstant,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,lat%o%lconstant)

        ! read the real-space vectors
        if (i_access(restf)) tios = findfirsttag(restf,"REAL-SPACE_VECTORS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: REAL-SPACE_VECTORS tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_double ; ndata = 9
          call readf(lat%o%vectors,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,lat%o%vectors)

        ! read the reciprocal-space vectors
        if (i_access(restf)) tios = findfirsttag(restf,"RECIPROCAL-SPACE_VECTORS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: RECIPROCAL-SPACE_VECTORS tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_double ; ndata = 9
          call readf(lat%o%ivectors,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,lat%o%ivectors)

        ! close the LATTICE block
100     if (i_access(restf)) call closeblock(restf)

200     call glean(thy(restf))

        if (error("Exit lattice_mod::constructor_lat_3")) continue

      end function

      subroutine my_lat(lat)
!doc$ subroutine my(lat)
        type(lattice_obj) :: lat

!cod$
        lat%ref = lat%ref + 1
        lat%o%ref = lat%o%ref + 1
      end subroutine

      subroutine my_new_lat(lati,lat)
!doc$ subroutine my(lati,lat)
        type(lattice_obj) :: lati, lat

!cod$
        lat%ref = 1
        lat%o => lati%o
        lat%o%ref = lat%o%ref + 1
      end subroutine

      function thy_lat(lat) result(lato)
!doc$ function thy(lat) result(lato)
        type(lattice_obj) :: lat, lato

!cod$
        lat%ref = lat%ref - 1
        lat%o%ref = lat%o%ref - 1
        lato%ref = lat%ref
        lato%o => lat%o
      end function

      subroutine glean_lat(lat)
!doc$ subroutine glean(lat)
        type(lattice_obj) :: lat

!cod$
        if (lat%o%ref < 1) then
          deallocate( lat%o )
        end if
      end subroutine

      subroutine bequeath_lat(lat)
!doc$ subroutine bequeath(lat)
        type(lattice_obj) :: lat

!cod$
        continue
      end subroutine

      subroutine assign_lat(lat,lat2)
!doc$ subroutine assignment(=)(lat,lat2)
        type(lattice_obj), intent(inout) :: lat
        type(lattice_obj), intent(in) :: lat2

!cod$
        type(lattice_obj) :: latt
        call my(lat2)
        latt%o => lat%o
        lat%o%ref = lat%o%ref - lat%ref
        lat%o => lat2%o
        lat%o%ref = lat%o%ref + lat%ref
        call glean(latt)
        call glean(thy(lat2))
      end subroutine

      function lat_ref(lat) result(r)
!doc$ function x_ref(lat) result(r)
        type(lattice_obj) :: lat
        integer, dimension(2) :: r
!       effects: Returns lat%ref and lat%o%ref.

!cod$
        r(1) = lat%ref
        r(2) = lat%o%ref
        call glean(lat)
      end function

      function lat_ghost(lat) result(g)
!doc$ function x_ghost(lat) result(g)
        type(lattice_obj) :: lat
        type(ghost) :: g
!       effects: Returns the ghost of lat.

!cod$
        call my(lat)
        g = lat%o%g
        call glean(thy(lat))
      end function

      function lat_cell_volume(lat) result(v)
!doc$ function x_cell_volume(lat) result(v)
        type(lattice_obj) :: lat
        real (double) :: v
!       effects: Returns the cell volume of lat.

!cod$
        call my(lat)
        v = abs(determinant(lat%o%vectors))
        call glean(thy(lat))
      end function

      function lat_lattice_constant(lat) result(c)
!doc$ function x_lattice_constant(lat) result(c)
        type(lattice_obj) :: lat
        real (double) :: c
!       effects: Returns the lattice constant of lat.

!cod$
        call my(lat)
        c = lat%o%lconstant
        call glean(thy(lat))
      end function

      function lat_lattice_vector(lat,i) result(v)
!doc$ function x_lattice_vector(lat,i) result(v)
        type(lattice_obj) :: lat
        integer, intent(in) :: i
        real (double), dimension(3) :: v
!       effects: Returns the i'th lattice vector.

!cod$
        call my(lat)
        v = lat%o%vectors(:,i)
        call glean(thy(lat))
      end function

      function r2lat_vector(lat,v) result(tv)
!doc$ function r2lat(lat,v) result(tv)
        type(lattice_obj) :: lat
        real(double), dimension(:), intent(in) :: v
        real(double), dimension(3) :: tv
!       requires: v be dimension(3).
!       effects: Transforms v from the {x,y,z} representation to the {a1,a2,a3} representation where
!                a1, a2, and a3 are direct lattice vectors.

!cod$
        call my(lat)
        tv(1) = lat%o%ivectors(1,1)*v(1) + lat%o%ivectors(1,2)*v(2) + lat%o%ivectors(1,3)*v(3)
        tv(2) = lat%o%ivectors(2,1)*v(1) + lat%o%ivectors(2,2)*v(2) + lat%o%ivectors(2,3)*v(3)
        tv(3) = lat%o%ivectors(3,1)*v(1) + lat%o%ivectors(3,2)*v(2) + lat%o%ivectors(3,3)*v(3)
        call glean(thy(lat))
      end function

      function r2lat_matrix(lat,m) result(tm)
!doc$ function r2lat(lat,m) result(tm)
        type(lattice_obj) :: lat
        real(double), dimension(:,:), intent(in) :: m
        real(double), dimension(3,3) :: tm
!       requires: m be dimension(3,3).
!       effects: Transforms m from the {x,y,z} representation to the {a1,a2,a3}
!                representation where a1, a2, and a3 are direct lattice vectors.

!cod$
        call my(lat)
        tm = matmul(lat%o%ivectors,matmul(m,lat%o%vectors))
        call glean(thy(lat))
      end function

      function r2lat_tensor(lat,t) result(tt)
!doc$ function r2lat_tensor(lat,t) result(tt)
        type(lattice_obj) :: lat
        real(double), dimension(:,:), intent(in) :: t
        real(double), dimension(3,3) :: tt
!       requires: t be dimension(3,3).
!       effects: Transforms t from the {x,y,z} representation to the {a1,a2,a3}
!                representation where a1, a2, and a3 are direct lattice vectors.

!cod$
        call my(lat)
        tt = matmul(lat%o%ivectors,matmul(t,transpose(lat%o%ivectors)))
        call glean(thy(lat))
      end function

      function lat2r_vector(lat,v) result(tv)
!doc$ function lat2r(lat,v) result(tv)
        type(lattice_obj) :: lat
        real(double), dimension(:), intent(in) :: v
        real(double), dimension(3) :: tv
!       requires: v be dimension(3).
!       effects: Transforms v from the {a1,a2,a3} representation to the {x,y,z}
!                representation where a1, a2, and a3 are direct lattice vectors.

!cod$
        call my(lat)
        tv(1) = lat%o%vectors(1,1)*v(1) + lat%o%vectors(1,2)*v(2) + lat%o%vectors(1,3)*v(3)
        tv(2) = lat%o%vectors(2,1)*v(1) + lat%o%vectors(2,2)*v(2) + lat%o%vectors(2,3)*v(3)
        tv(3) = lat%o%vectors(3,1)*v(1) + lat%o%vectors(3,2)*v(2) + lat%o%vectors(3,3)*v(3)
        call glean(thy(lat))
      end function

      function lat2r_matrix(lat,m) result(tm)
!doc$ function lat2r(lat,m) result(tm)
        type(lattice_obj) :: lat
        real(double), dimension(:,:), intent(in) :: m
        real(double), dimension(3,3) :: tm
!       requires: m be dimension(3,3).
!       effects: Transforms m from the {a1,a2,a3} representation to the {x,y,z}
!                representation where a1, a2, and a3 are direct lattice vectors.

!cod$
        call my(lat)
        tm = matmul(lat%o%vectors,matmul(m,lat%o%ivectors))
        call glean(thy(lat))
      end function

      function lat2r_tensor(lat,t) result(tt)
!doc$ function lat2r_tensor(lat,t) result(tt)
        type(lattice_obj) :: lat
        real(double), dimension(:,:), intent(in) :: t
        real(double), dimension(3,3) :: tt
!       requires: t be dimension(3,3).
!       effects: Transforms t from the {a1,a2,a3} representation to the {x,y,z}
!                representation where a1, a2, and a3 are direct lattice vectors.

!cod$
        call my(lat)
        tt = matmul(lat%o%vectors,matmul(t,transpose(lat%o%vectors)))
        call glean(thy(lat))
      end function

      function f2lat_vector(lat,v) result(tv)
!doc$ function f2lat(lat,v) result(tv)
        type(lattice_obj) :: lat
        real(double), dimension(:), intent(in) :: v
        real(double), dimension(3) :: tv
!       requires: v be dimension(3).
!       effects: Transforms v from the {gx,gy,gz} representation to the {b1,b2,b3}
!                representation where b1, b2, and b3 are reciprocal lattice vectors.

!cod$ 
        call my(lat)
        tv(1) = (v(1)*lat%o%vectors(1,1) + v(2)*lat%o%vectors(2,1) + v(3)*lat%o%vectors(3,1))/two_pi
        tv(2) = (v(1)*lat%o%vectors(1,2) + v(2)*lat%o%vectors(2,2) + v(3)*lat%o%vectors(3,2))/two_pi
        tv(3) = (v(1)*lat%o%vectors(1,3) + v(2)*lat%o%vectors(2,3) + v(3)*lat%o%vectors(3,3))/two_pi
        call glean(thy(lat))
      end function

      function f2lat_matrix(lat,m) result(tm)
!doc$ function f2lat(lat,m) result(tm)
        type(lattice_obj) :: lat
        real(double), dimension(:,:), intent(in) :: m
        real(double), dimension(3,3) :: tm
!       requires: m be dimension(3,3).
!       effects: Transforms m from the {gx,gy,gz} representation to the {b1,b2,b3}
!                representation where b1, b2, and b3 are reciprocal lattice vectors.

!cod$ 
        call my(lat)
        tm = matmul(transpose(lat%o%vectors),matmul(m,transpose(lat%o%ivectors)))
        call glean(thy(lat))
      end function

      function lat2f_vector(lat,v) result(tv)
!doc$ function lat2f(lat,v) result(tv)
        type(lattice_obj) :: lat
        real(double), dimension(:), intent(in) :: v
        real(double), dimension(3) :: tv
!       requires: v be dimension(3).
!       effects: Transforms v from the {b1,b2,b3} representation to the {gx,gy,gz}
!                representation where b1, b2, and b3 are reciprocal lattice vectors.

!cod$ 
        call my(lat)
        tv(1) = (v(1)*lat%o%ivectors(1,1) + v(2)*lat%o%ivectors(2,1) + v(3)*lat%o%ivectors(3,1))*two_pi
        tv(2) = (v(1)*lat%o%ivectors(1,2) + v(2)*lat%o%ivectors(2,2) + v(3)*lat%o%ivectors(3,2))*two_pi
        tv(3) = (v(1)*lat%o%ivectors(1,3) + v(2)*lat%o%ivectors(2,3) + v(3)*lat%o%ivectors(3,3))*two_pi
        call glean(thy(lat))
      end function

      function lat2f_matrix(lat,m) result(tm)
!doc$ function lat2f(lat,m) result(tm)
        type(lattice_obj) :: lat
        real(double), dimension(3,3), intent(in) :: m
        real(double), dimension(3,3) :: tm
!       requires: m be dimension(3,3).
!       effects: Transforms m from the {b1,b2,b3} representation to the {gx,gy,gz}
!                representation where b1, b2, and b3 are reciprocal lattice vectors.

!cod$ 
        call my(lat)
        tm = matmul(transpose(lat%o%ivectors),matmul(m,transpose(lat%o%vectors)))
        call glean(thy(lat))
      end function

      subroutine save_lat(lat,prefix)
!doc$ subroutine save(lat,prefix)
        type(lattice_obj) :: lat
        character(*), intent(in) :: prefix
!       effects : Writes lattice information from file prefix//"lat".
!       errors : Format problems or IO problems.

!cod$
        type(file_obj) :: f
        integer :: ios
        call my(lat)
        call my(file(prefix//"lat"),f);
        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100
        if (i_access(f)) then 
          write(x_unit(f),*) lat%o%vectors/lat%o%lconstant
          write(x_unit(f),*) lat%o%lconstant
          close(unit=x_unit(f))
        end if
100     call glean(thy(lat))
        call glean(thy(f))
        if (error("Exit lattice_mod::save_lat")) continue
      end subroutine

      function cubic_cell(lat) result(cc)
!doc$ function cubic_cell(lat) result(cc)
        type(lattice_obj) :: lat
        logical :: cc
!       effects: Returns .true. if cell is a cube and .false. otherwise.

!cod$
        real(double), parameter :: tol_nbhd = 1.0e-9_double
        real(double), dimension(3) :: a1, a2, a3

        call my(lat)

        a1 = lat2r(lat,real((/1,0,0/),double))
        a2 = lat2r(lat,real((/0,1,0/),double))
        a3 = lat2r(lat,real((/0,0,1/),double))

        cc = .true.

        cc = cc .and. (norm(a1) .in. nbhd(norm(a2),tol_nbhd))
        cc = cc .and. (norm(a1) .in. nbhd(norm(a3),tol_nbhd))
        cc = cc .and. (norm(a2) .in. nbhd(norm(a3),tol_nbhd))

        cc = cc .and. (dot_product(a1,a2) .in. nbhd(0.0_double,tol_nbhd))
        cc = cc .and. (dot_product(a1,a3) .in. nbhd(0.0_double,tol_nbhd))
        cc = cc .and. (dot_product(a2,a3) .in. nbhd(0.0_double,tol_nbhd))

        call glean(thy(lat))

      end function

      subroutine wigner_seitz_vectors(lat,mode,vws)
!doc$ subroutine wigner_seitz_vectors(lat,mode,vws)
        type(lattice_obj) :: lat
        integer :: mode
        real(double), dimension(:,:), pointer :: vws
!       requires: vws be nullified or allocated
!       effects:   Generates the set of lattice (or ilattice) vectors, vws, with shortest lengths such
!                  that vws(i)*vws(j)/|vws(j)|**2 < 0.5 for all i,j.  This set defines the Wigner-Seitz
!                  cell sometimes referred to as the first Brillouin zone when vws is in reciprocal space.
!       errors: mode not equal to 1 (position space) or 2 (reciprocal space)

!cod$
        real(double), parameter :: tol_nbhd = 1.0e-8_double
        logical :: complete, inside
        integer :: i, i1, i2, i3, iv, ivt, inc, j, nv, nv1, nv2, nv3, nvt, nv_old
        real(double) :: proj, radius, volume, rxmt
        real(double), dimension(3) :: v1, v2, v3, rt, rxt
        real(double), dimension(:), allocatable :: vxm, vxmt
        real(double), dimension(:,:), allocatable :: v, vx, vt, vxt

        call my(lat)

        select case (mode)
        case (PS)
          v1 = lat2r(lat,real((/1,0,0/),double))
          v2 = lat2r(lat,real((/0,1,0/),double))
          v3 = lat2r(lat,real((/0,0,1/),double))
        case (RS)
          v1 = lat2f(lat,real((/1,0,0/),double))
          v2 = lat2f(lat,real((/0,1,0/),double))
          v3 = lat2f(lat,real((/0,0,1/),double))
        case default
          if (error(.true.,"ERROR: unrecognized lattice mode")) goto 100
        end select
        volume = abs(dot_product(v1,cross_product(v2,v3)))
        radius = 0.0_double

        nv_old = 0
        complete = .false.
        do while (.not.complete)

          radius = radius + ((3.0_double*volume)/(four_pi))**(1.0_double/3.0_double)
          nv1 = ceiling( radius*norm(cross_product(v2,v3))/volume )
          nv2 = ceiling( radius*norm(cross_product(v3,v1))/volume )
          nv3 = ceiling( radius*norm(cross_product(v1,v2))/volume )

          nv = (2*nv1 + 1)*(2*nv2 + 1)*(2*nv3 + 1)
          if (allocated( vt )) deallocate( vt ) ; allocate( vt(3,nv) )
          if (allocated( vxt )) deallocate( vxt ) ; allocate( vxt(3,nv) )
          if (allocated( vxmt )) deallocate( vxmt ) ; allocate( vxmt(nv) )

          nv = 0
          do i1 = -nv1,+nv1
            do i2 = -nv2,+nv2
              do i3 = -nv3,+nv3
                if ((i1 == 0) .and. (i2 == 0) .and. (i3 == 0)) cycle
                nv = nv + 1
                vt(:,nv) = real((/i1,i2,i3/),double)
                select case (mode)
                case (PS)
                  vxt(:,nv) = lat2r(lat,vt(:,nv))
                case (RS)
                  vxt(:,nv) = lat2f(lat,vt(:,nv))
                end select
                vxmt(nv) = norm(vxt(:,nv))
              end do
            end do
          end do

          inc = 1
          do
            inc = 3*inc + 1
            if (inc > nv) exit
          end do
          do
            inc = inc/3
            do i = inc+1,nv
              rt = vt(:,i)
              rxt = vxt(:,i)
              rxmt = vxmt(i)
              j = i
              do
                if (vxmt(j-inc) <= rxmt) exit
                vt(:,j) = vt(:,j-inc)
                vxt(:,j) = vxt(:,j-inc)
                vxmt(j) = vxmt(j-inc)
                j = j - inc
                if (j <= inc) exit
              end do
              vt(:,j) = rt
              vxt(:,j) = rxt
              vxmt(j) = rxmt
            end do
            if (inc <= 1) exit
          end do

          if (allocated( v )) deallocate( v ) ; allocate( v(3,nv) )
          if (allocated( vx )) deallocate( vx ) ; allocate( vx(3,nv) )
          if (allocated( vxm )) deallocate( vxm ) ; allocate( vxm(nv) )

          nvt = nv
          nv = 0
          do ivt = 1,nvt
            inside = .true.
            do iv = 1,nv
              proj = dot_product(vxt(:,ivt),vx(:,iv))/vxm(iv)
              if ( (abs(proj) .in. nbhd(vxmt(ivt),tol_nbhd)) .and. (vxmt(ivt) .in. nbhd(vxm(iv),tol_nbhd)) ) exit
              if (abs(proj/vxm(iv)) > (1.0_double - tol_nbhd)) then
                inside = .false.
                exit
              end if
            end do
            if (inside) then
              nv = nv + 1
              v(:,nv) = vt(:,ivt)
              vx(:,nv) = vxt(:,ivt)
              vxm(nv) = vxmt(ivt)
            end if
          end do

          complete = .true.
          if (nv > nv_old) complete = .false.
          nv_old = nv

        end do

        if (associated( vws )) deallocate( vws ) ; allocate( vws(3,nv) )
        vws = v(:,1:nv)

100     if (allocated( vxm )) deallocate( vxm )
        if (allocated( vxmt )) deallocate( vxmt )
        if (allocated( vx )) deallocate( vx )
        if (allocated( vxt )) deallocate( vxt )
        if (allocated( v )) deallocate( v )
        if (allocated( vt )) deallocate( vt )

        call glean(thy(lat))

        if (error("Exit lattice_mod::wigner_seitz_vectors")) continue

      end subroutine

      function maximum_sphere_radius(lat) result(r)
!doc$ function maximum_sphere_radius(lat) result(r)
        type(lattice_obj) :: lat
        real(double) :: r
!       effects: Returns the radius of the largest sphere that does not fit
!                completely inside the parallelpiped defined by lat%o%vectors.

!cod$
        real(double) :: b1m, b2m, b3m
        call my(lat)
        b1m = norm(lat2f(lat,real((/1,0,0/),double)))
        b2m = norm(lat2f(lat,real((/0,1,0/),double)))
        b3m = norm(lat2f(lat,real((/0,0,1/),double)))
        r = pi/max(b1m,b2m,b3m)
        call glean(thy(lat))
      end function

      subroutine diary_lat(lat)
!doc$ subroutine diary(lat)
        type(lattice_obj) :: lat
!       effects: Writes lattice information to the diary file.

!cod$
        call my(lat)
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Primitive lattice vectors:")')
          write(x_unit(diaryfile()),'(/,t6,"a1:",3f17.10)') lat%o%vectors(:,1)
          write(x_unit(diaryfile()),'(  t6,"a2:",3f17.10)') lat%o%vectors(:,2)
          write(x_unit(diaryfile()),'(  t6,"a3:",3f17.10)') lat%o%vectors(:,3)
          write(x_unit(diaryfile()),'(/,t6,"b1:",3f17.10)') two_pi*lat%o%ivectors(1,:)
          write(x_unit(diaryfile()),'(  t6,"b2:",3f17.10)') two_pi*lat%o%ivectors(2,:)
          write(x_unit(diaryfile()),'(  t6,"b3:",3f17.10)') two_pi*lat%o%ivectors(3,:)
          write(x_unit(diaryfile()),'(/,t6,"cell volume = ",f0.5)') x_cell_volume(lat)
        end if
        call glean(thy(lat))
      end subroutine

      subroutine write_restart_lat(lat,nrestf)
!doc$ subroutine write_restart(lat,nrestf)
        type(lattice_obj) :: lat
        type(tagio_obj) :: nrestf
!       requires: nrestf be positioned inside the CRYSTAL block.
!       modifies: nrestf
!       effects: Writes lat information to nrestf.

!cod$
        integer(long) :: dsize, iosl, ndata

        call my(lat)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the LATTICE block
          call startblock(nrestf,"LATTICE")

          ! write the lattice constant
          call writetag(nrestf,"LATTICE_CONSTANT")
          dsize = sizeof_double ; ndata = 1
          call writef(lat%o%lconstant,dsize,ndata,x_tagfd(nrestf),iosl)

          ! write the real-space lattice vectors
          call writetag(nrestf,"REAL-SPACE_VECTORS")
          dsize = sizeof_double ; ndata = 9
          call writef(lat%o%vectors,dsize,ndata,x_tagfd(nrestf),iosl)

          ! write the reciprocal-space lattice vectors
          call writetag(nrestf,"RECIPROCAL-SPACE_VECTORS")
          dsize = sizeof_double ; ndata = 9
          call writef(lat%o%ivectors,dsize,ndata,x_tagfd(nrestf),iosl)

          ! end the LATTICE block
          call endblock(nrestf)

        end if

        call glean(thy(lat))
        call glean(thy(nrestf))

        if (error("Exit lattice_mod::write_restart_lat")) continue

      end subroutine

    end module

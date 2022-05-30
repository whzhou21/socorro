! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module atomic_potential_paw_mod
!doc$ module atomic_potential_paw_mod

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
      use axc_mod
      use diary_mod
      use math_mod
      use symmetry_mod
      use paw_data_mod
      use atomic_operators_paw_mod
      use atomic_density_paw_mod

!     One datatype is available here: type(atomic_potential_paw_obj).

!cod$
      implicit none
      private

      type :: gp_mat
        complex(double), dimension(:,:), pointer :: mat
      end type

      type :: indices
         integer :: ai
         integer :: ti
         integer :: tvi
      end type

      type :: atomic_potential_paw_rep
        integer :: ref
        type(ghost) :: g
        type(atomic_operators_paw_obj) :: ao            ! atomic operators object
        type(gp_mat), dimension(:), pointer :: dij      ! atomic potential matrix for an sgroup
      end type 

      type, public :: atomic_potential_paw_obj
        private
        integer :: ref
        type(atomic_potential_paw_rep), pointer :: o
      end type

!doc$
      public :: atomic_potential_paw
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_atomic_operators
      public :: atomic_hamiltonian
      public :: atomic_overlap
      public :: extract_potential
      public :: insert_potential
      public :: write_restart

!cod$
      interface atomic_potential_paw
        module procedure constructor_ap
      end interface
      interface update
        module procedure update_ap
      end interface
      interface my
        module procedure my_ap, my_new_ap
      end interface
      interface thy
        module procedure thy_ap
      end interface
      interface glean
        module procedure glean_ap
      end interface
      interface bequeath
        module procedure bequeath_ap
      end interface
      interface assignment(=)
        module procedure assign_ap
      end interface
      interface x_ref
        module procedure ap_ref
      end interface
      interface x_ghost
        module procedure ap_ghost
      end interface
      interface x_atomic_operators
        module procedure ap_atomic_operators
      end interface
      interface atomic_hamiltonian
        module procedure atomic_hamiltonian_1d_ap, atomic_hamiltonian_2d_ap
      end interface
      interface atomic_overlap
        module procedure atomic_overlap_1d_ap, atomic_overlap_2d_ap
      end interface
      interface extract_potential
        module procedure extract_potential_ap
      end interface
      interface insert_potential
        module procedure insert_potential_ap
      end interface
      interface write_restart
        module procedure write_restart_ap
      end interface

      contains

! public routines

      function constructor_ap(ad,hap,restf) result(ap)
!doc$ function atomic_potential_paw(ad,hap,restf) result(ap)
        type(atomic_density_paw_obj) :: ad
        type(grid_obj), optional :: hap
        type(tagio_obj), optional :: restf
        type(atomic_potential_paw_obj) :: ap
!       requires: restf pointer be positioned inside the FIELDS block. Either hap or restf be absent.
!       effects: Constructs a new ap.
!       errors: Tags not found. Passes errors.

!cod$
        logical :: np_ok
        character(1) :: tios
        integer :: i, i1, i2, ia, msg, n, na, nm, np, nsg
        integer, dimension(:), allocatable :: nps
        integer(long) :: dsize, iosl, ndata, s4
        integer(long), dimension(:), allocatable :: v4
        complex(double), dimension(:), allocatable :: c1

        call my(ad)
        if (present(hap)) call my(hap)
        if (present(restf)) call my(restf)

        ap%ref = 0
        allocate( ap%o )
        ap%o%ref = 0
        ap%o%g = x_ghost()

        call my(x_atomic_operators(ad),ap%o%ao)

        nullify( ap%o%dij )

        if (present(restf)) then

          nsg = mpi_nsgroups()
          msg = mpi_mysgroup()

          ! Open the ATOMIC_POTENTIAL block
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_POTENTIAL")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_POTENTIAL block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)

          ! Find the PAW tag
          if (i_access(restf)) tios = findfirsttag(restf,"PAW")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: PAW tag was not found")) goto 100

          ! Read the number of atoms
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_ATOMS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_ATOMS tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long
            ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            na = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,na)
          if (error(na /= x_n_atoms(ap%o%ao),"ERROR: different numbers of atoms")) goto 100

          ! Read the number of projectors
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_PROJECTORS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_PROJECTORS tag was not found")) goto 100
          allocate( nps(na), v4(na) )
          if (i_access(restf)) then
            dsize = sizeof_long
            ndata = na
            call readf(v4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            nps = v4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,nps)
          np_ok = .true.
          do ia = 1,na
            if (nps(ia) /= x_n_atom_projectors(ap%o%ao,ia)) then
              np_ok = .false.
              exit
            end if
          end do
          if (error(.not.np_ok,"ERROR: different numbers of projectors")) goto 100

          ! Set the file pointer at the beginning of the dij matrix elements
          if (i_access(restf)) tios = findfirsttag(restf,"DIJ_MATRIX_ELEMENTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: DIJ_MATRIX_ELEMENTS tag was not found")) goto 100

          ! Allocate space for the data
          n = 0
          do ia = 1,size(nps)
            n = n + nps(ia)*nps(ia)
          end do
          allocate( c1(n) )

          ! Read and distribute the data for spin group 1
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 2*n
            call readf(c1,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,c1)
          select case (msg)
          case (1)
            allocate( ap%o%dij(na) )
            i = 0
            do ia = 1,size(ap%o%dij)
              nm = nps(ia)
              allocate( ap%o%dij(ia)%mat(nm,nm) )
              do i2 = 1,nm
              do i1 = 1,nm
                i = i + 1
                ap%o%dij(ia)%mat(i1,i2) = c1(i)
              end do
              end do
            end do
          end select

          ! Read and distribute the data for spin group 2
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 2*n
            call readf(c1,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,c1)
          select case (msg)
          case (2)
            allocate( ap%o%dij(na) )
            i = 0
            do ia = 1,size(ap%o%dij)
              nm = nps(ia)
              allocate( ap%o%dij(ia)%mat(nm,nm) )
              do i2 = 1,nm
              do i1 = 1,nm
                i = i + 1
                ap%o%dij(ia)%mat(i1,i2) = c1(i)
              end do
              end do
            end do
          end select

          ! Close the ATOMIC_POTENTIAL block
100       if (i_access(restf)) call closeblock(restf)

        else

          ! Form the dij matrix elements
          call form_dij_i(ap%o,ad,hap)

        end if

200     if (allocated( nps )) deallocate( nps )
        if (allocated( v4 )) deallocate( v4 )
        if (allocated( c1 )) deallocate( c1 )

        call glean(thy(ad))
        if (present(hap)) call glean(thy(hap))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_potential_paw_mod::constructor_ap")) continue

      end function 

      subroutine update_ap(ap,ad,hap)
!doc$ subroutine update(ap,ad,hap)
        type(atomic_potential_paw_obj) :: ap
        type(atomic_density_paw_obj) :: ad
        type(grid_obj) :: hap ! grid hartree potential
!       effects: Updates ap.

!cod$
        call my(ap)
        call my(ad)
        call my(hap)

        call own_i(ap)
        ap%o%g = x_ghost()
        ap%o%ao = x_atomic_operators(ad)
        call form_dij_i(ap%o,ad,hap)

100     call glean(thy(ap))
        call glean(thy(ad))
        call glean(thy(hap))
        
        if (error("Exit atomic_potential_paw_mod::update_ap")) continue

      end subroutine

      subroutine my_ap(ap)
!doc$ subroutine my(ap)
        type(atomic_potential_paw_obj) :: ap

!cod$
        ap%ref = ap%ref + 1
        ap%o%ref = ap%o%ref + 1
      end subroutine
  
      subroutine my_new_ap(api,ap)
!doc$ subroutine my(api,ap)
        type(atomic_potential_paw_obj) :: api, ap

!cod$
        ap%ref = 1
        ap%o => api%o
        ap%o%ref = ap%o%ref + 1
      end subroutine

      function thy_ap(ap) result(apo)
!doc$ function thy(ap) result(apo)
        type(atomic_potential_paw_obj) :: ap, apo
        
!cod$
        ap%ref = ap%ref - 1
        ap%o%ref = ap%o%ref - 1
        apo%ref = ap%ref
        apo%o => ap%o
      end function 

      subroutine glean_ap(ap)
!doc$ subroutine glean(ap)
        type(atomic_potential_paw_obj) :: ap

!cod$
        integer :: i
        if (ap%o%ref < 1) then
          call glean(thy(ap%o%ao))
          if (associated( ap%o%dij )) then
            do i = 1,size(ap%o%dij)
              if (associated( ap%o%dij(i)%mat )) deallocate( ap%o%dij(i)%mat )
            end do
            deallocate( ap%o%dij )
          end if
          deallocate( ap%o )
        end if
    
      end subroutine

      subroutine bequeath_ap(ap)
!doc$ subroutine bequeath(ap)
        type(atomic_potential_paw_obj) :: ap

!cod$
        continue
      end subroutine

      subroutine assign_ap(ap,ap2)
!doc$ subroutine assignment(=)(ap,ap2)
        type(atomic_potential_paw_obj), intent(inout) :: ap
        type(atomic_potential_paw_obj), intent(in) :: ap2
        
!cod$
        type(atomic_potential_paw_obj) :: apt
        call my(ap2)
        apt%o => ap%o
        ap%o%ref = ap%o%ref - ap%ref
        ap%o => ap2%o
        ap%o%ref = ap%o%ref + ap%ref
        call glean(apt)
        call glean(thy(ap2))
      end subroutine

      function ap_ref(ap) result(r)
!doc$ function x_ref(ap) result(r)
        type(atomic_potential_paw_obj) :: ap
        integer, dimension(2) :: r
!       effects: Returns the reference counts of ap.

!cod$
        r(1) = ap%ref
        r(2) = ap%o%ref
        call glean(ap)
      end function 

      function ap_ghost(ap) result(g)
!doc$ function x_ghost(ap) result(g)
        type(atomic_potential_paw_obj) :: ap
        type(ghost) :: g
!       effects: Returns the ghost of ap.

!cod$ 
        call my(ap)
        g = ap%o%g
        call glean(thy(ap))
      end function 

      function ap_atomic_operators(ap) result(ao)
!doc$ function x_atomic_operators(ap) result(ao)
        type(atomic_potential_paw_obj) :: ap
        type(atomic_operators_paw_obj) :: ao
!       effects: Returns ap%o%ao.
  
!cod$
        call my(ap)
        call my(ap%o%ao,ao)
        call glean(thy(ap))
        call bequeath(thy(ao))
      end function

      subroutine atomic_hamiltonian_1d_ap(ap,pdots)
!doc$ subroutine atomic_hamiltonian(ap,pdots)
        type(atomic_potential_paw_obj) :: ap
        complex(double), dimension(:), intent(inout) :: pdots
!       modifies: pdots
!       effects: ap%o%dij()%mat*pdots --> pdots
!       errors: Passes errors.

!cod$
        character(1), parameter :: trans = 'n'
        integer :: bi, ia, nm
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:), allocatable :: d_pdots

        call my(ap)

        allocate( d_pdots(size(pdots)) )
        d_pdots = (0.0_double,0.0_double)
        do ia = 1,size(ap%o%dij)
          nm = size(ap%o%dij(ia)%mat,1)
          bi = x_atom_base(ap%o%ao,ia)
          call zgemv(trans,nm,nm,alpha,ap%o%dij(ia)%mat,nm,pdots(bi),1,beta,d_pdots(bi),1)
        end do
        pdots = d_pdots

100     if (allocated( d_pdots )) deallocate( d_pdots )

        call glean(thy(ap))

        if (error("Exit atomic_potential_paw_mod::atomic_hamiltonian_1d_ap")) continue

      end subroutine

      subroutine atomic_hamiltonian_2d_ap(ap,pdots)
!doc$ subroutine atomic_hamiltonian(ap,pdots)
        type(atomic_potential_paw_obj) :: ap
        complex(double), dimension(:,:), intent(inout) :: pdots
!       modifies: pdots
!       effects: ap%o%dij()%mat*pdots --> pdots
!       errors: Passes errors.

!cod$
        character(1), parameter :: transa = 'n', transb = 'n'
        integer :: bi, ia, nb, nm, np
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: d_pdots

        call my(ap)

        np = size(pdots,1)
        nb = size(pdots,2)
        allocate( d_pdots(np,nb) )
        d_pdots = (0.0_double,0.0_double)
        do ia = 1,size(ap%o%dij)
          nm = size(ap%o%dij(ia)%mat,1)
          bi = x_atom_base(ap%o%ao,ia)
          call zgemm(transa,transb,nm,nb,nm,alpha,ap%o%dij(ia)%mat,nm,pdots(bi,1),np,beta,d_pdots(bi,1),np)
        end do
        pdots = d_pdots

100     if (allocated( d_pdots )) deallocate( d_pdots )

        call glean(thy(ap))

        if (error("Exit atomic_potential_paw_mod::atomic_hamiltonian_2d_ap")) continue

      end subroutine

      subroutine atomic_overlap_1d_ap(ap,pdots)
!doc$ subroutine atomic_overlap(ap,pdots)
        type(atomic_potential_paw_obj) :: ap
        complex(double), dimension(:), intent(inout) :: pdots
!       modifies: pdots
!       effects: ap%o%ao%o%r_oij()%mat*pdots --> pdots
!       errors: Passes errors.

!cod$
        character(1), parameter :: trans = 'n'
        integer :: bi, ia, nm
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:), allocatable :: o_pdots
        complex(double), dimension(:,:), pointer :: a_oij

        call my(ap)

        allocate( o_pdots(size(pdots)) )
        o_pdots = (0.0_double,0.0_double)
        do ia = 1,size(ap%o%dij)
          call atom_r_oij(ap%o%ao,ia,a_oij) ; if (error()) goto 100
          nm = size(a_oij,1)
          bi = x_atom_base(ap%o%ao,ia)
          call zgemv(trans,nm,nm,alpha,a_oij,nm,pdots(bi),1,beta,o_pdots(bi),1)
        end do
        pdots = o_pdots

100     if (allocated( o_pdots )) deallocate( o_pdots )
        nullify( a_oij )

        call glean(thy(ap))

        if (error("Exit atomic_potential_paw_mod::atomic_overlap_1d_ap")) continue

      end subroutine

      subroutine atomic_overlap_2d_ap(ap,pdots)
!doc$ subroutine atomic_overlap(ap,pdots)
        type(atomic_potential_paw_obj) :: ap
        complex(double), dimension(:,:), intent(inout) :: pdots
!       modifies: pdots
!       effects: ap%o%ao%o%r_oij()%mat*pdots --> pdots
!       errors: Passes errors.

!cod$
        character(1), parameter :: transa = 'n', transb = 'n'
        integer :: bi, ia, nb, nm, np
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: o_pdots
        complex(double), dimension(:,:), pointer :: a_oij

        call my(ap)

        np = size(pdots,1)
        nb = size(pdots,2)
        allocate( o_pdots(np,nb) )
        o_pdots = (0.0_double,0.0_double)
        do ia = 1,size(ap%o%dij)
          call atom_r_oij(ap%o%ao,ia,a_oij) ; if (error()) goto 100
          nm = size(a_oij,1)
          bi = x_atom_base(ap%o%ao,ia)
          call zgemm(transa,transb,nm,nb,nm,alpha,a_oij,nm,pdots(bi,1),np,beta,o_pdots(bi,1),np)
        end do
        pdots = o_pdots

100     if (allocated( o_pdots )) deallocate( o_pdots )
        nullify( a_oij )

        call glean(thy(ap))

        if (error("Exit atomic_potential_paw_mod::atomic_overlap_2d_ap")) continue

      end subroutine

      subroutine extract_potential_ap(ap,c1d)
!doc$ subroutine extract_potential(ap,c1d)
        type(atomic_potential_paw_obj) :: ap
        complex(double), dimension(:), pointer :: c1d
!       requires: c1d be nullified or associated.
!       modifies: c1d
!       effects: Reshapes ap%o%dij into c1d.

!cod$
        integer :: i, i1, i2, ia, n1d

        call my(ap)

        n1d = 0
        do ia = 1,size(ap%o%dij)
          n1d = n1d + size(ap%o%dij(ia)%mat)
        end do
        if (associated( c1d )) deallocate( c1d ) ; allocate( c1d(n1d) )

        i = 0
        do ia = 1,size(ap%o%dij)
          do i2 = 1,size(ap%o%dij(ia)%mat,2)
          do i1 = 1,size(ap%o%dij(ia)%mat,1)
            i = i + 1
            c1d(i) = ap%o%dij(ia)%mat(i1,i2)
          end do
          end do
        end do

        call glean(thy(ap))

        if (error("Exit atomic_potential_paw_mod::extract_potential_ap")) continue

      end subroutine

      subroutine insert_potential_ap(c1d,ap)
!doc$ subroutine insert_potential(c1d,ap)
        complex(double), dimension(:), pointer :: c1d
        type(atomic_potential_paw_obj) :: ap
!       modifies: ap
!       requires: size(c1d) = size(ap%o%dij(*)%mat)
!       effects: Reshapes c1d into ap%o%dij.

!cod$
        integer :: i, i1, i2, ia

        call my(ap)

        call own_i(ap)
        ap%o%g = x_ghost()
        i = 0
        do ia = 1,size(ap%o%dij)
          do i2 = 1,size(ap%o%dij(ia)%mat,2)
          do i1 = 1,size(ap%o%dij(ia)%mat,1)
            i = i + 1
            ap%o%dij(ia)%mat(i1,i2) = c1d(i)
          end do
          end do
        end do

        call glean(thy(ap))

        if (error("Exit atomic_potential_paw_mod::insert_potential_ap")) continue

      end subroutine

      subroutine write_restart_ap(ap,nrestf)
!doc$ subroutine write_restart(ap,nrestf)
        type(atomic_potential_paw_obj) :: ap
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ap restart information to nrestf.

!cod$
        integer :: i, i1, i2, ia, msg, n, na, nsg
        integer(long) :: dsize, iosl, ndata, s4
        integer(long), dimension(:), allocatable :: v4
        complex(double), dimension(:), allocatable :: c1, c2

        call my(ap)
        call my(nrestf)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        ! Start the ATOMIC_POTENTIAL block
        if (i_access(nrestf)) call startblock(nrestf,"ATOMIC_POTENTIAL")

        ! Write the PAW tag
        if (i_access(nrestf)) call writetag(nrestf,"PAW")

        ! write the number of atoms
        if (i_access(nrestf)) then
          call writetag(nrestf,"NUMBER_OF_ATOMS")
          na = x_n_atoms(ap%o%ao)
          dsize = sizeof_long
          ndata = 1
          s4 = na
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Write the number of projectors
        if (i_access(nrestf)) then
          call writetag(nrestf,"NUMBER_OF_PROJECTORS")
          dsize = sizeof_long
          ndata = na
          allocate( v4(na) )
          do ia = 1,na
            v4(ia) = x_n_atom_projectors(ap%o%ao,ia)
          end do
          call writef(v4,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Reshape the matrix elements into a 1-D array
        n = 0
        do ia = 1,size(ap%o%dij)
          n = n + size(ap%o%dij(ia)%mat)
        end do
        allocate( c1(n) )
        c1 = (0.0_double,0.0_double)
        i = 0
        do ia = 1,size(ap%o%dij)
          do i2 = 1,size(ap%o%dij(ia)%mat,2)
          do i1 = 1,size(ap%o%dij(ia)%mat,1)
            i = i + 1
            c1(i) = ap%o%dij(ia)%mat(i1,i2)
          end do
          end do
        end do

        ! Write the spin group 1 data (residing on the FILE_SCOPE rank 0 process)
        if (i_access(nrestf)) then
          call writetag(nrestf,"DIJ_MATRIX_ELEMENTS")
          dsize = sizeof_double
          ndata = 2*n
          call writef(c1,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Reduce the spin group 2 data to the FILE_SCOPE rank 0 process and write
        select case (nsg)
        case (2)
          select case (msg)
          case (1)
            c1 = (0.0_double,0.0_double)
          end select
          allocate( c2(n) )
          call xcomm_reduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
          if (i_access(nrestf)) then
            dsize = sizeof_double
            ndata = 2*n
            call writef(c2,dsize,ndata,x_tagfd(nrestf),iosl)
          end if
        end select

        ! End the ATOMIC_POTENTIAL block
        if (i_access(nrestf)) call endblock(nrestf)

100     if (allocated( v4 )) deallocate( v4 )
        if (allocated( c1 )) deallocate( c1 )
        if (allocated( c2 )) deallocate( c2 )

        call glean(thy(ap))
        call glean(thy(nrestf))

        if (error("Exit atomic_potential_paw_mod::write_restart_ap")) continue
      
      end subroutine

! private routines

      subroutine form_dij_i(apr,ad,hap) 
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(grid_obj) :: hap
!       notes: Unless an "atomic_symmetry  off" line is present in argvf, this routine explicitly computes dij
!              contributions only for atoms that are unique with respect to symmetry. dij contributions from
!              symmetry-related atoms are then accumulated via symmetry operations.

        integer :: ia, it, iu, np, nu
        complex(double), dimension(:,:), pointer :: t_r2c, t_c2r
        type(gp_mat), dimension(:), pointer :: dij

        call my(ad)
        call my(hap)

        if (associated( apr%dij )) then
          do ia = 1,size(apr%dij)
            if (associated( apr%dij(ia)%mat )) deallocate( apr%dij(ia)%mat )
          end do
          deallocate( apr%dij )
        end if

        ! Create a temporary array for the dij contributions.
        nu = x_n_unique_atoms(apr%ao)
        allocate( dij(nu) )
        do iu = 1,nu
          ia = x_unique_atom(apr%ao,iu)
          np = x_n_atom_projectors(apr%ao,ia)
          allocate( dij(iu)%mat(np,np) )
          dij(iu)%mat = (0.0_double,0.0_double)
        end do

        ! Compute the dij contributions.
        call diagonal_terms_i(apr,dij)         ; if (error()) goto 100  
        call hartree_terms_i(apr,ad,dij)       ; if (error()) goto 100
        call xc_terms_i(apr,ad,dij)            ; if (error()) goto 100
        call hat_terms_i(apr,ad,dij)           ; if (error()) goto 100
        call multipole_terms_i(apr,ad,hap,dij) ; if (error()) goto 100

        ! Distribute the dij contributions to non-symmetry-unique atoms.
        call distribute_dij_i(apr,dij)

        ! Transform apr%dij to real-valued spherical harmonics.
        do ia = 1,size(apr%dij)
          it = x_atom_type(apr%ao,ia)
          call type_r2c(apr%ao,it,t_r2c)
          call type_c2r(apr%ao,it,t_c2r)
          apr%dij(ia)%mat = matmul(t_c2r,matmul(apr%dij(ia)%mat,t_r2c))
        end do

100     if (associated( dij )) then
          do ia = 1,size(dij)
            if (associated( dij(ia)%mat )) deallocate( dij(ia)%mat )
          end do
          deallocate( dij )
        end if
        nullify( t_r2c )
        nullify( t_c2r )

        call glean(thy(ad))
        call glean(thy(hap))

        if (error("Exit atomic_potential_paw_mod::form_dij_i")) continue

      end subroutine 

      subroutine diagonal_terms_i(apr,dij) ! PRB 55, 2005 (1997): Equations A10, A11, A12, and A31(1,2,3).
        type(atomic_potential_paw_rep) :: apr
        type(gp_mat), dimension(:), pointer :: dij
        
        integer :: i, ia, it, m
        integer :: bi, bj, li, lj, oi, oj, omi, omj
        integer :: iu, nu, nu_proc, first_u, last_u
        real(double) :: fn
        type(paw_data_obj), pointer :: pawd 
        
        nu = x_n_unique_atoms(apr%ao)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u
          ia = x_unique_atom(apr%ao,iu)
          it = x_unique_type(apr%ao,iu)
          call type_projector_data(apr%ao,it,pawd)
          i = 0 
          do bi = 1,x_n_basis(apr%ao,it)
            li = x_basis_l(apr%ao,it,bi)
            do bj = 1,bi
              lj = x_basis_l(apr%ao,it,bj) 
              if (lj == li) then 
                oi = x_basis_offset(apr%ao,it,bi)
                oj = x_basis_offset(apr%ao,it,bj)
                i = i + 1 
                fn = x_kinetic(pawd,i) + x_v_ion(pawd,i) 
                do m = -lj,+lj 
                  omi = oi + m
                  omj = oj + m
                  dij(iu)%mat(omi,omj) = dij(iu)%mat(omi,omj) + cmplx(fn,0,double)
                  if (bj /= bi) then 
                    dij(iu)%mat(omj,omi) = dij(iu)%mat(omj,omi) + cmplx(fn,0,double)
                  end if
                end do
              end if
            end do
          end do
        end do

        nullify( pawd )

        if (error("Exit atomic_potential_paw_mod::diagonal_terms_i")) continue

      end subroutine

      subroutine hartree_terms_i(apr,ad,dij) ! PRB 55, 2005 (1997): Equations A26 and A31 (6th term).
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(gp_mat), dimension(:), pointer :: dij
 
        integer :: i, ia, it
        integer :: bi, bj, bk, bl, mi, mj, mk, ml, omi, omj, omk, oml
        integer :: iu, nu, nu_proc, first_u, last_u
        complex(double) :: c1
        complex(double), dimension(:,:), pointer :: a_wijs
        type(paw_data_obj), pointer :: pawd

        call my(ad)
 
        nu = x_n_unique_atoms(apr%ao)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u
          ia = x_unique_atom(apr%ao,iu)
          it = x_unique_type(apr%ao,iu)
          call atom_wijs(ad,ia,a_wijs)
          call type_projector_data(apr%ao,it,pawd)
          do i = 1,x_cijkl_size(pawd)
            call cijkl_decode(pawd,i,bi,bj,bk,bl,mi,mj,mk)
            ml = mi - mj + mk
            omi = x_basis_offset(apr%ao,it,bi) + mi
            omj = x_basis_offset(apr%ao,it,bj) + mj
            omk = x_basis_offset(apr%ao,it,bk) + mk
            oml = x_basis_offset(apr%ao,it,bl) + ml
            c1 = x_cijkl(pawd,i)*a_wijs(omk,oml)
            dij(iu)%mat(omi,omj) = dij(iu)%mat(omi,omj) +  c1            
          end do
        end do

        nullify( a_wijs )
        nullify( pawd )

        call glean(thy(ad))

        if (error("Exit atomic_potential_paw_mod::hartree_terms_i")) continue

      end subroutine

      subroutine xc_terms_i(apr,ad,dij)
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(gp_mat), dimension(:), pointer :: dij

        call my(ad)

        if (uses_gradient(x_axc(apr%ao))) then
          call accum_xc_grad_i(apr,ad,dij)
        else
          call accum_xc_i(apr,ad,dij)
        end if

        call glean(thy(ad))

        if (error("Exit atomic_potential_paw_mod::xc_terms_i")) continue

      end subroutine

      subroutine accum_xc_i(apr,ad,dij)
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(gp_mat), dimension(:), pointer :: dij

        integer :: ia, it, iu, iv, nb, nr, nv
        integer :: bi, bj, li, lj, mi, mj, oi, oj, osi, osj
        real(double) :: fn, spin_factor, vwt
        real(double), dimension(:), allocatable :: r2dfxc_dn, r2dfxct_dn
        real(double), dimension(:), pointer :: n, del_n, lap_n, dfxc_dn, dfxc_dg, dfxc_dl
        real(double), dimension(:), pointer :: nt, del_nt, lap_nt, dfxct_dn, dfxct_dg, dfxct_dl
        real(double), dimension(:), pointer :: r2, wt, nc, nct
        real(double), dimension(:,:,:), pointer :: phir2ij, tphir2ij
        complex(double) :: cij
        complex(double), dimension(:,:), pointer :: wij, ylmij
        type(paw_data_obj), pointer :: pawd

        call my(ad)

        nullify( n, del_n, lap_n, dfxc_dn, dfxc_dg, dfxc_dl )
        nullify( nt, del_nt, lap_nt, dfxct_dn, dfxct_dg, dfxct_dl )

        allocate( n(0), dfxc_dn(0), r2dfxc_dn(0) )
        allocate( nt(0), dfxct_dn(0), r2dfxct_dn(0) )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        ! Compute dij contributions for symmetry-unique atoms.
        nv = x_n_vectors(apr%ao)
        do iv = 1,nv

          iu = x_vector_unique_atom(apr%ao,iv)
          ia = x_vector_atom(apr%ao,iv)
          it = x_vector_type(apr%ao,iv)

          call atom_wij(ad,ia,wij)
          call type_projector_data(apr%ao,it,pawd)

          nr = x_grid_size(pawd)
          if (size(n) /= nr) then
            deallocate( n, dfxc_dn, r2dfxc_dn ) ; allocate( n(nr), dfxc_dn(nr), r2dfxc_dn(nr) )
            deallocate( nt, dfxct_dn, r2dfxct_dn ) ; allocate( nt(nr), dfxct_dn(nr), r2dfxct_dn(nr) )
          end if
          n = 0.0_double
          nt = 0.0_double

          call type_r2(pawd,r2)
          call type_wt(pawd,wt)
          call type_core_density(pawd,nc)
          call type_coretail_density(pawd,nct)
          call type_phir2ij(pawd,phir2ij)
          call type_tphir2ij(pawd,tphir2ij)

          nb = x_n_basis(apr%ao,it)
          call vector_ylmij(apr%ao,iv,ylmij)
          vwt = x_vector_weight(apr%ao,iv)

          do bj = 1,nb
            lj = x_basis_l(apr%ao,it,bj)
            oj = x_basis_offset(apr%ao,it,bj)
            osj = lj**2 + lj + 1
            do bi = 1,nb
              li = x_basis_l(apr%ao,it,bi)
              oi = x_basis_offset(apr%ao,it,bi)
              osi = li**2 + li + 1
              cij = (0.0_double,0.0_double)
              do mj = -lj,+lj
              do mi = -li,+li
                cij = cij + wij(oi+mi,oj+mj)*ylmij(osi+mi,osj+mj)
              end do
              end do
              n = n + real(cij,double)*phir2ij(:,bi,bj)
              nt = nt + real(cij,double)*tphir2ij(:,bi,bj)
            end do
          end do
          n = n + spin_factor*nc
          nt = nt + spin_factor*nct

          call xc_derivatives(x_axc(apr%ao),n,del_n,lap_n,dfxc_dn,dfxc_dg,dfxc_dl)
          call xc_derivatives(x_axc(apr%ao),nt,del_nt,lap_nt,dfxct_dn,dfxct_dg,dfxct_dl)

          r2dfxc_dn = (r2*dfxc_dn*wt)*vwt
          r2dfxc_dn(nr) = 0.0_double
          r2dfxct_dn = (r2*dfxct_dn*wt)*vwt
          r2dfxct_dn(nr) = 0.0_double

          do bj = 1,nb
            lj = x_basis_l(apr%ao,it,bj)
            oj = x_basis_offset(apr%ao,it,bj)
            osj = lj**2 + lj + 1
            do bi = 1,nb
              li = x_basis_l(apr%ao,it,bi)
              oi = x_basis_offset(apr%ao,it,bi)
              osi = li**2 + li + 1
              fn = sum(r2dfxc_dn*phir2ij(:,bi,bj) - r2dfxct_dn*tphir2ij(:,bi,bj))
              do mj = -lj,+lj
              do mi = -li,+li
                dij(iu)%mat(oi+mi,oj+mj) = dij(iu)%mat(oi+mi,oj+mj) + fn*ylmij(osi+mi,osj+mj)
              end do
              end do
            end do
          end do

        end do

        if (associated( n )) then
          deallocate( n, dfxc_dn, r2dfxc_dn )
          deallocate( nt, dfxct_dn, r2dfxct_dn )
        end if
        nullify( r2 )
        nullify( wt)
        nullify( nc )
        nullify( nct )
        nullify( phir2ij )
        nullify( tphir2ij )
        nullify( wij )
        nullify( ylmij )
        nullify( pawd )

        call glean(thy(ad))

        if (error("Exit atomic_potential_paw_mod::accum_xc_i")) continue

      end subroutine 

      subroutine accum_xc_grad_i(apr,ad,dij)
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(gp_mat), dimension(:), pointer :: dij

        integer :: ia, ir, it, iu, iv, nb, nr, nv
        integer :: bi, bj, li, lj, mi, mj, oi, oj, osi, osj
        real(double) :: fn_r, fn_t, fn_p, spin_factor, vwt
        real(double), dimension(:), allocatable :: dn_dr, dn_dt, dn_dp, r2dfxc_dn, r2dfxc_dr, dfxc_dt, dfxc_dp
        real(double), dimension(:), allocatable :: dnt_dr, dnt_dt, dnt_dp, r2dfxct_dn, r2dfxct_dr, dfxct_dt, dfxct_dp
        real(double), dimension(:), pointer :: n, del_n, lap_n, dfxc_dn, dfxc_dg, dfxc_dl
        real(double), dimension(:), pointer :: nt, del_nt, lap_nt, dfxct_dn, dfxct_dg, dfxct_dl
        real(double), dimension(:), pointer :: r2, wt, nc, dnc_dr, nct, dnct_dr
        real(double), dimension(:,:,:), pointer :: phir2ij, dphir2ij_dr, tphir2ij, dtphir2ij_dr
        complex(double) :: cij, dcij_dt, dcij_dp
        complex(double), dimension(:,:), pointer :: wij, ylmij, dylmij_dt, dylmij_dp
        type(paw_data_obj), pointer :: pawd

        call my(ad)

        nullify( n, del_n, lap_n, dfxc_dn, dfxc_dg, dfxc_dl )
        nullify( nt, del_nt, lap_nt, dfxct_dn, dfxct_dg, dfxct_dl )

        allocate( n(0), dn_dr(0), dn_dt(0), dn_dp(0), del_n(0) )
        allocate( dfxc_dn(0), dfxc_dg(0) )
        allocate( r2dfxc_dn(0), r2dfxc_dr(0), dfxc_dt(0), dfxc_dp(0) )
        allocate( nt(0), dnt_dr(0), dnt_dt(0), dnt_dp(0), del_nt(0) )
        allocate( dfxct_dn(0), dfxct_dg(0) )
        allocate( r2dfxct_dn(0), r2dfxct_dr(0), dfxct_dt(0), dfxct_dp(0) )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        ! Compute dij contributions for symmetry-unique atoms.
        nv = x_n_vectors(apr%ao)
        do iv = 1,nv

          iu = x_vector_unique_atom(apr%ao,iv)
          ia = x_vector_atom(apr%ao,iv)
          it = x_vector_type(apr%ao,iv)

          call atom_wij(ad,ia,wij)
          call type_projector_data(apr%ao,it,pawd)

          nr = x_grid_size(pawd)
          if (size(n) /= nr) then
            deallocate( n, dn_dr, dn_dt, dn_dp, del_n ) ; allocate( n(nr), dn_dr(nr), dn_dt(nr), dn_dp(nr), del_n(nr) )
            deallocate( dfxc_dn, dfxc_dg ) ; allocate( dfxc_dn(nr), dfxc_dg(nr) )
            deallocate( r2dfxc_dn, r2dfxc_dr ) ; allocate( r2dfxc_dn(nr), r2dfxc_dr(nr) )
            deallocate( dfxc_dt, dfxc_dp ) ; allocate( dfxc_dt(nr), dfxc_dp(nr) )
            deallocate( nt, dnt_dr, dnt_dt, dnt_dp, del_nt ) ; allocate( nt(nr), dnt_dr(nr), dnt_dt(nr), dnt_dp(nr), del_nt(nr) )
            deallocate( dfxct_dn, dfxct_dg ) ; allocate( dfxct_dn(nr), dfxct_dg(nr) )
            deallocate( r2dfxct_dn, r2dfxct_dr ) ; allocate( r2dfxct_dn(nr), r2dfxct_dr(nr) )
            deallocate( dfxct_dt, dfxct_dp ) ; allocate( dfxct_dt(nr), dfxct_dp(nr) )
          end if
          n = 0.0_double
          dn_dr = 0.0_double
          dn_dt = 0.0_double
          dn_dp = 0.0_double
          nt = 0.0_double
          dnt_dr = 0.0_double
          dnt_dt = 0.0_double
          dnt_dp = 0.0_double

          call type_r2(pawd,r2)
          call type_wt(pawd,wt)
          call type_core_density(pawd,nc)
          call type_core_grad(pawd,dnc_dr)
          call type_coretail_density(pawd,nct)
          call type_coretail_grad(pawd,dnct_dr)
          call type_phir2ij(pawd,phir2ij)
          call type_dphir2ij_dr(pawd,dphir2ij_dr)
          call type_tphir2ij(pawd,tphir2ij)
          call type_dtphir2ij_dr(pawd,dtphir2ij_dr)

          nb = x_n_basis(apr%ao,it)
          call vector_ylmij(apr%ao,iv,ylmij)
          call vector_dylmij_dt(apr%ao,iv,dylmij_dt)
          call vector_dylmij_dp(apr%ao,iv,dylmij_dp)
          vwt = x_vector_weight(apr%ao,iv)

          do bj = 1,nb
            lj = x_basis_l(apr%ao,it,bj)
            oj = x_basis_offset(apr%ao,it,bj)
            osj = lj**2 + lj + 1
            do bi = 1,nb
              li = x_basis_l(apr%ao,it,bi)
              oi = x_basis_offset(apr%ao,it,bi)
              osi = li**2 + li + 1
              cij = (0.0_double,0.0_double)
              dcij_dt = (0.0_double,0.0_double)
              dcij_dp = (0.0_double,0.0_double)
              do mj = -lj,+lj
              do mi = -li,+li
                cij = cij + wij(oi+mi,oj+mj)*ylmij(osi+mi,osj+mj)
                dcij_dt = dcij_dt + wij(oi+mi,oj+mj)*dylmij_dt(osi+mi,osj+mj)
                dcij_dp = dcij_dp + wij(oi+mi,oj+mj)*dylmij_dp(osi+mi,osj+mj)
              end do
              end do
              n = n + real(cij,double)*phir2ij(:,bi,bj)
              dn_dr = dn_dr + real(cij,double)*dphir2ij_dr(:,bi,bj)
              dn_dt = dn_dt + real(dcij_dt,double)*phir2ij(:,bi,bj)
              dn_dp = dn_dp + real(dcij_dp,double)*phir2ij(:,bi,bj)
              nt = nt + real(cij,double)*tphir2ij(:,bi,bj)
              dnt_dr = dnt_dr + real(cij,double)*dtphir2ij_dr(:,bi,bj)
              dnt_dt = dnt_dt + real(dcij_dt,double)*tphir2ij(:,bi,bj)
              dnt_dp = dnt_dp + real(dcij_dp,double)*tphir2ij(:,bi,bj)
            end do
          end do
          n = n + spin_factor*nc
          dn_dr = dn_dr + spin_factor*dnc_dr
          nt = nt + spin_factor*nct
          dnt_dr = dnt_dr + spin_factor*dnct_dr

          del_n(1) = 1.0_double
          del_nt(1) = 1.0_double
          do ir = 2,nr
            del_n(ir) = sqrt(dn_dr(ir)**2 + dn_dt(ir)**2/r2(ir) + dn_dp(ir)**2/r2(ir))
            del_nt(ir) = sqrt(dnt_dr(ir)**2 + dnt_dt(ir)**2/r2(ir) + dnt_dp(ir)**2/r2(ir))
          end do

          call xc_derivatives(x_axc(apr%ao),n,del_n,lap_n,dfxc_dn,dfxc_dg,dfxc_dl)
          call xc_derivatives(x_axc(apr%ao),nt,del_nt,lap_nt,dfxct_dn,dfxct_dg,dfxct_dl)

          r2dfxc_dn = (r2*dfxc_dn*wt)*vwt
          r2dfxc_dn(1) = 0.0_double ; r2dfxc_dn(nr) = 0.0_double
          r2dfxc_dr = (r2*dfxc_dg*dn_dr*wt)*vwt
          r2dfxc_dr(1) = 0.0_double ; r2dfxc_dr(nr) = 0.0_double
          dfxc_dt = (dfxc_dg*dn_dt*wt)*vwt
          dfxc_dp = (dfxc_dg*dn_dp*wt)*vwt

          r2dfxct_dn = (r2*dfxct_dn*wt)*vwt
          r2dfxct_dn(nr) = 0.0_double
          r2dfxct_dr = (r2*dfxct_dg*dnt_dr*wt)*vwt
          r2dfxct_dr(nr) = 0.0_double
          dfxct_dt = (dfxct_dg*dn_dt*wt)*vwt
          dfxct_dp = (dfxct_dg*dn_dp*wt)*vwt

          do bj = 1,nb
            lj = x_basis_l(apr%ao,it,bj)
            oj = x_basis_offset(apr%ao,it,bj)
            osj = lj**2 + lj + 1
            do bi = 1,nb
              li = x_basis_l(apr%ao,it,bi)
              oi = x_basis_offset(apr%ao,it,bi)
              osi = li**2 + li + 1
              fn_r = sum(r2dfxc_dn*phir2ij(:,bi,bj) + r2dfxc_dr*dphir2ij_dr(:,bi,bj) &
                          - r2dfxct_dn*tphir2ij(:,bi,bj) - r2dfxct_dr*dtphir2ij_dr(:,bi,bj))
              fn_t = sum(dfxc_dt*phir2ij(:,bi,bj) - dfxct_dt*tphir2ij(:,bi,bj))
              fn_p = sum(dfxc_dp*phir2ij(:,bi,bj) - dfxct_dp*tphir2ij(:,bi,bj))
              do mj = -lj,+lj
              do mi = -li,+li
                dij(iu)%mat(oi+mi,oj+mj) = dij(iu)%mat(oi+mi,oj+mj) + fn_r*ylmij(osi+mi,osj+mj) &
                                                                    + fn_t*dylmij_dt(osi+mi,osj+mj) &
                                                                    + fn_p*dylmij_dp(osi+mi,osj+mj)
              end do
              end do
            end do
          end do

        end do

        if (associated( n )) then
          deallocate( n, dn_dr, dn_dt, dn_dp, del_n )
          deallocate( dfxc_dn, dfxc_dg )
          deallocate( r2dfxc_dn, r2dfxc_dr, dfxc_dt, dfxc_dp )
          deallocate( nt, dnt_dr, dnt_dt, dnt_dp, del_nt )
          deallocate( dfxct_dn, dfxct_dg )
          deallocate( r2dfxct_dn, r2dfxct_dr, dfxct_dt, dfxct_dp )
        end if
        nullify( r2 )
        nullify( wt )
        nullify( nc )
        nullify( dnc_dr )
        nullify( nct )
        nullify( dnct_dr )
        nullify( phir2ij )
        nullify( dphir2ij_dr )
        nullify( tphir2ij )
        nullify( dtphir2ij_dr )
        nullify( wij )
        nullify( ylmij )
        nullify( dylmij_dt )
        nullify( dylmij_dp )
        nullify( pawd )

        call glean(thy(ad))

        if (error("Exit atomic_potential_paw_mod::accum_xc_grad_i")) continue

      end subroutine 
  
      subroutine hat_terms_i(apr,ad,dij) ! PRB 55, 2005 (1997): Equations A23 and A31 (4th term).
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(gp_mat), dimension(:), pointer :: dij

        integer :: i, ia, it, l, m, qb
        integer :: bi, bj, mi, mj, omi, omj
        integer :: iu, nu, nu_proc, first_u, last_u
        complex(double) :: c1
        complex(double), dimension(:,:), pointer :: a_qlm
        type(paw_data_obj), pointer :: pawd

        call my(ad)

        nu = x_n_unique_atoms(apr%ao)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u

          ia = x_unique_atom(apr%ao,iu)
          it = x_unique_type(apr%ao,iu)
          call atom_qlm(ad,ia,a_qlm)
          call type_projector_data(apr%ao,it,pawd)
          qb = 2*x_type_l_max(apr%ao,it) + 1

          do i = 1,x_qvlm_size(pawd)
            call orbital_decode(pawd,i,bi,bj,mi,mj,l)
            m = mi - mj
            omi = x_basis_offset(apr%ao,it,bi) + mi
            omj = x_basis_offset(apr%ao,it,bj) + mj
            c1 = x_avlm(pawd,i)*a_qlm(l+1,qb+m)
            dij(iu)%mat(omi,omj) = dij(iu)%mat(omi,omj) - c1
          end do

        end do

        nullify( a_qlm )
        nullify( pawd )

        call glean(thy(ad))
 
        if (error("Exit atomic_potential_paw_mod::hat_terms_i")) continue

      end subroutine

      subroutine multipole_terms_i(apr,ad,hap,dij) ! PRB 55, 2005 (1997): Equation A27.
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(grid_obj) :: hap
        type(gp_mat), dimension(:), pointer :: dij

        integer :: i, ia, it, l, m, qb
        integer :: bi, bj, mi, mj, omi, omj
        integer :: iu, nu, nu_proc, first_u, last_u
        type(paw_data_obj), pointer :: pawd
        type(gp_mat), dimension(:), pointer :: de_dqlm
 
        call my(ad)
        call my(hap)
 
        call form_de_dqlm_i(apr,ad,hap,de_dqlm) ; if (error()) goto 100

        nu = x_n_unique_atoms(apr%ao)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u
          ia = x_unique_atom(apr%ao,iu)
          it = x_unique_type(apr%ao,iu)
          call type_projector_data(apr%ao,it,pawd)
          qb = 2*x_type_l_max(apr%ao,it) + 1
          do i = 1,x_qvlm_size(pawd)
            call orbital_decode(pawd,i,bi,bj,mi,mj,l)
            m = mj - mi
            omi = x_basis_offset(apr%ao,it,bi) + mi
            omj = x_basis_offset(apr%ao,it,bj) + mj
            dij(iu)%mat(omi,omj) = dij(iu)%mat(omi,omj) + x_aqlm(pawd,i)*de_dqlm(iu)%mat(l+1,qb+m)
          end do
        end do

100     if (associated( de_dqlm )) then
          do iu = 1,size(de_dqlm)
            if (associated( de_dqlm(iu)%mat )) deallocate( de_dqlm(iu)%mat )
          end do
          deallocate( de_dqlm )
        end if
        nullify( pawd )

        call glean(thy(ad))
        call glean(thy(hap))

        if (error("Exit atomic_potential_paw_mod::multipole_terms_i")) continue

      end subroutine 

      subroutine form_de_dqlm_i(apr,ad,hap,de_dqlm) ! PRB 55, 2005 (1997): Equation A28.
        type(atomic_potential_paw_rep) :: apr
        type(atomic_density_paw_obj) :: ad
        type(grid_obj) :: hap
        type(gp_mat), dimension(:), pointer :: de_dqlm

        integer :: i, ia, it, iu, l, m, ml, nu, qb
        integer :: bi, bj, mi, mj, omi, omj
        real(double) :: r1
        complex(double), dimension(:,:), pointer :: a_wijs, a_qlm
        type(paw_data_obj), pointer :: pawd

        call my(ad)
        call my(hap)
 
        nu = x_n_unique_atoms(apr%ao)

        allocate( de_dqlm(nu) )

        do iu = 1,nu

          ia = x_unique_atom(apr%ao,iu)
          it = x_unique_type(apr%ao,iu)
          call atom_wijs(ad,ia,a_wijs)
          call atom_qlm(ad,ia,a_qlm)
          call type_projector_data(apr%ao,it,pawd)
          ml = x_type_l_max(apr%ao,it)
          qb = 2*ml + 1

          allocate( de_dqlm(iu)%mat(size(a_qlm,1),size(a_qlm,2)) )
          de_dqlm(iu)%mat = (0.0_double,0.0_double)

          do i = 1,x_qvlm_size(pawd)
            call orbital_decode(pawd,i,bi,bj,mi,mj,l)
            m = mi - mj
            omi = x_basis_offset(apr%ao,it,bi) + mi
            omj = x_basis_offset(apr%ao,it,bj) + mj
            de_dqlm(iu)%mat(l+1,qb+m) = de_dqlm(iu)%mat(l+1,qb+m) - a_wijs(omi,omj)*x_avlm(pawd,i)
          end do

          do l = 0,2*ml
            r1 = 2.0_double*x_hat_self_energy(pawd,l)
            do m = -l,+l
              de_dqlm(iu)%mat(l+1,qb+m) = de_dqlm(iu)%mat(l+1,qb+m) - r1*conjg(a_qlm(l+1,qb+m))
            end do
          end do

          de_dqlm(iu)%mat(1,qb) = de_dqlm(iu)%mat(1,qb) - x_coretail_hatenergy(pawd)

        end do
 
        call hap_terms_i(apr%ao,hap,de_dqlm)

        nullify( a_wijs )
        nullify( a_qlm )
        nullify( pawd )

        call glean(thy(ad))
        call glean(thy(hap))
  
        if (error("Exit atomic_potential_paw_mod::form_de_dqlm_i")) continue

      end subroutine 

      subroutine hap_terms_i(ao,hap,de_dqlm) ! PRB 55, 2005 (1997): Equation A28 (1st term).
        type(atomic_operators_paw_obj) :: ao
        type(grid_obj) :: hap
        type(gp_mat), dimension(:), pointer :: de_dqlm

        integer :: i1, i2, i3, ia, it, iu, l, m, n, n0, nu, qb
        real(double) :: volume
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        real(double), dimension(:,:,:,:), pointer :: hdff
        complex(double) :: c0, c1, igr
        complex(double), dimension(13) :: ylm
        complex(double), dimension(:), allocatable :: c_sum1, c_sum2
        complex(double), dimension(:,:,:), pointer :: c_hap
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(layout_obj) :: lay

        call my(ao)
        call my(hap)

        call my(x_lattice(x_crystal(ao)),lat)
        call my(x_atoms(x_crystal(ao)),ats)
        call my(x_layout(hap),lay)

        nullify( gx, gy, gz, c_hap )

        if (error(x_ghost(lay) /= x_layout_ghost(ao),"ERROR: inconsistent layouts")) goto 100

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call take(c_hap,hap,CDF_KIND)

        nu = x_n_unique_atoms(ao)

        n = 0
        do iu = 1,nu
          it = x_unique_type(ao,iu)
          call hat_density_ff(ao,it,hdff)
          do l = 0,(size(hdff,4) - 1)
            n = n + 2*l + 1
          end do
        end do
        allocate( c_sum1(n), c_sum2(n) )
        c_sum1 = (0.0_double,0.0_double)

        n0 = 0
        do iu = 1,nu
          ia = x_unique_atom(ao,iu)
          it = x_unique_type(ao,iu)
          call hat_density_ff(ao,it,hdff)
          pos = lat2r(lat,x_position(ats,ia))
          do i3 = 1,size(gx,3)
          do i2 = 1,size(gx,2)
          do i1 = 1,size(gx,1)
            if (hdff(i1,i2,i3,1) == 0.0_double) cycle
            igr = (0.0_double,1.0_double)*(pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3))
            c0 = conjg(c_hap(i1,i2,i3))*exp(-igr)
            n = n0
            do l = 0,(size(hdff,4) - 1)
              c1 = c0*hdff(i1,i2,i3,l+1)*(0.0_double,-1.0_double)**l
              ylm = spharm(gx(i1,i2,i3),gy(i1,i2,i3),gz(i1,i2,i3),l,.true.)
              do m = -l,+l
                n = n + 1
                c_sum1(n) = c_sum1(n) + c1*ylm(l+m+1)
              end do
            end do
          end do
          end do
          end do
          n0 = n
        end do

        call allreduce(SGROUP,MPI_SUM,c_sum1,c_sum2)

        call put(c_hap,hap,CDF_KIND)

        volume = x_cell_volume(lat)

        n = 0
        do iu = 1,nu
          it = x_unique_type(ao,iu)
          call hat_density_ff(ao,it,hdff)
          qb = 2*x_type_l_max(ao,it) + 1  
          do l = 0,(size(hdff,4) - 1)
            do m = -l,+l
              n = n + 1
              de_dqlm(iu)%mat(l+1,qb+m) = de_dqlm(iu)%mat(l+1,qb+m) + c_sum2(n)*volume
            end do 
          end do
        end do

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (allocated( c_sum1 )) deallocate( c_sum1 )
        if (allocated( c_sum2 )) deallocate( c_sum2 )
        if (associated( c_hap )) deallocate( c_hap )

        nullify( hdff )

        call glean(thy(lat))
        call glean(thy(ats))
        call glean(thy(lay))
        
        call glean(thy(ao))
        call glean(thy(hap))

        if (error("Exit atomic_potential_paw_mod::hap_terms_i")) continue
        
      end subroutine

      subroutine distribute_dij_i(apr,dij)
        type(atomic_potential_paw_rep) :: apr
        type(gp_mat), dimension(:), pointer :: dij

        logical, dimension(:), allocatable :: filled
        integer :: i1d, i1, i2, ia, iau, it, iu, n1d, na, nb, np, nu
        integer :: bi, bj, li, lj, mi, mj, oi, oj
        complex(double), dimension(:), allocatable :: c1d, c1d_all
        complex(double), dimension(:,:,:), allocatable :: d_block
        type(space_group_obj) :: sg

        call my(x_space_group(apr%ao),sg)

        ! Sum the dij contributions.
        nu = x_n_unique_atoms(apr%ao)
        n1d = 0
        do iu = 1,nu
          n1d = n1d + size(dij(iu)%mat)
        end do
        allocate( c1d(n1d), c1d_all(n1d) )
        i1d = 0
        do iu = 1,nu
          do i2 = 1,size(dij(iu)%mat,2)
          do i1 = 1,size(dij(iu)%mat,1)
            i1d = i1d + 1
            c1d(i1d) = dij(iu)%mat(i1,i2)
          end do
          end do
        end do
        call allreduce(SGROUP,MPI_SUM,c1d,c1d_all)
        i1d = 0
        do iu = 1,nu
          do i2 = 1,size(dij(iu)%mat,2)
          do i1 = 1,size(dij(iu)%mat,1)
            i1d = i1d + 1
            dij(iu)%mat(i1,i2) = c1d_all(i1d)
          end do
          end do
        end do
        deallocate( c1d, c1d_all )

        na = x_n_atoms(apr%ao)

        ! Allocate space for apr%dij.
        allocate( apr%dij(na) )
        do ia = 1,na
          np = x_n_atom_projectors(apr%ao,ia)
          allocate( apr%dij(ia)%mat(np,np) )
        end do

        ! Distribute the dij values to apr%dij.
        if (nu == na) then
          do ia = 1,na
            apr%dij(ia)%mat = dij(ia)%mat
          end do
        else
          allocate( filled(na) )
          do iu = 1,nu
            iau = x_unique_atom(apr%ao,iu)
            it = x_atom_type(apr%ao,iau)
            nb = x_n_basis(apr%ao,it)
            do bi = 1,nb
              li = x_basis_l(apr%ao,it,bi)
              oi = x_basis_offset(apr%ao,it,bi)
              do bj = 1,nb
                lj = x_basis_l(apr%ao,it,bj)
                oj = x_basis_offset(apr%ao,it,bj)
                allocate( d_block(-li:+li,-lj:+lj,na))
                do mi = -li,+li
                do mj = -lj,+lj
                  d_block(mi,mj,iau) = dij(iu)%mat(oi+mi,oj+mj)
                end do
                end do
                call distribute_spherical_tensor(sg,iau,d_block,filled)
                do ia = 1,na
                  if (.not.filled(ia)) cycle
                  do mi = -li,+li
                  do mj = -lj,+lj
                    apr%dij(ia)%mat(oi+mi,oj+mj) = d_block(mi,mj,ia)
                  end do
                  end do
                end do
                deallocate( d_block )
              end do
            end do
          end do
        end if

        if (allocated( filled )) deallocate( filled )
        if (allocated( c1d )) deallocate( c1d )
        if (allocated( c1d_all )) deallocate( c1d_all )
        if (allocated( d_block )) deallocate( d_block )

        call glean(thy(sg))

        if (error("Exit atomic_density_paw_mod::distribute_dij_i")) continue

      end subroutine

      subroutine own_i(ap)
        type(atomic_potential_paw_obj) :: ap
        type(atomic_potential_paw_obj) :: apt

        integer :: i
        if (ap%ref < ap%o%ref) then
          allocate( apt%o )
          apt%o%ref = 0
          apt%o%g = ap%o%g
          call my(ap%o%ao,apt%o%ao)
          allocate( apt%o%dij(size(ap%o%dij)) )
          do i = 1,size(apt%o%dij)
            allocate( apt%o%dij(i)%mat(size(ap%o%dij(i)%mat,1),size(ap%o%dij(i)%mat,2)) )
            apt%o%dij(i)%mat = ap%o%dij(i)%mat
          end do
          ap%o%ref = ap%o%ref - ap%ref
          ap%o => apt%o
          ap%o%ref = ap%o%ref + ap%ref
        end if
      end subroutine

      end module

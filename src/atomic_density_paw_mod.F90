! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module atomic_density_paw_mod
!doc$ module atomic_density_paw_mod

      use kind_mod
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
      use axc_mod
      use math_mod
      use symmetry_mod
      use paw_data_mod
      use atomic_operators_paw_mod

!     One datatype is available here: type(atomic_density_paw_obj).

!cod$
      implicit none
      private

      type :: gp_mat
        complex(double), dimension(:,:), pointer :: mat
      end type

      type :: atomic_density_paw_rep
        integer :: ref
        type(ghost) :: g
        real(double) :: charge_state                    ! charge_state
        real(double) :: energy                          ! atomic energy
        type(atomic_operators_paw_obj) :: ao            ! atomic operators object
        type(gp_mat), dimension(:), pointer :: wij      ! projected occupation coefficients (for an sgroup)
        type(gp_mat), dimension(:), pointer :: wijs     ! projected occupation coefficients (sum of sgroup contributions)
        type(gp_mat), dimension(:), pointer :: qlm      ! compensation charge multipole moments (sum of sgroup contributions)
      end type

      type, public :: atomic_density_paw_obj
        private
        integer :: ref
        type(atomic_density_paw_rep), pointer :: o
      end type

!doc$
      public :: atomic_density_paw
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_atomic_operators
      public :: atom_wij
      public :: atom_wijs
      public :: atom_qlm
      public :: distance
      public :: symmetrize
      public :: merge_atomic_density
      public :: add_atomic_density
      public :: atomic_hartree_density
      public :: guess_density
      public :: atomic_energy
      public :: atomic_forces
      public :: extract_density
      public :: insert_density
      public :: get_normalization
      public :: write_restart

!cod$
      interface atomic_density_paw
        module procedure constructor_ad
      end interface
      interface update
        module procedure update_ad
      end interface
      interface my
        module procedure my_ad, my_new_ad
      end interface
      interface thy
        module procedure thy_ad
      end interface
      interface glean
        module procedure glean_ad
      end interface
      interface bequeath
        module procedure bequeath_ad
      end interface
      interface assignment(=)
        module procedure assign_ad
      end interface 
      interface x_ghost
        module procedure ad_ghost
      end interface
      interface x_ref
        module procedure ad_ref
      end interface
      interface x_atomic_operators
        module procedure ad_atomic_operators
      end interface
      interface distance
        module procedure distance_ad
      end interface
      interface symmetrize
        module procedure symmetrize_ad
      end interface
      interface merge_atomic_density
        module procedure merge_atomic_density_ad
      end interface
      interface add_atomic_density
        module procedure add_atomic_density_ad
      end interface
      interface atomic_hartree_density
        module procedure atomic_hartree_density_ad
      end interface
      interface guess_density
        module procedure guess_density_ad
      end interface
      interface atomic_energy
        module procedure atomic_energy_ad
      end interface
      interface atomic_forces
        module procedure atomic_forces_ad
      end interface
      interface extract_density
        module procedure extract_density_ad
      end interface
      interface insert_density
        module procedure insert_density_ad
      end interface
      interface get_normalization
        module procedure get_normalization_ad
      end interface
      interface write_restart
        module procedure write_restart_ad
      end interface

      contains

! public routines

      function constructor_ad(ao,restf,empty) result(ad)
!doc$ function atomic_density_paw(ao,restf,empty) result(ad)
        type(atomic_operators_paw_obj) :: ao
        type(tagio_obj), optional :: restf
        logical, optional :: empty
        type(atomic_density_paw_obj) :: ad
!       requires: restf pointer be positioned inside the FIELDS block.
!       effects: Constructs a new ad. If empty is present and .true., an empty ad is created.
!       errors: Passes errors.

!cod$
        logical :: found, np_ok, zero
        character(1) :: tios
        character(line_len) :: tag
        integer :: i, i1, i2, ia, msg, n, na, nm, nsg
        integer(long) :: dsize, iosl, ndata, s4
        integer, dimension(2) :: csr
        integer, dimension(:), allocatable :: nps
        integer(long), dimension(:), allocatable :: v4
        complex(double), dimension(:), allocatable :: c1

        call my(ao)
        if (present(restf)) call my(restf)

        ad%ref = 0
        allocate( ad%o )
        ad%o%ref = 0
        ad%o%g = x_ghost()

        call my(ao,ad%o%ao)

        if (present(restf)) then

          nsg = mpi_nsgroups()
          msg = mpi_mysgroup()

          ! Open the ATOMIC_DENSITY block
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_DENSITY")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_DENSITY block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)

          ! Find the PAW tag
          if (i_access(restf)) tios = findfirsttag(restf,"PAW")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: PAW tag was not found")) goto 100

          ! Read the charge state
          if (i_access(restf)) tios = findfirsttag(restf,"CHARGE_STATE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: CHARGE_STATE tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
            call readf(ad%o%charge_state,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,ad%o%charge_state)

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
          if (error(na /= x_n_atoms(ad%o%ao),"ERROR: different numbers of atoms")) goto 100

          ! Read the number of projectors.
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
            if (nps(ia) /= x_n_atom_projectors(ad%o%ao,ia)) then
              np_ok = .false.
              exit
            end if
          end do
          if (error(.not.np_ok,"ERROR: different numbers of projectors")) goto 100

          ! Set the file pointer at the beginning of the wij matrix elements
          if (i_access(restf)) tios = findfirsttag(restf,"WIJ_MATRIX_ELEMENTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: WIJ_MATRIX_ELEMENTS tag was not found")) goto 100

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
            allocate( ad%o%wij(na) )
            i = 0
            do ia = 1,size(ad%o%wij)
              nm = nps(ia)
              allocate( ad%o%wij(ia)%mat(nm,nm) )
              do i2 = 1,nm
              do i1 = 1,nm
                i = i + 1
                ad%o%wij(ia)%mat(i1,i2) = c1(i)
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
            allocate( ad%o%wij(na) )
            i = 0
            do ia = 1,size(ad%o%wij)
              nm = nps(ia)
              allocate( ad%o%wij(ia)%mat(nm,nm) )
              do i2 = 1,nm
              do i1 = 1,nm
                i = i + 1
                ad%o%wij(ia)%mat(i1,i2) = c1(i)
              end do
              end do
            end do
          end select

          zero = .false.

          ! Close the ATOMIC_DENSITY block
100       if (i_access(restf)) call closeblock(restf)
          if (error()) goto 200

        else

          ! Read the charge state
          call arglc("charge_state_mode",tag,found)
          if (.not.found) tag = "real_number"
          select case (trim(tag))
          case ("real_number")
            call arg("charge_state",ad%o%charge_state,found)
            if (.not.found) ad%o%charge_state = 0.0_double
          case ("integer_ratio")
            call arg("charge_state_ratio",csr,found)
            if (error(.not.found,"ERROR: charge_state_ratio was not found")) goto 200
            if (error(csr(2) == 0,"ERROR: denominator = 0")) goto 200
            ad%o%charge_state = real(csr(1),double)/real(csr(2),double)
          case default
            if (error(.true.,"ERROR: charge_state_mode was not recognized")) goto 200
          end select


          ! Form the wij matrix elements
          zero = .false.
          if (present(empty)) zero = empty
          call form_wij_i(ad%o,zero) ; if (error()) goto 200

        end if

        nullify( ad%o%wijs )
        nullify( ad%o%qlm )
        if (zero) then
          ad%o%energy = 0.0_double
        else
          call form_wijs_i(ad%o)    ; if (error()) goto 200
          call form_qlm_i(ad%o)     ; if (error()) goto 200
          call form_energy_i(ad%o)
        end if

200     if (allocated( nps )) deallocate( nps )
        if (allocated( v4 )) deallocate( v4 )
        if (allocated( c1 )) deallocate( c1 )

        call glean(thy(ao))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_density_paw_mod::constructor_ad")) continue

      end function 

      subroutine update_ad(ad,ao)
!doc$ subroutine update(ad,ao)
        type(atomic_density_paw_obj) :: ad
        type(atomic_operators_paw_obj) :: ao
!       requires: The number and types of atoms have not changed.
!       modifies: ad
!       effects: Updates ad with respect to ao.

!cod$
        call my(ad)
        call my(ao)

        call own_i(ad)
        ad%o%g = x_ghost()
        ad%o%ao = ao

        call glean(thy(ao))
        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::update_ad")) continue

      end subroutine

      subroutine my_ad(ad)
!doc$ subroutine my(ad)
        type(atomic_density_paw_obj) :: ad 
 
!cod$
        ad%ref = ad%ref + 1
        ad%o%ref = ad%o%ref + 1
      end subroutine

      subroutine my_new_ad(adi,ad)
!doc$ subroutine my(adi,ad)
        type(atomic_density_paw_obj) :: adi, ad

!cod$
        ad%ref = 1
        ad%o => adi%o
        ad%o%ref = ad%o%ref + 1
      end subroutine

      function thy_ad(ad) result(ado)
!doc$ function thy(ad) result(ado)
        type(atomic_density_paw_obj) :: ad, ado

!cod$
        ad%ref = ad%ref - 1
        ad%o%ref = ad%o%ref - 1
        ado%ref = ad%ref
        ado%o => ad%o
      end function

      subroutine glean_ad(ad)
!doc$ subroutine glean(ad)
        type(atomic_density_paw_obj) :: ad

!cod$
        integer :: i
        if (ad%o%ref < 1) then
          select case (mpi_nsgroups())
          case (1)
            nullify( ad%o%wijs )
          case (2)
            if (associated( ad%o%wijs )) then
              do i = 1,size(ad%o%wijs)
                if (associated( ad%o%wijs(i)%mat )) deallocate( ad%o%wijs(i)%mat )
              end do
              deallocate( ad%o%wijs )
            end if
          end select
          if (associated( ad%o%wij )) then
            do i = 1,size(ad%o%wij)
              if (associated( ad%o%wij(i)%mat )) deallocate( ad%o%wij(i)%mat )
            end do
            deallocate( ad%o%wij )
          end if
          if (associated( ad%o%qlm )) then
            do i = 1,size(ad%o%qlm)
              if (associated( ad%o%qlm(i)%mat )) deallocate( ad%o%qlm(i)%mat )
            end do
            deallocate( ad%o%qlm )
          end if
          call glean(thy(ad%o%ao))
          deallocate( ad%o )
        end if
      end subroutine

      subroutine bequeath_ad(ad)
!doc$ subroutine bequeath(ad)
        type(atomic_density_paw_obj) :: ad

!cod$
        continue
      end subroutine

      subroutine assign_ad(ad,ad2)
!doc$ subroutine assignment(=)(ad,ad2)
        type(atomic_density_paw_obj), intent(inout) :: ad
        type(atomic_density_paw_obj), intent(in) :: ad2
 
!cod$
        type(atomic_density_paw_obj) :: adt
        call my(ad2)
        adt%o => ad%o
        ad%o%ref = ad%o%ref - ad%ref
        ad%o => ad2%o
        ad%o%ref = ad%o%ref + ad%ref
        call glean(adt)
        call glean(thy(ad2))
      end subroutine

      function ad_ref(ad) result(r)
!doc$ function x_ref(ad) result(r)
        type(atomic_density_paw_obj) :: ad
        integer, dimension(2) :: r
!       effects: Returns the reference counts of ad.

!cod$
        r(1) = ad%ref
        r(2) = ad%o%ref
        call glean(ad)
      end function 

      function ad_ghost(ad) result(g)
!doc$ function x_ghost(ad) result(g)
        type(atomic_density_paw_obj) :: ad
        type(ghost) :: g
!       effects: Returns the ghost of ad.

!cod$ 
        call my(ad)
        g = ad%o%g
        call glean(thy(ad))
      end function 

      function ad_atomic_operators(ad) result(ao)
!doc$ function x_atomic_operators(ad) result(ao)
        type(atomic_density_paw_obj) :: ad
        type(atomic_operators_paw_obj) :: ao
!       effects: Returns ad%o%ao.
  
!cod$
        call my(ad)
        call my(ad%o%ao,ao)
        call glean(thy(ad))
        call bequeath(thy(ao))
      end function 

      subroutine atom_wij(ad,ia,a_wij)
!doc$ subroutine atom_wij(ad,ia,a_wij)
        type(atomic_density_paw_obj) :: ad
        integer, intent(in) :: ia
        complex(double), dimension(:,:), pointer :: a_wij
!       requires: ia be in range.
!       modifies: a_wij
!       effects: Points a_wij to ad%o%wij(ia)%mat

!cod$
        a_wij => ad%o%wij(ia)%mat
        call glean(ad)
      end subroutine

      subroutine atom_wijs(ad,ia,a_wijs)
!doc$ subroutine atom_wijs(ad,ia,a_wijs)
        type(atomic_density_paw_obj) :: ad
        integer, intent(in) :: ia
        complex(double), dimension(:,:), pointer :: a_wijs
!       requires: ia be in range.
!       modifies: a_wijs
!       effects: Points a_wijs to ad%o%wijs(ia)%mat

!cod$
        a_wijs => ad%o%wijs(ia)%mat
        call glean(ad)
      end subroutine

      subroutine atom_qlm(ad,ia,a_qlm)
!doc$ subroutine atom_qlm(ad,ia,a_qlm)
        type(atomic_density_paw_obj) :: ad
        integer, intent(in) :: ia
        complex(double), dimension(:,:), pointer :: a_qlm
!       requires: ia be in range.
!       modifies: a_qlm
!       effects: Points a_qlm to ad%o%qlm(ia)%mat

!cod$
        a_qlm => ad%o%qlm(ia)%mat
        call glean(ad)
      end subroutine

      function distance_ad(ad1,ad2) result(d)
!doc$ function distance(ad1,ad2) result(d)
        type(atomic_density_paw_obj) :: ad1, ad2
        real(double) :: d
!       requires: ad1 and ad2 have the same sizes.
!       effects: d = sart[sum{|ad1%o%wij1 - ad1%o%wij2|^2}].

!cod$
        integer :: i1, i2, ia, first_a, last_a, na, na_proc
        real(double) :: d2p, d2
        complex(double) :: c0

        na = size(ad1%o%wij)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,na,first_a,last_a,na_proc)
        d2p = 0.0_double
        do ia = first_a,last_a
          do i2 = 1,size(ad1%o%wij(ia)%mat,2)
          do i1 = 1,size(ad1%o%wij(ia)%mat,1)
            c0 = ad1%o%wij(ia)%mat(i1,i2) - ad2%o%wij(ia)%mat(i1,i2)
            d2p = d2p + real(c0*conjg(c0),double)
          end do
          end do
        end do
        call allreduce(SGROUP,MPI_SUM,d2p,d2)
        d = sqrt(d2)

        call glean(ad1)
        call glean(ad2)

      end function 

      subroutine symmetrize_ad(ad,sg)
!doc$ subroutine symmetrize(ad,sg)
        type(atomic_density_paw_obj) :: ad
        type(space_group_obj) :: sg
!       modifies: ad
!       effects: Symmetrizes ad with respect to sg.
!       errors: Passes errors.

!cod$
        logical, dimension(:), allocatable :: mask
        integer :: ia, it, na, nb, nt
        integer :: bi, bj, li, lj, mi, mj, oi, oj
        complex(double), dimension(:,:,:), allocatable :: w_block
    
        call my(ad)
        call my(sg)

        call own_i(ad)
        ad%o%g = x_ghost()

        na = x_n_atoms(ad%o%ao)
        nt = x_n_types(ad%o%ao)
        allocate( mask(na) )

        do it = 1,nt
          nb = x_n_basis(ad%o%ao,it)
          do bi = 1,nb
            li = x_basis_l(ad%o%ao,it,bi)
            oi = x_basis_offset(ad%o%ao,it,bi)
            do bj = 1,nb
              lj = x_basis_l(ad%o%ao,it,bj)
              oj = x_basis_offset(ad%o%ao,it,bj)
              allocate( w_block(-li:+li,-lj:+lj,na))
              mask = .false.
              do ia = 1,na
                if (x_atom_type(ad%o%ao,ia) == it) then
                  mask(ia) = .true.
                  do mj = -lj,+lj
                  do mi = -li,+li
                    w_block(mi,mj,ia) = ad%o%wij(ia)%mat(oi+mi,oj+mj)
                  end do
                  end do
                else
                  w_block(:,:,ia) = (0.0_double,0.0_double)
                end if
              end do
              call symmetrize_spherical_tensor(sg,w_block,mask) ; if (error()) goto 100
              do ia = 1,na
                if (mask(ia)) then
                  do mj = -lj,+lj
                  do mi = -li,+li
                    ad%o%wij(ia)%mat(oi+mi,oj+mj) = w_block(mi,mj,ia)
                  end do
                  end do
                end if
              end do
              deallocate( w_block )
            end do
          end do
        end do
  
        if (allocated( w_block )) deallocate( w_block )
        if (allocated( mask )) deallocate( mask )

        call form_wijs_i(ad%o)    ; if (error()) goto 100
        call form_qlm_i(ad%o)     ; if (error()) goto 100
        call form_energy_i(ad%o)

100     call glean(thy(ad))
        call glean(thy(sg))
  
        if (error("Exit atomic_density_paw_mod::symmetrize_ad")) continue        

      end subroutine 

      subroutine merge_atomic_density_ad(ad)
!doc$ subroutine merge_atomic_density(ad)
        type(atomic_density_paw_obj) :: ad
!       modifies: ad
!       effects: Merges ad from different kgroups.
!       errors: Passes errors.

!cod$
        integer :: i, i1, i2, ia, first_a, last_a, n, na, na_proc
        complex(double), dimension(:), allocatable :: c1, c2

        call my(ad)

        call own_i(ad)
        ad%o%g = x_ghost()

        na = size(ad%o%wij)

        ! Copy the non-zero matrix elements into c1.
        n = 0
        do ia = 1,size(ad%o%wij)
          n = n + size(ad%o%wij(ia)%mat)
        end do
        allocate( c1(n), c2(n) )
        call subdivide(mpi_myproc(KGROUP),mpi_nprocs(KGROUP),1,na,first_a,last_a,na_proc)
        i = 0
        do ia = 1,na
          if (ia == first_a) exit
          i = i + size(ad%o%wij(ia)%mat)
        end do
        c1 = (0.0_double,0.0_double)
        do ia = first_a,last_a
          do i2 = 1,size(ad%o%wij(ia)%mat,2)
          do i1 = 1,size(ad%o%wij(ia)%mat,1)
            i = i + 1
            c1(i) = ad%o%wij(ia)%mat(i1,i2)
          end do
          end do
        end do

        ! Combine the contributions to the matrix elements.
        call allreduce(SGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100

        ! Copy the matrix elements back into ad%o%wij.
        i = 0
        do ia = 1,size(ad%o%wij)
          do i2 = 1,size(ad%o%wij(ia)%mat,2)
          do i1 = 1,size(ad%o%wij(ia)%mat,1)
            i = i + 1
            ad%o%wij(ia)%mat(i1,i2) = c2(i)
          end do
          end do
        end do

100     if (allocated( c1 )) deallocate( c1 )
        if (allocated( c2 )) deallocate( c2 )

        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::merge_atomic_density_ad")) continue

      end subroutine

      subroutine add_atomic_density_ad(ad,pdots,weights)
!doc$ subroutine add_atomic_density(ad,pdots,weights)
        type(atomic_density_paw_obj) :: ad
        complex(double), dimension(:,:), intent(in) :: pdots
        real(double), dimension(:), intent(in) :: weights
!       modifies: ad
!       effects: Accumulates contributions from pdots into ad%o%wij.
!       errors: Passes errors.
!       notes: This routine is called from within the KGROUP scope.

!cod$
        character(1), parameter :: transa = 'n', transb = 't'
        integer :: bi, first_a, last_a, ia, ib, it, na, na_proc, nb, nm, np, nt
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: wcpdots
        complex(double), dimension(:,:), pointer :: t_tr2c, t_tc2r
        type(gp_mat), dimension(:), pointer :: twij

        call my(ad)

        call own_i(ad)
        ad%o%g = x_ghost()

        na = size(ad%o%wij)
        call subdivide(mpi_myproc(KGROUP),mpi_nprocs(KGROUP),1,na,first_a,last_a,na_proc)

        nt = x_n_types(ad%o%ao)
        allocate( twij(nt) )
        do it = 1,size(twij)
          nm = x_n_type_projectors(ad%o%ao,it)
          allocate( twij(it)%mat(nm,nm) )
        end do

        np = size(pdots,1)
        nb = size(pdots,2)
        allocate( wcpdots(np,nb) )
        do ib = 1,nb
          wcpdots(:,ib) = weights(ib)*conjg(pdots(:,ib))
        end do
        do ia = first_a,last_a
          nm = size(ad%o%wij(ia)%mat,1)
          bi = x_atom_base(ad%o%ao,ia)
          it = x_atom_type(ad%o%ao,ia)
          call type_tr2c(ad%o%ao,it,t_tr2c)
          call type_tc2r(ad%o%ao,it,t_tc2r)
          call zgemm(transa,transb,nm,nm,nb,alpha,wcpdots(bi,1),np,pdots(bi,1),np,beta,twij(it)%mat,nm)
          ad%o%wij(ia)%mat = ad%o%wij(ia)%mat + matmul(t_tc2r,matmul(twij(it)%mat,t_tr2c))
        end do

        if (allocated( wcpdots )) deallocate( wcpdots )
        if (associated( twij )) then
          do it = 1,size(twij)
            if (associated( twij(it)%mat )) deallocate( twij(it)%mat )
          end do
          deallocate( twij )
        end if
        nullify( t_tr2c )
        nullify( t_tc2r )

        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::add_atomic_density_ad")) continue

      end subroutine

      function atomic_hartree_density_ad(ad,lay) result(den)
!doc$ function atomic_hartree_density(ad,lay) result(den)
        type(atomic_density_paw_obj) :: ad
        type(layout_obj) :: lay
        type(grid_obj) :: den
!       effects: Returns the filtered hat density + filtered coretail density.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, it, l, m, qb
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        real(double), dimension(:,:,:,:), pointer :: hdff
        complex(double) :: c0, c1, igr
        complex(double), dimension(13) :: ylm
        complex(double), dimension(:,:,:), pointer :: c_den
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(ad)
        call my(lay)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)
        call my(x_atoms(x_crystal(ad%o%ao)),ats)
        call my(grid(lay,SGROUP),den)

        nullify( gx, gy, gz, c_den )

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)
        
        call alloc(c_den,lay,D_TYPE,SGROUP)
        c_den = (0.0_double,0.0_double)
        do it = 1,x_n_types(ad%o%ao)
          call hat_density_ff(ad%o%ao,it,hdff)
          qb = 2*x_type_l_max(ad%o%ao,it) + 1
          do ia = 1,x_n_atoms(ad%o%ao)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
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
                  c0 = c0 + c1*ylm(l+m+1)*ad%o%qlm(ia)%mat(l+1,qb+m)
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

        call add_coretail_density(ad%o%ao,den) ; if (error()) continue

        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )
        nullify( hdff )

        call glean(thy(lat))
        call glean(thy(ats))
        call bequeath(thy(den))

        call glean(thy(ad))
        call glean(thy(lay))

        if (error("Exit atomic_density_paw_mod::atomic_hartree_density_ad")) continue

      end function

      function guess_density_ad(ad,lay,ne) result(den)
!doc$ function guess_density(ad,lay,ne) result(den)
        type(atomic_density_paw_obj) :: ad
        type(layout_obj) :: lay
        real(double), optional :: ne
        type(grid_obj) :: den
!       effects: Returns a filtered approximation to the valence density.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, id, it
        integer :: bi, bj, l, li, lj, m, mi, mj, omi, omj
        real(double) :: ne_a, ne_g, spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:,:), pointer :: dmff
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        complex(double) :: c0, den_norm, igr, mil
        complex(double), dimension(13) :: ylm
        complex(double), dimension(:,:), pointer :: a_oij
        complex(double), dimension(:,:,:), pointer :: c_den
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(paw_data_obj), pointer :: pawd

        call my(ad)
        call my(lay)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)
        call my(x_atoms(x_crystal(ad%o%ao)),ats)
        call my(grid(lay,SGROUP),den)

        nullify( gx, gy, gz, c_den )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        if (present(ne)) then
          ne_g = ne
        else
          ne_a = 0.0_double
          do ia = 1,size(ad%o%wij)
            call atom_c_oij(ad%o%ao,ia,a_oij)
            ne_a = ne_a + real(sum(ad%o%wij(ia)%mat*a_oij),double)
          end do
          ne_g = spin_factor*(x_valence_electrons(ad%o%ao) - ad%o%charge_state) - ne_a
        end if
        den_norm = cmplx(ne_g,0,double)/x_cell_volume(lat)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call alloc(c_den,lay,D_TYPE,SGROUP)
        c_den = (0.0_double,0.0_double)
        do it = 1,x_n_types(ad%o%ao)
          call density_matrix_ff(ad%o%ao,it,dmff)
          call type_projector_data(ad%o%ao,it,pawd)
          do ia = 1,size(ad%o%wij)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
            pos = lat2r(lat,x_position(ats,ia))
            do id = 1,x_denmat_size(pawd)
              call denvhat_decode(pawd,id,bi,bj,l)
              mil = (0.0_double,-1.0_double)**l
              li = x_basis_l(ad%o%ao,it,bi)
              lj = x_basis_l(ad%o%ao,it,bj)
              do mi = -li,+li
              do mj = -lj,+lj
                m = mj - mi
                if (abs(m) <= l) then
                  omi = x_basis_offset(ad%o%ao,it,bi) + mi
                  omj = x_basis_offset(ad%o%ao,it,bj) + mj
                  c0 = ad%o%wij(ia)%mat(omi,omj)*gaunt_complex(l,m,li,mi,lj,mj)*mil
                  do i3 = 1,size(gx,3)
                  do i2 = 1,size(gx,2)
                  do i1 = 1,size(gx,1)
                    if (dmff(i1,i2,i3,1) == 0.0_double) cycle
                    igr = (0.0_double,1.0_double)*(pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3))
                    ylm = spharm(gx(i1,i2,i3),gy(i1,i2,i3),gz(i1,i2,i3),l,.true.)
                    c_den(i1,i2,i3) = c_den(i1,i2,i3) + c0*dmff(i1,i2,i3,id)*ylm(l+m+1)*exp(-igr)
                  end do
                  end do
                  end do
                end if
              end do
              end do
            end do
          end do
        end do
        call put(c_den,den,CDF_KIND)
        call set_normalization(den,den_norm)

        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )

        nullify( dmff )
        nullify( a_oij )
        nullify( pawd )

        call glean(thy(lat))
        call glean(thy(ats))
        call bequeath(thy(den))

        call glean(thy(ad))
        call glean(thy(lay))

        if (error("Exit atomic_density_paw_mod::guess_density_ad")) continue
        
      end function

      subroutine atomic_energy_ad(ad,e)
!doc$ subroutine atomic_energy(ad,e)
        type(atomic_density_paw_obj) :: ad
        real(double), intent(out) :: e
!       modifies: e
!       effects: Returns the atomic energy

!cod$
        call my(ad)
        e = ad%o%energy
        call glean(thy(ad))

      end subroutine

      subroutine atomic_forces_ad(ad,den,xcp,ahd,f)
!doc$ subroutine atomic_forces(ad,den,xcp,ahd,f)
        type(atomic_density_paw_obj) :: ad
        type(grid_obj) :: den
        type(grid_obj) :: xcp
        type(grid_obj) :: ahd
        real(double), dimension(:,:), intent(out) :: f
!       requires: f be dimensioned (3,x_n_atoms(ad%o%ao)).
!       modifies: f
!       effects: Returns unsymmetrized atomic contributions to the forces.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, it, l, m, qb
        real(double) :: r0, spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:,:), allocatable :: ft
        real(double), dimension(:,:,:), pointer :: gx, gy, gz, g2i
        real(double), dimension(:,:,:), pointer :: cdff, gpff
        real(double), dimension(:,:,:,:), pointer :: hdff
        complex(double) :: c0, c1, igr
        complex(double), dimension(13) :: ylm
        complex(double), dimension(:,:,:), pointer :: c_ahd, c_den, c_xcp, c_sum1, c_sum2
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(layout_obj) :: lay

        call my(ad)
        call my(den)
        call my(xcp)
        call my(ahd)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)
        call my(x_atoms(x_crystal(ad%o%ao)),ats)
        call my(x_layout(den),lay)

        nullify( gx, gy, gz, g2i, c_ahd, c_den, c_xcp, c_sum1, c_sum2 )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call fdelinv(g2i,lay,D_TYPE,SGROUP)

        call take(c_ahd,ahd,CDF_KIND)
        call take(c_den,den,CDF_KIND)
        call take(c_xcp,xcp,CDF_KIND)

        call alloc(c_sum1,lay,D_TYPE,SGROUP)
        call alloc(c_sum2,lay,D_TYPE,SGROUP)

        allocate( ft(size(f,1),size(f,2)) )
        ft = 0.0_double

        do it = 1,x_n_types(ad%o%ao)

          c_sum1 = eight_pi*g2i*conjg(c_den + spin_factor*c_ahd)

          call grid_potential_ff(ad%o%ao,it,gpff)
          call coretail_density_ff(ad%o%ao,it,cdff)
          c_sum2 = gpff*conjg(c_den) + cdff*(c_sum1 + spin_factor*conjg(c_xcp))

          call hat_density_ff(ad%o%ao,it,hdff)

          qb = 2*x_type_l_max(ad%o%ao,it) + 1

          do ia = 1,x_n_atoms(ad%o%ao)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
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
                  c0 = c0 + c1*ylm(l+m+1)*ad%o%qlm(ia)%mat(l+1,qb+m)
                end do
              end do
              igr = (0.0_double,1.0_double)*(pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3))
              r0 = aimag((c0*c_sum1(i1,i2,i3) + c_sum2(i1,i2,i3))*exp(-igr))
              ft(1,ia) = ft(1,ia) - gx(i1,i2,i3)*r0
              ft(2,ia) = ft(2,ia) - gy(i1,i2,i3)*r0
              ft(3,ia) = ft(3,ia) - gz(i1,i2,i3)*r0
            end do
            end do
            end do
          end do

        end do

        ft = x_cell_volume(lat)*ft
        call allreduce(CONFIG,MPI_SUM,ft,f)

100     call put(c_ahd,ahd,CDF_KIND)
        call put(c_den,den,CDF_KIND)
        call put(c_xcp,xcp,CDF_KIND)

        if (allocated( ft )) deallocate( ft )
        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( g2i )) deallocate( g2i )
        if (associated( c_ahd )) deallocate( c_ahd )
        if (associated( c_den )) deallocate( c_den )
        if (associated( c_xcp )) deallocate( c_xcp )
        if (associated( c_sum1 )) deallocate( c_sum1 )
        if (associated( c_sum2 )) deallocate( c_sum2 )

        nullify( cdff )
        nullify( gpff )
        nullify( hdff )

        call glean(thy(lat))
        call glean(thy(ats))
        call glean(thy(lay))
        
        call glean(thy(ad))
        call glean(thy(den))
        call glean(thy(xcp))
        call glean(thy(ahd))

        if (error("Exit atomic_density_paw_mod::atomic_forces_ad")) continue

      end subroutine

      subroutine extract_density_ad(ad,c1d)
!doc$ subroutine extract_density(ad,c1d)
        type(atomic_density_paw_obj) :: ad
        complex(double), dimension(:), pointer :: c1d
!       requires: c1d be nullified or associated.
!       modifies: c1d
!       effects: Reshapes ad%o%wij into c1d.

!cod$
        integer :: i, i1, i2, ia, first_a, last_a, n1d, na, na_proc

        call my(ad)

        na = x_n_atoms(ad%o%ao)
        call subdivide(mpi_myproc(SGROUP),mpi_nprocs(SGROUP),1,na,first_a,last_a,na_proc)

        n1d = 0
        do ia = first_a,last_a
          n1d = n1d + size(ad%o%wij(ia)%mat)
        end do
        if (associated( c1d )) deallocate( c1d ) ; allocate( c1d(n1d) )

        i = 0
        do ia = first_a,last_a
          do i2 = 1,size(ad%o%wij(ia)%mat,2)
          do i1 = 1,size(ad%o%wij(ia)%mat,1)
            i = i + 1
            c1d(i) = ad%o%wij(ia)%mat(i1,i2)
          end do
          end do
        end do

        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::extract_density_ad")) continue

      end subroutine

      subroutine insert_density_ad(c1d,ad)
!doc$ subroutine insert_density(c1d,ad)
        complex(double), dimension(:), pointer :: c1d
        type(atomic_density_paw_obj) :: ad
!       modifies: ad
!       requires: size(c1d) = size(ad%o%wij(*)%mat)
!       effects: Reshapes c1d into ad%o%wij.

!cod$
        integer :: i, i1, i2, ia, ip, ipt, n1d, n1d_all
        integer, dimension(:), allocatable :: count, disp
        complex(double), dimension(:), allocatable :: c1d_all

        call my(ad)

        n1d = size(c1d)

        n1d_all = 0
        do ia = 1,size(ad%o%wij)
          n1d_all = n1d_all + size(ad%o%wij(ia)%mat)
        end do

        allocate( c1d_all(n1d_all) )
        allocate( count(0:mpi_nprocs(SGROUP)-1), disp(0:mpi_nprocs(SGROUP)-1) )

        call allgather(SGROUP,n1d,count)
        disp = 0
        do ip = 0,(mpi_nprocs(SGROUP)-1)
          do ipt = 0,(ip-1)
            disp(ip) = disp(ip) + count(ipt)
          end do
        end do
        call allgatherv(SGROUP,c1d,c1d_all,count,disp)

        call own_i(ad)
        ad%o%g = x_ghost()
        i = 0
        do ia = 1,size(ad%o%wij)
          do i2 = 1,size(ad%o%wij(ia)%mat,2)
          do i1 = 1,size(ad%o%wij(ia)%mat,1)
            i = i + 1
            ad%o%wij(ia)%mat(i1,i2) = c1d_all(i)
          end do
          end do
        end do

        call form_wijs_i(ad%o)   ; if (error()) goto 100
        call form_qlm_i(ad%o)    ; if (error()) goto 100
        call form_energy_i(ad%o)

100     if (allocated( count )) deallocate( count )
        if (allocated( disp )) deallocate( disp )
        if (allocated( c1d_all )) deallocate( c1d_all )

        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::insert_density_ad")) continue

      end subroutine

      subroutine get_normalization_ad(ad,ne)
!doc$ subroutine get_normalization(ad,ne)
        type(atomic_density_paw_obj) :: ad
        real(double) :: ne
!       effects: Returns the number of electrons in ad.

!cod$
        integer :: ia
        complex(double), dimension(:,:), pointer :: a_oij

        call my(ad)

        ne = 0.0_double
        do ia = 1,size(ad%o%wij)
          call atom_c_oij(ad%o%ao,ia,a_oij)
          ne = ne + real(sum(ad%o%wij(ia)%mat*a_oij),double)
        end do

        nullify( a_oij )

        call glean(thy(ad))

      end subroutine

      subroutine write_restart_ad(ad,nrestf)
!doc$ subroutine write_restart(ad,nrestf)
        type(atomic_density_paw_obj) :: ad
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ad restart information to nrestf.

!cod$
        integer :: i, i1, i2, ia, msg, n, na, nsg
        integer(long) :: dsize, iosl, ndata, s4
        integer(long), dimension(:), allocatable :: v4
        complex(double), dimension(:), allocatable :: c1, c2

        call my(ad)
        call my(nrestf)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        ! Start the ATOMIC_DENSITY block
        if (i_access(nrestf)) call startblock(nrestf,"ATOMIC_DENSITY")

        ! Write the PAW tag
        if (i_access(nrestf)) call writetag(nrestf,"PAW")

        ! Write the charge state
        if (i_access(nrestf)) then
          call writetag(nrestf,"CHARGE_STATE")
          dsize = sizeof_double
          ndata = 1
          call writef(ad%o%charge_state,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Write the number of atoms
        if (i_access(nrestf)) then
          call writetag(nrestf,"NUMBER_OF_ATOMS")
          na = x_n_atoms(ad%o%ao)
          dsize = sizeof_long
          ndata = 1
          s4 = na
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Write the numbers of projectors
        if (i_access(nrestf)) then
          call writetag(nrestf,"NUMBER_OF_PROJECTORS")
          ndata = na
          allocate( v4(na) )
          do ia = 1,na
            v4(ia) = size(ad%o%wij(ia)%mat,1)
          end do
          call writef(v4,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Reshape the matrix elements into a 1-D array
        n = 0
        do ia = 1,size(ad%o%wij)
          n = n + size(ad%o%wij(ia)%mat)
        end do
        allocate( c1(n) )
        c1 = (0.0_double,0.0_double)
        i = 0
        do ia = 1,size(ad%o%wij)
          do i2 = 1,size(ad%o%wij(ia)%mat,2)
          do i1 = 1,size(ad%o%wij(ia)%mat,1)
            i = i + 1
            c1(i) = ad%o%wij(ia)%mat(i1,i2)
          end do
          end do
        end do

        ! Write the spin group 1 data (residing on the FILE_SCOPE rank 0 process)
        if (i_access(nrestf)) then
          call writetag(nrestf,"WIJ_MATRIX_ELEMENTS")
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

        ! End the ATOMIC_DENSITY block
        if (i_access(nrestf)) call endblock(nrestf)

100     if (allocated( v4 )) deallocate( v4 )
        if (allocated( c1 )) deallocate( c1 )
        if (allocated( c2 )) deallocate( c2 )

        call glean(thy(ad))
        call glean(thy(nrestf))

        if (error("Exit atomic_density_paw_mod::write_restart_ad")) continue
      
      end subroutine

! private routines

      subroutine form_wij_i(adr,zero)
        type(atomic_density_paw_rep) :: adr
        logical :: zero

        integer :: ia, ib, it, l, m, na, nb, np, o
        real(double) :: occ, spin_factor

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        na = x_n_atoms(adr%ao)
        allocate( adr%wij(na) )

        do ia = 1,na

          np = x_n_atom_projectors(adr%ao,ia)
          allocate( adr%wij(ia)%mat(np,np) )
          adr%wij(ia)%mat = (0.0_double,0.0_double)

          if (zero) cycle

          it = x_atom_type(adr%ao,ia)

          nb = x_n_basis(adr%ao,it)
          do ib = 1,nb
            l = x_basis_l(adr%ao,it,ib)
            o = x_basis_offset(adr%ao,it,ib)
            occ = spin_factor*(x_basis_occupation(adr%ao,it,ib))/real(2*l+1,double)
            do m = -l,+l
              adr%wij(ia)%mat(o+m,o+m) = cmplx(occ,0,double)
            end do
          end do

        end do

        if (error("Exit atomic_density_paw_mod::form_wij_i")) continue

      end subroutine

      subroutine form_wijs_i(adr)
        type(atomic_density_paw_rep) :: adr

        integer :: i, i1, i2, ia, n, n1, n2, na
        complex(double), dimension(:), allocatable :: c1, c2
        
        select case (mpi_nsgroups())
        case (1)

          adr%wijs => adr%wij

        case (2)

          if (.not.associated( adr%wijs )) then
            na = size(adr%wij)
            allocate( adr%wijs(na) )
            do ia = 1,na
              n1 = size(adr%wij(ia)%mat,1)
              n2 = size(adr%wij(ia)%mat,2)
              allocate( adr%wijs(ia)%mat(n1,n2) )
            end do
          end if

          n = 0
          do ia = 1,size(adr%wij)
            n = n + size(adr%wij(ia)%mat)
          end do
          allocate( c1(n), c2(n) )

          i = 0
          do ia = 1,size(adr%wij)
            do i2 = 1,size(adr%wij(ia)%mat,2)
            do i1 = 1,size(adr%wij(ia)%mat,1)
              i = i + 1
              c1(i) = adr%wij(ia)%mat(i1,i2)
            end do
            end do
          end do

          call xcomm_allreduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100

          i = 0
          do ia = 1,size(adr%wijs)
            do i2 = 1,size(adr%wijs(ia)%mat,2)
            do i1 = 1,size(adr%wijs(ia)%mat,1)
              i = i + 1
              adr%wijs(ia)%mat(i1,i2) = c2(i)
            end do
            end do
          end do

100       if (allocated( c1 )) deallocate( c1 )
          if (allocated( c2 )) deallocate( c2 )

        end select

        if (error("Exit atomic_density_paw_mod::form_wijs_i")) continue

      end subroutine 

      subroutine form_qlm_i(adr) ! PRB 55, 2005 (1997): Equation A22.
        type(atomic_density_paw_rep) :: adr

        integer :: i, ia, it, n1, n2, na, qb
        integer :: bi, bj, l, m, mi, mj, ml, omi, omj
        type(paw_data_obj), pointer :: pawd
        
        if (.not.associated( adr%qlm )) then
          na = size(adr%wijs)
          allocate( adr%qlm(na) )
          do ia = 1,na
            it = x_atom_type(adr%ao,ia)
            ml = x_type_l_max(adr%ao,it)
            n1 = 2*ml + 1
            n2 = 2*(2*ml) + 1
            allocate( adr%qlm(ia)%mat(n1,n2) )
          end do
        end if

        do ia = 1,size(adr%qlm)
          it = x_atom_type(adr%ao,ia)
          call type_projector_data(adr%ao,it,pawd)
          ml = x_type_l_max(adr%ao,it)
          qb = 2*ml + 1
          adr%qlm(ia)%mat = (0.0_double,0.0_double)
          adr%qlm(ia)%mat(1,qb) = cmplx(x_qeffion(pawd),0,double)
          do i = 1,x_qvlm_size(pawd)
            call orbital_decode(pawd,i,bi,bj,mi,mj,l)
            m = mj - mi
            omi = x_basis_offset(adr%ao,it,bi) + mi
            omj = x_basis_offset(adr%ao,it,bj) + mj
            adr%qlm(ia)%mat(l+1,qb+m) = adr%qlm(ia)%mat(l+1,qb+m) + x_aqlm(pawd,i)*adr%wijs(ia)%mat(omi,omj)
          end do
        end do

        nullify( pawd )

        if (error("Exit atomic_density_paw_mod::form_qlm_i")) continue

      end subroutine 

      subroutine form_energy_i(adr) 
        type(atomic_density_paw_rep) :: adr
!       notes: Unless an "atomic_symmetry  off" line is present in argvf, this routine explicitly computes energy
!              contributions only for atoms that are unique with respect to symmetry. Energy contributions from
!              symmetry-related atoms are then accumulated via symmetry operations.

        integer :: nu
        real(double), dimension(:), allocatable :: energy

        nu = x_n_unique_atoms(adr%ao)
        allocate( energy(nu) )
        energy = 0.0_double

        call diagonal_terms_i(adr,energy) ; if (error()) goto 100  
        call hartree_terms_i(adr,energy)  ; if (error()) goto 100
        call xc_terms_i(adr,energy)       ; if (error()) goto 100
        call hat_terms_i(adr,energy)      ; if (error()) goto 100

        call accumulate_energy_i(adr,energy)

100     if (allocated( energy )) deallocate( energy )

        if (error("Exit atomic_density_paw_mod::form_energy_i")) continue

      end subroutine 

      subroutine diagonal_terms_i(adr,energy) ! PRB 55, 2005 (1997): Equations A10, A11, A12, and A31(1,2,3).
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy

        integer :: i, ia, it, m
        integer :: bi, bj, li, lj, oi, oj, omi, omj
        integer :: iu, nu, nu_proc, first_u, last_u
        real(double) :: fn
        type(paw_data_obj), pointer :: pawd 

        nu = x_n_unique_atoms(adr%ao)
        call subdivide(mpi_myproc(CONFIG),mpi_nprocs(CONFIG),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u
          ia = x_unique_atom(adr%ao,iu)
          it = x_unique_type(adr%ao,iu)
          call type_projector_data(adr%ao,it,pawd)
          i = 0 
          do bi = 1,x_n_basis(adr%ao,it)
            li = x_basis_l(adr%ao,it,bi)
            do bj = 1,bi 
              lj = x_basis_l(adr%ao,it,bj) 
              if (lj == li) then 
                oi = x_basis_offset(adr%ao,it,bi)
                oj = x_basis_offset(adr%ao,it,bj)
                i = i + 1 
                fn = x_kinetic(pawd,i) + x_v_ion(pawd,i) 
                do m = -lj,+lj 
                  omi = oi + m
                  omj = oj + m
                  energy(iu) = energy(iu) + fn*real(adr%wijs(ia)%mat(omi,omj),double)
                  if (bj /= bi) then 
                    energy(iu) = energy(iu) + fn*real(adr%wijs(ia)%mat(omj,omi),double)
                  end if
                end do
              end if
            end do
          end do
        end do

        nullify( pawd )

        if (error("Exit atomic_density_paw_mod::diagonal_terms_i")) continue

      end subroutine

      subroutine hartree_terms_i(adr,energy) ! PRB 55, 2005 (1997): Equations A26 and A31 (6th term).
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy
  
        integer :: i, ia, it
        integer :: bi, bj, bk, bl, mi, mj, mk, ml, omi, omj, omk, oml
        integer :: iu, nu, nu_proc, first_u, last_u
        complex(double) :: c1
        type(paw_data_obj), pointer :: pawd

        nu = x_n_unique_atoms(adr%ao)
        call subdivide(mpi_myproc(CONFIG),mpi_nprocs(CONFIG),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u
          ia = x_unique_atom(adr%ao,iu)
          it = x_unique_type(adr%ao,iu)
          call type_projector_data(adr%ao,it,pawd)
          do i = 1,x_cijkl_size(pawd)
            call cijkl_decode(pawd,i,bi,bj,bk,bl,mi,mj,mk)
            ml = mi - mj + mk
            omi = x_basis_offset(adr%ao,it,bi) + mi
            omj = x_basis_offset(adr%ao,it,bj) + mj
            omk = x_basis_offset(adr%ao,it,bk) + mk
            oml = x_basis_offset(adr%ao,it,bl) + ml
            c1 = x_cijkl(pawd,i)*adr%wijs(ia)%mat(omk,oml)
            energy(iu) = energy(iu) + 0.5_double*real(c1*adr%wijs(ia)%mat(omi,omj),double)
          end do

        end do
        
        nullify( pawd )

        if (error("Exit atomic_density_paw_mod::hartree_terms_i")) continue

      end subroutine

      subroutine xc_terms_i(adr,energy)
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy

        if (uses_gradient(x_axc(adr%ao))) then
          call accum_xc_grad_i(adr,energy)
        else
          call accum_xc_i(adr,energy)
        end if

        if (error("Exit atomic_density_paw_mod::xc_terms_i")) continue

      end subroutine

      subroutine accum_xc_i(adr,energy)
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy

        integer :: ia, it, iu, iv, nb, nr, nv
        integer :: bi, bj, li, lj, mi, mj, oi, oj, osi, osj
        real(double) :: spin_factor, vwt
        real(double), dimension(:), pointer :: n, del_n, lap_n, exc
        real(double), dimension(:), pointer :: nt, del_nt, lap_nt, exct
        real(double), dimension(:), pointer :: r2, wt, nc, nct
        real(double), dimension(:,:,:), pointer :: phir2ij, tphir2ij
        complex(double) :: cij
        complex(double), dimension(:,:), pointer :: ylmij
        type(paw_data_obj), pointer :: pawd

        nullify( n, del_n, lap_n, exc )
        nullify( nt, del_nt, lap_nt, exct )

        allocate( n(0), exc(0) )
        allocate( nt(0), exct(0) )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        ! Compute energy contributions for symmetry-unique atoms.
        nv = x_n_vectors(adr%ao)
        do iv = 1,nv

          iu = x_vector_unique_atom(adr%ao,iv)
          ia = x_vector_atom(adr%ao,iv)
          it = x_vector_type(adr%ao,iv)

          call type_projector_data(adr%ao,it,pawd)

          nr = x_grid_size(pawd)
          if (size(n) /= nr) then
            deallocate( n, exc ) ; allocate( n(nr), exc(nr) )
            deallocate( nt, exct ) ; allocate( nt(nr), exct(nr) )
          end if
          n = 0.0_double
          nt = 0.0_double

          call type_r2(pawd,r2)
          call type_wt(pawd,wt)
          call type_core_density(pawd,nc)
          call type_coretail_density(pawd,nct)
          call type_phir2ij(pawd,phir2ij)
          call type_tphir2ij(pawd,tphir2ij)

          nb = x_n_basis(adr%ao,it)
          call vector_ylmij(adr%ao,iv,ylmij)
          vwt = x_vector_weight(adr%ao,iv)

          do bj = 1,nb
            lj = x_basis_l(adr%ao,it,bj)
            oj = x_basis_offset(adr%ao,it,bj)
            osj = lj**2 + lj + 1
            do bi = 1,nb
              li = x_basis_l(adr%ao,it,bi)
              oi = x_basis_offset(adr%ao,it,bi)
              osi = li**2 + li + 1
              cij = (0.0_double,0.0_double)
              do mj = -lj,+lj
              do mi = -li,+li
                cij = cij + adr%wij(ia)%mat(oi+mi,oj+mj)*ylmij(osi+mi,osj+mj)
              end do
              end do
              n = n + real(cij,double)*phir2ij(:,bi,bj)
              nt = nt + real(cij,double)*tphir2ij(:,bi,bj)
            end do
          end do
          n = n + spin_factor*nc
          nt = nt + spin_factor*nct

          call xc_energy_density(x_axc(adr%ao),n,del_n,lap_n,exc)
          call xc_energy_density(x_axc(adr%ao),nt,del_nt,lap_nt,exct)

          energy(iu) = energy(iu) + sum(r2*(exc*n - exct*nt)*wt)*vwt

        end do

        if (associated( n )) then
          deallocate( n, exc )
          deallocate( nt, exct )
        end if
        nullify( r2 )
        nullify( wt)
        nullify( nc )
        nullify( nct )
        nullify( phir2ij )
        nullify( tphir2ij )
        nullify( pawd )
        nullify( ylmij )

        if (error("Exit atomic_density_paw_mod::accum_xc_i")) continue

      end subroutine 

      subroutine accum_xc_grad_i(adr,energy)
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy

        integer :: ia, ir, it, iu, iv, nb, nr, nv
        integer :: bi, bj, li, lj, mi, mj, oi, oj, osi, osj
        real(double) :: spin_factor, vwt
        real(double), dimension(:), allocatable :: dn_dr, dn_dt, dn_dp
        real(double), dimension(:), allocatable :: dnt_dr, dnt_dt, dnt_dp
        real(double), dimension(:), pointer :: n, del_n, lap_n, exc
        real(double), dimension(:), pointer :: nt, del_nt, lap_nt, exct
        real(double), dimension(:), pointer :: r2, wt, nc, dnc_dr, nct, dnct_dr
        real(double), dimension(:,:,:), pointer :: phir2ij, dphir2ij_dr, tphir2ij, dtphir2ij_dr
        complex(double) :: cij, dcij_dt, dcij_dp
        complex(double), dimension(:,:), pointer :: ylmij, dylmij_dt, dylmij_dp
        type(paw_data_obj), pointer :: pawd

        nullify( n, del_n, lap_n, exc )
        nullify( nt, del_nt, lap_nt, exct )

        allocate( n(0), dn_dr(0), dn_dt(0), dn_dp(0), del_n(0) )
        allocate( exc(0) )
        allocate( nt(0), dnt_dr(0), dnt_dt(0), dnt_dp(0), del_nt(0) )
        allocate( exct(0) )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        ! Compute energy contributions for symmetry-unique atoms.
        nv = x_n_vectors(adr%ao)
        do iv = 1,nv

          iu = x_vector_unique_atom(adr%ao,iv)
          ia = x_vector_atom(adr%ao,iv)
          it = x_vector_type(adr%ao,iv)

          call type_projector_data(adr%ao,it,pawd)

          nr = x_grid_size(pawd)
          if (size(n) /= nr) then
            deallocate( n, dn_dr, dn_dt, dn_dp, del_n ) ; allocate( n(nr), dn_dr(nr), dn_dt(nr), dn_dp(nr), del_n(nr) )
            deallocate( exc ) ; allocate( exc(nr) )
            deallocate( nt, dnt_dr, dnt_dt, dnt_dp, del_nt ) ; allocate( nt(nr), dnt_dr(nr), dnt_dt(nr), dnt_dp(nr), del_nt(nr) )
            deallocate( exct ) ; allocate( exct(nr) )
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

          nb = x_n_basis(adr%ao,it)
          call vector_ylmij(adr%ao,iv,ylmij)
          call vector_dylmij_dt(adr%ao,iv,dylmij_dt)
          call vector_dylmij_dp(adr%ao,iv,dylmij_dp)
          vwt = x_vector_weight(adr%ao,iv)

          do bj = 1,nb
            lj = x_basis_l(adr%ao,it,bj)
            oj = x_basis_offset(adr%ao,it,bj)
            osj = lj**2 + lj + 1
            do bi = 1,nb
              li = x_basis_l(adr%ao,it,bi)
              oi = x_basis_offset(adr%ao,it,bi)
              osi = li**2 + li + 1
              cij = (0.0_double,0.0_double)
              dcij_dt = (0.0_double,0.0_double)
              dcij_dp = (0.0_double,0.0_double)
              do mj = -lj,+lj
              do mi = -li,+li
                cij = cij + adr%wij(ia)%mat(oi+mi,oj+mj)*ylmij(osi+mi,osj+mj)
                dcij_dt = dcij_dt + adr%wij(ia)%mat(oi+mi,oj+mj)*dylmij_dt(osi+mi,osj+mj)
                dcij_dp = dcij_dp + adr%wij(ia)%mat(oi+mi,oj+mj)*dylmij_dp(osi+mi,osj+mj)
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

          call xc_energy_density(x_axc(adr%ao),n,del_n,lap_n,exc)
          call xc_energy_density(x_axc(adr%ao),nt,del_nt,lap_nt,exct)

          energy(iu) = energy(iu) + sum(r2*(exc*n - exct*nt)*wt)*vwt

        end do

        if (associated( n )) then
          deallocate( n, dn_dr, dn_dt, dn_dp, del_n )
          deallocate( exc )
          deallocate( nt, dnt_dr, dnt_dt, dnt_dp, del_nt )
          deallocate( exct )
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
        nullify( pawd )
        nullify( ylmij )
        nullify( dylmij_dt )
        nullify( dylmij_dp )

        if (error("Exit atomic_density_paw_mod::accum_xc_grad_i")) continue

      end subroutine 
  
      subroutine hat_terms_i(adr,energy) ! PRB 55, 2005 (1997): Equations A23 and A31 (4th term).
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy
  
        integer :: i, ia, it, l, m, ml, qb
        integer :: bi, bj, mi, mj, omi, omj
        integer :: iu, nu, nu_proc, first_u, last_u
        real(double) :: r1
        complex(double) :: c1
        type(paw_data_obj), pointer :: pawd

        nu = x_n_unique_atoms(adr%ao)
        call subdivide(mpi_myproc(CONFIG),mpi_nprocs(CONFIG),1,nu,first_u,last_u,nu_proc)

        do iu = first_u,last_u
          ia = x_unique_atom(adr%ao,iu)
          it = x_unique_type(adr%ao,iu)
          call type_projector_data(adr%ao,it,pawd)
          ml = x_type_l_max(adr%ao,it)
          qb = 2*ml + 1
          do i = 1,x_qvlm_size(pawd)
            call orbital_decode(pawd,i,bi,bj,mi,mj,l)
            m = mi - mj
            omi = x_basis_offset(adr%ao,it,bi) + mi
            omj = x_basis_offset(adr%ao,it,bj) + mj
            c1 = x_avlm(pawd,i)*adr%qlm(ia)%mat(l+1,qb+m)
            energy(iu) = energy(iu) - real(c1*adr%wijs(ia)%mat(omi,omj),double)
          end do
          do l = 0,2*ml
            r1 = x_hat_self_energy(pawd,l)
            do m = -l,+l
              energy(iu) = energy(iu) - r1*real(adr%qlm(ia)%mat(l+1,qb+m)*conjg(adr%qlm(ia)%mat(l+1,qb+m)),double)
            end do
          end do
          energy(iu) = energy(iu) - x_coretail_selfenergy(pawd) - real(adr%qlm(ia)%mat(1,qb),double)*x_coretail_hatenergy(pawd)
        end do

        nullify( pawd )

        if (error("Exit atomic_density_paw_mod::hat_terms_i")) continue

      end subroutine

      subroutine own_i(ad)
        type(atomic_density_paw_obj) :: ad
        type(atomic_density_paw_obj) :: adt

        integer :: i
        if (ad%ref < ad%o%ref) then
          allocate( adt%o )
          adt%o%ref = 0
          adt%o%g = ad%o%g
          adt%o%charge_state = ad%o%charge_state
          adt%o%energy = ad%o%energy
          call my(ad%o%ao,adt%o%ao)
          allocate( adt%o%wij(size(ad%o%wij)) )
          do i = 1,size(ad%o%wij)
            allocate( adt%o%wij(i)%mat(size(ad%o%wij(i)%mat,1),size(ad%o%wij(i)%mat,2)) )
            adt%o%wij(i)%mat = ad%o%wij(i)%mat
          end do
          select case (mpi_nsgroups())
          case (1)
            adt%o%wijs => adt%o%wij
          case (2)
            if (associated( ad%o%wijs )) then
              allocate( adt%o%wijs(size(ad%o%wijs)) )
              do i = 1,size(ad%o%wijs)
                allocate( adt%o%wijs(i)%mat(size(ad%o%wijs(i)%mat,1),size(ad%o%wijs(i)%mat,2)) )
                adt%o%wijs(i)%mat = ad%o%wijs(i)%mat
              end do
            else
              nullify( adt%o%wijs )
            end if
          end select
          if (associated( ad%o%qlm )) then
            allocate( adt%o%qlm(size(ad%o%qlm)) )
            do i = 1,size(ad%o%qlm)
              allocate( adt%o%qlm(i)%mat(size(ad%o%qlm(i)%mat,1),size(ad%o%qlm(i)%mat,2)) )
              adt%o%qlm(i)%mat = ad%o%qlm(i)%mat
            end do
          else
            nullify( adt%o%qlm )
          end if
          ad%o%ref = ad%o%ref - ad%ref
          ad%o => adt%o
          ad%o%ref = ad%o%ref + ad%ref
        end if
      end subroutine

      subroutine accumulate_energy_i(adr,energy_pu)
        type(atomic_density_paw_rep) :: adr
        real(double), dimension(:) :: energy_pu

        integer :: ia, iu, na, nu
        real(double), dimension(:), allocatable :: energy_a, energy_u
        type(space_group_obj) :: sg

        call my(x_space_group(adr%ao),sg)

        na = x_n_atoms(adr%ao)
        nu = x_n_unique_atoms(adr%ao)

        ! Sum the processor energy contributions.
        allocate( energy_u(nu) )
        call allreduce(CONFIG,MPI_SUM,energy_pu,energy_u)

        ! Distribute the energy values to non-symmetry-unique atoms.
        allocate( energy_a(na) )
        do iu = 1,nu
          ia = x_unique_atom(adr%ao,iu)
          energy_a(ia) = energy_u(iu)
        end do
        if (nu < na) then
          do iu = 1,nu
            ia = x_unique_atom(adr%ao,iu)
            call distribute_energy(sg,ia,energy_a)
          end do
        end if

        ! Sum the atom energy contributions.
        adr%energy = sum(energy_a)

        if (allocated( energy_a )) deallocate( energy_a )
        if (allocated( energy_u )) deallocate( energy_u )

        call glean(thy(sg))

        if (error("Exit atomic_density_paw_mod::accumulate_energy_i")) continue

      end subroutine

      end module

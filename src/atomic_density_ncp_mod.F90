!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module atomic_density_ncp_mod
!doc$ module atomic_density_ncp_mod

      use kind_mod
      use path_mod
      use mpi_mod
      use diary_mod
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
      use ewald_mod
      use ncp_data_mod
      use atomic_operators_ncp_mod

!     One datatype is available here: type(atomic_density_ncp_obj).

!cod$
      implicit none
      private

      type :: gp_mat
        complex(double), dimension(:,:), pointer :: mat
      end type

      type :: atomic_density_ncp_rep
        integer :: ref
        type(ghost) :: g                               ! ghost
        real(double) :: charge_state                   ! charge state
        real(double) :: energy                         ! total wij energy (sum of sgroup contributions)
        type(gp_mat), dimension(:), pointer :: wij     ! projected occupation coefficients for an sgroup
        type(atomic_operators_ncp_obj) :: ao           ! atomic operators object
      end type

      type, public :: atomic_density_ncp_obj
        private
        integer :: ref
        type(atomic_density_ncp_rep), pointer :: o
      end type

!doc$
      public :: atomic_density_ncp
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_atomic_operators
!      public :: symmetrize
      public :: merge_atomic_density
      public :: add_atomic_density
      public :: atomic_hartree_density
      public :: guess_density
      public :: atomic_energy
      public :: atomic_forces
      public :: atomic_pressure
      public :: atomic_stress_tensor
      public :: write_restart

!cod$
      interface atomic_density_ncp
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
      interface x_ref
        module procedure ad_ref
      end interface
      interface x_ghost
        module procedure ad_ghost
      end interface
      interface x_atomic_operators
        module procedure ad_atomic_operators
      end interface
 !     interface symmetrize
  !       module procedure symmetrize_ad
   !   end interface
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
      interface atomic_pressure
        module procedure atomic_pressure_ad
      end interface
      interface atomic_stress_tensor
        module procedure atomic_stress_tensor_ad
      end interface
      interface write_restart
        module procedure write_restart_ad
      end interface

      contains

! public routines

      function constructor_ad(ao,restf) result(ad)
!doc$ function atomic_density_ncp(ao,restf) result(ad)
        type(atomic_operators_ncp_obj) :: ao
        type(tagio_obj), optional :: restf
        type(atomic_density_ncp_obj) :: ad
!       effects: Constructs a new ad with wij and energy set to 0.
!       errors: Passes errors.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: tag
        integer :: np, na, ia
        integer, dimension(2) :: csr
        integer(long) :: dsize, iosl, ndata

        call my(ao)
        if (present(restf)) call my(restf)

        ad%ref = 0
        allocate( ad%o )
        ad%o%ref = 0
        ad%o%g = x_ghost()

        call my(ao,ad%o%ao)

        ! open the ATOMIC_DENSITY block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_DENSITY")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_DENSITY block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)
        end if

        ! find the NCP tag
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"NCP")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NCP tag was not found")) goto 200
        end if

        ! read the charge state
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"CHARGE_STATE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: CHARGE_STATE tag was not found")) goto 200
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
            call readf(ad%o%charge_state,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,ad%o%charge_state)
        else
          call arglc("charge_state_mode",tag,found)
          if (.not.found) tag = "real_number"
          select case (trim(tag))
          case ("real_number")
            call arg("charge_state",ad%o%charge_state,found)
            if (.not.found) ad%o%charge_state = 0.0_double
            !if ((mpi_nconfigs() == 2) .and. (mpi_myconfig() == 2)) then
            !  call arglc("ts_method",tag,found)
            !  if (.not.found) tag = "none"
            !  select case (trim(tag))
            !  case ("electron_capture")
            !    ad%o%charge_state = ad%o%charge_state - 1.0_double
            !  case ("hole_capture")
            !    ad%o%charge_state = ad%o%charge_state + 1.0_double
            !  end select
            !end if
          case ("integer_ratio")
            call arg("charge_state_ratio",csr,found)
            if (error(.not.found,"ERROR: charge_state_ratio was not found")) goto 200
            if (error(csr(2) == 0,"ERROR: denominator = 0")) goto 200
            ad%o%charge_state = real(csr(1),double)/real(csr(2),double)
          case default
            if (error(.true.,"ERROR: charge_state_mode was not recognized")) goto 200
          end select
        end if

        ! initialize wij and the energy
        na = x_n_atoms(ad%o%ao)
        allocate ( ad%o%wij(na) )
        do ia = 1,na
           np = x_n_atom_projectors(ad%o%ao, ia)
           allocate ( ad%o%wij(ia)%mat(np,np) )
           ad%o%wij(ia)%mat = 0.0_double
        end do
        ad%o%energy = 0.0_double

        ! close the ATOMIC_DENSITY block
200     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

300     call glean(thy(ao))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_density_ncp_mod::constructor_ad")) continue

      end function

      subroutine update_ad(ad,ao)
!doc$ subroutine update(ad,ao)
        type(atomic_density_ncp_obj) :: ad
        type(atomic_operators_ncp_obj) :: ao
!       modifies: ad
!       effects: Updates ad.

!cod$
        logical :: ao_change

        call my(ad)
        call my(ao)

        ao_change = ( x_ghost(ad%o%ao) /= x_ghost(ao) )
        if (ao_change) then
          call own_i(ad)
          ad%o%g = x_ghost()
          ad%o%ao = ao
        end if

        call glean(thy(ao))
        call glean(thy(ad))

        if (error("Exit atomic_density_ncp_mod::update_ad")) continue

      end subroutine

      subroutine my_ad(ad)
!doc$ subroutine my(ad)
        type(atomic_density_ncp_obj) :: ad

!cod$
        ad%ref = ad%ref + 1
        ad%o%ref = ad%o%ref + 1
      end subroutine

      subroutine my_new_ad(adi,ad)
!doc$ subroutine my(adi,ad)
        type(atomic_density_ncp_obj) :: adi
        type(atomic_density_ncp_obj) :: ad

!cod$
        ad%ref = 1
        ad%o => adi%o
        ad%o%ref = ad%o%ref + 1
      end subroutine

      function thy_ad(ad) result(ado)
!doc$ function thy(ad) result(ado)
        type(atomic_density_ncp_obj) :: ad, ado

!cod$
        ad%ref = ad%ref - 1
        ad%o%ref = ad%o%ref - 1
        ado%ref = ad%ref
        ado%o => ad%o
      end function

      subroutine glean_ad(ad)
!doc$ subroutine glean(ad)
        type(atomic_density_ncp_obj) :: ad

!cod$
        integer :: i

        if (ad%o%ref < 1) then
           if (associated( ad%o%wij )) then
            do i = 1,size(ad%o%wij)
              if (associated( ad%o%wij(i)%mat )) deallocate( ad%o%wij(i)%mat )
            end do
            deallocate( ad%o%wij )
          end if
!          deallocate( ad%o )

        end if
      end subroutine

      subroutine bequeath_ad(ad)
!doc$ subroutine bequeath(ad)
        type(atomic_density_ncp_obj) :: ad

!cod$
        continue
      end subroutine
    
      subroutine assign_ad(ad,ad2)
!doc$ subroutine assignment(=)(ad,ad2)
        type(atomic_density_ncp_obj), intent(inout) :: ad
        type(atomic_density_ncp_obj), intent(in) :: ad2

!cod$
        type(atomic_density_ncp_obj) :: adt
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
        type(atomic_density_ncp_obj) :: ad
        integer, dimension(2) :: r
!       effects: Returns ad%ref and ad%o%ref.
      
!cod$ 
        r(1) = ad%ref
        r(2) = ad%o%ref
        call glean(ad)
      end function 

      function ad_ghost(ad) result(g)
!doc$ function x_ghost(ad) result(g)
        type(atomic_density_ncp_obj) :: ad
        type(ghost) :: g
!       effects: Returns the ghost of ad.

!cod$
        call my(ad)
        g = ad%o%g
        call glean(thy(ad))
      end function

      function ad_atomic_operators(ad) result(ao)
!doc$ function x_atomic_operators(ad) result(ao)
        type(atomic_density_ncp_obj) :: ad
        type(atomic_operators_ncp_obj) :: ao
!       effects: Returns the ao of ad.

!cod$
        call my(ad)
        call my(ad%o%ao,ao)
        call bequeath(thy(ao))
        call glean(thy(ad))
      end function


!      subroutine symmetrize_ad(ad,sg)
!doc$ subroutine symmetrize(ad,sg)
!        type(atomic_density_paw_obj) :: ad
!        type(space_group_obj) :: sg
!       modifies: ad
!       effects: Symmetrizes ad with respect to sg.
!       errors: Passes errors.

!cod$
!        logical, dimension(:), allocatable :: mask
!        integer :: ia, it, na, nb, nt
!        integer :: bi, bj, li, lj, mi, mj, oi, oj
!        complex(double), dimension(:,:,:), allocatable :: w_block

!        call my(ad)
!        call my(sg)

!        call own_i(ad)
!        ad%o%g = x_ghost()

!        na = x_n_atoms(ad%o%ao)
!        nt = x_n_types(ad%o%ao)
!        allocate( mask(na) )

!        do it = 1,nt
!          nb = x_n_basis(ad%o%ao,it)
!          do bi = 1,nb
!            li = x_basis_l(ad%o%ao,it,bi)
!            oi = x_basis_offset(ad%o%ao,it,bi)
!            do bj = 1,nb
!              lj = x_basis_l(ad%o%ao,it,bj)
!              oj = x_basis_offset(ad%o%ao,it,bj)
!              allocate( w_block(-li:+li,-lj:+lj,na))
!              mask = .false.
!              do ia = 1,na
!                if (x_atom_type(ad%o%ao,ia) == it) then
!                  mask(ia) = .true.
!                  do mj = -lj,+lj
!                  do mi = -li,+li
!                    w_block(mi,mj,ia) = ad%o%wij(ia)%mat(oi+mi,oj+mj)
!                  end do
!                  end do
!                else
!                  w_block(:,:,ia) = (0.0_double,0.0_double)
!                end if
!              end do
!              call symmetrize_spherical_tensor(sg,w_block,mask) ; if (error()) goto 100
!              do ia = 1,na
!                if (mask(ia)) then
!                  do mj = -lj,+lj
!                  do mi = -li,+li
!                    ad%o%wij(ia)%mat(oi+mi,oj+mj) = w_block(mi,mj,ia)
!                  end do
!                  end do
!                end if
!              end do
!              deallocate( w_block )
!            end do
!          end do
!        end do

!        if (allocated( w_block )) deallocate( w_block )
!        if (allocated( mask )) deallocate( mask )

!        call form_wijs_i(ad%o)    ; if (error()) goto 100
!        call form_qlm_i(ad%o)     ; if (error()) goto 100
!        call form_energy_i(ad%o)
!
!100     call glean(thy(ad))
!        call glean(thy(sg))

!        if (error("Exit atomic_density_paw_mod::symmetrize_ad")) continue

!      end subroutine


      subroutine merge_atomic_density_ad(ad)
!doc$ subroutine merge_atomic_density(ad)
        type(atomic_density_ncp_obj) :: ad
!       modifies: ad
!       effects: Merges ad from different kgroups and computes the energy.
!       errors: Passes errors.

!cod$
        integer :: i, i1, i2, ia, first_a, last_a, n, na, na_proc, ij
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

        ! do ia = 1,size(ad%o%wij)
        !   do ij = 1,size(ad%o%wij(ia)%mat,1)
        !      write(x_unit(diaryfile()), '(/,"Wii element, " , i3, " is ", f15.10)') ij, ad%o%wij(ia)%mat(ij,ij) ; if (error()) goto 100
        !   end do
        !end do


100     if (allocated( c1 )) deallocate( c1 )
        if (allocated( c2 )) deallocate( c2 )

        call form_energy_i(ad%o)

!        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"energy in atomic density is", f15.10)') ad%o%energy
        
        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::merge_atomic_density_ad")) continue

      end subroutine

      subroutine add_atomic_density_ad(ad,pdots,weights)
!doc$ subroutine add_atomic_density(ad,pdots,weights)
        type(atomic_density_ncp_obj) :: ad
        complex(double), dimension(:,:), intent(in) :: pdots
        real(double), dimension(:), intent(in) :: weights
!       modifies: ad
!       effects: Accumulates pdots contributions from different kgroups into ad%o%wii.
!       errors: Passes errors.

!cod$
        character(1), parameter :: transa = 'n', transb = 't'
        integer :: bi, first_a, last_a, ia, ib, it, na, na_proc, nb, nm, np, nt
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: wcpdots
        !complex(double), dimension(:,:), pointer :: t_tr2c, t_tc2r
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
          twij(it)%mat = (0.0_double,0.0_double)
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
          !call type_tr2c(ad%o%ao,it,t_tr2c)
          !call type_tc2r(ad%o%ao,it,t_tc2r)
          call zgemm(transa,transb,nm,nm,nb,alpha,wcpdots(bi,1),np,pdots(bi,1),np,beta,twij(it)%mat,nm)
          ad%o%wij(ia)%mat = ad%o%wij(ia)%mat + twij(it)%mat
        end do

        if (allocated( wcpdots )) deallocate( wcpdots )
        if (associated( twij )) then
          do it = 1,size(twij)
            if (associated( twij(it)%mat )) deallocate( twij(it)%mat )
          end do
          deallocate( twij )
        end if


        call glean(thy(ad))

        if (error("Exit atomic_density_paw_mod::add_atomic_density_ad")) continue

      end subroutine

      function atomic_hartree_density_ad(ad,lay) result(den)
!doc$ function atomic_hartree_density(ad,lay) result(den)
        type(atomic_density_ncp_obj):: ad
        type(layout_obj) :: lay
        type(grid_obj) :: den
!       effects: Returns an empty grid with SGROUP scope.

!cod$
        call my(ad)
        call my(lay)
        call my(grid(lay,SGROUP),den)
        call bequeath(thy(den))
        call glean(thy(ad))
        call glean(thy(lay))
        if (error("Exit atomic_density_ncp_mod::atomic_hartree_density_ad")) continue
      end function

      function guess_density_ad(ad,lay,ne) result(den)
!doc$ function guess_density(ad,lay,ne) result(den)
        type(atomic_density_ncp_obj) :: ad
        type(layout_obj) :: lay
        real(double), optional :: ne
        type(grid_obj) :: den
!       modifies: den
!       effects: Returns a filtered approximation to the valence density.

!cod$
        integer :: i1, i2, i3, ia, it
        real(double) :: ne_g, prefactor, spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        real(double), dimension(:,:,:), pointer :: vdff
        complex(double) :: den_norm, igr
        complex(double), dimension(:,:,:), pointer :: c_den
        type(lattice_obj) :: lat

        call my(ad)
        call my(lay)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)

        nullify( gx, gy, gz, c_den )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        if (present(ne)) then
          ne_g = ne
        else
          ne_g = spin_factor*(x_valence_electrons(ad%o%ao) - ad%o%charge_state)
        end if
        den_norm = cmplx(ne_g,0,double)/x_cell_volume(lat)

        call my(grid(lay,SGROUP),den)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call alloc(c_den,lay,D_TYPE,SGROUP)
        c_den = (0.0_double,0.0_double)
        do it = 1,x_n_types(ad%o%ao)
          call valence_density_ff(ad%o%ao,it,vdff)
          prefactor = spin_factor*x_type_valence(ad%o%ao,it)
          do ia = 1,x_n_atoms(ad%o%ao)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
            pos = lat2r(lat,x_position(x_atoms(x_crystal(ad%o%ao)),ia))
            do i3 = 1,size(gx,3)
            do i2 = 1,size(gx,2)
            do i1 = 1,size(gx,1)
              igr = (0.0_double,1.0_double)*( pos(1)*gx(i1,i2,i3) + pos(2)*gy(i1,i2,i3) + pos(3)*gz(i1,i2,i3) )
              c_den(i1,i2,i3) = c_den(i1,i2,i3) + prefactor*vdff(i1,i2,i3)*exp(-igr)
            end do
            end do
            end do
          end do
        end do
        call put(c_den,den,CDF_KIND)
        call set_normalization(den,den_norm)

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )
        nullify( vdff )

        call glean(thy(lat))
        call bequeath(thy(den))

        call glean(thy(ad))
        call glean(thy(lay))

        if (error("Exit atomic_density_ncp_mod::guess_density_ad")) continue

      end function

      subroutine atomic_energy_ad(ad,e)
!doc$ subroutine atomic_energy(ad,e)
        type(atomic_density_ncp_obj) :: ad
        real(double), intent(out) :: e
!       modifies: e
!       effects: Returns the atomic energy.

!cod$
        call my(ad)
        e = ad%o%energy + x_ewald_energy(ad%o%ao)
        call glean(thy(ad))

      end subroutine

      subroutine atomic_forces_ad(ad,den,xcp,ccd,f)
!doc$ subroutine atomic_forces(ad,den,xcp,ccd,f)
        type(atomic_density_ncp_obj) :: ad
        type(grid_obj) :: den, xcp, ccd
        real(double), dimension(:,:), intent(out) :: f
!       requires: f be dimension(3,x_n_atoms(ad%o%ao)).
!                 Input-grid scopes be SGROUP.
!       modifies: f
!       effects: Returns unsymmetrized atomic contributions to the forces.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, it
        real(double) :: r0, spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:), allocatable :: q
        real(double), dimension(:,:), allocatable :: ft
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        real(double), dimension(:,:,:), pointer :: cdff, gpff
        complex(double) :: igr
        complex(double), dimension(:,:,:), pointer :: c_ccd, c_den, c_xcp, c_sum
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(layout_obj) :: lay

        call my(ad)
        call my(den)
        call my(xcp)
        call my(ccd)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)
        call my(x_atoms(x_crystal(ad%o%ao)),ats)
        call my(x_layout(den),lay)

        nullify( gx, gy, gz, c_ccd, c_den, c_xcp, c_sum )
        nullify( cdff, gpff )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        if (x_type(ccd) /= EMPTY_KIND) then
          call take(c_ccd,ccd,CDF_KIND)
        end if
        call take(c_den,den,CDF_KIND)
        call take(c_xcp,xcp,CDF_KIND)

        call alloc(c_sum,lay,D_TYPE,SGROUP)

        allocate( ft(size(f,1),size(f,2)) )
        ft = 0.0_double

        do it = 1,x_n_types(ad%o%ao)

          call grid_potential_ff(ad%o%ao,it,gpff)
          c_sum = gpff*conjg(c_den)
          if (associated( c_ccd )) then
            c_sum = c_sum + spin_factor*gpff*conjg(c_ccd)
          end if

          call core_density_ff(ad%o%ao,it,cdff)
          if (associated( cdff )) then
            c_sum = c_sum + spin_factor*cdff*conjg(c_xcp)
          end if

          do ia = 1,x_n_atoms(ad%o%ao)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
            pos = lat2r(lat,x_position(ats,ia))
            do i3 = 1,size(gx,3)
            do i2 = 1,size(gx,2)
            do i1 = 1,size(gx,1)
              igr = (0.0_double,1.0_double)*(gx(i1,i2,i3)*pos(1) + gy(i1,i2,i3)*pos(2) + gz(i1,i2,i3)*pos(3))
              r0 = aimag(c_sum(i1,i2,i3)*exp(-igr))
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

        if (associated( c_ccd )) then
          call put(c_ccd,ccd,CDF_KIND)
        end if
        call put(c_den,den,CDF_KIND)
        call put(c_xcp,xcp,CDF_KIND)

        allocate( q(x_n_atoms(ad%o%ao)) )
        do ia = 1,size(q)
          q(ia) = x_atom_valence(ad%o%ao,ia)
        end do
        call ewald_forces(x_crystal(ad%o%ao),q,ft) ; if (error()) goto 100
        f = f + ft

100     if (allocated( q )) deallocate( q )
        if (allocated( ft )) deallocate( ft )
        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_ccd )) deallocate( c_ccd )
        if (associated( c_den )) deallocate( c_den )
        if (associated( c_xcp )) deallocate( c_xcp )
        if (associated( c_sum )) deallocate( c_sum )
        nullify( cdff )
        nullify( gpff )

        call glean(thy(lat))
        call glean(thy(ats))
        call glean(thy(lay))

        call glean(thy(ad))
        call glean(thy(den))
        call glean(thy(xcp))
        call glean(thy(ccd))

        if (error("Exit atomic_density_ncp_mod::atomic_forces")) continue

      end subroutine

      subroutine atomic_pressure_ad(ad,den,xcp,p)
!doc$ subroutine atomic_pressure(ad,den,xcp,p)
        type(atomic_density_ncp_obj) :: ad
        type(grid_obj) :: den, xcp
        real(double), intent(out) :: p
!       requires: Input-grid scopes be SGROUP.
!       effects: Returns atomic contributions to the pressure.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, it
        real(double) :: pt, spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:), allocatable :: q
        real(double), dimension(:,:,:), pointer :: g2, gx, gy, gz
        real(double), dimension(:,:,:), pointer :: scdff, sgpff
        complex(double) :: igr
        complex(double), dimension(:,:,:), pointer :: c_den, c_sum, c_xcp
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(layout_obj) :: lay

        call my(ad)
        call my(den)
        call my(xcp)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)
        call my(x_atoms(x_crystal(ad%o%ao)),ats)
        call my(x_layout(den),lay)

        nullify( g2, gx, gy, gz, c_den, c_sum, c_xcp )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)
        call alloc(g2,lay,D_TYPE,SGROUP)
        g2 = gx**2 + gy**2 + gz**2

        call take(c_den,den,CDF_KIND)
        if (core_density(ad%o%ao)) then
          call take(c_xcp,xcp,CDF_KIND)
        end if

        call alloc(c_sum,lay,D_TYPE,SGROUP)

        pt = 0.0_double

        do it = 1,x_n_types(ad%o%ao)

          call stress_grid_potential_ff(ad%o%ao,it,sgpff)
          c_sum = sgpff*c_den

          call stress_core_density_ff(ad%o%ao,it,scdff)
          if (associated( scdff )) then
            c_sum = c_sum + spin_factor*scdff*c_xcp
          end if

          do ia = 1,x_n_atoms(ad%o%ao)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
            pos = lat2r(lat,x_position(ats,ia))
            do i3 = 1,size(gx,3)
            do i2 = 1,size(gx,2)
            do i1 = 1,size(gx,1)
              igr = (0.0_double,1.0_double)*(gx(i1,i2,i3)*pos(1) + gy(i1,i2,i3)*pos(2) + gz(i1,i2,i3)*pos(3))
              pt = pt - g2(i1,i2,i3)*real(c_sum(i1,i2,i3)*exp(igr),double)
            end do
            end do
            end do
          end do

        end do

        pt = pt/3.0_double
        call allreduce(CONFIG,MPI_SUM,pt,p)

        call put(c_den,den,CDF_KIND)
        if (core_density(ad%o%ao)) then
          call put(c_xcp,xcp,CDF_KIND)
        end if

        allocate( q(x_n_atoms(ad%o%ao)) )
        do ia = 1,size(q)
          q(ia) = x_atom_valence(ad%o%ao,ia)
        end do
        call ewald_pressure(x_crystal(ad%o%ao),q,pt) ; if (error()) goto 100
        p = p + pt

100     if (associated( g2 )) deallocate( g2 )
        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )
        if (associated( c_sum )) deallocate( c_sum )
        if (associated( c_xcp )) deallocate( c_xcp )
        if (allocated( q )) deallocate( q )

        nullify( scdff )
        nullify( sgpff )

        call glean(thy(lat))
        call glean(thy(ats))
        call glean(thy(lay))

        call glean(thy(ad))
        call glean(thy(den))
        call glean(thy(xcp))

        if (error("Exit atomic_density_ncp_mod::atomic_pressure_ad")) continue

      end subroutine

      subroutine atomic_stress_tensor_ad(ad,den,xcp,s)
!doc$ subroutine atomic_stress_tensor(ad,den,xcp,s)
        type(atomic_density_ncp_obj) :: ad
        type(grid_obj) :: den, xcp
        real(double), dimension(:,:), intent(out) :: s
!       requires: s be dimension(3,3).
!                 Input-grid scopes be SGROUP.
!       effects: Returns atomic contributions to the stress tensor.
!       errors: Passes errors.

!cod$
        integer :: i1, i2, i3, ia, it
        real(double) :: r0, spin_factor
        real(double), dimension(3) :: pos
        real(double), dimension(:), allocatable :: q
        real(double), dimension(:,:), allocatable :: st
        real(double), dimension(:,:,:), pointer :: gx, gy, gz
        real(double), dimension(:,:,:), pointer :: scdff, sgpff
        complex(double) :: igr
        complex(double), dimension(:,:,:), pointer :: c_den, c_sum, c_xcp
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(layout_obj) :: lay

        call my(ad)
        call my(den)
        call my(xcp)

        call my(x_lattice(x_crystal(ad%o%ao)),lat)
        call my(x_atoms(x_crystal(ad%o%ao)),ats)
        call my(x_layout(den),lay)

        nullify( gx, gy, gz, c_den, c_sum, c_xcp )

        spin_factor = 1.0_double/real(mpi_nsgroups(),double)

        call fmesh(gx,gy,gz,lay,D_TYPE,SGROUP)

        call take(c_den,den,CDF_KIND)
        if (core_density(ad%o%ao)) then
          call take(c_xcp,xcp,CDF_KIND)
        end if

        call alloc(c_sum,lay,D_TYPE,SGROUP)

        allocate( st(3,3) )
        st = 0.0_double

        do it = 1,x_n_types(ad%o%ao)

          call stress_grid_potential_ff(ad%o%ao,it,sgpff)
          c_sum = sgpff*c_den

          call stress_core_density_ff(ad%o%ao,it,scdff)
          if (associated( scdff )) then
            c_sum = c_sum + spin_factor*scdff*c_xcp
          end if

          do ia = 1,x_n_atoms(ad%o%ao)
            if (x_atom_type(ad%o%ao,ia) /= it) cycle
            pos = lat2r(lat,x_position(ats,ia))
            do i3 = 1,size(gx,3)
            do i2 = 1,size(gx,2)
            do i1 = 1,size(gx,1)
              igr = (0.0_double,1.0_double)*(gx(i1,i2,i3)*pos(1) + gy(i1,i2,i3)*pos(2) + gz(i1,i2,i3)*pos(3))
              r0 = real(c_sum(i1,i2,i3)*exp(igr),double)
              st(1,1) = st(1,1) + r0*gx(i1,i2,i3)*gx(i1,i2,i3)
              st(1,2) = st(1,2) + r0*gx(i1,i2,i3)*gy(i1,i2,i3)
              st(1,3) = st(1,3) + r0*gx(i1,i2,i3)*gz(i1,i2,i3)
              st(2,2) = st(2,2) + r0*gy(i1,i2,i3)*gy(i1,i2,i3)
              st(2,3) = st(2,3) + r0*gy(i1,i2,i3)*gz(i1,i2,i3)
              st(3,3) = st(3,3) + r0*gz(i1,i2,i3)*gz(i1,i2,i3)
            end do
            end do
            end do
          end do

        end do

        st(2,1) = st(1,2)
        st(3,1) = st(1,3)
        st(3,2) = st(2,3)
        call allreduce(CONFIG,MPI_SUM,st,s)

        call put(c_den,den,CDF_KIND)
        if (core_density(ad%o%ao)) then
          call put(c_xcp,xcp,CDF_KIND)
        end if

        allocate( q(x_n_atoms(ad%o%ao)) )
        do ia = 1,size(q)
          q(ia) = x_atom_valence(ad%o%ao,ia)
        end do
        call ewald_stress_tensor(x_crystal(ad%o%ao),q,st) ; if (error()) goto 100
        s = s + st

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( c_den )) deallocate( c_den )
        if (associated( c_sum )) deallocate( c_sum )
        if (associated( c_xcp )) deallocate( c_xcp )
        if (allocated( st)) deallocate( st )
        if (allocated( q )) deallocate( q )
        nullify( scdff )
        nullify( sgpff )

        call glean(thy(lat))
        call glean(thy(ats))
        call glean(thy(lay))

        call glean(thy(ad))
        call glean(thy(den))
        call glean(thy(xcp))

        if (error("Exit atomic_density_ncp_mod::atomic_stress_tensor_ad")) continue

      end subroutine

      subroutine write_restart_ad(ad,nrestf)
!doc$ subroutine write_restart(ad,nrestf)
        type(atomic_density_ncp_obj) :: ad
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ad restart information to nrestf.

!cod$
        integer(long) :: dsize, iosl, ndata

        call my(ad)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the ATOMIC_DENSITY block
          call startblock(nrestf,"ATOMIC_DENSITY")

          ! write the NCP tag
          call writetag(nrestf,"NCP")

          ! write the charge state
          call writetag(nrestf,"CHARGE_STATE")
          dsize = sizeof_double ; ndata = 1
          call writef(ad%o%charge_state,dsize,ndata,x_tagfd(nrestf),iosl)

          ! end the ATOMIC_DENSITY block
          call endblock(nrestf)

        end if

        call glean(thy(ad))
        call glean(thy(nrestf))

        if (error("Exit atomic_density_ncp_mod::write_restart_ad")) continue

      end subroutine

! private routines

      subroutine own_i(ad)
        type(atomic_density_ncp_obj) :: ad
        type(atomic_density_ncp_obj) :: adt
        integer :: i
        if (ad%ref < ad%o%ref) then
          allocate( adt%o )
          adt%o%ref = 0
          adt%o%g = ad%o%g
          adt%o%charge_state = ad%o%charge_state
          adt%o%energy = ad%o%energy
          if (associated( ad%o%wij )) then
             allocate( adt%o%wij(size(ad%o%wij)) )
             do i = 1,size(ad%o%wij)
               if (associated( ad%o%wij(i)%mat )) then
                  allocate( adt%o%wij(i)%mat(size(ad%o%wij(i)%mat,1),size(ad%o%wij(i)%mat,2)) ) ; adt%o%wij(i)%mat = ad%o%wij(i)%mat
               else
                  nullify( adt%o%wij(i)%mat )
               end if
             end do
             else
                nullify( adt%o%wij )
          end if
          call my(ad%o%ao,adt%o%ao)
          ad%o%ref = ad%o%ref - ad%ref
          ad%o => adt%o
          ad%o%ref = ad%o%ref + ad%ref
        end if
      end subroutine

      subroutine form_energy_i(adr)
        type(atomic_density_ncp_rep) :: adr

        integer :: ip, na, ia
        real(double) :: energy_sg

        energy_sg = 0.0_double
        na = size(adr%wij)
        do ia = 1,na
          energy_sg = energy_sg + x_atom_energy(adr%ao, adr%wij(ia)%mat, ia) 
!          energy_sg = energy_sg + x_projector_kbf(adr%ao,ip)*adr%wii(ip)
        end do
        call xcomm_allreduce(XSGROUP,MPI_SUM,energy_sg,adr%energy)

        if (error("Exit atomic_density_ncp_mod::form_energy_i")) continue

      end subroutine

      end module

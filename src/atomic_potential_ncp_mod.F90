! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module atomic_potential_ncp_mod
!doc$ module atomic_potential_ncp_mod

      use kind_mod
      use arg_mod
      use path_mod
      use math_mod
      use mpi_mod
      use error_mod
      use diary_mod
      use io_mod
      use tagio_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use lattice_mod
      use atoms_mod
      use crystal_mod
      use ncp_data_mod
      use atomic_operators_ncp_mod
      use atomic_density_ncp_mod

!     One datatype is available here: type(atomic_potential_ncp_obj).

!cod$
      implicit none
      private

      type :: atomic_potential_ncp_rep
        integer :: ref
        type(ghost) :: g
        type(atomic_operators_ncp_obj) :: ao
      end type

      type, public :: atomic_potential_ncp_obj
        private
        integer :: ref
        type(atomic_potential_ncp_rep), pointer :: o
      end type

!doc$
      public :: atomic_potential_ncp
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
      public :: write_restart

!cod$
      interface atomic_potential_ncp
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
      interface write_restart
        module procedure write_restart_ap
      end interface

      contains

! public routines

      function constructor_ap(ad,restf) result(ap)
!doc$ function atomic_potential_ncp(ad,restf) result(ap)
        type(atomic_density_ncp_obj) :: ad
        type(tagio_obj), optional :: restf
        type(atomic_potential_ncp_obj) :: ap
!       effects: Creates a new ap.
!       errors: Passes errors.

!cod$
        character(1) :: tios

        call my(ad)
        if (present(restf)) call my(restf)

        ap%ref = 0
        allocate( ap%o )
        ap%o%ref = 0
        ap%o%g = x_ghost()

        call my(x_atomic_operators(ad),ap%o%ao)

        if (present(restf)) then

          ! open the ATOMIC_POTENTIAL block
          if (i_access(restf)) tios = findfirsttag(restf,"ATOMIC_POTENTIAL")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ATOMIC_POTENTIAL block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)

          ! find the NCP tag
          if (i_access(restf)) tios = findfirsttag(restf,"NCP")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NCP tag was not found")) goto 100

          ! close the ATOMIC_POTENTIAL block
100       if (i_access(restf)) call closeblock(restf)
          if (error()) goto 200

        end if

200     call glean(thy(ad))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit atomic_potential_ncp_mod::constructor_ap")) continue

      end function

      subroutine update_ap(ap,ad)
!doc$ subroutine update(ap,ad)
        type(atomic_potential_ncp_obj) :: ap
        type(atomic_density_ncp_obj) :: ad
!       effects: Updates ap.

!cod$
        logical :: ao_change

        call my(ap)
        call my(ad)

        ao_change = ( x_ghost(ap%o%ao) /= x_ghost(x_atomic_operators(ad)) )
        if (ao_change) then
          call own_i(ap)
          ap%o%g = x_ghost()
          ap%o%ao = x_atomic_operators(ad)
        end if

        call glean(thy(ap))
        call glean(thy(ad))

        if (error("Exit atomic_potential_ncp_mod::update_ap")) continue

      end subroutine

      subroutine my_ap(ap)
!doc$ subroutine my(ap)
        type(atomic_potential_ncp_obj) :: ap

!cod$
        ap%ref = ap%ref + 1
        ap%o%ref = ap%o%ref + 1
      end subroutine

      subroutine my_new_ap(api,ap)
!doc$ subroutine my(api,ap)
        type(atomic_potential_ncp_obj) :: api
        type(atomic_potential_ncp_obj) :: ap

!cod$
        ap%ref = 1
        ap%o => api%o
        ap%o%ref = ap%o%ref + 1
      end subroutine

      function thy_ap(ap) result(apo)
!doc$ function thy(ap) result(apo)
        type(atomic_potential_ncp_obj) :: ap, apo

!cod$
        ap%ref = ap%ref - 1
        ap%o%ref = ap%o%ref - 1
        apo%ref = ap%ref
        apo%o => ap%o
      end function

      subroutine glean_ap(ap)
!doc$ subroutine glean(ap)
        type(atomic_potential_ncp_obj) :: ap

!cod$
        if (ap%o%ref < 1) then
          call glean(thy(ap%o%ao))
          deallocate( ap%o )
        end if
      end subroutine

      subroutine bequeath_ap(ap)
!doc$ subroutine bequeath(ap)
        type(atomic_potential_ncp_obj) :: ap

!cod$
        continue
      end subroutine
    
      subroutine assign_ap(ap,ap2)
!doc$ subroutine assignment(=)(ap,ap2)
        type(atomic_potential_ncp_obj), intent(inout) :: ap
        type(atomic_potential_ncp_obj), intent(in) :: ap2

!cod$
        type(atomic_potential_ncp_obj) :: apt
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
        type(atomic_potential_ncp_obj) :: ap
        integer, dimension(2) :: r
!       effects: Returns ap%ref and ap%o%ref.

!cod$
        r(1) = ap%ref
        r(2) = ap%o%ref
        call glean(ap)
      end function

      function ap_ghost(ap) result(g)
!doc$ function x_ghost(ap) result(g)
        type(atomic_potential_ncp_obj) :: ap
        type(ghost) :: g
!       effects: Returns the ghost of ap.

!cod$
        call my(ap)
        g = ap%o%g
        call glean(thy(ap))
      end function

      function ap_atomic_operators(ap) result(ao)
!doc$ function x_atomic_operators(ap) result(ao)
        type(atomic_potential_ncp_obj) :: ap
        type(atomic_operators_ncp_obj) :: ao
!       effects: Returns the ao of ap.

!cod$
        call my(ap)
        call my(ap%o%ao,ao)
        call bequeath(thy(ao))
        call glean(thy(ap))
      end function

      subroutine atomic_hamiltonian_1d_ap(ap,pdots)
!doc$ subroutine atomic_hamiltonian(ap,pdots)
        type(atomic_potential_ncp_obj) :: ap
        complex(double), dimension(:), intent(inout) :: pdots
!       requires: size(pdots) = x_n_projectors(ap%o%ao).
!       modifies: pdots
!       effects: Multiplies pdots by the Kleinman-Bylander factors.

!cod$
        integer :: ip, i
        real(double) :: sf
        logical :: found

        call my(ap)
          call apply_projector_kbf(ap%o%ao, pdots)
          !call arg("pdots_scale_factor",sf,found)
          !if(.not.found) sf = 1.0_double
          !if (i_access(diaryfile())) write (x_unit(diaryfile()), '(/,"scaling factor is=", f15.10)') sf
          !pdots = sf*pdots
          call glean(thy(ap))

      end subroutine

      subroutine atomic_hamiltonian_2d_ap(ap,pdots)
!doc$ subroutine atomic_hamiltonian(ap,pdots)
        type(atomic_potential_ncp_obj) :: ap
        complex(double), dimension(:,:), intent(inout) :: pdots
!       requires: size(pdots,1) = x_n_projectors(ap%o%ao).
!       modifies: pdots
!       effects: Multiplies pdots by the Kleinman-Bylander factors.

!cod$
        integer :: ib,i
        complex(double), dimension(1:size(pdots,1)) :: pdots_p
        real(double) :: sf
        logical :: found

        call my(ap)
        do ib = 1,size(pdots,2)
          pdots_p = pdots(:,ib)
          call apply_projector_kbf(ap%o%ao, pdots_p)
          pdots(:,ib) = pdots_p
        end do
        !call arg("pdots_scale_factor",sf,found)
          !if(.not.found) sf = 1.0_double
          !if (i_access(diaryfile())) write (x_unit(diaryfile()), '(/,"scaling factor is=", f15.10)') sf
        !pdots = sf*pdots
        call glean(thy(ap))

      end subroutine

      subroutine atomic_overlap_1d_ap(ap,pdots)
!doc$ subroutine atomic_overlap(ap,pdots)
        type(atomic_potential_ncp_obj) :: ap
        complex(double), dimension(:), intent(inout) :: pdots
!       modifies: pdots
!       effects: Sets pdots to zero.

!cod$
        pdots = (0.0_double,0.0_double)
        call glean(ap)
      end subroutine

      subroutine atomic_overlap_2d_ap(ap,pdots)
!doc$ subroutine atomic_overlap(ap,pdots)
        type(atomic_potential_ncp_obj) :: ap
        complex(double), dimension(:,:), intent(inout) :: pdots
!       modifies: pdots
!       effects: Sets pdots to zero.

!cod$
        pdots = (0.0_double,0.0_double)
        call glean(ap)
      end subroutine

      subroutine write_restart_ap(ap,nrestf)
!doc$ subroutine write_restart(ap,nrestf)
        type(atomic_potential_ncp_obj) :: ap
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ap restart information to nrestf.

!cod$
        call my(ap)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the ATOMIC_POTENTIAL block
          call startblock(nrestf,"ATOMIC_POTENTIAL")

          ! write the NCP tag
          call writetag(nrestf,"NCP")

          ! end the ATOMIC_POTENTIAL block
          call endblock(nrestf)

        end if

        call glean(thy(ap))
        call glean(thy(nrestf))

        if (error("Exit atomic_potential_ncp_mod::write_restart_ap")) continue

      end subroutine

! private routines

      subroutine own_i(ap)
        type(atomic_potential_ncp_obj) :: ap
        type(atomic_potential_ncp_obj) :: apt
        if (ap%ref < ap%o%ref) then
          allocate( apt%o )
          apt%o%ref = 0
          apt%o%g = ap%o%g
          call my(ap%o%ao,apt%o%ao)
          ap%o%ref = ap%o%ref - ap%ref
          ap%o => apt%o
          ap%o%ref = ap%o%ref + ap%ref
        end if
      end subroutine

      end module

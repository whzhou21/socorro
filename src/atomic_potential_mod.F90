! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module atomic_potential_mod
!doc$ module atomic_potential_mod

      use kind_mod
      use error_mod
      use tagio_mod
      use ghost_mod
      use grid_mod
      use atomic_operators_mod
      use atomic_density_mod
      use atomic_potential_ncp_mod
      use atomic_potential_paw_mod

!     One datatype is defined here: type(atomic_potential_obj).

!cod$
      implicit none
      private

      type, public :: atomic_potential_obj
        private
        integer :: type
        type(atomic_potential_ncp_obj) :: ap_ncp
        type(atomic_potential_paw_obj) :: ap_paw
      end type

!doc$
      public :: atomic_potential
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: trivial
      public :: overlap_is_identity
      public :: atomic_hamiltonian
      public :: atomic_overlap
      public :: extract_potential
      public :: insert_potential
      public :: write_restart

!cod$
      interface atomic_potential
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
      interface trivial
        module procedure trivial_ap
      end interface
      interface overlap_is_identity
        module procedure overlap_is_identity_ap
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
!doc$ function atomic_potential(ad,hap,restf) result(ap)
        type(atomic_density_obj)  :: ad
        type(grid_obj), optional :: hap
        type(tagio_obj), optional :: restf
        type(atomic_potential_obj) :: ap
!       requires: Either hap or restf be absent.
!       effects: Constructs a new ap.

!cod$
        ap%type = x_type(ad)

        select case (ap%type)
        case (NCP)
          if (present(restf)) then
            call my(atomic_potential_ncp(x_atomic_density_ncp(ad),restf),ap%ap_ncp) ; if (error()) goto 100
          else
            call my(atomic_potential_ncp(x_atomic_density_ncp(ad)),ap%ap_ncp) ; if (error()) goto 100
            call glean(hap)
          end if
          call bequeath(thy(ap%ap_ncp))
        case (PAW)
          if (present(restf)) then
            call my(atomic_potential_paw(x_atomic_density_paw(ad),restf=restf),ap%ap_paw) ; if (error()) goto 100
          else
            call my(atomic_potential_paw(x_atomic_density_paw(ad),hap),ap%ap_paw) ; if (error()) goto 100
          end if
          call bequeath(thy(ap%ap_paw))
        end select

100     if (error("Exit atomic_potential_mod::constructor_ap")) continue

      end function

      subroutine update_ap(ap,ad,hap)
!doc$ subroutine update(ap,ad,hap)
        type(atomic_potential_obj) :: ap
        type(atomic_density_obj) :: ad
        type(grid_obj) :: hap
!       requires: ap%type /= ad%type.
!       effects: Updates ap.

!cod$
        select case (ap%type)
        case (NCP)
          call update(ap%ap_ncp,x_atomic_density_ncp(ad))
          call glean(hap)
        case (PAW)
          call update(ap%ap_paw,x_atomic_density_paw(ad),hap)
        end select
        if (error("Exit atomic_potential_mod::update_ap")) continue
      end subroutine

      subroutine my_ap(ap)
!doc$ subroutine my(ap)
        type(atomic_potential_obj) :: ap

!cod$
        select case (ap%type)
        case (NCP)
          call my(ap%ap_ncp)
        case (PAW)
          call my(ap%ap_paw)
        end select
      end subroutine

      subroutine my_new_ap(api,ap)
!doc$ subroutine my(api,ap)
        type(atomic_potential_obj) :: api, ap

!cod$
        ap%type = api%type
        select case (ap%type)
        case (NCP)
          call my(api%ap_ncp,ap%ap_ncp)
        case (PAW)
          call my(api%ap_paw,ap%ap_paw)
        end select
      end subroutine

      function thy_ap(ap) result(apo)
!doc$ function thy(ap) result(apo)
        type(atomic_potential_obj) :: ap, apo

!cod$
        apo%type = ap%type
        select case (ap%type)
        case (NCP)
          call my(thy(ap%ap_ncp),apo%ap_ncp)
          call bequeath(thy(apo%ap_ncp))
        case (PAW)
          call my(thy(ap%ap_paw),apo%ap_paw)
          call bequeath(thy(apo%ap_paw))
        end select
      end function

      subroutine glean_ap(ap)
!doc$ subroutine glean(ap)
        type(atomic_potential_obj) :: ap

!cod$
        select case (ap%type)
        case (NCP)
          call glean(ap%ap_ncp)
        case (PAW)
          call glean(ap%ap_paw)
        end select
      end subroutine

      subroutine bequeath_ap(ap)
!doc$ subroutine bequeath(ap)
        type(atomic_potential_obj) :: ap

!cod$
        select case (ap%type)
        case (NCP)
          call bequeath(ap%ap_ncp)
        case (PAW)
          call bequeath(ap%ap_paw)
        end select
      end subroutine

      subroutine assign_ap(ap,ap2)
!doc$ subroutine assignment(=)(ap,ap2)
        type(atomic_potential_obj), intent(inout) :: ap
        type(atomic_potential_obj), intent(in) :: ap2
!       requires: ap and ap2 have the same type.

!cod$
        select case (ap%type)
        case (NCP)
          ap%ap_ncp = ap2%ap_ncp
        case (PAW)
          ap%ap_paw = ap2%ap_paw
        end select
      end subroutine

      function ap_ref(ap) result(r)
!doc$ function x_ref(ap) result(r)
        type(atomic_potential_obj) :: ap
        integer, dimension(2) :: r
!       effects: Returns ap%ref and ap%o%ref.

!cod$
        select case (ap%type)
        case (NCP)
          r = x_ref(ap%ap_ncp)
        case (PAW)
          r = x_ref(ap%ap_paw)
        end select
      end function

      function ap_ghost(ap) result(g)
!doc$ function x_ghost(ap) result(g)
        type(ghost) :: g
        type(atomic_potential_obj) :: ap
!       effects: Returns the ghost of ap.

!cod$
        select case (ap%type)
        case (NCP)
          g = x_ghost(ap%ap_ncp)
        case (PAW)
          g = x_ghost(ap%ap_paw)
        end select
      end function

      function trivial_ap(ap) result(t)
!doc$ function trivial(ap) result(t)
        type(atomic_potential_obj) :: ap
        logical :: t
!       effects: Returns the triviality of ap.

!cod$
        select case (ap%type)
        case (NCP)
          t = .true.
        case (PAW)
          t = .false.
        end select
        call glean(ap)
      end function

      function overlap_is_identity_ap(ap) result(t)
!doc$ function overlap_is_identity(ap) result(t)
        type(atomic_potential_obj) :: ap
        logical :: t
!       effects: Returns true if the overlap is the identity

!cod$
        select case (ap%type)
        case (NCP)
          t = .true.
        case (PAW)
          t = .false.
        end select
        call glean(ap)
      end function

      subroutine atomic_hamiltonian_1d_ap(ap,pdots)
!doc$ subroutine atomic_hamiltonian(ap,pdots)
        type(atomic_potential_obj) :: ap
        complex(double), dimension(:), intent(inout) :: pdots
!       modifies: pdots
!       effects: Multiplies pdots by the atomic hamiltonian matrix.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          call atomic_hamiltonian(ap%ap_ncp,pdots)
        case (PAW)
          call atomic_hamiltonian(ap%ap_paw,pdots)
        end select
      end subroutine

      subroutine atomic_hamiltonian_2d_ap(ap,pdots)
!doc$ subroutine atomic_hamiltonian(ap,pdots)
        type(atomic_potential_obj) :: ap
        complex(double), dimension(:,:), intent(inout) :: pdots
!       modifies: pdots
!       effects: Multiplies pdots by the atomic hamiltonian matrix.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          call atomic_hamiltonian(ap%ap_ncp,pdots)
        case (PAW)
          call atomic_hamiltonian(ap%ap_paw,pdots)
        end select
      end subroutine

      subroutine atomic_overlap_1d_ap(ap,pdots)
!doc$ subroutine atomic_overlap(ap,pdots)
        type(atomic_potential_obj) :: ap
        complex(double), dimension(:), intent(inout) :: pdots
!       modifies: pdots
!       effects: Multiplies pdots by the atomic overlap matrix.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          call atomic_overlap(ap%ap_ncp,pdots)
        case (PAW)
          call atomic_overlap(ap%ap_paw,pdots)
        end select
      end subroutine

      subroutine atomic_overlap_2d_ap(ap,pdots)
!doc$ subroutine atomic_overlap(ap,pdot)s
        type(atomic_potential_obj) :: ap
        complex(double), dimension(:,:), intent(inout) :: pdots
!       modifies: pdots
!       effects: Multiplies pdots by the atomic overlap matrix.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          call atomic_overlap(ap%ap_ncp,pdots)
        case (PAW)
          call atomic_overlap(ap%ap_paw,pdots)
        end select
      end subroutine

      subroutine extract_potential_ap(ap,c1d)
!doc$ subroutine extract_potential(ap,c1d)
        type(atomic_potential_obj) :: ap
        complex(double), dimension(:), pointer :: c1d
!       requires: c1d be nullified or associated.
!       modifies: c1d
!       effects: Returns ap%o%dij in c1d.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          if (associated( c1d )) deallocate( c1d ) ; allocate( c1d(0) )
          call glean(ap)
        case (PAW)
          call extract_potential(ap%ap_paw,c1d)
        end select
      end subroutine

      subroutine insert_potential_ap(c1d,ap)
!doc$ subroutine insert_potential(c1d,ap)
        complex(double), dimension(:), pointer :: c1d
        type(atomic_potential_obj) :: ap
!       modifies: ap
!       effects: Inserts c1d into ap%o%dij.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          call glean(ap)
        case (PAW)
          call insert_potential(c1d,ap%ap_paw)
        end select
      end subroutine

      subroutine write_restart_ap(ap,nrestf)
!doc$ subroutine write_restart(ap,nrestf)
        type(atomic_potential_obj) :: ap
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes ap restart information to nrestf.
!       errors: Passes errors.

!cod$
        select case (ap%type)
        case (NCP)
          call write_restart(ap%ap_ncp,nrestf)
        case (PAW)
          call write_restart(ap%ap_paw,nrestf)
        end select
      end subroutine

      end module

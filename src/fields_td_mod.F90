!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module fields_td_mod
!doc$ module fields_td_mod

!     One datatype is available here: type(fields_td_obj)

!     fields_td_mod manages field quantities in calculations with a time-dependent hamiltonian.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use arg_mod
      use diary_mod
      use math_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use lattice_mod
      use atoms_mod
      use crystal_mod
      use symmetry_mod
      use external_mod
      use electrons_td_mod
      use xc_type_mod
      use xc_mod
      use gen_density_mod
      use gen_potential_mod
      use atomic_operators_mod
      use atomic_density_mod
      use atomic_potential_mod
      use timing_mod

!cod$
      implicit none
      private 

      integer, parameter :: MOD_SCOPE = CONFIG

      integer, parameter :: NONE = 0
      integer, parameter :: UBC  = 1
      integer, parameter :: LMCC = 2

      type :: fields_td_rep
        integer :: ref
        type(ghost) :: g
        type(ghost) :: g_ao                                ! atomic_operators ghost
        integer :: compensation_method                     ! compensation method
        real(double) :: exc_energy                         ! exchange-correlation energy
        real(double) :: har_energy                         ! Hartree energy
        real(double) :: agp_energy                         ! atomic grid potential energy
        real(double) :: atomic_energy                      ! atomic energy
        real(double) :: charge_state                       ! charge state of supercell
        real(double) :: lmcc_width                         ! width of the Gaussian charge distribution in the LMCC method
        real(double) :: lmcc_energy                        ! energy contribution from the LMCC method (independent of gden)
        real(double) :: lmcc_fixed_energy                  ! fixed part of lmcc_energy
        real(double) :: ccp_energy                         ! energy contribution from the LMCC method (dependent on gden)
        real(double), dimension(3) :: lmcc_site            ! site of the Gaussian charge distribution in the LMCC method
        type(grid_obj) :: gden                             ! grid density
        type(grid_obj) :: ahd                              ! atomic hartree density
        type(grid_obj) :: axcd                             ! atomic exchange-correlation density
        type(grid_obj) :: agp                              ! atomic grid potential
        type(grid_obj) :: hap                              ! hartree potential
        type(grid_obj) :: xcp                              ! exchange-correlation potential
        type(grid_obj) :: scp                              ! self-consistent potential
        type(grid_obj) :: ccp                              ! counter-charge potential in the LMCC method
        type(grid_obj) :: total                            ! total potential
        type(atomic_density_obj) :: aden                   ! atomic density object
        type(atomic_potential_obj) :: apot                 ! atomic potential object
        type(xc_obj) :: xc                                 ! exchange-correlation object
        real(double), dimension(:,:,:), pointer :: ck      ! Coulomb kernel (conformable to grid data with CDF_KIND)
      end type

      type, public :: fields_td_obj
        private
        integer :: ref
        type(fields_td_rep), pointer :: o
      end type

!doc$
      public :: fields_td
      public :: update
      public :: my
      public :: thy
      public :: bequeath
      public :: glean
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_energy
      public :: x_potential
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: put_smooth_density
      public :: put_local_potential
      public :: put_hartree_potential
      public :: put_xc_potential
      public :: put_atomic_xc_density
      public :: put_atomic_hartree_density
      public :: put_atomic_grid_potential
      public :: max_energy
      public :: diary_energies
      public :: write_restart

!cod$
      interface fields_td
        module procedure constructor_fd
      end interface
      interface update
        module procedure update_fd
      end interface
      interface my
        module procedure my_fd, my_new_fd
      end interface
      interface thy
        module procedure thy_fd
      end interface
      interface glean
        module procedure glean_fd
      end interface
      interface bequeath
        module procedure bequeath_fd
      end interface
      interface assignment(=)
        module procedure assign_fd
      end interface
      interface x_ref
        module procedure fd_ref
      end interface
      interface x_ghost
        module procedure fd_ghost
      end interface
      interface x_energy
        module procedure fd_energy
      end interface
      interface x_potential
        module procedure fd_potential
      end interface
      interface forces
        module procedure forces_fd
      end interface
      interface pressure
        module procedure pressure_fd
      end interface
      interface stress_tensor
        module procedure stress_tensor_fd
      end interface
      interface put_smooth_density
        module procedure put_smooth_density_fd
      end interface
      interface put_local_potential
        module procedure put_local_potential_fd
      end interface
      interface put_hartree_potential
        module procedure put_hartree_potential_fd
      end interface
      interface put_xc_potential
        module procedure put_xc_potential_fd
      end interface
      interface put_atomic_xc_density
        module procedure put_atomic_xc_density_fd
      end interface
      interface put_atomic_hartree_density
        module procedure put_atomic_hartree_density_fd
      end interface
      interface put_atomic_grid_potential
        module procedure put_atomic_grid_potential_fd
      end interface
      interface max_energy
        module procedure max_energy_fd
      end interface
      interface diary_energies
        module procedure diary_energies_fd
      end interface
      interface write_restart
        module procedure write_restart_fd
      end interface

      contains

! public routines

      function constructor_fd(ext,el,restf) result(fd)
!doc$ function fields_td(ext,el) result(fd)
        type(external_obj)               :: ext
        type(electrons_td_obj),optional  :: el
        type(tagio_obj), optional        :: restf
        type(fields_td_obj)              :: fd
!       requires: Either el or restf must be present but not both. If el present, must be consistent with ext.
!       effects: Constructs a new fd.
!       errors: Passes errors.

!cod$
        logical :: found
        real(double) :: e_local, e_global
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(lattice_obj) :: lat
        type(layout_obj) :: lay
        type(atomic_operators_obj) :: ao
        type(grid_obj) :: gr1

        call start_timer("fields_td: constructor")

        call my(ext)

        if (present(el) .and. present(restf)) then
           if (error(.true.,"Error: either electrons OR restartfile must be passed but not both")) goto 400
        endif

        if (present(el)) call my(el)
        if (present(restf)) call my(restf)

        nullify( c1, c2 )

        fd%ref = 0
        allocate( fd%o )
        fd%o%ref = 0
        fd%o%g = x_ghost()

        call my(x_lattice(x_crystal(ext)),lat)
        call my(x_layout(ext),lay)
        call my(x_atomic_operators(ext),ao)
        call my(grid(lay,SGROUP),gr1)

        fd%o%g_ao = x_ghost(ao)

        call my(xc(x_xc_type(ext),lay,x_space_group(ext)),fd%o%xc)

        ! Coulomb kernel
        nullify( fd%o%ck )
        call form_coulomb_kernel_i(fd%o,lay)

        if (present(restf)) then
          ! atomic density, grid density, total potential and scp potential
          call restart_portal_i(fd%o,ext,restf) ; if (error()) goto 100
        else
          ! Just atomic density and grid density
          call my(x_atomic_density(x_density(el)),fd%o%aden)

          call my(x_grid_density(x_density(el)),fd%o%gden)
        end if

        ! atomic energy
        call atomic_energy(fd%o%aden,fd%o%atomic_energy) ; if (error()) goto 100

        ! atomic exchange-correlation density
        call my(atomic_xc_density(ao,lay),fd%o%axcd) ; if (error()) goto 100

        ! exchange-correlation energy and potential
        gr1 = fd%o%gden
        call saxpby(1.0_double,gr1,1.0_double,fd%o%axcd)
        call my(xc_energy_and_potential(fd%o%xc,gr1,fd%o%gden,fd%o%exc_energy),fd%o%xcp) ; if (error()) goto 100

        ! atomic hartree density
        call my(atomic_hartree_density(fd%o%aden,lay),fd%o%ahd) ; if (error()) goto 100

        ! hartree potential and energy
        gr1 = fd%o%gden ; if (error ("Error - gr1 = fd%o%gden for hartree potential and energy")) goto 100
        call saxpby(1.0_double,gr1,1.0_double,fd%o%ahd)
        call take(c1,gr1,CDF_KIND)
        call alloc(c2,lay,D_TYPE,MOD_SCOPE)
        c2 = eight_pi*c1*fd%o%ck
        e_local = 0.5_double*real(sum(conjg(c1)*c2),double)
        call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
        fd%o%har_energy = e_global*x_cell_volume(lat)
        call put(c1,gr1,CDF_KIND)
        call my(grid(lay,MOD_SCOPE),fd%o%hap)
        call put(c2,fd%o%hap,CDF_KIND)

        ! atomic grid potential
        call my(atomic_grid_potential(ao,lay),fd%o%agp) ; if (error()) goto 100

        ! atomic grid potential energy
        call take(c1,fd%o%gden,CDF_KIND)
        call take(c2,fd%o%agp,CDF_KIND)
        e_local = real(sum(c1*conjg(c2)),double)
        call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
        fd%o%agp_energy = e_global*x_cell_volume(lat)
        call put(c1,fd%o%gden,CDF_KIND)
        call put(c2,fd%o%agp,CDF_KIND)

        ! charge compensation
        call form_compensation_i(fd%o,ext) ; if (error()) goto 100

        ! ccp energy
        if (x_type(fd%o%ccp) /= EMPTY_KIND) then
          call take(c1,fd%o%gden,CDF_KIND)
          call take(c2,fd%o%ccp,CDF_KIND)
          e_local = real(sum(c1*conjg(c2)),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          fd%o%ccp_energy = e_global*x_cell_volume(lat)
          call put(c1,fd%o%gden,CDF_KIND)
          call put(c2,fd%o%ccp,CDF_KIND)
        end if

        ! atomic potential
        call my(atomic_potential(fd%o%aden,fd%o%hap),fd%o%apot)

        ! self-consistent potential
        if (.not.present(restf)) then  !** scp is already set if restarting from file
           call my(fd%o%hap,fd%o%scp)
           call saxpby(1.0_double,fd%o%scp,1.0_double,fd%o%xcp)
        endif

        ! total potential 
        if (.not.present(restf)) then  !** total is already set if restarting from file
           call my(fd%o%scp,fd%o%total)
           call saxpby(1.0_double,fd%o%total,1.0_double,fd%o%agp)
           call saxpby(1.0_double,fd%o%total,1.0_double,fd%o%ccp)
           call transform(fd%o%total,RD_KIND) ; if (error()) goto 100
        end if

        call diary_construction_i(fd%o,ao)

100     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(gr1))

        call glean(thy(lay))
        call glean(thy(lat))
        call glean(thy(ao))

        call glean(thy(ext))
        if (present(el)) call glean(thy(el))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit fields_td_mod::constructor_fd")) continue

400     if (.not.error()) call stop_timer("fields_td: constructor")

      end function

      subroutine update_fd(fd,ext,el)
!doc$ subroutine update(fd,ext,el)
        type(fields_td_obj) :: fd
        type(external_obj) :: ext
        type(electrons_td_obj) :: el
!       requires: el be consistent with ext.
!       modifies: fd
!       effects: Updates fd.
!       errors: Passes errors.

!cod$
        logical :: ao_change
        real(double) :: e_local, e_global
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(lattice_obj) :: lat
        type(layout_obj) :: lay
        type(gen_density_obj) :: gd
        type(atomic_operators_obj) :: ao
        type(grid_obj) :: gr1

        call start_timer("fields_td: update")

        call my(fd)
        call my(ext)
        call my(el)

        nullify( c1, c2 )

        call my(x_lattice(x_crystal(ext)),lat)
        call my(x_layout(ext),lay)
        call my(x_atomic_operators(ext),ao)
        call my(x_density(el),gd)
        call my(grid(lay,SGROUP),gr1)

        call own_i(fd)
        fd%o%g = x_ghost()

        ao_change = ( x_ghost(ao) /= fd%o%g_ao )

        if (ao_change) then

          fd%o%g_ao = x_ghost(ao)

          ! atomic density
          fd%o%aden = x_atomic_density(gd)

          ! grid density
          fd%o%gden = x_grid_density(gd)

          call update(fd%o%aden,ao)

          ! atomic energy
          call atomic_energy(fd%o%aden,fd%o%atomic_energy) ; if (error()) goto 100

          ! atomic exchange-correlation density
          fd%o%axcd = atomic_xc_density(ao,lay) ; if (error()) goto 100

          ! exchange-correlation energy and potential
          gr1 = fd%o%gden
          call saxpby(1.0_double,gr1,1.0_double,fd%o%axcd)
          fd%o%xcp = xc_energy_and_potential(fd%o%xc,gr1,fd%o%gden,fd%o%exc_energy) ; if (error()) goto 100

          ! atomic hartree density
          fd%o%ahd = atomic_hartree_density(fd%o%aden,lay) ; if (error()) goto 100

          ! hartree potential and energy
          gr1 = fd%o%gden
          call saxpby(1.0_double,gr1,1.0_double,fd%o%ahd)
          call take(c1,gr1,CDF_KIND)
          call alloc(c2,lay,D_TYPE,MOD_SCOPE)
          c2 = eight_pi*c1*fd%o%ck
          e_local = 0.5_double*real(sum(conjg(c1)*c2),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          fd%o%har_energy = e_global*x_cell_volume(lat)
          call put(c1,gr1,CDF_KIND)
          call put(c2,fd%o%hap,CDF_KIND)

          ! atomic grid potential
          fd%o%agp = atomic_grid_potential(ao,lay) ; if (error()) goto 100

          ! atomic grid potential energy
          call take(c1,fd%o%gden,CDF_KIND)
          call take(c2,fd%o%agp,CDF_KIND)
          e_local = real(sum(c1*conjg(c2)),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          fd%o%agp_energy = e_global*x_cell_volume(lat)
          call put(c1,fd%o%gden,CDF_KIND)
          call put(c2,fd%o%agp,CDF_KIND)

          ! charge compensation
          call update_compensation_i(ext,fd%o) ; if (error()) goto 100

          ! ccp energy
          if (x_type(fd%o%ccp) /= EMPTY_KIND) then
            call take(c1,fd%o%gden,CDF_KIND)
            call take(c2,fd%o%ccp,CDF_KIND)
            e_local = real(sum(c1*conjg(c2)),double)
            call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
            fd%o%ccp_energy = e_global*x_cell_volume(lat)
            call put(c1,fd%o%gden,CDF_KIND)
            call put(c2,fd%o%ccp,CDF_KIND)
          end if

          ! atomic potential
          call update(fd%o%apot,fd%o%aden,fd%o%hap)

          ! self-consistent potential
          fd%o%scp = fd%o%hap
          call saxpby(1.0_double,fd%o%scp,1.0_double,fd%o%xcp)

          ! total potential
          fd%o%total = fd%o%scp
          call saxpby(1.0_double,fd%o%total,1.0_double,fd%o%agp)
          call saxpby(1.0_double,fd%o%total,1.0_double,fd%o%ccp)
          call transform(fd%o%total,RD_KIND)

        else

          ! atomic density
          fd%o%aden = x_atomic_density(gd)

          ! grid density
          fd%o%gden = x_grid_density(gd)

          ! atomic energy
          call atomic_energy(fd%o%aden,fd%o%atomic_energy) ; if (error()) goto 100

          ! exchange-correlation energy and potential
          gr1 = fd%o%gden
          !call my(fd%o%gden,gr1)
          call saxpby(1.0_double,gr1,1.0_double,fd%o%axcd)
          fd%o%xcp = xc_energy_and_potential(fd%o%xc,gr1,fd%o%gden,fd%o%exc_energy) ; if (error()) goto 100

          ! atomic hartree density
          fd%o%ahd = atomic_hartree_density(fd%o%aden,lay) ; if (error()) goto 100

          ! hartree potential and energy
          gr1 = fd%o%gden
          call saxpby(1.0_double,gr1,1.0_double,fd%o%ahd)
          call take(c1,gr1,CDF_KIND)
          call alloc(c2,lay,D_TYPE,MOD_SCOPE)
          c2 = eight_pi*c1*fd%o%ck
          e_local = 0.5_double*real(sum(conjg(c1)*c2),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          fd%o%har_energy = e_global*x_cell_volume(lat)
          call put(c1,gr1,CDF_KIND)
          call put(c2,fd%o%hap,CDF_KIND)

          ! atomic grid potential energy
          call take(c1,fd%o%gden,CDF_KIND)
          call take(c2,fd%o%agp,CDF_KIND)
          e_local = real(sum(c1*conjg(c2)),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          fd%o%agp_energy = e_global*x_cell_volume(lat)
          call put(c1,fd%o%gden,CDF_KIND)
          call put(c2,fd%o%agp,CDF_KIND)

          ! ccp energy
          if (x_type(fd%o%ccp) /= EMPTY_KIND) then
            call take(c1,fd%o%gden,CDF_KIND)
            call take(c2,fd%o%ccp,CDF_KIND)
            e_local = real(sum(c1*conjg(c2)),double)
            call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
            fd%o%ccp_energy = e_global*x_cell_volume(lat)
            call put(c1,fd%o%gden,CDF_KIND)
            call put(c2,fd%o%ccp,CDF_KIND)
          end if

          ! atomic potential
          call update(fd%o%apot,fd%o%aden,fd%o%hap)

          ! self-consistent potential
          fd%o%scp = fd%o%hap
          call saxpby(1.0_double,fd%o%scp,1.0_double,fd%o%xcp)

          ! total potential
          fd%o%total = fd%o%scp
          call saxpby(1.0_double,fd%o%total,1.0_double,fd%o%agp)
          call saxpby(1.0_double,fd%o%total,1.0_double,fd%o%ccp)
          call transform(fd%o%total,RD_KIND)

        end if

100     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(lat))
        call glean(thy(lay))
        call glean(thy(gd))
        call glean(thy(ao))
        call glean(thy(gr1))

        call glean(thy(fd))
        call glean(thy(ext))
        call glean(thy(el))

        if (error("Exit fields_td_mod::update_fd")) continue

        if (.not.error()) call stop_timer("fields_td: update")

      end subroutine

      subroutine my_fd(fd)
!doc$ subroutine my(fd)
        type(fields_td_obj) :: fd

!cod$
        fd%ref = fd%ref + 1
        fd%o%ref = fd%o%ref + 1
      end subroutine

      subroutine my_new_fd(fdi,fd)
!doc$ subroutine my(fdi,fd)
        type(fields_td_obj) :: fdi, fd

!cod$
        fd%ref = 1
        fd%o => fdi%o
        fd%o%ref = fd%o%ref + 1
      end subroutine

      function thy_fd(fd) result(fdo)
!doc$ function thy(fd) result(fdo)
        type(fields_td_obj) :: fd, fdo

!cod$
        fd%ref = fd%ref - 1
        fd%o%ref = fd%o%ref - 1
        fdo%ref = fd%ref
        fdo%o => fd%o
      end function

      subroutine glean_fd(fd)
!doc subroutine glean(fd)
        type(fields_td_obj) :: fd

!cod$
        if (fd%o%ref < 1) then
          call glean(thy(fd%o%gden))
          call glean(thy(fd%o%ahd))
          call glean(thy(fd%o%axcd))
          call glean(thy(fd%o%agp))
          call glean(thy(fd%o%hap))
          call glean(thy(fd%o%xcp))
          call glean(thy(fd%o%scp))
          call glean(thy(fd%o%ccp))
          call glean(thy(fd%o%total))
          call glean(thy(fd%o%aden))
          call glean(thy(fd%o%apot))
          call glean(thy(fd%o%xc))
          if (associated( fd%o%ck )) deallocate( fd%o%ck )
          deallocate( fd%o )
        end if
      end subroutine

      subroutine bequeath_fd(fd)
!doc$ subroutine bequeath(fd)
        type(fields_td_obj) :: fd

!cod$
        continue
      end subroutine

      subroutine assign_fd(fd,fd2)
!doc$ subroutine assign(fd,fd2)
        type(fields_td_obj), intent(inout) :: fd
        type(fields_td_obj), intent(in) :: fd2

!cod$
        type(fields_td_obj) :: fdt
        call my(fd2)
        fdt%o => fd%o
        fd%o%ref = fd%o%ref - fd%ref
        fd%o => fd2%o
        fd%o%ref = fd%o%ref + fd%ref
        call glean(fdt)
        call glean(thy(fd2))
      end subroutine

      function fd_ref(fd) result(r)
!doc$ function x_ref(fd) result(r)
        type(fields_td_obj) :: fd
        integer, dimension(2) :: r
!       effects: Returns fd%ref and fd%o%ref.

!cod$
        r(1) = fd%ref
        r(2) = fd%o%ref
        call glean(fd)
      end function

      function fd_ghost(fd) result(g)
!doc$ function x_ghost(fd) result(g)
        type(fields_td_obj) :: fd
        type(ghost) :: g
!       effects: Returns the ghost of fd.

!cod$
        call my(fd)
        g = fd%o%g
        call glean(thy(fd))
      end function

      function fd_energy(fd) result(e)
!doc$ function x_energy(fd) result(e)
        type(fields_td_obj) :: fd
        real(double) :: e
!       effects: Returns the fd energy.

!cod$
        call my(fd)
        e = fd%o%har_energy + fd%o%exc_energy + fd%o%agp_energy + fd%o%atomic_energy + fd%o%lmcc_energy + fd%o%ccp_energy
        call glean(thy(fd))
      end function

      function fd_potential(fd) result(gp)
!doc$ function x_potential(fd) result(gp)
        type(fields_td_obj) :: fd
        type(gen_potential_obj) :: gp
!       effects: Returns a gen_potential_obj constructed from fd.

!cod$
        call my(fd)
        call my(gen_potential(fd%o%total,fd%o%apot),gp)
        call glean(thy(fd))
        call bequeath(thy(gp))
      end function

      subroutine forces_fd(fd,f)
!doc$ subroutine forces(fd,f)
        type(fields_td_obj) :: fd
        real(double), dimension(:,:), intent(out) :: f
!       modifies: f
!       requires: f be dimension(3,number_of_atoms).
!       effects: Returns unsymmetrized forces due to fd and atomic_operators.
!       errors: Passes errors.

!cod$
        real(double), dimension(:,:,:), pointer :: g2
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: ccd

        call my(fd)

        call my(grid(x_layout(fd%o%ccp),MOD_SCOPE),ccd)

        nullify( g2, c1, c2 )

        if (x_type(fd%o%ccp) /= EMPTY_KIND) then
          call fdel(g2,x_layout(fd%o%ccp),D_TYPE,MOD_SCOPE)
          call take(c1,fd%o%ccp,CDF_KIND)
          call alloc(c2,x_layout(fd%o%ccp),D_TYPE,MOD_SCOPE)
          c2 = c1*g2/eight_pi
          call put(c1,fd%o%ccp,CDF_KIND)
          call put(c2,ccd,CDF_KIND)
          deallocate( g2 )
        end if
        call atomic_forces(fd%o%aden,fd%o%gden,fd%o%xcp,fd%o%ahd,ccd,f) ; if (error()) goto 100

100     if (associated( g2 )) deallocate( g2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(ccd))

        call glean(thy(fd))

        if (error("Exit fields_td_mod::forces_fd")) continue

      end subroutine

      subroutine pressure_fd(fd,p)
!doc$ subroutine pressure(fd,p)
        type(fields_td_obj) :: fd
        real(double), intent(out) :: p
!       effects: Returns pressure contributions due to fd.
!       errors: Passes errors.

!cod$
        real(double) :: p1, p2, p3, p4
        type(grid_obj) :: gr1

        call my(fd)

        call my(grid(x_layout(fd%o%gden),MOD_SCOPE),gr1)

        call hartree_pressure_i(fd%o,p1)                       ; if (error()) goto 100
        call atomic_grid_pressure_i(fd%o,p2)                   ; if (error()) goto 100
        call atomic_pressure(fd%o%aden,fd%o%gden,fd%o%xcp,p3)  ; if (error()) goto 100
        gr1 = fd%o%gden
        call saxpby(1.0_double,gr1,1.0_double,fd%o%axcd)
        call xc_grid_pressure(fd%o%xc,gr1,p4)                  ; if (error()) goto 100
        p = p1 + p2 + p3 + p4  
 
100     call glean(thy(gr1))

        call glean(thy(fd))

        if (error("Exit fields_td_mod::pressure_fd")) continue

      end subroutine

      subroutine stress_tensor_fd(fd,s)
!doc$ subroutine stress_tensor(fd,s)
        type(fields_td_obj) :: fd
        real(double), dimension(:,:), intent(out) :: s
!       requires: s be dimension(3,3).
!       effects: Returns unsymmetrized stress tensor contributions due to fd.
!       errors: Passes errors.

!cod$
        real(double), dimension(:,:), allocatable :: s1, s2, s3, s4
        type(grid_obj) :: gr1

        call my(fd)

        call my(grid(x_layout(fd%o%gden),MOD_SCOPE),gr1)

        allocate( s1(3,3), s2(3,3), s3(3,3), s4(3,3) )

        call hartree_stress_tensor_i(fd%o,s1)                       ; if (error()) goto 100
        call atomic_grid_stress_tensor_i(fd%o,s2)                   ; if (error()) goto 100
        call atomic_stress_tensor(fd%o%aden,fd%o%gden,fd%o%xcp,s3)  ; if (error()) goto 100
        gr1 = fd%o%gden
        call saxpby(1.0_double,gr1,1.0_double,fd%o%axcd)
        call xc_grid_stress_tensor(fd%o%xc,gr1,s4)                  ; if (error()) goto 100
        s = s1 + s2 + s3 + s4

100     if (allocated( s1 )) deallocate( s1 )
        if (allocated( s2 )) deallocate( s2 )
        if (allocated( s3 )) deallocate( s3 )
        if (allocated( s4 )) deallocate( s4 )

        call glean(thy(gr1))

        call glean(thy(fd))

        if (error("Exit fields_td_mod::stress_tensor_fd")) continue

      end subroutine

      subroutine put_smooth_density_fd(f,g)
!doc$ subroutine put_smooth_density(f,g)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects : Returns the grid density.
!       errors: Incompatible layouts.

!cod$ 
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%gden)),"ERROR: different layouts")) goto 100
        g = f%o%gden
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_smooth_density_fd")) continue
      end subroutine

      subroutine put_local_potential_fd(f,g)
!doc$ subroutine put_local_potential(f,g)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects: Returns the local potential.
!       errors: Incompatible layouts.

!cod$ 
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%total)),"ERROR: different layouts")) goto 100
        g = f%o%total
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_local_potential_fd")) continue
      end subroutine

      subroutine put_hartree_potential_fd(f,g)
!doc$ subroutine put_hartree_potential(f,pot)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects: Returns the hartree potential.
!       errors: Incompatible layouts.

!cod$ 
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%hap)),"ERROR: different layouts")) goto 100
        g = f%o%hap
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_hartree_potential_fd")) continue
      end subroutine

      subroutine put_xc_potential_fd(f,g)
!doc$ subroutine put_xc_potential(f,g)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects: Returns the exchange-correlation potential.
!       errors: Incompatible layouts.

!cod$ 
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%xcp)),"ERROR: different layouts")) goto 100
        g = f%o%xcp
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_xc_potential_fd")) continue
      end subroutine

      subroutine put_atomic_xc_density_fd(f,g)
!doc$ subroutine put_atomic_xc_density(f,g)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects: Returns the atomic exchange-correlation density.
!       errors: Incompatible layouts.

!cod$ 
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%axcd)),"ERROR: different layouts")) goto 100
        g = f%o%axcd
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_atomic_xc_density_fd")) continue
      end subroutine

      subroutine put_atomic_hartree_density_fd(f,g)
!doc$ subroutine put_atomic_hartree_density(f,g)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects: Returns the atomic Hartree density.
!       errors: Incompatible layouts.

!cod$
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%ahd)),"ERROR: different layouts")) goto 100
        g = f%o%ahd
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_atomic_hartree_density_fd")) continue
      end subroutine

      subroutine put_atomic_grid_potential_fd(f,g)
!doc$ subroutine put_atomic_grid_potential(f,g)
        type(fields_td_obj) :: f
        type(grid_obj) :: g
!       modifies: g
!       effects: Returns the atomic grid potential.
!       errors: Incompatible layouts.

!cod$ 
        call my(f)
        call my(g)
        if (error(x_ghost(x_layout(g)) /= x_ghost(x_layout(f%o%agp)),"ERROR: different layouts")) goto 100
        g = f%o%agp
100     call glean(thy(f))
        call glean(thy(g))
        if (error("Exit fields_td_mod::put_atomic_grid_potential_fd")) continue
      end subroutine

      function max_energy_fd(fd) result(e)
!doc$ function max_energy(fd) result(e)
        type(fields_td_obj) :: fd
        real(double) :: e
!       effects: Returns the maximum among the field energies.

        call my(fd)
        e = max(abs(fd%o%exc_energy),abs(fd%o%har_energy),abs(fd%o%agp_energy),abs(fd%o%atomic_energy), &
                & abs(fd%o%lmcc_energy + fd%o%ccp_energy))
        call glean(thy(fd))
      end function

      subroutine diary_energies_fd(fd,md)
!doc$ subroutine diary_energies(fd,md)
        type(fields_td_obj) :: fd
        integer, intent(in) :: md
!       modifies: Output stream
!       effects: Prints energies.

!cod$
        call my(fd)

        if (i_access(diaryfile())) then
          select case (md)
          case (12)
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f12.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f12.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f12.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f12.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f12.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          case (13)
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f13.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f13.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f13.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f13.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f13.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          case (14)
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f14.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f14.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f14.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f14.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f14.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          case (15)
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f15.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f15.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f15.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f15.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f15.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          case (16)
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f16.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f16.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f16.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f16.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f16.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          case (17)
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f17.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f17.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f17.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f17.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f17.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          case default
            write(x_unit(diaryfile()),'(t6,"Hartree               = ",f18.9)') fd%o%har_energy
            write(x_unit(diaryfile()),'(t6,"exchange-correlation  = ",f18.9)') fd%o%exc_energy
            write(x_unit(diaryfile()),'(t6,"local pseudopotential = ",f18.9)') fd%o%agp_energy
            write(x_unit(diaryfile()),'(t6,"atomic                = ",f18.9)') fd%o%atomic_energy
            if (fd%o%charge_state /= 0.0_double) then
              write(x_unit(diaryfile()),'(t6,"counter charge        = ",f18.9)') (fd%o%lmcc_energy + fd%o%ccp_energy)
            end if
          end select
        end if

        call glean(thy(fd))

      end subroutine


      subroutine write_restart_fd(fd,nrestf)
!doc$ subroutine write_restart(fd,nrestf)
        type(fields_td_obj) :: fd
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes fd restart information to nrestf.
!       errors: Passes errors.

!cod$
        character(line_len) :: tag

        call my(fd)
        call my(nrestf)

        if (i_access(nrestf)) call startblock(nrestf,"FIELDS")

        call write_restart(fd%o%aden,nrestf) ; if (error()) goto 100
        call write_restart(fd%o%apot,nrestf) ; if (error()) goto 100

        tag = "GRID_DENSITY"
        call write_restart(fd%o%gden,tag,nrestf) ; if (error()) goto 100

        tag = "GRID_POTENTIAL"
        call write_restart(fd%o%total,tag,nrestf) ; if (error()) goto 100

        tag = "SELF-CONSISTENT_POTENTIAL"
        call write_restart(fd%o%scp,tag,nrestf) ; if (error()) goto 100

100     if (i_access(nrestf)) call endblock(nrestf)

        call glean(thy(fd))
        call glean(thy(nrestf))

        if (error("Exit fields_td_mod::write_restart_fd")) continue

      end subroutine

! private routines

      subroutine form_coulomb_kernel_i(fdr,lay)
        type(fields_td_rep) :: fdr
        type(layout_obj) :: lay

        call my(lay)

        call fdelinv(fdr%ck,lay,D_TYPE,MOD_SCOPE)
        call filter(fdr%ck,lay,MOD_SCOPE)

        call glean(thy(lay))
        if (error("Exit fields_td_mod::form_coulomb_kernel_i")) continue

      end subroutine


      subroutine restart_portal_i(fdr,ext,restf)
        type(fields_td_rep) :: fdr
        type(external_obj) :: ext
        type(tagio_obj) :: restf

        logical :: found
        character(line_len) :: tag
        character(1) :: tios

        call my(ext)
        call my(restf)

        if (i_access(restf)) tios = findfirsttag(restf,"FIELDS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: FIELDS block was not found")) goto 200
        if (i_access(restf)) call openblock(restf)

        call my(atomic_density(x_atomic_operators(ext),restf),fdr%aden) ; if (error()) goto 100

        call my(grid(x_layout(ext),MOD_SCOPE),fdr%gden)
        tag = "GRID_DENSITY"
        call read_restart(fdr%gden,tag,restf) ; if (error()) goto 100

        call my(grid(x_layout(ext),MOD_SCOPE),fdr%total)
        tag = "GRID_POTENTIAL"
        call read_restart(fdr%total,tag,restf) ; if (error()) goto 100

        call my(grid(x_layout(ext),MOD_SCOPE),fdr%scp)
        tag = "SELF-CONSISTENT_POTENTIAL"
        call read_restart(fdr%scp,tag,restf) ; if (error()) goto 100

100     if (i_access(restf)) call closeblock(restf)

200     call glean(thy(ext))
        call glean(thy(restf))

        if (error("Exit fields_td_mod::restart_portal_i")) continue

      end subroutine

      subroutine form_compensation_i(fdr,ext)
        type(fields_td_rep) :: fdr
        type(external_obj) :: ext
!       effects: Forms the variables and arrays needed for charge-state calculations.

        logical :: found
        character(1) :: tios
        character(line_len) :: tag
        integer :: i1, i2, i3, ios, m1, m2, m3, n
        integer, dimension(2) :: csr
        integer, dimension(3) :: nd
        real(double), parameter :: tol_rm = 1.0e-10_double, tol_cc = 1.0e-9_double
        real(double) :: cc_local, cc_global, e0, e1, e2, e3, e4, e_local, e_global
        real(double) :: pfd, pfp, ratio, rmt, rsum_global, rsum_local
        real(double), dimension(3) :: s
        real(double), dimension(:,:,:), pointer :: ccd, ccp, ccp0, g2, rm, x, y, z
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(layout_obj) :: lay
        type(lattice_obj) :: lat
        type(crystal_obj) :: bulk_cr
        type(atomic_operators_obj) :: bulk_ao
        type(grid_obj) :: gr1
        type(tagio_obj) :: bulk_restf

        call my(ext)

        call arglc("charge_state_mode",tag,found)
        if (.not.found) tag = "real_number"
        select case (trim(tag))
        case ("real_number")
          call arg("charge_state",fdr%charge_state,found)
          if (.not.found) fdr%charge_state = 0.0_double
        case ("integer_ratio")
          call arg("charge_state_ratio",csr,found)
          if (error(.not.found,"ERROR: charge_state_ratio was not found")) goto 200
          if (error(csr(2) == 0,"ERROR: denominator = 0")) goto 200
          fdr%charge_state = real(csr(1),double)/real(csr(2),double)
        case default
          if (error(.true.,"ERROR: charge_state_mode was not recognized")) goto 200
        end select

        if (fdr%charge_state == 0.0_double) then
          fdr%compensation_method = NONE
          fdr%lmcc_width = 0.0_double
          fdr%lmcc_site = 0.0_double
          fdr%lmcc_energy = 0.0_double
          fdr%lmcc_fixed_energy = 0.0_double
          fdr%ccp_energy = 0.0_double
          call my(grid(x_layout(ext),MOD_SCOPE),fdr%ccp)
          goto 200
        end if

        call arglc("compensation",tag,found)
        if (.not.found) tag = "ub"
        select case (trim(tag))
        case ("ub")
          fdr%compensation_method = UBC
          fdr%lmcc_width = 0.0_double
          fdr%lmcc_site = 0.0_double
          fdr%lmcc_energy = 0.0_double
          fdr%lmcc_fixed_energy = 0.0_double
          fdr%ccp_energy = 0.0_double
          call my(grid(x_layout(ext),MOD_SCOPE),fdr%ccp)
          goto 200
        case ("lmcc")
          fdr%compensation_method = LMCC
          if (error(x_type(x_atomic_operators(ext)) == PAW,"ERROR: lmcc is not currently supported with paw")) goto 200
          call arg("lmcc_site",fdr%lmcc_site,found)
          if (error(.not.found,"ERROR: lmcc_site tag was not found")) goto 200
          s = fdr%lmcc_site
          call centralize_position_i(fdr%lmcc_site)
          if (any(fdr%lmcc_site /= s)) call warn("WARNING: the counter charge was translated into the central parallelpiped")
          if (error(.not.invariant_site(x_space_group(ext),fdr%lmcc_site),"ERROR: lmcc_site is not invariant")) goto 200
          call arg("lmcc_width",fdr%lmcc_width,found)
          if (.not.found) fdr%lmcc_width = 1.50_double
          if (error(fdr%lmcc_width < 0.50_double,"ERROR: lmcc_width is too small")) goto 200
        case default
          if (error(.true.,"ERROR: compensation tag was not recognized")) goto 200
        end select

        call my(x_layout(ext),lay)
        call my(x_lattice(x_crystal(ext)),lat)
        call my(grid(lay,MOD_SCOPE),gr1)

        nullify( x, y, z, rm, ccd, ccp, ccp0, g2, r1, r2, c1, c2 )

        call mesh(x,y,z,lay,D_TYPE,MOD_SCOPE)
        call alloc(rm,lay,D_TYPE,MOD_SCOPE)

        s = lat2r(lat,fdr%lmcc_site)                                     ! DISTANCES TO CC_SITE FROM CENTRAL CELL MESH POINTS
        do i3 = 1,size(rm,3)
        do i2 = 1,size(rm,2)
        do i1 = 1,size(rm,1)
          rm(i1,i2,i3) = norm((/x(i1,i2,i3),y(i1,i2,i3),z(i1,i2,i3)/)-s)
        end do
        end do
        end do

        do m1 = -1,+1                                                    ! DISTANCES TO THE NEAREST CC_SITE
        do m2 = -1,+1
        do m3 = -1,+1
          if ( all((/i1,i2,i3/) == 0)) cycle
          s = lat2r(lat,fdr%lmcc_site) + lat2r(lat,real((/m1,m2,m3/),double))
          do i3 = 1,size(rm,3)
          do i2 = 1,size(rm,2)
          do i1 = 1,size(rm,1)
            rmt = norm((/x(i1,i2,i3),y(i1,i2,i3),z(i1,i2,i3)/)-s)
            rm(i1,i2,i3) = min(rm(i1,i2,i3),rmt)
          end do
          end do
          end do
        end do
        end do
        end do
        deallocate( x, y, z )

        pfd = -fdr%charge_state/((pi**1.5_double)*(fdr%lmcc_width**3))   ! GAUSSIAN CHARGE DISTRIBUTION AT CC_SITE (ccd)
        pfp = -2.0_double*fdr%charge_state                               !  AND CORRESPONDING ERROR FUNCTION POTENTIAL (ccp)
        call alloc(ccd,lay,D_TYPE,MOD_SCOPE)
        call alloc(ccp,lay,D_TYPE,MOD_SCOPE)
        do i3 = 1,size(rm,3)
        do i2 = 1,size(rm,2)
        do i1 = 1,size(rm,1)
          ratio = rm(i1,i2,i3)/fdr%lmcc_width
          ccd(i1,i2,i3) = pfd*exp(-ratio**2)
          if (rm(i1,i2,i3) .in. nbhd(0.0_double,tol_rm)) then
            ccp(i1,i2,i3) = pfp*(2.0_double/(sqrt(pi)*fdr%lmcc_width))
          else
            ccp(i1,i2,i3) = pfp*(1.0_double - complementary_error(ratio))/rm(i1,i2,i3)
          end if
        end do
        end do
        end do

        nd = x_dims(lay)
        n = nd(1)*nd(2)*nd(3)

        cc_local = -sum(ccd)                                      ! CHECK THAT THE COUNTER CHARGE IS WITHIN THE SUPERCELL
        call allreduce(MOD_SCOPE,MPI_SUM,cc_local,cc_global)
        cc_global = cc_global*x_cell_volume(lat)/real(n,double)
        if (error(cc_global .out. nbhd(fdr%charge_state,tol_cc),"ERROR: counter charge extends outside the supercell")) goto 100

        rsum_local = sum(ccp)                                     ! ERROR FUNCTION POTENTIAL WITH AVERAGE VALUE = 0 (ccp0)
        call allreduce(MOD_SCOPE,MPI_SUM,rsum_local,rsum_global)  !  AND ASSOCIATED ENERGY CONTRIBUTION: -0.5*q*ccp0
        call alloc(ccp0,lay,D_TYPE,MOD_SCOPE)
        ccp0 = ccp - rsum_global/real(n,double)
        e0 = -0.5_double*fdr%charge_state*rsum_global/real(n,double)

        call my(grid(lay,MOD_SCOPE),fdr%ccp)                      ! DIFFERENCE POTENTIAL: fdr%ccp = ccp0 + eight_pi*del^-2(ccd)
        call put(ccp0,fdr%ccp,RD_KIND)                            !  WITH AVERAGE VALUE = 0
        call transform(fdr%ccp,CDF_KIND)
        call put(ccd,gr1,RD_KIND)                                 ! NOTE THAT ccp0 IS FOURIER TRANSFORMED HERE. THIS MAY NOT BE
        call take(c1,gr1,CDF_KIND)                                !  ADVISEABLE BECAUSE IT HAS A DERIVATIVE DISCONTINUITY AT THE
        call alloc(c2,lay,D_TYPE,MOD_SCOPE)                     !  BOUNDARIES OF A WIGNER-SEITZ CELL SURROUNDING THE GAUSSIAN.
        c2 = eight_pi*c1*fdr%ck
        call put(c2,gr1,CDF_KIND)
        call saxpby(1.0_double,fdr%ccp,-1.0_double,gr1)
        call filter(fdr%ccp)
        call put(c1,gr1,CDF_KIND)
        call take(ccd,gr1,RD_KIND)

        call put(ccd,gr1,RD_KIND)                                 ! SELF-INTERACTION ENERGY CORRECTION: -0.5*ccd*fdr%ccp
        call take(c1,gr1,CDF_KIND)
        gr1 = fdr%ccp
        call take(c2,gr1,CDF_KIND)
        e_local = real(sum(conjg(c1)*c2),double)
        call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
        e1 = -0.5_double*e_global*x_cell_volume(lat)
        call put(c1,gr1,CDF_KIND)
        call take(ccd,gr1,RD_KIND)
        deallocate( c2 )

        gr1 = fdr%agp                                             ! INTERACTION ENERGY OF THE DEFECT CELL ION DENSITY WITH THE
        call take(c1,gr1,CDF_KIND)                                !  WITH THE DIFFERENCE POTENTIAL: ion_density*fdr%ccp
        call fdel(g2,lay,D_TYPE,MOD_SCOPE)
        call alloc(c2,lay,D_TYPE,MOD_SCOPE)
        c2 = g2*c1/eight_pi
        call put(c1,gr1,CDF_KIND)
        gr1 = fdr%ccp
        call take(c1,gr1,CDF_KIND)
        e_local = real(sum(c1*conjg(c2)),double)
        call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
        e2 = e_global*x_cell_volume(lat)
        deallocate( g2, c1, c2 )

        call arglc("lmcc_bulk",tag,found)                         ! BULK DENSITIES AND ASSOCIATED ENERGY TERMS
        if (.not.found) tag = "on"
        select case (trim(tag))
        case ("off")
          e3 = 0.0_double
          e4 = 0.0_double
        case ("on")
          call my(tagio(trim(bulk_restart_path),TAGIO_READ,mkey,len(mkey)),bulk_restf)   ! BULK ELECTRON DENSITY: -density*fdr%ccp
          if (i_access(bulk_restf)) ios = x_tagfd(bulk_restf)
          if (i_comm(bulk_restf)) call broadcast(FILE_SCOPE,ios)
          if (error(ios == 0,"ERROR: bulk restart file was not found")) goto 100
          if (i_access(bulk_restf)) tios = findfirsttag(bulk_restf,"FIELDS")
          if (i_comm(bulk_restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: FIELDS block was not found")) goto 100
          if (i_access(bulk_restf)) call openblock(bulk_restf)
          tag = "GRID_DENSITY"
          call read_restart(gr1,tag,bulk_restf) ; if (error()) goto 100
          if (i_access(bulk_restf)) call closeblock(bulk_restf)
          call glean(thy(bulk_restf)) ; if (error()) goto 100
          call take(c1,gr1,CDF_KIND)
          gr1 = fdr%ccp
          call take(c2,gr1,CDF_KIND)
          e_local = real(sum(conjg(c1)*c2),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          e3 = -e_global*x_cell_volume(lat)
          deallocate( c1, c2 )
          call my(crystal(bulk_crystal_path),bulk_cr) ; if (error()) goto 100           ! BULK ION DENSITY: -ion_density*fdr%ccp
          call my(atomic_operators(bulk_cr,lay,x_space_group(ext),x_xc_type(ext)),bulk_ao) ; if (error()) goto 100
          gr1 = atomic_grid_potential(bulk_ao,lay) ; if (error()) goto 100
          call glean(thy(bulk_cr))
          call glean(thy(bulk_ao))
          call take(c1,gr1,CDF_KIND)
          call fdel(g2,lay,D_TYPE,MOD_SCOPE)
          call alloc(c2,lay,D_TYPE,MOD_SCOPE)
          c2 = g2*c1/eight_pi
          call put(c1,gr1,CDF_KIND)
          gr1 = fdr%ccp
          call take(c1,gr1,CDF_KIND)
          e_local = real(sum(c1*conjg(c2)),double)
          call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
          e4 = -e_global*x_cell_volume(lat)
          deallocate( g2, c1, c2 )
        end select

        fdr%lmcc_fixed_energy = e0 + e1 + e3 + e4
        fdr%lmcc_energy = e0 + e1 + e2 + e3 + e4

        call arglc("lmcc_energy",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on")
          if (i_access(diaryfile())) write(37,'("e0 (-0.5*charge_state*average_error_function_potential) = ",f15.10)') e0
          if (i_access(diaryfile())) write(37,'("e1 (-0.5*ccd*difference_potential)                      = ",f15.10)') e1
          if (i_access(diaryfile())) write(37,'("e2 (+ionic_density*difference_potential)                = ",f15.10)') e2
          if (i_access(diaryfile())) write(37,'("e3 (-bulk_density*difference_potential)                 = ",f15.10)') e3
          if (i_access(diaryfile())) write(37,'("e4 (-bulk_ionic_density*difference_potential)           = ",f15.10)') e4
          if (i_access(diaryfile())) write(37,'("lmcc_energy = ",f15.10)') fdr%lmcc_energy
        end select

100     if (associated( x )) deallocate( x )
        if (associated( y )) deallocate( y )
        if (associated( z )) deallocate( z )
        if (associated( rm )) deallocate( rm )
        if (associated( ccd )) deallocate( ccd )
        if (associated( ccp )) deallocate( ccp )
        if (associated( ccp0 )) deallocate( ccp0 )
        if (associated( g2 )) deallocate( g2 )
        if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(lay))
        call glean(thy(lat))
        call glean(thy(gr1))

200    call glean(thy(ext))

       if (error("Exit fields_td_mod::form_compensation_i")) continue

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

      subroutine diary_construction_i(fdr,ao)
        type(fields_td_rep) :: fdr
        type(atomic_operators_obj) :: ao

        call my(ao)

        if (i_access(diaryfile())) then

          write(x_unit(diaryfile()),'(/,"Fields object construction:")')

          write(x_unit(diaryfile()),'(/,t4,"Input-density initialization")')

          call diary_angular_mesh(ao)

          if (fdr%charge_state /= 0.0_double) then
            write(x_unit(diaryfile()),'(/,t4,sp,"Charge state = ",f0.6)') fdr%charge_state
            select case (fdr%compensation_method)
            case (UBC)
              write(x_unit(diaryfile()),'(/,t6,"Using the uniform-background-charge method")')
            case (LMCC)
              write(x_unit(diaryfile()),'(/,t6,"Using the local-moment-counter-charge method:")')
              write(x_unit(diaryfile()),'(/,t8,"position {a1,a2,a3} = ",3f14.10)') fdr%lmcc_site
              write(x_unit(diaryfile()),'(t8,"width = ",f10.6)') fdr%lmcc_width
            end select
          end if

        end if

      call glean(thy(ao))

      end subroutine

      subroutine update_compensation_i(ext,fdr)
        type(external_obj) :: ext
        type(fields_td_rep) :: fdr
!       effects: Updates fdr%lmcc_energy after a change in the atom positions.

        real(double) :: e2, e_local, e_global
        real(double), dimension(:,:,:), pointer :: g2
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(layout_obj) :: lay
        type(lattice_obj) :: lat
        type(grid_obj) :: gr1

        call my(ext)

        if (fdr%charge_state == 0.0_double) goto 200

        select case (fdr%compensation_method)
        case (UBC)
          goto 200
        end select

        call my(x_layout(ext),lay)
        call my(x_lattice(x_crystal(ext)),lat)
        call my(grid(lay,MOD_SCOPE),gr1)

        nullify( g2, c1, c2 )

        gr1 = fdr%agp                                             ! INTERACTION ENERGY OF THE DEFECT CELL ION DENSITY WITH THE
        call take(c1,gr1,CDF_KIND)                                !  WITH THE DIFFERENCE POTENTIAL: ion_density*fdr%ccp
        call fdel(g2,lay,D_TYPE,MOD_SCOPE)
        call alloc(c2,lay,D_TYPE,MOD_SCOPE)
        c2 = g2*c1/eight_pi
        call put(c1,gr1,CDF_KIND)
        gr1 = fdr%ccp
        call take(c1,gr1,CDF_KIND)
        e_local = real(sum(c1*conjg(c2)),double)
        call allreduce(MOD_SCOPE,MPI_SUM,e_local,e_global)
        e2 = e_global*x_cell_volume(lat)
        deallocate( g2, c1, c2 )

        fdr%lmcc_energy = fdr%lmcc_fixed_energy + e2

100     if (associated( g2 )) deallocate( g2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(lay))
        call glean(thy(lat))
        call glean(thy(gr1))

200    call glean(thy(ext))

        if (error("Exit fields_td_mod::update_compensation_i")) continue

      end subroutine

      subroutine own_i(fd)
        type(fields_td_obj) :: fd
        type(fields_td_obj) :: fdt
        if (fd%ref < fd%o%ref) then
          allocate( fdt%o )
          fdt%o%ref                 = 0
          fdt%o%g                   = fd%o%g
          fdt%o%g_ao                = fd%o%g_ao
          fdt%o%compensation_method = fd%o%compensation_method
          fdt%o%exc_energy          = fd%o%exc_energy
          fdt%o%har_energy          = fd%o%har_energy
          fdt%o%agp_energy          = fd%o%agp_energy
          fdt%o%atomic_energy       = fd%o%atomic_energy
          fdt%o%charge_state        = fd%o%charge_state
          fdt%o%lmcc_width          = fd%o%lmcc_width
          fdt%o%lmcc_site           = fd%o%lmcc_site
          fdt%o%lmcc_energy         = fd%o%lmcc_energy
          fdt%o%lmcc_fixed_energy   = fd%o%lmcc_fixed_energy
          fdt%o%ccp_energy          = fd%o%ccp_energy
          call my(fd%o%gden,fdt%o%gden)
          call my(fd%o%ahd,fdt%o%ahd)
          call my(fd%o%axcd,fdt%o%axcd)
          call my(fd%o%agp,fdt%o%agp)
          call my(fd%o%hap,fdt%o%hap)
          call my(fd%o%xcp,fdt%o%xcp)
          call my(fd%o%scp,fdt%o%scp)
          call my(fd%o%ccp,fdt%o%ccp)
          call my(fd%o%total,fdt%o%total)
          call my(fd%o%aden,fdt%o%aden)
          call my(fd%o%apot,fdt%o%apot)
          call my(fd%o%xc,fdt%o%xc)
          allocate( fdt%o%ck(size(fd%o%ck,1),size(fd%o%ck,2),size(fd%o%ck,3)) )
          fdt%o%ck = fd%o%ck
          fd%o%ref = fd%o%ref - fd%ref
          fd%o => fdt%o
          fd%o%ref = fd%o%ref + fd%ref
        end if
      end subroutine

      subroutine hartree_pressure_i(fdr,p)
        type(fields_td_rep) :: fdr
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the Hartree potential.

        real(double) :: p_local
        real(double), dimension(:,:,:), pointer :: r1
        complex(double), dimension(:,:,:), pointer :: c1
        type(layout_obj) :: lay

        nullify( r1, c1 )

        call my(x_layout(fdr%gden),lay)
        call take(c1,fdr%gden,CDF_KIND)
        call alloc(r1,lay,D_TYPE,MOD_SCOPE)
        r1 = real(c1*conjg(c1)*fdr%ck,double)
        call put(c1,fdr%gden,CDF_KIND)
        p_local = (four_pi/3.0_double)*sum(r1)
        call allreduce(MOD_SCOPE,MPI_SUM,p_local,p)

100     if (associated( r1 )) deallocate( r1 )
        if (associated( c1 )) deallocate( c1 )

        call glean(thy(lay))

        if (error("Exit fields_td_mod::hartree_pressure_i")) continue

      end subroutine

      subroutine atomic_grid_pressure_i(fdr,p)
        type(fields_td_rep) :: fdr
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the atomic-grid potential.

        real(double) :: p_local
        complex(double), dimension(:,:,:), pointer :: c1, c2

        nullify( c1, c2 )

        call take(c1,fdr%agp,CDF_KIND)
        call take(c2,fdr%gden,CDF_KIND)
        p_local = real(sum(c1*conjg(c2)),double)
        call put(c1,fdr%agp,CDF_KIND)
        call put(c2,fdr%gden,CDF_KIND)
        call allreduce(MOD_SCOPE,MPI_SUM,p_local,p)

100     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        if (error("Exit fields_td_mod::atomic_grid_pressure_i")) continue

      end subroutine

      subroutine hartree_stress_tensor_i(fdr,s)
        type(fields_td_rep) :: fdr
        real(double), dimension(:,:), intent(out) :: s
!       effects: Returns stress tensor contributions due to the Hartree potential.

        real(double) :: sum_r
        real(double), dimension(:,:), allocatable :: s_local
        real(double), dimension(:,:,:), pointer :: gx, gy, gz, r1
        complex(double), dimension(:,:,:), pointer :: c1
        type(layout_obj) :: lay

        call my(x_layout(fdr%gden),lay)

        nullify( gx, gy, gz, r1, c1 )

        call fmesh(gx,gy,gz,lay,D_TYPE,MOD_SCOPE)

        call alloc(r1,lay,D_TYPE,MOD_SCOPE)
        call take(c1,fdr%gden,CDF_KIND)
        r1 = real(c1*conjg(c1)*fdr%ck,double)
        call put(c1,fdr%gden,CDF_KIND)
        sum_r = sum(r1)
        r1 = r1*fdr%ck

        allocate( s_local(3,3) )
        s_local(1,1) = four_pi*(2.0_double*sum(r1*gx*gx) - sum_r)
        s_local(1,2) = four_pi*(2.0_double*sum(r1*gx*gy)        )
        s_local(1,3) = four_pi*(2.0_double*sum(r1*gx*gz)        )
        s_local(2,2) = four_pi*(2.0_double*sum(r1*gy*gy) - sum_r)
        s_local(2,3) = four_pi*(2.0_double*sum(r1*gy*gz)        )
        s_local(3,3) = four_pi*(2.0_double*sum(r1*gz*gz) - sum_r)
        s_local(2,1) = s_local(1,2)
        s_local(3,1) = s_local(1,3)
        s_local(3,2) = s_local(2,3)
        call allreduce(MOD_SCOPE,MPI_SUM,s_local,s)

100     if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )
        if (associated( r1 )) deallocate( r1 )
        if (associated( c1 )) deallocate( c1 )
        if (allocated( s_local )) deallocate( s_local )

        call glean(thy(lay))

        if (error("Exit fields_td_mod::hartree_stress_tensor_i")) continue

      end subroutine

      subroutine atomic_grid_stress_tensor_i(fdr,s)
        type(fields_td_rep) :: fdr
        real(double), dimension(:,:), intent(out) :: s
!       effects: Returns stress tensor contributions due to the atomic grid potential.

        real(double) :: sum_r
        real(double), dimension(:,:), allocatable :: s_local
        complex(double), dimension(:,:,:), pointer :: c1, c2

        nullify( c1, c2 )

        call take(c1,fdr%agp,CDF_KIND)
        call take(c2,fdr%gden,CDF_KIND)
        sum_r = -real(sum(c1*conjg(c2)),double)
        call put(c1,fdr%agp,CDF_KIND)
        call put(c2,fdr%gden,CDF_KIND)

        allocate( s_local(3,3) )
        s_local(1,1) = sum_r
        s_local(1,2) = 0.0_double
        s_local(1,3) = 0.0_double
        s_local(2,2) = sum_r
        s_local(2,3) = 0.0_double
        s_local(3,3) = sum_r
        s_local(2,1) = s_local(1,2)
        s_local(3,1) = s_local(1,3)
        s_local(3,2) = s_local(2,3)
        call allreduce(MOD_SCOPE,MPI_SUM,s_local,s)

100     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )
        if (allocated( s_local )) deallocate( s_local )

        if (error("Exit fields_td_mod::atomic_grid_stress_tensor_i")) continue

      end subroutine

      end module

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module fields_fh_mod
!doc$ module fields_fh_mod

!     One datatype is available here: type(fields_fh_obj)

!     fields_fh_mod manages field quantities in calculations with a fixed hamiltonian.

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
      use lattice_mod
      use crystal_mod
      use atoms_mod
      use external_mod
      use symmetry_mod
      use gen_potential_mod
      use grid_mod
      use atomic_operators_mod
      use atomic_density_mod
      use atomic_potential_mod
      use timing_mod

!cod$
      implicit none
      private 

      integer, parameter :: MOD_SCOPE = CONFIG

      integer, parameter :: NA   = 0
      integer, parameter :: UBC  = 1
      integer, parameter :: LMCC = 2

      type :: fields_fh_rep
        integer :: ref
        type(ghost) :: g
        integer :: compensation_method                     ! compensation method
        real(double) :: charge_state                       ! charge state of supercell
        real(double) :: lmcc_width                         ! width of the Gaussian charge distribution in the LMCC method
        real(double), dimension(3) :: lmcc_site            ! site of the Gaussian charge distribution in the LMCC method
        type(grid_obj) :: total                            ! total potential
        type(atomic_potential_obj) :: apot                 ! atomic potential object
      end type

      type, public :: fields_fh_obj
        private
        integer :: ref
        type(fields_fh_rep), pointer :: o
      end type

!doc$
      public :: fields_fh
      public :: my
      public :: thy
      public :: bequeath
      public :: glean
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_potential

!cod$
      interface fields_fh
        module procedure constructor_fd
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
      interface x_potential
        module procedure fd_potential
      end interface

      contains

! public routines

      function constructor_fd(ext,restf) result(fd)
!doc$ function fields_fh(ext,restf) result(fd)
        type(external_obj) :: ext
        type(tagio_obj) :: restf
        type(fields_fh_obj) :: fd
!       requires: ext and input files be consistent with restf.
!       effects: Constructs a new fd.
!       errors: Passes errors.

!cod$ 
        character(1) :: tios
        character(line_len) :: tag
        integer(long) :: dsize, iosl, ndata, s4
        type(atomic_density_obj) :: aden

        call start_timer("fields_fh: constructor")

        call my(ext)
        call my(restf)

        fd%ref = 0
        allocate( fd%o )
        fd%o%ref = 0
        fd%o%g = x_ghost()

        ! open the FIELDS block
        if (i_access(restf)) tios = findfirsttag(restf,"FIELDS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: FIELDS block was not found")) goto 300
        if (i_access(restf)) call openblock(restf)

        ! read the atomic potential
        call my(atomic_density(x_atomic_operators(ext),restf),aden) ; if (error()) goto 200
        call my(atomic_potential(aden,restf=restf),fd%o%apot)
        call glean(thy(aden)) ; if (error()) goto 200

        ! read the grid potential
        call my(grid(x_layout(ext),MOD_SCOPE),fd%o%total)
        tag = "GRID_POTENTIAL"
        call read_restart(fd%o%total,tag,restf) ; if (error()) goto 200

        ! open the COMPENSATION block
        if (i_access(restf)) tios = findfirsttag(restf,"COMPENSATION")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: COMPENSATION block was not found")) goto 200
        if (i_access(restf)) call openblock(restf)

        ! read the charge state
        if (i_access(restf)) tios = findfirsttag(restf,"CHARGE_STATE")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: CHARGE_STATE tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_double ; ndata = 1
          call readf(fd%o%charge_state,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,fd%o%charge_state)

        if (fd%o%charge_state == 0.0_double) then

          fd%o%compensation_method = NA
          fd%o%lmcc_width = 0.0_double
          fd%o%lmcc_site = 0.0_double

        else

          ! read the compensation method
          if (i_access(restf)) tios = findfirsttag(restf,"COMPENSATION_METHOD")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: COMPENSATION_METHOD tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long ; ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            fd%o%compensation_method = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,fd%o%compensation_method)

          select case (fd%o%compensation_method)
          case (UBC)

            fd%o%lmcc_site = 0.0_double
            fd%o%lmcc_width = 0.0_double

          case (LMCC)

            ! read the lmcc site
            if (i_access(restf)) tios = findfirsttag(restf,"LMCC_SITE")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: LMCC_SITE tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_double ; ndata = 3
              call readf(fd%o%lmcc_site,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,fd%o%lmcc_site)

            ! read the lmcc width
            if (i_access(restf)) tios = findfirsttag(restf,"LMCC_WIDTH")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: LMCC_WIDTH tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_double ; ndata = 1
              call readf(fd%o%lmcc_width,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,fd%o%lmcc_width)

          end select

        end if

        ! close the COMPENSATION block
100     if (i_access(restf)) call closeblock(restf)

        ! close the FIELDS block
200     if (i_access(restf)) call closeblock(restf)
        if (error()) goto 300

        call diary_construction_i(fd%o,x_atomic_operators(ext))

300     call glean(thy(ext))
        call glean(thy(restf))

        if (error("Exit fields_fh_mod::constructor_fd")) continue

        if (.not.error()) call stop_timer("fields_fh: constructor")

      end function

      subroutine my_fd(fd)
!doc$ subroutine my(fd)
        type(fields_fh_obj) :: fd

!cod$
        fd%ref = fd%ref + 1
        fd%o%ref = fd%o%ref + 1
      end subroutine

      subroutine my_new_fd(fdi,fd)
!doc$ subroutine my(fdi,fd)
        type(fields_fh_obj) :: fdi, fd

!cod$
        fd%ref = 1
        fd%o => fdi%o
        fd%o%ref = fd%o%ref + 1
      end subroutine

      function thy_fd(fd) result(fdo)
!doc$ function thy(fd) result(fdo)
        type(fields_fh_obj) :: fd, fdo

!cod$
        fd%ref = fd%ref - 1
        fd%o%ref = fd%o%ref - 1
        fdo%ref = fd%ref
        fdo%o => fd%o
      end function

      subroutine glean_fd(fd)
!doc subroutine glean(fd)
        type(fields_fh_obj) :: fd

!cod$
        if (fd%o%ref < 1) then
          call glean(thy(fd%o%total))
          call glean(thy(fd%o%apot))
          deallocate( fd%o )
        end if
      end subroutine

      subroutine bequeath_fd(fd)
!doc$ subroutine bequeath(fd)
        type(fields_fh_obj) :: fd

!cod$
        continue
      end subroutine

      subroutine assign_fd(fd,fd2)
!doc$ subroutine assign(fd,fd2)
        type(fields_fh_obj), intent(inout) :: fd
        type(fields_fh_obj), intent(in) :: fd2

!cod$
        type(fields_fh_obj) :: fdt
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
        type(fields_fh_obj) :: fd
        integer, dimension(2) :: r
!       effects: Returns fd%ref and fd%o%ref.

!cod$
        r(1) = fd%ref
        r(2) = fd%o%ref
        call glean(fd)
      end function

      function fd_ghost(fd) result(g)
!doc$ function x_ghost(fd) result(g)
        type(fields_fh_obj) :: fd
        type(ghost) :: g
!       effects: returns the ghost of fd.

!cod$
        call my(fd)
        g = fd%o%g
        call glean(thy(fd))
      end function

      function fd_potential(fd) result(gp)
!doc$ function x_potential(fd) result(gp)
        type(fields_fh_obj) :: fd
        type(gen_potential_obj) :: gp
!       effects: Returns a gen_potential_obj constructed from fd.

!cod$
        call my(fd)
        call my(gen_potential(fd%o%total,fd%o%apot),gp)
        call glean(thy(fd))
        call bequeath(thy(gp))
      end function

! private routines

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
        type(fields_fh_rep) :: fdr
        type(atomic_operators_obj) :: ao

        call my(ao)

        if (i_access(diaryfile())) then

          write(x_unit(diaryfile()),'(/,"Fields object construction:")')

          write(x_unit(diaryfile()),'(/,t4,"Restart-potential initialization")')

          call diary_angular_mesh(ao)

          if (fdr%charge_state /= 0.0_double) then
            write(x_unit(diaryfile()),'(/,t4,sp,"Charge state = ",f0.2)') fdr%charge_state
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

      end module

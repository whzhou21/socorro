! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module gen_density_mod
!doc$ module gen_density_mod

!     Defines a type(gen_density_obj), which encapsulates the density data involved in a calculation.  Data is "density data"
!     if it is obtained by summing a contribution from the wavefunctions at each k-point, symmetrized, mixed, and then used
!     to calculate the generalized potential. Currently, the gen_density_obj is a container for a type(grid_obj) and a
!     type(atomic_density_obj).

      use kind_mod
      use error_mod
      use ghost_mod
      use mpi_mod
      use grid_mod
      use layout_mod
      use lattice_mod
      use external_mod
      use symmetry_mod
      use atomic_density_mod

!cod$
      implicit none
      private

      type :: gen_density_rep
        integer :: ref
        type(ghost) :: g
        type(grid_obj) :: gden
        type(atomic_density_obj) :: aden
      end type

      type, public :: gen_density_obj
        private
        integer :: ref
        type(gen_density_rep), pointer :: o
      end type

!doc$
      public :: gen_density
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_grid_density
      public :: x_atomic_density
      public :: distance
      public :: spin_information
      public :: merge_density
      public :: symmetrize
      public :: filter
      public :: merge_symmetrize_filter

!cod$
      interface gen_density
        module procedure constructor_gd_1, constructor_gd_2
      end interface
      interface update
        module procedure update_gd
      end interface
      interface my
        module procedure my_gd, my_new_gd
      end interface
      interface thy
        module procedure thy_gd
      end interface
      interface glean
        module procedure glean_gd
      end interface
      interface bequeath
        module procedure bequeath_gd
      end interface
      interface assignment(=)
        module procedure assign_gd
      end interface
      interface x_ref
        module procedure gd_ref
      end interface
      interface x_ghost
        module procedure gd_ghost
      end interface
      interface x_grid_density
        module procedure gd_grid_density
      end interface
      interface x_atomic_density
        module procedure gd_atomic_density
      end interface
      interface distance
        module procedure distance_gd
      end interface
      interface spin_information
        module procedure spin_information_gd
      end interface
      interface merge_density
        module procedure merge_density_gd
      end interface
      interface symmetrize
        module procedure symmetrize_gd
      end interface
      interface filter
        module procedure filter_gd
      end interface
      interface merge_symmetrize_filter
        module procedure merge_symmetrize_filter_gd
      end interface

      contains

! public routines

      function constructor_gd_1(ext) result(gd)
!doc$ function gen_density(ext) result(gd)
        type(external_obj) :: ext
        type(gen_density_obj) :: gd
!       effects: Constructs a new zeroed gd with KGROUP scope.

!cod$
        logical, parameter :: empty = .true.
        real(double), dimension(:,:,:), pointer :: r1
        type(grid_obj) :: gden
        type(atomic_density_obj) :: aden

        call my(ext)

        gd%ref = 0
        allocate( gd%o )
        gd%o%ref = 0
        gd%o%g = x_ghost()

        nullify( r1 )

        ! construct a zeroed grid density
        call my(grid(x_layout(ext),KGROUP),gden)
        call alloc(r1,x_layout(gden),S_TYPE)
        r1 = 0.0_double
        call put(r1,gden,RS_KIND)

        ! construct a zeroed atomic density
        call my(atomic_density(x_atomic_operators(ext),empty=empty),aden)

        ! construct gen_density from gden and aden
        call my(gden,gd%o%gden)
        call my(aden,gd%o%aden)

100     if (associated( r1 )) deallocate( r1 )

        call glean(thy(gden))
        call glean(thy(aden))

        call glean(thy(ext))

        if (error("Exit gen_density_mod::constructor_gd_1")) continue

      end function

      function constructor_gd_2(gden,aden) result(gd)
!doc$ function gen_density(gden,aden) result(gd)
        type(grid_obj) :: gden
        type(atomic_density_obj) :: aden
        type(gen_density_obj) :: gd
!       effects: Constructs a new gd containing gden and aden.

!cod$
        call my(gden)
        call my(aden)

        gd%ref = 0
        allocate( gd%o )
        gd%o%ref = 0
        gd%o%g = x_ghost()

        call my(gden,gd%o%gden)
        call my(aden,gd%o%aden)

100     call glean(thy(gden))
        call glean(thy(aden))

        if (error("Exit gen_density_mod::constructor_gd_2")) continue

      end function

      subroutine update_gd(gd,gden,aden)
!doc$ subroutine update(gd,gden,aden)
        type(gen_density_obj) :: gd
        type(grid_obj), optional :: gden
        type(atomic_density_obj), optional :: aden
!       modifies: gd
!       effects: Updates gd with respect to gden and aden.

!cod$
        call my(gd)

        if (present(gden)) call my(gden)
        if (present(aden)) call my(aden)

        if (present(gden)) then
          if (x_ghost(gd%o%gden) /= x_ghost(gden) ) then
            call own_i(gd)
            gd%o%g = x_ghost()
            gd%o%gden = gden
          end if
        end if

        if (present(aden)) then
          if (x_ghost(gd%o%aden) /= x_ghost(aden) ) then
            call own_i(gd)
            gd%o%g = x_ghost()
            gd%o%aden = aden
          end if
        end if

        if (present(gden)) call glean(thy(gden))
        if (present(aden)) call glean(thy(aden))

        call glean(thy(gd))

      end subroutine

      subroutine my_gd(gd)
!doc$ subroutine my(gd)
        type(gen_density_obj) :: gd

!cod$
        gd%ref = gd%ref + 1
        gd%o%ref = gd%o%ref + 1
      end subroutine

      subroutine my_new_gd(gdi,gd)
!doc$ subroutine my(gdi,gd)
        type(gen_density_obj) :: gdi, gd

!cod$
        gd%ref = 1
        gd%o => gdi%o
        gd%o%ref = gd%o%ref + 1
      end subroutine

      function thy_gd(gd) result(gdo)
!doc$ function thy(gd) result(gdo)
        type(gen_density_obj) :: gd
        type(gen_density_obj) :: gdo

!cod$
        gd%ref = gd%ref - 1
        gd%o%ref = gd%o%ref - 1
        gdo%ref = gd%ref
        gdo%o => gd%o
      end function

      subroutine glean_gd(gd)
!doc$ subroutine glean(gd)
        type(gen_density_obj) :: gd

!cod$
        if (gd%o%ref < 1) then
          call glean(thy(gd%o%gden))
          call glean(thy(gd%o%aden))
          deallocate( gd%o )
        end if
      end subroutine

      subroutine bequeath_gd(gd)
!doc$ subroutine bequeath(gd)
        type(gen_density_obj) :: gd

!cod$
        continue
      end subroutine

      subroutine assign_gd(gd,gd2)
!doc$ subroutine assignment(=)(gd,gd2)
        type(gen_density_obj), intent(inout) :: gd
        type(gen_density_obj), intent(in) :: gd2

!cod$
        type(gen_density_obj) :: gdt

       call my(gd2)
        gdt%o => gd%o
        gd%o%ref = gd%o%ref - gd%ref
        gd%o => gd2%o
        gd%o%ref = gd%o%ref + gd%ref
        call glean(gdt)
        call glean(thy(gd2))
      end subroutine

      function gd_ref(gd) result(r)
!doc$ function x_ref(gd) result(r)
        type(gen_density_obj) :: gd
        integer, dimension(2) :: r
!       effects: Returns gd%ref and gd%o%ref.

!cod$
        r(1) = gd%ref
        r(2) = gd%o%ref
        call glean(gd)
      end function

      function gd_ghost(gd) result(g)
!doc$ function x_ghost(gd) result(g)
        type(gen_density_obj) :: gd
        type(ghost) :: g
!       effects: Returns the ghost of gd.
        
!cod$
        g = gd%o%g
        call glean(gd)
      end function

      function gd_grid_density(gd) result(gden)
!doc$ function x_grid_density(gd) result(gden)
        type(gen_density_obj) :: gd
        type(grid_obj) :: gden
!       effects: Returns the gden contained in gd.

!cod$
        call my(gd)
        call my(gd%o%gden,gden)
        call glean(thy(gd))
        call bequeath(thy(gden))
      end function

      function gd_atomic_density(gd) result(aden)
!doc$ function x_atomic_density(gd) result(aden)
        type(gen_density_obj) :: gd
        type(atomic_density_obj) :: aden
!       effects: Returns the aden contained in gd.

!cod$
        call my(gd)
        call my(gd%o%aden,aden)
        call glean(thy(gd))
        call bequeath(thy(aden))
      end function

      function distance_gd(gd1,gd2) result(d)
!doc$ function distance(gd1,gd2) result(d)
        type(gen_density_obj) :: gd1, gd2
        real(double) :: d
!       effects: Returns the rms difference between gd1 and gd2.

!cod$
        real(double) :: da, da_sg, dg, dg_sg

        call my(gd1)
        call my(gd2)

        da_sg = distance(gd1%o%aden,gd2%o%aden)              ; if (error()) goto 100
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,da_sg,da)  ; if (error()) goto 100
        dg_sg = distance(gd1%o%gden,gd2%o%gden)              ; if (error()) goto 100
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,dg_sg,dg)  ; if (error()) goto 100
        d = sqrt(da**2 + dg**2)

100     call glean(thy(gd1))
        call glean(thy(gd2))

        if (error("Exit gen_density_mod::distance_gd")) continue

      end function

      subroutine spin_information_gd(gd,na1,na2,ng1,ng2)
!doc$ subroutine spin_information(gd,na1,na2,ng1,ng2)
        type(gen_density_obj) :: gd
        real(double) :: na1, na2, ng1, ng2
!       requires: mpi_nsgroups() = 2.
!       effects: Returns the numbers of atomic and grid electrons in spin groups 1 and 2.
!       errors: x_scope(gd%o%gden) /= SGROUP.

!cod$
        integer :: msg
        real(double ) :: ne_a, ne_g
        real(double), dimension(2) :: ne_a_c, ne_a_sg, ne_g_c, ne_g_sg
        complex(double) :: gden_norm

        call my(gd)

        if (error(x_scope(gd%o%gden) /= SGROUP,"ERROR: grid density scope is not SGROUP")) goto 100

        msg = mpi_mysgroup()

        ! atomic contributions
        call get_normalization(gd%o%aden,ne_a)
        ne_a_sg = 0.0_double
        ne_a_sg(msg) = ne_a
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,ne_a_sg,ne_a_c) ; if (error()) goto 100
        na1 = ne_a_c(1)
        na2 = ne_a_c(2)

        ! grid contributions
        call get_normalization(gd%o%gden,gden_norm)
        ne_g = real(gden_norm,double)*x_cell_volume(x_lattice(x_layout(gd%o%gden)))
        ne_g_sg = 0.0_double
        ne_g_sg(msg) = ne_g
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,ne_g_sg,ne_g_c) ; if (error()) goto 100

        ng1 = ne_g_c(1)
        ng2 = ne_g_c(2)

100     call glean(thy(gd))

        if (error("Exit gen_density_mod::spin_information_gd")) continue

      end subroutine

      subroutine merge_density_gd(gd,ext)
!doc$ subroutine merge_density(gd,ext)
        type(gen_density_obj) :: gd
        type(external_obj) :: ext
!       requires: gd%o%gden have KGROUP scope.
!       effects: Merges gd.
!       modifies: gd
!       errors: Passes errors.
!       notes: gd%o%gden arrives with KGROUP scope and returns with SGROUP scope.

!cod$
        call my(gd)
        call my(ext)

        call own_i(gd)
        gd%o%g = x_ghost()

        call merge_grid_density(gd%o%gden) ; if (error()) goto 100
        call merge_atomic_density(gd%o%aden)

100     call glean(thy(gd))
        call glean(thy(ext))

        if (error("Exit gen_density_mod::merge_gd")) continue

      end subroutine

      subroutine symmetrize_gd(gd,ext)
!doc$ subroutine symmetrize(gd,ext)
        type(gen_density_obj) :: gd
        type(external_obj) :: ext
!       modifies: gd
!       effects: Symmetrizes gd.

!cod$
        call my(gd)
        call my(ext)

        call own_i(gd)
        gd%o%g = x_ghost()

        call symmetrize_grid(x_space_group(ext),gd%o%gden) ; if (error()) goto 100
        call symmetrize(gd%o%aden,x_space_group(ext))

100     call glean(thy(gd))
        call glean(thy(ext))

      end subroutine

      subroutine filter_gd(gd)
!doc$ subroutine filter(gd)
        type(gen_density_obj) :: gd
!       modifies: gd
!       effects: Filters gd%o%gden.

!cod$
        call my(gd)

        call own_i(gd)
        gd%o%g = x_ghost()

        call filter(gd%o%gden) ; if (error()) goto 100

100     call glean(thy(gd))

      end subroutine

      subroutine merge_symmetrize_filter_gd(gd,ext)
!doc$ subroutine merge_symmetrize_filter(gd,ext)
        type(gen_density_obj) :: gd
        type(external_obj) :: ext
!       requires: gd%o%gden have KGROUP scope.
!       effects: Merges and symmetrizes gd.
!       modifies: gd
!       errors: Passes errors.
!       notes: gd%o%gden arrives with KGROUP scope and returns with SGROUP scope.

!cod$
        call my(gd)
        call my(ext)

        call own_i(gd)
        gd%o%g = x_ghost()

        call merge_grid_density(gd%o%gden) ; if (error()) goto 100
        call symmetrize_grid(x_space_group(ext),gd%o%gden) ; if (error()) goto 100
        call filter(gd%o%gden) ; if (error()) goto 100

        call merge_atomic_density(gd%o%aden) ; if (error()) goto 100
        call symmetrize(gd%o%aden,x_space_group(ext)) ; if (error()) goto 100

100     call glean(thy(gd))
        call glean(thy(ext))

        if (error("Exit gen_density_mod::merge_symmetrize_filter_gd")) continue

      end subroutine

! local routines

      subroutine own_i(gd)
        type(gen_density_obj) :: gd

        type(gen_density_obj) :: gdt
        if (gd%ref < gd%o%ref) then
          allocate( gdt%o )
          gdt%o%ref = 0
          call my(gd%o%gden,gdt%o%gden)
          call my(gd%o%aden,gdt%o%aden)
          gdt%o%g = gd%o%g
          gd%o%ref = gd%o%ref - gd%ref
          gd%o => gdt%o
          gd%o%ref = gd%o%ref + gd%ref
        end if
      end subroutine

      end module

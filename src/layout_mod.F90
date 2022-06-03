!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module layout_mod
!doc$ module layout_mod

!     One datatype is available here: type(layout_obj)

!     layout_obj specifies the distributions of real- and reciprocal-space points.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use diary_mod
      use arg_mod
      use math_mod
      use fft_mod
      use ghost_mod
      use lattice_mod

!cod$
      implicit none
      private

      ! data distribution types
      integer, parameter :: S_TYPE = 1  ! serial
      integer, parameter :: D_TYPE = 2  ! distributed according to scope

      type :: distribution
        integer, dimension(3) :: locdims
        integer, dimension(3) :: bdims
        integer, dimension(3) :: base
        integer, dimension(3) :: block
        type(fft_distributed_plan) :: dplan
        logical, dimension(:,:,:), pointer :: filter
      end type

      type :: layout_rep
        integer :: ref
        type(ghost) :: g
        real(double) :: cutoff
        type(lattice_obj) :: lattice
        integer, dimension(3) :: dims
        integer, dimension(3) :: origin
        integer, dimension(3) :: fft_origin
        type(distribution) :: config
        type(distribution) :: sgroup
        type(distribution) :: kgroup
        type(fft_serial_plan) :: splan
      end type 

      type, public :: layout_obj
        private
        integer :: ref
        type(layout_rep), pointer :: o
      end type

!doc$
      public :: layout
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_cutoff
      public :: x_lattice
      public :: x_dims
      public :: x_origin
      public :: x_fft_origin
      public :: x_serial_plan
      public :: fft_serial
      public :: fft_distributed_sc
      public :: gather
      public :: scatter
      public :: alloc
      public :: consistent
      public :: mesh
      public :: a2rl
      public :: a2rc
      public :: rl2a
      public :: rc2a
      public :: ad2as
      public :: as2ad
      public :: fmesh
      public :: a2fl
      public :: a2fc
      public :: a2f2
      public :: fl2a
      public :: fc2a
      public :: fdel
      public :: fdelinv
      public :: filter
      public :: diary
      public :: write_restart
      public :: S_TYPE
      public :: D_TYPE
      public :: R_TO_Q ! defined in fft_mod
      public :: Q_TO_R ! defined in fft_mod

!cod$
      interface layout
        module procedure constructor_lay
      end interface
      interface my
        module procedure my_lay, my_new_lay
      end interface
      interface thy
        module procedure thy_lay
      end interface
      interface bequeath
        module procedure bequeath_lay
      end interface
      interface glean
        module procedure glean_lay
      end interface
      interface assignment(=)
        module procedure assign_lay
      end interface
      interface x_ref
        module procedure lay_ref
      end interface
      interface x_ghost
        module procedure lay_ghost
      end interface
      interface x_cutoff
        module procedure lay_cutoff
      end interface
      interface x_lattice
        module procedure lay_lattice
      end interface
      interface x_dims
        module procedure lay_dims
      end interface
      interface x_origin
        module procedure lay_origin
      end interface
      interface x_fft_origin
        module procedure lay_fft_origin
      end interface
      interface x_serial_plan
        module procedure lay_serial_plan
      end interface
      interface fft_serial
        module procedure fft_serial_lay
      end interface
      interface fft_distributed_sc
        module procedure fft_distributed_sc_lay
      end interface
      interface alloc
        module procedure alloc_l_lay, alloc_i_lay, alloc_r_lay, alloc_c_lay, alloc_mr_lay, alloc_mc_lay
      end interface
      interface consistent
        module procedure consistent_data_r, consistent_data_c
      end interface
      interface gather
        module procedure gather_i_lay, gather_r_lay, gather_c_lay
      end interface
      interface scatter
        module procedure scatter_i_lay, scatter_r_lay, scatter_c_lay
      end interface
      interface mesh
        module procedure mesh_lay
      end interface
      interface a2rl
        module procedure a2rl_lay
      end interface
      interface a2rc
        module procedure a2rc_lay
      end interface
      interface rl2a
        module procedure rl2a_lay
      end interface
      interface rc2a
        module procedure rc2a_lay
      end interface
      interface ad2as
        module procedure ad2as_lay
      end interface
      interface as2ad
        module procedure as2ad_lay
      end interface
      interface fmesh
        module procedure fmesh_lay
      end interface
      interface a2fl
        module procedure a2fl_lay
      end interface
      interface a2fc
        module procedure a2fc_lay
      end interface
      interface a2f2
        module procedure a2f2_lay
      end interface
      interface fl2a
        module procedure fl2a_lay
      end interface
      interface fc2a
        module procedure fc2a_lay
      end interface
      interface fdel
        module procedure fdel_lay
      end interface
      interface fdelinv
        module procedure fdelinv_lay
      end interface
      interface filter
        module procedure filter_r_lay, filter_c_lay
      end interface
      interface diary
        module procedure diary_lay
      end interface
      interface write_restart
        module procedure write_restart_lay
      end interface

      contains

! public routines

      function constructor_lay(lat,cutoff_in,restf) result(lay)
!doc$ function layout(lat,cutoff_in,restf) result(lay)
        type(lattice_obj) :: lat
        real(double), optional :: cutoff_in
        type(tagio_obj), optional :: restf
        type(layout_obj) :: lay
!       requires: Both cutoff_in and restf not be present. cutoff_in be > 0.
!       effects: Creates a new lay.

!cod$
        logical :: found
        character(1) :: tios
        integer :: ic
        integer, dimension(3) :: dims, n2, n3, n5, n7
        integer(long) :: dsize, iosl, ndata
        real(double) :: recip_space_volume
        real(double), dimension(3) :: b1, b2, b3

        if (error("  Error on entry")) then
          lay%ref = 0
          allocate( lay%o )
          lay%o%ref = 0
          goto 999
        end if

        call my(lat)
        if (present(restf)) call my(restf)

        lay%ref = 0
        allocate( lay%o )
        lay%o%ref = 0
        lay%o%g = x_ghost()

        call my(lat,lay%o%lattice)

        ! open the LAYOUT block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"LAYOUT")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: LAYOUT block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)
        end if

        ! determine the mesh cutoff
        if (present(cutoff_in)) then
          lay%o%cutoff = cutoff_in
        else
          call arg("den_cutoff",lay%o%cutoff,found)
          if (found) then
            if (error(lay%o%cutoff <= 0.0_double,"ERROR: den_cutoff <= 0")) goto 100
          else
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"CUTOFF")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR: CUTOFF tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_double ; ndata = 1
                call readf(lay%o%cutoff,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,lay%o%cutoff)
            else
              if (error(.true.,"ERROR: den_cutoff was not found")) goto 100
            end if
          end if
        end if

        ! compute the mesh dimensions
        b1 = lat2f(lay%o%lattice,real((/1,0,0/),double))
        b2 = lat2f(lay%o%lattice,real((/0,1,0/),double))
        b3 = lat2f(lay%o%lattice,real((/0,0,1/),double))
        recip_space_volume = abs( dot_product(b1,cross_product(b2,b3)) )
        dims(1) = 2*ceiling( sqrt(lay%o%cutoff)*norm(cross_product(b2,b3))/recip_space_volume ) + 1
        dims(2) = 2*ceiling( sqrt(lay%o%cutoff)*norm(cross_product(b3,b1))/recip_space_volume ) + 1
        dims(3) = 2*ceiling( sqrt(lay%o%cutoff)*norm(cross_product(b1,b2))/recip_space_volume ) + 1

        ! set dims so that they can be expressed in powers of 2, 3, 5, and 7
        do ic = 1,3
          call get_nf(dims(ic),n2(ic),n3(ic),n5(ic),n7(ic))
        end do

        ! compute the layout dimensions
        call set_layout_dims_i(lay%o,dims)

        ! form the cutoff filters
        call form_filters_i(lay%o)

        ! form the FFT plans
        call fft_create_serial_plan(lay%o%dims,lay%o%splan)
        call fft_create_distributed_plan(mpi_comm(CONFIG),lay%o%dims,lay%o%config%locdims,lay%o%config%base,lay%o%config%dplan)
        call fft_create_distributed_plan(mpi_comm(SGROUP),lay%o%dims,lay%o%sgroup%locdims,lay%o%sgroup%base,lay%o%sgroup%dplan)
        call fft_create_distributed_plan(mpi_comm(KGROUP),lay%o%dims,lay%o%kgroup%locdims,lay%o%kgroup%base,lay%o%kgroup%dplan)

100     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

200     call glean(thy(lat))
        if (present(restf)) call glean(thy(restf))

999     if (error("Exit layout_mod::constructor_lay")) continue

      end function

      subroutine my_lay(lay)
!doc$ subroutine my(lay)
        type(layout_obj) :: lay

!cod$
        lay%ref = lay%ref + 1
        lay%o%ref = lay%o%ref + 1
      end subroutine

      subroutine my_new_lay(layi,lay)
!doc$ subroutine my(layi,lay)
        type(layout_obj) :: layi, lay

!cod$
        lay%ref = 1
        lay%o => layi%o
        lay%o%ref = lay%o%ref + 1
      end subroutine

      function thy_lay(lay) result(layo)
!doc$ function thy(lay) result(layo)
        type(layout_obj) :: lay, layo

!cod$
        lay%ref = lay%ref - 1
        lay%o%ref = lay%o%ref - 1
        layo%ref = lay%ref
        layo%o => lay%o
      end function

      subroutine glean_lay(lay)
!doc$ subroutine glean(lay)
        type(layout_obj) :: lay

!cod$
        if (lay%o%ref < 1) then
          call glean(thy(lay%o%lattice))
          if (associated( lay%o%config%filter )) deallocate( lay%o%config%filter )
          if (associated( lay%o%sgroup%filter )) deallocate( lay%o%sgroup%filter )
          if (associated( lay%o%kgroup%filter )) deallocate( lay%o%kgroup%filter )
          call fft_destroy_serial_plan(lay%o%splan)
          call fft_destroy_distributed_plan(lay%o%config%dplan)
          call fft_destroy_distributed_plan(lay%o%sgroup%dplan)
          call fft_destroy_distributed_plan(lay%o%kgroup%dplan)
          deallocate( lay%o )
        end if
      end subroutine
      
      subroutine bequeath_lay(lay)
!doc$ subroutine bequeath(lay)
        type(layout_obj) :: lay

!cod$
        continue
      end subroutine

      subroutine assign_lay(lay,lay2)
!doc$ subroutine assign(lay,lay2)
        type(layout_obj), intent(inout) :: lay
        type(layout_obj), intent(in) :: lay2

!cod$
        type(layout_obj) :: layt
        call my(lay2)
        layt%o => lay%o
        lay%o%ref = lay%o%ref - lay%ref
        lay%o => lay2%o
        lay%o%ref = lay%o%ref + lay%ref
        call glean(layt)
        call glean(thy(lay2))
      end subroutine

      function lay_ref(lay) result(r)
!doc$ function x_ref(lay) result(r)
        type(layout_obj) :: lay
        integer, dimension(2) :: r
!       effects: Returns lay%ref and lay%o%ref.

!cod$
        r(1) = lay%ref
        r(2) = lay%o%ref
      end function

      function lay_ghost(lay) result(g)
!doc$ function x_ghost(lay) result(g)
        type(layout_obj) :: lay
        type(ghost) :: g
!       effects: Returns lay%o%g.

!cod$
        call my(lay)
        g = lay%o%g
        call glean(thy(lay))
      end function

      function lay_cutoff(lay) result(c)
!doc$ function x_cutoff(lay) result(c)
        type(layout_obj) :: lay
        real(double) :: c
!       effects: Returns lay%o%cutoff.

!cod$
        call my(lay)
        c = lay%o%cutoff
        call glean(thy(lay))
      end function

      function lay_lattice(lay) result(lat)
!doc$ function x_lattice(lay) result(lat)
        type(layout_obj) :: lay
        type(lattice_obj) :: lat
!       effects: Returns lay%o%lattice.

!cod$
        call my(lay)
        call my(lay%o%lattice,lat)
        call bequeath(thy(lat))
        call glean(thy(lay))
      end function

      function lay_dims(lay) result(dims)
!doc$ function x_dims(lay) result(dims)
        type(layout_obj) :: lay
        integer, dimension(3) :: dims
!       effects: Returns lay%o%dims.

!cod$
        call my(lay)
        dims = lay%o%dims
        call glean(thy(lay))
      end function

      function lay_origin(lay) result(origin)
!doc$ function x_origin(lay) result(origin)
        type(layout_obj) :: lay
        integer, dimension(3) :: origin
!       effects: Returns lay%o%origin.

!cod$
        call my(lay)
        origin = lay%o%origin
        call glean(thy(lay))
      end function

      function lay_fft_origin(lay,tp,sc) result(fft_origin)
!doc$ function x_fft_origin(lay,tp,sc) result(fft_origin)
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, dimension(3) :: fft_origin
        integer, intent(in), optional :: sc
!       requires: tp be S_TYPE or D_TYPE
!                 sc be present if tp = D_TYPE
!       effects: Returns the coordinates of the FFT origin if it is within the processor scope and zero otherwise.

!cod$
        call my(lay)
        select case (tp)
        case (S_TYPE)
          fft_origin = lay%o%fft_origin
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            fft_origin = lay%o%fft_origin - (lay%o%config%base - (/1,1,1/))
            if (any(fft_origin < 1) .or. any(fft_origin > lay%o%config%locdims)) fft_origin = 0
          case (SGROUP)
            fft_origin = lay%o%fft_origin - (lay%o%sgroup%base - (/1,1,1/))
            if (any(fft_origin < 1) .or. any(fft_origin > lay%o%sgroup%locdims)) fft_origin = 0
          case (KGROUP)
            fft_origin = lay%o%fft_origin - (lay%o%kgroup%base - (/1,1,1/))
            if (any(fft_origin < 1) .or. any(fft_origin > lay%o%kgroup%locdims)) fft_origin = 0
          end select
        end select
        call glean(thy(lay))
      end function

      function lay_serial_plan(lay) result(plan)
!doc$ function x_serial_plan(lay) result(plan)
        type(layout_obj) :: lay
        type(fft_serial_plan) :: plan
!       effects: Returns lay%o%splan.

!cod$
        call my(lay)
        plan = lay%o%splan
        call glean(thy(lay))
      end function

      subroutine alloc_l_lay(data,lay,tp,sc)
!doc$ subroutine alloc(data,lay,tp,sc)
        logical, dimension(:,:,:), pointer :: data
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: data be nullified.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE
!       modifies: data
!       effects: Allocates data.

!cod$
        call my(lay)
        select case(tp)
        case (S_TYPE)
          allocate( data(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            allocate( data(lay%o%config%locdims(1),lay%o%config%locdims(2),lay%o%config%locdims(3)) )
          case (SGROUP)
            allocate( data(lay%o%sgroup%locdims(1),lay%o%sgroup%locdims(2),lay%o%sgroup%locdims(3)) )
          case (KGROUP)
            allocate( data(lay%o%kgroup%locdims(1),lay%o%kgroup%locdims(2),lay%o%kgroup%locdims(3)) )
          end select
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::alloc_l_lay")) continue
      end subroutine

      subroutine alloc_i_lay(data,lay,tp,sc)
!doc$ subroutine alloc(data,lay,tp,sc)
        integer, dimension(:,:,:), pointer :: data
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: data be nullified.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       modifies: data
!       effects: Allocates data.

!cod$
        call my(lay)
        select case(tp)
        case (S_TYPE)
          allocate( data(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            allocate( data(lay%o%config%locdims(1),lay%o%config%locdims(2),lay%o%config%locdims(3)) )
          case (SGROUP)
            allocate( data(lay%o%sgroup%locdims(1),lay%o%sgroup%locdims(2),lay%o%sgroup%locdims(3)) )
          case (KGROUP)
            allocate( data(lay%o%kgroup%locdims(1),lay%o%kgroup%locdims(2),lay%o%kgroup%locdims(3)) )
          end select
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::alloc_i_lay")) continue
      end subroutine

      subroutine alloc_r_lay(data,lay,tp,sc)
!doc$ subroutine alloc(data,lay,tp,sc)
        real(double), dimension(:,:,:), pointer :: data
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: data be nullified.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       modifies: data
!       effects: Allocates data.

!cod$
        call my(lay)
        select case(tp)
        case (S_TYPE)
          allocate( data(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        case (D_TYPE)
          if (error(.not.present(sc),'Error: Must include scope in argument for D_TYPE alloc')) goto 100
          select case (sc)
          case (CONFIG)
            allocate( data(lay%o%config%locdims(1),lay%o%config%locdims(2),lay%o%config%locdims(3)) )
          case (SGROUP)
            allocate( data(lay%o%sgroup%locdims(1),lay%o%sgroup%locdims(2),lay%o%sgroup%locdims(3)) )
          case (KGROUP)
            allocate( data(lay%o%kgroup%locdims(1),lay%o%kgroup%locdims(2),lay%o%kgroup%locdims(3)) )
          end select
        end select
        call glean(thy(lay))
100     if (error("Exit layout_mod::alloc_r_lay")) continue
      end subroutine

      subroutine alloc_c_lay(data,lay,tp,sc)
!doc$ subroutine alloc(data,lay,tp,sc)
        complex(double), dimension(:,:,:), pointer :: data
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: data be nullified.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       modifies: data
!       effects: Allocates data.

!cod$
        call my(lay)
        select case(tp)
        case (S_TYPE)
          allocate( data(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            allocate( data(lay%o%config%locdims(1),lay%o%config%locdims(2),lay%o%config%locdims(3)) )
          case (SGROUP)
            allocate( data(lay%o%sgroup%locdims(1),lay%o%sgroup%locdims(2),lay%o%sgroup%locdims(3)) )
          case (KGROUP)
            allocate( data(lay%o%kgroup%locdims(1),lay%o%kgroup%locdims(2),lay%o%kgroup%locdims(3)) )
          end select
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::alloc_c_lay")) continue
      end subroutine

      subroutine alloc_mr_lay(data,lay,m,tp,sc)
!doc$ subroutine alloc(data,lay,m,tp,sc)
        real(double), dimension(:,:,:,:), pointer :: data
        type(layout_obj) :: lay
        integer, intent(in) :: m, tp
        integer, intent(in), optional :: sc
!       requires: data be nullified.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       modifies: data
!       effects: Allocates data.

!cod$
        call my(lay)
        select case(tp)
        case (S_TYPE)
          allocate( data(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3),m) )
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            allocate( data(lay%o%config%locdims(1),lay%o%config%locdims(2),lay%o%config%locdims(3),m) )
          case (SGROUP)
            allocate( data(lay%o%sgroup%locdims(1),lay%o%sgroup%locdims(2),lay%o%sgroup%locdims(3),m) )
          case (KGROUP)
            allocate( data(lay%o%kgroup%locdims(1),lay%o%kgroup%locdims(2),lay%o%kgroup%locdims(3),m) )
          end select
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::alloc_mr_lay")) continue
      end subroutine

      subroutine alloc_mc_lay(data,lay,m,tp,sc)
!doc$ subroutine alloc(data,lay,m,tp,sc)
        complex(double), dimension(:,:,:,:), pointer :: data
        type(layout_obj) :: lay
        integer, intent(in) :: m, tp
        integer, intent(in), optional :: sc
!       requires: data be nullified.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       modifies: data
!       effects: Allocates data.

!cod$
        call my(lay)
        select case(tp)
        case (S_TYPE)
          allocate( data(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3),m) )
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            allocate( data(lay%o%config%locdims(1),lay%o%config%locdims(2),lay%o%config%locdims(3),m) )
          case (SGROUP)
            allocate( data(lay%o%sgroup%locdims(1),lay%o%sgroup%locdims(2),lay%o%sgroup%locdims(3),m) )
          case (KGROUP)
            allocate( data(lay%o%kgroup%locdims(1),lay%o%kgroup%locdims(2),lay%o%kgroup%locdims(3),m) )
          end select
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::alloc_mc_lay")) continue
      end subroutine

      function consistent_data_r(data,lay,tp,sc) result(c)
!doc$ function consistent(data,lay,tp,sc) result(c)
        real(double), dimension(:,:,:) :: data
        type(layout_obj) :: lay
        integer, intent(in) :: sc, tp
        logical :: c
!       effects: Checks that shape(data) is consistent with tp and sc.

!cod$
        call my(lay)
        c = .false.
        select case (tp)
        case (S_TYPE)
          c = all(lay%o%dims == shape(data))
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            c = all(lay%o%config%locdims == shape(data))
          case (SGROUP)
            c = all(lay%o%sgroup%locdims == shape(data))
          case (KGROUP)
            c = all(lay%o%kgroup%locdims == shape(data))
          end select
        end select
        call glean(thy(lay))
      end function

      function consistent_data_c(data,lay,tp,sc) result(c)
!doc$ function consistent(data,lay,tp,sc) result(c)
        complex(double), dimension(:,:,:) :: data
        type(layout_obj) :: lay
        integer, intent(in) :: sc, tp
        logical :: c
!       effects: Checks that shape(data) is consistent with tp and sc.

!cod$
        call my(lay)
        c = .false.
        select case(tp)
        case (S_TYPE)
          c = all(lay%o%dims == shape(data))
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            c = all(lay%o%config%locdims == shape(data))
          case (SGROUP)
            c = all(lay%o%sgroup%locdims == shape(data))
          case (KGROUP)
            c = all(lay%o%kgroup%locdims == shape(data))
          end select
        end select
        call glean(thy(lay))
      end function

      subroutine gather_i_lay(lay,sc,idata,odata)
!doc$ subroutine gather(lay,sc,idata,odata)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        integer, dimension(:,:,:), intent(in) :: idata
        integer, dimension(:,:,:), intent(out) :: odata
!       requires: idata has scope sc. shape(odata) = x_dims(lay)
!       modifies: odata
!       effects: Constructs a S_TYPE representation of idata in odata.
!       errors: Passes errors.

!cod$
        integer, dimension(3) :: l, u
        integer, dimension(:,:,:), allocatable :: tmp
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = lay%o%config%base + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = lay%o%sgroup%base + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = lay%o%kgroup%base + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        allocate( tmp(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        tmp = 0
        tmp(l(1):u(1),l(2):u(2),l(3):u(3)) = idata
        call allreduce(sc,MPI_SUM,tmp,odata)
        if (allocated( tmp )) deallocate( tmp )
        call glean(thy(lay))
        if (error("Exit layout_mod::gather_i_lay")) continue
      end subroutine

      subroutine gather_r_lay(lay,sc,idata,odata)
!doc$ subroutine gather(lay,sc,idata,odata)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        real(double), dimension(:,:,:), intent(in) :: idata
        real(double), dimension(:,:,:), intent(out) :: odata
!       requires: idata has scope sc. shape(odata) = x_dims(lay)
!       modifies: odata
!       effects: Constructs a S_TYPE representation of idata in odata.
!       errors: Passes errors.

!cod$
        integer, dimension(3) :: l, u
        real(double), dimension(:,:,:), allocatable :: tmp
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = lay%o%config%base + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = lay%o%sgroup%base + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = lay%o%kgroup%base + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        allocate( tmp(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        tmp = 0.0_double
        tmp(l(1):u(1),l(2):u(2),l(3):u(3)) = idata
        call allreduce(sc,MPI_SUM,tmp,odata)
        if (allocated( tmp )) deallocate( tmp )
        call glean(thy(lay))
        if (error("Exit layout_mod::gather_r_lay")) continue
      end subroutine

      subroutine gather_c_lay(lay,sc,idata,odata)
!doc$ subroutine gather(lay,sc,idata,odata)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        complex(double), dimension(:,:,:), intent(in) :: idata
        complex(double), dimension(:,:,:), intent(out) :: odata
!       requires: idata has scope sc. shape(odata) = x_dims(lay)
!       modifies: odata
!       effects: Constructs a S_TYPE representation of idata in odata.
!       errors: Passes errors.

!cod$
        integer, dimension(3) :: l, u
        complex(double), dimension(:,:,:), allocatable :: tmp
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = lay%o%config%base + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = lay%o%sgroup%base + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = lay%o%kgroup%base + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        allocate( tmp(lay%o%dims(1),lay%o%dims(2),lay%o%dims(3)) )
        tmp = (0.0_double,0.0_double)
        tmp(l(1):u(1),l(2):u(2),l(3):u(3)) = idata
        call allreduce(sc,MPI_SUM,tmp,odata)
        if (allocated( tmp )) deallocate( tmp )
        call glean(thy(lay))
        if (error("Exit layout_mod::gather_c_lay")) continue
      end subroutine

      subroutine scatter_i_lay(lay,sc,idata,odata)
!doc$ subroutine scatter(lay,sc,idata,odata)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        integer, dimension(:,:,:), intent(in) :: idata
        integer, dimension(:,:,:), intent(out) :: odata
!       requires: odata has scope sc. shape(idata) = x_dims(lay)
!       modifies: odata
!       effects: Constructs a D_TYPE representation of idata in odata.

!cod$
        integer, dimension(3) :: l, u
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = lay%o%config%base + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = lay%o%sgroup%base + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = lay%o%kgroup%base + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        odata = idata(l(1):u(1),l(2):u(2),l(3):u(3))
        call glean(thy(lay))
      end subroutine

      subroutine scatter_r_lay(lay,sc,idata,odata)
!doc$ subroutine scatter(lay,sc,idata,odata)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        real(double), dimension(:,:,:), intent(in) :: idata
        real(double), dimension(:,:,:), intent(out) :: odata
!       requires: odata has scope sc. shape(idata) = x_dims(lay)
!       modifies: odata
!       effects: Constructs a D_TYPE representation of idata in odata.

!cod$
        integer, dimension(3) :: l, u
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = lay%o%config%base + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = lay%o%sgroup%base + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = lay%o%kgroup%base + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        odata = idata(l(1):u(1),l(2):u(2),l(3):u(3))
        call glean(thy(lay))
      end subroutine

      subroutine scatter_c_lay(lay,sc,idata,odata)
!doc$ subroutine scatter(lay,sc,idata,odata)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        complex(double), dimension(:,:,:), intent(in) :: idata
        complex(double), dimension(:,:,:), intent(out) :: odata
!       requires: odata has scope sc. shape(idata) = x_dims(lay)
!       modifies: odata
!       effects: Constructs a D_TYPE representation of idata in odata.

!cod$
        integer, dimension(3) :: l, u
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = lay%o%config%base + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = lay%o%sgroup%base + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = lay%o%kgroup%base + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        odata = idata(l(1):u(1),l(2):u(2),l(3):u(3))
        call glean(thy(lay))
      end subroutine

      subroutine fft_serial_lay(lay,data,dir)
!doc$ subroutine fft_serial(lay,data,dir)
        type(layout_obj) :: lay
        complex(double), dimension(:,:,:), intent(inout) :: data
        integer, intent(in) :: dir
!       requires: shape(data) = x_dims(lay). dir be R_TO_Q or Q_TO_R.
!       modifies: data
!       effects: data is overwritten with its serial Fourier transform.
!       usage: Meant to be called from grid_mod.
!       errors: Passes errors.

!cod$
        complex(double), dimension(:,:,:), pointer :: c1
        call my(lay)
        c1 => datalink(lay%o%splan)
        c1 = data
        call fft_serial(dir,lay%o%splan)
        data = c1
        nullify( c1 )
        call glean(thy(lay))
        if (error("Exit layout_mod::fft_serial_lay")) continue
      end subroutine

      subroutine fft_distributed_sc_lay(lay,sc,data,dir)
!doc$ subroutine fft_distributed_sc(lay,sc,data,dir)
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        complex(double), dimension(:,:,:), intent(inout) :: data
        integer, intent(in) :: dir
!       requires: data have scope sc. dir be R_TO_Q or Q_TO_R.
!       modifies: data
!       effects: data is overwritten with its Fourier transform.
!       errors: Passes errors.

!cod$
        call my(lay)
        select case (sc)
        case (CONFIG)
          call fft_distributed(data,dir,lay%o%config%dplan)
        case (SGROUP)
          call fft_distributed(data,dir,lay%o%sgroup%dplan)
        case (KGROUP)
          call fft_distributed(data,dir,lay%o%kgroup%dplan)
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::fft_distributed_sc_lay")) continue
      end subroutine

      subroutine mesh_lay(x,y,z,lay,tp,sc)
!doc$ subroutine mesh(x,y,z,lay,tp,sc)
        real(double), dimension(:,:,:), pointer :: x, y, z
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!                 x, y, z be nullified.
!       modifies: x, y, z
!       effects: Allocates x, y, z and returns real-space mesh values (cartesian representation).

!cod$
        integer :: i1, i2, i3
        integer, dimension(3) :: d, i, o, s
        real(double), dimension(3) :: rc, rl, rd

        if (error(associated(x),"Error: x data associated, must be nullified")) goto 100
        if (error(associated(y),"Error: y data associated, must be nullified")) goto 100
        if (error(associated(z),"Error: z data associated, must be nullified")) goto 100
        call my(lay)
        select case (tp)
        case (S_TYPE)
          d = lay%o%dims
          s = (/0,0,0/)
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            d = lay%o%config%locdims
            s = lay%o%config%base - (/1,1,1/)
          case (SGROUP)
            d = lay%o%sgroup%locdims
            s = lay%o%sgroup%base - (/1,1,1/)
          case (KGROUP)
            d = lay%o%kgroup%locdims
            s = lay%o%kgroup%base - (/1,1,1/)
          end select
        end select
        o = lay%o%origin
        rd = real(lay%o%dims,double)
        allocate( x(d(1),d(2),d(3)) )
        allocate( y(d(1),d(2),d(3)) )
        allocate( z(d(1),d(2),d(3)) )
        do i3 = 1,d(3)
        do i2 = 1,d(2)
        do i1 = 1,d(1)
          i = (/i1,i2,i3/) + s - o
          rl = real(i,double)/rd
          rc = lat2r(lay%o%lattice,rl)
          x(i1,i2,i3) = rc(1)
          y(i1,i2,i3) = rc(2)
          z(i1,i2,i3) = rc(3)
        end do
        end do
        end do
        call glean(thy(lay))
100     if (error("Exit layout_mod::mesh_lay")) continue
      end subroutine

      function a2rl_lay(a,lay,tp,sc) result(rl)
!doc$ function a2rl(a,lay,tp,sc) result(rl)
        integer, dimension(3), intent(in) :: a
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        real(double), dimension(3) :: rl
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       requires: (/1,1,1/) <= a <= tp mesh limits
!       effects: Returns real-space mesh values (lattice representation) corresponding to array indices a.
!       errors: none

!cod$
        integer, dimension(3) :: i, o, s
        real(double), dimension(3) :: rd
        call my(lay)
        select case (tp)
        case (S_TYPE)
          s = (/0,0,0/)
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            s = lay%o%config%base - (/1,1,1/)
          case (SGROUP)
            s = lay%o%sgroup%base - (/1,1,1/)
          case (KGROUP)
            s = lay%o%kgroup%base - (/1,1,1/)
          end select
        end select
        rd = real(lay%o%dims,double)
        o = lay%o%origin
        i = a + s - o
        rl = real(i,double)/rd
        call glean(thy(lay))
        if (error("Exit layout_mod::a2rl_lay")) continue
      end function

      function a2rc_lay(a,lay,tp,sc) result(rc)
!doc$ function a2rc(a,lay,tp,sc) result(rc)
        integer, dimension(3), intent(in) :: a
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        real(double), dimension(3) :: rc
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       requires: (/1,1,1/) <= a <= tp mesh limits
!       effects: Returns the real-space mesh values (cartesian representation) corresponding to array indices a.
!       errors: none

!cod$
        real(double), dimension(3) :: rl
        call my(lay)
        if (present(sc)) then
          rl = a2rl(a,lay,tp,sc)
        else
          rl = a2rl(a,lay,tp)
        end if
        rc = lat2r(lay%o%lattice,rl)
        call glean(thy(lay))
        if (error("Exit layout_mod::a2rc_lay")) continue
      end function

      function rl2a_lay(rl,lay,tp,sc) result(a)
!doc$ function rl2a(rl,lay,tp,sc) result(a)
        real(double), dimension(3), intent(in) :: rl
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        integer, dimension(3) :: a
!       requires: rl be a mesh value in the lattice representation.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       effects: Returns the array indices corresponding to real-space mesh value rl (lattice representation).
!       errors: none

!cod$
        integer, dimension(3) :: d, i, o, s
        real(double), dimension(3) :: rd
        call my(lay)
        select case (tp)
        case (S_TYPE)
          s = (/0,0,0/)
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            s = lay%o%config%base - (/1,1,1/)
          case (SGROUP)
            s = lay%o%sgroup%base - (/1,1,1/)
          case (KGROUP)
            s = lay%o%kgroup%base - (/1,1,1/)
          end select
        end select
        o = lay%o%origin
        d = lay%o%dims
        rd = real(d,double)
        i = nint(rl*rd) - s + o
        a = modulo(i,d)
        call glean(thy(lay))
        if (error("Exit layout_mod::rl2a_lay")) continue
      end function

      function rc2a_lay(rc,lay,tp,sc) result(a)
!doc$ function rc2a(rc,lay,tp,sc) result(a)
        real(double), dimension(3), intent(in) :: rc
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        integer, dimension(3) :: a
!       requires: rc be a mesh value in the cartesian representation.
!                 tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       effects: Returns the array indices corresponding to real-space mesh value rc (cartesian representation).

!cod$
        real(double), dimension(3) :: rl
        call my(lay)
        rl = r2lat(lay%o%lattice,rc)
        if (present(sc)) then
          a = rl2a(rl,lay,tp,sc)
        else
          a = rl2a(rl,lay,tp)
        end if
        call glean(thy(lay))
        if (error("Exit layout_mod::rc2a_lay")) continue
      end function

      function ad2as_lay(ad,lay,sc) result(as)
!doc$ function ad2as(ad,lay,sc) result(as)
        integer, dimension(3), intent(in) :: ad
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        integer, dimension(3) :: as
!       requires: ad be array indices in a D_TYPE layout.
!       effects: Returns S_TYPE array indices corresponding to D_TYPE array indices.

!cod$
        integer, dimension(3) :: s
        call my(lay)
        select case (sc)
        case (CONFIG)
          s = lay%o%config%base - (/1,1,1/)
        case (SGROUP)
          s = lay%o%sgroup%base - (/1,1,1/)
        case (KGROUP)
          s = lay%o%kgroup%base - (/1,1,1/)
        end select
        as = ad + s
        call glean(thy(lay))
      end function

      function as2ad_lay(as,lay,sc) result(ad)
!doc$ function as2ad(as,lay,sc) result(ad)
        integer, dimension(3), intent(in) :: as
        type(layout_obj) :: lay
        integer, intent(in) :: sc
        integer, dimension(3) :: ad
!       requires: as be array indices in a S_TYPE layout.
!       effects: Returns D_TYPE array indices corresponding to S-TYPE array indices as
!                 in the D-TYPE block and (/0,0,0/) if otherwise.

!cod$
        integer, dimension(3) :: l, u
        call my(lay)
        select case (sc)
        case (CONFIG)
          l = lay%o%config%base
          u = l + lay%o%config%locdims - (/1,1,1/)
        case (SGROUP)
          l = lay%o%sgroup%base
          u = l + lay%o%sgroup%locdims - (/1,1,1/)
        case (KGROUP)
          l = lay%o%kgroup%base
          u = l + lay%o%kgroup%locdims - (/1,1,1/)
        end select
        if (all(as >= l) .and. all(as <= u)) then
          ad = as - l + (/1,1,1/)
        else
          ad = (/0,0,0/)
        end if
        call glean(thy(lay))
      end function

      subroutine fmesh_lay(fx,fy,fz,lay,tp,sc)
!doc$ subroutine fmesh(fx,fy,fz,lay,tp,sc)
        real(double), dimension(:,:,:), pointer :: fx, fy, fz
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!                 fx, fy, fz be nullified.
!       modifies: fx, fy, fz
!       effects: Allocates fx, fy, fz and returns Fourier-space  mesh values (cartesian representation).

!cod$
        integer :: i1, i2, i3, ic
        integer, dimension(3) :: d, fl, i, m, n, o, s
        real(double), dimension(3) :: fc, rfl
        type(lattice_obj) :: lat
        call my(lay)
        call my(x_lattice(lay),lat)
        select case (tp)
        case (S_TYPE)
          d = lay%o%dims
          s = (/0,0,0/)
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            d = lay%o%config%locdims
            s = lay%o%config%base - (/1,1,1/)
          case (SGROUP)
            d = lay%o%sgroup%locdims
            s = lay%o%sgroup%base - (/1,1,1/)
          case (KGROUP)
            d = lay%o%kgroup%locdims
            s = lay%o%kgroup%base - (/1,1,1/)
          end select
        end select
        o = lay%o%fft_origin
        m = lay%o%dims
        do ic = 1,3
          if (odd(m(ic))) then
            n(ic) = (m(ic) - 1)/2
          else
            n(ic) = (m(ic) - 2)/2
          end if
        end do
        allocate( fx(d(1),d(2),d(3)) )
        allocate( fy(d(1),d(2),d(3)) )
        allocate( fz(d(1),d(2),d(3)) )
        do i3 = 1,d(3)
        do i2 = 1,d(2)
        do i1 = 1,d(1)
          i = (/i1,i2,i3/) + s - o
          fl = modulo((i + n),m) - n
          rfl = real(fl,double)
          fc = lat2f(lat,rfl)
          fx(i1,i2,i3) = fc(1)
          fy(i1,i2,i3) = fc(2)
          fz(i1,i2,i3) = fc(3)
        end do
        end do
        end do
        call glean(thy(lat))
        call glean(thy(lay))
        if (error("Exit layout_mod::fmesh_lay")) continue
      end subroutine

      function a2fl_lay(a,lay,tp,sc) result(fl)
!doc$ function a2fl(a,lay,tp,sc) result(fl)
        integer, dimension(3), intent(in) :: a
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        integer, dimension(3) :: fl
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       requires: (/1,1,1/) <= a <= tp mesh limits
!       effects: Returns the Fourier-space mesh value (lattice representation) corresponding to array index a.
!       errors: none

!cod$
        integer :: ic
        integer, dimension(3) :: i, m, n, o, s
        call my(lay)
        select case (tp)
        case (S_TYPE)
          s = (/0,0,0/)
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            s = lay%o%config%base - (/1,1,1/)
          case (SGROUP)
            s = lay%o%sgroup%base - (/1,1,1/)
          case (KGROUP)
            s = lay%o%kgroup%base - (/1,1,1/)
          end select
        end select
        o = lay%o%fft_origin
        m = lay%o%dims
        do ic = 1,3
          if (odd(m(ic))) then
            n(ic) = (m(ic) - 1)/2
          else
            n(ic) = (m(ic) - 2)/2
          end if
        end do
        i = a + s - o
        fl = modulo((i + n),m) - n
        call glean(thy(lay))
        if (error("Exit layout_mod::a2fl_lay")) continue
      end function

      function a2fc_lay(a,lay,tp,sc) result(fc)
!doc$ function a2fc(a,lay,tp,sc) result(fc)
        integer, dimension(3), intent(in) :: a
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        real(double), dimension(3) :: fc
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       effects: Returns the Fourier-space mesh value (cartesian representation) corresponding to array index a.
!       errors: none

!cod$
        integer, dimension(3) :: fl
        real(double), dimension(3) :: rfl
        call my(lay)
        if (present(sc)) then
          fl = a2fl(a,lay,tp,sc)
        else
          fl = a2fl(a,lay,tp)
        end if
        rfl = real(fl,double)
        fc = lat2f(lay%o%lattice,rfl)
        call glean(thy(lay))
        if (error("Exit layout_mod::a2fc_lay")) continue
      end function

      function a2f2_lay(a,lay,tp,sc) result(f2)
!doc$ function a2f2(a,lay,tp,sc) result(f2)
        integer, dimension(3), intent(in) :: a
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        real(double) :: f2
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       effects: Returns the squared Fourier-space mesh value corresponding to array index a.
!       errors: none

!cod$
        real(double), dimension(3) :: fc
        call my(lay)
        if (present(sc)) then
          fc = a2fc(a,lay,tp,sc)
        else
          fc = a2fc(a,lay,tp)
        end if
        f2 = fc(1)**2 + fc(2)**2 + fc(3)**2
        call glean(thy(lay))
        if (error("Exit layout_mod::a2f2_lay")) continue
      end function

      function fl2a_lay(fl,lay,tp,sc) result(a)
!doc$ function fl2a(fl,lay,tp,sc) result(a)
        integer, dimension(3), intent(in) :: fl
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        integer, dimension(3) :: a
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       effects: Returns the array index corresponding to Fourier-space mesh value fl (lattice representation).

!cod$
        integer, dimension(3) :: d, o, s
        call my(lay)
        select case (tp)
        case (S_TYPE)
          s = (/0,0,0/)
        case (D_TYPE)
          select case (sc)
          case (CONFIG)
            s = lay%o%config%base - (/1,1,1/)
          case (SGROUP)
            s = lay%o%sgroup%base - (/1,1,1/)
          case (KGROUP)
            s = lay%o%kgroup%base - (/1,1,1/)
          end select
        end select
        d = lay%o%dims
        o = lay%o%fft_origin
        a = modulo(fl,d) + o - s
        call glean(thy(lay))
        if (error("Exit layout_mod::fl2a_lay")) continue
      end function

      function fc2a_lay(fc,lay,tp,sc) result(a)
!doc$ function fc2a(fc,lay,tp,sc) result(a)
        real(double), dimension(3), intent(in) :: fc
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
        integer, dimension(3) :: a
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!       effects: Returns the array index corresponding to Fourier-space mesh value fc (cartesian representation).

!cod$
        integer, dimension(3) :: fl
        real(double), dimension(3) :: rfl
        call my(lay)
        rfl = f2lat(lay%o%lattice,fc)
        fl = nint(rfl)
        if (present(sc)) then
          a = fl2a(fl,lay,tp,sc)
        else
          a = fl2a(fl,lay,tp)
        end if
        call glean(thy(lay))
        if (error("Exit layout_mod::fc2a_lay")) continue
      end function

      subroutine fdel_lay(del,lay,tp,sc)
!doc$ subroutine fdel(del,lay,tp,sc)
        real(double), dimension(:,:,:), pointer :: del
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!                 del be nullified.
!       modifies: del
!       effects: Allocates del and returns the squared Fourier mesh values in it.

!cod$
        real(double), dimension(:,:,:), pointer :: x, y, z
        call my(lay)
        if (present(sc)) then
          call fmesh(x,y,z,lay,tp,sc)
          call alloc(del,lay,tp,sc)
        else
          call fmesh(x,y,z,lay,tp)
          call alloc(del,lay,tp)
        end if
        del = x**2 + y**2 + z**2
        deallocate( x, y, z )
        call glean(thy(lay))
        if (error("Exit layout_mod::fdel_lay")) continue
      end subroutine

      subroutine fdelinv_lay(deli,lay,tp,sc)
!doc$ subroutine fdelinv(deli,lay,tp,sc)
        real(double), dimension(:,:,:), pointer :: deli
        type(layout_obj) :: lay
        integer, intent(in) :: tp
        integer, intent(in), optional :: sc
!       requires: tp be S_TYPE or D_TYPE.
!                 sc be present if tp = D_TYPE.
!                 deli be nullified.
!       modifies: deli
!       effects: Allocates deli and returns the inverse of the squared Fourier mesh values in it.
!                Returns 0 where the mesh value is 0.

!cod$
        integer, dimension(3) :: ffto
        call my(lay)
        select case (tp)
        case (S_TYPE)
          call fdel(deli,lay,tp)
          ffto = lay%o%fft_origin
          deli(ffto(1),ffto(2),ffto(3)) = 1.0_double
          deli = 1.0_double/deli
          deli(ffto(1),ffto(2),ffto(3)) = 0.0_double
        case (D_TYPE)
          call fdel(deli,lay,tp,sc)
          select case (sc)
          case (CONFIG)
            ffto = lay%o%fft_origin - (lay%o%config%base - (/1,1,1/))
            if (all(ffto > 0) .and. all(ffto <= lay%o%config%locdims)) deli(ffto(1),ffto(2),ffto(3)) = 1.0_double
            deli = 1.0_double/deli
            if (all(ffto > 0) .and. all(ffto <= lay%o%config%locdims)) deli(ffto(1),ffto(2),ffto(3)) = 0.0_double
          case (SGROUP)
            ffto = lay%o%fft_origin - (lay%o%sgroup%base - (/1,1,1/))
            if (all(ffto > 0) .and. all(ffto <= lay%o%sgroup%locdims)) deli(ffto(1),ffto(2),ffto(3)) = 1.0_double
            deli = 1.0_double/deli
            if (all(ffto > 0) .and. all(ffto <= lay%o%sgroup%locdims)) deli(ffto(1),ffto(2),ffto(3)) = 0.0_double
          case (KGROUP)
            ffto = lay%o%fft_origin - (lay%o%kgroup%base - (/1,1,1/))
            if (all(ffto > 0) .and. all(ffto <= lay%o%kgroup%locdims)) deli(ffto(1),ffto(2),ffto(3)) = 1.0_double
            deli = 1.0_double/deli
            if (all(ffto > 0) .and. all(ffto <= lay%o%kgroup%locdims)) deli(ffto(1),ffto(2),ffto(3)) = 0.0_double
          end select
        end select
        call glean(thy(lay))
        if (error("Exit layout_mod::fdelinv_lay")) continue
      end subroutine

      subroutine filter_r_lay(data,lay,sc)
!doc$ subroutine filter(data,lay,sc)
        real(double), dimension(:,:,:), intent(inout) :: data
        type(layout_obj) :: lay
        integer, intent(in) :: sc
!       requires: data be in reciprocal space.
!                 data distribution be D_TYPE in scope sc.
!       modifies: data
!       effects: Sets data = 0 where lay%o%filter = .true.
!       notes: This is a special routine for real-valued data in a reciprocal-space representation.

!cod$
        call my(lay)
        select case (sc)
        case (CONFIG)
          where (lay%o%config%filter) data = 0.0_double
        case (SGROUP)
          where (lay%o%sgroup%filter) data = 0.0_double
        case (KGROUP)
          where (lay%o%kgroup%filter) data = 0.0_double
        end select
        call glean(thy(lay))
      end subroutine

      subroutine filter_c_lay(data,lay,sc)
!doc$ subroutine filter(data,lay,sc)
        complex(double), dimension(:,:,:), intent(inout) :: data
        type(layout_obj) :: lay
        integer, intent(in) :: sc
!       requires: data be in reciprocal space.
!                 data distribution be D_TYPE in scope sc.
!       modifies: data
!       effects: Sets data = 0 where lay%o%filter = .true.

!cod$
        call my(lay)
        select case (sc)
        case (CONFIG)
          where (lay%o%config%filter) data = (0.0_double,0.0_double)
        case (SGROUP)
          where (lay%o%sgroup%filter) data = (0.0_double,0.0_double)
        case (KGROUP)
          where (lay%o%kgroup%filter) data = (0.0_double,0.0_double)
        end select
        call glean(thy(lay))
      end subroutine

      subroutine diary_lay(lay)
!doc$ subroutine diary(lay)
        type(layout_obj) :: lay
!       modifies: Output stream.
!       effects: Writes a summary of layout information to the diary output.

!cod$
        integer :: i, nsize, bigslice, slice
        integer :: i1, i2, i3, npw
        integer, dimension(3) :: dims, bdims, n2, n3, n5, n7
        real(double), dimension(:,:,:), pointer :: g2

        call my(lay)

        nullify( g2 )

        if (i_access(diaryfile())) then

          write(x_unit(diaryfile()),'(/,t4,"Mesh:")')

          dims = lay%o%dims
          do i = 1,3
            call get_nf(dims(i),n2(i),n3(i),n5(i),n7(i))
          end do
          write(x_unit(diaryfile()),'(/,t6,"dimension 1 = ",i3,"  =  2^",i1," * 3^",i1," * 5^",i1," * 7^",i1)') &
              dims(1), n2(1), n3(1), n5(1), n7(1)
          write(x_unit(diaryfile()),'(  t6,"dimension 2 = ",i3,"  =  2^",i1," * 3^",i1," * 5^",i1," * 7^",i1)') &
              dims(2), n2(2), n3(2), n5(2), n7(2)
          write(x_unit(diaryfile()),'(  t6,"dimension 3 = ",i3,"  =  2^",i1," * 3^",i1," * 5^",i1," * 7^",i1)') &
              dims(3), n2(3), n3(3), n5(3), n7(3)
          write(x_unit(diaryfile()),'(/,t6,"Number of mesh points = ",i0)') product(dims)

          if ( (mpi_nkgroups() == 1) .and. (mpi_nsgroups() == 1) ) then

            write(x_unit(diaryfile()),'(/,t6,"Distribution:",/)')
            bdims = lay%o%config%bdims
            do i = 1,3
              if (bdims(i) == 1) then
                write(x_unit(diaryfile()),'(t8,"dimension ",i1,": on processor")') i
              else
                call subdivide(bdims(i),dims(i),nsize,slice,bigslice)
                if (slice == 1) then
                  write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," group with ",i0," planes")') i, slice, nsize
                else
                  if (nsize == 1) then
                    write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," plane per group")') i, slice,nsize
                  else
                    write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," planes per group")') i,slice,nsize
                  end if
                end if
                if (bigslice /= 0) then
                  if (bigslice == 1) then
                    write(x_unit(diaryfile()),'(t23,i0," group with ",i0," planes")') bigslice, (nsize + 1)
                  else
                    write(x_unit(diaryfile()),'(t21,i0," groups with ",i0," planes per group")') bigslice, (nsize + 1)
                  end if
                end if
              end if
            end do
            if (product(bdims) < mpi_nprocs(CONFIG)) then
              write(x_unit(diaryfile()),'(/,t6,i0," processes idle")') (mpi_nprocs(CONFIG) - product(bdims))
            end if

          else

            write(x_unit(diaryfile()),'(/,t6,"config distribution:",/)')
            bdims = lay%o%config%bdims
            do i = 1,3
              if (bdims(i) == 1) then
                write(x_unit(diaryfile()),'(t8,"dimension ",i1,": on processor")') i
              else
                call subdivide(bdims(i),dims(i),nsize,slice,bigslice)
                if (slice == 1) then
                  write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," group with ",i0," planes")') i, slice, nsize
                else
                  if (nsize == 1) then
                    write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," plane per group")') i, slice,nsize
                  else
                    write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," planes per group")') i,slice,nsize
                  end if
                end if
                if (bigslice /= 0) then
                  if (bigslice == 1) then
                    write(x_unit(diaryfile()),'(t21,i0," group with ",i0," planes")') bigslice, (nsize + 1)
                  else
                    write(x_unit(diaryfile()),'(t21,i0," groups with ",i0," planes per group")') bigslice, (nsize + 1)
                  end if
                end if
              end if
            end do
            if (product(bdims) < mpi_nprocs(CONFIG)) then
              write(x_unit(diaryfile()),'(/,t6,i0," processes idle")') (mpi_nprocs(CONFIG) - product(bdims))
            end if

            if (mpi_nsgroups() /= 1) then

              write(x_unit(diaryfile()),'(/,t6,"sgroup distribution:",/)')
              bdims = lay%o%sgroup%bdims
              dims = lay%o%dims
              do i = 1,3
                if (bdims(i) == 1) then
                  write(x_unit(diaryfile()),'(t8,"dimension ",i1,": on processor")') i
                else
                  call subdivide(bdims(i),dims(i),nsize,slice,bigslice)
                  if (slice == 1) then
                    write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," group with ",i0," planes")') i, slice, nsize
                  else
                    if (nsize == 1) then
                      write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," plane per group")') i, slice,nsize
                    else
                      write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," planes per group")') i,slice,nsize
                    end if
                  end if
                  if (bigslice /= 0) then
                    if (bigslice == 1) then
                      write(x_unit(diaryfile()),'(t21,i0," group with ",i0," planes")') bigslice, (nsize + 1)
                    else
                      write(x_unit(diaryfile()),'(t21,i0," groups with ",i0," planes per group")') bigslice, (nsize + 1)
                    end if
                  end if
                end if
              end do
              if (product(bdims) < mpi_nprocs(SGROUP)) then
                write(x_unit(diaryfile()),'(/,t6,i0," processes idle")') (mpi_nprocs(SGROUP) - product(bdims))
              end if

            end if

            if (mpi_nkgroups() /= 1) then

              write(x_unit(diaryfile()),'(/,t6,"kgroup distribution:",/)')
              bdims = lay%o%kgroup%bdims
              dims = lay%o%dims
              do i = 1,3
                if (bdims(i) == 1) then
                  write(x_unit(diaryfile()),'(t8,"dimension ",i1,": on processor")') i
                else
                  call subdivide(bdims(i),dims(i),nsize,slice,bigslice)
                  if (slice == 1) then
                    write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," group with ",i0," planes")') i, slice, nsize
                  else
                    if (nsize == 1) then
                      write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," plane per group")') i, slice,nsize
                    else
                      write(x_unit(diaryfile()),'(t8,"dimension ",i1,": ",i0," groups with ",i0," planes per group")') i,slice,nsize
                    end if
                  end if
                  if (bigslice /= 0) then
                    if (bigslice == 1) then
                      write(x_unit(diaryfile()),'(t21,i0," group with ",i0," planes")') bigslice, (nsize + 1)
                    else
                      write(x_unit(diaryfile()),'(t21,i0," groups with ",i0," planes per group")') bigslice, (nsize + 1)
                    end if
                  end if
                end if
              end do
              if (product(bdims) < mpi_nprocs(KGROUP)) then
                write(x_unit(diaryfile()),'(/,t6,i0," processes idle")') (mpi_nprocs(KGROUP) - product(bdims))
              end if

            end if

          end if

          call fdel(g2,lay,S_TYPE)
          npw = 0
          do i3 = 1,size(g2,3)
          do i2 = 1,size(g2,2)
          do i1 = 1,size(g2,1)
            if (g2(i1,i2,i3) <= lay%o%cutoff) npw = npw + 1
          end do
          end do
          end do
          deallocate( g2 )
          write(x_unit(diaryfile()),'(/,t6,"Plane wave cutoff energy = ",f0.2," Ryd")') lay%o%cutoff
          write(x_unit(diaryfile()),'(/,t8,"Number of plane waves = ",i0)') npw

        end if

100     if (associated( g2 )) deallocate( g2 )

        call glean(thy(lay))

        if (error("Exit layout_mod::diary_lay")) continue

      end subroutine

      subroutine write_restart_lay(lay,nrestf)
!doc$ subroutine write_restart(lay,nrestf)
        type(layout_obj) :: lay
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes lay information to nrestf.

!cod$      
        integer(long) :: dsize, iosl, ndata

        call my(lay)
        call my(nrestf)

        ! start the LAYOUT block
        if (i_access(nrestf)) call startblock(nrestf,"LAYOUT")

        ! write the CUTOFF
        if (i_access(nrestf)) then
          call writetag(nrestf,"CUTOFF")
          dsize = sizeof_double ; ndata = 1
          call writef(lay%o%cutoff,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! end the LAYOUT block
        if (i_access(nrestf)) call endblock(nrestf)

        call glean(thy(lay))
        call glean(thy(nrestf))

        if (error("Exit layout_mod::write_restart_lay")) continue

      end subroutine

! private routines

      subroutine set_layout_dims_i(layr,dims)
        type(layout_rep) :: layr
        integer, dimension(3), intent(in) :: dims
        integer :: i, j, myp, np
        integer, dimension(3) :: bdims, locdims, base, block

        layr%dims = dims
        layr%origin = (/1,1,1/)
        layr%fft_origin = (/1,1,1/)

        np = mpi_nprocs(CONFIG)
        myp = mpi_myproc(CONFIG)
        do
          if (np <= dims(3)) then
            bdims(1) = 1
            bdims(2) = 1
            bdims(3) = np
            exit
          elseif (np <= dims(2)) then
            bdims(1) = 1
            bdims(2) = np
            bdims(3) = 1
            exit
          elseif (np <= dims(1)) then
            bdims(1) = np
            bdims(2) = 1
            bdims(3) = 1
            exit
          elseif (np <= dims(2)*dims(3)) then
            call factor(np,i,j)
            bdims(1) = 1
            bdims(2) = i
            bdims(3) = j
            exit
          elseif (np <= dims(3)*dims(1)) then
            call factor(np,i,j)
            bdims(1) = i
            bdims(2) = 1
            bdims(3) = j
            exit
          elseif (np <= dims(1)*dims(2)) then
            call factor(np,i,j)
            bdims(1) = i
            bdims(2) = j
            bdims(3) = 1
            exit
          end if
          np = np - 1
        end do
        if (myp < product(bdims)) then
          block = map3d(bdims,myp+1) - (/1,1,1/)
          do i = 1,3
            call subdivide(block(i),bdims(i),1,dims(i),base(i),j,locdims(i))
          end do
        else
          call warn("mpi_myproc(CONFIG) >= product(bdims)")
          locdims = (/0,0,0/)
          base = (/0,0,0/)
          block = (/0,0,0/)
        end if
        layr%config%locdims = locdims
        layr%config%bdims = bdims
        layr%config%base = base
        layr%config%block = block

        np = mpi_nprocs(SGROUP)
        myp = mpi_myproc(SGROUP)
        do
          if (np <= dims(3)) then
            bdims(1) = 1
            bdims(2) = 1
            bdims(3) = np
            exit
          elseif (np <= dims(2)) then
            bdims(1) = 1
            bdims(2) = np
            bdims(3) = 1
            exit
          elseif (np <= dims(1)) then
            bdims(1) = np
            bdims(2) = 1
            bdims(3) = 1
            exit
          elseif (np <= dims(2)*dims(3)) then
            call factor(np,i,j)
            bdims(1) = 1
            bdims(2) = i
            bdims(3) = j
            exit
          elseif (np <= dims(3)*dims(1)) then
            call factor(np,i,j)
            bdims(1) = i
            bdims(2) = 1
            bdims(3) = j
            exit
          elseif (np <= dims(1)*dims(2)) then
            call factor(np,i,j)
            bdims(1) = i
            bdims(2) = j
            bdims(3) = 1
            exit
          end if
          np = np - 1
        end do
        if (myp < product(bdims)) then
          block = map3d(bdims,myp+1) - (/1,1,1/)
          do i = 1,3
            call subdivide(block(i),bdims(i),1,dims(i),base(i),j,locdims(i))
          end do
        else
          call warn("mpi_myproc(SGROUP) >= product(bdims)")
          locdims = (/0,0,0/)
          base = (/0,0,0/)
          block = (/0,0,0/)
        end if
        layr%sgroup%locdims = locdims
        layr%sgroup%bdims = bdims
        layr%sgroup%base = base
        layr%sgroup%block = block

        np = mpi_nprocs(KGROUP)
        myp = mpi_myproc(KGROUP)
        do
          if (np <= dims(3)) then
            bdims(1) = 1
            bdims(2) = 1
            bdims(3) = np
            exit
          elseif (np <= dims(2)) then
            bdims(1) = 1
            bdims(2) = np
            bdims(3) = 1
            exit
          elseif (np <= dims(1)) then
            bdims(1) = np
            bdims(2) = 1
            bdims(3) = 1
            exit
          elseif (np <= dims(2)*dims(3)) then
            call factor(np,i,j)
            bdims(1) = 1
            bdims(2) = i
            bdims(3) = j
            exit
          elseif (np <= dims(3)*dims(1)) then
            call factor(np,i,j)
            bdims(1) = i
            bdims(2) = 1
            bdims(3) = j
            exit
          elseif (np <= dims(1)*dims(2)) then
            call factor(np,i,j)
            bdims(1) = i
            bdims(2) = j
            bdims(3) = 1
            exit
          end if
          np = np - 1
        end do
        if (myp < product(bdims)) then
          block = map3d(bdims,myp+1) - (/1,1,1/)
          do i = 1,3
            call subdivide(block(i),bdims(i),1,dims(i),base(i),j,locdims(i))
          end do
        else
          call warn("mpi_myproc(KGROUP) >= product(bdims)")
          locdims = (/0,0,0/)
          base = (/0,0,0/)
          block = (/0,0,0/)
        end if
        layr%kgroup%locdims = locdims
        layr%kgroup%bdims = bdims
        layr%kgroup%base = base
        layr%kgroup%block = block

      end subroutine

      subroutine form_filters_i(layr)
        type(layout_rep) :: layr

        integer :: i1, i2, i3
        integer, dimension(3) :: d
        real(double), dimension(3) :: o, p, q, rd
        real(double), dimension(3) :: half, one

        rd = real(layr%dims,double)
        one = real((/1,1,1/),double)
        half = 0.5_double*one

        d = layr%config%locdims
        o = real(layr%fft_origin - (layr%config%base - (/1,1,1/)),double)
        allocate( layr%config%filter(layr%config%locdims(1),layr%config%locdims(2),layr%config%locdims(3)) )
        do i3 = 1,d(3)
        do i2 = 1,d(2)
        do i1 = 1,d(1)
          p = (real((/i1,i2,i3/),double) - o)/rd
          p = (mod(p+half,one) - half)*rd
          q = lat2f(layr%lattice,p)
          if ((q(1)**2 + q(2)**2 + q(3)**2) <= layr%cutoff) then
            layr%config%filter(i1,i2,i3) = .false.
          else
            layr%config%filter(i1,i2,i3) = .true.
          end if
        end do
        end do
        end do

        d = layr%sgroup%locdims
        o = real(layr%fft_origin - (layr%sgroup%base - (/1,1,1/)),double)
        allocate( layr%sgroup%filter(layr%sgroup%locdims(1),layr%sgroup%locdims(2),layr%sgroup%locdims(3)) )
        do i3 = 1,d(3)
        do i2 = 1,d(2)
        do i1 = 1,d(1)
          p = (real((/i1,i2,i3/),double) - o)/rd
          p = (mod(p+half,one) - half)*rd
          q = lat2f(layr%lattice,p)
          if ((q(1)**2 + q(2)**2 + q(3)**2) <= layr%cutoff) then
            layr%sgroup%filter(i1,i2,i3) = .false.
          else
            layr%sgroup%filter(i1,i2,i3) = .true.
          end if
        end do
        end do
        end do

        d = layr%kgroup%locdims
        o = real(layr%fft_origin - (layr%kgroup%base - (/1,1,1/)),double)
        allocate( layr%kgroup%filter(layr%kgroup%locdims(1),layr%kgroup%locdims(2),layr%kgroup%locdims(3)) )
        do i3 = 1,d(3)
        do i2 = 1,d(2)
        do i1 = 1,d(1)
          p = (real((/i1,i2,i3/),double) - o)/rd
          p = (mod(p+half,one) - half)*rd
          q = lat2f(layr%lattice,p)
          if ((q(1)**2 + q(2)**2 + q(3)**2) <= layr%cutoff) then
            layr%kgroup%filter(i1,i2,i3) = .false.
          else
            layr%kgroup%filter(i1,i2,i3) = .true.
          end if
        end do
        end do
        end do

      end subroutine

      end module

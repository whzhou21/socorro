!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module grid_mod
!doc$ module grid_mod

!     One datatype is available here: type(grid_obj).

!     grid_mod provides a common medium of exchange for mesh data and basic type conversions.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use diary_mod
      use math_mod
      use layout_mod
      use lattice_mod
      use ghost_mod

!cod$
      implicit none
      private

      ! Grid kinds
      integer, parameter :: EMPTY_KIND = 0  ! empty
      integer, parameter :: RS_KIND    = 1  ! real serial
      integer, parameter :: RD_KIND    = 2  ! real distributed
      integer, parameter :: CSP_KIND   = 3  ! complex serial position
      integer, parameter :: CDP_KIND   = 4  ! complex distributed position
      integer, parameter :: CSF_KIND   = 5  ! complex serial fourier
      integer, parameter :: CDF_KIND   = 6  ! complex distributed fourier

      ! File formats for output
      integer, parameter :: MATLAB     = 0  ! Matlab format
      integer, parameter :: AMIRA      = 1  ! Amira/Avizo structured grid format
      integer, parameter :: VTK        = 2  ! Visualization Toolkit format
      

      type, public :: grid_rep
        integer :: ref
        type(ghost) :: g
        integer :: scope
        integer :: type
        type(layout_obj) :: layout
        real(double), dimension(:,:,:), pointer :: rdata
        complex(double), dimension(:,:,:), pointer :: cdata
      end type 

      type, public :: grid_obj
        private
        integer :: ref
        type(grid_rep), pointer :: o
      end type

!doc$
      public :: grid
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_scope
      public :: x_type
      public :: x_layout
      public :: put
      public :: take
      public :: transform
      public :: norm
      public :: empty
      public :: filter
      public :: get_normalization
      public :: set_normalization
      public :: distance
      public :: saxpby
      public :: merge_grid_density
      public :: xsum
      public :: xdif
      public :: sgroup_to_kgroup
      public :: config_to_kgroup
      public :: read_restart
      public :: write_restart
      public :: write_to_file
      public :: EMPTY_KIND
      public :: RS_KIND
      public :: RD_KIND
      public :: CSP_KIND
      public :: CDP_KIND
      public :: CSF_KIND
      public :: CDF_KIND
      public :: MATLAB
      public :: AMIRA
      public :: VTK

!cod$
      interface grid
        module procedure constructor_grid
      end interface
      interface my
        module procedure my_grid, my_new_grid
      end interface
      interface thy
        module procedure thy_grid
      end interface
      interface glean
        module procedure glean_grid
      end interface
      interface bequeath
        module procedure bequeath_grid
      end interface
      interface assignment(=)
        module procedure assign_grid
      end interface
      interface x_ref
        module procedure grid_ref
      end interface
      interface x_ghost
        module procedure grid_ghost
      end interface
      interface x_scope
        module procedure grid_scope
      end interface
      interface x_type
        module procedure grid_type
      end interface
      interface x_layout
        module procedure grid_layout
      end interface
      interface put
        module procedure put_ptr_grid_r, put_ptr_grid_c
      end interface
      interface take
        module procedure take_ptr_grid_r, take_ptr_grid_c
      end interface
      interface transform
        module procedure transform_grid
      end interface
      interface norm
        module procedure norm_grid
      end interface
      interface empty
        module procedure empty_grid
      end interface
      interface filter
        module procedure filter_grid
      end interface
      interface get_normalization
        module procedure get_normalization_grid_c
      end interface
      interface set_normalization
        module procedure set_normalization_grid_r, set_normalization_grid_c
      end interface
      interface distance
        module procedure distance_grid
      end interface
      interface saxpby
        module procedure saxpby_grid
      end interface
      interface merge_grid_density
        module procedure merge_grid_density_grid
      end interface
      interface xsum
        module procedure xsum_grid
      end interface
      interface xdif
        module procedure xdif_grid
      end interface
      interface sgroup_to_kgroup
        module procedure sgroup_to_kgroup_grid
      end interface
      interface config_to_kgroup
        module procedure config_to_kgroup_grid
      end interface
      interface read_restart
        module procedure read_restart_grid
      end interface
      interface write_restart
        module procedure write_restart_grid
      end interface
      interface write_to_file
        module procedure write_to_file_grid
      end interface

      contains

! public routines

      function constructor_grid(lay,sc) result(g)
!doc$ function grid(lay,sc) result(g)
        type(grid_obj) :: g
        type(layout_obj) :: lay
        integer, intent(in) :: sc
!       effects : Creates a new grid_obj with type = EMPTY_KIND, and scope = sc.

!cod$
        g%ref = 0
        allocate( g%o )
        g%o%ref = 0
        g%o%g = x_ghost()
        g%o%scope = sc
        g%o%type = EMPTY_KIND
        call my(lay,g%o%layout)
        nullify( g%o%rdata )
        nullify( g%o%cdata )
      end function

      subroutine my_grid(g)
!doc$ subroutine my(g)
        type(grid_obj) :: g

!cod$
        g%ref = g%ref + 1
        g%o%ref = g%o%ref + 1
      end subroutine

      subroutine my_new_grid(gi,g)
!doc$ subroutine my(gi,g)
        type(grid_obj) :: gi, g

!cod$
        g%ref = 1
        g%o => gi%o
        g%o%ref = g%o%ref + 1
      end subroutine

      function thy_grid(g) result(go)
!doc$ function thy(g) result(go)
        type(grid_obj) :: g, go

!cod$
        g%ref = g%ref - 1
        g%o%ref = g%o%ref - 1
        go%ref = g%ref
        go%o => g%o
      end function

      subroutine glean_grid(g)
!doc$ subroutine glean(g)
        type(grid_obj) :: g

!cod$
        if (g%o%ref < 1) then
          select case (g%o%type)
          case (EMPTY_KIND)
            continue
          case (RD_KIND,RS_KIND)
            deallocate( g%o%rdata )
          case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
            deallocate( g%o%cdata )
          end select
          call glean(thy(g%o%layout))
          deallocate( g%o )
        end if
      end subroutine

      subroutine bequeath_grid(g)
!doc$ subroutine bequeath(g)
        type(grid_obj) :: g

!cod$
        continue
      end subroutine

      subroutine assign_grid(g,g2)
!doc$ subroutine assign(g,g2)
        type(grid_obj), intent(inout) :: g
        type(grid_obj), intent(in) :: g2
!       errors: g and g2 have different layouts. g and g2 have different scopes.

!cod$
        type(grid_obj) :: gt
        call my(g2)
        if (error(x_ghost(g%o%layout) /= x_ghost(g2%o%layout),"ERROR: incompatible layouts")) goto 100
        if (error(g%o%scope /= g2%o%scope,"ERROR: incompatible scopes")) goto 100
        gt%o => g%o
        g%o%ref = g%o%ref - g%ref
        g%o => g2%o
        g%o%ref = g%o%ref + g%ref
        call glean(gt)
100     call glean(thy(g2))
        if (error("Exit grid_mod::assign_grid")) continue
      end subroutine

      function grid_ref(g) result(r)
!doc$ function x_ref(g) result(r)
        type(grid_obj) :: g
        integer, dimension(2) :: r
!       effects: Returns g%ref and g%o%ref.

!cod$
        r(1) = g%ref
        r(2) = g%o%ref
        call glean(g)
      end function

      function grid_ghost(g) result(gg)
!doc$ function x_ghost(g) result(gg)
        type(grid_obj) :: g
        type(ghost) :: gg
!       effects: Returns g%o%g.

!cod$ 
        call my(g)
        gg = g%o%g
        call glean(thy(g))
      end function

      function grid_scope(g) result(s)
!doc$ function x_scope(g) result(s)
        type(grid_obj) :: g
        integer :: s
!       effects: Returns g%o%scope.

!cod$
        call my(g)
        s = g%o%scope
        call glean(thy(g))
      end function

      function grid_type(g) result(t)
!doc$ function x_type(g) result(t)
        type(grid_obj) :: g
        integer :: t
!       effects: Returns g%o%type.

!cod$
        call my(g)
        t = g%o%type
        call glean(thy(g))
      end function

      function grid_layout(g) result(lay)
!doc$ function x_layout(g) result(lay)
        type(grid_obj) :: g
        type(layout_obj) :: lay
!       effects: Returns g%o%layout.

!cod$
        call my(g)
        call my(g%o%layout,lay)
        call bequeath(thy(lay))
        call glean(thy(g))
      end function

      subroutine put_ptr_grid_r(data,g,kind)
!doc$ subroutine put(data,g,kind)
        real(double), dimension(:,:,:), pointer :: data
        type(grid_obj) :: g
        integer, intent(in) :: kind
!       modifies: data, g
!       effects: If (empty?(g)) x_type(g) = kind, put from data to grid and nullify(data).
!       errors: Incorrect data distribution. Passes errors.

!cod$
        integer :: tpo
        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp
        call my(g)
        select case (kind)
        case (RS_KIND)
          if (error(.not.consistent(data,g%o%layout,S_TYPE,g%o%scope),"ERROR: incorrect data distribution")) goto 100
        case (RD_KIND)
          if (error(.not.consistent(data,g%o%layout,D_TYPE,g%o%scope),"ERROR: incorrect data distribution")) goto 100
        end select
        call own_i(g)
        call empty_i(g)
        rtmp => data
        nullify( ctmp )
        nullify( data )
        tpo = kind
        call transfer_grid_i(g%o%layout,g%o%scope,g%o%type,g%o%rdata,g%o%cdata,tpo,rtmp,ctmp) ; if (error()) goto 100
        g%o%g = x_ghost()
100     call glean(thy(g))
        if (error("Exit grid_mod::put_ptr_grid_r")) continue
      end subroutine

      subroutine put_ptr_grid_c(data,g,kind)
!doc$ subroutine put(data,g,kind)
        complex(double), dimension(:,:,:), pointer :: data
        type(grid_obj) :: g
        integer, intent(in) :: kind
!       modifies: data, g
!       effects: If (empty?(g)) x_type(g) = kind, put from data to grid and nullify(data).
!       errors: Incorrect data distribution. Passes errors.

!cod$
        integer :: tpo
        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp
        call my(g)
        select case (kind)
        case (CSP_KIND,CSF_KIND)
          if (error(.not.consistent(data,g%o%layout,S_TYPE,g%o%scope),"ERROR: incorrect data distribution")) goto 100
        case (CDP_KIND,CDF_KIND)
          if (error(.not.consistent(data,g%o%layout,D_TYPE,g%o%scope),"ERROR: incorrect data distribution")) goto 100
        end select
        call own_i(g)
        call empty_i(g)
        ctmp => data
        nullify( rtmp )
        nullify( data )
        tpo = kind
        call transfer_grid_i(g%o%layout,g%o%scope,g%o%type,g%o%rdata,g%o%cdata,tpo,rtmp,ctmp) ; if (error()) goto 100
        g%o%g = x_ghost()
100     call glean(thy(g))
        if (error("Exit grid_mod::put_ptr_grid_c")) continue
      end subroutine

      subroutine take_ptr_grid_r(data,g,kind)
!doc$ subroutine take(data,g,kind)
        real(double), dimension(:,:,:), pointer :: data
        type(grid_obj) :: g
        integer, intent(in) :: kind
!       modifies: g, data
!       effects: Gets pointer to data from g (transforming if needed) and sets x_type(g) = EMPTY_KIND.
!       errors: Passes errors.

!cod$
        integer :: tpo 
        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp
        call my(g)
        call own_i(g)
        nullify( rtmp )
        nullify( ctmp )
        tpo = kind
        call transfer_grid_i(g%o%layout,g%o%scope,tpo,rtmp,ctmp,g%o%type,g%o%rdata,g%o%cdata) ; if (error()) goto 100
        g%o%g = x_ghost()
        data => rtmp
100     call glean(thy(g))
        if (error("Exit grid_mod::take_ptr_grid_r")) continue
      end subroutine

      subroutine take_ptr_grid_c(data,g,kind)
!doc$ subroutine take(data,g,kind)
        complex(double), dimension(:,:,:), pointer :: data
        type(grid_obj) :: g
        integer, intent(in) :: kind
!       modifies: g and data
!       effects: Gets pointer to data from g (transforming if needed) and sets x_type(g) = EMPTY_KIND.
!       errors: Passes errors.

!cod$
        integer :: tpo
        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp
        call my(g)
        call own_i(g)
        nullify( rtmp )
        nullify( ctmp )
        tpo = kind
        call transfer_grid_i(g%o%layout,g%o%scope,tpo,rtmp,ctmp,g%o%type,g%o%rdata,g%o%cdata) ; if (error()) goto 100
        g%o%g = x_ghost()
        data => ctmp
100     call glean(thy(g))
        if (error("Exit grid_mod::take_ptr_grid_c")) continue
      end subroutine

      subroutine transform_grid(g,kind)
!doc$ subroutine transform(g,kind)
        type(grid_obj) :: g
        integer, intent(in) :: kind
!       modifies: g
!       effects: Transforms g data such that x_type(g) = kind.
!       errors: Passes errors.

!cod$
        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp
        call my(g)
        select case (kind)
        case (EMPTY_KIND)
          call own_i(g)
          call empty_i(g)
          g%o%g = x_ghost()
        case (RS_KIND,RD_KIND)
          call take(rtmp,g,kind) ; if (error()) goto 100
          call put(rtmp,g,kind) ; if (error()) goto 100
        case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
          call take(ctmp,g,kind) ; if (error()) goto 100
          call put(ctmp,g,kind) ; if (error()) goto 100
        end select
100     call glean(thy(g))
        if (error("Exit grid_mod::transform_grid: ")) continue
      end subroutine

      function norm_grid(g) result(n)
!doc$ subroutine norm(g) result(n)
        type(grid_obj) :: g
        real(double) :: n
!       effects: Returns sqrt( sum( abs(grid data) ) ).
!       errors: Passes errors.

!cod$
        integer :: kind
        real(double) :: x_local, x_global
        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp
        call my(g)
        kind = x_type(g)
        select case (kind)
        case (EMPTY_KIND)
          n = 0.0_double
        case (RS_KIND)
          call take(rtmp,g,kind) ; if (error()) goto 100
          x_local = sum(rtmp**2)
          call put(rtmp,g,kind) ; if (error()) goto 100
          n = sqrt(x_local)
        case (RD_KIND)
          call take(rtmp,g,kind) ; if (error()) goto 100
          x_local = sum(rtmp**2)
          call put(rtmp,g,kind) ; if (error()) goto 100
          call allreduce(g%o%scope,MPI_SUM,x_local,x_global)
          n = sqrt(x_global)
        case (CSP_KIND,CSF_KIND)
          call take(ctmp,g,kind) ; if (error()) goto 100
          x_local = sum(abs(ctmp)**2)
          call put(ctmp,g,kind) ; if (error()) goto 100
          n = sqrt(x_local)
        case (CDP_KIND,CDF_KIND)
          call take(ctmp,g,kind) ; if (error()) goto 100
          x_local = sum(abs(ctmp)**2)
          call put(ctmp,g,kind) ; if (error()) goto 100
          call allreduce(g%o%scope,MPI_SUM,x_local,x_global)
          n = sqrt(x_global)
        end select
100     call glean(thy(g))
        if (error("Exit grid_mod::norm_grid: ")) continue
      end function

      subroutine empty_grid(g)
!doc$ subroutine empty(g)
        type(grid_obj) :: g
!       modifies: g
!       effects: makes g empty, whether or not g was empty before

!cod$
        call my(g)
        call own_i(g)
        call empty_i(g)
        if (error()) goto 100
        g%o%g = x_ghost()
100     call glean(thy(g))
        if (error("Exit grid_mod::empty_grid: ")) continue
      end subroutine

      subroutine filter_grid(g)
!doc$ subroutine filter(g)
        type(grid_obj) :: g
!       modifies: g
!       effects: Filters g data. If g%o%type = EMPTY_KIND, returns g with no modifications.

!cod$
        complex(double), dimension(:,:,:), pointer :: c1
        call my(g)
        if (g%o%type == EMPTY_KIND) goto 100
        call own_i(g)
        call take(c1,g,CDF_KIND)
        call filter(c1,g%o%layout,g%o%scope)
        call put(c1,g,CDF_KIND)
        g%o%g = x_ghost()
100     call glean(thy(g))
        if (error("Exit grid_mod::filter_grid: ")) continue
      end subroutine

      subroutine get_normalization_grid_c(g,n)
!doc$ subroutine get_normalization(g,n)
        type(grid_obj) :: g
        complex(double) :: n
!       effects: Returns the normalization of g. if g%o%type = EMPTY_KIND, returns 0.
!       notes: g may be changed if the mesh sizes permit aliasing (hence the call to own_i and the new ghost).

!cod$
        integer, dimension(3) :: ffto
        complex(double) :: ns
        complex(double), dimension(:,:,:), pointer :: c1
        call my(g)
        if (g%o%type == EMPTY_KIND) then
          n = (0.0_double,0.0_double)
          goto 100
        end if
        call own_i(g)
        call take(c1,g,CDF_KIND)
        ffto = x_fft_origin(g%o%layout,D_TYPE,g%o%scope)
        ns = (0.0_double,0.0_double)
        if (all(ffto /= 0)) ns = c1(ffto(1),ffto(2),ffto(3))
        call allreduce(g%o%scope,MPI_SUM,ns,n) ; if (error()) goto 100
        call put(c1,g,CDF_KIND)
        g%o%g = x_ghost()
100     call glean(thy(g))
      end subroutine

      subroutine set_normalization_grid_r(g,n)
!doc$ subroutine set_normalization(g,n)
        type(grid_obj) :: g
        real(double) :: n
!       effects: Sets the normalization of g to n. if g%o%type = EMPTY_KIND, returns with no changes.

!cod$
        integer, dimension(3) :: ffto
        real(double), dimension(:,:,:), pointer :: r1
        call my(g)
        if (g%o%type == EMPTY_KIND) goto 100
        call own_i(g)
        call take(r1,g,RD_KIND)
        ffto = x_fft_origin(g%o%layout,D_TYPE,g%o%scope)
        if (all(ffto /= 0)) r1(ffto(1),ffto(2),ffto(3)) = n
        call put(r1,g,RD_KIND)
        g%o%g = x_ghost()
100     call glean(thy(g))
      end subroutine

      subroutine set_normalization_grid_c(g,n)
!doc$ subroutine set_normalization(g,n)
        type(grid_obj) :: g
        complex(double) :: n
!       effects: Sets the normalization of g to n. if g%o%type = EMPTY_KIND, returns with no changes.

!cod$
        integer, dimension(3) :: ffto
        complex(double), dimension(:,:,:), pointer :: c1
        call my(g)
        if (g%o%type == EMPTY_KIND) goto 100
        call own_i(g)
        call take(c1,g,CDF_KIND)
        ffto = x_fft_origin(g%o%layout,D_TYPE,g%o%scope)
        if (all(ffto /= 0)) c1(ffto(1),ffto(2),ffto(3)) = n
        call put(c1,g,CDF_KIND)
        g%o%g = x_ghost()
100     call glean(thy(g))
      end subroutine

      function distance_grid(g1,g2) result(dist)
!doc$ function distance(g1,g2) result(dist)
        type(grid_obj) :: g1, g2
        real(double) :: dist
!       requires: g1 and g2 data are complex valued.
!       effects: Returns the rms difference between g1 and g2.
!       errors: g1 or g2 empty.

!cod$
        real(double) :: s_local, s_global
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: g1_tmp, g2_tmp

        call my(g1)
        call my(g2)

        nullify( c1, c2 )

        if (error(g1%o%type == EMPTY_KIND,"ERROR: g1 is empty")) goto 100
        if (error(g2%o%type == EMPTY_KIND,"ERROR: g2 is empty")) goto 100

        if ((g1%o%type == CDF_KIND) .and. (g2%o%type == CDF_KIND)) then
          s_local = sum(abs(g2%o%cdata - g1%o%cdata)**2)
          call allreduce(g1%o%scope,MPI_SUM,s_local,s_global)
          dist = sqrt(s_global)
        elseif ((g1%o%type == CSF_KIND) .and. (g2%o%type == CSF_KIND)) then
          dist = sqrt(sum(abs(g2%o%cdata - g1%o%cdata)**2))
        else
          call my(g1,g1_tmp)
          call my(g2,g2_tmp)
          call take(c1,g1_tmp,CDF_KIND)
          call take(c2,g2_tmp,CDF_KIND)
          s_local = sum(abs(c2 - c1)**2)
          call allreduce(g1%o%scope,MPI_SUM,s_local,s_global)
          dist = sqrt(s_global)
          call put(c1,g1_tmp,CDF_KIND)
          call put(c2,g2_tmp,CDF_KIND)
          call glean(thy(g1_tmp))
          call glean(thy(g2_tmp))
        end if

100     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(g1))
        call glean(thy(g2))

        if (error("Exit grid_mod::distance_grid")) continue

      end function

      subroutine saxpby_grid(a,g1,b,g2)
!doc$ subroutine saxpby(a,g1,b,g2)
        real(double), intent(in) :: a, b
        type(grid_obj) :: g1, g2
!       requires: g1 and g2 are compatible.
!       modifies: g1
!       effects: Returns g1 = a*g1 + b*g2.
!       notes: If g1 or g2 are empty, they are treated as if they are filled with zeros.

!cod$
        integer :: g1_type
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: g2_tmp

        call my(g1)
        call my(g2)

        nullify( r1, r2, c1, c2 )

        g1_type = g1%o%type
        select case (g1_type)
        case (EMPTY_KIND)
          select case (g2%o%type)
          case (RS_KIND)
            call alloc(r1,g2%o%layout,S_TYPE)
            r1 = b*g2%o%rdata
            call put(r1,g1,g2%o%type)
          case (RD_KIND)
            call alloc(r1,g2%o%layout,D_TYPE,g2%o%scope)
            r1 = b*g2%o%rdata
            call put(r1,g1,g2%o%type)
          case (CSP_KIND,CSF_KIND)
            call alloc(c1,g2%o%layout,S_TYPE)
            c1 = b*g2%o%cdata
            call put(c1,g1,g2%o%type)
          case (CDP_KIND,CDF_KIND)
            call alloc(c1,g2%o%layout,D_TYPE,g2%o%scope)
            c1 = b*g2%o%cdata
            call put(c1,g1,g2%o%type)
          end select
        case (RS_KIND,RD_KIND)
          call take(r1,g1,g1_type)
          if (g2%o%type == EMPTY_KIND) then
            r1 = a*r1
          elseif (g2%o%type == g1_type) then
            r1 = a*r1 + b*g2%o%rdata
          else
            call my(g2,g2_tmp)
            call take(r2,g2_tmp,g1_type)
            r1 = a*r1 + b*r2
            call put(r2,g2_tmp,g1_type)
            call glean(thy(g2_tmp))
          end if
          call put(r1,g1,g1_type)
        case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
          call take(c1,g1,g1_type)
          if (g2%o%type == EMPTY_KIND) then
            c1 = a*c1
          elseif (g2%o%type == g1_type) then
            c1 = a*c1 + b*g2%o%cdata
          else
            call my(g2,g2_tmp)
            call take(c2,g2_tmp,g1_type)
            c1 = a*c1 + b*c2
            call put(c2,g2_tmp,g1_type)
            call glean(thy(g2_tmp))
          end if
          call put(c1,g1,g1_type)
        end select

100     if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(g1))
        call glean(thy(g2))

        if (error("Exit grid_mod::saxpby_grid")) continue

      end subroutine

      subroutine merge_grid_density_grid(g)
!doc$ subroutine merge_grid_density(g)
        type(grid_obj) :: g
!       requires: g%o%type be the same on all processors.
!       modifies: g
!       effects: Merges g mesh values up one scoping level and places them in g.
!       errors: g%o%type /= RS_KIND, CSP_KIND, or CSF_KIND. Passes errors.

!cod$
        integer :: g_type
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:,:,:), pointer :: c1, c2

        call my(g)

        nullify( r1, r2, c1, c2 )

        select case (g%o%scope)
        case (KGROUP)
          g_type = g%o%type
          select case (g_type)
          case (RS_KIND)
            call take(r1,g,g_type)
            call alloc(r2,g%o%layout,S_TYPE)
            call xcomm_allreduce(XKGROUP,MPI_SUM,r1,r2) ; if (error()) goto 100
            call put(r2,g,g_type)
          case (CSP_KIND,CSF_KIND)
            call take(c1,g,g_type)
            call alloc(c2,g%o%layout,S_TYPE)
            call xcomm_allreduce(XKGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
            call put(c2,g,g_type)
          case default
            if (error(.true.,"ERROR: grid type is not RS_KIND, CSP_KIND, or CSF_KIND")) goto 100
          end select
          g%o%scope = SGROUP
        case (SGROUP)
          g_type = g%o%type
          select case (g_type)
          case (RS_KIND)
            call take(r1,g,g_type)
            call alloc(r2,g%o%layout,S_TYPE)
            call xcomm_allreduce(XSGROUP,MPI_SUM,r1,r2) ; if (error()) goto 100
            call put(r2,g,g_type)
          case (CSP_KIND,CSF_KIND)
            call take(c1,g,g_type)
            call alloc(c2,g%o%layout,S_TYPE)
            call xcomm_allreduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
            call put(c2,g,g_type)
          case default
            if (error(.true.,"ERROR: grid type is not RS_KIND, CSP_KIND, or CSF_KIND")) goto 100
          end select
          g%o%scope = CONFIG
        end select

        if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

100     call glean(thy(g))

        if (error("Exit grid_mod::merge_grid_density_grid")) continue

      end subroutine

      function xsum_grid(g) result(gs)
!doc$ function xsum(g) result(gs)
        type(grid_obj) :: g
        type(grid_obj) :: gs
!       effects: Returns the sum of g data across processes with the same rank in XSGROUP.
!       errors: g%o%scope /= SGROUP. Passes errors.

!cod$
        integer :: g_type
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:,:,:), pointer :: c1, c2

        call my(g)

        nullify( r1, r2, c1, c2 )

        call my(grid(g%o%layout,g%o%scope),gs)

        if (error(g%o%scope /= SGROUP,"ERROR: grid scope is not SGROUP")) goto 100

        g_type = g%o%type
        select case (g_type)
        case (RS_KIND)
          call take(r1,g,g_type)
          call alloc(r2,g%o%layout,S_TYPE)
          call xcomm_rank_allreduce(XSGROUP,MPI_SUM,r1,r2) ; if (error()) goto 100
          call put(r1,g,g_type)
          call put(r2,gs,g_type)
        case (RD_KIND)
          call take(r1,g,g_type)
          call alloc(r2,g%o%layout,D_TYPE,g%o%scope)
          call xcomm_rank_allreduce(XSGROUP,MPI_SUM,r1,r2) ; if (error()) goto 100
          call put(r1,g,g_type)
          call put(r2,gs,g_type)
        case (CSP_KIND,CSF_KIND)
          call take(c1,g,g_type)
          call alloc(c2,g%o%layout,S_TYPE)
          call xcomm_rank_allreduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
          call put(c1,g,g_type)
          call put(c2,gs,g_type)
        case (CDP_KIND,CDF_KIND)
          call take(c1,g,g_type)
          call alloc(c2,g%o%layout,D_TYPE,g%o%scope)
          call xcomm_rank_allreduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
          call put(c1,g,g_type)
          call put(c2,gs,g_type)
        end select

100     if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call bequeath(thy(gs))

        call glean(thy(g))

        if (error("Exit grid_mod::xsum_grid")) continue

      end function

      function xdif_grid(g) result(gd)
!doc$ function xdif(g) result(gd)
        type(grid_obj) :: g
        type(grid_obj) :: gd
!       effects: Returns the difference of g data across processes with the same rank in XSGROUP.
!                Returns an empty grid if mpi_nsgroups() = 1.
!       errors: g%o%scope /= SGROUP. Passes errors.

!cod$
        integer :: g_type, msg, nsg
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:,:,:), pointer :: c1, c2

        call my(g)

        nullify( r1, r2, c1, c2 )

        nsg = mpi_nsgroups()

        call my(grid(g%o%layout,g%o%scope),gd)

        if (error(g%o%scope /= SGROUP,"ERROR: grid scope is not SGROUP")) goto 100

        select case (nsg)
        case (2)
          msg = mpi_mysgroup()
          g_type = g%o%type
          select case (g_type)
          case (RS_KIND)
            call take(r1,g,g_type)
            call alloc(r2,g%o%layout,S_TYPE)
            if (msg == 2) r1 = -1.0_double*r1
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,r1,r2) ; if (error()) goto 100
            if (msg == 2) r1 = -1.0_double*r1
            call put(r1,g,g_type)
            call put(r2,gd,g_type)
          case (RD_KIND)
            call take(r1,g,g_type)
            call alloc(r2,g%o%layout,D_TYPE,g%o%scope)
            if (msg == 2) r1 = -1.0_double*r1
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,r1,r2) ; if (error()) goto 100
            if (msg == 2) r1 = -1.0_double*r1
            call put(r1,g,g_type)
            call put(r2,gd,g_type)
          case (CSP_KIND,CSF_KIND)
            call take(c1,g,g_type)
            call alloc(c2,g%o%layout,S_TYPE)
            if (msg == 2) c1 = -1.0_double*c1
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
            if (msg == 2) c1 = -1.0_double*c1
            call put(c1,g,g_type)
            call put(c2,gd,g_type)
          case (CDP_KIND,CDF_KIND)
            call take(c1,g,g_type)
            call alloc(c2,g%o%layout,D_TYPE,g%o%scope)
            if (msg == 2) c1 = -1.0_double*c1
            call xcomm_rank_allreduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
            if (msg == 2) c1 = -1.0_double*c1
            call put(c1,g,g_type)
            call put(c2,gd,g_type)
          end select
        end select

100     if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call bequeath(thy(gd))

        call glean(thy(g))

        if (error("Exit grid_mod::xdif_grid")) continue

      end function

      subroutine sgroup_to_kgroup_grid(g)
!doc$ subroutine sgroup_to_kgroup(g)
        type(grid_obj) :: g
!       requires: g%o%scope = SGROUP.
!       modifies: g
!       effects: Transforms g data from SGROUP to KGROUP distribution.
!       errors: Passes errors.

!cod$
        call my(g)

        call own_i(g)
        select case (g%o%type)
        case (RS_KIND)
          g%o%scope = KGROUP
        case (RD_KIND)
          call transform(g,RS_KIND)         ; if (error()) goto 100
          g%o%scope = KGROUP
          call transform(g,RD_KIND)         ; if (error()) goto 100
        case (CSP_KIND)
          g%o%scope = KGROUP
        case (CDP_KIND)
          call transform(g,CSP_KIND)        ; if (error()) goto 100
          g%o%scope = KGROUP
          call transform(g,CDP_KIND)        ; if (error()) goto 100
        case (CSF_KIND)
          g%o%scope = KGROUP
        case (CDF_KIND)
          call transform(g,CSF_KIND)        ; if (error()) goto 100
          g%o%scope = KGROUP
          call transform(g,CDF_KIND)        ; if (error()) goto 100
        end select

100     call glean(thy(g))

        if (error("Exit grid_mod::sgroup_to_kgroup_grid")) continue

      end subroutine

      subroutine config_to_kgroup_grid(g)
!doc$ subroutine config_to_kgroup(g)
        type(grid_obj) :: g
!       requires: g%o%scope = CONFIG.
!       modifies: g
!       effects: Transforms g data from CONFIG to KGROUP distribution.
!       errors: Passes errors.

!cod$
        call my(g)

        if (g%o%scope == KGROUP) goto 100

        call own_i(g)
        select case (g%o%type)
        case (RS_KIND)
          g%o%scope = KGROUP
        case (RD_KIND)
          call transform(g,RS_KIND)         ; if (error()) goto 100
          g%o%scope = KGROUP
          call transform(g,RD_KIND)         ; if (error()) goto 100
        case (CSP_KIND)
          g%o%scope = KGROUP
        case (CDP_KIND)
          call transform(g,CSP_KIND)        ; if (error()) goto 100
          g%o%scope = KGROUP
          call transform(g,CDP_KIND)        ; if (error()) goto 100
        case (CSF_KIND)
          g%o%scope = KGROUP
        case (CDF_KIND)
          call transform(g,CSF_KIND)        ; if (error()) goto 100
          g%o%scope = KGROUP
          call transform(g,CDF_KIND)        ; if (error()) goto 100
        end select

100     call glean(thy(g))

        if (error("Exit grid_mod::config_to_kgroup_grid")) continue

      end subroutine

      subroutine read_restart_grid(g,tag,restf)
!doc$ subroutine read_restart(g,tag,restf)
        type(grid_obj) :: g
        character(line_len) :: tag
        type(tagio_obj) :: restf
!       requires: All CONFIG processes be present. g%o%scope = SGROUP.
!       modifies: g
!       effects: Reads coefficients from restf and places them in g.
!       errors: Problems reading restf. Passes errors.

!cod$
        character(1) :: tios
        integer :: i1, i2, i3, j1, j2, j3, msg, nsg
        integer, dimension(3) :: ci, gi
        integer(long) :: dsize, iosl, ndata, s4
        real(double) :: g_cutoff, r_cutoff
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(layout_obj) :: g_layout, r_layout

        call my(g)
        call my(restf)

        nullify( c1 )
        nullify( c2 )

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        call own_i(g)
        call empty_i(g)

        call my(g%o%layout,g_layout)
        g_cutoff = x_cutoff(g_layout)

        ! Open the block
        if (i_access(restf)) tios = findfirsttag(restf,trim(tag))
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: grid block was not found")) goto 200
        if (i_access(restf)) call openblock(restf)

        ! Read the cutoff
        if (i_access(restf)) tios = findfirsttag(restf,"CUTOFF")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: CUTOFF tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_double
          ndata = 1
          call readf(r_cutoff,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,r_cutoff)

        ! Set the file pointer at the beginning of the data
        if (i_access(restf)) tios = findfirsttag(restf,"DATA")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: DATA tag was not found")) goto 100

        ! Allocate space for the data
        if (r_cutoff == g_cutoff) then
          call alloc(c1,g_layout,S_TYPE)
          select case (nsg)
          case (2)
            call alloc(c2,g_layout,S_TYPE)
          end select
        else
          call my(layout(x_lattice(g_layout),r_cutoff),r_layout)
          call alloc(c1,r_layout,S_TYPE)
          select case (nsg)
          case (2)
            call alloc(c2,r_layout,S_TYPE)
          end select
        end if

        ! Read the spin group 1 data
        if (i_access(restf)) then
          dsize = sizeof_double
          ndata = 2*size(c1)
          call readf(c1,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,c1)

        ! Read the spin group 2 data
        select case (nsg)
        case (2)
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 2*size(c2)
            call readf(c2,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,c2)
          select case (msg)
          case (2)
            c1 = c2
          end select
          deallocate( c2 )
        end select

        ! Transfer the data to g and convert to CDF_KIND
        if (r_cutoff == g_cutoff) then
          call put(c1,g,CSF_KIND)
        else
          call alloc(c2,g_layout,S_TYPE)
          if (r_cutoff < g_cutoff) then
            c2 = (0.0_double,0.0_double)
            do i3 = 1,size(c1,3)
            do i2 = 1,size(c1,2)
            do i1 = 1,size(c1,1)
              ci(1) = i1
              ci(2) = i2
              ci(3) = i3
              gi = a2fl(ci,r_layout,S_TYPE)
              ci = fl2a(gi,g_layout,S_TYPE)
              j1 = ci(1)
              j2 = ci(2)
              j3 = ci(3)
              c2(j1,j2,j3) = c1(i1,i2,i3)
            end do
            end do
            end do
            call put(c2,g,CSF_KIND)
          elseif (r_cutoff > g_cutoff) then
            do i3 = 1,size(c2,3)
            do i2 = 1,size(c2,2)
            do i1 = 1,size(c2,1)
              ci(1) = i1
              ci(2) = i2
              ci(3) = i3
              gi = a2fl(ci,g_layout,S_TYPE)
              ci = fl2a(gi,r_layout,S_TYPE)
              j1 = ci(1)
              j2 = ci(2)
              j3 = ci(3)
              c2(i1,i2,i3) = c1(j1,j2,j3)
            end do
            end do
            end do
            call put(c2,g,CSF_KIND)
            call filter(g)
          end if
        end if
        call transform(g,CDF_KIND)

        g%o%g = x_ghost()

        ! Close the block
100     if (i_access(restf)) call closeblock(restf)

200     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(g_layout))
        if (r_cutoff /= g_cutoff) call glean(thy(r_layout))

        call glean(thy(g))
        call glean(thy(restf))

        if (error("Exit grid_mod::read_restart_grid")) continue

      end subroutine

      subroutine write_restart_grid(g,tag,nrestf)
!doc$ subroutine write_restart(g,tag,nrestf)
        type(grid_obj) :: g
        character(line_len) :: tag
        type(tagio_obj) :: nrestf
!       requires: All CONFIG processes be present. g%o%scope = SGROUP.
!       modifies: nrestf
!       effects: Writes g coefficients to nrestf.
!       errors: g empty.

!cod$
        integer :: msg, nsg
        integer(long) :: dsize, iosl, ndata, s4
        complex(double), dimension(:,:,:), pointer :: c1, c2
        type(grid_obj) :: g_t

        call my(g)
        call my(nrestf)

        nullify( c1 )
        nullify( c2 )

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        ! Check that there is data on the grid
        if (error(g%o%type == EMPTY_KIND,"ERROR: grid is empty")) goto 200

        ! Start the block
        if (i_access(nrestf)) call startblock(nrestf,trim(tag))

        ! Write the cutoff (needed to restart with a different den_cutoff)
        if (i_access(nrestf)) then
          call writetag(nrestf,"CUTOFF")
          dsize = sizeof_double
          ndata = 1
          call writef(x_cutoff(g%o%layout),dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Make a copy of the grid
        call my(g,g_t)

        ! Get all the data on each process
        call take(c1,g_t,CSF_KIND)

        ! Write the spin group 1 data (residing on the FILE_SCOPE rank 0 process)
        if (i_access(nrestf)) then
          call writetag(nrestf,"DATA")
          dsize = sizeof_double
          ndata = 2*size(c1)
          call writef(c1,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Copy the spin group 2 data to the FILE_SCOPE rank 0 process and write
        select case (nsg)
        case (2)
          select case (msg)
          case (1)
            c1 = (0.0_double,0.0_double)
          end select
          call alloc(c2,g_t%o%layout,S_TYPE)
          call xcomm_reduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) goto 100
          if (i_access(nrestf)) then
            dsize = sizeof_double
            ndata = 2*size(c2)
            call writef(c2,dsize,ndata,x_tagfd(nrestf),iosl)
          end if
        end select

        ! End the block
        if (i_access(nrestf)) call endblock(nrestf)

100     if (associated( c1 )) deallocate( c1 )
        if (associated( c2 )) deallocate( c2 )

        call glean(thy(g_t))

200     call glean(thy(g))
        call glean(thy(nrestf))

        if (error("Exit grid_mod::write_restart_grid")) continue

      end subroutine



      subroutine write_to_file_grid(g,filename,file_format,kind,spin,kpt)
!doc$ subroutine write_to_file(g,filename,file_format,kind,spin,kpt)
        type(grid_obj) :: g
        character(line_len)         :: filename
        integer,intent(in)         :: file_format
        integer,intent(in)         :: kind
        integer,optional           :: spin
        real(double), optional     :: kpt(3)
!       modifies: nothing
!       requires: All CONFIG processes be present.
!       requires: scope must be SGROUP
!       requires: if present, kind must be serial (i.e. not distributed)
!       requires: if present, spin must be either 1 or 2.
!       effects:  Writes the mesh to a file called "filename" of type file_format
!       effects:    according to kind.  If given, will print out the appropriate spin
!       errors:   g empty. kind is distributed. filename must open a valid file. 
!                   passes errors.
!cod$
        ! Local Vars
        integer :: final_kind

        call my(g)

        ! Check that there is data on the grid
        if (error(g%o%type == EMPTY_KIND,"ERROR: grid is empty")) goto 200

        ! Check that the scope is SGROUP
        if (error(g%o%scope /= SGROUP, "ERROR: grid scope is not SGROUP")) goto 200

        ! Verify that kind is serial type, if not change it and throw a warning
        select case(kind)
        case (RS_KIND,CSP_KIND,CSF_KIND)
           ! Any of these are fine so press on...
           final_kind = kind
        case(RD_KIND)
           !  Reset to serial if kind was distributed
           final_kind = RS_KIND
           call warn("Warning: call to grid::write_to_file with kind=RD_KIND ")
           call warn("            setting    kind=RS_KIND")
        case(CDP_KIND)
           final_kind = CSP_KIND
           call warn("Warning: call to grid::write_to_file with kind=CDP_KIND")
           call warn("            setting    kind=CSF_KIND")
        case(CDF_KIND)
           final_kind = CSF_KIND
           call warn("Warning: call to grid::write_to_file with kind=CDF_KIND ")
           call warn("            setting    kind=CSF_KIND")
        case default
           if (error(.true.,'Error: unrecognized kind')) goto 200
        end select
        

        select case(file_format)
        case (MATLAB)
          if (present(spin)) then
            if (present(kpt)) then
              call write_grid_to_matlab_i(g,filename,final_kind,spin=spin,kpt=kpt)
            else
              call write_grid_to_matlab_i(g,filename,final_kind,spin=spin)
            endif
          else
            if (present(kpt)) then
              call write_grid_to_matlab_i(g,filename,final_kind,kpt=kpt)
            else
              call write_grid_to_matlab_i(g,filename,final_kind)
            endif
          end if
        case(AMIRA)
          if (present(spin)) then
            if (present(kpt)) then
              call write_grid_to_amreg_i(g,filename,final_kind,spin=spin,kpt=kpt)
            else
              call write_grid_to_amreg_i(g,filename,final_kind,spin=spin)
            endif
          else
            if (present(kpt)) then
              call write_grid_to_amreg_i(g,filename,final_kind,kpt=kpt)
            else
              call write_grid_to_amreg_i(g,filename,final_kind)
            endif
          end if
        case(VTK)
          if (present(spin)) then
            if (present(kpt)) then
              call write_grid_to_vtk_i(g,filename,final_kind,spin=spin,kpt=kpt)
            else
              call write_grid_to_vtk_i(g,filename,final_kind,spin=spin)
            endif
          else
            if (present(kpt)) then
              call write_grid_to_vtk_i(g,filename,final_kind,kpt=kpt)
            else
              call write_grid_to_vtk_i(g,filename,final_kind)
            endif
          end if                 
        case default
           if (error(.true.,"Error: unrecognized file_format")) goto 200
        end select
        

200     call glean(thy(g))

        if (error("Exit grid_mod::write_to_file_grid")) continue

      end subroutine



! private routines

      subroutine own_i(g)
        type(grid_obj) :: g
        type(grid_obj) :: gt
        if (g%ref < g%o%ref) then
          allocate( gt%o )
          gt%o%ref = 0
          call my(g%o%layout,gt%o%layout)
          gt%o%scope = g%o%scope
          gt%o%type = g%o%type
          select case (g%o%type)
          case (EMPTY_KIND)
            continue
          case (RS_KIND)
            call alloc(gt%o%rdata,g%o%layout,S_TYPE)
            gt%o%rdata = g%o%rdata
            nullify( gt%o%cdata )
          case (RD_KIND)
            call alloc(gt%o%rdata,g%o%layout,D_TYPE,g%o%scope)
            gt%o%rdata = g%o%rdata
            nullify( gt%o%cdata )
          case (CSP_KIND,CSF_KIND)
            call alloc(gt%o%cdata,g%o%layout,S_TYPE)
            gt%o%cdata = g%o%cdata
            nullify( gt%o%rdata )
          case (CDP_KIND,CDF_KIND)
            call alloc(gt%o%cdata,g%o%layout,D_TYPE,g%o%scope)
            gt%o%cdata = g%o%cdata
            nullify( gt%o%rdata )
          end select
          gt%o%g = g%o%g
          g%o%ref = g%o%ref - g%ref
          g%o => gt%o
          g%o%ref = g%o%ref + g%ref
        end if
      end subroutine

      subroutine empty_i(g)
        type(grid_obj) :: g
        select case (g%o%type)
        case (EMPTY_KIND)
          continue
        case (RD_KIND,RS_KIND)
          deallocate( g%o%rdata )
          nullify( g%o%rdata )
        case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
          deallocate( g%o%cdata )
          nullify( g%o%cdata )
        end select
        g%o%type = EMPTY_KIND
      end subroutine

      recursive subroutine transfer_grid_i(ly,sc,otype,ordata,ocdata,itype,irdata,icdata)
        type(layout_obj) :: ly
        integer, intent(in) :: sc
        integer :: otype, itype
        real(double), dimension(:,:,:), pointer :: ordata, irdata
        complex(double), dimension(:,:,:), pointer :: ocdata, icdata

!       requires: otype, ordata, ocdata be consistent with fields of type grid
!                 itype, irdata, icdata be consistent with fields of type grid
!       errors: if otype is empty, but passes errors
!       modifies: otype, ordata, ocdata, itype, irdata, icdata, performs the pointer swap and
!                 conversion required to transfer idata to odata.

        real(double), dimension(:,:,:), pointer :: rtmp
        complex(double), dimension(:,:,:), pointer :: ctmp

        integer :: i1, i2, i3

        call my(ly)

        if (error(itype == EMPTY_KIND,"ERROR: empty data source")) goto 100
        if ( (otype == EMPTY_KIND) .or. (otype == itype) ) then
          otype = itype
          itype = EMPTY_KIND
          ordata => irdata
          nullify( irdata )
          ocdata => icdata
          nullify( icdata )
          goto 100
        elseif ( (itype == CSP_KIND) .and. ( (otype == RS_KIND) .or. (otype == RD_KIND) ) ) then
          itype = RS_KIND
          call alloc(irdata,ly,S_TYPE)
          irdata = real(icdata,double)
          deallocate( icdata )
        elseif ( (itype == CDP_KIND) .and. ( (otype == RS_KIND) .or. (otype == RD_KIND) ) ) then
          itype = RD_KIND
          call alloc(irdata,ly,D_TYPE,sc)
          irdata = real(icdata,double)
          deallocate( icdata )
        elseif ( (itype == RS_KIND) .and. ( (otype == RD_KIND) .or. (otype == CDP_KIND) .or. (otype == CDF_KIND) ) ) then
          rtmp => irdata
          itype = RD_KIND
          call alloc(irdata,ly,D_TYPE,sc)
          call scatter(ly,sc,rtmp,irdata)
          deallocate( rtmp )
        elseif ( ( (itype == CSP_KIND) .or. (itype == CSF_KIND) ) .and. &
                 ( (otype == RD_KIND) .or. (otype == CDP_KIND) .or. (otype == CDF_KIND) ) ) then
          ctmp => icdata
          if (itype == CSP_KIND) itype = CDP_KIND
          if (itype == CSF_KIND) itype = CDF_KIND
          call alloc(icdata,ly,D_TYPE,sc)
          call scatter(ly,sc,ctmp,icdata)
          deallocate( ctmp )
        elseif ( (itype == RS_KIND) .and. ( (otype == CSP_KIND) .or. (otype == CSF_KIND) ) ) then
          itype = CSP_KIND
          call alloc(icdata,ly,S_TYPE)
          do i3 = 1, size(icdata, 3)
          do i2 = 1, size(icdata, 2)
          do i1 = 1, size(icdata, 1)
            icdata(i1,i2,i3) = cmplx(irdata(i1,i2,i3),0,double)
          end do
          end do
          end do
          deallocate( irdata )
        elseif ( (itype == RD_KIND) .and. &
                 ( (otype == CSP_KIND) .or. (otype == CSF_KIND) .or. (otype == CDP_KIND) .or. (otype == CDF_KIND) ) ) then
          itype = CDP_KIND
          call alloc(icdata,ly,D_TYPE,sc)
          icdata = cmplx(irdata,0,double)
          deallocate( irdata )
        elseif ( (itype == CSP_KIND) .and. (otype == CSF_KIND) ) then
          call warn("serial FFT")
          call fft_serial(ly,icdata,R_TO_Q)
          itype = CSF_KIND
        elseif ( (itype == CSF_KIND) .and. ( (otype == CSP_KIND) .or. (otype == RS_KIND) ) ) then
          call warn("serial FFT")
          call fft_serial(ly,icdata,Q_TO_R)
          itype = CSP_KIND
        elseif ( (itype == CDP_KIND) .and. ( (otype == CSF_KIND) .or. (otype == CDF_KIND) ) ) then
          call fft_distributed_sc(ly,sc,icdata,R_TO_Q)
          itype = CDF_KIND
        elseif ( (itype == CDF_KIND) .and. &
                 ( (otype == CSP_KIND) .or. (otype == CDP_KIND) .or. (otype == RS_KIND) .or. (otype == RD_KIND) ) ) then
          call fft_distributed_sc(ly,sc,icdata,Q_TO_R)
          itype = CDP_KIND
        elseif ( (itype == RD_KIND) .and. (otype == RS_KIND) ) then
          itype = RS_KIND
          rtmp => irdata; 
          call alloc(irdata,ly,S_TYPE)
          call gather(ly,sc,rtmp,irdata)
          deallocate( rtmp )
        elseif ( ( (itype == CDP_KIND) .or. (itype == CDF_KIND) ) .and. ( (otype == CSP_KIND) .or. (otype == CSF_KIND) ) ) then
          if (itype == CDP_KIND) itype = CSP_KIND
          if (itype == CDF_KIND) itype = CSF_KIND
          ctmp => icdata
          call alloc(icdata,ly,S_TYPE)
          call gather(ly,sc,ctmp,icdata)
          deallocate( ctmp )
        end if
        if (error()) goto 100

        call transfer_grid_i(ly,sc,otype,ordata,ocdata,itype,irdata,icdata)

100     call glean(thy(ly))

        if (error("Exit grid_mod::transfer_grid_i")) continue

      end subroutine
      

      !** Writes a grid to matlab format
      subroutine write_grid_to_matlab_i(g,fname,kind,spin,kpt)
        type(grid_obj)             :: g
        character(line_len)         :: fname
        integer, intent(in)        :: kind
        integer, optional          :: spin
        real(double), optional     :: kpt(3)

        ! Local Vars
        type(grid_obj)           :: g2
        type(file_obj)           :: f
        type(layout_obj)         :: lay
        real(double), pointer    :: x(:,:,:), y(:,:,:), z(:,:,:)
        real(double), pointer    :: rtmp(:,:,:)
        real(double), pointer    :: rtmp2(:,:,:)
        complex(double), pointer :: ctmp(:,:,:)
        complex(double), pointer :: ctmp2(:,:,:)
        integer :: px,py,pz,nx,ny,nz,ix,iy,iz
        integer :: ios, dims(3), nsg, msg

        logical ::debug, valid_spin

        !----------------------------------------------------------------------

        debug = .false.

        call my(g)
        call my(g,g2)
        call my(x_layout(g),lay)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        nullify(rtmp, rtmp2)
        nullify(ctmp, ctmp2)
        nullify(x,y,z)


        if (debug) call warn("write_grid_to_matlab_i:: starting")


        !** Pull the mesh out of the grid object
        select case(kind)
        case(RS_KIND)
           call take(rtmp,g2,kind) ; if (error()) goto 100
           px = size(rtmp,1)
           py = size(rtmp,2)
           pz = size(rtmp,3)

           ! If spin dependent, load rtmp with the appropriate spin
           select case (nsg)
           case (2)
              if (present(spin)) then
                 valid_spin = ((spin == 1) .or. (spin ==2))
                 if (error(.not.valid_spin,"Error: invalid spin (not 1 or 2)")) goto 100
                 ! set other component of spin to zero.
                 if (spin /= msg) then
                    rtmp = 0.0_double
                 end if
                 ! collect appropriate spin to proc 0
                 call alloc(rtmp2,g2%o%layout,S_TYPE)
                 call xcomm_reduce(XSGROUP,MPI_SUM,rtmp,rtmp2) ; if (error()) goto 100
                 rtmp = rtmp2
                 deallocate(rtmp2)
              end if
           end select

        case(CSP_KIND,CSF_KIND)
           call take(ctmp,g2,kind) ; if (error()) goto 100
           px = size(ctmp,1)
           py = size(ctmp,2)
           pz = size(ctmp,3)

           ! If spin dependent, load ctmp with the appropriate spin
           select case (nsg)
           case (2)
              if (present(spin)) then
                 valid_spin = ((spin == 1) .or. (spin ==2))
                 if (error(.not.valid_spin,"Error: invalid spin (not 1 or 2)")) goto 100
                 ! set other component of spin to zero.
                 if (spin /= msg) then
                    ctmp = 0.0_double
                 end if
                 ! collect appropriate spin to proc 0
                 call alloc(ctmp2,g2%o%layout,S_TYPE)
                 call xcomm_reduce(XSGROUP,MPI_SUM,ctmp,ctmp2) ; if (error()) goto 100
                 ctmp = ctmp2
                 deallocate(ctmp2)
              end if
           end select

        case default
           if (error(.true.,"Error: invalid kind")) goto 100
        end select


        if (debug) call warn("write_grid_to_matlab_i:: calling alloc")

!        call alloc(x,x_layout(g2),RS_KIND, S_TYPE) ; if (error()) goto 100
!        call alloc(y,x_layout(g2),RS_KIND, S_TYPE) ; if (error()) goto 100
!        call alloc(z,x_layout(g2),RS_KIND, S_TYPE) ; if (error()) goto 100

        if (debug) call warn("write_grid_to_matlab_i:: calling mesh")

        !** place the Spatial vector components on x, y, and z
        select case (kind)
        case(RS_KIND,CSP_KIND)
           call mesh(x,y,z,x_layout(g2),S_TYPE) ; if (error()) GOTO 100
        case(CSF_KIND)
           call fmesh(x,y,z,x_layout(g2),S_TYPE) ; if (error()) GOTO 100
        case default
           if (error(.true.,"Error: invalid kind")) goto 100
        end select


        nx = size(x,1)
        ny = size(x,2)
        nz = size(x,3)

        if (error((nx/=px) .or. (ny/=py) .or. (nz/=pz),"ERROR: n /= p")) goto 100

        if (debug) call warn("write_grid_to_matlab_i:: writing coords")


        !** Open the file
        call my(file(trim(fname)), f)
        if (i_access(f)) open(unit=x_unit(f), file=x_name(f), &
             form='formatted',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100


        if (i_access(f)) then

           write(x_unit(f),'(a)') '# Socorro grid written in Matlab .mat format'
           select case (kind)
           case(RS_KIND,CSP_KIND)
              write(x_unit(f),'(a)') '# Real space representation (units are Bohr radii)'   
           case(CSF_KIND)
              write(x_unit(f),'(a)') '# Fourier space representation (units are inverse Bohr radii)'
           case default
              if (error(.true.,"Error: invalid kind")) goto 100
           end select

           !** write out the spin component if spin polarized
           select case (nsg)
           case (1)
              write(x_unit(f),'(a)') '# Not spin dependent'
              write(x_unit(f),'(a)') '' 
           case (2)
              write(x_unit(f),'(a)') '# Spin polarized'
              write(x_unit(f),'(a)') ''
              write(x_unit(f),'(a)') '# name: spin'
              write(x_unit(f),'(a)') '# type: scalar'
              if (present(spin)) then
                 write(x_unit(f),'(I2)') spin
              else
                 write(x_unit(f),'(I2)') 1
              end if
              write(x_unit(f),'(a)') ''
           end select

           !** write out the k vector if present
           if (present(kpt)) then
              write(x_unit(f),'(a)') '# name: kpt'
              write(x_unit(f),'(a)') '# type: matrix'
              write(x_unit(f),'(a)') '# rows: 1'
              write(x_unit(f),'(a)') '# columns: 3'
              write(x_unit(f),'(3f14.8)') kpt
              write(x_unit(f),'(a)') ''
           end if

           !** write out the dimensions of the cell
           write(x_unit(f),'(a)') '# name: dims'
           write(x_unit(f),'(a)') '# type: matrix'
           write(x_unit(f),'(a)') '# rows: 1'
           write(x_unit(f),'(a)') '# columns: 3'
           write(x_unit(f),'(3i5)') nx, ny, nz
           write(x_unit(f),'(a)') ''

           !** write out the lattice vectors
           write(x_unit(f),'(a)') '# name: latvec'
           write(x_unit(f),'(a)') '# type: matrix'
           write(x_unit(f),'(a)') '# rows: 3'
           write(x_unit(f),'(a)') '# columns: 3'
           write(x_unit(f),'(3e18.8)') x_lattice_vector(x_lattice(lay),1)
           write(x_unit(f),'(3e18.8)') x_lattice_vector(x_lattice(lay),2)
           write(x_unit(f),'(3e18.8)') x_lattice_vector(x_lattice(lay),3)
           write(x_unit(f),'(a)') ''

           !** write out the data
           write(x_unit(f),'(a)') '# Data is in the following format:'
           write(x_unit(f),'(a)') '#     x              y             z         Re(data(xyz))     Im(data(xyz))'
                
           write(x_unit(f),'(a)') ''
           write(x_unit(f),'(a)') '# name: data'
           write(x_unit(f),'(a)') '# type: matrix'
           write(x_unit(f),'(a, i10)') '# rows: ', nx*ny*nz

           select case(kind)
           case (RS_KIND)   
              write(x_unit(f),'(a)') '# columns: 4'
              do ix=1,nx
                 do iy=1,ny
                    do iz=1,nz
                       write(x_unit(f),'(3e14.4,e18.8)') x(ix,iy,iz),&
                                                   y(ix,iy,iz),&
                                                   z(ix,iy,iz),&
                                                   rtmp(ix,iy,iz)
                    end do
                 end do
              end do
           case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
              write(x_unit(f),'(a)') '# columns: 5'
              do ix=1,nx
                 do iy=1,ny
                    do iz=1,nz
                       write(x_unit(f),'(3e14.4,2e18.8)') x(ix,iy,iz),&
                                                   y(ix,iy,iz),&
                                                   z(ix,iy,iz),&
                                                   real(ctmp(ix,iy,iz)),&
                                                   imag(ctmp(ix,iy,iz))
                    end do
                 end do
              end do
           end select

        end if

        if (debug) call warn("write_grid_to_matlab_i:: deallocating")

        !** Clean up the mess
        if (i_access(f)) close (x_unit(f))
        if (kind == RS_KIND) then
           deallocate(rtmp)
           nullify(ctmp,ctmp2)
        elseif ((kind == CSP_KIND) .or. (kind == CSF_KIND)) then
           deallocate(ctmp)
           nullify(rtmp,rtmp2)
        end if
        
        if (associated(x)) deallocate(x)
        if (associated(y)) deallocate(y)
        if (associated(z)) deallocate(z)

        nullify(x,y,z)

        call glean(thy(lay))
200     call glean(thy(f))
        call glean(thy(g2))
        call glean(thy(g))
100     if (error("Exit grid_mod:: write_grid_to_matlab_i: ")) continue

        if (debug) call warn("write_grid_to_matlab_i:: exiting")

      end subroutine




      !** Writes a grid to vtk format
      subroutine write_grid_to_vtk_i(g,fname,kind,spin,kpt)
        type(grid_obj)             :: g
        character(line_len)         :: fname
        integer, intent(in)        :: kind
        integer, optional          :: spin
        real(double), optional     :: kpt(3)

        ! Local Vars
        type(grid_obj)           :: g2
        type(file_obj)           :: f
        type(layout_obj)         :: Lay
        real(double), pointer    :: x(:,:,:), y(:,:,:), z(:,:,:)
        real(double), pointer    :: rtmp(:,:,:)
        real(double), pointer    :: rtmp2(:,:,:)
        complex(double), pointer :: ctmp(:,:,:)
        complex(double), pointer :: ctmp2(:,:,:)
        integer :: px,py,pz,nx,ny,nz,ix,iy,iz
        integer :: ios, dims(3), nsg, msg
        character(250)           :: title
        logical ::debug, valid_spin

        !----------------------------------------------------------------------

        debug = .false.

        call my(g)
        call my(g,g2)
        call my(x_layout(g),lay)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        nullify(rtmp, rtmp2)
        nullify(ctmp, ctmp2)

        nullify(x,y,z)

        if (debug) call warn("write_grid_to_vtk_i:: starting")


        !** Pull the mesh out of the grid object
        select case(kind)
        case(RS_KIND)
           call take(rtmp,g2,kind) ; if (error()) goto 100
           px = size(rtmp,1)
           py = size(rtmp,2)
           pz = size(rtmp,3)

           ! If spin dependent, load rtmp with the appropriate spin
           select case (nsg)
           case (2)
              if (present(spin)) then
                 valid_spin = ((spin == 1) .or. (spin ==2))
                 if (error(.not.valid_spin,"Error: invalid spin (not 1 or 2)")) goto 100
                 ! set other component of spin to zero.
                 if (spin /= msg) then
                    rtmp = 0.0_double
                 end if
                 ! collect appropriate spin to proc 0
                 call alloc(rtmp2,g2%o%layout,S_TYPE)
                 call xcomm_reduce(XSGROUP,MPI_SUM,rtmp,rtmp2) ; if (error()) goto 100
                 rtmp = rtmp2
                 deallocate(rtmp2)
              end if
           end select

        case(CSP_KIND,CSF_KIND)
           call take(ctmp,g2,kind) ; if (error()) goto 100
           px = size(ctmp,1)
           py = size(ctmp,2)
           pz = size(ctmp,3)

           ! If spin dependent, load ctmp with the appropriate spin
           select case (nsg)
           case (2)
              if (present(spin)) then
                 valid_spin = ((spin == 1) .or. (spin ==2))
                 if (error(.not.valid_spin,"Error: invalid spin (not 1 or 2)")) goto 100
                 ! set other component of spin to zero.
                 if (spin /= msg) then
                    ctmp = 0.0_double
                 end if
                 ! collect appropriate spin to proc 0
                 call alloc(ctmp2,g2%o%layout,S_TYPE)
                 call xcomm_reduce(XSGROUP,MPI_SUM,ctmp,ctmp2) ; if (error()) goto 100
                 ctmp = ctmp2
                 deallocate(ctmp2)
              end if
           end select

        case default
           if (error(.true.,"Error: invalid kind")) goto 100
        end select


        if (debug) call warn("write_grid_to_vtk_i:: calling mesh")

        !** place the Spatial vector components on x, y, and z
        select case (kind)
        case(RS_KIND,CSP_KIND)
           call mesh(x,y,z,x_layout(g2),S_TYPE) ; if (error()) GOTO 100
        case(CSF_KIND)
           call fmesh(x,y,z,x_layout(g2),S_TYPE) ; if (error()) GOTO 100
        case default
           if (error(.true.,"Error: invalid kind")) goto 100
        end select


        nx = size(x,1)
        ny = size(x,2)
        nz = size(x,3)

        if (error((nx/=px) .or. (ny/=py) .or. (nz/=pz),"ERROR: n /= p")) goto 100

        if (debug) call warn("write_grid_to_vtk_i:: writing coords")


        !** Open the file
        call my(file(trim(fname)), f)
        if (i_access(f)) open(unit=x_unit(f), file=x_name(f), &
             form='formatted',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100


        if (i_access(f)) then

           write(x_unit(f),'(a)') '# vtk DataFile Version 2.0'

           select case (kind)
           case(RS_KIND,CSP_KIND)
             select case (nsg)
             case(1)
               if (present(kpt)) then
                 write(x_unit(f),'(a,3f6.4)') &
                 'Socorro grid, Real space representation, Not spin dependent, kpt = ', kpt   
               else
                 write(x_unit(f),'(a)') &
                 'Socorro grid, Real space representation, Not spin dependent'
               end if
             case(2)
               if (present(kpt)) then
                 write(x_unit(f),'(a,3f6.4)') &
                 'Socorro grid, Real space representation, Spin polarized, kpt = ', kpt   
               else
                 write(x_unit(f),'(a)') &
                 'Socorro grid, Real space representation, Spin polarized'
               end if
             end select
           case(CSF_KIND)
             select case (nsg)
             case(1)
               if (present(kpt)) then
                 write(x_unit(f),'(a,3f6.4)') &
                 'Socorro grid, Fourier space representation, Not spin dependent, kpt = ', kpt   
               else
                 write(x_unit(f),'(a)') &
                 'Socorro grid, Fourier space representation, Not spin dependent'
               end if
             case(2)
               if (present(kpt)) then
                 write(x_unit(f),'(a,3f6.4)') &
                 'Socorro grid, Fourier space representation, Spin polarized, kpt = ', kpt   
               else
                 write(x_unit(f),'(a)') &
                 'Socorro grid, Fourier space representation, Spin polarized'
               end if
             end select
           case default
              if (error(.true.,"Error: invalid kind")) goto 100
           end select

           write(x_unit(f),'(a)') 'ASCII'
           write(x_unit(f),'(a)') 'DATASET STRUCTURED_GRID'
           write(x_unit(f),'(a,3i5)') 'DIMENSIONS ', nx, ny, nz
           write(x_unit(f),'(a,i10,a)') 'POINTS ', nx*ny*nz, ' FLOAT' 


           do ix=1,nx
             do iy=1,ny
               do iz=1,nz
                 write(x_unit(f),'(3f18.10)') x(ix,iy,iz),y(ix,iy,iz),z(ix,iy,iz)
               end do
             end do
           end do
           
           write(x_unit(f),'(a,i10,a)') 'POINT_DATA ', nx*ny*nz
           write(x_unit(f),*) ' '
           write(x_unit(f),'(a,a,a)') 'SCALARS data FLOAT'
           write(x_unit(f),'(a)') 'LOOKUP_TABLE default'
           

           if (debug) call warn("write_grid_to_vtk_i:: writing density")

           do ix=1,nx
             do iy=1,ny
               do iz=1,nz
                 select case(kind)
                 case (RS_KIND,RD_KIND)   
                   write(x_unit(f),'(f18.10)') rtmp(ix,iy,iz)
                 case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
                   write(x_unit(f),'(2f18.10)') ctmp(ix,iy,iz)
                 end select
               end do
             end do
           end do

        end if

        if (debug) call warn("write_grid_to_vtk_i:: deallocating")

        !** Clean up the mess
        if (i_access(f)) close (x_unit(f))
        if (kind == RS_KIND) then
           deallocate(rtmp)
           nullify(ctmp,ctmp2)
        elseif ((kind == CSP_KIND) .or. (kind == CSF_KIND)) then
           deallocate(ctmp)
           nullify(rtmp,rtmp2)
        end if

        if (associated(x)) deallocate(x)
        if (associated(y)) deallocate(y)
        if (associated(z)) deallocate(z)

        nullify(x,y,z)        

        call glean(thy(lay))
200     call glean(thy(f))
        call glean(thy(g2))
        call glean(thy(g))
100     if (error("Exit grid_mod:: write_grid_to_vtk_i: ")) continue

        if (debug) call warn("write_grid_to_vtk_i:: exiting")

      end subroutine

      !** Writes a grid to amira format
      subroutine write_grid_to_amreg_i(g,fname,kind,spin,kpt)
        type(grid_obj)             :: g
        character(line_len)         :: fname
        integer, intent(in)        :: kind
        integer, optional          :: spin
        real(double), optional     :: kpt(3)

        ! Local Vars
        type(grid_obj)           :: g2
        type(file_obj)           :: f
        type(layout_obj)         :: lay
        real(double)             :: xmin,xmax,ymin,ymax,zmin,zmax
        real(double), pointer    :: x(:,:,:), y(:,:,:), z(:,:,:)
        real(double), pointer    :: rtmp(:,:,:)
        real(double), pointer    :: rtmp2(:,:,:)
        complex(double), pointer :: ctmp(:,:,:)
        complex(double), pointer :: ctmp2(:,:,:)
        integer :: px,py,pz,nx,ny,nz,ix,iy,iz
        integer :: ios, dims(3), nsg, msg
        integer :: numx, numy, numz
        integer :: stepx = 1
        integer :: stepy = 1
        integer :: stepz = 1

        logical ::debug, valid_spin

        !----------------------------------------------------------------------

        debug = .false.

        call my(g)
        call my(g,g2)
        call my(x_layout(g),lay)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        nullify(rtmp, rtmp2)
        nullify(ctmp, ctmp2)
        nullify(x,y,z)

        if (debug) call warn("write_grid_to_amreg_i:: starting")


        !** Pull the mesh out of the grid object
        select case(kind)
        case(RS_KIND)
           call take(rtmp,g2,kind) ; if (error()) goto 100
           px = size(rtmp,1)
           py = size(rtmp,2)
           pz = size(rtmp,3)

           ! If spin dependent, load rtmp with the appropriate spin
           select case (nsg)
           case (2)
              if (present(spin)) then
                 valid_spin = ((spin == 1) .or. (spin ==2))
                 if (error(.not.valid_spin,"Error: invalid spin (not 1 or 2)")) goto 100
                 ! set other component of spin to zero.
                 if (spin /= msg) then
                    rtmp = 0.0_double
                 end if
                 ! collect appropriate spin to proc 0
                 call alloc(rtmp2,g2%o%layout,S_TYPE)
                 call xcomm_reduce(XSGROUP,MPI_SUM,rtmp,rtmp2) ; if (error()) goto 100
                 rtmp = rtmp2
                 deallocate(rtmp2)
              end if
           end select

        case(CSP_KIND,CSF_KIND)
           call take(ctmp,g2,kind) ; if (error()) goto 100
           px = size(ctmp,1)
           py = size(ctmp,2)
           pz = size(ctmp,3)

           ! If spin dependent, load ctmp with the appropriate spin
           select case (nsg)
           case (2)
              if (present(spin)) then
                 valid_spin = ((spin == 1) .or. (spin ==2))
                 if (error(.not.valid_spin,"Error: invalid spin (not 1 or 2)")) goto 100
                 ! set other component of spin to zero.
                 if (spin /= msg) then
                    ctmp = 0.0_double
                 end if
                 ! collect appropriate spin to proc 0
                 call alloc(ctmp2,g2%o%layout,S_TYPE)
                 call xcomm_reduce(XSGROUP,MPI_SUM,ctmp,ctmp2) ; if (error()) goto 100
                 ctmp = ctmp2
                 deallocate(ctmp2)
              end if
           end select

        case default
           if (error(.true.,"Error: invalid kind")) goto 100
        end select


        if (debug) call warn("write_grid_to_amreg_i:: calling mesh")

        !** place the Spatial vector components on x, y, and z
        select case (kind)
        case(RS_KIND,CSP_KIND)
           call mesh(x,y,z,x_layout(g2),S_TYPE) ; if (error()) GOTO 100
        case(CSF_KIND)
           call fmesh(x,y,z,x_layout(g2),S_TYPE) ; if (error()) GOTO 100
        case default
           if (error(.true.,"Error: invalid kind")) goto 100
        end select


        nx = size(x,1)
        ny = size(x,2)
        nz = size(x,3)

        if (error((nx/=px) .or. (ny/=py) .or. (nz/=pz),"ERROR: n /= p")) goto 100

        if (debug) call warn("write_grid_to_amreg_i:: writing coords")


        !** Open the file
        call my(file(trim(fname)), f)
        if (i_access(f)) open(unit=x_unit(f), file=x_name(f), &
             form='formatted',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100


        if (i_access(f)) then

           xmin = x(1,1,1)
           xmax = x(1,1,1)
           ymin = y(1,1,1)
           ymax = y(1,1,1)
           zmin = z(1,1,1)
           zmax = z(1,1,1)

           do ix=1,nx,stepx
             do iy=1,ny,stepy
               do iz=1,nz,stepz
                 if (x(ix,iy,iz) < xmin) xmin = x(ix,iy,iz)
                 if (xmax < x(ix,iy,iz)) xmax = x(ix,iy,iz)
                 
                 if (y(ix,iy,iz) < ymin) ymin = y(ix,iy,iz)
                 if (ymax < y(ix,iy,iz)) ymax = y(ix,iy,iz)
                 
                 if (z(ix,iy,iz) < zmin) zmin = z(ix,iy,iz)
                 if (zmax < z(ix,iy,iz)) zmax = z(ix,iy,iz)
               end do
             end do
           end do

           
           numx = ceiling(real(nx)/stepx)
           numy = ceiling(real(ny)/stepy)
           numz = ceiling(real(nz)/stepz)
           


           write(x_unit(f),'(a)')       '# AmiraMesh ASCII 1.0'
           write(x_unit(f),'(a,3i5)')   'define Lattice ', numx,numy,numz
!          write(x_unit(f),'(a,3i5)')   'define Lattice ', nz, ny, nx
           write(x_unit(f),'(a)')       'Parameters {'
           write(x_unit(f),'(a)')       '        CoordType "uniform"'
           write(x_unit(f),'(a,6e12.5)')'        BoundingBox ',xmin,xmax,ymin,ymax,zmin,zmax
           write(x_unit(f),'(a)')       '}'
           write(x_unit(f),'(a)')       'Lattice { float ScalarField } = @1'
           write(x_unit(f),'(a)')       '@1'

           if (debug) call warn("write_grid_to_am_i:: writing density")


           do iz=1,nz,stepx
              do iy=1,ny,stepy
                 do ix=1,nx,stepz
                    select case(kind)
                    case (RS_KIND,RD_KIND)   
                       write(x_unit(f),'(e18.8)') rtmp(ix,iy,iz)
                    case (CSP_KIND,CDP_KIND,CSF_KIND,CDF_KIND)
                       write(x_unit(f),'(2e18.8)') ctmp(ix,iy,iz)
                    end select
                 end do
              end do
           end do


        end if

        if (debug) call warn("write_grid_to_amreg_i:: deallocating")

        !** Clean up the mess
        if (i_access(f)) close (x_unit(f))
        if (kind == RS_KIND) then
           deallocate(rtmp)
           nullify(ctmp,ctmp2)
        elseif ((kind == CSP_KIND) .or. (kind == CSF_KIND)) then
           deallocate(ctmp)
           nullify(rtmp,rtmp2)
        end if
        
        if (associated(x)) deallocate(x)
        if (associated(y)) deallocate(y)
        if (associated(z)) deallocate(z)

        nullify(x,y,z)        

        call glean(thy(lay))
200     call glean(thy(f))
        call glean(thy(g2))
        call glean(thy(g))
100     if (error("Exit grid_mod:: write_grid_to_amreg_i: ")) continue

        if (debug) call warn("write_grid_to_amreg_i:: exiting")

      end subroutine





      end module

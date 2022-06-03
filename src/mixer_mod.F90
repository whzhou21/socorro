!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module mixer_mod
!doc$ module mixer_mod

!     One datatype is available here: type(mixer_obj)

!     mixer_mod applies mixing to field PAW atomic quantities in order to speed up the convergence of
!     a self-consistent electronic structure calculation.

      use kind_mod
      use arg_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use ghost_mod
      use math_mod
      use grid_mod
      use layout_mod
      use lattice_mod
      use atomic_density_mod
      use atomic_potential_mod
      use multivector_mod
      use dyad_kpoint_mod
      use dyad_mod

!cod$
      implicit none
      private

      ! atomic_status
      integer, parameter :: ATOMIC_TRIVIAL     = 0
      integer, parameter :: ATOMIC_NON_TRIVIAL = 1

      ! method
      integer, parameter :: SIMPLE   = 1
      integer, parameter :: PULAY    = 2

      ! damping_type and metric_type
      integer, parameter :: NONE   = 1
      integer, parameter :: UNITY  = 1
      integer, parameter :: KERKER = 2

      ! profiles
      integer, parameter :: FIXED   = 1
      integer, parameter :: RAMPED  = 2
      integer, parameter :: STEPPED = 3

      type :: mixer_rep
        type(ghost) :: g
        integer :: ref
        type(ghost) :: g_layout                               ! layout ghost
        type(file_obj) :: f                                   ! mixing report
        logical :: report                                     ! switch to print a report
        integer :: mxh_profile                                ! max history profile type
        integer :: wgt_profile                                ! weight profile type
        integer :: dgc_profile                                ! damping g crossover profile type
        integer :: max_steps                                  ! maximum number of steps (max_steps in config_sc_mod)
        integer :: step                                       ! current step
        integer :: atomic_status                              ! status of the atomic object argument
        integer :: method                                     ! mixing method
        integer :: damping_type                               ! damping type in field mixing
        integer :: metric_type                                ! metric type in field mixing (Pulay mixing)
        integer :: hs                                         ! history slot (Pulay method)
        integer :: nh                                         ! current number of histories (Pulay method)
        integer, dimension(:), pointer :: mxh                 ! step-dependent maximum histories (Pulay method)
        integer, dimension(:,:), pointer :: fmap              ! mapping between 3d and 1d data structures in field mixing
        real(double) :: metric_gc                             ! g crossover for Kerker metrics (Pulay method)
        real(double), dimension(:), pointer :: wgt            ! step-dependent weights
        real(double), dimension(:), pointer :: dgc            ! step-dependent g crossover (Kerker damping)
        real(double), dimension(:), pointer :: g2i            ! inverse (g**2) used (Kerker damping)
        real(double), dimension(:), pointer :: metric         ! g-dependent metrics (Pulay method)
        real(double), dimension(:,:), pointer :: mm           ! mixing matrix (Pulay method)
        complex(double), dimension(:), pointer :: fc_val      ! current values in field mixing
        complex(double), dimension(:), pointer :: fp_val      ! previous values in field mixing (Pulay method)
        complex(double), dimension(:), pointer :: fp_res      ! previous residuals in field mixing (Pulay method)
        complex(double), dimension(:,:), pointer :: fd_val    ! differences in values in field mixing (Pulay method)
        complex(double), dimension(:,:), pointer :: fd_res    ! differences in residuals in field mixing (Pulay method)
        complex(double), dimension(:), pointer :: ac_val      ! current values in atomic mixing
        complex(double), dimension(:), pointer :: ap_val      ! previous values in atomic mixing (Pulay method)
        complex(double), dimension(:), pointer :: ap_res      ! previous residuals in atomic mixing (Pulay method)
        complex(double), dimension(:,:), pointer :: ad_val    ! differences in values in atomic mixing (Pulay method)
        complex(double), dimension(:,:), pointer :: ad_res    ! differences in residuals in atomic mixing (Pulay method)
        real(double) :: wgt_dp                                ! mixing weight for the dyadic potential
        type(dyad_obj), pointer :: dp_val                     ! previous dyadic potential
      end type

      type, public :: mixer_obj
        private
        integer :: ref
        type(mixer_rep), pointer :: o
      end type

!doc$
      public :: mixer
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: mix
      public :: diary
      public :: x_ref
      public :: x_ghost

!cod$
      interface mixer
        module procedure constructor_mx
      end interface
      interface update
        module procedure update_mx
      end interface
      interface my
        module procedure my_mx, my_new_mx
      end interface
      interface thy
        module procedure thy_mx
      end interface
      interface glean
        module procedure glean_mx
      end interface
      interface bequeath
        module procedure bequeath_mx
      end interface
      interface assignment(=)
        module procedure assign_mx
      end interface
      interface mix
        module procedure mix_mx
      end interface
      interface mix
        module procedure mix_mx_dp
      end interface
      interface diary
        module procedure diary_mx
      end interface
      interface x_ref
        module procedure mx_ref
      end interface
      interface x_ghost
        module procedure mx_ghost
      end interface

      contains

! public routines

      function constructor_mx(fq,ad,ap) result(mx)
!doc$ function mixer(fq,ad.ap) result(mx)
        type(grid_obj) :: fq
        type(atomic_density_obj), optional :: ad
        type(atomic_potential_obj), optional :: ap
        type(mixer_obj) :: mx
!       requires: fq data be in D_TYPE distribution.
!                 Only one of ad and ap be present.
!       effects: Constructs a new mx.

!cod$
        logical :: found
        character(line_len) :: tag
        integer :: h1, h2, h3, i1, i2, i3, ig, ios, is, ng, nl, ramp_start, ramp_stop, s12, s23
        real(double) :: gci, gcf, ramp_rate, wti, wtf
        real(double), dimension(:,:,:), pointer :: g2, g2i
        type(layout_obj) :: lay

        call my(fq)
        if (present(ad)) call my(ad)
        if (present(ap)) call my(ap)

        mx%ref = 0
        allocate( mx%o )
        mx%o%ref = 0
        mx%o%g = x_ghost()

        nullify( g2, g2i )

        call my(x_layout(fq),lay)
        mx%o%g_layout = x_ghost(lay)

        ! nullify object arrays
        nullify( mx%o%mxh )
        nullify( mx%o%fmap )
        nullify( mx%o%wgt )
        nullify( mx%o%dgc )
        nullify( mx%o%g2i )
        nullify( mx%o%metric )
        nullify( mx%o%mm )
        nullify( mx%o%fc_val )
        nullify( mx%o%fp_val )
        nullify( mx%o%fp_res )
        nullify( mx%o%fd_val )
        nullify( mx%o%fd_res )
        nullify( mx%o%ac_val )
        nullify( mx%o%ap_val )
        nullify( mx%o%ap_res )
        nullify( mx%o%ad_val )
        nullify( mx%o%ad_res )
        nullify( mx%o%dp_val )

        ! determine the atomic status
        if (present(ad)) then
          mx%o%atomic_status = ATOMIC_NON_TRIVIAL
          if (trivial(ad)) mx%o%atomic_status = ATOMIC_TRIVIAL
        end if
        if (present(ap)) then
          mx%o%atomic_status = ATOMIC_NON_TRIVIAL
          if (trivial(ap)) mx%o%atomic_status = ATOMIC_TRIVIAL
        end if

        ! determine the reporting status
        call arglc("mix_report",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("off")
          mx%o%report = .false.
        case ("on")
          mx%o%report = .true.
        case default
          if (error(.true.,"WARNING: mix_report not recognized")) continue
        end select

        ! determine the size of field arrays
        call fdel(g2,lay,D_TYPE,SGROUP)
        ng = 0
        do i3 = 1,size(g2,3)
        do i2 = 1,size(g2,2)
        do i1 = 1,size(g2,1)
          if (g2(i1,i2,i3) > x_cutoff(lay)) cycle
          ng = ng + 1
        end do
        end do
        end do

        ! determine the 3d-1d mapping for field arrays
        allocate( mx%o%fmap(3,ng) )
        ig = 0
        do i3 = 1,size(g2,3)
        do i2 = 1,size(g2,2)
        do i1 = 1,size(g2,1)
          if (g2(i1,i2,i3) > x_cutoff(lay)) cycle
          ig = ig + 1
          mx%o%fmap(:,ig) = (/i1,i2,i3/)
        end do
        end do
        end do

        ! get the mixing method
        call arglc("mix_method",tag,found)
        if (.not.found) tag = "pulay"
        select case (trim(tag))
        case ("simple")
          mx%o%method = SIMPLE
        case ("pulay")
          mx%o%method = PULAY
        case default
          if (error(.true.,"ERROR: mix_method not recognized")) goto 100
        end select

        ! get the weight profile type
        call arglc("mix_wgt_profile",tag,found)
        if (.not.found) tag = "fixed"
        select case (trim(tag))
        case ("fixed")
          mx%o%wgt_profile = FIXED
        case ("ramped")
          mx%o%wgt_profile = RAMPED
        case default
          if (error(.true.,"ERROR: mix_wgt_profile not recognized")) goto 100
        end select

        ! get the damping type
        call arglc("mix_damping_type",tag,found)
        if (.not.found) tag = "none"
        select case (trim(tag))
        case ("none")
          mx%o%damping_type = NONE
        case ("kerker")
          mx%o%damping_type = KERKER
        case default
          if (error(.true.,"ERROR: mix_damping_type not recognized")) goto 100
        end select

        ! get the kerker damping profile
        select case (mx%o%damping_type)
        case (KERKER)
          call arglc("mix_dgc_profile",tag,found)
          if (.not.found) tag = "fixed"
          select case (trim(tag))
          case ("fixed")
            mx%o%dgc_profile = FIXED
          case ("ramped")
            mx%o%dgc_profile = RAMPED
          case default
            if (error(.true.,"ERROR: mix_dgc_profile not recognized")) goto 100
          end select
        end select

        ! get the history profile type
        select case (mx%o%method)
        case (PULAY)
          call arglc("mix_mxh_profile",tag,found)
          if (.not.found) tag = "fixed"
          select case (trim(tag))
          case ("fixed")
            mx%o%mxh_profile = FIXED
          case ("stepped")
            mx%o%mxh_profile = STEPPED
          case default
            if (error(.true.,"ERROR: mix_mxh_profile not recognized")) goto 100
          end select
        end select

        ! get the maximum number of mixing steps
        call arg("config_steps",mx%o%max_steps,found)
        if (.not.found) mx%o%max_steps = 40 ! Note: this default must be consistent with the one in config_sc_mod.f90

        ! set weight profile values
        allocate( mx%o%wgt(mx%o%max_steps) )
        call arg("mix_weight",wti,found)
          if (.not.found) wti = 0.8_double
          if (error(wti <= 0.0_double,"ERROR: mix_weight <= 0")) goto 100
          if (error(wti > 1.0_double,"ERROR: mix_weight > 1")) goto 100
        select case (mx%o%wgt_profile)
        case (FIXED)
          mx%o%wgt = wti
        case (RAMPED)
          call arg("mix_wgt_final",wtf,found)
            if (.not.found) wtf = wti + 0.2_double
            if (error(wtf <= 0.0_double,"ERROR: mix_wgt_final <= 0")) goto 100
            if (error(wtf > 1.0_double,"ERROR: mix_wgt_final > 1")) goto 100
            if (error(wtf == wti,"ERROR: mix_wgt_final = mix_weight")) goto 100
          call arg("mix_wgt_ramp_start",ramp_start,found)
            if (.not.found) ramp_start = 8
            if (error(ramp_start < 2,"ERROR: mix_wgt_ramp_start < 2")) goto 100
          call arg("mix_wgt_ramp_stop",ramp_stop,found)
            if (.not.found) ramp_stop = ramp_start + 6
            if (error(ramp_stop < ramp_start,"ERROR: mix_wgt_ramp_stop < mix_wgt_ramp_start")) goto 100
          ramp_rate = (wtf - wti)/real((ramp_stop - ramp_start + 1),double)
          do is = 1,mx%o%max_steps
            if (is < ramp_start) then
              mx%o%wgt(is) = wti
            elseif (is < ramp_stop) then
              mx%o%wgt(is) = mx%o%wgt(is-1) + ramp_rate
            else
              mx%o%wgt(is) = wtf
            end if
          end do
        end select

        ! set kerker damping values
        select case (mx%o%damping_type)
        case (KERKER)
          call arg("mix_damping_gc",gci,found)
            if (.not.found) gci = 0.8_double
            if (error(gci <= 0.0_double,"ERROR: mix_damping_gc <= 0")) goto 100
          allocate( mx%o%dgc(mx%o%max_steps) )
          select case (mx%o%dgc_profile)
          case (FIXED)
            mx%o%dgc = gci
          case (RAMPED)
            call arg("mix_dgc_final",gcf,found)
              if (.not.found) gcf = gci + 0.2_double
              if (error(gcf <= 0.0_double,"ERROR: mix_dgc_final <= 0")) goto 100
              if (error(gcf == gci,"ERROR: mix_dgc_final = mix_damping_gc")) goto 100
            call arg("mix_dgc_ramp_start",ramp_start,found)
              if (.not.found) ramp_start = 8
              if (error(ramp_start < 2,"ERROR: mix_dgc_ramp_start < 2")) goto 100
            call arg("mix_dgc_ramp_stop",ramp_stop,found)
              if (.not.found) ramp_stop = ramp_start + 6
              if (error(ramp_stop < ramp_start,"ERROR: mix_dgc_ramp_stop < mix_dgc_ramp_start")) goto 100
            ramp_rate = (gcf - gci)/real((ramp_stop - ramp_start + 1),double)
            do is = 1,mx%o%max_steps
              if (is < ramp_start) then
                mx%o%dgc(is) = gci
              elseif (is < ramp_stop) then
                mx%o%dgc(is) = mx%o%dgc(is-1) + ramp_rate
              else
                mx%o%dgc(is) = gcf
              end if
            end do
          end select
        end select

        ! determine mixing weight for the dyad potential
        call arg("mix_weight_dp",mx%o%wgt_dp,found)
        if (.not.found) mx%o%wgt_dp = 0.8_double
        if (error(mx%o%wgt_dp <= 0.0_double,"ERROR: mix_weight_dp <= 0")) goto 100
        if (error(mx%o%wgt_dp > 1.0_double,"ERROR: mix_weight_dp > 1")) goto 100

        ! set the inverse g^2 array
        select case (mx%o%damping_type)
        case (KERKER)
          call fdelinv(g2i,lay,D_TYPE,SGROUP)
          allocate( mx%o%g2i(ng) )
          ig = 0
          do i3 = 1,size(g2,3)
          do i2 = 1,size(g2,2)
          do i1 = 1,size(g2,1)
            if (g2(i1,i2,i3) > x_cutoff(lay)) cycle
            ig = ig + 1
            mx%o%g2i(ig) = g2i(i1,i2,i3)
          end do
          end do
          end do
        end select

        ! get history profile values
        select case (mx%o%method)
        case (PULAY)
          call arg("mix_history",h1,found)
            if (.not.found) h1 = 5
            if (error(h1 < 1,"ERROR: mix_history < 1")) goto 100
            if (error(h1 > 20,"ERROR: mix_history > 20")) goto 100
            if (error(h1 > mx%o%max_steps,"ERROR: mix_history > max_steps")) goto 100
          allocate( mx%o%mxh(mx%o%max_steps) )
          select case (mx%o%mxh_profile)
          case (FIXED)
            mx%o%mxh = h1
          case (STEPPED)
            call arg("mix_mxh_levels",nl,found)
              if (error(.not.found,"ERROR: mix_mxh_levels not found")) goto 100
              if (error(nl < 2,"ERROR: mix_mxh_levels < 2")) goto 100
              if (error(nl > 3,"ERROR: mix_mxh_levels > 3")) goto 100
            call arg("mix_mxh_step_1-2",s12,found)
              if (error(.not.found,"ERROR: mix_mxh_step_1-2 not found")) goto 100
              if (error(s12 < (h1 + 1),"ERROR: mix_mxh_step_1-2 < mix_history + 1")) goto 100
              if (error(s12 > (mx%o%max_steps - 1),"ERROR: mix_mxh_step_1-2 out of bounds")) goto 100
            call arg("mix_mxh_2",h2,found)
              if (error(.not.found,"ERROR: mix_mxh_2 not found")) goto 100
              if (error(h2 < 1,"ERROR: mix_mxh_2 < 1")) goto 100
              if (error(h2 > 20,"ERROR: mix_mxh_2 > 20")) goto 100
              if (error(h2 == h1,"ERROR: mix_mxh_2 = mix_history")) goto 100
            do is = 1,(s12-1)
              mx%o%mxh(is) = h1
            end do
            if (nl == 2) then
              do is = s12,mx%o%max_steps
                mx%o%mxh(is) = h2
              end do
            else
              call arg("mix_mxh_step_2-3",s23,found)
                if (error(.not.found,"ERROR: mix_mxh_step_2-3 not found")) goto 100
                if (error(s23 <= (s12 + h2),"ERROR: mix_mxh_step_2-3 <= mix_mxh_step_1-2 + mix_mxh_2")) goto 100
                if (error(s23 >= (mx%o%max_steps - 1),"ERROR: mix_mxh_step_2-3 out of bounds")) goto 100
              call arg("mix_mxh_3",h3,found)
                if (error(.not.found,"ERROR: mix_mxh_3 not found")) goto 100
                if (error(h3 < 1,"ERROR: mix_mxh_3 < 1")) goto 100
                if (error(h3 > 20,"ERROR: mix_mxh_3 > 20")) goto 100
                if (error(h3 == h2,"ERROR: mix_mxh_3 = mix_mxh_2")) goto 100
              do is = s12,(s23-1)
                mx%o%mxh(is) = h2
              end do
              do is = s23,mx%o%max_steps
                mx%o%mxh(is) = h3
              end do
            end if
          end select
        end select

        ! get metric information
        select case (mx%o%method)
        case (PULAY)
          call arglc("mix_metric_type",tag,found)
          if (.not.found) tag = "unity"
          select case (trim(tag))
          case ("unity")
            mx%o%metric_type = UNITY
          case ("kerker")
            mx%o%metric_type = KERKER
            call arg("mix_metric_gc",mx%o%metric_gc,found)
            if (.not.found) mx%o%metric_gc = 0.8_double
            if (error(mx%o%metric_gc <= 0.0_double,"ERROR: mix_metric_gc <= 0")) goto 100
          case default
            if (error(.true.,"ERROR: mix_metric_type not recognized")) goto 100
          end select
        end select

        ! construct the field metrics
        select case (mx%o%method)
        case (PULAY)
          allocate( mx%o%metric(ng) )
          select case (mx%o%metric_type)
          case (UNITY)
            mx%o%metric = 1.0_double
          case (KERKER)
            if (.not.associated( g2i )) then
              call fdelinv(g2i,lay,D_TYPE,SGROUP)
            end if
            ig = 0
            do i3 = 1,size(g2,3)
            do i2 = 1,size(g2,2)
            do i1 = 1,size(g2,1)
              if (g2(i1,i2,i3) > x_cutoff(lay)) cycle
              ig = ig + 1
              mx%o%metric(ig) = 1.0_double + g2i(i1,i2,i3)*mx%o%metric_gc**2
            end do
            end do
            end do
          end select
        end select

        ! set the field current values
        call grid_to_c1d_i(fq,mx%o%fmap,mx%o%fc_val) ; if (error()) goto 100

        ! set the atomic current values
        select case (mx%o%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          if (present(ad)) call extract_density(ad,mx%o%ac_val)
          if (present(ap)) call extract_potential(ap,mx%o%ac_val)
        end select


        ! allocate arrays for the Pulay method
        select case (mx%o%method)
        case (PULAY)
          allocate( mx%o%fp_val(size(mx%o%fc_val)) )
          allocate( mx%o%fp_res(size(mx%o%fc_val)) )
          allocate( mx%o%fd_val(size(mx%o%fc_val),mx%o%mxh(1)) )
          allocate( mx%o%fd_res(size(mx%o%fc_val),mx%o%mxh(1)) )
          select case (mx%o%atomic_status)
          case (ATOMIC_NON_TRIVIAL)
            allocate( mx%o%ap_val(size(mx%o%ac_val)) )
            allocate( mx%o%ap_res(size(mx%o%ac_val)) )
            allocate( mx%o%ad_val(size(mx%o%ac_val),mx%o%mxh(1)) )
            allocate( mx%o%ad_res(size(mx%o%ac_val),mx%o%mxh(1)) )
          end select
          allocate( mx%o%mm(mx%o%mxh(1),mx%o%mxh(1)) )
        end select

        ! initialize the step counter
        mx%o%step = 0

        ! initialize the history slot
        mx%o%hs = 0

        ! initialize the current number of histories
        mx%o%nh = 0

        ! open the report file
        if (mx%o%report) then
          call my(file(trim(mixer_report_path)),mx%o%f)
          if (i_access(mx%o%f)) open(x_unit(mx%o%f),file=x_name(mx%o%f),status='unknown',iostat=ios)
        end if

100     if (associated( g2 )) deallocate( g2 )
        if (associated( g2i )) deallocate( g2i )

        call glean(thy(lay))

        call glean(thy(fq))
        if (present(ad)) call glean(thy(ad))
        if (present(ap)) call glean(thy(ap))

        if (error("Exit mixer_mod::constructor_mx")) continue

      end function

      subroutine update_mx(mx,fq,ad,ap)
!doc$ subroutine update(mx,fq,ad,ap)
        type(mixer_obj) :: mx
        type(grid_obj) :: fq
        type(atomic_density_obj), optional :: ad
        type(atomic_potential_obj), optional :: ap
!       requires: fq data be in D_TYPE distribution.
!                 Only one of ad and ap be present.
!       effects: Updates mx.
!       errors: If layout has changed.

!cod$ 
        call my(mx)
        call my(fq)
        if (present(ad)) call my(ad)
        if (present(ap)) call my(ap)

        if (error(mx%o%g_layout /= x_ghost(x_layout(fq)),"ERROR: layout changes are not currently allowed")) goto 100

        call own_mx_i(mx)
        mx%o%g = x_ghost()

        call grid_to_c1d_i(fq,mx%o%fmap,mx%o%fc_val)

        select case (mx%o%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          if (present(ad)) call extract_density(ad,mx%o%ac_val)
          if (present(ap)) call extract_potential(ap,mx%o%ac_val)
        end select

        ! initialize the step counter
        mx%o%step = 0

        ! initialize the history slot
        mx%o%hs = 0

        ! initialize the current number of histories
        mx%o%nh = 0

100     call glean(thy(mx))
        call glean(thy(fq))
        if (present(ad)) call glean(thy(ad))
        if (present(ap)) call glean(thy(ap))

        if (error("Exit mixer_mod::update_mx")) continue

      end subroutine

      subroutine my_mx(mx)
!doc$ subroutine my(mx)
        type(mixer_obj) :: mx

!cod$
        mx%ref = mx%ref + 1
        mx%o%ref = mx%o%ref + 1
      end subroutine

      subroutine my_new_mx(mxi,mx)
!doc$ subroutine my(mxi,mx)
        type(mixer_obj) :: mxi, mx

!cod$
        mx%ref = 1
        mx%o => mxi%o
        mx%o%ref = mx%o%ref + 1
      end subroutine

      function thy_mx(mx) result(mxo)
!doc$ function thy(mx) result(mxo)
        type(mixer_obj) :: mx, mxo

!cod$
        mx%ref = mx%ref - 1
        mx%o%ref = mx%o%ref - 1
        mxo%ref = mx%ref
        mxo%o => mx%o
      end function

      subroutine glean_mx(mx)
!doc$ subroutine glean(mx)
        type(mixer_obj) :: mx

!cod$
        if (mx%o%ref < 1) then
          if (mx%o%report) then
            if (i_access(mx%o%f)) close(x_unit(mx%o%f))
            call glean(thy(mx%o%f))
          end if
          deallocate( mx%o%fmap )
          deallocate( mx%o%wgt )
          select case (mx%o%damping_type)
          case (KERKER)
            deallocate( mx%o%dgc )
            deallocate( mx%o%g2i )
          end select
          deallocate( mx%o%fc_val )
          select case (mx%o%atomic_status)
          case (ATOMIC_NON_TRIVIAL)
            deallocate( mx%o%ac_val )
          end select
          select case (mx%o%method)
          case (PULAY)
            deallocate( mx%o%mxh )
            deallocate( mx%o%fp_val )
            deallocate( mx%o%fp_res )
            deallocate( mx%o%fd_val )
            deallocate( mx%o%fd_res )
            deallocate( mx%o%mm )
            deallocate( mx%o%metric )
            select case (mx%o%atomic_status)
            case (ATOMIC_NON_TRIVIAL)
              deallocate( mx%o%ap_val )
              deallocate( mx%o%ap_res )
              deallocate( mx%o%ad_val )
              deallocate( mx%o%ad_res )
            end select
          end select
          if (associated(mx%o%dp_val)) then
             call glean(thy(mx%o%dp_val))
             deallocate(mx%o%dp_val)
          end if
          deallocate( mx%o )
        end if
      end subroutine

      subroutine bequeath_mx(mx)
!doc$ subroutine bequeath(mx)
        type(mixer_obj) :: mx

!cod$
        continue
      end subroutine

      subroutine assign_mx(mx1,mx2)
!doc$ subroutine assign(mx1,mx2)
        type(mixer_obj), intent(inout) :: mx1
        type(mixer_obj), intent(in) :: mx2

!cod$
        type(mixer_obj) :: mxt
        call my(mx2)
        mxt%o => mx1%o
        mx1%o%ref = mx1%o%ref - mx1%ref
        mx1%o => mx2%o
        mx1%o%ref = mx1%o%ref + mx1%ref
        call glean(mxt)
        call glean(thy(mx2))
      end subroutine

      subroutine mix_mx(mx,fq,ad,ap)
!doc$ subroutine mix(mx,fq,ad,ap)
        type(mixer_obj) :: mx
        type(grid_obj) :: fq
        type(atomic_density_obj), optional :: ad
        type(atomic_potential_obj) , optional :: ap
!       modifies: mx, fq, ad or ap
!       effects: Updates mx.
!       errors: Passes errors.

!cod$ 
        call my(mx)
        call my(fq)
        if (present(ad)) call my(ad)
        if (present(ap)) call my(ap)

        if (error(x_ghost(x_layout(fq)) .ne. mx%o%g_layout,"ERROR: inconsistent layouts")) goto 100

        call own_mx_i(mx)
        mx%o%g = x_ghost()

        mx%o%step = mx%o%step + 1
        if (error(mx%o%step > size(mx%o%wgt),"ERROR: step is out of bounds")) goto 100

        select case (mx%o%method)
        case (SIMPLE)
          if (present(ad)) call simple_mx_i(mx,fq,ad=ad)
          if (present(ap)) call simple_mx_i(mx,fq,ap=ap)
        case (PULAY)
          if (present(ad)) call pulay_mx_i(mx,fq,ad=ad)
          if (present(ap)) call pulay_mx_i(mx,fq,ap=ap)
        end select

100     call glean(thy(mx))
        call glean(thy(fq))
        if (present(ad)) call glean(thy(ad))
        if (present(ap)) call glean(thy(ap))

        if (error("Exit mixer_mod::mix_mx")) continue

      end subroutine

      subroutine mix_mx_dp(mx,dp)
!doc$ subroutine mix(mx,dp)
        type(mixer_obj) :: mx
        type(dyad_obj) :: dp
!       modifies: dp
!       effects: Updates mx.
!       errors: Passes errors.

!cod$
        integer :: ik, nk
        type(dyad_kpoint_obj) :: dk, dko
        type(multivector_obj) :: v, w, xow
        complex(double), parameter :: cp1 = (+1.0_double,0.0_double)
        complex(double), parameter :: cmh = (-0.5_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: e_mat
        call my(mx)
        call my(dp)

        call own_mx_i(mx)
        mx%o%g = x_ghost()

        if (associated(mx%o%dp_val)) then
           nk = x_n_kpoints(dp)
           do ik = 1, nk
              if (is_empty(x_dyad_kpoint(dp,ik))) cycle
              call my(x_dyad_kpoint(dp,ik),dk)
              call my(x_dyad_kpoint(mx%o%dp_val,ik),dko)
              if (error(is_empty(dko),"ERROR: Inconsistency between old and new dyad potentials")) goto 100
              call my(x_v(dk),v)
              call my(x_w(dk),w)
              call my(apply(dko,w),xow)
              allocate(e_mat(x_n_bands(w),x_n_bands(w)))
              call overlap(w,xow,e_mat)
              call transform(cp1,xow,cmh,w,e_mat)
              deallocate(e_mat)
              call combine(mx%o%wgt_dp,v,1.0_double - mx%o%wgt_dp,xow)
              call update(dk,v)
              call update(dp,ik,dk)
              call glean(thy(xow))
              call glean(thy(w))
              call glean(thy(v))
              call glean(thy(dko))
              call glean(thy(dk))
           end do
           mx%o%dp_val = dp
        else
           allocate(mx%o%dp_val)
           call my(dp,mx%o%dp_val)
        end if


100     call glean(thy(mx))
        call glean(thy(dp))

        if (error("Exit mixer_mod::mix_mx_dp")) continue

      end subroutine

      subroutine diary_mx(mx)
!doc$ subroutine diary(mx)
        type(mixer_obj) :: mx
!       effects: Writes mx information to the diary.

!cod$
        integer :: is

        call my(mx)

        if (i_access( diaryfile() )) then
          select case (mx%o%method)
          case (SIMPLE)
            write(x_unit(diaryfile()),'(/,t4,"Simple mixing:")')
          case (PULAY)
            write(x_unit(diaryfile()),'(/,t4,"Pulay mixing:")')
            select case (mx%o%mxh_profile)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t6,"history = ",i0)') mx%o%mxh(1)
            case (STEPPED)
              write(x_unit(diaryfile()),'(/,t6,"initial history = ",i0)') mx%o%mxh(1)
              do is = 2,size(mx%o%mxh)
                if (mx%o%mxh(is) /= mx%o%mxh(is-1)) then
                  write(x_unit(diaryfile()),'(t6,"new history     = ",i0," at step ",i0)') mx%o%mxh(is), is
                end if
              end do
            end select
            select case (mx%o%metric_type)
            case (KERKER)
              write(x_unit(diaryfile()),'(/,t6,"Kerker metric: g crossover = ",f5.3)') mx%o%metric_gc
            end select
          end select
          select case (mx%o%wgt_profile)
          case (FIXED)
            write(x_unit(diaryfile()),'(/,t6,"weight = ",f5.3)') mx%o%wgt(1)
          case (RAMPED)
            write(x_unit(diaryfile()),'(/,t6,"initial weight = ",f5.3)') mx%o%wgt(1)
            do is = 2,size(mx%o%wgt)
              if (mx%o%wgt(is) /= mx%o%wgt(is-1)) then
                write(x_unit(diaryfile()),'(t6,"new weight     = ",f5.3," at step ",i0)') mx%o%wgt(is), is
              end if
            end do
          end select
          select case (mx%o%damping_type)
          case (KERKER)
            select case (mx%o%dgc_profile)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t6,"g crossover = ",f5.3)') mx%o%dgc(1)
            case (RAMPED)
              write(x_unit(diaryfile()),'(/,t6,"initial G crossover = ",f5.3)') mx%o%dgc(1)
              do is = 2,size(mx%o%dgc)
                if (mx%o%dgc(is) /= mx%o%dgc(is-1)) then
                  write(x_unit(diaryfile()),'(t6,"new G crossover     = ",f5.3," at step ",i0)') mx%o%dgc(is), is
                end if
              end do
            end select
          end select
        end if

        call glean(thy(mx))

      end subroutine

      function mx_ref(mx) result(r)
!doc$ function x_ref(mx) result(r)
        type(mixer_obj) :: mx
        integer, dimension(2) :: r
!       effects: Returns mx%ref and mx%o%ref.

!cod$
        r(1) = mx%ref
        r(2) = mx%o%ref
        call glean(mx)
      end function

      function mx_ghost(mx) result(g)
!doc$ function x_ghost(mx) result(g)
        type(mixer_obj) :: mx
        type(ghost) :: g
!       effects: Returns the ghost of mx.

!cod$
        call my(mx)
        g = mx%o%g
        call glean(thy(mx))
      end function

! private routines

      subroutine simple_mx_i(mx,fq,ad,ap)
        type(mixer_obj) :: mx
        type(grid_obj) :: fq
        type(atomic_density_obj), optional :: ad
        type(atomic_potential_obj), optional :: ap

        real(double) :: a_wgt
        real(double), dimension(:), pointer :: f_wgt
        complex(double), dimension(:), pointer :: a_val, f_val
        complex(double), dimension(:), allocatable :: a_res, f_res

        call my(mx)
        call my(fq)
        if (present(ad)) call my(ad)
        if (present(ap)) call my(ap)

        nullify( f_wgt, a_val, f_val )

        if (mx%o%report) call report_parameters_i(mx%o)

        call get_weights_i(mx%o,a_wgt,f_wgt)

        ! field update
        call grid_to_c1d_i(fq,mx%o%fmap,f_val)
        allocate( f_res(size(f_val)) )
        f_res = f_val - mx%o%fc_val
        mx%o%fc_val = mx%o%fc_val + f_wgt*f_res
        call c1d_to_grid_i(mx%o%fc_val,mx%o%fmap,fq)

        ! atomic update
        select case (mx%o%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          if (present(ad)) call extract_density(ad,a_val)
          if (present(ap)) call extract_potential(ap,a_val)
          allocate( a_res(size(a_val)) )
          a_res = a_val - mx%o%ac_val
          mx%o%ac_val = mx%o%ac_val + a_wgt*a_res
          if (present(ad)) call insert_density(mx%o%ac_val,ad)
          if (present(ap)) call insert_potential(mx%o%ac_val,ap)
        end select

        if (mx%o%report) call report_residual_norms_i(mx%o,f_res,a_res)

100     if (associated( f_wgt )) deallocate( f_wgt )
        if (associated( a_val )) deallocate( a_val )
        if (associated( f_val )) deallocate( f_val )
        if (allocated( a_res )) deallocate( a_res )
        if (allocated( f_res )) deallocate( f_res )

        call glean(thy(mx))
        call glean(thy(fq))
        if (present(ad)) call glean(thy(ad))
        if (present(ap)) call glean(thy(ap))

        if (error("Exit mixer_mod::simple_mx_i")) continue

      end subroutine

      subroutine pulay_mx_i(mx,fq,ad,ap)

        type(mixer_obj) :: mx
        type(grid_obj) :: fq
        type(atomic_density_obj), optional :: ad
        type(atomic_potential_obj), optional :: ap

        integer :: is
        real(double) :: a_wgt
        real(double), dimension(:), pointer :: f_wgt
        real(double), dimension(:), allocatable :: alpha, dp, dp_l, dp_res
        real(double), dimension(:,:), allocatable :: amat, amat_inv
        complex(double), dimension(:), pointer :: a_val, f_val
        complex(double), dimension(:), allocatable :: a_res, f_res

        call my(mx)
        call my(fq)
        if (present(ad)) call my(ad)
        if (present(ap)) call my(ap)

        nullify( f_wgt, a_val, f_val )

        if (mx%o%step > 1) then
          if (mx%o%mxh(mx%o%step) /= mx%o%mxh(mx%o%step-1)) call reorder_i(mx%o)
        end if

        if (mx%o%report) call report_parameters_i(mx%o)

        call get_weights_i(mx%o,a_wgt,f_wgt)

        if (mx%o%nh > 0) then
          allocate( alpha(mx%o%nh) )
          allocate( dp(2*mx%o%nh) )
          allocate( dp_l(2*mx%o%nh) )
          allocate( dp_res(mx%o%nh) )
          allocate( amat(mx%o%nh,mx%o%nh) )
          allocate( amat_inv(mx%o%nh,mx%o%nh) )
        end if

        ! field mixing quantities
        call grid_to_c1d_i(fq,mx%o%fmap,f_val)
        allocate( f_res(size(f_val)) )
        f_res = f_val - mx%o%fc_val
        if (mx%o%hs > 0) then
          mx%o%fd_val(:,mx%o%hs) = mx%o%fc_val - mx%o%fp_val
          mx%o%fd_res(:,mx%o%hs) = f_res - mx%o%fp_res
        end if
        do is = 1,mx%o%nh
          dp_l(is) = sum(mx%o%metric*conjg(mx%o%fd_res(:,is))*mx%o%fd_res(:,mx%o%hs))
          dp_l(is+mx%o%nh) = sum(mx%o%metric*conjg(mx%o%fd_res(:,is))*f_res)
        end do

        ! atomic mixing quantities
        select case (mx%o%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          if (present(ad)) call extract_density(ad,a_val)
          if (present(ap)) call extract_potential(ap,a_val)
          allocate( a_res(size(a_val)) )
          a_res = a_val - mx%o%ac_val
          if (mx%o%hs > 0) then
            mx%o%ad_val(:,mx%o%hs) = mx%o%ac_val - mx%o%ap_val
            mx%o%ad_res(:,mx%o%hs) = a_res - mx%o%ap_res
          end if
          do is = 1,mx%o%nh
            dp_l(is) = dp_l(is) + sum(conjg(mx%o%ad_res(:,is))*mx%o%ad_res(:,mx%o%hs))
            dp_l(is+mx%o%nh) = dp_l(is+mx%o%nh) + sum(conjg(mx%o%ad_res(:,is))*a_res)
          end do
        end select

        ! mixing coefficients
        if (mx%o%nh > 0) then
          call allreduce(CONFIG,MPI_SUM,dp_l,dp)
          do is = 1,mx%o%nh
            mx%o%mm(is,mx%o%hs) = dp(is)
            mx%o%mm(mx%o%hs,is) = dp(is)
            dp_res(is) = dp(is+mx%o%nh)
          end do
          amat = mx%o%mm(1:mx%o%nh,1:mx%o%nh)
          call invert_i(amat,amat_inv,mx%o%nh) ; if (error()) goto 100
          alpha = -matmul(amat_inv,dp_res)
          if (mx%o%report) call report_alpha_i(mx%o,alpha)
!          if (mx%o%report) call report_amat_i(mx%o,amat)
!          if (mx%o%report) call report_amat_inv_i(mx%o,amat_inv)
!          if (mx%o%report) call report_dp_res_i(mx%o,dp_res)
        end if

        ! field updates
        mx%o%fp_val = mx%o%fc_val
        mx%o%fp_res = f_res
        mx%o%fc_val = mx%o%fc_val + f_wgt*f_res
        do is = 1,mx%o%nh
          mx%o%fc_val = mx%o%fc_val + alpha(is)*(mx%o%fd_val(:,is) + f_wgt*mx%o%fd_res(:,is))
        end do
        call c1d_to_grid_i(mx%o%fc_val,mx%o%fmap,fq)

        ! atomic updates
        select case (mx%o%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          mx%o%ap_val = mx%o%ac_val
          mx%o%ap_res = a_res
          mx%o%ac_val = mx%o%ac_val + a_wgt*a_res
          do is = 1,mx%o%nh
            mx%o%ac_val = mx%o%ac_val + alpha(is)*(mx%o%ad_val(:,is) + a_wgt*mx%o%ad_res(:,is))
          end do
          if (present(ad)) call insert_density(mx%o%ac_val,ad)
          if (present(ap)) call insert_potential(mx%o%ac_val,ap)
        end select

        if (mx%o%report) call report_residual_norms_i(mx%o,f_res,a_res)

        ! update history slot and number of histories
        mx%o%hs = mx%o%hs + 1
        if (mx%o%hs > mx%o%mxh(mx%o%step)) mx%o%hs = 1
        mx%o%nh = mx%o%nh + 1
        if (mx%o%nh > mx%o%mxh(mx%o%step)) mx%o%nh = mx%o%mxh(mx%o%step)

100     if (associated( f_wgt )) deallocate( f_wgt )
        if (allocated( alpha )) deallocate( alpha )
        if (allocated( dp )) deallocate( dp )
        if (allocated( dp_l )) deallocate( dp_l )
        if (allocated( dp_res )) deallocate( dp_res )
        if (allocated( amat )) deallocate( amat )
        if (allocated( amat_inv )) deallocate( amat_inv )
        if (associated( a_val )) deallocate( a_val )
        if (associated( f_val )) deallocate( f_val )
        if (allocated( a_res )) deallocate( a_res )
        if (allocated( f_res )) deallocate( f_res )

        call glean(thy(mx))
        call glean(thy(fq))
        if (present(ad)) call glean(thy(ad))
        if (present(ap)) call glean(thy(ap))

        if (error("Exit mixer_mod::pulay_mx_i")) continue

      end subroutine

      subroutine own_mx_i(mx)
        type(mixer_obj) :: mx
        type(mixer_obj) :: mxt

        if (mx%ref < mx%o%ref) then
          allocate( mxt%o )
          mxt%o%ref = 0
          nullify( mxt%o%mxh )
          nullify( mxt%o%fmap )
          nullify( mxt%o%wgt )
          nullify( mxt%o%dgc )
          nullify( mxt%o%g2i )
          nullify( mxt%o%metric )
          nullify( mxt%o%mm )
          nullify( mxt%o%fc_val )
          nullify( mxt%o%fp_val )
          nullify( mxt%o%fp_res )
          nullify( mxt%o%fd_val )
          nullify( mxt%o%fd_res )
          nullify( mxt%o%ac_val )
          nullify( mxt%o%ap_val )
          nullify( mxt%o%ap_res )
          nullify( mxt%o%ad_val )
          nullify( mxt%o%ad_res )
          mxt%o%g_layout = mx%o%g_layout
          mxt%o%report = mx%o%report
          if (mx%o%report) call my(mx%o%f,mxt%o%f)
          mxt%o%max_steps = mx%o%max_steps
          mxt%o%step = mx%o%step
          mxt%o%wgt_profile = mx%o%wgt_profile
          mxt%o%damping_type = mx%o%damping_type
          select case (mx%o%damping_type)
          case (KERKER)
            mxt%o%dgc_profile = mx%o%dgc_profile
            allocate( mxt%o%dgc(size(mx%o%dgc)) )  ; mxt%o%dgc = mx%o%dgc
            allocate( mxt%o%g2i(size(mx%o%g2i)) )  ; mxt%o%g2i = mx%o%g2i
          end select
          mxt%o%atomic_status = mx%o%atomic_status
          mxt%o%method = mx%o%method
          allocate( mxt%o%wgt(size(mx%o%wgt)) )                       ; mxt%o%wgt    = mx%o%wgt
          allocate( mxt%o%fmap(size(mx%o%fmap,1),size(mx%o%fmap,2)) ) ; mxt%o%fmap   = mx%o%fmap
          allocate( mxt%o%fc_val(size(mx%o%fc_val)) )                 ; mxt%o%fc_val = mx%o%fc_val
          select case (mx%o%atomic_status)
          case (ATOMIC_NON_TRIVIAL)
            allocate( mxt%o%ac_val(size(mx%o%ac_val)) )   ; mxt%o%ac_val = mx%o%ac_val
          end select
          select case (mxt%o%method)
          case (PULAY)
            mxt%o%hs = mx%o%hs
            mxt%o%nh = mx%o%nh
            mxt%o%mxh_profile = mx%o%mxh_profile
            mxt%o%metric_type = mx%o%metric_type
            select case (mx%o%metric_type)
            case (KERKER)
              mxt%o%metric_gc = mx%o%metric_gc
            end select
            allocate( mxt%o%mxh(size(mx%o%mxh)) )                                ; mxt%o%mxh     = mx%o%mxh
            allocate( mxt%o%metric(size(mx%o%metric)) )                        ; mxt%o%metric = mx%o%metric
            allocate( mxt%o%mm(size(mx%o%mm,1),size(mx%o%mm,2)) )              ; mxt%o%mm     = mx%o%mm
            allocate( mxt%o%fp_val(size(mx%o%fp_val)) )                        ; mxt%o%fp_val = mx%o%fp_val
            allocate( mxt%o%fp_res(size(mx%o%fp_res)) )                        ; mxt%o%fp_res = mx%o%fp_res
            allocate( mxt%o%fd_val(size(mx%o%fd_val,1),size(mx%o%fd_val,2)) )  ; mxt%o%fd_val = mx%o%fd_val
            allocate( mxt%o%fd_res(size(mx%o%fd_res,1),size(mx%o%fd_res,2)) )  ; mxt%o%fd_res = mx%o%fd_res
            select case (mx%o%atomic_status)
            case (ATOMIC_NON_TRIVIAL)
              allocate( mxt%o%ap_val(size(mx%o%ap_val)) )                        ; mxt%o%ap_val = mx%o%ap_val
              allocate( mxt%o%ap_res(size(mx%o%ap_res)) )                        ; mxt%o%ap_res = mx%o%ap_res
              allocate( mxt%o%ad_val(size(mx%o%ad_val,1),size(mx%o%ad_val,2)) )  ; mxt%o%ad_val = mx%o%ad_val
              allocate( mxt%o%ad_res(size(mx%o%ad_res,1),size(mx%o%ad_res,2)) )  ; mxt%o%ad_res = mx%o%ad_res
            end select
          end select
          mxt%o%wgt_dp = mx%o%wgt_dp
          if (associated(mx%o%dp_val)) then
             allocate(mxt%o%dp_val)
             call my(mx%o%dp_val,mxt%o%dp_val)
          else
             nullify(mxt%o%dp_val)
          end if
          mxt%o%g = mx%o%g
          mx%o%ref = mx%o%ref - mx%ref
          mx%o => mxt%o
          mx%o%ref = mx%o%ref + mx%ref
        end if

      end subroutine

      subroutine grid_to_c1d_i(g,m,c1d)
        type(grid_obj) :: g
        integer, dimension(:,:) :: m
        complex(double), dimension(:), pointer :: c1d
!       requires: c1d be nullified or associated.
!       effects: Copies g data to c1d according to m.
!       errors: Passes errors.

        integer :: i
        complex(double), dimension(:,:,:), pointer :: c3d

        call my(g)

        nullify( c3d )

        call take(c3d,g,CDF_KIND) ; if (error()) goto 100
        if (associated( c1d )) deallocate( c1d ) ; allocate( c1d(size(m,2)) )
        do i = 1,size(c1d)
          c1d(i) = c3d(m(1,i),m(2,i),m(3,i))
        end do
        call put(c3d,g,CDF_KIND) ; if (error()) goto 100

100     if (associated( c3d )) deallocate( c3d )

        call glean(thy(g))

        if (error("Exit mixer_mod::grid_to_c1d_i")) continue

      end subroutine

      subroutine c1d_to_grid_i(c1d,m,g)
        complex(double), dimension(:), pointer :: c1d
        integer, dimension(:,:) :: m
        type(grid_obj) :: g
!       modifies: c1d
!       requires: c1d be associated. size(c1d) = size(m,2).
!       effects: Moves c1d data to g according to m.
!       errors: Passes errors.

        integer :: i
        complex(double), dimension(:,:,:), pointer :: c3

        call my(g)

        nullify( c3 )

        call alloc(c3,x_layout(g),D_TYPE,SGROUP)
        c3 = cmplx(0,0,double)
        do i = 1,size(c1d)
          c3(m(1,i),m(2,i),m(3,i)) = c1d(i)
        end do
        call put(c3,g,CDF_KIND) ; if (error()) goto 100

100     if (associated( c3 )) deallocate( c3 )

        call glean(thy(g))

        if (error("Exit mixer_mod::c1d_to_grid_i")) continue

      end subroutine

      subroutine invert_i(a_in,b,n)
        real(double), dimension(:,:), intent(in) :: a_in
        real(double), dimension(:,:), intent(out) :: b
        integer, intent (in) :: n

        real(double), parameter :: condition_no = 1000

        integer :: i, ierr, lda, lwork
        real(double), dimension(size(a_in,1)) :: s
        real(double), dimension(size(a_in,1),size(a_in,2)) :: a, u, vt, s_mat
        real(double), dimension(size(a_in,1)*size(a_in,1)*10) :: work
        real(double) :: tol

        lda = size(a_in,1)
        lwork = 10*lda*lda

        b = a_in
        call dgesvd('A','A',n,n,b(1,1),lda,s(1),u(1,1),lda,vt(1,1),lda,work(1),lwork,ierr)
        if (error(ierr /= 0,"ERROR: dgesdd error condition")) goto 100

        s_mat = 0
        tol = abs(s(1)/condition_no)
        if (tol > machine_precision) then
          do i = 1,n
            if (abs(s(i)) > tol) s_mat(i,i) = 1.0_double/s(i)
          end do
        end If

        a(1:n,1:n) = matmul(s_mat(1:n,1:n),transpose(u(1:n,1:n)))
        b(1:n,1:n) = matmul(transpose(vt(1:n,1:n)),a(1:n,1:n))

100     if (error("Exit mixer_mod::invert_i")) continue

      end subroutine

      subroutine get_weights_i(mxr,a_wgt,f_wgt)
        type(mixer_rep) :: mxr
        real(double) :: a_wgt
        real(double), dimension(:), pointer :: f_wgt

        allocate( f_wgt(size(mxr%fc_val)) )
        select case (mxr%damping_type)
        case (NONE)
          f_wgt = mxr%wgt(mxr%step)
        case (KERKER)
          f_wgt = mxr%wgt(mxr%step)/(1.0_double + mxr%g2i*mxr%dgc(mxr%step)**2)
        end select
        select case (mxr%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          a_wgt = mxr%wgt(mxr%step)
        end select

      end subroutine

      subroutine report_parameters_i(mxr)
        type(mixer_rep) :: mxr

        if (i_access(mxr%f)) then
          if (mxr%step == 1) then
            write(x_unit(mxr%f),'("step ",i0,":")') mxr%step
          else
            write(x_unit(mxr%f),'(/,"step ",i0,":")') mxr%step
          end if
          write(x_unit(mxr%f),'(t4,"mixing parameters:")')
          write(x_unit(mxr%f),'(t6,"weight      = ",f0.3)') mxr%wgt(mxr%step)
          select case (mxr%damping_type)
          case (KERKER)
            write(x_unit(mxr%f),'(t6,"g crossover = ",f0.3)') mxr%dgc(mxr%step)
          end select
          select case (mxr%method)
          case (PULAY)
            write(x_unit(mxr%f),'(t6,"maximum histories    = ",i0)') mxr%mxh(mxr%step)
            write(x_unit(mxr%f),'(t6,"number of histories  = ",i0)') mxr%nh
            write(x_unit(mxr%f),'(t6,"history slot         = ",i0)') mxr%hs
          end select
        end if

      end subroutine

      subroutine report_alpha_i(mxr,alpha)
        type(mixer_rep) :: mxr
        real(double), dimension(:) :: alpha

        if (i_access(mxr%f)) then
          write(x_unit(mxr%f),'(/,t4,"alpha:")')
          write(x_unit(mxr%f),'(es14.2)') alpha
        end if

      end subroutine

      subroutine report_amat_i(mxr,amat)
        type(mixer_rep) :: mxr
        real(double), dimension(:,:) :: amat

        integer :: ir
        if (i_access(mxr%f)) then
          write(x_unit(mxr%f),'(/,t4,"amat:")')
          do ir = 1,mxr%nh
            write(x_unit(mxr%f),'(es14.2)') amat(ir,:)
          end do
        end if

      end subroutine

      subroutine report_amat_inv_i(mxr,amat_inv)
        type(mixer_rep) :: mxr
        real(double), dimension(:,:) :: amat_inv

        integer :: ir
        if (i_access(mxr%f)) then
          write(x_unit(mxr%f),'(/,t4,"amat_inv:")')
          do ir = 1,mxr%nh
            write(x_unit(mxr%f),'(es14.2)') amat_inv(ir,:)
          end do
        end if

      end subroutine

      subroutine report_dp_res_i(mxr,dp_res)
        type(mixer_rep) :: mxr
        real(double), dimension(:) :: dp_res

        integer :: ic
        if (i_access(mxr%f)) then
          write(x_unit(mxr%f),'(/,t4,"dp_res:")')
          do ic = 1,mxr%nh
            write(x_unit(mxr%f),'(es14.2)') dp_res(ic)
          end do
        end if

      end subroutine

      subroutine report_residual_norms_i(mxr,s_res,a_res)
        type(mixer_rep) :: mxr
        complex(double), dimension(:) :: s_res, a_res

        real(double) :: a_rn, s_rn, rn, rn_l

        rn_l = sum(abs(s_res)**2)
        call allreduce(SGROUP,MPI_SUM,rn_l,rn)
        s_rn = sqrt(rn)

        select case (mxr%atomic_status)
        case (ATOMIC_NON_TRIVIAL)
          rn_l = sum(abs(a_res)**2)
          call allreduce(SGROUP,MPI_SUM,rn_l,rn)
          a_rn = sqrt(rn)
        end select

        if (i_access(mxr%f)) then
          write(x_unit(mxr%f),'(/,t4,"residual norms:")')
          write(x_unit(mxr%f),'(t6,"field = ",es8.2)') s_rn
          select case (mxr%atomic_status)
          case (ATOMIC_NON_TRIVIAL)
            write(x_unit(mxr%f),'(t6,"atomic = ",es8.2)') a_rn
          end select
        end if

      end subroutine

      subroutine reorder_i(mxr)
        type(mixer_rep) :: mxr

        integer :: is, ist, nmh, omh, os, ost
        real(double), dimension(:,:), pointer :: r1
        complex(double), dimension(:,:), pointer :: c1, c2

        nmh = mxr%mxh(mxr%step)
        omh = mxr%mxh(mxr%step-1)

        if (nmh < omh) then
          allocate( c1(size(mxr%fd_val,1),nmh), c2(size(mxr%fd_res,1),nmh) )
          os = mxr%hs - nmh
          if (os < 1) os = os + omh
          do is = 1,nmh
            c1(:,is) = mxr%fd_val(:,os)
            c2(:,is) = mxr%fd_res(:,os)
            os = os + 1
            if (os > omh) os = 1
          end do
          deallocate( mxr%fd_val, mxr%fd_res )
          mxr%fd_val => c1
          mxr%fd_res => c2
          nullify( c1, c2 )
          select case (mxr%atomic_status)
          case (ATOMIC_NON_TRIVIAL)
            allocate( c1(size(mxr%ad_val,1),nmh), c2(size(mxr%ad_res,1),nmh) )
            os = mxr%hs - nmh
            if (os < 1) os = os + omh
            do is = 1,nmh
              c1(:,is) = mxr%ad_val(:,os)
              c2(:,is) = mxr%ad_res(:,os)
              os = os + 1
              if (os > omh) os = 1
            end do
            deallocate( mxr%ad_val, mxr%ad_res )
            mxr%ad_val => c1
            mxr%ad_res => c2
            nullify( c1, c2 )
          end select
          allocate( r1(nmh,nmh) )
          os = mxr%hs - nmh
          if (os < 1) os = os + omh
          do is = 1,nmh
            ost = mxr%hs - nmh
            if (ost < 1) ost = ost + omh
            do ist = 1,nmh
              r1(is,ist) = mxr%mm(os,ost)
              ost = ost + 1
              if (ost > omh) ost = 1
            end do
            os = os + 1
            if (os > omh) os = 1
          end do
          deallocate( mxr%mm )
          mxr%mm => r1
          nullify( r1 )
          mxr%hs = 1
          mxr%nh = nmh
        elseif (nmh > omh) then
          allocate( c1(size(mxr%fd_val,1),nmh), c2(size(mxr%fd_res,1),nmh) )
          os = mxr%hs
          do is = 1,omh
            c1(:,is) = mxr%fd_val(:,os)
            c2(:,is) = mxr%fd_res(:,os)
            os = os + 1
            if (os > omh) os = 1
          end do
          deallocate( mxr%fd_val, mxr%fd_res )
          mxr%fd_val => c1
          mxr%fd_res => c2
          nullify( c1, c2 )
          select case (mxr%atomic_status)
          case (ATOMIC_NON_TRIVIAL)
            allocate( c1(size(mxr%ad_val,1),nmh), c2(size(mxr%ad_res,1),nmh) )
            os = mxr%hs
            do is = 1,omh
              c1(:,is) = mxr%ad_val(:,os)
              c2(:,is) = mxr%ad_res(:,os)
              os = os + 1
              if (os > omh) os = 1
            end do
            deallocate( mxr%ad_val, mxr%ad_res )
            mxr%ad_val => c1
            mxr%ad_res => c2
            nullify( c1, c2 )
          end select
          allocate( r1(nmh,nmh) )
          os = mxr%hs
          do is = 1,omh
            ost = mxr%hs
            do ist = 1,nmh
              r1(is,ist) = mxr%mm(os,ost)
              ost = ost + 1
              if (ost > omh) ost = 1
            end do
            os = os + 1
            if (os > omh) os = 1
          end do
          deallocate( mxr%mm )
          mxr%mm => r1
          nullify( r1 )
          mxr%hs = omh + 1
          mxr%nh = omh + 1
        end if

      end subroutine

    end module

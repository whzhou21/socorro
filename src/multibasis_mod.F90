!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module multibasis_mod
!doc$ module multibasis_mod

!     One datatype is available here: type(multibasis_obj).

!     Multibasis_mod provides basic capabilities for operations performed in multivector.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use math_mod
      use fft_mod
      use arg_mod
      use lattice_mod
      use ghost_mod
      use layout_mod
      use timing_mod

!cod$
      implicit none
      private

      ! usage
      integer, parameter :: NORMAL    = 1
      integer, parameter :: AUXILIARY = 2

      ! remap_type
      integer, parameter :: MPI_REMAP = 1
      integer, parameter :: SJP_REMAP = 2

      integer, parameter, public :: COLLECTIVE       = 1
      integer, parameter, public :: POINT_TO_POINT   = 2
      integer, parameter, public :: POINT_TO_POINT_2 = 3

      type, public :: mpi_remap_plan
        integer, dimension(:), pointer :: c1, d1, c2, d2
      end type

      type, public :: sjp_remap_plan
        integer :: ilo_in, ihi_in, jlo_in, jhi_in
        integer :: ilo_out, ihi_out, jlo_out, jhi_out
        real(double) :: f, r
      end type

      type, public :: multibasis_rep
        integer :: ref
        type(ghost) :: g
        integer :: usage                                                     ! usage wrt kgroup
        real(double) :: cutoff                                               ! cutoff energy
        real(double) :: cellvol                                              ! unit cell volume
        real(double), dimension(3) :: kpt                                    ! k-point in cartesian representation
        real(double), dimension(:,:), pointer :: gpt                         ! g values (distributed) in cartesian representation
        integer, dimension(:,:), pointer :: gridmap                          ! map into the non-zero grid elements
        integer :: first_g, last_g                                           ! plane-wave division information
        integer :: nbands                                                    ! number of bands
        integer, dimension(:), pointer :: first_band, last_band, my_band     ! band division information
        logical, dimension(:), pointer :: band_participant                   ! band group participation indicator
        integer :: exx_comm_method                                           ! communication method for exx calculations
        integer, dimension(:), pointer :: first_spair, last_spair, my_spair  ! collective small-pair division information
        integer, dimension(:,:), pointer :: spair_index                      ! collective small-pair index
        logical, dimension(:), pointer :: spair_participant                  ! collective small-pair group participant indicator
        integer, dimension(:), pointer :: first_lpair, last_lpair, my_lpair  ! collective large-pair division information
        integer, dimension(:,:), pointer :: lpair_index                      ! collective large-pair index
        logical, dimension(:), pointer :: lpair_participant                  ! collective large-pair group participant indicator
        integer, dimension(:,:), pointer :: r2r_dmap                         ! point-to-point row-to-row distribution map
        integer, dimension(:,:), pointer :: r2c_dmap                         ! point-to-point row-to-column distribution map
        integer, dimension(:,:), pointer :: c2c_smap                         ! point-to-point column-to-column shift map
        logical, dimension(:), pointer :: c2c_ulst                           ! point-to-point column-to-column update list
        integer, dimension(:,:), pointer :: c2r_cmap                         ! point-to-point column-to-row collection map
        integer, dimension(:,:), pointer :: r2r_cmap                         ! point-to-point row-to-row collection map
        integer, dimension(:,:), pointer :: rc_band_index                    ! point-to-point row/column band indices
        integer :: remap_type                                                ! type of remap
        type(mpi_remap_plan) :: mpi_plan1, mpi_plan2, mpi_plan3, mpi_plan4   ! mpi remap plans
        type(sjp_remap_plan) :: sjp_plan1, sjp_plan2, sjp_plan3, sjp_plan4   ! sjp remap plans
        type(fft_serial_plan_i) :: isplan1                                   ! plan for indexed serial FFT's
        type(fft_serial_plan_i) :: isplan2                                   ! plan for indexed serial FFT's
        type(layout_obj) :: lay                                              ! layout
      end type

      type, public :: multibasis_obj
        private
        integer :: ref
        type(multibasis_rep), pointer :: o
      end type

!doc$
      public :: multibasis
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_usage
      public :: x_n_bands
      public :: x_cutoff
      public :: x_n_gvectors
      public :: x_layout
      public :: band_remap
      public :: spair_remap
      public :: lpair_remap
      public :: form_vmap
      public :: wormhole
      public :: diary

!cod$
      interface multibasis
        module procedure constructor_mb
      end interface
      interface update
        module procedure update_mb
      end interface
      interface my
        module procedure my_mb, my_new_mb
      end interface
      interface thy
        module procedure thy_mb
      end interface
      interface glean
        module procedure glean_mb
      end interface
      interface bequeath
        module procedure bequeath_mb
      end interface
      interface assignment(=)
        module procedure assign_mb
      end interface
      interface x_ref
        module procedure mb_ref
      end interface
      interface x_ghost
        module procedure mb_ghost
      end interface
      interface x_usage
        module procedure mb_usage
      end interface
      interface x_n_bands
        module procedure mb_n_bands
      end interface
      interface x_cutoff
        module procedure mb_cutoff
      end interface
      interface x_n_gvectors
        module procedure mb_n_gvectors
      end interface
      interface x_layout
        module procedure mb_layout
      end interface
      interface band_remap
        module procedure band_remap_1d_to_2d_mb, band_remap_2d_to_1d_mb
      end interface
      interface spair_remap
        module procedure spair_remap_1d_to_2d_mb, spair_remap_2d_to_1d_mb
      end interface
      interface lpair_remap
        module procedure lpair_remap_1d_to_2d_mb, lpair_remap_2d_to_1d_mb
      end interface
      interface form_vmap
        module procedure form_vmap_mb
      end interface
      interface wormhole
        module procedure mb_wormhole
      end interface
      interface diary
        module procedure diary_mb
      end interface

      contains

! public routines

      function constructor_mb(usage,cutoff,nbands,kpt,lay) result(mb)
!doc$ function multibasis(usage,cutoff,nbands,kpt,lay) result(mb)
        character(line_len), intent(in) :: usage
        real(double), intent(in) :: cutoff
        integer, intent(in) :: nbands
        real(double), dimension(3), intent(in) :: kpt
        type(layout_obj) :: lay
        type(multibasis_obj) :: mb
!       requires: usage be "normal" or "auxiliary"; kpt be in reciprocal-lattice representation.
!       effects: Returns a new multibasis object.
!       errors: None.
        
!cod$
        logical :: found
        logical, dimension(:,:,:), allocatable :: pattern
        character(line_len) :: tag
        integer :: i1, i2, i3, ic, ig, ip, is, isr, m, myp, my_ng, nb, nc, ng, ngroups, np, npairs, ns, size_g
        integer :: level, step_down, step_up
        integer :: block, ctag, length, sp, rp, row_band, col_band, tmp_band
        integer, dimension(3) :: nf
        integer, dimension(4) :: mark
        integer, dimension(:,:), allocatable :: fft_index
        real(double), dimension(:,:,:), pointer :: cp, fx, fy, fz

        call my(lay)

! INITIALIZE THE MULTIBASIS.

        mb%ref = 0
        allocate( mb%o )
        mb%o%ref = 0
        mb%o%g = x_ghost()

        select case (trim(usage))
        case ("normal")
          mb%o%usage = NORMAL
        case ("auxiliary")
          mb%o%usage = AUXILIARY
        end select

        mb%o%cutoff = cutoff
        mb%o%nbands = nbands
        call my(lay,mb%o%lay)
        mb%o%cellvol = x_cell_volume(x_lattice(mb%o%lay))
        mb%o%kpt = lat2f(x_lattice(mb%o%lay),kpt)

        nullify( cp, fx, fy, fz )

! GENERATE PACKING INFORMATION

        ! Generate full mesh
        call fmesh(fx,fy,fz,mb%o%lay,S_TYPE)
        call alloc(cp,mb%o%lay,S_TYPE)
        do i3 = 1,size(cp,3)
        do i2 = 1,size(cp,2)
        do i1 = 1,size(cp,1)
          cp(i1,i2,i3) = (fx(i1,i2,i3) + mb%o%kpt(1))**2 + (fy(i1,i2,i3) + mb%o%kpt(2))**2 + (fz(i1,i2,i3) + mb%o%kpt(3))**2
        end do
        end do
        end do
        ! Count n meshpoints within the cutoff
        ng = 0
        do i3 = 1,size(cp,3)
        do i2 = 1,size(cp,2)
        do i1 = 1,size(cp,1)
          if (cp(i1,i2,i3) < mb%o%cutoff) ng = ng + 1
        end do
        end do
        end do

        ! Subdivide the cutoff meshpoints
        call subdivide(mpi_myproc(KGROUP),mpi_nprocs(KGROUP),1,ng,mb%o%first_g,mb%o%last_g,size_g)
        allocate( mb%o%gridmap(3,ng), mb%o%gpt(size_g,3) )
        ! Form the per-processor gpt and gridmap lists
        ng = 0
        do i3 = 1,size(cp,3)
        do i2 = 1,size(cp,2)
        do i1 = 1,size(cp,1)
          if (cp(i1,i2,i3) < mb%o%cutoff) then
            ng = ng + 1
            mb%o%gridmap(1,ng) = i1
            mb%o%gridmap(2,ng) = i2
            mb%o%gridmap(3,ng) = i3
            if ((ng >= mb%o%first_g) .and. (ng <= mb%o%last_g)) then
              my_ng = ng - mb%o%first_g + 1
              mb%o%gpt(my_ng,1) = fx(i1,i2,i3)
              mb%o%gpt(my_ng,2) = fy(i1,i2,i3)
              mb%o%gpt(my_ng,3) = fz(i1,i2,i3)
            end if
          end if
        end do
        end do
        end do

! GENERATE BAND GROUPS FOR REDISTRIBUTING THE WAVEFUNCTIONS.

        np = mpi_nprocs(KGROUP)
        nb = mb%o%nbands

        ! Band groups
        ngroups = (nb + np - 1)/np
        allocate( mb%o%first_band(ngroups) )
        allocate( mb%o%last_band(ngroups) )
        allocate( mb%o%my_band(ngroups) )
        allocate( mb%o%band_participant(ngroups) )
        do ig = 1,ngroups
          mb%o%first_band(ig) = np*(ig-1) + 1
          mb%o%last_band(ig) = min(nb,np*ig)
          if ((mb%o%first_band(ig) + mpi_myproc(KGROUP)) <= nb) then
            mb%o%my_band(ig) = mb%o%first_band(ig) + mpi_myproc(KGROUP)
            mb%o%band_participant(ig) = .true.
          else
            mb%o%my_band(ig) = 0
            mb%o%band_participant(ig) = .false.
          end if
        end do

! GENERATE COMMUNICATION MAPS FOR CALCULATING EXCHANGE QUANTITIES.

        nullify( mb%o%first_spair )
        nullify( mb%o%last_spair )
        nullify( mb%o%my_spair )
        nullify( mb%o%spair_index )
        nullify( mb%o%spair_participant )
        nullify( mb%o%first_lpair )
        nullify( mb%o%last_lpair )
        nullify( mb%o%my_lpair )
        nullify( mb%o%lpair_index )
        nullify( mb%o%lpair_participant )

        nullify( mb%o%r2r_dmap )
        nullify( mb%o%r2c_dmap )
        nullify( mb%o%c2c_smap )
        nullify( mb%o%c2c_ulst )
        nullify( mb%o%c2r_cmap )
        nullify( mb%o%r2r_cmap )
        nullify( mb%o%rc_band_index )

        np = mpi_nprocs(KGROUP)
        nb = mb%o%nbands

        ! Communication method used in exx energy and derivative calculations
        call arglc("exx_comm_method",tag,found)
        if (.not.found) tag = "collective"
        select case (trim(tag))
        case ("collective")

          mb%o%exx_comm_method = COLLECTIVE

          ! Small-pair groups
          npairs = (nb*(nb + 1))/2
          allocate( mb%o%spair_index(2,npairs) )
          ip = 0
          do i1 = 1,nb
            do i2 = i1,nb
              ip = ip + 1
              mb%o%spair_index(1,ip) = i1
              mb%o%spair_index(2,ip) = i2
            end do
          end do
          ngroups = (npairs + np - 1)/np
          allocate( mb%o%first_spair(ngroups), mb%o%last_spair(ngroups), mb%o%my_spair(ngroups) )
          allocate( mb%o%spair_participant(ngroups) )
          do ig = 1,ngroups
            mb%o%first_spair(ig) = np*(ig-1) + 1
            mb%o%last_spair(ig) = min(npairs,np*ig)
            if ((mb%o%first_spair(ig) + mpi_myproc(KGROUP)) <= npairs) then
              mb%o%my_spair(ig) = mb%o%first_spair(ig) + mpi_myproc(KGROUP)
              mb%o%spair_participant(ig) = .true.
            else
              mb%o%my_spair(ig) = 0
              mb%o%spair_participant(ig) = .false.
            end if
          end do

          ! Large-pair groups
          npairs = nb*nb
          allocate( mb%o%lpair_index(2,npairs) )
          ip = 0
          do i1 = 1,nb
            do i2 = 1,nb
              ip = ip + 1
              mb%o%lpair_index(1,ip) = i1
              mb%o%lpair_index(2,ip) = i2
            end do
          end do
          ngroups = (npairs + np - 1)/np
          allocate( mb%o%first_lpair(ngroups), mb%o%last_lpair(ngroups), mb%o%my_lpair(ngroups) )
          allocate( mb%o%lpair_participant(ngroups) )
          do ig = 1,ngroups
            mb%o%first_lpair(ig) = np*(ig-1) + 1
            mb%o%last_lpair(ig) = min(npairs,np*ig)
            if ((mb%o%first_lpair(ig) + mpi_myproc(KGROUP)) <= npairs) then
              mb%o%my_lpair(ig) = mb%o%first_lpair(ig) + mpi_myproc(KGROUP)
              mb%o%lpair_participant(ig) = .true.
            else
              mb%o%my_lpair(ig) = 0
              mb%o%lpair_participant(ig) = .false.
            end if
          end do

        case ("point-to-point")

          mb%o%exx_comm_method = POINT_TO_POINT

          myp = mpi_myproc(KGROUP)

          ! Check constraints on np vs. nb
          if (error(np < nb,"ERROR: np < nb")) goto 100
          if (error(np > nb*(nb + 1)/2,"ERROR: np > nb(nb+1)/2")) goto 100
          if (error(mod(np,nb) /= 0,"ERROR: np is not an integer multiple of nb")) goto 100

          m = np/nb

          ! row-to-row distribution maps
          ns = 0
          do                                              ! Calculate the number of stages in the mapping.
            length = (2**ns)*nb
            if (length >= np) exit
            ns = ns + 1
          end do
          if (ns > 0) then                                ! Allocate and initialize the mapping array.
            allocate( mb%o%r2r_dmap(ns,2) )
            mb%o%r2r_dmap = -1
          end if
          if (associated( mb%o%r2r_dmap )) then           ! Calculate and store the mappings.
            do is = 1,size(mb%o%r2r_dmap,1)
              length = 2**(is - 1)*nb                     ! The algorithm works with two blocks (0 and 1) of
              block = int(myp/length)                     !   length "length", sending and receiving between them.
              if (block == 0) then                        ! Send to a process in block 1:
                sp = modulo(myp,length) + length          !   Process that myp sends to (if within range).
                if (sp > (np - 1)) cycle
                mb%o%r2r_dmap(is,1) = sp
              elseif (block == 1) then                    ! Receive from a process in block 0:
                mb%o%r2r_dmap(is,2) = modulo(myp,length)  !   Process that myp receives from.
              end if
            end do
          end if

          ! row-to-column distribution maps
          if (even(nb) .and. np == nb*(nb + 1)/2) then              ! Calculate the number of stages in the mapping.
            ns = 2                                                  ! For this case, column-to-column shifts are not needed.
          else
            ns = 1                                                  ! For this case, column-to-column shifts may be needed
          end if
          allocate( mb%o%r2c_dmap(ns,2) )                           ! Allocate and initialize the mapping array.
          mb%o%r2c_dmap = -1
          block = int(myp/nb)                                       ! Calculate and store the mappings.
          if (block == 0)  then                                     ! Processes in the first nb block.
            mb%o%r2c_dmap(1,1) = myp                                  ! Process that myp sends to in an nb block.
            mb%o%r2c_dmap(1,2) = myp                                  ! Process that myp receives from in an nb block.
            if (even(nb)) then
              if (np == nb*(nb + 1)/2 .and. (myp + 1) > nb/2 ) then
                mb%o%r2c_dmap(2,1) = myp + nb*(nb - 1)/2              ! Process that myp sends to in the half-nb block.
              end if
            end if
          elseif (block < int(np/nb)) then                          ! Processes in the remaining nb blocks.
            mb%o%r2c_dmap(1,1) = modulo(myp - block,nb) + block*nb    ! Process that myp sends to in a nb block.
            mb%o%r2c_dmap(1,2) = modulo(myp + block,nb) + block*nb    ! Process that myp receives from in an nb block.
          else                                                      ! Processes in the half-nb block (if there is one).
            mb%o%r2c_dmap(2,2) = mod(myp,nb) + nb/2                   ! Process that myp receives from in an nb block.
          end if

          ! column-to-column shift maps
          ns = 0                                                    ! Calculate the number of stages in the mapping.
          do                                                          ! Groups of length np (when np < nb(nb + 1)/2).
            if (np*(ns + 1) >= nb*(nb + 1)/2 ) exit
            ns = ns + 1
          end do
          if (ns > 0) then                                          ! Allocate and initialize the mapping array.
            allocate( mb%o%c2c_smap(ns,2) )
            mb%o%c2c_smap = -1
          end if
          if (associated( mb%o%c2c_smap )) then                     ! Calculate and store the mappings.
            do is = 1,size(mb%o%c2c_smap,1)
              if (odd(nb)) then
                if (myp + is*np < nb*(nb + 1)/2) then                        ! Processes within range of nb-length blocks.
                  mb%o%c2c_smap(is,1) = modulo(myp - m,nb) + int(myp/nb)*nb    ! Process that myp sends to.
                  mb%o%c2c_smap(is,2) = modulo(myp + m,nb) + int(myp/nb)*nb    ! Process that myp receives from.
                end if
              elseif (even(nb)) then
                if (myp + is*np < nb*(nb + 2)/2) then                        ! Processes within range of nb-length blocks.
                  mb%o%c2c_smap(is,1) = modulo(myp - m,nb) + int(myp/nb)*nb    ! Process that myp sends to.
                  mb%o%c2c_smap(is,2) = modulo(myp + m,nb) + int(myp/nb)*nb    ! Process that myp receives from.
                end if
              end if
            end do
          end if

          ! column-to-column update lists
          if (associated( mb%o%c2c_smap )) then
            ns = size(mb%o%c2c_smap,1)                              ! Allocate and initialize the update array.
            allocate( mb%o%c2c_ulst(ns) )
            mb%o%c2c_ulst = .false.
            do is = 1,size(mb%o%c2c_ulst,1)                         ! Calculate and store the updates.
              if (myp + is*np < nb*(nb + 1)/2) then                 ! Processes within range of row-column pairs.
                mb%o%c2c_ulst(is) = .true.
              end if
            end do
          end if

          ! column-to-row collection maps
          ns = size(mb%o%r2c_dmap,1)                                ! Allocate and initialize the mapping array.
          allocate( mb%o%c2r_cmap(ns,2) )
          mb%o%c2r_cmap = -1
          if (associated( mb%o%c2c_smap )) then                     ! For the case in which there are c2c maps.
            sp = myp                                                ! Process that myp sends to.
            do is = size(mb%o%c2c_smap,1),1,-1                        ! c2c_smap contribution.
              if (odd(nb)) then
                if (sp + is*np < nb*(nb + 1)/2) then
                  sp = modulo(sp + m,nb) + int(sp/nb)*nb
                end if
              elseif (even(nb)) then
                if (sp + is*np < nb*(nb + 2)/2) then
                  sp = modulo(sp + m,nb) + int(sp/nb)*nb
                end if
              end if
            end do
            block = int(sp/nb)                                        ! r2c_dmap contribution.
            if (block == 0) then
              mb%o%c2r_cmap(1,1) = sp
            elseif (block < int(np/nb)) then
              mb%o%c2r_cmap(1,1) = modulo(sp + block,nb) + block*nb
            end if
            rp = mb%o%r2c_dmap(1,1)                                 ! Process that myp receives from; r2c_dmap contribution.
            do is = 1,size(mb%o%c2c_smap,1)                           ! c2c_smap contribution.
              if (odd(nb)) then
                if (rp + is*np < nb*(nb + 1)/2) then
                  rp = modulo(rp - m,nb) + int(rp/nb)*nb
                end if
              elseif (even(nb)) then
                if (rp + is*np < nb*(nb + 2)/2) then
                  rp = modulo(rp - m,nb) + int(rp/nb)*nb
                end if
              end if
            end do
            mb%o%c2r_cmap(1,2) = rp
          else                                                      ! For the case in which there are no c2c maps.
            do is = 1,size(mb%o%c2r_cmap,1)
              mb%o%c2r_cmap(is,1) = mb%o%r2c_dmap(is,2)
              mb%o%c2r_cmap(is,2) = mb%o%r2c_dmap(is,1)
            end do
          end if

          ! row-to-row collection maps
          if (associated( mb%o%r2r_dmap )) then                     ! Allocate the mapping array.
            ns = size(mb%o%r2r_dmap,1)
            allocate( mb%o%r2r_cmap(ns,2) )
          end if
          if (associated( mb%o%r2r_cmap )) then                     ! Calculate and store the mappings.
            isr = ns
            do is = 1,ns
              mb%o%r2r_cmap(is,1) = mb%o%r2r_dmap(isr,2)              ! Process that myp sends to.
              mb%o%r2r_cmap(is,2) = mb%o%r2r_dmap(isr,1)              ! Process that myp receives from.
              isr = isr - 1
            end do
          end if

          ! initial row and column band indexes (after band remap)
          row_band = mb%o%my_band(1)

          if (associated( mb%o%r2r_dmap )) then   ! distribute row indexes to the rest of the processes
            do is = 1,size(mb%o%r2r_dmap,1)
              ctag = is
              sp = mb%o%r2r_dmap(is,1)
              if (sp /= -1) then
                call nonblocking_send(KGROUP,row_band,sp,ctag)  ; if (error()) goto 10
              end if
              rp = mb%o%r2r_dmap(is,2)
              if (rp /= -1) then
                call blocking_recv(KGROUP,row_band,rp,ctag)     ; if (error()) goto 10
              end if
10            call barrier(KGROUP)
              call sync_kgroup_process_errors() ; if (error()) goto 100
            end do
          end if

          do is = 1,size(mb%o%r2c_dmap,1)        ! distribute row indexes to the column indexes
            ctag = is
            sp = mb%o%r2c_dmap(is,1)
            rp = mb%o%r2c_dmap(is,2)
            if (sp /= -1 .and. rp /= -1) then
              if (sp == rp) then
                col_band = row_band
              else
                call nonblocking_send(KGROUP,row_band,sp,ctag)   ; if (error()) goto 20
                call blocking_recv(KGROUP,col_band,rp,ctag)      ; if (error()) goto 20
              end if
            elseif (sp /= -1 .or. rp /= -1) then
              if (sp /= -1) then
                call nonblocking_send(KGROUP,row_band,sp,ctag)   ; if (error()) goto 20
              end if
              if (rp /= -1) then
                call blocking_recv(KGROUP,col_band,rp,ctag)      ; if (error()) goto 20
              end if
            end if
20          call barrier(KGROUP)
            call sync_kgroup_process_errors() ; if (error()) goto 100
          end do

          nc = 1                                  ! initialize rc_band_index
          if (associated( mb%o%c2c_smap )) then
            nc = nc + size(mb%o%c2c_smap,1)
          end if
          allocate( mb%o%rc_band_index(nc,2) )
          ic = 1
          mb%o%rc_band_index(ic,1) = row_band
          mb%o%rc_band_index(ic,2) = col_band

          if (associated( mb%o%c2c_smap )) then  ! shift the column indexes
            do is = 1,size(mb%o%c2c_smap,1)
              ctag = is
              sp = mb%o%c2c_smap(is,1)
              rp = mb%o%c2c_smap(is,2)
              if (sp /= -1 .and. rp /= -1) then
                tmp_band = col_band
                call nonblocking_send(KGROUP,tmp_band,sp,ctag)   ; if (error()) goto 30
                call blocking_recv(KGROUP,col_band,rp,ctag)      ; if (error()) goto 30
              end if
              ic = ic + 1
              mb%o%rc_band_index(ic,1) = row_band
              mb%o%rc_band_index(ic,2) = col_band
30            call barrier(KGROUP)
              call sync_kgroup_process_errors() ; if (error()) goto 100
            end do
          end if

        case ("point-to-point_2")

          mb%o%exx_comm_method = POINT_TO_POINT_2

          myp = mpi_myproc(KGROUP)

          ! Check constraints on np.
          if (error(np < nb,"ERROR: np < nb")) goto 100
          if (error(np > nb*nb,"ERROR: np > nb*nb")) goto 100
          if (error(mod(np,nb) /= 0,"ERROR: np is not an integer multiple of nb")) goto 100

          ! multiplier
          m = np/nb

          ! row-to-row distribution maps
          ns = 0                                          ! Calculate the number of stages in the mappings.
          do
            length = (2**ns)*nb
            if (length >= np) exit
            ns = ns + 1
          end do
          if (ns > 0) then                                ! Allocate and initialize the mapping array.
            allocate( mb%o%r2r_dmap(ns,2) )
            mb%o%r2r_dmap = -1
          end if
          if (associated( mb%o%r2r_dmap )) then           ! Calculate and store the mappings.
            do is = 1,size(mb%o%r2r_dmap,1)
              length = 2**(is - 1)*nb                       ! The algorithm uses two blocks of processes
              block = int(myp/length)                       ! (0 and 1), sending and receiving between them.
              if (block == 0) then                            ! Process myp sends to in block 1.
                sp = modulo(myp,length) + length
                if (sp > (np - 1)) cycle
                mb%o%r2r_dmap(is,1) = sp
              elseif (block == 1) then                        ! Process myp receives from in block 0.
                rp = modulo(myp,length)
                mb%o%r2r_dmap(is,2) = rp
              end if
            end do
          end if

          ! row-to-column distribution maps
          ns = 1                                            ! The mappings are within nb blocks.
          allocate( mb%o%r2c_dmap(ns,2) )                   ! Allocate the mapping array.
          block = int(myp/nb)                                 ! nb-block number of process myp.
          sp = modulo(myp - block,nb) + block*nb              ! Process myp sends to in its nb block.
          rp = modulo(myp + block,nb) + block*nb              ! Process myp receives from in its nb block.
          mb%o%r2c_dmap(1,1) = sp
          mb%o%r2c_dmap(1,2) = rp

          ! column-to-column shift maps
          ns = 0                                            ! Calculate the number of shifts needed.
          do
            if (np*(ns + 1) >= nb*nb) exit
            ns = ns + 1
          end do
          if (ns > 0) then                                  ! Allocate and initialize the mapping array.
            allocate( mb%o%c2c_smap(ns,2) )
            mb%o%c2c_smap = -1
          end if
          if (associated( mb%o%c2c_smap )) then             ! Calculate and store the mappings.
            do is = 1,size(mb%o%c2c_smap,1)
              if (myp + is*np < nb*nb) then
                sp = modulo(myp - m,nb) + int(myp/nb)*nb      ! Process myp sends.
                rp = modulo(myp + m,nb) + int(myp/nb)*nb      ! Process myp receives from.
                mb%o%c2c_smap(is,1) = sp
                mb%o%c2c_smap(is,2) = rp
              end if
            end do
          end if

          ! column-to-column update lists
          if (associated( mb%o%c2c_smap )) then
            ns = size(mb%o%c2c_smap,1)                      ! Allocate and initialize the update array.
            allocate( mb%o%c2c_ulst(ns) )
            mb%o%c2c_ulst = .false.
            do is = 1,size(mb%o%c2c_ulst,1)                 ! Calculate and store the update directives.
              if (myp + is*np < nb*nb) then
                mb%o%c2c_ulst(is) = .true.
              end if
            end do
          end if

          ! row-to-row collection maps
          if (associated( mb%o%r2r_dmap )) then             ! Allocate the mapping array.
            ns = size(mb%o%r2r_dmap,1)
            allocate( mb%o%r2r_cmap(ns,2) )
          end if
          if (associated( mb%o%r2r_cmap )) then             ! Calculate and store the mappings.
            isr = ns
            do is = 1,ns
              sp = mb%o%r2r_dmap(isr,2)                       ! Process myp sends to.
              rp = mb%o%r2r_dmap(isr,1)                       ! Process myp receives from.
              mb%o%r2r_cmap(is,1) = sp
              mb%o%r2r_cmap(is,2) = rp
              isr = isr - 1
            end do
          end if

          ! row and column band indices
          row_band = mb%o%my_band(1)                        ! Initial row indices (after the band remap).

          if (associated( mb%o%r2r_dmap )) then             ! Distribute row indexes to other nb blocks.
            do is = 1,size(mb%o%r2r_dmap,1)
              ctag = is
              sp = mb%o%r2r_dmap(is,1)
              rp = mb%o%r2r_dmap(is,2)
              if (sp /= -1) then
                call nonblocking_send(KGROUP,row_band,sp,ctag)  ; if (error()) goto 40
              end if
              if (rp /= -1) then
                call blocking_recv(KGROUP,row_band,rp,ctag)     ; if (error()) goto 40
              end if
40            call barrier(KGROUP)
              call sync_kgroup_process_errors() ; if (error()) goto 100
            end do
          end if

          ctag = 1                                          ! Distribute column indexes in each nb block.
          sp = mb%o%r2c_dmap(1,1)
          rp = mb%o%r2c_dmap(1,2)
          if (sp == myp .and. rp == myp) then
            col_band = row_band
          else
            call nonblocking_send(KGROUP,row_band,sp,ctag)   ; if (error()) goto 50
            call blocking_recv(KGROUP,col_band,rp,ctag)      ; if (error()) goto 50
          end if
50        call barrier(KGROUP)
          call sync_kgroup_process_errors() ; if (error()) goto 100

          nc = 1                                            ! Allocate and initialize the band index array.
          if (associated( mb%o%c2c_smap )) then
            nc = nc + size(mb%o%c2c_smap,1)
          end if
          allocate( mb%o%rc_band_index(nc,2) )
          ic = 1
          mb%o%rc_band_index(ic,1) = row_band
          mb%o%rc_band_index(ic,2) = col_band

          if (associated( mb%o%c2c_smap )) then             ! Calculate and store the shifted column band indices.
            do is = 1,size(mb%o%c2c_smap,1)
              ctag = is
              sp = mb%o%c2c_smap(is,1)
              rp = mb%o%c2c_smap(is,2)
              if (sp /= -1 .and. rp /= -1) then
                tmp_band = col_band
                call nonblocking_send(KGROUP,tmp_band,sp,ctag)   ; if (error()) goto 60
                call blocking_recv(KGROUP,col_band,rp,ctag)      ; if (error()) goto 60
              end if
              ic = ic + 1
              mb%o%rc_band_index(ic,1) = row_band
              mb%o%rc_band_index(ic,2) = col_band
60            call barrier(KGROUP)
              call sync_kgroup_process_errors() ; if (error()) goto 100
            end do
          end if

        case default

          if (error(.true.,"ERROR: exx_comm_method was not recognized")) goto 100

        end select

! CREATE REMAP PLANS

        nullify( mb%o%mpi_plan1%c1 )
        nullify( mb%o%mpi_plan1%d1 )
        nullify( mb%o%mpi_plan1%c2 )
        nullify( mb%o%mpi_plan1%d2 )

        nullify( mb%o%mpi_plan2%c1 )
        nullify( mb%o%mpi_plan2%d1 )
        nullify( mb%o%mpi_plan2%c2 )
        nullify( mb%o%mpi_plan2%d2 )

        nullify( mb%o%mpi_plan3%c1 )
        nullify( mb%o%mpi_plan3%d1 )
        nullify( mb%o%mpi_plan3%c2 )
        nullify( mb%o%mpi_plan3%d2 )

        nullify( mb%o%mpi_plan4%c1 )
        nullify( mb%o%mpi_plan4%d1 )
        nullify( mb%o%mpi_plan4%c2 )
        nullify( mb%o%mpi_plan4%d2 )

        call arglc("remap_type",tag,found)
        if (.not.found) tag = "mpi"
        select case (trim(tag))
        case ("mpi")

          mb%o%remap_type = MPI_REMAP

          ! plan 1: band, small-pair, and large-pair groups where all processors are participants
          allocate( mb%o%mpi_plan1%c1(0:mpi_nprocs(KGROUP)-1) )
          allocate( mb%o%mpi_plan1%d1(0:mpi_nprocs(KGROUP)-1) )
          allocate( mb%o%mpi_plan1%c2(0:mpi_nprocs(KGROUP)-1) )
          allocate( mb%o%mpi_plan1%d2(0:mpi_nprocs(KGROUP)-1) )
          ng = size(mb%o%gridmap,2)
          do i1 = 0,mpi_nprocs(KGROUP)-1
            call subdivide(i1,mpi_nprocs(KGROUP),1,ng,i3,ig,mb%o%mpi_plan1%c1(i1))
            mb%o%mpi_plan1%d1(i1) = 0
            do i2 = 0,i1-1
              mb%o%mpi_plan1%d1(i1) = mb%o%mpi_plan1%d1(i1) + mb%o%mpi_plan1%c1(i2)
            end do
          end do
          my_ng = size(mb%o%gpt,1)
          do i1 = 0,mpi_nprocs(KGROUP)-1
            mb%o%mpi_plan1%c2(i1) = my_ng
            mb%o%mpi_plan1%d2(i1) = 0
            do i2 = 0,i1-1
              mb%o%mpi_plan1%d2(i1) = mb%o%mpi_plan1%d2(i1) + mb%o%mpi_plan1%c2(i2)
            end do
          end do

          ! plan 2: band group where some processors may not be participants
          allocate( mb%o%mpi_plan2%c1(0:mpi_nprocs(KGROUP)-1) )
          allocate( mb%o%mpi_plan2%d1(0:mpi_nprocs(KGROUP)-1) )
          allocate( mb%o%mpi_plan2%c2(0:mpi_nprocs(KGROUP)-1) )
          allocate( mb%o%mpi_plan2%d2(0:mpi_nprocs(KGROUP)-1) )
          mb%o%mpi_plan2%c1 = mb%o%mpi_plan1%c1
          mb%o%mpi_plan2%d1 = mb%o%mpi_plan1%d1
          ngroups = size(mb%o%my_band)
          if (.not.mb%o%band_participant(ngroups)) mb%o%mpi_plan2%c1 = 0
          mb%o%mpi_plan2%c2 = mb%o%mpi_plan1%c2
          mb%o%mpi_plan2%d2 = mb%o%mpi_plan1%d2
          do i1 = 0,mpi_nprocs(KGROUP)-1
            if ((mb%o%first_band(ngroups) + i1) > mb%o%last_band(ngroups)) mb%o%mpi_plan2%c2(i1) = 0
          end do

          select case (mb%o%exx_comm_method)
          case (COLLECTIVE)

            ! plan 3: small-pair group where some processors may not be participants
            allocate( mb%o%mpi_plan3%c1(0:mpi_nprocs(KGROUP)-1) )
            allocate( mb%o%mpi_plan3%d1(0:mpi_nprocs(KGROUP)-1) )
            allocate( mb%o%mpi_plan3%c2(0:mpi_nprocs(KGROUP)-1) )
            allocate( mb%o%mpi_plan3%d2(0:mpi_nprocs(KGROUP)-1) )
            mb%o%mpi_plan3%c1 = mb%o%mpi_plan1%c1
            mb%o%mpi_plan3%d1 = mb%o%mpi_plan1%d1
            ngroups = size(mb%o%my_spair)
            if (.not.mb%o%spair_participant(ngroups)) mb%o%mpi_plan3%c1 = 0
            mb%o%mpi_plan3%c2 = mb%o%mpi_plan1%c2
            mb%o%mpi_plan3%d2 = mb%o%mpi_plan1%d2
            do i1 = 0,mpi_nprocs(KGROUP)-1
              if ((mb%o%first_spair(ngroups) + i1) > mb%o%last_spair(ngroups)) mb%o%mpi_plan3%c2(i1) = 0
            end do

            ! plan 4: large-pair group where some processors may not be participants
            allocate( mb%o%mpi_plan4%c1(0:mpi_nprocs(KGROUP)-1) )
            allocate( mb%o%mpi_plan4%d1(0:mpi_nprocs(KGROUP)-1) )
            allocate( mb%o%mpi_plan4%c2(0:mpi_nprocs(KGROUP)-1) )
            allocate( mb%o%mpi_plan4%d2(0:mpi_nprocs(KGROUP)-1) )
            mb%o%mpi_plan4%c1 = mb%o%mpi_plan1%c1
            mb%o%mpi_plan4%d1 = mb%o%mpi_plan1%d1
            ngroups = size(mb%o%my_lpair)
            if (.not.mb%o%lpair_participant(ngroups)) mb%o%mpi_plan4%c1 = 0
            mb%o%mpi_plan4%c2 = mb%o%mpi_plan1%c2
            mb%o%mpi_plan4%d2 = mb%o%mpi_plan1%d2
            do i1 = 0,mpi_nprocs(KGROUP)-1
              if ((mb%o%first_lpair(ngroups) + i1) > mb%o%last_lpair(ngroups)) mb%o%mpi_plan4%c2(i1) = 0
            end do

          end select

        case ("sjp")

          mb%o%remap_type = SJP_REMAP

          ! plan 1: band, small-pair, and large-pair groups where all processors are participants
          mb%o%sjp_plan1%ilo_in = mb%o%first_g
          mb%o%sjp_plan1%ihi_in = mb%o%last_g
          mb%o%sjp_plan1%jlo_in = 1
          mb%o%sjp_plan1%jhi_in = mpi_nprocs(KGROUP)
          mb%o%sjp_plan1%ilo_out = 1
          mb%o%sjp_plan1%ihi_out = size(mb%o%gridmap,2)
          mb%o%sjp_plan1%jlo_out = mpi_myproc(KGROUP) + 1
          mb%o%sjp_plan1%jhi_out = mpi_myproc(KGROUP) + 1
          call create_sjp_remap_plan_i(mb%o%sjp_plan1)

          ! plan 2: band group where some processors may not be participants
          mb%o%sjp_plan2 = mb%o%sjp_plan1
          ngroups = size(mb%o%my_band)
          mb%o%sjp_plan2%jhi_in = mb%o%last_band(ngroups) - mb%o%first_band(ngroups) + 1
          if (any(.not.mb%o%band_participant)) mb%o%sjp_plan2%jhi_out = mb%o%sjp_plan2%jlo_out - 1
          call create_sjp_remap_plan_i(mb%o%sjp_plan2)

          select case (mb%o%exx_comm_method)
          case (COLLECTIVE)

            ! plan 3: small-pair group where some processors may not be participants
            mb%o%sjp_plan3 = mb%o%sjp_plan1
            ngroups = size(mb%o%my_spair)
            mb%o%sjp_plan3%jhi_in = mb%o%last_spair(ngroups) - mb%o%first_spair(ngroups) + 1
            if (any(.not.mb%o%spair_participant)) mb%o%sjp_plan3%jhi_out = mb%o%sjp_plan3%jlo_out - 1
            call create_sjp_remap_plan_i(mb%o%sjp_plan3)

            ! plan 4: large-pair group where some processors may not be participants
            mb%o%sjp_plan4 = mb%o%sjp_plan1
            ngroups = size(mb%o%my_lpair)
            mb%o%sjp_plan4%jhi_in = mb%o%last_lpair(ngroups) - mb%o%first_lpair(ngroups) + 1
            if (any(.not.mb%o%lpair_participant)) mb%o%sjp_plan4%jhi_out = mb%o%sjp_plan4%jlo_out - 1
            call create_sjp_remap_plan_i(mb%o%sjp_plan4)

          end select

        case default

          if (error(.true.,"ERROR: remap_type was not recognized")) goto 100

        end select

! GENERATE PLANS FOR PRUNED S_TYPE FFT'S.

        nf = x_dims(mb%o%lay)
        allocate( pattern(nf(1),nf(2),nf(3)) )
        pattern = .false.
        do ig = 1,size(mb%o%gridmap,2)
          pattern(mb%o%gridmap(1,ig),mb%o%gridmap(2,ig),mb%o%gridmap(3,ig)) = .true.
        end do

        allocate( fft_index(5,nf(3)) )
        do i3 = 1,nf(3)

          level = 0
          step_up   = 0
          step_down = 0
          mark = 0

          do i2 = 1,nf(2)-1

            select case (level)
            case (0)
              if ( any(pattern(:,i2,i3)) ) then
                level = 1
                step_up = step_up + 1
                mark(step_up + step_down) = i2
              end if
            case (1)
              if ( .not.any(pattern(:,i2,i3)) ) then
                level = 0
                step_down = step_down + 1
                mark(step_up + step_down) = i2 - 1
              end if
            end select

          end do

          select case (level)
          case (0)
            if ( any(pattern(:,nf(2),i3)) ) then
              level = 1
              step_up = step_up + 1
              mark(step_up + step_down) = nf(2)
              step_down = step_down + 1
              mark(step_up + step_down) = nf(2)
            end if
          case (1)
            if ( .not.any(pattern(:,nf(2),i3)) ) then
              level = 0
              step_down = step_down + 1
              mark(step_up + step_down) = nf(2) - 1
            else
              step_down = step_down + 1
              mark(step_up + step_down) = nf(2)
            end if
          end select

          fft_index(1,i3) = step_up
          fft_index(2:5,i3) = 0

          select case (step_up)
          case (1)
            fft_index(2,i3) = mark(1)
            fft_index(3,i3) = mark(2) - mark(1) + 1
          case (2)
            fft_index(2,i3) = mark(1)
            fft_index(3,i3) = mark(2) - mark(1) + 1
            fft_index(4,i3) = mark(3)
            fft_index(5,i3) = mark(4) - mark(3) + 1
          end select

        end do

        ! generate two plans for usage = NORMAL, one plan for usage = AUXILIARY
        call fft_create_serial_plan(nf,fft_index,mb%o%isplan1)
        select case (mb%o%usage)
        case (NORMAL)
          call fft_create_serial_plan(nf,fft_index,mb%o%isplan2)
        end select

100     if (associated( cp )) deallocate( cp )
        if (associated( fx )) deallocate( fx )
        if (associated( fy )) deallocate( fy )
        if (associated( fz )) deallocate( fz )
        if (allocated( pattern )) deallocate( pattern )
        if (allocated( fft_index )) deallocate( fft_index )

        call glean(thy(lay))

        if (error("Exit multibasis_mod::constructor_mb")) continue

      end function

      subroutine update_mb(mb,lay)
!doc$ subroutine update(mb,lay)
        type(multibasis_obj) :: mb
        type(layout_obj) :: lay
!       effects: Updates mb.

!cod$
        call my(mb)
        call my(lay)

        if (x_ghost(lay) /= x_ghost(mb%o%lay)) then
          if (error(.true.,"ERROR: layout changes are not currently supported")) continue
        end if

        call glean(thy(mb))
        call glean(thy(lay))

        if (error("Exit multibasis_mod::update_mb")) continue

      end subroutine

      subroutine my_mb(mb)
!doc$ subroutine my(mb)
        type(multibasis_obj) :: mb

!cod$
        mb%ref = mb%ref + 1
        mb%o%ref = mb%o%ref + 1
      end subroutine

      subroutine my_new_mb(mbi,mb)
!doc$ subroutine my(mbi,mb)
        type(multibasis_obj) :: mbi, mb

!cod$
        mb%ref = 1
        mb%o => mbi%o
        mb%o%ref = mb%o%ref + 1
      end subroutine

      function thy_mb(mb) result(mbo)
!doc$ function thy(mb) result(mbo)
        type(multibasis_obj) :: mb, mbo

!cod$
        mb%ref = mb%ref - 1
        mb%o%ref = mb%o%ref - 1
        mbo%ref = mb%ref
        mbo%o => mb%o
      end function
   
      subroutine glean_mb(mb)
!doc$ subroutine glean(mb)
        type(multibasis_obj) :: mb

!cod$
        if (mb%o%ref < 1) then
          deallocate( mb%o%gridmap )
          deallocate( mb%o%gpt )
          deallocate( mb%o%first_band )
          deallocate( mb%o%last_band )
          deallocate( mb%o%my_band )
          deallocate( mb%o%band_participant)
          select case (mb%o%exx_comm_method)
          case (COLLECTIVE)
            deallocate( mb%o%first_spair )
            deallocate( mb%o%last_spair )
            deallocate( mb%o%my_spair )
            deallocate( mb%o%spair_index )
            deallocate( mb%o%spair_participant)
            deallocate( mb%o%first_lpair )
            deallocate( mb%o%last_lpair )
            deallocate( mb%o%my_lpair )
            deallocate( mb%o%lpair_index )
            deallocate( mb%o%lpair_participant)
          case (POINT_TO_POINT)
            if (associated( mb%o%r2r_dmap )) deallocate( mb%o%r2r_dmap )
            if (associated( mb%o%r2c_dmap )) deallocate( mb%o%r2c_dmap )
            if (associated( mb%o%c2c_smap )) deallocate( mb%o%c2c_smap )
            if (associated( mb%o%c2c_ulst )) deallocate( mb%o%c2c_ulst )
            if (associated( mb%o%c2r_cmap )) deallocate( mb%o%c2r_cmap )
            if (associated( mb%o%r2r_cmap )) deallocate( mb%o%r2r_cmap )
            if (associated( mb%o%rc_band_index )) deallocate( mb%o%rc_band_index )
          case (POINT_TO_POINT_2)
            if (associated( mb%o%r2r_dmap )) deallocate( mb%o%r2r_dmap )
            if (associated( mb%o%r2c_dmap )) deallocate( mb%o%r2c_dmap )
            if (associated( mb%o%c2c_smap )) deallocate( mb%o%c2c_smap )
            if (associated( mb%o%c2c_ulst )) deallocate( mb%o%c2c_ulst )
            if (associated( mb%o%c2r_cmap )) deallocate( mb%o%c2r_cmap )
            if (associated( mb%o%r2r_cmap )) deallocate( mb%o%r2r_cmap )
            if (associated( mb%o%rc_band_index )) deallocate( mb%o%rc_band_index )
          end select
          select case (mb%o%remap_type)
          case (MPI_REMAP)
            deallocate( mb%o%mpi_plan1%c1 )
            deallocate( mb%o%mpi_plan1%d1 )
            deallocate( mb%o%mpi_plan1%c2 )
            deallocate( mb%o%mpi_plan1%d2 )
            deallocate( mb%o%mpi_plan2%c1 )
            deallocate( mb%o%mpi_plan2%d1 )
            deallocate( mb%o%mpi_plan2%c2 )
            deallocate( mb%o%mpi_plan2%d2 )
            select case (mb%o%exx_comm_method)
            case (COLLECTIVE)
              deallocate( mb%o%mpi_plan3%c1 )
              deallocate( mb%o%mpi_plan3%d1 )
              deallocate( mb%o%mpi_plan3%c2 )
              deallocate( mb%o%mpi_plan3%d2 )
              deallocate( mb%o%mpi_plan4%c1 )
              deallocate( mb%o%mpi_plan4%d1 )
              deallocate( mb%o%mpi_plan4%c2 )
              deallocate( mb%o%mpi_plan4%d2 )
            end select
          case (SJP_REMAP)
            call destroy_sjp_remap_plan_i(mb%o%sjp_plan1)
            call destroy_sjp_remap_plan_i(mb%o%sjp_plan2)
            select case (mb%o%exx_comm_method)
            case (COLLECTIVE)
              call destroy_sjp_remap_plan_i(mb%o%sjp_plan3)
              call destroy_sjp_remap_plan_i(mb%o%sjp_plan4)
            end select
          end select
          call fft_destroy_serial_plan(mb%o%isplan1)
          select case (mb%o%usage)
          case (NORMAL)
            call fft_destroy_serial_plan(mb%o%isplan2)
          end select
          call glean(thy(mb%o%lay))
          deallocate( mb%o )
        end if
      end subroutine

      subroutine bequeath_mb(mb)
!doc$ subroutine bequeath(mb)
        type(multibasis_obj) :: mb

!cod$
        continue
      end subroutine

      subroutine assign_mb(mb,mb2)
!doc$ subroutine assignment(=)(mb,mb2)
        type(multibasis_obj), intent(inout) :: mb
        type(multibasis_obj), intent(in) :: mb2

!cod$
        type(multibasis_obj) :: mbt
        call my(mb2)
        mbt%o => mb%o
        mb%o%ref = mb%o%ref - mb%ref
        mb%o => mb2%o
        mb%o%ref = mb%o%ref + mb%ref
        call glean(mbt)
        call glean(thy(mb2))
      end subroutine
 
      function mb_ref(mb) result(r)
!doc$ function x_ref(mb) result(r)
        type(multibasis_obj) :: mb
        integer, dimension(2) :: r
!       effects: Returns mb%ref and mb%o%ref.

!cod$
        r(1) = mb%ref
        r(2) = mb%o%ref
        call glean(mb)
      end function

      function mb_ghost(mb) result(g)
!doc$ function x_ghost(mb) result(g)
        type(multibasis_obj) :: mb
        type(ghost) :: g
!       effects: Returns mb%o%g.

!cod$
        call my(mb)
        g = mb%o%g
        call glean(thy(mb))
      end function

      function mb_usage(mb) result(u)
!doc$ function x_usage(mb) result(u)
        integer :: u
        type(multibasis_obj) :: mb
!       effects: Returns mb%o%usage.
        
!cod$
        call my(mb)
        u = mb%o%usage
        call glean(thy(mb))
      end function

      function mb_n_bands(mb) result(n)
!doc$ function x_n_bands(mb) result(n)
        integer :: n
        type(multibasis_obj) :: mb
!       effects: Returns mb%o%nbands.
        
!cod$
        call my(mb)
        n = mb%o%nbands
        call glean(thy(mb))
      end function

      function mb_cutoff(mb) result(c)
!doc$ function x_cutoff(mb) result(c)
        real(double) :: c
        type(multibasis_obj) :: mb
!       effects: Returns mb%o%cutoff.
        
!cod$
        call my(mb)
        c = mb%o%cutoff
        call glean(thy(mb))
      end function

      function mb_n_gvectors(mb) result(n)
!doc$ function x_n_gvectors(mb) result(n)
        integer :: n
        type(multibasis_obj) :: mb
!       effects: Returns the total number of G vectors in the basis.
        
!cod$
        call my(mb)
        n = size(mb%o%gridmap,2)
        call glean(thy(mb))
      end function

      function mb_layout(mb) result(l)
!doc$ function x_layout(mb) result(l)
        type(multibasis_obj) :: mb
        type(layout_obj) :: l
!       effects: Returns a pointer to mb%o%layout.
!doc$
        call my(mb)
        call my(mb%o%lay,l)
        call bequeath(thy(l))
        call glean(thy(mb))
      end function

      subroutine band_remap_2d_to_1d_mb(mb,grp,data_in,data_out)
!doc$ subroutine band_remap(mb,grp,data_in,data_out)
        type(multibasis_obj) :: mb
        integer, intent(in) :: grp
        complex(double), dimension(:,:), intent(in) :: data_in
        complex(double), dimension(:), intent(out) :: data_out
!       requires: data_in be dimensions (size(mb%o%gpt,1),mpi_nprocs(KGROUP)) and data_out be dimension (size(mb%o%gridmap,2)).
!       effects: Copies data_in to data_out across processors.
!       errors: Passes errors.

!cod$
        complex(double) :: c

        call start_timer("multibasis: band_remap")

        call my(mb)

        select case (mb%o%remap_type)
        case (MPI_REMAP)
          if (grp < size(mb%o%first_band)) then
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan1%c2,mb%o%mpi_plan1%d2,data_out,mb%o%mpi_plan1%c1,mb%o%mpi_plan1%d1)
          else
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan2%c2,mb%o%mpi_plan2%d2,data_out,mb%o%mpi_plan2%c1,mb%o%mpi_plan2%d1)
          end if
        case (SJP_REMAP)
          if (grp < size(mb%o%first_band)) then
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan1%f)
          else
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan2%f)
          end if
        end select

100     call glean(thy(mb))

        if (error("Exit multibasis_mod::band_remap_2d_to_1d_mb")) continue

        if (.not.error()) call stop_timer("multibasis: band_remap")

      end subroutine

      subroutine band_remap_1d_to_2d_mb(mb,grp,data_in,data_out)
!doc$ subroutine band_remap(mb,grp,data_in,data_out)
        type(multibasis_obj) :: mb
        integer, intent(in) :: grp
        complex(double), dimension(:), intent(in) :: data_in
        complex(double), dimension(:,:), intent(out) :: data_out
!       requires: data_in be dimension (size(mb%o%gridmap,2)) and data_out be dimensions (size(mb%o%gpt,1),mpi_nprocs(KGROUP)).
!       effects: Copies data_in to data_out across processors.
!       errors: Passes errors.

!cod$
        complex(double) :: c

        call start_timer("multibasis: band_remap")

        call my(mb)

        select case (mb%o%remap_type)
        case (MPI_REMAP)
          if (grp < size(mb%o%first_band)) then
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan1%c1,mb%o%mpi_plan1%d1,data_out,mb%o%mpi_plan1%c2,mb%o%mpi_plan1%d2)
          else
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan2%c1,mb%o%mpi_plan2%d1,data_out,mb%o%mpi_plan2%c2,mb%o%mpi_plan2%d2)
          end if
        case (SJP_REMAP)
          if (grp < size(mb%o%first_band)) then
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan1%r)
          else
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan2%r)
          end if
        end select

100     call glean(thy(mb))

        if (error("Exit multibasis_mod::band_remap_1d_to_2d_mb")) continue

        if (.not.error()) call stop_timer("multibasis: band_remap")

      end subroutine

      subroutine spair_remap_2d_to_1d_mb(mb,grp,data_in,data_out)
!doc$ subroutine spair_remap(mb,grp,data_in,data_out)
        type(multibasis_obj) :: mb
        integer, intent(in) :: grp
        complex(double), dimension(:,:), intent(in) :: data_in
        complex(double), dimension(:), intent(out) :: data_out
!       requires: data_in be dimensions (size(mb%o%gpt,1),mpi_nprocs(KGROUP)) and data_out be dimension (size(mb%o%gridmap,2)).
!       effects: Copies data_in to data_out across processors.
!       errors: Passes errors.

!cod$
        complex(double) :: c

        call start_timer("multibasis: spair_remap")

        call my(mb)

        select case (mb%o%remap_type)
        case (MPI_REMAP)
          if (grp < size(mb%o%first_spair)) then
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan1%c2,mb%o%mpi_plan1%d2,data_out,mb%o%mpi_plan1%c1,mb%o%mpi_plan1%d1)
          else
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan3%c2,mb%o%mpi_plan3%d2,data_out,mb%o%mpi_plan3%c1,mb%o%mpi_plan3%d1)
          end if
        case (SJP_REMAP)
          if (grp < size(mb%o%first_spair)) then
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan1%f)
          else
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan3%f)
          end if
        end select

100     call glean(thy(mb))

        if (error("Exit multibasis_mod::spair_remap_2d_to_1d_mb")) continue

        if (.not.error()) call stop_timer("multibasis: spair_remap")

      end subroutine

      subroutine spair_remap_1d_to_2d_mb(mb,grp,data_in,data_out)
!doc$ subroutine spair_remap(mb,grp,data_in,data_out)
        type(multibasis_obj) :: mb
        integer, intent(in) :: grp
        complex(double), dimension(:), intent(in) :: data_in
        complex(double), dimension(:,:), intent(out) :: data_out
!       requires: data_in be dimension (size(mb%o%gridmap,2)) and data_out be dimensions (size(mb%o%gpt,1),mpi_nprocs(KGROUP)).
!       effects: Copies data_in to data_out across processors.
!       errors: Passes errors.

!cod$
        complex(double) :: c

        call start_timer("multibasis: spair_remap")

        call my(mb)

        select case (mb%o%remap_type)
        case (MPI_REMAP)
          if (grp < size(mb%o%first_spair)) then
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan1%c1,mb%o%mpi_plan1%d1,data_out,mb%o%mpi_plan1%c2,mb%o%mpi_plan1%d2)
          else
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan3%c1,mb%o%mpi_plan3%d1,data_out,mb%o%mpi_plan3%c2,mb%o%mpi_plan3%d2)
          end if
        case (SJP_REMAP)
          if (grp < size(mb%o%first_spair)) then
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan1%r)
          else
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan3%r)
          end if
        end select

100     call glean(thy(mb))

        if (error("Exit multibasis_mod::spair_remap_1d_to_2d_mb")) continue

        if (.not.error()) call stop_timer("multibasis: spair_remap")

      end subroutine

      subroutine lpair_remap_2d_to_1d_mb(mb,grp,data_in,data_out)
!doc$ subroutine lpair_remap(mb,grp,data_in,data_out)
        type(multibasis_obj) :: mb
        integer, intent(in) :: grp
        complex(double), dimension(:,:), intent(in) :: data_in
        complex(double), dimension(:), intent(out) :: data_out
!       requires: data_in be dimensions (size(mb%o%gpt,1),mpi_nprocs(KGROUP)) and data_out be dimension (size(mb%o%gridmap,2)).
!       effects: Copies data_in to data_out across processors.
!       errors: Passes errors.

!cod$
        complex(double) :: c

        call start_timer("multibasis: lpair_remap")

        call my(mb)

        select case (mb%o%remap_type)
        case (MPI_REMAP)
          if (grp < size(mb%o%first_lpair)) then
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan1%c2,mb%o%mpi_plan1%d2,data_out,mb%o%mpi_plan1%c1,mb%o%mpi_plan1%d1)
          else
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan4%c2,mb%o%mpi_plan4%d2,data_out,mb%o%mpi_plan4%c1,mb%o%mpi_plan4%d1)
          end if
        case (SJP_REMAP)
          if (grp < size(mb%o%first_lpair)) then
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan1%f)
          else
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan4%f)
          end if
        end select

100     call glean(thy(mb))

        if (error("Exit multibasis_mod::lpair_remap_2d_to_1d_mb")) continue

        if (.not.error()) call stop_timer("multibasis: lpair_remap")

      end subroutine

      subroutine lpair_remap_1d_to_2d_mb(mb,grp,data_in,data_out)
!doc$ subroutine lpair_remap(mb,grp,data_in,data_out)
        type(multibasis_obj) :: mb
        integer, intent(in) :: grp
        complex(double), dimension(:), intent(in) :: data_in
        complex(double), dimension(:,:), intent(out) :: data_out
!       requires: data_in be dimension (size(mb%o%gridmap,2)) and data_out be dimensions (size(mb%o%gpt,1),mpi_nprocs(KGROUP)).
!       effects: Copies data_in to data_out across processors.
!       errors: Passes errors.

!cod$
        complex(double) :: c

        call start_timer("multibasis: lpair_remap")

        call my(mb)

        select case (mb%o%remap_type)
        case (MPI_REMAP)
          if (grp < size(mb%o%first_lpair)) then
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan1%c1,mb%o%mpi_plan1%d1,data_out,mb%o%mpi_plan1%c2,mb%o%mpi_plan1%d2)
          else
            call alltoallv(KGROUP,data_in,mb%o%mpi_plan4%c1,mb%o%mpi_plan4%d1,data_out,mb%o%mpi_plan4%c2,mb%o%mpi_plan4%d2)
          end if
        case (SJP_REMAP)
          if (grp < size(mb%o%first_lpair)) then
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan1%r)
          else
            call remap_2d(data_in,data_out,c,mb%o%sjp_plan4%r)
          end if
        end select

100     call glean(thy(mb))

        if (error("Exit multibasis_mod::lpair_remap_1d_to_2d_mb")) continue

        if (.not.error()) call stop_timer("multibasis: lpair_remap")

      end subroutine

      subroutine form_vmap_mb(mb,cutoff,vmap)
!doc$ subroutine_form_vmap(mb,cutoff,vmap)
        type(multibasis_obj) :: mb
        real(double) :: cutoff
        integer, dimension(:), pointer :: vmap
!       requires: cutoff be > x_cutoff(mb%o%lay). vmap be nullified.
!       effects: Returns a mapping of packed G-vectors from cutoff to mb%o%cutoff.
        
!cod$
        integer :: i1, i2, i3, ng, ng2
        real(double), dimension(:,:,:), pointer :: cp, fx, fy, fz

        call my(mb)

        nullify( cp, fx, fy, fz )

        call fmesh(fx,fy,fz,mb%o%lay,S_TYPE)
        call alloc(cp,mb%o%lay,S_TYPE)

        do i3 = 1,size(cp,3)
        do i2 = 1,size(cp,2)
        do i1 = 1,size(cp,1)
          cp(i1,i2,i3) = (fx(i1,i2,i3) + mb%o%kpt(1))**2 + (fy(i1,i2,i3) + mb%o%kpt(2))**2 + (fz(i1,i2,i3) + mb%o%kpt(3))**2
        end do
        end do
        end do

        ng = size(mb%o%gridmap,2)

        allocate( vmap(ng) )

        if (cutoff <= mb%o%cutoff) then

          ng = 0
          ng2 = 0
          do i3 = 1,size(cp,3)
          do i2 = 1,size(cp,2)
          do i1 = 1,size(cp,1)
            if (cp(i1,i2,i3) < mb%o%cutoff) then
              ng = ng + 1
              if (cp(i1,i2,i3) < cutoff) then
                ng2 = ng2 + 1
                vmap(ng) = ng2
              else
                vmap(ng) = 0
              end if
            end if
          end do
          end do
          end do

        elseif (cutoff > mb%o%cutoff) then

          ng = 0
          ng2 = 0
          do i3 = 1,size(cp,3)
          do i2 = 1,size(cp,2)
          do i1 = 1,size(cp,1)
            if (cp(i1,i2,i3) < cutoff) then
              ng2 = ng2 + 1
              if (cp(i1,i2,i3) < mb%o%cutoff) then
                ng = ng + 1
                vmap(ng) = ng2
              end if
            end if
          end do
          end do
          end do

        end if

100     if (associated( cp )) deallocate( cp )
        if (associated( fx )) deallocate( fx )
        if (associated( fy )) deallocate( fy )
        if (associated( fz )) deallocate( fz )

        call glean(thy(mb))

        if (error("Exit multibasis_mod::form_vmap_mb")) continue

      end subroutine

      function mb_wormhole(mb) result(hole)
!doc$ function wormhole(mb) result(hole)
        type(multibasis_obj) :: mb
        type(multibasis_rep), pointer :: hole      
!       effects: Points hole to the internal rep of mb.
!       errors: Wormhole is an implementation dependent thing and you should know what you're doing.

!cod$
        call my(mb)
        hole => mb%o
        call glean(thy(mb))
      end function

      subroutine diary_mb(mb)
!doc$ subroutine diary(mb)
        type(multibasis_obj) :: mb
!       effects: Writes mb information to the diary.

!cod$
        call my(mb)
        if (i_access(diaryfile())) then
          select case (mb%o%remap_type)
          case (MPI_REMAP)
            write(x_unit(diaryfile()),'(/,t4,"MPI remap")')
          case (SJP_REMAP)
            write(x_unit(diaryfile()),'(/,t4,"SJP remap")')
          end select
        end if
        call glean(thy(mb))
      end subroutine

! private routines

      subroutine own_i(mb)
        type(multibasis_obj) :: mb, mbt
        if (mb%ref < mb%o%ref) then
          allocate( mbt%o )
          mbt%o%ref = 0
          mbt%o%g = mb%o%g
          mbt%o%usage = mb%o%usage
          mbt%o%cutoff = mb%o%cutoff
          mbt%o%cellvol = mb%o%cellvol
          mbt%o%kpt = mb%o%kpt
          allocate( mbt%o%gpt(size(mb%o%gpt,1),size(mb%o%gpt,2)) )                    ;  mbt%o%gpt = mb%o%gpt
          allocate( mbt%o%gridmap(size(mb%o%gridmap,1),size(mb%o%gridmap,2)) )        ;  mbt%o%gridmap = mb%o%gridmap
          mbt%o%first_g = mb%o%first_g
          mbt%o%last_g = mb%o%last_g
          mbt%o%nbands = mb%o%nbands
          allocate( mbt%o%first_band(size(mb%o%first_band)) )                 ; mbt%o%first_band = mb%o%first_band
          allocate( mbt%o%last_band(size(mb%o%last_band)) )                   ; mbt%o%last_band = mb%o%last_band
          allocate( mbt%o%my_band(size(mb%o%my_band)) )                       ; mbt%o%my_band = mb%o%my_band
          allocate( mbt%o%band_participant(size(mb%o%band_participant)) )     ; mbt%o%band_participant = mb%o%band_participant
          mbt%o%exx_comm_method = mb%o%exx_comm_method
          select case (mbt%o%exx_comm_method)
          case (COLLECTIVE)
            allocate( mbt%o%first_spair(size(mb%o%first_spair)) )               ; mbt%o%first_spair = mb%o%first_spair
            allocate( mbt%o%last_spair(size(mb%o%last_spair)) )                 ; mbt%o%last_spair = mb%o%last_spair
            allocate( mbt%o%my_spair(size(mb%o%my_spair)) )                     ; mbt%o%my_spair = mb%o%my_spair
            allocate( mbt%o%spair_index(2,size(mb%o%spair_index)) )             ; mbt%o%spair_index = mb%o%spair_index
            allocate( mbt%o%spair_participant(size(mb%o%spair_participant)) )   ; mbt%o%spair_participant = mb%o%spair_participant
            allocate( mbt%o%first_lpair(size(mb%o%first_lpair)) )               ; mbt%o%first_lpair = mb%o%first_lpair
            allocate( mbt%o%last_lpair(size(mb%o%last_lpair)) )                 ; mbt%o%last_lpair = mb%o%last_lpair
            allocate( mbt%o%my_lpair(size(mb%o%my_lpair)) )                     ; mbt%o%my_lpair = mb%o%my_lpair
            allocate( mbt%o%lpair_index(2,size(mb%o%lpair_index)) )             ; mbt%o%lpair_index = mb%o%lpair_index
            allocate( mbt%o%lpair_participant(size(mb%o%lpair_participant)) )   ; mbt%o%lpair_participant = mb%o%lpair_participant
          case (POINT_TO_POINT)
            if (associated( mb%o%r2r_dmap )) then
              allocate( mbt%o%r2r_dmap(size(mb%o%r2r_dmap,1),2) )               ; mbt%o%r2r_dmap = mb%o%r2r_dmap
            end if
            if (associated( mb%o%r2c_dmap )) then
              allocate( mbt%o%r2c_dmap(size(mb%o%r2c_dmap,1),2) )               ; mbt%o%r2c_dmap = mb%o%r2c_dmap
            end if
            if (associated( mb%o%c2c_smap )) then
              allocate( mbt%o%c2c_smap(size(mb%o%c2c_smap,1),2) )               ; mbt%o%c2c_smap = mb%o%c2c_smap
              allocate( mbt%o%c2c_ulst(size(mb%o%c2c_ulst)) )                   ; mbt%o%c2c_ulst = mb%o%c2c_ulst
            end if
            if (associated( mb%o%c2r_cmap )) then
              allocate( mbt%o%c2r_cmap(size(mb%o%c2r_cmap,1),2) )               ; mbt%o%c2r_cmap = mb%o%c2r_cmap
            end if
            if (associated( mb%o%r2r_cmap )) then
              allocate( mbt%o%r2r_cmap(size(mb%o%r2r_cmap,1),2) )               ; mbt%o%r2r_cmap = mb%o%r2r_cmap
            end if
            if (associated( mb%o%rc_band_index )) then
              allocate( mbt%o%rc_band_index(size(mb%o%rc_band_index,1),2) )     ; mbt%o%rc_band_index = mb%o%rc_band_index
            end if
          case (POINT_TO_POINT_2)
            if (associated( mb%o%r2r_dmap )) then
              allocate( mbt%o%r2r_dmap(size(mb%o%r2r_dmap,1),2) )               ; mbt%o%r2r_dmap = mb%o%r2r_dmap
            end if
            if (associated( mb%o%r2c_dmap )) then
              allocate( mbt%o%r2c_dmap(size(mb%o%r2c_dmap,1),2) )               ; mbt%o%r2c_dmap = mb%o%r2c_dmap
            end if
            if (associated( mb%o%c2c_smap )) then
              allocate( mbt%o%c2c_smap(size(mb%o%c2c_smap,1),2) )               ; mbt%o%c2c_smap = mb%o%c2c_smap
              allocate( mbt%o%c2c_ulst(size(mb%o%c2c_ulst)) )                   ; mbt%o%c2c_ulst = mb%o%c2c_ulst
            end if
            if (associated( mb%o%c2r_cmap )) then
              allocate( mbt%o%c2r_cmap(size(mb%o%c2r_cmap,1),2) )               ; mbt%o%c2r_cmap = mb%o%c2r_cmap
            end if
            if (associated( mb%o%r2r_cmap )) then
              allocate( mbt%o%r2r_cmap(size(mb%o%r2r_cmap,1),2) )               ; mbt%o%r2r_cmap = mb%o%r2r_cmap
            end if
            if (associated( mb%o%rc_band_index )) then
              allocate( mbt%o%rc_band_index(size(mb%o%rc_band_index,1),2) )     ; mbt%o%rc_band_index = mb%o%rc_band_index
            end if
          end select
          mbt%o%remap_type = mb%o%remap_type
          select case (mbt%o%remap_type)
          case (MPI_REMAP)
            allocate( mbt%o%mpi_plan1%c1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan1%c1 = mb%o%mpi_plan1%c1
            allocate( mbt%o%mpi_plan1%d1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan1%d1 = mb%o%mpi_plan1%d1
            allocate( mbt%o%mpi_plan1%c2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan1%c2 = mb%o%mpi_plan1%c2
            allocate( mbt%o%mpi_plan1%d2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan1%d2 = mb%o%mpi_plan1%d2
            allocate( mbt%o%mpi_plan2%c1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan2%c1 = mb%o%mpi_plan2%c1
            allocate( mbt%o%mpi_plan2%d1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan2%d1 = mb%o%mpi_plan2%d1
            allocate( mbt%o%mpi_plan2%c2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan2%c2 = mb%o%mpi_plan2%c2
            allocate( mbt%o%mpi_plan2%d2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan2%d2 = mb%o%mpi_plan2%d2
            select case (mb%o%exx_comm_method)
            case (COLLECTIVE)
              allocate( mbt%o%mpi_plan3%c1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan3%c1 = mb%o%mpi_plan3%c1
              allocate( mbt%o%mpi_plan3%d1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan3%d1 = mb%o%mpi_plan3%d1
              allocate( mbt%o%mpi_plan3%c2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan3%c2 = mb%o%mpi_plan3%c2
              allocate( mbt%o%mpi_plan3%d2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan3%d2 = mb%o%mpi_plan3%d2
              allocate( mbt%o%mpi_plan4%c1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan4%c1 = mb%o%mpi_plan4%c1
              allocate( mbt%o%mpi_plan4%d1(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan4%d1 = mb%o%mpi_plan4%d1
              allocate( mbt%o%mpi_plan4%c2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan4%c2 = mb%o%mpi_plan4%c2
              allocate( mbt%o%mpi_plan4%d2(0:mpi_nprocs(KGROUP)-1) ) ; mbt%o%mpi_plan4%d2 = mb%o%mpi_plan4%d2
            end select
          case (SJP_REMAP)
            mbt%o%sjp_plan1 = mb%o%sjp_plan1
            mbt%o%sjp_plan2 = mb%o%sjp_plan2
            call create_sjp_remap_plan_i(mb%o%sjp_plan1)
            call create_sjp_remap_plan_i(mb%o%sjp_plan2)
            select case (mb%o%exx_comm_method)
            case (COLLECTIVE)
              mbt%o%sjp_plan3 = mb%o%sjp_plan3
              mbt%o%sjp_plan4 = mb%o%sjp_plan4
              call create_sjp_remap_plan_i(mb%o%sjp_plan3)
              call create_sjp_remap_plan_i(mb%o%sjp_plan4)
            end select
          end select
          mbt%o%isplan1 = mb%o%isplan1
          select case (mbt%o%usage)
          case (NORMAL)
            mbt%o%isplan2 = mb%o%isplan2
          end select
          call my(mb%o%lay,mbt%o%lay)
          mb%o%ref = mb%o%ref - mb%ref
          mb%o => mbt%o
          mb%o%ref = mb%o%ref + mb%ref
        end if
      end subroutine

      subroutine create_sjp_remap_plan_i(rm)
        type(sjp_remap_plan) :: rm
        integer, parameter :: nqty = 2, permute = 0, memory = 1, precision = 2
        integer :: comm
        comm = mpi_comm(KGROUP)
        call remap_2d_create_plan(comm,rm%ilo_in,rm%ihi_in,rm%jlo_in,rm%jhi_in, &
                                       rm%ilo_out,rm%ihi_out,rm%jlo_out,rm%jhi_out,nqty,permute,memory,precision,rm%f)
        call remap_2d_create_plan(comm,rm%ilo_out,rm%ihi_out,rm%jlo_out,rm%jhi_out, &
                                       rm%ilo_in,rm%ihi_in,rm%jlo_in,rm%jhi_in,nqty,permute,memory,precision,rm%r)
      end subroutine

      subroutine destroy_sjp_remap_plan_i(rm)
        type(sjp_remap_plan) :: rm
        call remap_2d_destroy_plan(rm%f)
        call remap_2d_destroy_plan(rm%r)
      end subroutine

      end module

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module config_fh_mod
!doc$ module config_fh_mod

!     One datatype is available here: type(config_fh_obj)

!     config_fh_mod encapsulates a fixed-hamiltonian electronic-structure solution for a configuration of atoms.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use arg_mod
      use tagio_mod
      use ghost_mod
      use diary_mod
      use electrons_fh_mod
      use crystal_mod
      use layout_mod  ! Needed for PGI compiler
      use lattice_mod
      use atoms_mod
      use external_mod
      use fields_fh_mod
      use interrupt_mod
      use timing_mod

!cod$
      implicit none
      private

      integer, parameter :: NONE          = 0
      integer, parameter :: WAVEFUNCTIONS = 1

      type :: config_fh_rep
        integer :: ref
        type(ghost) :: g
        integer :: cvg_mode                             ! mode used to determine convergence
        integer :: max_steps                            ! maximum number of steps in the convergence loop
        type(external_obj) :: external                  ! external object
        type(fields_fh_obj) :: fields                   ! fields object
        type(electrons_fh_obj) :: electrons             ! electrons object
      end type

      type, public :: config_fh_obj
        private
        integer :: ref
        type(config_fh_rep), pointer :: o
      end type

!doc$
      public :: config_fh
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_external
      public :: x_electrons
      public :: diary
      public :: decompose

!cod$
      interface config_fh
         module procedure constructor_cfg
      end interface
      interface my
         module procedure my_cfg, my_new_cfg
      end interface
      interface thy
         module procedure thy_cfg
      end interface
      interface glean
         module procedure glean_cfg
      end interface
      interface bequeath
         module procedure bequeath_cfg
      end interface
      interface assignment (=)
         module procedure assign_cfg
      end interface
      interface x_ref
         module procedure cfg_ref
      end interface
      interface x_ghost
         module procedure cfg_ghost
      end interface
      interface x_external
         module procedure cfg_external
      end interface
      interface x_electrons
         module procedure cfg_electrons
      end interface
      interface diary
         module procedure diary_cfg
      end interface
      interface decompose
         module procedure decompose_cfg
      end interface
      
      contains

! public routines

      function constructor_cfg() result(cfg)
!doc$ function config_fh() result(cfg)
        type(config_fh_obj) :: cfg
!       effects: Constructs a new cfg.
!       errors: Restart tag not recognized. Passes errors.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: mode, tag
        integer :: r_nsg
        integer(long) :: dsize, iosl, ndata, s4
        real(double) :: version
        type(tagio_obj) :: restf

        call start_timer("config_fh: constructor")

        cfg%ref = 0
        allocate( cfg%o )
        cfg%o%ref = 0
        cfg%o%g = x_ghost()

        ! read the electronic convergence criteria
        call arglc("config_convergence",tag,found)
        if (.not.found) tag = "wavefunctions"
        select case (trim(tag))
        case ("none")
          cfg%o%cvg_mode = NONE
        case ("wavefunctions")
          cfg%o%cvg_mode = WAVEFUNCTIONS
        case default
          if (error(.true.,"ERROR: config_convergence is not recognized")) goto 200
        end select
        call arg("config_steps",cfg%o%max_steps,found)
        if (.not.found) cfg%o%max_steps = 40
        if (error(cfg%o%max_steps < 0,"ERROR: config_steps < 0")) goto 200

        ! read the restart mode
        call arglc("restart",mode,found)
        if (error(.not.found,"ERROR: restart tag was not found")) goto 200

        ! open the restart file and check for the correct type (mkey)
        call my(tagio(trim(restart_path),TAGIO_READ,mkey,len(mkey)),restf)
        if (i_access(restf)) iosl = x_tagfd(restf)
        if (i_comm(restf)) call broadcast(FILE_SCOPE,iosl)
        if (error(iosl == 0,"ERROR: restart file was not found")) goto 200

        ! check the version of the restart file
        if (i_access(restf)) tios = findfirsttag(restf,"VERSION")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: VERSION tag was not found")) goto 100
        if (i_access(restf)) then
          dsize = sizeof_double ; ndata = 1
          call readf(version,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,version)
        if (error(version /= es_version,"ERROR: incorrect version of the restart file")) goto 100

        ! read the number of sgroups
        if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_SGROUPS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (tios == TAG_NOT_FOUND) then
          call warn("Warning: NUMBER_OF_SGROUPS tag was not found - using the input value")
          if (i_access(restf)) call rewind_tobeginning(restf)
        else
          if (i_access(restf)) then
            dsize = sizeof_long
            ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            r_nsg = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_nsg)
          if (error(r_nsg /= mpi_nsgroups(),"ERROR: different numbers of sgroups")) goto 100
        end if

        ! construct the external, fields, and electrons
        select case (trim(mode))
        case ("f")
          call my(external(),cfg%o%external)                    ; if (error()) goto 100
          call my(fields_fh(cfg%o%external,restf),cfg%o%fields) ; if (error()) goto 100
          call my(electrons_fh(cfg%o%external,x_potential(cfg%o%fields)),cfg%o%electrons)
        case ("fe")
          call my(external(),cfg%o%external)                    ; if (error()) goto 100
          call my(fields_fh(cfg%o%external,restf),cfg%o%fields) ; if (error()) goto 100
          call my(electrons_fh(cfg%o%external,x_potential(cfg%o%fields),restf),cfg%o%electrons)
        case default
          if (error(.true.,"ERROR: restart mode was not recognized")) continue
        end select

100     call glean(thy(restf))

        call sync_configuration_errors() ; if (error()) goto 200

        call diary_construction_i(cfg%o)

        ! iterate to convergence
        call iterator_i(cfg%o)
        call sync_configuration_errors()

200     if (error("Exit config_fh_mod::constructor_cfg")) continue

        if (.not.error()) call stop_timer("config_fh: constructor")

      end function

      subroutine my_cfg(cfg)
!doc$ subroutine my(cfg)
        type(config_fh_obj) :: cfg

!cod$
        cfg%ref = cfg%ref + 1
        cfg%o%ref = cfg%o%ref + 1
      end subroutine

      subroutine my_new_cfg(cfgi,cfg)
!doc$ subroutine my(cfgi,cfg)
        type(config_fh_obj) :: cfgi, cfg

!cod$
        cfg%ref = 1
        cfg%o => cfgi%o
        cfg%o%ref = cfg%o%ref + 1
      end subroutine

      function thy_cfg(cfg) result(cfgo)
!doc$ function thy(cfg) result(cfgo)
        type(config_fh_obj) :: cfg, cfgo

!cod$
        cfg%ref = cfg%ref - 1
        cfg%o%ref = cfg%o%ref - 1
        cfgo%ref = cfg%ref
        cfgo%o => cfg%o
      end function

      subroutine glean_cfg(cfg)
!doc$ subroutine glean(cfg)
        type(config_fh_obj) :: cfg

!cod$
        if (cfg%o%ref < 1) then
          call glean(thy(cfg%o%external))
          call glean(thy(cfg%o%fields))
          call glean(thy(cfg%o%electrons))
          deallocate( cfg%o )
        end if
      end subroutine

      subroutine bequeath_cfg(cfg)
!doc$ subroutine bequeath(cfg)
        type(config_fh_obj) :: cfg

!cod$
        continue
      end subroutine

      subroutine assign_cfg(cfg,cfg2)
!doc$ subroutine assignment(=)(cfg,cfg2)
        type(config_fh_obj), intent(inout) :: cfg
        type(config_fh_obj), intent(in) :: cfg2

!cod$
        type(config_fh_obj) :: cfgt
        call my(cfg2)
        cfgt%o => cfg%o
        cfg%o%ref = cfg%o%ref - cfg%ref
        cfg%o => cfg2%o
        cfg%o%ref = cfg%o%ref + cfg%ref
        call glean(cfgt)
        call glean(thy(cfg2))
      end subroutine

      function cfg_ref(cfg) result(r)
!doc$ function x_ref(cfg) result(r)
        type(config_fh_obj) :: cfg
        integer, dimension(2) :: r
!       effects: Returns cfg%ref and cfg%o%ref.

!cod$
        r(1) = cfg%ref
        r(2) = cfg%o%ref
        call glean(cfg)
      end function

      function cfg_ghost(cfg) result(g)
!doc$ function x_ghost(cfg) result(g)
        type(config_fh_obj) :: cfg
        type(ghost) :: g

!cod$
        call my(cfg)
        g = cfg%o%g
        call glean(thy(cfg))
      end function

      function cfg_external(cfg) result(ext)
!doc$ function x_external(cfg) result(ext)
        type(config_fh_obj) :: cfg
        type(external_obj) :: ext
!       effects: Returns the external component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%external,ext)
        call glean(thy(cfg))
        call bequeath(thy(ext))
      end function

      function cfg_electrons(cfg) result(el)
!doc$ function x_electrons(cfg) result(el)
        type(config_fh_obj) :: cfg
        type(electrons_fh_obj) :: el
!       effects: Returns the electrons component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%electrons,el)
        call glean(thy(cfg))
        call bequeath(thy(el))
      end function

      subroutine diary_cfg(cfg)
!doc$ subroutine diary(cfg)
        type(config_fh_obj) :: cfg
!       modifies: Output stream
!       effects: Prints cfg information.
!       errors: Passes errors.

!cod$
        call my(cfg)
        call diary(cfg%o%electrons)
        call glean(thy(cfg))
        if (error("Exit config_fh_mod::diary_cfg")) continue
      end subroutine

      subroutine decompose_cfg(cfg)
!doc$ subroutine decompose(cfg)
        type(config_fh_obj) :: cfg
!       effects: Decomposes the Kohn-Sham functions into s, p, & d spherical harmonics around user-defined and atom sites.
!       errors: decomposition tag not recognized. Passes errors.

!cod$
        logical :: exist_file, found
        character(tag_sz) :: type
        character(line_len) :: mode, switch
        character(14+tag_sz) :: dr_tag
        integer :: ia, ios, is, na, ns, nu
        real(double) :: radius
        real(double), dimension(3) :: pos_lat, pos_xyz
        real(double), dimension(:,:), allocatable :: site_data
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats
        type(file_obj) :: f

        call my(cfg)

        call arglc("decomposition",switch,found)
        if (.not.found) switch = "off"
        select case (trim(switch))
        case ("on")
          continue
        case ("off")
          goto 700
        case default
          if (error(.true.,"ERROR: decomposition tag is not recognized")) goto 700
        end select

        call start_timer("config_fh: decompose")

        call arglc("dcomp_mode",mode,found)
        if (.not.found) mode = "l"
        select case (trim(mode))
        case ("l","lm","xyz")
          continue
        case default
          if (error(.true.,"ERROR: dcomp_mode tag is not recognized")) goto 600
        end select

        call my(x_lattice(x_crystal(cfg%o%external)),lat)
        call my(x_atoms(x_crystal(cfg%o%external)),ats)

        call my(file(trim(dsites_path)),f)
        if (i_access(f)) inquire(file=x_name(f),exist=exist_file)
        if (i_comm(f)) call broadcast(FILE_SCOPE,exist_file)

        nu = 0
        if (exist_file) then
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to open dsites file")) goto 200
          if (i_access(f)) read(x_unit(f),*,iostat=ios) nu
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to read the number of sites")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,nu)
        end if
        na = x_n_atoms(ats)
        ns = na + nu
        allocate( site_data(4,ns) )

        do is = 1,nu
          if (i_access(f)) read(x_unit(f),*,iostat=ios) pos_lat, radius
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to read dsites data")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,pos_lat)
          site_data(1:3,is) = lat2r(lat,pos_lat)
          if (i_comm(f)) call broadcast(FILE_SCOPE,radius)
          site_data(4,is) = radius
        end do
        do is = nu+1,ns
          ia = is - nu
          site_data(1:3,is) = lat2r(lat,x_position(ats,ia))
          type = x_type(ats,ia)
          dr_tag = "dcomp_radius_"//type
          call arg(trim(dr_tag),radius,found)
          if (found) then
            site_data(4,is) = radius
          else
            call warn("dcomp_radius tag was not found - using the value 3.0")
            site_data(4,is) = 3.0_double
          end if
        end do

100     if (i_access(f)) close(x_unit(f))
200     call glean(thy(f)) ; if (error()) goto 500

        call my(file(trim(dcomp_path)),f)
        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open dcomp file")) goto 400

        if (i_access(f)) then
          select case (trim(mode))
          case ("L", "l")
            write(x_unit(f),'("Kohn-Sham function decomposition in l spherical harmonics:")')
          case ("LM", "lm")
            write(x_unit(f),'("Kohn-Sham function decomposition in lm spherical harmonics:")')
          case ("XYZ", "xyz")
            write(x_unit(f),'("Kohn-Sham function decomposition in xyz spherical harmonics:")')
          end select
          write(x_unit(f),'(/,t2,"Site information:")')
          write(x_unit(f),'(/,t4,"site",4x," type ",10x,"x",10x,"y",10x,"z",11x,"a1",9x,"a2",9x,"a3",10x,"radius")')
          write(x_unit(f),'(t3,14("-"),5x,32("-"),3x,31("-"),5x,9("-"))')
          do is = 1,nu
            type = "User"
            pos_xyz = site_data(1:3,is)
            pos_lat = r2lat(lat,pos_xyz)
            radius = site_data(4,is)
            write(x_unit(f),'(t5,i3,1x,a8,3x,2(1x,3(f11.5)),f14.5)') is, trim(type), pos_xyz, pos_lat, radius
          end do
          do is = nu+1,ns
            ia = is - nu
            type = x_type(ats,ia)
            pos_xyz = site_data(1:3,is)
            pos_lat = r2lat(lat,pos_xyz)
            radius = site_data(4,is)
            write(x_unit(f),'(t5,i3,1x,a8,3x,2(1x,3(f11.5)),f14.5)') is, trim(type), pos_xyz, pos_lat, radius
          end do
          write(x_unit(f),'(/,t2,"Decompositions:")')
        end if

        call decompose(cfg%o%electrons,site_data,mode,f) ; if (error()) goto 300

300     if (i_access(f)) close(x_unit(f))
400     call glean(thy(f)) ; if (error()) goto 500

500     if (allocated( site_data )) deallocate( site_data )
        call glean(thy(lat))
        call glean(thy(ats))

600     if (.not.error()) call stop_timer("config_fh: decompose")

700     call glean(thy(cfg))

        if (error("Exit config_fh_mod::decompose_cfg")) continue

      end subroutine

! private routines

      subroutine diary_construction_i(cfgr)
        type(config_fh_rep) :: cfgr

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,"Config object construction:")')
          write(x_unit(diaryfile()),'(/,t4,"Fixed-hamiltonian calculation")')
          select case (cfgr%cvg_mode)
          case (NONE)
            write(x_unit(diaryfile()),'(/,t4,"Convergence will not be checked")')
          case (WAVEFUNCTIONS)
            write(x_unit(diaryfile()),'(/,t4,"Convergence will be determined from the wavefunctions residual")')
          end select
        end if

      end subroutine

      subroutine iterator_i(cfgr)
        type(config_fh_rep) :: cfgr

        logical :: done
        integer :: lc

        lc = 0
        do
          if (error(user_abort(),"USER INITIATED ABORT")) goto 100
          lc = lc + 1
          call check_convergence_i(cfgr,lc,done)
          if (done) exit
          call update(cfgr%electrons) ; if (error()) goto 100
        end do

100     if (error("Exit config_fh_mod::iterator_i")) continue

      end subroutine

      subroutine check_convergence_i(cfgr,lc,done)
        type(config_fh_rep) :: cfgr
        integer, intent(in) :: lc
        logical, intent(out) :: done
        real(double) :: rn
        done = .false.
        select case (cfgr%cvg_mode)
        case (NONE)
          rn = x_residual_norm(cfgr%electrons)
          if (i_access(diaryfile())) then
            if (lc == 1) then
              write(x_unit(diaryfile()),'(/,t4,"Fixed-hamiltonian step  ",i1,":  wavefunctions residual = ",es10.4)') lc, rn
            elseif (lc < 10) then
              write(x_unit(diaryfile()),'(  t4,"Fixed-hamiltonian step  ",i1,":  wavefunctions residual = ",es10.4)') lc, rn
            else
              write(x_unit(diaryfile()),'(  t4,"Fixed-hamiltonian step ", i0,":  wavefunctions residual = ",es10.4)') lc, rn
            end if
          end if
          if (i_access(output)) then
            if (lc < 10) then
              write(x_unit(output),'("Fixed-hamiltonian step  ",i1,":  wavefunctions residual = ",es10.4)') lc, rn
            else
              write(x_unit(output),'("Fixed-hamiltonian step ", i0,":  wavefunctions residual = ",es10.4)') lc, rn
            end if
          end if
        case (WAVEFUNCTIONS)
          rn = x_residual_norm(cfgr%electrons)
          if (i_access(diaryfile())) then
            if (lc == 1) then
              write(x_unit(diaryfile()),'(/,t4,"Fixed-hamiltonian step  ",i1,":  wavefunctions residual = ",es10.4)') lc, rn
            elseif (lc < 10) then
              write(x_unit(diaryfile()),'(  t4,"Fixed-hamiltonian step  ",i1,":  wavefunctions residual = ",es10.4)') lc, rn
            else
              write(x_unit(diaryfile()),'(  t4,"Fixed-hamiltonian step ", i0,":  wavefunctions residual = ",es10.4)') lc, rn
            end if
          end if
          if (i_access(output)) then
            if (lc < 10) then
              write(x_unit(output),'("Fixed-hamiltonian step  ",i1,":  wavefunctions residual = ",es10.4)') lc, rn
            else
              write(x_unit(output),'("Fixed-hamiltonian step ", i0,":  wavefunctions residual = ",es10.4)') lc, rn
            end if
          end if
          done = x_converged(cfgr%electrons)
        end select
        if (lc >= cfgr%max_steps) done = .true.
      end subroutine

      end module

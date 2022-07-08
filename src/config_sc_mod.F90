!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module config_sc_mod
!doc$ module config_sc_mod

!     One datatype is available here: type(config_sc_obj)

!     Config_sc_mod encapsulates a self-consistent electronic structure solution for a configuration of atoms.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use arg_mod
      use tagio_mod
      use ghost_mod
      use diary_mod
      use electrons_sc_mod
      use crystal_mod
      use layout_mod
      use lattice_mod
      use atoms_mod
      use external_mod
      use fields_sc_mod
      use interrupt_mod
      use timing_mod

!cod$
      implicit none
      private

      integer, parameter :: NONE          = 0
      integer, parameter :: ENERGY        = 1
      integer, parameter :: DENSITY       = 2
      integer, parameter :: WAVEFUNCTIONS = 3

      type :: config_sc_rep
        integer :: ref
        type(ghost) :: g
        integer :: cvg_mode                               ! mode used to determine convergence
        integer :: max_steps                              ! maximum number of steps in the convergence loop
        real(double) :: energy_tol                        ! convergence tolerance when cvg_mode = ENERGY
        type(external_obj) :: external                    ! external object
        type(fields_sc_obj) :: fields                     ! fields object
        type(electrons_sc_obj) :: electrons               ! electrons object
        logical :: have_forces                            ! indicates whether or not the forces have been calculated
        logical :: have_pressure                          ! indicates whether or not the pressure has been calculated
        logical :: have_stress_tensor                     ! indicates whether or not the stress tensor has been calculated
        logical :: null_residual_force                    ! indicates whether or not to null the residual force
        real(double), dimension(3) :: residual_force      ! residual force
        real(double), dimension(:,:), pointer :: forces   ! forces in the Cartesian representation
        real(double) :: pressure                          ! pressure
        real(double), dimension(3,3) :: stress_tensor     ! stress tensor
      end type

      type, public :: config_sc_obj
        private
        integer :: ref
        type(config_sc_rep), pointer :: o
      end type

!doc$
      public :: config_sc
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_cvg_mode
      public :: x_max_steps
      public :: x_external
      public :: x_fields
      public :: x_electrons
      public :: x_forces
      public :: x_pressure
      public :: x_stress_tensor
      public :: x_cell_energy
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: decompose
      public :: diary
      public :: diary_energy
      public :: diary_forces
      public :: diary_pressure
      public :: diary_stress_tensor
      public :: write_els_potential
      public :: write_restart

!cod$
      interface config_sc
         module procedure constructor_cfg
      end interface
      interface update
         module procedure update_cfg
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
      interface x_cvg_mode
         module procedure cfg_cvg_mode
      end interface
      interface x_max_steps
         module procedure cfg_max_steps
      end interface
      interface x_external
         module procedure cfg_external
      end interface
      interface x_fields
         module procedure cfg_fields
      end interface
      interface x_electrons
         module procedure cfg_electrons
      end interface
      interface x_forces
         module procedure cfg_forces
      end interface
      interface x_pressure
         module procedure cfg_pressure
      end interface
      interface x_stress_tensor
         module procedure cfg_stress_tensor
      end interface
      interface x_cell_energy
         module procedure cfg_cell_energy
      end interface
      interface forces
         module procedure forces_cfg
      end interface
      interface pressure
         module procedure pressure_cfg
      end interface
      interface stress_tensor
         module procedure stress_tensor_cfg
      end interface
      interface decompose
         module procedure decompose_cfg
      end interface
      interface diary
         module procedure diary_cfg, diary_atom_step
      end interface
      interface diary_energy
         module procedure diary_energy_sc
      end interface
      interface diary_forces
         module procedure diary_forces_sc
      end interface
      interface diary_pressure
         module procedure diary_pressure_sc
      end interface
      interface diary_stress_tensor
         module procedure diary_stress_tensor_sc
      end interface
      interface write_els_potential
         module procedure write_els_potential_cfg
      end interface
      interface write_restart
         module procedure write_restart_cfg
      end interface
      
      contains

! public routines

      function constructor_cfg() result(cfg)
!doc$ function config_sc() result(cfg)
        type(config_sc_obj) :: cfg
!       effects: Constructs a new cfg.
!       errors: Passes errors.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: mode, tag
        integer :: r_nsg
        integer(long) :: dsize, iosl, ndata, s4
        real(double) :: version
        type(tagio_obj) :: restf

        call start_timer("config_sc: constructor")

        call sync_configuration_errors()
        if (error("  Error on entry")) then
          cfg%ref = 0
          allocate( cfg%o )
          cfg%o%ref = 0
          goto 999
        end if

        cfg%ref = 0
        allocate( cfg%o )
        cfg%o%ref = 0
        cfg%o%g = x_ghost()

        ! read electronic convergence criteria
        call arglc("config_convergence",tag,found)
        if (.not.found) tag = "density"
        select case (trim(tag))
        case ("none")
          cfg%o%cvg_mode = NONE
        case ("energy")
          cfg%o%cvg_mode = ENERGY
          call arg("energy_tolerance",cfg%o%energy_tol,found)
          if (.not.found) cfg%o%energy_tol = 1.0e-6_double
          if (error(cfg%o%energy_tol < 0.0_double,"ERROR: energy_tolerance < 0")) goto 300
        case ("density")
          cfg%o%cvg_mode = DENSITY
        case ("wavefunctions")
          cfg%o%cvg_mode = WAVEFUNCTIONS
        case default
          if (error(.true.,"ERROR: config_convergence was not recognized")) goto 300
        end select
        call arg("config_steps",cfg%o%max_steps,found)
        if (.not.found) cfg%o%max_steps = 40
        if (error(cfg%o%max_steps < 0,"ERROR: config_steps < 0")) goto 300

        ! read whether or not to null the residual force
        call arglc("null_residual_force",tag,found)
        if (.not.found) tag = "on"
        select case (trim(tag))
        case ("on")
          cfg%o%null_residual_force = .true.
        case ("off")
          cfg%o%null_residual_force = .false.
        case default
          if (error(.true.,"ERROR: null_residual_force tag was not recognized")) goto 300
        end select

        ! read the restart mode
        call arglc("restart",mode,found)
        if (.not.found) mode = "off"

        select case (trim(mode))
        case ("off")

          call my(external(),cfg%o%external)              ; if (error()) goto 200
          call my(fields_sc(cfg%o%external),cfg%o%fields) ; if (error()) goto 200
          call my(electrons_sc(cfg%o%external,x_potential(cfg%o%fields)),cfg%o%electrons)

        case ("e","f","ef","fe","efe")

          ! open the restart file and check for the correct type (mkey)
          call my(tagio(trim(restart_path),TAGIO_READ,mkey,len(mkey)),restf)
          if (i_access(restf)) iosl = x_tagfd(restf)
          if (i_comm(restf)) call broadcast(FILE_SCOPE,iosl)
          if (error(iosl == 0,"ERROR: restart file was not found")) goto 300

          ! check the version number
          if (i_access(restf)) tios = findfirsttag(restf,"VERSION")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: VERSION tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
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
          case ("e")
            call my(external(restf),cfg%o%external)         ; if (error()) goto 100
            call my(fields_sc(cfg%o%external),cfg%o%fields) ; if (error()) goto 100
            call my(electrons_sc(cfg%o%external,x_potential(cfg%o%fields)),cfg%o%electrons)
          case ("f")
            call my(external(),cfg%o%external)                    ; if (error()) goto 100
            call my(fields_sc(cfg%o%external,restf),cfg%o%fields) ; if (error()) goto 100
            call my(electrons_sc(cfg%o%external,x_potential(cfg%o%fields)),cfg%o%electrons)
          case ("ef")
            call my(external(restf),cfg%o%external)               ; if (error()) goto 100
            call my(fields_sc(cfg%o%external,restf),cfg%o%fields) ; if (error()) goto 100
            call my(electrons_sc(cfg%o%external,x_potential(cfg%o%fields)),cfg%o%electrons)
          case ("fe")
            call my(external(),cfg%o%external)                    ; if (error()) goto 100
            call my(fields_sc(cfg%o%external,restf),cfg%o%fields) ; if (error()) goto 100
            call my(electrons_sc(cfg%o%external,x_potential(cfg%o%fields),restf),cfg%o%electrons)
          case ("efe")
            call my(external(restf),cfg%o%external)               ; if (error()) goto 100
            call my(fields_sc(cfg%o%external,restf),cfg%o%fields) ; if (error()) goto 100
            call my(electrons_sc(cfg%o%external,x_potential(cfg%o%fields),restf),cfg%o%electrons)
          end select

100       call glean(thy(restf)) ; if (error()) goto 200

        case default

          if (error(.true.,"ERROR: restart tag was not recognized")) goto 300

        end select

        call sync_configuration_errors() ; if (error()) goto 300
        call diary_construction_i(cfg%o)

        ! iterate to convergence
        call iterator_i(cfg%o)
        call sync_configuration_errors() ; if (error()) goto 300

        ! initialize the forces, pressure, and stress tensor
        allocate( cfg%o%forces(3,x_n_atoms(x_atoms(x_crystal(cfg%o%external)))) )
        cfg%o%have_forces = .false.
        cfg%o%have_pressure = .false.
        cfg%o%have_stress_tensor = .false.

200     call sync_configuration_errors()

300     continue

999     if (error("Exit config_sc_mod::constructor_cfg")) continue

        if (.not.error()) call stop_timer("config_sc: constructor")

      end function

      subroutine update_cfg(cfg,ext)
!doc$ subroutine update(cfg,ext)
        type(config_sc_obj) :: cfg
        type(external_obj), optional :: ext
!       modifies: cfg
!       effects: Updates cfg.
!       errors: Change in the layout. Passes errors.

!cod$
        logical :: ext_change, layout_change

        call start_timer("config_sc: update")

        call my(cfg)
        if (present(ext)) call my(ext)

        ! determine what has changed
        if (present(ext)) then
          ext_change = ( x_ghost(cfg%o%external) /= x_ghost(ext) )
          if (ext_change) then
            layout_change = ( x_ghost(x_layout(ext)) /= x_ghost(x_layout(cfg%o%external)) )
!           N.M. The changes that I have made below hopefully enable changes in layout.
!                However, they restart the calculation of the electronic structure
!            if (error(layout_change,"ERROR: layout changes are not currently allowed")) goto 100
          end if
        else
          ext_change = .false.
        end if

        ! update cfg
        if (ext_change) then
           call warn("ext_change")
          call own_i(cfg)
          cfg%o%g = x_ghost()
          cfg%o%external = ext
          if (layout_change) then
             cfg%o%fields = fields_sc(cfg%o%external) ; if (error()) goto 100
             call warn("layout change fields")
             cfg%o%electrons = electrons_sc(cfg%o%external,x_potential(cfg%o%fields)) ; if (error()) goto 100
          else
             call warn("in else branch")
             call update(cfg%o%fields,cfg%o%external,cfg%o%electrons) ; if (error()) goto 100
             call warn("after electrons update")
             call update(cfg%o%electrons,cfg%o%external,x_potential(cfg%o%fields)) ; if (error()) goto 100
             call warn("after fields update")
          end if
!          call update(cfg%o%fields,cfg%o%external,cfg%o%electrons) ; if (error()) goto 100
!          call update(cfg%o%electrons,cfg%o%external,x_potential(cfg%o%fields)) ; if (error()) goto 100
          call iterator_i(cfg%o) ; if (error()) goto 100
          cfg%o%have_forces = .false.
          cfg%o%have_pressure = .false.
          cfg%o%have_stress_tensor = .false.
        end if

100     call sync_configuration_errors()

        call glean(thy(cfg))
        if (present(ext)) call glean(thy(ext))

        if (error("Exit config_sc_mod::update_cfg")) continue

        if (.not.error()) call stop_timer("config_sc: update")

      end subroutine

      subroutine my_cfg(cfg)
!doc$ subroutine my(cfg)
        type(config_sc_obj) :: cfg

!cod$
        cfg%ref = cfg%ref + 1
        cfg%o%ref = cfg%o%ref + 1
      end subroutine

      subroutine my_new_cfg(cfgi,cfg)
!doc$ subroutine my(cfgi,cfg)
        type(config_sc_obj) :: cfgi, cfg

!cod$
        cfg%ref = 1
        cfg%o => cfgi%o
        cfg%o%ref = cfg%o%ref + 1
      end subroutine

      function thy_cfg(cfg) result(cfgo)
!doc$ function thy(cfg) result(cfgo)
        type(config_sc_obj) :: cfg, cfgo

!cod$
        cfg%ref = cfg%ref - 1
        cfg%o%ref = cfg%o%ref - 1
        cfgo%ref = cfg%ref
        cfgo%o => cfg%o
      end function

      subroutine glean_cfg(cfg)
!doc$ subroutine glean(cfg)
        type(config_sc_obj) :: cfg

!cod$
        if (cfg%o%ref < 1) then
          call glean(thy(cfg%o%external))
          call glean(thy(cfg%o%fields))
          call glean(thy(cfg%o%electrons))
          if (associated( cfg%o%forces )) deallocate( cfg%o%forces )
          deallocate( cfg%o )
        end if
      end subroutine

      subroutine bequeath_cfg(cfg)
!doc$ subroutine bequeath(cfg)
        type(config_sc_obj) :: cfg

!cod$
        continue
      end subroutine

      subroutine assign_cfg(cfg,cfg2)
!doc$ subroutine assignment(=)(cfg,cfg2)
        type(config_sc_obj), intent(inout) :: cfg
        type(config_sc_obj), intent(in) :: cfg2

!cod$
        type(config_sc_obj) :: cfgt
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
        type(config_sc_obj) :: cfg
        integer, dimension(2) :: r
!       effects: Returns cfg%ref and cfg%o%ref.

!cod$
        r(1) = cfg%ref
        r(2) = cfg%o%ref
        call glean(cfg)
      end function

      function cfg_ghost(cfg) result(g)
!doc$ function x_ghost(cfg) result(g)
        type(config_sc_obj) :: cfg
        type(ghost) :: g

!cod$
        call my(cfg)
        g = cfg%o%g
        call glean(thy(cfg))
      end function

      function cfg_cvg_mode(cfg) result(m)
!doc$ function x_cvg_mode(cfg) result(m)
        type(config_sc_obj) :: cfg
        integer :: m

!cod$
        call my(cfg)
        m = cfg%o%cvg_mode
        call glean(thy(cfg))
      end function

      function cfg_max_steps(cfg) result(s)
!doc$ function x_max_steps(cfg) result(s)
        type(config_sc_obj) :: cfg
        integer :: s

!cod$
        call my(cfg)
        s = cfg%o%max_steps
        call glean(thy(cfg))
      end function

      function cfg_external(cfg) result(ext)
!doc$ function x_external(cfg) result(ext)
        type(config_sc_obj) :: cfg
        type(external_obj) :: ext
!       effects: Returns the external component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%external,ext)
        call glean(thy(cfg))
        call bequeath(thy(ext))
      end function

      function cfg_fields(cfg) result(fd)
!doc$ function x_fields(cfg) result(fd)
        type(config_sc_obj) :: cfg
        type(fields_sc_obj) :: fd
!       effects: Returns the fields component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%fields,fd)
        call glean(thy(cfg))
        call bequeath(thy(fd))
      end function

      function cfg_electrons(cfg) result(el)
!doc$ function x_electrons(cfg) result(el)
        type(config_sc_obj) :: cfg
        type(electrons_sc_obj) :: el
!       effects: Returns the electrons component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%electrons,el)
        call glean(thy(cfg))
        call bequeath(thy(el))
      end function

      function cfg_forces(cfg) result(f)
!doc$ function x_forces(cfg) result(f)
        type(config_sc_obj) :: cfg
        real(double), dimension(size(cfg%o%forces,1),size(cfg%o%forces,2)) :: f
!       modifies: cfg if cfg%o%have_forces is .false. (Note: this is an optimization).
!       requires: cfg%o%forces be allocated.
!       effects: Returns the forces in atomic units.
!       errors: Passes errors.

!cod$
        call my(cfg)
        if (.not.cfg%o%have_forces) then
          call forces_i(cfg%o) ; if (error()) goto 100
        end if
        f = cfg%o%forces
100     call glean(thy(cfg))
        if (error("Exit config_sc_mod::cfg_forces")) continue
      end function

      function cfg_pressure(cfg) result(p)
!doc$ function x_pressure(cfg) result(pressure)
        type(config_sc_obj) :: cfg
        real(double) :: p
!       modifies: cfg if cfg%o%pressure is .false. (Note: this is an optimization).
!       effects: Returns the pressure in atomic units.
!       errors: Passes errors.

!cod$
        call my(cfg)
        if (.not.cfg%o%have_pressure) then
          call pressure_i(cfg%o) ; if (error()) goto 100
        end if
        p = cfg%o%pressure
100     call glean(thy(cfg))
        if (error("Exit config_sc_mod::cfg_pressure")) continue
      end function

      function cfg_stress_tensor(cfg) result(s)
!doc$ function x_stress_tensor(cfg) result(s)
        type(config_sc_obj) :: cfg
        real(double), dimension(3,3) :: s
!       modifies: cfg if cfg%o%stress_tensor is .false. (Note: this is an optimization).
!       effects: Returns the stress tensor in atomic units.
!       errors: Passes errors.

!cod$
        call my(cfg)
        if (.not.cfg%o%have_stress_tensor) then
          call stress_tensor_i(cfg%o) ; if (error()) goto 100
        end if
        s = cfg%o%stress_tensor
100     call glean(thy(cfg))
        if (error("Exit config_sc_mod::cfg_stress_tensor")) continue
      end function
           
      function cfg_cell_energy(cfg) result(ce)
!doc$ function x_cell_energy(cfg) result(ce)
        type(config_sc_obj) :: cfg
        real(double) :: ce
!       effects:  Returns the cell energy.

!cod$
        call my(cfg)
        ce = x_energy(cfg%o%electrons) + x_energy(cfg%o%fields)
        call glean(thy(cfg))
      end function

      subroutine forces_cfg(cfg)
!doc$ subroutine forces(cfg)
        type(config_sc_obj) :: cfg
!       modifies: cfg (iff cfg%o%have_forces = .false.)
!       effects: Computes cfg%o%forces (iff cfg%o%have_forces = .false.).
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag

        call my(cfg)

        call arglc("forces",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          call forces_i(cfg%o)
        case ("off",".false.")
          continue
        case default
          if (error(.true.,"ERROR: forces tag was not recognized")) goto 100
        end select

100     call glean(thy(cfg))

        if (error("Exit config_sc_mod::forces")) continue

      end subroutine

      subroutine pressure_cfg(cfg)
!doc$ subroutine pressure(cfg)
        type(config_sc_obj) :: cfg
!       modifies: cfg (iff cfg%o%have_pressure = .false.)
!       effects: Computes cfg%o%pressure (iff cfg%o%have_pressure = .false.).
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag

        call my(cfg)

        call arglc("pressure",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          call pressure_i(cfg%o)
        case ("off",".false.")
          continue
        case default
          if (error(.true.,"ERROR: pressure tag was not recognized")) goto 100
        end select

100     call glean(thy(cfg))

        if (error("Exit config_sc_mod::pressure")) continue

      end subroutine

      subroutine stress_tensor_cfg(cfg)
!doc$ subroutine stress_tensor(cfg)
        type(config_sc_obj) :: cfg
!       modifies: cfg (iff cfg%o%have_stress_tensor = .false.)
!       effects: Computes cfg%o%stress_tensor (iff cfg%o%have_stress_tensor = .false.).
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag
        real(double), dimension(:,:), allocatable :: s, st

        call my(cfg)

        call arglc("stress_tensor",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          call stress_tensor_i(cfg%o)
        case ("off",".false.")
          continue
        case default
          if (error(.true.,"ERROR: stress_tensor tag was not recognized")) goto 100
        end select

100     if (allocated( s )) deallocate( s )
        if (allocated( st )) deallocate( st )

        call glean(thy(cfg))

        if (error("Exit config_sc_mod::stress_tensor")) continue

      end subroutine

      subroutine decompose_cfg(cfg)
!doc$ subroutine decompose(cfg)
        type(config_sc_obj) :: cfg
!       effects: Decomposes the Kohn-Sham functions into s, p, & d spherical harmonics around user-defined and atom sites.
!       errors: decomposition tag not recognized. Passes errors.

!cod$
        logical :: exist_file, found
        character(tag_sz) :: type
        character(line_len) :: mode, switch
        character(13+tag_sz) :: dr_tag
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

        call start_timer("config_sc: decompose")

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
            call warn("dcomp_radius tag was not found - using the value 2.5")
            site_data(4,is) = 2.5_double
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

600     if (.not.error()) call stop_timer("config_sc: decompose")

700     call glean(thy(cfg))

        if (error("Exit config_sc_mod::decompose_cfg")) continue

      end subroutine

      subroutine diary_cfg(cfg)
!doc$ subroutine diary(cfg)
        type(config_sc_obj) :: cfg
!       modifies: Output stream
!       effects: Prints cfg information.
!       errors: Passes errors.

!cod$
        call my(cfg)
        call diary(cfg%o%electrons) ; if (error()) goto 100
        call diary(cfg%o%fields) ; if (error()) goto 100
        call diary_energy(cfg) 

100     call glean(thy(cfg))
        if (error("Exit config_sc_mod::diary_cfg")) continue
      end subroutine

      subroutine diary_atom_step(cfg,step)
!doc$ subroutine diary(cfg,step)
        type(config_sc_obj) :: cfg
        integer, intent(inout) :: step
!       modifies: Output stream and cfg if cfg%o%have_forces is .false.
!       effects: Prints cfg information relevant to a particular step.
!       errors: Passes errors.

!cod$
        call my(cfg)
        call diary(cfg%o%external,step)     ; if (error()) goto 100
        call diary_forces(cfg)              ; if (error()) goto 100
        call diary_energy(cfg,all=.false.)  ; if (error()) goto 100
100     call glean(thy(cfg))
        if (error("Exit config_sc_mod::diary_atom_step")) continue
      end subroutine

      subroutine diary_energy_sc(cfg,all)
!doc$ subroutine diary_energy(cfg,all)
        type(config_sc_obj) :: cfg
        logical, intent(in), optional :: all
!       modifies: Output stream
!       effects: Prints energies.

!cod$
        logical :: print_all
        integer :: md
        real(double) :: cell_energy

        if (present(all)) then
          print_all = all
        else
          print_all = .true.
        end if

        call my(cfg)
 
        cell_energy = cell_energy_i(cfg%o)
        if (print_all) then
          md = 12 + log10(max(abs(max_energy(cfg%o%electrons)),abs(max_energy(cfg%o%fields)),abs(cell_energy)))
        end if

        if (i_access( diaryfile() )) then
          if (print_all) then
            write(x_unit(diaryfile()),'(/,t4,"Energy components (Ryd):")')
            call diary_energies(cfg%o%electrons,md)
            call diary_energies(cfg%o%fields,md)
            write(x_unit(diaryfile()),'(t6,"---------------------------------------")')
            select case (md)
            case (12)
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f12.9)') cell_energy
            case (13)
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f13.9)') cell_energy
            case (14)
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f14.9)') cell_energy
            case (15)
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f15.9)') cell_energy
            case (16)
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f16.9)') cell_energy
            case (17)
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f17.9)') cell_energy
            case default
              write(x_unit(diaryfile()),'(t6,"cell energy           = ",f18.9)') cell_energy
            end select
          else
            write(x_unit(diaryfile()),'(/,t4,"Cell energy = ",f0.9," Ryd")') cell_energy
          end if
        end if

        call glean(thy(cfg))

      end subroutine

      subroutine diary_forces_sc(cfg)
!doc$ subroutine diary_forces(cfg)
        type(config_sc_obj) :: cfg
!       modifies: Output stream (iff cfg%o%have_forces = .true.).
!       effects: Writes the forces to the diary file (iff cfg%o%have_forces = .true.).
!       errors: Passes errors.

!cod$
        integer :: ia, na

        call my(cfg)

        if (cfg%o%have_forces) then
          na = size(cfg%o%forces,2)
          if (i_access( diaryfile() )) then
            write(x_unit(diaryfile()),'(/,t4,"Atomic forces:")')
            write(x_unit(diaryfile()),'(20x,"atom",11x,"Fx",12x,"Fy",12x,"Fz")')
            write(x_unit(diaryfile()),'(19x,"---------------------------------------------------")')
            do ia = 1,na
              if (ia == 1) then
                write(x_unit(diaryfile()),'(19x,i4,3x,sp,3(4x,f10.6),3x,"Ryd/Bohr")') ia, cfg%o%forces(:,ia)
              else
                write(x_unit(diaryfile()),'(19x,i4,3x,sp,3(4x,f10.6))') ia, cfg%o%forces(:,ia)
              end if
            end do
            if (cfg%o%null_residual_force) then
              write(x_unit(diaryfile()),'(/,t4,"Residual force (removed) = ",sp,3(4x,f10.6))') cfg%o%residual_force
            else
              write(x_unit(diaryfile()),'(/,t4,"Residual force = ",sp,3(4x,f10.6))') cfg%o%residual_force
            end if
          end if
        end if

        call glean(thy(cfg))

        if (error("Exit config_sc_mod::diary_forces_sc")) continue

      end subroutine

      subroutine diary_pressure_sc(cfg)
!doc$ subroutine diary_pressure(cfg)
        type(config_sc_obj) :: cfg
!       modifies: Output stream (iff cfg%o%have_pressure = .true.)
!       effects: Writes the pressure to the diary file (iff cfg%o%have_pressure = .true.).
!       errors: Passes errors.

!cod$
        call my(cfg)

        if (cfg%o%have_pressure) then
          if (i_access( diaryfile() )) then
            write(x_unit(diaryfile()),'(/,t4,"Pressure: ",f14.4,3x,"kbar")') cfg%o%pressure*147105.164_double
          end if
        end if

        call glean(thy(cfg))

        if (error("Exit config_sc_mod::diary_pressure_sc")) continue

      end subroutine

      subroutine diary_stress_tensor_sc(cfg)
!doc$ subroutine diary_stress_tensor(cfg)
        type(config_sc_obj) :: cfg
!       modifies: Output stream (iff cfg%o%have_stress_tensor = .true.)
!       effects: Writes the stress tensor to the diary file (iff cfg%o%have_stress_tensor = .true.).
!       errors: Passes errors.

!cod$
        call my(cfg)

        if (cfg%o%have_stress_tensor) then
          if (i_access( diaryfile() )) then
            write(x_unit(diaryfile()),'(/,t4,"Stress Tensor:")')
            write(x_unit(diaryfile()),'(t8,3f12.4,3x,"kbar")') cfg%o%stress_tensor(1,:)*147105.164_double
            write(x_unit(diaryfile()),'(t8,3f12.4)')           cfg%o%stress_tensor(2,:)*147105.164_double
            write(x_unit(diaryfile()),'(t8,3f12.4)')           cfg%o%stress_tensor(3,:)*147105.164_double
          end if
        end if

        call glean(thy(cfg))

        if (error("Exit config_sc_mod::diary_stress_tensor_sc")) continue

      end subroutine

      subroutine write_els_potential_cfg(cfg)
!doc$ subroutine write_els_potential(cfg)
        type(config_sc_obj) :: cfg
!       modifies: Output stream
!       effects: Writes the sxdefectalign els_potential file.
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag
        type(lattice_obj) :: lat
        type(atoms_obj) :: atomss

        call my(cfg)

        call my(x_lattice(x_crystal(cfg%o%external)),lat)
        call my(x_atoms(x_crystal(cfg%o%external)),atomss)

        call arglc("write_els_potential",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          call write_els_potential(cfg%o%fields,lat,atomss) ; if (error()) goto 100
        case ("off",".false.")
          continue
        case default
          if (error(.true.,"ERROR: write_els_potential tag was not recognized")) goto 100
        end select

100     call glean(thy(lat))
        call glean(thy(atomss))

        call glean(thy(cfg))

        if (error("Exit config_sc_mod::write_els_potential_cfg")) continue

      end subroutine

      subroutine write_restart_cfg(cfg,mode_in)
!doc$ subroutine write_restart(cfg,mode_in)
        type(config_sc_obj) :: cfg
        character(line_len), intent(in), optional :: mode_in
!       effects: Writes restart information according to mode_in or write_restart tag in argvf.
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: mode
        integer(long) :: dsize, iosl, ndata, s4
        type(tagio_obj) :: nrestf

        call my(cfg)

        if (.not.present(mode_in)) then
          call arglc("write_restart",mode,found)
          if (.not.found) mode = "off"
        end if

        select case (trim(mode))
        case ("off")
          continue
        case default

          ! Open the restart file
          call my(tagio(trim(new_restart_path),TAGIO_WRITE,mkey,len(mkey)),nrestf)

          if (i_access(nrestf)) then

            ! Write the version number
            call writetag(nrestf,"VERSION")
            dsize = sizeof_double
            ndata = 1
            call writef(es_version,dsize,ndata,x_tagfd(nrestf),iosl)

            ! Write the number of spin groups
            call writetag(nrestf,"NUMBER_OF_SGROUPS")
            dsize = sizeof_long
            ndata = 1
            s4 = mpi_nsgroups()
            call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          end if

          ! Write the rest of the restart information
          select case (trim(mode))
          case ("e")
            call write_restart(cfg%o%external,nrestf)
          case ("ef","f")
            call write_restart(cfg%o%external,nrestf)
            call write_restart(cfg%o%fields,nrestf)
          case ("efe","fe")
            call write_restart(cfg%o%external,nrestf)
            call write_restart(cfg%o%fields,nrestf)
            call write_restart(cfg%o%electrons,nrestf)
          case default
            call warn("WARNING: write_restart tag is not recognized - writing a full restart file")
            call write_restart(cfg%o%external,nrestf)
            call write_restart(cfg%o%fields,nrestf)
            call write_restart(cfg%o%electrons,nrestf)
          end select

          ! Close the restart file
          call glean(thy(nrestf))

        end select

        call glean(thy(cfg))

        if (error("Exit config_sc_mod::write_restart_cfg")) continue

      end subroutine

! private routines

      subroutine iterator_i(cfgr)
        type(config_sc_rep) :: cfgr

        logical :: done
        integer :: lc

        lc = 0
        do
          if (error(user_abort(),"USER INITIATED ABORT")) goto 100
          lc = lc + 1
          !write(*,*) "Here 1",mpi_myproc(world)
          call update(cfgr%fields,cfgr%external,cfgr%electrons) ; if (error()) goto 100
          call check_convergence_i(cfgr,lc,done)
          if (done) exit
          !write(*,*) "Here 2",mpi_myproc(world)
          call update(cfgr%electrons,cfgr%external,x_potential(cfgr%fields)) ; if (error()) goto 100
        end do

100     if (error("Exit config_sc_mod::iterator_i")) continue

      end subroutine

      subroutine check_convergence_i(cfgr,lc,done)
        type(config_sc_rep) :: cfgr
        integer, intent(in) :: lc
        logical, intent(out) :: done
        real(double), save :: ce
        real(double) :: pce, rn
        done = .false.
        select case (cfgr%cvg_mode)
        case (NONE)
          ce = cell_energy_i(cfgr)
          rn = x_residual_norm(cfgr%fields)
          if (i_access(diaryfile())) then
            if (lc == 1) then
              write(x_unit(diaryfile()),'(/,t4,"Self-consistent step  ",i1,":  density residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            elseif (lc < 10) then
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step  ",i1,":  density residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            else
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step ", i0,":  density residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            end if
          end if
          if (i_access(output)) then
            if (mpi_nconfigs() == 1) then
              if (lc < 10) then
                write(x_unit(output),'("Self-consistent step  ",i1,":  density residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') lc, rn, ce
              else
                write(x_unit(output),'("Self-consistent step ", i0,":  density residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') lc, rn, ce
              end if
            else
              if (lc < 10) then
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step  ",i1, &
                                       & ":  density residual = ",es10.4,",  cell energy = ",f0.9)') &
                                       & mpi_myconfig(), lc, rn, ce
              else
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step  ",i0, &
                                       & ":  density residual = ",es10.4,",  cell energy = ",f0.9)') &
                                       & mpi_myconfig(), lc, rn, ce
              end if
            end if
          end if
        case (ENERGY)
          if (lc == 1) then
            pce = 0.0_double
          else
            pce = ce
          end if
          ce = cell_energy_i(cfgr)
          if (i_access(diaryfile())) then
            if (lc == 1) then
              write(x_unit(diaryfile()),'(/,t4,"Self-consistent step  ",i1,":  cell energy = ",f0.9)') lc, ce
            elseif (lc < 10) then
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step  ",i1,":  cell energy = ",f0.9, &
                                     & ",  energy change = ",es11.4 )') lc, ce, (ce - pce)
            else
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step ", i0,":  cell energy = ",f0.9, &
                                     & ",  energy change = ",es11.4 )') lc, ce, (ce - pce)
            end if
          end if
          if (i_access(output)) then
            if (mpi_nconfigs() == 1) then
              if (lc == 1) then
                write(x_unit(output),'("Self-consistent step  ",i1,":  cell energy = ",f0.9)') lc, ce
              elseif (lc < 10) then
                write(x_unit(output),'("Self-consistent step  ",i1,":  cell energy = ",f0.9, &
                                       & ",  energy change = ",es11.4 )') lc, ce, (ce - pce)
              else
                write(x_unit(output),'("Self-consistent step ", i0,":  cell energy = ",f0.9, &
                                       & ",  energy change = ",es11.4 )') lc, ce, (ce - pce)
              end if
            else
              if (lc == 1) then
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step  ",i1, &
                                       & ":  cell energy = ",f0.9)') mpi_myconfig(), lc, ce
              elseif (lc < 10) then
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step  ",i1, &
                                       & ":  cell energy = ",f0.9, &
                                       & ",  energy change = ",es11.4 )') mpi_myconfig(), lc, ce, (ce - pce)
              else
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step ", i0, &
                                       & ":  cell energy = ",f0.9, &
                                       & ",  energy change = ",es11.4 )') mpi_myconfig(), lc, ce, (ce - pce)
              end if
            end if
          end if
          if ( (lc > 1) .and. (abs(ce - pce) < cfgr%energy_tol) ) done = .true.
        case (DENSITY)
          ce = cell_energy_i(cfgr)
          rn = x_residual_norm(cfgr%fields)
          if (i_access(diaryfile())) then
            if (lc == 1) then
              write(x_unit(diaryfile()),'(/,t4,"Self-consistent step  ",i1,":  density residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            elseif (lc < 10) then
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step  ",i1,":  density residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            else
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step ", i0,":  density residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            end if
          end if
          if (i_access(output)) then
            if (mpi_nconfigs() == 1) then
              if (lc < 10) then
                write(x_unit(output),'("Self-consistent step  ",i1,":  density residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') lc, rn, ce
              else
                write(x_unit(output),'("Self-consistent step ", i0,":  density residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') lc, rn, ce
              end if
            else
              if (lc < 10) then
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step  ",i1, &
                                       & ":  density residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') mpi_myconfig(), lc, rn, ce
              else
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step ", i0, &
                                       & ":  density residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') mpi_myconfig(), lc, rn, ce
              end if
            end if
          end if
          done = x_converged(cfgr%fields)
        case (WAVEFUNCTIONS)
          ce = cell_energy_i(cfgr)
          rn = x_residual_norm(cfgr%electrons)
          if (i_access(diaryfile())) then
            if (lc == 1) then
              write(x_unit(diaryfile()),'(/,t4,"Self-consistent step  ",i1,":  wavefunctions residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            elseif (lc < 10) then
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step  ",i1,":  wavefunctions residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            else
              write(x_unit(diaryfile()),'(  t4,"Self-consistent step ", i0,":  wavefunctions residual = ",es10.4, &
                                     & ",  cell energy = ",f0.9)') lc, rn, ce
            end if
          end if
          if (i_access(output)) then
            if (mpi_nconfigs() == 1) then
              if (lc < 10) then
                write(x_unit(output),'("Self-consistent step  ",i1,":  wavefunctions residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') lc, rn, ce
              else
                write(x_unit(output),'("Self-consistent step ", i0,":  wavefunctions residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') lc, rn, ce
              end if
            else
              if (lc < 10) then
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step  ",i1, &
                                       & ":  wavefunctions residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') mpi_myconfig(), lc, rn, ce
              else
                write(x_unit(output),'("Configuration ",i2.2,";  Self-consistent step ", i0, &
                                       & ":  wavefunctions residual = ",es10.4, &
                                       & ",  cell energy = ",f0.9)') mpi_myconfig(), lc, rn, ce
              end if
            end if
          end if
          done = x_converged(cfgr%electrons)
        end select
        call flushbuf(diaryfile())
        if (lc >= cfgr%max_steps) done = .true.
      end subroutine

      subroutine diary_construction_i(cfgr)
        type(config_sc_rep) :: cfgr

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,"Config object construction:")')
          write(x_unit(diaryfile()),'(/,t4,"Self-consistent calculation")')
          select case (cfgr%cvg_mode)
          case (NONE)
            write(x_unit(diaryfile()),'(/,t4,"Convergence will not be checked")')
          case (ENERGY)
            write(x_unit(diaryfile()),'(/,t4,"Convergence will be determined from the total energy")')
          case (DENSITY)
            write(x_unit(diaryfile()),'(/,t4,"Convergence will be determined from the density residual")')
          case (WAVEFUNCTIONS)
            write(x_unit(diaryfile()),'(/,t4,"Convergence will be determined from the wavefunctions residual")')
          end select
        end if

      end subroutine

      subroutine own_i(cfg)
        type(config_sc_obj) :: cfg
        type(config_sc_obj) :: cfgt
        if (cfg%ref < cfg%o%ref) then
          allocate( cfgt%o )
          cfgt%o%ref = 0
          cfgt%o%g = cfg%o%g
          cfgt%o%cvg_mode = cfg%o%cvg_mode
          cfgt%o%max_steps = cfg%o%max_steps
          cfgt%o%energy_tol = cfg%o%energy_tol
          call my(cfg%o%external,cfgt%o%external)
          call my(cfg%o%fields,cfgt%o%fields)
          call my(cfg%o%electrons,cfgt%o%electrons)
          cfgt%o%null_residual_force = cfg%o%null_residual_force
          allocate( cfgt%o%forces(size(cfg%o%forces,1),size(cfg%o%forces,2)) )
          cfgt%o%have_forces = cfg%o%have_forces
          if (cfg%o%have_forces) then
            cfgt%o%forces = cfg%o%forces
            cfgt%o%residual_force = cfg%o%residual_force
          end if
          cfgt%o%have_pressure = cfg%o%have_pressure
          if (cfg%o%have_pressure) cfgt%o%pressure = cfg%o%pressure
          cfgt%o%have_stress_tensor = cfg%o%have_stress_tensor
          if (cfg%o%have_stress_tensor) cfgt%o%stress_tensor = cfg%o%stress_tensor
          cfg%o%ref = cfg%o%ref - cfg%ref
          cfg%o => cfgt%o
          cfg%o%ref = cfg%o%ref + cfg%ref
        end if
      end subroutine

      function cell_energy_i(cfgr) result(ce)
        type(config_sc_rep) :: cfgr
        real(double) :: ce
        ce = x_energy(cfgr%electrons) + x_energy(cfgr%fields)
      end function

      subroutine forces_i(cfgr)
        type(config_sc_rep) :: cfgr
!       requires: cfgr%forces be allocated.
        integer :: ia, na
        real(double), dimension(:,:), allocatable :: ft
        call start_timer("config_sc: forces_i")
        na = size(cfgr%forces,2)
        allocate( ft(3,na) )
        call forces(cfgr%fields,ft) ; if (error()) goto 100
        cfgr%forces = ft
        call forces(cfgr%electrons,ft) ; if (error()) goto 100
        cfgr%forces = cfgr%forces + ft
        do ia = 1,na
          ft(:,ia) = r2lat(x_lattice(x_crystal(cfgr%external)),cfgr%forces(:,ia))
        end do
        call symmetrize_vectors(cfgr%external,ft)
        do ia = 1,na
          cfgr%forces(:,ia) = lat2r(x_lattice(x_crystal(cfgr%external)),ft(:,ia))
        end do
        cfgr%residual_force = 0.0_double
        do ia = 1,na
          cfgr%residual_force = cfgr%residual_force + cfgr%forces(:,ia)
        end do
        if (cfgr%null_residual_force) then
          do ia = 1,na
            cfgr%forces(:,ia) = cfgr%forces(:,ia) - cfgr%residual_force/real(na,double)
          end do
        end if
        cfgr%have_forces = .true.
100     if (allocated( ft )) deallocate( ft )
        if (error("Exit config_sc_mod::forces_i")) continue
        call stop_timer("config_sc: forces_i")
      end subroutine

      subroutine pressure_i(cfgr)
        type(config_sc_rep) :: cfgr
        real(double) :: p, pt
        call start_timer("config_sc: pressure_i")
        p = 0.0_double
        call pressure(cfgr%fields,pt) ; if (error()) goto 100
        p = p + pt
        call pressure(cfgr%electrons,pt) ; if (error()) goto 100
        p = p + pt
        cfgr%pressure = p
        cfgr%have_pressure = .true.
100     if (error("Exit config_sc_mod::pressure_i")) continue
        if (.not.error()) call stop_timer("config_sc: pressure_i")
      end subroutine

      subroutine stress_tensor_i(cfgr)
        type(config_sc_rep) :: cfgr
        real(double), dimension(:,:), allocatable :: s, st
        call start_timer("config_sc: stress_tensor_i")
        allocate( s(3,3), st(3,3) )
        s = 0.0_double
        call stress_tensor(cfgr%fields,st) ; if (error()) goto 100
        s = s + st
        call stress_tensor(cfgr%electrons,st) ; if (error()) goto 100
        s = s + st
        st = r2lat_tensor(x_lattice(x_crystal(cfgr%external)),s)
        call symmetrize_tensor(cfgr%external,st) ; if (error()) goto 100
        s = lat2r_tensor(x_lattice(x_crystal(cfgr%external)),st)
        cfgr%stress_tensor = s
        cfgr%have_stress_tensor = .true.
100     if (allocated( s )) deallocate( s )
        if (allocated( st )) deallocate( st )
        if (error("Exit config_sc_mod::stress_tensor_i")) continue
        if (.not.error()) call stop_timer("config_sc: stress_tensor_i")
      end subroutine

      end module

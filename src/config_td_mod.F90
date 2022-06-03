!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module config_td_mod
!doc$ module config_td_mod

!     One datatype is available here: type(config_td_obj)

!     config_td_mod encapsulates a time-dependent electronic structure solution for a configuration of atoms.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use arg_mod
      use tagio_mod
      use ghost_mod
      use diary_mod
      use electrons_td_mod
      use electrons_sc_mod
      use crystal_mod
      use layout_mod
      use lattice_mod
      use atoms_mod
      use external_mod
      use fields_td_mod
      use config_sc_mod
      use fields_sc_mod
      use gen_potential_mod
      use grid_mod
      use interrupt_mod
      use timing_mod

!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = CONFIG
      integer, parameter :: SIMPLE = 0
      integer, parameter :: SELFCONSISTENT = 1

      real(double), parameter :: default_tol = 1.d-10

      type :: config_td_rep
        integer :: ref
        type(ghost) :: g
        type(external_obj) :: external                    ! external object
        type(fields_td_obj) :: fields                     ! fields object
        type(electrons_td_obj) :: electrons               ! electrons object
        type(config_sc_obj) :: cfg_sc                     ! a self-consistent config for the same external
        type(electrons_sc_obj) :: adiabatic_el        ! an electrons constructed with the current TDDFT Hamiltonian
        type(grid_obj) :: guess_pot                       ! extrapolated local potential
        integer :: propagation_type                       ! simple, self-consistent, etc...
        integer :: num_sc_steps                           ! number of self-conistent steps
        integer :: max_sc_steps                           ! max number of sc steps
        real(double) :: propagation_tol                   ! tolerance for the self-consistent prop
        logical :: compute_forces                         ! indicates whether or not to compute the forces
        logical :: compute_pressure                       ! indicates whether or not to compute the pressure
        logical :: compute_stress_tensor                  ! indicates whether or not to compute the stress tensor
        logical :: have_forces                            ! indicates whether or not the forces have been calculated
        logical :: have_pressure                          ! indicates whether or not the pressure has been calculated
        logical :: have_stress_tensor                     ! indicates whether or not the stress tensor has been calculated
        logical :: null_residual_force                    ! indicates whether or not to null the residual force
        real(double), dimension(:,:), pointer :: forces   ! forces in the Cartesian representation
        real(double) :: pressure                          ! pressure
        real(double), dimension(3,3) :: stress_tensor     ! stress tensor
      end type

      type, public :: config_td_obj
        private
        integer :: ref
        type(config_td_rep), pointer :: o
      end type

!doc$
      public :: config_td
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_external
      public :: x_fields
      public :: x_electrons
      public :: x_config_sc
      public :: x_adiabatic_states
      public :: x_forces
      public :: x_pressure
      public :: x_stress_tensor
      public :: x_cell_energy
      public :: x_num_sc_steps
      public :: diary
      public :: diary_energy
      public :: diary_forces
      public :: diary_pressure
      public :: diary_stress_tensor
      public :: decompose
      public :: write_restart

!cod$
      interface config_td
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
      interface x_external
         module procedure cfg_external
      end interface
      interface x_fields
         module procedure cfg_fields
      end interface
      interface x_electrons
         module procedure cfg_electrons
      end interface
      interface x_config_sc
         module procedure cfg_config_sc
      end interface
      interface x_adiabatic_states
         module procedure cfg_adiabatic_states
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
      interface x_num_sc_steps
         module procedure cfg_num_sc_steps
      end interface

      interface diary
         module procedure diary_cfg, diary_atom_step
      end interface
      interface diary_energy
         module procedure diary_energy_td
      end interface
      interface diary_forces
         module procedure diary_forces_td
      end interface
      interface diary_pressure
         module procedure diary_pressure_td
      end interface
      interface diary_stress_tensor
         module procedure diary_stress_tensor_td
      end interface
      interface decompose
         module procedure decompose_cfg
      end interface
      interface write_restart
         module procedure write_restart_cfg
      end interface
      
      contains

! public routines

      function constructor_cfg(restf) result(cfg)
!doc$ function config_td(restf) result(cfg)
        type(tagio_obj), optional :: restf
        type(config_td_obj) :: cfg
!       effects: Constructs a new cfg.
!       errors: Restart file not found. Restart tag not recognized. Passes errors.

!cod$
        logical :: found

        call start_timer("config_td: constructor")

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

        ! determine the initialization mode and branch accordingly
        if (.not.present(restf)) then
          call standard_portal_i(cfg%o)
        else
          call my(restf)
          call restart_portal_i(cfg%o,restf)
          call glean(thy(restf))
        end if

999     if (error("Exit config_td_mod::constructor_cfg")) continue

        if (.not.error()) call stop_timer("config_td: constructor")

      end function

      subroutine update_cfg(cfg,ext)
!doc$ subroutine update(cfg,ext)
        type(config_td_obj) :: cfg
        type(external_obj), optional :: ext
!       modifies: cfg
!       effects: Updates cfg.
!       errors: Change in the layout. Passes errors.

!cod$
        logical :: ext_change, layout_change

        call start_timer("config_td: update")

        call my(cfg)
        if (present(ext)) call my(ext)

        ! determine what has changed
        if (present(ext)) then
          ext_change = ( x_ghost(cfg%o%external) /= x_ghost(ext) )
          if (ext_change) then
             layout_change = ( x_ghost(x_layout(ext)) /= x_ghost(x_layout(cfg%o%external)) )
             if (error(layout_change,"ERROR: layout changes are not currently allowed")) goto 100
          end if
        else
          ext_change = .false.
        end if


        ! update cfg
        call own_i(cfg)
        cfg%o%g = x_ghost()
        
        !** Update external, if necessary
        if (ext_change) cfg%o%external = ext

        ! Propagate
        select case (cfg%o%propagation_type)
        case (SIMPLE)
           !** First propagate the electrons
           call update(cfg%o%electrons,cfg%o%external,x_potential(cfg%o%fields)) ; if (error()) goto 100

           !** Next update calculate the fields that are consistent with the propagated electrons.
           call update(cfg%o%fields,cfg%o%external,cfg%o%electrons) ; if (error()) goto 100

        case (SELFCONSISTENT)
           call take_sc_step_i(cfg%o)

        case default
           if (error(.true.,"ERROR: unrecognized propagation_type")) goto 100
        end select

        !** The update of cfg_sc is delayed until x_config_sc() is called

        !** The update of adiabatic_el is delayed until x_adiabatic_states() is called

        !** Forces, pressure, stress tensor
        deallocate( cfg%o%forces ) ; allocate( cfg%o%forces(3,x_n_atoms(x_atoms(x_crystal(cfg%o%external)))) )
        cfg%o%have_forces = .false.
        cfg%o%have_pressure = .false.
        cfg%o%have_stress_tensor = .false.

100     call sync_configuration_errors()

        call glean(thy(cfg))
        if (present(ext)) call glean(thy(ext))

        if (error("Exit config_td_mod::update_cfg")) continue

        if (.not.error()) call stop_timer("config_td: update")

      end subroutine

      subroutine my_cfg(cfg)
!doc$ subroutine my(cfg)
        type(config_td_obj) :: cfg

!cod$
        cfg%ref = cfg%ref + 1
        cfg%o%ref = cfg%o%ref + 1
      end subroutine

      subroutine my_new_cfg(cfgi,cfg)
!doc$ subroutine my(cfgi,cfg)
        type(config_td_obj) :: cfgi, cfg

!cod$
        cfg%ref = 1
        cfg%o => cfgi%o
        cfg%o%ref = cfg%o%ref + 1
      end subroutine

      function thy_cfg(cfg) result(cfgo)
!doc$ function thy(cfg) result(cfgo)
        type(config_td_obj) :: cfg, cfgo

!cod$
        cfg%ref = cfg%ref - 1
        cfg%o%ref = cfg%o%ref - 1
        cfgo%ref = cfg%ref
        cfgo%o => cfg%o
      end function

      subroutine glean_cfg(cfg)
!doc$ subroutine glean(cfg)
        type(config_td_obj) :: cfg

!cod$
        if (cfg%o%ref < 1) then
          call glean(thy(cfg%o%external))
          call glean(thy(cfg%o%fields))
          call glean(thy(cfg%o%electrons))
          call glean(thy(cfg%o%guess_pot))
          call glean(thy(cfg%o%cfg_sc))
          call glean(thy(cfg%o%adiabatic_el))
          if (associated( cfg%o%forces )) deallocate( cfg%o%forces )
          deallocate( cfg%o )
        end if
      end subroutine

      subroutine bequeath_cfg(cfg)
!doc$ subroutine bequeath(cfg)
        type(config_td_obj) :: cfg

!cod$
        continue
      end subroutine

      subroutine assign_cfg(cfg,cfg2)
!doc$ subroutine assignment(=)(cfg,cfg2)
        type(config_td_obj), intent(inout) :: cfg
        type(config_td_obj), intent(in) :: cfg2

!cod$
        type(config_td_obj) :: cfgt
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
        type(config_td_obj) :: cfg
        integer, dimension(2) :: r
!       effects: Returns cfg%ref and cfg%o%ref.

!cod$
        r(1) = cfg%ref
        r(2) = cfg%o%ref
        call glean(cfg)
      end function

      function cfg_ghost(cfg) result(g)
!doc$ function x_ghost(cfg) result(g)
        type(config_td_obj) :: cfg
        type(ghost) :: g

!cod$
        call my(cfg)
        g = cfg%o%g
        call glean(thy(cfg))
      end function

      function cfg_external(cfg) result(ext)
!doc$ function x_external(cfg) result(ext)
        type(config_td_obj) :: cfg
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
        type(config_td_obj) :: cfg
        type(fields_td_obj) :: fd
!       effects: Returns the fields component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%fields,fd)
        call glean(thy(cfg))
        call bequeath(thy(fd))
      end function

      function cfg_electrons(cfg) result(el)
!doc$ function x_electrons(cfg) result(el)
        type(config_td_obj) :: cfg
        type(electrons_td_obj) :: el
!       effects: Returns the electrons component of cfg.

!cod$
        call my(cfg)
        call my(cfg%o%electrons,el)
        call glean(thy(cfg))
        call bequeath(thy(el))
      end function

      function cfg_config_sc(cfg) result(cfg_sc)
!doc$ function x_config_sc(cfg) result(cfg_sc)
        type(config_td_obj) :: cfg
        type(config_sc_obj) :: cfg_sc
!       effects: Returns the self-consistent configuration component of cfg.

!cod$
        call my(cfg)
        call update(cfg%o%cfg_sc,cfg%o%external)  ; if (error()) goto 100
        call my(cfg%o%cfg_sc,cfg_sc)
100     call glean(thy(cfg))
        call bequeath(thy(cfg_sc))
        if (error("Exit config_td_mod::cfg_config_sc")) continue
      end function

      function cfg_adiabatic_states(cfg) result(el)
!doc$ function x_adiabatic_states(cfg) result(el)
        type(config_td_obj) :: cfg
        type(electrons_sc_obj) :: el
!       effects: Returns the component of cfg containing the current eigenstates.

!cod$

        logical :: done
        integer :: lc
        real(double) :: rn

        type(atoms_obj)   :: at
        type(crystal_obj) :: cr
        type(external_obj) :: ext

        call my(cfg)

!       N.M The following nonsense is needed to change the ghost on external
!           and force update(cfg%o%adiabatic_el,...) to restart the eigensolver
        call my(cfg%o%external,ext)
        call my(x_crystal(ext),cr)
        call my(x_atoms(cr),at)
        call move(at,x_position(at)); if (error()) goto 100
        call update(cr,atoms=at); if (error()) goto 100
        call update(ext,cr)

        if (i_access(diaryfile())) then
           write(x_unit(diaryfile()),'(/,t4,"Solving for adiabatic eigenstates")')
        end if
        if (i_access(output)) then
           write(x_unit(output),'("Solving for adiabatic eigenstates")')
        end if

        lc = 0
        done = .false.
        do while (.not. done)
          if (error(user_abort(),"USER INITIATED ABORT")) goto 100
          lc = lc + 1
          call update(cfg%o%adiabatic_el,ext,x_potential(cfg%o%fields)) ; if (error()) goto 100
          rn = x_residual_norm(cfg%o%adiabatic_el)
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
          done = x_converged(cfg%o%adiabatic_el)
!         N.M. The maximum number of steps should probably be a read from argvf instead of a fixed number
          if (lc >= 40) done = .true.
        end do

        call glean(thy(at))
        call glean(thy(cr))
        call glean(thy(ext))

        call my(cfg%o%adiabatic_el,el)
100     call glean(thy(cfg))
        call bequeath(thy(el))
        if (error("Exit config_td_mod::cfg_adiabatic_states")) continue
      end function

      function cfg_forces(cfg) result(f)
!doc$ function x_forces(cfg) result(f)
        type(config_td_obj) :: cfg
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
        if (error("Exit config_td_mod::cfg_forces")) continue
      end function

      function cfg_pressure(cfg) result(p)
!doc$ function x_pressure(cfg) result(pressure)
        type(config_td_obj) :: cfg
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
        if (error("Exit config_td_mod::cfg_pressure")) continue
      end function

      function cfg_stress_tensor(cfg) result(s)
!doc$ function x_stress_tensor(cfg) result(s)
        type(config_td_obj) :: cfg
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
        if (error("Exit config_td_mod::cfg_stress_tensor")) continue
      end function


      function cfg_cell_energy(cfg) result(ce)
!doc$ function x_cell_energy(cfg) result(ce)
        type(config_td_obj) :: cfg
        real(double) :: ce
!       effects:  Returns the cell energy.
!cod$
        call my(cfg)
        ce = x_energy(cfg%o%electrons) + x_energy(cfg%o%fields)
        call glean(thy(cfg))
      end function


      function cfg_num_sc_steps(cfg) result(num_steps)
!doc$ function x_num_sc_steps(cfg) result(num_steps)
        type(config_td_obj) :: cfg
        integer :: num_steps
!       effects:  Returns the number of self-consistent steps taken in the last time step
!cod$
        call my(cfg)
        num_steps = cfg%o%num_sc_steps
        call glean(thy(cfg))
      end function

           
      subroutine diary_cfg(cfg)
!doc$ subroutine diary(cfg)
        type(config_td_obj) :: cfg
!       modifies: Output stream
!       effects: Prints cfg information.
!       errors: Passes errors.
!cod$
        call my(cfg)
        if (cfg%o%have_forces) call diary_forces_i(cfg%o)
        if (cfg%o%have_pressure) call diary_pressure_i(cfg%o)
        if (cfg%o%have_stress_tensor) call diary_stress_tensor_i(cfg%o)
        call diary(cfg%o%electrons) ; if (error()) goto 100
        call diary_energy(cfg) ; if (error()) goto 100
100     call glean(thy(cfg))

        if (i_access(diaryfile())) call flushbuf(diaryfile())

        if (error("Exit config_td_mod::diary_cfg")) continue
      end subroutine

      subroutine diary_atom_step(cfg,step)
!doc$ subroutine diary(cfg,step)
        type(config_td_obj) :: cfg
        integer, intent(inout) :: step
!       modifies: Output stream and cfg if cfg%o%have_forces is .false.
!       effects: Prints cfg information relevant to a particular step.
!       errors: Passes errors.

!cod$
        call my(cfg)
        call diary(cfg%o%external,step) ; if (error()) goto 100
        call diary_forces(cfg) ; if (error()) goto 100
        call diary_energy(cfg,all=.false.) ; if (error()) goto 100
100     call glean(thy(cfg))
        if (error("Exit config_td_mod::diary_atom_step")) continue
      end subroutine

      subroutine diary_energy_td(cfg,all)
!doc$ subroutine diary_energy(cfg,all)
        type(config_td_obj) :: cfg
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
          call flushbuf(diaryfile())
        end if

        call glean(thy(cfg))

      end subroutine

      subroutine diary_forces_td(cfg)
!doc$ subroutine diary_forces(cfg)
        type(config_td_obj) :: cfg
!       modifies: Output stream and cfg if cfg%o%have_forces is .false. (Note: this is an optimization).
!       effects: Writes the forces to the diary file.
!       errors: Passes errors.

!cod$
        integer :: ia, na
        call my(cfg)
        if (.not.cfg%o%have_forces) call forces_i(cfg%o) ; if (error()) goto 100
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
        end if
100     call glean(thy(cfg))
        if (error("Exit config_td_mod::diary_forces")) continue
      end subroutine

      subroutine diary_pressure_td(cfg)
!doc$ subroutine diary_pressure(cfg)
        type(config_td_obj) :: cfg
!       modifies: Output stream and cfg if cfg%o%have_pressure is .false. (Note: this is an optimization).
!       effects: Writes the pressure to the diary file.
!       errors: Passes errors.

!cod$
        call my(cfg)
        if (.not.cfg%o%have_pressure) call pressure_i(cfg%o) ; if (error()) goto 100
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Pressure: ",f14.4,3x,"kbar")') cfg%o%pressure*147105.164_double
        end if
100     call glean(thy(cfg))
        if (error("Exit config_td_mod::diary_pressure")) continue
      end subroutine

      subroutine diary_stress_tensor_td(cfg)
!doc$ subroutine diary_stress_tensor(cfg)
        type(config_td_obj) :: cfg
!       modifies: Output stream and cfg if cfg%o%have_stress_tensor is .false. (Note: this is an optimization).
!       effects: Writes the stress tensor to the diary file.
!       errors: Passes errors.

!cod$
        call my(cfg)
        if (.not.cfg%o%have_stress_tensor) call stress_tensor_i(cfg%o) ; if (error()) goto 100
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Stress Tensor:")')
          write(x_unit(diaryfile()),'(t8,3f12.4,3x,"kbar")') cfg%o%stress_tensor(1,:)*147105.164_double
          write(x_unit(diaryfile()),'(t8,3f12.4)')           cfg%o%stress_tensor(2,:)*147105.164_double
          write(x_unit(diaryfile()),'(t8,3f12.4)')           cfg%o%stress_tensor(3,:)*147105.164_double
        end if
100     call glean(thy(cfg))
        if (error("Exit config_td_mod::diary_stress_tensor")) continue
      end subroutine

      subroutine decompose_cfg(cfg)
!doc$ subroutine decompose(cfg)
        type(config_td_obj) :: cfg
!       effects: Decomposes the wave functions into s, p, & d spherical harmonics around user-defined and atom sites.
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

        call start_timer("config_td: decompose")

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

600     if (.not.error()) call stop_timer("config_td: decompose")

700     call glean(thy(cfg))

        if (error("Exit config_td_mod::decompose_cfg")) continue

      end subroutine

      subroutine write_restart_cfg(cfg,nrestf)
!doc$ subroutine write_restart(cfg,mode_in)
        type(config_td_obj) :: cfg
        type(tagio_obj) :: nrestf
!       effects: Writes restart information to the nrestf
!       errors: Passes errors.
!cod$

        call my(cfg)
        call my(nrestf)

        call write_restart(cfg%o%external,nrestf) ; if (error()) goto 100
        call write_restart(cfg%o%fields,nrestf)   ; if (error()) goto 100
        call write_restart(cfg%o%electrons,nrestf)

100     call glean(thy(cfg))
        call glean(thy(nrestf))

        if (error("Exit config_td_mod::write_restart_cfg")) continue

      end subroutine

! private routines

      subroutine standard_portal_i(cfgr)
        type(config_td_rep) :: cfgr

        call read_instructions_i(cfgr)

        ! construct the self-consistent config
        call my(config_sc(),cfgr%cfg_sc) ; if (error()) goto 200

        ! construct the adiabatic eigenstates
        call my(x_electrons(cfgr%cfg_sc),cfgr%adiabatic_el)

        ! initialize external, electrons and fields
        call my(x_external(cfgr%cfg_sc),cfgr%external)
        call my(electrons_td(x_electrons(cfgr%cfg_sc),cfgr%external),cfgr%electrons) ; if (error()) goto 100
        call my(fields_td(cfgr%external,cfgr%electrons),cfgr%fields)
100     call sync_configuration_errors()

        call diary_construction_i()

        allocate( cfgr%forces(3,x_n_atoms(x_atoms(x_crystal(cfgr%external)))) )
        cfgr%have_forces = .false.
        cfgr%have_pressure = .false.
        cfgr%have_stress_tensor = .false.

        cfgr%num_sc_steps = 0

        call my(grid(x_layout(cfgr%external), CONFIG),cfgr%guess_pot); if (error()) goto 200

200     if (error("Exit config_td_mod::standard_portal_i")) continue

      end subroutine

      subroutine restart_portal_i(cfgr,restf)
        type(config_td_rep) :: cfgr
        type(tagio_obj) :: restf

        integer :: ios

        call my(restf)

        call read_instructions_i(cfgr) ; if (error()) goto 300

        ! read/construct external, fields and electrons
        call my(external(restf),cfgr%external) ; if (error()) goto 100
        call my(fields_td(cfgr%external,restf=restf),cfgr%fields) ; if (error()) goto 100
        call my(electrons_td(restf),cfgr%electrons)
        call my(config_sc(),cfgr%cfg_sc) ; if (error()) goto 100
        call my(x_electrons(cfgr%cfg_sc),cfgr%adiabatic_el)

100     call glean(thy(restf))

        call sync_configuration_errors() ; if (error()) goto 300

        call diary_construction_i()

        call sync_configuration_errors() ; if (error()) goto 300

        ! compute forces, pressure and stress tensor
        allocate( cfgr%forces(3,x_n_atoms(x_atoms(x_crystal(cfgr%external)))) )
        cfgr%have_forces = .false.
        if (cfgr%compute_forces) call forces_i(cfgr) ; if (error()) goto 200
        cfgr%have_pressure = .false.
        if (cfgr%compute_pressure) call pressure_i(cfgr) ; if (error()) goto 200
        cfgr%have_stress_tensor = .false.
        if (cfgr%compute_stress_tensor) call stress_tensor_i(cfgr) ; if (error()) goto 200
200     call sync_configuration_errors()

300     call glean(thy(restf))

        if (error("Exit config_td_mod::restart_portal_i")) continue

      end subroutine

      subroutine read_instructions_i(cfgr)
        type(config_td_rep) :: cfgr

        logical :: found
        character(line_len) :: tag
        real(double) :: tol
        integer :: max_steps

        call arglc("forces",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          cfgr%compute_forces = .true.
        case ("off",".false.")
          cfgr%compute_forces = .false.
        case default
          if (error(.true.,"ERROR: forces tag is not recognized")) goto 100
        end select

        call arglc("null_residual_force",tag,found)
        if (.not.found) tag = "on"
        select case (trim(tag))
        case ("on")
          cfgr%null_residual_force = .true.
        case ("off")
          cfgr%null_residual_force = .false.
        case default
          if (error(.true.,"ERROR: null_residual_force tag is not recognized")) goto 100
        end select

        call arglc("pressure",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          cfgr%compute_pressure = .true.
        case ("off",".false.")
          cfgr%compute_pressure = .false.
        case default
          if (error(.true.,"ERROR: pressure tag is not recognized")) goto 100
        end select

        call arglc("stress_tensor",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.")
          cfgr%compute_stress_tensor = .true.
        case ("off",".false.")
          cfgr%compute_stress_tensor = .false.
        case default
          if (error(.true.,"ERROR: stress tensor tag is not recognized")) goto 100
        end select

        call arglc("tddft_propagation_type",tag,found)
        if (.not.found) tag = "simple"
        select case (trim(tag))
        case ("simple")
          cfgr%propagation_type = SIMPLE
        case ("self-consistent","selfconsistent","sc")
          cfgr%propagation_type = SELFCONSISTENT
        case default
          if (error(.true.,"ERROR: propagation_type tag is not recognized")) goto 100
        end select

        call arg("tddft_propagation_tol",tol,found)
        if (.not.found) tol = default_tol
        cfgr%propagation_tol = tol

        call arg("tddft_max_sc_steps",max_steps,found)
        if (.not.found) max_steps = 10
        cfgr%max_sc_steps = max_steps

100     if (error("Exit config_sc_mod::read_instructions_i")) continue

      end subroutine



      subroutine diary_construction_i()

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,"Config object construction:")')
          write(x_unit(diaryfile()),'(/,t4,"Time-dependent calculation")')
        end if

      end subroutine

      subroutine own_i(cfg)
        type(config_td_obj) :: cfg
        type(config_td_obj) :: cfgt
        if (cfg%ref < cfg%o%ref) then
          allocate( cfgt%o )
          cfgt%o%ref = 0
          cfgt%o%g = cfg%o%g
          call my(cfg%o%external,cfgt%o%external)
          call my(cfg%o%fields,cfgt%o%fields)
          call my(cfg%o%electrons,cfgt%o%electrons)
          call my(cfg%o%cfg_sc,cfgt%o%cfg_sc)
          call my(cfg%o%adiabatic_el,cfgt%o%adiabatic_el)
          cfgt%o%have_forces = cfg%o%have_forces
          allocate( cfgt%o%forces(size(cfg%o%forces,1),size(cfg%o%forces,2)) )
          cfgt%o%forces = cfg%o%forces
          cfgt%o%have_pressure = cfg%o%have_pressure
          cfgt%o%pressure = cfg%o%pressure
          cfgt%o%have_stress_tensor = cfg%o%have_stress_tensor
          cfgt%o%stress_tensor = cfg%o%stress_tensor
          cfgt%o%propagation_type = cfg%o%propagation_type
          cfgt%o%num_sc_steps = cfg%o%num_sc_steps
          cfgt%o%propagation_tol = cfg%o%propagation_tol
          cfgt%o%max_sc_steps = cfg%o%max_sc_steps
          call my(cfg%o%guess_pot,cfgt%o%guess_pot)
          cfg%o%ref = cfg%o%ref - cfg%ref
          cfg%o => cfgt%o
          cfg%o%ref = cfg%o%ref + cfg%ref
        end if
      end subroutine

      function cell_energy_i(cfgr) result(ce)
        type(config_td_rep) :: cfgr
        real(double) :: ce
        ce = x_energy(cfgr%electrons) + x_energy(cfgr%fields)
      end function

      subroutine forces_i(cfgr)
        type(config_td_rep) :: cfgr
!       requires: cfgr%forces allocated.
        integer :: ia, na
        real(double), dimension(:,:), allocatable :: ft
        call start_timer("config_td: forces")
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
        cfgr%have_forces = .true.
100     if (allocated( ft )) deallocate( ft )
        if (error("Exit config_td_mod::forces_i")) continue
        if (.not.error()) call stop_timer("config_td: forces")
      end subroutine

      subroutine pressure_i(cfgr)
        type(config_td_rep) :: cfgr
        real(double) :: p, pt
        call start_timer("config_td: pressure")
        p = 0.0_double
        call pressure(cfgr%fields,pt) ; if (error()) goto 100
        p = p + pt
        call pressure(cfgr%electrons,pt) ; if (error()) goto 100
        p = p + pt
        cfgr%pressure = p
        cfgr%have_pressure = .true.
100     if (error("Exit config_td_mod::pressure_i")) continue
        if (.not.error()) call stop_timer("config_td: pressure")
      end subroutine

      subroutine stress_tensor_i(cfgr)
        type(config_td_rep) :: cfgr
        real(double), dimension(:,:), allocatable :: s, st
        call start_timer("config_td: stress tensor")
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
        if (error("Exit config_td_mod::stress_tensor_i")) continue
        if (.not.error()) call stop_timer("config_td: stress tensor")
      end subroutine

      subroutine diary_forces_i(cfgr)
        type(config_td_rep) :: cfgr
!       requires: cfgr%have_forces = .true.
        integer :: ia, na
        na = size(cfgr%forces,2)
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Atomic forces:")')
          write(x_unit(diaryfile()),'(20x,"atom",11x,"Fx",12x,"Fy",12x,"Fz")')
          write(x_unit(diaryfile()),'(19x,"---------------------------------------------------")')
          do ia = 1,na
            if (ia == 1) then
              write(x_unit(diaryfile()),'(19x,i4,3x,sp,3(4x,f10.6),3x,"Ryd/Bohr")') ia, cfgr%forces(:,ia)
            else
              write(x_unit(diaryfile()),'(19x,i4,3x,sp,3(4x,f10.6))') ia, cfgr%forces(:,ia)
            end if
          end do
        end if
      end subroutine

      subroutine diary_pressure_i(cfgr)
        type(config_td_rep) :: cfgr
!       requires: cfgr%have_pressure = .true.
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Pressure: ",f14.4,3x,"kbar")') cfgr%pressure*147105.164_double
        end if
      end subroutine

      subroutine diary_stress_tensor_i(cfgr)
        type(config_td_rep) :: cfgr
!       requires: cfgr%have_stress_tensor = .true.
        if (i_access( diaryfile() )) then
          write(x_unit(diaryfile()),'(/,t4,"Stress Tensor:")')
          write(x_unit(diaryfile()),'(t8,3f12.4,3x,"kbar")') cfgr%stress_tensor(1,:)*147105.164_double
          write(x_unit(diaryfile()),'(t8,3f12.4)')           cfgr%stress_tensor(2,:)*147105.164_double
          write(x_unit(diaryfile()),'(t8,3f12.4)')           cfgr%stress_tensor(3,:)*147105.164_double
        end if
      end subroutine

      subroutine take_sc_step_i(cfgr)
        type(config_td_rep)  :: cfgr
        ! Local Vars  ---------------------------------------
        type(gen_potential_obj) :: gen_pot_start
        type(gen_potential_obj) :: gen_pot_half
!        type(gen_potential_obj) :: gen_pot_end
        type(grid_obj)          :: grid_pot_start
        type(grid_obj)          :: grid_pot_half
        type(grid_obj)          :: grid_pot_end

        type(electrons_td_obj)  :: elec_start

        real(double),pointer    :: pot_start(:,:,:)
        real(double),pointer    :: pot_half(:,:,:)
        real(double),pointer    :: pot_end(:,:,:)
        real(double),pointer    :: guess_pot(:,:,:)

        real(double) :: external_vecpot(3)
        real(double) :: resid_local
        real(double) :: resid_global
        real(double) :: sumpot_local
        real(double) :: sumpot_global
        integer      :: istep

        !====================================================

        !RMH  Need to fix this once the external vecpot is back in.
        !** Obtain the external vector potential
        !external_vecpot = x_external_vecpot( &
        !     x_tprop(x_stepper(x_wavefunctions(cfg%o%cfg_static%o%elec,1))))

       !** Cache the electrons and potential objects at t=t0
        call my(cfgr%electrons,elec_start) ; if (error()) goto 200
        call my(x_potential(cfgr%fields),gen_pot_start) ; if (error()) goto 400

        !** Initialize the various potential related quantities at the half_time_step
        call my(gen_pot_start,gen_pot_half)

        !** Extract the guess at the potential for this time step.
        if (x_type(cfgr%guess_pot) /= EMPTY_KIND) then
           !** update the general potential at the half time step.
           call update(gen_pot_half,gpot=cfgr%guess_pot) ; if (error()) goto 400
        end if

        !** Extract the grid potential.
        call my(x_grid_potential(gen_pot_half),grid_pot_half)

!        !** Update the gen_pot_half object with the guess potential.
!        call update(gen_pot_half,gpot=grid_pot_half) ; if (error()) goto 400

        !** Time step the underlying wavefunctions using the guess potential.
        call update(cfgr%electrons,cfgr%external,gen_pot_half) ; if (error()) goto 100

        !** Update the underlying fields with respect to the wavefunctions at t=t0+dt
        call update(cfgr%fields,cfgr%external,cfgr%electrons) ; if (error()) goto 100

        !** Extract the grid potential back at t=t0
        call my(x_grid_potential(gen_pot_start),grid_pot_start) ; if (error()) goto 600

        !** Extract the grid potential at t=t0+dt
        call my(x_grid_potential(x_potential(cfgr%fields)),grid_pot_end); if (error()) goto 600

        !** Calculate the potential at the half step
        call take(pot_start, grid_pot_start, RD_KIND) ; if (error()) goto 700
        call take(pot_half,  grid_pot_half,  RD_KIND) ; if (error()) goto 700
        call take(pot_end,   grid_pot_end,   RD_KIND) ; if (error()) goto 700

        !** Initialize sumpot which is a variable that stores the average magnitude of the potential
        sumpot_local = sum(abs(pot_start))
        call allreduce(CONFIG,MPI_SUM,sumpot_local,sumpot_global)


        !** Calculate the residual error between the half potential and the start potential
        resid_local = sum(abs((pot_start+pot_end)/2.0_double - pot_half))/sumpot_global
        call allreduce(CONFIG,MPI_SUM,resid_local,resid_global)

        pot_half = (pot_start + pot_end)/2.0_double

!** Calculate the residual error between the half potential and the start potential
!resid_local = sum(abs(pot_start - pot_half))/sumpot_global
!call allreduce(MPI_SUM,resid_local,resid_global)
!write(*,*) 'initial residual (pot_start - pot_half) = ', resid_global

!        call put(pot_start,grid_pot_start,RD_KIND) ; if (error()) goto 700

        istep = 1

        do while ((istep < cfgr%max_sc_steps) .and. (cfgr%propagation_tol < resid_global))

           call put(pot_half, grid_pot_half, RD_KIND) ; if (error()) goto 700
           call put(pot_end,  grid_pot_end,  RD_KIND) ; if (error()) goto 700

           !** update the general potential at the half time step.
           call update(gen_pot_half,gpot=grid_pot_half) ; if (error()) goto 700

           !** Reset the electrons object to the electrons object at t=t0
           cfgr%electrons = elec_start

           !** Time step the underlying wavefunctions using the potential at t=t0+dt/2
           call update(cfgr%electrons,cfgr%external,gen_pot_half) ; if (error()) goto 100

           !** Update the underlying fields with respect to the wavefunctions at t=t0+dt
           call update(cfgr%fields,cfgr%external,cfgr%electrons) ; if (error()) goto 100

!           !** Extract the potential from the fields object at t=t0+dt
!           call update(gen_pot_end,x_grid_potential(x_potential(cfgr%fields)))
!           if (error()) goto 400

           !** Extract the grid potential at t=t0+dt
           grid_pot_end = x_grid_potential(x_potential(cfgr%fields)); if (error()) goto 100

           !call glean(thy(grid_pot_end))
           !call my(x_grid_potential(gen_pot_end), grid_pot_end ) ; if (error()) goto 700

           call take(pot_end,grid_pot_end,RD_KIND)
           call take(pot_half,grid_pot_half,RD_KIND)


!write(*,*) 'istep ', istep
!write(*,*) '  pot_start(3,3,3) = ', pot_start(3,3,3)
!write(*,*) '   pot_half(3,3,3) = ', pot_half(3,3,3)
!write(*,*) '    pot_end(3,3,3) = ', pot_end(3,3,3)
!write(*,*) '      resid(3,3,3) = ', 0.5*(pot_start(3,3,3) + pot_end(3,3,3)) - pot_half(3,3,3)

           !** Calculate the residual error between the half potential and the start potential
           resid_local = sum(abs((pot_start+pot_end)/2.0_double - pot_half))/sumpot_global
           call allreduce(CONFIG,MPI_SUM,resid_local,resid_global)
!write(*,*) '    residual (pot_start - pot_half) = ', resid_global

           pot_half = (pot_start + pot_end)/2.0_double

           istep = istep + 1

        end do

        !** return the number of time propagation operations (add one to account for
        !      the intial propagation operation step before the loop.
        cfgr%num_sc_steps = istep

        !** Calculate and store a guess for the potential at the next half time step
        if (x_type(cfgr%guess_pot) == EMPTY_KIND) then
           call alloc(guess_pot,x_layout(grid_pot_half),RD_KIND,CONFIG)
           if (error()) goto 100
           guess_pot = pot_end
!           guess_pot = 2.0_double*pot_end - pot_half
!           guess_pot = 3.0_double*pot_end - 3.0_double*pot_half + pot_start
        else  !** guess_pot not empty
           call take(guess_pot,cfgr%guess_pot,RD_KIND)
           guess_pot = pot_end
!           guess_pot = 2.0_double*pot_end - pot_half
!           guess_pot = 3.0_double*pot_end - 3.0_double*pot_half + pot_start
        end if

        call put(guess_pot,cfgr%guess_pot,RD_KIND)

        call put(pot_half, grid_pot_half, RD_KIND) ; if (error()) goto 700
        call put(pot_end,  grid_pot_end,  RD_KIND) ; if (error()) goto 700
        call put(pot_start,grid_pot_start,RD_KIND) ; if (error()) goto 700

        if (associated(pot_start)) deallocate(pot_start)
        if (associated(pot_half)) deallocate(pot_half)
        if (associated(pot_end)) deallocate(pot_end)
        nullify(pot_start,pot_half,pot_end)
        call glean(thy(grid_pot_half))
700     call glean(thy(grid_pot_end))
600     call glean(thy(grid_pot_start))
        call glean(thy(gen_pot_half))
!500     call glean(thy(gen_pot_end))
400     call glean(thy(gen_pot_start))
200     call glean(thy(elec_start))
100     if (error("config_td_mod::take_sc_step_i - Exiting")) continue

      end subroutine take_sc_step_i

      end module

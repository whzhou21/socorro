!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module transition_state_mod
!doc$ module transition_state_mod

!     transition_state_mod implements procedures for determining transitions states, namely the dimer
!     method and the nudged-elastic-band method as developed by Graeme Henkelman, University of Texas,
!     Austin, and an algorithm to find transitions states for carrier capture developed by Normand Modine
!     at Sandia National Laboratories. The portal to these methods is the transition_state function.
!     Note that dimer_obj, neb_obj and ccts_obj are private and thus can only be used from within this module.

      use kind_mod
      use ghost_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use arg_mod
      use path_mod
      use math_mod
      use interrupt_mod
      use config_sc_mod
      use external_mod
      use crystal_mod
      use lattice_mod
      use shortcuts_mod

!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = WORLD

      integer, parameter :: ALL_IMAGES = 1
      integer, parameter :: CLIMBING_IMAGE =  2

      type :: dimer_rep
        integer :: ref
        type(ghost) :: g
        real(double), dimension(:,:), pointer :: pos_cart_01, pos_cart_02   ! atom positions of the configurations
        real(double), dimension(:,:), pointer :: force_01, force_02         ! atomic forces of the configurations
        real(double), dimension(:,:), pointer :: pos_cart_center            ! position of the dimer center
        real(double), dimension(:,:), pointer :: force_center               ! force on the dimer center
        real(double), dimension(:,:), pointer :: torque                     ! torque on the dimer
        real(double), dimension(:,:), pointer :: eff_force
        real(double), dimension(:,:), pointer :: separation_dir             ! normalized separation of the configurations
        real(double), dimension(:,:), pointer :: torque_dir                 ! normalized direction of the torque
        real(double) :: cell_energy_01, cell_energy_02                      ! energies of the configurations
        real(double) :: energy_center                                       ! energy of the dimer center
        real(double) :: force_center_mag                                    ! magnitude of the dimer force
        real(double) :: torque_mag                                          ! magnitude of the dimer torque
        real(double) :: separation_mag                                      ! user supplied dimer separation
        real(double) :: min_curvature
        real(double) :: force_tol                                           ! tolerance for translation optimization
        real(double) :: torque_tol                                          ! tolerance for rotation optimization
        integer :: n_steps                                                  ! number of steps taken during optimization
        integer :: max_steps                                                ! maximum number of steps taken for optimization
      end type
      
      type :: dimer_obj
        integer :: ref
        type(ghost) :: g
        type(dimer_rep), pointer :: o
      end type

      type :: neb_rep
        integer :: ref
        type(ghost) :: g
        real(double), dimension(:,:,:), pointer :: pos_cart_end ! atom positions of the end configurations
        real(double), dimension(:,:,:), pointer :: pos_cart     ! atom positions of the moving configurations
        real(double), dimension(:,:,:), pointer :: force        ! atomic forces of the mofing configurations
        real(double), dimension(:), pointer :: cell_energy      ! energies of the moving configurations
        real(double) :: spring                                  ! spring constant
        real(double) :: force_tol                               ! tolerance for the optimization
        integer :: n_steps                                      ! number of steps taken during optimization
        integer :: max_steps                                    ! maximum number of steps for optimization
        integer :: force_test                                   ! optimization test criterion
        logical :: climbing_image                               ! switch to use the climbing image method
        integer :: max_energy_image                             ! configuration with the highest energy
        real(double) :: max_energy                              ! energy of configuration max_energy_image
        real(double) :: rms_force                               ! rms force
        integer :: max_rms_force_image                          ! configuration with the largest rms force
        real(double) :: max_rms_force                           ! rms force of configuration max_rms_force_image
      end type
      
      type :: neb_obj
        integer :: ref
        type(ghost) :: g
        type(neb_rep), pointer :: o
      end type

      type :: ccts_rep
        integer :: ref
        type(ghost) :: g
        real(double), dimension(:,:), pointer :: pos_cart_01, pos_cart_02   ! atom positions of the configurations
        real(double), dimension(:,:), pointer :: force_01, force_02 ! atomic forces of the configurations
        real(double), dimension(:,:), pointer :: last_pos_cart
        real(double), dimension(:,:), pointer :: force_sum ! sum of the forces
        real(double), dimension(:,:), pointer :: force_diff ! difference of the forces
        real(double), dimension(:,:), pointer :: eff_force ! effective force on the configuration
        real(double) :: cell_energy_01, cell_energy_02 ! energies of the configurations
        real(double) :: level ! current difference of the energies
        real(double) :: target_level ! target difference of the energies
        real(double) :: level_err ! error in the current difference in energies
        real(double) :: force_sum_mag_sq ! magnitude of the force sum
        real(double) :: force_diff_mag_sq ! magnitude of the force difference
        real(double) :: eff_force_rms ! rms magnitude of the effective force components
        real(double) :: force_tol ! tolerance for contrained relaxation at a fixed energy difference
        real(double) :: level_tol ! tolerance for maintaining the energy difference during relaxation
        real(double) :: step ! force scaling factor for steepest descent relaxation
        real(double) :: newton_scaling ! scaling factor for Newton step enforcing level constraint
        real(double) :: move_length ! the initial distance between configurations in the path
        real(double) :: move_noise ! a random component added to each coordinate before each relaxation
        real(double) :: path_length ! the total distance traveled along the path
        integer :: n_steps ! number of steps taken during optimization
        integer :: max_steps ! maximum number of steps taken for optimization
        integer :: n_moves ! number of moves taken along the path
        integer :: max_moves ! maximum number of moves used to determine the path

      end type

      type :: ccts_obj
        integer :: ref
        type(ghost) :: g
        type(ccts_rep), pointer :: o
      end type

!doc$
      public :: transition_state

!cod$
      contains

      function transition_state(cfg) result(changed)
!doc$ function transition_state(cfg) result(changed)
        type(config_sc_obj) :: cfg
        logical :: changed
!       effects: Finds the transition state using the specified method. changed indicates if cfg has changed.
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag
        type(dimer_obj) :: dmr
        type(neb_obj) :: neb
        type(ccts_obj) :: ccts

        call my(cfg)
        
        changed = .false.

        call arglc("ts_method",tag,found)
        if (.not.found) tag = "none"
        select case (trim(tag))
        case ("dimer","d")
          changed = .true.
          call my_new_dmr(constructor_dmr(cfg),dmr) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),"('dimer construction complete')")
          if (i_access(output)) write(x_unit(output),"('dimer construction complete')")
          call glean_dmr(thy_dmr(dmr))
        case ("neb","n")
          changed = .true.
          call my_new_neb(constructor_neb_i(cfg),neb) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),"('neb construction complete')")
          if (i_access(output)) write(x_unit(output),"('neb construction complete')")
          call glean_neb(thy_neb(neb))
        case ("electron_capture","hole_capture")
          changed = .true.
          call my_new_ccts(constructor_ccts(cfg),ccts) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),"('ccts construction complete')")
          if (i_access(output)) write(x_unit(output),"('ccts construction complete')")
          call glean_ccts(thy_ccts(ccts))
        end select

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::transition_state")) continue
        
      end function

      subroutine my_dmr(dmr)
        type(dimer_obj) :: dmr
        dmr%ref = dmr%ref + 1
        dmr%o%ref = dmr%o%ref + 1
      end subroutine

      subroutine my_new_dmr(dmri,dmr)
      type(dimer_obj) :: dmri, dmr
        dmr%ref = 1
        dmr%o => dmri%o
        dmr%o%ref = dmr%o%ref + 1
      end subroutine

      function thy_dmr(dmr) result(dmro)
        type(dimer_obj) :: dmr, dmro
        dmr%ref = dmr%ref - 1
        dmr%o%ref = dmr%o%ref - 1
        dmro%ref = dmr%ref
        dmro%o => dmr%o
      end function

      subroutine glean_dmr(dmr)
        type(dimer_obj) :: dmr
        if (dmr%o%ref < 1) then
          if (associated( dmr%o%pos_cart_01 )) deallocate( dmr%o%pos_cart_01 )
          if (associated( dmr%o%pos_cart_02 )) deallocate( dmr%o%pos_cart_02 )
          if (associated( dmr%o%pos_cart_center )) deallocate( dmr%o%pos_cart_center )
          if (associated( dmr%o%separation_dir )) deallocate( dmr%o%separation_dir )
          if (associated( dmr%o%force_01 )) deallocate( dmr%o%force_01 )
          if (associated( dmr%o%force_02 )) deallocate( dmr%o%force_02 )
          if (associated( dmr%o%force_center )) deallocate( dmr%o%force_center )
          if (associated( dmr%o%eff_force )) deallocate( dmr%o%eff_force )
          if (associated( dmr%o%torque )) deallocate( dmr%o%torque )
          if (associated( dmr%o%torque_dir )) deallocate( dmr%o%torque_dir )
          deallocate( dmr%o )
        end if
      end subroutine

      subroutine bequeath_dmr(dmr)
        type(dimer_obj) :: dmr
      end subroutine

      function constructor_dmr(cfg) result(dmr)
        type(config_sc_obj) :: cfg
        type(dimer_obj) :: dmr
!       effects: Constructs and optimizes a dmr.

        logical :: found
        character(line_len) :: tag
        integer :: na

        call my(cfg)

        dmr%ref = 0
        allocate( dmr%o )
        dmr%o%g = x_ghost()
        dmr%o%ref = 0

        nullify( dmr%o%pos_cart_01, dmr%o%pos_cart_02, dmr%o%pos_cart_center, dmr%o%separation_dir )
        nullify( dmr%o%force_01, dmr%o%force_02, dmr%o%force_center, dmr%o%eff_force, dmr%o%torque, dmr%o%torque_dir )

        if (error(mpi_nconfigs() /= 2,"ERROR: number of configurations /= 2")) goto 100

        call arg("dimer_force_tol",dmr%o%force_tol,found)
        if (.not.found) dmr%o%force_tol = 0.001_double
        if (error(dmr%o%force_tol < 0.0_double,"ERROR: force_tol < 0")) goto 100

        call arg("dimer_torque_tol",dmr%o%torque_tol,found)
        if (.not.found) dmr%o%torque_tol = 0.0001_double
        if (error(dmr%o%torque_tol < 0.0_double,"ERROR: torque_tol < 0")) goto 100

        call arg("dimer_separation",dmr%o%separation_mag,found)
        if (.not.found) dmr%o%separation_mag = 0.01_double
        if (error(dmr%o%separation_mag < 0.0_double,"ERROR: dimer_separation < 0")) goto 100

        call arg("dimer_max_steps",dmr%o%max_steps,found)
        if (.not.found) dmr%o%max_steps = 100
        if (error(dmr%o%max_steps < 0,"ERROR: max_steps < 0")) goto 100
        dmr%o%n_steps = 0

        na = x_n_atoms(cfg)
        allocate( dmr%o%pos_cart_01(3,na), dmr%o%pos_cart_02(3,na), dmr%o%pos_cart_center(3,na), dmr%o%separation_dir(3,na) )
        allocate( dmr%o%force_01(3,na), dmr%o%force_02(3,na), dmr%o%force_center(3,na), dmr%o%eff_force(3,na) )
        allocate( dmr%o%torque(3,na), dmr%o%torque_dir(3,na) )

        call initialize_dimer_i(dmr%o,cfg) ; if (error()) goto 100

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Completed dimer initialization")')
        if (i_access(output)) write(x_unit(output),'("Completed dimer initialization")')

        call diary_dimer_i(dmr%o)

        call arglc("dimer_optimizer",tag,found)
        if (.not.found) tag = "conjugate_gradient"
        select case (trim(tag))
        case ("steepest_descent","sd")
          call optimize_dimer_sd_i(dmr%o,cfg) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'("  Steepest-descent optimization")')
          if (i_access(output)) write(x_unit(output),'("  Steepest-descent optimization")')
        case ("conjugate_gradient","cg")
          call optimize_dimer_cg_i(dmr%o,cfg) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'("  Conjugate-gradient optimization")')
          if (i_access(output)) write(x_unit(output),'("  Conjugate-gradient optimization")')
        case default
          if (error(.true.,"ERROR: unrecognized dimer_optimizer")) goto 100
        end select

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Completed optimization of dimer")')
        if (i_access(output)) write(x_unit(output),'("Completed optimization of dimer")')

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::constructor_dmr")) continue
        
      end function

      subroutine initialize_dimer_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        real(double) :: separation_mag_init
        type(lattice_obj) :: lat

        call my(cfg)
        call my(x_lattice(x_crystal(x_external(cfg))),lat)
        
        call share_dmr_pfe_i(dmrr,cfg) ; if (error()) goto 100

        dmrr%separation_dir = dmrr%pos_cart_01 - dmrr%pos_cart_02
        call reduce_vectors_i(lat,dmrr%separation_dir)
        separation_mag_init = sqrt(atom_dot_product(dmrr%separation_dir))
        dmrr%separation_dir = dmrr%separation_dir/separation_mag_init

        dmrr%pos_cart_center = dmrr%pos_cart_01 - 0.5_double*separation_mag_init*dmrr%separation_dir
        dmrr%pos_cart_01 = dmrr%pos_cart_center + dmrr%separation_mag*dmrr%separation_dir
        dmrr%pos_cart_02 = dmrr%pos_cart_center - dmrr%separation_mag*dmrr%separation_dir

        call update_dimer_i(dmrr,cfg) ; if (error()) goto 100

100     call glean(thy(lat))

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::initialize_dimer_i")) continue

      end subroutine

      subroutine optimize_dimer_sd_i(dmrr,cfg) 
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        logical :: found
        integer :: na
        real(double) :: alpha, delta
        real(double), dimension(:,:), allocatable :: direction

        call my(cfg)

        call arg("dimer_sd_alpha",alpha,found)
        if (.not.found) alpha = 3.5_double
        if (error(alpha < 0.0_double,"ERROR: dimer_sd_alpha < 0")) goto 100

        na = size(dmrr%eff_force,2)
        allocate( direction(3,na) )

        call optimize_rotation_i(dmrr,cfg) ; if (error()) goto 100

        if (dmrr%n_steps > dmrr%max_steps) then
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of dimer terminated: n_steps = max_steps")')
          if (i_access(output)) write(x_unit(output),'("Optimization of dimer terminated: n_steps = max_steps")')
          goto 100
        end if

        do

          if (finish_test_dmr_i(dmrr)) exit

          call barrier(MOD_SCOPE)
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            goto 100
          end if

          if (i_access(output)) write(x_unit(output),'("Line minimization")')

          direction = dmrr%eff_force/sqrt(atom_dot_product(dmrr%eff_force))
          delta = alpha*sqrt(atom_dot_product(dmrr%eff_force))
          call translate_dimer_i(dmrr,cfg,direction,delta) ; if (error()) goto 100

          if (dmrr%n_steps > dmrr%max_steps) then
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of dimer terminated: n_steps = max_steps")')
            if (i_access(output)) write(x_unit(output),'("Optimization of dimer terminated: n_steps = max_steps")')
            goto 100
          end if

          call optimize_rotation_i(dmrr,cfg) ; if (error()) goto 100

          if (dmrr%n_steps > dmrr%max_steps) then
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of dimer terminated: n_steps = max_steps")')
            if (i_access(output)) write(x_unit(output),'("Optimization of dimer terminated: n_steps = max_steps")')
            goto 100
          end if

        end do

100     if (allocated( direction )) deallocate( direction )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_dimer_sd_i")) continue

      end subroutine

      subroutine optimize_dimer_cg_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        logical :: found
        integer :: na
        real(double) :: a, b, curvature, delta, finite_diff, gam, proj_f1, proj_f2
        real(double) :: max_move, step_size
        real(double), dimension(:,:), allocatable :: direction, direction_norm, direction_old, force_old

        na = size(dmrr%eff_force,2)
        allocate( direction(3,na) )
        allocate( direction_norm(3,na) )
        allocate( direction_old(3,na) )
        allocate( force_old(3,na) )

        call my(cfg)

        call arg("dimer_cg_finite_diff",finite_diff,found)
        if (.not.found) finite_diff = 0.01_double
        if (error(finite_diff < 0.0_double,"ERROR: dimer_cg_finite_diff < 0")) goto 100

        call arg("dimer_cg_max_move",max_move,found)
        if (.not.found) max_move = 6.0_double
        if (error(max_move < 0.0_double,"ERROR: max_move < 0")) goto 100

        direction = 0.0_double
        direction_old = 0.0_double
        force_old = 0.0_double

        call optimize_rotation_i(dmrr,cfg) ; if (error()) goto 100

        if (dmrr%n_steps > dmrr%max_steps) then
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of dimer terminated: n_steps = max_steps")')
          if (i_access(output)) write(x_unit(output),'("Optimization of dimer terminated: n_steps = max_steps")')
          goto 100
        end if

        do

          if (finish_test_dmr_i(dmrr)) exit

          call barrier(MOD_SCOPE)
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            goto 100
          end if

          if (i_access(output)) write(x_unit(output),'("Line minimization")')

          a = abs(atom_dot_product(dmrr%eff_force,force_old))
          b = atom_dot_product(force_old)

          if (a < 0.5_double*b) then
            gam = atom_dot_product(dmrr%eff_force,dmrr%eff_force - force_old)/b
          else
            gam = 0.0_double
          endif

          direction = dmrr%eff_force + gam*direction_old
          direction_old = direction
          force_old = dmrr%eff_force
          direction_norm = direction/sqrt(atom_dot_product(direction))          

          delta = finite_diff
          call translate_dimer_i(dmrr,cfg,direction,delta) ; if (error()) goto 100

          proj_f1 = atom_dot_product(force_old,direction_norm)
          proj_f2 = atom_dot_product(dmrr%eff_force,direction_norm)
          curvature = (proj_f1 - proj_f2)/finite_diff

          step_size = max_move
          if (curvature > 0.0_double) then
            step_size = proj_f1/curvature
          endif

          if (step_size > max_move) then
            step_size = max_move
          endif

          delta = step_size - finite_diff
          call translate_dimer_i(dmrr,cfg,direction,delta) ; if (error()) goto 100

          if (dmrr%n_steps > dmrr%max_steps) then
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of dimer terminated: n_steps = max_steps")')
            if (i_access(output)) write(x_unit(output),'("Optimization of dimer terminated: n_steps = max_steps")')
            goto 100
          end if

          call optimize_rotation_i(dmrr,cfg) ; if (error()) goto 100

          if (dmrr%n_steps > dmrr%max_steps) then
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of dimer terminated: n_steps = max_steps")')
            if (i_access(output)) write(x_unit(output),'("Optimization of dimer terminated: n_steps = max_steps")')
            goto 100
          end if

        end do

100     if (allocated( direction )) deallocate( direction )
        if (allocated( direction_norm )) deallocate( direction_norm )
        if (allocated( direction_old )) deallocate( direction_old )
        if (allocated( force_old )) deallocate( force_old )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_dimer_cg_i")) continue

      end subroutine

      subroutine update_dimer_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        real(double) :: dp, dp_a, dp_b

        call my(cfg)
        
        select case (mpi_myconfig())
        case (1)
          call update_config(cfg,dmrr%pos_cart_01)
        case (2)
          call update_config(cfg,dmrr%pos_cart_02)
        end select
        call sync_configuration_errors() ; if (error()) goto 100

        call share_dmr_fe_i(dmrr,cfg) ; if (error()) goto 100

        dmrr%force_center = 0.5_double*(dmrr%force_01 + dmrr%force_02)
        dmrr%force_center_mag = sqrt(atom_dot_product(dmrr%force_center))

        dp_a = atom_dot_product(dmrr%force_01,dmrr%separation_dir)
        dp_b = atom_dot_product(dmrr%force_02,dmrr%separation_dir)

        dmrr%torque = (dmrr%force_01 - dp_a*dmrr%separation_dir) - (dmrr%force_02 - dp_b*dmrr%separation_dir)

        dmrr%torque_dir = dmrr%torque
        dmrr%torque_mag = sqrt(atom_dot_product(dmrr%torque))
        dmrr%torque_dir = dmrr%torque/dmrr%torque_mag

        dp = atom_dot_product(dmrr%force_01-dmrr%force_02,dmrr%separation_dir)
        
        dmrr%energy_center = 0.5_double*(dmrr%cell_energy_01 + dmrr%cell_energy_02) + 0.25_double*dmrr%separation_mag*dp
        dmrr%min_curvature = -0.5_double*dp/dmrr%separation_mag

        dp = atom_dot_product(dmrr%force_center,dmrr%separation_dir)

        if (dmrr%min_curvature > 0.0_double) then
          dmrr%eff_force = -dp*dmrr%separation_dir
        else
          dmrr%eff_force = dmrr%force_center - 2.0_double*dp*dmrr%separation_dir
        end if

        dmrr%n_steps = dmrr%n_steps + 1

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::update_dimer_i")) continue

      end subroutine

      subroutine optimize_rotation_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        call my(cfg)

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'("optimize rotation")')
        if (i_access(output)) write(x_unit(output),'("optimize rotation")')

        call optimize_rotation_in_plane_i(dmrr,cfg) ; if (error()) goto 100
        if (dmrr%n_steps > dmrr%max_steps) goto 100
        do while (rotate_test_i(dmrr))
          call optimize_rotation_in_plane_i(dmrr,cfg) ; if (error()) goto 100
          if (dmrr%n_steps > dmrr%max_steps) goto 100
        end do

        call diary_dimer_i(dmrr)
        
100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_rotation_i")) continue

      end subroutine

      subroutine translate_dimer_i(dmrr,cfg,direction,delta)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg
        real(double), dimension(:,:), intent(in) :: direction
        real(double), intent(in) :: delta

        call my(cfg)

        if (i_access(diaryfile())) then
           write(x_unit(diaryfile()),"('translate center of the dimer by',g15.5)") delta
           write(x_unit(diaryfile()),"('in the direction:')") 
           write(x_unit(diaryfile()),"(2x,3g15.5)") direction
        end if
        if (i_access(output)) then
           write(x_unit(output),"('translate center of the dimer by',g15.5)") delta
           write(x_unit(output),"('in the direction:')") 
           write(x_unit(output),"(2x,3g15.5)") direction
        end if

        dmrr%pos_cart_center = dmrr%pos_cart_center + delta*direction
        dmrr%pos_cart_01 = dmrr%pos_cart_center + dmrr%separation_mag*dmrr%separation_dir
        dmrr%pos_cart_02 = dmrr%pos_cart_center - dmrr%separation_mag*dmrr%separation_dir
        
        call update_dimer_i(dmrr,cfg) ; if (error()) goto 100
        if (dmrr%n_steps > dmrr%max_steps) goto 100

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::translate_dimer_i")) continue

      end subroutine

      subroutine optimize_rotation_in_plane_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        integer :: na
        real(double), dimension(:,:), allocatable :: dir1, dir2
        real(double), parameter :: delta = 0.1_double
        real(double) :: f0, fdelta, theta0

        call my(cfg)

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'("optimize rotation in plane")')
        if (i_access(output)) write(x_unit(output),'("optimize rotation in plane")')
        
        na = size(dmrr%separation_dir,2)
        allocate( dir1(3,na), dir2(3,na) )

        call sync_configuration_errors() ; if (error()) goto 100

        dir1 = dmrr%separation_dir
        dir2 = dmrr%torque_dir

        f0 = atom_dot_product(dmrr%torque,dir2)

        call rotate_dimer_i(dmrr,cfg,dir1,dir2,delta) ; if (error()) goto 100
        if (dmrr%n_steps > dmrr%max_steps) goto 100

        fdelta = atom_dot_product(dmrr%torque,dir2)

        theta0 = 0.5_double*(delta - atan(delta*(f0+fdelta)/(fdelta-f0)))

        call rotate_dimer_i(dmrr,cfg,dir1,dir2,theta0) ; if (error()) goto 100
        if (dmrr%n_steps > dmrr%max_steps) goto 100

        call diary_dimer_i(dmrr)

        call sync_configuration_errors() ; if (error()) goto 100

100     if (allocated( dir1 )) deallocate( dir1 )
        if (allocated( dir2 )) deallocate( dir2 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_rotation_in_plane_i")) continue
        
      end subroutine

      subroutine rotate_dimer_i(dmrr,cfg,dir1,dir2,delta)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg
        real(double), dimension(:,:), intent(in) :: dir1, dir2
        real(double), intent(in) :: delta

        call my(cfg)

        if (i_access(diaryfile())) then
           write(x_unit(diaryfile()),"('rotate the dimer orientation by ',g15.5)") delta
           write(x_unit(diaryfile()),"('from the direction:')")
           write(x_unit(diaryfile()),"(2x,3g15.5)") dir1
           write(x_unit(diaryfile()),"('towards the direction:')")
           write(x_unit(diaryfile()),"(2x,3g15.5)") dir2
        end if

        dmrr%separation_dir = cos(delta)*dir1 + sin(delta)*dir2
        dmrr%pos_cart_01 = dmrr%pos_cart_center + dmrr%separation_mag*dmrr%separation_dir
        dmrr%pos_cart_02 = dmrr%pos_cart_center - dmrr%separation_mag*dmrr%separation_dir
        
        call update_dimer_i(dmrr,cfg) ; if (error()) goto 100
        if (dmrr%n_steps > dmrr%max_steps) goto 100

100     call glean(thy(cfg))

        if (error("Exit: transition_state_mod::rotate_dimer_i")) continue

      end subroutine

      function rotate_test_i(dmrr) result(go)
        type(dimer_rep) :: dmrr
        logical :: go

        integer :: na

        if (error()) then
          go = .false.
          goto 100
        end if
        
        go = .true.

        na = size(dmrr%torque,2)
        if (sqrt(atom_dot_product(dmrr%torque,dmrr%torque)/(3.0_double*real(na,double))) < dmrr%torque_tol ) go = .false.

100     if (error("Exit transition_state_mod::rotate_test_i")) continue

      end function

      function finish_test_dmr_i(dmrr) result(done)
        type(dimer_rep) :: dmrr
        logical :: done
        
        integer :: na
        real(double) :: test
        
        done = .false.
        na = size(dmrr%force_center,2)
        test = sqrt(atom_dot_product(dmrr%eff_force)/(3.0_double*real(na,double)))
        if ((test < dmrr%force_tol) .and. (dmrr%min_curvature < 0.0_double)) done = .true.
        
      end function
        
      function translate_test_i(dmrr,direction) result(done)
        type(dimer_rep) :: dmrr
        real(double), dimension(:,:), intent(in) :: direction
        logical :: done

        integer :: na
        real(double) :: test

        done = .false.
        na = size(dmrr%force_center,2)
        test = atom_dot_product(dmrr%eff_force,direction)/sqrt(3.0_double*real(na,double))
        if ( test < dmrr%force_tol ) done = .true.
 
      end function
           
      subroutine diary_dimer_i(dmrr)
        type(dimer_rep) :: dmrr

        integer :: ia, na

        if ( i_access(diaryfile()) ) then
          write(x_unit(diaryfile()),'(/,"dimer current update, maximum updates: ",i0,", ",i0)') dmrr%n_steps, dmrr%max_steps
          write(x_unit(diaryfile()),'("dimer 01 energy, 02 energy, center energy: ",f0.8,", ",f0.8,", ",f0.8)') &
                                                              & dmrr%cell_energy_01, dmrr%cell_energy_02, dmrr%energy_center
          na = size(dmrr%eff_force,2)
          write(x_unit(diaryfile()),'("dimer effective-force magnitude, force tolerance: ",f0.6,", ",f0.6)') &
                                      & sqrt(atom_dot_product(dmrr%eff_force)/(3.0_double*real(na,double))), dmrr%force_tol
          write(x_unit(diaryfile()),'("dimer center:")')
          write(x_unit(diaryfile()),'(9x,"atom",12x,"x",13x,"y",13x,"z")')
          write(x_unit(diaryfile()),'(8x,"---------------------------------------------------")')
          do ia = 1,size(dmrr%pos_cart_center,2)
            write(x_unit(diaryfile()),'(8x,i3,3x,sp,3(4x,f10.6))') ia, dmrr%pos_cart_center(:,ia)
          end do
          write(x_unit(diaryfile()),'("dimer separation: ",f0.6)') dmrr%separation_mag
          write(x_unit(diaryfile()),'("dimer orientation:")')
          write(x_unit(diaryfile()),'(9x,"atom",12x,"x",13x,"y",13x,"z")')
          write(x_unit(diaryfile()),'(8x,"---------------------------------------------------")')
          do ia = 1,size(dmrr%separation_dir,2)
            write(x_unit(diaryfile()),'(9x,i3,3x,sp,3(4x,f10.6))') ia, dmrr%separation_dir(:,ia)
          end do
          write(x_unit(diaryfile()),'("dimer curvature along orientation: ",f0.6)') dmrr%min_curvature
          write(x_unit(diaryfile()),'("dimer center force magnitude: ",f0.6)') dmrr%force_center_mag
          write(x_unit(diaryfile()),'("dimer center force:")')
          write(x_unit(diaryfile()),'(9x,"atom",11x,"Fx",12x,"Fy",12x,"Fz")')
          write(x_unit(diaryfile()),'(8x,"---------------------------------------------------")')
          do ia = 1,size(dmrr%force_center,2)
            write(x_unit(diaryfile()),'(9x,i3,3x,sp,3(4x,f10.6))') ia, dmrr%force_center(:,ia)
          end do
          write(x_unit(diaryfile()),'("dimer effective force:")')
          write(x_unit(diaryfile()),'(9x,"atom",11x,"Fx",12x,"Fy",12x,"Fz")')
          write(x_unit(diaryfile()),'(8x,"---------------------------------------------------")')
          do ia = 1,size(dmrr%eff_force,2)
            write(x_unit(diaryfile()),'(9x,i3,3x,sp,3(4x,f10.6))') ia, dmrr%eff_force(:,ia)
          end do
          write(x_unit(diaryfile()),'("dimer effective-torque magnitude, torque tolerance: ",f0.6,", ",f0.6)') &
                                      & sqrt(atom_dot_product(dmrr%torque)/(3.0_double*real(na,double))), dmrr%torque_tol
          write(x_unit(diaryfile()),'("dimer torque:")')
          write(x_unit(diaryfile()),'(9x,"atom",11x,"Tx",12x,"Ty",12x,"Tz")')
          write(x_unit(diaryfile()),'(8x,"---------------------------------------------------")')
          do ia = 1,size(dmrr%torque,2)
            write(x_unit(diaryfile()),'(9x,i3,3x,sp,3(4x,f10.6))') ia, dmrr%torque(:,ia)
          end do
        end if

      end subroutine

      subroutine share_dmr_pfe_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        integer :: na
        real(double) :: tmp0
        real(double), dimension(2) :: tmp1
        real(double), dimension(:,:), allocatable :: tmp2
        real(double), dimension(:,:,:), allocatable :: tmp3

        call my(cfg)
        
        na = x_n_atoms(cfg)
        allocate( tmp2(3,na), tmp3(3,na,2) )

        call x_cart_positions(cfg,tmp2)
        call xcomm_allgather(XCONFIG,tmp2,tmp3)
        if (error()) goto 100
        dmrr%pos_cart_01 = tmp3(:,:,1)
        dmrr%pos_cart_02 = tmp3(:,:,2)

        tmp2 = x_forces(cfg)
        call xcomm_allgather(XCONFIG,tmp2,tmp3)
        if (error()) goto 100
        dmrr%force_01 = tmp3(:,:,1)
        dmrr%force_02 = tmp3(:,:,2)

        tmp0 = x_cell_energy(cfg)
        call xcomm_allgather(XCONFIG,tmp0,tmp1)
        if (error()) goto 100
        dmrr%cell_energy_01 = tmp1(1)
        dmrr%cell_energy_02 = tmp1(2)

100     if (allocated( tmp2 )) deallocate( tmp2 )
        if (allocated( tmp3 )) deallocate( tmp3 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::share_dmr_pfe_i")) continue

      end subroutine

      subroutine share_dmr_fe_i(dmrr,cfg)
        type(dimer_rep) :: dmrr
        type(config_sc_obj) :: cfg

        integer :: na
        real(double) :: tmp0
        real(double), dimension(2) :: tmp1
        real(double), dimension(:,:), allocatable :: tmp2
        real(double), dimension(:,:,:), allocatable :: tmp3

        call my(cfg)
        
        na = x_n_atoms(cfg)
        allocate( tmp2(3,na), tmp3(3,na,2) )

        tmp2 = x_forces(cfg)
        call xcomm_allgather(XCONFIG,tmp2,tmp3)
        if (error()) goto 100
        dmrr%force_01 = tmp3(:,:,1)
        dmrr%force_02 = tmp3(:,:,2)

        tmp0 = x_cell_energy(cfg)
        call xcomm_allgather(XCONFIG,tmp0,tmp1)
        if (error()) goto 100
        dmrr%cell_energy_01 = tmp1(1)
        dmrr%cell_energy_02 = tmp1(2)

100     if (allocated( tmp2 )) deallocate( tmp2 )
        if (allocated( tmp3 )) deallocate( tmp3 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::share_dmr_fe_i")) continue

      end subroutine

      subroutine my_neb(neb)
        type(neb_obj) :: neb
        neb%ref = neb%ref + 1
        neb%o%ref = neb%o%ref + 1
      end subroutine

      subroutine my_new_neb(nebi,neb)
      type(neb_obj) :: nebi, neb
        neb%ref = 1
        neb%o => nebi%o
        neb%o%ref = neb%o%ref + 1
      end subroutine

      function thy_neb(neb) result(nebo)
        type(neb_obj) :: neb, nebo
        neb%ref = neb%ref - 1
        neb%o%ref = neb%o%ref - 1
        nebo%ref = neb%ref
        nebo%o => neb%o
      end function

      subroutine glean_neb(neb)
        type(neb_obj) :: neb
        if (neb%o%ref < 1) then
          if (associated( neb%o%pos_cart_end )) deallocate( neb%o%pos_cart_end )
          if (associated( neb%o%pos_cart )) deallocate( neb%o%pos_cart )
          if (associated( neb%o%force )) deallocate( neb%o%force )
          if (associated( neb%o%cell_energy )) deallocate( neb%o%cell_energy )
          deallocate( neb%o )
        end if
      end subroutine

      subroutine bequeath_neb(neb)
        type(neb_obj) :: neb
      end subroutine

      function constructor_neb_i(cfg) result(neb)
        type(config_sc_obj) :: cfg
        type(neb_obj) :: neb
!       effects: Constructs and optimizes an neb.
!       errors: Incorrect tag values. Unrecognized tags. Passes errors.

        logical :: found
        character(line_len) :: tag
        integer :: na, nc
        type(crystal_obj) :: cr_00, cr_nn

        call my(cfg)

        neb%ref = 0
        allocate( neb%o )
        neb%o%g = x_ghost()
        neb%o%ref = 0

        nullify( neb%o%pos_cart_end, neb%o%pos_cart, neb%o%force, neb%o%cell_energy )

        call arg("neb_force_tol",neb%o%force_tol,found)
        if (.not.found) neb%o%force_tol = 1.0e-3_double
        if (error(neb%o%force_tol < 0.0_double,"ERROR: force_tol < 0")) goto 100

        call arg("neb_climbing_image",neb%o%climbing_image,found)
        if (.not.found) neb%o%climbing_image = .true.

        call arg("neb_spring",neb%o%spring,found)
        if (.not.found) neb%o%spring = 1.0e-2_double
        if (error(neb%o%spring < 0.0_double,"ERROR: neb_spring < 0")) goto 100

        call arg("neb_max_steps",neb%o%max_steps,found)
        if (.not.found) neb%o%max_steps = 100
        if (error(neb%o%max_steps < 0,"ERROR: max_steps < 0")) goto 100
        neb%o%n_steps = 0

        call arg("neb_force_test",tag,found)
        if (.not.found) tag = "ALL_IMAGES"
        select case (trim(tag))
        case ("ALL_IMAGES","All_Images","all_images","AI","ai","A","a")
          neb%o%force_test = ALL_IMAGES
        case ("CLIMBING_IMAGE","Climbing_image","climbing_image","CI","ci","C","c")
          neb%o%force_test = CLIMBING_IMAGE
        case default
          if (error(.true.,"ERROR: unrecognized neb_force_test")) goto 100
        end select

        na = x_n_atoms(cfg)
        nc = mpi_nconfigs()

        call my(crystal(neb_crystal_00_path),cr_00)
        call my(crystal(neb_crystal_nn_path),cr_nn)
        allocate( neb%o%pos_cart_end(3,na,2) )
        call x_cart_positions(cr_00,neb%o%pos_cart_end(:,:,1))
        call x_cart_positions(cr_nn,neb%o%pos_cart_end(:,:,2))
        call glean(thy(cr_00))
        call glean(thy(cr_nn))

        allocate( neb%o%pos_cart(3,na,nc), neb%o%force(3,na,nc), neb%o%cell_energy(nc) )
        call share_neb_pfe_i(neb%o,cfg) ; if (error()) goto 100

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'("Completed neb initialization:")')
          if (neb%o%climbing_image) then
            write(x_unit(diaryfile()),'("neb climbing image on")')
          else
            write(x_unit(diaryfile()),'("neb climbing image off")')
          end if
        end if
        if (i_access(output)) then
          write(x_unit(output),'("Completed neb initialization:")')
          if (neb%o%climbing_image) then
            write(x_unit(output),'("neb climbing image on")')
          else
            write(x_unit(output),'("neb climbing image off")')
          end if
        end if

        call arglc("neb_optimizer",tag,found)
        if (.not.found) tag = "qm"
        select case (trim(tag))
        case ("steepest_descent","sd")
          call optimize_neb_sd_i(neb%o,cfg) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'("  Steepest-descent optimization")')
          if (i_access(output)) write(x_unit(output),'("  Steepest-descent optimization")')
        case ("quickmin","qm")
          call optimize_neb_qm_i(neb%o,cfg) ; if (error()) goto 100
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'("  Quickmin optimization")')
          if (i_access(output)) write(x_unit(output),'("  Quickmin optimization")')
        case default
          if (error(.true.,"ERROR: unrecognized neb_optimizer")) goto 100
        end select

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::constructor_neb_i")) continue
        
      end function

      subroutine optimize_neb_sd_i(nebr,cfg) 
        type(neb_rep) :: nebr
        type(config_sc_obj) :: cfg

        logical :: found
        integer :: ic
        real(double) :: sd_factor

        call my(cfg)

        call arg("neb_sd_factor",sd_factor,found)
        if (.not.found) sd_factor = 2.50_double

        do
          call project_force_i(nebr,cfg)
          call diary_neb_i(nebr)
          if (finish_test_neb_i(nebr)) exit
          if (nebr%n_steps == nebr%max_steps) then
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of neb terminated: n_steps = max_steps")')
            if (i_access(output)) write(x_unit(output),'("Optimization of neb terminated: n_steps = max_steps")')
            goto 100
          end if
          call barrier(MOD_SCOPE)
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            goto 100
          end if
          do ic = 1,size(nebr%pos_cart,3)
            nebr%pos_cart(:,:,ic) = sd_factor*nebr%force(:,:,ic) + nebr%pos_cart(:,:,ic)
          end do
          call update_neb_i(nebr,cfg) ; if (error()) goto 100
        end do

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Completed optimization of neb")')
        if (i_access(output)) write(x_unit(output),'("Completed optimization of neb")')

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_neb_sd_i")) continue

      end subroutine
      
      subroutine optimize_neb_qm_i(nebr,cfg) 
        type(neb_rep) :: nebr
        type(config_sc_obj) :: cfg

        logical :: found
        integer :: ic, na
        real(double) :: max_move, rdr, time_step, vdv
        real(double), dimension(:,:), allocatable :: v, dv

        call my(cfg)

        call arg("neb_qm_time_step",time_step,found)
        if (.not.found) time_step = 1.50_double

        call arg("neb_qm_max_move",max_move,found)
        if (.not.found) max_move = 1.0_double

        na = size(nebr%pos_cart,2)
        allocate( v(3,na), dv(3,na) )
        v = 0.0_double

        do
          call project_force_i(nebr,cfg)
          call diary_neb_i(nebr)
          if (finish_test_neb_i(nebr)) exit
          if (nebr%n_steps == nebr%max_steps) then
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of neb terminated: n_steps = max_steps")')
            if (i_access(output)) write(x_unit(output),'("Optimization of neb terminated: n_steps = max_steps")')
            goto 100
          end if
          call barrier(MOD_SCOPE)
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            goto 100
          end if
          do ic = 1,size(nebr%pos_cart,3)
            dv = nebr%force(:,:,ic)*time_step
            where (dv*v < 0.0_double) v = 0.0_double
            vdv = sum(v*dv)
            v = dv*(1.0_double + vdv/sum(dv**2))
            rdr = sqrt(sum((v*time_step)**2))
            if (rdr > max_move) then
              nebr%pos_cart(:,:,ic) = nebr%pos_cart(:,:,ic) + max_move*v/sum(v**2)
            else
              nebr%pos_cart(:,:,ic) = nebr%pos_cart(:,:,ic) + v*time_step
            end if
          end do
          call update_neb_i(nebr,cfg) ; if (error()) goto 100
        end do

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Completed optimization of neb")')
        if (i_access(output)) write(x_unit(output),'("Completed optimization of neb")')

100     if (allocated( v )) deallocate( v )
        if (allocated( dv )) deallocate( dv )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_neb_qm_i")) continue

      end subroutine
      
      subroutine update_neb_i(nebr,cfg)
        type(neb_rep) :: nebr
        type(config_sc_obj) :: cfg

        call my(cfg)
        
        call update_config(cfg,nebr%pos_cart(:,:,mpi_myconfig()))
        call sync_configuration_errors() ; if (error()) goto 100

        call share_neb_fe_i(nebr,cfg) ; if (error()) goto 100

        nebr%n_steps = nebr%n_steps + 1

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::update_neb_i")) continue

      end subroutine

      function finish_test_neb_i(nebr) result(done)
        type(neb_rep) :: nebr
        logical :: done
          
        integer :: ic
        real(double) :: scale, test
          
        scale = 1.0_double/(3.0_double*real(size(nebr%pos_cart,2),double))

        done = .true.

        select case (nebr%force_test)
        case (ALL_IMAGES)
          do ic = 1,size(nebr%force,3)
            test = sqrt(scale*atom_dot_product(nebr%force(:,:,ic)))
            if ( test >= nebr%force_tol ) then
              done = .false.
              exit
            end if
          end do
        case (CLIMBING_IMAGE)
          ic = nebr%max_energy_image
          test = sqrt(scale*atom_dot_product(nebr%force(:,:,ic)))
          if ( test >= nebr%force_tol ) done = .false.
        end select

      end function

      subroutine project_force_i(nebr,cfg)
        type(neb_rep) :: nebr
        type(config_sc_obj) :: cfg

        integer :: ic, na, nc
        logical :: higher_prev, higher_next
        real(double) :: rms_force, scale, umin, umax
        real(double) :: energy_to_next, energy_to_prev
        real(double) :: dist_to_next, dist_to_prev, force_proj
        real(double), dimension(:,:), allocatable :: tangent, vec_diff
        type(lattice_obj) :: lat

        call my(cfg)

        call my(x_lattice(x_crystal(x_external(cfg))),lat)

        na = size(nebr%pos_cart,2)
        nc = size(nebr%pos_cart,3)
        allocate( tangent(3,na), vec_diff(3,na) )

        do ic = 1,nc
 
          if (ic == 1) then                                               ! FIND THE TANGENT TO THE PATH AT EACH IMAGE.
            tangent = nebr%pos_cart(:,:,ic) - nebr%pos_cart_end(:,:,1)
            call reduce_vectors_i(lat,tangent)
          elseif (ic == nc) then
            tangent = nebr%pos_cart_end(:,:,2) - nebr%pos_cart(:,:,ic)
            call reduce_vectors_i(lat,tangent)
          else                                                              ! ARE THE NEIGHBORING IMAGES HIGHER IN ENERGY?
            higher_prev = nebr%cell_energy(ic-1) > nebr%cell_energy(ic)
            higher_next = nebr%cell_energy(ic+1) > nebr%cell_energy(ic)
            if (higher_prev .neqv. higher_next) then                        ! IF AT AN EXTREMA OF ENERGY, THE TANGENT IS THE
              if (higher_prev) then                                         ! VECTOR TO THE LOWER ENERGY NEIGHBOR.
                tangent = nebr%pos_cart(:,:,ic) - nebr%pos_cart(:,:,ic-1)
                call reduce_vectors_i(lat,tangent)
              else
                tangent = nebr%pos_cart(:,:,ic+1) - nebr%pos_cart(:,:,ic)
                call reduce_vectors_i(lat,tangent)
              end if
            else                                                              ! AT AN EXTREMA OF ENERGY, INTERPOLATE THE TANGENT  
              energy_to_prev = nebr%cell_energy(ic-1) - nebr%cell_energy(ic)  ! LINEARLY WITH THE ENERGY.
              energy_to_next = nebr%cell_energy(ic+1) - nebr%cell_energy(ic)
              umin = min(abs(energy_to_prev),abs(energy_to_next))
              umax = max(abs(energy_to_prev),abs(energy_to_next))
              if (energy_to_prev > energy_to_next) then
                tangent = nebr%pos_cart(:,:,ic+1) - nebr%pos_cart(:,:,ic)
                call reduce_vectors_i(lat,tangent)
                tangent = tangent*umin
                vec_diff = nebr%pos_cart(:,:,ic) - nebr%pos_cart(:,:,ic-1)
                call reduce_vectors_i(lat,vec_diff)
                tangent = tangent + vec_diff*umax
              else
                tangent = nebr%pos_cart(:,:,ic+1) - nebr%pos_cart(:,:,ic)
                call reduce_vectors_i(lat,tangent)
                tangent = tangent*umax
                vec_diff = nebr%pos_cart(:,:,ic) - nebr%pos_cart(:,:,ic-1)
                call reduce_vectors_i(lat,vec_diff)
                tangent = tangent + vec_diff*umin
              end if
            end if
          end if
          tangent = tangent/sqrt(sum(tangent**2))

          force_proj = sum(nebr%force(:,:,ic)*tangent)                               ! PROJECT THE FORCE ALONG THE TANGENT, UNLESS
          if (nebr%climbing_image .and. ic == nebr%max_energy_image) then            ! THIS IS THE CLIMBING IMAGE, FOR WHICH THE
            nebr%force(:,:,ic) = nebr%force(:,:,ic) - 2.0_double*tangent*force_proj  ! FORCE ALONG THE TANGENT IS INVERTED.
          else
            nebr%force(:,:,ic) = nebr%force(:,:,ic) - tangent*force_proj
          end if

          if (ic == 1) then                                              ! ADD SPRING FORCES ALONG THE BAND IF THIS IS NOT THE
            vec_diff = nebr%pos_cart(:,:,ic) - nebr%pos_cart_end(:,:,1)  ! CLIMBING IMAGE.
            call reduce_vectors_i(lat,vec_diff)
            dist_to_prev = sqrt(sum(vec_diff**2))
          else
            vec_diff = nebr%pos_cart(:,:,ic) - nebr%pos_cart(:,:,ic-1)
            call reduce_vectors_i(lat,vec_diff)
            dist_to_prev = sqrt(sum(vec_diff**2))
          end if
          if (ic == nc) then
            vec_diff = nebr%pos_cart(:,:,ic) - nebr%pos_cart_end(:,:,2)
            call reduce_vectors_i(lat,vec_diff)
            dist_to_next = sqrt(sum(vec_diff**2))
          else
            vec_diff = nebr%pos_cart(:,:,ic) - nebr%pos_cart(:,:,ic+1)
            call reduce_vectors_i(lat,vec_diff)
            dist_to_next = sqrt(sum(vec_diff**2))
          end if
          if (nebr%climbing_image .and. ic == nebr%max_energy_image) then
            continue
          else
            nebr%force(:,:,ic) = nebr%force(:,:,ic) + nebr%spring*tangent*(dist_to_next - dist_to_prev)
          end if
        end do

        scale = 1.0_double/(3.0_double*real(size(nebr%pos_cart,2),double))  ! FIND THE IMAGE WITH THE LARGEST PROJECTED FORCE
        nebr%max_rms_force = 0.0_double
        nebr%max_rms_force_image = 0
        do ic = 1,size(nebr%force,3)
          rms_force = sqrt(scale*atom_dot_product(nebr%force(:,:,ic)))
          if (ic == mpi_myconfig()) nebr%rms_force = rms_force
          if (rms_force > nebr%max_rms_force) then
            nebr%max_rms_force = rms_force
            nebr%max_rms_force_image = ic
          end if
        end do

        if (allocated( tangent )) deallocate( tangent )
        if (allocated( vec_diff )) deallocate( vec_diff )

        call glean(thy(lat))

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::project_force_i")) continue

      end subroutine

      subroutine diary_neb_i(nebr)
        type(neb_rep) :: nebr

        integer :: ia, ic

        ic = mpi_myconfig()

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,"neb current update, maximum updates: ",i0,", ",i0)') nebr%n_steps, nebr%max_steps
          write(x_unit(diaryfile()),'("neb energies:")')
          write(x_unit(diaryfile()),'("neb image, energy: ",4x,i0,", ",f0.8)') ic, nebr%cell_energy(ic)
          write(x_unit(diaryfile()),'("neb max image, energy: ",i0,", ",f0.8)') nebr%max_energy_image, nebr%max_energy
          write(x_unit(diaryfile()),'("neb RMS forces:")')
          write(x_unit(diaryfile()),'("neb image, force: ",4x,i0,", ",f0.6)') ic, nebr%rms_force
          write(x_unit(diaryfile()),'("neb max image, force: ",i0,", ",f0.6)') nebr%max_rms_force_image, nebr%max_rms_force
          write(x_unit(diaryfile()),'("atomic forces:")')
          write(x_unit(diaryfile()),'(19x,"atom",11x,"Fx",12x,"Fy",12x,"Fz")')
          write(x_unit(diaryfile()),'(18x,"---------------------------------------------------")')
          do ia = 1,size(nebr%force,2)
            write(x_unit(diaryfile()),'(19x,i3,3x,sp,3(4x,f10.6))') ia, nebr%force(:,ia,ic)
          end do
        end if
          
      end subroutine

      subroutine share_neb_pfe_i(nebr,cfg)
        type(neb_rep) :: nebr
        type(config_sc_obj) :: cfg

        integer :: ic, na
        real(double) :: tmp0
        real(double), dimension(:,:), allocatable :: tmp2

        call my(cfg)
        
        na = size(nebr%pos_cart,2)
        allocate( tmp2(3,na) )

        call x_cart_positions(cfg,tmp2)
        call xcomm_allgather(XCONFIG,tmp2,nebr%pos_cart) ; if (error()) goto 100

        tmp2 = x_forces(cfg)
        call xcomm_allgather(XCONFIG,tmp2,nebr%force) ; if (error()) goto 100

        tmp0 = x_cell_energy(cfg)
        call xcomm_allgather(XCONFIG,tmp0,nebr%cell_energy) ; if (error()) goto 100

        nebr%max_energy_image = 1
        nebr%max_energy = nebr%cell_energy(nebr%max_energy_image)
        do ic = 2,size(nebr%cell_energy)
          if (nebr%cell_energy(ic) <= nebr%max_energy) cycle
          nebr%max_energy_image = ic
          nebr%max_energy = nebr%cell_energy(ic)
        end do

100     if (allocated( tmp2 )) deallocate( tmp2 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::share_neb_pfe_i")) continue

      end subroutine

      subroutine share_neb_fe_i(nebr,cfg)
        type(neb_rep) :: nebr
        type(config_sc_obj) :: cfg

        integer :: ic, na
        real(double) :: tmp0
        real(double), dimension(:,:), allocatable :: tmp2

        call my(cfg)
        
        na = x_n_atoms(cfg)
        allocate( tmp2(3,na) )

        tmp2 = x_forces(cfg)
        call xcomm_allgather(XCONFIG,tmp2,nebr%force) ; if (error()) goto 100

        tmp0 = x_cell_energy(cfg)
        call xcomm_allgather(XCONFIG,tmp0,nebr%cell_energy) ; if (error()) goto 100

        nebr%max_energy_image = 1
        nebr%max_energy = nebr%cell_energy(nebr%max_energy_image)
        do ic = 2,size(nebr%cell_energy)
          if (nebr%cell_energy(ic) <= nebr%max_energy) cycle
          nebr%max_energy_image = ic
          nebr%max_energy = nebr%cell_energy(ic)
        end do

100     if (allocated( tmp2 )) deallocate( tmp2 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::share_neb_fe_i")) continue

      end subroutine

      subroutine my_ccts(ccts)
        type(ccts_obj) :: ccts
        ccts%ref = ccts%ref + 1
        ccts%o%ref = ccts%o%ref + 1
      end subroutine

      subroutine my_new_ccts(cctsi,ccts)
      type(ccts_obj) :: cctsi, ccts
        ccts%ref = 1
        ccts%o => cctsi%o
        ccts%o%ref = ccts%o%ref + 1
      end subroutine

      function thy_ccts(ccts) result(cctso)
        type(ccts_obj) :: ccts, cctso
        ccts%ref = ccts%ref - 1
        ccts%o%ref = ccts%o%ref - 1
        cctso%ref = ccts%ref
        cctso%o => ccts%o
      end function

      subroutine glean_ccts(ccts)
        type(ccts_obj) :: ccts
        if (ccts%o%ref < 1) then
          if (associated( ccts%o%pos_cart_01 )) deallocate( ccts%o%pos_cart_01 )
          if (associated( ccts%o%pos_cart_02 )) deallocate( ccts%o%pos_cart_02 )
          if (associated( ccts%o%force_01 )) deallocate( ccts%o%force_01)
          if (associated( ccts%o%force_02 )) deallocate( ccts%o%force_02)
          if (associated( ccts%o%last_pos_cart )) deallocate( ccts%o%last_pos_cart )
          if (associated( ccts%o%force_sum )) deallocate( ccts%o%force_sum )
          if (associated( ccts%o%force_diff )) deallocate( ccts%o%force_diff )
          if (associated( ccts%o%eff_force )) deallocate( ccts%o%eff_force )

          deallocate( ccts%o )
        end if
      end subroutine

      subroutine bequeath_ccts(ccts)
        type(ccts_obj) :: ccts
      end subroutine

      function constructor_ccts(cfg) result(ccts)
        type(config_sc_obj) :: cfg
        type(ccts_obj) :: ccts
!       effects: Constructs and optimizes a ccts (carrier capture transition state finder.

        logical :: found
        character(line_len) :: tag
        integer :: na

        call my(cfg)

        ccts%ref = 0
        allocate( ccts%o )
        ccts%o%g = x_ghost()
        ccts%o%ref = 0

        nullify( ccts%o%pos_cart_01, ccts%o%pos_cart_02,ccts%o%force_01, ccts%o%force_02)
        nullify( ccts%o%last_pos_cart, ccts%o%force_sum, ccts%o%force_diff, ccts%o%eff_force)

        if (error(mpi_nconfigs() /= 2,"ERROR: number of configurations /= 2")) goto 100

        call arg("ccts_force_tol",ccts%o%force_tol,found)
        if (.not.found) ccts%o%force_tol = 0.001_double
        if (error(ccts%o%force_tol < 0.0_double,"ERROR: ccts_force_tol < 0")) goto 100

        call arg("ccts_energy_tol",ccts%o%level_tol,found)
        if (.not.found) ccts%o%level_tol = 0.0001_double
        if (error(ccts%o%level_tol < 0.0_double,"ERROR: ccts_energy_tol < 0")) goto 100

        call arg("ccts_relax_step",ccts%o%step,found)
        if (.not.found) ccts%o%step = 0.1_double
        if (error(ccts%o%step < 0.0_double,"ERROR: ccts_relax_step < 0")) goto 100

        call arg("ccts_newton_scaling",ccts%o%newton_scaling,found)
        if (.not.found) ccts%o%newton_scaling = 1.0_double
        if (error(ccts%o%newton_scaling < 0.0_double,"ERROR: ccts_newton_scaling < 0")) goto 100

        call arg("ccts_max_steps",ccts%o%max_steps,found)
        if (.not.found) ccts%o%max_steps = 100
        if (error(ccts%o%max_steps < 0,"ERROR: ccts_max_steps < 0")) goto 100
        ccts%o%n_steps = 0

        call arg("ccts_move_length",ccts%o%move_length,found)
        if (.not.found) ccts%o%move_length = 0.1_double

        call arg("ccts_max_moves",ccts%o%max_moves,found)
        if (.not.found) ccts%o%max_moves = 10
        if (error(ccts%o%max_moves < 0,"ERROR: ccts_max_moves < 0")) goto 100
        ccts%o%n_moves = 0
        ccts%o%path_length = 0.0_double

        call arg("ccts_move_noise",ccts%o%move_noise,found)
        if (.not.found) ccts%o%move_noise = 0.0_double

        na = x_n_atoms(cfg)
        allocate( ccts%o%pos_cart_01(3,na), ccts%o%pos_cart_02(3,na), ccts%o%force_01(3,na), ccts%o%force_02(3,na) )
        allocate( ccts%o%last_pos_cart(3,na), ccts%o%force_sum(3,na), ccts%o%force_diff(3,na), ccts%o%eff_force(3,na) )

        call initialize_ccts_i(ccts%o,cfg) ; if (error()) goto 100

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Completed ccts initialization")')
        if (i_access(output)) write(x_unit(output),'(/,"Completed ccts initialization")')

!       call diary_ccts_i(ccts%o)
        if ( i_access(diaryfile()) ) then
          write(x_unit(diaryfile()),'(" Initial 01 energy, 02 energy: ",f0.8,", ",f0.8)') ccts%o%cell_energy_01, ccts%o%cell_energy_02
          write(x_unit(diaryfile()),'(" Initial target level set to: ",f0.8,", ",f0.8)') ccts%o%target_level
        end if


        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Beginning optimization of the ccts")')
        if (i_access(output)) write(x_unit(output),'(/,"Beginning optimization of the ccts")')

        call arglc("ccts_optimizer",tag,found)
        if (.not.found) tag = "steepest_descent"
        select case (trim(tag))
        case ("steepest_descent","sd")
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(" Using steepest-descent optimization")')
          if (i_access(output)) write(x_unit(output),'(" Using steepest-descent optimization")')
          call optimize_ccts_sd_i(ccts%o,cfg) ; if (error()) goto 100
        case default
          if (error(.true.,"ERROR: unrecognized ccts_optimizer")) goto 100
        end select

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Completed optimization of ccts")')
        if (i_access(output)) write(x_unit(output),'(/,"Completed optimization of ccts")')

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::constructor_ccts")) continue
        
      end function

      subroutine initialize_ccts_i(cctsr,cfg)
        type(ccts_rep) :: cctsr
        type(config_sc_obj) :: cfg

        real(double) :: pos_diff_sq_norm
        integer :: na

        call my(cfg)

        call share_ccts_pfe_i(cctsr,cfg) ; if (error()) goto 100

        na = size(cctsr%eff_force,2)

        pos_diff_sq_norm = sum((cctsr%pos_cart_01 - cctsr%pos_cart_02)**2)
        if (error(pos_diff_sq_norm.gt.machine_precision,"ERROR: 01/crystal and 02/crystal must be the same")) goto 100

        cctsr%pos_cart_02 = cctsr%pos_cart_01

        cctsr%force_sum = cctsr%force_01 + cctsr%force_02
        cctsr%force_sum_mag_sq = atom_dot_product(cctsr%force_sum)

        cctsr%force_diff = cctsr%force_01 - cctsr%force_02
        cctsr%force_diff_mag_sq = atom_dot_product(cctsr%force_diff)

        cctsr%eff_force = cctsr%force_01  &
                - cctsr%force_diff * atom_dot_product(cctsr%force_diff,cctsr%force_01)/cctsr%force_diff_mag_sq
        cctsr%eff_force_rms = sqrt(atom_dot_product(cctsr%eff_force)/(3.0_double*real(na,double)))

        cctsr%level = cctsr%cell_energy_02 - cctsr%cell_energy_01
        cctsr%target_level = cctsr%level
        cctsr%level_err = cctsr%level - cctsr%target_level

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::initialize_ccts_i")) continue

      end subroutine

      subroutine optimize_ccts_sd_i(cctsr,cfg) 
        type(ccts_rep) :: cctsr
        type(config_sc_obj) :: cfg

        integer :: na
        real(double), dimension(:,:), allocatable :: move

        na = size(cctsr%eff_force,2)
        allocate( move(3,na))

        call my(cfg)

        do while (cctsr%n_moves .le. cctsr%max_moves)

           if (i_access(diaryfile())) write(x_unit(diaryfile()),"(/,'Begin constrained minimization')")
           if (i_access(output)) write(x_unit(output),'("/,Begin constrained minimization")')

!          Add a random component to the coordinates to break any symmetries
           if (cctsr%move_noise > 0.0_double) then

             if (i_access(diaryfile())) write(x_unit(diaryfile()),"(/,'Perturbing coordinates to break symmetries')")
             if (i_access(output)) write(x_unit(output),'("/,Perturbing coordinates to break symmetries")')

             call random_number(move)
             move = cctsr%move_noise*(2.0_double*move - 1.0_double)
             move = move - cctsr%force_diff * atom_dot_product(cctsr%force_diff,move)/cctsr%force_diff_mag_sq
             call broadcast(WORLD,move) ; if (error()) goto 100

             if (i_access(diaryfile())) then
               write(x_unit(diaryfile()),"(' move the ccts by:')")
               write(x_unit(diaryfile()),"(2x,3g15.5)") move
             end if
             if (i_access(output)) then
               write(x_unit(output),"(' move the ccts by:')")
               write(x_unit(output),"(2x,3g15.5)") move
             end if

             cctsr%pos_cart_01 = cctsr%pos_cart_01 + move
             cctsr%pos_cart_02 = cctsr%pos_cart_01

             call update_ccts_i(cctsr,cfg) ; if (error()) goto 100

             if ( i_access(diaryfile()) ) then
               write(x_unit(diaryfile()),'(/," 01 energy, 02 energy: ",f0.8,", ",f0.8)') cctsr%cell_energy_01, cctsr%cell_energy_02
               write(x_unit(diaryfile()),'(" level, target level: ",f0.8,", ",f0.8)') &
                                                              & cctsr%level, cctsr%target_level
               write(x_unit(diaryfile()),'(" level error, level error tolerance: ",f0.8,", ",f0.8)') &
                                                              & abs(cctsr%level_err), cctsr%level_tol
               write(x_unit(diaryfile()),'(" rms effective-force, force tolerance: ",f0.6,", ",f0.6)') &
                                      & cctsr%eff_force_rms, cctsr%force_tol
             end if

           end if

           do while ( (cctsr%eff_force_rms .gt. cctsr%force_tol) .or. &
                      (abs(cctsr%level_err) .gt. cctsr%level_tol) )
   
             call barrier(MOD_SCOPE)
             if (user_stop()) then
               call warn("WARNING: USER INITIATED STOP")
               goto 100
             end if
   
             if (i_access(diaryfile())) then
                write(x_unit(diaryfile()),"(/,'Constrained minimization step')")
                write(x_unit(diaryfile()),'(" current step, maximum steps: ",i0,", ",i0)') cctsr%n_steps, cctsr%max_steps
             end if
             if (i_access(output)) write(x_unit(output),'("/,Constrained minimization step")')
   
             move = cctsr%step*cctsr%eff_force - cctsr%newton_scaling*cctsr%level_err*cctsr%force_diff/cctsr%force_diff_mag_sq

             if (i_access(diaryfile())) then
               write(x_unit(diaryfile()),"(' move the ccts by:')")
               write(x_unit(diaryfile()),"(2x,3g15.5)") move
             end if
             if (i_access(output)) then
               write(x_unit(output),"(' move the ccts by:')")
               write(x_unit(output),"(2x,3g15.5)") move
             end if
   
             cctsr%pos_cart_01 = cctsr%pos_cart_01 + move
             cctsr%pos_cart_02 = cctsr%pos_cart_01
   
             call update_ccts_i(cctsr,cfg) ; if (error()) goto 100
   
             if ( i_access(diaryfile()) ) then
               write(x_unit(diaryfile()),'(/," 01 energy, 02 energy: ",f0.8,", ",f0.8)') cctsr%cell_energy_01, cctsr%cell_energy_02
               write(x_unit(diaryfile()),'(" level, target level: ",f0.8,", ",f0.8)') &
                                                              & cctsr%level, cctsr%target_level
               write(x_unit(diaryfile()),'(" level error, level error tolerance: ",f0.8,", ",f0.8)') &
                                                              & abs(cctsr%level_err), cctsr%level_tol
               write(x_unit(diaryfile()),'(" rms effective-force, force tolerance: ",f0.6,", ",f0.6)') &
                                      & cctsr%eff_force_rms, cctsr%force_tol
             end if
             
             if (cctsr%n_steps > cctsr%max_steps) then
               if (i_access(diaryfile())) write(x_unit(diaryfile()),'("Optimization of ccts terminated: n_steps = max_steps")')
               if (i_access(output)) write(x_unit(output),'("Optimization of ccts terminated: n_steps = max_steps")')
               goto 100
             end if
   
           end do

!          call diary_ccts_i(cctsr)

           if (cctsr%n_moves /= 0) then
              move = cctsr%pos_cart_01 - cctsr%last_pos_cart
              cctsr%path_length = cctsr%path_length + sqrt(atom_dot_product(move))
           end if

           if ( i_access(diaryfile()) ) then
              write(x_unit(diaryfile()),'(/,"Constrained minimization completed")')
              write(x_unit(diaryfile()),'(" path length, 01 energy, 02 energy, level: ",f0.8,", ",f0.8,", ",f0.8,", ",f0.8)') &
                                      cctsr%path_length, cctsr%cell_energy_01, cctsr%cell_energy_02, cctsr%level
           end if
           if (i_access(output)) write(x_unit(output),'(/,"Constrained minimization completed")')

           cctsr%n_moves = cctsr%n_moves + 1

           if (i_access(diaryfile())) then
              write(x_unit(diaryfile()),"(/,'Moving along the reaction path')")
              write(x_unit(diaryfile()),'(" current move, maximum moves: ",i0,", ",i0)') cctsr%n_moves, cctsr%max_moves
           end if
           if (i_access(output)) write(x_unit(output),'(/,"Moving along the reaction path")') 

           if (cctsr%n_moves == 1) then
              move = cctsr%move_length*cctsr%force_diff / sqrt(cctsr%force_diff_mag_sq)
           else
              move = cctsr%pos_cart_01 - cctsr%last_pos_cart
              move = move - cctsr%force_diff * atom_dot_product(cctsr%force_diff,move)/cctsr%force_diff_mag_sq
              move = move + cctsr%move_length * cctsr%force_diff / sqrt(cctsr%force_diff_mag_sq)
           end if
           move = abs(cctsr%move_length) * move / sqrt(atom_dot_product(move))
!          move = cctsr%move_length * cctsr%force_diff / sqrt(cctsr%force_diff_mag_sq)


           if (i_access(diaryfile())) then
              write(x_unit(diaryfile()),"(' move the ccts by:')")
              write(x_unit(diaryfile()),"(2x,3g15.5)") move
           end if
           if (i_access(output)) then
              write(x_unit(output),"(' move the ccts by:')")
              write(x_unit(output),"(2x,3g15.5)") move
           end if

           cctsr%last_pos_cart = cctsr%pos_cart_01
           cctsr%pos_cart_01 = cctsr%pos_cart_01 + move
           cctsr%pos_cart_02 = cctsr%pos_cart_01

           cctsr%n_steps = 0
           call update_ccts_i(cctsr,cfg) ; if (error()) goto 100

           cctsr%target_level = cctsr%level
           cctsr%level_err = cctsr%level - cctsr%target_level

            if ( i_access(diaryfile()) ) then
               write(x_unit(diaryfile()),'(/," new 01 energy, 02 energy: ",f0.8,", ",f0.8)') cctsr%cell_energy_01, cctsr%cell_energy_02
               write(x_unit(diaryfile()),'(" new target level: ",f0.8)') cctsr%target_level
             end if

        end do
   
100     if (allocated( move )) deallocate( move )
        call glean(thy(cfg))

        if (error("Exit transition_state_mod::optimize_ccts_sd_i")) continue

      end subroutine

      subroutine update_ccts_i(cctsr,cfg)
        type(ccts_rep) :: cctsr
        type(config_sc_obj) :: cfg

        integer :: na

        call my(cfg)

        na = size(cctsr%eff_force,2)
        
        select case (mpi_myconfig())
        case (1)
          call update_config(cfg,cctsr%pos_cart_01)
        case (2)
          call update_config(cfg,cctsr%pos_cart_02)
        end select
        call sync_configuration_errors() ; if (error()) goto 100

        call share_ccts_fe_i(cctsr,cfg) ; if (error()) goto 100

        cctsr%force_sum = cctsr%force_01 + cctsr%force_02
        cctsr%force_sum_mag_sq = atom_dot_product(cctsr%force_sum)

        cctsr%force_diff = cctsr%force_01 - cctsr%force_02
        cctsr%force_diff_mag_sq = atom_dot_product(cctsr%force_diff)

        cctsr%eff_force = cctsr%force_01  &
                - cctsr%force_diff * atom_dot_product(cctsr%force_diff,cctsr%force_01)/cctsr%force_diff_mag_sq
        cctsr%eff_force_rms = sqrt(atom_dot_product(cctsr%eff_force)/(3.0_double*real(na,double)))

        cctsr%level = cctsr%cell_energy_02 - cctsr%cell_energy_01
        cctsr%level_err = cctsr%level - cctsr%target_level

        cctsr%n_steps = cctsr%n_steps + 1

100     call glean(thy(cfg))

        if (error("Exit transition_state_mod::update_ccts_i")) continue

      end subroutine

      subroutine diary_ccts_i(cctsr)
        type(ccts_rep) :: cctsr

        integer :: ia

        if ( i_access(diaryfile()) ) then
          write(x_unit(diaryfile()),'("ccts 01 energy, 02 energy: ",f0.8,", ",f0.8)') cctsr%cell_energy_01, cctsr%cell_energy_02
          write(x_unit(diaryfile()),'("ccts level, target level: ",f0.8,", ",f0.8)') &
                                                              & cctsr%level, cctsr%target_level
        end if

      end subroutine

      subroutine share_ccts_pfe_i(cctsr,cfg)
        type(ccts_rep) :: cctsr
        type(config_sc_obj) :: cfg

        integer :: na
        real(double) :: tmp0
        real(double), dimension(2) :: tmp1
        real(double), dimension(:,:), allocatable :: tmp2
        real(double), dimension(:,:,:), allocatable :: tmp3

        call my(cfg)
        
        na = x_n_atoms(cfg)
        allocate( tmp2(3,na), tmp3(3,na,2) )

        call x_cart_positions(cfg,tmp2)
        call xcomm_allgather(XCONFIG,tmp2,tmp3)
        if (error()) goto 100
        cctsr%pos_cart_01 = tmp3(:,:,1)
        cctsr%pos_cart_02 = tmp3(:,:,2)

        tmp2 = x_forces(cfg)
        call xcomm_allgather(XCONFIG,tmp2,tmp3)
        if (error()) goto 100
        cctsr%force_01 = tmp3(:,:,1)
        cctsr%force_02 = tmp3(:,:,2)

        tmp0 = x_cell_energy(cfg)
        call xcomm_allgather(XCONFIG,tmp0,tmp1)
        if (error()) goto 100
        cctsr%cell_energy_01 = tmp1(1)
        cctsr%cell_energy_02 = tmp1(2)

100     if (allocated( tmp2 )) deallocate( tmp2 )
        if (allocated( tmp3 )) deallocate( tmp3 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::share_ccts_pfe_i")) continue

      end subroutine

      subroutine share_ccts_fe_i(cctsr,cfg)
        type(ccts_rep) :: cctsr
        type(config_sc_obj) :: cfg

        integer :: na
        real(double) :: tmp0
        real(double), dimension(2) :: tmp1
        real(double), dimension(:,:), allocatable :: tmp2
        real(double), dimension(:,:,:), allocatable :: tmp3

        call my(cfg)
        
        na = x_n_atoms(cfg)
        allocate( tmp2(3,na), tmp3(3,na,2) )

        tmp2 = x_forces(cfg)
        call xcomm_allgather(XCONFIG,tmp2,tmp3)
        if (error()) goto 100
        cctsr%force_01 = tmp3(:,:,1)
        cctsr%force_02 = tmp3(:,:,2)

        tmp0 = x_cell_energy(cfg)
        call xcomm_allgather(XCONFIG,tmp0,tmp1)
        if (error()) goto 100
        cctsr%cell_energy_01 = tmp1(1)
        cctsr%cell_energy_02 = tmp1(2)

100     if (allocated( tmp2 )) deallocate( tmp2 )
        if (allocated( tmp3 )) deallocate( tmp3 )

        call glean(thy(cfg))

        if (error("Exit transition_state_mod::share_ccts_fe_i")) continue

      end subroutine

      subroutine reduce_vectors_i(lat,vecx)
        type(lattice_obj) :: lat
        real(double), dimension(:,:), intent(inout) :: vecx
!       requires: vecx be in the cartesian representation.
!       effects: Modifies vecx such that all vectors end within the Wigner-Seitz cell.
!       errors: Passes errors.

        logical :: inside
        integer :: iv, iws, m
        real(double), parameter :: tol_nbhd = 1.0e-9_double
        real(double) :: dp
        real(double), dimension(3) :: vecl
        real(double), dimension(:,:), pointer :: ws
        real(double), dimension(:,:), allocatable :: wsx

        call my(lat)

        m = 1
        nullify( ws )
        call wigner_seitz_vectors(lat,m,ws)

        allocate( wsx(size(ws,1),size(ws,2)) )
        do iws = 1,size(ws,2)
          wsx(:,iws) = lat2r(lat,ws(:,iws))
          wsx(:,iws) = wsx(:,iws)/dot_product(wsx(:,iws),wsx(:,iws))
        end do

        do iv = 1,size(vecx,2)
          vecl = r2lat(lat,vecx(:,iv))
          inside = .false.
          do while (.not.inside)
            do iws = 1,size(ws,2)
              dp = dot_product(vecx(:,iv),wsx(:,iws))
              if ((dp < 0.5_double) .or. (dp .in. nbhd(0.5_double,tol_nbhd))) then
                inside = .true.
              else
                inside = .false.
                vecl = vecl - ws(:,iws)
                vecx(:,iv) = lat2r(lat,vecl)
                exit
              end if
            end do
          end do
        end do

100     if (associated( ws )) deallocate( ws )
        if (allocated( wsx )) deallocate( wsx )

        call glean(thy(lat))

        if (error("Exit transition_state_mod::reduce_vectors_i")) continue

      end subroutine
 
    end module

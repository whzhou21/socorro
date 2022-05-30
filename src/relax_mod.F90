! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module relax_mod
!doc$ module relax_mod

!     relax_mod defines relax_obj which defines the type of relaxation
!     performed for lattice vectors, as well as, a method and related
!     parameters used to optimize the atomic positions in a crystal.
!     The current choices for lattice vector relxation are
!        1. none
!        2. cell shape optimization
!        3. cell shape optimization with quenched motion
!        4. cell shape optimization along single axis(z) with quenched motion
!     Three methods are currently implemented to optimize atomic positions:
!        1. steepest descent
!        2. conjugate gradient
!        3. quench minimization

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use diary_mod
      use arg_mod
      use ghost_mod
      use config_sc_mod
      use external_mod
      use crystal_mod
      use lattice_mod
      use shortcuts_mod
      use math_mod
      use interrupt_mod
      use atomic_operators_mod
      use layout_mod
      use atoms_mod
      use symmetry_mod
      use xc_type_mod

!cod$
      implicit none
      private

      integer, parameter :: NONE = 0
      integer, parameter :: STEEPEST__DESCENT = 1
      integer, parameter :: CONJUGATE__GRADIENT = 2
      integer, parameter :: QUENCH__MINIMIZATION = 3
      integer, parameter :: SHAPE_LINEAR = 4
      integer, parameter :: SHAPE_QUENCH = 5
      integer, parameter :: SHAPE_QUENCH_Z = 6

      real(double), parameter :: AMU_2_ELECTRON_MASS = 9.1144424e+02_double
      real(double) :: stress_tol

      type :: relax_rep
         integer :: ref
         type(ghost) :: g
         integer :: lattice_relaxation, method, max_steps, step
         real(double) :: force_nrm, force_tol
         real(double) :: sd_con
         real(double) :: time_step
      end type

      type, public :: relax_obj
         private
         integer :: ref
         type(relax_rep), pointer :: o
      end type

!doc$
      public :: relax
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: optimize_crystal
      public :: optimize_positions
      public :: steepest_descent
      public :: conjugate_gradient
      public :: quench_minimization
      public :: diary

!cod$
      interface relax
         module procedure constructor_rx_1, constructor_rx_2
      end interface
      interface update
         module procedure update_relax
      end interface
      interface my
         module procedure my_relax, my_new_relax
      end interface
      interface thy
         module procedure thy_relax
      end interface
      interface glean
         module procedure glean_relax
      end interface
      interface bequeath
         module procedure bequeath_relax
      end interface
      interface assignment(=)
         module procedure assign_relax
      end interface
      interface optimize_crystal
         module procedure relax_structure, relax_structure_rx
      end interface
      interface optimize_positions
         module procedure relax_brancher, relax_brancher_rx
      end interface
      interface steepest_descent
         module procedure relax_sd
      end interface
      interface conjugate_gradient
         module procedure relax_cg
      end interface
      interface quench_minimization
         module procedure relax_qm
      end interface
      interface diary
         module procedure diary_relax
      end interface
      
      contains

      function constructor_rx_1() result(rx)
!doc$ function relax() result(rx)
        type(relax_obj) :: rx
!       effects: Constructs rx using a tag and parameters read from an arguments file.
!       errors: Unrecognized tag or out-of-range parameters.
!       arguments:
!         lattice_relaxation: none, shape, shape_quench
!         relax_method: none, steepest_descent, conjugate_gradient, quench_minimization
!         For relax_method = steepest_descent, conjugate_gradient, or quench_minimization:
!            relax_steps: > 0  (default = 100)
!            relax_tol: > 0  (default = 1e^-3)
!         For relax_method = steepest_descent:
!            relax_prefactor: > 0  (default = 1)
!         For relax_method = quench_minimization:
!            relax_time_step: > 0  (default = 1)

!cod$
        character(line_len) :: tag
        logical :: found

        rx%ref = 0
        allocate( rx%o )
        rx%o%ref = 0
        rx%o%g = x_ghost()

        call arg("lattice_relaxation",tag,found)
        if (.not.found) tag = "NONE"
        select case (tag(1:len_trim(tag)))
        case ( "none", "NONE", "None" )
          rx%o%lattice_relaxation = NONE
        case ( "shape", "SHAPE", "Shape" )
          rx%o%lattice_relaxation = SHAPE_LINEAR
        case ( "shape_quench", "SHAPE_QUENCH", "Shape_quench" )
          rx%o%lattice_relaxation = SHAPE_QUENCH
        case ( "shape_quench_z", "SHAPE_QUENCH_Z", "Shape_quench_z")
           rx%o%lattice_relaxation = SHAPE_QUENCH_Z
        case default
          if (error(.true.,"ERROR: unrecognized lattice_relaxation")) goto 100
        end select

        select case (rx%o%lattice_relaxation)
        case (SHAPE_LINEAR, SHAPE_QUENCH, SHAPE_QUENCH_Z)
           call arg("lattice_relax_tol",stress_tol,found)
           if (.not.found) stress_tol = 1.0_double
           if (error(stress_tol < 0.0_double, "ERROR: lattice_relax_tol < 0")) goto 100
        end select   

        call arg("relax_method",tag,found)
        if (.not.found) tag = "NONE"
        select case (tag(1:len_trim(tag)))
        case ( "none", "NONE", "None" )
          rx%o%method = NONE
        case ( "steepest_descent", "STEEPEST_DESCENT", "Steepest_descent", "sd", "SD" )
          rx%o%method = STEEPEST__DESCENT
        case ( "conjugate_gradient", "CONJUGATE_GRADIENT", "Conjugate_gradient", "cg", "CG" )
          rx%o%method = CONJUGATE__GRADIENT
        case ( "quench_minimization", "QUENCH_MINIMIZATION", "Quench_minimization", "qm", "QM" )
          rx%o%method = QUENCH__MINIMIZATION
        case default
          if (error(.true.,"ERROR: unrecognized relax_method")) goto 100
        end select

        select case (rx%o%method)
        case (STEEPEST__DESCENT,CONJUGATE__GRADIENT,QUENCH__MINIMIZATION)
          call arg("relax_steps",rx%o%max_steps,found)
          if (.not.found) rx%o%max_steps = 100
          if (error(rx%o%max_steps < 0,"ERROR: relax_steps < 0")) goto 100
          call arg("relax_tol",rx%o%force_tol,found)
          if (.not.found) rx%o%force_tol = 1.0e-3_double
          if (error(rx%o%force_tol < 0.0_double,"ERROR: relax_tol < 0")) goto 100
          rx%o%step = 0
          rx%o%force_nrm = 1.0e10_double
        end select

        select case (rx%o%method)
        case (STEEPEST__DESCENT)
          call arg("relax_prefactor",rx%o%sd_con,found)
          if (.not.found) rx%o%sd_con = 1.0_double
          if (error(rx%o%sd_con <= 0.0_double,"ERROR: relax_prefactor <= 0")) goto 100
        case (QUENCH__MINIMIZATION)
          call arg("relax_time_step",rx%o%time_step,found)
          if (.not.found) rx%o%time_step = 100.0_double
          if (error(rx%o%time_step <= 0.0_double,"ERROR: relax_time_step <= 0")) goto 100
        end select

100     if (error("Exit relax_mod::constructor_rx_1")) continue

      end function
      
      function constructor_rx_2(lattice_relaxation,method,max_steps,force_tol,sd_con,time_step) result(rx)
!doc$ function relax(lattice_relaxation,method,max_steps,force_tol,sd_con,time_step) result(rx)
        character(line_len), intent(in) :: lattice_relaxation
        character(line_len), intent(in) :: method
        integer, intent(in), optional :: max_steps
        real(double), intent(in), optional :: force_tol, sd_con, time_step
        type(relax_obj) :: rx
!       effects: Constructs rx using method and optional parameters.
!       errors: Unrecognized method or out-of-range parameters.
        logical :: found

!cod$
        rx%ref = 0
        allocate( rx%o )
        rx%o%ref = 0
        rx%o%g = x_ghost()

        select case (lattice_relaxation(1:len_trim(lattice_relaxation)))
        case ( "none", "NONE", "None" )
          rx%o%lattice_relaxation = NONE
        case ( "shape", "SHAPE", "Shape" )
          rx%o%lattice_relaxation = SHAPE_LINEAR
        case ( "shape_quench", "SHAPE_QUENCH", "Shape_quench" )
          rx%o%lattice_relaxation = SHAPE_QUENCH
        case ( "shape_quench_z", "SHAPE_QUENCH_Z", "Shape_quench_z")
           rx%o%lattice_relaxation = SHAPE_QUENCH_Z
        case default
          if (error(.true.,"ERROR: unrecognized lattice_relaxation")) goto 100
        end select

        select case (rx%o%lattice_relaxation)
        case (SHAPE_LINEAR, SHAPE_QUENCH, SHAPE_QUENCH_Z)
           call arg("lattice_relax_tol",stress_tol,found)
           if (.not.found) stress_tol = 1.0_double
           if (error(stress_tol < 0.0_double, "ERROR: lattice_relax_tol < 0")) goto 100
        end select

        select case (method(1:len_trim(method)))
        case ( "none", "NONE", "None" )
          rx%o%method = NONE
        case ( "steepest_descent", "STEEPEST_DESCENT", "Steepest_descent", "sd", "SD" )
          rx%o%method = STEEPEST__DESCENT
        case ( "conjugate_gradient", "CONJUGATE_GRADIENT", "Conjugate_gradient", "cg", "CG" )
          rx%o%method = CONJUGATE__GRADIENT
        case ( "quench_minimization", "QUENCH_MINIMIZATION", "Quench_minimization", "qm", "QM" )
          rx%o%method = QUENCH__MINIMIZATION
        case default
          if (error(.true.,"ERROR: unrecognized method")) goto 100
        end select

        select case (rx%o%method)
        case (STEEPEST__DESCENT,CONJUGATE__GRADIENT,QUENCH__MINIMIZATION)
          rx%o%max_steps = 100
          if (present(max_steps)) rx%o%max_steps = max_steps
          if (error(rx%o%max_steps < 1,"ERROR: max_steps < 1")) goto 100
          rx%o%force_tol = 1.0e-3_double
          if (present(force_tol)) rx%o%force_tol = force_tol
          if (error(rx%o%force_tol < 0.0_double,"ERROR: force_tol < 0")) goto 100
          rx%o%step = 1
          rx%o%force_nrm = 1.0e10_double
        end select

        select case (rx%o%method)
        case (STEEPEST__DESCENT)
          rx%o%sd_con = 1.0_double
          if (present(sd_con)) rx%o%sd_con = sd_con
          if (error(rx%o%sd_con <= 0.0_double,"ERROR: sd_con <= 0")) goto 100
        case (QUENCH__MINIMIZATION)
          rx%o%time_step = 1.0_double
          if (present(time_step)) rx%o%time_step = time_step
          if (error(rx%o%time_step <= 0.0_double,"ERROR: time_step <= 0")) goto 100
        end select

100     if (error("Exit relax_mod::constructor_rx_2")) continue

      end function
      
      subroutine update_relax(rx,lattice_relaxation,method,max_steps,force_tol,sd_con,time_step)
!doc$ subroutine update(rx,lattice_relaxation,method,max_steps,force_tol,sd_con,time_step)
        type(relax_obj), intent(inout) :: rx
        character(line_len), intent(in), optional :: lattice_relaxation
        character(line_len), intent(in), optional :: method
        integer, intent(in), optional :: max_steps
        real(double), intent(in), optional :: force_tol, sd_con, time_step
!       modifies: rx
!       effects: Updates rx. Resets rx%o%step and rx%o%force_nrm.
!       errors: Unrecognized lattice_relaxation or method or out-of-range parameters.

!cod$
        logical :: lattice_relaxation_change, method_change, found
        integer :: old_lattice_relaxation, old_method

        call my(rx)

        call own_i(rx)
        rx%o%g = x_ghost()

        old_lattice_relaxation = rx%o%lattice_relaxation
        if (present(lattice_relaxation)) then
          select case (lattice_relaxation(1:len_trim(lattice_relaxation)))
          case ( "none", "NONE", "None" )
             rx%o%lattice_relaxation = NONE
          case ( "shape", "SHAPE", "Shape" )
             rx%o%lattice_relaxation = SHAPE_LINEAR
          case ( "shape_quench", "SHAPE_QUENCH", "Shape_quench" )
             rx%o%lattice_relaxation = SHAPE_QUENCH
        case ( "shape_quench_z", "SHAPE_QUENCH_Z", "Shape_quench_z")
           rx%o%lattice_relaxation = SHAPE_QUENCH_Z
          case default
            if (error(.true.,"ERROR: unrecognized lattice_relaxation")) goto 100
          end select
        end if
        lattice_relaxation_change = ( rx%o%lattice_relaxation /= old_lattice_relaxation )

        select case (rx%o%lattice_relaxation)
        case (SHAPE_LINEAR, SHAPE_QUENCH, SHAPE_QUENCH_Z)
           call arg("lattice_relax_tol",stress_tol,found)
           if (.not.found) stress_tol = 1.0_double
           if (error(stress_tol < 0.0_double, "ERROR: lattice_relax_tol < 0")) goto 100
        end select

        old_method = rx%o%method
        if (present(method)) then
          select case (method(1:len_trim(method)))
          case ( "none", "NONE", "None" )
            rx%o%method = NONE
          case ( "steepest_descent", "STEEPEST_DESCENT", "Steepest_descent", "sd", "SD" )
            rx%o%method = STEEPEST__DESCENT
          case ( "conjugate_gradient", "CONJUGATE_GRADIENT", "Conjugate_gradient", "cg", "CG" )
            rx%o%method = CONJUGATE__GRADIENT
          case ( "quench_minimization", "QUENCH_MINIMIZATION", "Quench_minimization", "qm", "QM" )
            rx%o%method = QUENCH__MINIMIZATION
          case default
            if (error(.true.,"ERROR: unrecognized method")) goto 100
          end select
        end if
        method_change = ( rx%o%method /= old_method )

        select case (rx%o%method)
        case (STEEPEST__DESCENT,CONJUGATE__GRADIENT,QUENCH__MINIMIZATION)
          if (method_change) rx%o%max_steps = 100
          if (present(max_steps)) rx%o%max_steps = max_steps
          if (error(rx%o%max_steps < 1,"ERROR: max_steps < 1")) goto 100
          if (method_change) rx%o%force_tol = 1.0e-3_double
          if (present(force_tol)) rx%o%force_tol = force_tol
          if (error(rx%o%force_tol < 0.0_double,"ERROR: force_tol < 0")) goto 100
          rx%o%step = 1
          rx%o%force_nrm = 1.0e10_double
        end select

        select case (rx%o%method)
        case (STEEPEST__DESCENT)
          if (method_change) rx%o%sd_con = 1.0_double
          if (present(sd_con)) rx%o%sd_con = sd_con
          if (error(rx%o%sd_con <= 0.0_double,"ERROR: sd_con <= 0")) goto 100
        case (QUENCH__MINIMIZATION)
          if (method_change) rx%o%time_step = 1.0_double
          if (present(time_step)) rx%o%time_step = time_step
          if (error(rx%o%time_step <= 0.0_double,"ERROR: time_step <= 0")) goto 100
        end select

100     call glean(thy(rx))

        if (error("Exit relax_mod::update_relax")) continue

      end subroutine

      function relax_structure(cfg) result(changed)
!doc$ function optimize_crystal(cfg) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        logical :: changed
!       modifies: cfg
!       effects: Constructs rx and uses it to relax atom positions and/or
!          lattice vectors in cfg.
!       errors: Passes errors.

!cod$
        type(relax_obj) :: rx

        call my(cfg)
        call my(relax(),rx) ; if (error()) goto 100

        changed = optimize_crystal(cfg,rx)

100     call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::relax_structure")) continue

      end function

      function relax_structure_rx(cfg,rx) result(changed)
!doc$ function optimize_crystal(cfg,rx) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        logical :: changed
!       modifies: cfg and rx
!       effects: Uses rx to relax atom positions and/or lattice vectors in cfg.
!       errors: Passes errors.

!cod$
        real(double), dimension(3,3) :: s, v, step
        type(external_obj) :: ext
        type(crystal_obj) :: cr
        type(lattice_obj) :: lat
        real(double) :: latc, check
        integer :: stress_step, max_steps
        real(double) :: stress_norm
        integer, parameter :: max_stress_steps = 50
        real(double), parameter :: au_2_kbar = 147105.164_double
        real(double), parameter :: lattice_prefactor = 1.0e-4_double
        integer :: i, j, k

        call my(cfg)
        call my(rx)

        call my(x_external(cfg),ext)  ; if (error()) goto 100
        call my(x_crystal(ext),cr)  ; if (error()) goto 200
        call my(x_lattice(cr),lat)  ; if (error()) goto 300

        SELECT CASE (rx%o%lattice_relaxation)
           CASE (NONE)
              changed = optimize_positions(cfg,rx) ; if (error()) goto 400
          CASE (SHAPE_QUENCH_Z)
             !j sets the columns of the stress tensor we work with.  3 corresponds to z-axis values
              j = 3

              max_steps = rx%o%max_steps
              rx%o%max_steps = 5
              changed = optimize_positions(cfg,rx) ; if (error()) goto 400

              step = 0.0_double
              v(1,:) = x_lattice_vector(lat,1)  ; if (error()) goto 400
              v(2,:) = x_lattice_vector(lat,2)  ; if (error()) goto 400
              v(3,:) = x_lattice_vector(lat,3)  ; if (error()) goto 400
              latc = x_lattice_constant(lat)  ; if (error()) goto 400

              stress_step = 0
              s = x_stress_tensor(cfg) ; if (error()) goto 400
              call diary_stress_tensor(cfg) ; if (error()) goto 400
              !calculate stress contributions that effect the z-axis
              stress_norm = au_2_kbar*sqrt(sum(s(:,3)**2)/3.0_double)
              if (i_access( diaryfile() )) then
                 write(x_unit(diaryfile()),'(/,"Current RMS stress (z components only) = ",f0.6,":")') stress_norm
              end if

              do while (stress_norm >= stress_tol)
                stress_step = stress_step + 1
                if (i_access( diaryfile() )) then
                  if (stress_step < 10) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i1,":")') stress_step
                  elseif (stress_step < 100) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i2,":")') stress_step
                  elseif (stress_step < 1000) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i3,":")') stress_step
                  else
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i4,":")') stress_step
                  end if
                end if

                ! compare to previous steps direction to speed up minimization procedure,
                !   only need 3rd column for z axis motion (j=3) 
                do i = 1,3
                   check = step(i,j)*s(i,j)
                  if (check < 0) then
                      step(i,j) = s(i,j)
                  else
                     step(i,j) = step(i,j)+s(i,j)
                  end if
                end do

                ! update positions
                v = v - au_2_kbar*lattice_prefactor*matmul(v,step)
                call broadcast(CONFIG,v) ; if (error()) goto 400
                lat = lattice(v(1,:)/latc,v(2,:)/latc,v(3,:)/latc,latc) ; if (error()) goto 400
                cr = crystal(x_atoms(cr),lat,x_name(cr)) ; if (error()) goto 400
                call save(cr); if (error()) goto 400
                ext = external(cr = cr) ; if (error()) goto 400
                changed = .TRUE.
                call update(cfg,ext)  ; if (error()) goto 400

                !  reset atomic relaxation step number
                rx%o%step = 1
                changed = optimize_positions(cfg,rx) .or. changed ; if (error()) goto 400
                call diary_energy(cfg,all=.false.)

                s = x_stress_tensor(cfg) ; if (error()) goto 400
                call diary_stress_tensor(cfg) ; if (error()) goto 400
               !calculate stress contributions that effect the z-axis
                stress_norm = au_2_kbar*sqrt(sum(s(:,3)**2)/3.0_double)
                if (i_access( diaryfile() )) then
                   write(x_unit(diaryfile()),'(/,"Current RMS stress (z components only) = ",f0.6,":")') stress_norm
                end if

                if (stress_step > max_stress_steps) exit

              end do

!perform full atomic relaxation
              rx%o%max_steps = max_steps
              !  reset atomic relaxation step number
              rx%o%step = 1
              changed = optimize_positions(cfg,rx) .or. changed; if (error()) goto 400

           CASE (SHAPE_QUENCH)

              max_steps = rx%o%max_steps
              rx%o%max_steps = 5
              changed = optimize_positions(cfg,rx) ; if (error()) goto 400

              step = 0.0_double
              v(1,:) = x_lattice_vector(lat,1)  ; if (error()) goto 400
              v(2,:) = x_lattice_vector(lat,2)  ; if (error()) goto 400
              v(3,:) = x_lattice_vector(lat,3)  ; if (error()) goto 400
              latc = x_lattice_constant(lat)  ; if (error()) goto 400

              stress_step = 0
              s = x_stress_tensor(cfg) ; if (error()) goto 400
              call diary_stress_tensor(cfg) ; if (error()) goto 400
              stress_norm = au_2_kbar*sqrt(sum(s**2)/9.0_double)
              if (i_access( diaryfile() )) then
                 write(x_unit(diaryfile()),'(/,"Current RMS stress = ",f0.6,":")') stress_norm
              end if

              do while (stress_norm >= stress_tol)
                stress_step = stress_step + 1
                if (i_access( diaryfile() )) then
                  if (stress_step < 10) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i1,":")') stress_step
                  elseif (stress_step < 100) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i2,":")') stress_step
                  elseif (stress_step < 1000) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i3,":")') stress_step
                  else
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i4,":")') stress_step
                  end if
                end if

                ! compare to previous steps direction to speed up minimization procedure
                do i = 1,3
                   do j = 1,3
                   check = step(i,j)*s(i,j)
                  if (check < 0) then
                      step(i,j) = s(i,j)
                  else
                     step(i,j) = step(i,j)+s(i,j)
                  end if
                  end do
                end do

                ! update lattice vectors
                v = v - au_2_kbar*lattice_prefactor*matmul(v,step)
                call broadcast(CONFIG,v) ; if (error()) goto 400
                lat = lattice(v(1,:)/latc,v(2,:)/latc,v(3,:)/latc,latc) ; if (error()) goto 400
                cr = crystal(x_atoms(cr),lat,x_name(cr)) ; if (error()) goto 400
                call save(cr); if (error()) goto 400
                ext = external(cr = cr) ; if (error()) goto 400
                changed = .TRUE.
                call update(cfg,ext)  ; if (error()) goto 400

                !  reset atomic relaxation step number
                rx%o%step = 1
                changed = optimize_positions(cfg,rx) .or. changed; if (error()) goto 400
                call diary_energy(cfg,all=.false.)

                s = x_stress_tensor(cfg) ; if (error()) goto 400
                call diary_stress_tensor(cfg) ; if (error()) goto 400
                stress_norm = au_2_kbar*sqrt(sum(s**2)/9.0_double)
                if (i_access( diaryfile() )) then
                   write(x_unit(diaryfile()),'(/,"Current RMS stress = ",f0.6,":")') stress_norm
                end if

                if (stress_step > max_stress_steps) exit

              end do

!perform full atomic relaxation
              rx%o%max_steps = max_steps
              !  reset atomic relaxation step number
              rx%o%step = 1
              changed = optimize_positions(cfg,rx) .or. changed; if (error()) goto 400

! linear optimization of lattice
           CASE (SHAPE_LINEAR)

              changed = optimize_positions(cfg,rx) ; if (error()) goto 400

              v(1,:) = x_lattice_vector(lat,1)  ; if (error()) goto 400
              v(2,:) = x_lattice_vector(lat,2)  ; if (error()) goto 400
              v(3,:) = x_lattice_vector(lat,3)  ; if (error()) goto 400
              latc = x_lattice_constant(lat)  ; if (error()) goto 400

              stress_step = 0
              s = x_stress_tensor(cfg) ; if (error()) goto 400
              call diary_stress_tensor(cfg) ; if (error()) goto 400
              stress_norm = au_2_kbar*sqrt(sum(s**2)/9.0_double)
              if (i_access( diaryfile() )) then
                 write(x_unit(diaryfile()),'(/,"Current RMS stress = ",f0.6,":")') stress_norm
              end if

              do while (stress_norm >= stress_tol)
                stress_step = stress_step + 1
                if (i_access( diaryfile() )) then
                  if (stress_step < 10) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i1,":")') stress_step
                  elseif (stress_step < 100) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i2,":")') stress_step
                  elseif (stress_step < 1000) then
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i3,":")') stress_step
                  else
                     write(x_unit(diaryfile()),'(/,"Cell shape relaxation step #",i4,":")') stress_step
                  end if
                end if

                ! update lattice vectors
                v = v - au_2_kbar*lattice_prefactor*matmul(v,s)
                call broadcast(CONFIG,v) ; if (error()) goto 400
                lat = lattice(v(1,:)/latc,v(2,:)/latc,v(3,:)/latc,latc) ; if (error()) goto 400
                cr = crystal(x_atoms(cr),lat,x_name(cr)) ; if (error()) goto 400
                call save(cr); if (error()) goto 400
                ext = external(cr = cr) ; if (error()) goto 400
                changed = .TRUE.
                call update(cfg,ext)  ; if (error()) goto 400

                !reset atomic relaxation step number
                rx%o%step = 1
                changed = optimize_positions(cfg,rx) .or. changed; if (error()) goto 400
                call diary_energy(cfg,all=.false.)

                s = x_stress_tensor(cfg) ; if (error()) goto 400
                call diary_stress_tensor(cfg) ; if (error()) goto 400
                stress_norm = au_2_kbar*sqrt(sum(s**2)/9.0_double)
                if (i_access( diaryfile() )) then
                   write(x_unit(diaryfile()),'(/,"Current RMS stress = ",f0.6,":")') stress_norm
                end if

                if (stress_step > max_stress_steps) exit

             end do

        END SELECT

400     call glean(thy(lat))
300     call glean(thy(cr))
200     call glean(thy(ext))

100     call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::relax_structure_rx")) continue

      end function


      function relax_brancher(cfg) result(changed)
!doc$ function optimize_positions(cfg) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        logical :: changed
!       modifies: cfg
!       effects: Constructs rx and uses it to relax atom positions in cfg.
!       errors: Passes errors.

!cod$        
        type(relax_obj) :: rx
        
        call my(cfg)
        call my(relax(),rx) ; if (error()) goto 100
        
        branch: SELECT CASE (rx%o%method)
           CASE (NONE)
              changed = .false.
           CASE (STEEPEST__DESCENT)
              changed = steepest_descent(cfg,rx)
           CASE (CONJUGATE__GRADIENT)
              changed = conjugate_gradient(cfg,rx)
           CASE (QUENCH__MINIMIZATION)
              changed = quench_minimization(cfg,rx)
           END SELECT branch

100     call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::relax_brancher")) continue

      end function

      function relax_brancher_rx(cfg,rx) result(changed)
!doc$ function optimize_positions(cfg,rx) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        logical :: changed
!       modifies: cfg and rx
!       effects: Uses rx to relax atom positions in cfg.
!       errors: Passes errors.

!cod$        
        call my(cfg)
        call my(rx)
        
        branch: SELECT CASE (rx%o%method)
           CASE (NONE)
              changed = .false.
           CASE (STEEPEST__DESCENT)
              changed = steepest_descent(cfg,rx)     ; if (error()) goto 100
           CASE (CONJUGATE__GRADIENT)
              changed = conjugate_gradient(cfg,rx)   ; if (error()) goto 100
           CASE (QUENCH__MINIMIZATION)
              changed = quench_minimization(cfg,rx)  ; if (error()) goto 100
           END SELECT branch

100     call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::relax_brancher_rx")) continue

      end function

      function relax_sd(cfg,rx) result(changed)
!doc$ function steepest_descent(cfg,rx) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        logical :: changed
!       modifies: cfg and rx
!       effects: Performs steepest descent relaxation of atom positions in cfg.
!       errors: Passes errors.

!cod$
        integer :: na
        real(double), dimension(:,:), allocatable :: frcs, pos_cart
        
        call my(cfg)
        call my(rx)

        call own_i(rx)
        rx%o%g = x_ghost()

        na = x_n_atoms(cfg)
        allocate ( pos_cart(3,na) )
        allocate ( frcs(3,na) )

        changed = .false.
        do
          frcs = x_forces(cfg)  ; if (error()) goto 100
          call diary(cfg,rx%o%step)
          if (.not.continue_relax_i(cfg,rx)) exit
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            exit
          end if
          call x_cart_positions(cfg,pos_cart)
          pos_cart = pos_cart + rx%o%sd_con*frcs
          call update_config(cfg,pos_cart)             ; if (error()) goto 100
          changed = .true.
          rx%o%step = rx%o%step + 1
        end do

100     if (allocated( frcs )) deallocate ( frcs )
        if (allocated( pos_cart )) deallocate ( pos_cart )

        call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::relax_sd")) continue

      end function

      function relax_cg(cfg,rx) result(changed)
!doc$ function conjugate_gradient(cfg,rx) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        logical :: changed
!       modifies: cfg and rx
!       effects: Performs conjugate gradient relaxation of atom positions in cfg.
!       errors: Passes errors.

!cod$
        integer :: na, n_linmin
        real(double) :: fp, gg, dgg, gam
        real(double), dimension(:,:), allocatable :: xi, g, h

        call my(cfg)
        call my(rx)

        call own_i(rx)
        rx%o%g = x_ghost()

        na = x_n_atoms(cfg)
        allocate ( xi(3,na), g(3,na), h(3,na) )

        fp = x_cell_energy(cfg)
        xi = -x_forces(cfg)
        g = -xi
        h = g
        xi = h
        n_linmin = 0

        call diary(cfg,rx%o%step)

        changed = .false.
        do
          if (.not.continue_relax_i(cfg,rx)) exit
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            exit
          end if
          call linmin_i(cfg,rx,xi) ; if (error()) goto 100
          n_linmin = n_linmin + 1
          changed = .true.
          fp = x_cell_energy(cfg)
          xi = -x_forces(cfg)
          if (mod(n_linmin,5) /= 0) then  ! periodically restart conjugate gradient process
            gg = atom_dot_product(g)
!            dgg = atom_dot_product(xi+g,g)
            dgg = atom_dot_product(xi)
            gam = dgg/gg
            g = -xi
            h = g + gam*h
            xi = h
          else
            g = -xi
            h = g
            xi = h
          end if
        end do

100     deallocate ( xi, g, h )
        
        call glean(thy(cfg))
        call glean(thy(rx))
        
        if (error("Exit relax_mod::relax_cg")) continue

      end function

      function relax_qm(cfg,rx) result(changed)
!doc$ function quench_minimization(cfg,rx) result(changed)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        logical :: changed
!       modifies: cfg and rx
!       effects: Performs quench minimization relaxation of atom positions in cfg.
!       errors: Passes errors.

!cod$
        logical :: found, match
        integer :: ia, iat, na
        character(line_len) :: tag
        character(tag_sz), dimension(:), allocatable :: tags
        real(double) :: fdotv, mass
        real(double), dimension(:), allocatable :: masses
        real(double), dimension(:,:), allocatable :: accelerations, forces, positions, velocities ! cartesian representation

        call my(cfg)
        call my(rx)

        call own_i(rx)
        rx%o%g = x_ghost()

        na = x_n_atoms(cfg)
        allocate( tags(na) )
        allocate( masses(na) )
        allocate( accelerations(3,na) )
        allocate( forces(3,na) )
        allocate( positions(3,na) )
        allocate( velocities(3,na) )

        ! read the atom masses
        do ia = 1,na
          tags(ia) = x_type(cfg,ia)
          tag = "atom_mass_"//tags(ia)
          call arg(trim(tag),mass,found)
          if (error(.not.found,"ERROR: atom_mass was not found")) goto 100
          if (error(mass <= 0.0_double,"ERROR: atom_mass <= 0")) goto 100
          masses(ia) = AMU_2_ELECTRON_MASS*mass
        end do

        ! diary the atom masses
        if ( i_access(diaryfile()) ) then
          do ia = 1,na
            if (ia == 1) then
              write(x_unit(diaryfile()),"(/,t4,'Mass for atom',i6,' =',2g14.5)") ia, masses(ia), masses(ia)/AMU_2_ELECTRON_MASS
            else
              match = .false.
              do iat = 1,(ia-1)
                if (tags(ia) == tags(iat)) then
                  match = .true.
                  exit
                end if
              end do
              if (.not.match) then
                write(x_unit(diaryfile()),"(t4,'Mass for atom',i6,' =',2g14.5)") ia, masses(ia), masses(ia)/AMU_2_ELECTRON_MASS
              end if
            end if
          end do
        end if
        deallocate( tags )

        ! get/compute the initial positions, forces, velocities, and accelerations
        call x_cart_positions(cfg,positions)
        forces = x_forces(cfg)
        do ia = 1,na
          velocities(:,ia) = 0.0_double
          accelerations(:,ia) = forces(:,ia)/masses(ia)
        end do

        changed = .false.
        do

          ! diary step results
          call diary(cfg,rx%o%step)

          ! check for tolerance
          if (.not.continue_relax_i(cfg,rx)) exit

          ! check for user-initiated stop
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            exit
          end if

          ! update cfg and associated quantities
          positions = positions + rx%o%time_step*(velocities + (rx%o%time_step/2.0_double)*accelerations)
          velocities = velocities + (rx%o%time_step/2.0_double)*accelerations
          call update_config(cfg,positions) ; if (error()) goto 100
          changed = .true.
          rx%o%step = rx%o%step + 1
          forces = x_forces(cfg)
          do ia = 1,na
            accelerations(:,ia) = forces(:,ia)/masses(ia)
          end do

          ! update velocities
          velocities = velocities + (rx%o%time_step/2.0_double)*accelerations
          do ia = 1,na
            fdotv = dot_product(forces(:,ia),velocities(:,ia))
            if (fdotv <= 0.0_double) then
              velocities(:,ia) = 0.0_double
            else
              velocities(:,ia) = fdotv*forces(:,ia)/dot_product(forces(:,ia),forces(:,ia))
            end if
          end do
           
        end do

100     if (allocated( tags )) deallocate( tags )
        if (allocated( masses )) deallocate( masses )
        if (allocated( accelerations )) deallocate( accelerations )
        if (allocated( forces )) deallocate( forces )
        if (allocated( positions )) deallocate( positions )
        if (allocated( velocities )) deallocate( velocities )

        call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::relax_qm")) continue

      end function

      subroutine my_relax(rx)
!doc$ subroutine my(rx)
        type(relax_obj) :: rx

!cod$
        rx%ref = rx%ref + 1
        rx%o%ref = rx%o%ref + 1
      end subroutine

      subroutine my_new_relax(rxi,rx)
!doc$ subroutine my(rxi,rx)
        type(relax_obj) :: rxi, rx

!cod$
        rx%ref = 1
        rx%o => rxi%o
        rx%o%ref = rx%o%ref + 1
      end subroutine

      function thy_relax(rx) result(rxo)
!doc$ function thy(rx) result(rxo)
        type(relax_obj) :: rx, rxo

!cod$
        rx%ref = rx%ref - 1
        rx%o%ref = rx%o%ref - 1
        rxo%ref = rx%ref
        rxo%o => rx%o
      end function

      subroutine glean_relax(rx)
!doc$ subroutine glean(rx)
        type(relax_obj) :: rx

!cod$
        if (rx%o%ref < 1) then
           deallocate ( rx%o )
        end if
      end subroutine

      subroutine bequeath_relax(rx)
!doc$ subroutine bequeath(rx)
        type(relax_obj) :: rx

!cod$
        continue
      end subroutine

      subroutine assign_relax(rx,rx2)
!doc$ subroutine assignment(=)(rx,rx2)
        type(relax_obj), intent(inout) :: rx
        type(relax_obj), intent(in) :: rx2

!cod$
        type(relax_obj) :: rxt
        call my(rx2)
        rxt%o => rx%o
        rx%o%ref = rx%o%ref - rx%ref
        rx%o => rx2%o
        rx%o%ref = rx%o%ref + rx%ref
        call glean(rxt)
        call glean(thy(rx2))
      end subroutine

      subroutine diary_relax(rx)
!doc$ subroutine diary(rx)
        type(relax_obj) :: rx
!       modifies: diary stream
!       effects: Writes information about the optimization method.

!cod$
        call my(rx)

        if (i_access( diaryfile() )) then

          select case (rx%o%lattice_relaxation)
             case (NONE)
               continue
!               write(x_unit(diaryfile()),'(/,"No optimization of lattice vectors")')
             case (SHAPE_LINEAR)
               write(x_unit(diaryfile()),'(/,"The cell shape is optimized")')
             case (SHAPE_QUENCH)
                write(x_unit(diaryfile()),'(/,"The cell shape is optimized with quenched motion")')
             case (SHAPE_QUENCH_Z)
                write(x_unit(diaryfile()),'(/,"The cell shape is optimized in z-direction only with quenched motion")')
           end select 


           wr: select case(rx%o%method)
              case (NONE)
                 write(x_unit(diaryfile()),'(/,"No optimization of structure")')

              case (STEEPEST__DESCENT)
                 write(x_unit(diaryfile()),'(/,t4,"Steepest Descents with constant: ",es12.5)') rx%o%sd_con
                 write(x_unit(diaryfile()),'(t4,"Maximum number of optimization steps:",i4)') rx%o%max_steps 
                 write(x_unit(diaryfile()),'(t4,"Number of optimization steps taken:  ",i4)') rx%o%step 
                 write(x_unit(diaryfile()),'(t4,"Tolerance for RMS force value:",es11.4)') rx%o%force_tol
                 write(x_unit(diaryfile()),'(t4,"Current RMS force value:      ",es11.4)') rx%o%force_nrm
              case (CONJUGATE__GRADIENT)
                 write(x_unit(diaryfile()),'(/,t4,"Conjugate gradient optimization of structure")') 
                 write(x_unit(diaryfile()),'(t4,"Maximum number of optimization steps:",i4)') rx%o%max_steps 
                 write(x_unit(diaryfile()),'(t4,"Number of force calls made:          ",i4)') rx%o%step 
                 write(x_unit(diaryfile()),'(t4,"Tolerance for RMS force value:",es11.4)') rx%o%force_tol
                 write(x_unit(diaryfile()),'(t4,"Current RMS force value:      ",es11.4)') rx%o%force_nrm
              case (QUENCH__MINIMIZATION)
                 write(x_unit(diaryfile()),'(/,t4,"Quench Minimization with time_step: ",es11.5)') rx%o%time_step
                 write(x_unit(diaryfile()),'(t4,"Maximum number of optimization steps:",i4)') rx%o%max_steps 
                 write(x_unit(diaryfile()),'(t4,"Number of optimization steps taken:  ",i4)') rx%o%step 
                 write(x_unit(diaryfile()),'(t4,"Tolerance for RMS force value:",es11.4)') rx%o%force_tol
                 write(x_unit(diaryfile()),'(t4,"Current RMS force value:      ",es11.4)') rx%o%force_nrm
              end select wr

        end if

        call glean(thy(rx))

      end subroutine

      subroutine own_i(rx)
        type(relax_obj), intent(inout) :: rx
        type(relax_obj) :: rxtmp
        if (rx%ref < rx%o%ref) then
          allocate ( rxtmp%o )
          rxtmp%o%ref = 0
          rxtmp%o%lattice_relaxation = rx%o%lattice_relaxation
          rxtmp%o%method = rx%o%method
          select case (rxtmp%o%method)
          case (STEEPEST__DESCENT,CONJUGATE__GRADIENT,QUENCH__MINIMIZATION)
            rxtmp%o%max_steps = rx%o%max_steps
            rxtmp%o%force_tol = rx%o%force_tol
            rxtmp%o%step = rx%o%step
            rxtmp%o%force_nrm = rx%o%force_nrm
          end select
          select case (rxtmp%o%method)
          case (STEEPEST__DESCENT)
            rxtmp%o%sd_con = rx%o%sd_con
          case (QUENCH__MINIMIZATION)
            rxtmp%o%time_step = rx%o%time_step
          end select
          rxtmp%o%g = rx%o%g
          rx%o%ref = rx%o%ref - rx%ref
          rx%o => rxtmp%o
          rx%o%ref = rx%o%ref + rx%ref
        end if
      end subroutine

      subroutine linmin_i(cfg,rx,displace)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        real(double), dimension(:,:), intent(inout) :: displace

        logical :: found, bracket
        integer :: na, status
        real(double) :: alpha_a, alpha_b, alpha_c
        real(double) :: energy_a, energy_c
        real(double) :: grad_a, grad_b, grad_c, grad_init
        real(double) :: tst
        real(double), dimension(:,:), allocatable :: pos_cart,pos_cart_start
        real(double), dimension(:,:), allocatable :: direction

        call my(cfg)
        call my(rx)

        na = x_n_atoms(cfg)

        allocate(pos_cart(3,na),pos_cart_start(3,na),direction(3,na),STAT=status)
        if (error(status /= 0,"limnin_i: allocate error")) goto 100

        if (error(size(displace,1) /= 3 .or. size(displace,2) /= na,&
             &"linmin_i: displace array wrong dimensions")) goto 100

        if (i_access(output)) write(x_unit(output),'(/,"Beginning line minimization")')
        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Beginning line minimization")')

        ! positions and direction for the initial cfg
        direction = displace
        grad_init = sqrt(atom_dot_product(direction))
        direction = direction/grad_init
        grad_init = atom_dot_product(x_forces(cfg),direction)
        call x_cart_positions(cfg,pos_cart_start)
        
        ! initial values for point a
        alpha_a = 0.0_double
        grad_a = atom_dot_product(x_forces(cfg),direction)
        energy_a = x_cell_energy(cfg)

        ! value of alpha_c on the other side of the minimum
        alpha_c = grad_init

        do
          pos_cart = pos_cart_start + alpha_c*direction
          call update_config(cfg,pos_cart) ; if (error()) goto 100
          call warn('out of update config')
          rx%o%step = rx%o%step + 1
          call diary(cfg,rx%o%step)
          grad_c = atom_dot_product(x_forces(cfg),direction)
          energy_c = x_cell_energy(cfg)
          tst = grad_a*grad_c
          if (.not.continue_relax_i(cfg,rx,direction)) then
            found = .true.
            alpha_b = alpha_c
            exit
          end if
          if (tst < 0.0_double) then
            found = .false.
            bracket = .true.
            exit
          end if
          alpha_b = estimate_minimum_i(alpha_a,alpha_c,grad_a,grad_c)
          if (energy_c < energy_a) then
            alpha_a = alpha_c
            grad_a = grad_c
            energy_a = energy_c
          end if
          alpha_c = alpha_b
          if (i_access(output)) then
             write(x_unit(output),'(/,"Increase step to bracket minimum")')
             write(x_unit(output),'(t2,"alpha_a,alpha_c:",2f10.3)') alpha_a, alpha_c
          end if
          if (i_access(diaryfile())) then
             write(x_unit(diaryfile()),'(/,"Increase step to bracket minimum")')
             write(x_unit(diaryfile()),'(t2,"alpha_a,alpha_c:",2f10.3)') alpha_a, alpha_c
          end if
        end do

        if (.not.found) then
          do
            alpha_b = alpha_a + grad_a*(alpha_a - alpha_c)/(grad_c - grad_a)
            pos_cart = pos_cart_start + alpha_b*direction
            call update_config(cfg,pos_cart) ; if (error()) goto 100
            rx%o%step = rx%o%step + 1
            call diary(cfg,rx%o%step)
            if (.not.continue_relax_i(cfg,rx,direction)) exit
            grad_b = atom_dot_product(x_forces(cfg),direction)
            if (grad_a*grad_b < 0.0_double) then
              grad_c = grad_b
              alpha_c = alpha_b
            else
              grad_a = grad_b
              alpha_a = alpha_b
            end if
          end do
        end if

        displace = alpha_b*direction

        if (i_access(output)) write(x_unit(output),'(/,"Completed line minimization")')
        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Completed line minimization")')

100     deallocate(pos_cart,pos_cart_start,direction)
        call glean(thy(cfg))
        call glean(thy(rx))

        if (error("Exit relax_mod::linmin_i")) continue

      end subroutine

      function continue_relax_i(cfg,rx,direction) result(go)
        type(config_sc_obj), intent(inout) :: cfg
        type(relax_obj), intent(inout) :: rx
        real(double), dimension(:,:), intent(in), optional :: direction
        logical :: go

        integer :: na
        real(double) :: sum, line_grad

        call my(cfg)
        call my(rx)

        na = x_n_atoms(cfg)

        sum = atom_dot_product(x_forces(cfg))
        sum = sqrt(sum/(3.0_double*dble(na)))
        rx%o%force_nrm = sum

        call diary(rx)
        if (i_access(output)) then
          write(x_unit(output),'(t3,"RMS force value = ",es11.4,";",3x,"Tolerance = ",es11.4)') &
                                rx%o%force_nrm, rx%o%force_tol
          write(x_unit(output),'(t3,"Current number of force call = ",i5,";",3x,"Limit = ",i5)') &
                                rx%o%step, rx%o%max_steps
        end if

        if (.not.present(direction)) then
          go = ( rx%o%force_nrm >= rx%o%force_tol )
        else
          line_grad = atom_dot_product(x_forces(cfg),direction)/sqrt(3.0_double*na)
          line_grad = abs(line_grad)
          if (i_access(output)) write(x_unit(output),'(t3,"RMS force along direction:",es12.5)') line_grad
          go = ( line_grad >= rx%o%force_tol )
        end if
        go = ( go .and. (rx%o%step < rx%o%max_steps) )

        if (error('Exit relax_mod:continue_relax_i')) goto 100

100     call glean(thy(cfg))
        call glean(thy(rx))

      end function

      function estimate_minimum_fancy_i(alpha_a,alpha_c,energy_a,energy_c,grad_a,grad_c) result(alpha)
        real(double), intent(in) :: alpha_a, alpha_c, energy_a, energy_c, grad_a, grad_c

        real(double) :: alpha
        real(double) :: ec, ga, gc, xc
        real(double) a, b, c, root1, root2, test

        ga = -grad_a
        gc = -grad_c
        xc = alpha_c - alpha_a
        ec = energy_c - energy_a

        c = ga

        b = (3.0_double*ec - (gc + 2.0_double*ga)*xc)/xc**2

        a = (ec - c*xc - b*xc**2)/xc**3

        test = b**2 - 3.0_double*a*c
        if (test <= 0.0_double) then          
           root1 = -ga*xc/(gc - ga)
           alpha = alpha_a + root1
        else
           root1 = (1.0_double/(3.0_double*a))*(-b + sqrt(test))
           root2 = (1.0_double/(3.0_double*a))*(-b - sqrt(test))
           if (min(root1,root2) > xc) then
              alpha = alpha_a + min(root1,root2)
           else if (max(root1,root2) > xc) then
              alpha = alpha_c + max(root1,root2)
           else
              alpha = alpha_a + 2.0_double*(alpha_c - alpha_a)
           end if
        end if

      end function

      function estimate_minimum_i(alpha_a,alpha_c,grad_a,grad_c) result(alpha)
        real(double), intent(in) :: alpha_a, alpha_c, grad_a, grad_c

        real(double) :: alpha
        real(double) :: ga, gc, root1, xc

        ga = -grad_a
        gc = -grad_c
        xc = alpha_c - alpha_a
        root1 = -ga*xc/(gc - ga)
        alpha = alpha_a + root1

      end function

    end module

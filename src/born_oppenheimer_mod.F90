!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

      module born_oppenheimer_mod
!doc$ module born_oppenheimer_mod

!     born_oppenheimer_mod provides procedures for performing molecular dynamics simulations.

      use kind_mod
      use ghost_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use diary_mod
      use arg_mod
      use math_mod
      use interrupt_mod
      use config_sc_mod
      use shortcuts_mod

!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = CONFIG

      real(double), parameter :: K_BOLTZMAN = 6.3333e-06_double
      real(double), parameter :: AMU_2_ELECTRON_MASS = 9.1144424e+02_double
      real(double), parameter :: FS_2_ARU = 1.0_double/4.8377687e-02_double

      integer, parameter :: NONE = 0
      integer, parameter :: NVE = 1
      integer, parameter :: NVT_RESCALE = 2
      integer, parameter :: NVT_ANDERSON = 3
      integer, parameter :: NVT_HOOVER = 4

      type, public :: bo_ave_obj
         integer :: n_ave_steps
         real(double) :: sum_temperature,sum_temperature_sq
         real(double) :: sum_kin_energy,sum_kin_e_sq
         real(double) :: sum_pot_energy,sum_pot_e_sq
         real(double) :: sum_tot_energy,sum_tot_e_sq
         logical :: ave_pressure, ave_stress_tensor
         real(double) :: sum_pressure, sum_pressure_sq
         real(double), dimension(3,3) :: sum_stress_tensor, sum_stress_tensor_sq
      end type bo_ave_obj

      type, public :: bo_state_obj
         integer :: natoms
         real(double) :: time
         real(double) :: kin_energy,temperature
         real(double) :: pot_energy
         real(double) :: tot_energy
         real(double), dimension(:,:), pointer :: pos_cart
         real(double), dimension(:,:), pointer :: velocities
         real(double), dimension(:), pointer :: masses
         real(double) :: hoover_mass
         real(double) :: hoover_xi
         real(double) :: hoover_s
      end type bo_state_obj

      type :: bo_rep
         integer :: ref
         type(ghost) :: g
         integer :: md_method
         integer :: n_md_steps
         integer :: n_skip_steps
         integer :: max_md_steps
         integer :: temp_mod_freq
         integer :: generate_velocities
         real(double) :: time_step
         real(double) :: init_temp
         real(double) :: desired_temp
         type(file_obj) :: md_file
         type(file_obj) :: init_vel_file
         type(file_obj) :: new_vel_file
         type(bo_state_obj), pointer :: current
         type(bo_ave_obj), pointer :: ave
      end type bo_rep

      type, public :: bo_obj
         private
         integer :: ref
         type(bo_rep), pointer :: o
      end type bo_obj

!doc$
      public :: born_oppenheimer
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: born_oppenheimer_dynamics
      public :: bo_velocity_verlet
      public :: bo_rescale
      public :: bo_anderson
      public :: diary
      
!cod$

      interface born_oppenheimer
         module procedure constructor_bo
      end interface
      interface my
         module procedure my_bo, my_new_bo
      end interface
      interface thy
         module procedure thy_bo
      end interface
      interface glean
         module procedure glean_bo
      end interface
      interface bequeath
         module procedure bequeath_bo
      end interface
      interface assignment(=)
         module procedure assign_bo
      end interface
      interface born_oppenheimer_dynamics
         module procedure bo_brancher
      end interface
      interface diary
         module procedure diary_bo, diary_bo_state, diary_bo_ave
      end interface

      contains

      subroutine my_bo(bo)
!doc$ subroutine my(bo)
        type(bo_obj) :: bo

!cod$
        bo%ref = bo%ref + 1
        bo%o%ref = bo%o%ref + 1

      end subroutine

      subroutine my_new_bo(boi,bo)
!doc$ subroutine my(boi,bo)
      type(bo_obj) :: boi, bo

!cod$
        bo%ref = 1
        bo%o => boi%o
        bo%o%ref = bo%o%ref + 1
      end subroutine

      function thy_bo(bo) result(boo)
!doc$ function thy(bo) result(boo)
        type(bo_obj) :: bo, boo

!cod$
        bo%ref = bo%ref - 1
        bo%o%ref = bo%o%ref - 1
        boo%ref = bo%ref
        boo%o => bo%o
      end function

      subroutine glean_bo(bo)
!doc$ subroutine glean(bo)
        type(bo_obj) :: bo

!cod$
        if (bo%o%ref < 1) then
          if (associated( bo%o%current%pos_cart )) deallocate( bo%o%current%pos_cart )
          if (associated( bo%o%current%velocities )) deallocate( bo%o%current%velocities )
          if (associated( bo%o%current%masses )) deallocate( bo%o%current%masses )
          deallocate( bo%o%current )
          deallocate( bo%o%ave )
          deallocate( bo%o )
        end if
      end subroutine

      subroutine bequeath_bo(bo)
!doc$ subroutine bequeath(bo)
        type(bo_obj) :: bo

!cod$
        continue
      end subroutine

      subroutine assign_bo(bo,bo2)
!doc$ subroutine assignment(=) (bo,bo2)
        type(bo_obj), intent(inout) :: bo
        type(bo_obj), intent(in) :: bo2
       
!cod$
        type(bo_obj) :: bot
        call my(bo2)
        bot%o => bo%o
        bo%o%ref = bo%o%ref - bo%ref
        bo%o => bo2%o
        bo%o%ref = bo%o%ref + bo%ref
        call glean(bot)
        call glean(thy(bo2))
      end subroutine

      function constructor_bo(cfg,time_step) result(bo)
!doc$ function born_oppenheimer(cfg,time_step) result(bo)
        type(config_sc_obj) :: cfg
        real(double), optional :: time_step
        type(bo_obj) :: bo
!       effects: sets the structural optimization parameters

!cod$
        logical :: exist_file, found, match
        character(line_len) :: tag
        character(tag_sz), dimension(:), allocatable :: tags
        integer :: i, j, na, ios, status
        real(double) :: mass_temp

        call my(cfg)

        bo%ref = 0
        na = x_n_atoms(cfg)
        allocate( bo%o, STAT=status)
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate(bo%o%ave, STAT=status)
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate(bo%o%current, STAT=status)
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate(bo%o%current%pos_cart(3,na), STAT=status)
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate(bo%o%current%velocities(3,na), STAT=status)
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate(bo%o%current%masses(na), STAT=status )
        if (error(status /= 0,"ERROR: allocate failure")) goto 100

        bo%o%g = x_ghost()
        bo%o%ref = 0

        bo%o%n_md_steps = 0

        ! read the parameter defining the method
        !  =0: no molecular dynamics
        !  =1: NVE
        !  =2: NVT_RESCALE
        !  =3: NVT_ANDERSON        
        !  =4: NVT_HOOVER
        
        call arglc("md_method",tag,found)
        if (.not.found) tag = "none"
        select case (tag(1:len_trim(tag)))
        case ("none")
           bo%o%md_method = NONE
        case ("nve","verlet")
           bo%o%md_method = NVE
        case ("nvt_rescale","rescale")
           bo%o%md_method = NVT_RESCALE
        case ("nvt_bo_anderson","bo_anderson","stochastic")
           bo%o%md_method = NVT_ANDERSON
        case ("nvt_hoover","hoover")
           bo%o%md_method = NVT_HOOVER
        case default
           if (error(.true.,"ERROR: unrecognized md_method")) goto 100
        end select
              
        if (bo%o%md_method .ne. NONE) then

           ! read md parameters 
           call arg("md_steps",bo%o%max_md_steps,found)
           if (.not.found) bo%o%max_md_steps = 0
           if (error(bo%o%max_md_steps<0,"ERROR: max_md_steps < 0")) goto 100

           call arg("md_skip_steps",bo%o%n_skip_steps,found)
           if (.not.found) bo%o%n_skip_steps=0
           if (error(bo%o%n_skip_steps<0,"ERROR: n_skip_steps < 0")) goto 100

           if (present(time_step)) then
              bo%o%time_step = time_step
           else
              call arg("md_time_step",bo%o%time_step,found)
              if (.not.found) bo%o%time_step = 100.
           end if
           
           if (error(bo%o%time_step<0,"ERROR: time_step < 0")) goto 100

           call arglc("md_gen_velocities",tag,found)
           if (.not.found) tag = "yes"
           select case (tag(1:len_trim(tag)))
           case ("yes")
              bo%o%generate_velocities = 1
           case ("no")
              bo%o%generate_velocities = 0
           case default
              if (error(.true.,"ERROR: md_gen_velocities")) goto 100
           end select

           call arg("md_init_temp",bo%o%init_temp,found)
           if (.not.found) bo%o%init_temp = 0.0_double
           if (error(bo%o%init_temp<0,"ERROR: init_temp < 0")) goto 100
           
           call arg("md_desired_temp",bo%o%desired_temp,found)
           if (.not.found) bo%o%desired_temp = bo%o%init_temp
           if (error(bo%o%desired_temp<0,"ERROR: desired_temp < 0")) goto 100

           call arg("md_temp_freq",bo%o%temp_mod_freq,found)
           if (.not.found) bo%o%temp_mod_freq=5
           if (error(bo%o%temp_mod_freq<1,"ERROR: temp_mod_freq < 1")) goto 100

           call arg("md_hoover_mass",bo%o%current%hoover_mass,found)
           if (.not.found) bo%o%current%hoover_mass=1000.
           if (error(bo%o%current%hoover_mass <= 0.0_double,"ERROR: hoover_mass < 0")) goto 100

           call arglc("pressure",tag,found)
           if (.not.found) tag = "off"
           select case (trim(tag))
           case ("on",".true.")
             bo%o%ave%ave_pressure = .true.
           case ("off",".false.")
             bo%o%ave%ave_pressure = .false.
           case default
             if (error(.true.,"ERROR: pressure tag was not recognized")) goto 100
           end select
           call arglc("stress_tensor",tag,found)
           if (.not.found) tag = "off"
           select case (trim(tag))
           case ("on",".true.")
             bo%o%ave%ave_stress_tensor = .true.
           case ("off",".false.")
             bo%o%ave%ave_stress_tensor = .false.
           case default
             if (error(.true.,"ERROR: pressure tag was not recognized")) goto 100
           end select

           call my(file(trim(md_trajectory_path)),bo%o%md_file)
           if (i_access(bo%o%md_file)) open(x_unit(bo%o%md_file),file=x_name(bo%o%md_file))
           
           call my(file(trim(new_velocity_path)),bo%o%new_vel_file)

           bo%o%current%natoms = x_n_atoms(cfg)
           
           ! read the atom masses
           allocate( tags(bo%o%current%natoms) )
           do i = 1,bo%o%current%natoms
             tags(i) = x_type(cfg,i)
             tag = "atom_mass_"//tags(i)
             call arg(trim(tag),mass_temp,found)
             if (error(.not.found,"ERROR: atom_mass was not found")) goto 100
             if (error(mass_temp <= 0.0_double,"ERROR: atom_mass <= 0")) goto 100
             bo%o%current%masses(i) = AMU_2_ELECTRON_MASS*mass_temp
           end do

           ! diary the atom masses
           if ( i_access(diaryfile()) ) then
             do i = 1,bo%o%current%natoms
               if (i == 1) then
                 write(x_unit(diaryfile()),"(/,t4,'Mass for atom',i6,' =',2g14.5)") &
                                                  & i, bo%o%current%masses(i), bo%o%current%masses(i)/AMU_2_ELECTRON_MASS
               else
                 match = .false.
                 do j = 1,i
                   if (tags(i) == tags(j)) then
                     match = .true.
                     exit
                   end if
                 end do
                 if (.not.match) then
                     write(x_unit(diaryfile()),"(t4,'Mass for atom',i6,' =',2g14.5)") &
                                                  & i, bo%o%current%masses(i), bo%o%current%masses(i)/AMU_2_ELECTRON_MASS
                 end if
               end if
             end do
           end if
           deallocate( tags )

           call x_cart_positions(cfg,bo%o%current%pos_cart)

           if (bo%o%generate_velocities == 0) then
              call my(file(trim(velocity_path)),bo%o%init_vel_file)
              if (i_access(bo%o%init_vel_file)) &
                   &inquire(file=x_name(bo%o%init_vel_file),exist=exist_file)
              if (i_comm(bo%o%init_vel_file)) call broadcast(FILE_SCOPE,exist_file)
              if (error(.not.exist_file,"ERROR: velocity file does not exist")) goto 100

              if (i_access(bo%o%init_vel_file)) &
                   &open(unit=x_unit(bo%o%init_vel_file), &
                   &file=x_name(bo%o%init_vel_file), iostat=ios)
              if (i_comm(bo%o%init_vel_file)) call broadcast(FILE_SCOPE,ios)
              if (error(ios /= 0,"ERROR: unable to open velocity file")) goto 100

              do i = 1,x_n_atoms(cfg)
                 if (i_access(bo%o%init_vel_file)) &
                      &read(x_unit(bo%o%init_vel_file),*,iostat=ios) &
                      &  bo%o%current%velocities(:,i)
                 if (i_comm(bo%o%init_vel_file)) call broadcast(FILE_SCOPE,ios)
                 if (error(ios /= 0,"ERROR: read error in velocity file")) goto 100
                 if (i_comm(bo%o%init_vel_file)) call broadcast(FILE_SCOPE,bo%o%current%velocities(:,i))
              end do
           else
              call init_velocities(bo%o%current%velocities,bo%o%current%masses,bo%o%init_temp)
           end if
           bo%o%current%hoover_s = 0.0_double
           bo%o%current%hoover_xi = 0.0_double

           bo%o%ave%n_ave_steps = 0
           bo%o%ave%sum_temperature = 0.0_double
           bo%o%ave%sum_temperature_sq = 0.0_double
           bo%o%ave%sum_kin_energy = 0.0_double
           bo%o%ave%sum_kin_e_sq = 0.0_double
           bo%o%ave%sum_pot_energy = 0.0_double
           bo%o%ave%sum_pot_e_sq = 0.0_double
           bo%o%ave%sum_tot_energy = 0.0_double
           bo%o%ave%sum_tot_e_sq = 0.0_double
           bo%o%ave%sum_pressure = 0.0_double
           bo%o%ave%sum_pressure_sq = 0.0_double
           do i = 1,3
              do j = 1,3
                 bo%o%ave%sum_stress_tensor(i,j) = 0.0_double
                 bo%o%ave%sum_stress_tensor_sq(i,j) = 0.0_double
              end do
           end do
           
           call update_mdinfo_i(cfg,bo%o)
           
        end if

100     if (allocated( tags )) deallocate( tags )

        call glean(thy(cfg))

        if (error("Exit born_oppenheimer_mod::constructor_bo")) continue

      end function
      
      subroutine diary_bo(bo)
!doc$ subroutine diary(bo)
        type(bo_obj) :: bo
!       effects: Writes information about the optimization

!cod$

        if (i_access( diaryfile() )) then

           wr: select case(bo%o%md_method)
              case (NONE)
                 write(x_unit(diaryfile()),'(/,"No molecular dynamics requested")')

              case (NVE)
                 write(x_unit(diaryfile()),'(/,"NVE MD (velocity Verlet) with time step: ",g12.5)') bo%o%time_step
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') bo%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') bo%o%n_skip_steps
                 if (bo%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') bo%o%n_md_steps
                 end if
                 call diary(bo%o%current)
                 call diary(bo%o%ave)
              case (NVT_RESCALE)
                 write(x_unit(diaryfile()),'(/,"NVT MD (velocity rescaling) with time step: ",g12.5)') bo%o%time_step
                 write(x_unit(diaryfile()),'(/," velocity rescaled every",i4," time steps")') bo%o%temp_mod_freq
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') bo%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') bo%o%n_skip_steps
                 if (bo%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') bo%o%n_md_steps
                 end if
                 call diary(bo%o%current)
                 call diary(bo%o%ave)
              case (NVT_ANDERSON)
                 write(x_unit(diaryfile()), &
                     &'(/,"NVT MD (stochastic velocities ala Anderson) with time step: ",g12.5)') &
                     &bo%o%time_step
                 write(x_unit(diaryfile()), &
                     &'(/," atom velocity rescaled at each time step with probability ", g12.5)') &
                     &1.0_double/bo%o%temp_mod_freq
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') bo%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') bo%o%n_skip_steps
                 if (bo%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') bo%o%n_md_steps
                 end if
                 call diary(bo%o%current)
                 call diary(bo%o%ave)
              case (NVT_HOOVER)
                 write(x_unit(diaryfile()), &
                     &'(/,"NVT MD (Hoover implementation of Nose thermostat) with time step: ",g12.5)') &
                     & bo%o%time_step
                 write(x_unit(diaryfile()),'(/," fictitious temperature-coupling mass ", g12.5)') bo%o%current%hoover_mass
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') bo%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') bo%o%n_skip_steps
                 if (bo%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') bo%o%n_md_steps
                 end if
                 call diary(bo%o%current)
                 write(x_unit(diaryfile()),'("Value of effective Nose-Hoover Hamiltonian",g18.8)') &
                      &bo%o%current%kin_energy + bo%o%current%pot_energy + &
                      &bo%o%current%hoover_xi**2*bo%o%current%hoover_mass/2.0_double + &
                      &3.0_double*bo%o%current%natoms*bo%o%current%hoover_s*bo%o%desired_temp*K_BOLTZMAN

                 call diary(bo%o%ave)
              case default
                 write(x_unit(diaryfile()),'(/,"Unknown md method:",i6)') bo%o%md_method
              end select wr

        end if

      end subroutine

      subroutine diary_bo_state(bo_state)
!doc$ subroutine diary(bo_state)
        type(bo_state_obj) :: bo_state
!       effects: Writes information about the current md state

!cod$

        if (i_access( diaryfile() )) then

           write(x_unit(diaryfile()),'(/,"MD time:",g18.8)') bo_state%time
           write(x_unit(diaryfile()),'(/,"Temperature     :",g18.8)') bo_state%temperature
           write(x_unit(diaryfile()),'(/,"Kinetic   Energy:",g18.8)') bo_state%kin_energy
           write(x_unit(diaryfile()),'(/,"Potential Energy:",g18.8)') bo_state%pot_energy
           write(x_unit(diaryfile()),'(/,"Total     Energy:",g18.8)') bo_state%tot_energy

        end if

      end subroutine

      subroutine diary_bo_ave(bo_ave)
!doc$ subroutine diary(bo_ave)
        type(bo_ave_obj) :: bo_ave
!       effects: Writes information about the current md averages

!cod$
        real(double) :: ave,sd,ave1,ave2,ave3,sd1,sd2,sd3
        real(double) :: press2kbar = 1.471076D+05
        integer :: i

        if (bo_ave%n_ave_steps .gt.0 .and. i_access( diaryfile() )) then

           write(x_unit(diaryfile()),'(/,"Number of Averaging Steps:",i5)') bo_ave%n_ave_steps
           ave = bo_ave%sum_temperature/bo_ave%n_ave_steps
           sd = std_dev(bo_ave%sum_temperature,bo_ave%sum_temperature_sq,bo_ave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Temperature     :",g18.8,"(",g12.4,")")') &
             & ave,sd
           ave = bo_ave%sum_kin_energy/bo_ave%n_ave_steps
           sd = std_dev(bo_ave%sum_kin_energy,bo_ave%sum_kin_e_sq,bo_ave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Kinetic   Energy:",g18.8,"(",g12.4,")")') &
             & ave,sd
           ave = bo_ave%sum_pot_energy/bo_ave%n_ave_steps
           sd = std_dev(bo_ave%sum_pot_energy,bo_ave%sum_pot_e_sq,bo_ave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Potential Energy:",g18.8,"(",g12.4,")")') &
             & ave,sd
           ave = bo_ave%sum_tot_energy/bo_ave%n_ave_steps
           sd = std_dev(bo_ave%sum_tot_energy,bo_ave%sum_tot_e_sq,bo_ave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Total     Energy:",g18.8,"(",g12.4,")")') &
             & ave,sd
           if (bo_ave%ave_pressure) then
              ave = bo_ave%sum_pressure/bo_ave%n_ave_steps
              sd = std_dev(bo_ave%sum_pressure,bo_ave%sum_pressure_sq,bo_ave%n_ave_steps)
              write(x_unit(diaryfile()), &
                   &'(/,"Average Pressure (kbar, excludes thermal):",g18.8,"(",g12.4,")")') &
                   & press2kbar*ave,press2kbar*sd
           end if
           if (bo_ave%ave_stress_tensor) then
              write(x_unit(diaryfile()), &
                 &'(/,"Average Stress Tensor (kbar, excludes thermal):")')
              do i = 1,3
                 ave1 = bo_ave%sum_stress_tensor(i,1)/bo_ave%n_ave_steps
                 sd1 = std_dev(bo_ave%sum_stress_tensor(i,1), &
                         &bo_ave%sum_stress_tensor_sq(i,1),bo_ave%n_ave_steps)
                 ave2 = bo_ave%sum_stress_tensor(i,2)/bo_ave%n_ave_steps
                 sd2 = std_dev(bo_ave%sum_stress_tensor(i,2), &
                         &bo_ave%sum_stress_tensor_sq(i,2),bo_ave%n_ave_steps)
                 ave3 = bo_ave%sum_stress_tensor(i,3)/bo_ave%n_ave_steps
                 sd3 = std_dev(bo_ave%sum_stress_tensor(i,3), &
                         &bo_ave%sum_stress_tensor_sq(i,3),bo_ave%n_ave_steps)
                 write(x_unit(diaryfile()),&
                      &'(5x,3(g15.5,"(",g10.2,")"))') &
                      &press2kbar*ave1,press2kbar*sd1, &
                      &press2kbar*ave2,press2kbar*sd2, &
                      &press2kbar*ave3,press2kbar*sd3
              end do
           end if
        end if

      end subroutine

      function std_dev(s_f,s_f_sq,n) result(sd)
        real(double), intent(in) :: s_f, s_f_sq
        integer, intent(in) :: n
        real(double) :: sd

        sd = (s_f_sq - s_f*s_f/n)/n
        sd = sqrt(sd)

      end function

      function bo_brancher(cfg) result(changed)
!doc$ function born_oppenheimer_dynamics(cfg)
!  initial crude brancher for choosing relaxation method etc.
        type(config_sc_obj) :: cfg
        logical :: changed

!cod$        

        type(bo_obj) :: bo
        
        call my(cfg)
        call my(born_oppenheimer(cfg),bo)
        
        branch: SELECT CASE (bo%o%md_method)
           CASE (NONE)
              changed = .false.
           CASE (NVE)
              changed = bo_velocity_verlet(cfg,bo)
           CASE (NVT_RESCALE)
              changed = bo_rescale(cfg,bo)
           CASE (NVT_ANDERSON)
              changed = bo_anderson(cfg,bo)
           CASE (NVT_HOOVER)
              changed = hoover_md(cfg,bo)
           CASE default
              if (error(.true.,"ERROR: md_method not defined")) call diary(bo)
              changed = .false.
           END SELECT branch

        call glean(thy(bo))
        call glean(thy(cfg))

      end function

      function bo_velocity_verlet(cfg,bo) result(changed)
!doc$ function bo_velocity_verlet(cfg,bo)
        type(config_sc_obj) :: cfg
        type(bo_obj) :: bo
        logical :: changed
! perform NVE molecular dynamics with the 'velocity Verlet' algorithm

!cod$

! define a set of objects that will be used for the modification of the config
        integer :: natoms,i
        integer :: status
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double) :: time_step

        call my(cfg)
        call my(bo)

        natoms = x_n_atoms(cfg)

        changed = .false.

        time_step = bo%o%time_step

        allocate( frcs(3,natoms), accel(3,natoms), STAT=status )
        if (error(status /= 0,"ERROR: allocate failure")) goto 100

        call x_cart_positions(cfg,bo%o%current%pos_cart)

        frcs = x_forces(cfg)
        do i = 1,3
          accel(i,:) = frcs(i,:)/bo%o%current%masses(:)
        end do

        ! start the molecular dynamics loop

        do
          call diary(cfg,bo%o%n_md_steps)
          if (bo%o%n_md_steps .ge. bo%o%max_md_steps) exit
          if (user_stop()) then
            call warn("WARNING: USER INITIATED STOP")
            exit
          end if
          bo%o%current%pos_cart = bo%o%current%pos_cart + &
                                  time_step*(bo%o%current%velocities + (time_step/2.0_double)*accel)
          bo%o%current%velocities = bo%o%current%velocities + (time_step/2.0_double)*accel
          call update_config(cfg,bo%o%current%pos_cart) ; if (error()) goto 100
          changed = .true.
          frcs = x_forces(cfg)
          do i = 1,3
            accel(i,:) = frcs(i,:)/bo%o%current%masses(:)
          end do
          bo%o%current%velocities = bo%o%current%velocities + (time_step/2.0_double)*accel
          bo%o%n_md_steps = bo%o%n_md_steps + 1
          call update_mdinfo_i(cfg,bo%o)
          call diary(bo%o%current)
          call diary(bo%o%ave)
        end do

100     call glean(thy(bo))
        call glean(thy(cfg))

        if (error("Exit born_oppenheimer_mod::bo_velocity_verlet")) continue

      end function

      function bo_rescale(cfg,bo) result(changed)
!doc$ function bo_rescale(cfg,bo)
        type(config_sc_obj) :: cfg
        type(bo_obj) :: bo
        logical :: changed
! perform NVE molecular dynamics with the 'velocity Verlet' algorithm

!cod$

! define a set of objects that will be used for the modification of the config
        integer :: i
        integer :: status
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double) :: time_step
        
        call my(cfg)
        call my(bo)
        
        changed = .false.

        time_step = bo%o%time_step
        
        allocate( frcs(3,bo%o%current%natoms), &
             & accel(3,bo%o%current%natoms), STAT=status )
        
        if (error(status /= 0,"ERROR: allocate failure")) goto 100

        call x_cart_positions(cfg,bo%o%current%pos_cart)

        frcs = x_forces(cfg)
        do i = 1,3
           accel(i,:) = frcs(i,:)/bo%o%current%masses(:)
        end do
        
        
! start the molecular dynamics loop

        do
           call diary(cfg,bo%o%n_md_steps)
           if (bo%o%n_md_steps .ge. bo%o%max_md_steps) exit
           if (user_stop()) then
             call warn("WARNING: USER INITIATED STOP")
             exit
           end if

           bo%o%current%pos_cart = bo%o%current%pos_cart + &
                         time_step*(bo%o%current%velocities + (time_step/2.0_double)*accel)

           bo%o%current%velocities = bo%o%current%velocities + (time_step/2.0_double)*accel

           call update_config(cfg,bo%o%current%pos_cart) ; if (error()) goto 100

           changed = .true.

           frcs = x_forces(cfg)
           do i = 1,3
              accel(i,:) = frcs(i,:)/bo%o%current%masses(:)
           end do

           bo%o%current%velocities = bo%o%current%velocities + (time_step/2.0_double)*accel

           if (mod(bo%o%n_md_steps,bo%o%temp_mod_freq) == 0) then
              call rescale_velocity(bo,bo%o%desired_temp)
           end if
           
           bo%o%n_md_steps = bo%o%n_md_steps + 1
           
           call update_mdinfo_i(cfg,bo%o)

           call diary(bo%o%current)

           call diary(bo%o%ave)

           
        end do

100     call glean(thy(bo))
        call glean(thy(cfg))

        if (error("Exit born_oppenheimer_mod::bo_rescale")) continue

      end function

      function bo_anderson(cfg,bo) result(changed)
!doc$ function bo_anderson(cfg,bo)
        type(config_sc_obj) :: cfg
        type(bo_obj) :: bo
        logical :: changed
! perform NVT molecular dynamics with the 'velocity Verlet' algorithm coupled with 
! the bo_anderson thermostat (coupling to stochastic bath)

!cod$

        integer :: i
        integer :: seed, status
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double), allocatable, dimension(:) :: test_numbers
        real(double) :: r1, time_step
        
        call my(cfg)
        call my(bo)
        
        changed = .false.

        time_step = bo%o%time_step
        
        allocate( frcs(3,bo%o%current%natoms), STAT=status)
        if (error(status /= 0,"ERROR: allocate failure - forces")) goto 100
        allocate( accel(3,bo%o%current%natoms), STAT=status)
        if (error(status /= 0,"ERROR: allocate failure - accel")) goto 100
        allocate( test_numbers(bo%o%current%natoms), STAT=status )
        if (error(status /= 0,"ERROR: allocate failure - test")) goto 100

        call x_cart_positions(cfg,bo%o%current%pos_cart)

        frcs = x_forces(cfg)
        
        do i = 1,3
           accel(i,:) = frcs(i,:)/bo%o%current%masses(:)
        end do

! start the molecular dynamics loop   
        do
           call diary(cfg,bo%o%n_md_steps)
           if (bo%o%n_md_steps .ge. bo%o%max_md_steps) exit
           if (user_stop()) then
             call warn("WARNING: USER INITIATED STOP")
             exit
           end if

           bo%o%current%pos_cart = bo%o%current%pos_cart + &
                         time_step*(bo%o%current%velocities + (time_step/2.0_double)*accel)

           bo%o%current%velocities = bo%o%current%velocities + (time_step/2.0_double)*accel

           call update_config(cfg,bo%o%current%pos_cart) ; if (error()) goto 100

           changed = .true.

           frcs = x_forces(cfg)
           do i = 1,3
              accel(i,:) = frcs(i,:)/bo%o%current%masses(:)
           end do

           bo%o%current%velocities = bo%o%current%velocities + (time_step/2.0_double)*accel

!          The following generates the same 'random' numbers at every pass through the loop
!          if (mpi_first(MOD_SCOPE)) then
!            seed = 101
!            do i = 1,53
!              r1 = random(seed)
!            end do
!            do i = 1,size(test_numbers)
!              test_numbers(i) = random(seed)
!            end do
!          end if

           if (mpi_first(MOD_SCOPE)) call random_number(test_numbers)
           call broadcast(MOD_SCOPE,test_numbers)
           do i = 1,bo%o%current%natoms
              if (test_numbers(i) < 1.0_double/bo%o%temp_mod_freq) then
                 bo%o%current%velocities(:,i) = sqrt(bo%o%desired_temp*K_BOLTZMAN/bo%o%current%masses(i))*rgauss(3)
              end if
           end do
           
           bo%o%n_md_steps = bo%o%n_md_steps + 1
           
           call update_mdinfo_i(cfg,bo%o)

           call diary(bo%o%current)

           call diary(bo%o%ave)
           
        end do

100     call glean(thy(bo))
        call glean(thy(cfg))

        if (error("Exit born_oppenheimer_mod::bo_anderson")) continue

      end function

      function hoover_md(cfg,bo) result(changed)
!doc$ function hoover_md(cfg,bo)
        type(config_sc_obj) :: cfg
        type(bo_obj) :: bo
        logical :: changed
! perform NVT molecular dynamics with the Hoover implementation of Nose Thermostat algorithm
! follows algorithm in "Understanding Molecular Simulation", by Daan Frenkel and Berend Smit

!cod$

        real(double), parameter :: SOLVE_TOL = 1.d-6
        integer, parameter :: MAX_SOLVE_TRIES = 12

        integer :: natoms,i
        integer :: status
        integer :: n_solve_tries
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double), allocatable, dimension(:,:) :: b,c,h
        real(double), allocatable, dimension(:,:) :: vel_prime
        real(double), allocatable, dimension(:,:) :: v1, v2
        real(double) :: a,d,hplus1,del_xi,xi1,xi2
        real(double) :: diff_ke
        real(double) :: test_val
        real(double) :: xi_prime
        real(double) :: time_step
        logical :: solved
        
        call my(cfg)
        call my(bo)
        
        natoms = x_n_atoms(cfg)

        changed = .false.

        time_step = bo%o%time_step
        
        allocate( frcs(3,natoms), accel(3,natoms), STAT=status )
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate( b(3,natoms), c(3,natoms), h(3,natoms), STAT=status )
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        allocate( v1(3,natoms), v2(3,natoms), vel_prime(3,natoms), STAT=status )
        if (error(status /= 0,"ERROR: allocate failure")) goto 100
        
        call x_cart_positions(cfg,bo%o%current%pos_cart)

        frcs = x_forces(cfg)
        
! start the molecular dynamics loop

        do
           call diary(cfg,bo%o%n_md_steps)
           if (bo%o%n_md_steps .ge. bo%o%max_md_steps) exit
           if (user_stop()) then
             call warn("WARNING: USER INITIATED STOP")
             exit
           end if

           do i = 1,3
              accel(i,:) = frcs(i,:)/bo%o%current%masses(:) - &
                   & bo%o%current%hoover_xi*bo%o%current%velocities(i,:)
           end do
        
           bo%o%current%pos_cart = bo%o%current%pos_cart + &
                         time_step*(bo%o%current%velocities + (time_step/2.0_double)*accel)

           vel_prime = bo%o%current%velocities + (time_step/2.0_double)*accel

           diff_ke = 0.0_double
           do i = 1,natoms
              diff_ke = diff_ke + &
                & bo%o%current%masses(i)*&
                & dot_product(bo%o%current%velocities(:,i),bo%o%current%velocities(:,i))
           end do
           diff_ke = diff_ke - 3.0_double*natoms*K_BOLTZMAN*bo%o%desired_temp
           
           bo%o%current%hoover_s = bo%o%current%hoover_s + &
                & bo%o%current%hoover_xi*time_step + &
                & diff_ke*time_step**2/(2.0_double*bo%o%current%hoover_mass)

           xi_prime = bo%o%current%hoover_xi + &
                & diff_ke*time_step/(2.0_double*bo%o%current%hoover_mass)

           call update_config(cfg,bo%o%current%pos_cart) ; if (error()) goto 100

           changed = .true.

           frcs = x_forces(cfg)

           solved = .false.
           n_solve_tries = 0

           v1 = vel_prime
           xi1 = xi_prime

           do while (.not.solved)

              xi2 = xi1
              v2 = v1
              
              do i = 1,3
                 accel(i,:) = frcs(i,:)/bo%o%current%masses(:) - &
                      & xi2*v2(i,:)
                 b(i,:) = bo%o%current%masses(:)*v2(i,:)*&
                      &time_step/bo%o%current%hoover_mass
              end do

              diff_ke = 0.0_double
              do i = 1,natoms
                 diff_ke = diff_ke + &
                      & bo%o%current%masses(i)*dot_product(v2(:,i),v2(:,i))
              end do
              diff_ke = diff_ke - 3.0_double*natoms*K_BOLTZMAN*bo%o%desired_temp

              h = vel_prime + accel*(time_step/2.0_double) - v2
              hplus1 = xi_prime - xi2 + &
                  & diff_ke*time_step/(2.0_double*bo%o%current%hoover_mass)
              
              c = -v2*time_step/2.0_double
              a = -1.0_double
              d = -xi2*time_step/2.0_double - 1.0_double

              del_xi = (hplus1*d - atom_dot_product(h,b))/&
                   &(-a*d + atom_dot_product(b,c))
              xi1 = xi2 + del_xi
              v1 = v2 + (-h - del_xi*c)/d

              test_val = max(maxval(abs(v1-v2)/v1),abs(xi1-xi2)/xi1)
              if (test_val < SOLVE_TOL) solved = .true.

              n_solve_tries = n_solve_tries + 1
              if (error(n_solve_tries > MAX_SOLVE_TRIES,"ERROR: solve failed")) goto 100

           end do
           
           bo%o%current%velocities = v1
           bo%o%current%hoover_xi = xi1

           bo%o%n_md_steps = bo%o%n_md_steps + 1
           
           call update_mdinfo_i(cfg,bo%o)

           call diary(bo)

        end do

100     deallocate(frcs,accel,b,c,h,v1,v2,vel_prime)
        
        call glean(thy(bo))
        call glean(thy(cfg))

        if (error("Exit born_oppenheimer_mod::hoover_md")) continue

      end function

      subroutine update_mdinfo_i(cfg,bor)
        type(config_sc_obj) :: cfg
        type(bo_rep) :: bor

        integer :: i,j
        real(double) :: p
        real(double), dimension(3,3) :: st

        call my(cfg)

        bor%current%natoms = x_n_atoms(cfg)

        bor%current%pot_energy = x_cell_energy(cfg)
        bor%current%kin_energy = 0.0_double
        do i = 1,bor%current%natoms
           bor%current%kin_energy = bor%current%kin_energy + &
             0.5_double*bor%current%masses(i)*dot_product(bor%current%velocities(:,i),bor%current%velocities(:,i))
        end do
        bor%current%temperature = compute_temperature_i(bor)
        bor%current%tot_energy = bor%current%pot_energy + bor%current%kin_energy
        bor%current%time = bor%n_md_steps * bor%time_step

        if (bor%n_md_steps > bor%n_skip_steps) then
           bor%ave%n_ave_steps = bor%ave%n_ave_steps + 1
           bor%ave%sum_temperature = bor%ave%sum_temperature + bor%current%temperature
           bor%ave%sum_temperature_sq = bor%ave%sum_temperature_sq + bor%current%temperature**2
           bor%ave%sum_kin_energy = bor%ave%sum_kin_energy + bor%current%kin_energy
           bor%ave%sum_kin_e_sq = bor%ave%sum_kin_e_sq + bor%current%kin_energy**2
           bor%ave%sum_pot_energy = bor%ave%sum_pot_energy + bor%current%pot_energy
           bor%ave%sum_pot_e_sq = bor%ave%sum_pot_e_sq + bor%current%pot_energy**2
           bor%ave%sum_tot_energy = bor%ave%sum_tot_energy + bor%current%tot_energy
           bor%ave%sum_tot_e_sq = bor%ave%sum_tot_e_sq + bor%current%tot_energy**2
           if (bor%ave%ave_pressure) then
              p = x_pressure(cfg)
              if (i_access( diaryfile() )) then
                write(x_unit(diaryfile()),'(/,"   current pressure:",f14.4,3x,"kbar")') p*147105.164_double
              end if
              bor%ave%sum_pressure = bor%ave%sum_pressure + p
              bor%ave%sum_pressure_sq = bor%ave%sum_pressure_sq + p**2
           end if
           if (bor%ave%ave_stress_tensor) then
              st = x_stress_tensor(cfg)
              if (i_access( diaryfile() )) then
                write(x_unit(diaryfile()),'(/,t4,"stress tensor:")')
                write(x_unit(diaryfile()),'(t8,3f12.4,3x,"kbar")') st(1,:)*147105.164_double
                write(x_unit(diaryfile()),'(t8,3f12.4)')           st(2,:)*147105.164_double
                write(x_unit(diaryfile()),'(t8,3f12.4)')           st(3,:)*147105.164_double
              end if
              do i = 1,3
                 do j = 1,3
                    bor%ave%sum_stress_tensor(i,j) = bor%ave%sum_stress_tensor(i,j) + st(i,j)
                    bor%ave%sum_stress_tensor_sq(i,j) = bor%ave%sum_stress_tensor_sq(i,j) + st(i,j)**2
                 end do
              end do
           end if
        end if
           
        if (i_access(bor%md_file)) then
           write(x_unit(bor%md_file),"(3g18.8)") bor%current%time, bor%current%temperature, bor%current%tot_energy
           write(x_unit(bor%md_file),"(3g18.8)") x_lattice_vector(cfg,1)
           write(x_unit(bor%md_file),"(3g18.8)") x_lattice_vector(cfg,2)
           write(x_unit(bor%md_file),"(3g18.8)") x_lattice_vector(cfg,3)
           do i = 1,bor%current%natoms
              write(x_unit(bor%md_file),"(a8,3g18.8,/,8x,3g18.8)") &
                   x_type(cfg,i), bor%current%pos_cart(:,i),bor%current%velocities(:,i)
           end do
           call flushbuf(bor%md_file)
        end if
              
        if (i_access(bor%new_vel_file)) then
           open(unit=x_unit(bor%new_vel_file),file=x_name(bor%new_vel_file))
           do i = 1,bor%current%natoms
              write(x_unit(bor%new_vel_file),"(3g18.8)") bor%current%velocities(:,i)
           end do
           close(unit=x_unit(bor%new_vel_file))
        end if

        call glean(thy(cfg))

      end subroutine
      
      subroutine init_velocities(vel,masses,temperature)
        real(double), dimension(:,:), intent(inout) :: vel
        real(double), dimension(:), intent(in) :: masses
        real(double), intent(in) :: temperature

        integer :: i, na, nd
        real(double) :: ekin
        real(double), dimension(size(vel,1)) :: pcm
        real(double) :: kbt

        nd = size(vel,1)
        na = size(vel,2)

        if (error(size(vel,2) /= size(masses),"ERROR: size incompatability")) goto 100

        do i = 1,na
           if (error(masses(i) <= 0.,"ERROR: mass <= 0")) goto 100
        end do
           
! here vel is really momentum
        if (temperature <= 0.) then
           vel = 0.0_double
        else
           kbt = K_BOLTZMAN*temperature
           do i = 1,na
              vel(:,i) = sqrt(kbt*masses(i))*rgauss(nd)
           end do
           pcm = sum(vel,2)/na
           do i = 1,na
              vel(:,i) = vel(:,i) - pcm
           end do
           ekin = 0.0_double
           do i = 1,na
              ekin = ekin + dot_product(vel(:,i),vel(:,i))/(2.0_double*masses(i))
           end do
           vel = vel*sqrt((1.50_double*(na-1)*kbt)/ekin)
! now convert momentums back to velocities
           do i = 1,na
              vel(:,i) = vel(:,i)/masses(i)
              if ( i_access(diaryfile()) ) then
                write(x_unit(diaryfile()),"('Initial velocity for atom',i6,' =',3g14.5)") i,vel(1,i),vel(2,i),vel(3,i)
              end if
           end do
        end if

100     if (error("Exit born_oppenheimer_mod::init_velocities")) continue

      end subroutine
      
      function rgauss(n1) result(grnd)
        integer, intent(in) :: n1
        real(double), dimension(n1) :: grnd

        integer :: i, seed
        real(double) :: r1
        real(double), dimension(n1) :: a1, a2

!       The following generates the same 'random' numbers every time rgauss is called
!       if (mpi_first(MOD_SCOPE)) then
!         seed = 37
!         do i = 1,55
!           r1 = random(seed)
!         end do
!         do i = 1,size(a1)
!           a1(i) = random(seed)
!         end do
!       end if
!       call broadcast(MOD_SCOPE,a1)
!       if (mpi_first(MOD_SCOPE)) then
!         seed = 39
!         do i = 1,57
!           r1 = random(seed)
!         end do
!         do i = 1,size(a2)
!           a2(i) = random(seed)
!         end do
!       end if
!       call broadcast(MOD_SCOPE,a2)

        if (mpi_first(MOD_SCOPE)) call random_number(a1)
        call broadcast(MOD_SCOPE,a1)
        if (mpi_first(MOD_SCOPE)) call random_number(a2)
        call broadcast(MOD_SCOPE,a2)
        grnd = sqrt(-2.0_double*log(a1))*cos(two_pi*a2)

      end function rgauss

      subroutine rescale_velocity(bo,new_temperature)
        type(bo_obj) :: bo
        real(double), intent(in) :: new_temperature

        real(double) :: ratio
        real(double) :: current_t

        call my(bo)

        current_t = compute_temperature_i(bo%o)
        
        if (current_t <= 0.0_double) then
           call init_velocities(bo%o%current%velocities,bo%o%current%masses,new_temperature)
           bo%o%current%temperature = new_temperature
           bo%o%current%kin_energy = 1.50_double*(bo%o%current%natoms-1)*K_BOLTZMAN
           bo%o%current%tot_energy = bo%o%current%pot_energy + bo%o%current%kin_energy
        else
           if (new_temperature <= 0.0_double) then
              bo%o%current%velocities = 0.0_double
              bo%o%current%kin_energy = 0.0_double
              bo%o%current%tot_energy = bo%o%current%pot_energy
              bo%o%current%temperature = 0.0_double
           else
              ratio = new_temperature/current_t
              bo%o%current%velocities = sqrt(ratio)*bo%o%current%velocities
              bo%o%current%temperature = new_temperature
              bo%o%current%kin_energy = ratio*bo%o%current%kin_energy
              bo%o%current%tot_energy = bo%o%current%pot_energy + bo%o%current%kin_energy
           end if
        end if

        call glean(thy(bo))

      end subroutine

      function compute_temperature_i(bor) result(temp)
        type(bo_rep) :: bor
        real(double) :: temp

        real(double) :: sum_ke
        integer :: i

        sum_ke = 0.0_double
        do i = 1,bor%current%natoms
           sum_ke = sum_ke + 0.50_double*bor%current%masses(i)* &
                & dot_product(bor%current%velocities(:,i),bor%current%velocities(:,i))
        end do

        select case (bor%md_method)
        case (NVE,NVT_RESCALE)
           temp = sum_ke/(1.50_double*(bor%current%natoms-1)*K_BOLTZMAN)
        case (NVT_ANDERSON,NVT_HOOVER)
           temp = sum_ke/(1.50_double*bor%current%natoms*K_BOLTZMAN)
        case default
           if (error(.true.,"ERROR: unknown md_method")) temp = 0.0_double
        end select

      end function

      end module

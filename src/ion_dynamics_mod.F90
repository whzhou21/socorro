!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module ion_dynamics_mod
!doc$ module ion_dynamics_mod

!     ion_dynamics_mod provides procedures for moving the ions.

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
      use config_td_mod
      use shortcuts_mod
      use grid_mod
      use gen_density_mod

!cod$
      implicit none
      private 

      real(double), parameter :: K_BOLTZMAN = 6.3333d-06
      real(double), parameter :: AMU_2_ELECTRON_MASS = 1.83615d+03
      real(double), parameter :: AMU_2_ARU = 911.4430774971
      real(double), parameter :: FS_2_ARU = 1/0.048377687

      integer, parameter :: NONE = 0
      integer, parameter :: NVE = 1
      integer, parameter :: NVT_RESCALE = 2
      integer, parameter :: NVT_ANDERSON = 3
      integer, parameter :: NVT_HOOVER = 4
! Nuclear position file format
      integer, parameter :: PDB = 1
      integer, parameter :: XYZ = 2

      integer, parameter         :: MOD_SCOPE = CONFIG

      type, public :: iondyn_ave_obj
         integer :: n_ave_steps
         real(double) :: sum_temperature,sum_temperature_sq
         real(double) :: sum_kin_energy,sum_kin_e_sq
         real(double) :: sum_pot_energy,sum_pot_e_sq
         real(double) :: sum_tot_energy,sum_tot_e_sq
         logical :: ave_pressure, ave_stress_tensor
         real(double) :: sum_pressure, sum_pressure_sq
         real(double), dimension(3,3) :: sum_stress_tensor, sum_stress_tensor_sq
      end type iondyn_ave_obj

      type, public :: iondyn_state_obj
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
      end type iondyn_state_obj

      type :: ion_dynamics_rep
         integer          :: ref
         type(ghost)      :: g
         integer          :: md_method
         integer          :: n_md_steps
         integer          :: n_skip_steps
         integer          :: max_md_steps
         integer          :: temp_mod_freq
         integer          :: generate_velocities
         integer          :: nucpos_file_type
         integer          :: nucpos_interval
         integer          :: mesh_file_type
         integer          :: density_interval
         real(double)     :: time_step
         real(double)     :: init_temp
         real(double)     :: desired_temp
         type(file_obj)   :: md_file
         type(file_obj)   :: init_vel_file
         type(file_obj)   :: new_vel_file
         type(iondyn_state_obj), pointer :: current
         type(iondyn_ave_obj), pointer   :: ave
      end type ion_dynamics_rep

      type, public :: ion_dynamics_obj
         private
         integer :: ref
         type(ion_dynamics_rep), pointer :: o
      end type ion_dynamics_obj

       
!doc$
      public :: ion_dynamics, move_atoms
!      public :: v_verlet,rescale_md,anderson_md
      public :: my, thy, glean, bequeath, assignment(=)
      public :: write_restart
      public :: diary
      public :: x_step_size, x_time, x_method, x_istep, x_kinetic_energy
      public :: x_pos_cart
public :: write_info
interface write_info
   module procedure write_info_ion_dynamics
end interface
      
!cod$

      interface ion_dynamics
         module procedure cons_tddft_tags
      end interface

      interface my
         module procedure my_ion_dynamics, my_new_ion_dynamics
      end interface
      interface thy
         module procedure thy_ion_dynamics
      end interface
      interface glean
         module procedure glean_ion_dynamics
      end interface
      interface bequeath
         module procedure bequeath_ion_dynamics
      end interface
      interface assignment(=)
         module procedure assign_ion_dynamics
      end interface

      interface diary
         module procedure diary_ion_dynamics, diary_iondyn_state, diary_mdave
      end interface
      interface write_restart
         module procedure write_restart_ion_dynamics
      end interface
      
!      interface run_ion_dynamics
!         module procedure ion_dynamics_brancher
!      end interface

      interface move_atoms
         module procedure ion_dynamics_move_atoms
      end interface

      interface x_step_size
         module procedure ion_dynamics_step_size
      end interface
      interface x_time
         module procedure ion_dynamics_time
      end interface
      interface x_method
         module procedure ion_dynamics_method
      end interface
      interface x_istep
         module procedure ion_dynamics_istep
      end interface
      interface x_kinetic_energy
         module procedure ion_dynamics_kinetic_energy
      end interface
      interface x_pos_cart
         module procedure ion_dynamics_pos_cart,ion_dynamics_pos_cart_ia
      end interface
      
      contains

      subroutine my_ion_dynamics(iondyn)
!doc$ subroutine my<type(ion_dynamics_obj)>(iondyn)
        type(ion_dynamics_obj) :: iondyn

!cod$
        iondyn%ref = iondyn%ref + 1
        iondyn%o%ref = iondyn%o%ref + 1
      end subroutine

      subroutine my_new_ion_dynamics(iondyni,iondyn)
!doc$ subroutine my<type(ion_dynamics_obj)>(iondyni,iondyn)
      type(ion_dynamics_obj) :: iondyni,iondyn

!cod$

        iondyn%ref = 1
        iondyn%o => iondyni%o
        iondyn%o%ref = iondyn%o%ref + 1
      end subroutine

      function thy_ion_dynamics(iondyn) result(thy_iondyn)
!doc$ function thy<type(ion_dynamics_obj)>(iondyn) result(thy_iondyn)
        type(ion_dynamics_obj) :: iondyn, thy_iondyn

!cod$
        iondyn%ref = iondyn%ref - 1
        iondyn%o%ref = iondyn%o%ref - 1
        thy_iondyn%ref = iondyn%ref
        thy_iondyn%o => iondyn%o
      end function thy_ion_dynamics

      subroutine glean_ion_dynamics(iondyn)
!doc$ subroutine glean<type(ion_dynamics_obj)>(iondyn)
        type(ion_dynamics_obj) :: iondyn

!cod$
        if (iondyn%o%ref < 1) then
           if (associated( iondyn%o%current%pos_cart )) deallocate( iondyn%o%current%pos_cart )
           if (associated( iondyn%o%current%velocities )) deallocate( iondyn%o%current%velocities )
           if (associated( iondyn%o%current%masses )) deallocate( iondyn%o%current%masses )
           deallocate( iondyn%o%current )
           deallocate( iondyn%o%ave )
           deallocate( iondyn%o )
        end if
      end subroutine glean_ion_dynamics

      subroutine bequeath_ion_dynamics(iondyn)
!doc$ subroutine bequeath<type(ion_dynamics_obj)>(iondyn)
        type(ion_dynamics_obj) :: iondyn
!cod$
        continue
      end subroutine bequeath_ion_dynamics

      subroutine assign_ion_dynamics(iondyn,iondyn2)
!doc$ subroutine assignment(=) (iondyn,iondyn2)
        type(ion_dynamics_obj), intent(inout) :: iondyn
        type(ion_dynamics_obj), intent(in) :: iondyn2

!cod$
        type(ion_dynamics_obj) :: iondynt
        call my(iondyn2)
        iondynt%o => iondyn%o
        iondyn%o%ref = iondyn%o%ref - iondyn%ref
        iondyn%o => iondyn2%o
        iondyn%o%ref = iondyn%o%ref + iondyn%ref
        call glean(iondynt)
        call glean(thy(iondyn2))
      end subroutine assign_ion_dynamics

      subroutine diary_ion_dynamics(iondyn)
!doc$ subroutine diary(iondyn)
        type(ion_dynamics_obj) :: iondyn
!       effects: Writes information about the optimization

!cod$

        call my(iondyn)

        if (i_access( diaryfile() )) then

           wr: select case(iondyn%o%md_method)
              case (NONE)
                 write(x_unit(diaryfile()),'(/,"No molecular dynamics requested")')

              case (NVE)
                 write(x_unit(diaryfile()),'(/,"NVE MD (velocity Verlet) with time step: ",g12.5)') iondyn%o%time_step
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') iondyn%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') iondyn%o%n_skip_steps
                 if (iondyn%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') iondyn%o%n_md_steps
                 end if
                 call diary(iondyn%o%current)
                 call diary(iondyn%o%ave)
              case (NVT_RESCALE)
                 write(x_unit(diaryfile()),'(/,"NVT MD (velocity rescaling) with time step: ",g12.5)') iondyn%o%time_step
                 write(x_unit(diaryfile()),'(/," velocity rescaled every",i4," time steps")') iondyn%o%temp_mod_freq
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') iondyn%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') iondyn%o%n_skip_steps
                 if (iondyn%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') iondyn%o%n_md_steps
                 end if
                 call diary(iondyn%o%current)
                 call diary(iondyn%o%ave)
              case (NVT_ANDERSON)
                 write(x_unit(diaryfile()), &
                     &'(/,"NVT MD (stochastic velocities ala Anderson) with time step: ",g12.5)') &
                     &iondyn%o%time_step
                 write(x_unit(diaryfile()), &
                     &'(/," atom velocity rescaled at each time step with probability ", g12.5)') &
                     &1.0_double/iondyn%o%temp_mod_freq
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') iondyn%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') iondyn%o%n_skip_steps
                 if (iondyn%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') iondyn%o%n_md_steps
                 end if
                 call diary(iondyn%o%current)
                 call diary(iondyn%o%ave)
              case (NVT_HOOVER)
                 write(x_unit(diaryfile()), &
                     &'(/,"NVT MD (Hoover implementation of Nose thermostat) with time step: ",g12.5)') &
                     & iondyn%o%time_step
                 write(x_unit(diaryfile()),'(/," fictitious temperature-coupling mass ", g12.5)') iondyn%o%current%hoover_mass
                 write(x_unit(diaryfile()),'("Maximum number of md steps:",i4)') iondyn%o%max_md_steps
                 write(x_unit(diaryfile()),'("Number of equilibration md steps:",i4)') iondyn%o%n_skip_steps
                 if (iondyn%o%n_md_steps .gt. 0 ) then
                    write(x_unit(diaryfile()),'("Number of md steps taken:",i4)') iondyn%o%n_md_steps
                 end if
                 call diary(iondyn%o%current)
                 write(x_unit(diaryfile()),'("Value of effective Nose-Hoover Hamiltonian",g18.8)') &
                      &iondyn%o%current%kin_energy + iondyn%o%current%pot_energy + &
                      &iondyn%o%current%hoover_xi**2*iondyn%o%current%hoover_mass/2.0_double + &
                      &3.0_double*iondyn%o%current%natoms*iondyn%o%current%hoover_s*iondyn%o%desired_temp*K_BOLTZMAN

                 call diary(iondyn%o%ave)
              case default
                 write(x_unit(diaryfile()),'(/,"Unknown md method:",i6)') iondyn%o%md_method
              end select wr

        end if

        call glean(thy(iondyn))

      end subroutine diary_ion_dynamics

      subroutine diary_iondyn_state(mdst)
!doc$ subroutine diary(mdst)
        type(iondyn_state_obj), intent(in) :: mdst
!       effects: Writes information about the current md state

!cod$

        if (i_access( diaryfile() )) then

           write(x_unit(diaryfile()),'(/,"MD time:",g18.8)') mdst%time
           write(x_unit(diaryfile()),'(/,"Temperature     :",g18.8)') mdst%temperature
           write(x_unit(diaryfile()),'(/,"Kinetic   Energy:",g18.8)') mdst%kin_energy
           write(x_unit(diaryfile()),'(/,"Potential Energy:",g18.8)') mdst%pot_energy
           write(x_unit(diaryfile()),'(/,"Total     Energy:",g18.8)') mdst%tot_energy

        end if

      end subroutine diary_iondyn_state

      subroutine diary_mdave(mdave)
!doc$ subroutine diary(mdave)
        type(iondyn_ave_obj), intent(in) :: mdave
!       effects: Writes information about the current md averages

!cod$
        real(double), parameter :: press2kbar = 1.471076E+05_double
        real(double) :: ave, sd, ave1, ave2, ave3, sd1, sd2, sd3
        integer :: i, j

        if (mdave%n_ave_steps .gt.0 .and. i_access( diaryfile() )) then

           write(x_unit(diaryfile()),'(/,"Number of Averaging Steps:",i5)') mdave%n_ave_steps
           ave = mdave%sum_temperature/mdave%n_ave_steps
           sd = std_dev(mdave%sum_temperature,mdave%sum_temperature_sq,mdave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Temperature     :",g18.8,"(",g12.4,")")') &
             & ave,sd
           ave = mdave%sum_kin_energy/mdave%n_ave_steps
           sd = std_dev(mdave%sum_kin_energy,mdave%sum_kin_e_sq,mdave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Kinetic   Energy:",g18.8,"(",g12.4,")")') &
             & ave,sd
           ave = mdave%sum_pot_energy/mdave%n_ave_steps
           sd = std_dev(mdave%sum_pot_energy,mdave%sum_pot_e_sq,mdave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Potential Energy:",g18.8,"(",g12.4,")")') &
             & ave,sd
           ave = mdave%sum_tot_energy/mdave%n_ave_steps
           sd = std_dev(mdave%sum_tot_energy,mdave%sum_tot_e_sq,mdave%n_ave_steps)
           write(x_unit(diaryfile()), &
             &'(/,"Average Total     Energy:",g18.8,"(",g12.4,")")') &
             & ave,sd
           if (mdave%ave_pressure) then
              ave = mdave%sum_pressure/mdave%n_ave_steps
              sd = std_dev(mdave%sum_pressure,mdave%sum_pressure_sq,mdave%n_ave_steps)
              write(x_unit(diaryfile()), &
                   &'(/,"Average Pressure (kbar, excludes thermal):",g18.8,"(",g12.4,")")') &
                   & press2kbar*ave, press2kbar*sd
           end if
           if (mdave%ave_stress_tensor) then
              write(x_unit(diaryfile()), &
                 &'(/,"Average Stress Tensor (kbar, excludes thermal):")')
              do i = 1,3
                 ave1 = mdave%sum_stress_tensor(i,1)/mdave%n_ave_steps
                 sd1 = std_dev(mdave%sum_stress_tensor(i,1), &
                         &mdave%sum_stress_tensor_sq(i,1),mdave%n_ave_steps)
                 ave2 = mdave%sum_stress_tensor(i,2)/mdave%n_ave_steps
                 sd2 = std_dev(mdave%sum_stress_tensor(i,2), &
                         &mdave%sum_stress_tensor_sq(i,2),mdave%n_ave_steps)
                 ave3 = mdave%sum_stress_tensor(i,3)/mdave%n_ave_steps
                 sd3 = std_dev(mdave%sum_stress_tensor(i,3), &
                         &mdave%sum_stress_tensor_sq(i,3),mdave%n_ave_steps)
                 write(x_unit(diaryfile()),&
                      &'(5x,3(g15.5,"(",g10.2,")"))') &
                      & press2kbar*ave1, press2kbar*sd1, &
                      & press2kbar*ave2, press2kbar*sd2, &
                      & press2kbar*ave3, press2kbar*sd3
              end do
           end if
        end if

      end subroutine diary_mdave

      function std_dev(s_f,s_f_sq,n) result(sd)
        real(double), intent(in) :: s_f, s_f_sq
        integer, intent(in) :: n
        real(double) :: sd

        sd = (s_f_sq - s_f*s_f/n)/n
        sd = sqrt(sd)

      end function std_dev


      subroutine update_mdinfo_i(cfg,iondynr)
        type(config_td_obj), intent(inout) :: cfg
        type(ion_dynamics_rep), intent(inout) :: iondynr
        ! Local Vars ---------------------------------------
        integer            :: i,j
        real(double)       :: p
        real(double)       :: st(3,3)
        type(grid_obj)     :: g
        character(line_len) :: filename
        logical            :: debug
        !===================================================

        debug = .false.

        if (debug) call warn('update_mdinfo::starting')

        call my(cfg)

        iondynr%current%natoms = x_n_atoms(cfg)

        iondynr%current%pot_energy = x_cell_energy(cfg)
        iondynr%current%kin_energy = 0.0_double
        do i = 1,iondynr%current%natoms
           iondynr%current%kin_energy = iondynr%current%kin_energy + &
             0.5_double*iondynr%current%masses(i)*dot_product(iondynr%current%velocities(:,i),iondynr%current%velocities(:,i))
        end do
        iondynr%current%temperature = compute_temperature_i(iondynr)
        iondynr%current%tot_energy = iondynr%current%pot_energy + iondynr%current%kin_energy
        iondynr%current%time = iondynr%n_md_steps * iondynr%time_step

        if (iondynr%n_md_steps > iondynr%n_skip_steps) then
           iondynr%ave%n_ave_steps = iondynr%ave%n_ave_steps + 1
           iondynr%ave%sum_temperature = iondynr%ave%sum_temperature + iondynr%current%temperature
           iondynr%ave%sum_temperature_sq = iondynr%ave%sum_temperature_sq + iondynr%current%temperature**2
           iondynr%ave%sum_kin_energy = iondynr%ave%sum_kin_energy + iondynr%current%kin_energy
           iondynr%ave%sum_kin_e_sq = iondynr%ave%sum_kin_e_sq + iondynr%current%kin_energy**2
           iondynr%ave%sum_pot_energy = iondynr%ave%sum_pot_energy + iondynr%current%pot_energy
           iondynr%ave%sum_pot_e_sq = iondynr%ave%sum_pot_e_sq + iondynr%current%pot_energy**2
           iondynr%ave%sum_tot_energy = iondynr%ave%sum_tot_energy + iondynr%current%tot_energy
           iondynr%ave%sum_tot_e_sq = iondynr%ave%sum_tot_e_sq + iondynr%current%tot_energy**2
           if (iondynr%ave%ave_pressure) then
              p = x_pressure(cfg)
              if (i_access( diaryfile() )) then
                write(x_unit(diaryfile()),'(/,"   current pressure:",f14.4,3x,"kbar")') p*147105.164_double
              end if
              iondynr%ave%sum_pressure = iondynr%ave%sum_pressure + p
              iondynr%ave%sum_pressure_sq = iondynr%ave%sum_pressure_sq + p**2
           end if
           if (iondynr%ave%ave_stress_tensor) then
              st = x_stress_tensor(cfg)
              if (i_access( diaryfile() )) then
                write(x_unit(diaryfile()),'(/,t4,"stress tensor:")')
                write(x_unit(diaryfile()),'(t8,3f12.4,3x,"kbar")') st(1,:)*147105.164_double
                write(x_unit(diaryfile()),'(t8,3f12.4)')           st(2,:)*147105.164_double
                write(x_unit(diaryfile()),'(t8,3f12.4)')           st(3,:)*147105.164_double
              end if
              do i = 1,3
                 do j = 1,3
                    iondynr%ave%sum_stress_tensor(i,j) = iondynr%ave%sum_stress_tensor(i,j) + st(i,j)
                    iondynr%ave%sum_stress_tensor_sq(i,j) = iondynr%ave%sum_stress_tensor_sq(i,j) + st(i,j)**2
                 end do
              end do
           end if
        end if
           
        if (i_access(iondynr%md_file)) then
           write(x_unit(iondynr%md_file),"(3g18.8)") &
                iondynr%current%time, iondynr%current%temperature, iondynr%current%tot_energy
           write(x_unit(iondynr%md_file),"(3g18.8)") x_lattice_vector(cfg,1)
           write(x_unit(iondynr%md_file),"(3g18.8)") x_lattice_vector(cfg,2)
           write(x_unit(iondynr%md_file),"(3g18.8)") x_lattice_vector(cfg,3)
           do i = 1,iondynr%current%natoms
              write(x_unit(iondynr%md_file),"(a8,3g18.8,/,8x,3g18.8)") &
                   x_type(cfg,i), iondynr%current%pos_cart(:,i),iondynr%current%velocities(:,i)
           end do
           call flushbuf(iondynr%md_file)
        end if
              
        if (i_access(iondynr%new_vel_file)) then
           open(unit=x_unit(iondynr%new_vel_file),file=x_name(iondynr%new_vel_file))
           do i = 1,iondynr%current%natoms
              write(x_unit(iondynr%new_vel_file),"(3g18.8)") iondynr%current%velocities(:,i)
           end do
           close(unit=x_unit(iondynr%new_vel_file))
        end if

        if (debug) call warn('update_mdinfo::checking to see if the nuclear positions need to be written out')

        !** Write the nuclear positions out to file if indicated by argvf parameters
        if (0 < iondynr%nucpos_interval) then
           if (mod(iondynr%n_md_steps,iondynr%nucpos_interval)==0) then
              select case (iondynr%nucpos_file_type)
              case (PDB)
                 call write_pdb_file_i(cfg,iondynr)
              case (XYZ)
                 call write_xyz_file_i(cfg,iondynr)
              end select
           end if
        end if
        
!        if (debug) call warn('update_mdinfo::writing density out if indicating by argvf')
!        !** Write out the density to file if indicated by argvf parameters.
!        if (0 < iondynr%density_interval) then
!           if (mod(iondynr%n_md_steps,iondynr%density_interval) == 0) then
!              call my(x_grid_density(x_density(x_electrons(cfg))),g)
!              if (error()) goto 100
!              if (iondynr%mesh_file_type == AMIRA_REGULAR) then
!                 call get_density_am_filename_i(filename,iondynr%n_md_steps)
!                 call write_to_file(g,filename,AMIRA_REGULAR,RS_KIND) ; if (error()) goto 100
!              else if (iondynr%mesh_file_type == VTK) then
!                 call get_density_vtk_filename_i(filename,iondynr%n_md_steps)
!                 call write_to_file(g,filename,VTK,RS_KIND) ; if (error()) goto 100
!              end if
!              call glean(thy(g))
!           end if
!        end if

        if (debug) call warn('update_mdinfo::gleaning')

100     call glean(thy(cfg))

        if (debug) call warn('update_mdinfo::exiting')


      end subroutine update_mdinfo_i
      
      subroutine init_velocities(vel,masses,temperature)
        real(double), dimension(:,:), intent(inout) :: vel
        real(double), dimension(:) , intent(in) :: masses
        real(double), intent(in) :: temperature

        integer :: i,na,nd
        real(double) :: ekin
        real(double), dimension(size(vel,1)) :: pcm
        real(double) :: kbt

        nd = size(vel,1)
        na = size(vel,2)

        if (error(size(vel,2) /= size(masses),"init_velocities: size incompatability")) goto 100

        do i = 1,na
           if (error(masses(i) <= 0.,"init_velocity: mass <= 0")) goto 100
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
           vel = vel*sqrt((1.5_double*(na-1)*kbt)/ekin)
! now convert momentums back to velocities
           do i = 1,na
              vel(:,i) = vel(:,i)/masses(i)
           end do
        end if

100     if (error("Exit ion_dynamics_mod:init_velocities")) continue

      end subroutine init_velocities
      
      function rgauss(n1) result(grnd)
        integer, intent(in) :: n1
        real(double),dimension(n1) :: grnd

        real(double), dimension(n1) :: a1,a2
        real(double) :: twopi = 6.28318530717958647690_double

        if (mpi_isroot(MOD_SCOPE)) call random_number(a1)
        call broadcast(MOD_SCOPE,a1)
        if (mpi_isroot(MOD_SCOPE)) call random_number(a2)
        call broadcast(MOD_SCOPE,a2)
        grnd = sqrt(-2.*log(a1))*cos(twopi*a2)
      end function rgauss

      subroutine rescale_velocity(iondyn,new_temperature)
        type(ion_dynamics_obj), intent(inout) :: iondyn
        real(double), intent(in) :: new_temperature

        real(double) :: ratio
        real(double) :: current_t
        integer :: i

        call my(iondyn)

        current_t = compute_temperature_i(iondyn%o)
        
        if (current_t <= 0.0_double) then
           call init_velocities(iondyn%o%current%velocities,iondyn%o%current%masses, &
                &new_temperature)
           iondyn%o%current%temperature = new_temperature
           iondyn%o%current%kin_energy = 1.5_double*(iondyn%o%current%natoms-1)*K_BOLTZMAN
           iondyn%o%current%tot_energy = &
                &iondyn%o%current%pot_energy + iondyn%o%current%kin_energy
        else
           if (new_temperature <= 0.0_double) then
              iondyn%o%current%velocities = 0.0_double
              iondyn%o%current%kin_energy = 0.0_double
              iondyn%o%current%tot_energy = iondyn%o%current%pot_energy
              iondyn%o%current%temperature = 0.0_double
           else
              ratio = new_temperature/current_t
              iondyn%o%current%velocities = sqrt(ratio)*iondyn%o%current%velocities
              iondyn%o%current%temperature = new_temperature
              iondyn%o%current%kin_energy = ratio*iondyn%o%current%kin_energy
              iondyn%o%current%tot_energy = iondyn%o%current%pot_energy + &
                   & iondyn%o%current%kin_energy
           end if
        end if

        call glean(thy(iondyn))

      end subroutine rescale_velocity

      function compute_temperature_i(iondynr) result(temp)
        type(ion_dynamics_rep), intent(in) :: iondynr
        real(double) :: temp

        real(double) :: sum_ke
        integer :: i

        sum_ke = 0.0_double
        do i = 1,iondynr%current%natoms
           sum_ke = sum_ke + 0.5_double*iondynr%current%masses(i)*&
                &dot_product(iondynr%current%velocities(:,i),iondynr%current%velocities(:,i))
        end do

        select case (iondynr%md_method)
        case (NVE,NVT_RESCALE)
           temp = sum_ke/(1.5_double*(iondynr%current%natoms-1)*K_BOLTZMAN)
        case (NVT_ANDERSON,NVT_HOOVER)
           temp = sum_ke/(1.5_double*iondynr%current%natoms*K_BOLTZMAN)
        case default
           if (error(.true.,"compute_temperature: unknown md_method")) temp = 0.0_double
        end select

      end function compute_temperature_i


      subroutine write_restart_ion_dynamics(iondyn,nrestf)
!doc$ subroutine write_restart(iondyn,nrestf)
        type(ion_dynamics_obj) :: iondyn
        type(tagio_obj), intent(inout) :: nrestf
!       modifies: nrestf
!       effects: Writes iondyn restart information to nrestf.
!       errors: Passes errors.

!cod$
        integer(long) :: dsize, ios, ndata, s4
        integer :: ia, na

        call my(iondyn)
        call my(nrestf)

        if (i_access(nrestf) .and. iondyn%o%md_method /= NONE) then
           call startblock(nrestf,"ION_DYNAMICS")

           call writetag(nrestf,"N_MD_STEPS")
           s4 = iondyn%o%n_md_steps ; dsize = sizeof_long ; ndata = 1
           call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)

           call writetag(nrestf,"CURRENT")
           dsize = sizeof_double ; ndata = 1
           call writef(iondyn%o%current%time,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%current%kin_energy,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%current%pot_energy,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%current%tot_energy,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%current%temperature,dsize,ndata,x_tagfd(nrestf),ios)

           na = size(iondyn%o%current%pos_cart,1)
           ndata = size(iondyn%o%current%pos_cart,2)
           do ia=1,na
              call writef(iondyn%o%current%pos_cart(ia,:),dsize,ndata,x_tagfd(nrestf),ios)
           end do

           na = size(iondyn%o%current%velocities,1)
           dsize = sizeof_double ; ndata = size(iondyn%o%current%velocities,2)
           do ia=1,na
              call writef(iondyn%o%current%velocities(ia,:),dsize,ndata,x_tagfd(nrestf),ios)
           end do

           dsize = sizeof_double ; ndata = 1
           call writef(iondyn%o%current%hoover_xi,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%current%hoover_s,dsize,ndata,x_tagfd(nrestf),ios)

           call writetag(nrestf,"AVE")
           s4 = iondyn%o%ave%n_ave_steps ; dsize = sizeof_long ; ndata = 1
           call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)

           dsize = sizeof_double ; ndata = 1
           call writef(iondyn%o%ave%sum_temperature,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_temperature_sq,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_kin_energy,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_kin_e_sq,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_pot_energy,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_pot_e_sq,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_tot_energy,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_tot_e_sq,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_pressure,dsize,ndata,x_tagfd(nrestf),ios)
           call writef(iondyn%o%ave%sum_pressure_sq,dsize,ndata,x_tagfd(nrestf),ios)

           ndata = 3
           do ia=1,3
              call writef(iondyn%o%ave%sum_stress_tensor(ia,:),dsize,ndata,x_tagfd(nrestf),ios)
           end do

           do ia=1,3
              call writef(iondyn%o%ave%sum_stress_tensor_sq(ia,:),dsize,ndata,x_tagfd(nrestf),ios)
           end do

           call endblock(nrestf) ; if (error()) goto 100
        end if

100     call glean(thy(iondyn))
        call glean(thy(nrestf))

        if (error("Exit ion_dynamics_mod:: write_restart_ion_dynamics")) continue

      end subroutine


      subroutine read_restart_i(iondynr,restf)
        type(ion_dynamics_rep) :: iondynr
        type(tagio_obj) :: restf

        ! Local Vars
        integer(long)    :: dsize, ios, ndata, s4
        integer          :: i, na
        character(1)     :: tios


        call my(restf)

        if (i_access(restf)) tios = findfirsttag(restf,"ION_DYNAMICS")
        if (i_comm(restf)) call broadcast(MOD_SCOPE,tios)
        if (tios /= TAG_START_BLOCK) then
           call warn("WARNING: ION_DYNAMICS block was not found - reverting to a standard construction")
           goto 100
        end if
        
        if (i_access(restf)) call openblock(restf)
        
        !** Read in the number of md steps
        if (i_access(restf)) then
           tios = findfirsttag(restf,"N_MD_STEPS")
           dsize = sizeof_long ; ndata = 1
           call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
        end if
        if (i_comm(restf)) call broadcast(MOD_SCOPE,s4)
        iondynr%n_md_steps = s4
        
        !** Read in the info in the iondynr%current data structure
        if (i_access(restf)) then
           tios = findfirsttag(restf,"CURRENT")
           dsize = sizeof_double ; ndata = 1
           call readf(iondynr%current%time,       dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%current%kin_energy, dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%current%pot_energy, dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%current%tot_energy, dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%current%temperature,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           
           na = size(iondynr%current%pos_cart,1)
           ndata = size(iondynr%current%pos_cart,2)
           do i=1,na
              call readf(iondynr%current%pos_cart(i,:),dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           end do
           
           ndata = size(iondynr%current%velocities,2)
           do i=1,na
              call readf(iondynr%current%velocities(i,:),dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           end do
           
           ndata = 1
           call readf(iondynr%current%hoover_xi,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%current%hoover_s, dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
        end if
        if (i_comm(restf)) then
           call broadcast(MOD_SCOPE,iondynr%current%time)
           call broadcast(MOD_SCOPE,iondynr%current%kin_energy)
           call broadcast(MOD_SCOPE,iondynr%current%pot_energy)
           call broadcast(MOD_SCOPE,iondynr%current%tot_energy)
           call broadcast(MOD_SCOPE,iondynr%current%temperature)
           call broadcast(MOD_SCOPE,iondynr%current%pos_cart)
           call broadcast(MOD_SCOPE,iondynr%current%velocities)
           call broadcast(MOD_SCOPE,iondynr%current%hoover_xi)
           call broadcast(MOD_SCOPE,iondynr%current%hoover_s)
        end if
        
        
        !** Read in the info in the iondynr%ave data structure
        if (i_access(restf)) then
           tios = findfirsttag(restf,"AVE")
           dsize = sizeof_long ; ndata = 1
           call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
        end if
        iondynr%ave%n_ave_steps = s4
        if (i_comm(restf)) call broadcast(MOD_SCOPE,iondynr%ave%n_ave_steps)
        
        if (i_access(restf)) then
           
           dsize = sizeof_double ; ndata = 1
           call readf(iondynr%ave%sum_temperature,   dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_temperature_sq,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_kin_energy,    dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_kin_e_sq,      dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_pot_energy,    dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_pot_e_sq,      dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_tot_energy,    dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_tot_e_sq,      dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_pressure,      dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           call readf(iondynr%ave%sum_pressure_sq,   dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           
           ndata = 3
           do i=1,3
              call readf(iondynr%ave%sum_stress_tensor(i,:),dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           end do
           
           ndata = 3
           do i=1,3
              call readf(iondynr%ave%sum_stress_tensor_sq(i,:),dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           end do
           
        end if
        if (i_comm(restf)) then
           call broadcast(MOD_SCOPE,iondynr%ave%sum_temperature)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_temperature_sq)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_kin_energy)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_kin_e_sq)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_pot_energy)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_pot_e_sq)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_tot_energy)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_tot_e_sq)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_stress_tensor)
           call broadcast(MOD_SCOPE,iondynr%ave%sum_stress_tensor_sq)
        end if
        
        if (i_access(restf)) call closeblock(restf)
        
100     call glean(thy(restf))
        
      end subroutine

      function ion_dynamics_step_size(iondyn) result(step_size)
!doc$ function x_step_size(iondyn) result(step_size)
        type(ion_dynamics_obj) :: iondyn
        real(double) :: step_size

!     effects: Returns the time interval (in ARU's) between md steps.
!cod$

        call my(iondyn)
        step_size = iondyn%o%time_step
        call glean(thy(iondyn))

      end function 

      function ion_dynamics_time(iondyn) result(time)
!doc$ function x_time(iondyn) result(time)
        type(ion_dynamics_obj) :: iondyn
        real(double) :: time

!     effects: Returns the time (in ARU's) of the last md step
!cod$

        call my(iondyn)
        time = iondyn%o%current%time
        call glean(thy(iondyn))

      end function 

      function ion_dynamics_method(iondyn) result(method)
!doc$ function x_method(iondyn) result(method)
        type(ion_dynamics_obj) :: iondyn
        integer :: method

!     effects: Returns an integer representing the method used for this calc.
!cod$

        call my(iondyn)
        method = iondyn%o%md_method
        call glean(thy(iondyn))

      end function 

      function ion_dynamics_istep(iondyn) result(istep)
!doc$ function x_istep(iondyn) result(istep)
        type(ion_dynamics_obj) :: iondyn
        integer :: istep

!     effects: Returns the current md step
!cod$

        call my(iondyn)
        istep = iondyn%o%n_md_steps
        call glean(thy(iondyn))

      end function 

      function ion_dynamics_kinetic_energy(iondyn) result(ke)
!doc$ function x_kinetic_energy(iondyn) result(ke)
        type(ion_dynamics_obj) :: iondyn
        real(double) :: ke

!     effects: Returns the kinetic energy of the nuclei
!cod$

        call my(iondyn)
        ke = iondyn%o%current%kin_energy
        call glean(thy(iondyn))

      end function 

      function ion_dynamics_pos_cart(iondyn) result(pos_cart)
!doc$ function x_pos_cart(iondyn) result(pos_cart)
        type(ion_dynamics_obj) :: iondyn
        real(double), dimension(3,iondyn%o%current%natoms) :: pos_cart
!       effects: Returns the atom positions in the cartesian representation.

!cod$
        integer :: ia
        call my(iondyn)
        do ia = 1,iondyn%o%current%natoms
          pos_cart(:,ia) = iondyn%o%current%pos_cart(:,ia)
        end do
        call glean(thy(iondyn))
      end function

      function ion_dynamics_pos_cart_ia(iondyn,ia) result(pos_cart)
!doc$ function x_pos_cart(iondyn) result(pos_cart)
        type(ion_dynamics_obj) :: iondyn
        integer :: ia
        real(double) :: pos_cart(3)
!       effects: Returns the atom positions in the cartesian representation.

!cod$
        call my(iondyn)
        pos_cart = iondyn%o%current%pos_cart(:,ia)
        call glean(thy(iondyn))
      end function


!************************************************************************************
!************************************************************************************
!*******    Below are the routines directly related to the TDDFT-based MD     *******
!****                                                                            ****
!**  Note that there are many similarities but that only a single move takes       **
!**  place in a given call.  The electronic system propagates by time evolving     **
!**  the Kohn-Sham wavefunctions explicitly                                        **
!****                                                                            ****
!************************************************************************************
      
      function cons_tddft_tags(cfg, time_step, restf) result(iondyn)
!doc$ function ion_dynamics(cfg) result(iondyn)
        type(config_td_obj)       :: cfg
        real(double)              :: time_step
        type(tagio_obj), optional :: restf
        type(ion_dynamics_obj)    :: iondyn
!       effects: sets the structural optimization parameters

!cod$
        ! Local Vars -----------------------------------------------
        integer          :: i,j,na
        integer          :: status
        integer          :: ios
        logical          :: exist_file
        logical          :: debug
        !===========================================================

        debug = .false.

        if (debug) call warn('cons_tddft_ion_dynamics::starting')

        call my(cfg)

        iondyn%ref = 0
        na = x_n_atoms(cfg)
        allocate( iondyn%o, STAT=status)
        if (error(status /= 0,"cons_tddft_ion_dynamics: allocate failure")) goto 100
        allocate(iondyn%o%ave, STAT=status)
        if (error(status /= 0,"cons_tddft_ion_dynamics: allocate failure")) goto 100
        allocate(iondyn%o%current, STAT=status)
        if (error(status /= 0,"cons_tddft_ion_dynamics: allocate failure")) goto 100
        allocate(iondyn%o%current%pos_cart(3,na), STAT=status)
        if (error(status /= 0,"cons_tddft_ion_dynamics: allocate failure")) goto 100
        allocate(iondyn%o%current%velocities(3,na), STAT=status)
        if (error(status /= 0,"cons_tddft_ion_dynamics: allocate failure")) goto 100
        allocate(iondyn%o%current%masses(na), STAT=status )
        if (error(status /= 0,"cons_tddft_ion_dynamics: allocate failure")) goto 100

        iondyn%o%g = x_ghost()
        iondyn%o%ref = 0

        iondyn%o%n_md_steps = 0

        if (debug) call warn('cons_tddft_ion_dynamics::calling read_params_i')

        if (error(time_step < 0,"cons_tddft_ion_dynamics: time_step < 0")) goto 100
        
        iondyn%o%time_step = time_step
        
        call read_params_i(iondyn%o,cfg) ; if (error()) goto 100

        if (debug) call warn('cons_tddft_ion_dynamics::setting constants')

        iondyn%o%current%hoover_s = 0.0_double
        iondyn%o%current%hoover_xi = 0.0_double
        
        iondyn%o%ave%n_ave_steps = 0
        iondyn%o%ave%sum_temperature = 0.0_double
        iondyn%o%ave%sum_temperature_sq = 0.0_double
        iondyn%o%ave%sum_kin_energy = 0.0_double
        iondyn%o%ave%sum_kin_e_sq = 0.0_double
        iondyn%o%ave%sum_pot_energy = 0.0_double
        iondyn%o%ave%sum_pot_e_sq = 0.0_double
        iondyn%o%ave%sum_tot_energy = 0.0_double
        iondyn%o%ave%sum_tot_e_sq = 0.0_double
        iondyn%o%ave%sum_pressure = 0.0_double
        iondyn%o%ave%sum_pressure_sq = 0.0_double
        do i = 1,3
           do j = 1,3
              iondyn%o%ave%sum_stress_tensor(i,j) = 0.0_double
              iondyn%o%ave%sum_stress_tensor_sq(i,j) = 0.0_double
           end do
        end do
        
        if (iondyn%o%md_method .ne. NONE) then

           if (debug) call warn('cons_tddft_ion_dynamics::initializing md_file')

           call my(file(trim(tddft_md_trajectory_path)),iondyn%o%md_file)
           if (i_access(iondyn%o%md_file)) open(x_unit(iondyn%o%md_file),file=x_name(iondyn%o%md_file))
           
           if (debug) call warn('cons_tddft_ion_dynamics::initializing new_vel_file')

           call my(file(trim(new_velocity_path)),iondyn%o%new_vel_file)

           if (debug) call warn('cons_tddft_ion_dynamics::calling x_cart_positions')

           call x_cart_positions(cfg,iondyn%o%current%pos_cart)

           if (debug) call warn('cons_tddft_ion_dynamics::if(generate_velocities)...')

           if (iondyn%o%generate_velocities == 0) then
              call my(file(trim(velocity_path)),iondyn%o%init_vel_file)
              if (i_access(iondyn%o%init_vel_file)) &
                   &inquire(file=x_name(iondyn%o%init_vel_file),exist=exist_file)
              if (i_comm(iondyn%o%init_vel_file)) call broadcast(MOD_SCOPE,exist_file)
              if (error(.not.exist_file,"ERROR: velocity file does not exist")) goto 100

              if (i_access(iondyn%o%init_vel_file)) &
                   &open(unit=x_unit(iondyn%o%init_vel_file), &
                   &file=x_name(iondyn%o%init_vel_file), iostat=ios)
              if (i_comm(iondyn%o%init_vel_file)) call broadcast(MOD_SCOPE,ios)
              if (error(ios /= 0,"ERROR: unable to open velocity file")) goto 100

              do i = 1,x_n_atoms(cfg)
                 if (i_access(iondyn%o%init_vel_file)) &
                      &read(x_unit(iondyn%o%init_vel_file),*,iostat=ios) &
                      &  iondyn%o%current%velocities(:,i)
                 if (i_comm(iondyn%o%init_vel_file)) call broadcast(MOD_SCOPE,ios)
                 if (error(ios /= 0,"ERROR: read error in velocity file")) goto 100
                 if (i_comm(iondyn%o%init_vel_file)) call broadcast(MOD_SCOPE,iondyn%o%current%velocities(:,i))
              end do
           else
              call init_velocities(iondyn%o%current%velocities,iondyn%o%current%masses,iondyn%o%init_temp)
           end if


           if (debug) call warn('cons_tddft_ion_dynamics::checking for restart file')
           
           if (present(restf)) then
              call my(restf)
              call read_restart_i(iondyn%o,restf)
              call glean(thy(restf))

              !RMH
              !call update_config(cfg,iondyn%o%current%pos_cart)
              !if (error()) goto 100
           else ! Only call update mdinfo if NOT restarting!!!
              if (debug) call warn('cons_tddft_ion_dynamics::no restart file, calling update_mdinfo')
              call update_mdinfo_i(cfg,iondyn%o)
           end if
   
        end if

        if (debug) call warn('cons_tddft_ion_dynamics::gleaning')

100     call glean(thy(cfg))

        if (error("Exit ion_dynamics_mod::cons_tddft_ion_dynamics")) continue

        if (debug) call warn('cons_tddft_ion_dynamics::exiting')

      end function cons_tddft_tags


      subroutine read_params_i(iondynr,cfg)
        type(ion_dynamics_rep) :: iondynr        
        type(config_td_obj)    :: cfg

        ! Local vars
        character(line_len) :: tag
        logical            :: fnd,exist_file
        real(double)       :: mass_temp
        integer            :: i
        integer            :: interval

        ! read the parameter defining the method
        !  =0: no molecular dynamics
        !  =1: NVE
        !  =2: NVT_RESCALE
        !  =3: NVT_ANDERSON        
        !  =4: NVT_HOOVER

        call my(cfg)
        
        !** Read in the type of method used to move the nuclei.
        call arg("tddft_md_method",tag,fnd)
        if (.not.fnd) tag = "none"
        call lowercase(tag)
        select case (tag(1:len_trim(tag)))
        case ("none")
           iondynr%md_method = NONE
           goto 100
        case ("nve","verlet")
           iondynr%md_method = NVE
        case ("nvt_rescale", "rescale")
           iondynr%md_method = NVT_RESCALE
        case ("nvt_anderson","anderson","stochastic")
           iondynr%md_method = NVT_ANDERSON
        case ("nvt_hoover","hoover")
           iondynr%md_method = NVT_HOOVER
        case default
           if (error(.true.,"cons_tddft_ion_dynamics: unrecognized md_method")) goto 100
        end select
              
        ! read md parameters 
        call arg("tddft_md_steps",iondynr%max_md_steps,fnd)
        if (.not.fnd) iondynr%max_md_steps = 0
        if (error(iondynr%max_md_steps<0,"cons_tddft_ion_dynamics: max_md_steps < 0")) goto 100

        call arg("tddft_md_skip_steps",iondynr%n_skip_steps,fnd)
        if (.not.fnd) iondynr%n_skip_steps=0
        if (error(iondynr%n_skip_steps<0,"cons_tddft_ion_dynamics: n_skip_steps < 0")) goto 100

        
!        call arg("tddft_md_dt",iondynr%time_step,fnd)
!        if (.not.fnd) iondynr%time_step = 100.
!        iondynr%time_step = iondynr%time_step*FS_2_ARU
           
!        if (error(iondynr%time_step<0,"cons_tddft_ion_dynamics: time_step < 0")) goto 100

        call arg("tddft_md_gen_velocities",tag,fnd)
        if (.not.fnd) tag = "yes"
        call lowercase(tag)
        select case (tag(1:len_trim(tag)))
        case ("yes","true","on")
           iondynr%generate_velocities = 1
        case ("no","false","off")
           iondynr%generate_velocities = 0
        case default
           if (error(.true.,"unrecognized response: tddft_md_gen_velocities")) goto 100
        end select
        
        call arg("tddft_md_init_temp",iondynr%init_temp,fnd)
        if (.not.fnd) iondynr%init_temp = 0.0_double
        if (error(iondynr%init_temp<0,"cons_tddft_ion_dynamics: init_temp < 0")) goto 100
        
        call arg("tddft_md_desired_temp",iondynr%desired_temp,fnd)
        if (.not.fnd) iondynr%desired_temp = iondynr%init_temp
        if (error(iondynr%desired_temp<0,"cons_tddft_ion_dynamics: desired_temp < 0")) goto 100
        
        call arg("tddft_md_temp_freq",iondynr%temp_mod_freq,fnd)
        if (.not.fnd) iondynr%temp_mod_freq=5
        if (error(iondynr%temp_mod_freq<1,"cons_tddft_ion_dynamics temp_mod_freq < 1")) goto 100
        
        call arg("tddft_md_hoover_mass",iondynr%current%hoover_mass,fnd)
        if (.not.fnd) iondynr%current%hoover_mass=1000.
        if (error(iondynr%current%hoover_mass <= 0.0_double,"cons_tddft_ion_dynamics hoover_mass < 0")) goto 100
 
        call arg("nucpos_file_type",tag,fnd)
        if (.not.fnd) tag = 'pdb'
        call lowercase(tag)
        select case (trim(tag))
        case ('pdb')
           iondynr%nucpos_file_type = PDB
        case ('xyz')
           iondynr%nucpos_file_type = XYZ
        case ('yes','true','.true.','on')
           iondynr%nucpos_file_type = PDB
        case ('no','false','.false.','off','none')
           iondynr%nucpos_file_type = NONE
        case default
           call warn('Unrecognized option for nucpos_file_type.')
           call warn('    not writing out separate MD files.')
           iondynr%nucpos_file_type = NONE
        end select
        
        call arg("nucpos_interval",interval,fnd)
        if (.not.fnd) interval = -1
        iondynr%nucpos_interval = interval
                
        call arg("mesh_file_type",tag,fnd)
        if (.not.fnd) tag = 'amira'
        call lowercase(tag)
        select case (trim(tag))
        case ('am','amira')
           iondynr%mesh_file_type = AMIRA
        case ('vtk')
           iondynr%mesh_file_type = VTK
        case default
           call warn('Unrecognized option for mesh_file_type.')
           call warn('        Not writing out seperate Ion_dynamics density files')
        end select

!        call arg("md_density_interval",interval,fnd)
!        if (.not.fnd) interval = -1
!        iondynr%density_interval = interval
        iondynr%density_interval = -1
        
        call arg("pressure",tag,fnd)
        if (.not.fnd) tag = "off"
        call lowercase(tag)
        select case (trim(tag))
        case ("on",".true.","true")
           iondynr%ave%ave_pressure = .true.
        case ("off",".false.","false")
           iondynr%ave%ave_pressure = .false.
        case default
           if (error(.true.,"ERROR: pressure tag was not recognized")) goto 100
        end select
        call arg("stress_tensor",tag,fnd)
        if (.not.fnd) tag = "OFF"
        call lowercase(tag)
        select case (trim(tag))
        case ("on",".true.","true")
           iondynr%ave%ave_stress_tensor = .true.
        case ("off",".false.","false")
           iondynr%ave%ave_stress_tensor = .false.
        case default
           if (error(.true.,"ERROR: pressure tag was not recognized")) goto 100
        end select
          
        iondynr%current%natoms = x_n_atoms(cfg)

        do i = 1,iondynr%current%natoms
           tag = "atom_mass_"//x_type(cfg,i)
           call arg(tag(1:len_trim(tag)),mass_temp,fnd)
           if (error(.not.fnd,"ERROR: mass not found")) goto 100
!           iondynr%current%masses(i) = AMU_2_ELECTRON_MASS*mass_temp
           iondynr%current%masses(i) = AMU_2_ARU*mass_temp
           if (error(iondynr%current%masses(i) <= 0.0_double,"ERROR: mass <= 0")) goto 100
           if ( i_access(diaryfile()) ) then
!              write(x_unit(diaryfile()),"('mass for atom',i6,' =',2g14.5)") &
!                   & i,iondynr%current%masses(i),iondynr%current%masses(i)/AMU_2_ELECTRON_MASS
              write(x_unit(diaryfile()),"('mass for atom',i6,' =',2g14.5)") &
                   & i,iondynr%current%masses(i),iondynr%current%masses(i)/AMU_2_ARU
           end if
        end do
        
100     call glean(thy(cfg))
        if (error("ERROR:ion_dynamics_mod:read_params_i exiting")) continue
        
      end subroutine read_params_i


      function ion_dynamics_move_atoms(iondyn,cfg) result (changed)
!doc$ function              move_atoms(cfg)
!  Crude routine to move the atoms using underlying MD propagation algorithms.
        type(ion_dynamics_obj) :: iondyn
        type(config_td_obj)    :: cfg
        logical                :: changed
!cod$        

        call my(cfg)
        call my(iondyn) 
        
        branch: SELECT CASE (iondyn%o%md_method)
        CASE (NONE)
           changed = .false.
        CASE (NVE)
           changed = v_verlet_tddft(cfg,iondyn); if (error()) goto 100
        CASE (NVT_RESCALE)
           changed = rescale_tddft(cfg,iondyn); if (error()) goto 100
        CASE (NVT_ANDERSON)
           changed = anderson_tddft(cfg,iondyn); if (error()) goto 100
        CASE (NVT_HOOVER)
           changed = hoover_tddft(cfg,iondyn); if (error()) goto 100
        CASE default
           if (error(.true.,"tddft_md_method not defined")) call diary(iondyn)
           changed = .false.
        END SELECT branch
           
100     call glean(thy(cfg))
        call glean(thy(iondyn))

        if (error("Exit ion_dynamics_mod::move_atoms")) continue

      end function 


      function v_verlet_tddft(cfg,iondyn) result(changed)
!doc$ function v_verlet_tddft(cfg,iondyn)
        type(config_td_obj)    :: cfg
        type(ion_dynamics_obj) :: iondyn
        logical                :: changed
! perform single move using NVE molecular dynamics with the 'velocity Verlet' algorithm

!cod$

! define a set of objects that will be used for the modification of the config
        ! Local Vars ---------------------------------------------
        integer :: natoms,i
        integer :: status
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double) :: time_step
        real(double) :: fdotv,fnrm
        logical :: found, debug
        !---------------------------------------------------------
        
        debug = .false.

        call my(cfg)
        call my(iondyn)
        
        if (debug) call warn('ion_dynamics::v_verlet_tddft - starting')

        natoms = x_n_atoms(cfg)

        changed = .false.

        time_step = iondyn%o%time_step

        allocate( frcs(3,natoms), accel(3,natoms), STAT=status )
        if (error(status /= 0,"v_verlet: allocate failure")) goto 100

        call x_cart_positions(cfg,iondyn%o%current%pos_cart)
 
        frcs = x_forces(cfg)
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:)
        end do

!RMH        call diary(cfg_static,iondyn%o%n_md_steps)

        iondyn%o%current%pos_cart = iondyn%o%current%pos_cart + &
             time_step*(iondyn%o%current%velocities + (time_step/2.0_double)*accel)

        iondyn%o%current%velocities = iondyn%o%current%velocities + (time_step/2.0_double)*accel

        if (debug) call warn('ion_dynamics::v_verlet_tddft - about to call update()')

!        call update_config(cfg,iondyn%o%current%pos_cart)
        if (error()) goto 100

        changed = .true.
        
        frcs = x_forces(cfg)
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:)
        end do

        iondyn%o%current%velocities = iondyn%o%current%velocities + (time_step/2.0_double)*accel

        iondyn%o%n_md_steps = iondyn%o%n_md_steps + 1

        if (debug) call warn('ion_dynamics::v_verlet_tddft - call update_mdinfo_i()')

        call update_mdinfo_i(cfg,iondyn%o)

        call diary(iondyn%o%current)

        call diary(iondyn%o%ave)


100     call glean(thy(iondyn))
        call glean(thy(cfg))

        if (error("Exit ion_dynamics_mod:v_verlet_tddft")) continue

        if (debug) call warn('ion_dynamics::v_verlet_tddft - exiting...')


      end function


      function rescale_tddft(cfg,iondyn) result(changed)
!doc$ function rescale_tddft(cfg,iondyn)
        type(config_td_obj) :: cfg
        type(ion_dynamics_obj)       :: iondyn
        logical :: changed
! perform a single move using NVE molecular dynamics with the 'velocity Verlet' algorithm

!cod$

! define a set of objects that will be used for the modification of the config
        integer :: natoms,i
        integer :: status
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double) :: time_step
        real(double) :: fdotv,fnrm
        logical :: found
        
        call my(cfg)
        call my(iondyn)
        
        changed = .false.

        time_step = iondyn%o%time_step
        
        allocate( frcs(3,iondyn%o%current%natoms), &
             & accel(3,iondyn%o%current%natoms), STAT=status )
        
        if (error(status /= 0,"relax_md: allocate failure")) goto 100

        call x_cart_positions(cfg,iondyn%o%current%pos_cart)

        frcs = x_forces(cfg)
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:)
        end do
                
!RMH fix        call diary(cfg_static,iondyn%o%n_md_steps)
        
        iondyn%o%current%pos_cart = iondyn%o%current%pos_cart + &
             time_step*(iondyn%o%current%velocities + (time_step/2.0_double)*accel)
        
        iondyn%o%current%velocities = iondyn%o%current%velocities + (time_step/2.0_double)*accel
        
!RMH
!        call update_config(cfg,iondyn%o%current%pos_cart)
        if (error()) goto 100
        
        changed = .true.
        
        frcs = x_forces(cfg)
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:)
        end do
        
        iondyn%o%current%velocities = iondyn%o%current%velocities + (time_step/2.0_double)*accel
        
        if (mod(iondyn%o%n_md_steps,iondyn%o%temp_mod_freq) == 0) then
           call rescale_velocity(iondyn,iondyn%o%desired_temp)
        end if
        
        iondyn%o%n_md_steps = iondyn%o%n_md_steps + 1
        
        call update_mdinfo_i(cfg,iondyn%o)
        
        call diary(iondyn%o%current)
        
        call diary(iondyn%o%ave)
        

100     call glean(thy(iondyn))
        call glean(thy(cfg))

        if (error("Exit ion_dynamics_mod:rescale_tddft")) continue

      end function 


      function anderson_tddft(cfg,iondyn) result(changed)
!doc$ function anderson_tddft(cfg,iondyn)
        type(config_td_obj) :: cfg
        type(ion_dynamics_obj)       :: iondyn
        logical :: changed
! perform a single move using NVT molecular dynamics with the 'velocity Verlet' algorithm coupled with 
! the anderson thermostat (coupling to stochastic bath)

!cod$

        integer :: natoms,i
        integer :: status
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double), allocatable, dimension(:) :: test_numbers
        real(double) :: time_step
        real(double) :: fdotv,fnrm
        logical :: found
        
        call my(cfg)
        call my(iondyn)
        
        changed = .false.

        time_step = iondyn%o%time_step
        
        allocate( frcs(3,iondyn%o%current%natoms), STAT=status)
        if (error(status /= 0,"relax_md: allocate failure - forces")) goto 100
        allocate( accel(3,iondyn%o%current%natoms), STAT=status)
        if (error(status /= 0,"relax_md: allocate failure - accel")) goto 100
        allocate( test_numbers(iondyn%o%current%natoms), STAT=status )
        if (error(status /= 0,"relax_md: allocate failure - test")) goto 100

        call x_cart_positions(cfg,iondyn%o%current%pos_cart)

        frcs = x_forces(cfg)
        
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:)
        end do

!RMH         call diary(cfg_static,iondyn%o%n_md_steps)

        iondyn%o%current%pos_cart = iondyn%o%current%pos_cart + &
             time_step*(iondyn%o%current%velocities + (time_step/2.0_double)*accel)
        
        iondyn%o%current%velocities = iondyn%o%current%velocities + (time_step/2.0_double)*accel
        
        !RMH
        !call update_config(cfg,iondyn%o%current%pos_cart)
        !if (error()) goto 100
        
        changed = .true.
        
        frcs = x_forces(cfg)
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:)
        end do
        
        iondyn%o%current%velocities = iondyn%o%current%velocities + (time_step/2.0_double)*accel
        
        if (mpi_isroot(MOD_SCOPE)) call random_number(test_numbers)
        call broadcast(MOD_SCOPE,test_numbers)
        do i = 1,iondyn%o%current%natoms
           if (test_numbers(i) < 1.0_double/iondyn%o%temp_mod_freq) then
              iondyn%o%current%velocities(:,i) = &
                   & sqrt(iondyn%o%desired_temp*K_BOLTZMAN/iondyn%o%current%masses(i))&
                   & *rgauss(3)
           end if
        end do
        
        iondyn%o%n_md_steps = iondyn%o%n_md_steps + 1
        
        call update_mdinfo_i(cfg,iondyn%o)
        
        call diary(iondyn%o%current)
        
        call diary(iondyn%o%ave)
           
100     call glean(thy(iondyn))
        call glean(thy(cfg))

        if (error("Exit ion_dynamics_mod:anderson_tddft")) continue

      end function 


      function hoover_tddft(cfg,iondyn) result(changed)
!doc$ function hoover_tddft(cfg,iondyn)
        type(config_td_obj)    :: cfg
        type(ion_dynamics_obj) :: iondyn
        logical :: changed
! perform a single move using NVT molecular dynamics with the Hoover implementation of Nose Thermostat 
! algorithm.  Follows algorithm in "Understanding Molecular Simulation", by Daan Frenkel and Berend Smit

!cod$

        real(double), parameter :: SOLVE_TOL = 1.d-6
        integer, parameter :: MAX_SOLVE_TRIES = 12

        integer :: natoms,i
        integer :: status
        integer :: n_solve_tries
        real(double), allocatable, dimension(:,:) :: frcs, accel
        real(double), allocatable, dimension(:,:) :: b,c,h
        real(double), allocatable, dimension(:,:) :: vel_prime
        real(double), allocatable, dimension(:,:) :: v1,v2
        real(double) :: a,d,hplus1,del_xi,xi1,xi2
        real(double) :: diff_ke
        real(double) :: test_val
        real(double) :: xi_prime
        real(double) :: time_step
        real(double) :: fdotv,fnrm
        logical :: found,solved
        
        call my(cfg)
        call my(iondyn)
        
        natoms = x_n_atoms(cfg)

        changed = .false.

        time_step = iondyn%o%time_step
        
        allocate( frcs(3,natoms), accel(3,natoms), STAT=status )
        if (error(status /= 0,"hoover_md: allocate failure")) goto 100
        allocate( b(3,natoms), c(3,natoms), h(3,natoms), STAT=status )
        if (error(status /= 0,"hoover_md: allocate failure")) goto 100
        allocate( v1(3,natoms), v2(3,natoms), vel_prime(3,natoms), STAT=status )
        if (error(status /= 0,"hoover_md: allocate failure")) goto 100
        
        call x_cart_positions(cfg,iondyn%o%current%pos_cart)

        frcs = x_forces(cfg)
        
        ! start the molecular dynamics loop

!RMH        call diary(cfg_static,iondyn%o%n_md_steps)
        
        do i = 1,3
           accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:) - &
                & iondyn%o%current%hoover_xi*iondyn%o%current%velocities(i,:)
        end do
        
        iondyn%o%current%pos_cart = iondyn%o%current%pos_cart + &
             time_step*(iondyn%o%current%velocities + (time_step/2.0_double)*accel)
        
        vel_prime = iondyn%o%current%velocities + (time_step/2.0_double)*accel
        
        diff_ke = 0.0_double
        do i = 1,natoms
           diff_ke = diff_ke + &
                & iondyn%o%current%masses(i)*&
                & dot_product(iondyn%o%current%velocities(:,i),iondyn%o%current%velocities(:,i))
        end do
        diff_ke = diff_ke - 3.0_double*natoms*K_BOLTZMAN*iondyn%o%desired_temp
        
        iondyn%o%current%hoover_s = iondyn%o%current%hoover_s + &
             & iondyn%o%current%hoover_xi*time_step + &
             & diff_ke*time_step**2/(2.0_double*iondyn%o%current%hoover_mass)
        
        xi_prime = iondyn%o%current%hoover_xi + &
             & diff_ke*time_step/(2.0_double*iondyn%o%current%hoover_mass)
        
        !RMH
        !call update_config(cfg,iondyn%o%current%pos_cart)
        !if (error()) goto 100
        
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
              accel(i,:) = frcs(i,:)/iondyn%o%current%masses(:) - &
                   & xi2*v2(i,:)
              b(i,:) = iondyn%o%current%masses(:)*v2(i,:)*&
                   &time_step/iondyn%o%current%hoover_mass
           end do
           
           diff_ke = 0.0_double
           do i = 1,natoms
              diff_ke = diff_ke + &
                   & iondyn%o%current%masses(i)*dot_product(v2(:,i),v2(:,i))
           end do
           diff_ke = diff_ke - 3.0_double*natoms*K_BOLTZMAN*iondyn%o%desired_temp
           
           h = vel_prime + accel*(time_step/2.0_double) - v2
           hplus1 = xi_prime - xi2 + &
                & diff_ke*time_step/(2.0_double*iondyn%o%current%hoover_mass)
           
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
           if (error(n_solve_tries > MAX_SOLVE_TRIES,"hoover_md: solve failed")) goto 100
           
        end do
        
        iondyn%o%current%velocities = v1
        iondyn%o%current%hoover_xi = xi1
        
        iondyn%o%n_md_steps = iondyn%o%n_md_steps + 1
        
        call update_mdinfo_i(cfg,iondyn%o)
        
        call diary(iondyn)
        
100     deallocate(frcs,accel,b,c,h,v1,v2,vel_prime)
        
        call glean(thy(iondyn))
        call glean(thy(cfg))

        if (error("Exit ion_dynamics_mod: hoover_tddft")) continue

      end function


      !** Write the positions of the atoms out to a file in PDB file format
      subroutine write_pdb_file_i(cfg,iondynr)
        type(config_td_obj)    :: cfg
        type(ion_dynamics_rep) :: iondynr
        ! Local Vars ------------------------------------------------------------
        type(file_obj)     :: f
        character(line_len) :: filename
        character(line_len) :: atom_type
        integer            :: ia, na
        integer            :: ios
        logical            :: debug
        !========================================================================
        
        debug = .false.

        if (debug) call warn('ion_dynamics::write_pdb_file_i - starting')

        call my(cfg)

        !** Construct a file name
        call get_nucpos_filename_i(iondynr%n_md_steps,filename,iondynr%nucpos_file_type)

        if (debug) call warn('ion_dynamics::write_pdb_file_i - initializing the file object')

        !** Initialize the file object
        call my(file(trim(filename)),f)

        if (debug) call warn('ion_dynamics::write_pdb_file_i - open the file')
        
        !** Open the file
        if (i_access(f)) then
           open(unit=x_unit(f), &
                   file=x_name(f), &
                   iostat=ios)
        end if
        if (i_comm(f)) call broadcast(MOD_SCOPE,ios)
        if (error(ios /= 0,"ERROR: error opening file")) goto 100

        if (debug) call warn('ion_dynamics::write_pdb_file_i - writing the file')

        !** Write the postions of the atoms to a file in PDB file format
        if (i_access(f)) then
           na = size(iondynr%current%pos_cart,2)
           do ia=1,na
              atom_type = x_type(cfg,ia)
              write(x_unit(f),'(a6,i5,a3,a6,f18.3,2f8.3)') &
                   "HETATM", ia, atom_type(1:2), "MOL", &
                   iondynr%current%pos_cart(:,ia)
           end do

           write(x_unit(f),'(a3)') "TER"

           !** Close the file
           close(unit=x_unit(f))

        end if

100     call glean(thy(cfg))
        call glean(thy(f))
        if (error("ion_dynamics::write_pdb_file_i")) continue

        if (debug) call warn('ion_dynamics::write_pdb_file_i - exiting')

      end subroutine write_pdb_file_i


      !** Write the positions of the atoms out to a file in XYZ file format
      subroutine write_xyz_file_i(cfg,iondynr)
        type(config_td_obj) :: cfg
        type(ion_dynamics_rep) :: iondynr
        ! Local Vars ------------------------------------------------------------
        type(file_obj)     :: f
        character(line_len) :: filename
        character(line_len) :: atom_type
        integer            :: ia, na
        integer            :: ios
        !========================================================================
        
        call my(cfg)

        !** Construct a file name
        call get_nucpos_filename_i(iondynr%n_md_steps,filename,iondynr%nucpos_file_type)
        if (error()) goto 100

        !** Initialize the file object
        call my(file(trim(filename)),f) ; if (error()) goto 100
        
        !** Open the file
        if (i_access(f)) then
           open(unit=x_unit(f), &
                   file=x_name(f), &
                   iostat=ios)
        end if
        if (i_comm(f)) call broadcast(MOD_SCOPE,ios)
        if (error(ios /= 0,"ERROR: error opening file")) goto 100

        !** Write the atom postions to a file in xyz format
        if (i_access(f)) then
           na = size(iondynr%current%pos_cart,2)
           write(x_unit(f),'(i6)') na
           write(x_unit(f),'(a14)') "Some Material"
           do ia=1,na
              atom_type = x_type(cfg,ia)
              write(x_unit(f),'(a3,3f10.3)') atom_type(1:2), iondynr%current%pos_cart(:,ia)
           end do

           !** Close the file
           close(unit=x_unit(f))

        end if

100     call glean(thy(cfg))
        call glean(thy(f))
        if (error("ion_dynamics::write_xyz_file_i")) continue

      end subroutine


      !** Return the name of a file that is to contain the nuclear positions
      subroutine get_nucpos_filename_i(istep,filename,filetype)
        integer            :: istep
        character(line_len) :: filename
        integer            :: filetype

        ! Local Vars ------------------------------------------------------------
        character(line_len) :: num
        integer            :: j, pos, a
        logical            :: debug
        !========================================================================

        debug = .false.

        if (debug) write(*,*) 'ion_dynamics::get_nucpos_filename_i - starting'
        write(num,"(I8)") istep

        if (filetype == PDB) then
           filename = "nucpos_        .pdb"
        else if (filetype == XYZ) then
           filename = "nucpos_        .xyz"
        else
           if (error(.true.,"Error - unrecognized filetype")) goto 100
        end if

        pos = 8

        do j=1,8
           a = iachar(num(j:j) )
           if (a == 32) then
              filename(pos:pos) = "0"
              pos = pos + 1
           else
              filename(pos:pos) = num(j:j)
              pos = pos + 1
           end if
        end do

100     if (error("ion_dynamics::get_nucpos_filename_i - Exiting")) continue

        if (debug) write(*,*) 'ion_dynamics::get_nucpos_filename_i - exiting'

      end subroutine

        
      subroutine get_density_vtk_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i
        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "density_        .vtk"
        pos = 9
        !write(*,*) "in get_density_filename_i"
        do j=1,8
           a = iachar(num(j:j) )
           !write(*,*) "  char(j)", j, " ", a
           if (a == 32) then
              filename(pos:pos) = "0"
              pos = pos + 1
           else
              filename(pos:pos) = num(j:j)
              pos = pos + 1
           end if
        end do

        !write(*,*) "num =", num
        !write(*,*) "filename = ", filename

      end subroutine


      subroutine get_density_am_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i
        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "density_        .am"
        pos =9 
        !write(*,*) "in get_density_filename_i"
        do j=1,8
           a = iachar(num(j:j) )
           !write(*,*) "  char(j)", j, " ", a
           if (a == 32) then
              filename(pos:pos) = "0"
              pos = pos + 1
           else
              filename(pos:pos) = num(j:j)
              pos = pos + 1
           end if
        end do

        !write(*,*) "num =", num
        !write(*,*) "filename = ", filename

      end subroutine




      !*****************************************************************************
      !
      ! Utilities
      !
      !*****************************************************************************
      
      subroutine uppercase(str)
        character(*), intent(INOUT) :: str
        integer  :: i, j, k
        
        j = len(Str)
        do i=1, j
           k = iachar(str(i:i))
           if ((k>96) .and. (k<123)) str(i:i) = achar(k-32)
        end do
        
        return
      end subroutine 
      
      subroutine lowercase(str)
        character(*), intent(INOUT) :: str
        
        integer  :: i, j, k
        
        j = len(Str)
        
        do i=1, j
           k = iachar(str(i:i))
           if ((k>64) .and. (k<91)) str(i:i) = achar(k+32)
        end do
        
        return
      end subroutine
      



subroutine write_info_ion_dynamics(iondyn)
  type(ion_dynamics_obj) :: iondyn

  call my(iondyn)

  write(*,*) 'ion_dynamics::write_info starting'
  write(*,*) 'ghost = ', iondyn%o%g
  write(*,*) 'md_method = ', iondyn%o%md_method
  write(*,*) 'n_md_steps = ', iondyn%o%n_md_steps
  write(*,*) 'n_skip_steps = ', iondyn%o%n_skip_steps
  write(*,*) 'max_md_steps = ', iondyn%o%max_md_steps
  write(*,*) 'temp_mod_freq = ', iondyn%o%temp_mod_freq
  write(*,*) 'generate_velocities = ', iondyn%o%generate_velocities
  write(*,*) 'time_step = ', iondyn%o%time_step
  write(*,*) 'init_temp = ', iondyn%o%init_temp
  write(*,*) 'desired_temp = ', iondyn%o%desired_temp
  write(*,*) 'init_temp = ', iondyn%o%init_temp
  write(*,*) 'n_ave_steps = ', iondyn%o%ave%n_ave_steps
  write(*,*) 'sum_temperature = ', iondyn%o%ave%sum_temperature
  write(*,*) 'sum_temperature_sq = ', iondyn%o%ave%sum_temperature_sq
  write(*,*) 'sum_kin_energy = ', iondyn%o%ave%sum_kin_energy
  write(*,*) 'sum_kin_e_sq = ', iondyn%o%ave%sum_kin_e_sq
  write(*,*) 'sum_pot_energy = ', iondyn%o%ave%sum_pot_energy
  write(*,*) 'sum_pot_e_sq = ', iondyn%o%ave%sum_pot_e_sq
  write(*,*) 'sum_tot_energy = ', iondyn%o%ave%sum_tot_energy
  write(*,*) 'sum_tot_e_sq = ', iondyn%o%ave%sum_tot_e_sq
  write(*,*) 'sum_pressure = ', iondyn%o%ave%sum_pressure
  write(*,*) 'sum_pressure_sq = ', iondyn%o%ave%sum_pressure_sq
  write(*,*) 'natoms = ', iondyn%o%current%natoms
  write(*,*) 'time = ', iondyn%o%current%time
  write(*,*) 'kin_energy = ', iondyn%o%current%kin_energy
  write(*,*) 'pot_energy = ', iondyn%o%current%pot_energy
  write(*,*) 'tot_energy = ', iondyn%o%current%tot_energy
  write(*,*) 'pos_cart = ', iondyn%o%current%pos_cart
  write(*,*) 'velocities = ', iondyn%o%current%velocities
  write(*,*) 'masses = ', iondyn%o%current%masses

  call glean(thy(iondyn))

end subroutine





    end module

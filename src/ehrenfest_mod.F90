!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module ehrenfest_mod
!doc$ module ehrenfest_mod

!     There are no publicly available data types here.
!
!     ehrenfest_mod is responsible for running a time dependent simulation where the electrons
!     and the atoms are allowed to propagate in time (the atoms do not have to move in time).
!     Tags are used to read in the following parameters which control the simulation:
!
!     tddft_step_size dt  
!        dt should be entered in units of fs (10E-15 seconds)
!
!
!     ....
!
!

      use arg_mod
      use atomic_operators_mod
      use config_td_mod
      use electrons_td_mod
      use error_mod
      use atoms_mod
      use lattice_mod
      use crystal_mod
      use external_mod
      use fields_td_mod
      use gen_density_mod
      use gen_potential_mod
      use grid_mod
      use kind_mod
      use kpoints_mod
      use interrupt_mod
      use io_mod
      use layout_mod
      use ion_dynamics_mod
      use mpi_mod
      use multibasis_mod
      use multivector_mod
      use shortcuts_mod
      use tagio_mod
      use wavefunctions_td_mod
      use config_sc_mod
      use electrons_sc_mod
      use wavefunctions_es_mod

!      use strings
!      use config_sc_mod
!      use config_fh_mod
!      use grid_math_mod

!use hamiltonian_mod
!use crystal_mod
!use atoms_mod
!use lattice_mod
!cod$
      implicit none
      private

      !** State variables (not read in from argvf).
      type ehrenfest_state
         real(double)  :: time                 !** Simulation time for electrons
         integer       :: istep                !** current electronic time step 
         integer       :: atoms_istep          !** current md time step
         real(double)  :: last_partial_charge  !** 
         real(double)  :: last_partial_time    !** 
         integer       :: num_sc_steps         !** number of iterations from last electronic time step
         type(ion_dynamics_obj) :: iondyn      !** iondyn holds an object that can be used to propagate the ions
         type(config_td_obj)    :: cfg_td      !** config_td object holds the progated config
         type(external_obj)     :: ext         !** ext holds the external potentials including atom 
                                               !   positions and free fields
      end type  

      !** The following parameters are read in from argvf to control the time dependent simulation
      type ehrenfest_params
         real(double)   :: runtime      ! Total time of the calculation.
         real(double)   :: dt           ! Size of each electronic time step
         real(double)   :: ions_dt     ! Size of atomic time step
         real(double)   :: ke_cutoff    ! Kinetic energy cutoff (used for numerical stability)
         real(double)   :: den_cutoff   ! Density cutoff (used for numerical stability)
         real(double)   :: den_beta     ! smoothing parameter for cutoffs.
         real(double)   :: surface_origin(3)
         real(double)   :: surface_direction(3)
         character(line_len) :: restart_filename     ! Name of the restart file
         character(line_len) :: occupations_filename ! Name of the file containing the new occupations.
!         logical        :: move_atoms
         logical        :: selfconsistent
         integer        :: mesh_file_type       ! Mesh format for any write to file
         integer        :: move_atoms_interval  ! number of electronic time steps per MD time step 
                                                !   (NB this is not directly read in)
         integer        :: restf_interval       ! num el time steps per write to restart file
         integer        :: screen_dump_interval ! num el time steps per write to screen (STDOUT)
         integer        :: file_dump_interval   ! num el time steps per write of tddft quantities to *out files
         integer        :: den_interval         ! num el time steps per write of charge density to file
         integer        :: currden_interval     ! num el time steps per write of current density to file
         integer        :: pot_interval         ! num el time steps per write of total local potential to file
         integer        :: hap_interval         ! num el time steps per write of hartree potential to file
         integer        :: wf_interval          ! num el time steps per write of wavefunctions to file
         integer        :: proj_interval        ! num el time steps per write of projections to file
         integer        :: orth_interval        ! num el time steps per write of orthogonalization to file
         integer        :: em_energy_interval   ! num el time steps per write of EM energy terms to file
         integer        :: joule_power_interval ! num el time steps per write of Joule power terms to file
         integer        :: eigval_interval      ! num el time steps per write of eigenvalues to file
         integer        :: partial_charge_interval  ! num el time steps per write of partial charge to file
      end type

      !** Files that track time dependent quantities of interest
      type ehrenfest_files
         type(file_obj)  :: tddft_out
         type(file_obj)  :: fields_out
         type(file_obj)  :: energy_out
         type(file_obj)  :: proj_out
         type(file_obj)  :: proj_curr_out
         type(file_obj)  :: partial_charge_out
         type(file_obj)  :: eigval_out
         type(file_obj)  :: orth_out
      end type 

      !** Global parameters for this module
      integer, parameter         :: NONE = 0
      integer, parameter         :: MOD_SCOPE = CONFIG
      ! Grid file type
!      integer, parameter :: VTK = 1            ! VTK format
!      integer, parameter :: AMIRA_REGULAR = 2  ! Regular mesh in amira format
      
      !** Conversion Constants
      real(double), parameter    :: FS_2_ARU = 1.0_double/0.048377687
      real(double), parameter    :: ARU_2_FS = 0.048377687  ! convert from ARU time units to femtoseconds
      real(double), parameter    :: ARU_2_VoltsPerAngstrom = 36.360903 ! convert from ARU to Volts per Angstrom
      real(double), parameter    :: PHZ_2_ARU = 0.048377687
      real(double), parameter    :: EV_2_ARU = 1.0_double/13.605691
      real(double), parameter    :: VoltsPerAngstrom_2_ARU = 1.0_double/36.360903      
      real(double), parameter    :: Volts_2_ARU = 1.0_double/19.241363

      !** Default parameters for this module
      real(double), parameter :: default_dt = 1.0d-3    ! In units of fs


!doc$
      public :: ehrenfest_dynamics
      
!cod$

    contains
      
!******************************************************************************************
      
      subroutine ehrenfest_dynamics()
!doc$ subroutine ehrenfest_dynamics()
!effects : manages the ehrenfest run at the highest level
!cod$   

        ! Local Vars
        type(ehrenfest_params) :: params
        type(ehrenfest_state)  :: state
        type(ehrenfest_files)  :: files
        real(double),parameter :: tol = 1.0d-5
        integer                :: ik,nk
        logical                :: debug
        logical                :: fnd

        !-----------------------------------------------------------------------------

        debug = .false.

        if (debug) call warn("ehrenfest_dynamics starting")

        
        !** Initialize the ehrenfest params 
        if (debug) call warn("ehrenfest_dynamics: calling init_params_i")
        call init_params_i( params )      ; if (error()) goto 100

        !** Initialize the ehrenfest state variables
        if (debug) call warn("ehrenfest_dynamics: calling init_state_i")
        call init_state_i( state, params) ; if (error()) goto 100

        !** Initialize the ehrenfest output files
        if (debug) call warn("ehrenfest_dynamics: calling init_files_i")
        call init_files_i( files, params ); if (error()) goto 100

        !** print out initial values of various quantities to output files
        if (debug) call warn("ehrenfest_dynamics: dump output ")
        call dump_output_i(state,params,files,disable_restart=.true.)  ; if (error()) goto 100

        
        !*********************************************************************
        !** Main TDDFT loop. *************************************************
        if (debug) call warn("ehrenfest_dynamics: before do while(time < runtime)")
        do while ( state%time < params%runtime )

           !** Update the external potential (move atoms and update any free fields)
           if (debug) call warn("ehrenfest_dynamics:main_loop: before update_ext_i ")
           call update_ext_i(state, params)

           !** This update will result in a single electronic time step.
           if (debug) call warn("ehrenfest_dynamics:main_loop: before propagate ")
!           call update(state%cfg_td,state%ext,state%time,params%dt,params%selfconsistent,state%num_sc_steps)
           call update(state%cfg_td, state%ext)
           if (error()) goto 100
           
           !** update the time
           state%time = state%time + params%dt
           state%istep = state%istep + 1

           !** update the number of self-consistent steps
           state%num_sc_steps = x_num_sc_steps(state%cfg_td)
           
           !** Call a routine to dump various info out to file and or console
           call dump_output_i(state,params,files)  ; if (error()) goto 100

           if (debug) call warn("ehrenfest_dynamics::main loop - calling user_stop()")
              
           if (user_stop()) then
              call warn("WARNING: USER INITIATED STOP")
              exit
           end if

        end do
        !---------------------------------------------------------------------
        !---------------------------------------------------------------------

        !** Clean up the mess
        call close_files_i(files,params)
        call close_state_i(state,params)

100     if (error("Exit ehrenfest_mod::ehrenfest_dynamics")) continue

        if (debug) call warn ("ehrenfest_mod::ehrenfest_dynamics() closed state, exiting")

      end subroutine


      !** Private routines

      !******************************************************************************
      !**  init_params_i
      !** 
      !** This routine is responsible for reading input parameters from argvf
      subroutine init_params_i( params )

        type(ehrenfest_params) :: params

        !local vars
        character(line_len) :: tag
        logical            :: fnd
        integer            :: interval
        real(double)       :: cutoff, beta
        real(double)       :: orig(3)
        real(double)       :: dir(3)

        
        !** Get the total time of the calculation
        call arg("tddft_runtime",params%runtime,fnd)
        if (.not.fnd) params%runtime=1.0_double  ! Default is 1fs
        ! Convert from femtoseconds to ARU time units
        params%runtime = params%runtime*FS_2_ARU
        if (error(params%runtime<0,"ERROR: run_time < 0")) goto 100
        
        ! Size of each time step - N.B. the step size is read in assuming units of fs and
        !    then converted to atomic rydberg time units (4.84E-17s)
        call arg("tddft_step_size",params%dt,fnd)
        if (.not.fnd) params%dt = default_dt
        if (error(params%dt<0,"ERROR! step_size < 0")) goto 100
        params%dt = params%dt*FS_2_ARU

        !** Set the interval over which the atoms will be moved. Must be an integral multiple
        !     of the electronic time step.
        call arg("tddft_md_dt",params%ions_dt,fnd)
        if (.not.fnd) params%ions_dt = default_dt*FS_2_ARU
        if (error(params%ions_dt<0,"ERROR! step_size < 0")) goto 100
        params%ions_dt = params%ions_dt*FS_2_ARU

        if (params%ions_dt < params%dt) then
           params%move_atoms_interval = 1
        else
           params%move_atoms_interval = floor(params%ions_dt/params%dt)
        end if
        ! rescale so the time steps sync correctly
        params%ions_dt = params%ions_dt*(params%move_atoms_interval*1.0_double)

        !** ke_cutoff is the energy cutoff for the kinetic energy operator
        call arg("tddft_ke_cutoff", cutoff,fnd)
        if (.not.fnd) cutoff = -1.0_double
        params%ke_cutoff = cutoff

        !** den_cutoff is the energy cutoff for the density
        call arg("tddft_den_cutoff", cutoff,fnd)
        if (.not.fnd) cutoff = -1.0_double
        params%den_cutoff = cutoff

        !** den_beta is the smoothing factor that is used to cutoff the density
        call arg("tddft_den_beta",beta,fnd)
        if (.not.fnd) beta = 1.0_double
        params%den_beta = beta

        !** surface origin is a point that along with the surface direction defines a
        !    directed surface that separates the supercell into two regions.
        call arg("surface_origin",orig,fnd)
        if (.not.fnd) orig = (/0.0,0.0,0.0/)
        params%surface_origin = orig

        !** surface direction is a point that along with the surface origin defines a
        !    directed surface that separates the supercell into two regions.
        call arg("surface_direction",dir,fnd)
        if (.not.fnd) dir = (/0.0,0.0,0.0/)
        params%surface_direction = dir

        !** Check for a restart file.
        call arg("tddft_restart",tag,fnd)
        if (.not.fnd) tag = "off"
        params%restart_filename = trim(tag)

        !** Determine if new occupations are to be read in and then read them in 
        !    if indicated
        call arg("tddft_occupations",tag,fnd)
        if (.not.fnd) tag = "off"
        params%occupations_filename = trim(tag)

        !** Read in whether the time steps should be made self-consistently
        call arglc("tddft_self_consistent",tag,fnd)
        if (.not.fnd) tag = "false"
        select case (tag(1:len_trim(tag)))
        case ("yes",".true.","true","y","on")
           params%selfconsistent = .true.
        case default
           params%selfconsistent = .false.
        end select

        !** The mesh_file_type is the file format to which any 3d mesh quantities will be written
        call arglc("mesh_file_type",tag,fnd)
        if (.not.fnd) tag = "amira"
        select case(tag(1:len_trim(tag)))
        case("vtk")
           params%mesh_file_type = VTK
        case default
           params%mesh_file_type = AMIRA
        end select

        !** Read in whether the atoms should move
        call arglc("tddft_move_atoms",tag,fnd)
        if (.not.fnd) tag = "false"
        select case (tag(1:len_trim(tag)))
        case ("yes",".true.","true","y","on")
           params%selfconsistent = .true.
        case default
           params%selfconsistent = .false.
        end select

        !** The params%restf_interval is the number of time steps per write of a restart file
        call arg("tddft_write_restart_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%restf_interval = interval

        !** The screen dump interval is the number of time steps per write of various info to console
        call arg("tddft_screen_dump_interval",interval,fnd)
        if (.not.fnd) interval = 100
        params%screen_dump_interval = interval

        !** The tddft file dump interval is the number of time steps per write 
        !     current, efield, vecpot etc... to file
        call arg("tddft_file_dump_interval",interval,fnd)
        if (.not.fnd) interval = 10
        params%file_dump_interval = interval

        !** The den_interval is the number of time steps per write of a density file
        call arg("tddft_write_den_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%den_interval = interval

        !** The currden_interval is the number of time steps per write of current density file
        call arg("tddft_write_currden_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%currden_interval = interval

        !** The pot_interval is the number of time steps per write of a total potential file
        call arg("tddft_write_pot_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%pot_interval = interval

        !** The hap_interval is the number of time steps per write of the hartree potential file
        call arg("tddft_write_hap_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%hap_interval = interval

        !** The wf_interval is the number of time steps per write of the
        !     wavefunctions to files
        call arg("tddft_write_wf_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%wf_interval = interval

        !** The projection_interval is the number of time steps per calculation and
        !     write out of the projection of the time evolved wavefunctions on the 
        !     original ground state wavefunctions.
        call arg("tddft_projection_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%proj_interval = interval

        !** The orth interval determines how often the internal overlap of the 
        !     wavefunctions are calculated and written out to file
        call arg("tddft_orth_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%orth_interval = interval

        !** The em_energy_interval is the number of time steps per write of the
        !     electromagnetic energy to file
        call arg("tddft_em_energy_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%em_energy_interval = interval

        !** The joule_power_interval is the number of time steps per write of the
        !     joule power to file
        call arg("tddft_joule_power_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%joule_power_interval = interval

        !** The eigenvalue interval determines how often the eigenvalues
        !     are calculated and written out to file.
        call arg("tddft_eigval_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%eigval_interval = interval

        !** The partial_charge interval determines how often the amount of charge
        !    on the positive side of the surface defined by surface_origin and
        !    surface_direction is calculated.
        call arg("tddft_partial_charge_interval",interval,fnd)
        if (.not.fnd) interval = -1
        params%partial_charge_interval = interval

100     if (error("Exit ehrenfest_mod::init_params_i")) continue

      end subroutine init_params_i

      
      !********************************************************************
      !** 
      !**  init_state_i
      !**
      !** This routine is responsible for initializing the state including
      !**   restarting from a file if indicated.
      subroutine init_state_i( state, params )
        type(ehrenfest_state)  :: state
        type(ehrenfest_params) :: params

        !local vars --------------
        integer(long)        :: dsize, s4, ndata, ios
        type(tagio_obj)      :: restf
        logical              :: debug
        integer              :: ik,nk,nb
        real(double),pointer :: new_occs(:)
        real(double)         :: cutoff, beta
        character(line_len)   :: usage
        type(multibasis_obj) :: mb
        character(1)         :: tios
        !-----------------------------------------------------

        debug = .false.
        
        if (debug) call warn("init_state_i:: starting ")


        !*********************************************************************************************
        !** Check to see if this is a restart
        select case(params%restart_filename(1:len_trim(params%restart_filename)))
        !---------------------------------------------------------------------------------------------
        !-------------- Start at time = 0 ------------------------------------------------------------
        case ("none","off","false") !** Starting (i.e. not restarting)
        
           state%time = 0.0_double                !** Simulation time
           state%istep = 0                        !** current electronic time step 
           state%atoms_istep = 0                  !** current md time step
           state%last_partial_charge = 0.0_double !** total charge in some subsection of supercell from last step
           state%last_partial_time = 0.0_double   !** time from last sampling of the partial charge
           state%num_sc_steps = 1                 !** number of iterations from last electronic time step

           call my(config_td(),state%cfg_td)                      ; if (error()) goto 100
           call diary(state%cfg_td)
           call my(x_external(state%cfg_td),state%ext)            ; if (error()) goto 100

           call my(ion_dynamics(state%cfg_td, params%ions_dt),state%iondyn) ; if (error()) goto 100

           nk = x_n_kpoints(x_kpoints(x_electrons(state%cfg_td)))  ; if (error()) goto 100
           nb = x_n_bands(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),1))) ;  if (error()) goto 100

        !---------------------------------------------------------------------------------------------
        !-------------- Restart from a file ----------------------------------------------------------
        case default   

           call my(tagio(trim(params%restart_filename),TAGIO_READ,mkey,len(mkey)),restf) ; if (error()) goto 100
           if (i_access(restf)) ios = x_tagfd(restf)
           if (i_comm(restf)) call broadcast(MOD_SCOPE,ios)
           if (error(ios == 0,"ERROR: restart file was not found")) goto 100


           if (i_access(restf)) tios = findfirsttag(restf,"EHRENFEST")
           if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
           if (error(tios /= TAG_START_BLOCK,"ERROR: EHRENFEST block was not found")) goto 100
           
           if (i_access(restf)) call openblock(restf)

           !** Read in the state data
           if (i_access(restf)) tios = findfirsttag(restf,"STATE")
           if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
           if (error(tios == TAG_NOT_FOUND,"ERROR: STATE tag was not found in restart file")) goto 100

           ! Read the state
           if (i_access(restf)) then
              
              dsize = sizeof_double ; ndata = 1
              call readf(state%time,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)

              dsize = sizeof_long ; ndata = 1
              call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
              state%istep = s4
              call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
              state%atoms_istep = s4

              dsize = sizeof_double ; ndata = 1
              call readf(state%last_partial_charge,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
              call readf(state%last_partial_time,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)

              dsize = sizeof_long ; ndata = 1
              call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
              state%num_sc_steps = s4
              
           end if


           !** Broadcast the state info to the rest of the processors
           if (i_comm(restf)) call broadcast(FILE_SCOPE,state%time)
           if (i_comm(restf)) call broadcast(FILE_SCOPE,state%istep)
           if (i_comm(restf)) call broadcast(FILE_SCOPE,state%atoms_istep)
           if (i_comm(restf)) call broadcast(FILE_SCOPE,state%last_partial_charge)
           if (i_comm(restf)) call broadcast(FILE_SCOPE,state%last_partial_time)
           if (i_comm(restf)) call broadcast(FILE_SCOPE,state%num_sc_steps)

           !** Construct a TDDFT stepper object from the restart file
           call my(config_td(restf=restf),state%cfg_td)             ; if (error()) goto 100

!           call my(ion_dynamics(state%cfg_td, params%ions_dt, restf),state%iondyn) ; if (error()) goto 100
           !call my(x_external(state%cfg_td),state%ext)        ; if (error()) goto 100

!           call my(external(restf),state%ext)        ; if (error()) goto 100


           nk = x_n_kpoints(x_kpoints(x_electrons(state%cfg_td)))  ; if (error()) goto 100
           nb = x_n_bands(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),1))) ;  if (error()) goto 100


           if (i_access(restf)) call closeblock(restf)



           call glean(thy(restf))


        end select
        !*********************************************************************************************
        !*********************************************************************************************

        if (debug) call warn("init_state_i:: call update(cfg_td...) ")

        !** update the config object with the cutoffs. Note that if a cutoff
        !      is negative, nothing will be changed in the underlying hamiltonian
        !      or in the construction of the density.
!RMH - disabled for now. Need to fix...
!        call update(state%cfg_td,params%ke_cutoff,params%den_cutoff,params%den_beta) 
        if (error()) goto 100

100     if (error("Exit ehrenfest_mod::init_state_i")) continue


      end subroutine init_state_i


      !********************************************************************
      !** 
      !**  init_files_i
      !**
      !** This routine is responsible for initializing the files used for
      !**  outputting info
      subroutine init_files_i( files, params )
        type(ehrenfest_files) :: files
        type(ehrenfest_params) :: params

        !local vars --------------
        logical            :: debug
        logical            :: file_exists
        !-----------------------------------------------------

        debug = .false.
        
        if (debug) call warn("init_files_i:: starting ")

        !*********************************************************************************************
        !** Check to see if this is a restart
!        select case(params%restart_filename(1:len_trim(params%restart_filename)))
        select case(trim(params%restart_filename))
        case ("none","off","false") !** Starting (i.e. not restarting)
        
           !------------------------------------------------------------------------------------------
           !-------------- Initialize files from scratch ---------------------------------------------
           !** Write the current, dipole moment and normalization to tddft_out
           call my(file('tddft_out'), files%tddft_out)
           file_exists = open_file_i(files%tddft_out,append=.false.)
           if (error("Error opening the tddft_out file")) goto 100

           if (i_access(files%tddft_out)) then
              write(x_unit(files%tddft_out),"(a6,a16,a13,5a20,a26)") &
                   "% Step","time(fs)","Ix","Iy","Iz","px","py","pz","Normalization"
           end if
           
           !** Write some field quantities to 'fields_out'
           if (debug) call warn("init_ehrenfest_i:: write out the vector potential, etc... ")
           call my(file('fields_out'), files%fields_out)
           file_exists = open_file_i(files%fields_out,append=.false.)
           if (error("Error opening the fields_out file")) goto 100

           if (i_access(files%fields_out)) then
              write(x_unit(files%fields_out),"(a6,a16,a11,5a16)") "% Step","time(fs)", &
                   " Aext_x ", " Aext_y ", " Aext_z ", " Eext_x ", " Eext_y ", " Eext_z "
           end if

           !** Write the energy terms to 'energy_out'
           call my(file('energy_out'), files%energy_out)
           file_exists = open_file_i(files%energy_out,append=.false.)
           if (error("Error opening the energy_out file")) goto 100

           if (i_access(files%energy_out)) then
              write(x_unit(files%energy_out),"(a8,a14,a13,a18,a22)") &
                   "%   Step","time(fs)","Total","Electronic(KS)","Ion Kinetic"
           end if
!           write(x_unit(files%energy_out),"(i10,f12.5,3e18.10)") &
!                state%istep, state%time, cell_energy + ion_kinetic_energy, cell_energy, ion_kinetic_energy

           !** Initialize the projection output file
           if (0 < params%proj_interval) then
              call my(file('projections_out'), files%proj_out)
              file_exists = open_file_i(files%proj_out,append=.false.)
              if (error("Error opening the projections_out file")) goto 100

              if (i_access(files%proj_out)) then
                 write(x_unit(files%proj_out),"(a6,a16,a16)") "% Step","time(fs)", "Num Electrons"
              end if

              call my(file('projections_curr_out'), files%proj_curr_out)
              file_exists = open_file_i(files%proj_curr_out,append=.false.)
              if (error("Error opening the projections_curr_out file")) goto 100

              if (i_access(files%proj_curr_out)) then
                 write(x_unit(files%proj_curr_out),"(a6,a16,a16)") "% Step","time(fs)", "Num Electrons"
              end if
           end if

           !** Initialize the partial_charge output file
           if (0 < params%partial_charge_interval) then
              call my(file('partial_charge_out'), files%partial_charge_out)
              file_exists = open_file_i(files%partial_charge_out,append=.false.)
              if (error("Error opening the partial_charge_out file")) goto 100

              if (i_access(files%partial_charge_out)) then
                 write(x_unit(files%partial_charge_out),"(a6,a16,2a16)") "% Step","time(fs)", "Charge", "Current (dQ/dt)"
              end if
           end if

           !** Initialize the eigenvalue output file
           if (0 < params%eigval_interval) then
              call my(file('eigval_out'), files%eigval_out)
              file_exists = open_file_i(files%eigval_out,append=.false.)
              if (error("Error opening the eigval_out file")) goto 100

              if (i_access(files%eigval_out)) then
                 write(x_unit(files%eigval_out),"(a6,a16,a16)") "% Step","time(fs)", "Eigenvalues"
              end if
           end if


           !** Initialize the orthogonalization output file
           if (0 < params%orth_interval) then
              call my(file('orth_out'), files%orth_out)
              file_exists = open_file_i(files%orth_out,append=.false.)
              if (error("Error opening the orth_out file")) goto 100

              if (i_access(files%orth_out)) then
                 write(x_unit(files%orth_out),"(a6,a16,a16)") "% Step","time(fs)", "Non-orthogonal overlap"
              end if
           end if

        !---------------------------------------------------------------------------------------------
        !-------------- Restart from a file ----------------------------------------------------------
        case default   

          !** Open tddft_out in append mode
           call my(file('tddft_out'), files%tddft_out)
           file_exists = open_file_i(files%tddft_out,append=.true.)
           if (error("Error opening the tddft_out file in append mode")) goto 100
           
           !** open fields_out in append mode
           call my(file('fields_out'), files%fields_out)
           file_exists = open_file_i(files%fields_out,append=.true.)
           if (error("Error opening the fields_out file in append mode")) goto 100

           !** open energy_out in append mode
           call my(file('energy_out'), files%energy_out)
           file_exists = open_file_i(files%energy_out,append=.true.)
           if (error("Error opening the energy_out file in append mode")) goto 100

           !** open proj_out in append mode
           if (0 < params%proj_interval) then
              call my(file('projections_out'), files%proj_out)
              file_exists = open_file_i(files%proj_out,append=.true.)
              if (error("Error opening the proj_out file in append mode")) goto 100
              call my(file('projections_curr_out'), files%proj_curr_out)
              file_exists = open_file_i(files%proj_curr_out,append=.true.)
              if (error("Error opening the proj_curr_out file in append mode")) goto 100
           end if

           !** open partial_charge_out in append mode
           if (0 < params%partial_charge_interval) then
              call my(file('partial_charge_out'), files%partial_charge_out)
              file_exists = open_file_i(files%partial_charge_out,append=.true.)
              if (error("Error opening the partial_charge_out file in append mode")) goto 100
           end if

           !** open eigval_out in append mode
           if (0 < params%eigval_interval) then
              call my(file('eigval_out'), files%eigval_out)
              file_exists = open_file_i(files%eigval_out,append=.true.)
              if (error("Error opening the eigval_out file in append mode")) goto 100
           end if

           !** open orth_out in append mode
           if (0 < params%orth_interval) then
              call my(file('orth_out'), files%orth_out)
              file_exists = open_file_i(files%orth_out,append=.true.)
              if (error("Error opening the orth_out file in append mode")) goto 100
           end if

        end select
        !*********************************************************************************************
        !*********************************************************************************************


100     if (error("Exit ehrenfest_mod::init_files_i")) continue

      end subroutine 


      !** Close down the files in the files object
      subroutine close_files_i(files,params)
        type(ehrenfest_files)  :: files
        type(ehrenfest_params) :: params
        logical :: debug
        !-----------------------

        debug = .false.

        if (debug) call warn("ehrenfest::close_files_i - starting")
        if (i_access(files%tddft_out)) close(x_unit(files%tddft_out))

        if (i_access(files%fields_out)) close(x_unit(files%fields_out))
        if (i_access(files%energy_out)) close(x_unit(files%energy_out))
        if (0 < params%proj_interval) then
           if (i_access(files%proj_out)) close(x_unit(files%proj_out))
           if (i_access(files%proj_curr_out)) close(x_unit(files%proj_curr_out))
        end if

        if (0 < params%partial_charge_interval) then
           if (i_access(files%partial_charge_out)) close(x_unit(files%partial_charge_out))
        end if
        if (0 < params%eigval_interval) then
           if (i_access(files%eigval_out)) close(x_unit(files%eigval_out))
        end if
        if (0 < params%orth_interval) then
           if (i_access(files%orth_out)) close(x_unit(files%orth_out))
        end if

        if (debug) call warn("ehrenfest::close_files_i - gleaning")

        call glean(thy(files%tddft_out))   ; if (error()) goto 100
        call glean(thy(files%fields_out))  ; if (error()) goto 100
        call glean(thy(files%energy_out))  ; if (error()) goto 100
        if (0 < params%proj_interval) call glean(thy(files%proj_out))    ; if (error()) goto 100
        if (0 < params%proj_interval) call glean(thy(files%proj_curr_out)) ; if (error()) goto 100
        if (0 < params%partial_charge_interval) call glean(thy(files%partial_charge_out)) ; if (error()) goto 100
        if (0 < params%eigval_interval) call glean(thy(files%eigval_out))  ; if (error()) goto 100
        if (0 < params%orth_interval) call glean(thy(files%orth_out))    ; if (error()) goto 100

100     if (error("Exit ehrenfest_mod::close_files_i")) continue

        if (debug) call warn("ehrenfest::close_files_i - exiting")

      end subroutine close_files_i


      subroutine close_state_i(state,params)
        type(ehrenfest_state)  :: state
        type(ehrenfest_params) :: params

        integer :: ik

        call glean(thy(state%cfg_td))
        call glean(thy(state%iondyn))
        call glean(thy(state%ext))

100     if (error("Exit ehrenfest_mod::close_state_i")) continue

      end subroutine


      !** This routine updates the total external potential including both the atomic positions and 
      !       any external free fields that depend explicitly on time.
      subroutine update_ext_i(state, params)
        type(ehrenfest_state)  :: state
        type(ehrenfest_params) :: params
        ! Local vars
        !------------------
        integer :: na, ia
        real(double)      :: pos_cart(3)
        real(double)      :: pos_lat(3)
        type(atoms_obj)   :: at
        type(lattice_obj) :: lat
        type(crystal_obj) :: cr

        call my(x_crystal(state%ext),cr)
        call my(x_atoms(cr),at)
        call my(x_lattice(cr),lat)

        ! Move atoms
        if (x_method(state%iondyn) /= NONE) then
           if (mod(state%istep,params%move_atoms_interval) == 0) then
              if (move_atoms(state%iondyn, state%cfg_td)) then
                 if (error("Move atoms failed")) goto 100
!                 if (mpi_first(MOD_SCOPE)) write(*,*) '**md_step ', x_istep(state%iondyn)
                 na = x_n_atoms(at); if (error()) goto 100

                 do ia=1,na
                    pos_cart = x_pos_cart(state%iondyn,ia); if (error()) goto 100
                    pos_lat = r2lat(lat,pos_cart); if (error()) goto 100
                    call move(at,ia,pos_lat); if (error()) goto 100
                 end do

                 call update(cr,atoms=at); if (error()) goto 100
                 call update(state%ext,cr)

              end if
           end if
           if (error("Move atoms failed")) goto 100
        end if


        call glean(thy(at))
        call glean(thy(lat))
        call glean(thy(cr))
100     return

      end subroutine


      !** dump_output_i is responsible for writing out tddft related info 
      !     to all the various files as well as to the console as the 
      !     calculation progresses.
      subroutine dump_output_i(state,params,files,disable_restart)
        type(ehrenfest_state)  :: state
        type(ehrenfest_params) :: params
        type(ehrenfest_files) :: files
        logical, optional      :: disable_restart

        ! Local Vars ------------------------------------------------------------
        type(tagio_obj)        :: nrestf
        type(grid_obj)         :: g
        type(grid_obj)         :: gx,gy,gz
        type(electrons_td_obj) :: el
        character(line_len)     :: filename
        real(double)           :: norm
        integer                :: num_hpsis
        logical                :: restart_enabled
        logical                :: debug
        !========================================================================

        debug = .false.

        restart_enabled = .true.
        if (present(disable_restart)) restart_enabled = .not.disable_restart

        !** Dump data to tddft_out and fields_out
        if (debug) call warn("dump_output_i:: dump_tddft_fields_out_i ")
        if (0 < params%file_dump_interval) then
           if (mod(state%istep,params%file_dump_interval) == 0) then
              call dump_tddft_fields_out_i(state,params,files) ; if (error()) goto 100
           end if
        end if

        if (debug) call warn("dump_output_i:: printing progress out to screen ")
        
        !** print out progress to console
        if (0 < params%screen_dump_interval) then
           if (mod(state%istep,params%screen_dump_interval) == 0) then
              if (debug) call warn("dump_output_i:: calling get_norm")
              norm = get_norm(x_electrons(state%cfg_td))
              if (debug) call warn("dump_output_i:: setting num_hpsis")

              num_hpsis = num_hamiltonian_ops(x_electrons(state%cfg_td))

              if (mpi_first(MOD_SCOPE)) write(x_unit(output),"(a8,i8,a6,f9.4,a5,f20.16,a7,i4,a9,i4)") &
                   "td_step ", state%istep, &
                   "  time", state%time*ARU_2_FS, &
                   "  |v|", norm, &
                   " H|v>s:", num_hpsis, "sc steps", state%num_sc_steps

              if (debug) call warn("dump_output_i:: finished write(x_unit)")
           end if
        end if

        if (debug) call warn("dump_output_i:: writing out a restart file (maybe) ")
        
        !** Write out a restart file if enough time steps have been taken
        if (0 < params%restf_interval .and. restart_enabled) then
           if (mod(state%istep,params%restf_interval) == 0) then
              call write_restart_i(state, params)
              !** Flush the i/o buffers so that the results get written to file immediately
              call flushbuf(files%tddft_out) 
              call flushbuf(files%fields_out)
              call flushbuf(files%energy_out)
           end if
        end if

!!!!RMH this is broken
        !** Write out the grid density to file if enough steps have been taken
!        if (0 < params%den_interval) then
!           if (mod(state%istep,params%den_interval) == 0) then
!              call get_density_filename_i(filename,params%mesh_file_type,state%istep)
!              if (overlap_is_identity(x_atomic_potential(x_potential(x_fields(state%cfg_td))))) then
!                 call my(x_grid_density(x_density(x_electrons(state%cfg_td))),g)
!              else
!                 call my(all_electron_grid_density(x_electrons(state%cfg_td)),g)
!              end if
!              if (error()) goto 100
!              call write_to_file(g,filename,params%mesh_file_type,RS_KIND) ; if (error()) goto 100
!              call glean(thy(g))
!           end if
!        end if
        
!!!!RMH this is broken
        !** Write out the grid current density to file if enough steps have been taken
!        if (0 < params%currden_interval) then
!           if (mod(state%istep,params%currden_interval) == 0) then
!              call get_current_density_filename_i(filename,params%mesh_file_type,state%istep)

!              call my(x_grid_current_density(x_density(x_electrons(state%cfg_td)),1),gx)
!              call my(x_grid_current_density(x_density(x_electrons(state%cfg_td)),2),gy)
!              call my(x_grid_current_density(x_density(x_electrons(state%cfg_td)),3),gz)
!              if (error()) goto 100
!              call write_to_file(gx,gy,gz,filename,params%mesh_file_type,RS_KIND) ; if (error()) goto 100
!              call glean(thy(gx))
!              call glean(thy(gy))
!              call glean(thy(gz))
!           end if
!        end if

        if (debug) call warn("dump_output_i:: writing out local potential ")

        !** Write out the local potential to file if enough steps have been taken
!        if (0 < params%pot_interval) then
!           if (mod(state%istep,params%pot_interval) == 0) then
!              call get_potential_filename_i(filename,params%mesh_file_type,state%istep)
!              call my(local_potential(x_fields(state%cfg_td)),g)
!              if (error()) goto 100
!              call write_to_file(g,filename,params%mesh_file_type,RS_KIND) ; if (error()) goto 100
!              call glean(thy(g))
!           end if
!        end if


!!!!RMH This is broken - it's not printing out the harpot is printing out the total local potential !!!!        
        !** Write out the hartree potential to file if enough steps have been taken
!        if (0 < params%hap_interval) then
!           if (mod(state%istep,params%hap_interval) == 0) then
!              call get_harpot_filename_i(filename,params%mesh_file_type,state%istep)
!              call my(grid(x_layout(local_potential(x_fields(state%cfg_td)))),g)
!              call put_hap(x_fields(state%cfg_td),g)
!              if (error()) goto 100
!              call write_to_file(g,filename,params%mesh_file_type,RS_KIND) ; if (error()) goto 100
!              call glean(thy(g))
!           end if
!        end if

        !** Write out the wavefunctions to file if enough steps have been taken
!        if (0 < params%wf_interval) then
!           if (mod(state%istep,params%wf_interval) == 0) then
!              call my(x_electrons(state%cfg_td), el)
!              call write_wavefunctions(el, params%mesh_file_type, state%istep ) 
!              if (error()) goto 100
!              call glean(thy(el))
!           end if
!        end if

        if (debug) call warn("dump_output_i:: writing out the projection of the propagated wavefunctions... ")

        !** calculate and write out the projection of the propagated wavefuncs with the original wavefuncs
        if (0 < params%proj_interval) then
           if (mod(state%istep,params%proj_interval) == 0) then
              call write_proj_out_i(state,files)
           end if
        end if

        if (debug) call warn("dump_output_i:: calculating partial charge ")

        !** Find the amount of charge on the positive side of the surface defined by surface_origin
        !    and surface_direction tags.  Also calculate the current.
        if (0 < params%partial_charge_interval) then
           if (mod(state%istep,params%partial_charge_interval) == 0) then
              call write_partial_charge_out_i(state,params,files)
           end if
        end if

        if (debug) call warn("dump_output_i:: dumping eigenvalues ")

        !** Calculate and write out the eigenvalues 
        if (0 < params%eigval_interval) then
           if (mod(state%istep,params%eigval_interval) == 0) then
              call write_eigvals_i(state,files)
           end if
        end if

        if (debug) call warn("dump_output_i:: dumping orthogonalization ")
        !** Calculate and write out the eigenvalues 
        if (0 < params%orth_interval) then
           if (mod(state%istep,params%orth_interval) == 0) then
              call write_orth_i(state,files)
           end if
        end if

100     if (error("Exit tddft_mod::dump_output_i")) continue

      end subroutine dump_output_i
      

      subroutine dump_tddft_fields_out_i(state,params,files)
        type(ehrenfest_state)        :: state
        type(ehrenfest_params)       :: params
        type(ehrenfest_files)       :: files
        
        ! Local Vars
        real(double) :: current(3)
        real(double) :: electric_field(3)
        real(double) :: dEdt(3)
        real(double) :: avg_vp_tot(3)
        real(double) :: rms_vp_tot(3)
        real(double) :: avg_vp_har(3)
        real(double) :: rms_vp_har(3)
        real(double) :: avg_vp_xc(3)
        real(double) :: rms_vp_xc(3)
        real(double) :: A_ext(3)
        real(double) :: dipole_mom(3)
        real(double) :: norm

        real(double) :: kinetic_energy
        real(double) :: har_scalar_energy
        real(double) :: har_vecpot_energy
        real(double) :: xc_energy
        real(double) :: scp_energy
        real(double) :: el_ion_energy
        real(double) :: ion_ion_energy
        real(double) :: cell_energy
        real(double) :: eig_energy
        real(double) :: ion_kinetic_energy
        real(double) :: electromagnetic_energy
        real(double) :: joule_p
        logical      :: debug 

        debug = .false.

        call my(files%tddft_out)
        call my(files%fields_out)
        call my(files%energy_out)

        if (debug) call warn("dump_fields_tddft_i:: starting ")

        !** Calculate various quantities we may be interested in.
!        current = total_current(x_electrons(state%cfg_td),state%time) ; if (error()) goto 100
!        current = total_current(x_electrons(state%cfg_td)) ; if (error()) goto 100
!        dipole_mom = dipole_moment(x_density(x_electrons(state%cfg_td))) ; if (error()) goto 100

!!!RMH
electric_field = 0.0_double
!        electric_field = x_electric_field(x_external(state%cfg_td),state%time)*ARU_2_VoltsPerAngstrom 
        if (error()) goto 100

!!!RMH
electric_field = 0.0_double
!        dEdt = x_dEdt(x_external(state%cfg_td),state%time)*ARU_2_VoltsPerAngstrom/ARU_2_FS  
        if (error()) goto 100

!        avg_vp_tot = avg_vecpot_tot(cfg)*ARU_2_VoltsPerAngstrom*ARU_2_FS ; if (error()) goto 100
!        rms_vp_tot = rms_vecpot_tot(cfg)*ARU_2_VoltsPerAngstrom*ARU_2_FS ; if (error()) goto 100
!        avg_vp_har = avg_vecpot_har(cfg)*ARU_2_VoltsPerAngstrom*ARU_2_FS ; if (error()) goto 100
!        rms_vp_har = rms_vecpot_har(cfg)*ARU_2_VoltsPerAngstrom*ARU_2_FS ; if (error()) goto 100
!        avg_vp_xc  = avg_vecpot_xc(cfg)*ARU_2_VoltsPerAngstrom*ARU_2_FS ; if (error()) goto 100
!        rms_vp_xc  = rms_vecpot_xc(cfg)*ARU_2_VoltsPerAngstrom*ARU_2_FS ; if (error()) goto 100

!!!RMH
A_ext = 0.0_double
!        A_ext = x_external_vecpot(x_external(state%cfg_td),state%time)*ARU_2_VoltsPerAngstrom*ARU_2_FS 
        if (error()) goto 100


        norm = get_norm(x_electrons(state%cfg_td)) ; if (error()) goto 100
        cell_energy = x_cell_energy(state%cfg_td) ; if (error()) goto 100
        if (x_method(state%iondyn) /= NONE) then
           ion_kinetic_energy = x_kinetic_energy(state%iondyn) 
        else
           ion_kinetic_energy = 0.0_double
        end if


!!!!!RMH - Can add additional energy terms if necessary
!        kinetic_energy = x_kinetic_energy(state%cfg_td) ; if (error()) goto 100
!        har_scalar_energy = x_hartree_energy(state%cfg_td) ; if (error()) goto 100
!        har_vecpot_energy = hartree_vecpot_energy(state%cfg_td) ; if (error()) goto 100
!        xc_energy = x_exc_energy(state%cfg_td)          ; if (error()) goto 100
!        scp_energy = x_scp_energy(state%cfg_td)          ; if (error()) goto 100
!        el_ion_energy = x_el_ion_energy(state%cfg_td)   ; if (error()) goto 100
!        ion_ion_energy = x_ion_ion_energy(state%cfg_td) ; if (error()) goto 100
        
!        eig_energy = x_eig_energy(state%cfg_td) ; if (error()) goto 100
!        electromagnetic_energy = em_energy(state%cfg_td) ; if (error()) goto 100
        
!        joule_p =g joule_power(state%cfg_td) ; if (error()) goto 100


        !** Dump the current, dipole moment and normalization to files%tddft_out
        if (debug) call warn("dump_fields_tddft_i:: dumping current, dipole moment and norm")
current = 0.0
dipole_mom = 0.0
        if (i_access(files%tddft_out)) then
           write(x_unit(files%tddft_out),"(i10,f12.5,6e20.11,f20.16)") &
                state%istep,state%time*ARU_2_FS,current,dipole_mom,norm
        end if

        if (debug) call warn("dump_fields_tddft_i:: dumping vecpot elec field ... ")
        
        !** Dump the Vector Potential, Electric Field and Derivative of the E field to files%fields_out
        if (i_access(files%fields_out)) then
           write(x_unit(files%fields_out),"(i10,f12.5,6e16.7)") state%istep,state%time*ARU_2_FS, &
                A_ext, electric_field !,rms_vp_tot,ext_vp,avg_vp_har,rms_vp_har,avg_vp_xc,rms_vp_xc
        end if

        !** Dump the various energy terms to energy out
        if (i_access(files%energy_out)) then
           write(x_unit(files%energy_out),"(i10,f12.5,3e18.10)") &
                state%istep, state%time*ARU_2_FS, cell_energy + ion_kinetic_energy, cell_energy, ion_kinetic_energy
        end if
        !!!!RYAN Can write out all the energy terms if desired
!        if (i_access(files%energy_out)) then
!           write(x_unit(files%energy_out),"(i10,f12.5,10e18.10)") &
!                state%istep, state%time, cell_energy + ion_kinetic_energy, kinetic_energy,&
!                har_scalar_energy, xc_energy, el_ion_energy, ion_ion_energy, &
!                cell_energy, ion_kinetic_energy, scp_energy, eig_energy
!        end if


        !** Flush the i/o buffers so that the results get written to file immediately
        call flushbuf(files%tddft_out)
        call flushbuf(files%fields_out)
        call flushbuf(files%energy_out)

!write(*,*) 'called flushbuf'

100     call glean(thy(files%tddft_out))
        call glean(thy(files%fields_out))
        call glean(thy(files%energy_out))


        if (error('Exiting tddft::dump_tddft_fields_out_i')) continue
!write(*,*) 'exiting dump_tddft_fields_out_i...'

      end subroutine dump_tddft_fields_out_i


      !** Write out a restart file.
      subroutine write_restart_i(state,params)
        type(ehrenfest_state)  :: state
        type(ehrenfest_params) :: params
        !**  Local vars
        type(tagio_obj)        :: nrestf
        character(line_len)     :: filename
        integer                :: ik, nk
        integer(long)          :: dsize, ios, ndata, s4
        !------------------------------------

        !** Construct the restart file
        call get_restart_filename_i(filename,state%istep)
        call my(tagio(trim(filename),TAGIO_WRITE,mkey,len(mkey)),nrestf)

        !** Define a tddft block
        if (i_access(nrestf)) call startblock(nrestf,"EHRENFEST") 

        !** Write out the contents of the tddft state
        if (i_access(nrestf)) then
           call writetag(nrestf,"STATE")

           dsize = sizeof_double ; ndata = 1
           call writef(state%time,dsize,ndata,x_tagfd(nrestf),ios)           

           dsize = sizeof_long ; ndata = 1
           s4 = state%istep
           call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)   
           s4 = state%atoms_istep
           call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)   

           dsize = sizeof_double ; ndata = 1
           call writef(state%last_partial_charge,dsize,ndata,x_tagfd(nrestf),ios)           
           call writef(state%last_partial_time,dsize,ndata,x_tagfd(nrestf),ios)           

           dsize = sizeof_long ; ndata = 1
           s4 = state%num_sc_steps
           call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)   

        end if
           
        call write_restart(state%cfg_td,nrestf) ; if (error("Error writing cfg_td restart")) goto 100
        call write_restart(state%iondyn,nrestf) ; if (error("Error writing iondyn restart")) goto 100
        call write_restart(state%ext,nrestf) ; if (error("Error writing ext restart")) goto 100
        
        if (i_access(nrestf)) call endblock(nrestf) 

100     call glean(thy(nrestf))
        
      end subroutine


      subroutine get_restart_filename_i(fullname, i)
        character(line_len) :: fullname
        integer :: i

        character(line_len) :: num
        character(line_len) :: filename
        character(line_len) :: path
        integer :: j, pos, a
        
        write(num,"(I8)") i
                                  
        !path = "/scratch/rhatcher/"
        filename = "tddft_restartf_        "
        pos = 16
        !write(*,*) "in get_restart_filename_i"
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

!        fullname = trim(path)//trim(filename)
        fullname = trim(filename)

      end subroutine


      subroutine get_density_filename_i(filename,file_type,i)
        character(line_len) :: filename
        integer            :: file_type
        integer            :: i

        select case(file_type)
        case (AMIRA)
           call get_density_am_filename_i(filename,i)
        case(VTK)
           call get_density_vtk_filename_i(filename,i)
        case default
           call warn("tddft::get_density_filename_i - unrecognized file_type, reverting to AMIRA")
           call get_density_am_filename_i(filename,i)
        end select

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

      subroutine get_current_density_filename_i(filename,file_type,i)
        character(line_len) :: filename
        integer            :: file_type
        integer            :: i

        select case(file_type)
        case (AMIRA)
           call get_currden_am_filename_i(filename,i)
        case(VTK)
           call get_currden_vtk_filename_i(filename,i)
        case default
           call warn("tddft::get_current_density_filename_i - unrecognized file_type, reverting to AMIRA")
           call get_currden_am_filename_i(filename,i)
        end select

      end subroutine

      subroutine get_currden_am_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "currentdensity_        .am"
        pos = 16
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

      subroutine get_currden_vtk_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "currentdensity_        .vtk"
        pos = 16
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


      subroutine get_potential_filename_i(filename,file_type,i)
        character(line_len) :: filename
        integer            :: file_type
        integer            :: i

        select case(file_type)
        case (AMIRA)
           call get_potential_am_filename_i(filename,i)
        case(VTK)
           call get_potential_vtk_filename_i(filename,i)
        case default
           call warn("tddft::get_potential_filename_i - unrecognized file_type, reverting to AMIRA")
           call get_potential_am_filename_i(filename,i)
        end select

      end subroutine

      subroutine get_potential_am_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "potential_        .am"
        pos = 11
        !write(*,*) "in get_potential_filename_i"
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

      subroutine get_potential_vtk_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "potential_        .vtk"
        pos = 11 
        !write(*,*) "in get_potential_filename_i"
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

      subroutine get_harpot_filename_i(filename,file_type,i)
        character(line_len) :: filename
        integer            :: file_type
        integer            :: i

        select case(file_type)
        case (AMIRA)
           call get_harpot_am_filename_i(filename,i)
        case(VTK)
           call get_harpot_vtk_filename_i(filename,i)
        case default
           call warn("tddft::get_harpot_filename_i - unrecognized file_type, reverting to AMIRA")
           call get_harpot_am_filename_i(filename,i)
        end select

      end subroutine

      subroutine get_harpot_am_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "hartree_pot_        .am"
        pos = 13
        !write(*,*) "in get_potential_filename_i"
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

      subroutine get_harpot_vtk_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

        filename = "hartree_pot_        .vtk"
        pos = 13
        !write(*,*) "in get_potential_filename_i"
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

      subroutine read_occupations_i(filename,nk,nb,new_occs)
        character(line_len)    :: filename
        integer               :: nk
        integer               :: nb
        real(double), pointer :: new_occs(:)
        
        ! local vars
        type(file_obj) :: occ_file
        integer        :: ik, ib, ios
        integer        :: start, finish
        
        call my(file(trim(filename)), occ_file)
        
        if (i_access(occ_file)) open(unit=x_unit(occ_file), &
                file=x_name(occ_file),iostat=ios)
        
        if (i_comm(occ_file)) call broadcast(MOD_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open occupation file")) goto 200
        
        new_occs = 0.0_double

        if (i_access(occ_file)) then
           do ik=1,nk
              start = (ik-1)*nb + 1
              finish = ik*nb
              read(x_unit(occ_file),*) new_occs(start:finish)
           end do
        end if

        if (i_comm(occ_file)) call broadcast(MOD_SCOPE,new_occs)
        
        if (i_access(occ_file)) close(unit=x_unit(occ_file))

100     call glean(thy(occ_file))

200     if (error('tddft_mod::read_occupations_i')) continue

      end subroutine

      subroutine write_proj_out_i(state,files)
        type(ehrenfest_state)    :: state
        type(ehrenfest_files)    :: files
        ! Local Vars -------------------------------
        integer                  :: ik,nk
        integer                  :: ib,nb
        integer                  :: beg
        real(double)             :: time
        real(double)             :: num_electrons
        real(double), pointer    :: coeffs(:)
        real(double), pointer    :: occs(:,:)
        real(double), pointer    :: eigenvalues(:,:)
        
        type(config_sc_obj)   :: cfg_sc
        complex(double), pointer :: cmat(:,:)
        type(electrons_sc_obj)   :: adiabatic_el

        !===========================================
        
        call my(state%cfg_td)
        
        !** Find the number of k-points and the number of bands
        nk = x_n_kpoints(x_kpoints(x_electrons(state%cfg_td)))
        nb = x_n_bands(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),1)))
        
        !** Allocate cmat, occs, eigenvalues, and coeffs
        allocate(cmat(nb,nb),occs(nk,nb),eigenvalues(nk,nb),coeffs(nb*nk))
        
        !** Get the occupations
        occs = x_occupations(x_electrons(state%cfg_td))
        
        !** Get the self-consistent config for the current external
        call my(x_config_sc(state%cfg_td),cfg_sc)

        eigenvalues = x_eigenvalues(x_electrons(cfg_sc))

        !** Calc the overlap and multiply by the occupation.
        do ik=1,nk
           call overlap(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),ik)), &
                        x_multivector(x_wavefunctions(x_electrons(cfg_sc),ik)),cmat) ; if (error()) goto 100
           beg = (ik-1)*nb 
           do ib=1,nb
              coeffs(beg+ib) = sum(occs(ik,:)*conjg(cmat(ib,:))*cmat(ib,:))
           end do
        end do

        !** clean up the self-consistent config
        call glean(thy(cfg_sc))

        num_electrons = sum(coeffs)

        !** Write the coefficients out to the proj_out file
        time = state%time*ARU_2_FS
        call my(files%proj_out) ; if (error()) goto 200

        if (i_access(files%proj_out)) then
           write(x_unit(files%proj_out),*) state%istep,time,num_electrons
           do ik=1,nk
              beg = (ik-1)*nb
              do ib=1,nb
                 write(x_unit(files%proj_out),*) eigenvalues(ik,ib), coeffs(beg+ib)
              end do
           end do
        end if

        call flushbuf(files%proj_out)

        !** Get the electrons containing the eigenstates for the current TDDFT Hamiltonian
        call my(x_adiabatic_states(state%cfg_td),adiabatic_el)

        eigenvalues = x_eigenvalues(adiabatic_el)

        !** Calc the overlap and multiply by the occupation.
        do ik=1,nk
           call overlap(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),ik)), &
                        x_multivector(x_wavefunctions(adiabatic_el,ik)),cmat) ; if (error()) goto 100
           beg = (ik-1)*nb 
           do ib=1,nb
              coeffs(beg+ib) = sum(occs(ik,:)*conjg(cmat(ib,:))*cmat(ib,:))
           end do
        end do

        !** clean up the adiabatic electrons
        call glean(thy(adiabatic_el))

        num_electrons = sum(coeffs)

        !** Write the coefficients out to the proj_curr_out file
        time = state%time*ARU_2_FS
        call my(files%proj_curr_out) ; if (error()) goto 300
        if (i_access(files%proj_curr_out)) then
           write(x_unit(files%proj_curr_out),*) state%istep,time,num_electrons
           do ik=1,nk
              beg = (ik-1)*nb
              do ib=1,nb
                 write(x_unit(files%proj_curr_out),*) eigenvalues(ik,ib), coeffs(beg+ib)
              end do
           end do
        end if

        call flushbuf(files%proj_curr_out)

        !** Deallocate cmat, occs, eigenvalues, and coeffs
        deallocate(cmat, occs, eigenvalues, coeffs)

300     call glean(thy(files%proj_out))
200     call glean(thy(files%proj_curr_out))

100     if (error("Exit tddft_mod::write_proj_out_i")) continue
        
      end subroutine write_proj_out_i



      subroutine write_orth_i(state,files)
        type(ehrenfest_state)    :: state
        type(ehrenfest_files)    :: files
        ! Local Vars -------------------------------
        integer                  :: ik,nk,j,a,pos
        integer                  :: ib,jb,nb
        integer                  :: beg
        integer                  :: ehrenfest_step
        real(double)             :: time
        real(double)             :: total_overlap
        real(double), pointer    :: coeffs(:)
!        real(double), pointer    :: occs(:,:)
        
        type(multivector_obj)    :: mvec
        complex(double), pointer :: cmat(:,:)
        character(line_len)       :: fmt
        character(line_len)       :: num
        !===========================================
        
        !** Find the number of k-points and the number of bands
        nk = x_n_kpoints(x_kpoints(x_electrons(state%cfg_td)))
        nb = x_n_bands(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),1)))
        
        !** Allocate cmat, occs, and coeffs
        allocate(cmat(nb,nb),coeffs(nb*nk))
        
        coeffs = 0.0_double
        !** Calc the overlap and multiply by the occupation.
        do ik=1,nk
           call my(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),ik)),mvec)
           call overlap(mvec,cmat) ; if (error()) goto 100
           beg = (ik-1)*nb 
           do ib=1,nb
              do jb=1,nb
                 if (jb /= ib) coeffs(beg+ib) = coeffs(beg+ib) + conjg(cmat(ib,jb))*cmat(ib,jb)
              end do
           end do
           call glean(thy(mvec))
        end do

        total_overlap = sum(coeffs)*0.5_double

        !** Write the coefficients out to the orth_out file
        time = state%time*ARU_2_FS
        call my(files%orth_out) ; if (error()) goto 200

        !** Construct the format statement
        fmt = "(i10,f10.3,         e12.5)"
        write(num,"(I8)") nb*nk + 1
        pos = 13
        do j=1,8
           a = iachar(num(j:j) )
           if (a /= 32) fmt(pos:pos) = num(j:j)
           pos = pos + 1
        end do

!write(*,*) 'write_ortho, fmt = ', fmt

        if (i_access(files%orth_out)) then
           write(x_unit(files%orth_out),fmt) state%istep,time,total_overlap,coeffs
        end if

        call flushbuf(files%orth_out)

200     call glean(thy(files%orth_out))

100     if (error("Exit tddft_mod::write_orth_i")) continue
        
      end subroutine write_orth_i


      subroutine write_eigvals_i(state,files)

        type(ehrenfest_state)    :: state
        type(ehrenfest_files)    :: files
        ! Local Vars -------------------------------
        integer                  :: ik,nk,j,a,pos
        integer                  :: ib,nb
        integer                  :: beg
        real(double)             :: time
        real(double)             :: num_electrons
        real(double), pointer    :: evals(:)
        real(double), pointer    :: all_evals(:)
        character(line_len)       :: fmt
        character(line_len)       :: num
        logical                  :: debug
        !===========================================
        
        debug = .false.
        
        if (debug) call warn('starting write_eigvals_i')

        !** Find the number of k-points and the number of bands
        nk = x_n_kpoints(x_kpoints(x_electrons(state%cfg_td)))
        nb = x_n_bands(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),1)))
        
        !** Allocate cmat, occs, and coeffs
        allocate(evals(nb),all_evals(nb*nk))
        

        if (debug) call warn('write_eigvals_i:: getting evals..')

        !** Calc the overlap and multiply by the occupation.
        do ik=1,nk
!!!RMH
           evals = 0.0_double
!           evals = x_eigenvalues(x_wavefunctions(x_electrons(state%cfg_td),ik))

           beg = (ik-1)*nb 
           do ib=1,nb
              all_evals(beg+ib) = evals(ib)
           end do
        end do

        !** Write the coefficients out to the eigval_out file
        time = state%time*ARU_2_FS

        if (debug) call warn('write_eigvals_i:: call my(files%eigval_out)')
        call my(files%eigval_out) ; if (error()) goto 200

        !** Construct the format statement
        fmt = "(i10,f10.3,         e11.3)"
        write(num,"(I8)") nb*nk
        pos = 13
        do j=1,8
           a = iachar(num(j:j) )
           if (a /= 32) fmt(pos:pos) = num(j:j)
           pos = pos + 1
        end do

        if (debug) call warn('write_eigvals_i write(eigval_out,fmt)')

        if (i_access(files%eigval_out)) then
           write(x_unit(files%eigval_out),fmt) state%istep,time,all_evals
        end if

        if (debug) call warn('write_eigvals_i:: call flushbuf')
        call flushbuf(files%eigval_out) 

        if (debug) call warn('write_eigvals_i:: gleaning')

200     call glean(thy(files%eigval_out))

100     deallocate(evals,all_evals)

        if (error("Exit tddft_mod::write_eigvals_i")) continue
        
        if (debug) call warn('write_eigvals_i:: exiting')

      end subroutine write_eigvals_i


      subroutine write_partial_charge_out_i(state,params,files)
        type(ehrenfest_state)       :: state
        type(ehrenfest_params)      :: params
        type(ehrenfest_files)       :: files
        ! Local Vars -------------------------------
        real(double)             :: time
        real(double)             :: dt
        real(double)             :: partial_q
        real(double)             :: current
        real(double)             :: orig(3), dir(3)
        real(double),parameter   :: tol = 1.0d-6
        type(grid_obj)           :: g
        
        !===========================================
        

        !** Write the coefficients out to the proj_out file
        time = state%time*ARU_2_FS

!!!RMH Fix this
!        if (overlap_is_identity(x_atomic_potential(gen_potential(x_fields(state%cfg_td))))) then
!           call my(x_grid_density(x_density(x_electrons(state%cfg_td))),g)
!        else  ! If PAW need to call a different routine to calculate the all electron density
!           call my(all_electron_grid_density(x_electrons(state%cfg_td)),g) 
!        end if
        
        orig = params%surface_origin
        dir = params%surface_direction
!!!!RMH
!        partial_q = partial_charge(g,orig,dir) ; if (error()) goto 200

        dt = time - state%last_partial_time
        if (dt < tol) then
           dt = tol
           current = 0.0
        else
!           current = (partial_q - state%last_partial_charge)/dt
        end if

!** NB maybe we should convert to some convenient units here...        

        call my(files%partial_charge_out) ; if (error()) goto 200
        if (i_access(files%partial_charge_out)) then
!           write(x_unit(files%partial_charge_out),*) state%istep,time,partial_q,current
!           write(x_unit(files%partial_charge_out),'(i6,f10.4,f12.6,f12.6)') state%istep,time,partial_q,current
        end if

!        state%last_partial_charge = partial_q
        state%last_partial_time = time

        call flushbuf(files%partial_charge_out)

        call glean(thy(g))
200     call glean(thy(files%partial_charge_out))

100     if (error("Exit tddft_mod::write_proj_out_i")) continue
        
      end subroutine 



      !** Opens the file associated with the file_obj.  If append is true
      !     and the file exists, it will be opened in append mode.  Otherwise
      !     the file will be written over.
      function open_file_i(f, append) result(file_exists)
        type(file_obj) :: f
        logical        :: append
        logical        :: file_exists

        ! Local Vars -----------------------------------------
        integer :: ios
        !=====================================================

        call my(f)
        
        if (i_access(f)) then
           inquire(file=trim(x_name(f)),exist=file_exists)
           if (file_exists .and. append) then
              open(unit=x_unit(f), &
                   file=x_name(f), &
                   status='old', &
                   position='append', &
                   iostat=ios)
           else ! The file 'tddft_out' does not exist or we are not to append
              open(unit=x_unit(f), &
                   file=x_name(f), &
                   iostat=ios)
           end if
        end if
        if (i_comm(f)) call broadcast(MOD_SCOPE,ios)
        if (i_comm(f)) call broadcast(MOD_SCOPE,file_exists)

        if (error(ios /= 0,"ERROR: error opening file")) goto 100

100     call glean(thy(f))
        
        if (error("tddft_mod::open_file_i exiting")) continue

      end function

!      subroutine test_current(cfg)
!        type(config_tddft_obj) :: cfg
!        
!        ! local vars
!        type(grid_obj) :: gx
!        type(grid_obj) :: gy
!        type(grid_obj) :: gz
!        real(double),pointer :: rtmpx(:,:,:)
!        real(double),pointer :: rtmpy(:,:,:)
!        real(double),pointer :: rtmpz(:,:,:)
!       
!        real(double) :: Ix,Iy,Iz
!        real(double) :: current(3)
!        
!write(*,*) 'test_current starting'
!        call my(x_grid_current_density(x_density(x_electrons(state%cfg_td)), 1),gx)
!!        call my(x_grid_current_density(x_density(x_electrons(state%cfg_td)), 2),gy)
!        call my(x_grid_current_density(x_density(x_electrons(state%cfg_td)), 3),gz)
!        
!        call take(rtmpx,gx,RS_KIND) ; if (error()) goto 100
!        call take(rtmpy,gy,RS_KIND) ; if (error()) goto 100
!        call take(rtmpz,gz,RS_KIND) ; if (error()) goto 100
!        
!        Ix = sum(rtmpx)
!        Iy = sum(rtmpy)
!        Iz = sum(rtmpz)
!        
!        current = total_current(cfg)
!        if (mpi_myproc() == 0) write(*,*) ' '
!        if (mpi_myproc() == 0) write(*,*) '  Ix,Iy,Iz = ', Ix, Iy, Iz  
!        if (mpi_myproc() == 0) write(*,*) '   current = ', current(1),current(2),current(3)
!       
!        deallocate (rtmpx,rtmpy,rtmpz)
!        
!100     call glean(thy(gx))
!        call glean(thy(gy))
!        call glean(thy(gz))
!       
!        if (error('tddft_mod::test_current exiting')) continue
!      
!      end subroutine

!      subroutine write_atom_pos_i(cfg)
!        type(config_tddft_obj) :: cfg

!        type(atoms_obj)   :: at
!        type(lattice_obj) :: lat
!        integer               :: ia,na
!        real(double), pointer :: pos(:,:)

!        call my(cfg)
!        call my(x_atoms(x_crystal(x_external(state%cfg_td))),at)
!        call my(x_lattice(x_crystal(x_external(state%cfg_td))),lat)

!        na = x_natoms(at)
!        allocate(pos(3,na))

!        write(*,*) lat2r(lat,x_position(at,1)), lat2r(lat,x_position(at,2))

!        call glean(thy(cfg))
!        call glean(thy(at))
!        call glean(thy(lat))
        
!      end subroutine



!      subroutine write_wfns(cfg)
!        type(config_tddft_obj) :: cfg

!        type(multivector_obj) :: mvec
!        type(multivector_rep), pointer :: worm
!        integer :: ik,nk

!        nk = x_n_kpoints(x_kpoints(x_electrons(state%cfg_td)))

!        do ik=1,nk
!           call my(x_multivector(x_wavefunctions(x_electrons(state%cfg_td),ik)),mvec)
          
!           worm => wormhole(mvec)
           
!           write(*,*) '    tddft::wfns ', worm%mat(2,2), worm%mat(4,2), worm%mat(55,4)
!           write(*,*) 'wfns:: ', worm%mat(2,2), worm%mat(4,2), worm%mat(55,4)
!           call glean(thy(mvec))
!        end do

!      end subroutine


!      subroutine write_wfn_slice_i(cfg)
!        type(config_obj) :: cfg

!        type(multivector_obj) :: mvec
!        type(multivector_rep), pointer :: worm
!        type(grid_obj) :: g
!        complex(double), pointer :: wf(:,:,:)
!        integer :: ik,nk, ib,nb
!        character(line_len) :: filename

!write(*,*) 'write_wfn_slice_i starting...'

!        nk = x_n_kpoints(x_kpoints(x_electrons(cfg)))

!write(*,*) 'write_wfn_slice_i num kpoints = ', nk


!        do ik=1,nk

!write(*,*) 'write_wfn_slice_i calling my(mvec)'

!           call my(x_multivector(x_wavefunctions(x_electrons(cfg),ik)),mvec)

!           worm => wormhole(mvec)

!write(*,*) 'write_wfn_slice_i size(mvec) = ', size(worm%mat,1), size(worm%mat,2)

!           nb = x_n_bands(mvec)

!write(*,*) 'write_wfn_slice_i nb = ',nb

!           call my(grid(x_layout(x_multibasis(mvec))),g)

!           do ib=1,nb
!write(*,*) 'write_wfn_slice_i ib = ', ib

!              call put(mvec,ib,g)
!              if (error('error after call put(mvec,ib,g)')) goto 100
!              call get_wfn_filename_i(filename,ib)

!write(*,*) 'write_wfn_slice_i filename = ', filename

!              call transform(g,CSP_KIND)

!              call write_to_file(g,filename,AMIRA,CSP_KIND)

!              if (error('error after write_to_file')) goto 100
              
!           end do

!100        call glean(thy(mvec))
!           call glean(thy(g))
!        end do

!      end subroutine

      subroutine get_wfn_filename_i(filename, i)
        character(line_len) :: filename
        integer :: i

        character(line_len) :: num
        integer :: j, pos, a
        
        write(num,"(I8)") i

!        filename = "density_        .am"
        filename = "wfn_        .am"
        pos =5 
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

    end module

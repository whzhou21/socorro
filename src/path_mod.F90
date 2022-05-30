! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module path_mod
!doc$ module path_mod

!     path_mod: - sets and provides the default range of file units.
!               - sets directories, names, and paths for socorro files.
!               - users can customize their installation by modifying the character parameters.

!     Note: The range of file units must include 5 and 6.

      use kind_mod

!cod$
      implicit none
      private

! file units

      integer, parameter :: first_unit = 1
      integer, parameter :: last_unit  = 100

! directories

      character(line_len), parameter :: input_directory = ""  ! change to "" if Sandia build
      character(line_len), parameter :: potential_directory = ""  ! change to "" if Sandia build
      character(line_len), parameter :: output_directory = ""

! run directory file names

      character(line_len), parameter :: arg_name  = "argvf"
      character(line_len), parameter :: stop_name = "stopf"

! potential file names and file paths

      character(line_len), parameter :: ncp_name = "NCP."
      character(line_len), parameter :: paw_name = "PAW."

      character(line_len) :: ncp_path
      character(line_len) :: paw_path

! input file names and file paths

      character(line_len), parameter :: crystal_name       = "crystal"
      character(line_len), parameter :: kpoints_name       = "kpoints"
      character(line_len), parameter :: dsites_name        = "dsites"
      character(line_len), parameter :: velocity_name      = "initial_velocity"
      character(line_len), parameter :: lattice_group_name = "lattice_group"
      character(line_len), parameter :: restart_name       = "restartf"
      character(line_len), parameter :: bulk_crystal_name  = "bulk_crystal"
      character(line_len), parameter :: bulk_restart_name  = "bulk_restartf"

      character(line_len) :: crystal_path
      character(line_len) :: kpoints_path
      character(line_len) :: dsites_path
      character(line_len) :: velocity_path
      character(line_len) :: lattice_group_path
      character(line_len) :: restart_path
      character(line_len) :: neb_crystal_00_path
      character(line_len) :: neb_crystal_nn_path
      character(line_len) :: bulk_crystal_path
      character(line_len) :: bulk_restart_path

! output file names and file paths

      character(line_len), parameter :: p_error_name             = "errorf_"
      character(line_len), parameter :: f_error_name             = "errorf"
      character(line_len), parameter :: diary_name               = "diaryf"
      character(line_len), parameter :: dcomp_name               = "dcompf"
      character(line_len), parameter :: new_crystal_name         = "new_crystal"
      character(line_len), parameter :: new_kpoints_name         = "new_kpoints"
      character(line_len), parameter :: new_velocity_name        = "new_velocity"
      character(line_len), parameter :: new_lattice_group_name   = "new_lattice_group"
      character(line_len), parameter :: new_point_group_name     = "new_point_group"
      character(line_len), parameter :: new_space_group_name     = "new_space_group"
      character(line_len), parameter :: new_restart_name         = "new_restartf"
      character(line_len), parameter :: md_trajectory_name       = "md_trajectory"
      character(line_len), parameter :: tddft_md_trajectory_name = "tddft_md_trajectory"
      character(line_len), parameter :: eigenvalues_name         = "eigenvalues"
      character(line_len), parameter :: gks_eigenvalues_name     = "gks_eigenvalues"
      character(line_len), parameter :: band_structure_name      = "band_structure"
      character(line_len), parameter :: els_potential_name       = "els_potential"
      character(line_len), parameter :: mixer_report_name        = "mixer_report"
      character(line_len), parameter :: solver_report_name       = "solver_report"

      character(line_len) :: p_error_path
      character(line_len) :: f_error_path
      character(line_len) :: diary_path
      character(line_len) :: dcomp_path
      character(line_len) :: new_crystal_path
      character(line_len) :: new_kpoints_path
      character(line_len) :: new_velocity_path
      character(line_len) :: new_lattice_group_path
      character(line_len) :: new_point_group_path
      character(line_len) :: new_space_group_path
      character(line_len) :: new_restart_path
      character(line_len) :: md_trajectory_path
      character(line_len) :: tddft_md_trajectory_path
      character(line_len) :: eigenvalues_path
      character(line_len) :: gks_eigenvalues_path
      character(line_len) :: band_structure_path
      character(line_len) :: els_potential_path
      character(line_len) :: mixer_report_path
      character(line_len) :: solver_report_path

!doc$
      public :: first_unit
      public :: last_unit
      public :: arg_name
      public :: stop_name
      public :: ncp_path
      public :: paw_path
      public :: crystal_path
      public :: kpoints_path
      public :: dsites_path
      public :: velocity_path
      public :: lattice_group_path
      public :: restart_path
      public :: neb_crystal_00_path
      public :: neb_crystal_nn_path
      public :: bulk_crystal_path
      public :: bulk_restart_path
      public :: p_error_path
      public :: f_error_path
      public :: diary_path
      public :: dcomp_path
      public :: new_crystal_path
      public :: new_kpoints_path
      public :: new_velocity_path
      public :: new_lattice_group_path
      public :: new_point_group_path
      public :: new_space_group_path
      public :: new_restart_path
      public :: md_trajectory_path
      public :: tddft_md_trajectory_path
      public :: eigenvalues_path
      public :: gks_eigenvalues_path
      public :: band_structure_path
      public :: els_potential_path
      public :: mixer_report_path
      public :: solver_report_path
      public :: set_paths

!cod$

      contains

      subroutine set_paths(nc,myc,c_myp)
!doc$ subroutine set_paths(nc,myc,c_myp)
!       requires: 1 >= nc < 99. Directories along paths exist in the file system.
!       effects: Defines paths to files.
        integer :: nc, myc, c_myp

!cod$
        character(3) :: config_dir
        character(line_len) :: input_dir, output_dir

        if (nc == 1) then
          input_dir = trim(input_directory)
          output_dir = trim(output_directory)
        else
          write(config_dir,'(i2.2,"/")') myc
          input_dir = trim(input_directory)//trim(config_dir)
          output_dir = trim(output_directory)//trim(config_dir)
        end if

        ncp_path           = trim(potential_directory)//trim(ncp_name)
        paw_path           = trim(potential_directory)//trim(paw_name)

        crystal_path       = trim(input_dir)//trim(crystal_name)
        kpoints_path       = trim(input_dir)//trim(kpoints_name)
        dsites_path        = trim(input_dir)//trim(dsites_name)
        velocity_path      = trim(input_dir)//trim(velocity_name)
        lattice_group_path = trim(input_dir)//trim(lattice_group_name)
        restart_path       = trim(input_dir)//trim(restart_name)
        bulk_crystal_path  = trim(input_dir)//trim(bulk_crystal_name)
        bulk_restart_path  = trim(input_dir)//trim(bulk_restart_name)

        write(p_error_path,'(a,i4.4)') trim(output_dir)//trim(p_error_name), c_myp
        f_error_path             = trim(output_dir)//trim(f_error_name)
        diary_path               = trim(output_dir)//trim(diary_name)
        dcomp_path               = trim(output_dir)//trim(dcomp_name)
        new_crystal_path         = trim(output_dir)//trim(new_crystal_name)
        new_kpoints_path         = trim(output_dir)//trim(new_kpoints_name)
        new_velocity_path        = trim(output_dir)//trim(new_velocity_name)
        new_lattice_group_path   = trim(output_dir)//trim(new_lattice_group_name)
        new_point_group_path     = trim(output_dir)//trim(new_point_group_name)
        new_space_group_path     = trim(output_dir)//trim(new_space_group_name)
        new_restart_path         = trim(output_dir)//trim(new_restart_name)
        md_trajectory_path       = trim(output_dir)//trim(md_trajectory_name)
        tddft_md_trajectory_path = trim(output_dir)//trim(tddft_md_trajectory_name)
        eigenvalues_path         = trim(output_dir)//trim(eigenvalues_name)
        gks_eigenvalues_path     = trim(output_dir)//trim(gks_eigenvalues_name)
        els_potential_path       = trim(output_dir)//trim(els_potential_name)
        band_structure_path      = trim(output_dir)//trim(band_structure_name)
        mixer_report_path        = trim(output_dir)//trim(mixer_report_name)
        solver_report_path       = trim(output_dir)//trim(solver_report_name)

        if (nc == 1) then
          write(config_dir,'(i2.2,"/")') 0
          neb_crystal_00_path = trim(input_directory)//trim(config_dir)//trim(crystal_name)
          neb_crystal_nn_path = trim(input_directory)//trim(config_dir)//trim(crystal_name)
        else
          write(config_dir,'(i2.2,"/")') 0
          neb_crystal_00_path = trim(input_directory)//trim(config_dir)//trim(crystal_name)
          write(config_dir,'(i2.2,"/")') (nc + 1)
          neb_crystal_nn_path = trim(input_directory)//trim(config_dir)//trim(crystal_name)
        end if

      end subroutine

      end module

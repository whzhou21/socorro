README: - describes the executables that can be compiled.
        - describes what is needed to compile.
        - describes conventions.
        - describes the input and output files.
        - describes how to run a test calculation.
        - descriptions the tags used to control a run.

Executables:

  socorro: The main executable. Used to perform electronic structure calculations.
  check_symmetry: Used to determine the symmetry of a set of atom coordinates in a parallelpiped
                  and, optionally, to symmetrize the atom coordinates to machine precision.
  check_kpoints: Used to determine the symmetry of a set of atom coordinates in a parallelpiped
                  and the Monkhorst-Pack special k-points for a set of Monkhorst-Pack parameters.
  taginfo: Used to view the data blocks in a restart file.

What is needed to compile:

  1) Fortan 90/95 compiler
  2) C compiler
  3) FFT library routines
  4) BLAS library
  5) LAPACK library
  6) BLACS library
  7) ScaLAPACK library
  8) Gnu make or equivalent
  9) MPI library
 10) CMMF library (optional)

  Detailed instructions can be found in /socorro/INSTALL.

Conventions:

  Atomic units are used wherein energies are in Ryd, lengths are in Bohr, and e^2 = 2 Ryd*Bohr.
  Calculations are always performed using periodic boundary conditions and Kohn-Sham functions
  (wave functions) are represented either on a real-space mesh or a reciprocal-space mesh. One
  representation can be transformed into the other via a Fourier transform and a wave function
  is normalized such that:

  sum(phi_i*conjg(phi_i)) = N, where the sum is over the N real-space mesh points in a periodic cell

  sum(phi_i*conjg(phi_i)) = 1, where the sum is over the reciprocal-space mesh points

  With these conventions, the real-space to reciprocal-space Fourier transform is scaled by 1/N.

Input and output files.

   Socorro requires a set of input files and a unix directory system, and produces a set of output
   files. The numbers of input and output files depend on the type of calculation being performed.
   The names of the input and output files can be customized by editing /socorro/src/path_mod.f90.
   With two exceptions, the names of the input and output directories can likewise be customized by
   editing /socorro/src/path_mod.f90. The exceptions are the control and interrupt files which must
   be located in the run directory. Brief descriptions of the default input and output files for a
   single configuration run are given below. For calculations with multiple configurations (e.g. dimer
   method), a configuration number (01, 02,...) is inserted after the input or output directory.
   
   Input files:

     argvf: Control file containing parameters used to specify and control a run. The parameters are
            denoted with the tags described below.

     stopf (optional): Interrupt file containing a keyword used to halt a run:
                         "abort" causes the run to stop quickly using the error handling system.
                         "stop" cause the run to stop gracefully, performing post-processing.

     data/crystal: Lattice vectors, atom types, and atom positions. The format is:
                   <identifier_string>
                   <lattice constant>
                   <lattice vector 1>
                   <lattice vector 2>
                   <lattice vector 3>
                   one of ["cartesian", "lattice", or specific variations of these]
                   <number_of_atoms>
                   ATOM_TAG <position in cartesian or lattice coordinates>
                   ATOM_TAG <position in cartesian or lattice coordinates>
                   .....
                   (EOF)

     data/NCP.ATOM_TAG: Data for the ATOM_TAG norm-conserving pseudopotential.

     data/PAW.ATOM_TAG: Data for the ATOM_TAG projector-augmented-wave functions.

     data/kpoints: (optional) Sampling points in the Brillouin zone. The format is:
                   <number of sampling points>
                   one of ["cartesian", "lattice", or specific variations of these]
                   <point in cartesian or lattice coordinates> <degeneracy>
                   .....
                   (EOF)

     data/dsites:  (optional) Sites and radii for decomposition of the Kohn-Sham functions. The format is:
                   <number of sites>
                   <site in lattice coordinates> <radius>
                   <site in lattice coordinates> <radius>
                   .....
                   (EOF)

     data/initial_velocity: (optional) Initial velocities for an MD run. The format is:
                            <velocity components for atom 1>
                            <velocity components for atom 2>
                            .....
                            (EOF)

     data/lgroup: (optional) Lattice point-group operations. The format is:
                  <number of point operators> <number of translations per point operation>
                  <space>
                  <3x3 matrix denoting a point-group operator in lattice representation>
                  .....
                  (EOF)

     data/restartf: (optional) Restart file.

     data/bulk_crystal: (optional) Bulk crystal file used in the LMCC method.

     data/bulk_restartf: (optional) Bulk restart file used in the LMCC method.

   Output files:

     errorf_xxx: Error and warning messages for process xxx.

     diaryf: Diary of the run.

     dcompf: (optional) Decomposition of the Kohn-Sham functions.

     new_crystal: The sequence of new crystal files during relaxation or molecular-dynamics (MD) runs.

     new velocity: New atom velocities after an MD update.

     new_kpoints: (optional) Brillouin zone sampling points generated using the Monkhorst-Pack scheme.

     new_lgroup: (optional) Lattice point-group operations.

     new_pgroup: (optional) Point-group operations of the space group.

     new sgroup: (optional) Space-group operations.

     new_restartf: (optional) Restart file.

     sym_crystal: (optional) Symmetrized crystal file output by check_symmetry.

     md_trajectory: Atom positions during an MD run.

     eigenvalues: Eigenvalues and occupations

     solver_report: (optional) Solver parameters and residual tolerances.

     mixer_report (optional) Mixer parameters, coefficients, and residual tolerances.

   check_symmetry requires two input files and a unix directory system, and produces a set of output
   files with the number of output files depending on the type of calculation being performed. The
   names of the input and output files can be customized by editing /socorro/src/path_mod.f90. With
   one exception, the names of the input and output directories can likewise be customized by editing
   /socorro/src/path_mod.f90. The exceptions is the control file which must be located in the run
   directory

How to run a test calculation:

   To run the test in /socorro/testdata/si_ncp:
     # cd run
     # ./linktest si_ncp
     # ./socorro

Tips for dealing with symmetry problems.

    When there are not as many point-group symmetry operations as expected, check that the magnitudes of
    lattice vectors related by symmetry are equal to each other. If this is not possible, increase the
    number of digits in the lattice vectors to make the them equal to within 1.0e-8 Bohr (the tolerance
    used in symmetry_mod.f90).

Tips for dealing with electronic-structure convergence problems.

    When performing calculations using large supercells, it is not uncommon to encounter slow convergence of
    the electronic structure. When this occurs, try increasing the number of bands (nbands tag), increasing
    the number of solver directions (solver_directions tag), and using kerker weighting in the mixing algorithm
    (mix_weight_tp  tag). You may also try reducing the value of mix_weight_pf and lowering mix_weight_g to
    damp out oscillations in the electron density from one iteration to the next.

Tips for dealing with memory problems.

    When performing calculations for large numbers of atoms and/or large cutoffs, the memory requirements may
    exceed the amount of RAM resulting in poor performance due to swapping to disk or in a crash. The remedy
    for this is to reduce the size of the calculation by lowering the energy cutoff or increase the memory by
    increasing the number of processors. If the problem occurs in the decomposition of the Kohn-Sham functions,
    you can decrease the memory requirement by setting the dcomp_memory tag to s.

Tips when performing transition-state calculations.

    When performing a nudged-elastic-band calculation, it is necessary to null the residual atomic forces (the
    default) and it is necessary that the prior calculations to obtain the atomic structures at the end points
    also were performed with the residual forces nulled.

    When performing a dimer calculation, it is necessary to null the residual atomic forces (the default).

Schematic diagram of program flow (when relaxing atom positions):

     ______________________________ Loop to relax atom positions:
    |                                 Atom positions are updated according to relax_... tags.
    |
    |       _______________________ Loop to attain self-consistency for a configuration of atoms:
    |      |                          New density mixed with earlier densities to form new hamiltonian.
    |      |                          Mixing is controlled via mix_... tags.
    |      |
    |      |      _________________ Loop to update wavefunctions wrt. the current hamiltonian:
    |      |     |                    Wavefunctions are updated according to solver_... tags.
    |      |     |                    New density is constructed from the final wavefunctions.
    |      |     |
    |      |     |
    |      |     |
    |      |     |
    |      |     |
    |      |     |_________________
    |      |
    |      |
    |      |
    |      |_______________________
    |
    |
    |______________________________


Tags and values used to control a run:

   The lines below explain tags and their associated values which can be used to control a Socorro calculation.  To
   use a tag, simply give it on a separate line in argvf followed by the tag value or values. For example, to specify
   a wavefunction cutoff of 20 Ryd add the following line to argvf:

   wf_cutoff  20.0d0

   Note that: Tags can be written with any combination of upper and lower case.
              Logical parameters can be written with any combination of upper and lower case.
              Character values can be written with any combination of upper and lower case.
              Energy values must be given in Rydberg units.
              Distances must be given in atomic units (Bohr radii).

   Tag descriptions are grouped according to code areas, and the format for the descriptions is as follows:

    tag: Meaning
        type; valid values or range of numerical values; default value
        modules where tag is sought; object where value is stored; name of variable holding the value
        first option if tag is not found, second option if tag is not found, ...
        Note #1
        Note #2
        .
        .
        .

Error handling related tags:

    error_file_mode: Determines how many error files to write.
        character; all, first, none; first
        sys_mod, tag ; error_mod, file_status, file_mode
        run aborts
        notes; all (the default) causes all processes to write error files (errorf_xxxx)
               first causes each rank 0 processes in an MPI CONFIG scope write an errof file (errorf)
               none causes no errorf files to be written

Electrons related tags:

    wf_cutoff: Determines the number of plane waves used to expand the Kohn-Sham functions.
        real number; >= 0.0; no default
        multibasis_mod, layout_mod (optional); protobasis_obj; pb%cutoff
        run aborts
        Expansion includes  plane waves for which (G + k)^2 < wf_cutoff, where G is a wave vector and k is a sampling point

    wf_init: Determines how the Kohn-Sham functions are initialized.
        character; RANDOM, DIAGNOSTIC; RANDOM
        multivector_mod; NA; NA
        default is used
        RANDOM: Fourier coefficients are set to random numbers between -0.5 and +0.5
        DIAGNOSTIC: Fourier coefficients are set to Ross Lippert's diagnostic values

    nbands: Number of bands (wavefunctions at a Brillouin zone sampling point).
        integer; >= total charge/2; no default
        multibasis_mod; protobasis_obj; pb%nbands
        run aborts

    write_band_occupations: Controls output of band occupations to diary file.
        character; on, off; off
        electrons_sc_mod; NA; NA
        default used

    write_eigenvalues: Controls output of eigenvalues and occupations to the eigenvalues file.
        character; on, off; off
        electrons_sc_mod; NA; NA
        default used

    kpoints: Method of obtaining sampling points (k-points) in the Brillouin zone.
        character; GMP, GBP, USP; GMP
        kpoints_mod; kpoints_obj; kp%o%mode
        default is used
        For GMP (Generate Monkhorst Pack), a mesh is generated from mpparams (see below). These points are closed using the
           lattice group and then reduced using the point group of the space group (augmented with inversion if it is not
           already present in the space group).
        For GBP (Generate Baldereschi Point), the point at (0.25,0.25,0.25) is used for sampling. This mode  is meant to be
           used only in molecular dynamics calculations and only with a cubic cell and the C_1 space group.
        For USP (Read Special Points), sampling points are read from file /data/kpoints. The format is
           ------------------------------------------
           n
           rep
           k1(1) k1(2) k1(3) d(1)
           k2(1) k2(2) k2(3) d(2)
           .
           .
           .
           kn(1) kn(2) kn(3) d(n)
           ------------------------------------------
           where n is the number of sampling points, rep is the representation (lattice or cartesian), ki(j) is the j'th
           coordinate of the i'th sampling point, and d(i) is the degeneracy. The points are assumed to consistent with the
           lattice group and space group.
        
    kgroups: The number of kgroups (groups of k-points) to be run independently
        integer; >= 1; 1
        mpi_mod; the_mpi; nkgroups
        default is used
        kgroups must divide the number of sgroup processes equally
        kgroups must be >= the number of k-points
        For efficiency, it is recommended that kgroups divide the number of k-points equally, although this is not required

    mpparams: Numbers of Brillouin zone sampling points along the three reciprocal lattice vectors.
        3 integers; -inf < i < +inf; no defaults
        kpoints_mod; kpoints_obj; kp%o%mpp
        run aborts
        Monkhorst-Pack parameters with shifts of the mesh are accomplished by giving one or more negative parameters.
        Values along two reciprocal lattice vectors related by a point group operation must be the same.
        The point at zero for a reciprocal lattice vector is included if the value for that vector is < 0.
        Values 0 and -1 yield the same sampling points as the value 1.
        For a hexagonal lattice, sampling points fall on the line through the origin and parallel to the c-axis.

    mp_closure: Turns on/off closure of the Monkhorst-Pack sampling points using the lattice group.
        character; ON, OFF; OFF
        kpoints_mod
        default is used
        Note: ON is left for backward compatibility but should not be used otherwise

    mp_reduction: Turns on/off reduction of the Monkhorst-Pack sampling points using the double group.
        character; ON, OFF; ON
        kpoints_mod
        default is used

    save_kpoints; .true., .false.; .false.
        kpoints_mod; NA; NA
        default is used
        The file is saved in the run directory with the name "new_kpoints". The format is as given above with rep = "lattice".
        Saves are performed only when the kpoints_mode is GMP (the default setting).

    occupation_method: Determines the method used to occupy the Kohn-Sham states.
        character; THERMAL, UNIFORM; THERMAL
        electrons_xx_mod; electrons_xx_obj; el%o%occupation_method
        default is used
        THERMAL: Fermi function is used to determine band occupations with kt defining the energy broadening
        UNIFORM: k-points are occupied uniformly within each band

    kt: Energy broadening in the THERMAL occupation_method.
        real number; > 0.0; 5.0e-3 (Ryd)
        electrons_xx_mod; electrons_xx_obj; el%o%kt
        default is used

    free_energy: Determines whether or not the Mermin Free Energy expression is used.
        character; ON, OFF; OFF
        electrons_xx_mod; electrons_xx_obj; el%o%use_free_energy
        default is used
        Note: ON is only valid with THERMAL occupation_method

    sgroups: The number of sgroups (spin groups) to be run independently
        integer; 1, 2; 1
        mpi_mod; the_mpi; nsgroups
        default is used
        sgroups must divide the number of configuration processes equally

    polarization_method: Method for determining the spin polarization.
        character; VARIABLE, FIXED; VARIABLE
        electrons_xx_mod; electrons_xx_obj; el%o%polarization_method
        default is used
        Note: VARIABLE determines the spin polarization as part of the electrons solution
        Note: FIXED produces a solution with a fixed (input) spin polarization

    spin_polarization: Spin polarization in the FIXED polarization method
        real number; < half the total number of electrons; 0.0
        electrons_xx_mod; electrons_xx_obj; el%o%spin_polarization
        default is used

    starting_spin_polarization: Starting spin polarization in the VARIABLE polarization method
        real number; < half the total number of electrons; 0.0
        electrons_xx_mod; electrons_xx_obj; el%o%spin_polarization
        default is used
        Note: The spin polarization is fixed at this value during the first solution for the occupations
              and then allowed to vary in subsequent solutions

    decomposition: Switch determining whether or not to decompose the Kohn-Sham functions into s, p, & d spherical harmonics.
        character; ON, OFF; OFF
        config_xx_mod; NA; NA
        default is used
        ON: Decomposes the Kohn-Sham functions into s, p, & d spherical harmonics.
        OFF: Does not decompose the Kohn-Sham functions into s, p, & d spherical harmonics.
        Note: Decomposition is performed for atom sites and for user-supplied sites in file dsites (format given above).
        Note: Results are written to file dcompf (see above).

    dcomp_mode: Type of spherical-harmonic decomposition to perform when decomposition = ON.
        character; L, LM, XYZ; L
        config_xx_mod; NA; NA
        default is used
        L: Reports the l (0, 1, 2) decomposition of the Kohn-Sham functions within the specified radius.
        LM: Reports the l (0, 1, 2) and m decomposition of the Kohn-Sham functions within the specified radius.
        XYZ: Reports the l (0, 1, 2) and x, y, z decomposition of the Kohn-Sham functions within the specified radius.

    dcomp_range: Energy range for decomposition of the Kohn-Sham functions.
        real number; >= 0; full range of Kohn-Sham function eigenvalues
        electrons_xx_mod; NA; NA
        default is used
        Note: Energy range is in Ryd and is distributed half above and half below the Fermi energy.

    dcomp_memory: Memory size for decomposition of the Kohn-Sham functions.
        character; S, M, L; M
        multivector_mod; NA; NA
        default is used
        S: Small memory (rank 2 array) is used.
        M: Medium memory (rank 3 array) is used.
        L: Large memory (rank 4 array) is used.
        Note: Use M or S when the number of atoms > 250.

    dcomp_radius_xxx: Radius within which to perform a spherical-harmonic decomposition for  atom species xxx.
        real number; > 0; 3.0d0
        config_xx_mod; NA; NA
        default is used

Mesh related tags:

    den_cutoff: Determines the number of plane waves in the fields (e.g. electron density) expansions.
        real number; >= 0.0; no default
        layout_mod; layout_obj; lay%o%cutoff
        run aborts
        Expansion includes plane waves for which G^2 < den_cutoff where G is a wave vector

Iterative solver related tags:

    config_type: Type of electronic structure calculation.
        character; SELF-CONSISTENT, FIXED-HAMILTONIAN; SELF-CONSISTENT
        socorro; NA; NA
        default value used
        SELF-CONSISTENT invokes a calculation wherein the wavefunctions and the electron density are iterated to consistency
        FIXED-HAMILTONIAN invokes a calculation wherein the wavefunctions are determined for a fixed hamiltonian

    config_convergence: Quantity used to determine when the configuration (electronic structure) is fully constructed.
        character; NONE, ENERGY, DENSITY, WAVEFUNCTIONS; DENSITY (config_type = SELF-CONSISTENT)
        character; NONE, WAVEFUNCTION; WAVEFUNCTIONS; WAVEFUNCTIONS (config_type = FIXED-HAMILTONIAN)
        config_xx_mod; config_xx_obj; cfg%o%cvg_mode
        default value used
        NONE causes config_steps iterations to be taken in the electronic-structure loop
        ENERGY stops iterating when the cell energy is converged to the energy_tolerance tag value
        DENSITY stops iterating when the electron density residual is converged to the density_tolerance tag value
        WAVEFUNCTIONS stops iterating when the wavefunction residuals are converged to the wavefunctions_tolerance tag value

    config_steps: Maximum number of iterations used to converge the electronic structure.
        integer; >= 0; 40
        config_xx_mod; config_xx_obj; cfg%o%max_steps
        default value used

    energy_tolerance: Tolerance for determining convergence when cfg%o%cvg_mode = ENERGY.
        real number; > 0.0; 1.0d-6
        config_sc_mod; config_sc_obj; cfg%o%energy_tol
        default value used
        The iterative solver is terminated when the cell energy is converged to cfg%o%energy_tol

    density_tolerance: Tolerance for determining convergence when cfg%o%cvg_mode = DENSITY.
        real number; > 0.0; 1.0d-8
        fields_sc_mod; fields_sc_obj; fd%o%res_norm_tol
        default value used
        The iterative solver is terminated when the electron-density residual is converged to cfg%o%density_tol

    wavefunctions_tolerance: Tolerance for determining convergence when cfg%o%cvg_mode = WAVEFUNCTIONS.
        real number; > 0.0; 1.0d-4 (xx = sc) or 1.0d-7 (xx = fh)
        electrons_xx_mod; electrons_xx_obj; el%o%res_norm_tol
        default value used
        The iterative solver is terminated when the wavefunction residuals are converged to el%o%res_norm_tol

Charged-defect related tags:

    charge_state: Charge deviation from a neutral system.
        real number; total charge in the system > 0; 0
        atomic_density_xxx_mod, electrons_xx_mod, and fields_xx_mod;
          atomic_density_xxx_obj, electrons_xx_obj and fields_xx_obj;
          ad%o%charge_state, es%o%charge_state, and fd%o%charge_state
        default is used
        Note: cannot be changed in a restart

    charge_state_mode: mode for inputing the charge state
        character; real_number, integer_ratio; real_number
        atomic_density_xxx_mod, electrons_xx_mod, and fields_xx_mod
        default is used

    charge_state_ratio: ratio (i1/i2) defines charge deviation from a neutral system when charge_state_mode = integer_ratio
        integers; total charge in the system > 0; no defaults
        atomic_density_xxx_mod, electrons_xx_mod, and fields_xx_mod;
          atomic_density_xxx_obj, electrons_xx_obj and fields_xx_obj;
          ad%o%charge_state, es%o%charge_state, and fd%o%charge_state
        run aborts

    compensation: Method used to deal with the divergent electrostatic potential.
        character; UBC, LMCC; UBC
        fields_xx_mod; fields_xx_obj; fd%o%ccompensation
        default is used
        UBC invokes the use of the uniform compensating-background-charge method.
        LMCC invokes the local-moment-counter-charge method. [Physical Review B 60, 1551 (1999)]
                                                             [Physical Review Letters 84, 1942 (2000)]
        Note: The lmcc method is currently available only for ncp calculations.
        Note: cannot be changed in a restart

    lmcc_site: Coordinates (lattice representation) of the counter charge in the lmcc method.
        3 real numbers; [0,1) enforced by program; no defaults
        fields_xx_mod; fields_xx_obj; fd%o%lmcc_site
        program terminates
        Note: cannot be changed in a restart

    lmcc_width: Width of the Gaussian counter charge in the lmcc method.
        real number; > 0.5 and small enough so that the entire charge is contained in the cell; 1.50
        fields_xx_mod; fields_xx_obj; fd%o%lmcc_width
        default is used
        Note: cannot be changed in a restart

    lmcc_bulk: Indicates that the cell is to be surrounded by bulk material in the lmcc method.
        character; ON, OFF; ON
        fields_xx_mod; NA; NA
        default is used

Eigensolver related tags:

    solver_report: Switch indicating whether or not to print a solver report.
        character; ON, OFF; OFF
        eigensolver_mod; eigensolver_obj; es%o%solver_report
        default is used
        default file name is solver_report

    solver_rf_mode: Switch indicating whether or not to use special solver parameters in the first wavefunctions solve.
        character; ON, OFF; ON
        eigensolver_mod; eigensolver_obj; es%o%rf_mode
        default is used
        This capability is only used in a self-consistent restart run with either "ef" or "f" mode.

    solver_method: Type of eigensolver to use.
        character; CONJUGATE-GRADIENTS (CG), GRASSMAN-CONJUGATE-GRADIENTS (GCG), BLOCKED-DAVIDSON (BD); BD
        eigensolver_mod; eigensolver_obj; es%o%solver_method
        default is used

    solver_dir_profile: Profile of the step-dependent maximum number of wavefunction updates (directions).
        character; FIXED, STEPPED; FIXED
        eigensolver_mod; eigensolver_obj; es%o%dir_profile
        default is used

    solver_rfm_directions: Maximum number of directions in the first wavefunctions solve for a restart in "f" or "fe" mode.
        integer; > 0; config_type = self-consistent: 15, 15, 15  (es%o%method = CG, GCG, BD)
                    ; config_type = fixed-hamiltonian: not invoked
        eigensolver_mod; eigensolver_obj; es%o%rfm_dir
        default is used
        This capability is only used when rf_mode is on.

    solver_directions: Maximum number of directions in a FIXED profile or in level 1 of a STEPPED profile.
        integer; > 0; config_type = self-consistent: 3, 10, 3  (es%o%method = CG, GCG, BD)
                    ; config_type = fixed-hamiltonian: 5, 12, 5  (es%o%method = CG, GCG, BD)
        eigensolver_mod; eigensolver_obj; es%o%dir
        default is used

    solver_dir_levels: Number of levels when using a STEPPED directions profile.
        integer; 2 or 3
        run halts

    solver_dir_step_1-2: Step at which the first change takes place in a STEPPED directions profile.
        integer; > 1 and < es%o%max_steps
        run halts

    solver_dir_2: Maximum number of directions in level 2 of a STEPPED profile.
        integer; > 0
        run halts

    solver_dir_step_2-3: Step at which the second change takes place in a STEPPED directions profile.
        integer; > solver_dir_step_1-2 and < es%o%max_steps
        run halts

    solver_dir_3: Maximum number of directions in level 3 of a STEPPED profile.
        integer; > 0
        run halts

    solver_tol_profile: Profile of the step-dependent tolerance for wavefunctions updates.
        character; FIXED, STEPPED; FIXED
        eigensolver_mod; eigensolver_obj; es%o%tol_profile
        default is used

    solver_rfm_tolerance: Tolerance in the first wavefunctions solve for a restart in "f" or "fe" mode.
        real number; >= 0.0; config_type = self-consistent: 1.0e-5, 1.0e-8, 1.0e-5 (es%o%method = CG, GCG, BD)
                    ; config_type = fixed-hamiltonian: not invoked
        eigensolver_mod; eigensolver_obj; es%o%rfm_tol
        default is used
        This capability is only used when rf_mode is on.

    solver_tolerance: Tolerance in a FIXED profile or in level 1 of a STEPPED profile.
        real number; >= 0.0; config_type = self-consistent: 1.0e-4, 1.0e-5, 1.0e-4 (es%o%method = CG, GCG, BD)
                          ; config_type = fixed-hamiltonian: 5.0e-6, 1.0e-8, 5.0e-6 (es%o%method = CG, GCG, BD)
        eigensolver_mod; eigensolver_obj; es%o%tol
        default is used
        Tolerance refers to the RMS wavefunction residual

    solver_tol_levels: Number of levels when using a STEPPED tolerance profile.
        integer; 2 or 3
        run halts

    solver_tol_step_1-2: Step at which the first change takes place in a STEPPED tolerance profile.
        integer; > 1 and < es%o%max_steps
        run halts

    solver_tol_2: Tolerance in level 2 of a STEPPED profile.
        real number; >= 0.0
        run halts

    solver_tol_step_2-3: Step at which the second change takes place in a STEPPED tolerance profile.
        integer; > solver_tol_step_1-2 and < es%o%max_steps
        run halts

    solver_tol_3: Tolerance in level 3 of a STEPPED profile.
        real number; >= 0.0
        run halts

    remap_type: Type of routine used to remap wavefunctions for fast-Fourier transforms.
        character; SJP, MPI; SJP
        multibasis_mod; multibasis_obj; mb%o%remap_type
        default is used
        SJP invokes the C remap_2d routine written by Steve J. Plimpton
        MPI invokes the MPI_ALLTOALLV routine

    diagonalization_method: Type of diagonalization routine to use.
        character; LAPACK, SCALAPACK, CMMF; LAPACK
        eigensolver_mod; eigensolver_obj; es%o%diagonalization_method
        default is used
        LAPACK invokes LAPACK diagonalization routines
        SCALAPACK invokes ScaLAPACK diagonalization routines
        CMMF invokes Mark Sears' CMMF diagonalization routines (not available for public release)
        Note: to use the CMMF diagonalization routine, do the following before compiling:
          1. Remove the # in the following line in make.conf: #USE_CMMF = 1
          2. Search for the phrase "Sandia build" in eigensolver_mod.f90 and follow the instructions.
          3. Copy cmmf.c to the src directory.

    blocksize: Blocksize for use in scalapack routines.
        integer; > 0; 30
        eigensolver_mod; NA; NA
        default is used

    memory_factor: Multiplier for work space in scalapack generalized diagonalization routines
        integer; > 0; 100
        eigensolver_mod; NA; NA
        default is used

Projector related tags:

    projector_type: Type of non-local projectors to use.
        character; RECIPROCAL, REAL; RECIPROCAL
        operators_mod; h_common_obj; hc%o%projector_type
        default value used
        RECIPROCAL invokes reciprocal-spaced projectors
        REAL invokes real-space projectors constructed using the scheme proposed by King-Smith et al.

    projector_radius_xxx: Radius used to optimize real-space projectors for atom type xxx.
        real number; < radius of sphere which will fit inside the supercell; 4.0
        ncp_data_mod; ncp_data_obj; pd%o%r_opt
        run aborts

    pdots_blocksize: Blocksize used to reduce memory usage in calculations of pdots and grad pdots.
        integer; >= 1 and <= number-of-bands; number-of-bands
        operators_mod; h_common_obj; hc%o%pbs
        default value used
        Note that this is only implemented with Fourier-space projectors

    optimization_points: Number of points used to optimize the real-space projectors.
        integer; > 0; 171
        ncp_data_mod; NA; NA
        default is used
        RESET THIS VALUE ONLY IF YOU UNDERSTAND THE OPTIMIZATION ROUTINE!

Self-consistency related tags:

    self-consistency_method: Method used to achieve self-consistency
        character; MIXING, OPTIMIZED-EFFECTIVE-POTENTIAL; MIXING (xctype_dependence = density or hybrid)
                                                          OEP (xctype_dependence = orbital)
        fields_sc_mod; fields_sc_obj; fd%o%self_consistency_method
        default is used

    mix_report: Switch indicating whether or not to print a mixer report.
        character; ON, OFF; OFF
        mixer_mod; mixer_obj; mx%o%mixer_report
        default is used
        default file name is mixer_report

    mix_type: Type of mixing.
        character; DENSITY, POTENTIAL; DENSITY
        fields_sc_mod; fields_sc_obj; fd%o%mixing_type
        default is used

    mix_method: Mixing method to be used.
        character; SIMPLE, PULAY; PULAY
        mixer_mod; mixer_obj; fmx%o%method
        default is used
        Note: To turn mixing completely off, use: mix_method = SIMPLE, mix_weight = 1.0, and mix_damping_type = NONE.

    mix_wgt_profile: Profile of the step-dependent weight.
        character; FIXED, RAMPED; FIXED
        mixer_mod; mixer_obj; mx%o%wgt_profile
        default is used

    mix_weight: Mixing weight in a FIXED profile and in the initial phase of a RAMPED profile.
        real number; > 0.0 and <= 1.0; 0.8
        mixer_mod; mixer_obj; fmx%o%wgt
        default is used

    mix_weight_dp: Mixing weight for the dyad potential in a generalized Kohn-Sham calculation.
        real number; > 0.0 and <= 1.0; 0.8
        mixer_mod; mixer_obj; fmx%o%wgt_dp
        default is used

    mix_wgt_final: Mixing weight in the final phase of a RAMPED profile.
        real number; > 0.0 and <= 1.0; mix_weight + 0.2
        default is used

    mix_wgt_ramp_start: Step at which the ramp is started in a RAMPED profile.
        integer: > 1; 8
        default is used

    mix_wgt_ramp_stop: Step at which the ramp is stopped in a RAMPED profile.
        integer: > mix_wgt_ramp_start; mix_wgt_ramp_start + 6
        default is used

    mix_damping_type: Type of damping used with field quantities.
        character; NONE, KERKER; NONE
        mixer_mod; mixer_obj; mx%o%damping_type
        default is used

    mix_dgc_profile: Profile of the step-dependent KERKER damping G crossover.
        character; FIXED, RAMPED; FIXED
        mixer_mod; mixer_obj; mx%o%dgc_profile
        default is used

    mix_damping_gc: G crossover in a FIXED profile and in the initial phase of a RAMPED profile when using KERKER damping.
        real number; > 0.0; 0.8
        mixer_mod; mixer_obj; mx%o%dgc
        default is used
        For non-zero wave vector G, dgc = 1.0/(1.0 + dgc^2/G^2). for G = 0, dgc = 1.0.

    mix_dgc_final: G crossover in the final phase of a RAMPED profile when using KERKER damping.
        real number; > 0.0; mix_damping_gc + 0.2
        default is used

    mix_dgc_ramp_start: Step at which the ramp is started in a RAMPED profile.
        integer: > 1; 8
        default is used

    mix_dgc_ramp_stop: Step at which the ramp is stopped in a RAMPED profile.
        integer: > mix_dgc_ramp_start; mix_dgc_ramp_start + 6
        default is used

    mix_mxh_profile: Profile of the step-dependent PULAY history.
        character; FIXED, STEPPED; FIXED
        mixer_mod; mixer_obj; mx%o%mxh_profile
        default is used

    mix_history: Maximum number of saved mixing residuals (history) in PULAY mixing.
        integer; > 0 and < 21; 5
        mixer_mod; mixer_obj; mx%o%mxh
        default is used

    mix_mxh_levels: Number of levels when using a STEPPED history profile.
        integer; 2 or 3
        run halts

    mix_mxh_step_1-2: Step at which the first change takes place in a STEPPED history profile.
        integer; > mix_history and < mx%o%max_steps
        run halts

    mix_mxh_2: Number of histories in level 2 of a STEPPED history profile.
        integer; > 0 and < 21
        run halts

    mix_mxh_step_2-3: Step at which the second change takes place in a STEPPED history profile.
        integer; > (mix_mxh_step_1-2 + mix_mxh_2) and < mx%o%max_steps
        run halts

    mix_mxh_3: Number of histories in level 3 of a STEPPED history profile.
        integer; > 0 and < 21
        run halts

    mix_metric_type: Type of metric used in PULAY mixing.
        character; UNITY, KERKER; UNITY
        mixer_mod; mixer_obj; mx%o%metric_type
        default is used
        For mix_metric_tp = UNITY, weight is 1.0 independent of wave vector.
        For mix_metric_tp = KERKER, weight depends on wave vector according to mix_metric_gc.

    mix_metric_gc: G crossover used in KERKER metric.
        real number; >= 0.0; 0.8
        mixer_mod; mixer_obj; mx%o%metric_gc
        default is used
        For non-zero wave vector G, metric = 1.0/(1.0 + metric_gc^2/G^2). for G = 0, metric = 1.0.

    oep_optimizer: Type of optimizer used to obtain the OEP
        character; CHEBYSHEV, SIMPLE; CHEBYSHEV
        oep_mod; oep_obj; oep%o%oep_optimizer
        default is used

    oep_chebyshev_norm: Norm used to initialize the Chebyshev optimizer.
        real; > 0.0; 1.0
        oep_mod; oep_obj; oep%o%cs_norm
        default is used

    oep_chebyshev_cond: Condition number used to initialize the Chebyshev optimizer.
        real; > 0.0; 1000.0
        oep_mod; oep_obj; oep%o%cs_cond
        default is used

    oep_chebyshev_period: Reset period used in the Chebyshev optimizerer.
        integer; > 0; 1 + floor(2.0*sqrt(oep%o%cs_cond))
        oep_mod; oep_obj; oep%o%cs_period
        default is used

    oep_step_size: Amount of the (preconditioned) gradient added to the self-consistent potential in the SIMPLE optimizer.
        real; > 0.0; 8.0
        oep_mod; oep_obj; oep%o%step_size
        default is used
        scp --> scp + step_size*(preconditioner*de_dv)

    oep_preconditioner: Type of preconditioner to use with the SIMPLE OEP optimizer
        character; NONE, PC1, PC2; PC2
        oep_mod; oep_obj; oep%o%preconditioner
        default is used
        PC1: de_dv -> de_dv(1  + pc_factor*|G|^pc_power)
        PC2: de_dv -> de_dv(pc_factor*|G|^pc_power)

    oep_pc_factor: Multiplying factor in PC preconditioners
        real; >= 0; 1.0
        oep_mod; oep_obj, oep%o%pc_factor
        default is used

    oep_pc_power: Exponent of |G| in the form of the PC preconditioners
        real; >= 0; 2.0
        oep_mod; oep_obj, oep%o%pc_power
        default is used

    oep_normalization: Type of nNormalization applied to the OEP exchange-correlation potential.
        character; ZERO, DXCP; ZERO (wavefunction-dependent functional), DXCP (density-dependent functional)
        oep_mod; oep_obj; oep%o%normalization
        default is used
        note: DXCP will cause an error with wavefunction-dependent functionals
        note: ZERO will give a warning with density-dependent functionals

    oep_solver_dir: Number of directions used in the conjugate-gradients oep solver.
        integer; > 0; 50
        oep_mod; oep_obj; oep%o%solver_dir
        default is used

    oep_solver_tol: Tolerance used in the conjugate-gradients oep solver.
        real; > 0.0; 1.0e-8
        oep_mod; oep_obj; oep%o%solver_tol
        default is used

    oep_cutoff: Cutoff energy used to filter the OEP exchange-correlation potential.
        real; > 0.0; wf_cutoff
        oep_mod; oep_obj; oep%o%cutoff
        default is used

Crystal related tags:

    neighbor_range: Determines the range used to construct a neighbor list.
        real; > 0; neighbor_range_default = 5.0
        crystal_mod; crystal_obj; cr%o%neighbor_range
        default is used

Symmetry related tags:

    lattice_symmetry: Determines how the lattice group is obtained.
        character; AUTO, USER; AUTO
        symmetry_mod; point_group_obj; pg%o%mode (AUTO or USER)
        default is used
        For AUTO, the lattice group is generated from the lattice.
        For USER, the lattice group is read from file lgroup with format given below.
        The lattice group is used to close the Monkhorst-Pack k-point mesh.

    symmetry: Determines whether or not space-group symmetry is used.
        character; FULL, AUTO, OFF; AUTO
        symmetry_mod; space_group_obj; sg%o%mode
        default is used
        For FULL, the space group is generated from the lattice group and atom positions. A new space group is generated
            when the lattice group or atom positions change.
        For AUTO, the space group is generated from the lattice group and atom positions. A new space group is generated
            when the lattice group changes and when the old space group is not compatible with new atom positions.
        For OFF, the trivial space group (C_1) is generated.
        The space group (augmented by inversion if it is not already present) is used to reduce the closed Monkhorst-Pack points.

    sg_tolerance: Tolerance used to determine valid atom permutations.
        real; > sg_tolerance/minimum-lattice-vector_magnitude > 1.0e-4; tolerance/minimum-lattice-vector-magnitude = 1.0e-8
        symmetry_mod; space_group_obj; sg%o%tol_aperm
        default is used
        Note: If a run stops with the message "ERROR: unknown point-group type" in the error files, this indicates that the
              atom positions are slightly off of symmetry positions. The user should either add more digits to the atom
              positions in the crystal file or loosen sg_tolerance. Likewise if the space-group symmetry is not as large as
              expected.

    symmetrize_atoms: Symmetrizes the starting atom positions.
        character; .true., .false.; .false.
        check_symmetry_mod, external_mod; NA; NA
        default is used
        
    list_lattice_group: Prints the lattice-group operations to diaryf.
        logical; .TRUE., .FALSE.; .FALSE.
        symmetry_mod; NA; NA
        default is used

    list_space_group: Prints the space-group operations to diaryf.
        logical; .TRUE., .FALSE.; .FALSE.
        symmetry_mod; NA; NA
        default is used
        Note: The matrices are printed in the lattice representation. To compare these with the matrices given in
              Burns and Glazer, the printed matrices must be transformed to the cartesian representation.

    save_lattice_group; .true., .false.; .false.
        symmetry_mod; NA; NA
        default is used
        The file is saved in the run directory with the name "new_lgroup". The format is as follows:
           number_of_point_operators number_of_translations_per_point_operation
           <space>
           3x3 matrix denoting the first point operator (in lattice representation)
           <space>
           3x3 matrix denoting the next point operator (in lattice representation)
           <space>
           .
           .
    save_space_group; .true., .false.; .false.
        symmetry_mod; NA; NA
        default is used
        The file is saved in the run directory with the name "new_sgroup". The format is as follows:
           # of point operators # of translations per point operation
           <space>
           3x3 matrix denoting the first point operator (in lattice representation)
           first translation for the first point operator
           second translation for the first point operator
           .
           .
           last translation for the first point operator
           <space>
           3x3 matrix denoting the next point operator (in lattice representation)
           first translation for the next point operator
           second translation for the next point operator
           .
           .
           last translation for the next point operator
           <space>
           .
           .

Exchange-Correlation-Type related tags:
  - specifies the exchange-correlation functional
  - prefix is xctype_

    xctype_dependence: Dependence of the functional.
         character; density, orbital, hybrid; density
         xc_type_mod; xc_type_obj; xct%o%dependence
         default is used (if no restart file is found)
         Note: density: density-dependent (DD, semilocal) functional
               orbital: exact exchange (EXX) functional
               hybrid: mixed EXX and DD functional

    xctype_ddf_source: Source of the density-dependent functionals
         character; native, libxc; native
         xc_type_mod; xc_type_obj; xct%o%ddf_source
         default is used (if no restart file is found)
         Note: Native means internal Socorro routines.
               Native can be used with density-dependent functionals.
               Libxc can be used with hybrid and density-dependent functionals.
         Note: There are two categories of hybrids:
               1. Named hybrids in which the EXX fraction is fixed (and provided by libxc).
               2. Custom hybrids in which the EXX fraction is set by the user.

    xctype_functional: Type of exchange-correlation functional.
        character; (native) lda, lsda, pw91, pbe, blyp, am05;
        integer; (libxc) see the file libxc_funcs.f90 in the libxc-2.x.x src directory (types having _XC_ in the name)
        xc_type_mod; xc_type_obj; xc%o%functional
        checks in a restart file and then for the legacy tag (functional)
        Note: (native) Used to designate an exchange plus its default correlation (LDA + PZ).
        Note: (libxc)  Used to designate a combined exchange-correlation functional.
        Note: Named hybrid functionals are designated with xctype_functional.

    xctype_exchange: Type of exchange functional.
        character; (native) lda, lsda, pw91, pbe, blyp, am05;
        integer; (libxc) see the file libxc_funcs.f90 in the libxc-2.x.x src directory (types having _X_ in the name)
        xc_type_mod; xc_type_obj; xc%o%exchange
        checks in a restart file and then for the legacy tag (exchange)
        Note: Only checked if xctype_functional is not present.
        Note: Custom hybrid functionals are designated with xctype_exchange and xctype_correlation.

    xctype_correlation: Type of correlation functional.
        character; (native) PZ, PW, LYP;
        integer; (libxc) see the file libxc_funcs.f90 in the libxc-2.x.x src directory (types having _C_ in the name)
        xc_type_mod; xc_type_obj; xc%o%correlation
        checks in a restart file and then for the legacy tag (correlation)
        Note: Only checked if xctype_functional is not present.

    xctype_derivative_method: Method used to calculate (functional) derivatives of density-dependent functionals.
        character; analytical, numerical; depends on functional type
        xc_type_mod; xc_type_obj; xc%o%derivative_method
        checks in a restart file and then for the legacy tag (derivative_method)

    xctype_hybrid_mixing: exact-exchange fraction of a hybrid functional.
        real; > 0.0 and < 1.0;
        xc_type_mod; xc_type_obj; xc%o%hybrid_mixing
        Named hybrid: xc%o%hybrid_mixing is obtained from libxc (user input ignored)
        Custom hybrid: Run aborts if xctype_hybrid_mixing is not found.

    xctype_omega_orb: Screening length for orbital-dependent part of hybrid functionals
        real; > 0; 0.11
        xc_type_mod; xc_type_obj; xct%o%sd_aux_form
        default is used
        Note: Specifying this will override the default screening length in named hybrids.

    xctype_omega_den: Screening length for density-dependent part of hybrid functionals
        real; > 0; 0.11
        xc_type_mod; xc_type_obj; xct%o%sd_aux_form
        default is used
        Note: Specifying this will override the default screening length in named hybrids.

    xctype_nlcc: Switch to use non-linear core correction with hybrid functionals
        character; on, off; off
        xc_type_mod; xc_type_obj; xct%o%sd_aux_form
        run aborts

Exact -Exchange related tags:
  - provides instructions and parameters for calculating the exchange energy and derivative
  - prefix is exx_

    exx_coulomb_kernel: Type of Coulomb-potential kernel to use in EXX calculations
        character; normal, attenuated, screened; none
        xc_type_mod; xc_type_obj; xc%o%coulomb_kernel
        run aborts

    exx_auxiliary_type: Type of auxiliary function used to treat integrable divergences.
        character; legacy, structure_dependent; legacy
        xc_type_mod; xc_type_obj; xct%o%auxiliary_type
        checks for it in a restart file and uses the default if a restart file is not found

    exx_sd_aux_form: Form of the struture-dependent auxiliary type
        character; sc, bcc, fcc; none
        xc_type_mod; xc_type_obj; xct%o%sd_aux_form
        run aborts

    exx_si_aux_form: Form of the struture-independent auxiliary type
        character; sfh, crg; none
        xc_type_mod; xc_type_obj; xct%o%si_aux_form
        run aborts
        Note: This form is not currently available

    exx_comm_method: Communication method used to compute the energy and derivative
        character; collective, point-to-point, point-to-point_2; collective
        multibasis_mod; multibasis_obj; mb%o%exx_comm_method
        default is used

Atomic-representation related tags:

    atomic_representation: Method for representing atoms.
        character; NCP, PAW; NCP
        atomic_operators_mod; atomic_operators_obj; ao%type
        default is used
        Note: cannot be changed in a restart

    atomic_symmetry: Use of symmetry in computing atomic energy and dij matrices.
        character; ON, OFF; ON
        atomic_operators_mod
        default is used
        Note: Use of atomic symmetry may produce small changes in the dij matrix elements and the atomic energy. These
              changes can be systematically reduced to zero by increasing the number of theta_points and phi_points.

    theta_points_xxx: Number of points in the angular mesh for angle theta and atom xxx.
        integer; theta_points > 0; 12
        atomic_operators_paw_mod
        default is used
        Note: for optimal performance, theta_points should be adjusted for each set of calculations.

    phi_points_xxx: Number of points in the angular mesh for angle phi and atom xxx.
        integer; phi_points > 0; 12
        atomic_operators_paw_mod
        default is used
        Note: for optimal performance, phi_points should be adjusted for each set of calculations.

     atomic_occupations_xxx: Initial occupations of the wij matrix for atom type xxx in non-spin-polarized calculations.
        real(number of basis functions); >= 0; INITOCC values in PAW.xxx
        atomic_operators_paw_mod
        default is used

     atomic_occupations_sgx_yyy: Initial occupations of the wij matrix for spin group x and atom type yyy.
        real(number of basis functions); >= 0; INITOCC values in PAW.xxx
        atomic_operators_paw_mod
        default is used

Post-processing related tags:

    forces: Controls whether or not the forces are computed automatically
        character; ON, OFF, .TRUE., .FALSE.; OFF
        config_sc_mod; NA; NA
        default is used
        Note: Even when the forces tag = OFF, the forces will be computed when x_forces(cfg)
              or diary_forces(cfg) is called.

    null_residual_force: Switch indicating whether or not to null residual forces.
        character; ON, OFF; ON
        config_sc_mod; config_sc_obj; cfg%o%null_residual_force
        default is used
        The residual force is the sum of all atom forces.

    pressure: Controls whether or not the pressure is computed automatically
        character; ON, OFF, .TRUE., .FALSE.; OFF
        config_sc_mod; NA; NA
        default is used
        Note: Even when the pressure tag = OFF, the pressure will be computed when x_pressure(cfg)
              or diary_pressure(cfg) is called.

    stress_tensor: Controls whether or not the stress tensor is computed automatically
        character; ON, OFF, .TRUE., .FALSE.; ON
        config_sc_mod; NA; NA
        default is used
        Note: The stress tensor calculation is not currently implemented for the PAW method.
        Note: Even when the stress_tensor tag = OFF, the stress_tensor will be computed when x_stress_tensor(cfg)
              or diary_stress_tensor(cfg) is called.

    write_els_potential: Determines whether or not to write a sxdefectalign els_potential file
        character; ON, OFF; OFF
        config_sc_mod; NA; NA
        default is used
        ON: els_potential file will be written
        OFF: els_potential file will not be written
        Note: The els_potential file is written at the end of a run.

Restart related tags:

    restart: Determines whether or not to restart and how much information to read
        character; OFF, E, EF, EFE, F, FE; OFF
        config_xx_mod; NA; NA
        default is used
        OFF: restart information will not be read
        E: external (crystal) information will be read
        EF: external and fields (atomic and grid densities) information will be read
        EFE: external, fields and electrons information will be read
        F: fields information will be read
        FE: fields and electrons (Kohn-Sham functions) information will be read
        Note: Information about the exchange-correlation functionals are in external.

    write_restart: Determines whether or not to write a restart file and how much information to write
        character; OFF, E, EF, EFE, F, FE; OFF
        config_sc_mod; NA; NA
        default is used
        OFF: restart information will not be written
        E: external (crystal) information will be written
        EF and F: external and fields (atomic and grid densities) information will be written
        EFE and FE: external, fields and electrons information will be written
        Note: The restart file is written at the end of a run.

Transition-state related tags:
  - provides general parameters related to transition-state-finding methods

    configurations: The number of configurations (electronic structure calculations) to run simultaneously.
        integer; > 0, < 99; 1
        mpi_mod; the_mpi; nconfigs
        default is used
        configurations must divide the number of process equally

    ts_method: Determines which transition-state finding method to use.
        character; DIMER, NEB, NONE; NONE
        transition_state_mod; dimer_obj, neb_obj; NA
        default is used
        DIMER: dimer method will be used
        NEB: nudged-elastic-band method will be used
        NONE: transition-state calculation will not be performed

Dimer-method related tags:
  - provides instructions and parameters related to dimer-method calculations
  - prefix is dimer_

    dimer_force_tol: force tolerance used to determine when the dimer translation is finished
        real; > 0; 0.001
        transition_state_mod; dimer_obj; force_tol
        default is used

    dimer_torque_tol: force tolerance used to determine when the dimer rotation is finished
        real; > 0; 0.0001
        transition_state_mod; dimer_obj; torque_tol
        default is used

    dimer_separation: enforced separation of the dimer configurations
        real; > 0; 0.01
        transition_state_mod; dimer_obj; separation_mag
        default is used
        Note: If the supplied configurations do not have this separation, they are modified to have it.

    dimer_max_steps: maximum number of configuration updates in the dimer calculation
        integer; > 0; 100
        transition_state_mod; dimer_obj; max_steps
        default is used

    dimer_optimizer: optimization method to be used in the dimer calculation
        character; sd, cg; cg
        transition_state_mod; NA; NA
        default is used
        Note: sd = steepest-descent, cg = conjugate-gradient

    dimer_sd_alpha: alpha parameter in a dimer sd calculation
        real; > 0; 3.5
        transition_state_mod; optimize_dimer_sd_i routine
        default is used

    dimer_cg_finite_diff: finite_diff parameter in a dimer cg calculation
        real; > 0; 0.01
        transition_state_mod; optimize_dimer_cg_i routine
        default is used

    dimer_cg_max_move: max_move parameter in a cg dimer calculation
        real; > 0; 6.0
        transition_state_mod; optimize_dimer_cg_i routine
        default is used

Nudged-elastic-method related tags:
  - provides instructions and parameters related to nudged-elastic-method calculations
  - prefix is neb_

    neb_force_tol: force tolerance used to determine when the neb calculation is finished
        real; > 0; 0.001
        transition_state_mod; neb_obj; force_tol
        default is used

    neb_climbing_image: switch indicating whether or not the climbing-image method is used
        logical; .true., .false; .true.
        transition_state_mod; neb_obj; climbing_image
        default is used

    neb_spring: spring constant used in the neb calculation
        real; > 0; 0.01
        transition_state_mod; neb_obj; spring
        default is used

    neb_max_steps: maximum number of configuration updates in the neb calculation
        integer; > 0; 100
        transition_state_mod; neb_obj; max_steps
        default is used

    neb_force_test: indicates how to apply force_tol
        character; ALL_IMAGES, CLIMBING_IMAGE; ALL_IMAGES
        transition_state_mod; neb_obj; force_test
        default is used

    neb_optimizer: optimization method to be used in the neb calculation
        character; SD, QM; QM
        transition_state_mod; NA; NA
        default is used
        Note: SD = steepest-descent, QM = quickmin

    neb_sd_factor: prefactor in steepest-descent optimizer
        real; > 0; 2.50
        transition_state_mod; NA; NA
        default is used

    neb_qm_time_step: time step used in the quickmin optimizer
        real; > 0; 1.5
        transition_state_mod; NA; NA
        default is used

    neb_qm_max_move: maximum move allowed in the quickmin optimizer
        real; > 0; 1.0
        transition_state_mod; NA; NA
        default is used

Structural-optimization related tags:
  - provides instructions and parameters related to structural optimization
  - prefix is relax_

    relax_method
        character; none, steepest_descent, conjugate_gradient, quench_minimization; none
  
    relax_tol
        real number; ; 1.d-3
  
    relax_max_steps
        integer; 0 < relax_max_steps ; 100

    relax_prefactor
        real number; 0 < relax_prefactor; 1.0

    relax_time_step
        real number; 0 < relax_time_step; 1.0

    lattice_relaxation
	character; none, shape, shape_quench, shape_quench_z
 	NONE default does not relax lattice
	SHAPE relaxes lattice by minimizing the stress tensor down to LATTICE_RELAX_TOL with a linear algorithm
	SHAPE_QUENCH performs lattice relaxation with a quenched motion
	SHAPE_QUENCH_Z performs the same lattice relaxation as shape_quench except it only relaxes along a single direction

    lattice_relax_tol
	real number; 0 < lattice_relax_tol; 1.0	
	sets the threshold for the various lattice_relaxation options

    extrapolation_method: method used to extrapolate the smooth density after atom moves (also used in MD calculations)
        character; NONE, GUESS_DENSITY; GUESS_DENSITY
        fields_sc_mod; fields_sc_obj; extrapolation_method
        default is used

ENERGY MINIMIZATION:

Fixed volume optimization of the atomic positions. All atoms are moved
- there is currently no method for constraining specific atoms. 
This is called by the statement "if (optimize_lattice(cfg)) call
diary(cfg)"  in socorro.f90.  
optimize_lattice returns a logical that indicates whether any
modifications to cfg have occurred.  The cfg that is returned
corresponds to the optimized positions.

The relevant options in argvf are as follows

relax_method: Specify the relaxation method to be used.

 NONE default Do not perform a structural optimization

 STEEPEST_DESCENTS, SD  'steepest descents' (I have seen
  different definitions of 'steepest descents'.) What is implemented
  here is the simple approach of changing the coordinates by             
  F*relax_prefactor where F is the current force and relax_prefactor 
  is a parameter that you can set (see below).  This is usually NOT 
  the best way to get to a minimum.  It was included because it was the 
  obvious first thing to code and is a standard 'brute force' approach to
  optimization.

 CONJUGATE_GRADIENT, CG  'conjugate gradients' 
  Implements a conjugate gradient search for the minimum
  energy structure.  (See Press, et.al, 'Numerical Recipes' for 
  a description of this algorithm.)  This is typically the best 
  way to go when the initial positions are close to the minimum.

 QUENCHED_MINIMIZATION, QM 'quenched MD'  
  This implements an algorithm described in Della Valle and Andersen, 
  J. Chem. Phys. 97, 2682 (1992).  It performs a molecular dynamics
  simulation using the 'velocity Verlet' algorithm with the following
  modifications.  At each time step and for each particle, the velocity 
  is reset as follows.  If the projection of the force along the velocity 
  is positive, the velocity is replaced by the projection of the velocity
  along the direction of the force.  If the projection is negative, the   
  velocity is set to zero. This will quench the dynamics to the minimum.   
  This approach seems to be best suited for getting close to the minimum   
  when the initial guess may be poor.

relax_tol:  This is the stopping criteria for all methods.  
  The code stops when the root mean square value of the force components is
  less than this value.  default: 1.d-3.

relax_steps:  This sets a maximum number of force calls that can be made to
  attempt to find the minimum.  Note that in some cases, the actual
  number of calls may exceed this slightly since it will always try
  to finish the line minimizations in the conjugate gradient
  approach.  default: 100

relax_prefactor: This is the constant used in relax_method=SD (see above)
  Note that the efficiency and stability of this method depends on this 
  choice.  If the value is too large, the optimization will become
  unstable.  If it is too small, a large number of force calls is
  required to get to the minimum.  default: 1.0

relax_time_step: This is used in relax_method=QM.  It sets the fixed time
  step for the MD simulation.  default: 1.0 (probably too small for most
  cases.)

atom_mass_xxx This is the mass of the atom in amu used for 
  quench_minimization.  Here xxx is replaced by the tag used in the 
  crystal file.  There should be one line like this for each type.  
  The default is to assume a mass of 1 amu.



Molecular Dynamics related tags:
  - provide instructions and parameters for molecular-dynamics runs
  - pertains to born_oppenheimer_mod, ion_dynamics_mod and ehrenfest_mod
  - prefix is md_

    md_method
         character; none, nve, nvt_rescale, nvt_anderson, nvt_hoover; none
         xc_type_mod; xc_type_obj; xct%o%dependence
         run aborts

    md_time_step
        real number; 0 < md_time_step; 100.

    md_steps
        integer; 0 <= md_steps; 0

    md_skip_steps
        integer; 0 <= md_skip_steps <= md_steps; 0

    md_init_temp
        real number; 0 <= md_init_temp; 0.

    md_desired_temp
        real number; 0 <= md_desired_temp; md_init_temp

    md_temp_freq
        integer; 0 < md_temp_freq; 1

    md_hoover_mass
        real number; 0 < md_hoover_mass; 1000.

    md_gen_velocities
        character; YES, NO; YES

The code can currently perform either a NVE (constant number,
volume and energy) micro-canonical simulation or a NVT (constant
number, volume and temperature) simulation.  For the later case, a
variety of standard thermostat methods are implemented.  The
initial velocities can be either read from the file
"data/initial_velocity" or are set randomly based on an input
temperature.  The integration is performed via the 'velocity
Verlet' algorithm (see any book on MD simulations) with a fixed
time step.  Intermediate atomic positions are output to the file
'md_output' at every time step.  The masses are set with the arg
parameter atom_mass_xxx (see below).

This is called by the statement "if (run_moldyn(cfg)) call
diary(cfg)" in socorro.f90.  run_moldyn returns a logical that
indicates whether any modifications to cfg have occurred.  If an
MD simulation is performed, the cfg returned is that for the final
time step.

The relevant parameters are

md_method: Specify the molecular dynamics method to be used.

 NONE (default) - do no MD.

 NVE perform a NVE simulation using a fixed time-step velocity
  Verlet algorithm.

 NVT_RESCALE perform a NVT (isochoric,isothermal) simulation using a 
  fixed time-step.  The temperature control is through periodic 
  rescaling of the velocities to achieve the desired temperature.

 NVT_ANDERSON perform a NVT (isochoric, isothermal) simulation using 
  a fixed time-step.  The temperature is controlled via a stochastic 
  method due to Anderson (see H C Anderson, J. Chem. Phys. 72, 2384
 (1980)).  At each time step and for each atom, a random velocity 
  from a Maxwell-Boltzman distribution replaces the velocity with a 
  probability give by 1/temp_freq.

 NVT_HOOVER perform a NVT (isochoric, isothermal) simulation using a 
  fixed time-step.  The temperature is controlled using the Hoover
  implementation of the Nosé thermostat. (See, for example, 
  "Understanding Molecular Simulations" by Frenkel and Smit).  The 
  rate of energy flow between the ions and the heat bath is controlled 
  by md_hoover_mass.

md_time_step Time step used for the MD. See note below regarding units.   
  default: 100.

md_steps Number of time steps to integrate the equation of motion. 
  default: 0

md_skip_steps  Number of time steps to ignore before starting to 
  compute averages - in other words (skip_steps)*(time_step) is an 
  equilibration time.  default: 0

md_init_temp  Temperature used to define the initial velocity 
  distributions.  The velocities are selected from a Maxwell-Boltzman
  distribution.  Then they are adjusted to give zero total momentum 
  and rescaled to give the exact temperature requested. default: 0

md_desired_temp Target temperature for the isothermal simulation methods.
  default: md_init_temp

md_temp_freq Parameter that determines the frequency of velocity 
  modifications for the isothermal simulation methods.  For NVT\_RESCALE, 
  the velocity is rescaled every md\_temp\_freq time steps.  For 
  NVT\_ANDERSON, it give the inverse of the probability that an 
  atom will get a random velocity in a given time step.

atom_mass_xxx This is the mass of the atom in amu.  Here xxx is 
  replaced by the tag used in the crystal file.  There should be 
  one line like this for each type.  The default is to assume a 
  mass of 1 amu.

md_hoover_mass Used by NVT\_HOOVER. It is the effective mass 
  associated with the additional coordinate added to control the 
  temperature.  It may need to be adjusted by trial and error.  
  default: 1000.

md_gen_velocities  Specify the method for initializing velocities

 YES Determine the initial velocities based on md_init_temp.
  (default)

 NO Read the initial velocities from the file 'initial_velocity' in the 
  run directory.  This file contains a line with the x, y, and z 
  velocity for each atom on a separate line.  The order of the atoms 
  is assumed to be the same as in the crystal file.


NOTE ON UNITS:

The convention in the code is that energies are in Rydbergs and that 
distances are in Bohrs.  I have made the choice to have the code work 
with the nuclear masses in units of the electron mass - keep in the spirit
of atomic units. For convenience, when masses are entered, they are assumed
to be in amu (atomic mass units) and are converted in the code 
to electron masses.  Having made this choice for the mass unit, the choice 
of the time unit is now fixed.  The unit of time is 3.421E-17 sec.  Since 
a typical MD time step is on the order of a few femtoseconds (fs), 
the typical time steps will be on the order of 100 in the units used here.

There is a subtlety in the determination of the temperature.  In classical 
thermodynamics (MD is classical in the treatment of the ionic motion), 
each degree of freedom has a kinetic energy of kT/2. The issue is related to
the number of degrees of freedom.  If the MD method used conserves 
the total momentum, then the number of degrees of freedom is 3(N-1).  
If the total momentum is not conserved, the number of degrees of freedom 
is 3N.  The NVE and NVT_RESCALE methods conserve the total momentum.
(Actually, for NVT_RESCALE the total momentum remains zero if it starts 
out zero. The initial velocities are generated in the code to have
zero total momentum.)  For these methods, the code uses 3(N-1) degrees of
freedom to determine the temperature.  For NVT\_ANDERSON, the total 
momentum is not conserved, so the code uses 3N degrees of freedom to
compute the temperature.


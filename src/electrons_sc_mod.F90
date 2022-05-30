! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module electrons_sc_mod
!doc$ module electrons_sc_mod

!     One datatype is available here: type(electrons_sc_obj)

!     electrons_sc_mod creates and maintains a set of electrons in calculations with a self-consistent hamiltonian.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use arg_mod
      use diary_mod
      use ghost_mod
      use grid_mod
      use math_mod
      use layout_mod
      use lattice_mod
      use crystal_mod
      use external_mod
      use xc_type_mod
      use kpoints_mod
      use wavefunctions_es_mod
      use operators_mod
      use multivector_mod
      use multibasis_mod
      use gen_density_mod
      use gen_potential_mod
      use dyad_mod
      use dyad_kpoint_mod
      use timing_mod

!cod$
      implicit none
      private

      integer, parameter :: NA = 0

      ! occupation method
      integer, parameter :: THERMAL = 1
      integer, parameter :: UNIFORM = 2

      ! polarization method
      integer, parameter :: FIXED    = 1
      integer, parameter :: VARIABLE = 2

      type, public :: electrons_sc_rep
        integer :: ref
        type(ghost) :: g                                           ! ghost
        type(ghost) :: g_external                                  ! ghost of external object used in construction
        logical :: use_free_energy                                 ! switch for using the free energy
        logical :: res_norm_cvg                                    ! convergence status
        integer :: occupation_method                               ! occupation method
        integer :: polarization_method                             ! spin polarization method
        real(double) :: res_norm                                   ! wavefunctions residual norm
        real(double) :: res_norm_tol                               ! tolerance for wavefunctions residual norm
        real(double) :: total_charge                               ! total electronic charge in the supercell
        real(double) :: charge_state                               ! charge state of supercell
        real(double) :: spin_polarization                          ! spin polarization
        real(double) :: fermi_level                                ! Fermi level
        real(double) :: kt                                         ! kT
        real(double) :: cutoff                                     ! wavefunctions cutoff energy
        real(double) :: kinetic_energy                             ! kinetic energy
        real(double) :: mermin_energy                              ! Mermin energy
        integer, dimension(:), pointer :: kgroup_index             ! mapping of k-points to kgroups
        real(double), dimension(:,:), pointer :: eigs              ! eigenvalues
        real(double), dimension(:,:), pointer :: occs              ! occupations
        real(double), dimension(:,:), pointer :: oep_gksc          ! eigenvalue gks corrections from oep_mod
        type(kpoints_obj) :: kpoints                               ! k-points object
        type(gen_density_obj) :: density                           ! generalized density object
        type(h_common_obj) :: hc                                   ! common hamiltonian object
        type(wavefunctions_es_obj), dimension(:), pointer :: wf    ! wavefunction objects
      end type

      type, public :: electrons_sc_obj
        private
        integer :: ref
        type(electrons_sc_rep), pointer :: o
      end type

!doc$
      public :: electrons_sc
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: wormhole
      public :: x_ref
      public :: x_ghost
      public :: x_fermi_level
      public :: x_energy
      public :: x_kinetic_energy
      public :: x_mermin_energy
      public :: x_wavefunctions
      public :: x_density
      public :: x_kpoints
      public :: x_kt
      public :: x_n_bands
      public :: x_kgroup_index
      public :: x_eigenvalue
      public :: x_eigenvalues
      public :: x_occupation
      public :: x_occupations
      public :: x_residual_norm
      public :: x_converged
      public :: s_oep_gksc
      public :: thermal_occupations
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: diary
      public :: max_energy
      public :: diary_energies
      public :: decompose
      public :: distribute
      public :: release
      public :: write_restart

!cod$
      interface electrons_sc
        module procedure constructor_el
      end interface
      interface update
        module procedure update_el
      end interface
      interface my
        module procedure my_el, my_new_el
      end interface
      interface thy
        module procedure thy_el
      end interface
      interface glean
        module procedure glean_el
      end interface
      interface bequeath
        module procedure bequeath_el
      end interface
      interface assignment(=)
        module procedure assign_el
      end interface
      interface wormhole
        module procedure wormhole_el
      end interface
      interface x_ref
        module procedure el_ref
      end interface
      interface x_ghost
        module procedure el_ghost
      end interface
      interface x_fermi_level
        module procedure el_fermi_level
      end interface
      interface x_energy
        module procedure el_energy
      end interface
      interface x_kinetic_energy
        module procedure el_kinetic_energy
      end interface
      interface x_mermin_energy
        module procedure el_mermin_energy
      end interface
      interface x_wavefunctions
         module procedure el_wavefunctions
      end interface
      interface x_density
        module procedure el_density
      end interface
      interface x_kpoints
         module procedure el_kpoints
      end interface
      interface x_kt
         module procedure el_kt
      end interface
      interface x_n_bands
         module procedure el_n_bands
      end interface
      interface x_kgroup_index
         module procedure el_kgroup_index
      end interface
      interface x_eigenvalue
         module procedure el_eigenvalue
      end interface
      interface x_eigenvalues
         module procedure el_eigenvalues
      end interface
      interface x_occupation
         module procedure el_occupation
      end interface
      interface x_occupations
         module procedure el_occupations
      end interface
      interface x_residual_norm
         module procedure el_residual_norm
      end interface
      interface x_converged
        module procedure el_converged
      end interface
      interface s_oep_gksc
        module procedure el_oep_gksc
      end interface
      interface thermal_occupations
        module procedure thermal_occupations_el
      end interface
      interface forces
        module procedure forces_el
      end interface
      interface pressure
        module procedure pressure_el
      end interface
      interface stress_tensor
        module procedure stress_tensor_el
      end interface
      interface diary
        module procedure diary_el
      end interface
      interface max_energy
        module procedure max_energy_el
      end interface
      interface diary_energies
        module procedure diary_energies_el
      end interface
      interface decompose
        module procedure decompose_el
      end interface
      interface distribute
        module procedure distribute_el
      end interface
      interface release
        module procedure release_el
      end interface
      interface write_restart
        module procedure write_restart_el
      end interface

      contains

! public routines

      function constructor_el(ext,gp,restf) result(el)
!doc$ function electrons_sc(ext,gp,restf) result (el)
        type(external_obj) :: ext
        type(gen_potential_obj) :: gp
        type(tagio_obj), optional :: restf
        type(electrons_sc_obj) :: el
!       requires: Consistent ext and gp.
!       effects: Constructs a new el.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: tag, usage
        integer :: ik, nb, nk, nsg
        integer(long) :: dsize, iosl, ndata, s4
        integer, dimension(2) :: csr
        real(double), dimension(3) :: kpt
        type(layout_obj) :: lay
        type(dyad_kpoint_obj) :: dk

        call start_timer("electrons_sc: constructor")

        if (error("  Error on entry")) then
          el%ref = 0
          allocate( el%o )
          el%o%ref = 0
          goto 999
        end if

        call my(ext)
        call my(gp)
        if (present(restf)) call my(restf)

        el%ref = 0
        allocate( el%o )
        el%o%ref = 0
        el%o%g = x_ghost()

        el%o%g_external = x_ghost(ext)

        call my(x_layout(ext),lay)

        nsg = mpi_nsgroups()

        ! Open the ELECTRONS block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"ELECTRONS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ELECTRONS block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)
        end if

        ! Determine the k-points
        if (present(restf)) then
          call my(kpoints(ext,restf),el%o%kpoints) ; if (error()) goto 200
        else
          call my(kpoints(ext),el%o%kpoints) ; if (error()) goto 200
        end if
        nk = x_n_kpoints(el%o%kpoints)

        ! Open the PARAMETERS block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"PARAMETERS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: PARAMETERS block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)
        end if

        ! Determine the wavefunctions cutoff - amenable to change during a restart
        call arg("wf_cutoff",el%o%cutoff,found)
        if (found) then
          if (error(el%o%cutoff <= 0.0_double,"ERROR: wf_cutoff <= 0")) goto 100
          if (error(el%o%cutoff > 0.5_double*x_cutoff(lay),"ERROR: wf_cutoff is too large")) goto 100
        else
          if (present(restf)) then
            if (i_access(restf)) tios = findfirsttag(restf,"CUTOFF")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: CUTOFF tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_double
              ndata = 1
              call readf(el%o%cutoff,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%cutoff)
          else
            if (error(.true.,"ERROR: wf_cutoff not found")) goto 100
          end if
        end if

        ! Determine the number of bands - amenable to change during a restart
        call arg("nbands",nb,found)
        if (found) then
          if (error(nb <= 0,"ERROR: nbands <= 0")) goto 100
        else
          if (present(restf)) then
            if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_BANDS")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_BANDS tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_long
              ndata = 1
              call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              nb = s4
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,nb)
          else
            if (error(.not.found,"ERROR: nbands was not found")) goto 100
          end if
        end if

        ! Determine the charge state
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"CHARGE_STATE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: CHARGE_STATE tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
            call readf(el%o%charge_state,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%charge_state)
        else
          call arglc("charge_state_mode",tag,found)
          if (.not.found) tag = "real_number"
          select case (trim(tag))
          case ("real_number")
            call arg("charge_state",el%o%charge_state,found)
            if (.not.found) el%o%charge_state = 0.0_double  ! default value
          case ("integer_ratio")
            call arg("charge_state_ratio",csr,found)
            if (error(.not.found,"ERROR: charge_state_ratio was not found")) goto 100
            if (error(csr(2) == 0,"ERROR: denominator = 0")) goto 100
            el%o%charge_state = real(csr(1),double)/real(csr(2),double)
          case default
            if (error(.true.,"ERROR: charge_state_mode was not recognized")) goto 100
          end select
        end if

        ! Compute the total number of electrons
        el%o%total_charge = valence_electrons(ext) - el%o%charge_state
        if (error(el%o%total_charge <= 0.0_double,"ERROR: total charge <= 0")) goto 100

        ! Check the number of bands
        if (error(nb < el%o%total_charge/2.0_double,"ERROR: not enough bands")) goto 100

        ! Determine the occupation method - amenable to change during a restart
        call arglc("occupation_method",tag,found)
        if (found) then
          select case (trim(tag))
          case ("thermal")
            el%o%occupation_method = THERMAL
          case ("uniform")
            el%o%occupation_method = UNIFORM
            if (error(nsg == 2,"ERROR: UNIFORM occupation_method is not compatible with spin")) goto 100
          case default
            if (error(.true.,"ERROR: occupation_method was not recognized")) goto 100
          end select
        else
          el%o%occupation_method = THERMAL  ! default value
          if (present(restf)) then
            if (i_access(restf)) tios = findfirsttag(restf,"OCCUPATION_METHOD")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: OCCUPATION_METHOD tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_long
              ndata = 1
              call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              el%o%occupation_method = s4
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%occupation_method)
          end if
        end if

        select case(el%o%occupation_method)
        case (THERMAL)

          ! Determine kT - amenable to change during a restart
          call arg("kt",el%o%kt,found)
          if (found) then
            if (error(el%o%kt <= 0.0_double,"ERROR: kT <= 0.0")) goto 100
          else
            el%o%kt = 5.0e-3_double  ! default value
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"KT")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR: KT tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_double
                ndata = 1
                call readf(el%o%kt,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%kt)
            end if
          end if

          ! Determine the free energy status - amenable to change during a restart
          call arglc("free_energy",tag,found)
          if (found) then
            select case (trim(tag))
            case ("on","true")
              el%o%use_free_energy = .true.
            case ("off","false")
              el%o%use_free_energy = .false.
            case default
              if (error(.true.,"ERROR: free_energy tag was not recognized")) goto 100
            end select
          else
            el%o%use_free_energy = .false.  ! default value
            if (present(restf)) then
              if (i_access(restf)) tios = findfirsttag(restf,"FREE_ENERGY_STATUS")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR: FREE_ENERGY_STATUS tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_long
                ndata = 1
                call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                select case(s4)
                case (1)
                  el%o%use_free_energy = .true.
                case (0)
                  el%o%use_free_energy = .false.
                end select
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%use_free_energy)
            end if
          end if
          select case (x_functional_dependence(x_xc_type(ext)))
          case (FD_ORBITAL)
            if (error(.not.el%o%use_free_energy,"ERROR: free_energy must be on for an orbital-dependent functional")) goto 100
          case (FD_HYBRID)
            if (error(.not.el%o%use_free_energy,"ERROR: free_energy must be on for a hybrid functional")) goto 100
          end select

          select case (nsg)
          case (1)

            ! Set null values for irrelevant elements
            el%o%polarization_method = NA
            el%o%spin_polarization = 0.0_double

          case (2)

            ! Determine the polarization method - amenable to change during a restart
            call arglc("polarization_method",tag,found)
            if (found) then
              select case (trim(tag))
              case ("fixed")
                el%o%polarization_method = FIXED
              case ("variable")
                el%o%polarization_method = VARIABLE
              case default
                if (error(.true.,"ERROR: polarization_method was not recognized")) goto 100
              end select
            else
              el%o%polarization_method = VARIABLE  ! default value
              if (present(restf)) then
                if (i_access(restf)) tios = findfirsttag(restf,"POLARIZATION_METHOD")
                if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                if (error(tios == TAG_NOT_FOUND,"ERROR: POLARIZATION_METHOD tag was not found")) goto 100
                if (i_access(restf)) then
                  dsize = sizeof_long
                  ndata = 1
                  call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  el%o%polarization_method = s4
                end if
                if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%polarization_method)
              end if
            end if

            ! Determine the fixed/starting spin polarization - amenable to change during restart
            select case (el%o%polarization_method)
            case (FIXED)
              call arg("spin_polarization",el%o%spin_polarization,found)
              if (.not.found) then
                el%o%spin_polarization = 0.0_double  ! default value
                if (present(restf)) then
                  if (i_access(restf)) tios = findfirsttag(restf,"SPIN_POLARIZATION")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (error(tios == TAG_NOT_FOUND,"ERROR: SPIN_POLARIZATION tag was not found")) goto 100
                  if (i_access(restf)) then
                    dsize = sizeof_double
                    ndata = 1
                    call readf(el%o%spin_polarization,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%spin_polarization)
                end if
              end if
            case (VARIABLE)
              call arg("starting_spin_polarization",el%o%spin_polarization,found)
              if (.not.found) then
                el%o%spin_polarization = 0.0_double  ! default value
                if (present(restf)) then
                  if (i_access(restf)) tios = findfirsttag(restf,"SPIN_POLARIZATION")
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
                  if (error(tios == TAG_NOT_FOUND,"ERROR: SPIN_POLARIZATION tag was not found")) goto 100
                  if (i_access(restf)) then
                    dsize = sizeof_double
                    ndata = 1
                    call readf(el%o%spin_polarization,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
                  end if
                  if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%spin_polarization)
                end if
              end if
            end select
            if (error(el%o%spin_polarization < 0.0_double,"ERROR: spin_polarization < 0")) goto 100

          end select

        case (UNIFORM)

          ! Set null values for irrelevant elements
          el%o%use_free_energy = .false.
          el%o%polarization_method = NA
          el%o%spin_polarization = 0.0_double

        end select

        ! Determine the wavefunctions tolerance - amenable to change during a restart
        call arg("wavefunctions_tolerance",el%o%res_norm_tol,found)
        if (found) then
          if (error(el%o%res_norm_tol < 0.0_double,"ERROR: wavefunctions_tolerance < 0")) goto 100
        else
          el%o%res_norm_tol = 1.0e-4_double  ! default value
          if (present(restf)) then
            if (i_access(restf)) tios = findfirsttag(restf,"WAVEFUNCTIONS_TOLERANCE")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: WAVEFUNCTIONS_TOLERANCE tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_double
              ndata = 1
              call readf(el%o%res_norm_tol,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%res_norm_tol)
          end if
        end if

        ! Close the PARAMETERS block
100     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if
        if (error()) goto 200

        ! Divide the k-points among processes
        if (error(nk < mpi_nkgroups(),"ERROR: nk is less than kgroups")) then
          call notify("Number of k-points = ",nk)
          goto 200
        end if
        if (mod(nk,mpi_nkgroups()) /= 0) then
          call warn("WARNING: non-equal division of k-points among kgroups")
        end if
        allocate( el%o%kgroup_index(nk) )
        do ik = 1,nk
          el%o%kgroup_index(ik) = mod(ik-1,mpi_nkgroups()) + 1
        end do

        ! Form the common hamiltonian
        call my(common_hamiltonian(ext,gp,el%o%cutoff,nb),el%o%hc) ; if (error()) goto 200

        ! Determine the starting wavefunctions
        allocate( el%o%wf(nk) )
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"WAVEFUNCTIONS_START")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: WAVEFUNCTIONS_START tag was not found")) goto 200
        end if
        if (has_dyadic_potential(gp)) call my(dyad_kpoint(),dk)
        do ik = 1,nk
          kpt = x_kpoint(el%o%kpoints,ik) ; if (error()) exit
          if (mpi_mykgroup() == el%o%kgroup_index(ik)) then
            usage = "normal"
          else
            usage = "auxiliary"
          end if
          if (present(restf)) then
            call my(wavefunctions_es(el%o%hc,multibasis(usage,el%o%cutoff,nb,kpt,lay),restf=restf),el%o%wf(ik)) ; if (error()) exit
          else
            if (has_dyadic_potential(gp)) then
               dk = x_dyad_kpoint(x_dyadic_potential(gp),ik)
               call my(wavefunctions_es(el%o%hc,multibasis(usage,el%o%cutoff,nb,kpt,lay),dk),el%o%wf(ik)) ; if (error()) exit
            else
               call my(wavefunctions_es(el%o%hc,multibasis(usage,el%o%cutoff,nb,kpt,lay)),el%o%wf(ik)) ; if (error()) exit
            end if
          end if
        end do
        if (has_dyadic_potential(gp)) call glean(thy(dk))

        call sync_config_process_errors() ; if (error()) goto 200

        ! Eigenvalues, occupations, and density
        allocate( el%o%eigs(nk,nb) )
        nullify( el%o%occs )
        call eigenvalues_i(el%o)                ; if (error()) goto 200
        call occupations_i(el%o)                ; if (error()) goto 200
        call my(gen_density(ext),el%o%density)  ; if (error()) goto 200
        call accumulate_density_i(el%o,ext)     ; if (error()) goto 200

        ! OEP eigenvalue corrections
        nullify( el%o%oep_gksc )

        ! Residual norm and energy contributions
        call residual_norm_i(el%o)
        call kinetic_energy_i(el%o)
        call mermin_energy_i(el%o)

        if (present(restf)) then
          call diary_construction_i(el%o,restf)
        else
          call diary_construction_i(el%o)
        end if

        ! Close the ELECTRONS block
200     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

300     call glean(thy(lay))

        call glean(thy(ext))
        call glean(thy(gp))
        if (present(restf)) call glean(thy(restf))

999     if (error("Exit electrons_sc_mod::constructor_el")) continue

        if (.not.error()) call stop_timer("electrons_sc: constructor")

      end function

      subroutine update_el(el,ext,gp)
!doc$ subroutine update(el,ext,gp)
        type(electrons_sc_obj) :: el
        type(external_obj) :: ext
        type(gen_potential_obj) :: gp
!       requires: Consistent ext and gp.
!       modifies: el
!       effects: Updates el with respect to its dependencies.
!       errors: Passes errors.

!cod$
        logical :: e_change, ge_change
        integer :: ik
        type(ghost) :: oldg
        type(dyad_kpoint_obj) :: dk

        call start_timer("electrons_sc: update")

        call my(el)
        call my(ext)
        call my(gp)

        call own_i(el)

        if (x_ghost(ext) /= el%o%g_external) el%o%g_external = x_ghost(ext)

        call update(el%o%hc,ext,gp)
        if (has_dyadic_potential(gp)) call my(dyad_kpoint(),dk)
        do ik = 1,size(el%o%wf)
          oldg = x_ghost(el%o%wf(ik))
          if (has_dyadic_potential(gp)) then
             dk = x_dyad_kpoint(x_dyadic_potential(gp),ik)
             call update(el%o%wf(ik),el%o%hc,dk = dk) ; if (error()) exit
          else
             call update(el%o%wf(ik),el%o%hc) ; if (error()) exit
          end if
          e_change = (oldg /= x_ghost(el%o%wf(ik)))
        end do
        if (has_dyadic_potential(gp)) call glean(thy(dk))

        call sync_config_process_errors() ; if (error()) goto 100

        call xcomm_allreduce(XKGROUP,MPI_LOR,e_change,ge_change)
        if (ge_change) then
          el%o%g = x_ghost()
          call eigenvalues_i(el%o)             ; if (error()) goto 100
          call occupations_i(el%o)             ; if (error()) goto 100
          el%o%density = gen_density(ext)      ; if (error()) goto 100
          call accumulate_density_i(el%o,ext)  ; if (error()) goto 100
        end if
        call residual_norm_i(el%o)
        call kinetic_energy_i(el%o)
        call mermin_energy_i(el%o)

100     call glean(thy(el))
        call glean(thy(ext))
        call glean(thy(gp))

        if (error("Exit electrons_sc_mod::update_el")) continue

        if (.not.error()) call stop_timer("electrons_sc: update")

      end subroutine

      subroutine my_el(el)
!doc$ subroutine my(el)
        type(electrons_sc_obj) :: el

!cod$
        el%ref = el%ref + 1
        el%o%ref = el%o%ref + 1
      end subroutine

      subroutine my_new_el(eli,el)
!doc$ subroutine my(eli,el)
        type(electrons_sc_obj) :: eli, el

!cod$
        el%ref = 1
        el%o => eli%o
        el%o%ref = el%o%ref + 1
      end subroutine

      function thy_el(el) result(elo)
!doc$ function thy(el) result(elo)
        type(electrons_sc_obj) :: el, elo

!cod$
        el%ref = el%ref - 1
        el%o%ref = el%o%ref - 1
        elo%ref = el%ref
        elo%o => el%o
      end function

      subroutine glean_el(el)
!doc$ subroutine glean(el)
        type(electrons_sc_obj) :: el

!cod$
        integer :: ik
        if (el%o%ref < 1) then
          if (associated( el%o%kgroup_index )) deallocate( el%o%kgroup_index )
          if (associated( el%o%eigs )) deallocate( el%o%eigs )
          if (associated( el%o%occs )) deallocate( el%o%occs )
          if (associated( el%o%oep_gksc )) deallocate( el%o%oep_gksc )
          if (associated( el%o%wf )) then
            do ik = 1,size(el%o%wf)
              call glean(thy(el%o%wf(ik)))
            end do
            deallocate( el%o%wf )
          end if
          call glean(thy(el%o%kpoints))
          call glean(thy(el%o%density))
          call glean(thy(el%o%hc))
          deallocate( el%o )
        end if
      end subroutine

      subroutine bequeath_el(el)
!doc$ subroutine bequeath(el)
        type(electrons_sc_obj) :: el

!cod$
        continue
      end subroutine

      subroutine assign_el(el,el2)
!doc$ subroutine assign(el,el2)
        type(electrons_sc_obj), intent(inout) :: el
        type(electrons_sc_obj), intent(in) :: el2

!cod$
        type(electrons_sc_obj) :: elt
        call my(el2)
        elt%o => el%o
        el%o%ref = el%o%ref - el%ref
        el%o => el2%o
        el%o%ref = el%o%ref + el%ref
        call glean(elt)
        call glean(thy(el2))
      end subroutine

      function wormhole_el(el) result(elr)
!doc$ function wormhole(el) result(elr)
        type(electrons_sc_obj) :: el
        type(electrons_sc_rep), pointer :: elr
!       effects: Points elr at el%o
!       errors: Wormhole is an implementation dependent thing and you should know what you're doing.

!cod$
        call my(el)
        elr => el%o
        call glean(thy(el))
      end function

      function el_ref(el) result(r)
!doc$ function x_ref(el) result(r)
        type(electrons_sc_obj) :: el
        integer, dimension(2) :: r
!       effects: Returns el%ref and el%o%ref.

!cod$
        r(1) = el%ref
        r(2) = el%o%ref
        call glean(el)
      end function

      function el_ghost(el) result(g)
!doc$ function x_ghost(el) result(g)
        type(electrons_sc_obj) :: el
        type(ghost) :: g
!       effects: Returns ghost of el.

!cod$
        call my(el)
        g = el%o%g
        call glean(thy(el))
      end function

      function el_fermi_level(el) result(ef)
!doc$ function x_fermi_level(el) result(ef)
        type(electrons_sc_obj) :: el
        real(double) :: ef
!       effects: Returns the Fermi level of el.

!cod$
        call my(el)
        ef = el%o%fermi_level
        call glean(thy(el))
      end function

      function el_energy(el) result(e)
!doc$ function x_energy(el) result(e)
        type(electrons_sc_obj) :: el
        real(double) :: e
!       effects: Returns the electrons energy.

!cod$
        call my(el)
        e = el%o%kinetic_energy
        if (el%o%use_free_energy) e = e + el%o%mermin_energy
        call glean(thy(el))
      end function

      function el_kinetic_energy(el) result(e)
!doc$ function x_kinetic_energy(el) result(e)
        type(electrons_sc_obj) :: el
        real(double) :: e
!       effects: Returns the kinetic energy.

!cod$
        call my(el)
        e = el%o%kinetic_energy
        call glean(thy(el))
      end function

      function el_mermin_energy(el) result(e)
!doc$ function x_mermin_energy(el) result(e)
        type(electrons_sc_obj) :: el
        real(double) :: e
!       effects: Returns the entropy contribution to the Mermin free energy

!cod$
        call my(el)
        e = el%o%mermin_energy
        call glean(thy(el))
      end function

      function el_wavefunctions(el,ik) result(wf)
!doc$ function x_wavefunctions(el,ik) result(wf)
        type(electrons_sc_obj) :: el
        integer, intent(in) :: ik
        type(wavefunctions_es_obj) :: wf
!       effects: Returns the ik set of wavefunctions.
!       errors: ik out of range.

!cod$
        call my(el)
        if ( error((ik < 1) .or. (ik > size(el%o%wf)),"ERROR: ik is out of range")) goto 100
        call my(el%o%wf(ik),wf)
        call bequeath(thy(wf))
100     call glean(thy(el))
        if (error("Exit electrons_sc_mod::el_wavefunctions")) continue
      end function 

      function el_density(el) result(dens)
!doc$ function  x_density(el) result(dens)
        type(electrons_sc_obj) :: el
        type(gen_density_obj) :: dens
!       effects: Returns the density object
!cod$
        call my(el)
        call my(el%o%density,dens)
        call bequeath(thy(dens))
        call glean(thy(el))
      end function 

      function el_kpoints(el) result(kp)
!doc$ function x_kpoints(el) result(kp)
        type(electrons_sc_obj) :: el
        type(kpoints_obj) :: kp
!       effects: Returns the kpoints_obj of el.

!cod$
        call my(el)
        call my(el%o%kpoints,kp)
        call glean(thy(el))
        call bequeath(thy(kp))
      end function
 
      function el_kt(el) result(kt)
!doc$ function x_kt(el) result(kt)
        type(electrons_sc_obj) :: el
        real(double) :: kt
!       effects: Returns the temperature (in Rydbergs) used to calculate thermal occupations.

!cod$
        call my(el)
        kt = el%o%kt
        call glean(thy(el))
      end function

      function el_n_bands(el) result(n)
!doc$ function x_n_bands(el) result(n)
        type(electrons_sc_obj) :: el
        integer :: n
!       effects: Returns the number of bands.

!cod$
        call my(el)
        n = size(el%o%eigs,2)
        call glean(thy(el))
      end function

      function el_kgroup_index(el,ik) result(ikg)
!doc$ function x_kgroup_index(el,ik) result(ikg)
        type(electrons_sc_obj) :: el
        integer, intent(in) :: ik
        integer :: ikg
!       effects: Returns the kgroup index of ik.

!cod$
        call my(el)
        ikg = el%o%kgroup_index(ik)
        call glean(thy(el))
      end function

      function el_eigenvalue(el,ik,ib) result(ev)
!doc$ function x_eigenvalue(el,ik,ib) result(ev)
        type(electrons_sc_obj) :: el
        integer :: ik, ib
        real(double) :: ev
!       effects: Returns the eigenvalue for k-point ik and band ib.
!       errors: ik or ib out of range.

!cod$
        call my(el)
        if (error((ik < 1) .or. (ik > size(el%o%eigs,1)),"ERROR: ik is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(el%o%eigs,2)),"ERROR: ib is out of range")) goto 100
        ev = el%o%eigs(ik,ib)
100     call glean(thy(el))
        if (error("Exit electrons_sc_mod::el_eigenvalue")) continue
      end function

      function el_eigenvalues(el) result(evs)
!doc$ function x_eigenvalues(el) result(evs)
        type(electrons_sc_obj) :: el
        real(double), dimension(size(el%o%eigs,1),size(el%o%eigs,2)) :: evs
!       effects: Returns the eigenvalues.

!cod$
        call my(el)
        evs = el%o%eigs
        call glean(thy(el))
      end function

      function el_occupation(el,ik,ib) result(occ)
!doc$ function x_occupation(el,ik,ib) result(occ)
        type(electrons_sc_obj) :: el
        integer :: ik, ib
        real(double) :: occ
!       effects: Returns the occupation for k-point ik and band ib.
!       errors: ik or ib out of range.

!cod$
        call my(el)
        if (error((ik < 1) .or. (ik > size(el%o%occs,1)),"ERROR: ik is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(el%o%occs,2)),"ERROR: ib is out of range")) goto 100
        occ = el%o%occs(ik,ib)
100     call glean(thy(el))
        if (error("Exit electrons_sc_mod::el_occupation")) continue
      end function

      function el_occupations(el) result(occs)
!doc$ function x_occupations(el) result(occs)
        type(electrons_sc_obj) :: el
        real(double), dimension(size(el%o%occs,1),size(el%o%occs,2)) :: occs
!       effects: Returns the occupations.

!cod$
        call my(el)
        occs = el%o%occs
        call glean(thy(el))
      end function

      function el_residual_norm(el) result(rn)
!doc$ function x_residual_norm(el) result(rn)
        type(electrons_sc_obj) :: el
        real(double) :: rn
!       effects: Returns the weighted sum of the wavefunction residual norms.

!cod$
        call my(el)
        rn = el%o%res_norm
        call glean(thy(el))
      end function

      function el_converged(el) result(cvg)
!doc$ function x_converged(el) result(cvg)
        type(electrons_sc_obj) :: el
        logical :: cvg
!       effects: Returns the convergence status of el with respect to res_norm and res_norm_tol.

!cod$
        call my(el)
        cvg = el%o%res_norm_cvg
        call glean(thy(el))
      end function

      subroutine el_oep_gksc(el,ec)
!doc$ subroutine s_oep_gksc(el,ec)
        type(electrons_sc_obj) :: el
        real(double), dimension(size(el%o%eigs,1),size(el%o%eigs,2)) :: ec
!       effects: Sets the eigenvalue corrections from oep_mod.
!       requires: ec dimensions be the same as el%o%eigs.

!cod$
        call my(el)
        if (.not.associated( el%o%oep_gksc )) allocate( el%o%oep_gksc(size(el%o%eigs,1),size(el%o%eigs,2)) )
        el%o%oep_gksc = ec
        call glean(thy(el))
      end subroutine

      function thermal_occupations_el(el) result(tho)
!doc$ function thermal_occupations(el) result(tho)
        type(electrons_sc_obj) :: el
        logical :: tho
!       effects: Returns .true. if thermal occupations are being used and false otherwise

!cod$
        call my(el)
        tho = .false.
        select case (el%o%occupation_method)
        case (THERMAL)
          tho = .true.
        end select
        call glean(thy(el))
      end function

      subroutine forces_el(el,f)
!doc$ subroutine forces(el,f)
        type(electrons_sc_obj) :: el
        real(double), dimension(:,:), intent(out) :: f
!       modifies: f
!       effects: Returns atomic forces due to electrons.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double), dimension(:), allocatable :: wts
        real(double), dimension(:,:), allocatable :: f_k, f_kg, f_sg

        call my(el)

        allocate( wts(size(el%o%eigs,2)) )
        allocate( f_k(size(f,1),size(f,2)), f_kg(size(f,1),size(f,2)), f_sg(size(f,1),size(f,2)) )

        f_kg = 0.0_double
        do ik = 1,size(el%o%wf)
          wts = el%o%occs(ik,:)*x_kweight(el%o%kpoints,ik)
          call forces(el%o%wf(ik),wts,f_k) ; if (error()) goto 100
          f_kg = f_kg + f_k
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,f_kg,f_sg) ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_SUM,f_sg,f)    ; if (error()) goto 100

100     if (allocated( wts )) deallocate( wts )
        if (allocated( f_k ))  deallocate( f_k )
        if (allocated( f_kg )) deallocate( f_kg )
        if (allocated( f_sg )) deallocate( f_sg )

        call glean(thy(el))

        if (error("Exit electrons_sc_mod::forces_el")) continue

      end subroutine

      subroutine pressure_el(el,p)
!doc$ subroutine pressure(el,p)
        type(electrons_sc_obj) :: el
        real(double), intent(out) :: p
!       effects: Returns pressure contributions due to electrons.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double) :: p_k, p_kg, p_sg
        real(double), dimension(:), allocatable :: wts

        call my(el)

        allocate( wts(size(el%o%occs,2)) )

        p_kg = 0.0_double
        do ik = 1,size(el%o%wf)
          wts = el%o%occs(ik,:)*x_kweight(el%o%kpoints,ik)
          call pressure(el%o%wf(ik),wts,p_k) ; if (error()) goto 100
          p_kg = p_kg + p_k
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,p_kg,p_sg) ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_SUM,p_sg,p)    ; if (error()) goto 100

100     if (allocated( wts )) deallocate( wts )

        call glean(thy(el))

        if (error("Exit electrons_sc_mod::pressure_el")) continue

      end subroutine

      subroutine stress_tensor_el(el,s)
!doc$ subroutine stress_tensor(el,s)
        type(electrons_sc_obj) :: el
        real(double), dimension(:,:), intent(out) :: s
!       requires: s be dimension(3,3)
!       effects: Returns stress tensor contributions due to electrons.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double), dimension(:), allocatable :: wts
        real(double), dimension(3,3) :: s_k, s_kg, s_sg

        call my(el)

        allocate( wts(size(el%o%occs,2)) )

        s_kg = 0.0_double
        do ik = 1,size(el%o%wf)
          wts = el%o%occs(ik,:)*x_kweight(el%o%kpoints,ik)
          call stress_tensor(el%o%wf(ik),wts,s_k) ; if (error()) goto 100
          s_kg = s_kg + s_k
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,s_kg,s_sg) ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_SUM,s_sg,s)    ; if (error()) goto 100

100     if (allocated( wts )) deallocate( wts )

        call glean(thy(el))

        if (error("Exit electrons_sc_mod::stress_tensor_el")) continue

      end subroutine

      subroutine diary_el(el)
!doc$ subroutine diary(el)
        type(electrons_sc_obj) :: el
!       effects: Writes el information to the diary.

!cod$
        logical :: found
        character(line_len) :: tag
        integer :: ib, ik, ios, is, msg, nb, nk, nsg, file_format
        integer :: band_range(2), brl, bru, kpoint_range(2),krl, kru
        real(double) :: bocc1, bocc2, ne1, ne2, kpt(3), cv
        real(double), dimension(:), allocatable :: fl_c, fl_sg
        real(double), dimension(:,:), allocatable :: eigs_gks
        real(double), dimension(:,:,:), allocatable :: eigs_c, eigs_sg
        real(double), dimension(:,:,:), allocatable :: occs_c, occs_sg
        real(double), dimension(:,:,:), allocatable :: eigs_gks_c, eigs_gks_sg
        real(double), pointer    :: rtmp(:,:,:)
        complex(double), pointer :: ctmp(:,:,:)
        type(grid_obj) :: g
        type(file_obj) :: f
        !type(multivector_obj) :: mvec
        character(line_len)    :: filename
        integer               :: kind

        call my(el)

        nk = size(el%o%eigs,1)
        nb = size(el%o%eigs,2)
        nsg = mpi_nsgroups()

        ! Compute the generalized Kohn-Sham eigenvalues
        if (associated( el%o%oep_gksc )) then
          allocate( eigs_gks(nk,nb) )
          eigs_gks = el%o%eigs + el%o%oep_gksc
        end if

        ! Gather spin information to the CONFIG rank 0 process.
        select case (nsg)
        case (2)
          msg = mpi_mysgroup()
          ! Fermi level
          allocate( fl_c(nsg), fl_sg(nsg) )
          fl_sg = 0.0_double
          fl_sg(msg) = el%o%fermi_level
          call xcomm_reduce(XSGROUP,MPI_SUM,fl_sg,fl_c)                                             ; if (error()) goto 300
          ! eigenvalues
          allocate( eigs_c(nsg,nk,nb), eigs_sg(nsg,nk,nb) )
          eigs_sg = 0.0_double
          eigs_sg(msg,:,:) = el%o%eigs
          call xcomm_reduce(XSGROUP,MPI_SUM,eigs_sg,eigs_c)                                         ; if (error()) goto 300
          ! occupations
          allocate( occs_c(nsg,nk,nb), occs_sg(nsg,nk,nb) )
          occs_sg = 0.0_double
          occs_sg(msg,:,:) = el%o%occs
          call xcomm_reduce(XSGROUP,MPI_SUM,occs_sg,occs_c)                                         ; if (error()) goto 300
          ! gks eigenvalues
          if (associated( el%o%oep_gksc )) then
            allocate( eigs_gks_c(nsg,nk,nb), eigs_gks_sg(nsg,nk,nb) )
            eigs_gks_sg = 0.0_double
            eigs_gks_sg(msg,:,:) = el%o%eigs + el%o%oep_gksc
            call xcomm_reduce(XSGROUP,MPI_SUM,eigs_gks_sg,eigs_gks_c)                               ; if (error()) goto 300
          end if
        end select

        ! Write information to the diary file
        if (i_access( diaryfile() )) then

          ! Kohn-Sham eigenvalues and occupations
          select case (nsg)
          case (1)
            select case (el%o%occupation_method)
            case (THERMAL)
              write(x_unit(diaryfile()),'(/,t4,"Fermi level = ",f14.10," Ryd")') el%o%fermi_level
            end select
            do ik = 1,size(el%o%eigs,1)
              write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
              write(x_unit(diaryfile()),'(/,21x,"band",8x,"eigenvalue (Ryd)",7x,"occupation")')
              write(x_unit(diaryfile()),'(19x,49("-"))')
              do ib = 1,size(el%o%eigs,2)
                write(x_unit(diaryfile()),'(21x,i4,9x,f14.10,11x,f6.4)') ib, el%o%eigs(ik,ib), el%o%occs(ik,ib)
              end do
            end do
          case (2)
            select case (el%o%occupation_method)
            case (THERMAL)
              write(x_unit(diaryfile()),'(/,t4,"Fermi level:",10x,f14.10,15x,f14.10," Ryd")') fl_c(1), fl_c(2)
            end select
            do ik = 1,size(eigs_c,2)
              write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
              write(x_unit(diaryfile()),'(/,11x,"band",8x,"eigenvalue (Ryd)",7x,"occupation", &
                                                    & 10x,"eigenvalue (Ryd)",7x,"occupation")')
              write(x_unit(diaryfile()),'(9x,92("-"))')
              do ib = 1,size(eigs_c,3)
                write(x_unit(diaryfile()),'(11x,i4,9x,f14.10,11x,f6.4,13x,f14.10,11x,f6.4)') ib, eigs_c(1,ik,ib), occs_c(1,ik,ib), &
                                                                                               & eigs_c(2,ik,ib), occs_c(2,ik,ib)
              end do
            end do
          end select

          ! Kohn-Sham band occupations
          call arglc("write_band_occupations",tag,found)
          if (.not.found) tag = "off"
          select case (trim(tag))
          case ("on")
            write(x_unit(diaryfile()),'(/,t6,"Band occupations:")')
            select case (nsg)
            case (1)
              write(x_unit(diaryfile()),'(/,21x,"band",8x,"occupation")')
              write(x_unit(diaryfile()),'(19x,26("-"))')
              do ib = 1,size(el%o%occs,2)
                bocc1 = 0.0_double
                do ik = 1,size(el%o%occs,1)
                  bocc1 = bocc1 + x_kweight(el%o%kpoints,ik)*el%o%occs(ik,ib)
                end do
                write(x_unit(diaryfile()),'(21x,i4,9x,f6.4)') ib, bocc1
              end do
            case (2)
              write(x_unit(diaryfile()),'(/,11x,"band",8x,"occupation",10x,"occupation")')
              write(x_unit(diaryfile()),'(9x,46("-"))')
              ne1 = 0.0_double
              ne2 = 0.0_double
              do ib = 1,size(occs_c,3)
                bocc1 = 0.0_double
                bocc2 = 0.0_double
                do ik = 1,size(el%o%occs,1)
                  bocc1 = bocc1 + x_kweight(el%o%kpoints,ik)*occs_c(1,ik,ib)
                  bocc2 = bocc2 + x_kweight(el%o%kpoints,ik)*occs_c(2,ik,ib)
                end do
                ne1 = ne1 + bocc1
                ne2 = ne2 + bocc2
                write(x_unit(diaryfile()),'(11x,i4,9x,f6.4,14x,f6.4)') ib, bocc1, bocc2
              end do
              write(x_unit(diaryfile()),'(/,t6,"Total occupations:")')
              write(x_unit(diaryfile()),'(/,23x,"spin group 1",8x,"spin group 2",9x,"|difference|")')
              write(x_unit(diaryfile()),'(21x,57("-"))')
              write(x_unit(diaryfile()),'(21x,f9.4,11x,f9.4,13x,f9.4)') ne1, ne2, abs(ne1 - ne2)
            end select
          end select

          ! Generalized Kohn-Sham eigenvalues
          if (allocated( eigs_gks )) then
            write(x_unit(diaryfile()),'(/,t4,"Generalized Kohn-Sham eigenvalues")')
            select case (nsg)
            case (1)
              do ik = 1,size(eigs_gks,1)
                write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
                write(x_unit(diaryfile()),'(/,21x,"band",8x,"eigenvalue (Ryd)")')
                write(x_unit(diaryfile()),'(19x,31("-"))')
                do ib = 1,size(eigs_gks,2)
                  write(x_unit(diaryfile()),'(21x,i4,9x,f14.10)') ib, eigs_gks(ik,ib)
                end do
              end do
            case (2)
              do ik = 1,size(eigs_gks_c,2)
                write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
                write(x_unit(diaryfile()),'(/,11x,"band",8x,"eigenvalue (Ryd)",10x,"eigenvalue (Ryd)")')
                write(x_unit(diaryfile()),'(9x,58("-"))')
                do ib = 1,size(eigs_gks_c,3)
                  write(x_unit(diaryfile()),'(11x,i4,9x,f14.10,13x,f14.10)') ib, eigs_gks_c(1,ik,ib), eigs_gks_c(2,ik,ib)
                end do
              end do
            end select
          end if

        end if

        ! Write the eigenvalues and occupations to a file
        call arglc("write_eigenvalues",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on")
          call my(file(trim(eigenvalues_path)),f)
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to open eigenvalues file")) goto 100
          if (i_access(f)) then
            select case (nsg)
            case (1)
              do ik = 1,size(el%o%eigs,1)
                write(x_unit(f),'(/,"Special k-point #",i0,":",/)') ik
                do ib = 1,size(el%o%eigs,2)
                  write(x_unit(f),'(1x,i5,10x,f15.10,10x,f12.10)') ib, el%o%eigs(ik,ib), el%o%occs(ik,ib)
                end do
              end do
            case ( 2)
              do is = 1,size(eigs_c,1)
                write(x_unit(f),'(/,"Spin #",i0,":",/)') is
                do ik = 1,size(eigs_c,2)
                  write(x_unit(f),'(/,"Special k-point #",i0,":",/)') ik
                  do ib = 1,size(eigs_c,3)
                    write(x_unit(f),'(1x,i5,10x,f15.10,10x,f12.10)') ib, eigs_c(is,ik,ib), occs_c(is,ik,ib)
                  end do
                end do
              end do
            end select
            close(x_unit(f))
          end if
100       call glean(thy(f))
        end select

        ! Write the gks eigenvalues to a file
        call arglc("write_eigenvalues",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on")
          call my(file(trim(gks_eigenvalues_path)),f)
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to open gks_eigenvalues file")) goto 200
          if (i_access(f)) then
            select case (nsg)
            case (1)
              do ik = 1,size(eigs_gks,1)
                write(x_unit(f),'(/,"Special k-point #",i0,":",/)') ik
                do ib = 1,size(eigs_gks,2)
                  write(x_unit(f),'(1x,i5,10x,f15.10)') ib, eigs_gks(ik,ib)
                end do
              end do
            case ( 2)
              do is = 1,size(eigs_gks_c,1)
                write(x_unit(f),'(/,"Spin #",i0,":",/)') is
                do ik = 1,size(eigs_gks_c,2)
                  write(x_unit(f),'(/,"Special k-point #",i0,":",/)') ik
                  do ib = 1,size(eigs_gks_c,3)
                    write(x_unit(f),'(1x,i5,10x,f15.10)') ib, eigs_gks_c(is,ik,ib)
                  end do
                end do
              end do
            end select
            close(x_unit(f))
          end if
200       call glean(thy(f))
        end select

        ! Write the band structure to a file
        call arglc("write_band_structure",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on")
          call my(file(trim(band_structure_path)),f)
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to open band_structure file")) goto 400
          if (i_access(f)) then
            select case (nsg)
            case (1)
              write(x_unit(f),*) '% ik   kx   ky   kz  Efermi  eval1  eval2  eval3 ...'
              do ik = 1,size(el%o%eigs,1)
                kpt = x_kpoint(el%o%kpoints,ik)
                write(x_unit(f),*) ik,kpt(1),kpt(2),kpt(3),el%o%fermi_level, el%o%eigs(ik,:)
              end do
            case ( 2)
              do is = 1,size(eigs_c,1)
                write(x_unit(f),'(/,"Spin #",i0,":",/)') is
                write(x_unit(f),*) '% ik   kx   ky   kz  Efermi  eval1  eval2  eval3 ...'
                do ik = 1,size(eigs_c,2)
                  kpt = x_kpoint(el%o%kpoints,ik)
                  write(x_unit(f),*) ik,kpt(1),kpt(2),kpt(3),el%o%fermi_level, eigs_c(is,ik,:)
                end do
              end do
            end select
            close(x_unit(f))
          end if
400       call glean(thy(f))
        end select

        ! Read in the range of bands that should be written out
        call arg("write_band_range",band_range,found)
        if (found) then
           brl = band_range(1)
           bru = band_range(2)
           if (error((brl<1),"Error: write_band_range(1) is not valid (<1)")) goto 300
           if (error((nb<brl),"Error: write_band_range(1) is not valid (>nb)")) goto 300
           if (error((bru<1),"Error: write_band_range(2) is not valid (<1)")) goto 300
           if (error((nb<bru),"Error: write_band_range(2) is not valid (>nb)")) goto 300
           if (error((bru<brl),"Error: write_band_range is not valid (bru<brl)")) goto 300
        else
           brl = 1
           bru = nb
        end if
        
        ! Read in the range of kpoint that should be written out
        call arg("write_kpoint_range",kpoint_range,found)
        if (found) then
           krl = kpoint_range(1)
           kru = kpoint_range(2)
           if (error((krl<1),"Error: write_band_range(1) is not valid (<1)")) goto 300
           if (error((nk<krl),"Error: write_band_range(1) is not valid (>nk)")) goto 300
           if (error((kru<1),"Error: write_band_range(2) is not valid (<1)")) goto 300
           if (error((nk<kru),"Error: write_band_range(2) is not valid (>nk)")) goto 300
           if (error((kru<krl),"Error: write_band_range is not valid (kru<krl)")) goto 300
        else
           krl = 1
           kru = nk
        end if

        ! Read in the file format for any output grid files
        call arglc("write_file_format",tag,found)
        if (.not.found) tag = "matlab"
        select case (trim(tag))
        case ("matlab","mat","m")
           file_format = MATLAB
        case("amira","am","avizo","av","a")
           file_format = AMIRA
        case("vtk","visualizationtoolkit","v")
           file_format = VTK
        case default
           call warn("Unrecognized entry for write_file_format, setting to MATLAB")
           file_format = MATLAB
        end select
          
          ! Write out (a subset of) the wavefunctions to file
        call arglc("write_wavefunctions",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.","t","yes","true")


          ! Read in whether the wavefunction should be written out in real or reciprocal space
          call arglc("write_wavefunctions_rep",tag,found)
          if (.not.found) tag = "realspace"
          select case (trim(tag))
          case("realspace","real","position","pos","p")
             kind = CSP_KIND
          case("reciprocalspace","reciprocal","fourier","f")
             kind = CSF_KIND
          case default
             if (error(.true.,"Error: unrecognized write_wavefunctions_rep tag")) goto 300
          end select

          ! Write out wavefunctions
          call my(grid(x_layout(x_grid_density(el%o%density)),KGROUP),g)
          do ik=krl,kru
!            call my(x_multivector(el%o%wf(ik)),mvec) ; if (error()) goto 300
            kpt = x_kpoint(el%o%kpoints,ik)
            do ib=brl,bru

              select case (nsg)
              case (1)
                call get_wf_filename_i(filename,file_format,ik,nk,ib,nb)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                !call put(mvec,ib,g)
                call put(x_multivector(el%o%wf(ik)),ib,g)
                call merge_grid_density(g)

                ! Write to file
                call write_to_file(g,filename,file_format,kind,kpt=kpt) 
                if (error()) goto 300

              case (2)
                ! First write out spin 1
                call get_wf_filename_i(filename,file_format,ik,nk,ib,nb,is=1)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                !call put(mvec,ib,g)
                call put(x_multivector(el%o%wf(ik)),ib,g)
                call merge_grid_density(g)

                call write_to_file(g,filename,file_format,kind,spin=1,kpt=kpt) 
                if (error()) goto 300

                ! Next write out spin 2
                call get_wf_filename_i(filename,file_format,ik,nk,ib,nb,is=2)

                call write_to_file(g,filename,file_format,kind,spin=2,kpt=kpt) 
                if (error()) goto 300
              end select
              ! Reset scope to kgroup for next iteration
              call sgroup_to_kgroup(g)
            end do
            !call glean(thy(mvec))
          end do
          call glean(thy(g))

        end select



        ! Write out (a subset of) the state_densities to file
        call arglc("write_state_densities",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.","t","yes","true")
          ! Read in whether the wavefunction should be written out in real or reciprocal space
          call arglc("write_state_densities_rep",tag,found)
          if (.not.found) tag = "realspace"
          select case (trim(tag))
          case("realspace","real","position","pos","p")
             kind = RS_KIND
          case("reciprocalspace","reciprocal","fourier","f")
             kind = CSF_KIND
          case default
             if (error(.true.,"Error: unrecognized write_state_densities_rep tag")) goto 300
          end select

          ! Write out state_densities
          call my(grid(x_layout(x_grid_density(el%o%density)),KGROUP),g)
          cv = x_cell_volume(x_lattice(x_layout(g)))
          call alloc(rtmp, x_layout(g), S_TYPE)

          do ik=krl,kru
            !call my(x_multivector(el%o%wf(ik)),mvec) ; if (error()) goto 300
            kpt = x_kpoint(el%o%kpoints,ik)
            do ib=brl,bru

              select case (nsg)
              case (1)
                call get_sd_filename_i(filename,file_format,ik,nk,ib,nb)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                !call put(mvec,ib,g)
                call put(x_multivector(el%o%wf(ik)),ib,g)
                call merge_grid_density(g)

                ! extract the wfn, multiply it by its cconjg to obtain the state density
                call take(ctmp,g,CSP_KIND)
                rtmp = conjg(ctmp)*ctmp !/cv
                if (associated(ctmp)) deallocate(ctmp)
                nullify(ctmp)

                ! put the state density back in the grid object
                call put(rtmp,g,RS_KIND)

                ! Write to file
                call write_to_file(g,filename,file_format,kind,kpt=kpt) ; if (error()) goto 300

              case (2)
                ! First write out spin 1
                call get_sd_filename_i(filename,file_format,ik,nk,ib,nb,is=1)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                !call put(mvec,ib,g)
                call put(x_multivector(el%o%wf(ik)),ib,g)
                call merge_grid_density(g)

                ! extract the wfn, multiply it by its cconjg to obtain the state density
                call take(ctmp,g,CSP_KIND)
                rtmp = real(conjg(ctmp)*ctmp) !/cv
                if (associated(ctmp)) deallocate(ctmp)
                nullify(ctmp)


                ! put the state density back in the grid object
                call put(rtmp,g,RS_KIND)

                call write_to_file(g,filename,file_format,kind,spin=1,kpt=kpt) ; if (error()) goto 300

                ! Next write out spin 2
                call get_sd_filename_i(filename,file_format,ik,nk,ib,nb,is=2)

                call write_to_file(g,filename,file_format,kind,spin=2,kpt=kpt) ; if (error()) goto 300
              end select

              ! Reset scope to KGROUP for next iteration
              call sgroup_to_kgroup(g)
              ! Take rtmp only to keep from having to re-alloc
              call take(rtmp,g,RS_KIND)
            end do
            !call glean(thy(mvec))
          end do
          call put(rtmp,g,RS_KIND)
          call glean(thy(g))

        end select


        ! Write out the electronic density to file
        call arglc("write_density",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on",".true.","t","yes","true")
          ! Read in whether the density should be written out in real or reciprocal space
          call arglc("write_density_rep",tag,found)
          if (.not.found) tag = "realspace"
          select case (trim(tag))
          case("realspace","real","position","pos","p")
             kind = RS_KIND
          case("reciprocalspace","reciprocal","fourier","f")
             kind = CSF_KIND
          case default
             if (error(.true.,"Error: unrecognized write_density_rep tag")) goto 300
          end select

          ! Grab the grid density
          call my(x_grid_density(el%o%density),g)
           
          select case (nsg)
          case (1)
             select case(file_format)
             case(MATLAB)
                filename = trim("density.mat")
             case(AMIRA)
                filename = trim("density.am")
             case(VTK)
                filename = trim("density.vtk")
             end select
             
             ! Write to file
             call write_to_file(g,filename,file_format,kind) ; if (error()) goto 300
             
          case (2)
             ! Gen filename for spin 1 (call it spin up)
             select case(file_format)
             case(MATLAB)
                filename = trim("density_up.mat")
             case(AMIRA)
                filename = trim("density_up.am")
             case(VTK)
                filename = trim("density_up.vtk")
             end select

             ! Write out spin 1
             call write_to_file(g,filename,file_format,kind,spin=1) ; if (error()) goto 300
             
             ! Gen filename for spin 2 (spin down)
             select case(file_format)
             case(MATLAB)
                filename = trim("density_dn.mat")
             case(AMIRA)
                filename = trim("density_dn.am")
             case(VTK)
                filename = trim("density_dn.vtk")
             end select
             
             ! Next write out spin 2
             call write_to_file(g,filename,file_format,kind,spin=2) ; if (error()) goto 300

          end select
          call glean(thy(g))
        end select
       

 
300     call barrier(CONFIG)

        if (allocated( fl_c ))  deallocate( fl_c )
        if (allocated( fl_sg ))  deallocate( fl_sg )
        if (allocated( eigs_gks ))  deallocate( eigs_gks )
        if (allocated( eigs_c ))  deallocate( eigs_c )
        if (allocated( eigs_sg )) deallocate( eigs_sg )
        if (allocated( occs_c ))  deallocate( occs_c )
        if (allocated( occs_sg )) deallocate( occs_sg )
        if (allocated( eigs_gks_c ))  deallocate( eigs_gks_c )
        if (allocated( eigs_gks_sg )) deallocate( eigs_gks_sg )

        call glean(thy(el))

        if (error("Exit electrons_sc_mod::diary_el")) continue

      end subroutine

      function max_energy_el(el) result(e)
!doc$ function max_energy(el) result(e)
        type(electrons_sc_obj) :: el
        real(double) :: e
!       effects: Returns the maximum among the el energies.

        call my(el)
        if (el%o%use_free_energy) then
          e = max(abs(el%o%kinetic_energy),abs(el%o%mermin_energy))
        else
          e = abs(el%o%kinetic_energy)
        end if
        call glean(thy(el))
      end function

      subroutine diary_energies_el(el,md)
!doc$ subroutine diary_energies(el,md)
        type(electrons_sc_obj) :: el
        integer, intent(in) :: md
!       modifies: Output stream
!       effects: Prints energies.

!cod$
        call my(el)

        if (i_access(diaryfile())) then
          select case (md)
          case (12)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f12.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f12.9)') el%o%mermin_energy
          case (13)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f13.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f13.9)') el%o%mermin_energy
          case (14)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f14.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f14.9)') el%o%mermin_energy
          case (15)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f15.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f15.9)') el%o%mermin_energy
          case (16)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f16.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f16.9)') el%o%mermin_energy
          case (17)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f17.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f17.9)') el%o%mermin_energy
          case default
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f18.9)') el%o%kinetic_energy
            if (el%o%use_free_energy) write(x_unit(diaryfile()),'(t6,"Mermin                = ",f18.9)') el%o%mermin_energy
          end select
        end if

        call glean(thy(el))

      end subroutine

      subroutine decompose_el(el,site_data,mode,f)
!doc$ subroutine decompose(el,site_data,mode,f)
        type(electrons_sc_obj) :: el
        real(double), dimension(:,:), intent(in) :: site_data
        character(line_len), intent(in) :: mode
        type(file_obj) :: f
!       requires: x_unit(f) be open.
!       effects: Decomposes the Kohn-Sham functions into s, p, & d spherical harmonics around sites.
!                Outputs these results to file dcompf in one of three formats.
!       errors: Passes errors.

!cod$
        logical :: found
        integer :: b1, b1_sg, b2, b2_sg, bt
        integer :: ib, ik, is, msg, nb, nk, ns, nsg
        real(double) :: eig, eig1, eig2, occ, occ1, occ2
        real(double) :: sum_s, sum_s1, sum_s2, sum_p, sum_p1, sum_p2, sum_d, sum_d1, sum_d2
        real(double) :: emax, emin, kwt, range
        real(double), dimension(3) :: kpt
        real(double), dimension(:,:,:), allocatable :: eigs_c, eigs_sg, occs_c, occs_sg
        real(double), dimension(:,:,:), allocatable :: rsa, rsa_kg
        real(double), dimension(:,:,:,:), allocatable :: rsa_c, rsa_sg

        call my(el)
        call my(f)

        nk = size(el%o%eigs,1)
        nb = size(el%o%eigs,2)
        ns = size(site_data,2)
        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        ! Read the eigenvalue range to decompose.
        call arg("dcomp_range",range,found)
        if (found) then
          if (error(range < 0.0_double,"ERROR: dcomp_range < 0")) goto 100
          emin = el%o%fermi_level - 0.5_double*range
          emax = el%o%fermi_level + 0.5_double*range
        else
          emin = minval(el%o%eigs)
          emax = maxval(el%o%eigs)
        end if

        ! Convert the eigenvalue range to a band range
        b1_sg = nb
        b2_sg = 1
        do ik = 1,nk
          bt = 1
          do ib = 1,nb
            if (el%o%eigs(ik,ib) >= emin) exit
            bt = ib
          end do
          b1_sg = min(b1_sg,bt)
          bt = nb
          do ib = nb,1,-1
            if (el%o%eigs(ik,ib) <= emax) exit
            bt = ib
          end do
          b2_sg = max(b2_sg,bt)
        end do
        call xcomm_allreduce(XSGROUP,MPI_MIN,b1_sg,b1) ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_MAX,b2_sg,b2) ; if (error()) goto 100
        nb = b2 - b1 + 1
        if (nb < 0) then
          call warn("WARNING: aborting decomposition because nb < 0")
          goto 100
        end if

        allocate(    rsa(9,nb,ns) )
        allocate( rsa_kg(9,nb,ns) )
        select case (nsg)
        case (2)
          allocate(  eigs_c(nsg,nk,nb) )
          allocate( eigs_sg(nsg,nk,nb) )
          allocate(  occs_c(nsg,nk,nb) )
          allocate( occs_sg(nsg,nk,nb) )
          allocate( rsa_sg(nsg,9,nb,ns) )
          allocate(  rsa_c(nsg,9,nb,ns) )
        end select

        do ik = 1,size(el%o%wf)

          ! Decompose the wavefunctions
          call decompose(el%o%wf(ik),site_data,mode,rsa_kg,b1) ; if (error()) goto 100

          ! Gather the decomposition information to the rank 0 SGROUP processes
          call xcomm_reduce(XKGROUP,MPI_SUM,rsa_kg,rsa) ; if (error()) goto 100

          select case (nsg)
          case (1)

            ! Output the results
            if (i_access(f)) then
              kpt = x_kpoint(el%o%kpoints,ik)
              kwt = x_kweight(el%o%kpoints,ik)
              write(x_unit(f),'(/,t2,"k-point #",i0,":",3f9.5,"; weight = ",f7.5)') ik, kpt, kwt
              do ib = 1,nb
                bt = ib - 1 + b1
                eig = el%o%eigs(ik,bt)
                occ = el%o%occs(ik,bt)
                write(x_unit(f),'(/,t3,"band #",i0,": eigenvalue = ",f9.5,"; occupation = ",f7.5)') bt, eig, occ
                select case (mode)
                case ("L", "l")
                  write(x_unit(f),'(/,t4,"site",7x,"s",9x,"p",9x,"d",7x,"total")')
                  write(x_unit(f),'(t3,46("-"))')
                  sum_s = 0.0_double
                  sum_p = 0.0_double
                  sum_d = 0.0_double
                end select
                do is = 1,ns
                  select case (mode)
                  case ("L", "l")
                    write(x_unit(f),'(t3,i5,4f10.4)') is, rsa(1,ib,is), sum(rsa(2:4,ib,is)), sum(rsa(5:9,ib,is)), sum(rsa(:,ib,is))
                    sum_s = sum_s + rsa(1,ib,is)
                    sum_p = sum_p + sum(rsa(2:4,ib,is))
                    sum_d = sum_d + sum(rsa(5:9,ib,is))
                  case ("LM", "lm")
                    write(x_unit(f),'(/,t6,"site #",i0)') is
                    write(x_unit(f),'(t9,"l = 0, m =  0:",2x,f6.4)') rsa(1,ib,is)
                    write(x_unit(f),'(t9,"l = 1, m = -1:",4x,f6.4)') rsa(2,ib,is)
                    write(x_unit(f),'(t16,      "m =  0:",4x,f6.4)') rsa(3,ib,is)
                    write(x_unit(f),'(t16,      "m = +1:",4x,f6.4)') rsa(4,ib,is)
                    write(x_unit(f),'(t9,"l = 2, m = -2:",6x,f6.4)') rsa(5,ib,is)
                    write(x_unit(f),'(t16,      "m = -1:",6x,f6.4)') rsa(6,ib,is)
                    write(x_unit(f),'(t16,      "m =  0:",6x,f6.4)') rsa(7,ib,is)
                    write(x_unit(f),'(t16,      "m = +1:",6x,f6.4)') rsa(8,ib,is)
                    write(x_unit(f),'(t16,      "m = +2:",6x,f6.4)') rsa(9,ib,is)
                  case ("XYZ", "xyz")
                    write(x_unit(f),'(/,t6,"site #",i0)') is
                    write(x_unit(f),'(t9,"l = 0:",14x,f6.4)') rsa(1,ib,is)
                    write(x_unit(f),'(t9,"l = 1, x:",13x,f6.4)') rsa(2,ib,is)
                    write(x_unit(f),'(t16,      "y:",13x,f6.4)') rsa(3,ib,is)
                    write(x_unit(f),'(t16,      "z:",13x,f6.4)') rsa(4,ib,is)
                    write(x_unit(f),'(t9,"l = 2, xy:",14x,f6.4)') rsa(5,ib,is)
                    write(x_unit(f),'(t16,      "xz:",14x,f6.4)') rsa(6,ib,is)
                    write(x_unit(f),'(t16,      "yz:",14x,f6.4)') rsa(7,ib,is)
                    write(x_unit(f),'(t16,      "(x**2 - y**2):",3x,f6.4)') rsa(8,ib,is)
                    write(x_unit(f),'(t16,      "(3z**2 - r**2):",2x,f6.4)') rsa(9,ib,is)
                  end select
                end do
                select case (mode)
                case ("L", "l")
                  write(x_unit(f),'(t11,38("="))')
                  write(x_unit(f),'(t8,4f10.4)') sum_s, sum_p, sum_d, (sum_s + sum_p + sum_d)
                end select
                write(x_unit(f),'(" ")')
              end do
            end if

          case (2)

            ! Gather all decomposition results to the rank 0 CONFIG process
            rsa_sg = 0.0_double
            rsa_sg(msg,:,:,:) = rsa
            call xcomm_reduce(XSGROUP,MPI_SUM,rsa_sg,rsa_c) ; if (error()) goto 100

            ! Gather all eigenvalues and occupations to the rank 0 CONFIG process
            eigs_sg = 0.0_double
            eigs_sg(msg,:,:) = el%o%eigs
            call xcomm_reduce(XSGROUP,MPI_SUM,eigs_sg,eigs_c) ; if (error()) goto 100
            occs_sg = 0.0_double
            occs_sg(msg,:,:) = el%o%occs
            call xcomm_reduce(XSGROUP,MPI_SUM,occs_sg,occs_c) ; if (error()) goto 100

            ! Output the results
            if (i_access(f)) then
              kpt = x_kpoint(el%o%kpoints,ik)
              kwt = x_kweight(el%o%kpoints,ik)
              write(x_unit(f),'(/,t2,"k-point #",i0,":",3f9.5,"; weight = ",f7.5)') ik, kpt, kwt
              do ib = 1,nb
                bt = ib - 1 + b1
                eig1 = eigs_c(1,ik,bt)
                eig2 = eigs_c(2,ik,bt)
                occ1 = occs_c(1,ik,bt)
                occ2 = occs_c(2,ik,bt)
                write(x_unit(f),'(/,t3,"band #",i0,": eigenvalue = ",f9.5,"; occupation = ",f7.5, &
                                                & 6x,"eigenvalue = ",f9.5,"; occupation = ",f7.5)') bt, eig1, occ1, eig2, occ2
                select case (mode)
                case ("L", "l")
                  write(x_unit(f),'(/,t4,"site",7x,"s",9x,"p",9x,"d",7x,"total",17x,"s",9x,"p",9x,"d",7x,"total")')
                  write(x_unit(f),'(t3,96("-"))')
                  sum_s1 = 0.0_double
                  sum_s2 = 0.0_double
                  sum_p1 = 0.0_double
                  sum_p2 = 0.0_double
                  sum_d1 = 0.0_double
                  sum_d2 = 0.0_double
                end select
                do is = 1,ns
                  select case (mode)
                  case ("L", "l")
                    write(x_unit(f),'(t3,i5,4f10.4,10x,4f10.4)') is, &
                                     & rsa_c(1,1,ib,is), sum(rsa_c(1,2:4,ib,is)), sum(rsa_c(1,5:9,ib,is)), sum(rsa_c(1,:,ib,is)), &
                                     & rsa_c(2,1,ib,is), sum(rsa_c(2,2:4,ib,is)), sum(rsa_c(2,5:9,ib,is)), sum(rsa_c(2,:,ib,is))
                    sum_s1 = sum_s1 + rsa_c(1,1,ib,is)
                    sum_s2 = sum_s2 + rsa_c(2,1,ib,is)
                    sum_p1 = sum_p1 + sum(rsa_c(1,2:4,ib,is))
                    sum_p2 = sum_p2 + sum(rsa_c(2,2:4,ib,is))
                    sum_d1 = sum_d1 + sum(rsa_c(1,5:9,ib,is))
                    sum_d2 = sum_d2 + sum(rsa_c(2,5:9,ib,is))
                  case ("LM", "lm")
                    write(x_unit(f),'(/,t6,"site #",i0)') is
                    write(x_unit(f),'(t9,"l = 0, m =  0:",2x,f6.4,14x,f6.4)') rsa_c(1,1,ib,is), rsa_c(2,1,ib,is)
                    write(x_unit(f),'(t9,"l = 1, m = -1:",4x,f6.4,14x,f6.4)') rsa_c(1,2,ib,is), rsa_c(2,2,ib,is)
                    write(x_unit(f),'(t16,      "m =  0:",4x,f6.4,14x,f6.4)') rsa_c(1,3,ib,is), rsa_c(2,3,ib,is)
                    write(x_unit(f),'(t16,      "m = +1:",4x,f6.4,14x,f6.4)') rsa_c(1,4,ib,is), rsa_c(2,4,ib,is)
                    write(x_unit(f),'(t9,"l = 2, m = -2:",6x,f6.4,14x,f6.4)') rsa_c(1,5,ib,is), rsa_c(2,5,ib,is)
                    write(x_unit(f),'(t16,      "m = -1:",6x,f6.4,14x,f6.4)') rsa_c(1,6,ib,is), rsa_c(2,6,ib,is)
                    write(x_unit(f),'(t16,      "m =  0:",6x,f6.4,14x,f6.4)') rsa_c(1,7,ib,is), rsa_c(2,7,ib,is)
                    write(x_unit(f),'(t16,      "m = +1:",6x,f6.4,14x,f6.4)') rsa_c(1,8,ib,is), rsa_c(2,8,ib,is)
                    write(x_unit(f),'(t16,      "m = +2:",6x,f6.4,14x,f6.4)') rsa_c(1,9,ib,is), rsa_c(2,9,ib,is)
                  case ("XYZ", "xyz")
                    write(x_unit(f),'(/,t6,"site #",i0)') is
                    write(x_unit(f),'(t9,"l = 0:",14x,f6.4,14x,f6.4)')                rsa_c(1,1,ib,is), rsa_c(2,1,ib,is)
                    write(x_unit(f),'(t9,"l = 1, x:",13x,f6.4,14x,f6.4)')             rsa_c(1,2,ib,is), rsa_c(2,2,ib,is)
                    write(x_unit(f),'(t16,      "y:",13x,f6.4,14x,f6.4)')             rsa_c(1,3,ib,is), rsa_c(2,3,ib,is)
                    write(x_unit(f),'(t16,      "z:",13x,f6.4,14x,f6.4)')             rsa_c(1,4,ib,is), rsa_c(2,4,ib,is)
                    write(x_unit(f),'(t9,"l = 2, xy:",14x,f6.4,14x,f6.4)')            rsa_c(1,5,ib,is), rsa_c(2,5,ib,is)
                    write(x_unit(f),'(t16,      "xz:",14x,f6.4,14x,f6.4)')            rsa_c(1,6,ib,is), rsa_c(2,6,ib,is)
                    write(x_unit(f),'(t16,      "yz:",14x,f6.4,14x,f6.4)')            rsa_c(1,7,ib,is), rsa_c(2,7,ib,is)
                    write(x_unit(f),'(t16,      "(x**2 - y**2):",3x,f6.4,14x,f6.4)')  rsa_c(1,8,ib,is), rsa_c(2,8,ib,is)
                    write(x_unit(f),'(t16,      "(3z**2 - r**2):",2x,f6.4,14x,f6.4)') rsa_c(1,9,ib,is), rsa_c(2,9,ib,is)
                  end select
                end do
                select case (mode)
                case ("L", "l")
                  write(x_unit(f),'(t11,88("="))')
                  write(x_unit(f),'(t8,4f10.4,10x,4f10.4)') sum_s1, sum_p1, sum_d1, (sum_s1 + sum_p1 + sum_d1), &
                                                          & sum_s2, sum_p2, sum_d2, (sum_s2 + sum_p2 + sum_d2)
                end select
                write(x_unit(f),'(" ")')
              end do
            end if

          end select

        end do

100     if (allocated( rsa ))    deallocate( rsa )
        if (allocated( rsa_kg )) deallocate( rsa_kg )
        if (allocated( eigs_c ))  deallocate( eigs_c )
        if (allocated( eigs_sg )) deallocate( eigs_sg )
        if (allocated( occs_c ))  deallocate( occs_c )
        if (allocated( occs_sg )) deallocate( occs_sg )
        if (allocated( rsa_sg ))  deallocate( rsa_sg )
        if (allocated( rsa_c ))   deallocate( rsa_c )

        call glean(thy(el))
        call glean(thy(f))

        if (error("Exit electrons_sc_mod::decompose_el")) continue

      end subroutine

      subroutine distribute_el(el)
!doc$ subroutine distribute(el)
        type(electrons_sc_obj) :: el
!       requires: mpi_nkgroups() > 1. Must be followed by a call to release_el with no intervening modification of el.
!       modifies: el%o%wf
!       effects: Distributes wf data among all kgroups.
!       errors: Passes errors.

!cod$
        integer :: ik, rank

        call my(el)

        do ik = 1,size(el%o%wf)
          rank = (el%o%kgroup_index(ik) - 1)
          call distribute(el%o%wf(ik),rank) ; if (error()) goto 100
        end do

100     call glean(thy(el))

        if (error("Exit electrons_sc_mod::distribute_el")) continue

      end subroutine

      subroutine release_el(el)
!doc$ subroutine release(el)
        type(electrons_sc_obj) :: el
!       requires: mpi_nkgroups() > 1. Must be preceeded by a call to distribute_el with no intervening modification of el.
!       modifies: el%o%wf
!       effects: Releases unneeded wf memory.
!       errors: Passes errors.

!cod$
        integer :: ik

        call my(el)

        do ik = 1,size(el%o%wf)
          call release(el%o%wf(ik)) ; if (error()) exit
        end do

        call glean(thy(el))

        if (error("Exit electrons_sc_mod::release_el")) continue

      end subroutine

      subroutine write_restart_el(el,nrestf)
!doc$ subroutine write_restart(el,nrestf)
        type(electrons_sc_obj) :: el
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes el restart information to nrestf.
!       errors: Passes errors.

!cod$
        integer :: ik
        integer(long) :: dsize, iosl, ndata, s4

        call my(el)
        call my(nrestf)

        ! Start the ELECTRONS block
        if (i_access(nrestf)) call startblock(nrestf,"ELECTRONS")

        ! Write the k-point information
        call write_restart(el%o%kpoints,nrestf)

        if (i_access(nrestf)) then

          ! Start the PARAMETERS block
          call startblock(nrestf,"PARAMETERS")

          ! Write the cutoff
          call writetag(nrestf,"CUTOFF")
          dsize = sizeof_double
          ndata = 1
          call writef(el%o%cutoff,dsize,ndata,x_tagfd(nrestf),iosl)

          ! Write the number of bands
          call writetag(nrestf,"NUMBER_OF_BANDS")
          s4 = size(el%o%eigs,2)
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          ! Write the charge state
          call writetag(nrestf,"CHARGE_STATE")
          dsize = sizeof_double
          ndata = 1
          call writef(el%o%charge_state,dsize,ndata,x_tagfd(nrestf),iosl)

          ! Write the occupation method
          call writetag(nrestf,"OCCUPATION_METHOD")
          s4 = el%o%occupation_method
          dsize = sizeof_long
          ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          select case(el%o%occupation_method)
          case (THERMAL)

            ! Write kT
            call writetag(nrestf,"KT")
            dsize = sizeof_double
            ndata = 1
            call writef(el%o%kt,dsize,ndata,x_tagfd(nrestf),iosl)

            ! Write the free energy status
            call writetag(nrestf,"FREE_ENERGY_STATUS")
            if (el%o%use_free_energy) then
              s4 = 1
            else
              s4 = 0
            end if
            dsize = sizeof_long
            ndata = 1
            call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

            select case (mpi_nsgroups())
            case (2)

              ! Write the polarization method
              call writetag(nrestf,"POLARIZATION_METHOD")
              s4 = el%o%polarization_method
              dsize = sizeof_long
              ndata = 1
              call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

              ! Write the spin polarization
              call writetag(nrestf,"SPIN_POLARIZATION")
              dsize = sizeof_double
              ndata = 1
              call writef(el%o%spin_polarization,dsize,ndata,x_tagfd(nrestf),iosl)

            end select

          end select

          ! Write the wavefunctions tolerance
          call writetag(nrestf,"WAVEFUNCTIONS_TOLERANCE")
          dsize = sizeof_double
          ndata = 1
          call writef(el%o%res_norm_tol,dsize,ndata,x_tagfd(nrestf),iosl)

          ! End the PARAMETERS block
          call endblock(nrestf)

        end if

        ! Write the wavefunctions information
        if (i_access(nrestf)) call writetag(nrestf,"WAVEFUNCTIONS_START")
        do ik = 1,size(el%o%wf)
          call write_restart(el%o%wf(ik),nrestf) ; if (error()) exit
        end do

        ! End the ELECTRONS block
        if (i_access(nrestf)) call endblock(nrestf)

        call glean(thy(el))
        call glean(thy(nrestf))

        if (error("Exit electrons_sc_mod::write_restart_el")) continue

      end subroutine

! private routines

      subroutine accumulate_density_i(elr,ext)
        type(electrons_sc_rep) :: elr
        type(external_obj) :: ext
!       requires: elr%density exist, have KGROUP scope, and be zeroed.
!       notes: elr%density enters with KGROUP scope and exits with CONFIG scope.

        integer :: ik
        real(double), dimension(:), allocatable :: wts

        call my(ext)

        allocate( wts(size(elr%occs,2)) )
        do ik = 1,size(elr%wf)
          wts = elr%occs(ik,:)*x_kweight(elr%kpoints,ik) ; if (error()) goto 100
          call add_density(elr%wf(ik),wts,elr%density)   ; if (error()) goto 100
        end do
        call merge_symmetrize_filter(elr%density,ext)

100     if (allocated( wts )) deallocate( wts )

        call glean(thy(ext))

        if (error("Exit electrons_sc_mod::accumulate_density_i")) continue

      end subroutine

      subroutine eigenvalues_i(elr)
        type(electrons_sc_rep) :: elr

        integer :: ik
        real(double), dimension(:), allocatable :: eigs_k
        real(double), dimension(:,:), allocatable :: eigs_kg

        allocate( eigs_k(size(elr%eigs,2)) )
        allocate( eigs_kg(size(elr%eigs,1),size(elr%eigs,2)) )
        do ik = 1,size(elr%wf)
          call eigenvalues(elr%wf(ik),eigs_k)
          eigs_kg(ik,:) = eigs_k
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,eigs_kg,elr%eigs)

        if (allocated( eigs_k)) deallocate( eigs_k )
        if (allocated( eigs_kg )) deallocate( eigs_kg )

        if (error("Exit electrons_sc_mod::eigenvalues_i")) continue

      end subroutine

      subroutine occupations_i(elr)
        type(electrons_sc_rep) :: elr

        real(double), parameter :: eftol = 1.0e-15_double
        real(double), parameter :: qtol = 1.0e-10_double
        logical, dimension(:), allocatable :: kmap
        integer :: ib, ik, ikt, it, msg, nb, nk, ns, nsg
        real(double) :: d, emax, emin, q, total_charge, sp, spin_degeneracy, x
        real(double) :: emin1, emin2, emax1, emax2, fsg1, fsg2, qsg1, qsg2
        real(double) :: rns, sum_occs
        real(double), dimension(:), allocatable :: kwts
        real(double), dimension(:,:), allocatable :: occs_us
        real(double), dimension(:,:,:), allocatable :: eigs_c, eigs_sg

        nk = size(elr%eigs,1)
        nb = size(elr%eigs,2)

        allocate( kwts(nk) )
        kwts = x_kweights(elr%kpoints)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()
        spin_degeneracy = 2.0_double/real(nsg,double)

        select case (elr%occupation_method)
        case (THERMAL)

          select case (nsg)
          case (1)  ! no spin

            ! Find the Fermi level
            emin = minval(elr%eigs)
            emax = maxval(elr%eigs)
            d = emax - emin
            emin = emin - 0.5_double*d
            emax = emax + 0.5_double*d
            it = 0
            do while ( (it < 1000) .and. ((emax - emin) > eftol) )
              it = it + 1
              elr%fermi_level = (emax + emin)/2.0_double
              q = 0.0_double
              do ik = 1,nk
                do ib = 1,nb
                  x = (elr%eigs(ik,ib) - elr%fermi_level)
                  q = q + spin_degeneracy*kwts(ik)*fermi_distr(x,elr%kt)
                end do
              end do
              if (q .in. nbhd(elr%total_charge,qtol)) exit
              if (q > elr%total_charge) then
                emax = elr%fermi_level
              elseif (q < elr%total_charge) then
                emin = elr%fermi_level
              end if
            end do
            if (error(.not.(q .in. nbhd(elr%total_charge,qtol)), "ERROR: Fermi level was not found")) then
              call notify_occupations_i(elr,it,emin,emax,elr%fermi_level,q,spin_degeneracy,kwts)
              goto 100
            end if

          case (2)  ! spin

            select case (elr%polarization_method)
            case (FIXED)

              ! Set the desired spin group charges
              select case (msg)
              case (1)
                total_charge = (elr%total_charge - elr%spin_polarization)/2.0_double
              case (2)
                total_charge = (elr%total_charge + elr%spin_polarization)/2.0_double
              end select

              ! Find the Fermi levels
              emin = minval(elr%eigs)
              emax = maxval(elr%eigs)
              d = emax - emin
              emin = emin - 0.5_double*d
              emax = emax + 0.5_double*d
              it = 0
              do while ( (it < 1000) .and. ((emax - emin) > eftol) )
                it = it + 1
                elr%fermi_level = (emax + emin)/2.0_double
                q = 0.0_double
                do ik = 1,nk
                  do ib = 1,nb
                    x = (elr%eigs(ik,ib) - elr%fermi_level)
                    q = q + spin_degeneracy*kwts(ik)*fermi_distr(x,elr%kt)
                  end do
                end do
                if (q .in. nbhd(total_charge,qtol)) exit
                if (q > total_charge) then
                  emax = elr%fermi_level
                elseif (q < total_charge) then
                  emin = elr%fermi_level
                end if
              end do
              if (error(.not.(q .in. nbhd(total_charge,qtol)), "ERROR: Fermi level was not found")) then
                call notify_occupations_i(elr,it,emin,emax,elr%fermi_level,q,spin_degeneracy,kwts)
              end if

              ! Synchronize the CONFIG error states
              call sync_config_process_errors() ; if (error()) goto 100

            case (VARIABLE)

              if (.not.associated(elr%occs)) then  ! Same as for FIXED polarization method

                ! Set the desired spin group charges
                select case (msg)
                case (1)
                  total_charge = (elr%total_charge - elr%spin_polarization)/2.0_double
                case (2)
                  total_charge = (elr%total_charge + elr%spin_polarization)/2.0_double
                end select

                ! Find the Fermi levels
                emin = minval(elr%eigs)
                emax = maxval(elr%eigs)
                d = emax - emin
                emin = emin - 0.5_double*d
                emax = emax + 0.5_double*d
                it = 0
                do while ( (it < 1000) .and. ((emax - emin) > eftol) )
                  it = it + 1
                  elr%fermi_level = (emax + emin)/2.0_double
                  q = 0.0_double
                  do ik = 1,nk
                    do ib = 1,nb
                      x = (elr%eigs(ik,ib) - elr%fermi_level)
                      q = q + spin_degeneracy*kwts(ik)*fermi_distr(x,elr%kt)
                    end do
                  end do
                  if (q .in. nbhd(total_charge,qtol)) exit
                  if (q > total_charge) then
                    emax = elr%fermi_level
                  elseif (q < total_charge) then
                    emin = elr%fermi_level
                  end if
                end do
                if (error(.not.(q .in. nbhd(total_charge,qtol)), "ERROR: Fermi level was not found")) then
                  call notify_occupations_i(elr,it,emin,emax,elr%fermi_level,q,spin_degeneracy,kwts)
                end if

                ! Synchronize CONFIG error states
                call sync_config_process_errors() ; if (error()) goto 100

              else

                ! Allgather the eigenvalues to the CONFIG processes
                allocate( eigs_c(nsg,nk,nb), eigs_sg(nsg,nk,nb) )
                eigs_sg = 0.0_double
                eigs_sg(msg,:,:) = elr%eigs
                call xcomm_allreduce(XSGROUP,MPI_SUM,eigs_sg,eigs_c) ; if (error()) goto 100
                deallocate( eigs_sg )

                ! Set a constraint on the spin polarization
                sp = elr%spin_polarization

                ! Find the Fermi levels
                emin1 = minval(eigs_c(1,:,:))
                emax1 = maxval(eigs_c(1,:,:))
                d = emax1 - emin1
                emin1 = emin1 - 0.5_double*d
                emax1 = emax1 + 0.5_double*d
                emin2 = minval(eigs_c(2,:,:))
                emax2 = maxval(eigs_c(2,:,:))
                d = emax2 - emin2
                emin2 = emin2 - 0.5_double*d
                emax2 = emax2 + 0.5_double*d
                it = 0
                do while ( (it < 1000) .and. ((emax1 - emin1) > eftol) .and. ((emax2 - emin2) > eftol) )
                  it = it + 1
                  fsg1 = (emax1 + emin1)/2.0_double
                  qsg1 = 0.0_double - sp/2.0_double
                  do ik = 1,nk
                    do ib = 1,nb
                      x = (eigs_c(1,ik,ib) - fsg1)
                      qsg1 = qsg1 + spin_degeneracy*kwts(ik)*fermi_distr(x,elr%kt)
                    end do
                  end do
                  fsg2 = (emax2 + emin2)/2.0_double
                  qsg2 = 0.0_double + sp/2.0_double
                  do ik = 1,nk
                    do ib = 1,nb
                      x = (eigs_c(2,ik,ib) - fsg2)
                      qsg2 = qsg2 + spin_degeneracy*kwts(ik)*fermi_distr(x,elr%kt)
                    end do
                  end do
                  q = qsg1 + qsg2
                  if (q .in. nbhd(elr%total_charge,qtol)) exit
                  if (q > elr%total_charge) then
                    emax1 = fsg1
                    emax2 = fsg2
                  elseif (q < elr%total_charge) then
                    emin1 = fsg1
                    emin2 = fsg2
                  end if
                  sp = 0.0_double  ! Release the constraint on the spin polarization
                end do
                if (error(.not.(q .in. nbhd(elr%total_charge,qtol)), "ERROR: Fermi level was not found")) then
                  call notify_occupations_i(elr,it,emin1,emax1,fsg1,q,spin_degeneracy,kwts)
                  goto 100
                end if
                deallocate( eigs_c )

                ! Set the Fermi level (and q for temporary notification purposes)
                select case (msg)
                case (1)
                  elr%fermi_level = fsg1
                  q = qsg1
                case (2)
                  elr%fermi_level = fsg2
                  q = qsg2
                end select

                ! Reset the spin polarization
                elr%spin_polarization = qsg2 - qsg1

              end if

            end select  ! polarization method

          end select  ! nsg

          ! Compute the unsymmetrized occupations
          allocate( occs_us(nk,nb) )
          do ik = 1,nk
            do ib = 1,nb
              x = elr%eigs(ik,ib) - elr%fermi_level
              occs_us(ik,ib) = spin_degeneracy*fermi_distr(x,elr%kt)
            end do
          end do

          ! Allocate space for the occupations
          if (.not.associated(elr%occs)) allocate( elr%occs(nk,nb) )

          ! Symmetrize the occupations
          allocate( kmap(nk) )
          do ik = 1,nk
            kmap = x_kmap(elr%kpoints,ik) ; if (error()) goto 100
            ns = 0
            do ikt = 1,nk
              if (kmap(ikt)) ns = ns + 1
            end do
            rns = real(ns,double)
            do ib = 1,nb
              sum_occs = 0.0_double
              do ikt = 1,nk
                if (kmap(ikt)) sum_occs = sum_occs + occs_us(ikt,ib)
              end do
              elr%occs(ik,ib) = sum_occs/rns
            end do
          end do

          deallocate( occs_us )
          deallocate( kmap )

        case (UNIFORM)

          ! Allocate space for the occupations
          if (.not.associated(elr%occs)) allocate( elr%occs(nk,nb) )

          ! Compute the occupations
          q = elr%total_charge
          do ib = 1,nb
            if (q >= spin_degeneracy) then
              elr%occs(:,ib) = spin_degeneracy
              q = q - spin_degeneracy
            else
              elr%occs(:,ib) = q
              q = 0.0_double
            end if
          end do

          ! Set the Fermi level to the highest eigenvalue in the highest occupied band
          ib = 0                             ! (used in decomposing the wavefunctions) 
          do
            ib = ib + 1
            if (elr%occs(1,ib) == 0.0_double) exit
          end do
          elr%fermi_level = maxval(elr%eigs(:,ib-1))

        end select

100     if (allocated( kwts )) deallocate( kwts )
        if (allocated( eigs_c )) deallocate( eigs_c )
        if (allocated( eigs_sg )) deallocate( eigs_sg )

        if (error("Exit electrons_sc_mod::occupations_i")) continue

      end subroutine

      subroutine notify_occupations_i(elr,it,emin,emax,fl,q,sd,kwts)
        type(electrons_sc_rep) :: elr
        integer :: it
        real(double) :: emin, emax, fl, q, sd
        real(double), dimension(:) :: kwts

        integer :: ib, ik
        real(double) :: x

        call warn(" ")
        call notify("it",it)
        call notify("emin",emin)
        call notify("emax",emax)
        call notify("fermi_level",fl)
        call notify("total_charge",elr%total_charge)
        call notify("q",q)
        q = 0.0_double
        do ik = 1,size(elr%eigs,1)
          call warn(" ")
          call notify("    ik",ik)
          do ib = 1,size(elr%eigs,2)
            x = (elr%eigs(ik,ib) - fl)
            q = q + sd*kwts(ik)*fermi_distr(x,elr%kt)
            call notify("    e",elr%eigs(ik,ib))
            call notify("    o",sd*fermi_distr(x,elr%kt))
            call notify("    q",q)
          end do
        end do

      end subroutine

      subroutine residual_norm_i(elr)
        type(electrons_sc_rep) :: elr

        integer :: ik
        real(double) :: rn_k, rn_kg, rn_sg

        rn_kg = 0.0_double
        do ik = 1,size(elr%wf)
          call residual_norm(elr%wf(ik),rn_k)
          rn_kg = rn_kg + rn_k*x_kweight(elr%kpoints,ik)
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,rn_kg,rn_sg)        ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_SUM,rn_sg,elr%res_norm) ; if (error()) goto 100

        elr%res_norm_cvg = (elr%res_norm < elr%res_norm_tol)

100     if (error("Exit electrons_sc_mod::residual_norm_i")) continue

      end subroutine

      subroutine kinetic_energy_i(elr)
        type(electrons_sc_rep) :: elr

        integer :: ik
        real(double) :: ke_kg, ke_sg
        real(double), dimension(:), allocatable :: wts
        
        allocate( wts(size(elr%occs,2)) )
        ke_kg = 0.0_double
        do ik = 1,size(elr%wf)
          wts = elr%occs(ik,:)*x_kweight(elr%kpoints,ik)
          ke_kg = ke_kg + kinetic_energy(elr%wf(ik),wts)
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,ke_kg,ke_sg) ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_SUM,ke_sg,elr%kinetic_energy)

100     if (allocated( wts )) deallocate( wts )

        if (error("Exit electrons_sc_mod::kinetic_energy_i")) continue

      end subroutine

      subroutine mermin_energy_i(elr)
        type(electrons_sc_rep) :: elr

        integer :: ib, ik
        real(double) :: me_sg, spin_degeneracy, x

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        me_sg = 0.0_double
        do ik = 1,size(elr%eigs,1)
          do ib = 1,size(elr%eigs,2)
            x = elr%eigs(ik,ib) - elr%fermi_level
            me_sg = me_sg - spin_degeneracy*elr%kt*fermi_entropy(x,elr%kt)*x_kweight(elr%kpoints,ik)
          end do
        end do
        call xcomm_allreduce(XSGROUP,MPI_SUM,me_sg,elr%mermin_energy)

        if (error("Exit electrons_sc_mod::mermin_energy_i")) continue

      end subroutine

      subroutine diary_construction_i(elr,restf)
        type(electrons_sc_rep) :: elr
        type(tagio_obj), optional :: restf

        integer :: ik, nk, nkg, nsg

        if (present(restf)) call my(restf)

        nsg = mpi_nsgroups()

        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,"Electrons object construction:")')
          if (present(restf)) then
            write(x_unit(diaryfile()),'(/,t4,"Restart-file initialization")')
          end if
          write(x_unit(diaryfile()),'(/,t4,"Plane wave cutoff energy = ",f0.2," Ryd")') elr%cutoff
          write(x_unit(diaryfile()),'(/,t4,"Number of bands = ",i0)') size(elr%eigs,2)
        end if
        call diary(elr%kpoints)
        nk = x_n_kpoints(elr%kpoints)
        nkg = mpi_nkgroups()
        if (nkg == 1) then
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Plane wave expansion:")')
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t6,"k-point     plane waves",/)')
          do ik = 1,nk
            if (i_access(diaryfile())) write(x_unit(diaryfile()),'(t6,i5,9x,i7)') ik, x_n_gvectors(elr%wf(ik))
          end do
        else
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Plane wave expansion:")')
          if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t6,"k-point     kgroup     plane waves",/)')
          if (i_access(diaryfile())) then
            do ik = 1,nk
              write(x_unit(diaryfile()),'(t6,i5,8x,i4,8x,i7)') ik, (mod(ik-1,nkg) + 1), x_n_gvectors(elr%wf(ik))
            end do
          end if
        end if
        if (i_access(diaryfile())) then
          select case (elr%occupation_method)
          case (THERMAL)
            write(x_unit(diaryfile()),'(/,t4,"Fermi occupations with kT = ",es9.3," Ryd")') elr%kt
            if (elr%use_free_energy) write(x_unit(diaryfile()),'(/,t4,"Using Mermin free-energy")')
          case (UNIFORM)
            write(x_unit(diaryfile()),'(/,t4,"Uniform occupations")')
          end select
        end if
        if (i_access(diaryfile())) then
          select case (nsg)
          case (2)
            select case (elr%polarization_method)
            case (FIXED)
              write(x_unit(diaryfile()),'(/,t4,"Using fixed spin polarization = ",f0.3)') elr%spin_polarization
            case (VARIABLE)
              write(x_unit(diaryfile()),'(/,t4,"Using variable spin polarization with starting value = ",f0.3)') elr%spin_polarization
            end select
          end select
        end if
        call diary(elr%hc)
        do ik = 1,nk
          if (mpi_mykgroup() /= elr%kgroup_index(ik)) cycle
          call diary(elr%wf(ik))
          exit
        end do

        if (present(restf)) call glean(thy(restf))

        if (error("Exit electrons_sc_mod::diary_construction_i")) continue

      end subroutine

      subroutine own_i(el)
        type(electrons_sc_obj) :: el
        type(electrons_sc_obj) :: elt
        integer :: ik
        if (el%ref < el%o%ref) then
          allocate( elt%o )
          elt%o%ref = 0
          elt%o%g = el%o%g
          elt%o%g_external = el%o%g_external
          elt%o%use_free_energy = el%o%use_free_energy
          elt%o%res_norm_cvg = el%o%res_norm_cvg
          elt%o%occupation_method = el%o%occupation_method
          elt%o%polarization_method = el%o%polarization_method
          elt%o%res_norm = el%o%res_norm
          elt%o%res_norm_tol = el%o%res_norm_tol
          elt%o%total_charge = el%o%total_charge
          elt%o%charge_state = el%o%charge_state
          elt%o%spin_polarization = el%o%spin_polarization
          elt%o%fermi_level = el%o%fermi_level
          elt%o%kt = el%o%kt
          elt%o%cutoff = el%o%cutoff
          elt%o%kinetic_energy = el%o%kinetic_energy
          elt%o%mermin_energy = el%o%mermin_energy
          allocate( elt%o%kgroup_index(size(el%o%kgroup_index)) )
          elt%o%kgroup_index = el%o%kgroup_index
          allocate( elt%o%eigs(size(el%o%eigs,1),size(el%o%eigs,2)) )
          elt%o%eigs = el%o%eigs
          allocate( elt%o%occs(size(el%o%occs,1),size(el%o%occs,2)) )
          elt%o%occs = el%o%occs
          if (associated(el%o%oep_gksc)) then
             allocate( elt%o%oep_gksc(size(el%o%oep_gksc,1),size(el%o%oep_gksc,2)) )
             elt%o%oep_gksc = el%o%oep_gksc
          else
             nullify(elt%o%oep_gksc)
          end if
          call my(el%o%kpoints,elt%o%kpoints)
          call my(el%o%density,elt%o%density)
          call my(el%o%hc,elt%o%hc)
          allocate( elt%o%wf(size(el%o%wf)) )
          do ik = 1,size(elt%o%wf)
            call my(el%o%wf(ik),elt%o%wf(ik))
          end do
          el%o%ref = el%o%ref - el%ref
          el%o => elt%o
          el%o%ref = el%o%ref + el%ref
        end if
      end subroutine



      subroutine get_wf_filename_i(fname,ftype,ik,nk,ib,nb,is,it)
        character(line_len) :: fname
        integer :: ftype
        integer :: ik
        integer :: nk
        integer :: ib
        integer :: nb
        integer, optional :: is   ! Spin index
        integer, optional :: it   ! Time step

        ! Local Vars
        integer            :: pos
        integer            :: j
        integer            :: a
        integer            :: kpt_offset
        integer            :: band_offset
        integer            :: spin_offset
        integer            :: time_offset
        character(line_len) :: it_char
        character(line_len) :: ik_char
        character(line_len) :: ib_char
        character(4)       :: kpt_fmt
        character(4)       :: band_fmt


        !** calculate offsets
        if (error(nk<0,"Error:get_wf_filename_i - invalid nk, nk<0")) goto 100
        if (nk < 10) then
           kpt_offset = 1
           kpt_fmt = "(I1)"
        elseif (nk < 100) then
           kpt_offset = 2
           kpt_fmt = "(I2)"
        elseif (nk < 1000) then
           kpt_offset = 3
           kpt_fmt = "(I3)"
        elseif (nk < 10000) then
           kpt_offset = 4
           kpt_fmt = "(I4)"
        elseif (nk < 100000) then
           kpt_offset = 5
           kpt_fmt = "(I5)"
        elseif (nk < 1000000) then
           kpt_offset = 6
           kpt_fmt = "(I6)"
        elseif (nk < 10000000) then
           kpt_offset = 7
           kpt_fmt = "(I7)"
        elseif (nk < 100000000) then
           kpt_offset = 8
           kpt_fmt = "(I8)"
        else
           if (error(.true.,"Error:get_wf_filename_i - invalid nk, nk too large")) goto 100
        end if

        if (error(nb<0,"Error:get_wf_filename_i - invalid nb, nb<0")) goto 100
        if (nb < 10) then
           band_offset = 1
           band_fmt = "(I1)"
        elseif (nb < 100) then
           band_offset = 2
           band_fmt = "(I2)"
        elseif (nb < 1000) then
           band_offset = 3
           band_fmt = "(I3)"
        elseif (nb < 10000) then
           band_offset = 4
           band_fmt = "(I4)"
        elseif (nb < 100000) then
           band_offset = 5
           band_fmt = "(I5)"
        elseif (nb < 1000000) then
           band_offset = 6
           band_fmt = "(I6)"
        elseif (nb < 10000000) then
           band_offset = 7
           band_fmt = "(I7)"
        elseif (nb < 100000000) then
           band_offset = 8
           band_fmt = "(I8)"
        else
           if (error(.true.,"Error:get_wf_filename_i - invalid nb, nb too large")) goto 100
        end if


        !** Initialize the filename
        fname = "wfn_"
        pos = 5


        
        !** Write the time step if present
        if (present(it)) then
          if (error(it<0,"Error:get_wf_filename_i - invalid time step, it<0")) goto 100
          !   Convert it which is an integer to a character string, it_char
          fname(pos:(pos+2)) = "it_"
          pos = pos + 3
          write(it_char,"(I8)") it
          do j=1,8
             a = iachar(it_char(j:j) )
             if (a == 32) then
                fname(pos:pos) = "0"
             else
                fname(pos:pos) = it_char(j:j)
             end if
             pos = pos + 1
          end do
          fname(pos:pos) = "_"
          pos = pos + 1
        end if

        !** Write the spin, if present
        if (present(is)) then
          select case(is)
          case (1)
            fname(pos:(pos+2)) = "up_"
          case (2)
            fname(pos:(pos+2)) = "dn_"
          case default
             if (error(.true.,'Error:get_wf_filename_i - invalid spin index')) goto 100
          end select
 
          pos = pos + 3
        end if

        
        !** Write the k point
        !   Convert ik which is an integer to a character string, ik_char
        fname(pos:(pos+2)) = "ik_"
        pos = pos + 3
        write(ik_char,kpt_fmt) ik
        do j=1,kpt_offset
           a = iachar(ik_char(j:j) )
           if (a == 32) then
              fname(pos:pos) = "0"
           else
              fname(pos:pos) = ik_char(j:j)
           end if
           pos = pos + 1
        end do
        fname(pos:pos) = "_"
        pos = pos + 1



        !** Write the band index
        !   Convert ib which is an integer to a character string, ib_char
        fname(pos:(pos+2)) = "ib_"
        pos = pos + 3
        write(ib_char,band_fmt) ib
        do j=1,band_offset
           a = iachar(ib_char(j:j) )
           if (a == 32) then
              fname(pos:pos) = "0"
           else
              fname(pos:pos) = ib_char(j:j)
           end if
           pos = pos + 1
        end do


        !** Write the suffix
        select case (ftype)
        case ( MATLAB )
           fname(pos:(pos+3)) = ".mat"
        case ( AMIRA )
           fname(pos:(pos+2)) = ".am"
        case ( VTK )
           fname(pos:(pos+3)) = ".vtk"
        case default
           if (error(.true.,"Error! Unrecognized file type")) goto 100
        end select


100     if (error("electrons_mod::get_wf_filename_i - Exiting")) continue

      end subroutine 


      subroutine get_sd_filename_i(fname,ftype,ik,nk,ib,nb,is,it)
        character(line_len) :: fname
        integer :: ftype
        integer :: ik
        integer :: nk
        integer :: ib
        integer :: nb
        integer, optional :: is   ! Spin index
        integer, optional :: it   ! Time step

        ! Local Vars
        integer            :: pos
        integer            :: j
        integer            :: a
        integer            :: kpt_offset
        integer            :: band_offset
        integer            :: spin_offset
        integer            :: time_offset
        character(line_len) :: it_char
        character(line_len) :: ik_char
        character(line_len) :: ib_char
        character(4)       :: kpt_fmt
        character(4)       :: band_fmt


        !** calculate offsets
        if (error(nk<0,"Error:get_sd_filename_i - invalid nk, nk<0")) goto 100
        if (nk < 10) then
           kpt_offset = 1
           kpt_fmt = "(I1)"
        elseif (nk < 100) then
           kpt_offset = 2
           kpt_fmt = "(I2)"
        elseif (nk < 1000) then
           kpt_offset = 3
           kpt_fmt = "(I3)"
        elseif (nk < 10000) then
           kpt_offset = 4
           kpt_fmt = "(I4)"
        elseif (nk < 100000) then
           kpt_offset = 5
           kpt_fmt = "(I5)"
        elseif (nk < 1000000) then
           kpt_offset = 6
           kpt_fmt = "(I6)"
        elseif (nk < 10000000) then
           kpt_offset = 7
           kpt_fmt = "(I7)"
        elseif (nk < 100000000) then
           kpt_offset = 8
           kpt_fmt = "(I8)"
        else
           if (error(.true.,"Error:get_sd_filename_i - invalid nk, nk too large")) goto 100
        end if

        if (error(nb<0,"Error:get_sd_filename_i - invalid nb, nb<0")) goto 100
        if (nb < 10) then
           band_offset = 1
           band_fmt = "(I1)"
        elseif (nb < 100) then
           band_offset = 2
           band_fmt = "(I2)"
        elseif (nb < 1000) then
           band_offset = 3
           band_fmt = "(I3)"
        elseif (nb < 10000) then
           band_offset = 4
           band_fmt = "(I4)"
        elseif (nb < 100000) then
           band_offset = 5
           band_fmt = "(I5)"
        elseif (nb < 1000000) then
           band_offset = 6
           band_fmt = "(I6)"
        elseif (nb < 10000000) then
           band_offset = 7
           band_fmt = "(I7)"
        elseif (nb < 100000000) then
           band_offset = 8
           band_fmt = "(I8)"
        else
           if (error(.true.,"Error:get_sd_filename_i - invalid nb, nb too large")) goto 100
        end if


        !** Initialize the filename
        fname = "stateden_"
        pos = 10


        
        !** Write the time step if present
        if (present(it)) then
          if (error(it<0,"Error:get_sd_filename_i - invalid time step, it<0")) goto 100
          !   Convert it which is an integer to a character string, it_char
          fname(pos:(pos+2)) = "it_"
          pos = pos + 3
          write(it_char,"(I8)") it
          do j=1,8
             a = iachar(it_char(j:j) )
             if (a == 32) then
                fname(pos:pos) = "0"
             else
                fname(pos:pos) = it_char(j:j)
             end if
             pos = pos + 1
          end do
          fname(pos:pos) = "_"
          pos = pos + 1
        end if

        !** Write the spin, if present
        if (present(is)) then
          select case(is)
          case (1)
            fname(pos:(pos+2)) = "up_"
          case (2)
            fname(pos:(pos+2)) = "dn_"
          case default
             if (error(.true.,'Error:get_sd_filename_i - invalid spin index')) goto 100
          end select
 
          pos = pos + 3
        end if

        
        !** Write the k point
        !   Convert ik which is an integer to a character string, ik_char
        fname(pos:(pos+2)) = "ik_"
        pos = pos + 3
        write(ik_char,kpt_fmt) ik
        do j=1,kpt_offset
           a = iachar(ik_char(j:j) )
           if (a == 32) then
              fname(pos:pos) = "0"
           else
              fname(pos:pos) = ik_char(j:j)
           end if
           pos = pos + 1
        end do
        fname(pos:pos) = "_"
        pos = pos + 1



        !** Write the band index
        !   Convert ib which is an integer to a character string, ib_char
        fname(pos:(pos+2)) = "ib_"
        pos = pos + 3
        write(ib_char,band_fmt) ib
        do j=1,band_offset
           a = iachar(ib_char(j:j) )
           if (a == 32) then
              fname(pos:pos) = "0"
           else
              fname(pos:pos) = ib_char(j:j)
           end if
           pos = pos + 1
        end do


        !** Write the suffix
        select case (ftype)
        case ( MATLAB )
           fname(pos:(pos+3)) = ".mat"
        case ( AMIRA )
           fname(pos:(pos+2)) = ".am"
        case ( VTK )
           fname(pos:(pos+3)) = ".vtk"
        case default
           if (error(.true.,"Error! Unrecognized file type")) goto 100
        end select


100     if (error("electrons_mod::get_sd_filename_i - Exiting")) continue

      end subroutine 







      end module

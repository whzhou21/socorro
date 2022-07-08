!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module electrons_fh_mod
!doc$ module electrons_fh_mod

!     One datatype is available here: type(electrons_fh_obj)

!     electrons_fh_mod creates and maintains a set of electrons in calculations with a fixed hamiltonian.

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
      use kpoints_mod
      use wavefunctions_es_mod
      use operators_mod
      use multibasis_mod
      use multivector_mod
      use gen_potential_mod
      use timing_mod

!cod$
      implicit none
      private

      ! occupation method
      integer, parameter :: THERMAL = 1
      integer, parameter :: UNIFORM = 2

      type :: electrons_fh_rep
        integer :: ref
        type(ghost) :: g
        logical :: res_norm_cvg                                    ! convergence status
        integer :: occupation_method                               ! occupation method
        real(double) :: res_norm                                   ! wavefunctions residual norm
        real(double) :: res_norm_tol                               ! tolerance for wavefunctions residual norm
        real(double) :: total_charge                               ! total electronic charge in supercell
        real(double) :: spin_polarization                          ! spin polarization
        real(double) :: charge_state                               ! charge state of supercell
        real(double) :: fermi_level                                ! Fermi energy
        real(double) :: kt                                         ! kT
        real(double) :: cutoff                                     ! wavefunctions cutoff energy
        integer, dimension(:), pointer :: kgroup_index             ! mapping of k-points to kgroups
        real(double), dimension(:,:), pointer :: eigs              ! eigenvalues
        real(double), dimension(:,:), pointer :: occs              ! occupations
        type(kpoints_obj) :: kpoints                               ! k-points object
        type(h_common_obj) :: hc                                   ! common hamiltonian object
        type(wavefunctions_es_obj), dimension(:), pointer :: wf    ! wavefunction objects
      end type

      type, public :: electrons_fh_obj
        private
        integer :: ref
        type(electrons_fh_rep), pointer :: o
      end type

!doc$
      public :: electrons_fh
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_wavefunctions
      public :: x_n_bands
      public :: x_eigenvalue
      public :: x_eigenvalues
      public :: x_residual_norm
      public :: x_converged
      public :: diary
      public :: decompose

!cod$
      interface electrons_fh
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
      interface x_ref
        module procedure el_ref
      end interface
      interface x_ghost
        module procedure el_ghost
      end interface
      interface x_wavefunctions
         module procedure el_wavefunctions
      end interface
      interface x_n_bands
         module procedure el_n_bands
      end interface
      interface x_eigenvalue
         module procedure el_eigenvalue
      end interface
      interface x_eigenvalues
         module procedure el_eigenvalues
      end interface
      interface x_residual_norm
         module procedure el_residual_norm
      end interface
      interface x_converged
        module procedure el_converged
      end interface
      interface diary
        module procedure diary_el
      end interface
      interface decompose
        module procedure decompose_el
      end interface

      contains

! public routines

      function constructor_el(ext,gp,restf) result(el)
!doc$ function electrons_fh(ext,gp,restf) result(el)
        type(external_obj) :: ext
        type(gen_potential_obj) :: gp
        type(tagio_obj), optional :: restf
        type(electrons_fh_obj) :: el
!       requires: Consistent ext and gp.
!       effects: Constructs a new el.
!       errors: Passes errors.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: tag, usage
        integer :: ik, nb, nk, nsg, r_nb
        integer(long) :: dsize, iosl, ndata, s4
        integer, dimension(2) :: csr
        real(double), dimension(3) :: kpt
        type(layout_obj) :: lay

        if (error("  Error on entry")) then
          el%ref = 0
          allocate( el%o )
          el%o%ref = 0
          goto 999
        end if

        call start_timer("electrons_fh: constructor")

        call my(ext)
        call my(gp)
        if (present(restf)) call my(restf)

        el%ref = 0
        allocate( el%o )
        el%o%ref = 0
        el%o%g = x_ghost()

        call my(x_layout(ext),lay)

        nsg = mpi_nsgroups()

        if (present(restf)) then

          ! Open the ELECTRONS block
          if (i_access(restf)) tios = findfirsttag(restf,"ELECTRONS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: ELECTRONS block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)

          ! Read the k-points
          call my(kpoints(ext,restf),el%o%kpoints) ; if (error()) goto 200
          nk = x_n_kpoints(el%o%kpoints)

          ! Open the PARAMETERS block
          if (i_access(restf)) tios = findfirsttag(restf,"PARAMETERS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: PARAMETERS block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)

          ! Read the wavefunctions cutoff
          if (i_access(restf)) tios = findfirsttag(restf,"CUTOFF")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: CUTOFF tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
            call readf(el%o%cutoff,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%cutoff)

          ! Read the number of bands
          call arg("nbands",nb,found)
          if (error(.not.found,"ERROR: nbands tag was not found")) goto 100
          if (error(nb <= 0,"ERROR: nbands <= 0")) goto 100

          ! Check that the number of bands is greater than the number of bands in the restart file
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_BANDS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_BANDS tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long
            ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            r_nb = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_nb)
          if (error(r_nb <= nb,"ERROR: nb <= the number of bands in the restart file")) goto 100

          ! Read the charge state
          if (i_access(restf)) tios = findfirsttag(restf,"CHARGE_STATE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: CHARGE_STATE tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
            call readf(el%o%charge_state,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%charge_state)

          ! Compute the total number of electrons
          el%o%total_charge = valence_electrons(ext) - el%o%charge_state
          if (error(el%o%total_charge <= 0.0_double,"ERROR: total charge <= 0")) goto 100

          ! Read the occupation method
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

          select case(el%o%occupation_method)
          case (THERMAL)

            ! Read kT
            if (i_access(restf)) tios = findfirsttag(restf,"KT")
            if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
            if (error(tios == TAG_NOT_FOUND,"ERROR: KT tag was not found")) goto 100
            if (i_access(restf)) then
              dsize = sizeof_double
              ndata = 1
              call readf(el%o%kt,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%kt)

            select case (nsg)
            case (1)

              ! Set null values for irrelevant elements
              el%o%spin_polarization = 0.0_double

            case (2)

              ! Read the spin polarization
              if (i_access(restf)) tios = findfirsttag(restf,"SPIN_POLARIZATION")
              if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
              if (error(tios == TAG_NOT_FOUND,"ERROR: SPIN_POLARIZATION tag was not found")) goto 100
              if (i_access(restf)) then
                dsize = sizeof_double
                ndata = 1
                call readf(el%o%spin_polarization,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,el%o%spin_polarization)

            end select

          case (UNIFORM)

            ! Set null values for irrelevant elements
            el%o%spin_polarization = 0.0_double

          end select

          ! Read the wavefunctions tolerance
          call arg("wavefunctions_tolerance",el%o%res_norm_tol,found)
          if (.not.found) el%o%res_norm_tol = 1.0e-7_double  ! default value
          if (error(el%o%res_norm_tol < 0.0_double,"ERROR: wavefunctions_tolerance < 0")) continue

          ! Close the PARAMETERS block
100       if (i_access(restf)) call closeblock(restf)
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

          ! Form the wavefunctions
          allocate( el%o%wf(nk) )
          do ik = 1,nk
            kpt = x_kpoint(el%o%kpoints,ik) ; if (error()) exit
            if (mpi_mykgroup() == el%o%kgroup_index(ik)) then
              usage = "normal"
            else
              usage = "auxiliary"
            end if
            call my(wavefunctions_es(el%o%hc,multibasis(usage,el%o%cutoff,nb,kpt,lay),restf=restf),el%o%wf(ik)) ; if (error()) exit
          end do

          call sync_config_process_errors() ; if (error()) goto 200

          ! Eigenvalues and occupations
          allocate( el%o%eigs(nk,nb), el%o%occs(nk,nb) )
          call eigenvalues_i(el%o) ; if (error()) goto 200
          call occupations_i(el%o) ; if (error()) goto 200

          ! Residual norm
          call residual_norm_i(el%o)

          call diary_construction_i(el%o,restf)

          ! Close the ELECTRONS block
200       if (i_access(restf)) call closeblock(restf)

        else

          ! Generate the k-points
          call my(kpoints(ext),el%o%kpoints) ; if (error()) goto 300
          nk = x_n_kpoints(el%o%kpoints)

          ! Read the wavefunctions cutoff: should be the same as what was used to generate the restart file
          call arg("wf_cutoff",el%o%cutoff,found)
          if (error(.not.found,"ERROR: wf_cutoff was not found")) goto 300
          if (error(el%o%cutoff <= 0.0_double,"ERROR: wf_cutoff <= 0")) goto 300
          if (error(el%o%cutoff > 0.5_double*x_cutoff(lay),"ERROR: wf_cutoff is too large")) goto 300

          ! Read the number of bands
          call arg("nbands",nb,found)
          if (error(.not.found,"ERROR: nbands tag was not found")) goto 300
          if (error(nb <= 0,"ERROR: nbands <= 0")) goto 300

          ! Read the charge state: should be the same as what was used to generate the restart file
          call arglc("charge_state_mode",tag,found)
          if (.not.found) tag = "real_number"
          select case (trim(tag))
          case ("real_number")
            call arg("charge_state",el%o%charge_state,found)
            if (.not.found) el%o%charge_state = 0.0_double  ! default value
          case ("integer_ratio")
            call arg("charge_state_ratio",csr,found)
            if (error(.not.found,"ERROR: charge_state_ratio was not found")) goto 300
            if (error(csr(2) == 0,"ERROR: denominator = 0")) goto 300
            el%o%charge_state = real(csr(1),double)/real(csr(2),double)
          case default
            if (error(.true.,"ERROR: charge_state_mode was not recognized")) goto 300
          end select

          ! Compute the total number of electrons
          el%o%total_charge = valence_electrons(ext) - el%o%charge_state

          ! Check the number of bands
          if (error(nb < el%o%total_charge/2.0_double,"ERROR: not enough bands")) goto 300

          ! Read the occupation method: should be the same as what was used to generate the restart file
          call arglc("occupation_method",tag,found)
          if (.not.found) tag = "thermal"  ! default value
          select case (trim(tag))
          case ("thermal")
            el%o%occupation_method = THERMAL
          case ("uniform")
            el%o%occupation_method = UNIFORM
            if (error(nsg == 2,"ERROR: UNIFORM occupation_method is not compatible with spin")) goto 300
          case default
            if (error(.true.,"ERROR: occupation_method was not recognized")) goto 300
          end select

          select case(el%o%occupation_method)
          case (THERMAL)

            ! Read kT: should be the same as what was used to generate the restart file
            call arg("kt",el%o%kt,found)
            if (.not.found) el%o%kt = 5.0e-3_double  ! default value
            if (error(el%o%kt <= 0.0_double,"ERROR: kT <= 0.0")) goto 300

            select case (nsg)
            case (1)

              ! Set null values for irrelevant elements
              el%o%spin_polarization = 0.0_double

            case (2)

              ! Read the spin polarization: should be the same as what was used to generate the restart file
              call arg("spin_polarization",el%o%spin_polarization,found)
              if (.not.found) el%o%spin_polarization = 0.0_double  ! default value

            end select

          case (UNIFORM)

            ! Set null values for irrelevant elements
            el%o%spin_polarization = 0.0_double

          end select

          ! Read the wavefunctions tolerance
          call arg("wavefunctions_tolerance",el%o%res_norm_tol,found)
          if (.not.found) el%o%res_norm_tol = 1.0e-7_double  ! default value
          if (error(el%o%res_norm_tol < 0.0_double,"ERROR: wavefunctions_tolerance < 0")) goto 300

          ! Divide the k-points among processes
          if (error(nk < mpi_nkgroups(),"ERROR: nk is less than kgroups")) then
            call notify("Number of k-points = ",nk)
            goto 300
          end if
          if (mod(nk,mpi_nkgroups()) /= 0) then
            call warn("WARNING: non-equal division of k-points among kgroups")
          end if
          allocate( el%o%kgroup_index(nk) )
          do ik = 1,nk
            el%o%kgroup_index(ik) = mod(ik-1,mpi_nkgroups()) + 1
          end do

          ! Form the common hamiltonian
          call my(common_hamiltonian(ext,gp,el%o%cutoff,nb),el%o%hc) ; if (error()) goto 300

          ! Form the wavefunctions
          allocate( el%o%wf(nk) )
          do ik = 1,nk
            kpt = x_kpoint(el%o%kpoints,ik) ; if (error()) exit
            if (mpi_mykgroup() == el%o%kgroup_index(ik)) then
              usage = "normal"
            else
              usage = "auxiliary"
            end if
            call my(wavefunctions_es(el%o%hc,multibasis(usage,el%o%cutoff,nb,kpt,lay)),el%o%wf(ik)) ; if (error()) exit
          end do

          call sync_config_process_errors() ; if (error()) goto 300

          ! Eigenvalues and occupations
          allocate( el%o%eigs(nk,nb), el%o%occs(nk,nb) )
          call eigenvalues_i(el%o) ; if (error()) goto 300
          call occupations_i(el%o) ; if (error()) goto 300

          ! Residual norm
          call residual_norm_i(el%o)

          call diary_construction_i(el%o)

        end if

300     call glean(thy(lay))

        call glean(thy(ext))
        call glean(thy(gp))
        if (present(restf)) call glean(thy(restf))

999     if (error("Exit electrons_fh_mod::constructor_el")) continue

        if (.not.error()) call stop_timer("electrons_fh: constructor")

      end function

      subroutine update_el(el)
!doc$ subroutine update(el)
        type(electrons_fh_obj) :: el
!       modifies: el
!       effects: Updates el.
!       errors: Passes errors.

!cod$
        integer :: ik

        call start_timer("electrons_fh: update")

        call my(el)

        call own_i(el)

        do ik = 1,size(el%o%wf)
          call update(el%o%wf(ik),el%o%hc) ; if (error()) exit
        end do

        call sync_config_process_errors() ; if (error()) goto 100

        call eigenvalues_i(el%o) ; if (error()) goto 100
        call occupations_i(el%o) ; if (error()) goto 100

        call residual_norm_i(el%o)

100     call glean(thy(el))

        if (error("Exit electrons_fh_mod::update_el")) continue

        if (.not.error()) call stop_timer("electrons_fh: update")

      end subroutine

      subroutine my_el(el)
!doc$ subroutine my(el)
        type(electrons_fh_obj) :: el

!cod$
        el%ref = el%ref + 1
        el%o%ref = el%o%ref + 1
      end subroutine

      subroutine my_new_el(eli,el)
!doc$ subroutine my(eli,el)
        type(electrons_fh_obj) :: eli, el

!cod$
        el%ref = 1
        el%o => eli%o
        el%o%ref = el%o%ref + 1
      end subroutine

      function thy_el(el) result(elo)
!doc$ function thy(el) result(elo)
        type(electrons_fh_obj) :: el, elo

!cod$
        el%ref = el%ref - 1
        el%o%ref = el%o%ref - 1
        elo%ref = el%ref
        elo%o => el%o
      end function

      subroutine glean_el(el)
!doc$ subroutine glean(el)
        type(electrons_fh_obj) :: el

!cod$
        integer :: ik
        if (el%o%ref < 1) then
          if (associated( el%o%kgroup_index )) deallocate( el%o%kgroup_index )
          if (associated( el%o%eigs )) deallocate( el%o%eigs )
          if (associated( el%o%occs )) deallocate( el%o%occs )
          if (associated( el%o%wf )) then
            do ik = 1,size(el%o%wf)
              call glean(thy(el%o%wf(ik)))
            end do
            deallocate( el%o%wf )
          end if
          call glean(thy(el%o%kpoints))
          call glean(thy(el%o%hc))
          deallocate( el%o )
        end if
      end subroutine

      subroutine bequeath_el(el)
!doc$ subroutine bequeath(el)
        type(electrons_fh_obj) :: el

!cod$
        continue
      end subroutine

      subroutine assign_el(el,el2)
!doc$ subroutine assign(el,el2)
        type(electrons_fh_obj), intent(inout) :: el
        type(electrons_fh_obj), intent(in) :: el2

!cod$
        type(electrons_fh_obj) :: elt
        call my(el2)
        elt%o => el%o
        el%o%ref = el%o%ref - el%ref
        el%o => el2%o
        el%o%ref = el%o%ref + el%ref
        call glean(elt)
        call glean(thy(el2))
      end subroutine

      function el_ref(el) result(r)
!doc$ function x_ref(el) result(r)
        type(electrons_fh_obj) :: el
        integer, dimension(2) :: r
!       effects: Returns el%ref and el%o%ref.

!cod$
        r(1) = el%ref
        r(2) = el%o%ref
        call glean(el)
      end function

      function el_ghost(el) result(g)
!doc$ function x_ghost(el) result(g)
        type(electrons_fh_obj) :: el
        type(ghost) :: g
!       effects: Returns ghost of el.

!cod$
        call my(el)
        g = el%o%g
        call glean(thy(el))
      end function

      function el_wavefunctions(el,ik) result(wf)
!doc$ function x_wavefunctions(el,ik) result(wf)
        type(electrons_fh_obj) :: el
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
        if (error("Exit electrons_fh_mod::el_wavefunctions")) continue
      end function

      function el_n_bands(el) result(n)
!doc$ function x_n_bands(el) result(n)
        type(electrons_fh_obj) :: el
        integer :: n
!       effects: Returns the number of bands.

!cod$
        call my(el)
        n = size(el%o%eigs,2)
100     call glean(thy(el))
        if (error("Exit electrons_fh_mod::el_n_bands")) continue
      end function

      function el_eigenvalue(el,ik,ib) result(ev)
!doc$ function x_eigenvalue(el,ik,ib) result(ev)
        type(electrons_fh_obj) :: el
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
        if (error("Exit electrons_fh_mod::el_eigenvalue")) continue
      end function

      function el_eigenvalues(el) result(evs)
!doc$ function x_eigenvalues(el) result(evs)
        type(electrons_fh_obj) :: el
        real(double), dimension(size(el%o%eigs,1),size(el%o%eigs,2)) :: evs
!       effects: Returns the eigenvalues.

!cod$
        call my(el)
        evs = el%o%eigs
        call glean(thy(el))
      end function

      function el_residual_norm(el) result(rn)
!doc$ function x_residual_norm(el) result(rn)
        type(electrons_fh_obj) :: el
        real(double) :: rn
!       effects: Returns the weighted sum of the wavefunction residual norms.

!cod$
        call my(el)
        rn = el%o%res_norm
        call glean(thy(el))
      end function

      function el_converged(el) result(cvg)
!doc$ function x_converged(el) result(cvg)
        type(electrons_fh_obj) :: el
        logical :: cvg
!       effects: Returns the convergence status of el with respect to res_norm and res_norm_tol.

!cod$
        call my(el)
        cvg = el%o%res_norm_cvg
        call glean(thy(el))
      end function

      subroutine diary_el(el)
!doc$ subroutine diary(el)
        type(electrons_fh_obj) :: el
!       effects: Writes el information to the diary.

!cod$
        logical :: found
        character(line_len) :: tag
        integer :: ib, ik, ios, is, msg, nb, nk, nsg
        integer :: band_range(2), brl, bru, kpoint_range(2),krl, kru, file_format
        real(double) :: bocc1, bocc2, ne1, ne2, kpt(3)
        real(double), dimension(:), allocatable :: fl_c, fl_sg
        real(double), dimension(:,:,:), allocatable :: eigs_c, eigs_sg, occs_c, occs_sg
        real(double), pointer    :: rtmp(:,:,:)
        complex(double), pointer :: ctmp(:,:,:)
        type(grid_obj) :: g
        type(file_obj) :: f
        type(multivector_obj) :: mvec
        character(line_len)    :: filename
        integer               :: kind


        call my(el)

        nk = size(el%o%eigs,1)
        nb = size(el%o%eigs,2)
        nsg = mpi_nsgroups()

        ! Gather spin information to the CONFIG rank 0 process.
        select case (nsg)
        case (2)
          msg = mpi_mysgroup()
          allocate( fl_c(nsg), fl_sg(nsg) )
          fl_sg = 0.0_double
          fl_sg(msg) = el%o%fermi_level
          call xcomm_reduce(XSGROUP,MPI_SUM,fl_sg,fl_c)                                             ; if (error()) goto 200
          allocate( eigs_c(nsg,nk,nb), eigs_sg(nsg,nk,nb), occs_c(nsg,nk,nb), occs_sg(nsg,nk,nb) )
          eigs_sg = 0.0_double
          eigs_sg(msg,:,:) = el%o%eigs
          call xcomm_reduce(XSGROUP,MPI_SUM,eigs_sg,eigs_c)                                         ; if (error()) goto 200
          occs_sg = 0.0_double
          occs_sg(msg,:,:) = el%o%occs
          call xcomm_reduce(XSGROUP,MPI_SUM,occs_sg,occs_c)                                         ; if (error()) goto 200
        end select

        ! Write information to the diary file
        if (i_access( diaryfile() )) then

          ! Eigenvalues and occupations
          select case (nsg)
          case (1)
            select case (el%o%occupation_method)
            case (THERMAL)
              write(x_unit(diaryfile()),'(/,t4,"Fermi level = ",f13.10," Ryd")') el%o%fermi_level
            end select
            do ik = 1,size(el%o%eigs,1)
              write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
              write(x_unit(diaryfile()),'(/,21x,"band",8x,"eigenvalue (Ryd)",7x,"occupation")')
              write(x_unit(diaryfile()),'(19x,49("-"))')
              do ib = 1,size(el%o%eigs,2)
                write(x_unit(diaryfile()),'(21x,i4,9x,f13.10,11x,f6.4)') ib, el%o%eigs(ik,ib), el%o%occs(ik,ib)
              end do
            end do
          case (2)
            select case (el%o%occupation_method)
            case (THERMAL)
              write(x_unit(diaryfile()),'(/,t4,"Fermi level:",10x,f13.10,15x,f13.10," Ryd")') fl_c(1), fl_c(2)
            end select
            do ik = 1,size(eigs_c,2)
              write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
              write(x_unit(diaryfile()),'(/,11x,"band",8x,"eigenvalue (Ryd)",7x,"occupation", &
                                                    & 10x,"eigenvalue (Ryd)",7x,"occupation")')
              write(x_unit(diaryfile()),'(9x,92("-"))')
              do ib = 1,size(eigs_c,3)
                write(x_unit(diaryfile()),'(11x,i4,9x,f13.10,11x,f6.4,13x,f13.10,11x,f6.4)') ib, eigs_c(1,ik,ib), occs_c(1,ik,ib), &
                                                                                               & eigs_c(2,ik,ib), occs_c(2,ik,ib)
              end do
            end do
          end select

          ! Band occupations
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
          call my(grid(x_layout(el%o%hc),KGROUP),g)
          do ik=krl,kru
            call my(x_multivector(el%o%wf(ik)),mvec) ; if (error()) goto 300
            kpt = x_kpoint(el%o%kpoints,ik)
            do ib=brl,bru

              select case (nsg)
              case (1)
                call get_wf_filename_i(filename,file_format,ik,nk,ib,nb)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                call put(mvec,ib,g)
                call merge_grid_density(g)

                ! Write to file
                call write_to_file(g,filename,file_format,kind,kpt=kpt) 
                if (error()) goto 300

              case (2)
                ! First write out spin 1
                call get_wf_filename_i(filename,file_format,ik,nk,ib,nb,is=1)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                call put(mvec,ib,g)
                call merge_grid_density(g)

                call write_to_file(g,filename,file_format,kind,spin=1,kpt=kpt) 
                if (error()) goto 300

                ! Next write out spin 2
                call get_wf_filename_i(filename,file_format,ik,nk,ib,nb,is=2)

                call write_to_file(g,filename,file_format,kind,spin=2,kpt=kpt) 
                if (error()) goto 300
              end select
              call sgroup_to_kgroup(g)
            end do
            call glean(thy(mvec))
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
          call my(grid(x_layout(el%o%hc),KGROUP),g)
          call alloc(rtmp, x_layout(g), S_TYPE)
          do ik=krl,kru
            call my(x_multivector(el%o%wf(ik)),mvec) ; if (error()) goto 300
            kpt = x_kpoint(el%o%kpoints,ik)
            do ib=brl,bru

              select case (nsg)
              case (1)
                call get_sd_filename_i(filename,file_format,ik,nk,ib,nb)

                ! Extract the i'th band from the underlying multivector and put it on a grid
                call empty(g)
                call put(mvec,ib,g)
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
                call put(mvec,ib,g)
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
              call sgroup_to_kgroup(g)
              ! Take rtmp only to keep from having to re-alloc
              call take(rtmp,g,RS_KIND)
            end do
            call glean(thy(mvec))
          end do
          call put(rtmp,g,RS_KIND)
          call glean(thy(g))

        end select
       
200     call barrier(CONFIG)

        if (allocated( fl_c ))  deallocate( fl_c )
        if (allocated( fl_sg ))  deallocate( fl_sg )
        if (allocated( eigs_c ))  deallocate( eigs_c )
        if (allocated( eigs_sg )) deallocate( eigs_sg )
        if (allocated( occs_c ))  deallocate( occs_c )
        if (allocated( occs_sg )) deallocate( occs_sg )

300     call glean(thy(el))

        if (error("Exit electrons_fh_mod::diary_el")) continue

      end subroutine

      subroutine decompose_el(el,site_data,mode,f)
!doc$ subroutine decompose(el,site_data,mode,f)
        type(electrons_fh_obj) :: el
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

        if (error("Exit electrons_fh_mod::decompose_el")) continue

      end subroutine

! private routines

      subroutine own_i(el)
        type(electrons_fh_obj) :: el
        type(electrons_fh_obj) :: elt
        integer :: ik
        if (el%ref < el%o%ref) then
          allocate( elt%o )
          elt%o%ref = 0
          elt%o%g = el%o%g
          elt%o%res_norm_cvg = el%o%res_norm_cvg
          elt%o%occupation_method = el%o%occupation_method
          elt%o%res_norm = el%o%res_norm
          elt%o%res_norm_tol = el%o%res_norm_tol
          elt%o%total_charge = el%o%total_charge
          elt%o%spin_polarization = el%o%spin_polarization
          elt%o%charge_state = el%o%charge_state
          elt%o%fermi_level = el%o%fermi_level
          elt%o%kt = el%o%kt
          elt%o%cutoff = el%o%cutoff
          allocate( elt%o%kgroup_index(size(el%o%kgroup_index)) )
          elt%o%kgroup_index = el%o%kgroup_index
          allocate( elt%o%eigs(size(el%o%eigs,1),size(el%o%eigs,2)) )
          elt%o%eigs = el%o%eigs
          allocate( elt%o%occs(size(el%o%occs,1),size(el%o%occs,2)) )
          elt%o%occs = el%o%occs
          call my(el%o%kpoints,elt%o%kpoints)
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

      subroutine diary_construction_i(elr,restf)
        type(electrons_fh_rep) :: elr
        type(tagio_obj), optional :: restf

        integer :: ik, nk, nkg

        if (present(restf)) call my(restf)

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
          case (UNIFORM)
            write(x_unit(diaryfile()),'(/,t4,"Uniform occupations")')
          end select
        end if
        call diary(elr%hc)
        do ik = 1,nk
          if (mpi_mykgroup() /= elr%kgroup_index(ik)) cycle
          call diary(elr%wf(ik))
          exit
        end do

        if (present(restf)) call glean(thy(restf))

        if (error("Exit electrons_fh_mod::diary_construction_i")) continue

      end subroutine

      subroutine eigenvalues_i(elr)
        type(electrons_fh_rep) :: elr

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

        if (error("Exit electrons_fh_mod::eigenvalues_i")) continue

      end subroutine

      subroutine occupations_i(elr)
        type(electrons_fh_rep) :: elr

        real(double), parameter :: eftol = 1.0e-15_double
        real(double), parameter :: qtol = 1.0e-10_double
        integer :: ib, ik, it, msg, nb, nk, nsg
        real(double) :: d, emax, emin, q, total_charge, spin_degeneracy, x
        real(double), dimension(:), allocatable :: kwts
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

          end select

          ! Compute the occupations
          do ik = 1,nk
            do ib = 1,nb
              x = elr%eigs(ik,ib) - elr%fermi_level
              elr%occs(ik,ib) = spin_degeneracy*fermi_distr(x,elr%kt)
            end do
          end do

        case (UNIFORM)

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

        if (error("Exit electrons_fh_mod::occupations_i")) continue

      end subroutine

      subroutine notify_occupations_i(elr,it,emin,emax,fl,q,sd,kwts)
        type(electrons_fh_rep) :: elr
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
        type(electrons_fh_rep) :: elr

        integer :: ik
        real(double) :: rn_k, rn_kg, rn_sg

        rn_kg = 0.0_double
        do ik = 1,size(elr%wf)
          call residual_norm(elr%wf(ik),rn_k)
          rn_kg = rn_kg + rn_k*x_kweight(elr%kpoints,ik)
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,rn_kg,rn_sg)        ; if (error()) goto 100
        call xcomm_allreduce(XSGROUP,MPI_SUM,rn_sg,elr%res_norm) ; if (error()) goto 100
!        elr%res_norm = elr%res_norm/real(mpi_nsgroups(),double)

        elr%res_norm_cvg = (elr%res_norm < elr%res_norm_tol)

100     if (error("Exit electrons_fh_mod::residual_norm_i")) continue

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


100     if (error("electrons_fh_mod::get_wf_filename_i - Exiting")) continue

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


100     if (error("electrons_fh_mod::get_sd_filename_i - Exiting")) continue

      end subroutine 





      end module

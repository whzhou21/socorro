!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module electrons_td_mod
!doc$ module electrons_td_mod

!     One datatype is available here: type(electrons_td_obj)

!     electrons_sc_mod creates and maintains a set of electrons in calculations with a time-dependent hamiltonian..

      use kind_mod
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
      use symmetry_mod
      use atomic_density_mod
      use atomic_potential_mod
      use wavefunctions_td_mod
      use operators_mod
      use multibasis_mod
      use gen_density_mod
      use gen_potential_mod
      use electrons_sc_mod
      use time_propagator_mod
      use timing_mod
!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = CONFIG

      type :: electrons_td_rep
        integer :: ref
        type(ghost) :: g
!        type(ghost) :: g_external                                  ! ghost of external object used in construction
        real(double) :: kinetic_energy                             ! kinetic energy
        real(double) :: cutoff                                     ! wavefunctions cutoff energy
        integer, dimension(:), pointer :: kgroup_index             ! mapping of k-points to kgroups
        real(double), dimension(:,:), pointer :: exps              ! expectation values
        real(double), dimension(:,:), pointer :: occs              ! occupations
        type(kpoints_obj) :: kpoints                               ! k-points object
        type(gen_density_obj) :: density                           ! generalized density object
        type(h_common_obj) :: hc                                   ! common hamiltonian object
        type(wavefunctions_td_obj), dimension(:), pointer :: wf    ! wavefunction objects
      end type

      type, public :: electrons_td_obj
        private
        integer :: ref
        type(electrons_td_rep), pointer :: o
      end type

!doc$
      public :: electrons_td
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_energy
      public :: x_kinetic_energy
      public :: x_wavefunctions
      public :: x_density
      public :: x_kpoints
      public :: x_n_bands
      public :: x_expectation_value
      public :: x_expectation_values
      public :: x_occupation
      public :: x_occupations
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: diary
      public :: max_energy
      public :: diary_energies
      public :: decompose
      public :: write_restart
      public :: get_norm
      public :: num_hamiltonian_ops

!cod$
      interface electrons_td
        module procedure constructor_el, constructor_restart_el
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
      interface x_energy
        module procedure el_energy
      end interface
      interface x_kinetic_energy
        module procedure el_kinetic_energy
      end interface
      interface x_wavefunctions
         module procedure el_wavefunctions
      end interface
      interface x_kpoints
         module procedure el_kpoints
      end interface
      interface x_n_bands
         module procedure el_n_bands
      end interface
      interface x_expectation_value
         module procedure el_expectation_value
      end interface
      interface x_expectation_values
         module procedure el_expectation_values
      end interface
      interface x_occupation
         module procedure el_occupation
      end interface
      interface x_occupations
         module procedure el_occupations
      end interface
      interface x_density
        module procedure el_density
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
      interface write_restart
        module procedure write_restart_el
      end interface
      interface get_norm
         module procedure get_norm_el
      end interface
      interface num_hamiltonian_ops
         module procedure num_hamiltonian_ops_el
      end interface


      contains

! public routines

      function constructor_el(el_sc,ext) result(el)
!doc$ function electrons_td(el_sc,ext) result (el)
        type(electrons_sc_obj) :: el_sc
        type(external_obj) :: ext
        type(electrons_td_obj) :: el
!       requires: el_sc must have been constructed, el_sc should be consistent with ext
!       effects: Constructs a new el for time t = 0.
!       errors: Passes errors.
!       notes: This routine does not propagate the electrons so it does not receive a gen_potential.

!cod$
        call start_timer("electrons_td: constructor")

        if (error("  Error on entry")) goto 999

        call my(el_sc)
        call my(ext)

        el%ref = 0
        allocate( el%o )
        el%o%ref = 0
        el%o%g = x_ghost()

        call standard_portal_i(el%o,el_sc,ext) ; if (error()) goto 100

100     call glean(thy(ext))
        call glean(thy(el_sc))

999     if (error("Exit electrons_td_mod::constructor_el")) continue

        if (.not.error()) call stop_timer("electrons_td: constructor")

      end function

      function constructor_restart_el(restf) result(el)
!doc$ function electrons_td(restf) result (el)
        type(tagio_obj)        :: restf
        type(electrons_td_obj) :: el
!       requires: valid restart file restf
!       effects: Constructs a new el from file
!       errors: Passes errors.
!       notes: This routine does not propagate the electrons so it does not receive a gen_potential.

!cod$
        call start_timer("electrons_td: restart constructor")

        if (error("  Error on entry")) goto 999

        call my(restf)

        el%ref = 0
        allocate( el%o )
        el%o%ref = 0
        el%o%g = x_ghost()

        call restart_portal_i(el%o,restf) ; if (error()) goto 100

100     call glean(thy(restf))

999     if (error("Exit electrons_td_mod::constructor_el")) continue

        if (.not.error()) call stop_timer("electrons_td: restart constructor")

      end function


      subroutine update_el(el,ext,gp)
!doc$ subroutine update(el,ext,gp)
        type(electrons_td_obj) :: el
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

        call start_timer("electrons_td: update")

        call my(el)
        call my(ext)
        call my(gp)

        call own_i(el)

!        if ( x_ghost(ext) /= el%o%g_external ) el%o%g_external = x_ghost(ext)

        call update(el%o%hc,ext,gp)
        do ik = 1,size(el%o%wf)
          if (mpi_mykgroup() /= el%o%kgroup_index(ik)) cycle
          oldg = x_ghost(el%o%wf(ik))
          call update(el%o%wf(ik),el%o%hc) ; if (error()) goto 100
          e_change = (oldg /= x_ghost(el%o%wf(ik)))
        end do

        call sync_config_process_errors() ; if (error()) goto 100

        call xcomm_allreduce(XKGROUP,MPI_LOR,e_change,ge_change)
        if (ge_change) then
          el%o%g = x_ghost()
          el%o%density = gen_density(ext) ; if (error()) goto 100
          call accumulate_density_i(el%o,ext) ; if (error()) goto 100
        end if
        call kinetic_energy_i(el%o)

100     call glean(thy(el))
        call glean(thy(ext))
        call glean(thy(gp))

        if (error("Exit electrons_td_mod::update_el")) continue

        if (.not.error()) call stop_timer("electrons_td: update")

      end subroutine

      subroutine my_el(el)
!doc$ subroutine my(el)
        type(electrons_td_obj) :: el

!cod$
        el%ref = el%ref + 1
        el%o%ref = el%o%ref + 1
      end subroutine

      subroutine my_new_el(eli,el)
!doc$ subroutine my(eli,el)
        type(electrons_td_obj) :: eli, el

!cod$
        el%ref = 1
        el%o => eli%o
        el%o%ref = el%o%ref + 1
      end subroutine

      function thy_el(el) result(elo)
!doc$ function thy(el) result(elo)
        type(electrons_td_obj) :: el, elo

!cod$
        el%ref = el%ref - 1
        el%o%ref = el%o%ref - 1
        elo%ref = el%ref
        elo%o => el%o
      end function

      subroutine glean_el(el)
!doc$ subroutine glean(el)
        type(electrons_td_obj) :: el

!cod$
        integer :: ik
        if (el%o%ref < 1) then
          if (associated( el%o%kgroup_index )) deallocate( el%o%kgroup_index )
          if (associated( el%o%exps )) deallocate( el%o%exps )
          if (associated( el%o%occs )) deallocate( el%o%occs )
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
        type(electrons_td_obj) :: el

!cod$
        continue
      end subroutine

      subroutine assign_el(el,el2)
!doc$ subroutine assign(el,el2)
        type(electrons_td_obj), intent(inout) :: el
        type(electrons_td_obj), intent(in) :: el2

!cod$
        type(electrons_td_obj) :: elt
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
        type(electrons_td_obj) :: el
        integer, dimension(2) :: r
!       effects: Returns el%ref and el%o%ref.

!cod$
        r(1) = el%ref
        r(2) = el%o%ref
        call glean(el)
      end function

      function el_ghost(el) result(g)
!doc$ function x_ghost(el) result(g)
        type(electrons_td_obj) :: el
        type(ghost) :: g
!       effects: Returns ghost of el.

!cod$
        call my(el)
        g = el%o%g
        call glean(thy(el))
      end function

      function el_energy(el) result(e)
!doc$ function x_energy(el) result(e)
        type(electrons_td_obj) :: el
        real(double) :: e
!       effects: Returns the kinetic energy.

!cod$
        call my(el)
        e = el%o%kinetic_energy
        call glean(thy(el))
        if (error("Exit electrons_td_mod::el_energy")) continue

      end function

      function el_kinetic_energy(el) result(e)
!doc$ function x_kinetic_energy(el) result(e)
        type(electrons_td_obj) :: el
        real(double) :: e
!       effects: Returns kinetic energy.

!cod$
        call my(el)
        e = el%o%kinetic_energy
        call glean(thy(el))
      end function

      function el_wavefunctions(el,ik) result(wf)
!doc$ function x_wavefunctions(el,ik) result(wf)
        type(electrons_td_obj) :: el
        integer, intent(in) :: ik
        type(wavefunctions_td_obj) :: wf
!       effects: Returns the ik set of wavefunctions.
!       errors: ik out of range.

!cod$
        call my(el)
        if ( error((ik < 1) .or. (ik > size(el%o%wf)),"ERROR: ik is out of range")) goto 100
        call my(el%o%wf(ik),wf)
        call bequeath(thy(wf))
100     call glean(thy(el))
        if (error("Exit electrons_td_mod::el_wavefunctions")) continue
      end function 

      function el_density(el) result(dens)
!doc$ function  x_density(el) result(dens)
        type(electrons_td_obj) :: el
        type(gen_density_obj) :: dens
!       effects: Returns the density object
!cod$
        call my(el)
        call my(el%o%density,dens)
        call bequeath(thy(dens))
100     call glean(thy(el))
        if (error("Exit electrons_td_mod::el_density")) continue
      end function 

      function el_kpoints(el) result(kp)
!doc$ function x_kpoints(el) result(kp)
        type(electrons_td_obj) :: el
        type(kpoints_obj) :: kp
!       effects: Returns the kpoints_obj of el.

!cod$
        call my(el)
        call my(el%o%kpoints,kp)
        call glean(thy(el))
        call bequeath(thy(kp))
      end function
 
      function el_n_bands(el) result(n)
!doc$ function x_n_bands(el) result(n)
        type(electrons_td_obj) :: el
        integer :: n
!       effects: Returns the number of bands.

!cod$
        call my(el)
        n = size(el%o%exps,2)
100     call glean(thy(el))
        if (error("Exit electrons_td_mod::el_n_bands")) continue
      end function

      function el_expectation_value(el,ik,ib) result(ev)
!doc$ function x_expectation_value(el,ik,ib) result(ev)
        type(electrons_td_obj) :: el
        integer :: ik, ib
        real(double) :: ev
!       effects: Returns the expectation value for k-point ik and band ib.
!       errors: ik or ib out of range.

!cod$
        call my(el)
        if (error((ik < 1) .or. (ik > size(el%o%exps,1)),"ERROR: ik is out of range")) goto 100
        if (error((ib < 1) .or. (ib > size(el%o%exps,2)),"ERROR: ib is out of range")) goto 100
        ev = el%o%exps(ik,ib)
100     call glean(thy(el))
        if (error("Exit electrons_td_mod::el_expectation_value")) continue
      end function

      function el_expectation_values(el) result(evs)
!doc$ function x_expectation_values(el) result(evs)
        type(electrons_td_obj) :: el
        real(double), dimension(size(el%o%exps,1),size(el%o%exps,2)) :: evs
!       effects: Returns the expectation_values.

!cod$
        call my(el)
        evs = el%o%exps
        call glean(thy(el))
      end function

      function el_occupation(el,ik,ib) result(occ)
!doc$ function x_occupation(el,ik,ib) result(occ)
        type(electrons_td_obj) :: el
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
        if (error("Exit electrons_td_mod::el_occupation")) continue
      end function

      function el_occupations(el) result(occs)
!doc$ function x_occupations(el) result(occs)
        type(electrons_td_obj) :: el
        real(double), dimension(size(el%o%occs,1),size(el%o%occs,2)) :: occs
!       effects: Returns the occupations.

!cod$
        call my(el)
        occs = el%o%occs
        call glean(thy(el))
      end function

      subroutine forces_el(el,f)
!doc$ subroutine forces(el,f)
        type(electrons_td_obj) :: el
        real(double), dimension(:,:), intent(out) :: f
!       modifies: f
!       effects: Returns atomic forces due to electrons.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double), dimension(:), allocatable :: wts
        real(double), dimension(:,:), allocatable :: f1, f2

        call my(el)

        allocate( wts(size(el%o%exps,2)) )
        allocate( f1(size(f,1),size(f,2)), f2(size(f,1),size(f,2)) )

        f1 = 0.0_double
        do ik = 1,size(el%o%wf)
          if (mpi_mykgroup() /= el%o%kgroup_index(ik)) cycle
          wts = el%o%occs(ik,:)*x_kweight(el%o%kpoints,ik)      ; if (error()) goto 100
          call forces(el%o%wf(ik),wts,f2)                       ; if (error()) goto 100
          f1 = f1 + f2
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,f1,f)              ; if (error()) goto 100

        if (allocated( wts )) deallocate( wts )
        if (allocated( f1 )) deallocate( f1 )
        if (allocated( f2 )) deallocate( f2 )

100     call glean(thy(el))

        if (error("Exit electrons_td_mod::forces_el")) continue

      end subroutine

      subroutine pressure_el(el,p)
!doc$ subroutine pressure(el,p)
        type(electrons_td_obj) :: el
        real(double), intent(out) :: p
!       effects: Returns pressure contributions due to electrons.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double) :: p1, p2
        real(double), dimension(:), allocatable :: wts

        call my(el)

        allocate( wts(size(el%o%occs,2)) )

        p1 = 0.0_double
        do ik = 1,size(el%o%wf)
          if (mpi_mykgroup() /= el%o%kgroup_index(ik)) cycle
          wts = el%o%occs(ik,:)*x_kweight(el%o%kpoints,ik)      ; if (error()) goto 100
          call pressure(el%o%wf(ik),wts,p2)                     ; if (error()) goto 100
          p1 = p1 + p2
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,p1,p)              ; if (error()) goto 100

        if (allocated( wts )) deallocate( wts )

100     call glean(thy(el))

        if (error("Exit electrons_td_mod::pressure_el")) continue

      end subroutine

      subroutine stress_tensor_el(el,s)
!doc$ subroutine stress_tensor(el,s)
        type(electrons_td_obj) :: el
        real(double), dimension(:,:), intent(out) :: s
!       requires: s be dimension(3,3)
!       effects: Returns stress tensor contributions due to electrons.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double), dimension(:), allocatable :: wts
        real(double), dimension(3,3) :: s1, s2

        call my(el)

        allocate( wts(size(el%o%occs,2)) )

        s1 = 0.0_double
        do ik = 1,size(el%o%wf)
          if (mpi_mykgroup() /= el%o%kgroup_index(ik)) cycle
          wts = el%o%occs(ik,:)*x_kweight(el%o%kpoints,ik)     ; if (error()) goto 100
          call stress_tensor(el%o%wf(ik),wts,s2)               ; if (error()) goto 100
          s1 = s1 + s2
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,s1,s)             ; if (error()) goto 100

        if (allocated( wts )) deallocate( wts )

100     call glean(thy(el))

        if (error("Exit electrons_td_mod::stress_tensor_el")) continue

      end subroutine

      subroutine diary_el(el)
!doc$ subroutine diary(el)
        type(electrons_td_obj) :: el
!       effects: Writes el information to the diary.

!cod$
        integer :: ib, ik

        call my(el)

        if (i_access( diaryfile() )) then
           do ik = 1,size(el%o%exps,1)
             write(x_unit(diaryfile()),'(/,t6,"Special k-point #",i0,":")') ik
             write(x_unit(diaryfile()),'(/,21x,"band",8x,"expectation value (Ryd)",7x,"occupation")')
             write(x_unit(diaryfile()),'(19x,"-------------------------------------------------")')
             do ib = 1,size(el%o%exps,2)
               write(x_unit(diaryfile()),'(21x,i4,13x,f13.10,14x,f6.4)') ib, el%o%exps(ik,ib), el%o%occs(ik,ib)
             end do
           end do
        end if

        call glean(thy(el))

      end subroutine

      function max_energy_el(el) result(e)
!doc$ function max_energy(el) result(e)
        type(electrons_td_obj) :: el
        real(double) :: e
!       effects: Returns the kinetic energy.

        call my(el)
        e = abs(el%o%kinetic_energy)
        call glean(thy(el))
      end function

      subroutine diary_energies_el(el,md)
!doc$ subroutine diary_energies(el,md)
        type(electrons_td_obj) :: el
        integer, intent(in) :: md
!       modifies: Output stream
!       effects: Prints energies.

!cod$
        call my(el)

        if (i_access(diaryfile())) then
          select case (md)
          case (12)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f12.9)') el%o%kinetic_energy
          case (13)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f13.9)') el%o%kinetic_energy
          case (14)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f14.9)') el%o%kinetic_energy
          case (15)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f15.9)') el%o%kinetic_energy
          case (16)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f16.9)') el%o%kinetic_energy
          case (17)
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f17.9)') el%o%kinetic_energy
          case default
            write(x_unit(diaryfile()),'(t6,"kinetic               = ",f18.9)') el%o%kinetic_energy
          end select
          
          call flushbuf(diaryfile())
        end if

        call glean(thy(el))

      end subroutine

      function get_norm_el(el) result(norm)
!doc$ function get_norm(el) result(norm)
        type(electrons_td_obj) :: el
        real(double) :: norm
!       effects: Returns the largest normalization discrepancy in el.
!       errors: Passes errors.

!cod$
        integer :: ik
        real(double) :: err, tmp_err, tmp_norm

        call my(el)

        norm = get_norm(el%o%wf(1)) ; if (error()) goto 100

        err = abs(norm - 1.0_double)
        do ik = 2,size(el%o%wf)
          tmp_norm = get_norm(el%o%wf(ik)) ; if (error()) goto 100
          tmp_err = abs(tmp_norm - 1.0_double)
          if (err < tmp_err) then
            norm = tmp_norm
            err = tmp_err
          end if
        end do

100     call glean(thy(el))
        if (error("Exit electrons_mod::get_norm_el")) continue

      end function

      function num_hamiltonian_ops_el(el) result(num_h_ops)
!doc$ function num_hamiltonian_ops(el) result(num_h_ops)
        type(electrons_td_obj) :: el
        integer :: num_h_ops
!       effects: Returns the largest number of hamiltonian operations for the most recent time propagation operation
!       errors: Passes errors.
!cod$
        integer :: ik
        integer :: num_h_ops_tmp
        real(double) :: err, tmp_err, tmp_norm

        call my(el)

        num_h_ops = x_num_hamiltonian_ops(x_time_propagator(el%o%wf(1))) ; if (error()) goto 100

        do ik = 2,size(el%o%wf)
           num_h_ops_tmp = x_num_hamiltonian_ops(x_time_propagator(el%o%wf(ik))) ; if (error()) goto 100
           if (num_h_ops < num_h_ops_tmp) num_h_ops = num_h_ops_tmp
        end do
100     call glean(thy(el))
        if (error("Exit electrons_mod::num_hamiltonian_ops")) continue

      end function

      subroutine decompose_el(el,site_data,mode,f)
!doc$ subroutine decompose(el,site_data,mode,f)
        type(electrons_td_obj) :: el
        real(double), dimension(:,:), intent(in) :: site_data
        character(line_len), intent(in) :: mode
        type(file_obj) :: f
!       requires: x_unit(f) be open.
!       effects: Decomposes the Kohn-Sham functions into s, p, & d spherical harmonics around sites.
!       errors: Passes errors.

!cod$
        logical :: found
        integer :: b1, b2, bt, ib, ik, is, nb, nk, ns
        real(double) :: exp, emax, emin, emid, kwt, occ, range, sum_s, sum_p, sum_d
        real(double), dimension(3) :: kpt
        real(double), dimension(:,:,:), allocatable :: rsa, rsa_kg

        call my(el)
        call my(f)

        nk = size(el%o%exps,1)
        nb = size(el%o%exps,2)
        ns = size(site_data,2)

        call arg("dcomp_range",range,found)
        if (found) then
          if (error(range < 0.0_double,"ERROR: dcomp_range < 0")) goto 100
          do ib = 1,nb
            if (el%o%occs(1,ib) <= 0.5_double) then
              emid = el%o%exps(1,ib)
              exit
            end if
          end do
          emin = emid - 0.5_double*range
          emax = emid + 0.5_double*range
        else
          emin = minval(el%o%exps)
          emax = maxval(el%o%exps)
        end if

        b1 = nb
        b2 = 1
        do ik = 1,nk
          bt = 1
          do ib = 1,nb
            if (el%o%exps(ik,ib) >= emin) exit
            bt = ib
          end do
          b1 = min(b1,bt)
          bt = nb
          do ib = nb,1,-1
            if (el%o%exps(ik,ib) <= emax) exit
            bt = ib
          end do
          b2 = max(b2,bt)
        end do
        nb = b2 - b1 + 1
        if (nb < 0) then
          call warn("WARNING: aborting decomposition because nb < 0")
          goto 100
        end if

        allocate( rsa(9,nb,ns), rsa_kg(9,nb,ns) )

        do ik = 1,size(el%o%wf)
          if (mpi_mykgroup() == el%o%kgroup_index(ik)) then
            call decompose(el%o%wf(ik),site_data,mode,rsa_kg,b1) ; if (error()) goto 100
          else
            rsa_kg = 0.0_double
          end if
          call xcomm_reduce(XKGROUP,MPI_SUM,rsa_kg,rsa) ; if (error()) goto 100
          if (i_access(f)) then
            kpt = x_kpoint(el%o%kpoints,ik) ; kwt = x_kweight(el%o%kpoints,ik)
            write(x_unit(f),'(/,t2,"k-point #",i0,":",3f9.5,"; weight = ",f7.5)') ik, kpt, kwt
            do ib = 1,nb
              bt = ib - 1 + b1
              exp = el%o%exps(ik,bt)
              occ = el%o%occs(ik,bt)
              write(x_unit(f),'(/,t3,"band #",i0,": expectation value = ",f9.5,"; occupation = ",f7.5)') bt, exp, occ
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
        end do

        if (allocated( rsa )) deallocate( rsa )
        if (allocated( rsa_kg )) deallocate( rsa_kg )

100     call glean(thy(el))
        call glean(thy(f))

        if (error("Exit electrons_td_mod::decompose_el")) continue

      end subroutine


      subroutine write_restart_el(el,nrestf)
!doc$ subroutine write_restart(el,nrestf)
        type(electrons_td_obj) :: el
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes el restart information to nrestf.
!       errors: Passes errors.

!cod$
        integer :: ik, nb
        integer(long) :: dsize, ios, ndata, s4
        real(double), dimension(:), allocatable :: exps, occs
        type(gen_potential_obj)    :: hc_genpot
        type(external_obj)         :: hc_ext

        call my(el)
        call my(nrestf)

        if (i_access(nrestf)) call startblock(nrestf,"ELECTRONS")

! Ryan: Not sure why the call to write external is here instead of in config_td_mod.
!       Likewise not sure why the calls to write out densities and potentials are
!        here instead of in fields_td_mod.

!        if (i_access(nrestf)) call writetag(nrestf,"HC_EXT")
!        call my(external(x_crystal(el%o%hc)),hc_ext)
!        call write_restart(hc_ext,nrestf)


!        !** Write out components needed for the density member object
!        if (i_access(nrestf)) call writetag(nrestf,"DENSITY_GDEN")
!        call write_restart(x_grid_density(el%o%density),nrestf)

!        if (i_access(nrestf)) call writetag(nrestf,"DENSITY_ADEN")
!        call write_restart(x_atomic_density(el%o%density),nrestf)
        
!        !** Continue writing out quantities needed to construct the hc object
!        ! gpot
!        call my(get_genpot(el%o%hc),hc_genpot)

!        if (i_access(nrestf)) call writetag(nrestf,"HC_GRIDPOT")
!        call write_restart(x_grid_potential(hc_genpot),nrestf)

!        if (i_access(nrestf)) call writetag(nrestf,"HC_APOT")
!        call write_restart(x_atomic_potential(hc_genpot),nrestf)

        if (i_access(nrestf)) then
          call writetag(nrestf,"HC_PROTOBASIS")
          dsize = sizeof_double ; ndata = 1
          call writef(el%o%cutoff,dsize,ndata,x_tagfd(nrestf),ios)
          nb = size(el%o%exps,2)
          s4 = nb ; dsize = sizeof_long ; ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),ios)
        end if

        !** Write out kpoints
        call write_restart(el%o%kpoints,nrestf) ; if (error()) goto 100

        !** Wavefunctions
        if (i_access(nrestf)) call writetag(nrestf,"EL_WAVEFUNCS")
        do ik = 1,size(el%o%wf)
          call write_restart(el%o%wf(ik),nrestf) ; if (error()) goto 100
        end do

        !** Expectation values
        if (i_access(nrestf)) call writetag(nrestf,"EXPECTATIONVALUES")
        allocate( exps(size(el%o%exps,2)) )
        do ik = 1,size(el%o%wf)
           exps = el%o%exps(ik,:)
           if (i_access(nrestf)) then
              dsize = sizeof_double ; ndata = size(exps)
              call writef(exps,dsize,ndata,x_tagfd(nrestf),ios)
           end if
        end do

        !** Occupations
        if (i_access(nrestf)) call writetag(nrestf,"OCCUPATIONS")
        allocate( occs(size(el%o%occs,2)) )
        do ik = 1,size(el%o%wf)
           occs = el%o%occs(ik,:)
           if (i_access(nrestf)) then
              dsize = sizeof_double ; ndata = size(occs)
              call writef(occs,dsize,ndata,x_tagfd(nrestf),ios)
           end if
        end do


100     if (i_access(nrestf)) call endblock(nrestf) ; if (error()) goto 200

        if (allocated( exps )) deallocate( exps )
        if (allocated( occs )) deallocate( occs )

200     call glean(thy(el))
        call glean(thy(hc_genpot))
        call glean(thy(hc_ext))
        call glean(thy(nrestf))

        if (error("Exit electrons_td_mod::write_restart_el")) continue

      end subroutine

! private routines

      subroutine standard_portal_i(elr,el_sc,ext)
        type(electrons_td_rep) :: elr
        type(electrons_sc_obj) :: el_sc
        type(external_obj) :: ext

        integer :: ik, ib, nb, nk
        real(double) :: neb, nek
        type(electrons_sc_rep), pointer :: elr_sc
        real(double), dimension(:,:), pointer :: exps_local
        character(line_len) :: tag
        logical :: fnd

        call my(el_sc)
        call my(ext)

        elr_sc => wormhole(el_sc)

        nk = size(elr_sc%eigs,1)
        nb = size(elr_sc%eigs,2)

        
        if (error((elr_sc%g_external /= x_ghost(ext)),"ERROR: el_sc is not consistent with ext")) goto 100
!        elr%g_external = elr_sc%g_external

        ! copy the common part of the hamiltonian
        call my(elr_sc%hc,elr%hc)

        ! copy the k-points
        call my(elr_sc%kpoints,elr%kpoints)

        ! copy the kgroup information
        allocate( elr%kgroup_index(nk) )
        elr%kgroup_index = elr_sc%kgroup_index

        ! copy the wavefunctions
        allocate( elr%wf(nk) )
        do ik = 1,nk
          call my(wavefunctions_td(elr_sc%wf(ik)),elr%wf(ik))
        end do

        ! get the expectation values
        allocate( elr%exps(nk,nb), exps_local(nk,nb) )
        exps_local = 0.0_double
        do ik = 1,nk
          if (mpi_mykgroup() /= elr%kgroup_index(ik)) cycle
          exps_local(ik,:) = x_expectation_values(elr%wf(ik))
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,exps_local,elr%exps) ; if (error()) goto 100
        deallocate(exps_local)

        ! get the occupations
        allocate( elr%occs(nk,nb) )

        call arg("tddft_occupations", tag, fnd)
        if (.not.fnd) tag = "mermin-wagner"
        select case(tag(1:len_trim(tag)))
        case ("mermin-wagner","Mermin-Wagner")

           elr%occs = elr_sc%occs

        case ("runge-gross","Runge-Gross")

           neb = 2.0_double/real(mpi_nsgroups(),double)
           do ik = 1,nk
              nek = sum(elr_sc%occs(ik,:))
              do ib = 1, nb
                if (nek > neb) then
                   elr%occs(ik,ib) = neb
                   nek = nek - neb
                else
                   elr%occs(ik,ib) = nek
                   nek = 0.0_double
                endif
              end do
           end do

        case default

           if (error(.true.,"ERROR: tddft_occupations tag value not recognized")) goto 100    

        end select

        ! get the electronic density
        call my(gen_density(ext),elr%density) ; if (error()) goto 100
        call accumulate_density_i(elr,ext) ; if (error()) goto 100

        ! get the kinetic energy
        call kinetic_energy_i(elr)

        ! copy the cutoff
        elr%cutoff = elr_sc%cutoff

        call diary_standard_construction_i()

        nullify( elr_sc )

100     call glean(thy(ext))
        call glean(thy(el_sc))

        if (error("Exit electrons_td_mod::standard_portal_i")) continue

      end subroutine


      subroutine restart_portal_i(elr,restf)
        type(electrons_td_rep) :: elr
        type(tagio_obj) :: restf

        character(1)                :: tios
        character(line_len)          :: usage
        integer                     :: ik, nb, nk
        integer(long)               :: dsize, ios, ndata, s4
        real(double), dimension(3)  :: kpt
        real(double), allocatable   :: exps(:), occs(:)
        type(layout_obj)            :: lay
        type(external_obj)          :: ext
        type(grid_obj)              :: gridden
        type(atomic_density_obj)    :: adens
        type(grid_obj)              :: gridpot
        type(atomic_potential_obj)  :: apot
        type(gen_potential_obj)     :: genpot


        call my(restf)

        if (i_access(restf)) tios = findfirsttag(restf,"ELECTRONS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: ELECTRONS block was not found")) goto 100

        if (i_access(restf)) call openblock(restf)

!        !** First read in the external object that was stored from the hc object
!        if (i_access(restf)) tios = findfirsttag(restf,"HC_EXT")
!        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
!        if (error(tios == TAG_NOT_FOUND,"ERROR: HC_EXT tag was not found")) goto 100
        call my(external(restf),ext)
        
        call my(x_layout(ext),lay)

!        !** Read in quantities needed to build the density object
!        ! read in the grid density object
!        if (i_access(restf)) tios = findfirsttag(restf,"DENSITY_GDEN")
!        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
!        if (error(tios == TAG_NOT_FOUND,"ERROR: DENSITY_GDEN tag was not found")) goto 100
!        call my(grid(x_layout(ext),MOD_SCOPE),gridden)
!        call read_restart(gridden,restf) ; if (error()) goto 100
        
!        ! read in the atomic density object
!        if (i_access(restf)) tios = findfirsttag(restf,"DENSITY_ADEN")
!        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
!        if (error(tios == TAG_NOT_FOUND,"ERROR: DENSITY_ADEN tag was not found")) goto 100
!        call my(atomic_density(x_atomic_operators(ext),restf),adens)

!        ! construct the density object from the grid density and atomic density
!        call my(gen_density(gridden,adens),elr%density)
        
!        !** Next read in the rest of the quantities needed to build the hc member object
!        !** Read in quantities needed to construct the general potential
!        !  grid potential
!        if (i_access(restf)) tios = findfirsttag(restf,"HC_GRIDPOT")
!        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
!        if (error(tios == TAG_NOT_FOUND,"ERROR: HC_GRIDPOT tag was not found")) goto 100
!        call my(grid(x_layout(ext),MOD_SCOPE),gridpot)
!        call read_restart(gridpot,restf) ; if (error()) goto 100

!        !  atomic potential
!        if (i_access(restf)) tios = findfirsttag(restf,"HC_APOT")
!        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
!        if (error(tios == TAG_NOT_FOUND,"ERROR: HC_APOT tag was not found")) goto 100
!        call my(atomic_potential(adens,restf=restf),apot)

!        ! construct generalized potential object
!        call my(gen_potential(gridpot,apot),genpot)

        if (i_access(restf)) tios = findfirsttag(restf,"HC_PROTOBASIS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: HC_PROTOBASIS tag was not found")) goto 100

        ! read the wavefunctions cutoff and number of bands, and then form the protobasis
        if (i_access(restf)) then
          dsize = sizeof_double ; ndata = 1
          call readf(elr%cutoff,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,elr%cutoff)

        if (i_access(restf)) then
          dsize = sizeof_long ; ndata = 1
          call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
          nb = s4
        end if
        if (i_comm(restf)) call broadcast(FILE_SCOPE,nb)

        ! Construct the hc object
        call my(common_hamiltonian(ext,genpot,elr%cutoff,nb),elr%hc)
        
        !** Read in the kpoints object
        call my(kpoints(ext,restf),elr%kpoints) ; if (error()) goto 100
        nk = x_n_kpoints(elr%kpoints)

        ! Divide the k-points among processors
        if (error(nk < mpi_nkgroups(),"ERROR: nk is less than kgroups")) then
          call notify("Number of k-points = ",nk)
          goto 100
        end if
        if (mod(nk,mpi_nkgroups()) /= 0) then
          call warn("WARNING: non-equal division of k-points among kgroups")
        end if
        allocate( elr%kgroup_index(nk) )
        do ik = 1,nk
          elr%kgroup_index(ik) = mod(ik-1,mpi_nkgroups()) + 1
        end do

        !** Read in the wavefunctions
        if (i_access(restf)) tios = findfirsttag(restf,"EL_WAVEFUNCS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: EL_WAVEFUNCS tag was not found")) goto 100
        ! initialize the wavefunctions
        allocate( elr%wf(nk) )
        do ik = 1,nk
           kpt = x_kpoint(elr%kpoints,ik) ; if (error()) goto 100
           usage = "normal"
           if (mpi_mykgroup() /= elr%kgroup_index(ik)) usage = "auxiliary"
           call my(wavefunctions_td(elr%hc,multibasis(usage,elr%cutoff,nb,kpt,lay),restf),elr%wf(ik)) ; if (error()) goto 100
        end do

        !** Read in the expectation values
        if (i_access(restf)) tios = findfirsttag(restf,"EXPECTATIONVALUES")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: EXPECTATIONVALUES tag was not found")) goto 100
        allocate(elr%exps(nk,nb), exps(nb))
        do ik = 1,nk
          if (i_access(restf)) then
            dsize = sizeof_double ; ndata = nb
            call readf(exps,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,exps)
          elr%exps(ik,:) = exps(:)
        end do
        deallocate(exps)
        

        !** Read in the occupations
        if (i_access(restf)) tios = findfirsttag(restf,"OCCUPATIONS")
        if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
        if (error(tios == TAG_NOT_FOUND,"ERROR: OCCUPATIONS tag was not found")) goto 100
        allocate(elr%occs(nk,nb), occs(nb))
        do ik=1,nk
           if (i_access(restf)) then
              dsize = sizeof_double ; ndata = nb
              call readf(occs,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           end if
           if (i_comm(restf)) call broadcast(FILE_SCOPE,occs)
           elr%occs(ik,:) = occs(:)
        end do
        deallocate(occs)
        

        call sync_config_process_errors() ; if (error()) goto 100

        call diary_restart_construction_i(elr,lay)


        if (i_access(restf)) call closeblock(restf)


        call glean(thy(lay))
        call glean(thy(ext))
        call glean(thy(gridden))
        call glean(thy(adens))
        call glean(thy(gridpot))
        call glean(thy(apot))
        call glean(thy(genpot))
        call glean(thy(restf))

100     if (error("Exit electrons_td_mod::restart_portal_i")) continue

      end subroutine



      subroutine accumulate_density_i(elr,ext)
        type(electrons_td_rep) :: elr
        type(external_obj) :: ext

        integer :: ik
        real(double), dimension(:), allocatable :: wts

        call my(ext)

        allocate( wts(size(elr%occs,2)) )
        do ik = 1,size(elr%wf)
          if (mpi_mykgroup() /= elr%kgroup_index(ik)) cycle
          wts = elr%occs(ik,:)*x_kweight(elr%kpoints,ik)     ; if (error()) goto 100
          call add_density(elr%wf(ik),wts,elr%density)       ; if (error()) goto 100
        end do
        call merge_symmetrize_filter(elr%density,ext)        ; if (error()) goto 100

        if (allocated( wts )) deallocate( wts )

100     call glean(thy(ext))

        if (error("Exit electrons_td_mod::accumulate_density_i")) continue

      end subroutine

      subroutine kinetic_energy_i(elr)
        type(electrons_td_rep) :: elr

        integer :: ik
        real(double) :: e
        real(double), dimension(:), allocatable :: wts

        allocate( wts(size(elr%occs,2)) )
        e = 0.0_double
        do ik = 1,size(elr%wf)
          if (mpi_mykgroup() /= elr%kgroup_index(ik)) cycle
          wts = elr%occs(ik,:)*x_kweight(elr%kpoints,ik)                 ; if (error()) goto 100
          e = e + kinetic_energy(elr%wf(ik),wts)                         ; if (error()) goto 100
        end do
        call xcomm_allreduce(XKGROUP,MPI_SUM,e,elr%kinetic_energy)       ; if (error()) goto 100

        if (allocated( wts )) deallocate( wts )

100     if (error("Exit electrons_td_mod::kinetic_energy_i")) continue

      end subroutine

      subroutine diary_standard_construction_i()

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Electrons object construction:")')
        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Fixed occupations")')

      end subroutine

      subroutine diary_restart_construction_i(elr,lay)
        type(electrons_td_rep) :: elr
        type(layout_obj) :: lay

        integer :: ik, nk, ng

        call my(lay)

        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,"Electrons object construction:")')
        call diary(elr%kpoints)
        nk = x_n_kpoints(elr%kpoints)
        ng = mpi_nkgroups()
        if (ng == 1) then
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
              write(x_unit(diaryfile()),'(t6,i5,8x,i4,8x,i7)') ik, (mod(ik-1,ng) + 1), x_n_gvectors(elr%wf(ik))
            end do
          end if
        end if
        if (i_access(diaryfile())) write(x_unit(diaryfile()),'(/,t4,"Fixed occupations")')
        call diary(elr%hc)
        do ik = 1,nk
          if (mpi_mykgroup() /= elr%kgroup_index(ik)) cycle
          call diary(elr%wf(ik))
          exit
        end do

100     call glean(thy(lay))

        if (error("Exit electrons_td_mod::diary_restart_construction_i")) continue

      end subroutine

      subroutine own_i(el)
        type(electrons_td_obj) :: el
        type(electrons_td_obj) :: elt
        integer :: ik
        if (el%ref < el%o%ref) then
          allocate( elt%o )
          elt%o%ref = 0
          elt%o%g = el%o%g
!          elt%o%g_external = el%o%g_external
          elt%o%kinetic_energy = el%o%kinetic_energy
          elt%o%cutoff = el%o%cutoff
          allocate( elt%o%kgroup_index(size(el%o%kgroup_index)) )
          elt%o%kgroup_index = el%o%kgroup_index
          allocate( elt%o%exps(size(el%o%exps,1),size(el%o%exps,2)) )
          elt%o%exps = el%o%exps
          allocate( elt%o%occs(size(el%o%occs,1),size(el%o%occs,2)) )
          elt%o%occs = el%o%occs
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

      end module

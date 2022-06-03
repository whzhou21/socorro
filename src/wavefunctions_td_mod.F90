!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module wavefunctions_td_mod
!doc$ module wavefunctions_td_mod

!     One datatype is available here: type(wavefunctions_td_obj).

!     wavefunctions_td_mod creates and maintains a set of wavefunctions in calculations involving time propagation.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod
      use tagio_mod
      use arg_mod
      use diary_mod
      use crystal_mod
      use operators_mod
      use multivector_mod
      use multibasis_mod
      use time_propagator_mod
      use ghost_mod
      use gen_density_mod
      use wavefunctions_es_mod
      use math_mod

!cod$
      implicit none
      private
   
      integer, parameter :: MOD_SCOPE = KGROUP

      ! usage
      integer, parameter :: NORMAL    = 1
      integer, parameter :: AUXILIARY = 2

      type :: wavefunctions_td_rep
        integer :: ref
        type(ghost) :: g
!        type(ghost) :: g_crystal                        ! crystal ghost
        integer :: usage                                ! usage wrt kgroup
        type(time_propagator_obj) :: tp                 ! time propagator
        type(multivector_obj) :: mv                     ! multivector
        real(double), dimension(:), pointer :: exps     ! expectation values
        type(h_kpoint_obj) :: hk                        ! k-point dependent part of the hamiltonian
      end type

      type, public :: wavefunctions_td_obj
        private
        integer :: ref
        type(wavefunctions_td_rep), pointer :: o
      end type

!doc$
      public :: wavefunctions_td
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_multivector
      public :: x_h_kpoint
      public :: x_n_gvectors
      public :: x_expectation_value
      public :: x_expectation_values
      public :: x_time_propagator
      public :: kinetic_energy
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: add_density
      public :: get_norm
      public :: decompose
      public :: diary
      public :: write_restart

!cod$
      interface wavefunctions_td
        module procedure constructor_wf, constructor_restart_wf
      end interface
      interface update
        module procedure update_wf
      end interface
      interface my
        module procedure my_wf, my_new_wf
      end interface
      interface thy
        module procedure thy_wf
      end interface
      interface glean
        module procedure glean_wf
      end interface
      interface bequeath
        module procedure bequeath_wf
      end interface
      interface assignment(=)
        module procedure assign_wf
      end interface
      interface x_ref
        module procedure wf_ref
      end interface
      interface x_ghost
        module procedure wf_ghost
      end interface
      interface x_multivector
        module procedure wf_multivector
      end interface
      interface x_h_kpoint
        module procedure wf_h_kpoint
      end interface
      interface x_n_gvectors
        module procedure wf_n_gvectors
      end interface
      interface x_expectation_value
        module procedure wf_expectation_value
      end interface
      interface x_expectation_values
        module procedure wf_expectation_values
      end interface
      interface x_time_propagator
        module procedure wf_time_propagator
      end interface
      interface kinetic_energy
        module procedure kinetic_energy_wf
      end interface
      interface forces
        module procedure forces_wf
      end interface
      interface pressure
        module procedure pressure_wf
      end interface
      interface stress_tensor
        module procedure stress_tensor_wf
      end interface
      interface add_density
        module procedure add_density_wf
      end interface
      interface get_norm
        module procedure get_norm_wf
      end interface
      interface decompose
        module procedure decompose_wf
      end interface
      interface diary
        module procedure diary_wf
      end interface
      interface write_restart
        module procedure write_restart_wf
      end interface

      contains

! public routines

      function constructor_wf(wf_es) result(wf)
!doc$ function wavefunctions_td(wf_es) result(wf)
        type(wavefunctions_es_obj) ::wf_es
        type(wavefunctions_td_obj) :: wf
!       effects: Constructs a new wf for time t = 0.
!       errors: Passes errors.

!cod$
        type(wavefunctions_es_rep), pointer :: wfr_es
        character(line_len) :: tag
        logical :: found
        integer :: nb, ib, jb

        integer, parameter :: quad = SELECTED_REAL_KIND(32)
        real(quad) :: kt
        real(quad), dimension(:), pointer :: eigs, weights
        complex(quad), dimension(:,:), pointer :: coeff, overlap
        complex(quad) :: prod
        real(quad) :: norm, moe_qk

        real(double) :: ktd, moe_dk, moe
        real(double), dimension(:,:), pointer :: phases
        complex(double), dimension(:,:), pointer :: coefficients


        if (error("  Error on entry")) then
          wf%ref = 0
          allocate( wf%o )
          wf%o%ref = 0
          goto 999
        end if

        call my(wf_es)
        
        wf%ref = 0
        allocate( wf%o )
        wf%o%ref = 0
        wf%o%g = x_ghost()

        wfr_es => wormhole(wf_es)
        wf%o%usage = wfr_es%usage

        select case (wf%o%usage)
        case (NORMAL)
!          wfr%g_crystal = wfr_es%g_crystal
          call my(time_propagator(),wf%o%tp) ; if (error()) goto 100
          call my(wfr_es%mv,wf%o%mv)

          call arg("tddft_thermalize_wavefunctions", tag, found)
          if (.not.found) tag = "no"
          select case(tag(1:len_trim(tag)))
          case ("yes","on","true")
             nb = size(wfr_es%eigs)
             allocate( coefficients(nb,nb) )

             if (mpi_first(KGROUP)) then

                allocate( eigs(nb), weights(nb), phases(nb,nb), coeff(nb,nb), overlap(nb,nb) )

                call arg("tddft_thermalize_kt",ktd,found)
                if (found) then
                   if (error(ktd <= 0.0_double,"ERROR: tddft_thermalize_ktd <= 0.0")) goto 100
                else
                   ktd = 0.05_double
                end if
                kt = ktd  !  Convert from double to quad
                eigs = wfr_es%eigs  !  Convert from double to quad
                weights = exp(-eigs/(2.0_quad*kt))


                call random_number(phases)
                phases = two_pi*phases
                coeff = cmplx(cos(phases),sin(phases)) !  Convert from double to quad
                do ib = 1, nb
                   coeff(:,ib) = coeff(:,ib)*weights
                end do
                do ib = 1, nb
                   do jb = 1, ib - 1
                      prod = dot_product(coeff(:,jb),coeff(:,ib))
                      coeff(:,ib) = coeff(:,ib) - prod*coeff(:,jb)
                   end do
                   norm = sqrt(dot_product(coeff(:,ib),coeff(:,ib)))
                   coeff(:,ib) = coeff(:,ib) / norm
                end do

                overlap = matmul(coeff,conjg(transpose(coeff)))

                moe_qk = 0.0_quad
                do ib = 1, nb
                   do jb = 1, nb
                      if ((ib /= jb) .and. (abs(overlap(ib,jb)) .gt. moe_qk)) then
                         moe_qk = abs(overlap(ib,jb))
                      end if
                   end do
                end do

                moe_dk = moe_qk  ! convert from quad to double
                
                coefficients = coeff  ! convert from quad to double

                deallocate( eigs, weights, phases, coeff, overlap)

             end if
             call broadcast(KGROUP,coefficients)
             call broadcast(KGROUP,moe_dk)

             call xcomm_allreduce(XKGROUP,MPI_MAX,moe_dk,moe) ; if (error()) goto 100
             if (i_access( diaryfile() )) then
                write(x_unit(diaryfile()),'(/,t4,"Thermalized TDDFT wavefunction initialization")')
                write(x_unit(diaryfile()),'(/,t6,"Maximum orthogonalization error = ",es7.1)') moe
             end if
             if (error(moe > 1.0d-10,"ERROR: thermal wavefunctions are not sufficiently orthogonal")) goto 100

             call transform(wf%o%mv,coefficients)

             allocate( wf%o%exps(nb) )
             do ib = 1, nb
                wf%o%exps(ib) = sum(conjg(coefficients(:,ib))*coefficients(:,ib)*wfr_es%eigs)
             end do

             call my(wfr_es%hk,wf%o%hk)

             deallocate(coefficients)

          case default

             allocate( wf%o%exps(size(wfr_es%eigs)) )
             wf%o%exps = wfr_es%eigs
             call my(wfr_es%hk,wf%o%hk)

          end select

        case (AUXILIARY)
          call my(wfr_es%mv,wf%o%mv)
        end select

100     nullify( wfr_es )

        call glean(thy(wf_es))

999     if (error("Exit wavefunctions_td_mod::constructor_wf")) continue

      end function

      function constructor_restart_wf(hc,mb,restf) result(wf)
!doc$ function wavefunctions_td(hc,mb,restf) result(wf)
        type(h_common_obj)         :: hc
        type(multibasis_obj)       :: mb
        type(tagio_obj)            :: restf
        type(wavefunctions_td_obj) :: wf
!       requires: 
!       effects: Constructs a wf from a restart file
!       errors: Passes errors.
!cod$
        integer                     :: nb, ng
        character(1)                :: tios
        integer(long)               :: dsize, ios, ndata, s4

        if (error("  Error on entry")) then
          wf%ref = 0
          allocate( wf%o )
          wf%o%ref = 0
          goto 999
        end if

        call my(hc)
        call my(mb)
        call my(restf)

        wf%ref = 0
        allocate( wf%o )
        wf%o%ref = 0
        wf%o%g = x_ghost()

        if (i_access(restf)) tios = findnexttag(restf,"WAVEFUNCTIONS_TD")
        if (i_comm(restf)) call broadcast(MOD_SCOPE,tios)
        if (error(tios /= TAG_START_BLOCK,"ERROR: WAVEFUNCTIONS_TD block was not found")) goto 100

        if (i_access(restf)) call openblock(restf)

        select case (x_usage(mb))
        case (NORMAL)

           wf%o%usage = NORMAL

           !** Construct the k-point dependent hamiltonian
           call my(kpoint_hamiltonian(hc,mb),wf%o%hk)         ; if (error()) goto 100
           
           !** construct the time propagator object
           call my(time_propagator(mb,restf),wf%o%tp)         ; if (error()) goto 100

           !** Read in the expectation values
           if (i_access(restf)) tios = findfirsttag(restf,"EXPECTATIONVALUES")
           if (i_comm(restf)) call broadcast(MOD_SCOPE,tios)
           if (error(tios == TAG_NOT_FOUND,"ERROR: EXPECTATIONVALUES tag was not found")) goto 100
           nb = x_n_bands(mb)
           allocate( wf%o%exps(nb) )
           if (i_access(restf)) then
              dsize = sizeof_double ; ndata = nb
              call readf(wf%o%exps,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),ios)
           end if
           if (i_comm(restf)) call broadcast(MOD_SCOPE,wf%o%exps)

           !** Read in the multivector
           call my(multivector(mb,restf=restf),wf%o%mv) ; if (error()) goto 100

        case (AUXILIARY)

           wf%o%usage = AUXILIARY
           nullify(wf%o%exps)
           call my(multivector(mb,restf=restf),wf%o%mv) ; if (error()) goto 100

        end select

        if (i_access(restf)) call closeblock(restf)

100     call glean(thy(restf))

999     if (error("Exit wavefunctions_td_mod::constructor_restart_wf")) continue

      end function

      subroutine update_wf(wf,hc,mb)
!doc$ subroutine update(wf,hc,mb)
        type(wavefunctions_td_obj) :: wf
        type(h_common_obj) :: hc
        type(multibasis_obj), optional :: mb
!       requires: wf%o%usage = NORMAL
!       modifies: wf
!       effects: Updates wf with respect to its dependencies.
!       errors: mb change. Passes errors.

!cod$
!        logical :: crystal_change
        logical :: mb_change
        logical :: debug
        !-----------------------------------
        
        debug = .false.
        
        if (debug) call warn("wavefunctions_td::update_wf starting")

        call my(wf)
        call my(hc)
        if (present(mb)) call my(mb)

!        crystal_change = ( x_ghost(x_crystal(hc)) /= wf%o%g_crystal )
        mb_change = .false.

        if (present(mb)) then
          mb_change = ( x_ghost(x_multibasis(wf%o%mv)) /= x_ghost(mb) )
          if (error(mb_change,"ERROR: multibasis change is not currently allowed")) goto 100
        end if

        call own_i(wf)
        wf%o%g = x_ghost()
        call update(wf%o%hk,hc) ; if (error()) goto 100

        if (debug) call warn("wavefunctions_td::update_wf calling propagate")

        call propagate(wf%o%tp,hamiltonian(wf%o%hk),wf%o%mv,wf%o%exps) ; if (error()) goto 100

100     call glean(thy(wf))
        call glean(thy(hc))
        if (present(mb)) call glean(thy(mb))

        if (error("Exit wavefunctions_td_mod::update_wf")) continue

        if (debug) call warn("wavefunctions_td::update_wf exiting")

      end subroutine

      subroutine my_wf(wf)
!doc$ subroutine my(wf)
        type(wavefunctions_td_obj) :: wf

!cod$
        wf%ref = wf%ref + 1
        wf%o%ref = wf%o%ref + 1
      end subroutine

      subroutine my_new_wf(wfi,wf)
!doc$ subroutine my(wfi,wf)
        type(wavefunctions_td_obj) :: wfi, wf

!cod$
        wf%ref = 1
        wf%o => wfi%o
        wf%o%ref = wf%o%ref + 1
      end subroutine

      function thy_wf(wf) result(wfo)
!doc$ function thy(wf) result(wfo)
        type(wavefunctions_td_obj) :: wf, wfo

!cod$
        wf%ref = wf%ref - 1
        wf%o%ref = wf%o%ref - 1
        wfo%ref = wf%ref
        wfo%o => wf%o
      end function

      subroutine glean_wf(wf)
!doc$ subroutine glean(wf)
        type(wavefunctions_td_obj) :: wf

!cod$
        if (wf%o%ref < 1) then
          select case (wf%o%usage)
          case (NORMAL)
             call glean(thy(wf%o%tp))
             call glean(thy(wf%o%mv))
             if (associated( wf%o%exps )) deallocate( wf%o%exps )
             call glean(thy(wf%o%hk))
          case (AUXILIARY)
             call glean(thy(wf%o%mv))
          end select
          deallocate( wf%o )
        end if
      end subroutine

      subroutine bequeath_wf(wf)
!doc$ subroutine bequeath(wf)
        type(wavefunctions_td_obj) :: wf

!cod$
        continue
      end subroutine

      subroutine assign_wf(wf,wf2)
!doc$ subroutine assignment(=)(wf,wf2)
        type(wavefunctions_td_obj), intent(inout) :: wf
        type(wavefunctions_td_obj), intent(in) :: wf2

!cod$
        type(wavefunctions_td_obj) :: wft
        call my(wf2)
        wft%o => wf%o
        wf%o%ref = wf%o%ref - wf%ref
        wf%o => wf2%o
        wf%o%ref = wf%o%ref + wf%ref
        call glean(wft)
        call glean(thy(wf2))
      end subroutine

      function wf_ref(wf) result(r)
!doc$ function x_ref(wf) result(r)
        type(wavefunctions_td_obj) :: wf
        integer, dimension(2) :: r
!       effects: Returns wf%ref and wf%o%ref.

!cod$
        r(1) = wf%ref
        r(2) = wf%o%ref
        call glean(wf)
      end function

      function wf_ghost(wf) result(g)
!doc$ function x_ghost(wf) result(g)
        type(wavefunctions_td_obj) :: wf
        type(ghost) :: g
!       effects: Returns wf%o%g.

!cod$
        call my(wf)
        g = wf%o%g
        call glean(thy(wf))
      end function

      function wf_multivector(wf) result(mv)
!doc$ function x_multivector(wf) result(mv)
        type(wavefunctions_td_obj) :: wf
        type(multivector_obj) :: mv
!       effects: Returns wf%o%mv.

!cod$
        call my(wf)
        call my(wf%o%mv,mv)
        call bequeath(thy(mv))
        call glean(thy(wf))
      end function

      function wf_h_kpoint(wf) result(hk)
!doc$ function x_h_kpoint(wf) result(hk)
        type(wavefunctions_td_obj) :: wf
        type(h_kpoint_obj) :: hk
!       requires: wf%o%usage = NORMAL
!       effects: Returns wf%o%hk.

!cod$
        call my(wf)
        call my(wf%o%hk,hk)
        call bequeath(thy(hk))
        call glean(thy(wf))
      end function

      function wf_n_gvectors(wf) result(n)
!doc$ function x_n_gvectors(wf) result(n)
        type(wavefunctions_td_obj) :: wf
        integer :: n
!       effects: Returns the number of G vectors used in wf.

!cod$
        call my(wf)
        n = x_n_gvectors(wf%o%mv)
        call glean(thy(wf))
      end function

      function wf_expectation_value(wf,ib) result(ev)
!doc$ function x_expectation_value(wf,ib) result(ev)
        type(wavefunctions_td_obj) :: wf
        integer, intent(in) :: ib
        real(double) :: ev
!       requires: wf%o%usage = NORMAL
!       requires: 1 <= ib <= size(wf%o%exps)
!       effects: Returns wf%o%exps(ib).

!cod$
        call my(wf)
        ev = wf%o%exps(ib)
        call glean(thy(wf))
      end function

      function wf_expectation_values(wf) result(ev)
!doc$ function x_expectation_values(wf) result(ev)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(size(wf%o%exps)) :: ev
!       requires: wf%o%usage = NORMAL
!       effects: Returns wf%o%exps.

!cod$
        call my(wf)
        ev = wf%o%exps
        call glean(thy(wf))
      end function

      function wf_time_propagator(wf) result(tp)
!doc$ function x_time_propagator(wf) result(tp)
        type(wavefunctions_td_obj) :: wf
        type(time_propagator_obj) :: tp
!       effects: Returns wf%o%tp.

!cod$
        call my(wf)
        call my(wf%o%tp,tp)
        call bequeath(thy(tp))
        call glean(thy(wf))
      end function

      function kinetic_energy_wf(wf,wts) result(ke)
!doc$ function kinetic_energy(wf,wts) result(ke)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(:), intent(in) :: wts
        real(double) :: ke
!       requires: wf%o%usage = NORMAL
!       effects: Returns the kinetic energy of wf.

!cod$
        call my(wf)
        ke = kinetic_energy(wf%o%mv,wf%o%hk,wts)
        call glean(thy(wf))
        if (error("Exit wavefunctions_td_mod::kinetic_energy_wf")) continue
      end function

      subroutine forces_wf(wf,wts,fwf)
!doc$ subroutine forces(wf,wts,fwf)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(:) :: wts
        real(double), dimension(:,:), intent(out) :: fwf
!       requires: wf%o%usage = NORMAL
!       modifies: fwf
!       effects: Returns atomic force contributions due to wavefunctions.
!       errors: Passes errors.

!cod$
        call my(wf)
        call forces(wf%o%mv,hamiltonian(wf%o%hk),wts,wf%o%exps,fwf) ; if (error()) goto 100
100     call glean(thy(wf))
        if (error("Exit wavefunctions_td_mod::forces_wf")) continue
      end subroutine

      subroutine pressure_wf(wf,wts,p)
!doc$ subroutine pressure(wf,wts,p)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(:), intent(in) :: wts
        real(double), intent(out) :: p
!       requires: wf%o%usage = NORMAL
!       effects: Returns pressure contributions due to wavefunctions.
!       errors: Passes errors.

!cod$
        call my(wf)
        call pressure(wf%o%mv,hamiltonian(wf%o%hk),wts,p) ; if (error()) goto 100
100     call glean(thy(wf))
        if (error("Exit wavefunctions_td_mod::pressure_wf")) continue
      end subroutine

      subroutine stress_tensor_wf(wf,wts,s)
!doc$ subroutine stress_tensor(wf,wts,s)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(:), intent(in) :: wts
        real(double), dimension(:,:), intent(out) :: s
!       requires: s be dimension(3,3). wf%o%usage = NORMAL
!       modifies: s
!       effects: Returns stress tensor contributions due to wavefunctions.
!       errors: Passes errors.

!cod$
        call my(wf)
        call stress_tensor(wf%o%mv,hamiltonian(wf%o%hk),wts,s) ; if (error()) goto 100
100     call glean(thy(wf))
        if (error("Exit wavefunctions_td_mod::stress_tensor_wf")) continue
      end subroutine

      subroutine add_density_wf(wf,weights,den)
!doc$ subroutine add_density(wf,weights,den)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(:), intent(in) :: weights
        type(gen_density_obj) :: den
!       requires: wf%o%usage = NORMAL
!       modifies: den
!       effects: Adds wf band densities weighted by weights to den.

!cod$
        call my(wf)
        call my(den)
        call add_density(hamiltonian(wf%o%hk),wf%o%mv,weights,den) ; if (error()) goto 100
100     call glean(thy(wf))
        call glean(thy(den))
        if (error("Exit wavefunctions_td_mod::add_density_wf")) continue
      end subroutine

      function get_norm_wf(wf) result(norm)
!doc$ function get_norm(wf) result(norm)
        type(wavefunctions_td_obj) :: wf
        real(double) :: norm
!       requires: wf%o%usage = NORMAL
!       effects: Returns the largest normalization discrepancy in wf.
!       errors: Passes errors.

!cod$
        character(line_len) :: init
        integer :: ib, nb
        real(double) :: err, tmp_err, tmp_norm
        real(double), dimension(:), allocatable :: norms
        type(multibasis_obj) :: mb
        type(multivector_obj) :: ov

        call my(wf)
        call my(x_multibasis(wf%o%mv),mb)

        nb = x_n_bands(wf%o%mv)
        allocate( norms(nb) )

        init = "zeros"
        call my(multivector(mb,init=init),ov)

        call apply_overlap(wf%o%mv,ov,hamiltonian(wf%o%hk)) ! O*v --> ov
        
        call multiply(wf%o%mv,ov,norms)
        norm = norms(1)
        err = abs(norm - 1.0_double)
        do ib = 2,nb
          tmp_norm = norms(ib)
          tmp_err = abs(norms(ib) - 1.0_double)
          if (err < tmp_err) then
            norm = tmp_norm
            err = tmp_err
          end if
        end do

        call glean(thy(ov))
        call glean(thy(mb))
        
100     if (allocated( norms )) deallocate( norms )

        call glean(thy(wf))

        if (error("Exit wavefunctions_mod::get_norm_wf")) continue

      end function

      subroutine decompose_wf(wf,site_data,mode,rsa,b1)
!doc$ subroutine decompose(wf,site_data,mode,rsa)
        type(wavefunctions_td_obj) :: wf
        real(double), dimension(:,:), intent(in) :: site_data
        character(line_len), intent(in) :: mode
        real(double), dimension(:,:,:), intent(out) :: rsa
        integer :: b1
!       requires: wf%o%usage = NORMAL
!       modifies: rsa
!       effects: Returns the decomposition of the Kohn-Sham functions into s, p, & d spherical harmonics.
!       errors: Passes errors.

!cod$
        call my(wf)
        call decompose(wf%o%mv,site_data,mode,rsa,b1) ; if (error()) goto 100
100     call glean(thy(wf))
        if (error("Exit wavefunctions_td_mod::decompose_wf")) continue

      end subroutine

      subroutine diary_wf(wf)
!doc$ subroutine diary(wf)
        type(wavefunctions_td_obj) :: wf
!       requires: wf%o%usage = 	NORMAL
!       effects: Writes wf information to the diary.

!cod$
        call my(wf)
        call diary(wf%o%mv)
!       Ryan: reconcile this
!        call diary(wf%o%tp)
        call glean(thy(wf))
      end subroutine


      subroutine write_restart_wf(wf,nrestf)
!doc$ subroutine write_restart(wf,nrestf)
        type(wavefunctions_td_obj) :: wf
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes restart information to nrestf.
!       errors: Passes errors.

!cod$
        integer(long) :: dsize, ios, ndata, s4
        

        call my(wf)
        call my(nrestf)

        if (i_access(nrestf)) call startblock(nrestf,"WAVEFUNCTIONS_TD")

        !** Time propagator object
        call write_restart(wf%o%tp,nrestf)

        !** Expectation values
        if (i_access(nrestf)) then
           call writetag(nrestf,"EXPECTATIONVALUES")
           dsize = sizeof_double ; ndata = size(wf%o%exps)
           call writef(wf%o%exps,dsize,ndata,x_tagfd(nrestf),ios)
        end if

        call write_restart(wf%o%mv,nrestf) ; if (error()) goto 100

        if (i_access(nrestf)) call endblock(nrestf) 

100     call glean(thy(wf))
        call glean(thy(nrestf))

        if (error("Exit wavefunctions_td_mod::write_restart_wf")) continue

      end subroutine

! local routines

      subroutine own_i(wf)
        type(wavefunctions_td_obj) :: wf, wft
        if (wf%ref < wf%o%ref) then
          call warn("WARNING: wavefunctions mutation called: massive copy going on")
          allocate( wft%o )
          wft%o%ref = 0
          wft%o%g = wf%o%g
!          wft%o%g_crystal = wf%o%g_crystal
          wft%o%usage = wf%o%usage
          select case (wf%o%usage)
          case (NORMAL)
             call my(wf%o%tp,wft%o%tp)
             call my(wf%o%mv,wft%o%mv)
             allocate( wft%o%exps(size(wf%o%exps)) )
             wft%o%exps = wf%o%exps
             call my(wf%o%hk,wft%o%hk)
          case (AUXILIARY)
             call my(wf%o%mv,wft%o%mv)
          end select
          wf%o%ref = wf%o%ref - wf%ref
          wf%o => wft%o
          wf%o%ref = wf%o%ref + wf%ref
        end if
      end subroutine

      end module

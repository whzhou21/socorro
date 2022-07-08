!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module wavefunctions_es_mod
!doc$ module wavefunctions_es_mod

!     One datatype is available here: type(wavefunctions_es_obj).

!     wavefunctions_es_mod creates and maintains a set of wavefunctions in calculations involving eigensolutions.

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
      use eigensolver_mod
      use ghost_mod
      use gen_density_mod
      use dyad_kpoint_mod

!cod$
      implicit none
      private
   
      ! usage
      integer, parameter :: NORMAL    = 1
      integer, parameter :: AUXILIARY = 2

      type, public :: wavefunctions_es_rep
        integer :: ref
        type(ghost) :: g
        type(ghost) :: g_crystal                        ! crystal ghost
        integer :: usage                                ! usage wrt kgroup
        type(eigensolver_obj) :: es                     ! eigensolver
        type(multivector_obj) :: mv                     ! multivector
        real(double), dimension(:), pointer :: eigs     ! eigenvalues
        real(double) :: res_norm                        ! mv residual norm
        type(h_kpoint_obj) :: hk                        ! k-point dependent part of the hamiltonian
      end type

      type, public :: wavefunctions_es_obj
        private
        integer :: ref
        type(wavefunctions_es_rep), pointer :: o
      end type

!doc$
      public :: wavefunctions_es
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: wormhole
      public :: x_ref
      public :: x_ghost
      public :: x_multivector
      public :: x_h_kpoint
      public :: x_n_gvectors
      public :: eigenvalues
      public :: residual_norm
      public :: kinetic_energy
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: add_density
      public :: decompose
      public :: distribute
      public :: release
      public :: diary
      public :: write_restart

!cod$
      interface wavefunctions_es
        module procedure constructor_wf
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
      interface wormhole
        module procedure wormhole_wf
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
      interface eigenvalues
        module procedure eigenvalues_wf
      end interface
      interface residual_norm
        module procedure residual_norm_wf
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
      interface decompose
        module procedure decompose_wf
      end interface
      interface distribute
        module procedure distribute_wf
      end interface
      interface release
        module procedure release_wf
      end interface
      interface diary
        module procedure diary_wf
      end interface
      interface write_restart
        module procedure write_restart_wf
      end interface

      contains

! public routines

      function constructor_wf(hc,mb,dk,restf) result(wf)
!doc$ function wavefunctions_es(hc,mb,dk,restf) result(wf)
        type(h_common_obj) :: hc
        type(multibasis_obj) :: mb
        type(dyad_kpoint_obj), optional :: dk
        type(tagio_obj), optional :: restf
        type(wavefunctions_es_obj) :: wf
!       requires: usage be either "normal" or "auxiliary".
!       effects: Constructs a new wf.
!       errors: Error on entry. Passes errors.

!cod$
        logical :: found
        character(line_len) :: init

        if (error("  Error on entry")) then
          wf%ref = 0
          allocate( wf%o )
          wf%o%ref = 0
          goto 999
        end if

        call my(hc)
        call my(mb)
        if (present(dk)) call my(dk)
        if (present(restf)) call my(restf)

        wf%ref = 0
        allocate( wf%o )
        wf%o%ref = 0
        wf%o%g = x_ghost()

        wf%o%usage = x_usage(mb)

        if (present(restf)) then
          select case (wf%o%usage)
          case (NORMAL)
            wf%o%g_crystal = x_ghost(x_crystal(hc))
            call my(eigensolver(),wf%o%es)                                                 ; if (error()) goto 100
            call my(multivector(mb,restf=restf),wf%o%mv)                                   ; if (error()) goto 100
            call my(kpoint_hamiltonian(hc,mb),wf%o%hk)                                     ; if (error()) goto 100
            call orthonormalize(wf%o%mv,hamiltonian(wf%o%hk))                              ; if (error()) goto 100
            allocate( wf%o%eigs(x_n_bands(mb)) )
            call eigensolve(wf%o%es,hamiltonian(wf%o%hk),wf%o%mv,wf%o%eigs,wf%o%res_norm)  ; if (error()) goto 100
          case (AUXILIARY)
            call my(multivector(mb,restf=restf),wf%o%mv)                                   ; if (error()) goto 100
          end select
        else
          select case (wf%o%usage)
          case (NORMAL)
            wf%o%g_crystal = x_ghost(x_crystal(hc))
            call my(eigensolver(),wf%o%es)                                                 ; if (error()) goto 100
            call arglc("wf_init",init,found)
            if (.not.found) init = "random"
            select case (trim(init))
            case ("random","diagnostic")
              continue
            case default
              if (error(.true.,"ERROR: wf_init was not recognized")) goto 100
            end select
            call my(multivector(mb,init),wf%o%mv)                                          ; if (error()) goto 100
            if (present(dk)) then
               call my(kpoint_hamiltonian(hc,mb,dk),wf%o%hk)
            else
               call my(kpoint_hamiltonian(hc,mb),wf%o%hk)                                  ; if (error()) goto 100
            end if
            call orthonormalize(wf%o%mv,hamiltonian(wf%o%hk))                              ; if (error()) goto 100
            allocate( wf%o%eigs(x_n_bands(mb)) )
            call eigensolve(wf%o%es,hamiltonian(wf%o%hk),wf%o%mv,wf%o%eigs,wf%o%res_norm)  ; if (error()) goto 100
          case (AUXILIARY)
            call my(multivector(mb),wf%o%mv)                                               ; if (error()) goto 100
          end select
        end if

100     call glean(thy(hc))
        call glean(thy(mb))
        if (present(dk)) call glean(thy(dk))
        if (present(restf)) call glean(thy(restf))

999     if (error("Exit wavefunctions_es_mod::constructor_wf")) continue

      end function

      subroutine update_wf(wf,hc,mb,dk)
!doc$ subroutine update(wf,hc,mb,dk)
        type(wavefunctions_es_obj) :: wf
        type(h_common_obj) :: hc
        type(multibasis_obj), optional :: mb
        type(dyad_kpoint_obj), optional :: dk

!       modifies: wf
!       effects: Updates wf with respect to its dependencies.
!       errors: mb change. Passes errors.

!cod$
        logical :: crystal_change, mb_change

        call my(wf)
        call my(hc)
        if (present(mb)) call my(mb)
        if (present(dk)) call my(dk)

        crystal_change = ( x_ghost(x_crystal(hc)) /= wf%o%g_crystal )
        mb_change = .false.
        if (present(mb)) then
          mb_change = ( x_ghost(x_multibasis(wf%o%mv)) /= x_ghost(mb) )
          if (error(mb_change,"ERROR: multibasis change is not currently allowed")) goto 100
        end if

        if (wf%o%usage == AUXILIARY) goto 100

        call own_i(wf)
        wf%o%g = x_ghost()

        if (present(dk)) then
           call update(wf%o%hk,hc,dk = dk) ; if (error()) goto 100
        else
           call update(wf%o%hk,hc) ; if (error()) goto 100
        end if
        if (crystal_change) then
          wf%o%g_crystal = x_ghost(x_crystal(hc))
          call update(wf%o%es)
          if (.not.overlap_is_identity(hc)) call orthonormalize(wf%o%mv,hamiltonian(wf%o%hk)) ; if (error()) goto 100
        end if
        call eigensolve(wf%o%es,hamiltonian(wf%o%hk),wf%o%mv,wf%o%eigs,wf%o%res_norm) ; if (error()) goto 100

100     call glean(thy(wf))
        call glean(thy(hc))
        if (present(mb)) call glean(thy(mb))
        if (present(dk)) call glean(thy(dk))

        if (error("Exit wavefunctions_es_mod::update_wf")) continue

      end subroutine

      subroutine my_wf(wf)
!doc$ subroutine my(wf)
        type(wavefunctions_es_obj) :: wf

!cod$
        wf%ref = wf%ref + 1
        wf%o%ref = wf%o%ref + 1
      end subroutine

      subroutine my_new_wf(wfi,wf)
!doc$ subroutine my(wfi,wf)
        type(wavefunctions_es_obj) :: wfi, wf

!cod$
        wf%ref = 1
        wf%o => wfi%o
        wf%o%ref = wf%o%ref + 1
      end subroutine

      function thy_wf(wf) result(wfo)
!doc$ function thy(wf) result(wfo)
        type(wavefunctions_es_obj) :: wf, wfo

!cod$
        wf%ref = wf%ref - 1
        wf%o%ref = wf%o%ref - 1
        wfo%ref = wf%ref
        wfo%o => wf%o
      end function

      subroutine glean_wf(wf)
!doc$ subroutine glean(wf)
        type(wavefunctions_es_obj) :: wf

!cod$
        if (wf%o%ref < 1) then
          select case (wf%o%usage)
          case (NORMAL)
            call glean(thy(wf%o%es))
            call glean(thy(wf%o%mv))
            if (associated( wf%o%eigs )) deallocate( wf%o%eigs )
            call glean(thy(wf%o%hk))
          case (AUXILIARY)
            call glean(thy(wf%o%mv))
          end select
          deallocate( wf%o )
        end if
      end subroutine

      subroutine bequeath_wf(wf)
!doc$ subroutine bequeath(wf)
        type(wavefunctions_es_obj) :: wf

!cod$
        continue
      end subroutine

      subroutine assign_wf(wf,wf2)
!doc$ subroutine assignment(=)(wf,wf2)
        type(wavefunctions_es_obj), intent(inout) :: wf
        type(wavefunctions_es_obj), intent(in) :: wf2

!cod$
        type(wavefunctions_es_obj) :: wft
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
        type(wavefunctions_es_obj) :: wf
        integer, dimension(2) :: r
!       effects: Returns wf%ref and wf%o%ref.

!cod$
        r(1) = wf%ref
        r(2) = wf%o%ref
        call glean(wf)
      end function

      function wormhole_wf(wf) result(wfr)
!doc$ function wormhole(wf) result(wfr)
        type(wavefunctions_es_obj) :: wf
        type(wavefunctions_es_rep), pointer :: wfr
!       effects: Points wfr at wf%o
!       errors: Wormhole is an implementation dependent thing and you should know what you're doing.

!cod$
        call my(wf)
        wfr => wf%o
        call glean(thy(wf))
      end function

      function wf_ghost(wf) result(g)
!doc$ function x_ghost(wf) result(g)
        type(wavefunctions_es_obj) :: wf
        type(ghost) :: g
!       effects: Returns wf%o%g.

!cod$
        call my(wf)
        g = wf%o%g
        call glean(thy(wf))
      end function

      function wf_multivector(wf) result(mv)
!doc$ function x_multivector(wf) result(mv)
        type(wavefunctions_es_obj) :: wf
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
        type(wavefunctions_es_obj) :: wf
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
        type(wavefunctions_es_obj) :: wf
        integer :: n
!       effects: Returns the number of G vectors used in wf.

!cod$
        call my(wf)
        n = x_n_gvectors(wf%o%mv)
        call glean(thy(wf))
      end function

      subroutine eigenvalues_wf(wf,ev)
!doc$ subroutine eigenvalues(wf,ev)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:) :: ev
!       effects: Returns the eigenvalues of a NORMAL usage wf and 0 otherwise.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          ev = wf%o%eigs
        case (AUXILIARY)
          ev = 0.0_double
        end select
        call glean(thy(wf))
      end subroutine

      subroutine residual_norm_wf(wf,rn)
!doc$ subroutine residual_norm(wf,rn)
        type(wavefunctions_es_obj) :: wf
        real(double) :: rn
!       effects: Returns the residual norm of a NORMAL usage wf and 0 otherwise.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          rn = wf%o%res_norm
        case (AUXILIARY)
          rn = 0.0_double
        end select
        call glean(thy(wf))
      end subroutine

      function kinetic_energy_wf(wf,wts) result(ke)
!doc$ function kinetic_energy(wf,wts) result(ke)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:), intent(in) :: wts
        real(double) :: ke
!       effects: Returns the kinetic energy of a NORMAL usage wf and 0 otherwise.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          ke = kinetic_energy(wf%o%mv,wf%o%hk,wts)
        case (AUXILIARY)
          ke = 0.0_double
        end select
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::kinetic_energy_wf")) continue
      end function

      subroutine forces_wf(wf,wts,f)
!doc$ subroutine forces(wf,wts,f)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:) :: wts
        real(double), dimension(:,:), intent(out) :: f
!       requires: f be dimension(3,number-of-atoms)
!       modifies: f
!       effects: Returns atomic force contributions due to wavefunctions.
!       errors: Passes errors.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          call forces(wf%o%mv,hamiltonian(wf%o%hk),wts,wf%o%eigs,f)
        case (AUXILIARY)
          f = 0.0_double
        end select
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::forces_wf")) continue
      end subroutine

      subroutine pressure_wf(wf,wts,p)
!doc$ subroutine pressure(wf,wts,p)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:), intent(in) :: wts
        real(double), intent(out) :: p
!       effects: Returns pressure contributions due to wavefunctions.
!       errors: Passes errors.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          call pressure(wf%o%mv,hamiltonian(wf%o%hk),wts,p)
        case (AUXILIARY)
          p = 0.0_double
        end select
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::pressure_wf")) continue
      end subroutine

      subroutine stress_tensor_wf(wf,wts,s)
!doc$ subroutine stress_tensor(wf,wts,s)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:), intent(in) :: wts
        real(double), dimension(:,:), intent(out) :: s
!       requires: s be dimension(3,3)
!       modifies: s
!       effects: Returns stress tensor contributions due to wavefunctions.
!       errors: Passes errors.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          call stress_tensor(wf%o%mv,hamiltonian(wf%o%hk),wts,s)
        case (AUXILIARY)
          s = 0.0_double
        end select
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::stress_tensor_wf")) continue
      end subroutine

      subroutine add_density_wf(wf,weights,den)
!doc$ subroutine add_density(wf,weights,den)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:), intent(in) :: weights
        type(gen_density_obj) :: den
!       modifies: den
!       effects: Adds wf band densities weighted by weights to den.

!cod$
        call my(wf)
        call my(den)
        select case (wf%o%usage)
        case (NORMAL)
          call add_density(hamiltonian(wf%o%hk),wf%o%mv,weights,den) ; if (error()) goto 100
        end select
100     call glean(thy(wf))
        call glean(thy(den))
        if (error("Exit wavefunctions_es_mod::add_density_wf")) continue
      end subroutine

      subroutine decompose_wf(wf,site_data,mode,rsa,b1)
!doc$ subroutine decompose(wf,site_data,mode,rsa)
        type(wavefunctions_es_obj) :: wf
        real(double), dimension(:,:), intent(in) :: site_data
        character(line_len), intent(in) :: mode
        real(double), dimension(:,:,:), intent(out) :: rsa
        integer :: b1
!       requires: rsa be dimension(9,*)
!       modifies: rsa
!       effects: Returns the decomposition of the Kohn-Sham functions into s, p, & d spherical harmonics.
!       errors: Passes errors.

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (NORMAL)
          call decompose(wf%o%mv,site_data,mode,rsa,b1)
        case (AUXILIARY)
          rsa = 0.0_double
        end select
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::decompose_wf")) continue

      end subroutine

      subroutine distribute_wf(wf,rank)
!doc$ subroutine distribute(wf,rank)
        type(wavefunctions_es_obj) :: wf
        integer, intent(in) :: rank
!       requires: Must be followed by a call to release_wf with no intervening modification of wf.
!       effects: Distributes wf%o%mv data among all kgroups.
!       errors:

!cod$
        call my(wf)
        call distribute(wf%o%mv,rank)
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::distribute_wf")) continue

      end subroutine

      subroutine release_wf(wf)
!doc$ subroutine release(wf)
        type(wavefunctions_es_obj) :: wf
!       requires: Must be preceeded by a call to distribute_wf with no intervening modification of wf.
!       effects: Releases memory used by wf%o%mv.
!       errors:

!cod$
        call my(wf)
        select case (wf%o%usage)
        case (AUXILIARY)
          call release(wf%o%mv)
        end select
        call glean(thy(wf))
        if (error("Exit wavefunctions_es_mod::release_wf")) continue

      end subroutine

      subroutine diary_wf(wf)
!doc$ subroutine diary(wf)
        type(wavefunctions_es_obj) :: wf
!       requires: wf%o%usage = NORMAL.
!       effects: Writes wf information to the diary.

!cod$
        call my(wf)

        call diary(wf%o%mv)
        call diary(wf%o%es)

        call glean(thy(wf))

      end subroutine

      subroutine write_restart_wf(wf,nrestf)
!doc$ subroutine write_restart(wf,nrestf)
        type(wavefunctions_es_obj) :: wf
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes wf information to nrestf.
!       errors: Passes errors.

!cod$
        call my(wf)
        call my(nrestf)

        call write_restart(wf%o%mv,nrestf) ; if (error()) goto 100

100     call glean(thy(wf))
        call glean(thy(nrestf))

        if (error("Exit wavefunctions_es_mod::write_restart_wf")) continue

      end subroutine

! private routines

      subroutine own_i(wf)
        type(wavefunctions_es_obj) :: wf, wft
        if (wf%ref < wf%o%ref) then
          call warn("WARNING: wavefunctions mutation called: massive copy going on")
          allocate( wft%o )
          wft%o%ref = 0
          wft%o%g = wf%o%g
          wft%o%usage = wf%o%usage
          select case (wf%o%usage)
          case (NORMAL)
            wft%o%g_crystal = wf%o%g_crystal
            call my(wf%o%es,wft%o%es)
            call my(wf%o%mv,wft%o%mv)
            allocate( wft%o%eigs(size(wf%o%eigs)) )
            wft%o%eigs = wf%o%eigs
            wft%o%res_norm = wf%o%res_norm
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

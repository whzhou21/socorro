!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module paw_data_mod
!doc$ module paw_data_mod

!     One datatype is defined here: type(paw_data_obj). paw_data_obj contains information about paw
!     functions for a given atom type.

      use kind_mod
      use mpi_mod
      use error_mod
      use io_mod 
      use arg_mod
      use diary_mod
      use math_mod

!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = CONFIG

      integer, parameter :: LINEAR = 1
      integer, parameter :: LOG    = 2

      integer :: w_ok, w_eol, w_eof

      type :: paw_data_rep
        integer :: ref                                              ! reference count
        character(tag_sz) :: name                                   ! name of this data type
        real(double) :: valence_electrons                           ! number of valence electrons
        real(double) :: matching_radius                             ! matching radius for smooth basis functions
        character(tag_sz) :: xc_type                                ! exchange-correlation type
        real(double) :: mass                                        ! mass (a.u.)
        real(double), dimension(:), pointer :: r_sg                 ! standard grid: values
        real(double), dimension(:), pointer :: r2_sg                !                values^2
        real(double), dimension(:), pointer :: wt_sg                !                integration weights
        real(double), dimension(:), pointer :: r_cg                 ! coretail grid: values
        integer, dimension(:), pointer :: l_value                   ! basis functions: l values
        real(double), dimension(:), pointer :: occupation           ! basis functions: initial occupations
        integer, dimension(:), pointer :: nl_base                   !                  starting indices
        integer, dimension(:,:), pointer :: lut_nlm                 !                  n, l, and m
        real(double), dimension(:,:,:), pointer :: phir2ij          !                  phi_i*phi_j/r^2 
        real(double), dimension(:,:,:), pointer :: dphir2ij_dr      !                  d(phi_i*phi_j/r^2)/dr
        real(double), dimension(:,:,:), pointer :: tphir2ij         !                  phi~_i*phi~_j/r^2 
        real(double), dimension(:,:,:), pointer :: dtphir2ij_dr     !                  d(phi~_i*phi~_j/r^2)/dr
        integer :: theta_points                                     !                  angular grid theta parameter
        integer :: phi_points                                       !                  angular grid phi parameter
        real(double), dimension(:,:), pointer :: tp                 ! projector functions
        real(double) :: r_opt                                       ! real-space projectors: optimization radius
        real(double) :: g_icut                                      !                        g_vector inner cutoff
        real(double) :: g_ocut                                      !                        g-vector outer cutoff
        real(double), dimension(:), pointer :: g                    !                        g-vector grid
        real(double), dimension(:), pointer :: w_max                !                        maximum errors
        real(double), dimension(:,:), pointer :: tpfo               !                        tp optimized Fourier coefficients
        real(double), dimension(:), pointer :: kinetic              ! kinetic energy matrix
        real(double), dimension(:), pointer :: v_ion                ! v ion matrix
        real(double), dimension(:), pointer :: vlocal               ! local potential
        real(double), dimension(:,:), pointer :: oij                ! overlap matrix
        real(double), dimension(:,:), pointer :: hatden             ! hat density grid values
        real(double), dimension(:), pointer :: hat_selfenergy       ! hat self-energy 
        real(double), dimension(:,:), pointer :: denmat             ! density matrix grid values
        integer, dimension(:), pointer :: lut_denvhat               ! look-up-table for the density and hat potential
        real(double), dimension(:), pointer :: cijkl                ! coefficients of w^a_kl
        integer, dimension(:), pointer :: lut_cijkl                 ! encoded indices for cij
        real(double), dimension(:), pointer :: aqlm                 ! G^(LM)_limiljmj,lkmkllml
        real(double), dimension(:), pointer :: avlm                 ! avlm
        integer, dimension(:), pointer :: lut_orbital               ! encoded indices for aqlm and avlm
        real(double) :: qeffion                                     ! effective ion charge
        real(double), dimension(:), pointer :: core_density         ! core density: values
        real(double), dimension(:), pointer :: core_grad            !               gradient magnitude
        real(double), dimension(:), pointer :: coretail_density_cg  ! coretail density: values (coretail grid)
        real(double), dimension(:), pointer :: coretail_density_sg  !                   values (standard grid)
        real(double), dimension(:), pointer :: coretail_grad_sg     !                   gradient magnitude (standard grid)
        real(double) :: coretail_selfenergy                         !                   self-energy
        real(double) :: coretail_hatenergy                          !                   interaction energy with hat density
      end type 

      type, public :: paw_data_obj
        private
        integer :: ref
        type(paw_data_rep), pointer :: o
      end type

      type :: word_context
        character :: comment, literal
        character(30) :: delimiter
        character(200) :: text
        integer :: column, length, error, last_error
        type(file_obj) :: f
      end type

!doc$
      public :: paw_data
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_name
      public :: x_valence_electrons
      public :: x_matching_radius
      public :: x_xc_type
      public :: x_mass
      public :: x_grid_size
      public :: type_r2
      public :: type_wt
      public :: x_l_value
      public :: x_occupation
      public :: x_nl_base
      public :: x_nlm_size
      public :: x_nlm
      public :: x_basis_size
      public :: type_phir2ij
      public :: type_dphir2ij_dr
      public :: type_tphir2ij
      public :: type_dtphir2ij_dr
      public :: x_theta_points
      public :: x_phi_points
      public :: projector_f_value
      public :: projector_f_values
      public :: x_radius
      public :: x_inner_cutoff
      public :: x_outer_cutoff
      public :: x_error
      public :: projector_r_value
      public :: projector_r_gradients
      public :: diary_rs_projectors
      public :: x_kinetic
      public :: x_v_ion
      public :: x_oij
      public :: local_f_value
      public :: local_f_values
      public :: x_hatden_size
      public :: hatden_f_value
      public :: hatden_f_values
      public :: x_hat_self_energy
      public :: x_denmat_size
      public :: denmat_f_value
      public :: denmat_f_values
      public :: denvhat_decode
      public :: gaunt_complex
      public :: x_cijkl_size
      public :: x_cijkl
      public :: cijkl_decode
      public :: x_qvlm_size
      public :: x_aqlm
      public :: x_avlm
      public :: orbital_decode
      public :: x_qeffion
      public :: type_core_density
      public :: type_core_grad
      public :: type_coretail_density
      public :: type_coretail_grad
      public :: coretail_f_value
      public :: coretail_f_values
      public :: x_coretail_hatenergy
      public :: x_coretail_selfenergy

!cod$
      interface paw_data
        module procedure constructor_pd
      end interface
      interface my
        module procedure my_pd, my_new_pd
      end interface
      interface thy
        module procedure thy_pd
      end interface
      interface glean
        module procedure glean_pd
      end interface
      interface bequeath
        module procedure bequeath_pd
      end interface
      interface assignment(=)
        module procedure assign_pd
      end interface
      interface x_ref
        module procedure pd_ref
      end interface
      interface x_name
        module procedure pd_name
      end interface
      interface x_valence_electrons
        module procedure pd_valence_electrons
      end interface
      interface x_matching_radius
        module procedure pd_matching_radius
      end interface
      interface x_xc_type
        module procedure pd_xc_type
      end interface
      interface x_mass
        module procedure pd_mass
      end interface
      interface x_grid_size
        module procedure pd_grid_size
      end interface
      interface x_l_value
        module procedure pd_l_value
      end interface
      interface x_occupation
        module procedure pd_occupation
      end interface
      interface x_nl_base
        module procedure pd_nl_base
      end interface
      interface x_nlm_size
        module procedure pd_nlm_size
      end interface
      interface x_nlm
        module procedure pd_nlm
      end interface
      interface x_basis_size
        module procedure pd_basis_size
      end interface
      interface x_theta_points
        module procedure pd_theta_points
      end interface
      interface x_phi_points
        module procedure pd_phi_points
      end interface
      interface projector_f_value
        module procedure pd_projector_f_value
      end interface
      interface projector_f_values
        module procedure pd_projector_f_values
      end interface
      interface x_radius
        module procedure pd_radius
      end interface
      interface x_inner_cutoff
        module procedure pd_inner_cutoff
      end interface
      interface x_outer_cutoff
        module procedure pd_outer_cutoff
      end interface
      interface x_error
        module procedure pd_error
      end interface
      interface projector_r_value
        module procedure projector_r_value_pd
      end interface
      interface projector_r_gradients
        module procedure projector_r_gradients_pd
      end interface
      interface diary_rs_projectors
        module procedure diary_rs_projectors_pd
      end interface
      interface x_kinetic
        module procedure pd_kinetic
      end interface
      interface x_v_ion
        module procedure pd_v_ion
      end interface
      interface x_oij
        module procedure pd_oij
      end interface
      interface local_f_value
        module procedure local_f_value_pd
      end interface
      interface local_f_values
        module procedure local_f_values_pd
      end interface
      interface x_hatden_size
        module procedure pd_hatden_size
      end interface
      interface hatden_f_value
        module procedure hatden_f_value_pd
      end interface
      interface hatden_f_values
        module procedure hatden_f_values_pd
      end interface
      interface x_hat_self_energy
        module procedure pd_hat_self_energy
      end interface
      interface x_denmat_size
        module procedure pd_denmat_size
      end interface
      interface denmat_f_value
        module procedure denmat_f_value_pd
      end interface
      interface denmat_f_values
        module procedure denmat_f_values_pd
      end interface
      interface x_cijkl_size
        module procedure pd_cijkl_size
      end interface
      interface x_cijkl
        module procedure pd_cijkl
      end interface
      interface x_qvlm_size
        module procedure pd_qvlm_size
      end interface
      interface x_aqlm
        module procedure pd_aqlm
      end interface
      interface x_avlm
        module procedure pd_avlm
      end interface
      interface x_qeffion
        module procedure pd_qeffion
      end interface
      interface coretail_f_value
        module procedure coretail_f_value_pd
      end interface
      interface coretail_f_values
        module procedure coretail_f_values_pd
      end interface
      interface x_coretail_hatenergy
        module procedure pd_coretail_hatenergy
      end interface
      interface x_coretail_selfenergy
        module procedure pd_coretail_selfenergy
      end interface

      contains

      function constructor_pd(f) result(pd)
!doc$ function paw_data(f) result(pd)
        type(file_obj) :: f   
        type(paw_data_obj) :: pd
!       requires: f be open and pointing to paw data.
!       effects: Creates a new pd.
!       errors: Format errors.

!cod$
        logical :: found
        character(tag_sz) :: shape_type
        character(17+tag_sz) :: pr_tag
        character(line_len) :: token
        integer :: base_i, base_j, i, j, k, l, m, n, nili, njlj, grid_type
        integer :: basis_size, denvhat_size, hartree_size, lcao_size, grid_size_sg, grid_size_cg, nlm_size, overlap_size
        integer, dimension(:), pointer :: lut_v_hartree
        real(double) :: atom_energy, total_electrons, core_electrons, log_grid_pf, grid_step, lcao_step, pf
        real(double), dimension(:), allocatable :: r2i_sg, drdu_sg, r2_cg, r2i_cg, wt_cg, drdu_cg
        real(double), dimension(:), allocatable :: density, v_hat
        real(double), dimension(:), pointer :: input_hatden, input_oij, r1d, v_hartree
        real(double), dimension(:,:), allocatable :: phi, tphi
        real(double), dimension(:,:), pointer :: r2d
        type(word_context) :: wc   

        call my(f)

        pd%ref = 0
        allocate( pd%o )
        pd%o%ref = 0

        nullify( lut_v_hartree, input_hatden, input_oij, r1d, v_hartree, r2d )

        nullify( pd%o%r_sg )
        nullify( pd%o%r2_sg )
        nullify( pd%o%wt_sg )
        nullify( pd%o%r_cg )
        nullify( pd%o%l_value )
        nullify( pd%o%nl_base )
        nullify( pd%o%lut_nlm )
        nullify( pd%o%phir2ij )
        nullify( pd%o%dphir2ij_dr )
        nullify( pd%o%tphir2ij )
        nullify( pd%o%dtphir2ij_dr )
        nullify( pd%o%tp )
        nullify( pd%o%w_max )
        nullify( pd%o%g )
        nullify( pd%o%tpfo )
        nullify( pd%o%kinetic )
        nullify( pd%o%v_ion )
        nullify( pd%o%vlocal )
        nullify( pd%o%oij )
        nullify( pd%o%hatden )
        nullify( pd%o%hat_selfenergy )
        nullify( pd%o%denmat )
        nullify( pd%o%lut_denvhat )
        nullify( pd%o%cijkl )
        nullify( pd%o%lut_cijkl )
        nullify( pd%o%aqlm )
        nullify( pd%o%avlm )
        nullify( pd%o%lut_orbital )
        nullify( pd%o%core_density )
        nullify( pd%o%core_grad )
        nullify( pd%o%coretail_density_cg )
        nullify( pd%o%coretail_density_sg )
        nullify( pd%o%coretail_grad_sg )

        grid_type = LINEAR

        call open_word_i(f,wc) ; if (error()) goto 100

        call get_next_word_i(wc,token)
        do
          if (token == "ATOMTYPE") then
            call get_next_word_i(wc,pd%o%name)
            if (error(len_trim(pd%o%name) > 6,"ERROR: Atom name has more than six characters")) goto 100
          elseif (token == "MASS") then
            call get_real_i(wc,pd%o%mass) 
          elseif (token == "ATOMXCTYPE") then
            call get_next_word_i(wc,pd%o%xc_type)
          elseif (token == "ATOMIC_CHARGE") then
            call get_integer_i(wc,i)
            total_electrons = real(i,double)
          elseif (token == "CORE_CHARGE") then
            call get_real_i(wc,core_electrons)
          elseif (token == "RC") then
            call get_real_i(wc,pd%o%matching_radius)
          elseif (token == "SHAPE_TYPE") then
            call get_next_word_i(wc,shape_type)
          elseif (token == "BASIS_SIZE") then
            call get_integer_i(wc,basis_size)
          elseif (token == "ORBITALS") then
            allocate( pd%o%l_value(basis_size) )
            call get_integers_i(wc,pd%o%l_value)
          elseif (token == "INITOCC") then
            allocate( pd%o%occupation(basis_size) )
            call get_reals_i(wc,pd%o%occupation)
          elseif (token == "MESH_SIZE") then
            call get_integer_i(wc,grid_size_sg) 
          elseif (token == "MESH_STEP") then
            call get_real_i(wc,grid_step)
          else if (token == "LOG_GRID") then
            grid_type = LOG
            call get_real_i(wc,log_grid_pf)
          elseif (token == "CORETAIL_POINTS") then
            call get_integer_i(wc,grid_size_cg) 
          elseif (token == "LCAO_SIZE") then
            call get_integer_i(wc,lcao_size) 
          elseif (token == "LCAO_STEP") then
            call get_real_i(wc,lcao_step)
          elseif (token == "CORE_DENSITY") then
            allocate( pd%o%core_density(grid_size_sg) )
            call get_reals_i(wc,pd%o%core_density)
          elseif (token == "CORETAIL_DENSITY") then
            allocate( pd%o%coretail_density_cg(grid_size_cg) )
            call get_reals_i(wc,pd%o%coretail_density_cg)
          elseif (token == "SHAPE_FUNC") then
            allocate( input_hatden(grid_size_sg) )
            call get_reals_i(wc,input_hatden)
          elseif (token == "PSEUDO_VALENCE_DENSITY") then
            call get_integer_i(wc,j)
            allocate( r1d(j) )
            call get_reals_i(wc,r1d,l=.false.)
            deallocate( r1d )
          elseif (token == "VLOCFUN") then
            allocate( pd%o%vlocal(grid_size_sg) )
            call get_reals_i(wc,pd%o%vlocal)
          elseif (token == "VLOCION") then
            call get_integer_i(wc,j)
            allocate( r1d(j) )
            call get_reals_i(wc,r1d,l=.false.)
            deallocate( r1d )
          elseif (token == "VLOCION_NOHAT") then
            call get_integer_i(wc,j)
            allocate( r1d(j) )
            call get_reals_i(wc,r1d,l=.false.)
            deallocate( r1d )
          elseif (token == "TPROJECTOR") then
            call get_integer_i(wc,j)
            if (.not.associated( pd%o%tp )) allocate( pd%o%tp(grid_size_sg,basis_size) )
            call get_reals_i(wc,pd%o%tp(:,j))    
          elseif (token == "PHI") then
            call get_integer_i(wc,j)
            if (.not.allocated( phi )) allocate( phi(grid_size_sg,basis_size) )
            call get_reals_i(wc,phi(:,j))
          elseif (token == "TPHI") then
            call get_integer_i(wc,j)
            if (.not.allocated( tphi )) allocate( tphi(grid_size_sg,basis_size) )
            call get_reals_i(wc,tphi(:,j))
          elseif (token == "TPHI_LCAO") then
            call get_integer_i(wc,j)
            if (.not.associated( r2d )) allocate( r2d(lcao_size,basis_size) )
            call get_reals_i(wc,r2d(:,j),l=.false.)
            deallocate( r2d )
          elseif (token == "OVERLAP_SIZE") then
            call get_integer_i(wc,overlap_size) 
          elseif (token == "OVERLAP_MATRIX") then
            allocate( input_oij(overlap_size) )
            call get_reals_i(wc,input_oij)
          elseif (token == "KINETIC_ENERGY_MATRIX") then
            allocate( pd%o%kinetic(overlap_size) )
            call get_reals_i(wc,pd%o%kinetic)
          elseif (token == "V_ION_MATRIX") then
            allocate( pd%o%v_ion(overlap_size) )
            call get_reals_i(wc,pd%o%v_ion)
          elseif (token == "DENVHAT_SIZE") then
            call get_integer_i(wc,denvhat_size) 
          elseif (token == "DENSITY") then
            allocate( density(denvhat_size) )
            call get_den_vhat_i(wc,density,pd%o%lut_denvhat) ; if (error()) goto 100
          elseif (token == "V_HAT") then
            allocate( v_hat(denvhat_size) )
            call get_den_vhat_i(wc,v_hat,pd%o%lut_denvhat) ; if (error()) goto 100
          elseif (token == "HARTREE_SIZE") then
            call get_integer_i(wc,hartree_size)
          elseif (token == "V_HARTREE") then
            allocate( v_hartree(hartree_size), lut_v_hartree(hartree_size) )
            call get_vhartree_i(wc,v_hartree,lut_v_hartree)
          elseif (token == "HAT_SELF-ENERGY") then
            call get_integer_i(wc,i)
            allocate( pd%o%hat_selfenergy(0:i) )
            do j = 0,i
              call get_integer_i(wc,l)
              call get_real_i(wc,pd%o%hat_selfenergy(l))
            end do
          elseif (token == "CORETAILSELFENERGY") then
            call get_real_i(wc,pd%o%coretail_selfenergy) 
          elseif (token == "CORETAILHATENERGY") then
            call get_real_i(wc,pd%o%coretail_hatenergy) 
          elseif (token == "ENERGY") then
            call get_real_i(wc,atom_energy)
          elseif (token == "THETA_POINTS") then
            call get_integer_i(wc,pd%o%theta_points) 
          elseif (token == "PHI_POINTS") then
            call get_integer_i(wc,pd%o%phi_points) 
          elseif (token == "END") then
            continue
          else 
            call warn("WARNING: unrecognized token in PAW file: ")
          end if
          call get_next_word_i(wc,token)
          if (wc%error == w_eof) exit
        end do

        call close_word_i(wc)

        pd%o%mass = 0.0_double
        pd%o%valence_electrons = total_electrons - core_electrons

        nlm_size = 0
        do i = 1,size(pd%o%l_value)
          nlm_size = nlm_size + 2*pd%o%l_value(i) + 1
        end do
    
        allocate( pd%o%nl_base(basis_size), pd%o%lut_nlm(3,nlm_size) )
        k = 0
        do n = 1,size(pd%o%l_value)
          pd%o%nl_base(n) = k + 1
          l = pd%o%l_value(n)
          do m = -l,+l
            k = k + 1
            pd%o%lut_nlm(:,k) = (/n,l,m/)
          end do
        end do
    
        allocate( pd%o%oij(nlm_size,nlm_size) )
        pd%o%oij = 0.0_double
        k = 0
        do i = 1,size(pd%o%l_value)
          do j = 1,i
            if (pd%o%l_value(i) == pd%o%l_value(j)) then
              k = k + 1
              l = pd%o%l_value(i)           
              base_i = pd%o%nl_base(i) + l
              base_j = pd%o%nl_base(j) + l
              do m = -l,+l
                pd%o%oij(base_i+m,base_j+m) = input_oij(k)
                pd%o%oij(base_j+m,base_i+m) = input_oij(k)
              end do
            end if
          end do
        end do

        call orbital_matrix_i(pd%o,density,v_hat) 

        call coulomb_matrix_i(pd%o,v_hartree,lut_v_hartree)

        allocate( pd%o%r_sg(grid_size_sg), drdu_sg(grid_size_sg) )
        select case (grid_type)
        case (LINEAR)
          do i = 1,size(pd%o%r_sg)
            pd%o%r_sg(i) = grid_step*real(i-1,double)
            drdu_sg(i) = 1.0_double
          end do
        case (LOG)
          do i = 1,size(pd%o%r_sg)
            pd%o%r_sg(i) = log_grid_pf*(exp(grid_step*real(i-1,double)) - 1.0_double)
            drdu_sg(i) = log_grid_pf*exp(grid_step*real(i-1,double))
          end do
        end select
        allocate( pd%o%r2_sg(grid_size_sg), r2i_sg(grid_size_sg) )
        pd%o%r2_sg(1) = pd%o%r_sg(1)**2
        r2i_sg(1) = 0.0_double
        do i = 2,size(r2i_sg)
          pd%o%r2_sg(i) = pd%o%r_sg(i)**2
          r2i_sg(i) = 1.0_double/pd%o%r2_sg(i)
        end do

        n = basis_size**2
        allocate( pd%o%phir2ij(grid_size_sg,basis_size,basis_size), pd%o%tphir2ij(grid_size_sg,basis_size,basis_size) )
        allocate( pd%o%dphir2ij_dr(grid_size_sg,basis_size,basis_size), pd%o%dtphir2ij_dr(grid_size_sg,basis_size,basis_size) )

        do i = 1,basis_size
          do j = 1,basis_size
            pd%o%phir2ij(:,i,j) = phi(:,i)*phi(:,j)*r2i_sg
            pd%o%phir2ij(1,i,j) = 3.0_double*(pd%o%phir2ij(2,i,j) - pd%o%phir2ij(3,i,j)) + pd%o%phir2ij(4,i,j)
            pd%o%tphir2ij(:,i,j) = tphi(:,i)*tphi(:,j)*r2i_sg
            pd%o%tphir2ij(1,i,j) = 3.0_double*(pd%o%tphir2ij(2,i,j) - pd%o%tphir2ij(3,i,j)) + pd%o%tphir2ij(4,i,j)
            call nderiv_i(grid_step,pd%o%phir2ij(:,i,j),pd%o%dphir2ij_dr(:,i,j))
            pd%o%dphir2ij_dr(:,i,j) = pd%o%dphir2ij_dr(:,i,j)/drdu_sg
            call nderiv_i(grid_step,pd%o%tphir2ij(:,i,j),pd%o%dtphir2ij_dr(:,i,j))
            pd%o%dtphir2ij_dr(:,i,j) = pd%o%dtphir2ij_dr(:,i,j)/drdu_sg
          end do
        end do
        pd%o%theta_points = 12
        pd%o%phi_points = 12

        pd%o%core_density = r2i_sg*pd%o%core_density/four_pi
        pd%o%core_density(1)= 3.0_double*(pd%o%core_density(2) - pd%o%core_density(3)) + pd%o%core_density(4)

        allocate( pd%o%core_grad(grid_size_sg) )
        call nderiv_i(grid_step,pd%o%core_density,pd%o%core_grad)
        pd%o%core_grad = pd%o%core_grad/drdu_sg

        allocate( pd%o%wt_sg(grid_size_sg) )
        if (odd(grid_size_sg)) then
          pd%o%wt_sg(1)                  = (1.0_double/3.0_double)*drdu_sg(1)*grid_step
          pd%o%wt_sg(2:grid_size_sg-1:2) = (4.0_double/3.0_double)*drdu_sg(2:grid_size_sg-1:2)*grid_step
          pd%o%wt_sg(3:grid_size_sg-1:2) = (2.0_double/3.0_double)*drdu_sg(3:grid_size_sg-1:2)*grid_step
          pd%o%wt_sg(grid_size_sg)       = (1.0_double/3.0_double)*drdu_sg(grid_size_sg)*grid_step
        else
          pd%o%wt_sg(1)                  = (1.0_double/3.0_double)*drdu_sg(1)*grid_step
          pd%o%wt_sg(2:grid_size_sg-2:2) = (4.0_double/3.0_double)*drdu_sg(2:grid_size_sg-2:2)*grid_step
          pd%o%wt_sg(3:grid_size_sg-2:2) = (2.0_double/3.0_double)*drdu_sg(3:grid_size_sg-2:2)*grid_step
          pd%o%wt_sg(grid_size_sg-1)     = (5.0_double/6.0_double)*drdu_sg(grid_size_sg-1)*grid_step
          pd%o%wt_sg(grid_size_sg)       = (1.0_double/2.0_double)*drdu_sg(grid_size_sg)*grid_step
        end if

        pd%o%vlocal = four_pi*pd%o%vlocal*pd%o%r2_sg*pd%o%wt_sg

        n = 2*maxval(pd%o%l_value) + 1
        allocate( pd%o%hatden(grid_size_sg,n) )
        allocate( r1d(grid_size_sg) )
        do i = 1,n
          l = i - 1
          r1d = pd%o%r_sg**(2*l+2)*input_hatden*pd%o%wt_sg
          pf = sqrt(four_pi)/sum(r1d)
          pd%o%hatden(:,i) = pf*input_hatden*pd%o%r_sg**(l+2)*pd%o%wt_sg
        end do
        deallocate( r1d )

        allocate( pd%o%denmat(grid_size_sg,denvhat_size) )
        do i = 1,denvhat_size
          call denvhat_decode_i(pd%o,i,nili,njlj,l)
          pf = sqrt(four_pi)
          pd%o%denmat(:,i) = pf*tphi(:,nili)*tphi(:,njlj)*pd%o%wt_sg
        end do

        do i = 1,basis_size
          pd%o%tp(:,i) = pd%o%tp(:,i)*pd%o%r_sg*pd%o%wt_sg
        end do
        pr_tag = "projector_radius_"//pd%o%name
        call arg(trim(pr_tag),pd%o%r_opt,found)
        if (.not.found) pd%o%r_opt = 4.0_double


        if (error(grid_size_cg < grid_size_sg,"ERROR: grid_size_cg < grid_size_sg")) goto 100
        allocate( pd%o%r_cg(grid_size_cg), drdu_cg(grid_size_cg) )
        select case (grid_type)
        case (LINEAR)
          do i = 1,size(pd%o%r_cg)
            pd%o%r_cg(i) = grid_step*real(i-1,double)
            drdu_cg(i) = 1.0_double
          end do
        case (LOG)
          do i = 1,size(pd%o%r_cg)
            pd%o%r_cg(i) = log_grid_pf*(exp(grid_step*real(i-1,double)) - 1.0_double)
            drdu_cg(i) = log_grid_pf*exp(grid_step*real(i-1,double))
          end do
        end select
        allocate( r2_cg(grid_size_cg), r2i_cg(grid_size_cg) )
        r2_cg(1) = pd%o%r_cg(1)**2
        r2i_cg(1) = 0.0_double
        do i = 2,size(r2i_cg)
          r2_cg(i) = pd%o%r_cg(i)**2
          r2i_cg(i) = 1.0_double/r2_cg(i)
        end do

        pd%o%coretail_density_cg = r2i_cg*pd%o%coretail_density_cg/four_pi
        pd%o%coretail_density_cg(1) = 3.0_double*(pd%o%coretail_density_cg(2) - pd%o%coretail_density_cg(3)) &
                                       + pd%o%coretail_density_cg(4)
        allocate( pd%o%coretail_density_sg(grid_size_sg) )
        pd%o%coretail_density_sg = pd%o%coretail_density_cg(1:grid_size_sg)
        allocate( r1d(grid_size_sg) )
        r1d = four_pi*pd%o%r2_sg*(pd%o%core_density - pd%o%coretail_density_sg)*drdu_sg
        pd%o%qeffion = -total_electrons + simpson_integral(r1d,grid_step)
        deallocate( r1d )

        allocate( pd%o%coretail_grad_sg(grid_size_sg) )
        call nderiv_i(grid_step,pd%o%coretail_density_sg,pd%o%coretail_grad_sg)
        pd%o%coretail_grad_sg = pd%o%coretail_grad_sg/drdu_sg

        allocate( wt_cg(grid_size_cg) )
        if (odd(grid_size_cg)) then
          wt_cg(1)                  = (1.0_double/3.0_double)*drdu_cg(1)*grid_step
          wt_cg(2:grid_size_cg-1:2) = (4.0_double/3.0_double)*drdu_cg(2:grid_size_cg-1:2)*grid_step
          wt_cg(3:grid_size_cg-1:2) = (2.0_double/3.0_double)*drdu_cg(3:grid_size_cg-1:2)*grid_step
          wt_cg(grid_size_cg)       = (1.0_double/3.0_double)*drdu_cg(grid_size_cg)*grid_step
        else
          wt_cg(1)                  = (1.0_double/3.0_double)*drdu_cg(1)*grid_step
          wt_cg(2:grid_size_cg-2:2) = (4.0_double/3.0_double)*drdu_cg(2:grid_size_cg-2:2)*grid_step
          wt_cg(3:grid_size_cg-2:2) = (2.0_double/3.0_double)*drdu_cg(3:grid_size_cg-2:2)*grid_step
          wt_cg(grid_size_cg-1)     = (5.0_double/6.0_double)*drdu_cg(grid_size_cg-1)*grid_step
          wt_cg(grid_size_cg)       = (1.0_double/2.0_double)*drdu_cg(grid_size_cg)*grid_step
        end if

        pd%o%coretail_density_cg = four_pi*pd%o%coretail_density_cg*r2_cg*wt_cg

100     if (associated( lut_v_hartree )) deallocate( lut_v_hartree )
        if (allocated( r2i_sg )) deallocate( r2i_sg )
        if (allocated( drdu_sg )) deallocate( drdu_sg )
        if (allocated( r2_cg )) deallocate( r2_cg )
        if (allocated( r2i_cg )) deallocate( r2i_cg )
        if (allocated( wt_cg )) deallocate( wt_cg )
        if (allocated( drdu_cg )) deallocate( drdu_cg )
        if (allocated( density )) deallocate( density )
        if (allocated( v_hat )) deallocate( v_hat )
        if (associated( input_hatden )) deallocate( input_hatden )
        if (associated( input_oij )) deallocate( input_oij )
        if (associated( r1d )) deallocate( r1d )
        if (associated( v_hartree )) deallocate( v_hartree )
        if (allocated( phi )) deallocate( phi )
        if (allocated( tphi )) deallocate( tphi )
        if (associated( r2d )) deallocate( r2d )

        call glean(thy(f))

        if (error("Exit paw_data_mod::constructor_pd")) continue

      end function 

      subroutine my_pd(pd)
!doc$ subroutine my(pd)
        type(paw_data_obj) :: pd

!cod$
        pd%ref = pd%ref + 1
        pd%o%ref = pd%o%ref + 1
      end subroutine

      subroutine my_new_pd(pdi,pd)
!doc$ subroutine my(pdi,pd)
        type(paw_data_obj) :: pdi
        type(paw_data_obj) :: pd

!cod$ 
        pd%ref = 1
        pd%o => pdi%o
        pd%o%ref = pd%o%ref + 1
      end subroutine

      function thy_pd(pd) result(pdo)
!doc$ function thy(pd) result(pdo)
        type(paw_data_obj) :: pd, pdo

!cod$
        pd%ref = pd%ref-1
        pd%o%ref = pd%o%ref-1
        pdo%ref = pd%ref
        pdo%o => pd%o
      end function

      subroutine glean_pd(pd)
!doc$ subroutine glean(pd)
        type(paw_data_obj) :: pd

!cod$
        if (pd%o%ref < 1) then 
          if (associated( pd%o%r_sg )) deallocate( pd%o%r_sg )
          if (associated( pd%o%r2_sg )) deallocate( pd%o%r2_sg )
          if (associated( pd%o%wt_sg )) deallocate( pd%o%wt_sg )
          if (associated( pd%o%r_cg )) deallocate( pd%o%r_cg )
          if (associated( pd%o%l_value )) deallocate( pd%o%l_value )
          if (associated( pd%o%occupation )) deallocate( pd%o%occupation )
          if (associated( pd%o%nl_base )) deallocate( pd%o%nl_base )
          if (associated( pd%o%lut_nlm )) deallocate( pd%o%lut_nlm )
          if (associated( pd%o%phir2ij )) deallocate( pd%o%phir2ij )
          if (associated( pd%o%dphir2ij_dr )) deallocate( pd%o%dphir2ij_dr )
          if (associated( pd%o%tphir2ij)) deallocate( pd%o%tphir2ij )
          if (associated( pd%o%dtphir2ij_dr )) deallocate( pd%o%dtphir2ij_dr )
          if (associated( pd%o%tp )) deallocate( pd%o%tp )
          if (associated( pd%o%w_max )) deallocate( pd%o%w_max )
          if (associated( pd%o%g )) deallocate( pd%o%g )
          if (associated( pd%o%tpfo )) deallocate( pd%o%tpfo )
          if (associated( pd%o%kinetic )) deallocate( pd%o%kinetic )
          if (associated( pd%o%v_ion )) deallocate( pd%o%v_ion )
          if (associated( pd%o%vlocal )) deallocate( pd%o%vlocal )
          if (associated( pd%o%oij )) deallocate( pd%o%oij )
          if (associated( pd%o%hatden )) deallocate( pd%o%hatden )
          if (associated( pd%o%hat_selfenergy )) deallocate( pd%o%hat_selfenergy )
          if (associated( pd%o%denmat )) deallocate( pd%o%denmat )
          if (associated( pd%o%lut_denvhat )) deallocate( pd%o%lut_denvhat )
          if (associated( pd%o%cijkl )) deallocate( pd%o%cijkl )
          if (associated( pd%o%lut_cijkl )) deallocate( pd%o%lut_cijkl )
          if (associated( pd%o%aqlm )) deallocate( pd%o%aqlm )
          if (associated( pd%o%avlm )) deallocate( pd%o%avlm )
          if (associated( pd%o%lut_orbital )) deallocate( pd%o%lut_orbital )
          if (associated( pd%o%core_density )) deallocate( pd%o%core_density )
          if (associated( pd%o%core_grad )) deallocate( pd%o%core_grad )
          if (associated( pd%o%coretail_density_cg )) deallocate( pd%o%coretail_density_cg )
          if (associated( pd%o%coretail_density_sg )) deallocate( pd%o%coretail_density_sg )
          if (associated( pd%o%coretail_grad_sg )) deallocate( pd%o%coretail_grad_sg )
          deallocate( pd%o )
        end if
        
      end subroutine 

      subroutine bequeath_pd(pd)
!doc$ subroutine bequeath(pd)
        type(paw_data_obj) :: pd

!cod$
        continue
      end subroutine

      subroutine assign_pd(pd,pd2)
!doc$ subroutine assignment(=)(pd,pd2)
        type(paw_data_obj), intent(inout) :: pd
        type(paw_data_obj), intent(in) :: pd2

!cod$
        type(paw_data_obj) :: pdt
        call my(pd2)
        pdt%o => pd%o
        pd%o%ref = pd%o%ref - pd%ref
        pd%o => pd2%o
        pd%o%ref = pd%o%ref + pd%ref
        call glean(pdt)
        call glean(thy(pd2))
      end subroutine

      function pd_ref(pd) result(r)
!doc$ function x_ref(pd) result(r)
        type(paw_data_obj) :: pd
        integer, dimension(2) :: r
!       effects: Returns pd%ref and pd%o%ref.

!cod$
        r(1) = pd%ref
        r(2) = pd%o%ref
        call glean(pd)
      end function

      function pd_name(pd) result(c)
!doc$ function x_name(pd) result(c)
        type(paw_data_obj) :: pd
        character(10) :: c

!cod$
        call my(pd)
        c = pd%o%name
        call glean(thy(pd))
      end function 

      function pd_valence_electrons(pd) result(ve)
!doc$ function x_valence_electrons(pd) result(ve)
        type(paw_data_obj) :: pd
        real(double) :: ve

!cod$
        call my(pd)
        ve = pd%o%valence_electrons
        call glean(thy(pd)) 
      end function 

      function pd_matching_radius(pd) result(mr)
!doc$ function x_matching_radius(pd) result(mr)
        type(paw_data_obj) :: pd
        real(double) :: mr

!cod$
        call my(pd)
        mr = pd%o%matching_radius
        call glean(thy(pd)) 
      end function 

      function pd_xc_type(pd) result(t)
!doc$ function x_xc_type(pd) result(t)
        type(paw_data_obj) :: pd
        character(tag_sz) :: t

!cod$
        call my(pd)
        t = pd%o%xc_type
        call glean(thy(pd))
      end function 

      function pd_mass(pd) result(m)
!doc$ function x_mass(pd) result(m)
        type(paw_data_obj) :: pd
        real(double) :: m

!cod$
        call my(pd)
        m = pd%o%mass
        call glean(thy(pd))
      end function 

      function pd_grid_size(pd) result(n)
!doc$ function x_grid_size(pd) result(n)
        type(paw_data_obj) :: pd
        integer :: n

!cod$  
        call my(pd)
        n = size(pd%o%r_sg)
        call glean(thy(pd))
      end function 

      subroutine type_r2(pd,p)
!doc$ subroutine type_r2(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%r2_sg.

!cod$
        call my(pd)
        p => pd%o%r2_sg
        call glean(thy(pd))
      end subroutine

      subroutine type_wt(pd,p)
!doc$ subroutine type_wt(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%wt_sg.

!cod$
        call my(pd)
        p => pd%o%wt_sg
        call glean(thy(pd))
      end subroutine

      function pd_l_value(pd,i) result(l)
!doc$ function x_l_value(pd,i) result(l)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        integer :: l

!cod$
        call my(pd)
        if (error(i < lbound(pd%o%l_value,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%l_value,1),"ERROR: i is above the upper bound")) goto 100
        l = pd%o%l_value(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_l_value")) continue
      end function

      function pd_occupation(pd,i) result(o)
!doc$ function x_occupation(pd,i) result(o)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: o

!cod$
        call my(pd)
        if (error(i < lbound(pd%o%occupation,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%occupation,1),"ERROR: i is above the upper bound")) goto 100
        o = pd%o%occupation(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_occupation")) continue
      end function

      function pd_nl_base(pd,i) result(n)
!doc$ function x_nl_base(pd,i) result(n)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        integer :: n

!cod$
        call my(pd)
        if (error(i < lbound(pd%o%nl_base,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%nl_base,1),"ERROR: i is above the upper bound")) goto 100
        n = pd%o%nl_base(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_nl_base")) continue
      end function

      function pd_nlm_size(pd) result(n)
!doc$ function x_nlm_size(pd) result(n)
        type(paw_data_obj) :: pd
        integer :: n

!cod$
        call my(pd)        
        n = size(pd%o%lut_nlm,2)
        call glean(thy(pd))
      end function

      function pd_nlm(pd,i) result(nlm)
!doc$ function x_nlm(pd,i) result(nlm)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        integer, dimension(3) :: nlm

!cod$
        call my(pd)
        if (error(i < lbound(pd%o%lut_nlm,2),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%lut_nlm,2),"ERROR: i is above the upper bound")) goto 100
        nlm = pd%o%lut_nlm(:,i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_lut_nlm")) continue
      end function 

      function pd_basis_size(pd) result(n)
!doc$ function x_basis_size(pd) result(n)
        type(paw_data_obj) :: pd 
        integer :: n

!cod$
        call my(pd)
        n = size(pd%o%l_value)
        call glean(thy(pd))
      end function 

      subroutine type_phir2ij(pd,p)
!doc$ subroutine type_phir2ij(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:,:,:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%phir2ij.

!cod$
        call my(pd)
        p => pd%o%phir2ij
        call glean(thy(pd))
      end subroutine

      subroutine type_dphir2ij_dr(pd,p)
!doc$ subroutine type_dphir2ij_dr(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:,:,:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%dphir2ij_dr.

!cod$
        call my(pd)
        p => pd%o%dphir2ij_dr
        call glean(thy(pd))
      end subroutine

      subroutine type_tphir2ij(pd,p)
!doc$ subroutine type_tphir2ij(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:,:,:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%tphir2ij.

!cod$
        call my(pd)
        p => pd%o%tphir2ij
        call glean(thy(pd))
      end subroutine

      subroutine type_dtphir2ij_dr(pd,p)
!doc$ subroutine type_dtphir2ij_dr(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:,:,:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%dtphir2ij_dr.

!cod$
        call my(pd)
        p => pd%o%dtphir2ij_dr
        call glean(thy(pd))
      end subroutine

      function pd_theta_points(pd) result(p)
!doc$ function x_theta_points(pd) result(p)
        type(paw_data_obj) :: pd
        integer :: p
!       effects: Returns the number of theta points.

!cod$
        call my(pd)
        p = pd%o%theta_points
        call glean(thy(pd))
      end function

      function pd_phi_points(pd) result(p)
!doc$ function x_phi_points(pd) result(p)
        type(paw_data_obj) :: pd
        integer :: p
!       effects: Returns the number of phi points.

!cod$
        call my(pd)
        p = pd%o%phi_points
        call glean(thy(pd))
      end function

      function pd_projector_f_value(pd,g,ip) result(pfv)
!doc$ function projector_f_value(pd,g,ip) result(pfv)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: g
        integer, intent(in) :: ip
        real(double) :: pfv
!       effect: Returns the g Fourier coefficient of projector ip.
!       error: ip out of range. Passes errors.

!cod$
        logical :: switch
        integer :: l
        real(double), dimension(:), allocatable :: f
        if (error((ip < 1) .or. (ip > size(pd%o%tp,2)),"ERROR: ip is out of range")) goto 100
        allocate( f(size(pd%o%r_sg)) )
        l = pd%o%l_value(ip)
        switch = .true.
        f = four_pi*pd%o%tp(:,ip)*pd%o%r_sg**l
        pfv = radial_f_value_i(f,pd%o%r_sg,g,l,switch)
100     if (allocated( f )) deallocate( f )
        if (error("Exit paw_data_mod::pd_projector_f_value")) continue
      end function

      function pd_projector_f_values(pd,g,ip) result(pfv)
!doc$ function projector_f_values(pd,g,ip) result(pfv)
        type(paw_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        integer, intent(in) :: ip
        real(double), dimension(size(g)) :: pfv
!       effect: Returns the g Fourier coefficients of projector ip.
!       error: ip out of range. Passes errors.

!cod$
        logical :: switch
        integer :: ig, l
        real(double), dimension(:), allocatable :: f
        if (error((ip < 1) .or. (ip > size(pd%o%tp,2)),"ERROR: ip is out of range")) goto 100
        allocate( f(size(pd%o%r_sg)) )
        l = pd%o%l_value(ip)
        switch = .true.
        f = four_pi*pd%o%tp(:,ip)*pd%o%r_sg**l
        do ig = 1,size(g)
          pfv(ig) = radial_f_value_i(f,pd%o%r_sg,g(ig),l,switch)
        end do
100     if (allocated( f )) deallocate( f )
        if (error("Exit paw_data_mod::pd_projector_f_value")) continue
      end function

      function pd_radius(pd) result(r)
!doc$ function x_radius(pd) result(r)
        type(paw_data_obj) :: pd
        real(double) :: r
!       effects: Returns pd%o%r_opt in Bohr.

!cod$
        call my(pd)
        r = pd%o%r_opt
        call glean(thy(pd))
      end function

      function pd_inner_cutoff(pd) result(c)
!doc$ function x_inner_cutoff(pd) result(c)
        type(paw_data_obj) :: pd
        real(double) :: c
!       effects: Returns pd%o%g_icut in sqrt(Ryd).

!cod$
        call my(pd)
        c = pd%o%g_icut
        call glean(thy(pd))
      end function

      function pd_outer_cutoff(pd) result(c)
!doc$ function x_outer_cutoff(pd) result(c)
        type(paw_data_obj) :: pd
        real(double) :: c
!       effects: Returns pd%o%g_ocut in sqrt(Ryd).

!cod$
        call my(pd)
        c = pd%o%g_ocut
        call glean(thy(pd))
      end function

      function pd_error(pd,i) result(e)
!doc$ function x_error(pd,i) result(e)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: e
!       effects: Returns the i'th tpfo error in pd.
!       errors: If pd%o%w_max is not associated. If i is out of bounds.

!cod$
        call my(pd)
        if (error(.not.associated( pd%o%w_max ),"ERROR: w_max is not associated")) goto 100
        if (error((i < 1) .or. (i > size(pd%o%w_max)),"ERROR: i is out of bounds")) goto 100
        e = pd%o%w_max(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_error")) continue
      end function

      function projector_r_value_pd(pd,r,ipb,g_icut,g_ocut) result(prv)
!doc$ function projector_r_value(pd,r,ipb,g_icut,g_ocut) result(prv)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: r
        integer, intent(in) :: ipb
        real(double), intent(in) :: g_icut, g_ocut
        real(double) :: prv
!       effects: Returns the value of an optimized projector at r. Initializes pd%o%tpfo if necessary.
!       errors: r > pd%o%r_opt. ipb out of range. Passes errors.

!cod$
        logical :: g_icut_change, g_ocut_change
        integer :: l

        call my(pd)

        if (error(r > pd%o%r_opt,"ERROR: r > r_opt")) goto 100
        if (error((ipb < 1) .or. (ipb > size(pd%o%l_value)),"ERROR: ipb is out of range")) goto 100

        if (.not.associated( pd%o%tpfo )) then
          pd%o%g_icut = g_icut
          pd%o%g_ocut = g_ocut
          call optimize_nlp_i(pd%o) ; if (error()) goto 100
        else
          g_icut_change = ( g_icut /= pd%o%g_icut )
          g_ocut_change = ( g_ocut /= pd%o%g_ocut )
          if (g_icut_change .or. g_ocut_change) then
            if (g_icut_change) pd%o%g_icut = g_icut
            if (g_ocut_change) pd%o%g_ocut = g_ocut
            if (associated( pd%o%tpfo )) deallocate( pd%o%tpfo )
            if (associated( pd%o%g )) deallocate( pd%o%g )
            if (associated( pd%o%w_max )) deallocate( pd%o%w_max )
            call warn("WARNING: real-space projectors are being re-optimized")
            call optimize_nlp_i(pd%o) ; if (error()) goto 100
          end if
        end if

        l = pd%o%l_value(ipb)
        select case(l)
        case(0)
          prv = sum(pd%o%tpfo(:,ipb)*spherical_bessel(pd%o%g*r,l,.true.)) ; if (error()) goto 100
        case default
          prv = sum(pd%o%tpfo(:,ipb)*pd%o%g**l*spherical_bessel(pd%o%g*r,l,.true.)) ; if (error()) goto 100
        end select

100     call glean(thy(pd))

        if (error("Exit paw_data_mod::projector_r_value_pd")) continue

      end function

      function projector_r_gradients_pd(pd,r,ipb) result(prg)
!doc$ function projector_r_gradients(pd,r,ipb) result(prg)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: r
        integer, intent(in) :: ipb
        real(double), dimension(2) :: prg
!       effects: Returns the gradients of the ipb optimized real-space projector at r.
!       errors: r > pd%o%r_opt. ipb out of range. pd%o%tpfo not associated. Passes errors.

!cod$
        integer :: l

        call my(pd)

        if (error(r > pd%o%r_opt,"ERROR: r > r_opt")) goto 100
        if (error((ipb < 1) .or. (ipb > size(pd%o%l_value)),"ERROR: ipb is out of range")) goto 100
        if (error(.not.associated( pd%o%tpfo ),"ERROR: real-space projectors are not yet optimized")) goto 100

        l = pd%o%l_value(ipb)
        select case(l)
        case(0)
          prg(1) = 0.0_double
        case default
          prg(1) = sum(pd%o%tpfo(:,ipb)*pd%o%g**(l)*spherical_bessel(pd%o%g*r,l,.true.)) ; if (error()) goto 100
        end select
        prg(2) = sum(pd%o%tpfo(:,ipb)*pd%o%g**(l+2)*spherical_bessel(pd%o%g*r,l+1,.true.)) ; if (error()) goto 100

100     call glean(thy(pd))

        if (error("Exit paw_data_mod::projector_r_gradients_pd")) continue

      end function

      subroutine diary_rs_projectors_pd(pd)
!doc$ subroutine diary_rs_projectors(pd)
        type(paw_data_obj) :: pd
!       requires: Real-space projectors be in use.
!       effects: Writes real-space projector information to the diary.

!cod$
        integer :: ib
        call my(pd)
        if (i_access( diaryfile() )) then
          select case (len_trim(pd%o%name))
          case (1)
            write(x_unit(diaryfile()),'(/,t6,a1,":")') trim(pd%o%name)
          case (2)
            write(x_unit(diaryfile()),'(/,t6,a2,":")') trim(pd%o%name)
          case (3)
            write(x_unit(diaryfile()),'(/,t6,a3,":")') trim(pd%o%name)
          case (4)
            write(x_unit(diaryfile()),'(/,t6,a4,":")') trim(pd%o%name)
          case (5)
            write(x_unit(diaryfile()),'(/,t6,a5,":")') trim(pd%o%name)
          case default
            write(x_unit(diaryfile()),'(/,t6,a6,":")') trim(pd%o%name)
          end select
          write(x_unit(diaryfile()),'(t8,"optimization radius = ",f4.2," Bohr")') pd%o%r_opt
          write(x_unit(diaryfile()),'(t8,"inner cutoff energy = ",f7.2," Ryd")') pd%o%g_icut**2
          write(x_unit(diaryfile()),'(t8,"outer cutoff energy = ",f7.2," Ryd")') pd%o%g_ocut**2
          do ib = 1,size(pd%o%l_value)
            write(x_unit(diaryfile()),'(t8,"l = ",i1,": maximum Fourier error = ",es8.2)') pd%o%l_value(ib), pd%o%w_max(ib)
          end do
        end if
        call glean(thy(pd))

      end subroutine

      function pd_kinetic(pd,i) result(k)
!doc$ function x_kinetic(pd,i) result(k)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: k

!cod$
        call my(pd)
        if (error(i < lbound(pd%o%kinetic,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%kinetic,1),"ERROR: i is above the upper bound")) goto 100
        k = pd%o%kinetic(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_kinetic")) continue
      end function

      function pd_v_ion(pd,i) result(v)
!doc$ function x_v_ion(pd,i) result(v)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: v

!cod$
        call my(pd)  
        if (error(i < lbound(pd%o%v_ion,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%v_ion,1),"ERROR: i is above the upper bound")) goto 100
        v = pd%o%v_ion(i)
100     call glean(thy(pd)) 
        if (error("Exit paw_data_mod::pd_v_ion")) continue
      end function 

      function pd_oij(pd) result(r)
!doc$ function x_oij(pd) result(r)
        type(paw_data_obj) :: pd
        real(double), dimension(size(pd%o%oij,1),size(pd%o%oij,2)) :: r

!cod$
        call my(pd)
        r = pd%o%oij
        call glean(thy(pd))
      end function 

      function local_f_value_pd(pd,g) result(f)
!doc$ function local_f_value(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: g
        real(double) :: f
!       effects: Returns the Fourier coefficient of the local pseudopotential.
!       errors: Passes errors.

!cod$
        integer :: l
        call my(pd)
        l = 0
        f = radial_f_value_i(pd%o%vlocal,pd%o%r_sg,g,l)
        call glean(thy(pd))
        if (error("Exit paw_data_mod::local_f_value_pd")) continue
      end function

      function local_f_values_pd(pd,g) result(f)
!doc$ function local_f_values(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        real(double), dimension(size(g)) :: f
!       effects: Returns the Fourier coefficients of the local pseudopotential.
!       errors: Passes errors.

!cod$
        integer :: ig, l
        call my(pd)
        l = 0
        do ig = 1,size(f)
          f(ig) = radial_f_value_i(pd%o%vlocal,pd%o%r_sg,g(ig),l)
        end do
        call glean(thy(pd))
        if (error("Exit paw_data_mod::local_f_values_pd")) continue
      end function

      function pd_hatden_size(pd) result(n)
!doc$ function x_hatden_size(pd) result(n)
        type(paw_data_obj) :: pd
        integer :: n

!cod$
        call my(pd)  
        n = size(pd%o%hatden,2)
        call glean(thy(pd)) 
      end function

      function hatden_f_value_pd(pd,g) result(f)
!doc$ function hatden_f_value(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: g
        real(double), dimension(size(pd%o%hatden,2)) :: f
!       effects: Returns the Fourier coefficient of the hat density.
!       errors: Passes errors.

!cod$
        integer :: i, l
        call my(pd)
        do i = 1,size(f)
          l = i - 1
          f(i) = radial_f_value_i(pd%o%hatden(:,i),pd%o%r_sg,g,l)
        end do
        call glean(thy(pd))
        if (error("Exit paw_data_mod::hatden_f_value_pd")) continue
      end function

      function hatden_f_values_pd(pd,g) result(f)
!doc$ function hatden_f_values(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        real(double), dimension(size(g),size(pd%o%hatden,2)) :: f
!       effects: Returns the Fourier coefficients of the hat density.
!       errors: Passes errors.

!cod$
        integer :: i, ig, l
        call my(pd)
        do i = 1,size(f,2)
          l = i - 1
          do ig = 1,size(f,1)
            f(ig,i) = radial_f_value_i(pd%o%hatden(:,i),pd%o%r_sg,g(ig),l)
          end do
        end do
        call glean(thy(pd))
        if (error("Exit paw_data_mod::hatden_f_values_pd")) continue
      end function

      function pd_hat_self_energy(pd,i) result(r)
!doc$ function x_hat_self_energy(pd,i) result(r)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: r

!cod$
        call my(pd)
        if (error(i < lbound(pd%o%hat_selfenergy,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%hat_selfenergy,1),"ERROR: i is above the upper bound")) goto 100
        r = pd%o%hat_selfenergy(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_hat_self_energy")) continue
      end function

      function pd_denmat_size(pd) result(n)
!doc$ function x_denmat_size(pd) result(n)
        type(paw_data_obj) :: pd
        integer :: n

!cod$
        call my(pd)
        n = size(pd%o%denmat,2)
        call glean(thy(pd))
      end function 

      function denmat_f_value_pd(pd,g) result(f)
!doc$ function denmat_f_value(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: g
        real(double), dimension(size(pd%o%denmat,2)) :: f
!       effects: Returns the Fourier coefficient of the density matrix.
!       errors: Passes errors.

!cod$
        integer :: i, l, nili, njlj
        call my(pd)
        do i = 1,size(f)
          call denvhat_decode_i(pd%o,i,nili,njlj,l)
          f(i) = radial_f_value_i(pd%o%denmat(:,i),pd%o%r_sg,g,l)
        end do
        call glean(thy(pd))
        if (error("Exit paw_data_mod::denmat_f_value_pd")) continue
      end function

      function denmat_f_values_pd(pd,g) result(f)
!doc$ function denmat_f_values(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        real(double), dimension(size(g),size(pd%o%denmat,2)) :: f
!       effects: Returns the Fourier coefficients of the density matrix.
!       errors: Passes errors.

!cod$
        integer :: i, ig, l, nili, njlj
        call my(pd)
        do i = 1,size(f,2)
          call denvhat_decode_i(pd%o,i,nili,njlj,l)
          do ig = 1,size(f,1)
            f(ig,i) = radial_f_value_i(pd%o%denmat(:,i),pd%o%r_sg,g(ig),l)
          end do
        end do
        call glean(thy(pd))
        if (error("Exit paw_data_mod::denmat_f_values_pd")) continue
      end function

      subroutine denvhat_decode(pd,i,nili,njlj,l)
!doc$ subroutine denvhat_decode(pd,i,nili,njlj,l)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        integer, intent(out) :: nili, njlj, l
!       effects: Decodes the denvhat indices.
!       errors: i out of range.

!cod$
        integer, parameter :: dd_p8 = 7936
        integer, parameter :: dd_p3 = 248
        integer :: n
        call my(pd)
        if (error(i < lbound(pd%o%lut_denvhat,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%lut_denvhat,1),"ERROR: i is above the upper bound")) goto 100
        n = pd%o%lut_denvhat(i)
        nili = ishft(iand(dd_p8,n),-8)
        njlj = ishft(iand(dd_p3,n),-3)
        l = iand(7,n)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::denvhat_decode")) continue
      end subroutine

      function gaunt_complex(l,m,l1,m1,l2,m2) result(gc)
!doc$ function gaunt_complex(l,m,l1,m1,l2,m2) result(gc)
        integer, intent(in) :: l, m, l1, l2, m1, m2
        real(double) :: gc

!cod$
        logical :: ok
        integer :: j, j1, j2, j3, j1half, j2half, j3half, j_half, k1, k2, i1, i2, n1, n2
        real(double) :: xx,yy, sum, sign, argument
        gc = 0.0_double
        sum = 0.0_double
        ok = .true.
        if ((m2 - m1 - m) /= 0) ok = .false.
        if (abs(m1) > l1) ok = .false.
        if (abs(m) > l) ok = .false.
        if (abs(m2) > l2) ok = .false.
        j = l1 + l + l2
        if (mod(j,2) /= 0) ok = .false.
        j1 = j - 2*l2
        j2 = j - 2*l
        j3 = j - 2*l1
        if ((j1 < 0) .or. (j2 < 0) .or. (j3 < 0)) ok = .false.
        if (ok) then
          xx = (2*l1 + 1)*(2*l + 1)*(2*l2 + 1)
          j1half = j1/2
          j2half = j2/2
          j3half = j3/2
          j_half = j/2
          gc = (-1)**j1half*sqrt(xx)
          gc = gc*factorial(j2)*factorial(j3)/factorial(j+1)
          gc = gc*factorial(j_half)/(factorial(j1half)*factorial(j2half)*factorial(j3half))
          yy = factorial(l2+m2)*factorial(l2-m2)
          if (m >= 0) then
            yy = yy*permutations(l+m,2*m)
          else
            yy = yy/permutations(l-m,-2*m)
          end if
          if (m1 >= 0) then
            yy = yy/permutations(l1+m1,2*m1)
          else
            yy = yy*permutations(l1-m1,-2*m1)
          end if
          gc = gc*sqrt(yy)
          i1 = l2 - l - m1
          i2 = l2 - l1 + m
          k1 = -min(0,i1,i2)
          n1 = l1 + m1
          n2 = l - m
          k2 = min(j1,n1,n2)
          sign = 1.0_double
          if(k1 > 0) sign = (-1)**k1
          argument = sign*permutations(n1,k1)/factorial(k1)
          argument = argument*permutations(n2,k1)/factorial(i1+k1)
          argument = argument*permutations(j1,k1)/factorial(i2+k1)
          sum = sum + argument
          sign = -sign
          k1 = k1 + 1
          do while(k1 <= k2)
            argument = sign*permutations(n1,k1)/factorial(k1)
            argument = argument*permutations(n2,k1)/factorial(i1+k1)
            argument = argument*permutations(j1,k1)/factorial(i2+k1)
            sum = sum + argument
            sign = -sign
            k1 = k1 + 1
          end do
        end if
        gc = gc*sum
      end function

      function pd_cijkl_size(pd) result(n)
!doc$ function x_cijkl_size(pd) result(n)
        type(paw_data_obj) :: pd
        integer :: n

!cod$
        call my(pd)
        n = size(pd%o%cijkl)
        call glean(thy(pd))
      end function

      function pd_cijkl(pd, i) result(r)
!doc$ function x_cijkl(pd, i) result(r)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: r

!cod$
        call my(pd) 
        if (error(i < lbound(pd%o%cijkl,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%cijkl,1),"ERROR: i is above the upper bound")) goto 100
        r = pd%o%cijkl(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_cijkl")) continue
      end function

      subroutine cijkl_decode(pd,i,nili,njlj,nklk,nlll,mi,mj,mk)
!doc$ subroutine cijkl_decode(pd,i,nili,njlj,nklk,nlll,mi,mj,mk)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        integer, intent(out) :: nili, njlj, nklk, nlll, mi, mj, mk
!       effects: Decodes the cijkl indices.
!       errors: i out of range.

!cod$
        integer, parameter :: cd_p24 = 520093696
        integer, parameter :: cd_p19 = 16252928
        integer, parameter :: cd_p14 = 507904
        integer, parameter :: cd_p9 = 15872
        integer, parameter :: cd_p6 = 448
        integer, parameter :: cd_p3 = 56
        integer :: n
        call my(pd)
        if (error(i < lbound(pd%o%lut_cijkl,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%lut_cijkl,1),"ERROR: i is above the upper bound")) goto 100
        n = pd%o%lut_cijkl(i)
        nili = ishft(iand(cd_p24,n),-24)
        njlj = ishft(iand(cd_p19,n),-19)
        nklk = ishft(iand(cd_p14,n),-14)
        nlll = ishft(iand(cd_p9,n),-9)
        mi = ishft(iand(cd_p6,n),-6) - 3
        mj = ishft(iand(cd_p3,n),-3) - 3
        mk = iand(7,n) - 3
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::cijkl_decode")) continue
      end subroutine

      function pd_qvlm_size(pd) result(n)
!doc$ function x_qvlm_size(pd) result(n)
        type(paw_data_obj) :: pd
        integer :: n

!cod$
        call my(pd) 
        n = size(pd%o%aqlm)
        call glean(thy(pd)) 
      end function

      function pd_aqlm(pd,i) result(r)
!doc$ function x_aqlm(pd,i) result(r)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: r

!cod$  
        call my(pd)
        if (error(i < lbound(pd%o%aqlm,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%aqlm,1),"ERROR: i is above the upper bound")) goto 100
        r = pd%o%aqlm(i)
100     call glean(thy(pd))  
        if (error("Exit paw_data_mod::pd_aqlm")) continue
      end function

      function pd_avlm(pd,i) result(r)
!doc$ function x_avlm(pd,i) result(r)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: r

!cod$
        call my(pd) 
        if (error(i < lbound(pd%o%avlm,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%avlm,1),"ERROR: i is above the upper bound")) goto 100
        r = pd%o%avlm(i)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::pd_avlm")) continue
      end function

      subroutine orbital_decode(pd,i,nili,njlj,mi,mj,l)
!doc$ subroutine orbital_decode(pd,i,nili,njlj,mi,mj,l)
        type(paw_data_obj) :: pd
        integer, intent(in) :: i
        integer, intent(out) :: nili, njlj, mi, mj, l
!       effects: Decodes the orbital indices.
!       errors: i out of range.

!cod$
        integer, parameter :: od_p14 = 507904
        integer, parameter :: od_p9 = 15872
        integer, parameter :: od_p6 = 448
        integer, parameter :: od_p3 = 56
        integer :: n
        call my(pd)
        if (error(i < lbound(pd%o%lut_orbital,1),"ERROR: i is below the lower bound")) goto 100
        if (error(i > ubound(pd%o%lut_orbital,1),"ERROR: i is above the upper bound")) goto 100
        n = pd%o%lut_orbital(i)
        nili = ishft(iand(od_p14,n),-14)
        njlj = ishft(iand(od_p9,n),-9)
        mi = ishft(iand(od_p6,n),-6) - 3
        mj = ishft(iand(od_p3,n),-3) - 3
        l  = iand(7,n)
100     call glean(thy(pd))
        if (error("Exit paw_data_mod::orbital_decode")) continue
      end subroutine

      function pd_qeffion(pd) result(r)
!doc$ function x_qeffion(pd) result(r)
        type(paw_data_obj) :: pd
        real(double) :: r

!cod$
        call my(pd) 
        r = pd%o%qeffion
        call glean(thy(pd))
      end function

      subroutine type_core_density(pd,p)
!doc$ subroutine type_core_density(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%core_density.

!cod$
        call my(pd)
        p => pd%o%core_density
        call glean(thy(pd))
      end subroutine

      subroutine type_core_grad(pd,p)
!doc$ subroutine type_core_grad(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%core_grad.

!cod$
        call my(pd)
        p => pd%o%core_grad
        call glean(thy(pd))
      end subroutine

      subroutine type_coretail_density(pd,p)
!doc$ subroutine type_coretail_density(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%coretail_density_sg.

!cod$
        call my(pd)
        p => pd%o%coretail_density_sg
        call glean(thy(pd))
      end subroutine

      subroutine type_coretail_grad(pd,p)
!doc$ subroutine type_coretail_grad(pd,p)
        type(paw_data_obj) :: pd
        real(double), dimension(:), pointer :: p
!       modifies: p
!       effects: Points p to pd%o%coretail_grad_sg.

!cod$
        call my(pd)
        p => pd%o%coretail_grad_sg
        call glean(thy(pd))
      end subroutine

      function coretail_f_value_pd(pd,g) result(f)
!doc$ function coretail_f_value(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), intent(in) :: g
        real(double) :: f
!       effects: Returns the Fourier coefficient of the coretail density.
!       errors: Passes errors.

!cod$
        integer :: l
        call my(pd)
        l = 0
        f = radial_f_value_i(pd%o%coretail_density_cg,pd%o%r_cg,g,l)
        call glean(thy(pd))
        if (error("Exit paw_data_mod::coretail_f_value_pd")) continue
      end function

      function coretail_f_values_pd(pd,g) result(f)
!doc$ function coretail_f_values(pd,g) result(f)
        type(paw_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        real(double), dimension(size(g)) :: f
!       effects: Returns the Fourier coefficients of the coretail density.
!       errors: Passes errors.

!cod$
        integer :: ig, l
        call my(pd)
        l = 0
        do ig = 1,size(f)
          f(ig) = radial_f_value_i(pd%o%coretail_density_cg,pd%o%r_cg,g(ig),l)
        end do
        call glean(thy(pd))
        if (error("Exit paw_data_mod::coretail_f_values_pd")) continue
      end function

      function pd_coretail_hatenergy(pd) result(r)
!doc$ function x_coretail_hatenergy(pd) result(r)
        type(paw_data_obj) :: pd
        real(double) :: r

!cod$
        call my(pd) 
        r = pd%o%coretail_hatenergy
        call glean(thy(pd))
      end function

      function pd_coretail_selfenergy(pd) result(r)
!doc$ function x_coretail_selfenergy(pd) result(r)
        type(paw_data_obj) :: pd
        real(double) :: r

!cod$
        call my(pd) 
        r = pd%o%coretail_selfenergy
        call glean(thy(pd))
      end function

! local routines

      subroutine nderiv_i(h,y,z)
        real(double), intent(in) :: h
        real(double), dimension(:), intent(in) :: y
        real(double), dimension(:), intent(out) :: z

        integer :: i, n
        real(double) :: hh, a, b, c
        n = size(y)
        hh = 0.08333333333333333_double/h
        b = hh*(-25.0_double*y(1) + 48.0_double*y(2) - 36.0_double*y(3) + 16.0_double*y(4) - 3.0_double*y(5))
        c = hh*(-3.0_double*y(1) - 10.0_double*y(2) + 18.0_double*y(3) - 6.0_double*y(4) + y(5))
        do i = 5,n
          a = b
          b = c
          c = hh*(y(i-4) - y(i) + 8.0_double*(y(i-1) - y(i-3)))
          z(i-4) = a
        end do
        a = hh*(-y(n-4) + 6.0_double*y(n-3) - 18.0_double*y(n-2) + 10.0_double*y(n-1) + 3.0_double*y(n))
        z(n) = hh*(3.0_double*y(n-4) - 16.0_double*y(n-3) + 36.0_double*y(n-2) - 48.0_double*y(n-1) + 25.0_double*y(n))
        z(n-1) = a
        z(n-2) = c
        z(n-3) = b
      end subroutine

      subroutine get_den_vhat_i(wc,a,lut)
        type(word_context) :: wc
        real(double), dimension(:) :: a
        integer, dimension(:), pointer :: lut

        logical :: exact
        integer :: i, j, k
        integer, dimension(3) :: ind
        if (.not.associated( lut )) then
          allocate( lut(size(a)) )
          do i = 1,size(a)
            call get_integers_i(wc,ind)
            j = denvhat_encode_i(ind(1),ind(2),ind(3))
            lut(i) = j
            call get_real_i(wc,a(i))
          end do
        else
          do i = 1,size(a)
            call get_integers_i(wc,ind)
            j = denvhat_encode_i(ind(1),ind(2),ind(3))
            k = linear_search_i(lut,j,exact)
            if (error(.not.exact,"ERROR: match was not found")) goto 100
            call get_real_i(wc,a(k))
          end do
        end if
100     if (error("Exit paw_data_mod::get_den_vhat_i")) continue

      end subroutine 

      subroutine get_vhartree_i(wc,a,lut)
        type(word_context) :: wc
        real(double), dimension(:) :: a
        integer, dimension(:),  pointer :: lut

        integer :: i
        integer, dimension(5) :: ind
        do i = 1,size(a)
          call get_integers_i(wc,ind)
          lut(i) = vhartree_encode_i(ind(1),ind(2),ind(3),ind(4),ind(5))
          call get_real_i(wc,a(i))
        end do
      end subroutine

      subroutine coulomb_matrix_i(pdr,v_hartree,lut_v_hartree)
        type(paw_data_rep) :: pdr
        real(double), dimension(:) :: v_hartree
        integer, dimension(:) :: lut_v_hartree

        logical :: ok
        integer :: i, ih, j, n, nh
        integer :: nili, nili2, njlj, njlj2, nklk, nklk2, nlll, nlll2
        integer :: l, li, li2, lj, lj2, lk, lk2, ll, ll2
        integer :: m, mi, mj, mk, ml
        integer :: loop
        integer, dimension(:), allocatable :: i1d, lut, l_list
        real(double) :: gc1, gc2, msign, term
        real(double), dimension(:), allocatable :: r1d
        nh = size(v_hartree)
        allocate( lut(nh), i1d(nh), r1d(nh) )
        call insertion_sort_i(lut_v_hartree,lut)
        i1d = lut_v_hartree(lut)
        lut_v_hartree = i1d
        r1d = v_hartree(lut)
        v_hartree = r1d
        deallocate( lut, i1d, r1d )
        m = 1
        ih = 1
        do while (ih <= nh)
          call vhartree_decode_i(lut_v_hartree(ih),nili,njlj,nklk,nlll,l)
          nili2 = nili
          njlj2 = njlj
          nklk2 = nklk
          nlll2 = nlll
          i = 0
          do while ((nili2 == nili) .and. (njlj2 == njlj) .and. (nklk2 == nklk) .and. (nlll2 == nlll) .and. ((ih+i) < nh))
            i = i + 1
            call vhartree_decode_i(lut_v_hartree(ih+i),nili2,njlj2,nklk2,nlll2,l)
          end do
          m = max(m,i)
          ih = ih + 1
        end do
        allocate( l_list(m) )
        n = 0
        ih = 1
        do while (ih <= nh)
          call vhartree_decode_i(lut_v_hartree(ih),nili,njlj,nklk,nlll,l)
          nili2 = nili
          njlj2 = njlj
          nklk2 = nklk
          nlll2 = nlll
          i = 0
          do while ((nili2 == nili) .and. (njlj2 == njlj) .and. (nklk2 == nklk) .and. (nlll2 == nlll) .and. ((ih+i) < nh))
            i = i + 1
            l_list(i) = l
            call vhartree_decode_i(lut_v_hartree(ih+i),nili2,njlj2,nklk2,nlll2,l)
          end do
          if (i == 0) then
            i = 1
            l_list(1) = l
          end if
          nili2 = nili ; li = pdr%l_value(nili) ; li2 = li
          njlj2 = njlj ; lj = pdr%l_value(njlj) ; lj2 = lj
          nklk2 = nklk ; lk = pdr%l_value(nklk) ; lk2 = lk
          nlll2 = nlll ; ll = pdr%l_value(nlll) ; ll2 = ll
          do loop = 1,4
            ok = .false.
            if (loop == 1) then
              ok = .true.
            elseif ((loop == 2) .and. (nili /= njlj)) then
              ok = .true.
              nili = njlj2 ; li = lj2
              njlj = nili2 ; lj = li2
              nklk = nklk2 ; lk = lk2
              nlll = nlll2 ; ll = ll2
            elseif ((loop == 3) .and. (nklk /= nlll)) then
              ok = .true.
              nili = nili2 ; li = li2
              njlj = njlj2 ; lj = lj2
              nklk = nlll2 ; lk = ll2
              nlll = nklk2 ; ll = lk2
            elseif ((loop == 4) .and. (nili /= njlj) .and. (nklk /= nlll)) then
              ok = .true.
              nili = njlj2 ; li = lj2
              njlj = nili2 ; lj = li2
              nklk = nlll2 ; lk = ll2
              nlll = nklk2 ; ll = lk2
            end if
            if (ok) then
              do mi = -li,+li
                do mj = -lj,+lj
                  m = mi - mj        
                  msign = (-1.0_double)**m
                  do mk = -lk,+lk
                    do ml = -ll,+ll
                      if (m == (ml-mk)) then    
                        term = 0.0_double
                        do j = 1,i
                          if (abs(m) <= l_list(j)) then
                            gc1 = gaunt_complex(l_list(j),-m,li,mi,lj,mj)
                            gc2 = gaunt_complex(l_list(j),m,lk,mk,ll,ml)
                            term = term + v_hartree(ih+j-1)*gc1*gc2
                          end if
                        end do
                        if (abs(term) > 1.0e-20_double) n = n + 1
                      end if
                    end do
                  end do
                end do
              end do
            end if
          end do
          ih = ih + i
        end do
        if (associated( pdr%cijkl )) deallocate( pdr%cijkl ) ; allocate( pdr%cijkl(n) )
        if (associated( pdr%lut_cijkl )) deallocate( pdr%lut_cijkl ) ; allocate( pdr%lut_cijkl(n) )
        n = 0
        ih = 1
        do while (ih <= nh)
          call vhartree_decode_i(lut_v_hartree(ih),nili,njlj,nklk,nlll,l)
          nili2 = nili
          njlj2 = njlj
          nklk2 = nklk
          nlll2 = nlll
          i = 0
          do while ((nili2 == nili) .and. (njlj2 == njlj) .and. (nklk2 == nklk) .and. (nlll2 == nlll) .and. ((ih+i) < nh))
            i = i + 1
            l_list(i) = l
            call vhartree_decode_i(lut_v_hartree(ih+i),nili2,njlj2,nklk2,nlll2,l)
          end do
          if (i == 0) then
            i = 1
            l_list(1) = l
          end if
          nili2 = nili ; li = pdr%l_value(nili) ; li2 = li
          njlj2 = njlj ; lj = pdr%l_value(njlj) ; lj2 = lj
          nklk2 = nklk ; lk = pdr%l_value(nklk) ; lk2 = lk
          nlll2 = nlll ; ll = pdr%l_value(nlll) ; ll2 = ll
          do loop = 1,4
            ok = .false.
            if (loop == 1) then
              ok = .true.
            elseif ((loop == 2) .and. (nili /= njlj)) then
              ok = .true.
              nili = njlj2 ; li = lj2
              njlj = nili2 ; lj = li2
              nklk = nklk2 ; lk = lk2
              nlll = nlll2 ; ll = ll2
            elseif ((loop == 3) .and. (nklk /= nlll)) then
              ok = .true.
              nili = nili2 ; li = li2
              njlj = njlj2 ; lj = lj2
              nklk = nlll2 ; lk = ll2
              nlll = nklk2 ; ll = lk2
            elseif ((loop == 4) .and. (nili /= njlj) .and. (nklk /= nlll)) then
              ok = .true.
              nili = njlj2 ; li = lj2
              njlj = nili2 ; lj = li2
              nklk = nlll2 ; lk = ll2
              nlll = nklk2 ; ll = lk2
            end if
            if (ok) then
              do mi = -li,+li
                do mj = -lj,+lj
                  m = mi - mj        
                  msign = (-1.0_double)**m
                  do mk = -lk,+lk
                    do ml = -ll,+ll
                      if (m == (ml-mk)) then    
                        term = 0.0_double
                        do j = 1,i
                          if (abs(m) <= l_list(j)) then
                            gc1 = gaunt_complex(l_list(j),-m,li,mi,lj,mj)
                            gc2 = gaunt_complex(l_list(j),m,lk,mk,ll,ml)
                            term = term + v_hartree(ih+j-1)*gc1*gc2
                          end if
                        end do
                        if (abs(term) > 1.0e-20_double) then
                          n = n + 1
                          pdr%cijkl(n) = term*msign
                          pdr%lut_cijkl(n) = cijkl_encode_i(nili,njlj,nklk,nlll,mi,mj,mk)
                        end if
                      end if
                    end do
                  end do
                end do
              end do
            end if
          end do
          ih = ih + i
        end do
100     if (allocated( i1d )) deallocate( i1d )
        if (allocated( lut )) deallocate( lut )
        if (allocated( l_list )) deallocate( l_list )
        if (allocated( r1d )) deallocate( r1d )
      end subroutine 
 
      subroutine orbital_matrix_i(pdr,density,v_hat)
        type(paw_data_rep) :: pdr
        real(double), dimension(:) :: density, v_hat
!       requires: pdr%aqlm, pdr%avlm, and pdr%lut_orbital be nullified or associated.
  
        integer :: ic, n
        integer :: i, j, li, lj, mi, mj, l, m, nili, njlj 
        real(double) :: gc
        n = 0
        do ic = 1,size(pdr%lut_denvhat)
          call denvhat_decode_i(pdr,ic,nili,njlj,l)
          if (nili == njlj) then
            j = 1
          else
            j = 2
          end if
          do i = 1,j
            if (i == 2) then
              m = nili ; nili = njlj ; njlj = m
            end if
            li = pdr%l_value(nili)
            lj = pdr%l_value(njlj)
            do mi = -li,+li
              do mj = -lj,+lj
                m = mj - mi
                if (abs(m) <= l) n = n + 1
              end do
            end do
          end do
        end do
        if (associated( pdr%aqlm )) deallocate( pdr%aqlm ) ; allocate( pdr%aqlm(n) )
        if (associated( pdr%avlm )) deallocate( pdr%avlm ) ; allocate( pdr%avlm(n) )
        if (associated( pdr%lut_orbital )) deallocate( pdr%lut_orbital ) ; allocate( pdr%lut_orbital(n) )
        n = 0
        do ic = 1,size(pdr%lut_denvhat)
          call denvhat_decode_i(pdr,ic,nili,njlj,l)
          if (nili == njlj) then
            j = 1
          else
            j = 2
          end if
          do i = 1,j
            if (i == 2) then
              m = nili ; nili = njlj ; njlj = m
            end if
            li = pdr%l_value(nili)
            lj = pdr%l_value(njlj)
            do mi = -li,+li
              do mj = -lj,+lj
                m = mj - mi
                if (abs(m) <= l) then
                  n = n + 1
                  pdr%lut_orbital(n) = orbital_encode_i(nili,njlj,mi,mj,l)
                  gc = gaunt_complex(l,m,li,mi,lj,mj)
                  pdr%aqlm(n) = gc*density(ic)
                  pdr%avlm(n) = gc*v_hat(ic)
                  if (odd(m)) pdr%avlm(n) = -pdr%avlm(n)
                end if
              end do
            end do
          end do
        end do
      end subroutine

      function cijkl_encode_i(nili,njlj,nklk,nlll,mi,mj,mk) result(n)
        integer, intent(in) :: nili, njlj, nklk, nlll, mi, mj, mk
        integer :: n

        n = nili
        n = ishft(n,5) + njlj
        n = ishft(n,5) + nklk
        n = ishft(n,5) + nlll
        n = ishft(n,3) + mi + 3
        n = ishft(n,3) + mj + 3
        n = ishft(n,3) + mk + 3
      end function
 
      function denvhat_encode_i(nili,njlj,l) result(n)
        integer, intent(in) :: nili, njlj, l
        integer :: n

        n = nili
        n = ishft(n,5) + njlj
        n = ishft(n,3) + l
      end function

      function vhartree_encode_i(nili,njlj,nklk,nlll,l) result(n)
        integer, intent(in) :: nili, njlj, nklk, nlll, l
        integer :: n

        n = nili
        n = ishft(n,5) + njlj
        n = ishft(n,5) + nklk
        n = ishft(n,5) + nlll
        n = ishft(n,3) + l
      end function

      function orbital_encode_i(nili,njlj,mi,mj,l) result(n)
        integer, intent(in) :: nili, njlj, mi, mj, l
        integer :: n
        n = nili

        n = ishft(n,5) + njlj
        n = ishft(n,3) + mi + 3
        n = ishft(n,3) + mj + 3
        n = ishft(n,3) + l
      end function
 
      subroutine vhartree_decode_i(n,nili,njlj,nklk,nlll,l)
        integer, intent(in) :: n
        integer, intent(out) :: nili, njlj, nklk, nlll, l

        integer, parameter :: vd_p18 = 8126464
        integer, parameter :: vd_p13 = 253952
        integer, parameter :: vd_p8 = 7936
        integer, parameter :: vd_p3 = 248
        nili = ishft(iand(vd_p18,n),-18)
        njlj = ishft(iand(vd_p13,n),-13)
        nklk = ishft(iand(vd_p8,n),-8)
        nlll = ishft(iand(vd_p3,n),-3)
        l = iand(7,n)
      end subroutine

      function linear_search_i(a,match_value,exact) result(best_index)
        integer, dimension(:), intent(in) :: a
        integer, intent(in) :: match_value
        logical, optional, intent(out) :: exact
        integer :: best_index

        logical :: found
        integer :: i, best_value, best_miss
        best_value = a(1)
        best_index = 1
        best_miss = abs(best_value - match_value)
        i = 2
        found = .false.
        do while ((i <= size(a)) .and. (best_value /= match_value))
          if (abs(a(i) - match_value) < best_miss) then
            best_index = i
            best_value = a(i)
            best_miss = abs(best_value - match_value)
          end if
          i = i + 1
        end do
        if (present(exact)) then
          if (best_value == match_value) then
            exact = .true.
          else
            exact = .false.
          end if
        end if
      end function

      subroutine insertion_sort_i(a,lut)
        integer, dimension(:), intent(in) :: a
        integer, dimension(:), intent(inout) :: lut

        integer :: i, j, cv
        do i = 1,size(lut)
          lut(i) = i
        end do
        do i = 2,size(a)
          cv = a(lut(i))
          j = i - 1
          do while ((cv < a(lut(j))) .and. (j > 1))
            lut(j+1) = lut(j)
            j = j - 1
          end do
          if (cv < a(lut(j))) then
            lut(j+1) = lut(j)
            j = j - 1
          end if
          lut(j+1) = i
        end do
      end subroutine

      subroutine open_word_i(f,wc)
        type(file_obj) :: f
        type(word_context), intent(out) :: wc

        w_ok = IOSTAT_OK
        w_eol = IOSTAT_EOL
        w_eof = IOSTAT_EOF
        wc%error = w_ok
        wc%delimiter = " ()[]{},"
        wc%comment = "#"
        wc%literal = "'"
        wc%column = 0
        wc%text = " "
        wc%length = 0
        call my(f,wc%f)
        if (error("Exit paw_data_mod::open_word_i")) continue
      end subroutine

      subroutine close_word_i(wc)
        type(word_context), intent(out) :: wc

        call glean(thy(wc%f))
      end subroutine

      subroutine get_next_word_i(wc,token)
        type(word_context), intent(inout) :: wc
        character(*), intent(inout) :: token

        integer :: length
        if (i_access(wc%f)) then
          call get_word_i(wc,token,length)
          do while ((length == 0) .and. (wc%error /= w_eof))
            call get_word_i(wc,token,length)
          end do
          if (wc%error == w_eol) wc%error = w_ok
          token = trim(token)
        end if
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,token)
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,wc%error)
      end subroutine

      subroutine get_next_word_i_nb(wc,token,length)
        type(word_context), intent(inout) :: wc
        character(*), intent(inout) :: token
        integer, intent(out) :: length

        if (i_access(wc%f)) then
          call get_word_i(wc,token,length)
          do while ((length == 0) .and. (wc%error /= w_eof))
            call get_word_i(wc,token,length)
          end do
          if (wc%error == w_eol) wc%error = w_ok
          token = trim(token)
        end if
      end subroutine

      subroutine get_word_i(wc,token,length)
        type(word_context), intent(inout) :: wc
        character(*), intent(inout) :: token
        integer, intent(out) :: length

        logical :: finished
        character :: c
        character(len(token)) :: temp
        integer :: k
        call next_word_i(wc,token,length)
        if (wc%error /= w_eof) then
          finished = .false.
          do while (.not.finished)
            finished = .true.
            temp = token
            if (index(token,wc%comment) /= 0) then
              do while (index(token,wc%comment) /= 0)
                do while (wc%error == w_ok)
                  call get_c_i(wc,c)
                end do
                k = index(token,wc%comment)
                if (k < 2) then
                  call next_word_i(wc,token,length)
                else
                  token(k:length) = " "
                  length = length - k
                end if
                finished = .false.
              end do
            end if
          end do
        end if
      end subroutine

      subroutine next_word_i(wc,token,length)
        type(word_context), intent(inout) :: wc
        character(*), intent(inout) :: token
        integer, intent(out) :: length

        character :: endchar
        wc%error = w_ok
        length = 0
        do while ((length == 0) .and. (wc%error /= w_eof))
          length = get_token_i(wc,token,endchar)
        end do
      end subroutine

      function get_token_i(wc,token,endchar) result(t)
        type(word_context), intent(inout) :: wc
        character(*), intent(inout) :: token
        character(*), intent(inout) :: endchar
        integer :: t

        character :: c
        integer :: pos

        pos = 0
        token = " "
        call get_c_i(wc,c)
        do while ((scan(wc%delimiter,c) /= 0) .and. (wc%literal /= c) .and. (wc%error == w_ok))
          call get_c_i(wc,c)
        end do
        if (wc%error == w_ok) then
          if (wc%literal == c) then
            call get_c_i(wc,c)
            do while ((wc%literal /= c) .and. (wc%error == w_ok))
              pos = pos + 1
              token(pos:pos) = c
              call get_c_i(wc,c)
            end do
          else
            do while ((scan(wc%delimiter,c) == 0) .and. (wc%error == w_ok))
              pos = pos + 1
              token(pos:pos) = c
              call get_c_i(wc,c)
            end do
          end if
        end if
        endchar = c
        t = pos
      end function

      subroutine get_c_i(wc,c)
        type(word_context), intent(inout) :: wc
        character, intent(out) :: c

        wc%error = w_ok
        if (wc%column >= wc%length) then
          call load_i(x_unit(wc%f),wc%text,wc%length,wc%last_error)
          wc%column = 0
        end if
        wc%column = wc%column + 1
        c = wc%text(wc%column:wc%column)
        if (wc%column >= wc%length) wc%error = wc%last_error
      end subroutine

      subroutine load_i(unit,text,length,last_error)
        integer, intent(in) :: unit
        character(*), intent(out) :: text
        integer, intent(out) :: length, last_error
        character :: c

        length = 0
        last_error = w_ok
        do while (last_error == w_ok)
          read(unit,'(a1)',advance='no',iostat=last_error) c
          length = length + 1
          text(length:length) = c
        end do
      end subroutine

      subroutine get_integer_i(wc,a)
        type(word_context), intent(inout) :: wc
        integer, intent(out) :: a

        character(line_len) :: token
        integer :: length
        call get_next_word_i_nb(wc,token,length)
        if (length > 0) read(token,*,iostat=wc%error) a
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,a)
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,wc%error)
      end subroutine

      subroutine get_integers_i(wc,a)
        type(word_context), intent(inout) :: wc
        integer, dimension(:), intent(out) :: a

        character(line_len) :: token
        integer :: i, length
        do i = 1,size(a)
          call get_next_word_i_nb(wc,token,length)
          if (length > 0) read(token,*,iostat=wc%error) a(i)
        end do
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,a)
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,wc%error)
      end subroutine

      subroutine get_real_i(wc,a)
        type(word_context), intent(inout) :: wc
        real(double), intent(out) :: a

        character(line_len) :: token
        integer :: length
        call get_next_word_i_nb(wc,token,length)
        if (length > 0) read(token,*,iostat=wc%error) a
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,a)
        if (i_comm(wc%f)) call broadcast(FILE_SCOPE,wc%error)
      end subroutine

      subroutine get_reals_i(wc,a,l)
        type(word_context), intent(inout) :: wc
        real(double), dimension(:), intent(out) :: a
        logical, optional :: l

        character(line_len) :: token
        integer :: i, length
        do i = 1,size(a)
          call get_next_word_i_nb(wc,token,length)
          if (length > 0) read(token,*,iostat=wc%error) a(i)
        end do
        if (present(l)) then
          if (l) then
            if (i_comm(wc%f)) call broadcast(FILE_SCOPE,a)
            if (i_comm(wc%f)) call broadcast(FILE_SCOPE,wc%error)
          end if
        else
          if (i_comm(wc%f)) call broadcast(FILE_SCOPE,a)
          if (i_comm(wc%f)) call broadcast(FILE_SCOPE,wc%error)
        end if
      end subroutine

      subroutine denvhat_decode_i(pdr,i,nili,njlj,l)
        type(paw_data_rep) :: pdr
        integer, intent(in) :: i
        integer, intent(out) :: nili, njlj, l
!       requires: i be in range.

        integer, parameter :: dd_p8 = 7936
        integer, parameter :: dd_p3 = 248
        integer :: n
        n = pdr%lut_denvhat(i)
        nili = ishft(iand(dd_p8,n),-8)
        njlj = ishft(iand(dd_p3,n),-3)
        l = iand(7,n)
      end subroutine

      function radial_f_value_i(f,r,g,l,switch) result(rfv)
        real(double), dimension(:), intent(in) :: f, r
        real(double), intent(in) :: g
        integer, intent(in) :: l
        logical, intent(in), optional :: switch
        real(double) :: rfv

        if (present(switch)) then
          rfv = sum(f*spherical_bessel(g*r,l,switch))
        else
          rfv = sum(f*spherical_bessel(g*r,l))
        end if
        if (error("Exit paw_data_mod::radial_f_value_i")) continue
      end function

      function radial_r_value_i(f,g,r,l) result(rrv)
        real(double), dimension(:), intent(in) :: f, g
        real(double), intent(in) :: r
        integer, intent(in) :: l
        real(double) :: rrv

        rrv = sum(f*spherical_bessel(g*r,l))
        if (error("Exit paw_data_mod::radial_r_value_i")) continue
      end function

      subroutine optimize_nlp_i(pdr)
        type(paw_data_rep) :: pdr

        integer, parameter :: nr = 1000
        logical :: found
        integer :: ib, ic, ig, ir, l, ll, nb, ng_d, ng_i, ng_o, ns
        real(double) :: dg, dr, r0, r1
        real(double), dimension(:), allocatable :: b, b0, b1, f, gr, r, wg, wt
        real(double), dimension(:,:), allocatable :: a, t

        call arg("optimization_points",ng_o,found)
        if (.not.found) ng_o = 181

        allocate( pdr%g(ng_o) )  ! GENERATE A RADIAL G-VECTOR_MESH.
        dg = pdr%g_ocut/real(ng_o,double)
        do ig = 1,ng_o
          pdr%g(ig) = real(ig-1,double)*dg
        end do

        ng_i = int(pdr%g_icut/dg)
        ng_d = ng_o - ng_i

        nb = size(pdr%l_value)

        allocate( pdr%tpfo(ng_o,nb) )

        do ib = 1,nb  ! COMPUTE FOURIER COMPONENTS UP TO GI_CUT
          l = pdr%l_value(ib)
          do ig = 1,ng_i
            pdr%tpfo(ig,ib) = radial_f_value_i(pdr%tp(:,ib),pdr%r_sg,pdr%g(ig),l)
          end do
        end do

        allocate( a(ng_d,ng_d), b(ng_d), t(ng_o,ng_o) )  ! COMPUTE FOURIER COMPONENTS BETWEEN GI_CUT AND GO_CUT.
        allocate( gr(ng_o), b0(ng_o), b1(ng_o) )
        gr = pdr%g*pdr%r_opt
        r0 = 1.0_double/pdr%r_opt
        do ib = 1,nb
          l = pdr%l_value(ib)
          r1 = real(2*l+1,double)
          b0 = spherical_bessel(gr,l,.true.)
          b1 = spherical_bessel(gr,l+1,.true.)
          do ic = 1,ng_o
            do ir = 1,ic-1
              t(ir,ic) = r0*(gr(ir)*gr(ic))**(l+2)/(gr(ir)**2 - gr(ic)**2)*(gr(ir)**2*b1(ir)*b0(ic) - gr(ic)**2*b0(ir)*b1(ic))
              t(ic,ir) = t(ir,ic)
            end do
            t(ic,ic) = 0.5_double*r0*gr(ic)**(2*l+4)*(b0(ic)**2 + (gr(ic)*b1(ic))**2 - r1*b0(ic)*b1(ic))
          end do
          do ir = 1,ng_d
            b(ir) = sum(pdr%tpfo(1:ng_i,ib)*t(1:ng_i,ng_i+ir)*dg)
            a(ir,:) = -t(ng_i+1:ng_o,ng_i+ir)*dg
            a(ir,ir) = a(ir,ir) + 0.5_double*pi*pdr%g(ng_i+ir)**2
          end do
          pdr%tpfo(ng_i+1:ng_o,ib) = solve(a,b) ; if (error()) goto 100
        end do
        deallocate( a, b, t )
        deallocate( gr, b0, b1 )

        allocate( t(ng_o,nb) )  ! CONSTRUCT INTEGRANDS USED TO COMPUTE THE REAL-SPACE PROJECTORS.
        allocate( wt(ng_o) )
        if (mod(ng_o,2) == 0) then
          ns = 4
          wt(1) = 3.0_double/8.0_double
          wt(2) = 9.0_double/8.0_double
          wt(3) = 9.0_double/8.0_double
          wt(4) = 3.0_double/8.0_double + 1.0_double/3.0_double
        elseif (mod(ng_o,2) == 1) then
          ns = 1
          wt(1) = 1.0_double/3.0_double
        end if
        ll = 4
        do ig = ns+1,ng_o-1
          wt(ig) = real(ll,double)/3.0_double
          ll = 6 - ll
        end do
        wt(ng_o) = 1.0_double/3.0_double
        do ib = 1,nb
          t(1:ng_i,ib) = pdr%tpfo(1:ng_i,ib)
          pdr%tpfo(:,ib) = (2.0_double/pi)*pdr%tpfo(:,ib)*(pdr%g**2)*wt*dg
        end do
        deallocate( wt )

        allocate( r(nr), wt(nr) )  ! COMPUTE THE MAXIMUM ERRORS FOR THE REAL-SPACE PROJECTORS.
        dr = pdr%r_opt/real(nr,double)
        do ir = 1,nr
          r(ir) = dr*real(ir-1,double)
        end do
        if (mod(nr,2) == 0) then
          ns = 4
          wt(1) = 3.0_double/8.0_double
          wt(2) = 9.0_double/8.0_double
          wt(3) = 9.0_double/8.0_double
          wt(4) = 3.0_double/8.0_double + 1.0_double/3.0_double
        elseif (mod(nr,2) == 1) then
          ns = 1
          wt(1) = 1.0_double/3.0_double
        end if
        ll = 4
        do ir = ns+1,nr-1
          wt(ir) = real(ll,double)/3.0_double
          ll = 6 - ll
        end do
        wt(nr) = 1.0_double/3.0_double
        allocate( pdr%w_max(nb) )
        allocate( f(nr), wg(ng_i) )
        do ib = 1,nb
          do ir = 1,nr
            f(ir) = radial_r_value_i(pdr%tpfo(:,ib),pdr%g,r(ir),pdr%l_value(ib))*r(ir)**2*wt(ir)*dr ; if (error()) goto 100
          end do
          do ig = 1,ng_i
            wg(ig) = radial_f_value_i(f,r,pdr%g(ig),pdr%l_value(ib)) - t(ig,ib) ; if (error()) goto 100
          end do
          pdr%w_max(ib) = maxval(abs(wg))
        end do
        deallocate( f, wg )
        deallocate( r, wt )
        deallocate( t )

100     if (allocated( f )) deallocate( f )
        if (allocated( r )) deallocate( r )
        if (allocated( a )) deallocate( a )
        if (allocated( b )) deallocate( b )
        if (allocated( t )) deallocate( t )
        if (allocated( gr )) deallocate( gr )
        if (allocated( b0 )) deallocate( b0 )
        if (allocated( b1 )) deallocate( b1 )
        if (allocated( wt )) deallocate( wt )
        if (allocated( wg )) deallocate( wg )

        if (error("Exit paw_data_mod::optimize_nlp_i")) continue

      end subroutine

      end module

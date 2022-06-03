!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module ncp_data_mod
!doc$ module ncp_data_mod

!     One datatype is available here: type(ncp_data_obj).

!     This module encapsulates radial functions supplied by the various types of norm-conserving
!     pseudopotentials, converting them into forms required by atomic_operators_ncp_mod.

      use kind_mod
      use path_mod
      use mpi_mod
      use error_mod
      use io_mod
      use math_mod
      use arg_mod
      use diary_mod

!cod$
      implicit none
      private 

      integer, parameter :: MOD_SCOPE = CONFIG

      real(double), parameter :: e2 = 2.0_double
      real(double), parameter :: tol_g_zero = 1.0e-10_double

      type :: nonlocal_hgh_data
         integer :: c                                            ! number of l channels
         real(double) :: r                                            ! l-indexed radius for projectors
         real(double), dimension(:,:), pointer :: h                            ! projector coefficients for SR piece
         real(double), dimension(:,:), pointer :: k                            ! projector coefficients for SO piece
      end type


      type :: ncp_data_rep
        integer :: ref
        character(tag_sz) :: name                     ! name associated with data
        real(double), dimension(:), pointer :: r      ! real-space radial grid
        real(double) :: exponent                      ! exponent for the long-ranged part of the local pseudopotential
        real(double), dimension(:), pointer :: lp     ! short-ranged part of the local pseudopotential on the real-space grid
        real(double), dimension(:), pointer :: vd     ! valence density on the real-space grid (optional)
        real(double), dimension(:), pointer :: cd     ! core density on the real-space grid (optional)
        real(double) :: r_cut                         ! cutoff radius for the non-local pseudopotentials
        integer, dimension(:), pointer :: l           ! angular momenta for the non-local pseudopotentials
        real(double), dimension(:), pointer :: kb     ! Kleinman-Bylander parameters for the non-local pseudopotentials
        real(double), dimension(:,:), pointer :: nlp  ! non-local pseudopotentials on the real-space grid
        real(double) :: r_opt                         ! radius for the optimized non-local pseudopotentials
        real(double) :: gi_cut                        ! inner g_vector cutoff for the optimized non-local pseudopotentials
        real(double) :: go_cut                        ! outer g-vector cutoff for the optimized non-local pseudopotentials
        real(double), dimension(:), pointer :: w_max  ! maximum error for the optimized non-local pseudopotentials
        real(double), dimension(:), pointer :: g      ! reciprocal-space radial grid 
        real(double), dimension(:,:), pointer :: nlpo ! Fourier coefficients for the optimized non-local pseudopotentials
        real(double) :: valence                       ! valence of atom
        type(nonlocal_hgh_data), dimension(:), pointer :: nl           ! nonlocal index for hgh pseudopotentials
        integer :: nnonloc                                            ! number of nonlocal(l) index
        real(double) :: rloc                                           ! local radius for hgh
        real(double), dimension(:), pointer :: cn                               ! expansion coefficients for local hgh potential
        integer :: zatom                                               ! total atomic charge from hgh potential file
        integer :: pspcod                                              ! pspcod info from hgh potential file
        integer :: pspxc                                               ! pspxc info from hgh potential file
        integer :: lmax                                                ! lmax info from hgh potential file
        character(len("Goedecker")) :: id                       ! id tag to determine if reading in hgh file
     end type

      type, public :: ncp_data_obj
        private
        integer ref
        type(ncp_data_rep), pointer :: o
      end type

      character(len("simple")), parameter    :: SIMPLE_STRING    = "simple"
      character(len("linear")), parameter    :: LINEAR_STRING    = "linear"
      character(len("linear_cc")), parameter :: LINEAR_CC_STRING = "linear_cc"
      character(len("log")), parameter       :: LOG_STRING       = "log"
      character(len("log_cc")), parameter    :: LOG_CC_STRING    = "log_cc"
      character(len("log_wt")), parameter    :: LOG_WT_STRING    = "log_wt"
      character(len("log_wt_cc")), parameter :: LOG_WT_CC_STRING = "log_wt_cc"

!doc$
      public :: ncp_data
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_name
      public :: x_nlcomps
      public :: x_lcomp
      public :: x_kbparam
      public :: x_l_channels
      public :: x_valence
      public :: x_valence_electrons
      public :: x_radius
      public :: x_inner_cutoff
      public :: x_outer_cutoff
      public :: x_error
      public :: core_density
      public :: diary_rs_projectors
      public :: valence_f_value
      public :: local_f_value
      public :: core_f_value
      public :: local_stress_f_value
      public :: core_stress_f_value
      public :: apply_projector_kbf
      public :: projector_f_values
      public :: projector_stress_f_values
      public :: projector_r_value
      public :: projector_r_gradients

!cod$
      interface ncp_data
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
      interface x_nlcomps
        module procedure pd_nlcomps
      end interface
      interface x_lcomp
        module procedure pd_lcomp, pd_lcomps
      end interface
      interface x_kbparam
        module procedure pd_kbparam, pd_kbparams
      end interface
      interface x_l_channels
         module procedure pd_l_channels
      end interface
      interface x_valence
        module procedure pd_valence
      end interface
      interface x_valence_electrons
         module procedure pd_valence_electrons
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
      interface core_density
        module procedure pd_core_density
      end interface
      interface diary_rs_projectors
        module procedure diary_rs_projectors_pd
      end interface
      interface valence_f_value
        module procedure valence_f_value_pd
      end interface
      interface local_f_value
        module procedure local_f_value_pd
      end interface
      interface core_f_value
        module procedure core_f_value_pd
      end interface
      interface local_stress_f_value
        module procedure local_stress_f_value_pd
      end interface
      interface core_stress_f_value
        module procedure core_stress_f_value_pd
      end interface
      interface apply_projector_kbf
         module procedure pd_apply_projector_kbf, pd_apply_projector_kbf_wij
      end interface
      interface projector_f_values
        module procedure projector_f_values_pd
      end interface
      interface projector_stress_f_values
        module procedure projector_stress_f_values_pd
      end interface
      interface projector_r_value
        module procedure projector_r_value_pd
      end interface
      interface projector_r_gradients
        module procedure projector_r_gradients_pd
      end interface

      contains

      function constructor_pd(tag) result(pd)
!doc$ function ncp_data(f,tag) result(pd)
        character(tag_sz), intent(in) :: tag
        type(ncp_data_obj) :: pd
!        requires: f be open and pointing to a ncp data block.
!        effects: Creates a new pd.
!        errors: Format errors.

!cod$

        type(file_obj) :: f
        logical :: found
        character(line_len) :: format_string, first_string
        character(17+tag_sz) :: pr_tag
        integer :: ios,i

        !call my(f)

        pd%ref = 0
        allocate( pd%o )
        pd%o%ref = 0
        

        pd%o%name = tag
        call my(file(trim(ncp_path)//pd%o%name),f)

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: problem opening file = "//x_name(f))) goto 100
        if (i_access(f)) read(x_unit(f), *, iostat=ios) pd%o%id
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to pp id line")) goto 100
        if (i_comm(f)) call broadcast(FILE_SCOPE,pd%o%id)
        if (i_access( diaryfile() )) then
        !write(x_unit(diaryfile()), '(/,"Valence of atom is =", a12)') pd%o%id
        end if
        if (pd%o%id == 'Goedecker') then

          nullify ( pd%o%l )
          nullify ( pd%o%r )
          nullify ( pd%o%lp )
          nullify ( pd%o%vd )
          nullify ( pd%o%cd )
          nullify ( pd%o%l )
          nullify ( pd%o%kb )
          nullify ( pd%o%nlp )
          nullify ( pd%o%w_max )
          nullify ( pd%o%g )
          nullify ( pd%o%nlpo )

          call read_hgh_i(pd%o,f); if (error()) goto 100

        else

          nullify ( pd%o%cn )
          nullify ( pd%o%nl )

          read (pd%o%id, *, iostat=ios) pd%o%valence
          if (error(ios /= 0, "ERROR: failed string conversion to double")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,pd%o%valence)

          if (i_access(f)) read(unit=x_unit(f),fmt=*,iostat=ios) format_string
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to read format_string")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,format_string)

          select case (trim(format_string))
          case (SIMPLE_STRING)
            call read_simple_i(pd%o,f) ; if (error()) goto 100
          case (LINEAR_STRING)
            call read_linear_i(pd%o,f) ; if (error()) goto 100
          case (LINEAR_CC_STRING)
            call read_linear_cc_i(pd%o,f) ; if (error()) goto 100
          case (LOG_STRING)
            call read_log_i(pd%o,f) ; if (error()) goto 100
          case (LOG_CC_STRING)
            call read_log_cc_i(pd%o,f) ; if (error()) goto 100
          case (LOG_WT_STRING)
            call read_log_wt_i(pd%o,f) ; if (error()) goto 100
          case (LOG_WT_CC_STRING)
            call read_log_wt_cc_i(pd%o,f) ; if (error()) goto 100
          case default
            if (error(.true.,"ERROR: unrecognized format_string")) goto 100
          end select

        end if

        pd%o%name = tag
        if (error(len_trim(pd%o%name) > 6,"ERROR: Atom name has more than six characters")) goto 100

        pr_tag = "projector_radius_"//pd%o%name
        call arg(trim(pr_tag),pd%o%r_opt,found)
        if (.not.found) pd%o%r_opt = 4.0_double
        
        if (i_access(f)) close(x_unit(f))

100     call glean(thy(f))

        if (error("Exit ncp_data_mod::constructor_pd")) continue

      end function

      subroutine my_pd(pd)
!doc$ subroutine my(pd)
        type(ncp_data_obj) :: pd

!cod$
        pd%ref = pd%ref + 1
        pd%o%ref = pd%o%ref + 1
      end subroutine

      subroutine my_new_pd(pdi,pd)
!doc$ subroutine my(pdi,pd)
        type(ncp_data_obj) :: pdi, pd

!cod$
        pd%ref = 1
        pd%o => pdi%o
        pd%o%ref = pd%o%ref + 1
      end subroutine

      function thy_pd(pd)
!doc$ function thy(pd)
        type(ncp_data_obj) :: pd, thy_pd

!cod$
        pd%ref = pd%ref - 1
        pd%o%ref = pd%o%ref - 1
        thy_pd%ref = pd%ref
        thy_pd%o => pd%o
      end function

      subroutine glean_pd(pd)
!doc$ subroutine glean(pd)
        type(ncp_data_obj) :: pd

!cod$
        integer :: i
        if (pd%o%ref < 1) then
          if (associated( pd%o%nl )) then
             do i = 1,size( pd%o%nl )
                 if (associated( pd%o%nl(i)%h )) deallocate( pd%o%nl(i)%h )
                 if (associated( pd%o%nl(i)%k )) deallocate( pd%o%nl(i)%k )
             end do
             deallocate( pd%o%nl )
          end if
          if (associated( pd%o%cn )) deallocate( pd%o%cn )   
          if (associated( pd%o%r )) deallocate( pd%o%r )
          if (associated( pd%o%lp )) deallocate( pd%o%lp )
          if (associated( pd%o%vd )) deallocate( pd%o%vd )
          if (associated( pd%o%cd )) deallocate( pd%o%cd )
          if (associated( pd%o%l )) deallocate( pd%o%l )
          if (associated( pd%o%kb )) deallocate( pd%o%kb )
          if (associated( pd%o%nlp )) deallocate( pd%o%nlp )
          if (associated( pd%o%w_max )) deallocate( pd%o%w_max )
          if (associated( pd%o%g )) deallocate( pd%o%g )
          if (associated( pd%o%nlpo )) deallocate( pd%o%nlpo )
          deallocate( pd%o )
        end if
      end subroutine
      
      subroutine bequeath_pd(pd)
!doc$ subroutine bequeath(pd)
        type(ncp_data_obj) :: pd

!cod$
        continue
      end subroutine

      subroutine assign_pd(pd,pd2)
!doc$ subroutine assignment(=)(pd,pd2)
        type(ncp_data_obj), intent(inout) :: pd
        type(ncp_data_obj), intent(in) :: pd2

!cod$
        type(ncp_data_obj) :: pdt
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
        type(ncp_data_obj) :: pd
        integer, dimension(2) :: r
!       effects: Returns pd%ref and pd%o%ref.

!cod$
        r(1) = pd%ref
        r(2) = pd%o%ref
        call glean(pd)
      end function

      function pd_name(pd) result(n)
!doc$ function x_name(pd) result(n)
        type(ncp_data_obj) :: pd
        character(tag_sz) :: n
!       effects: Returns the name associated with pd.

!cod$
        call my(pd)
        n = pd%o%name
        call glean(thy(pd))
      end function

      function pd_nlcomps(pd) result(nl)
!doc$ function x_nlcomps(pd) result(nl)
        type(ncp_data_obj) :: pd
        integer :: nl
!       effects: Returns the number of angular momenta in pd.

!cod$
        call my(pd)
        if (associated( pd%o%l )) then
          nl = size(pd%o%l)
        else if (associated( pd%o%nl )) then
           nl = size(pd%o%nl)
        else
          nl = 0
        end if
        call glean(thy(pd))
      end function

      function pd_lcomp(pd,i) result(l)
!doc$ function x_lcomp(pd,i) result(l)
        type(ncp_data_obj) :: pd
        integer, intent(in) :: i
        integer :: l
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns the i'th angular momentum in pd.
!       errors: If i is out of bounds.

!cod$
        call my(pd)
        if (pd%o%id == 'Goedecker') then
           if (error((i < 1) .or. (i > size(pd%o%nl)),"ERROR: i is out of bounds"))then
           goto 100
           else
              l = i-1
         end if
        else
         if (error((i < 1) .or. (i > size(pd%o%l)),"ERROR: i is out of bounds"))then 
              goto 100
         else 
          l = pd%o%l(i)
         end if
        end if
        
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::pd_lcomp")) continue
      end function

      function pd_lcomps(pd) result(l)
!doc$ function x_lcomp(pd) result(l)
        type(ncp_data_obj) :: pd
        integer, dimension(size(pd%o%l)) :: l
!       requires: x_nlcomps (pd) be > 0.
!       effects: Returns the angular momenta in pd.

!cod$
        call my(pd)
        l = pd%o%l
        call glean(thy(pd))
      end function

      function pd_l_channels(pd,i) result(c)
!doc$ function x_l_channels(pd,i) result(nl)
        type(ncp_data_obj) :: pd
        integer, intent(in) :: i
        integer :: c
!
!       effects: Returns the number of projector channels with the same l

        call my(pd)
        if (pd%o%id == 'Goedecker') then
           if (error((i < 1) .or. (i > size(pd%o%nl)), "ERROR: i is out of bounds")) then
              goto 100
              else
                 c = pd%o%nl(i)%c
           end if
        else
           c = 1
        end if
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::pd_l_channels")) continue
      end function

      function pd_kbparam(pd,i) result(kb)
!doc$ function x_kbparam(pd,i) result(kb)
        type(ncp_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: kb
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns the i'th Kleinman-Bylander parameter in pd.
!       errors: If i is out of bounds.

!cod$
        call my(pd)
        if (pd%o%id == 'Goedecker') then
           kb = 0.0_double
        else
        if (error((i < 1) .or. (i > size(pd%o%l)),"ERROR: i is out of bounds")) then
           goto 100
        end if
        kb = pd%o%kb(i)
        end if
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::pd_kbparam")) continue
      end function

      function pd_kbparams(pd) result(kb)
!doc$ function x_kbparam(pd) result(kb)
        type(ncp_data_obj) :: pd
        real(double), dimension(size(pd%o%kb)) :: kb
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns the Kleinman-Bylander parameters in pd.

!cod$
        call my(pd)
        kb = pd%o%kb
        call glean(thy(pd))
      end function

      function pd_valence(pd) result(v)
!doc$ function x_valence(pd) result(v)
        type(ncp_data_obj) :: pd
        logical :: v
!       effects: Returns .true. if pd contains valence density.

!cod$
        call my(pd)
        v = associated( pd%o%vd )
        call glean(thy(pd))
      end function

      function pd_valence_electrons(pd) result(ve)
!doc$ function x_valence_electrons(pd) result(ve)
        type(ncp_data_obj) :: pd
        real(double) :: ve
!       effects: Returns the valence of the atom, ve
        call my(pd)
        ve = pd%o%valence
        call glean(thy(pd))
      end function

      function pd_radius(pd) result(r)
!doc$ function x_radius(pd) result(r)
        type(ncp_data_obj) :: pd
        real(double) :: r
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns pd%o%r_opt in Bohr.

!cod$
        call my(pd)
        r = pd%o%r_opt
        call glean(thy(pd))
      end function

      function pd_inner_cutoff(pd) result(c)
!doc$ function x_inner_cutoff(pd) result(c)
        type(ncp_data_obj) :: pd
        real(double) :: c
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns pd%o%gi_cut in sqrt(Ryd).

!cod$
        call my(pd)
        c = pd%o%gi_cut
        call glean(thy(pd))
      end function

      function pd_outer_cutoff(pd) result(c)
!doc$ function x_outer_cutoff(pd) result(c)
        type(ncp_data_obj) :: pd
        real(double) :: c
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns pd%o%go_cut in sqrt(Ryd).

!cod$
        call my(pd)
        c = pd%o%go_cut
        call glean(thy(pd))
      end function

      function pd_error(pd,i) result(e)
!doc$ function x_error(pd,i) result(e)
        type(ncp_data_obj) :: pd
        integer, intent(in) :: i
        real(double) :: e
!       requires: x_nlcomps(pd) be > 0.
!       effects: Returns the i'th nlpo error in pd.
!       errors: If pd%o%w_max is not associated. If i is out of bounds.

!cod$
        call my(pd)
        if (error(.not.associated( pd%o%w_max ),"ERROR: w_max is not associated")) goto 100
        if (error((i < 1) .or. (i > size(pd%o%w_max)),"ERROR: i is out of bounds")) goto 100
        e = pd%o%w_max(i)
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::pd_error")) continue
      end function

      subroutine diary_rs_projectors_pd(pd)
!doc$ subroutine diary_rs_projectors(pd)
        type(ncp_data_obj) :: pd
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
          write(x_unit(diaryfile()),'(t8,"inner cutoff energy = ",f7.2," Ryd")') pd%o%gi_cut**2
          write(x_unit(diaryfile()),'(t8,"outer cutoff energy = ",f7.2," Ryd")') pd%o%go_cut**2
          do ib = 1,size(pd%o%l)
            write(x_unit(diaryfile()),'(t8,"l = ",i1,": maximum Fourier error = ",es8.2)') pd%o%l(ib), pd%o%w_max(ib)
          end do
        end if
        call glean(thy(pd))
      end subroutine

      function pd_core_density(pd) result(c)
!doc$ function core_density(pd) result(c)
        type(ncp_data_obj) :: pd
        logical :: c
!       effects: Returns .true. if pd contains a core density.

!cod$
        call my(pd)
        c = associated( pd%o%cd )
        call glean(thy(pd))
      end function

      function valence_f_value_pd(pd,g) result(f)
!doc$ function valence_f_value(pd,g) result(f)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: g
        real(double) :: f
!       requires: x_valence(pd) = .true.
!       effects: Returns the Fourier coefficient of a valence density.
!       errors: Passes errors.

!cod$
        call my(pd)
        f = radial_f_value_i(pd%o,pd%o%r*pd%o%vd,g,0)
        call glean(thy(pd))
        if (error("Exit ncp_data_mod::valence_f_value_pd")) continue
      end function

      function local_f_value_pd(pd,g,z) result(f)
!doc$ function local_f_value(pd,g,z) result(f)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: g, z
        real(double) :: f
!       effects: Returns the Fourier coefficient of the local pseudopotential.
!       errors: Passes errors.

!cod$
        real(double) :: c, sf
        integer :: i
        logical :: found

        call my(pd)
        if (pd%o%id == 'Goedecker') then
           if (g < tol_g_zero) then
              f = four_pi*z*pd%o%rloc**2/2
           else
              f = -four_pi*z*exp(-(g*pd%o%rloc)**2/2)/g**2
          end if
          c = sqrt(8*pi**3)*pd%o%rloc**3*exp(-(g*pd%o%rloc)**2/2)
           do i = 1,size(pd%o%cn)
              select case(i)
                 case(1)
                    f = f + c*pd%o%cn(1)
                 case(2)
                    f = f + c*pd%o%cn(2)*(3-g**2*pd%o%rloc**2)
                 case(3)
                    f = f + c*pd%o%cn(3)*(15-10*(g*pd%o%rloc)**2+(g*pd%o%rloc)**4)
                 case(4)
                    f = f + c*pd%o%cn(4)*(105-105*(g*pd%o%rloc)**2+21*(g*pd%o%rloc)**4-(g*pd%o%rloc)**6)
              end select      
           end do
          f = 2.0_double*f
        else
        if (g < tol_g_zero) then
          f = four_pi*z*e2/(4.0_double*pd%o%exponent)
        else
          f = -four_pi*z*e2*exp(-g**2/(4.0_double*pd%o%exponent))/g**2
        end if
        if (associated( pd%o%lp ))then
           f = f + radial_f_value_i(pd%o,pd%o%r*pd%o%lp,g,0)
        end if
        end if
        call glean(thy(pd))
        if (error("Exit ncp_data_mod::local_f_value_pd")) continue
      end function

      function core_f_value_pd(pd,g) result(f)
!doc$ function core_f_value(pd,g) result(f)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: g
        real(double) :: f
!       requires: x_core(pd) = .true.
!       effects: Returns the Fourier coefficient of a core density.
!       errors: Passes errors.

!cod$
        call my(pd)
        f = radial_f_value_i(pd%o,pd%o%r*pd%o%cd,g,0)
        call glean(thy(pd))
        if (error("Exit ncp_data_mod::core_f_value_pd")) continue
      end function

      function local_stress_f_value_pd(pd,g,z) result(f)
!doc$ function local_stress_f_value(pd,g,z) result(f)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: g, z
        real(double) :: f
!       effects: Returns the Fourier coefficient of a local pseudopotential contribution to the stress tensor.
!       errors: Passes errors.

!cod$
        real(double) :: r1, c
        integer :: i
        call my(pd)
         if (pd%o%id == 'Goedecker') then
          if (g < tol_g_zero) then
              f = 0.0_double
           else
              f = -four_pi*z*exp(-(g*pd%o%rloc)**2/2)*(2.0_double+g**2*pd%o%rloc**2)/g**4
          end if
          c = sqrt(8*pi**3)*pd%o%rloc**5*exp(-(g*pd%o%rloc)**2/2)
           do i = 1,size(pd%o%cn)
              select case(i)
                 case(1)
                    f = f + c*pd%o%cn(1)
                 case(2)
                    f = f + c*pd%o%cn(2)*(5-g**2*pd%o%rloc**2)
                 case(3)
                    f = f + c*pd%o%cn(3)*(35-14*(g*pd%o%rloc)**2+(g*pd%o%rloc)**4)
                 case(4)
                    f = f + c*pd%o%cn(4)*(315-189*(g*pd%o%rloc)**2+27*(g*pd%o%rloc)**4-(g*pd%o%rloc)**6)
              end select
           end do
           f = 2.0_double*f
        else
        if (g < tol_g_zero) then
          f = 0.0_double
        else
          r1 = g**2/(4.0_double*pd%o%exponent)
          f = -eight_pi*z*e2*(1.0_double + r1)*exp(-r1)/g**4
          if (associated( pd%o%lp )) f = f + radial_f_value_i(pd%o,pd%o%r**2*pd%o%lp,g,1)/g
        end if
        end if
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::local_stress_f_value_pd")) continue
      end function

      function core_stress_f_value_pd(pd,g) result(f)
!doc$ function core_stress_f_value(pd,g) result(f)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: g
        real(double) :: f
!       requires: x_core(pd) = .true.
!       effects: Returns the Fourier coefficient of a core density contribution to the stress tensor.
!       errors: Passes errors.

!cod$
        call my(pd)
        f = sum(pd%o%r**3*pd%o%cd*spherical_bessel(g*pd%o%r,1,.true.))
        call glean(thy(pd))
        if (error("Exit ncp_data_mod::core_stress_f_value_pd")) continue
      end function

      subroutine pd_apply_projector_kbf(pd, pdots, il)
!doc$ subroutine apply_projector_kbf(pd, pdots_kbf) result(kbf_pdots)
        type(ncp_data_obj) :: pd
        complex(double), dimension(:), intent(inout) :: pdots
        integer, intent(in) :: il
!
!
!

!cod$
        call my(pd)
        if (pd%o%id == 'Goedecker') then
           pdots = 1/2.0_double*matmul(pd%o%nl(il)%h,pdots)
        else
           pdots = pd%o%kb(il)*pdots
        end if
        call glean(thy(pd))
      end subroutine

      subroutine pd_apply_projector_kbf_wij(pd, wij_kbf, il)
!doc$ subroutine apply_projector_kbf(pd, pdots_kbf) result(kbf_pdots)
        type(ncp_data_obj) :: pd
        complex(double), dimension(:, :), intent(inout) :: wij_kbf
        integer, intent(in) :: il
!
!  matrix multiplication of the projector coefficients with the wij matrix
!

!cod$
        call my(pd)
        if (pd%o%id == 'Goedecker') then
           wij_kbf = 1/2.0_double*matmul(pd%o%nl(il)%h,wij_kbf)
        else
           wij_kbf = pd%o%kb(il)*wij_kbf
        end if
        call glean(thy(pd))
      end subroutine


      function projector_f_values_pd(pd,g,l, c) result(f)
!doc$ function projector_f_values(pd,g,l) result(f)
        type(ncp_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        integer, intent(in) :: l, c
        real(double), dimension(size(g)) :: f 
!       requires: x_nlcomps(pd) be > 0. l be one of pd%o%l(i).
!       effects: Returns the g Fourier coefficients of pd%o%nlp for channel c of angular momentum l.
!       errors: Passes errors.

!doc$
        integer :: i
        real(double) :: f_sum, f_ssum
        call my(pd)
        if (pd%o%id == 'Goedecker') then
           do i = 1,size(g)
              f(i) = 2.0_double*hgh_nonlocal_radial_f_value_i(pd%o,g(i),l,c)
           end do
        else   
        do i = 1,size(g)
          f(i) = nonlocal_radial_f_value_i(pd%o,g(i),l)
        end do
        end if
        call glean(thy(pd))
        if (error("Exit ncp_data_mod::projector_f_values_pd")) continue
      end function

      function projector_stress_f_values_pd(pd,g,l,c) result(f)
!doc$ function projector_stress_f_values(pd,g,l) result(f)
        type(ncp_data_obj) :: pd
        real(double), dimension(:), intent(in) :: g
        integer, intent(in) :: l, c
        real(double), dimension(2,size(g)) :: f
!       requires: x_nlcomps(pd) be > 0. l be one of pd%o%l(i).
!       effects: Returns the g Fourier coefficients of the projector contribution to the stress tensor.
!       errors: Passes errors.

!doc$
        integer :: i
        call my(pd)
        if (pd%o%id == 'Goedecker') then
           do i = 1,size(g)
            f(1,i) = 2.0_double*hgh_nonlocal_radial_f_value_i(pd%o,g(i),l,c)  
            f(2,i) = 2.0_double*hgh_nonlocal_stress_radial_f_value_i(pd%o,g(i),l,c)
           end do
        else
        do i = 1,size(g)
          f(:,i) = nonlocal_stress_radial_f_valu_i(pd%o,g(i),l)
        end do
        end if
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::projector_stress_f_values_pd")) continue
      end function

      function projector_r_value_pd(pd,r,l,gi_cut,go_cut) result(v)
!doc$ function projector_r_value(pd,r,l,gi_cut,go_cut) result(v)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: r
        integer, intent(in) :: l
        real(double), intent(in) :: gi_cut, go_cut
        real(double) :: v
!       requires: x_nlcomps(pd) be > 0. l be one of pd%o%l(i). gi_cut < go_cut.
!       effects: Returns the value of the l'th component of the optimized non-local pseudopotential at r.
!                Initializes pd%o%nlpo(:,:) if this has not already been done.
!       errors: r > pd%o%r_opt. Passes errors

!doc$
        logical :: gi_cut_change, go_cut_change
        call my(pd)
        if (error(r > pd%o%r_opt,"ERROR: r > r_opt")) goto 100
        if (.not.associated( pd%o%nlpo )) then
          pd%o%gi_cut = gi_cut
          pd%o%go_cut = go_cut
          call optimize_nlp_i(pd%o) ; if (error()) goto 100
        else
          gi_cut_change = ( gi_cut /= pd%o%gi_cut )
          go_cut_change = ( go_cut /= pd%o%go_cut )
          if (gi_cut_change .or. go_cut_change) then
            if (gi_cut_change) pd%o%gi_cut = gi_cut
            if (go_cut_change) pd%o%go_cut = go_cut
            if (associated( pd%o%nlpo )) deallocate( pd%o%nlpo )
            if (associated( pd%o%g )) deallocate( pd%o%g )
            if (associated( pd%o%w_max )) deallocate( pd%o%w_max )
            call warn("WARNING: real-space projectors are being re-optimized")
            call optimize_nlp_i(pd%o) ; if (error()) goto 100
          end if
        end if
        v = nonlocal_radial_r_value_i(pd%o,r,l) ; if (error()) goto 100
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::projector_r_value_pd")) continue
      end function

      function projector_r_gradients_pd(pd,r,l) result(g)
!doc$ function projector_r_gradients(pd,r,l) result(g)
        type(ncp_data_obj) :: pd
        real(double), intent(in) :: r
        integer, intent(in) :: l
        real(double), dimension(2) :: g
!       requires: x_nlcomps(pd) be > 0. l be one of pd%o%l(i).
!       effects: Returns the gradients of the l'th component of the optimized non-local pseudopotential at r.
!       errors: r > pd%o%r_opt. pd%o%nlpo not associated. Passes errors

!doc$
        call my(pd)
        if (error(r > pd%o%r_opt,"ERROR: r > r_opt")) goto 100
        if (error(.not.associated( pd%o%nlpo ),"ERROR: real-space projectors not yet optimized")) goto 100
        g = nonlocal_radial_r_gradients_i(pd%o,r,l) ; if (error()) goto 100
100     call glean(thy(pd))
        if (error("Exit ncp_data_mod::projector_r_gradients_pd")) continue
      end function

! local routines

      subroutine read_hgh_i(pdr,f)
        type(ncp_data_rep) :: pdr
        type(file_obj) :: f

        integer :: nloc, iloc, inonloc, il, jl, al, void
        integer :: ios
        if (i_access(f)) read(x_unit(f),*, IOSTAT=ios) pdr%zatom, pdr%valence
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: problem opening file = "//x_name(f))) goto 100
        if (i_comm(f)) then
           call broadcast(FILE_SCOPE,pdr%zatom)
           call broadcast(FILE_SCOPE,pdr%valence)
        end if
        if (i_access(f)) read(x_unit(f),*) pdr%pspcod, pdr%pspxc, pdr%lmax
        if (i_comm(f)) then
           call broadcast(FILE_SCOPE,pdr%pspcod)
           call broadcast(FILE_SCOPE,pdr%pspxc)
           call broadcast(FILE_SCOPE,pdr%lmax)
        end if
        if (i_access(f)) read(x_unit(f),'(F15.8, I5)', ADVANCE='NO', IOSTAT=ios) pdr%rloc, nloc
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: reading non advancing input = "//x_name(f))) goto 100 
        if (i_comm(f)) then
           call broadcast(FILE_SCOPE,pdr%rloc)
           call broadcast(FILE_SCOPE,nloc)
        end if
        allocate ( pdr%cn(nloc) )
        if (i_access(f)) read(x_unit(f),*) pdr%cn
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%cn)
        if (i_access(f)) read(x_unit(f),*) pdr%nnonloc
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%nnonloc)
        allocate(pdr%nl(pdr%nnonloc))
        do inonloc = 1,pdr%nnonloc      
           if (i_access(f)) read(x_unit(f),'(F15.8, I5)', ADVANCE='NO', IOSTAT=ios) pdr%nl(inonloc)%r, pdr%nl(inonloc)%c
           if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
           if (error(ios /= 0,"ERROR: reading non advancing input of r and nl = "//x_name(f))) goto 100
              if (i_comm(f)) then
                call broadcast(FILE_SCOPE,pdr%nl(inonloc)%r)
                call broadcast(FILE_SCOPE,pdr%nl(inonloc)%c)
              end if
           if (inonloc == 1) then
              allocate(pdr%nl(inonloc)%h(pdr%nl(inonloc)%c,pdr%nl(inonloc)%c))
              nullify(pdr%nl(inonloc)%k)
!read in h parameter matrix for s orbitals
              do il = 1,pdr%nl(inonloc)%c
                 do jl = il,pdr%nl(inonloc)%c
                    if (jl /= pdr%nl(inonloc)%c) then
       !advance through blank space in file                
                       if (il /= 1) then
                           do al = 1,il
                            if (al==1) then
                            if (i_access(f)) then
                            read(x_unit(f),'(I20)', ADVANCE='no') void
                            end if
                            else
                            if (i_access(f)) then
                            read(x_unit(f), '(I15)', ADVANCE='no') void
                            end if
                            end if
                           end do
                       end if
      !resume normal read in
                      if (i_access(f)) then 
                      read(x_unit(f),'(F15.8)', ADVANCE='no') pdr%nl(inonloc)%h(il,jl)
                      end if
                    else
                       if (i_access(f)) then
                       read(x_unit(f),*) pdr%nl(inonloc)%h(il, jl)
                       end if
                    end if
!fill in symmetric elements
                    if (jl /= il) pdr%nl(inonloc)%h(jl,il)=pdr%nl(inonloc)%h(il,jl)
                  end do     
              end do
              if (i_comm(f)) then
                call broadcast(FILE_SCOPE,pdr%nl(inonloc)%h)
              end if
           else
            allocate(pdr%nl(inonloc)%h(pdr%nl(inonloc)%c,pdr%nl(inonloc)%c))
            allocate(pdr%nl(inonloc)%k(pdr%nl(inonloc)%c,pdr%nl(inonloc)%c))
!read in h parameter matrix for higher non-s states
              do il = 1,pdr%nl(inonloc)%c
                 do jl = il,pdr%nl(inonloc)%c
                    if (jl /= pdr%nl(inonloc)%c) then
                     !advance through blank space in document
                       if (il /= 1) then
                           do al = 1,il
                            if (al==1) then
                            if (i_access(f)) then
                            read(x_unit(f),'(I20)', ADVANCE='no') void
                            end if
                            else
                            if (i_access(f)) then
                            read(x_unit(f), '(I15)', ADVANCE='no') void
                            end if
                            end if
                           end do
                        end if
!resume read in
                       if (i_access(f)) then
                      read(x_unit(f),'(F15.8)', ADVANCE='no') pdr%nl(inonloc)%h(il,jl)
                      end if
                    else
                       if (i_access(f)) then
                       read(x_unit(f),*) pdr%nl(inonloc)%h(il, jl)
                       end if
                    end if
!fill in symmetric elements
                    if (jl /= il) pdr%nl(inonloc)%h(jl,il)=pdr%nl(inonloc)%h(il,jl)
                  end do
              end do
              if (i_comm(f)) then
                 call broadcast(FILE_SCOPE,pdr%nl(inonloc)%h)
              end if
!read in  k parameter matrix for non-s states
              do il = 1,pdr%nl(inonloc)%c
                 do jl = il,pdr%nl(inonloc)%c
                    if (jl /= pdr%nl(inonloc)%c) then
        !advance through blank space in file
                      do al = 1,il
                            if (al==1) then
                            if (i_access(f)) then
                            read(x_unit(f),'(I20)', ADVANCE='no') void
                            end if
                            else
                            if (i_access(f)) then
                            read(x_unit(f), '(I15)', ADVANCE='no') void
                            end if
                            end if
                      end do
       !resume normal read in
                      if (i_access(f)) then
                      read(x_unit(f),'(F15.8)', ADVANCE='no') pdr%nl(inonloc)%k(il,jl)
                      end if
                    else
                       if (i_access(f)) then
                       read(x_unit(f),*) pdr%nl(inonloc)%k(il, jl)
                       end if
                    end if
!fill in symmetric elements
                    if (jl /= il) pdr%nl(inonloc)%k(jl,il)=pdr%nl(inonloc)%k(il,jl)
                  end do
              end do
              if (i_comm(f)) then
                 call broadcast(FILE_SCOPE,pdr%nl(inonloc)%k)
              end if
              end if
     end do
100      end subroutine


      subroutine read_simple_i(pdr,f)                         !------------(begin block)
        type(ncp_data_rep) :: pdr                             ! simple
        type(file_obj) :: f                                   ! <exponent>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent       !------------(end block)
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent)
        nullify( pdr%l )
        nullify( pdr%kb )
        nullify( pdr%r )
        nullify( pdr%lp )
        nullify( pdr%nlp )
        nullify( pdr%vd )
        nullify( pdr%cd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
100     if (error("Exit ncp_data_mod::read_simple_i")) continue
      end subroutine

      subroutine read_linear_i(pdr,f)                          !--------------------------------------------(begin block)
        type(ncp_data_rep) :: pdr                              ! linear
        type(file_obj) :: f                                    ! <exponent>
!       requires: Data be on an equal-spaced radial mesh.      ! <numl>
        integer :: il, ir, nl, nr                              ! <l(1)> <l(2)> .... <l(numl)>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent        ! <kbparm(1)> <kbparm(2)> .... <kbparm(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent) ! <numrad>
        if (i_access(f)) read(x_unit(f),*) nl                  ! <gridval> <loc> <nloc(1)> ... <nloc(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,nl)           ! <gridval> <loc> <nloc(1)> ... <nloc(numl)>
        if (nl == 0) then                                      ! ..........................................
          nullify( pdr%l )                                     ! <gridval> <loc> <nloc(1)> ... <nloc(numl)>
          nullify( pdr%kb )                                    !--------------------------------------------(end block)
        else
          allocate( pdr%l(nl) )
          allocate( pdr%kb(nl) )
          if (i_access(f)) then
            read(x_unit(f),*) pdr%l
            read(x_unit(f),*) pdr%kb
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%l)
            call broadcast(FILE_SCOPE,pdr%kb)
          end if
        end if
        if (i_access(f)) read(x_unit(f),*) nr
        if (i_comm(f)) call broadcast(FILE_SCOPE,nr)
        if (nl == 0) then
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          nullify( pdr%nlp )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
          end if
        else
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          allocate( pdr%nlp(nr,nl) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir), (pdr%nlp(ir,il),il=1,nl)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%nlp)
          end if
        end if
        nullify( pdr%vd )
        nullify( pdr%cd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
        call canonicalize_linear_i(pdr) ; if (error()) goto 100
100     if (error("Exit ncp_data_mod::read_linear_i")) continue
      end subroutine

      subroutine read_linear_cc_i(pdr,f)                       !--------------------------------------------------(begin block)
        type(ncp_data_rep) :: pdr                              ! linear_cc
        type(file_obj) :: f                                    ! <exponent>
!       requires: Data be on an equal-spaced radial mesh.      ! <numl>
        integer :: il, ir, nl, nr                              ! <l(1)> <l(2)> .... <l(numl)>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent        ! <kbparm(1)> <kbparm(2)> .... <kbparm(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent) ! <numrad>
        if (i_access(f)) read(x_unit(f),*) nl                  ! <gridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
        if (i_comm(f)) call broadcast(FILE_SCOPE,nl)           ! <gridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
        if (nl == 0) then                                      ! .................................................
          nullify( pdr%l )                                     ! <gridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
          nullify( pdr%kb )                                    !--------------------------------------------------(end block)
        else
          allocate( pdr%l(nl) )
          allocate( pdr%kb(nl) )
          if (i_access(f)) then
            read(x_unit(f),*) pdr%l
            read(x_unit(f),*) pdr%kb
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%l)
            call broadcast(FILE_SCOPE,pdr%kb)
          end if
        end if
        if (i_access(f)) read(x_unit(f),*) nr
        if (i_comm(f)) call broadcast(FILE_SCOPE,nr)
        if (nl == 0) then
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          nullify( pdr%nlp )
          allocate( pdr%cd(nr) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir), pdr%cd(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%cd)
          end if
        else
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          allocate( pdr%nlp(nr,nl) )
          allocate( pdr%cd(nr) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir), (pdr%nlp(ir,il),il=1,nl), pdr%cd(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%nlp)
            call broadcast(FILE_SCOPE,pdr%cd)
          end if
        end if
        nullify( pdr%vd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
        call canonicalize_linear_i(pdr) ; if (error()) goto 100
100     if (error("Exit ncp_data_mod::read_linear_cc_i")) continue
      end subroutine

      subroutine canonicalize_linear_i(pdr)
        type(ncp_data_rep) :: pdr
!       requires: Data be on an equal-spaced radial mesh.

        integer :: il, nl, nr
        real(double) :: dr
        real(double), dimension(:), allocatable :: wt
        real(double), dimension(:,:), allocatable :: f

        nr = size(pdr%r)
        dr = pdr%r(2) - pdr%r(1)
        allocate( wt(nr) )
        wt(1:nr:8) =   989.0_double
        wt(2:nr:8) =  5888.0_double
        wt(3:nr:8) =  -928.0_double
        wt(4:nr:8) = 10496.0_double
        wt(5:nr:8) = -4540.0_double
        wt(6:nr:8) = 10496.0_double
        wt(7:nr:8) =  -928.0_double
        wt(8:nr:8) =  5888.0_double
        if (mod(size(wt)-1,8) == 0) then
          wt(9:size(wt)-1:8) = wt(9:size(wt)-1:8)*2.0_double
        else
          wt(9:size(wt):8) = wt(9:size(wt):8)*2.0_double
        end if
        wt = 4.0_double*wt/14175.0_double
        pdr%lp = eight_pi*pdr%r*pdr%lp*wt*dr
        if (associated( pdr%vd )) pdr%vd = pdr%r*pdr%vd*wt*dr
        if (associated( pdr%cd )) pdr%cd = four_pi*pdr%r*pdr%cd*wt*dr
        if (associated( pdr%l )) then
          nl = size(pdr%l)
          allocate( f(nr,nl) )
          do il = 1,nl
            f(1,il) = sqrt(2.0_double)*pdr%nlp(1,il)
            f(2:nr,il) = sqrt(2.0_double)*pdr%nlp(2:nr,il)/pdr%r(2:nr)**2
          end do
          pdr%r_cut = pdr%r(cutoff_point_i(f)) ; if (error()) goto 100
          do il = 1,nl
            pdr%nlp(:,il) = sqrt(2.0_double)*four_pi*pdr%nlp(:,il)*wt*dr
          end do
        end if

100     if (allocated( f )) deallocate( f )
        if (allocated( wt )) deallocate( wt )

        if (error("Exit ncp_data_mod::canonicalize_linear_i")) continue

      end subroutine

      subroutine read_log_i(pdr,f)                             !------------------------------------------------------(begin block)
        type(ncp_data_rep) :: pdr                              ! log
        type(file_obj) :: f                                    ! <exponent>
        integer :: il, ir, nl, nr                              ! <numl>
        real(double), dimension(:), allocatable :: dr          ! <l(1)> <l(2)> .... <l(numl)>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent        ! <kbparm(1)> <kbparm(2)> .... <kbparm(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent) ! <numrad>
        if (i_access(f)) read(x_unit(f),*) nl                  ! <gridval> <dgridval> <loc> <nloc(1)> ... <nloc(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,nl)           ! <gridval> <dgridval> <loc> <nloc(1)> ... <nloc(numl)>
        if (nl == 0) then                                      ! .....................................................
          nullify( pdr%l )                                     ! <gridval> <dgridval> <loc> <nloc(1)> ... <nloc(numl)>
          nullify( pdr%kb )                                    !------------------------------------------------------(end block)
        else
          allocate( pdr%l(nl) )
          allocate( pdr%kb(nl) )
          if (i_access(f)) then
            read(x_unit(f),*) pdr%l
            read(x_unit(f),*) pdr%kb
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%l)
            call broadcast(FILE_SCOPE,pdr%kb)
          end if
        end if
        if (i_access(f)) read(x_unit(f),*) nr
        if (i_comm(f)) call broadcast(FILE_SCOPE,nr)
        if (nl == 0) then
          allocate( pdr%r(nr) )
          allocate( dr(nr) )
          allocate( pdr%lp(nr) )
          nullify( pdr%nlp )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), dr(ir), pdr%lp(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,dr)
            call broadcast(FILE_SCOPE,pdr%lp)
          end if
        else
          allocate( pdr%r(nr) )
          allocate( dr(nr) )
          allocate( pdr%lp(nr) )
          allocate( pdr%nlp(nr,nl) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), dr(ir), pdr%lp(ir), (pdr%nlp(ir,il),il=1,nl)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,dr)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%nlp)
          end if
        end if
        nullify( pdr%vd )
        nullify( pdr%cd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
        call canonicalize_log_i(pdr,dr) ; if (error()) goto 100
100     if (allocated( dr )) deallocate( dr )
        if (error("Exit ncp_data_mod::read_log_i")) continue
      end subroutine

      subroutine read_log_cc_i(pdr,f)                          !-------------------------------------------------------(begin block)
        type(ncp_data_rep) :: pdr                              ! log_cc
        type(file_obj) :: f                                    ! <exponent>
        integer :: il, ir, nl, nr                              ! <numl>
        real(double), dimension(:), allocatable :: dr          ! <l(1)> <l(2)> .... <l(numl)>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent        ! <kbparm(1)> <kbparm(2)> .... <kbparm(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent) ! <numrad>
        if (i_access(f)) read(x_unit(f),*) nl                  ! <gridval> <dgridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
        if (i_comm(f)) call broadcast(FILE_SCOPE,nl)           ! <gridval> <dgridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
        if (nl == 0) then                                      ! ............................................................
          nullify( pdr%l )                                     ! <gridval> <dgridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
          nullify( pdr%kb )                                    !---------------------------------------------------------(end block)
        else
          allocate( pdr%l(nl) )
          allocate( pdr%kb(nl) )
          if (i_access(f)) then
            read(x_unit(f),*) pdr%l
            read(x_unit(f),*) pdr%kb
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%l)
            call broadcast(FILE_SCOPE,pdr%kb)
          end if
        end if
        if (i_access(f)) read(x_unit(f),*) nr
        if (i_comm(f)) call broadcast(FILE_SCOPE,nr)
        if (nl == 0) then
          allocate( pdr%r(nr) )
          allocate( dr(nr) )
          allocate( pdr%lp(nr) )
          nullify( pdr%nlp )
          allocate( pdr%cd(nr) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), dr(ir), pdr%lp(ir), pdr%cd(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,dr)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%cd)
          end if
        else
          allocate( pdr%r(nr) )
          allocate( dr(nr) )
          allocate( pdr%lp(nr) )
          allocate( pdr%nlp(nr,nl) )
          allocate( pdr%cd(nr) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), dr(ir), pdr%lp(ir), (pdr%nlp(ir,il),il=1,nl), pdr%cd(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,dr)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%nlp)
            call broadcast(FILE_SCOPE,pdr%cd)
          end if
        end if
        nullify( pdr%vd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
        call canonicalize_log_i(pdr,dr) ; if (error()) goto 100
100     if (allocated( dr )) deallocate( dr )
        if (error("Exit ncp_data_mod::read_log_cc_i")) continue
      end subroutine

      subroutine canonicalize_log_i(pdr,dr)
        type(ncp_data_rep) :: pdr
        real(double), dimension(:), intent(in) :: dr

        integer :: il, ir, ll, nl, nr, ns
        real(double), dimension(:), allocatable :: wt
        real(double), dimension(:,:), allocatable :: f

        nr = size(pdr%r)
        allocate( wt(nr) )
        if (mod(nr,2) == 0) then
          ns = 4
          wt(1) = 3.0_double/8.0_double
          wt(2) = 9.0_double/8.0_double
          wt(3) = 9.0_double/8.0_double
          wt(4) = 3.0_double/8.0_double + 1.0_double/3.0_double
        else
          ns = 1
          wt(1) = 1.0_double/3.0_double
        end if
        ll = 4
        do ir = ns+1,nr-1
          wt(ir) = real(ll,double)/3.0_double
          ll = 6 - ll
        end do
        wt(nr) = 1.0_double/3.0_double
        pdr%lp = four_pi*pdr%r*pdr%lp*wt*dr
        if (associated( pdr%vd )) pdr%vd = pdr%r*pdr%vd*wt*dr
        if (associated( pdr%cd )) pdr%cd = four_pi*pdr%r*pdr%cd*wt*dr
        if (associated( pdr%l )) then
          nl = size(pdr%l)
          allocate( f(nr,nl) )
          f = pdr%nlp
          pdr%r_cut = pdr%r(cutoff_point_i(f)) ; if (error()) goto 100
          do il = 1,nl
            pdr%nlp(:,il) = four_pi*pdr%r**2*pdr%nlp(:,il)*wt*dr
          end do
        end if

100     if (allocated( f )) deallocate( f )
        if (allocated( wt )) deallocate( wt )

        if (error("Exit ncp_data_mod::canonicalize_log_i")) continue

      end subroutine

      subroutine read_log_wt_i(pdr,f)                          !-------------------------------------------(begin block) 
        type(ncp_data_rep) :: pdr                              ! log_wt
        type(file_obj) :: f                                    ! <exponent>
        integer :: il, ir, nl, nr                              ! <numl>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent        ! <l(1)> <l(2)> .... <l(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent) ! <kbparm(1)> <kbparm(2)> .... <kbparm(numl)>
        if (i_access(f)) read(x_unit(f),*) nl                  ! <numrad>
        if (i_comm(f)) call broadcast(FILE_SCOPE,nl)           ! <gridval> <loc> <nloc(1)> ... <nloc(numl)>
        if (nl == 0) then                                      ! <gridval> <loc> <nloc(1)> ... <nloc(numl)>
          nullify( pdr%l )                                     ! ..........................................
          nullify( pdr%kb )                                    ! <gridval> <loc> <nloc(1)> ... <nloc(numl)>
        else                                                   !-------------------------------------------(end block)
          allocate( pdr%l(nl) )
          allocate( pdr%kb(nl) )
          if (i_access(f)) then
            read(x_unit(f),*) pdr%l
            read(x_unit(f),*) pdr%kb
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%l)
            call broadcast(FILE_SCOPE,pdr%kb)
          end if
        end if
        if (i_access(f)) read(x_unit(f),*) nr
        if (i_comm(f)) call broadcast(FILE_SCOPE,nr)
        if (nl == 0) then
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          nullify( pdr%nlp )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
          end if
        else
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          allocate( pdr%nlp(nr,nl) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir), (pdr%nlp(ir,il),il=1,nl)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%nlp)
          end if
        end if
        nullify( pdr%vd )
        nullify( pdr%cd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
        call canonicalize_log_wt_i(pdr) ; if (error()) goto 100
100     if (error("Exit ncp_data_mod::read_log_wt_i")) continue
      end subroutine

      subroutine read_log_wt_cc_i(pdr,f)                       !--------------------------------------------------(begin block) 
        type(ncp_data_rep) :: pdr                              ! log_wt_cc
        type(file_obj) :: f                                    ! <exponent>
        integer :: il, ir, nl, nr                              ! <numl>
        if (i_access(f)) read(x_unit(f),*) pdr%exponent        ! <l(1)> <l(2)> .... <l(numl)>
        if (i_comm(f)) call broadcast(FILE_SCOPE,pdr%exponent) ! <kbparm(1)> <kbparm(2)> .... <kbparm(numl)>
        if (i_access(f)) read(x_unit(f),*) nl                  ! <numrad>
        if (i_comm(f)) call broadcast(FILE_SCOPE,nl)           ! <gridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
        if (nl == 0) then                                      ! <gridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
          nullify( pdr%l )                                     ! .................................................
          nullify( pdr%kb )                                    ! <gridval> <loc> <nloc(1)> ... <nloc(numl)> <core>
        else                                                   !--------------------------------------------------(end block)
          allocate( pdr%l(nl) )
          allocate( pdr%kb(nl) )
          if (i_access(f)) then
            read(x_unit(f),*) pdr%l
            read(x_unit(f),*) pdr%kb
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%l)
            call broadcast(FILE_SCOPE,pdr%kb)
          end if
        end if
        if (i_access(f)) read(x_unit(f),*) nr
        if (i_comm(f)) call broadcast(FILE_SCOPE,nr)
        if (nl == 0) then
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          nullify( pdr%nlp )
          allocate( pdr%cd(nr) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir), pdr%cd(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%cd)
          end if
        else
          allocate( pdr%r(nr) )
          allocate( pdr%lp(nr) )
          allocate( pdr%nlp(nr,nl) )
          allocate( pdr%cd(nr) )
          if (i_access(f)) then
            do ir = 1,nr
              read(x_unit(f),*) pdr%r(ir), pdr%lp(ir), (pdr%nlp(ir,il),il=1,nl), pdr%cd(ir)
            end do
          end if
          if (i_comm(f)) then
            call broadcast(FILE_SCOPE,pdr%r)
            call broadcast(FILE_SCOPE,pdr%lp)
            call broadcast(FILE_SCOPE,pdr%nlp)
            call broadcast(FILE_SCOPE,pdr%cd)
          end if
        end if
        nullify( pdr%vd )
        nullify( pdr%w_max )
        nullify( pdr%g )
        nullify( pdr%nlpo )
        call canonicalize_log_wt_i(pdr) ; if (error()) goto 100
100     if (error("Exit ncp_data_mod::read_log_wt_cc_i")) continue
      end subroutine

      subroutine canonicalize_log_wt_i(pdr)
        type(ncp_data_rep) :: pdr

        integer :: il, ir, ll, nl, nr, ns
        real(double) :: b
        real(double), dimension(:), allocatable :: dr, wt
        real(double), dimension(:,:), allocatable :: f

        if (associated( pdr%l )) then
          nr = size(pdr%r)
          allocate( dr(nr) )
          b = log(pdr%r(nr)/pdr%r(nr-1))
          dr = b*pdr%r
          allocate( wt(nr) )
          if (mod(nr,2) == 0) then
            ns = 4
            wt(1) = 3.0_double/8.0_double
            wt(2) = 9.0_double/8.0_double
            wt(3) = 9.0_double/8.0_double
            wt(4) = 3.0_double/8.0_double + 1.0_double/3.0_double
          else
            ns = 1
            wt(1) = 1.0_double/3.0_double
          end if
          ll = 4
          do ir = ns+1,nr-1
            wt(ir) = real(ll,double)/3.0_double
            ll = 6 - ll
          end do
          wt(nr) = 1.0_double/3.0_double
          nl = size(pdr%l)
          allocate( f(nr,nl) )
          do il = 1,nl
            f(1,il) = pdr%nlp(1,il)
            f(2:nr,il) = pdr%nlp(2:nr,il)/(four_pi*pdr%r(2:nr)**2*wt(2:nr)*dr(2:nr))
          end do
          pdr%r_cut = pdr%r(cutoff_point_i(f)) ; if (error()) goto 100
        end if

100     if (allocated( dr )) deallocate( dr )
        if (allocated( f )) deallocate( f )
        if (allocated( wt )) deallocate( wt )

        if (error("Exit ncp_data_mod::canonicalize_log_wt_i")) continue

      end subroutine

      function radial_f_value_i(pdr,pdc,g,l) result(f)
        type(ncp_data_rep) :: pdr
        real(double), dimension(:), intent(in) :: pdc
        real(double), intent(in) :: g
        integer, intent(in) :: l
        real(double) :: f
        f = sum(pdc*spherical_bessel(g*pdr%r,l))
        if (error("Exit ncp_data_mod::radial_f_value_i")) continue
      end function

      function nonlocal_radial_f_value_i(pdr,g,l) result(f)
        type(ncp_data_rep) :: pdr
        real(double), intent(in) :: g
        integer, intent(in) :: l
        real(double) :: f

        integer :: i, il
        do i = 1,size(pdr%l)
          if (pdr%l(i) == l) then
            il = i
            exit
          end if
        end do
        select case(l)
        case(0)
          f = sum(pdr%nlp(:,il)*spherical_bessel(g*pdr%r,l))
        case default
          f = sum(pdr%nlp(:,il)*pdr%r**l*spherical_bessel(g*pdr%r,l,.true.))
        end select
      end function
      
      function hgh_nonlocal_radial_f_value_i(pdr,g,l,c) result(f)
        type(ncp_data_rep) :: pdr
        real(double), intent(in) :: g
        integer, intent(in) :: l, c
        real(double) :: f


        logical :: found
        real(double) :: sf
        select case(l)
           case(0)
              select case(c)
                 case(1)
                    f = 4*sqrt(2.0_double)*sqrt(pdr%nl(1)%r**3)*pi**(1.25000)/exp((g*pdr%nl(1)%r)**2/2)
                 case(2)
                    f = 8*sqrt(2*(pdr%nl(1)%r**3))*sqrt(0.06666667)*pi**(1.250000)*(3-(g*pdr%nl(1)%r)**2)/exp((g*pdr%nl(1)%r)**2/2)
                 case(3)
                    f = 16*sqrt(2*(pdr%nl(1)%r**3)/105)*pi**(1.250000)*(15-10*g**2*pdr%nl(1)%r**2+g**4*pdr%nl(1)%r**4)/(3*exp((g*pdr%nl(1)%r)**2/2))
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select      
           case(1)
              select case(c)
                 case(1)
                    f = 8*sqrt((pdr%nl(2)%r**5))*sqrt(0.3333333333)*pi**(1.250000)/exp((g*pdr%nl(2)%r)**2/2)
                 case(2)
                    f = 16*sqrt((pdr%nl(2)%r**5)/105)*pi**(1.250000)*(5-g**2*pdr%nl(2)%r**2)/exp((g*pdr%nl(2)%r)**2/2)
                 case(3)
                    f = 32*sqrt((pdr%nl(2)%r**5)/1155)*pi**(1.250000)*(35-14*g**2*pdr%nl(2)%r**2+g**4*pdr%nl(2)%r**4)/(3*exp((g*pdr%nl(2)%r)**2/2))
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select
           case(2)
              select case(c)
                 case(1)
                    f = 8*sqrt((2*pdr%nl(3)%r**7)/15)*pi**(1.250000)/exp((g*pdr%nl(3)%r)**2/2)
                 case(2)
                    f = 16*sqrt((2*pdr%nl(3)%r**7)/105)*pi**(1.250000)*(7-g**2*pdr%nl(3)%r**2)/(3*exp((g*pdr%nl(3)%r)**2/2))
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select
           case(3)
              select case(c)
                 case(1)
                    f = 16*sqrt((pdr%nl(4)%r**9)/105)*pi**(1.250000)/exp((g*pdr%nl(4)%r)**2/2)
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select      
        end select      
100  end function  

      function nonlocal_stress_radial_f_valu_i(pdr,g,l) result(f)
        type(ncp_data_rep) :: pdr
        real(double), intent(in) :: g
        integer, intent(in) :: l
        real(double), dimension(2) :: f

        integer :: i, il
        do i = 1,size(pdr%l)
          if (pdr%l(i) == l) then
            il = i
            exit
          end if
        end do
        select case(l)
        case(0)
          f(1) = sum(pdr%nlp(:,il)*spherical_bessel(g*pdr%r,l))
        case default
          f(1) = sum(pdr%nlp(:,il)*pdr%r**l*spherical_bessel(g*pdr%r,l,.true.))
        end select
        f(2) = sum(pdr%nlp(:,il)*pdr%r**(l+2)*spherical_bessel(g*pdr%r,l+1,.true.))
        if (error("Exit ncp_data_mod::nonlocal_stress_radial_f_valu_i")) continue
      end function

      function hgh_nonlocal_stress_radial_f_value_i(pdr,g,l,c) result(f)
        type(ncp_data_rep) :: pdr
        real(double), intent(in) :: g
        integer, intent(in) :: l, c
        real(double) :: f


        logical :: found
        real(double) :: sf
        select case(l)
           case(0)
              select case(c)
                 case(1)
                    f = 4*sqrt(2.0_double)*sqrt(pdr%nl(1)%r**7)*pi**(1.25000)/exp((g*pdr%nl(1)%r)**2/2)
                 case(2)
                    f = 8*sqrt(2*(pdr%nl(1)%r**7))*sqrt(0.06666667)*pi**(1.250000)*(5-(g*pdr%nl(1)%r)**2)/exp((g*pdr%nl(1)%r)**2/2)
                 case(3)
                    f = 16*sqrt(2*(pdr%nl(1)%r**7)/105)*pi**(1.250000)*(35-14*g**2*pdr%nl(1)%r**2+g**4*pdr%nl(1)%r**4)/(3*exp((g*pdr%nl(1)%r)**2/2))
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select
           case(1)
              select case(c)
                 case(1)
                    f = 8*sqrt((pdr%nl(2)%r**9))*sqrt(0.3333333333)*pi**(1.250000)/exp((g*pdr%nl(2)%r)**2/2)
                 case(2)
                    f = 16*sqrt((pdr%nl(2)%r**9)/105)*pi**(1.250000)*(7-g**2*pdr%nl(2)%r**2)/exp((g*pdr%nl(2)%r)**2/2)
                 case(3)
                    f = 32*sqrt((pdr%nl(2)%r**9)/1155)*pi**(1.250000)*(63-18*g**2*pdr%nl(2)%r**2+g**4*pdr%nl(2)%r**4)/(3*exp((g*pdr%nl(2)%r)**2/2))
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select
           case(2)
              select case(c)
                 case(1)
                    f = 8*sqrt((2*pdr%nl(3)%r**11)/15)*pi**(1.250000)/exp((g*pdr%nl(3)%r)**2/2)
                 case(2)
                    f = 16*sqrt((2*pdr%nl(3)%r**11)/105)*pi**(1.250000)*(9-g**2*pdr%nl(3)%r**2)/(3*exp((g*pdr%nl(3)%r)**2/2))
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select
           case(3)
              select case(c)
                 case(1)
                    f = 16*sqrt((pdr%nl(4)%r**13)/105)*pi**(1.250000)/exp((g*pdr%nl(4)%r)**2/2)
                 case default
                    if (error(.true., "ERROR: Invalid angular momentum channel")) goto 100
              end select
        end select
100  end function


      function nonlocal_radial_r_value_i(pdr,r,l) result(v)
        type(ncp_data_rep) :: pdr
        real(double), intent(in) :: r
        integer, intent(in) :: l
        real(double) :: v

        integer :: i, il
        do i = 1,size(pdr%l)
          if (pdr%l(i) == l) then
            il = i
            exit
          end if
        end do
        select case(l)
        case(0)
          v = sum(pdr%nlpo(:,il)*spherical_bessel(pdr%g*r,l,.true.)) ; if (error()) goto 100
        case default
          v = sum(pdr%nlpo(:,il)*pdr%g**(l)*spherical_bessel(pdr%g*r,l,.true.)) ; if (error()) goto 100
        end select
100     if (error("Exit ncp_data_mod::nonlocal_radial_r_value_i")) continue
      end function

      function nonlocal_radial_r_gradients_i(pdr,r,l) result(g)
        type(ncp_data_rep) :: pdr
        real(double), intent(in) :: r
        integer, intent(in) :: l
        real(double), dimension(2) :: g

        integer :: i, il
        do i = 1,size(pdr%l)
          if (pdr%l(i) == l) then
            il = i
            exit
          end if
        end do
        select case(l)
        case(0)
          g(1) = 0.0_double
        case default
          g(1) = sum(pdr%nlpo(:,il)*pdr%g**(l)*spherical_bessel(pdr%g*r,l,.true.)) ; if (error()) goto 100
        end select
        g(2) = sum(pdr%nlpo(:,il)*pdr%g**(l+2)*spherical_bessel(pdr%g*r,l+1,.true.)) ; if (error()) goto 100
100     if (error("Exit ncp_data_mod::nonlocal_radial_r_gradients_i")) continue
      end function

      function radial_r_value_i(pdr,pdc,r,l) result(f)
        type(ncp_data_rep) :: pdr
        real(double), dimension(:), intent(in) :: pdc
        real(double), intent(in) :: r
        integer, intent(in) :: l
        real(double) :: f
        f = sum(pdc*spherical_bessel(pdr%g*r,l))
        if (error("Exit ncp_data_mod::radial_r_value_i")) continue
      end function

      subroutine optimize_nlp_i(pdr)
        type(ncp_data_rep) :: pdr

        integer, parameter :: nr = 1000

        logical :: found
        integer :: ib, ic, ig, ir, l, ll, nb, ng_d, ng_i, ng_o, ns
        real(double) :: dg, dr, r0, r1
        real(double), dimension(:), allocatable :: b, b0, b1, f, gr, r, wg, wt
        real(double), dimension(:,:), allocatable :: a, t

        call arg("optimization_points",ng_o,found)
        if (.not.found) ng_o = 181

        allocate( pdr%g(ng_o) )  ! GENERATE A RADIAL G-VECTOR_MESH.
        dg = pdr%go_cut/real(ng_o,double)
        do ig = 1,ng_o
          pdr%g(ig) = dg*real(ig-1,double)
        end do

        ng_i = int(pdr%gi_cut/dg)
        ng_d = ng_o - ng_i

        nb = size(pdr%l)

        allocate( pdr%nlpo(ng_o,nb) )

        allocate( f(size(pdr%r)) )  ! COMPUTE FOURIER COMPONENTS UP TO GI_CUT.
        do ib = 1,nb
          f = pdr%nlp(:,ib)/four_pi
          do ig = 1,ng_i
            pdr%nlpo(ig,ib) = radial_f_value_i(pdr,f,pdr%g(ig),pdr%l(ib)) ; if (error()) goto 100
          end do
        end do
        deallocate( f )

        allocate( a(ng_d,ng_d), b(ng_d), t(ng_o,ng_o) )  ! COMPUTE FOURIER COMPONENTS BETWEEN GI_CUT AND GO_CUT.
        allocate( gr(ng_o), b0(ng_o), b1(ng_o) )
        gr = pdr%g*pdr%r_opt
        r0 = 1.0_double/pdr%r_opt
        do ib = 1,nb
          l = pdr%l(ib)
          r1 = real(2*l+1,double)
          b0 = spherical_bessel(gr,l,.true.) ; if (error()) goto 100
          b1 = spherical_bessel(gr,l+1,.true.) ; if (error()) goto 100
          do ic = 1,ng_o
            do ir = 1,ic-1
              t(ir,ic) = r0*(gr(ir)*gr(ic))**(l+2)/(gr(ir)**2 - gr(ic)**2)*(gr(ir)**2*b1(ir)*b0(ic) - gr(ic)**2*b0(ir)*b1(ic))
              t(ic,ir) = t(ir,ic)
            end do
            t(ic,ic) = 0.5_double*r0*gr(ic)**(2*l+4)*(b0(ic)**2 + (gr(ic)*b1(ic))**2 - r1*b0(ic)*b1(ic))
          end do
          do ir = 1,ng_d
            b(ir) = sum(pdr%nlpo(1:ng_i,ib)*t(1:ng_i,ng_i+ir)*dg)
            a(ir,:) = -t(ng_i+1:ng_o,ng_i+ir)*dg
            a(ir,ir) = a(ir,ir) + 0.5_double*pi*pdr%g(ng_i+ir)**2
          end do
          pdr%nlpo(ng_i+1:ng_o,ib) = solve(a,b) ; if (error()) goto 100
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
          t(1:ng_i,ib) = pdr%nlpo(1:ng_i,ib)
          pdr%nlpo(:,ib) = (2.0_double/pi)*pdr%nlpo(:,ib)*(pdr%g**2)*wt*dg
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
            f(ir) = radial_r_value_i(pdr,pdr%nlpo(:,ib),r(ir),pdr%l(ib))*r(ir)**2*wt(ir)*dr ; if (error()) goto 100
          end do
          do ig = 1,ng_i
            wg(ig) = sum(f*spherical_bessel(pdr%g(ig)*r,pdr%l(ib))) - t(ig,ib) ; if (error()) goto 100
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

        if (error("Exit ncp_data_mod::optimize_nlp_i")) continue

      end subroutine

      function cutoff_point_i(f) result(cp)
        real(double), dimension(:,:), intent(in) :: f
        integer :: cp
        real(double), parameter :: fl = 1.0e-10_double

        cp = size(f,1) - 1
        do while (cp > 0)
          if (maxval(abs(f(cp,:))) > fl) exit
          cp = cp - 1
        end do
        cp = cp + 1
        if (error(cp == 1,"ERROR: cutoff point was not found")) goto 100

100     if (error("Exit ncp_data_mod::cutoff_point_i")) continue

      end function

      end module

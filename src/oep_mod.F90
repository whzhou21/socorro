! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module oep_mod
!doc$ module oep_mod

!     One datatype is available here: type(oep_obj).
!     oep_mod encapsulates an iterative method for finding the Optimized Effective Potential.

      use kind_mod
      use mpi_mod
      use io_mod
      use tagio_mod
      use error_mod
      use diary_mod
      use arg_mod
      use ghost_mod
      use math_mod
      use layout_mod
      use grid_mod
      use symmetry_mod
      use xc_type_mod
      use xc_mod
      use electrons_sc_mod
      use wavefunctions_es_mod
      use multibasis_mod
      use multivector_mod
      use operators_mod
      use kpoints_mod
      use timing_mod

!cod$
      implicit none
      private

      ! normalization type
      integer, parameter :: ZERO = 0
      integer, parameter :: DXCP = 1

      ! preconditioner
      integer, parameter :: NONE = 0
      integer, parameter :: PC0  = 1
      integer, parameter :: PC1  = 2
      integer, parameter :: PC2  = 3

      ! optimizer
      integer, parameter :: CHEBYSHEV = 1
      integer, parameter :: SIMPLE    = 2

      type :: oep_rep
        integer :: ref
        type(ghost) :: g                                                ! ghost
        integer :: normalization                                        ! normalization type
        integer :: preconditioner                                       ! preconditioner
          real(double) :: pc_factor                                     ! preconditioner factor
          real(double) :: pc_power                                      ! preconditioner power factor
          real(double), dimension(:,:,:), pointer :: pc_array           ! preconditioning array
        integer :: optimizer                                            ! method used to optimize the oep
          type(chebyshev_state) :: cs                                   ! Chebyshev state (CHEBYSHEV)
          real(double) :: cs_cond                                       ! condition number (CHEBYSHEV)
          real(double) :: cs_norm                                       ! norm (CHEBYSHEV)
          integer :: cs_period                                          ! period (CHEBYSHEV)
          type(grid_obj) :: velocity                                    ! velocity (CHEBYSHEV)
          real(double) :: step_size                                     ! step size (SIMPLE)
        real(double) :: cutoff                                          ! xcp cutoff energy
        logical, dimension(:,:,:), pointer :: filter                    ! xcp cutoff filter (conformable to grid data with CDF_KIND)
        type(xc_obj) :: xc                                              ! exchange-correlation interface
        type(grid_obj) :: scp                                           ! self-consistent potential
        integer :: solver_dir                                           ! solver directions
        real(double) :: solver_tol                                      ! solver tolerance
        logical, dimension(:), pointer :: initialized                   ! indicates whether or not solutions is initialized
        type(multivector_obj), dimension(:), pointer :: solutions       ! solution multivectors
      end type

      type, public :: oep_obj
        private
        integer :: ref
        type(oep_rep), pointer :: o
      end type

!doc$
      public :: oep
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: step
      public :: diary

!cod$
      interface oep
        module procedure constructor_oep
      end interface
      interface update
        module procedure update_oep
      end interface
      interface my
        module procedure my_oep, my_new_oep
      end interface
      interface thy
        module procedure thy_oep
      end interface
      interface glean
        module procedure glean_oep
      end interface
      interface bequeath
        module procedure bequeath_oep
      end interface
      interface assignment(=)
        module procedure assign_oep
      end interface
      interface step
        module procedure step_oep
      end interface
      interface diary
        module procedure diary_oep
      end interface

      contains

! public routines

      function constructor_oep(xc,n_tot,n_val,hap,xcp,scp,restf) result(oep)
!doc$ function oep(xc,n_tot,n_val,hap,xcp,scp,restf) result(oep)
        type(xc_obj) :: xc
        type(grid_obj) :: n_tot, n_val, hap, xcp, scp
        type(tagio_obj), optional :: restf
        type(oep_obj) :: oep
!       requires: restf pointer be in FIELDS block.
!       modifies: xcp and scp
!       effects: Constructs a new oep. Initializes, xcp and scp.
!       errors: Solver parameters out of bounds. Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag
        real(double), parameter :: m1 = -1.0_double
        real(double), parameter :: p1 = +1.0_double
        type(layout_obj) :: lay

        call my(xc)
        call my(n_tot)
        call my(n_val)
        call my(hap)
        call my(xcp)
        call my(scp)
        if (present(restf)) call my(restf)

        oep%ref = 0
        allocate( oep%o )
        oep%o%ref = 0
        oep%o%g = x_ghost()

        call my(xc,oep%o%xc)

        call my(x_layout(scp),lay)

        ! form xcp filter
        call form_filter_i(oep%o,lay) ; if (error()) goto 100

        ! form xcp and scp
        select case(x_functional_dependence(oep%o%xc))
        case (FD_DENSITY)
          xcp = xc_potential(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
          call filter_i(oep%o,xcp)
          scp = hap
          call saxpby(p1,scp,p1,xcp)
        case (FD_ORBITAL)
          if (present(restf)) then
            tag = "SELF-CONSISTENT_POTENTIAL"
            call read_restart(scp,tag,restf) ; if (error()) goto 100
            xcp = scp
            call saxpby(p1,xcp,m1,hap)
            call filter_i(oep%o,xcp)
          else
            call empty(xcp)
            scp = hap
          end if
        case (FD_HYBRID)
          if (present(restf)) then
            tag = "SELF-CONSISTENT_POTENTIAL"
            call read_restart(scp,tag,restf) ; if (error()) goto 100
            xcp = scp
            call saxpby(p1,xcp,m1,hap)
            call filter_i(oep%o,xcp)
          else
            xcp = xc_potential(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
            call filter_i(oep%o,xcp)
            scp = hap
            call saxpby(p1,scp,p1,xcp)
          end if
        end select
        call my(scp,oep%o%scp)

        ! set normalization type
        select case(x_functional_dependence(oep%o%xc))
        case (FD_DENSITY)
          oep%o%normalization = DXCP
        case (FD_ORBITAL)
          oep%o%normalization = ZERO
        case (FD_HYBRID)
          oep%o%normalization = ZERO
        end select
        call arglc("oep_normalization",tag,found)
        if (found) then
          select case (trim(tag))
          case ("zero")
            if (oep%o%normalization == DXCP) then
              call warn("WARNING: Using ZERO oep_normalization with a density-dependent functional")
              oep%o%normalization = ZERO
            end if
          case ("dxcp")
            if (error(oep%o%normalization == ZERO,"ERROR: DXCP oep_normalization is not permitted")) goto 100
          case default
            if (error(.true.,"ERROR: oep_normalization was not recognized")) goto 100
          end select
        end if

        ! get optimization method and initialize
        call arglc("oep_optimizer",tag,found)
        if (.not.found) tag = "chebyshev"
        select case (trim(tag))
        case ("chebyshev","c")
          oep%o%optimizer = CHEBYSHEV
          call chebyshev_initialize_i(oep%o,lay) ; if (error()) goto 100
        case ("simple","s")
          oep%o%optimizer = SIMPLE
          call simple_initialize_i(oep%o) ; if (error()) goto 100
        case default
          if (error(.true.,"ERROR: oep_optimizer tag was not recognized")) goto 100
        end select

        ! get the oep preconditioner
        call form_preconditioner_i(oep%o,lay) ; if (error()) goto 100

        ! get solver parameters and initialize workspace
        call arg("oep_solver_dir",oep%o%solver_dir,found)
        if (.not.found) oep%o%solver_dir = 50
        if (error(oep%o%solver_dir < 1,"ERROR: oep_solver_dir < 1")) goto 100
        call arg("oep_solver_tol",oep%o%solver_tol,found)
        if (.not.found) oep%o%solver_tol = 1.0e-8_double
        if (error(oep%o%solver_tol < 0.0_double,"ERROR: oep_solver_tol < 0")) goto 100
        nullify( oep%o%initialized, oep%o%solutions )

100     call glean(thy(lay))

        call glean(thy(xc))
        call glean(thy(n_tot))
        call glean(thy(n_val))
        call glean(thy(hap))
        call glean(thy(xcp))
        call glean(thy(scp))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit oep_mod::constructor_oep")) continue

      end function

      subroutine update_oep(oep,xc,n_tot,n_val,hap,xcp,scp)
!doc$ subroutine update(oep,xc,n_tot,n_val,hap,xcp,scp)
        type(oep_obj) :: oep
        type(xc_obj) :: xc
        type(grid_obj) :: n_tot, n_val, hap, xcp, scp
!       modifies: oep, xcp, and scp
!       effects: Updates oep. Initializes xcp and scp.
!       errors: Passes errors.
!       requires: xc changes are confined to its layout and space group.

!cod$
        integer :: i
        real(double), parameter :: m1 = -1.0_double
        real(double), parameter :: p1 = 1.0_double

        call my(oep)
        call my(xc)
        call my(n_tot)
        call my(n_val)
        call my(hap)
        call my(xcp)
        call my(scp)

        call own_i(oep)
        oep%o%g = x_ghost()

        ! form xcp and scp
        select case(x_functional_dependence(oep%o%xc))
        case (FD_DENSITY)
          xcp = xc_potential(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
          scp = hap
          call saxpby(p1,scp,p1,xcp)
          call filter_i(oep%o,xcp)
        case (FD_ORBITAL)
          call empty(xcp)
          scp = hap
        case (FD_HYBRID)
          xcp = xc_potential(oep%o%xc,n_tot,n_val); if (error()) goto 100
          scp = hap
          call saxpby(p1,scp,p1,xcp)
          call filter_i(oep%o,xcp)
        end select
        oep%o%scp = scp

        select case (oep%o%optimizer)
        case (CHEBYSHEV)
          call chebyshev_initialize(oep%o%cs,oep%o%cs_norm,oep%o%cs_cond,oep%o%cs_period)
          call empty(oep%o%velocity)
        end select

        ! initialize work space
        do i = 1,size(oep%o%initialized)
          if (oep%o%initialized(i)) call glean(thy(oep%o%solutions(i)))
        end do
        deallocate( oep%o%initialized, oep%o%solutions )

100     call glean(thy(oep))
        call glean(thy(xc))
        call glean(thy(n_tot))
        call glean(thy(n_val))
        call glean(thy(hap))
        call glean(thy(xcp))
        call glean(thy(scp))

        if (error("Exit oep_mod::update_oep")) continue

      end subroutine

      subroutine my_oep(oep)
!doc$ subroutine my(oep)
        type(oep_obj) :: oep

!cod$
        oep%ref = oep%ref + 1
        oep%o%ref = oep%o%ref + 1
      end subroutine

      subroutine my_new_oep(oepi,oep)
!doc$ subroutine my(oepi,oep)
        type(oep_obj) :: oepi, oep

!cod$
        oep%ref = 1
        oep%o => oepi%o
        oep%o%ref = oep%o%ref + 1
      end subroutine

      function thy_oep(oep) result(oepo)
!doc$  function thy(oep) result(oepo)
        type(oep_obj) :: oep, oepo

!cod$
        oep%ref = oep%ref - 1
        oep%o%ref = oep%o%ref - 1
        oepo%ref = oep%ref
        oepo%o => oep%o
      end function
      
      subroutine glean_oep(oep)
!doc$ subroutine glean(oep)
        type(oep_obj) :: oep

!cod$
        integer :: i
        if (oep%o%ref < 1) then
          if (associated( oep%o%pc_array )) deallocate( oep%o%pc_array )
          select case (oep%o%optimizer)
          case (CHEBYSHEV)
            call glean(thy(oep%o%velocity))
          end select
          if (associated( oep%o%filter )) deallocate( oep%o%filter )
          call glean(thy(oep%o%xc))
          call glean(thy(oep%o%scp))
          if (associated( oep%o%initialized )) then
            do i = 1,size(oep%o%initialized)
              if (oep%o%initialized(i)) call glean(thy(oep%o%solutions(i)))
            end do
            deallocate( oep%o%initialized, oep%o%solutions )
          end if
          deallocate( oep%o )
        end if
      end subroutine

      subroutine bequeath_oep(oep)
!doc$ subroutine bequeath(oep)
        type(oep_obj) :: oep

!cod$
        continue
      end subroutine
    
      subroutine assign_oep(oep,oep2) 
!doc$ subroutine assign(oep,oep2) 
        type(oep_obj), intent(inout) :: oep
        type(oep_obj), intent(in) :: oep2

!cod$
        type(oep_obj) :: oept
        call my(oep2)
        oept%o => oep%o
        oep%o%ref = oep%o%ref - oep%ref
        oep%o => oep2%o
        oep%o%ref = oep%o%ref + oep%ref
        call glean(oept)
        call glean(thy(oep2))
      end subroutine

      subroutine step_oep(oep,el,n_tot,n_val,hap,xcp,scp,e)
!doc$ subroutine step(oep,el,n_tot,n_val,hap,xcp,scp,e)
        type(oep_obj) :: oep
        type(electrons_sc_obj) :: el
        type(grid_obj) :: n_tot, n_val, hap, xcp, scp
        real(double), intent(out) :: e
!       modifies: oep, xcp, scp, and e
!       errors: Passes errors.

!cod$
        real(double), parameter :: m1 = -1.0_double
        real(double), parameter :: p1 = +1.0_double
        real(double) :: e_span
        complex(double) :: scp_norm
        complex(double), dimension(:,:,:), pointer :: c1
        type(layout_obj) :: lay
        type(grid_obj) :: de_drho, de_dv

        call my(oep)
        call my(el)
        call my(n_tot)
        call my(n_val)
        call my(hap)
        call my(xcp)
        call my(scp)

        nullify( c1 )

        if (error(.not.thermal_occupations(el),"ERROR: thermal occupations must be used with the OEP method")) goto 200

        call own_i(oep)
        oep%o%g = x_ghost()

        call my(x_layout(scp),lay)

        ! form initial e
        select case(x_functional_dependence(oep%o%xc))
        case (FD_DENSITY)
          e = xc_energy(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
        case (FD_ORBITAL)
          e = 0.0_double
        case (FD_HYBRID)
          e = xc_energy(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
        end select

        ! form initial xcp
        select case(x_functional_dependence(oep%o%xc))
        case (FD_DENSITY)
          xcp = xc_potential(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
        case (FD_ORBITAL)
          call empty(xcp)
        case (FD_HYBRID)
          xcp = xc_potential(oep%o%xc,n_tot,n_val) ; if (error()) goto 100
        end select

        call filter_i(oep%o,xcp)

        ! form initial scp
        scp = hap
        call saxpby(p1,scp,p1,xcp)

        ! get initial normalization
        select case (oep%o%normalization)
        case (ZERO)
          scp_norm = (0.0_double,0.0_double)
        case (DXCP)
          call get_normalization(oep%o%scp,scp_norm)
        end select

        ! form de_drho
        call my(scp,de_drho)
        call saxpby(p1,de_drho,m1,oep%o%scp)
        call sgroup_to_kgroup(de_drho)

        ! form de_dv and augment e
        call my(grid(lay,KGROUP),de_dv)
        call solver_i(oep%o,el,de_drho,xcp,de_dv,e) ; if (error()) goto 100

        ! precondition de_dv
        select case (oep%o%optimizer)
        case (CHEBYSHEV)
          e_span = x_fermi_level(el) - minval(x_eigenvalues(el))
          call take(c1,de_dv,CDF_KIND)
          c1 = c1*oep%o%pc_array/e_span
          call put(c1,de_dv,CDF_KIND)
        case (SIMPLE)
          select case (oep%o%preconditioner)
          case (PC1,PC2)
            call take(c1,de_dv,CDF_KIND)
            c1 = c1*oep%o%pc_array
            call put(c1,de_dv,CDF_KIND)
          end select
        end select

        ! update scp
        select case (oep%o%optimizer)
        case (CHEBYSHEV)
          call chebyshev_step_i(oep%o,de_dv) ; if (error()) goto 100
        case (SIMPLE)
          call simple_step_i(oep%o,de_dv)    ; if (error()) goto 100
        end select

        ! reset scp normalization
        call set_normalization(oep%o%scp,scp_norm)

        ! form final xcp
        xcp = oep%o%scp
        call saxpby(p1,xcp,m1,hap)
        call filter_i(oep%o,xcp)

        ! form final scp
        scp = hap
        call saxpby(p1,scp,p1,xcp)
        oep%o%scp = scp

        if (associated( c1 )) deallocate( c1 )

100     call glean(thy(lay))
        call glean(thy(de_drho))
        call glean(thy(de_dv))

200     call glean(thy(oep))
        call glean(thy(el))
        call glean(thy(n_tot))
        call glean(thy(n_val))
        call glean(thy(hap))
        call glean(thy(xcp))
        call glean(thy(scp))

        if (error("Exit oep_mod::step_oep")) continue

      end subroutine

      subroutine diary_oep(oep)
!doc$ subroutine diary(oep)
        type(oep_obj) :: oep
!       requires: Valid oep.
!       effects: Writes oep information to the diary.

!cod$
        call my(oep)

        if (i_access(diaryfile())) then
          select case (oep%o%optimizer)
          case (CHEBYSHEV)
            write(x_unit(diaryfile()),'(t8,"Chebyshev acceleration for OEP solution:")')
            write(x_unit(diaryfile()),'(t10,"condition number = ",f0.2)') oep%o%cs_cond
            write(x_unit(diaryfile()),'(t10,"norm = ",f0.2)') oep%o%cs_norm
            write(x_unit(diaryfile()),'(t10,"reset period = ",i0)') oep%o%cs_period
          case (SIMPLE)
            write(x_unit(diaryfile()),'(t8,"Simple updates for OEP solution:")')
            write(x_unit(diaryfile()),'(t10,"step size = ",f0.3)') oep%o%step_size
          end select
          select case (oep%o%preconditioner)
          case (PC1)
            write(x_unit(diaryfile()),'(t8,"Using the pc1 preconditioner:")')
            write(x_unit(diaryfile()),'(t10,"pc factor = ",f0.2)') oep%o%pc_factor
            write(x_unit(diaryfile()),'(t10,"pc power  = ",f0.2)') oep%o%pc_power
          case (PC2)
            write(x_unit(diaryfile()),'(t8,"Using the pc2 preconditioner:")')
            write(x_unit(diaryfile()),'(t10,"pc factor = ",f0.2)') oep%o%pc_factor
            write(x_unit(diaryfile()),'(t10,"pc power  = ",f0.2)') oep%o%pc_power
          end select
          select case (oep%o%normalization)
          case (ZERO)
            write(x_unit(diaryfile()),'(t8,"Normalizing the self-consistent potential to 0")')
          case (DXCP)
            write(x_unit(diaryfile()),'(t8,"Normalizing the self-consistent potential to the density-derived XC potential")')
          end select
          write(x_unit(diaryfile()),'(t8,"Filtering the exchange-correlation potential at ",f0.2," Ryd")') oep%o%cutoff
          write(x_unit(diaryfile()),'(t8,"Conjugate-gradients solver for dE/dV:")')
          write(x_unit(diaryfile()),'(t10,"maximum number of directions = ",i0)') oep%o%solver_dir
          write(x_unit(diaryfile()),'(t10,"tolerance for residual norm = ",es8.2)') oep%o%solver_tol
        end if

        call glean(thy(oep))

      end subroutine

! private routines

      subroutine own_i(oep)
        type(oep_obj), intent(inout) :: oep
!       notes: oept%o%cs is not set here and therefore must be done in the calling routine.
        type(oep_obj) :: oept
        integer :: i, n
        if (oep%ref < oep%o%ref) then
          allocate( oept%o )
          oept%o%ref = 0
          oept%o%g = oep%o%g
          oept%o%normalization = oep%o%normalization
          oept%o%preconditioner = oep%o%preconditioner
          oept%o%pc_factor = oep%o%pc_factor
          oept%o%pc_power = oep%o%pc_power
          select case(oep%o%preconditioner)
          case (PC0,PC1,PC2)
            allocate( oept%o%pc_array(size(oep%o%pc_array,1),size(oep%o%pc_array,2),size(oep%o%pc_array,3)) )
            oept%o%pc_array = oep%o%pc_array
          end select
          oept%o%optimizer = oep%o%optimizer
          select case (oept%o%optimizer)
          case (CHEBYSHEV)
            call my(oep%o%velocity,oept%o%velocity)
          case (SIMPLE)
            oept%o%step_size = oep%o%step_size
          end select
          oept%o%cutoff        = oep%o%cutoff
          allocate( oept%o%filter(size(oep%o%filter,1),size(oep%o%filter,2),size(oep%o%filter,3)) )
          oept%o%filter        = oep%o%filter
          call my(oep%o%xc,oept%o%xc)
          call my(oep%o%scp,oept%o%scp)
          oept%o%solver_dir    = oep%o%solver_dir
          oept%o%solver_tol    = oep%o%solver_tol
          if (associated( oep%o%initialized )) then
            n = size(oep%o%initialized)
            allocate( oept%o%initialized(n) )
            oept%o%initialized = oep%o%initialized
            allocate( oept%o%solutions(n) )
            do i = 1,n
              if (oep%o%initialized(i)) call my(oep%o%solutions(i),oept%o%solutions(i))
            end do
          else
            nullify( oept%o%initialized )
            nullify( oept%o%solutions )
          end if
          oep%o%ref = oep%o%ref - oep%ref
          oep%o => oept%o
          oep%o%ref = oep%o%ref + oep%ref
        end if
      end subroutine

      subroutine solver_i(oepr,el,de_drho,xcp,de_dv,e)
        type(oep_rep) :: oepr
        type(electrons_sc_obj) :: el
        type(grid_obj) :: de_drho, de_dv, xcp
        real(double), intent(inout) :: e

        character(line_len), parameter :: init = "zeros"
        real(double), parameter :: p1 = 1.0_double, m1 = -1.0_double
        complex(double), parameter :: cm1 = (-1.0_double,0.0_double)
        complex(double), parameter :: c0  = ( 0.0_double,0.0_double)
        complex(double), parameter :: cp1 = (+1.0_double,0.0_double)
        integer :: ib, ib1, ib2, ik, iks, it, nb, nk, nks
        real(double) :: alpha, ek, kt, residual_norm, kg_ek, kg_tj, kg_tje, tj, tje, tj_ratio, weight, x
        real(double), dimension(:), allocatable :: me_vec, occ_vec, p1_vec, rnorm_vec
        real(double), dimension(:,:), allocatable :: eigs, gksc, kg_gksc, occs, def
        complex(double) :: alpha_c, cm1w, cp1w
        complex(double), dimension(:), allocatable :: a_vec, b_vec, bn_vec, bd_vec, cp1_vec, wrk_vec
        complex(double), dimension(:,:), allocatable :: e_mat, j_mat
        type(kpoints_obj) :: kp
        type(hamiltonian_obj) :: h
        type(multibasis_obj) :: mb
        type(multivector_obj) :: e_mv, t_mv, r_mv, pr_mv, sd_mv, wrk_mv

        call my(el)
        call my(de_drho)
        call my(xcp)
        call my(de_dv)

        call my(x_kpoints(el),kp)

        nk = x_n_kpoints(kp)
        nb = x_n_bands(el)

        allocate( eigs(nk,nb), occs(nk,nb), def(nk,nb) )
        eigs = x_eigenvalues(el)
        occs = x_occupations(el)
        def = eigs - x_fermi_level(el)

        allocate( gksc(nk,nb), kg_gksc(nk,nb) )

        kt = x_kt(el)

        allocate( me_vec(nb), occ_vec(nb), p1_vec(nb), rnorm_vec(nb) )
        p1_vec = p1

        allocate( a_vec(nb), b_vec(nb), bn_vec(nb), bd_vec(nb), cp1_vec(nb), wrk_vec(nb) )
        cp1_vec = cp1

        allocate( e_mat(nb,nb), j_mat(nb,nb) )

        if (.not.associated( oepr%initialized )) then
          nks = 0
          do ik = 1,nk
            if (mpi_mykgroup() == x_kgroup_index(el,ik)) nks = nks + 1
          end do
          allocate( oepr%initialized(nks), oepr%solutions(nks) )
          oepr%initialized = .false.
        end if

        kg_tj = 0.0_double
        kg_tje = 0.0_double

        kg_gksc = 0.0_double

        kg_ek = 0.0_double

        iks = 0

        select case(x_functional_dependence(oepr%xc))
        case (FD_ORBITAL,FD_HYBRID)
          if (mpi_nkgroups() > 1) call distribute(el)
        end select

        do ik = 1,nk

          if (mpi_mykgroup() /= x_kgroup_index(el,ik)) cycle

          weight = x_kweight(kp,ik)

          me_vec = -eigs(ik,:)
          occ_vec = occs(ik,:)

          call my(x_multivector(x_wavefunctions(el,ik)),e_mv)
          call my(hamiltonian(x_h_kpoint(x_wavefunctions(el,ik))),h)

          call my(x_multibasis(e_mv),mb)

          ! init ='zeros' so t_mv = 0
          call my(multivector(mb,init),t_mv)

          select case(x_functional_dependence(oepr%xc))
          case (FD_DENSITY)
            alpha = 0.0_double
            call apply_field(c0,t_mv,cp1,e_mv,de_drho)
          case (FD_ORBITAL)
            alpha = 1.0_double
            call xc_energy_and_derivative(oepr%xc,el,ik,kg_ek,t_mv)
            call apply_field(cp1,t_mv,cp1,e_mv,de_drho)
          case(FD_HYBRID)
            alpha = x_hybrid_mixing(oepr%xc)
            alpha_c = cmplx(alpha,0.0,double)
            call xc_energy_and_derivative(oepr%xc,el,ik,kg_ek,t_mv)
            call apply_field(alpha_c,t_mv,cp1,e_mv,de_drho)
          end select

          call project(e_mv,t_mv,e_mat)

          do ib1 = 1,nb
            do ib2 = 1,nb
              x = fermi_distr_divdiff(def(ik,ib1),def(ik,ib2),kt)
              j_mat(ib1,ib2) = x*e_mat(ib1,ib2)
            end do
            x = weight*fermi_distr_divdiff(def(ik,ib1),def(ik,ib1),kt)
            kg_tj = kg_tj + x
            kg_tje = kg_tje + x*real(e_mat(ib1,ib1),double)
            kg_gksc(ik,ib1) = real(e_mat(ib1,ib1),double)
          end do

          ! initialize the solutions and the residual
          iks = iks + 1
          if (.not.oepr%initialized(iks)) then
            call my(multivector(mb,init),oepr%solutions(iks))           ; if (error()) goto 100
            oepr%initialized(iks) = .true.
            call my(multivector(mb,init),r_mv)                          ; if (error()) goto 100
          else
            call project(e_mv,oepr%solutions(iks))
            call my(multivector(mb,init),r_mv)                          ; if (error()) goto 100
            call apply_hamiltonian(oepr%solutions(iks),r_mv,h)
            call combine(p1_vec,r_mv,me_vec,oepr%solutions(iks))
          end if
          call combine(m1,r_mv,p1,t_mv)
          call project(e_mv,r_mv)

          call my(multivector(mb,init),pr_mv)  ; if (error()) goto 100
          call my(multivector(mb,init),sd_mv)  ; if (error()) goto 100
          call my(multivector(mb,init),wrk_mv) ; if (error()) goto 100

          it = 0
          do

            it = it + 1

            pr_mv = r_mv
            call precondition(pr_mv,t_mv,h)

            call project(e_mv,pr_mv)

            call multiply(r_mv,pr_mv,bn_vec)
            if (it == 1) then
              b_vec = (0.0_double,0.0_double)
            else
              b_vec = bn_vec/bd_vec
            end if
            call combine(b_vec,sd_mv,cp1_vec,pr_mv)
            bd_vec = bn_vec

            call apply_hamiltonian(sd_mv,wrk_mv,h)
            call combine(p1_vec,wrk_mv,me_vec,sd_mv)
            call project(e_mv,wrk_mv)
            call multiply(sd_mv,wrk_mv,wrk_vec)

            a_vec = bn_vec/wrk_vec
            call combine(cp1_vec,oepr%solutions(iks),a_vec,sd_mv)
            a_vec = -bn_vec/wrk_vec
            call combine(cp1_vec,r_mv,a_vec,wrk_mv)

            call multiply(r_mv,rnorm_vec)
            residual_norm = sum(occ_vec*sqrt(rnorm_vec))/real(nb,double)

            if ((it == oepr%solver_dir) .or. (residual_norm <= oepr%solver_tol)) exit

          end do

          wrk_mv = oepr%solutions(iks)
          call portion(occ_vec,wrk_mv)

          cm1w = cm1*weight
          cp1w = cp1*weight
          call transform(cm1w,wrk_mv,cp1w,e_mv,j_mat)

          call add_grid_density(wrk_mv,e_mv,de_dv) ; if (error()) goto 100

          call glean(thy(h))
          call glean(thy(mb))
          call glean(thy(e_mv))
          call glean(thy(t_mv))
          call glean(thy(r_mv))
          call glean(thy(pr_mv))
          call glean(thy(sd_mv))
          call glean(thy(wrk_mv)) ; if (error()) goto 100

        end do

        select case(x_functional_dependence(oepr%xc))
        case (FD_ORBITAL,FD_HYBRID)
          if (mpi_nkgroups() > 1) call release(el)
        end select

        call xcomm_allreduce(XKGROUP,MPI_SUM,kg_ek,ek)
        e = e + alpha*ek

        call xcomm_allreduce(XKGROUP,MPI_SUM,kg_tj,tj)
        call xcomm_allreduce(XKGROUP,MPI_SUM,kg_tje,tje)
        tj_ratio = tje/tj
        j_mat = (0.0_double,0.0_double)

        call xcomm_allreduce(XKGROUP,MPI_SUM,kg_gksc,gksc)
        call s_oep_gksc(el,gksc)

        do ik = 1,nk

          if (mpi_mykgroup() /= x_kgroup_index(el,ik)) cycle

          weight = x_kweight(kp,ik)

          do ib = 1,nb
            x = fermi_distr_divdiff(def(ik,ib),def(ik,ib),kt)
            j_mat(ib,ib) = cmplx(-x*tj_ratio,0,double)
          end do

          call my(x_multivector(x_wavefunctions(el,ik)),e_mv)    ; if (error()) goto 100
          call my(multivector(x_multibasis(e_mv),init),wrk_mv)   ; if (error()) goto 100

          cp1w = cp1*weight
          call transform(c0,wrk_mv,cp1w,e_mv,j_mat)

          call add_grid_density(wrk_mv,e_mv,de_dv) ; if (error()) goto 100

          call glean(thy(e_mv))
          call glean(thy(wrk_mv)) ; if (error()) goto 100

        end do

!       if (i_access(diaryfile())) write(x_unit(diaryfile()),'("norm de_dv = ",f10.6)') norm(de_dv)

        call merge_grid_density(de_dv)
        call symmetrize_grid(x_space_group(oepr%xc),de_dv)
        call filter_i(oepr,de_dv)

        if (allocated( me_vec )) deallocate( me_vec )
        if (allocated( occ_vec )) deallocate( occ_vec )
        if (allocated( p1_vec )) deallocate( p1_vec )
        if (allocated( rnorm_vec )) deallocate( rnorm_vec )
        if (allocated( eigs )) deallocate( eigs )
        if (allocated( gksc )) deallocate( gksc )
        if (allocated( kg_gksc )) deallocate( kg_gksc )
        if (allocated( occs )) deallocate( occs )
        if (allocated( def )) deallocate( def )
        if (allocated( a_vec )) deallocate( a_vec )
        if (allocated( b_vec )) deallocate( b_vec )
        if (allocated( bn_vec )) deallocate( bn_vec )
        if (allocated( bd_vec )) deallocate( bd_vec )
        if (allocated( cp1_vec )) deallocate( cp1_vec )
        if (allocated( wrk_vec )) deallocate( wrk_vec )
        if (allocated( e_mat )) deallocate( e_mat )
        if (allocated( j_mat )) deallocate( j_mat )

100     call glean(thy(kp))

200     call glean(thy(el))
        call glean(thy(de_drho))
        call glean(thy(xcp))
        call glean(thy(de_dv))

        if (error("Exit oep_mod::solver_i")) continue

      end subroutine

      subroutine form_filter_i(oepr,lay)
        type(oep_rep) :: oepr
        type(layout_obj) :: lay

        logical :: found
        real(double), dimension(:,:,:), pointer :: g2

        call my(lay)

        nullify( g2 )

        call arg("oep_cutoff",oepr%cutoff,found)
        if (.not.found) then
          call arg("wf_cutoff",oepr%cutoff,found)
          if (error(.not.found,"ERROR: oep_cutoff was not identified")) goto 100
        end if

        call alloc(oepr%filter,lay,D_TYPE,SGROUP)
        call fdel(g2,lay,D_TYPE,SGROUP)
        where (g2 <= oepr%cutoff)
          oepr%filter = .false.
        elsewhere
          oepr%filter = .true.
        end where

        if (associated( g2 )) deallocate( g2 )

100     call glean(thy(lay))

        if (error("Exit oep_mod::form_filter_i")) continue

      end subroutine

      subroutine filter_i(oepr,g)
        type(oep_rep) :: oepr
        type(grid_obj) :: g

        complex(double), dimension(:,:,:), pointer :: c1

        call my(g)

        if (x_type(g) /= EMPTY_KIND) then
          call take(c1,g,CDF_KIND)
          where (oepr%filter) c1 = (0.0_double,0.0_double)
          call put(c1,g,CDF_KIND)
        end if

        call glean(thy(g))

      end subroutine

      subroutine form_preconditioner_i(oepr,lay)
        type(oep_rep) :: oepr
        type(layout_obj) :: lay

        logical :: found
        character(line_len) :: tag
        real(double) :: p
        real(double), dimension(:,:,:), pointer :: r1

        call my(lay)

        nullify( r1 )

        nullify( oepr%pc_array )

        ! get the oep preconditioner
        select case (oepr%optimizer)
        case (CHEBYSHEV)
          oepr%preconditioner = PC0
          oepr%pc_factor = 0.0_double
          oepr%pc_power = 0.0_double
          call fdel(r1,lay,D_TYPE,SGROUP)
          call filter(r1,lay,SGROUP)
          call alloc(oepr%pc_array,lay,D_TYPE,SGROUP)
          oepr%pc_array = 0.75_double*r1
        case (SIMPLE)
          call arglc("oep_preconditioner",tag,found)
          if (.not.found) tag = "pc1"
          select case (trim(tag))
          case ("none")
            oepr%preconditioner = NONE
            oepr%pc_factor = 0.0_double
            oepr%pc_power = 0.0_double
          case ("pc1")
            oepr%preconditioner = PC1
            call arg("oep_pc_factor",oepr%pc_factor,found)
            if (.not.found) oepr%pc_factor = 1.0_double
            if (error(oepr%pc_factor <= 0.0_double,"ERROR: oep_pc_factor <= 0")) goto 100
            call arg("oep_pc_power",oepr%pc_power,found)
            if (.not.found) oepr%pc_power = 2.0_double
            if (error(oepr%pc_power <= 0.0_double,"ERROR: oep_pc_power <= 0")) goto 100
            call fdel(r1,lay,D_TYPE,SGROUP)
            call filter(r1,lay,SGROUP)
            p = oepr%pc_power/2.0_double
            r1 = r1**p
            call alloc(oepr%pc_array,lay,D_TYPE,SGROUP)
            oepr%pc_array = 1.0_double + oepr%pc_factor*r1
          case ("pc2")
            oepr%preconditioner = PC2
            call arg("oep_pc_factor",oepr%pc_factor,found)
            if (.not.found) oepr%pc_factor = 1.0_double
            if (error(oepr%pc_factor <= 0.0_double,"ERROR: oep_pc_factor <= 0")) goto 100
            call arg("oep_pc_power",oepr%pc_power,found)
            if (.not.found) oepr%pc_power = 2.0_double
            if (error(oepr%pc_power <= 0.0_double,"ERROR: oep_pc_power <= 0")) goto 100
            call fdel(r1,lay,D_TYPE,SGROUP)
            call filter(r1,lay,SGROUP)
            p = oepr%pc_power/2.0_double
            r1 = r1**p
            call alloc(oepr%pc_array,lay,D_TYPE,SGROUP)
            oepr%pc_array = oepr%pc_factor*r1
          case default
            if (error(.true.,"ERROR: oep_preconditioner tag was not recognized")) goto 100
          end select
        end select

100     if (associated( r1 )) deallocate( r1 )

        call glean(thy(lay))

        if (error("Exit oep_mod::form_preconditioner_i")) continue

      end subroutine

      subroutine chebyshev_initialize_i(oepr,lay)
        type(oep_rep) :: oepr
        type(layout_obj) :: lay

        logical :: found

        call my(lay)

        call arg("oep_chebyshev_norm",oepr%cs_norm,found)
        if (.not.found) oepr%cs_norm = 1.0_double
        call arg("oep_chebyshev_cond",oepr%cs_cond,found)
        if (.not.found) oepr%cs_cond = 1000.0_double
        call arg("oep_chebyshev_period",oepr%cs_period,found)
        if (.not.found) oepr%cs_period = 1 + floor(2.0_double*sqrt(oepr%cs_cond))
        call chebyshev_initialize(oepr%cs,oepr%cs_norm,oepr%cs_cond,oepr%cs_period)
        call my(grid(lay,SGROUP),oepr%velocity)

        call glean(thy(lay))

        if (error("Exit oep_mod::chebyshev_initialize_i")) continue

      end subroutine

      subroutine chebyshev_step_i(oepr,de_dv)
        type(oep_rep) :: oepr
        type(grid_obj) :: de_dv
        ! requires: de_dv be filtered

        real(double), parameter :: p1 = 1.0_double
        real(double) :: mp1, mp2

        call chebyshev_accelerate(oepr%cs,mp1,mp2)
!       if (i_access(diaryfile())) write(x_unit(diaryfile()),'(t7,"mp1 = ",f10.7,"; mp2 = ",f10.7)') mp1, mp2
        call saxpby(mp1,oepr%velocity,mp2,de_dv)
        call saxpby(p1,oepr%scp,p1,oepr%velocity)

        if (error("Exit oep_mod::chebyshev_step_i")) continue

      end subroutine

      subroutine simple_initialize_i(oepr)
        type(oep_rep) :: oepr

        logical :: found

        call arg("oep_step_size",oepr%step_size,found)
        if (.not.found) oepr%step_size = 4.0_double
        if (error(oepr%step_size <= 0.0_double,"ERROR: oep_step_size <= 0")) goto 100

100     if (error("Exit oep_mod::simple_initialize_i")) continue

      end subroutine

      subroutine simple_step_i(oepr,de_dv)
        type(oep_rep) :: oepr
        type(grid_obj) :: de_dv
        ! requires: de_dv be filtered

        real(double), parameter :: p1 = +1.0_double
        real(double) :: step

        step = -oepr%step_size
        call saxpby(p1,oepr%scp,step,de_dv)

        if (error("Exit oep_mod::simple_step_i")) continue

      end subroutine

      end module

! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module xc_orbital_mod
!doc$ module xc_orbital_mod

!     One datatype is available here: type(xc_orbital_obj)

!     xc_orbital_mod encapsulates routines for computing exchange-correlation quantities
!     using functionals that depend explicitly on the Kohn-Sham orbitals.

      use electrons_sc_mod
      use error_mod
      use ghost_mod
      use kind_mod
      use kpoints_mod
      use math_mod
      use mpi_mod
      use multivector_mod
      use lattice_mod
      use layout_mod
      use operators_mod
      use path_mod
      use timing_mod
      use symmetry_mod
      use wavefunctions_es_mod
      use xc_type_mod

!cod$
      implicit none
      private

      type, public :: xc_orbital_rep
        integer :: ref
        type(ghost) :: g                                           ! ghost
        integer :: coulomb_kernel                                  ! EXX Coulomb kernel type
        integer :: auxiliary_type                                  ! EXX auxiliary type
        integer :: sd_aux_form                                     ! EXX structure-dependent auxiliary form
        integer :: si_aux_form                                     ! EXX structure-independent auxiliary form
        type(xc_type_obj) :: xct                                   ! xc_type object
        type(layout_obj) :: lay                                    ! layout object
        type(space_group_obj) :: sg                                ! space group object
      end type

      type, public :: xc_orbital_obj
        private
        integer :: ref
        type(xc_orbital_rep), pointer :: o
      end type

!doc$
      public :: xc_orbital
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_xc_type
      public :: x_layout
      public :: x_space_group
      public :: exx_energy
      public :: exx_derivative
      public :: exx_energy_and_derivative

!cod$
      interface xc_orbital
        module procedure constructor_xco
      end interface
      interface update
        module procedure update_xco
      end interface
      interface my
        module procedure my_xco, my_new_xco
      end interface
      interface thy
        module procedure thy_xco
      end interface
      interface glean
        module procedure glean_xco
      end interface
      interface bequeath
        module procedure bequeath_xco
      end interface
      interface assignment(=)
        module procedure assign_xco
      end interface
      interface x_ref
        module procedure xco_ref
      end interface
      interface x_ghost
        module procedure xco_ghost
      end interface
      interface x_xc_type
        module procedure xco_xc_type
      end interface
      interface x_layout
        module procedure xco_layout
      end interface
      interface x_space_group
        module procedure xco_space_group
      end interface
      interface exx_energy
        module procedure exx_energy_xco
      end interface
      interface exx_derivative
        module procedure exx_derivative_xco
      end interface
      interface exx_energy_and_derivative
        module procedure exx_energy_and_derivative_xco
      end interface

      contains

! public routines

      function constructor_xco(xct,lay,sg) result(xco)
!doc$ function xc_orbital(xct,lay,sg) result (xco)
        type(xc_type_obj) :: xct
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
        type(xc_orbital_obj) :: xco
!       requires:
!       effects: Constructs a new xco.

!cod$
        call my(xct)
        call my(lay)
        call my(sg)

        xco%ref = 0
        allocate( xco%o )
        xco%o%ref = 0
        xco%o%g = x_ghost()

        call my(xct,xco%o%xct)
        call my(lay,xco%o%lay)
        call my(sg,xco%o%sg)

        ! get parameters needed for EXX calculations
        xco%o%coulomb_kernel = x_coulomb_kernel(xct)
        xco%o%auxiliary_type = x_auxiliary_type(xct)
        xco%o%sd_aux_form = x_sd_aux_form(xct)
        xco%o%si_aux_form = x_si_aux_form(xct)

        call glean(thy(xct))
        call glean(thy(lay))
        call glean(thy(sg))

        if (error("Exit xc_orbital_mod::constructor_xco")) continue

      end function

      subroutine update_xco(xco,lay,sg)
!doc$ subroutine update(xco,lay,sg)
        type(xc_orbital_obj) :: xco
        type(layout_obj) :: lay
        type(space_group_obj) :: sg
!       modifies: xco
!       effects: Currently updates lay and sg, but not xct.
!       errors: Passes errors.

!cod$
        call my(xco)
        call my(lay)
        call my(sg)

        xco%o%lay = lay
        xco%o%sg = sg

100     call glean(thy(xco))
        call glean(thy(lay))
        call glean(thy(sg))

        if (error("Exit xc_orbital_mod::update_xco")) continue

      end subroutine

      subroutine my_xco(xco)
!doc$ subroutine my(xco)
        type(xc_orbital_obj) :: xco

!cod$
        xco%ref = xco%ref + 1
        xco%o%ref = xco%o%ref + 1
      end subroutine

      subroutine my_new_xco(xcoi,xco)
!doc$ subroutine my(xcoi,xco)
        type(xc_orbital_obj) :: xcoi, xco

!cod$
        xco%ref = 1
        xco%o => xcoi%o
        xco%o%ref = xco%o%ref + 1
      end subroutine

      function thy_xco(xco) result(xcoo)
!doc$ function thy(xco) result(xcoo)
        type(xc_orbital_obj) :: xco, xcoo

!cod$
        xco%ref = xco%ref - 1
        xco%o%ref = xco%o%ref - 1
        xcoo%ref = xco%ref
        xcoo%o => xco%o
      end function

      subroutine glean_xco(xco)
!doc$ subroutine glean(xco)
        type(xc_orbital_obj) :: xco

!cod$
        integer :: ik
        if (xco%o%ref < 1) then
          call glean(thy(xco%o%xct))
          call glean(thy(xco%o%lay))
          call glean(thy(xco%o%sg))
          deallocate( xco%o )
        end if
      end subroutine

      subroutine bequeath_xco(xco)
!doc$ subroutine bequeath(xco)
        type(xc_orbital_obj) :: xco

!cod$
        continue
      end subroutine

      subroutine assign_xco(xco,xco2)
!doc$ subroutine assign(xco,xco2)
        type(xc_orbital_obj), intent(inout) :: xco
        type(xc_orbital_obj), intent(in) :: xco2

!cod$
        type(xc_orbital_obj) :: xcot
        call my(xco2)
        xcot%o => xco%o
        xco%o%ref = xco%o%ref - xco%ref
        xco%o => xco2%o
        xco%o%ref = xco%o%ref + xco%ref
        call glean(xcot)
        call glean(thy(xco2))
      end subroutine

      function xco_ref(xco) result(r)
!doc$ function x_ref(xco) result(r)
        type(xc_orbital_obj) :: xco
        integer, dimension(2) :: r
!       effects: Returns xco%ref and xco%o%ref.

!cod$
        r(1) = xco%ref
        r(2) = xco%o%ref
        call glean(xco)
      end function

      function xco_ghost(xco) result(g)
!doc$ function x_ghost(xco) result(g)
        type(xc_orbital_obj) :: xco
        type(ghost) :: g
!       effects: Returns the ghost of xco.

!cod$
        call my(xco)
        g = xco%o%g
        call glean(thy(xco))
      end function

      function xco_xc_type(xco) result(xct)
!doc$ function x_xc_type(xco) result(xct)
        type(xc_orbital_obj) :: xco
        type(xc_type_obj) :: xct
!       effects: Returns the xc_type of xco.

!cod$
        call my(xco)
        call my(xco%o%xct,xct)
        call glean(thy(xco))
        call bequeath(thy(xct))
      end function

      function xco_layout(xco) result(lay)
!doc$ function x_layout(xco) result(lay)
        type(xc_orbital_obj) :: xco
        type(layout_obj) :: lay
!       effects: Returns the lay of xco.

!cod$
        call my(xco)
        call my(xco%o%lay,lay)
        call glean(thy(xco))
        call bequeath(thy(lay))
      end function

      function xco_space_group(xco) result(sg)
!doc$ function x_space_group(xco) result(sg)
        type(xc_orbital_obj) :: xco
        type(space_group_obj) :: sg
!       effects: Returns the sg of xco.

!cod$
        call my(xco)
        call my(xco%o%sg,sg)
        call glean(thy(xco))
        call bequeath(thy(sg))
      end function

      subroutine exx_energy_xco(xco,el,e)
!doc$ subroutine exx_energy(xco,el,e)
        type(xc_orbital_obj) :: xco
        type(electrons_sc_obj) :: el
        real(double) :: e
!       effects: Returns the exx energy.
!       errors: Passes errors.

!cod$
        integer :: ik1, ik2, nb, nk
        real(double) :: afi, afv, cell_volume, e1, pfe, rc
        real(double), dimension(:), allocatable :: wt1, wt2
        type(electrons_sc_rep), pointer :: elr
        type(multivector_obj) :: mv1, mv2

        call my(xco)
        call my(el)

        elr => wormhole(el)

        nk = size(elr%occs,1)
        nb = size(elr%occs,2)

        cell_volume = x_cell_volume(x_lattice(x_layout(elr%hc)))

        allocate( wt1(nb), wt2(nb) )

        e1 = 0.0_double

        ! auxiliary function values
        select case (xco%o%coulomb_kernel)
        case (CK_NORMAL)
          afi = aux_function_integral_i(xco%o)                             ; if (error()) goto 300
          select case (mpi_nsgroups())
          case (1)
            pfe = two_pi/cell_volume
          case (2)
            pfe = four_pi/cell_volume
          end select
        end select

        ! attenuation radius
        rc = 0.0_double
        select case (xco%o%coulomb_kernel)
        case (CK_ATTENUATED)
          rc = ((3.0_double/four_pi)*real(nk,double)*cell_volume)**(1.0_double/3.0_double)
        end select

        do ik1 = 1,nk
          wt1 = x_kweight(elr%kpoints,ik1)*elr%occs(ik1,:)
          call my(x_multivector(elr%wf(ik1)),mv1)                          ; if (error()) goto 300
          call exx_energy(mv1,wt1,xco%o%xct,rc,e1)                         ; if (error()) goto 200
          select case (xco%o%coulomb_kernel)
          case (CK_NORMAL)
            e1 = e1 - pfe*afi*sum(elr%occs(ik1,:)*wt1)
          end select
          do ik2 = 1,nk
            if (ik2 == ik1) cycle
            wt2 = x_kweight(elr%kpoints,ik2)*elr%occs(ik2,:)
            call my(x_multivector(elr%wf(ik2)),mv2)                        ; if (error()) goto 200
            call exx_energy(mv1,mv2,wt1,wt2,xco%o%xct,rc,e1)               ; if (error()) goto 100
            select case (xco%o%coulomb_kernel)
            case (CK_NORMAL)
              afv = aux_function_value_i(xco%o,elr,ik1,ik2)                ; if (error()) goto 100
              e1 = e1 + pfe*afv*sum(wt1*wt1)
            end select
100         call glean(thy(mv2))                                           ; if (error()) exit
          end do
200       call glean(thy(mv1))                                             ; if (error()) exit
        end do

        call xcomm_allreduce(XSGROUP,MPI_SUM,e1,e)

300     if (allocated( wt1 )) deallocate( wt1 )
        if (allocated( wt2 )) deallocate( wt2 )

        nullify( elr )

        call glean(thy(xco))
        call glean(thy(el))

        if (error("Exit xc_orbital_mod::exx_energy_xco")) continue

      end subroutine

      subroutine exx_derivative_xco(xco,el,ik1,mvo)
!doc$ subroutine exx_derivative(xco,el,ik1,mvo)
        type(xc_orbital_obj) :: xco
        type(electrons_sc_obj) :: el
        integer, intent(in) :: ik1
        type(multivector_obj) :: mvo
!       requires: ik1 be a valid k-point index in el.
!       effects: Augments mvo with ik1 contributions to the exx functional derivative.
!       errors: Passes errors.

!cod$
        integer :: ik2, nb, nk
        real(double) :: afi, afv, cell_volume, pfm, rc
        real(double), dimension(:), allocatable :: rv1, rvo, wt1, wt2
        type(electrons_sc_rep), pointer :: elr
        type(multivector_obj) :: mv1, mv2

        call my(xco)
        call my(el)
        call my(mvo)

        elr => wormhole(el)

        nk = size(elr%occs,1)
        nb = size(elr%occs,2)

        cell_volume = x_cell_volume(x_lattice(x_layout(elr%hc)))

        allocate( wt1(nb), wt2(nb) )

        ! auxiliary function values
        select case (xco%o%coulomb_kernel)
        case (CK_NORMAL)
          allocate( rv1(nb), rvo(nb) )
          afi = aux_function_integral_i(xco%o)                    ; if (error()) goto 300
          select case (mpi_nsgroups())
          case (1)
            pfm = four_pi/cell_volume
          case (2)
            pfm = eight_pi/cell_volume
          end select
          rvo = 1.0_double
        end select

        ! attenuation radius
        rc = 0.0_double
        select case (xco%o%coulomb_kernel)
        case (CK_ATTENUATED)
          rc = ((3.0_double/four_pi)*real(nk,double)*cell_volume)**(1.0_double/3.0_double)
        end select

        wt1 = x_kweight(elr%kpoints,ik1)*elr%occs(ik1,:)
        call my(x_multivector(elr%wf(ik1)),mv1)                   ; if (error()) goto 300
        call exx_derivative(mv1,wt1,xco%o%xct,rc,mvo)             ; if (error()) goto 200
        select case (xco%o%coulomb_kernel)
        case (CK_NORMAL)
          rv1 = -pfm*afi*elr%occs(ik1,:)
          call combine(rvo,mvo,rv1,mv1)
        end select
        do ik2 = 1,nk
          if (ik2 == ik1) cycle
          wt2 = x_kweight(elr%kpoints,ik2)*elr%occs(ik2,:)
          call my(x_multivector(elr%wf(ik2)),mv2)                 ; if (error()) goto 200
          call exx_derivative(mv1,mv2,wt2,xco%o%xct,rc,mvo)       ; if (error()) goto 100
          select case (xco%o%coulomb_kernel)
          case (CK_NORMAL)
            afv = aux_function_value_i(xco%o,elr,ik1,ik2)         ; if (error()) goto 100
            rv1 = pfm*afv*wt1
            call combine(rvo,mvo,rv1,mv1)
          end select
100       call glean(thy(mv2))                                    ; if (error()) exit
        end do

200      call glean(thy(mv1))

300     if (allocated( rv1 )) deallocate( rv1 )
        if (allocated( rvo )) deallocate( rvo )
        if (allocated( wt1 )) deallocate( wt1 )
        if (allocated( wt2 )) deallocate( wt2 )

        nullify( elr )

        call glean(thy(xco))
        call glean(thy(el))
        call glean(thy(mvo))

        if (error("Exit xco::exx_derivative_xco")) continue

      end subroutine

      subroutine exx_energy_and_derivative_xco(xco,el,ik1,e,mvo)
!doc$ subroutine exx_energy_and_derivative(xco,el,ik1,e,mvo)
        type(xc_orbital_obj) :: xco
        type(electrons_sc_obj) :: el
        integer, intent(in) :: ik1
        real(double), intent(inout) :: e
        type(multivector_obj) :: mvo
!       requires: ik1 be a valid k-point index in el.
!       effects: Augments e and mvo with ik1 contributions to the exx functional and its derivative.
!       errors: Passes errors.

!cod$
        integer :: ik2, nb, nk
        real(double) :: afi, afv, cell_volume, e1, pfe, pfm, rc
        real(double), dimension(:), allocatable :: rv1, rvo, wt1, wt2
        type(electrons_sc_rep), pointer :: elr
        type(multivector_obj) :: mv1, mv2

        call my(xco)
        call my(el)
        call my(mvo)

        elr => wormhole(el)

        nk = size(elr%occs,1)
        nb = size(elr%occs,2)

        cell_volume = x_cell_volume(x_lattice(x_layout(elr%hc)))

        allocate( wt1(nb), wt2(nb) )

        e1 = 0.0_double

        ! auxiliary function values
        select case (xco%o%coulomb_kernel)
        case (CK_NORMAL)
          allocate( rv1(nb), rvo(nb) )
          select case (mpi_nsgroups())
          case (1)
            pfe = two_pi/cell_volume
            pfm = four_pi/cell_volume
          case (2)
            pfe = four_pi/cell_volume
            pfm = eight_pi/cell_volume
          end select
          rvo = 1.0_double
        end select

        ! attenuation radius
        rc = 0.0_double
        select case (xco%o%coulomb_kernel)
        case (CK_ATTENUATED)
          rc = ((3.0_double/four_pi)*real(nk,double)*cell_volume)**(1.0_double/3.0_double)
        end select

        wt1 = x_kweight(elr%kpoints,ik1)*elr%occs(ik1,:)
        call my(x_multivector(elr%wf(ik1)),mv1)                                ; if (error()) goto 300
        call exx_energy_and_derivative(mv1,wt1,xco%o%xct,rc,e1,mvo)            ; if (error()) goto 200
        select case (xco%o%coulomb_kernel)
        case (CK_NORMAL)
          afi = aux_function_integral_i(xco%o)                                 ; if (error()) goto 200
          e1 = e1 - pfe*afi*sum(elr%occs(ik1,:)*wt1)
          rv1 = -pfm*afi*elr%occs(ik1,:)
          call combine(rvo,mvo,rv1,mv1)
        end select
        do ik2 = 1,nk
          if (ik2 == ik1) cycle
          wt2 = x_kweight(elr%kpoints,ik2)*elr%occs(ik2,:)
          call my(x_multivector(elr%wf(ik2)),mv2)                              ; if (error()) goto 200
          call exx_energy_and_derivative(mv1,mv2,wt1,wt2,xco%o%xct,rc,e1,mvo)  ; if (error()) goto 100
          select case (xco%o%coulomb_kernel)
          case (CK_NORMAL)
            afv = aux_function_value_i(xco%o,elr,ik1,ik2)                      ; if (error()) goto 100
            e1 = e1 + pfe*afv*sum(wt1*wt1)
            rv1 = pfm*afv*wt1
            call combine(rvo,mvo,rv1,mv1)
          end select
100       call glean(thy(mv2))                                                 ; if (error()) exit
        end do

        e = e + e1

200     call glean(thy(mv1))

300     if (allocated( rv1 )) deallocate( rv1 )
        if (allocated( rvo )) deallocate( rvo )
        if (allocated( wt1 )) deallocate( wt1 )
        if (allocated( wt2 )) deallocate( wt2 )

        nullify( elr )

        call glean(thy(xco))
        call glean(thy(el))
        call glean(thy(mvo))

        if (error("Exit xc_orbital_mod::exx_energy_and_derivative_xco")) continue

      end subroutine

! private routines

      subroutine own_i(xco)
        type(xc_orbital_obj) :: xco
        type(xc_orbital_obj) :: xcot
        if (xco%ref < xco%o%ref) then
          allocate( xcot%o )
          xcot%o%ref            = 0
          xcot%o%g              = xco%o%g
          xcot%o%coulomb_kernel = xco%o%coulomb_kernel
          xcot%o%auxiliary_type = xco%o%auxiliary_type
          xcot%o%sd_aux_form    = xco%o%sd_aux_form
          xcot%o%si_aux_form    = xco%o%si_aux_form
          call my(xco%o%xct,xcot%o%xct)
          call my(xco%o%lay,xcot%o%lay)
          call my(xco%o%sg,xcot%o%sg)
          xco%o%ref = xco%o%ref - xco%ref
          xco%o => xcot%o
          xco%o%ref = xco%o%ref + xco%ref
        end if
      end subroutine

      function aux_function_integral_i(xcor) result(afi)
        type(xc_orbital_rep) :: xcor
        real(double) :: afi
!       effects: Returns the analytic integral of the exx auxiliary function.

!        real(double) :: alpha
        real(double) :: a, cv
        type(lattice_obj) :: lat

        call my(x_lattice(xcor%lay),lat)

        select case (xcor%auxiliary_type)
        case (AT_LEGACY)
          afi = 0.0_double
        case (AT_STRUCTURE_DEPENDENT)  ! Wenzien, Cappellini, Bechstedt: PRB 51, 14701 (1995)
          cv = x_cell_volume(lat)
          select case(xcor%sd_aux_form)
          case(SDF_SC)
            a = cv**(1.0_double/3.0_double)
            afi = 9.9774204_double*(a/two_pi)**2
          case(SDF_BCC)
            a = (2.0_double*cv)**(1.0_double/3.0_double)
            afi = 6.7172039_double*(a/two_pi)**2
          case(SDF_FCC)
            a = (4.0_double*cv)**(1.0_double/3.0_double)
            afi = 4.4237580_double*(a/two_pi)**2
          end select
!        case (AT_STRUCTURE_INDEPENDENT)
!          cv = x_cell_volume(lat)
!          select case (xcor%si_aux_form)
!          case (SIF_SFH)  ! Sorouri, Foulkes, Hine: JCP 124, 064105 (2006)
!            A value of alpha must be provided here
!            afi = (cv/two_pi**3)*two_pi*sqrt(pi/alpha)  ! I don't know if the cv/two_pi**3 factor should be here??
!          case (SIF_CRG)  ! Carrier, Rohra, Gorling: PRB 75. 205126 (2007)
!            if (error(.true.,"ERROR: crg exx_si_aux_form is not yet available")) goto 100
!          end select
        end select

100     call glean(thy(lat))

        if (error("Exit xc_orbital_mod::aux_function_integral_i")) continue

      end function

      function aux_function_value_i(xcor,elr,ik1,ik2) result(afv)
        type(xc_orbital_rep) :: xcor
        type(electrons_sc_rep) :: elr
        integer, intent(in) :: ik1, ik2
        real(double) :: afv
!       effects: Returns the value of the exx auxiliary function at the difference of k-points corresponding to ik1 and ik2.

!        integer :: i1, i2, i3
!        real(double) :: afv_local, alpha, kg2
        real(double) :: a, cv, dkx, dky, dkz
        real(double), dimension(3) :: kpt1, kpt2
        real(double), dimension(:,:,:), pointer :: g2, gx, gy, gz
        type(lattice_obj) :: lat

        nullify( g2 )
        nullify( gx )
        nullify( gy )
        nullify( gz )

        call my(x_lattice(xcor%lay),lat)

        kpt1 = lat2f(lat,x_kpoint(elr%kpoints,ik1))
        kpt2 = lat2f(lat,x_kpoint(elr%kpoints,ik2))

        dkx = kpt1(1) - kpt2(1)
        dky = kpt1(2) - kpt2(2)
        dkz = kpt1(3) - kpt2(3)

        select case (xcor%auxiliary_type)
        case (AT_LEGACY)
          afv = 0.0_double
        case (AT_STRUCTURE_DEPENDENT)  ! Wenzien, Cappellini, Bechstedt: PRB 51, 14701 (1995)
          cv = x_cell_volume(lat)
          select case(xcor%sd_aux_form)
          case(SDF_SC)
            a = cv**(1.0_double/3.0_double)
            afv = (a**2/2.0_double)/(3.0_double - cos(a*dkx) - cos(a*dky) - cos(a*dkz))
          case(SDF_BCC)
            a = (2.0_double*cv)**(1.0_double/3.0_double)
            afv = (a**2/8.0_double)/(1.0_double - cos(a*dkx/2.0_double) - cos(a*dky/2.0_double) - cos(a*dkz/2.0_double))
          case(SDF_FCC)
            a = (4.0_double*cv)**(1.0_double/3.0_double)
            afv = (a**2/4.0_double)/(3.0_double - cos(a*dkx/2.0_double)*cos(a*dky/2.0_double) &
                                                - cos(a*dkx/2.0_double)*cos(a*dkz/2.0_double) &
                                                - cos(a*dky/2.0_double)*cos(a*dkz/2.0_double))
          end select
!        case (AT_STRUCTURE_INDEPENDENT)
!          cv = x_cell_volume(lat)
!          select case (xcor%si_aux_form)
!          case (SIF_SFH)  ! Sorouri, Foulkes, Hine: JCP 124, 064105 (2006)
!            A value of alpha must be provided here
!            call fdel(g2,xcor%lay,D_TYPE,CONFIG)
!            call fmesh(gx,gy,gz,xcor%lay,D_TYPE,CONFIG)
!            afv_local = 0.0_double
!            do i3 = 1,size(g2,3)
!            do i2 = 1,size(g2,2)
!            do i1 = 1,size(g2,1)
!              if (g2(i1,i2,i3) > x_cutoff(xcor%lay)) cycle
!              kg2 = (dkx - gx(i1,i2,i3))**2 + (dky - gy(i1,i2,i3))**2 + (dkz - gz(i1,i2,i3))**2
!              afv_local = afv_local + exp(-alpha*kg2)/kg2
!            end do
!            end do
!            end do
!            call allreduce(CONFIG,MPI_SUM,afv_local,afv)
!          case (SIF_CRG)  ! Carrier, Rohra, Gorling: PRB 75. 205126 (2007)
!            if (error(.true.,"ERROR: crg exx_si_aux_form is not yet available")) goto 100
!          end select
        end select

100     if (associated( g2 )) deallocate( g2 )
        if (associated( gx )) deallocate( gx )
        if (associated( gy )) deallocate( gy )
        if (associated( gz )) deallocate( gz )

        call glean(thy(lat))

        if (error("Exit xc_orbital_mod::aux_function_value_i")) continue

      end function

      end module

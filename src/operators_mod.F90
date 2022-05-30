! Copyright 2011 Sandia Corporation. 

! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

      module operators_mod
!doc$ module operators_mod

!     Three datatypes are available here: type(h_common_obj)
!                                         type(h_kpoint_obj)
!                                         type(hamiltonian_obj)

!     operators_mod generates or processes components of operators and applies them to multivectors.

      use kind_mod
      use mpi_mod
      use arg_mod
      use error_mod
      use io_mod
      use diary_mod
      use ghost_mod
      use layout_mod
      use grid_mod
      use multibasis_mod
      use math_mod
      use fft_mod
      use multivector_mod
      use external_mod
      use crystal_mod
      use lattice_mod
      use atoms_mod
      use gen_density_mod
      use gen_potential_mod
      use point_blas_mod
      use atomic_operators_mod
      use atomic_density_mod
      use atomic_potential_mod
      use dyad_kpoint_mod
      use timing_mod

!cod$
      implicit none
      private
   
      integer, parameter :: FOURIER_SPACE = 1
      integer, parameter :: REAL_SPACE    = 2

      type :: hc_replica_data
        integer :: sr_index                                      ! index into summing (srdots and gsrdots) arrays
        real(double), dimension(3) :: pos                        ! replica position (cartesian)
      end type

      type :: hc_site_data
        logical :: new                                           ! indicates if the projector is at a new site
        real(double), dimension(3) :: pos                        ! projector position in the {0,0,0} parallelpiped (cartesian)
        type(hc_replica_data), dimension(:), pointer :: replica  ! replica information for this projector
      end type

      type :: h_common_rep                   ! hamiltonian information common to all k-points
        integer :: ref
        type(ghost) :: g 
        type(ghost) :: g_gpot
        type(crystal_obj) :: crystal                           ! crystal
        type(layout_obj) :: layout                             ! layout
        type(atomic_operators_obj) :: ao                       ! atomic operators
        type(atomic_potential_obj) :: apot                     ! atomic potential
        real(double), dimension(:,:,:), pointer :: lp          ! local potential
        integer :: projector_type                              ! type of non-local projectors
        real(double) :: cutoff                                 ! cutoff
        integer :: pbs                                         ! pdots blocksize (Fourier-space representation)
        type(hc_site_data), dimension(:), pointer :: site      ! real-space projector site data
        integer, dimension(:), pointer :: rsi                  ! real-space projector indexes
        integer, dimension(:,:,:), pointer :: rsc1, rsc2       ! real-space projector counters
        real(double), dimension(:), pointer :: rsp             ! real-space projectors
        real(double), dimension(:,:), pointer :: grsp          ! real-space projector gradients
      end type
   
      type, public :: h_common_obj
        private         
        integer :: ref
        type(h_common_rep), pointer :: o
      end type

      type :: hk_type_info
        real(double), dimension(:,:), pointer :: form_factor         ! form factors
      end type
   
      type :: hk_replica_data
        complex(double) :: phase                                     ! phase factor for the parallelpiped containing the replica
      end type

      type :: hk_site_data
        type(hk_replica_data), dimension(:), pointer :: replica      ! replica information for this site
      end type

      type :: h_kpoint_rep             ! hamiltonian information for a particular k-point
        integer :: ref
        type(ghost) :: g 
        type(multibasis_obj) :: mb                         ! multibasis
        type(h_common_obj) :: hc                           ! common part of the full hamiltonian
        type(dyad_kpoint_obj), pointer :: dk               ! dyadic potential for this k-point
        real(double), dimension(:), pointer :: ke          ! kinetic energy operator
        logical :: need_gk                                 ! indicates whether the (g + k) vector is needed
        real(double), dimension(:), pointer :: gkx         ! x component of the (g + k) vector
        real(double), dimension(:), pointer :: gky         ! y component of the (g + k) vector
        real(double), dimension(:), pointer :: gkz         ! z component of the (g + k) vector
        type(hk_type_info), dimension(:), pointer :: type  ! type-specific information for reciprocal-space projectors
        type(hk_site_data), dimension(:), pointer :: site  ! real-space projector site data
      end type
   
      type, public :: h_kpoint_obj
        private
        integer :: ref
        type(h_kpoint_rep), pointer :: o
      end type

      type :: hamiltonian_rep           ! full hamiltonian
        integer :: ref
        type(ghost) :: g
        type(h_kpoint_obj) :: hk
        complex(double), dimension(:,:), pointer :: vnl      ! expanded reciprocal-space projectors
        complex(double), dimension(:,:,:), pointer :: phase  ! phase factors used with the real-space projectors
      end type

      type, public :: hamiltonian_obj
        private
        integer :: ref
        type(hamiltonian_rep), pointer :: o
      end type
   
      type :: point_data
        integer, dimension(:), pointer :: i
        real(double), dimension(:), pointer :: p
        real(double), dimension(:,:), pointer :: pg
      end type

      type :: type_info
        real(double), dimension(:,:), pointer :: pff
        real(double), dimension(:,:,:), pointer :: sff
      end type

      type :: atom_info
        complex(double), dimension(:), pointer :: sf
      end type

!doc$
      public :: common_hamiltonian
      public :: kpoint_hamiltonian
      public :: hamiltonian
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_crystal
      public :: x_layout
      public :: x_cutoff
      public :: x_h_common
      public :: x_multibasis
      public :: lp_minimum
      public :: lp_maximum
      public :: apply_hamiltonian
      public :: apply_overlap
      public :: add_density
      public :: orthonormalize
      public :: precondition
      public :: kinetic_energy
      public :: forces
      public :: pressure
      public :: stress_tensor
      public :: overlap_is_identity
      public :: diary

!cod$

      interface common_hamiltonian
        module procedure constructor_hc
      end interface
      interface kpoint_hamiltonian
        module procedure constructor_hk
      end interface
      interface hamiltonian
        module procedure constructor_h
      end interface
      interface update
        module procedure update_hc, update_hk
      end interface
      interface my
        module procedure my_hc, my_new_hc, my_hk, my_new_hk, my_h, my_new_h
      end interface
      interface thy
        module procedure thy_hc, thy_hk, thy_h
      end interface
      interface glean
        module procedure glean_hc, glean_hk, glean_h
      end interface
      interface bequeath
        module procedure bequeath_hc, bequeath_hk, bequeath_h
      end interface
      interface assignment(=)
        module procedure assign_hc, assign_hk, assign_h
      end interface
      interface x_ref
        module procedure hc_ref, hk_ref, h_ref
      end interface
      interface x_ghost
        module procedure hc_ghost, hk_ghost, h_ghost
      end interface
      interface x_crystal
        module procedure hc_crystal
      end interface
      interface x_layout
        module procedure hc_layout
      end interface
      interface x_cutoff
        module procedure hc_cutoff, hk_cutoff, h_cutoff
      end interface
      interface x_h_common
        module procedure hk_h_common, h_h_common
      end interface
      interface x_multibasis
        module procedure hk_multibasis
      end interface
      interface lp_minimum
        module procedure lp_minimum_hc
      end interface
      interface lp_maximum
        module procedure lp_maximum_hc
      end interface
      interface apply_hamiltonian
        module procedure apply_hamiltonian_h
      end interface
      interface apply_overlap
        module procedure apply_overlap_h
      end interface
      interface add_density
        module procedure add_density_h
      end interface
      interface orthonormalize
        module procedure orthonormalize_h
      end interface
      interface precondition
        module procedure precondition_h
      end interface
      interface kinetic_energy
        module procedure kinetic_energy_hk
      end interface
      interface forces
        module procedure forces_h
      end interface
      interface pressure
        module procedure pressure_h
      end interface
      interface stress_tensor
        module procedure stress_tensor_h
      end interface
      interface overlap_is_identity
        module procedure overlap_is_identity_hc, overlap_is_identity_h
      end interface
      interface diary
        module procedure diary_hc
      end interface
   
      contains

! public hc routines

      function constructor_hc(ext,gpot,cutoff,nb) result(hc)
!doc$ function common_hamiltonian(ext,gpot,cutoff,nb) result(hc)
        type(external_obj) :: ext
        type(gen_potential_obj) :: gpot
        real(double), intent(in) :: cutoff
        integer, intent(in) :: nb
        type(h_common_obj) :: hc
!       effects: Creates a new hc.
!       errors: Projector_type not recognized. Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag
        type(grid_obj) :: lpg

        if (error("  Error on entry")) then
          hc%ref = 0
          allocate( hc%o )
          hc%o%ref = 0
          goto 999
        end if

        call my(ext)
        call my(gpot)

        hc%ref = 0
        allocate( hc%o )
        hc%o%ref = 0
        hc%o%g = x_ghost()

        hc%o%g_gpot = x_ghost(gpot)

        call my(x_crystal(ext),hc%o%crystal)
        call my(x_layout(ext),hc%o%layout)
        call my(x_atomic_operators(ext),hc%o%ao)
        call my(x_atomic_potential(gpot),hc%o%apot)

        call my(x_grid_potential(gpot),lpg) ! EXTRACT THE LOCAL POTENTIAL
        call take(hc%o%lp,lpg,RS_KIND)
        call glean(thy(lpg))

        call arglc("projector_type",tag,found) ! DETERMINE THE TYPE OF PROJECTORS AND DIARY THAT INFORMATION
        if (.not.found) tag = "reciprocal"
        select case (trim(tag))
        case ("reciprocal")
          hc%o%projector_type = FOURIER_SPACE
        case ("real")
          hc%o%projector_type = REAL_SPACE
        case default
          if (error(.true.,"ERROR: projector_type not recognized")) goto 100
        end select

        hc%o%cutoff = cutoff

        select case (hc%o%projector_type)
        case (FOURIER_SPACE)
          call arg("pdots_blocksize",hc%o%pbs,found)
          if (.not.found) hc%o%pbs = nb
          if (error(hc%o%pbs < 1,"ERROR: hc%o%pbs < 1")) goto 100
          if (error(hc%o%pbs > nb,"ERROR: hc%o%pbs > nb")) goto 100
        case (REAL_SPACE)
          hc%o%pbs = 0
        end select

        nullify( hc%o%site, hc%o%rsi, hc%o%rsc1, hc%o%rsc2, hc%o%rsp, hc%o%grsp )
        select case (hc%o%projector_type) ! FORM THE REAL-SPACE PROJECTORS
        case (REAL_SPACE)
          call form_rs_projectors_i(hc%o) ; if (error()) goto 100
        end select

100     call glean(thy(ext))
        call glean(thy(gpot))

999     if (error("Exit operators_mod::constructor_hc")) continue

      end function

      subroutine update_hc(hc,ext,gpot,cutoff)
!doc$ subroutine update(hc,ext,gpot,cutoff)
        type(h_common_obj) :: hc
        type(external_obj) :: ext
        type(gen_potential_obj) :: gpot
        real(double), intent(in), optional :: cutoff
!       modifies: hc
!       effects: Updates hc.
!       errors: Passes errors.

!cod$
        logical :: cutoff_change, crystal_change, layout_change, ao_change, gpot_change, any_change
        type(grid_obj) :: lpg

        if (error("  Error on entry")) goto 999

        call my(hc)
        call my(ext)
        call my(gpot)

        if (present(cutoff)) then
          cutoff_change = (cutoff /= hc%o%cutoff)
          if (error(cutoff_change,"ERROR: cutoff changes are not currently allowed")) goto 100
        end if

        crystal_change = (x_ghost(x_crystal(ext)) /= x_ghost(hc%o%crystal))

        layout_change = (x_ghost(x_layout(ext)) /= x_ghost(hc%o%layout))

        ao_change = (x_ghost(x_atomic_operators(ext)) /= x_ghost(hc%o%ao))

        gpot_change = (x_ghost(gpot) /= hc%o%g_gpot)

        any_change = (crystal_change .or. layout_change .or. ao_change .or. gpot_change)

        if (any_change) then
          call own_hc_i(hc)
          hc%o%g = x_ghost()
          if (gpot_change) then
            hc%o%g_gpot = x_ghost(gpot)
            hc%o%apot = x_atomic_potential(gpot)
            call my(x_grid_potential(gpot),lpg)
            deallocate( hc%o%lp )
            call take(hc%o%lp,lpg,RS_KIND)
            call glean(thy(lpg))
          end if
          if (crystal_change) hc%o%crystal = x_crystal(ext)
          if (layout_change) hc%o%layout = x_layout(ext)
          if (ao_change) hc%o%ao = x_atomic_operators(ext)
          if (ao_change .or. crystal_change .or. layout_change) then
            select case (hc%o%projector_type)
            case (REAL_SPACE)
              call form_rs_projectors_i(hc%o) ; if (error()) goto 100
            end select
          end if
        end if

100     call glean(thy(hc))
        call glean(thy(ext))
        call glean(thy(gpot))

999     if (error("Exit operators_mod::update_hc")) continue

      end subroutine

      subroutine my_hc(hc)
!doc$ subroutine my(hc)
        type(h_common_obj) :: hc

!cod$
        hc%ref = hc%ref + 1
        hc%o%ref = hc%o%ref + 1
      end subroutine

      subroutine my_new_hc(hci,hc)
!doc$ subroutine my(hci,hc)
        type(h_common_obj) :: hci, hc

!cod$
        hc%ref = 1
        hc%o => hci%o
        hc%o%ref = hc%o%ref + 1
      end subroutine

      function thy_hc(hc) result(hco)
!doc$ function thy(hc) result(hco)
        type(h_common_obj) :: hc, hco

!cod$
        hc%ref = hc%ref - 1
        hc%o%ref = hc%o%ref - 1
        hco%ref = hc%ref
        hco%o => hc%o
      end function
   
      subroutine glean_hc(hc)
!doc$ subroutine glean(hc)
        type(h_common_obj) :: hc

!cod$
        integer :: ip
        if (hc%o%ref < 1) then
          call glean(thy(hc%o%crystal))
          call glean(thy(hc%o%layout))
          call glean(thy(hc%o%ao))
          call glean(thy(hc%o%apot))
          if (associated( hc%o%lp )) deallocate( hc%o%lp )
          if (associated( hc%o%site )) then
            do ip = 1,size(hc%o%site)
              if (associated( hc%o%site(ip)%replica )) deallocate( hc%o%site(ip)%replica )
            end do
            deallocate( hc%o%site )
          end if
          if (associated( hc%o%rsi )) deallocate( hc%o%rsi )
          if (associated( hc%o%rsc1 )) deallocate( hc%o%rsc1 )
          if (associated( hc%o%rsc2 )) deallocate( hc%o%rsc2 )
          if (associated( hc%o%rsp )) deallocate( hc%o%rsp )
          if (associated( hc%o%grsp )) deallocate( hc%o%grsp )
          deallocate( hc%o )
        end if
      end subroutine

      subroutine bequeath_hc(hc)
!doc$ subroutine bequeath(hc)
        type(h_common_obj) :: hc

!cod$
        continue
      end subroutine

      subroutine assign_hc(hc,hc2)
!doc$ subroutine assignment(=)(hc,hc2)
        type(h_common_obj), intent(inout) :: hc
        type(h_common_obj), intent(in) :: hc2

!cod$
        type(h_common_obj) :: hct
        call my(hc2)
        hct%o => hc%o
        hc%o%ref = hc%o%ref - hc%ref
        hc%o => hc2%o
        hc%o%ref = hc%o%ref + hc%ref
        call glean(hct)
        call glean(thy(hc2))
      end subroutine
 
      function hc_ref(hc) result(r)
!doc$ function x_ref(hc) result(r)
        type(h_common_obj) :: hc
        integer, dimension(2) :: r
!       effects: Returns hc%ref and hc%o%ref.

!cod$
        r(1) = hc%ref
        r(2) = hc%o%ref
        call glean(hc)
      end function

      function hc_ghost(hc) result(g)
!doc$ function x_ghost(hc) result(g)
        type(h_common_obj) :: hc
        type(ghost) :: g
!       effects: Returns hc%o%g.

!cod$
        call my(hc)
        g = hc%o%g
        call glean(thy(hc))
      end function

      function hc_crystal(hc) result(crys)
!doc$ function x_crystal(hc) result(crys)
        type(h_common_obj) :: hc
        type(crystal_obj) :: crys
!       effects: Returns hc%o%crystal.

!cod$
        call my(hc)
        call my(hc%o%crystal,crys)
        call bequeath(thy(crys))
        call glean(thy(hc))
      end function

      function hc_layout(hc) result(lay)
!doc$ function x_layout(hc) result(lay)
        type(h_common_obj) :: hc
        type(layout_obj) :: lay
!       effects: Returns hc%o%layout.

!cod$
        call my(hc)
        call my(hc%o%layout,lay)
        call bequeath(thy(lay))
        call glean(thy(hc))
      end function

      function hc_cutoff(hc) result(c)
!doc$ function x_cutoff(hc) result(c)
        type(h_common_obj) :: hc
        real(double) :: c
!       effects: Returns hc%o%cutoff.

!cod$
        call my(hc)
        c = hc%o%cutoff
        call glean(thy(hc))
      end function

      function lp_minimum_hc(hc) result(m)
!doc$ function lp_minimum(hc) result(m)
        type(h_common_obj) :: hc
        real(double) :: m
!       effects: Returns the minimum value of hc%o%lp.
!cod$
        call my(hc)
        m = minval(hc%o%lp)
        call glean(thy(hc))
      end function

      function lp_maximum_hc(hc) result(m)
!doc$ function lp_maximum(hc) result(m)
        type(h_common_obj) :: hc
        real(double) :: m
!       effects: Returns the maximum value of hc%o%lp.
!cod$
        call my(hc)
        m = maxval(hc%o%lp)
        call glean(thy(hc))
      end function

      function overlap_is_identity_hc(hc) result(l)
!doc$ function overlap_is_identity(hc) result(l)
        type(h_common_obj) :: hc
        logical :: l
!       effects: Returns true if the overlap operator is the identity.

!cod$
        call my(hc)
        l = overlap_is_identity(hc%o%apot)
        call glean(thy(hc))
      end function

      subroutine diary_hc(hc)
!doc$ subroutine diary(hc)
        type(h_common_obj) :: hc
!       effects: Writes hc information to the diary.

!cod$
        call my(hc)

        if (x_n_projectors(hc%o%ao) /= 0) then
          select case (hc%o%projector_type)
          case (FOURIER_SPACE)
            if (i_access( diaryfile() )) write(x_unit(diaryfile()),'(/,t4,"Reciprocal-space projectors")')
            if (i_access( diaryfile() )) write(x_unit(diaryfile()),'(t6,"using pdots blocksize = ",i0)') hc%o%pbs
          case (REAL_SPACE)
            if (i_access( diaryfile() )) write(x_unit(diaryfile()),'(/,t4,"Real-space projectors:")')
            call diary_rs_projectors(hc%o%ao)
          end select
        end if

        call glean(thy(hc))

      end subroutine

! public hk routines

      function constructor_hk(hc,mb,dk) result(hk)
!doc$ function kpoint_hamiltonian(hc,dk,mb) result(hk)
        type(h_common_obj) :: hc
        type(multibasis_obj) :: mb
        type(dyad_kpoint_obj), optional :: dk
        type(h_kpoint_obj) :: hk
!       effects: Creates a new hk.
!       errors: Passes errors.

!cod$
        if (error("  Error on entry")) then
          hk%ref = 0
          allocate( hk%o )
          hk%o%ref = 0
          goto 999
        end if

        call my(hc)
        call my(mb)
        if (present(dk)) call my(dk)

        hk%ref = 0
        allocate( hk%o )
        hk%o%ref = 0
        hk%o%g = x_ghost()

        call my(hc,hk%o%hc)
        if (present(dk)) then
           allocate(hk%o%dk)
           call my(dk,hk%o%dk)
        else
           nullify(hk%o%dk)
        end if
        call my(mb,hk%o%mb)

        hk%o%need_gk = .false. ! FORM THE KINETIC-ENERGY OPERATOR
        nullify( hk%o%ke, hk%o%gkx, hk%o%gky, hk%o%gkz )
        call form_kinetic_energy_i(hk%o) ; if (error()) goto 100

        nullify( hk%o%type, hk%o%site )
        select case (hk%o%hc%o%projector_type) ! FORM THE FOURIER-SPACE PROJECTORS OR THE REAL-SPACE PROJECTOR REPLICA PHASES
        case (FOURIER_SPACE)
          call form_fs_type_projectors_i(hk%o) ; if (error()) goto 100
        case (REAL_SPACE)
          call form_rs_replica_phases_i(hk%o) ; if (error()) goto 100
        end select

100     call glean(thy(hc))
        call glean(thy(mb))
        if (present(dk)) call glean(thy(dk))

999     if (error("Exit operators_mod::constructor_hk")) continue

      end function

      subroutine update_hk(hk,hc,mb,dk)
!doc$ subroutine update(hk,hc,mb,dk)
        type(h_kpoint_obj) :: hk
        type(h_common_obj), optional :: hc
        type(multibasis_obj), optional :: mb
        type(dyad_kpoint_obj), optional :: dk
!       modifies: hk
!       effects: Updates hk.
!       errors: Change in projector type. Passes errors.

!cod$
        logical :: hc_change, mb_change, crystal_change, layout_change, lattice_change, dk_change

        if (error("  Error on entry")) goto 999

        call my(hk)
        if (present(hc)) call my(hc)
        if (present(mb)) call my(mb)
        if (present(dk)) call my(dk)

        hc_change = .false.
        mb_change = .false.
        crystal_change = .false.
        layout_change = .false.
        lattice_change = .false.

        if (present(hc)) then
          hc_change = (x_ghost(hc) /= x_ghost(hk%o%hc))
          if (hc_change) then
            crystal_change = (x_ghost(x_crystal(hc)) /= x_ghost(x_crystal(hk%o%hc)))
            layout_change = (x_ghost(x_layout(hc)) /= x_ghost(x_layout(hk%o%hc)))
            if (crystal_change) lattice_change = (x_ghost(x_lattice(x_crystal(hc))) /= x_ghost(x_lattice(x_crystal(hk%o%hc))))
            if (error(hc%o%projector_type /= hk%o%hc%o%projector_type,"ERROR: projector type has changed")) goto 100
          end if
        end if
        if (present(mb)) mb_change = (x_ghost(mb) /= x_ghost(hk%o%mb))
        if (present(dk)) then
          if (associated(hk%o%dk)) then
             dk_change = (x_ghost(dk) /= x_ghost(hk%o%dk))
          else
             dk_change = .true.
          end if
        else
          dk_change = .false.
        end if

        if (hc_change .or. dk_change .or. mb_change) then
          call own_hk_i(hk)
          hk%o%g = x_ghost()
          if (hc_change) hk%o%hc = hc
          if (dk_change) then
             if (associated(hk%o%dk)) then
                hk%o%dk = dk
             else
                allocate(hk%o%dk)
                call my(dk,hk%o%dk)
             end if
          end if
          if (mb_change) then
            hk%o%mb = mb
            call form_kinetic_energy_i(hk%o) ; if (error()) goto 100
          end if
          select case (hk%o%hc%o%projector_type)
          case (FOURIER_SPACE)
            if (lattice_change .or. mb_change) then
              call form_fs_type_projectors_i(hk%o) ; if (error()) goto 100
            end if
          case (REAL_SPACE)
            if (crystal_change .or. layout_change .or. mb_change) then
              call form_rs_replica_phases_i(hk%o) ; if (error()) goto 100
            end if
          end select
        end if

100     call glean(thy(hk))
        if (present(hc)) call glean(thy(hc))
        if (present(mb)) call glean(thy(mb))
        if (present(dk)) call glean(thy(dk))

999     if (error("Exit operators_mod::update_hk")) continue

      end subroutine

      subroutine my_hk(hk)
!doc$ subroutine my(hk)
        type(h_kpoint_obj) :: hk

!cod$
        hk%ref = hk%ref + 1
        hk%o%ref = hk%o%ref + 1
      end subroutine

      subroutine my_new_hk(hki,hk)
!doc$ subroutine my(hki,hk)
        type(h_kpoint_obj) :: hki, hk

!cod$
        hk%ref = 1
        hk%o => hki%o
        hk%o%ref = hk%o%ref + 1
      end subroutine

      function thy_hk(hk) result(hko)
!doc$ function thy(hk) result(hko)
        type(h_kpoint_obj) :: hk, hko

!cod$
        hk%ref = hk%ref - 1
        hk%o%ref = hk%o%ref - 1
        hko%ref = hk%ref
        hko%o => hk%o
      end function
   
      subroutine glean_hk(hk)
!doc$ subroutine glean(hk)
        type(h_kpoint_obj) :: hk

!cod$
        integer :: ip, it
        if (hk%o%ref < 1) then
          if (associated( hk%o%ke )) deallocate( hk%o%ke )
          if (associated( hk%o%gkx )) deallocate( hk%o%gkx )
          if (associated( hk%o%gky )) deallocate( hk%o%gky )
          if (associated( hk%o%gkz )) deallocate( hk%o%gkz )
          call glean(thy(hk%o%hc))
          if (associated(hk%o%dk)) then
             call glean(thy(hk%o%dk))
             deallocate(hk%o%dk)
          end if
          call glean(thy(hk%o%mb))
          if (associated( hk%o%type )) then
            do it = 1,size(hk%o%type)
              if (associated( hk%o%type(it)%form_factor )) deallocate( hk%o%type(it)%form_factor )
            end do
            deallocate( hk%o%type )
          end if
          if (associated( hk%o%site )) then
            do ip = 1,size(hk%o%site)
              if (associated( hk%o%site(ip)%replica )) deallocate( hk%o%site(ip)%replica )
            end do
            deallocate( hk%o%site )
          end if
          deallocate( hk%o )
        end if
      end subroutine

      subroutine bequeath_hk(hk)
!doc$ subroutine bequeath(hk)
        type(h_kpoint_obj) :: hk

!cod$
        continue
      end subroutine

      subroutine assign_hk(hk,hk2)
!doc$ subroutine assignment(=)(hk,hk2)
        type(h_kpoint_obj), intent(inout) :: hk
        type(h_kpoint_obj), intent(in) :: hk2

!cod$
        type(h_kpoint_obj) :: hkt
        call my(hk2)
        hkt%o => hk%o
        hk%o%ref = hk%o%ref - hk%ref
        hk%o => hk2%o
        hk%o%ref = hk%o%ref + hk%ref
        call glean(hkt)
        call glean(thy(hk2))
      end subroutine
 
      function hk_ref(hk) result(r)
!doc$ function x_ref(hk) result(r)
        type(h_kpoint_obj) :: hk
        integer, dimension(2) :: r
!       effects: Returns hk%ref and hk%o%ref.

!cod$
        r(1) = hk%ref
        r(2) = hk%o%ref
        call glean(hk)
      end function

      function hk_ghost(hk) result(g)
!doc$ function x_ghost(hk) result(g)
        type(h_kpoint_obj) :: hk
        type(ghost) :: g
!       effects: Returns hk%o%g.

!cod$
        call my(hk)
        g = hk%o%g
        call glean(thy(hk))
      end function

      function hk_cutoff(hk) result(c)
!doc$ function x_cutoff(hk) result(c)
        type(h_kpoint_obj) :: hk
        real(double) :: c
!       effects: Returns hk%o%hc%o%cutoff.

!cod$
        call my(hk)
        c = x_cutoff(hk%o%hc)
        call glean(thy(hk))
      end function

      function hk_h_common(hk) result(hc)
!doc$ function x_h_common(hk) result(hc)
        type(h_kpoint_obj) :: hk
        type(h_common_obj) :: hc
!       effects: Returns hk%o%hc.

!cod$
        call my(hk)
        call my(hk%o%hc,hc)
        call bequeath(thy(hc))
        call glean(thy(hk))
      end function

      function hk_multibasis(hk) result(mb)
!doc$ function x_multibasis(hk) result(mb)
        type(h_kpoint_obj) :: hk
        type(multibasis_obj) :: mb
!       effects: Returns the hk multibasis.

!cod$
        call my(hk)
        call my(hk%o%mb,mb)
        call bequeath(thy(mb))
        call glean(thy(hk))
      end function

      function kinetic_energy_hk(v,hk,wts) result(ke)
!doc$ function kinetic_energy(v,hk,wts) result(ke)
        type(multivector_obj) :: v
        type(h_kpoint_obj) :: hk
        real(double), dimension(:), intent(in) :: wts
        real(double) :: ke
!       requires: size(wts) = size(v%o%mat,2)
!       effects: Returns the kinetic energy of v.
!       errors: v and hk multibases are different.

!cod$
        integer :: nb, ng
        type(multivector_rep), pointer :: worm_v

        call my(v)
        call my(hk)

        worm_v => wormhole(v)

        if (error(x_ghost(hk%o%mb) /= x_ghost(worm_v%mb),"ERROR: multibases are different")) goto 100

        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)

        call kernel_kinetic_energy_i(ng,nb,hk%o%ke,worm_v%mat,wts,ke)

100     nullify( worm_v )

        call glean(thy(v))
        call glean(thy(hk))

        if (error("Exit operators_mod::kinetic_energy_hk")) continue

      end function

! public h routines

      function constructor_h(hk) result(h)
!doc$ function hamiltonian(hk) result(h)
        type(h_kpoint_obj) :: hk
        type(hamiltonian_obj) :: h
!       effects: Constructs a new hamiltonian_obj.
!       errors: Passes errors.

!cod$
        if (error("  Error on entry")) then
          h%ref = 0
          allocate( h%o )
          h%o%ref = 0
          goto 999
        end if

        call my(hk)

        h%ref = 0
        allocate( h%o )
        h%o%ref = 0
        h%o%g = x_ghost()

        call my(hk,h%o%hk)

        nullify( h%o%vnl, h%o%phase )
        select case (h%o%hk%o%hc%o%projector_type) ! EXPAND THE FOURIER-SPACE PROJECTORS OR FORM THE REAL-SPACE PROJECTOR PHASE
        case (FOURIER_SPACE)
          call expand_fs_projectors_i(h%o) ; if (error()) goto 100
        case (REAL_SPACE)
          call form_rs_global_phases_i(h%o) ; if (error()) goto 100
        end select

100     call glean(thy(hk))

999     if (error("Exit operators_mod::constructor_h")) continue

      end function

      subroutine my_h(h)
!doc$ subroutine my(h)
        type(hamiltonian_obj) :: h

!cod$
        h%ref = h%ref + 1
        h%o%ref = h%o%ref + 1
      end subroutine

      subroutine my_new_h(hi,h)
!doc$ subroutine my(hi,h)
        type(hamiltonian_obj) :: hi, h

!cod$
        h%ref = 1
        h%o => hi%o
        h%o%ref = h%o%ref + 1
      end subroutine

      function thy_h(h) result(ho)
!doc$ function thy(h) result(ho)
        type(hamiltonian_obj) :: h, ho

!cod$
        h%ref = h%ref - 1
        h%o%ref = h%o%ref - 1
        ho%ref = h%ref
        ho%o => h%o
      end function
   
      subroutine glean_h(h)
!doc$ subroutine glean(h)
        type(hamiltonian_obj) :: h

!cod$
        if (h%o%ref < 1) then
          call glean(thy(h%o%hk))
          if (associated( h%o%vnl )) deallocate( h%o%vnl )
          if (associated( h%o%phase )) deallocate( h%o%phase )
          deallocate( h%o )
        end if
      end subroutine

      subroutine bequeath_h(h)
!doc$ subroutine bequeath(h)
        type(hamiltonian_obj) :: h

!cod$
        continue
      end subroutine

      subroutine assign_h(h,h2)
!doc$ subroutine assignment(=)(h,h2)
        type(hamiltonian_obj), intent(inout) :: h
        type(hamiltonian_obj), intent(in) :: h2

!cod$
        type(hamiltonian_obj) :: ht
        call my(h2)
        ht%o => h%o
        h%o%ref = h%o%ref - h%ref
        h%o => h2%o
        h%o%ref = h%o%ref + h%ref
        call glean(ht)
        call glean(thy(h2))
      end subroutine

      function h_ref(h) result(r)
!doc$ function x_ref(h) result(r)
        type(hamiltonian_obj) :: h
        integer, dimension(2) :: r
!       effects: Returns h%ref and h%o%ref.

!cod$
        r(1) = h%ref
        r(2) = h%o%ref
        call glean(h)
      end function

      function h_ghost(h) result(g)
!doc$ function x_ghost(h) result(g)
        type(hamiltonian_obj) :: h
        type(ghost) :: g
!       effects: Returns h%o%g.

!cod$
        call my(h)
        g = h%o%g
        call glean(thy(h))
      end function

      function h_cutoff(h) result(c)
!doc$ function x_cutoff(h) result(c)
        type(hamiltonian_obj) :: h
        real(double) :: c
!       effects: Returns h%o%hk%o%hc%o%cutoff.

!cod$
        call my(h)
        c = x_cutoff(h%o%hk)
        call glean(thy(h))
      end function

      function h_h_common(h) result(hc)
!doc$ function x_h_common(h) result(hc)
        type(hamiltonian_obj) :: h
        type(h_common_obj) :: hc
!       effects: Returns h%o%hk%o%hc.

!cod$
        call my(h)
        call my(h%o%hk%o%hc,hc)
        call bequeath(thy(hc))
        call glean(thy(h))
      end function

      subroutine apply_hamiltonian_h(v,hv,h)
!doc$ subroutine apply_hamiltonian(v,hv,h)
        type(multivector_obj) :: v
        type(multivector_obj) :: hv
        type(hamiltonian_obj) :: h
!       effects: Performs the operation h*v => hv.
!       modifies: hv.
!       errors: Passes errors.

!cod$
        call start_timer("operators: apply_hamiltonian")

        call my(v)
        call my(hv)
        call my(h)

        if (associated(h%o%hk%o%dk)) then
           hv = apply(h%o%hk%o%dk,v)
        else
           call zero_mv_i(hv)
        end if
        call apply_kinetic_energy_i(v,hv,h%o%hk%o) ; if (error()) goto 100
        call apply_potentials_i(v,hv,h%o) ; if (error()) goto 100

100     call glean(thy(v))
        call glean(thy(hv))
        call glean(thy(h))

        if (error("Exit operators_mod::apply_hamiltonian_h")) continue

        if (.not.error()) call stop_timer("operators: apply_hamiltonian")

      end subroutine

      subroutine apply_overlap_h(v,ov,h)
!doc$ subroutine apply_overlap(v,ov,h)
        type(multivector_obj) :: v
        type(multivector_obj) :: ov
        type(hamiltonian_obj) :: h
!       modifies: ov
!       errors: Inconsistent multibases. Passes errors.

!cod$
        call my(v)
        call my(ov)
        call my(h)

        if (error(x_ghost(h%o%hk%o%mb) /= x_ghost(x_multibasis(v)),"ERROR: inconsistent multibases 1")) goto 100
        if (error(x_ghost(x_multibasis(v)) /= x_ghost(x_multibasis(ov)),"ERROR: inconsistent multibases 2")) goto 100

        ov = v
        if (x_n_projectors(h%o%hk%o%hc%o%ao) /= 0) then
          select case (h%o%hk%o%hc%o%projector_type)
          case (FOURIER_SPACE)
            call apply_fs_overlap_i(v,ov,h%o) ; if (error()) goto 100
          case (REAL_SPACE)
            call apply_rs_overlap_i(v,ov,h%o) ; if (error()) goto 100
          end select
        end if

100     call glean(thy(v))
        call glean(thy(ov))
        call glean(thy(h))

        if (error("Exit operators_mod::apply_overlap_h")) continue

      end subroutine
 
      subroutine add_density_h(h,v,wts,gen_den)
!doc$ subroutine add_density(h,v,wts,gen_den)
        type(hamiltonian_obj) :: h
        type(multivector_obj) :: v
        real(double), dimension(:), intent(in) :: wts
        type(gen_density_obj) :: gen_den
!       modifies: gen_den
!       effects: Updates gen_den with respect to h and v.
!       errors: Passes errors.

!cod$
        complex(double), dimension(:,:), pointer :: pdots
        type(grid_obj) :: grid_den
        type(atomic_density_obj) :: atomic_den

        call my(h)
        call my(v)
        call my(gen_den)

        nullify( pdots )

        call my(x_grid_density(gen_den),grid_den)
        call my(x_atomic_density(gen_den),atomic_den)

        call add_grid_density(v,wts,grid_den)            ; if (error()) goto 100
        if (x_n_projectors(h%o%hk%o%hc%o%ao) /= 0) then
          call form_pdots_i(v,h%o,pdots)                 ; if (error()) goto 100
          call add_atomic_density(atomic_den,pdots,wts)  ; if (error()) goto 100
        end if
        call update(gen_den,grid_den,atomic_den)         ; if (error()) goto 100

100     if (associated( pdots )) deallocate( pdots )
        call glean(thy(grid_den))
        call glean(thy(atomic_den))

        call glean(thy(h))
        call glean(thy(v))
        call glean(thy(gen_den))

        if (error("Exit operators_mod::add_density_h")) continue

      end subroutine

      subroutine orthonormalize_h(mv,h)
!doc$ subroutine orthonormalize(mv,h)
        type(multivector_obj) :: mv
        type(hamiltonian_obj) :: h
!       modifies: mv
!       effects: Orthonormalizes mv.
!       errors: Passes errors.

!cod$
        integer :: nb
        complex(double) :: cmv, ctmv
        complex(double), dimension(:,:), allocatable :: s
        type(multivector_obj) :: tmv

        call my(mv)
        call my(h)
        call my(mv,tmv)

        nb = x_n_bands(mv)
        allocate( s(nb,nb) )
        s = (0.0_double,0.0_double)
        if (overlap_is_identity(h)) then
          call overlap(mv,s)
        else
          call apply_overlap(mv,tmv,h) ; if (error()) goto 100
          call overlap(mv,tmv,s)
          tmv = mv
        end if

        call inverse_cholesky(s) ; if (error()) goto 100
        cmv = cmplx(0,0,double) ; ctmv = cmplx(1,0,double)
        call transform(cmv,mv,ctmv,tmv,s)

100     if (allocated( s )) deallocate( s )

        call glean(thy(tmv))
        call glean(thy(h))
        call glean(thy(mv))

        if (error("Exit operators_mod::orthonormalize_h")) continue

      end subroutine

      subroutine precondition_h(r,v,h)
!doc$ subroutine precondition(r,v,h)
        type(multivector_obj) :: r
        type(multivector_obj) :: v
        type(hamiltonian_obj) :: h
!       modifies: r
!       effects: r <-- precondition(r)
!       errors: Multibases are different for v and h, or v and r.

!cod$
        integer :: nb, ng
        type(multivector_rep), pointer :: worm_v, worm_r

        call my(r)
        call my(v)
        call my(h)

        worm_v => wormhole(v)
        worm_r => wormhole(r)

        if (error(x_ghost(worm_v%mb) /= x_ghost(h%o%hk%o%mb),"ERROR: multibases are different 1")) goto 100
        if (error(x_ghost(worm_v%mb) /= x_ghost(worm_r%mb),"ERROR: multibases are different 2")) goto 100

        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)

        call kernel_precondition_h_i(ng,nb,h%o%hk%o%ke,worm_v%mat,worm_r%mat)
        worm_r%g = x_ghost()

100     nullify( worm_v )
        nullify( worm_r )

        call glean(thy(r))
        call glean(thy(v))
        call glean(thy(h))

        if (error("Exit operators_mod::precondition_h")) continue

      end subroutine

      subroutine forces_h(v,h,wts,eigs,f)
!doc$ subroutine forces(v,h,wts,eigs,f)
        type(multivector_obj) :: v
        type(hamiltonian_obj) :: h
        real(double), dimension(:) :: wts
        real(double), dimension(:) :: eigs
        real(double), dimension(:,:), intent(out) :: f
!       modifies: f
!       requires: size(wts) = size(v%o%mat,2)
!       effects: Returns the forces on atoms due to projectors.
!       errors: Passes errors.

!cod$
        integer :: ia, ib, ic, ip, nb, np
        complex(double), dimension(:,:), allocatable :: o_pdots
        complex(double), dimension(:,:), pointer :: pdots, gpdots

        call my(v)
        call my(h)

        nullify( pdots )
        nullify( gpdots )

        f = 0.0_double
        if (x_n_projectors(h%o%hk%o%hc%o%ao) /= 0) then

          call form_pdots_i(v,h%o,pdots) ; if (error()) goto 100
          np = size(pdots,1)
          nb = size(pdots,2)
          if (overlap_is_identity(h)) then
            call atomic_hamiltonian(h%o%hk%o%hc%o%apot,pdots)
            do ib = 1,nb
              pdots(:,ib) = wts(ib)*pdots(:,ib)
            end do
          else
            allocate( o_pdots(np,nb) )
            o_pdots = pdots
            call atomic_overlap(h%o%hk%o%hc%o%apot,o_pdots)
            call atomic_hamiltonian(h%o%hk%o%hc%o%apot,pdots)
            do ib = 1,nb
              pdots(:,ib) = wts(ib)*pdots(:,ib) - wts(ib)*eigs(ib)*o_pdots(:,ib)
            end do
            deallocate( o_pdots )
          end if

          do ic = 1,3
            call form_grad_pdots_i(v,h%o,ic,gpdots) ; if (error()) goto 100
            do ip = 1,np
              ia = x_projector_atom(h%o%hk%o%hc%o%ao,ip) ; if (error()) goto 100
              f(ic,ia) = f(ic,ia) - 2.0_double*real(sum(conjg(gpdots(ip,:))*pdots(ip,:))) 
            end do
          end do

        end if

100     if (allocated( o_pdots )) deallocate( o_pdots )
        if (associated( pdots )) deallocate( pdots )
        if (associated( gpdots )) deallocate( gpdots )

        call glean(thy(v))
        call glean(thy(h))

        if (error("Exit operators_mod::forces_h")) continue

      end subroutine

      subroutine pressure_h(v,h,wts,p)
!doc$ subroutine pressure(v,h,wts,p)
        type(multivector_obj) :: v
        type(hamiltonian_obj) :: h
        real(double), dimension(:), intent(in) :: wts
        real(double), intent(out) :: p
!       requires: v and h use the same multibasis_obj.
!       effects: Returns hamiltonian contributions to the pressure.
!       errors: Passes errors.

!cod$
        real(double) :: p1, p2

        call my(v)
        call my(h)

        call kinetic_energy_pressure_i(v,h%o%hk,wts,p1) ; if (error()) goto 100
        select case (h%o%hk%o%hc%o%projector_type)
        case (FOURIER_SPACE)
          call projector_pressure_fs_i(v,h%o,wts,p2) ; if (error()) goto 100
        case (REAL_SPACE)
          call warn("WARNING: pressure due to real-space projectors is not yet implemented")
          p2 = 0.0_double
        end select
        p = p1 + p2

100     call glean(thy(v))
        call glean(thy(h))

        if (error("Exit operators_mod::pressure_h")) continue

      end subroutine

      subroutine stress_tensor_h(v,h,wts,s)
!doc$ subroutine stress_tensor(v,h,wts,s)
        type(multivector_obj) :: v
        type(hamiltonian_obj) :: h
        real(double), dimension(:), intent(in) :: wts
        real(double), dimension(:,:), intent(out) :: s
!       requires: v and h use the same multibasis_obj. s be dimension(3,3).
!       effects: Returns hamiltonian contributions to the stress tensor.
!       errors: Passes errors.

!cod$
        real(double), dimension(:,:), allocatable :: s1, s2

        call my(v)
        call my(h)

        allocate( s1(3,3), s2(3,3) )

        call kinetic_energy_stress_tensor_i(v,h%o%hk%o,wts,s1) ; if (error()) goto 100
        select case (h%o%hk%o%hc%o%projector_type)
        case (FOURIER_SPACE)
          call projector_stress_tensor_fs_i(v,h%o,wts,s2) ; if (error()) goto 100
        case (REAL_SPACE)
          call warn("WARNING: stress tensor due to real-space projectors is not yet implemented")
          s2 = 0.0_double
        end select
        s = s1 + s2

100     if (allocated( s1 )) deallocate( s1 )
        if (allocated( s2 )) deallocate( s2 )

        call glean(thy(v))
        call glean(thy(h))

        if (error("Exit operators_mod::stress_tensor_h")) continue

      end subroutine

      function overlap_is_identity_h(h) result(l)
!doc$ function overlap_is_identity(h) result(l)
        type(hamiltonian_obj) :: h
        logical :: l
!       effects: Returns true if the overlap operator is the identity.

!cod$
        call my(h)
        l = overlap_is_identity(h%o%hk%o%hc%o%apot)
        call glean(thy(h))
      end function

! private hc routines

      subroutine own_hc_i(hc)
        type(h_common_obj) :: hc
        type(h_common_obj) :: hct
        integer :: ip, ir
        if (hc%ref < hc%o%ref) then
          allocate( hct%o )
          hct%o%ref = 0
          hct%o%g = hc%o%g
          hct%o%g_gpot = hc%o%g_gpot
          call my(hc%o%crystal,hct%o%crystal)
          call my(hc%o%layout,hct%o%layout)
          call my(hc%o%ao,hct%o%ao)
          call my(hc%o%apot,hct%o%apot)
          allocate( hct%o%lp(size(hc%o%lp,1),size(hc%o%lp,2),size(hc%o%lp,3)) )
          hct%o%lp = hc%o%lp
          hct%o%projector_type = hc%o%projector_type
          hct%o%cutoff = hc%o%cutoff
          hct%o%pbs = hc%o%pbs
          nullify( hct%o%site )
          if (associated( hc%o%site )) then
            allocate( hct%o%site(size(hc%o%site)) )
            do ip = 1,size(hc%o%site)
              hct%o%site(ip)%new = hc%o%site(ip)%new
              hct%o%site(ip)%pos = hc%o%site(ip)%pos
              allocate( hct%o%site(ip)%replica(size(hc%o%site(ip)%replica)) )
              do ir = 1, size(hc%o%site(ip)%replica)
                hct%o%site(ip)%replica(ir)%sr_index = hc%o%site(ip)%replica(ir)%sr_index
                hct%o%site(ip)%replica(ir)%pos = hc%o%site(ip)%replica(ir)%pos
              end do
            end do
          end if
          nullify( hct%o%rsi )
          if (associated( hc%o%rsi )) then
            allocate( hct%o%rsi(size(hc%o%rsi)) )
            hct%o%rsi = hc%o%rsi
          end if
          nullify( hct%o%rsc1 )
          if (associated( hc%o%rsc1 )) then
            allocate( hct%o%rsc1(size(hc%o%rsc1,1),size(hc%o%rsc1,2),size(hc%o%rsc1,3)) )
            hct%o%rsc1 = hc%o%rsc1
          end if
          nullify( hct%o%rsc2 )
          if (associated( hc%o%rsc2 )) then
            allocate( hct%o%rsc2(size(hc%o%rsc2,1),size(hc%o%rsc2,2),size(hc%o%rsc2,3)) )
            hct%o%rsc2 = hc%o%rsc2
          end if
          nullify( hct%o%rsp )
          if (associated( hc%o%rsp )) then
            allocate( hct%o%rsp(size(hc%o%rsp)) )
            hct%o%rsp = hc%o%rsp
          end if
          nullify( hct%o%grsp )
          if (associated( hc%o%grsp )) then
            allocate( hct%o%grsp(size(hc%o%grsp,1),size(hc%o%grsp,2)) )
            hct%o%grsp = hc%o%grsp
          end if
          hc%o%ref = hc%o%ref - hc%ref
          hc%o => hct%o
          hc%o%ref = hc%o%ref + hc%ref
        end if
      end subroutine

      subroutine form_rs_projectors_i(hcr)
        type(h_common_rep) :: hcr
!       effects: Forms the real-space projectors.

        if (x_n_projectors(hcr%ao) /= 0) then
          call get_site_data_i(hcr) ; if (error()) goto 100
          call get_point_data_i(hcr) ; if (error()) goto 100
        end if

100     if (error("Exit operators_mod::form_rs_projectors_i")) continue

      end subroutine

      subroutine get_site_data_i(hcr)
        type(h_common_rep) :: hcr
!       requires: At least one projector exists. hcr%site and underlying pointers be associated or nullified.
!       effects: Identifies replicas (atoms which project into the central parallelpiped) and gathers data about them.
!       errors: If a replica sphere is larger than the parallelpiped.
!       documents: rs_ncp_projectors.

        integer :: i, j, i1, i2, i3, ia, ip, ipr, ir, n, np, nr, start, stop
        integer, dimension(3) :: dim
        integer, dimension(:), allocatable :: count, disp, gnr, lnr
        real(double) :: radius, last_radius
        real(double), dimension(3) :: pos, last_pos, pos_trial, v1, v2, v3, va1, va2, va3
        real(double), dimension(:,:), allocatable :: r
        real(double), dimension(:,:,:), allocatable :: gpr, lpr
        type(atoms_obj) :: ats
        type(lattice_obj) :: lat

        if (associated( hcr%site )) then
          do i = 1,size(hcr%site)
            if (associated( hcr%site(i)%replica )) deallocate( hcr%site(i)%replica )
          end do
          deallocate( hcr%site )
        end if

        call my(x_atoms(hcr%crystal),ats)
        call my(x_lattice(hcr%crystal),lat)

        dim = x_dims(hcr%layout) ! IDENTIFY MESH POINTS BOUNDING THE CENTRAL PARALLELPIPED
        v1 = lat2r(lat,real((/1,0,0/),double))/real(dim(1),double)
        v2 = lat2r(lat,real((/0,1,0/),double))/real(dim(2),double)
        v3 = lat2r(lat,real((/0,0,1/),double))/real(dim(3),double)
        va1 = v1*real(dim(1)-1,double)
        va2 = v2*real(dim(2)-1,double)
        va3 = v3*real(dim(3)-1,double)
        nr = 2*dim(1)*dim(2) + 2*dim(1)*(dim(3)-2) + 2*(dim(2)-2)*(dim(3)-2)
        allocate( r(3,nr) )
        ir = 0
        do i1 = 0,dim(1)-1
          do i2 = 0,dim(2)-1
            r(:,ir+1) = v1*real(i1,double) + v2*real(i2,double)
            r(:,ir+2) = r(:,ir+1) + va3
            ir = ir + 2
         end do
        end do
        do i1 = 0,dim(1)-1
          do i3 = 1,dim(3)-2
            r(:,ir+1) = v1*real(i1,double) + v3*real(i3,double)
            r(:,ir+2) = r(:,ir+1) + va2
            ir = ir + 2
          end do
        end do
        do i2 = 1,dim(2)-2
          do i3 = 1,dim(3)-2
            r(:,ir+1) = v2*real(i2,double) + v3*real(i3,double)
            r(:,ir+2) = r(:,ir+1) + va1
            ir = ir + 2
          end do
        end do

        np = x_n_projectors(hcr%ao)

        call subdivide(mpi_myproc(CONFIG),mpi_nprocs(CONFIG),1,np,start,stop,n) ! IDENTIFY REPLICAS AND GATHER REPLICA DATA
        allocate( lnr(n), gnr(np), lpr(3,8,n), gpr(3,8,np) )
        last_pos = real((/0,0,0/),double)
        last_radius = -1.0_double
        do i = 1,size(lnr)
          ip = start + i - 1
          ia = x_projector_atom(hcr%ao,ip) ; if (error()) goto 100
          pos = lat2r(lat,modulo(x_position(ats,ia),real((/1,1,1/),double)))
          radius = x_projector_radius(hcr%ao,ip) ; if (error()) goto 100
          if (all(pos == last_pos) .and. (radius == last_radius)) then
            lnr(i) = lnr(i-1)
            lpr(:,1:lnr(i),i) = lpr(:,1:lnr(i-1),i-1)
          else
            if (error(radius >= maximum_sphere_radius(lat),"ERROR: radius is too large")) goto 100
            lnr(i) = 0
            do i1 = -1,+1
              do i2 = -1,+1
                do i3 = -1,+1
                  if ( all((/i1,i2,i3/) == 0)) then
                    lnr(i) = lnr(i) + 1
                    lpr(:,lnr(i),i) = pos
                    cycle
                  end if
                  pos_trial = lat2r(lat,real((/i1,i2,i3/),double)) + pos
                  do ir = 1,nr
                    if (norm(r(:,ir) - pos_trial) <= radius) then
                      lnr(i) = lnr(i) + 1
                      lpr(:,lnr(i),i) = pos_trial
                      exit
                    end if
                  end do
                end do
              end do
            end do
            last_pos = pos
            last_radius = radius
          end if
        end do
        allocate( count(0:mpi_nprocs(CONFIG)-1), disp(0:mpi_nprocs(CONFIG)-1) )
        do i = 0,mpi_nprocs(CONFIG)-1
          call subdivide(i,mpi_nprocs(CONFIG),1,np,start,stop,n)
          count(i) = n
          disp(i) = 0
          do j = 0,i-1
            disp(i) = disp(i) + count(j)
          end do
        end do
        call allgatherv(CONFIG,lnr,gnr,count,disp) ; if (error()) goto 100
        count = 3*8*count ; disp = 3*8*disp
        call allgatherv(CONFIG,lpr,gpr,count,disp) ; if (error()) goto 100

        allocate( hcr%site(np) ) ! STORE THE REPLICA DATA
        last_pos = real((/0,0,0/),double)
        last_radius = -1.0_double
        do ip = 1,np
          ia = x_projector_atom(hcr%ao,ip) ; if (error()) goto 100
          pos = lat2r(lat,modulo(x_position(ats,ia),real((/1,1,1/),double)))
          radius = x_projector_radius(hcr%ao,ip) ; if (error()) goto 100
          if (all(pos == last_pos) .and. (radius == last_radius)) then
            hcr%site(ip)%new = .false.
            hcr%site(ip)%pos = hcr%site(ip-1)%pos
            allocate( hcr%site(ip)%replica(size(hcr%site(ip-1)%replica)) )
            do ir = 1,size(hcr%site(ip)%replica)
              hcr%site(ip)%replica(ir)%pos = hcr%site(ip-1)%replica(ir)%pos
            end do
          else
            hcr%site(ip)%new = .true.
            hcr%site(ip)%pos = pos
            allocate( hcr%site(ip)%replica(gnr(ip)) )
            do ir = 1,size(hcr%site(ip)%replica)
              hcr%site(ip)%replica(ir)%pos = gpr(:,ir,ip)
            end do
            last_pos = pos
            last_radius = radius
          end if
        end do
        ipr = 0
        do ip = 1,size(hcr%site)
          do ir = 1,size(hcr%site(ip)%replica)
            ipr = ipr + 1
            hcr%site(ip)%replica(ir)%sr_index = ipr
          end do
        end do

100     if (allocated( r )) deallocate( r )
        if (allocated( lnr )) deallocate( lnr )
        if (allocated( lpr )) deallocate( lpr )
        if (allocated( gnr )) deallocate( gnr )
        if (allocated( gpr )) deallocate( gpr )
        if (allocated( count )) deallocate( count )
        if (allocated( disp )) deallocate( disp )

        call glean(thy(ats))
        call glean(thy(lat))

        if (error("Exit operators_mod::get_site_data_i")) continue

      end subroutine

      subroutine get_point_data_i(hcr)
        type(h_common_rep) :: hcr
!       requires: At least one projector exists.
!                 hcr%rsi, hcr%rsc1, hcr%rsc2, hcr%rsp, hcr%grsp be associated or nullified.
!                 Projector m values for a given l value be in increasing order.
!       effects: Computes projector and projector gradient values at real-space mesh points.
!       documents: rs_ncp_projectors.
!       notes: CONFIG processes work as a team in this routine.

        logical :: hit
        integer :: first, i, ip, ir, j, last, last_l, last_m, n, n_dim
        integer, dimension(3) :: c, dim
        integer, dimension(:), allocatable :: count, disp, gnp, l, m, lnp, rsi_loc
        integer, dimension(:,:), allocatable :: pr_index
        real(double) :: gi_cut, go_cut, prv, rn
        real(double), dimension(2) :: prg
        real(double), dimension(3) :: r, rm
        real(double), dimension(:), allocatable :: radius, rsp_loc
        real(double), dimension(:,:), allocatable :: grsp_loc
        type(point_data), dimension(:), pointer :: point

        nullify( point )

        if (associated( hcr%rsi ))  deallocate( hcr%rsi )
        if (associated( hcr%rsc1 )) deallocate( hcr%rsc1 )
        if (associated( hcr%rsc2 )) deallocate( hcr%rsc2 )
        if (associated( hcr%rsp ))  deallocate( hcr%rsp )
        if (associated( hcr%grsp )) deallocate( hcr%grsp )

        dim = x_dims(hcr%layout)
        n_dim = product(dim)

        allocate( count(0:mpi_nprocs(CONFIG)-1), disp(0:mpi_nprocs(CONFIG)-1) )

        gi_cut = sqrt(hcr%cutoff)
        go_cut = 2.0_double*sqrt(x_cutoff(hcr%layout)) - gi_cut

        n = x_n_projectors(hcr%ao)
        allocate( l(n), m(n), radius(n), pr_index(2,n) )
        do i = 1,n
          l(i) = x_projector_l(hcr%ao,i) ; if (error()) goto 100
          m(i) = x_projector_m(hcr%ao,i) ; if (error()) goto 100
          radius(i) = x_projector_radius(hcr%ao,i) ; if (error()) goto 100
        end do

        call subdivide(mpi_myproc(CONFIG),mpi_nprocs(CONFIG),1,n_dim,first,last,n)
        allocate( lnp(n), point(n) )

        do i = 1,size(lnp)

          c = map3d(dim,i+first-1)                 ! IDENTIFY REPLICA MESH POINTS INSIDE THE CENTRAL PARALLELPIPED
          rm = a2rc(c,hcr%layout,S_TYPE)
          hit = .false.
          n = 0
          do ip = 1,size(hcr%site)
            if (hcr%site(ip)%new) then
              do ir = 1,size(hcr%site(ip)%replica)
                r = rm - hcr%site(ip)%replica(ir)%pos
                if (norm(r) <= radius(ip)) then
                  n = n + 1
                  pr_index(:,n) = (/ip,ir/)
                  hit = .true.
                  exit
                else
                  hit = .false.
                end if
              end do
            else
              if (hit) then
                n = n + 1
                pr_index(1,n) = ip
                pr_index(2,n) = pr_index(2,n-1)
              end if
            end if
          end do
          lnp(i) = n

          allocate( point(i)%i(lnp(i)) )   ! COMPUTE PROJECTOR AND PROJECTOR GRADIENT VALUES AT THE REPLICA MESH POINTS
          allocate( point(i)%p(lnp(i)) )
          allocate( point(i)%pg(3,lnp(i)) )
          last_l = -1
          last_m = -100
          do j = 1,lnp(i)
            ip = pr_index(1,j)
            ir = pr_index(2,j)
            r = rm - hcr%site(ip)%replica(ir)%pos
            rn = norm(r)
            if ((l(ip) /= last_l) .or. (m(ip) /= (last_m+1))) then
              prv = projector_r_value(hcr%ao,ip,rn,gi_cut,go_cut) ; if (error()) goto 100
              prg = projector_r_gradients(hcr%ao,ip,rn) ; if (error()) goto 100
            end if
            point(i)%i(j) = hcr%site(ip)%replica(ir)%sr_index
            point(i)%p(j) = prv*real_spherical_harmonic(l(ip),m(ip),r)
            point(i)%pg(:,j) = prg(2)*r*real_spherical_harmonic(l(ip),m(ip),r) - prg(1)*grad_real_spherical_harmonic(l(ip),m(ip),r)
            last_l = l(ip)
            last_m = m(ip)
          end do

        end do

        allocate( gnp(n_dim) )                     ! STORE THE PROJECTOR VALUES
        do i = 0,mpi_nprocs(CONFIG)-1
          call subdivide(i,mpi_nprocs(CONFIG),1,n_dim,first,last,n)
          count(i) = n
          disp(i) = 0
          do j = 0,i-1
            disp(i) = disp(i) + count(j)
          end do
        end do
        call allgatherv(CONFIG,lnp,gnp,count,disp) ; if (error()) goto 100

        allocate( hcr%rsc1(dim(1),dim(2),dim(3)), hcr%rsc2(dim(1),dim(2),dim(3)) )
        n = 0
        do i = 1,size(gnp)
          c = map3d(dim,i)
          hcr%rsc1(c(1),c(2),c(3)) = n + 1
          hcr%rsc2(c(1),c(2),c(3)) = n + gnp(i)
          n = n + gnp(i)
        end do

        n = sum(lnp)
        allocate( rsi_loc(n), rsp_loc(n) )
        n = 0
        do i = 1,size(lnp)
          do j = 1,lnp(i)
            n = n + 1
            rsi_loc(n) = point(i)%i(j)
            rsp_loc(n) = point(i)%p(j)
          end do
        end do

        n = sum(gnp)
        allocate( hcr%rsi(n), hcr%rsp(n) )
        do i = 0,mpi_nprocs(CONFIG)-1
          call subdivide(i,mpi_nprocs(CONFIG),1,n_dim,first,last,n)
          count(i) = sum(gnp(first:last))
          disp(i) = 0
          do j = 0,i-1
            disp(i) = disp(i) + count(j)
          end do
        end do
        call allgatherv(CONFIG,rsi_loc,hcr%rsi,count,disp) ; if (error()) goto 100
        call allgatherv(CONFIG,rsp_loc,hcr%rsp,count,disp) ; if (error()) goto 100

        n = sum(lnp)                               ! STORE THE PROJECTOR GRADIENT VALUES
        allocate( grsp_loc(3,n) )
        n = 0
        do i = 1,size(lnp)
          do j = 1,lnp(i)
            n = n + 1
            grsp_loc(:,n) = point(i)%pg(:,j)
          end do
        end do

        do i = 0,mpi_nprocs(CONFIG)-1
          call subdivide(i,mpi_nprocs(CONFIG),1,n_dim,first,last,n)
          count(i) = 0
          do j = first,last
            c = map3d(dim,j)
            count(i) = count(i) + hcr%rsc2(c(1),c(2),c(3)) - hcr%rsc1(c(1),c(2),c(3)) + 1
          end do
          disp(i) = 0
          do j = 0,i-1
            disp(i) = disp(i) + count(j)
          end do
        end do
        allocate( hcr%grsp(3,sum(count)) )
        count = 3*count
        disp = 3*disp
        call allgatherv(CONFIG,grsp_loc,hcr%grsp,count,disp) ; if (error()) goto 100

100     if (allocated( pr_index )) deallocate( pr_index )
        if (allocated( l )) deallocate( l )
        if (allocated( m )) deallocate( m )
        if (allocated( radius )) deallocate( radius )
        if (allocated( lnp )) deallocate( lnp )
        if (allocated( gnp )) deallocate( gnp )
        if (allocated( count )) deallocate( count )
        if (allocated( disp )) deallocate( disp )
        if (associated( point )) then
          do i = 1,size(point)
            if (associated( point(i)%i )) deallocate( point(i)%i )
            if (associated( point(i)%p )) deallocate( point(i)%p )
            if (associated( point(i)%pg )) deallocate( point(i)%pg )
          end do
          deallocate( point )
        end if
        if (allocated( rsi_loc )) deallocate( rsi_loc )
        if (allocated( rsp_loc )) deallocate( rsp_loc )
        if (allocated( grsp_loc )) deallocate( grsp_loc )

        if (error("Exit operators_mod::get_point_data_i")) continue

      end subroutine

! private hk routines

      subroutine own_hk_i(hk)
        type(h_kpoint_obj) :: hk
        type(h_kpoint_obj) :: hkt
        integer :: ip, ir, it
        if (hk%ref < hk%o%ref) then
          allocate( hkt%o )
          hkt%o%ref = 0
          hkt%o%g = hk%o%g
          call my(hk%o%hc,hkt%o%hc)
          call my(hk%o%mb,hkt%o%mb)
          if (associated( hk%o%dk )) then
             allocate(hkt%o%dk)
             call my(hk%o%dk, hkt%o%dk)
          else
             nullify(hkt%o%dk)
          end if
          nullify( hkt%o%type )
          if (associated( hk%o%type )) then
            allocate( hkt%o%type(size(hk%o%type)) )
            do it = 1,size(hk%o%type)
              nullify( hkt%o%type(it)%form_factor )
              if (associated( hk%o%type(it)%form_factor )) then
                allocate( hkt%o%type(it)%form_factor(size(hk%o%type(it)%form_factor,1),size(hk%o%type(it)%form_factor,2)) )
                hkt%o%type(it)%form_factor = hk%o%type(it)%form_factor
              end if
            end do
          end if
          nullify( hkt%o%site )
          if (associated( hk%o%site )) then
            allocate( hkt%o%site(size(hk%o%site)) )
            do ip = 1,size(hk%o%site)
              allocate( hkt%o%site(ip)%replica(size(hk%o%site(ip)%replica)) )
              do ir = 1,size(hk%o%site(ip)%replica)
                hkt%o%site(ip)%replica(ir)%phase = hk%o%site(ip)%replica(ir)%phase
              end do
            end do
          end if
          allocate( hkt%o%ke(size(hk%o%ke)) )
          hkt%o%ke = hk%o%ke
          hkt%o%need_gk = hk%o%need_gk
          nullify( hkt%o%gkx )
          if (associated(hk%o%gkx)) then
             allocate( hkt%o%gkx(size(hk%o%gkx)) )
             hkt%o%gkx = hk%o%gkx
          end if
          nullify( hkt%o%gky )
          if (associated(hk%o%gky)) then
             allocate( hkt%o%gky(size(hk%o%gky)) )
             hkt%o%gky = hk%o%gky
          end if
          nullify( hkt%o%gkz )
          if (associated(hk%o%gkz)) then
             allocate( hkt%o%gkz(size(hk%o%gkz)) )
             hkt%o%gkz = hk%o%gkz
          end if
          hk%o%ref = hk%o%ref - hk%ref
          hk%o => hkt%o
          hk%o%ref = hk%o%ref + hk%ref
        end if
      end subroutine

      subroutine apply_kinetic_energy_i(v,hv,hkr)
        type(multivector_obj) :: v
        type(multivector_obj) :: hv
        type(h_kpoint_rep) :: hkr
!       requires: hv contain valid data.
!       modifies: hv

        type(multivector_rep), pointer :: worm_v, worm_hv

        call my(v)
        call my(hv)

        worm_v => wormhole(v)
        worm_hv => wormhole(hv)

        call point_mamxv(worm_hv%mat,worm_v%mat,hkr%ke)
        worm_hv%g = x_ghost()

        nullify( worm_v )
        nullify( worm_hv )

        call glean(thy(v))
        call glean(thy(hv))

      end subroutine

      subroutine kinetic_energy_pressure_i(v,hk,wts,p)
        type(multivector_obj) :: v
        type(h_kpoint_obj) :: hk
        real(double), dimension(:), intent(in) :: wts
        real(double), intent(out) :: p
!       effects: Returns the pressure contribution due to the kinetic energy.
!       errors: Passes errors.

        real(double) :: cell_volume

        call my(v)
        call my(hk)

        cell_volume = x_cell_volume(x_lattice(hk%o%hc%o%crystal))
        p = (2.0_double/(3.0_double*cell_volume))*kinetic_energy(v,hk,wts) ; if (error()) goto 100

100     call glean(thy(v))
        call glean(thy(hk))

        if (error("Exit operators_mod::kinetic_energy_pressure_i")) continue

      end subroutine

      subroutine kinetic_energy_stress_tensor_i(v,hkr,wts,s)
        type(multivector_obj) :: v
        type(h_kpoint_rep) :: hkr
        real(double), dimension(:), intent(in) :: wts
        real(double), dimension(:,:), intent(out) :: s
!       effects: Returns kinetic energy contributions to the stress tensor.

        integer :: ib
        real(double) :: cell_volume
        real(double), dimension(:), allocatable :: gkx, gky, gkz, tmp
        real(double), dimension(:,:), allocatable :: s_local
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb

        call my(v)

        nullify( worm_v )
        nullify( worm_mb )

        cell_volume = x_cell_volume(x_lattice(hkr%hc%o%crystal))

        worm_v => wormhole(v)
        worm_mb => wormhole(hkr%mb)
        allocate( tmp(size(worm_v%mat,1)) )
        allocate( gkx(size(worm_mb%gpt,1)), gky(size(worm_mb%gpt,1)), gkz(size(worm_mb%gpt,1)) )
        gkx = worm_mb%gpt(:,1) + worm_mb%kpt(1)
        gky = worm_mb%gpt(:,2) + worm_mb%kpt(2)
        gkz = worm_mb%gpt(:,3) + worm_mb%kpt(3)

        allocate( s_local(3,3) )
        s_local = 0.0_double
        do ib = 1,size(worm_v%mat,2)
          tmp = (2.0_double*wts(ib)/cell_volume)*(worm_v%mat(:,ib)*conjg(worm_v%mat(:,ib)))
          s_local(1,1) = s_local(1,1) - sum(tmp*gkx*gkx)
          s_local(1,2) = s_local(1,2) - sum(tmp*gkx*gky)
          s_local(1,3) = s_local(1,3) - sum(tmp*gkx*gkz)
          s_local(2,2) = s_local(2,2) - sum(tmp*gky*gky)
          s_local(2,3) = s_local(2,3) - sum(tmp*gky*gkz)
          s_local(3,3) = s_local(3,3) - sum(tmp*gkz*gkz)
        end do
        s_local(2,1) = s_local(1,2)
        s_local(3,1) = s_local(1,3)
        s_local(3,2) = s_local(2,3)
        call allreduce(KGROUP,MPI_SUM,s_local,s)

100     if (allocated( gkx )) deallocate( gkx )
        if (allocated( gky )) deallocate( gky )
        if (allocated( gkz )) deallocate( gkz )
        if (allocated( tmp )) deallocate( tmp )
        if (allocated( s_local )) deallocate( s_local )
        nullify( worm_v )
        nullify( worm_mb )

        call glean(thy(v))

        if (error("Exit operators_mod::kinetic_energy_stress_tensor_i")) continue

      end subroutine

      subroutine form_kinetic_energy_i(hkr)
        type(h_kpoint_rep) :: hkr
!       requires: hkr%ke, hkr%gkx, hkr%gky, hkr%gkz be nullified or associated.
!       effects: Forms hkr%ke, hkr%gkx, hkr%gky, and hkr%gkz.

        if (associated( hkr%ke )) deallocate( hkr%ke )
        if (associated( hkr%gkx )) deallocate( hkr%gkx )
        if (associated( hkr%gky )) deallocate( hkr%gky )
        if (associated( hkr%gkz )) deallocate( hkr%gkz )

        call form_gk_i(hkr) ; if (error()) goto 100
        allocate( hkr%ke(size(hkr%gkx)) )
        hkr%ke = (hkr%gkx**2 + hkr%gky**2 + hkr%gkz**2)
        if (.not.hkr%need_gk) then
          deallocate( hkr%gkx )
          deallocate( hkr%gky )
          deallocate( hkr%gkz )
        end if

100     if (error("Exit operators_mod::form_kinetic_energy_i")) continue

      end subroutine

      subroutine form_gk_i(hkr)
        type(h_kpoint_rep) :: hkr
!       requires: hkr%gkx, hkr%gky%, hkr%gkz be nullified or associated.
!       effects: Forms hkr%gkx, hkr%gky, and hkr%gkz.

        integer :: ng
        type(multibasis_rep), pointer :: worm_mb

        nullify( worm_mb )

        if (associated( hkr%gkx )) deallocate( hkr%gkx )
        if (associated( hkr%gky )) deallocate( hkr%gky )
        if (associated( hkr%gkz )) deallocate( hkr%gkz )

        worm_mb => wormhole(hkr%mb)
      
        ng = size(worm_mb%gpt,1)
        allocate( hkr%gkx(ng) )
        allocate( hkr%gky(ng) )
        allocate( hkr%gkz(ng) )

        hkr%gkx = worm_mb%gpt(:,1) + worm_mb%kpt(1)
        hkr%gky = worm_mb%gpt(:,2) + worm_mb%kpt(2)
        hkr%gkz = worm_mb%gpt(:,3) + worm_mb%kpt(3)

100     nullify( worm_mb )

        if (error("Exit operators_mod::form_gk_i")) continue

      end subroutine

      subroutine form_fs_type_projectors_i(hkr)
        type(h_kpoint_rep) :: hkr
!       requires: hkr%type be nullified or associated. Projector m values for a given l value be in increasing order.
!       effects: Forms hkr%type and underlying data structures.
!       documents: ncp_projectors.

        integer :: ig, it, itp, l, last_l, last_m, m, n, ntp
        real(double) :: rsh, sr_icv
        real(double), dimension(:), allocatable :: gkm, pfv
        real(double), dimension(:,:), allocatable :: gk
        type(multibasis_rep), pointer :: worm_mb

        nullify( worm_mb )

        if (associated( hkr%type )) then
          do it = 1,size(hkr%type)
            if (associated( hkr%type(it)%form_factor )) deallocate( hkr%type(it)%form_factor )
          end do
          deallocate( hkr%type )
        end if

        if (x_n_projectors(hkr%hc%o%ao) /= 0) then

          sr_icv = sqrt(1.0_double/x_cell_volume(x_lattice(hkr%hc%o%crystal)))

          worm_mb => wormhole(hkr%mb)
          n = size(worm_mb%gpt,1)
          allocate( gk(3,n), gkm(n), pfv(n) )
          do ig = 1,size(gk,2)
            gk(1,ig) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2,ig) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3,ig) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = norm(gk(:,ig))
          end do
          nullify( worm_mb )

          n = x_n_types(hkr%hc%o%ao)
          allocate( hkr%type(n) )

          do it = 1,size(hkr%type)
            nullify( hkr%type(it)%form_factor )
            ntp = x_n_type_projectors(hkr%hc%o%ao,it) ; if (error()) goto 100
            if (ntp == 0) cycle
            allocate( hkr%type(it)%form_factor(size(pfv),ntp) )
            last_l = -1
            last_m = 0
            do itp = 1,ntp
              l = x_projector_l(hkr%hc%o%ao,it,itp) ; if (error()) goto 100
              m = x_projector_m(hkr%hc%o%ao,it,itp) ; if (error()) goto 100
              !if ((l /= last_l) .or. (m /= last_m + 1)) then
                pfv = projector_f_values(hkr%hc%o%ao,gkm,it,itp) ; if (error()) goto 100
              !end if
              do ig = 1,size(gk,2)
                rsh = real_spherical_harmonic(l,m,gk(:,ig))
                hkr%type(it)%form_factor(ig,itp) = pfv(ig)*rsh*sr_icv
              end do
              last_l = l
              last_m = m
            end do
          end do
          deallocate( gk, gkm, pfv )

        end if

100     if (allocated( gk )) deallocate( gk )
        if (allocated( gkm )) deallocate( gkm )
        if (allocated( pfv )) deallocate( pfv )

        nullify( worm_mb )

        if (error("Exit operators_mod::form_fs_type_projectors_i")) continue

      end subroutine

      subroutine form_rs_replica_phases_i(hkr)
        type(h_kpoint_rep) :: hkr
!       requires: hkr%site and its underlying pointers not be associated.
!       effects: Forms hkr%site and its underlying data structures.

        integer :: ip, ir
        real(double), dimension(3) :: k, r
        complex(double), parameter :: i = (0.0_double,1.0_double)
        type(multibasis_rep), pointer :: worm_mb

        nullify( worm_mb )

        if (associated( hkr%site )) then
          do ip = 1,size(hkr%site)
            if (associated( hkr%site(ip)%replica )) deallocate( hkr%site(ip)%replica )
          end do
          deallocate( hkr%site )
        end if

        if (x_n_projectors(hkr%hc%o%ao) /= 0) then

          worm_mb => wormhole(hkr%mb)
          k = worm_mb%kpt
          nullify( worm_mb )

          allocate( hkr%site(size(hkr%hc%o%site)) )
          do ip = 1,size(hkr%site)
            allocate( hkr%site(ip)%replica(size(hkr%hc%o%site(ip)%replica)) )
            do ir = 1,size(hkr%site(ip)%replica)
              r = hkr%hc%o%site(ip)%replica(ir)%pos - hkr%hc%o%site(ip)%pos
              hkr%site(ip)%replica(ir)%phase = exp(-i*(k(1)*r(1) + k(2)*r(2) + k(3)*r(3)))
            end do
          end do

        end if

100     if (error("Exit operators_mod::form_rs_replica_phases_i")) continue

      end subroutine

! private h routines

      subroutine zero_mv_i(v)
        type(multivector_obj) :: v

        type(multivector_rep), pointer :: worm_v

        call my(v)

        worm_v => wormhole(v)
        worm_v%mat = (0.0_double,0.0_double)
        worm_v%g = x_ghost()
        nullify( worm_v )

        call glean(thy(v))

      end subroutine

      subroutine apply_potentials_i(v,hv,hr)
        type(multivector_obj) :: v
        type(multivector_obj) :: hv
        type(hamiltonian_rep) :: hr
!       requires: hv contain valid data.

        call my(v)
        call my(hv)

        if (x_n_projectors(hr%hk%o%hc%o%ao) == 0) then
          call apply_local_i(v,hv,hr) ; if (error()) goto 100
        else
          select case (hr%hk%o%hc%o%projector_type)
          case (FOURIER_SPACE)
            call apply_local_i(v,hv,hr) ; if (error()) goto 100
            call apply_fs_projectors_i(v,hv,hr) ; if (error()) goto 100
          case (REAL_SPACE)
            call apply_local_and_rs_projectors_i(v,hv,hr) ; if (error()) goto 100
          end select
        end if

100     call glean(thy(v))
        call glean(thy(hv))

        if (error("Exit operators_mod::apply_potentials_i")) continue

      end subroutine

      subroutine apply_local_i(v,hv,hr)
        type(multivector_obj) :: v
        type(multivector_obj) :: hv
        type(hamiltonian_rep) :: hr

        integer :: grp, nb, ng, ngt, npr
        integer, dimension(3) :: nd
        complex(double), dimension(:), allocatable :: vec
        complex(double), dimension(:,:), allocatable :: mat
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multivector_rep), pointer :: worm_v, worm_hv
        type(multibasis_rep), pointer :: worm_mb
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        call my(v)
        call my(hv)

        worm_v => wormhole(v)
        worm_hv => wormhole(hv)
        worm_mb => wormhole(worm_v%mb)

        hcr => hr%hk%o%hc%o
        hkr => hr%hk%o

        wf1 => datalink(worm_mb%isplan1)

        nd = x_dims(hcr%layout)
        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)
        ngt = size(worm_mb%gridmap,2)
        npr = mpi_nprocs(KGROUP)
        allocate( vec(ngt), mat(ng,npr) )

        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_v%mat)
          call band_remap(hkr%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call kernel_mm_cr_i(nd,wf1,hcr%lp)
            call fft_serial(R_TO_Q,worm_mb%isplan1)
            call gather_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
          end if
          call band_remap(hkr%mb,grp,vec,mat) ; if (error()) goto 100
          call augment_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_hv%mat)
        end do
        worm_hv%g = x_ghost()

100     if (allocated( vec )) deallocate( vec )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )
        nullify( worm_v, worm_hv )
        nullify( worm_mb )
        nullify( hcr )
        nullify( hkr )

        call glean(thy(v))
        call glean(thy(hv))

        if (error("Exit operators_mod::apply_local_i")) continue

      end subroutine

      subroutine kernel_mm_cr_i(nd,c1,r1)
        integer, dimension(3) :: nd
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1
        real(double), dimension(nd(1),nd(2),nd(3)) :: r1

        integer :: i1, i2, i3

        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              c1(i1,i2,i3) = c1(i1,i2,i3)*r1(i1,i2,i3)
            end do
          end do
        end do

      end subroutine

      subroutine apply_fs_projectors_i(v,hv,hr)
        type(multivector_obj) :: v
        type(multivector_obj) :: hv
        type(hamiltonian_rep) :: hr
!       requires: At least one projector exists.

        character(1), parameter :: transa = 'n', transb = 'n'
        integer :: nb, np, nw
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (1.0_double,0.0_double)
        complex(double), dimension(:,:), pointer :: pdots
        type(multivector_rep), pointer :: worm_hv

        call my(v)
        call my(hv)

        nullify( pdots )
        nullify( worm_hv )

        call form_pdots_i(v,hr,pdots) ; if (error()) goto 100

        call atomic_hamiltonian(hr%hk%o%hc%o%apot,pdots)

        np = size(pdots,1)
        nb = size(pdots,2)
        nw = size(hr%vnl,1)

        worm_hv => wormhole(hv)

        call start_timer("operators: zgemm")
        call zgemm(transa,transb,nw,nb,np,alpha,hr%vnl,nw,pdots,np,beta,worm_hv%mat,nw)
        call stop_timer("operators: zgemm")
        worm_hv%g = x_ghost()

100     if (associated( pdots )) deallocate( pdots )
        nullify( worm_hv )

        call glean(thy(v))
        call glean(thy(hv))

        if (error("Exit operators_mod::apply_fs_projectors_i")) continue

      end subroutine

      subroutine apply_local_and_rs_projectors_i(v,hv,hr)
        type(multivector_obj) :: v
        type(multivector_obj) :: hv
        type(hamiltonian_rep) :: hr
!       requires: At least one projector exists.

        integer :: grp, nb, ng, ngt, np, npr, nr, ns, nsr
        integer, dimension(3) :: nd
        real(double) :: mesh_volume
        complex(double), dimension(:), allocatable :: pdots, srdots, vec
        complex(double), dimension(:,:), allocatable :: mat
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multivector_rep), pointer :: worm_v, worm_hv
        type(multibasis_rep), pointer :: worm_mb
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        call my(v)
        call my(hv)

        worm_v => wormhole(v)
        worm_hv => wormhole(hv)
        worm_mb => wormhole(worm_v%mb)

        hcr => hr%hk%o%hc%o
        hkr => hr%hk%o

        nd = x_dims(hcr%layout)
        mesh_volume = x_cell_volume(x_lattice(hcr%crystal))/real(product(nd),double)

        wf1 => datalink(worm_mb%isplan1)

        ns = size(hcr%site)
        nr = size(hcr%site(ns)%replica)
        nsr = hcr%site(ns)%replica(nr)%sr_index
        allocate( srdots(nsr), pdots(ns) )

        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)
        ngt = size(worm_mb%gridmap,2)
        np = size(hcr%rsi)
        npr = mpi_nprocs(KGROUP)
        allocate( vec(ngt), mat(ng,npr) )

        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_v%mat)
          call band_remap(hkr%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call form_srdots_i(nd,wf1,hr%phase,hcr%rsc1,hcr%rsc2,np,hcr%rsi,hcr%rsp,nsr,srdots)
            call combine_srdots_i(srdots,pdots,hcr,hkr)
            pdots = mesh_volume*pdots
            call atomic_hamiltonian(hcr%apot,pdots)
            call disperse_pdots_i(srdots,pdots,hcr,hkr)
            call distribute_srdots_and_lp_i(nd,wf1,hr%phase,hcr%lp,hcr%rsc1,hcr%rsc2,np,hcr%rsi,hcr%rsp,nsr,srdots)
            call fft_serial(R_TO_Q,worm_mb%isplan1)
            call gather_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
          end if
          call band_remap(hkr%mb,grp,vec,mat) ; if (error()) goto 100
          call augment_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_hv%mat)
        end do
        worm_hv%g = x_ghost()

100     if (allocated( pdots )) deallocate( pdots )
        if (allocated( srdots )) deallocate( srdots )
        if (allocated( vec )) deallocate( vec )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )
        nullify( worm_v, worm_hv )
        nullify( worm_mb )
        nullify( hcr )
        nullify( hkr )

        call glean(thy(v))
        call glean(thy(hv))

        if (error("Exit operators_mod::apply_local_and_rs_projectors_i")) continue

      end subroutine

      subroutine apply_fs_overlap_i(v,ov,hr)
        type(multivector_obj) :: v
        type(multivector_obj) :: ov
        type(hamiltonian_rep) :: hr
!       requires: At least one projector exist.

        character(1), parameter :: transa = 'n', transb = 'n'
        integer :: nb, np, nw
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (1.0_double,0.0_double)
        complex(double), dimension(:,:), pointer :: pdots
        type(multivector_rep), pointer :: worm_ov

        call my(v)
        call my(ov)

        nullify( pdots )
        nullify( worm_ov )

        call form_pdots_i(v,hr,pdots) ; if (error()) goto 100

        call atomic_overlap(hr%hk%o%hc%o%apot,pdots)

        np = size(pdots,1)
        nb = size(pdots,2)
        nw = size(hr%vnl,1)

        worm_ov => wormhole(ov)

        call start_timer("operators: zgemm")
        call zgemm(transa,transb,nw,nb,np,alpha,hr%vnl,nw,pdots,np,beta,worm_ov%mat,nw)
        call stop_timer("operators: zgemm")
        worm_ov%g = x_ghost()

100     if (associated( pdots )) deallocate( pdots )
        nullify( worm_ov )

        call glean(thy(v))
        call glean(thy(ov))

        if (error("Exit operators_mod::apply_fs_overlap_i")) continue

      end subroutine

      subroutine apply_rs_overlap_i(v,ov,hr)
        type(multivector_obj) :: v
        type(multivector_obj) :: ov
        type(hamiltonian_rep) :: hr
!       requires: At least one projector exist.
!       documents: rs_ncp_projectors.

        integer :: grp, nb, ng, ngt, np, npr, nr, ns, nsr
        integer, dimension(3) :: nd
        real(double) :: mesh_volume
        complex(double), dimension(:), allocatable :: pdots, srdots, vec
        complex(double), dimension(:,:), allocatable :: mat
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multivector_rep), pointer :: worm_v, worm_ov
        type(multibasis_rep), pointer :: worm_mb
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        call my(v)
        call my(ov)

        worm_v => wormhole(v)
        worm_ov => wormhole(ov)
        worm_mb => wormhole(worm_v%mb)

        hcr => hr%hk%o%hc%o
        hkr => hr%hk%o

        nd = x_dims(hcr%layout)
        mesh_volume = x_cell_volume(x_lattice(hcr%crystal))/real(product(nd),double)

        wf1 => datalink(worm_mb%isplan1)

        ns = size(hcr%site)
        nr = size(hcr%site(ns)%replica)
        nsr = hcr%site(ns)%replica(nr)%sr_index
        allocate( srdots(nsr), pdots(ns) )

        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)
        ngt = size(worm_mb%gridmap,2)
        np = size(hcr%rsi)
        npr = mpi_nprocs(KGROUP)
        allocate( vec(ngt), mat(ng,npr) )

        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_v%mat)
          call band_remap(hkr%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call form_srdots_i(nd,wf1,hr%phase,hcr%rsc1,hcr%rsc2,np,hcr%rsi,hcr%rsp,nsr,srdots)
            call combine_srdots_i(srdots,pdots,hcr,hkr)
            pdots = mesh_volume*pdots
            call atomic_overlap(hcr%apot,pdots)
            call disperse_pdots_i(srdots,pdots,hcr,hkr)
            call distribute_srdots_i(nd,wf1,hr%phase,hcr%rsc1,hcr%rsc2,np,hcr%rsi,hcr%rsp,nsr,srdots)
            call fft_serial(R_TO_Q,worm_mb%isplan1)
            call gather_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
          end if
          call band_remap(hkr%mb,grp,vec,mat) ; if (error()) goto 100
          call augment_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_ov%mat)
        end do
        worm_ov%g = x_ghost()

100     if (allocated( pdots )) deallocate( pdots )
        if (allocated( srdots )) deallocate( srdots )
        if (allocated( vec )) deallocate( vec )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )
        nullify( worm_v, worm_ov )
        nullify( worm_mb )
        nullify( hcr )
        nullify( hkr )

        call glean(thy(v))
        call glean(thy(ov))

        if (error("Exit operators_mod::apply_rs_overlap_i")) continue

      end subroutine

      subroutine kernel_precondition_h_i(ng,nb,keo,v_mat,r_mat)
        integer :: nb, ng
        real(double), dimension(ng) :: keo
        complex(double), dimension(ng,nb) :: r_mat, v_mat
        real(double), parameter :: min_avg_ke = 0.01_double

        integer :: ib, ig
        real(double) :: x
        real(double), dimension(:), allocatable :: avg_ke_local, avg_ke

        allocate( avg_ke_local(nb), avg_ke(nb) )
        do ib = 1,nb
          avg_ke_local(ib) = 0.0_double
          do ig = 1,ng
            avg_ke_local(ib) = avg_ke_local(ib) + keo(ig)*conjg(v_mat(ig,ib))*v_mat(ig,ib)
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,avg_ke_local,avg_ke)
        do ib = 1,nb
          avg_ke(ib) = max(avg_ke(ib),min_avg_ke)
        end do
        do ib = 1,nb
          do ig = 1,ng
            x = keo(ig)/(1.5_double*avg_ke(ib))
            r_mat(ig,ib) = r_mat(ig,ib)*(1.0_double + x*(1.0_double + x*(1.0_double + x)))/ &
                                (1.0_double + x*(1.0_double + x*(1.0_double + x*(1.0_double + 3.0_double*x))))
          end do
        end do
        deallocate( avg_ke_local, avg_ke )

      end subroutine
      
      subroutine kernel_kinetic_energy_i(ng,nb,keo,v_mat,wts,ke)
        integer :: nb, ng
        real(double) :: ke
        real(double), dimension(ng) :: keo
        real(double), dimension(nb) :: wts
        complex(double), dimension(ng,nb) :: v_mat

        integer :: ib, ig
        real(double) :: ke_local, sum_ke

        ke_local = 0.0_double
        do ib = 1,nb
          sum_ke = 0.0_double
          do ig = 1,ng
            sum_ke = sum_ke + keo(ig)*conjg(v_mat(ig,ib))*v_mat(ig,ib)
          end do
          ke_local = ke_local + wts(ib)*sum_ke
        end do
        call allreduce(KGROUP,MPI_SUM,ke_local,ke)

      end subroutine
      
      subroutine projector_pressure_fs_i(v,hr,wts,p)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        real(double), dimension(:), intent(in) :: wts
        real(double), intent(out) :: p
!       requires:  Projector m values for a given l value be in increasing order.
!       effects: Returns the contribution to the pressure due to the  Fourier-space projectors.
!       errors: Passes errors.
!       documents: ncp_projector_stress.

        character(1) :: transa, transb
        integer :: ia, ib, ig, ip, it, itp
        integer :: na, nap, nb, ng, np, nt, ntp
        integer :: l, last_l, m, last_m
        real(double) :: cell_volume, r0, rsh, sr_icv
        real(double), dimension(3) :: pos
        real(double), dimension(:), allocatable :: gkm
        real(double), dimension(:,:), allocatable :: gk, psfv
        complex(double), parameter :: i = (0.0_double,1.0_double)
        complex(double) :: alpha, beta
        complex(double), dimension(:,:), allocatable :: pdots, pdots_local, d_vnl
        complex(double), dimension(:,:), allocatable :: d_pdots, d_pdots_local
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb
        type(type_info), dimension(:), pointer :: atom_type
        type(atom_info), dimension(:), pointer :: atom

        call my(v)

        nullify( worm_v )
        nullify( worm_mb )
        nullify( atom_type )
        nullify( atom )

        p = 0.0_double

        if (x_n_projectors(hr%hk%o%hc%o%ao) /= 0) then

          worm_v => wormhole(v)
          worm_mb => wormhole(hr%hk%o%mb)

          np = x_n_projectors(hr%hk%o%hc%o%ao)
          nt = x_n_types(hr%hk%o%hc%o%ao)
          na = x_n_atoms(hr%hk%o%hc%o%ao)
          ng = size(worm_v%mat,1)
          nb = size(worm_v%mat,2)

          cell_volume = x_cell_volume(x_lattice(hr%hk%o%hc%o%crystal))
          sr_icv = sqrt(1.0_double/cell_volume)

          allocate( gk(3,ng), gkm(ng), psfv(2,ng) )
          do ig = 1,ng
            gk(1,ig) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2,ig) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3,ig) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = norm(gk(:,ig))
          end do

          allocate( atom_type(nt) )
          do it = 1,nt
            ntp = x_n_type_projectors(hr%hk%o%hc%o%ao,it) ; if (error()) goto 100
            if (ntp == 0) then
              nullify( atom_type(it)%pff )
              cycle
            else
              allocate( atom_type(it)%pff(ng,ntp) )
            end if
            last_l = -1
            last_m = 0
            do itp = 1,ntp
              l = x_projector_l(hr%hk%o%hc%o%ao,it,itp) ; if (error()) goto 100
              m = x_projector_m(hr%hk%o%hc%o%ao,it,itp) ; if (error()) goto 100
              r0 = real(l,double) + 1.5_double
              !if ((l /= last_l) .or. (m /= last_m + 1)) then
                psfv = projector_stress_f_values(hr%hk%o%hc%o%ao,gkm,it,itp) ; if (error()) goto 100
              !end if
              do ig = 1,ng
                rsh = real_spherical_harmonic(l,m,gk(:,ig))
                atom_type(it)%pff(ig,itp) = (r0*psfv(1,ig) - gkm(ig)**2*psfv(2,ig))*rsh*sr_icv
              end do
              last_l = l
              last_m = m
            end do
          end do
          deallocate( gk, gkm, psfv )

          allocate( pdots(np,nb), pdots_local(np,nb) )
          transa = 'c' ; transb = 'n'
          alpha = cmplx(1,0,double) ; beta = cmplx(0,0,double)
          call start_timer("operators: zgemm")
          call zgemm(transa,transb,np,nb,ng,alpha,hr%vnl,ng,worm_v%mat,ng,beta,pdots_local,np)
          call stop_timer("operators: zgemm")
          call allreduce(KGROUP,MPI_SUM,pdots_local,pdots)
          call atomic_hamiltonian(hr%hk%o%hc%o%apot,pdots)
          deallocate( pdots_local )

          allocate( atom(na) )
          do ia = 1,na
            nap = x_n_atom_projectors(hr%hk%o%hc%o%ao,ia) ; if (error()) goto 100
            if (nap == 0) then
              nullify( atom(ia)%sf )
            else
              allocate( atom(ia)%sf(ng) )
              pos = lat2r(x_lattice(hr%hk%o%hc%o%crystal),x_position(x_atoms(hr%hk%o%hc%o%crystal),ia))
              atom(ia)%sf = exp(i*(pos(1)*worm_mb%gpt(:,1) + pos(2)*worm_mb%gpt(:,2) + pos(3)*worm_mb%gpt(:,3)))
            end if
          end do

          allocate( d_pdots(np,nb), d_pdots_local(np,nb), d_vnl(ng,np) )
          transa = 't' ; transb = 'n'
          alpha = cmplx(1,0,double) ; beta = cmplx(0,0,double)
          do ip = 1,np
            ia = x_projector_atom(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            it = x_projector_type(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            itp = x_projector_index_in_type(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            l = x_projector_l(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            d_vnl(:,ip) = i**l*atom(ia)%sf*atom_type(it)%pff(:,itp)
          end do
          call start_timer("operators: zgemm")
          call zgemm(transa,transb,np,nb,ng,alpha,d_vnl,ng,worm_v%mat,ng,beta,d_pdots_local,np)
          call stop_timer("operators: zgemm")
          call allreduce(KGROUP,MPI_SUM,d_pdots_local,d_pdots)
          deallocate( d_pdots_local, d_vnl )
          do it = 1,nt
            if (associated( atom_type(it)%pff )) deallocate( atom_type(it)%pff )
          end do
          deallocate( atom_type )
          do ia = 1,na
            if (associated( atom(ia)%sf )) deallocate( atom(ia)%sf )
          end do
          deallocate( atom )

          do ib = 1,nb
            pdots(:,ib) = wts(ib)*conjg(pdots(:,ib))
          end do
          p = (2.0_double/(3.0_double*cell_volume))*real(sum(pdots*d_pdots(:,:)),double)
          deallocate( pdots, d_pdots )

        end if

100     if (allocated( gk )) deallocate( gk )
        if (allocated( gkm )) deallocate( gkm )
        if (allocated( psfv )) deallocate( psfv )
        if (allocated( pdots )) deallocate( pdots )
        if (allocated( pdots_local )) deallocate( pdots_local )
        if (allocated( d_pdots )) deallocate( d_pdots )
        if (allocated( d_pdots_local )) deallocate( d_pdots_local )
        if (allocated( d_vnl )) deallocate( d_vnl )
        nullify( worm_v )
        nullify( worm_mb )
        if (associated( atom_type )) then
          do it = 1,nt
            if (associated( atom_type(it)%pff )) deallocate( atom_type(it)%pff )
          end do
          deallocate( atom_type )
        end if
        if (associated( atom )) then
          do ia = 1,na
            if (associated( atom(ia)%sf )) deallocate( atom(ia)%sf )
          end do
          deallocate( atom )
        end if

        call glean(thy(v))

        if (error("Exit operators_mod::projector_pressure_fs_i")) continue

      end subroutine

      subroutine projector_stress_tensor_fs_i(v,hr,wts,s)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        real(double), dimension(:), intent(in) :: wts
        real(double), dimension(:,:), intent(out) :: s
!       effects: Returns contributions to the stress tensor from Fourier-space projectors.
!       errors: Passes errors.
!       documents: ncp_projector_stress.

        character(1) :: transa, transb
        integer :: ia, ib, ic, ig, im, ip, ir, it, itp
        integer :: na, nap, nb, ng, np, nt, ntp
        integer :: l, last_l, m, last_m
        real(double) :: cell_volume, rsh, sr_icv
        real(double), dimension(3) :: grsh, pos
        real(double), dimension(:), allocatable :: gkm
        real(double), dimension(:,:), allocatable :: gk, psfv
        complex(double), parameter :: i = (0.0_double,1.0_double)
        complex(double) :: alpha, beta
        complex(double), dimension(:,:), allocatable :: pdots, pdots_local, d_vnl
        complex(double), dimension(:,:,:), allocatable :: d_pdots, d_pdots_local
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb
        type(type_info), dimension(:), pointer :: atom_type
        type(atom_info), dimension(:), pointer :: atom

        call my(v)

        nullify( worm_v )
        nullify( worm_mb )
        nullify( atom_type )
        nullify( atom )

        s = 0.0_double

        if (x_n_projectors(hr%hk%o%hc%o%ao) /= 0) then

          worm_v => wormhole(v)
          worm_mb => wormhole(hr%hk%o%mb)

          np = x_n_projectors(hr%hk%o%hc%o%ao)
          nt = x_n_types(hr%hk%o%hc%o%ao)
          na = x_n_atoms(hr%hk%o%hc%o%ao)
          ng = size(worm_v%mat,1)
          nb = size(worm_v%mat,2)

          cell_volume = x_cell_volume(x_lattice(hr%hk%o%hc%o%crystal))
          sr_icv = sqrt(1.0_double/cell_volume)

          allocate( gk(3,ng), gkm(ng), psfv(2,ng) )
          do ig = 1,ng
            gk(1,ig) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2,ig) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3,ig) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = norm(gk(:,ig))
          end do

          allocate( atom_type(nt) )
          do it = 1,nt
            ntp = x_n_type_projectors(hr%hk%o%hc%o%ao,it) ; if (error()) goto 100
            if (ntp == 0) then
              nullify( atom_type(it)%sff )
              cycle
            else
              allocate( atom_type(it)%sff(ng,ntp,6) )
            end if
            last_l = -1
            last_m = 0
            do itp = 1,ntp
              l = x_projector_l(hr%hk%o%hc%o%ao,it,itp) ; if (error()) goto 100
              m = x_projector_m(hr%hk%o%hc%o%ao,it,itp) ; if (error()) goto 100
              !if ((l /= last_l) .or. (m /= last_m + 1)) then
                psfv = projector_stress_f_values(hr%hk%o%hc%o%ao,gkm,it,itp) ; if (error()) goto 100
              !end if
              do ig = 1,ng
                rsh = real_spherical_harmonic(l,m,gk(:,ig))
                grsh = grad_real_spherical_harmonic(l,m,gk(:,ig))
                im = 0
                do ic = 1,3
                  do ir = 1,ic
                    im = im + 1
                    atom_type(it)%sff(ig,itp,im) = (gk(ir,ig)*gk(ic,ig)*psfv(2,ig)*rsh - gk(ir,ig)*psfv(1,ig)*grsh(ic))*sr_icv
                    if (ir == ic) atom_type(it)%sff(ig,itp,im) = atom_type(it)%sff(ig,itp,im) - 0.5_double*psfv(1,ig)*rsh*sr_icv
                  end do
                end do
              end do
              last_l = l
              last_m = m
            end do
          end do
          deallocate( gk, gkm, psfv )

          allocate( pdots(np,nb), pdots_local(np,nb) )
          transa = 'c' ; transb = 'n'
          alpha = cmplx(1,0,double) ; beta = cmplx(0,0,double)
          call start_timer("operators: zgemm")
          call zgemm(transa,transb,np,nb,ng,alpha,hr%vnl,ng,worm_v%mat,ng,beta,pdots_local,np)
          call stop_timer("operators: zgemm")
          call allreduce(KGROUP,MPI_SUM,pdots_local,pdots)
          deallocate( pdots_local )
          call atomic_hamiltonian(hr%hk%o%hc%o%apot,pdots)

          allocate( atom(na) )
          do ia = 1,na
            nap = x_n_atom_projectors(hr%hk%o%hc%o%ao,ia) ; if (error()) goto 100
            if (nap == 0) then
              nullify( atom(ia)%sf )
            else
              allocate( atom(ia)%sf(ng) )
              pos = lat2r(x_lattice(hr%hk%o%hc%o%crystal),x_position(x_atoms(hr%hk%o%hc%o%crystal),ia))
              atom(ia)%sf = exp(i*(pos(1)*worm_mb%gpt(:,1) + pos(2)*worm_mb%gpt(:,2) + pos(3)*worm_mb%gpt(:,3)))
            end if
          end do

          allocate( d_pdots(np,nb,6), d_pdots_local(np,nb,6), d_vnl(ng,np) )
          transa = 't' ; transb = 'n'
          alpha = cmplx(1,0,double) ; beta = cmplx(0,0,double)
          im = 0
          do ic = 1,3
            do ir = 1,ic
              im = im + 1
              do ip = 1,np
                ia = x_projector_atom(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
                it = x_projector_type(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
                itp = x_projector_index_in_type(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
                l = x_projector_l(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
                d_vnl(:,ip) = i**l*atom(ia)%sf*atom_type(it)%sff(:,itp,im)
              end do
              call start_timer("operators: zgemm")
              call zgemm(transa,transb,np,nb,ng,alpha,d_vnl,ng,worm_v%mat,ng,beta,d_pdots_local(:,:,im),np)
              call stop_timer("operators: zgemm")
            end do
          end do
          call allreduce(KGROUP,MPI_SUM,d_pdots_local,d_pdots)
          deallocate( d_pdots_local, d_vnl )
          do it = 1,nt
            if (associated( atom_type(it)%sff )) deallocate( atom_type(it)%sff )
          end do
          deallocate( atom_type )
          do ia = 1,na
            if (associated( atom(ia)%sf )) deallocate( atom(ia)%sf )
          end do
          deallocate( atom )

          do ib = 1,nb
            pdots(:,ib) = wts(ib)*conjg(pdots(:,ib))
          end do
          im = 0
          do ic = 1,3
            do ir = 1,ic
              im = im + 1
              s(ir,ic) = (2.0_double/cell_volume)*real(sum(pdots*d_pdots(:,:,im)),double)
              if (ic /= ir) s(ic,ir) = s(ir,ic)
            end do
          end do
          deallocate( pdots, d_pdots )

        end if

100     if (allocated( gk )) deallocate( gk )
        if (allocated( gkm )) deallocate( gkm )
        if (allocated( psfv )) deallocate( psfv )
        if (associated( atom_type )) then
          do it = 1,nt
            if (associated( atom_type(it)%sff )) deallocate( atom_type(it)%sff )
          end do
          deallocate( atom_type )
        end if
        if (associated( atom )) then
          do ia = 1,na
            if (associated( atom(ia)%sf )) deallocate( atom(ia)%sf )
          end do
          deallocate( atom )
        end if
        if (allocated( pdots )) deallocate( pdots )
        if (allocated( pdots_local )) deallocate( pdots_local )
        if (allocated( d_pdots )) deallocate( d_pdots )
        if (allocated( d_pdots_local )) deallocate( d_pdots_local )
        if (allocated( d_vnl )) deallocate( d_vnl )
        nullify( worm_v )
        nullify( worm_mb )

        call glean(thy(v))

        if (error("Exit operators_mod::projector_stress_tensor_fs_i")) continue

      end subroutine

      subroutine expand_fs_projectors_i(hr)
        type(hamiltonian_rep) :: hr
!       requires: hr%vnl be nullified.
!       effects: Forms hr%vnl.
!       errors: Passes errors.
!       documents: ncp_projectors.

        integer :: ia, ic, ip, it, l, last_ia, np, nw
        real(double), dimension(3) :: pos
        complex(double), parameter :: i = (0.0_double,1.0_double)
        complex(double), dimension(:), allocatable :: phase
        type(multibasis_rep), pointer :: worm_mb

        nullify( worm_mb )

        if (x_n_projectors(hr%hk%o%hc%o%ao) /= 0) then

          worm_mb => wormhole(hr%hk%o%mb)
          np = x_n_projectors(hr%hk%o%hc%o%ao)
          nw = size(worm_mb%gpt,1)
          allocate( hr%vnl(nw,np) )

          allocate( phase(nw) )
          last_ia = -1
          do ip = 1,np
            ia = x_projector_atom(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            it = x_projector_type(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            ic = x_projector_index_in_type(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            l = x_projector_l(hr%hk%o%hc%o%ao,ip) ; if (error()) goto 100
            if (ia /= last_ia) then
              pos = lat2r(x_lattice(hr%hk%o%hc%o%crystal),x_position(x_atoms(hr%hk%o%hc%o%crystal),ia))
              phase = exp( -i*(pos(1)*worm_mb%gpt(:,1) + pos(2)*worm_mb%gpt(:,2) + pos(3)*worm_mb%gpt(:,3)) )
            end if
            hr%vnl(:,ip) = (-i)**l*phase*hr%hk%o%type(it)%form_factor(:,ic)
            last_ia = ia
          end do

        end if

100     if (allocated( phase )) deallocate( phase )
        nullify( worm_mb )

        if (error("Exit operators_mod::expand_fs_projectors_i")) continue

      end subroutine

      subroutine form_rs_global_phases_i(hr)
        type(hamiltonian_rep) :: hr
!       requires: hr%phase be nullified.
!       effects: Forms hr%phase.
!       errors: Passes errors.

        real(double), dimension(3) :: k
        real(double), dimension(:,:,:), pointer :: x, y, z
        complex(double), parameter :: i = (0.0_double,1.0_double)
        type(multibasis_rep), pointer :: worm_mb

        nullify( x )
        nullify( y )
        nullify( z )
        nullify( worm_mb )

        if (x_n_projectors(hr%hk%o%hc%o%ao) /= 0) then

          worm_mb => wormhole(hr%hk%o%mb)
          k = worm_mb%kpt
          nullify( worm_mb )

          call mesh(x,y,z,hr%hk%o%hc%o%layout,S_TYPE)
          call alloc(hr%phase,hr%hk%o%hc%o%layout,S_TYPE)
          hr%phase = exp(i*(k(1)*x + k(2)*y + k(3)*z))

        end if

100     if (associated( x )) deallocate( x )
        if (associated( y )) deallocate( y )
        if (associated( z )) deallocate( z )

        if (error("Exit operators_mod::form_rs_global_phases_i")) continue

      end subroutine

      subroutine form_pdots_i(v,hr,pdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        complex(double), dimension(:,:), pointer :: pdots
!       requires: At least one projector exist.

        call my(v)

        select case (hr%hk%o%hc%o%projector_type)
        case (FOURIER_SPACE)
          call form_fs_pdots_i(v,hr,pdots)
        case (REAL_SPACE)
          call form_rs_pdots_i(v,hr,pdots)
        end select

        call glean(thy(v))

        if (error("Exit operators_mod::form_pdots_i")) continue

      end subroutine

      subroutine form_fs_pdots_i(v,hr,pdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        complex(double), dimension(:,:), pointer :: pdots
!       requires: At least one projector exist.

        character(1), parameter :: transa = 'c', transb = 'n'
        integer :: nb, nb1, nb2, nbs, np, nw
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: pdots_local
        type(multivector_rep), pointer :: worm_v

        call my(v)

        nb = x_n_bands(v)
        nw = size(hr%vnl,1)
        np = size(hr%vnl,2)

        worm_v => wormhole(v)

        if (associated( pdots )) deallocate( pdots )
        allocate( pdots(np,nb) )

        nbs = hr%hk%o%hc%o%pbs
        allocate( pdots_local(np,nbs) )
        nb2 = 0
        do while (nb2 < nb)
          nb1 = nb2 + 1
          nb2 = nb2 + nbs
          if (nb2 > nb) then
            nb2 = nb
            nbs = nb2 - nb1 + 1
            deallocate( pdots_local )
            allocate( pdots_local(np,nbs) )
          end if
          call start_timer("operators: zgemm")
          call zgemm(transa,transb,np,nbs,nw,alpha,hr%vnl,nw,worm_v%mat(1,nb1),nw,beta,pdots_local,np)
          call stop_timer("operators: zgemm")
          call allreduce(KGROUP,MPI_SUM,pdots_local,pdots(:,nb1:nb2))
        end do

100     if (allocated( pdots_local )) deallocate( pdots_local )

        nullify( worm_v )

        call glean(thy(v))

        if (error("Exit operators_mod::form_fs_pdots_i")) continue

      end subroutine

      subroutine form_rs_pdots_i(v,hr,all_pdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        complex(double), dimension(:,:), pointer :: all_pdots
!       requires: At least one projector exist.

        integer :: grp, nb, ng, ngt, np, npr, nr, ns, nsr
        integer, dimension(3) :: nd
        real(double) :: sqrt_mesh_volume
        complex(double), dimension(:), allocatable :: band_pdots, srdots, vec
        complex(double), dimension(:,:), allocatable :: mat, pdots
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        call my(v)

        worm_v => wormhole(v)
        worm_mb => wormhole(worm_v%mb)

        hcr => hr%hk%o%hc%o
        hkr => hr%hk%o

        nd = x_dims(hcr%layout)
        sqrt_mesh_volume = sqrt(x_cell_volume(x_lattice(hcr%crystal)))/real(product(nd),double)

        wf1 => datalink(worm_mb%isplan1)

        ns = size(hcr%site)
        nr = size(hcr%site(ns)%replica)
        nsr = hcr%site(ns)%replica(nr)%sr_index
        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)
        ngt = size(worm_mb%gridmap,2)
        np = size(hcr%rsi)
        npr = mpi_nprocs(KGROUP)

        allocate( band_pdots(ns) )
        allocate( srdots(nsr) )
        allocate( vec(ngt) )
        allocate( mat(ng,npr) )
        allocate( pdots(ns,nb) )

        pdots = (0.0_double,0.0_double)
        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_v%mat)
          call band_remap(hkr%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call form_srdots_i(nd,wf1,hr%phase,hcr%rsc1,hcr%rsc2,np,hcr%rsi,hcr%rsp,nsr,srdots)
            call combine_srdots_i(srdots,band_pdots,hcr,hkr)
            pdots(:,worm_mb%my_band(grp)) = band_pdots
          end if
        end do
        pdots = sqrt_mesh_volume*pdots

        if (associated( all_pdots )) deallocate( all_pdots )
        allocate( all_pdots(ns,nb) )
        call allreduce(KGROUP,MPI_SUM,pdots,all_pdots)

100     if (allocated( band_pdots )) deallocate( band_pdots )
        if (allocated( srdots )) deallocate( srdots )
        if (allocated( vec )) deallocate( vec )
        if (allocated( pdots )) deallocate( pdots )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )
        nullify( worm_v )
        nullify( worm_mb )
        nullify( hcr )
        nullify( hkr )

        call glean(thy(v))

        if (error("Exit operators_mod::form_rs_pdots_i")) continue

      end subroutine

      subroutine form_grad_pdots_i(v,hr,ic,gpdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        integer :: ic
        complex(double), dimension(:,:), pointer :: gpdots
!       requires: At least one projector exist.

        call my(v)

        select case (hr%hk%o%hc%o%projector_type)
        case (FOURIER_SPACE)
          call form_fs_grad_pdots_i(v,hr,ic,gpdots)
        case (REAL_SPACE)
          call form_rs_grad_pdots_i(v,hr,ic,gpdots)
        end select

        call glean(thy(v))

        if (error("Exit operators_mod::form_grad_pdots_i")) continue

      end subroutine

      subroutine form_fs_grad_pdots_i(v,hr,ic,gpdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        integer :: ic
        complex(double), dimension(:,:), pointer :: gpdots
!       requires: At least one projector exist.

        character(1), parameter :: transa = 'c', transb = 'n'
        integer :: ib, nb, nb1, nb2, nbs, np, nw
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), parameter :: i = (0.0_double,1.0_double)
        complex(double), dimension(:,:), allocatable :: gmat, gpdots_local
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb

        call my(v)

        nb = x_n_bands(v)
        nw = size(hr%vnl,1)
        np = size(hr%vnl,2)

        worm_v => wormhole(v)
        worm_mb => wormhole(hr%hk%o%mb)

        allocate( gmat(nw,nb) )
        do ib = 1,nb
          gmat(:,ib) = i*worm_mb%gpt(:,ic)*worm_v%mat(:,ib)
        end do

        if (associated( gpdots )) deallocate( gpdots )
        allocate( gpdots(np,nb) )

        nbs = hr%hk%o%hc%o%pbs
        allocate( gpdots_local(np,nbs) )
        nb2 = 0
        do while (nb2 < nb)
          nb1 = nb2 + 1
          nb2 = nb2 + nbs
          if (nb2 > nb) then
            nb2 = nb
            nbs = nb2 - nb1 + 1
            deallocate( gpdots_local )
            allocate( gpdots_local(np,nbs) )
          end if
          call start_timer("operators: zgemm")
          call zgemm(transa,transb,np,nbs,nw,alpha,hr%vnl,nw,gmat(1,nb1),nw,beta,gpdots_local,np)
          call stop_timer("operators: zgemm")
          call allreduce(KGROUP,MPI_SUM,gpdots_local,gpdots(:,nb1:nb2))
        end do

100     if (allocated( gmat )) deallocate( gmat )
        if (allocated( gpdots_local )) deallocate( gpdots_local )

        nullify( worm_v )
        nullify( worm_mb )

        call glean(thy(v))

        if (error("Exit operators_mod::form_fs_grad_pdots_i")) continue

      end subroutine

      subroutine form_rs_grad_pdots_i(v,hr,ic,all_gpdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        integer :: ic
        complex(double), dimension(:,:), pointer :: all_gpdots
!       requires: At least one projector exist.

        integer :: grp, nb, ng, ngt, np, npr, nr, ns, nsr
        integer, dimension(3) :: nd
        real(double) :: sqrt_mesh_volume
        complex(double), dimension(:), allocatable :: band_gpdots, gsrdots, vec
        complex(double), dimension(:,:), allocatable :: gpdots, mat
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        call my(v)

        worm_v => wormhole(v)
        worm_mb => wormhole(worm_v%mb)

        hcr => hr%hk%o%hc%o
        hkr => hr%hk%o

        nd = x_dims(hcr%layout)
        sqrt_mesh_volume = sqrt(x_cell_volume(x_lattice(hcr%crystal)))/real(product(nd),double)

        wf1 => datalink(worm_mb%isplan1)

        ns = size(hcr%site)
        nr = size(hcr%site(ns)%replica)
        nsr = hcr%site(ns)%replica(nr)%sr_index
        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)
        ngt = size(worm_mb%gridmap,2)
        np = size(hcr%rsi)
        npr = mpi_nprocs(KGROUP)

        allocate( band_gpdots(ns) )
        allocate( gsrdots(nsr) )
        allocate( vec(ngt) )
        allocate( gpdots(ns,nb) )
        allocate( mat(ng,npr) )

        gpdots = (0.0_double,0.0_double)
        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_v%mat)
          call band_remap(hkr%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call form_gsrdots_i(nd,wf1,hr%phase,hcr%rsc1,hcr%rsc2,np,hcr%rsi,ic,hcr%grsp,nsr,gsrdots)
            call combine_gsrdots_1_i(gsrdots,band_gpdots,hcr,hkr)
            gpdots(:,worm_mb%my_band(grp)) = band_gpdots
          end if
        end do
        gpdots = sqrt_mesh_volume*gpdots

        if (associated( all_gpdots )) deallocate( all_gpdots )
        allocate( all_gpdots(ns,nb) )
        call allreduce(KGROUP,MPI_SUM,gpdots,all_gpdots)

100     if (allocated( band_gpdots )) deallocate( band_gpdots )
        if (allocated( gsrdots )) deallocate( gsrdots )
        if (allocated( vec )) deallocate( vec )
        if (allocated( gpdots )) deallocate( gpdots )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )
        nullify( worm_v )
        nullify( worm_mb )
        nullify( hcr )
        nullify( hkr )

        call glean(thy(v))

        if (error("Exit operators_mod::form_rs_grad_pdots_i")) continue

      end subroutine

      subroutine form_pdots_and_grad_pdots_i(v,hr,pdots,gpdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        complex(double), dimension(:,:), pointer :: pdots
        complex(double), dimension(:,:,:), pointer :: gpdots
!       requires: At least one projector exist.

        call my(v)

        select case (hr%hk%o%hc%o%projector_type)
        case (FOURIER_SPACE)
          call form_fs_pdots_and_grad_pdots_i(v,hr,pdots,gpdots)
        case (REAL_SPACE)
          call form_rs_pdots_and_grad_pdots_i(v,hr,pdots,gpdots)
        end select

        call glean(thy(v))

        if (error("Exit operators_mod::form_pdots_and_grad_pdots_i")) continue

      end subroutine

      subroutine form_fs_pdots_and_grad_pdots_i(v,hr,pdots,gpdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        complex(double), dimension(:,:), pointer :: pdots
        complex(double), dimension(:,:,:), pointer :: gpdots
!       requires: At least one projector exist.

        character(1), parameter :: transa = 'c', transb = 'n'
        integer :: ib, ic, nb, np, nw
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), parameter :: i = (0.0_double,1.0_double)
        complex(double), dimension(:,:), allocatable :: gmat, pdots_local
        complex(double), dimension(:,:,:), allocatable :: gpdots_local
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb

        call my(v)

        nb = x_n_bands(v)
        nw = size(hr%vnl,1)
        np = size(hr%vnl,2)

        worm_v => wormhole(v)
        worm_mb => wormhole(hr%hk%o%mb)

        if (associated( pdots )) deallocate( pdots )
        allocate( pdots(np,nb) )

        allocate( pdots_local(np,nb) )
        call start_timer("operators: zgemm")
        call zgemm(transa,transb,np,nb,nw,alpha,hr%vnl,nw,worm_v%mat,nw,beta,pdots_local,np)
        call stop_timer("operators: zgemm")
        call allreduce(KGROUP,MPI_SUM,pdots_local,pdots)
        deallocate( pdots_local )

        if (associated( gpdots )) deallocate( gpdots )
        allocate( gpdots(np,nb,3) )

        allocate( gmat(nw,nb) )
        allocate( gpdots_local(np,nb,3) )
        do ic = 1,3
          do ib = 1,nb
            gmat(:,ib) = i*worm_mb%gpt(:,ic)*worm_v%mat(:,ib)
          end do
          call start_timer("operators: zgemm")
          call zgemm(transa,transb,np,nb,nw,alpha,hr%vnl,nw,gmat,nw,beta,gpdots_local(1,1,ic),np)
          call stop_timer("operators: zgemm")
        end do
        call allreduce(KGROUP,MPI_SUM,gpdots_local,gpdots)

100     if (allocated( gmat )) deallocate( gmat )
        if (allocated( pdots_local )) deallocate( pdots_local )
        if (allocated( gpdots_local )) deallocate( gpdots_local )

        nullify( worm_v )
        nullify( worm_mb )

        call glean(thy(v))

        if (error("Exit operators_mod::form_fs_pdots_and_grad_pdots_i")) continue

      end subroutine

      subroutine form_rs_pdots_and_grad_pdots_i(v,hr,all_pdots,all_gpdots)
        type(multivector_obj) :: v
        type(hamiltonian_rep) :: hr
        complex(double), dimension(:,:), pointer :: all_pdots
        complex(double), dimension(:,:,:), pointer :: all_gpdots
!       requires: At least one projector exist.

        integer :: grp, nb, ng, ngt, np, npr, nr, ns, nsr
        integer, dimension(3) :: nd
        real(double) :: sqrt_mesh_volume
        complex(double), dimension(:), allocatable :: band_pdots, srdots, vec
        complex(double), dimension(:,:), allocatable :: band_gpdots, gsrdots, mat, pdots
        complex(double), dimension(:,:,:), allocatable :: gpdots
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multivector_rep), pointer :: worm_v
        type(multibasis_rep), pointer :: worm_mb
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        call my(v)

        worm_v => wormhole(v)
        worm_mb => wormhole(worm_v%mb)

        hcr => hr%hk%o%hc%o
        hkr => hr%hk%o

        nd = x_dims(hcr%layout)
        sqrt_mesh_volume = sqrt(x_cell_volume(x_lattice(hcr%crystal)))/real(product(nd),double)

        wf1 => datalink(worm_mb%isplan1)

        ns = size(hcr%site)
        nr = size(hcr%site(ns)%replica)
        nsr = hcr%site(ns)%replica(nr)%sr_index
        ng = size(worm_v%mat,1)
        nb = size(worm_v%mat,2)
        ngt = size(worm_mb%gridmap,2)
        np = size(hcr%rsi)
        npr = mpi_nprocs(KGROUP)

        if (allocated( pdots )) deallocate( pdots )
        allocate( pdots(ns,nb) )

        if (allocated( gpdots )) deallocate( gpdots )
        allocate( gpdots(ns,nb,3) )

        allocate( srdots(nsr), band_pdots(ns) )
        allocate( gsrdots(3,nsr), band_gpdots(ns,3) )
        allocate( vec(ngt), mat(ng,npr) )

        pdots = (0.0_double,0.0_double)
        gpdots = (0.0_double,0.0_double)
        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,worm_v%mat)
          call band_remap(hkr%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call form_srdots_and_gsrdots_i(nd,wf1,hr%phase,hcr%rsc1,hcr%rsc2,np,hcr%rsi,hcr%rsp,hcr%grsp,nsr,srdots,gsrdots)
            call combine_srdots_i(srdots,band_pdots,hcr,hkr)
            call combine_gsrdots_3_i(gsrdots,band_gpdots,hcr,hkr)
            pdots(:,worm_mb%my_band(grp)) = band_pdots
            gpdots(:,worm_mb%my_band(grp),:) = band_gpdots
          end if
        end do
        pdots = sqrt_mesh_volume*pdots
        gpdots = sqrt_mesh_volume*gpdots

        if (associated( all_pdots )) deallocate( all_pdots )
        allocate( all_pdots(ns,nb) )

        if (associated( all_gpdots )) deallocate( all_gpdots )
        allocate( all_gpdots(ns,nb,3) )

        call allreduce(KGROUP,MPI_SUM,pdots,all_pdots)
        call allreduce(KGROUP,MPI_SUM,gpdots,all_gpdots)

100     if (allocated( band_pdots )) deallocate( band_pdots )
        if (allocated( srdots )) deallocate( srdots )
        if (allocated( vec )) deallocate( vec )
        if (allocated( band_gpdots )) deallocate( band_gpdots )
        if (allocated( gsrdots )) deallocate( gsrdots )
        if (allocated( pdots )) deallocate( pdots )
        if (allocated( mat )) deallocate( mat )
        if (allocated( gpdots )) deallocate( gpdots )
        nullify( wf1 )
        nullify( worm_v )
        nullify( worm_mb )
        nullify( hcr )
        nullify( hkr )

        call glean(thy(v))

        if (error("Exit operators_mod::form_rs_pdots_and_grad_pdots_i")) continue

      end subroutine

      subroutine form_srdots_i(nd,wf,phase,rsc1,rsc2,np,rsi,rsp,nsr,srdots)
        integer :: np, nsr
        integer, dimension(3) :: nd
        integer, dimension(np) :: rsi
        integer, dimension(nd(1),nd(2),nd(3)) :: rsc1, rsc2
        real(double), dimension(np) :: rsp
        complex(double), dimension(nsr) :: srdots
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf, phase

        integer :: i1, i2, i3, ip
        complex(double) :: c1

        srdots = (0.0_double,0.0_double)
        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              c1 = wf(i1,i2,i3)*phase(i1,i2,i3)
              do ip = rsc1(i1,i2,i3),rsc2(i1,i2,i3)
                srdots(rsi(ip)) = srdots(rsi(ip)) + rsp(ip)*c1
              end do
            end do
          end do
        end do

      end subroutine

      subroutine form_srdots_and_gsrdots_i(nd,wf,phase,rsc1,rsc2,np,rsi,rsp,grsp,nsr,srdots,gsrdots)
        integer :: np, nsr
        integer, dimension(3) :: nd
        integer, dimension(np) :: rsi
        integer, dimension(nd(1),nd(2),nd(3)) :: rsc1, rsc2
        real(double), dimension(np) :: rsp
        real(double), dimension(3,np) :: grsp
        complex(double), dimension(nsr) :: srdots
        complex(double), dimension(3,nsr) :: gsrdots
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf, phase

        integer :: i1, i2, i3, ip
        complex(double) :: c1

        srdots = (0.0_double,0.0_double)
        gsrdots = (0.0_double,0.0_double)
        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              c1 = wf(i1,i2,i3)*phase(i1,i2,i3)
              do ip = rsc1(i1,i2,i3),rsc2(i1,i2,i3)
                srdots(rsi(ip)) = srdots(rsi(ip)) + rsp(ip)*c1
                gsrdots(:,rsi(ip)) = gsrdots(:,rsi(ip)) + grsp(:,ip)*c1
              end do
            end do
          end do
        end do

      end subroutine

      subroutine distribute_srdots_i(nd,wf,phase,rsc1,rsc2,np,rsi,rsp,nsr,srdots)
        integer :: np, nsr
        integer, dimension(3) :: nd
        integer, dimension(np) :: rsi
        integer, dimension(nd(1),nd(2),nd(3)) :: rsc1, rsc2
        real(double), dimension(np) :: rsp
        complex(double), dimension(nsr) :: srdots
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf, phase

        integer :: i1, i2, i3, ip
        complex(double) :: c1

        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              c1 = (0.0_double,0.0_double)
              do ip = rsc1(i1,i2,i3),rsc2(i1,i2,i3)
                c1 = c1 + srdots(rsi(ip))*rsp(ip)
              end do
              wf(i1,i2,i3) = c1*conjg(phase(i1,i2,i3))
            end do
          end do
        end do

      end subroutine

      subroutine form_gsrdots_i(nd,wf,phase,rsc1,rsc2,np,rsi,ic,grsp,nsr,gsrdots)
        integer :: np, ic, nsr
        integer, dimension(3) :: nd
        integer, dimension(np) :: rsi
        integer, dimension(nd(1),nd(2),nd(3)) :: rsc1, rsc2
        real(double), dimension(3,np) :: grsp
        complex(double), dimension(nsr) :: gsrdots
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf, phase

        integer :: i1, i2, i3, ip
        complex(double) :: c1

        gsrdots = (0.0_double,0.0_double)
        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              c1 = wf(i1,i2,i3)*phase(i1,i2,i3)
              do ip = rsc1(i1,i2,i3),rsc2(i1,i2,i3)
                gsrdots(rsi(ip)) = gsrdots(rsi(ip)) + grsp(ic,ip)*c1
              end do
            end do
          end do
        end do

      end subroutine

      subroutine distribute_srdots_and_lp_i(nd,wf,phase,lp,rsc1,rsc2,np,rsi,rsp,nsr,srdots)
        integer :: np, nsr
        integer, dimension(3) :: nd
        integer, dimension(np) :: rsi
        integer, dimension(nd(1),nd(2),nd(3)) :: rsc1, rsc2
        real(double), dimension(np) :: rsp
        real(double), dimension(nd(1),nd(2),nd(3)) :: lp
        complex(double), dimension(nsr) :: srdots
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf, phase

        integer :: i1, i2, i3, ip
        complex(double) :: c1

        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              c1 = (0.0_double,0.0_double)
              do ip = rsc1(i1,i2,i3),rsc2(i1,i2,i3)
                c1 = c1 + srdots(rsi(ip))*rsp(ip)
              end do
              wf(i1,i2,i3) = wf(i1,i2,i3)*lp(i1,i2,i3) + c1*conjg(phase(i1,i2,i3))
            end do
          end do
        end do

      end subroutine

      subroutine combine_srdots_i(srdots,pdots,hcr,hkr)
        complex(double), dimension(:) :: srdots, pdots
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        integer :: ir, is

        do is = 1,size(hcr%site)
          pdots(is) = (0.0_double,0.0_double)
          do ir = 1,size(hcr%site(is)%replica)
            pdots(is) = pdots(is) + srdots(hcr%site(is)%replica(ir)%sr_index)*hkr%site(is)%replica(ir)%phase
          end do
        end do

      end subroutine

      subroutine combine_gsrdots_3_i(gsrdots,gpdots,hcr,hkr)
        complex(double), dimension(:,:) :: gsrdots, gpdots
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        integer :: ir, is

        do is = 1,size(hcr%site)
          gpdots(is,:) = (0.0_double,0.0_double)
          do ir = 1,size(hcr%site(is)%replica)
            gpdots(is,:) = gpdots(is,:) + gsrdots(:,hcr%site(is)%replica(ir)%sr_index)*hkr%site(is)%replica(ir)%phase
          end do
        end do

      end subroutine

      subroutine combine_gsrdots_1_i(gsrdots,gpdots,hcr,hkr)
        complex(double), dimension(:) :: gsrdots, gpdots
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        integer :: ir, is

        do is = 1,size(hcr%site)
          gpdots(is) = (0.0_double,0.0_double)
          do ir = 1,size(hcr%site(is)%replica)
            gpdots(is) = gpdots(is) + gsrdots(hcr%site(is)%replica(ir)%sr_index)*hkr%site(is)%replica(ir)%phase
          end do
        end do

      end subroutine

      subroutine disperse_pdots_i(srdots,pdots,hcr,hkr)
        complex(double), dimension(:) :: srdots, pdots
        type(h_common_rep), pointer :: hcr
        type(h_kpoint_rep), pointer :: hkr

        integer :: ir, is

        do is = 1,size(hcr%site)
          do ir = 1,size(hcr%site(is)%replica)
            srdots(hcr%site(is)%replica(ir)%sr_index) = pdots(is)*conjg(hkr%site(is)%replica(ir)%phase)
          end do
        end do

      end subroutine

      subroutine gather_wf_i(nd,wf,ng,map,wfv)
        integer :: ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        complex(double), dimension(ng) :: wfv
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf

        integer :: ig

        do ig = 1,ng
          wfv(ig) = wf(map(1,ig),map(2,ig),map(3,ig))
        end do

      end subroutine

      subroutine scatter_wf_i(nd,wf,ng,map,wfv)
        integer :: ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        complex(double), dimension(ng) :: wfv
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf

        integer :: i1, i2, i3, ig

        do i3 = 1,nd(3)
          do i2 = 1,nd(2)
            do i1 = 1,nd(1)
              wf(i1,i2,i3) = (0.0_double,0.0_double)
            end do
          end do
        end do
        do ig = 1,ng
          wf(map(1,ig),map(2,ig),map(3,ig)) = wfv(ig)
        end do

      end subroutine

      subroutine extract_vectors_i(fb,lb,ng,np,mat,nb,bmat)
        integer :: fb, lb, nb, ng, np
        complex(double), dimension(ng,np) :: mat
        complex(double), dimension(ng,nb) :: bmat

        integer :: ib, ig, ip

        ip = 0
        do ib = fb,lb
          ip = ip + 1
          do ig = 1,ng
            mat(ig,ip) = bmat(ig,ib)
          end do
        end do
        if (np > nb) then
          do ip = nb+1,np
            do ig = 1,ng
              mat(ig,ip) = (0.0_double,0.0_double)
            end do
          end do
        end if

      end subroutine

      subroutine augment_vectors_i(fb,lb,ng,np,mat,nb,bmat)
        integer :: fb, lb, nb, ng, np
        complex(double), dimension(ng,np) :: mat
        complex(double), dimension(ng,nb) :: bmat

        integer :: ib, ig, ip

        ip = 0
        do ib = fb,lb
          ip = ip + 1
          do ig = 1,ng
            bmat(ig,ib) = bmat(ig,ib) + mat(ig,ip)
          end do
        end do

      end subroutine

      end module 

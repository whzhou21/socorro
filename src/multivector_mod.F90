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

      module multivector_mod
!doc$ module multivector_mod

      use arg_mod
      use error_mod
      use fft_mod
      use ghost_mod
      use grid_mod
      use io_mod
      use kind_mod
      use lattice_mod
      use layout_mod
      use math_mod
      use mpi_mod
      use multibasis_mod
      use point_blas_mod
      use tagio_mod
      use timing_mod
      use xc_type_mod

!     One datatype is available here: type(multivector_obj).

!     multivector_mod encapsulates the set of wave functions at a specific sampling point (k-point).

!     Notes: type(multivector_obj) is not lazy-copied meaning that each multivector_obj has a distinct
!            multivector_rep. This allows the user to mutate %mat data in a wormhole without fear of
!            something counterintuitive hapening. As such the user is responsible for changing the ghost
!            after altering the multivector_rep data via a wormhole.

!cod$
      implicit none
      private

      ! usage
      integer, parameter :: NORMAL    = 1
      integer, parameter :: AUXILIARY = 2

      type, public :: multivector_rep
        type(ghost) :: g
        type(ghost) :: g_storage
        integer :: usage                                  ! usage wrt kgroup
        type(multibasis_obj) :: mb                        ! multibasis
        complex(double), dimension(:,:), pointer :: mat   ! multivector plane-wave coefficients
      end type

      type, public :: multivector_obj
        integer :: ref
        type(multivector_rep), pointer :: o
      end type

!doc$
      public :: multivector
      public :: update
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: wormhole
      public :: x_ref
      public :: x_ghost
      public :: x_multibasis
      public :: x_n_bands
      public :: x_n_gvectors
      public :: put
      public :: portion
      public :: combine
      public :: transform
      public :: project
      public :: orthonormalize
      public :: residual
      public :: decompose
      public :: multiply
      public :: overlap
      public :: apply_filter
      public :: apply_field
      public :: add_grid_density
      public :: distribute
      public :: release
      public :: exx_energy
      public :: exx_derivative
      public :: exx_energy_and_derivative
      public :: diary
      public :: write_restart

!cod$
      interface multivector
        module procedure constructor_mv
      end interface
      interface update
        module procedure update_mv
      end interface
      interface my
        module procedure my_mv, my_new_mv
      end interface
      interface thy
        module procedure thy_mv
      end interface
      interface glean
        module procedure glean_mv
      end interface
      interface bequeath
        module procedure bequeath_mv
      end interface
      interface assignment(=)
        module procedure assign_mv
      end interface
      interface wormhole
        module procedure wormhole_mv
      end interface
      interface x_ref
        module procedure mv_ref
      end interface
      interface x_ghost
        module procedure mv_ghost
      end interface
      interface x_multibasis
        module procedure mv_multibasis
      end interface
      interface x_n_bands
        module procedure mv_n_bands
      end interface
      interface x_n_gvectors
        module procedure mv_n_gvectors
      end interface
      interface put
        module procedure put_mv
      end interface
      interface portion
        module procedure portion_mv_real_scalar, portion_mv_complex_scalar, &
                         portion_mv_real_vector, portion_mv_complex_vector
      end interface
      interface combine
        module procedure combine_2mv_real_scalar, combine_2mv_complex_scalar, &
                         combine_2mv_real_vector, combine_2mv_complex_vector, &
                         combine_3mv_real_vector, combine_3mv_complex_vector
      end interface
      interface transform
        module procedure transform_mv_complex, transform_2mv_complex
      end interface
      interface project
        module procedure project_2mv, project_2mv_mat
      end interface
      interface orthonormalize
        module procedure orthonormalize_mv
      end interface
      interface residual
        module procedure residual_mv
      end interface
      interface decompose
        module procedure decompose_mv
      end interface
      interface multiply
        module procedure multiply_mv_real, multiply_2mv_real, multiply_2mv_complex
      end interface
      interface overlap
        module procedure overlap_mv, overlap_2mv
      end interface
      interface apply_filter
        module procedure filter_mv
      end interface
      interface apply_field
        module procedure scale_mv
      end interface
      interface add_grid_density
        module procedure add_grid_density_mv, add_grid_density_2mv
      end interface
      interface distribute
        module procedure distribute_mv
      end interface
      interface release
        module procedure release_mv
      end interface
      interface exx_energy
        module procedure exx_energy_mv, exx_energy_2mv
      end interface
      interface exx_derivative
        module procedure exx_derivative_mv, exx_derivative_2mv
      end interface
      interface exx_energy_and_derivative
        module procedure exx_energy_and_derivative_mv, exx_energy_and_derivative_2mv
      end interface
      interface diary
        module procedure diary_mv
      end interface
      interface write_restart
        module procedure write_restart_mv
      end interface

 contains

! public routines

      function constructor_mv(mb,init,restf) result(mv)
!doc$ function multivector(mb,init,restf) result(mv)
        type(multibasis_obj) :: mb
        character(line_len), intent(in), optional :: init
        type(tagio_obj), optional :: restf
        type(multivector_obj) :: mv
!       requires: init and restf not both be present; init be appropriate.
!       effects: Constructs a new mv.
!       errors: init not recognized. Problems reading restf.  Mismatched parameters. Passes errors

!cod$
        character(1) :: tios
        integer :: ib, ig, msg, nb, nsg, r_nb, r_ngt
        integer, dimension(:), pointer :: vmap
        integer(long) :: dsize, iosl, ndata, s4
        real(double), parameter :: tol_nbhd = 1.0e-10_double
        real(double) :: r_cutoff
        real(double), dimension(3) :: kpt
        complex(double), dimension(:), allocatable :: c1
        type(multibasis_rep), pointer :: worm_mb

        call my(mb)
        if (present(restf)) call my(restf)

        mv%ref = 0
        allocate( mv%o )
        mv%o%g = x_ghost()
        mv%o%g_storage = x_ghost()

        call my(mb,mv%o%mb)

        nullify( vmap )

        worm_mb => wormhole(mv%o%mb)

        mv%o%usage = worm_mb%usage

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        if (present(restf)) then

          ! Open the WAVEFUNCTIONS block
          if (i_access(restf)) tios = findnexttag(restf,"WAVEFUNCTIONS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: WAVEFUNCTIONS block was not found")) goto 400
          if (i_access(restf)) call openblock(restf)

          ! Open the PARAMETERS block
          if (i_access(restf)) tios = findfirsttag(restf,"PARAMETERS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: PARAMETERS block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)

          ! Read the k-point
          if (i_access(restf)) tios = findfirsttag(restf,"K-POINT")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: K-POINT tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 3
            call readf(kpt,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,kpt)
          if (error(.not.(kpt .in. nbhd(worm_mb%kpt,tol_nbhd)),"ERROR: k-points do not match")) goto 100

          ! Read the number of bands
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

          ! Read the cutoff
          if (i_access(restf)) tios = findfirsttag(restf,"CUTOFF")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: CUTOFF tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_double
            ndata = 1
            call readf(r_cutoff,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_cutoff)
          if (error(r_cutoff > x_cutoff(worm_mb%lay),"ERROR: restart wavefunctions cutoff > mesh cutoff")) goto 100

          ! Read the number of g-points
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_G-POINTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_G-POINTS tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long
            ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            r_ngt = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,r_ngt)

          ! Close the PARAMETERS block
100       if (i_access(restf)) call closeblock(restf)
          if (error()) goto 300

          ! Form the map of the packed wavefunction data
          call form_vmap(mb,r_cutoff,vmap)

          ! Initialize mv%o%mat
          select case (mv%o%usage)
          case (NORMAL)
            call random_init_i(mv%o) ; if (error()) goto 300
          case (AUXILIARY)
            nullify( mv%o%mat )
          end select

          ! Open the COEFFICIENTS block
          if (i_access(restf)) tios = findfirsttag(restf,"COEFFICIENTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: COEFFICIENTS block was not found")) goto 300
          if (i_access(restf)) call openblock(restf)

          ! Allocate space for the coefficients
          allocate( c1(r_ngt) )

          nb = x_n_bands(mv%o%mb)
          do ib = 1,nb

            ! Read and distribute the data for spin group 1
            if (i_access(restf)) then
              dsize = sizeof_double
              ndata = 2*r_ngt
              call readf(c1,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            end if
            if (i_comm(restf)) call broadcast(FILE_SCOPE,c1)
            select case (msg)
            case (1)
              select case (mv%o%usage)
              case (NORMAL)
                do ig = 1,size(mv%o%mat,1)
                  if (vmap(worm_mb%first_g-1+ig) == 0) cycle
                  mv%o%mat(ig,ib) = c1(vmap(worm_mb%first_g-1+ig))
                end do
              end select
            end select

            ! Read and distribute the data for spin group 2
            select case (nsg)
            case (2)
              if (i_access(restf)) then
                dsize = sizeof_double
                ndata = 2*r_ngt
                call readf(c1,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              end if
              if (i_comm(restf)) call broadcast(FILE_SCOPE,c1)
              select case (msg)
              case (2)
                select case (mv%o%usage)
                case (NORMAL)
                  do ig = 1,size(mv%o%mat,1)
                    if (vmap(worm_mb%first_g-1+ig) == 0) cycle
                    mv%o%mat(ig,ib) = c1(vmap(worm_mb%first_g-1+ig))
                  end do
                end select
              end select
            end select

            ! Exit if the number of bands read is equal to the requested number of bands
            if (ib == r_nb) exit

          end do

          ! Close the COEFFICIENTS block
200       if (i_access(restf)) call closeblock(restf)

          ! Close the WAVEFUNCTIONS block
300       if (i_access(restf)) call closeblock(restf)

        else

          ! Initialize mv%o%mat
          select case (mv%o%usage)
          case (NORMAL)
            select case (trim(init))
            case ("zeros")
              call zeros_init_i(mv%o)
            case ("random")
              call random_init_i(mv%o)
            case ("diagnostic")
              call diagnostic_init_i(mv%o)
            end select
          case (AUXILIARY)
            nullify( mv%o%mat )
          end select

        end if

400     if (associated( vmap )) deallocate( vmap )
        if (allocated( c1 ))   deallocate( c1 )

        nullify( worm_mb )

        call glean(thy(mb))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit multivector_mod::constructor_mv")) continue

      end function

      subroutine update_mv(mv,mb)
!doc$ subroutine update(mv,mb)
        type(multivector_obj) :: mv
        type(multibasis_obj) :: mb
!       modifies: mv
!       effects: Updates mv.
!       errors: Changes in the number of G vectors or the number of bands.

!cod$        
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)
        call my(mb)

        worm_mb => wormhole(mb)
        if (associated(mv%o%mat)) then
           if (error( (size(worm_mb%gpt,1) /= size(mv%o%mat,1)) .or. &
            (x_n_bands(mb) /= size(mv%o%mat,2)), &
            "ERROR: multivector unable to update a dimension change in the basis")) goto 100
        end if
        if (x_ghost(mb) /= x_ghost(mv%o%mb)) then
          mv%o%mb = mb
          mv%o%g = x_ghost()
        end if

100     nullify( worm_mb )

        call glean(thy(mv))
        call glean(thy(mb))

        if (error("Exit multivector_mod::update_mv")) continue

      end subroutine

      subroutine my_mv(mv)
!doc$ subroutine my(mv)
        type(multivector_obj) :: mv

!cod$
        mv%ref = mv%ref + 1
      end subroutine

      subroutine my_new_mv(mvi,mv)
!doc$ subroutine my(mvi,mv)
        type(multivector_obj) :: mvi, mv
!       errors: Not able to allocate space.

!cod$
        mv%ref = 1
        if (mvi%ref < 1) then
          mv%o => mvi%o
        else
          allocate( mv%o )
          mv%o%g = mvi%o%g
          mv%o%g_storage = x_ghost()
          mv%o%usage = mvi%o%usage
          call my(mvi%o%mb,mv%o%mb)
          if (associated(mvi%o%mat)) then
             allocate( mv%o%mat(size(mvi%o%mat,1),size(mvi%o%mat,2)) )
             mv%o%mat = mvi%o%mat
          else
             nullify(mv%o%mat)
          end if

        end if
        if (error("Exit multivector_mod:: my_new_mv")) continue
      end subroutine

      function thy_mv(mv) result(mvo)
!doc$ function thy(mv) result(mvo)
        type(multivector_obj) :: mv, mvo

!cod$
        mv%ref = mv%ref - 1
        mvo%ref = mv%ref
        mvo%o => mv%o
      end function
   
      subroutine glean_mv(mv)
!doc$ subroutine glean(mv)
        type(multivector_obj) :: mv

!cod$
        if (mv%ref < 1) then
          call glean(thy(mv%o%mb))
          if (associated( mv%o%mat )) deallocate( mv%o%mat )
          deallocate( mv%o )
        end if
      end subroutine

      subroutine bequeath_mv(mv)
!doc$ subroutine bequeath(mv)
        type(multivector_obj) :: mv

!cod$
        continue
      end subroutine

      subroutine assign_mv(mv,mv2)
!doc$ subroutine assignment(=)(mv,mv2)
        type(multivector_obj), intent(inout) :: mv
        type(multivector_obj), intent(in) :: mv2
!       requires: mv%o%usge = mv2%o%usage.

!cod$
        type(multivector_obj) :: mvt

        if (error(mv%o%usage /= mv2%o%usage,"Error. mv usage /= to mv2 usage")) goto 100 
        
        call my(mv2)
        if (x_ghost(mv%o%mb) == x_ghost(mv2%o%mb)) then
          mv%o%g = mv2%o%g
          if (associated(mv2%o%mat)) then
             mv%o%mat = mv2%o%mat
          else
             nullify(mv%o%mat)
          end if
        else
          call glean(thy(mv%o%mb))
          if (associated( mv%o%mat )) deallocate( mv%o%mat )
          deallocate( mv%o )
          call my(mv2,mvt)
          mv%o => mvt%o
        end if
100     call glean(thy(mv2))   
      end subroutine

      function wormhole_mv(mv) result(wh)
!doc$ function wormhole(mv) result(wh)
        type(multivector_obj) :: mv
        type(multivector_rep), pointer :: wh
!       effects: Points wh at the wormhole mv underlying.
!       errors: Wormhole is an implementation dependent thing and you should know what you're doing.

!cod$
        call my(mv)
        wh => mv%o
        call glean(thy(mv))
      end function

      function mv_ref(mv) result(r)
!doc$ function x_ref(mv) result(r)
        type(multivector_obj) :: mv
        integer :: r
!       effects: Returns the reference count of mv.

!cod$
        r = mv%ref
        call glean(mv)
      end function

      function mv_ghost(mv) result(g)
!doc$ function x_ghost(mv) result(g)
        type(multivector_obj) :: mv
        type(ghost) :: g
!       effects: Returns the ghost of mv.

!cod$
        call my(mv)
        g = mv%o%g
        call glean(thy(mv))
      end function

      function mv_multibasis(mv) result(mb)
!doc$ function x_multibasis(mv) result(mb)
        type(multivector_obj) :: mv
        type(multibasis_obj) :: mb
!       effects: Returns the mv multibasis.

!cod$
        call my(mv)
        call my(mv%o%mb,mb)
        call bequeath(thy(mb))
        call glean(thy(mv))
      end function

      function mv_n_bands(mv) result(n)
!doc$ function x_n_bands(mv) result(n)
        type(multivector_obj) :: mv
        integer :: n
!       effects: Returns the number of bands in mv%o%multibasis.

!cod$
        call my(mv)
        n = x_n_bands(mv%o%mb)
        call glean(thy(mv))
      end function
      
      function mv_n_gvectors(mv) result(n)
!doc$ function x_n_gvectors(mv) result(n)
        type(multivector_obj) :: mv
        integer :: n
!       effects: Returns the total number of G vectors in mv%o%mb.

!cod$
        call my(mv)
        n = x_n_gvectors(mv%o%mb)
        call glean(thy(mv))
      end function

      subroutine put_mv(mv,ib,g)
!doc$ subroutine put(mv,ib,g)
        type(multivector_obj) :: mv
        integer, intent(in) :: ib
        type(grid_obj) :: g
!       requires: g have SGROUP scope.
!       modifies: g
!       effects: Copies mv%o%mat(:,ib) to g with kind CSF_KIND.
!       errors: ib out of range.

!cod$
        integer :: grp, nb, ng, ngt, npr
        integer, dimension(3) :: nd
        real(double), dimension(:), allocatable :: wts
        complex(double), dimension(:), allocatable :: vec
        complex(double), dimension(:,:), allocatable :: mat
        complex(double), dimension(:,:,:), pointer :: c1, c2
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multibasis_rep), pointer :: worm_mb
        type(grid_obj) :: gr1

        call my(mv)
        call my(g)

        nullify( c1 )
        nullify( c2 )

        worm_mb => wormhole(mv%o%mb)

        if (error( ((ib<1).or.(x_n_bands(mv%o%mb)<ib)),"ERROR: number of weights larger than number of bands")) goto 100
        if (error(x_ghost(x_layout(g)) /= x_ghost(worm_mb%lay),"ERROR: mismatched layouts")) goto 100

        allocate(wts(x_n_bands(mv%o%mb)))
        wts = 0.0_double
        wts(ib) = 1.0_double
        
        nd = x_dims(worm_mb%lay)
        ng = size(mv%o%mat,1)
        nb = size(mv%o%mat,2)
        ngt = size(worm_mb%gridmap,2)
        npr = mpi_nprocs(KGROUP)
        allocate( vec(ngt), mat(ng,npr) )
        
        wf1 => datalink(worm_mb%isplan1)

        call alloc(c1,worm_mb%lay,S_TYPE)
        c1 = (0.0_double,0.0_double)
        
        select case (mv%o%usage)
        case (NORMAL)
          do grp = 1,size(worm_mb%first_band)
            call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,mv%o%mat)
            call band_remap(mv%o%mb,grp,mat,vec) ; if (error()) goto 100
            if (worm_mb%band_participant(grp)) then
              call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
              call fft_serial(Q_TO_R,worm_mb%isplan1)
              call point_mxs(wf1,(wts(worm_mb%my_band(grp))/sqrt(worm_mb%cellvol)))
              c1 = c1 + wf1
            end if
          end do
        end select

        call my(grid(worm_mb%lay,KGROUP),gr1)
        call alloc(c2,worm_mb%lay,S_TYPE)
        call allreduce(KGROUP,MPI_SUM,c1,c2)
        call put(c2,gr1,CSP_KIND)
          
        call saxpby(1.0_double,g,1.0_double,gr1)
          
100     if (allocated( wts )) deallocate( wts )
        if (allocated( vec )) deallocate( vec )
        if (allocated( mat )) deallocate( mat )
        if (associated( c2 )) deallocate( c2 )
        if (associated( c1 )) deallocate( c1 )
        nullify( wf1 )
        nullify( worm_mb )

        call glean(thy(gr1))

        call glean(thy(mv))
        call glean(thy(g))

        if (error("Exit multivector_mod::put_mv")) continue

      end subroutine

      subroutine portion_mv_real_scalar(r,mv)
!doc$ subroutine portion(r,mv)
        real(double), intent(in) :: r
        type(multivector_obj) :: mv
!       modifies: mv
!       effects: mv <--- r*mv

!cod$
        integer :: nb, ng

        call my(mv)
        
        if (associated(mv%o%mat)) then
          ng = size(mv%o%mat,1)
          nb = size(mv%o%mat,2)
          call kernel_portion_mv_rs_i(ng,nb,r,mv%o%mat)
        end if
        mv%o%g = x_ghost()

        call glean(thy(mv))

      end subroutine

      subroutine portion_mv_complex_scalar(c,mv)
!doc$ subroutine portion(c,mv)
        complex(double), intent(in) :: c
        type(multivector_obj) :: mv
!       modifies: mv
!       effects: mv <--- c*mv

!cod$
        integer :: nb, ng

        call my(mv)

        if (associated(mv%o%mat)) then
          ng = size(mv%o%mat,1)
          nb = size(mv%o%mat,2)

          call kernel_portion_mv_cs_i(ng,nb,c,mv%o%mat)
        end if
        mv%o%g = x_ghost()
 
        call glean(thy(mv))

      end subroutine

      subroutine portion_mv_real_vector(rv,mv)
!doc$ subroutine portion(rv,mv)
        real(double), dimension(:), intent(in) :: rv
        type(multivector_obj) :: mv
!       modifies: mv
!       effects: mv <--- rv*mv
!       errors: rv is not of length n_bands.

!cod$
        integer :: nb, ng

        call my(mv)

        if (associated(mv%o%mat)) then
          ng = size(mv%o%mat,1)
          nb = size(mv%o%mat,2)

          if (error((size(rv) /= nb),"ERROR: improper size for rv")) goto 100

          call kernel_portion_mv_rv_i(ng,nb,rv,mv%o%mat)
        end if
        mv%o%g = x_ghost()

100     call glean(thy(mv))

        if (error("Exit multivector_mod::ERROR: portion_mv_real_vector")) continue

      end subroutine

      subroutine portion_mv_complex_vector(cv,mv)
!doc$ subroutine portion(cv,mv)
        complex(double), dimension(:), intent(in) :: cv
        type(multivector_obj) :: mv
!       modifies: mv
!       effects: mv <--- cv*mv
!       errors: cv is not of length n_bands.

!cod$
        integer :: nb, ng

        call my(mv)

        if (associated(mv%o%mat)) then
           ng = size(mv%o%mat,1)
           nb = size(mv%o%mat,2)

           if (error((size(cv) /= nb),"ERROR: improper size for cv")) goto 100

           call kernel_portion_mv_cv_i(ng,nb,cv,mv%o%mat)
        end if
        mv%o%g = x_ghost()


100     call glean(thy(mv))

        if (error("Exit multivector_mod::ERROR: portion_mv_complex_vector")) continue

      end subroutine

      subroutine combine_2mv_real_scalar(r1,mv1,r2,mv2)
!doc$ subroutine combine(r1,mv1,r2,mv2)
        real(double), intent(in) :: r1, r2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
!       modifies: mv1
!       effects:  mv1 <--- r1*mv1 + r2*mv2
!       errors: mv1 and mv2 have different multibases.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        if (associated(mv1%o%mat)) then
           ng = size(mv1%o%mat,1)
           nb = size(mv1%o%mat,2)
           if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100

           if (associated(mv2%o%mat)) then
              call kernel_combine_2mv_rs_i(ng,nb,r1,mv1%o%mat,r2,mv2%o%mat)
           else
              call kernel_portion_mv_rs_i(ng,nb,r1,mv1%o%mat)
           end if
        end if
        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::combine_2mv_real_scalar")) continue

      end subroutine

      subroutine combine_2mv_complex_scalar(c1,mv1,c2,mv2)
!doc$ subroutine combine(c1,mv1,c2,mv2)
        complex(double), intent(in) :: c1, c2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
!       modifies: mv1
!       effects: mv1 <--- c1*mv1 + c2*mv2
!       errors: mv1 and mv2 have different multibases

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        if (associated(mv1%o%mat)) then
           ng = size(mv1%o%mat,1)
           nb = size(mv1%o%mat,2)
           if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100

           if (associated(mv2%o%mat)) then
              call kernel_combine_2mv_cs_i(ng,nb,c1,mv1%o%mat,c2,mv2%o%mat)
           else
              call kernel_portion_mv_cs_i(ng,nb,c1,mv1%o%mat)
           end if
        end if

        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::combine_2mv_complex_scalar")) continue

      end subroutine

      subroutine combine_2mv_real_vector(rv1,mv1,rv2,mv2)
!doc$ subroutine combine(rv1,mv1,rv2,mv2)
        real(double), dimension(:), intent(in) :: rv1, rv2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
!       modifies: mv1
!       effects:  mv1 <--- rv1*mv1 + rv2*mv2
!       errors: mv1 and mv2 have different multibases. rv1 or rv2 is not of length n_bands.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        if (associated(mv1%o%mat)) then
           ng = size(mv1%o%mat,1)
           nb = size(mv1%o%mat,2)

           if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), &
                "ERROR: multibases are different")) goto 100
           if (error((size(rv1) /= nb),"ERROR: improper size for rv1")) goto 100
           if (error((size(rv2) /= nb),"ERROR: improper size for rv2")) goto 100

           if (associated(mv2%o%mat)) then
              call kernel_combine_2mv_rv_i(ng,nb,rv1,mv1%o%mat,rv2,mv2%o%mat)
           else
              call warn("Warning! mvec::combine_2mv_rv. mv2%o%mat not associated")
              call kernel_portion_mv_rv_i(ng,nb,rv1,mv1%o%mat)
           end if
        else
           call warn("Warning! mvec::combine_2mv_rv. mv%o%mat not associated")
        end if

        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::combine_2mv_real_vector")) continue

      end subroutine

      subroutine combine_2mv_complex_vector(cv1,mv1,cv2,mv2)
!doc$ subroutine combine(cv1,mv1,cv2,mv2)
        complex(double), dimension(:), intent(in) :: cv1, cv2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
!       modifies: mv1
!       effects: mv1 <--- cv1*mv1 + cv2*mv2
!       errors: mv1 and mv2 have different multibases. cv1 or cv2 is not of length n_bands.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100
        if (error((size(cv1) /= nb),"ERROR: improper size for cv1")) goto 100
        if (error((size(cv2) /= nb),"ERROR: improper size for cv2")) goto 100

        call kernel_combine_2mv_cv_i(ng,nb,cv1,mv1%o%mat,cv2,mv2%o%mat)
        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::combine_2mv_complex_vector")) continue

      end subroutine

      subroutine combine_3mv_real_vector(rv1,mv1,rv2,mv2,mv3)
!doc$ subroutine combine(rv1,mv1,rv2,mv2,mv3)
        real(double), dimension(:), intent(in) :: rv1, rv2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        type(multivector_obj) :: mv3
!       modifies: mv3
!       effects:  mv3 <--- rv1*mv1 + rv2*mv2
!       errors: mv1, mv2 and mv3 have different multibases. rv1 or rv2 is not of length n_bands.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)
        call my(mv3)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different 1")) goto 100
        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv3%o%mb), "ERROR: multibases are different 2")) goto 100
        if (error((size(rv1) /= nb),"ERROR: improper size for rv1")) goto 100
        if (error((size(rv2) /= nb),"ERROR: improper size for rv2")) goto 100

        call kernel_combine_3mv_rv_i(ng,nb,rv1,mv1%o%mat,rv2,mv2%o%mat,mv3%o%mat)
        mv3%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(mv3))

        if (error("Exit multivector_mod::combine_3mv_real_vector")) continue

      end subroutine

      subroutine combine_3mv_complex_vector(cv1,mv1,cv2,mv2,mv3)
!doc$ subroutine combine(cv1,mv1,cv2,mv2,mv3)
        complex(double), dimension(:), intent(in) :: cv1, cv2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        type(multivector_obj) :: mv3
!       modifies: mv3
!       effects:  mv3 <--- cv1*mv1 + cv2*mv2
!       errors: mv1, mv2 and mv3 have different multibases. cv1 or cv2 is not of length n_bands.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)
        call my(mv3)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different 1")) goto 100
        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv3%o%mb), "ERROR: multibases are different 2")) goto 100
        if (error((size(cv1) /= nb),"ERROR: improper size for cv1")) goto 100
        if (error((size(cv2) /= nb),"ERROR: improper size for cv2")) goto 100

        call kernel_combine_3mv_cv_i(ng,nb,cv1,mv1%o%mat,cv2,mv2%o%mat,mv3%o%mat)
        mv3%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(mv3))

        if (error("Exit multivector_mod::combine_3mv_complex_vector")) continue

      end subroutine

      subroutine transform_mv_complex(mv,cmat)
!doc$ subroutine transform(mv,cmat)
        type(multivector_obj) :: mv
        complex(double), dimension(:,:), intent(in) :: cmat
!       modifies: mv
!       effects: mv <-- mv*cmat
!       requires: cmat is stored in contiguous memory, e.g., be careful passing sections
!       errors: cmat is not of size n_bands X n_bands.

!cod$
        character(1), parameter :: transa = 'n'
        character(1), parameter :: transb = 'n'
        integer :: nb, nw
        complex(double), parameter :: c0 = (0.0_double,0.0_double), c1 = (1.0_double,0.0_double)
        complex(double), dimension(:,:), pointer :: newmat, tmpmat 

        call my(mv)

        nw = size(mv%o%mat,1)
        nb = size(mv%o%mat,2)

        if (error((size(cmat,1) /= nb) .or. (size(cmat,2) /= nb),"ERROR: improper size for cmat")) goto 100

        allocate( newmat(nw,nb) )
        newmat = (0.0_double,0.0_double)
        call start_timer("multivector: zgemm")
        call zgemm(transa,transb,nw,nb,nb,c1,mv%o%mat,nw,cmat,nb,c0,newmat,nw)
        call stop_timer("multivector: zgemm")
        tmpmat => mv%o%mat
        mv%o%mat => newmat
        mv%o%g = x_ghost()
        mv%o%g_storage = x_ghost()

        if (associated( tmpmat )) deallocate( tmpmat )

100     nullify( newmat )

        call glean(thy(mv))

        if (error("Exit multivector_mod::transform_mv_complex")) continue

      end subroutine

      subroutine transform_2mv_complex(c1,mv1,c2,mv2,cmat)
!doc$ subroutine transform(c1,mv1,c2,mv2,cmat)
        complex(double), intent(in) :: c1, c2
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        complex(double), dimension(:,:), intent(in) :: cmat
!       modifies: mv1
!       effects: mv1 <-- c1*mv1 + c2*mv2*cmat
!       requires: cmat is stored in contiguous memory, e.g., be careful passing sections
!       errors: mv1 and mv2 have different multibases. cmat is not of size n_bands X n_bands.

!cod$
        character(1), parameter :: transa = 'n'
        character(1), parameter :: transb = 'n'
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are diferent")) goto 100
        if (error((size(cmat,1) /= nb) .or. (size(cmat,2) /= nb),"ERROR: improper size for cmat")) goto 100

        if (mv1%o%g_storage == mv2%o%g_storage) then
          call warn("WARNING: Using matmul to perform mv operation")
          mv1%o%mat = c1*mv1%o%mat + c2*matmul(mv2%o%mat,cmat)
        else
          call start_timer("multivector: zgemm")
          call zgemm(transa,transb,ng,nb,nb,c2,mv2%o%mat,ng,cmat,nb,c1,mv1%o%mat,ng)
          call stop_timer("multivector: zgemm")
        end if
        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::transform_2mv_complex")) continue

      end subroutine

      subroutine project_2mv(mv1,mv2)
!doc$ subroutine project(mv1,mv2)
        type(multivector_obj) :: mv1, mv2
!       modifies: mv2
!       effects: mv2 <-- (1,0)*mv2 - (1,0)*mv1*(transpose(conjg(mv1%o%mat))*mv2%o%mat).
!       errors: mv1 and mv2 have different multibases.

!cod$
        character(1) :: transa
        character(1), parameter :: transb = 'n'
        integer :: nb, ng
        complex(double), parameter :: c0 = (0.0_double,0.0_double), cp1 = (1.0_double,0.0_double), cm1 = (-1.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: mat_global, mat_local

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100

        allocate( mat_global(nb,nb), mat_local(nb,nb) )
        if (mv1%o%g_storage == mv2%o%g_storage) then
          call warn("WARNING: Using matmul to perform mv operation")
          mat_local = matmul(transpose(conjg(mv1%o%mat)),mv2%o%mat)
          call allreduce(KGROUP,MPI_SUM,mat_local,mat_global)
          mv2%o%mat = cp1*mv2%o%mat + cm1*matmul(mv1%o%mat,mat_global)
        else
          mat_local = (0.0_double,0.0_double)
          transa = 'c'
          call start_timer("multivector: zgemm")
          call zgemm(transa,transb,nb,nb,ng,cp1,mv1%o%mat,ng,mv2%o%mat,ng,c0,mat_local,nb)
          call stop_timer("multivector: zgemm")
          call allreduce(KGROUP,MPI_SUM,mat_local,mat_global)
          transa = 'n'
          call start_timer("multivector: zgemm")
          call zgemm(transa,transb,ng,nb,nb,cm1,mv1%o%mat,ng,mat_global,nb,cp1,mv2%o%mat,ng)
          call stop_timer("multivector: zgemm")
        end if
        deallocate( mat_global, mat_local )
        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::project_2mv")) continue

      end subroutine

      subroutine project_2mv_mat(mv1,mv2,mat)
!doc$ subroutine project(mv1,mv2,mat)
        type(multivector_obj) :: mv1, mv2
        complex(double), dimension(:,:), intent(out) :: mat
!       requires: mat be of size n_bands x n_bands.
!       modifies: mv2 and mat
!       effects: mat <-- transpose(conjg(mv1%o%mat))*mv2%o%mat. mv2 <-- (1,0)*mv2 - (1,0)*mv1*mat.
!       errors: mv1 and mv2 have different multibases.

!cod$
        character(1) :: transa
        character(1), parameter :: transb = 'n'
        integer :: nb, ng
        complex(double), parameter :: c0 = (0.0_double,0.0_double), cp1 = (1.0_double,0.0_double), cm1 = (-1.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: mat_local

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100

        allocate( mat_local(nb,nb) )
        if (mv1%o%g_storage == mv2%o%g_storage) then
          call warn("WARNING: Using matmul to perform mv operation")
          mat_local = matmul(transpose(conjg(mv1%o%mat)),mv2%o%mat)
          call allreduce(KGROUP,MPI_SUM,mat_local,mat)
          mv2%o%mat = cp1*mv2%o%mat + cm1*matmul(mv1%o%mat,mat)
        else
          mat_local = (0.0_double,0.0_double)
          transa = 'c'
          call start_timer("multivector: zgemm")
          call zgemm(transa,transb,nb,nb,ng,cp1,mv1%o%mat,ng,mv2%o%mat,ng,c0,mat_local,nb)
          call stop_timer("multivector: zgemm")
          call allreduce(KGROUP,MPI_SUM,mat_local,mat)
          transa = 'n'
          call start_timer("multivector: zgemm")
          call zgemm(transa,transb,ng,nb,nb,cm1,mv1%o%mat,ng,mat,nb,cp1,mv2%o%mat,ng)
          call stop_timer("multivector: zgemm")
        end if
        deallocate( mat_local )
        mv1%o%g = x_ghost()

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::project_2mv_mat")) continue

      end subroutine

      subroutine orthonormalize_mv(mv)
!doc$ subroutine orthonormalize(mv)
        type(multivector_obj) :: mv
!       modifies: mv
!       effects: applies the transformation t(mv) such that transpose(conjg(t(mv)))*t(mv) = unit matrix

!cod$
        integer :: nb
        complex(double), parameter :: c0 = (0.0_double,0.0_double), c1 = (1.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: s
        type(multivector_obj) :: tmp_mv

        call my(mv)

        call my(mv,tmp_mv)

        nb = size(mv%o%mat,2)
        allocate( s(nb,nb) )

        s = (0.0_double,0.0_double)
        call overlap(mv,s)
        call inverse_cholesky(s)
        call transform(c0,mv,c1,tmp_mv,s)

        if (allocated( s )) deallocate( s )

100     call glean(thy(tmp_mv))

        call glean(thy(mv))

        if (error("Exit multivector_mod::orthonormalize_mv")) continue

      end subroutine

      subroutine residual_mv(v,hv,e,r)
!doc$ subroutine residual(v,hv,e,r)
        type(multivector_obj) :: v, hv
        type(multivector_obj) :: r
        complex(double), dimension(:), intent(in) :: e
!       modifies: r
!       effects:  r <--- e*v - hv
!       errors: v and hv have different multibases. v and r have different multibases. e is not of size n_bands.

!cod$
        integer :: nb, ng

        call my(v)
        call my(hv)
        call my(r)

        ng = size(v%o%mat,1)
        nb = size(v%o%mat,2)

        if (error(x_ghost(v%o%mb) /= x_ghost(hv%o%mb), "ERROR: multibases are diferent 1")) goto 100
        if (error(x_ghost(v%o%mb) /= x_ghost(r%o%mb), "ERROR: multibases are diferent 2")) goto 100
        if (error( (size(e) /= nb),"ERROR: improper size for e")) goto 100

        call kernel_residual_mv_i(ng,nb,v%o%mat,hv%o%mat,e,r%o%mat)
        r%o%g = x_ghost()

100     call glean(thy(v))
        call glean(thy(hv))
        call glean(thy(r))

        if (error("Exit multivector_mod::residual_mv")) continue

      end subroutine

      subroutine decompose_mv(mv,sd,mode,rsa,b1)
!doc$ subroutine decompose(mv,sd,mode,rsa,b1)
        type(multivector_obj) :: mv
        real(double), dimension(:,:), intent(in) :: sd
        character(line_len), intent(in) :: mode
        real(double), dimension(:,:,:), intent(out) :: rsa
        integer :: b1
!       requires: rsa dimensions be 9 x number-of-bands x number-of-sites (the second dimension of sd).
!       modifies: rsa
!       effects: Returns the decomposition of the Kohn-Sham functions into s, p, & d spherical harmonics.

!cod$
        logical :: found
        character(line_len) :: tag

        call my(mv)

        call arglc("dcomp_memory",tag,found)
        if (.not.found) tag = "m"
        select case (trim(tag))
        case ("s")
          call decompose_small_i(mv,sd,mode,rsa,b1)
        case ("m")
          call decompose_medium_i(mv,sd,mode,rsa,b1)
        case ("l")
          call decompose_large_i(mv,sd,mode,rsa,b1)
        case default
          if (error(.true.,"ERROR: dcomp_memory tag was not recognized")) goto 100
        end select

100     call glean(thy(mv))

        if (error("Exit multivector_mod::decompose_mv")) continue

      end subroutine

      subroutine multiply_mv_real(mv,rv)
!doc$ subroutine multiply(mv,rv)
        type(multivector_obj) :: mv
        real(double), dimension(:), intent(out) :: rv
!       modifies: rv
!       effects: rv <-- dot_product(mv,mv)
!       errors: rv is not of size n_bands.

!cod$
        integer :: nb, ng

        call my(mv)

        ng = size(mv%o%mat,1)
        nb = size(mv%o%mat,2)

        if (error( (size(rv) /= nb),"ERROR: target vector is not of size n_bands")) goto 100

        call kernel_multiply_mv_r_i(ng,nb,mv%o%mat,rv)

100     call glean(thy(mv))

        if (error("Exit multivector_mod::multiply_mv_real")) continue

      end subroutine

      subroutine multiply_2mv_real(mv1,mv2,rv)
!doc$ subroutine multiply(mv1,mv2,rv)
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        real(double), dimension(:), intent(out) :: rv
!       modifies: rv
!       effects:  rv <-- real( dot_product(mv1,mv2) )
!       errors: mv1 and mv2 have different multibases. rv is not of size n_bands.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are diferent")) goto 100
        if (error( (size(rv) /= nb),"ERROR: improper size for rv")) goto 100

        call kernel_multiply_2mv_r_i(ng,nb,mv1%o%mat,mv2%o%mat,rv)

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::multiply_2mv_real")) continue

      end subroutine

      subroutine multiply_2mv_complex(mv1,mv2,cv)
!doc$ subroutine multiply(mv1,mv2,cv)
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        complex(double), dimension(:), intent(out) :: cv
!       modifies: cv
!       effects:  cv <-- dot_product(mv1,mv2)
!       errors : mv1 and mv2 have different multibases. cv is not of size n_bands.

!cod$
        integer :: nb, ng

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100
        if (error( (size(cv) /= nb),"ERROR: improper size for cv")) goto 100

        call kernel_multiply_2mv_c_i(ng,nb,mv1%o%mat,mv2%o%mat,cv)

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::multiply_2mv_complex")) continue

      end subroutine

      subroutine overlap_mv(mv,cmat)
!doc$ subroutine overlap(mv,cmat)
        type(multivector_obj) :: mv
        complex(double), dimension(:,:), intent(out) :: cmat
!       modifies: cmat
!       effects: cmat <-- transpose(conjg(mv%o%mat))*mv%o%mat).
!       errors: cmat is not of size n_bands X n_bands.

!cod$
        character(1), parameter :: trans = 'c', uplo = 'u'
        complex(double), parameter :: c0 = (0.0_double,0.0_double), c1 = (1.0_double,0.0_double)
        integer :: i, j, n, nb, ng
        complex(double), dimension(:), allocatable :: r1_local, r1_global

        call my(mv)

        ng = size(mv%o%mat,1)
        nb = size(mv%o%mat,2)

        if (error( (size(cmat,1) /= nb) .or. (size(cmat,2) /= nb), "ERROR: improper size for cmat")) goto 100

        cmat = (0.0_double,0.0_double)
        call zherk(uplo,trans,nb,ng,c1,mv%o%mat,ng,c0,cmat,nb)
        allocate( r1_local(nb*(nb+1)/2), r1_global(nb*(nb+1)/2) )
        n = 0
        do i = 1,nb
          n = n + 1
          r1_local(n) = cmat(i,i)
          do j = 1,(i - 1)
            n = n + 1
            r1_local(n) = cmat(j,i)
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,r1_local,r1_global)
        n = 0
        do i = 1,nb
          n = n + 1
          cmat(i,i) = r1_global(n)
          do j = 1,(i - 1)
            n = n + 1
            cmat(j,i) = r1_global(n)
            cmat(i,j) = conjg(r1_global(n))
          end do
        end do

        if (allocated( r1_local )) deallocate( r1_local )
        if (allocated( r1_global )) deallocate( r1_global )

100     call glean(thy(mv))

        if (error("Exit multivector_mod::overlap_mv")) continue

      end subroutine

      subroutine overlap_2mv(mv1,mv2,cmat)
!doc$ subroutine overlap(mv1,mv2,cmat)
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        complex(double), dimension(:,:), intent(out) :: cmat
!       modifies: cmat
!       effects: cmat <-- transpose(conjg(mv1%o%mat))*mv2%o%mat).
!       errors: mv1 and mv2 have different multibases. cmat is not of size n_bands X n_bands.

!cod$
        character(1), parameter :: trans1 = 'c', trans2 = 'n'
        integer :: nb, ng
        complex(double), parameter :: c0 = (0.0_double,0.0_double), c1 = (1.0_double,0.0_double)
        complex(double), dimension(:,:), allocatable :: cmat_local

        call my(mv1)
        call my(mv2)

        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb), "ERROR: multibases are different")) goto 100
        if (error( (size(cmat,1) /= nb) .or. (size(cmat,2) /= nb), "ERROR: improper size for cmat")) goto 100

        allocate( cmat_local(nb,nb) )
        if (mv1%o%g_storage == mv2%o%g_storage) then
          call warn("WARNING: Using matmul to perform mv operation")
          cmat_local = matmul(transpose(conjg(mv1%o%mat)),mv2%o%mat)
        else
          cmat_local = (0.0_double,0.0_double)
          call start_timer("multivector: zgemm")
          call zgemm(trans1,trans2,nb,nb,ng,c1,mv1%o%mat,ng,mv2%o%mat,ng,c0,cmat_local,nb)
          call stop_timer("multivector: zgemm")
        end if
        call allreduce(KGROUP,MPI_SUM,cmat_local,cmat)

        if (allocated( cmat_local )) deallocate( cmat_local )

100     call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::overlap_2mv")) continue

      end subroutine

      subroutine filter_mv(alpha,mv1,beta,mv2,f)
!doc$ subroutine apply_filter(alpha,mv1,beta,mv2,f)
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        type(grid_obj) :: f
        complex(double), intent(in) :: alpha, beta
!       modifies: mv1
!       effects: mv1 <--- alpha*mv1 + beta*(f.*.mv2) where .*. applies fourier filter in f to a multivector.
!       errors: mv1 and mv2 have different multibases. f incompatible with mv multibases.

!cod$
        integer :: nb, ng
        integer, dimension(3) :: nd
        complex(double), dimension(:,:,:), pointer :: filter
        type(grid_obj) :: myf
        type(multibasis_rep), pointer :: worm_mb

        call my(mv1)
        call my(mv2)
        call my(f)

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb),"ERROR: multibases are not the same")) goto 100
        if (error(x_ghost(x_layout(f)) /= x_ghost(x_layout(mv1%o%mb)),"ERROR: mismatched layouts")) goto 100

        worm_mb => wormhole(mv1%o%mb)

        nd = x_dims(worm_mb%lay)
        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)

        call my(f,myf)
        call take(filter,myf,CSF_KIND)
        call glean(thy(myf))

        call kernel_filter_mv_i(worm_mb%first_g,worm_mb%last_g,ng,nb,nd,alpha,beta,mv1%o%mat,mv2%o%mat,worm_mb%gridmap,filter)
        mv1%o%g = x_ghost()

        if (associated( filter )) deallocate( filter )

100     nullify( worm_mb )

        call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(f))

        if (error("Exit multivector_mod::filter_mv")) continue

      end subroutine

      subroutine scale_mv(alpha,mv1,beta,mv2,f)
!doc$ subroutine apply_field(alpha,mv1,beta,mv2,f)
        type(multivector_obj) :: mv1
        type(multivector_obj) :: mv2
        type(grid_obj) :: f
        complex(double), intent(in) :: alpha, beta
!       modifies: mv1
!       effects: mv1 <--- alpha*mv1 + beta*(f.*.mv2) where .*. applies a field to a mutlivector.
!       errors: mv1 and mv2 have different multibases. f and mv1 have different layouts.

!cod$
        integer :: grp, nb, ng, ngt, npr
        integer, dimension(3) :: nd
        real(double), dimension(:,:,:), pointer :: field
        complex(double), dimension(:), allocatable :: vec
        complex(double), dimension(:,:), allocatable :: mat
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multibasis_rep), pointer :: worm_mb
        type(grid_obj) :: myf

        call my(mv1)
        call my(mv2)
        call my(f)

        nullify( field )

        if (error(x_ghost(mv1%o%mb) /= x_ghost(mv2%o%mb),"ERROR: mutibases are different")) goto 200
        if (error(x_ghost(x_layout(f)) /= x_ghost(x_layout(mv1%o%mb)),"ERROR: mismatched layouts")) goto 200

        worm_mb => wormhole(mv1%o%mb)

        wf1 => datalink(worm_mb%isplan1)

        call my(f,myf)
        call take(field,myf,RS_KIND)

        nd = x_dims(worm_mb%lay)
        ng = size(mv1%o%mat,1)
        nb = size(mv1%o%mat,2)
        ngt = size(worm_mb%gridmap,2)
        npr = mpi_nprocs(KGROUP)
        allocate( vec(ngt), mat(ng,npr) )

        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,mv2%o%mat)
          call band_remap(mv1%o%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call kernel_mm_cr_i(nd,wf1,field)
            call fft_serial(R_TO_Q,worm_mb%isplan1)
            call gather_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
          end if
          call band_remap(mv1%o%mb,grp,vec,mat) ; if (error()) goto 100
          call augment_vectors_ab_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,mv1%o%mat,alpha,beta)
        end do
        mv1%o%g = x_ghost()

100     call glean(thy(myf))

        if (associated( field )) deallocate( field )
        if (allocated( vec )) deallocate( vec )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )

200     nullify( worm_mb )

        call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(f))

        if (error("Exit multivector_mod::apply_field")) continue

      end subroutine

      subroutine add_grid_density_mv(mv,wts,den)
!doc$ subroutine add_grid_density(mv,wts,den)
        type(multivector_obj) :: mv
        real(double), dimension(:) :: wts
        type(grid_obj) :: den
!       modifies: den
!       effects: Adds the densities from vectors in mv to den with weights wts.
!       errors: size(wts) > x_n_bands(mv). den layout /= mv layout.

!cod$
        integer :: grp, nb, ng, ngt, npr
        integer, dimension(3) :: nd
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:), allocatable :: vec
        complex(double), dimension(:,:), allocatable :: mat
        complex(double), dimension(:,:,:), pointer :: wf1
        type(multibasis_rep), pointer :: worm_mb
        type(grid_obj) :: gr1

        call my(mv)
        call my(den)

        nullify( r1, r2 )

        worm_mb => wormhole(mv%o%mb)

        if (error(size(wts) > x_n_bands(mv%o%mb),"ERROR: number of weights larger than number of bands")) goto 100
        if (error(x_ghost(x_layout(den)) /= x_ghost(worm_mb%lay),"ERROR: mismatched layouts")) goto 100

        nd = x_dims(worm_mb%lay)
        ng = size(mv%o%mat,1)
        nb = size(mv%o%mat,2)
        ngt = size(worm_mb%gridmap,2)
        npr = mpi_nprocs(KGROUP)
        allocate( vec(ngt), mat(ng,npr) )

        wf1 => datalink(worm_mb%isplan1)

        call alloc(r1,worm_mb%lay,S_TYPE)
        r1 = 0.0_double

        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat,nb,mv%o%mat)
          call band_remap(mv%o%mb,grp,mat,vec) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call point_accumden(r1,wf1,(wts(worm_mb%my_band(grp))/worm_mb%cellvol))
          end if
        end do

        call my(grid(worm_mb%lay,KGROUP),gr1)
        call alloc(r2,worm_mb%lay,S_TYPE)
        call allreduce(KGROUP,MPI_SUM,r1,r2)
        call put(r2,gr1,RS_KIND)

        call saxpby(1.0_double,den,1.0_double,gr1)

100     if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (allocated( vec )) deallocate( vec )
        if (allocated( mat )) deallocate( mat )
        nullify( wf1 )
        nullify( worm_mb )

        call glean(thy(gr1))

        call glean(thy(mv))
        call glean(thy(den))

        if (error("Exit multivector_mod::add_grid_density_mv")) continue

      end subroutine

      subroutine add_grid_density_2mv(mvec1,mvec2,den)
!doc$ subroutine add_grid_density(mvec1,mvec2,den)
        type(multivector_obj) :: mvec1
        type(multivector_obj) :: mvec2
        type(grid_obj) :: den
!       modifies: den
!       effects: Adds (conjugate(mvec1(r))*mvec2(r) + conjugate(mvec2(r))*mvec1(r)) to den(r). If den is empty it is zeroed.
!       errors: mvec1 and mvec2 have different multibases.  den layout /= mvec1 layout.

!cod$
        integer :: grp, nb, ng, ngt, npr
        integer, dimension(3) :: nd
        real(double), dimension(:,:,:), pointer :: r1, r2
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: wf1, wf2
        type(multibasis_rep), pointer :: worm_mb
        type(grid_obj) :: gr1

        call my(mvec1)
        call my(mvec2)
        call my(den)

        nullify( r1, r2 )

        worm_mb => wormhole(mvec1%o%mb)

        if (error(x_ghost(mvec1%o%mb) /= x_ghost(mvec2%o%mb),"ERROR: mutibases are different")) goto 100
        if (error(x_ghost(x_layout(den)) /= x_ghost(worm_mb%lay),"ERROR: mismatched layouts")) goto 100

        nd = x_dims(worm_mb%lay)
        ng = size(mvec1%o%mat,1)
        nb = size(mvec1%o%mat,2)
        ngt = size(worm_mb%gridmap,2)
        npr = mpi_nprocs(KGROUP)
        allocate( vec1(ngt), vec2(ngt), mat1(ng,npr), mat2(ng,npr) )

        wf1 => datalink(worm_mb%isplan1)
        wf2 => datalink(worm_mb%isplan2)

        call alloc(r1,worm_mb%lay,S_TYPE)
        r1 = 0.0_double

        do grp = 1,size(worm_mb%first_band)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat1,nb,mvec1%o%mat)
          call extract_vectors_i(worm_mb%first_band(grp),worm_mb%last_band(grp),ng,npr,mat2,nb,mvec2%o%mat)
          call band_remap(mvec1%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call band_remap(mvec2%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb%band_participant(grp)) then
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt,worm_mb%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call fft_serial(Q_TO_R,worm_mb%isplan2)
            r1 = r1 + dble(conjg(wf1)*wf2 + wf1*conjg(wf2))/worm_mb%cellvol
          end if
        end do

        call alloc(r2,worm_mb%lay,S_TYPE)
        call allreduce(KGROUP,MPI_SUM,r1,r2)
        call my(grid(worm_mb%lay,KGROUP),gr1)
        call put(r2,gr1,RS_KIND)

        call saxpby(1.0_double,den,1.0_double,gr1)

100     if (associated( r1 )) deallocate( r1 )
        if (associated( r2 )) deallocate( r2 )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        nullify( wf1, wf2 )
        nullify( worm_mb )

        call glean(thy(gr1))

        call glean(thy(mvec1))
        call glean(thy(mvec2))
        call glean(thy(den))

        if (error("Exit multivector_mod::add_grid_density_2mv")) continue

      end subroutine

      subroutine distribute_mv(mv,rank)
!doc$ subroutine distribute(mv,rank)
        type(multivector_obj) :: mv
        integer, intent(in) :: rank
!       requires: Must be followed by a call to release_mv (if mv%o%mat is nullified) with no intervening modification of mv.
!       modifies: mv iff mv%o%mat is nullified.
!       effects: Distributes mv%o%mat data to all kgroups.

!cod$
        integer :: nb, ng
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)

        worm_mb => wormhole(mv%o%mb)

        nb = x_n_bands(mv%o%mb)
        ng = size(worm_mb%gpt,1)
        if (.not.associated( mv%o%mat )) allocate( mv%o%mat(ng,nb) )

        call xcomm_broadcast(XKGROUP,rank,mv%o%mat) ; if (error()) goto 100

        nullify( worm_mb )

100     call glean(thy(mv))

        if (error("Exit multivector_mod::distribute_mv")) continue

      end subroutine

      subroutine release_mv(mv)
!doc$ subroutine release(mv)
        type(multivector_obj) :: mv
!       requires: mv%o%mat be associated. Must be preceeded by a call to distribute_mv with no intervening modification of mv.
!       modifies: mv
!       effects: Deallocates mv%o%mat.

!cod$
        call my(mv)

        deallocate( mv%o%mat )

        call glean(thy(mv))

        if (error("Exit multivector_mod::release_mv")) continue

      end subroutine

      subroutine exx_energy_mv(mv,wt,xct,rc,e)
!doc$ subroutine exx_energy(mv,wt,xct,rc,e)
        type(multivector_obj) :: mv
        real(double), dimension(:), intent(in) :: wt
        type(xc_type_obj)             :: xct
        real(double), intent(in)      :: rc
        real(double), intent(inout)   :: e
!       modifies: e
!       effects: Computes non-singular contributions to e from a single BZ sampling point.
!       errors: Passes errors

!cod$
        integer :: ib1, ib2, grp, nb, ngp, ngt, np, npairs
        integer, dimension(3) :: nd
        integer, dimension(:), allocatable :: index1, index2
        real(double) :: e1, e2, pfe, r0, spin_degeneracy, dk(3)
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: c1, wf1, wf2
        type(multibasis_rep), pointer :: worm_mb
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g

        call my(mv)
        call my(xct)

        nullify( ck )

        worm_mb => wormhole(mv%o%mb)

        call my(worm_mb%lay,lay)
        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = 0.0_double
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv%o%mat,2)

        ngp = size(mv%o%mat,1)
        ngt = size(worm_mb%gridmap,2)
        allocate( vec1(ngt), vec2(ngt) )
        allocate( mat1(ngp,np), mat2(ngp,np) )

        c1 => datalink(x_serial_plan(lay))

        wf1 => datalink(worm_mb%isplan1)
        wf2 => datalink(worm_mb%isplan2)

        pfe = (four_pi/worm_mb%cellvol)/spin_degeneracy

        npairs = size(worm_mb%spair_index,2)
        allocate( index1(npairs), index2(npairs) )
        index1 = worm_mb%spair_index(1,:)
        index2 = worm_mb%spair_index(2,:)

        ! compute the energy contributions
        e1 = 0.0_double
        do grp = 1,size(worm_mb%first_spair)
          call extract_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index1,ngp,np,mat1,nb,mv%o%mat)
          call extract_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index2,ngp,np,mat2,nb,mv%o%mat)
          call spair_remap(mv%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call spair_remap(mv%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb%spair_participant(grp)) then
            ib1 = index1(worm_mb%my_spair(grp))
            ib2 = index2(worm_mb%my_spair(grp))
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt,worm_mb%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call fft_serial(Q_TO_R,worm_mb%isplan2)
            call kernel_mmm_ccc_i(nd,c1,wf1,wf2)
            call fft_serial(R_TO_Q,x_serial_plan(lay))
            call kernel_mm_cr_i(nd,c1,ck)
            call fft_serial(Q_TO_R,x_serial_plan(lay))
            call kernel_rsum_mmm_ccc_i(nd,c1,wf1,wf2,r0)
            if (ib1 == ib2) then
              e1 = e1 - pfe*wt(ib1)*wt(ib1)*r0
            else
              e1 = e1 - 2.0_double*pfe*wt(ib1)*wt(ib2)*r0
            end if
          end if
        end do
        call allreduce(KGROUP,MPI_SUM,e1,e2)
        e = e + e2

100     if (allocated( index1 )) deallocate( index1 )
        if (allocated( index2 )) deallocate( index2 )
        if (associated( ck )) deallocate( ck )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        nullify( c1 )
        nullify( wf1 )
        nullify( wf2 )
        nullify( worm_mb )

        call glean(thy(lay))
        call glean(thy(ck_g))

        call glean(thy(xct))
        call glean(thy(mv))

        if (error("Exit multivector_mod::exx_energy_mv")) continue

      end subroutine

      subroutine exx_energy_2mv(mv1,mv2,wt1,wt2,xct,rc,e)
!doc$ subroutine exx_energy(mv1,mv2,wt1,wt2,xct,rc,e)
        type(multivector_obj) :: mv1, mv2
        real(double), dimension(:), intent(in) :: wt1, wt2
        type(xc_type_obj) :: xct
        real(double), intent(in) :: rc
        real(double), intent(inout) :: e
!       requires: mv1 and mv2 have different kpts.
!       modifies: e
!       effects: Computes contributions to e from distinct BZ sampling points.
!       errors: Passes errors

!cod$
        integer :: ib1, ib2, grp, nb, ngp1, ngp2, ngt1, ngt2, np, npairs
        integer, dimension(3) :: nd
        integer, dimension(:), allocatable :: index1, index2
        real(double) :: dk(3), e1, e2, pfe, r0, spin_degeneracy
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: c1, wf1, wf2
        type(multibasis_rep), pointer :: worm_mb1, worm_mb2
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g

        call my(mv1)
        call my(mv2)
        call my(xct)

        nullify( ck )

        worm_mb1 => wormhole(mv1%o%mb)
        worm_mb2 => wormhole(mv2%o%mb)

        call my(worm_mb1%lay,lay)
        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = worm_mb1%kpt - worm_mb2%kpt
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv1%o%mat,2)

        ngp1 = size(mv1%o%mat,1)
        ngt1 = size(worm_mb1%gridmap,2)
        allocate( vec1(ngt1), mat1(ngp1,np) )

        ngp2 = size(mv2%o%mat,1)
        ngt2 = size(worm_mb2%gridmap,2)
        allocate( vec2(ngt2), mat2(ngp2,np) )

        npairs = size(worm_mb1%lpair_index,2)
        allocate( index1(npairs), index2(npairs) )
        index1 = worm_mb1%lpair_index(1,:)
        index2 = worm_mb1%lpair_index(2,:)

        c1 => datalink(x_serial_plan(lay))

        wf1 => datalink(worm_mb1%isplan1)
        wf2 => datalink(worm_mb2%isplan1)

        pfe = (four_pi/worm_mb1%cellvol)/spin_degeneracy

        ! compute the energy contributions
        e1 = 0.0_double
        do grp = 1,size(worm_mb1%first_lpair)
          call extract_vectors_indexed_i(worm_mb1%first_lpair(grp),worm_mb1%last_lpair(grp),npairs,index1,ngp1,np,mat1,nb,mv1%o%mat)
          call extract_vectors_indexed_i(worm_mb2%first_lpair(grp),worm_mb2%last_lpair(grp),npairs,index2,ngp2,np,mat2,nb,mv2%o%mat)
          call lpair_remap(mv1%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call lpair_remap(mv2%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb1%lpair_participant(grp)) then
            ib1 = index1(worm_mb1%my_lpair(grp))
            ib2 = index2(worm_mb2%my_lpair(grp))
            call scatter_wf_i(nd,wf1,ngt1,worm_mb1%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt2,worm_mb2%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb1%isplan1)
            call fft_serial(Q_TO_R,worm_mb2%isplan1)
            call kernel_mmm_ccc_i(nd,c1,wf1,wf2)
            call fft_serial(R_TO_Q,x_serial_plan(lay))
            call kernel_mm_cr_i(nd,c1,ck)
            call fft_serial(Q_TO_R,x_serial_plan(lay))
            call kernel_rsum_mmm_ccc_i(nd,c1,wf1,wf2,r0)
            e1 = e1 - pfe*wt1(ib1)*wt2(ib2)*r0
          end if
        end do
        call allreduce(KGROUP,MPI_SUM,e1,e2)
        e = e + e2

100     if (allocated( index1 )) deallocate( index1 )
        if (allocated( index2 )) deallocate( index2 )
        if (associated( ck )) deallocate( ck )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        nullify( c1 )
        nullify( wf1 )
        nullify( wf2 )
        nullify( worm_mb1 )
        nullify( worm_mb2 )

        call glean(thy(lay))
        call glean(thy(ck_g))

        call glean(thy(xct))
        call glean(thy(mv1))
        call glean(thy(mv2))

        if (error("Exit multivector_mod::exx_energy_2mv")) continue

      end subroutine

      subroutine exx_derivative_mv(mv,wt,xct,rc,mvo)
!doc$ subroutine exx_derivative(mv,wt,xct,rc,mvo)
        type(multivector_obj) :: mv
        real(double), dimension(:), intent(in) :: wt
        type(xc_type_obj)                      ::xct
        real(double), intent(in) :: rc
        type(multivector_obj) :: mvo
!       requires: mv and mvo have the same kpt.
!       modifies: mvo
!       effects: Computes non-singular contributions to mvo from a single BZ sampling point.
!       errors: Passes errors

!cod$
        integer :: ib1, ib2, grp, nb, ngp, ngt, np, npairs
        integer, dimension(3) :: nd
        integer, dimension(:), allocatable :: index1, index2
        real(double) :: pfm, r0, spin_degeneracy, dk(3)
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: c1, wf1, wf2
        type(multibasis_rep), pointer :: worm_mb
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g

        call my(mv)
        call my(mvo)
        call my(xct)

        nullify( ck )

        worm_mb => wormhole(mv%o%mb)

        call my(worm_mb%lay,lay)
        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = 0.0_double
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv%o%mat,2)

        ngp = size(mv%o%mat,1)
        ngt = size(worm_mb%gridmap,2)
        allocate( vec1(ngt), vec2(ngt) )
        allocate( mat1(ngp,np), mat2(ngp,np) )

        c1 => datalink(x_serial_plan(lay))

        wf1 => datalink(worm_mb%isplan1)
        wf2 => datalink(worm_mb%isplan2)

        pfm = (eight_pi/worm_mb%cellvol)/spin_degeneracy

        npairs = size(worm_mb%spair_index,2)
        allocate( index1(npairs), index2(npairs) )
        index1 = worm_mb%spair_index(1,:)
        index2 = worm_mb%spair_index(2,:)

        ! compute the derivative contributions
        do grp = 1,size(worm_mb%first_spair)
          call extract_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index1,ngp,np,mat1,nb,mv%o%mat)
          call extract_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index2,ngp,np,mat2,nb,mv%o%mat)
          call spair_remap(mv%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call spair_remap(mv%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb%spair_participant(grp)) then
            ib1 = index1(worm_mb%my_spair(grp))
            ib2 = index2(worm_mb%my_spair(grp))
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt,worm_mb%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call fft_serial(Q_TO_R,worm_mb%isplan2)
            call kernel_mmm_ccc_i(nd,c1,wf1,wf2)
            call fft_serial(R_TO_Q,x_serial_plan(lay))
            call kernel_mm_cr_i(nd,c1,ck)
            call fft_serial(Q_TO_R,x_serial_plan(lay))
            call kernel_mm_cc_i(nd,wf2,c1)
            call fft_serial(R_TO_Q,worm_mb%isplan2)
            call gather_wf_i(nd,wf2,ngt,worm_mb%gridmap,vec1)
            r0 = -pfm*wt(ib2)
            call kernel_portion_v_rs_i(ngt,r0,vec1)
            if (ib1 == ib2) then
              vec2 = (0.0_double,0.0_double)
            else
              call kernel_mm_ccs_i(nd,wf1,c1)
              call fft_serial(R_TO_Q,worm_mb%isplan1)
              call gather_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec2)
              r0 = -pfm*wt(ib1)
              call kernel_portion_v_rs_i(ngt,r0,vec2)
            end if
          end if
          call spair_remap(mv%o%mb,grp,vec1,mat1) ; if (error()) goto 100
          call spair_remap(mv%o%mb,grp,vec2,mat2) ; if (error()) goto 100
          call augment_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index1,ngp,np,mat1,nb,mvo%o%mat)
          call augment_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index2,ngp,np,mat2,nb,mvo%o%mat)
        end do
        mvo%o%g = x_ghost()

100     if (allocated( index1 )) deallocate( index1 )
        if (allocated( index2 )) deallocate( index2 )
        if (associated( ck )) deallocate( ck )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        nullify( c1 )
        nullify( wf1 )
        nullify( wf2 )
        nullify( worm_mb )

        call glean(thy(lay))
        call glean(thy(ck_g))

        call glean(thy(mv))
        call glean(thy(mvo))
        call glean(thy(xct))

        if (error("Exit multivector_mod::exx_derivative_mv")) continue

      end subroutine

      subroutine exx_derivative_2mv(mv1,mv2,wt2,xct,rc,mvo)
!doc$ subroutine exx_derivative(mv1,mv2,wt2,xct,rc,mvo)
        type(multivector_obj)                  :: mv1, mv2
        real(double), dimension(:), intent(in) :: wt2
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        type(multivector_obj)                  :: mvo
!       requires: mv1 and mv2 have different kpts. mv1 and mvo have the same kpt.
!       modifies: mvo
!       effects: Computes contributions to mvo from distinct BZ sampling points.
!       errors: Passes errors

!cod$
        integer :: ib1, ib2, grp, nb, ngp1, ngp2, ngt1, ngt2, np, npairs
        integer, dimension(3) :: nd
        integer, dimension(:), allocatable :: index1, index2
        real(double) :: dk(3), pfm, r0, spin_degeneracy
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: c1, wf1, wf2
        type(multibasis_rep), pointer :: worm_mb1, worm_mb2
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g

        call my(mv1)
        call my(mv2)
        call my(mvo)
        call my(xct)

        nullify( ck )

        worm_mb1 => wormhole(mv1%o%mb)
        worm_mb2 => wormhole(mv2%o%mb)

        call my(worm_mb1%lay,lay)
        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = worm_mb1%kpt - worm_mb2%kpt
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv1%o%mat,2)

        ngp1 = size(mv1%o%mat,1)
        ngt1 = size(worm_mb1%gridmap,2)
        allocate( vec1(ngt1), mat1(ngp1,np) )

        ngp2 = size(mv2%o%mat,1)
        ngt2 = size(worm_mb2%gridmap,2)
        allocate( vec2(ngt2), mat2(ngp2,np) )

        npairs = size(worm_mb1%lpair_index,2)
        allocate( index1(npairs), index2(npairs) )
        index1 = worm_mb1%lpair_index(1,:)
        index2 = worm_mb1%lpair_index(2,:)

        c1 => datalink(x_serial_plan(lay))

        wf1 => datalink(worm_mb1%isplan1)
        wf2 => datalink(worm_mb2%isplan1)

        pfm = (eight_pi/worm_mb1%cellvol)/spin_degeneracy

        ! compute the derivative contributions
        do grp = 1,size(worm_mb1%first_lpair)
          call extract_vectors_indexed_i(worm_mb1%first_lpair(grp),worm_mb1%last_lpair(grp),npairs,index1,ngp1,np,mat1,nb,mv1%o%mat)
          call extract_vectors_indexed_i(worm_mb2%first_lpair(grp),worm_mb2%last_lpair(grp),npairs,index2,ngp2,np,mat2,nb,mv2%o%mat)
          call lpair_remap(mv1%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call lpair_remap(mv2%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb1%lpair_participant(grp)) then
            ib1 = index1(worm_mb1%my_lpair(grp))
            ib2 = index2(worm_mb2%my_lpair(grp))
            call scatter_wf_i(nd,wf1,ngt1,worm_mb1%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt2,worm_mb2%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb1%isplan1)
            call fft_serial(Q_TO_R,worm_mb2%isplan1)
            call kernel_mmm_ccc_i(nd,c1,wf1,wf2)
            call fft_serial(R_TO_Q,x_serial_plan(lay))
            call kernel_mm_cr_i(nd,c1,ck)
            call fft_serial(Q_TO_R,x_serial_plan(lay))
            call kernel_mm_cc_i(nd,c1,wf2)
            wf1 = c1
            call fft_serial(R_TO_Q,worm_mb1%isplan1)
            call gather_wf_i(nd,wf1,ngt1,worm_mb1%gridmap,vec1)
            r0 = -pfm*wt2(ib2)
            call kernel_portion_v_rs_i(ngt1,r0,vec1)
          end if
          call lpair_remap(mv1%o%mb,grp,vec1,mat1) ; if (error()) goto 100
          call augment_vectors_indexed_i(worm_mb1%first_lpair(grp),worm_mb1%last_lpair(grp),npairs,index1,ngp1,np,mat1,nb,mvo%o%mat)
        end do
        mvo%o%g = x_ghost()

100     if (allocated( index1 )) deallocate( index1 )
        if (allocated( index2 )) deallocate( index2 )
        if (associated( ck )) deallocate( ck )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        nullify( c1 )
        nullify( wf1 )
        nullify( wf2 )
        nullify( worm_mb1 )
        nullify( worm_mb2 )

        call glean(thy(lay))
        call glean(thy(ck_g))

        call glean(thy(xct))
        call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(mvo))

        if (error("Exit multivector_mod::exx_derivative_2mv")) continue

      end subroutine

      subroutine exx_energy_and_derivative_mv(mv,wt,xct,rc,e,mvo)
!doc$ subroutine    exx_energy_and_derivative(mv,wt,xct,rc,e,mvo)
        type(multivector_obj)                  :: mv
        real(double), dimension(:), intent(in) :: wt
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        real(double), intent(inout)            :: e
        type(multivector_obj)                  :: mvo
!       requires: mv and mvo have the same kpt.
!       modifies: e and mvo
!       effects: Computes contributions to e and mvo from the same k-point as mv.
!       errors: Passes errors

!cod$
        integer :: nb, ng
        real(double) :: e1
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)
        call my(mvo)
        call my(xct)

        worm_mb => wormhole(mv%o%mb)

        ! apply the exchange operator to mv
        call start_timer("multivector: exchange_operator_mv")
        select case (worm_mb%exx_comm_method)
        case (COLLECTIVE)
          call exchange_operator_spair_mv_i(mv,wt,xct,rc,mvo)   ; if (error()) goto 100
        case (POINT_TO_POINT)
          call exchange_operator_p2p_mv_i(mv,wt,xct,rc,mvo)     ; if (error()) goto 100
        case (POINT_TO_POINT_2)
          call exchange_operator_p2p_2_mv_i(mv,wt,xct,rc,mvo)   ; if (error()) goto 100
        end select
        call stop_timer("multivector: exchange_operator_mv")

        ! accumulate the energy
        ng = size(mv%o%mat,1)
        nb = size(mv%o%mat,2)
        call kernel_trace_2mv_wt_i(ng,nb,mv%o%mat,mvo%o%mat,wt,e1)
        e = e + e1

100     nullify( worm_mb )

        call glean(thy(mv))
        call glean(thy(mvo))
        call glean(thy(xct))

        if (error("Exit multivector_mod::exx_energy_and_derivative_mv")) continue

      end subroutine

      subroutine exx_energy_and_derivative_2mv(mv1,mv2,wt1,wt2,xct,rc,e,mvo)
!doc$ subroutine     exx_energy_and_derivative(mv1,mv2,wt1,wt2,xct,rc,e,mvo)
        type(multivector_obj)                  :: mv1, mv2
        real(double), dimension(:), intent(in) :: wt1, wt2
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        real(double), intent(inout)            :: e
        type(multivector_obj)                  :: mvo
!       requires: mv1 and mv2 have different kpts. mv1 and mvo have the same kpt.
!       modifies: e and mvo
!       effects: Computes contributions to e and mvo from a different k-point than that of mv1.
!       errors: Passes errors

!cod$
        type(multibasis_rep), pointer :: worm_mb

        call my(mv1)
        call my(mv2)
        call my(mvo)
        call my(xct)

        worm_mb => wormhole(mv1%o%mb)

        ! apply the exchange operator to mv
        call start_timer("multivector: exchange_operator_2mv")
        select case (worm_mb%exx_comm_method)
        case (COLLECTIVE)
          call exchange_operator_lpair_2mv_i(mv1,mv2,wt1,wt2,xct,rc,e,mvo)   ; if (error()) goto 100
        case (POINT_TO_POINT_2)
          call exchange_operator_p2p_2_2mv_i(mv1,mv2,wt1,wt2,xct,rc,e,mvo)   ; if (error()) goto 100
        end select
        call stop_timer("multivector: exchange_operator_2mv")

100     nullify( worm_mb )

        call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(mvo))
        call glean(thy(xct))

        if (error("Exit multivector_mod::exx_energy_and_derivative_2mv")) continue

      end subroutine

      subroutine diary_mv(mv)
!doc$ subroutine diary(mv)
        type(multivector_obj) :: mv
!       effects: Writes mv information to the diary.

!cod$
        call my(mv)

        call diary(mv%o%mb)

        call glean(thy(mv))

      end subroutine

      subroutine write_restart_mv(mv,nrestf)
!doc$ subroutine write_restart(mv,nrestf)
        type(multivector_obj) :: mv
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes restart information to nrestf.

!cod$
        integer :: ib, nb, msg, ngt, nsg
        integer(long) :: dsize, iosl, ndata, s4
        complex(double), dimension(:), allocatable :: c1, c2
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)
        call my(nrestf)

        nsg = mpi_nsgroups()
        msg = mpi_mysgroup()

        worm_mb => wormhole(mv%o%mb)

        ! Start the WAVEFUNCTIONS block
        if (i_access(nrestf)) call startblock(nrestf,"WAVEFUNCTIONS")

        ! Start the PARAMETERS block
        if (i_access(nrestf)) call startblock(nrestf,"PARAMETERS")

        ! Write the k-point
        if (i_access(nrestf)) then
          call writetag(nrestf,"K-POINT")
          dsize = sizeof_double ; ndata = 3
          call writef(worm_mb%kpt,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Write the number of bands (needed to restart with a different number of bands)
        nb = x_n_bands(mv%o%mb)
        if (i_access(nrestf)) then
          call writetag(nrestf,"NUMBER_OF_BANDS")
          dsize = sizeof_long
          ndata = 1
          s4 = nb
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Write the cutoff (needed to restart with a different wf_cutoff))
        if (i_access(nrestf)) then
          call writetag(nrestf,"CUTOFF")
          dsize = sizeof_double
          ndata = 1
          call writef(worm_mb%cutoff,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! Write the number of g-points (needed to restart with a different wf_cutoff)
        ngt = x_n_gvectors(mv%o%mb)
        if (i_access(nrestf)) then
          call writetag(nrestf,"NUMBER_OF_G-POINTS")
          dsize = sizeof_long
          ndata = 1
          s4 = ngt
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)
        end if

        ! End the PARAMETERS block
        if (i_access(nrestf)) call endblock(nrestf)

        ! Start the COEFFICIENTS block
        if (i_access(nrestf)) call startblock(nrestf,"COEFFICIENTS")

        ! Allocate space for the coefficients
        allocate( c1(ngt) )
        allocate( c2(ngt) )

        do ib = 1,nb

          c1 = (0.0_double,0.0_double)

          ! Gather the coefficients to the KGROUP rank 0 processes
          select case (mv%o%usage)
          case (NORMAL)
            c1(worm_mb%first_g:worm_mb%last_g) = mv%o%mat(:,ib)
          end select
          call reduce(KGROUP,MPI_SUM,c1,c2)

          ! Reduce the coefficients to the SGROUP rank 0 processes (c1 contains one set because of the mv%o%usage construct)
          call xcomm_reduce(XKGROUP,MPI_SUM,c2,c1) ; if (error()) exit

          ! Write the spin group 1 coefficients (residing on the FILE_SCOPE rank 0 process)
          if (i_access(nrestf)) then
            dsize = sizeof_double
            ndata = 2*size(c1)
            call writef(c1,dsize,ndata,x_tagfd(nrestf),iosl)
          end if

          ! Reduce the spin group 2 coefficents to the FILE_SCOPE rank 0 process and write
          select case (nsg)
          case (2)
            select case (msg)
            case (1)
              c1 = (0.0_double,0.0_double)
            end select
            call xcomm_reduce(XSGROUP,MPI_SUM,c1,c2) ; if (error()) exit
            if (i_access(nrestf)) then
              dsize = sizeof_double
              ndata = 2*size(c2)
              call writef(c2,dsize,ndata,x_tagfd(nrestf),iosl)
            end if
          end select

        end do
        if (error()) goto 100

        ! End the COEFFICIENTS block
        if (i_access(nrestf)) call endblock(nrestf)

        ! End the WAVEFUNCTIONS block
        if (i_access(nrestf)) call endblock(nrestf)

100     if (allocated( c1 )) deallocate( c1 )
        if (allocated( c2 )) deallocate( c2 )

        nullify( worm_mb )

        call glean(thy(mv))
        call glean(thy(nrestf))

        if (error("Exit multivector_mod::write_restart_mv")) continue

      end subroutine

! private routines

      subroutine exchange_operator_spair_mv_i(mv,wt,xct,rc,mvo)
        type(multivector_obj)                  :: mv
        real(double), dimension(:), intent(in) :: wt
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        type(multivector_obj)                  :: mvo
!       method: Uses collective communication.
!       errors: Passes errors

!cod$
        integer :: ib1, ib2, grp, nb, ngp, ngt, np, npairs
        integer, dimension(3) :: nd
        integer, dimension(:), allocatable :: index1, index2
        real(double) :: pfm, r0, spin_degeneracy, dk(3)
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: c1, wf1, wf2
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)
        call my(mvo)
        call my(xct)

        nullify( ck )

        worm_mb => wormhole(mv%o%mb)

        call my(worm_mb%lay,lay)
        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = 0.0_double
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv%o%mat,2)

        ngp = size(mv%o%mat,1)
        ngt = size(worm_mb%gridmap,2)
        allocate( vec1(ngt), vec2(ngt), mat1(ngp,np), mat2(ngp,np) )

        c1 => datalink(x_serial_plan(lay))

        wf1 => datalink(worm_mb%isplan1)
        wf2 => datalink(worm_mb%isplan2)

        pfm = (eight_pi/worm_mb%cellvol)/spin_degeneracy

        npairs = size(worm_mb%spair_index,2)
        allocate( index1(npairs), index2(npairs) )
        index1 = worm_mb%spair_index(1,:)
        index2 = worm_mb%spair_index(2,:)

        ! compute the energy and derivative contributions
        do grp = 1,size(worm_mb%first_spair)
          call extract_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index1,ngp,np,mat1,nb,mv%o%mat)
          call extract_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index2,ngp,np,mat2,nb,mv%o%mat)
          call spair_remap(mv%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call spair_remap(mv%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb%spair_participant(grp)) then
            ib1 = index1(worm_mb%my_spair(grp))
            ib2 = index2(worm_mb%my_spair(grp))
            call scatter_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt,worm_mb%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb%isplan1)
            call fft_serial(Q_TO_R,worm_mb%isplan2)
            call kernel_mmm_ccc_i(nd,c1,wf1,wf2)
            call fft_serial(R_TO_Q,x_serial_plan(lay))
            call kernel_mm_cr_i(nd,c1,ck)
            call fft_serial(Q_TO_R,x_serial_plan(lay))
            call kernel_mm_cc_i(nd,wf2,c1)
            call fft_serial(R_TO_Q,worm_mb%isplan2)
            call gather_wf_i(nd,wf2,ngt,worm_mb%gridmap,vec1)
            r0 = -pfm*wt(ib2)
            call kernel_portion_v_rs_i(ngt,r0,vec1)
            if (ib1 == ib2) then
              vec2 = (0.0_double,0.0_double)
            else
              call kernel_mm_ccs_i(nd,wf1,c1)
              call fft_serial(R_TO_Q,worm_mb%isplan1)
              call gather_wf_i(nd,wf1,ngt,worm_mb%gridmap,vec2)
              r0 = -pfm*wt(ib1)
              call kernel_portion_v_rs_i(ngt,r0,vec2)
            end if
          end if
          call spair_remap(mv%o%mb,grp,vec1,mat1) ; if (error()) goto 100
          call spair_remap(mv%o%mb,grp,vec2,mat2) ; if (error()) goto 100
          call augment_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index1,ngp,np,mat1,nb,mvo%o%mat)
          call augment_vectors_indexed_i(worm_mb%first_spair(grp),worm_mb%last_spair(grp),npairs,index2,ngp,np,mat2,nb,mvo%o%mat)
        end do
        mvo%o%g = x_ghost()

100     if (allocated( index1 )) deallocate( index1 )
        if (allocated( index2 )) deallocate( index2 )
        if (associated( ck )) deallocate( ck )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        nullify( c1 )
        nullify( wf1 )
        nullify( wf2 )
        nullify( worm_mb )

        call glean(thy(lay))
        call glean(thy(ck_g))

        call glean(thy(mv))
        call glean(thy(mvo))
        call glean(thy(xct))

200     if (error("Exit multivector_mod::exchange_operator_spair_mv_i")) continue

      end subroutine

      subroutine exchange_operator_lpair_2mv_i(mv1,mv2,wt1,wt2,xct,rc,e,mvo)
        type(multivector_obj)                  :: mv1, mv2
        real(double), dimension(:), intent(in) :: wt1, wt2
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        real(double), intent(inout)            :: e
        type(multivector_obj)                  :: mvo
!       method: Uses collective communication.
!       errors: Passes errors

!cod$
        integer :: ib1, ib2, grp, nb, ngp1, ngp2, ngt1, ngt2, np, npairs
        integer, dimension(3) :: nd
        integer, dimension(:), allocatable :: index1, index2
        real(double) :: dk(3), e1, e2, pfe, pfm, r0, spin_degeneracy
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:), allocatable :: vec1, vec2
        complex(double), dimension(:,:), allocatable :: mat1, mat2
        complex(double), dimension(:,:,:), pointer :: c1, wf1, wf2
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g
        type(multibasis_rep), pointer :: worm_mb1, worm_mb2

        call my(mv1)
        call my(mv2)
        call my(mvo)
        call my(xct)

        nullify( ck )

        worm_mb1 => wormhole(mv1%o%mb)
        worm_mb2 => wormhole(mv2%o%mb)

        call my(worm_mb1%lay,lay)
        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = worm_mb1%kpt - worm_mb2%kpt
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv1%o%mat,2)

        ngp1 = size(mv1%o%mat,1)
        ngt1 = size(worm_mb1%gridmap,2)
        allocate( vec1(ngt1), mat1(ngp1,np) )

        ngp2 = size(mv2%o%mat,1)
        ngt2 = size(worm_mb2%gridmap,2)
        allocate( vec2(ngt2), mat2(ngp2,np) )

        c1 => datalink(x_serial_plan(lay))
        wf1 => datalink(worm_mb1%isplan1)
        wf2 => datalink(worm_mb2%isplan1)

        pfe = (four_pi/worm_mb1%cellvol)/spin_degeneracy
        pfm = (eight_pi/worm_mb1%cellvol)/spin_degeneracy

        npairs = size(worm_mb1%lpair_index,2)
        allocate( index1(npairs), index2(npairs) )
        index1 = worm_mb1%lpair_index(1,:)
        index2 = worm_mb1%lpair_index(2,:)

        ! compute the energy and derivative contributions
        e1 = 0.0_double
        do grp = 1,size(worm_mb1%first_lpair)
          call extract_vectors_indexed_i(worm_mb1%first_lpair(grp),worm_mb1%last_lpair(grp),npairs,index1,ngp1,np,mat1,nb,mv1%o%mat)
          call extract_vectors_indexed_i(worm_mb2%first_lpair(grp),worm_mb2%last_lpair(grp),npairs,index2,ngp2,np,mat2,nb,mv2%o%mat)
          call lpair_remap(mv1%o%mb,grp,mat1,vec1) ; if (error()) goto 100
          call lpair_remap(mv2%o%mb,grp,mat2,vec2) ; if (error()) goto 100
          if (worm_mb1%lpair_participant(grp)) then
            ib1 = index1(worm_mb1%my_lpair(grp))
            ib2 = index2(worm_mb2%my_lpair(grp))
            call scatter_wf_i(nd,wf1,ngt1,worm_mb1%gridmap,vec1)
            call scatter_wf_i(nd,wf2,ngt2,worm_mb2%gridmap,vec2)
            call fft_serial(Q_TO_R,worm_mb1%isplan1)
            call fft_serial(Q_TO_R,worm_mb2%isplan1)
            call kernel_mmm_ccc_i(nd,c1,wf1,wf2)
            call fft_serial(R_TO_Q,x_serial_plan(lay))
            call kernel_mm_cr_i(nd,c1,ck)
            call fft_serial(Q_TO_R,x_serial_plan(lay))
            call kernel_rsum_mmm_ccc_i(nd,c1,wf1,wf2,r0)
            e1 = e1 - pfe*wt1(ib1)*wt2(ib2)*r0
            call kernel_mm_cc_i(nd,c1,wf2)
            wf1 = c1
            call fft_serial(R_TO_Q,worm_mb1%isplan1)
            call gather_wf_i(nd,wf1,ngt1,worm_mb1%gridmap,vec1)
            r0 = -pfm*wt2(ib2)
            call kernel_portion_v_rs_i(ngt1,r0,vec1)
          end if
          call lpair_remap(mv1%o%mb,grp,vec1,mat1) ; if (error()) goto 100
          call augment_vectors_indexed_i(worm_mb1%first_lpair(grp),worm_mb1%last_lpair(grp),npairs,index1,ngp1,np,mat1,nb,mvo%o%mat)
        end do
        call allreduce(KGROUP,MPI_SUM,e1,e2)
        e = e + e2
        mvo%o%g = x_ghost()

100     if (allocated( index1 )) deallocate( index1 )
        if (allocated( index2 )) deallocate( index2 )
        if (associated( ck )) deallocate( ck )
        nullify( c1 )
        nullify( wf1 )
        nullify( wf2 )
        if (allocated( vec1 )) deallocate( vec1 )
        if (allocated( vec2 )) deallocate( vec2 )
        if (allocated( mat1 )) deallocate( mat1 )
        if (allocated( mat2 )) deallocate( mat2 )
        call glean(thy(lay))
        call glean(thy(ck_g))
        nullify( worm_mb1 )
        nullify( worm_mb2 )

        call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(mvo))
        call glean(thy(xct))

        if (error("Exit multivector_mod::exchange_operator_lpair_2mv_i")) continue

      end subroutine

      subroutine exchange_operator_p2p_mv_i(mv,wt,xct,rc,mvo)
        type(multivector_obj)                  :: mv
        real(double), dimension(:), intent(in) :: wt
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        type(multivector_obj)                  :: mvo
!       method: Uses point-to-point communication except to accumulate e.
!       notes: Currently has constraints on np: np = m*nbands where m is an integer .and. np <= nbands*(nbands + 1)/2.
!       errors: Passes errors

!cod$
        integer :: ic, is
        integer :: myp, np, sp, rp, tag
        integer :: grp, fb, lb, nb
        integer :: ngp, ngt
        integer :: row_band, col_band
        integer :: f1, l1, f2, l2
        integer, dimension(3) :: nd
        real(double) :: pfm, r0, spin_degeneracy
        real(double), dimension(3) :: dk
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:,:,:), pointer :: c1
        complex(double), dimension(:,:,:), pointer :: row_wf, col_wf
        complex(double), dimension(:,:,:), pointer :: row_owf
        complex(double), dimension(:), allocatable :: row_vec, col_vec
        complex(double), dimension(:), allocatable :: row_ovec, col_ovec, tmp_ovec
        complex(double), dimension(:), allocatable :: col_dble, tmp_dble
        complex(double), dimension(:,:), allocatable :: mat
        type(multibasis_rep), pointer :: worm_mb
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g

        call my(mv)
        call my(mvo)
        call my(xct)

        nullify( ck )
        nullify( row_owf )

        worm_mb => wormhole(mv%o%mb)

        call my(worm_mb%lay,lay)
        nd = x_dims(lay)

        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = 0.0_double
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        myp = mpi_myproc(KGROUP)
        np = mpi_nprocs(KGROUP)

        nb = size(mv%o%mat,2)

        ngp = size(mv%o%mat,1)
        ngt = size(worm_mb%gridmap,2)

        allocate( row_vec(ngt) )
        allocate( col_vec(ngt) )
        allocate( col_ovec(ngt) )

        c1 => datalink(x_serial_plan(lay))
        row_wf => datalink(worm_mb%isplan1)
        col_wf => datalink(worm_mb%isplan2)
        call alloc(row_owf,lay,S_TYPE)

        pfm = (eight_pi/worm_mb%cellvol)/spin_degeneracy

        ! remap the vectors to the first nb processes
        grp = 1
        call band_remap(mv%o%mb,grp,mv%o%mat,row_vec) ; if (error()) goto 100

        ! distribute row vectors to the rest of the processes
        if (associated( worm_mb%r2r_dmap )) then
          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb%r2r_dmap,1)
            tag = is
            sp = worm_mb%r2r_dmap(is,1)
            rp = worm_mb%r2r_dmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_vec,sp,tag) ; if (error()) goto 10
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,row_vec,rp,tag)    ; if (error()) goto 10
            end if
10          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")
        end if

        ! distribute column vectors to the processes
        call start_timer("multivector: p2p communication")
        do is = 1,size(worm_mb%r2c_dmap,1)
          tag = is
          sp = worm_mb%r2c_dmap(is,1)
          rp = worm_mb%r2c_dmap(is,2)
          if (sp /= -1 .and. rp /= -1) then
            if (sp == rp) then
              col_vec = row_vec
            else
              call nonblocking_send(KGROUP,row_vec,sp,tag)   ; if (error()) goto 20
              call blocking_recv(KGROUP,col_vec,rp,tag)      ; if (error()) goto 20
            end if
          elseif (sp /= -1 .or. rp /= -1) then
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_vec,sp,tag)   ; if (error()) goto 20
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,col_vec,rp,tag)      ; if (error()) goto 20
            end if
          end if
20        call barrier(KGROUP)
!          call sync_kgroup_process_errors() ; if (error()) goto 100
        end do
        call stop_timer("multivector: p2p communication")

        ! initialize the calculation counter
        ic = 1

        ! calculate results
        call start_timer("multivector: p2p calculation")
        row_band = worm_mb%rc_band_index(ic,1)
        col_band = worm_mb%rc_band_index(ic,2)
        call scatter_wf_i(nd,row_wf,ngt,worm_mb%gridmap,row_vec)
        call fft_serial(Q_TO_R,worm_mb%isplan1)
        if (col_band == row_band) then
          col_wf = row_wf
        else
          call scatter_wf_i(nd,col_wf,ngt,worm_mb%gridmap,col_vec)
          call fft_serial(Q_TO_R,worm_mb%isplan2)
        end if
        call kernel_mmm_ccc_i(nd,c1,row_wf,col_wf)
        call fft_serial(R_TO_Q,x_serial_plan(lay))
        call kernel_mm_cr_i(nd,c1,ck)
        call fft_serial(Q_TO_R,x_serial_plan(lay))
        r0 = -pfm*wt(col_band)
        call kernel_portion_mmm_ccc_rs_i(nd,row_owf,r0,col_wf,c1)
        if (row_band == col_band) then
          col_ovec = (0.0_double,0.0_double)
        else
          call kernel_mmm_ccc_i(nd,col_wf,row_wf,c1)
          call fft_serial(R_TO_Q,worm_mb%isplan2)
          r0 = -pfm*wt(row_band)
          call gather_wf_portion_i(nd,col_wf,ngt,worm_mb%gridmap,r0,col_ovec)
        end if
        call stop_timer("multivector: p2p calculation")

        deallocate( row_vec )

        ! iterate over column shifts
        if (associated( worm_mb%c2c_smap )) then

          ! temporary storage
          f1 = 1
          l1 = ngt
          f2 = f1 + ngt
          l2 = l1 + ngt
          allocate( col_dble(f1:l2) )
          allocate( tmp_dble(f1:l2) )

          ! shift the column vectors
          do is = 1,size(worm_mb%c2c_smap,1)

            call start_timer("multivector: p2p communication")
            tag = is
            sp = worm_mb%c2c_smap(is,1)
            rp = worm_mb%c2c_smap(is,2)
            if (sp /= -1 .and. rp /= -1) then
              tmp_dble(f1:l1) = col_vec
              tmp_dble(f2:l2) = col_ovec
              call nonblocking_send(KGROUP,tmp_dble,sp,tag)    ; if (error()) goto 30
              call blocking_recv(KGROUP,col_dble,rp,tag)       ; if (error()) goto 30
              col_vec  = col_dble(f1:l1)
              col_ovec = col_dble(f2:l2)
            end if
30          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
            call stop_timer("multivector: p2p communication")

            ! calculate and accumulate the results
            call start_timer("multivector: p2p calculation")
            if (worm_mb%c2c_ulst(is)) then
              ic = ic + 1
              row_band = worm_mb%rc_band_index(ic,1)
              col_band = worm_mb%rc_band_index(ic,2)
              call scatter_wf_i(nd,col_wf,ngt,worm_mb%gridmap,col_vec)
              call fft_serial(Q_TO_R,worm_mb%isplan2)
              call kernel_mmm_ccc_i(nd,c1,row_wf,col_wf)
              call fft_serial(R_TO_Q,x_serial_plan(lay))
              call kernel_mm_cr_i(nd,c1,ck)
              call fft_serial(Q_TO_R,x_serial_plan(lay))
              r0 = -pfm*wt(col_band)
              call kernel_csum_portion_mmm_ccc_rs_i(nd,row_owf,r0,col_wf,c1)
              call kernel_mmm_ccc_i(nd,col_wf,row_wf,c1)
              call fft_serial(R_TO_Q,worm_mb%isplan2)
              r0 = -pfm*wt(row_band)
              call gather_wf_sum_portion_i(nd,col_wf,ngt,worm_mb%gridmap,r0,col_ovec)
            end if
            call stop_timer("multivector: p2p calculation")
            call start_timer("multivector: p2p calculation barrier")
            call barrier(KGROUP)
            call stop_timer("multivector: p2p calculation barrier")

          end do

          deallocate( col_dble )
          deallocate( tmp_dble )

        end if

        call barrier(KGROUP)

        deallocate( ck )
        deallocate( col_vec )

        ! Fourier transform row_owf and gather to row_ovec
        call start_timer("multivector: p2p calculation")
        row_wf = row_owf
        call fft_serial(R_TO_Q,worm_mb%isplan1)
        allocate( row_ovec(ngt) )
        call gather_wf_i(nd,row_wf,ngt,worm_mb%gridmap,row_ovec)
        call stop_timer("multivector: p2p calculation")

        allocate( tmp_ovec(ngt) )

        ! collect the column results to the rows
        call start_timer("multivector: p2p communication")
        do is = 1,size(worm_mb%c2r_cmap,1)
          tag = is
          sp = worm_mb%c2r_cmap(is,1)
          rp = worm_mb%c2r_cmap(is,2)
          if (sp /= -1 .and. rp /= -1) then
            if (myp /= rp) then
              call nonblocking_send(KGROUP,col_ovec,sp,tag)    ; if (error()) goto 40
              call blocking_recv(KGROUP,tmp_ovec,rp,tag)       ; if (error()) goto 40
              row_ovec = row_ovec + tmp_ovec
            end if
          elseif (sp /= -1 .or. rp /= -1) then
            if (sp /= -1) then
              call nonblocking_send(KGROUP,col_ovec,sp,tag)    ; if (error()) goto 40
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,tmp_ovec,rp,tag)       ; if (error()) goto 40
              row_ovec = row_ovec + tmp_ovec
            end if
          end if
40        call barrier(KGROUP)
!          call sync_kgroup_process_errors() ; if (error()) goto 100
        end do
        call stop_timer("multivector: p2p communication")

        call barrier(KGROUP)

        deallocate( col_ovec )

        ! collect the row results to the first nb rows
        if (associated( worm_mb%r2r_cmap )) then
          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb%r2r_cmap,1)
            tag = is
            sp = worm_mb%r2r_cmap(is,1)
            rp = worm_mb%r2r_cmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_ovec,sp,tag)    ; if (error()) goto 50
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,tmp_ovec,rp,tag)       ; if (error()) goto 50
              row_ovec = row_ovec + tmp_ovec
            end if
50          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")
        end if

        deallocate( tmp_ovec )

        ! remap the row_ovec's of the first nb processes and add them to mvo
        grp = 1
        allocate( mat(ngp,nb) )
        call band_remap(mv%o%mb,grp,row_ovec,mat) ; if (error()) goto 100
        fb = worm_mb%first_band(grp)
        lb = worm_mb%last_band(grp)
        call augment_vectors_i(fb,lb,ngp,nb,mat,nb,mvo%o%mat)

        ! update the mvo ghost
        mvo%o%g = x_ghost()

100     if (associated( ck )) deallocate( ck )
        nullify( c1 )
        nullify( row_wf )
        nullify( col_wf )
        if (associated( row_owf )) deallocate( row_owf )
        if (allocated( row_vec ))  deallocate( row_vec )
        if (allocated( col_vec ))  deallocate( col_vec )
        if (allocated( row_ovec )) deallocate( row_ovec )
        if (allocated( col_ovec )) deallocate( col_ovec )
        if (allocated( tmp_ovec )) deallocate( tmp_ovec )
        if (allocated( col_dble )) deallocate( col_dble )
        if (allocated( tmp_dble )) deallocate( tmp_dble )
        if (allocated( mat )) deallocate( mat )
        call glean(thy(lay))
        call glean(thy(ck_g))
        nullify( worm_mb )

        call glean(thy(mv))
        call glean(thy(mvo))
        call glean(thy(xct))

        if (error("Exit multivector_mod::exchange_operator_p2p_mv_i")) continue

      end subroutine

      subroutine exchange_operator_p2p_2_mv_i(mv,wt,xct,rc,mvo)
        type(multivector_obj)                  :: mv
        real(double), dimension(:), intent(in) :: wt
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        type(multivector_obj)                  :: mvo
!       method: Uses point-to-point communication to accumulate mvo.
!       notes: Has constraints on np: np = m*nb where m is an integer .and. np <= nb*nb.
!       errors: Passes errors

!cod$
        integer :: ic, is
        integer :: myp, np, sp, rp, tag
        integer :: grp, fb, lb, nb
        integer :: ngp, ngt
        integer :: row_band, col_band
        integer, dimension(3) :: nd
        real(double) :: pfm, r0, spin_degeneracy
        real(double), dimension(3) :: dk
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:,:,:), pointer :: c1
        complex(double), dimension(:,:,:), pointer :: row_wf
        complex(double), dimension(:,:,:), pointer :: col_wf
        complex(double), dimension(:,:,:), pointer :: row_owf
        complex(double), dimension(:), allocatable :: row_vec
        complex(double), dimension(:), allocatable :: col_vec
        complex(double), dimension(:), allocatable :: row_ovec
        complex(double), dimension(:), allocatable :: tmp_vec
        complex(double), dimension(:), allocatable :: tmp_ovec
        complex(double), dimension(:,:), allocatable :: mat
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)
        call my(xct)
        call my(mvo)

        nullify( ck )
        nullify( row_owf )

        worm_mb => wormhole(mv%o%mb)

        call my(worm_mb%lay,lay)
        nd = x_dims(lay)

        call my(grid(lay,KGROUP),ck_g)

        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        ! construct the Coulomb kernel
        dk = 0.0_double
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        myp = mpi_myproc(KGROUP)
        np = mpi_nprocs(KGROUP)

        nb = size(mv%o%mat,2)

        ngp = size(mv%o%mat,1)
        ngt = size(worm_mb%gridmap,2)

        c1 => datalink(x_serial_plan(lay))
        row_wf => datalink(worm_mb%isplan1)
        col_wf => datalink(worm_mb%isplan2)
        call alloc(row_owf,lay,S_TYPE)

        pfm = (eight_pi/worm_mb%cellvol)/spin_degeneracy

        allocate( row_vec(ngt) )
        allocate( col_vec(ngt) )

        ! remap the vectors to the first nb processes
        grp = 1
        call band_remap(mv%o%mb,grp,mv%o%mat,row_vec) ; if (error()) goto 100

        ! distribute row vectors to the rest of the processes
        if (associated( worm_mb%r2r_dmap )) then
          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb%r2r_dmap,1)
            tag = is
            sp = worm_mb%r2r_dmap(is,1)
            rp = worm_mb%r2r_dmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_vec,sp,tag) ; if (error()) goto 10
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,row_vec,rp,tag)    ; if (error()) goto 10
            end if
10          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")
        end if

        ! distribute column vectors to the processes
        call start_timer("multivector: p2p communication")
        tag = 1
        sp = worm_mb%r2c_dmap(1,1)
        rp = worm_mb%r2c_dmap(1,2)
        if (sp == myp .and. rp == myp) then
          col_vec = row_vec
        else
          call nonblocking_send(KGROUP,row_vec,sp,tag)   ; if (error()) goto 20
          call blocking_recv(KGROUP,col_vec,rp,tag)      ; if (error()) goto 20
        end if
20      call barrier(KGROUP)
!        call sync_kgroup_process_errors() ; if (error()) goto 100
        call stop_timer("multivector: p2p communication")

        ! initialize the calculation counter
        ic = 1

        ! calculate results
        call start_timer("multivector: p2p calculation")
        row_band = worm_mb%rc_band_index(ic,1)
        col_band = worm_mb%rc_band_index(ic,2)
        call scatter_wf_i(nd,row_wf,ngt,worm_mb%gridmap,row_vec)
        call fft_serial(Q_TO_R,worm_mb%isplan1)
        if (col_band == row_band) then
          col_wf = row_wf
        else
          call scatter_wf_i(nd,col_wf,ngt,worm_mb%gridmap,col_vec)
          call fft_serial(Q_TO_R,worm_mb%isplan2)
        end if
        call kernel_mmm_ccc_i(nd,c1,row_wf,col_wf)
        call fft_serial(R_TO_Q,x_serial_plan(lay))
        call kernel_mm_cr_i(nd,c1,ck)
        call fft_serial(Q_TO_R,x_serial_plan(lay))
        r0 = -pfm*wt(col_band)
        call kernel_portion_mmm_ccc_rs_i(nd,row_owf,r0,col_wf,c1)
        call stop_timer("multivector: p2p calculation")

        deallocate( row_vec )

        ! iterate over column shifts
        if (associated( worm_mb%c2c_smap )) then

          allocate( tmp_vec(ngt) )

          do is = 1,size(worm_mb%c2c_smap,1)

            ! shift the column vectors
            call start_timer("multivector: p2p communication")
            tag = is
            sp = worm_mb%c2c_smap(is,1)
            rp = worm_mb%c2c_smap(is,2)
            tmp_vec = col_vec
            if (sp /= -1 .and. rp /= -1) then
              call nonblocking_send(KGROUP,tmp_vec,sp,tag)    ; if (error()) goto 30
              call blocking_recv(KGROUP,col_vec,rp,tag)       ; if (error()) goto 30
            end if
30          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
            call stop_timer("multivector: p2p communication")

            ! calculate and accumulate results
            call start_timer("multivector: p2p calculation")
            if (worm_mb%c2c_ulst(is)) then
              ic = ic + 1
              row_band = worm_mb%rc_band_index(ic,1)
              col_band = worm_mb%rc_band_index(ic,2)
              call scatter_wf_i(nd,col_wf,ngt,worm_mb%gridmap,col_vec)
              call fft_serial(Q_TO_R,worm_mb%isplan2)
              call kernel_mmm_ccc_i(nd,c1,row_wf,col_wf)
              call fft_serial(R_TO_Q,x_serial_plan(lay))
              call kernel_mm_cr_i(nd,c1,ck)
              call fft_serial(Q_TO_R,x_serial_plan(lay))
              r0 = -pfm*wt(col_band)
              call kernel_csum_portion_mmm_ccc_rs_i(nd,row_owf,r0,col_wf,c1)
            end if
            call stop_timer("multivector: p2p calculation")
            call start_timer("multivector: p2p calculation barrier")
            call barrier(KGROUP)
            call stop_timer("multivector: p2p calculation barrier")

          end do

          deallocate( tmp_vec )

        end if

        call barrier(KGROUP)

        deallocate( ck )
        deallocate( col_vec )

        ! Fourier transform row_owf and gather to row_ovec
        call start_timer("multivector: p2p calculation")
        row_wf = row_owf
        call fft_serial(R_TO_Q,worm_mb%isplan1)
        allocate( row_ovec(ngt) )
        call gather_wf_i(nd,row_wf,ngt,worm_mb%gridmap,row_ovec)
        call stop_timer("multivector: p2p calculation")

        ! collect the row results to the first nb rows
        if (associated( worm_mb%r2r_cmap )) then

          allocate( tmp_ovec(ngt) )

          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb%r2r_cmap,1)
            tag = is
            sp = worm_mb%r2r_cmap(is,1)
            rp = worm_mb%r2r_cmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_ovec,sp,tag)    ; if (error()) goto 40
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,tmp_ovec,rp,tag)       ; if (error()) goto 40
              row_ovec = row_ovec + tmp_ovec
            end if
40          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")

          deallocate( tmp_ovec )

        end if

        allocate( mat(ngp,nb) )

        ! remap the row_ovec's of the first nb processes and add them to mvo
        grp = 1
        call band_remap(mv%o%mb,grp,row_ovec,mat) ; if (error()) goto 100
        fb = worm_mb%first_band(grp)
        lb = worm_mb%last_band(grp)
        call augment_vectors_i(fb,lb,ngp,nb,mat,nb,mvo%o%mat)

        ! update the mvo ghost
        mvo%o%g = x_ghost()

100     if (associated( ck )) deallocate( ck )
        nullify( c1 )
        nullify( row_wf )
        nullify( col_wf )
        if (associated( row_owf )) deallocate( row_owf )
        if (allocated( row_vec ))  deallocate( row_vec )
        if (allocated( col_vec ))  deallocate( col_vec )
        if (allocated( row_ovec )) deallocate( row_ovec )
        if (allocated( tmp_vec )) deallocate( tmp_vec )
        if (allocated( tmp_ovec )) deallocate( tmp_ovec )
        if (allocated( mat )) deallocate( mat )
        call glean(thy(lay))
        call glean(thy(ck_g))
        nullify( worm_mb )

        call glean(thy(mv))
        call glean(thy(xct))
        call glean(thy(mvo))

        if (error("Exit multivector_mod::exchange_operator_p2p_2_mv_i")) continue

      end subroutine

      subroutine exchange_operator_p2p_2_2mv_i(mv1,mv2,wt1,wt2,xct,rc,e,mvo)
        type(multivector_obj)                  :: mv1, mv2
        real(double), dimension(:), intent(in) :: wt1, wt2
        type(xc_type_obj)                      :: xct
        real(double), intent(in)               :: rc
        real(double), intent(inout)            :: e
        type(multivector_obj)                  :: mvo
!       method: Uses point-to-point communication to accumulate e and mvo.
!       notes: Has constraints on np: np = m*nb where m is an integer .and. np <= nb*nb.
!       notes: For a given (mv1,mv2) pair, a single set of processors with communicator KGROUP
!              processes the code below.
!       errors: Passes errors

!cod$
        integer :: ic, is
        integer :: myp, np, sp, rp, tag
        integer :: grp, fb, lb, nb
        integer :: ngp1, ngp2
        integer :: ngt1, ngt2
        integer :: row_band, col_band
        integer, dimension(3) :: nd
        real(double) :: e1, e2, pfe
        real(double) :: pfm, r0, r1, spin_degeneracy
        real(double), dimension(3) :: dk
        real(double), dimension(:,:,:), pointer :: ck
        complex(double), dimension(:,:,:), pointer :: c1
        complex(double), dimension(:,:,:), pointer :: row_wf
        complex(double), dimension(:,:,:), pointer :: col_wf
        complex(double), dimension(:,:,:), pointer :: row_owf
        complex(double), dimension(:), allocatable :: row_vec
        complex(double), dimension(:), allocatable :: col_vec
        complex(double), dimension(:), allocatable :: row_ovec
        complex(double), dimension(:), allocatable :: tmp_vec
        complex(double), dimension(:), allocatable :: tmp_ovec
        complex(double), dimension(:,:), allocatable :: mat
        type(layout_obj) :: lay
        type(grid_obj) :: ck_g
        type(multibasis_rep), pointer :: worm_mb1
        type(multibasis_rep), pointer :: worm_mb2

        call my(mv1)
        call my(mv2)
        call my(xct)
        call my(mvo)

        nullify( ck )

        worm_mb1 => wormhole(mv1%o%mb)
        worm_mb2 => wormhole(mv2%o%mb)

        call my(worm_mb1%lay,lay)

        myp = mpi_myproc(KGROUP)

        ! parameters that are the same for both k-points
        nd = x_dims(lay)
        np = mpi_nprocs(KGROUP)
        nb = size(mv1%o%mat,2)

        ! paramaters that may be different for the two k-points
        ngp1 = size(mv1%o%mat,1)
        ngp2 = size(mv2%o%mat,1)
        ngt1 = size(worm_mb1%gridmap,2)
        ngt2 = size(worm_mb2%gridmap,2)

        ! allocate pointers for the on-processor data
        c1 => datalink(x_serial_plan(lay))
        row_wf =>  datalink(worm_mb1%isplan1)
        col_wf =>  datalink(worm_mb2%isplan1)
        row_owf => datalink(worm_mb1%isplan2)

        ! construct the Coulomb kernel
        dk = worm_mb1%kpt - worm_mb2%kpt
        call my(grid(lay,KGROUP),ck_g)
        call construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g)  ; if (error()) goto 100
        call take(ck,ck_g,RS_KIND)                           ; if (error()) goto 100

        ! a.k.a. maximum occupation of each band
        spin_degeneracy = 2.0_double/real(mpi_nsgroups(),double)

        pfe = (four_pi/worm_mb1%cellvol)/spin_degeneracy
        pfm = (eight_pi/worm_mb1%cellvol)/spin_degeneracy

        ! allocate space for the row vectors
        allocate( row_vec(ngt1) )

        ! remap the k-point 1 distributed vectors to serial (row) vectors on the first nb processes
        grp = 1
        call band_remap(mv1%o%mb,grp,mv1%o%mat,row_vec)  ; if (error()) goto 100

        ! distribute copies of the row vectors to the rest of the processes
        if (associated( worm_mb1%r2r_dmap )) then
          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb1%r2r_dmap,1)
            tag = is
            sp = worm_mb1%r2r_dmap(is,1)
            rp = worm_mb1%r2r_dmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_vec,sp,tag)  ; if (error()) goto 10
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,row_vec,rp,tag)     ; if (error()) goto 10
            end if
10          call barrier(KGROUP)
!            call sync_kgroup_process_errors()  ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")
        end if

        ! allocate space for the column vectors
        allocate( col_vec(ngt2) )
        allocate( tmp_vec(ngt2) )

        ! remap the k-point 2 distributed vectors to serial (column) vectors on the first nb processes
        grp = 1
        call band_remap(mv2%o%mb,grp,mv2%o%mat,tmp_vec)  ; if (error()) goto 100

        ! distribute copies of the column vectors to the rest of the processes
        if (associated( worm_mb2%r2r_dmap )) then
          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb2%r2r_dmap,1)
            tag = is
            sp = worm_mb2%r2r_dmap(is,1)
            rp = worm_mb2%r2r_dmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,tmp_vec,sp,tag)  ; if (error()) goto 15
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,tmp_vec,rp,tag)     ; if (error()) goto 15
            end if
15          call barrier(KGROUP)
!            call sync_kgroup_process_errors()  ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")
        end if

        ! arrange the column vectors on the nb blocks of processes
        call start_timer("multivector: p2p communication")
        tag = 1
        sp = worm_mb2%r2c_dmap(1,1)
        rp = worm_mb2%r2c_dmap(1,2)
        if (sp == myp .and. rp == myp) then
          col_vec = tmp_vec
        else
          call nonblocking_send(KGROUP,tmp_vec,sp,tag)   ; if (error()) goto 20
          call blocking_recv(KGROUP,col_vec,rp,tag)      ; if (error()) goto 20
        end if
20      call barrier(KGROUP)
!        call sync_kgroup_process_errors() ; if (error()) goto 100
        call stop_timer("multivector: p2p communication")

        ! deallocate the temporary space
        deallocate( tmp_vec )

        ! initialize the process energy
        e1 = 0.0_double

        ! initialize the calculation counter
        ic = 1

        ! calculate results
        call start_timer("multivector: p2p calculation")
        row_band = worm_mb1%rc_band_index(ic,1)
        col_band = worm_mb2%rc_band_index(ic,2)
        call scatter_wf_i(nd,row_wf,ngt1,worm_mb1%gridmap,row_vec)
        call scatter_wf_i(nd,col_wf,ngt2,worm_mb2%gridmap,col_vec)
        call fft_serial(Q_TO_R,worm_mb1%isplan1)
        call fft_serial(Q_TO_R,worm_mb2%isplan1)
        call kernel_mmm_ccc_i(nd,c1,row_wf,col_wf)
        call fft_serial(R_TO_Q,x_serial_plan(lay))
        call kernel_mm_cr_i(nd,c1,ck)
        call fft_serial(Q_TO_R,x_serial_plan(lay))
        call kernel_rsum_mmm_ccc_i(nd,c1,row_wf,col_wf,r0)
        e1 = e1 - pfe*wt1(row_band)*wt2(col_band)*r0
        r1 = -pfm*wt2(col_band)
        call kernel_portion_mmm_ccc_rs_i(nd,row_owf,r1,col_wf,c1)
        call stop_timer("multivector: p2p calculation")

        deallocate( row_vec )

        ! iterate over column shifts
        if (associated( worm_mb2%c2c_smap )) then

          ! allocate temporary space
          allocate( tmp_vec(ngt2) )

          do is = 1,size(worm_mb2%c2c_smap,1)

            ! shift the column vectors
            call start_timer("multivector: p2p communication")
            tag = is
            sp = worm_mb2%c2c_smap(is,1)
            rp = worm_mb2%c2c_smap(is,2)
            tmp_vec = col_vec
            if (sp /= -1 .and. rp /= -1) then
              call nonblocking_send(KGROUP,tmp_vec,sp,tag)    ; if (error()) goto 30
              call blocking_recv(KGROUP,col_vec,rp,tag)       ; if (error()) goto 30
            end if
30          call barrier(KGROUP)
!            call sync_kgroup_process_errors() ; if (error()) goto 100
            call stop_timer("multivector: p2p communication")

            ! calculate and accumulate results
            call start_timer("multivector: p2p calculation")
            if (worm_mb1%c2c_ulst(is)) then
              ic = ic + 1
              row_band = worm_mb1%rc_band_index(ic,1)
              col_band = worm_mb2%rc_band_index(ic,2)
              call scatter_wf_i(nd,col_wf,ngt2,worm_mb2%gridmap,col_vec)
              call fft_serial(Q_TO_R,worm_mb2%isplan1)
              call kernel_mmm_ccc_i(nd,c1,row_wf,col_wf)
              call fft_serial(R_TO_Q,x_serial_plan(lay))
              call kernel_mm_cr_i(nd,c1,ck)
              call fft_serial(Q_TO_R,x_serial_plan(lay))
              call kernel_rsum_mmm_ccc_i(nd,c1,row_wf,col_wf,r0)
              e1 = e1 - pfe*wt1(row_band)*wt2(col_band)*r0
              r1 = -pfm*wt2(col_band)
              call kernel_csum_portion_mmm_ccc_rs_i(nd,row_owf,r1,col_wf,c1)
            end if
            call stop_timer("multivector: p2p calculation")
            call start_timer("multivector: p2p calculation barrier")
            call barrier(KGROUP)
            call stop_timer("multivector: p2p calculation barrier")

          end do

          deallocate( tmp_vec )

        end if

        call barrier(KGROUP)

        deallocate( ck )
        deallocate( col_vec )

        ! accumulate the energy contributions
        call allreduce(KGROUP,MPI_SUM,e1,e2)
        e = e + e2

        ! Fourier transform row_owf and gather to row_ovec
        call start_timer("multivector: p2p calculation")
        allocate( row_ovec(ngt1) )
        call fft_serial(R_TO_Q,worm_mb1%isplan2)
        call gather_wf_i(nd,row_owf,ngt1,worm_mb1%gridmap,row_ovec)
        call stop_timer("multivector: p2p calculation")

        ! accumulate the row-vector results to the first nb processes
        if (associated( worm_mb1%r2r_cmap )) then

          allocate( tmp_ovec(ngt1) )

          call start_timer("multivector: p2p communication")
          do is = 1,size(worm_mb1%r2r_cmap,1)
            tag = is
            sp = worm_mb1%r2r_cmap(is,1)
            rp = worm_mb1%r2r_cmap(is,2)
            if (sp /= -1) then
              call nonblocking_send(KGROUP,row_ovec,sp,tag)  ; if (error()) goto 40
            end if
            if (rp /= -1) then
              call blocking_recv(KGROUP,tmp_ovec,rp,tag)     ; if (error()) goto 40
              row_ovec = row_ovec + tmp_ovec
            end if
40          call barrier(KGROUP)
!            call sync_kgroup_process_errors()  ; if (error()) goto 100
          end do
          call stop_timer("multivector: p2p communication")

          deallocate( tmp_ovec )

        end if

        ! remap the serial row_ovec's on the first nb processes and add them to mvo
        allocate( mat(ngp1,nb) )
        grp = 1
        call band_remap(mv1%o%mb,grp,row_ovec,mat)  ; if (error()) goto 100
        fb = worm_mb1%first_band(grp)
        lb = worm_mb1%last_band(grp)
        call augment_vectors_i(fb,lb,ngp1,nb,mat,nb,mvo%o%mat)
        deallocate( mat )

        ! update the mvo ghost
        mvo%o%g = x_ghost()

100     if (associated( ck )) deallocate( ck )
        nullify( c1 )
        nullify( row_wf )
        nullify( col_wf )
        nullify( row_owf )
        if (allocated( row_vec ))  deallocate( row_vec )
        if (allocated( col_vec ))  deallocate( col_vec )
        if (allocated( row_ovec )) deallocate( row_ovec )
        if (allocated( tmp_vec ))  deallocate( tmp_vec )
        if (allocated( tmp_ovec )) deallocate( tmp_ovec )
        if (allocated( mat )) deallocate( mat )
        call glean(thy(lay))
        call glean(thy(ck_g))
        nullify( worm_mb1 )
        nullify( worm_mb2 )

        call glean(thy(mv1))
        call glean(thy(mv2))
        call glean(thy(xct))
        call glean(thy(mvo))

        if (error("Exit multivector_mod::exchange_operator_p2p_2_2mv_i")) continue

      end subroutine

      subroutine construct_coulomb_kernel_i(lay,xct,rc,dk,ck_g) !, ck_norm)
        type(layout_obj)              :: lay
        type(xc_type_obj)             :: xct
        real(double), intent(in)      :: rc
        real(double), intent(in)      :: dk(3)
        type(grid_obj)                :: ck_g
        real(double), dimension(:,:,:), pointer :: ck, g2, gdk2, gx, gy, gz
        real(double)                            :: four_omega2
        real(double)                            :: ck_norm
        real(double)                            :: dkx, dky, dkz
        real(double), parameter                 :: tol = 1.0e-10_double

        call my(lay)
        call my(xct)
        call my(ck_g)
        
        nullify(ck,g2,gdk2,gx,gy,gz)

        call fmesh(gx,gy,gz,lay,D_TYPE,KGROUP)
        call alloc(g2,lay,D_TYPE,KGROUP)
        g2 = gx**2 + gy**2 + gz**2
        dkx = dk(1)
        dky = dk(2)
        dkz = dk(3)
        call alloc(gdk2,lay,D_TYPE,KGROUP)
        gdk2 = ((gx + dkx)**2 + (gy + dky)**2 + (gz + dkz)**2)

        call alloc(ck,lay,D_TYPE,KGROUP)

        select case (x_coulomb_kernel(xct))
        case (CK_NORMAL)
           where (gdk2 < tol) 
              ck = 0.0_double
           elsewhere (x_cutoff(lay) < g2)
              ck = 0.0_double
           elsewhere
              ck = 1.0_double/gdk2
           end where
           ck_norm = 0.0_double
        case (CK_ATTENUATED)
           where (gdk2 < tol) 
              ck = 0.0_double
           elsewhere (x_cutoff(lay) < g2)
              ck = 0.0_double
           elsewhere
              ck = (1.0_double - cos(sqrt(gdk2)*rc))/gdk2
           end where
           ck_norm = 0.5_double*rc**2
        case(CK_SCREENED)
           four_omega2 = 4.0_double*x_omega_orb(xct)**2

           ! In this limit omega->0, the short range contribution to Ex is zero.
           if (four_omega2 < tol) then
              ck = 0.0_double
              ck_norm = 0.0_double
           else  ! omega > 0 so contruct the screened coulomb kernel.
              where (gdk2 < tol) 
                 ck = 1.0/four_omega2
              elsewhere (x_cutoff(lay) < g2)
                 ck = 0.0_double
              elsewhere
                 ck = (1.0_double - exp(-gdk2/four_omega2))/gdk2
              end where
              ck_norm = 1.0_double/four_omega2
           end if
        case default 
           if (error(.true.,"Unrecognized coulomb kernel type")) goto 100
        end select

        deallocate( g2, gdk2, gx, gy, gz )

        call put(ck,ck_g,RD_KIND)

        if (dk .in. nbhd((/0.0_double,0.0_double,0.0_double/),tol)) then
           call set_normalization(ck_g,ck_norm);  if (error()) goto 100
        end if

        call glean(thy(ck_g))
        call glean(thy(xct))
        call glean(thy(lay))

100     if (error("Exit multivector_mod::construct_coulomb_kernel_i")) continue 

      end subroutine

      subroutine zeros_init_i(mvr)
        type(multivector_rep) :: mvr

        integer :: nb, ng
        type(multibasis_rep), pointer :: worm_mb

        worm_mb => wormhole(mvr%mb)

        nb = x_n_bands(mvr%mb)
        ng = size(worm_mb%gpt,1)

        allocate( mvr%mat(ng,nb) )
        mvr%mat = (0.0_double,0.0_double)

        nullify( worm_mb )

        if (error("Exit multivector_mod::zeros_init_i")) continue

      end subroutine

      subroutine random_init_i(mvr)
        type(multivector_rep) :: mvr

        integer :: i, ib, ig, nb, ng, seed
        real(double), parameter :: tol_nbhd = 1.0e-10_double
        real(double) :: gk2, r1, r2
        real(double), dimension(:), allocatable :: gk2i
        type(multibasis_rep), pointer :: worm_mb

        worm_mb => wormhole(mvr%mb)

        ng = size(worm_mb%gpt,1)
        nb = x_n_bands(mvr%mb)

        allocate( gk2i(ng) )
        if (worm_mb%kpt .in. nbhd((/0.0_double,0.0_double,0.0_double/),tol_nbhd)) then
          do ig = 1,ng
            gk2 = (worm_mb%gpt(ig,1) + worm_mb%kpt(1))**2 + &
                  (worm_mb%gpt(ig,2) + worm_mb%kpt(2))**2 + &
                  (worm_mb%gpt(ig,3) + worm_mb%kpt(3))**2
            if (gk2 .in. nbhd(0.0_double,tol_nbhd)) then
              gk2i(ig) = 1.0_double
            else
              gk2i(ig) = 1.0_double/gk2
            end if
          end do
        else
          do ig = 1,ng
            gk2 = (worm_mb%gpt(ig,1) + worm_mb%kpt(1))**2 + &
                  (worm_mb%gpt(ig,2) + worm_mb%kpt(2))**2 + &
                  (worm_mb%gpt(ig,3) + worm_mb%kpt(3))**2
            gk2i(ig) = 1.0_double/gk2
          end do
        end if

        seed = mpi_myproc(KGROUP) + 1
        do i = 1,50
          r1 = random(seed)
        end do

        allocate( mvr%mat(ng,nb) )
        do ib = 1,nb
          do ig = 1,ng
            r1 = 2.0_double*random(seed) - 1.0_double
            r2 = 2.0_double*random(seed) - 1.0_double
            mvr%mat(ig,ib) = gk2i(ig)*cmplx(r1,r2,double)
          end do
        end do

        if (allocated( gk2i )) deallocate( gk2i )

        nullify( worm_mb )

        if (error("Exit multivector_mod::random_init_i")) continue

      end subroutine

      subroutine diagnostic_init_i(mvr)
        type(multivector_rep) :: mvr

        integer :: ib, ig, nb, ng
        type(multibasis_rep), pointer :: worm_mb

        worm_mb => wormhole(mvr%mb)

        nb = x_n_bands(mvr%mb)
        ng = size(worm_mb%gpt,1)

        allocate( mvr%mat(ng,nb) )
        do ib = 1,nb
          do ig = 1,ng
            mvr%mat(ig,ib) = cmplx(exp(real(ib+2,double))*worm_mb%gpt(ig,1)**1  &
                                            + sqrt(real(ib+4,double))*worm_mb%gpt(ig,2)**2  &
                                            + sin(real(ib+1,double))*worm_mb%gpt(ig,3)**3, &
                                       6*sin(real(ib-5,double))*worm_mb%gpt(ig,1)**2  &
                                            + 2*exp(real(ib+2,double))*worm_mb%gpt(ig,2)**3  &
                                            + 11*sqrt(real(ib,double))*worm_mb%gpt(ig,3)**1,double)
          end do
        end do

        nullify( worm_mb )

        if (error("Exit multivector_mod::diagnostic_init_i")) continue

      end subroutine

      subroutine kernel_portion_v_rs_i(ng,r,v)
        integer :: ng
        real(double) :: r
        complex(double), dimension(ng) :: v

        integer :: ig

        do ig = 1,ng
          v(ig) = r*v(ig)
        end do

      end subroutine

      subroutine kernel_portion_mv_rs_i(ng,nb,r,mat)
        integer :: nb, ng
        real(double) :: r
        complex(double), dimension(ng,nb) :: mat

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat(ig,ib) = r*mat(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_portion_mv_cs_i(ng,nb,c,mat)
        integer :: nb, ng
        complex(double) :: c
        complex(double), dimension(ng,nb) :: mat

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat(ig,ib) = c*mat(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_portion_mv_rv_i(ng,nb,rv,mat)
        integer :: nb, ng
        real(double), dimension(nb) :: rv
        complex(double), dimension(ng,nb) :: mat

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat(ig,ib) = rv(ib)*mat(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_portion_mv_cv_i(ng,nb,cv,mat)
        integer :: nb, ng
        complex(double), dimension(nb) :: cv
        complex(double), dimension(ng,nb) :: mat

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat(ig,ib) = cv(ib)*mat(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_combine_2mv_rs_i(ng,nb,r1,mat1,r2,mat2)
        integer :: nb, ng
        real(double) :: r1, r2
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat1(ig,ib) = r1*mat1(ig,ib) + r2*mat2(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_combine_2mv_cs_i(ng,nb,c1,mat1,c2,mat2)
        integer :: nb, ng
        complex(double) :: c1, c2
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat1(ig,ib) = c1*mat1(ig,ib) + c2*mat2(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_combine_2mv_rv_i(ng,nb,rv1,mat1,rv2,mat2)
        integer :: nb, ng
        real(double), dimension(nb) :: rv1, rv2
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat1(ig,ib) = rv1(ib)*mat1(ig,ib) + rv2(ib)*mat2(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_combine_2mv_cv_i(ng,nb,cv1,mat1,cv2,mat2)
        integer :: nb, ng
        complex(double), dimension(nb) :: cv1, cv2
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat1(ig,ib) = cv1(ib)*mat1(ig,ib) + cv2(ib)*mat2(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_combine_3mv_rv_i(ng,nb,rv1,mat1,rv2,mat2,mat3)
        integer :: nb, ng
        real(double), dimension(nb) :: rv1, rv2
        complex(double), dimension(ng,nb) :: mat1, mat2, mat3

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat3(ig,ib) = rv1(ib)*mat1(ig,ib) + rv2(ib)*mat2(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_combine_3mv_cv_i(ng,nb,cv1,mat1,cv2,mat2,mat3)
        integer :: nb, ng
        complex(double), dimension(nb) :: cv1, cv2
        complex(double), dimension(ng,nb) :: mat1, mat2, mat3

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            mat3(ig,ib) = cv1(ib)*mat1(ig,ib) + cv2(ib)*mat2(ig,ib)
          end do
        end do

      end subroutine

      subroutine kernel_residual_mv_i(ng,nb,v_mat,hv_mat,e,r_mat)
        integer :: nb, ng
        complex(double), dimension(nb) :: e
        complex(double), dimension(ng,nb) :: v_mat, hv_mat, r_mat

        integer :: ib, ig

        do ib = 1,nb
          do ig = 1,ng
            r_mat(ig,ib) = e(ib)*v_mat(ig,ib) - hv_mat(ig,ib)
          end do
        end do

      end subroutine

      subroutine decompose_small_i(mv,sd,mode,rsa,b1)
        type(multivector_obj) :: mv
        real(double), dimension(:,:), intent(in) :: sd
        character(line_len), intent(in) :: mode
        real(double), dimension(:,:,:), intent(out) :: rsa
        integer :: b1

        character(1), parameter :: transa = 't'
        character(1), parameter :: transb = 'n'
        integer, parameter :: nc = 9, nr = 20
        integer :: ig, ir, is, nb, ng, ns
        real(double) :: dr, r, wt
        real(double), dimension(3) :: gk
        real(double), dimension(0:3) :: rm
        real(double), allocatable, dimension(:) :: gkm, gkr, sbr
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), allocatable, dimension(:) :: c1
        complex(double), allocatable, dimension(:,:) :: cm2, cm2l, tmp, ylm
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)

        worm_mb => wormhole(mv%o%mb)

        nb = size(rsa,2)
        ns = size(rsa,3)
        ng = size(worm_mb%gpt,1)

        allocate ( gkm(ng), ylm(ng,nc) )
        select case (mode)
        case ("L", "l", "LM", "lm")
          do ig = 1,ng
            gk(1) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = sqrt(gk(1)**2 + gk(2)**2 + gk(3)**2)
            ylm(ig,1) = spherical_harmonic(0, 0,gk)
            ylm(ig,2) = spherical_harmonic(1,-1,gk)
            ylm(ig,3) = spherical_harmonic(1, 0,gk)
            ylm(ig,4) = spherical_harmonic(1,+1,gk)
            ylm(ig,5) = spherical_harmonic(2,-2,gk)
            ylm(ig,6) = spherical_harmonic(2,-1,gk)
            ylm(ig,7) = spherical_harmonic(2, 0,gk)
            ylm(ig,8) = spherical_harmonic(2,+1,gk)
            ylm(ig,9) = spherical_harmonic(2,+2,gk)
          end do
        case ("XYZ", "xyz")
          do ig = 1,ng
            gk(1) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = sqrt(gk(1)**2 + gk(2)**2 + gk(3)**2)
            ylm(ig,1) = real_spherical_harmonic(0, 0,gk)
            ylm(ig,2) = real_spherical_harmonic(1,+1,gk)
            ylm(ig,3) = real_spherical_harmonic(1,-1,gk)
            ylm(ig,4) = real_spherical_harmonic(1, 0,gk)
            ylm(ig,5) = real_spherical_harmonic(2,-2,gk)
            ylm(ig,6) = real_spherical_harmonic(2,+1,gk)
            ylm(ig,7) = real_spherical_harmonic(2,-1,gk)
            ylm(ig,8) = real_spherical_harmonic(2,+2,gk)
            ylm(ig,9) = real_spherical_harmonic(2, 0,gk)
          end do
        end select

        rm(0) = (2.0_double/45.0_double)*32.0_double*(four_pi)**2/worm_mb%cellvol    ! Bode's rule: nr must be 4 x integer
        rm(1) = (2.0_double/45.0_double)*14.0_double*(four_pi)**2/worm_mb%cellvol
        rm(2) = (2.0_double/45.0_double)*32.0_double*(four_pi)**2/worm_mb%cellvol
        rm(3) = (2.0_double/45.0_double)*12.0_double*(four_pi)**2/worm_mb%cellvol
        allocate ( c1(ng), gkr(ng), sbr(ng), tmp(ng,nc) )
        allocate ( cm2(nc,nb), cm2l(nc,nb) )
        rsa = 0.0_double
        do is = 1,ns
          c1 = exp((0.0_double,1.0_double)*(worm_mb%gpt(:,1)*sd(1,is) + worm_mb%gpt(:,2)*sd(2,is) + worm_mb%gpt(:,3)*sd(3,is)))
          dr = sd(4,is)/real(nr,double)
          do ir = 1,nr
            r = dr*real(ir,double)
            gkr = gkm*r
            tmp(:,1) = ylm(:,1)*c1*spherical_bessel(gkr,0,.true.)
            sbr = spherical_bessel(gkr,1,.true.)*r
            tmp(:,2) = ylm(:,2)*c1*sbr
            tmp(:,3) = ylm(:,3)*c1*sbr
            tmp(:,4) = ylm(:,4)*c1*sbr
            sbr = spherical_bessel(gkr,2,.true.)*r**2
            tmp(:,5) = ylm(:,5)*c1*sbr
            tmp(:,6) = ylm(:,6)*c1*sbr
            tmp(:,7) = ylm(:,7)*c1*sbr
            tmp(:,8) = ylm(:,8)*c1*sbr
            tmp(:,9) = ylm(:,9)*c1*sbr
            call start_timer("multivector: zgemm")
            call zgemm(transa,transb,nc,nb,ng,alpha,tmp(1,1),ng,mv%o%mat(1,b1),ng,beta,cm2l(1,1),nc)
            call stop_timer("multivector: zgemm")
            call allreduce(KGROUP,MPI_SUM,cm2l,cm2)
            wt = rm(mod(ir+1,4))*dr*r**2
            if (ir == nr) wt = wt/2.0_double
            rsa(:,:,is) = rsa(:,:,is) + wt*real(cm2*conjg(cm2),double)
          end do
        end do

        if (allocated( gkm )) deallocate( gkm )
        if (allocated( ylm )) deallocate( ylm )
        if (allocated( c1 )) deallocate( c1 )
        if (allocated( gkr )) deallocate( gkr )
        if (allocated( sbr )) deallocate( sbr )
        if (allocated( tmp )) deallocate( tmp )
        if (allocated( cm2 )) deallocate( cm2 )
        if (allocated( cm2l )) deallocate( cm2l )

        nullify( worm_mb )

        call glean(thy(mv))

        if (error("Exit multivector_mod::decompose_small_i")) continue

      end subroutine

      subroutine decompose_medium_i(mv,sd,mode,rsa,b1)
        type(multivector_obj) :: mv
        real(double), dimension(:,:), intent(in) :: sd
        character(line_len), intent(in) :: mode
        real(double), dimension(:,:,:), intent(out) :: rsa
        integer :: b1

        character(1), parameter :: transa = 't'
        character(1), parameter :: transb = 'n'
        integer, parameter :: nc = 9, nr = 20
        integer :: ig, ir, is, nb, ng, ns
        real(double) :: dr, r, wt
        real(double), dimension(3) :: gk
        real(double), dimension(0:3) :: rm
        real(double), allocatable, dimension(:) :: gkm, gkr, sbr
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), allocatable, dimension(:) :: c1
        complex(double), allocatable, dimension(:,:) :: cm2l, tmp, ylm
        complex(double), allocatable, dimension(:,:,:) :: cm3, cm3l
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)

        worm_mb => wormhole(mv%o%mb)

        nb = size(rsa,2)
        ns = size(rsa,3)
        ng = size(worm_mb%gpt,1)

        allocate ( gkm(ng), ylm(ng,nc) )
        select case (mode)
        case ("L", "l", "LM", "lm")
          do ig = 1,ng
            gk(1) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = sqrt(gk(1)**2 + gk(2)**2 + gk(3)**2)
            ylm(ig,1) = spherical_harmonic(0, 0,gk)
            ylm(ig,2) = spherical_harmonic(1,-1,gk)
            ylm(ig,3) = spherical_harmonic(1, 0,gk)
            ylm(ig,4) = spherical_harmonic(1,+1,gk)
            ylm(ig,5) = spherical_harmonic(2,-2,gk)
            ylm(ig,6) = spherical_harmonic(2,-1,gk)
            ylm(ig,7) = spherical_harmonic(2, 0,gk)
            ylm(ig,8) = spherical_harmonic(2,+1,gk)
            ylm(ig,9) = spherical_harmonic(2,+2,gk)
          end do
        case ("XYZ", "xyz")
          do ig = 1,ng
            gk(1) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = sqrt(gk(1)**2 + gk(2)**2 + gk(3)**2)
            ylm(ig,1) = real_spherical_harmonic(0, 0,gk)
            ylm(ig,2) = real_spherical_harmonic(1,+1,gk)
            ylm(ig,3) = real_spherical_harmonic(1,-1,gk)
            ylm(ig,4) = real_spherical_harmonic(1, 0,gk)
            ylm(ig,5) = real_spherical_harmonic(2,-2,gk)
            ylm(ig,6) = real_spherical_harmonic(2,+1,gk)
            ylm(ig,7) = real_spherical_harmonic(2,-1,gk)
            ylm(ig,8) = real_spherical_harmonic(2,+2,gk)
            ylm(ig,9) = real_spherical_harmonic(2, 0,gk)
          end do
        end select

        rm(0) = (2.0_double/45.0_double)*32.0_double*(four_pi)**2/worm_mb%cellvol    ! Bode's rule: nr must be 4 x integer
        rm(1) = (2.0_double/45.0_double)*14.0_double*(four_pi)**2/worm_mb%cellvol
        rm(2) = (2.0_double/45.0_double)*32.0_double*(four_pi)**2/worm_mb%cellvol
        rm(3) = (2.0_double/45.0_double)*12.0_double*(four_pi)**2/worm_mb%cellvol
        allocate ( c1(ng), gkr(ng), sbr(ng), tmp(ng,nc) )
        allocate ( cm2l(nc,nb), cm3(nc,nb,nr), cm3l(nc,nb,nr) )
        rsa = 0.0_double
        do is = 1,ns
          c1 = exp((0.0_double,1.0_double)*(worm_mb%gpt(:,1)*sd(1,is) + worm_mb%gpt(:,2)*sd(2,is) + worm_mb%gpt(:,3)*sd(3,is)))
          dr = sd(4,is)/real(nr,double)
          do ir = 1,nr
            r = dr*real(ir,double)
            gkr = gkm*r
            tmp(:,1) = ylm(:,1)*c1*spherical_bessel(gkr,0,.true.)
            sbr = spherical_bessel(gkr,1,.true.)*r
            tmp(:,2) = ylm(:,2)*c1*sbr
            tmp(:,3) = ylm(:,3)*c1*sbr
            tmp(:,4) = ylm(:,4)*c1*sbr
            sbr = spherical_bessel(gkr,2,.true.)*r**2
            tmp(:,5) = ylm(:,5)*c1*sbr
            tmp(:,6) = ylm(:,6)*c1*sbr
            tmp(:,7) = ylm(:,7)*c1*sbr
            tmp(:,8) = ylm(:,8)*c1*sbr
            tmp(:,9) = ylm(:,9)*c1*sbr
            call start_timer("multivector: zgemm")
            call zgemm(transa,transb,nc,nb,ng,alpha,tmp(1,1),ng,mv%o%mat(1,b1),ng,beta,cm2l(1,1),nc)
            call stop_timer("multivector: zgemm")
            cm3l(:,:,ir) = cm2l
          end do
          call allreduce(KGROUP,MPI_SUM,cm3l,cm3)
          do ir = 1,nr
            r = dr*real(ir,double)
            wt = rm(mod(ir+1,4))*dr*r**2
            if (ir == nr) wt = wt/2.0_double
            rsa(:,:,is) = rsa(:,:,is) + wt*real(cm3(:,:,ir)*conjg(cm3(:,:,ir)),double)
          end do
        end do

        if (allocated( gkm )) deallocate( gkm )
        if (allocated( ylm )) deallocate( ylm )
        if (allocated( c1 )) deallocate( c1 )
        if (allocated( gkr )) deallocate( gkr )
        if (allocated( sbr )) deallocate( sbr )
        if (allocated( tmp )) deallocate( tmp )
        if (allocated( cm2l )) deallocate( cm2l )
        if (allocated( cm3 )) deallocate( cm3 )
        if (allocated( cm3l )) deallocate( cm3l )

        nullify( worm_mb )

        call glean(thy(mv))

        if (error("Exit multivector_mod::decompose_medium_i")) continue

      end subroutine

      subroutine decompose_large_i(mv,sd,mode,rsa,b1)
        type(multivector_obj) :: mv
        real(double), dimension(:,:), intent(in) :: sd
        character(line_len), intent(in) :: mode
        real(double), dimension(:,:,:), intent(out) :: rsa
        integer :: b1

        character(1), parameter :: transa = 't'
        character(1), parameter :: transb = 'n'
        integer, parameter :: nc = 9, nr = 20
        integer :: ig, ir, is, nb, ng, ns
        real(double) :: dr, r, wt
        real(double), dimension(3) :: gk
        real(double), dimension(0:3) :: rm
        real(double), allocatable, dimension(:) :: gkm, gkr, sbr
        complex(double), parameter :: alpha = (1.0_double,0.0_double), beta = (0.0_double,0.0_double)
        complex(double), allocatable, dimension(:) :: c1
        complex(double), allocatable, dimension(:,:) :: cm2l, tmp, ylm
        complex(double), allocatable, dimension(:,:,:,:) :: cm4, cm4l
        type(multibasis_rep), pointer :: worm_mb

        call my(mv)

        worm_mb => wormhole(mv%o%mb)

        nb = size(rsa,2)
        ns = size(rsa,3)
        ng = size(worm_mb%gpt,1)

        allocate ( gkm(ng), ylm(ng,nc) )
        select case (mode)
        case ("L", "l", "LM", "lm")
          do ig = 1,ng
            gk(1) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = sqrt(gk(1)**2 + gk(2)**2 + gk(3)**2)
            ylm(ig,1) = spherical_harmonic(0, 0,gk)
            ylm(ig,2) = spherical_harmonic(1,-1,gk)
            ylm(ig,3) = spherical_harmonic(1, 0,gk)
            ylm(ig,4) = spherical_harmonic(1,+1,gk)
            ylm(ig,5) = spherical_harmonic(2,-2,gk)
            ylm(ig,6) = spherical_harmonic(2,-1,gk)
            ylm(ig,7) = spherical_harmonic(2, 0,gk)
            ylm(ig,8) = spherical_harmonic(2,+1,gk)
            ylm(ig,9) = spherical_harmonic(2,+2,gk)
          end do
        case ("XYZ", "xyz")
          do ig = 1,ng
            gk(1) = worm_mb%gpt(ig,1) + worm_mb%kpt(1)
            gk(2) = worm_mb%gpt(ig,2) + worm_mb%kpt(2)
            gk(3) = worm_mb%gpt(ig,3) + worm_mb%kpt(3)
            gkm(ig) = sqrt(gk(1)**2 + gk(2)**2 + gk(3)**2)
            ylm(ig,1) = real_spherical_harmonic(0, 0,gk)
            ylm(ig,2) = real_spherical_harmonic(1,+1,gk)
            ylm(ig,3) = real_spherical_harmonic(1,-1,gk)
            ylm(ig,4) = real_spherical_harmonic(1, 0,gk)
            ylm(ig,5) = real_spherical_harmonic(2,-2,gk)
            ylm(ig,6) = real_spherical_harmonic(2,+1,gk)
            ylm(ig,7) = real_spherical_harmonic(2,-1,gk)
            ylm(ig,8) = real_spherical_harmonic(2,+2,gk)
            ylm(ig,9) = real_spherical_harmonic(2, 0,gk)
          end do
        end select

        allocate ( c1(ng), gkr(ng), sbr(ng), tmp(ng,nc) )
        allocate ( cm2l(nc,nb), cm4(nc,nb,nr,ns), cm4l(nc,nb,nr,ns) )
        do is = 1,ns
          dr = sd(4,is)/real(nr,double)
          c1 = exp((0.0_double,1.0_double)*(worm_mb%gpt(:,1)*sd(1,is) + worm_mb%gpt(:,2)*sd(2,is) + worm_mb%gpt(:,3)*sd(3,is)))
          do ir = 1,nr
            r = dr*real(ir,double)
            gkr = gkm*r
            tmp(:,1) = ylm(:,1)*c1*spherical_bessel(gkr,0,.true.)
            sbr = spherical_bessel(gkr,1,.true.)*r
            tmp(:,2) = ylm(:,2)*c1*sbr
            tmp(:,3) = ylm(:,3)*c1*sbr
            tmp(:,4) = ylm(:,4)*c1*sbr
            sbr = spherical_bessel(gkr,2,.true.)*r**2
            tmp(:,5) = ylm(:,5)*c1*sbr
            tmp(:,6) = ylm(:,6)*c1*sbr
            tmp(:,7) = ylm(:,7)*c1*sbr
            tmp(:,8) = ylm(:,8)*c1*sbr
            tmp(:,9) = ylm(:,9)*c1*sbr
            call start_timer("multivector: zgemm")
            call zgemm(transa,transb,nc,nb,ng,alpha,tmp(1,1),ng,mv%o%mat(1,b1),ng,beta,cm2l(1,1),nc)
            call stop_timer("multivector: zgemm")
            cm4l(:,:,ir,is) = cm2l
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,cm4l,cm4)

        rm(0) = (2.0_double/45.0_double)*32.0_double*(four_pi)**2/worm_mb%cellvol    ! Bode's rule: nr must be 4 x integer
        rm(1) = (2.0_double/45.0_double)*14.0_double*(four_pi)**2/worm_mb%cellvol
        rm(2) = (2.0_double/45.0_double)*32.0_double*(four_pi)**2/worm_mb%cellvol
        rm(3) = (2.0_double/45.0_double)*12.0_double*(four_pi)**2/worm_mb%cellvol
        rsa = 0.0_double
        do is = 1,ns
          dr = sd(4,is)/real(nr,double)
          do ir = 1,nr
            r = dr*real(ir,double)
            wt = rm(mod(ir+1,4))*dr*r**2
            if (ir == nr) wt = wt/2.0_double
            rsa(:,:,is) = rsa(:,:,is) + wt*real(cm4(:,:,ir,is)*conjg(cm4(:,:,ir,is)),double)
          end do
        end do

        if (allocated( gkm )) deallocate( gkm )
        if (allocated( ylm )) deallocate( ylm )
        if (allocated( c1 )) deallocate( c1 )
        if (allocated( gkr )) deallocate( gkr )
        if (allocated( sbr )) deallocate( sbr )
        if (allocated( tmp )) deallocate( tmp )
        if (allocated( cm2l )) deallocate( cm2l )
        if (allocated( cm4 )) deallocate( cm4 )
        if (allocated( cm4l )) deallocate( cm4l )

        nullify( worm_mb )

        call glean(thy(mv))

        if (error("Exit multivector_mod::decompose_large_i")) continue

      end subroutine

      subroutine kernel_multiply_mv_r_i(ng,nb,mat,rv)
        integer :: nb, ng
        real(double), dimension(nb) :: rv
        complex(double), dimension(ng,nb) :: mat

        integer :: ib, ig
        real(double), dimension(:), allocatable :: rv_local

        allocate( rv_local(nb) )
        do ib = 1,nb
          rv_local(ib) = 0.0_double
          do ig = 1,ng
            rv_local(ib) = rv_local(ib) + real(mat(ig,ib)*conjg(mat(ig,ib)),double)
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,rv_local,rv)
        deallocate( rv_local )

      end subroutine

      subroutine kernel_multiply_2mv_r_i(ng,nb,mat1,mat2,rv)
        integer :: nb, ng
        real(double), dimension(nb) :: rv
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig
        real(double), dimension(:), allocatable :: rv_local

        allocate( rv_local(nb) )
        do ib = 1,nb
          rv_local(ib) = 0.0_double
          do ig = 1,ng
            rv_local(ib) = rv_local(ib) + real(mat1(ig,ib)*conjg(mat2(ig,ib)),double)
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,rv_local,rv)
        deallocate( rv_local )

      end subroutine

      subroutine kernel_multiply_2mv_c_i(ng,nb,mat1,mat2,cv)
        integer :: nb, ng
        complex(double), dimension(nb) :: cv
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig
        complex(double), dimension(:), allocatable :: cv_local

        allocate( cv_local(nb) )
        do ib = 1,nb
          cv_local(ib) = (0.0_double,0.0_double)
          do ig = 1,ng
            cv_local(ib) = cv_local(ib) + conjg(mat1(ig,ib))*mat2(ig,ib)
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,cv_local,cv)
        deallocate( cv_local )

      end subroutine

      subroutine kernel_trace_2mv_wt_i(ng,nb,mat1,mat2,wt,r)
        integer :: nb, ng
        real(double) :: r
        real(double), dimension(nb) :: wt
        complex(double), dimension(ng,nb) :: mat1, mat2

        integer :: ib, ig
        real(double) :: r_local

        r_local = 0.0_double
        do ib = 1,nb
          do ig = 1,ng
            r_local = r_local + 0.5_double*wt(ib)*real(mat1(ig,ib)*conjg(mat2(ig,ib)),double)
          end do
        end do
        call allreduce(KGROUP,MPI_SUM,r_local,r)

      end subroutine

      subroutine kernel_filter_mv_i(fg,lg,ng,nb,nd,alpha,beta,mat1,mat2,map,filter)
        integer :: fg, lg, nb, ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        complex(double) :: alpha, beta
        complex(double), dimension(ng,nb) :: mat1, mat2
        complex(double), dimension(nd(1),nd(2),nd(3)) :: filter

        integer :: ib, ig, igp
        complex(double), dimension(:), allocatable :: bp_filter

        allocate( bp_filter(ng) )
        igp = 0
        do ig = fg,lg
          igp = igp + 1
          bp_filter(igp) = beta*filter(map(1,ig),map(2,ig),map(3,ig))
        end do
        if (alpha == (0.0_double,0.0_double)) then
          do ib = 1,nb
            do ig = 1,ng
              mat1(ig,ib) = bp_filter(ig)*mat2(ig,ib)
            end do
          end do
        else
          do ib = 1,nb
            do ig = 1,ng
              mat1(ig,ib) = alpha*mat1(ig,ib) + bp_filter(ig)*mat2(ig,ib)
            end do
          end do
        end if
        deallocate( bp_filter )

      end subroutine

      subroutine gather_wf_i(nd,wf,ng,map,wfv)
        integer :: ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        complex(double), dimension(ng) :: wfv
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf

        integer :: ig

        call start_timer("multivector: gather_wf")
        do ig = 1,ng
          wfv(ig) = wf(map(1,ig),map(2,ig),map(3,ig))
        end do
        call stop_timer("multivector: gather_wf")

      end subroutine

      subroutine gather_wf_portion_i(nd,wf,ng,map,r,wfv)
        integer :: ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        real(double) :: r
        complex(double), dimension(ng) :: wfv
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf

        integer :: ig

        do ig = 1,ng
          wfv(ig) = r*wf(map(1,ig),map(2,ig),map(3,ig))
        end do

      end subroutine

      subroutine gather_wf_sum_portion_i(nd,wf,ng,map,r,wfv)
        integer :: ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        real(double) :: r
        complex(double), dimension(ng) :: wfv
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf

        integer :: ig

        do ig = 1,ng
          wfv(ig) = wfv(ig) + r*wf(map(1,ig),map(2,ig),map(3,ig))
        end do

      end subroutine

      subroutine scatter_wf_i(nd,wf,ng,map,wfv)
        integer :: ng
        integer, dimension(3) :: nd
        integer, dimension(3,ng) :: map
        complex(double), dimension(ng) :: wfv
        complex(double), dimension(nd(1),nd(2),nd(3)) :: wf

        integer :: i1, i2, i3, ig

        call start_timer("multivector: scatter_wf")
!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          wf(i1,i2,i3) = (0.0_double,0.0_double)
        end do
        end do
        end do
!$omp end parallel do
!$omp parallel do
        do ig = 1,ng
          wf(map(1,ig),map(2,ig),map(3,ig)) = wfv(ig)
        end do
!$omp end parallel do
        call stop_timer("multivector: scatter_wf")

      end subroutine

      subroutine extract_vector_i(ib2,n1,n2,mat,nb2,bmat)
        integer :: ib2, n1, n2, nb2
        complex(double), dimension(n1,n2) :: mat
        complex(double), dimension(n1,nb2) :: bmat

        integer :: i1, i2

        do i2 = 1,n2
          do i1 = 1,n1
            mat(i1,i2) = bmat(i1,ib2)
          end do
        end do

      end subroutine

      subroutine extract_vectors_i(fb2,lb2,n1,n2,mat,nb2,bmat)
        integer :: fb2, lb2, n1, n2, nb2
        complex(double), dimension(n1,n2) :: mat
        complex(double), dimension(n1,nb2) :: bmat

        integer :: i1, i2, ib2

        i2 = 0
        do ib2 = fb2,lb2
          i2 = i2 + 1
          do i1 = 1,n1
            mat(i1,i2) = bmat(i1,ib2)
          end do
        end do

      end subroutine

      subroutine extract_vectors_indexed_i(fi,li,ni,index,n1,n2,mat,nb2,bmat)
        integer :: fi, li, ni, n1, n2, nb2
        integer, dimension(ni) :: index
        complex(double), dimension(n1,n2) :: mat
        complex(double), dimension(n1,nb2) :: bmat

        integer :: ii, i1, i2, ib2

        i2 = 0
        do ii = fi,li
          i2 = i2 + 1
          ib2 = index(ii)
          do i1 = 1,n1
            mat(i1,i2) = bmat(i1,ib2)
          end do
        end do

      end subroutine

      subroutine augment_vector_i(fb2,lb2,tb2,n1,n2,mat,nb2,bmat)
        integer :: fb2, lb2, tb2, n1, n2, nb2
        complex(double), dimension(n1,n2) :: mat
        complex(double), dimension(n1,nb2) :: bmat

        integer :: i1, i2, ib2

        i2 = 0
        do ib2 = fb2,lb2
          i2 = i2 + 1
          do i1 = 1,n1
            bmat(i1,tb2) = bmat(i1,tb2) + mat(i1,i2)
          end do
        end do

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

      subroutine augment_vectors_ab_i(fb2,lb2,n1,n2,mat,nb2,bmat,alpha,beta)
        integer :: fb2, lb2, n1, n2, nb2
        complex(double) :: alpha, beta
        complex(double), dimension(n1,n2) :: mat
        complex(double), dimension(n1,nb2) :: bmat

        integer :: i1, i2, ib2

        i2 = 0
        do ib2 = fb2,lb2
          i2 = i2 + 1
          do i1 = 1,n1
            bmat(i1,ib2) = alpha*bmat(i1,ib2) + beta*mat(i1,i2)
          end do
        end do

      end subroutine

      subroutine augment_vectors_indexed_i(fi,li,ni,index,n1,n2,mat,nb2,bmat)
        integer :: fi, li, ni, n1, n2, nb2
        integer, dimension(ni) :: index
        complex(double), dimension(n1,n2) :: mat
        complex(double), dimension(n1,nb2) :: bmat

        integer :: ii, i1, i2, ib2

        i2 = 0
        do ii = fi,li
          i2 = i2 + 1
          ib2 = index(ii)
          do i1 = 1,n1
            bmat(i1,ib2) = bmat(i1,ib2) + mat(i1,i2)
          end do
        end do

      end subroutine

      subroutine kernel_mm_cr_i(nd,c1,r1)
        integer, dimension(3) :: nd
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1
        real(double), dimension(nd(1),nd(2),nd(3)) :: r1

        integer :: i1, i2, i3

        call start_timer("multivector: kernel_mm_cr")
!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          c1(i1,i2,i3) = c1(i1,i2,i3)*r1(i1,i2,i3)
        end do
        end do
        end do
!$omp end parallel do
        call stop_timer("multivector: kernel_mm_cr")

      end subroutine

      subroutine kernel_mm_cc_i(nd,c1,c2)
        integer, dimension(3) :: nd
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1, c2

        integer :: i1, i2, i3

!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          c1(i1,i2,i3) = c1(i1,i2,i3)*c2(i1,i2,i3)
        end do
        end do
        end do
!$omp end parallel do

      end subroutine

      subroutine kernel_mm_ccs_i(nd,c1,c2)
        integer, dimension(3) :: nd
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1, c2

        integer :: i1, i2, i3

!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          c1(i1,i2,i3) = c1(i1,i2,i3)*conjg(c2(i1,i2,i3))
        end do
        end do
        end do
!$omp end parallel do

      end subroutine

      subroutine kernel_mmm_ccc_i(nd,c1,c2,c3)
        integer, dimension(3) :: nd
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1, c2, c3

        integer :: i1, i2, i3

        call start_timer("multivector: kernel_mmm_ccc")
!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          c1(i1,i2,i3) = c2(i1,i2,i3)*conjg(c3(i1,i2,i3))
        end do
        end do
        end do
!$omp end parallel do
        call stop_timer("multivector: kernel_mmm_ccc")

      end subroutine

      subroutine kernel_rsum_mmm_ccc_i(nd,c1,c2,c3,s)
        integer, dimension(3) :: nd
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1, c2, c3
        real(double) :: s

        integer :: i1, i2, i3
        complex(double) :: cs

        cs = (0.0_double,0.0_double)
!$omp parallel do reduction(+:cs)
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          cs = cs + c1(i1,i2,i3)*conjg(c2(i1,i2,i3))*c3(i1,i2,i3)
        end do
        end do
        end do
!$omp end parallel do
        s = real(cs,double)/real(nd(1)*nd(2)*nd(3),double)

      end subroutine

      subroutine kernel_portion_mmm_ccc_rs_i(nd,c1,r,c2,c3)
        integer, dimension(3) :: nd
        real(double) :: r
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1, c2, c3

        integer :: i1, i2, i3

!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          c1(i1,i2,i3) = r*c2(i1,i2,i3)*c3(i1,i2,i3)
        end do
        end do
        end do
!$omp end parallel do

      end subroutine

      subroutine kernel_csum_portion_mmm_ccc_rs_i(nd,c1,r,c2,c3)
        integer, dimension(3) :: nd
        real(double) :: r
        complex(double), dimension(nd(1),nd(2),nd(3)) :: c1, c2, c3

        integer :: i1, i2, i3

        call start_timer("multivector: kernel_csum_portion_mmm_ccc_rs")
!$omp parallel do
        do i3 = 1,nd(3)
        do i2 = 1,nd(2)
        do i1 = 1,nd(1)
          c1(i1,i2,i3) = c1(i1,i2,i3) + r*c2(i1,i2,i3)*c3(i1,i2,i3)
        end do
        end do
        end do
!$omp end parallel do
        call stop_timer("multivector: kernel_csum_portion_mmm_ccc_rs")

      end subroutine

      end module

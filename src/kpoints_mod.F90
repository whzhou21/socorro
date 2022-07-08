!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module kpoints_mod
!doc$ module kpoints_mod

!     One datatype is available here: type(kpoints_obj).

!     kpoints_mod creates and maintains sampling points (k-points) in the Brillouin zone.

      use kind_mod
      use path_mod
      use mpi_mod
      use arg_mod
      use error_mod
      use io_mod
      use diary_mod
      use tagio_mod
      use ghost_mod
      use math_mod
      use lattice_mod
      use crystal_mod
      use external_mod
      use symmetry_mod

!cod$
      implicit none
      private

      integer, parameter :: GMP = 1   ! Generate Special Points
      integer, parameter :: GBP = 2   ! Generate Baldereschi Point
      integer, parameter :: USP = 3   ! User Supplied Points

      type :: set
        real(double), dimension(3) :: sp        ! sampling point in the lattice representation
        real(double), dimension(3) :: spx       ! sampling point in the cartesian representation
        logical, dimension(:), pointer :: map   ! logical map to symmetry-related sampling points
        integer :: deg                          ! degeneracy wrt. symmetry reduction
        real(double) :: wgt                     ! weight wrt. symmetry reduction
      end type

      type :: kpoints_rep
        integer :: ref
        type(ghost) :: g
        type(ghost) :: g_lattice                         ! lattice ghost
        type(ghost) :: g_double_group                    ! double-group ghost (used only with GMP mode)
        integer :: mode                                  ! mode for obtaining sampling points
        integer, dimension(3) :: mpp                     ! Monkhorst-Pack parameters (used only with GMP mode)
        integer, dimension(3) :: mpq                     ! Monkhorst-Pack q parameters (used only with GMP mode)
        integer, dimension(3) :: mps                     ! Monkhorst-Pack s parameters (used only with GMP mode)
        type(set), dimension(:), pointer :: set          ! sampling point and related items
      end type

      type, public :: kpoints_obj
        private
        integer :: ref
        type(kpoints_rep), pointer :: o
      end type

!doc$
      public :: kpoints
      public :: update
      public :: change
      public :: my
      public :: thy
      public :: glean
      public :: bequeath
      public :: assignment(=)
      public :: x_ref
      public :: x_ghost
      public :: x_n_kpoints
      public :: x_kpoint
      public :: x_kpoints
      public :: x_kweight
      public :: x_kweights
      public :: x_kmap
      !public :: x_kmaps  ! Issue with implicit logical to double conversion
      public :: diary
      public :: save
      public :: write_restart

!cod$
      interface kpoints
        module procedure constructor_kp, constructor_kpc
      end interface
      interface update
        module procedure update_kp
      end interface
      interface change
        module procedure change_kp
      end interface
      interface my
        module procedure my_kp, my_new_kp
      end interface
      interface thy
        module procedure thy_kp
      end interface
      interface glean
        module procedure glean_kp
      end interface
      interface bequeath
        module procedure bequeath_kp
      end interface
      interface assignment(=)
        module procedure assign_kp
      end interface
      interface x_ref
        module procedure kp_ref
      end interface
      interface x_ghost
        module procedure kp_ghost
      end interface
      interface x_n_kpoints
         module procedure kp_n_kpoints
      end interface
      interface x_kpoint
         module procedure kp_kpoint
      end interface
      interface x_kpoints
         module procedure kp_kpoints
      end interface
      interface x_kweight
         module procedure kp_kweight
      end interface
      interface x_kweights
         module procedure kp_kweights
      end interface
      interface x_kmap
         module procedure kp_kmap
      end interface
      !interface x_kmaps
      !   module procedure kp_kmaps
      !end interface
      interface diary
        module procedure diary_kp
      end interface
      interface save
        module procedure save_kp
      end interface
      interface write_restart
        module procedure write_restart_kp
      end interface

      contains

! public routines

      function constructor_kp(ext,restf) result(kp)
!doc$ function kpoints(ext,restf) result (kp)
        type(external_obj) :: ext
        type(tagio_obj), optional :: restf
        type(kpoints_obj) :: kp
!       effects: Constructs a new kp.
!       errors: Passes errors.

!cod$
        logical :: found
        character(1) :: tios
        character(line_len) :: tag
        integer :: mode
        integer(long) :: dsize, iosl, ndata, s4

        call my(ext)
        if (present(restf)) call my(restf)

        kp%ref = 0
        allocate( kp%o )
        kp%o%ref = 0
        kp%o%g = x_ghost()

        ! open the K-POINTS block
        if (present(restf)) then
          if (i_access(restf)) tios = findfirsttag(restf,"K-POINTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios /= TAG_START_BLOCK,"ERROR: K-POINTS block was not found")) goto 200
          if (i_access(restf)) call openblock(restf)
        end if

        if (present(restf)) then

          ! read the mode
          if (i_access(restf)) tios = findfirsttag(restf,"MODE")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: MODE tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long ; ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            mode = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,mode)

          ! read or generate k-point information
          select case (mode)
          case (GMP)
            call gmp_i(ext,kp%o,restf) ; if (error()) goto 100
          case (GBP)
            call gbp_i(ext,kp%o) ; if (error()) goto 100
          case (USP)
            call usp_i(ext,kp%o,restf) ; if (error()) goto 100
          end select

        else

          ! read the mode
          call arglc("kpoints",tag,found)
          if (.not.found) tag = "gmp"

          ! read or generate k-point information
          select case (trim(tag))
          case ("gmp")
            call gmp_i(ext,kp%o) ; if (error()) goto 200
          case ("gbp")
            call gbp_i(ext,kp%o) ; if (error()) goto 200
          case ("usp")
            call usp_i(ext,kp%o) ; if (error()) goto 200
          case default
            if (error(.true.,"ERROR: kpoints tag was not recognized")) goto 100
          end select

        end if

        ! close the K-POINTS block
100     if (present(restf)) then
          if (i_access(restf)) call closeblock(restf)
        end if

200     call glean(thy(ext))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit kpoints_mod::constructor_kp")) continue

      end function

      function constructor_kpc(lat,lat_grp,dbl_grp) result(kp)
!doc$ function kpoints(lat,lat_grp,dbl_grp) result (kp)
        type(lattice_obj) :: lat
        type(point_group_obj) :: lat_grp, dbl_grp
        type(kpoints_obj) :: kp
!       effects: Constructs a new kp.
!       errors: Passes errors.

!cod$
        logical :: found
        character(line_len) :: tag

        call my(lat)
        call my(lat_grp)
        call my(dbl_grp)

        kp%ref = 0
        allocate( kp%o )
        kp%o%ref = 0
        kp%o%g = x_ghost()

        call arglc("kpoints",tag,found)
        if (.not.found) tag = "gmp"
        select case (trim(tag))
        case ("gmp")
          call gmp_ck_i(lat,lat_grp,dbl_grp,kp%o) ; if (error()) goto 100
        case default
          if (error(.true.,"ERROR: kpoints tag was not recognized")) goto 100
        end select

100     call glean(thy(lat))
        call glean(thy(lat_grp))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::constructor_kpc")) continue

      end function

      subroutine update_kp(kp,ext)
!doc$ subroutine update(kp,ext)
        type(kpoints_obj) :: kp
        type(external_obj), optional :: ext
!       modifies: kp
!       effects: Updates kp.
!       errors: Passes errors.

!cod$
        logical :: remake

        call my(kp)
        if (present(ext)) call my(ext)

        remake = .false.
        if (present(ext)) remake = remake .or. change(kp,ext)

        if (remake) then
          call own_i(kp)
          kp%o%g = x_ghost()
          select case (kp%o%mode)
          case (GMP)
            call monkhorst_pack_i(ext,kp%o) ; if (error()) goto 100
          case (GBP)
            call baldereschi_point_i(ext,kp%o) ; if (error()) goto 100
          case (USP)
            call update_usp_i(ext,kp%o) ; if (error()) goto 100
          end select
        end if

100     call glean(thy(kp))
        if (present(ext)) call glean(thy(ext))

        if (error("Exit kpoints_mod::update_kp")) continue

      end subroutine

      function change_kp(kp,ext) result(ckp)
!doc$ function change(kp,ext) result(ckp)
        type(kpoints_obj) :: kp
        type(external_obj) :: ext
        logical :: ckp
!       effects: Determines whether or not kp will change due to ext.
!       errors: Passes errors.

!cod$
        call my(kp)
        call my(ext)

        ckp = .false.
        ckp = ckp .or. (kp%o%g_lattice /= x_ghost(x_lattice(x_crystal(ext))))
        select case (kp%o%mode)
          case (GMP)
            ckp = ckp .or. (kp%o%g_double_group /= x_ghost(x_double_group(ext)))
          case (GBP)
            ckp = ckp .or. .not.trivial(x_space_group(ext))
        end select

100     call glean(thy(kp))
        call glean(thy(ext))

        if (error("Exit kpoints_mod::change_kp")) continue

      end function

      subroutine my_kp(kp)
!doc$ subroutine my(kp)
        type(kpoints_obj) :: kp

!cod$
        kp%ref = kp%ref + 1
        kp%o%ref = kp%o%ref + 1
      end subroutine

      subroutine my_new_kp(kpi,kp)
!doc$ subroutine my(kpi,kp)
        type(kpoints_obj) :: kpi, kp

!cod$
        kp%ref = 1
        kp%o => kpi%o
        kp%o%ref = kp%o%ref + 1
      end subroutine

      function thy_kp(kp) result(kpo)
!doc$ function thy(kp) result(kpo)
        type(kpoints_obj) :: kp, kpo

!cod$
        kp%ref = kp%ref - 1
        kp%o%ref = kp%o%ref - 1
        kpo%ref = kp%ref
        kpo%o => kp%o
      end function

      subroutine glean_kp(kp)
!doc$ subroutine glean(kp)
        type(kpoints_obj) :: kp

!cod$
        integer :: ik
        if (kp%o%ref < 1) then
          if (associated( kp%o%set )) then
            do ik = 1,size(kp%o%set)
              if (associated( kp%o%set(ik)%map )) deallocate( kp%o%set(ik)%map )
            end do
            deallocate( kp%o%set )
          end if
          deallocate( kp%o )
        end if
      end subroutine

      subroutine bequeath_kp(kp)
!doc$ subroutine bequeath(kp)
        type(kpoints_obj) :: kp

!cod$
        continue
      end subroutine

      subroutine assign_kp(kp,kp2)
!doc$ subroutine assign(kp,kp2)
        type(kpoints_obj), intent(inout) :: kp
        type(kpoints_obj), intent(in) :: kp2
!       effects: kp = kp2

!cod$
        type(kpoints_obj) :: kpt
        call my(kp2)
        kpt%o => kp%o
        kp%o%ref = kp%o%ref - kp%ref
        kp%o => kp2%o
        kp%o%ref = kp%o%ref + kp%ref
        call glean(kpt)
        call glean(thy(kp2))
      end subroutine

      function kp_ref(kp) result(r)
!doc$ function x_ref(kp) result(r)
        type(kpoints_obj) :: kp
        integer, dimension(2) :: r
!       effects: Returns kp%ref and kp%o%ref.

!cod$
        r(1) = kp%ref
        r(2) = kp%o%ref
        call glean(kp)
      end function

      function kp_ghost(kp) result(g)
!doc$ function x_ghost(kp) result(g)
        type(kpoints_obj) :: kp
        type(ghost) :: g
!       effects: Returns ghost of kp.

!cod$
        call my(kp)
        g = kp%o%g
        call glean(thy(kp))
      end function

      function kp_n_kpoints(kp) result(n)
!doc$ function x_n_kpoints(kp) result(n)
        type(kpoints_obj) :: kp
        integer :: n
!       effects: Returns the number of k-points.

!cod$
        call my(kp)
        n = size(kp%o%set)
        call glean(thy(kp))
      end function

      function kp_kpoint(kp,ik) result(sp)
!doc$ function x_kpoint(kp,ik) result(sp)
        type(kpoints_obj) :: kp
        integer :: ik
        real(double), dimension(3) :: sp
!       effects: Returns the ik k-point in the lattice representation.
!       errors: ik out of range.

!cod$
        call my(kp)
        if (error((ik < 1) .or. (ik > size(kp%o%set)),"ERROR: index is out of range")) goto 100
        sp = kp%o%set(ik)%sp
100     call glean(thy(kp))
        if (error("Exit kpoints_mod::kp_kpoint")) continue
      end function

      function kp_kpoints(kp) result(sps)
!doc$ function x_kpoints(kp) result(sps)
        type(kpoints_obj) :: kp
        real(double), dimension(size(kp%o%set(1)%sp),size(kp%o%set)) :: sps
!       effects: Returns the k-points in the lattice representation.

!cod$
        integer :: ik
        call my(kp)
        do ik = 1,size(kp%o%set)
          sps(:,ik) = kp%o%set(ik)%sp
        end do
        call glean(thy(kp))
      end function

      function kp_kweight(kp,ik) result(wt)
!doc$ function x_kweight(kp,ik) result(wt)
        type(kpoints_obj) :: kp
        integer :: ik
        real(double) :: wt
!       effects: Returns the ik k-point weight.
!       errors: ik out of range.

!cod$
        call my(kp)
        if (error((ik < 1) .or. (ik > size(kp%o%set)),"ERROR: ik is out of range")) goto 100
        wt = kp%o%set(ik)%wgt
100     call glean(thy(kp))
        if (error("Exit kpoints_mod::kp_kweight")) continue
      end function

      function kp_kweights(kp) result(wts)
!doc$ function x_kweights(kp) result(wts)
        type(kpoints_obj) :: kp
        real(double), dimension(size(kp%o%set)) :: wts
!       effects: Returns the k-point weights.

!cod$
        integer :: ik
        call my(kp)
        do ik = 1,size(kp%o%set)
          wts(ik) = kp%o%set(ik)%wgt
        end do
        call glean(thy(kp))
      end function

      function kp_kmap(kp,ik) result(map)
!doc$ function x_kmap(kp,ik) result(map)
        type(kpoints_obj) :: kp
        integer :: ik
        logical, dimension(size(kp%o%set)) :: map
!       effects: Returns the ik k-point map.
!       errors: ik out of range.

!cod$
        call my(kp)
        if (error((ik < 1) .or. (ik > size(kp%o%set)),"ERROR: index is out of range")) goto 100
        map = kp%o%set(ik)%map
100     call glean(thy(kp))
        if (error("Exit kpoints_mod::kp_kmap")) continue
      end function

      function kp_kmaps(kp) result(maps)
!doc$ function x_kmaps(kp) result(maps)
        type(kpoints_obj) :: kp
        real(double), dimension(size(kp%o%set),size(kp%o%set)) :: maps
!       effects: Returns the k-point maps.

!cod$
        integer :: ik
        call my(kp)
        do ik = 1,size(kp%o%set)
        !  maps(ik,:) = kp%o%set(ik)%map ! CHANGED
        end do
        call glean(thy(kp))
      end function

      subroutine diary_kp(kp)
!doc$ subroutine diary(kp)
        type(kpoints_obj) :: kp
!       effects: Writes kp information to the diary.

!cod$
        integer :: ik
        call my(kp)
        if (i_access(diaryfile())) then
          write(x_unit(diaryfile()),'(/,t4,"Brillouin-zone sampling:")')
          select case (kp%o%mode)
          case (GMP)
            write(x_unit(diaryfile()),'(/,t6,"Monkhorst-Pack method:")')
            if (kp%o%mps(1) == 0) then
              write(x_unit(diaryfile()),'(/,t8,"q1: ",i0," --> ",i0)') kp%o%mpp(1), kp%o%mpq(1)
            else
              write(x_unit(diaryfile()),'(/,t8,"q1: ",i0," --> ",i0," + shift")') kp%o%mpp(1), kp%o%mpq(1)
            end if
            if (kp%o%mps(2) == 0) then
              write(x_unit(diaryfile()),'(t8,"q2: ",i0," --> ",i0)') kp%o%mpp(2), kp%o%mpq(2)
            else
              write(x_unit(diaryfile()),'(t8,"q2: ",i0," --> ",i0," + shift")') kp%o%mpp(2), kp%o%mpq(2)
            end if
            if (kp%o%mps(3) == 0) then
              write(x_unit(diaryfile()),'(t8,"q3: ",i0," --> ",i0)') kp%o%mpp(3), kp%o%mpq(3)
            else
              write(x_unit(diaryfile()),'(t8,"q3: ",i0," --> ",i0," + shift")') kp%o%mpp(3), kp%o%mpq(3)
            end if
          case (GBP)
            write(x_unit(diaryfile()),'(/,t6,"Baldereschi method:")')
          case (USP)
            write(x_unit(diaryfile()),'(/,t6,"User-supplied points:")')
          end select
          write(x_unit(diaryfile()),'(/,t8,"k-point",14x,"{b1,b2,b3}",26x,"{x,y,z}",18x,"weight",/)')
          do ik = 1,size(kp%o%set)
            write(x_unit(diaryfile()),'(t8,i5,5x,sp,3f10.6,5x,3f10.6,4x,ss,f10.6)') &
                                        ik, kp%o%set(ik)%sp, kp%o%set(ik)%spx, kp%o%set(ik)%wgt
          end do
        end if
        call glean(thy(kp))
      end subroutine

      subroutine save_kp(kp,rep)
!doc$ subroutine save(kp,rep)
        type(kpoints_obj) :: kp
        integer, intent(in), optional :: rep
!       modifies: file system
!       effects: Saves k-points and degeneracies in file new_kpoints_path. rep = 1 writes the k-points in the
!                (default) lattice representation. rep = 2 writes the k-points in the cartesian representation.
!       errors: File I/O problems.

!cod$
        integer :: ik, ios, local_rep, nk
        type(file_obj) :: f

        call my(kp)

        call my(file(trim(new_kpoints_path)),f)

        local_rep = 1
        if (present(rep)) then
          select case (rep)
          case (2)
            local_rep = 2
          end select
        end if

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100

        nk = size(kp%o%set)
        if (i_access(f)) then
          write(x_unit(f),'(i0)') nk
          select case (local_rep)
          case (1)
            write(x_unit(f),'(a)') "lattice"
            do ik = 1,nk
              write(x_unit(f),'(3f16.10,i5)') kp%o%set(ik)%sp, kp%o%set(ik)%deg
            end do
          case (2)
            write(x_unit(f),'(a)') "cartesian"
            do ik = 1,nk
              write(x_unit(f),'(3f16.10,i5)') kp%o%set(ik)%spx, kp%o%set(ik)%deg
            end do
          end select
        end if

        if (i_access(f)) close(x_unit(f))

100     call glean(thy(kp))
        call glean(thy(f))

        if (error("Exit kpoints_mod::save_kp")) continue

      end subroutine

      subroutine write_restart_kp(kp,nrestf)
!doc$ subroutine write_restart(kp,nrestf)
        type(kpoints_obj) :: kp
        type(tagio_obj) :: nrestf
!       modifies: nrestf
!       effects: Writes kp restart information to nrestf.

!cod$
        integer :: ik
        integer(long) :: dsize, iosl, ndata, s4
        integer(long), dimension(3) :: v4

        call my(kp)
        call my(nrestf)

        if (i_access(nrestf)) then

          ! start the K-POINTS block
          call startblock(nrestf,"K-POINTS")

          ! write the mode
          call writetag(nrestf,"MODE")
          s4 = kp%o%mode ; dsize = sizeof_long ; ndata = 1
          call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)

          select case (kp%o%mode)
          case (GMP)
            call writetag(nrestf,"MP_PARAMETERS")
            v4 = kp%o%mpp ; dsize = sizeof_long ; ndata = 3
            call writef(v4,dsize,ndata,x_tagfd(nrestf),iosl)
          case (USP)
            call writetag(nrestf,"NUMBER_OF_K-POINTS")
            s4 = size(kp%o%set) ; dsize = sizeof_long ; ndata = 1
            call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)
            call writetag(nrestf,"K-POINT_DATA")
            do ik = 1,size(kp%o%set)
              dsize = sizeof_double ; ndata = 3
              call writef(kp%o%set(ik)%sp(1),dsize,ndata,x_tagfd(nrestf),iosl)
              s4 = kp%o%set(ik)%deg ; dsize = sizeof_long ; ndata = 1
              call writef(s4,dsize,ndata,x_tagfd(nrestf),iosl)
            end do
          end select

          ! end the K-POINTS block
          call endblock(nrestf)

        end if

        call glean(thy(kp))
        call glean(thy(nrestf))

        if (error("Exit kpoints_mod::write_restart_kp")) continue

      end subroutine

! private routines

      subroutine gmp_i(ext,kpr,restf)
        type(external_obj) :: ext
        type(kpoints_rep) :: kpr
        type(tagio_obj), optional :: restf

        logical :: found
        character(1) :: tios
        integer(long) :: dsize, iosl, ndata
        integer(long), dimension(3) :: v4

        call my(ext)
        if (present(restf)) call my(restf)

        kpr%mode = GMP

        nullify( kpr%set )

        ! read the Monkhorst-Pack parameters
        if (present(restf)) then

          if (i_access(restf)) tios = findfirsttag(restf,"MP_PARAMETERS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: MP_PARAMETERS tag was not found")) goto 100
          if (i_access(restf)) then
            dsize = sizeof_long ; ndata = 3
            call readf(v4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            kpr%mpp = v4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,kpr%mpp)

        else

          call arg("mpparams",kpr%mpp,found)
          if (error(.not.found,"ERROR: mpparams tag was not found")) goto 100

        end if

        call monkhorst_pack_i(ext,kpr) ; if (error()) goto 100

100     call glean(thy(ext))
        if (present(restf)) call glean(thy(restf))

        if (error("Exit kpoints_mod::gmp_i")) continue

      end subroutine

      subroutine gmp_ck_i(lat,lat_grp,dbl_grp,kpr)
        type(lattice_obj) :: lat
        type(point_group_obj) :: lat_grp, dbl_grp
        type(kpoints_rep) :: kpr

        logical :: found

        call my(lat)
        call my(lat_grp)
        call my(dbl_grp)

        kpr%mode = GMP

        nullify( kpr%set )

        call arg("mpparams",kpr%mpp,found)
        if (error(.not.found,"ERROR: mpparams tag was not found")) goto 100

        call monkhorst_pack_ck_i(lat,lat_grp,dbl_grp,kpr) ; if (error()) goto 100

100     call glean(thy(lat))
        call glean(thy(lat_grp))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::gmp_ck_i")) continue

      end subroutine

      subroutine monkhorst_pack_i(ext,kpr)
        type(external_obj) :: ext
        type(kpoints_rep) :: kpr

        real(double), parameter :: tol_nbhd = 1.0e-9_double

        logical :: found, inside, lsave
        logical, dimension(:,:), pointer :: tmap
        character(line_len) :: tag
        integer :: i1, i2, i3, ib, ik, m, m1, m2, m3, nb, nk, q1, q2, q3, s1, s2, s3, sum_deg
        integer, dimension(:), pointer :: tdeg
        real(double) :: dp
        real(double), dimension(3) :: kx
        real(double), dimension(:,:), pointer :: b, tsp
        real(double), dimension(:,:), allocatable :: bx
        type(lattice_obj) :: lat
        type(point_group_obj) :: lat_grp, dbl_grp

        nullify( tmap )
        nullify( tdeg )
        nullify( b )
        nullify( tsp )

        call my(ext)
        call my(x_lattice(x_crystal(ext)),lat)
        call my(x_lattice_group(ext),lat_grp)
        call my(x_double_group(ext),dbl_grp)

        if (associated( kpr%set )) then
          do ik = 1,size(kpr%set)
            if (associated( kpr%set(ik)%map )) deallocate( kpr%set(ik)%map )
          end do
          deallocate( kpr%set )
        end if

        kpr%g_lattice = x_ghost(lat)
        kpr%g_double_group = x_ghost(dbl_grp)

        q1 = kpr%mpp(1)
        q2 = kpr%mpp(2)
        q3 = kpr%mpp(3)
        call check_monkhorst_pack(lat_grp,q1,q2,q3) ; if (error()) goto 100

        s1 = 0
        if (q1 < 0) then
          q1 = abs(q1)
          if (even(q1)) s1 = 1
        end if
        if (q1 == 0) q1 = 1
        m1 = mod(q1,2)
        s2 = 0
        if (q2 < 0) then
          q2 = abs(q2)
          if (even(q2)) s2 = 1
        end if
        if (q2 == 0) q2 = 1
        m2 = mod(q2,2)
        s3 = 0
        if (q3 < 0) then
          q3 = abs(q3)
          if (even(q3)) s3 = 1
        end if
        if (q3 == 0) q3 = 1
        m3 = mod(q3,2)

        kpr%mpq = (/q1,q2,q3/)
        kpr%mps = (/s1,s2,s3/)

        if ((q1 == 1) .and. (q2 == 1) .and. (q3 == 1)) then
          allocate( kpr%set(1) )
          kpr%set(1)%sp = (/0.0_double,0.0_double,0.0_double/)
          kpr%set(1)%spx = lat2f(lat,kpr%set(1)%sp)
          allocate( kpr%set(1)%map(1) )
          kpr%set(1)%map(1) = .true.
          kpr%set(1)%deg = 1
          kpr%set(1)%wgt = 1.0_double
          goto 100
        end if

        m = 2
        call wigner_seitz_vectors(lat,m,b)

        nb = size(b,2)
        allocate( bx(3,nb) )
        do ib = 1,nb
          bx(:,ib) = lat2f(lat,b(:,ib))
          bx(:,ib) = bx(:,ib)/dot_product(bx(:,ib),bx(:,ib))
        end do

        allocate( tsp(3,q1*q2*q3) )
        nk = 0
        do i1 = 1,2*q1,2
          do i2 = 1,2*q2,2
            do i3 = 1,2*q3,2

              nk = nk + 1
              tsp(1,nk) = real(i1-m1-s1,double)/real(2*q1,double)
              if (tsp(1,nk) > 0.5_double) tsp(1,nk) = tsp(1,nk) - 1.0_double
              tsp(2,nk) = real(i2-m2-s2,double)/real(2*q2,double)
              if (tsp(2,nk) > 0.5_double) tsp(2,nk) = tsp(2,nk) - 1.0_double
              tsp(3,nk) = real(i3-m3-s3,double)/real(2*q3,double)
              if (tsp(3,nk) > 0.5_double) tsp(3,nk) = tsp(3,nk) - 1.0_double

              inside = .false.
              do while (.not.inside)
                kx = lat2f(lat,tsp(:,nk))
                do ib = 1,nb
                  dp = dot_product(kx,bx(:,ib))
                  if ((dp < 0.5_double) .or. (dp .in. nbhd(0.5_double,tol_nbhd))) then
                    inside = .true.
                  else
                    inside = .false.
                    tsp(:,nk) = tsp(:,nk) - b(:,ib)
                    exit
                  end if
                end do
              end do

            end do
          end do
        end do

        call arglc("mp_closure",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on")
          call closure(lat_grp,lat,tsp)
        case ("off")
          continue
        case default
          call warn("WARNING: mp_closure tag was not recognized - not using closure")
        end select

        call arglc("mp_reduction",tag,found)
        if (.not.found) tag = "on"
        select case (trim(tag))
        case ("on")
          call reduction(dbl_grp,lat,tsp,tdeg) ; if (error()) goto 100
        case ("off")
          allocate( tdeg(size(tsp,2)) )
          tdeg = 1
        case default
          call warn("WARNING: mp_reduction tag was not recognized - using reduction")
          call reduction(dbl_grp,lat,tsp,tdeg) ; if (error()) goto 100
        end select

        call kpoint_map(dbl_grp,lat,tsp,tmap) ; if (error()) goto 100

        nk = size(tdeg)
        allocate( kpr%set(nk) )
        sum_deg = sum(tdeg)
        do ik = 1,nk
          kpr%set(ik)%sp = tsp(:,ik)
          kpr%set(ik)%spx = lat2f(lat,tsp(:,ik))
          allocate( kpr%set(ik)%map(nk) )
          kpr%set(ik)%map = tmap(ik,:)
          kpr%set(ik)%deg = tdeg(ik)
          kpr%set(ik)%wgt = real(tdeg(ik),double)/real(sum_deg,double)
        end do

        call arg("save_kpoints",lsave,found)
        if (.not.found) lsave = .false.
        if (lsave) call save_i(kpr)

100     if (associated( tmap )) deallocate( tmap )
        if (associated( tdeg )) deallocate( tdeg )
        if (associated( tsp )) deallocate( tsp )
        if (associated( b )) deallocate( b )
        if (allocated( bx )) deallocate( bx )

        call glean(thy(ext))

        call glean(thy(lat))
        call glean(thy(lat_grp))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::monkhorst_pack_i")) continue

      end subroutine

      subroutine monkhorst_pack_ck_i(lat,lat_grp,dbl_grp,kpr)
        type(lattice_obj) :: lat
        type(point_group_obj) :: lat_grp, dbl_grp
        type(kpoints_rep) :: kpr

        real(double), parameter :: tol_nbhd = 1.0e-9_double

        logical :: found, inside, lsave
        logical, dimension(:,:), pointer :: tmap
        character(line_len) :: tag
        integer :: i1, i2, i3, ib, ik, m, m1, m2, m3, nb, nk, q1, q2, q3, s1, s2, s3, sum_deg
        integer, dimension(:), pointer :: tdeg
        real(double) :: dp
        real(double), dimension(3) :: kx
        real(double), dimension(:,:), pointer :: b, tsp
        real(double), dimension(:,:), allocatable :: bx

        nullify( tmap )
        nullify( tdeg )
        nullify( b )
        nullify( tsp )

        call my(lat)
        call my(lat_grp)
        call my(dbl_grp)

        kpr%g_lattice = x_ghost(lat)
        kpr%g_double_group = x_ghost(dbl_grp)

        q1 = kpr%mpp(1)
        q2 = kpr%mpp(2)
        q3 = kpr%mpp(3)
        call check_monkhorst_pack(lat_grp,q1,q2,q3) ; if (error()) goto 100

        s1 = 0
        if (q1 < 0) then
          q1 = abs(q1)
          if (even(q1)) s1 = 1
        end if
        if (q1 == 0) q1 = 1
        m1 = mod(q1,2)
        s2 = 0
        if (q2 < 0) then
          q2 = abs(q2)
          if (even(q2)) s2 = 1
        end if
        if (q2 == 0) q2 = 1
        m2 = mod(q2,2)
        s3 = 0
        if (q3 < 0) then
          q3 = abs(q3)
          if (even(q3)) s3 = 1
        end if
        if (q3 == 0) q3 = 1
        m3 = mod(q3,2)

        kpr%mpq = (/q1,q2,q3/)
        kpr%mps = (/s1,s2,s3/)

        if ((q1 == 1) .and. (q2 == 1) .and. (q3 == 1)) then
          allocate( kpr%set(1) )
          kpr%set(1)%sp = (/0.0_double,0.0_double,0.0_double/)
          kpr%set(1)%spx = lat2f(lat,kpr%set(1)%sp)
          allocate( kpr%set(1)%map(1) )
          kpr%set(1)%map(1) = .true.
          kpr%set(1)%deg = 1
          kpr%set(1)%wgt = 1.0_double
          goto 100
        end if

        m = 2
        call wigner_seitz_vectors(lat,m,b)

        nb = size(b,2)
        allocate( bx(3,nb) )
        do ib = 1,nb
          bx(:,ib) = lat2f(lat,b(:,ib))
          bx(:,ib) = bx(:,ib)/dot_product(bx(:,ib),bx(:,ib))
        end do

        allocate( tsp(3,q1*q2*q3) )
        nk = 0
        do i1 = 1,2*q1,2
          do i2 = 1,2*q2,2
            do i3 = 1,2*q3,2

              nk = nk + 1
              tsp(1,nk) = real(i1-m1-s1,double)/real(2*q1,double)
              if (tsp(1,nk) > 0.5_double) tsp(1,nk) = tsp(1,nk) - 1.0_double
              tsp(2,nk) = real(i2-m2-s2,double)/real(2*q2,double)
              if (tsp(2,nk) > 0.5_double) tsp(2,nk) = tsp(2,nk) - 1.0_double
              tsp(3,nk) = real(i3-m3-s3,double)/real(2*q3,double)
              if (tsp(3,nk) > 0.5_double) tsp(3,nk) = tsp(3,nk) - 1.0_double

              inside = .false.
              do while (.not.inside)
                kx = lat2f(lat,tsp(:,nk))
                do ib = 1,nb
                  dp = dot_product(kx,bx(:,ib))
                  if ((dp < 0.5_double) .or. (dp .in. nbhd(0.5_double,tol_nbhd))) then
                    inside = .true.
                  else
                    inside = .false.
                    tsp(:,nk) = tsp(:,nk) - b(:,ib)
                    exit
                  end if
                end do
              end do

            end do
          end do
        end do

        call arglc("mp_closure",tag,found)
        if (.not.found) tag = "off"
        select case (trim(tag))
        case ("on")
          call closure(lat_grp,lat,tsp)
        case ("off")
          continue
        case default
          call warn("WARNING: mp_closure tag was not recognized - not using closure")
        end select

        call arglc("mp_reduction",tag,found)
        if (.not.found) tag = "on"
        select case (trim(tag))
        case ("on")
          call reduction(dbl_grp,lat,tsp,tdeg) ; if (error()) goto 100
        case ("off")
          allocate( tdeg(size(tsp,2)) )
          tdeg = 1
        case default
          call warn("WARNING: mp_reduction tag was not recognized - using reduction")
          call reduction(dbl_grp,lat,tsp,tdeg) ; if (error()) goto 100
        end select

        call kpoint_map(dbl_grp,lat,tsp,tmap) ; if (error()) goto 100

        nk = size(tdeg)
        allocate( kpr%set(nk) )
        sum_deg = sum(tdeg)
        do ik = 1,nk
          kpr%set(ik)%sp = tsp(:,ik)
          kpr%set(ik)%spx = lat2f(lat,tsp(:,ik))
          allocate( kpr%set(ik)%map(nk) )
          kpr%set(ik)%map = tmap(ik,:)
          kpr%set(ik)%deg = tdeg(ik)
          kpr%set(ik)%wgt = real(tdeg(ik),double)/real(sum_deg,double)
        end do

        call arg("save_kpoints",lsave,found)
        if (.not.found) lsave = .false.
        if (lsave) call save_i(kpr)

100     if (associated( tmap )) deallocate( tmap )
        if (associated( tdeg )) deallocate( tdeg )
        if (associated( b )) deallocate( b )
        if (allocated( bx )) deallocate( bx )

        call glean(thy(lat))
        call glean(thy(lat_grp))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::monkhorst_pack_ck_i")) continue

      end subroutine

      subroutine gbp_i(ext,kpr)
        type(external_obj) :: ext
        type(kpoints_rep) :: kpr

        call my(ext)

        kpr%mode = GBP

        nullify( kpr%set )

        call baldereschi_point_i(ext,kpr) ; if (error()) goto 100

100     call glean(thy(ext))

        if (error("Exit kpoints_mod::gbp_i")) continue

      end subroutine

      subroutine baldereschi_point_i(ext,kpr)
        type(external_obj) :: ext
        type(kpoints_rep) :: kpr

        logical :: found, lsave
        type(lattice_obj) :: lat
        type(point_group_obj) :: dbl_grp

        call my(ext)

        call my(x_lattice(x_crystal(ext)),lat)
        call my(x_double_group(ext),dbl_grp)

        kpr%g_lattice = x_ghost(lat)
        kpr%g_double_group = x_ghost(dbl_grp)

        if (associated( kpr%set )) deallocate( kpr%set )

        if (error(.not.cubic_cell(lat),"ERROR: GBP mode can be used only with a cubic cell")) goto 100
        if (error(.not.trivial(x_space_group(ext)),"ERROR: GBP mode can be used only with C1 symmetry")) goto 100

        allocate( kpr%set(1) )
        kpr%set(1)%sp = (/0.25_double,0.25_double,0.25_double/)
        kpr%set(1)%spx = lat2f(lat,kpr%set(1)%sp)
        allocate( kpr%set(1)%map(1) )
        kpr%set(1)%map(1) = .true.
        kpr%set(1)%deg = 1
        kpr%set(1)%wgt = 1.0_double

        call arg("save_kpoints",lsave,found)
        if (.not.found) lsave = .false.
        if (lsave) call save_i(kpr)

100     call glean(thy(ext))

        call glean(thy(lat))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::baldereschi_point_i")) continue

      end subroutine

      subroutine usp_i(ext,kpr,restf)
        type(external_obj) :: ext
        type(kpoints_rep) :: kpr
        type(tagio_obj), optional :: restf

        logical :: exist_file, found, lattice_rep, lsave
        logical, dimension(:,:), pointer :: tmap
        character(1) :: tios
        character(line_len) :: rep
        integer :: ik, ios, nk, sum_deg
        integer, dimension(:), pointer :: tdeg
        integer(long) :: dsize, iosl, ndata, s4
        real(double), dimension(3) :: ssp
        real(double), dimension(:,:), pointer :: tsp
        type(lattice_obj) :: lat
        type(point_group_obj) :: dbl_grp
        type(file_obj) :: f

        nullify( tmap )
        nullify( tdeg )
        nullify( tsp )

        call my(ext)
        if (present(restf)) call my(restf)

        call my(x_lattice(x_crystal(ext)),lat)
        call my(x_double_group(ext),dbl_grp)

        kpr%g_lattice = x_ghost(lat)
        kpr%g_double_group = x_ghost(dbl_grp)

        kpr%mode = USP

        nullify( kpr%set )

        if (present(restf)) then

          ! read the number of k-points
          if (i_access(restf)) tios = findfirsttag(restf,"NUMBER_OF_K-POINTS")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: NUMBER_OF_K-POINTS tag was not found")) goto 200
          if (i_access(restf)) then
            dsize = sizeof_long ; ndata = 1
            call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
            nk = s4
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,nk)

          ! read the k-point data
          allocate( tsp(3,nk), tdeg(nk) )
          if (i_access(restf)) tios = findfirsttag(restf,"K-POINT_DATA")
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tios)
          if (error(tios == TAG_NOT_FOUND,"ERROR: K-POINT_DATA tag was not found")) goto 200
          if (i_access(restf)) then
            do ik = 1,nk
              dsize = sizeof_double ; ndata = 3
              call readf(tsp(1,ik),dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              dsize = sizeof_long ; ndata = 1
              call readf(s4,dsize,ndata,x_tagfd(restf),x_swapbytes(restf),iosl)
              tdeg(ik) = s4
            end do
          end if
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tsp)
          if (i_comm(restf)) call broadcast(FILE_SCOPE,tdeg)

        else

          call my(file(trim(kpoints_path)),f)

          ! find the k-points file
          if (i_access(f)) inquire(file=x_name(f),exist=exist_file)
          if (i_comm(f)) call broadcast(FILE_SCOPE,exist_file)
          if (error(.not.exist_file,"ERROR: k-points file was not found")) goto 200

          ! open the k-points file
          if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='old',iostat=ios)
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to open file")) goto 200

          ! read the number of k-points
          if (i_access(f)) read(x_unit(f),*,iostat=ios) nk
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to read the number of k-points")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,nk)
          allocate( kpr%set(nk) )

          ! read the k-points representation
          if (i_access(f)) read(x_unit(f),'(a)',iostat=ios) rep
          if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
          if (error(ios /= 0,"ERROR: unable to read the k-points representation")) goto 100
          if (i_comm(f)) call broadcast(FILE_SCOPE,rep)
          select case (trim(rep))
          case ("lattice","Lattice","LATTICE","lat","Lat","LAT","l","L")
            lattice_rep = .true.
          case ("cartesian","Cartesian","CARTESIAN","cart","Cart","CART","c","C")
            lattice_rep = .false.
          case default
            if (error(.true.,"ERROR: k-points representation was not recognized")) goto 100
          end select

          ! read the k-points data
          allocate( tsp(3,nk), tdeg(nk) )
          if (i_access(f)) then
            do ik = 1,nk
              read(x_unit(f),*) tsp(:,ik), tdeg(ik)
            end do
          end if
          if (i_comm(f)) call broadcast(FILE_SCOPE,tsp)
          if (i_comm(f)) call broadcast(FILE_SCOPE,tdeg)

          ! convert to the lattice representation
          if (lattice_rep) then
            continue
          else
            do ik = 1,nk
              ssp = tsp(:,ik)
              tsp(:,ik) = f2lat(lat,ssp)
            end do
          end if

100       if (i_access(f)) close(x_unit(f))

          call glean(thy(f))
          if (error()) goto 200

        end if

        call kpoint_map(dbl_grp,lat,tsp,tmap) ; if (error()) goto 200

        sum_deg = sum(tdeg)
        allocate( kpr%set(nk) )
        do ik = 1,nk
          kpr%set(ik)%sp = tsp(:,ik)
          kpr%set(ik)%spx = lat2f(lat,tsp(:,ik))
          allocate( kpr%set(ik)%map(nk) )
          kpr%set(ik)%map = tmap(ik,:)
          kpr%set(ik)%deg = tdeg(ik)
          kpr%set(ik)%wgt = real(tdeg(ik),double)/real(sum_deg,double)
        end do

        call arg("save_kpoints",lsave,found)
        if (.not.found) lsave = .false.
        if (lsave) call save_i(kpr)

200     if (associated( tmap )) deallocate( tmap )
        if (associated( tdeg )) deallocate( tdeg )
        if (associated( tsp )) deallocate( tsp )

        call glean(thy(ext))
        if (present(restf)) call glean(thy(restf))

        call glean(thy(lat))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::usp_i")) continue

      end subroutine

      subroutine update_usp_i(ext,kpr)
        type(external_obj) :: ext
        type(kpoints_rep) :: kpr

        logical, dimension(:,:), pointer :: tmap
        integer :: ik, nk
        real(double), dimension(:,:), pointer :: tsp
        type(lattice_obj) :: lat
        type(point_group_obj) :: dbl_grp

        nullify( tmap )
        nullify( tsp )

        call my(ext)

        call my(x_lattice(x_crystal(ext)),lat)
        call my(x_double_group(ext),dbl_grp)

        kpr%g_lattice = x_ghost(lat)
        kpr%g_double_group = x_ghost(dbl_grp)

        nk = size(kpr%set)
        allocate( tsp(3,nk) )
        do ik = 1,nk
          tsp(:,ik) = kpr%set(ik)%sp
        end do

        call kpoint_map(dbl_grp,lat,tsp,tmap) ; if (error()) goto 100

        nk = size(kpr%set)
        do ik = 1,nk
          kpr%set(ik)%spx = lat2f(lat,kpr%set(ik)%sp)
          kpr%set(ik)%map = tmap(ik,:)
        end do

100     if (associated( tmap )) deallocate( tmap )
        if (associated( tsp )) deallocate( tsp )

        call glean(thy(ext))

        call glean(thy(lat))
        call glean(thy(dbl_grp))

        if (error("Exit kpoints_mod::update_usp_i")) continue

      end subroutine

      subroutine save_i(kpr)
        type(kpoints_rep) :: kpr

        integer :: ik, ios, nk
        type(file_obj) :: f

        call my(file(trim(new_kpoints_path)),f)

        if (i_access(f)) open(unit=x_unit(f),file=x_name(f),status='unknown',iostat=ios)
        if (i_comm(f)) call broadcast(FILE_SCOPE,ios)
        if (error(ios /= 0,"ERROR: unable to open file")) goto 100

        nk = size(kpr%set)
        if (i_access(f)) then
          write(x_unit(f),'(i0)') nk
          write(x_unit(f),'(a)') "lattice"
          do ik = 1,nk
            write(x_unit(f),'(3f16.10,i5)') kpr%set(ik)%sp, kpr%set(ik)%deg
          end do
        end if

        if (i_access(f)) close(x_unit(f))

100     call glean(thy(f))

        if (error("Exit kpoints_mod::save_i")) continue

      end subroutine

      subroutine own_i(kp)
        type(kpoints_obj) :: kp
        type(kpoints_obj) :: kpt
        integer :: ik
        if (kp%ref < kp%o%ref) then
          allocate( kpt%o )
          kpt%o%ref = 0
          kpt%o%g = kp%o%g
          kpt%o%g_lattice = kp%o%g_lattice
          kpt%o%mode = kp%o%mode
          select case (kp%o%mode)
          case (GMP)
            kpt%o%mpp = kp%o%mpp
            kpt%o%mpq = kp%o%mpq
            kpt%o%mps = kp%o%mps
            kpt%o%g_double_group = kp%o%g_double_group
          end select
          allocate( kpt%o%set(size(kp%o%set)) )
          do ik = 1,size(kp%o%set)
            kpt%o%set(ik)%sp = kp%o%set(ik)%sp
            kpt%o%set(ik)%spx = kp%o%set(ik)%spx
            allocate( kpt%o%set(ik)%map(size(kp%o%set)) )
            kpt%o%set(ik)%map = kp%o%set(ik)%map
            kpt%o%set(ik)%deg = kp%o%set(ik)%deg
            kpt%o%set(ik)%wgt = kp%o%set(ik)%wgt
          end do
          kp%o%ref = kp%o%ref - kp%ref
          kp%o => kpt%o
          kp%o%ref = kp%o%ref + kp%ref
        end if
      end subroutine

      end module

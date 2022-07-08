!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module ewald_mod
!doc$ module ewald_mod

!     ewald_mod encapsulates procedures involving Ewald sums.

      use kind_mod
      use mpi_mod
      use error_mod
      use diary_mod
      use math_mod
      use crystal_mod
      use lattice_mod
      use atoms_mod

!cod$
      implicit none
      private

      integer, parameter :: MOD_SCOPE = CONFIG

!doc$
      public :: ewald_energy
      public :: ewald_forces
      public :: ewald_pressure
      public :: ewald_stress_tensor

!cod$
      contains

      subroutine ewald_energy(cr,q,eii)
!doc$ subroutine ewald_energy(cr,q,eii)
        type(crystal_obj) :: cr
        real(double), dimension(:), intent(in) :: q
        real(double), intent(out) :: eii
!       effects: Returns the energy due to ion-ion interactions in the presence of a uniform compensating electron density.
!       errors: If the calculation does not converge.

!cod$
        integer, parameter :: maxit = 5
        real(double), parameter :: tol_eii = 1.0e-9_double

        integer :: at1, at2, i1, i2, i3, ia, ip, it, na
        integer :: start, stop, length
        integer, dimension(3) :: maxg, maxr, ng, nr
        real(double) :: cell_volume, eta, eta2, p2, pda, pm, sum_gg, sum_rr
        real(double) :: eii_start, eii_local, eii_old
        real(double), dimension(3) :: amag, da, p
        real(double), dimension(:,:), allocatable :: pos
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(cr)

        call my(x_lattice(cr),lat)
        call my(x_atoms(cr),ats)

        amag(1) = norm(lat2r(lat,real((/1,0,0/),double)))
        amag(2) = norm(lat2r(lat,real((/0,1,0/),double)))
        amag(3) = norm(lat2r(lat,real((/0,0,1/),double)))
        cell_volume = x_cell_volume(lat)

        na = x_n_atoms(ats)
        allocate( pos(3,na) )
        do ia = 1,na
          pos(:,ia) = lat2r(lat,x_position(ats,ia))
        end do

        eta = sqrt(pi)/( product(amag)**(1.0_double/3.0_double) )
        eta2 = eta**2

        maxr = max(int(5.0_double/(eta*amag)),(/1,1,1/))
        maxg = max(int(5.0_double*eta*amag/pi),(/1,1,1/))

        call subdivide(mpi_myproc(MOD_SCOPE),mpi_nprocs(MOD_SCOPE),1,na**2,start,stop,length)

        eii_start = -2.0_double*sum(q**2)*(eta/sqrt(pi))/mpi_nprocs(MOD_SCOPE)
        do it = 1,maxit

          nr = it*maxr - (/1,1,1/)
          ng = it*maxg - (/1,1,1/)

          eii_local = eii_start
          do ip = start,stop
            at1 = int((ip-1)/na) + 1
            at2 = mod((ip-1),na) + 1
            da = pos(:,at1) - pos(:,at2)
            sum_rr = 0.0_double
            do i1 = -nr(1),+nr(1)
            do i2 = -nr(2),+nr(2)
            do i3 = -nr(3),+nr(3)
              if (all((/i1,i2,i3/) == 0)) then
                if (at1 /= at2) then
                  pm = norm(da)
                  sum_rr = sum_rr + complementary_error(eta*pm)/pm
                end if
              else
                p = lat2r(lat,real((/i1,i2,i3/),double)) + da
                pm = norm(p)
                sum_rr = sum_rr + complementary_error(eta*pm)/pm
              end if
            end do
            end do
            end do
            sum_gg = 0.0_double
            do i1 = -ng(1),+ng(1)
            do i2 = -ng(2),+ng(2)
            do i3 = -ng(3),+ng(3)
              if (all((/i1,i2,i3/) == 0)) cycle
              p = lat2f(lat,real((/i1,i2,i3/),double))
              p2 = dot_product(p,p)
              pda = dot_product(p,da)
              sum_gg = sum_gg + cos(pda)*exp(-p2/(4.0_double*eta2))/p2
            end do
            end do
            end do
            sum_gg = (four_pi/cell_volume)*sum_gg
            eii_local = eii_local + q(at1)*q(at2)*(sum_rr + sum_gg - pi/(cell_volume*eta2))
          end do
          call allreduce(MOD_SCOPE,MPI_SUM,eii_local,eii)

          if (it > 1) then
            if (eii .in. nbhd(eii_old,tol_eii)) exit
          end if
          eii_old = eii
          if (error(it == maxit,"ERROR: Ewald energy did not converge")) goto 100

        end do

100     if (allocated( pos)) deallocate( pos )

        call glean(thy(lat))
        call glean(thy(ats))

        call glean(thy(cr))

        if (error("Exit ewald_mod::ewald_energy")) continue

      end subroutine

      subroutine ewald_forces(cr,q,fii)
!doc$ subroutine ewald_forces(cr,q,fii)
        type(crystal_obj) :: cr
        real(double), dimension(:), intent(in) :: q
        real(double), dimension(:,:), intent(out) :: fii
!       modifies: fii
!       requires: fii have dimensions (3,x_n_atoms(x_atoms(cr))).
!       effects: Returns the forces due to ion-ion interactions in the presence of a uniform compensating electron density.
!       errors: If the calculation does not converge.

!cod$
        integer, parameter :: maxit = 5
        real(double), parameter :: tol_fii = 1.0e-6_double

        integer :: at1, at2, i1, i2, i3, ia, ip, it, na
        integer :: start, stop, length
        integer, dimension(3) :: maxg, maxr, ng, nr
        real(double) :: cell_volume, dif_max, eta, eta2, p2, pda, pm, r1
        real(double), dimension(3) :: amag, da, p, sum_gg, sum_rr
        real(double), dimension(:,:), allocatable :: pos, fii_local, fii_old
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(cr)

        call my(x_lattice(cr),lat)
        call my(x_atoms(cr),ats)

        amag(1) = norm(lat2r(lat,real((/1,0,0/),double)))
        amag(2) = norm(lat2r(lat,real((/0,1,0/),double)))
        amag(3) = norm(lat2r(lat,real((/0,0,1/),double)))
        cell_volume = x_cell_volume(lat)

        na = x_n_atoms(ats)
        allocate( pos(3,na) )
        do ia = 1,na
          pos(:,ia) = lat2r(lat,x_position(ats,ia))
        end do

        eta = sqrt(pi)/( product(amag)**(1.0_double/3.0_double) )
        eta2 = eta**2

        maxr = max(int(5.0_double/(eta*amag)),(/1,1,1/))
        maxg = max(int(5.0_double*eta*amag/pi),(/1,1,1/))

        call subdivide(mpi_myproc(MOD_SCOPE),mpi_nprocs(MOD_SCOPE),1,na**2,start,stop,length)

        allocate( fii_local(3,na), fii_old(3,na) )

        do it = 1,maxit

          nr = it*maxr - (/1,1,1/)
          ng = it*maxg - (/1,1,1/)

          fii_local = 0.0_double
          do ip = start,stop
            at1 = int((ip-1)/na) + 1
            at2 = mod((ip-1),na) + 1
            if (at1 == at2) cycle
            da = pos(:,at1) - pos(:,at2)
            sum_rr = 0.0_double
            do i1 = -nr(1),+nr(1)
            do i2 = -nr(2),+nr(2)
            do i3 = -nr(3),+nr(3)
              p = lat2r(lat,real((/i1,i2,i3/),double)) + da
              pm = norm(p)
              r1 = 2.0_double*eta*exp(-(eta*pm)**2)/(sqrt(pi)*pm**2) + complementary_error(eta*pm)/pm**3
              sum_rr = sum_rr + r1*p
            end do
            end do
            end do
            sum_gg = 0.0_double
            do i1 = -ng(1),+ng(1)
            do i2 = -ng(2),+ng(2)
            do i3 = -ng(3),+ng(3)
              if ( all((/i1,i2,i3/) == 0)) cycle
              p = lat2f(lat,real((/i1,i2,i3/),double))
              p2 = dot_product(p,p)
              pda = dot_product(p,da)
              r1 = sin(pda)*exp(-p2/(4.0_double*eta2))/p2
              sum_gg = sum_gg + r1*p
            end do
            end do
            end do
            sum_gg = (four_pi/cell_volume)*sum_gg
            fii_local(:,at1) = fii_local(:,at1) + 2.0_double*q(at1)*q(at2)*(sum_rr + sum_gg)
          end do
          call allreduce(MOD_SCOPE,MPI_SUM,fii_local,fii)

          if (it > 1) then
            dif_max = maxval(abs(fii - fii_old))
            if (dif_max .in. nbhd(0.0_double,tol_fii)) exit
          end if
          fii_old = fii
          if (error(it == maxit,"ERROR: Ewald forces did not converge")) goto 100

        end do

100     if (allocated( pos )) deallocate( pos )
        if (allocated( fii_local )) deallocate( fii_local )
        if (allocated( fii_old )) deallocate( fii_old )

        call glean(thy(lat))
        call glean(thy(ats))

        call glean(thy(cr))

        if (error("Exit ewald_mod::ewald_forces")) continue

      end subroutine

      subroutine ewald_pressure(cr,q,pii)
!doc$ subroutine ewald_pressure(cr,q,pii)
        type(crystal_obj) :: cr
        real(double), dimension(:) , intent(in) :: q
        real(double), intent(out) :: pii
!       effects: Returns the pressure due to ion-ion interactions in the presence of a uniform compensating electron density.
!       errors: If the calculation does not converge.

!cod$
        integer, parameter :: maxit = 5
        real(double), parameter :: tol_pii = 1.0e-6_double

        integer :: at1, at2, i1, i2, i3, ia, ip, it, na
        integer :: start, stop, length
        integer, dimension(3) :: maxg, maxr, ng, nr
        real(double) :: cell_volume, eta, p2, pm, r1, r2, r3
        real(double) :: pii_local, pii_old
        real(double), dimension(3) :: amag, da, p
        real(double), dimension(:,:), allocatable :: pos
        complex(double) :: c1
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(cr)

        call my(x_lattice(cr),lat)
        call my(x_atoms(cr),ats)

        amag(1) = norm(lat2r(lat,real((/1,0,0/),double)))
        amag(2) = norm(lat2r(lat,real((/0,1,0/),double)))
        amag(3) = norm(lat2r(lat,real((/0,0,1/),double)))
        cell_volume = x_cell_volume(lat)

        eta = sqrt(pi)/( product(amag)**(1.0_double/3.0_double) )

        na = x_n_atoms(ats)
        allocate( pos(3,na) )
        do ia = 1,na
          pos(:,ia) = lat2r(lat,x_position(ats,ia))
        end do

        maxr = max(int(5.0_double/(eta*amag)),(/1,1,1/))
        maxg = max(int(5.0_double*eta*amag/pi),(/1,1,1/))

        call subdivide(mpi_myproc(MOD_SCOPE),mpi_nprocs(MOD_SCOPE),1,na**2,start,stop,length)

        r3 = (pi*sum(q)**2)/(cell_volume*eta)**2

        do it = 1,maxit

          nr = it*maxr - (/1,1,1/)
          ng = it*maxg - (/1,1,1/)

          pii_local = 0.0_double
          do ip = start,stop
            at1 = int((ip-1)/na) + 1
            at2 = mod((ip-1),na) + 1
            da = pos(:,at1) - pos(:,at2)
            do i1 = -nr(1),+nr(1)
            do i2 = -nr(2),+nr(2)
            do i3 = -nr(3),+nr(3)
              if ((at1 == at2) .and. (all((/i1,i2,i3/) == 0))) cycle
              p = lat2r(lat,real((/i1,i2,i3/),double)) + da
              p2 = dot_product(p,p)
              pm = norm(p)
              r1 = q(at1)*q(at2)*((2.0_double*eta/sqrt(pi))*exp(-(pm*eta)**2) + complementary_error(eta*pm)/pm)/(p2*cell_volume)
              pii_local = pii_local + r1*p2
            end do
            end do
            end do
          end do
          call allreduce(MOD_SCOPE,MPI_SUM,pii_local,pii)

          do i1 = -ng(1),+ng(1)
          do i2 = -ng(2),+ng(2)
          do i3 = -ng(3),+ng(3)
            if (all((/i1,i2,i3/) == 0)) cycle
            p = lat2f(lat,real((/i1,i2,i3/),double))
            p2 = dot_product(p,p)
            c1 = cmplx(0,0,double)
            do ia = 1,na
              c1 = c1 + q(ia)*exp(cmplx(0,1,double)*(p(1)*pos(1,ia) + p(2)*pos(2,ia) + p(3)*pos(3,ia)))
            end do
            r1 = (four_pi/cell_volume)*(exp(-p2/(4.0_double*eta**2))/p2)*(c1*conjg(c1))/cell_volume
            r2 = r1*(2.0_double/p2)*(p2/(4.0_double*eta**2) + 1.0_double)
            pii = pii + 3.0_double*r1 - r2*p2
          end do
          end do
          end do

          pii = pii/3.0_double - r3

          if (it > 1) then
            if (pii .in. nbhd(pii_old,tol_pii)) exit
          end if
          pii_old = pii
          if (error(it == maxit,"ERROR: Ewald pressure did not converge")) goto 100

        end do

100     if (allocated( pos )) deallocate( pos )

        call glean(thy(lat))
        call glean(thy(ats))

        call glean(thy(cr))

        if (error("Exit ewald_mod::ewald_pressure")) continue

      end subroutine

      subroutine ewald_stress_tensor(cr,q,sii)
!doc$ subroutine ewald_stress_tensor(cr,q,sii)
        type(crystal_obj) :: cr
        real(double), dimension(:) , intent(in) :: q
        real(double), dimension(3,3), intent(out) :: sii
!       requires: sii be dimension(3,3).
!       effects: Returns the stress tensor due to ion-ion interactions in the presence of a uniform compensating electron density.
!       errors: If the calculation does not converge.

!cod$
        integer, parameter :: maxit = 5
        real(double), parameter :: tol_sii = 1.0e-6_double

        integer :: at1, at2, i1, i2, i3, ia, ip, it, na
        integer :: start, stop, length
        integer, dimension(3) :: maxg, maxr, ng, nr
        real(double) :: cell_volume, dif_max, eta, p2, pm, r1, r2, r3
        real(double), dimension(3) :: amag, da, p
        real(double), dimension(:,:), allocatable :: pos, sii_old, sii_local
        complex(double) :: c1
        type(lattice_obj) :: lat
        type(atoms_obj) :: ats

        call my(cr)

        call my(x_lattice(cr),lat)
        call my(x_atoms(cr),ats)

        amag(1) = norm(lat2r(lat,real((/1,0,0/),double)))
        amag(2) = norm(lat2r(lat,real((/0,1,0/),double)))
        amag(3) = norm(lat2r(lat,real((/0,0,1/),double)))
        cell_volume = x_cell_volume(lat)

        eta = sqrt(pi)/( product(amag)**(1.0_double/3.0_double) )

        na = x_n_atoms(ats)
        allocate( pos(3,na) )
        do ia = 1,na
          pos(:,ia) = lat2r(lat,x_position(ats,ia))
        end do

        maxr = max(int(5.0_double/(eta*amag)),(/1,1,1/))
        maxg = max(int(5.0_double*eta*amag/pi),(/1,1,1/))

        call subdivide(mpi_myproc(MOD_SCOPE),mpi_nprocs(MOD_SCOPE),1,na**2,start,stop,length)

        allocate( sii_local(3,3), sii_old(3,3) )

        r3 = (pi*sum(q)**2)/(cell_volume*eta)**2

        do it = 1,maxit

          nr = it*maxr - (/1,1,1/)
          ng = it*maxg - (/1,1,1/)

          sii_local = 0.0_double
          do ip = start,stop
            at1 = int((ip-1)/na) + 1
            at2 = mod((ip-1),na) + 1
            da = pos(:,at1) - pos(:,at2)
            do i1 = -nr(1),+nr(1)
            do i2 = -nr(2),+nr(2)
            do i3 = -nr(3),+nr(3)
              if ((at1 == at2) .and. (all((/i1,i2,i3/) == 0))) cycle
              p = lat2r(lat,real((/i1,i2,i3/),double)) + da
              p2 = dot_product(p,p)
              pm = norm(p)
              r1 = q(at1)*q(at2)*((2.0_double*eta/sqrt(pi))*exp(-(pm*eta)**2) + complementary_error(eta*pm)/pm)/(p2*cell_volume)
              sii_local(1,1) = sii_local(1,1) - r1*p(1)*p(1)
              sii_local(1,2) = sii_local(1,2) - r1*p(1)*p(2)
              sii_local(1,3) = sii_local(1,3) - r1*p(1)*p(3)
              sii_local(2,2) = sii_local(2,2) - r1*p(2)*p(2)
              sii_local(2,3) = sii_local(2,3) - r1*p(2)*p(3)
              sii_local(3,3) = sii_local(3,3) - r1*p(3)*p(3)
            end do
            end do
            end do
          end do
          call allreduce(MOD_SCOPE,MPI_SUM,sii_local,sii)

          do i1 = -ng(1),+ng(1)
          do i2 = -ng(2),+ng(2)
          do i3 = -ng(3),+ng(3)
            if (all((/i1,i2,i3/) == 0)) cycle
            p = lat2f(lat,real((/i1,i2,i3/),double))
            p2 = dot_product(p,p)
            c1 = cmplx(0,0,double)
            do ia = 1,na
              c1 = c1 + q(ia)*exp(cmplx(0,1,double)*(p(1)*pos(1,ia) + p(2)*pos(2,ia) + p(3)*pos(3,ia)))
            end do
            r1 = (four_pi/cell_volume)*(exp(-p2/(4.0_double*eta**2))/p2)*(c1*conjg(c1))/cell_volume
            r2 = r1*(2.0_double/p2)*(p2/(4.0_double*eta**2) + 1.0_double)
            sii(1,1) = sii(1,1) + r2*p(1)*p(1) - r1
            sii(1,2) = sii(1,2) + r2*p(1)*p(2)
            sii(1,3) = sii(1,3) + r2*p(1)*p(3)
            sii(2,2) = sii(2,2) + r2*p(2)*p(2) - r1
            sii(2,3) = sii(2,3) + r2*p(2)*p(3)
            sii(3,3) = sii(3,3) + r2*p(3)*p(3) - r1
          end do
          end do
          end do

          sii(1,1) = sii(1,1) + r3
          sii(2,2) = sii(2,2) + r3
          sii(3,3) = sii(3,3) + r3

          sii(2,1) = sii(1,2)
          sii(3,1) = sii(1,3)
          sii(3,2) = sii(2,3)

          if (it > 1) then
            dif_max = maxval(abs(sii - sii_old))
            if (dif_max .in. nbhd(0.0_double,tol_sii)) exit
          end if
          sii_old = sii
          if (error(it == maxit,"ERROR: Ewald stress tensor did not converge")) goto 100

        end do

100     if (allocated( pos )) deallocate( pos )
        if (allocated( sii_local )) deallocate( sii_local )
        if (allocated( sii_old )) deallocate( sii_old )

        call glean(thy(lat))
        call glean(thy(ats))

        call glean(thy(cr))

        if (error("Exit ewald_mod::ewald_stress_tensor")) continue

      end subroutine

      end module

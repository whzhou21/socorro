!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module math_mod
!doc$ module math_mod

!     Math_mod is a repository for numerical constants and mathematics procedures.

      use kind_mod
      use error_mod
      use diary_mod
      use timing_mod

!cod$
      implicit none
      private

      real(double), parameter :: pi       =            3.1415926535897932_double
      real(double), parameter :: two_pi   = 2.0_double*3.1415926535897932_double
      real(double), parameter :: four_pi  = 4.0_double*3.1415926535897932_double
      real(double), parameter :: eight_pi = 8.0_double*3.1415926535897932_double

      real(double) :: machine_precision
      real(double) :: machine_zero
      real(double) :: machine_infinity
      real(double) :: minlog
      real(double) :: maxlog
      real(double) :: minexp
      real(double) :: maxexp
      real(double) :: minlogarg
      real(double) :: maxlogarg
      real(double) :: minexparg
      real(double) :: maxexparg

      type, public :: real_nbhd
        private
        real(double) :: tol
        real(double) :: point
      end type

      type, public :: real_vector_nbhd
        private
        real(double) :: tol
        real(double), dimension(3) :: point
      end type

      type, public :: real_vector_zone
        private
        real(double) :: tol
        real(double), dimension(3) :: point
      end type

      type, public :: real_matrix_nbhd
        private
        real(double) :: tol
        real(double), dimension(3,3) :: point
      end type

      type, public :: complex_nbhd
        private
        real(double) :: tol
        complex(double) :: point
      end type

      type, public :: chebyshev_state
         integer :: phase
         integer :: period
         real(double) :: alpha
         real(double) :: beta
         real(double) :: tau0
         real(double) :: tau1
      end type 

!doc$
      public :: pi
      public :: two_pi
      public :: four_pi
      public :: eight_pi
      public :: machine_precision
      public :: machine_zero
      public :: machine_infinity
      public :: minlog
      public :: maxlog
      public :: minexp
      public :: maxexp
      public :: minlogarg
      public :: maxlogarg
      public :: minexparg
      public :: maxexparg
      public :: subdivide
      public :: factor
      public :: map3d
      public :: imap3d
      public :: nbhd
      public :: zone
      public :: operator(.in.)
      public :: operator(.out.)
      public :: even
      public :: odd
      public :: norm
      public :: cross_product
      public :: determinant
      public :: inverse
      public :: trace
      public :: solve
      public :: complementary_error
      public :: spherical_bessel
      public :: fermi_distr
      public :: fermi_distr_divdiff
      public :: fermi_entropy
      public :: random
      public :: inverse_cholesky
      public :: permutations
      public :: factorial
      public :: double_factorial
      public :: spharm
      public :: gradylm
      public :: spherical_harmonic
      public :: grad_spherical_harmonic
      public :: real_spherical_harmonic
      public :: grad_real_spherical_harmonic
      public :: simpson_integral
      public :: gauss_legendre
      public :: atom_dot_product
      public :: init_machine_constants
      public :: chebyshev_initialize
      public :: chebyshev_accelerate
      public :: besjn

!cod$
      interface subdivide
        module procedure subdivide_mine, subdivide_prime
      end interface
      interface nbhd
        module procedure nbhd_real, nbhd_complex, nbhd_vector, nbhd_matrix
      end interface
      interface zone
        module procedure zone_vector
      end interface
      interface operator(.in.)
        module procedure in_real, in_complex, in_vector, in_matrix, in_vector_zone
      end interface
      interface operator(.out.)
        module procedure out_real, out_complex, out_vector, out_matrix, out_vector_zone
      end interface
      interface norm
        module procedure norm_vector
      end interface
      interface cross_product
        module procedure cross_product_vector
      end interface
      interface determinant
        module procedure determinant_matrix
      end interface
      interface inverse
        module procedure inverse_matrix
      end interface
      interface trace
        module procedure trace_real_matrix, trace_complex_matrix
      end interface
      interface solve
        module procedure solve_r
      end interface
      interface random
        module procedure random_sjp
      end interface
      interface simpson_integral
        module procedure simpson_integral_real, simpson_integral_complex
      end interface
      interface atom_dot_product
        module procedure all_atom_dot_product_1, all_atom_dot_product_2
      end interface
      interface init_machine_constants
        module procedure init_constants
      end interface
      interface besjn
        module procedure bessj
      end interface

      contains

      subroutine subdivide_mine(myrank,mygroup,start,stop,mystart,mystop,mylength)
!doc$ subroutine subdivide(myrank,mygroup,start,stop,mystart,mystop,mylength)
        integer, intent(in) :: myrank, mygroup, start, stop
        integer, intent(out) :: mystart, mystop, mylength
!       modifies: mystart, mystop, mylength
!       effects: Takes the interval of integers [start,stop], and divides it into ngroup disjoint intervals.
!                If 0 <= myrank < ngroup, mystart, mystop, and mylength (== mystop - mystart + 1) are the
!                (myrank+1)^th interval characteristics. If myrank < 0, mystart = start and mystop = start-1.
!                If myrank >= ngroup, mystart=stop and mystop=stop-1
!doc$
        integer :: total, sl, bsl, sz
        total = stop - start + 1
        call subdivide_prime(mygroup,total,sz,sl,bsl)
        if (myrank < 0) then
          mylength = 0
          mystart = start
          mystop = start - 1
          return
        end if
        if (myrank >= mygroup) then
          mylength = 0
          mystart = stop
          mystop = stop - 1
          return
        end if
        if (myrank < bsl) then
          mylength = sz + 1
          mystart = start + myrank*(sz + 1)
          mystop = mystart + sz
        else
          mylength = sz
          mystart = start + myrank*sz + bsl
          mystop = mystart + sz - 1
        end if
      end subroutine

      subroutine subdivide_prime(group,total,size,slices,bigslices)
!doc$ subroutine subdivide(group,total,size,slices,bigslices)
        integer, intent(in) :: group, total
        integer, intent(out) :: size, slices, bigslices
!       modifies: size, slices, bigslices
!       effects: Divides total into group slices, such that size*slices + (size+1)*bigslices = total,
!                and slices+bigslices = group
!
!cod$
        bigslices = mod(total,group)
        slices = group - bigslices
        size = total/group
      end subroutine

      subroutine factor(n,p,q)
!doc$ subroutine factor(n,p,q)
        integer, intent(in) :: n
        integer, intent(inout) :: p, q
!       modifies: p and q
!       effects: Finds the factorization p*q = n such that p+q is minimized.
!
!cod$
        integer :: bestper, i
        p = 1
        q = n
        bestper = p + q
        do i = 1,n/2
          if (mod(n,i) == 0) then
            if ((i + n/i) < bestper) then
              bestper = i + n/i
              p = i
              q = n/i
            end if
          end if
        end do
      end subroutine

      function map3d(d,p) result(f)
!doc$ function map3d(d,p) result(f)
        integer, dimension(3), intent(in) :: d
        integer, intent(in) :: p
        integer, dimension(3) :: f
!       requires: 1 <= p <= product(d)
!       effects: Returns the (p)^th column ordered coordinates in a
!                rank 3 array with dimension i running from 1 to d(i).

!cod$
        f(1) = mod((p-1),d(1)) + 1
        f(2) = mod(int((p-1)/d(1)),d(2)) + 1
        f(3) = mod(int((p-1)/(d(1)*d(2))),d(3)) + 1
      end function

      function imap3d(d,i) result(f)
!doc$ function imap3d(d,i) result(f)
        integer, dimension(3), intent(in) :: d, i
        integer :: f
!       requires: all(1 <= i) .and. all(i <= d)
!       effects: Returns the ordinal number (from 1..product(d)) of the vector
!                 i in a rank three array with dimensions d in column ordering.
!cod$
        f = 1 + (i(1) - 1) + ((i(2) - 1) + (i(3) - 1)*d(2))*d(1)
      end function

      function nbhd_real(r,tol) result(f)
!doc$ function nbhd(r,tol) result(f)
        real(double), intent(in) :: r
        real(double), intent(in) :: tol
!       effects: Returns the tol neighborhood around r.

!cod$
        type(real_nbhd) :: f
        f%point = r
        f%tol = tol
      end function
      
      function nbhd_complex(c,tol) result(f)
!doc$ function nbhd(c,tol) result(f)
        complex(double), intent(in) :: c
        real(double), intent(in) :: tol
        type(complex_nbhd) :: f
!       effects: Returns the tol neighborhood around c.

!cod$
        f%point = c
        f%tol = tol
      end function
      
      function nbhd_vector(v,tol) result(f)
!doc$ function nbhd(v,tol) result(f)
        real(double), dimension(3), intent(in) :: v
        real(double), intent(in) :: tol
        type(real_vector_nbhd) :: f
!       effects: Returns the tol neighborhood around v.

!cod$
        f%point = v
        f%tol = tol
      end function
      
      function nbhd_matrix(m,tol) result(f)
!doc$ function nbhd(m,tol) result(f)
        real(double), dimension(3,3), intent(in) :: m
        real(double), intent(in) :: tol
        type(real_matrix_nbhd) :: f
!       effects: Returns the tol neighborhood around m.

!cod$
        f%point = m
        f%tol = tol
      end function
      
      function zone_vector(v,tol) result(f)
!doc$ function zone(v,tol) result(f)
        real(double), dimension(3), intent(in) :: v
        real(double), intent(in) :: tol
        type(real_vector_zone) :: f
!       effects: Returns the tol neighborhood around v and its integer translates.

!cod$
        f%point = v
        f%tol = tol
      end function
      
      function in_real(r,b) result(f)
!doc$ function in(r,b) result(f)
        real(double), intent(in) :: r
        type(real_nbhd), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff r is inside the b nbhd.

!cod$
        f = (abs(r-b%point) < b%tol)
      end function

      function in_complex(c,b) result(f)
!doc$ function in(c,b) result(f)
        complex(double), intent(in) :: c
        type(complex_nbhd), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff c is inside the b nbhd.

!cod$
        f = ( (abs(real(c)-real(b%point)) < b%tol) .and. (abs(aimag(c)-aimag(b%point)) < b%tol) )
      end function

      function in_vector(v,b) result(f)
!doc$ function in(v,b) result(f)
        real(double), dimension(3), intent(in) :: v
        type(real_vector_nbhd), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff v is inside the b nbhd.

!cod$
        f = all(abs(v-b%point) < b%tol)
      end function

      function in_matrix(m,b) result(f)
!doc$ function in(m,b) result(f)
        type(real_matrix_nbhd), intent(in) :: b
        real(double), intent(in) :: m(3,3)
        logical :: f
!       effects: Returns .true. iff m is inside the b nbhd.

!cod$
        logical, dimension(9) :: tmp
        tmp = reshape(abs(m-b%point) < b%tol,(/9/))
        f = all(tmp)
      end function

      function in_vector_zone(v,b) result(f)
!doc$ function in(v,b) result(f)
        real(double), dimension(3), intent(in) :: v
        type(real_vector_zone), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff v is inside the b zone.

!cod$
        real(double), dimension(3) :: tmp
        tmp = v - b%point
        f = all(abs(tmp-nint(tmp)) < b%tol)
      end function

      function out_real(r,b) result(f)
!doc$ function out(r,b) result(f)
        real(double), intent(in) :: r
        type(real_nbhd), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff r is on or outside the b nbhd.

!cod$
        f = (abs(r-b%point) >= b%tol)
      end function

      function out_complex(c,b) result(f)
!doc$ function out(c,b) result(f)
        complex(double), intent(in) :: c
        type(complex_nbhd), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff c is on or outside the b nbhd.

!cod$
        f = ( (abs(real(c)-real(b%point)) >= b%tol) .and. (abs(aimag(c)-aimag(b%point)) >= b%tol) )
      end function

      function out_vector(v,b) result(f)
!doc$ function out(v,b) result(f)
        real(double), dimension(3), intent(in) :: v
        type(real_vector_nbhd), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff v is on or outside the b nbhd.

!cod$
        f = all(abs(v-b%point) >= b%tol)
      end function

      function out_matrix(m,b) result(f)
!doc$ function out(m,b) result(f)
        type(real_matrix_nbhd), intent(in) :: b
        real(double), intent(in) :: m(3,3)
        logical :: f
!       effects: Returns .true. iff m is on or outside the b nbhd.

!cod$
        logical, dimension(9) :: tmp
        tmp = reshape(abs(m-b%point) >= b%tol,(/9/))
        f = all(tmp)
      end function

      function out_vector_zone(v,b) result(f)
!doc$ function out(v,b) result(f)
        real(double), dimension(3), intent(in) :: v
        type(real_vector_zone), intent(in) :: b
        logical :: f
!       effects: Returns .true. iff v is on or outside the b zone.

!cod$
        real(double), dimension(3) :: tmp
        tmp = v - b%point
        f = all(abs(tmp-nint(tmp)) >= b%tol)
      end function

      function even(n) result(f)
!doc$ function even(n) result(f)
        integer, intent(in) :: n
        logical :: f
!       effects: Returns .true. iff n is even.

!cod$
        f = (mod(n,2) == 0)
      end function
 
      function odd(n) result(f)
!doc$ function odd(n) result(f)
        integer, intent(in) :: n
        logical :: f
!       effects: Returns .true. iff n is odd.

!cod$
        f = (mod(n,2) /= 0)
      end function
 
      function norm_vector(v) result(f)
!doc$ function norm(v) result(f)
        real(double), dimension(3), intent(in) :: v
        real(double) :: f
!       effects: Returns sqrt(v(1)**2 + v(2)**2 + v(3)**2).

!cod$
        f = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
      end function
 
      function cross_product_vector(v1,v2) result(f)
!doc$ function cross_product(v1,v2) result(f)
        real(double), dimension(3), intent(in) :: v1, v2
        real(double), dimension(3) :: f
!       effects: Returns the cross product of v1 and v2.

!cod$
        f(1) = v1(2)*v2(3) - v1(3)*v2(2)
        f(2) = v1(3)*v2(1) - v1(1)*v2(3)
        f(3) = v1(1)*v2(2) - v1(2)*v2(1)
      end function

      function determinant_matrix(m) result(f)
!doc$ function determinant(m) result(f)
        real(double), dimension(3,3), intent(in) :: m
        real(double) :: f
!       effects: Returns the determinant of m.

!cod$
        f = m(1,1)*(m(2,2)*m(3,3) - m(2,3)*m(3,2)) + &
            m(1,2)*(m(2,3)*m(3,1) - m(2,1)*m(3,3)) + &
            m(1,3)*(m(2,1)*m(3,2) - m(2,2)*m(3,1))
      end function

      function inverse_matrix(m) result(f)
!doc$ function inverse(m) result(f)
        real(double), dimension(3,3), intent(in) :: m
        real(double), dimension(3,3) :: f
!       requires: Determinant of m not equal to zero.
!       effects: Returns the inverse of m.

!cod$
        real(double) :: det
        det = determinant(m)
        f(1,1) = (-m(2,3)*m(3,2) + m(2,2)*m(3,3))/det
        f(1,2) = ( m(1,3)*m(3,2) - m(1,2)*m(3,3))/det
        f(1,3) = (-m(1,3)*m(2,2) + m(1,2)*m(2,3))/det
        f(2,1) = ( m(2,3)*m(3,1) - m(2,1)*m(3,3))/det
        f(2,2) = (-m(1,3)*m(3,1) + m(1,1)*m(3,3))/det
        f(2,3) = ( m(1,3)*m(2,1) - m(1,1)*m(2,3))/det
        f(3,1) = (-m(2,2)*m(3,1) + m(2,1)*m(3,2))/det
        f(3,2) = ( m(1,2)*m(3,1) - m(1,1)*m(3,2))/det
        f(3,3) = (-m(1,2)*m(2,1) + m(1,1)*m(2,2))/det
      end function

      function trace_real_matrix(m) result(f)
!doc$ function trace(m) result(f)
        real(double), dimension(:,:), intent(in) :: m
        real(double) :: f
!       effects: Returns the trace of m.

!cod$
        integer :: i, n
        f = 0.0_double
        n = minval(shape(m))
        do i = 1,n
          f = f + m(i,i)
        end do
      end function
 
      function trace_complex_matrix(m) result(f)
!doc$ function trace(m) result(f)
        complex(double), dimension(:,:), intent(in) :: m
        complex(double) :: f
!       effects: Returns the trace of m.

!cod$
        integer :: i, n
        f = cmplx(0,0,double)
        n = minval(shape(m))
        do i = 1,n
          f = f + m(i,i)
        end do
      end function
 
      function solve_r(m,x) result(y)
!doc$ function solve(m,x) result(y)
        real(double), dimension(:,:), intent(in) :: m
        real(double), dimension(:) :: x
        real(double), dimension(size(m,2)) :: y
!       requires: m be nonsingular.
!       effects: Returns the solution y, of x = matmul(m,y).
!       errors: If m not square or m and x incompatible dimensions. Passes errors.
 
!cod$
	integer, dimension(size(m,2)) :: ipiv
        integer :: info
        real(double), dimension(size(m,1),size(m,2)) :: tmp
        if (error(size(m,1) /= size(m,2),"ERROR: matrix m must be square")) goto 100
        if (error(size(x,1) /= size(m,1),"ERROR: dimensions of x and m are incompatable")) goto 100
	y = x
        tmp = m
        call dgesv(size(m,1),1,tmp,size(m,1),ipiv,y,size(m,1),info)
        if (error(info /= 0,"ERROR: non-zero variable info in LAPACK call")) goto 100
100     if (error("Exit math_mod::solve_r")) continue
      end function
      
      function complementary_error(x) result(f)
!doc$ function complementary_error(x) result(f)
        real(double), intent(in) :: x
        real(double) :: f
!       effects: Returns the complementary error function evaluated at x.

!cod$
        real(double), parameter, dimension(13) :: erfcs = &
                      (/ -0.049046121234691808_double,     -0.142261205103713640_double, &
                          0.010035582187599796_double,     -0.000576876469976748_double, &
                          0.000027419931252196_double,     -0.000001104317550734_double, &
                          0.000000038488755420_double,     -0.000000001180858253_double, &
                          0.000000000032334215_double,     -0.000000000000799101_double, &
                          0.000000000000017990_double,     -0.000000000000000371_double, &
                          0.000000000000000007_double /)
        real(double), parameter, dimension(23) :: erc2cs = &
                      (/ -0.069601346602309501_double,     -0.041101339362620893_double, &
                          0.003914495866689626_double,     -0.000490639565054897_double, &
                          0.000071574790013770_double,     -0.000011530716341312_double, &
                          0.000001994670590201_double,     -0.000000364266647159_double, &
                          0.000000069443726100_double,     -0.000000013712209021_double, &
                          0.000000002788389661_double,     -0.000000000581416472_double, &
                          0.000000000123892049_double,     -0.000000000026906391_double, &
                          0.000000000005942614_double,     -0.000000000001332386_double, &
                          0.000000000000302804_double,     -0.000000000000069666_double, &
                          0.000000000000016208_double,     -0.000000000000003809_double, &
                          0.000000000000000904_double,     -0.000000000000000216_double, &
                          0.000000000000000052_double /)
        real(double), parameter, dimension(24) :: erfccs = &
                      (/  0.071517931020292500_double,     -0.026532434337606719_double, &
                          0.001711153977920853_double,     -0.000163751663458512_double, &
                          0.000019871293500549_double,     -0.000002843712412769_double, &
                          0.000000460616130901_double,     -0.000000082277530261_double, &
                          0.000000015921418724_double,     -0.000000003295071356_double, &
                          0.000000000722343973_double,     -0.000000000166485584_double, &
                          0.000000000040103931_double,     -0.000000000010048164_double, &
                          0.000000000002608272_double,     -0.000000000000699105_double, &
                          0.000000000000192946_double,     -0.000000000000054704_double, &
                          0.000000000000015901_double,     -0.000000000000004729_double, &
                          0.000000000000001432_double,     -0.000000000000000439_double, &
                          0.000000000000000138_double,     -0.000000000000000048_double /)
        integer :: nterf, nterfc, nterc2
        real(double) :: eta, sqeps, xmax, xsml, y
        real(double) :: eps = epsilon(1.0_double)
        real(double) :: smallest = tiny(1.0_double)
        eta = 1.0_double*eps
        nterf = inits_i(erfcs,eta)
        nterfc = inits_i(erfccs,eta)
        nterc2 = inits_i(erc2cs,eta)
        xsml = -sqrt(-log(sqrt(pi)*eps))
        xmax = sqrt(-log(sqrt(pi)*smallest))
        xmax = xmax - 0.5_double*log(xmax)/xmax - 0.01_double
        sqeps = sqrt(2.0_double*eps)
        if (x <= xsml) then
          f = 2.0_double
          return
        end if
        if (x <= xmax) then
          y = abs(x)
          if (y <= 1.0_double) then
            if (y < sqeps) f = 1.0_double - 2.0_double*x/sqrt(pi)
            if (y >= sqeps) then
              f = 1.0_double - x*(1.0_double + csevl_i(2.0_double*x*x-1.0_double,erfcs,nterf)) ; if (error()) goto 100
            end if
            return
          else
            y = y*y
            if (y <= 4.0_double) then
              f = exp(-y)/abs(x)*(0.5_double + csevl_i((8.0_double/y-5.0_double)/3.0_double,erc2cs,nterc2)) ; if (error()) goto 100
            end if
            if (y > 4.0_double) then
              f = exp(-y)/abs(x)*(0.5_double + csevl_i(8.0_double/y-1.0_double,erfccs,nterfc)) ; if (error()) goto 100
            end if
            if (x < 0.0_double) f = 2.0_double - f
            return
          end if
        end if
        f = 0.0_double
100     if (error("Exit math_mod::complementary_error")) continue
      end function

      function inits_i(os,eta) result(f)
        real(double), dimension(:), intent(in) :: os
        real(double), intent(in) :: eta
        integer :: f
!       warns: If eta is too small.
        integer :: i, ii, nos
        real(double) :: err
        nos = size(os)
        err = 0.0_double
        do ii = 1,nos
          i = nos + 1 - ii
          err = err + abs(os(i))
          if (err > eta) exit
        end do
        if (i == nos) call warn("eta may be too small")
 100    f = i
      end function

      function csevl_i(x,cs,n) result(f)
        real(double), intent(in) :: x
        real(double), dimension(:), intent(in) :: cs
        integer, intent(in) :: n
        real(double) :: f
!       requires: n <= size(cs)
!       errors: -1 <= x =< 1.
        integer :: i, ni
        real(double) :: b0, b1, b2, twox
        f = 0.0_double
        if (error(x < -1.0_double,"ERROR: x < -1")) goto 100
        if (error(x > +1.0_double,"ERROR: x > +1")) goto 100
        b1 = 0.0_double
        b0 = 0.0_double
        twox = 2.0_double*x
        do i = 1,n
          b2 = b1
          b1 = b0
          ni = n + 1 - i
          b0 = twox*b1 - b2 + cs(ni)
        end do
        f = 0.5_double*(b0 - b2)
100     if (error("Exit math_mod::csevl_i")) continue
      end function

      function spherical_bessel(x,l,switch) result(spb)
!doc$ function spherical_bessel(x,l,switch) result(spb)
        real(double), dimension(:), intent(in) :: x
        integer, intent(in) :: l
        logical, intent(in), optional :: switch
        real(double), dimension(size(x)) :: spb
!       requires: All x >= 0. l >= 0.
!       effects: Returns the l^th spherical Bessel function evaluated at values x.
!       errors: If tolerance is not reached.

!cod$
        integer, parameter :: n_limit = 500
        real(double), parameter :: tol_t = 1.0e-20_double
        logical :: divide
        integer :: i, m, n
        real(double) :: bs, bs0, bs1, x2, t

        divide = .false.
        if (present(switch)) divide = switch
        do i = 1,size(x)
          if (x(i) <= 0.5_double) then
            x2 = x(i)**2
            m = 2*l + 3
            t = -x2/real(2*m,double)
            bs = 1.0_double + t
            do  n = 2,n_limit
              m = m + 2
              t = -t*(x2/real(2*n*m,double))
              bs = bs + t
              if (abs(t) <= tol_t) exit
            end do
            if (error(n == n_limit,"ERROR: n_limit reached")) goto 100
            if (l > 0) then  ! Needed to properly handle the case x(i) = 0 when l = 0.
              do n = 1,l
                bs = bs/real(2*n+1,double)
              end do
              if (.not.divide) bs = bs*x(i)**l
            end if
          else
            bs = sin(x(i))/x(i)
            if (l > 0) then    
              bs0 = bs
              bs = (bs0 - cos(x(i)))/x(i)
              if (l > 1) then
                bs1 = bs
                do n = 2,l
                  bs = (real(2*n-1,double)*bs1)/x(i) - bs0
                  bs0 = bs1
                  bs1 = bs
                end do
              end if
              if (divide) bs = bs/x(i)**l
            end if
          end if       
          spb(i) = bs
        end do 

100     if (error("Exit math_mod::spherical_bessel")) continue

      end function

      function fermi_distr(x,kt) result(f)
!doc$ function fermi_distr(x,kt) result(f)
        real(double), intent(in) :: x
        real(double), intent(in), optional :: kt
        real(double) :: f
!       effects: Returns the Fermi function evaluated at x.
!                This is mathematically equivalent to 1/(1+exp(x))
!                or 1/(1+exp(x/kt))
!                This implementation will avoid overflows
!cod$
        real(double) :: xkt
        if (present(kt)) then
	   xkt = x/kt
	else
	   xkt = x
	endif
        if (xkt <= 0.0_double) then
           f = 1.0_double/(1.0_double + exp(xkt))
        else
           f = exp(-xkt)/(1.0_double + exp(-xkt))
        endif
      end function


      function fermi_entropy(x,kt) result(f)
!doc$ function fermi_entropy(x,kt) result(f)
        real(double), intent(in) :: x
        real(double), intent(in), optional :: kt
        real(double) :: f
!       effects: Returns the Fermi entropy evaluated at x.
!                This is mathematically equivalent to wlog(w) + (1-w)log(1-w),
!                where w = fermi_distr(x,kt)
!                This implementation will avoid overflows
!cod$   
        real(double) :: xkt, omw
        if (present(kt)) then
           xkt = x/kt 
        else
           xkt = x
        endif
        xkt = abs(xkt)
        omw = 1.0_double/(1.0_double+exp(-xkt))
        f = xkt*(1.0_double-omw)-log(omw)
      end function


      function fermi_distr_divdiff(x,y,kt) result(f)
!doc$ function fermi_distr_divdiff(x,y,kt) result(f)
        real(double), intent(in) :: x,y
        real(double), intent(in) :: kt
        real(double) :: f
!       effects: Returns the divided difference of the Fermi function
!                evaluated at x and y at temperature kt.
!                This is mathematically equivalent to
!                        (1/(1+exp(x/kt)) - 1/(1+exp(y/kt)))/(x-y)
!                Avoids overflows and inaccuracy when x ~= y
!cod$
! (1/(1+exp(x)) - 1/(1+exp(y)))/(x-y)
! == 1/((1+exp(x))(1+exp(y))) * (exp(y) - exp(x))/(x-y)
! xy = sort([x y]);
! xkt = xy(1)/kt; ykt = xy(2)/kt;
! (fermi_distr(xkt)-fermi_distr(ykt))/(kt*(xkt-ykt))
! f = -1.0/kt;
! if (xkt <= 0)
!   f = f*(1.0/(1.0+exp(xkt)));
! else
!   f = f*(exp(-xkt)/(1.0+exp(-xkt)));
! endif
! if (ykt <= 0)
!   f = f*(exp(ykt)/(1.0+exp(ykt)));
! else
!   f = f*(1.0/(1.0+exp(-ykt)));
! endif
! dkt = xkt - ykt;
! if (dkt < -1)
!   f = f*(exp(dkt)-1.0)/dkt;
! elseif (dkt < 0)
!   f = f*2*exp(dkt/2)*sinh(dkt/2)/dkt;
! endif
	real(double) :: xkt,ykt,dkt
	if (x <= y) then
	   xkt = x/kt
	   ykt = y/kt
	else
           xkt = y/kt
           ykt = x/kt
        end if

	f = -1.0/kt

	if (xkt <= 0.0_double) then
	   f = f*(1.0_double/(1.0_double+exp(xkt)))
	else
	   f = f*(exp(-xkt)/(1.0_double+exp(-xkt)))
	endif

	if (ykt <= 0.0_double) then
	   f = f*(exp(ykt)/(1.0_double+exp(ykt)))
	else
	   f = f*(1.0_double/(1.0_double+exp(-ykt)))
	endif

	dkt = xkt - ykt
	if (dkt < -1) then
	   f = f*((exp(dkt)-1.0_double)/dkt)
	elseif (dkt < 0) then
	   f = f*(2.0_double*exp(dkt/2.0_double)*sinh(dkt/2.0_double)/dkt)
	endif
      end function

      function random_sjp(seed) result(f)
!doc$ function random(seed) result(f)
        integer, intent(inout) :: seed
        real(double) :: f
!       modifies: seed
!       effects: Returns a random number based on seed.

!cod$
        real(double), parameter :: aa = 16807.0_double
        real(double), parameter :: mm = 2147483647.0_double
        real(double) :: real_seed
        real_seed = mod(aa*real(seed,double),mm)
        f = real_seed/mm
        seed = int(real_seed)
      end function

      subroutine inverse_cholesky(ans)
!doc$ subroutine inverse_cholesky(ans)
        complex(double), dimension(:,:), intent(inout) :: ans
!       modifies: ans
!       effects: Replaces ans with its inverse cholesky factor where inv(ans_in) = D*D', ans_out <-- D.

!cod$
        character(1) :: diag, uplo
        integer :: i, ierr, j, n

        n = size(ans,1)
        do i = 1,n
          do j = 1,(i-1)
            ans(i,j) = cmplx(0,0,double)
          end do
        end do
        uplo = 'u'
        call zpotrf(uplo,n,ans,n,ierr)
        if (error(ierr /= 0,"ERROR: zpotrf ierr /= 0")) goto 100
        diag = 'n'
        call ztrtri(uplo,diag,n,ans,n,ierr)
        if (error(ierr /= 0,"ERROR: ztrtri ierr /= 0")) goto 100

100     if (error("Exit math_mod::inverse_cholesky")) continue

      end subroutine

      function permutations(n,k) result(f)
!doc$ function permutations(n,k) result(f)
        integer, intent(in) :: n, k
        real(double) :: f
!       effects: Returns n!/(n-k)! if n >= 0 and n-k >= 0. Returns 0 otherwise.

!cod$
        integer :: i
        if ((n >= 0) .and. ((n-k) >= 0)) then
          f = 1.0_double
          do i = (n-k+1),n
            f = f*real(i,double)
          end do
        else
          f = 0.0_double
        end if
      end function

      function factorial(n) result(f)
!doc$ function factorial(n) result(f)
        integer, intent(in) :: n
        real(double) :: f
!       effects: Returns n!

!cod$
        integer :: i
        f = 1.0_double
        do i = 2,n
          f = f*real(i,double)
        end do
      end function

      function double_factorial(n) result(f)
!doc$ function double_factorial(n) result(f)
        integer, intent(in) :: n
        real(double) :: f
!       effects: Returns n!!

!cod$
        integer :: i
        f = 1.0_double
        do i = n,1,-2
          f = f*real(i,double)
        end do
      end function

      function spharm(x,y,z,l,do_norm) result(f)
!doc$ function spharm(x,y,z,l,do_norm) result(f)
        real(double), intent(in) :: x, y, z
        integer, intent(in) :: l
        logical, intent(in), optional :: do_norm
        complex(double), dimension(13) :: f
!       requires: size(f) = 13, l >= 0.
!       effects: Returns spherical harmonics for l in the order m = -l,...0,...+l
!                and in the Condon-Shortley convention. If do_norm = .false., the
!                results are multiplied by r**l where r = sqrt(x**2 + y**2 + z**2).
!       errors: If l > 6.

!cod$
        real(double), parameter :: tol = 1.0e-13_double
        integer :: m
        real(double) :: ylm_norm, r, costheta, cosphi, sinphi, tmp
        complex(double) :: phase

        if (error(l > 6,"ERROR: l > 6")) goto 100

        r = sqrt(x*x + y*y + z*z)
        if (present(do_norm)) then
          if (do_norm) then
            ylm_norm = sqrt(real(2*l+1,double)/four_pi)
          else
            ylm_norm = (r**l)*sqrt(real(2*l+1,double)/four_pi)
          end if
        end if

        f = cmplx(0,0,double)
        if (r < tol) then
          select case (l)
          case (0)
            f(1) = ylm_norm*cmplx(1,0,double)
          end select
          return
        end if
  
        select case (l)
        case (0)
          f(1) = ylm_norm*cmplx(1,0,double)
        case (1:6)
          costheta = z/r
          cosphi = 1.0_double
          sinphi = 0.0_double
          tmp = sqrt( (1.0_double - costheta)*(1.0_double + costheta) )
          if (tmp > tol) then
            cosphi = x/(r*tmp)
            sinphi = y/(r*tmp)
          end if
          f(l+1) = ylm_norm*assosp_i(l,0,costheta)
          phase = cmplx(cosphi,sinphi,double)
          tmp = 1.0_double/real(l*(l+1),double)
          do m = 1,l
            f(l+1+m) = ylm_norm*sqrt(tmp)*assosp_i(l,m,costheta)*(phase**m)
            f(l+1-m) = conjg(f(l+1+m))*(-1.0_double)**m
            if (m < l) tmp = tmp/real((l+m+1)*(l-m),double)
          end do
        end select

100     if (error("Exit math_mod::spharm")) continue

      end function

      subroutine gradylm(x,y,z,l,dtheta,dphi)
!doc$ subroutine gradylm(x,y,z,l,dtheta,dphi)
        real(double), intent(in) :: x, y, z
        integer, intent(in) :: l
        complex(double), dimension(13), intent(out):: dtheta, dphi
!       requires: size(dtheta) = 13, size(dphi) = 13, l >= 0.
!       effects: Returns gradients of spherical harmonics for l in the order m = -l,...0,...+l
!                and in the Condon-Shortley convention.
!       errors: If l > 6.

!cod$
        real(double), parameter :: tol = 1.0e-13_double
        integer :: m
        real(double) :: ylm_norm, r, costheta, cosphi, sinphi, tmp
        complex(double) :: phase, iphase

        if (error(l > 6,"ERROR: l > 6")) goto 100

        dtheta = cmplx(0,0,double)
        dphi = cmplx(0,0,double)

        r = sqrt(x*x + y*y + z*z)
        If (r < tol) return

        select case (l)
        case (0)
          return
        case (1:6)
          ylm_norm= sqrt(real(2*l+1,double)/four_pi)
          costheta = z/r
          cosphi = 1.0_double
          sinphi = 0.0_double
          tmp = sqrt( (1.0_double - costheta)*(1.0_double + costheta) )
          if (tmp > tol) then
            cosphi = x/(r*tmp)
            sinphi = y/(r*tmp)
          end if
          dtheta(l+1) = ylm_norm*der_theta_p_i(l,0,costheta)
          dphi(l+1) = cmplx(0,0,double)
          phase = cmplx(cosphi,sinphi,double)
          iphase = cmplx(-sinphi,cosphi,double)
          tmp = 1.0_double/real(l*(l+1),double)
          do m = 1,l
            dtheta(l+1+m) = ylm_norm*sqrt(tmp)*der_theta_p_i(l,m,costheta)*(phase**m)
            dtheta(l+1-m) = real((-1)**m,double)*conjg(dtheta(l+1+m))
            dphi(l+1+m) = ylm_norm*sqrt(tmp)*der_phi_p_i(l,m,costheta)*iphase*(phase**(m-1))
            dphi(l+1-m) = real((-1)**m,double)*conjg(dphi(l+1+m))
            if (m < l) tmp = tmp/real((l+m+1)*(l-m),double)
          end do
        end select

100     if (error("Exit math_mod::gradylm")) continue
  
      end subroutine

      function der_theta_p_i(l,m,x) result(f)
        integer, intent(in) :: l, m
        real(double), intent(in) :: x
        real(double) :: f
!       requires: 0 <= m <= l and |x| <= 1.
!       origin: Numerical Recipies, 2nd edition, page 247. Modified to return diff(P^m_l(x),x),theta)

        integer :: i, fact, ll
        real(double) :: dpll, pll, dpmm, pmm, dpmmp1, pmmp1, somx2, dosomx2

        pmm = 1.0_double
        dpmm = 1.0_double
        dosomx2 = 1.0_double
        somx2 = sqrt( (1.0_double - x)*(1.0_double + x) )
        if (m == 0) then
          dpmm = 0.0_double
        elseif (m > 0) then
          fact = 1
          do i = 1,m
            pmm = -pmm*real(fact,double)*somx2
            dpmm = -dpmm*real(fact,double)
            fact = fact + 2
          end do
           if (m > 1) then
            do i = 2,m
              dosomx2 = dosomx2*somx2
            end do
          end if  
          dpmm = dpmm*real(m,double)*x*dosomx2 
        end if

        if (l == m) then
          f = dpmm
        else
          pmmp1 = x*real(2*m+1,double)*pmm
          dpmmp1 = -real(2*m+1,double)*somx2*pmm + x*real(2*m+1,double)*dpmm
          if (l == (m+1)) then
            f = dpmmp1
          else
            do ll = (m+2),l
              pll = (x*real(2*ll-1,double)*pmmp1 - real(ll+m-1,double)*pmm)/real(ll-m,double)
              dpll = (-somx2*real(2*ll-1,double)*pmmp1 + x*real(2*ll-1,double)*dpmmp1 - real(ll+m-1,double)*dpmm)/real(ll-m,double)
              pmm = pmmp1
              pmmp1 = pll
              dpmm = dpmmp1
              dpmmp1 = dpll
            end do
            f = dpll
          end if
        end if

      end function

      function der_phi_p_i(l,m,x) result(f)
        integer, intent(in) :: l, m
        real(double), intent(in) :: x
        real(double) :: f
!       requires: 0 <= m <= l and |x| <= 1.
!       origin: Numerical Recipies, 2nd edition, page 247. Modified to return m*P_l(x)/sqrt((1-x^2)).

        integer :: i, fact, ll
        real(double) :: pll, pmm, pmmp1, somx2, dosomx2

        f = 0.0_double
        if (m == 0) return

        pmm = 1.0_double
        dosomx2 = 1.0_double
        if (m > 0) then
          somx2 = sqrt( (1.0_double - x)*(1.0_double + x) )
          fact = 1
          do i = 1,m
            pmm = -pmm*real(fact,double)
            fact = fact + 2
          end do
          if (m > 1) then
            do i = 2,m
              dosomx2 = somx2*dosomx2
            end do
          end if
          pmm = pmm*dosomx2
        end if

        if (l == m) then
          f = pmm*real(m,double)
        else
          pmmp1 = x*real(2*m+1,double)*pmm
          if (l == (m+1)) then
            f = pmmp1*real(m,double)
          else
            do ll = (m+2),l
              pll = (x*real(2*ll-1,double)*pmmp1 - real(ll+m-1,double)*pmm)/real(ll-m,double)
              pmm = pmmp1
              pmmp1 = pll
            end do
            f = pll*real(m,double)
          end if
        end if

      end function

      function assosp_i(l,m,x) result(f)
        integer, intent(in) :: l, m
        real(double), intent(in) :: x
        real(double) :: f
!       requires: 0 <= m <= l and |x| <= 1.
!       origin: Numerical Recipies, 2nd edition, page 247.

        integer :: i, fact, ll
        real(double) :: pll, pmm, pmmp1, somx2

        pmm = 1.0_double
        if (m > 0) then
          somx2 = sqrt( (1.0_double - x)*(1.0_double + x) )
          fact = 1
          do i = 1,m
            pmm = -pmm*real(fact,double)*somx2
            fact = fact + 2
          end do
        end if

        if (l == m) then 
          f = pmm
        else
          pmmp1 = x*real(2*m+1,double)*pmm
          if (l == (m+1)) then
            f = pmmp1
          else
            do ll = (m+2),l 
              pll = (x*real(2*ll-1,double)*pmmp1-real(ll+m-1,double)*pmm)/real(ll-m,double)
              pmm = pmmp1
              pmmp1 = pll
            end do
            f = pll
          end if
        end if

      end function

      function simpson_integral_real(f,h) result(s)
!doc$ function simpson_integral(f,h) result(s)
        real(double), dimension(:), intent(in) :: f
        real(double), intent(in) :: h
        real(double) :: s
!       requires: f be tabulated at points separated by equal intervals h.
!       effects: Returns the integral of f using Simpson's method.
 
!cod$
        integer :: i, j, n
  
        n = size(f)
        select case (n)
        case (0)
          s = 0.0_double
        case (1)
          s = h*f(1)
        case (2)
          s = h*(f(1) + f(2))/2.0_double
        case default
          s = f(1) + 4.0_double*f(2) + f(3)
          j = 2*((n-1)/2) + 1
          if (j > 3) then
            do i = 5,n,2
              s = s + f(i-2) + 4.0_double*f(i-1) + f(i)
            end do
          end if
          s = s*(h/3.0_double)
          if (n > j) s = s + h*(f(n-1) + f(n))/2.0_double 
        end select

      end function

      function simpson_integral_complex(f,h) result(s)
!doc$ function simpson_integral(f,h) result(s)
        complex(double), dimension(:), intent(in) :: f
        real(double), intent(in) :: h
        complex(double) :: s
!       requires: f be tabulated at points separated by equal intervals h.
!       effects: Returns the integral of f using Simpson's method.
 
!cod$
        integer :: i, j, n
  
        n = size(f)
        select case (n)
        case (0)
          s = 0.0_double
        case (1)
          s = h*f(1)
        case (2)
          s = h*(f(1) + f(2))/2.0_double
        case default
          s = f(1) + 4.0_double*f(2) + f(3)
          j = 2*((n-1)/2) + 1
          if (j > 3) then
            do i = 5,n,2
              s = s + f(i-2) + 4.0_double*f(i-1) + f(i)
            end do
          end if
          s = s*(h/3.0_double)
          if (n > j) s = s + h*(f(n-1) + f(n))/2.0_double 
        end select

      end function

      subroutine gauss_legendre(x1,x2,x,w)
!doc$ subroutine gauss_legendre(x1,x2,x,w)
        real(double), intent(in) :: x1, x2
        real(double), dimension(:), intent(out) :: x
        real(double), dimension(:), intent(out) :: w
!       requires: size(x) = size(w)
!       effects: Returns Gauss-Legendre integration points (x) and weights (w) on the interval [-1,1].
!       source: Numerical Recipies.

!cod$
        real(double), parameter :: eps = 3.0e-14_double
        real(double) :: p1, p2, p3, pp, xl, xm, z, z1
        integer :: i, j, m, n

        n = size(x)
        m = (n+1)/2
        xl = (x2 - x1)/2.0_double
        xm = (x2 + x1)/2.0_double
        do i = 1,m
          z = cos(pi*(real(i,double) - 0.25_double)/(real(n,double) + 0.5_double))
          z1 = z + 10.0_double*eps
          do while (abs(z-z1) > eps)
            p1 = 1.0_double
            p2 = 0.0_double
            do j = 1,n
              p3 = p2
              p2 = p1
              p1 = (z*p2*(2.0_double*real(j,double) - 1.0_double) - p3*(real(j,double) - 1.0_double))/real(j,double)
            end do
            pp = real(n,double)*(z*p1 - p2)/(z*z - 1.0_double)
            z1 = z
            z = z1 - p1/pp
          end do
          x(i) = xm - xl*z
          x(n+1-i) = xm + xl*z
          w(i) = 2.0_double*xl/((1.0_double  - z*z)*pp*pp)
          w(n+1-i) = w(i)
        end do

      end subroutine


      function all_atom_dot_product_1(fa) result(f)
!doc$ function atom_dot_product(fa) result(f)
        real(double), dimension(:,:), intent(in) :: fa
        real(double) :: f
!       effects: Returns the sum of the dot_products of the atom positions.
  
!cod$
        integer :: ia,na
        na = size(fa,2)
        f = 0.0_double
        do ia = 1,na
          f = f + dot_product(fa(:,ia),fa(:,ia))
        end do
      end function

      function all_atom_dot_product_2(fa,fb) result(f)
!doc$ function atom_dot_product(fa,fb) result(f)
        real(double), dimension(:,:), intent(in) :: fa, fb
        real(double) :: f
!       effects: Returns the sum of dot products between fa and fb.
!       errors: If fa and fb have incompatible sizes.
  
!cod$
        integer :: ia, na
        f = 0.0_double
        if (error(size(fa,2) /= size(fb,2),"ERROR: incompatable array sizes")) goto 100
        na = size(fa,2)
        do ia = 1,na
          f = f + dot_product(fa(:,ia),fb(:,ia))
        end do
100     if (error("Exit math_mod::all_atom_dot_product_2")) continue  
      end function


      subroutine init_constants()
!doc$ subroutine init_machine_constants()
!       effects: ?

!cod$
        real(double) :: a1, a2, a3

        machine_precision = 0.0_double
        a1 = 4.0_double/3.0_double
        do while (machine_precision == 0.0_double)
          a2 = a1 - 1.0_double
          a3 = a2 + a2 + a2
          machine_precision = abs(a3 - 1.0_double)
        end do

        machine_zero = machine_precision**4
        machine_infinity = 1.0_double/machine_zero
        minlogarg = machine_precision
        minlog = log(minlogarg)
        maxlogarg = 1.0_double/machine_precision
        maxlog = log(maxlogarg)
        minexparg = log(machine_precision)
        minexp = 0.0_double
        maxexparg = -log(machine_precision)
        maxexp = exp(maxexparg)

      end subroutine

      function spherical_harmonic(l,m,v) result(f)
!doc$ function spherical_harmonic(l,m,v) result(f)
        integer, intent(in) :: l, m
        real(double), dimension(3), intent(in) :: v
        complex(double) :: f
!       requires: -l <= m <= +l.
!       effects: Returns (|v|^l)*Ylm.
!       documents: spherical_harmonics.pdf

!cod$
        integer :: am
        real(double) :: pf
 
        am = abs(m)
        pf = sqrt(real(2*l+1,double)/four_pi)*sqrt(factorial(l-am)/factorial(l+am))
        If (m < 0) then
          f = pf*plm_i(l,am,v)*cmplx(-v(1),v(2),double)**am
        elseif (m == 0) then
          f = pf*plm_i(l,am,v)*cmplx(1,0,double)
        else
          f = pf*plm_i(l,am,v)*cmplx(v(1),v(2),double)**am
        end if
 
      end function
 
      function grad_spherical_harmonic(l,m,v) result(f)
!doc$ function grad_spherical_harmonic(l,m,v) result(f)
        integer, intent(in) :: l, m
        real(double), dimension(3), intent(in) :: v
        complex(double), dimension(3) :: f
!       requires: -l <= m <= +l.
!       effects: Returns grad{(|v|^l)*Ylm}.
!       documents: spherical_harmonics.pdf

!cod$
        integer :: am
        real(double) :: pf
 
        if (l == 0) then
          f = cmplx(0,0,double)
        else
          am = abs(m)
          pf = sqrt(real(2*l+1,double)/four_pi)*sqrt(factorial(l-am)/factorial(l+am))
          if (m < 0) then
            f = pf*(grad_plm_i(l,am,v)*cmplx(-v(1),v(2),double)**am &
                    - cmplx((/1,0,0/),(/0,-1,0/),double)*plm_i(l,am,v)*real(am,double)*cmplx(-v(1),v(2),double)**(am-1))
          elseif (m == 0) then
            f = pf*grad_plm_i(l,am,v)*cmplx(1,0,double)
          else
            f = pf*(grad_plm_i(l,am,v)*cmplx(v(1),v(2),double)**am &
                    + cmplx((/1,0,0/),(/0,1,0/),double)*plm_i(l,am,v)*real(am,double)*cmplx(v(1),v(2),double)**(am-1))
          end if
        end if

      end function

      function real_spherical_harmonic(l,m,v) result(f)
!doc$ function real_spherical_harmonic(l,m,v) result(f)
        integer, intent(in) :: l, m
        real(double), dimension(3), intent(in) :: v
        real(double) :: f
!       requires: -l <= m <= +l.
!       effects: Returns (|v|^l)*Klm.
!       documents: spherical_harmonics.pdf

!cod$
        integer :: am

        if (m < 0) then
          am = abs(m)
          f = sqrt(2.0_double)*aimag(spherical_harmonic(l,am,v))
        elseif (m == 0) then
          f = real(spherical_harmonic(l,m,v))
        else
          f = sqrt(2.0_double)*real(spherical_harmonic(l,m,v))
        end if

      end function

      function grad_real_spherical_harmonic(l,m,v) result(f)
!doc$ function grad_real_spherical_harmonic(l,m,v) result(f)
        integer, intent(in) :: l, m
        real(double), dimension(3), intent(in) :: v
        real(double), dimension(3) :: f
!       requires: -l <= m <= +l.
!       effects: Returns grad{(|v|^l)*Klm}.
!       documents: spherical_harmonics.pdf

!cod$
        integer :: am

        if (m < 0) then
          am = abs(m)
          f = sqrt(2.0_double)*aimag(grad_spherical_harmonic(l,am,v))
        elseif (m == 0) then
          f = real(grad_spherical_harmonic(l,m,v))
        else
          f = sqrt(2.0_double)*real(grad_spherical_harmonic(l,m,v))
        end if

      end function

      function plm_i(l,m,v) result(f)
        integer, intent(in) :: l, m
        real(double), dimension(3), intent(in) :: v
        real(double) :: f
!       requires: 0 <= m <= l.
!       effects: Returns (|v|^l)*Plm.
!       documents: spherical_harmonics.pdf

        integer :: ll
        real(double) :: pll, pmm, pmmp1, vn

        pmm = 1.0_double
        if (m > 0) pmm = real((-1)**m,double)*double_factorial(2*m-1)

        if (l == m) then 
          f = pmm
        else
          pmmp1 = v(3)*real(2*m+1,double)*pmm
          if (l == m+1) then
            f = pmmp1
          else
            vn = norm(v)
            do ll = m+2,l 
              pll = (v(3)*real(2*ll-1,double)*pmmp1 - vn**2*real(ll+m-1,double)*pmm)/real(ll-m,double)
              pmm = pmmp1
              pmmp1 = pll
            end do
            f = pll
          end if
        end if

      end function

      function grad_plm_i(l,m,v) result(f)
        integer, intent(in) :: l, m
        real(double), dimension(3), intent(in) :: v
        real(double), dimension(3) :: f
!       requires: 0 <= m <= l.
!       effects: Returns grad{(|v|^l)*Plm}.
!       documents: spherical_harmonics.pdf

        integer :: ll
        real(double) :: pll, pmm, pmmp1, vn
        real(double), dimension(3) :: dpmm, dpmmp1, dpll

        pmm = 1.0_double
        if (m > 0) pmm = real((-1)**m,double)*double_factorial(2*m-1)
        dpmm = real((/0,0,0/),double)

        if (l == m) then 
          f = dpmm
        else
          pmmp1 = v(3)*real(2*m+1,double)*pmm
          dpmmp1 = real((/0,0,1/),double)*real(2*m+1,double)*pmm
          if (l == m+1) then
            f = dpmmp1
          else
            vn = norm(v)
            do ll = m+2,l
              pll = (v(3)*real(2*ll-1,double)*pmmp1 - vn**2*real(ll+m-1,double)*pmm)/real(ll-m,double)
              dpll = (v(3)*real(2*ll-1,double)*dpmmp1 - vn**2*real(ll+m-1,double)*dpmm &
                      + real((/0,0,1/),double)*real(2*ll-1,double)*pmmp1 &
                      - 2.0_double*(/v(1),v(2),v(3)/)*real(ll+m-1,double)*pmm)/real(ll-m,double)
              pmm = pmmp1
              pmmp1 = pll
              dpmm = dpmmp1
              dpmmp1 = dpll
            end do
            f = dpll
          end if
        end if

      end function

      subroutine chebyshev_initialize(c,norm,cond,per)
!doc$ subroutine chebyshev_initialize(c,norm,cond,per)
        type(chebyshev_state) :: c
        real(double), intent(in) :: norm, cond
        integer, intent(in) :: per
!       requires: norm > 0, cond > 0, per > 0
!       effects: Initializes the state, c, of a chebyshev acceleration. norm and cond should be estimates of the norm
!                and condition number of the matrix being inverted, while per is a periodic reset parameter.  A good 
!                default for per is 2*sqrt(cond).

!cod$
        c%beta = ((cond - 1.0_double)/(cond + 1.0_double))**2
        c%alpha = 2.0_double/(norm + norm/cond)
        c%period = per
        c%phase = 0
      end subroutine

      subroutine chebyshev_accelerate(c,v_wt,g_wt)
!doc$ subroutine chebyshev_accelerate(c,v_wt,g_wt)
        type(chebyshev_state) :: c
        real(double), intent(out) :: v_wt, g_wt
!       requires: c to have been initialized by chebyshev_initialize.
!       effects: Returns two weights, v_wt and g_wt such that the update
!         v_i = v_wt * v_{i-1} + g_wt * (A * x_{i-1} - b)
!         x_i = x_{i-1} + v_i
!         makes x_i a chebyshev accelerated sequence of minimizers of the quadratic function: x'*A*x - 2*x'*b.

!cod$
        real(double) :: tau2
        if (c%phase < 0) c%phase = 0
        c%phase = mod(c%phase,c%period) + 1
        if (c%phase == 1) then
           c%tau0 = 1.0_double
           c%tau1 = 1.0_double
           v_wt = 0.0_double
           g_wt = -c%alpha
        else
           tau2 = 2.0_double*c%tau1 - c%tau0*c%beta
           v_wt = c%tau0*c%beta/tau2
           g_wt = -2.0_double*c%tau1*c%alpha/tau2
           c%tau0 = c%tau1
           c%tau1 = tau2
        end if
      end subroutine

      function bessj(n,x) result(r)
!doc$ function besjn(n,x) result(r)
        integer, intent(in) :: n
        real(double), intent(in) :: x
        real(double) :: r

!cod$
        integer, parameter :: iacc = 40
        real(double), parameter :: bigno = 1.0e+10_double
        real(double), parameter :: bigni = 1.0e-10_double
        integer :: j, jsum, m
        real(double) :: bj, bjm, bjp, sum, tox

        if (n == 0) then
          r = bessj0(x)
          return
        end if
        if (n == 1) then
          r = bessj1(x)
          return
        end if
        if (x == 0.0_double) then
          r = 0.0_double
          return
        end if
        tox = 2.0_double/x
        if (x > real(n,double)) then
          bjm = bessj0(x)
          bj  = bessj1(x)
          do j = 1,n-1
            bjp = j*tox*bj - bjm
            bjm = bj
            bj  = bjp
          end do
          r = bj
        else
          m = 2*((n + int(sqrt(real(iacc*n,double))))/2)
          r = 0.0_double
          jsum = 0
          sum = 0.0_double
          bjp = 0.0_double
          bj  = 1.0_double
          do j = m,1,-1
            bjm = j*tox*bj - bjp
            bjp = bj
            bj  = bjm
            if (abs(bj) > bigno) then
              bj  = bj*bigni
              bjp = bjp*bigni
              r = r*bigni
              sum = sum*bigni
            end if
            if (jsum /= 0) sum = sum + bj
            jsum = 1 - jsum
            if (j == n) r = bjp
          end do
          sum = 2.0_double*sum - bj
          r = r/sum
        end if
      end function

      function bessj0(x) result(bj0)
        real(double), intent(in) :: x
        real(double) :: bj0

        integer :: k, k0
        real(double) :: x2, r, t1, cu, p0, q0
        real(double), dimension(12) :: a, b
        
        
        x2 = x*x
        if (x == 0.0_double) then
          bj0 = 1.0_double
          return
        end if
        if (x <= 12.0_double) then
          bj0 = 1.0_double
          r = 1.0_double
          do k = 1,30
            r = -0.25_double*r*x2/(k*k)
            bj0 = bj0 + r
            if (abs(r)< abs(bj0)*1.0e-15_double) exit
          end do
          r = 1.0_double
        else
          data a/ -0.7031250000000000e-01_double, 0.1121520996093750e+00_double, &
     &            -0.5725014209747314e+00_double, 0.6074042001273483e+01_double, &
     &            -0.1100171402692467e+03_double, 0.3038090510922384e+04_double, &
     &            -0.1188384262567832e+06_double, 0.6252951493434797e+07_double, &
     &            -0.4259392165047669e+09_double, 0.3646840080706556e+11_double, &
     &            -0.3833534661393944e+13_double, 0.4854014686852901e+15_double /
          data b/  0.7324218750000000e-01_double, -0.2271080017089844e+00_double, &
     &             0.1727727502584457e+01_double, -0.2438052969955606e+02_double, &
     &             0.5513358961220206e+03_double, -0.1825775547429318e+05_double, &
     &             0.8328593040162893e+06_double, -0.5006958953198893e+08_double, &
     &             0.3836255180230433e+10_double, -0.3649010818849833e+12_double, &
     &             0.4218971570284096e+14_double, -0.5827244631566907e+16_double /
          k0 = 12
          if (x >= 35.0_double) k0 = 10
          if (x >= 50.0_double) k0 = 8
          t1 = x - 0.25_double*pi
          p0 = 1.0_double
          q0 = -0.125_double/x
          do k = 1,k0
            p0 = p0 + a(k)*x**(-2*k)
            q0 = q0 + b(k)*x**(-2*k-1)
          end do
          cu = sqrt(2.0_double/(pi*x))
          bj0 = cu*(p0*cos(t1) - q0*sin(t1))
        end if

      end function

      function bessj1(x) result(bj1)
        real(double), intent(in) :: x
        real(double) :: bj1

        integer :: k, k0
        real(double) :: x2, r, t2, cu, p1, q1
        real(double), dimension(12) :: a1(12), b1(12)
          
        x2 = x*x
        if (x == 0.0_double) then
          bj1 = 0.0_double
          return
        end if
        if (x <= 12.0_double) then
          bj1 = 1.0_double
          r = 1.0_double
          do k = 1,30
            r = -0.25_double*r*x2/(k*(k + 1.0_double))
            bj1 = bj1 + r
            if (abs(r)< abs(bj1)*1.0e-15_double) exit
          end do
          bj1 = 0.5_double*x*bj1
        else
          data a1/ 0.1171875000000000e+00_double, -0.1441955566406250e+00_double, &
     &             0.6765925884246826e+00_double, -0.6883914268109947e+01_double, &
     &             0.1215978918765359e+03_double, -0.3302272294480852e+04_double, &
     &             0.1276412726461746e+06_double, -0.6656367718817688e+07_double, &
     &             0.4502786003050393e+09_double, -0.3833857520742790e+11_double, &
     &             0.4011838599133198e+13_double, -0.5060568503314727e+15_double /
          data b1/ -0.1025390625000000e+00_double, 0.2775764465332031e+00_double, &
     &             -0.1993531733751297e+01_double, 0.2724882731126854e+02_double, &
     &             -0.6038440767050702e+03_double, 0.1971837591223663e+05_double, &
     &             -0.8902978767070678e+06_double, 0.5310411010968522e+08_double, &
     &             -0.4043620325107754e+10_double, 0.3827011346598605e+12_double, &
     &             -0.4406481417852278e+14_double, 0.6065091351222699e+16_double /
          k0 = 12
          if (x >= 35.0_double) k0 = 10
          if (x >= 50.0_double) k0 = 8
          t2 = x - 0.75_double*pi
          p1 = 1.0_double
          q1 = 0.375_double/x
          do k = 1,k0
            p1 = p1 + a1(k)*x**(-2*k)
            q1 = q1 + b1(k)*x**(-2*k-1)
          end do
          cu = sqrt(2_double/(pi*X))
          bj1 = cu*(p1*cos(t2) - q1*sin(t2))
        end if

      end function

      end module 

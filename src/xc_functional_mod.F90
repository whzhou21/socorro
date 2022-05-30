! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

      module xc_functional_mod
!doc$ module xc_functional_mod

!     functional_mod contains routines for computing native exchange-correlation quantities.

      use kind_mod
      use math_mod
      use error_mod
      use mpi_mod
      use xc_type_mod

!cod$
      implicit none
      private

!doc$
      public :: xcq_analytic
      public :: xcd_analytic
      public :: xce_analytic

!cod$
      interface xcq_analytic
        module procedure xcq_analytic_1d, xcq_analytic_3d
      end interface
      interface xcd_analytic
        module procedure xcd_analytic_1d, xcd_analytic_3d
      end interface
      interface xce_analytic
        module procedure xce_analytic_1d, xce_analytic_3d
      end interface

      contains

! public routines

      subroutine xcq_analytic_1d(etype,ctype,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xcq_analytic(etype,ctype,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        integer, intent(in) :: ctype, etype
        real(double), dimension(:), pointer :: n, dn, lp, exc, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Relevant pointers be allocated and correctly sized.
!       effects: Computes the exchange-correlation density and its derivatives using analytic formulas.

!cod$
        integer :: np
        real(double), dimension(:), pointer :: ns, zeta

        nullify( ns, zeta )

        np = size(n)

        select case (etype)
        case (E_LDA)
          select case (ctype)
          case (C_PZ)
            call xcq_lda_pz_1d_i(np,n,exc,dfxcdn)
          case (C_PW)
            call xcq_lda_pw_1d_i(np,n,exc,dfxcdn)
          end select
        case (E_LSDA)
          call spin_functions_1d_i(n,ns,zeta) ; if (error()) goto 100
          select case (ctype)
          case (C_PW)
            call xcq_lsda_pw_1d_i(np,ns,zeta,exc,dfxcdn)
          end select
        case (E_PW91)
          select case (ctype)
          case (C_PZ)
            call xcq_pw91_pz_1d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcq_pw91_pw_1d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          end select
        case (E_PBE)
          select case (ctype)
          case (C_PZ)
            call xcq_pbe_pz_1d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcq_pbe_pw_1d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          end select
        case (E_AM05)
          select case (ctype)
          case (C_PZ)
            call xcq_am05_pz_1d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcq_am05_pw_1d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          end select
        case (E_BLYP)
          call xcq_blyp_1d_i(np,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select

100     if (associated( ns )) deallocate( ns )
        if (associated( zeta )) deallocate( zeta )

        if (error("Exit functional1d_mod::xcq_analytic_1d")) continue

      end subroutine

      subroutine xcd_analytic_1d(etype,ctype,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xcd_analytic(etype,ctype,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
        integer, intent(in) :: ctype, etype
        real(double), dimension(:), pointer :: n, dn, lp, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Relevant pointers be allocated and correctly sized.
!       effects: Computes the exchange-correlation density derivatives using analytic formulas.

!cod$
        integer :: np
        real(double), dimension(:), pointer :: ns, zeta

        nullify( ns, zeta )

        np = size(n)

        select case (etype)
        case (E_LDA)
          select case (ctype)
          case (C_PZ)
            call xcd_lda_pz_1d_i(np,n,dfxcdn)
          case (C_PW)
            call xcd_lda_pw_1d_i(np,n,dfxcdn)
          end select
        case (E_LSDA)
          call spin_functions_1d_i(n,ns,zeta) ; if (error()) goto 100
          select case (ctype)
          case (C_PW)
            call xcd_lsda_pw_1d_i(np,ns,zeta,dfxcdn)
          end select
        case (E_PW91)
          select case (ctype)
          case (C_PZ)
            call xcd_pw91_pz_1d_i(np,n,dn,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcd_pw91_pw_1d_i(np,n,dn,dfxcdn,dfxcddnodn)
          end select
        case (E_PBE)
          select case (ctype)
          case (C_PZ)
            call xcd_pbe_pz_1d_i(np,n,dn,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcd_pbe_pw_1d_i(np,n,dn,dfxcdn,dfxcddnodn)
          end select
        case (E_AM05)
          select case (ctype)
          case (C_PZ)
            call xcd_am05_pz_1d_i(np,n,dn,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcd_am05_pw_1d_i(np,n,dn,dfxcdn,dfxcddnodn)
          end select
        case (E_BLYP)
          call xcd_blyp_1d_i(np,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select

100     if (associated( ns )) deallocate( ns )
        if (associated( zeta )) deallocate( zeta )

        if (error("Exit functional1d_mod::xcd_analytic_1d")) continue

      end subroutine

      subroutine xce_analytic_1d(etype,ctype,n,dn,lp,exc)
!doc$ subroutine xce_analytic(etype,ctype,n,dn,lp,exc)
        integer, intent(in) :: ctype, etype
        real(double), dimension(:), pointer :: n, dn, lp, exc
!       requires: Relevant pointers be allocated and correctly sized.
!       effect: Computes the exchange-correlation density using analytic formulas.

!cod$
        integer :: np
        real(double), dimension(:), pointer :: ns, zeta

        nullify( ns, zeta )

        np = size(n)

        select case (etype)
        case (E_LDA)
          select case (ctype)
          case (C_PZ)
            call xce_lda_pz_1d_i(np,n,exc)
          case (C_PW)
            call xce_lda_pw_1d_i(np,n,exc)
          end select
        case (E_LSDA)
          call spin_functions_1d_i(n,ns,zeta) ; if (error()) goto 100
          select case (ctype)
          case (C_PW)
            call xce_lsda_pw_1d_i(np,ns,zeta,exc)
          end select
        case (E_PW91)
          select case (ctype)
          case (C_PZ)
            call xce_pw91_pz_1d_i(np,n,dn,exc)
          case (C_PW)
            call xce_pw91_pw_1d_i(np,n,dn,exc)
          end select
        case (E_PBE)
          select case (ctype)
          case (C_PZ)
            call xce_pbe_pz_1d_i(np,n,dn,exc)
          case (C_PW)
            call xce_pbe_pw_1d_i(np,n,dn,exc)
          end select
        case (E_AM05)
          select case (ctype)
          case (C_PZ)
            call xce_am05_pz_1d_i(np,n,dn,exc)
          case (C_PW)
            call xce_am05_pw_1d_i(np,n,dn,exc)
          end select
        case (E_BLYP)
          call xce_blyp_1d_i(np,n,dn,lp,exc)
        end select

100     if (associated( ns )) deallocate( ns )
        if (associated( zeta )) deallocate( zeta )

        if (error("Exit functional1d_mod::xce_analytic_1d")) continue

      end subroutine

      subroutine xcq_analytic_3d(etype,ctype,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xcq_analytic(etype,ctype,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        integer, intent(in) :: ctype, etype
        real(double), dimension(:,:,:), pointer :: n, dn, lp, exc, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Relevant pointers be allocated and correctly sized.
!       effects: Computes the exchange-correlation density and its derivatives using analytic formulas.

!cod$
        integer, dimension(3) :: np
        real(double), dimension(:,:,:), pointer :: ns, zeta

        nullify( ns, zeta )

        np(1) = size(n,1) 
        np(2) = size(n,2) 
        np(3) = size(n,3)

        select case (etype)
        case (E_LDA)
          select case (ctype)
          case (C_PZ)
            call xcq_lda_pz_3d_i(np,n,exc,dfxcdn)
          case (C_PW)
            call xcq_lda_pw_3d_i(np,n,exc,dfxcdn)
          end select
        case (E_LSDA)
          call spin_functions_3d_i(n,ns,zeta) ; if (error()) goto 100
          select case (ctype)
          case (C_PW)
            call xcq_lsda_pw_3d_i(np,ns,zeta,exc,dfxcdn)
          end select
        case (E_PW91)
          select case (ctype)
          case (C_PZ)
            call xcq_pw91_pz_3d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcq_pw91_pw_3d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          end select
        case (E_PBE)
          select case (ctype)
          case (C_PZ)
            call xcq_pbe_pz_3d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcq_pbe_pw_3d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          end select
        case (E_AM05)
          select case (ctype)
          case (C_PZ)
            call xcq_am05_pz_3d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcq_am05_pw_3d_i(np,n,dn,exc,dfxcdn,dfxcddnodn)
          end select
        case (E_BLYP)
          call xcq_blyp_3d_i(np,n,dn,lp,exc,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select

100     if (associated( ns )) deallocate( ns )
        if (associated( zeta )) deallocate( zeta )

        if (error("Exit functional_mod::xcq_analytic_3d")) continue

      end subroutine

      subroutine xcd_analytic_3d(etype,ctype,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
!doc$ subroutine xcd_analytic(etype,ctype,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
        integer, intent(in) :: ctype, etype
        real(double), dimension(:,:,:), pointer :: n, dn, lp, dfxcdn, dfxcddnodn, dfxcdlapl
!       requires: Relevant pointers be allocated and correctly sized.
!       effect: Computes the exchange-correlation density derivatives using analytic formulas.

!cod$
        integer, dimension(3) :: np
        real(double), dimension(:,:,:), pointer :: ns, zeta

        nullify( ns, zeta )

        np(1) = size(n,1) 
        np(2) = size(n,2) 
        np(3) = size(n,3)

        select case (etype)
        case (E_LDA)
          select case (ctype)
          case (C_PZ)
            call xcd_lda_pz_3d_i(np,n,dfxcdn)
          case (C_PW)
            call xcd_lda_pw_3d_i(np,n,dfxcdn)
          end select
        case (E_LSDA)
          call spin_functions_3d_i(n,ns,zeta) ; if (error()) goto 100
          select case (ctype)
          case (C_PW)
            call xcd_lsda_pw_3d_i(np,ns,zeta,dfxcdn)
          end select
        case (E_PW91)
          select case (ctype)
          case (C_PZ)
            call xcd_pw91_pz_3d_i(np,n,dn,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcd_pw91_pw_3d_i(np,n,dn,dfxcdn,dfxcddnodn)
          end select
        case (E_PBE)
          select case (ctype)
          case (C_PZ)
            call xcd_pbe_pz_3d_i(np,n,dn,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcd_pbe_pw_3d_i(np,n,dn,dfxcdn,dfxcddnodn)
          end select
        case (E_AM05)
          select case (ctype)
          case (C_PZ)
            call xcd_am05_pz_3d_i(np,n,dn,dfxcdn,dfxcddnodn)
          case (C_PW)
            call xcd_am05_pw_3d_i(np,n,dn,dfxcdn,dfxcddnodn)
          end select
        case (E_BLYP)
          call xcd_blyp_3d_i(np,n,dn,lp,dfxcdn,dfxcddnodn,dfxcdlapl)
        end select

100     if (associated( ns )) deallocate( ns )
        if (associated( zeta )) deallocate( zeta )

        if (error("Exit functional_mod::xcd_analytic_3d")) continue

      end subroutine

      subroutine xce_analytic_3d(etype,ctype,n,dn,lp,exc)
!doc$ subroutine xce_analytic(etype,ctype,n,dn,lp,exc)
        integer, intent(in) :: ctype, etype
        real(double), dimension(:,:,:), pointer :: n, dn, lp, exc
!       requires: Relevant pointers be allocated and correctly sized.
!       effect: Computes the exchange-correlation density using analytic formulas.

!cod$
        integer, dimension(3) :: np
        real(double), dimension(:,:,:), pointer :: ns, zeta

        nullify( ns, zeta )

        np(1) = size(n,1) 
        np(2) = size(n,2) 
        np(3) = size(n,3)

        select case (etype)
        case (E_LDA)
          select case (ctype)
          case (C_PZ)
            call xce_lda_pz_3d_i(np,n,exc)
          case (C_PW)
            call xce_lda_pw_3d_i(np,n,exc)
          end select
        case (E_LSDA)
          call spin_functions_3d_i(n,ns,zeta) ; if (error()) goto 100
          select case (ctype)
          case (C_PW)
            call xce_lsda_pw_3d_i(np,ns,zeta,exc)
          end select
        case (E_PW91)
          select case (ctype)
          case (C_PZ)
            call xce_pw91_pz_3d_i(np,n,dn,exc)
          case (C_PW)
            call xce_pw91_pw_3d_i(np,n,dn,exc)
          end select
        case (E_PBE)
          select case (ctype)
          case (C_PZ)
            call xce_pbe_pz_3d_i(np,n,dn,exc)
          case (C_PW)
            call xce_pbe_pw_3d_i(np,n,dn,exc)
          end select
        case (E_AM05)
          select case (ctype)
          case (C_PZ)
            call xce_am05_pz_3d_i(np,n,dn,exc)
          case (C_PW)
            call xce_am05_pw_3d_i(np,n,dn,exc)
          end select
        case (E_BLYP)
          call xce_blyp_3d_i(np,n,dn,lp,exc)
        end select

100     if (associated( ns )) deallocate( ns )
        if (associated( zeta )) deallocate( zeta )

        if (error("Exit functional_mod::xce_analytic_3d")) continue

      end subroutine

! private routines

      subroutine xcq_lda_pz_1d_i(np,n,e,d)
        integer :: np
        real(double), dimension(np) :: n, e, d

        integer :: i
        real(double), parameter :: p01 = -0.09600_double   ! c1_ca
        real(double), parameter :: p02 = +0.06220_double   ! c2_ca
        real(double), parameter :: p03 = -0.02320_double   ! c3_ca
        real(double), parameter :: p04 = +0.00400_double   ! c4_ca
        real(double), parameter :: p05 = -0.28460_double   ! d1_ca
        real(double), parameter :: p06 = +1.05290_double   ! d2_ca
        real(double), parameter :: p07 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: rs, rsln, rssq

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p01 + p03*rs + (p02 + p04*rs)*rsln
            ndecdn = -one_third*(p02 + (p03 + p04)*rs + p04*rs*rsln)
          else
            ec = p05/(1.0_double + p06*rssq + p07*rs)
            ndecdn = one_sixth*p05*(p06*rssq + 2.0_double*p07*rs)/(1.0_double + p06*rssq + p07*rs)**2
          end if

          e(i) = ex + ec
          d(i) = four_thirds*ex + ec + ndecdn

        end do

      end subroutine

      subroutine xcq_lda_pw_1d_i(np,n,e,d)
        integer :: np
        real(double), dimension(np) :: n, e, d

        integer :: i
        real(double), parameter :: p01 = +0.0310907_double   ! a
        real(double), parameter :: p02 = +0.21370_double     ! a1
        real(double), parameter :: p03 = +7.59570_double     ! b1
        real(double), parameter :: p04 = +3.58760_double     ! b2
        real(double), parameter :: p05 = +1.63820_double     ! b3
        real(double), parameter :: p06 = +0.49294_double     ! b4
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: c1, c2, c3, c4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p01*c1*(p03 + c1*(p04 + c1*(p05 + p06*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p01*(1.0_double + p02*c1**2)*c3
          c4 = -(1.0_double + p02*c1**2)*c1*(p03 + c1*(2.0_double*p04 + c1*(3.0_double*p05 + 4.0_double*p06*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p01*(p02*c1**2*c3 + p01*c4)

          e(i) = ex + ec
          d(i) = four_thirds*ex + ec + ndecdn

        end do

      end subroutine

      subroutine xcq_lsda_pw_1d_i(np,ns,zeta,e,d)
        integer :: np
        real(double), dimension(np) :: ns, zeta, e, d

        integer :: i
        real(double), parameter :: p10 =  +0.0310907_double   ! a_0
        real(double), parameter :: p20 =  +0.2137000_double   ! a1_0
        real(double), parameter :: p30 =  +7.5957000_double   ! b1_0
        real(double), parameter :: p40 =  +3.5876000_double   ! b2_0
        real(double), parameter :: p50 =  +1.6382000_double   ! b3_0
        real(double), parameter :: p60 =  +0.4929400_double   ! b4_0
        real(double), parameter :: p11 =  +0.0155450_double   ! a_1
        real(double), parameter :: p21 =  +0.2054800_double   ! a1_1
        real(double), parameter :: p31 = +14.1189000_double   ! b1_1
        real(double), parameter :: p41 =  +6.1977000_double   ! b2_1
        real(double), parameter :: p51 =  +3.3662000_double   ! b3_1
        real(double), parameter :: p61 =  +0.6251700_double   ! b4_1
        real(double), parameter :: p12 =  +0.0168870_double   ! a_a
        real(double), parameter :: p22 =  +0.1112500_double   ! a1_a
        real(double), parameter :: p32 = +10.3570000_double   ! b1_a
        real(double), parameter :: p42 =  +3.6231000_double   ! b2_a
        real(double), parameter :: p52 =  +0.8802600_double   ! b3_a
        real(double), parameter :: p62 =  +0.4967100_double   ! b4_a
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: f01, f02, f03
        real(double) :: c1, c1_2, c20, c21, c22, c30, c31, c32
        real(double) :: dc20drs, dc21drs, dc22drs
        real(double) :: f, dfdz, d2fdz2
        real(double) :: ex, dexdrs, dexdz
        real(double) :: ec, decdrs, decdz
        real(double) :: ec0, ec1, ec2, dec0drs, dec1drs, dec2drs
        real(double) :: omz, opz, s
        real(double) :: z, z3, z4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 2.0_double**four_thirds - 2.0_double

        d2fdz2 = (8.0_double/9.0_double)/f03  ! d2f/dz2 | z=0

        if (mpi_mysgroup() == 1) s = -1.0_double
        if (mpi_mysgroup() == 2) s = +1.0_double

        do i = 1,np

          z = zeta(i)
          z3 = z**3
          z4 = z**4

          opz = 1.0_double + z
          omz = 1.0_double - z

          f = (opz**four_thirds + omz**four_thirds - 2.0_double)/f03
          dfdz =  four_thirds*(opz**one_third - omz**one_third)/f03

          c1 = f02/ns(i)**one_sixth  ! sqrt(rs)
          c1_2 = c1**2               ! rs

          ! exchange
          ex = 0.50_double*(f01*ns(i)**one_third)*(opz**four_thirds + omz**four_thirds )
          dexdrs = -ex/c1_2
          dexdz = 0.50_double*(f01*ns(i)**one_third)*four_thirds*(opz**one_third - omz**one_third )

          ! correlation
          c20 = 2.0_double*p10*c1*(p30 + c1*(p40 + c1*(p50 + p60*c1)))
          c21 = 2.0_double*p11*c1*(p31 + c1*(p41 + c1*(p51 + p61*c1)))
          c22 = 2.0_double*p12*c1*(p32 + c1*(p42 + c1*(p52 + p62*c1)))
          c30 = log(1.0_double + 1.0_double/c20)
          c31 = log(1.0_double + 1.0_double/c21)
          c32 = log(1.0_double + 1.0_double/c22)
          ec0 = -4.0_double*p10*(1.0_double + p20*c1_2)*c30
          ec1 = -4.0_double*p11*(1.0_double + p21*c1_2)*c31
          ec2 = -4.0_double*p12*(1.0_double + p22*c1_2)*c32
          ec = ec0 + ec2*(f/d2fdz2)*(1.0_double - z4) + (ec1 - ec0)*f*z4
          dc20drs = (p10/c1)*(p30 + c1*(2.0_double*p40 + c1*(3.0_double*p50 + c1*4.0_double*p60)))
          dc21drs = (p11/c1)*(p31 + c1*(2.0_double*p41 + c1*(3.0_double*p51 + c1*4.0_double*p61)))
          dc22drs = (p12/c1)*(p32 + c1*(2.0_double*p42 + c1*(3.0_double*p52 + c1*4.0_double*p62)))
          dec0drs = 4.0_double*p10*(dc20drs*(1.0_double + p20*c1_2)/(c20*(1.0_double + c20)) - p20*c30)
          dec1drs = 4.0_double*p11*(dc21drs*(1.0_double + p21*c1_2)/(c21*(1.0_double + c21)) - p21*c31)
          dec2drs = 4.0_double*p12*(dc22drs*(1.0_double + p22*c1_2)/(c22*(1.0_double + c22)) - p22*c32)
          decdrs = dec0drs + dec2drs*f*(1.0_double - z4)/d2fdz2 + (dec1drs - dec0drs)*f*z4
          decdz = ec2*(dfdz*(1.0_double - z4) + f*(1.0_double - 4.0_double*z3))/d2fdz2 + (ec1 - ec0)*(dfdz*z4 + 4.0_double*f*z3)

          ! combine
          e(i) = ex + ec
          d(i) = ex + ec - (dexdrs + decdrs)*one_third*c1_2 + (dexdz + decdz)*(1.0_double + s*z)

        end do

100     if (error("Exit functional1d_mod::xcq_lsda_pw_1d_i")) continue

      end subroutine

      subroutine xcq_pw91_pz_1d_i(np,n,dn,e,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, e, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = -0.09600_double       ! c1_ca
        real(double), parameter :: p08 = +0.06220_double       ! c2_ca
        real(double), parameter :: p09 = -0.02320_double       ! c3_ca
        real(double), parameter :: p10 = +0.00400_double       ! c4_ca
        real(double), parameter :: p11 = -0.28460_double       ! d1_ca
        real(double), parameter :: p12 = +1.05290_double       ! d2_ca
        real(double), parameter :: p13 = +0.33340_double       ! d3_ca
        real(double), parameter :: p14 = +0.002568_double      ! c0
        real(double), parameter :: p15 = +0.023266_double      ! c1
        real(double), parameter :: p16 = +7.389e-6_double      ! c2
        real(double), parameter :: p17 = +8.723_double         ! c3
        real(double), parameter :: p18 = +0.472_double         ! c4
        real(double), parameter :: p19 = +7.389e-2_double      ! c5
        real(double), parameter :: p20 = -0.001667212_double   ! cx
        real(double), parameter :: p21 = +0.004235_double      ! cc0
        real(double), parameter :: p22 = +0.09_double          ! alfa
        real(double), parameter :: p23 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: rs, rsln, rssq, s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p21*p23                                                 ! beta
        f06 = 2.0_double*p22/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p07 + p09*rs + (p08 + p10*rs)*rsln
            ndecdn = -one_third*(p08 + (p09 + p10)*rs + p10*rs*rsln)
          else
            ec = p11/(1.0_double + p12*rssq + p13*rs)
            ndecdn = one_sixth*p11*(p12*rssq + 2.0_double*p13*rs)/(1.0_double + p12*rssq + p13*rs)**2
          end if
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i)**one_third
          c2 = 100.0_double*t2*f09/n(i)**one_third
          c3 = exp(-c2)
          c4 = -p20 + (p14 + p15*c1 + p16*c1**2)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)
          h1 = 2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*t2*c3
          c5 = 2.0_double*(p14*p18 - p16) + c1*(p15*p18 - p16*p17 + p19*(3.0_double*p14 + c1*(2.0_double*p15 + p16*c1)))
          c5 = c1*(p14*p17 - p15 + c1*c5)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p23*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*c3*(1.0_double - c2)

          e(i) = ex*fx + ec + h0 + h1
          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i) = f07*ex/n(i)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcq_pw91_pw_1d_i(np,n,dn,e,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, e, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = + 0.0310907_double    ! a
        real(double), parameter :: p08 = +0.21370_double       ! a1
        real(double), parameter :: p09 = +7.59570_double       ! b1
        real(double), parameter :: p10 = +3.58760_double       ! b2
        real(double), parameter :: p11 = +1.63820_double       ! b3
        real(double), parameter :: p12 = +0.49294_double       ! b4
        real(double), parameter :: p13 = +0.002568_double      ! c0
        real(double), parameter :: p14 = +0.023266_double      ! c1
        real(double), parameter :: p15 = +7.389e-6_double      ! c2
        real(double), parameter :: p16 = +8.723_double         ! c3
        real(double), parameter :: p17 = +0.472_double         ! c4
        real(double), parameter :: p18 = +7.389e-2_double      ! c5
        real(double), parameter :: p19 = -0.001667212_double   ! cx
        real(double), parameter :: p20 = +0.004235_double      ! cc0
        real(double), parameter :: p21 = +0.09_double          ! alfa
        real(double), parameter :: p22 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p20*p22                                                 ! beta
        f06 = 2.0_double*p21/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p07*c1*(p09 + c1*(p10 + c1*(p11 + p12*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p08*c1**2)*c1
          c4 = c4*(p09 + c1*(2.0_double*p10 + c1*(3.0_double*p11 + 4.0_double*p12*c1)))/(c2*(1.0_double + c2))
          ec = -4.0_double*p07*(1.0_double + p08*c1**2)*c3
          ndecdn = four_thirds*p07*(p08*c1**2*c3 + p07*c4)
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i)**one_third
          c2 = 100.0_double*t2*f09/n(i)**one_third
          c3 = exp(-c2)
          c4 = -p19 + (p13 + p14*c1 + p15*c1**2)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)
          h1 = 2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*t2*c3
          c5 = 2.0_double*(p13*p17 - p15) + c1*(p14*p17 - p15*p16 + p18*(3.0_double*p13 + c1*(2.0_double*p14 + p15*c1)))
          c5 = c1*(p13*p16 - p14 + c1*c5)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p22*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*c3*(1.0_double - c2)

          e(i) = ex*fx + ec + h0 + h1
          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i) = f07*ex/n(i)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcq_pbe_pz_1d_i(np,n,dn,e,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, e, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.804_double      ! kappa
        real(double), parameter :: p02 = +0.21951_double    ! mu
        real(double), parameter :: p03 = +0.031091_double   ! gamma
        real(double), parameter :: p04 = +0.066725_double   ! beta
        real(double), parameter :: p05 = -0.09600_double    ! c1_ca
        real(double), parameter :: p06 = +0.06220_double    ! c2_ca
        real(double), parameter :: p07 = -0.02320_double    ! c3_ca
        real(double), parameter :: p08 = +0.00400_double    ! c4_ca
        real(double), parameter :: p09 = -0.28460_double    ! d1_ca
        real(double), parameter :: p10 = +1.05290_double    ! d2_ca
        real(double), parameter :: p11 = +0.33340_double    ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: rs, rsln, rssq, s, s2, t, t2, c1, c2, c3, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p05 + p07*rs + (p06 + p08*rs)*rsln
            ndecdn = -one_third*(p06 + (p07 + p08)*rs + p08*rs*rsln)
          else
            ec = p09/(1.0_double + p10*rssq + p11*rs)
            ndecdn = one_sixth*p09*(p10*rssq + 2.0_double*p11*rs)/(1.0_double + p10*rssq + p11*rs)**2
          end if
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c2 = p04/p03/(c1 - 1.0_double)*t2
          c3 = (1.0_double + c2*(1.0_double + c2))**2 + p04/p03*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndhdadadn = -c1*c2**3*ndecdn*(2.0_double + c2)/c3
          dhdtot = p04*(4.0_double + 8.0_double*c2)/c3

          e(i) = ex*fx + ec + h
          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i) = f05*ex/n(i)**five_thirds*dfxdsos + f06*dhdtot/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcq_pbe_pw_1d_i(np,n,dn,e,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, e, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.804_double       ! kappa
        real(double), parameter :: p02 = +0.21951_double     ! mu
        real(double), parameter :: p03 = +0.031091_double    ! gamma
        real(double), parameter :: p04 = +0.066725_double    ! beta
        real(double), parameter :: p05 = +0.0310907_double   ! a
        real(double), parameter :: p06 = +0.21370_double     ! a1
        real(double), parameter :: p07 = +7.59570_double     ! b1
        real(double), parameter :: p08 = +3.58760_double     ! b2
        real(double), parameter :: p09 = +1.63820_double     ! b3
        real(double), parameter :: p10 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: s, s2, t, t2, c1, c2, c3, c4, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p05*c1*(p07 + c1*(p08 + c1*(p09 + p10*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p06*c1**2)*c1*(p07 + c1*(2.0_double*p08 + c1*(3.0_double*p09 + 4.0_double*p10*c1)))/ &
                 (c2*(1.0_double + c2))
          ec = -4.0_double*p05*(1.0_double + p06*c1**2)*c3
          ndecdn = four_thirds*p05*(p06*c1**2*c3 + p05*c4)
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c3 = p04/p03/(c1 - 1.0_double)*t2
          c4 = (1.0_double + c3*(1.0_double + c3))**2 + p04/p03*t2*(1.0_double + c3)*(1.0_double + c3*(1.0_double + c3))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c3)/(1.0_double + c3*(1.0_double + c3)))
          ndhdadadn = -c1*c3**3*ndecdn*(2.0_double + c3)/c4
          dhdtot = p04*(4.0_double + 8.0_double*c3)/c4

          e(i) = ex*fx + ec + h
          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i) = f05*ex/n(i)**five_thirds*dfxdsos + f06*dhdtot/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcq_am05_pz_1d_i(np,n,dn,e,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, e, d, dd

        integer :: i, j
        real(double), parameter :: p01 = +2.804_double     ! alpha
        real(double), parameter :: p02 = +0.8098_double    ! gamma
        real(double), parameter :: p03 = +0.7168_double    ! c
        real(double), parameter :: p04 = -0.09600_double   ! c1_ca
        real(double), parameter :: p05 = +0.06220_double   ! c2_ca
        real(double), parameter :: p06 = -0.02320_double   ! c3_ca
        real(double), parameter :: p07 = +0.00400_double   ! c4_ca
        real(double), parameter :: p08 = -0.28460_double   ! d1_ca
        real(double), parameter :: p09 = +1.05290_double   ! d2_ca
        real(double), parameter :: p10 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: rs, rsln, rssq, s, s2, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do j = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p04 + p06*rs + (p05 + p07*rs)*rsln
            ndecdn = -one_third*(p05 + (p06 + p07)*rs + p07*rs*rsln)
          else
            ec = p08/(1.0_double + p09*rssq + p10*rs)
            ndecdn = one_sixth*p08*(p09*rssq + 2.0_double*p10*rs)/(1.0_double + p09*rssq + p10*rs)**2
          end if

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          e(i) = ex*hx + ec*hc
          d(i) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i))

        end do

200     call sync_config_process_errors()
        if (error("Exit functional1d_mod::xcq_am05_pz_1d_i")) continue

      end subroutine

      subroutine xcq_am05_pw_1d_i(np,n,dn,e,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, e, d, dd

        integer :: i, j
        real(double), parameter :: p01 = +2.804_double       ! alpha
        real(double), parameter :: p02 = +0.8098_double      ! gamma
        real(double), parameter :: p03 = +0.7168_double      ! c
        real(double), parameter :: p04 = +0.0310907_double   ! a
        real(double), parameter :: p05 = +0.21370_double     ! a1
        real(double), parameter :: p06 = +7.59570_double     ! b1
        real(double), parameter :: p07 = +3.58760_double     ! b2
        real(double), parameter :: p08 = +1.63820_double     ! b3
        real(double), parameter :: p09 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: s, s2, c1, c2, c3, c4, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do j = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p04*c1*(p06 + c1*(p07 + c1*(p08 + p09*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p04*(1.0_double + p05*c1**2)*c3
          c4 = -(1.0_double + p05*c1**2)*c1*(p06 + c1*(2.0_double*p07 + c1*(3.0_double*p08 + 4.0_double*p09*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p04*(p05*c1**2*c3 + p04*c4)

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          e(i) = ex*hx + ec*hc
          d(i) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i))

        end do

200     call sync_config_process_errors()
        if (error("Exit functional1d_mod::xcq_am05_pw_1d_i")) continue

      end subroutine

      subroutine xcq_blyp_1d_i(np,n,dn,lp,e,d,dd,ddd)
        integer :: np
        real(double), dimension(np) :: n, dn, lp, e, d, dd, ddd

        integer :: i
        real(double), parameter :: p01 = +0.27430_double    ! c1
        real(double), parameter :: p02 = +0.19645_double    ! c2
        real(double), parameter :: p03 = +7.79560_double    ! c3
        real(double), parameter :: p04 = +0.04918_double    ! a
        real(double), parameter :: p05 = +0.132_double      ! b
        real(double), parameter :: p06 = +0.2533_double     ! c
        real(double), parameter :: p07 = +0.349_double      ! d
        real(double), parameter :: p08 = +2.871234_double   ! cf
        real(double) :: f01, f02, f03
        real(double) :: one_eighth, one_third, two_thirds, four_thirds, five_thirds, seventeen_thirds
        real(double) :: ex, fx, tb, dfxdsos, ec, dfcdn, dfcddnodn, dfcdd2n
        real(double) :: s, s2, c1, c2, c3, c4, c5, x1, x2, x3

        one_eighth = 1.0_double/8.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        seventeen_thirds = 17.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f03 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f02*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          fx = 1.0_double + p01*s2/x2
          x3 = 1.0_double - p02*p03*s2/x1
          tb = 1.0_double - p01*s2*x3/x2**2
          dfxdsos = p01*(x3/x2 + 1.0_double)/x2

          ! correlation
          c1 = dn(i)/n(i)**four_thirds
          c2 = lp(i)/n(i)**five_thirds
          c3 = n(i)**one_third
          c4 = one_eighth*(c1**2 - c2)
          c5 = (p08*p05 + p05/9.0_double*(c2/2.0_double - 17.0_double*c4))*exp(-p06/c3)
          ec = -2.0_double*p04/(1.0_double + p07/c3)*(1.0_double + c5)
          c4 = seventeen_thirds*(p06*p07 + c3*(p06 - 4.0_double*p07 - 5.0_double*c3))*c1**2
          c4 = c4 - 7.0_double*(p06*p07 + c3*(p06 - p07 - 2.0_double*c3))*c2
          c5 = p05*exp(-p06/c3)
          c4 = c4*c5/24.0_double
          c5 = c5*p08*(p06*p07 + c3*(p06 + 4.0_double*p07 + 3.0_double*c3))
          c4 = c4 - c5 - c3*(4.0_double*p07 + 3.0_double*c3)
          dfcdn = 2.0_double*p04*c4/(3.0_double*(p07 + c3)**2)
          c4 = p05/c3**5/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcddnodn = 34.0_double*p04/36.0_double*c4
          c4 = p05/c3**2/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcdd2n = -14.0_double*p04/24.0_double*c4

          e(i) = ex*fx + ec
          d(i) = four_thirds*ex*tb + dfcdn
          dd(i) = f03*ex/n(i)**five_thirds*dfxdsos + dfcddnodn
          ddd(i) = dfcdd2n

        end do

      end subroutine

      subroutine xcd_lda_pz_1d_i(np,n,d)
        integer :: np
        real(double), dimension(np) :: n, d

        integer :: i
        real(double), parameter :: p01 = -0.09600_double   ! c1_ca
        real(double), parameter :: p02 = +0.06220_double   ! c2_ca
        real(double), parameter :: p03 = -0.02320_double   ! c3_ca
        real(double), parameter :: p04 = +0.00400_double   ! c4_ca
        real(double), parameter :: p05 = -0.28460_double   ! d1_ca
        real(double), parameter :: p06 = +1.05290_double   ! d2_ca
        real(double), parameter :: p07 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: rs, rsln, rssq

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p01 + p03*rs + (p02 + p04*rs)*rsln
            ndecdn = -one_third*(p02 + (p03 + p04)*rs + p04*rs*rsln)
          else
            ec = p05/(1.0_double + p06*rssq + p07*rs)
            ndecdn = one_sixth*p05*(p06*rssq + 2.0_double*p07*rs)/(1.0_double + p06*rssq + p07*rs)**2
          end if

          d(i) = four_thirds*ex + ec + ndecdn

        end do

      end subroutine

      subroutine xcd_lda_pw_1d_i(np,n,d)
        integer :: np
        real(double), dimension(np) :: n, d

        integer :: i
        real(double), parameter :: p01 = +0.0310907_double   ! a
        real(double), parameter :: p02 = +0.21370_double     ! a1
        real(double), parameter :: p03 = +7.59570_double     ! b1
        real(double), parameter :: p04 = +3.58760_double     ! b2
        real(double), parameter :: p05 = +1.63820_double     ! b3
        real(double), parameter :: p06 = +0.49294_double     ! b4
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: c1, c2, c3, c4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p01*c1*(p03 + c1*(p04 + c1*(p05 + p06*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p01*(1.0_double + p02*c1**2)*c3
          c4 = -(1.0_double + p02*c1**2)*c1*(p03 + c1*(2.0_double*p04 + c1*(3.0_double*p05 + 4.0_double*p06*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p01*(p02*c1**2*c3 + p01*c4)

          d(i) = four_thirds*ex + ec + ndecdn

        end do

      end subroutine

      subroutine xcd_lsda_pw_1d_i(np,ns,zeta,d)
        integer :: np
        real(double), dimension(np) :: ns, zeta, d

        integer :: i
        real(double), parameter :: p10 =  +0.0310907_double   ! a_0
        real(double), parameter :: p20 =  +0.2137000_double   ! a1_0
        real(double), parameter :: p30 =  +7.5957000_double   ! b1_0
        real(double), parameter :: p40 =  +3.5876000_double   ! b2_0
        real(double), parameter :: p50 =  +1.6382000_double   ! b3_0
        real(double), parameter :: p60 =  +0.4929400_double   ! b4_0
        real(double), parameter :: p11 =  +0.0155450_double   ! a_1
        real(double), parameter :: p21 =  +0.2054800_double   ! a1_1
        real(double), parameter :: p31 = +14.1189000_double   ! b1_1
        real(double), parameter :: p41 =  +6.1977000_double   ! b2_1
        real(double), parameter :: p51 =  +3.3662000_double   ! b3_1
        real(double), parameter :: p61 =  +0.6251700_double   ! b4_1
        real(double), parameter :: p12 =  +0.0168870_double   ! a_a
        real(double), parameter :: p22 =  +0.1112500_double   ! a1_a
        real(double), parameter :: p32 = +10.3570000_double   ! b1_a
        real(double), parameter :: p42 =  +3.6231000_double   ! b2_a
        real(double), parameter :: p52 =  +0.8802600_double   ! b3_a
        real(double), parameter :: p62 =  +0.4967100_double   ! b4_a
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: f01, f02, f03
        real(double) :: c1, c1_2, c20, c21, c22, c30, c31, c32
        real(double) :: dc20drs, dc21drs, dc22drs
        real(double) :: f, dfdz, d2fdz2
        real(double) :: ex, dexdrs, dexdz
        real(double) :: ec, decdrs, decdz
        real(double) :: ec0, ec1, ec2, dec0drs, dec1drs, dec2drs
        real(double) :: omz, opz, s
        real(double) :: z, z3, z4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 2.0_double**four_thirds - 2.0_double

        d2fdz2 = (8.0_double/9.0_double)/f03  ! d2f/dz2 | z=0

        if (mpi_mysgroup() == 1) s = -1.0_double
        if (mpi_mysgroup() == 2) s = +1.0_double

        do i = 1,np

          z = zeta(i)
          z3 = z**3
          z4 = z**4

          opz = 1.0_double + z
          omz = 1.0_double - z

          f = (opz**four_thirds + omz**four_thirds - 2.0_double)/f03
          dfdz =  four_thirds*(opz**one_third - omz**one_third)/f03

          c1 = f02/ns(i)**one_sixth  ! sqrt(rs)
          c1_2 = c1**2               ! rs

          ! exchange
          ex = 0.50_double*(f01*ns(i)**one_third)*(opz**four_thirds + omz**four_thirds )
          dexdrs = -ex/c1_2
          dexdz = 0.50_double*(f01*ns(i)**one_third)*four_thirds*(opz**one_third - omz**one_third )

          ! correlation
          c20 = 2.0_double*p10*c1*(p30 + c1*(p40 + c1*(p50 + p60*c1)))
          c21 = 2.0_double*p11*c1*(p31 + c1*(p41 + c1*(p51 + p61*c1)))
          c22 = 2.0_double*p12*c1*(p32 + c1*(p42 + c1*(p52 + p62*c1)))
          c30 = log(1.0_double + 1.0_double/c20)
          c31 = log(1.0_double + 1.0_double/c21)
          c32 = log(1.0_double + 1.0_double/c22)
          ec0 = -4.0_double*p10*(1.0_double + p20*c1_2)*c30
          ec1 = -4.0_double*p11*(1.0_double + p21*c1_2)*c31
          ec2 = -4.0_double*p12*(1.0_double + p22*c1_2)*c32
          ec = ec0 + ec2*(f/d2fdz2)*(1.0_double - z4) + (ec1 - ec0)*f*z4
          dc20drs = (p10/c1)*(p30 + c1*(2.0_double*p40 + c1*(3.0_double*p50 + c1*4.0_double*p60)))  
          dc21drs = (p11/c1)*(p31 + c1*(2.0_double*p41 + c1*(3.0_double*p51 + c1*4.0_double*p61)))  
          dc22drs = (p12/c1)*(p32 + c1*(2.0_double*p42 + c1*(3.0_double*p52 + c1*4.0_double*p62)))  
          dec0drs = 4.0_double*p10*(dc20drs*(1.0_double + p20*c1_2)/(c20*(1.0_double + c20)) - p20*c30)  
          dec1drs = 4.0_double*p11*(dc21drs*(1.0_double + p21*c1_2)/(c21*(1.0_double + c21)) - p21*c31)  
          dec2drs = 4.0_double*p12*(dc22drs*(1.0_double + p22*c1_2)/(c22*(1.0_double + c22)) - p22*c32)  
          decdrs = dec0drs + dec2drs*f*(1.0_double - z4)/d2fdz2 + (dec1drs - dec0drs)*f*z4
          decdz = ec2*(dfdz*(1.0_double - z4) + f*(1.0_double - 4.0_double*z3))/d2fdz2 + (ec1 - ec0)*(dfdz*z4 + 4.0_double*f*z3)

          ! combine
          d(i) = ex + ec - (dexdrs + decdrs)*one_third*c1_2 + (dexdz + decdz)*(1.0_double + s*z)

        end do

100     if (error("Exit functional1d_mod::xcd_lsda_pw_1d_i")) continue

      end subroutine

      subroutine xcd_pw91_pz_1d_i(np,n,dn,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = -0.09600_double       ! c1_ca
        real(double), parameter :: p08 = +0.06220_double       ! c2_ca
        real(double), parameter :: p09 = -0.02320_double       ! c3_ca
        real(double), parameter :: p10 = +0.00400_double       ! c4_ca
        real(double), parameter :: p11 = -0.28460_double       ! d1_ca
        real(double), parameter :: p12 = +1.05290_double       ! d2_ca
        real(double), parameter :: p13 = +0.33340_double       ! d3_ca
        real(double), parameter :: p14 = +0.002568_double      ! c0
        real(double), parameter :: p15 = +0.023266_double      ! c1
        real(double), parameter :: p16 = +7.389e-6_double      ! c2
        real(double), parameter :: p17 = +8.723_double         ! c3
        real(double), parameter :: p18 = +0.472_double         ! c4
        real(double), parameter :: p19 = +7.389e-2_double      ! c5
        real(double), parameter :: p20 = -0.001667212_double   ! cx
        real(double), parameter :: p21 = +0.004235_double      ! cc0
        real(double), parameter :: p22 = +0.09_double          ! alfa
        real(double), parameter :: p23 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: rs, rsln, rssq, s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p21*p23                                                 ! beta
        f06 = 2.0_double*p22/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p07 + p09*rs + (p08 + p10*rs)*rsln
            ndecdn = -one_third*(p08 + (p09 + p10)*rs + p10*rs*rsln)
          else
            ec = p11/(1.0_double + p12*rssq + p13*rs)
            ndecdn = one_sixth*p11*(p12*rssq + 2.0_double*p13*rs)/(1.0_double + p12*rssq + p13*rs)**2
          end if
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i)**one_third
          c2 = 100.0_double*t2*f09/n(i)**one_third
          c3 = exp(-c2)
          c4 = -p20 + (p14 + p15*c1 + p16*c1**2)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)
          h1 = 2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*t2*c3
          c5 = 2.0_double*(p14*p18 - p16) + c1*(p15*p18 - p16*p17 + p19*(3.0_double*p14 + c1*(2.0_double*p15 + p16*c1)))
          c5 = c1*(p14*p17 - p15 + c1*c5)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p23*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*c3*(1.0_double - c2)

          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i) = f07*ex/n(i)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcd_pw91_pw_1d_i(np,n,dn,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = + 0.0310907_double    ! a
        real(double), parameter :: p08 = +0.21370_double       ! a1
        real(double), parameter :: p09 = +7.59570_double       ! b1
        real(double), parameter :: p10 = +3.58760_double       ! b2
        real(double), parameter :: p11 = +1.63820_double       ! b3
        real(double), parameter :: p12 = +0.49294_double       ! b4
        real(double), parameter :: p13 = +0.002568_double      ! c0
        real(double), parameter :: p14 = +0.023266_double      ! c1
        real(double), parameter :: p15 = +7.389e-6_double      ! c2
        real(double), parameter :: p16 = +8.723_double         ! c3
        real(double), parameter :: p17 = +0.472_double         ! c4
        real(double), parameter :: p18 = +7.389e-2_double      ! c5
        real(double), parameter :: p19 = -0.001667212_double   ! cx
        real(double), parameter :: p20 = +0.004235_double      ! cc0
        real(double), parameter :: p21 = +0.09_double          ! alfa
        real(double), parameter :: p22 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p20*p22                                                 ! beta
        f06 = 2.0_double*p21/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p07*c1*(p09 + c1*(p10 + c1*(p11 + p12*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p08*c1**2)*c1
          c4 = c4*(p09 + c1*(2.0_double*p10 + c1*(3.0_double*p11 + 4.0_double*p12*c1)))/(c2*(1.0_double + c2))
          ec = -4.0_double*p07*(1.0_double + p08*c1**2)*c3
          ndecdn = four_thirds*p07*(p08*c1**2*c3 + p07*c4)
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i)**one_third
          c2 = 100.0_double*t2*f09/n(i)**one_third
          c3 = exp(-c2)
          c4 = -p19 + (p13 + p14*c1 + p15*c1**2)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)
          h1 = 2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*t2*c3
          c5 = 2.0_double*(p13*p17 - p15) + c1*(p14*p17 - p15*p16 + p18*(3.0_double*p13 + c1*(2.0_double*p14 + p15*c1)))
          c5 = c1*(p13*p16 - p14 + c1*c5)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p22*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*c3*(1.0_double - c2)

          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i) = f07*ex/n(i)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcd_pbe_pz_1d_i(np,n,dn,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.804_double      ! kappa
        real(double), parameter :: p02 = +0.21951_double    ! mu
        real(double), parameter :: p03 = +0.031091_double   ! gamma
        real(double), parameter :: p04 = +0.066725_double   ! beta
        real(double), parameter :: p05 = -0.09600_double    ! c1_ca
        real(double), parameter :: p06 = +0.06220_double    ! c2_ca
        real(double), parameter :: p07 = -0.02320_double    ! c3_ca
        real(double), parameter :: p08 = +0.00400_double    ! c4_ca
        real(double), parameter :: p09 = -0.28460_double    ! d1_ca
        real(double), parameter :: p10 = +1.05290_double    ! d2_ca
        real(double), parameter :: p11 = +0.33340_double    ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: rs, rsln, rssq, s, s2, t, t2, c1, c2, c3, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p05 + p07*rs + (p06 + p08*rs)*rsln
            ndecdn = -one_third*(p06 + (p07 + p08)*rs + p08*rs*rsln)
          else
            ec = p09/(1.0_double + p10*rssq + p11*rs)
            ndecdn = one_sixth*p09*(p10*rssq + 2.0_double*p11*rs)/(1.0_double + p10*rssq + p11*rs)**2
          end if
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c2 = p04/p03/(c1 - 1.0_double)*t2
          c3 = (1.0_double + c2*(1.0_double + c2))**2 + p04/p03*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndhdadadn = -c1*c2**3*ndecdn*(2.0_double + c2)/c3
          dhdtot = p04*(4.0_double + 8.0_double*c2)/c3

          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i) = f05*ex/n(i)**five_thirds*dfxdsos + f06*dhdtot/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcd_pbe_pw_1d_i(np,n,dn,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, d, dd

        integer :: i
        real(double), parameter :: p01 = +0.804_double       ! kappa
        real(double), parameter :: p02 = +0.21951_double     ! mu
        real(double), parameter :: p03 = +0.031091_double    ! gamma
        real(double), parameter :: p04 = +0.066725_double    ! beta
        real(double), parameter :: p05 = +0.0310907_double   ! a
        real(double), parameter :: p06 = +0.21370_double     ! a1
        real(double), parameter :: p07 = +7.59570_double     ! b1
        real(double), parameter :: p08 = +3.58760_double     ! b2
        real(double), parameter :: p09 = +1.63820_double     ! b3
        real(double), parameter :: p10 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: s, s2, t, t2, c1, c2, c3, c4, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p05*c1*(p07 + c1*(p08 + c1*(p09 + p10*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p06*c1**2)*c1*(p07 + c1*(2.0_double*p08 + c1*(3.0_double*p09 + 4.0_double*p10*c1)))/ &
                 (c2*(1.0_double + c2))
          ec = -4.0_double*p05*(1.0_double + p06*c1**2)*c3
          ndecdn = four_thirds*p05*(p06*c1**2*c3 + p05*c4)
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c3 = p04/p03/(c1 - 1.0_double)*t2
          c4 = (1.0_double + c3*(1.0_double + c3))**2 + p04/p03*t2*(1.0_double + c3)*(1.0_double + c3*(1.0_double + c3))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c3)/(1.0_double + c3*(1.0_double + c3)))
          ndhdadadn = -c1*c3**3*ndecdn*(2.0_double + c3)/c4
          dhdtot = p04*(4.0_double + 8.0_double*c3)/c4

          d(i) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i) = f05*ex/n(i)**five_thirds*dfxdsos + f06*dhdtot/n(i)**four_thirds

        end do

      end subroutine

      subroutine xcd_am05_pz_1d_i(np,n,dn,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, d, dd

        integer :: i, j
        real(double), parameter :: p01 = +2.804_double     ! alpha
        real(double), parameter :: p02 = +0.8098_double    ! gamma
        real(double), parameter :: p03 = +0.7168_double    ! c
        real(double), parameter :: p04 = -0.09600_double   ! c1_ca
        real(double), parameter :: p05 = +0.06220_double   ! c2_ca
        real(double), parameter :: p06 = -0.02320_double   ! c3_ca
        real(double), parameter :: p07 = +0.00400_double   ! c4_ca
        real(double), parameter :: p08 = -0.28460_double   ! d1_ca
        real(double), parameter :: p09 = +1.05290_double   ! d2_ca
        real(double), parameter :: p10 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: rs, rsln, rssq, s, s2, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do j = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p04 + p06*rs + (p05 + p07*rs)*rsln
            ndecdn = -one_third*(p05 + (p06 + p07)*rs + p07*rs*rsln)
          else
            ec = p08/(1.0_double + p09*rssq + p10*rs)
            ndecdn = one_sixth*p08*(p09*rssq + 2.0_double*p10*rs)/(1.0_double + p09*rssq + p10*rs)**2
          end if

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          d(i) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i))

        end do

200     call sync_config_process_errors()
        if (error("Exit functional1d_mod::xcd_am05_pz_1d_i")) continue

      end subroutine

      subroutine xcd_am05_pw_1d_i(np,n,dn,d,dd)
        integer :: np
        real(double), dimension(np) :: n, dn, d, dd

        integer :: i, j
        real(double), parameter :: p01 = +2.804_double       ! alpha
        real(double), parameter :: p02 = +0.8098_double      ! gamma
        real(double), parameter :: p03 = +0.7168_double      ! c
        real(double), parameter :: p04 = +0.0310907_double   ! a
        real(double), parameter :: p05 = +0.21370_double     ! a1
        real(double), parameter :: p06 = +7.59570_double     ! b1
        real(double), parameter :: p07 = +3.58760_double     ! b2
        real(double), parameter :: p08 = +1.63820_double     ! b3
        real(double), parameter :: p09 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: s, s2, c1, c2, c3, c4, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do j = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p04*c1*(p06 + c1*(p07 + c1*(p08 + p09*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p04*(1.0_double + p05*c1**2)*c3
          c4 = -(1.0_double + p05*c1**2)*c1*(p06 + c1*(2.0_double*p07 + c1*(3.0_double*p08 + 4.0_double*p09*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p04*(p05*c1**2*c3 + p04*c4)

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          d(i) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i))

        end do

200     call sync_config_process_errors()
        if (error("Exit functional1d_mod::xcd_am05_pw_1d_i")) continue

      end subroutine

      subroutine xcd_blyp_1d_i(np,n,dn,lp,d,dd,ddd)
        integer :: np
        real(double), dimension(np) :: n, dn, lp, d, dd, ddd

        integer :: i
        real(double), parameter :: p01 = +0.27430_double    ! c1
        real(double), parameter :: p02 = +0.19645_double    ! c2
        real(double), parameter :: p03 = +7.79560_double    ! c3
        real(double), parameter :: p04 = +0.04918_double    ! a
        real(double), parameter :: p05 = +0.132_double      ! b
        real(double), parameter :: p06 = +0.2533_double     ! c
        real(double), parameter :: p07 = +0.349_double      ! d
        real(double), parameter :: p08 = +2.871234_double   ! cf
        real(double) :: f01, f02, f03
        real(double) :: one_eighth, one_third, two_thirds, four_thirds, five_thirds, seventeen_thirds
        real(double) :: ex, tb, dfxdsos, dfcdn, dfcddnodn, dfcdd2n
        real(double) :: s, s2, c1, c2, c3, c4, c5, x1, x2, x3

        one_eighth = 1.0_double/8.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        seventeen_thirds = 17.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f03 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f02*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = 1.0_double - p02*p03*s2/x1
          tb = 1.0_double - p01*s2*x3/x2**2
          dfxdsos = p01*(x3/x2 + 1.0_double)/x2

          ! correlation
          c1 = dn(i)/n(i)**four_thirds
          c2 = lp(i)/n(i)**five_thirds
          c3 = n(i)**one_third
          c4 = one_eighth*(c1**2 - c2)
          c5 = (p08*p05 + p05/9.0_double*(c2/2.0_double - 17.0_double*c4))*exp(-p06/c3)
          c4 = seventeen_thirds*(p06*p07 + c3*(p06 - 4.0_double*p07 - 5.0_double*c3))*c1**2
          c4 = c4 - 7.0_double*(p06*p07 + c3*(p06 - p07 - 2.0_double*c3))*c2
          c5 = p05*exp(-p06/c3)
          c4 = c4*c5/24.0_double
          c5 = c5*p08*(p06*p07 + c3*(p06 + 4.0_double*p07 + 3.0_double*c3))
          c4 = c4 - c5 - c3*(4.0_double*p07 + 3.0_double*c3)
          dfcdn = 2.0_double*p04*c4/(3.0_double*(p07 + c3)**2)
          c4 = p05/c3**5/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcddnodn = 34.0_double*p04/36.0_double*c4
          c4 = p05/c3**2/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcdd2n = -14.0_double*p04/24.0_double*c4

          d(i) = four_thirds*ex*tb + dfcdn
          dd(i) = f03*ex/n(i)**five_thirds*dfxdsos + dfcddnodn
          ddd(i) = dfcdd2n

        end do

      end subroutine

      subroutine xce_lda_pz_1d_i(np,n,e)
        integer :: np
        real(double), dimension(np) :: n, e

        integer :: i
        real(double), parameter :: p01 = -0.09600_double   ! c1_ca
        real(double), parameter :: p02 = +0.06220_double   ! c2_ca
        real(double), parameter :: p03 = -0.02320_double   ! c3_ca
        real(double), parameter :: p04 = +0.00400_double   ! c4_ca
        real(double), parameter :: p05 = -0.28460_double   ! d1_ca
        real(double), parameter :: p06 = +1.05290_double   ! d2_ca
        real(double), parameter :: p07 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02
        real(double) :: one_third
        real(double) :: ex, ec
        real(double) :: rs, rsln, rssq

        one_third = 1.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p01 + p03*rs + (p02 + p04*rs)*rsln
          else
            ec = p05/(1.0_double + p06*rssq + p07*rs)
          end if

          e(i) = ex + ec

        end do

      end subroutine

      subroutine xce_lda_pw_1d_i(np,n,e)
        integer :: np
        real(double), dimension(np) :: n, e

        integer :: i
        real(double), parameter :: p01 = +0.0310907_double   ! a
        real(double), parameter :: p02 = +0.21370_double     ! a1
        real(double), parameter :: p03 = +7.59570_double     ! b1
        real(double), parameter :: p04 = +3.58760_double     ! b2
        real(double), parameter :: p05 = +1.63820_double     ! b3
        real(double), parameter :: p06 = +0.49294_double     ! b4
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third
        real(double) :: ex, ec
        real(double) :: c1, c2, c3

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p01*c1*(p03 + c1*(p04 + c1*(p05 + p06*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p01*(1.0_double + p02*c1**2)*c3

          e(i) = ex + ec

        end do

      end subroutine

      subroutine xce_lsda_pw_1d_i(np,ns,zeta,e)
        integer :: np
        real(double), dimension(np) :: ns, zeta, e

        integer :: i
        real(double), parameter :: p10 =  +0.0310907_double   ! a_0
        real(double), parameter :: p20 =  +0.2137000_double   ! a1_0
        real(double), parameter :: p30 =  +7.5957000_double   ! b1_0
        real(double), parameter :: p40 =  +3.5876000_double   ! b2_0
        real(double), parameter :: p50 =  +1.6382000_double   ! b3_0
        real(double), parameter :: p60 =  +0.4929400_double   ! b4_0
        real(double), parameter :: p11 =  +0.0155450_double   ! a_1
        real(double), parameter :: p21 =  +0.2054800_double   ! a1_1
        real(double), parameter :: p31 = +14.1189000_double   ! b1_1
        real(double), parameter :: p41 =  +6.1977000_double   ! b2_1
        real(double), parameter :: p51 =  +3.3662000_double   ! b3_1
        real(double), parameter :: p61 =  +0.6251700_double   ! b4_1
        real(double), parameter :: p12 =  +0.0168870_double   ! a_a
        real(double), parameter :: p22 =  +0.1112500_double   ! a1_a
        real(double), parameter :: p32 = +10.3570000_double   ! b1_a
        real(double), parameter :: p42 =  +3.6231000_double   ! b2_a
        real(double), parameter :: p52 =  +0.8802600_double   ! b3_a
        real(double), parameter :: p62 =  +0.4967100_double   ! b4_a
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: f01, f02, f03
        real(double) :: c1, c1_2, c20, c21, c22, c30, c31, c32
        real(double) :: f, d2fdz2
        real(double) :: ex
        real(double) :: ec, ec0, ec1, ec2
        real(double) :: omz, opz
        real(double) :: z, z4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 2.0_double**four_thirds - 2.0_double

        d2fdz2 = (8.0_double/9.0_double)/f03  ! d2f/dz2 | z=0

        do i = 1,np

          z = zeta(i)
          z4 = z**4

          opz = 1.0_double + z
          omz = 1.0_double - z

          f = (opz**four_thirds + omz**four_thirds - 2.0_double)/f03

          c1 = f02/ns(i)**one_sixth  ! sqrt(rs)
          c1_2 = c1**2               ! rs

          ! exchange
          ex = 0.50_double*(f01*ns(i)**one_third)*(opz**four_thirds + omz**four_thirds )

          ! correlation
          c20 = 2.0_double*p10*c1*(p30 + c1*(p40 + c1*(p50 + p60*c1)))
          c21 = 2.0_double*p11*c1*(p31 + c1*(p41 + c1*(p51 + p61*c1)))
          c22 = 2.0_double*p12*c1*(p32 + c1*(p42 + c1*(p52 + p62*c1)))
          c30 = log(1.0_double + 1.0_double/c20)
          c31 = log(1.0_double + 1.0_double/c21)
          c32 = log(1.0_double + 1.0_double/c22)
          ec0 = -4.0_double*p10*(1.0_double + p20*c1_2)*c30
          ec1 = -4.0_double*p11*(1.0_double + p21*c1_2)*c31
          ec2 = -4.0_double*p12*(1.0_double + p22*c1_2)*c32
          ec = ec0 + ec2*(f/d2fdz2)*(1.0_double - z4) + (ec1 - ec0)*f*z4

          ! combine
          e(i) = ex + ec

        end do

100     if (error("Exit functional1d_mod::xce_lsda_pw_1d_i")) continue

      end subroutine

      subroutine xce_pw91_pz_1d_i(np,n,dn,e)
        integer :: np
        real(double), dimension(np) :: n, dn, e

        integer :: i
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = -0.09600_double       ! c1_ca
        real(double), parameter :: p08 = +0.06220_double       ! c2_ca
        real(double), parameter :: p09 = -0.02320_double       ! c3_ca
        real(double), parameter :: p10 = +0.00400_double       ! c4_ca
        real(double), parameter :: p11 = -0.28460_double       ! d1_ca
        real(double), parameter :: p12 = +1.05290_double       ! d2_ca
        real(double), parameter :: p13 = +0.33340_double       ! d3_ca
        real(double), parameter :: p14 = +0.002568_double      ! c0
        real(double), parameter :: p15 = +0.023266_double      ! c1
        real(double), parameter :: p16 = +7.389e-6_double      ! c2
        real(double), parameter :: p17 = +8.723_double         ! c3
        real(double), parameter :: p18 = +0.472_double         ! c4
        real(double), parameter :: p19 = +7.389e-2_double      ! c5
        real(double), parameter :: p20 = -0.001667212_double   ! cx
        real(double), parameter :: p21 = +0.004235_double      ! cc0
        real(double), parameter :: p22 = +0.09_double          ! alfa
        real(double), parameter :: p23 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f08, f09
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h0, h1
        real(double) :: rs, rsln, rssq, s, s2, s4, t, t2, c1, c2, c3, c4, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p21*p23                                                 ! beta
        f06 = 2.0_double*p22/f05                                      ! aob
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p07 + p09*rs + (p08 + p10*rs)*rsln
          else
            ec = p11/(1.0_double + p12*rssq + p13*rs)
          end if
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          c1 = f08/n(i)**one_third
          c2 = 100.0_double*t2*f09/n(i)**one_third
          c3 = exp(-c2)
          c4 = -p20 + (p14 + p15*c1 + p16*c1**2)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)
          h1 = 2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*t2*c3

          e(i) = ex*fx + ec + h0 + h1

        end do

      end subroutine

      subroutine xce_pw91_pw_1d_i(np,n,dn,e)
        integer :: np
        real(double), dimension(np) :: n, dn, e

        integer :: i
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = + 0.0310907_double    ! a
        real(double), parameter :: p08 = +0.21370_double       ! a1
        real(double), parameter :: p09 = +7.59570_double       ! b1
        real(double), parameter :: p10 = +3.58760_double       ! b2
        real(double), parameter :: p11 = +1.63820_double       ! b3
        real(double), parameter :: p12 = +0.49294_double       ! b4
        real(double), parameter :: p13 = +0.002568_double      ! c0
        real(double), parameter :: p14 = +0.023266_double      ! c1
        real(double), parameter :: p15 = +7.389e-6_double      ! c2
        real(double), parameter :: p16 = +8.723_double         ! c3
        real(double), parameter :: p17 = +0.472_double         ! c4
        real(double), parameter :: p18 = +7.389e-2_double      ! c5
        real(double), parameter :: p19 = -0.001667212_double   ! cx
        real(double), parameter :: p20 = +0.004235_double      ! cc0
        real(double), parameter :: p21 = +0.09_double          ! alfa
        real(double), parameter :: p22 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f08, f09
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h0, h1
        real(double) :: s, s2, s4, t, t2, c1, c2, c3, c4, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p20*p22                                                 ! beta
        f06 = 2.0_double*p21/f05                                      ! aob
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p07*c1*(p09 + c1*(p10 + c1*(p11 + p12*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p07*(1.0_double + p08*c1**2)*c3
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          c1 = f08/n(i)**one_third
          c2 = 100.0_double*t2*f09/n(i)**one_third
          c3 = exp(-c2)
          c4 = -p19 + (p13 + p14*c1 + p15*c1**2)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)
          h1 = 2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*t2*c3

          e(i) = ex*fx + ec + h0 + h1

        end do

      end subroutine

      subroutine xce_pbe_pz_1d_i(np,n,dn,e)
        integer :: np
        real(double), dimension(np) :: n, dn, e

        integer :: i
        real(double), parameter :: p01 = +0.804_double      ! kappa
        real(double), parameter :: p02 = +0.21951_double    ! mu
        real(double), parameter :: p03 = +0.031091_double   ! gamma
        real(double), parameter :: p04 = +0.066725_double   ! beta
        real(double), parameter :: p05 = -0.09600_double    ! c1_ca
        real(double), parameter :: p06 = +0.06220_double    ! c2_ca
        real(double), parameter :: p07 = -0.02320_double    ! c3_ca
        real(double), parameter :: p08 = +0.00400_double    ! c4_ca
        real(double), parameter :: p09 = -0.28460_double    ! d1_ca
        real(double), parameter :: p10 = +1.05290_double    ! d2_ca
        real(double), parameter :: p11 = +0.33340_double    ! d3_ca
        real(double) :: f01, f02, f03, f04
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h
        real(double) :: rs, rsln, rssq, s, s2, t, t2, c1, c2, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p05 + p07*rs + (p06 + p08*rs)*rsln
          else
            ec = p09/(1.0_double + p10*rssq + p11*rs)
          end if
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c2 = p04/p03/(c1 - 1.0_double)*t2
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))

          e(i) = ex*fx + ec + h

        end do

      end subroutine

      subroutine xce_pbe_pw_1d_i(np,n,dn,e)
        integer :: np
        real(double), dimension(np) :: n, dn, e

        integer :: i
        real(double), parameter :: p01 = +0.804_double       ! kappa
        real(double), parameter :: p02 = +0.21951_double     ! mu
        real(double), parameter :: p03 = +0.031091_double    ! gamma
        real(double), parameter :: p04 = +0.066725_double    ! beta
        real(double), parameter :: p05 = +0.0310907_double   ! a
        real(double), parameter :: p06 = +0.21370_double     ! a1
        real(double), parameter :: p07 = +7.59570_double     ! b1
        real(double), parameter :: p08 = +3.58760_double     ! b2
        real(double), parameter :: p09 = +1.63820_double     ! b3
        real(double), parameter :: p10 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h
        real(double) :: s, s2, t, t2, c1, c2, c3, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p05*c1*(p07 + c1*(p08 + c1*(p09 + p10*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p05*(1.0_double + p06*c1**2)*c3
          t = f04*dn(i)/n(i)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c3 = p04/p03/(c1 - 1.0_double)*t2
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c3)/(1.0_double + c3*(1.0_double + c3)))

          e(i) = ex*fx + ec + h

        end do

      end subroutine

      subroutine xce_am05_pz_1d_i(np,n,dn,e)
        integer :: np
        real(double), dimension(np) :: n, dn, e

        integer :: i, j
        real(double), parameter :: p01 = +2.804_double     ! alpha
        real(double), parameter :: p02 = +0.8098_double    ! gamma
        real(double), parameter :: p03 = +0.7168_double    ! c
        real(double), parameter :: p04 = -0.09600_double   ! c1_ca
        real(double), parameter :: p05 = +0.06220_double   ! c2_ca
        real(double), parameter :: p06 = -0.02320_double   ! c3_ca
        real(double), parameter :: p07 = +0.00400_double   ! c4_ca
        real(double), parameter :: p08 = -0.28460_double   ! d1_ca
        real(double), parameter :: p09 = +1.05290_double   ! d2_ca
        real(double), parameter :: p10 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08
        real(double) :: one_sixth, one_quarter, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, xs, hx, hc
        real(double) :: rs, rsln, rssq, s, s2, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do j = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          rsln = one_third*(f02 - log(n(i)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p04 + p06*rs + (p05 + p07*rs)*rsln
          else
            ec = p08/(1.0_double + p09*rssq + p10*rs)
          end if

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02

          e(i) = ex*hx + ec*hc

        end do

200     call sync_config_process_errors()
        if (error("Exit functional1d_mod::xce_am05_pz_1d_i")) continue

      end subroutine

      subroutine xce_am05_pw_1d_i(np,n,dn,e)
        integer :: np
        real(double), dimension(np) :: n, dn, e

        integer :: i, j
        real(double), parameter :: p01 = +2.804_double       ! alpha
        real(double), parameter :: p02 = +0.8098_double      ! gamma
        real(double), parameter :: p03 = +0.7168_double      ! c
        real(double), parameter :: p04 = +0.0310907_double   ! a
        real(double), parameter :: p05 = +0.21370_double     ! a1
        real(double), parameter :: p06 = +7.59570_double     ! b1
        real(double), parameter :: p07 = +3.58760_double     ! b2
        real(double), parameter :: p08 = +1.63820_double     ! b3
        real(double), parameter :: p09 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08
        real(double) :: one_sixth, one_quarter, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, xs, hx, hc
        real(double) :: s, s2, c1, c2, c3, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f03*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do j = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          c1 = f02/n(i)**one_sixth
          c2 = 2.0_double*p04*c1*(p06 + c1*(p07 + c1*(p08 + p09*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p04*(1.0_double + p05*c1**2)*c3

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02

          e(i) = ex*hx + ec*hc

        end do

200     call sync_config_process_errors()
        if (error("Exit functional1d_mod::xce_am05_pw_1d_i")) continue

      end subroutine

      subroutine xce_blyp_1d_i(np,n,dn,lp,e)
        integer :: np
        real(double), dimension(np) :: n, dn, lp, e

        integer :: i
        real(double), parameter :: p01 = +0.27430_double    ! c1
        real(double), parameter :: p02 = +0.19645_double    ! c2
        real(double), parameter :: p03 = +7.79560_double    ! c3
        real(double), parameter :: p04 = +0.04918_double    ! a
        real(double), parameter :: p05 = +0.132_double      ! b
        real(double), parameter :: p06 = +0.2533_double     ! c
        real(double), parameter :: p07 = +0.349_double      ! d
        real(double), parameter :: p08 = +2.871234_double   ! cf
        real(double) :: f01, f02
        real(double) :: one_eighth, one_third, four_thirds, five_thirds
        real(double) :: ex, fx, ec
        real(double) :: s, s2, c1, c2, c3, c4, c5, x1, x2

        one_eighth = 1.0_double/8.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1

        do i = 1,np

          ! exchange
          ex = f01*n(i)**one_third
          s = f02*dn(i)/n(i)**four_thirds
          s2 = s**2
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          fx = 1.0_double + p01*s2/x2

          ! correlation
          c1 = dn(i)/n(i)**four_thirds
          c2 = lp(i)/n(i)**five_thirds
          c3 = n(i)**one_third
          c4 = one_eighth*(c1**2 - c2)
          c5 = (p08*p05 + p05/9.0_double*(c2/2.0_double - 17.0_double*c4))*exp(-p06/c3)
          ec = -2.0_double*p04/(1.0_double + p07/c3)*(1.0_double + c5)

          e(i) = ex*fx + ec

        end do

      end subroutine

      subroutine spin_functions_1d_i(n,ns,zeta)
        real(double), dimension(:) :: n
        real(double), dimension(:), pointer :: ns, zeta

        allocate( ns(size(n)) )
        allocate( zeta(size(n)) )

        ! sum of the spin densities
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,n,ns)

        ! fractional difference of the spin densities
        if (mpi_mysgroup() == 2) n = -1.0_double*n
        call xcomm_pair_allreduce(XSGROUP,MPI_SUM,n,zeta)
        zeta = zeta/ns
        if (mpi_mysgroup() == 2) n = -1.0_double*n

      end subroutine

      subroutine xcq_lda_pz_3d_i(np,n,e,d)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, e, d

        integer :: i1, i2, i3
        real(double), parameter :: p01 = -0.09600_double   ! c1_ca
        real(double), parameter :: p02 = +0.06220_double   ! c2_ca
        real(double), parameter :: p03 = -0.02320_double   ! c3_ca
        real(double), parameter :: p04 = +0.00400_double   ! c4_ca
        real(double), parameter :: p05 = -0.28460_double   ! d1_ca
        real(double), parameter :: p06 = +1.05290_double   ! d2_ca
        real(double), parameter :: p07 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: rs, rsln, rssq

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p01 + p03*rs + (p02 + p04*rs)*rsln
            ndecdn = -one_third*(p02 + (p03 + p04)*rs + p04*rs*rsln)
          else
            ec = p05/(1.0_double + p06*rssq + p07*rs)
            ndecdn = one_sixth*p05*(p06*rssq + 2.0_double*p07*rs)/(1.0_double + p06*rssq + p07*rs)**2
          end if

          e(i1,i2,i3) = ex + ec
          d(i1,i2,i3) = four_thirds*ex + ec + ndecdn

        end do
        end do
        end do

      end subroutine

      subroutine xcq_lda_pw_3d_i(np,n,e,d)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, e, d

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.0310907_double   ! a
        real(double), parameter :: p02 = +0.21370_double     ! a1
        real(double), parameter :: p03 = +7.59570_double     ! b1
        real(double), parameter :: p04 = +3.58760_double     ! b2
        real(double), parameter :: p05 = +1.63820_double     ! b3
        real(double), parameter :: p06 = +0.49294_double     ! b4
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: c1, c2, c3, c4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p01*c1*(p03 + c1*(p04 + c1*(p05 + p06*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p01*(1.0_double + p02*c1**2)*c3
          c4 = -(1.0_double + p02*c1**2)*c1*(p03 + c1*(2.0_double*p04 + c1*(3.0_double*p05 + 4.0_double*p06*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p01*(p02*c1**2*c3 + p01*c4)

          e(i1,i2,i3) = ex + ec
          d(i1,i2,i3) = four_thirds*ex + ec + ndecdn

        end do
        end do
        end do

      end subroutine

      subroutine xcq_lsda_pw_3d_i(np,ns,zeta,e,d)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: ns, zeta, e, d

        integer :: i1, i2, i3
        real(double), parameter :: p10 =  +0.0310907_double   ! a_0
        real(double), parameter :: p20 =  +0.2137000_double   ! a1_0
        real(double), parameter :: p30 =  +7.5957000_double   ! b1_0
        real(double), parameter :: p40 =  +3.5876000_double   ! b2_0
        real(double), parameter :: p50 =  +1.6382000_double   ! b3_0
        real(double), parameter :: p60 =  +0.4929400_double   ! b4_0
        real(double), parameter :: p11 =  +0.0155450_double   ! a_1
        real(double), parameter :: p21 =  +0.2054800_double   ! a1_1
        real(double), parameter :: p31 = +14.1189000_double   ! b1_1
        real(double), parameter :: p41 =  +6.1977000_double   ! b2_1
        real(double), parameter :: p51 =  +3.3662000_double   ! b3_1
        real(double), parameter :: p61 =  +0.6251700_double   ! b4_1
        real(double), parameter :: p12 =  +0.0168870_double   ! a_a
        real(double), parameter :: p22 =  +0.1112500_double   ! a1_a
        real(double), parameter :: p32 = +10.3570000_double   ! b1_a
        real(double), parameter :: p42 =  +3.6231000_double   ! b2_a
        real(double), parameter :: p52 =  +0.8802600_double   ! b3_a
        real(double), parameter :: p62 =  +0.4967100_double   ! b4_a
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: f01, f02, f03
        real(double) :: c1, c1_2, c20, c21, c22, c30, c31, c32
        real(double) :: dc20drs, dc21drs, dc22drs
        real(double) :: f, dfdz, d2fdz2
        real(double) :: ex, dexdrs, dexdz
        real(double) :: ec, decdrs, decdz
        real(double) :: ec0, ec1, ec2, dec0drs, dec1drs, dec2drs
        real(double) :: omz, opz, s
        real(double) :: z, z3, z4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 2.0_double**four_thirds - 2.0_double

        d2fdz2 = (8.0_double/9.0_double)/f03  ! d2f/dz2 | z=0

        if (mpi_mysgroup() == 1) s = -1.0_double
        if (mpi_mysgroup() == 2) s = +1.0_double

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          z = zeta(i1,i2,i3)
          z3 = z**3
          z4 = z**4

          opz = 1.0_double + z
          omz = 1.0_double - z

          f = (opz**four_thirds + omz**four_thirds - 2.0_double)/f03
          dfdz =  four_thirds*(opz**one_third - omz**one_third)/f03

          c1 = f02/ns(i1,i2,i3)**one_sixth  ! sqrt(rs)
          c1_2 = c1**2                      ! rs

          ! exchange
          ex = 0.50_double*(f01*ns(i1,i2,i3)**one_third)*(opz**four_thirds + omz**four_thirds )
          dexdrs = -ex/c1_2
          dexdz = 0.50_double*(f01*ns(i1,i2,i3)**one_third)*four_thirds*(opz**one_third - omz**one_third )

          ! correlation
          c20 = 2.0_double*p10*c1*(p30 + c1*(p40 + c1*(p50 + p60*c1)))
          c21 = 2.0_double*p11*c1*(p31 + c1*(p41 + c1*(p51 + p61*c1)))
          c22 = 2.0_double*p12*c1*(p32 + c1*(p42 + c1*(p52 + p62*c1)))
          c30 = log(1.0_double + 1.0_double/c20)
          c31 = log(1.0_double + 1.0_double/c21)
          c32 = log(1.0_double + 1.0_double/c22)
          ec0 = -4.0_double*p10*(1.0_double + p20*c1_2)*c30
          ec1 = -4.0_double*p11*(1.0_double + p21*c1_2)*c31
          ec2 = -4.0_double*p12*(1.0_double + p22*c1_2)*c32
          ec = ec0 + ec2*(f/d2fdz2)*(1.0_double - z4) + (ec1 - ec0)*f*z4
          dc20drs = (p10/c1)*(p30 + c1*(2.0_double*p40 + c1*(3.0_double*p50 + c1*4.0_double*p60)))
          dc21drs = (p11/c1)*(p31 + c1*(2.0_double*p41 + c1*(3.0_double*p51 + c1*4.0_double*p61)))
          dc22drs = (p12/c1)*(p32 + c1*(2.0_double*p42 + c1*(3.0_double*p52 + c1*4.0_double*p62)))
          dec0drs = 4.0_double*p10*(dc20drs*(1.0_double + p20*c1_2)/(c20*(1.0_double + c20)) - p20*c30)
          dec1drs = 4.0_double*p11*(dc21drs*(1.0_double + p21*c1_2)/(c21*(1.0_double + c21)) - p21*c31)
          dec2drs = 4.0_double*p12*(dc22drs*(1.0_double + p22*c1_2)/(c22*(1.0_double + c22)) - p22*c32)
          decdrs = dec0drs + dec2drs*f*(1.0_double - z4)/d2fdz2 + (dec1drs - dec0drs)*f*z4
          decdz = ec2*(dfdz*(1.0_double - z4) + f*(1.0_double - 4.0_double*z3))/d2fdz2 + (ec1 - ec0)*(dfdz*z4 + 4.0_double*f*z3)

          ! combine
          e(i1,i2,i3) = ex + ec
          d(i1,i2,i3) = ex + ec - (dexdrs + decdrs)*one_third*c1_2 + (dexdz + decdz)*(1.0_double + s*z)

        end do
        end do
        end do

100     if (error("Exit functional_mod::xcq_lsda_pw_3d_i")) continue

      end subroutine

      subroutine xcq_pw91_pz_3d_i(np,n,dn,e,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = -0.09600_double       ! c1_ca
        real(double), parameter :: p08 = +0.06220_double       ! c2_ca
        real(double), parameter :: p09 = -0.02320_double       ! c3_ca
        real(double), parameter :: p10 = +0.00400_double       ! c4_ca
        real(double), parameter :: p11 = -0.28460_double       ! d1_ca
        real(double), parameter :: p12 = +1.05290_double       ! d2_ca
        real(double), parameter :: p13 = +0.33340_double       ! d3_ca
        real(double), parameter :: p14 = +0.002568_double      ! c0
        real(double), parameter :: p15 = +0.023266_double      ! c1
        real(double), parameter :: p16 = +7.389e-6_double      ! c2
        real(double), parameter :: p17 = +8.723_double         ! c3
        real(double), parameter :: p18 = +0.472_double         ! c4
        real(double), parameter :: p19 = +7.389e-2_double      ! c5
        real(double), parameter :: p20 = -0.001667212_double   ! cx
        real(double), parameter :: p21 = +0.004235_double      ! cc0
        real(double), parameter :: p22 = +0.09_double          ! alfa
        real(double), parameter :: p23 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: rs, rsln, rssq, s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p21*p23                                                 ! beta
        f06 = 2.0_double*p22/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p07 + p09*rs + (p08 + p10*rs)*rsln
            ndecdn = -one_third*(p08 + (p09 + p10)*rs + p10*rs*rsln)
          else
            ec = p11/(1.0_double + p12*rssq + p13*rs)
            ndecdn = one_sixth*p11*(p12*rssq + 2.0_double*p13*rs)/(1.0_double + p12*rssq + p13*rs)**2
          end if
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i1,i2,i3)**one_third
          c2 = 100.0_double*t2*f09/n(i1,i2,i3)**one_third
          c3 = exp(-c2)
          c4 = -p20 + (p14 + p15*c1 + p16*c1**2)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)
          h1 = 2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*t2*c3
          c5 = 2.0_double*(p14*p18 - p16) + c1*(p15*p18 - p16*p17 + p19*(3.0_double*p14 + c1*(2.0_double*p15 + p16*c1)))
          c5 = c1*(p14*p17 - p15 + c1*c5)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p23*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*c3*(1.0_double - c2)

          e(i1,i2,i3) = ex*fx + ec + h0 + h1
          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i1,i2,i3) = f07*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcq_pw91_pw_3d_i(np,n,dn,e,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = + 0.0310907_double    ! a
        real(double), parameter :: p08 = +0.21370_double       ! a1
        real(double), parameter :: p09 = +7.59570_double       ! b1
        real(double), parameter :: p10 = +3.58760_double       ! b2
        real(double), parameter :: p11 = +1.63820_double       ! b3
        real(double), parameter :: p12 = +0.49294_double       ! b4
        real(double), parameter :: p13 = +0.002568_double      ! c0
        real(double), parameter :: p14 = +0.023266_double      ! c1
        real(double), parameter :: p15 = +7.389e-6_double      ! c2
        real(double), parameter :: p16 = +8.723_double         ! c3
        real(double), parameter :: p17 = +0.472_double         ! c4
        real(double), parameter :: p18 = +7.389e-2_double      ! c5
        real(double), parameter :: p19 = -0.001667212_double   ! cx
        real(double), parameter :: p20 = +0.004235_double      ! cc0
        real(double), parameter :: p21 = +0.09_double          ! alfa
        real(double), parameter :: p22 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p20*p22                                                 ! beta
        f06 = 2.0_double*p21/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p07*c1*(p09 + c1*(p10 + c1*(p11 + p12*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p08*c1**2)*c1
          c4 = c4*(p09 + c1*(2.0_double*p10 + c1*(3.0_double*p11 + 4.0_double*p12*c1)))/(c2*(1.0_double + c2))
          ec = -4.0_double*p07*(1.0_double + p08*c1**2)*c3
          ndecdn = four_thirds*p07*(p08*c1**2*c3 + p07*c4)
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i1,i2,i3)**one_third
          c2 = 100.0_double*t2*f09/n(i1,i2,i3)**one_third
          c3 = exp(-c2)
          c4 = -p19 + (p13 + p14*c1 + p15*c1**2)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)
          h1 = 2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*t2*c3
          c5 = 2.0_double*(p13*p17 - p15) + c1*(p14*p17 - p15*p16 + p18*(3.0_double*p13 + c1*(2.0_double*p14 + p15*c1)))
          c5 = c1*(p13*p16 - p14 + c1*c5)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p22*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*c3*(1.0_double - c2)

          e(i1,i2,i3) = ex*fx + ec + h0 + h1
          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i1,i2,i3) = f07*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcq_pbe_pz_3d_i(np,n,dn,e,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.804_double      ! kappa
        real(double), parameter :: p02 = +0.21951_double    ! mu
        real(double), parameter :: p03 = +0.031091_double   ! gamma
        real(double), parameter :: p04 = +0.066725_double   ! beta
        real(double), parameter :: p05 = -0.09600_double    ! c1_ca
        real(double), parameter :: p06 = +0.06220_double    ! c2_ca
        real(double), parameter :: p07 = -0.02320_double    ! c3_ca
        real(double), parameter :: p08 = +0.00400_double    ! c4_ca
        real(double), parameter :: p09 = -0.28460_double    ! d1_ca
        real(double), parameter :: p10 = +1.05290_double    ! d2_ca
        real(double), parameter :: p11 = +0.33340_double    ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: rs, rsln, rssq, s, s2, t, t2, c1, c2, c3, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p05 + p07*rs + (p06 + p08*rs)*rsln
            ndecdn = -one_third*(p06 + (p07 + p08)*rs + p08*rs*rsln)
          else
            ec = p09/(1.0_double + p10*rssq + p11*rs)
            ndecdn = one_sixth*p09*(p10*rssq + 2.0_double*p11*rs)/(1.0_double + p10*rssq + p11*rs)**2
          end if
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c2 = p04/p03/(c1 - 1.0_double)*t2
          c3 = (1.0_double + c2*(1.0_double + c2))**2 + p04/p03*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndhdadadn = -c1*c2**3*ndecdn*(2.0_double + c2)/c3
          dhdtot = p04*(4.0_double + 8.0_double*c2)/c3

          e(i1,i2,i3) = ex*fx + ec + h
          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i1,i2,i3) = f05*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f06*dhdtot/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcq_pbe_pw_3d_i(np,n,dn,e,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.804_double       ! kappa
        real(double), parameter :: p02 = +0.21951_double     ! mu
        real(double), parameter :: p03 = +0.031091_double    ! gamma
        real(double), parameter :: p04 = +0.066725_double    ! beta
        real(double), parameter :: p05 = +0.0310907_double   ! a
        real(double), parameter :: p06 = +0.21370_double     ! a1
        real(double), parameter :: p07 = +7.59570_double     ! b1
        real(double), parameter :: p08 = +3.58760_double     ! b2
        real(double), parameter :: p09 = +1.63820_double     ! b3
        real(double), parameter :: p10 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: s, s2, t, t2, c1, c2, c3, c4, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p05*c1*(p07 + c1*(p08 + c1*(p09 + p10*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p06*c1**2)*c1*(p07 + c1*(2.0_double*p08 + c1*(3.0_double*p09 + 4.0_double*p10*c1)))/ &
                 (c2*(1.0_double + c2))
          ec = -4.0_double*p05*(1.0_double + p06*c1**2)*c3
          ndecdn = four_thirds*p05*(p06*c1**2*c3 + p05*c4)
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c3 = p04/p03/(c1 - 1.0_double)*t2
          c4 = (1.0_double + c3*(1.0_double + c3))**2 + p04/p03*t2*(1.0_double + c3)*(1.0_double + c3*(1.0_double + c3))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c3)/(1.0_double + c3*(1.0_double + c3)))
          ndhdadadn = -c1*c3**3*ndecdn*(2.0_double + c3)/c4
          dhdtot = p04*(4.0_double + 8.0_double*c3)/c4

          e(i1,i2,i3) = ex*fx + ec + h
          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i1,i2,i3) = f05*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f06*dhdtot/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcq_am05_pz_3d_i(np,n,dn,e,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e, d, dd

        integer :: i, i1, i2, i3
        real(double), parameter :: p01 = +2.804_double     ! alpha
        real(double), parameter :: p02 = +0.8098_double    ! gamma
        real(double), parameter :: p03 = +0.7168_double    ! c
        real(double), parameter :: p04 = -0.09600_double   ! c1_ca
        real(double), parameter :: p05 = +0.06220_double   ! c2_ca
        real(double), parameter :: p06 = -0.02320_double   ! c3_ca
        real(double), parameter :: p07 = +0.00400_double   ! c4_ca
        real(double), parameter :: p08 = -0.28460_double   ! d1_ca
        real(double), parameter :: p09 = +1.05290_double   ! d2_ca
        real(double), parameter :: p10 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: rs, rsln, rssq, s, s2, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do i = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p04 + p06*rs + (p05 + p07*rs)*rsln
            ndecdn = -one_third*(p05 + (p06 + p07)*rs + p07*rs*rsln)
          else
            ec = p08/(1.0_double + p09*rssq + p10*rs)
            ndecdn = one_sixth*p08*(p09*rssq + 2.0_double*p10*rs)/(1.0_double + p09*rssq + p10*rs)**2
          end if

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i1,i2,i3)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          e(i1,i2,i3) = ex*hx + ec*hc
          d(i1,i2,i3) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i1,i2,i3) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i1,i2,i3))

        end do
        end do
        end do

200     call sync_config_process_errors()
        if (error("Exit functional_mod::xcq_am05_pz_3d_i")) continue

      end subroutine

      subroutine xcq_am05_pw_3d_i(np,n,dn,e,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e, d, dd

        integer :: i, i1, i2, i3
        real(double), parameter :: p01 = +2.804_double       ! alpha
        real(double), parameter :: p02 = +0.8098_double      ! gamma
        real(double), parameter :: p03 = +0.7168_double      ! c
        real(double), parameter :: p04 = +0.0310907_double   ! a
        real(double), parameter :: p05 = +0.21370_double     ! a1
        real(double), parameter :: p06 = +7.59570_double     ! b1
        real(double), parameter :: p07 = +3.58760_double     ! b2
        real(double), parameter :: p08 = +1.63820_double     ! b3
        real(double), parameter :: p09 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: s, s2, c1, c2, c3, c4, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do i = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p04*c1*(p06 + c1*(p07 + c1*(p08 + p09*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p04*(1.0_double + p05*c1**2)*c3
          c4 = -(1.0_double + p05*c1**2)*c1*(p06 + c1*(2.0_double*p07 + c1*(3.0_double*p08 + 4.0_double*p09*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p04*(p05*c1**2*c3 + p04*c4)

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i1,i2,i3)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          e(i1,i2,i3) = ex*hx + ec*hc
          d(i1,i2,i3) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i1,i2,i3) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i1,i2,i3))

        end do
        end do
        end do

200     call sync_config_process_errors()
        if (error("Exit functional_mod::xcq_am05_pw_3d_i")) continue

      end subroutine

      subroutine xcq_blyp_3d_i(np,n,dn,lp,e,d,dd,ddd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, lp, e, d, dd, ddd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double    ! c1
        real(double), parameter :: p02 = +0.19645_double    ! c2
        real(double), parameter :: p03 = +7.79560_double    ! c3
        real(double), parameter :: p04 = +0.04918_double    ! a
        real(double), parameter :: p05 = +0.132_double      ! b
        real(double), parameter :: p06 = +0.2533_double     ! c
        real(double), parameter :: p07 = +0.349_double      ! d
        real(double), parameter :: p08 = +2.871234_double   ! cf
        real(double) :: f01, f02, f03
        real(double) :: one_eighth, one_third, two_thirds, four_thirds, five_thirds, seventeen_thirds
        real(double) :: ex, fx, tb, dfxdsos, ec, dfcdn, dfcddnodn, dfcdd2n
        real(double) :: s, s2, c1, c2, c3, c4, c5, x1, x2, x3

        one_eighth = 1.0_double/8.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        seventeen_thirds = 17.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f03 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f02*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          fx = 1.0_double + p01*s2/x2
          x3 = 1.0_double - p02*p03*s2/x1
          tb = 1.0_double - p01*s2*x3/x2**2
          dfxdsos = p01*(x3/x2 + 1.0_double)/x2

          ! correlation
          c1 = dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          c2 = lp(i1,i2,i3)/n(i1,i2,i3)**five_thirds
          c3 = n(i1,i2,i3)**one_third
          c4 = one_eighth*(c1**2 - c2)
          c5 = (p08*p05 + p05/9.0_double*(c2/2.0_double - 17.0_double*c4))*exp(-p06/c3)
          ec = -2.0_double*p04/(1.0_double + p07/c3)*(1.0_double + c5)
          c4 = seventeen_thirds*(p06*p07 + c3*(p06 - 4.0_double*p07 - 5.0_double*c3))*c1**2
          c4 = c4 - 7.0_double*(p06*p07 + c3*(p06 - p07 - 2.0_double*c3))*c2
          c5 = p05*exp(-p06/c3)
          c4 = c4*c5/24.0_double
          c5 = c5*p08*(p06*p07 + c3*(p06 + 4.0_double*p07 + 3.0_double*c3))
          c4 = c4 - c5 - c3*(4.0_double*p07 + 3.0_double*c3)
          dfcdn = 2.0_double*p04*c4/(3.0_double*(p07 + c3)**2)
          c4 = p05/c3**5/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcddnodn = 34.0_double*p04/36.0_double*c4
          c4 = p05/c3**2/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcdd2n = -14.0_double*p04/24.0_double*c4

          e(i1,i2,i3) = ex*fx + ec
          d(i1,i2,i3) = four_thirds*ex*tb + dfcdn
          dd(i1,i2,i3) = f03*ex/n(i1,i2,i3)**five_thirds*dfxdsos + dfcddnodn
          ddd(i1,i2,i3) = dfcdd2n

        end do
        end do
        end do

      end subroutine

      subroutine xcd_lda_pz_3d_i(np,n,d)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, d

        integer :: i1, i2, i3
        real(double), parameter :: p01 = -0.09600_double   ! c1_ca
        real(double), parameter :: p02 = +0.06220_double   ! c2_ca
        real(double), parameter :: p03 = -0.02320_double   ! c3_ca
        real(double), parameter :: p04 = +0.00400_double   ! c4_ca
        real(double), parameter :: p05 = -0.28460_double   ! d1_ca
        real(double), parameter :: p06 = +1.05290_double   ! d2_ca
        real(double), parameter :: p07 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: rs, rsln, rssq

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p01 + p03*rs + (p02 + p04*rs)*rsln
            ndecdn = -one_third*(p02 + (p03 + p04)*rs + p04*rs*rsln)
          else
            ec = p05/(1.0_double + p06*rssq + p07*rs)
            ndecdn = one_sixth*p05*(p06*rssq + 2.0_double*p07*rs)/(1.0_double + p06*rssq + p07*rs)**2
          end if

          d(i1,i2,i3) = four_thirds*ex + ec + ndecdn

        end do
        end do
        end do

      end subroutine

      subroutine xcd_lda_pw_3d_i(np,n,d)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, d

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.0310907_double   ! a
        real(double), parameter :: p02 = +0.21370_double     ! a1
        real(double), parameter :: p03 = +7.59570_double     ! b1
        real(double), parameter :: p04 = +3.58760_double     ! b2
        real(double), parameter :: p05 = +1.63820_double     ! b3
        real(double), parameter :: p06 = +0.49294_double     ! b4
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: ex, ec, ndecdn
        real(double) :: c1, c2, c3, c4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p01*c1*(p03 + c1*(p04 + c1*(p05 + p06*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p01*(1.0_double + p02*c1**2)*c3
          c4 = -(1.0_double + p02*c1**2)*c1*(p03 + c1*(2.0_double*p04 + c1*(3.0_double*p05 + 4.0_double*p06*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p01*(p02*c1**2*c3 + p01*c4)

          d(i1,i2,i3) = four_thirds*ex + ec + ndecdn

        end do
        end do
        end do

      end subroutine

      subroutine xcd_lsda_pw_3d_i(np,ns,zeta,d)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: ns, zeta, d

        integer :: i1, i2, i3
        real(double), parameter :: p10 =  +0.0310907_double   ! a_0
        real(double), parameter :: p20 =  +0.2137000_double   ! a1_0
        real(double), parameter :: p30 =  +7.5957000_double   ! b1_0
        real(double), parameter :: p40 =  +3.5876000_double   ! b2_0
        real(double), parameter :: p50 =  +1.6382000_double   ! b3_0
        real(double), parameter :: p60 =  +0.4929400_double   ! b4_0
        real(double), parameter :: p11 =  +0.0155450_double   ! a_1
        real(double), parameter :: p21 =  +0.2054800_double   ! a1_1
        real(double), parameter :: p31 = +14.1189000_double   ! b1_1
        real(double), parameter :: p41 =  +6.1977000_double   ! b2_1
        real(double), parameter :: p51 =  +3.3662000_double   ! b3_1
        real(double), parameter :: p61 =  +0.6251700_double   ! b4_1
        real(double), parameter :: p12 =  +0.0168870_double   ! a_a
        real(double), parameter :: p22 =  +0.1112500_double   ! a1_a
        real(double), parameter :: p32 = +10.3570000_double   ! b1_a
        real(double), parameter :: p42 =  +3.6231000_double   ! b2_a
        real(double), parameter :: p52 =  +0.8802600_double   ! b3_a
        real(double), parameter :: p62 =  +0.4967100_double   ! b4_a
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: f01, f02, f03
        real(double) :: c1, c1_2, c20, c21, c22, c30, c31, c32
        real(double) :: dc20drs, dc21drs, dc22drs
        real(double) :: f, dfdz, d2fdz2
        real(double) :: ex, dexdrs, dexdz
        real(double) :: ec, decdrs, decdz
        real(double) :: ec0, ec1, ec2, dec0drs, dec1drs, dec2drs
        real(double) :: omz, opz, s
        real(double) :: z, z3, z4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 2.0_double**four_thirds - 2.0_double

        d2fdz2 = (8.0_double/9.0_double)/f03  ! d2f/dz2 | z=0

        if (mpi_mysgroup() == 1) s = -1.0_double
        if (mpi_mysgroup() == 2) s = +1.0_double

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          z = zeta(i1,i2,i3)
          z3 = z**3
          z4 = z**4

          opz = 1.0_double + z
          omz = 1.0_double - z

          f = (opz**four_thirds + omz**four_thirds - 2.0_double)/f03
          dfdz =  four_thirds*(opz**one_third - omz**one_third)/f03

          c1 = f02/ns(i1,i2,i3)**one_sixth  ! sqrt(rs)
          c1_2 = c1**2                      ! rs

          ! exchange
          ex = 0.50_double*(f01*ns(i1,i2,i3)**one_third)*(opz**four_thirds + omz**four_thirds )
          dexdrs = -ex/c1_2
          dexdz = 0.50_double*(f01*ns(i1,i2,i3)**one_third)*four_thirds*(opz**one_third - omz**one_third )

          ! correlation
          c20 = 2.0_double*p10*c1*(p30 + c1*(p40 + c1*(p50 + p60*c1)))
          c21 = 2.0_double*p11*c1*(p31 + c1*(p41 + c1*(p51 + p61*c1)))
          c22 = 2.0_double*p12*c1*(p32 + c1*(p42 + c1*(p52 + p62*c1)))
          c30 = log(1.0_double + 1.0_double/c20)
          c31 = log(1.0_double + 1.0_double/c21)
          c32 = log(1.0_double + 1.0_double/c22)
          ec0 = -4.0_double*p10*(1.0_double + p20*c1_2)*c30
          ec1 = -4.0_double*p11*(1.0_double + p21*c1_2)*c31
          ec2 = -4.0_double*p12*(1.0_double + p22*c1_2)*c32
          ec = ec0 + ec2*(f/d2fdz2)*(1.0_double - z4) + (ec1 - ec0)*f*z4
          dc20drs = (p10/c1)*(p30 + c1*(2.0_double*p40 + c1*(3.0_double*p50 + c1*4.0_double*p60)))
          dc21drs = (p11/c1)*(p31 + c1*(2.0_double*p41 + c1*(3.0_double*p51 + c1*4.0_double*p61)))
          dc22drs = (p12/c1)*(p32 + c1*(2.0_double*p42 + c1*(3.0_double*p52 + c1*4.0_double*p62)))
          dec0drs = 4.0_double*p10*(dc20drs*(1.0_double + p20*c1_2)/(c20*(1.0_double + c20)) - p20*c30)
          dec1drs = 4.0_double*p11*(dc21drs*(1.0_double + p21*c1_2)/(c21*(1.0_double + c21)) - p21*c31)
          dec2drs = 4.0_double*p12*(dc22drs*(1.0_double + p22*c1_2)/(c22*(1.0_double + c22)) - p22*c32)
          decdrs = dec0drs + dec2drs*f*(1.0_double - z4)/d2fdz2 + (dec1drs - dec0drs)*f*z4
          decdz = ec2*(dfdz*(1.0_double - z4) + f*(1.0_double - 4.0_double*z3))/d2fdz2 + (ec1 - ec0)*(dfdz*z4 + 4.0_double*f*z3)

          ! combine
          d(i1,i2,i3) = ex + ec - (dexdrs + decdrs)*one_third*c1_2 + (dexdz + decdz)*(1.0_double + s*z)

        end do
        end do
        end do

100     if (error("Exit functional_mod::xcd_lsda_pw_3d_i")) continue

      end subroutine

      subroutine xcd_pw91_pz_3d_i(np,n,dn,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = -0.09600_double       ! c1_ca
        real(double), parameter :: p08 = +0.06220_double       ! c2_ca
        real(double), parameter :: p09 = -0.02320_double       ! c3_ca
        real(double), parameter :: p10 = +0.00400_double       ! c4_ca
        real(double), parameter :: p11 = -0.28460_double       ! d1_ca
        real(double), parameter :: p12 = +1.05290_double       ! d2_ca
        real(double), parameter :: p13 = +0.33340_double       ! d3_ca
        real(double), parameter :: p14 = +0.002568_double      ! c0
        real(double), parameter :: p15 = +0.023266_double      ! c1
        real(double), parameter :: p16 = +7.389e-6_double      ! c2
        real(double), parameter :: p17 = +8.723_double         ! c3
        real(double), parameter :: p18 = +0.472_double         ! c4
        real(double), parameter :: p19 = +7.389e-2_double      ! c5
        real(double), parameter :: p20 = -0.001667212_double   ! cx
        real(double), parameter :: p21 = +0.004235_double      ! cc0
        real(double), parameter :: p22 = +0.09_double          ! alfa
        real(double), parameter :: p23 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: rs, rsln, rssq, s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p21*p23                                                 ! beta
        f06 = 2.0_double*p22/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p07 + p09*rs + (p08 + p10*rs)*rsln
            ndecdn = -one_third*(p08 + (p09 + p10)*rs + p10*rs*rsln)
          else
            ec = p11/(1.0_double + p12*rssq + p13*rs)
            ndecdn = one_sixth*p11*(p12*rssq + 2.0_double*p13*rs)/(1.0_double + p12*rssq + p13*rs)**2
          end if
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i1,i2,i3)**one_third
          c2 = 100.0_double*t2*f09/n(i1,i2,i3)**one_third
          c3 = exp(-c2)
          c4 = -p20 + (p14 + p15*c1 + p16*c1**2)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)
          h1 = 2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*t2*c3
          c5 = 2.0_double*(p14*p18 - p16) + c1*(p15*p18 - p16*p17 + p19*(3.0_double*p14 + c1*(2.0_double*p15 + p16*c1)))
          c5 = c1*(p14*p17 - p15 + c1*c5)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p23*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*c3*(1.0_double - c2)

          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i1,i2,i3) = f07*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcd_pw91_pw_3d_i(np,n,dn,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = + 0.0310907_double    ! a
        real(double), parameter :: p08 = +0.21370_double       ! a1
        real(double), parameter :: p09 = +7.59570_double       ! b1
        real(double), parameter :: p10 = +3.58760_double       ! b2
        real(double), parameter :: p11 = +1.63820_double       ! b3
        real(double), parameter :: p12 = +0.49294_double       ! b4
        real(double), parameter :: p13 = +0.002568_double      ! c0
        real(double), parameter :: p14 = +0.023266_double      ! c1
        real(double), parameter :: p15 = +7.389e-6_double      ! c2
        real(double), parameter :: p16 = +8.723_double         ! c3
        real(double), parameter :: p17 = +0.472_double         ! c4
        real(double), parameter :: p18 = +7.389e-2_double      ! c5
        real(double), parameter :: p19 = -0.001667212_double   ! cx
        real(double), parameter :: p20 = +0.004235_double      ! cc0
        real(double), parameter :: p21 = +0.09_double          ! alfa
        real(double), parameter :: p22 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09, f10
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h0, h1, ndh0dadadn, ndh1drsdrsdn, dh0dtot, dh1dtot
        real(double) :: s, s2, s4, t, t2, c1, c2, c3, c4, c5, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p20*p22                                                 ! beta
        f06 = 2.0_double*p21/f05                                      ! aob
        f07 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third
        f10 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4
          x2 = x3 - p06*s2
          x3 = -p02*p03*x2*s2/x1 + x2*(1.0_double - 3.0_double*p06*s4)
          dfxdsos = (x3/x4 + x2 - 2.0_double*s2*(p06 - p05*(p01 - p06*s2 - x2)))/x4

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p07*c1*(p09 + c1*(p10 + c1*(p11 + p12*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p08*c1**2)*c1
          c4 = c4*(p09 + c1*(2.0_double*p10 + c1*(3.0_double*p11 + 4.0_double*p12*c1)))/(c2*(1.0_double + c2))
          ec = -4.0_double*p07*(1.0_double + p08*c1**2)*c3
          ndecdn = four_thirds*p07*(p08*c1**2*c3 + p07*c4)
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          c3 = ((1.0_double + c2*(1.0_double + c2))**2 + f06*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2)))
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndh0dadadn = -c2**3*c1*ndecdn*(2.0_double + c2)/c3
          dh0dtot = 2.0_double*f05*(2.0_double + 4.0_double*c2)/c3
          c1 = f08/n(i1,i2,i3)**one_third
          c2 = 100.0_double*t2*f09/n(i1,i2,i3)**one_third
          c3 = exp(-c2)
          c4 = -p19 + (p13 + p14*c1 + p15*c1**2)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)
          h1 = 2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*t2*c3
          c5 = 2.0_double*(p13*p17 - p15) + c1*(p14*p17 - p15*p16 + p18*(3.0_double*p13 + c1*(2.0_double*p14 + p15*c1)))
          c5 = c1*(p13*p16 - p14 + c1*c5)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)**2
          ndh1drsdrsdn = one_third*(c2*h1 + 2.0_double*p22*t2*c3*c5)
          dh1dtot = 2.0_double*2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*c3*(1.0_double - c2)

          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + &
                        ec + ndecdn + h0 + h1 + ndh0dadadn + ndh1drsdrsdn - seven_sixths*t2*(dh0dtot + dh1dtot)
          dd(i1,i2,i3) = f07*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f10*(dh0dtot + dh1dtot)/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcd_pbe_pz_3d_i(np,n,dn,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.804_double      ! kappa
        real(double), parameter :: p02 = +0.21951_double    ! mu
        real(double), parameter :: p03 = +0.031091_double   ! gamma
        real(double), parameter :: p04 = +0.066725_double   ! beta
        real(double), parameter :: p05 = -0.09600_double    ! c1_ca
        real(double), parameter :: p06 = +0.06220_double    ! c2_ca
        real(double), parameter :: p07 = -0.02320_double    ! c3_ca
        real(double), parameter :: p08 = +0.00400_double    ! c4_ca
        real(double), parameter :: p09 = -0.28460_double    ! d1_ca
        real(double), parameter :: p10 = +1.05290_double    ! d2_ca
        real(double), parameter :: p11 = +0.33340_double    ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: rs, rsln, rssq, s, s2, t, t2, c1, c2, c3, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p05 + p07*rs + (p06 + p08*rs)*rsln
            ndecdn = -one_third*(p06 + (p07 + p08)*rs + p08*rs*rsln)
          else
            ec = p09/(1.0_double + p10*rssq + p11*rs)
            ndecdn = one_sixth*p09*(p10*rssq + 2.0_double*p11*rs)/(1.0_double + p10*rssq + p11*rs)**2
          end if
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c2 = p04/p03/(c1 - 1.0_double)*t2
          c3 = (1.0_double + c2*(1.0_double + c2))**2 + p04/p03*t2*(1.0_double + c2)*(1.0_double + c2*(1.0_double + c2))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          ndhdadadn = -c1*c2**3*ndecdn*(2.0_double + c2)/c3
          dhdtot = p04*(4.0_double + 8.0_double*c2)/c3

          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i1,i2,i3) = f05*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f06*dhdtot/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcd_pbe_pw_3d_i(np,n,dn,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, d, dd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.804_double       ! kappa
        real(double), parameter :: p02 = +0.21951_double     ! mu
        real(double), parameter :: p03 = +0.031091_double    ! gamma
        real(double), parameter :: p04 = +0.066725_double    ! beta
        real(double), parameter :: p05 = +0.0310907_double   ! a
        real(double), parameter :: p06 = +0.21370_double     ! a1
        real(double), parameter :: p07 = +7.59570_double     ! b1
        real(double), parameter :: p08 = +3.58760_double     ! b2
        real(double), parameter :: p09 = +1.63820_double     ! b3
        real(double), parameter :: p10 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06
        real(double) :: one_sixth, seven_sixths, one_third, two_thirds, four_thirds, five_thirds
        real(double) :: ex, fx, dfxdsos, ec, ndecdn, h, ndhdadadn, dhdtot
        real(double) :: s, s2, t, t2, c1, c2, c3, c4, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)
        f06 = 1.0_double/(16.0_double*(3.0_double/pi)**one_third)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1
          dfxdsos = 2.0_double*p02/x1**2

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p05*c1*(p07 + c1*(p08 + c1*(p09 + p10*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          c4 = -(1.0_double + p06*c1**2)*c1*(p07 + c1*(2.0_double*p08 + c1*(3.0_double*p09 + 4.0_double*p10*c1)))/ &
                 (c2*(1.0_double + c2))
          ec = -4.0_double*p05*(1.0_double + p06*c1**2)*c3
          ndecdn = four_thirds*p05*(p06*c1**2*c3 + p05*c4)
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c3 = p04/p03/(c1 - 1.0_double)*t2
          c4 = (1.0_double + c3*(1.0_double + c3))**2 + p04/p03*t2*(1.0_double + c3)*(1.0_double + c3*(1.0_double + c3))
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c3)/(1.0_double + c3*(1.0_double + c3)))
          ndhdadadn = -c1*c3**3*ndecdn*(2.0_double + c3)/c4
          dhdtot = p04*(4.0_double + 8.0_double*c3)/c4

          d(i1,i2,i3) = four_thirds*ex*(fx - s2*dfxdsos) + ec + ndecdn + h + ndhdadadn - seven_sixths*t2*dhdtot
          dd(i1,i2,i3) = f05*ex/n(i1,i2,i3)**five_thirds*dfxdsos + f06*dhdtot/n(i1,i2,i3)**four_thirds

        end do
        end do
        end do

      end subroutine

      subroutine xcd_am05_pz_3d_i(np,n,dn,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, d, dd

        integer :: i, i1, i2, i3
        real(double), parameter :: p01 = +2.804_double     ! alpha
        real(double), parameter :: p02 = +0.8098_double    ! gamma
        real(double), parameter :: p03 = +0.7168_double    ! c
        real(double), parameter :: p04 = -0.09600_double   ! c1_ca
        real(double), parameter :: p05 = +0.06220_double   ! c2_ca
        real(double), parameter :: p06 = -0.02320_double   ! c3_ca
        real(double), parameter :: p07 = +0.00400_double   ! c4_ca
        real(double), parameter :: p08 = -0.28460_double   ! d1_ca
        real(double), parameter :: p09 = +1.05290_double   ! d2_ca
        real(double), parameter :: p10 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: rs, rsln, rssq, s, s2, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do i = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p04 + p06*rs + (p05 + p07*rs)*rsln
            ndecdn = -one_third*(p05 + (p06 + p07)*rs + p07*rs*rsln)
          else
            ec = p08/(1.0_double + p09*rssq + p10*rs)
            ndecdn = one_sixth*p08*(p09*rssq + 2.0_double*p10*rs)/(1.0_double + p09*rssq + p10*rs)**2
          end if

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i1,i2,i3)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          d(i1,i2,i3) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i1,i2,i3) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i1,i2,i3))

        end do
        end do
        end do

200     call sync_config_process_errors()
        if (error("Exit functional_mod::xcd_am05_pz_3d_i")) continue

      end subroutine

      subroutine xcd_am05_pw_3d_i(np,n,dn,d,dd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, d, dd

        integer :: i, i1, i2, i3
        real(double), parameter :: p01 = +2.804_double       ! alpha
        real(double), parameter :: p02 = +0.8098_double      ! gamma
        real(double), parameter :: p03 = +0.7168_double      ! c
        real(double), parameter :: p04 = +0.0310907_double   ! a
        real(double), parameter :: p05 = +0.21370_double     ! a1
        real(double), parameter :: p06 = +7.59570_double     ! b1
        real(double), parameter :: p07 = +3.58760_double     ! b2
        real(double), parameter :: p08 = +1.63820_double     ! b3
        real(double), parameter :: p09 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08, f09
        real(double) :: one_sixth, one_quarter, three_quarters, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, ndecdn, xs, kf, hx, hxsos, hc, hcsos, fsos, xsos
        real(double) :: s, s2, c1, c2, c3, c4, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        three_quarters = 3.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4
        f09 = (3.0_double*pi**2)**one_third

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do i = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p04*c1*(p06 + c1*(p07 + c1*(p08 + p09*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p04*(1.0_double + p05*c1**2)*c3
          c4 = -(1.0_double + p05*c1**2)*c1*(p06 + c1*(2.0_double*p07 + c1*(3.0_double*p08 + 4.0_double*p09*c1)))/ &
                    (c2*(1.0_double + c2))
          ndecdn = four_thirds*p04*(p05*c1**2*c3 + p04*c4)

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02
          kf = f09*n(i1,i2,i3)**one_third
          xsos = -2.0_double*p01*xs**2
          fsos = p03/x5**2*(2.0_double - x3*((1.0_double - p03*s2)*(1.0_double + x4)**one_quarter + &
                   (1.0_double + p03*s2)*(1.0_double + 1.5_double*x4)/(1.0_double + x4)**three_quarters/(1.0_double + x2)))
          hxsos = (1.0_double - xs)*fsos - (fx - 1.0_double)*xsos
          hcsos = xsos*(1.0_double - p02)

          d(i1,i2,i3) = four_thirds*ex*hx + (ec + ndecdn)*hc - four_thirds*s2*(ex*hxsos + ec*hcsos)
          dd(i1,i2,i3) = (ex*hxsos + ec*hcsos)/((2.0_double*kf)**2*n(i1,i2,i3))

        end do
        end do
        end do

200     call sync_config_process_errors()
        if (error("Exit functional_mod::xcd_am05_pw_3d_i")) continue

      end subroutine

      subroutine xcd_blyp_3d_i(np,n,dn,lp,d,dd,ddd)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, lp, d, dd, ddd

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double    ! c1
        real(double), parameter :: p02 = +0.19645_double    ! c2
        real(double), parameter :: p03 = +7.79560_double    ! c3
        real(double), parameter :: p04 = +0.04918_double    ! a
        real(double), parameter :: p05 = +0.132_double      ! b
        real(double), parameter :: p06 = +0.2533_double     ! c
        real(double), parameter :: p07 = +0.349_double      ! d
        real(double), parameter :: p08 = +2.871234_double   ! cf
        real(double) :: f01, f02, f03
        real(double) :: one_third, two_thirds, four_thirds, five_thirds, seventeen_thirds
        real(double) :: ex, tb, dfxdsos, dfcdn, dfcddnodn, dfcdd2n
        real(double) :: s, s2, c1, c2, c3, c4, c5, x1, x2, x3

        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        seventeen_thirds = 17.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f03 = 1.0_double/(4.0_double*(3.0_double*pi**2)**two_thirds)

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f02*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = 1.0_double - p02*p03*s2/x1
          tb = 1.0_double - p01*s2*x3/x2**2
          dfxdsos = p01*(x3/x2 + 1.0_double)/x2

          ! correlation
          c1 = dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          c2 = lp(i1,i2,i3)/n(i1,i2,i3)**five_thirds
          c3 = n(i1,i2,i3)**one_third
          c4 = seventeen_thirds*(p06*p07 + c3*(p06 - 4.0_double*p07 - 5.0_double*c3))*c1**2
          c4 = c4 - 7.0_double*(p06*p07 + c3*(p06 - p07 - 2.0_double*c3))*c2
          c5 = p05*exp(-p06/c3)
          c4 = c4*c5/24.0_double
          c5 = c5*p08*(p06*p07 + c3*(p06 + 4.0_double*p07 + 3.0_double*c3))
          c4 = c4 - c5 - c3*(4.0_double*p07 + 3.0_double*c3)
          dfcdn = 2.0_double*p04*c4/(3.0_double*(p07 + c3)**2)
          c4 = p05/c3**5/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcddnodn = 34.0_double*p04/36.0_double*c4
          c4 = p05/c3**2/(1.0_double + p07/c3)*exp(-p06/c3)
          dfcdd2n = -14.0_double*p04/24.0_double*c4

          d(i1,i2,i3) = four_thirds*ex*tb + dfcdn
          dd(i1,i2,i3) = f03*ex/n(i1,i2,i3)**five_thirds*dfxdsos + dfcddnodn
          ddd(i1,i2,i3) = dfcdd2n

        end do
        end do
        end do

      end subroutine

      subroutine xce_lda_pz_3d_i(np,n,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = -0.09600_double   ! c1_ca
        real(double), parameter :: p02 = +0.06220_double   ! c2_ca
        real(double), parameter :: p03 = -0.02320_double   ! c3_ca
        real(double), parameter :: p04 = +0.00400_double   ! c4_ca
        real(double), parameter :: p05 = -0.28460_double   ! d1_ca
        real(double), parameter :: p06 = +1.05290_double   ! d2_ca
        real(double), parameter :: p07 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02
        real(double) :: one_third
        real(double) :: ex, ec
        real(double) :: rs, rsln, rssq

        one_third = 1.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p01 + p03*rs + (p02 + p04*rs)*rsln
          else
            ec = p05/(1.0_double + p06*rssq + p07*rs)
          end if

          e(i1,i2,i3) = ex + ec

        end do
        end do
        end do

      end subroutine

      subroutine xce_lda_pw_3d_i(np,n,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.0310907_double   ! a
        real(double), parameter :: p02 = +0.21370_double     ! a1
        real(double), parameter :: p03 = +7.59570_double     ! b1
        real(double), parameter :: p04 = +3.58760_double     ! b2
        real(double), parameter :: p05 = +1.63820_double     ! b3
        real(double), parameter :: p06 = +0.49294_double     ! b4
        real(double) :: f01, f02
        real(double) :: one_sixth, one_third
        real(double) :: ex, ec
        real(double) :: c1, c2, c3

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p01*c1*(p03 + c1*(p04 + c1*(p05 + p06*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p01*(1.0_double + p02*c1**2)*c3

          e(i1,i2,i3) = ex + ec

        end do
        end do
        end do

      end subroutine

      subroutine xce_lsda_pw_3d_i(np,ns,zeta,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: ns, zeta, e

        integer :: i1, i2, i3
        real(double), parameter :: p10 =  +0.0310907_double   ! a_0
        real(double), parameter :: p20 =  +0.2137000_double   ! a1_0
        real(double), parameter :: p30 =  +7.5957000_double   ! b1_0
        real(double), parameter :: p40 =  +3.5876000_double   ! b2_0
        real(double), parameter :: p50 =  +1.6382000_double   ! b3_0
        real(double), parameter :: p60 =  +0.4929400_double   ! b4_0
        real(double), parameter :: p11 =  +0.0155450_double   ! a_1
        real(double), parameter :: p21 =  +0.2054800_double   ! a1_1
        real(double), parameter :: p31 = +14.1189000_double   ! b1_1
        real(double), parameter :: p41 =  +6.1977000_double   ! b2_1
        real(double), parameter :: p51 =  +3.3662000_double   ! b3_1
        real(double), parameter :: p61 =  +0.6251700_double   ! b4_1
        real(double), parameter :: p12 =  +0.0168870_double   ! a_a
        real(double), parameter :: p22 =  +0.1112500_double   ! a1_a
        real(double), parameter :: p32 = +10.3570000_double   ! b1_a
        real(double), parameter :: p42 =  +3.6231000_double   ! b2_a
        real(double), parameter :: p52 =  +0.8802600_double   ! b3_a
        real(double), parameter :: p62 =  +0.4967100_double   ! b4_a
        real(double) :: one_sixth, one_third, four_thirds
        real(double) :: f01, f02, f03
        real(double) :: c1, c1_2, c20, c21, c22, c30, c31, c32
        real(double) :: f, d2fdz2
        real(double) :: ex
        real(double) :: ec, ec0, ec1, ec2
        real(double) :: omz, opz, z, z4

        one_sixth = 1.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 2.0_double**four_thirds - 2.0_double

        d2fdz2 = (8.0_double/9.0_double)/f03  ! d2f/dz2 | z=0

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          z = zeta(i1,i2,i3)
          z4 = z**4

          opz = 1.0_double + z
          omz = 1.0_double - z

          f = (opz**four_thirds + omz**four_thirds - 2.0_double)/f03

          c1 = f02/ns(i1,i2,i3)**one_sixth  ! sqrt(rs)
          c1_2 = c1**2                      ! rs

          ! exchange
          ex = 0.50_double*(f01*ns(i1,i2,i3)**one_third)*(opz**four_thirds + omz**four_thirds )

          ! correlation
          c20 = 2.0_double*p10*c1*(p30 + c1*(p40 + c1*(p50 + p60*c1)))
          c21 = 2.0_double*p11*c1*(p31 + c1*(p41 + c1*(p51 + p61*c1)))
          c22 = 2.0_double*p12*c1*(p32 + c1*(p42 + c1*(p52 + p62*c1)))
          c30 = log(1.0_double + 1.0_double/c20)
          c31 = log(1.0_double + 1.0_double/c21)
          c32 = log(1.0_double + 1.0_double/c22)
          ec0 = -4.0_double*p10*(1.0_double + p20*c1_2)*c30
          ec1 = -4.0_double*p11*(1.0_double + p21*c1_2)*c31
          ec2 = -4.0_double*p12*(1.0_double + p22*c1_2)*c32
          ec = ec0 + ec2*(f/d2fdz2)*(1.0_double - z4) + (ec1 - ec0)*f*z4

          ! combine
          e(i1,i2,i3) = ex + ec

        end do
        end do
        end do

100     if (error("Exit functional_mod::xce_lsda_pw_3d_i")) continue

      end subroutine

      subroutine xce_pw91_pz_3d_i(np,n,dn,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = -0.09600_double       ! c1_ca
        real(double), parameter :: p08 = +0.06220_double       ! c2_ca
        real(double), parameter :: p09 = -0.02320_double       ! c3_ca
        real(double), parameter :: p10 = +0.00400_double       ! c4_ca
        real(double), parameter :: p11 = -0.28460_double       ! d1_ca
        real(double), parameter :: p12 = +1.05290_double       ! d2_ca
        real(double), parameter :: p13 = +0.33340_double       ! d3_ca
        real(double), parameter :: p14 = +0.002568_double      ! c0
        real(double), parameter :: p15 = +0.023266_double      ! c1
        real(double), parameter :: p16 = +7.389e-6_double      ! c2
        real(double), parameter :: p17 = +8.723_double         ! c3
        real(double), parameter :: p18 = +0.472_double         ! c4
        real(double), parameter :: p19 = +7.389e-2_double      ! c5
        real(double), parameter :: p20 = -0.001667212_double   ! cx
        real(double), parameter :: p21 = +0.004235_double      ! cc0
        real(double), parameter :: p22 = +0.09_double          ! alfa
        real(double), parameter :: p23 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f08, f09
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h0, h1
        real(double) :: rs, rsln, rssq, s, s2, s4, t, t2, c1, c2, c3, c4, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p21*p23                                                 ! beta
        f06 = 2.0_double*p22/f05                                      ! aob
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p07 + p09*rs + (p08 + p10*rs)*rsln
          else
            ec = p11/(1.0_double + p12*rssq + p13*rs)
          end if
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          c1 = f08/n(i1,i2,i3)**one_third
          c2 = 100.0_double*t2*f09/n(i1,i2,i3)**one_third
          c3 = exp(-c2)
          c4 = -p20 + (p14 + p15*c1 + p16*c1**2)/(1.0_double + p17*c1 + p18*c1**2 + p19*c1**3)
          h1 = 2.0_double*p23*(c4 - p21 - 3.0_double/7.0_double*p20)*t2*c3

          e(i1,i2,i3) = ex*fx + ec + h0 + h1

        end do
        end do
        end do

      end subroutine

      subroutine xce_pw91_pw_3d_i(np,n,dn,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double       ! c1
        real(double), parameter :: p02 = +0.19645_double       ! c2
        real(double), parameter :: p03 = +7.79560_double       ! c3
        real(double), parameter :: p04 = +0.15084_double       ! c4
        real(double), parameter :: p05 = +100.0_double         ! c5
        real(double), parameter :: p06 = +0.004_double         ! c6
        real(double), parameter :: p07 = + 0.0310907_double    ! a
        real(double), parameter :: p08 = +0.21370_double       ! a1
        real(double), parameter :: p09 = +7.59570_double       ! b1
        real(double), parameter :: p10 = +3.58760_double       ! b2
        real(double), parameter :: p11 = +1.63820_double       ! b3
        real(double), parameter :: p12 = +0.49294_double       ! b4
        real(double), parameter :: p13 = +0.002568_double      ! c0
        real(double), parameter :: p14 = +0.023266_double      ! c1
        real(double), parameter :: p15 = +7.389e-6_double      ! c2
        real(double), parameter :: p16 = +8.723_double         ! c3
        real(double), parameter :: p17 = +0.472_double         ! c4
        real(double), parameter :: p18 = +7.389e-2_double      ! c5
        real(double), parameter :: p19 = -0.001667212_double   ! cx
        real(double), parameter :: p20 = +0.004235_double      ! cc0
        real(double), parameter :: p21 = +0.09_double          ! alfa
        real(double), parameter :: p22 = +15.75592_double      ! nu
        real(double) :: f01, f02, f03, f04, f05, f06, f08, f09
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h0, h1
        real(double) :: s, s2, s4, t, t2, c1, c2, c3, c4, x1, x2, x3, x4

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1
        f05 = p20*p22                                                 ! beta
        f06 = 2.0_double*p21/f05                                      ! aob
        f08 = (3.0_double/(4.0_double*pi))**one_third
        f09 = 4.0_double/pi/(3.0_double*pi**2)**one_third

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          s4 = s**4
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          x3 = p01 - p04*exp(-p05*s2)
          x4 = x2 + p06*s4
          fx = (x2 + x3*s2)/x4

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p07*c1*(p09 + c1*(p10 + c1*(p11 + p12*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p07*(1.0_double + p08*c1**2)*c3
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-f06/f05*ec/2.0_double)
          c2 = f06/(c1 - 1.0_double)*t2
          h0 = 2.0_double*f05/f06*log(1.0_double + f06*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))
          c1 = f08/n(i1,i2,i3)**one_third
          c2 = 100.0_double*t2*f09/n(i1,i2,i3)**one_third
          c3 = exp(-c2)
          c4 = -p19 + (p13 + p14*c1 + p15*c1**2)/(1.0_double + p16*c1 + p17*c1**2 + p18*c1**3)
          h1 = 2.0_double*p22*(c4 - p20 - 3.0_double/7.0_double*p19)*t2*c3

          e(i1,i2,i3) = ex*fx + ec + h0 + h1

        end do
        end do
        end do

      end subroutine

      subroutine xce_pbe_pz_3d_i(np,n,dn,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.804_double      ! kappa
        real(double), parameter :: p02 = +0.21951_double    ! mu
        real(double), parameter :: p03 = +0.031091_double   ! gamma
        real(double), parameter :: p04 = +0.066725_double   ! beta
        real(double), parameter :: p05 = -0.09600_double    ! c1_ca
        real(double), parameter :: p06 = +0.06220_double    ! c2_ca
        real(double), parameter :: p07 = -0.02320_double    ! c3_ca
        real(double), parameter :: p08 = +0.00400_double    ! c4_ca
        real(double), parameter :: p09 = -0.28460_double    ! d1_ca
        real(double), parameter :: p10 = +1.05290_double    ! d2_ca
        real(double), parameter :: p11 = +0.33340_double    ! d3_ca
        real(double) :: f01, f02, f03, f04
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h
        real(double) :: rs, rsln, rssq, s, s2, t, t2, c1, c2, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p05 + p07*rs + (p06 + p08*rs)*rsln
          else
            ec = p09/(1.0_double + p10*rssq + p11*rs)
          end if
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c2 = p04/p03/(c1 - 1.0_double)*t2
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c2)/(1.0_double + c2*(1.0_double + c2)))

          e(i1,i2,i3) = ex*fx + ec + h

        end do
        end do
        end do

      end subroutine

      subroutine xce_pbe_pw_3d_i(np,n,dn,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.804_double       ! kappa
        real(double), parameter :: p02 = +0.21951_double     ! mu
        real(double), parameter :: p03 = +0.031091_double    ! gamma
        real(double), parameter :: p04 = +0.066725_double    ! beta
        real(double), parameter :: p05 = +0.0310907_double   ! a
        real(double), parameter :: p06 = +0.21370_double     ! a1
        real(double), parameter :: p07 = +7.59570_double     ! b1
        real(double), parameter :: p08 = +3.58760_double     ! b2
        real(double), parameter :: p09 = +1.63820_double     ! b3
        real(double), parameter :: p10 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04
        real(double) :: one_sixth, seven_sixths, one_third, four_thirds
        real(double) :: ex, fx, ec, h
        real(double) :: s, s2, t, t2, c1, c2, c3, x1

        one_sixth = 1.0_double/6.0_double
        seven_sixths = 7.0_double/6.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 1.0_double/(4.0_double*(3.0_double/pi)**one_sixth)      ! t1

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = 1.0_double + p02*s2/p01
          fx = 1.0_double + p01 - p01/x1

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p05*c1*(p07 + c1*(p08 + c1*(p09 + p10*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p05*(1.0_double + p06*c1**2)*c3
          t = f04*dn(i1,i2,i3)/n(i1,i2,i3)**seven_sixths
          t2 = t**2
          c1 = exp(-ec/2.0_double/p03)
          c3 = p04/p03/(c1 - 1.0_double)*t2
          h = 2.0_double*p03*log(1.0_double + p04/p03*t2*(1.0_double + c3)/(1.0_double + c3*(1.0_double + c3)))

          e(i1,i2,i3) = ex*fx + ec + h

        end do
        end do
        end do

      end subroutine

      subroutine xce_am05_pz_3d_i(np,n,dn,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e

        integer :: i, i1, i2, i3
        real(double), parameter :: p01 = +2.804_double     ! alpha
        real(double), parameter :: p02 = +0.8098_double    ! gamma
        real(double), parameter :: p03 = +0.7168_double    ! c
        real(double), parameter :: p04 = -0.09600_double   ! c1_ca
        real(double), parameter :: p05 = +0.06220_double   ! c2_ca
        real(double), parameter :: p06 = -0.02320_double   ! c3_ca
        real(double), parameter :: p07 = +0.00400_double   ! c4_ca
        real(double), parameter :: p08 = -0.28460_double   ! d1_ca
        real(double), parameter :: p09 = +1.05290_double   ! d2_ca
        real(double), parameter :: p10 = +0.33340_double   ! d3_ca
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08
        real(double) :: one_sixth, one_quarter, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, xs, hx, hc
        real(double) :: rs, rsln, rssq, s, s2, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = log((3.0_double/(4.0_double*pi)))
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do i = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          rsln = one_third*(f02 - log(n(i1,i2,i3)))
          rssq = exp(0.5_double*rsln)
          rs = rssq**2
          if (rs < 1.0_double) then
            ec = p04 + p06*rs + (p05 + p07*rs)*rsln
          else
            ec = p08/(1.0_double + p09*rssq + p10*rs)
          end if

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02

          e(i1,i2,i3) = ex*hx + ec*hc

        end do
        end do
        end do

200     call sync_config_process_errors()
        if (error("Exit functional_mod::xce_am05_pz_3d_i")) continue

      end subroutine

      subroutine xce_am05_pw_3d_i(np,n,dn,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, e

        integer :: i, i1, i2, i3
        real(double), parameter :: p01 = +2.804_double       ! alpha
        real(double), parameter :: p02 = +0.8098_double      ! gamma
        real(double), parameter :: p03 = +0.7168_double      ! c
        real(double), parameter :: p04 = +0.0310907_double   ! a
        real(double), parameter :: p05 = +0.21370_double     ! a1
        real(double), parameter :: p06 = +7.59570_double     ! b1
        real(double), parameter :: p07 = +3.58760_double     ! b2
        real(double), parameter :: p08 = +1.63820_double     ! b3
        real(double), parameter :: p09 = +0.49294_double     ! b4
        real(double) :: f01, f02, f03, f04, f05, f06, f07, f08
        real(double) :: one_sixth, one_quarter, one_third, two_thirds, four_thirds
        real(double) :: ex, fx, ec, xs, hx, hc
        real(double) :: s, s2, c1, c2, c3, x1, x2, x3, x4, x5

        one_sixth = 1.0_double/6.0_double
        one_quarter = 1.0_double/4.0_double
        one_third = 1.0_double/3.0_double
        two_thirds = 2.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = (3.0_double/(4.0_double*pi))**one_sixth
        f03 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1
        f04 = 2.0_double*sqrt(6.0_double)
        f05 = 1.45_double*exp(1.0_double) - 1.0_double
        f06 = 2.0_double*exp(1.0_double)
        f07 = 2.0_double*3.0_double**one_third
        f08 = (729.0_double/1024.0_double)/pi**4

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f03*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = s**1.5_double/f04
          if (x1 < 1.0e-20_double) then
            x2 = x1
          else
            if (x1 > f05) then
              x2 = log(x1)
              x2 = x2 - log(x2)
            else
              x2 = sqrt(f06*x1 + 2.0_double) - 1.0_double
            end if
            do i = 1,10
              if (x2 == -1.0_double) goto 100
              x3 = exp(x2)
              x4 = x2*x3 - x1
              x5 = x4/(x3*(x2 + 1.0_double) - 0.5_double*(x2 + 2.0_double)*x4/(x2 + 1.0_double))
              x2 = x2 - x5
              if (abs(x5) < 2.48e-14_double*(1.0_double + abs(x2))) goto 100
            end do
            if (error(abs(x5) >= 2.48e-14_double*(1.0_double + abs(x2)),"ERROR: AM05 iteration failed")) goto 200
          end if
100       if (s < 1.0e-14_double) then
            x3 = 1.0_double
          else
            x3 = f07*x2**two_thirds/s
          end if
          x4 = s2*f08*x3**2
          x5 = 1.0_double + p03*s2*x3*(1.0_double + x4)**one_quarter
          fx = (p03*s2 + 1.0_double)/x5

          ! correlation
          c1 = f02/n(i1,i2,i3)**one_sixth
          c2 = 2.0_double*p04*c1*(p06 + c1*(p07 + c1*(p08 + p09*c1)))
          c3 = log(1.0_double + 1.0_double/c2)
          ec = -4.0_double*p04*(1.0_double + p05*c1**2)*c3

          xs = 1.0_double - p01*s2/(1.0_double + p01*s2)
          hx = xs + (1.0_double - xs)*fx
          hc = xs + (1.0_double - xs)*p02

          e(i1,i2,i3) = ex*hx + ec*hc

        end do
        end do
        end do

200     call sync_config_process_errors()
        if (error("Exit functional_mod::xce_am05_pw_3d_i")) continue

      end subroutine

      subroutine xce_blyp_3d_i(np,n,dn,lp,e)
        integer, dimension(3) :: np
        real(double), dimension(np(1),np(2),np(3)) :: n, dn, lp, e

        integer :: i1, i2, i3
        real(double), parameter :: p01 = +0.27430_double    ! c1
        real(double), parameter :: p02 = +0.19645_double    ! c2
        real(double), parameter :: p03 = +7.79560_double    ! c3
        real(double), parameter :: p04 = +0.04918_double    ! a
        real(double), parameter :: p05 = +0.132_double      ! b
        real(double), parameter :: p06 = +0.2533_double     ! c
        real(double), parameter :: p07 = +0.349_double      ! d
        real(double), parameter :: p08 = +2.871234_double   ! cf
        real(double) :: f01, f02
        real(double) :: one_eighth, one_third, four_thirds, five_thirds
        real(double) :: ex, fx, ec
        real(double) :: s, s2, c1, c2, c3, c4, c5, x1, x2

        one_eighth = 1.0_double/8.0_double
        one_third = 1.0_double/3.0_double
        four_thirds = 4.0_double/3.0_double
        five_thirds = 5.0_double/3.0_double
        f01 = -3.0_double/2.0_double*(3.0_double/pi)**one_third
        f02 = 1.0_double/(2.0_double*(3.0_double*pi**2)**one_third)   ! s1

        do i3 = 1,np(3)
        do i2 = 1,np(2)
        do i1 = 1,np(1)

          ! exchange
          ex = f01*n(i1,i2,i3)**one_third
          s = f02*dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          s2 = s**2
          x1 = sqrt(1.0_double + (p03*s)**2)
          x2 = 1.0_double + p02*s*log(p03*s + x1)
          fx = 1.0_double + p01*s2/x2

          ! correlation
          c1 = dn(i1,i2,i3)/n(i1,i2,i3)**four_thirds
          c2 = lp(i1,i2,i3)/n(i1,i2,i3)**five_thirds
          c3 = n(i1,i2,i3)**one_third
          c4 = one_eighth*(c1**2 - c2)
          c5 = (p08*p05 + p05/9.0_double*(c2/2.0_double - 17.0_double*c4))*exp(-p06/c3)
          ec = -2.0_double*p04/(1.0_double + p07/c3)*(1.0_double + c5)

          e(i1,i2,i3) = ex*fx + ec

        end do
        end do
        end do

      end subroutine

      subroutine spin_functions_3d_i(n,ns,zeta)
        real(double), dimension(:,:,:) :: n
        real(double), dimension(:,:,:), pointer :: ns, zeta

        allocate( ns(size(n,1),size(n,2),size(n,3)) )
        allocate( zeta(size(n,1),size(n,2),size(n,3)) )

        ! sum of the spin densities
        call xcomm_rank_allreduce(XSGROUP,MPI_SUM,n,ns)

        ! fractional difference of the spin densities
        if (mpi_mysgroup() == 2) n = -1.0_double*n
        call xcomm_rank_allreduce(XSGROUP,MPI_SUM,n,zeta)
        zeta = zeta/ns
        if (mpi_mysgroup() == 2) n = -1.0_double*n

      end subroutine

      end module xc_functional_mod

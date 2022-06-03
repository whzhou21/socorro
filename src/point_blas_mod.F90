!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module point_blas_mod
!doc$ module point_blas_mod

!*******************************************************************
!
! point_blas_mod.f90 - contains misc routines designed to get around
!                  the F90 pointer vs. allocatable problem.
!
!  The module routines are provided so that you can still have type
!  checking.  The actual routines that do the work are listed 
!  after the end of the module definition.  These are the routines
!  that start with Kernel.  There is a slight performance
!  degradation associated with this 2-level method but I think it's
!  worth for the type checking ability.  The Kernel routines
!  are not in the scope of the point_blas module for a reason. 
!  This forces F90 to pass the arrays using F77 syntax and also
!  compiles the kernel routines assuming they are statically
!  created giving the best performance.  This is similar to 
!  typecasting a data type in C/C++.  
!
!  NOTE:  If the arrays are passed into the point_blas routines
!         using F90 array subscripting unpredictable things will
!         happen.           
!
!*******************************************************************

Use kind_mod

Implicit NONE!
Private

Integer, PARAMETER, PUBLIC :: POINT_F77 = 0  !** F77 style
Integer, PARAMETER, PUBLIC :: POINT_F90 = 1  !** F90 style


Logical :: DoF77  !** Do F77 type operations


Public :: point_MxM, Point_MxS, Point_MAMxV, Point_AccumDen, Point_Mode

Interface Point_MxM  
   Module Procedure Point_MxM_Z3_D3 
End Interface

Interface Point_MxS
   Module Procedure Point_MxS_Z3_D
   Module Procedure Point_MxS_D3_D
End Interface

Interface Point_MAMxV
   Module Procedure Point_MAMxV_Z2Z2D
End Interface

!*******************************************************************
Contains
!*******************************************************************

!*******************************************************************

Subroutine Point_Mode(mode)
  Integer, Intent(IN) :: mode

  DoF77 = (mode == POINT_F77)

  Return
End Subroutine
 
!*******************************************************************

Subroutine Point_MxM_Z3_D3(A,B)
  Complex(DOUBLE), Intent(INOUT) :: A(:,:,:)
  Real(DOUBLE),    Intent(IN)    :: B(:,:,:)

  if (DoF77) then
     Call Kernel_MxM_ZD(A,B,size(A))
  else
     A = A*B
  end if   

End Subroutine

!*******************************************************************

Subroutine Point_MxS_Z3_D(A,scale)
  Complex(DOUBLE), Intent(INOUT) :: A(:,:,:)
  Real(DOUBLE),    Intent(IN)    :: Scale

  If (DoF77) then
     Call Kernel_MxS_ZD(A,scale,size(A))
  else
     A = A*scale
  End if   

End Subroutine

!*******************************************************************

Subroutine Point_MxS_D3_D(A,scale)
  Real(DOUBLE), Intent(INOUT) :: A(:,:,:)
  Real(DOUBLE),    Intent(IN)    :: Scale

  If (DoF77) then  
     Call Kernel_MxS_DD(A,scale,size(A))
  else
     A = A*scale 
  End If

End Subroutine

!*******************************************************************

Subroutine Point_MAMxV_Z2Z2D(A,B,V)
  Complex(DOUBLE), Intent(INOUT) :: A(:,:)
  Complex(DOUBLE), Intent(IN)    :: B(:,:)
  Real(DOUBLE),    Intent(IN)    :: V(:)

  integer :: i, cols

  if (DoF77) then  
     Call Kernel_MAMxV_Z2Z2D(A,B,V, size(A,1), size(A,2))
  else
     cols = size(A,2)
     Do i=1, Cols
        A(:,i) = A(:,i) + B(:,i)*V
     End Do
  End If

End Subroutine

!*******************************************************************

Subroutine Point_AccumDen(A,B, scale)
  Real(DOUBLE),    Intent(INOUT) :: A(:,:,:)
  Complex(DOUBLE), Intent(IN)    :: B(:,:,:)
  Real(DOUBLE),    Intent(IN)    :: scale
    
  if (DoF77) then  
     Call Kernel_AccumDen(A,B,scale,size(A))
  else
     A = A + scale * (abs(B))**2
  End If
   
End Subroutine

End Module

!*******************************************************************
!*******************************************************************

Subroutine Kernel_MxM_ZD(A,B,N)
  Use kind_mod
  Integer         :: N 
  Complex(Double) :: A(N)
  Real(DOUBLE)    :: B(N)
  
  A(1:N) = A(1:N)*B(1:N)
End Subroutine

!*******************************************************************

Subroutine Kernel_MxS_ZD(A,Scale,N)
  Use kind_mod
  Integer         :: N 
  Complex(Double) :: A(N)
  Real(DOUBLE)    :: Scale
  
  A(1:N) = A(1:N)*scale
End Subroutine

!*******************************************************************

Subroutine Kernel_MxS_DD(A,Scale,N)
  Use kind_mod
  Integer         :: N 
  Real(Double) :: A(N)
  Real(DOUBLE)    :: Scale
  
  A(1:N) = A(1:N)*scale
End Subroutine

!*******************************************************************

Subroutine Kernel_MAMxV_Z2Z2D(A,B,V, rows, Cols)
  Use kind_mod
  Integer         :: Rows
  Integer         :: Cols
  Complex(Double) :: A(Rows,Cols)
  Complex(Double) :: B(rows, Cols)
  Real(DOUBLE)    :: V(Rows)
  
  Integer :: i

  Do i=1, Cols
     A(1:rows,i) = A(1:rows,i) + B(1:rows,i)*V(1:rows)
  End Do   

End Subroutine  

!*******************************************************************

Subroutine Kernel_AccumDen(A,B, scale, N)
  USe kind_mod
  Integer          :: N
  Real(DOUBLE)     :: A(N)
  Complex(DOUBLE)  :: B(N)
  Real(DOUBLE)     :: scale

  A(1:N) = A(1:N) + scale * (abs(B(1:N))**2)

End Subroutine

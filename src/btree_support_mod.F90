! Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms
! of Contract DE-NA0003525 with NTESS, the U.S. Government retains certains rights to this software.

!******************************************************************************
!
! File : btree_support_mod.f90
!   by : Alan Tackett
!   on : 07/12/01
!  for : Socorro Timing routines
!
!  Btree data structures and support routines
!
!******************************************************************************

Module btree_support_mod

Use kind_mod

Implicit NONE!!!!!!!!!!!!!

Type Timer_Type                 !** Timer definition
   Integer          :: state    !** Timer state
   Integer          :: N        !** Number of invocations
   Real(double)     :: elapsed  !** Elapsed time
   Real(double)     :: start_time !** Start time of last invocation
   Character(len=40):: Name     !** Timer name (Duplicated in watch)
End Type


Type BinaryTreeData   !** Define the information stored
     Type (BinaryTreeData), Pointer :: Level_Next !*Link to next point/level

     character*80               :: Name        !* Name of the timer
     Type (Timer_type), Pointer :: timer       !* Timer pointer
End type

Type BinaryTreeNode   !** Define the btree node
     Type (BinaryTreeNode), Pointer :: LeftChild      !* children
     Type (BinaryTreeNode), Pointer :: RightChild
     Type (BinaryTreeData), Pointer :: Info           !* pointer to data
     Integer                        :: BalanceFactor  !* relative branch ht
End Type


!******************************************************************************
Contains
!******************************************************************************


!******************************************************************************
!
!  Compare - Compares two data items and determines their relationship.
!              The ordering is determined as |G|,Z,Y,X.
!
!    P1, P2  - Data items to compare
!    
!    Return Values
!         -1 -- P1<P2
!          0 -- P1==P2
!          1 -- P1>P2
!
!******************************************************************************

Integer Function Compare(pd1, pd2)
  Type (BinaryTreeData) :: pd1, pd2


  if (PD1%name < PD2%name) then
     Compare = -1
  elseif (PD1%name == PD2%name) then
     Compare = 0
  else
     Compare = 1
  End If

  Return
End Function

End Module



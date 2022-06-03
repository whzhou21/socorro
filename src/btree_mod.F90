!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

!******************************************************************************
!
! File : btree_mod.f90
!   by : Alan Tackett
!   on : 08/06/98
!  for : WAE program
!
!  Module for implementing a height balanced binary tree.
!
!  The code used here was based on that given in "Introduction to data 
!  structures and algorithm analysis" by Thomas L. Naps in chapter 7
!
!******************************************************************************

Module btree_mod

Use btree_support_mod   !** This module contains the data types and 
                        !** support routines for the specific application

Implicit NONE!!!!!!!!!!!!
Private

!doc$
Public :: FreeTree, SearchBtree, InsertNode, CreateBtree

!cod$

!******************************************************************************
Contains
!******************************************************************************

!******************************************************************************
!
!  FreeTree - Frees the tree 
!
!     Root     - Binary tree root           
!
!******************************************************************************

Recursive Subroutine FreeTree(Root, FreeData)
  Type (BinaryTreeNode), Pointer  :: Root
  Logical,             Intent(IN) :: FreeData

  If (Associated(Root)) then
     Call FreeTree(Root%LeftChild, FreeData)
     Call FreeTree(Root%RightChild, FreeData)
     If (FreeData) then
        Deallocate(Root%Info%timer)
	DeAllocate(Root%Info)
     End If
     DeAllocate(Root)
  EndIf

  Return
End Subroutine

!******************************************************************************
!
!  SearchBTree - Searches the b-tree for a given node
!
!    Root - Binary tree root
!    Item - Data to match
!    Result  - Pointer to matching node if found
!    Success - Determines if item is successfully inserted
!
!******************************************************************************

Subroutine SearchBTree(Root, Item, Result, Success)
  Type (BinaryTreeNode), Pointer     :: Root
  Type (BinaryTreeData), TARGET      :: Item
  Type (BinaryTreeNode), Pointer     :: Result
  Logical,               Intent(OUT) :: Success

  Integer :: CmpResult
  
  Success = .FALSE.

  If (.NOT. Associated(Root)) then
     Nullify(Result)
     RETURN
  End If

  Result => Root

  Do While (Associated(Result) .AND. (.NOT. Success))
     CmpResult = Compare(Item, Result%Info)

     Select Case (CmpResult)
        Case (-1)
           Result => Result%LeftChild
        Case (0)
           Success = .TRUE.
        Case (1)
           Result => Result%RightChild
     End Select
  End Do

  Return
End Subroutine
  
!******************************************************************************
!
!  LeftOfLeft - AVL rotation case 1 (p351)
!
!    Pivot - Pivot point
!
!******************************************************************************

Subroutine LeftOfLeft(Pivot)
  Type (BinaryTreeNode), Pointer :: Pivot

  Type (BinaryTreeNode), Pointer :: P, Q

  P => Pivot%LeftChild
  Q => P%RightChild

  P%RightChild => Pivot
  Pivot%LeftChild => Q
  Pivot => P

  Pivot%BalanceFactor = 0
  Pivot%RightChild%BalanceFactor = 0

  Return
End Subroutine

!******************************************************************************
!
!  RightOfright - AVL rotation case 2 (p353)
!
!    Pivot - Pivot point
!
!******************************************************************************

Subroutine RightOfRight(Pivot)
  Type (BinaryTreeNode), Pointer :: Pivot

  Type (BinaryTreeNode), Pointer :: P, Q

  P => Pivot%RightChild
  Q => P%LeftChild

  P%LeftChild => Pivot
  Pivot%RightChild => Q
  Pivot => P

  Pivot%BalanceFactor = 0
  Pivot%LeftChild%BalanceFactor = 0

  Return
End Subroutine

!******************************************************************************
!
!  RightOfLeft - AVL rotation case 3 (p354)
!
!    Pivot - Pivot point
!
!******************************************************************************

Subroutine RightOfLeft(Pivot)
  Type (BinaryTreeNode), Pointer :: Pivot

  Type (BinaryTreeNode), Pointer :: X, Y

  X => Pivot%LeftChild
  Y => X%RightChild

If (.NOT. Associated(Y)) Write(*,*) 'Y => NULL'

  Pivot%LeftChild => Y%RightChild
  X%RightChild => Y%LeftChild
  Y%LeftChild => X
  Y%RightChild => Pivot
  Pivot => Y


  Select Case (Pivot%BalanceFactor)
     Case (0)
        Pivot%LeftChild%BalanceFactor = 0
        Pivot%RightChild%BalanceFactor = 0
     Case (1)
        Pivot%BalanceFactor = 0
        Pivot%LeftChild%BalanceFactor = 0
        Pivot%RightChild%BalanceFactor = -1
     Case Default
        Pivot%BalanceFactor = 0
        Pivot%LeftChild%BalanceFactor = 1
        Pivot%RightChild%BalanceFactor = 0
   End Select

   Return
End Subroutine

!******************************************************************************
!
!  LeftOfRight - AVL rotation case 4 (p355)
!
!    Pivot - Pivot point
!
!******************************************************************************

Subroutine LeftOfRight(Pivot)
  Type (BinaryTreeNode), Pointer :: Pivot

  Type (BinaryTreeNode), Pointer :: X, Y

  X => Pivot%RightChild
  Y => X%LeftChild

  Pivot%RightChild => Y%LeftChild
  X%LeftChild => Y%RightChild
  Y%RightChild => X
  Y%LeftChild => Pivot
  Pivot => Y

  Select Case (Pivot%BalanceFactor)
     Case (0)
        Pivot%LeftChild%BalanceFactor = 0
        Pivot%RightChild%BalanceFactor = 0
     Case (-1)
        Pivot%BalanceFactor = 0
        Pivot%LeftChild%BalanceFactor = 1
        Pivot%RightChild%BalanceFactor = 0
     Case Default
        Pivot%BalanceFactor = 0
        Pivot%LeftChild%BalanceFactor = 0
        Pivot%RightChild%BalanceFactor = -1
   End Select

   Return
End Subroutine

!******************************************************************************
!
!  InsertNode - Inserts a node in the b-tree
!
!    Root - Binary tree root
!    Item - Data to be inserted
!    Success - Determines if item is successfully inserted
!
!******************************************************************************

Subroutine InsertNode(Root, Item, Success)
  Type (BinaryTreeNode), Pointer     :: Root
  Type (BinaryTreeData), Pointer     :: Item
!**  Type (BinaryTreeData), TARGET     :: Item
  Logical,               Intent(OUT) :: Success

  Type (BinaryTreeNode), Pointer :: P, Piv, PivParent, InP, InParent, Q
  Integer :: err
  
  Success = .TRUE.

  Allocate(P, STAT=err)
  If (err /= 0) then
    Write(*,*) 'InsertNode: Alloc Error=',err
    STOP
  End If

  P%Info => Item
  P%BalanceFactor = 0
  Nullify(P%LeftChild)
  Nullify(P%RightChild)

  If (.NOT. Associated(Root)) then
     Root => P
     RETURN
  End If

  InP => Root
  Piv => Root
  Nullify(InParent)
  Nullify(PivParent)

  !*** Search for insertion point and pivot ****
  Do While (Associated(InP))
     If (InP%BalanceFactor /= 0) then
        Piv => InP
        PivParent => InParent
     End If

     InParent => InP

     If (Compare(Item, InP%Info) == -1) then
        InP => InP%LeftChild
     else
        InP => InP%RightChild
     End If
  End Do

  !*** Insert the node as a child of InParent ***
  If (Compare(Item, InParent%Info) == -1) then
     InParent%LeftChild => P
  else
     InParent%RightChild => P
  End If

  !*** Now recompute the balance factors between Piv and InParent ***
  Q => Piv
  If (Compare(Item, Q%Info) == -1) then
     Q%BalanceFactor = Q%BalanceFactor + 1
     Q => Q%LeftChild
  else
     Q%BalanceFactor = Q%BalanceFactor - 1
     Q => Q%RightChild
  End If

  Do While (.NOT. Associated(Q, P))
     If (Compare(Item, Q%Info) == -1) then
        Q%BalanceFactor = Q%BalanceFactor + 1
        Q => Q%LeftChild
     else
        Q%BalanceFactor = Q%BalanceFactor - 1
        Q => Q%RightChild
     End If
  End Do

  !*** Check to see if an AVL rotation is needed ****
  If ((Piv%BalanceFactor < -1) .OR. (Piv%BalanceFactor > 1)) then
     If (Compare(Item, Piv%Info) == -1) then
        If (Compare(Item, Piv%LeftChild%Info) == -1) then
           If (Associated(Piv, Root)) then
              Call LeftOfLeft(Root)
           else If (Associated(Piv, PivParent%LeftChild)) then
              Call LeftOfLeft(PivParent%LeftChild)
           else
              Call LeftOfLeft(PivParent%RightChild)
           End If
        else if (Associated(Piv, Root)) then
           Call RightOfLeft(Root)
        else If (Associated(Piv, PivParent%LeftChild)) then
           Call RightOfLeft(PivParent%LeftChild)
        else
           Call RightOfLeft(PivParent%RightChild)
        End If
     else If (Compare(Item, Piv%RightChild%Info) /= -1) then
        If (Associated(Root, Piv)) then
           Call RightOfRight(Root)
        else If (Associated(Piv, PivParent%LeftChild)) then
           Call RightOfRight(PivParent%LeftChild)
        else
           Call RightOfRight(PivParent%RightChild)
        End If
     else if (Associated(Piv, Root)) then
        Call LeftOfRight(Root)
     else if (Associated(Piv, PivParent%LeftChild)) then
        Call LeftOfRight(PivParent%LeftChild)
     else
        Call LeftOfRight(PivParent%RightChild)
     End If
  End If

  Return
End Subroutine


!******************************************************************************
!
!  CreateBtree - Creates an empty binary tree
!
!    Root - Binary tree root
!
!******************************************************************************

Subroutine CreateBtree(root)
  Type (BinaryTreeNode), Pointer     :: Root

  Nullify(root)

  Return
End Subroutine

End Module


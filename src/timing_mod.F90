!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

!******************************************************************************
!
! File : timing_mod.f90
!   by : Alan Tackett
!   on : 07/12/01
!  for : Socorro Timing routines
!
!  This module contains the stopwatch definitions and initializes them
!
!******************************************************************************

Module timing_mod

Use btree_support_mod
Use btree_mod
use kind_mod
Use mpi_mod
Use error_mod
Use io_mod
Use diary_mod
use utils_mod

implicit NONE!!!!!!!!
Private

Type (BinaryTreeNode), Pointer :: Root            !* Btree root node

!doc$
Public :: Print_ALL_Timers, InitTimers, Kill_All_Timers
Public :: Start_Timer, Stop_Timer, Pause_Timer, Resume_Timer, Reset_Timer

!**** Timer states *****
Integer, PRIVATE, PARAMETER :: STATE_STOPPED = 0
Integer, PRIVATE, PARAMETER :: STATE_RUNNING = 1
Integer, PRIVATE, PARAMETER :: STATE_PAUSED  = 2

!cod$

!******************************************************************************
Contains
!******************************************************************************

!******************************************************************************
! time_now - Returns the current time
!******************************************************************************

Function time_now() result(now)
  Real(double) :: now
  Real(double), external :: MPI_WTime
  
  now = MPI_WTime()
End Function

!******************************************************************************
! timer_uppercase - Just converts the timer name to uppercase
!******************************************************************************

Subroutine Timer_UpperCase(str)
  Character*(*), Intent(INOUT) :: str
  
  integer, PARAMETER :: a = iachar('a')
  Integer, PARAMETER :: z = iachar('z')
  Integer, PARAMETER :: shift = (iachar('A') - iachar('a'))
  
  integer i, lstr, c
  
  lstr = len(TRIM(str))
  
  do i = 1,lstr
     c = iachar(str(i:i))
     if ((a <= c) .and. (c <= z)) then
        str(i:i) = ACHAR(c + shift)
    else if (c == iachar('-')) then
       str(i:i) = "^"
    End If
  end do
  
end subroutine 
  
!******************************************************************************
!
!  WriteTimer - Does the actual writing of the timer info to the screen
!
!    Timer - Timer to print
!
!******************************************************************************

Subroutine WriteTimer(timer)

   type (timer_Type), Intent(IN) :: Timer

   integer :: ncalls, j, ii
   real(double) :: my_time, wall_time, cpu_time
   character(line_len) :: name
   character(line_len) :: fmt = '(a45,3x,i15,3x,2(f15.2,3x))'

   my_time = timer%elapsed
   cpu_time = 0.0_double
   call reduce(CONFIG,MPI_SUM, my_time, cpu_time)

   ncalls = 0
   j = timer%N  !* Don't want to mess up potential timer call
   call reduce(CONFIG,MPI_SUM, j, ncalls)

   wall_time = cpu_time/real(mpi_nprocs(CONFIG),double)

   if ( i_access( diaryfile() ) ) then
      name = timer%name
      do ii = 1, len(trimstr(name))
         if ( name(ii:ii) == "-" ) name(ii:ii) = " "
      end do
      write(x_unit(diaryfile()),fmt) name,ncalls,cpu_time,wall_time
   end if

End Subroutine

!******************************************************************************
!
!  PrintTimerInfo - Performs an in order traversal of the binary tree to print
!               the timer info
!
!     Root     - Binary tree root
!
!******************************************************************************

Recursive Subroutine PrintTimerInfo(Root)
  Type (BinaryTreeNode), Pointer  :: Root

  If (Associated(Root)) then
     Call PrintTimerInfo(Root%LeftChild)

     Call Writetimer(Root%Info%Timer)     

     Call PrintTimerInfo(Root%RightChild)
  EndIf

  Return
End Subroutine


!******************************************************************************
!
!  Print_ALL_timers - Prints the timing inforation to the specified unit
!      
!******************************************************************************

Subroutine Print_ALL_Timers()

   integer :: unit
   character(line_len) :: fmt = '(a45,3(3x,a15))'

   if ( i_access( diaryfile() ) ) then
      unit = x_unit(diaryfile())
      write(unit,'(/)')
      write(unit,fmt) "Timer Breakdown                              ",          "Calls",   "CPU Time (s)",  "Wall Time (s)"
      write(unit,fmt) "---------------------------------------------","---------------","---------------","---------------"
   end if

   call PrintTimerInfo(Root)
   call flush(unit)

End Subroutine

!******************************************************************************
!
!  Add_Timer - Adds a new timer to the list and expands it if needed.
!              The new timer index is returned.
!
!     Name - Name of new timer
!
!******************************************************************************

Function Add_timer(Name)  Result(Timer)
  Character*(*),     Intent(IN) :: Name
  Type (timer_type), Pointer    :: Timer

  Type (BinaryTreeData), Pointer     :: Item
  Logical                            :: Success


  Allocate(Timer)
  Call Create_timer(timer, Name)

  Allocate(Item)
  Item%Name = Name
  Call Timer_UpperCase(Item%name)  !*Ignore case in compares
  Item%timer => Timer
  Call InsertNode(Root, Item, Success)

  If (.NOT. Success) then
    Write(*,*) 'Add_timer: Error!!!!!!!!!!!!!!!!!!!!!!!!!'
    !* Should complain and exit here
  end If

  Return
End Function

!******************************************************************************
!
! Get_Timer - Finds the timer and optionally creates it if needed.
!             Returning the index in timers or -1 if not found or created.
!
!    Name -Name of the timer
!    CreateIfNeeded - Allows the timer to be created if it doesn't exist!
!        this is an optional parameter and is FALSE by default.
!
!******************************************************************************

Function Get_timer(Name, CreateIfNeeded)
  Character*(*),     Intent(IN) :: Name
  Logical, OPTIONAL, Intent(IN) :: CreateIfNeeded
  Type (timer_type), Pointer :: Get_Timer

  Type (BinaryTreeData) :: item
  Type (BinaryTreeNode), Pointer :: Result
  Logical :: success, DoCreate

  Item%Name = name
  Call Timer_UpperCase(Item%Name)
  Call SearchBTree(Root, Item, Result, Success)

  If (Success) then
     Get_Timer => Result%Info%Timer
  else
     DoCreate = .FALSE.
     If (Present(CreateIfNeeded)) DoCreate = CreateIfNeeded

     If (.NOT. Docreate) then
        Nullify(Get_timer)  !** Exit early since I'm not allowed to create it
        RETURN
     else               !** Create a new timer
         Get_Timer => Add_timer(Name)
     End If     
  End If

  Return
End Function

!******************************************************************************
!  Start_Timer - Starts the timer
!******************************************************************************

Subroutine Start_Timer(name)
  character*(*), Intent(IN) :: Name

  Type (timer_type), Pointer :: Timer
  Logical :: test
  
  timer => Get_timer(name, .TRUE.)
  
  test = timer%state /= STATE_STOPPED
  if (error(test, "timing::start_timer: ERROR: Timer NOT stopped! timer:" // TRIM(name))) goto 9

  timer%state = STATE_RUNNING
  timer%start_time = time_now()

9 continue

  Return
End Subroutine

!******************************************************************************
!  Stop_Timer - Stops the timer
!******************************************************************************

Subroutine Stop_Timer(name)
  character*(*), Intent(IN) :: Name

  Type (timer_type), Pointer :: Timer
  Logical :: test

  timer => Get_timer(name)

  test = .NOT. Associated(timer)
  if (error(test,"ERROR: unknown timer "//trim(name))) goto 9
  test = timer%state /= STATE_RUNNING
  if (error(test, "timing::stop_timer: ERROR: Timer not running! timer:" // TRIM(name))) goto 9

  timer%state = STATE_STOPPED
  timer%elapsed = timer%elapsed + time_now() - timer%start_time 
  timer%n = timer%n + 1

9 continue

  Return
End Subroutine


!******************************************************************************
!  Pause_Timer - Pauses the timer
!******************************************************************************

Subroutine Pause_Timer(Name)
  character*(*), Intent(IN) :: Name

  Type (timer_type), Pointer :: Timer
  Logical :: test

  timer => Get_timer(name)
   
  test = .NOT. Associated(timer)
  if (error(test,"ERROR: Unknown Timer "//TRIM(Name))) goto 9
  test = timer%state /= STATE_RUNNING
  if (error(test, "timing::pause_timer: ERROR: Timer not PAUSED! timer:" // TRIM(name))) goto 9

  timer%state = STATE_PAUSED
  timer%elapsed = timer%elapsed + time_now() - timer%start_time 

9 continue

  Return
End Subroutine

!******************************************************************************
!  Resume_Timer - Resumes the timer after a pause
!******************************************************************************

Subroutine Resume_Timer(Name)
  character*(*), Intent(IN) :: Name

  Type (timer_type), Pointer :: Timer
  logical :: not_associated
  logical :: not_paused
  
  timer => Get_timer(name)
   
  not_associated = .NOT. Associated(timer)
  if (error(not_associated,"timing::resume_timer: ERROR: Unknown Timer "//TRIM(Name))) goto 9
  not_paused = (timer%state == STATE_PAUSED) != STATE_PAUSED
  if (error(not_paused, "timing::resume_timer: ERROR: Timer not PAUSED! timer:" // TRIM(name))) goto 9

  timer%state = STATE_RUNNING
  timer%start_time = time_now()

9 continue

  Return
End Subroutine

!******************************************************************************
!  Reset_Timer - Resets the timer
!******************************************************************************

Subroutine Reset_Timer(Name)
  character*(*), Intent(IN) :: Name

  Type (timer_type), Pointer :: Timer
  Logical :: NotFound
  
  timer => Get_timer(name, .FALSE.)
   
  NotFound = .NOT. Associated(timer)
  if (error(NotFound,"ERROR: Unknown Timer "//TRIM(Name))) goto 9

  timer%N = 0
  timer%start_time = 0.0_double
  timer%state = STATE_STOPPED
  timer%elapsed = 0.0_double

9 if (error("Exit timing::Reset_Timer")) continue

  Return
End Subroutine

!******************************************************************************
!  Create_Timer - Create the timer
!******************************************************************************

Subroutine Create_Timer(Timer, timer_name)
  Type (Timer_Type), Intent(INOUT) :: Timer
  Character(len=*),  Intent(IN)    :: timer_name

  timer%Name = Timer_Name
  timer%N = 0
  timer%start_time = 0.0_double
  timer%state = STATE_STOPPED
  timer%elapsed = 0.0_double

  Return
End Subroutine

!******************************************************************************
!
!  InitTimers - Initializes all the timers for use.
!
!******************************************************************************

Subroutine InitTimers()

  Call CreateBtree(Root)

  Return
End Subroutine


!******************************************************************************
!
!  Kill_All_Timers - recursively deletes all the timers in the tree.
!
!******************************************************************************

Subroutine Kill_All_Timers()

  Call FreeTree(Root,.true.)

  Return
End Subroutine


End Module

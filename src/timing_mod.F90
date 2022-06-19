!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of contract DE-NA0003525     !
!  with NTESS, the United States Government retains certain rights to this software. This software is distributed uner the         !
!  modified Berkeley Software Distribution (BSD) License.                                                                          !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

module timing_mod

   use diary_mod
   use error_mod
   use io_mod
   use kind_mod
   use mpi_mod
   use utils_mod

   implicit none ; private

   !* Tags for tracking the state of a timer

   integer, parameter :: state_stopped = 0
   integer, parameter :: state_running = 1
   integer, parameter :: state_paused = 2

   !* Derived type that encapsulates a timer

   type stopwatch_obj

      private

      character(line_len) :: name
      integer :: state
      integer :: ncalls
      real(double) :: start_time
      real(double) :: elapsed_time

   end type stopwatch_obj

   !* Derived type that encapsulates the binary tree data

   type tree_data_obj

      private

      character(line_len)  :: name
      type(stopwatch_obj), pointer :: timer
      type(tree_data_obj), pointer :: level

   end type tree_data_obj

   !* Derived type that encapsulates the binary tree node

   type tree_node_obj

      private

      integer  :: balance_factor
      type(tree_data_obj), pointer :: info
      type(tree_node_obj), pointer :: left_child
      type(tree_node_obj), pointer :: right_child

   end type tree_node_obj

   !* Binary tree object containing the root node

   type(tree_node_obj), pointer :: root

   !* Publicly available parameters and procedures

   public :: timer_ctor
   public :: timer_dtor

   public :: start_timer
   public :: stop_timer
   public :: pause_timer
   public :: resume_timer
   public :: reset_timer
   public :: write_timers

   !* Interfaces for the publicly available procedures

   interface timer_ctor
      module procedure timer_ctor_
   end interface

   interface timer_dtor
      module procedure timer_dtor_
   end interface

   interface start_timer
      module procedure start_timer_
   end interface

   interface stop_timer
      module procedure stop_timer_
   end interface

   interface pause_timer
      module procedure pause_timer_
   end interface

   interface resume_timer
      module procedure resume_timer_
   end interface

   interface reset_timer
      module procedure reset_timer_
   end interface

   interface write_timers
      module procedure write_timers_
   end interface

contains


   !* Method to construct the timer system

   subroutine timer_ctor_()

      call binary_tree_ctor_(root)

   end subroutine timer_ctor_


   !* Method to destroy the timer system

   subroutine timer_dtor_()

      call binary_tree_dtor_(root,free_data=.true.)

   end subroutine timer_dtor_


   !* Method to start a given timer

   subroutine start_timer_( name )

      character(*) :: name

      logical :: test
      type(stopwatch_obj), pointer :: timer

      timer => get_timer_(name, create_timer = .true.)

      test = ( timer%state /= state_stopped )
      if ( error(test,"The timer '"//trimstr(name)//"' was not stopped before starting") ) goto 90

      timer%state = state_running
      timer%start_time = current_time_()

90    if ( error("Exiting timing_mod::start_timer_()") ) continue

   end subroutine start_timer_


   !* Method to stop a given timer

   subroutine stop_timer_( name )

      character(*) :: name

      logical :: test
      type(stopwatch_obj), pointer :: timer

      timer => get_timer_(name, create_timer = .false.)

      test = ( .not.associated( timer ) )
      if ( error(test,"Unknown timer "//trimstr(name)) ) goto 90

      test = ( timer%state /= state_running )
      if ( error(test,"The timer '"//trimstr(name)//"' was not running before stopping") ) goto 90

      timer%state = state_stopped
      timer%ncalls = timer%ncalls + 1
      timer%elapsed_time = timer%elapsed_time - timer%start_time + current_time_()

90    if ( error("Exiting timing_mod::stop_timer_()") ) continue

   end subroutine stop_timer_


   !* Method to pause a given timer

   subroutine pause_timer_( name )

      character*(*), intent(in) :: name

      logical :: test
      type(stopwatch_obj), pointer :: timer

      timer => get_timer_(name, create_timer = .false.)

      test = ( .not.associated( timer ) )
      if ( error(test,"Unknown timer "//trimstr(name)) ) goto 90

      test = ( timer%state /= state_running )
      if ( error(test,"The timer '"//trimstr(name)//"' was not running before being paused") ) goto 90

      timer%state = state_paused
      timer%elapsed_time = timer%elapsed_time - timer%start_time + current_time_()

90    if ( error("Exiting timing_mod::pause_timer_()") ) continue

   end subroutine


   !* Method to resume a given timer

   subroutine resume_timer_( name )

      character*(*), intent(in) :: name

      logical :: test
      type(stopwatch_obj), pointer :: timer

      timer => get_timer_(name, create_timer = .false.)

      test = ( .not.associated( timer ) )
      if ( error(test,"Unknown timer "//trimstr(name)) ) goto 90

      test = ( timer%state /= state_paused )
      if ( error(test,"The timer '"//trimstr(name)//"' was not paused before resuming") ) goto 90

      timer%state = state_running
      timer%elapsed_time = current_time_()

90    if ( error("Exiting timing_mod::resume_timer_()") ) continue

   end subroutine


   !* Method to reset a given timer

   subroutine reset_timer_( name )

      character*(*), intent(in) :: name

      logical :: test
      type(stopwatch_obj), pointer :: timer

      timer => get_timer_(name, create_timer = .false.)

      test = ( .not.associated( timer ) )
      if ( error(test,"Unknown timer "//trimstr(name)) ) goto 90

      timer%state = state_stopped
      timer%ncalls = 0
      timer%start_time = 0.0d0
      timer%elapsed_time = 0.0d0

90    if ( error("Exiting timing_mod::reset_timer_()") ) continue

   end subroutine


   !* Method to write all timer information to the diary file

   subroutine write_timers_()

      integer :: unit
      character(line_len) :: fmt = '(a45,3(3x,a15))'

      if ( i_access( diaryfile() ) ) then
         unit = x_unit(diaryfile())
         write(unit,'(/,"Runtime task breakdown:",/)')
         write(unit,fmt) "Timer                                        ",          "Calls",   "CPU Time (s)",  "Wall Time (s)"
         write(unit,fmt) "---------------------------------------------","---------------","---------------","---------------"
      end if

      call write_timer_tree_(root)
      call flush(unit)

   end subroutine write_timers_


   !* Method to retrieve a timer from the timer tree

   function get_timer_( name , create_timer ) result( timer )

      character(*), intent(in) :: name
      logical, intent(in), optional :: create_timer
      type(stopwatch_obj), pointer :: timer

      logical :: success, create
      type(tree_data_obj) :: item
      type(tree_node_obj), pointer :: result

      item%name = name
      call binary_tree_search_(root,item,result,success)

      if ( success ) then
         timer => result%info%timer
      else
         create = .false.
         if ( present( create_timer ) ) create = create_timer
         if ( .not.create ) then
            nullify( timer )
         else
            timer => add_timer_(name)
         end if
      end if

90    if ( error("Exiting timing_mod::get_timer_()") ) continue

   end function get_timer_


   !* Method to add a timer to the timer tree

   function add_timer_( name ) result( timer )

      character(*), intent(in) :: name
      type(stopwatch_obj), pointer :: timer

      logical :: success
      type(tree_data_obj), pointer :: item

      allocate( timer )

      timer%name = name
      timer%state = state_stopped
      timer%ncalls = 0
      timer%start_time = 0.0d0
      timer%elapsed_time = 0.0d0

      allocate( item )

      item%name = name
      item%timer => timer

      call binary_tree_insert_(root,item,success)

      if ( error((.not.success),"There was a problem inserting the '"//trimstr(name)//"' timer") ) continue
90    if ( error("Exiting timing_mod::add_timer_()") ) continue

   end function add_timer_


   !* Method to retrieve the current elapsed time

   function current_time_() result( now )

      real(double) :: now

      real(double) :: mpi_wtime

      now = mpi_wtime()

   end function current_time_


   !* Method to construct the binary tree

   subroutine binary_tree_ctor_( node )

      type(tree_node_obj), pointer :: node

      nullify( node )

   end subroutine binary_tree_ctor_


   !* Method to destroy the binary tree

   recursive subroutine binary_tree_dtor_( node , free_data )

      type(tree_node_obj), pointer :: node
      logical, intent(in) :: free_data

      if ( associated( node ) ) then
         call binary_tree_dtor_(node%left_child,free_data)
         call binary_tree_dtor_(node%right_child,free_data)
         if ( free_data ) then
            deallocate( node%info%timer )
            deallocate( node%info )
         end if
         deallocate( node )
      end if

   end subroutine binary_tree_dtor_


   !*

   subroutine binary_tree_search_( node , item , result , success )

      type(tree_node_obj), pointer :: node
      type(tree_data_obj), target :: item
      type(tree_node_obj), pointer :: result
      logical, intent(out) :: success

      integer :: comparison

      success = .false.

      if ( .not.associated( node ) ) then
         nullify( result )
         return
      end if

      result => node

      do while ( .not.success .and. associated( result ) )
         comparison = compare_tree_data_(item,result%info)
         select case ( comparison )
            case ( -1 )
               result => result%left_child
            case (  0 )
               success = .true.
            case (  1 )
               result => result%right_child
         end select
      end do

   end subroutine binary_tree_search_


   !*

   subroutine binary_tree_insert_( node , item , success )

      type(tree_node_obj), pointer :: node
      type(tree_data_obj), pointer :: item
      logical, intent(out) :: success

      integer  :: error
      type(tree_node_obj), pointer :: p, piv, piv_parent, in_p, in_parent, q

      success = .true.

      allocate( p )

      p%info => item
      p%balance_factor = 0

      nullify( p%left_child )
      nullify( p%right_child )

      if ( .not.associated( node ) ) then
         node => p
         return
      end if

      in_p => node
      piv => node

      nullify( in_parent )
      nullify( piv_parent )

      ! Search for an insertion point and pivot

      do while ( associated( in_p ) )
         if ( in_p%balance_factor /= 0 ) then
            piv => in_p
            piv_parent => in_parent
         end if
         in_parent => in_p
         if ( compare_tree_data_(item,in_p%info) == -1 ) then
            in_p => in_p%left_child
         else
            in_p => in_p%right_child
         end if
      end do

      ! Insert the node as a child of in_parent

      if ( compare_tree_data_(item,in_parent%info) == -1 ) then
         in_parent%left_child => p
      else
         in_parent%right_child => p
      end if

      ! Compute the balance factors between piv and in_parent

      q => piv

      if ( compare_tree_data_(item,q%info) == -1 ) then
         q%balance_factor = q%balance_factor + 1
         q => q%left_child
      else
         q%balance_factor = q%balance_factor - 1
         q => q%right_child
      end if

      do while ( .not.associated( q , p) )
         if ( compare_tree_data_(item, q%info) == -1 ) then
            q%balance_factor = q%balance_factor + 1
            q => q%left_child
         else
            q%balance_factor = q%balance_factor - 1
            q => q%right_child
         end if
      end do

      ! Determine of AVL rotations are needed

      if ( (piv%balance_factor < -1) .or. (piv%balance_factor > 1) ) then
         if ( compare_tree_data_(item,piv%info) == -1 ) then
            if ( compare_tree_data_(item, piv%left_child%info) == -1 ) then
               if ( associated( piv , node ) ) then
                  call left_of_left(node)
               else if ( associated(piv, piv_parent%left_child) ) then
                  call left_of_left(piv_parent%left_child)
               else
                  call left_of_left(piv_parent%right_child)
               end if
            else if ( associated(piv, node) ) then
               call right_of_left(node)
            else if ( associated(piv, piv_parent%left_child) ) then
               call right_of_left(piv_parent%left_child)
            else
               call right_of_left(piv_parent%right_child)
            end if
         else if ( compare_tree_data_(item, piv%right_child%info) /= -1 ) then
            if ( associated(node, piv) ) then
               call right_of_right(node)
            else if ( associated(piv, piv_parent%left_child) ) then
               call right_of_right(piv_parent%left_child)
            else
               call right_of_right(piv_parent%right_child)
            end if
         else if ( associated(piv, node) ) then
            call left_of_right(node)
         else if ( associated(piv, piv_parent%left_child) ) then
            call left_of_right(piv_parent%left_child)
         else
            call left_of_right(piv_parent%right_child)
         end if
      end if

   end subroutine binary_tree_insert_


   !*

   subroutine left_of_left(pivot)

      type(tree_node_obj), pointer :: pivot
      
      type(tree_node_obj), pointer :: p, q
      
      p => pivot%left_child
      q => p%right_child
      
      p%right_child => pivot
      pivot%left_child => q
      pivot => p
      
      pivot%balance_factor = 0
      pivot%right_child%balance_factor = 0
      
   end subroutine left_of_left


   !*

   subroutine right_of_right(pivot)
      
      type(tree_node_obj), pointer :: pivot
      
      type(tree_node_obj), pointer :: p, q
      
      p => pivot%right_child
      q => p%left_child
      
      p%left_child => pivot
      pivot%right_child => q
      pivot => p
      
      pivot%balance_factor = 0
      pivot%left_child%balance_factor = 0
      
   end subroutine right_of_right


   !*

   subroutine left_of_right(pivot)
      
      type(tree_node_obj), pointer :: pivot
      
      type(tree_node_obj), pointer :: x, y
      
      x => pivot%right_child
      y => x%left_child
      
      pivot%right_child => y%left_child
      x%left_child => y%right_child
      y%right_child => x
      y%left_child => pivot
      pivot => y
      
      select case ( pivot%balance_factor )
         case ( -1 )
            pivot%balance_factor = 0
            pivot%left_child%balance_factor = 1
            pivot%right_child%balance_factor = 0
         case (  0 )
            pivot%left_child%balance_factor = 0
            pivot%right_child%balance_factor = 0
         case default
            pivot%balance_factor = 0
            pivot%left_child%balance_factor = 0
            pivot%right_child%balance_factor = -1
      end select
      
   end subroutine left_of_right


   !*

   subroutine right_of_left(pivot)

      type(tree_node_obj), pointer :: pivot

      type(tree_node_obj), pointer :: x, y

      x => pivot%left_child
      y => x%right_child

      pivot%left_child => y%right_child
      x%right_child => y%left_child
      y%left_child => x
      y%right_child => pivot
      pivot => y

      select case ( pivot%balance_factor )
         case ( 0 )
            pivot%left_child%balance_factor = 0
            pivot%right_child%balance_factor = 0
         case ( 1 )
            pivot%balance_factor = 0
            pivot%left_child%balance_factor = 0
            pivot%right_child%balance_factor = -1
         case default
            pivot%balance_factor = 0
            pivot%left_child%balance_factor = 1
            pivot%right_child%balance_factor = 0
      end select

   end subroutine right_of_left


   !* 

   function compare_tree_data_( td1 , td2 ) result( val )

      type(tree_data_obj) :: td1, td2
      integer :: val

      if ( td1%name < td2%name ) then
         val = -1
      else if ( td1%name == td2%name ) then
         val = 0
      else
         val = 1
      end if

   end function compare_tree_data_


   !* Method to write the timer tree to the diary file

   recursive subroutine write_timer_tree_( node )

      type (tree_node_obj), pointer :: node

      if ( associated( node ) ) then
         call write_timer_tree_(node%right_child)
         call write_timer_data_(node%info%timer)
         call write_timer_tree_(node%left_child)
      end if

   end subroutine write_timer_tree_


   !* Method to write the timer data to the diary file

   subroutine write_timer_data_( timer )

      type(stopwatch_obj), intent(in) :: timer

      integer :: ii, ncalls
      real(double) :: cpu_time, wall_time
      character(line_len) :: name
      character(line_len) :: fmt = '(a45,3x,i15,3x,2(f15.2,3x))'

      ncalls = 0
      call reduce(CONFIG,MPI_SUM, timer%ncalls, ncalls)

      cpu_time = 0.0d0
      call reduce(CONFIG,MPI_SUM, timer%elapsed_time, cpu_time)
      wall_time = ( cpu_time / real(mpi_nprocs(config), double) )

      if ( i_access( diaryfile() ) ) then
         name = timer%name
         do ii = 1, len(trimstr(name))
            if ( name(ii:ii) == "-" ) name(ii:ii) = " "
         end do
         write(x_unit(diaryfile()),fmt) name,ncalls,cpu_time,wall_time
      end if

   end subroutine write_timer_data_


end module timing_mod

!* ------------------------------------------------------------------------------------------------------------------------------ *!
!  Socorro is a plane-wave density functional theory code for solid-state electronic structure calculations.                       !
!  See the README file in the top-level directory.                                                                                 !
!                                                                                                                                  !
!  Copyright 2011 National Technology and Engineering Solutions of Sandia, LLC (NTESS).                                            !
!  This software is distributed uner the modified Berkeley Software Distribution (BSD) License.                                    !
!  Under the terms of contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights to this software.               !
!* ------------------------------------------------------------------------------------------------------------------------------ *!

#include "macros.h"

      module tagio_mod
!doc$ module tagio_mod

!     One datatype is available here: type(tagio_obj).

!     Tagio_mod gives I/O abstractions to F90/95 processes.

      use path_mod
      use kind_mod
      use intstack_mod
      use io_mod

!cod$
      implicit none
      private

      type tag_obj
         integer(LONGLONG)  :: next_tag
         character(line_len) :: tag
         character(1)       :: tag_type
      end type

      type, public :: tagio_obj
         private
         type (file_obj)   :: f            !* File object
         character(line_len):: name         !* File name
         integer           :: ref          !* Ref count for SADR
         integer           :: mode         !* Access mode
         integer(CPTR)     :: fd           !* C file descriptor
         integer           :: filesex      !* Endian-nes of the file
         integer           :: cpusex       !* Endian-ness of processor
         integer(LONGLONG) :: fpos         !* Current File position
         integer(LONGLONG) :: curr_tag_start !* Current tag start position
         integer(LONGLONG) :: curr_tag_offset !* Current tag offset position
         integer(LONGLONG) :: curr_tag_end !* Current tag end position
         integer(LONGLONG) :: first_tag_pos!* First tag position
         integer(LONGLONG) :: data_start   !* beginning of file
         integer(LONGLONG) :: Bounds(3)    !* Block Bounds
         integer(LONGLONG) :: data_seg     !* Amount of data between tags
         type(tag_obj)     :: tag          !* Current tag
         type(intstack_obj), pointer :: head !* Tag stack
         Logical           :: BlockOpened
         Logical           :: End_Scope
      end type

!doc$
      character(11), public, parameter :: mkey = "Socorro_1.0"
      real(double), public, parameter  :: es_version = 2.0_double

      character(1), public, PARAMETER :: TAG_NOT_FOUND   = "U"
      character(1), public, PARAMETER :: TAG_NORMAL      = "N"
      character(1), public, PARAMETER :: TAG_START_BLOCK = "B"
      character(1), public, PARAMETER :: TAG_END_BLOCK   = "E"
      character(1), public, PARAMETER :: TAG_END_FILE    = "F"
      character(1), public, PARAMETER :: TAG_END_SCOPE   = "S"

      integer, public, parameter :: TAGIO_READ  = 1
      integer, public, parameter :: TAGIO_WRITE = 2

      integer, public, PARAMETER :: UNKNOWN_ENDIAN = 0
      integer, public, PARAMETER :: BIG_ENDIAN     = 1
      integer, public, PARAMETER :: LITTLE_ENDIAN  = 2
      integer, public, PARAMETER :: MIDDLE_ENDIAN  = 3

!cod$
      integer(CPTR),  external :: openf
      integer,        external :: closef
      integer,        external, public :: machinebyteorder
      integer(LONGLONG), external, public :: cfilesize

      External :: swapbytes;   public   :: swapbytes
      External :: getposf;     public   :: getposf
      External :: setposf;     public   :: setposf
      External :: readf;       public   :: readf
      External :: writef;      public   :: writef

      integer, PARAMETER :: B_START = 1
      integer, PARAMETER :: B_END   = 2
      integer, PARAMETER :: B_TAG_START   = 3

!doc$
      public :: i_access, i_comm, open_tagio, close_tagio
      public :: x_tagfd, x_datasize, x_end_scope, str_ftoc, str_ctof
      public :: x_swapbytes, x_noswap, x_filesex, x_cpusex
      public :: writetag, startblock, endblock
      public :: rewind_tofirsttag, rewind_tobeginning, rewind_toblockstart
      public :: findfirsttag, findnexttag, getnexttag
      public :: openblock, closeblock
      public :: tagio_readmagic, tagio_openerror
      public :: my, thy, glean, bequeath
      public :: tagio, x_ref

!cod$
      interface i_access
         module procedure tagio_i_access
      end interface
      interface i_comm
         module procedure tagio_i_comm
      end interface

      interface my
         module procedure my_tagio, my_new_tagio
      end interface
      interface thy
         module procedure thy_tagio
      end interface
      interface glean
         module procedure glean_tagio
      end interface
      interface bequeath
         module procedure bequeath_tagio
      end interface

      interface tagio
         module procedure open_tagio
      end interface

      interface x_ref
         module procedure tagio_ref
      end interface

      contains

!******************************************************************        

      subroutine my_tagio(tagio)
!doc$ subroutine my(tagio)
        type (tagio_obj) :: tagio

!cod$
        tagio%ref = tagio%ref + 1

      end subroutine

!******************************************************************        

      subroutine my_new_tagio(tin,tout)
!doc$ subroutine my(tin,tout)
        type (tagio_obj), Intent(IN)  :: tin
        type (tagio_obj)              :: tout

!cod$
        tout = tin
        tout%ref = tout%ref + 1

      end subroutine

!******************************************************************        

      function thy_tagio(tagio) result(thy_obj)
!doc$ function thy(tagio) result(thy_obj)
        type (tagio_obj), Intent(INOUT)  :: tagio
        type (tagio_obj)                 :: thy_obj

!cod$
        tagio%ref = tagio%ref - 1
        thy_obj = tagio        

      end function

!******************************************************************        

      subroutine glean_tagio(tagio)
!doc$ subroutine glean(tagio)
        type (tagio_obj) :: tagio

!cod$
        if (tagio%ref < 1) call close_tagio(tagio)

      end subroutine

!******************************************************************        

      subroutine bequeath_tagio(tagio)
!doc$ subroutine bequeath(tagio)
        type (tagio_obj), Intent(INOUT) :: tagio

!cod$
        continue

      end subroutine

!******************************************************************        

      Subroutine tagio_readmagic(name, mkey, klen) 
        character(*), intent(in)  :: name
        character(*), intent(out) :: mkey
        Integer,      intent(out) :: klen

        character(line_len)   :: fname
        character(1)         :: c
        integer(long)        :: datasize, numdata  
        integer(long)              :: err
        type(tagio_obj)      :: tagio
        
        call str_ftoc("r", mkey)
        call str_ftoc(name,fname)
        tagio%fd = openf(fname, mkey)

        numdata = 1; datasize = 1;
!        Call readf(c, datasize, numdata, tagio%fd, x_noswap(tagio), err)
        Call readf(c, datasize, numdata, tagio%fd, x_noswap(), err)
        klen = ichar(c)
        
        numdata = klen; datasize = 1;
!        Call readf(mkey, datasize, numdata, tagio%fd, x_noswap(tagio), err)
        Call readf(mkey, datasize, numdata, tagio%fd, x_noswap(), err)
        
      end subroutine

!******************************************************************        

      function open_tagio(name, mode, magic_key, key_len) result(tagio)
        type (tagio_obj)          :: tagio
        character(*), intent(in) :: name
        integer,       intent(in) :: mode
        character(*), intent(in) :: magic_key
        integer,       intent(in) :: key_len
        
        character(line_len)   :: mkey, fname
        character(1)         :: c
        integer(long)        :: datasize, numdata, err
        integer              :: j
        type (tag_obj)       :: tag
        logical              :: e_ok

        tagio%fpos = -1
        tagio%curr_tag_start = -1
        tagio%curr_tag_end   = -1
        tagio%curr_tag_offset = -1
        tagio%bounds = 0
        tagio%name = name
        tagio%ref = 0

        Call createstack(tagio%head)

        tagio%mode = mode
        Call my(file(name),tagio%f)
        tagio%cpusex = machinebyteorder()
!write(*,*) 'cpusex=',tagio%cpusex

        !** Exit early if I'm not supposedto access anything
        if (.NOT. i_access(tagio%f)) then
!Write(*,*) 'open_tagio:  Exiting early'
           RETURN
        End If

!Write(*,*) 'open_tagio: Opening file'
        if (mode == TAGIO_READ) then
           Inquire(FILE=trim(name), EXIST=e_ok)
           if (.NOT.  e_ok) then
              tagio%fd = 0
              RETURN
           end if

           call str_ftoc("r", mkey)
           call str_ftoc(name,fname)
           tagio%fd = openf(fname, mkey)
!write(*,*) 'FD=',tagio%fd 

           numdata = 1; datasize = 1;
!           Call readf(c, datasize, numdata, tagio%fd, x_noswap(tagio), err)
           Call readf(c, datasize, numdata, tagio%fd, x_noswap(), err)
           j = ichar(c)

           numdata = key_len; datasize = 1;
!           Call readf(mkey, datasize, numdata, tagio%fd, x_noswap(tagio), err)
           Call readf(mkey, datasize, numdata, tagio%fd, x_noswap(), err)

           if (mkey(1:key_len) /= magic_key(1:key_len)) then
              write(*,*) 'tagio: Magic key wrong!'
              RETURN
           end if

           !** Read the file byte ordering
           numdata = 1; datasize = 1;
           Call readf(c, datasize, numdata, tagio%fd, x_swapbytes(tagio), err)
           tagio%filesex = ichar(c)

           !** Read the 1st tag for setting the position
           call readtag(tagio, tag)
!write(*,*) 'open_tagio: tag=',trim(tag%tag)

           tagio%first_tag_pos = tag%next_tag
           Call getposf(tagio%fd, tagio%data_start)

           tagio%bounds(B_START) = tagio%data_start
           call str_ftoc(x_name(tagio%f), fname)
           tagio%bounds(B_END) = cfilesize(fname)
           tagio%bounds(B_TAG_START) = tagio%first_tag_pos

!write(*,*) 'FD=',tagio%fd !**, ' * DS=',tagio%data_start, ' * FT=',tagio%first_tag_pos
        else     !** Init the file by writing the magic word
           call str_ftoc("w+", mkey)
           call str_ftoc(name,fname)
           tagio%fd = openf(fname, mkey)

!write(*,*) 'FD=',tagio%fd
           numdata = 1; datasize = 1;
           c = achar(key_len)
           Call  writef(c, datasize, numdata, tagio%fd, err)

           numdata = key_len; datasize = 1;
           Call writef(magic_key, datasize, numdata, tagio%fd, err)

           !** Set and store the byte ordering
           tagio%filesex = tagio%cpusex
           numdata = 1; datasize = 1;
           c = achar(tagio%cpusex)
           Call writef(c, datasize, numdata, tagio%fd, err)

           !** Lastly store the 1st tag
           call writetag(tagio, "BASETAG!!!!!!")           
        end if

        tagio%fpos = key_len

        return
      end function 


!******************************************************************        

      subroutine close_tagio(tagio)
        type (tagio_obj) :: tagio

        integer :: result

        call deletestack(tagio%head)

        if (i_access(tagio%f)) result = closef(tagio%fd)

	Call glean(thy(tagio%f))

        tagio%fd = 0
      end subroutine

!******************************************************************        

      integer(CPTR) function x_tagfd(tagio)
        type (tagio_obj) :: tagio

        x_tagfd = tagio%fd
      end function

!******************************************************************        

!      integer(CPTR) function x_noswap(tagio)
!        type (tagio_obj) :: tagio

      integer(CPTR) function x_noswap()

        x_noswap = -1

      end function

!******************************************************************        

      integer(CPTR) function x_swapbytes(tagio)
        type (tagio_obj) :: tagio

        if (tagio%cpusex /= tagio%filesex) then
           x_swapbytes = 1
        else
           x_swapbytes = -1
        end If

      end function

!******************************************************************        

      integer function x_filesex(tagio)
        type (tagio_obj) :: tagio

        x_filesex = tagio%filesex

      end function

!******************************************************************        

      integer function x_cpusex(tagio)
        type (tagio_obj) :: tagio

        x_cpusex = tagio%cpusex

      end function
      
!******************************************************************        

      logical function tagio_openerror(tagio)
        type (tagio_obj) :: tagio

        if (tagio%fd == 0) then
           tagio_openerror = .TRUE.
        else
           tagio_openerror = .FALSE.
        end if

      end function tagio_openerror
!******************************************************************        

!doc$ logical function i_access(tagio)
!doc$ type(tagio_obj) :: tagio
!doc$ effects : returns (not x_first_only) or (mpi_isroot())
!doc$           to determine if the process should call read or write
!******************************************************************

      logical function tagio_i_access(tagio)
        type(tagio_obj) :: tagio

        tagio_i_access = i_access(tagio%f)
      end function 

!******************************************************************
!doc$ logical function i_comm(tagio)
!doc$ type(tagio_obj) :: tagio
!doc$ effects : returns x_first_only
!doc$           to determine if communication is required
!doc$           before writes or after reads.
!******************************************************************

      logical function tagio_i_comm(tagio)
        type(tagio_obj) :: tagio

        tagio_i_comm = i_comm(tagio%f)
      end function


!**********************************************************

      integer(LONGLONG) function x_datasize(tagio)
        type (tagio_obj) :: tagio

        x_datasize = tagio%data_seg

        return
      end function x_datasize

!**********************************************************

      Logical function x_end_scope(tagio)
        type (tagio_obj) :: tagio

        x_end_scope = tagio%end_scope

        return
      end function

!**********************************************************

      subroutine str_ftoc(fstr, cstr)
        character(*), intent(IN) :: fstr
        character(*), intent(OUT) :: cstr
        
        integer :: i
        
        cstr = fstr
        i = len(trim(fstr))+1
        cstr(i:i) = achar(0)          
        
        return
      end subroutine str_ftoc

!**********************************************************
      
      subroutine str_ctof(cstr, fstr)
        character(*), intent(IN) :: cstr
        character(*), intent(OUT) :: fstr
        
        integer :: i
        character :: null
        
        null = achar(0)
        
        i = index(cstr, null)
        if (i==0) then
           fstr = cstr
        else
           fstr = cstr(1:i-1)
        end if
        
        return
      end subroutine str_ctof

!**********************************************************

      subroutine writetag(tagio, tag, blockinfo)
        type(tagio_obj) :: tagio
        character(*),    intent(in)    :: tag
        character(1), optional, intent(in)  :: blockinfo

        integer(longlong) :: offset
        integer(long)     :: datasize, ndata, i
        character(1)      :: Block
        character(line_len):: ctag

        !** Get the current position
        Call getposf(tagio%fd, tagio%curr_tag_start)

        !** write the block info
        Block = TAG_NORMAL
        if (present(blockinfo)) then
           Block = blockinfo
        end if
        datasize = 1; ndata = 1;
        Call writef(block, datasize, ndata, tagio%fd, i)
        
        !** write the tag
        ndata = len(trim(tag)) + 1;  datasize = 1;
        call str_ftoc(tag, ctag)
        call writef(ctag, datasize, ndata, tagio%fd, i)

        !** If needed update the previous tag's link **
        if (tagio%curr_tag_offset > 0) then
           Call endtag(tagio)
        End If
           

        !** store the offset placeholder for the current tag
        Call getposf(tagio%fd, tagio%curr_tag_offset)
        ndata = 1; datasize = sizeof_longlong
        offset = 0
        Call writef(offset, datasize, ndata, tagio%fd, i)

!        !** Store the start of the tag
!        tagio%curr_tag_start = currpos
!write(*,*) 'writetag: tstart=',tagio%curr_tag_start


      end subroutine

!**********************************************************

      subroutine endtag(tagio)
        type(tagio_obj) :: tagio

        integer(longlong) :: currpos
        integer(long)     :: i, ndata, dsize

        if (tagio%curr_tag_offset < 0) then
           write(*,*) 'closetag: No tag opened!!!!!!!!!!!!!'
           STOP
        end if

        !** Get the current position
        Call  getposf(tagio%fd, currpos)

        !** Now move back to the offset plac
        !** Now move back to the offset placeholder and eholder and 
        !** store the offset
        Call setposf(tagio%fd, tagio%curr_tag_offset, i)

        ndata=1; dsize = sizeof_longlong
        Call writef(tagio%curr_tag_start, dsize, ndata, tagio%fd, i)

        !** finally move back to the end
        Call setposf(tagio%fd, currpos, i)

        !** and clear the tag info
        tagio%curr_tag_offset = -1
        tagio%curr_tag_start = -1

      end subroutine

!**********************************************************

      subroutine startblock(tagio, tag)
        type(tagio_obj) :: tagio
        character(*),    intent(in)    :: tag

        call writetag(tagio, tag, TAG_START_BLOCK)

      end subroutine

!**********************************************************

      subroutine endblock(tagio)
        type(tagio_obj) :: tagio


        character(20)   :: tag

        tag = "BLOCKEND"
        call writetag(tagio, tag, TAG_END_BLOCK)
        
        !** reset this since we don't really have a new tag.
        tagio%curr_tag_start = -1  

      end subroutine

!**********************************************************

!**********************************************************

      subroutine rewind_tofirsttag(tagio)
        type(tagio_obj) :: tagio

        Integer(LONG) :: i
        Character(line_len) :: fname

        tagio%bounds(B_START) = tagio%data_start
        call str_ftoc(x_name(tagio%f), fname)
        tagio%bounds(B_END) = cfilesize(fname)
        tagio%bounds(B_TAG_START) = tagio%first_tag_pos

        call deletestack(tagio%head)
        call createstack(tagio%head)

        tagio%curr_tag_start = tagio%first_tag_pos
        tagio%tag%tag_type = TAG_NORMAL
        tagio%tag%next_tag = tagio%first_tag_pos

        Call setposf(tagio%fd, tagio%first_tag_pos, i)
      end subroutine

!**********************************************************

      subroutine rewind_tobeginning(tagio)
        type(tagio_obj) :: tagio

        integer(long) :: i
        Character(line_len) :: fname

        tagio%bounds(B_START) = tagio%data_start
        call str_ftoc(x_name(tagio%f), fname)
        tagio%bounds(B_END) = cfilesize(fname)
        tagio%bounds(B_TAG_START) = tagio%first_tag_pos

        call deletestack(tagio%head)
        call createstack(tagio%head)

        tagio%curr_tag_start = tagio%first_tag_pos
        tagio%tag%tag_type = TAG_NORMAL
        tagio%tag%next_tag = tagio%first_tag_pos

        Call setposf(tagio%fd, tagio%data_start, i)
      end subroutine

!**********************************************************

      subroutine rewind_toblockstart(tagio)
        type(tagio_obj) :: tagio

        integer(long) :: i

!write(*,*) 'rewind_toblockstart: BS=',tagio%bounds(B_START)
        Call setposf(tagio%fd, tagio%bounds(B_START), i)
        tagio%curr_tag_start = tagio%bounds(B_TAG_START)
        tagio%tag%next_tag = tagio%bounds(B_TAG_START)

      end subroutine


!**********************************************************

      subroutine readtag(tagio, tag)
        type(tagio_obj) :: tagio
        type (tag_obj) :: tag

        integer(long)     :: i, ndata, dsize, result
        character(1)      :: c

        tag%tag = ""

        Call getposf(tagio%fd, tagio%curr_tag_start)

        ndata = 1;  dsize = 1
!        Call readf(tag%tag_type, dsize, ndata, tagio%fd, x_noswap(tagio), result)
        Call readf(tag%tag_type, dsize, ndata, tagio%fd, x_noswap(), result)
!write(*,*) 'readtag: tstart=',tagio%curr_tag_start, ' * type=',tag%tag_type
    
        i = 1;
!        Call readf(c, dsize, ndata, tagio%fd, x_noswap(tagio), result)
        Call readf(c, dsize, ndata, tagio%fd, x_noswap(), result)
        do while (ichar(c) /= 0)
!write(*,*) 'readtag: c=',c
           if (i <= len(tag%tag)) tag%tag(i:i) = c
           i = i + 1;
!           Call readf(c, dsize, ndata, tagio%fd, x_noswap(tagio), result)
           Call readf(c, dsize, ndata, tagio%fd, x_noswap(), result)
        end do

        ndata = 1;  dsize = sizeof_longlong
        Call readf(tag%next_tag, dsize, ndata, tagio%fd, x_swapbytes(tagio), result)

!write(*,*) 'readtag: tag=',trim(tag%tag), ' * type=',tag%tag_type, ' * next=',tag%next_tag, ' * swap=',x_swapbytes(tagio)
        Call getposf(tagio%fd, tagio%curr_tag_end)

        tagio%data_seg = tag%next_tag - tagio%curr_tag_end

        tagio%tag = tag

        tagio%BlockOpened = .FALSE.

        tagio%end_scope = .FALSE.
        if ((tag%next_tag > tagio%bounds(B_END)) .OR.  &
            (tag%next_tag == 0)) tagio%end_scope = .TRUE.

      end subroutine readtag


!**********************************************************

      subroutine internal_nexttag(tagio, tag)
        type(tagio_obj) :: tagio
        type (tag_obj) :: tag
        
        integer(long) :: err

!write(*,*) 'internal_nexttag: nexT_tag = ',tagio%tag%next_tag

        if (tagio%tag%next_tag /= 0) then
           Call setposf(tagio%fd, tagio%tag%next_tag, err)
           
           call readtag(tagio, tag)
        else
           tag%tag_type = TAG_END_FILE
           tag%next_tag = 0
        end if

      end subroutine

!**********************************************************

      subroutine skipblock(tagio)
        type(tagio_obj) :: tagio

        type (tag_obj) :: tag
        type (intstack_obj), pointer :: head
        Logical :: ok, pop_ok
        integer(longlong) :: pos

!write(*,*) 'skipblock: start'
        call createstack(head)

        ok = .FALSE.
        do while (.NOT. Ok)
           call internal_nexttag(tagio, tag)

           select case (tag%tag_type)
              case (TAG_START_BLOCK)
                 call getposf(tagio%fd, pos)
                 call push(head, pos)
              case (TAG_END_BLOCK)
                 pos = pop(head, pop_Ok)
                 Ok = .NOT. pop_ok
              case (TAG_END_FILE)
                 Ok = .TRUE.
           end select
        end do

        call deletestack(head)
!write(*,*) 'skipblock: end'
      end subroutine skipblock

!**********************************************************

      character(1) function scantag(tagio, key, movetostart)
        type(tagio_obj) :: tagio
        character(*),    intent(IN)    :: key
        logical,         intent(in)    :: movetostart

        integer(long) :: result
        type (tag_obj) :: tag

        if (movetostart) then
!write(*,*) 'scantag: Bpos=',tagio%bounds(B_TAG_START)
           Call setposf(tagio%fd, tagio%bounds(B_TAG_START), result)
        else
!write(*,*) 'scantag: Cpos=',tagio%curr_tag_start
           Call setposf(tagio%fd, tagio%tag%next_tag, result)
        end if

        call readtag(tagio, tag)        
        do while ((trim(tag%tag) /= trim(key)) .and. (tag%next_tag /= 0) &
             .and. (tag%next_tag <= (tagio%bounds(B_END)-1)))
!write(*,*) 'next1=',tag%next_tag
           if (tag%tag_type == TAG_START_BLOCK) Call skipblock(tagio)
!write(*,*) 'next2=',tag%next_tag
           call internal_nexttag(tagio, tag)
!write(*,*) 'trim test=', (trim(tag%tag) /= trim(key))
!write(*,*) 'next_tag test=',(tag%next_tag /= 0) , ' * nt=',tag%next_tag
!write(*,*) 'bounds=', (tag%next_tag <= (tagio%bounds(B_END)-1))

        end do

!write(*,*) 'scantag: tag=',tag%tag
        scantag = TAG_NOT_FOUND
        if (trim(tag%tag) == trim(key)) scantag = tag%tag_type

      end function 

!**********************************************************

      character(1) function findfirsttag(tagio, key)
        type(tagio_obj) :: tagio
        character(*),    intent(IN)    :: key

        findfirsttag = scantag(tagio, key, .TRUE.)

      end function findfirsttag

!**********************************************************

      character(1) function findnexttag(tagio, key)
        type(tagio_obj) :: tagio
        character(*),    intent(IN)    :: key

        findnexttag = scantag(tagio, key, .FALSE.)

      end function findnexttag

!**********************************************************

      character(1) function getnexttag(tagio, key)
        type(tagio_obj) :: tagio
        character(*),    intent(out)    :: key

        type (tag_obj) :: tag

        if (tagio%curr_tag_start >= tagio%bounds(B_END)) then 
           key = "UNKNOWN"
           getnexttag = TAG_END_BLOCK
           tagio%end_scope = .TRUE.
           RETURN
        elseif (tagio%tag%next_tag <= 0) then
           key = "UNKNOWN"
           getnexttag = TAG_END_FILE
           tagio%end_scope = .TRUE.
           RETURN
        end if

        if (tagio%tag%tag_type == TAG_START_BLOCK) then
           if (.NOT. tagio%BlockOpened) then
              Call Skipblock(tagio)
           end If
        End IF

        call internal_nexttag(tagio, tag)
        
           
        key = tag%tag

        getnexttag = tag%tag_type

      end function

!**********************************************************

      subroutine openblock(tagio)
        type(tagio_obj) :: tagio

        integer(longlong) :: currpos
        type (tag_obj)    :: tag
        integer(long)     :: err

!write(*,*) 'openblock:  Old Bounds=',tagio%bounds
        call push(tagio%head, tagio%bounds(B_START))
        call push(tagio%head, tagio%bounds(B_END))
        call push(tagio%head, tagio%bounds(B_TAG_START))

        currpos = tagio%curr_tag_start

        tagio%bounds(B_START) = tagio%curr_tag_end
        tagio%bounds(B_TAG_START) = tagio%tag%next_tag
        call skipblock(tagio)
        tagio%bounds(B_END) = tagio%curr_tag_start

        Call setposf(tagio%fd, currpos, err)
        Call readtag(tagio, tag)

        tagio%tag%next_tag = tagio%bounds(B_TAG_START)

!write(*,*) 'openblock:  Bounds=',tagio%bounds

        tagio%BlockOpened = .TRUE.

      end subroutine

!**********************************************************

      subroutine closeblock(tagio)
        type(tagio_obj) :: tagio

        type(tag_obj) :: tag
        integer(long) :: result
        logical       :: Ok

        if (.NOT. Associated(tagio%head)) then
           write(*,*) 'closeblock:  Error!!! No block opened!'
           STOP
        end if

        Call setposf(tagio%fd, tagio%bounds(B_END), result)

        call readtag(tagio, tag)

        tagio%bounds(B_TAG_START) = pop(tagio%head, Ok)
        tagio%bounds(B_END) = pop(tagio%head, Ok)
        tagio%bounds(B_START) = pop(tagio%head, Ok)
        tagio%end_scope = .FALSE.

!write(*,*) 'closeblock: tag=',trim(tag%tag), ' * Bnds=',tagio%bounds

      end subroutine

!**********************************************************

      function tagio_ref(tagio) result(r)
!doc$ function x_ref(tagio) result(r)
        type(tagio_obj) :: tagio
        integer :: r
!       effects: Returns tagio%ref

!cod$
        r = tagio%ref
        call glean(tagio)

      end function

end module

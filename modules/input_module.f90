!
! ParaGauss, a program package for high-performance computations
! of molecular systems
! Copyright (C) 2014
! T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
! M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
! A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
! T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
! M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
! M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License version 2 as published
! by the Free Software Foundation [1].
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
!===============================================================
! Public interface of module
!===============================================================
module input_module
!---------------------------------------------------------------
!
!  Purpose: performs input preprocessing. read routine should
!           access the input file via this modules by calling
!           the subroutine input_read_to_intermediate() which
!           writes one line of preprocessed input as a single
!           record to an intermediate file. To access this
!           data use the function input_intermediate_unit()
!           which return the unit number of the intermediate
!           file. Opening and closing is done automatically.
!
!           Counts line numbers.
!
!           Offers the subroutine input_error(message) to print
!           errormessage showing line number file name and
!           input line.
!
!           Offers the subroutines input_open(filename),
!           input_close() and input_close_one() to open and
!           close input file.
!           Input_open() may be called an arbitrary number of
!           times (as long as units are available) to open new
!           files while other files are still open. The states
!           of theses files are saved an the last opened file
!           is read till the end. After that, reading is continued
!           in the file open before. input_close() closes all open
!           files and input_close_one() only the latest file.
!
!           Does the following Preprocessing to input file:
!             + removing everything after "#" in line, allowing
!               comments beginning with this letter.
!             + removes trailing and preceeding blanks
!             + removes empty lines
!             + allows continuing lines: the preceeding line must
!               contain "%". Everything ater this character is ignored
!               and the next line is read.
!             + "#", "%" and \ can be masked with an \ preceding it.
!             + allows to include recursively other files:
!               "~" at the beginning of line followed by filename
!
!
!
!  Module called by: everything reading input
!             the read routines of the various modules are called in
!             subroutinr read_input
!
!
!
!  Author: TB
!  Date: 10/95
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: UB
! Date:   7/97
! Description: a) Modified to allow for arbitrary length input lines.
!              b) Preceeding blanks in a line are generally ignored now.
!
! Modification (Please copy before editing)
! Author: AS
! Date:   9/2003
! Description: new procedure to read in data from single intermediate
!              input file was included
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!----------------------------------------------------------------
#include "def.h"
use type_module ! type specification parameters
use iounitadmin_module, only : get_iounit,openget_iounit,returnclose_iounit ! to get input unit
use filename_module, only: strlen=>filename_namelengthmax
implicit none
save
private         ! by default, all names are private

!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------

public :: input_open
public :: input_close
public :: input_intermediate_unit
public :: input_read_to_intermediate
public :: input_which_namelist
public :: input_error
public :: input_line_is_namelist
public :: input_end_of_file

public :: truefalse ! truefalse('True'/'False') retruns true/false



!================================================================
! End of public interface of module
!================================================================

!------------ Declaration of constants and variables ----

type store_type
   integer :: unit, line_nbr
   character(len=strlen) :: filename
   type(store_type), pointer :: last
end type store_type

integer(kind=i4_kind), parameter :: block_length = 132
character(len=7)     , parameter :: format = "(132A1)"
type(store_type), pointer            :: store, first_store
integer(kind=i4_kind)                :: intermediate_unit
integer(kind=i4_kind)                :: new_inp_unit
character(len=1), pointer            :: raw_line(:)
character(len=1), pointer            :: processed_line(:)
integer(kind=i4_kind)                :: i_processed
integer(kind=i4_kind)                :: length_raw
logical :: intermediate_open, line_read, store_associated=.false., &
     question_mode, is_namelist, namelist_open
integer, parameter                   :: namelistname_max_length=32
character(len=namelistname_max_length) :: namelistname



!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

   !*************************************************************
   integer function input_intermediate_unit()
   !  Purpose: returns unit number of intermediate file and opens
   !           file if necessary
   !** End of interface *****************************************
#ifndef WITH_OLD_INPUT
      input_intermediate_unit = new_inp_unit
#else
      if ( intermediate_unit .eq. 0 ) intermediate_unit = get_iounit()
      input_intermediate_unit = intermediate_unit
#endif
   end function input_intermediate_unit
   !*************************************************************

   !*************************************************************
   subroutine input_read_to_intermediate()
   !  Purpose: reads input and writes one preprocessed line
   !  to intermediate file, rewinding it. Opens intermediate
   !  file if necessary.
   !** End of interface *****************************************
   use filename_module, only: tmpfile
   implicit none
#ifndef WITH_OLD_INPUT
   ! do nothing:
   RETURN
#else
   integer(kind=i4_kind) :: i

   DPRINT 'input_read_to_intermediate: entered'
   if ( intermediate_unit .eq. 0 ) intermediate_unit = get_iounit()
   if (.not. intermediate_open ) then
      DPRINT 'input_read_to_intermediate: open ',trim(tmpfile("input_intermediate"))
      open(intermediate_unit, err=777, &
          file = trim(tmpfile("input_intermediate")), &
          action = "readwrite", status = "replace", recl = 65536 )
   else
      DPRINT 'input_read_to_intermediate: rewind(intermediate_unit, err=779)'
      rewind(intermediate_unit, err=779)
   endif
   if ( .not. line_read ) call input_readline()
!:VPP_BUG[
!  ! The trailing blank in processed_line is skipped here
!  do i=1,i_processed-1,block_length
!     if (i_processed > i+block_length) then
!        write(intermediate_unit,format,advance='NO',err=778) &
!            processed_line(i:i+block_length-1)
!     else
!        write(intermediate_unit,format,advance='YES',err=778) &
!             processed_line(i:i_processed-1)
!     endif
!  enddo
   ! Due to a compiler bug on VPP the entire processed line must be written
   do i=1,size(processed_line),block_length
      if (size(processed_line) > i+block_length-1) then
         DPRINT 'input_read_to_intermediate: WRITE>>',processed_line(i:i+block_length-1),'<<'
         write(intermediate_unit,format,advance='NO',err=778) &
             processed_line(i:i+block_length-1)
      else
         DPRINT 'input_read_to_intermediate: WRITE>>',processed_line(i:),'<<'
         write(intermediate_unit,format,advance='YES',err=778) &
              processed_line(i:)
      endif
   enddo
!:VPP_BUG]
   DPRINT 'input_read_to_intermediate: rewind intermediate_unit'
   rewind(intermediate_unit, err=779)
   line_read = .false.
   intermediate_open = .true.
   return
777 call error_handler("input_read_to_intermediate: open failed")
778 call error_handler("input_read_to_intermediate: writing failed")
779 call error_handler("input_read_to_intermediate: rewind failed")
#endif
   end subroutine input_read_to_intermediate
   !*************************************************************


   !*************************************************************
   integer function input_linenbr()
   !  Purpose: returns actual line number in input file
   !** End of interface *****************************************
   input_linenbr = store%line_nbr
   end function input_linenbr
   !*************************************************************


   !*************************************************************
   logical function input_end_of_file()
   !  Purpose: .TRUE. if either no input file has been opened so
   !                  far or the EOF mark of the input is reached
   !  UB /9/97
   !** End of interface *****************************************

#ifndef WITH_OLD_INPUT
   character(len=300) :: line

!!$      input_end_of_file=eof(new_inp_unit)
      read(new_inp_unit,*,end=100) line
      input_end_of_file = .false.
      backspace new_inp_unit
      return
100   input_end_of_file = .true.
#else
      if ( .not. line_read ) then
         if ( store_associated ) then
            question_mode = .true.
            call input_readline()
            question_mode = .false.
         endif
      endif
      input_end_of_file = .not. line_read
#endif
   end function input_end_of_file
   !*************************************************************

   !*************************************************************
   function string(array)
     !  Purpose: converts an array of type character(len=1) into
     !           a single string of thype cvharacter(len=*)
     !------------ Declaration of formal parameters ---------------
     character(len=1), intent(in) :: array(:)
     character(len=size(array))   :: string   ! automatic variable
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     integer(kind=i4_kind)     :: pos ! , indent, length
     !------------ Executable Statements --------------------------
     do pos=1,size(array)
        string(pos:pos) = array(pos)
     enddo
   end function string
   !*************************************************************

   !*************************************************************
   logical function input_line_is_namelist(string)
   !
   !  Purpose: returns true if input line is
   !           namelist string
   !
   !------------ Declaration of formal parameters ---------------
   character(len=*), intent(in) :: string
   !** End of interface *****************************************

   character(len=LEN(string)) :: string_lower
   character(len=LEN(string)) :: namelist_name
   logical :: is_namelist

   string_lower = string

   call lowcase(string_lower)

   !
   ! Returns current namelist name in lower case, or a failure:
   !
   is_namelist = input_which_namelist(namelist_name)

   input_line_is_namelist = is_namelist .and. namelist_name == string_lower
   end function input_line_is_namelist
   !*************************************************************


   !*************************************************************
   logical function input_which_namelist(namelist_name)
   !  Purpose: gives name ot current input namelist.
   !  Returns .true. if current input line is namelist
   !  and .false. if it is no namelist or end of file is reached.
   !------------ Declaration of formal parameters ---------------
   character(len=*), intent(out) :: namelist_name ! in lower case letters
   !** End of interface *****************************************

#ifndef WITH_OLD_INPUT
   character(len=200) :: line
   character(len=200) :: line1
   integer(i4_kind) :: i,j

      input_which_namelist = .false.
      namelist_name = ""
!!$      if (eof(new_inp_unit)) return
      if (input_end_of_file()) return

      read(new_inp_unit,'(a200)') line
      backspace new_inp_unit

      call lowcase(line)
      if (check_string(line,'&')) then
         input_which_namelist = .true.
         i=0
         l1: do
            i=i+1
            if (line(i:i) /= '&') cycle l1
            j=0
            l2: do
               j=j+1
               if (line(i+j:i+j) == ' ') exit l1
               line1(j:j) = line(i+j:i+j)
            end do l2
         end do l1
         namelist_name=line1(1:j-1)
      end if
#else
      if ( .not. line_read ) then
         if ( store_associated ) then
            question_mode = .true.
            call input_readline()
            question_mode = .false.
         endif
         if ( .not. store_associated ) then
            input_which_namelist= .false.
            return
         endif
      endif
      if ( is_namelist ) then
         namelist_name = trim(namelistname)
      else
         namelist_name = ""
      endif
      input_which_namelist = is_namelist
#endif
   end function input_which_namelist
   !*************************************************************


#ifdef WITH_OLD_INPUT
   !*************************************************************
   subroutine input_readline()
   !  Purpose: reads in raw_line, does preprocessing and
   !           writes to processed_line.
   !           Errors and end of file are detected and handled
   !           by writting appropriate error message and
   !           terminating program
   !** End of interface *****************************************
   integer :: i_raw, indent_raw, i_namelistname
   logical :: escape, finished, add_blank, only_blanks, &
        is_namelistname, set_is_namelistname
   character :: c
!:VPP_BUG[
!  integer :: status
!  ! reset the input buffer processed_line
!  if (size(processed_line) /= block_length) then
!     deallocate(processed_line,stat=status)
!     if (status /= 0) call error_handler( &
!          "input_module(input_readline): deallocate (1) failed" )
!     allocate(processed_line(block_length),stat=status)
!     if (status /= 0) call error_handler( &
!          "input_module(input_readline): allocate (1) failed" )
!  endif
   ! don`t reset processed_line! (because of a compiler bug on VPP)
!:VPP_BUG]
   ! processed_line always contains one trailing blank!
   processed_line = " "
   i_processed = 1
   i_namelistname = 1
   finished = .false.
   escape = .false.
   namelist_open = .false.
   add_blank = .false.
   only_blanks = .true.
   is_namelist = .false.
   is_namelistname = .false.
   set_is_namelistname = .false.

   call read_rawline()
   do while ( .not. finished )
      call get_char()
      if ( escape ) then
         call put_char()
         escape = .false.
      else
         select case(c)
         case (achar(92)) ! escape
            escape = .true.
         case ('#') ! comment
            if (only_blanks .or. namelist_open) then
               call read_rawline()
            else
               finished = .true.
            endif
         case ('&') ! namelist begin
            if (only_blanks) then
               namelist_open = .true.
               is_namelist = .true.
               only_blanks = .false.
               set_is_namelistname = .true.
            endif
            call put_char()
            if (set_is_namelistname) then
               set_is_namelistname = .false.
               is_namelistname = .true.
               namelistname = repeat(" ",namelistname_max_length)
            endif
         case ('/') ! namelist end
            if (namelist_open) then
               namelist_open = .false.
               finished = .true.
            endif
            call put_char()
         case ('%') ! continue line
            call read_rawline()
         case (' ')
            is_namelistname = .false.
            call put_char()
         case default
            only_blanks = .false.
            call put_char()
         end select
      endif
   enddo

   line_read = store_associated

   contains

     subroutine get_line(iostat)
       integer, intent(out) :: iostat
       integer :: i
       !  Purpose: read of one record from the active raw input unit
       character(len=1), pointer :: raw_line_read_in(:)
       integer                   :: alloc_stat, chars_found
       !
       ! reset the input buffer raw_line
       if (size(raw_line) /= block_length) then
          deallocate(raw_line,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): deallocate (1) failed" )
          allocate(raw_line(block_length),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): allocate (1) failed" )
       endif
       length_raw = 0
       indent_raw = 1
       iostat = 1
       do while (.true.)
          read(int(store%unit), &
               format,advance='NO',size=chars_found,eor=770,end=772,err=773) &
               (raw_line(i),i=length_raw+1,length_raw+block_length)
          ! end of physical input line (EOR) not yet reached
          length_raw = length_raw + block_length
          ! enlarge the input buffer raw_line
          raw_line_read_in => raw_line
          allocate(raw_line(length_raw+block_length),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): allocate (2) failed" )
          raw_line(1:length_raw) = raw_line_read_in
          deallocate(raw_line_read_in,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): deallocate (2) failed" )
       end do
       ! end of physical input line (EOR) reached
       770 length_raw = length_raw + chars_found
       iostat = 0
       do indent_raw=1,length_raw
         if (raw_line(indent_raw) /= " ") return
       enddo
       indent_raw = length_raw + 1
       return
       !
       772 iostat = -1
       773 return
     end subroutine get_line

     subroutine read_rawline()
       !  Purpose: reads in raw input line, removing empty lines
       !** End of interface *****************************************
       integer :: unit, status
       !
       i_raw = 0
       do while ( i_raw == 0 .and. store_associated )
          store%line_nbr = store%line_nbr + 1
          call get_line(iostat=status)
          if (status > 0) call input_error("read error")
          if (status < 0) then
             call input_close_one()
             finished = .not. store_associated
          elseif (indent_raw <= length_raw) then
             if (raw_line(indent_raw) == "~") then
                unit = input_open_old &
                     (trim (adjustl (string (raw_line(indent_raw+1:length_raw)))))
             else
                if (raw_line(indent_raw) /= "#" .and. &
                    raw_line(indent_raw) /= "%" ) i_raw = indent_raw
             endif
          endif
       enddo
       return
     end subroutine read_rawline

     subroutine get_char()
       if (add_blank) then
          add_blank = .false.
          c = " "
          return
       endif
       c = raw_line(i_raw)
       if ( i_raw .eq. length_raw) then
          if (namelist_open) then
             if ( c /= "/" .and. c /= "%" ) then
                if ( c /= "#" ) call read_rawline()
                add_blank = .true.
             endif
          else
             finished = ( c /= "#" .or. .not. only_blanks ) .and. c /= "%"
          endif
       else
          i_raw = i_raw + 1
       endif
     end subroutine get_char

     subroutine put_char()
       integer :: status
       integer, parameter :: to_lower = ichar('a') - ichar('A')
       character(len=1), pointer :: already_processed_line(:)
       if (i_processed == size(processed_line)) then
          ! enlarge the input buffer processed_line
          already_processed_line => processed_line
          allocate(processed_line(i_processed+block_length),stat=status)
          if (status /= 0 ) call error_handler( &
               "input_module(put_char): allocate (1) failed" )
          processed_line(:i_processed) = already_processed_line
          processed_line(i_processed+1:) = " "
          deallocate(already_processed_line,stat=status)
          if (status /= 0 ) call error_handler( &
               "input_module(put_char): deallocate (1) failed" )
       endif
       processed_line(i_processed) = c
       if (is_namelistname) then
          if (i_namelistname .gt. namelistname_max_length) &
               call input_error("namelist name too long")
          if ('A' <= c .and. c <= 'Z') then
             namelistname(i_namelistname:i_namelistname) = char(ichar(c)+to_lower)
          else
             namelistname(i_namelistname:i_namelistname) = c
          endif
          i_namelistname = i_namelistname + 1
       endif
       i_processed = i_processed + 1
     end subroutine put_char

   end subroutine input_readline
   !*************************************************************
#endif


   !*************************************************************
   subroutine input_error(message)
   !  Purpose: Prints error message including message and actual
   !           input line number and line.
   !           Terminates calculation.
   !------------ Declaration of formal parameters ---------------
   character(len=*), intent(in), optional :: message
   !** End of interface *****************************************
   character(len=10)           :: nbr
   !------------ Executable code --------------------------------
   if ( store_associated ) then
      write(nbr,'(i10)') store%line_nbr
      if ( present(message) ) then
         call error_handler( "input_error: " // trim(message) // &
              " at line " // trim(nbr) // " of " // trim(store%filename) // &
              " : " // string(raw_line(:length_raw)) )
      else
         call error_handler( "input_error: at line " //  &
              trim(nbr) // " of " // trim(store%filename) //  &
              " : " // string(raw_line(:length_raw)) )
      endif
   else
      if ( present(message) ) then
         call error_handler( "input_error: " // trim(message) )
      else
         call error_handler( "input_error" )
      endif
   endif
   end subroutine input_error
   !*************************************************************

   subroutine input_open(filename)
   implicit none
   character(len=*), intent(in) :: filename
   ! *** end of interface ***

   integer :: old_inp_unit

#ifdef WITH_OLD_INPUT
   old_inp_unit = input_open_old (filename)
#else
   old_inp_unit = input_open_old (filename)

   ! open a file $TMP_DIR/new_input:
   new_inp_unit = input_open_new("new_input")

   ! dump the preprocessed input content into $TMP_DIR/new_input:
   call preprocess(old_inp_unit, new_inp_unit)

   call input_close_old()
#endif
   end subroutine input_open

   !*************************************************************
   function input_open_old (filename) result (unit)
   !
   ! Purpose: opens  new file and stores  unit if other  input file is
   ! open its state is stored
   !
   implicit none
   character (len=*), intent (in) :: filename
   integer :: unit
   !** End of interface *****************************************

   type (store_type), pointer :: store_intermediate
   integer :: status
   if ( store_associated ) then
      store_intermediate => store
      allocate( store, stat = status)
      if ( status .ne. 0 ) call error_handler( &
           "input_open: allocate (1) failed" )
      store%last => store_intermediate
   else
      allocate( store, stat = status)
      if ( status .ne. 0 ) call error_handler( &
           "input_open: allocate (2) failed" )
      store_associated = .true.
      first_store => store
      intermediate_unit = 0
      intermediate_open = .false.
      line_read = .false.
      question_mode = .false.
      namelist_open = .false.
      ! init the input buffers raw_line and processed_line
      allocate(raw_line(block_length),stat=status)
      if ( status .ne. 0 ) call error_handler( &
           "input_open: allocate (3) failed" )
      length_raw = 0
      ! processed_line always contains one trailing blank!
         allocate(processed_line(block_length),stat=status)
         if (status /= 0 ) call error_handler( &
              "input_readline: allocate (1) failed" )
         processed_line = " "
         i_processed = 1
   endif

   !
   ! Return the unit number. Status="old" was not there originally, so
   ! that  including  a  non-existent  file  by  ~file  directive  was
   ! silently creating a new empty file.
   !
   unit = openget_iounit (file=filename, action="read", status="old")

   store%unit = unit
   store%filename = filename
   store%line_nbr = 0
   end function input_open_old
   !*************************************************************


   !*************************************************************
   subroutine input_close_one()
   !  Purpose: closes file and restores state of outer file
   !           if there is no outer file, the program is
   !           terminated with error: unexpected end of file
   !** End of interface *****************************************
   type (store_type), pointer :: store_intermediate
   integer :: status
   if ( store_associated ) then
      call returnclose_iounit(store%unit)
      if ( .not. associated(store,first_store) ) then
         store_intermediate => store%last
         deallocate( store, stat = status)
         if ( status .ne. 0 ) call error_handler( &
              "input_close_one: deallocate (1) failed" )
         store => store_intermediate
      else
#ifndef WITH_OLD_INPUT
         question_mode = .true. !!!!!!!!!!!!!!!!!!AS
#endif
         if ( question_mode .and. .not. namelist_open ) then
            deallocate( store, stat = status)
            if ( status .ne. 0 ) call error_handler( &
                 "input_close_one: deallocate (2) failed" )
            store_associated = .false.
            ! close the input buffers processed_line and raw_line as well
#ifdef WITH_OLD_INPUT
               deallocate(processed_line,stat=status)
               if ( status .ne. 0 ) call error_handler( &
                    "input_close: deallocate (3) failed" )
#endif
            deallocate(raw_line,stat=status)
            if ( status .ne. 0 ) call error_handler( &
                 "input_close: deallocate (4) failed" )
         else
            call input_error("unexpected end of file")
         endif
      endif
   else
      call error_handler("unexpected end of file")
   endif
   end subroutine input_close_one
   !*************************************************************

   subroutine input_close()
   !  Purpose: closes all input file and if necessary intermediate file
   !** End of interface *****************************************

#ifdef WITH_OLD_INPUT
   call input_close_old()
#else
   call input_close_new()
#endif
   end subroutine input_close

   !*************************************************************
   subroutine input_close_old()
   !  Purpose: closes all input file and if necessary intermediate file
   !** End of interface *****************************************
   type (store_type), pointer :: store_intermediate
   integer :: status
   do while ( store_associated )
      call returnclose_iounit(store%unit)
      if ( .not. associated(store,first_store) ) then
         store_intermediate => store%last
         deallocate( store, stat = status)
         if ( status .ne. 0 ) call error_handler( &
              "input_close: deallocate (1) failed" )
         store => store_intermediate
      else
         deallocate( store, stat = status)
         if ( status .ne. 0 ) call error_handler( &
              "input_close: deallocate (2) failed" )
         store_associated = .false.
         ! close the input buffers processed_line and raw_line as well
            deallocate(processed_line,stat=status)
            if ( status .ne. 0 ) call error_handler( &
                 "input_close: deallocate (3) failed" )
         deallocate(raw_line,stat=status)
         if ( status .ne. 0 ) call error_handler( &
              "input_close: deallocate (4) failed" )
      endif
   enddo

   if ( intermediate_open ) &
        call returnclose_iounit(intermediate_unit, status="delete")
   intermediate_open = .false.
   intermediate_unit = 0
   end subroutine input_close_old
   !*************************************************************


#ifndef WITH_OLD_INPUT
   !*************************************************************
   subroutine preprocess(old, new)
    !
    ! Does the magic with input files
    !
    implicit none
    integer, intent(in) :: old, new ! unit numbers
    ! *** end of interface ***

#if 0
    character(len=1024) :: line

    !
    ! Pass input (almost) as is, trimming the trailing
    ! spaces and cutting lines to at most 1024 characters.
    ! Beware long numerical blocks that may extend beyond
    ! column 1024. FIXME: for a better "cat" implementation,
    ! use F2003 stream IO.
    !
    do
        read(old, '(A)', END=100) line
        write(new, '(A)') trim(line)
    enddo
100 CONTINUE

    !
    ! Seek to the begining of the file (legacy):
    !
    rewind new
#else
     integer :: i_raw, indent_raw, i_namelistname
     logical :: escape, finished, add_blank, only_blanks, &
          is_namelistname, set_is_namelistname
     character :: c

     do while (store_associated)
        processed_line = " "
        i_processed = 1
        i_namelistname = 1
        finished = .false.
        escape = .false.
        namelist_open = .false.
        add_blank = .false.
        only_blanks = .true.
        is_namelist = .false.
        is_namelistname = .false.
        set_is_namelistname = .false.

        call read_rawline()
        do while ( .not. finished )
           call get_char()
           if ( escape ) then
              call put_char()
              escape = .false.
           else
              select case(c)
              case (achar(92)) ! escape
                 escape = .true.
              case ('#') ! comment
                 if (only_blanks .or. namelist_open) then
                    call read_rawline()
                 else
                    finished = .true.
                 endif
              case ('&') ! namelist begin
                 if (only_blanks) then
                    namelist_open = .true.
                    is_namelist = .true.
                    only_blanks = .false.
                    set_is_namelistname = .true.
                 endif
                 call put_char()
                 if (set_is_namelistname) then
                    set_is_namelistname = .false.
                    is_namelistname = .true.
                    namelistname = repeat(" ",namelistname_max_length)
                 endif
              case ('/') ! namelist end
                 if (namelist_open) then
                    namelist_open = .false.
                    finished = .true.
                 endif
                 call put_char()
              case ('%') ! continue line
                 call read_rawline()
              case (' ')
                 is_namelistname = .false.
                 call put_char()
              case default
                 only_blanks = .false.
                 call put_char()
              end select
           endif
        enddo

        line_read = store_associated
        if(store_associated) call write_new_inp(new_inp_unit, processed_line)
        !write(new_inp_unit,*) processed_line

     end do
     rewind new_inp_unit

   contains

     subroutine write_new_inp(unit,chars)
       ! trims and writes
       integer(i4_kind), intent(in) :: unit
       character(len=1), intent(in) :: chars(:)
       ! *** end of interface ***

       integer(i4_kind) :: i,len

       do i=size(chars),1,-1
          if(chars(i)/=' ')then
             len = i
             exit ! loop
          endif
       enddo
       write(unit,*) chars(:len)
     end subroutine write_new_inp

     subroutine get_line(iostat)
       integer, intent(out) :: iostat
       integer :: i
       !  Purpose: read of one record from the active raw input unit
       character(len=1), pointer :: raw_line_read_in(:)
       integer                   :: alloc_stat, chars_found
       !
       ! reset the input buffer raw_line
       if (size(raw_line) /= block_length) then
          deallocate(raw_line,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): deallocate (1) failed" )
          allocate(raw_line(block_length),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): allocate (1) failed" )
       endif
       length_raw = 0
       indent_raw = 1
       iostat = 1
       do while (.true.)
          read(int(store%unit), &
               format,advance='NO',size=chars_found,eor=770,end=772,err=773) &
               (raw_line(i),i=length_raw+1,length_raw+block_length)
          ! end of physical input line (EOR) not yet reached
          length_raw = length_raw + block_length
          ! enlarge the input buffer raw_line
          raw_line_read_in => raw_line
          allocate(raw_line(length_raw+block_length),stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): allocate (2) failed" )
          raw_line(1:length_raw) = raw_line_read_in
          deallocate(raw_line_read_in,stat=alloc_stat)
          if (alloc_stat /= 0) call error_handler( &
               "input_module(get_line): deallocate (2) failed" )
       end do
       ! end of physical input line (EOR) reached
770    length_raw = length_raw + chars_found
       iostat = 0
       do indent_raw=1,length_raw
          if (raw_line(indent_raw) /= " ") return
       enddo
       indent_raw = length_raw + 1
       return
       !
772    iostat = -1
773    return
     end subroutine get_line

     subroutine read_rawline()
       !  Purpose: reads in raw input line, removing empty lines
       !** End of interface *****************************************
       integer :: unit, status
       !
       i_raw = 0
       do while ( i_raw == 0 .and. store_associated )
          store%line_nbr = store%line_nbr + 1
          call get_line(iostat=status)
          if (status > 0) call input_error("read error")
          if (status < 0) then
             call input_close_one()
             finished = .not. store_associated
          elseif (indent_raw <= length_raw) then
             if (raw_line(indent_raw) == "~") then
                unit = input_open_old &
                     (trim ( adjustl (string (raw_line(indent_raw+1:length_raw)) )))
             else
                if (raw_line(indent_raw) /= "#" .and. &
                     raw_line(indent_raw) /= "%" ) i_raw = indent_raw
             endif
          endif
       enddo
       return
     end subroutine read_rawline

     subroutine get_char()
       if (add_blank) then
          add_blank = .false.
          c = " "
          return
       endif
       c = raw_line(i_raw)
       if ( i_raw .eq. length_raw) then
          if (namelist_open) then
             if ( c /= "/" .and. c /= "%" ) then
                if ( c /= "#" ) call read_rawline()
                add_blank = .true.
             endif
          else
             finished = ( c /= "#" .or. .not. only_blanks ) .and. c /= "%"
          endif
       else
          i_raw = i_raw + 1
       endif
     end subroutine get_char

     subroutine put_char()
       integer :: status
       integer, parameter :: to_lower = ichar('a') - ichar('A')
       character(len=1), pointer :: already_processed_line(:)
       if (i_processed == size(processed_line)) then
          ! enlarge the input buffer processed_line
          already_processed_line => processed_line
          allocate(processed_line(i_processed+block_length),stat=status)
          if (status /= 0 ) call error_handler( &
               "input_module(put_char): allocate (1) failed" )
          processed_line(:i_processed) = already_processed_line
          processed_line(i_processed+1:) = " "
          deallocate(already_processed_line,stat=status)
          if (status /= 0 ) call error_handler( &
               "input_module(put_char): deallocate (1) failed" )
       endif
       processed_line(i_processed) = c
       if (is_namelistname) then
          if (i_namelistname .gt. namelistname_max_length) &
               call input_error("namelist name too long")
          if ('A' <= c .and. c <= 'Z') then
             namelistname(i_namelistname:i_namelistname) = char(ichar(c)+to_lower)
          else
             namelistname(i_namelistname:i_namelistname) = c
          endif
          i_namelistname = i_namelistname + 1
       endif
       i_processed = i_processed + 1
     end subroutine put_char
#endif
   end subroutine preprocess
   !*************************************************************


   !*************************************************************
   function input_open_new(filename) result(unit)
     use filename_module, only: tmpfile
     implicit none
     !------------ Declaration of formal parameters -------------
     character(len=*), intent(in) :: filename
     integer                      :: unit
     !** End of interface ***************************************

     DPRINT 'input_open_new: NEW_INPUT=', tmpfile('new_input')

     unit = openget_iounit(file=tmpfile(filename))
   end function input_open_new
   !*************************************************************


   !*************************************************************
   subroutine input_close_new()
     !------------ Declaration of formal parameters -------------
     !** End of interface ***************************************

     call returnclose_iounit(new_inp_unit,'delete')
   end subroutine input_close_new
   !*************************************************************
#endif


   !******************************************************************
   subroutine lowcase(string)

     character(len=*), intent(inout) :: string
     integer(kind=i4_kind) :: i,ln,ich

     ln = len(trim(string))
     do i=1,ln
        ich=iachar(string(i:i))
        if(ich>=65 .and. ich<=90) then
           ich=ich+32
           string(i:i)=achar(ich)
        end if
     end do

   end subroutine lowcase
   !******************************************************************


   !******************************************************************
   function check_string(string,word)

     character(len=*), intent(in) :: string
     character(len=*), intent(in) :: word

     logical :: check_string

     if(index(string,word) /= 0) then
        check_string = .true.
     else
        check_string = .false.
     end if

   end function check_string
   !******************************************************************

   function truefalse(word,stat) result(yes)
     implicit none
     character(len=*), intent(in)  :: word
     integer(i4_kind), intent(out) :: stat
     optional stat
     logical                       :: yes ! result
     ! *** end of interface ***

     integer(i4_kind) :: istat

     istat = 0
     select case( word )
     case ('FALSE','false','False','F','f','no')
       yes = .false.
     case ('TRUE','true','True','T','t','yes')
       yes = .true.
     case default
       ! default to false on failure:
       yes = .false.
       istat = 1
     end select
     if( present(stat) )then
       stat = istat
     else if( istat /= 0 )then
       ! the default behaviour on error:
       print*,'ERROR: Invalid value of logical flag=',word
       print*,'       Valid values are:'
       print*,'       FALSE, false, False, F, f, no'
       print*,'       TRUE, true, True, T, t, yes'
       ABORT('no such case, see tty')
     endif
   end function truefalse

!--------------- End of module ----------------------------------
end module input_module

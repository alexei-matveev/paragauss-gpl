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
!================================================================
! Public interface of module
!================================================================
module echo_input_module
  !---------------------------------------------------------------
  !
  !  Purpose: This module provides a couple of utilities for
  !           echoing the input into output files
  !
  !  Author: UB
  !  Date: 8/97
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   5/01
  ! Description: real_d and intg_d subroutines has been written
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module, only: r8_kind, i4_kind
  implicit none
  private
  save

  !== Interrupt end of public interface of module =================

  !------------ public parameters ---------------------------------
  integer(i4_kind), public, parameter :: echo_level_short   = 0, &
                                         echo_level_default = 1, &
                                         echo_level_full    = 2

  !------------ public variables ----------------------------------

  integer, parameter :: MAX_FMT = 35
  character(len=MAX_FMT), public, protected :: flag_format  = '("    ",a," = ",a20:" # ",a)'
  character(len=MAX_FMT), public, protected :: intg_format  = '("    ",a," = ",i20:" # ",a)'
  character(len=MAX_FMT), public, protected :: real_format  = '("    ",a," = ",g20.14:" # ",a)'

  character(len=MAX_FMT), public :: &
       real_format1 = '("    ",a," = ",f10.8:" # ",a)', &
                            & real_format2 = '("    ",a2," = ",f10.8:" # ",a2)', &
                            & real_format3 = '("    ",a7," = ",es10.3:" # ",a7)', &
                            & real_format4 = '("    ",a ," = ",E20.10:" # ",a)', &
                            & word_format  = '("    ",a," = ",a10  :" # ",a)', &
                            & string_format= '("    ",a," = ",a    :" # ",a)'

  !------------ public functions and subroutines ------------------
  public :: start, flag, intg, real, word, stop, echo
  public :: strng, intg_d, real_d

!================================================================
! End of public interface of module
!================================================================

  !------------ private variables ----------------------------------
  character(len=78), private :: list_name, prog_name
  integer(i4_kind) , private :: wrte_unit, echo_mode
  logical          , private :: written

  contains

  subroutine start(name,prog,unit,level)
    ! Start writing a namelist
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in) :: name, prog
    integer(i4_kind), intent(in) :: unit
    integer(i4_kind), intent(in) :: level
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer  status
    external error_handler

    ! first check input strings
    if (len_trim(name) > len(list_name)) call error_handler &
         (trim(prog)//": namelist name '"//trim(name)//"' too long")
    if (len_trim(prog) > len(prog_name)) call error_handler &
         (trim(prog)//": subroutine name '"//trim(prog)//"' too long")
    if (truncated(flag_format)) call error_handler &
         (trim(prog)//": flag_format '"//trim(flag_format)//"' too long")
    if (truncated(intg_format)) call error_handler &
         (trim(prog)//": intg_format '"//trim(intg_format)//"' too long")
    if (truncated(real_format)) call error_handler &
         (trim(prog)//": real_format '"//trim(real_format)//"' too long")
    if (truncated(word_format)) call error_handler &
         (trim(prog)//": word_format '"//trim(word_format)//"' too long")
    if (truncated(real_format1)) call error_handler &
         (trim(prog)//": real_format1 '"//trim(real_format1)//"' too long")
    if (truncated(real_format2)) call error_handler &
         (trim(prog)//": real_format2 '"//trim(real_format2)//"' too long")
    if (truncated(real_format3)) call error_handler &
         (trim(prog)//": real_format3 '"//trim(real_format3)//"' too long")
    list_name = trim(name)
    prog_name = trim(prog)
    wrte_unit = unit
    echo_mode = level
    written   = .false.

    write(wrte_unit,'(" &",a )',iostat=status)trim(list_name)

    if (status /= 0) call error_handler( trim(prog_name)// &
         ": write namelist header '"//trim(list_name)//"' failed")

    contains

    logical function truncated(format)
       character(len=*) :: format
       character(len=1) :: char
       integer(i4_kind) :: pos, level

       level = 0
       do pos=1,len(format)
          char = format(pos:pos)
          if (char == "(") level = level + 1
          if (char == ")") level = level - 1
       enddo

       truncated = level /= 0
    end function truncated

  end subroutine start

  !***************************************************************

  subroutine stop(keep_namelist,empty_line,drop_namelist)
    ! Finish writing a namelist
    logical, optional :: keep_namelist ! default = echo_level controlled
    logical, optional :: empty_line    ! default = .true.
    logical, optional :: drop_namelist ! default = echo_level controlled
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    character(len=9) :: format, comment*33
    integer     :: status, pos
    logical     :: keep, drop
    external error_handler

    if (present(keep_namelist)) then
       keep = keep_namelist
    else
       keep = .false.
    endif
    if (present(drop_namelist)) then
       drop = drop_namelist
    else
       drop = .false.
    endif
    format = '(" /",a/)'
    if (present(empty_line)) then
       if (.not.empty_line) format = '(" /",a)'
    endif

    if (written .or. ( keep .and. echo_mode == echo_level_short ) ) then
       pos = index(trim(list_name),' ')
       if (pos == 0) then
          write(wrte_unit,format,iostat=status)trim(list_name)
       else
          write(wrte_unit,format,iostat=status)trim(list_name(:pos-1))
       endif
       if (status /= 0) call error_handler( trim(prog_name)// &
            ": write namelist '"//trim(list_name)//"' terminator failed")
    elseif (echo_mode == echo_level_default) then ! .not.written
       if (keep) then
          comment = " (empty namelist must remain)"
       else
          if (drop) then
             comment = " (empty namelist must be dropped)"
          else
             comment = " (namelist may be dropped)"
          endif
       endif
       pos = index(trim(list_name),' ')
       if (pos == 0) then
          write(wrte_unit,format,iostat=status)trim(list_name)// &
                trim(comment)
       else
          write(wrte_unit,format,iostat=status)trim(list_name(:pos-1))// &
                trim(comment)
       endif
       if (status /= 0) call error_handler( trim(prog_name)// &
            ": write namelist '"//trim(list_name)//"' terminator failed")
    else
       backspace(wrte_unit,iostat=status)
       if (status /= 0) call error_handler( trim(prog_name)// &
            ": skipping namelist header '"//trim(list_name)//"' failed")
    endif

  end subroutine stop

  !***************************************************************

  logical function echo(keep_namelist)
    ! if called after subrouitine stop it tells whether the last
    ! namelist has been echoed to the output or not.
    logical, optional :: keep_namelist ! default = echo_level controlled
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    logical     :: keep

    if (present(keep_namelist)) then
       keep = keep_namelist
    else
       keep = .false.
    endif
    echo = written .or. ( keep .and. echo_mode == echo_level_short ) .or. &
           echo_mode == echo_level_default
  end function echo

  !***************************************************************

  subroutine flag(name,value,default,comment)
     ! Write a logical input parameter
     !------------ Declaration of formal parameters ---------------
     character(len=*),           intent(in) :: name
     logical         ,           intent(in) :: value, default
     character(len=*), optional, intent(in) :: comment
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     integer(i4_kind)                       :: status
     character(len=32)                      :: commentstr

     commentstr = ' '
     if ( present(comment) ) commentstr = comment

     status = 0
     if (echo_mode == echo_level_full .or. (value .neqv. default)) then
        if (value) then
           write( wrte_unit, flag_format, iostat=status )                      &
           name," TRUE",trim(commentstr)
        else
           write( wrte_unit, flag_format, iostat=status )                      &
           name,"FALSE",trim(commentstr)
        endif
        written = .true.
     elseif (echo_mode == echo_level_default) then ! value .eqv. default
        if (value) then
           write(wrte_unit,flag_format,iostat=status)name," TRUE", &
                trim(commentstr)//"(the default)"
        else
           write(wrte_unit,flag_format,iostat=status)name,"FALSE", &
                trim(commentstr)//"(the default)"
        endif
     endif
     if (status /= 0) call error_handler( trim(prog_name)// &
          ": write logical parameter '"//trim(name)//"' failed")

  end subroutine flag

  !***************************************************************

  subroutine intg(name, value, default)
     ! Write an integer input parameter
     !------------ Declaration of formal parameters ---------------
     character(len=*), intent(in) :: name
     integer(i4_kind), intent(in) :: value, default
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     integer(i4_kind)             :: status

     status = 0
     if (echo_mode == echo_level_full .or. (value /= default)) then
        write(wrte_unit, intg_format, iostat=status) name, value
        written = .true.
     elseif (echo_mode == echo_level_default) then ! value == default
        write(wrte_unit, intg_format, iostat=status) name, value, &
              "(the default)"
     endif
     if (status /= 0) call error_handler( trim(prog_name)// &
          ": write integer parameter '"//trim(name)//"' failed")

  end subroutine intg

  !***************************************************************

  subroutine intg_d(name,ivalues,idefaults,length)
     ! AS 02/2001
     ! Write a integer input parameters
     !------------ Declaration of formal parameters ---------------
     integer(kind=i4_kind), intent(in) :: length
     character(len=*),      intent(in) :: name
     integer(kind=i4_kind), intent(in) :: ivalues(length), idefaults(length)
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     character(len=256) :: i_format,i_format1
!!$     character(len=38) :: i_format,i_format1
     character(len=2) :: buf
     integer(i4_kind) :: status,length1,i

     length1=length
     if(length == 0) length1=1
     status = 0
     write(buf,'(i2)') length1
!!$     i_format='("    ",a," = ",'//trim(buf)//'(1x,i5)" # ",a)'
!!$     i_format1='("    ",a," = ",'//trim(buf)//'(1x,i5))'
     i_format='("    ",a," = ",'
     i_format1='("    ",a," = ",'
     do i=1,length1
        i_format=trim(i_format)//'i5,","'
        i_format1=trim(i_format1)//'i5,","'
     end do
     i_format=trim(i_format)//'," # ",a)'
     i_format1=trim(i_format1)//')'
     if(echo_mode == echo_level_full .or. .not.check_values_i()) then
        write(wrte_unit,trim(i_format1),iostat=status)name,ivalues(1:length1)
     elseif(echo_mode == echo_level_default) then ! values == defaults
        write(wrte_unit,trim(i_format),iostat=status)name,ivalues(1:length1), "(the default)"
     endif
     if (status /= 0) call error_handler( trim(prog_name)// &
          ": write integer parameterS '"//trim(name)//"' failed")

   contains
     logical function check_values_i()
       !------------ Declaration of local variables---------------
       integer(kind=i4_kind) :: i

       check_values_i=.true.
       do i=1,length
          if(ivalues(i) /= idefaults(i)) check_values_i=.false.
       enddo
     end function check_values_i
   end subroutine intg_d

   !******************************************************************

  subroutine real(name,value,default,format,fmt,comment)
     ! Write a real input parameter
     !------------ Declaration of formal parameters ---------------
     character(len=*), intent(in)           :: name
     real(r8_kind)   , intent(in)           :: value, default
     integer(i4_kind), intent(in), optional :: format
     character(len=*), intent(in), optional :: fmt, comment
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     integer(i4_kind)                       :: status

     character(len=MAX_FMT) :: outfmt
     character(len=MAX_FMT) :: commentstr

     commentstr = ' '
     if ( present(comment) ) commentstr = comment

     if(present(format).and.present(fmt))then
        ABORT('two format specs')
     endif

     outfmt = real_format
     if( present(fmt) )then
        outfmt = fmt
     elseif (present(format)) then
        select case(format)
        case (1)
           outfmt = real_format1
        case (2)
           outfmt = real_format2
        case (3)
           outfmt = real_format3
        case (4)
           outfmt = real_format4
        case default
           ABORT('no such format')
        end select
     endif

     status = 0
     if (echo_mode == echo_level_full .or. (value /= default)) then
        write(wrte_unit,outfmt,iostat=status) name, value, trim(commentstr)
        written = .true.
     elseif (echo_mode == echo_level_default) then ! value == default
        write(wrte_unit,outfmt,iostat=status) name, value                      &
                                            , trim(commentstr)//"(the default)"
     endif
     if (status /= 0) call error_handler( trim(prog_name)// &
          ": write real parameter '"//trim(name)//"' failed")

  end subroutine real

  !***************************************************************

  subroutine real_d(name,values,defaults,length)
     ! AS 02/2001
     ! Write a real input parameters
     !------------ Declaration of formal parameters ---------------
     integer(kind=i4_kind), intent(in) :: length
     character(len=*),      intent(in) :: name
     real(kind=r8_kind),    intent(in) :: values(length), defaults(length)
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     character(len=256) :: r_format,r_format1
!!$     character(len=38) :: r_format,r_format1
     character(len=2) :: buf
     integer(i4_kind) :: status,length1,i

     length1=length
     if(length == 0) length1=1
     status = 0
     write(buf,'(i2)') length1
!!$     r_format='("    ",a," = ",'//trim(buf)//'(1x,f15.10)" # ",a)'
!!$     r_format1='("    ",a," = ",'//trim(buf)//'(1x,f15.10))'
     r_format='("    ",a," = ",'
     r_format1='("    ",a," = ",'
     do i=1,length1
        r_format=trim(r_format)//'f15.10,","'
        r_format1=trim(r_format1)//'f15.10,","'
     end do
     r_format=trim(r_format)//'," # ",a)'
     r_format1=trim(r_format1)//')'
     if(echo_mode == echo_level_full .or. .not.check_values_r()) then
        write(wrte_unit,trim(r_format1),iostat=status)name,values(1:length1)
     elseif(echo_mode == echo_level_default) then ! values == defaults
        write(wrte_unit,trim(r_format),iostat=status)name,values(1:length1), "(the default)"
     endif
     if (status /= 0) call error_handler( trim(prog_name)// &
          ": write real parameterS '"//trim(name)//"' failed")

   contains
     logical function check_values_r()
       !------------ Declaration of local variables---------------
       integer(kind=i4_kind) :: i

       check_values_r=.true.
       do i=1,length
          if(values(i) /= defaults(i)) check_values_r=.false.
       enddo
     end function check_values_r
   end subroutine real_d
  !***************************************************************

  subroutine word(name,value,default)
     ! Write an integer input parameter
     !** End of interface *****************************************
     !------------ Declaration of formal parameters ---------------
     character(len=*), intent(in) :: name
     character(len=*), intent(in) :: value, default
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     integer(i4_kind)             :: status

     status = 0
     if (echo_mode == echo_level_full .or. (value /= default)) then
        write(wrte_unit,word_format,iostat=status)name,'"'//value//'"'
        written = .true.
     elseif (echo_mode == echo_level_default) then ! value == default
        write(wrte_unit,word_format,iostat=status)name,'"'//value//'"', &
              "(the default)"
     endif
     if (status /= 0) call error_handler( trim(prog_name)// &
          ": write character parameter '"//trim(name)//"' failed")

  end subroutine word

  subroutine strng(name,value,default)
    ! Write an integer input parameter
    !** End of interface *****************************************
    !------------ Declaration of formal parameters ---------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value, default
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(i4_kind)             :: status

    status = 0
    if (echo_mode == echo_level_full .or. (value /= default)) then
       write(wrte_unit,string_format,iostat=status)name,'"'//trim(value)//'"'
       written = .true.
       elseif (echo_mode == echo_level_default) then ! value == default
       write(wrte_unit,string_format,iostat=status)name,'"'//trim(value)//'"', &
            "(the default)"
    endif
    if (status /= 0) call error_handler( trim(prog_name)// &
         ": write string parameter '"//trim(name)//"' failed")

  end subroutine strng

  !***************************************************************

end module echo_input_module

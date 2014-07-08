!
! ParaGauss,  a program package  for high-performance  computations of
! molecular systems
!
! Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
! F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
! A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
! D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
! S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
! A. Nikodem, T. Soini, M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify
! it under  the terms of the  GNU General Public License  version 2 as
! published by the Free Software Foundation [1].
!
! This program is distributed in the  hope that it will be useful, but
! WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
! MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
  !===================================================================
! Public interface of module
  !===================================================================
module readwriteblocked_module
!
! Purpose: high performance buffered IO of
!          real(kind=r8_kind) arrays.
!          There is a machine dependent ideal blocklength
!          for read and write operations. This module
!          offers routines both for reading and writing.
!          When writing, the arrays are stored in a buffer
!          with length blocklength till the buffer is full
!          and then the buffer is written in one sweep.
!          When reading, one buffer is read at one sweep and
!          the arrays are then fetched from the buffer.
!
!  To allow several read and write operations to be performed
!  at once, a handle technique is introduced:
!  The routines of this module operate on handles associated
!  with files passed as arguments to them.
!
!  Three modes of operation exist:
!  1. Standard mode:
!       the last buffer is written / read completely
!       even if it is not full. Fine for large files.
!       Do not use optional parameters variable_length
!       and total_length of subroutines readwriteblocked_startwrite,
!       readwriteblocked_stopwrite and readwriteblocked_startread.
!       (appending of data not possible in standard mode)
!  2. Smart mode:
!       All buffers are written together with their sizes.
!       Only the filled part of the last buffer is written.
!       The price is lower reading performance. For smaller files.
!       Set variable_length = .true. when invoking
!       readwriteblocked_startwrite and readwriteblocked_startread.
!       (appending of data possible)
!  3. Optimized mode:
!       The total size of the file is calculated when writing
!       and has to be spezified when reading. The programmer
!       has to save this total size in between.
!       Only the filled part of the last buffer is written.
!       Optimal performance for files of any size.
!       Inquire size by optional parameter total_length to
!       readwriteblocked_stopwrite and pass this parameter to
!       readwriteblocked_startread.
!       (appending of data possible if total_length is provided)
!
! Module called by: predominantly everything which needs
!          3-Z-Integrals, such as ham_calc_module, int_send_2cob3c_module
!
! Author: Folke Noertemann, translated from f77 to f90
! Date  : 8/95
!
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: TB
! Date:   2.10.95
! Description: Changed to widget technique and introduced
!              type readwriteblocked_tapehandle.
!              changed name of module from readtape_module
!              and writetape_module and combined them.
!              Changed names of subroutines to include
!              module name.
!              Moved parameter definition of lenblk into
!              the module.
!
! Modification (Please copy before editing)
! Author: TB
! Date:   2/97
! Description: Introduced smart and optimized mode and adjustable
!              blocklength. The start / stop routines now also
!              open / close the files and allocate / free the
!              internal buffer associated with th.
!
!---------------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!-------------------------------------------------

!#define FPP_TIMERS 2
#include "def.h"
use type_module ! type specification parameters
use iounitadmin_module
#ifndef FPP_OPTIMIZER
use machineparameters_module, only : lenblk
#endif
implicit none
private         ! by default, all names are private
save

!--- declaration of private variables ------------
#ifdef FPP_OPTIMIZER
integer(kind=i4_kind), parameter ,private  :: lenblk= 8*1024  !-> machineparameters_module
      ! This is the maximum block size that the
      ! NAG-f90 Compiler  allows to be read in or
      ! written out.
#endif

!== Interrupt end of public interface of module =================



!------------ Declaration of types ------------------------------
type readwriteblocked_tapehandle
   private
   !== Interrupt of public interface of module =====================
   integer(kind=i4_kind)  :: left,i_start,io_unit,lenblk,lentot,lenfin
   integer(kind=i4_kind)  :: rec
   real(kind=r8_kind), pointer :: buffer(:)
   character(len=150) :: filename
   logical :: variable_length
   !== Interrupt end of public interface of module =================
end type readwriteblocked_tapehandle
public readwriteblocked_tapehandle


!------------ public functions and subroutines ------------------
public readwriteblocked_read, &
       readwriteblocked_startread, &
       readwriteblocked_stopread, &
       readwriteblocked_skipread, &
       readwriteblocked_write, &
       readwriteblocked_startwrite, &
       readwriteblocked_stopwrite, &
       readwriteblocked_blocklength, &
       readwriteblocked_returnclose &
      ,readwriteblocked_rec,readwriteblocked_setrec &
      ,readwriteblocked_filename


  !===================================================================
  ! End of public interface of module
  !===================================================================



FPP_TIMER_DECL(rea)
FPP_TIMER_DECL(wrt)

!---------------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine readwriteblocked_read(item,th,eof,rec)
    !  Purpose: reads in item(length) from th
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    real(kind=r8_kind),                intent(out)   :: item(:)
    type(readwriteblocked_tapehandle), intent(inout) :: th
    logical, optional,                 intent(out)   :: eof
    ! returns .TRUE. if end-of-file is reached during reading
    ! not allowed in standard mode
    integer(kind=i4_kind), intent(inout), optional :: rec
    !** End of interface *****************************************
    ! --- declaration of local variables  ----------------------
    integer(kind=i4_kind) :: ns,nn,length
    !------------ Executable code ------------------------------
    FPP_TIMER_START(rea)
    if (present(eof)) then
       eof = .false.
       if (.not.th%variable_length .and. th%lenfin == -1) &
            call error_handler("readwriteblocked_read: &
            &argument 'eof' not allowed in standard mode")
    endif
    ns = 0
    nn = int(ubound(item,1),kind=i4_kind)
!       print*, th%io_unit,th%left, th%i_start,'th%io_unit th%left th%i_start'
    do
       if (th%left .ge. nn) then
          item(ns+1:ns+nn) = th%buffer(th%i_start+1:th%i_start+nn)
          th%i_start = th%i_start + nn
          th%left = th%left - nn
          th%lentot = th%lentot + nn
    FPP_TIMER_STOP(rea)
          return
       else
          item(ns+1:ns+th%left) = th%buffer(th%i_start+1:th%i_start+th%left)
          ns = ns + th%left
          nn = nn - th%left
          th%lentot = th%lentot + th%left
          sm_mod: if ( th%variable_length ) then
             ! smart mode
             read(int(th%io_unit),ERR=99,END=98) length, th%buffer(1:length)
             th%left = length
          elseif ( th%lentot .gt. th%lenfin - th%lenblk .and. &
               th%lenfin .ne. -1 ) then
             ! optimized mode, last block
             length = th%lenfin - th%lentot
             if (length <= 0) goto 98
           if(present(rec)) then
!           print*,'read_sub opt mod'//trim(th%filename),rec
              call read_sub(th%buffer,length,rec=rec)
           else
              call read_sub(th%buffer,length)
           endif
             th%left = length
          else sm_mod
             ! standard mode and optimized mode before last block
            if(present(rec)) then
!            print*,'standart read_sub'//trim(th%filename),rec
             call read_sub(th%buffer,th%lenblk,rec=rec)
!            print*,th%buffer(1),th%buffer(th%lenblk)
            else
             call read_sub(th%buffer,th%lenblk)
            endif
             th%left = th%lenblk
          endif sm_mod
          th%i_start = 0
       endif
    enddo

    FPP_TIMER_STOP(rea)
    return

98  if (present(eof)) then
       eof = .true.
       th%left = 0
    FPP_TIMER_STOP(rea)
       return
    endif
    call error_handler("readwriteblocked_read: end of file "//trim(th%filename))
99  call error_handler("readwriteblocked_read: read error for file "//trim(th%filename))

  contains

    subroutine read_sub(array,length,rec)
      integer(kind=i4_kind), intent(in) :: length
      real(kind=r8_kind), intent(out) :: array(length)
      integer(kind=i4_kind), intent(inout), optional :: rec
      if(present(rec)) then
        read(int(th%io_unit),ERR=999,rec=rec) array
        rec=rec+1
        else
        read(int(th%io_unit),ERR=999,END=998) array
        endif
      return
998   call error_handler("readwriteblocked_read: read_sub: end of file "//trim(th%filename))
999   call error_handler("readwriteblocked_read: read_sub: read error for file "//trim(th%filename))
    end subroutine read_sub

  end subroutine readwriteblocked_read
  !*************************************************************


  !*************************************************************
  subroutine readwriteblocked_startread(filename,th,blocklength, &
       variable_length,total_length,rec,unit)
    !  Purpose: initialises th
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    character(len=*),                  intent(in)    :: filename
    type(readwriteblocked_tapehandle), intent(inout) :: th
    integer(kind=i4_kind), optional,   intent(in)    :: blocklength
    ! default can be obtained by readwriteblocked_blocklength()
    logical,               optional,   intent(in)    :: variable_length
    ! default: false
    ! Set true for smart mode:
    ! the records are read together with their length, allowing
    ! the last record to be shorter. Must be set to the same value
    ! at readwriteblocked_startwrite
    integer(kind=i4_kind), optional,   intent(in)    :: total_length
    ! total number of floats contained in  file. In case this argument
    ! If given optimized mode is choosen:
    ! the last record is not read totally and the same
    ! argument has to be suplied at readwriteblocked_stopwrite.
    integer(kind=i4_kind), optional,   intent(inout) :: rec
    integer(kind=i4_kind), optional,   intent(inout) :: unit
    !** End of interface *****************************************
    integer :: status
    !------------ Executable code ------------------------------
    FPP_TIMER_START(rea)
    th%filename=filename
    th%i_start = 0
    th%left = 0
    th%lentot = 0
    if ( present(total_length) ) then
       th%lenfin = total_length
    else
       th%lenfin = -1
    endif
    if ( present(blocklength) ) then
       th%lenblk = blocklength
    else
       th%lenblk = lenblk
    endif
    if ( present(variable_length) ) then
       th%variable_length = variable_length
    else
       th%variable_length = .false.
    endif
    if(present(rec)) then
     if(rec.eq.0) then
      th%io_unit = openget_iounit(filename, status='OLD', &
         form='UNFORMATTED', action='READ', access='DIRECT', &
         recl=th%lenblk * bit_size(r8_kind) )
      if(present(unit)) unit=th%io_unit
      rec=1
      else
      if(present(unit)) th%io_unit=unit
     endif
    else
     th%io_unit = openget_iounit(filename, status='OLD', &
         form='UNFORMATTED', action='READ', position='REWIND', &
         recl=th%lenblk * bit_size(r8_kind) )
     if ( th%variable_length ) th%lenblk = th%lenblk - 1
    endif
    allocate(th%buffer(th%lenblk),stat=status)
    if ( status .ne. 0 ) call error_handler( &
         "readwriteblocked_startread: allocate failed")
    FPP_TIMER_STOP(rea)
  end subroutine readwriteblocked_startread
  !*************************************************************


  !*************************************************************
  subroutine readwriteblocked_stopread(th,status,rec)
    !  Purpose: resets th and closes the corresponding io-unit
    !           ONLY if the file was opened at all.
    !           This was done to prevent error messages in the case
    !           that this routine is called also if e.g. the array
    !           of fitcoeff`s for this host is empty.
    use iounitadmin_module, only: return_iounit
    implicit none
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    character(len=*), optional, intent(in) :: status ! for close
    integer(kind=i4_kind),intent(in),optional:: rec
    !** End of interface *****************************************
    ! --- decalaration of local variables ---------------------
    integer :: stat
    !------------ Executable code ------------------------------
    FPP_TIMER_START(rea)
    if (th%io_unit.ne.0) then
       if(.not.present(rec))  &
                 call returnclose_iounit(int(th%io_unit),status)
       deallocate(th%buffer,stat=stat)
       if ( stat .ne. 0 ) call error_handler( &
            "readwriteblocked_stopread: deallocate failed")
    endif
    FPP_TIMER_STOP(rea)
    DPRINT  'RWB:  READ time: ',FPP_TIMER_VALUE(rea)
  end subroutine readwriteblocked_stopread
  !*************************************************************


  !*************************************************************
  subroutine readwriteblocked_skipread(length,th,eof)
    !  Purpose: skips length real(kind=r8_kind) numbers from th
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    integer(kind=i4_kind),             intent(in)    :: length
    logical, optional,                 intent(out)   :: eof
    ! returns .TRUE. if end-of-file is reached during skipping
    ! not allowed in standard mode
    !** End of interface *****************************************
    ! --- declaration of local variables  ----------------------
    integer(kind=i4_kind) nn, ll
    !------------ Executable code ------------------------------
    FPP_TIMER_START(rea)
    if (present(eof)) then
       eof = .false.
       if (.not.th%variable_length .and. th%lenfin == -1) &
            call error_handler("readwriteblocked_skipread: &
            &argument 'eof' not allowed in standard mode")
    endif

    ! in smart mode each record may have a different length, thus:
    nn = length
    do
       if (nn.le.th%left) then
          th%i_start = th%i_start + nn
          th%left = th%left - nn
          th%lentot = th%lentot + nn
    FPP_TIMER_STOP(rea)
          return
       else
          nn = nn - th%left
          th%lentot = th%lentot + th%left
          if ( th%variable_length ) then
             ! smart mode
             if (nn > th%lenblk) then
                read(int(th%io_unit),ERR=99,END=98) ll
             else
                read(int(th%io_unit),ERR=99,END=98) ll, th%buffer(1:ll)
             endif
             th%left = ll
          elseif( th%lentot .gt. th%lenfin - th%lenblk .and. &
               th%lenfin .ne. -1 ) then
             ! optimized mode last block
             ll = th%lenfin - th%lentot
             if (ll <= 0) goto 98
             call skip_sub(th%buffer,ll)
             th%left = ll
          else
             ! standard mode and optimized mode before last block
             call skip_sub(th%buffer,th%lenblk)
             th%left = th%lenblk
          endif
          th%i_start = 0
       endif
    enddo

    FPP_TIMER_STOP(rea)
    return

98  if (present(eof)) then
       eof = .true.
       th%left = 0
    FPP_TIMER_STOP(rea)
       return
    endif
    call error_handler("readwriteblocked_skipread: end of file "//trim(th%filename))
99  call error_handler("readwriteblocked_skipread: read error for " &
         // trim(th%filename) )

  contains

    subroutine skip_sub(array,length)
      integer(kind=i4_kind), intent(in) :: length
      real(kind=r8_kind), intent(out) :: array(length)
      if (nn >= length) then
         read(int(th%io_unit),ERR=999,END=998)
      else
         read(int(th%io_unit),ERR=999,END=998) array
      endif
      return
998   call error_handler("readwriteblocked_skipread: skip_sub: end of file "//trim(th%filename))
999   call error_handler("readwriteblocked_skipread: skip_sub: read error for file "//trim(th%filename))
   end subroutine skip_sub

  end subroutine readwriteblocked_skipread
  !*************************************************************


  !*************************************************************
  subroutine readwriteblocked_write(item,th,rec)
    !  Purpose: writes item(:) to th
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    real(kind=r8_kind),                intent(in)    :: item(:)
    type(readwriteblocked_tapehandle), intent(inout) :: th
    integer(kind=i4_kind), intent(inout), optional   :: rec
    !** End of interface *****************************************
    ! --- declaration of local variables  ----------------------
    integer(kind=i4_kind) :: ns,nn,io_stat
    logical               :: virtual_open
    character(len=255)    :: err_msg,iostat_str
    !------------ Executable code ------------------------------
    FPP_TIMER_START(wrt)
    ns = 0
    nn = int(ubound(item,1),kind=i4_kind)
    th%lentot = th%lentot + nn
    virtual_open = th%io_unit == -1
    do
       if (th%left.gt.nn) then
          th%buffer(th%i_start+1:th%i_start+nn) = item(ns+1:ns+nn)
          th%i_start = th%i_start + nn
          th%left = th%left - nn
          exit
       else
          th%buffer(th%i_start+1:th%i_start+th%left) = item(ns+1:ns+th%left)
          ns = ns + th%left
          nn = nn - th%left
          th%i_start = 0
          th%left = th%lenblk
          if (th%io_unit == -1)call readwriteblocked_openwrite(th)
          if ( th%variable_length ) then
             write(int(th%io_unit),iostat=io_stat) th%lenblk,th%buffer
             if(io_stat/=0) then
                write(iostat_str,'(i5)') io_stat
                err_msg="readwriteblocked_write: write error 1 (iostat="
                err_msg=trim(err_msg)//trim(iostat_str)//") for"
                err_msg=trim(err_msg)//" "//trim(th%filename)
                call error_handler(trim(err_msg))
             end if
          else
           if(present(rec)) then
!          print*,'write sub'//trim(th%filename),th%buffer(1),th%buffer(size(th%buffer)),size(th%buffer),rec
             write(int(th%io_unit),iostat=io_stat,rec=rec) th%buffer
             rec=rec+1
           else
             write(int(th%io_unit),iostat=io_stat) th%buffer
           endif
             if(io_stat/=0) then
                write(iostat_str,'(i5)') io_stat
                err_msg="readwriteblocked_write: write error 2 (iostat="
                err_msg=trim(err_msg)//trim(iostat_str)//") for"
                err_msg=trim(err_msg)//" "//trim(th%filename)
                call error_handler(trim(err_msg))
             end if
          endif
          if (nn.eq.0) exit
       endif
    enddo
    if (virtual_open .and. th%io_unit /= -1)call readwriteblocked_closewrite(th)
    FPP_TIMER_STOP(wrt)
    return
  end subroutine readwriteblocked_write
  !*************************************************************

  !*************************************************************
  subroutine readwriteblocked_closewrite(th)
    !  Purpose: closes file without destroying the information in
    !           the tapehandle and the item buffer. This is usefull
    !           if you want to close a file temporarily and re-open
    !           it later on to append more data
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************
    ! --- declaration of local variables  ----------------------
    !------------ Executable code ------------------------------
    FPP_TIMER_START(wrt)
    call returnclose_iounit(th%io_unit,status='KEEP')
    th%io_unit = -1
    FPP_TIMER_STOP(wrt)
  end subroutine readwriteblocked_closewrite
  !*************************************************************

  !*************************************************************
  subroutine readwriteblocked_openwrite(th)
    !  Purpose: open file which has been cloes with
    !           readwriteblocked_closewrite before
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !** End of interface *****************************************
    ! --- declaration of local variables  ----------------------
    !------------ Executable code ------------------------------
    FPP_TIMER_START(wrt)
    th%io_unit=openget_iounit(th%filename,status='OLD',&
        form='UNFORMATTED', action='WRITE', position='APPEND', &
         recl=th%lenblk*bit_size(r8_kind))
    FPP_TIMER_STOP(wrt)
  end subroutine readwriteblocked_openwrite
  !*************************************************************

  !*************************************************************
  subroutine readwriteblocked_startwrite(filename,th,blocklength,&
             variable_length,append,total_length,virtual_open,rec,unit)
    !  Purpose: initialises th
    !------------ Modules used ---------------------------------
    use type_module
    implicit none
    ! --- declaration of formal parameters ---------------------
    character(len=*),                  intent(in)    :: filename
    type(readwriteblocked_tapehandle), intent(inout) :: th
    integer(kind=i4_kind), optional,   intent(in)    :: blocklength
    integer(kind=i4_kind), optional,   intent(out)   :: rec
    integer(kind=i4_kind), optional,   intent(inout) :: unit
    ! default can be obtained by readwriteblocked_blocklength()
    logical,               optional,   intent(in)    :: variable_length
    ! default: false
    ! Set true for smart mode:
    ! the records are written together with their length, allowing
    ! the last record to be shorter. Must be set to the same value
    ! at readwriteblocked_startread
    logical,               optional,   intent(in)    :: append
    ! default: false
    ! Set true for appending (not allowed in standard mode)
    integer(kind=i4_kind), optional,   intent(in)    :: total_length
    ! Only required for appending in optimized mode
    logical,               optional,   intent(in)    :: virtual_open
    ! default: false
    ! Set true if the file should be re-opened for each buffer
    ! flushing and closed again immediately afterwords. This
    ! mode is far less effective than the defauklt mode and
    ! should only be used if IO units are getting rare.
    ! The virtual_open option can be changed after starting by
    ! readwriteblocked_closewrite => virtual_open = .true. and
    ! readwriteblocked_openwrite => virtual_open = .false.
    !** End of interface *****************************************
    integer :: status
    logical :: recreate
    !------------ Executable code ------------------------------
    FPP_TIMER_START(wrt)
    th%filename=filename
    th%i_start = 0
    th%lentot = 0
    if ( present(total_length) ) then
       th%lenfin = total_length
    else
       th%lenfin = -1
    endif
    if ( present(blocklength) ) then
       th%lenblk = blocklength
    else
       th%lenblk = lenblk
    endif
    th%left = th%lenblk
    if ( present(variable_length) ) then
       th%variable_length = variable_length
    else
       th%variable_length = .false.
    endif
    if (present(append)) then
       if (.not.th%variable_length .and. th%lenfin == -1) &
            call error_handler("readwriteblocked_startwrite: &
            &argument 'append' not allowed standard mode")
       recreate = .not.append
    else
       recreate = .true.
    endif
    if (recreate) then
       if(present(rec)) then
        if(rec.lt.1) then
        th%io_unit = openget_iounit(filename, status='REPLACE', &
            form='UNFORMATTED', action='WRITE', &
            recl=th%lenblk * bit_size(r8_kind),access='DIRECT' )
            rec=1
            th%rec=1
            if(present(unit))  unit=th%io_unit
        else
         if(present(unit)) th%io_unit = unit
        endif
       else
       th%io_unit = openget_iounit(filename, status='REPLACE', &
            form='UNFORMATTED', action='WRITE', position='REWIND', &
            recl=th%lenblk * bit_size(r8_kind) )
       endif
    else
       th%io_unit = openget_iounit(filename, form='UNFORMATTED', &
            action='READWRITE', status='OLD', position='APPEND',&
            recl=th%lenblk * bit_size(r8_kind) )
    endif
    if (th%variable_length) then
       th%lenblk = th%lenblk - 1
       th%left = th%left - 1
    endif

    if(present(rec)) then
!    if(rec.eq.0) then
    allocate(th%buffer(th%lenblk),stat=status)
    if ( status .ne. 0 ) call error_handler( &
         "readwriteblocked_startwrite: allocate failed")
!     endif
     else
    allocate(th%buffer(th%lenblk),stat=status)
    if ( status .ne. 0 ) call error_handler( &
         "readwriteblocked_startwrite: allocate failed")
    endif

    if (.not.recreate .and. .not.th%variable_length) then
       ! read last (partial) record of optimized mode file
       th%i_start = mod(th%lenfin,th%lenblk)
       if (th%i_start > 0) then
          th%left = th%left - th%i_start
          backspace(int(th%io_unit))
          read(int(th%io_unit),ERR=99,END=98) th%buffer(1:th%i_start)
          backspace(int(th%io_unit))
       endif
       th%lentot = th%lenfin
    endif
    if (present(virtual_open)) then
       if (virtual_open) then
          endfile(th%io_unit) ! to ensure overwriting of old data
          call readwriteblocked_closewrite(th)
       endif
    endif

    FPP_TIMER_STOP(wrt)
    return

98  call error_handler("readwriteblocked_startwrite: &
         &end of file "//trim(th%filename))
99  call error_handler("readwriteblocked_startwrite: &
         &read error for file "//trim(th%filename))

  end subroutine readwriteblocked_startwrite
  !*************************************************************


  !*************************************************************
  subroutine readwriteblocked_stopwrite(th,total_length,rec)
    !  Purpose: closes file associate with th
    implicit none
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    integer(kind=i4_kind), intent(out), optional :: total_length
    integer(kind=i4_kind), intent(inout), optional :: rec
    ! total number of floats written to file. In case this argument
    ! is given optimized mode is invoked:
    ! The last record is not written totally and the same
    ! argument has to be suplied at readwriteblocked_startread.
    !** End of interface *****************************************
    ! --- decalaration of local variables ---------------------
    integer :: stat
    !------------ Executable code ------------------------------
    FPP_TIMER_START(wrt)
    if ( present(total_length) ) total_length = th%lentot
    if (th%io_unit > 0 .or. th%io_unit == -1) then
       if (th%left.ne.th%lenblk) then
          if (th%io_unit == -1)call readwriteblocked_openwrite(th)
          if ( th%variable_length ) then
             ! smart mode
             write(int(th%io_unit),ERR=99) th%i_start, th%buffer(1:th%i_start)
          elseif ( present(total_length) ) then
             ! optimized mode
            if(present(rec)) then
!            print*,'i_start ',th%i_start,size(th%buffer),associated(th%buffer), &
!                              th%buffer(1),th%buffer(th%i_start),int(th%io_unit)
             write(int(th%io_unit),rec=rec) th%buffer(1:th%i_start)
             rec=rec+1
            else
             write(int(th%io_unit),ERR=99) th%buffer(1:th%i_start)
            endif
          else
             ! standard mode
            if(present(rec)) then
             write(int(th%io_unit),rec=rec) th%buffer
             rec=rec+1
            else
             write(int(th%io_unit)) th%buffer
            endif
          endif
       endif
       if(present(rec)) then
!        print*,'skiped to next record buffer deallocated'
       else
        if (th%io_unit /= -1)call returnclose_iounit(th%io_unit,status='KEEP')
       endif

!       if(.not.present(rec)) then
         deallocate(th%buffer,stat=stat)
          if ( stat .ne. 0 ) call error_handler( &
            "readwriteblocked_stopwrite: deallocate failed")
!       endif
    endif
    FPP_TIMER_STOP(wrt)
    DPRINT  'RWB: WRITE time: ',FPP_TIMER_VALUE(wrt)
    return
99  call error_handler("readwriteblocked_stopwrite: write error for " &
         // trim(th%filename) )
  end subroutine readwriteblocked_stopwrite
  !*************************************************************

  subroutine readwriteblocked_returnclose(th)
   implicit none
   type(readwriteblocked_tapehandle), intent(in) :: th
   call returnclose_iounit(th%io_unit,status='KEEP')
  end subroutine readwriteblocked_returnclose

  !*************************************************************
  integer function  readwriteblocked_give_iounit(th)
    !  Purpose: returns iounit of th
    implicit none
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(in) :: th
    !** End of interface *****************************************
    !------------ Executable code ------------------------------
    readwriteblocked_give_iounit = th%io_unit
  end function readwriteblocked_give_iounit
  !*************************************************************

  function  readwriteblocked_filename(th)
    !  Purpose: returns iounit of th
    implicit none
    character(len=150) :: readwriteblocked_filename
    ! --- declaration of formal parameters ---------------------
    type(readwriteblocked_tapehandle), intent(in) :: th
    !** End of interface *****************************************
    !------------ Executable code ------------------------------
    readwriteblocked_filename = th%filename
  end function readwriteblocked_filename


  !*************************************************************
  integer(kind=i4_kind) function  readwriteblocked_blocklength()
    !  Purpose: returns standard blocklength
    !** End of interface *****************************************
    readwriteblocked_blocklength = lenblk
  end function readwriteblocked_blocklength
  !*************************************************************

  integer(kind=i4_kind) function  readwriteblocked_rec(th)
    !  Purpose: returns rec number
   implicit none
   ! --- declaration of formal parameters ---------------------
   type(readwriteblocked_tapehandle), intent(in) :: th
    !** End of interface *****************************************
    readwriteblocked_rec = th%rec
  end function readwriteblocked_rec

 subroutine readwriteblocked_setrec(th,rec)
 implicit none
 type(readwriteblocked_tapehandle), intent(inout) :: th
 integer(kind=i4_kind),intent(in) :: rec
 th%rec=rec
 end subroutine readwriteblocked_setrec

end module readwriteblocked_module

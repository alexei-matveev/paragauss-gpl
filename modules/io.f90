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
!=====================================================================
! Public interface of module
!=====================================================================
module io
  !-------------------------------------------------------------------
  !
  !  Purpose: Yet another IO interface...
  !           Fortran IO is very fragile, and one often needs
  !           to workaround "big" record sizes or compiler bugs.
  !           It would be good to workaround those at one place.
  !
  !  Usage:   call file_open(iou,file,'w') ! open for 'w'rite, 'r'ead or 'a'ppend
  !           call file_write(iou,...)
  !           call file_close(iou,'k')     ! 'k'eep or 'd'elete
  !
  !  Note:    One should not assume that ``iou'' is the *Fortran* unit number,
  !           though normally it is ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
# include "def.h"
  use type_module, only:&
       & IK => i4_kind,&
       & RK => r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  integer(IK), parameter, public ::&
       & IO_OK            =  0, &
       & IO_EOR           = -1, &
       & IO_EOF           = -2, &
       & IO_LINE_TOO_LONG = -100, &
         IO_ABORT         =  1  , &
         IO_WARN          =  2  , &
         IO_SILENT        =  3


  !------------ Interface statements ---------------------------------

  interface file_open
     module procedure file_open          ! (iou,file-name,rwa-mode)
  end interface

  interface file_close
     module procedure file_close         ! (iou,kd-mode)
  end interface

  interface file_write
     module procedure file_write_int           ! (iou,int)
     module procedure file_write_int_1D        ! (iou,buf(:))
     module procedure file_write_int_2D        ! (iou,buf(:,:))
     module procedure file_write_int_3D        ! (iou,buf(:,:,:))
     module procedure file_write_real          ! (iou,real)
     module procedure file_write_real_1D       ! (iou,buf(:))
     module procedure file_write_real_2D       ! (iou,buf(:,:))
     module procedure file_write_real_3D       ! (iou,buf(:,:,:))
     ! shorthands for open/write/close:
     module procedure file_write_file_real_1D  ! (file-name,buf(:))
     module procedure file_write_file_real_2D  ! (file-name,buf(:,:))
     module procedure file_write_file_real_3D  ! (file-name,buf(:,:,:))
  end interface

  interface file_read
     module procedure file_read_int           ! (iou,int)
     module procedure file_read_int_1D        ! (iou,buf(:))
     module procedure file_read_int_2D        ! (iou,buf(:,:))
     module procedure file_read_int_3D        ! (iou,buf(:,:,:))
     module procedure file_read_real          ! (iou,real)
     module procedure file_read_real_1D       ! (iou,buf(:))
     module procedure file_read_real_2D       ! (iou,buf(:,:))
     module procedure file_read_real_3D       ! (iou,buf(:,:,:))
     ! shorthands for open/read/close:
     module procedure file_read_file_real_1D  ! (file-name,buf(:))
     module procedure file_read_file_real_2D  ! (file-name,buf(:,:))
     module procedure file_read_file_real_3D  ! (file-name,buf(:,:,:))
  end interface

  interface readline
     module procedure readline_line
     module procedure readline_real
     module procedure readline_int
     module procedure readline_char
  end interface

  interface write_buffer ! MUST DIE!
!    module procedure file_write_real_1D       ! (io,buf(:))
!    module procedure file_write_real_2D       ! (io,buf(:,:))
!    module procedure file_write_real_3D       ! (io,buf(:,:,:))
     module procedure file_write_file_real_1D  ! (file-name,buf(:))
     module procedure file_write_file_real_2D  ! (file-name,buf(:,:))
     module procedure file_write_file_real_3D  ! (file-name,buf(:,:,:))
  end interface

  interface read_buffer ! MUST DIE!
!    module procedure file_read_real_1D       ! (io,buf(:))
!    module procedure file_read_real_2D       ! (io,buf(:,:))
!    module procedure file_read_real_3D       ! (io,buf(:,:,:))
     module procedure file_read_file_real_1D  ! (file-name,buf(:))
     module procedure file_read_file_real_2D  ! (file-name,buf(:,:))
     module procedure file_read_file_real_3D  ! (file-name,buf(:,:,:))
  end interface

  !------------ public functions and subroutines ---------------------

  public :: file_open
  public :: file_close
  public :: file_write
  public :: file_read

  public :: write_buffer ! MUST DIE!
  public :: read_buffer  ! MUST DIE!

  public :: readline

  public :: io_set_error_handler
  public :: io_status

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  integer(IK), parameter :: line_length = 256
  character(line_length) :: line

  integer(IK)            :: on_error      = IO_ABORT
  integer(IK)            :: on_error_save = 0
  integer(IK)            :: iostat = 0



  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  function io_status() result(stat)
    implicit none
    integer(IK) :: stat
    ! *** end of interface ***

    stat = iostat
  end function io_status

  subroutine io_set_error_handler(action)
    implicit none
    integer(IK),intent(in),optional :: action
    ! *** end of interface ***

    if ( present(action) ) then
      on_error = action
    else
      on_error = IO_ABORT
    endif
  end subroutine io_set_error_handler

#define SETERR(stat)     call seterr(stat)
  subroutine seterr(stat)
    implicit none
    integer(IK), intent(out), optional :: stat
    ! *** end of interface ***

    if(present(stat))then
      ! errors are not fatal:
      on_error_save = on_error
      on_error = IO_SILENT
    else
      on_error_save = on_error
      on_error = IO_ABORT
    endif
  end subroutine seterr

#define GETERR(stat)     call geterr(stat)
  subroutine geterr(stat)
    implicit none
    integer(IK), intent(out), optional :: stat
    ! *** end of interface ***

    if(present(stat))then
      ! errors are not fatal:
      stat = iostat
    else
      ! catch errors earlier:
      ASSERT(iostat==0)
    endif
    ! restore default error handler:
    on_error = on_error_save
  end subroutine geterr

  subroutine readline_char(io,res,stat)
    implicit none
    integer(IK),intent(in)  :: io
    character,  intent(out) :: res
    integer(IK),intent(out) :: stat
    ! *** end of interface ***

    stat = IO_OK

    read(io,'(A)', advance='no',EOR=1001,END=1002) res
    return

1001   continue ! EOR exeption
    DPRINT 'io::readline_char: EOR!'
    stat = IO_EOR
    return

1002   continue ! EOF exeption
    DPRINT 'io::readline_char: EOF!'
    stat = IO_EOF
    return
  end subroutine readline_char

  subroutine readline_line(io,line,n,stat)
    implicit none
    integer(IK),intent(in)        :: io   ! unit
    character(len=*), intent(out) :: line ! storage
    integer(IK),intent(inout)     :: n ! in=want,out=got
    integer(IK),intent(out)       :: stat
    ! *** end of interface ***

    integer(IK) :: strlen,i,nn
    character   :: char

    strlen = LEN(line)
    ASSERT(n<=strlen)

    stat = IO_OK
    nn = 0
    do i=1,n
       call readline_char(io,char,stat)
       if(stat/=IO_OK)then
          exit ! loop
       endif
       line(i:i) = char
       nn = nn + 1
    enddo
    n = nn
  end subroutine readline_line

  subroutine readline_real(io,res,stat)
    use strings, only: isblank
    implicit none
    integer(IK),intent(in)  :: io
    real(RK),   intent(out) :: res
    integer(IK),intent(out) :: stat
    ! *** end of interface ***

    character   :: char
    integer(IK) :: n

    DPRINT 'io::readline_real: entered; unit=',io

    stat = IO_OK

    do while(stat.eq.IO_OK)
!!$       read(io,'(A)', advance='no',EOR=1001,END=1002) char
       call readline_char(io,char,stat)
       if(.not. isblank(char)) exit
    enddo
    n = 0
    scan: do while(stat.eq.IO_OK)
       n = n + 1
       if(n>line_length)then
          DPRINT 'io::readline_real: line too long!'
          stat = IO_LINE_TOO_LONG
          exit scan
       endif
       line(n:n) = char

!!$       read(io,'(A)', advance='no',EOR=1001,END=1002) char
       call readline_char(io,char,stat)
       if(isblank(char)) exit scan
       cycle scan

    enddo scan

    if(n>0) read(line(1:n),*) res
    DPRINT 'io::readline_real: ', res, ' stat=',stat
  end subroutine readline_real

  subroutine readline_int(io,res,stat)
    use strings, only: isblank
    implicit none
    integer(IK),intent(in)  :: io
    integer(IK),intent(out) :: res
    integer(IK),intent(out) :: stat
    ! *** end of interface ***

    character   :: char
    integer(IK) :: n

    DPRINT 'io::readline_int: entered; unit=',io

    stat = IO_OK

    do while(stat.eq.IO_OK)
!!$       read(io,'(A)', advance='no',EOR=1001,END=1002) char
       call readline_char(io,char,stat)
       if(.not.isblank(char)) exit
    enddo
    n = 0
    scan: do while(stat.eq.IO_OK)
       n = n + 1
       if(n>line_length)then
          DPRINT 'io::readline_int: line too long!'
          stat = IO_LINE_TOO_LONG
          exit scan
       endif
       line(n:n) = char

!!$       read(io,'(A)', advance='no',EOR=1001,END=1002) char
       call readline_char(io,char,stat)
       if(isblank(char)) exit scan
       cycle scan

    enddo scan

    if(n>0) read(line(1:n),*) res
    DPRINT 'io::readline_int: ', res, ' stat=',stat
  end subroutine readline_int

  subroutine file_write_real(io,r)
    implicit none
    integer(IK),intent(in) :: io
    real(RK),   intent(in) :: r
    ! *** end of interface ***

    call file_write_real_buf(io,(/r/),1)
  end subroutine file_write_real

  subroutine file_write_real_1D(io,buf)
    implicit none
    integer(IK),intent(in) :: io
    real(RK),   intent(in) :: buf(:)
    ! *** end of interface ***

    call file_write_real_buf(io,buf,size(buf))
  end subroutine file_write_real_1D

  subroutine file_read_real(io,r)
    implicit none
    integer(IK),intent(in)  :: io
    real(RK),   intent(out) :: r
    ! *** end of interface ***

    real(RK) :: buf(1)

    call file_read_real_buf(io,buf,1)
    r = buf(1)
  end subroutine file_read_real

  subroutine file_read_real_1D(io,buf)
    implicit none
    integer(IK),intent(in)  :: io
    real(RK),   intent(out) :: buf(:)
    ! *** end of interface ***

    call file_read_real_buf(io,buf,size(buf))
  end subroutine file_read_real_1D

  subroutine file_write_real_2D(io,buf)
    implicit none
    integer(IK),intent(in) :: io
    real(RK),   intent(in) :: buf(:,:)
    ! *** end of interface ***

    call file_write_real_buf(io,buf,size(buf))
  end subroutine file_write_real_2D

  subroutine file_read_real_2D(io,buf)
    implicit none
    integer(IK),intent(in)  :: io
    real(RK),   intent(out) :: buf(:,:)
    ! *** end of interface ***

    call file_read_real_buf(io,buf,size(buf))
  end subroutine file_read_real_2D

  subroutine file_write_real_3D(io,buf)
    implicit none
    integer(IK),intent(in) :: io
    real(RK),   intent(in) :: buf(:,:,:)
    ! *** end of interface ***

    call file_write_real_buf(io,buf,size(buf))
  end subroutine file_write_real_3D

  subroutine file_read_real_3D(io,buf)
    implicit none
    integer(IK),intent(in)  :: io
    real(RK),   intent(out) :: buf(:,:,:)
    ! *** end of interface ***

    call file_read_real_buf(io,buf,size(buf))
  end subroutine file_read_real_3D

  subroutine file_write_int(io,i)
    implicit none
    integer(IK),intent(in) :: io
    integer(IK),intent(in) :: i
    ! *** end of interface ***

    call file_write_int_buf(io,(/i/),1)
  end subroutine file_write_int

  subroutine file_write_int_1D(io,buf)
    implicit none
    integer(IK),intent(in) :: io
    integer(IK),intent(in) :: buf(:)
    ! *** end of interface ***

    call file_write_int_buf(io,buf,size(buf))
  end subroutine file_write_int_1D

  subroutine file_write_int_2D(io,buf)
    implicit none
    integer(IK),intent(in) :: io
    integer(IK),intent(in) :: buf(:,:)
    ! *** end of interface ***

    call file_write_int_buf(io,buf,size(buf))
  end subroutine file_write_int_2D

  subroutine file_write_int_3D(io,buf)
    implicit none
    integer(IK),intent(in) :: io
    integer(IK),intent(in) :: buf(:,:,:)
    ! *** end of interface ***

    call file_write_int_buf(io,buf,size(buf))
  end subroutine file_write_int_3D

  subroutine file_read_int(io,i,stat)
    implicit none
    integer(IK),intent(in)  :: io
    integer(IK),intent(out) :: i
    integer(IK),intent(out),optional :: stat
    ! *** end of interface ***

    integer(IK) :: buf(1)

    SETERR(stat)
    call file_read_int_buf(io,buf,1)
    i = buf(1)
    GETERR(stat)
  end subroutine file_read_int

  subroutine file_read_int_1D(io,buf,stat)
    implicit none
    integer(IK),intent(in)  :: io
    integer(IK),intent(out) :: buf(:)
    integer(IK),intent(out),optional :: stat
    ! *** end of interface ***

    SETERR(stat)
    call file_read_int_buf(io,buf,size(buf))
    GETERR(stat)
  end subroutine file_read_int_1D

  subroutine file_read_int_2D(io,buf,stat)
    implicit none
    integer(IK),intent(in)  :: io
    integer(IK),intent(out) :: buf(:,:)
    integer(IK),intent(out),optional :: stat
    ! *** end of interface ***

    SETERR(stat)
    call file_read_int_buf(io,buf,size(buf))
    GETERR(stat)
  end subroutine file_read_int_2D

  subroutine file_read_int_3D(io,buf,stat)
    implicit none
    integer(IK),intent(in)  :: io
    integer(IK),intent(out) :: buf(:,:,:)
    integer(IK),intent(out),optional :: stat
    ! *** end of interface ***

    SETERR(stat)
    call file_read_int_buf(io,buf,size(buf))
    GETERR(stat)
  end subroutine file_read_int_3D


  subroutine file_write_real_buf(io,buf,siz)
    ! DISK INTERFACE: working horse for writing
    implicit none
    integer(IK),intent(in) :: io, siz
    real(RK),   intent(in) :: buf(*)
#ifdef WITH_FXDR
    include 'fxdr.inc'
#endif
    ! *** end of interface ***

#ifdef WITH_FXDR
    ! FXDR call to read/write DOUBLE arrays:
    iostat = IXDRDMAT(io,siz,buf)
    ! NOTE: read or write mode is defined when opening...
#else
#ifndef _IO_SEQ_FAILSAFE
    ! variable RECL = buf(:siz)
    DPRINT 'file_write_real_buf(',siz,')'
    WRITE(io,IOSTAT=iostat) buf(:siz)
#else
    ! constant RECL = 1 REAL(r8_kind)
    ! may consume twice the space:
    integer(IK)            :: i
    do i=1,siz
       WRITE(io) buf(i)
    enddo
#endif
#endif
    if( iostat /= 0 )then
      select case ( on_error )
      case ( IO_ABORT )
        print *,'file_write_real_buf: IO error: got IOSTAT=',iostat,' while writing unit ',io
        ABORT('file_write_real_buf: IO error')
      case ( IO_WARN )
        print *,'file_write_real_buf: ignore IO error: got IOSTAT=',iostat,' while writing unit ',io
        WARN('file_write_real_buf: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif
  end subroutine file_write_real_buf

  subroutine file_read_real_buf(io,buf,siz)
    ! DISK INTERFACE: working horse for reading
    implicit none
    integer(IK),intent(in)  :: io, siz
    real(RK),   intent(out) :: buf(*)
#ifdef WITH_FXDR
    include 'fxdr.inc'
#endif
    ! *** end of interface ***

#ifdef WITH_FXDR
    ! FXDR call to read/write DOUBLE arrays:
    iostat = IXDRDMAT(io,siz,buf)
    ! NOTE: read or write mode is defined when opening...
#else
#ifndef _IO_SEQ_FAILSAFE
    ! variable RECL = buf(:siz)
    DPRINT 'file_read_real_buf(',siz,')'
    READ(io,IOSTAT=iostat) buf(:siz)
#else
    ! constant RECL = 1 REAL(r8_kind)
    ! may consume twice the space:
    integer(IK)             :: i
    do i=1,siz
       READ(io) buf(i)
    enddo
#endif
#endif
    if( iostat /= 0 )then
      select case ( on_error )
      case ( IO_ABORT )
        print *,'file_read_real_buf: IO error: got IOSTAT=',iostat,' while reading unit ',io
        ABORT('file_read_real_buf: IO error')
      case ( IO_WARN )
        print *,'file_read_real_buf: ignore IO error: got IOSTAT=',iostat,' while reading unit ',io
        WARN('file_read_real_buf: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif
  end subroutine file_read_real_buf

  subroutine file_write_int_buf(io,buf,siz)
    ! DISK INTERFACE: working horse for writing
    implicit none
    integer(IK),intent(in) :: io, siz
    integer(IK),intent(in) :: buf(*)
#ifdef WITH_FXDR
    include 'fxdr.inc'
#endif
    ! *** end of interface ***

#ifdef WITH_FXDR
    ! FXDR call to read/write INTEGER*4 arrays:
    iostat = IXDRIMAT(io,siz,buf)
    ! NOTE: read or write mode is defined when opening...
#else
#ifndef _IO_SEQ_FAILSAFE
    ! variable RECL = buf(:siz)
    DPRINT 'file_write_int_buf(',siz,')'
    WRITE(io,IOSTAT=iostat) buf(:siz)
#else
    ! constant RECL = 1
    ! may consume twice the space:
    integer(IK)            :: i
    do i=1,siz
       WRITE(io) buf(i)
    enddo
#endif
#endif
    if( iostat /= 0 )then
      select case ( on_error )
      case ( IO_ABORT )
        print *,'file_write_int_buf: IO error: got IOSTAT=',iostat,' while writing unit ',io
        ABORT('file_write_int_buf: IO error')
      case ( IO_WARN )
        print *,'file_write_int_buf: ignore IO error: got IOSTAT=',iostat,' while writing unit ',io
        WARN('file_write_int_buf: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif
  end subroutine file_write_int_buf

  subroutine file_read_int_buf(io,buf,siz)
    ! DISK INTERFACE: working horse for reading
    implicit none
    integer(IK),intent(in)  :: io, siz
    integer(IK),intent(out) :: buf(*)
#ifdef WITH_FXDR
    include 'fxdr.inc'
#endif
    ! *** end of interface ***

#ifdef WITH_FXDR
    ! FXDR call to read/write INTEGER*4 arrays:
    iostat = IXDRIMAT(io,siz,buf)
    ! NOTE: read or write mode is defined when opening...
#else
#ifndef _IO_SEQ_FAILSAFE
    ! variable RECL = buf(:siz)
    DPRINT 'file_read_int_buf(',siz,')'
    READ(io,IOSTAT=iostat) buf(:siz)
#else
    ! constant RECL = 1
    ! may consume twice the space:
    integer(IK)             :: i
    do i=1,siz
       READ(io) buf(i)
    enddo
#endif
#endif
    if( iostat /= 0 )then
      select case ( on_error )
      case ( IO_ABORT )
        print *,'file_read_int_buf: IO error: got IOSTAT=',iostat,' while reading unit ',io
        ABORT('file_read_int_buf: IO error')
      case ( IO_WARN )
        print *,'file_read_int_buf: ignore IO error: got IOSTAT=',iostat,' while reading unit ',io
        WARN('file_read_int_buf: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif
  end subroutine file_read_int_buf

  subroutine file_write_file_real_1D(fn,buf)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(in)  :: buf(:)
    ! *** end of interface ***

    call file_write_file_real_buf(fn,buf,size(buf))
  end subroutine file_write_file_real_1D

  subroutine file_read_file_real_1D(fn,buf)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(out) :: buf(:)
    ! *** end of interface ***

    call file_read_file_real_buf(fn,buf,size(buf))
  end subroutine file_read_file_real_1D

  subroutine file_write_file_real_2D(fn,buf)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(in)  :: buf(:,:)
    ! *** end of interface ***

    call file_write_file_real_buf(fn,buf,size(buf))
  end subroutine file_write_file_real_2D

  subroutine file_read_file_real_2D(fn,buf)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(out) :: buf(:,:)
    ! *** end of interface ***

    call file_read_file_real_buf(fn,buf,size(buf))
  end subroutine file_read_file_real_2D

  subroutine file_write_file_real_3D(fn,buf)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(in)  :: buf(:,:,:)
    ! *** end of interface ***

    call file_write_file_real_buf(fn,buf,size(buf))
  end subroutine file_write_file_real_3D

  subroutine file_read_file_real_3D(fn,buf)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(out) :: buf(:,:,:)
    ! *** end of interface ***

    call file_read_file_real_buf(fn,buf,size(buf))
  end subroutine file_read_file_real_3D

#ifdef WITH_FXDR
  subroutine file_open(iou,fn,mode)
    implicit none
    integer(IK)     , intent(out) :: iou
    character(len=*), intent(in)  :: fn
    character(len=1), intent(in)  :: mode ! 'r', 'w', or 'a'
    ! external subs and return codes:
    include 'fxdr.inc'
    ! *** end of interface ***

    DPRINT 'file_open_fxdr: ',trim(fn),' mode=',mode
    select case ( mode )
    case ( 'r','w','a' )
      ! accept
    case default
      ABORT('file_open_fxdr: action?')
    end select

    ! FXDR call:
    iou = INITXDR( fn, mode, .TRUE. )
    ! .TRUE. means return error codes, dont abort...

    if( iou <= 0 )then
      print *,'file_open_fxdr: IO error: got iou=',iou,' while opening ',trim(fn)
      select case ( on_error )
      case ( IO_ABORT )
        ABORT('file_open_fxdr: IO error')
      case ( IO_WARN )
        WARN('file_open_fxdr: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif
  end subroutine file_open

  subroutine file_close(iou,mode)
    use iounitadmin_module, only: return_iounit
    implicit none
    integer(IK)     , intent(inout) :: iou
    character(len=1), intent(in)    :: mode ! 'k', or 'd'
    ! external subs and return codes:
    include 'fxdr.inc'
    ! *** end of interface ***

    integer(IK) :: stat

    DPRINT 'file_close_fxdr: unit',iou,' mode=',mode
    select case ( mode )
    case ( 'k' )
      ! accept
    case ( 'd' )
      WARN('file_close_fxdr: dont know how to delete!')
    case default
      ABORT('file_close_fxdr: action?')
    end select

    stat = IXDRCLOSE(iou)

    if( stat /= 0 )then
      print *,'file_close_fxdr: IO error: got IOSTAT=',stat,' while closing unit',iou
      select case ( on_error )
      case ( IO_ABORT )
        ABORT('file_close_fxdr: IO error')
      case ( IO_WARN )
        WARN('file_close_fxdr: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif
    ! return zero or error code:
    iou = stat
  end subroutine file_close
#else
  subroutine file_open(iou,fn,mode)
    use iounitadmin_module, only: get_iounit, return_iounit
    implicit none
    integer(IK)     , intent(out) :: iou
    character(len=*), intent(in)  :: fn
    character(len=1), intent(in)  :: mode ! 'r', 'w', or 'a'
    ! *** end of interface ***

    character(len=12) :: action
    character(len=12) :: status
    character(len=12) :: position

    DPRINT 'file_open_fortran: ',trim(fn),' mode=',mode
    select case ( mode )
    case ( 'r' )
      action   = 'read'
      status   = 'old'
      position = 'rewind'
    case ( 'w' )
      action   = 'write'
      status   = 'unknown'
      position = 'rewind'
    case ( 'a' )
      action   = 'write'
      status   = 'old'
      position = 'append'
    case default
      ABORT('file_open_fortran: action?')
    end select

    iou = get_iounit()
    OPEN(  iou                      &
         , FILE     = trim(fn)      &
         , FORM     = 'unformatted' &
         , ACTION   = action        &
         , POSITION = position      &
         , STATUS   = status        &
         , IOSTAT   = iostat        &
         )

    if( iostat /= 0 )then
      print *,'file_open_fortran: IO error: got IOSTAT=',iostat,' while opening ',trim(fn)
      select case ( on_error )
      case ( IO_ABORT )
        ABORT('file_open_fortran: IO error')
      case ( IO_WARN )
        WARN('file_open_fortran: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
      ! return io unit:
      call return_iounit(iou)
      iou = -1
    endif
  end subroutine file_open

  subroutine file_close(iou,mode)
    use iounitadmin_module, only: return_iounit
    implicit none
    integer(IK)     , intent(inout) :: iou
    character(len=1), intent(in)    :: mode ! 'k', or 'd'
    ! *** end of interface ***

    character(len=12) :: status

    DPRINT 'file_close_fortran: unit',iou,' mode=',mode
    select case ( mode )
    case ( 'k' )
      status   = 'keep'
    case ( 'd' )
      status   = 'delete'
    case default
      ABORT('file_close_fortran: action?')
    end select

    CLOSE( iou             &
         , STATUS = status &
         , IOSTAT = iostat &
         )

    if( iostat /= 0 )then
      print *,'file_close_fortran: IO error: got IOSTAT=',iostat,' while closing unit',iou
      select case ( on_error )
      case ( IO_ABORT )
        ABORT('file_close_fortran: IO error')
      case ( IO_WARN )
        WARN('file_close_fortran: ignore IO error')
      case ( IO_SILENT )
        ! continue silently
      case default
        ABORT('no such case')
      end select
    endif

    ! return io unit:
    call return_iounit(iou)
    iou = iostat
  end subroutine file_close
#endif

  subroutine file_write_file_real_buf(fn,buf,siz)
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(in)  :: buf(*)
    integer(IK)     , intent(in)  :: siz
    ! *** end of interface ***

    integer(IK) :: iou

    DPRINT 'file_write_file_real_buf: ',trim(fn),siz

    if(siz.eq.0) RETURN

    call file_open(iou,fn,'w')

    if( iou < 0 )then
      print *,'file_write_file_real_buf: IO error: got iou=',iou,' while opening ',trim(fn)
      select case ( on_error )
      case ( IO_ABORT )
        ABORT('file_write_file_real_buf: IO error')
      case ( IO_WARN )
        WARN('file_write_file_real_buf: ignore IO error')
        GOTO 999
      case ( IO_SILENT )
        ! continue silently
        GOTO 999
      case default
        ABORT('no such case')
      end select
    endif

    call file_write_real_buf(iou,buf,siz)
    call file_close(iou,'k')
999 CONTINUE
  end subroutine file_write_file_real_buf

  subroutine file_read_file_real_buf(fn,buf,siz)
    use iounitadmin_module, only: get_iounit, return_iounit
    implicit none
    character(len=*), intent(in)  :: fn
    real(RK)        , intent(out) :: buf(*)
    integer(IK)     , intent(in)  :: siz
    ! *** end of interface ***

    integer(IK) :: iou

    DPRINT 'file_read_file_real_buf : <',trim(fn),siz

    if(siz.eq.0) RETURN

    call file_open(iou,fn,'r')

    if( iou < 0 )then
      print *,'file_read_file_real_buf: IO error: got iou=',iou,' while opening ',trim(fn)
      select case ( on_error )
      case ( IO_ABORT )
        ABORT('file_read_file_real_buf: IO error')
      case ( IO_WARN )
        WARN('file_read_file_real_buf: ignore IO error')
        GOTO 999
      case ( IO_SILENT )
        ! continue silently
        GOTO 999
      case default
        ABORT('no such case')
      end select
    endif

    call file_read_real_buf(iou,buf,siz)
    call file_close(iou,'k')
999 CONTINUE
  end subroutine file_read_file_real_buf

  !--------------- End of module -------------------------------------
end module io

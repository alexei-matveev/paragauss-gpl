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
!===============================================================
! Public interface of module
!===============================================================
module error_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
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
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private


  !------------ Interface statements ------------------------------
  interface error
     module procedure error_wrapper
     !
     module procedure error_logical
     module procedure error_i4_kind
  end interface

  interface warn
     module procedure warn_message
  end interface

  !------------ public functions and subroutines ------------------

  public init_error_module,&
       & done_error_module,&
       & error,&
       & warn


  integer(i4_kind)        :: MyIndex = -1
  character(len=6),public :: MyID    = "[????]"

  logical                 :: initialized = .false.


contains

  subroutine init_error_module(Rank)
    !use comm_module
    implicit none
    integer, intent(in) :: Rank
    ! *** end of interface ***

    character(len=4) :: buf

    !yIndex = comm_myindex()
    MyIndex = Rank

    write(buf,'(I4)') MyIndex

    MyID = "["//buf//"]"

    initialized = .true.
  end subroutine init_error_module

  subroutine done_error_module()
    implicit none
    ! *** end of interface ***

    initialized = .false.
  end subroutine done_error_module

  subroutine warn_message(msg)
    implicit none
    character(LEN=*),intent(in) :: msg
    ! *** end of interface ***

    print *,MyID,msg
  end subroutine warn_message

  subroutine error_wrapper(msg,file,line)
    implicit none
    external error_handler
    character(LEN=*), intent(in), optional :: msg
    character(LEN=*), intent(in), optional :: file
    integer(i4_kind), intent(in), optional :: line
    ! *** end of interface ***

    character(LEN=256) :: buf

    buf = repeat(' ',LEN(buf))
    if ( present(line) ) then
       write(buf,'(I8)') line
       buf = " line: "//trim(adjustl(buf))
    endif
    if ( present(file) ) then
       buf = " file: "//trim(file)//trim(buf)
    endif
    if ( present(msg) ) then
       buf = trim(msg)//trim(buf)
    endif

    call error_handler(MyID//"Error: "//trim(buf))
  end subroutine error_wrapper

  subroutine error_logical(err,msg)
    implicit none
    logical,intent(in)          :: err
    character(LEN=*),intent(in) :: msg
    ! *** end of interface ***

    if(err)call error(MyID//"errm/el: "//msg)
  end subroutine error_logical

  subroutine error_i4_kind(err,msg)
    implicit none
    integer(I4_kind),intent(in) :: err
    character(LEN=*),intent(in) :: msg
    ! *** end of interface ***

    if(err/=0)call error(MyID//"errm/ei: "//msg)
  end subroutine error_i4_kind

end module error_module

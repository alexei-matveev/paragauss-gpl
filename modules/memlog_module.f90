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
module memlog_module
  !-------------------------------------------------------------------
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
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: meminc,memset
  public :: memusage ! returns module variable ``usage''
  public :: memshow  ! dumps stats on tty

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------
 
  !------------ Declaration of constants and variables ---------------
  integer(IK) ::&
       usage          = 0, &
       max_usage      = 0, &
       max_usage_glob = 0

  integer(IK) :: last_usage = 0

  real(RK), parameter :: interval = 120.0
  real(RK)            :: time
  integer(IK)         :: next_count=-2

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains


  !*************************************************************
  subroutine meminc(inc)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use options_module, only: options_debug_key ! <- backdoor
    use time_module, only: clktime
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in) :: inc
    !** End of interface *****************************************

    time = clktime()

    usage = usage + inc
    if( usage > max_usage )then
       max_usage = usage
       ! print *, 'MEMLOG: max_usage=',max_usage*8,' b'
    endif

    if ( options_debug_key(2**30) /= 0 ) then
      call memshow(1)
    else
      if( time > next_count * interval )then
        next_count = INT( time / interval ) + 1
        call memshow(2)
      endif
    endif
  end subroutine meminc
  !*************************************************************

  function memusage() result(u)
    implicit none
    integer(IK) :: u ! result
    ! *** end of interface ***

    u = usage
  end function memusage

  subroutine memshow(level)
    use error_module, only: MyID
    use time_module, only: clktime
    implicit none
    integer(IK), intent(in) :: level
    ! *** end of interface ***

    if( level > 0 )then
      print *, MyID, 'MEMLOG(STATUS):'       &
           , 'usage=',  usage*8              &
           , 'change=', (usage-last_usage)*8 &
           , 'peak =',  max_usage*8          &
           , 'b'                             &
           , ' time=',clktime()
      last_usage = usage
    endif

    if( level > 1 )then
      call proc_self_status()
    endif

    if( level > 2 )then
      print *, MyID, 'MEMLOG(RESET):'          &
           , 'usage=',  usage*8                &
           , 'peak(loc) =',  max_usage*8       &
           , 'peak(glob) =',  max_usage_glob*8 &
           , 'b'
    endif
  end subroutine memshow

  !*************************************************************
  subroutine memset(val)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use error_module, only: MyID
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in) :: val
    !** End of interface *****************************************

    max_usage_glob = MAX(max_usage_glob,max_usage)

!   call display(3)

!   usage     = val
!   max_usage = val
!   call meminc(0) ! just for consistency
  end subroutine memset
  !*************************************************************

  !*************************************************************
  subroutine proc_self_status()
#ifdef _LINUX
     !  Purpose: cat /proc/self/status:
     !    VmSize:     2500 kB
     !    VmLck:         0 kB
     !    VmRSS:       536 kB
     !    VmData:      152 kB
     !    VmStk:         8 kB
     !    VmExe:        16 kB
     !    VmLib:      1152 kB
     !------------ Modules used ------------------- ---------------
     use error_module, only: MyID
     use iounitadmin_module, only: get_iounit, return_iounit
     implicit none
     !------------ Declaration of formal parameters ---------------
     !** End of interface *****************************************

     integer            :: u
     integer            :: n
     character(len=256) :: line
     integer            :: stat

     u = get_iounit()

     OPEN(  u                        &
          , FILE     = '/proc/self/status' &
          , FORM     = 'formatted'   &
          , ACTION   = 'read'        &
          , STATUS   = 'old'         &
          , POSITION = 'rewind'      &
          )
     do
       read(u,'(A)',iostat=stat) line
       if( stat /= 0 ) exit ! do-loop
       if( line(1:2) == 'Vm' )then
         n = index(line,':') - 1

         select case( line(1:n) )
         case ( 'VmSize', 'VmRSS', 'VmData' ) !, 'VmStk' )
           print *, MyID, 'MEMLOG(STATUS):',trim(line)
         case default
           ! do nothing
         end select
       endif
     enddo
     CLOSE(u)
     call return_iounit(u)
#endif
  end subroutine proc_self_status
  !*************************************************************

  !--------------- End of module -------------------------------------
end module memlog_module

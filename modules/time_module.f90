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
module time_module
!
!  Purpose: timing of program
!
!  Type timer_type is defined to time one specific task
!  within program. If this task is executed several times,
!  timer_type also keeps track of sum of time required for
!  all executions of this task, number of executions and
!  average time required per execution. Define an instance
!  of type timer_type ("timer") for all tasks you want to time.
!  Call subroutine init_timer to initialize timer.
!  Call subroutine start_timer before execution of task
!  and subroutine stop_timer after execution. Use subroutine
!  print_timer to print timing of task.
!
!  Use subroutine print_time to print time since start of program.
!
!  Subroutine init_time must be called at start of program.
!
!  All timers defined should be contained in timers_module
!
!
!  Author: TB
!  Date: 9/95
!
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Author: TG, TB
! Date:   2/96
! Description: Complete reorganisation of timing concept and module
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

#ifdef FPP_AIX_XLF
# define ETIME  etime_
#endif

#if defined(_SGI) || defined(_SR8000) || defined(_LINUX) || defined(_DEC) || defined(_COMPAC_FORTRAN)
# ifdef _UNDERSCORE
#  define SECNDS secnds_
#  define ETIME etime_
# else
#  define SECNDS secnds
#  define ETIME  etime
# endif
#endif

#if defined(_SGI) || defined(_SR8000) || defined(_DEC) || defined(_COMPAC_FORTRAN)
# define SYSTEM_CLOCK
# define ITIME
#endif

!------------ Modules used --------------------------------------
#define FPP_TIMERS 2
#include "def.h"
use type_module ! type specification parameters

implicit none
save
private

!== Interrupt end of public interface of module ================

!------------ Declaration of public types -----------------------

type, public :: timer_type
   ! to contain all time information for one task to be timed
   real(kind=r8_kind) :: real_timediff, usr_timediff, sys_timediff, &
                         real_starttime, usr_starttime, sys_starttime, &
                         real_timesum, usr_timesum, sys_timesum, &
                         real_averagetime, usr_averagetime, sys_averagetime, &
                         real_absoluttime, usr_absoluttime, sys_absoluttime
   ! timediff: timing from last call to start_timer to last call to stop_timer
   ! timesum: sum of all timediffs at last call to stop_timer
   ! averagetime: average off all timediffs at last call to stop_timer
   ! absoluttime: timing since call to init_time at last call to stop_timer
   ! starttime: timing at last call of start_timer, do not modify
   ! all timings are given in secounds
   integer :: nbr_of_calls ! to start_timer, do not modify
   logical :: running ! for error checks, do not modify
   character(len=64)  :: name
end type timer_type


!------------ public functions and subroutines ------------------

public :: usrtime ! () -> real(r8_kind), since time_setup()
public :: clktime ! () -> real(r8_kind), since time_setup()

public :: &
     time_setup, &
     init_timer,&
     start_timer,&
     stop_timer, &
     print_timer,&
     print_time,&
     timer_small_to_large,&
     add_timer,&
     divide_timer,&
     pack_timer,&
     unpack_timer,&
     unpack_add_timer,&
     switch_timer,&
     init_timer_v2


!================================================================
! End of public interface of module
!================================================================


!------------ Declaration of private variables ----------------

#if defined(_SGI) || defined(_SR8000) || defined(_DEC) || defined(_COMPAC_FORTRAN)
integer, private :: count, old_count, count_rate, count_max, &
     old_st, st, old_time(3), new_time(3), iday,iyear,imon, &
     iday_old, iyear_old, imon_old
real(kind=r8_kind), private :: st_per_cycle, real_count_rate
#endif

real(kind=r8_kind), private :: real_time, usr_time, sys_time
real(kind=r8_kind), private :: real_time_at_start

! see ../include/def.h for macros:
FPP_USR_TIMER_DECL(usr)
FPP_CLK_TIMER_DECL(clk)



!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

   function usrtime() result(time)
      implicit none
      real(r8_kind) :: time
      ! *** end of interface ***

      ! see ../include/def.h for macros:
      FPP_USR_TIMER_STOP(usr)
      time = FPP_USR_TIMER_VALUE(usr)
      FPP_USR_TIMER_START(usr)
   end function usrtime

   function clktime() result(time)
      implicit none
      real(r8_kind) :: time
      ! *** end of interface ***

      ! see ../include/def.h for macros:
      FPP_CLK_TIMER_STOP(clk)
      time = FPP_CLK_TIMER_VALUE(clk)
      FPP_CLK_TIMER_START(clk)
   end function clktime

   !*************************************************************
   subroutine time_setup()
   !  Purpose: should be called at start of program to start clock
   !** End of interface *****************************************
   !------------ Modules used -----------------------------------
#ifdef FPP_AIX_XLF
   use xlfutility
#endif
   implicit none
#if defined(_SGI) || defined(_SR8000) || defined(_DEC) || defined(_COMPAC_FORTRAN)
   external time
   integer :: time
#endif
#ifndef INTRINSIC_SECNDS
    real, external :: SECNDS
#endif
   !------------ Declaration of local variables -----------------
   !------------ Executable code --------------------------------
#ifdef FPP_AIX_XLF
    real_time      = timef()
#endif
#ifdef SECNDS
    real_time          = real(SECNDS(0.0),r8_kind)
    real_time_at_start = real_time
    real_time = real_time - real_time_at_start ! ZERO
#endif
#ifdef SYSTEM_CLOCK
    call system_clock(old_count, count_rate, count_max)
    real_count_rate = real(count_rate ,kind=r8_kind)
    st_per_cycle = count_max / real_count_rate
#endif
#ifdef ITIME
    call itime(old_time)
    call idate(imon,iday,iyear)
#endif

    ! see ../include/def.h for macros:
    FPP_USR_TIMER_START(usr)
    FPP_CLK_TIMER_START(clk)
   end subroutine time_setup
   !*************************************************************


   !*************************************************************
   subroutine fetch_timing()
   !  Purpose: fetches timing from system clock and stores it in
   !           private vars real_time, usr_time, sys_time
   !** End of interface *****************************************
#ifdef FPP_AIX_XLF
   use xlfutility
#endif
   implicit none
#if defined(_SGI) || defined(_SR8000) || defined(_DEC) || defined(_COMPAC_FORTRAN)
   external itime, idate
   integer :: warp
#endif
#ifndef INTRINSIC_SECNDS
   real, external :: SECNDS
#endif
#ifndef INTRINSIC_ETIME
   real, external :: ETIME
#endif
   !------------ Declaration of local variables -----------------
#ifdef FPP_AIX_XLF
   type(tb_type) buffer
#endif
#if defined(_SGI) || defined(_SR8000) || defined(_DEC) || defined(_COMPAC_FORTRAN)
   real :: real_sec
#endif
#if defined(_SGI) || defined(_LINUX) || defined(_SR8000)
   real*4 :: tsum, tarray(2)
#endif
   !------------ Executable code --------------------------------
#ifdef FPP_AIX_XLF
   real_time      = real(ETIME(buffer) ,kind=r8_kind)
   real_time      = timef() * 1.0E-3_r8_kind
   usr_time       = real(buffer%usrtime ,kind=r8_kind)
   sys_time       = real(buffer%systime ,kind=r8_kind)
#endif
#ifdef SYSTEM_CLOCK
   call system_clock(count)
   real_time = real(count-old_count,kind=r8_kind) / real_count_rate
#endif
#ifdef ITIME
   call itime(new_time)
   call idate(imon,iday,iyear)
   real_sec=real((new_time(1)-old_time(1))*3600+(new_time(2)-old_time(2))*60+&
        (iday-iday_old)*86400+new_time(3)-old_time(3))

   ! what is this magic for?:
   warp=int(real(real_sec-real_time)/st_per_cycle+0.5)
   real_time=real_time+warp*st_per_cycle
#endif
#ifdef SECNDS
    real_time=real(SECNDS(0.0),r8_kind) - real_time_at_start
#endif
#ifdef ETIME
   tsum = ETIME(tarray)
   usr_time = real(tarray(1) ,kind=r8_kind)
   sys_time = real(tarray(2) ,kind=r8_kind)
#endif

      ! update timers (to avoid overruns in CLK)
      ! see ../include/def.h for macros:
      FPP_USR_TIMER_STOP(usr)
      FPP_USR_TIMER_START(usr)
      FPP_CLK_TIMER_STOP(clk)
      FPP_CLK_TIMER_START(clk)
   end subroutine fetch_timing
   !*************************************************************


   !*************************************************************
   subroutine stop_timer( tt )
   !  Purpose: stops timer and returns current timing
   !------------ Modules used -----------------------------------
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type), intent(inout) :: tt
   !** End of interface *****************************************
   DPRINT 'stop_timer: '//trim(tt%name)
   if ( .not. tt%running ) call error_handler( &
        "stop_timer: timer not running: "//trim(tt%name) )
   call fetch_timing()
   tt%real_timediff  = real_time - tt%real_starttime
   tt%sys_timediff  = sys_time - tt%sys_starttime
   tt%usr_timediff  = usr_time - tt%usr_starttime
   tt%real_absoluttime = real_time
   tt%sys_absoluttime  = sys_time
   tt%usr_absoluttime  = usr_time
   tt%real_timesum  = tt%real_timesum + tt%real_timediff
   tt%usr_timesum  = tt%usr_timesum + tt%usr_timediff
   tt%sys_timesum  = tt%sys_timesum + tt%sys_timediff
   tt%real_averagetime  = tt%real_timesum / tt%nbr_of_calls
   tt%usr_averagetime  = tt%usr_timesum / tt%nbr_of_calls
   tt%sys_averagetime  = tt%sys_timesum / tt%nbr_of_calls
   tt%running          = .false.
   end subroutine stop_timer
   !*************************************************************


   !*************************************************************
   subroutine start_timer( tt )
   !  Purpose: starts timer
   !------------ Modules used -----------------------------------
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type), intent(out) :: tt
   !** End of interface *****************************************
   DPRINT 'start_timer: '//trim(tt%name)
   if ( tt%running ) call error_handler( &
        "start_timer: timer already running: "//trim(tt%name) )
   call fetch_timing()
   tt%real_starttime = real_time
   tt%usr_starttime  = usr_time
   tt%sys_starttime  = sys_time
   tt%nbr_of_calls   = tt%nbr_of_calls + 1
   tt%running        = .true.
   end subroutine start_timer
   !*************************************************************


   !*************************************************************
   subroutine switch_timer( tt_stop, tt_start )
   !  Purpose: stops timer tt_stop and starts timer tt_start
   !------------ Modules used -----------------------------------
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type), intent(inout) :: tt_stop, tt_start
   !** End of interface *****************************************
   DPRINT 'switch_timer: '//trim(tt_stop%name)//"<>"//trim(tt_start%name)
   if ( .not. tt_stop%running ) call error_handler( &
        "stop_timer: timer not running" )
   call fetch_timing()
   tt_stop%real_timediff  = real_time - tt_stop%real_starttime
   tt_stop%sys_timediff  = sys_time - tt_stop%sys_starttime
   tt_stop%usr_timediff  = usr_time - tt_stop%usr_starttime
   tt_stop%real_absoluttime = real_time
   tt_stop%sys_absoluttime  = sys_time
   tt_stop%usr_absoluttime  = usr_time
   tt_stop%real_timesum  = tt_stop%real_timesum + tt_stop%real_timediff
   tt_stop%usr_timesum  = tt_stop%usr_timesum + tt_stop%usr_timediff
   tt_stop%sys_timesum  = tt_stop%sys_timesum + tt_stop%sys_timediff
   tt_stop%real_averagetime  = tt_stop%real_timesum / tt_stop%nbr_of_calls
   tt_stop%usr_averagetime  = tt_stop%usr_timesum / tt_stop%nbr_of_calls
   tt_stop%sys_averagetime  = tt_stop%sys_timesum / tt_stop%nbr_of_calls
   tt_stop%running          = .false.
   if ( tt_start%running ) call error_handler( &
        "switch_timer: timer already running: "//trim(tt_start%name) )
   tt_start%real_starttime = real_time
   tt_start%usr_starttime  = usr_time
   tt_start%sys_starttime  = sys_time
   tt_start%nbr_of_calls   = tt_start%nbr_of_calls + 1
   tt_start%running        = .true.
   end subroutine switch_timer
   !*************************************************************


   !*************************************************************
   subroutine print_timer( tt, io_unit, message, &
        diff, absolut, start, sum, average, nbr, empty )
   !  Purpose: prints out information in tt and message
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),        intent(in) :: tt
   integer(kind=i4_kind),   intent(in) :: io_unit
   character*(*), optional, intent(in) :: message
   !  to decide witch times to print:
   logical, optional, intent(in) :: diff    ! default: .true.
   logical, optional, intent(in) :: absolut ! default: .true.
   logical, optional, intent(in) :: start   ! default: .false.
   logical, optional, intent(in) :: sum     ! default: .false.
   logical, optional, intent(in) :: average ! default: .false.
   logical, optional, intent(in) :: nbr ! of calls, default: .false.
   logical, optional, intent(in) :: empty ! default: .false.
   ! print if tt%nbr_of_calls .eq. 0
   !** End of interface *****************************************
   logical :: print_diff, print_absolut, print_start, print_sum, &
        print_average, print_nbr, print_empty
   !------------ Executable code --------------------------------
   DPRINT 'print_timer: '//trim(tt%name)
   if ( present(empty) ) then
      print_empty = empty
   else
      print_empty = .false.
   endif
   if ( .not. print_empty .and. tt%nbr_of_calls .eq. 0 ) return

   if (tt % running) then
      call error_handler ("print_timer: timer running " // tt % name)
   endif

   if ( present(diff) ) then
      print_diff = diff
   else
      print_diff = .true.
   endif
   if ( present(absolut) ) then
      print_absolut = absolut
   else
      print_absolut = .true.
   endif
   if ( present(start) ) then
      print_start = start
   else
      print_start = .false.
   endif
   if ( present(sum) ) then
      print_sum = sum
   else
      print_sum = .false.
   endif
   if ( present(average) ) then
      print_average = average
   else
      print_average = .false.
   endif
   if ( present(nbr) ) then
      print_nbr = nbr
   else
      print_nbr = .false.
   endif

   write(io_unit,*)
   if ( present(message) ) then
      write(io_unit,*) 'Timing:  ', message
   else
      write(io_unit,*) 'Timing'
   endif
   if ( print_absolut ) &
        write(io_unit,'("   since start of program           :  real  ",F10.3,"  usr  ",F10.3,"  sys  ",F10.3)') &
        tt%real_absoluttime, tt%usr_absoluttime, tt%sys_absoluttime
   if ( print_diff ) &
        write(io_unit,'("   for task (present execution)     :  real  ",F10.3,"  usr  ",F10.3,"  sys  ",F10.3)') &
        tt%real_timediff, tt%usr_timediff, tt%sys_timediff
   if ( print_start ) &
        write(io_unit,'("   at start of task                 :  real  ",F10.3,"  usr  ",F10.3,"  sys  ",F10.3)') &
        tt%real_starttime, tt%usr_starttime, tt%sys_starttime
   if ( print_sum ) &
        write(io_unit,'("   for task (all executions)        :  real  ",F10.3,"  usr  ",F10.3,"  sys  ",F10.3)') &
        tt%real_timesum, tt%usr_timesum, tt%sys_timesum
   if ( print_average ) &
        write(io_unit,'("   for task (average per executions):  real  ",F10.3,"  usr  ",F10.3,"  sys  ",F10.3)') &
        tt%real_averagetime, tt%usr_averagetime, tt%sys_averagetime
   if ( print_nbr ) &
        write(io_unit,'("   number of executions of task     :  ",I5)') tt%nbr_of_calls
   write(io_unit,*)
   end subroutine print_timer
   !*************************************************************


   !*************************************************************
   subroutine init_timer( tt )
   !  Purpose: sets times in tt to 0.0
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(inout) :: tt
   !** End of interface *****************************************
   !------------ Executable code --------------------------------
   tt%real_timediff    = 0.0_r8_kind
   tt%usr_timediff     = 0.0_r8_kind
   tt%sys_timediff     = 0.0_r8_kind
   tt%real_absoluttime = 0.0_r8_kind
   tt%sys_absoluttime  = 0.0_r8_kind
   tt%usr_absoluttime  = 0.0_r8_kind
   tt%real_timesum     = 0.0_r8_kind
   tt%usr_timesum      = 0.0_r8_kind
   tt%sys_timesum      = 0.0_r8_kind
   tt%real_averagetime = 0.0_r8_kind
   tt%usr_averagetime  = 0.0_r8_kind
   tt%sys_averagetime  = 0.0_r8_kind
   tt%real_starttime   = 0.0_r8_kind
   tt%usr_starttime    = 0.0_r8_kind
   tt%sys_starttime    = 0.0_r8_kind
   tt%nbr_of_calls     = 0
   tt%running          = .false.
   tt%name             = "undefined"
   end subroutine init_timer
   !*************************************************************

   !*************************************************************
   subroutine init_timer_v2( tt , name)
     !  Purpose: sets times in tt to 0.0
     implicit none
     !------------ Declaration of formal parameters ---------------
     type(timer_type),       intent(inout) :: tt
     character(len=*), intent(in)          :: name
     !** End of interface *****************************************
     !------------ Executable code --------------------------------
     DPRINT 'init_timer_v2: '//trim(name)
     tt%real_timediff    = 0.0_r8_kind
     tt%usr_timediff     = 0.0_r8_kind
     tt%sys_timediff     = 0.0_r8_kind
     tt%real_absoluttime = 0.0_r8_kind
     tt%sys_absoluttime  = 0.0_r8_kind
     tt%usr_absoluttime  = 0.0_r8_kind
     tt%real_timesum     = 0.0_r8_kind
     tt%usr_timesum      = 0.0_r8_kind
     tt%sys_timesum      = 0.0_r8_kind
     tt%real_averagetime = 0.0_r8_kind
     tt%usr_averagetime  = 0.0_r8_kind
     tt%sys_averagetime  = 0.0_r8_kind
     tt%real_starttime   = 0.0_r8_kind
     tt%usr_starttime    = 0.0_r8_kind
     tt%sys_starttime    = 0.0_r8_kind
     tt%nbr_of_calls     = 0
     tt%running          = .false.
     tt%name             = name
   end subroutine init_timer_v2
   !*************************************************************


   !*************************************************************
   subroutine timer_small_to_large( tt_small, tt_large )
   !  Purpose: Intended for situations when multiple calls of
   !  timer start and stop, as within grid loop, measured with timer
   !  tt_small, should be transfered as one call to timer tt_large.
   !  Makes the sum timings off tt_small to one timing call to
   !  timer tt_large, modifying all its contents
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(in) :: tt_small
   type(timer_type),       intent(inout) :: tt_large
   !** End of interface *****************************************
   !------------ Executable code --------------------------------
   if ( tt_small%running ) call error_handler( &
        "small_to_large: timer small running" )
   if ( tt_large%running ) call error_handler( &
        "small_to_large: timer large running" )
   tt_large%real_timediff    = tt_small%real_timesum
   tt_large%usr_timediff     = tt_small%usr_timesum
   tt_large%sys_timediff     = tt_small%sys_timesum
   tt_large%real_absoluttime = tt_small%real_absoluttime
   tt_large%sys_absoluttime  = tt_small%sys_absoluttime
   tt_large%usr_absoluttime  = tt_small%usr_absoluttime
   tt_large%real_timesum     = tt_large%real_timesum + tt_small%real_timesum
   tt_large%usr_timesum      = tt_large%usr_timesum + tt_small%usr_timesum
   tt_large%sys_timesum      = tt_large%sys_timesum + tt_small%sys_timesum
   tt_large%nbr_of_calls     = tt_large%nbr_of_calls + 1
   tt_large%real_averagetime = tt_large%real_timesum / tt_large%nbr_of_calls
   tt_large%usr_averagetime  = tt_large%usr_timesum / tt_large%nbr_of_calls
   tt_large%sys_averagetime  = tt_large%sys_timesum / tt_large%nbr_of_calls
   tt_large%real_starttime   = 0.0_r8_kind
   tt_large%usr_starttime    = 0.0_r8_kind
   tt_large%sys_starttime    = 0.0_r8_kind
   end subroutine timer_small_to_large
   !*************************************************************


   !*************************************************************
   subroutine add_timer( tt_summand, tt_sum )
   !  Purpose: Adds the timesums and updates tt_sum%nbr_of_calls
   !           and tt_sum%xxx_averagetime accordingley.
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(in)    :: tt_summand
   type(timer_type),       intent(inout) :: tt_sum
   !** End of interface *****************************************
   !------------ Executable code --------------------------------
   DPRINT 'add_timer: '//trim(tt_sum%name)//"+="//trim(tt_summand%name)
   if ( tt_summand%running ) call error_handler( &
        "timer_add: timer summand running" )
   if ( tt_sum%running ) call error_handler( &
        "timer_add: timer sum running" )
   tt_sum%real_timesum     = tt_sum%real_timesum + tt_summand%real_timesum
   tt_sum%usr_timesum      = tt_sum%usr_timesum + tt_summand%usr_timesum
   tt_sum%sys_timesum      = tt_sum%sys_timesum + tt_summand%sys_timesum
   tt_sum%nbr_of_calls     = tt_sum%nbr_of_calls + tt_summand%nbr_of_calls
   if ( tt_sum%nbr_of_calls .gt. 0 ) then
      tt_sum%real_averagetime = tt_sum%real_timesum / tt_sum%nbr_of_calls
      tt_sum%usr_averagetime  = tt_sum%usr_timesum / tt_sum%nbr_of_calls
      tt_sum%sys_averagetime  = tt_sum%sys_timesum / tt_sum%nbr_of_calls
   endif
   end subroutine add_timer
   !*************************************************************


   !*************************************************************
   subroutine divide_timer( tt, divisor )
   !  Purpose: Divides the timesums and tt_sum%nbr_of_calls
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(inout) :: tt
   real(kind=r8_kind),     intent(in)    :: divisor
   !** End of interface *****************************************
   !------------ Executable code --------------------------------
   DPRINT 'divide_timer: '//trim(tt%name)
   if ( tt%running ) call error_handler( &
        "timer_add: timer running" )
   tt%real_timesum     = tt%real_timesum / divisor
   tt%usr_timesum      = tt%usr_timesum / divisor
   tt%sys_timesum      = tt%sys_timesum / divisor
   tt%nbr_of_calls     = int( tt%nbr_of_calls / divisor, i4_kind )
   end subroutine divide_timer
   !*************************************************************


   !*************************************************************
   subroutine print_time( io_unit, message )
   !  Purpose: prints out time since start of program and message
   implicit none
   !------------ Declaration of formal parameters ---------------
   integer(kind=i4_kind),   intent(in) :: io_unit
   character*(*), optional, intent(in) :: message
   !** End of interface *****************************************
   !------------ Executable code --------------------------------
   call fetch_timing()
   write(io_unit,*)
   if ( present(message) ) then
      write(io_unit,*) 'Timing:  ', message
   else
      write(io_unit,*) 'Timing'
   endif
   write(io_unit,'("   since start of program           :  real  ",F10.3,"  usr  ",F10.3,"  sys  ",F10.3)') &
        real_time, usr_time, sys_time
   write(io_unit,*)
   end subroutine print_time
   !*************************************************************


   !*************************************************************
   subroutine pack_timer( tt )
   !  Purpose: packs contents of tt
   use comm_module, only: commpack
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(inout) :: tt
   !** End of interface *****************************************
   integer(kind=i4_kind) :: info
   !------------ Executable code --------------------------------
   call commpack(tt%real_timediff,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%real_timediff failed")
   call commpack(tt%usr_timediff,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%usr_timediff failed")
   call commpack(tt%sys_timediff,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%sys_timediff failed")
   call commpack(tt%real_absoluttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%real_absoluttime failed")
   call commpack(tt%sys_absoluttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%sys_absoluttime failed")
   call commpack(tt%usr_absoluttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%usr_absoluttime failed")
   call commpack(tt%real_timesum,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%real_timesum failed")
   call commpack(tt%usr_timesum,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%usr_timesum failed")
   call commpack(tt%sys_timesum,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%sys_timesum failed")
   call commpack(tt%real_averagetime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%real_averagetime failed")
   call commpack(tt%usr_averagetime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%usr_averagetime failed")
   call commpack(tt%sys_averagetime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%sys_averagetime failed")
   call commpack(tt%real_starttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%real_starttime failed")
   call commpack(tt%usr_starttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%usr_starttime failed")
   call commpack(tt%sys_starttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%sys_starttime failed")
   call commpack(tt%nbr_of_calls,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%nbr_of_calls failed")
   call commpack(tt%running,info)
   if ( info .ne. 0 ) call error_handler( &
        "pack_timer: packing of tt%running failed")
   end subroutine pack_timer
   !*************************************************************


   !*************************************************************
   subroutine unpack_timer( tt )
   !  Purpose: unpacks contents of tt
   use comm_module, only: communpack
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(inout) :: tt
   !** End of interface *****************************************
   integer(kind=i4_kind) :: info
   !------------ Executable code --------------------------------
   call communpack(tt%real_timediff,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%real_timediff failed")
   call communpack(tt%usr_timediff,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%usr_timediff failed")
   call communpack(tt%sys_timediff,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%sys_timediff failed")
   call communpack(tt%real_absoluttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%real_absoluttime failed")
   call communpack(tt%sys_absoluttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%sys_absoluttime failed")
   call communpack(tt%usr_absoluttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%usr_absoluttime failed")
   call communpack(tt%real_timesum,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%real_timesum failed")
   call communpack(tt%usr_timesum,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%usr_timesum failed")
   call communpack(tt%sys_timesum,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%sys_timesum failed")
   call communpack(tt%real_averagetime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%real_averagetime failed")
   call communpack(tt%usr_averagetime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%usr_averagetime failed")
   call communpack(tt%sys_averagetime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%sys_averagetime failed")
   call communpack(tt%real_starttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%real_starttime failed")
   call communpack(tt%usr_starttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%usr_starttime failed")
   call communpack(tt%sys_starttime,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%sys_starttime failed")
   call communpack(tt%nbr_of_calls,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%nbr_of_calls failed")
   call communpack(tt%running,info)
   if ( info .ne. 0 ) call error_handler( &
        "unpack_timer: unpacking of tt%running failed")
   end subroutine unpack_timer
   !*************************************************************


   !*************************************************************
   subroutine unpack_add_timer( tt )
   !  Purpose: unpacks timer and adds it to timer tt
   implicit none
   !------------ Declaration of formal parameters ---------------
   type(timer_type),       intent(inout) :: tt
   !** End of interface *****************************************
   type(timer_type) :: tt_local
   !------------ Executable code --------------------------------
   call unpack_timer( tt_local )
   call add_timer( tt_local, tt )
   end subroutine unpack_add_timer
   !*************************************************************


!--------------- End of module ----------------------------------
end module time_module

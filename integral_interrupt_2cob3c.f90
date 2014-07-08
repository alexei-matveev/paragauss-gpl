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
subroutine  integral_interrupt_2cob3c
!---------------------------------------------------------------------
!
!  Purpose: This routine is the interrupt routine for the
!           integral calculation of one quadrupel
!           ( unique atom 1, l 1, unique atom 2, l 2 )
!           for 2 center orbital and three center integral
!           calculation.
!           the routine tries to receive from any host the
!           following messages:
!             msgtag_int_2cob3c_result
!
!  Subroutine called by: integral_main_2cob3c, 
!     integral_calc_quad_2cob3c and routines called there
!
!  References: Publisher Document: Concepts of Integral Part
! 
!  Author: TB
!  Date: 6/96
!
  !===================================================================
  ! End of public interface of module
  !===================================================================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification DLB instead of master/slave for 2cob3c integrals
! Author: AN
! Date:   4/11
! Description: Integrals are now scheduled with DLB
!              Interrupt does not need to appear that often
!              no need to wait for messages about task requesting
!!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...

!---------------------------------------------------------------------

  !------------ Modules used ------------------------------------
  use type_module ! type specification parameters
  use comm_module,                 only: comm_parallel                         &
                                       , comm_all_other_hosts                  &
                                       , comm_save_recv_nonblocking
  use msgtag_module,               only: msgtag_int_2cob3c_result
  use iounitadmin_module, only: write_to_output_units
  use output_module, only: output_slaveoperations
  use timer_module,                only: timer_int_calc_2cob3c                 &
                                       , timer_int_interrupt_2cob3c
  use time_module, only: start_timer , stop_timer
  use integralpar_module, only: integralpar_i_int_part
  use quadrupel_module,            only: quadrupel_unpack                      &
                                       , quadrupel_type
  use int_send_2cob3c_module,      only: int_send_2cob3c_receive
  use int_send_2cob3c_spor_module, only: int_send_2cob3c_spor_receive
  use options_module,              only: options_spin_orbit

  implicit none

  !------------ Declaration of local variables ------------------
  logical                     :: restart_timer, timer_running

  if ( .not. comm_parallel() ) return

  timer_running = .false.


  ! try to receive message containing results
  do while (comm_save_recv_nonblocking(comm_all_other_hosts,msgtag_int_2cob3c_result))

     if (.not. timer_running) call start_interrupt_timing()

     if ( output_slaveoperations ) &
          call write_to_output_units("integral_interrupt_2cob3c: integral_2cob_result")
     if (options_spin_orbit) then
        !
        ! SPIN ORBIT
        !
        call int_send_2cob3c_spor_receive()
     else
        call int_send_2cob3c_receive()
     endif
  enddo
  
  if (timer_running) call stop_interrupt_timing()

contains

  subroutine start_interrupt_timing()
    if ( timer_int_calc_2cob3c(integralpar_i_int_part)%running ) then
       call stop_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
       restart_timer = .true.
    else
       restart_timer = .false.
    endif
    call start_timer(timer_int_interrupt_2cob3c(integralpar_i_int_part))
    timer_running = .true.
  end subroutine start_interrupt_timing


  subroutine stop_interrupt_timing()
    call stop_timer(timer_int_interrupt_2cob3c(integralpar_i_int_part))
    if ( restart_timer ) call start_timer(timer_int_calc_2cob3c(integralpar_i_int_part))
  end subroutine stop_interrupt_timing


end subroutine integral_interrupt_2cob3c

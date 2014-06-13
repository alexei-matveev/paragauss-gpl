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
subroutine  integral_main_2cff
!----------------------------------------------------------------
!
!  Purpose: This is the main routine of the 2 center
!           fitfunction part of the integral part
!
!  Subroutine called by: main_integral
!
!  References: Publisher Document: Concepts of Integral Part
! 
!
!  Author: TB
!  Date: 5/96
!
!
!================================================================
! End of public interface of module
!================================================================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   6/98
! Description: ...
!
! Modification Master/Slave concept to DLB
! Author: AN
! Date:   4/11
! Description: for scheduling the 2cff integrals DLB is used
!              it replaces the master/slave concept
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
#define FPP_TIMERS 2
#include "def.h"
use type_module          ! type specification parameters
use output_module        ! defines amount of output
use iounitadmin_module   ! to open output units
use comm_module           ! comm related information and routines
use msgtag_module
use integralpar_module   ! steering information for integral part
use int_distribute_module, only: int_distribute_start, int_distribute_next_job
use int_distribute_module, only: int_distribute_setup, int_distribute_2cff_run
use int_distribute_module, only: int_distribute_shutdown
use int_data_2cff_module, only: quadrupel
use int_send_2cff_module
use timer_module
use time_module
use fit_coeff_module, only: fit_coeff_bcast
use  calc3c_switches
implicit none

integer(kind=i4_kind) :: n_quads

!----------------------------------------------------------------
!------------ Executable code -----------------------------------

! ================ CONTEXT: EVERY PROC ================

if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cff: start")

if( comm_i_am_master() )then
! ================ CONTEXT: MASTER ONLY ================
call start_timer(timer_int_2cff(integralpar_i_int_part))
call start_timer(timer_int_idle_2cff(integralpar_i_int_part))


! prepare distribution of integral quadrupels
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cff: int_distribute_setup")
endif ! i am master
! ================ CONTEXT: EVERY PROC ================
! do setup work, starting the idle timer
n_quads = int_distribute_setup(int_distribute_2cff_run)

if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cff: setup")
call integral_setup_2cff()

if( comm_i_am_master() )then
! ================ CONTEXT: MASTER ONLY ================
! if integrals have to be sent ( in normal run, but not in gradient run )
! prepare sending
if(integralpar_send_2c) then
   ! prepare receiving of results
   if (output_int_detailedprogress) call write_to_output_units( &
        "main_integral_2cff: int_send_2cff_setup")
   call int_send_2cff_setup(n_quads)
end if

! ================ CONTEXT: EVERY PROC ================
! start calculation
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cff: int_distribute_start")
endif ! i am master
call int_distribute_start()

! Main loop, distribution of jobs will be done by int_distribute_module
! integral_calc_quad_2cff will calculate the quadrupel stored in
! int_data_2cff_module, which will be set by int_distribute_next_job
! if int_distribute_next_job returns False quadrupel is invalid
do while ( int_distribute_next_job(quadrupel))
    call integral_calc_quad_2cff()
    if (comm_i_am_master() .and. integralpar_send_2c) then
        ! master tries to collect results by a nonblocking routine
        call try_get_results()
    endif
enddo
! Finished calculating quadrupels

DPRINT 't_2c_total',FPP_TIMER_VALUE(t_2c_total)
DPRINT 't_calc_2c_dervs',FPP_TIMER_VALUE(t_calc_2c_dervs)
DPRINT 't_contract_2c_dervs',FPP_TIMER_VALUE(t_contract_2c_dervs)
DPRINT 't_cpks_grad_fitmat',FPP_TIMER_VALUE(t_cpks_grad_fitmat)

if( comm_i_am_master() )then
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cff: int_distribute_shutdown")
endif
call int_distribute_shutdown()

! ================ CONTEXT: MASTER ONLY ================
if( comm_i_am_master() )then
! shutdown int_send_module if necessary
if(integralpar_send_2c) then
   if (output_int_detailedprogress) call write_to_output_units( &
        "main_integral_2cff: int_send_2cff_shutdown")
   ! after the next routine master should have all results sended to him
   call int_send_2cff_shutdown()
end if
endif ! i am master

! stop timer
call stop_timer(timer_int_idle_2cff(integralpar_i_int_part))
call stop_timer(timer_int_2cff(integralpar_i_int_part))

! ================ CONTEXT: EVERY PROC ================

! now broadcast the norms which were stored by calling 
!    fit_coeff_store_norm()
! from
!    int_send_2cff_shutdown()
call fit_coeff_bcast()

if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cff: done")

!Module contains only explicit messages (send and receive given msgtag)
! does it is not needed to have a barrier here
contains

  subroutine try_get_results()
    !Purpose: master tries to get the results from the slaves
    !         uses nonblocking mpi
    !         Slaves send their results at the end of integral_calc_quad_2cff
    use comm_module, only: comm_save_recv_nonblocking
    implicit none
    do while (comm_save_recv_nonblocking &
                (comm_all_other_hosts, msgtag_int_2cff_result))
       call int_send_2cff_receive()
    enddo
  end subroutine try_get_results

end subroutine integral_main_2cff

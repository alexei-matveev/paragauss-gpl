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
subroutine  integral_main_2cob3c
!----------------------------------------------------------------
!
!  Purpose: This is the main routine of the 2 center orbital
!           and 3 center part of the integral part
!           executed by the master.
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
! Date:   6/97
! Description: ...
!
! Modification master/slave -> DLB
! Author: AN
! Date:   4/11
! Description: use DLB for scheduling of 2cob3c integrals
!              master/slave use the same code to run
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
! define FPP_TIMERS 3
#include "def.h"
use type_module          ! type specification parameters
use output_module        ! defines amount of output
use iounitadmin_module   ! to open output units
use comm_module           ! comm related information and routines
use integralpar_module   ! steering information for integral part
use int_distribute_module
use timer_module
use time_module
use int_data_2cob3c_module, only: quadrupel
use int_send_2cob3c_module, only: int_send_2cob3c_receive_all
#ifdef FPP_TIMERS
 use dlb, only: dlb_print_statistics
 use comm, only: comm_rank, comm_size, comm_reduce
 USE_MPI, only: MPI_WTIME
#endif

implicit none
integer :: n_quads
#ifdef FPP_TIMERS
    real(kind=r8_kind), allocatable :: times(:, :)
    integer(kind=i4_kind) :: ierr
    integer(kind=i4_kind) :: np, rank
#endif
FPP_TIMER_DECL(tot)
FPP_TIMER_DECL(loop)
FPP_TIMER_DECL(work)
FPP_TIMER_DECL(sched)

!----------------------------------------------------------------
!------------ Executable code -----------------------------------


! ================ CONTEXT: EVERY PROC ================
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: start")

if( comm_i_am_master() )then
! ================ CONTEXT: MASTER ONLY ================
call start_timer(timer_int_2cob3c(integralpar_i_int_part))
call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))

! prepare distribution of integral quadrupels
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: int_distribute_setup")
endif ! i am master
! ================ CONTEXT: EVERY PROC ================
! do setup work, starting the idle timer
n_quads = int_distribute_setup(int_distribute_2cob3c_run)

if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: setup")
call integral_setup_2cob3c(n_quads)

if( comm_i_am_master() )then
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: int_distribute_start")
endif ! i am master

! perform calculation
! last preparations for starting the calculation
call int_distribute_start()

#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif
FPP_TIMER_START(tot)

! Main loop over the calculations, int_distribute_next_job is true
! as long the given processor has something to calculate, it puts
! the quadrupel to calculate into quadrupel, which should be the
! one stored in int_data_2cob3c_module, as this one is calculated
! in integral_calc_quad_2cob3c, if int_distribute_next_job returns
! false quadrupel contains invalid input
FPP_TIMER_START(loop)
FPP_TIMER_START(sched)
do while ( int_distribute_next_job(quadrupel))
    FPP_TIMER_STOP(sched)
    FPP_TIMER_START(work)
    call integral_calc_quad_2cob3c()
    if (integralpar_send_3c) then
        ! with this option the processors share their results, they
        ! send it at the end of integral_calc_quad_2cob3c, here
        ! they check (nonblocking) if anything of that kind has already
        ! arrived
        call integral_interrupt_2cob3c()
    endif
    FPP_TIMER_STOP(work)
    FPP_TIMER_START(sched)
enddo
#if FPP_TIMERS
    call comm_barrier()  ! for ensuring clean time measurments
#endif

FPP_TIMER_STOP(sched)
FPP_TIMER_STOP(loop)
! Here all results are received, if no results are exchanged, it
! will notice it by having not to wait
call int_send_2cob3c_receive_all()
FPP_TIMER_STOP(tot)

! do shutdown work
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: shutdown")
call integral_shutdown_2cob3c()

if( comm_i_am_master() )then
if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: int_distribute_shutdown")
endif ! i am master

call int_distribute_shutdown()


if( comm_i_am_master() )then
call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
call stop_timer(timer_int_2cob3c(integralpar_i_int_part))
endif ! i am master

#ifdef FPP_TIMERS
    np = comm_size()
    rank = comm_rank()
    allocate(times(0:np-1, 6), stat=ierr)
    ASSERT(ierr == 0)

    times = 0
    times(rank, 1) = FPP_TIMER_VALUE(tot)
    times(rank, 2) = FPP_TIMER_VALUE(loop)
    times(rank, 3) = FPP_TIMER_VALUE(work)
    times(rank, 4) = FPP_TIMER_VALUE(sched)
    times(rank, 5) = FPP_TIMER_SLICE(sched)

    call comm_reduce(times)
    call comm_barrier()  ! for ensuring clean output

    if(comm_i_am_master()) then
        print *, ""
        print *, "TIMINGS STATISTICS FOR LOAD BALANCING 2COB3C"
        print *, "All times are in seconds"
        print *, "Elapsed times all          :", maxval(times(:, 1)) * np ! total
        print *, "|- Elapsed times for loop  :", maxval(times(:, 2)) * np ! loop1
        print *, " |- Work in loop           :", sum(times(:, 3)) ! work2
        print *, " |- Scheduling in loop     :", sum(times(:, 4) - times(:, 5)) ! sched2 w/o last interval
        print *, " |- Rest at end of loop    :", sum(times(:, 5)) - np * minval(times(:, 5)) !sched2 last interval- term
        print *, "  \- Varianz of loop end   :", (maxval(times(:, 5)) - minval(times(:, 5))) * np ! sched2 last ineterval max -min
        print *, " |- Termination after loop :", minval(times(:, 5)) * np ! sched2 last interval
        print *, " Estimated efficiency      :", sum(times(:, 3)) / (maxval(times(:, 2)) * np) ! work2 / loop2
        print *, ""

    endif

    deallocate(times, stat = ierr)
    ASSERT(ierr == 0)
    call dlb_print_statistics(0)
    call comm_barrier()  ! for ensuring clean output
#endif



if (output_int_detailedprogress) call write_to_output_units( &
     "main_integral_2cob3c: done")

!As in this routine messages are received only by explicit message tags
! it is not needed to wait here for each other with the fear to get a
! message to early
end subroutine integral_main_2cob3c

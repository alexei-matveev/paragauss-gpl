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
module dlb_impl

!---------------------------------------------------------------
!
!  Takes care  about dynamical load balancing. Is  the static routine,
!  meaning, not dynamical, but  for comparision with them with exactly
!  the same interface.
!
!  INTERFACE to others:
!
!  call dlb_init()
!
!  Called once before first acces.
!
!  call dlb_setup(init_job)
!
!  Called once every time a dlb should be used, init_job should be the
!  part of the current proc of an initial job distribution.
!
!  dlb_give_more(n, jobs)
!
!  n should  be the number  of jobs requested  at once. The  next time
!  dlb_give_more()  is called  again,  all jobs  from  jobs should  be
!  finished. Jobs are at most n, they are just the number of the jobs,
!  which still have to be transformed between each other. It should be
!  done the error slice from jobs(0) + 1 to jobs(1) if dlb_give_more()
!  returns jobs(0) == jobs(1) there are on no proc any jobs left
!
!  Module called by: ...
!
!
!  References:
!
!  Author: AN
!  Date: 09/10
!
!
!---------------------------------------------------------------------
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!---------------------------------------------------------------------
# include "dlb.h"
USE_MPI, only: MPI_THREAD_SINGLE
USE_MPI, only: MPI_Wtime
use dlb_common, only: ik, lk, JLENGTH
use dlb_common, only: timer_give_more, timer_give_more_last
implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private

public :: dlb_init, dlb_finalize, dlb_setup, dlb_give_more !for using the module

! Program from outside might want to know the thread-safety-level required form DLB
integer (ik), parameter, public :: DLB_THREAD_REQUIRED = MPI_THREAD_SINGLE

!== Interrupt end of public interface of module =================

integer (lk) :: job_storage(JLENGTH) ! store all the jobs, belonging to this processor

contains

  subroutine dlb_init(world)
    ! Initalization of needed stuff
    use dlb_common, only: set_empty_job, dlb_common_init
    use dlb_common, only: my_rank, time_stamp
    implicit none
    integer, intent(in) :: world
    ! *** end of interface ***

    job_storage = set_empty_job()

    ! this also sets my_rank in dlb_common:
    call dlb_common_init(world)

    if ( my_rank == 0 ) then
        call time_stamp("dlb_init: using variant 'static'", output_level=0)
    endif
  end subroutine dlb_init

  subroutine dlb_finalize()
    ! Cleaning up everything, after last call
    use dlb_common, only: dlb_common_finalize
    implicit none
    ! *** end of interface ***

    call dlb_common_finalize()
  end subroutine dlb_finalize

  subroutine dlb_give_more(n, my_job)
    !
    ! Returns next  bunch of up to  n jobs, if  my_job(1) <= my_job(2)
    ! there are no more jobs there,  else returns the jobs done by the
    ! procs should  be my_job(1) + 1  to my_job(2) in  the related job
    ! list
    !
    use dlb_common, only: steal_local, JLEFT, JRIGHT, L_JOB
    use dlb_common, only: empty, time_stamp
    implicit none
    integer (lk), intent(in)  :: n
    integer (lk), intent(out) :: my_job(:)

    !** End of interface *****************************************
    double precision              :: start_timer_gm
    integer (lk) :: jobs(JLENGTH), remaining(JLENGTH)

    !
    ! Return of an empty job interval will be interpreted as
    ! "no jobs left", thus refuse requests for zero jobs:
    !
    start_timer_gm = MPI_Wtime() ! for debugging
    ASSERT(n>0)
    ASSERT(size(my_job)==L_JOB)

    if ( steal_local(n, job_storage, remaining, jobs) ) then
        job_storage = remaining
    endif

    !
    ! Named constants are not known outside:
    !
    my_job = jobs(:L_JOB)

    if ( empty(jobs) ) then
        ! output for debugging, only reduced inforamtions needed compared to
        ! dynamical cases
        timer_give_more_last = MPI_Wtime() - start_timer_gm ! for debugging
        call time_stamp("dlb_give_more: no jobs left", output_level=1)
    else
        timer_give_more = MPI_Wtime() - start_timer_gm ! for debugging
    endif
  end subroutine dlb_give_more

  subroutine dlb_setup(job)
    !
    ! Initialization  of a  DLB run,  each  proc should  call it  with
    ! inital jobs. The inital jobs  should be a static distribution of
    ! all available jobs, each job should be given to one and only one
    ! of the  procs jobs  should be  given as a  job range  (STP, EP),
    ! where all  jobs should  be the numbers  from START to  END, with
    ! START <= STP <= EP <= END
    !
    use dlb_common, only: JLEFT, JRIGHT, JOWNER, my_rank
    use dlb_common, only: time_stamp, L_JOB
    implicit none
    integer (lk), intent(in) :: job(:)
    !** End of interface *****************************************

    ASSERT(size(job)==L_JOB)

    job_storage(JLEFT) = job(JLEFT)
    job_storage(JRIGHT) = job(JRIGHT)
    job_storage(JOWNER) = my_rank ! FIXME: is this field ever used?

    call time_stamp("dlb_setup: done", output_level=1)
  end subroutine dlb_setup

  !--------------- End of module -------------------------------------
end module dlb_impl

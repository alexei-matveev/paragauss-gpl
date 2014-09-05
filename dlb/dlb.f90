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
module dlb
!---------------------------------------------------------------
!
!  Purpose: takes care about dynamical load balancing, there are two
!           different cases possible: one where all jobs are equal in
!           name and a second one where each job has a special color
!           It is assumed that several succeeding jobs have the same
!           color, thus the color is wanted in ranges, It is assured
!           that all the jobs given back have the same color, color
!           should be given as an integer There are several DLB
!           routines possible, which could be used to generate the
!           dynamical balancing themselves
!
!  Interface:
!
!           call dlb_init(world) - once before first access
!                 world is an integer specifying an MPI_WORLD
!
!           call dlb_finalize() - once after last access
!
!           There are two choices for the actual DLB run, they can be
!           used alternately but routines of them may not be mixed up.
!           First WITHOUT COLORS (e.g. all jobs are equal):
!
!           call dlb_setup(N) - once every time a DLB should
!              be used, N should be the number of jobs
!
!           dlb_give_more(n, jobs) - n should be the number of jobs
!              requested at once, the next time dlb_give_more is
!              called again, all jobs from jobs should be finished,
!              jobs are at most n, they are just the number of the
!              jobs, which still have to be transformed between each
!              other, it should be done the error slice from jobs(0)
!              +1 to jobs(1) if dlb_give_more returns jobs(0)==jobs(1)
!              there is no chance that this proc will get any job
!              during this DLB-round
!
!           Second WITH COLORS (e.g. each job has an integer number
!           attached to it, which gives the "color" of the job. It is
!           supposed that several succeeding jobs have the same color,
!           thus the color is given for job slices)
!
!           call dlb_setup_color(distr) - once every time a DLB should
!              be used, distr is an array containing elements like:
!              (/color, startnumber, endnumber/), it is possible to
!              have several not succeeding jobs with the same color,
!              then for each one the number has to be given in its own
!              array element startnumber and endnumber of each color
!              are independent of each other, they may have
!              overlapping intervals, in this case one has to ensure,
!              that the right jobs are done, as the DLB routine gives
!              back the numbers of the array
!
!           dlb_give_more_color(n, color, jobs) - the same as
!              dlb_give_more, but for the color case gives back the
!              numbers between startnumber and endnumber of the
!              element color, it only gives jobs of the same color,
!              even if it still has some left but from another color,
!              but it also will only give back empty jobs if all is
!              finished
!
!           Third WITH ROUND_ROBIN (e.g. all processes start with jobs
!           with small numbers. This interface was designed for tasks,
!           which are ordered after size (or expected cost), starting with
!           the large one. The interval is changed, as it provides also
!           a value for the stride.
!
!           call dlb_setup_rr(N) - once every time a DLB should
!              be used, N should be the number of jobs, belongs to
!              distribution with dlb_give_more_rr. It is illegal to
!              intermix the _rr routines with the other two interfaces.
!
!           dlb_give_more_rr(n, jobs) - the same as dlb_give_more, but
!              be aware that jobs has three elements. The jobslice it
!              provdides is jobs(0) +1: jobs(1): jobs(2), where jobs(2)
!              is the stride. jobs(0) +1 and jobs(2) are elements of the
!              interval. Be aware that for n = 1, jobs(2) is not the same
!              than jobs(1) + 1.
!
!  Module called by: ...
!
!
!  References:
!
!  Author: AN
!  Date: 10/10
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

!
! There are different implementations of DLB available, the interface
! looks all the same. The linking will define which one is used.
! See "use dlb_impl" instances in text, these are the references to
! the actual implementation.
!

# include "dlb.h"
! Need here some stuff, that is already defined elsewhere
use dlb_common, only: L_JOB
use dlb_common, only: lk, ik
use dlb_impl, only: DLB_THREAD_REQUIRED
use dlb_common, only: dlb_timers

implicit none
save            ! save all variables defined in this module
private         ! by default, all names are private
!== Interrupt end of public interface of module =================
 !------------ Declaration of types ------------------------------

public dlb_init, dlb_finalize, dlb_setup, dlb_setup_color, dlb_give_more
public dlb_give_more_color
public dlb_give_more_rr, dlb_setup_rr
public DLB_THREAD_REQUIRED
public dlb_print_statistics

! storage for the distribution of jobs over the colors, should hold exactly
! the distr one has given in dlb_setup_color, start_color is a helper variable
! for linking the intern job-numbers (succeeding ones) on the job-ids of the
! job distribution
integer, public, parameter :: idlb_kind = lk
integer (lk), allocatable :: start_color(:)

! in the color case one may not be able to hand over all jobs at once, thus store them here
integer (lk)              :: current_jobs(L_JOB)

! for the rr variant it is required to know the stride for the interval
integer (lk)              :: stride
contains

  subroutine dlb_init(world)
    !  Purpose: initialization of needed stuff
    !------------ Modules used ------------------- ---------------
    use dlb_impl, only: dlb_impl_init => dlb_init
    implicit none
    integer (ik), intent(in) :: world
    !** End of interface *****************************************

    call dlb_impl_init(world)
    current_jobs = 0
    stride = -1
  end subroutine dlb_init

  subroutine dlb_finalize()
    !  Purpose: cleaning up everything, after last call
    !------------ Modules used ------------------- ---------------
    use dlb_impl, only: dlb_impl_finalize => dlb_finalize
    use dlb_common, only: ik
    implicit none
    ! *** end of interface ***

    integer (ik)              :: ierr

    !
    ! Only if colored-version was in use:
    !
    if (allocated(start_color)) then
      deallocate(start_color, stat = ierr)
      ASSERT(ierr==0)
    endif

    call dlb_impl_finalize()
  end subroutine dlb_finalize

  subroutine dlb_setup(N)
    !  Purpose: initialization of a DLB run, each proc should call
    !           it with the number of jobs altogether. This is the
    !           version without color distinguishing, thus this information
    !           is enough
    !------------ Modules used ------------------- ---------------
    use dlb_common, only: distribute_jobs, n_procs, my_rank
    use dlb_impl, only: dlb_impl_setup => dlb_setup
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: N
    ! *** end of interface ***

    !
    ! FIXME: why my_rank as an argument?
    !        Should the input to dlb_impl_setup(...)
    !        depend on where it is executed?
    !        dlb_impl_setup gets on every process the tasks for this
    !        process.
    !
    call dlb_impl_setup(distribute_jobs(N, n_procs, my_rank))
  end subroutine dlb_setup

  subroutine dlb_setup_rr(N)
    !  Purpose: initialization of a DLB run, each proc should call
    !           it with the number of jobs altogether. This is the
    !           the version, where a round robin iteration is done
    !           over the jobs. Every process gets all tasks. It is
    !           also required to set the stride
    !------------ Modules used ------------------- ---------------
    use dlb_common, only: n_procs, my_rank, L_JOB, JRIGHT, JLEFT, JOWNER
    use dlb_impl, only: dlb_impl_setup => dlb_setup
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: N
    integer (lk)                :: jobs(L_JOB)
    ! *** end of interface ***
    jobs(JLEFT) = min(my_rank, N)
    jobs(JRIGHT) = N
    jobs(JOWNER) = my_rank
    stride = n_procs
    call dlb_impl_setup(jobs)
  end subroutine dlb_setup_rr

  subroutine dlb_setup_color(N)
    !  Purpose: initialization of a DLB run, each proc should call
    !           it with the number of jobs altogether. This is the
    !           version with color distinguishing, thus the distribution
    !           of the jobs over the color have to be given
    !------------ Modules used ------------------- ---------------
    use dlb_common, only: ik
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: N(:)
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer (ik)              :: ierr
    integer (lk)                :: i, number_jobs

    if (allocated(start_color)) then
      deallocate(start_color, stat = ierr)
      ASSERT(ierr==0)
    endif

    !
    ! There are size(N) many colors:
    !
    allocate(start_color(size(N) + 1), stat = ierr)
    ASSERT(ierr==0)

    ! start_color just keeps for every internal color the information
    ! with which of the internal job-numbers this color will
    ! start. The last entry point gives the number of all jobs
    ! altogether
    start_color(1) = 0
    do i = 2, size(N) + 1
      start_color(i) = start_color(i - 1) + N(i - 1)
    enddo

    !
    ! Internally jobs are treated as equal by DLB:
    !
    number_jobs = sum(N)
    call dlb_setup(number_jobs)
  end subroutine dlb_setup_color

  logical function dlb_give_more(n, my_job)
    !  Purpose: Returns next bunch of up to n jobs, if jobs(JRIGHT)<=
    !  jobs(JLEFT) there are no more jobs there, else returns the jobs
    !  done by the procs should be jobs(JLEFT) + 1 to jobs(JRIGHT) in
    !  the related job list
    !------------ Modules used ------------------- ---------------
    use dlb_impl, only: dlb_impl_give_more => dlb_give_more
    use dlb_common, only: JLEFT, JRIGHT, L_JOB
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: n
    integer (lk), intent(out  ) :: my_job(2)
    !** End of interface *****************************************
    integer (lk)                :: my_job_raw(L_JOB)
    ASSERT (stride < 0)

    call dlb_impl_give_more(n, my_job_raw)
    dlb_give_more = (my_job_raw(JLEFT) < my_job_raw(JRIGHT))
    my_job(1) = my_job_raw(JLEFT)
    my_job(2) = my_job_raw(JRIGHT)
  end function dlb_give_more

  logical function dlb_give_more_rr(n, my_job)
    !  Purpose: Returns next bunch of up to n jobs
    !           returns false if there are no more jobs left
    !           else it returns the job interval a+1:b:s as
    !           a = element left to the leftborder of interval
    !           b = right border of interval (inclusive)
    !           and s as stride.
    !
    !------------ Modules used ------------------- ---------------
    use dlb_impl, only: dlb_impl_give_more => dlb_give_more
    use dlb_common, only: JLEFT, JRIGHT, L_JOB, JOWNER
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: n
    integer (lk), intent(out  ) :: my_job(3)
    !** End of interface *****************************************
    integer (lk)                 :: my_job_raw(L_JOB), rest
    logical :: not_empty !, is_right
    ASSERT (stride > 0)
    not_empty = .false.

    ! cycle till the interval contains some tasks (or there is global termination).
    do while ( .not. not_empty ) ! second break if there are no more task around
        call dlb_impl_give_more(n * stride, my_job_raw)
        ! the global termination is reached if the interval is empty
        dlb_give_more_rr = (my_job_raw(JLEFT) < my_job_raw(JRIGHT))
        if (.not. dlb_give_more_rr) exit ! second break, no work left

        ! move the interval borders to values with the correct modulo (but still inside the interval).
        ! Consider: first JOWNER is 0, first task is 1, but we get left task -1 for historical reasons.
        rest = modulo(my_job_raw(JLEFT), stride) - my_job_raw(JOWNER)
        if (rest > 0 ) my_job_raw(JLEFT) = my_job_raw(JLEFT) + stride - rest
        if (rest < 0 ) my_job_raw(JLEFT) = my_job_raw(JLEFT) - rest
        ! Here we get the last valid right task. As those with modulo 1 should be mapped to
        ! owner number 0, we need to decrease the job ID here.
        rest = - modulo(my_job_raw(JRIGHT) - 1, stride) + my_job_raw(JOWNER)
        if (rest < 0) my_job_raw(JRIGHT) = my_job_raw(JRIGHT) + rest
        if (rest > 0 ) my_job_raw(JRIGHT) = my_job_raw(JRIGHT) - stride + rest
        not_empty = (my_job_raw(JLEFT) < my_job_raw(JRIGHT))

!       is_right = ( modulo(my_job_raw(JLEFT), stride) == my_job_raw(JOWNER))
!       ASSERT (is_right )
!       is_right = ( modulo(my_job_raw(JRIGHT) - 1, stride) == my_job_raw(JOWNER))
!       ASSERT (is_right )
    enddo
    my_job(1) = my_job_raw(JLEFT)
    my_job(2) =  my_job_raw(JRIGHT)
    my_job(3) = stride

    if (not_empty) then
        ASSERT (dlb_give_more_rr)
    endif
    ! next setup might be one without rr, thus reset stride to not get confused.
    if (.not. dlb_give_more_rr) stride = -1
  end function dlb_give_more_rr

  function dlb_give_more_color(n, color, my_job) result(more)
    !  Purpose: Returns next bunch of up to n jobs, if jobs(JRIGHT)<=
    !  jobs(JLEFT) there are no more jobs there, else returns the jobs
    !  done by the procs should be jobs(JLEFT) + 1 to jobs(JRIGHT) in
    !  the related job list
    !  this version gives also back the color of the jobs and ensures
    !  that the jobs given back all have the same color
    !  it keeps other colored jobs in its own storage
    !------------ Modules used ------------------- ---------------
    use dlb_common, only: ik
    use dlb_common, only: JLEFT, JRIGHT
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: n
    integer (lk), intent(out  ) :: my_job(2), color
    logical                              :: more ! result
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer (lk)                :: i,  w, jobs_all, jobs_color
    integer (ik)              :: ierr
    ASSERT (stride < 0)

    if (current_jobs(JLEFT) < current_jobs(JRIGHT)) then
        ! some are left over from the last time:
        more = .true.
    else
        ! only if the own storage is empty, refill:
        more = dlb_give_more(n, current_jobs)
    endif

    !
    ! ATTENTION: this if statement contains a return
    !
    if ( .not. more ) then
        ! got empty job from DLB, thus all done, quit
        my_job(1) = current_jobs(JLEFT)
        my_job(2) = current_jobs(JRIGHT)
        color = 0

        if (allocated(start_color)) then
            deallocate(start_color, stat = ierr)
            ASSERT(ierr==0)
        endif

        RETURN
    endif

    !
    ! Find out which color the first current job has:
    !
    do i = 2, size(start_color)
      if (start_color(i) > current_jobs(JLEFT)) then
        color = i - 1
        exit
      endif
    enddo

    ! will want
    jobs_all = current_jobs(JRIGHT) - current_jobs(JLEFT) ! how many jobs do I have
    jobs_color = start_color(color + 1) - current_jobs(JLEFT) ! how many jobs are there left
    ! of the color
    w = min(jobs_all, jobs_color)
    ! now share the own storage with the calling program
    my_job(1) = current_jobs(JLEFT) - start_color(color)
    my_job(2) = my_job(JLEFT) + w
    current_jobs(JLEFT) = current_jobs(JLEFT) + w
  end function dlb_give_more_color


  subroutine dlb_print_statistics(output_level)
    ! Purpose: Collects and prints some statistics of how the current DLB
    !          run has been performed
    !          For output_level = 0 exits without doing anything
    !          For output_level = 1 time spend in dlb_give_more (split)
    !          For output_level = 2 + wait for new tasks
    !          For output_level = 3 + last work slice
    !          For output_level = 4 + work length/task sizes + complete program
    !          ATTENTION: gives always the statistics of the last DLB run
    !                     if all statistics from several runs are wanted, has
    !                     to be called before the next dlb_setup
    !          The static backend provides only support for up to output_level 1
    !          Output for higher output level will be start values (0 or infinity)
    !------------ Modules used ------------------- ---------------
    use dlb_common, only: ik
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (ik), intent(in)  :: output_level
    !------------ Declaration of local variables ---------------------
    call dlb_timers(output_level)
  end subroutine dlb_print_statistics

end module dlb

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
module dlb_common
  !-------------------------------------------------------------------
  !
  !  Purpose: code common for all DLB implementations
  !
  !  Module called by: ...
  !
  !
  !  References:
  !
  !  Author: AN
  !  Date: 08/10->09/10
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
# include "dlb.h"
  USE_MPI!, only: MPI_STATUS_SIZE, MPI_SUCCESS, MPI_COMM_NULL, &
       !MPI_DOUBLE_PRECISION, MPI_WTIME, MPI_DATATYPE_NULL
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================
  !------------ public functions and subroutines ---------------------

  public :: distribute_jobs

  public :: steal_local!(m, local, remaining, stolen) result(ok)

  public :: length!(jobs) -> integer
  public :: empty!(jobs) -> logical

  public :: split_at! (C, AB, AC, CB)

  public :: time_stamp
  public :: dlb_common_init, dlb_common_finalize
  public :: set_start_job, set_empty_job

  public :: dlb_timers

  ! For those integers which will  stay 4 bytes no matter what changes
  ! around  (for example  for statistics).   One is  better  off using
  ! default  integers here as  this is  what MPI  likely uses  for its
  ! object handlers and counts. FIXME: maybe use kind(0) instead?
  integer, parameter, public :: ik = selected_int_kind (9)

  ! Integer with 8 bytes, or rather with range of 18 decimal digits:
  integer, parameter, public :: i8_kind = selected_int_kind(18)

  ! For job IDs, the actual data that DLB serves, use integers with at
  ! least that many decimal digits:
  integer, parameter, private :: kind_of_lk = 18
  integer, parameter, public :: lk = selected_int_kind (kind_of_lk)

  ! MPI has its own convention for  types, it needs its own version of
  ! the  integer  to  send.  This   will  be  set  to  something  more
  ! meaningfull in dlb_common_init():
  integer, public, protected :: lk_mpi = MPI_DATATYPE_NULL


  integer,  public, protected :: comm_world

  integer (ik), parameter, public  :: OUTPUT_BORDER = FPP_OUTPUT_BORDER

  integer (ik), parameter, public  :: DONE_JOB = 1, NO_WORK_LEFT = 2, RESP_DONE = 3 ! message tags
  integer (ik), parameter, public  :: WORK_REQUEST = 4, WORK_DONAT = 5 ! messages for work request
  integer (ik), parameter, public  :: JLENGTH = 3 ! Length of a single job in interface
  integer (ik), parameter, public  :: L_JOB = 3  ! Length of job to give back from inner interface (dlb_imp)
  integer (ik), parameter, public  :: JOWNER = 3 ! Number in job, where rank of origin proc is stored
  integer (ik), parameter, public  :: JLEFT = 1 ! Number in job, where stp (start point) is stored
  integer (ik), parameter, public  :: JRIGHT = 2 ! Number in job, where ep (end point) is stored

  ! Variables need for debug and efficiency testing
  ! They are the timers, needed in all variants, in multithreaded variants only
  ! MAIN is allowed to acces them
  double precision, public  :: main_wait_all, main_wait_max, main_wait_last
  double precision, public  :: max_work, last_work, average_work, num_jobs
  double precision, public  :: dlb_time, min_work, second_last_work
  double precision, public  :: timer_give_more, timer_give_more_last

  integer (ik), public, protected :: my_rank, n_procs ! some synonyms, They will be initialized
                                        !once and afterwards be read only

  !
  ! Termination master is the  process that gathers completion reports
  ! and  tells everyone to terminate.
  !
  integer (ik), public, protected :: termination_master

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !
  ! This scalar  holds the  number of jobs  initially assigned  to the
  ! process owning this variable.  FIXME: make it an array (1:n_procs)
  ! to keep  full info about initial assignment.  Currently not doable
  ! as  the  setup  procedure  is  provided only  one  assignment  (to
  ! myself).
  !
  integer (lk) :: responsibility = -1

  !
  ! This array holds  the counts of jobs that  were initially assigned
  ! to the process owning this structure and later reported as done by
  ! one  of  the workers.   Note  that  sum(reported_by(:))  is to  be
  ! compared  with with  the total  number  of jobs  from the  initial
  ! assignment kept in the variable "responsibility".
  !
  integer (lk), allocatable :: reported_by(:) ! (0:n_procs-1)

  !
  ! The next  array holds  the counts of  jobs that were  delivered to
  ! "userspace" of  the process owning  this structure (aka  jobs that
  ! were "done") and reported to the respective job owner according to
  ! initial assignment.   Note that sum(reported_to(:))  is the number
  ! of jobs "done" by the process owning this structure. We dont abuse
  ! message buffers to keep this data anymore.
  !
  integer (lk), allocatable :: reported_to(:) ! (0:n_procs-1)

  integer (ik), allocatable :: req_dj(:) ! need to store the
                                               ! messages for
                                               ! DONE_JOB,

  ! Need to store the messages for  DONE_JOB, as they may still be not
  ! finished, when the subroutine for generating them is finshed.
  integer (lk), allocatable :: messages(:,:)

  ! There may  be a lot of  them, message_on_way keeps  track, to whom
  ! they  are already  on their  way  the requests  are handeled  also
  ! separately, as it is needed to know, which one has finsished:
  logical, allocatable :: message_on_way(:) ! which messages are
                                            ! already sent

  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  double precision :: time_offset = -1.0
  logical, allocatable              :: all_done(:)
    ! only allocated on termination_master, stores which procs
    ! jobs are finished. If with masterserver: which proc has terminated
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

! ONLY FOR DEBBUGING (WITHOUT PARAGAUSS)
  subroutine time_stamp(msg, output_level)
    implicit none
    character(len=*), intent(in) :: msg
    integer (ik), intent(in) :: output_level
    ! *** end of interface ***

    double precision :: time
    character(len=28) :: prefix
    ! ATTENTION: write into character and printing might be not thread save
    ! for all compilers. Thus keep write and print in the if-clause allowing
    ! to get rid of them by setting OUTPUT_BORDER=0.

    time = MPI_Wtime()
    if ( time_offset < 0.0 ) then
      time_offset = time
      if(output_level < OUTPUT_BORDER) then
         write(prefix, '("#", I3, G20.10)') my_rank, time_offset
         print *, prefix, "(BASE TIME)"
      endif
    endif

    if(output_level < OUTPUT_BORDER) then
      write(prefix, '("#", I3, G20.10)') my_rank, time - time_offset

      print *, prefix, msg
    endif
  end subroutine time_stamp
! END ONLY FOR DEBUGGING

  subroutine dlb_timers(output_level)
    ! Purpose: Collects and prints some statistics of how the current dlb
    !          run has been performed
    !          For output_level = 0 exits without doing anything
    !          For output_level = 1 time spend in dlb_give_more (split)
    !          For output_level = 2 + wait for new tasks
    !          For output_level = 3 + last work slice
    !          For output_level = 4 + work length/task sizes + complete program
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (ik), intent(in)  :: output_level
    !------------ Declaration of local variables ---------------------
    double precision, allocatable :: times(:,:), help_arr(:)
    double precision              :: time_singles(12)
    integer (ik)              :: ierr
    ! *** end of interface ***
    ! Return if none output is wanted, then there is also no need to
    ! send the output to the master
    if (output_level == 0) RETURN

    ! Master needs some place to store all results
    if (my_rank == 0) then
        allocate(times(12, n_procs), help_arr(n_procs), stat = ierr )
        ASSERT (ierr == 0)
        times = 0
    else
        ! will not contain any information, but better have it also on the
        ! slaves
        allocate(times(1, 1), stat = ierr )
        ASSERT (ierr == 0)
        times = 0
    endif

    time_singles(1) = timer_give_more
    time_singles(2) = timer_give_more_last
    time_singles(3) = main_wait_all
    time_singles(4) = main_wait_max
    time_singles(5) = main_wait_last
    time_singles(6) = dlb_time
    time_singles(7) = max_work
    time_singles(8) = min_work
    time_singles(9) = average_work
    time_singles(10) = last_work
    time_singles(11) = second_last_work
    time_singles(12) = num_jobs
    ASSERT (size(time_singles) == 12)

    call MPI_GATHER(time_singles, size(time_singles), MPI_DOUBLE_PRECISION,  &
          times, size(time_singles), MPI_DOUBLE_PRECISION, 0, comm_world, ierr)
    ASSERT (ierr == 0)

    if (my_rank == 0) then
        print *, "--DLB-statistics--DLB-statistics--DLB-statistics--DLB-statistics--"
        print *, "Statistics for last DLB loop (times in seconds):"
        print *, "1. Time spend in DLB"
        print *, "  DLB loop time (without terminating) = ", sum(times(1,:))
        print *, "  DLB termination                     = ", minval(times(2,:)) * n_procs
        print *, "  DLB last call - termination         = ", sum(times(2,:)) - n_procs * minval(times(2,:))
        if (output_level > 1) then
            print *, "2. Waiting times for new task"
            print *, "  Maximum wait  min/max               =", minval(times(4,:)), maxval(times(4,:))
            print *, "  Total wait                          =", sum(times(3,:))
            print *, "  Total wait times without last       =", sum(times(3,:)) - sum(times(5,:))
            print *, "  Number of tasks on proc.  min/max   =", minval(times(12,:)) , maxval(times(12,:))
        endif
        if (output_level > 2) then
           print *, "3. Times of last work slices (time between 2 dlb_give_more calls)"
           print *, "  Last work slice min/max             =", minval(times(10,:)) , maxval(times(10,:))
           print *, "  Second last work slice min/max      =", minval(times(11,:)) , maxval(times(11,:))
        endif
        if (output_level > 3) then
           print *, "  Maximum task on proc min/max        =", minval(times(7,:)) , maxval(times(7,:))
           print *, "  Minimal task on proc min/max        =", minval(times(8,:)) , maxval(times(8,:))
           help_arr = times(9,:) / times(12,:)
           print *, "4. Average work load (time between the dlb calls / Number of tasks)"
           print *, "  Average work load all/min/max       =", sum(times(9,:)) / sum(times(12,:)), minval(help_arr) &
                                                           , maxval(help_arr)
           print *, "5. Complete Time (between dlb_setup and termination)"
           print *, "  Complete Time all                   =", sum(times(6,:))
           print *, "  Time for single proc min/max        =", minval(times(6,:)), maxval(times(6,:))
        endif
        print *, "----DLB-statistics-end--DLB-statistics-end--DLB-statistics-end----"
        print *, ""
        if (allocated(times)) then
            deallocate(times, stat = ierr)
            ASSERT (ierr == 0)
        endif
    endif
  end subroutine dlb_timers


  subroutine dlb_common_init(world)
    ! Intialization of common stuff, needed by all routines
    implicit none
    integer (ik) ,  intent(in) :: world
    ! *** end of interface ***

    integer (ik) :: ierr, alloc_stat

    !
    ! Set global communicator as a DUP of the world:
    !
    call MPI_COMM_DUP(world, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    call MPI_COMM_RANK(comm_world, my_rank, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    call MPI_COMM_SIZE(comm_world, n_procs, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    !
    ! Create a  new MPI_Datatype, here  lk_mpi. It is not  a kind
    ! (like  lk), it  is not  a  number of  decimal digits  (like
    ! kind_of_lk) but a real  MPI type, represented as an integer
    ! as everything else in Fortran bindings.
    !
    ! FIXME: even with the guard, valgrind reports a memory leak here.
    ! But at least one MPI interpretation says:
    !
    !   "An  MPI_Datatype  returned  by  this  subroutine  is  already
    !    committed.  It cannot be freed with MPI_TYPE_FREE()."
    !
    ! Here kind_of_lk is a  PARAMETER and does not change between
    ! invocations. So do it just once:
    !
    if (lk_mpi == MPI_DATATYPE_NULL) then
       call MPI_TYPE_CREATE_F90_INTEGER (kind_of_lk, lk_mpi, ierr)
       ASSERT(ierr==MPI_SUCCESS)
    endif

    !
    ! Anyone  can  be  the  termination  master. So  choose  the  last
    ! processor  to  be  the  termination  master.   It  should  exist
    ! irrespective  what  the  number  of  processors  is.  The  first
    ! processor is most likely playing the role of master for the real
    ! calculation,  this  way the  additional  load  should be  better
    ! balanced.
    !
    termination_master = n_procs - 1

    !
    ! Only on termination master:
    !
    if (my_rank == termination_master) then
      allocate(all_done(n_procs), stat = alloc_stat)
      ASSERT(alloc_stat==0)
    endif
  end subroutine dlb_common_init

  subroutine dlb_common_finalize()
    ! shut down the common stuff, as needed
    implicit none
    ! *** end of interface ***

    integer (ik) :: ierr, alloc_stat

    !
    ! Only on termination master:
    !
    if (allocated(all_done)) then
      deallocate(all_done, stat=alloc_stat)
      ASSERT(alloc_stat==0)
    endif

    call MPI_COMM_FREE(comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
    ASSERT(comm_world==MPI_COMM_NULL)
  end subroutine dlb_common_finalize

  ! To make it easier to add for example new job informations
  ! The direct setting of a comlete job is always done here
  ! There are two functions needed: setting of the inital job
  ! for the own rank and setting an empty, easy as that to
  ! recognise job
  ! it is supposed, that all the informations needed so far, will be
  ! still there, but that new informations will be handled on too
  ! (without this two functions, the informations will just copied
  ! and special ones will be changed)

  function set_start_job(job) result(job_data)
    !Purpose: gives a complete starting job description
    !         after gotten just the job range
    implicit none
    integer (lk), intent(in) :: job(L_JOB)
    integer (lk)             :: job_data(JLENGTH)
    ! *** end of interface ***

    job_data(:L_JOB) = job
    ASSERT (job_data(JOWNER) == my_rank)
  end function set_start_job

  function set_empty_job() result(job_data)
    !Purpose: gives a complete starting job description
    !         after gotten just the job range
    implicit none
    integer (lk) :: job_data(JLENGTH)
    ! *** end of interface ***

    job_data(JRIGHT) = 0
    job_data(JLEFT) = 0
    job_data(JOWNER) = -1
  end function set_empty_job

  function steal_local(m, local, remaining, stolen) result(ok)
    !
    ! "Stealing" from myself is implemented as splitting the interval
    !
    !     (A, B] = local(1:2)
    !
    ! into
    !
    !     stolen(1:2) = (A, C] and remaining(1:2) = (C, B]
    !
    ! at the point
    !
    !     C = A + min(B - A, m)
    !
    ! NOTE: This function has to adhere to the interface of modify(...)
    ! argument in try_read_modify_write(...)
    !
    ! use dlb_common, only: lk, JLENGTH, split_at
    implicit none
    integer (lk), intent(in)  :: m
    integer (lk), intent(in)  :: local(:) ! (JLENGTH)
    integer (lk), intent(out) :: remaining(:) ! (JLENGTH)
    integer (lk), intent(out) :: stolen(:) ! (JLENGTH)
    logical                       :: ok ! result
    ! *** end of interface ***

    integer (lk) :: work, c

    ASSERT(size(local)==JLENGTH)
    ASSERT(size(remaining)==JLENGTH)
    ASSERT(size(stolen)==JLENGTH)
    ASSERT(m>=0)

    ! The give_grid function needs only up to m jobs at once:
    work = min(length(local), m)
    ! work = max(work, 0)
    ASSERT(work>=0)

    ! split here:
    c = local(JLEFT) + work

    !
    ! Note the order of (stolen, remaining) --- left interval is stolen, right
    ! interval is remaining:
    !
    call split_at(c, local, stolen, remaining)

    ok = .not. empty(stolen)
  end function steal_local

  logical function empty(jobs)
    implicit none
    integer (lk), intent(in) :: jobs(:)
    ! *** end of interface ***

    empty = length(jobs) == 0
  end function empty

  function length(jobs) result(n)
    implicit none
    integer (lk), intent(in) :: jobs(:)
    integer (lk)             :: n ! result
    ! *** end of interface ***

    ASSERT(size(jobs)==JLENGTH)

    n = max(jobs(JRIGHT) - jobs(JLEFT), 0)
  end function length

  subroutine split_at(C, AB, AC, CB)
    !
    ! Split (A, B] into (A, C] and (C, B]
    !
    implicit none
    integer (lk), intent(in)  :: C, AB(:)
    integer (lk), intent(out) :: AC(:), CB(:)
    ! *** end of interface ***

    ASSERT(size(AB)==JLENGTH)
    ASSERT(size(AC)==JLENGTH)
    ASSERT(size(CB)==JLENGTH)

    ASSERT(C>=AB(JLEFT))
    ASSERT(C<=AB(JRIGHT))

    ! copy trailing posiitons, if any:
    AC(:) = AB(:)
    CB(:) = AB(:)

    AC(JLEFT)  = AB(JLEFT)
    AC(JRIGHT) = C

    CB(JLEFT)  = C
    CB(JRIGHT) = AB(JRIGHT)
  end subroutine split_at

  function distribute_jobs(N, procs, rank) result(my_jobs)
    !  Purpose: given the number of jobs alltogether, decides how many
    !           will be done on each proc and where is start and endpoint
    !           in the global job range, this will be fed directly in
    !           the dlb_setup of the dlb routine
    !
    !           The last processes will get fewer jobs, if it is not
    !           dividable equally.
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer (lk), intent(in   ) :: N
    integer (ik), intent(in ) :: procs, rank
    integer (lk)                :: my_jobs(L_JOB)
    !** End of interface *****************************************

    !------------ Declaration of local variables ---------------------
    integer (lk) :: n_procs, my_rank
    integer (lk) :: jobs_per_proc
    integer (lk) :: rest
    n_procs = procs
    my_rank = rank

    ! jobs_per_proc = minimum number of jobs per processes
    ! first rest jobs will get one more job
    jobs_per_proc = N / n_procs
    ! N = jobs_per_proc * n_procs + rest; rest < n_procs
    rest = N - jobs_per_proc * n_procs

    if (my_rank == termination_master .and. 1 < OUTPUT_BORDER) then
        print *, "Distributing of jobs on the processors"
        print *, "There are ", N, " jobs altogether"
        print *, "each processor will get approximatly", jobs_per_proc
    endif

    ! Starting point of own interval:
    my_jobs(JLEFT) = jobs_per_proc * my_rank + min(rest, my_rank)

    ! if this proc has to do one more job tell him so
    ! changes job_per_proc from minimum value
    !to jobs_per_proc for me
    if (my_rank < rest) jobs_per_proc = jobs_per_proc + 1

    my_jobs(JRIGHT) = my_jobs(JLEFT) + jobs_per_proc
    my_jobs(JOWNER) = my_rank
  end function distribute_jobs

end module dlb_common

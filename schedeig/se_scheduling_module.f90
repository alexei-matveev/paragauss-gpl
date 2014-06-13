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
module se_scheduling_module
  !---------------------------------------------------------------
  !
  !  Purpose: Generates a scheduling for the diagonalization of
  !           a set of matrices on a set of processors. The runtime
  !           of the eigensolvers employed for the diagonalization
  !           is modelled by a cost function.
  !
  !
  !  Module called by: se_eigen_module.f90
  !
  !
  !  References: W. Ludwig and P. Tiwari: Scheduling malleable and
  !              nonmalleable parallel tasks. SODA ’94: Proceedings
  !              of the fifth annual ACM-SIAM symposium on Discrete
  !              algorithms, 167–176 (1994)
  !
  !
  !
  !  Author: Martin Roderus
  !  Date: 4/2010
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

# include "def.h"
  ! FIXME:  original  code  use  single  precision  for  real  numbers
  ! (r4_kind), but  comm primitives are not implemented  for them. Add
  ! them to comm.f90 or use MPI directly.
  use type_module, only: IP => i4_kind, WP => r8_kind ! kind params
  use se_timefunction_module, only: se_timefunction_type
  implicit none
  private


  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !
  ! The  optimal   scheduling  is  stored  in   the  following  nested
  ! structure.  The  eigensolver schedules the tasks  as defined here,
  ! so it accesses the structure in read-only fashion.
  !
  ! Times  are not  necessary for  the eigensolver.   Required  is the
  ! field  %tasks.  In  the struct  array %tasks,  the  only necessary
  ! fields  are %allotedProcessors,  %id.   The field  %size is  used,
  ! among other things to verify if the schedule can be re-used in the
  ! next iteration in case the problem sizes did not change.
  !
  ! The local constraints on redundant data include
  !
  !     size(%allotedProcessors) == %slots
  !
  type, public :: taskType
     integer(IP), allocatable, public :: allotedProcessors(:)

     !
     ! Each  entry of  the array  of  taskType structs  refers to  the
     ! original task.  These are  their IDs (sequence numbers) and the
     ! corresponding  problem size. FIXME:  %size duplicates  the info
     ! that can  be derived  by the caller  from the task  ID. Problem
     ! sizes were the input to this module, obviousely:
     !
     integer(IP), public :: id = -1
     integer(IP), public :: size = -1
  end type taskType

  !------------ public functions and subroutines ------------------
  public :: se_scheduling_run

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of interface blocks --------------------
  interface sort
    module procedure sort_r, sort_i
  end interface sort

  type :: npts_task
     !
     ! We will be "permuting" the tasks. These are their IDs (sequence
     ! numbers), number of workers for  the task to be executed on and
     ! an (approximate) time span:
     !
     integer(IP) :: id = -1
     integer(IP) :: slots = -1 ! = size(%allotedProcessors)
     real(WP) :: span = -1. ! = getExecTime(%slots, %size)
  end type npts_task

  !------------ Declaration of constants and variables ------------
  type(se_timefunction_type), allocatable :: timeFunctions(:)
  real(WP), parameter :: STD_TIME = 0.001_WP


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine se_scheduling_run(tasks, task_sizes, np)
    !
    ! Purpose:  public  interface.  Performs  the  necessary steps  to
    ! generate  a  scheduling and  stores  it  into  array of  structs
    ! tasks(:). FIXME: intent(out) arg in the first position?
    !
    use comm, only: comm_rank
    implicit none
    type(taskType), intent(out), allocatable :: tasks(:) !(size(task_sizes)), the (reordered) schedule info
    integer(IP), intent(in) :: task_sizes(:) ! Contains information about the task sizes
    integer(IP), intent(in) :: np ! Indicates how many processors are available for the scheduling
    ! *** end of interface ***

    type(npts_task) :: npts_tasks(size(task_sizes))
    integer :: i

    !
    ! Find  for each  matrix a  number of  processors on  which  it is
    ! processed.  This converts an  MPTS problem  to the  NPTS problem
    ! that is next solved by a combinatoric approach:
    !
    call findOmega(task_sizes, np, npts_tasks)

    ASSERT(all(npts_tasks%id>0))

    !
    ! Generate scheduling.  Reorderes npts_tasks so that makeSpan() is
    ! minimized.  Check their %id to get the new sequence:
    !
    call combinatoricSchedule(npts_tasks, np)

    !
    ! Finally, write  the scheduling data  into the tasks(:)  array of
    ! structs understood outside of this module.
    !
    ! FIXME: The  time estimates invoked  in recCombinatoricSchedule()
    ! and  the  detailed   final  schedule  representation  should  be
    ! consistently  generated.  Currently   the  logic  is  duplicated
    ! between heavily used makeSpan(), O(N!) uses, and makeSchedule(),
    ! only one use here:
    !
    call makeSchedule(npts_tasks, np, tasks)

    !
    ! FIXME: code outside  relies on task sizes stored  here too. Note
    ! that  combinatoricSchedule() returns  NPTS-tasks  in a  permuted
    ! sequence  such  that  a  makeSchedule() generates  an  "optimal"
    ! schedule:
    !
    tasks%size = task_sizes(tasks%id)

    if (comm_rank() == 0) then
       do i = 1, size(tasks)
          print *, "SCHEDULE: taks(", tasks(i)%id, ") ", &
               "size=", tasks(i)%size, &
               "span=", npts_tasks(i)%span, &
               "cpus=", tasks(i)%allotedProcessors
       enddo
    endif
  end subroutine se_scheduling_run

#ifdef WITH_GUILE
  subroutine findOmega(task_sizes, np, tasks)
    !
    ! Convert an MPTS problem to NPTS problem by deciding for a number
    ! of workers a mallable task  should be executed on, and providing
    ! an estimate for the time  required.  Note that sum of all worker
    ! counts does not necessarily yield total number of workers, NP.
    !
    use scm, only: scm_t, scm_from, scm_call, scm_car, scm_cdr, &
         assignment(=)
    use scheme, only: scheme_make_list, scheme_lookup
    implicit none
    integer(IP), intent(in) :: task_sizes(:)
    integer, intent(in) :: np
    type(npts_task), intent(out) :: tasks(:)
    ! *** end of interface ***

    type(scm_t) :: proc, pairs, pair
    integer :: i, slots
    real(WP) :: span

    !
    ! We delegate the task to  convert an MPTS problem to NPTS problem
    ! to  the Scheme  procedure.  It  is  an error  when the  function
    ! cannot be found:
    !
    proc = scheme_lookup ("qm-mpts->npts", mod="guile scheduling")

    !
    ! Return value is a list of pairs (slots . span) for each task:
    !
    pairs = scm_call (proc, scheme_make_list (task_sizes), scm_from (np))

    do i = 1, size (task_sizes)
       pair = scm_car (pairs)
       slots = scm_car (pair)   ! integer
       span = scm_cdr (pair)    ! double

       ! Type constructor here:
       tasks (i) = npts_task (i, slots, span)

       ! Advance the head "pointer":
       pairs = scm_cdr (pairs)
    enddo
  end subroutine findOmega
#else
  subroutine findOmega(task_sizes, np, tasks)
    !
    ! Purpose: find  the smallest time  value which yields  a feasible
    ! scheduling.
    !
    use se_timefunction_module, only: se_timefunction_init
    implicit none
    integer(IP), intent(in) :: task_sizes(:)
    integer, intent(in) :: np
    type(npts_task), intent(out) :: tasks(:)
    !** End of interface *****************************************

    real(WP)              :: meanWork, meanWork_min, omega, omega_min, tau, tau_min
    real(WP), allocatable :: timeMatrix(:)
    integer(IP)           :: i, j, mTasks, maxproc_index

    !
    ! Fetch the coeffs of the cost function:
    !
    call se_timefunction_init(timeFunctions)

    mTasks = size(task_sizes)

    !
    ! Find index for costfunction for processor count <= max available processors,
    !
    maxproc_index = -1
    do i = 1, size(timeFunctions)
       if( timeFunctions(i)%procCount > np ) exit
       maxproc_index = i
    end do

    !
    ! If there is none, return early. FIXME: does the "short-circuit" branch
    ! allocate/set all of the "output" fields?
    !
    if( maxproc_index < 1 ) then
       ! Assign standard values and return
       return ! *** RETURN POINT ***
    end if


    ! Matrix with execution times of each task with each possible
    ! processor count. Tasks are in rows, processor count in columns
    allocate(  timeMatrix( mTasks*maxproc_index )  )


    do i=1, mTasks
       do j=1, maxproc_index
          timeMatrix( i + (j-1)*mTasks ) = getExecTime(timeFunctions(j)%procCount, task_sizes(i))
       end do
    end do

    call sort(timeMatrix)


    omega_min = timeMatrix(size(timeMatrix))        ! The highest possible value for omega
    tau_min = omega_min
    meanWork_min = 0.
    do i=1, mTasks*maxproc_index
       if( timeMatrix(i) .gt. 0. ) then
          tau = timeMatrix(i)
          if (findProcessorMap(task_sizes, np, tau, tasks)) then

             ! Compute mean work of allotment
             meanWork=0.
             do j=1, mTasks
                meanWork = meanWork + tasks(j)%span * tasks(j)%slots
             end do
             meanWork = meanWork / np

             ! Write tau, mean work and omega into the scheduling type
             ! Get omega = max( W(tau), tau )
             omega = max(tau, meanWork)
             if( omega .lt. omega_min ) then
                omega_min = omega
                tau_min = tau
                meanWork_min = meanWork
             end if
          end if
       end if
    end do

    if (.not. findProcessorMap(task_sizes, np, tau_min, tasks)) then
       ! print fallback warning
    end if

    deallocate( timeMatrix )

    !
    ! Cleaning up internally allocated  memory --- deallocate a nested
    ! struct:
    !
    deallocate(timeFunctions)
  end subroutine findOmega
  !*************************************************************


  !*************************************************************
  logical function findProcessorMap(task_sizes, np, tau, tasks)
    !  Purpose: ..
    ! Called by: findOmega
    implicit none
    integer(IP), intent(in) :: task_sizes(:)
    integer, intent(in) :: np
    real(WP), intent(in) :: tau ! time threshold
    type(npts_task), intent(out) :: tasks(:)
    !** End of interface *****************************************

    integer(IP)          :: i, j, taskCount, maxproc_index
    real(WP)             :: timeDummy

    findProcessorMap = .true.

    taskCount = size(task_sizes)

    !
    ! Find index for costfunction for processor count <= max available processors,
    !
    maxproc_index = -1
    do i = 1, size(timeFunctions)
       if( timeFunctions(i)%procCount > np ) exit
       maxproc_index = i
    end do

    !
    ! If there is none, return early. FIXME: does the "short-circuit" branch
    ! allocate/set all of the "output" fields?
    !
    if( maxproc_index < 1 ) then
       ! Generate fallback solution and return
       do i = 1, size(task_sizes)
          tasks(i) = npts_task(i, 1, STD_TIME)
       end do
       return ! *** RETURN POINT ***
    end if

    do i = 1, taskCount
       do j=1, maxproc_index
          timeDummy = getExecTime(timeFunctions(j)%procCount, task_sizes(i))
          if( timeDummy .le. tau  .and.  timeDummy .ge. 0 ) then
             tasks(i) = npts_task(i, timeFunctions(j)%procCount, timeDummy)
             exit
          else if( timeDummy .lt. 0) then
             ! Time function did not yield a sensible time. Allot one processor and a standard time constant.
             tasks(i) = npts_task(i, 1, STD_TIME)
             findProcessorMap = .false.
             exit
          end if
       end do

       if (j > maxproc_index) then
          !
          ! No time  below "tau" was found. Allot  P processors.  This
          ! is not a dead code, it is enetered in case the loop over j
          ! = 1, maxproc_index runs through.
          !
          tasks(i) = npts_task(i, &
               timeFunctions(maxproc_index)%procCount, &
               getExecTime(tasks(i)%slots, task_sizes(i)))
          findProcessorMap = .false.
       end if
    end do
  end function findProcessorMap
  !*************************************************************


  !*************************************************************
  real(WP) function getExecTime(np, task_size)
    !
    ! Purpose: returns a  predicted execution time of a  task on given
    ! 'np' processors. The cost function which is imployed to make the
    ! prediction is a polynomial function; its coefficients are stored
    ! in 'timeFunctions'.
    !
    implicit none
    integer(IP), intent(in) :: np
    integer(IP), intent(in) :: task_size
    !** End of interface *****************************************

    integer(IP) :: i, polyDegree, minS, iProc = -1
    real(WP) :: scaleFactor

    do i=1, size(timeFunctions)
       if( timeFunctions(i)%procCount .eq. np ) then
          iProc = i
          exit
       end if
    end do

    if( iProc .eq. -1 ) then
       getExecTime = -1.
       return
    end if

    polyDegree = size(timeFunctions(iProc)%coefficients, 1)
    scaleFactor = timeFunctions(iProc)%scaleFactor

    minS = timeFunctions(iProc)%minSize
    if( task_size .lt. timeFunctions(iProc)%minSize ) then
       getExecTime = -1_WP
       return
    end if

    getExecTime = 0.
    do i=1, polyDegree
       getExecTime = getExecTime + timeFunctions(iProc)%coefficients(i) * (task_size/scaleFactor)**(i-1)
    end do

    if( getExecTime .le. 0 ) then
       getExecTime = -1_WP
    end if
  end function getExecTime
#endif /* ifdef WITH_GUILE */


  !*************************************************************
  subroutine combinatoricSchedule(tasks, np)
    !
    ! Purpose:  initialization for  recursive combinatoric  routine to
    ! find the scheduling with the best makespan 'minExecTime'.
    !
    ! Most accesses  to the tasks(:) array  in se_eigen_module occured
    ! in this form: i  -> tasks(task_sequence(i)).  So for convenience
    ! we order the tasks in their "processing" sequence instead of the
    ! original  one   here  to  address  them  more   directly:  i  ->
    ! tasks(i). The original task IDs (sequence indices as provided by
    ! the caller) are stored in tasks(:)%id.
    !
    use comm, only: comm_size, comm_rank, comm_minloc, comm_bcast
    implicit none
    type(npts_task), intent(inout) :: tasks(:)
    integer, intent(in) :: np
    !** End of interface *****************************************

    type(npts_task) :: permute(size(tasks))
    integer(IP) :: i, j, k, rank, counter, loc(1)
    real(WP) :: minExecTime

    ASSERT(np==comm_size())
    rank = comm_rank()
    counter = 0

    !
    ! Fill local  array "permute" with an  initial scheduling sequence
    ! of tasks:
    !
    j = 1
    k = size(tasks)
    do i = 1, size(tasks)
       if (tasks(i)%slots == np) then
          ! If the task has alloted all available processors, schedule
          ! it at the beginning:
          permute(j) = tasks(i)
          j = j + 1
       else
          ! Otherwise at the end:
          permute(k) = tasks(i)
          k = k - 1
       end if
    end do

    !
    ! These    two   arguments   to    recCombinatoricSchedule()   are
    ! intent(inout):
    !
    minExecTime = huge(minExecTime) ! some very high number
    tasks = permute

    !
    ! Start the  recursive combinatoric procedure. Only  the tasks are
    ! considered  which have allotted  less processors  than available
    ! (j:size(tasks)). The first array is a working array, its entries
    ! are  permuted  in  all  possible combinations  and  tested  with
    ! makeSpan().   The last  argument  is the  "optimal" sequence  on
    ! output  (or unchanged  if no  better  was found).   I think  the
    ! constrain is  that initial values of "permute"  and "tasks" need
    ! to be equal:
    !
    call recCombinatoricSchedule(permute(j:), 1, np, rank, counter, &
         minExecTime, tasks(j:))

    !
    ! Each worker processed a subset of all permutations, they need to
    ! agree on the best result  now. Here loc(1) will contain the rank
    ! of the worker holding (one of) the minimal value of minExecTime:
    !
    loc = comm_minloc(minExecTime)

    !
    ! FIXME: non-continous arrays here, will cause copy in/out on call
    ! to MPI_BCAST:
    !
    call comm_bcast(tasks(j:)%id, root=loc(1))
    call comm_bcast(tasks(j:)%slots, root=loc(1))
    call comm_bcast(tasks(j:)%span, root=loc(1))
  end subroutine combinatoricSchedule
  !*************************************************************


  !*************************************************************
  recursive subroutine recCombinatoricSchedule(tasks, i, np, rank, counter, minExecTime, optimal)
    !
    ! Purpose:  Recursive  routine to  find  the  scheduling with  the
    ! minimum makespan  'minExecTime'.  A scheduling is  made from all
    ! possible permutations of 'tasks'. The permutations are generated
    ! recursively.
    !
    ! The procedure searches only a subset of the permutation sequence
    ! such that rank = modulo(counter, np).
    !
    ! FIXME: at least  of O(N!) complexity with N  = size(tasks). Call
    ! stack of O(N) depth.
    !
    ! FIXME: "np" is both, the size of the "virtual" machine for which
    ! the schedule is optimized AND the  size of the MPI wolrd used to
    ! share the work.
    !
    implicit none

    type(npts_task), intent(inout) :: tasks(:)
    integer(IP), intent(in) :: i ! only tasks(i:) to be permuted
    integer(IP), intent(in) :: np, rank
    integer(IP), intent(inout) :: counter ! modulo np
    type(npts_task), intent(inout) :: optimal(:)
    real(WP), intent(inout) :: minExecTime
    ! *** end of interface ***

    integer(IP) :: j
    real(WP) :: execTime

    if (i < size(tasks)) then
       !
       ! "Branch node" of the  recursion.  Initiates size(tasks) - i +
       ! 1  recursive calls:
       !
       do j = i, size(tasks)
          call swap(tasks(i), tasks(j))
          call recCombinatoricSchedule(tasks, i + 1, np, rank, counter, minExecTime, optimal)
          call swap(tasks(i), tasks(j))
       end do
    else
       !
       ! "Leaf node"  of the  recursion over permutations.   The first
       ! part is the poor man  work sharing construct --- skip some of
       ! the permutations if the remainder does not match my rank. The
       ! caller will  need to reduce  the result from all  ranks (Note
       ! that if  called with  initial value of  counter = 0  then the
       ! first permutation will be processed  by rank = 1, not by rank
       ! = 0):
       !
       counter = modulo(counter + 1, np)
       if (counter /= rank) then
          ! skip computation  of the makespan and proceed  to the next
          ! permutation (return):
          return ! *** RETURN POINT ***
       endif

       execTime = makeSpan(tasks, np)

       if (execTime < minExecTime) then
          !
          ! Found a  sequence for the  tasks that, when  combined with
          ! the  reservation procedure coded  in makeSpan/makeSchedule
          ! results in shorter time estimate:
          !
          minExecTime = execTime
          optimal = tasks
       end if
    end if

    contains
      subroutine swap(v1, v2)
        !  Purpose: helper function to swap two elements
        implicit none
        type(npts_task), intent(inout) :: v1, v2
        !** End of interface *****************************************

        type(npts_task) :: v_temp

        v_temp = v1
        v1 = v2
        v2 = v_temp
      end subroutine swap
  end subroutine recCombinatoricSchedule


  subroutine makeSchedule(npts_tasks, np, tasks)
    !
    ! Purpose: Generates a scheduling for given tasks on all available
    ! processors.   Implemented is  a LIST-scheduling  according  to a
    ! certain sequence  of tasks as  given. Each task is  scheduled as
    ! soon as the required number of processors is available.
    !
    ! FIXME: O(N*P^2) complexity with N  = size(tasks) and P being the
    ! typical number of processors per task.
    !
    ! FIXME: the logic has to be identical to that in makeSpan().
    !
    implicit none
    type(npts_task), intent(in) :: npts_tasks(:)
    integer(IP), intent(in) :: np
    type(taskType), allocatable, intent(out) :: tasks(:)
    ! *** end of interface ***

    integer(IP) :: i, k, allotedProcCount

    ! Array which stores the blocked times for each processor. Compare
    ! with the array times(:) in  makeSpan(). Here each bank is exatly
    ! one  slot, because  we  are  supposed to  return  an info  about
    ! identities of the reserved processors:
    real(WP) :: procTimeArray(np)
    real(WP) :: time, execTime

    ! Contains the processor numbers which are available
    ! first. Resulting array for MINLOC-function:
    integer(IP) :: chosenProcessors(np)

    allocate(tasks(size(npts_tasks)))

    !
    ! Fill in the  sizes, code of this module  prefers to operate with
    ! taskType structs:
    !
    tasks%id = npts_tasks%id
    tasks%size = -1

    ! Initialization
    procTimeArray = 0.
    execTime = 0.

    do i = 1, size(npts_tasks)
       ! Copy number of alloted processors for relevant task to local variable
       allotedProcCount = npts_tasks(i)%slots

       !
       ! Get  the indices  of  processors which  are available  first,
       ! ordered by the time they  are available. Think of (take (sort
       ! list) n).  Thus,  in the last position we  will find the time
       ! point  for the  "task(i)" to  have all  of "allotedProcCount"
       ! processors available to start the execution of the task:
       !
       call minlocs(procTimeArray, chosenProcessors(:allotedProcCount))

       !
       ! At this moment the "task(i)" will start execution:
       !
       time = procTimeArray(chosenProcessors(allotedProcCount))

       !
       ! At this moment the "task(i)" will be completed:
       !
       time = time + npts_tasks(i)%span

       !
       ! Update 'procTimeArray'. Up  to this point in time  all of the
       ! chosen processrs will be considered blocked:
       !
       do k = 1, allotedProcCount
          procTimeArray(chosenProcessors(k)) = time
       end do

       !
       ! Check if the new scheduled task exceeds the makespan. If yes,
       ! update it:
       !
       if (time > execTime) execTime = time

       allocate(tasks(i)%allotedProcessors(allotedProcCount))
       tasks(i)%allotedProcessors = chosenProcessors(:allotedProcCount)
       call sort(tasks(i)%allotedProcessors)
    end do

    ! FIXME: this is consistency checking:
    execTime = makeSpan(npts_tasks, np) - execTime
    ASSERT(execTime==0.0)

  contains
    subroutine minlocs(a, locs)
      !
      ! Find size(locs) smallest values  in array "a" and return their
      ! indices.
      !
      ! O(N*M) version  with N  = size(a) and  M =  size(locs). FIXME:
      ! does an optimization here make sense?
      !
      implicit none
      real(WP), intent(in) :: a(:)
      integer, intent(out) :: locs(:)
      ! *** end of interface ***

      integer :: i
      logical :: mask(size(a))

      mask = .true.

      do i = 1, size(locs)
         locs(i:i) = minloc(a, mask)
         mask(locs(i)) = .false.
      enddo
    end subroutine minlocs
  end subroutine makeSchedule


  pure function makeSpan(tasks, np) result(span)
    !
    ! Purpose:  Estimates   the  makespan   for  given  tasks   on  NP
    ! processors.  Implemented is  a LIST-scheduling of tasks(:), that
    ! is the tasks  in a particular sequence.  Each  task is scheduled
    ! as soon as the required number of processors is available.
    !
    ! FIXME:  the algorithm  here should  be consistent  with  the one
    ! which generates detailed schedule description in makeSchedule().
    !
    ! FIXME: O(N^2) complexity with  N = size(tasks). This function is
    ! called N! times.
    !
    implicit none
    type(npts_task), intent(in) :: tasks(:)
    integer, intent(in) :: np
    real(WP) :: span ! total time
    ! *** end of interface ***

    integer(IP) :: i, k, P, Q, loc(1)

    ! Two arrays that store the  time and the number of slots becoming
    ! available. Each position in those  two array will be referred to
    ! as a  "bank" below. The theory  is one does not  need more banks
    ! than an the number of  tasks --- each completed task potentially
    ! corresponds to at most one bank of some number of slots becoming
    ! available:
    real(WP) :: times(size(tasks))
    integer(IP) :: slots(size(tasks))
    real(WP) :: time

    span = 0.0

    !
    ! We represent empty banks in this way:
    !
    times = huge(time)
    slots = 0

    !
    ! At the beginning  all slots are available at  the same time 0.0,
    ! we put  them in the first  bank (albeit after we  enter the loop
    ! over tasks):
    !
    P = np
    time = 0.0

    do i = 1, size(tasks)
       !
       ! Create a  new bank with P  = "P" slots  becoming available at
       ! the time = "time" (time of  the previous task, or just 0.0 in
       ! the first iteration).  FIXME: Up to this point in time all of
       ! the  chosen processors  will be  considered blocked,  even if
       ! some of them are idle inbetween.  For that find an empty bank
       ! and update it:
       !
       k = 1
       do while (slots(k) /= 0)
          k = k + 1
          ! The  theory is  that the  number of  non-empty  banks will
          ! increase by at most one for each task. So that the initial
          ! allocation of size(seq) will suffice:
       enddo
       slots(k) = P
       times(k) = time

       !
       ! Now   that   the  banks   were   updated,   proceed  to   the
       ! reservation. Number of slots to be reserved for this task:
       !
       P = tasks(i)%slots

       !
       ! Make a  reservation by taking slots from  the banks available
       ! at an earlierst  point.  Thus, in the last  iteration we will
       ! find  the time  point  for the  "task(i)"  to have  all of  P
       ! processors available to start the execution.
       !
       ! This will be decremented as we reserve more and more slots:
       !
       Q = P
       do while (Q > 0)
          !
          ! Find the slots to be available as soon as possible, loc is
          ! a 1-sized  array. FIXME: dynamically adjust  the number of
          ! banks, maybe?
          !
          loc = minloc(times)

          ! this is an index into arrays slots(:) and times(:):
          k = loc(1)

          ! update the estimate for the earliest time the task may
          ! start:
          time = times(k)

          if (Q < slots(k)) then
             ! reserve Q slots from that bank, bank remains non-empty,
             ! reservation is finished:
             slots(k) = slots(k) - Q
             Q = 0
             ! times(k) for remaining slots(k) is unchanged ...
          else
             ! reserve all slots from the bank, bank becomes empty,
             ! reservation may or may not be finished:
             Q = Q - slots(k)
             slots(k) = 0
             times(k) = huge(time)
          endif
       enddo

       !
       ! The "task(i)" will start  execution at "time" computed during
       ! the reservation above and finish here:
       !
       time = time + tasks(i)%span

       !
       ! Check  if  the completion  time  of  the  new scheduled  task
       ! exceeds the makespan. If yes, update the makespan:
       !
       if (time > span) span = time

       !
       ! Put P processors  back to new CPU bank  with the availability
       ! time computed above. See the beginning of the loop body ...
       !
    end do
  end function makeSpan


  !*************************************************************
  recursive subroutine sort_r(A)
    !
    ! Purpose: Implementation of QuickSort for Real arrays.
    !
    implicit none
    real(WP), intent(inout) :: A(:)
    !** End of interface *****************************************

    integer(IP) :: iq

    !
    ! End of recursion, zero-sized and length-one array are sorted:
    !
    if (size(A) <= 1) return

    call partition_r(A, iq)
    call sort_r(A(:iq-1))
    call sort_r(A(iq:))
  end subroutine sort_r
  !*************************************************************


  !*************************************************************
  subroutine partition_r(A, marker)
    !  Purpose: Helper routine for 'sort_r'
   implicit none
    !------------ Declaration of formal parameters ---------------
    real(WP),    intent(inout) :: A(:)
    integer(IP), intent(out)   :: marker
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(IP) :: i, j
    real(WP)    :: temp
    real(WP)    :: x      ! pivot point
    !------------ Executable code --------------------------------

    x = A(1)
    i= 0
    j= size(A) + 1

    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine partition_r
  !*************************************************************


  !*************************************************************
  recursive subroutine sort_i(A)
    !
    ! Purpose: Implementation of QuickSort for Integer arrays.
    !
    implicit none
    integer(IP), intent(inout) :: A(:)
    !** End of interface *****************************************

    integer(IP) :: iq

    !
    ! End of recursion, zero-sized and length-one array are sorted:
    !
    if (size(A) <= 1) return

    call partition_i(A, iq)
    call sort_i(A(:iq-1))
    call sort_i(A(iq:))
  end subroutine sort_i
  !*************************************************************


  !*************************************************************
  subroutine partition_i(A, marker)
    !  Purpose: Helper routine for 'sort_i'
   implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IP), intent(inout) :: A(:)
    integer(IP), intent(out)   :: marker
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(IP) :: i, j
    integer(IP) :: temp
    integer(IP) :: x      ! pivot point
    !------------ Executable code --------------------------------

    x = A(1)
    i= 0
    j= size(A) + 1

    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine partition_i
  !*************************************************************


end module se_scheduling_module

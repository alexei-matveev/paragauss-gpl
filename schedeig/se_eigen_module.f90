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
module se_eigen_module
  !---------------------------------------------------------------
  !
  !  Purpose:  Computes  all eigenvalues  and  -vectors  of the  given
  !  matrices from 'ham_tot' according  to the scheduling given in the
  !  structure  'scheduling'.   LAPACK   and  ScaLAPACK  routines  are
  !  employed to perform the eigenvalue computations.
  !
  !
  !
  !  Module called by: eigen_data_module
  !
  !
  !  Author: Martin Roderus
  !  Date: 04/2010
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
! define FPP_TIMERS 1
# include "def.h"
  use type_module ! type specification parameters
  use datatype
  USE_MPI
  use se_scheduling_module, only: taskType

  implicit none
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  public :: se_eigen_compeigs
  public :: se_eigen_close!()

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of interface blocks --------------------

  !------------ Declaration of constants and variables ----
  integer(i4_kind), parameter :: BLOCKSIZE = 64 ! ScaLAPACK block size

  ! MPI tags used in this module:
  integer(i4_kind), parameter :: TAG_A_BLOCK = 499
  integer(i4_kind), parameter :: TAG_RESULTS = 699

  integer(i4_kind) :: mpi_comm  ! MPI communicator
  integer(i4_kind) :: myrank    ! MPI rank of running process
  integer(i4_kind) :: nprocs    ! All available MPI processes

  !
  ! Private matrix/vector arrays:
  !

  ! Local array of all (partial) REAL matrices (irreps), in scheduling
  ! order:
  type(arrmat2), allocatable :: locHamTot(:)

  ! Local array of all (partial) REAL matrices (irreps), in scheduling
  ! order:
  type(arrmat2), allocatable :: locOverlapReal(:)

  ! Local array of all eigenvalues, in scheduling order
  type(arrmat1), allocatable :: locEigValArray(:)

  ! Local  array of  all  (partial) REAL  eigenvectors, in  scheduling
  ! order:
  type(arrmat2), allocatable :: locEigVecArray(:)

  !
  ! Cache schedule info between SCF/Geo iterations here:
  !
  type(taskType), allocatable, private :: oldSchedule(:)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  subroutine se_eigen_compeigs (H, S, e, V, mpi_communicator)
    !
    !  Purpose: public interface. Called from both, master and slave.
    !
    use se_scheduling_module, only: se_scheduling_run
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(arrmat3), intent(in) :: H(:)
    type(arrmat2), intent(in) :: S(:)
    type(arrmat2), intent(inout) :: e(:)
    type(arrmat3), intent(inout) :: V(:)
    integer(i4_kind), intent(in) :: mpi_communicator
    !** End of interface *****************************************

    integer(i4_kind) :: i, dims(size(H))
    integer(i4_kind) :: ierr
    FPP_TIMER_DECL(sched)

    ! Set global MPI communicator variable
    mpi_comm = mpi_communicator

    ! Query global MPI variables 'myrank' and 'nprocs'
    call MPI_Comm_rank (mpi_comm, myrank, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    call MPI_Comm_size (mpi_comm, nprocs, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    !
    ! Extract dimensions (overlap is not allocated on slaves?):
    !
    do i = 1, size (dims)
       dims(i) = size (H(i) % m, 1)
    enddo
    DPRINT "se_eigen_compeigs: dims=", dims, "on", myrank

    !
    ! First  propose a  schedule.  Combinatoric  approach  to schedule
    ! takes  significant   time  so  implement   caching  (FIXME:  and
    ! parallelize?).    Problem   sizes    can   be   retrieved   from
    ! oldSchedule(:)%size (?)
    !

    !
    ! If sizes differ, invalidate the cache:
    !
    if (allocated (oldSchedule)) then
       if (size (oldSchedule) /= size(dims)) then
          deallocate (oldSchedule)
       endif
    endif

    if (allocated (oldSchedule)) then
       if (any (oldSchedule(:) % size /= dims(oldSchedule(:) % id))) then
          deallocate (oldSchedule)
       endif
    endif

    !
    ! If there is no schedule, compute one:
    !
    if (.not. allocated (oldSchedule)) then
       !
       ! Prepare scheduling strategy for blocked Hamiltonian.  This is
       ! currently a serial code executed on all workers:
       !
       TRACE("se_scheduling_run/calling")
       DPRINT "se_eigen_compeigs: call se_scheduling_run(...)"
       FPP_TIMER_START(sched)
       call se_scheduling_run (oldSchedule, dims, nprocs)
       FPP_TIMER_STOP(sched)
       DPRINT "se_eigen_compeigs: time generating schedule =", FPP_TIMER_SLICE(sched)
       TRACE("se_scheduling_run/returned")
    endif

    !
    ! From here on a valid schedule  is supposed to exists in a global
    ! module    variable   (FIXME!)     named    (by   now    somewhat
    ! inappropriately) oldSchedule.  From that struct these fields are
    ! definitely used here, anything else?
    !
    !   oldSchedule(:)%id
    !   oldSchedule(:)%size
    !   oldSchedule(:)%allotedProcessors(:)
    !

    !
    ! Does the real work:
    !
    DPRINT "se_eigen_compeigs: call compEigSched"
    call compEigSched (oldSchedule, H, S, e, V)

    DPRINT "se_eigen_compeigs: call deallocateData()"
    ! Final deallocation
    call deallocateData()
    DPRINT "se_eigen_compeigs: exit"
  end subroutine se_eigen_compeigs
  !*************************************************************

  subroutine se_eigen_close()
    !
    ! Clean up the module.
    !
    implicit none
    ! *** end of interface ***

    !
    ! In  some cases,  eg. SO,  the parallel  eigensolver is  not used
    ! because  this  module  does  not  yet  support  solving  complex
    ! hermitean problems.  On the other hand, this procedure is called
    ! unconditionally   from  modules/initialization.f90.    So  check
    ! before attempting to deallocate:
    !
    if (allocated(oldSchedule)) deallocate (oldSchedule)
  end subroutine se_eigen_close


  !*************************************************************
  subroutine compEigSched (tasks, H, S, e, V)
    !
    !  Purpose: Main routine, executed in parallel context.
    !
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(taskType), intent(in) :: tasks(:)
    type(arrmat3), intent(in) :: H(:)
    type(arrmat2), intent(in) :: S(:)
    type(arrmat2), intent(inout) :: e(:)
    type(arrmat3), intent(inout) :: V(:)
    !** End of interface *****************************************

    integer(i4_kind) :: k, spin, n_spin
    integer :: ctx(size (tasks))

    n_spin = size (H(1) % m, 3)
    DPRINT "compEigSched: rank=", myrank, "n_spin=", n_spin

    !
    ! This fills global variables with meaningful data ...
    !
    call allocateData (size (tasks))

    DPRINT "compEigSched: rank=", myrank, "generate BLACS handles"
    !
    ! This sets  the context  for every task  (not all of  them are
    ! valid BLACS contexts):
    !
    do k = 1, size (tasks)
       ctx(k) = make_context (tasks(k))
    enddo

    !
    ! FIXME: prefer SPMD-style, (single program for all workers).
    !        This Master/Slave paradigm is twice the code size, and
    !        orders of magnitude more confusing.
    !
    do spin = 1, n_spin

       DPRINT "compEigSched: rank=", myrank, "call scatterData()"
       do k = 1, size (tasks)

          ! FIXME: are they required or just "for the case"?
          call barrier ()

          call scatterData (H, S, V, k, tasks(k), ctx(k), spin)
       enddo

       DPRINT "compEigSched: rank=", myrank, "call execEigSolver()"
       do k = 1, size (tasks)
          if (myrank == 0) then
             call execEigSolverMaster (H, S, e, V, k, tasks(k), ctx(k), spin)
          else
             call execEigSolverSlave (H, S, e, V, k, tasks(k), ctx(k), spin)
          endif
       enddo

       DPRINT "compEigSched: rank=", myrank, "call gatherData()"
       do k = 1, size (tasks)
          call gatherData (e(tasks(k) % id) % m(:, spin), V(tasks(k) % id) % m(:, :, spin), &
               k, tasks(k), ctx(k))
       enddo
       DPRINT "compEigSched: rank=", myrank, "done spin", spin
    enddo

    ! BLACS context is not needed anymore, it is freed
    do k = 1, size (tasks)
       call destroy_context (ctx(k))
    enddo
  end subroutine compEigSched


  function make_context (task) result(ctx)
    !
    ! Generate  an BLACS handle  from information  about how  many and
    ! which  processors  are  alloted  to each  task.   Collective  on
    ! mpi_comm.  Appears  to return negative context  on those workers
    ! that are not part of the BLACS grid.
    !
    use f77_scalapack, only: blacs_gridmap
    implicit none
    type(taskType), intent(in) :: task
    integer :: ctx
    ! *** end of interface ***

    integer :: i, j, nrow, ncol
    integer, allocatable :: ranks(:, :)

    ! Grid shape for this task:
    call grid_shape (task, nrow, ncol)

    ! Put ranks  of those workers  that will be  in BLACS grid  into a
    ! rectangular  array. Rows/columns  are base-0.   Same  as reshape
    ! (task % allotedProcessors(:nrow * ncol), [nrow, ncol]) - 1
    allocate (ranks(0:nrow-1, 0:ncol-1))
    do j = 0, ncol - 1
       do i = 0, nrow - 1
          ranks(i, j) = grid_rank (task, i, j)
       enddo
    enddo

    ! Context   is  formally   of   intent(inout)  in   the  call   to
    ! blacs_gridmap()  --- on  input it  is a  "system"  (or default?)
    ! context, on  output it is the  grid context. It  also looks like
    ! all processes in mpi_comm have  to participate. Also I dont find
    ! an  exact quote but  it seems  that the  context is  negative on
    ! return for those workers which are not part of the grid:
    ctx = mpi_comm
    call blacs_gridmap (ctx, ranks, size (ranks, 1), nrow, ncol)

    ! We  use this  apparent feature  to screen  workers which  do not
    ! participate in a task:
    if (any (ranks == myrank)) then
       ASSERT(ctx>=0)
    else
       ASSERT(ctx<0)
    endif
  end function make_context

  subroutine destroy_context (ctx)
    use f77_scalapack, only: blacs_gridexit
    implicit none
    integer, intent(in) :: ctx
    ! *** end of interface ***

    ! Negative context means I am not part of that grid:
    if (ctx > 0) then
       call blacs_gridexit (ctx)
    else
       ! nothing
    endif
  end subroutine destroy_context

  subroutine scatter (task, ctx, A, B)
    !
    ! Distribute a (big) matrix A(:, :) from rank = 0 to local B(:, :)
    ! on all workers participating in the task.
    !
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    type(taskType), intent(in) :: task
    integer, intent(in) :: ctx           ! BLACS context
    real(r8_kind), intent(in) :: A(:, :) ! global
    real(r8_kind), allocatable, intent(out) :: B(:, :) ! local
    ! *** end of interface ***

    ! BLACS variables
    integer(i4_kind) :: myrow, mycol
    integer(i4_kind) :: nprow, npcol

    integer :: m, n             ! global dimensions

    ! Problem dimensions:
    m = size (A, 1)
    n = size (A, 2)

    ! Allocate   local  partial   matrix,  this   has  an   effect  of
    ! deallocating it for ctx < 0:
    call reallocate2 (ctx, m, n, B)

    ! I am not involved:
    if (ctx < 0 .and. myrank /= 0) return

    if (ctx >= 0) then  ! I am involved, allocate local matrix
       ! Get my row- and column index in the BLACS process grid
       call blacs_gridinfo (ctx, nprow, npcol, myrow, mycol)
    else
       ! The case  for master  when it does  not participate  in task.
       ! Not all of np workers form a BLACS grid:
       call grid_shape (task, nprow, npcol)
       myrow = -1
       mycol = -1
    end if

    if (nprow == 1 .and. npcol == 1) then
       !
       ! "Optimization" was  there in the original  version.  Send all
       ! in one block.  In fact, there is no  evidence whatsoever that
       ! this is better/faster.
       !
       call block_cyclic (m, n)
    else
       call block_cyclic (BLOCKSIZE, BLOCKSIZE)
    endif

  contains

    subroutine block_cyclic (mbs, nbs)
      implicit none
      integer, intent(in) :: mbs, nbs ! block sizes
      ! *** end of interface ***

      integer :: mb, nb           ! number of blocks
      integer :: i, j             ! block indices
      integer :: ia, ja           ! block corner in A
      integer :: ib, jb           ! block corner in B
      integer :: ni, nj           ! block dimensions
      integer :: row, col         ! block owner

      ! Number  of blocks  along each  dimension, including  the last,
      ! eventually smaller one:
      mb = (m + mbs - 1) / mbs
      nb = (n + nbs - 1) / nbs

      ! Send cyclic blocks over the grid
      do j = 0, nb - 1
         do i = 0, mb - 1
            !
            ! Compute base-0  offsets into array A, section  B and the
            ! block sizes. First block corner in A:
            !
            ia = i * mbs
            ja = j * nbs

            ! Block corner in B:
            ib = i / nprow * mbs
            jb = j / npcol * nbs

            ! Block sizes. Last blocks along  i and j may have smaller
            ! sizes:
            ni = min (mbs, m - ia)
            nj = min (nbs, n - ja)

            ! (row, col)  are the grid coordinates of  the worker this
            ! block should be sent to:
            row = mod (i, nprow)
            col = mod (j, npcol)

            if (row == myrow .and. col == mycol) then
               !
               ! Mine block, fill the local partial matrix B:
               !
               if (myrank == 0) then
                  ! copy
                  B(ib + 1:ib + ni, jb + 1:jb + nj) = &
                       A(ia + 1:ia + ni, ja + 1:ja + nj)
               else
                  ! recv
                  call recv (B(ib + 1:ib + ni, jb + 1:jb + nj), &
                       0, TAG_A_BLOCK)
               endif
            else if (myrank == 0) then
               !
               ! Not mine block, send a section of A from master:
               !
               call send (A(ia + 1:ia + ni, ja + 1:ja + nj), &
                    grid_rank (task, row, col), TAG_A_BLOCK)
            else
               ! nothing to do on slaves for blocks that do not belong
               ! to them
            end if
         end do
      end do   ! Send cyclic blocks over the grid
    end subroutine block_cyclic
  end subroutine scatter


  subroutine gather (task, ctx, B, A)
    !
    ! Gather a (big) matrix A(:, :) on  rank = 0 from local B(:, :) on
    ! all workers participating in the task.
    !
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    type(taskType), intent(in) :: task
    integer, intent(in) :: ctx           ! BLACS context
    real(r8_kind), allocatable, intent(in) :: B(:, :) ! local
    real(r8_kind), intent(out) :: A(:, :)             ! global
    ! *** end of interface ***

    ! BLACS variables
    integer(i4_kind) :: myrow, mycol
    integer(i4_kind) :: nprow, npcol

    integer :: m, n             ! global dimensions

    ! Problem dimensions:
    m = size (A, 1)
    n = size (A, 2)

    ! I am not involved:
    if (ctx < 0 .and. myrank /= 0) return

    if (ctx >= 0) then  ! I am involved, allocate local matrix
       ! Get my row- and column index in the BLACS process grid
       call blacs_gridinfo (ctx, nprow, npcol, myrow, mycol)
    else
       ! The case  for master  when it does  not participate  in task.
       ! Not all of np workers form a BLACS grid:
       call grid_shape (task, nprow, npcol)
       myrow = -1
       mycol = -1
    end if

    if (nprow == 1 .and. npcol == 1) then
       !
       ! "Optimization" was  there in the original  version.  Send all
       ! in one block.  In fact, there is no  evidence whatsoever that
       ! this is better/faster.
       !
       call block_cyclic (m, n)
    else
       call block_cyclic (BLOCKSIZE, BLOCKSIZE)
    endif

  contains

    subroutine block_cyclic (mbs, nbs)
      implicit none
      integer, intent(in) :: mbs, nbs ! block sizes
      ! *** end of interface ***

      integer :: mb, nb           ! number of blocks
      integer :: i, j             ! block indices
      integer :: ia, ja           ! block corner in A
      integer :: ib, jb           ! block corner in B
      integer :: ni, nj           ! block dimensions
      integer :: row, col         ! block owner

      ! Number  of blocks  along each  dimension, including  the last,
      ! eventually smaller one:
      mb = (m + mbs - 1) / mbs
      nb = (n + nbs - 1) / nbs

      ! Send cyclic blocks over the grid
      do j = 0, nb - 1
         do i = 0, mb - 1
            !
            ! Compute base-0  offsets into array A, section  B and the
            ! block sizes. First block corner in A:
            !
            ia = i * mbs
            ja = j * nbs

            ! Block corner in B:
            ib = i / nprow * mbs
            jb = j / npcol * nbs

            ! Block sizes. Last blocks along  i and j may have smaller
            ! sizes:
            ni = min (mbs, m - ia)
            nj = min (nbs, n - ja)

            ! (row, col)  are the grid coordinates of  the worker this
            ! block should be sent to:
            row = mod (i, nprow)
            col = mod (j, npcol)

            if (row == myrow .and. col == mycol) then
               !
               ! Mine block, fill the local partial matrix B:
               !
               if (myrank == 0) then
                  ! copy
                  A(ia + 1:ia + ni, ja + 1:ja + nj) = &
                       B(ib + 1:ib + ni, jb + 1:jb + nj)
               else
                  ! send
                  call ssend (B(ib + 1:ib + ni, jb + 1:jb + nj), &
                       0, TAG_A_BLOCK)
               endif
            else if (myrank == 0) then
               !
               ! Not mine block, recv a section of A from a worker:
               !
               call recv (A(ia + 1:ia + ni, ja + 1:ja + nj), &
                    grid_rank (task, row, col), TAG_A_BLOCK)
            else
               ! nothing to do on slaves for blocks that do not belong
               ! to them
            end if
         end do
      end do   ! Send cyclic blocks over the grid
    end subroutine block_cyclic
  end subroutine gather


  subroutine scatterData (H, S, V, k, task, ctx, spin)
    !  Purpose: scatter matrices in a block-cyclic way over the processors
    implicit none
    type(arrmat3), intent(in) :: H(:)
    type(arrmat2), intent(in) :: S(:)
    type(arrmat3), intent(inout) :: V(:)
    integer(i4_kind), intent(in) :: k ! task index, 1-based
    type(taskType), intent(in) :: task
    integer, intent(in) :: ctx ! BLACS context
    integer(i4_kind), intent(in) :: spin
    !** End of interface *****************************************

    ! BLACS variables
    integer(i4_kind) :: N
    integer(i4_kind) :: np

    ! Problem size:
    N = task % size

    ! Number of workers reserved for this task:
    np = size (task % allotedProcessors)

    if (np == 1 .and. myrank == 0 .and. ctx >= 0) then
       !
       ! FIXME: Ugly special case: master does that task alone:
       !
       call reallocate2 (ctx, N, N, locOverlapReal(k) % m)
       locOverlapReal(k) % m = S(task % id) % m

       V(task % id) % m(:, :, spin) = H(task % id) % m(:, :, spin)
    else   ! More than one alloted process
       !
       ! Allocate and fill  the data with my blocks,  if I am involved
       ! in the task:
       !
       call scatter (task, ctx, H(task % id) % m(:, :, spin), locHamTot(k) % m)
       call scatter (task, ctx, S(task % id) % m(:, :), locOverlapReal(k) % m)

       ! FIXME: not sure if these are used in all corner cases:
       call reallocate2 (ctx, N, N, locEigVecArray(k) % m)
       call reallocate1 (ctx, N, locEigValArray(k) % m)
    end if  ! Only one alloted process
  end subroutine scatterData

  subroutine execEigSolverMaster (H, S, e, V, k, task, ctx, spin)
    !
    ! Purpose: execute the LAPACK-  and ScaLAPACK solvers according to
    ! the given scheduling.
    !
    use matrix_eigenval, only: dsygv90
    implicit none
    type(arrmat3), intent(in) :: H(:)
    type(arrmat2), intent(in) :: S(:)
    type(arrmat2), intent(inout) :: e(:)
    type(arrmat3), intent(inout) :: V(:)
    integer(i4_kind), intent(in) :: k ! task index, 1-based
    type(taskType), intent(in) :: task
    integer, intent(in) :: ctx ! BLACS context
    integer(i4_kind), intent(in)  :: spin
    !** End of interface *****************************************

    integer(i4_kind) :: n
    integer(i4_kind) :: np

    !
    ! Do  I participate  in  diagonalization of  this  block? If  not,
    ! nothing to be done:
    !
    if (ctx < 0) return

    n = task % size

    ! Number of workers reserved for this task:
    np = size (task % allotedProcessors)

    if (np == 1) then
       !
       ! Only one alloted process => sequential LAPACK routine:
       !
       call dsygv90 (V(task % id) % m(:, :, spin), &
            locOverlapReal(k) % m, &
            e(task % id) % m(:, spin))

       deallocate (locOverlapReal(k) % m)
    else
       !
       ! Multiple alloted processes => parallel ScaLAPACK routine
       !
       call geigs_scalapack (ctx, n, &
            locHamTot(k) % m, &
            locOverlapReal(k) % m, &
            e(task % id) % m(:, spin), &
            locEigVecArray(k) % m)

       deallocate (locHamTot(k) % m, locOverlapReal(k) % m)
    end if
  end subroutine execEigSolverMaster
  !*************************************************************


  !*************************************************************
  subroutine execEigSolverSlave (H, S, e, V, k, task, ctx, spin)
    !
    ! Purpose: execute the LAPACK-  and ScaLAPACK solvers according to
    ! the given scheduling.
    !
    use matrix_eigenval, only: dsygv90
    implicit none
    type(arrmat3), intent(in) :: H(:)
    type(arrmat2), intent(in) :: S(:)
    type(arrmat2), intent(inout) :: e(:)
    type(arrmat3), intent(inout) :: V(:)
    integer(i4_kind), intent(in) :: k ! task index, 1-based
    type(taskType), intent(in) :: task ! unused
    integer, intent(in) :: ctx ! BLACS context
    integer(i4_kind), intent(in)  :: spin ! unused
    !** End of interface *****************************************

    integer(i4_kind) :: n
    integer(i4_kind) :: np

    !
    ! Do  I participate  in  diagonalization of  this  block? If  not,
    ! nothing to be done:
    !
    if (ctx < 0) return

    n = task % size

    ! Number of workers reserved for this task:
    np = size (task % allotedProcessors)

    if (np == 1) then
       !
       ! Only one alloted process => sequential LAPACK routine:
       !
       call dsygv90 (locHamTot(k) % m, &
            locOverlapReal(k) % m, &
            locEigValArray(k) % m)

       deallocate (locOverlapReal(k) % m)
    else
       !
       ! Multiple alloted processes => parallel ScaLAPACK routine:
       !
       call geigs_scalapack (ctx, n, &
            locHamTot(k) % m, &
            locOverlapReal(k) % m, &
            locEigValArray(k) % m, &
            locEigVecArray(k) % m)

       deallocate (locHamTot(k) % m, locOverlapReal(k) % m)
    end if
  end subroutine execEigSolverSlave
  !*************************************************************

  subroutine geigs_scalapack (context, n, A, B, e, V)
    !
    ! Purpose: execute  the ScaLAPACK  solvers according to  the given
    ! scheduling
    !
    use f77_scalapack, only: numroc, blacs_gridinfo, descinit, pdsygvx
    implicit none
    integer, intent(in) :: context, n
    real(r8_kind), intent(inout) :: A(:, :), B(:, :) ! ScaLAPACK overwrites ...
    real(r8_kind), intent(out) :: e(:), V(:, :)
    ! *** end of interface ***

    integer(i4_kind)              :: info

    ! BLACS variables
    integer(i4_kind)              :: myrow, mycol   ! BLACS grid row- and column rank
    integer(i4_kind)              :: nprow, npcol   ! BLACS grid row- and column total number

    ! LAPACK (DSYEV) attributes
    integer(i4_kind), parameter   :: itype = 1
    real(r8_kind),    allocatable :: work(:)
    real(r8_kind)                 :: workdummy(10)
    character, parameter :: jobz = 'V', uplo = 'L'
    integer(i4_kind) :: lda, ldb, lwork = -1

    ! Additional ScaLAPACK parameters and variables
    integer(i4_kind), parameter :: ia=1, ja=1, ib=1, jb=1, iz=1, jz=1, il=1, iu=1
    character, parameter :: range='A'

    !  ORFAC  (global input)  Specifies which  eigenvectors  should be
    !          reorthogonalized.   Eigenvectors   that  correspond  to
    !          eigenvalues which are  within tol=ORFAC*norm(A) of each
    !          other  are  to be  reorthogonalized.   However, if  the
    !          workspace  is  insufficient  (see  LWORK), tol  may  be
    !          decreased until all eigenvectors to be reorthogonalized
    !          can be  stored in one  process.  No reorthogonalization
    !          will be done if ORFAC  equals zero.  A default value of
    !          10^-3 is  used if ORFAC  is negative.  ORFAC  should be
    !          identical on all processes.
    !
    ! FIXME: literal constants for eigenvalue spacing and "accidental"
    ! degeneracies here:
    !
    double precision, parameter :: orfac = 1.0d-12
    integer, parameter :: cluster_size = 7

    double precision, parameter :: vl = 0.0, vu = 0.0
    real(r8_kind),    parameter   :: abstol=0._r8_kind
    real(r8_kind),    allocatable :: gap(:)
    integer(i4_kind)              :: desca(9), descb(9), descz(9)
    integer(i4_kind)              :: liwork, iworkdummy(10), m, nz
    integer(i4_kind), allocatable :: iwork(:), ifail(:), iclustr(:)

    !
    ! Multiple alloted processes => parallel ScaLAPACK routine
    !

    ! Get my row- and column index in the BLACS process grid
    call blacs_gridinfo (context, nprow, npcol, myrow, mycol)
    lda = numroc (n, BLOCKSIZE, myrow, 0, nprow)
    ldb = lda

    ! Initialize ScaLAPACK matrix descriptors
    call descinit (desca, n, n, BLOCKSIZE, BLOCKSIZE, 0, 0, context, lda, info)
    call descinit (descb, n, n, BLOCKSIZE, BLOCKSIZE, 0, 0, context, ldb, info)
    call descinit (descz, n, n, BLOCKSIZE, BLOCKSIZE, 0, 0, context, lda, info)

    ! First run for workspace query
    lwork = -1
    liwork = -1

    allocate (iclustr(2*nprow*npcol), gap(nprow*npcol))
    allocate (ifail(n))

    call pdsygvx (itype, jobz, range, uplo, n, A, ia, ja, desca, B, &
         & ib, jb, descb, vl, vu, il, iu, abstol, m, nz, e, orfac, V, &
         & iz, jz, descz, workdummy, lwork, iworkdummy, liwork, ifail, iclustr, gap, info)

    if (info /= 0) then
       print *, "geigs_scalapack: info=", info
       ABORT("pdsygvx failed (dry run)")
    endif

    !
    ! Quote  from  http://www.netlib.org/scalapack/double/pdsyevx.f
    ! (reference implementation). On  output "if JOBZ='V' WORK(1) =
    ! optimal amount  of workspace required  to compute eigenvalues
    ! and   eigenvectors   efficiently   with  *no   guarantee   on
    ! orthogonality*."   Emphasis mine.  So the  value  returned in
    ! workdummy(1) may not be what one needs.
    !
    ! Another  quote:   "The  computed  eigenvectors   may  not  be
    ! orthogonal if the minimal  workspace is supplied and ORFAC is
    ! too small.   If you want  to guarantee orthogonality  (at the
    ! cost  of potentially  poor  performance) you  should add  the
    ! following to LWORK: (CLUSTERSIZE-1)*N"
    !
    lwork = int (workdummy(1), kind (lwork)) + (cluster_size - 1) * n
    liwork = iworkdummy(1)

    ! Allocate workspace
    allocate (work(lwork), iwork(liwork))

    ! Second run for actual diagonalisation
    call pdsygvx(itype, jobz, range, uplo, n, A, ia, ja, desca, B, &
         & ib, jb, descb, vl, vu, il, iu, abstol, m, nz, e, orfac, V, &
         & iz, jz, descz, work, lwork, iwork, liwork, ifail, iclustr, gap, info)

    if (info /= 0) then
       print *, "geigs_scalapack: info=", info
       ABORT("pdsygvx failed")
    endif
  end subroutine geigs_scalapack


  subroutine gatherData (e, V, k, task, ctx)
    !  Purpose: Collect computed eigenvalues and -vectors from slaves
    use f77_scalapack, only: blacs_gridinfo
    implicit none
    real(r8_kind), intent(out) :: e(:)
    real(r8_kind), intent(out) :: V(:, :)
    integer, intent(in) :: k
    type(taskType), intent(in) :: task
    integer, intent(in) :: ctx ! BLACS context
    !** End of interface *****************************************

    ! BLACS variables
    integer(i4_kind) :: myrow, mycol   ! BLACS grid row- and column rank
    integer(i4_kind) :: nprow, npcol   ! BLACS grid row- and column total number
    integer(i4_kind) :: np
    logical :: blacs_master

    ! Number of workers reserved for this task:
    np = size (task % allotedProcessors)

    if (ctx >= 0) then
       call blacs_gridinfo (ctx, nprow, npcol, myrow, mycol)
    else
       myrow = -1
       mycol = -1
    endif

    ! This one will be sending eigenvalues:
    blacs_master = (myrow == 0 .and. mycol == 0)

    !
    ! Eigenvalues:
    !
    if (myrank == 0) then
       if (.not. blacs_master) then
          call recv1 (e, MPI_ANY_SOURCE, TAG_RESULTS)
       else
          ! FIXME:  again  a  sigularity,  master  does  not  allocate
          ! locEigValArray(k)  %  m,  at  all.   Uses  final  location
          ! directly.  No copy needed.  No deallocation needed. Why do
          ! I feel so dirty?
       endif
    else
       if (blacs_master) then
          call ssend1 (locEigValArray(k) % m, 0, TAG_RESULTS)
       endif

       if (ctx >= 0) then
          deallocate (locEigValArray(k) % m)
       endif
    endif

    call barrier()

    !
    ! Eigenvectors:
    !
    if (np <= 1) then
       if (myrank == 0 .and. ctx >= 0) then
          ! FIXME:  when master  does  the work  alone  the result  is
          ! already  in   V.   locEigVecArray(k)  %  m   is  not  even
          ! allocated, argh!
       else
          ! FIXME:  LAPACK  vs.   ScaLAPACK singularity,  with  LAPACK
          ! eigenvectors on output are in locHamTot(k). I hate my job!
          call gather (task, ctx, locHamTot(k) % m, V)
       endif
    else
       call gather (task, ctx, locEigVecArray(k) % m, V)
    end if

    call barrier()
  end subroutine gatherData


  !*********************************

  !****     Helper Routines     ****

  !*********************************


  !*************************************************************

  subroutine reallocate2 (ctx, M, N, buf)
    !
    ! Reallocate buf(:,  :) to  hold the  local section of  the M  x N
    ! matrix  according  to  block-cyclic  distribution  algorithm  of
    ! ScaLAPACK.  On those workers  that supply invalid BLACS context,
    ! ctx < 0 this  subroutine effectively deallocates the intent(out)
    ! argument.
    !
    use f77_scalapack, only: blacs_gridinfo, numroc
    implicit none
    integer, intent(in) :: ctx  ! BLACS context
    integer, intent(in) :: M, N ! global dimensions
    real(r8_kind), intent(out), allocatable :: buf(:, :)
    ! *** end of interface ***

    integer :: nprow, npcol, myrow, mycol
    integer :: locM, locN

    ! Doing nothing has an effect of deallocating intent(out) buf:
    if (ctx < 0) return

    ! Get my row- and column index in the BLACS process grid
    call blacs_gridinfo (ctx, nprow, npcol, myrow, mycol)

    ! Query local row- and column number
    locM = numroc (M, BLOCKSIZE, myrow, 0, nprow)
    locN = numroc (N, BLOCKSIZE, mycol, 0, npcol)

    ! Allocate local partial matrix
    allocate (buf(locM, locN))
  end subroutine reallocate2

  subroutine reallocate1 (ctx, N, buf)
    !
    ! See reallocate2() above.
    !
    implicit none
    integer, intent(in) :: ctx  ! BLACS context
    integer, intent(in) :: N    ! global dimension
    real(r8_kind), intent(out), allocatable :: buf(:)
    ! *** end of interface ***

    ! Doing nothing has an effect of deallocating intent(out) buf:
    if (ctx < 0) return

    ! Allocate diagonal matrix in full:
    allocate (buf(N))
  end subroutine reallocate1


#ifdef WITH_GUILE
  subroutine grid_shape (task, nrow, ncol)
    !
    ! Returns  the grid  shape nrow  x ncol.   Input: total  number of
    ! processes p reserved  for the task.  Output: number  of row- and
    ! column processes of the BLACS process grid.
    !
    use scm, only: scm_t, scm_call, scm_car, scm_cdr, assignment(=)
    use scheme, only: scheme_lookup
    implicit none
    type (taskType), intent (in) :: task
    integer, intent (out) :: nrow, ncol
    ! *** end of interface ***

    type (scm_t) :: proc, np, pair

    !
    ! We   delegate   the   task    to   the   Scheme   procedure   in
    ! guile/scheduling.scm. It is an error when the function cannot be
    ! found:
    !
    proc = scheme_lookup ("grid-shape", mod="guile scheduling")

    ! Number of workers for this task (SCM int):
    np = size (task % allotedProcessors)

    !
    ! The Scheme  procedure is  supposed to return  a pair  of integer
    ! values (nrow . ncol):
    !
    pair = scm_call (proc, np)

    ! Defined assignment here (int <- SCM int):
    nrow = scm_car (pair)
    ncol = scm_cdr (pair)
  end subroutine grid_shape
#else
  subroutine grid_shape (task, nrow, ncol)
    !
    ! Returns  the grid  shape nrow  x ncol.   Input: total  number of
    ! processes p reserved  for the task.  Output: number  of row- and
    ! column processes of the BLACS process grid.  If p < 9, nprow = 1
    ! and npcol = p. Else nprow = npcol = floor(sqrt(p))
    !
    implicit none
    type (taskType), intent (in) :: task
    integer, intent (out) :: nrow, ncol
    ! *** end of interface ***

    integer :: np

    ! number of workers for this task:
    np = size (task % allotedProcessors)

    ! This is our convention for rectangular grids.  NOTE: nrow * ncol
    ! <= np, not  exactly equal. Not all of  the reserved workers will
    ! be a part of the BLACS grid:
    if (np < 9) then
       ! Less than 9 processors alloted to the task, generate 1-dimensional
       ! process grid
       nrow = 1
       ncol = np
    else
       ! 9 processors or more alloted, generate square process grid
       nrow = int (sqrt (real (np, r8_kind)), kind (nrow))
       ncol = nrow
    end if
  end subroutine grid_shape
#endif


  function grid_rank (task, row, col) result (rank)
    !
    ! Encodes the grid map.
    !
    implicit none
    integer, intent(in) :: row, col ! base-0
    type(taskType), intent(in) :: task
    integer :: rank
    ! *** end of interface ***

    integer :: nrow, ncol

    call grid_shape (task, nrow, ncol)

    ASSERT(row>=0)
    ASSERT(col>=0)
    ASSERT(row<nrow)
    ASSERT(col<ncol)

    ! base-0 ranks, from base-1 worker IDs:
    rank = task % allotedProcessors(1 + row + nrow * col) - 1
  end function grid_rank


  subroutine allocateData(taskCount)
    implicit none
    integer(i4_kind), intent(in) :: taskCount
    ! *** end of interface ***

    allocate (locHamTot(taskCount))
    allocate (locOverlapReal(taskCount))
    allocate (locEigVecArray(taskCount))
    allocate (locEigValArray(taskCount))
  end subroutine allocateData


  subroutine deallocateData

    deallocate (locHamTot)
    deallocate (locOverlapReal)
    deallocate (locEigVecArray)
    deallocate (locEigValArray)
  end subroutine deallocateData


  subroutine barrier ()
    implicit none
    ! *** end of interface ***

    integer :: ierr

    call MPI_Barrier(mpi_comm, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine barrier


  subroutine send (buf, dest, tag)
    implicit none
    real(r8_kind), intent(in) :: buf(:, :)
    integer, intent(in) :: dest, tag
    ! *** end of interface ***

    integer :: ierr

    call MPI_Send (buf, size (buf), MPI_DOUBLE_PRECISION, dest, tag, mpi_comm, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine send


  subroutine ssend (buf, dest, tag)
    !
    ! FIXME: check if send() isnt enough.
    !
    implicit none
    real(r8_kind), intent(in) :: buf(:, :)
    integer, intent(in) :: dest, tag
    ! *** end of interface ***

    integer :: ierr

    call MPI_Ssend (buf, size (buf), MPI_DOUBLE_PRECISION, dest, tag, mpi_comm, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine ssend


  subroutine ssend1 (buf, dest, tag)
    implicit none
    real(r8_kind), intent(in) :: buf(:)
    integer, intent(in) :: dest, tag
    ! *** end of interface ***

    integer :: ierr

    call MPI_Ssend (buf, size (buf), MPI_DOUBLE_PRECISION, dest, tag, mpi_comm, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine ssend1


  subroutine recv (buf, src, tag)
    implicit none
    real(r8_kind), intent(out) :: buf(:, :)
    integer, intent(in) :: src, tag
    ! *** end of interface ***

    integer :: ierr, status(MPI_STATUS_SIZE)

    call MPI_Recv (buf, size (buf), MPI_DOUBLE_PRECISION, src, tag, mpi_comm, status, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine recv


  subroutine recv1 (buf, src, tag)
    implicit none
    real(r8_kind), intent(out) :: buf(:)
    integer, intent(in) :: src, tag
    ! *** end of interface ***

    integer :: ierr, status(MPI_STATUS_SIZE)

    call MPI_Recv (buf, size (buf), MPI_DOUBLE_PRECISION, src, tag, mpi_comm, status, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine recv1

end module se_eigen_module

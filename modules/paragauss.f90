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
module paragauss
  !---------------------------------------------------------------
  !
  !  Purpose:  holds what  was previousely  in main.f90.   This  is to
  !  allow calling ParaGauss as a subroutine from more than one place.
  !
  !  Call sequence:
  !
  !     world = qm_init()       ! call once
  !     call qm_run(world)      ! call zero or more times
  !     call qm_finalize(world) ! call once
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
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: qm_init             ! () -> world
  public :: qm_run              ! (world)
  public :: qm_finalize         ! (world)
  public :: toggle_legacy_mode  ! ()

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine qm_run(world) bind(C)
    !
    !  Purpose:   this is the main routine that only enrolls to comm
    !             and decides if this is a master or slave run,
    !             calling main_master or main_slave
    !
    !
    !  Author: TB, FN
    !  Date: 10/95
    !
    use iso_c_binding
    use comm, only: comm_world_push, comm_barrier, comm_bcast
    use comm_module, only: comm_enroll, comm_end ! legacy comm routines
    use initialization, only: initialize_environment, finalize
    use dlb, only: dlb_init, dlb_finalize
    use energy_calc_module, only: get_energy
#ifdef WITH_GUILE
    use scheme, only: scheme_define
#endif
    implicit none
    integer(C_INT), intent(in), value :: world
    ! *** end of interface ***

    real(r8_kind) :: energy

    TRACE("qm_run/entered")
    !
    ! Communication layer has to use this communicator from now on:
    !
    call comm_world_push(world)

    !
    ! This initializes  the legacy communicaiton  layer, including the
    ! C-layer with the communicator converted to C-representation:
    !
    call comm_enroll(world)

    !
    ! DLB uses a DUP of the world saved in static module variable:
    !
    call dlb_init(world)

    !
    ! Actual work done within initialize_environment() includes:
    !
    ! - setting the prefix for error messages
    ! - initializing clock and timers
    ! - reading environment variables,
    ! - setting up directory paths,
    ! - opening output units.
    !
    call initialize_environment()

    !
    ! This runs  on all  workers. Thinks of  it as of  master-plan for
    ! everyone:
    !
    call main_master()

    !
    ! For the moment, this value is only meaningful on master, let him
    ! broadcast the number:
    !
    call get_energy(energy)
    call comm_bcast(energy)
#ifdef WITH_GUILE
    call scheme_define ("*energy*", energy)
#endif

    ! do cleaning up
    call finalize()

    call dlb_finalize()

    !
    ! Clean  up   the  legacy  communicaiton  layer   (does  NOT  call
    ! MPI_Finalize() if you wonder:
    !
    call comm_end()

    !
    ! In  scripts,  qm_run()  may  be  followed by  cleanup  or  other
    ! modifications of the filesystem, so in case someone is not ready
    ! with, say, closing output  files wait here. Alternative would be
    ! to require a barrier after each call to qm_run():
    !
    call comm_barrier()
    TRACE("qm_run/return")
  end subroutine qm_run

  function qm_init() result(world) bind(C)
    !
    ! Calls  MPI_INIT and  returns communicator  to be  used, normally
    ! that is MPI_COMM_WORLD.
    !
    use iso_c_binding
    use comm, only: comm_init
    use dlb, only: dlb_init, DLB_THREAD_REQUIRED
    USE_MPI, only: MPI_THREAD_SERIALIZED, MPI_THREAD_MULTIPLE
    USE_DEBUG, only: print_pids_and_sleep
    implicit none
    integer(C_INT) :: world ! usually MPI_COMM_WORLD
    ! *** end of interface ***

    integer :: req

    ! Calculate  the required  thread safety  level, so  far  only DLB
    ! needs to be considered therein:
    req = DLB_THREAD_REQUIRED

    ! As there will  be also some message passing  of the main program
    ! while   DLB  distributes  the   jobs,  in   case  of   DLB  with
    ! MPI_THREAD_SERIALIZED (only the  separated helper thread will do
    ! MPI), there  is need of  level MPI_THREAD_MULTIPLE, so  that the
    ! message passing of the main program will not cause problems:
    if (req == MPI_THREAD_SERIALIZED) req = MPI_THREAD_MULTIPLE

    !
    ! Start MPI (consider required  thread savety level), return world
    ! communicator (e.g. MPI_COMM_WORLD):
    !
    world = comm_init(req)

    ! To give you some time to attach the debugger by "gdb -p PID":
    DCALL print_pids_and_sleep(20) ! ... seconds

#ifdef WITH_GUILE
    !
    ! Make a few procedure available to the Scheme interpreter:
    !
    call qm_init_scheme()
#endif
  end function qm_init

  subroutine qm_finalize(world) bind(C)
    use iso_c_binding
    use comm, only: comm_finalize
    implicit none
    integer(C_INT), intent(in), value :: world
    ! *** end of interface ***

    ! finish MPI
    call comm_finalize()
  end subroutine qm_finalize


  function toggle_legacy_mode() result (flag)
    !
    ! Slaves enter the receive/dispatch loop in main_slave() and leave
    ! upon  receiving of  a message  sent  from here.  On slaves,  the
    ! return value is unconditionally  false.  On master, the function
    ! returns true/false  on odd/even calls, respectively.  To be used
    ! in constructs like this:
    !
    !   do while (toggle_legacy_mode())
    !     !
    !     ! Master  does some work,  eventually sending  slaves orders
    !     ! that  are handled  in main_slave().  Master  executes this
    !     ! block exactly once, slaves do not enter here.
    !     !
    !     ...
    !   enddo
    !
    ! The idea  is to  migrate to common  control flow on  all workers
    ! only  occasionally  entering  legacy  master/slave  mode  during
    ! transition period.
    !
    ! I have a dream that one  day all of the workers will execute the
    ! same code irrespective of their rank. ABOLISH THE SLAVERY!
    !
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, comm_all_other_hosts
    use msgtag_module, only: msgtag_finito
    implicit none
    logical :: flag
    ! *** end of interface ***

    ! On master, flips on every call:
    logical, save :: legacy_mode = .false.

    if (comm_rank() /= 0) then
       legacy_mode = .not. legacy_mode ! flip!
       !
       ! Exit upon receival of a msgtag_finito:
       !
       call main_slave()
       legacy_mode = .not. legacy_mode ! flip!
    else
       if (legacy_mode) then
          !
          ! Second call.  If we are already in  master/slave mode then
          ! tell   the  slaves   to   leave  that   mode  by   exiting
          ! main_slave():
          !
          call comm_init_send (comm_all_other_hosts, msgtag_finito)
          call comm_send()
       endif
       legacy_mode = .not. legacy_mode ! flip!
    endif

    ! Slaves flip it twice, thus always return false:
    flag = legacy_mode
  end function toggle_legacy_mode


#ifdef WITH_GUILE
  subroutine qm_init_scheme() bind(c)
    !
    ! Exports a  few auxiliary methods  to the Scheme  interpreter. At
    ! the moment qm_init_scheme() is  run from qm_init().  So that the
    ! gsubrs  it   exposes  only  become   available  after  executing
    ! (qm-init)  in  Scheme interpreter  OR  upon (use-modules  (guile
    ! paragauss)).
    !
    use iso_c_binding, only: c_funptr, c_funloc
    use scm, only: scm_t, scm_define_gsubr
    use xc_func, only: qm_xc
    use vdw_dft, only: vdw_dft_phi
    use grid_module, only: guile_pople_radius, guile_slater_radius, &
         guile_ionic_radius
#ifdef WITH_BGY3D_NON_GPL
    use bgy3d, only: bgy3d_init_scheme
#endif
    implicit none
    ! *** end of interface ***

    type(c_funptr) :: fun
    type(scm_t) :: proc

#ifdef  WITH_MATRIX_PARALLEL
    !
    ! FIXME: GFortran 4.3 fails if  you inline call to c_funloc(), use
    ! a temp var as a workaround:
    !
    fun = c_funloc (qm_test_eigensolver)
    proc = scm_define_gsubr ("qm-test-eigensolver", 1, 0, 0, fun)
#endif

    fun = c_funloc (qm_xc)
    proc = scm_define_gsubr ("qm-xc", 6, 0, 0, fun)

    fun = c_funloc (vdw_dft_phi)
    proc = scm_define_gsubr ("qm-vdw-dft-phi", 2, 0, 0, fun)

    fun = c_funloc (guile_pople_radius)
    proc = scm_define_gsubr ("pople-radius", 1, 0, 0, fun)

    fun = c_funloc (guile_slater_radius)
    proc = scm_define_gsubr ("slater-radius", 1, 0, 0, fun)

    fun = c_funloc (guile_ionic_radius)
    proc = scm_define_gsubr ("ionic-radius", 1, 0, 0, fun)

#ifdef WITH_BGY3D_NON_GPL
    call bgy3d_init_scheme()
#endif
  end subroutine qm_init_scheme

#ifdef  WITH_MATRIX_PARALLEL
  function qm_test_eigensolver (n) result (time) bind (c)
    !
    ! Scheme procedure to test the parallel eigensolver:
    !
    ! (qm-test-eigensolver n) -> time
    !
    use scm, only: scm_t, scm_to_int, assignment(=)
    type(scm_t), intent(in), value :: n ! SCM int
    type(scm_t) :: time ! SCM double
    ! *** end of interface ***

    time = test (scm_to_int (n)) ! defined assignment
  contains
    function test (n) result (time)
      USE_MPI, only: MPI_WTIME
      use matrix_parallel, only: rmatrix, rdmatrix, matrix, geigs
      implicit none
      integer, intent(in) :: n
      double precision :: time
      ! *** end of interface ***

      real(r8_kind) :: A_(n, n), B_(n, n)
      type(rmatrix) :: A, B, V
      type(rdmatrix) :: e
      integer :: i

      call random_number (A_)
      call random_number (B_)
      A_ = A_ + transpose (A_)

      forall (i=1:n) B_(i, i) = 1000.0

      A = matrix (A_)
      B = matrix (B_)

      time = MPI_WTIME ()
      call geigs (A, B, e, V)
      time = MPI_WTIME () - time
    end function test
  end function qm_test_eigensolver
#endif /* WITH_MATRIX_PARALLEL */
#endif /* WITH_GUILE */

  !--------------- End of module ----------------------------------
end module paragauss

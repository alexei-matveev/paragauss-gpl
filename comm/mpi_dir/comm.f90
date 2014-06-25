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
module comm
  !---------------------------------------------------------------
  !
  ! Hides the ugliest pieces of MPI behind a slimer interface.
  !
  ! Copyright (c) 2010-2013 Alexei Matveev
  ! Copyright (c) 2010-2011 Astrid Nikodem
  ! Copyright (c) 2011 Thomas Soini
  ! Copyright (c) 2012 Bo Li
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
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
#define COMM_MPI
#ifdef COMM_MPI
  USE_MPI, only: MPI_COMM_NULL, MPI_ANY_SOURCE, MPI_STATUS_SIZE
#endif
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  interface comm_send
    module procedure comm_send_int
    module procedure comm_send_int1D
    module procedure comm_send_double1D
  end interface

  public :: comm_send! (auto data, integer destination)

  interface comm_recv
    module procedure comm_recv_int
    module procedure comm_recv_int1D
    module procedure comm_recv_double1D
  end interface

  public :: comm_recv! (auto data, integer source)

  interface comm_bsend
    module procedure comm_bsend_int
    module procedure comm_bsend_int1D
    module procedure comm_bsend_double1D
  end interface

  public :: comm_bsend! (auto data, integer destination)

  interface comm_bcast
    module procedure comm_bcast_int
    module procedure comm_bcast_int1D
    module procedure comm_bcast_int2D
    module procedure comm_bcast_int3D
    module procedure comm_bcast_double
    module procedure comm_bcast_double1D
    module procedure comm_bcast_double2D
    module procedure comm_bcast_double3D
    module procedure comm_bcast_double4D
    module procedure comm_bcast_string
    module procedure comm_bcast_logical
  end interface

  public :: comm_bcast! (auto data, integer root)

  interface comm_allreduce
    module procedure comm_allreduce_int
    module procedure comm_allreduce_double
    module procedure comm_allreduce_double1D
    module procedure comm_allreduce_double2D
    module procedure comm_allreduce_logical
  end interface

  public :: comm_allreduce! (auto data)

  interface comm_minloc
    module procedure comm_minloc_double
  end interface

  public :: comm_minloc! (auto data) -> integer array

  interface comm_reduce
    module procedure comm_reduce_int
    module procedure comm_reduce_int1D
    module procedure comm_reduce_double
    module procedure comm_reduce_double1D
    module procedure comm_reduce_double2D
    module procedure comm_reduce_double3D
    module procedure comm_reduce_double4D
  end interface

  public :: comm_reduce! (auto data, integer root?, character op?)

  interface comm_gather
     module procedure comm_gatherv_double1D!(send, recv, counts)
     module procedure comm_gatherv_in_place_double1D!(buf, counts)
     module procedure comm_gather_double
  end interface

  public :: comm_gather

  interface comm_scatter
     module procedure comm_scatterv
  end interface

  public :: comm_scatter

  interface comm_same
     !
     ! For debug only, returns true if the argument is the same on all
     ! workers:
     !
     module procedure comm_same_int!(int) -> logical
     module procedure comm_same_logical!(logical) -> logical
     module procedure comm_same_double!(real) -> logical
     module procedure comm_same_double2D!(array) -> logical
  end interface

  public :: comm_same

  public :: comm_init!()
  public :: comm_finalize!()
  public :: comm_size!() -> integer np
  public :: comm_rank!() -> integer rank
  public :: comm_barrier!()
  public :: comm_parallel!() -> logical
  public :: comm_source!() -> integer src

  public :: comm_world_push!(integer world), DO NOT ABUSE!

#ifdef COMM_MPI
  !
  ! Set by comm_world_push() from comm_init(), do not use before:
  !
  integer, public, protected :: comm_world = MPI_COMM_NULL

  integer, parameter, public :: COMM_ANY_SOURCE = MPI_ANY_SOURCE
#endif
  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

#ifdef COMM_MPI
  ! status of the last receive:
  integer,            private :: stat(MPI_STATUS_SIZE) = -1
  integer, parameter, private :: default_tag = 123

  !
  ! If we call  MPI_INIT() from comm_init() we also  set this variable
  ! so that comm_finalize() can decide  whether to call or not to call
  ! MPI_FINALIZE():
  !
  logical, private :: we_called_mpi_init = .false.
#endif

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function comm_parallel() result(yes)
    implicit none
    logical :: yes ! result
    ! *** end of interface ***

    yes = comm_size() > 1
  end function comm_parallel

#ifdef COMM_MPI
  function comm_init (req) result(world)
    !
    ! Initializes  MPI and returns  MPI_COMM_WORLD. Idempotent  --- if
    ! MPI  is  already initialized  this  subroutine  only checks  the
    ! thread  safety  level  and  resets  the  world  communicator  to
    ! MPI_COMM_WORLD again.
    !
    USE_MPI, only: MPI_SUCCESS, MPI_COMM_WORLD, MPI_THREAD_SINGLE
!   USE_MPI, only: MPI_ERRORS_RETURN !, MPI_ERRHANDLER_SET
    implicit none
    integer(IK), intent(in) :: req
    integer :: world
    ! *** end of interface ***

    integer(IK) :: ierr
    integer(IK) :: prov
    logical :: initialized

    call MPI_INITIALIZED (initialized, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    if (req /= MPI_THREAD_SINGLE) then
       !
       ! Some parts of the code might  want to use threads this is for
       ! example  the case  for some  implementations of  DLB. Request
       ! MPI_INIT_THREAD() to give us the thread safety level:
       !
       if (.not. initialized) then
          call MPI_INIT_THREAD (req, prov, ierr)
       else
          call MPI_QUERY_THREAD (prov, ierr)
       endif
       ASSERT(ierr==MPI_SUCCESS)

       !
       ! We might get more or less than we expected if we got less the
       ! program is not  supposed to be stable thus  abort.  If we got
       ! more than we  wanted, just give a warning  but this should be
       ! all right,  as it  only means that  the special  thread level
       ! does not exist in the current MPI implementation:
       !
       if (prov < req) then
          print *, "Wanted", req, "but got only", prov
          ABORT("MPI thread safety level too low!")
       else if(prov > req) then
          WARN("MPI thread safety level higher than required!")
          print *, "Wanted", req, "but got instead", prov
       endif
    else
       !
       ! Stay  conservative, MPI_INIT_THREAD()  is  more fragile,  and
       ! segfaults on SuperMUC (as of 10/2011) in violation of the MPI
       ! standard saying that
       !
       !   "A  call to  MPI_INIT  has the  same  effect as  a call  to
       !    MPI_INIT_THREAD with a required = MPI_THREAD_SINGLE."
       !
       if (.not. initialized) then
          call MPI_INIT (ierr)
          ASSERT(ierr==MPI_SUCCESS)
       endif
    endif

    if (.not. initialized) then
       ! We  initialized MPI,  so  it is  also  our responsibility  to
       ! finalize it:
       we_called_mpi_init = .true.
    endif

    !
    ! To   eventually  continue   after  an   MPI  error   and  invoke
    ! error_handler(),  e.g. on assertion  failure, uncomment  the two
    ! lines below (and import above).
    !
!   call MPI_ERRHANDLER_SET(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
!   ASSERT(ierr==MPI_SUCCESS)

    !
    ! The  rest   of  the  program  and  libraries   should  use  this
    ! communicator, return to the outside world:
    !
    world = MPI_COMM_WORLD

    !
    ! To make sure the procedures in this module also use the
    ! same communicator set the private global variable:
    !
    call comm_world_push(world)
  end function comm_init

  subroutine comm_world_push(world)
    !
    ! Sets global module variable comm_world.
    !
    implicit none
    integer, intent(in) :: world
    ! *** end of interface ***

    comm_world = world
  end subroutine comm_world_push

  subroutine comm_finalize()
    USE_MPI, only: MPI_SUCCESS
    implicit none
    ! *** end of interface ***

    integer(IK) :: ierr
    logical :: finalized

    ! We initialized MPI, so it is also our responsibility to finalize
    ! it:
    if (we_called_mpi_init) then
       call MPI_FINALIZED (finalized, ierr)
       ASSERT(ierr==MPI_SUCCESS)

       if (.not. finalized) then
          call MPI_FINALIZE (ierr)
          ASSERT(ierr==MPI_SUCCESS)
       else
          ! somebody else already done that:
          WARN("MPI already finalized!")
       endif

       ! restore status quo:
       we_called_mpi_init = .false.
    endif
  end subroutine comm_finalize

  function comm_rank() result(rank)
    USE_MPI, only: MPI_SUCCESS
    implicit none
    integer(IK) :: rank ! result
    ! *** end of interface ***

    integer(IK) :: ierr

    call MPI_COMM_RANK( comm_world, rank, ierr )
    ASSERT(ierr==MPI_SUCCESS)
  end function comm_rank

  function comm_size() result(np)
    USE_MPI, only: MPI_SUCCESS
    implicit none
    integer(IK) :: np ! result
    ! *** end of interface ***

    integer(IK)           :: ierr

    call MPI_COMM_SIZE( comm_world, np, ierr )
    ASSERT(ierr==MPI_SUCCESS)
  end function comm_size

  function comm_source() result(src)
    USE_MPI, only: MPI_SOURCE
    implicit none
    integer(IK) :: src ! result
    ! *** end of interface ***

    src = stat(MPI_SOURCE)

    ! only meaningfull after a receive:
    ASSERT(src>=0)
  end function comm_source

  subroutine comm_barrier()
    USE_MPI, only: MPI_SUCCESS
    implicit none
    ! *** end of interface ***

    integer :: ierr

    call MPI_BARRIER(comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_barrier

  subroutine comm_send_int(i, dst)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(in) :: i, dst
    ! *** end of interface ***

    integer :: ierr

    call MPI_SEND(i, 1, MPI_INTEGER4, dst, default_tag, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_send_int

  subroutine comm_recv_int(i, src)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(out) :: i
    integer(IK), intent(in)  :: src
    ! *** end of interface ***

    integer :: ierr

    ! upodates global stat:
    call MPI_RECV(i, 1, MPI_INTEGER4, src, default_tag, comm_world, stat, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_recv_int

  subroutine comm_bsend_int(i, dst)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(in) :: i, dst
    ! *** end of interface ***

    integer :: ierr

    call MPI_BSEND(i, 1, MPI_INTEGER4, dst, default_tag, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bsend_int

  subroutine comm_send_int_buf(v, n, dst)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    integer :: ierr

    call MPI_SEND(v, n, MPI_INTEGER4, dst, default_tag, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_send_int_buf

  subroutine comm_recv_int_buf(v, n, src)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(out) :: v(*) ! (n)
    integer(IK), intent(in)  :: n, src
    ! *** end of interface ***

    integer :: ierr

    ! upodates global stat:
    call MPI_RECV(v, n, MPI_INTEGER4, src, default_tag, comm_world, stat, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_recv_int_buf

  subroutine comm_bsend_int_buf(v, n, dst)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    integer :: ierr

    call MPI_BSEND(v, n, MPI_INTEGER4, dst, default_tag, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bsend_int_buf

  subroutine comm_send_double_buf(v, n, dst)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK),    intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    integer :: ierr

    call MPI_SEND(v, n, MPI_DOUBLE_PRECISION, dst, default_tag, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_send_double_buf

  subroutine comm_recv_double_buf(v, n, src)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK),    intent(out) :: v(*) ! (n)
    integer(IK), intent(in)  :: n, src
    ! *** end of interface ***

    integer :: ierr

    ! upodates global stat:
    call MPI_RECV(v, n, MPI_DOUBLE_PRECISION, src, default_tag, comm_world, stat, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_recv_double_buf

  subroutine comm_bsend_double_buf(v, n, dst)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK),    intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    integer :: ierr

    call MPI_BSEND(v, n, MPI_DOUBLE_PRECISION, dst, default_tag, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bsend_double_buf

  subroutine comm_bcast_int(v, root)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(inout) :: v
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_BCAST(v, 1, MPI_INTEGER4, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bcast_int

  subroutine comm_bcast_double(v, root)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK), intent(inout) :: v
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_BCAST(v, 1, MPI_DOUBLE_PRECISION, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bcast_double

  subroutine comm_bcast_int_buf(v, n, root)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS
    implicit none
    integer(IK), intent(inout) :: v(*)
    integer(IK), intent(in) :: n, root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_BCAST(v, n, MPI_INTEGER4, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bcast_int_buf

  subroutine comm_bcast_double_buf(v, n, root)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK), intent(inout) :: v(*)
    integer(IK), intent(in) :: n, root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_BCAST(v, n, MPI_DOUBLE_PRECISION, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bcast_double_buf

  subroutine comm_bcast_string(s, root)
    USE_MPI, only: MPI_CHARACTER, MPI_SUCCESS
    implicit none
    character(len=*), intent(inout) :: s
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_BCAST(s, len(s), MPI_CHARACTER, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bcast_string

  subroutine comm_bcast_logical (s, root)
    USE_MPI, only: MPI_LOGICAL, MPI_SUCCESS
    implicit none
    logical, intent(inout) :: s
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_BCAST(s, 1, MPI_LOGICAL, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_bcast_logical

  function opcode(op) result(iop)
    USE_MPI, only: MPI_SUM, MPI_MAX, MPI_MIN, MPI_MAXLOC, MPI_MINLOC, &
         MPI_LAND, MPI_LOR, MPI_LXOR
    implicit none
    character(len=*), intent(in), optional :: op
    integer :: iop
    ! *** end of interface ***

    if (.not. present(op)) then
       iop = MPI_SUM
    else
       iop = -1
       select case (op)
       case ('sum')
          iop = MPI_SUM
       case ('max')
          iop = MPI_MAX
       case ('min')
          iop = MPI_MIN
       case ('maxloc')
          iop = MPI_MAXLOC
       case ('minloc')
          iop = MPI_MINLOC
       case ('and')
          iop = MPI_LAND
       case ('or')
          iop = MPI_LOR
       case ('xor')
          iop = MPI_LXOR
       case default
          ABORT('not implemented')
       end select
    endif
  end function opcode

  !
  ! MPI_IN_PLACE   in  reduce  operations   is  apparently   an  MPI-2
  ! feature. An MPI-1 implementation seems to require a temporary copy
  ! for reduce operations.
  !

#ifdef FPP_NO_MPI_IN_PLACE

  subroutine comm_allreduce_double_buf(v, n, op)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK), intent(inout)    :: v(*)
    integer(IK), intent(in)    :: n
    character(len=*), intent(in), optional :: op
    ! *** end of interface ***

    integer :: ierr
    real(RK) :: recv(n)

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    call MPI_ALLREDUCE (v, recv, n, MPI_DOUBLE_PRECISION, opcode(op), comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    v(1:n) = recv(1:n)
  end subroutine comm_allreduce_double_buf

  subroutine comm_reduce_int(v, root)
    USE_MPI, only: MPI_INTEGER4, MPI_SUM, MPI_SUCCESS
    implicit none
    integer(IK), intent(inout) :: v
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    integer(IK) :: root_rank, ierr
    integer(IK) :: recv ! FIXME: used only on root_rank

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_REDUCE(v, recv, 1, MPI_INTEGER4, MPI_SUM, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    if( comm_rank() == root_rank ) v = recv
  end subroutine comm_reduce_int

  subroutine comm_reduce_double(v, root, op)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK), intent(inout) :: v
    integer(IK), intent(in) :: root
    character(len=*), intent(in) :: op
    optional :: root, op
    ! *** end of interface ***

    integer(IK) :: root_rank, ierr, iop
    real(RK) :: recv ! FIXME: used only on root_rank

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    iop = opcode('sum')
    if(present(op)) iop = opcode(op)

    call MPI_REDUCE(v, recv, 1, MPI_DOUBLE_PRECISION, iop, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    if( comm_rank() == root_rank ) v = recv
  end subroutine comm_reduce_double

  subroutine comm_reduce_int_buf(v, n, root)
    USE_MPI, only: MPI_INTEGER4, MPI_SUM, MPI_SUCCESS
    implicit none
    integer(IK), intent(inout) :: v(*)
    integer(IK), intent(in)    :: n, root
    optional :: root
    ! *** end of interface ***

    integer(IK) :: root_rank, ierr
    integer(IK) :: recv(n) ! FIXME: used only on root_rank

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_REDUCE(v, recv, n, MPI_INTEGER4, MPI_SUM, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    if( comm_rank() == root_rank ) v(1:n) = recv(1:n)
  end subroutine comm_reduce_int_buf

  subroutine comm_reduce_double_buf(v, n, root)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_SUCCESS
    implicit none
    real(RK), intent(inout)    :: v(*)
    integer(IK), intent(in)    :: n, root
    optional :: root
    ! *** end of interface ***

    integer(IK) :: root_rank, ierr
    real(RK) :: recv(n) ! FIXME: used only on root_rank

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    root_rank = 0
    if(present(root)) root_rank = root

    call MPI_REDUCE(v, recv, n, MPI_DOUBLE_PRECISION, MPI_SUM, root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    if( comm_rank() == root_rank ) v(1:n) = recv(1:n)
  end subroutine comm_reduce_double_buf

#else /* that is we do have working MPI_IN_PLACE */

  subroutine comm_allreduce_double_buf(v, n, op)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUCCESS
    implicit none
    real(RK), intent(inout)    :: v(*)
    integer(IK), intent(in)    :: n
    character(len=*), intent(in), optional :: op
    ! *** end of interface ***

    integer :: ierr

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    call MPI_ALLREDUCE (MPI_IN_PLACE, v, n, MPI_DOUBLE_PRECISION, opcode(op), comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_allreduce_double_buf

  subroutine comm_reduce_int(v, root)
    USE_MPI, only: MPI_INTEGER4, MPI_SUM, MPI_IN_PLACE, MPI_SUCCESS
    implicit none
    integer(IK), intent(inout) :: v
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    if ( comm_rank() == root_rank ) then
      call MPI_REDUCE( MPI_IN_PLACE, v, 1, MPI_INTEGER4, MPI_SUM, root_rank, comm_world, ierr)
    else
      call MPI_REDUCE( v, MPI_IN_PLACE, 1, MPI_INTEGER4, MPI_SUM, root_rank, comm_world, ierr)
    endif
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_reduce_int

  subroutine comm_reduce_double(v, root, op)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUCCESS
    implicit none
    real(RK), intent(inout) :: v
    integer(IK), intent(in) :: root
    character(len=*), intent(in) :: op
    optional :: root, op
    ! *** end of interface ***

    integer :: root_rank, ierr, iop

    if( .not. comm_parallel() ) return

    root_rank = 0
    if(present(root)) root_rank = root

    iop = opcode('sum')
    if(present(op)) iop = opcode(op)

    if ( comm_rank() == root_rank ) then
      call MPI_REDUCE( MPI_IN_PLACE, v, 1, MPI_DOUBLE_PRECISION, iop, root_rank, comm_world, ierr)
    else
      call MPI_REDUCE( v, MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, iop, root_rank, comm_world, ierr)
    endif
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_reduce_double

  subroutine comm_reduce_int_buf(v, n, root)
    USE_MPI, only: MPI_INTEGER4, MPI_SUM, MPI_IN_PLACE, MPI_SUCCESS
    implicit none
    integer(IK), intent(inout) :: v(*)
    integer(IK), intent(in)    :: n, root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    root_rank = 0
    if(present(root)) root_rank = root

    if ( comm_rank() == root_rank ) then
      call MPI_REDUCE( MPI_IN_PLACE, v, n, MPI_INTEGER4, MPI_SUM, root_rank, comm_world, ierr)
    else
      call MPI_REDUCE( v, MPI_IN_PLACE, n, MPI_INTEGER4, MPI_SUM, root_rank, comm_world, ierr)
    endif
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_reduce_int_buf

  subroutine comm_reduce_double_buf(v, n, root)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE, MPI_SUCCESS
    implicit none
    real(RK), intent(inout)    :: v(*)
    integer(IK), intent(in)    :: n, root
    optional :: root
    ! *** end of interface ***

    integer :: root_rank, ierr

    ! dont call MPI with zero sized arrays:
    if(.not. comm_parallel() .or. n <= 0) return

    root_rank = 0
    if(present(root)) root_rank = root

    if ( comm_rank() == root_rank ) then
      call MPI_REDUCE( MPI_IN_PLACE, v, n, MPI_DOUBLE_PRECISION, MPI_SUM, root_rank, comm_world, ierr)
    else
      call MPI_REDUCE( v, MPI_IN_PLACE, n, MPI_DOUBLE_PRECISION, MPI_SUM, root_rank, comm_world, ierr)
    endif
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_reduce_double_buf

  function comm_same_int (v) result (yes)
    !
    ! Returns a true if v is the same on all procs. For debug only.
    !
    implicit none
    integer (IK), intent (in) :: v
    logical :: yes
    ! *** end of interface ***

    integer (IK) :: maxv

    maxv = v
    call comm_allreduce (maxv, "max")

    yes = (v == maxv)

    call comm_allreduce (yes, "and")
  end function comm_same_int

  function comm_same_logical (v) result (yes)
    !
    ! Returns a true if v is the same on all procs. For debug only.
    !
    implicit none
    logical, intent (in) :: v
    logical :: yes
    ! *** end of interface ***

    if (v) then
       yes = comm_same (1)
    else
       yes = comm_same (0)
    endif
  end function comm_same_logical

  function comm_same_double (v) result (yes)
    !
    ! Returns a true if v is the same on all procs. For debug only.
    !
    implicit none
    real(RK), intent(in) :: v
    logical :: yes
    ! *** end of interface ***

    real(RK) :: maxv

    maxv = v
    call comm_allreduce (maxv, "max")

    yes = (v == maxv)

    call comm_allreduce (yes, "and")
  end function comm_same_double

  function comm_same_double2D (v) result (yes)
    !
    ! Returns a true if v is the same on all procs. For debug only.
    !
    implicit none
    real(RK), intent(in) :: v(:, :)
    logical :: yes
    ! *** end of interface ***

    real(RK) :: maxv(size(v, 1), size(v, 2))

    maxv = v
    call comm_allreduce (maxv, "max")

    yes = all (v == maxv)

    call comm_allreduce (yes, "and")
  end function comm_same_double2D

  function comm_minloc_double(v) result(loc)
    !
    ! Returns a rank  of a process holding the minimal  value "v" in a
    ! size-1 array.
    !
    USE_MPI, only: MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_IN_PLACE, MPI_SUCCESS
    implicit none
    real(RK), intent(in) :: v
    integer :: loc(1)
    ! *** end of interface ***

    real(RK) :: pair(2)
    integer :: ierr

    pair(1) = v
    ! cast to double here:
    pair(2) = comm_rank()

    call MPI_ALLREDUCE (MPI_IN_PLACE, pair, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)

    ! Cast to integer here:
    loc(1) = int (pair(2), kind = kind (loc))
  end function comm_minloc_double

#endif

  subroutine comm_scatterv_double_buf(sendbuf, recvbuf, counts)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK), intent(in)    :: sendbuf(*)
    real(RK), intent(out)   :: recvbuf(*)
    integer(IK), intent(in) :: counts(:) ! (comm_size)
    ! *** end of interface ***

    integer :: root_rank = 0, rank, ierr
    integer :: displs(size(counts)), i

    rank = comm_rank()

    displs(1) = 0
    do i = 2, size(displs)
        displs(i) = displs(i-1) + counts(i-1)
    enddo

    call MPI_SCATTERV( sendbuf, counts, displs, MPI_DOUBLE_PRECISION &
                     , recvbuf, counts(1+rank), MPI_DOUBLE_PRECISION &
                     , root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_scatterv_double_buf

  subroutine comm_gatherv_double_buf(sendbuf, recvbuf, counts)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    real(RK), intent(in)    :: sendbuf(*)
    real(RK), intent(out)   :: recvbuf(*)
    integer(IK), intent(in) :: counts(:) ! (comm_size)
    ! *** end of interface ***

    integer :: root_rank = 0, rank, ierr
    integer :: displs(size(counts)), i

    rank = comm_rank()

    displs(1) = 0
    do i = 2, size(displs)
        displs(i) = displs(i-1) + counts(i-1)
    enddo

    call MPI_GATHERV( sendbuf, counts(1+rank), MPI_DOUBLE_PRECISION &
                    , recvbuf, counts, displs, MPI_DOUBLE_PRECISION &
                    , root_rank, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_gatherv_double_buf

  subroutine comm_gatherv_in_place_double_buf(buf, counts)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS, MPI_IN_PLACE
    implicit none
    real(RK), intent(inout) :: buf(*) ! intent(in) on non-root
    integer(IK), intent(in) :: counts(:) ! (comm_size)
    ! *** end of interface ***

    integer :: root_rank = 0, rank, ierr
    integer :: displs(size(counts)), i

    rank = comm_rank()

    displs(1) = 0
    do i = 2, size(displs)
        displs(i) = displs(i-1) + counts(i-1)
    enddo

    if ( rank == root_rank ) then
        ! send buf, send count and send type are not significant:
        call MPI_GATHERV( MPI_IN_PLACE, -1, MPI_DOUBLE_PRECISION &
                        , buf, counts, displs, MPI_DOUBLE_PRECISION &
                        , root_rank, comm_world, ierr)
    else
        ! recv buf, counts, displs, and recv type are not significant:
        call MPI_GATHERV( buf, counts(1+rank), MPI_DOUBLE_PRECISION &
                        , MPI_IN_PLACE, counts, displs, MPI_DOUBLE_PRECISION &
                        , root_rank, comm_world, ierr)
    endif
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_gatherv_in_place_double_buf

  subroutine comm_gather_double(SENDBUF, RECVBUF)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS
    implicit none
    !Note: Buffers are allocated by the calling programs
    real(RK), intent(in)    :: SENDBUF(:)
    real(RK), intent(inout) :: RECVBUF(:)
    ! *** end of interface ***

    integer :: ierr

    call MPI_GATHER(SENDBUF, 1, MPI_DOUBLE_PRECISION, RECVBUF, 1, MPI_DOUBLE_PRECISION, 0, comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_gather_double
#else
  function comm_init(req) result(world)
    implicit none
    integer(IK), intent(in) :: req
    integer :: world
    ! *** end of interface ***

    world = -1 ! MPI_COMM_NULL would be better
  end subroutine comm_init

  subroutine comm_finalize()
    implicit none
    ! *** end of interface ***

  end subroutine comm_finalize

  function comm_rank() result(rank)
    implicit none
    integer(IK) :: rank ! result
    ! *** end of interface ***

    rank = 0
  end function comm_rank

  function comm_size() result(np)
    implicit none
    integer(IK) :: np ! result
    ! *** end of interface ***

    np = 1
  end function comm_size

  function comm_source() result(src)
    implicit none
    integer(IK) :: src ! result
    ! *** end of interface ***

    src = 0
    ABORT("FIXME: serial version")
  end function comm_source

  subroutine comm_barrier()
    implicit none
    ! *** end of interface ***

  end subroutine comm_barrier

  subroutine comm_send_int(i, dst)
    implicit none
    integer(IK), intent(in) :: i, dst
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_send_int

  subroutine comm_recv_int(i, src)
    implicit none
    integer(IK), intent(out) :: i
    integer(IK), intent(in)  :: src
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_recv_int

  subroutine comm_bsend_int(i, dst)
    implicit none
    integer(IK), intent(in) :: i, dst
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_bsend_int

  subroutine comm_send_int_buf(v, n, dst)
    implicit none
    integer(IK), intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_send_int_buf

  subroutine comm_recv_int_buf(v, n, src)
    implicit none
    integer(IK), intent(out) :: v(*) ! (n)
    integer(IK), intent(in)  :: n, src
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_recv_int_buf

  subroutine comm_bsend_int_buf(v, n, dst)
    implicit none
    integer(IK), intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_bsend_int_buf

  subroutine comm_send_double_buf(v, n, dst)
    implicit none
    real(RK),    intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_send_double_buf

  subroutine comm_recv_double_buf(v, n, src)
    implicit none
    real(RK),    intent(out) :: v(*) ! (n)
    integer(IK), intent(in)  :: n, src
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_recv_double_buf

  subroutine comm_bsend_double_buf(v, n, dst)
    implicit none
    real(RK),    intent(in) :: v(*) ! (n)
    integer(IK), intent(in) :: n, dst
    ! *** end of interface ***

    ABORT("FIXME: serial version")
  end subroutine comm_bsend_double_buf

  subroutine comm_bcast_int(v, root)
    implicit none
    integer(IK), intent(inout) :: v
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_bcast_int

  subroutine comm_bcast_double(v, root)
    implicit none
    real(RK), intent(inout) :: v
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_bcast_double

  subroutine comm_bcast_int_buf(v, n, root)
    implicit none
    integer(IK), intent(inout) :: v(*)
    integer(IK), intent(in) :: n, root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_bcast_int_buf

  subroutine comm_bcast_double_buf(v, n, root)
    implicit none
    real(RK), intent(inout) :: v(*)
    integer(IK), intent(in) :: n, root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_bcast_double_buf

  subroutine comm_bcast_string(s, root)
    implicit none
    character(len=*), intent(inout) :: s
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_bcast_string

  subroutine comm_bcast_logical (s, root)
    implicit none
    logical, intent(inout) :: s
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_bcast_logical

  subroutine comm_reduce_int(v, root)
    implicit none
    integer(IK), intent(inout) :: v
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_reduce_int

  subroutine comm_reduce_double(v, root)
    implicit none
    real(RK), intent(inout) :: v
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_reduce_double

  subroutine comm_reduce_int_buf(v, n, root)
    implicit none
    integer(IK), intent(inout) :: v(*)
    integer(IK), intent(in)    :: n, root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_reduce_int_buf

  subroutine comm_reduce_double_buf(v, n, root)
    implicit none
    real(RK), intent(inout)    :: v(*)
    integer(IK), intent(in)    :: n, root
    optional :: root
    ! *** end of interface ***

  end subroutine comm_reduce_double_buf

  subroutine comm_reduce_maximum (double)
    !
    ! FIXME: better yet --- generalize comm_reduce
    !
    implicit none
    real(RK), intent(inout) :: double
    ! *** end of interface ***

  end subroutine comm_reduce_maximum

  subroutine comm_gather_double (SENDBUF, RECVBUF)
    implicit none
    real(RK),  intent(in)    :: SENDBUF(:)
    real(RK),  intent(inout) :: RECVBUF(:)
    ! *** end of interface ***

  end subroutine comm_gather_double
#endif

  subroutine comm_send_int1D(v, dst)
    implicit none
    integer(IK), intent(in) :: v(:)
    integer(IK), intent(in) :: dst
    ! *** end of interface ***

    call comm_send_int_buf(v, size(v), dst)
  end subroutine comm_send_int1D

  subroutine comm_recv_int1D(v, src)
    implicit none
    integer(IK), intent(out) :: v(:)
    integer(IK), intent(in)  :: src
    ! *** end of interface ***

    call comm_recv_int_buf(v, size(v), src)
  end subroutine comm_recv_int1D

  subroutine comm_bsend_int1D(v, dst)
    implicit none
    integer(IK), intent(in) :: v(:)
    integer(IK), intent(in) :: dst
    ! *** end of interface ***

    call comm_bsend_int_buf(v, size(v), dst)
  end subroutine comm_bsend_int1D

  subroutine comm_send_double1D(v, dst)
    implicit none
    real(RK),    intent(in) :: v(:)
    integer(IK), intent(in) :: dst
    ! *** end of interface ***

    call comm_send_double_buf(v, size(v), dst)
  end subroutine comm_send_double1D

  subroutine comm_recv_double1D(v, src)
    implicit none
    real(RK),    intent(out) :: v(:)
    integer(IK), intent(in)  :: src
    ! *** end of interface ***

    call comm_recv_double_buf(v, size(v), src)
  end subroutine comm_recv_double1D

  subroutine comm_bsend_double1D(v, dst)
    implicit none
    real(RK),    intent(in) :: v(:)
    integer(IK), intent(in) :: dst
    ! *** end of interface ***

    call comm_bsend_double_buf(v, size(v), dst)
  end subroutine comm_bsend_double1D

  subroutine comm_bcast_int1D(v, root)
    implicit none
    integer(IK), intent(inout) :: v(:)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_int_buf(v, size(v), root)
  end subroutine comm_bcast_int1D

  subroutine comm_bcast_int2D(v, root)
    implicit none
    integer(IK), intent(inout) :: v(:, :)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_int_buf(v, size(v), root)
  end subroutine comm_bcast_int2D

  subroutine comm_bcast_int3D(v, root)
    implicit none
    integer(IK), intent(inout) :: v(:, :, :)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_int_buf(v, size(v), root)
  end subroutine comm_bcast_int3D

  subroutine comm_bcast_double1D(v, root)
    implicit none
    real(RK), intent(inout) :: v(:)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_double_buf(v, size(v), root)
  end subroutine comm_bcast_double1D

  subroutine comm_bcast_double2D(v, root)
    implicit none
    real(RK), intent(inout) :: v(:,:)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_double_buf(v, size(v), root)
  end subroutine comm_bcast_double2D

  subroutine comm_bcast_double3D(v, root)
    implicit none
    real(RK), intent(inout) :: v(:,:,:)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_double_buf(v, size(v), root)
  end subroutine comm_bcast_double3D

  subroutine comm_bcast_double4D(v, root)
    implicit none
    real(RK), intent(inout) :: v(:,:,:,:)
    integer(IK), intent(in) :: root
    optional :: root
    ! *** end of interface ***

    call comm_bcast_double_buf(v, size(v), root)
  end subroutine comm_bcast_double4D

  subroutine comm_allreduce_int (v, op)
    USE_MPI, only: MPI_INTEGER4, MPI_SUCCESS, MPI_IN_PLACE
    implicit none
    integer (IK), intent (inout) :: v
    character (len=*), intent (in), optional :: op
    ! *** end of interface ***

    integer :: ierr

    call MPI_ALLREDUCE (MPI_IN_PLACE, v, 1, MPI_INTEGER4, opcode(op), comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_allreduce_int

  subroutine comm_allreduce_double (v, op)
    USE_MPI, only: MPI_DOUBLE_PRECISION, MPI_SUCCESS, MPI_IN_PLACE
    implicit none
    real(RK), intent(inout) :: v
    character(len=*), intent(in), optional :: op
    ! *** end of interface ***

    integer :: ierr

    call MPI_ALLREDUCE (MPI_IN_PLACE, v, 1, MPI_DOUBLE_PRECISION, opcode(op), comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_allreduce_double

  subroutine comm_allreduce_logical (v, op)
    USE_MPI, only: MPI_LOGICAL, MPI_SUCCESS, MPI_IN_PLACE
    implicit none
    logical, intent(inout) :: v
    character(len=*), intent(in), optional :: op
    ! *** end of interface ***

    integer :: ierr

    call MPI_ALLREDUCE (MPI_IN_PLACE, v, 1, MPI_LOGICAL, opcode(op), comm_world, ierr)
    ASSERT(ierr==MPI_SUCCESS)
  end subroutine comm_allreduce_logical

  subroutine comm_allreduce_double1D(v, op)
    implicit none
    real(RK), intent(inout) :: v(:)
    character(len=*), intent(in), optional :: op
    ! *** end of interface ***

    call comm_allreduce_double_buf(v, size(v), op)
  end subroutine comm_allreduce_double1D

  subroutine comm_allreduce_double2D(v, op)
    implicit none
    real(RK), intent(inout)    :: v(:, :)
    character(len=*), intent(in), optional :: op
    ! *** end of interface ***

    call comm_allreduce_double_buf(v, size(v), op)
  end subroutine comm_allreduce_double2D

  subroutine comm_reduce_int1D(v, root)
    implicit none
    integer(IK), intent(inout) :: v(:)
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    call comm_reduce_int_buf(v, size(v), root)
  end subroutine comm_reduce_int1D

  subroutine comm_reduce_double1D(v, root)
    implicit none
    real(RK), intent(inout)    :: v(:)
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    call comm_reduce_double_buf(v, size(v), root)
  end subroutine comm_reduce_double1D

  subroutine comm_reduce_double2D(v,root)
    implicit none
    real(RK), intent(inout)    :: v(:,:)
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    call comm_reduce_double_buf(v, size(v), root)
  end subroutine comm_reduce_double2D

  subroutine comm_reduce_double3D(v,root)
    implicit none
    real(RK), intent(inout)    :: v(:,:,:)
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    call comm_reduce_double_buf(v, size(v), root)
  end subroutine comm_reduce_double3D

  subroutine comm_reduce_double4D(v,root)
    implicit none
    real(RK), intent(inout)    :: v(:,:,:,:)
    integer(IK), intent(in)    :: root
    optional :: root
    ! *** end of interface ***

    call comm_reduce_double_buf(v, size(v), root)
  end subroutine comm_reduce_double4D

  subroutine comm_scatterv(send, recv, counts)
    implicit none
    real(RK), intent(in)     :: send(:)
    real(RK), intent(out)    :: recv(:)
    integer(IK), intent(in)  :: counts(:)
    ! *** end of interface ***

    integer :: rank

    rank = comm_rank()
    if ( rank == 0 ) then
      ASSERT(size(send)==sum(counts))
    endif
    ASSERT(size(recv)==counts(1+rank))

    call comm_scatterv_double_buf(send, recv, counts)
  end subroutine comm_scatterv

  subroutine comm_gatherv_double1D(send, recv, counts)
    implicit none
    real(RK), intent(in)     :: send(:)
    real(RK), intent(out)    :: recv(:)
    integer(IK), intent(in)  :: counts(:)
    ! *** end of interface ***

    integer :: rank

    WARN('untested')
    rank = comm_rank()
    if ( rank == 0 ) then
      ASSERT(size(recv)==sum(counts))
    endif
    ASSERT(size(send)==counts(1+rank))

    call comm_gatherv_double_buf(send, recv, counts)
  end subroutine comm_gatherv_double1D

  subroutine comm_gatherv_in_place_double1D(buf, counts)
    implicit none
    real(RK), intent(inout) :: buf(:) ! inout on root
    integer(IK), intent(in) :: counts(:)
    ! *** end of interface ***

    integer :: rank

    rank = comm_rank()
    if ( rank == 0 ) then
        ASSERT(size(buf)==sum(counts))
    else
        ASSERT(size(buf)==counts(1+rank))
    endif

    call comm_gatherv_in_place_double_buf(buf, counts)
  end subroutine comm_gatherv_in_place_double1D

  !--------------- End of module ----------------------------------
end module comm

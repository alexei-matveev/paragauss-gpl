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
module comm_module
  !---------------------------------------------------------------
  !comm_module.f90
  !
  !  Purpose: This  module contains all the necessary entities
  !           to deal with MPI in PARAGAU. Data on all running threads are
  !           kept und methods are defined that make the MPI user
  !           interface more convenient, taking away the necessity
  !           of external variables and error checks, for instance
  !           comm_send_master() and comm_mcast_slaves(). Direct calls
  !           to MPI routines are no longer necessary.
  !           The handling of error conditions with the necessity
  !           to display the error on the master and to terminate
  !           all running threads is done by the external subroutine
  !           error_handler that cooperates with the internal
  !           routine comm_save_recv. This routine should be called
  !           for any receive.
  !           See the list of public subroutines and functions below!
  !
  !  Contents:
  !
  !     - Public interfaces:  commpack  -> packing routines
  !                           communpack->  dito for unpacking
  !          (both defined in [mpi|pvm]pack_module)
  !
  !  Prerequesits for using the module:
  !      The program should enrol to the message passing demon by calling
  !      comm_enroll().
  !      The information on the virtual machine must be present in the slave
  !      threads
  !      as well, thus it must be send from master to slave; use
  !      commpack and communpack. Only after these steps, the other
  !      methods of the module will work.
  !      The slaves should have a main loop that tries to receive
  !      any message from the master and executes different tasks
  !      depending on the received message tag and exits if message
  !      comm_msgtag_finito id received.
  !
  !
  !
  !  Module called by: all routines dealing with MPI routines
  !
  !
  !  References: MPI man pages etc.
  !
  !  Author: TG
  !  Date: 8/97
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification
  ! Author: AS
  ! Date:   8/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !---------------------------------------------------------------

  !------------ Modules used --------------------------------------
# include "def.h"
#ifdef _COMPAC_FORTRAN
#  define MPI_COMM_RANK MPI_Comm_rank_
#  define MPI_COMM_SIZE mpi_comm_size_
#  define MPI_ABORT MPI_Abort_
#  define MPI_BARRIER MPI_Barrier_
#else
#  define MPI_COMM_RANK MPI_Comm_rank
#  define MPI_COMM_SIZE mpi_comm_size
#  define MPI_ABORT MPI_Abort
#  define MPI_BARRIER MPI_Barrier
#endif
  use type_module ! type specification parameters
  ! includes standard MPI fortran header file:
  USE_MPI, only: MPI_COMM_NULL, MPI_MAX_PROCESSOR_NAME
  use commpack_module ! interface for MPI packing and unpacking of data
  implicit none
  save ! save all variables defined in this module
  private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of constants public variables----------
  integer(kind=i4_kind), parameter, public :: comm_master_host = 1
  integer(kind=i4_kind), parameter, public :: comm_all_other_hosts = -1
  integer(kind=i4_kind), parameter, public :: comm_any_message = -1

  !------------ public functions and subroutines ------------------

  !
  ! All public nams should be prefixed with "comm_"
  !

  ! re-export commpack/communpack from commpack_module:
  public :: commpack
  public :: communpack

  public comm_enroll, &
       comm_abort, &
       comm_get_n_processors, &
       comm_i_am_master, &
       comm_hostname, comm_myindex, &
       comm_parallel, &
       comm_end, comm_print_conf, &
       comm_init_send, comm_save_recv_nonblocking

  public :: comm_save_recv

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----

  ! number of processors
  integer(kind=i4_kind), private  :: n_procs_cvm

  integer(i4_kind), parameter, private :: MAX_PROCESSOR_NAME = MPI_MAX_PROCESSOR_NAME

  ! hostnames of processors; name(1) master
  character(len=MAX_PROCESSOR_NAME), allocatable, private :: name_cvm(:)

  ! my own index
  integer(kind=i4_kind),private :: my_index_cvm

  logical, private :: i_am_master_cvm

#ifdef _DPRINT
  logical, private :: verbose = .false.
#endif

  !------------ Declaration of private functions and subroutines --

  !----------------------------------------------------------------

  !---------------- Interface statements ---------------------------
  ! The routines described by the following interfaces are defined in comm.c

  ! re-export them under new names:
  public :: comm_msgtag
  public :: comm_sendinghost
  public :: comm_send

  interface comm_msgtag
     integer(c_int) function comm_msgtag_buf() bind(C)
        use iso_c_binding
     end function comm_msgtag_buf
  end interface

  interface comm_sendinghost
     integer(c_int) function comm_sending_host_buf() bind(C)
        use iso_c_binding
     end function comm_sending_host_buf
  end interface

  interface
     subroutine comm_send() bind(C)
        use iso_c_binding
     end subroutine comm_send
  end interface

  !
  ! These are used within the module, not exported:
  !
  interface
     subroutine comm_init_buffer_data(world) bind(C)
        use iso_c_binding
        integer(c_int) :: world ! FIXME: should it be rather MPI_Fint?
     end subroutine comm_init_buffer_data

     subroutine comm_init_send_buf(target, msgtag) bind(C)
        use iso_c_binding
        integer(c_int) :: target, msgtag
     end subroutine comm_init_send_buf

     subroutine comm_save_recv_c(source, msgtag, err) bind(C)
        use iso_c_binding
        integer(c_int) :: source, msgtag, err
     end subroutine comm_save_recv_c

     integer(c_int) function comm_save_recv_nonblocking_buf(source, msgtag) bind(C)
        use iso_c_binding
        integer(c_int) :: source, msgtag
     end function comm_save_recv_nonblocking_buf
  end interface


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine comm_enroll(world)
    !
    !  Purpose: enroll to MPI.
    !
    !      This subroutine must be called before module
    !      is beeing used !!!!
    !
    ! Executed by all workers, initializes legacy communication
    ! layer.
    !
    ! Some of the code inside already uses the legacy "comm_send"
    ! interface!
    !
    implicit none
    integer, intent(in) :: world
    ! *** end of interface ***

    integer(kind=i4_kind)    :: ierror

    call MPI_COMM_RANK(world, my_index_cvm, ierror)
    ASSERT(ierror==0)

    ! Naturally, MPI uses indices ranging from 0 to n_procs-1.
    ! Cause the index scheme should be consistent to PVM,
    ! we shift to 1...n_procs.
    ! This is also assumed in comm_init_buffer_data():
    my_index_cvm = my_index_cvm + 1

    ! set the global private var "n_procs_cvm" here:
    call MPI_COMM_SIZE(world, n_procs_cvm, ierror)
    ASSERT(ierror==0)

    i_am_master_cvm = my_index_cvm .eq. 1

    ! Init communication buffers and basic variables in the
    ! C-module. Tell them to use the same communicator as we do:
    call comm_init_buffer_data(world)

    !
    ! From this point on the comm-interface is working ...
    !

    ! ask the system for my hostname and report that to the rest,
    ! only used for error reporting as far as I understand:
    call determine_hostnames()

    ! Send/recv  global  module  data  that was  initialized  only  on
    ! master.  Currently there is  no such  data, so  this is  just an
    ! example for early error detection:
    call send_recv_data()
  end subroutine comm_enroll

  !*************************************************************

  subroutine comm_print_conf(io_u)
    !  Purpose: prints mpi configuration
    !           determine_hostnames has to be called before
    !------------ Declaration of formal parameters ---------------
    implicit none
    integer(kind=i4_kind),intent(in)  :: io_u
    !** End of interface *****************************************
    !------------ Declaration of local variables ----------------
    integer(i4_kind) :: i, n
    !------------ Executable code --------------------------------
    write(io_u,*)
    write(io_u,*) 'MPI configuration:'
    write(io_u,*)
    write(io_u,*)'Number of processors is :', n_procs_cvm
    write(io_u,*)'Processor names are:'
    i = 1
    do while (i <= size(name_cvm))
      !
      ! Count repeated entries:
      !
      n = 0
      do while ( i + n <= size(name_cvm))
        if ( name_cvm(i) == name_cvm(i+n) ) then
          n = n + 1
        else
          exit ! inner while loop
        endif
      enddo

      !
      ! Compress output if there are repeated entries:
      !
      write(io_u, '(2X,(A), ":", I2)') trim(name_cvm(i)), n

      !
      ! Advance:
      !
      i = i + n
    enddo
    write(io_u,*)
  end subroutine comm_print_conf

  !*************************************************************

  subroutine send_recv_data()
    !
    ! Sending private information to  slaves. Currently only for early
    ! testing.
    !
    use msgtag_module, only: msgtag_commdata
    implicit none
    !** End of interface *****************************************

    integer (i4_kind) :: info, test

    if (comm_i_am_master()) then
      test = 123
      call comm_init_send (comm_all_other_hosts, msgtag_commdata)
      call commpack (test, info)
      ASSERT(info==0)
      call comm_send()
    else
      test = 0
      call comm_save_recv (comm_master_host, msgtag_commdata)
      call communpack (test, info)
      ASSERT(info==0)
    endif
    ASSERT(test==123)
  end subroutine send_recv_data

  !*************************************************************

  integer(kind=i4_kind) function comm_get_n_processors()
     ! purpose : returns number of processors
     !** End of interface *****************************************
     comm_get_n_processors = n_procs_cvm
  end function comm_get_n_processors

  !*************************************************************

  subroutine determine_hostnames()
    !  Purpose: name_cvm(:) is filled on ALL processors.
    !** End of interface *****************************************
    integer(kind=i4_kind) :: info, i
    integer(kind=i4_kind) :: strlen
    character(len=MAX_PROCESSOR_NAME) :: my_name

    ! each procs will have all hostnames:
    allocate(name_cvm(n_procs_cvm), STAT=info)
    ASSERT(info==0)

    call MPI_GET_PROCESSOR_NAME(my_name, strlen, info)
    if(info/=0)then
      print *,'determine_hostnames: error from mpi_get_processor_name, info=',info
      stop
    endif

    ! overwrite zeros (^@) returned:
    my_name(strlen+1:MAX_PROCESSOR_NAME) = ' ' ! FIXME: still needed?

    ! fill my array entry:
    name_cvm(my_index_cvm) = my_name

    !
    ! FIXME: I guess we expect sends to be non-blocking below.
    !

    ! send my index and my name to all other hosts:
    call comm_init_send(comm_all_other_hosts, 101)
    call commpack(name_cvm(my_index_cvm), info)
    ASSERT(info==0)
    call comm_send()

    ! and receive names of others:
    do i = 1, n_procs_cvm
       if ( i == my_index_cvm ) cycle

       ! recv from proc indexed by "p":
       call comm_save_recv(i, 101)
       call communpack(name_cvm(i), info)
       ASSERT(info==0)
    enddo
  end subroutine determine_hostnames

  !***************************************************************

  subroutine comm_init_send (target, msgtag)
    !
    ! Initialises send. First argument is only a dummy which is needed
    ! to keep the PVM counterpart happy.
    !
#ifdef _DPRINT
    use msgtag_module, only: msgtag_name
#endif
    implicit none
    integer (i4_kind), intent (in) :: target, msgtag
    !** End of interface *****************************************

#ifdef _DPRINT
    if (verbose) then
    print *, comm_myindex(), &
         'comm_init_send: tag ', trim (msgtag_name (msgtag)), &
         ' target=', target
    endif
#endif
    call comm_init_send_buf (target, msgtag)
  end subroutine comm_init_send

  !***************************************************************

  subroutine comm_abort(errorcode)
    !
    ! purpose: kill all MPI tasks belonging to the run
    !
    USE_MPI, only: MPI_COMM_WORLD
    implicit none
    integer(i4_kind), intent(in) :: errorcode
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------
    integer(kind=i4_kind) :: ierror

    ! we will keep MPI_COMM_WORLD here as for abort it might
    ! be more reasonable choice:
    call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierror)
    ABORT('returned from MPI_ABORT!')
  end subroutine comm_abort

  !***************************************************************

  subroutine comm_end()
    !
    ! Cleans up legacy communication layer (this module and C-layer)
    !
    implicit none
    !** End of interface *****************************************

    integer :: info

    deallocate(name_cvm, STAT=info)
    ASSERT(info==0)
  end subroutine comm_end

  !***************************************************************

  subroutine comm_save_recv (source, msgtag, info)
    !
    ! Receive a message, block until a matching one arrives
    !
    use msgtag_module, only: msgtag_error
#ifdef _DPRINT
    use msgtag_module, only: msgtag_name
#endif
    implicit none
    integer(kind=i4_kind),intent(in) :: source, msgtag
    integer(kind=i4_kind), optional, intent(out) :: info
    !** End of interface *****************************************

    integer(kind=i4_kind) :: err
    character(len=500) :: message

#ifdef _DPRINT
    if (verbose) then
       print *, comm_myindex(), &
            'comm_save_recv: waiting for tag ', trim (msgtag_name (msgtag)), &
            ' from=',source
    endif
#endif

    call comm_save_recv_c (source, msgtag, err)
    if (present (info)) info = err

#ifdef _DPRINT
    if (verbose) then
       print *, comm_myindex(), &
            'comm_save_recv: received tag ', trim (msgtag_name (comm_msgtag())), &
            ' from=', source
    endif
#endif

    ! FIXME: error handling code, possibly needs revision:
    if (msgtag == comm_any_message) then
       if (comm_msgtag() == msgtag_error) then
          call communpack(message, err)
          if (err /= 0) call error_handler &
               ("comm_save_recv: unpacking of error message failed")
          call error_handler (message)
       endif
    endif
  end subroutine comm_save_recv

  !***************************************************************

  logical function comm_save_recv_nonblocking (source, msgtag, info)
    !
    ! Try  to  receive,   but  do  not  block  if   there  is  nothing
    ! matching. Just return "false" in that case.
    !
    use msgtag_module, only: msgtag_error
#ifdef _DPRINT
    use msgtag_module, only: msgtag_name
#endif
    implicit none
    integer(kind=i4_kind), intent(in) :: source, msgtag
    integer(kind=i4_kind),optional,intent(out) :: info
    !** End of interface *****************************************

    integer(kind=i4_kind) :: recv, err
    character(len=500) :: message

    recv = comm_save_recv_nonblocking_buf (source, msgtag)
    ASSERT(recv==0.or.recv==1)
    if (present (info)) info = 0

    if (recv == 0) comm_save_recv_nonblocking = .false.

    if(recv == 1) then
       comm_save_recv_nonblocking = .true.

#ifdef _DPRINT
       if (verbose) then
         print *, comm_myindex(), &
              'comm_save_recv_nonblocking: received tag ', trim (msgtag_name (msgtag)), &
              ' from=', source
       endif
#endif

       ! FIXME: error handling code, possibly needs revision:
       if (msgtag == comm_any_message) then
          if (comm_msgtag() == msgtag_error) then
             call communpack (message, err)
             if (err /= 0) call error_handler &
                  ("comm_save_recv_nonblocking: unpacking of error message failed")
             call error_handler (message)
          endif
       endif
    endif
  end function comm_save_recv_nonblocking

  !***************************************************************

  logical function comm_i_am_master()
    ! purpose:  returns .true. if master
    !** End of interface *****************************************
    comm_i_am_master = i_am_master_cvm
  end function comm_i_am_master


  logical function comm_parallel()
    ! purpose:  returns .true. if parallel run
    !** End of interface *****************************************
    comm_parallel = ( n_procs_cvm .gt. 1 )
  end function comm_parallel

  !***************************************************************

  function comm_hostname() result(hostname)
    ! purpose:  returns hostname
    implicit none
    character(len=MAX_PROCESSOR_NAME) :: hostname
    !** End of interface *****************************************

    if ( allocated(name_cvm) ) then
       hostname = name_cvm(my_index_cvm)
    else
       hostname = ' '
    endif
  end function comm_hostname

  !***************************************************************

  integer function comm_myindex()
    !
    ! Purpose:   returns  index   of  present   host  between   1  and
    ! comm_get_n_processors(). In other words 1 + MPI_comm_rank().
    !
    implicit none
    !** End of interface *****************************************

    comm_myindex = my_index_cvm
  end function comm_myindex


#ifdef _DPRINT
  subroutine comm_set_verbose()
    implicit none
    !** End of interface *****************************************
    verbose = .true.
  end subroutine comm_set_verbose
#endif

  !--------------- End of module ----------------------------------
end module comm_module

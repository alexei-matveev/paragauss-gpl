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
module comm_agent
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: AM
  !  Date: 04/2009
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
  ! so far using the default fortran integers, might be better
  ! for interoperability with MPI:
  use type_module, only: IK=>i4_kind ! type specification parameters
  USE_MPI ! is made available by OpenMPI at least
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: comm_agent_fork!( world ) -> new_world or MPI_COMM_NULL
  public :: comm_agent_default!()
  public :: comm_agent_wait!()

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  integer            :: glob_universe = MPI_COMM_NULL ! == MPI_COMM_WORLD
  integer            :: glob_world    = MPI_COMM_NULL ! differs in two worlds
  integer            :: glob_pipe     = MPI_COMM_NULL ! inter-communicator

  integer, parameter :: AGENT_EXIT = 888

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function  comm_agent_fork( universe ) result( world )
    !
    ! If possible, return a new world smaller by one
    ! and save the inter-communicator in global module state.
    ! On agent the procedure returns MPI_COMM_NULL.
    ! If universe consists of one process return the universe
    ! itself.
    !
    implicit none
    integer, intent(in)  :: universe ! input communicator
    integer              :: world ! new intra-comm, MPI_COMM_NULL on agent
    ! *** end of interface ***

    integer :: color

    ! FIXME: universes cannot be nested so far:
    ASSERT(glob_universe==MPI_COMM_NULL)

    ! save the universe for the case of cases:
    glob_universe = universe

    ! workers return color == 1, the agent color == 0,
    ! save the inter-comm in global variable:
    color = fork( glob_universe, glob_world, glob_pipe )

    ! workers get new smaller world returned, agent nothing:
    if( color == 1 )then
      world = glob_world
    else
      world = MPI_COMM_NULL
    endif
  end function comm_agent_fork

  !*************************************************************

  subroutine comm_agent_wait()
    !
    ! called by workers, tells agent to return
    ! not really "wait", rather order ...
    !
    implicit none
    ! *** end of interface ***

    integer :: rank, ierr
    integer :: msg=999

    ! dont do anything if agent is not running:
    if( glob_pipe == MPI_COMM_NULL ) return

    ! this should act as intra-comm for our side:
    call MPI_COMM_RANK( glob_pipe, rank, ierr )
    ASSERT(ierr==0)

    ! send only one message:
    if( rank == 0 )then
      print *,'comm_agent_wait: terminating agent'
      call MPI_SEND( msg, 1, MPI_INTEGER, 0, AGENT_EXIT, glob_pipe, ierr )
      ASSERT(ierr==0)
      print *,'comm_agent_wait: agent terminated?'

      ! reset the global vars:
      call reset()
    endif
  end subroutine comm_agent_wait

  !*************************************************************

  subroutine comm_agent_default()
    implicit none
    ! *** end of interface ***

    integer :: ierr
    logical :: yes
    integer :: stat(MPI_STATUS_SIZE)
    integer :: msg, src, tag

    print *,'default agent started'
    do ! forever
      call MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, glob_pipe, yes, stat, ierr )
      ASSERT(ierr==0)

      if( .not. yes )then
        print *,'agent is idle'
        call sleep(1)
        cycle ! the probe loop
      endif

      print *,'agent sees a message, stat =',stat
      src = stat(MPI_SOURCE)
      tag = stat(MPI_TAG)

      call MPI_RECV( msg, 1, MPI_INTEGER, src, tag, glob_pipe, stat, ierr )
      ASSERT(ierr==0)

      select case( tag )
      case ( AGENT_EXIT )
        print *,'agent received AGENT_EXIT, exiting, BTW msg was',msg

        ! reset the global vars:
        call reset()

        exit ! the probe loop and return
      case default
        print *,'agent received tag =',tag,'msg =',msg,'need some work to do!'
      end select
    enddo
  end subroutine comm_agent_default

  !*************************************************************

  subroutine reset()
    implicit none
    ! *** end of interface ***

    ! FIXME: shouldnt we free the communicators here?
    glob_universe = MPI_COMM_NULL
    glob_world    = MPI_COMM_NULL
    glob_pipe     = MPI_COMM_NULL
  end subroutine reset

  !*************************************************************

  function fork( world, new_world, pipe) result(color)
    implicit none
    integer, intent(in)  :: world
    integer, intent(out) :: new_world
    integer, intent(out) :: pipe
    integer              :: color ! 1 for workers, 0 for agent
    ! *** end of interface ***

    integer :: n, p, peer, ierr

    call MPI_COMM_SIZE( world, n, ierr )
    ASSERT(ierr==0)

    if( n == 1 )then ! there is no other side
      color     = 1
      new_world = world
      pipe      = MPI_COMM_NULL
      RETURN
    endif

    call MPI_COMM_RANK( world, p, ierr )
    ASSERT(ierr==0)

    ! the processor with the highest rank will be the agent
    ! having different color, if there is more than one, that is:
    color = 1
    if( p == n - 1 ) color = 0

    call MPI_COMM_SPLIT( world, color, p, new_world, ierr )
    ASSERT(ierr==0)

    ! the first and the last processes are the respective
    ! peers for the two worlds:
    if( color == 1 )then
      peer = n - 1
    else
      peer = 0
    endif

    ! create the inter-world communicator:
    call MPI_INTERCOMM_CREATE( new_world, 0, world, peer, 999, pipe, ierr )
    ASSERT(ierr==0)
  end function fork

  !*************************************************************

  !--------------- End of module ----------------------------------
end module comm_agent

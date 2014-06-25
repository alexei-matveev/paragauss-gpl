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
  !comm_module.f90 - version for serial work
  !
  !    This module contains all necessary changes in communicated
  !    and information subroutines and functions which allow to
  !    use directly Paragaus for a work on single processor machine.
  !
  !  Author: AS, TB
  !  Date: 7/98
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

  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use commpack_module ! interface for packing and unpacking of comm data
  use msgtag_module

  implicit none
  save ! save all variables defined in this module

  !== Interrupt end of public interface of module =================



  integer(kind=i4_kind), parameter, public :: comm_master_host = 1
  integer(kind=i4_kind), parameter, public :: comm_all_other_hosts = -1
  integer(kind=i4_kind), parameter, public :: comm_all_slave_hosts = -2
  integer(kind=i4_kind), parameter, public :: comm_any_message = -1

  !------------ public functions and subroutines ------------------
  public comm_enroll, &
       comm_abort, &
       comm_get_n_processors, &
       comm_i_am_master, comm_hostname, comm_myindex, &
       comm_msgtag, &
       comm_parallel, &
       comm_end, comm_print_conf, &
       comm_write, comm_sendinghost, &
       comm_myname, &
       comm_init_send, comm_send, comm_save_recv, &
       comm_save_recv_nonblocking

  public :: comm_barrier

  !================================================================
  ! End of public interface of module
  !================================================================
  !------------ Declaration of private constants and variables ----
  integer(kind=i4_kind), private :: n_procs


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************

  function comm_enroll() result(return_color)
    !  Purpose: enrol to pvm and set private tid vars.
    !
    !      This subroutine must be called before module
    !      is beeing used !!!!
    !
    use filename_module, only: filename_set_hostindex
    implicit none
    integer :: return_color ! 1 for master 0 for slaves
    !** End of interface *****************************************

    return_color = 1
  end function comm_enroll


  !*************************************************************

  subroutine comm_print_conf(io_u)
    !  Purpose: prints pvm configuration
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(in)  :: io_u
    !** End of interface *****************************************
    !------------ Declaration of local variables ----------------
    !** End of interface *****************************************
    !------------ Executable code --------------------------------
    write(io_u,*)
    write(io_u,*) 'This is a seriall run:'
    write(io_u,*)
  end subroutine comm_print_conf

  !*************************************************************

  integer(kind=i4_kind) function comm_get_n_processors()
  ! purpose : returns number of processors
  !** End of interface *****************************************
  comm_get_n_processors = 1
  !** End of interface *****************************************
  end function comm_get_n_processors


 subroutine comm_init_send(index_of_host,msgtag)
   ! purpose: pass to "comm_send" message`s address and tag and
   ! call "pvmfinitsend" for initialization of send buffer
   !  checking for errors
   !------------ Declaration of formal parameters ---------------
   integer(kind=i4_kind), intent(in) :: index_of_host, msgtag
   !** End of interface *****************************************
   !print '(I3)', msgtag
  end subroutine comm_init_send

 !***************************************************************

 subroutine comm_send
  ! purpose: send data from any of hosts to any of one or all
  ! except itself, checking for errors
  !** End of interface *****************************************
 end subroutine comm_send

 !***************************************************************

 subroutine comm_save_recv(index_of_host,msgtag,bufid_out)
  !purpose: receive dummy
   integer(kind=i4_kind), intent(in) :: index_of_host, msgtag
   integer(kind=i4_kind), optional,intent(out) ::bufid_out
   !** End of interface *****************************************
   print '(I3)', msgtag
 end subroutine comm_save_recv

 !***************************************************************

 logical function comm_save_recv_nonblocking(index_of_host,msgtag,bufid_out)
   ! purpose: Trying to receive msgtag and then in case of being
   ! master trying to receive error message (tag msgtag_error).
   ! If an error_handler is received an errormessage is written
   ! and the calculation is terminatedon all hosts.
   ! Note that the configuration is not checked (i.e. it is not
   ! controlled if all hosts are alive). Call comm_check_config()
   ! manually in case you want that.
   ! Returns if msgtag was received.
   integer(kind=i4_kind), intent(in) :: index_of_host, msgtag
   integer(kind=i4_kind), optional,intent(out) ::bufid_out
   !** End of interface *****************************************
  end function comm_save_recv_nonblocking

 !***************************************************************

 subroutine comm_abort(errorcode)
   ! purpose: kill all other pvm processes belonging to the run
   implicit none
   integer(i4_kind), intent(in) :: errorcode
   !** End of interface *****************************************

   ABORT('why calling?')
 end subroutine comm_abort

 !***************************************************************

 subroutine comm_end()
   ! purpose: exits pvm
   !** End of interface *****************************************
 end subroutine comm_end

 !***************************************************************

 logical function comm_i_am_master()
   ! purpose:  returns .true. if master
   !** End of interface *****************************************
   comm_i_am_master = .true.
 end function comm_i_am_master


 logical function comm_parallel()
   ! purpose:  returns .true. if parallel run
   !** End of interface *****************************************
   comm_parallel = .false.
 end function comm_parallel


 !***************************************************************


 integer(kind=i4_kind) function comm_msgtag(bufid_in)
   ! purpose:  returns msgtag of currently received buffer by
   ! comm_save_recv or of buffer bufid_in
   !------------ Declaration of formal parameters ---------------
   integer(kind=i4_kind), optional, intent(in) ::bufid_in
   !** End of interface *****************************************
   comm_msgtag=0
 end function comm_msgtag

 !***************************************************************

 integer(kind=i4_kind) function comm_sendinghost(bufid_in)
   ! purpose:  returns index of sending host currently received
   ! buffer by comm_save_recv or of buffer bufid_in
   !------------ Declaration of formal parameters ---------------
   integer(kind=i4_kind), optional, intent(in) ::bufid_in
   !** End of interface *****************************************
   comm_sendinghost=0
 end function comm_sendinghost

 !***************************************************************

 character*30 function comm_hostname()
   ! purpose:  returns hostname
   !** End of interface *****************************************
   comm_hostname=' '
 end function comm_hostname

 !***************************************************************

 integer function comm_myindex()
    !
    ! Purpose:   returns  index   of  present   host  between   1  and
    ! comm_get_n_processors(). In other words always 1.
    !
    implicit none
    !** End of interface *****************************************

   comm_myindex = 1
 end function comm_myindex


 character*12 function comm_myname()
   ! same as comm_hostname, obsolete
   !** End of interface *****************************************
   comm_myname=' '
 end function comm_myname

  subroutine comm_barrier()
    implicit none
    !** End of interface *****************************************
  end subroutine comm_barrier

!--------------- End of module ----------------------------------
end module comm_module




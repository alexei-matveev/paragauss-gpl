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
subroutine error_handler(message)
!
!  Purpose: printing an error message and terminating
!           program. The message is preceded by:
!            "Error on <hostname> : "
!           All differences in treatment of the
!           error between master and slave are hidden
!           within this routine.
!
!  Author: TB
!  Date: 10/95
!
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!---------------------------------------------------------------------

!------------ Modules used --------------------------------------
#include "def.h"
use comm_module
use msgtag_module
use type_module
use iounitadmin_module, only: write_to_output_units
use time_module, only: usrtime, clktime
use error_module, only: MyID
USE_MEMLOG, only: memshow
implicit none
!== Interrupt end of public interface of module ================
!------------ Declaration of formal parameters ------------------
character(LEN=*) :: message
  !===================================================================
  ! End of public interface of module
  !===================================================================
!------------ Declaration of local variables --------------------
logical, save                       :: already_called = .false.
!------------ Executable code -----------------------------------
if ( already_called ) then
   print *, MyID//'recursive call of error_handler'
   call write_to_stderr(MyID//'recursive call of error_handler' )
   goto 999 ! and dont even try any fancy things ...
endif
already_called = .true.

if ( comm_i_am_master() ) then
   print *,MyID//'error on master'
else
   print *,MyID//'error on slave'
endif

call write_to_output_units( MyID//'Error on '//trim(comm_hostname())//': '//trim(message) )
call write_to_stderr( MyID//'Error on '//trim(comm_hostname())//': '//trim(message) )

#ifdef WITH_MEMLOG
  print *,MyID//'error_handler: memory stats before abort:'
  call memshow(2)
#endif

! print timing:
print *,MyID//'error_handler: clock time =',clktime()
print *,MyID//'error_handler:  user time =',usrtime()

!
! This should try to terminate all connected processes,
! comment this if you want a traceback output of TRACEBACKQQ:
!
call comm_abort(1)

999 CONTINUE
! dont try to call anything that may recursively
! invoke error_handler from this point on:
call terminate()

contains

  subroutine terminate
    ! stops calculation with return status 1
#ifdef FPP_TRACEBACKQQ /* Intel provides it */
    print *,'error_handler::terminate: call TRACEBACKQQ()'
    call TRACEBACKQQ() ! with Intel should provide the traceback
#else
    ! HP: use f90_unix, only: exit
    ! HP, SGI: call exit(1)
    ! IBM: call exit_(1_i4_kind)
    ! VPP: stop 1

    !
    ! To get s backtrace (with gfrotran) uncomment this:
    !
    ! call abort()

    stop ! just in case
#endif
  end subroutine terminate

  subroutine write_to_stderr(message)
    character(LEN=*), intent(in) :: message
    integer, parameter :: err_unit=0
    write(err_unit,*) message
  end subroutine write_to_stderr

end subroutine error_handler

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
subroutine  integral_shutdown_2cob3c()
!----------------------------------------------------------------
!
!  Purpose: Contains calls to shutdown routines to be executed
!           at end of 2 Center orbital integral and 3 Center
!           3 Center integral calculation. Is executed both
!           by master and slave.
!
!
!  Subroutine called by: integral_main_2cob3c, main_slave
!
!
!  Author: TB
!  Date: 5/96
!
!
!----------------------------------------------------------------
  !===================================================================
  ! End of public interface of module
  !===================================================================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: AS
! Date:   7/98
! Description: ...
!
! Modification (Please copy before editing)
! Author: MM
! Date:   9/97
! Description: Modifications for spin orbit run
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
#include <def.h>
#ifdef FPP_DEBUG
use error_module, only: MyID
#endif
use type_module ! type specification parameters
use output_module
use iounitadmin_module   ! to open output units
use timer_module
use time_module
use integralpar_module, only: integralpar_i_int_part,integralpar_send_3c
use integralpar_module, only: integralpar_pot_for_secderiv !!!!!!!!!!!AS
use options_module, only: options_spin_orbit
use int_send_2cob3c_module
use int_send_2cob3c_spor_module
use comm_module
use msgtag_module
use unique_atom_module
use spin_orbit_module,   only: is_on,op_FitTrafo
use int_send_aux_module, only: int_send_aux_done=>done
implicit none


!----------------------------------------------------------------
!------------ Executable code -----------------------------------


if (output_int_detailedprogress) call write_to_output_units( &
     "integral_shutdown_2cob3c:  begin")

! if necesarry ( in normal integral run), shutdown int_send_2cob3c_module
if(integralpar_send_3c .and. .not. integralpar_pot_for_secderiv) then !!!!!!!!!
   if (output_int_detailedprogress) call write_to_output_units( &
        "integral_shutdown_2cob3c: int_send_2cob3c_shutdown")
   if (options_spin_orbit) then
      !
      ! SPIN ORBIT
      !
      DPRINT MyID,'integral_shutdown_2cob3c: call int_send_2cob3c_spor_shutdown()'
      call int_send_2cob3c_spor_shutdown()
      DPRINT MyID,'integral_shutdown_2cob3c: .'
   else
      call int_send_2cob3c_shutdown()
   endif
end if

if(options_spin_orbit)then

   if(is_on(op_FitTrafo))then
      DPRINT MyID,'integral_shutdown_2cob3c: call int_send_aux_done()'
      call int_send_aux_done()
      DPRINT MyID,'integral_shutdown_2cob3c: .'
   endif
endif

if ( .not. comm_i_am_master() ) then
   call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
   call stop_timer(timer_int_2cob3c(integralpar_i_int_part))
endif


if (output_int_detailedprogress) call write_to_output_units( &
     "integral_shutdown_2cob3c:  end")

end subroutine integral_shutdown_2cob3c

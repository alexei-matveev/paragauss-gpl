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
subroutine  integral_setup_2cob3c(n_quad)
!----------------------------------------------------------------
!
!  Purpose: Contains calls to setup-routines to be executed
!           at beginning of 2 Center orbital and 3 Center
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
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

!------------ Modules used --------------------------------------
#include "def.h"
use type_module ! type specification parameters
use output_module, only: output_int_detailedprogress
use iounitadmin_module   ! to open output units
use timer_module
use time_module
use integralpar_module
use int_send_2cob3c_module, only: int_send_2cob3c_setup
use int_send_2cob3c_spor_module, only: int_send_2cob3c_spor_setup
use comm_module, only: comm_i_am_master
use spin_orbit_module,   only: is_on,op_FitTrafo
use int_send_aux_module, only: int_send_aux_init=>init
use options_module, only: options_spin_orbit

implicit none
integer(kind=i4_kind), intent(in) :: n_quad

!----------------------------------------------------------------
!------------ Executable code -----------------------------------

if (output_int_detailedprogress) call write_to_output_units( &
     "integral_setup_2cob3c:  begin")


if ( .not. comm_i_am_master() ) then
   call start_timer(timer_int_2cob3c(integralpar_i_int_part))
   call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
endif

! if necesarry, setup int_send_2cob3c_module
if(integralpar_send_3c .and. .not. integralpar_pot_for_secderiv) then 
   if (.not.integralpar_spor) then
      if (output_int_detailedprogress) call write_to_output_units( &
           "integral_setup_2cob3c: int_send_2cob3c_setup")
      call int_send_2cob3c_setup(n_quad)
   else
      if (output_int_detailedprogress) call write_to_output_units( &
           "integral_setup_2cob3c: int_send_2cob3c_spor_setup")
      call int_send_2cob3c_spor_setup(n_quad)
   endif
end if

if(options_spin_orbit)then
   if(is_on(op_FitTrafo))then
      call int_send_aux_init()
   endif
endif

if (output_int_detailedprogress) call write_to_output_units( &
     "integral_setup_2cob3c:  end")
end subroutine integral_setup_2cob3c

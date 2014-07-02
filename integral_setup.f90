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
subroutine  integral_setup()
!----------------------------------------------------------------
!
!  Purpose: Contains calls to setup-routines to be executed
!           at beginning of integral part.
!
!
!  Subroutine called by: main_integral, main_slave
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
use comm_module,   only: comm_parallel
use output_module, only: output_int_progress
use iounitadmin_module   ! to open output units
use integralpar_module, only: integralpar_send_receive, integralpar_i_int_part,integralpar_setup
use time_module          ! timing routines
use timer_module         ! timing database
use gamma_module, only: gamma_setup
implicit none

!----------------------------------------------------------------
!------------ Executable code -----------------------------------


if (output_int_progress) call write_to_output_units( &
     "integral_setup:  begin")

!if (output_int_progress) call write_to_output_units( &
!     "integral_setup:  call integralpar_setup()")
! moved to integralpar_set():
! call integralpar_setup()

if (output_int_progress) call write_to_output_units( &
     "integral_setup:  call integralpar_send_receive()")
 call integralpar_send_receive()
 call integralpar_setup()
! moved to integralpar_setup():
! integralpar_dervs=integralpar_dervs.and.integralpar_gradients

if (output_int_progress) call write_to_output_units( &
     "integral_setup:  gamma_setup")
call gamma_setup(16) ! cp. default numj=17 in old vers. of gamma_module

call start_timer(timer_int(integralpar_i_int_part))

if (output_int_progress) call write_to_output_units( &
     "integral_setup:  end")


end subroutine integral_setup

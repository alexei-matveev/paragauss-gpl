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
subroutine  integral_shutdown()
!----------------------------------------------------------------
!
!  Purpose: Contains calls to shutdown routines to be executed
!           at end of integral part.
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
use type_module ! type specification parameters
use output_module, only: output_int_progress
use iounitadmin_module   ! to open output units
use time_module          ! timing routines
use timer_module, only: timer_gather_slave_int_timing, timer_int
use integralpar_module, only: integralpar_i_int_part, &
     integralpar_relativistic
use gamma_module, only: gamma_close
use relgrads_store, only: rg_close
implicit none


!----------------------------------------------------------------
!------------ Executable code -----------------------------------


if (output_int_progress) call write_to_output_units( &
     "integral_shutdown:  begin")

if( integralpar_relativistic )then
  ! close the storage for relativistic matrices/gradients/derivatives:
  call rg_close('k') ! keep
endif

if (output_int_progress) call write_to_output_units( &
     "integral_shutdown:  gamma_close")

call gamma_close()

! stop timer
call stop_timer(timer_int(integralpar_i_int_part))


! send/receive timing of integral part to master
if (output_int_progress) call write_to_output_units( &
     "integral_shutdown: send/recv slave_int_timing")

call timer_gather_slave_int_timing(integralpar_i_int_part)

if (output_int_progress) call write_to_output_units( &
     "integral_shutdown:  end")


end subroutine integral_shutdown

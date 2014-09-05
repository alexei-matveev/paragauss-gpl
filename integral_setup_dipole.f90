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
subroutine  integral_setup_dipole()
  !-------------------------------------------------------------------
  !
  !  Purpose: Contains calls to setup-routines to be executed
  !           at beginning of dipole integral calculation.
  !
  !  Subroutine called by: integral_main_dipole, main_slave
  !
  !  Author: TB
  !  Date: 7/97
  !
  !
  !-------------------------------------------------------------------
  !===================================================================
  ! End of public interface of module
  !===================================================================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  !------------ Modules used -----------------------------------------
  use type_module ! type specification parameters
  use output_module, only: output_int_detailedprogress
  use iounitadmin_module
  use timer_module
  use time_module
  use integralpar_module
  use comm_module, only: comm_i_am_master

  implicit none


  !-------------------------------------------------------------------
  !------------ Executable code -----------------------------------


  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_setup_dipole:  begin")


  if ( .not. comm_i_am_master() ) then
     call start_timer(timer_int_2cob3c(integralpar_i_int_part))
     call start_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
  endif


  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_setup_dipole:  end")


end subroutine integral_setup_dipole

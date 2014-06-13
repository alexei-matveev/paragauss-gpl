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
subroutine  integral_shutdown_dipole()
  !----------------------------------------------------------------
  !
  !  Purpose: Contains calls to shutdown routines to be executed
  !           at end of dipole integral calculation.
  !
  !
  !  Subroutine called by: integral_main_dipole, main_slave
  !
  !
  !  Author: TB
  !  Date: 7/97
  !
  !
  !----------------------------------------------------------------
  !================================================================
  ! End of public interface of module
  !================================================================
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
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use output_module, only: output_int_detailedprogress
  use timer_module
  use time_module
  use comm_module
  use msgtag_module
  use integralpar_module
  use iounitadmin_module

  implicit none


  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------


  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_shutdown_dipole:  begin")

  if (.not. comm_i_am_master()) then
     call stop_timer(timer_int_idle_2cob3c(integralpar_i_int_part))
     call stop_timer(timer_int_2cob3c(integralpar_i_int_part))
     if (output_int_detailedprogress) call write_to_output_units( &
          "integral_shutdown_dipole: reporting back to master")
     call comm_init_send(comm_master_host,msgtag_int_dipole_rep_back) 
     call comm_send()  
  endif

  if (output_int_detailedprogress) call write_to_output_units( &
       "integral_shutdown_dipole:  end")


end subroutine integral_shutdown_dipole

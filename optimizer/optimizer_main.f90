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
program optimizer_main
  !----------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Subroutine called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: ...
  !  Date: ...
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

  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use optimizer, only: main_opt
  use opt_data_module, only: filename_setup_opt

  implicit none

  !== Interrupt end of public interface of module =================
  !------------ Declaration of formal parameters ------------------

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of subroutines ------------------------

  !------------ Declaration of local constants --------------------

  !------------ Declaration of local variables --------------------
  logical :: converged


  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------


  print *, 'optimizer_main: entered'
  print *, 'optimizer_main: call filename_setup_opt()'
  call filename_setup_opt(abortif=.false.,optonly=.true.)
  print *, 'optimizer_main: call main_opt(GeoOpt, converged)'
  call main_opt("GeoOpt", converged)
  ! The last line should be kept like this if there is no reason
  ! to change. At least in one case an external wrapper uses this
  ! to decide if geometry is converged:
  print *, 'optimizer_main: converged=',converged
end program optimizer_main

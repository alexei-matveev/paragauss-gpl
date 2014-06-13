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
subroutine build_mol_surfaces()
  !----------------------------------------------------------------
  !
  !  Purpose: build molecular surfaces for different solvation tasks
  !           and calculate cavitation and disp-rep energies
  !
  !  Subroutine called by: main_master
  !
  !
  !  References: ...
  ! 
  !
  !  Author: AS
  !  Date: 12.07.2007
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
  use solv_cavity_module
  use solv_electrostat_module

  implicit none

  !== Interrupt end of public interface of module =================
  !------------ Declaration of formal parameters ------------------
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of subroutines ------------------------
  !------------ Declaration of local constants --------------------
  !------------ Declaration of local variables --------------------
  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------

  do_gradients=.false.
  if(disp_rep_energy) then
     do_cavitation=.false.
     do_disp_rep=.true.
     call disp_rep_wrap()
  endif
  if(cavitation_energy) then
     do_cavitation=.true.
     do_disp_rep=.false.
     call points_on_cavity_surface()
     call energy_and_grad_of_cavity()
  endif
  call correction_param()
  do_cavitation=.false.
  do_disp_rep=.false.
  call points_on_cavity_surface()
  call matrix_generation()
  call solv_poten_transfer_data()

end subroutine build_mol_surfaces

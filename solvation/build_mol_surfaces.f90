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
subroutine build_mol_surfaces()
  !----------------------------------------------------------------
  !
  ! Build  molecular  surfaces   for  different  solvation  tasks  and
  ! calculate cavitation and disp-rep energies.
  !
  ! Subroutine called by: main_master() and runs on all workers.
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
# include "def.h"
  use solv_cavity_module, only: cavitation_energy, &
       disp_rep_energy, do_cavitation, do_disp_rep, do_gradients, &
       points_on_cavity_surface, correction_param
  use solv_electrostat_module, only: matrix_generation, &
       solv_poten_transfer_data
  use paragauss, only: toggle_legacy_mode
  use comm, only: comm_same
  implicit none
  ! *** end of interface ***

  do_gradients = .false.
  ASSERT(comm_same(disp_rep_energy))
  if (disp_rep_energy) then
     do_cavitation = .false.
     do_disp_rep = .true.
     do while (toggle_legacy_mode())
        call disp_rep_wrap()
     enddo
  endif

  ASSERT(comm_same(cavitation_energy))
  if (cavitation_energy) then
     do_cavitation = .true.
     do_disp_rep = .false.
     do while (toggle_legacy_mode())
        call points_on_cavity_surface()
        call energy_and_grad_of_cavity()
     enddo
  endif
  do while (toggle_legacy_mode())
     call correction_param()
  enddo
  do_cavitation = .false.
  do_disp_rep = .false.

  do while (toggle_legacy_mode())
     call points_on_cavity_surface()
     call matrix_generation()
     call solv_poten_transfer_data()
  enddo
end subroutine build_mol_surfaces

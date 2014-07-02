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
  ! calculate cavitation and disp-rep energies.  Subroutine called by:
  ! main_master()  and runs  on all  workers. Though  the bulk  of the
  ! actual  work is  done on  master  only. Does  not seem  to do  any
  ! communication.
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
  use solv_cavity_module, only: cavitation_energy, &
       disp_rep_energy, do_cavitation, do_disp_rep, do_gradients, &
       points_on_cavity_surface, correction_param
  use solv_electrostat_module, only: matrix_generation, &
       solv_poten_transfer_data
  use comm, only: comm_rank
  implicit none
  ! *** end of interface ***

  integer :: rank

  rank = comm_rank()

  ! Flags used in the conditionals seem to have the same values on all
  ! workers.
  do_gradients = .false.
  if (disp_rep_energy) then
     do_cavitation = .false.
     do_disp_rep = .true.
     if (rank == 0) then
        call disp_rep_wrap()    ! no comm?
     endif
  endif

  if (cavitation_energy) then
     do_cavitation = .true.
     do_disp_rep = .false.
     if (rank == 0) then
        call points_on_cavity_surface()  ! no comm?
        call energy_and_grad_of_cavity() ! no comm
     endif
  endif
  ! This one  sets a few  parameters, I think  there is no harm  to do
  ! that on all workers:
  call correction_param()       ! no comm
  do_cavitation = .false.
  do_disp_rep = .false.

  if (rank == 0) then
     call points_on_cavity_surface() ! no comm?
     call matrix_generation()        ! no comm
     call solv_poten_transfer_data() ! no comm
  endif
end subroutine build_mol_surfaces

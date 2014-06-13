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
subroutine grad_solv_calculate(task)
!
!  This is the warp procedure to calculate
!  matrix elements (integrals) of gradients of electronic part of 
!  electrostatic term of the free solvation energy
!
!  Author: AS
!  Date: 5/00
!
!================================================================
! End of public interface of module
!================================================================
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

  use type_module          ! contains standard data types
  use iounitadmin_module
  use output_module        ! defines amount of output
  use integralpar_module, only: integralpar_set
  use solv_cavity_module, only: send_receive_geom_grad
  use comm_module
  use gradient_data_module
  use interfaces, only: main_integral
  use integralpar_module, only: integralpar_cpks_contribs
#ifdef WITH_EFP
  use unique_atom_module, only: N_unique_atoms
#endif
  implicit none

  character(len=*), intent(in) :: task
  logical :: Q_deriv

  ! *** end of interface ***

  Q_deriv=.false.
  if(trim(task) == 'Q_SolvGrads') Q_deriv=.true.

  ! this sub now runs in a parallel context!

  if(.not.integralpar_cpks_contribs .and. .not.Q_deriv) then
     call send_receive_geom_grad(gradient_index)
#ifdef WITH_EFP
     if (N_unique_atoms == 0) return
#endif
  end if
 
  if( output_int_progress ) call write_to_output_units &
       ("grad_solv_calculate: Starting the Integral Part")
  if( output_int_progress ) call write_to_output_units &
       ("grad_solv_calculate: calling grad_solv_calculate ")
  call write_to_trace_unit("grad_solv: calling grad_solv_calculate ")
  call integralpar_set(task)

  call main_integral ()

  if( output_int_progress ) call write_to_output_units &
       ("grad_solv: done with the Integral Part")

  call write_to_trace_unit("grad_solv: finising grad_solv_calculate ")

end subroutine grad_solv_calculate

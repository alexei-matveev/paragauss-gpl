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
subroutine potential_calculate(task)
!
!  This is the warp procedure to calculate
!  matrix elements (integrals) of electronic part
!  of electrostatic potential in some special
!  points (e.g. surface points of a cavity)
!
!  Author: AS
!  Date: 11/99
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
  use integralpar_module, only: integralpar_pot_for_secderiv
  use potential_module, only: bounds_calc_poten, bounds_send_poten
  use comm_module, only: comm_parallel
  use interfaces, only: main_integral, IMAST,IPARA
  implicit none
  character(len=*), intent(in) :: task ! 'Potential', what else?
  ! *** end of interface ***

  if( output_int_progress ) call write_to_output_units &
       ("potential_calculate: Starting the Integral Part")
  if( output_int_progress ) call write_to_output_units &
       ("potential_calculate: calling potential_calculate ")
  call write_to_output_units("potential_calculate: calling potential_calculate ")

  ! FIXME: I dont really get trough the logic, FIX THIS:
  select case(task)
  case ('Solvation','Potential')
    ! is Solvation the same as 'Potential?', there are three calls in main_master!
    call integralpar_set('Potential')
  case default
    call integralpar_set(task) ! also sets integralpar_pot_for_secderiv
  end select

  ! FIXME: does this sub run in master-only context?
  !    AM: NO! it is called from main_master AND from parallel
  !        context in main_gradient
  if(integralpar_pot_for_secderiv) then
     call main_integral(IPARA)
  else
     call main_integral(IMAST)
  end if

  if( output_int_progress ) call write_to_output_units &
       ("potential_calculate: done with the Integral Part")

  if(.not.integralpar_pot_for_secderiv) then
     call bounds_calc_poten
     if (comm_parallel()) call bounds_send_poten
  end if

  call write_to_output_units("potential_calculate: finising potential_calculate ")
end subroutine potential_calculate

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
subroutine field_calculate
!
!  This is the warp procedure to calculate matrix elements (integrals)
!  of electronic  part of electrostatic  field at some  special points
!  (e.g. surface points of a solute cavity)
!
!  Runs in parallel context.
!
!  Author: AS
!  Date: 03/00
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
  use elec_static_field_module, only: bounds_calc_field
  use comm_module, only: comm_parallel
  use interfaces, only: main_integral

  implicit none

  if (output_int_progress) call write_to_output_units &
       ("field_calculate: Starting the Integral Part")
  if (output_int_progress) call write_to_output_units &
       ("field_calculate: calling field_calculate ")
  call write_to_trace_unit ("field_calculate: calling field_calculate")
  call integralpar_set ('Field')

  call main_integral ()

  if (output_int_progress) call write_to_output_units &
       ("field_calculate: done with the Integral Part")

  call bounds_calc_field ()

  call write_to_trace_unit ("field_calculate: finising field_calculate")
end subroutine field_calculate


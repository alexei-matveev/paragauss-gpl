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
logical function  precision_check(org,new)
!----------------------------------------------------------------
!
!  Purpose: check if relative difference between org and new
!           is below machineperameters_checkpreci from
!           machineparameters_module (an input parameter).
!
!  Subroutine called by: routines in integral part when
!           operations_integraltest is selected
!
!  Author: TB
!  Date:   10/96
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
  use machineparameters_module, only: machineparameters_checkpreci
  implicit none

!== Interrupt end of public interface of module =================
  !------------ Declaration of formal parameters ------------------
  real(kind=r8_kind), intent(in) :: org, new
!================================================================
! End of public interface of module
!================================================================

  if ( abs(org) .gt. machineparameters_checkpreci ) then
     precision_check = abs((org-new)/org) .gt. machineparameters_checkpreci
  else
     precision_check = abs(org-new) .gt. machineparameters_checkpreci
  endif

end function precision_check

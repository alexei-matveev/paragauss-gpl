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
module type_module
!
!  Purpose: Spezification of machine independent type parameters
!           for use in variable declarations and constants
!
!  Example: subroutine xx()
!           use type_module
!           real (kind=r8_kind) :: x
!           integer(kind=i4_kind) :: i
!           x = 1.0E1_r8_kind
!           i = int(x,i4_kind) ! intrinsic function
!           end subroutine xx
!            
!
!  Called by: everything
!
!  Author: TB
!  Date: 06/27/95
!
!== Interrupt of public interface of module ====================
!
!---------------------------------------------------------------
! Modifications
!---------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!--------------------------------------------------------------
!== Interrupt end of public interface of module ===============

implicit none
save

! real with 8 bytes, precision 15 decimal digits
integer, parameter :: double_precision_kind = selected_real_kind(15)
integer, parameter :: r8_kind = double_precision_kind

! real with 4 bytes, precision 6 decimal digits
integer, parameter :: single_precision_kind = selected_real_kind(6)
integer, parameter :: r4_kind = single_precision_kind

! complex with 16 bytes, precision 15 decimal digits
integer, parameter :: c16_kind = r8_kind

! complex with 8 bytes, precision 6 decimal digits
integer, parameter :: c8_kind = r4_kind

! integer with 8 bytes, range 9 decimal digits
integer, parameter :: i8_kind = 8

! integer with 4 bytes, range 9 decimal digits
integer, parameter :: integer_kind = selected_int_kind(9)
integer, parameter :: i4_kind = integer_kind

! integer with 2 bytes, range 4 decimal digits
integer, parameter :: i2_kind = selected_int_kind(4)

! integer with 1 byte, range 2 decimal digits
integer, parameter :: i1_kind = selected_int_kind(2)

  !===================================================================
  ! End of public interface of module
  !===================================================================

end module type_module

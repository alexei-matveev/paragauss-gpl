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
module <name>
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
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

  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------
  type, public ::  <name>_
  end type <name>_

  !------------ Declaration of constants and variables ------------
  integer(i4_kind), parameter, public :: <name>_
  real(r8_kind), parameter, public :: <name>_
  logical, parameter, public :: <name>_
  character, parameter, public :: <name>_
  integer(i4_kind), public :: <name>_
  real(r8_kind), public :: <name>_
  logical, public :: <name>_
  character, public :: <name>_


  !------------ Interface statements ------------------------------
  interface <name>_
  end interface<name>_
  public :: <name>_

  !------------ public functions and subroutines ------------------
  public :: <name>_

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------
  type
  end type

  !------------ Declaration of constants and variables ----
  integer(i4_kind), parameter ::
  real(r8_kind), parameter ::
  logical, parameter ::
  character(len=...), parameter ::
  integer(i4_kind) ::
  real(r8_kind) ::
  logical ::
  character(len=...) ::


  !------------ Subroutines ---------------------------------------
contains

  subroutine <name>_
    !
    ! Purpose ...
    !
    use
    implicit none
    integer(i4_kind), intent( ) ::
    real(r8_kind), intent( ) ::
    logical, intent( ) ::
    character, intent( ) ::
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    integer(i4_kind) ::
    real(r8_kind) ::
    logical ::
    character(len=...) ::

    !------------ Executable code --------------------------------


  end subroutine <name>_

  !--------------- End of module ----------------------------------
end module <name>

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
module matrix_xpack
  !-------------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  ! TO BE DONE:
  ! 1) ...
  ! 2) ...

# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  use matrix_types
  implicit none
  save

  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Interface statements ---------------------------------

  interface cpack
     module procedure pack_cmatrix
     module procedure pack_rmatrix
  end interface

  interface cunpack
     module procedure unpack_cmatrix
     module procedure unpack_rmatrix
  end interface

  !------------ public functions and subroutines ---------------------
  public &
       & cpack,&
       & cunpack

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables -----


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  ! FIXME: so far for square only

  subroutine pack_cmatrix(a)
    use xpack
    implicit none
    type(cmatrix),intent(in) :: a
    ! *** end of interface ***

    integer(IK) :: n(2)

    n = shape(a%re)

    call pck(n)
    call pck(a%re)
    call pck(a%im)
  end subroutine pack_cmatrix

  subroutine unpack_cmatrix(a)
    use xpack
    use matrix_methods, only: alloc
    implicit none
    type(cmatrix),intent(inout) :: a
    ! *** end of interface ***

    integer(IK) :: n(2)

    call upck(n)
    call alloc(n(1), n(2), a)
    call upck(a%re)
    call upck(a%im)
  end subroutine unpack_cmatrix

  subroutine pack_rmatrix(a)
    use xpack
    implicit none
    type(rmatrix),intent(in) :: a
    ! *** end of interface ***

    integer(IK) :: n(2)

    n = shape(a%m)

    call pck(n)
    call pck(a%m)
  end subroutine pack_rmatrix

  subroutine unpack_rmatrix(a)
    use xpack
    use matrix_methods, only: alloc
    implicit none
    type(rmatrix),intent(inout) :: a
    ! *** end of interface ***

    integer(IK) :: n(2)

    call upck(n)
    call alloc(n(1), n(2), a)
    call upck(a%m)
  end subroutine unpack_rmatrix

  !--------------- End of module -------------------------------------
end module matrix_xpack

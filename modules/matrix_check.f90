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
module matrix_check
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
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
#include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind,&
       & CK=>c16_kind ! type specification parameters
  use matrix_types
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  interface square
     module procedure square_r
     module procedure square_c

     module procedure square_2
  end interface

  interface conform
     module procedure conform_c_c
     module procedure conform_r_r
     module procedure conform_d_d

     module procedure conform_1_1
     module procedure conform_2_2
  end interface

  interface mconform
     module procedure mconform_r_r
     module procedure mconform_c_c
     module procedure mconform_d_c
     module procedure mconform_c_d
     module procedure mconform_r_c
     module procedure mconform_c_r

     module procedure mconform_2_2
     module procedure mconform_2_1
     module procedure mconform_1_2
     module procedure mconform_i1_i1_t
  end interface

  !------------ public functions and subroutines ------------------
  public :: square, conform, mconform ! is a matrix <...> ?

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !
  ! CHECK EVERYTHING:
  !

  function square_c(a) result(yes)
    implicit none
    type(cmatrix), intent(in) :: a
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(a%re,a%im))

    yes = square(a%re)
  end function square_c

  function square_r(a) result(yes)
    implicit none
    type(rmatrix), intent(in) :: a
    logical :: yes
    ! *** end of interface ***

    yes = square(a%m)
  end function square_r

  function square_2(a) result(yes)
    implicit none
    real(RK), intent(in) :: a(:, :)
    logical :: yes
    ! *** end of interface ***

    yes = size(a, 1) .eq. size(a, 2)
  end function square_2

  !
  ! are two matrices conform for summation (multiplication)? >>>
  !

  function conform_c_c(a, b) result(yes)
    implicit none
    type(cmatrix), intent(in) :: a, b
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(a%re,a%im))
    ASSERT(conform(b%re,b%im))

    yes = conform(a%re, b%re)
  end function conform_c_c

  function conform_r_r(a, b) result(yes)
    implicit none
    type(rmatrix), intent(in) :: a, b
    logical :: yes
    ! *** end of interface ***

    yes = conform(a%m, b%m)
  end function conform_r_r

  function conform_d_d(a, b) result(yes)
    implicit none
    type(rdmatrix), intent(in) :: a, b
    logical :: yes
    ! *** end of interface ***

    yes = conform(a%d, b%d)
  end function conform_d_d

  function conform_1_1(a, b) result(yes)
    implicit none
    real(RK), intent(in) :: a(:), b(:)
    logical :: yes
    ! *** end of interface ***

    yes = size(a).eq.size(b)
  end function conform_1_1

  function conform_2_2(a, b) result(yes)
    implicit none
    real(RK), intent(in) :: a(:,:), b(:,:)
    logical :: yes
    ! *** end of interface ***

    yes = all(shape(a).eq.shape(b))
  end function conform_2_2

  function mconform_c_c(a, b, trans) result(yes)
    implicit none
    type(cmatrix), intent(in) :: a, b
    character(len=2), intent(in), optional :: trans
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(a%re,a%im))
    ASSERT(conform(b%re,b%im))

    yes = mconform(a%re, b%re, trans)
  end function mconform_c_c

  function mconform_r_r(a, b, trans) result(yes)
    implicit none
    type(rmatrix), intent(in) :: a, b
    character(len=2), intent(in), optional :: trans
    logical :: yes
    ! *** end of interface ***

    yes = mconform(a%m, b%m, trans)
  end function mconform_r_r

  function mconform_2_2(a, b, trans) result(yes)
    implicit none
    real(RK), intent(in) :: a(:, :), b(:, :)
    character(len=2), intent(in), optional :: trans
    logical :: yes
    ! *** end of interface ***

    yes = mconform(shape(a), shape(b), trans)
  end function mconform_2_2

  function mconform_i1_i1_t(a, b, trans) result(yes)
    implicit none
    integer(IK), intent(in) :: a(:), b(:) ! shapes
    character(len=2), intent(in), optional :: trans
    logical :: yes
    ! *** end of interface ***

    integer(IK) :: a2, b1
    character(len=2) :: trans_

    ASSERT(size(a)==2)
    ASSERT(size(b)==2)

    if(present(trans))then
       trans_ = trans
    else
       trans_ = "nn"
    endif

    if(trans_(1:1).eq."n")then
       a2 = a(2)
    else
       a2 = a(1)
    endif
    if(trans_(2:2).eq."n")then
       b1 = b(1)
    else
       b1 = b(2)
    endif

    yes = (a2.eq.b1)
  end function mconform_i1_i1_t

  function mconform_d_c(a, b) result(yes)
    implicit none
    type(rdmatrix), intent(in) :: a
    type(cmatrix), intent(in) :: b
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(b%re,b%im))

    yes = mconform(a%d, b%re)
  end function mconform_d_c

  function mconform_c_d(a, b) result(yes)
    implicit none
    type(cmatrix), intent(in) :: a
    type(rdmatrix), intent(in) :: b
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(a%re,a%im))

    yes = mconform(a%re, b%d)
  end function mconform_c_d

  function mconform_1_2(a, b) result(yes)
    implicit none
    real(RK), intent(in) :: a(:)
    real(RK), intent(in) :: b(:,:)
    logical :: yes
    ! *** end of interface ***

    yes = size(a) == size(b, 1)
  end function mconform_1_2

  function mconform_2_1(a, b) result(yes)
    implicit none
    real(RK), intent(in) :: a(:,:)
    real(RK), intent(in) :: b(:)
    logical :: yes
    ! *** end of interface ***

    yes = size(a, 2) == size(b)
  end function mconform_2_1

  function mconform_r_c(a, b) result(yes)
    implicit none
    type(rmatrix), intent(in) :: a
    type(cmatrix), intent(in) :: b
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(b%re,b%im))

    yes = mconform(a%m, b%re)
  end function mconform_r_c

  function mconform_c_r(a, b) result(yes)
    implicit none
    type(cmatrix), intent(in) :: a
    type(rmatrix), intent(in) :: b
    logical :: yes
    ! *** end of interface ***

    ASSERT(conform(a%re,a%im))

    yes = mconform(a%re, b%m)
  end function mconform_c_r

  !--------------- End of module ----------------------------------
end module matrix_check

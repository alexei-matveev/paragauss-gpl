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
module matrix_matmult
  !-------------------------------------------------------------------
  !
  !  Purpose: We have a fundamental problem, intrinsic MATMUL
  !           having access to array descriptors may be able to
  !           handle non-continous sections and transposition
  !           properly (e.g. without creating continous temporaries).
  !           On the other hand it may be less efficient than
  !           platform optimized BLAS. Also in the past some
  !           MATMUL versions didnt handle matrix * vector
  !           multiplicaitons properly which actualy was one of
  !           the reasons of creating the MATMULT fallback interface.
  !
  !           Which one, MATMUL or MATMULT, to use may thus depend
  !           on situation.
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

# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind,&
       & CK=>c16_kind ! type specification parameters
  use matrix_types, only: verbose
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  real(RK), parameter :: one  = 1.0_rk
  real(RK), parameter :: zero = 0.0_rk

  !------------ Interface statements ---------------------------------

  interface matmult
     !
     ! try to replace every instance 
     ! of MATMUL with MATMULT
     !
     module procedure matmult_2_2
     module procedure matmult_2_2_t
     module procedure matmult_dgemm
     module procedure matmult_row_rr
     module procedure matmult_col_rr
     module procedure matmult_col_cc
  end interface


  !------------ public functions and subroutines ---------------------
  public matmult

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  function matmult_2_2(A, B) result(C)
    !
    ! Plain matrix multiplicaiton
    !
    implicit none
    real(RK), intent(in) :: A(:, :), B(:, :)
    real(RK) :: C(size(A, 1), size(B, 2))
    !** End of interface ***

    ASSERT(size(A,2)==size(B,1))

    ! FIXME: is intrinsic MATMUL good enough?
    C = matmult(size(C, 1), size(C, 2), size(A, 2), A, B, "nn")
  end function matmult_2_2

  function matmult_2_2_t(A, B, trans) result(C)
    !
    ! USE WITH CARE, IT CANNOT BE GENERALIZED TO ARBITRARY SHAPES!
    !
    ! Matrix multiplicaiton, pass transposiiton info separately
    !
    ! FIXME: we cannot adapt size of C depending on
    !        transa, transb. Assume a matrix that is to
    !        be transposed is square.
    !
    !        Use more flexible matmult_dgemm(m, n, k, ...)
    !        directly.
    !
    implicit none
    real(RK), intent(in) :: A(:, :), B(:, :)
    character(len=2), intent(in) :: trans
    real(RK) :: C(size(A, 1), size(B, 2))
    !** End of interface ***

    integer(IK) :: m, n, k

    if ( trans(1:1) /= 'n' ) then
        ASSERT(trans(1:1)=='t')
        ! see size(A, 1) in decl of C:
        ASSERT(size(A,1)==size(A,2))
    endif

    if ( trans(2:2) /= 'n' ) then
        ASSERT(trans(2:2)=='t')
        ! see size(B, 2) in decl of C:
        ASSERT(size(B,1)==size(B,2))
    endif

    ! just as in decl of C:
    m = size(A, 1)
    n = size(B, 2)

    ! if B is transposed, then it is square anyway:
    k = size(B, 1) ! or size(A, 2)

    ! FIXME: is intrinsic MATMUL good enough?
    C = matmult(m, n, k, A, B, trans)
  end function matmult_2_2_t

  function matmult_dgemm(m, n, k, A, B, trans) result(C)
    use f77_blas, only: dgemm
    implicit none
    integer(IK), intent(in) :: m, n, k
    real(RK), intent(in) :: A(:, :), B(:, :)
    character(len=2), intent(in) :: trans
    real(RK) :: C(m, n)
    !** End of interface ***

    if(verbose)print *, 'mm/matmult_dgemm: entered'

    if ( trans(1:1) == "n" ) then
        ASSERT(m<=size(A,1))
        ASSERT(k<=size(A,2))
    else
        ASSERT(trans(1:1)=='t')
        ASSERT(m<=size(A,2))
        ASSERT(k<=size(A,1))
    endif

    if ( trans(2:2) == "n" ) then
        ASSERT(n<=size(B,2))
        ASSERT(k<=size(B,1))
    else
        ASSERT(trans(2:2)=='t')
        ASSERT(n<=size(B,1))
        ASSERT(k<=size(B,2))
    endif

    !
    ! FIXME: more "physical" matrix multiplication than purely mathematical
    !
    if ( m == 0 ) then
       if ( verbose ) print *, 'mm/matmult_dgemm: case 0xKxN'
       return
    endif
    if ( n == 0 ) then
       if ( verbose ) print *, 'mm/matmult_dgemm: case MxKx0'
       return
    endif
    if ( k == 0 ) then
       C = zero
       if ( verbose ) print *, 'mm/matmult_dgemm: case Mx0xN'
       return
    endif

    call dgemm(trans(1:1), trans(2:2), m, n, k, one, &
        A, size(A, 1), &
        B, size(B, 1), zero, &
        C, size(C, 1))
  end function matmult_dgemm

  !
  ! These were meant as a workaround for the buggy intrinsic
  ! MATMUL not handling properly matrix-vector multiplications.
  !

  function matmult_row_rr(s, m) result(c)
    !
    ! Real row vector * real matrix
    !
    implicit none
    real(RK), intent(in) :: s(:), m(:, :)
    real(RK) :: c(size(m, 2))
    ! *** end of interface ***

    integer(IK) :: j, k

    ASSERT(size(s)==size(m,1))
    do k = 1, size(m, 2)
       c(k) = 0.0_RK
       do j = 1, size(s)
          c(k) = c(k) + s(j) * m(j, k)
       enddo
    enddo
  end function matmult_row_rr

  function matmult_col_rr(m, c) result(s)
    !
    ! Real matrix * real column vector
    !
    implicit none
    real(RK), intent(in) :: c(:), m(:,:)
    real(RK) :: s(size(m, 1))
    ! *** end of interface ***

    integer(IK) :: i, j

    ASSERT(size(c)==size(m,2))
    do i = 1, size(m, 1)
       s(i) = 0.0_RK
       do j = 1, size(c)
          s(i) = s(i) + m(i, j) * c(j)
       enddo
    enddo
  end function matmult_col_rr

  function matmult_col_cc(m, c) result(s)
    !
    ! Complex matrix * complex column vector
    !
    implicit none
    complex(CK), intent(in) :: c(:), m(:,:)
    complex(CK) :: s(size(m, 1))
    ! *** end of interface ***

    integer(IK) :: i, j

    ASSERT(size(c)==size(m,2))
    do i = 1, size(m, 1)
       s(i) = 0.0_CK
       do j = 1, size(c)
          s(i) = s(i) + m(i, j) * c(j)
       enddo
    enddo
  end function matmult_col_cc

  !--------------- End of module -------------------------------------
end module matrix_matmult





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
module matrix_methods
  !-------------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Names are mangled by suffixes indicating the types
  !  of accepted arguments:
  !
  !     0       double precision scalar
  !     1       double precision array
  !             interpreted as diagonal matrix
  !     2       double precision array with two axes
  !             interpreted as a rectangular matrix
  !     r       type(rmatrix)
  !     c       type(cmatrix)
  !     d       type(rdmatrix)
  !     h       type(chmatrix)
  !     i       integer
  !
  !  Return value, when indicated, follows the double underscore.
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
  ! 1) remove check for hermiteanity(?) from pck- and from_h__c
  ! 2) sim() can fail for non-square matrices
  ! 3) rdmatrix ar[iy]thmetics ?
  ! 4) ...
# include "def.h"
  use type_module, only: IK => i4_kind, RK => r8_kind ! kinds
  use matrix_types, only: rmatrix, cmatrix, &
       rdmatrix, chmatrix, verbose
  use matrix_check
  use matrix_matmult, only: matmult
  implicit none

  save
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Interface statements ---------------------------------

  interface matrix
     !
     ! Convert array into dense- or diagonal matrix:
     !
     module procedure matrix_2
     module procedure matrix_1
  end interface

  public :: matrix

  interface array
     !
     ! Convert a matrix to a plain array:
     !
     module procedure array_r
     module procedure array_d
  end interface

  public :: array

  interface assignment(=)
     ! rdmatrix = array(:)
     module procedure assign_r_1
  end interface

  interface convert
     ! pack and unpack complex hermitean mtrx:
     module procedure from_c__h
     module procedure from_h__c
  end interface

  interface diag
     ! extracting diagonal
     module procedure diag_r
  end interface

  interface operator(*)
     ! matrix * number:
     module procedure mult_0_r
     module procedure mult_r_0
     module procedure mult_0_c
     module procedure mult_c_0

     ! matrix * matrix:
     module procedure mult_r_r
     module procedure mult_c_c

     ! matrix * diagonal:
     module procedure mult_r_d
     module procedure mult_d_r

     ! matrix * rdmatrix:
     module procedure mult_d_c
     module procedure mult_c_d

     ! diagonal matrix * diagonal matrix:
     module procedure mult_d_d
     module procedure mult_d_1 ! MUSTDIE!
     module procedure mult_1_d ! MUSTDIE!

     ! diagonal matrix * number:
     module procedure mult_d_0
     module procedure mult_0_d
  end interface

  interface operator(**)
     module procedure pow_d_i
  end interface

  interface operator(+)
     module procedure plus_r_r
     module procedure plus_c_c
     module procedure plus_d_d
  end interface

  interface operator(-)
     module procedure minus_r
     module procedure minus_r_r
     module procedure minus_c
     module procedure minus_c_c
     module procedure minus_d_d
  end interface

  interface tr
     !
     ! transpose ( hermitean conjugation )
     !
     module procedure tr_r
     module procedure tr_2
     module procedure tr_c
     module procedure tr_c_t
  end interface

  interface mult
     module procedure mult_0_r
     module procedure mult_r_0
     module procedure mult_0_c
     module procedure mult_c_0
     module procedure mult_r_r
     module procedure mult_c_c
     module procedure mult_r_d
     module procedure mult_d_r
     !
     module procedure mult_2_2
     module procedure mult_2_2_t
     module procedure mult_1_2
     module procedure mult_2_1
     !
     ! matrix * rdmatrix:
     module procedure mult_d_c
     module procedure mult_c_d
     !
     ! three factor products:
     module procedure mult3_d_r_d
     module procedure mult3_d_c_d
     module procedure mult3_1_2_1
     !
     ! diagonal matrix * diagonal matrix:
     module procedure mult_d_d
     module procedure mult_d_1
     module procedure mult_1_d
  end interface

  interface sim
     !
     ! similarity transformation U(+)* V * U
     ! function sim(V,U)
     !
     module procedure sim_r_r
     module procedure sim_c_c
     module procedure sim_2_2
  end interface

  interface bsim
     !
     ! [back] similarity transformation U * V * U(+)
     ! function bsim(U,V)
     !
     module procedure bsim_r_r
     module procedure bsim_c_c
  end interface

  interface alloc
     module procedure alloc_i_r
     module procedure alloc_i_i_r
     module procedure alloc_i_r_r_plus ! for two or more
     module procedure alloc_i_c
     module procedure alloc_i_i_c
     module procedure alloc_i_c_c_plus ! for two or more
     !
     ! rdmatrix:
     module procedure alloc_i_d
     module procedure alloc_i_d_d_plus ! for two or more
     !
     ! chmatrix:
     module procedure alloc_i_h
  end interface

  interface pgfree
     module procedure free_r
     module procedure free_r_r_plus ! for two or more
     module procedure free_c
     module procedure free_c_c_plus ! for two or more
     !
     ! rdmatrix:
     module procedure free_d
     module procedure free_d_d_plus
     !
     ! chmatrix:
     module procedure free_h
  end interface

  interface size
     module procedure size_d!(rdmatrix)
     module procedure size_r_i!(rmatrix, axis)
     module procedure size_c_i!(cmatrix, axis)
  end interface

  !------------ public functions and subroutines ---------------------
  public :: assignment(=)                               ! [cr]matrix = array
  public :: operator(*), mult                           ! multiplication
  public :: operator(**)                                ! exponentiation
  public :: operator(+), operator(-)                    ! ...
  public :: tr                                          ! TRANSPOSED, NOT TRACE !!
  public :: sim, bsim                                   ! similarity trafo
  public :: convert                                     ! c <-> h conversion
  public :: alloc, pgfree                               ! memory managment
  public :: diag                                        ! extracting diagonal
  public :: size                                        ! dimension

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables -----

  real(RK), parameter :: zero = 0.0_rk
  real(RK), parameter :: one = 1.0_rk

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

#define AT_ENTRY(proc)  ! if ( verbose ) print *, __STRING(proc), "entered"

  function matrix_2(array) result(matrix)
    implicit none
    real(RK), intent(in) :: array(:,:)
    type(rmatrix) :: matrix
    ! *** end of interface ***

    AT_ENTRY(matrix_2)

    matrix = rmatrix(array)
  end function matrix_2

  function array_r(matrix) result(array)
    implicit none
    type(rmatrix), intent(in) :: matrix
    real(RK) :: array(size(matrix, 1), size(matrix, 2))
    ! *** end of interface ***

    AT_ENTRY(array_r)

    array = matrix%m
  end function array_r

  function matrix_1(array) result(matrix)
    implicit none
    real(RK), intent(in) :: array(:)
    type(rdmatrix) :: matrix
    ! *** end of interface ***

    AT_ENTRY(matrix_1)

    matrix = rdmatrix(array)
  end function matrix_1

  function array_d(matrix) result(array)
    implicit none
    type(rdmatrix), intent(in) :: matrix
    real(RK) :: array(size(matrix))
    ! *** end of interface ***

    AT_ENTRY(array_d)

    array = matrix%d
  end function array_d

  subroutine assign_r_1(a, b)
    !
    ! FIXME:  some  legacy  code  contains  assignments  of  the  type
    ! rdmatrix = array.
    !
    implicit none
    type(rdmatrix), intent(out) :: a
    real(RK), intent(in) :: b(:)
    !** End of interface ***

    AT_ENTRY(assign_r_1)

    a = rdmatrix(b)
  end subroutine assign_r_1

  function tr_r(a) result(c)
    !
    ! Transpose of rmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a
    type(rmatrix) :: c
    ! *** end of interface **

    AT_ENTRY(tr_r)

#ifndef FPP_GFORTRAN_BUGS /* 4.6 produces NaNs here */
    ! type constructor here:
    c = rmatrix(transpose(a%m))
#else
    call alloc(size(a, 2), size(a, 1), c)
    c%m = transpose(a%m)
#endif
  end function tr_r

  function tr_2(a) result(c)
    implicit none
    real(RK),dimension(:,:),intent(in)      :: a
    real(RK),dimension(size(a,2),size(a,1)) :: c
    ! *** end of interface **

    AT_ENTRY(tr_2)

    c = transpose(a)
  end function tr_2

  function tr_c(a) result(c)
    implicit none
    type(cmatrix), intent(in) :: a
    type(cmatrix) :: c
    ! *** end of interface **

    AT_ENTRY(tr_c)

#ifndef FPP_GFORTRAN_BUGS /* 4.6 produces NaNs here */
    ! type constructor here:
    c = cmatrix(transpose(a%re), -transpose(a%im))
#else
    call alloc(size(a, 2), size(a, 1), c)
    c%re = transpose(a%re)
    c%im = -transpose(a%im)
#endif
  end function tr_c

  function tr_c_t(a, op) result(c)
    implicit none
    type(cmatrix), intent(in) :: a
    character(len=1), intent(in) :: op
    type(cmatrix) :: c
    ! *** end of interface **

    AT_ENTRY(tr_c_t)

    select case(op)
    case ("h")
       c = tr(a)
       WARN("call tr(A) instead")
    case ("t")
       ! type constructor here:
       c = cmatrix(transpose(a%re), transpose(a%im))
    case default
       ! to make compiler happy:
       c = cmatrix(NULL(), NULL())
       ABORT("mm/tr_c_t: no such case")
    end select
  end function tr_c_t

  function sim_2_2(v, u) result(w)
    ! returns U(+) * V * U
    implicit none
    real(RK), intent(in) :: v(:,:) ! (n, n)
    real(RK), intent(in) :: u(:,:) ! (n, m)
    real(RK), dimension(size(u,2), size(u,2)) :: w ! (m, m)
    ! *** end of interface **

    integer(IK) :: m, n

    AT_ENTRY(sim_2_2)

    m = size(u, 2)
    n = size(u, 1)

    ASSERT(n==size(v,1))
    ASSERT(n==size(v,2))

    ! FIXME: need a rectangular temporary otherwise:
    ASSERT(m==n)

    ! w  = tr(u) * v * u
    w = matmult(n, m, n, v, u, "nn")
    w = matmult(m, m, n, u, w, "tn")
  end function sim_2_2

  function sim_r_r(v, u) result(w)
    !
    ! Simialrity transformation:
    !
    !    T
    !   U  *  V  *  U
    !
    implicit none
    type(rmatrix), intent(in) :: v, u
    type(rmatrix) :: w
    ! *** end of interface **

    AT_ENTRY(sim_r_r)

    ASSERT(mconform(v,u))
    ASSERT(mconform(u,v,'tn'))

    w  = tr(u) * v * u
  end function sim_r_r

  function sim_c_c(v, u) result(w)
    implicit none
    type(cmatrix), intent(in) :: v, u
    type(cmatrix) :: w
    ! *** end of interface **

    AT_ENTRY(sim_c_c)

    ASSERT(mconform(v,u))
    ASSERT(mconform(u,v,'tn'))

    w  = tr(u) * v * u
  end function sim_c_c

  function bsim_r_r(u, v) result(w)
    !
    ! Backwards simialrity transformation:
    !
    !                T
    !   U  *  V  *  U
    !
    implicit none
    type(rmatrix), intent(in) :: v, u
    type(rmatrix) :: w
    ! *** end of interface **

    AT_ENTRY(bsim_r_r)

    ASSERT(mconform(u,v))
    ASSERT(mconform(v,u,'nt'))

    w  = u * (v * tr(u))
  end function bsim_r_r

  function bsim_c_c(u, v) result(w)
    implicit none
    type(cmatrix), intent(in) :: v, u
    type(cmatrix) :: w
    ! *** end of interface **

    AT_ENTRY(bsim_c_c)

    ASSERT(mconform(u,v))
    ASSERT(mconform(v,u,'nt'))

    w  = u * (v * tr(u))
  end function bsim_c_c

  function mult3_1_2_1(d1, a, d2) result(c)
    implicit none
    real(RK), dimension(:,:), intent(in) :: a
    real(RK), dimension(:)  , intent(in) :: d1 ,d2
    real(RK), dimension(size(a,1), size(a,2)) :: c
    ! *** end of interface **

    integer(IK) :: n1,n2

    AT_ENTRY(mult3_1_2_1)

    n1 = size(a,1)
    n2 = size(a,2)

    ASSERT(n1==size(d1))
    ASSERT(n2==size(d2))

    c(:,:) = spread(d1,2,n2) * spread(d2,1,n1) ! <<< used as tmp storage !!!

    c(:,:) = a * c(:,:)
  end function mult3_1_2_1

  function mult3_d_r_d(d1, a, d2) result(c)
    !
    ! rdmatrix * rmatrix * rdmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a
    type(rdmatrix), intent(in) :: d1, d2
    type(rmatrix) :: c
    ! *** end of interface **

    AT_ENTRY(mult3_d_r_d)

    ! type constructor here:
    c = rmatrix(mult(d1%d, a%m, d2%d))
  end function mult3_d_r_d

  function mult3_d_c_d(d1, a, d2) result(c)
    implicit none
    type(cmatrix), intent(in) :: a
    type(rdmatrix), intent(in) :: d1, d2
    type(cmatrix) :: c
    ! *** end of interface **

    AT_ENTRY(mult3_d_c_d)

    ! type constructor here:
    c = cmatrix(mult(d1%d, a%re, d2%d), mult(d1%d, a%im, d2%d))
  end function mult3_d_c_d

  !
  ! PLUS - MINUS OPERATORS:
  !

  function plus_r_r(a, b) result(c)
    !
    ! rmatrix + rmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a, b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(plus_r_r)

    ASSERT(conform(a,b))

    ! type constructor here:
    c = rmatrix(a%m + b%m)
  end function plus_r_r

  function minus_r_r(a, b) result(c)
    !
    ! rmatrix - rmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a, b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(minus_r_r)

    ASSERT(conform(a,b))

    ! type constructor here:
    c = rmatrix(a%m - b%m)
  end function minus_r_r

  function minus_r(a) result(c)
    !
    ! - rmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(minus_r)

    ! type constructor here:
    c = rmatrix(-a%m)
  end function minus_r

  function plus_c_c(a, b) result(c)
    implicit none
    type(cmatrix), intent(in) :: a, b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(plus_c_c)

    ASSERT(conform(a,b))

    c = cmatrix(a%re + b%re, a%im + b%im)
  end function plus_c_c

  function minus_c_c(a, b) result(c)
    implicit none
    type(cmatrix), intent(in) :: a, b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(minus_c_c)

    ASSERT(conform(a,b))

    ! type constructor here:
    c = cmatrix(a%re - b%re, a%im - b%im)
  end function minus_c_c

  function minus_c(a) result(c)
    implicit none
    type(cmatrix), intent(in) :: a
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(minus_c)

    ! type constructor here:
    c = cmatrix(-a%re, -a%im)
  end function minus_c

  function plus_d_d(a, b) result(c)
    !
    ! rdmatrix + rdmatrix
    !
    implicit none
    type(rdmatrix), intent(in) :: a, b
    type(rdmatrix) :: c
    !** End of interface ***

    AT_ENTRY(plus_d_d)

    ASSERT(conform(a,b))

    ! type constructor here:
    c = rdmatrix(a%d + b%d)
  end function plus_d_d

  function minus_d_d(a, b) result(c)
    !
    ! rdmatrix - rdmatrix
    !
    implicit none
    type(rdmatrix), intent(in) :: a, b
    type(rdmatrix) :: c
    !** End of interface ***

    AT_ENTRY(minus_d_d)

    ASSERT(conform(a,b))

    ! type constructor here:
    c = rdmatrix(a%d - b%d)
  end function minus_d_d

  !
  ! MATRIX MULTIPLICATIONS:
  !

  !
  ! first multiplications with numbers:
  !

  function mult_0_r(a, b) result(c)
    !
    ! Number * rmatrix
    !
    implicit none
    real(RK), intent(in) :: a
    type(rmatrix), intent(in) :: b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_0_r)

    ! type constructor here:
    c = rmatrix(a * b%m)
  end function mult_0_r

  function mult_r_0(a, b) result(c)
    !
    ! rmatrix * number
    !
    implicit none
    type(rmatrix), intent(in) :: a
    real(RK), intent(in) :: b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_r_0)

    ! type constructor here:
    c = rmatrix(a%m * b)
  end function mult_r_0

  function mult_0_c(a,b) result(c)
    implicit none
    real(RK), intent(in) :: a
    type(cmatrix), intent(in) :: b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_0_c)

    ! type constructor here:
    c = cmatrix(a * b%re, a * b%im)
  end function mult_0_c

  function mult_c_0(a,b) result(c)
    implicit none
    type(cmatrix), intent(in) :: a
    real(RK), intent(in) :: b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_c_0)

    ! type constructor here:
    c = cmatrix(a%re * b, a%im * b)
  end function mult_c_0

  !
  ! now MATRIX * MATRIX:
  !

  function mult_r_r(a, b) result(c)
    !
    ! rmatrix * rmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a, b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_r_r)

    ASSERT(mconform(a,b))

    ! type constructor here:
    c = rmatrix(matmult(a%m, b%m))
  end function mult_r_r

  !
  ! plain matrix multiplication:
  !

  function mult_2_2(a, b) result(c)
    implicit none
    real(RK), dimension(:,:), intent(in) :: a, b
    real(RK), dimension(size(a,1), size(b,2)) :: c
    !** End of interface ***

    AT_ENTRY(mult_2_2)

    c = matmult(a, b)
  end function mult_2_2

  function mult_2_2_t(a, b, trans) result(c)
    implicit none
    real(RK), dimension(:,:), intent(in) :: a, b
    character(len=2), intent(in) :: trans ! "nn", "nt", etc.
    real(RK), dimension(size(a,1), size(b,2)) :: c
    !** End of interface ***

    AT_ENTRY(mult_2_2_t)

    c = matmult(a, b, trans)
  end function mult_2_2_t

  function mult_c_c(a,b) result(c)
    implicit none
    type(cmatrix), intent(in) :: a, b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_c_c)

    ASSERT(mconform(a,b))

    ! type constructor here:
    c = cmatrix(matmult(a%re, b%re) - matmult(a%im, b%im), &
                matmult(a%re, b%im) + matmult(a%im, b%re))
  end function mult_c_c

  !
  ! now MATRIX * DIAGONAL MATRIX:
  !

  function mult_1_2(a,b) result(c)
    implicit none
    real(RK), dimension(:), intent(in) :: a
    real(RK), dimension(:,:), intent(in) :: b
    real(RK), dimension(size(b,1),size(b,2)) :: c
    !** End of interface ***

    integer(IK) :: n1,n2,i

    AT_ENTRY(mult_1_2)

    ASSERT(mconform(a,b))

    n1 = size(b,1)
    n2 = size(b,2)
    do i=1,n1
       c(i,1:n2) = a(i) * b(i,1:n2)
    enddo
  end function mult_1_2

  function mult_2_1(a, b) result(c)
    implicit none
    real(RK), dimension(:,:), intent(in) :: a
    real(RK), dimension(:), intent(in) :: b
    real(RK), dimension(size(a,1),size(a,2)) :: c
    !** End of interface ***

    integer(IK) :: n1,n2,i

    AT_ENTRY(mult_2_1)

    ASSERT(mconform(a,b))

    n1 = size(a,1)
    n2 = size(a,2)
    do i=1,n2
       c(1:n1,i) = a(1:n1,i) * b(i)
    enddo
  end function mult_2_1

  function mult_r_d(a,b) result(c)
    !
    ! rmatrx * rdmatrix
    !
    implicit none
    type(rmatrix), intent(in) :: a
    type(rdmatrix), intent(in) :: b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_r_d)

    ! type constructor here:
    c = rmatrix(mult(a%m, b%d))
  end function mult_r_d

  function mult_d_r(a,b) result(c)
    !
    ! rmatrx * rdmatrix
    !
    implicit none
    type(rdmatrix), intent(in) :: a
    type(rmatrix), intent(in) :: b
    type(rmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_d_r)

    ! type constructor here:
    c = rmatrix(mult(a%d, b%m))
  end function mult_d_r

  function mult_d_c(a, b) result(c)
    implicit none
    type(rdmatrix), intent(in) :: a
    type(cmatrix), intent(in) :: b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_d_c)

    ! type constructor here:
    c = cmatrix(mult(a%d, b%re), mult(a%d, b%im))
  end function mult_d_c

  function mult_c_d(a, b) result(c)
    implicit none
    type(cmatrix), intent(in) :: a
    type(rdmatrix), intent(in) :: b
    type(cmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_c_d)

    ! type constructor here:
    c = cmatrix(mult(a%re, b%d), mult(a%im, b%d))
  end function mult_c_d

  !
  ! now DIAGONAL MATRIX * DIAGONAL MATRIX:
  !

  function mult_d_d(a,b) result(c)
    implicit none
    type(rdmatrix), intent(in) :: a
    type(rdmatrix), intent(in) :: b
    real(RK), dimension(size(a%d)) :: c
    !** End of interface ***

    AT_ENTRY(mult_d_d)
    ASSERT(size(a%d)==size(b%d))

    c = a%d * b%d
  end function mult_d_d

  function pow_d_i(a, p) result(ap)
    !
    !           p
    ! Computes A  for diagonal A.
    !
    implicit none
    type(rdmatrix), intent(in) :: a
    integer, intent(in) :: p
    type(rdmatrix) :: ap
    ! *** end of interface ***

    AT_ENTRY(pow_d_i)

    ! type constructor here:
    ap = rdmatrix(a%d**p)
  end function pow_d_i

  function mult_d_1(a,b) result(c)
    implicit none
    type(rdmatrix),intent(in)        :: a
    real(RK)      ,intent(in)        :: b(:)
    real(RK), dimension(size(b))     :: c !<<<result
    !** End of interface ***

    AT_ENTRY(mult_d_1)
    ASSERT(size(a%d)==size(b))

    c = a%d * b
  end function mult_d_1

  function mult_1_d(a,b) result(c)
    implicit none
    real(RK)      ,intent(in)        :: a(:)
    type(rdmatrix),intent(in)        :: b
    real(RK), dimension(size(a))     :: c !<<<result
    !** End of interface ***

    AT_ENTRY(mult_1_d)
    ASSERT(size(b%d)==size(a))

    c = a * b%d
  end function mult_1_d

  function mult_d_0(a, b) result(c)
    implicit none
    type(rdmatrix), intent(in) :: a
    real(RK), intent(in) :: b
    type(rdmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_d_0)

    ! type constructor here:
    c = rdmatrix(a%d * b)
  end function mult_d_0

  function mult_0_d(a, b) result(c)
    implicit none
    real(RK), intent(in) :: a
    type(rdmatrix), intent(in) :: b
    type(rdmatrix) :: c
    !** End of interface ***

    AT_ENTRY(mult_0_d)

    c = rdmatrix(a * b%d)
  end function mult_0_d

  !
  ! SOME ROUTINES FOR TYPE CHMATRIX:
  !

  function from_h__c(h, UPLO) result(cm)
    implicit none
    type(chmatrix), intent(in) :: h
    character, optional, intent(in) :: UPLO
    type(cmatrix) :: cm
    ! *** end of interface ***

    integer(IK) :: m, n, mn
    character   :: uplo_

    AT_ENTRY(from_h__c)

    uplo_ = 'U'
    if(present(UPLO)) uplo_ = UPLO

    call alloc(h%n, cm)

    select case(uplo_)
    case('U','u')
       mn = 0
       do n = 1, h%n
          do m = 1, n
             mn = mn + 1
             cm%re(m, n) =  h%re(mn)
             cm%im(m, n) =  h%im(mn)
             cm%re(n, m) =  h%re(mn)
             cm%im(n, m) = -h%im(mn)
          enddo
       enddo
    case('L','l')
       mn = 0
       do n = 1, h%n
          do m = 1, n
             mn = mn + 1
             cm%re(n, m) =  h%re(mn)
             cm%im(n, m) =  h%im(mn)
             cm%re(m, n) =  h%re(mn)
             cm%im(m, n) = -h%im(mn)
          enddo
       enddo
    case default
       ABORT('no such order: '//uplo_)
    end select
  end function from_h__c

  function from_c__h(c, UPLO) result(h)
    implicit none
    type(cmatrix), intent(in) :: c
    character, optional, intent(in) :: UPLO
    type(chmatrix) :: h
    ! *** end of interface ***

    integer(IK) :: m, n, mn
    character   :: uplo_

    AT_ENTRY(from_c__h)

    ASSERT(square(c))

    uplo_ = 'U'
    if(present(UPLO)) uplo_ = UPLO

    call alloc(size(c%re, 1), h)

    select case(uplo_)
    case('U','u')
       mn = 0
       do n = 1, h%n
          do m = 1, n
             mn = mn + 1
             h%re(mn) = c%re(m, n)
             h%im(mn) = c%im(m, n)
          enddo
       enddo
    case('L','l')
       mn = 0
       do n = 1, h%n
          do m = 1, n
             mn = mn + 1
             h%re(mn) = c%re(n, m)
             h%im(mn) = c%im(n, m)
          enddo
       enddo
    case default
       ABORT('no such order: '//uplo_)
    end select
  end function from_c__h

  !
  ! EXTRACTING DIAGONAL:
  !

  function diag_r(a) result (d)
    implicit none
    type(rmatrix), intent(in) :: a ! NxN
    real(RK), dimension(MIN(size(a%m,1), size(a%m,2))) :: d ! N
    ! *** end of interface ***

    integer(IK) :: i

    AT_ENTRY(diag_r)

    do i = 1, size(d)
       d(i) = a%m(i,i)
    enddo
  end function diag_r

  !
  ! MEMORY MANAGMENT PROCEDURES:
  !

  subroutine alloc_i_r(n, a)
    implicit none
    integer(IK),intent(in)      :: n
    type(rmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(alloc_i_r)

    allocate(a%m(n,n),STAT=memstat)
    ASSERT(memstat==0)
  end subroutine alloc_i_r

  subroutine alloc_i_i_r(n1, n2, a)
    implicit none
    integer(IK),intent(in)      :: n1,n2
    type(rmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(alloc_i_i_r)

    allocate(a%m(n1,n2),STAT=memstat)
    ASSERT(memstat==0)
  end subroutine alloc_i_i_r

  subroutine free_r(a)
    implicit none
    type(rmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(free_r)

    deallocate(a%m,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine free_r

  subroutine alloc_i_c(n,a)
    implicit none
    integer(IK),intent(in)      :: n
    type(cmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(alloc_i_c)

    allocate(a%re(n, n), a%im(n, n), STAT=memstat)
    ASSERT(memstat==0)

    ! FIXME: why?
    a%re =  zero
    a%im =  zero
  end subroutine alloc_i_c

  subroutine alloc_i_i_c(n1,n2,a)
    implicit none
    integer(IK),intent(in)      :: n1,n2
    type(cmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(alloc_i_i_c)

    allocate(a%re(n1, n2), a%im(n1, n2), STAT=memstat)
    ASSERT(memstat==0)

    ! FIXME: why?
    a%re =  zero
    a%im =  zero
  end subroutine alloc_i_i_c

  subroutine free_c(a)
    implicit none
    type(cmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(free_c)

    deallocate(a%re, a%im, STAT=memstat)
    ASSERT(memstat==0)
  end subroutine free_c

  subroutine alloc_i_c_c_plus(n,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
    implicit none
    integer(IK),intent(in)      :: n
    type(cmatrix),intent(inout) ::&
         & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
    optional a2,a3,a4,a5,a6,a7,a8,a9
    !** End of interface ***

    AT_ENTRY(alloc_i_c_c_plus)

    call alloc(n,a0)
    call alloc(n,a1)
    if(present(a2)) call alloc(n,a2)
    if(present(a3)) call alloc(n,a3)
    if(present(a4)) call alloc(n,a4)
    if(present(a5)) call alloc(n,a5)
    if(present(a6)) call alloc(n,a6)
    if(present(a7)) call alloc(n,a7)
    if(present(a8)) call alloc(n,a8)
    if(present(a9)) call alloc(n,a9)
  end subroutine alloc_i_c_c_plus

  subroutine alloc_i_r_r_plus(n,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
    implicit none
    integer(IK),intent(in)      :: n
    type(rmatrix),intent(inout) ::&
         & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
    optional a2,a3,a4,a5,a6,a7,a8,a9
    !** End of interface ***

    AT_ENTRY(alloc_i_r_r_plus)

    call alloc(n,a0)
    call alloc(n,a1)
    if(present(a2)) call alloc(n,a2)
    if(present(a3)) call alloc(n,a3)
    if(present(a4)) call alloc(n,a4)
    if(present(a5)) call alloc(n,a5)
    if(present(a6)) call alloc(n,a6)
    if(present(a7)) call alloc(n,a7)
    if(present(a8)) call alloc(n,a8)
    if(present(a9)) call alloc(n,a9)
  end subroutine alloc_i_r_r_plus

  subroutine free_c_c_plus(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
    implicit none
    type(cmatrix),intent(inout) ::&
         & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
    optional a2,a3,a4,a5,a6,a7,a8,a9
    !** End of interface ***

    AT_ENTRY(free_c_c_plus)

    call pgfree(a0)
    call pgfree(a1)
    if(present(a2)) call pgfree(a2)
    if(present(a3)) call pgfree(a3)
    if(present(a4)) call pgfree(a4)
    if(present(a5)) call pgfree(a5)
    if(present(a6)) call pgfree(a6)
    if(present(a7)) call pgfree(a7)
    if(present(a8)) call pgfree(a8)
    if(present(a9)) call pgfree(a9)
  end subroutine free_c_c_plus

  subroutine free_r_r_plus(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
    implicit none
    type(rmatrix),intent(inout) ::&
         & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
    optional a2,a3,a4,a5,a6,a7,a8,a9
    !** End of interface ***

    AT_ENTRY(free_r_r_plus)

    call pgfree(a0)
    call pgfree(a1)
    if(present(a2)) call pgfree(a2)
    if(present(a3)) call pgfree(a3)
    if(present(a4)) call pgfree(a4)
    if(present(a5)) call pgfree(a5)
    if(present(a6)) call pgfree(a6)
    if(present(a7)) call pgfree(a7)
    if(present(a8)) call pgfree(a8)
    if(present(a9)) call pgfree(a9)
  end subroutine free_r_r_plus

  subroutine alloc_i_d(n,a)
    implicit none
    integer(IK),intent(in)       :: n
    type(rdmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(alloc_i_d)

    allocate(a%d(n),STAT=memstat)
    ASSERT(memstat==0)

    ! FIXME: why?
    a%d = zero
  end subroutine alloc_i_d

  subroutine alloc_i_d_d_plus(n,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
    implicit none
    integer(IK),intent(in)      :: n
    type(rdmatrix),intent(inout) ::&
         & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
    optional a2,a3,a4,a5,a6,a7,a8,a9
    !** End of interface ***

    AT_ENTRY(alloc_i_d_d_plus)

    call alloc(n,a0)
    call alloc(n,a1)
    if(present(a2)) call alloc(n,a2)
    if(present(a3)) call alloc(n,a3)
    if(present(a4)) call alloc(n,a4)
    if(present(a5)) call alloc(n,a5)
    if(present(a6)) call alloc(n,a6)
    if(present(a7)) call alloc(n,a7)
    if(present(a8)) call alloc(n,a8)
    if(present(a9)) call alloc(n,a9)
  end subroutine alloc_i_d_d_plus

  subroutine free_d(a)
    implicit none
    type(rdmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(free_d)

    deallocate(a%d,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine free_d

  subroutine free_d_d_plus(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
    implicit none
    type(rdmatrix),intent(inout) ::&
         & a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
    optional a2,a3,a4,a5,a6,a7,a8,a9
    !** End of interface ***

    AT_ENTRY(free_d_d_plus)

    call pgfree(a0)
    call pgfree(a1)
    if(present(a2)) call pgfree(a2)
    if(present(a3)) call pgfree(a3)
    if(present(a4)) call pgfree(a4)
    if(present(a5)) call pgfree(a5)
    if(present(a6)) call pgfree(a6)
    if(present(a7)) call pgfree(a7)
    if(present(a8)) call pgfree(a8)
    if(present(a9)) call pgfree(a9)
  end subroutine free_d_d_plus

  subroutine alloc_i_h(n,a)
    implicit none
    integer(IK),intent(in)       :: n
    type(chmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(alloc_i_h)

    if(verbose)then
       if(a%n>=0)then
          print *,'n=',a%n
          WARN('alloc_i_h: already alloced?')
       endif
    endif

    a%n    = n
    a%n_ij = n*(n+1)/2

    allocate(a%re(n*(n+1)/2), a%im(n*(n+1)/2), STAT=memstat)
    ASSERT(memstat==0)

    ! FIXME: why?
    a%re = zero
    a%im = zero
  end subroutine alloc_i_h

  subroutine free_h(a)
    implicit none
    type(chmatrix),intent(inout) :: a
    !** End of interface ***

    integer :: memstat

    AT_ENTRY(free_h)

    a%n    = -1
    a%n_ij = -1

    deallocate(a%re, a%im, STAT=memstat)
    ASSERT(memstat==0)
  end subroutine free_h

  pure function size_d(a) result(n)
    !
    ! Size of a diagonal matrix, here length of a diagonal.
    !
    implicit none
    type(rdmatrix), intent(in) :: a
    integer :: n
    ! *** end of interface ***

    n = size(a%d)
  end function size_d

  pure function size_r_i(a, axis) result(n)
    !
    ! Dimension of an rmatrix along the axis.
    !
    implicit none
    type(rmatrix), intent(in) :: a
    integer, intent(in) :: axis
    integer :: n
    ! *** end of interface ***

    !
    ! This  is a  pure function,  we  cannot abort  here, leave  error
    ! checking to intrinsic size:
    !
    n = size(a%m, axis)
  end function size_r_i

  pure function size_c_i(a, axis) result(n)
    !
    ! Dimension of an cmatrix along the axis.
    !
    implicit none
    type(cmatrix), intent(in) :: a
    integer, intent(in) :: axis
    integer :: n
    ! *** end of interface ***

    !
    ! This  is a  pure function,  we  cannot abort  here, leave  error
    ! checking to intrinsic size:
    !
    n = size(a%re, axis)
  end function size_c_i

  !--------------- End of module -------------------------------------
end module matrix_methods

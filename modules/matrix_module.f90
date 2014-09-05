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
module matrix_module
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

  use matrix_types
  use matrix_check
  use matrix_matmult, only: matmult
  use matrix_methods
  use matrix_eigenval
  use matrix_xpack
  use matrix_linsolve, only: linsolve
  implicit none

  save
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  ! from matrix_types:
  public :: rmatrix     ! real matrix
  public :: cmatrix     ! complex matrix
  public :: rdmatrix    ! real diagonal matrix
  public :: chmatrix    ! complex hermitean matrix

  !------------ public functions and subroutines ---------------------


  ! from matrix_types
  public :: init

  ! from matrix_check
  public :: square

  ! from matrix_matmult
  public :: matmult

  ! from matrix_methods:
  public :: matrix                                      ! array -> matrix
  public :: array                                       ! matrix -> array
  public :: assignment(=)                               ! [cr]matrix = array
  public :: operator(*), mult                           ! multiplication
  public :: operator(**)                                ! exponentiation
  public :: operator(+), operator(-)                    ! ...
  public :: tr                                          ! TRANSPOSED, NOT TRACE !!
  public :: sim,bsim                                    ! similarity trafo
  public :: convert                                     ! cm <-> chm conversion
  public :: alloc, pgfree                               ! memory managment
  public :: diag                                        ! extract diagonal
  public :: size                                        ! dimension

  ! from matrix_xpack
  public :: cpack
  public :: cunpack

  ! from matrix_eigenval
  public :: eigs,geigs   ! eigensolver
  public :: svd          ! SVD decomposition

  ! from matrix_linsolve
  public :: linsolve

  !--------------- End of module -------------------------------------
end module matrix_module

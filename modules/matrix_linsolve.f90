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
module matrix_linsolve
  !-------------------------------------------------------------------
  !
  !  Purpose: comfortable call of lapack routine DGESV for solving linear equations
  !          
  !
  !
  !  Module called by: matrix_module
  !
  !
  !  References: ...
  ! 
  !
  !  Author: AN
  !  Date: 07.04.2009
  !
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
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

 interface linsolve
   ! solves A*X = B, A is a n*n matrix,B and X are the same size
   ! B and X can be vectors or matrices
   module procedure linsolve1
   module procedure linsolveN
 end interface
 
  public linsolve

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  function linsolve1(A,B) result (X)
    !
    !  Purpose: solves A*X=B with the routine DGESV; A is a matrix,
    !           b, x are vectors
    !
    !------------ Modules used ------------------- ---------------
    use f77_lapack, only: DGESV
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),    intent(in)    :: A(:,:),B(:)
    real(kind=r8_kind)                   :: X(size(B))
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: n, info, ipiv(size(A,1))
    real(kind=r8_kind)                   :: AC(size(A,1),size(A,2))

    !------------ Executable code ------------------------------------

    ASSERT(size(A,1).eq.size(A,2))
    ASSERT(size(A,1).eq.size(B))

    ! needed variable by DGESV
    n = size(A,1)
    ! copy A and B so they won't be overwriten
    AC = A
    ! B will be overwriten with the result, so x will get the right values this way
    X = B

    call DGESV(n,1,AC,n,ipiv,X,n,info)
    ASSERT(info.eq.0)
  end function linsolve1

  function linsolveN(A,B) result (X)
    !
    !  Purpose: solves A*X=B with the routine DGESV; A,B,X are matrices
    !
    !------------ Modules used ------------------- ---------------
    use f77_lapack, only: DGESV
    implicit none
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),    intent(in)    :: A(:,:),B(:,:)
    real(kind=r8_kind)                   :: X(size(B,1),size(B,2))
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: n, n_b, info, ipiv(size(A,1))
    real(kind=r8_kind)                   :: AC(size(A,1),size(A,2))

    !------------ Executable code ------------------------------------

    ASSERT(size(A,1).eq.size(A,2))
    ASSERT(size(A,1).eq.size(B,1))

    ! needed variables by DGESv
    n = size(A,1)
    n_b = size(B,2)
    ! copy A and B so they won't be overwriten
    AC = A
    ! (B will be overwriten with the result, so x will get the right values this way)
    X = B

    call DGESV(n,n_b,AC,n,ipiv,X,n,info)
    ASSERT(info.eq.0)
  end function linsolveN

  !--------------- End of module -------------------------------------
end module matrix_linsolve

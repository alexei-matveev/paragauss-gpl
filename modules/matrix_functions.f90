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
module matrix_functions
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

# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  implicit none

  save
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: funm!(M,fun)

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------


  !------------ Declaration of constants and variables -----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function funm(M,fun) result(FM)
    !
    ! Computes the user-supplied function fun(M) of a 
    !
    !    SQUARE and SYMMETRIC matrix
    !
    use matrix_eigenval, only: eigs
    use matrix_methods,  only: mult
    implicit none
    real(RK), intent(in) :: M(:,:)
    ! real(RK)           :: fun ! function fun :: real -> real
    ! actual declaration in this interface:
    interface
      function fun(x) result(f)
        use type_module, only: RK => r8_kind
        real(RK), intent(in) :: x
        real(RK)             :: f ! result
      end function fun
    end interface
    real(RK)             :: FM(size(M,1),size(M,2)) ! result
    !** End of interface ***

    real(RK), dimension(size(M,1))           :: md,fd
    real(RK), dimension(size(M,1),size(M,2)) :: U

    integer(IK) :: n,i

    n = size(M,1)
    ASSERT(n==size(M,2))

    ! compute eigenvalues (and eigenvectors):
    call eigs(M, md, U)

    ! apply function to diagonal represenation:
    do i = 1, n
      fd(i) = fun( md(i) )
    enddo

    ! back-transform from diagonal representation:
    !
    !       FM = U * fun(md) * U^T
    !
    ! using the fact that U^T = U^-1
    !
    FM = mult( mult( U, fd ), U, 'nt' )
  end function funm

  !--------------- End of module ----------------------------------
end module matrix_functions

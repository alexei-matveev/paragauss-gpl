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
module matrix_types
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

  use type_module, only:&
       & IK => i4_kind,&
       & RK => r8_kind,&
       & CK => c16_kind ! type specification parameters
  use datatype, only: rmatrix => arrmat2
  implicit none

  save
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

! type, public :: rmatrix
!    real(RK), allocatable :: m(:, :) ! m(n1, n2)
! end type rmatrix
  public :: rmatrix

  type, public :: cmatrix
     real(RK), allocatable :: re(:, :) ! re(n1, n2)
     real(RK), allocatable :: im(:, :) ! im(n1, n2)
  end type cmatrix

  type, public :: rdmatrix ! real diagonal matrix
     real(RK), allocatable :: d(:) ! d(n)
  end type rdmatrix

  type, public :: chmatrix ! complex hermitean matrix
     integer(IK) :: n,n_ij
     real(RK), allocatable :: re(:) ! re(n*(n+1)/2)
     real(RK), allocatable :: im(:) ! im(n*(n+1)/2)
  end type chmatrix

  !------------ Declaration of types ------------------------------
  !
  ! INTERMEDIATE  REPRESENTATION OF MATRICES:
  !
  ! real matrix               <-> real(n1,n2)
  ! complex matrix            <-> real(n1,n2,2)
  ! real diagonal matrix      <-> real(n)
  ! ...................................................
  ! these are not currently supported >>>
  !
  ! real hermitean matrix     <-> real(n*(n+1)/2)
  !                              (always longer (>=) than
  !                               real diagonal matrix!!!)
  ! complex hermitean matrix  <-> real(2, n*(n+1)/2)

  !------------ Declaration of types ------------------------------
  ! added sparse matrices:

  type, public :: SparseRMatrix
     integer(IK)         :: n
     real(RK),pointer    :: c(:)
     integer(IK),pointer :: i(:),j(:)
     ! c(n),i(n),j(n) == n non-zero Cij elements
  end type SparseRMatrix

  type, public :: SparseCMatrix
     integer(IK)         :: n
     complex(CK),pointer :: z(:)
     real(RK),pointer    :: re(:),im(:)
     integer(IK),pointer :: i(:),j(:)
     ! z(n)==(re(n),im(n)),i(n),j(n) == n non-zero Zij elements
  end type SparseCMatrix

  logical, public, protected :: verbose= .false.
  logical, public, protected :: debug  = .true.

  public :: init!(verb, deb)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine init(verb, deb)
    implicit none
    logical,intent(in),optional :: verb,deb
    ! *** end of interface ***

    if(present(verb)) verbose = verb
    if(present(deb))  debug   = deb
  end subroutine init

  !--------------- End of module ----------------------------------
end module matrix_types

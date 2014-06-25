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
module matrix_sparse
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
  use matrix_types, only:&
       & SparseRMatrix,&
       & SparseCMatrix
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  public SparseRMatrix,SparseCMatrix

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  interface alloc ! public
     module procedure alloc_SparseRMatrix
     module procedure alloc_SparseCMatrix
  end interface

  interface dealloc
     module procedure dealloc_SparseRMatrix
     module procedure dealloc_SparseCMatrix
  end interface

  interface free ! public
     module procedure dealloc_SparseRMatrix
     module procedure dealloc_SparseCMatrix
  end interface

  interface spack ! public
     module procedure pack_SparseRMatrix
     module procedure pack_SparseCMatrix
  end interface


  !------------ public functions and subroutines ------------------
  public alloc,&
       & free,&
       & spack

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----
  real(RK),parameter :: AFewDigits = 10000.0_rk

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine pack_SparseRMatrix(a,pa)
    implicit none
    real(RK),intent(in)               :: a(:,:)
    type(SparseRMatrix),intent(inout) :: pa
    ! *** end of interface ***

    integer(IK) :: nn,n,m,i,j
    real(RK)    :: eps
    logical     :: mask(size(a,1),size(a,2))

    eps = AFewDigits * epsilon(AFewDigits)

    n = size(a,1)
    m = size(a,2)

    mask(:,:) = (abs(a(:,:)) > eps)
    ! a(:,:).ne.0.0_rk
       
    nn = count(mask(:,:))
    call alloc(nn,pa)
    if(nn==0)return !<<< for HP

    pa%c = pack(a(:,:),mask(:,:))
    pa%i = pack(spread((/(i,i=1,n)/), 2, m),mask(:,:))
    pa%j = pack(spread((/(j,j=1,m)/), 1, n),mask(:,:))
  end subroutine pack_SparseRMatrix

  subroutine pack_SparseCMatrix(a,pa)
    implicit none
    complex(CK),intent(in)            :: a(:,:)
    type(SparseCMatrix),intent(inout) :: pa
    ! *** end of interface ***

    integer(IK) :: nn,n,m,i,j
    real(RK)    :: eps
    logical     :: mask(size(a,1),size(a,2))

    eps = AFewDigits * epsilon(AFewDigits)

    n = size(a,1)
    m = size(a,2)

    mask(:,:) = (abs(a(:,:)) > eps)
    ! a(:,:).ne.0.0_rk
       
    nn = count(mask(:,:))
    call alloc(nn,pa)
    if(nn==0)return !<<< for HP

    pa%z = pack(a(:,:),mask(:,:))
    pa%re =  REAL(pa%z,KIND=RK)
    pa%im = AIMAG(pa%z)
    pa%i = pack(spread((/(i,i=1,n)/), 2, m),mask(:,:))
    pa%j = pack(spread((/(j,j=1,m)/), 1, n),mask(:,:))
  end subroutine pack_SparseCMatrix

  subroutine alloc_SparseRMatrix(n,pa)
    use error_module
    implicit none
    integer(IK),intent(in)            :: n
    type(SparseRMatrix),intent(inout) :: pa
    ! *** end of interface ***

    integer(IK) :: memstat

    pa%n = n
    allocate(pa%c(n),pa%i(n),pa%j(n),STAT=memstat)
    if(memstat.ne.0) call error("ms/alloc_SparseRMatrix: alloc failed")
  end subroutine alloc_SparseRMatrix

  subroutine dealloc_SparseRMatrix(pa)
    use error_module
    implicit none
    type(SparseRMatrix),intent(inout) :: pa
    ! *** end of interface ***

    integer(IK) :: memstat

    deallocate(pa%c,pa%i,pa%j,STAT=memstat)
    if(memstat.ne.0) call error("ms/dealloc_SparseRMatrix: dealloc failed")
    pa%n = -1
  end subroutine dealloc_SparseRMatrix

  subroutine alloc_SparseCMatrix(n,pa)
    use error_module
    implicit none
    integer(IK),intent(in)            :: n
    type(SparseCMatrix),intent(inout) :: pa
    ! *** end of interface ***

    integer(IK) :: memstat

    pa%n = n
    allocate(pa%z(n),pa%re(n),pa%im(n),pa%i(n),pa%j(n),STAT=memstat)
    if(memstat.ne.0) call error("ms/alloc_SparseCMatrix: alloc failed")
  end subroutine alloc_SparseCMatrix

  subroutine dealloc_SparseCMatrix(pa)
    use error_module
    implicit none
    type(SparseCMatrix),intent(inout) :: pa
    ! *** end of interface ***

    integer(IK) :: memstat

    deallocate(pa%z,pa%re,pa%im,pa%i,pa%j,STAT=memstat)
    if(memstat.ne.0) call error("ms/dealloc_SparseCMatrix: dealloc failed")
    pa%n = -1
  end subroutine dealloc_SparseCMatrix

  !--------------- End of module ----------------------------------
end module matrix_sparse

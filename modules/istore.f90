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
module istore
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

# include "def.h"
  use type_module, only: IK => i4_kind, RK => r8_kind ! type specification parameters
  use readwriteblocked_module, only: readwriteblocked_tapehandle
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  type, private :: file_t
    character(len=16) :: base_name
    !
    ! In spin-orbit case data is usually written in two files
    ! for real and imaginary parts:
    !
    type(readwriteblocked_tapehandle) :: re, im
  end type

  type, public :: istore_t
    private ! iplementation details are subject to change ...
    integer(IK) :: m, n
    real(RK), pointer :: d(:) => NULL() ! column major (m, n)

    !
    ! These hold the state while iterating over the storage:
    !
    integer(IK) :: p ! in [0, n), progress along the n-axis
    type(file_t) :: f ! open while iterating
  end type

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  interface istore_new
    module procedure new
  end interface

  interface istore_apply_and_advance
    module procedure left_apply_and_advance_1
    module procedure left_apply_and_advance_2
    module procedure left_apply_and_advance_3
  end interface

  !------------ public functions and subroutines ------------------

  public :: istore_new
  public :: istore_apply_and_advance

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function new(base_name, m, n, d) result(store)
    !
    ! We dont export intrinsic structure constructor,
    ! see private attribute of type fields.
    !
    implicit none
    character(len=*), intent(in) :: base_name
    integer(IK), intent(in) :: m, n
    real(RK), intent(in), target :: d(:)
    type(istore_t) :: store
    ! *** end of interface ***

    type(readwriteblocked_tapehandle) :: re, im ! in udefined state?

    store = istore_t(m, n, d, 0, file_t(base_name, re, im))
    !print *, "new: store%p =", store%p
  end function new

  subroutine left_apply_and_advance_1(a, b, c)
    !
    ! The case of one rhs b(:)
    !
    implicit none
    type(istore_t), intent(inout) :: a ! (m, n), inout because of position
    real(RK), intent(in) :: b(:) ! (k), k <= n - p
    real(RK), intent(out) :: c(:) ! (m)
    ! *** end of interface ***

    !print *, "left_apply_and_advance_1: store%p =", a%p
    ASSERT(size(b)<=a%n-a%p)
    ASSERT(size(c)==a%m)

    !
    ! 1d array will be viewed as (size(b), 1) matrix:
    !
    call left_apply_and_advance(size(b), 1, a, b, c)
  end subroutine left_apply_and_advance_1

  subroutine left_apply_and_advance_2(a, b, c)
    !
    ! The case of many rhs, b = b(:, 1:nb)
    !
    implicit none
    type(istore_t), intent(inout) :: a ! (m, n), inout because of position
    real(RK), intent(in) :: b(:, :) ! (k, nb), k <= n - p
    real(RK), intent(out) :: c(:, :) ! (m, nb)
    ! *** end of interface ***

    !print *, "left_apply_and_advance_2: store%p =", a%p
    ASSERT(size(b,1)<=a%n-a%p)
    ASSERT(size(c,1)==a%m)
    ASSERT(size(b,2)==size(c,2))

    !
    ! 2d array is a matrix:
    !
    call left_apply_and_advance(size(b,1), size(b,2), a, b, c)
  end subroutine left_apply_and_advance_2

  subroutine left_apply_and_advance_3(a, b, c)
    !
    ! The case of many rhs, with spin-orbit, b = b(1:2, :, 1:nb)
    !
    implicit none
    type(istore_t), intent(inout) :: a ! (m, 2*n), inout because of position
    real(RK), intent(in) :: b(:, :, :) ! (2, k, nb), k <= n - p
    real(RK), intent(out) :: c(:, :) ! (m, nb)
    ! *** end of interface ***

    !print *, "left_apply_and_advance_3: store%p =", a%p
    ASSERT(size(b,1)==2)
    ASSERT(2*size(b,2)<=a%n-a%p)
    ASSERT(size(c,1)==a%m)
    ASSERT(size(b,3)==size(c,2))

    !
    ! 2d array is a matrix:
    !
    call left_apply_and_advance(2*size(b,2), size(b,3), a, b, c)
  end subroutine left_apply_and_advance_3

  subroutine left_apply_and_advance(k, nb, a, b, c)
    !
    !  View arguments in proper shape and dispatch
    !  in-core/on-file verision.
    !
    use options_module, only: options_spin_orbit
    implicit none
    integer(IK), intent(in) :: k ! n or d*(d+1)/2 with irrep dimension d
    integer(IK), intent(in) :: nb ! 1 or more
    type(istore_t), intent(inout) :: a ! (m, n), inout only because of a%p
    real(RK), intent(in) :: b(k, nb)
    real(RK), intent(out) :: c(a%m, nb) ! e.g. (m, 1)
    ! *** end of interface ***

    !print *, "left_apply_and_advance: store%p =", a%p
    if ( associated(a%d) ) then
        ASSERT(size(a%d)==a%m*a%n)
        !print *, "left_apply_and_advance: in core"

        ! just a matrix multiply:
        call mm_in_core(a%m, a%n, k, nb, a%p, a%d, b, c)
    else
        ! need to get the matrix from file:
        if ( options_spin_orbit ) then
            !print *, "left_apply_and_advance: on two files"
            call mm_on_two_files(a%m, a%n, k, nb, a%p, a%f, b, c)
        else
            !print *, "left_apply_and_advance: on file"
            call mm_on_file(a%m, a%n, k, nb, a%p, a%f, b, c)
        endif
    endif

    ! advance the storage pointer (e.g. to next irrep):
    a%p = a%p + k

    ! wrap around:
    a%p = mod(a%p, a%n)
  end subroutine left_apply_and_advance

  subroutine mm_in_core(m, n, k, nb, p, a, b, c)
    !
    ! Basically a reshape and a matrix multiply.
    !
    use matrix_matmult, only: matmult
    implicit none
    integer(IK), intent(in) :: m
    integer(IK), intent(in) :: n
    integer(IK), intent(in) :: k
    integer(IK), intent(in) :: nb
    integer(IK), intent(in) :: p ! current position along n-axis
    real(RK), intent(in) :: a(m, n) ! integrals
    real(RK), intent(in) :: b(k, nb) ! e.g. density matrices
    real(RK), intent(out) :: c(m, nb)
    ! *** end of interface ***

    !print *, "mm_in_core: m, n, k, nb, p =", m, n, k, nb, p

    c = matmult(a(:, p+1:p+k), b)
  end subroutine mm_in_core

  subroutine mm_on_file(m, n, k, nb, p, f, b, c)
    !
    ! Basically a reshape and a matrix multiply,
    ! while fetching the matrix from a file.
    !
    use filename_module, only: tmpfile
    use readwriteblocked_module, only: readwriteblocked_startread, &
        readwriteblocked_read, readwriteblocked_stopread
    implicit none
    integer(IK), intent(in) :: m
    integer(IK), intent(in) :: n
    integer(IK), intent(in) :: k
    integer(IK), intent(in) :: nb ! e.g. 1 or more
    integer(IK), intent(in) :: p ! current position along n-axis
    type(file_t), intent(inout) :: f
    real(RK), intent(in) :: b(k, nb) ! e.g. density matrix
    real(RK), intent(out) :: c(m, nb)
    ! *** end of interface ***

    integer(IK) :: i, j

    ! FIXME: test if reading more than one column works
    real(RK) :: col(m)

    !print *, "mm_on_file: m, n, k, nb, p =", m, n, k, nb, p

    ! open on first iteration:
    if ( p == 0 ) then
        !print *, "mm_on_file: open ", f%base_name
        ! legacy open:
        call readwriteblocked_startread(trim(tmpfile(trim(f%base_name)//'.dat')), f%re)
    endif

    c = 0.0
    do i = 1, k
        ! reading column (p+k) with a legacy read:
        call readwriteblocked_read(col, f%re)

        do j = 1, nb
            c(:, j) = c(:, j) + col(:) * b(i, j)
        enddo
    enddo

    ! close when reaching the end:
    if ( p + k == n ) then
        !print *, "mm_on_file: close ", f%base_name
        ! legacy close:
        call readwriteblocked_stopread(f%re)
    endif
    !
    ! NOTE: pointer p is advanced in the caller code common for both
    !       in-core and on-file versions.
    !
  end subroutine mm_on_file

  subroutine mm_on_two_files(m, n, k, nb, p, f, b, c)
    !
    ! Basically a reshape and a matrix multiply,
    ! while fetching the matrix from a file.
    !
    use filename_module, only: tmpfile
    use readwriteblocked_module, only: readwriteblocked_startread, &
        readwriteblocked_read, readwriteblocked_stopread
    implicit none
    integer(IK), intent(in) :: m
    integer(IK), intent(in) :: n
    integer(IK), intent(in) :: k
    integer(IK), intent(in) :: nb ! e.g. 1 or more
    integer(IK), intent(in) :: p ! current position along n-axis
    type(file_t), intent(inout) :: f
    real(RK), intent(in) :: b(k, nb) ! e.g. density matrix
    real(RK), intent(out) :: c(m, nb)
    ! *** end of interface ***

    integer(IK) :: i, j

    ! FIXME: test if reading more than one column works
    real(RK) :: col(m)

    !print *, "mm_on_two_files: m, n, k, nb, p =", m, n, k, nb, p
    ASSERT(2*(k/2)==k)

    ! open on first iteration:
    if ( p == 0 ) then
        !print *, "mm_on_two_files: open ", f%base_name
        ! legacy open:
        call readwriteblocked_startread(trim(tmpfile(trim(f%base_name)//'_real.dat')), f%re, &
                                        total_length=m*n/2)
        call readwriteblocked_startread(trim(tmpfile(trim(f%base_name)//'_imag.dat')), f%im, &
                                        total_length=m*n/2)
    endif

    c = 0.0
    do i = 1, k, 2
        ! reading column (p+k) with a legacy read:
        call readwriteblocked_read(col, f%re)
        do j = 1, nb
            c(:, j) = c(:, j) + col(:) * b(i, j)
        enddo
    enddo
    do i = 2, k, 2
        ! reading column (p+k) with a legacy read:
        call readwriteblocked_read(col, f%im)
        do j = 1, nb
            c(:, j) = c(:, j) + col(:) * b(i, j)
        enddo
    enddo

    ! close when reaching the end:
    if ( p + k == n ) then
        !print *, "mm_on_two_files: close ", f%base_name
        ! legacy close:
        call readwriteblocked_stopread(f%re)
        call readwriteblocked_stopread(f%im)
    endif
    !
    ! NOTE: pointer p is advanced in the caller code common for both
    !       in-core and on-file versions.
    !
  end subroutine mm_on_two_files

  !--------------- End of module ----------------------------------
end module istore

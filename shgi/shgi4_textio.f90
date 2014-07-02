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
module shgi4_textio
  !-------------------------------------------------------------------
  !
  ! Prepares the input for two-electron SHGI code.
  !
  ! Copyright (c) 2010 Alexei Matveev
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
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: shgi4_write
  public :: shgi4_read

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi4_write(ibasis, iposition, lvalues, exponents, contractions, positions, iou)
    !
    ! Shell structure is represented by indexing into arrays parametrizing
    ! bases and positions of atoms. Specifically, for a shell number "ishl"
    !
    !   ibas = ibasis(ishl)
    !
    ! is an index into arrays (lvalues, exponents, contractions)
    ! to access L-value, exponents and contractions for the shell "ishl"
    !
    !   L = lvalues(ibas)
    !   E = exponents(ibas)%m(:)
    !   C = contractions(ibas)%m(:, :)
    !
    ! Similarly,
    !
    !   ipos = iposition(ishl)
    !
    ! is an index into array of positions:
    !
    !   X = positions(:, ipos)
    !
    use shgi4_driver, only: array1, array2
    implicit none
    integer(IK),  intent(in) :: ibasis(:), iposition(:) ! (nshl)
    integer(IK),  intent(in) :: lvalues(:)              ! (nbas)
    type(array1), intent(in) :: exponents(:)            ! (nbas)%m(nexp)
    type(array2), intent(in) :: contractions(:)         ! (nbas)%m(nexp, ncnt)
    real(RK),     intent(in) :: positions(:, :)         ! (3, natm)
    integer(IK),  intent(in) :: iou
    ! *** end of interface ***

    integer(IK) :: nshl, nbas, nexp, ncnt, natm
    integer(IK) :: i, j

    ! numer of shells:
    nshl = size(ibasis)
    ASSERT(nshl==size(iposition))

    ! dump shell info:
    write(iou, "(I6, ' NSHL')") nshl
    do i = 1, nshl
      write(iou, "(2I6)") ibasis(i), iposition(i)
    enddo

    ! number of different basis sets:
    nbas = size(exponents)
    ASSERT(nbas==size(lvalues))
    ASSERT(nbas==size(contractions))

    ! dump basis info:
    write(iou, "(I6, ' NBAS')") nbas
    do i = 1, nbas
      write(iou, "(I6, ' LVAL')") lvalues(i)

      nexp = size(exponents(i)%m)

      write(iou, "(I6, ' NEXP')") nexp
      write(iou, "(1P, 4E25.16E3)") exponents(i)%m

      ncnt = size(contractions(i)%m, 2)

      write(iou, "(I6, ' NCNT')") ncnt
      do j = 1, ncnt
        write(iou, "(1P, 4E25.16E3)") contractions(i)%m(:, j)
      enddo
    enddo

    ! dump positions:
    natm = size(positions, 2)

    write(iou, "(I6, ' NPOS')") natm
    do i = 1, natm
      write(iou, "(1P, 3E25.16E3)") positions(:, i)
    enddo
  end subroutine shgi4_write

  subroutine shgi4_read(ibasis, iposition, lvalues, exponents, contractions, positions, iou)
    !
    ! Allocates storage and reads the info from IO unit.
    !
    ! Requires allocatable dummy arguments feature (TR15581)
    !
    use shgi4_driver, only: array1, array2
    implicit none
    !
    ! arguments unallocated on entry, rhs comments reflect array shapes on exit:
    !
    integer(IK),  allocatable, intent(out) :: ibasis(:), iposition(:) ! (nshl)
    integer(IK),  allocatable, intent(out) :: lvalues(:)              ! (nbas)
    type(array1), allocatable, intent(out) :: exponents(:)            ! (nbas)%m(nexp)
    type(array2), allocatable, intent(out) :: contractions(:)         ! (nbas)%m(nexp, ncnt)
    real(RK),     allocatable, intent(out) :: positions(:, :)         ! (3, natm)
    integer(IK),  intent(in)  :: iou
    ! *** end of interface ***

    integer(IK) :: nshl, nbas, nexp, ncnt, natm
    integer(IK) :: i, j

    ! read shell info:
    read(iou, *) nshl
    allocate(ibasis(nshl), iposition(nshl))

    do i = 1, nshl
      read(iou, *) ibasis(i), iposition(i)
    enddo

    ! read basis info:
    read(iou, *) nbas
    allocate(lvalues(nbas), exponents(nbas), contractions(nbas))

    do i = 1, nbas
      read(iou, *) lvalues(i)

      read(iou, *) nexp
      allocate(exponents(i)%m(nexp))

      read(iou, *) exponents(i)%m

      read(iou, *) ncnt
      allocate(contractions(i)%m(nexp, ncnt))

      do j = 1, ncnt
        read(iou, *) contractions(i)%m(:, j)
      enddo
    enddo

    ! read positions:
    read(iou, *) natm
    allocate(positions(3, natm))

    do i = 1, natm
      read(iou, *) positions(:, i)
    enddo

!   ! for testing purposes dump the info on stdout:
!   print *, 'shgi4_read: info read from iou=', iou, 'below:'
!   call shgi4_write(ibasis, iposition, lvalues, exponents, contractions, positions, iou=6)
  end subroutine shgi4_read

  !--------------- End of module -------------------------------------
end module shgi4_textio

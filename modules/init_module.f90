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
module init_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: provide an initialization to zero
  !           for a wide variety of data types
  !
  ! The generic subroutine init(xxx) works for
  ! xxx of the following types: 
  !   type(arrmat2)         :: arr(:)
  !   type(arrmat3)         :: arr(:)
  !   type(arrmat2int)      :: arr(:)
  !   real(kind=r8_kind)    :: arr(:)
  !   integer(kind=i4_kind) :: arr(:)
  !   real(kind=r8_kind)    :: arr(:,:)
  !   integer(kind=i4_kind) :: arr(:,:)
  !   real(kind=r8_kind)    :: arr(:,:,:)
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
!================================================================
! End of public interface of module
!================================================================
  !
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

  !------------ Modules used --------------------------------------
  use type_module ! type specification parameters
  use datatype    ! user defined stuff

  implicit none
  private
  save

  !------------ Declaration of public interfaces ------------------
  interface init
     module procedure init_ham
     module procedure init_arrmat2_vec
     module procedure init_arrmat2int_vec
     module procedure init_arr_1dim
     module procedure init_arr_1dim_int
     module procedure init_arr_2dim
     module procedure init_arr_3dim 
     module procedure init_arr_2dim_int
 end interface

!------------ Declaration of private constants and variables ----
 real(kind=r8_kind),    parameter, private :: zero = 0.0_r8_kind
 integer(kind=i4_kind), parameter, private :: zero_int = 0_i4_kind


public :: init,arrmat3,arrmat2,arrmat2int


!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine init_ham(ham)
    !------------ Declaration of formal parameters ---------------
    type(arrmat3),intent(inout) :: ham(:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)      :: i_gamma
    !------------ Executable code --------------------------------
    do i_gamma = 1,ubound(ham,1)
       ham(i_gamma)%m = zero
    enddo
  end subroutine init_ham
  !*************************************************************

  !*************************************************************
  subroutine init_arr_1dim(arr)
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),    intent(out ) :: arr(:)
    !------------ Executable code --------------------------------
    arr = zero
  end subroutine init_arr_1dim
  !*************************************************************

  !*************************************************************
  subroutine init_arr_1dim_int(arr)
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),    intent(out ) :: arr(:)
    !------------ Executable code --------------------------------
    arr = zero_int
  end subroutine init_arr_1dim_int
  !*************************************************************

  !*************************************************************
  subroutine init_arr_2dim(arr)
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(out)   :: arr(:,:)
    !------------ Executable code --------------------------------
    arr = zero
  end subroutine init_arr_2dim
  !*************************************************************

  !*************************************************************
  subroutine init_arr_2dim_int(arr)
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),intent(out)   :: arr(:,:)
    !------------ Executable code --------------------------------
    arr = zero_int
  end subroutine init_arr_2dim_int
  !*************************************************************

  !*************************************************************
  subroutine init_arr_3dim(arr)
    !------------ Declaration of formal parameters ---------------
    real(kind=r8_kind),intent(out)   :: arr(:,:,:)
    !------------ Executable code --------------------------------
    arr = zero
  end subroutine init_arr_3dim
  !*************************************************************

  !*************************************************************
  subroutine init_arrmat2_vec(arr)
    !------------ Declaration of formal parameters ---------------
    type(arrmat2),intent(inout)       :: arr(:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)             :: i
    !------------ Executable code --------------------------------
    do i=lbound(arr,1),ubound(arr,1)
       arr(i)%m = zero
    enddo
  end subroutine init_arrmat2_vec
  !*************************************************************

  !*************************************************************
  subroutine init_arrmat2int_vec(arr)
    !------------ Declaration of formal parameters ---------------
    type(arrmat2int),intent(inout)    :: arr(:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)             :: i
    !------------ Executable code --------------------------------
    do i=lbound(arr,1),ubound(arr,1)
       arr(i)%m = zero_int
    enddo
  end subroutine init_arrmat2int_vec
  !*************************************************************

end module init_module
     

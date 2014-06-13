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
module direct_scf
  !----------------------------------------------------------------------------+
  !                                                                            !
  !  Purpose: ...                                                              !
  !                                                                            !
  !                                                                            !
  !  Module called by: ...                                                     !
  !                                                                            !
  !                                                                            !
  !  References: ...                                                           !
  !                                                                            !
  !                                                                            !
  !  Author: ...                                                               !
  !  Date: ...                                                                 !
  !                                                                            !
  !                                                                            !
  !----------------------------------------------------------------------------+
  ! Modifications                                                              !
  !----------------------------------------------------------------------------+
  !                                                                            !
  ! Modification (Please copy before editing)                                  !
  ! Author: ...                                                                !
  ! Date:   ...                                                                !
  ! Description: ...                                                           !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
# include "def.h"
  use type_module,         only: IK=>i4_kind, RK=>r8_kind                      !
  !                                                                            !
  use datatype, only: arrmat2, arrmat3
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  !============================================================================+
  ! DEFINITION OF SAVED QUANTITIES                                             !
  !============================================================================+
  !                                                                            !
  type(arrmat2), allocatable, save :: Jlast(:)                                 !
  type(arrmat3), allocatable, save :: Klast(:)                                 !
  type(arrmat3), allocatable, save :: Plast(:)                                 !
  !                                                                            !
  interface allocate_as                                                        !
    module procedure allocate_as2                                              !
    module procedure allocate_as3                                              !
    module procedure allocate_as32                                             !
  end interface                                                                !
  !                                                                            !
  public :: allocate_as                                                        !
  public :: get_diff_P                                                         !
  public :: get_full_K                                                         !
  public :: get_full_J                                                         !
  !                                                                            !
  contains                                                                     !
  !                                                                            !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine get_diff_P( P_act, Pdiff )                                        !
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    ! Store density matrix P_act and calculate differential Pdiff              !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    type(arrmat3), intent(in)               :: P_act(:)                        !
    !                                                                          !
    type(arrmat3), allocatable, intent(out) :: Pdiff(:)                        !
    !                                                                          !
    integer(IK)                             :: irr                             !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    if (.not. allocated( Plast )) then                                         !
      call allocate_as(P_act, Plast )                                          !
    end if                                                                     !
    if (.not. allocated( Pdiff )) then                                         !
      call allocate_as(P_act, Pdiff)                                           !
    end if                                                                     !
    !                                                                          !
    do irr = 1, size( Plast )                                                  !
      Pdiff(irr)%m = P_act(irr)%m - Plast(irr)%m                               !
      Plast(irr)%m = P_act(irr)%m                                              !
    enddo                                                                      !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
  end subroutine get_diff_P                                                    !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine get_full_K( K_mat )                                               !
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    ! Increment old exchange matrix by differential, store and deliver result  !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    type(arrmat3), intent(inout) :: K_mat(:)                                   !
    !                                                                          !
    !                                                                          !
    integer(IK)                  :: irr                                        !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    if (.not. allocated( Klast )) then                                         !
      call allocate_as(K_mat, Klast)                                           !
    end if                                                                     !
    !                                                                          !
    do irr = 1, size( Klast)                                                   !
      Klast(irr)%m = Klast(irr)%m + K_mat(irr)%m                               !
      K_mat(irr)%m = Klast(irr)%m                                              !
    enddo                                                                      !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
  end subroutine get_full_K                                                    !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine get_full_J( J_mat )                                               !
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    ! Increment old coulomb matrix by differential, store and deliver result   !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    type(arrmat2), intent(inout) :: J_mat(:)                                   !
    !                                                                          !
    !                                                                          !
    integer(IK)                  :: irr                                        !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    if (.not. allocated( Jlast )) then                                         !
      call allocate_as(J_mat, Jlast)                                           !
    end if                                                                     !
    !                                                                          !
    do irr = 1, size( Jlast )                                                  !
      Jlast(irr)%m = Jlast(irr)%m + J_mat(irr)%m                               !
      J_mat(irr)%m = Jlast(irr)%m                                              !
    enddo                                                                      !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
  end subroutine get_full_J                                                    !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
#if 0
  subroutine store_ham_exx( mat )                                              !
    !--------------------------------------------------------------------------+
    ! save ham_exx for later use                                               !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(arrmat3), intent(in) :: mat(:)                                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)            :: irr                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    if (.not. allocated(ham_exx_saved)) then                                   !
      call allocate_as( mat, ham_exx_saved )                                   !
    end if                                                                     !
    !                                                                          !
    do irr = 1, size( ham_exx_saved )                                          !
      !                                                                        !
      ham_exx_saved(irr)%m = mat(irr)%m                                        !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine store_ham_exx                                                 !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine update_densmat_saved_to_densmat_diff(densmat, c_EXX)              !
    !--------------------------------------------------------------------------+
    ! save difference for later use                                            !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),      intent(in) :: c_EXX                                         !
    !                                                                          !
    type(arrmat3), intent(in) :: densmat(:)                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)            :: irr                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    if (.not. allocated(densmat_saved)) then                                   !
      call allocate_as(densmat, densmat_saved)                                 !
    end if                                                                     !
    !                                                                          !
    do irr = 1, size(densmat_saved)                                            !
      !                                                                        !
      densmat_saved(irr)%m = (densmat(irr)%m - densmat_saved(irr)%m) * c_EXX   !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine update_densmat_saved_to_densmat_diff                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine update_densmat_saved(densmat)                                     !
    !--------------------------------------------------------------------------+
    ! save difference for later use                                            !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(arrmat3), intent(in) :: densmat(:)                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)            :: irr                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    if (.not. allocated(densmat_saved)) then                                   !
      call allocate_as(densmat, densmat_saved)                                 !
    end if                                                                     !
    !                                                                          !
    do irr = 1, size(densmat_saved)                                            !
      !                                                                        !
      densmat_saved(irr)%m = densmat(irr)%m                                    !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine update_densmat_saved                                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine ham_coul_diff_to_total(ham_coul)                                  !
    !--------------------------------------------------------------------------+
    ! Expects differential ham_coul, saves and gives back total ham_coul       !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(arrmat2), intent(inout) :: ham_coul(:)                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)            :: irr                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    if (.not. allocated(ham_coul_saved)) then                                  !
      call allocate_as(ham_coul, ham_coul_saved)                               !
    endif                                                                      !
    !                                                                          !
    do irr = 1, size(ham_coul)                                                 !
      !                                                                        !
      ham_coul_saved(irr)%m = ham_coul_saved(irr)%m + ham_coul(irr)%m          !
      ham_coul(irr)%m       = ham_coul_saved(irr)%m                            !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine ham_coul_diff_to_total                                        !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine ham_exx_diff_to_total(ham_exx)                                    !
    !--------------------------------------------------------------------------+
    ! Expects differential ham_exx, saves and gives back total ham_exx         !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(arrmat3), intent(inout) :: ham_exx(:)                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)            :: irr                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    if (.not. allocated(ham_exx_saved)) then                                   !
      call allocate_as(ham_exx, ham_exx_saved)                                 !
    endif                                                                      !
    !                                                                          !
    do irr = 1, size(ham_exx)                                                  !
      !                                                                        !
      ham_exx_saved(irr)%m = ham_exx_saved(irr)%m + ham_exx(irr)%m             !
      ham_exx(irr)%m       = ham_exx_saved(irr)%m                              !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine ham_exx_diff_to_total                                         !
#endif
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine allocate_as2(a, b)                                                !
    !--------------------------------------------------------------------------+
    ! Allocate b to look like a                                                !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(arrmat2), intent(in)  :: a(:)                                         !
    type(arrmat2), allocatable :: b(:)                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK) :: size1, size2, i                                             !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ASSERT(.not.allocated(b))
    allocate(b(size(a)))                                                       !
    !                                                                          !
    do i = 1, size(a)                                                          !
      size1 = size(a(i)%m, 1)                                                  !
      size2 = size(a(i)%m, 2)                                                  !
      !                                                                        !
      allocate(b(i)%m(size1, size2))                                           !
      !                                                                        !
      ! FIXME: initialize as well, see *diff_to_total                          !
      b(i)%m = 0.0                                                             !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine allocate_as2                                                  !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine allocate_as3(a, b)                                                !
    !--------------------------------------------------------------------------+
    ! Allocate b to look like a                                                !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(arrmat3), intent(in)  :: a(:)                                         !
    type(arrmat3), allocatable :: b(:)                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK) :: size1, size2, size3, i                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ASSERT(.not.allocated(b))
    allocate(b(size(a)))                                                       !
    !                                                                          !
    do i = 1, size(a)                                                          !
      size1 = size(a(i)%m, 1)                                                  !
      size2 = size(a(i)%m, 2)                                                  !
      size3 = size(a(i)%m, 3)                                                  !
      !                                                                        !
      allocate(b(i)%m(size1, size2, size3))                                    !
      !                                                                        !
      ! FIXME: initialize as well, see *diff_to_total                          !
      b(i)%m = 0.0                                                             !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine allocate_as3                                                  !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine allocate_as32(a, b)                                               !
    ! Allocate b to look like a                                                !
    implicit none                                                              !
    !                                                                          !
    type(arrmat3), intent(in)  :: a(:)                                         !
    type(arrmat2), allocatable :: b(:)                                         !
    !                                                                          !
    integer(IK)                :: i                                            !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    allocate(b(size(a)))                                                       !
    !                                                                          !
    do i = 1, size(a)                                                          !
      !                                                                        !
      allocate( b(i)%m(size( a(i)%m, 1 ), size(a(i)%m, 2) ) )                  !
      b(i)%m = 0.0                                                             !
      !                                                                        !
    enddo                                                                      !
    !                                                                          !
  end subroutine allocate_as32                                                 !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module direct_scf

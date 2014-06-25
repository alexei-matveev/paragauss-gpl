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
module cpks_common
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
  use type_module, only: &
       IK => i4_kind, &
       RK => r8_kind ! type specification parameters
  use datatype, only: arrmat2, arrmat3
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  ! (yet another) density matrix, and energy-weighted density matrix
  ! that are used in gradient integral part
  type(arrmat3), allocatable, public :: cpks_p0(:) ! (n_irr)%m(dim_irr,dim_irr,n_spn)
  type(arrmat2), allocatable, public :: cpks_w0(:) ! (n_irr)%m(dim_irr,dim_irr)

  ! Gradients of the density matrix, and energy-weighted density matrix
  ! that are used in  the second derivative integral part
  type(arrmat3), allocatable, public :: cpks_p1(:,:) ! (n_irr,n_gra)%m(dim_irr,dim_irr,n_spn)
  type(arrmat2), allocatable, public :: cpks_w1(:,:) ! (n_irr,n_gra)%m(dim_irr,dim_irr)


  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: cpks_alloc_p0w0
  public :: cpks_free_p0w0

  public :: cpks_alloc_p1w1
  public :: cpks_free_p1w1

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine cpks_alloc_p0w0(dims,n_spn)
    !  Purpose: allocates storage for P1/W1
    implicit none
    integer(IK), intent(in) :: dims(:) ! (n_irreps)
    integer(IK), intent(in) :: n_spn   ! 1 or 2 (for MDA only)
    ! *** end of interface ***

    integer(IK) :: n_irr,irr,dim
    integer(IK) :: memstat

    n_irr = size(dims)

    ASSERT(.not.allocated(cpks_p0))
    ASSERT(.not.allocated(cpks_w0))

    allocate( cpks_p0(n_irr) &
            , cpks_w0(n_irr) &
            , STAT=memstat   &
            )
    ASSERT(memstat==0)

    do irr=1,n_irr
       dim = dims(irr)
       allocate( cpks_p0(irr)%m(dim,dim,n_spn) &
               , cpks_w0(irr)%m(dim,dim)       &
               , STAT=memstat                  &
               )
       ASSERT(memstat==0)
    enddo
  end subroutine cpks_alloc_p0w0

  subroutine cpks_free_p0w0()
    !  Purpose: deallocates storage for P1/W1
    implicit none
    ! *** end of interface ***

    integer(IK) :: n_irr,irr
    integer(IK) :: memstat

    ASSERT(allocated(cpks_p0))
    ASSERT(allocated(cpks_w0))

    n_irr = size(cpks_p0)

    do irr=1,n_irr
       deallocate( cpks_p0(irr)%m &
                 , cpks_w0(irr)%m &
                 , STAT=memstat   &
                 )
       ASSERT(memstat==0)
    enddo

    deallocate( cpks_p0      &
              , cpks_w0      &
              , STAT=memstat &
              )
    ASSERT(memstat==0)
  end subroutine cpks_free_p0w0

  subroutine cpks_alloc_p1w1(dims,n_spn,n_gra)
    !  Purpose: allocates storage for P1/W1
    implicit none
    integer(IK), intent(in) :: dims(:) ! (n_irreps)
    integer(IK), intent(in) :: n_spn   ! 1 or 2 (for MDA only)
    integer(IK), intent(in) :: n_gra
    ! *** end of interface ***

    integer(IK) :: n_irr,irr,dim,igr
    integer(IK) :: memstat

    print *,'cpks_alloc_p1w1:  dims=',dims
    print *,'cpks_alloc_p1w1: n_spn=',n_spn
    print *,'cpks_alloc_p1w1: n_gra=',n_gra

    n_irr = size(dims)

    ASSERT(.not.allocated(cpks_p1))
    ASSERT(.not.allocated(cpks_w1))

    allocate( cpks_p1(n_irr,n_gra) &
            , cpks_w1(n_irr,n_gra) &
            , STAT=memstat         &
            )
    ASSERT(memstat==0)

    do igr=1,n_gra
    do irr=1,n_irr
       dim = dims(irr)
       allocate( cpks_p1(irr,igr)%m(dim,dim,n_spn) &
               , cpks_w1(irr,igr)%m(dim,dim)       &
               , STAT=memstat                      &
               )
       ASSERT(memstat==0)
    enddo
    enddo
  end subroutine cpks_alloc_p1w1

  subroutine cpks_free_p1w1()
    !  Purpose: deallocates storage for P1/W1
    implicit none
    ! *** end of interface ***

    integer(IK) :: n_irr,irr,n_gra,igr
    integer(IK) :: memstat

    ASSERT(allocated(cpks_p1))
    ASSERT(allocated(cpks_w1))

    n_irr = size(cpks_p1,1)
    n_gra = size(cpks_p1,2)

    do igr=1,n_gra
    do irr=1,n_irr
       deallocate( cpks_p1(irr,igr)%m &
                 , cpks_w1(irr,igr)%m &
                 , STAT=memstat       &
                 )
       ASSERT(memstat==0)
    enddo
    enddo

    deallocate( cpks_p1      &
              , cpks_w1      &
              , STAT=memstat &
              )
    ASSERT(memstat==0)
  end subroutine cpks_free_p1w1

  !--------------- End of module ----------------------------------
end module cpks_common

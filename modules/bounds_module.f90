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
module bounds_module
  !
  !  Purpose: provide the bounds when sending the splitted
  !           array of fit coefficients to the slaves
  !           (and also the master if comm_master_work is TRUE)
  !
  !  Contents:
  !     - Public types: fit_bounds
  !
  !     - Routines:
  !       (i) PUBLIC
  !                   bounds_calc   -> calculate the bounds to split up the
  !                                    arrays of fit_coeffs when they are to
  !                                    be sent to the slave processors.
  !                   get_bounds    -> access to BOUNDS_CH or BOUNDS_XC
  !
  !                   bounds_send   -> send the bounds information to the
  !                                    slaves
  !
  !  Module called by: build_hamiltonian
  !                    chargefit
  !                    split_tapes in 'convert'
  !
  !
  !  Author: Folke Noertemann
  !  Date: 10/95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   8/97
  ! Description: bounds_free introduced to release the private
  !              pointer bounds_{ch,xc}%item_arr(:)
  !
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: treatment of exchange integral bounds adapted
  !              to the use of the extended model density approach
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use fit_coeff_module, only: fit
  implicit none
  private
  save
  !== Interrupt of public interface of module =====================

  !------------ Declaration of types ------------------------------

  integer(i4_kind), allocatable, public, protected :: bounds_ch(:) ! (NPROCS)
  integer(i4_kind), allocatable, public, protected :: bounds_xc(:) ! (NPROCS)

  !------------ public functions and subroutines ------------------
  public :: bounds_calc!(n_fit_struct) allocates or re-allocates public arrays
  public :: bounds_free!()

!================================================================
! End of public interface of module
!================================================================

  !--- public types of formal parameters of subroutines ---------
  public :: fit

  !------------ Declaration of private constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function work_shares(n, np) result(work)
    !
    ! Legacy algorithm for sharing N jobs among NP workers
    !
    implicit none
    integer(i4_kind), intent(in) :: n, np
    integer(i4_kind)             :: work(np) ! sum(work) == n
    ! *** end of interface ***

    integer(i4_kind) :: i, workers
    real(r8_kind) :: rem

    workers = np
    rem = real(n, kind=r8_kind)

    do i = 1, np
       ! assign a fraction:
       work(i) = ceiling(rem / workers)

       ! the remainder to be distributed among the others:
       workers = workers - 1
       rem = rem - real(work(i), kind=r8_kind)
    enddo
    ASSERT(sum(work)==n)
  end function work_shares

  subroutine bounds_calc(n_fit)
    !  Purpose: Calculate the bounds of a vector that
    !           is to be splitted up and sent to the
    !           slave processors.
    !  Input parameter (not modified on output):
    !  nfit             fit basis inforamtion -> fit_coeff_module
    !  Subroutine called by: build_hamiltonian,
    !                        split_tapes ( see : convert.f90)
    !                        chargefit
    implicit none
    !------------ Declaration of formal parameters ---------------
    type(fit), intent(in) :: n_fit
    !** End of interface *****************************************

    ASSERT(n_fit%n_ch>=0)
    call bounds_calc1(n_fit%n_ch, bounds_ch)

    ASSERT(n_fit%n_xc>=0)
    call bounds_calc1(n_fit%n_xc, bounds_xc)
  end subroutine bounds_calc

  subroutine bounds_calc1(n, bounds)
    !  Purpose: Calculate the bounds of a vector that
    !           is to be splitted up and sent to the
    !           slave processors.
    !  Input parameter (not modified on output):
    !  n        fit basis dimension
    !
    use comm, only: comm_size
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: n
    integer(i4_kind), allocatable, intent(inout) :: bounds(:)
    !** End of interface *****************************************

    integer (i4_kind) :: workers, alloc_stat

    if (allocated (bounds)) then
      WARN ("already done, repeating")
      deallocate (bounds, STAT=alloc_stat)
      ASSERT(alloc_stat.eq.0)
    endif

    workers = comm_size()

    allocate (bounds(workers), STAT=alloc_stat)
    ASSERT(alloc_stat.eq.0)

    bounds(:) = work_shares (n, workers)
  end subroutine bounds_calc1

  subroutine bounds_free()
    !  Purpose: Release the private bounds pointer
    !  Input parameter (not modified on output):
    !  nfit             fit basis inforamtion -> fit_coeff_module
    !  Subroutine called by: build_hamiltonian,
    !                        split_tapes ( see : convert.f90)
    !                        chargefit
    implicit none
    !** End of interface *****************************************

    integer :: alloc_stat

    if(allocated(bounds_ch)) then
       deallocate(bounds_ch, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if

    if(allocated(bounds_xc)) then
       deallocate(bounds_xc, STAT=alloc_stat)
       ASSERT(alloc_stat.eq.0)
    end if
  end subroutine bounds_free

!--------------- End of module ----------------------------------
end module bounds_module

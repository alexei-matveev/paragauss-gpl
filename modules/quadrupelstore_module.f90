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
module  quadrupelstore_module
!---------------------------------------------------------------
!
!  Purpose: A two-way associative store type to hold qudruples and
!           methods for it are also defined. Elements in the store
!           are sorted by ascending priorities. There are methods
!           to construct and destruct a store, to add new
!           quadrupels and to return and remove the first
!           and last element of store. See below.
!
!  Module called by: int_distribute_module
!
!  References: Publisher Document: Concepts of Integral Part
!
!  Author: TB
!  Date: 5/96
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
  use type_module ! type specification parameters
  use quadrupel_module

  implicit none

  private         ! by default, all names are private
  save

!== Interrupt end of public interface of module =================


  !------------ Declaration of types ----------------------------
  type, private :: cell_type
     ! one cell of store
     private
     !== Interrupt of public interface of module ================
     real(kind=r8_kind)       :: priority
     type(quadrupel_type)     :: q
     type(cell_type), pointer :: pred => NULL(), next => NULL()
     !== Interrupt end of public interface of module ============
  end type cell_type

  type, public :: quadrupelstore_type
     private
     !== Interrupt of public interface of module ================
     type(cell_type), pointer                :: first => NULL(), last => NULL(), &
                                                actual => NULL(), putactual => NULL()
     integer(kind=i4_kind)                   :: n_elements
     logical                                 :: is_fitfunction
     logical                                 :: read_done
     !== Interrupt end of public interface of module ============
  end type quadrupelstore_type

  public quadrupel_type
  ! export type from quadrupel_module necessary for formal
  ! parameter of subroutines below

  !------------ Interface statements ----------------------------

  !------------ public functions and subroutines ----------------
  public quadrupelstore_construct, &
         quadrupelstore_give_num

!================================================================
! End of public interface of module
!================================================================

  ! for knowing which one number belongs to the current job
  integer(kind=i4_kind) :: current_num

  !--------------------------------------------------------------
  !------------ Subroutines -------------------------------------
contains



  !**************************************************************
  subroutine quadrupelstore_construct(sa, lmin, lmax, quadrupels_triangular,  is_fitfunction )
    !  Purpose: allocates and initialices store s
    use unique_atom_module, only: N_unique_atoms
    implicit none
    !------------ Declaration of formal parameters -------------

    type(quadrupel_type), allocatable   :: sa(:)
    logical, intent(in) :: is_fitfunction ! default .false.
    logical, intent(in) :: quadrupels_triangular
    integer(kind=i4_kind), intent(in) :: lmin, lmax(:)
    !** End of interface ****************************************
    type(quadrupelstore_type)    :: s
    integer(kind=i4_kind) :: n_quads, ierr
    integer(kind=i4_kind) :: i_ua1, i_l1,  i_ua2, i_l2
    integer(kind=i4_kind) :: uaL1, uaL2

    ! First build object of old quadrupelstore_type, which can be sorted easily
    s%n_elements = 0
    s%is_fitfunction = is_fitfunction
    s%read_done = .true.

    ! fill quadrupel stack directly here, not in int_distribute_module:
    ! reset quadrupel count:
    n_quads = 0

    ! loop over shells 1:
    uaL1 = 0
    do i_ua1 = 1, N_unique_atoms
       do i_l1 = lmin, lmax(i_ua1)
          uaL1 = uaL1 + 1

          ! loop over shells 2:
          uaL2 = 0
          do i_ua2 = 1, N_unique_atoms
             do i_l2 = lmin, lmax(i_ua2)
                uaL2 = uaL2 + 1

                if ( quadrupels_triangular .and. uaL2  > uaL1 ) cycle

                ! increment quadrupel count:
                n_quads = n_quads + 1

                ! store quadrupel description:
                call quadrupelstore_put( s, &
                     i_ua1, i_l1, i_ua2, i_l2 )

             enddo
          enddo ! loop over second shell

       enddo
    enddo ! loop over first shell

    ! set global job counter on first element
    current_num = 1

    ! move quadrupelstore_type content in array:
    allocate(sa(n_quads), stat = ierr)
    ASSERT(ierr==0)

    do i_ua1 = 1, n_quads - 1
       sa(i_ua1) = s%last%q
       s%last => s%last%pred
    enddo
    ! last element has no valid %pred
    sa(n_quads) = s%last%q
    ! quadruplelstore object is not needed any more
    call quadrupelstore_destruct(s)
  end subroutine quadrupelstore_construct
  !**************************************************************


  !**************************************************************
  subroutine quadrupelstore_put(s,ua1,l1,ua2,l2)
    !  Purpose: adds quadrupel q to store
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupelstore_type)          :: s
    integer(kind=i4_kind), intent(in)  :: ua1,l1,ua2,l2
    !** End of interface ****************************************

    call put_quadrupel( s, quadrupel_type(ua1,l1,ua2,l2) )
  end subroutine quadrupelstore_put
  !**************************************************************


  !**************************************************************
  subroutine put_quadrupel(s,q)
    !  Purpose: adds quadrupel q to store
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupelstore_type)        :: s
    type(quadrupel_type), intent(in) :: q
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    type(cell_type),  pointer :: c
    integer(kind=i4_kind)     :: status
    !------------ Executable code -------------------------------
    allocate( c, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "put_quadrupel: allocate failed")
    c%priority = quadrupel_priority(q,s%is_fitfunction)
    c%q = q
    call put_cell(s,c)
  end subroutine put_quadrupel
  !**************************************************************


  !**************************************************************
  subroutine put_cell(s,c)
    !  Purpose: puts quadrupel q in actual position of store
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupelstore_type) :: s
    type(cell_type), pointer  :: c
    !** End of interface ****************************************
    !------------ Executable code -------------------------------
    if ( s%n_elements .eq. 0 ) then
       s%first => c
       s%last => c
       s%actual => c
    else
       if ( s%putactual%priority .lt. c%priority ) then
          do while ( .not. associated(s%putactual,s%last) )
             if (  s%putactual%next%priority .ge. c%priority ) exit
             s%putactual => s%putactual%next
          enddo
          ! putactual should point to cell after which to add
          if ( .not. associated(s%putactual,s%last) ) then
             c%next => s%putactual%next
             c%next%pred => c
          else
             s%last => c
          endif
          c%pred => s%putactual
          s%putactual%next => c
       else
          do while ( .not. associated(s%putactual,s%first) )
             if (  s%putactual%pred%priority .lt. c%priority ) exit
             s%putactual => s%putactual%pred
          enddo
          ! putactual should point to cell before which to add
          if ( .not. associated(s%putactual,s%first) ) then
             c%pred => s%putactual%pred
             c%pred%next => c
          else
             s%first => c
          endif
          c%next => s%putactual
          s%putactual%pred => c
       endif
    endif
    s%putactual => c
    s%n_elements = s%n_elements + 1
  end subroutine put_cell
  !**************************************************************


  !**************************************************************
  function quadrupelstore_give_num(s, num) result(q)
    !  Purpose: gives quadrupel of number num back
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type), intent(in) :: s(:)
    integer(kind=i4_kind), intent(in) :: num
    type(quadrupel_type) :: q
    !** End of interface ****************************************
    !------------ Executable code -------------------------------
    q = s(num)
  end function quadrupelstore_give_num
  !**************************************************************



  !**************************************************************
  subroutine quadrupelstore_destruct(s)
    !  Purpose: deallocates store s and all quadruples contained
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupelstore_type) :: s
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    type(cell_type), pointer :: cell, next
    !------------ Executable code -------------------------------

     cell => s%first
     do while ( associated(cell) )
        ! save pointer to next before deallocation:
        next => cell%next

        deallocate(cell)

        ! advance:
        cell => next
     enddo
  end subroutine quadrupelstore_destruct
  !**************************************************************

end module quadrupelstore_module

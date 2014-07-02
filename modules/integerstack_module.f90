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
module  integerstack_module
!---------------------------------------------------------------
!
!  Purpose: stack to hold integer numbers i
!           add new number i to stack with push(i)
!           and get latest number i and delete it from
!           stack with pop(i).
!
!
!  Module called by: int_distribute_module
!
! 
!  Author: TB
!  Date: 5/96
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

  use type_module ! type specification parameters

  implicit none

  private         ! by default, all names are private
  save

!== Interrupt end of public interface of module =================


  !------------ Declaration of types ----------------------------
  !== Interrupt of public interface of module ================
  type, private :: integerstack_cell_type
     type(integerstack_cell_type), pointer :: pred
     integer(kind=i4_kind)                 :: content
  end type integerstack_cell_type
  !== Interrupt end of public interface of module ============

  type, public :: integerstack_type
     private
     !== Interrupt of public interface of module ================
     type(integerstack_cell_type), pointer :: first,last
     integer(kind=i4_kind)                   :: n_elements
     !== Interrupt end of public interface of module ============
  end type integerstack_type
        

  !------------ public functions and subroutines ----------------
  public integerstack_construct, &
         integerstack_destruct, &
         integerstack_push, &
         integerstack_pop, &
         integerstack_is_empty


  !===================================================================
  ! End of public interface of module
  !===================================================================



  !--------------------------------------------------------------
  !------------ Subroutines -------------------------------------
contains



  !**************************************************************
  subroutine integerstack_construct(s)
    !  Purpose: allocates and initialices stack s
    implicit none
    !------------ Declaration of formal parameters -------------
    type(integerstack_type), pointer :: s
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind)                :: status
    !------------ Executable code -------------------------------
    allocate( s, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integerstack_construct: allocate failed")
    s%n_elements = 0
  end subroutine integerstack_construct
  !**************************************************************


  !**************************************************************
  subroutine integerstack_destruct(s)
    !  Purpose: deallocates stack s
    implicit none
    !------------ Declaration of formal parameters -------------
    type(integerstack_type), pointer :: s
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    type(integerstack_cell_type), pointer :: c
    integer(kind=i4_kind)                 :: status
    !------------ Executable code -------------------------------
    if ( s%n_elements .gt. 0 ) then
       do while ( .not. associated(s%first,s%last) )
          c => s%last
          s%last => s%last%pred
          deallocate( c, stat=status )
          if ( status .ne. 0 ) call error_handler( &
               "integerstack_dectruct: deallocate 1 of cell failed")
       enddo
       c => s%first
       deallocate( c, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "integerstack_dectruct: deallocate 2 of cell failed")
    endif
    deallocate( s, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integerstack_destruct: deallocate failed")
  end subroutine integerstack_destruct
  !**************************************************************


  !**************************************************************
  subroutine integerstack_push(s,i)
    !  Purpose: puts i on top of stack
    implicit none
    !------------ Declaration of formal parameters --------------
    type(integerstack_type), pointer  :: s
    integer(kind=i4_kind), intent(in) :: i
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    type(integerstack_cell_type), pointer :: c
    integer(kind=i4_kind)                 :: status
    !------------ Executable code -------------------------------
    allocate( c, stat=status )
    if ( status .ne. 0 ) call error_handler( &
         "integerstack_push: allocate of cell failed")
    c%content = i
    if ( s%n_elements .eq. 0 ) then
       s%first => c
    else
       c%pred => s%last
    endif
    s%last => c
    s%n_elements = s%n_elements + 1
  end subroutine integerstack_push
  !**************************************************************


  !**************************************************************
  logical function integerstack_pop(s,i)
    !  Purpose: gives i from top of stack, removing it from stack 
    !           and returns if there was any element in stack.
    implicit none
    !------------ Declaration of formal parameters --------------
    type(integerstack_type), pointer   :: s
    integer(kind=i4_kind), intent(out) :: i
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    type(integerstack_cell_type), pointer :: c
    integer(kind=i4_kind)                 :: status
    !------------ Executable code -------------------------------
    integerstack_pop = s%n_elements .gt. 0
    if ( integerstack_pop ) then
       i = s%last%content
       c => s%last
       if ( .not. associated(s%last,s%first) ) s%last => s%last%pred
       deallocate( c, stat=status )
       if ( status .ne. 0 ) call error_handler( &
            "integerstack_pop: deallocate of cell failed")
       s%n_elements = s%n_elements - 1
    endif
  end function integerstack_pop
  !**************************************************************


  !**************************************************************
  logical function integerstack_is_empty(s)
    !  Purpose: returns if there are no element in stack.
    implicit none
    !------------ Declaration of formal parameters --------------
    type(integerstack_type), pointer   :: s
    !** End of interface ****************************************
    !------------ Executable code -------------------------------
    integerstack_is_empty = s%n_elements .eq. 0
  end function integerstack_is_empty
  !**************************************************************


end module integerstack_module

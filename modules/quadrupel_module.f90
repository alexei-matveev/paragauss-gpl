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
module  quadrupel_module
!---------------------------------------------------------------
!
!  Purpose: type definition and comm pack and unpack routines
!           for quadrupels
!           ( unique atom 1, l 1, unique atom 2, l 2 )
!           as used in the integral part.
!
!           A priority function ua_l_quadrupel_priority is also
!           defined, see below.
!
!           An order on quadrupels is defined by saying that
!           the later the quadrupel appears in the following
!           loop, the bigger it is:
!             do ua1
!                do l1
!                   do ua2
!                      do l2
!                         quadrupel(ua1,l1,ua2,l2) 
!                      enddo
!                   enddo
!                enddo
!             enddo
!           logical operators .gt. .ge. .lt. .le. .eq. are
!           defined.
!
!           Note the interpretation of l1/l2 for fitfunctions:
!             -1   : l = 0
!             0    : r**2
!             else : l = l1/l2
!
!  Module called by: quadrupelstack_module, int_distribute_module
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
! Modification 
! Author: AM
! Date:   01.07.99
! Description: quadrupel_transpose function added
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

# include "def.h"
  use type_module ! type specification parameters

  implicit none

  private         ! by default, all names are private
  save

!== Interrupt end of public interface of module =================


  !------------ Declaration of types ----------------------------
  type, public ::  quadrupel_type
     integer(kind=i4_kind) :: ua1, l1 , ua2, l2
  end type quadrupel_type



  !------------ Interface statements ----------------------------
  interface operator(.eq.)
     module procedure quadrupel_equal
  end interface
  public :: operator(.eq.)


  !------------ public functions and subroutines ----------------
  public :: quadrupel_priority
  public :: quadrupel_transpose
  public :: quadrupel_pack
  public :: quadrupel_unpack
  public :: quadrupel_write


  !===================================================================
  ! End of public interface of module
  !===================================================================

!------------ Declaration of private constants and variables ----

  !--------------------------------------------------------------
  !------------ Subroutines -------------------------------------
contains


  !**************************************************************
  subroutine quadrupel_pack(quadrupel)
    !  Purpose: comm packing of quadrupel
    !------------ Modules used ----------------------------------
    use comm_module
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type),intent(in)    :: quadrupel
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind)                :: info
    !------------ Executable code -------------------------------
    call commpack(quadrupel%ua1,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_pack: ua1")
    call commpack(quadrupel%l1,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_pack: l1")
    call commpack(quadrupel%ua2,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_pack: ua2")
    call commpack(quadrupel%l2,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_pack: l2")
  end subroutine quadrupel_pack
  !**************************************************************


  !**************************************************************
  subroutine quadrupel_unpack(quadrupel)
    !  Purpose: comm unpacking of quadrupel
    !------------ Modules used ----------------------------------
    use comm_module
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type),intent(out)    :: quadrupel
    !** End of interface ****************************************
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind)                :: info
    !------------ Executable code -------------------------------
    call communpack(quadrupel%ua1,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_unpack: ua1")
    call communpack(quadrupel%l1,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_unpack: l1")
    call communpack(quadrupel%ua2,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_unpack: ua2")
    call communpack(quadrupel%l2,info)
    if ( info .ne. 0 ) call error_handler("quadrupel_unpack: l2")
  end subroutine quadrupel_unpack
  !**************************************************************



  !**************************************************************
  logical function quadrupel_equal(q1,q2)
    !  Purpose: returns if q1 and q2 are equal
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type), intent(in) :: q1, q2
    !** End of interface ****************************************
    !------------ Executable code -------------------------------
    quadrupel_equal = q1%ua1 .eq. q2%ua1 .and. &
                           q1%ua2 .eq. q2%ua2 .and. &
                           q1%l1  .eq. q2%l1  .and. &
                           q1%l2  .eq. q2%l2
  end function quadrupel_equal
  !**************************************************************


  !mdf>>>
  !**************************************************************
  function quadrupel_transpose(q) result(tq)
    !  Purpose: exchange 2 and 1
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type), intent(in) :: q
    type(quadrupel_type)             :: tq
    !** End of interface ****************************************
    !------------ Executable code -------------------------------

    tq%ua2 = q%ua1
    tq%L2  = q%L1
    tq%ua1 = q%ua2
    tq%L1  = q%L2
  end function quadrupel_transpose
  !**************************************************************

  !**************************************************************
  subroutine quadrupel_write(q,output_unit,number,message,priority)
    !  Purpose: writes quadrupel q to output_unit
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type),            intent(in) :: q
    integer(kind=i4_kind),           intent(in) :: output_unit
    integer(kind=i4_kind), optional, intent(in) :: number
    real(kind=r8_kind),    optional, intent(in) :: priority
    character(len=*),      optional, intent(in) :: message
    !** End of interface ****************************************
    !------------ Executable code -------------------------------
    if ( present(number) .and. present(priority) &
         .and. present(message) ) then
       write(output_unit,'(I4,3X,A20,"  priority ",F10.0,4(3X,A3,I4))') &
            number, message, priority,&
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( present(number) .and. .not. present(priority) &
             .and. present(message) ) then
       write(output_unit,'(I4,3X,A20,21X,4(3X,A3,I4))') &
            number, message, &
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( present(number) .and. present(priority) &
             .and. .not. present(message) ) then
       write(output_unit,'(I4,23X,"  priority ",F10.0,4(3X,A3,I4))') &
            number, priority, &
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( present(number) .and. .not. present(priority) &
             .and. .not. present(message) ) then
       write(output_unit,'(I4,44X,4(3X,A3,I4))') &
            number, &
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( .not. present(number) .and. present(priority) &
             .and. present(message) ) then
       write(output_unit,'(7X,A20,"  priority ",F10.0,4(3X,A3,I4))') &
            message, priority,&
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( .not. present(number) .and. present(priority) &
             .and. .not. present(message) ) then
       write(output_unit,'(27X,"  priority ",F10.0,4(3X,A3,I4))') &
            priority,&
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( .not. present(number) .and. .not. present(priority) &
             .and. present(message) ) then
       write(output_unit,'(7X,A20,21X,4(3X,A3,I4))') &
            message,&
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    elseif ( .not. present(number) .and. present(message) ) then
       write(output_unit,'(A20,28X,4(3X,A3,I4))') &
            message, &
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    else
       write(output_unit,'(48X,4(3X,A3,I4))') &
            "ua1", q%ua1, "l1 ", q%l1, "ua2", q%ua2, "l2 ", q%l2
    end if
  end subroutine quadrupel_write
  !**************************************************************



  !**************************************************************
  real(kind=r8_kind) function quadrupel_priority(q,is_ff)
    !  Purpose: returns a priority number for given quadrupel
    !           depending upon size and l of basis sets involved.
    !           The priority number should be a meassure for the
    !           amount of work connected with the quadrupel:
    !           The higher the priority, the more work is involved.
    !           Thus, diagonal quadrupels with lowest priority
    !           should be calculated first and then non-diagonal
    !           quadrupels with highest priority.
    !------------ Modules used ----------------------------------
    use unique_atom_module
    implicit none
    !------------ Declaration of formal parameters -------------
    type(quadrupel_type), intent(in) :: q
    logical, optional, intent(in)    :: is_ff
      ! set .true. if you want priority for fitfunction basis,
      ! default: .false.
    !** End of interface ****************************************
    !------------ Declaration of formal parameters -------------
    logical                       :: fitfunction
    integer(kind=i4_kind)         :: l1_recursion, l2_recursion, &
         n1_exp, n2_exp, n_equiv_dist, n1_m, n2_m
    real(kind=r8_kind), parameter :: wight_recursion = 2.0_r8_kind, &
                                     wight_linear = 1.0_r8_kind
    !------------ Executable code ------------------------------

    if ( present(is_ff) ) then
       fitfunction = is_ff
    else
       fitfunction = .false.
    endif

    ! error check
    if ( q%ua1 .gt. N_unique_atoms .or.  q%ua1 .lt. 1 ) &
         call error_handler("quadrupel_priotity: wrong ua1")
    if ( q%ua2 .gt. N_unique_atoms .or.  q%ua2 .lt. 1 ) &
         call error_handler("quadrupel_priotity: wrong ua2")

    ! number of different sym.equiv. distances
    n_equiv_dist = unique_atom_symequiv(q%ua1,q%ua2)%n_equiv_dist

    ! recursion depth for primitive integrals
    if (q%l1 .eq. -1 ) then
       l1_recursion = 1
    else
       l1_recursion = q%l1 + 1
    endif
    if (q%l2 .eq. -1 ) then
       l2_recursion = 1
    else
       l2_recursion = q%l2 + 1
    endif

    ! number of magnetical quantum numbers
    n1_m = 2 * l1_recursion - 1
    n2_m = 2 * l2_recursion - 1
    
    ! number of exponents
    if (fitfunction ) then
       if (q%l1 .eq. -1 ) then
          n1_exp = &
               unique_atoms(q%ua1)%l_ch(0)%N_exponents + &
               unique_atoms(q%ua1)%l_xc(0)%N_exponents
       elseif (q%l1 .eq. 0 ) then
          n1_exp = &
               unique_atoms(q%ua1)%r2_ch%N_exponents + &
               unique_atoms(q%ua1)%r2_xc%N_exponents
       elseif ( ( q%l1 .gt. unique_atoms(q%ua1)%lmax_ch .and. &
                  q%l1 .gt. unique_atoms(q%ua1)%lmax_xc ) &
                 .or. q%l1 .lt. -1 ) then
          call error_handler("quadrupel_priotity: wrong l1")
       else
          if ( q%l1 .le. unique_atoms(q%ua1)%lmax_ch ) then
             n1_exp = unique_atoms(q%ua1)%l_ch(q%l1)%N_exponents
          else
             n1_exp = 0
          endif
          if ( q%l1 .le. unique_atoms(q%ua1)%lmax_xc ) then
             n1_exp = n1_exp + unique_atoms(q%ua1)%l_xc(q%l1)%N_exponents
          endif
       endif

       if (q%l2 .eq. -1 ) then
          n2_exp = &
               unique_atoms(q%ua2)%l_ch(0)%N_exponents + &
               unique_atoms(q%ua2)%l_xc(0)%N_exponents
       elseif (q%l2 .eq. -1 ) then
          n2_exp = &
               unique_atoms(q%ua2)%r2_ch%N_exponents + &
               unique_atoms(q%ua2)%r2_xc%N_exponents
       elseif ( ( q%l2 .gt. unique_atoms(q%ua2)%lmax_ch .and. &
                  q%l2 .gt. unique_atoms(q%ua2)%lmax_xc ) &
                 .or. q%l2 .lt. -1 ) then
          call error_handler("quadrupel_priotity: wrong l2")
       else
          if ( q%l2 .le. unique_atoms(q%ua2)%lmax_ch ) then
             n2_exp = unique_atoms(q%ua2)%l_ch(q%l2)%N_exponents
          else
             n2_exp = 0
          endif
          if ( q%l2 .le. unique_atoms(q%ua2)%lmax_xc ) then
             n2_exp = n2_exp + unique_atoms(q%ua2)%l_xc(q%l2)%N_exponents
          endif
       endif
    else ! orbitals
       n1_exp = unique_atoms(q%ua1)%l_ob(q%l1)%N_exponents
       n2_exp = unique_atoms(q%ua2)%l_ob(q%l2)%N_exponents
    endif

    quadrupel_priority = n_equiv_dist * &
         ( l1_recursion * wight_recursion + wight_linear ) * n1_exp * n1_m * &
         ( l2_recursion * wight_recursion + wight_linear ) * n2_exp * n2_m

  end function quadrupel_priority
  !**************************************************************



end module quadrupel_module

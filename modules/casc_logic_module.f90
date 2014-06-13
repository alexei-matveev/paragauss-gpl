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
module casc_logic_module
  !
  ! Purpose: Contains routines which provide steering logic
  !   for cascadic summation of data.
  ! 
  ! References: ~ttfs/DOCUMENTATION/DOK/kaskade.html
  !
  ! Module called by: ham_calc_module, xc_hamiltonian_module
  !
  ! Author: TG
  ! Date: 12/96
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
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/98
  ! Description: ...
  !
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------
  use type_module
  use comm_module, only: comm_get_n_processors

  implicit none
  private
  save

  !== Interrupt end of public interface of module =================
  
  !------------ public functions and subroutines ------------------
  public :: comp_rpn, comp_vpn

!================================================================
! End of public interface of module
!================================================================

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  !*************************************************************
  subroutine comp_rpn(rpn)
    !  Purpose: provides the important information for
    !           virtual processor number (VPN) to
    !           real processor number (RPN) mapping.
    !    
    !  Output parameter:
    !  rpn(2*scomm%n_procs-1)     rpn(i) is the RPN to vpn i
    !
    !  called by: cascadic_ham_calc
    !
    !  Author: TG, 9/96
    !-------------------------------------------------------------
    implicit none
    !------------ Declaration of formal parameters ----------------
    integer(kind=i4_kind) :: rpn(:)
    !** End of interface *****************************************  
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: i,n
    !------------------ Executable code ---------------------------
    n=comm_get_n_processors()
    rpn(1)=1
    do i=2,2*n-1
       if(mod(i,2).eq.0_i4_kind) then
          rpn(i)=rpn(i/2_i4_kind)
       else
          rpn(i)=(i+1_i4_kind)/2_i4_kind
       endif
    enddo
  end subroutine comp_rpn



  subroutine comp_vpn(vpn)
    !  Purpose: provides the important information for
    !           real processor number (RPN) to
    !           virtual processor number (VPN) mapping.
    !       
    !  Output parameter:
    !  vpn(scomm%n_procs)         vpn(i) is the highest VPN to RPN i
    !
    !  called by: build_hamiltonian
    !
    !  Author: TG, 9/96
    !---------------------------------------------------------------
    implicit none
    !------------ Declaration of formal parameters ----------------
    integer(kind=i4_kind) :: vpn(:)
    !** End of interface *****************************************  
    !------------ Declaration of local variables --------------------
    integer(kind=i4_kind) :: i,n
    integer(kind=i4_kind),allocatable :: rpn(:)
    !------------------ Executable code ---------------------------
    n=comm_get_n_processors()
    allocate(rpn(2*n-1))
    rpn(1)=1
    do i=2,2*n-1
       if(mod(i,2).eq.0_i4_kind) then
          rpn(i)=rpn(i/2_i4_kind)
       else
          rpn(i)=(i+1_i4_kind)/2_i4_kind
       endif
    enddo
    do i=1,2*n-1
       vpn(rpn(i))=i
    enddo
    deallocate(rpn)
  end subroutine comp_vpn
  !*************************************************************

!--------------- End of module ----------------------------------
end module casc_logic_module

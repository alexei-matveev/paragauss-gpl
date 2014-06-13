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
! Public interface of module
!================================================================
module symmetry_element
!-------------------------------------------------------------------------
!
!  Module stores symmetrieelements for (nearly) all pointgroups
!
!  The symmetry elements are stored as follows:
!  For every rotation the axis and the order of axis is stored.
!  For every mirror plane a normal vector on the plane is stored.
!  This datas are used in subroutine get_symmetry in grid_module.
! 
!  Author: MS
!  Date: 11/95
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
  use type_module
  use datatype

  implicit none
  private         ! by default, all names are private
  save

  !== Interrupt end of public interface of module =================


  type, public :: s_elements
     real(kind=r8_kind),dimension(:,:),pointer :: rotaxis,sigma
     integer(kind=i4_kind),dimension(:),pointer :: n_axis  !order of rotaxis
     integer(kind=i4_kind)   ::  num_axis,num_sigma
     logical  :: inversion
  end type s_elements

  type(s_elements), public :: sym_element, sym_element_local
  logical, private         :: allocated_sym_element_local = .false. ! SAVE by default

  !------------ public functions and subroutines ------------------
  public set_sym_element


  !================================================================
  ! End of public interface of module
  !================================================================

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine set_sym_element(ptgrp)
    implicit none
    character(len=4),intent(in)  :: ptgrp
    !** End of interface *****************************************


    real(kind=r8_kind),parameter :: pi=3.1415926535897932368_r8_kind
    real(kind=r8_kind),parameter :: g2r=pi/180.0_r8_kind

    real(kind=r8_kind) :: sin60,cos60,sin144,cos144,sin72,cos72,sqrt3m1,&
         sqrt2m1,zero,one,sin30,cos30,sin225,cos225,cos120,sin120
    
    integer :: allocstat

    sin30  = sin(g2r*30.0_r8_kind)
    cos30  = cos(g2r*30.0_r8_kind)
    sin60  = sin(g2r*60.0_r8_kind)
    cos60  = cos(g2r*60.0_r8_kind)
    cos72  = cos(g2r*72.0_r8_kind)
    sin72  = sin(g2r*72.0_r8_kind)
    cos144 = cos(g2r*144.0_r8_kind)
    sin144 = sin(g2r*144.0_r8_kind)
    cos225 = cos(g2r*22.5_r8_kind)
    sin225 = sin(g2r*22.5_r8_kind)
    sin120 = sin(g2r*120.0_r8_kind)
    cos120 = cos(g2r*120.0_r8_kind)
    zero=0.0_r8_kind
    one=1.0_r8_kind
    sqrt2m1=1/sqrt(2.0_r8_kind)
    sqrt3m1=1/sqrt(3.0_r8_kind)

    sym_element%inversion=.false.
    allocstat = 0

    select case(ptgrp)

    case('C1')
      if(.not.associated(sym_element%rotaxis)) then
       allocate(sym_element%rotaxis(3,0),sym_element%sigma(3,0), &
            STAT=allocstat)
       ASSERT(allocstat.eq.0)
      endif

    case('C1H')
       allocate(sym_element%rotaxis(3,0),sym_element%sigma(3,1), &
            STAT=allocstat)          
       sym_element%sigma(:,1)=(/zero,zero,one/)
      

    case('C2')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler &
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=2

    case('C2H')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,1), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=2
       sym_element%sigma(:,1)=(/zero,zero,one/)
       sym_element%inversion=.true.

    case('S6')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=3
       sym_element%inversion=.true.

    case('C3')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=3

    case('C3H')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,1), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=3
       sym_element%sigma(:,1)=(/zero,zero,one/)

    case('C4')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=4

    case('C4H')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,1), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=4
       sym_element%sigma(:,1)=(/zero,zero,one/)
       sym_element%inversion=.true.

    case('C5')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=5

    case('C5H')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,1), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=5
       sym_element%sigma(:,1)=(/zero,zero,one/)

    case('C6')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=6

    case('C6H')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,1), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=6
       sym_element%sigma(:,1)=(/zero,zero,one/)
       sym_element%inversion=.true.

    case('C7')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU &set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=7

    case('C8')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,0), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=8

    case('C8H')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,1), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/0.0_r8_kind,0.0_r8_kind,1.0_r8_kind/)
       sym_element%n_axis(1)=8
       sym_element%sigma(:,1)=(/zero,zero,one/)
       sym_element%inversion=.true.

    case('C2V')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,2), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=2
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/one,zero,zero/)

    case('C3V')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,3), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=3
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/sin60,cos60,zero/)
       sym_element%sigma(:,3)=(/sin60,-cos60,zero/)
       
    case('C4V')
       if(.not.associated(sym_element%rotaxis)) then
        allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,4), &
            sym_element%n_axis(1),STAT=allocstat)
        if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       endif

       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=4
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/one,zero,zero/)
       sym_element%sigma(:,3)=(/sqrt2m1,sqrt2m1,zero/)
       sym_element%sigma(:,4)=(/sqrt2m1,-sqrt2m1,zero/)

    case('C5V')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,5), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=5
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/sin144,cos144,zero/)
       sym_element%sigma(:,3)=(/-sin144,cos144,zero/)
       sym_element%sigma(:,4)=(/sin72,cos72,zero/)
       sym_element%sigma(:,5)=(/-sin72,cos72,zero/)


    case('C6V')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,6), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=6
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/one,zero,zero/)
       sym_element%sigma(:,3)=(/sin30,cos30,zero/)
       sym_element%sigma(:,4)=(/sin60,cos60,zero/)
       sym_element%sigma(:,5)=(/sin30,-cos30,zero/)
       sym_element%sigma(:,6)=(/sin60,-cos60,zero/)
!!$       sym_element%sigma(:,7)=(/sin225,-cos225,zero/)
!!$       sym_element%sigma(:,8)=(/-sin225,-cos225,zero/)



    case('C8V')
       allocate(sym_element%rotaxis(3,1),sym_element%sigma(3,8), &
            sym_element%n_axis(1),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in SU set_sym_element')
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=8
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/one,zero,zero/)
       sym_element%sigma(:,3)=(/sqrt2m1,sqrt2m1,zero/)
       sym_element%sigma(:,4)=(/sqrt2m1,-sqrt2m1,zero/)
       sym_element%sigma(:,5)=(/sin225,cos225,zero/)
       sym_element%sigma(:,6)=(/-sin225,cos225,zero/)
       sym_element%sigma(:,7)=(/sin225,-cos225,zero/)
       sym_element%sigma(:,8)=(/-sin225,-cos225,zero/)

    case('D2','D2H')
       if(ptgrp=='D2') then
          allocate(sym_element%rotaxis(3,3),sym_element%sigma(3,0), &
               sym_element%n_axis(3),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       else 
          allocate(sym_element%rotaxis(3,3),sym_element%sigma(3,3), &
               sym_element%n_axis(3),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       endif
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=2
       sym_element%rotaxis(:,2)=(/zero,one,zero/)
       sym_element%n_axis(2)=2
       sym_element%rotaxis(:,3)=(/one,zero,zero/)
       sym_element%n_axis(3)=2 
       if (ptgrp=='D2H') then
          sym_element%sigma(:,1)=(/zero,one,zero/)
          sym_element%sigma(:,2)=(/one,zero,zero/)
          sym_element%sigma(:,3)=(/one,one,zero/)
       endif

    case ('D2D','D4D','D3D','D5D','D6D') 

    case('D3','D3H')
       if (ptgrp=='D3') then
          allocate(sym_element%rotaxis(3,4),sym_element%sigma(3,0), &
               sym_element%n_axis(4),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       else
          allocate(sym_element%rotaxis(3,4),sym_element%sigma(3,4), &
               sym_element%n_axis(4),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       endif
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=3
       sym_element%rotaxis(:,2)=(/one,zero,zero/)
       sym_element%n_axis(2)=2
       sym_element%rotaxis(:,3)=(/cos120,sin120,zero/)
       sym_element%n_axis(3)=2
       sym_element%rotaxis(:,4)=(/cos120,-sin120,zero/)
       sym_element%n_axis(4)=2
       if (ptgrp=='D3H') then
          sym_element%sigma(:,1)=(/zero,one,zero/)
          sym_element%sigma(:,2)=(/sin60,cos60,zero/)
          sym_element%sigma(:,3)=(/sin60,-cos60,zero/)
          sym_element%sigma(:,4)=(/zero,zero,one/)
       endif

    case('D4','D4H')
       if (ptgrp=='D4') then
          allocate(sym_element%rotaxis(3,5),sym_element%sigma(3,0), &
               sym_element%n_axis(5),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       else
          allocate(sym_element%rotaxis(3,5),sym_element%sigma(3,5), &
               sym_element%n_axis(5),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       endif
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=4
       sym_element%rotaxis(:,2)=(/one,zero,zero/)
       sym_element%n_axis(2)=2
       sym_element%rotaxis(:,3)=(/zero,one,zero/)
       sym_element%n_axis(3)=2
       sym_element%rotaxis(:,4)=(/sqrt2m1,sqrt2m1,zero/)
       sym_element%n_axis(4)=2
       sym_element%rotaxis(:,5)=(/-sqrt2m1,-sqrt2m1,zero/)
       sym_element%n_axis(5)=2
       if(ptgrp=='D4H') then
          sym_element%sigma(:,1)=(/zero,one,zero/)
          sym_element%sigma(:,2)=(/one,zero,zero/)
          sym_element%sigma(:,3)=(/sqrt2m1,sqrt2m1,zero/)
          sym_element%sigma(:,4)=(/sqrt2m1,-sqrt2m1,zero/)
          sym_element%sigma(:,5)=(/zero,zero,one/)
       endif

    case('D5','D5H')
       if (ptgrp=='D5') then
          allocate(sym_element%rotaxis(3,6),sym_element%sigma(3,0), &
               sym_element%n_axis(6),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       else
          allocate(sym_element%rotaxis(3,6),sym_element%sigma(3,6), &
               sym_element%n_axis(6),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       endif
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=5
       sym_element%rotaxis(:,2)=(/one,zero,zero/)
       sym_element%n_axis(2)=2
       sym_element%rotaxis(:,3)=(/cos72,sin72,zero/)
       sym_element%n_axis(3)=2
       sym_element%rotaxis(:,4)=(/cos72,-sin72,zero/)
       sym_element%n_axis(4)=2
       sym_element%rotaxis(:,5)=(/cos144,sin144,zero/)
       sym_element%n_axis(5)=2
       sym_element%rotaxis(:,6)=(/cos144,-sin144,zero/)
       sym_element%n_axis(6)=2
       if (ptgrp=='D5H') then
          sym_element%sigma(:,1)=(/zero,one,zero/)
          sym_element%sigma(:,2)=(/sin144,cos144,zero/)
          sym_element%sigma(:,3)=(/-sin144,cos144,zero/)
          sym_element%sigma(:,4)=(/sin72,cos72,zero/)
          sym_element%sigma(:,5)=(/-sin72,cos72,zero/)
          sym_element%sigma(:,6)=(/zero,zero,one/)
       endif

    case('D6','D6H')

       if(ptgrp=='D6') then

          allocate(sym_element%rotaxis(3,7),sym_element%sigma(3,0), &
               sym_element%n_axis(7),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       else
          allocate(sym_element%rotaxis(3,7),sym_element%sigma(3,7), &
               sym_element%n_axis(7),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       endif
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=6
       sym_element%rotaxis(:,2)=(/one,zero,zero/)
       sym_element%n_axis(2)=2
       sym_element%rotaxis(:,3)=(/zero,one,zero/)
       sym_element%n_axis(3)=2
       sym_element%rotaxis(:,4)=(/cos30,sin30,zero/)
       sym_element%n_axis(4)=2
       sym_element%rotaxis(:,5)=(/cos30,-sin30,zero/)
       sym_element%n_axis(5)=2
       sym_element%rotaxis(:,6)=(/cos60,sin60,zero/)
       sym_element%n_axis(6)=2
       sym_element%rotaxis(:,7)=(/cos60,-sin60,zero/)
       sym_element%n_axis(7)=2

       if (ptgrp=='D6H') then
       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/one,zero,zero/)
       sym_element%sigma(:,3)=(/sin30,cos30,zero/)
       sym_element%sigma(:,4)=(/sin60,cos60,zero/)
       sym_element%sigma(:,5)=(/sin30,-cos30,zero/)
       sym_element%sigma(:,6)=(/sin60,-cos60,zero/)
       sym_element%sigma(:,7)=(/zero,zero,one/)
       endif

    case('D8','D8H')
       if(ptgrp=='D8') then
          allocate(sym_element%rotaxis(3,9),sym_element%sigma(3,0), &
               sym_element%n_axis(9),STAT=allocstat)
          if(allocstat/=0) call error_handler&
               ('Allocation failed in SU set_sym_element')
       else
          allocate(sym_element%rotaxis(3,9),sym_element%sigma(3,9), &
               sym_element%n_axis(9),STAT=allocstat)
          if(allocstat/=0) call  error_handler&
               ('Allocation failed in SU set_sym_element')          
       endif
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=8
       sym_element%rotaxis(:,2)=(/one,zero,zero/)
       sym_element%n_axis(2)=2
       sym_element%rotaxis(:,3)=(/zero,one,zero/)
       sym_element%n_axis(3)=2
       sym_element%rotaxis(:,4)=(/sqrt2m1,sqrt2m1,zero/)
       sym_element%n_axis(4)=2
       sym_element%rotaxis(:,5)=(/sqrt2m1,-sqrt2m1,zero/)
       sym_element%n_axis(5)=2
       sym_element%rotaxis(:,6)=(/cos225,sin225,zero/)
       sym_element%n_axis(6)=2
       sym_element%rotaxis(:,7)=(/cos225,-sin225,zero/)
       sym_element%n_axis(7)=2
       sym_element%rotaxis(:,8)=(/-cos225,sin225,zero/)
       sym_element%n_axis(8)=2
       sym_element%rotaxis(:,9)=(/sin225,cos225,zero/)
       sym_element%n_axis(9)=2
       if(ptgrp=='D8H') then
          sym_element%sigma(:,1)=(/zero,one,zero/)
          sym_element%sigma(:,2)=(/one,zero,zero/)
          sym_element%sigma(:,3)=(/sqrt2m1,sqrt2m1,zero/)
          sym_element%sigma(:,4)=(/sqrt2m1,-sqrt2m1,zero/)
          sym_element%sigma(:,5)=(/sin225,cos225,zero/)
          sym_element%sigma(:,6)=(/-sin225,cos225,zero/)
          sym_element%sigma(:,7)=(/sin225,-cos225,zero/)
          sym_element%sigma(:,8)=(/-sin225,-cos225,zero/)
          sym_element%sigma(:,9)=(/zero,zero,one/)
       endif

    case('IH')

    case('CS')

    case('TD')

    case('OH')
       allocate(sym_element%rotaxis(3,7),sym_element%sigma(3,9), &
            sym_element%n_axis(7),STAT=allocstat)
       if(allocstat/=0) call error_handler&
            ('Allocation failed in su set_sym_element')
       sym_element%rotaxis(:,1)=(/zero,zero,one/)
       sym_element%n_axis(1)=4
       sym_element%rotaxis(:,2)=(/one,zero,zero/)
       sym_element%n_axis(2)=4
       sym_element%rotaxis(:,3)=(/zero,one,zero/)
       sym_element%n_axis(3)=4
       sym_element%rotaxis(:,4)=(/sqrt3m1,sqrt3m1,sqrt3m1/)
       sym_element%n_axis(4)=3
       sym_element%rotaxis(:,5)=(/-sqrt3m1,sqrt3m1,sqrt3m1/)
       sym_element%n_axis(5)=3
       sym_element%rotaxis(:,6)=(/sqrt3m1,-sqrt3m1,sqrt3m1/)
       sym_element%n_axis(6)=3
       sym_element%rotaxis(:,7)=(/sqrt3m1,sqrt3m1,-sqrt3m1/)
       sym_element%n_axis(7)=3

       sym_element%sigma(:,1)=(/zero,one,zero/)
       sym_element%sigma(:,2)=(/one,zero,zero/)
       sym_element%sigma(:,3)=(/zero,zero,one/)
       sym_element%sigma(:,4)=(/sqrt2m1,sqrt2m1,zero/)
       sym_element%sigma(:,5)=(/sqrt2m1,-sqrt2m1,zero/)
       sym_element%sigma(:,6)=(/sqrt2m1,zero,sqrt2m1/)
       sym_element%sigma(:,7)=(/sqrt2m1,zero,-sqrt2m1/)
       sym_element%sigma(:,8)=(/zero,sqrt2m1,sqrt2m1/)
       sym_element%sigma(:,9)=(/zero,sqrt2m1,-sqrt2m1/)
       
    case default
      write(*,*) 'set_sym_element: warning, no symmetry information available for the desired pointgroup'

    end select
    if(allocstat/=0) then
       call error_handler('Allocation failed in su get_symmetrie_element')
    endif

    if( .not. allocated_sym_element_local )then
      ! this is used in grid_module to store (temporary?) local symmetries of 
      ! atomic sites:
      allocate( sym_element_local%rotaxis(3,10) &
              , sym_element_local%sigma(3,10)   &
              , sym_element_local%n_axis(10)    &
              , STAT=allocstat)
      ASSERT(allocstat==0)
      allocated_sym_element_local = .true.
    endif

  end subroutine set_sym_element

end module symmetry_element

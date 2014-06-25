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
module shgi_common
  !---------------------------------------------------------------
  !
  ! Holds global variables used in several modules.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
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
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  PUBLIC
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  real(RK)              :: XD(3), XD2

  real(RK), allocatable :: W(:,:), WC(:,:)             ! (NAB,3)
  real(RK), allocatable :: WDA(:), WDB(:), W2(:)       ! (NAB)
  real(RK), allocatable :: WDAL(:,:), WDBL(:,:)        ! (NAB,1+LA/B) for WDA**L
  logical,  allocatable :: CUTOFF(:,:)                 ! (NEA,NEB)
  real(RK), allocatable :: LAMBDA(:), ZETA(:)          ! (NAB)
  real(RK), allocatable :: NORM(:,:)                   ! (NAB,3)
  real(RK), allocatable :: LAMBDACH(:,:), NORMCH(:,:)  ! (NAB,NECS )
  real(RK), allocatable :: LAMBDAR2(:,:), NORMR2(:,:)  ! (NAB,NECR2)
  real(RK), allocatable :: F132R2(:,:)                 ! (NAB,NECR2)

  real(RK), allocatable :: X5(            :,        :,        :,        :,        :,       :)
  !                        X5(    (LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+SUM(L))

  real(RK), allocatable :: S5(  :,        :,        :,        :,        :,        :)
  !                        S5(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)

  real(RK), allocatable :: K4(  :,        :,        :,                  :,        :)
  !                        K4(NAB,(LA+1)**2,(LB+1)**2,          (LD+1)**2,(LE+1)**2)

  real(RK), allocatable :: YL(  :,        :,        :,        :,        :,        :,       :)
  !                        YL(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+SUM(L))

  real(RK), allocatable :: YS(:,:,:,:,:)               ! (NAB,2*LA+1,2*LB+1,1+LA+LB+LC,CXX:CYY)

  ! angular factor for GW and GD gradients:
  real(RK), allocatable :: GS(:,:,:,:,:)               ! (NAB,2*LA+1,2*LB+1,1+LA+LB+LC+LD,6+)

  real(RK)              :: XA(3), XB(3)                ! for PP
  real(RK), allocatable :: AEXP(:), BEXP(:)            ! (NEA), (NEB) a copy of exponents

  !--------------- End of module ----------------------------------
end module shgi_common

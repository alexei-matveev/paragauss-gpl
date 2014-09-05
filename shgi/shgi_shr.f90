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
module shgi_shr
  !-------------------------------------------------------------------
  !
  !  1) Computing harmonic derivatives of F(w^2/2), given tailor
  !     series FL for F:
  !
  !  * Example:
  !
  !      call SHR_D3Fv(N, LA, LB, LC, W, FL, D3F)
  !                    |-------  in -------| ^- out
  !
  !    executed with N vectors W(N, 3) and Tailor expansion FL(N,
  !    1+LA+LB+LC) gives triple harmonic derivatives wrt W:
  !
  !      D3F(:N,lma,lmb,lmc) = D(lma)D(lmb)D(lmc) F(W^2/2)
  !
  !    for all "lma", "lmb", and "lmc" up to respective limits
  !    (LA+1)^2, (LB+1)^2, and (LC+1)^2
  !
  !  * Call Tree:
  !
  !    SHR_DxFv           (x = 2, 3, 4, 5, 6)
  !       |----- SHR_DyFv (y = x - 1, called eventualy if some of
  !       |      angular momenta are zero )
  !       |----- XDxFv
  !                |------ SHR_DyC (precomputes derivatives of SH)
  !                |------ ZDxFv (actually computes harmonic derivatives)
  !                          |------ ZDyFv ( y = x - 1 )
  !                        .....
  !
  !  2) Computing the angular part of harmonic derivatives. If output
  !     is later convoluted with a tailor expansion for a specific
  !     F(w^2/2) the result is the harmonic derivative as given by
  !     "SHR_DxFv".
  !
  !  * Example
  !
  !      call SHR_D3Av(N, LA, LB, LC, W, D3A)
  !
  !    so that now one can apply it to some specific F:
  !
  !                          1+LA+LB+LC
  !                            ____
  !                            \
  !      D3F(:N,lma,lmb,lmc) = /    FL(:N,L) * D3A(:N,lma,lmb,lmc,L)
  !                           /___
  !                            L=1
  !
  !  * Call Tree:
  !
  !    SHR_DxAv           ( x = 2, 3, 4, 5 )
  !       |
  !      ... in complete analogy to SHR_DxFv.
  !
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
  ! Copyright (c) 2007-2011 Alexey Shor
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

!# define FPP_TIMERS 2
# include "def.h"
  use type_module, only: &
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  ! FIXME: move it to subs:
  use solhrules_module, only: &
       solhrules_differential, &
       solhrules_product,solhrules_l_and_m_of_lm
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: SHR_YLSv

  public :: SHR_D2Av
  public :: SHR_D3Av
  public :: SHR_D4Av
  public :: SHR_D5Av

  public :: SHR_X2Ps
  public :: SHR_PR2v
  public :: SHR_PR2Ls
  public :: SHR_doPR2Ls
  public :: SHR_doPR3Ls

  public :: SHR_D2Fv
  public :: SHR_D3Fv
  public :: SHR_D4Fv
  public :: SHR_D5Fv
  public :: SHR_D6Fv

  !======================================================
  ! These subs compute derivatives of spherical harmonics
  !
  !   DnC(:,lma,lmb,lmc,...) == D(lmb)D(lmc)... C(lma)
  !             \_    _/        \___    ___/
  !                \/                \/
  !              n-times           n-times
  !
  ! and used mainly internaly from within this module
  !
  public :: SHR_D1C !(N,LA,LB,            W,D1C)
  public :: SHR_D2C !(N,LA,LB,LC,         W,D2C)
  public :: SHR_D3C !(N,LA,LB,LC,LD,      W,D3C)
  public :: SHR_D4C !(N,LA,LB,LC,LD,LE,   W,D4C)
  public :: SHR_D5C !(N,LA,LB,LC,LD,LE,LF,W,D5C)


  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  integer(IK)   :: L_,M_
  integer(IK), parameter :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)
  integer(IK), parameter :: mof( (MAXL+1)**2 ) = (/((M_,M_=1,2*L_+1),L_=0,MAXL)/)
  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

   function SHR_X2Ps(LA,LB,LP,CA,X3A,CB,X3B) result(X2P)
     ! missing
     implicit none
     integer(IK), intent(in) :: LA,LB,LP
     real(RK)   , intent(in) :: CA(:), CB(:)           ! ( >= (LX+1)**2 ), X=A,B
     real(RK)   , intent(in) :: X3A(:,:,:), X3B(:,:,:) ! ((LA+1)**2,(LB+1)**2,1+LA+LB)
     real(RK)                :: X2P(2*LA+1,2*LB+1,1+LA,1+LB,1+LP)
     ! *** end of interface ***

     real(RK), parameter :: zero = 0.0_rk
#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: XX((LA+1)**2,(LB+1)**2)
#else
     real(RK), allocatable :: XX(:,:)
#endif
     integer(IK) :: ma,mb,lma,lmb
     integer(IK) :: sa,sb,lma1,lma2,lmb1,lmb2,la1,lb1
     integer(IK) :: ilp,mp,lmp
     integer(IK) :: LLA, LLB
     real(RK)    :: cf,cfa,cfb

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     allocate(XX((LA+1)**2,(LB+1)**2))
#endif

     ! FIXME: not all of it is used:
     X2P = zero

     do ilp=0,LP
        LLA = min(LA,LP)
        LLB = min(LB,LP)

        ! FIXME: not all of it is used:
        XX = zero
        do lmb=1,(LLB+1)**2
           do lma=1,(LLA+1)**2

              XX(lma,lmb) = zero
              ! sum over PP magnetic number:
              do mp=1,2*ilp+1
                 lmp = mp + (ilp)**2 ! index of {ilp,mp}
                 XX(lma,lmb) = XX(lma,lmb) + X3A(lma,lmp,1+ilp) * X3B(lmb,lmp,1+ilp)
              enddo
           enddo
        enddo

        do mb=1,2*LB+1
           lmb = mb + LB**2
           do ma=1,2*LA+1
              lma = ma + LA**2
              do sb=1,solhrules_product(lmb)%n_summands
                 cfb  = solhrules_product(lmb)%coef  (sb)
                 lmb1 = solhrules_product(lmb)%lm_sh1(sb)
                 lmb2 = solhrules_product(lmb)%lm_sh2(sb)
                 lb1  = lof(lmb1)
              do sa=1,solhrules_product(lma)%n_summands
                 cfa  = solhrules_product(lma)%coef  (sa)
                 lma1 = solhrules_product(lma)%lm_sh1(sa)
                 lma2 = solhrules_product(lma)%lm_sh2(sa)
                 la1  = lof(lma1)

                 cf = cfa * cfb
                 X2P(ma,mb,1+la1,1+lb1,1+ilp) = X2P(ma,mb,1+la1,1+lb1,1+ilp) &
                      + cf * CA(lma1) * CB(lmb1) * XX(lma2,lmb2)
              enddo
              enddo
           enddo
        enddo
     enddo ! LP

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(XX)
#endif

   end function SHR_X2Ps

   ! ***************************   D3A  *******************************
   !
   ! D3A subs compute angular factors
   !
   !   D3A(lma,lmb,lmc,L) L<=LA+LB+LC (1-off! -> L=0 at position 1)
   !
   ! of all TRIPLE derivatives
   !
   !   D(lma)D(lmb)D(lmc) S(x2/2) = SUM_L S^L(x2/2) * D3A(lma,lmb,lmc,L)
   !
   ! NOTE: D2A is a special case of D3A
   !   D3A(LA,LB,0) ~ D3A(0,LB,LC) ~ D3A(0,LB,LC) ~ D2A(L1,L2)
   !
   ! FIXME: it is highly symmetric in (lma,lmb,lmc) indices
   !        (because derivatives commute)
   !        This fact is not used yet in implementation

   subroutine ZD2Av(N,LA,LB,D1C,D2A)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB
     real(RK)   , intent(in)  :: D1C(:,:,:)
     ! D1C(N,>=(LA+1)**2,>=(LB+1)**2)
     real(RK)   , intent(out) :: &
          D2A(N,(LA+1)**2,(LB+1)**2, * ) ! * = 1:1+LA+LB
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,LA1
     integer(IK) :: lmb1,lmb2,LB1
     integer(IK) :: L,LLO,LHI
     integer(IK) :: sb
     real(RK)    :: cfb

     ! assume LA >= LB
     ASSERT(LA>=LB)
     ASSERT(N==size(D1C,1))
     ASSERT((LA+1)**2<=size(D1C,2))
     ASSERT((LB+1)**2<=size(D1C,3))

     ! FIXME: write and call ZD2Av()
     ! ================================================================
     ! D(lmb) D(lma) F =
     !              PR(lmb,lmb1,lmb2)
     !  ( D(lmb1) * F(LA) ) * ( D(lmb2) * C(lma) )
     !        \_____________________/

     D2A(:, :, :, :1 + LA + LB) = ZERO

     ! Compute symmetric D2A(:, lma, lmb, :)
     ! only for lma >= lmb:

     do lmb = 1, (LB + 1)**2 ! FIXME: skip L=0
     ! PR(lmb, lmb1, lmb2):
     do sb=1,  solhrules_product(lmb)%n_summands
        cfb  = solhrules_product(lmb)%coef(sb)
        lmb1 = solhrules_product(lmb)%lm_sh1(sb)
        lmb2 = solhrules_product(lmb)%lm_sh2(sb)
        LB1 = lof(lmb1)

     do lma = lmb, (LA + 1)**2    ! FIXME: skip L=0
        LA1 = lof(lma)
        ! in general (see below):
        !   LLO = LA1 + max(LB1)
        !   LHI = LA1 +     LB1
        !   do L=LLO,LHI
        !     ...
        ! but in this case only:
        L   = LA1 + LB1 ! so that L - LA1 == LB1
        D2A(:, lma, lmb, 1 + L) =  D2A(:, lma, lmb, 1 + L) &
             + cfb * D1C(:, lmb1, c00) * D1C(:, lma, lmb2)
        !                    \________________________/
     enddo
     enddo
     enddo

     ! ==================
     do lmb = 1, (LB + 1)**2
     do lma = 1, lmb - 1
              LA1 = lof(lma)
              LB1 = lof(lmb)
              LLO = 1 + max(LA1, LB1)
              LHI = 1 +     LA1 + LB1

        D2A(:, lma, lmb, LLO:LHI) = D2A(:, lmb, lma, LLO:LHI)
     enddo
     enddo
   end subroutine ZD2Av

   subroutine ZD2Fv(N,LA,LB,D1C,FL,F2)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB
     real(RK)   , intent(in)  :: D1C(:,:,:)
     ! D1C(N,>=(LA+1)**2,>=(LB+1)**2)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB)
     real(RK)   , intent(out) :: F2(N,(LA+1)**2,(LB+1)**2)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,LA1
     integer(IK) :: lmb1,lmb2,LB1
     integer(IK) :: L
     integer(IK) :: sb
     real(RK)    :: cfb

     ! assume LA >= LB
     ASSERT(LA>=LB)
     ASSERT(N==size(D1C,1))
     ASSERT((LA+1)**2<=size(D1C,2))
     ASSERT((LB+1)**2<=size(D1C,3))

     ! FIXME: write and call ZD2Av()
     ! ================================================================
     ! D(lmb) D(lma) F =
     !              PR(lmb,lmb1,lmb2)
     !  ( D(lmb1) * F(LA) ) * ( D(lmb2) * C(lma) )
     !        \_____________________/

     do lmb = 1, (LB + 1)**2 ! FIXME: skip L=0
        F2(:, :, lmb) = ZERO
     ! PR(lmb,lmb1,lmb2):
     do sb=1,  solhrules_product(lmb)%n_summands
        cfb  = solhrules_product(lmb)%coef(sb)
        lmb1 = solhrules_product(lmb)%lm_sh1(sb)
        lmb2 = solhrules_product(lmb)%lm_sh2(sb)
        LB1 = lof(lmb1)

     do lma = lmb, (LA + 1)**2    ! FIXME: skip L=0
        LA1 = lof(lma)
        ! in general (see below):
        !   LLO = LA1 + max(LB1)
        !   LHI = LA1 +     LB1
        !   do L=LLO,LHI
        !     ...
        ! but in this case only:
        L   = LA1 + LB1 ! so that L-LA1 == LB1
        F2(:, lma, lmb) =  F2(:, lma, lmb) + FL(:, 1 + L) &
             * cfb * D1C(:, lmb1, c00) * D1C(:, lma, lmb2)
        !                    \________________________/
     enddo
     enddo
     enddo

     ! ==================
     do lmb = 1, (LB + 1)**2
     do lma = 1, lmb - 1
        F2(:, lma, lmb) = F2(:, lmb, lma)
     enddo
     enddo
   end subroutine ZD2Fv

   subroutine XD2Av(N,LA,LB,W,D2A)
     ! prepare 2-derivatives of CLM
     ! then call ZD3Av to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: &
          D2A(N,(LA+1)**2,(LB+1)**2, * ) ! * = 1:1+SUM(L)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D1C(N,(LA+1)**2,(LB+1)**2)
#else
     real(RK), allocatable :: D1C(:,:,:)

     allocate(D1C(N,(LA+1)**2,(LB+1)**2))
#endif
     ! ================================================================
     ! precompute diffs of C(lm)
     ! D1C(lma,lmb) := D(lmb) * C(lma)
     call SHR_D1C(N,LA,LB,W,D1C)

     call ZD2Av(N,LA,LB,D1C,D2A)
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D1C)
#endif
   end subroutine XD2Av

   subroutine XD2Fv(N,LA,LB,W,FL,F2)
     ! prepare 2-derivatives of CLM
     ! then call ZD3Av to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB
     real(RK)   , intent(in)  :: W(N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB)
     real(RK)   , intent(out) :: F2(N,(LA+1)**2,(LB+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D1C(N,(LA+1)**2,(LB+1)**2)
#else
     real(RK), allocatable ::D1C(:,:,:)

     allocate(D1C(N,(LA+1)**2,(LB+1)**2))
#endif
     ! ================================================================
     ! precompute diffs of C(lm)
     ! D1C(lma,lmb) := D(lmb) * C(lma)
     call SHR_D1C(N,LA,LB,W,D1C)

     call ZD2Fv(N,LA,LB,D1C,FL,F2)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D1C)
#endif

   end subroutine XD2Fv

   subroutine ZD3Av(N,LA,LB,LC,D2C,D3A)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC
     real(RK)   , intent(in)  :: D2C(:,:,:,:)
     ! D2C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2)
     real(RK)   , intent(out) :: &
          D3A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2, * ) ! * = 1:1+LA+LB+LC
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc,LA1
     integer(IK) :: lmb1,lmb2,LB1,LB2
     integer(IK) :: lmc1,lmc2,LC1,LC2
     integer(IK) :: lm(3)
     integer(IK) :: L,LLO,LHI
     integer(IK) :: sb,sc
     real(RK)    :: cfb,cfc,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: X2A(N,(LB+1)**2,(LC+1)**2,1+LB+LC)
#else
     real(RK), allocatable :: X2A(:,:,:,:)

     allocate(X2A(N,(LB+1)**2,(LC+1)**2,1+LB+LC))
#endif

     FPP_TIMER_DECL(x2a)
     FPP_TIMER_DECL(x3a)

     ! assume LA >= LB >= LC
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(N==size(D2C,1))
     ASSERT((LA+1)**2<=size(D2C,2))
     ASSERT((LB+1)**2<=size(D2C,3))
     ASSERT((LC+1)**2<=size(D2C,4))

     FPP_TIMER_START(x2a)

     call ZD2Av(N,LB,LC,D2C(:,:,:,c00),X2A)

     FPP_TIMER_STOP(x2a)
     ! from now on:
     ! D3A(0,lmb,lmc,LLO:LHI) are (nonzero) ``Tailor'' coeffs of
     ! D3(0,lmb,lmc) * F. Here LLO=max(LC,LB), LHI=LC+LB

     FPP_TIMER_START(x3a)
     ! ================================================================
     ! D3(lma,lmb,lmc) F =
     !                     PR(lmb,lmb1,lmb2)
     !                     PR(lmc,lmc1,lmc2)
     !  ( D(lmc1) D(lmb1) * F(LA) ) * ( D(lmc2) D(lmb2) * C(lma) )
     !        \_______\____________________/       /
     !                 \__________________________/
     !
     !X3A = ZERO
     D3A(:,:,:,:,:1+LA+LB+LC) = ZERO

     ! Compute symmetric D3A(:,lma,lmb,lmc,:)
     ! only for lma >= lmb >= lmc:

        do lmc=1,(LC+1)**2 ! FIXME: skip L=0
        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC1  = lof(lmc1)
           LC2  = lof(lmc2)

        do lmb=lmc,(LB+1)**2 ! FIXME: skip L=0
        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB1  = lof(lmb1)
           LB2  = lof(lmb2)

           cf = cfb * cfc

!    do lma=1+(LB2+LC2)**2,(LA+1)**2 ! LA1 < LB2+LC2 => D(lmc2)D(lmb2)*C(lma) == 0
     do lma=lmb,(LA+1)**2
        LA1 = lof(lma)
        if( LA1 < LB2+LC2 ) cycle ! => D(lmc2)D(lmb2)*C(lma) == 0

           LLO  = LA1 + max(LB1,LC1) ! >= max(LA1,LB0,LC0)
           LHI  = LA1 +     LB1+LC1  ! <=     LA1+LB0+LC0
           do L = LLO,LHI            ! L-LC1 <= LB1+LA1
              D3A(:,lma,lmb,lmc,1+L) =  D3A(:,lma,lmb ,lmc ,1+L    ) &
                                 + cf * X2A(:,    lmb1,lmc1,1+L-LA1) * D2C(:,lma,lmb2,lmc2)
              !                                       \____\_______________________/    /
              !                                             \__________________________/
           enddo
     enddo    ! lma
        enddo ! sb
        enddo ! lmb
        enddo ! sc
        enddo ! lmc

     ! ==================
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc ) cycle
              LA1 = lof(lma)
              LB1 = lof(lmb)
              LC1 = lof(lmc)
              LLO = 1+ max(LA1,LB1,LC1)
              LHI = 1+     LA1+LB1+LC1

        lm = srt3(lma,lmb,lmc)
        D3A(:,lma,lmb,lmc,LLO:LHI) = D3A(:,lm(1),lm(2),lm(3),LLO:LHI)
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(X2A)
#endif

     FPP_TIMER_STOP(x3a)

     DPRINT 'XD3Av: X2A: ', FPP_TIMER_VALUE(x2a)
     DPRINT 'XD3Av: X3A: ', FPP_TIMER_VALUE(x3a)
   end subroutine ZD3Av

     function srt3(lm1,lm2,lm3) result(lm)
       ! returns reordered lm(1) >= lm(2) >= lm(3)
       implicit none
       integer(IK), intent(in) :: lm1,lm2,lm3 ! positive!
       integer(IK)             :: lm(3)
       ! *** end of interface ***

       integer(IK) :: jm(3), i, j , jx

       jm = (/ lm1, lm2, lm3 /)

       jx = 1        ! guess position of max
       do i=1,3
       do j=1,3
          if( jm(jx) < jm(j) ) jx = j
       enddo
          lm(i)  = jm(jx)
          jm(jx) = -1 ! erase current max
       enddo
     end function srt3

   subroutine ZD3Fv(N,LA,LB,LC,D2C,FL,F3)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC
     real(RK)   , intent(in)  :: D2C(:,:,:,:)
     ! D2C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC)
     real(RK)   , intent(out) :: F3(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc
     integer(IK) :: LA1,ma
     integer(IK) :: lmb1,lmb2,LB2
     integer(IK) :: lmc1,lmc2,LC2
     integer(IK) :: lm(3)
     integer(IK) :: sb,sc
     real(RK)    :: cfb,cfc,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: F2(N,(LB+1)**2,(LC+1)**2)
#else
     real(RK), allocatable :: F2(:,:,:)

     allocate(F2(N,(LB+1)**2,(LC+1)**2))
#endif

     ! assume LA >= LB >= LC
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(N==size(D2C,1))
     ASSERT((LA+1)**2<=size(D2C,2))
     ASSERT((LB+1)**2<=size(D2C,3))
     ASSERT((LC+1)**2<=size(D2C,4))

     ! ================================================================
     ! D3(lma,lmb,lmc) F =
     !                     PR(lmb,lmb1,lmb2)
     !                     PR(lmc,lmc1,lmc2)
     !  ( D(lmc1) D(lmb1) * F(LA) ) * ( D(lmc2) D(lmb2) * C(lma) )
     !        \_______\____________________/       /
     !                 \__________________________/
     !

     ! Compute symmetric F3(:,lma,lmb,lmc)
     ! only for lma >= lmb >= lmc:

     lma = 0
     do LA1=0,LA

        ! Recursive call to compute D(lmb)D(lmc) * F(LA1):
        call ZD2Fv(N,LB,LC,D2C(:,:,:,c00),FL(:,1+LA1:),F2)

     do ma=1,2*LA1+1
        lma = lma + 1
     do lmb=1,MIN((LB+1)**2,lma)
     do lmc=1,MIN((LC+1)**2,lmb)

        F3(:,lma,lmb,lmc) =  ZERO

        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB2  = lof(lmb2)

        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC2  = lof(lmc2)

           if( LA1 < LB2+LC2 ) cycle ! => D(lmc2)D(lmb2)*C(lma) == 0

           cf = cfb * cfc

           F3(:,lma,lmb,lmc) =  F3(:,lma,lmb,lmc) &
                         + cf * F2(:,lmb1,lmc1) * D2C(:,lma,lmb2,lmc2)
           !                            \____\_______________/    /
           !                                  \_________________/
        enddo ! sc
        enddo ! sb
     enddo ! lmc
     enddo ! lmb
     enddo !  ma
     enddo !  LA1

     ! ==================
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc ) cycle

        lm = srt3(lma,lmb,lmc)
        F3(:,lma,lmb,lmc) = F3(:,lm(1),lm(2),lm(3))
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(F2)
#endif

   end subroutine ZD3Fv

   subroutine XD3Av(N,LA,LB,LC,W,D3A)
     ! prepare 2-derivatives of CLM
     ! then call ZD3Av to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: &
          D3A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2, * ) ! * = 1:1+SUM(L)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D2C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
#else
     real(RK), allocatable :: D2C(:,:,:,:)

     allocate(D2C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2))
#endif

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd) := D(lmd) D(lmc) D(lmb) * C(lma)
     call SHR_D2C(N,LA,LB,LC,W,D2C)

     call ZD3Av(N,LA,LB,LC,D2C,D3A)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D2C)
#endif

   end subroutine XD3Av

   subroutine XD3Fv(N,LA,LB,LC,W,FL,F3)
     ! prepare 2-derivatives of CLM
     ! then call ZD3Fv to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC
     real(RK)   , intent(in)  :: W(N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC)
     real(RK)   , intent(out) :: F3(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D2C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
#else
     real(RK), allocatable :: D2C(:,:,:,:)

     allocate(D2C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2))
#endif
     ! ================================================================
     ! precompute diffs of C(lm)
     ! D2C(lma,lmb,lmc) := D(lmc) D(lmb) * C(lma)
     call SHR_D2C(N,LA,LB,LC,W,D2C)

     call ZD3Fv(N,LA,LB,LC,D2C,FL,F3)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D2C)
#endif
   end subroutine XD3Fv

   subroutine ZD4Av(N,LA,LB,LC,LD,D3C,D4A)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD
     real(RK)   , intent(in)  :: D3C(:,:,:,:,:)
     ! D3C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2,>=(LD+1)**2)
     real(RK)   , intent(out) :: &
          D4A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2, * ) ! * = 1:1+SUM(L)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc,lmd,LA1
     integer(IK) :: lmb1,lmb2,LB1,LB2
     integer(IK) :: lmc1,lmc2,LC1,LC2
     integer(IK) :: lmd1,lmd2,LD1,LD2
     integer(IK) :: lm(4)
     integer(IK) :: L,LLO,LHI
     integer(IK) :: sb,sc,sd
     real(RK)    :: cfb,cfc,cfd,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: X3A(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,1+LB+LC+LD)
#else
     real(RK), allocatable :: X3A(:,:,:,:,:)

     allocate(X3A(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,1+LB+LC+LD))
#endif
     ! FIXME: re-use D3A storage instead!

     ! assume LA >= LB >= LC >= LD
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(LC>=LD)
     ASSERT(N==size(D3C,1))
     ASSERT((LA+1)**2<=size(D3C,2))
     ASSERT((LB+1)**2<=size(D3C,3))
     ASSERT((LC+1)**2<=size(D3C,4))
     ASSERT((LD+1)**2<=size(D3C,5))

     call ZD3Av(N,LB,LC,LD,D3C(:,:,:,:,c00),X3A)

     ! ================================================================
     ! D3(lma,lmb,lmc) F =
     !                              PR(lmb,lmb1,lmb2)
     !                              PR(lmc,lmc1,lmc2)
     !                              PR(lmd,lmd1,lmd2)
     !  ( D(lmd1) D(lmc1) D(lmb1) * F(LA) ) * ( D(lmd2) D(lmc2) D(lmb2) * C(lma) )
     !        \_______\________\___________________/       /       /
     !                 \_________\_______________________/       /
     !                             \___________________________/
     ! FIXME: move lma loop inside? Do PR(b) and PR(c) sequentially

     D4A(:,:,:,:,:,:1+LA+LB+LC+LD) = ZERO

     ! Compute symmetric D4A(:,lma,lmb,lmc,lmd,:)
     ! only for lma >= lmb >= lmc >= lmd:

        ! PR(lmd,lmd1,lmd2):
        do lmd=1,(LD+1)**2 ! FIXME: skip L=0
        do sd=1,  solhrules_product(lmd)%n_summands
           cfd  = solhrules_product(lmd)%coef(sd)
           lmd1 = solhrules_product(lmd)%lm_sh1(sd)
           lmd2 = solhrules_product(lmd)%lm_sh2(sd)
           LD1  = lof(lmd1)
           LD2  = lof(lmd2)

        do lmc=lmd,(LC+1)**2 ! FIXME: skip L=0
        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC1  = lof(lmc1)
           LC2  = lof(lmc2)

        do lmb=lmc,(LB+1)**2 ! FIXME: skip L=0
        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB1  = lof(lmb1)
           LB2  = lof(lmb2)

           cf = cfb * cfc * cfd

!    do lma=1+(LB2+LC2+LD2)**2,(LA+1)**2 ! LA1 < LB2+LC2+LD2 => D(lmd2)D(lmc2)D(lmb2)*C(lma) == 0
     do lma=lmb,(LA+1)**2
        LA1 = lof(lma)
        if( LA1 < LB2+LC2+LD2 ) cycle ! => D(lmd2)D(lmc2)D(lmb2)*C(lma) == 0

           LLO  = LA1 + max(LB1,LC1,LD1) ! >= max(LA1,LB0,LC0,LD0)
           LHI  = LA1 +     LB1+LC1+LD1  ! <=     LA1+LB0+LC0+LD0
           do L = LLO,LHI                ! L-LC1 <= LB1+LA1+LD1
              D4A(:,lma,lmb,lmc,lmd,1+L) =  D4A(:,lma,lmb ,lmc ,lmd ,1+L    ) &
                                     + cf * X3A(:,    lmb1,lmc1,lmd1,1+L-LA1) * D3C(:,lma,lmb2,lmc2,lmd2)
              !                                         \____\_____\_______________________/    /    /
              !                                                \_____\________________________/    /
              !                                                        \_________________________/
           enddo
     enddo    ! lma
        enddo ! sb
        enddo ! lmb
        enddo ! sc
        enddo ! lmc
        enddo ! sd
        enddo ! lmd

     ! ==================
     do lmd=  1,(LD+1)**2
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc .and. lmc>=lmd ) cycle
              LA1 = lof(lma)
              LB1 = lof(lmb)
              LC1 = lof(lmc)
              LD1 = lof(lmd)
              LLO = 1+ max(LA1,LB1,LC1,LD1)
              LHI = 1+     LA1+LB1+LC1+LD1

        lm = srt4(lma,lmb,lmc,lmd)
        D4A(:,lma,lmb,lmc,lmd,LLO:LHI) = D4A(:,lm(1),lm(2),lm(3),lm(4),LLO:LHI)
     enddo
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(X3A)
#endif

   end subroutine ZD4Av

     function srt4(lm1,lm2,lm3,lm4) result(lm)
       ! returns reordered lm(1) >= lm(2) >= lm(3)
       implicit none
       integer(IK), intent(in) :: lm1,lm2,lm3,lm4 ! positive!
       integer(IK)             :: lm(4)
       ! *** end of interface ***

       integer(IK) :: jm(4), i, j , jx

       jm = (/ lm1, lm2, lm3, lm4 /)

       jx = 1        ! guess position of max
       do i=1,4
       do j=1,4
          if( jm(jx) < jm(j) ) jx = j
       enddo
          lm(i)  = jm(jx)
          jm(jx) = -1 ! erase current max
       enddo
     end function srt4

   subroutine ZD4Fv(N,LA,LB,LC,LD,D3C,FL,F4)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD
     real(RK)   , intent(in)  :: D3C(:,:,:,:,:)
     ! D2C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2,>=(LD+1)**2)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD)
     real(RK)   , intent(out) :: F4(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc,lmd
     integer(IK) :: LA1,ma
     integer(IK) :: lmb1,lmb2,LB2
     integer(IK) :: lmc1,lmc2,LC2
     integer(IK) :: lmd1,lmd2,LD2
     integer(IK) :: lm(4)
     integer(IK) :: sb,sc,sd
     real(RK)    :: cfb,cfc,cfd,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: F3(N,(LB+1)**2,(LC+1)**2,(LD+1)**2)
#else
     real(RK), allocatable :: F3(:,:,:,:)

     allocate(F3(N,(LB+1)**2,(LC+1)**2,(LD+1)**2))
#endif

     ! assume LA >= LB >= LC >= LD
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(LC>=LD)
     ASSERT(N==size(D3C,1))
     ASSERT((LA+1)**2<=size(D3C,2))
     ASSERT((LB+1)**2<=size(D3C,3))
     ASSERT((LC+1)**2<=size(D3C,4))
     ASSERT((LD+1)**2<=size(D3C,5))

     ! ================================================================
     ! D4(lma,lmb,lmc,lmd) F =
     !                     PR(lmb,lmb1,lmb2)
     !                     PR(lmc,lmc1,lmc2)
     !                     PR(lmd,lmd1,lmd2)
     !  ( D(lmd1) D(lmc1) D(lmb1) * F(LA) ) * ( D(lmd2) D(lmc2) D(lmb2) * C(lma) )
     !        \_______\_______\___________________/       /      /
     !                 \_______\________________________/      /
     !                          \____________________________/
     !

     ! Compute symmetric F4(:,lma,lmb,lmc,lmd)
     ! only for lma >= lmb >= lmc >= lmd:

     lma = 0
     do LA1=0,LA

        ! Recursive call to compute D(lmb)D(lmc)D(lmd) * F(LA1):
        call ZD3Fv(N,LB,LC,LD,D3C(:,:,:,:,c00),FL(:,1+LA1:),F3)

     do ma=1,2*LA1+1
        lma = lma + 1
     do lmb=1,MIN((LB+1)**2,lma)
     do lmc=1,MIN((LC+1)**2,lmb)
     do lmd=1,MIN((LD+1)**2,lmc)

        F4(:,lma,lmb,lmc,lmd) =  ZERO

        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB2  = lof(lmb2)

        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC2  = lof(lmc2)

        ! PR(lmd,lmd1,lmd2):
        do sd=1,  solhrules_product(lmd)%n_summands
           cfd  = solhrules_product(lmd)%coef(sd)
           lmd1 = solhrules_product(lmd)%lm_sh1(sd)
           lmd2 = solhrules_product(lmd)%lm_sh2(sd)
           LD2  = lof(lmd2)

           if( LA1 < LB2+LC2+LD2 ) cycle ! => D(lmd2)D(lmc2)D(lmb2)*C(lma) == 0

           cf = cfb * cfc * cfd

           F4(:,lma,lmb,lmc,lmd) =  F4(:,lma,lmb,lmc,lmd) &
                         + cf * F3(:,lmb1,lmc1,lmd1) * D3C(:,lma,lmb2,lmc2,lmd2)
           !                            \____\___\________________/    /    /
           !                                  \___\__________________/    /
           !                                       \____________________/
        enddo ! sd
        enddo ! sc
        enddo ! sb
     enddo ! lmd
     enddo ! lmc
     enddo ! lmb
     enddo !  ma
     enddo !  LA1

     ! ==================
     do lmd=  1,(LD+1)**2
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc .and. lmc>=lmd) cycle

        lm = srt4(lma,lmb,lmc,lmd)
        F4(:,lma,lmb,lmc,lmd) = F4(:,lm(1),lm(2),lm(3),lm(4))
     enddo
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(F3)
#endif

   end subroutine ZD4Fv

   subroutine XD4Av(N,LA,LB,LC,LD,W,D4A)
     ! prepare 3-derivatives of CLM
     ! then call ZD4Av to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD
     real(RK)   , intent(in)  :: W(N,3)   ! (N,3)
     real(RK)   , intent(out) :: &
          D4A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2, * ) ! * = 1:1+SUM(L)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D3C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
#else
     real(RK), allocatable :: D3C(:,:,:,:,:)

     allocate(D3C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2))
#endif

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd) := D(lmd) D(lmc) D(lmb) * C(lma)
     call SHR_D3C(N,LA,LB,LC,LD,W,D3C)

     call ZD4Av(N,LA,LB,LC,LD,D3C,D4A)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D3C)
#endif

   end subroutine XD4Av

   subroutine XD4Fv(N,LA,LB,LC,LD,W,FL,F4)
     ! prepare 4-derivatives of CLM
     ! then call ZD4Fv to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD
     real(RK)   , intent(in)  :: W(N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD)
     real(RK)   , intent(out) :: F4(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D3C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
#else
     real(RK), allocatable :: D3C(:,:,:,:,:)

     allocate(D3C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2))
#endif
     ! ================================================================
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd) := D(lmd) D(lmc) D(lmb) * C(lma)
     call SHR_D3C(N,LA,LB,LC,LD,W,D3C)

     call ZD4Fv(N,LA,LB,LC,LD,D3C,FL,F4)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D3C)
#endif

   end subroutine XD4Fv

   subroutine ZD5Av(N,LA,LB,LC,LD,LE,D4C,D5A)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE
     real(RK)   , intent(in)  :: D4C(:,:,:,:,:,:)
     ! D4C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2,>=(LD+1)**2,>=(LE+1)**2)
     real(RK)   , intent(out) :: &
          D5A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2, * ) ! * = 1:1+SUM(L)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc,lmd,lme,LA1
     integer(IK) :: lmb1,lmb2,LB1,LB2
     integer(IK) :: lmc1,lmc2,LC1,LC2
     integer(IK) :: lmd1,lmd2,LD1,LD2
     integer(IK) :: lme1,lme2,LE1,LE2
     integer(IK) :: lm(5)
     integer(IK) :: L,LLO,LHI
     integer(IK) :: sb,sc,sd,se
     real(RK)    :: cfb,cfc,cfd,cfe,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: X4A(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+LB+LC+LD+LE)
#else
     real(RK),allocatable :: X4A(:,:,:,:,:,:)

     allocate(X4A(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+LB+LC+LD+LE))
#endif
     ! FIXME: re-use D5A storage instead!

     ! assume LA >= LB >= LC >= LD
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(LC>=LD)
     ASSERT(LD>=LE)
     ASSERT(N==size(D4C,1))
     ASSERT((LA+1)**2<=size(D4C,2))
     ASSERT((LB+1)**2<=size(D4C,3))
     ASSERT((LC+1)**2<=size(D4C,4))
     ASSERT((LD+1)**2<=size(D4C,5))
     ASSERT((LE+1)**2<=size(D4C,6))

     call ZD4Av(N,LB,LC,LD,LE,D4C(:,:,:,:,:,c00),X4A)

     ! ================================================================
     ! D5(lma,lmb,lmc,lmd,lme) F =
     !                              PR(lmb,lmb1,lmb2)
     !                                    ...
     !                              PR(lme,lme1,lme2)
     !  ( D(lme1)   ...   D(lmb1) * F(LA) ) * ( D(lme2)   ...   D(lmb2) * C(lma) )
     !        \                \___________________/_______________/
     !         \__________________________________/
     !

     D5A(:,:,:,:,:,:,:1+LA+LB+LC+LD+LE) = ZERO

     ! Compute symmetric D5A(:,lma,lmb,lmc,lmd,lme,:)
     ! only for lma >= lmb >= lmc >= lmd >= lme:

        ! PR(lme,lme1,lme2):
        do lme=1,(LE+1)**2 ! FIXME: skip L=0
        do se=1,  solhrules_product(lme)%n_summands
           cfe  = solhrules_product(lme)%coef(se)
           lme1 = solhrules_product(lme)%lm_sh1(se)
           lme2 = solhrules_product(lme)%lm_sh2(se)
           LE1  = lof(lme1)
           LE2  = lof(lme2)

        ! PR(lmd,lmd1,lmd2):
        do lmd=lme,(LD+1)**2 ! FIXME: skip L=0
        do sd=1,  solhrules_product(lmd)%n_summands
           cfd  = solhrules_product(lmd)%coef(sd)
           lmd1 = solhrules_product(lmd)%lm_sh1(sd)
           lmd2 = solhrules_product(lmd)%lm_sh2(sd)
           LD1  = lof(lmd1)
           LD2  = lof(lmd2)

        do lmc=lmd,(LC+1)**2 ! FIXME: skip L=0
        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC1  = lof(lmc1)
           LC2  = lof(lmc2)

        do lmb=lmc,(LB+1)**2 ! FIXME: skip L=0
        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB1  = lof(lmb1)
           LB2  = lof(lmb2)

           cf = cfb * cfc * cfd * cfe

     do lma=lmb,(LA+1)**2
        LA1 = lof(lma)
        if( LA1 < LB2+LC2+LD2+LE2 ) cycle ! => D(lme2)...D(lmb2)*C(lma) == 0

           LLO  = LA1 + max(LB1,LC1,LD1,LE1) ! >= max(LA1,LB0,LC0,LD0,LE0)
           LHI  = LA1 +     LB1+LC1+LD1+LE1  ! <=     LA1+LB0+LC0+LD0+LE0
           do L = LLO,LHI
              D5A(:,lma,lmb,lmc,lmd,lme,1+L) =  D5A(:,lma,lmb ,lmc ,lmd ,lme ,1+L    ) &
                                         + cf * X4A(:,    lmb1,lmc1,lmd1,lme1,1+L-LA1) &
                                              * D4C(:,lma,lmb2,lmc2,lmd2,lme2)

           enddo
     enddo    ! lma
        enddo ! sb
        enddo ! lmb
        enddo ! sc
        enddo ! lmc
        enddo ! sd
        enddo ! lmd
        enddo ! se
        enddo ! lme

     ! ==================
     do lme=  1,(LE+1)**2
     do lmd=  1,(LD+1)**2
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc .and. lmc>=lmd .and. lmd>=lme ) cycle
              LA1 = lof(lma)
              LB1 = lof(lmb)
              LC1 = lof(lmc)
              LD1 = lof(lmd)
              LE1 = lof(lme)
              LLO = 1+ max(LA1,LB1,LC1,LD1,LE1)
              LHI = 1+     LA1+LB1+LC1+LD1+LE1

        lm = srt5(lma,lmb,lmc,lmd,lme)
        D5A(:,lma,lmb,lmc,lmd,lme,LLO:LHI) = D5A(:,lm(1),lm(2),lm(3),lm(4),lm(5),LLO:LHI)
     enddo
     enddo
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(X4A)
#endif

   end subroutine ZD5Av

     function srt5(lm1,lm2,lm3,lm4,lm5) result(lm)
       ! returns reordered lm(1) >= lm(2) >= lm(3)
       implicit none
       integer(IK), intent(in) :: lm1,lm2,lm3,lm4,lm5 ! positive!
       integer(IK)             :: lm(5)
       ! *** end of interface ***

       integer(IK) :: jm(5), i, j , jx

       jm = (/ lm1, lm2, lm3, lm4, lm5 /)

       jx = 1        ! guess position of max
       do i=1,5
       do j=1,5
          if( jm(jx) < jm(j) ) jx = j
       enddo
          lm(i)  = jm(jx)
          jm(jx) = -1 ! erase current max
       enddo
     end function srt5

     function srt6(lm1,lm2,lm3,lm4,lm5,lm6) result(lm)
       ! returns reordered lm(1) >= lm(2) >= lm(3)
       implicit none
       integer(IK), intent(in) :: lm1,lm2,lm3,lm4,lm5,lm6 ! positive!
       integer(IK)             :: lm(6)
       ! *** end of interface ***

       integer(IK) :: jm(6), i, j , jx

       jm = (/ lm1, lm2, lm3, lm4, lm5, lm6 /)

       jx = 1        ! guess position of max
       do i=1,6
       do j=1,6
          if( jm(jx) < jm(j) ) jx = j
       enddo
          lm(i)  = jm(jx)
          jm(jx) = -1 ! erase current max
       enddo
     end function srt6

   subroutine ZD5Fv(N,LA,LB,LC,LD,LE,D4C,FL,F5)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE
     real(RK)   , intent(in)  :: D4C(:,:,:,:,:,:)
     ! D2C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2,>=(LD+1)**2,>=(LE+1)**2)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD+LE)
     real(RK)   , intent(out) :: F5(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc,lmd,lme
     integer(IK) :: LA1,ma
     integer(IK) :: lmb1,lmb2,LB2
     integer(IK) :: lmc1,lmc2,LC2
     integer(IK) :: lmd1,lmd2,LD2
     integer(IK) :: lme1,lme2,LE2
     integer(IK) :: lm(5)
     integer(IK) :: sb,sc,sd,se
     real(RK)    :: cfb,cfc,cfd,cfe,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: F4(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
#else
     real(RK), allocatable :: F4(:,:,:,:,:)

     allocate(F4(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2))
#endif

     ! assume LA >= LB >= LC >= LD >= LE
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(LC>=LD)
     ASSERT(LD>=LE)
     ASSERT(N==size(D4C,1))
     ASSERT((LA+1)**2<=size(D4C,2))
     ASSERT((LB+1)**2<=size(D4C,3))
     ASSERT((LC+1)**2<=size(D4C,4))
     ASSERT((LD+1)**2<=size(D4C,5))
     ASSERT((LE+1)**2<=size(D4C,6))

     ! ================================================================
     ! D5(lma,lmb,lmc,lmd,lme) F =
     !                              PR(lmb,lmb1,lmb2)
     !                                    ...
     !                              PR(lme,lme1,lme2)
     !  ( D(lme1)   ...   D(lmb1) * F(LA) ) * ( D(lme2)   ...   D(lmb2) * C(lma) )
     !        \                \___________________/_______________/
     !         \__________________________________/
     !

     ! Compute symmetric F5(:,lma,lmb,lmc,lmd,lme)
     ! only for lma >= lmb >= lmc >= lmd >= lme:

     lma = 0
     do LA1=0,LA

        ! Recursive call to compute D(lmb)D(lmc)D(lmd)D(lme) * F(LA1):
        call ZD4Fv(N,LB,LC,LD,LE,D4C(:,:,:,:,:,c00),FL(:,1+LA1:),F4)

     do ma=1,2*LA1+1
        lma = lma + 1
     do lmb=1,MIN((LB+1)**2,lma)
     do lmc=1,MIN((LC+1)**2,lmb)
     do lmd=1,MIN((LD+1)**2,lmc)
     do lme=1,MIN((LE+1)**2,lmd)

        F5(:,lma,lmb,lmc,lmd,lme) =  ZERO

        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB2  = lof(lmb2)

        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC2  = lof(lmc2)

        ! PR(lmd,lmd1,lmd2):
        do sd=1,  solhrules_product(lmd)%n_summands
           cfd  = solhrules_product(lmd)%coef(sd)
           lmd1 = solhrules_product(lmd)%lm_sh1(sd)
           lmd2 = solhrules_product(lmd)%lm_sh2(sd)
           LD2  = lof(lmd2)

        ! PR(lme,lme1,lme2):
        do se=1,  solhrules_product(lme)%n_summands
           cfe  = solhrules_product(lme)%coef(se)
           lme1 = solhrules_product(lme)%lm_sh1(se)
           lme2 = solhrules_product(lme)%lm_sh2(se)
           LE2  = lof(lme2)

           if( LA1 < LB2+LC2+LD2+LE2 ) cycle ! => D(lme2)D(lmd2)D(lmc2)D(lmb2)*C(lma) == 0

           cf = cfb * cfc * cfd * cfe

           F5(:,lma,lmb,lmc,lmd,lme) =  F5(:,lma,lmb,lmc,lmd,lme) &
                         + cf * F4(:,lmb1,lmc1,lmd1,lme1) * D4C(:,lma,lmb2,lmc2,lmd2,lme2)
           !                            \____\___\____\____ prb _______/    /    /   /
           !                                  \___\____\____ prc _________/    /   /
           !                                       \____\____ prd ___________/   /
        enddo ! se                                       \____ pre ____________/
        enddo ! sd
        enddo ! sc
        enddo ! sb
     enddo ! lme
     enddo ! lmd
     enddo ! lmc
     enddo ! lmb
     enddo !  ma
     enddo !  LA1

     ! ==================
     do lme=  1,(LE+1)**2
     do lmd=  1,(LD+1)**2
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc .and. lmc>=lmd .and. lmd>=lme ) cycle

        lm = srt5(lma,lmb,lmc,lmd,lme)
        F5(:,lma,lmb,lmc,lmd,lme) = F5(:,lm(1),lm(2),lm(3),lm(4),lm(5))
     enddo
     enddo
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(F4)
#endif

   end subroutine ZD5Fv

   subroutine ZD6Fv(N,LA,LB,LC,LD,LE,LF,D5C,FL,F6)
     use constants, only: ZERO, ONE
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE,LF
     real(RK)   , intent(in)  :: D5C(:,:,:,:,:,:,:)
     ! D5C(N,>=(LA+1)**2,>=(LB+1)**2,>=(LC+1)**2,>=(LD+1)**2,>=(LE+1)**2,>=(LF+1)**2)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD+LE+LF)
     real(RK)   , intent(out) :: F6(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
     ! *** end of interface ***

     integer(IK), parameter :: c00 = 1
     integer(IK) :: lma,lmb,lmc,lmd,lme,lmf
     integer(IK) :: LA1,ma
     integer(IK) :: lmb1,lmb2,LB2
     integer(IK) :: lmc1,lmc2,LC2
     integer(IK) :: lmd1,lmd2,LD2
     integer(IK) :: lme1,lme2,LE2
     integer(IK) :: lmf1,lmf2,LF2
     integer(IK) :: lm(6)
     integer(IK) :: sb,sc,sd,se,sf
     real(RK)    :: cfb,cfc,cfd,cfe,cff,cf

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: F5(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
#else
     real(RK), allocatable :: F5(:,:,:,:,:,:)

     allocate(F5(N,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2))
#endif

     ! assume LA >= LB >= LC >= LD >= LE >= LF
     ASSERT(LA>=LB)
     ASSERT(LB>=LC)
     ASSERT(LC>=LD)
     ASSERT(LD>=LE)
     ASSERT(LE>=LF)
     ASSERT(N==size(D5C,1))
     ASSERT((LA+1)**2<=size(D5C,2))
     ASSERT((LB+1)**2<=size(D5C,3))
     ASSERT((LC+1)**2<=size(D5C,4))
     ASSERT((LD+1)**2<=size(D5C,5))
     ASSERT((LE+1)**2<=size(D5C,6))
     ASSERT((LF+1)**2<=size(D5C,7))

     ! ================================================================
     ! D5(lma,lmb,lmc,lmd,lme,lmf) F =
     !                              PR(lmb,lmb1,lmb2)
     !                                    ...
     !                              PR(lmf,lmf1,lmf2)
     !  ( D(lmf1)   ...   D(lmb1) * F(LA) ) * ( D(lmf2)   ...   D(lmb2) * C(lma) )
     !        \                \___________________/_______________/
     !         \__________________________________/
     !

     ! Compute symmetric F5(:,lma,lmb,lmc,lmd,lme,lmf)
     ! only for lma >= lmb >= lmc >= lmd >= lme >= lmf:

     lma = 0
     do LA1=0,LA

        ! Recursive call to compute D(lmb)D(lmc)D(lmd)D(lme)D(lmf) * F(LA1):
        call ZD5Fv(N,LB,LC,LD,LE,LF,D5C(:,:,:,:,:,:,c00),FL(:,1+LA1:),F5)

     do ma=1,2*LA1+1
        lma = lma + 1
     do lmb=1,MIN((LB+1)**2,lma)
     do lmc=1,MIN((LC+1)**2,lmb)
     do lmd=1,MIN((LD+1)**2,lmc)
     do lme=1,MIN((LE+1)**2,lmd)
     do lmf=1,MIN((LF+1)**2,lme)

        F6(:,lma,lmb,lmc,lmd,lme,lmf) =  ZERO

        ! PR(lmb,lmb1,lmb2):
        do sb=1,  solhrules_product(lmb)%n_summands
           cfb  = solhrules_product(lmb)%coef(sb)
           lmb1 = solhrules_product(lmb)%lm_sh1(sb)
           lmb2 = solhrules_product(lmb)%lm_sh2(sb)
           LB2  = lof(lmb2)

        ! PR(lmc,lmc1,lmc2):
        do sc=1,  solhrules_product(lmc)%n_summands
           cfc  = solhrules_product(lmc)%coef(sc)
           lmc1 = solhrules_product(lmc)%lm_sh1(sc)
           lmc2 = solhrules_product(lmc)%lm_sh2(sc)
           LC2  = lof(lmc2)

        ! PR(lmd,lmd1,lmd2):
        do sd=1,  solhrules_product(lmd)%n_summands
           cfd  = solhrules_product(lmd)%coef(sd)
           lmd1 = solhrules_product(lmd)%lm_sh1(sd)
           lmd2 = solhrules_product(lmd)%lm_sh2(sd)
           LD2  = lof(lmd2)

        ! PR(lme,lme1,lme2):
        do se=1,  solhrules_product(lme)%n_summands
           cfe  = solhrules_product(lme)%coef(se)
           lme1 = solhrules_product(lme)%lm_sh1(se)
           lme2 = solhrules_product(lme)%lm_sh2(se)
           LE2  = lof(lme2)

        ! PR(lmf,lmf1,lmf2):
        do sf=1,  solhrules_product(lmf)%n_summands
           cff  = solhrules_product(lmf)%coef(sf)
           lmf1 = solhrules_product(lmf)%lm_sh1(sf)
           lmf2 = solhrules_product(lmf)%lm_sh2(sf)
           LF2  = lof(lmf2)

           if( LA1 < LB2+LC2+LD2+LE2+LF2 ) cycle ! => D(lmf2)D(lme2)D(lmd2)D(lmc2)D(lmb2)*C(lma) == 0

           cf = cfb * cfc * cfd * cfe * cff

           F6(:,lma,lmb,lmc,lmd,lme,lmf) =  F6(:,lma,lmb,lmc,lmd,lme,lmf) &
                         + cf * F5(:,lmb1,lmc1,lmd1,lme1,lmf1) * D5C(:,lma,lmb2,lmc2,lmd2,lme2,lmf2)
           !                            \____\___\____\____\____ prb _______/    /    /   /    /
           !                                  \___\____\____\____ prc __________/    /   /    /
           !                                       \____\____\____ prd _____________/   /    /
           !                                             \____\_____ pre ______________/    /
           !                                                   \_____ prf _________________/
        enddo ! sf
        enddo ! se
        enddo ! sd
        enddo ! sc
        enddo ! sb
     enddo ! lmf
     enddo ! lme
     enddo ! lmd
     enddo ! lmc
     enddo ! lmb
     enddo !  ma
     enddo !  LA1

     ! ==================
     do lmf=  1,(LF+1)**2
     do lme=  1,(LE+1)**2
     do lmd=  1,(LD+1)**2
     do lmc=  1,(LC+1)**2
     do lmb=  1,(LB+1)**2
     do lma=  1,(LA+1)**2
        if( lma>=lmb .and. lmb>=lmc .and. lmc>=lmd .and. lmd>=lme .and. lme>=lmf ) cycle

        lm = srt6(lma,lmb,lmc,lmd,lme,lmf)
        F6(:,lma,lmb,lmc,lmd,lme,lmf) = F6(:,lm(1),lm(2),lm(3),lm(4),lm(5),lm(6))
     enddo
     enddo
     enddo
     enddo
     enddo
     enddo

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(F5)
#endif

   end subroutine ZD6Fv

   subroutine XD5Av(N,LA,LB,LC,LD,LE,W,D5A)
     ! prepare 4-derivatives of CLM
     ! then call ZD5Av to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: &
          D5A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2, * ) ! * = 1:1+SUM(L)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D4C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
#else
     real(RK), allocatable :: D4C(:,:,:,:,:,:)

     allocate(D4C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2))
#endif

     call SHR_D4C(N,LA,LB,LC,LD,LE,W,D4C)

     call ZD5Av(N,LA,LB,LC,LD,LE,D4C,D5A)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D4C)
#endif

   end subroutine XD5Av

   subroutine XD5Fv(N,LA,LB,LC,LD,LE,W,FL,F5)
     ! prepare 4-derivatives of CLM
     ! then call ZD5Fv to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD+LE)
     real(RK)   , intent(out) :: F5(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D4C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
#else
     real(RK), allocatable :: D4C(:,:,:,:,:,:)

     allocate(D4C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2))
#endif

     call SHR_D4C(N,LA,LB,LC,LD,LE,W,D4C)

     call ZD5Fv(N,LA,LB,LC,LD,LE,D4C,FL,F5)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D4C)
#endif
   end subroutine XD5Fv

   subroutine XD6Fv(N,LA,LB,LC,LD,LE,LF,W,FL,F6)
     ! prepare 5-derivatives of CLM
     ! then call ZD6Fv to do the real work
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE,LF
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD+LE+LF)
     real(RK)   , intent(out) :: F6(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: D5C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
#else
     real(RK), allocatable :: D5C(:,:,:,:,:,:,:)

     allocate(D5C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2))
#endif

     call SHR_D5C(N,LA,LB,LC,LD,LE,LF,W,D5C)

     call ZD6Fv(N,LA,LB,LC,LD,LE,LF,D5C,FL,F6)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(D5C)
#endif
   end subroutine XD6Fv

   subroutine SHR_D1C(N,LA,LB,W,D1C)
     ! prepare 1-derivatives of CLM
     use solid_harmonics_module, only: &
         solid_harmonics_calc
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: D1C(N,(LA+1)**2,(LB+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: CLM(N,(LA+1)**2) ! (N,(LMAX+1)**2)
#else
     real(RK), allocatable :: CLM(:,:)

     allocate(CLM(N,(LA+1)**2))
#endif

     CLM = solid_harmonics_calc(LA,W)

     call XD1C(N,LA,LB,CLM,D1C)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(CLM)
#endif

   end subroutine SHR_D1C

   subroutine SHR_D2C(N,LA,LB,LC,W,D2C)
     ! prepare 2-derivatives of CLM
     use solid_harmonics_module, only: &
         solid_harmonics_calc
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: D2C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: CLM(N,(LA+1)**2) ! (N,(LMAX+1)**2)
#else
     real(RK), allocatable :: CLM(:,:)

     allocate(CLM(N,(LA+1)**2))
#endif

     CLM = solid_harmonics_calc(LA,W)

     call XD2C(N,LA,LB,LC,CLM,D2C)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(CLM)
#endif


   end subroutine SHR_D2C

   subroutine SHR_D3C(N,LA,LB,LC,LD,W,D3C)
     ! prepare 3-derivatives of CLM
     use solid_harmonics_module, only: &
         solid_harmonics_calc
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: D3C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: CLM(N,(LA+1)**2) ! (N,(LMAX+1)**2)
#else
     real(RK), allocatable :: CLM(:,:)

     allocate(CLM(N,(LA+1)**2))
#endif

     CLM = solid_harmonics_calc(LA,W)

     call XD3C(N,LA,LB,LC,LD,CLM,D3C)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(CLM)
#endif

   end subroutine SHR_D3C

   subroutine SHR_D4C(N,LA,LB,LC,LD,LE,W,D4C)
     ! prepare 4-derivatives of CLM
     use solid_harmonics_module, only: &
         solid_harmonics_calc
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: D4C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: CLM(N,(LA+1)**2) ! (N,(LMAX+1)**2)
#else
     real(RK), allocatable :: CLM(:,:)

     allocate(CLM(N,(LA+1)**2))
#endif

     CLM = solid_harmonics_calc(LA,W)

     call XD4C(N,LA,LB,LC,LD,LE,CLM,D4C)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(CLM)
#endif

   end subroutine SHR_D4C

   subroutine SHR_D5C(N,LA,LB,LC,LD,LE,LF,W,D5C)
     ! prepare 5-derivatives of CLM
     use solid_harmonics_module, only: &
         solid_harmonics_calc
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE,LF
     real(RK)   , intent(in)  :: W(N,3) ! (N,3)
     real(RK)   , intent(out) :: D5C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
     ! *** end of interface ***

#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
     real(RK)    :: CLM(N,(LA+1)**2) ! (N,(LMAX+1)**2)
#else
     real(RK), allocatable :: CLM(:,:)

     allocate(CLM(N,(LA+1)**2))
#endif

     CLM = solid_harmonics_calc(LA,W)

     call XD5C(N,LA,LB,LC,LD,LE,LF,CLM,D5C)

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
     deallocate(CLM)
#endif

   end subroutine SHR_D5C

   subroutine XD1C(N,LA,LB,CLM,D1C)
     ! precompute diffs of C(lm)
     ! D1C(lma,lmb) := D(lmb) * C(lma)
     use constants, only: ZERO
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB
     real(RK)   , intent(in)  :: CLM(N,*) ! (N,(LMAX+1)**2)
     real(RK)   , intent(out) :: &
          D1C(N,(LA+1)**2,(LB+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb
     integer(IK) :: lm1
     integer(IK) :: sb
     real(RK)    :: cfb

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D1C(lma,lmb) := D(lmb) * C(lma)
     D1C = ZERO
     do lma=1,(LA+1)**2
        do lmb=1,(LB+1)**2
           if( lof(lmb) > lof(lma) ) cycle
           do sb=1, solhrules_differential(lmb,lma)%n_summands
              cfb = solhrules_differential(lmb,lma)%coef(sb)
              lm1 = solhrules_differential(lmb,lma)%lm_sh(sb)

              D1C(:,lma,lmb) = D1C(:,lma,lmb) &
                   + cfb * CLM(:,lm1)
        enddo !  sb
        enddo ! lmb
     enddo    ! lma
   end subroutine XD1C

   subroutine XD2C(N,LA,LB,LC,CLM,D2C)
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc) := D(lmc) D(lmb) * C(lma)
     use constants, only: ZERO
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC
     real(RK)   , intent(in)  :: CLM(N,*) ! (N,(LMAX+1)**2)
     real(RK)   , intent(out) :: &
          D2C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc
     integer(IK) :: lm1,lm2
     integer(IK) :: sb,sc
     real(RK)    :: cfb,cfc

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc) := D(lmc) D(lmb) * C(lma)
     D2C = ZERO
     do lma=1,(LA+1)**2
        do lmb=1,(LB+1)**2
           if( lof(lmb) > lof(lma) ) cycle
           do sb=1, solhrules_differential(lmb,lma)%n_summands
              cfb = solhrules_differential(lmb,lma)%coef(sb)
              lm1 = solhrules_differential(lmb,lma)%lm_sh(sb)
        do lmc=1,(LC+1)**2
           if( lof(lmc) > lof(lm1) ) cycle
           do sc=1, solhrules_differential(lmc,lm1)%n_summands
              cfc = solhrules_differential(lmc,lm1)%coef(sc)
              lm2 = solhrules_differential(lmc,lm1)%lm_sh(sc)

              D2C(:,lma,lmb,lmc) = D2C(:,lma,lmb,lmc) &
                   + cfb * cfc * CLM(:,lm2)
        enddo !  sc
        enddo ! lmc
        enddo !  sb
        enddo ! lmb
     enddo    ! lma
   end subroutine XD2C

   subroutine XD3C(N,LA,LB,LC,LD,CLM,D3C)
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd) := D(lmd) D(lmc) D(lmb) * C(lma)
     use constants, only: ZERO
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD
     real(RK)   , intent(in)  :: CLM(N,*) ! (N,(LMAX+1)**2)
     real(RK)   , intent(out) :: &
          D3C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd
     integer(IK) :: lm1,lm2,lm3
     integer(IK) :: sb,sc,sd
     real(RK)    :: cfb,cfc,cfd

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd) := D(lmd) D(lmc) D(lmb) * C(lma)
     D3C = ZERO
     do lma=1,(LA+1)**2
        do lmb=1,(LB+1)**2
           if( lof(lmb) > lof(lma) ) cycle
           do sb=1, solhrules_differential(lmb,lma)%n_summands
              cfb = solhrules_differential(lmb,lma)%coef(sb)
              lm1 = solhrules_differential(lmb,lma)%lm_sh(sb)
        do lmc=1,(LC+1)**2
           if( lof(lmc) > lof(lm1) ) cycle
           do sc=1, solhrules_differential(lmc,lm1)%n_summands
              cfc = solhrules_differential(lmc,lm1)%coef(sc)
              lm2 = solhrules_differential(lmc,lm1)%lm_sh(sc)
        do lmd=1,(LD+1)**2
           if( lof(lmd) > lof(lm2) ) cycle
           do sd=1, solhrules_differential(lmd,lm2)%n_summands
              cfd = solhrules_differential(lmd,lm2)%coef(sd)
              lm3 = solhrules_differential(lmd,lm2)%lm_sh(sd)

              D3C(:,lma,lmb,lmc,lmd) = D3C(:,lma,lmb,lmc,lmd) &
                   + cfb * cfc * cfd * CLM(:,lm3)
        enddo !  sd
        enddo ! lmd
        enddo !  sc
        enddo ! lmc
        enddo !  sb
        enddo ! lmb
     enddo    ! lma
   end subroutine XD3C

   subroutine XD4C(N,LA,LB,LC,LD,LE,CLM,D4C)
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd,lme) := D(lme) D(lmd) D(lmc) D(lmb) * C(lma)
     use constants, only: ZERO
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE
     real(RK)   , intent(in)  :: CLM(N,*) ! (N,(LMAX+1)**2)
     real(RK)   , intent(out) :: &
          D4C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd,lme
     integer(IK) :: lm1,lm2,lm3,lm4
     integer(IK) :: sb,sc,sd,se
     real(RK)    :: cfb,cfc,cfd,cfe

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D3C(lma,lmb,lmc,lmd,lme) := D(lme) D(lmd) D(lmc) D(lmb) * C(lma)
     D4C = ZERO
     do lma=1,(LA+1)**2
        do lmb=1,(LB+1)**2
           if( lof(lmb) > lof(lma) ) cycle
           do sb=1, solhrules_differential(lmb,lma)%n_summands
              cfb = solhrules_differential(lmb,lma)%coef(sb)
              lm1 = solhrules_differential(lmb,lma)%lm_sh(sb)
        do lmc=1,(LC+1)**2
           if( lof(lmc) > lof(lm1) ) cycle
           do sc=1, solhrules_differential(lmc,lm1)%n_summands
              cfc = solhrules_differential(lmc,lm1)%coef(sc)
              lm2 = solhrules_differential(lmc,lm1)%lm_sh(sc)
        do lmd=1,(LD+1)**2
           if( lof(lmd) > lof(lm2) ) cycle
           do sd=1, solhrules_differential(lmd,lm2)%n_summands
              cfd = solhrules_differential(lmd,lm2)%coef(sd)
              lm3 = solhrules_differential(lmd,lm2)%lm_sh(sd)
        do lme=1,(LE+1)**2
           if( lof(lme) > lof(lm3) ) cycle
           do se=1, solhrules_differential(lme,lm3)%n_summands
              cfe = solhrules_differential(lme,lm3)%coef(se)
              lm4 = solhrules_differential(lme,lm3)%lm_sh(se)

              D4C(:,lma,lmb,lmc,lmd,lme) = D4C(:,lma,lmb,lmc,lmd,lme) &
                   + cfb * cfc * cfd * cfe * CLM(:,lm4)
        enddo !  se
        enddo ! lme
        enddo !  sd
        enddo ! lmd
        enddo !  sc
        enddo ! lmc
        enddo !  sb
        enddo ! lmb
     enddo    ! lma
   end subroutine XD4C

   subroutine XD5C(N,LA,LB,LC,LD,LE,LF,CLM,D5C)
     ! precompute diffs of C(lm)
     ! D5C(lma,lmb,lmc,lmd,lme,lmf) := D(lmf) D(lme) D(lmd) D(lmc) D(lmb) * C(lma)
     use constants, only: ZERO
     implicit none
     integer(IK), intent(in)  :: N
     integer(IK), intent(in)  :: LA,LB,LC,LD,LE,LF
     real(RK)   , intent(in)  :: CLM(N,*) ! (N,(LMAX+1)**2)
     real(RK)   , intent(out) :: &
          D5C(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd,lme,lmf
     integer(IK) :: lm1,lm2,lm3,lm4,lm5
     integer(IK) :: sb,sc,sd,se,sf
     real(RK)    :: cfb,cfc,cfd,cfe,cff

     ! ================================================================
     ! precompute diffs of C(lm)
     ! D5C(lma,lmb,lmc,lmd,lme,lmf) := D(lf) D(lme) D(lmd) D(lmc) D(lmb) * C(lma)
     D5C = ZERO
     do lma=1,(LA+1)**2
        do lmb=1,(LB+1)**2
           if( lof(lmb) > lof(lma) ) cycle
           do sb=1, solhrules_differential(lmb,lma)%n_summands
              cfb = solhrules_differential(lmb,lma)%coef(sb)
              lm1 = solhrules_differential(lmb,lma)%lm_sh(sb)
        do lmc=1,(LC+1)**2
           if( lof(lmc) > lof(lm1) ) cycle
           do sc=1, solhrules_differential(lmc,lm1)%n_summands
              cfc = solhrules_differential(lmc,lm1)%coef(sc)
              lm2 = solhrules_differential(lmc,lm1)%lm_sh(sc)
        do lmd=1,(LD+1)**2
           if( lof(lmd) > lof(lm2) ) cycle
           do sd=1, solhrules_differential(lmd,lm2)%n_summands
              cfd = solhrules_differential(lmd,lm2)%coef(sd)
              lm3 = solhrules_differential(lmd,lm2)%lm_sh(sd)
        do lme=1,(LE+1)**2
           if( lof(lme) > lof(lm3) ) cycle
           do se=1, solhrules_differential(lme,lm3)%n_summands
              cfe = solhrules_differential(lme,lm3)%coef(se)
              lm4 = solhrules_differential(lme,lm3)%lm_sh(se)
        do lmf=1,(LF+1)**2
           if( lof(lmf) > lof(lm4) ) cycle
           do sf=1, solhrules_differential(lmf,lm4)%n_summands
              cff = solhrules_differential(lmf,lm4)%coef(sf)
              lm5 = solhrules_differential(lmf,lm4)%lm_sh(sf)

              D5C(:,lma,lmb,lmc,lmd,lme,lmf) = D5C(:,lma,lmb,lmc,lmd,lme,lmf) &
                   + cfb * cfc * cfd * cfe * cff * CLM(:,lm5)
        enddo !  sf
        enddo ! lmf
        enddo !  se
        enddo ! lme
        enddo !  sd
        enddo ! lmd
        enddo !  sc
        enddo ! lmc
        enddo !  sb
        enddo ! lmb
     enddo    ! lma
   end subroutine XD5C

   subroutine SHR_D2Av(N,LA,LB,W,D2A)
     ! Computes angular part of double derivative
     !      F2(lma,lmb) = D(lma) D(lmb) * F
     !
     ! (vector version)
     implicit none
     integer(IK), intent(in) :: N
     integer(IK), intent(in) :: LA,LB
     real(RK), intent(in)    :: W(N,3) ! (N,3)
     real(RK), intent(out)   :: &
          D2A(N,(LA+1)**2,(LB+1)**2,*) ! * == 1+LA+LB)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,LA1,LB1
     integer(IK) :: llo,lhi
     real(RK), allocatable :: T2A(:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,1+LA+LB)

     if( LA < LB )then
        allocate(T2A(N,(LB+1)**2,(LA+1)**2,1+LA+LB))

        call XD2Av(N,LB,LA,W,T2A)

        D2A(:,:,:,:1+LA+LB) = 0.0_rk

        do lmb=1,(LB+1)**2
           LB1 = lof(lmb)
        do lma=1,(LA+1)**2
           LA1 = lof(lma)
           llo = 1+ max(LA1,LB1)
           lhi = 1+     LA1+LB1
           D2A(:,lma,lmb,llo:lhi) = T2A(:,lmb,lma,llo:lhi)
        enddo
        enddo

        deallocate(T2A)
     else
        call XD2Av(N,LA,LB,W,D2A)
     endif
   end subroutine SHR_D2Av

   subroutine SHR_D2Fv(N,LA,LB,W,FL,F2)
     ! Computes double derivative
     !      F2(lma,lmb) = D(lma) D(lmb) * F
     !
     implicit none
     integer(IK), intent(in)  :: N,LA,LB
     real(RK)   , intent(in)  :: W(N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB)
     real(RK)   , intent(out) :: F2(N,(LA+1)**2,(LB+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb
     real(RK), allocatable :: T2(:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2)

     if( LA < LB )then

        allocate(T2(N,(LB+1)**2,(LA+1)**2))

        call XD2Fv(N,LB,LA,W,FL,T2)

        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           F2(:,lma,lmb) = T2(:,lmb,lma)
        enddo
        enddo

        deallocate(T2)
     else
        call XD2Fv(N,LA,LB,W,FL,F2)
     endif
   end subroutine SHR_D2Fv

   subroutine SHR_D3Av(N,LA,LB,LC,W,D3A)
     ! Computes angular part of triple derivative
     !      F3(lma,lmb,lmc) = D(lma) D(lmb) D(lmc) * F
     !
     ! (vector version)
     implicit none
     integer(IK), intent(in) :: N
     integer(IK), intent(in) :: LA,LB,LC
     real(RK), intent(in)    :: W(N,3) ! (N,3)
     real(RK), intent(out)   :: &
          D3A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,*) ! * == 1+LA+LB+LC)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc
     integer(IK) :: L(3),lm(3)
     integer(IK) :: llo,lhi
     real(RK), allocatable :: T3A(:,:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2,1+LA+LB+LC)

     if( LC == 0 )then
        call SHR_D2Av(N,LA,LB,   W,D3A)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D2Av(N,LA,   LC,W,D3A)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D2Av(N,   LB,LC,W,D3A)
        RETURN
     endif

     L = srt3(LA,LB,LC)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) )then

        allocate(T3A(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2,1+SUM(L)))

        call XD3Av(N,L(1),L(2),L(3),W,T3A)

        D3A(:,:,:,:,:1+SUM(L)) = 0.0_rk

        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt3(lma,lmb,lmc)
           L  = lof(lm)
           llo = 1+ maxval(L)
           lhi = 1+    sum(L)
           D3A(:,lma,lmb,lmc,llo:lhi) = T3A(:,lm(1),lm(2),lm(3),llo:lhi)
        enddo
        enddo
        enddo

        deallocate(T3A)
     else
        call XD3Av(N,L(1),L(2),L(3),W,D3A)
     endif
   end subroutine SHR_D3Av

   subroutine SHR_D3Fv(N,LA,LB,LC,W,FL,F3)
     ! Computes 3-derivative
     !      F3(lma,lmb,lmc) = D(lma) D(lmb) D(lmc) * F
     !
     implicit none
     integer(IK), intent(in)  :: N,LA,LB,LC
     real(RK)   , intent(in)  :: W(N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC)
     real(RK)   , intent(out) :: F3(N,(LA+1)**2,(LB+1)**2,(LC+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc
     integer(IK) :: L(3),lm(3)
     real(RK), allocatable :: T3(:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2)

     if( LC == 0 )then
        call SHR_D2Fv(N,LA,LB,   W,FL,F3)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D2Fv(N,LA,   LC,W,FL,F3)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D2Fv(N,   LB,LC,W,FL,F3)
        RETURN
     endif

     L = srt3(LA,LB,LC)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) )then

        allocate(T3(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2))

        call XD3Fv(N,L(1),L(2),L(3),W,FL,T3)

        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt3(lma,lmb,lmc)
           F3(:,lma,lmb,lmc) = T3(:,lm(1),lm(2),lm(3))
        enddo
        enddo
        enddo

        deallocate(T3)
     else
        call XD3Fv(N,L(1),L(2),L(3),W,FL,F3)
     endif
   end subroutine SHR_D3Fv

   subroutine SHR_D4Av(N,LA,LB,LC,LD,W,D4A)
     ! Computes angular part of 4-derivative
     !      F4(lma,lmb,lmc,lmd) = D(lma)D(lmb)D(lmc)D(lmd) * F
     !
     ! (vector version)
     implicit none
     integer(IK), intent(in) :: N
     integer(IK), intent(in) :: LA,LB,LC,LD
     real(RK), intent(in)    :: W(N,3) ! (N,3)
     real(RK), intent(out)   :: &
          D4A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,*) ! * == 1+LA+LB+LC+LD)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd
     integer(IK) :: L(4),lm(4)
     integer(IK) :: llo,lhi
     real(RK), allocatable :: T4A(:,:,:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,1+LA+LB+LC+LD)

     if( LD == 0 )then
        call SHR_D3Av(N,LA,LB,LC,   W,D4A)
        RETURN
     endif
     if( LC == 0 )then
        call SHR_D3Av(N,LA,LB,   LD,W,D4A)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D3Av(N,LA,   LC,LD,W,D4A)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D3Av(N,   LB,LC,LD,W,D4A)
        RETURN
     endif

     L = srt4(LA,LB,LC,LD)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) .or. LD/=L(4) )then

        allocate(T4A(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2,(L(4)+1)**2,1+SUM(L)))

        call XD4Av(N,L(1),L(2),L(3),L(4),W,T4A)

        D4A(:,:,:,:,:,:1+SUM(L)) = 0.0_rk

        do lmd=1,(LD+1)**2
        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt4(lma,lmb,lmc,lmd)
           L  = lof(lm)
           llo = 1+ maxval(L)
           lhi = 1+    sum(L)
           D4A(:,lma,lmb,lmc,lmd,llo:lhi) = T4A(:,lm(1),lm(2),lm(3),lm(4),llo:lhi)
        enddo
        enddo
        enddo
        enddo

        deallocate(T4A)
     else
        call XD4Av(N,L(1),L(2),L(3),L(4),W,D4A)
     endif
   end subroutine SHR_D4Av

   subroutine SHR_D4Fv(N,LA,LB,LC,LD,W,FL,F4)
     ! Computes 4-derivative
     !      F4(lma,lmb,lmc,lmd) = D(lma)D(lmb)D(lmc)D(lmd) * F
     !
     implicit none
     integer(IK), intent(in)  :: N,LA,LB,LC,LD
     real(RK)   , intent(in)  :: W(N,3)
     real(RK)   , intent(in)  :: FL(N,1+LA+LB+LC+LD)
     real(RK)   , intent(out) :: F4(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd
     integer(IK) :: L(4),lm(4)
     real(RK), allocatable :: T4(:,:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)

     if( LD == 0 )then
        call SHR_D3Fv(N,LA,LB,LC,   W,FL,F4)
        RETURN
     endif
     if( LC == 0 )then
        call SHR_D3Fv(N,LA,LB,   LD,W,FL,F4)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D3Fv(N,LA,   LC,LD,W,FL,F4)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D3Fv(N,   LB,LC,LD,W,FL,F4)
        RETURN
     endif

     L = srt4(LA,LB,LC,LD)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) .or. LD/=L(4) )then

        allocate(T4(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2,(L(4)+1)**2))

        call XD4Fv(N,L(1),L(2),L(3),L(4),W,FL,T4)

        do lmd=1,(LD+1)**2
        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt4(lma,lmb,lmc,lmd)
           F4(:,lma,lmb,lmc,lmd) = T4(:,lm(1),lm(2),lm(3),lm(4))
        enddo
        enddo
        enddo
        enddo

        deallocate(T4)
     else
        call XD4Fv(N,L(1),L(2),L(3),L(4),W,FL,F4)
     endif
   end subroutine SHR_D4Fv

   subroutine SHR_D5Av(N,LA,LB,LC,LD,LE,W,D5A)
     ! Computes angular part of 5-derivative
     !      F5(lma,..,lme) = D(lma)..D(lme) * F
     !
     ! (vector version)
     implicit none
     integer(IK), intent(in) :: N
     integer(IK), intent(in) :: LA,LB,LC,LD,LE
     real(RK), intent(in)    :: W(N,3) ! (N,3)
     real(RK), intent(out)   :: &
          D5A(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,*) ! * == 1+SUM(L)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd,lme
     integer(IK) :: L(5),lm(5)
     integer(IK) :: llo,lhi
     real(RK), allocatable :: T5A(:,:,:,:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+SUM(L))

     if( LE == 0 )then
        call SHR_D4Av(N,LA,LB,LC,LD,   W,D5A)
        RETURN
     endif
     if( LD == 0 )then
        call SHR_D4Av(N,LA,LB,LC,   LE,W,D5A)
        RETURN
     endif
     if( LC == 0 )then
        call SHR_D4Av(N,LA,LB,   LD,LE,W,D5A)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D4Av(N,LA,   LC,LD,LE,W,D5A)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D4Av(N,   LB,LC,LD,LE,W,D5A)
        RETURN
     endif

     L = srt5(LA,LB,LC,LD,LE)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) .or. LD/=L(4) .or. LE/=L(5) )then

        allocate(T5A(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2,(L(4)+1)**2,(L(5)+1)**2,1+SUM(L)))

        call XD5Av(N,L(1),L(2),L(3),L(4),L(5),W,T5A)

        D5A(:,:,:,:,:,:,:1+SUM(L)) = 0.0_rk

        do lme=1,(LE+1)**2
        do lmd=1,(LD+1)**2
        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt5(lma,lmb,lmc,lmd,lme)
           L  = lof(lm)
           llo = 1+ maxval(L)
           lhi = 1+    sum(L)
           D5A(:,lma,lmb,lmc,lmd,lme,llo:lhi) = T5A(:,lm(1),lm(2),lm(3),lm(4),lm(5),llo:lhi)
        enddo
        enddo
        enddo
        enddo
        enddo

        deallocate(T5A)
     else
        call XD5Av(N,L(1),L(2),L(3),L(4),L(5),W,D5A)
     endif
   end subroutine SHR_D5Av

   subroutine SHR_D5Fv(N,LA,LB,LC,LD,LE,W,FL,F5)
     ! Computes the 5-derivative
     !      F5(lma,..,lme) = D(lma)..D(lme) * F
     !
     ! (vector version)
     implicit none
     integer(IK), intent(in) :: N
     integer(IK), intent(in) :: LA,LB,LC,LD,LE
     real(RK), intent(in)    :: W(N,3) ! (N,3)
     real(RK), intent(in)    :: FL(N,1+LA+LB+LC+LD+LE) !?
     real(RK), intent(out)   :: F5(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd,lme
     integer(IK) :: L(5),lm(5)
     real(RK), allocatable :: T5(:,:,:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)

     if( LE == 0 )then
        call SHR_D4Fv(N,LA,LB,LC,LD,   W,FL,F5)
        RETURN
     endif
     if( LD == 0 )then
        call SHR_D4Fv(N,LA,LB,LC,   LE,W,FL,F5)
        RETURN
     endif
     if( LC == 0 )then
        call SHR_D4Fv(N,LA,LB,   LD,LE,W,FL,F5)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D4Fv(N,LA,   LC,LD,LE,W,FL,F5)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D4Fv(N,   LB,LC,LD,LE,W,FL,F5)
        RETURN
     endif

     L = srt5(LA,LB,LC,LD,LE)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) .or. LD/=L(4) .or. LE/=L(5) )then

        allocate(T5(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2,(L(4)+1)**2,(L(5)+1)**2))

        call XD5Fv(N,L(1),L(2),L(3),L(4),L(5),W,FL,T5)

        do lme=1,(LE+1)**2
        do lmd=1,(LD+1)**2
        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt5(lma,lmb,lmc,lmd,lme)
           F5(:,lma,lmb,lmc,lmd,lme) = T5(:,lm(1),lm(2),lm(3),lm(4),lm(5))
        enddo
        enddo
        enddo
        enddo
        enddo

        deallocate(T5)
     else
        call XD5Fv(N,L(1),L(2),L(3),L(4),L(5),W,FL,F5)
     endif
   end subroutine SHR_D5Fv

   subroutine SHR_D6Fv(N,LA,LB,LC,LD,LE,LF,W,FL,F6)
     ! Computes the 6-derivative
     !      F5(lma,..,lmf) = D(lma)..D(lmf) * F
     !
     ! (vector version)
     implicit none
     integer(IK), intent(in) :: N
     integer(IK), intent(in) :: LA,LB,LC,LD,LE,LF
     real(RK), intent(in)    :: W(N,3) ! (N,3)
     real(RK), intent(in)    :: FL(N,1+LA+LB+LC+LD+LE+LF)
     real(RK), intent(out)   :: F6(N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lmc,lmd,lme,lmf
     integer(IK) :: L(6),lm(6)
     real(RK), allocatable :: T6(:,:,:,:,:,:,:)
     ! (N,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LE+1)**2)

     if( LF == 0 )then
        call SHR_D5Fv(N,LA,LB,LC,LD,LE,   W,FL,F6)
        RETURN
     endif
     if( LE == 0 )then
        call SHR_D5Fv(N,LA,LB,LC,LD,  LF, W,FL,F6)
        RETURN
     endif
     if( LD == 0 )then
        call SHR_D5Fv(N,LA,LB,LC,   LE,LF,W,FL,F6)
        RETURN
     endif
     if( LC == 0 )then
        call SHR_D5Fv(N,LA,LB,   LD,LE,LF,W,FL,F6)
        RETURN
     endif
     if( LB == 0 )then
        call SHR_D5Fv(N,LA,   LC,LD,LE,LF,W,FL,F6)
        RETURN
     endif
     if( LA == 0 )then
        call SHR_D5Fv(N,   LB,LC,LD,LE,LF,W,FL,F6)
        RETURN
     endif

     L = srt6(LA,LB,LC,LD,LE,LF)

     if( LA/=L(1) .or. LB/=L(2) .or. LC/=L(3) .or. LD/=L(4) .or. LE/=L(5) .or. LF/=L(6) )then

        allocate(T6(N,(L(1)+1)**2,(L(2)+1)**2,(L(3)+1)**2,(L(4)+1)**2,(L(5)+1)**2,(L(6)+1)**2))

        call XD6Fv(N,L(1),L(2),L(3),L(4),L(5),L(6),W,FL,T6)

        do lmf=1,(LF+1)**2
        do lme=1,(LE+1)**2
        do lmd=1,(LD+1)**2
        do lmc=1,(LC+1)**2
        do lmb=1,(LB+1)**2
        do lma=1,(LA+1)**2
           lm = srt6(lma,lmb,lmc,lmd,lme,lmf)
           F6(:,lma,lmb,lmc,lmd,lme,lmf) = T6(:,lm(1),lm(2),lm(3),lm(4),lm(5),lm(6))
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        deallocate(T6)
     else
        call XD6Fv(N,L(1),L(2),L(3),L(4),L(5),L(6),W,FL,F6)
     endif
   end subroutine SHR_D6Fv

   subroutine SHR_PR2v(N,LA,LB,F,S,FS,FACT)
     ! couples double derivatives
     !    F(lm1,ln1)
     ! and
     !    S(lm2,ln2)
     ! according to the (double) SH product rule
     !
     !    FS(ma,mb) += PR(lma,lm1,lm2) PR(lmb,ln1,ln2)
     !                  *  F(lm1,ln1) * S(lm2,ln2)
     !
     ! WARNING: ADDS RESULT TO FS, DONT FORGET TO CLEAR FS!
     !
     implicit none
     integer(IK), intent(in)    :: N
     integer(IK), intent(in)    :: LA,LB
     real(RK)   , intent(in)    :: F(:,:,:)  ! (N,(LA+1)**2,(LB+1)**2)
     real(RK)   , intent(in)    :: S(:,:,:)  ! (N,(LA+1)**2,(LB+1)**2)
     real(RK)   , intent(inout) :: FS(:,:,:) ! (N,2*LA+1,2*LB+1)
     real(RK)   , intent(in)    :: FACT(:)   ! (N)
     optional :: FACT
     ! *** end of interface ***

     integer(IK) :: ma,mb,lma,lmb,lm1,lm2,ln1,ln2,ta,tb
     real(RK)    :: cfa,cfb,cf
     real(RK)    :: dFS(N)

     ! not zeroing FS!, we increment it instead!
     ASSERT(N==size(F,1))
     ASSERT(N==size(S,1))
     ASSERT(N==size(FS,1))

     if( .not.present(FACT) )then
       do mb=1,2*LB+1
          lmb = LB**2 + mb
       do ma=1,2*LA+1
          lma = LA**2 + ma

       dFS = 0.0_rk

       do ta=1, solhrules_product(lma)%n_summands
          cfa = solhrules_product(lma)%coef(ta)
          lm1 = solhrules_product(lma)%lm_sh1(ta)
          lm2 = solhrules_product(lma)%lm_sh2(ta)
       do tb=1, solhrules_product(lmb)%n_summands
          cfb = solhrules_product(lmb)%coef(tb)
          ln1 = solhrules_product(lmb)%lm_sh1(tb)
          ln2 = solhrules_product(lmb)%lm_sh2(tb)

          cf = cfa * cfb

          dFS = dFS + cf * F(:,lm1,ln1) * S(:,lm2,ln2)
       enddo ! tb
       enddo ! ta
          FS(:,ma,mb) = FS(:,ma,mb) + dFS
       enddo ! ma
       enddo ! mb
     else
       ! same code as above but see NOTE below:
       do mb=1,2*LB+1
          lmb = LB**2 + mb
       do ma=1,2*LA+1
          lma = LA**2 + ma

       dFS = 0.0_rk

       do ta=1, solhrules_product(lma)%n_summands
          cfa = solhrules_product(lma)%coef(ta)
          lm1 = solhrules_product(lma)%lm_sh1(ta)
          lm2 = solhrules_product(lma)%lm_sh2(ta)
       do tb=1, solhrules_product(lmb)%n_summands
          cfb = solhrules_product(lmb)%coef(tb)
          ln1 = solhrules_product(lmb)%lm_sh1(tb)
          ln2 = solhrules_product(lmb)%lm_sh2(tb)

          cf = cfa * cfb

          dFS = dFS + cf * F(:,lm1,ln1) * S(:,lm2,ln2)
       enddo ! tb
       enddo ! ta
          FS(:,ma,mb) = FS(:,ma,mb) + dFS * FACT ! NOTE the FACT factor!
       enddo ! ma
       enddo ! mb
     endif
   end subroutine SHR_PR2v

   subroutine SHR_doPR2Ls(LMAS,LMAE,LMBS,LMBE,Y,S,YLS)
     ! Double Product Rule preserving orders Y(L)
     ! for the range lma=LMAS:LMAE, lmb=LMBS:LMBE
     !
     ! couples "angular" factors
     !    Y(lm1,ln1,L)
     ! and
     !    S(lm2,ln2)
     ! performing double (L)-preserving product rule:
     !
     !    YLS(lma,lmb,L) = PR(lma,lm1,lm2) PR(lmb,ln1,ln2)
     !                   *  Y(lm1,ln1,L) * S(lm2,ln2)
     !
     implicit none
     integer(IK), intent(in) :: LMAS,LMAE
     integer(IK), intent(in) :: LMBS,LMBE
     real(RK), intent(in)    :: Y(:,:,:) ! ((LA+1)**2,(LB+1)**2,1+LA+LB)
     real(RK), intent(in)    :: S(:,:)   ! ((LA+1)**2,(LB+1)**2        )
     real(RK), intent(inout) :: YLS(LMAS:,LMBS:,:)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lm1,lm2,ln1,ln2,ta,tb,L
     real(RK)    :: cfa,cfb
     integer(IK) :: LA1,LB1

     lma=size(Y,1)
     ASSERT(LMAE-LMAS+1<=lma)
     lmb=size(Y,2)
     ASSERT(LMBE-LMBS+1<=lmb)
     lma = lof(LMAE)
     lmb = lof(LMBE)
     ASSERT(1+lma+lmb<=size(Y,3))

     ! FIXME: dont forget to YLS = 0.0_rk
     do lmb=LMBS,LMBE
        do lma=LMAS,LMAE
              do ta=1,solhrules_product(lma)%n_summands
                 cfa  = solhrules_product(lma)%coef(ta)
                 lm1 = solhrules_product(lma)%lm_sh1(ta)
                 lm2 = solhrules_product(lma)%lm_sh2(ta)
                 LA1 = lof(lm1)
              do tb=1,solhrules_product(lmb)%n_summands
                 cfb  = solhrules_product(lmb)%coef(tb)
                 ln1 = solhrules_product(lmb)%lm_sh1(tb)
                 ln2 = solhrules_product(lmb)%lm_sh2(tb)
                 LB1 = lof(ln1)

                 do L=1+max(LA1,LB1),1+LA1+LB1
                    YLS(lma,lmb,L) = YLS(lma,lmb,L) &
                         + cfa * cfb * Y(lm1,ln1,L) * S(lm2,ln2)
                 enddo
              enddo ! tb
              enddo ! ta
        enddo
     enddo
   end subroutine SHR_doPR2Ls

   function SHR_PR2Ls(LA,LB,Y,S) result(YLS)
     ! Double Product Rule preserving orders Y(L):
     !
     ! couples "angular" factors
     !    Y(lm1,ln1,L)
     ! and
     !    S(lm2,ln2)
     ! performing double (L)-preserving product rule:
     !
     !    YLS(ma,mb,L) = PR(lma,lm1,lm2) PR(lmb,ln1,ln2)
     !                   *  Y(lm1,ln1,L) * S(lm2,ln2)
     !
     implicit none
     integer(IK), intent(in) :: LA,LB
     real(RK), intent(in)    :: Y(:,:,:) ! ((LA+1)**2,(LB+1)**2,1+LA+LB)
     real(RK), intent(in)    :: S(:,:)   ! ((LA+1)**2,(LB+1)**2        )
     real(RK)                :: YLS(2*LA+1,2*LB+1,1+LA+LB)
     ! *** end of interface ***

     YLS = 0.0_rk
     call SHR_doPR2Ls(LA**2+1,(LA+1)**2,LB**2+1,(LB+1)**2,Y,S,YLS)
   end function SHR_PR2Ls

   function SHR_YLSv(LA,LB,NG,Y,S,FAC) result(YLS)
     ! couples "angular" factors
     !    Y(lm1,ln1,L)
     ! and
     !    S(lm2,ln2)
     ! to get the 3rd-center exponent independent factor
     !
     !    YLS(ma,mb,L) = PR(lma,lm1,lm2) PR(lmb,ln1,ln2)
     !                   *  Y(lm1,ln1,L) * S(lm2,ln2)
     !
     implicit none
     integer(IK), intent(in) :: LA,LB
     integer(IK), intent(in) :: NG         ! number of extra ``gradients''
     real(RK), intent(in)    :: Y(:,:,:,:) ! (NAB,(LA+1)**2,(LB+1)**2,1+LA+LB+x)
     real(RK), intent(in)    :: S(:,:,:)   ! (NAB,(LA+1)**2,(LB+1)**2        )
     real(RK), intent(in)    :: FAC(:)     ! (NAB)
     optional                :: FAC
     real(RK)                :: YLS(size(Y,1),2*LA+1,2*LB+1,size(Y,4))
     ! *** end of interface ***

     integer(IK) :: ma,mb,lma,lmb,lm1,lm2,ln1,ln2,ta,tb,L
     real(RK)    :: cfa,cfb
     integer(IK) :: LA1,LB1
     integer(IK) :: LC

     LC = 0
     if( NG > 0 ) LC = 1

     YLS = 0.0_rk
     do mb=1,2*LB+1
        lmb = LB**2 + mb
        do ma=1,2*LA+1
           lma = LA**2 + ma
              do ta=1,solhrules_product(lma)%n_summands
                 cfa  = solhrules_product(lma)%coef(ta)
                 lm1 = solhrules_product(lma)%lm_sh1(ta)
                 lm2 = solhrules_product(lma)%lm_sh2(ta)
                 LA1 = lof(lm1)
                 do tb=1,solhrules_product(lmb)%n_summands
                    cfb  = solhrules_product(lmb)%coef(tb)
                    ln1 = solhrules_product(lmb)%lm_sh1(tb)
                    ln2 = solhrules_product(lmb)%lm_sh2(tb)
                    LB1 = lof(ln1)
                    if( .not.present(FAC) )then
                    do L=1+max(LA1,LB1,LC),1+LA1+LB1+(LC*NG)
                       YLS(:,ma,mb,L) = YLS(:,ma,mb,L) &
                            + cfa * cfb * Y(:,lm1,ln1,L) * S(:,lm2,ln2)
                    enddo
                    else
                    do L=1+max(LA1,LB1,LC),1+LA1+LB1+(LC*NG)
                       YLS(:,ma,mb,L) = YLS(:,ma,mb,L) &
                            + cfa * cfb * Y(:,lm1,ln1,L) * S(:,lm2,ln2) &
                            * FAC
                    enddo
                    endif
                 enddo
              enddo
        enddo
     enddo
   end function SHR_YLSv

   subroutine SHR_doPR3Ls(LMDS,LMDE,LMAS,LMAE,LMBS,LMBE,Y,S,YLS)
     ! Double Product Rule preserving orders Y(L)
     ! for the range lma=LMAS:LMAE, lmb=LMBS:LMBE
     !
     ! couples "angular" factors
     !    Y(lm1,ln1,L)
     ! and
     !    S(lm2,ln2)
     ! performing double (L)-preserving product rule:
     !
     !    YLS(lma,lmb,L) = PR(lma,lm1,lm2) PR(lmb,ln1,ln2)
     !                   *  Y(lm1,ln1,L) * S(lm2,ln2)
     !
     implicit none
     integer(IK), intent(in) :: LMDS,LMDE
     integer(IK), intent(in) :: LMAS,LMAE
     integer(IK), intent(in) :: LMBS,LMBE
     real(RK), intent(in)    :: Y(:,:,:,:) ! ((LD+1)**2,(LA+1)**2,(LB+1)**2,1+LD+LA+LB)
     real(RK), intent(in)    :: S(:,:,:)   ! ((LD+1)**2,(LA+1)**2,(LB+1)**2        )
     real(RK), intent(inout) :: YLS(LMDS:,LMAS:,LMBS:,:)
     ! *** end of interface ***

     integer(IK) :: lma,lmb,lm1,lm2,ln1,ln2,ta,tb,L
     integer(IK) :: lmd,lk1,lk2,td
     real(RK)    :: cfd,cfa,cfb
     integer(IK) :: LD1,LA1,LB1

     lmd=size(Y,1)
     ASSERT(LMDE-LMDS+1<=lmd)
     lma=size(Y,2)
     ASSERT(LMAE-LMAS+1<=lma)
     lmb=size(Y,3)
     ASSERT(LMBE-LMBS+1<=lmb)
     lmd = solhrules_l_and_m_of_lm(1,LMDE)
     lma = solhrules_l_and_m_of_lm(1,LMAE)
     lmb = solhrules_l_and_m_of_lm(1,LMBE)
     ASSERT(1+lma+lmb+lmd<=size(Y,4))

     ! FIXME: dont forget to YLS = 0.0_rk
     do lmb=LMBS,LMBE
        do lma=LMAS,LMAE
          do lmd=LMDS,LMDE
              do td=1,solhrules_product(lmd)%n_summands
                 cfd  = solhrules_product(lmd)%coef(td)
                 lk1 = solhrules_product(lmd)%lm_sh1(td)
                 lk2 = solhrules_product(lmd)%lm_sh2(td)
                 LD1 = solhrules_l_and_m_of_lm(1,lk1)
              do ta=1,solhrules_product(lma)%n_summands
                 cfa  = solhrules_product(lma)%coef(ta)
                 lm1 = solhrules_product(lma)%lm_sh1(ta)
                 lm2 = solhrules_product(lma)%lm_sh2(ta)
                 LA1 = solhrules_l_and_m_of_lm(1,lm1)
              do tb=1,solhrules_product(lmb)%n_summands
                 cfb  = solhrules_product(lmb)%coef(tb)
                 ln1 = solhrules_product(lmb)%lm_sh1(tb)
                 ln2 = solhrules_product(lmb)%lm_sh2(tb)
                 LB1 = solhrules_l_and_m_of_lm(1,ln1)

                 do L=1+max(LA1,LB1,LD1),1+LA1+LB1+LD1
                    YLS(lmd,lma,lmb,L) = YLS(lmd,lma,lmb,L) &
                         + cfd * cfa * cfb * Y(lk1,lm1,ln1,L) * S(lk2,lm2,ln2)
                 enddo ! L
              enddo ! tb
              enddo ! ta
              enddo ! td
         enddo ! lmd
        enddo
     enddo
   end subroutine SHR_doPR3Ls

   function SHR_PR3Ls(LD,LA,LB,Y,S) result(YLS)
     ! Double Product Rule preserving orders Y(L):
     !
     ! couples "angular" factors
     !    Y(lm1,ln1,L)
     ! and
     !    S(lm2,ln2)
     ! performing double (L)-preserving product rule:
     !
     !    YLS(ma,mb,L) = PR(lma,lm1,lm2) PR(lmb,ln1,ln2)
     !                   *  Y(lm1,ln1,L) * S(lm2,ln2)
     !
     implicit none
     integer(IK), intent(in) :: LD,LA,LB
     real(RK), intent(in)    :: Y(:,:,:,:) ! ((LD+1)**2,(LA+1)**2,(LB+1)**2,1+LD+LA+LB)
     real(RK), intent(in)    :: S(:,:,:)   ! ((LD+1)**2,(LA+1)**2,(LB+1)**2        )
     real(RK)                :: YLS(2*LD+1,2*LA+1,2*LB+1,1+LD+LA+LB)
     ! *** end of interface ***

     YLS = 0.0_rk
     call SHR_doPR3Ls(LD**2+1,(LD+1)**2,LA**2+1,(LA+1)**2,LB**2+1,(LB+1)**2,Y,S,YLS)
   end function SHR_PR3Ls

#ifdef NO_BETTER_PLACE_FOR_IT
   ! find a better place for it:
   RECURSIVE function asort(X) result(Y)
     ! sorts in ascending order
     implicit none
     integer(IK), intent(in) :: X(:)
     integer(IK)             :: Y(size(X))
     ! *** end of interface ***

     integer(IK) :: n,m

     n = size(X)
     select case (n)
     case (0)
        ! do nothing
     case (1)
        ! return input
        Y = X
     case default
        m = n/2
        Y = ajoin( asort( X(:m) ), asort( X(m+1:) ) )
     end select

   contains

     function ajoin(A,B) result (Y)
       ! joins sorted arrays in to one sorted
       implicit none
       integer(IK), intent(in) :: A(:), B(:)
       integer(IK)             :: Y(size(A)+size(B))
       ! *** end of interface ***

       integer(IK) :: n,na,nb

       n  = 1
       na = 1
       nb = 1

       do while( na <= size(A) .and. nb <= size(B) )
          if( A(na) > B(nb) )then
             Y(n) = A(na)
             na = na + 1
          else
             Y(n) = B(nb)
             nb = nb + 1
          endif
          n = n + 1
       end do

       do while( na <= size(A) )
             Y(n) = A(na)
             na = na + 1
             n  = n  + 1
       end do

       do while( nb <= size(B) )
             Y(n) = B(nb)
             nb = nb + 1
             n  = n  + 1
       end do
     end function ajoin
   end function asort

   function perm(A,B) result(PM)
     !
     ! computes permutation suitable for reshape:
     !
     !  AMAT = reshape( BMAT, SHAPE(AMAT), order=PERM )
     !
     ! e.g.
     ! unsorted=  20  20  30  10 (A)
     !   sorted=  10  20  20  30 (B)
     !
     ! PMM: . T T .
     ! PMM: . T T .
     ! PMM: . . . T
     ! PMM: T . . .
     !
     ! perm    =   4   1   2   3  (NOT: 3   2   4   1)
     !
     ! { the element BMAT(2,1,1,1) goes to AMAT(1,1,1,2) }
     !
     ! Prefers more ``diagonal'' left to more ``unity'' right:
     !     PM: \ T . .            PM: \ . T .
     !     PM: . \ T .            PM: . T . .
     !     PM: . . \ T            PM: . . \ T
     !     PM: T . . \            PM: T . . \
     implicit none
     integer(IK), intent(in) :: A(:),B(:)
     integer(IK)             :: PM(size(B))
     ! *** end of interface ***

     integer(IK) :: i,j,NA,NB
     logical     :: PMM(size(A),size(B)), mask(size(A)), col(size(A))

     NA = size(A)
     NB = size(B)
     ! ASSERT(NA==NB)

     PM   = 0
     mask = .true.

     do i=1,NA
        do j=1,NB
           PMM(i,j) = A(i) .eq. B(j)
        enddo
     enddo

     do j=1,NB
        col = PMM(:,j) .and. mask
        ! find perm.:
        do i=1,NA
           if(.not.col(i)) cycle
           PM(j)   = i
           mask(i) = .false.
           exit ! this loop
        enddo
     enddo
   end function perm
#endif

  !--------------- End of module -------------------------------------
end module shgi_shr

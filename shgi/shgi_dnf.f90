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
module shgi_dnf
  !---------------------------------------------------------------
  !
  ! Computes harmonic derivatives
  !
  !   D(lma) ... D(lme) S(x^2/2)
  !
  ! of a scalar function.
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
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: doD2F
  public :: doD3F
  public :: doD4F
  public :: doD5F
  public :: doD2S
! public :: doD3S
  public :: doD4S
  public :: doD5S

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  ! ALL INTEGER CONSTANTS, LIKE GAX,GAY,GAZ ARE NOW IN
  !                     shgi_cntrl.f90
  ! THIS HAS BEEN DONE TO ENABLE ITS USE IN OTHER MODULES

  ! MOST GLOBAL VARIABLES HOLDING ``ANGULAR'' FACTORS WERE MOVED TO
  !                     shgi_common.f90
  ! THIS HAS BEEN DONE TO SPLIT THIS FILE INTO PARTS LATER

  ! a copy from shgi_shr.f90:
  integer(IK)   :: L_,M_
  integer(IK), parameter :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)
  integer(IK), parameter :: mof( (MAXL+1)**2 ) = (/((M_,M_=1,2*L_+1),L_=0,MAXL)/)

  !------------ Subroutines ---------------------------------------
contains

  subroutine doD2F(LA,LB,FL,YL,DF)
    !
    !   D2F = D(lma)D(lmb) F
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB
    real(RK)   , intent(in)  :: FL(:,:)     ! (NAB,1+SUM(L))
    real(RK)   , intent(in)  :: YL(:,:,:,:) ! (NAB,(LA+1)**2,(LB+1)**2,1+SUM(L))
    real(RK)   , intent(out) :: DF(:,:,:)   ! (NAB,(LA+1)**2,(LB+1)**2)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb

    ASSERT(1+LA+LB<=size(FL,2))
    ASSERT(1+LA+LB<=size(YL,4))

    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DF(:,lma,lmb) = 0.0
       do L=1+max(ila,ilb),1+ila+ilb
          DF(:,lma,lmb) =  DF(:,lma,lmb) &
              +      FL(:,L) * YL(:,lma,lmb,L)
       enddo
    enddo
    enddo
  end subroutine doD2F

  subroutine doD3F(LA,LB,LC,FL,YL,DF)
    !
    !   D3F = D(lma)D(lmb)D(lmc) F
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC
    real(RK)   , intent(in)  :: FL(:,:)       ! (NAB,1+SUM(L))
    real(RK)   , intent(in)  :: YL(:,:,:,:,:) ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,1+SUM(L))
    real(RK)   , intent(out) :: DF(:,:,:,:)   ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc

    ASSERT(1+LA+LB+LC<=size(FL,2))
    ASSERT(1+LA+LB+LC<=size(YL,5))

    do lmc=1,(LC+1)**2
       ilc=lof(lmc)
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DF(:,lma,lmb,lmc) = 0.0
       do L=1+max(ila,ilb,ilc),1+ila+ilb+ilc
          DF(:,lma,lmb,lmc) =  DF(:,lma,lmb,lmc) &
              +      FL(:,L) * YL(:,lma,lmb,lmc,L)
       enddo
    enddo
    enddo
    enddo
  end subroutine doD3F

  subroutine doD4F(LA,LB,LC,LD,FL,YL,DF)
    !
    !   D4F = D(lma)D(lmb)D(lmc)D(lmd) F
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC,LD
    real(RK)   , intent(in)  :: FL(:,:)         ! (NAB,1+SUM(L))
    real(RK)   , intent(in)  :: YL(:,:,:,:,:,:) ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,1+SUM(L))
    real(RK)   , intent(out) :: DF(:,:,:,:,:)   ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc
    integer(IK) :: lmd,ild

    ASSERT(1+LA+LB+LC+LD<=size(FL,2))
    ASSERT(1+LA+LB+LC+LD<=size(YL,6))

    do lmd=1,(LD+1)**2
       ild=lof(lmd)
    do lmc=1,(LC+1)**2
       ilc=lof(lmc)
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DF(:,lma,lmb,lmc,lmd) = 0.0
       do L=1+max(ila,ilb,ilc,ild),1+ila+ilb+ilc+ild
          DF(:,lma,lmb,lmc,lmd) =  DF(:,lma,lmb,lmc,lmd) &
              +          FL(:,L) * YL(:,lma,lmb,lmc,lmd,L)
       enddo
    enddo
    enddo
    enddo
    enddo
  end subroutine doD4F

  subroutine doD5F(LA,LB,LC,LD,LE,FL,YL,DF)
    !
    !   D5F = D(lma)D(lmb)D(lmc)D(lme) F
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC,LD,LE
    real(RK)   , intent(in)  :: FL(:,:)           ! (NAB,1+SUM(L))
    real(RK)   , intent(in)  :: YL(:,:,:,:,:,:,:) ! (NAB,(LA+1)**2,...,(LE+1)**2,1+SUM(L))
    real(RK)   , intent(out) :: DF(:,:,:,:,:,:)   ! (NAB,(LA+1)**2,...,(LE+1)**2)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc
    integer(IK) :: lmd,ild
    integer(IK) :: lme,ile

    ASSERT(1+LA+LB+LC+LD+LE<=size(FL,2))
    ASSERT(1+LA+LB+LC+LD+LE<=size(YL,7))

    do lme=1,(LE+1)**2
       ile=lof(lme)
    do lmd=1,(LD+1)**2
       ild=lof(lmd)
    do lmc=1,(LC+1)**2
       ilc=lof(lmc)
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DF(:,lma,lmb,lmc,lmd,lme) = 0.0
       do L=1+max(ila,ilb,ilc,ild,ile),1+ila+ilb+ilc+ild+ile
          DF(:,lma,lmb,lmc,lmd,lme) =  DF(:,lma,lmb,lmc,lmd,lme) &
                  +          FL(:,L) * YL(:,lma,lmb,lmc,lmd,lme,L)
       enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  end subroutine doD5F

  subroutine doD2S(LA,LB,SL,X,DS)
    !
    !   D2S = (-)^lb * D(lma)D(lmb) S
    !
    ! the (-)^lb factor may be interpreted as
    !
    !   Da(lma)Db(lmb)...
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB
    real(RK)   , intent(in)  :: SL(:,:)
    real(RK)   , intent(in)  :: X(:,:,:)
    real(RK)   , intent(out) :: DS(:,:,:)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb

    ASSERT(1+LA+LB<=size(SL,2))
    ASSERT(1+LA+LB<=size(X,3))

    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DS(:,lma,lmb) = 0.0
       do L=1+max(ila,ilb),1+ila+ilb
          DS(:,lma,lmb) =  DS(:,lma,lmb) &
              + (-1)**ilb * SL(:,L) * X(lma,lmb,L)
       enddo
    enddo
    enddo
  end subroutine doD2S

  subroutine doD3S(LA,LB,LC,SL,X,DS)
    !
    !   D3S = (-)^lb * D(lma)D(lmb)D(lmc) S
    !
    ! the (-)^lb factor may be interpreted as
    !
    !   Da(lma)Db(lmb)...
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC
    real(RK)   , intent(in)  :: SL(:,:)
    real(RK)   , intent(in)  :: X(:,:,:,:)
    real(RK)   , intent(out) :: DS(:,:,:,:)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc

    ASSERT(1+LA+LB+LC<=size(SL,2))
    ASSERT(1+LA+LB+LC<=size(X,4))

    do lmc=1,(LC+1)**2
       ilc=lof(lmc)
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DS(:,lma,lmb,lmc) = 0.0
       do L=1+max(ila,ilb,ilc),1+ila+ilb+ilc
          DS(:,lma,lmb,lmc) =  DS(:,lma,lmb,lmc) &
              + (-1)**ilb * SL(:,L) * X(lma,lmb,lmc,L)
       enddo
    enddo
    enddo
    enddo
  end subroutine doD3S

  subroutine doD4S(LA,LB,LC,LD,SL,X,DS)
    !
    !   D4S = (-)^lb * D(lma)D(lmb)D(lmc)D(lmd) S
    !
    ! the (-)^lb factor may be interpreted as
    !
    !   Da(lma)Db(lmb)...
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC,LD
    real(RK)   , intent(in)  :: SL(:,:)
    real(RK)   , intent(in)  :: X(:,:,:,:,:)
    real(RK)   , intent(out) :: DS(:,:,:,:,:)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc
    integer(IK) :: lmd,ild

    ASSERT(1+LA+LB+LC+LD<=size(SL,2))
    ASSERT(1+LA+LB+LC+LD<=size(X,5))

    do lmd=1,(LD+1)**2
       ild=lof(lmd)
    do lmc=1,(LC+1)**2
       ilc=lof(lmc)
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DS(:,lma,lmb,lmc,lmd) = 0.0
       do L=1+max(ila,ilb,ilc,ild),1+ila+ilb+ilc+ild
          DS(:,lma,lmb,lmc,lmd) =  DS(:,lma,lmb,lmc,lmd) &
              + (-1)**ilb * SL(:,L) * X(lma,lmb,lmc,lmd,L)
       enddo
    enddo
    enddo
    enddo
    enddo
  end subroutine doD4S

  subroutine doD5S(LA,LB,LC,LD,LE,SL,X,DS)
    !
    !   D5S = (-)^lb * D(lma)D(lmb)D(lmc)D(lmd)D(lme) S
    !
    ! the (-)^lb factor may be interpreted as
    !
    !   Da(lma)Db(lmb)...
    !
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC,LD,LE
    real(RK)   , intent(in)  :: SL(:,:)
    real(RK)   , intent(in)  :: X(:,:,:,:,:,:)
    real(RK)   , intent(out) :: DS(:,:,:,:,:,:)
    ! *** end of interface ***

    integer(IK) :: L
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc
    integer(IK) :: lmd,ild
    integer(IK) :: lme,ile

    ASSERT(1+LA+LB+LC+LD+LE<=size(SL,2))
    ASSERT(1+LA+LB+LC+LD+LE<=size(X,6))

    do lme=1,(LE+1)**2
       ile=lof(lme)
    do lmd=1,(LD+1)**2
       ild=lof(lmd)
    do lmc=1,(LC+1)**2
       ilc=lof(lmc)
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DS(:,lma,lmb,lmc,lmd,lme) = 0.0
       do L=1+max(ila,ilb,ilc,ild,ile),1+ila+ilb+ilc+ild+ile
          DS(:,lma,lmb,lmc,lmd,lme) =  DS(:,lma,lmb,lmc,lmd,lme) &
                + (-1)**ilb * SL(:,L) * X(lma,lmb,lmc,lmd,lme,L)
       enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  end subroutine doD5S

  !--------------- End of module ----------------------------------
end module shgi_dnf

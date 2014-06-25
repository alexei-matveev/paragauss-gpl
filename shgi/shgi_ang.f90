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
module shgi_ang
  !---------------------------------------------------------------
  !
  ! For  integrals  involving the  third  center  C precompute  common
  ! (angular) quantities.
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
! use CPU_TIME for timers:
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind, & ! type specification parameters
       I8K=>i8_kind
  use constants
  use shgi_cntrl
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: shgi_set_lcde
  public :: shgi_set_c

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

  !----------------------------------------------------------------
  ! a copy from shgi_shr.f90:
  integer(IK)   :: L_,M_
  integer(IK), parameter :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)
  integer(IK), parameter :: mof( (MAXL+1)**2 ) = (/((M_,M_=1,2*L_+1),L_=0,MAXL)/)
  !------------ Subroutines ---------------------------------------
contains

  subroutine shgi_set_lcde(ILC,ILD,ILE,ILF)
    ! (re)sets global LC and
    ! (re)allocates storage for the angular factors YL, YS, GS
    use shgi_common
    implicit none
    integer(IK), intent(in)  :: ILC,ILD,ILE
    integer(IK), intent(in), optional  :: ILF
    ! *** end of interface ***

    integer(IK) :: memstat
    integer(IK) :: CXX,CYY

    DPRINT  'SHGI: shgi_set_lcde(',ILC,ILD,ILE,')'

    ! set globals:
    LC   = ILC
    LD   = ILD
    LE   = ILE
    if(present(ILF)) LF   = ILF
    LMAX = max(LA,LB,LC,LD,LE,LF)

    if( allocated(YL) ) deallocate(YL)
    ! Allocation of the "angular" component
    allocate(YL(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+LA+LB+LC+LD+LE),&
         stat=memstat)
    ASSERT(memstat==0)

    ! YS for "special" harmonic types for SR and SO
    ! the last index of YS goes from minus something
    ! to 1=C00 (maybe later even higher)
    !
    ! The default range is, just C00:C00, is enough
    ! for (NR) nuclear, s- and r2-fit integrals
    !
    CXX = C00
    CYY = C00

    if(is_on(IYSNU+IYSSR+IYSSO))then
    if(is_on(IYSSR))then
       ! FIXME: in case you change the order
       CXX = min(CP2,CR2)
       CYY = C00
    endif

    if(is_on(IYSSO))then
       ! FIXME: in case you change the order
       CXX = min(CVX,CVY,CVZ,CP2,CR2)
       CYY = C00
    endif

    DPRINT 'SHGI: allocating YS from ',CXX,' to ',CYY

    if( allocated(YS) ) deallocate(YS)
    allocate(YS(NAB,2*LA+1,2*LB+1,1+LA+LB+LC,CXX:CYY),stat=memstat)
    ASSERT(memstat==0)
    endif

    if(is_on(IGSNU+IGSSR))then
       ASSERT(LC>=1)
       CXX = min(GWX,GWY,GWZ,GDX,GDY,GDZ)
       CYY = max(GWX,GWY,GWZ,GDX,GDY,GDZ)
       if(is_on(IGSSR))then
          ASSERT(LD>=1)
          ! raise the upper bound:
          CYY = CYY + 12
       endif
       if( allocated(GS) ) deallocate(GS)
       allocate(GS(NAB,2*LA+1,2*LB+1,1+LA+LB+LC+LD,CXX:CYY),stat=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine shgi_set_lcde

  subroutine shgi_set_c(xc)
    use shgi_common, only: W, WC, W2, YL, YS, GS, S5, K4
    use shgi_common, only: WDAL, WDBL
    use shgi_shr, only: SHR_D5Av, SHR_YLSv
    implicit none
    real(RK)   , intent(in) :: xc(3)
    ! *** end of interface ***

    !eal(RK) :: YL(NAB,(LA+1)**2,...,(LE+1)**2,1+SUM(L))

    integer(IK) :: iab
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc,ilc
    integer(IK) :: lmd,ild
    integer(IK) :: lme,ile
    integer(IK) :: L

    FPP_TIMER_START(tsc)

    DPRINT 'shgi_set_c:',xc

    do iab=1,NAB
       WC(iab,:) =  W(iab,:) - xc
       W2(iab  ) = WC(iab,1)**2 + WC(iab,2)**2 + WC(iab,3)**2
    enddo

!   YL = HUGE(1.0_rk) ! FIXME: see if propagates
!   YL = NaN()        ! FIXME: see if propagates
    FPP_TIMER_START(td3a)
    call SHR_D5Av(NAB,LA,LB,LC,LD,LE,WC,YL)
    FPP_TIMER_STOP(td3a)

    ! D(lma)D(lmb) D(lmc)I -> DA(lma)DB(lmb) D(lmc)I
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

       do L=1+max(ila,ilb,ilc,ild,ile),1+ila+ilb+ilc+ild+ile

          YL(:,lma,lmb,lmc,lmd,lme,L) = YL(:,lma,lmb,lmc,lmd,lme,L) &
               * WDAL(:,1+ila) * WDBL(:,1+ilb)
       enddo
    enddo
    enddo
    enddo
    enddo
    enddo

    FPP_TIMER_START(tyls)

    if(is_on(IYSNU))then
    ! needed for NUC, and FIT integrals
    YS(:,:,:,:,C00) = SHR_YLSv(LA,LB,0,YL(:,:,:,C0,C0,C0,:),S5(:,:,:,C0,C0,C0))
    endif ! IYSNU

    if(is_on(IYSSR))then
    ASSERT(LC>=1)
    ! only needed for SR integrals:
    ! PR(r2) : F(x)*S(x) + F(y)*S(y) + F(z)*S(z)
    YS(:,:,:,:,CR2) = SHR_YLSv(LA,LB,1,YL(:,:,:,CX,C0,C0,:),S5(:,:,:,CX,C0,C0)) &
                    + SHR_YLSv(LA,LB,1,YL(:,:,:,CY,C0,C0,:),S5(:,:,:,CY,C0,C0)) &
                    + SHR_YLSv(LA,LB,1,YL(:,:,:,CZ,C0,C0,:),S5(:,:,:,CZ,C0,C0))
    ! PR(p2) : F(0)*Sp2(0)
    YS(:,:,:,:,CP2) = SHR_YLSv(LA,LB,0,YL(:,:,:,C0,C0,C0,:),K4(:,:,:,C0,C0)   )
    endif ! IYSSR

    if(is_on(IYSSO))then
    ASSERT(LC>=1)
    ! only needed for SO integrals:
    ! PR(sox) : F(y)*S(z) - F(z)*S(y)
    ! PR(soy) : F(z)*S(x) - F(x)*S(z)
    ! PR(soz) : F(x)*S(y) - F(y)*S(x)
    YS(:,:,:,:,CVX) = SHR_YLSv(LA,LB,1,YL(:,:,:,CY,C0,C0,:),S5(:,:,:,CZ,C0,C0)) &
                    - SHR_YLSv(LA,LB,1,YL(:,:,:,CZ,C0,C0,:),S5(:,:,:,CY,C0,C0))

    YS(:,:,:,:,CVY) = SHR_YLSv(LA,LB,1,YL(:,:,:,CZ,C0,C0,:),S5(:,:,:,CX,C0,C0)) &
                    - SHR_YLSv(LA,LB,1,YL(:,:,:,CX,C0,C0,:),S5(:,:,:,CZ,C0,C0))

    YS(:,:,:,:,CVZ) = SHR_YLSv(LA,LB,1,YL(:,:,:,CX,C0,C0,:),S5(:,:,:,CY,C0,C0)) &
                    - SHR_YLSv(LA,LB,1,YL(:,:,:,CY,C0,C0,:),S5(:,:,:,CX,C0,C0))
    endif ! IYSSO

    if(is_on(IGSNU))then
    ASSERT(LC>=1)
!   GS = HUGE(1.0_rk) ! FIXME: see if it propagates
!   GS = NaN()        ! FIXME: see if it propagates

    ! GS will contain angular factors for the two independent gradients
    !
    !                         !!!  GW and GD !!!
    !
    ! once obtained GW and GD grads may be converted to GA, GB, and GC
    ! by a call to shgi_gr_wd_to_ab(gwd(..,:6),gab(..,:9))

    ! GR(wx) : F(x)*S(0)
    ! GR(wy) : F(y)*S(0)
    ! GR(wz) : F(z)*S(0)
    GS(:,:,:,:,   GWX) = SHR_YLSv( LA,LB,1, YL(:,:,:,CX,C0,C0,:), S5(:,:,:,C0,C0,C0))
    GS(:,:,:,:,   GWY) = SHR_YLSv( LA,LB,1, YL(:,:,:,CY,C0,C0,:), S5(:,:,:,C0,C0,C0))
    GS(:,:,:,:,   GWZ) = SHR_YLSv( LA,LB,1, YL(:,:,:,CZ,C0,C0,:), S5(:,:,:,C0,C0,C0))

    ! GR(dx) : F(0)*S(x)
    ! GR(dy) : F(0)*S(y)
    ! GR(dz) : F(0)*S(z)
    GS(:,:,:,:,   GDX) = SHR_YLSv( LA,LB,0, YL(:,:,:,C0,C0,C0,:), S5(:,:,:,CX,C0,C0))
    GS(:,:,:,:,   GDY) = SHR_YLSv( LA,LB,0, YL(:,:,:,C0,C0,C0,:), S5(:,:,:,CY,C0,C0))
    GS(:,:,:,:,   GDZ) = SHR_YLSv( LA,LB,0, YL(:,:,:,C0,C0,C0,:), S5(:,:,:,CZ,C0,C0))
    endif ! IGSNU

    if(is_on(IGSSR))then
    ASSERT(LD>=1)
    ! PR(r2) : F(x)*S(x) + F(y)*S(y) + F(z)*S(z)
    GS(:,:,:,:, 6+GWX) = SHR_YLSv( LA,LB,2, YL(:,:,:,CX,CX,C0,:), S5(:,:,:,CX,C0,C0)) &
                       + SHR_YLSv( LA,LB,2, YL(:,:,:,CY,CX,C0,:), S5(:,:,:,CY,C0,C0)) &
                       + SHR_YLSv( LA,LB,2, YL(:,:,:,CZ,CX,C0,:), S5(:,:,:,CZ,C0,C0))

    GS(:,:,:,:, 6+GWY) = SHR_YLSv( LA,LB,2, YL(:,:,:,CX,CY,C0,:), S5(:,:,:,CX,C0,C0)) &
                       + SHR_YLSv( LA,LB,2, YL(:,:,:,CY,CY,C0,:), S5(:,:,:,CY,C0,C0)) &
                       + SHR_YLSv( LA,LB,2, YL(:,:,:,CZ,CY,C0,:), S5(:,:,:,CZ,C0,C0))

    GS(:,:,:,:, 6+GWZ) = SHR_YLSv( LA,LB,2, YL(:,:,:,CX,CZ,C0,:), S5(:,:,:,CX,C0,C0)) &
                       + SHR_YLSv( LA,LB,2, YL(:,:,:,CY,CZ,C0,:), S5(:,:,:,CY,C0,C0)) &
                       + SHR_YLSv( LA,LB,2, YL(:,:,:,CZ,CZ,C0,:), S5(:,:,:,CZ,C0,C0))

    GS(:,:,:,:, 6+GDX) = SHR_YLSv( LA,LB,1, YL(:,:,:,CX,C0,C0,:), S5(:,:,:,CX,CX,C0)) &
                       + SHR_YLSv( LA,LB,1, YL(:,:,:,CY,C0,C0,:), S5(:,:,:,CY,CX,C0)) &
                       + SHR_YLSv( LA,LB,1, YL(:,:,:,CZ,C0,C0,:), S5(:,:,:,CZ,CX,C0))

    GS(:,:,:,:, 6+GDY) = SHR_YLSv( LA,LB,1, YL(:,:,:,CX,C0,C0,:), S5(:,:,:,CX,CY,C0)) &
                       + SHR_YLSv( LA,LB,1, YL(:,:,:,CY,C0,C0,:), S5(:,:,:,CY,CY,C0)) &
                       + SHR_YLSv( LA,LB,1, YL(:,:,:,CZ,C0,C0,:), S5(:,:,:,CZ,CY,C0))

    GS(:,:,:,:, 6+GDZ) = SHR_YLSv( LA,LB,1, YL(:,:,:,CX,C0,C0,:), S5(:,:,:,CX,CZ,C0)) &
                       + SHR_YLSv( LA,LB,1, YL(:,:,:,CY,C0,C0,:), S5(:,:,:,CY,CZ,C0)) &
                       + SHR_YLSv( LA,LB,1, YL(:,:,:,CZ,C0,C0,:), S5(:,:,:,CZ,CZ,C0))

    ! PR(p2) : F(0)*Sp2(0)
    GS(:,:,:,:,12+GWX) = SHR_YLSv( LA,LB,1, YL(:,:,:,C0,CX,C0,:), K4(:,:,:,C0,C0)   )
    GS(:,:,:,:,12+GWY) = SHR_YLSv( LA,LB,1, YL(:,:,:,C0,CY,C0,:), K4(:,:,:,C0,C0)   )
    GS(:,:,:,:,12+GWZ) = SHR_YLSv( LA,LB,1, YL(:,:,:,C0,CZ,C0,:), K4(:,:,:,C0,C0)   )
    ASSERT(size(K4,4)>=4)
    GS(:,:,:,:,12+GDX) = SHR_YLSv( LA,LB,0, YL(:,:,:,C0,C0,C0,:), K4(:,:,:,CX,C0)   )
    GS(:,:,:,:,12+GDY) = SHR_YLSv( LA,LB,0, YL(:,:,:,C0,C0,C0,:), K4(:,:,:,CY,C0)   )
    GS(:,:,:,:,12+GDZ) = SHR_YLSv( LA,LB,0, YL(:,:,:,C0,C0,C0,:), K4(:,:,:,CZ,C0)   )
    endif ! IGSSR

    FPP_TIMER_STOP(tyls)
    FPP_TIMER_STOP(tsc)
  end subroutine shgi_set_c

  !--------------- End of module ----------------------------------
end module shgi_ang

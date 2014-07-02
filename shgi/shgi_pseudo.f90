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
module shgi_pseudo
  !-------------------------------------------------------------------
  !
  ! One-electron  integrals for  local and  non-local pseudo-potential
  ! terms.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
  ! Copyright (c) 2006 Vladimir Nasluzov
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
! use CPU_TIME for timers:
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  use constants
  use shgi_cntrl
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: shgi_pseu_set_ab
  public :: shgi_pseu_set_abc
  public :: shgi_pseu_close_abc
  public :: shgi_pseu
  public :: shgi_gr_pseu
  public :: shgi_sd_pseu

! public :: AEXP, BEXP, XA, XB

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  ! ALL INTEGER CONSTANTS, LIKE GAX,GAY,GAZ ARE NOW IN
  !                     shgi_cntrl.f90
  ! THIS HAS BEEN DONE TO ENABLE ITS USE IN OTHER MODULES

  ! MOST GLOBAL VARIABLES HOLDING ``ANGULAR'' FACTORS WERE MOVED TO
  !                     shgi_common.f90
  ! THIS HAS BEEN DONE TO SPLIT THIS FILE INTO PARTS LATER

  real(RK)              ::        XAC(3), XAC2         !               for PP
  real(RK)              ::        XBC(3), XBC2         !               for PP

! real(RK), allocatable :: AEXP(:), BEXP(:)            ! (NEA), (NEB), for PP

  ! only for PP:
  real(RK), allocatable :: X2P(:,:,:,:,:)              ! (    2*LA+1,2*LB+1,1+LA,1+LB,1+LP)
  ! for PP gradients:
  real(RK), allocatable :: G2P(:,:,:,:,:,:)            ! (    2*LA+1,2*LB+1,2+LA,2+LB,1+LP,9)
  real(RK), allocatable :: D2P(:,:,:,:,:,:,:)          ! (    2*LA+1,2*LB+1,3+LA,3+LB,1+LP,9,9)

  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi_pseu_set_ab()
    use bessel_module, only: bessel_setup
    implicit none
    ! *** end of interface ***

    call bessel_setup()
    ! bessel_close() does nothing important
  end subroutine shgi_pseu_set_ab

  subroutine shgi_pseu_set_abc(LP,order)
    implicit none
    integer(IK), intent(in) :: LP
    integer(IK), intent(in) :: order
    ! *** end of interface ***

    integer(IK) :: memstat

    select case(order)

    case (0) ! Energy
      ! deallocate in shgi_pseu_close_abc(),
      ! LP maybe different for different C, or not required at all
      allocate(X2P(2*LA+1,2*LB+1,1+LA,1+LB,1+LP),stat=memstat)
      ASSERT(memstat==0)

    case (1) ! Gradients
      allocate(G2P(2*LA+1,2*LB+1,2+LA,2+LB,1+LP,9),stat=memstat)
      ASSERT(memstat==0)

    case (2) ! SecDers
      allocate(D2P(2*LA+1,2*LB+1,3+LA,3+LB,1+LP,9,9),stat=memstat)
      ASSERT(memstat==0)

    case default
       ABORT('no such case')
    end select
  end subroutine shgi_pseu_set_abc

  subroutine shgi_pseu_close_abc()
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat

    ! allocated in shgi_pseu_set_abc(),
    ! LP maybe different for different C, or not required at all
    if( allocated(X2P) )then
       deallocate(X2P,stat=memstat)
       ASSERT(memstat==0)
    endif

    ! for gradients:
    if( allocated(G2P) )then
       deallocate(G2P,stat=memstat)
       ASSERT(memstat==0)
    endif

    ! for derivatives:
    if(allocated(D2P) )then
      deallocate(D2P,stat=memstat)
      ASSERT(memstat==0)
    endif
  end subroutine shgi_pseu_close_abc

  subroutine shgi_pseu_set_c(xc,LP,LG)
    use shgi_shr, only: SHR_D3Av, SHR_PR2Ls
    use shgi_common, only: XA,XB
    implicit none
    real(RK)   , intent(in) :: xc(3)
    integer(IK), intent(in) :: LP
    integer(IK), intent(in) :: LG ! LG=1 for gradients
    ! *** end of interface ***

    integer(IK) :: LAP,LBP
    real(RK)    :: X3A((LG+1)**2,(LA+1)**2,(LP+1)**2,1+LG+LA+LP)
    real(RK)    :: X3B((LG+1)**2,(LB+1)**2,(LP+1)**2,1+LG+LB+LP)
    real(RK)    :: PRA(2*LG+1   ,2*LA+1             ,1+LG+LA   )
    real(RK)    :: PRB(2*LG+1   ,2*LB+1             ,1+LG+LB   )
    integer(IK) :: llp,mp,lmp
    integer(IK) :: lla,ma
    integer(IK) :: llb,mb

    FPP_TIMER_START(tpsc)

    DPRINT 'shgi_pseu_set_c: LP=',LP,' XC=',xc

    XAC  = XA - xc
    XBC  = XB - xc
    XAC2 = XAC(1)**2 + XAC(2)**2 + XAC(3)**2
    XBC2 = XBC(1)**2 + XBC(2)**2 + XBC(3)**2

    LAP = max(LG,LA,LP)
    LBP = max(LG,LB,LP)

    ! FIXME: only X3[AB](:,:,LPMP,1+LP) is needed:
    call SHR_D3Av(1,LG,LA,LP,XAC,X3A)
    call SHR_D3Av(1,LG,LB,LP,XBC,X3B)

    ! allocated by shgi_pseu_set_abc(LP), once for
    ! all EAs:
    ASSERT(allocated(X2P))

    X2P = zero
    lmp = 0
    do llp=0,LP
       do mp=1,2*llp+1
          lmp = lmp + 1
          ! FIXME: if llp < LG+LA it is zero anyway.
          ! FIXME: X3[AB](:,:,C00,:1+LG+L[AB]) as the first arg.
          PRA(:,:,:) = SHR_PR2Ls(LG,LA,X3A(:,:,C00,:),X3A(:,:,lmp,1+llp))
          PRB(:,:,:) = SHR_PR2Ls(LG,LB,X3B(:,:,C00,:),X3B(:,:,lmp,1+llp))

          ! for energy calc (LG=0):
          forall(ma=1:2*LA+1, mb=1:2*LB+1, lla=0:LA, llb=0:LB)
             X2P(ma,mb,1+lla,1+llb,1+llp) = X2P(ma,mb,1+lla,1+llb,1+llp) &
                  + PRA(C00,ma,1+lla) * PRB(C00,mb,1+llb)
          end forall
       enddo
    enddo

#if _DPRINT
    DPRINT 'shgi_pseu_set_c: === X2P ==='
    DPRINT 'shgi_pseu_set_c: *** upto la`=',LA,' lb`=',LB,' LP=',LP
    do llp=0,LP
       do llb=0,LB
          do lla=0,LA
             DPRINT 'shgi_pseu_set_c: +++ X(',lla,llb,')[lp=',llp,']'
             do ma=1,2*LA+1
                DPRINT 'shgi_pseu_set_c: ',X2P(ma,:,1+lla,1+llb,1+llp)
             enddo
          enddo
       enddo
    enddo
#endif
    FPP_TIMER_STOP(tpsc)
  end subroutine shgi_pseu_set_c

  subroutine shgi_gr_pseu_set_c(xc,LP)
    use shgi_shr, only: SHR_D3Av, SHR_doPR2Ls
    use shgi_common, only: XA,XB
    USE_DEBUG, only: show
    implicit none
    real(RK)   , intent(in) :: xc(3)
    integer(IK), intent(in) :: LP
    ! *** end of interface ***

    integer(IK), parameter :: LG=1
    integer(IK) :: LAP,LBP
    real(RK)    :: X3A((LG+1)**2,(LA+1)**2,(LP+1)**2,1+LG+LA+LP)
    real(RK)    :: X3B((LG+1)**2,(LB+1)**2,(LP+1)**2,1+LG+LB+LP)
    real(RK)    :: PRA((LG+1)**2,2*LA+1             ,1+LG+LA   )
    real(RK)    :: PRB((LG+1)**2,2*LB+1             ,1+LG+LB   )
    integer(IK) :: ilp,mp,llp,lmp
    integer(IK) :: ila,ma
    integer(IK) :: ilb,mb

    FPP_TIMER_START(tpsc)

    DPRINT  'shgi_gr_pseu_set_c: LP=',LP,' XC=',xc

    XAC  = XA - xc
    XBC  = XB - xc
    XAC2 = XAC(1)**2 + XAC(2)**2 + XAC(3)**2
    XBC2 = XBC(1)**2 + XBC(2)**2 + XBC(3)**2

    LAP = max(LG,LA,LP)
    LBP = max(LG,LB,LP)

    ! FIXME: only X3[AB](:,:,LPMP,1+LP) is needed:
    call SHR_D3Av(1,LG,LA,LP,XAC,X3A)
    call SHR_D3Av(1,LG,LB,LP,XBC,X3B)

    ! allocated by shgi_pseu_set_abc(LP), once for all EAs:
    ASSERT(allocated(G2P))
    G2P = zero

    lmp = 0
    do LLP=0,LP
       do mp=1,2*LLP+1
          lmp = lmp + 1
          ilp = LLP + 1 ! index in storage
          ! FIXME: if LLP < LG+LA it is zero anyway.
          ! FIXME: X3[AB](:,:,C00,:1+LG+L[AB]) as the first arg.
          PRA = zero
          PRB = zero
          call SHR_doPR2Ls(1,(LG+1)**2,LA**2+1,(LA+1)**2,X3A(:,:,C00,:),X3A(:,:,lmp,ilp),PRA)
          call SHR_doPR2Ls(1,(LG+1)**2,LB**2+1,(LB+1)**2,X3B(:,:,C00,:),X3B(:,:,lmp,ilp),PRB)

          ! if moving A:
          forall(ma=1:2*LA+1, mb=1:2*LB+1, ila=1:LA+2, ilb=1:LB+1  )
             G2P(ma,mb,ila,ilb,ilp,GAX) = G2P(ma,mb,ila,ilb,ilp,GAX) + PRA(C1X,ma,ila) * PRB(C00,mb,ilb)
             G2P(ma,mb,ila,ilb,ilp,GAY) = G2P(ma,mb,ila,ilb,ilp,GAY) + PRA(C1Y,ma,ila) * PRB(C00,mb,ilb)
             G2P(ma,mb,ila,ilb,ilp,GAZ) = G2P(ma,mb,ila,ilb,ilp,GAZ) + PRA(C1Z,ma,ila) * PRB(C00,mb,ilb)
          end forall

          ! if moving B:
          forall(ma=1:2*LA+1, mb=1:2*LB+1, ila=1:LA+1  , ilb=1:LB+2)
             G2P(ma,mb,ila,ilb,ilp,GBX) = G2P(ma,mb,ila,ilb,ilp,GBX) + PRA(C00,ma,ila) * PRB(C1X,mb,ilb)
             G2P(ma,mb,ila,ilb,ilp,GBY) = G2P(ma,mb,ila,ilb,ilp,GBY) + PRA(C00,ma,ila) * PRB(C1Y,mb,ilb)
             G2P(ma,mb,ila,ilb,ilp,GBZ) = G2P(ma,mb,ila,ilb,ilp,GBZ) + PRA(C00,ma,ila) * PRB(C1Z,mb,ilb)
          end forall

       enddo ! sum over MP within LP
          ! if moving C:
          G2P(:,:,:,:,ilp,GCX) = - G2P(:,:,:,:,ilp,GAX) - G2P(:,:,:,:,ilp,GBX)
          G2P(:,:,:,:,ilp,GCY) = - G2P(:,:,:,:,ilp,GAY) - G2P(:,:,:,:,ilp,GBY)
          G2P(:,:,:,:,ilp,GCZ) = - G2P(:,:,:,:,ilp,GAZ) - G2P(:,:,:,:,ilp,GBZ)
    enddo

    FPP_TIMER_STOP(tpsc)
  end subroutine shgi_gr_pseu_set_c

  subroutine shgi_sd_pseu_set_c(xc,LP)
    use shgi_shr, only: SHR_doPR2Ls, SHR_D4Av, SHR_doPR3Ls
    use shgi_common, only: XA,XB
    USE_DEBUG, only: show
    implicit none
    real(RK)   , intent(in) :: xc(3)
    integer(IK), intent(in) :: LP
    ! *** end of interface ***

    integer(IK), parameter :: LG=1
    integer(IK), parameter :: LD=1
    integer(IK) :: LAP,LBP
    real(RK)    :: X4A((LD+1)**2,(LG+1)**2,(LA+1)**2,(LP+1)**2,1+LD+LG+LA+LP)
    real(RK)    :: X4B((LD+1)**2,(LG+1)**2,(LB+1)**2,(LP+1)**2,1+LD+LG+LB+LP)
    real(RK)    :: PRA((LD+1)**2,(LG+1)**2,2*LA+1             ,1+LD+LG+LA   )
    real(RK)    :: PRB((LD+1)**2,(LG+1)**2,2*LB+1             ,1+LD+LG+LB   )
    integer(IK) :: ilp,mp,llp,lmp
    integer(IK) :: ila,ma
    integer(IK) :: ilb,mb

    FPP_TIMER_START(tpsc)

    DPRINT 'shgi_sd_pseu_set_c: LP=',LP,' XC=',xc

    XAC  = XA - xc
    XBC  = XB - xc
    XAC2 = XAC(1)**2 + XAC(2)**2 + XAC(3)**2
    XBC2 = XBC(1)**2 + XBC(2)**2 + XBC(3)**2

    LAP = max(LD,LG,LA,LP)
    LBP = max(LD,LG,LB,LP)

    ! FIXME: only X4[AB](:,:,LPMP,1+LP) is needed:
    call SHR_D4Av(1,LD,LG,LA,LP,XAC,X4A)
    call SHR_D4Av(1,LD,LG,LB,LP,XBC,X4B)

    ! allocated by shgi_pseu_set_abc(LP), once for all EAs:
    ASSERT(allocated(D2P))
    D2P = zero

    lmp = 0
    do LLP=0,LP
       do mp=1,2*LLP+1
          lmp = lmp + 1
          ilp = LLP + 1 ! index in storage
          ! FIXME: if LLP < LG+LA it is zero anyway.
          ! FIXME: X3[AB](:,:,C00,:1+LG+L[AB]) as the first arg.
          PRA = zero
          PRB = zero
          call SHR_doPR3Ls(1,(LD+1)**2,1,(LG+1)**2,LA**2+1,(LA+1)**2,X4A(:,:,:,C00,:),X4A(:,:,:,lmp,ilp),PRA)
          call SHR_doPR3Ls(1,(LD+1)**2,1,(LG+1)**2,LB**2+1,(LB+1)**2,X4B(:,:,:,C00,:),X4B(:,:,:,lmp,ilp),PRB)

          ! if moving A:
          forall(ma=1:2*LA+1, mb=1:2*LB+1, ila=1:LA+3, ilb=1:LB+1  )
             D2P(ma,mb,ila,ilb,ilp,GAX,GAX) = D2P(ma,mb,ila,ilb,ilp,GAX,GAX) + PRA(C1X,C1X,ma,ila) * PRB(C00,C00,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAX,GAY) = D2P(ma,mb,ila,ilb,ilp,GAX,GAY) + PRA(C1X,C1Y,ma,ila) * PRB(C00,C00,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAX,GAZ) = D2P(ma,mb,ila,ilb,ilp,GAX,GAZ) + PRA(C1X,C1Z,ma,ila) * PRB(C00,C00,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAY,GAY) = D2P(ma,mb,ila,ilb,ilp,GAY,GAY) + PRA(C1Y,C1Y,ma,ila) * PRB(C00,C00,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAY,GAZ) = D2P(ma,mb,ila,ilb,ilp,GAY,GAZ) + PRA(C1Y,C1Z,ma,ila) * PRB(C00,C00,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAZ,GAZ) = D2P(ma,mb,ila,ilb,ilp,GAZ,GAZ) + PRA(C1Z,C1Z,ma,ila) * PRB(C00,C00,mb,ilb)
          end forall

          ! if moving B:
          forall(ma=1:2*LA+1, mb=1:2*LB+1, ila=1:LA+1  , ilb=1:LB+3)
             D2P(ma,mb,ila,ilb,ilp,GBX,GBX) = D2P(ma,mb,ila,ilb,ilp,GBX,GBX) + PRA(C00,C00,ma,ila) * PRB(C1X,C1X,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GBX,GBY) = D2P(ma,mb,ila,ilb,ilp,GBX,GBY) + PRA(C00,C00,ma,ila) * PRB(C1X,C1Y,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GBX,GBZ) = D2P(ma,mb,ila,ilb,ilp,GBX,GBZ) + PRA(C00,C00,ma,ila) * PRB(C1X,C1Z,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GBY,GBY) = D2P(ma,mb,ila,ilb,ilp,GBY,GBY) + PRA(C00,C00,ma,ila) * PRB(C1Y,C1Y,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GBY,GBZ) = D2P(ma,mb,ila,ilb,ilp,GBY,GBZ) + PRA(C00,C00,ma,ila) * PRB(C1Y,C1Z,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GBZ,GBZ) = D2P(ma,mb,ila,ilb,ilp,GBZ,GBZ) + PRA(C00,C00,ma,ila) * PRB(C1Z,C1Z,mb,ilb)
          end forall

          forall(ma=1:2*LA+1, mb=1:2*LB+1, ila=1:LA+2  , ilb=1:LB+2)
             D2P(ma,mb,ila,ilb,ilp,GAX,GBX) = D2P(ma,mb,ila,ilb,ilp,GAX,GBX) + PRA(C00,C1X,ma,ila) * PRB(C00,C1X,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAX,GBY) = D2P(ma,mb,ila,ilb,ilp,GAX,GBY) + PRA(C00,C1X,ma,ila) * PRB(C00,C1Y,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAX,GBZ) = D2P(ma,mb,ila,ilb,ilp,GAX,GBZ) + PRA(C00,C1X,ma,ila) * PRB(C00,C1Z,mb,ilb)

             D2P(ma,mb,ila,ilb,ilp,GAY,GBX) = D2P(ma,mb,ila,ilb,ilp,GAY,GBX) + PRA(C00,C1Y,ma,ila) * PRB(C00,C1X,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAY,GBY) = D2P(ma,mb,ila,ilb,ilp,GAY,GBY) + PRA(C00,C1Y,ma,ila) * PRB(C00,C1Y,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAY,GBZ) = D2P(ma,mb,ila,ilb,ilp,GAY,GBZ) + PRA(C00,C1Y,ma,ila) * PRB(C00,C1Z,mb,ilb)

             D2P(ma,mb,ila,ilb,ilp,GAZ,GBX) = D2P(ma,mb,ila,ilb,ilp,GAZ,GBX) + PRA(C00,C1Z,ma,ila) * PRB(C00,C1X,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAZ,GBY) = D2P(ma,mb,ila,ilb,ilp,GAZ,GBY) + PRA(C00,C1Z,ma,ila) * PRB(C00,C1Y,mb,ilb)
             D2P(ma,mb,ila,ilb,ilp,GAZ,GBZ) = D2P(ma,mb,ila,ilb,ilp,GAZ,GBZ) + PRA(C00,C1Z,ma,ila) * PRB(C00,C1Z,mb,ilb)
          end forall

       enddo ! sum over MP within LP
    enddo
    FPP_TIMER_STOP(tpsc)
  end subroutine shgi_sd_pseu_set_c


  subroutine shgi_pseu(ZCORE,LP,pps,xc,PIPSEU)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    USE_DEBUG, only: show
    implicit none
    real(RK)   , intent(in)    :: ZCORE
    integer(IK), intent(in)    :: LP
    type(ppt)  , intent(in)    :: pps(0:) ! (0:LP+1)
    real(RK)   , intent(in)    :: xc(3)
    real(RK)   , intent(inout) :: PIPSEU(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: LPLOW,LPHIGH,llp

    ASSERT(size(pps)==2+LP)

    ! first the nonlocal part:
    LPLOW  = 0
    LPHIGH = LP ! in general all PPs contribute

    if(is_on(IAEQC))then
       ! A==C then only LA==LP may contibute:
       LPLOW  = max(LA,LPLOW)
       LPHIGH = min(LA,LPHIGH)
    endif

    if(is_on(IBEQC))then
       ! B==C then only LB==LP may contibute:
       LPLOW  = max(LB,LPLOW)
       LPHIGH = min(LB,LPHIGH)
    endif
    ! if AEQC and BEQC then
    ! LPLOW  == max(LA,LB)
    ! LPHIGH == min(LA,LB,LP)
    ! which assures that only LA==LB <= LP
    ! are computed
    DPRINT 'SHGI: only doing LPLOW=',LPLOW,' .. LPHIGH=',LPHIGH

    DCALL show('pp+0'  ,PIPSEU(1,:,:))
    if( LPLOW <= LPHIGH )then
       call shgi_pseu_set_c(xc,LP,0) ! LG=0 i.e. no gradients
       do llp=LPHIGH,LPLOW,-1
          ! loop range 0:LP in general, or LA:LA, or LB:LB, or not at all
          ! sum the smallest terms first, from HIGH LP to LOW LP
          call shgi_pseuL(llp,pps(llp),PIPSEU)
          DCALL show('pp+lp',PIPSEU(1,:,:))
       enddo
    else
       DPRINT 'SHGI: only local contribution, LA=',LA,' LB=',LB,' LP=',LP
    endif

    ! second the local part, as maybe the largest contribution?:
    call shgi_pseu0(ZCORE,xc,pps(LP+1),PIPSEU)
    DCALL show('pp+loc',PIPSEU(1,:,:))
  end subroutine shgi_pseu

  subroutine shgi_gr_pseu(ZCORE,LP,pps,xc,PIPSEU)
    use type_module, only: i8_kind
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    USE_DEBUG, only: show
    implicit none
    real(RK)   , intent(in)    :: ZCORE
    integer(IK), intent(in)    :: LP
    type(ppt)  , intent(in)    :: pps(0:) ! (0:LP+1)
    real(RK)   , intent(in)    :: xc(3)
    real(RK)   , intent(inout) :: PIPSEU(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    ! *** end of interface ***

    integer(IK) :: LPLOW,LPHIGH,llp
    integer (i8_kind) :: same

    ASSERT(size(pps)==2+LP)

    ! first the nonlocal part:
    same = whatis(IXEQC)
    select case(same)
    case (0) ! all centers differ
       LPLOW  = 0
       LPHIGH = LP ! in general all PPs contribute
    case (IAEQC) ! A=C
       ! A==C then only LA==LP may contibute:
       LPLOW  = LA
       LPHIGH = LA
    case (IBEQC) ! B=C
       ! B==C then only LB==LP may contibute:
       LPLOW  = LB
       LPHIGH = LB
    case (IXEQC) ! A=B=C
       ! same center, no forces at all
       WARN('A=B=C, why calling?')
       return !?
    case default
       ABORT('what?')
    end select
    DPRINT 'SHGI: only doing LPLOW=',LPLOW,' .. LPHIGH=',LPHIGH

!!$    DCALL show('pp+0'  ,PIPSEU(1,:,:))
    if( LPLOW <= LP )then
       call shgi_gr_pseu_set_c(xc,LP) ! calc G2P
       do llp=LPHIGH,LPLOW,-1
          ! loop range 0:LP in general, or LA:LA, or LB:LB, or not at all
          ! sum the smallest terms first, from HIGH LP to LOW LP
          DPRINT 'call shgi_gr_pseuL(',llp,',pps(llp),PIPSEU)'
          call shgi_gr_pseuL(llp,pps(llp),PIPSEU)
!!$          DCALL show('pp+lp',PIPSEU(1,:,:))
       enddo
    else
       DPRINT 'SHGI: only local contribution, LA=',LA,' LB=',LB,' LP=',LP
    endif

    ! second the local part, as maybe the largest contribution?:
    DPRINT 'call shgi_gr_pseu0(',ZCORE,'pps(LP+1),PIPSEU)'
    call shgi_gr_pseu0(ZCORE,xc,pps(LP+1),PIPSEU)
!!$    DCALL show('pp+loc',PIPSEU(1,:,:))
  end subroutine shgi_gr_pseu

 subroutine shgi_sd_pseu(ZCORE,LP,pps,xc,PIPSEU)
   use type_module, only: i8_kind
   use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
   USE_DEBUG, only: show
   implicit none
   real(RK)   , intent(in)    :: ZCORE
   integer(IK), intent(in)    :: LP
   type(ppt)  , intent(in)    :: pps(0:) ! (0:LP+1)
   real(RK)   , intent(in)    :: xc(3)
   real(RK)   , intent(inout) :: PIPSEU(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9,9)
   ! *** end of interface ***

   integer(IK) :: LPLOW,LPHIGH,llp
   integer (i8_kind) :: same

   ASSERT(size(pps)==2+LP)

   ! first the nonlocal part:
   same = whatis(IXEQC)
   select case(same)
   case (0) ! all centers differ
      LPLOW  = 0
      LPHIGH = LP ! in general all PPs contribute
   case (IAEQC) ! A=C
      ! A==C then only LA==LP may contibute:
      LPLOW  = LA
      LPHIGH = LA
   case (IBEQC) ! B=C
      ! B==C then only LB==LP may contibute:
      LPLOW  = LB
      LPHIGH = LB
   case (IXEQC) ! A=B=C
      ! same center, no forces at all
      WARN('A=B=C, why calling?')
      return !?
   case default
      ABORT('what?')
   end select
   DPRINT 'SHGI: only doing LPLOW=',LPLOW,' .. LPHIGH=',LPHIGH

   if( LPLOW <= LP )then
      call shgi_sd_pseu_set_c(xc,LP) ! calc D2P
      do llp=LPHIGH,LPLOW,-1
         ! loop range 0:LP in general, or LA:LA, or LB:LB, or not at all
         ! sum the smallest terms first, from HIGH LP to LOW LP
         call shgi_sd_pseuL(llp,pps(llp),PIPSEU) !!! (1)
      enddo
   else
      DPRINT 'SHGI: only local contribution, LA=',LA,' LB=',LB,' LP=',LP
   endif

   ! second the local part, as maybe the largest contribution?:
   call shgi_sd_pseu0(ZCORE,xc,pps(LP+1),PIPSEU)
 end subroutine shgi_sd_pseu

#if SHGI_USE_ANGULAR
  subroutine shgi_pseu0(Z,pp,nuc)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: W2, NORM, LAMBDA, YS
    use shgi_rad, only: doIL, doPQ1L
    use shgi_utils, only: doRadAng !!! (1)
    implicit none
    real(RK)   , intent(in)    :: Z ! Z - ZCORE
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: nuc(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,1+LA+LB)

    ! do the nuclear attraction of Z=Z-ZCORE,
    ! sets IL, ignores the state of IL:
    call doIL  (LA+LB,W2,LAMBDA, Z * NORM(:,3), IL)
    ! adds local PP contribution:
    call doPQ1L(LA+LB,W2,LAMBDA, pp, NORM(:,1), IL)
    call doRadAng(LA+LB,IL,YS(:,:,:,:,C00),nuc)
  end subroutine shgi_pseu0
#else
  subroutine shgi_pseu0(Z,x,pp,nuc)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: W, NORM, LAMBDA, S5
    use shgi_rad, only: doIL, doPQ1L
!   use shgi_utils, only: doD2F
    use shgi_utils, only: shgi_dww_to_dab
    use shgi_shr, only: SHR_D2Fv, SHR_PR2v
    implicit none
    real(RK)   , intent(in)    :: Z ! Z - ZCORE
    real(RK)   , intent(in)    :: x(3)
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: nuc(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK) :: IL(NAB,1+LA+LB)
    real(RK) :: F2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK) :: wc(NAB,3), w2(NAB)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    ! do the nuclear attraction of Z=Z-ZCORE,
    ! sets IL, ignores the state of IL:
    call doIL  (LA+LB  ,w2,LAMBDA, Z * NORM(:,3), IL)

    ! adds local PP contribution:
    call doPQ1L(LA+LB  ,w2,LAMBDA, pp, NORM(:,1), IL)

    ! FIXME: maybe to use YL, if available:
!   call doD2F(LA,LB,  IL,YL(:,:,:,C0,C0,C0,:),F2)

    ! D(lma)D(lmb) *  F:
    call SHR_D2Fv(NAB,LA,LB,wc,IL,F2)

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,1,F2)

    FPP_TIMER_START(tpr2)
    call SHR_PR2v(NAB,LA,LB,F2(:,:,:),S5(:,:,:,C0,C0,C0),nuc)
    FPP_TIMER_STOP(tpr2)
  end subroutine shgi_pseu0
#endif

#if SHGI_USE_ANGULAR
  subroutine shgi_gr_pseu0(Z,pp,nuc)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: W2, NORM, LAMBDA, GS
    use shgi_rad, only: doIL, doPQ1L
    use shgi_utils, only: doRadAng, shgi_gr_wd_to_abc
    implicit none
    real(RK)   , intent(in)    :: Z ! Z - ZCORE
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: nuc(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    ! *** end of interface ***

    real(RK) :: IL(NAB,1+LA+LB+1)
    real(RK) :: gr(NAB,2*LA+1,2*LB+1,6)

    ! do the nuclear attraction of Z=Z-ZCORE,
    ! sets IL, ignores the state of IL:
    call doIL  (LA+LB+1,W2,LAMBDA, Z * NORM(:,3), IL)
    ! adds local PP contribution:
    call doPQ1L(LA+LB+1,W2,LAMBDA, pp, NORM(:,1), IL)

    gr = 0.0_rk
    ! GW:
    call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GWX),gr(:,:,:,GWX))
    call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GWY),gr(:,:,:,GWY))
    call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GWZ),gr(:,:,:,GWZ))

    ! GD:
    call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GDX),gr(:,:,:,GDX))
    call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GDY),gr(:,:,:,GDY))
    call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GDZ),gr(:,:,:,GDZ))

    ! (GW,GD) -> (GA,GB,GC)
    call shgi_gr_wd_to_abc(gr,nuc,+1)
  end subroutine shgi_gr_pseu0
#else
  subroutine shgi_gr_pseu0(Z,x,pp,nuc)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: W, NORM, LAMBDA, S5
    use shgi_rad, only: doIL, doPQ1L
    use shgi_shr, only: SHR_D3Fv
    use shgi_utils, only: grFS, shgi_gr_wd_to_abc
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    real(RK)   , intent(in)    :: Z ! Z - ZCORE
    real(RK)   , intent(in)    :: x(3)
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: nuc(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    ! *** end of interface ***

    real(RK) :: IL(NAB,1+LA+LB+1)
    real(RK) :: gr(NAB,2*LA+1,2*LB+1,6)
    real(RK) :: F3(NAB,(LA+1)**2,(LB+1)**2,4) ! 4 = (P+1)**2
    real(RK) :: wc(NAB,3), w2(NAB)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    ! do the nuclear attraction of Z=Z-ZCORE,
    ! sets IL, ignores the state of IL:
    call doIL  (LA+LB+1,W2,LAMBDA, Z * NORM(:,3), IL)
    ! adds local PP contribution:
    call doPQ1L(LA+LB+1,W2,LAMBDA, pp, NORM(:,1), IL)

    ! FIXME: maybe to use YL, if available:
!   call doD3F(LA,LB,1,IL,YL(:,:,:,:,C0,C0,:),F3)

    ! D(lma)D(lmb) *  F:
    call SHR_D3Fv(NAB,LA,LB,1,wc,IL,F3)

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,4,F3) ! 4 = (P+1)**2

    ! F(g)*S(0), F(0)*S(g):
    gr = 0.0_rk
    call grFS(F3(:,:,:,:),S5(:,:,:,:,C0,C0),gr)

    ! (GW,GD) -> (GA,GB,GC)
    call shgi_gr_wd_to_abc(gr,nuc,+1)
  end subroutine shgi_gr_pseu0
#endif

  subroutine shgi_sd_pseu0(Z,x,pp,nuc)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: W, NORM, LAMBDA, S5
    use shgi_rad, only: doIL, doPQ1L
    use shgi_shr, only: SHR_D4Fv
    use shgi_utils, only: shgi_sd_wd_to_abc, sdFS, sdSYM
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    real(RK)   , intent(in)    :: Z ! Z - ZCORE
    real(RK)   , intent(in)    :: x(3)
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: nuc(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9,9)
    ! *** end of interface ***

    real(RK) :: IL(NAB,1+LA+LB+2)
    real(RK) :: sd(NAB,2*LA+1,2*LB+1,6,6)
    real(RK) :: F4(NAB,(LA+1)**2,(LB+1)**2,4,4) ! 4 = (P+1)**2
    real(RK) :: wc(NAB,3), w2(NAB)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    ! do the nuclear attraction of Z=Z-ZCORE,
    ! sets IL, ignores the state of IL:
    call doIL  (LA+LB+2,W2,LAMBDA, Z * NORM(:,3), IL)

    ! adds local PP contribution:
    call doPQ1L(LA+LB+2,W2,LAMBDA, pp, NORM(:,1), IL)

    ! FIXME: maybe to use YL, if available:
!   call doD4F(LA,LB,1,1,IL,YL(:,:,:,:,:,C0,:),F4)

    ! D(lma)D(lmb) *  F:
    call SHR_D4Fv(NAB,LA,LB,1,1,wc,IL,F4)

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,4*4,F4) ! 4 = (P+1)**2

    ! F(g)*S(0), F(0)*S(g):
    sd = 0.0_rk
    call sdFS(F4(:,:,:,:,:),S5(:,:,:,:,:,C0),sd)

    call sdSYM(sd)

    call shgi_sd_wd_to_abc(sd,nuc,+1)
  end subroutine shgi_sd_pseu0

  !
  !  non-local part (bessel_2):
  !

  subroutine shgi_pseuL(LP,pp,nuc)
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: &
         XD2,NORM,WDA,WDB,LAMBDA,CUTOFF
    use shgi_common, only: AEXP,BEXP
    use shgi_rad, only: doPQ2L
    implicit none
    integer(IK), intent(in)    :: LP
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: nuc(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,1+LA,1+LB)


    if(is_on(IAEQC))then
       ASSERT(LP==LA)
    endif
    if(is_on(IBEQC))then
       ASSERT(LP==LB)
    endif

    call doPQ2L(LP,LA,LB,XAC2,XBC2,XD2,AEXP,BEXP,pp,NORM(:,1),IL,WDA,WDB,LAMBDA,CUTOFF)
    call doRadAng2(LA,LB,IL,X2P(:,:,:,:,1+LP),nuc)
    DPRINT 'SHGI: shgi_pseuL(nuc,IL): LP=',LP,sum(nuc),sum(IL),sum(X2P),shape(X2P)

  end subroutine shgi_pseuL

  subroutine shgi_gr_pseuL(LP,pp,inuc)
    use type_module, only: i8_kind
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use shgi_common, only: &
         XD2,NORM,WDA,WDB,LAMBDA,CUTOFF
    use shgi_common, only: AEXP,BEXP
    use shgi_rad, only: doPQ2L
    implicit none
    integer(IK), intent(in)    :: LP
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(inout) :: inuc(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,1+LA+1,1+LB+1)
    integer (i8_kind) :: same
    real(RK)    :: nuc(NAB,2*LA+1,2*LB+1,9)

    DPRINT 'shgi_gr_pseuL: LP=',LP

    if(is_on(IAEQC))then
       ASSERT(LP==LA)
    endif
    if(is_on(IBEQC))then
       ASSERT(LP==LB)
    endif

    same = whatis(IXEQC)
    if( same == IXEQC )then
       WARN('A=B=C, why calling?')
       return
    endif

    nuc = ZERO

    DPRINT 'call doPQ2L(',LP,LA+1,LB+1,XAC2,XBC2,'...)'
    call doPQ2L(LP,LA+1,LB+1,XAC2,XBC2,XD2,AEXP,BEXP,pp,NORM(:,1),IL,WDA,WDB,LAMBDA,CUTOFF)

    select case(same)
    case (0)
       ! FIXME: only if moving A:
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAX),nuc(:,:,:,GAX))
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAY),nuc(:,:,:,GAY))
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAZ),nuc(:,:,:,GAZ))
       DPRINT  'shgi_gr_pseuL: doRadAng2(A)',sum(nuc(:,:,:,GAX)),sum(nuc(:,:,:,GAY)),sum(nuc(:,:,:,GAZ))
       ! FIXME: only if moving B:
       DPRINT  'shgi_gr_pseuL: doRadAng2(B)'
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBX),nuc(:,:,:,GBX))
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBY),nuc(:,:,:,GBY))
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBZ),nuc(:,:,:,GBZ))
       ! FIXME: only if moving C:
!!$       DPRINT  'shgi_gr_pseuL: doRadAng2(C)'
!!$       call doRadAng2(LA+1,LB+1,IL,G2P(:,:,:,:,1+LP,GCX),nuc(:,:,:,GCX))
!!$       call doRadAng2(LA+1,LB+1,IL,G2P(:,:,:,:,1+LP,GCY),nuc(:,:,:,GCY))
!!$       call doRadAng2(LA+1,LB+1,IL,G2P(:,:,:,:,1+LP,GCZ),nuc(:,:,:,GCZ))
       DPRINT  'shgi_gr_pseuL: doRadAng2(-AB)'
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAX),nuc(:,:,:,GCX),-ONE)
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAY),nuc(:,:,:,GCY),-ONE)
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAZ),nuc(:,:,:,GCZ),-ONE)
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBX),nuc(:,:,:,GCX),-ONE)
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBY),nuc(:,:,:,GCY),-ONE)
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBZ),nuc(:,:,:,GCZ),-ONE)
    case (IAEQC)
       ! FIXME: only if moving B:
       DPRINT  'shgi_gr_pseuL: doRadAng2(B,A=C)'
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBX),nuc(:,:,:,GBX))
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBY),nuc(:,:,:,GBY))
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBZ),nuc(:,:,:,GBZ))
       ! FIXME: only if moving C:
       DPRINT  'shgi_gr_pseuL: doRadAng2(A=C,B)'
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBX),nuc(:,:,:,GCX),-ONE)
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBY),nuc(:,:,:,GCY),-ONE)
       call doRadAng2(LA  ,LB+1,IL,G2P(:,:,:,:,1+LP,GBZ),nuc(:,:,:,GCZ),-ONE)
    case (IBEQC)
       ! FIXME: only if moving A:
       DPRINT  'shgi_gr_pseuL: doRadAng2(A,B=C)'
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAX),nuc(:,:,:,GAX))
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAY),nuc(:,:,:,GAY))
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAZ),nuc(:,:,:,GAZ))
       ! FIXME: only if moving C:
       DPRINT  'shgi_gr_pseuL: doRadAng2(B=C,A)'
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAX),nuc(:,:,:,GCX),-ONE)
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAY),nuc(:,:,:,GCY),-ONE)
       call doRadAng2(LA+1,LB  ,IL,G2P(:,:,:,:,1+LP,GAZ),nuc(:,:,:,GCZ),-ONE)
    end select
    inuc = inuc + nuc
  end subroutine shgi_gr_pseuL

   subroutine shgi_sd_pseuL(LP,pp,inuc)
   use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
   use shgi_rad, only: doPQ2L
   use shgi_common
   implicit none
   integer(IK), intent(in)    :: LP
   type(ppt)  , intent(in)    :: pp
   real(RK)   , intent(inout) :: inuc(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
   ! *** end of interface ***

   real(RK)    :: IL(NAB,3+LA,3+LB) !dervs
   real(RK)    :: nuc(NAB,2*LA+1,2*LB+1,9,9)

   DPRINT 'shgi_sd_pseuL: LP=',LP

   if(is_on(IAEQC))then
      ASSERT(LP==LA)
   endif
   if(is_on(IBEQC))then
      ASSERT(LP==LB)
   endif

   if (whatis (IXEQC) == IXEQC) then
      WARN('A=B=C, why calling?')
      return
   endif

   nuc = ZERO

   DPRINT 'call doPQ2L(',LP,LA+2,LB+2,XAC2,XBC2,'...)'
   call doPQ2L(LP,LA+2,LB+2,XAC2,XBC2,XD2,AEXP,BEXP,pp,NORM(:,1),IL,WDA,WDB,LAMBDA,CUTOFF)
   ! calculate IL of shape dependent on LA+2 and LB+2

      !    AxA block
      call doRadAng2(LA+2,LB  ,IL,D2P(:,:,:,:,1+LP,GAX,GAX),nuc(:,:,:,GAX,GAX))
      call doRadAng2(LA+2,LB  ,IL,D2P(:,:,:,:,1+LP,GAX,GAY),nuc(:,:,:,GAX,GAY))
      call doRadAng2(LA+2,LB  ,IL,D2P(:,:,:,:,1+LP,GAX,GAZ),nuc(:,:,:,GAX,GAZ))
      call doRadAng2(LA+2,LB  ,IL,D2P(:,:,:,:,1+LP,GAY,GAY),nuc(:,:,:,GAY,GAY))
      call doRadAng2(LA+2,LB  ,IL,D2P(:,:,:,:,1+LP,GAY,GAZ),nuc(:,:,:,GAY,GAZ))
      call doRadAng2(LA+2,LB  ,IL,D2P(:,:,:,:,1+LP,GAZ,GAZ),nuc(:,:,:,GAZ,GAZ))
      nuc(:,:,:,GAY,GAX)=nuc(:,:,:,GAX,GAY)
      nuc(:,:,:,GAZ,GAX)=nuc(:,:,:,GAX,GAZ)
      nuc(:,:,:,GAZ,GAY)=nuc(:,:,:,GAY,GAZ)

      ! BxB block
      call doRadAng2(LA  ,LB+2,IL,D2P(:,:,:,:,1+LP,GBX,GBX),nuc(:,:,:,GBX,GBX))
      call doRadAng2(LA  ,LB+2,IL,D2P(:,:,:,:,1+LP,GBX,GBY),nuc(:,:,:,GBX,GBY))
      call doRadAng2(LA  ,LB+2,IL,D2P(:,:,:,:,1+LP,GBX,GBZ),nuc(:,:,:,GBX,GBZ))
      call doRadAng2(LA  ,LB+2,IL,D2P(:,:,:,:,1+LP,GBY,GBY),nuc(:,:,:,GBY,GBY))
      call doRadAng2(LA  ,LB+2,IL,D2P(:,:,:,:,1+LP,GBY,GBZ),nuc(:,:,:,GBY,GBZ))
      call doRadAng2(LA  ,LB+2,IL,D2P(:,:,:,:,1+LP,GBZ,GBZ),nuc(:,:,:,GBZ,GBZ))
      nuc(:,:,:,GBY,GBX)=nuc(:,:,:,GBX,GBY)
      nuc(:,:,:,GBZ,GBX)=nuc(:,:,:,GBX,GBZ)
      nuc(:,:,:,GBZ,GBY)=nuc(:,:,:,GBY,GBZ)

      ! AxB block
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAX,GBX),nuc(:,:,:,GAX,GBX))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAX,GBY),nuc(:,:,:,GAX,GBY))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAX,GBZ),nuc(:,:,:,GAX,GBZ))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAY,GBX),nuc(:,:,:,GAY,GBX))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAY,GBY),nuc(:,:,:,GAY,GBY))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAY,GBZ),nuc(:,:,:,GAY,GBZ))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAZ,GBX),nuc(:,:,:,GAZ,GBX))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAZ,GBY),nuc(:,:,:,GAZ,GBY))
      call doRadAng2(LA+1  ,LB+1,IL,D2P(:,:,:,:,1+LP,GAZ,GBZ),nuc(:,:,:,GAZ,GBZ))

     !  BxA block
     nuc(:,:,:,GBX,GAX)=nuc(:,:,:,GAX,GBX)
     nuc(:,:,:,GBY,GAX)=nuc(:,:,:,GAX,GBY)
     nuc(:,:,:,GBZ,GAX)=nuc(:,:,:,GAX,GBZ)
     nuc(:,:,:,GBX,GAY)=nuc(:,:,:,GAY,GBX)
     nuc(:,:,:,GBY,GAY)=nuc(:,:,:,GAY,GBY)
     nuc(:,:,:,GBZ,GAY)=nuc(:,:,:,GAY,GBZ)
     nuc(:,:,:,GBX,GAZ)=nuc(:,:,:,GAZ,GBX)
     nuc(:,:,:,GBY,GAZ)=nuc(:,:,:,GAZ,GBY)
     nuc(:,:,:,GBZ,GAZ)=nuc(:,:,:,GAZ,GBZ)


      ! AxC block
      nuc(:,:,:,GAX,GCX:GCZ)=-nuc(:,:,:,GAX,GAX:GAZ)-nuc(:,:,:,GAX,GBX:GBZ) ! 3 terms
      nuc(:,:,:,GAY:GAZ,GCX)=-nuc(:,:,:,GAY:GAZ,GAX)-nuc(:,:,:,GAY:GAZ,GBX) ! 2 terms
      nuc(:,:,:,GAY,GCY:GCZ)=-nuc(:,:,:,GAY,GAY:GAZ)-nuc(:,:,:,GAY,GBY:GBZ) ! 2 terms
      nuc(:,:,:,GAZ,GCY:GCZ)=-nuc(:,:,:,GAZ,GAY:GAZ)-nuc(:,:,:,GAZ,GBY:GBZ) ! 2 terms

      ! CxA block
      nuc(:,:,:,GCX,GAX:GAZ)=nuc(:,:,:,GAX:GAZ,GCX) !3
      nuc(:,:,:,GCY:GCZ,GAX)=nuc(:,:,:,GAX,GCY:GCZ) !2
      nuc(:,:,:,GCY,GAY:GAZ)=nuc(:,:,:,GAY:GAZ,GCY) !2
      nuc(:,:,:,GCZ,GAY:GAZ)=nuc(:,:,:,GAY:GAZ,GCZ) !2

      ! BxC block bc= -bb - ba
      nuc(:,:,:,GBX,GCX:GCZ)=-nuc(:,:,:,GBX,GBX:GBZ)-nuc(:,:,:,GBX,GAX:GAZ) ! 3 terms
      nuc(:,:,:,GBY:GBZ,GCX)=-nuc(:,:,:,GBY:GBZ,GBX)-nuc(:,:,:,GBY:GBZ,GAX) ! 2 terms
      nuc(:,:,:,GBY,GCY:GCZ)=-nuc(:,:,:,GBY,GBY:GBZ)-nuc(:,:,:,GBY,GAY:GAZ) ! 2 terms
      nuc(:,:,:,GBZ,GCY:GCZ)=-nuc(:,:,:,GBZ,GBY:GBZ)-nuc(:,:,:,GBZ,GAY:GAZ) ! 2 terms
      ! CxB block
      nuc(:,:,:,GCX,GBX:GBZ)=nuc(:,:,:,GBX:GBZ,GCX) !3
      nuc(:,:,:,GCY:GCZ,GBX)=nuc(:,:,:,GBX,GCY:GCZ) !2
      nuc(:,:,:,GCY,GBY:GBZ)=nuc(:,:,:,GBY:GBZ,GCY) !2
      nuc(:,:,:,GCZ,GBY:GBZ)=nuc(:,:,:,GBY:GBZ,GCZ) !2

      ! CxC block
      nuc(:,:,:,GCX,GCX:GCZ)=nuc(:,:,:,GAX,GAX:GAZ)+nuc(:,:,:,GBX,GBX:GBZ)+nuc(:,:,:,GAX,GBX:GBZ)+nuc(:,:,:,GBX,GAX:GAZ)
      nuc(:,:,:,GCY,GCX)=nuc(:,:,:,GCX,GCY)
      nuc(:,:,:,GCZ,GCX)=nuc(:,:,:,GCX,GCZ)
      nuc(:,:,:,GCY,GCY:GCZ)=nuc(:,:,:,GAY,GAY:GAZ)+nuc(:,:,:,GBY,GBY:GBZ)+nuc(:,:,:,GAY,GBY:GBZ)+nuc(:,:,:,GBY,GAY:GAZ)
      nuc(:,:,:,GCZ,GCY)=nuc(:,:,:,GCY,GCZ)
      nuc(:,:,:,GCZ,GCZ:GCZ)=nuc(:,:,:,GAZ,GAZ:GAZ)+nuc(:,:,:,GBZ,GBZ:GBZ)+nuc(:,:,:,GAZ,GBZ:GBZ)+nuc(:,:,:,GBZ,GAZ:GAZ)

    inuc = inuc + nuc
  end subroutine shgi_sd_pseuL


  subroutine doRadAng2(LLA,LLB,IL,Y,NUC,fact)
    ! couples two radial derivatives w.r.t. A and B
    ! with corresponding "angular" part.
    ! Adds integral to NUC
    ! DONT FORGET TO CLEAR NUC!
    implicit none
    integer(IK), intent(in) :: LLA,LLB
    real(RK), intent(in)    :: IL(:,:,:)  ! (NAB,              1+LLA,1+LLB)
    real(RK), intent(in)    :: Y(:,:,:,:) ! (    2*LA+1,2*LB+1,1+LLA,1+LLB)
    real(RK), intent(inout) :: NUC(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK), intent(in)    :: fact
    optional :: fact
    ! *** end of interface ***

    integer(IK) :: ma,mb
    integer(IK) :: ila,ilb

    FPP_TIMER_START(tra2)
    if(.not.present(fact))then
    do ilb=1,1+LLB    ! Energy: LLB=LB, Gradients: LLB=LB+1
       do ila=1,1+LLA ! Energy: LLA=LA, Gradients: LLA=LA+1
          do mb=1,2*LB+1
             do ma=1,2*LA+1
                NUC(:,ma,mb) =  NUC(:,ma,mb) &
                     + IL(:,ila,ilb) * Y(ma,mb,ila,ilb)
             enddo
          enddo
       enddo
    enddo
    else
    do ilb=1,1+LLB    ! Energy: LLB=LB, Gradients: LLB=LB+1
       do ila=1,1+LLA ! Energy: LLA=LA, Gradients: LLA=LA+1
          do mb=1,2*LB+1
             do ma=1,2*LA+1
                NUC(:,ma,mb) =  NUC(:,ma,mb) &
                     + IL(:,ila,ilb) * Y(ma,mb,ila,ilb) * fact
             enddo
          enddo
       enddo
    enddo
    endif
    FPP_TIMER_STOP(tra2)
  end subroutine doRadAng2

  !--------------- End of module -------------------------------------
end module shgi_pseudo


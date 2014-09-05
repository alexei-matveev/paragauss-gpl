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
module shgi_pcm
  !-------------------------------------------------------------------
  !
  ! Here _pcm  == Point Charge  Module. Computes the integrals  of the
  ! field of point charge collection.
  !
  ! Copyright (c) 2006-2009 Alexey Shor
  ! Copyright (c) 2006-2013 Alexei Matveev
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
       IK=>i4_kind, RK=>r8_kind, & ! type specification parameters
       I8K=>i8_kind
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
  public :: shgi_pc
  public :: shgi_gr_pc
  public :: shgi_sd_pc
  public :: shgi_pcs
  public :: shgi_gr_pcs
  public :: shgi_pc_grad
  public :: shgi_sd_pcs
#ifdef WITH_EFP
  public :: shgi_PCs_wrap
#endif

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains
#ifdef WITH_EFP
  subroutine shgi_PCs_wrap(IU1,IE1,IL1,IU2,IE2,IL2, uas,NUCL)
    use unique_atom_module   , only: uat=>unique_atom_type
    use pointcharge_module, only: pcs=>pointcharge_array
    use options_module, only: options_integral_expmax
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: upack2c, shgi_timing
    integer(IK) , intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)   , intent(in)  :: uas(:)   ! array of unique atoms, normally all of them
    real(RK)    , intent(out) :: NUCL(:,:,:,:)     ! (N2E,N1E,2*L2+1,2*L1+1  )
    ! *** end of interface ***

    real(RK), allocatable :: PINUCL(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    integer(IK) :: memstat


    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(0,0,0)

    ! compute overlap S5:
    call shgi_set_ovrl(LA,LB,0,0,0) !(2) shgi_pot_drv

    allocate(PINUCL(NAB,2*LA+1,2*LB+1),stat=memstat)
    ASSERT(memstat==0)
    PINUCL=0.0_RK

    call shgi_pcs(pcs,PINUCL)

    call upack2c(PINUCL,NUCL,-1)

    deallocate(PINUCL,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()

    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_PCs_wrap
#endif

  subroutine shgi_pcs(pcs,nu)
    use datatype, only: pct=>pointcharge_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_utils, only: shgi_dww_to_dab
    use shgi_shr, only: SHR_PR2v
    implicit none
    type(pct)  , intent(in)    :: pcs(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: f2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: p2(NAB,(LA+1)**2,(LB+1)**2)
#ifdef WITH_EFP
    real(RK)    :: LAM1(NAB),LAM2(NAB),ETA(NAB),NRM(NAB),a,c
#endif

    DPRINT  'SHGI: shgi_pcs(',size(pcs),'nu)'

    FPP_TIMER_START(pcs)

    p2 = 0.0_rk

    do ua=1,size(pcs)
      z = pcs(ua)%Z
      do ea=1,pcs(ua)%n_equal_charges

        x = pcs(ua)%position(:,ea)

        wc(:,1) = W(:,1) - x(1)
        wc(:,2) = W(:,2) - x(2)
        wc(:,3) = W(:,3) - x(3)

        ! computes radial dervs -> haramonic dervs:
        call doD2Fw(LA,LB,wc,LAMBDA,z*NORM(:,3),f2)

        p2 = p2 + f2

#ifdef WITH_EFP
        if(pcs(ua)%c /= 0.0_RK) then
           a=pcs(ua)%A
           c=pcs(ua)%C
           ETA=LAMBDA/(LAMBDA+a)
           NRM=two*sqrt(LAMBDA/pi)*ETA*c*NORM(:,2)
           LAM1=ETA*a
           LAM2=ETA*LAMBDA

           call doD2Fw_s(LA,LB,wc,LAM1,LAM2,z*NRM,f2)

           p2 = p2 - f2
        end if
#endif
      enddo
    enddo

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,1,p2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,p2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)

    FPP_TIMER_STOP(pcs)
  end subroutine shgi_pcs

  subroutine shgi_gr_pcs(pcs,nu)
    ! WARNING: returns 6 grads w.r.t. W and D
    use datatype, only: pct=>pointcharge_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_utils, only: grFS
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(pct)  , intent(in)    :: pcs(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    integer(IK), parameter :: LC=1
    real(RK)    :: wc(NAB,3)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
#ifdef WITH_EFP
    real(RK)    :: LAM1(NAB),LAM2(NAB),ETA(NAB),NRM(NAB),a,c
#endif

    DPRINT  'SHGI: shgi_gr_pcs(',size(pcs),'nu)'

    FPP_TIMER_START(pcs)

    p3 = 0.0_rk

    do ua=1,size(pcs)
      z = pcs(ua)%Z
      do ea=1,pcs(ua)%n_equal_charges

        x = pcs(ua)%position(:,ea)

        wc(:,1) = W(:,1) - x(1)
        wc(:,2) = W(:,2) - x(2)
        wc(:,3) = W(:,3) - x(3)

        ! computes radial dervs -> haramonic dervs:
        call doD3Fw(LA,LB,LC,wc,LAMBDA,z*NORM(:,3),f3)

        p3 = p3 + f3

#ifdef WITH_EFP
        if(pcs(ua)%c /= 0.0_RK) then
           a=pcs(ua)%A
           c=pcs(ua)%C
           ETA=LAMBDA/(LAMBDA+a)
           NRM=two*sqrt(LAMBDA/pi)*ETA*c*NORM(:,2)
           LAM1=ETA*a
           LAM2=ETA*LAMBDA

           call doD3Fw_s(LA,LB,LC,wc,LAM1,LAM2,z*NRM,f3)

           p3 = p3 - f3
        end if
#endif
      enddo
    enddo

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

    FPP_TIMER_STOP(pcs)
  end subroutine shgi_gr_pcs

  subroutine shgi_sd_pcs(pcs,nu)
    ! WARNING: returns 6 grads w.r.t. W and D
    use datatype, only: pct=>pointcharge_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_utils, only: sdFS, sdSYM
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(pct)  , intent(in)    :: pcs(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    integer(IK), parameter :: LC=1,LD=1
    real(RK)    :: wc(NAB,3)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)    :: p4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)

    DPRINT  'SHGI: shgi_sd_pcs(',size(pcs),'nu)'

    FPP_TIMER_START(pcs)

    p4 = 0.0_rk

    do ua=1,size(pcs)
      z = pcs(ua)%Z
      do ea=1,pcs(ua)%n_equal_charges

        x = pcs(ua)%position(:,ea)

        wc(:,1) = W(:,1) - x(1)
        wc(:,2) = W(:,2) - x(2)
        wc(:,3) = W(:,3) - x(3)

        ! computes radial dervs -> haramonic dervs:
        call doD4Fw(LA,LB,LC,LD,wc,LAMBDA,z*NORM(:,3),f4)

        p4 = p4 + f4
      enddo
    enddo

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2*(LD+1)**2,p4)

    ! F(i,j)*S(0), F(i,0)*S(0,j), F(0)*S(i,j):
    call sdFS(p4,S5(:,:,:,:,:,C0),nu)

    call sdSYM(nu)

    FPP_TIMER_STOP(pcs)
  end subroutine shgi_sd_pcs

  subroutine shgi_pc(z,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_utils, only: shgi_dww_to_dab
    use shgi_shr, only: SHR_PR2v
    implicit none
    real(RK)   , intent(in)    :: z, x(3)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: f2(NAB,(LA+1)**2,(LB+1)**2)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)

    ! computes radial dervs -> haramonic dervs:
    call doD2Fw(LA,LB,wc,LAMBDA,z*NORM(:,3),f2)

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,1,f2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,f2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)
  end subroutine shgi_pc

  subroutine shgi_gr_pc(z,x,nu,c,a)
    ! WARNING: returns 6 grads w.r.t. W and D
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_utils, only: grFS
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    real(RK)   , intent(in)    :: z, x(3)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in), optional :: c,a
    ! *** end of interface ***

    integer(IK), parameter :: LC=1
    real(RK)    :: wc(NAB,3)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
#ifdef WITH_EFP
    real(RK)    :: ff3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: LAM1(NAB), LAM2(NAB), ETA(NAB), NRM(NAB)
#endif
    real(RK)    :: c1,a1

    c1=0.0_RK; a1=0.0_RK
    if(present(c) .and. present(a)) then
       c1=c; a1=a
    end if

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)

    ! computes radial dervs -> haramonic dervs:
    call doD3Fw(LA,LB,LC,wc,LAMBDA,z*NORM(:,3),f3)

#ifdef WITH_EFP
    if(c1 /= 0.0_RK) then
       ETA=LAMBDA/(LAMBDA+a1)
       NRM=two*sqrt(LAMBDA/pi)*ETA*c1*NORM(:,2)
       LAM1=ETA*a1
       LAM2=ETA*LAMBDA

       call doD3Fw_s(LA,LB,LC,wc,LAM1,LAM2,z*NRM,ff3)

       f3 = f3 - ff3
    end if
#endif

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,f3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(f3,S5(:,:,:,:,C0,C0),nu)
  end subroutine shgi_gr_pc

  subroutine shgi_pc_grad(IU1,IE1,IL1,IU2,IE2,IL2,uas,pcs,PCGR)
    use datatype          , only: arrmat4
    use unique_atom_module, only: uat=>unique_atom_type
    use datatype          , only: pct=>pointcharge_type
    use options_module, only: options_integral_expmax
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: shgi_X_wd_store, shgi_timing
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)       :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)       :: uas(:)  ! array of unique atoms, normally all of them
    type(pct)  , intent(in)       :: pcs(:)  ! array of point charges
    type(arrmat4) , intent(inout) :: PCGR(:) ! (num of tot-sym grads)
    ! *** end of interface ***
    !------------ Declaration of local variables --------------------
    integer(IK) :: up,ep,n_equal
    real(RK)    :: x(3)
    real(RK)    :: z,c,a
    ! temp storage for pc grads
    real(RK), allocatable :: PIPCGR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    integer(IK) :: memstat
    !------------ Executable code -----------------------------------

    DPRINT 'SHGI: shgi_pc_grad: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)
    FPP_TIMER_START(pcs)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(1,0,0)

    ! compute overlap S5:
    call shgi_set_ovrl(LA,LB,1,0,0)

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    allocate(PIPCGR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    do up=1,size(pcs)
       n_equal=pcs(up)%n_equal_charges
       z=pcs(up)%Z
       do ep=1,n_equal
          PIPCGR = 0.0_rk
          x=pcs(up)%position(:,ep)
          c=pcs(up)%C
          a=pcs(up)%A
          call shgi_gr_pc(z,x,PIPCGR,c,a)
          call shgi_X_wd_store(PC,up,ep,PIPCGR,PCGR,-1)
       end do
    end do

    deallocate(PIPCGR,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()

    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_pc_grad

  subroutine shgi_sd_pc(z,x,nu)
    ! WARNING: returns 6 grads w.r.t. W and D
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_utils, only: sdFS, sdSYM
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    real(RK)   , intent(in)    :: z, x(3)
    real(RK)   , intent(inout) :: nu(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    ! *** end of interface ***

    integer(IK), parameter :: LC=1,LD=1
    real(RK)    :: wc(NAB,3)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)

    ! computes radial dervs -> haramonic dervs:
    call doD4Fw(LA,LB,LC,LD,wc,LAMBDA,z*NORM(:,3),f4)

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2*(LD+1)**2,f4)

    ! F(i,j)*S(0), F(i,0)*S(0,j), F(0)*S(i,j):
    call sdFS(f4,S5(:,:,:,:,:,C0),nu)

    call sdSYM(nu)
  end subroutine shgi_sd_pc

  subroutine doD2Fw(LA,LB,w,lambda,norm,DF)
    !
    !   D2F(w) = D(lma)D(lmb) F(w2/2)
    !
    !         ( D(lm) == Dw(lm) )
    !
#ifdef FPP_TIMERS
    use shgi_cntrl, only: FPP_TIMER_VARS(tpcr),FPP_TIMER_VARS(td3f)
#endif
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_D2Fv
    implicit none
    integer(IK), intent(in)  :: LA,LB
    real(RK)   , intent(in)  :: w(:,:)        ! (NAB,3)
    real(RK)   , intent(in)  :: lambda(:)     ! (NAB)
    real(RK)   , intent(in)  :: norm(:)       ! (NAB)
    real(RK)   , intent(out) :: DF(:,:,:)     ! (NAB,(LA+1)**2,(LB+1)**2)
    ! *** end of interface ***

    real(RK)    :: w2(size(w,1))
    real(RK)    :: fL(size(w,1),1+LA+LB)

    w2(:)   = w(:,1)**2 + w(:,2)**2 + w(:,3)**2

    FPP_TIMER_START(tpcr)
    ! radial part:
    ASSERT(LC==0)
    call doIL(LA+LB,w2,lambda,norm,fL)
    FPP_TIMER_STOP(tpcr)

    FPP_TIMER_START(td3f)
    ! Double derivatives D(lma)D(lmb) * F:
    call SHR_D2Fv(size(w,1),LA,LB,w,fL,DF)
    FPP_TIMER_STOP(td3f)
  end subroutine doD2Fw

  subroutine doD2Fw_s(LA,LB,w,lambda1,lambda2,norm,DF)
    !
    !   D2F(w) = D(lma)D(lmb) F(w2/2)
    !
    !         ( D(lm) == Dw(lm) )
    !
#ifdef FPP_TIMERS
    use shgi_cntrl, only: FPP_TIMER_VARS(tpcr),FPP_TIMER_VARS(td3f)
#endif
    use shgi_rad, only: doEIL
    use shgi_shr, only: SHR_D2Fv
    implicit none
    integer(IK), intent(in)  :: LA,LB
    real(RK)   , intent(in)  :: w(:,:)        ! (NAB,3)
    real(RK)   , intent(in)  :: lambda1(:)    ! (NAB)
    real(RK)   , intent(in)  :: lambda2(:)    ! (NAB)
    real(RK)   , intent(in)  :: norm(:)       ! (NAB)
    real(RK)   , intent(out) :: DF(:,:,:)     ! (NAB,(LA+1)**2,(LB+1)**2)
    ! *** end of interface ***

    real(RK)    :: w2(size(w,1))
    real(RK)    :: fL(size(w,1),1+LA+LB)

    w2(:)   = w(:,1)**2 + w(:,2)**2 + w(:,3)**2

    FPP_TIMER_START(tpcr)
    ! radial part:
    call doEIL(LA+LB,w2,lambda1,lambda2,norm,fL)
    FPP_TIMER_STOP(tpcr)

    FPP_TIMER_START(td3f)
    ! Double derivatives D(lma)D(lmb) * F:
    call SHR_D2Fv(size(w,1),LA,LB,w,fL,DF)
    FPP_TIMER_STOP(td3f)
  end subroutine doD2Fw_s

  subroutine doD3Fw_s(LA,LB,LC,w,lambda1,lambda2,norm,DF)
    !
    !   D3F(w) = D(lma)D(lmb)D(lmc) F(w2/2)
    !
    !         ( D(lm) == Dw(lm) )
    !
#ifdef FPP_TIMERS
    use shgi_cntrl, only: FPP_TIMER_VARS(tpcr),FPP_TIMER_VARS(td3f)
#endif
    use shgi_rad, only: doEIL
    use shgi_shr, only: SHR_D3Fv
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC
    real(RK)   , intent(in)  :: w(:,:)        ! (NAB,3)
    real(RK)   , intent(in)  :: lambda1(:)    ! (NAB)
    real(RK)   , intent(in)  :: lambda2(:)    ! (NAB)
    real(RK)   , intent(in)  :: norm(:)       ! (NAB)
    real(RK)   , intent(out) :: DF(:,:,:,:)   ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    ! *** end of interface ***

    real(RK)    :: w2(size(w,1))
    real(RK)    :: fL(size(w,1),1+LA+LB+LC)

    w2(:)   = w(:,1)**2 + w(:,2)**2 + w(:,3)**2

    FPP_TIMER_START(tpcr)
    ! radial part:
    call doEIL(LA+LB+LC,w2,lambda1,lambda2,norm,fL)
    FPP_TIMER_STOP(tpcr)

    FPP_TIMER_START(td3f)
    ! Double derivatives D(lma)D(lmb) * F:
    call SHR_D3Fv(size(w,1),LA,LB,LC,w,fL,DF)
    FPP_TIMER_STOP(td3f)
  end subroutine doD3Fw_s

  subroutine doD3Fw(LA,LB,LC,w,lambda,norm,DF)
    !
    !   D3F(w) = D(lma)D(lmb)D(lmc) F(w2/2)
    !
    !         ( D(lm) == Dw(lm) )
    !
#ifdef FPP_TIMERS
    use shgi_cntrl, only: FPP_TIMER_VARS(tpcr),FPP_TIMER_VARS(td3f)
#endif
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_D3Fv
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC
    real(RK)   , intent(in)  :: w(:,:)        ! (NAB,3)
    real(RK)   , intent(in)  :: lambda(:)     ! (NAB)
    real(RK)   , intent(in)  :: norm(:)       ! (NAB)
    real(RK)   , intent(out) :: DF(:,:,:,:)   ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    ! *** end of interface ***

    real(RK)    :: w2(size(w,1))
    real(RK)    :: fL(size(w,1),1+LA+LB+LC)

    w2(:)   = w(:,1)**2 + w(:,2)**2 + w(:,3)**2

    FPP_TIMER_START(tpcr)
    ! radial part:
    call doIL(LA+LB+LC,w2,lambda,norm,fL)
    FPP_TIMER_STOP(tpcr)

    FPP_TIMER_START(td3f)
    ! Double derivatives D(lma)D(lmb) * F:
    call SHR_D3Fv(size(w,1),LA,LB,LC,w,fL,DF)
    FPP_TIMER_STOP(td3f)
  end subroutine doD3Fw

  subroutine doD4Fw(LA,LB,LC,LD,w,lambda,norm,DF)
    !
    !   D4F(w) = D(lma)D(lmb)D(lmc)D(lmd) F(w2/2)
    !
    !         ( D(lm) == Dw(lm) )
    !
#ifdef FPP_TIMERS
    use shgi_cntrl, only: FPP_TIMER_VARS(tpcr),FPP_TIMER_VARS(td3f)
#endif
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_D4Fv
    implicit none
    integer(IK), intent(in)  :: LA,LB,LC,LD
    real(RK)   , intent(in)  :: w(:,:)        ! (NAB,3)
    real(RK)   , intent(in)  :: lambda(:)     ! (NAB)
    real(RK)   , intent(in)  :: norm(:)       ! (NAB)
    real(RK)   , intent(out) :: DF(:,:,:,:,:) ! (NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    ! *** end of interface ***

    real(RK)    :: w2(size(w,1))
    real(RK)    :: fL(size(w,1),1+LA+LB+LC+LD)

    w2(:)   = w(:,1)**2 + w(:,2)**2 + w(:,3)**2

    FPP_TIMER_START(tpcr)
    ! radial part:
    call doIL(LA+LB+LC+LD,w2,lambda,norm,fL)
    FPP_TIMER_STOP(tpcr)

    FPP_TIMER_START(td3f)
    ! Double derivatives D(lma)D(lmb) * F:
    call SHR_D4Fv(size(w,1),LA,LB,LC,LD,w,fL,DF)
    FPP_TIMER_STOP(td3f)
  end subroutine doD4Fw

  !--------------- End of module -------------------------------------
end module shgi_pcm

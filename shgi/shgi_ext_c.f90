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
module shgi_ext_c
  !---------------------------------------------------------------
  !
  ! Here _ext_c == external centers module (dipoles, octupoles,
  ! quadrupoles, induced dipoles)
  !
  ! Copyright (c) 2007-2011 Alexey Shor
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
  public :: shgi_ext,shgi_ext_gr,shgi_X_grad,shgi_X_torq
#ifdef WITH_EFP
  public :: shgi_X_wrap
#endif

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
#ifdef WITH_EFP
  subroutine shgi_X_wrap(IU1,IE1,IL1,IU2,IE2,IL2, uas,NUCL)
    use unique_atom_module   , only: uat=>unique_atom_type
    use options_module, only: options_integral_expmax
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: upack2c, shgi_timing
    integer(IK) , intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)   , intent(in)  :: uas(:)   ! array of unique atoms, normally all of them
    real(RK)    , intent(out) :: NUCL(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
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

    call shgi_ext(PINUCL)

    call upack2c(PINUCL,NUCL,-1)

    deallocate(PINUCL,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()

    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_X_wrap
#endif
!*************************************************************
  subroutine shgi_ext(nu)
    use shgi_ang, only: shgi_set_lcde
    use point_dqo_module, only: pd_array,N_pd
    use point_dqo_module, only: pq_array,N_pq
    use point_dqo_module, only: po_array,N_po
    use point_dqo_module, only: rc_array,N_rc
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    if(N_pd > 0) then
       call shgi_set_lcde(1,0,0)
       call shgi_pd(pd_array,nu)
    end if

    if(N_pq > 0) then
       call shgi_set_lcde(1,1,0)
       call shgi_pq(pq_array,nu)
    end if

    if(N_po > 0) then
       call shgi_set_lcde(1,1,1)
       call shgi_po(po_array,nu)
    end if

    if(N_rc > 0) then
       call shgi_set_lcde(0,0,0)
       call shgi_rc(rc_array,nu)
    end if

  end subroutine shgi_ext
!*************************************************************
  subroutine shgi_ext_gr(nu)
    use shgi_ang, only: shgi_set_lcde
    use point_dqo_module, only: pd_array,N_pd
    use point_dqo_module, only: pq_array,N_pq
    use point_dqo_module, only: po_array,N_po
    use point_dqo_module, only: rc_array,N_rc
    use induced_dipoles_module, only: ipd_array,N_ipd
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    if(N_pd > 0) then
       call shgi_set_lcde(1,1,0)
       call shgi_pd_gr(pd_array,nu)
    end if

    if(N_pq > 0) then
       call shgi_set_lcde(1,1,1)
       call shgi_pq_gr(pq_array,nu)
    end if

    if(N_po > 0) then
       call shgi_set_lcde(1,1,1,1)
       call shgi_po_gr(po_array,nu)
    end if

    if(N_ipd > 0) then
       call shgi_set_lcde(1,1,0)
       call shgi_ipd_gr(ipd_array,nu)
    end if

    if(N_rc > 0) then
       call shgi_set_lcde(1,0,0)
       call shgi_rc_gr(rc_array,nu)
    end if

  end subroutine shgi_ext_gr
!*************************************************************
  subroutine shgi_X_grad(IU1,IE1,IL1,IU2,IE2,IL2, &
       uas,   &
       pds,   &
       pqs,   &
       pos,   &
       rcs,   &
       ipds,  &
       PDGR,  &
       PQGR,  &
       POGR,  &
       PRCGR, &
       PIDGR)
    use datatype             , only: arrmat4
    use options_module, only: options_integral_expmax
    use unique_atom_module   , only: uat=>unique_atom_type
    use point_dqo_module     , only: pdt=>pointdipole_type,N_pd
    use point_dqo_module     , only: pqt=>pointquadrupole_type,N_pq
    use point_dqo_module     , only: pot=>pointoctopole_type,N_po
    use point_dqo_module     , only: rct=>repulsion_type,N_rc
    use induced_dipoles_module, only: ipdt=>ind_dipole_type,N_ipd
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: shgi_X_wd_store, shgi_timing
    !------------ Declaration of formal parameters ------------------
    integer(IK) , intent(in)       :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)   , intent(in)       :: uas(:)   ! array of unique atoms, normally all of them
    type(pdt)   , intent(in)       :: pds(:)   ! array of point dipoles
    type(pqt)   , intent(in)       :: pqs(:)   ! array of point quadrupoles
    type(pot)   , intent(in)       :: pos(:)   ! array of point octopoles
    type(rct)   , intent(in)       :: rcs(:)   ! array of repulsive centers
    type(ipdt)  , intent(in)       :: ipds(:)  ! array of point dipoles
    type(arrmat4)  , intent(inout) :: PDGR(:)  ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: PQGR(:)  ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: POGR(:)  ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: PRCGR(:) ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: PIDGR(:) ! (num of tot-sym grads)
    ! *** end of interface ***
    !------------ Declaration of local variables --------------------
    integer(IK) :: up,ep,n_equal
    real(RK)    :: x(3),c,a
    real(RK)    :: d(3),q(3,3),o(3,3,3)
    ! temp storage for pc grads
    real(RK), allocatable :: PIPDGR(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIPQGR(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIPOGR(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIRCGR(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIPIDGR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    integer(IK) :: memstat
    real(RK),parameter    :: small=1.0e-8_RK
    !------------ Executable code -----------------------------------

    DPRINT 'SHGI: shgi_X_grad: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

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

    call shgi_set_lcde(1,0,0)

    call shgi_set_ovrl(LA,LB,1,0,0)

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    if(N_pd > 0) then
       FPP_TIMER_START(pd)

       call shgi_set_lcde(1,1,0)
       allocate(PIPDGR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_pd
          n_equal=pds(up)%n_equal_dipoles
          do ep=1,n_equal
             PIPDGR = 0.0_rk
             d=pds(up)%dipole(:,ep)
             if(sum(d*d) <= small) cycle
             x=pds(up)%position(:,ep)
             call shgi_gr_pd(d,x,PIPDGR)
             call shgi_X_wd_store(PD,up,ep,PIPDGR,PDGR,-1)
          end do
       end do

       deallocate(PIPDGR,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pd)
    end if

    if(N_pq > 0) then
       FPP_TIMER_START(pq)

       call shgi_set_lcde(1,1,1)
       allocate(PIPQGR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_pq
          n_equal=pqs(up)%n_equal_qpoles
          do ep=1,n_equal
             PIPQGR = 0.0_rk
             q=pqs(up)%quadrupole(:,:,ep)
             if(sum(q*q) <= small) cycle
             x=pqs(up)%position(:,ep)
             call shgi_gr_pq(q,x,PIPQGR)
             call shgi_X_wd_store(PQ,up,ep,PIPQGR,PQGR,-1)
          end do
       end do

       deallocate(PIPQGR,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pq)
    end if

    if(N_po > 0) then
       FPP_TIMER_START(po)

       call shgi_set_lcde(1,1,1,1)
       allocate(PIPOGR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_po
          n_equal=pos(up)%n_equal_opoles
          do ep=1,n_equal
             PIPOGR = 0.0_rk
             o=pos(up)%octopole(:,:,:,ep)
             if(sum(o*o) <= small) cycle
             x=pos(up)%position(:,ep)
             call shgi_gr_po(o,x,PIPOGR)
             call shgi_X_wd_store(PO,up,ep,PIPOGR,POGR,-1)
          end do
       end do

       deallocate(PIPOGR,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(po)
    end if

    if(N_ipd > 0) then
       FPP_TIMER_START(pd)

       call shgi_set_lcde(1,1,0)
       allocate(PIPIDGR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_ipd
          n_equal=ipds(up)%n_equal_dipoles
          do ep=1,n_equal
             PIPIDGR = 0.0_rk
             d=-(ipds(up)%idipole(:,ep)+ipds(up)%idipole1(:,ep))*0.5_RK
             x=ipds(up)%position(:,ep)
             call shgi_gr_pd(d,x,PIPIDGR)
             call shgi_X_wd_store(IPD,up,ep,PIPIDGR,PIDGR,-1)
          end do
       end do

       deallocate(PIPIDGR,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pd)
    end if

    if(N_rc > 0) then
       FPP_TIMER_START(pd)

       call shgi_set_lcde(1,0,0)
       allocate(PIRCGR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_rc
          n_equal=rcs(up)%n_equal_centers
          do ep=1,n_equal
             PIRCGR = 0.0_rk
             c=rcs(up)%C
             a=rcs(up)%A
             x=rcs(up)%position(:,ep)
             call shgi_gr_rc(c,a,x,PIRCGR)
             call shgi_X_wd_store(RC,up,ep,PIRCGR,PRCGR,-1)
          end do
       end do

       deallocate(PIRCGR,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pd)
    end if

    call shgi_close_ab()

    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_X_grad
!*************************************************************
  subroutine shgi_X_torq(IU1,IE1,IL1,IU2,IE2,IL2, &
       uas,    &
       pds,    &
       pqs,    &
       pos,    &
       ipds,   &
       PDTRQ,  &
       PQTRQ,  &
       POTRQ,  &
       PIDTRQ)
    use datatype             , only: arrmat4
    use options_module, only: options_integral_expmax
    use unique_atom_module   , only: uat=>unique_atom_type
    use point_dqo_module     , only: pdt=>pointdipole_type,N_pd
    use point_dqo_module     , only: pqt=>pointquadrupole_type,N_pq
    use point_dqo_module     , only: pot=>pointoctopole_type,N_po
    use induced_dipoles_module, only: ipdt=>ind_dipole_type,N_ipd
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: shgi_X_wd_store, shgi_timing
    !------------ Declaration of formal parameters ------------------
    integer(IK) , intent(in)       :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)   , intent(in)       :: uas(:)    ! array of unique atoms, normally all of them
    type(pdt)   , intent(in)       :: pds(:)    ! array of point dipoles
    type(pqt)   , intent(in)       :: pqs(:)    ! array of point quadrupoles
    type(pot)   , intent(in)       :: pos(:)    ! array of point octopoles
    type(ipdt)  , intent(in)       :: ipds(:)   ! array of point dipoles
    type(arrmat4)  , intent(inout) :: PDTRQ(:)  ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: PQTRQ(:)  ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: POTRQ(:)  ! (num of tot-sym grads)
    type(arrmat4)  , intent(inout) :: PIDTRQ(:) ! (num of tot-sym grads)
    ! *** end of interface ***
    !------------ Declaration of local variables --------------------
    integer(IK) :: up,ep,n_equal
    real(RK)    :: x(3)!,c,a
    real(RK)    :: d(3),q(3,3),o(3,3,3)
    ! temp storage for pc grads
    real(RK), allocatable :: PIPDTRQ(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIPQTRQ(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIPOTRQ(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PIPIDTRQ(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    integer(IK) :: memstat
    real(RK),parameter    :: small=1.0e-8_RK
    !------------ Executable code -----------------------------------

    DPRINT 'SHGI: shgi_X_torq: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

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

    call shgi_set_lcde(1,0,0)

    call shgi_set_ovrl(LA,LB,1,0,0)

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    if(N_pd > 0) then
       FPP_TIMER_START(pd)

       call shgi_set_lcde(1,0,0)
       allocate(PIPDTRQ(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_pd
          n_equal=pds(up)%n_equal_dipoles
          do ep=1,n_equal
             PIPDTRQ = 0.0_rk
             d=pds(up)%dipole(:,ep)
             if(sum(d*d) <= small) cycle
             x=pds(up)%position(:,ep)
             call shgi_trq_pd(d,x,PIPDTRQ)
             call shgi_X_wd_store(PD,up,ep,PIPDTRQ,PDTRQ,-1)
          end do
       end do

       deallocate(PIPDTRQ,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pd)
    end if

    if(N_pq > 0) then
       FPP_TIMER_START(pq)

       call shgi_set_lcde(1,1,0)
       allocate(PIPQTRQ(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_pq
          n_equal=pqs(up)%n_equal_qpoles
          do ep=1,n_equal
             PIPQTRQ = 0.0_rk
             q=pqs(up)%quadrupole(:,:,ep)
             if(sum(q*q) <= small) cycle
             x=pqs(up)%position(:,ep)
             call shgi_trq_pq(q,x,PIPQTRQ)
             call shgi_X_wd_store(PQ,up,ep,PIPQTRQ,PQTRQ,1)
          end do
       end do

       deallocate(PIPQTRQ,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pq)
    end if

    if(N_po > 0) then
       FPP_TIMER_START(po)

       call shgi_set_lcde(1,1,1)
       allocate(PIPOTRQ(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_po
          n_equal=pos(up)%n_equal_opoles
          do ep=1,n_equal
             PIPOTRQ = 0.0_rk
             o=pos(up)%octopole(:,:,:,ep)
             if(sum(o*o) <= small) cycle
             x=pos(up)%position(:,ep)
             call shgi_trq_po(o,x,PIPOTRQ)
             call shgi_X_wd_store(PO,up,ep,PIPOTRQ,POTRQ,-1)
          end do
       end do

       deallocate(PIPOTRQ,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(po)
    end if

    if(N_ipd > 0) then
       FPP_TIMER_START(pd)

       call shgi_set_lcde(1,0,0)
       allocate(PIPIDTRQ(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)

       do up=1,N_ipd
          n_equal=ipds(up)%n_equal_dipoles
          do ep=1,n_equal
             PIPIDTRQ = 0.0_rk
             d=-(ipds(up)%idipole(:,ep)+ipds(up)%idipole1(:,ep))*0.5_RK
             x=ipds(up)%position(:,ep)
             call shgi_trq_pd(d,x,PIPIDTRQ)
             call shgi_X_wd_store(IPD,up,ep,PIPIDTRQ,PIDTRQ,-1)
          end do
       end do

       deallocate(PIPIDTRQ,stat=memstat)
       ASSERT(memstat==0)

       FPP_TIMER_STOP(pd)
    end if

    call shgi_close_ab()

    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_X_torq
!*************************************************************
  subroutine shgi_pd(pd,nu)
    use point_dqo_module, only: pdt=>pointdipole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_PR2v
    use shgi_shr, only: SHR_D3Fv
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(pdt)  , intent(in)    :: pd(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: p2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_pd(',size(pd),'nu)'

    FPP_TIMER_START(pd)

    p2 = 0.0_rk

    do ua=1,size(pd)
       z = 1.0_rk
       do ea=1,pd(ua)%n_equal_dipoles

          if(sum(pd(ua)%dipole(:,ea)*pd(ua)%dipole(:,ea)) <= small) cycle
          x = pd(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpd)
          call doIL(LA+LB+LC,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D3Fv(NAB,LA,LB,LC,wc,fL,f3)
          FPP_TIMER_STOP(td3f)

          p2 = p2 - &
               f3(:,:,:,CX)*pd(ua)%dipole(1,ea)- &
               f3(:,:,:,CY)*pd(ua)%dipole(2,ea)- &
               f3(:,:,:,CZ)*pd(ua)%dipole(3,ea)
       enddo
    enddo

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,1,p2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,p2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)

    FPP_TIMER_STOP(pd)
  end subroutine shgi_pd
!*************************************************************
#if 0
  subroutine shgi_pd1(pd,nu)
    use point_dqo_module, only: pdt=>pointdipole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_PR2v
    use shgi_shr, only: SHR_D3Fv
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(pdt)  , intent(in)    :: pd(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC)
    real(RK)    :: nu1(NAB,2*LA+1,2*LB+1,3)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_pd(',size(pd),'nu)'

    FPP_TIMER_START(pd)

    do ua=1,size(pd)
       z = 1.0_rk
       do ea=1,pd(ua)%n_equal_dipoles

          if(sum(pd(ua)%dipole(:,ea)*pd(ua)%dipole(:,ea)) <= small) cycle
          x = pd(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpd)
          call doIL(LA+LB+LC,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D3Fv(NAB,LA,LB,LC,wc,fL,f3)
          FPP_TIMER_STOP(td3f)

          ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
          call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,f3)

          nu1=zero
          FPP_TIMER_START(tpr2)
          ! Double Product Rule:
          call SHR_PR2v(NAB,LA,LB,f3(:,:,:,CX),S5(:,:,:,C0,C0,C0),nu1(:,:,:,GWX))
          call SHR_PR2v(NAB,LA,LB,f3(:,:,:,CY),S5(:,:,:,C0,C0,C0),nu1(:,:,:,GWY))
          call SHR_PR2v(NAB,LA,LB,f3(:,:,:,CZ),S5(:,:,:,C0,C0,C0),nu1(:,:,:,GWZ))
          FPP_TIMER_STOP(tpr2)

          nu(:,:,:)=nu(:,:,:)- &
               nu1(:,:,:,GWX)*pd(ua)%dipole(1,ea)- &
               nu1(:,:,:,GWY)*pd(ua)%dipole(2,ea)- &
               nu1(:,:,:,GWZ)*pd(ua)%dipole(3,ea)
       enddo
    enddo

  end subroutine shgi_pd1
#endif
!*************************************************************
#if 0
  subroutine shgi_ipd(ipd,nu)
    use induced_dipoles_module, only: ipdt=>ind_dipole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_PR2v
    use shgi_shr, only: SHR_D3Fv
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(ipdt)  , intent(in)    :: ipd(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: p2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC)

    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb
    integer(IK) :: lmc

    DPRINT  'SHGI: shgi_ipd(',size(pd),'nu)'

    FPP_TIMER_START(pd)

    p2 = 0.0_rk

    do ua=1,size(ipd)
       z = 1.0_rk
       do ea=1,ipd(ua)%n_equal_dipoles

          x = ipd(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpd)
          call doIL(LA+LB+LC,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D3Fv(NAB,LA,LB,LC,wc,fL,f3)
          FPP_TIMER_STOP(td3f)

          p2 = p2 - &
               f3(:,:,:,CX)*ipd(ua)%idipole(1,ea)- &
               f3(:,:,:,CY)*ipd(ua)%idipole(2,ea)- &
               f3(:,:,:,CZ)*ipd(ua)%idipole(3,ea)
       enddo
    enddo

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,1,p2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,p2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)

    FPP_TIMER_STOP(pd)
  end subroutine shgi_ipd
#endif
  !*************************************************************
  subroutine shgi_pd_gr(pd,nu)
    use point_dqo_module, only: pdt=>pointdipole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D4Fv
    implicit none
    type(pdt)  , intent(in)    :: pd(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_pd_gr(',size(pd),'nu)'

    FPP_TIMER_START(pd)

    p3 = 0.0_rk

    do ua=1,size(pd)
       z = 1.0_rk
       do ea=1,pd(ua)%n_equal_dipoles

          if(sum(pd(ua)%dipole(:,ea)*pd(ua)%dipole(:,ea)) <= small) cycle
          x = pd(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpd)
          call doIL(LA+LB+LC+LD,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D4Fv(NAB,LA,LB,LC,LD,wc,fL,f4)
          FPP_TIMER_STOP(td3f)

          p3 = p3 - &
               f4(:,:,:,:,CX)*pd(ua)%dipole(1,ea)- &
               f4(:,:,:,:,CY)*pd(ua)%dipole(2,ea)- &
               f4(:,:,:,:,CZ)*pd(ua)%dipole(3,ea)
       enddo
    enddo

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

    FPP_TIMER_STOP(pd)
  end subroutine shgi_pd_gr
  !*************************************************************
  subroutine shgi_ipd_gr(ipd,nu)
    use induced_dipoles_module, only: ipdt=>ind_dipole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D4Fv
    implicit none
    type(ipdt)  , intent(in)    :: ipd(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD)

    DPRINT  'SHGI: shgi_pd_gr(',size(pd),'nu)'

    FPP_TIMER_START(pd)

    p3 = 0.0_rk

    do ua=1,size(ipd)
       z = 0.5_rk
       ! NOTE: overall factor 1/2!
       do ea=1,ipd(ua)%n_equal_dipoles

          x = ipd(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpd)
          call doIL(LA+LB+LC+LD,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D4Fv(NAB,LA,LB,LC,LD,wc,fL,f4)
          FPP_TIMER_STOP(td3f)

          p3 = p3 + &
               f4(:,:,:,:,CX)*(ipd(ua)%idipole(1,ea)+ipd(ua)%idipole1(1,ea))+ &
               f4(:,:,:,:,CY)*(ipd(ua)%idipole(2,ea)+ipd(ua)%idipole1(2,ea))+ &
               f4(:,:,:,:,CZ)*(ipd(ua)%idipole(3,ea)+ipd(ua)%idipole1(3,ea))
      enddo
    enddo

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

    FPP_TIMER_STOP(pd)
  end subroutine shgi_ipd_gr
  !*************************************************************
  subroutine shgi_gr_pd(d,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D4Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: d(3),x(3)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD)

    p3 = 0.0_rk

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    FPP_TIMER_START(tpd)
    call doIL(LA+LB+LC+LD,w2,LAMBDA,NORM(:,3),fL)
    FPP_TIMER_STOP(tpd)

    FPP_TIMER_START(td3f)
    call SHR_D4Fv(NAB,LA,LB,LC,LD,wc,fL,f4)
    FPP_TIMER_STOP(td3f)

    p3 = p3 - &
         f4(:,:,:,:,CX)*d(1)- &
         f4(:,:,:,:,CY)*d(2)- &
         f4(:,:,:,:,CZ)*d(3)

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_gr_pd
  !*************************************************************
  subroutine shgi_trq_pd(d,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D3Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: d(3),x(3)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC)

    p3 = 0.0_rk

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    FPP_TIMER_START(tpd)
    call doIL(LA+LB+LC,w2,LAMBDA,NORM(:,3),fL)
    FPP_TIMER_STOP(tpd)

    FPP_TIMER_START(td3f)
    call SHR_D3Fv(NAB,LA,LB,LC,wc,fL,f3)
    FPP_TIMER_STOP(td3f)

    p3(:,:,:,CX)=p3(:,:,:,CX)+(d(2)*f3(:,:,:,CZ)-d(3)*f3(:,:,:,CY))
    p3(:,:,:,CY)=p3(:,:,:,CY)+(d(3)*f3(:,:,:,CX)-d(1)*f3(:,:,:,CZ))
    p3(:,:,:,CZ)=p3(:,:,:,CZ)+(d(1)*f3(:,:,:,CY)-d(2)*f3(:,:,:,CX))

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_trq_pd
  !*************************************************************
  subroutine shgi_pq(pq,nu)
    use point_dqo_module, only: pqt=>pointquadrupole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_PR2v
    use shgi_shr, only: SHR_D4Fv
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(pqt)  , intent(in)    :: pq(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)    :: p2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_pq(',size(pq),'nu)'

    FPP_TIMER_START(pq)

    p2 = 0.0_rk

    do ua=1,size(pq)
       z = -1.0_rk/3.0_rk
       ! NOTE: overall factor -1/3!
       do ea=1,pq(ua)%n_equal_qpoles

          if(sum(pq(ua)%quadrupole(:,:,ea)*pq(ua)%quadrupole(:,:,ea)) <= small) cycle
          x = pq(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpq)
          call doIL(LA+LB+LC+LD,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpq)


          FPP_TIMER_START(td3f)
          call SHR_D4Fv(NAB,LA,LB,LC,LD,wc,fL,f4)
          FPP_TIMER_STOP(td3f)

          p2 = p2 - &
               f4(:,:,:,CX,CX)*pq(ua)%quadrupole(1,1,ea)- &
               f4(:,:,:,CX,CY)*pq(ua)%quadrupole(1,2,ea)- &
               f4(:,:,:,CX,CZ)*pq(ua)%quadrupole(1,3,ea)- &
               f4(:,:,:,CY,CX)*pq(ua)%quadrupole(2,1,ea)- &
               f4(:,:,:,CY,CY)*pq(ua)%quadrupole(2,2,ea)- &
               f4(:,:,:,CY,CZ)*pq(ua)%quadrupole(2,3,ea)- &
               f4(:,:,:,CZ,CX)*pq(ua)%quadrupole(3,1,ea)- &
               f4(:,:,:,CZ,CY)*pq(ua)%quadrupole(3,2,ea)- &
               f4(:,:,:,CZ,CZ)*pq(ua)%quadrupole(3,3,ea)
       enddo
    enddo

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,1,p2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,p2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)

    FPP_TIMER_STOP(pq)
  end subroutine shgi_pq
  !*************************************************************
  subroutine shgi_pq_gr(pq,nu)
    use point_dqo_module, only: pqt=>pointquadrupole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D5Fv
    implicit none
    type(pqt)  , intent(in)    :: pq(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f5(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD+LE)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_pq_gr(',size(pq),'nu)'

    FPP_TIMER_START(pq)

    p3 = 0.0_rk

    do ua=1,size(pq)
       z = -1.0_rk/3.0_rk
       ! NOTE: overall factor -1/3!
       do ea=1,pq(ua)%n_equal_qpoles

          if(sum(pq(ua)%quadrupole(:,:,ea)*pq(ua)%quadrupole(:,:,ea)) <= small) cycle
          x = pq(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpq)
          call doIL(LA+LB+LC+LD+LE,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpq)

          FPP_TIMER_START(td3f)
          call SHR_D5Fv(NAB,LA,LB,LC,LD,LE,wc,fL,f5)
          FPP_TIMER_STOP(td3f)

          p3 = p3 - &
               f5(:,:,:,:,CX,CX)*pq(ua)%quadrupole(1,1,ea)- &
               f5(:,:,:,:,CX,CY)*pq(ua)%quadrupole(1,2,ea)- &
               f5(:,:,:,:,CX,CZ)*pq(ua)%quadrupole(1,3,ea)- &
               f5(:,:,:,:,CY,CX)*pq(ua)%quadrupole(2,1,ea)- &
               f5(:,:,:,:,CY,CY)*pq(ua)%quadrupole(2,2,ea)- &
               f5(:,:,:,:,CY,CZ)*pq(ua)%quadrupole(2,3,ea)- &
               f5(:,:,:,:,CZ,CX)*pq(ua)%quadrupole(3,1,ea)- &
               f5(:,:,:,:,CZ,CY)*pq(ua)%quadrupole(3,2,ea)- &
               f5(:,:,:,:,CZ,CZ)*pq(ua)%quadrupole(3,3,ea)
       enddo
    enddo

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)
    FPP_TIMER_STOP(pq)
  end subroutine shgi_pq_gr
  !*************************************************************
  subroutine shgi_gr_pq(q,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D5Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: q(3,3),x(3)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f5(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD+LE)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    FPP_TIMER_START(tpq)
    call doIL(LA+LB+LC+LD+LE,w2,LAMBDA,-NORM(:,3)/3.0_RK,fL)
    ! NOTE: overall factor -1/3!
    FPP_TIMER_STOP(tpq)

    FPP_TIMER_START(td3f)
    call SHR_D5Fv(NAB,LA,LB,LC,LD,LE,wc,fL,f5)
    FPP_TIMER_STOP(td3f)

    p3= -f5(:,:,:,:,CX,CX)*q(1,1)- &
         f5(:,:,:,:,CX,CY)*q(1,2)- &
         f5(:,:,:,:,CX,CZ)*q(1,3)- &
         f5(:,:,:,:,CY,CX)*q(2,1)- &
         f5(:,:,:,:,CY,CY)*q(2,2)- &
         f5(:,:,:,:,CY,CZ)*q(2,3)- &
         f5(:,:,:,:,CZ,CX)*q(3,1)- &
         f5(:,:,:,:,CZ,CY)*q(3,2)- &
         f5(:,:,:,:,CZ,CZ)*q(3,3)

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_gr_pq
  !*************************************************************
  subroutine shgi_trq_pq(q,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D4Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: q(3,3),x(3)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD)

    p3 = 0.0_rk

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    FPP_TIMER_START(tpd)
    call doIL(LA+LB+LC+LD,w2,LAMBDA,2.0_RK*NORM(:,3)/3.0_RK,fL)
    ! NOTE: overall factor 2/3!
    FPP_TIMER_STOP(tpd)

    FPP_TIMER_START(td3f)
    call SHR_D4Fv(NAB,LA,LB,LC,LD,wc,fL,f4)
    FPP_TIMER_STOP(td3f)

    p3(:,:,:,CX)=p3(:,:,:,CX)+(q(2,1)*f4(:,:,:,CX,CZ)+   &
                               q(2,2)*f4(:,:,:,CY,CZ)+   &
                               q(2,3)*f4(:,:,:,CZ,CZ)) - &
                              (q(3,1)*f4(:,:,:,CX,CY)+   &
                               q(3,2)*f4(:,:,:,CY,CY)+   &
                               q(3,3)*f4(:,:,:,CZ,CY))

    p3(:,:,:,CY)=p3(:,:,:,CY)+(q(3,1)*f4(:,:,:,CX,CX)+   &
                               q(3,2)*f4(:,:,:,CY,CX)+   &
                               q(3,3)*f4(:,:,:,CZ,CX)) - &
                              (q(1,1)*f4(:,:,:,CX,CZ)+   &
                               q(1,2)*f4(:,:,:,CY,CZ)+   &
                               q(1,3)*f4(:,:,:,CZ,CZ))

    p3(:,:,:,CZ)=p3(:,:,:,CZ)+(q(1,1)*f4(:,:,:,CX,CY)+   &
                               q(1,2)*f4(:,:,:,CY,CY)+   &
                               q(1,3)*f4(:,:,:,CZ,CY)) - &
                              (q(2,1)*f4(:,:,:,CX,CX)+   &
                               q(2,2)*f4(:,:,:,CY,CX)+   &
                               q(2,3)*f4(:,:,:,CZ,CX))

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_trq_pq
  !*************************************************************
  subroutine shgi_po(po,nu)
    use point_dqo_module, only: pot=>pointoctopole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_PR2v
    use shgi_shr, only: SHR_D5Fv
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(pot)  , intent(in)    :: po(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f5(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
    real(RK)    :: p2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD+LE)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_po(',size(po),'nu)'

    FPP_TIMER_START(po)

    p2 = 0.0_rk

    do ua=1,size(po)
       z = 1.0_rk/15.0_rk
       ! NOTE: overall factor 1/15!
       do ea=1,po(ua)%n_equal_opoles

          if(sum(po(ua)%octopole(:,:,:,ea)*po(ua)%octopole(:,:,:,ea)) <= small) cycle
          x = po(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpo)
          call doIL(LA+LB+LC+LD+LE,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpo)

          FPP_TIMER_START(td3f)
          call SHR_D5Fv(NAB,LA,LB,LC,LD,LE,wc,fL,f5)
          FPP_TIMER_STOP(td3f)

          p2 = p2 - &
               f5(:,:,:,CX,CX,CX)*po(ua)%octopole(1,1,1,ea)- &
               f5(:,:,:,CX,CY,CX)*po(ua)%octopole(1,2,1,ea)- &
               f5(:,:,:,CX,CZ,CX)*po(ua)%octopole(1,3,1,ea)- &
               f5(:,:,:,CY,CX,CX)*po(ua)%octopole(2,1,1,ea)- &
               f5(:,:,:,CY,CY,CX)*po(ua)%octopole(2,2,1,ea)- &
               f5(:,:,:,CY,CZ,CX)*po(ua)%octopole(2,3,1,ea)- &
               f5(:,:,:,CZ,CX,CX)*po(ua)%octopole(3,1,1,ea)- &
               f5(:,:,:,CZ,CY,CX)*po(ua)%octopole(3,2,1,ea)- &
               f5(:,:,:,CZ,CZ,CX)*po(ua)%octopole(3,3,1,ea)- &
               f5(:,:,:,CX,CX,CY)*po(ua)%octopole(1,1,2,ea)- &
               f5(:,:,:,CX,CY,CY)*po(ua)%octopole(1,2,2,ea)- &
               f5(:,:,:,CX,CZ,CY)*po(ua)%octopole(1,3,2,ea)- &
               f5(:,:,:,CY,CX,CY)*po(ua)%octopole(2,1,2,ea)- &
               f5(:,:,:,CY,CY,CY)*po(ua)%octopole(2,2,2,ea)- &
               f5(:,:,:,CY,CZ,CY)*po(ua)%octopole(2,3,2,ea)- &
               f5(:,:,:,CZ,CX,CY)*po(ua)%octopole(3,1,2,ea)- &
               f5(:,:,:,CZ,CY,CY)*po(ua)%octopole(3,2,2,ea)- &
               f5(:,:,:,CZ,CZ,CY)*po(ua)%octopole(3,3,2,ea)- &
               f5(:,:,:,CX,CX,CZ)*po(ua)%octopole(1,1,3,ea)- &
               f5(:,:,:,CX,CY,CZ)*po(ua)%octopole(1,2,3,ea)- &
               f5(:,:,:,CX,CZ,CZ)*po(ua)%octopole(1,3,3,ea)- &
               f5(:,:,:,CY,CX,CZ)*po(ua)%octopole(2,1,3,ea)- &
               f5(:,:,:,CY,CY,CZ)*po(ua)%octopole(2,2,3,ea)- &
               f5(:,:,:,CY,CZ,CZ)*po(ua)%octopole(2,3,3,ea)- &
               f5(:,:,:,CZ,CX,CZ)*po(ua)%octopole(3,1,3,ea)- &
               f5(:,:,:,CZ,CY,CZ)*po(ua)%octopole(3,2,3,ea)- &
               f5(:,:,:,CZ,CZ,CZ)*po(ua)%octopole(3,3,3,ea)
       end do
    end do

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,1,p2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,p2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)

    FPP_TIMER_STOP(po)
  end subroutine shgi_po
  !*************************************************************
  subroutine shgi_po_gr(po,nu)
    use point_dqo_module, only: pot=>pointoctopole_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D6Fv
    implicit none
    type(pot)  , intent(in)    :: po(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z, x(3)

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f6(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD+LE+LF)

    real(RK),parameter    :: small=1.0e-8_RK

    DPRINT  'SHGI: shgi_po_gr(',size(po),'nu)'

    FPP_TIMER_START(po)

    p3 = 0.0_rk

    do ua=1,size(po)
       z = 1.0_rk/15.0_rk
       ! NOTE: overall factor 1/15!
       do ea=1,po(ua)%n_equal_opoles

          if(sum(po(ua)%octopole(:,:,:,ea)*po(ua)%octopole(:,:,:,ea)) <= small) cycle
          x = po(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          FPP_TIMER_START(tpo)
          call doIL(LA+LB+LC+LD+LE+LF,w2,LAMBDA,z*NORM(:,3),fL)
          FPP_TIMER_STOP(tpo)


          FPP_TIMER_START(td3f)
          call SHR_D6Fv(NAB,LA,LB,LC,LD,LE,LF,wc,fL,f6)
          FPP_TIMER_STOP(td3f)

          p3 = p3 - &
               f6(:,:,:,:,CX,CX,CX)*po(ua)%octopole(1,1,1,ea)- &
               f6(:,:,:,:,CX,CY,CX)*po(ua)%octopole(1,2,1,ea)- &
               f6(:,:,:,:,CX,CZ,CX)*po(ua)%octopole(1,3,1,ea)- &
               f6(:,:,:,:,CY,CX,CX)*po(ua)%octopole(2,1,1,ea)- &
               f6(:,:,:,:,CY,CY,CX)*po(ua)%octopole(2,2,1,ea)- &
               f6(:,:,:,:,CY,CZ,CX)*po(ua)%octopole(2,3,1,ea)- &
               f6(:,:,:,:,CZ,CX,CX)*po(ua)%octopole(3,1,1,ea)- &
               f6(:,:,:,:,CZ,CY,CX)*po(ua)%octopole(3,2,1,ea)- &
               f6(:,:,:,:,CZ,CZ,CX)*po(ua)%octopole(3,3,1,ea)- &
               f6(:,:,:,:,CX,CX,CY)*po(ua)%octopole(1,1,2,ea)- &
               f6(:,:,:,:,CX,CY,CY)*po(ua)%octopole(1,2,2,ea)- &
               f6(:,:,:,:,CX,CZ,CY)*po(ua)%octopole(1,3,2,ea)- &
               f6(:,:,:,:,CY,CX,CY)*po(ua)%octopole(2,1,2,ea)- &
               f6(:,:,:,:,CY,CY,CY)*po(ua)%octopole(2,2,2,ea)- &
               f6(:,:,:,:,CY,CZ,CY)*po(ua)%octopole(2,3,2,ea)- &
               f6(:,:,:,:,CZ,CX,CY)*po(ua)%octopole(3,1,2,ea)- &
               f6(:,:,:,:,CZ,CY,CY)*po(ua)%octopole(3,2,2,ea)- &
               f6(:,:,:,:,CZ,CZ,CY)*po(ua)%octopole(3,3,2,ea)- &
               f6(:,:,:,:,CX,CX,CZ)*po(ua)%octopole(1,1,3,ea)- &
               f6(:,:,:,:,CX,CY,CZ)*po(ua)%octopole(1,2,3,ea)- &
               f6(:,:,:,:,CX,CZ,CZ)*po(ua)%octopole(1,3,3,ea)- &
               f6(:,:,:,:,CY,CX,CZ)*po(ua)%octopole(2,1,3,ea)- &
               f6(:,:,:,:,CY,CY,CZ)*po(ua)%octopole(2,2,3,ea)- &
               f6(:,:,:,:,CY,CZ,CZ)*po(ua)%octopole(2,3,3,ea)- &
               f6(:,:,:,:,CZ,CX,CZ)*po(ua)%octopole(3,1,3,ea)- &
               f6(:,:,:,:,CZ,CY,CZ)*po(ua)%octopole(3,2,3,ea)- &
               f6(:,:,:,:,CZ,CZ,CZ)*po(ua)%octopole(3,3,3,ea)
       end do
    end do

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)
    FPP_TIMER_STOP(po)
  end subroutine shgi_po_gr
  !*************************************************************
  subroutine shgi_gr_po(o,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D6Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: o(3,3,3),x(3)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f6(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,(LF+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD+LE+LF)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    FPP_TIMER_START(tpo)
    call doIL(LA+LB+LC+LD+LE+LF,w2,LAMBDA,NORM(:,3)/15.0_RK,fL)
    ! NOTE: overall factor 1/15!
    FPP_TIMER_STOP(tpo)


    FPP_TIMER_START(td3f)
    call SHR_D6Fv(NAB,LA,LB,LC,LD,LE,LF,wc,fL,f6)
    FPP_TIMER_STOP(td3f)

    p3= -f6(:,:,:,:,CX,CX,CX)*o(1,1,1)- &
         f6(:,:,:,:,CX,CY,CX)*o(1,2,1)- &
         f6(:,:,:,:,CX,CZ,CX)*o(1,3,1)- &
         f6(:,:,:,:,CY,CX,CX)*o(2,1,1)- &
         f6(:,:,:,:,CY,CY,CX)*o(2,2,1)- &
         f6(:,:,:,:,CY,CZ,CX)*o(2,3,1)- &
         f6(:,:,:,:,CZ,CX,CX)*o(3,1,1)- &
         f6(:,:,:,:,CZ,CY,CX)*o(3,2,1)- &
         f6(:,:,:,:,CZ,CZ,CX)*o(3,3,1)- &
         f6(:,:,:,:,CX,CX,CY)*o(1,1,2)- &
         f6(:,:,:,:,CX,CY,CY)*o(1,2,2)- &
         f6(:,:,:,:,CX,CZ,CY)*o(1,3,2)- &
         f6(:,:,:,:,CY,CX,CY)*o(2,1,2)- &
         f6(:,:,:,:,CY,CY,CY)*o(2,2,2)- &
         f6(:,:,:,:,CY,CZ,CY)*o(2,3,2)- &
         f6(:,:,:,:,CZ,CX,CY)*o(3,1,2)- &
         f6(:,:,:,:,CZ,CY,CY)*o(3,2,2)- &
         f6(:,:,:,:,CZ,CZ,CY)*o(3,3,2)- &
         f6(:,:,:,:,CX,CX,CZ)*o(1,1,3)- &
         f6(:,:,:,:,CX,CY,CZ)*o(1,2,3)- &
         f6(:,:,:,:,CX,CZ,CZ)*o(1,3,3)- &
         f6(:,:,:,:,CY,CX,CZ)*o(2,1,3)- &
         f6(:,:,:,:,CY,CY,CZ)*o(2,2,3)- &
         f6(:,:,:,:,CY,CZ,CZ)*o(2,3,3)- &
         f6(:,:,:,:,CZ,CX,CZ)*o(3,1,3)- &
         f6(:,:,:,:,CZ,CY,CZ)*o(3,2,3)- &
         f6(:,:,:,:,CZ,CZ,CZ)*o(3,3,3)

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_gr_po
  !*************************************************************
  subroutine shgi_trq_po(o,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doIL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D5Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: o(3,3,3),x(3)
    ! *** end of interface ***

    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f5(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC+LD+LE)

    p3 = 0.0_rk

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    FPP_TIMER_START(tpo)
    call doIL(LA+LB+LC+LD+LE,w2,LAMBDA,NORM(:,3)/15.0_RK,fL)
    ! NOTE: overall factor 1/15!
    FPP_TIMER_STOP(tpo)

    FPP_TIMER_START(td3f)
    call SHR_D5Fv(NAB,LA,LB,LC,LD,LE,wc,fL,f5)
    FPP_TIMER_STOP(td3f)

    p3(:,:,:,CX)=p3(:,:,:,CX)+o(1,1,2)*  f5(:,:,:,CX,CX,CZ)                                        &
                             -o(1,1,3)*  f5(:,:,:,CX,CX,CY)                                        &
                             +o(1,2,1)*  f5(:,:,:,CX,CZ,CX)                                        &
                             +o(1,2,2)*( f5(:,:,:,CX,CZ,CY)+f5(:,:,:,CX,CY,CZ))                    &
                             +o(1,2,3)*( f5(:,:,:,CX,CZ,CZ)-f5(:,:,:,CX,CY,CY))                    &
                             -o(1,3,1)*  f5(:,:,:,CX,CY,CX)                                        &
                             +o(1,3,2)*(-f5(:,:,:,CX,CY,CY)+f5(:,:,:,CX,CZ,CZ))                    &
                             +o(1,3,3)*(-f5(:,:,:,CX,CY,CZ)-f5(:,:,:,CX,CZ,CY))                    &
                             +o(2,1,1)*  f5(:,:,:,CZ,CX,CX)                                        &
                             +o(2,1,2)*( f5(:,:,:,CZ,CX,CY)+f5(:,:,:,CY,CX,CZ))                    &
                             +o(2,1,3)*( f5(:,:,:,CZ,CX,CZ)-f5(:,:,:,CY,CX,CY))                    &
                             +o(2,2,1)*( f5(:,:,:,CZ,CY,CX)+f5(:,:,:,CY,CZ,CX))                    &
                             +o(2,2,2)*( f5(:,:,:,CZ,CY,CY)+f5(:,:,:,CY,CZ,CY)+f5(:,:,:,CY,CY,CZ)) &
                             +o(2,2,3)*( f5(:,:,:,CZ,CY,CZ)+f5(:,:,:,CY,CZ,CZ)-f5(:,:,:,CY,CY,CY)) &
                             +o(2,3,1)*( f5(:,:,:,CZ,CZ,CX)-f5(:,:,:,CY,CY,CX))                    &
                             +o(2,3,2)*( f5(:,:,:,CZ,CZ,CY)-f5(:,:,:,CY,CY,CY)+f5(:,:,:,CY,CZ,CZ)) &
                             +o(2,3,3)*( f5(:,:,:,CZ,CZ,CZ)-f5(:,:,:,CY,CY,CZ)-f5(:,:,:,CY,CZ,CY)) &
                             -o(3,1,1)*  f5(:,:,:,CY,CX,CX)                                        &
                             +o(3,1,2)*(-f5(:,:,:,CY,CX,CY)+f5(:,:,:,CZ,CX,CZ))                    &
                             +o(3,1,3)*(-f5(:,:,:,CY,CX,CZ)-f5(:,:,:,CZ,CX,CY))                    &
                             +o(3,2,1)*(-f5(:,:,:,CY,CY,CX)+f5(:,:,:,CZ,CZ,CX))                    &
                             +o(3,2,2)*(-f5(:,:,:,CY,CY,CY)+f5(:,:,:,CZ,CZ,CY)+f5(:,:,:,CZ,CY,CZ)) &
                             +o(3,2,3)*(-f5(:,:,:,CY,CY,CZ)+f5(:,:,:,CZ,CZ,CZ)-f5(:,:,:,CZ,CY,CY)) &
                             +o(3,3,1)*(-f5(:,:,:,CY,CZ,CX)-f5(:,:,:,CZ,CY,CX))                    &
                             +o(3,3,2)*(-f5(:,:,:,CY,CZ,CY)-f5(:,:,:,CZ,CY,CY)+f5(:,:,:,CZ,CZ,CZ)) &
                             +o(3,3,3)*(-f5(:,:,:,CY,CZ,CZ)-f5(:,:,:,CZ,CY,CZ)-f5(:,:,:,CZ,CZ,CY))

    p3(:,:,:,CY)=p3(:,:,:,CY)+o(1,1,1)*(-f5(:,:,:,CZ,CX,CX)-f5(:,:,:,CX,CZ,CX)-f5(:,:,:,CX,CX,CZ)) &
                             +o(1,1,2)*(-f5(:,:,:,CZ,CX,CY)-f5(:,:,:,CX,CZ,CY))                    &
                             +o(1,1,3)*(-f5(:,:,:,CZ,CX,CZ)-f5(:,:,:,CX,CZ,CZ)+f5(:,:,:,CX,CX,CX)) &
                             +o(1,2,1)*(-f5(:,:,:,CZ,CY,CX)-f5(:,:,:,CX,CY,CZ))                    &
                             -o(1,2,2)*  f5(:,:,:,CZ,CY,CY)                                        &
                             +o(1,2,3)*(-f5(:,:,:,CZ,CY,CZ)+f5(:,:,:,CX,CY,CX))                    &
                             +o(1,3,1)*(-f5(:,:,:,CZ,CZ,CX)+f5(:,:,:,CX,CX,CX)-f5(:,:,:,CX,CZ,CZ)) &
                             +o(1,3,2)*(-f5(:,:,:,CZ,CZ,CY)+f5(:,:,:,CX,CX,CY))                    &
                             +o(1,3,3)*(-f5(:,:,:,CZ,CZ,CZ)+f5(:,:,:,CX,CX,CZ)+f5(:,:,:,CX,CZ,CX)) &
                             +o(2,1,1)*(-f5(:,:,:,CY,CZ,CX)-f5(:,:,:,CY,CX,CZ))                    &
                             -o(2,1,2)*  f5(:,:,:,CY,CZ,CY)                                        &
                             +o(2,1,3)*(-f5(:,:,:,CY,CZ,CZ)+f5(:,:,:,CY,CX,CX))                    &
                             -o(2,2,1)*  f5(:,:,:,CY,CY,CZ)                                        &
                             +o(2,2,3)*  f5(:,:,:,CY,CY,CX)                                        &
                             +o(2,3,1)*( f5(:,:,:,CY,CX,CX)-f5(:,:,:,CY,CZ,CZ))                    &
                             +o(2,3,2)*  f5(:,:,:,CY,CX,CY)                                        &
                             +o(2,3,3)*( f5(:,:,:,CY,CX,CZ)+f5(:,:,:,CY,CZ,CX))                    &
                             +o(3,1,1)*( f5(:,:,:,CX,CX,CX)-f5(:,:,:,CZ,CZ,CX)-f5(:,:,:,CZ,CX,CZ)) &
                             +o(3,1,2)*( f5(:,:,:,CX,CX,CY)-f5(:,:,:,CZ,CZ,CY))                    &
                             +o(3,1,3)*( f5(:,:,:,CX,CX,CZ)-f5(:,:,:,CZ,CZ,CZ)+f5(:,:,:,CZ,CX,CX)) &
                             +o(3,2,1)*( f5(:,:,:,CX,CY,CX)-f5(:,:,:,CZ,CY,CZ))                    &
                             +o(3,2,2)*  f5(:,:,:,CX,CY,CY)                                        &
                             +o(3,2,3)*( f5(:,:,:,CX,CY,CZ)+f5(:,:,:,CZ,CY,CX))                    &
                             +o(3,3,1)*( f5(:,:,:,CX,CZ,CX)+f5(:,:,:,CZ,CX,CX)-f5(:,:,:,CZ,CZ,CZ)) &
                             +o(3,3,2)*( f5(:,:,:,CX,CZ,CY)+f5(:,:,:,CZ,CX,CY))                    &
                             +o(3,3,3)*( f5(:,:,:,CX,CZ,CZ)+f5(:,:,:,CZ,CX,CZ)+f5(:,:,:,CZ,CZ,CX))

    p3(:,:,:,CZ)=p3(:,:,:,CZ)+o(1,1,1)*( f5(:,:,:,CY,CX,CX)+f5(:,:,:,CX,CY,CX)+f5(:,:,:,CX,CX,CY)) &
                             +o(1,1,2)*( f5(:,:,:,CY,CX,CY)+f5(:,:,:,CX,CY,CY)-f5(:,:,:,CX,CX,CX)) &
                             +o(1,1,3)*( f5(:,:,:,CY,CX,CZ)+f5(:,:,:,CX,CY,CZ))                    &
                             +o(1,2,1)*( f5(:,:,:,CY,CY,CX)-f5(:,:,:,CX,CX,CX)+f5(:,:,:,CX,CY,CY)) &
                             +o(1,2,2)*( f5(:,:,:,CY,CY,CY)-f5(:,:,:,CX,CX,CY)-f5(:,:,:,CX,CY,CX)) &
                             +o(1,2,3)*( f5(:,:,:,CY,CY,CZ)-f5(:,:,:,CX,CX,CZ))                    &
                             +o(1,3,1)*( f5(:,:,:,CY,CZ,CX)+f5(:,:,:,CX,CZ,CY))                    &
                             +o(1,3,2)*( f5(:,:,:,CY,CZ,CY)-f5(:,:,:,CX,CZ,CX))                    &
                             +o(1,3,3)*  f5(:,:,:,CY,CZ,CZ)                                        &
                             +o(2,1,1)*(-f5(:,:,:,CX,CX,CX)+f5(:,:,:,CY,CY,CX)+f5(:,:,:,CY,CX,CY)) &
                             +o(2,1,2)*(-f5(:,:,:,CX,CX,CY)+f5(:,:,:,CY,CY,CY)-f5(:,:,:,CY,CX,CX)) &
                             +o(2,1,3)*(-f5(:,:,:,CX,CX,CZ)+f5(:,:,:,CY,CY,CZ))                    &
                             +o(2,2,1)*(-f5(:,:,:,CX,CY,CX)-f5(:,:,:,CY,CX,CX)+f5(:,:,:,CY,CY,CY)) &
                             +o(2,2,2)*(-f5(:,:,:,CX,CY,CY)-f5(:,:,:,CY,CX,CY)-f5(:,:,:,CY,CY,CX)) &
                             +o(2,2,3)*(-f5(:,:,:,CX,CY,CZ)-f5(:,:,:,CY,CX,CZ))                    &
                             +o(2,3,1)*(-f5(:,:,:,CX,CZ,CX)+f5(:,:,:,CY,CZ,CY))                    &
                             +o(2,3,2)*(-f5(:,:,:,CX,CZ,CY)-f5(:,:,:,CY,CZ,CX))                    &
                             -o(2,3,3)*  f5(:,:,:,CX,CZ,CZ)                                        &
                             +o(3,1,1)*( f5(:,:,:,CZ,CY,CX)+f5(:,:,:,CZ,CX,CY))                    &
                             +o(3,1,2)*( f5(:,:,:,CZ,CY,CY)-f5(:,:,:,CZ,CX,CX))                    &
                             +o(3,1,3)*  f5(:,:,:,CZ,CY,CZ)                                        &
                             +o(3,2,1)*(-f5(:,:,:,CZ,CX,CX)+f5(:,:,:,CZ,CY,CY))                    &
                             +o(3,2,2)*(-f5(:,:,:,CZ,CX,CY)-f5(:,:,:,CZ,CY,CX))                    &
                             -o(3,2,3)*  f5(:,:,:,CZ,CX,CZ)                                        &
                             +o(3,3,1)*  f5(:,:,:,CZ,CZ,CY)                                        &
                             -o(3,3,2)*  f5(:,:,:,CZ,CZ,CX)

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_trq_po
  !*************************************************************
  subroutine shgi_rc(rc,nu)
    use point_dqo_module, only: rct=>repulsion_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doEL
    use shgi_shr, only: SHR_PR2v
    use shgi_shr, only: SHR_D2Fv
    use shgi_utils, only: shgi_dww_to_dab
    implicit none
    type(rct)  , intent(in)    :: rc(:)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z,g,x(3)

    real(RK)    :: ny(NAB)
    real(RK)    :: LAMD(NAB)
    real(RK)    :: wc(NAB,3)
    real(RK)    :: NRM(NAB)
    real(RK)    :: w2(NAB)
    real(RK)    :: f2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: p2(NAB,(LA+1)**2,(LB+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB)

    DPRINT  'SHGI: shgi_po(',size(rc),'nu)'

    FPP_TIMER_START(pd)

    p2 = 0.0_rk

    do ua=1,size(rc)
       z=rc(ua)%C
       g=rc(ua)%A
       do ea=1,rc(ua)%n_equal_centers

          x = rc(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          ny=LAMBDA/(LAMBDA+g)
          LAMD=ny*g
          NRM=NORM(:,2)*ny**(3.0_RK/2.0_RK)

          FPP_TIMER_START(tpd)
          call doEL(LA+LB,w2,LAMD,NRM*z,fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D2Fv(NAB,LA,LB,wc,fL,f2)
          FPP_TIMER_STOP(td3f)

          p2 = p2 - f2
       end do
    end do

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,1,p2)

    FPP_TIMER_START(tpr2)
    ! Double Product Rule:
    call SHR_PR2v(NAB,LA,LB,p2,S5(:,:,:,C0,C0,C0),nu)
    FPP_TIMER_STOP(tpr2)

    FPP_TIMER_STOP(pd)
  end subroutine shgi_rc
  !*************************************************************
  subroutine shgi_rc_gr(rc,nu)
    use point_dqo_module, only: rct=>repulsion_type
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doEL
    use shgi_shr, only: SHR_D3Fv
    use shgi_utils, only: grFS, shgi_dww_to_dab
    implicit none
    type(rct)  , intent(in)    :: rc(:)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    integer(IK) :: ua,ea
    real(RK)    :: z,g,x(3)

    real(RK)    :: ny(NAB)
    real(RK)    :: LAMD(NAB)
    real(RK)    :: wc(NAB,3)
    real(RK)    :: NRM(NAB)
    real(RK)    :: w2(NAB)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: p3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC)

    DPRINT  'SHGI: shgi_po(',size(rc),'nu)'

    FPP_TIMER_START(pd)

    p3 = 0.0_rk

    do ua=1,size(rc)
       z=rc(ua)%C
       g=rc(ua)%A
       do ea=1,rc(ua)%n_equal_centers

          x = rc(ua)%position(:,ea)

          wc(:,1) = W(:,1) - x(1)
          wc(:,2) = W(:,2) - x(2)
          wc(:,3) = W(:,3) - x(3)
          w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

          ny=LAMBDA/(LAMBDA+g)
          LAMD=ny*g
          NRM=NORM(:,2)*ny**(3.0_RK/2.0_RK)

          FPP_TIMER_START(tpd)
          call doEL(LA+LB+LC,w2,LAMD,NRM*z,fL)
          FPP_TIMER_STOP(tpd)

          FPP_TIMER_START(td3f)
          call SHR_D3Fv(NAB,LA,LB,LC,wc,fL,f3)
          FPP_TIMER_STOP(td3f)

          p3 = p3 - f3
       end do
    end do

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,p3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(p3,S5(:,:,:,:,C0,C0),nu)
    FPP_TIMER_STOP(pd)
  end subroutine shgi_rc_gr
  !*************************************************************
  subroutine shgi_gr_rc(c,a,x,nu)
    use shgi_common, only: W, LAMBDA, NORM, S5
    use shgi_rad, only: doEL
    use shgi_utils, only: grFS, shgi_dww_to_dab
    use shgi_shr, only: SHR_D3Fv
    implicit none
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: c,a,x(3)
    ! *** end of interface ***

    real(RK)    :: ny(NAB)
    real(RK)    :: LAMD(NAB)
    real(RK)    :: NRM(NAB)
    real(RK)    :: wc(NAB,3)
    real(RK)    :: w2(NAB)
    real(RK)    :: f3(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)    :: fL(NAB,1+LA+LB+LC)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    ny=LAMBDA/(LAMBDA+a)
    LAMD=ny*a
    NRM = - NORM(:,2)*ny**(3.0_RK/2.0_RK)
    ! NOTE: minus for overall sign!

    FPP_TIMER_START(tpd)
    call doEL(LA+LB+LC,w2,LAMD,NRM*c,fL)
    FPP_TIMER_STOP(tpd)

    FPP_TIMER_START(td3f)
    call SHR_D3Fv(NAB,LA,LB,LC,wc,fL,f3)
    FPP_TIMER_STOP(td3f)

    ! D(lma)D(lmb) ... F => DA(lma)DB(lmb) ... F:
    call shgi_dww_to_dab(NAB,LA,LB,(LC+1)**2,f3)

    ! F(g)*S(0), F(0)*S(g):
    call grFS(f3,S5(:,:,:,:,C0,C0),nu)

  end subroutine shgi_gr_rc
  !*************************************************************
  !--------------- End of module ----------------------------------
end module shgi_ext_c

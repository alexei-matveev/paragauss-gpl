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
module shgi_dip
  !-------------------------------------------------------------------
  !
  ! For  dipole integral of  LL, LS  and SS  type.  Module  called by:
  ! integral_calc_quad_dipole.f90
  !
  ! Copyright (c) 2009 Gopal Dixit
  ! Copyright (c) 2009 Alexei Matveev
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
! USE debug
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: shgi_dip_drv
  public :: shgi_gr_efield_drv

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

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  ! a copy from shgi_shr.f90:
  integer(IK)             :: L_,M_
  integer(IK), parameter  :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter  :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)

  !------------ Subroutines ------------------------------------------
contains

  !****************************************************************
  !****************** DRIVER FOR DIPOLE INTEGRALS *****************
  !****************************************************************

  subroutine shgi_dip_drv(IU1,IE1,IL1,IU2,IE2,IL2, uas, DIPL)
    use unique_atom_module, only: uat=>unique_atom_type
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: upack2c, shgi_timing
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)  :: uas(:)            ! array of unique atoms, normally all of them
    real(RK)   , intent(out) :: DIPL(:,:,:,:,:)   ! (N2E,N1E,    2*L2+1,2*L1+1,3)
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------

    real(RK), allocatable :: PIDIPL(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,3)
    integer(IK)           :: i

    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

!   print*, 'SHGI: shgi_dip_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2
    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)

    call setif(-1_i8k,.false.) ! zero all bits

    ! For some reason even conservative screening of
    ! integrals (MAXEXP==50) causes the MO overlap matrix
    ! to be non-positive definite -- dpptrf() returns INFO>0.
    ! So in normal integral run dont use any screening:
    call shgi_set_maxexp(HUGE(1.0_rk))
    ! OTOH, in "property" runs, eg. gradients, second derivatives
    ! it may be safer to apply integral screening.
    ! See shgi_gr_drv(), shgi_sd_drv(), etc ...

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    allocate(PIDIPL(NAB,2*LA+1,2*LB+1,3))

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(0,0,0)

    ! 2c overlap and kinetic integrals:
    call shgi_set_ovrl(LA,LB,0,0,0)

    ! actually compute the dipole matrix elements
    call shgi_dipole(PIDIPL)

    ! unpack dipole:
    do i=1,3
      call upack2c(PIDIPL(:,:,:,i),DIPL(:,:,:,:,i))
    enddo

!   call show('dip', DIPL(1,1,:,:,3))
    deallocate(PIDIPL)

    call shgi_close_ab()
    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()
    DPRINT  'SHGI: shgi_dip_drv: exit'

  end subroutine shgi_dip_drv

  subroutine shgi_dipole(dip)
    ! for all type dipole intergral, i.e. ll, ls and ss type
    use shgi_common, only: W, S5, NORM
    use shgi_shr, only: SHR_D2C, SHR_PR2v
    use shgi_common, only: WDAL, WDBL
    implicit none
    !------------ Declaration of formal parameters ------------------
    real(RK)   , intent(out) :: dip(:,:,:,:)   ! (NAB,2*LA+1,2*LB+1,3)
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------
    integer(IK), parameter  :: LW = 1
    real(RK)                :: D2C(NAB,(LW+1)**2,(LA+1)**2,(LB+1)**2)
    integer(IK)             :: lm
    integer(IK)             :: lma,ila
    integer(IK)             :: lmb,ilb
    integer(IK)             :: i,ma,mb
    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    ! compute the derivatives D(lma)D(lmb) * C(lmw), LW == 1:
    call SHR_D2C(NAB,LW,LA,LB,W,D2C)

    ! D(lma)D(lmb) -> DA(lma)DB(lmb):
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)
       if( ila + ilb > LW ) CYCLE ! high derivatives of vector are anyway zero

       do lm =1,(LW+1)**2
         D2C(:,lm,lma,lmb) = D2C(:,lm,lma,lmb) &
                           * WDAL(:,1+ila) * WDBL(:,1+ilb)
       enddo
    enddo
    enddo

    ! SHR_PR2v() increments the integrals, so initialize first:
    dip = 0.0
    ! switch from (zxy) to (xyz) here:
    call SHR_PR2v(NAB,LA,LB,D2C(:,CZ,:,:),S5(:,:,:,C0,C0,C0),dip(:,:,:,VZ))
    call SHR_PR2v(NAB,LA,LB,D2C(:,CX,:,:),S5(:,:,:,C0,C0,C0),dip(:,:,:,VX))
    call SHR_PR2v(NAB,LA,LB,D2C(:,CY,:,:),S5(:,:,:,C0,C0,C0),dip(:,:,:,VY))

    ! renorm the integrals (S5 was not normalized):
    do i=1,3
    do mb=1,2*LB+1
    do ma=1,2*LA+1
        dip(:,ma,mb,i) = dip(:,ma,mb,i) * NORM(:,2)
    enddo
    enddo
    enddo
  end subroutine shgi_dipole

  !****************************************************************
  !*************** DRIVER FOR GRADIENT OF DIPOLE INTEGRALS ********
  !****************************************************************

  subroutine shgi_gr_efield_drv(IU1,IE1,IL1,IU2,IE2,IL2,uas,densmat,&
                        & efield,VECUAGR)
    use unique_atom_module, only: uat=>unique_atom_type
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: upack2c, shgi_timing, shgi_gr_wd_to_ab
    use shgi_utils , only: shgi_ve_abc_store
    use shgi_common, only: CUTOFF
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)    :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)    :: uas(:)            ! array of unique atoms, normally all of them
    real(RK)   , intent(in)    :: densmat(:,:,:,:)  ! (NEA,NEB,2*LA+1,2*LB+1), section of the density matrix
    real(RK)   , intent(in)    :: efield(3) ! Electric field strength
    real(RK)   , intent(inout) :: VECUAGR(:)
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------

    real(RK), allocatable :: dip_gr(:,:,:,:)  ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: gr(:,:,:,:)      ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), allocatable :: PAB(:,:,:)       ! (NAB,2*LA+1,2*LB+1)   -- packed densmat
    integer(IK)           :: i,ma,mb
    real(RK)              :: VAB(6)           !                   (6) -- wrt A, B trace(PAB,dip_gr)
    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

!   print*, 'SHGI: shgi_dip_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2
    DPRINT 'SHGI: shgi_gr_efield_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2
    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)

    call setif(-1_i8k,.false.) ! zero all bits

    ! For some reason even conservative screening of
    ! integrals (MAXEXP==50) causes the MO overlap matrix
    ! to be non-positive definite -- dpptrf() returns INFO>0.
    ! So in normal integral run dont use any screening:
    call shgi_set_maxexp(HUGE(1.0_rk))
    ! OTOH, in "property" runs, eg. gradients, second derivatives
    ! it may be safer to apply integral screening.
    ! See shgi_gr_drv(), shgi_sd_drv(), etc ...

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    allocate(gr(NAB,2*LA+1,2*LB+1,6))
    allocate(dip_gr(NAB,2*LA+1,2*LB+1,6))

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(1,0,0)

    ! S5 overlap integrals:
    call shgi_set_ovrl(LA,LB,1,0,0)

    allocate(PAB(NAB,2*LA+1,2*LB+1))

    ! pack density matrix using CUTOFF that was set in shgi_set_ab():
    do mb=1,2*LB+1
    do ma=1,2*LA+1
      PAB(:,ma,mb) = pack(densmat(:,:,ma,mb),CUTOFF)
    enddo
    enddo

    ! actually compute the gradient of dipole integral elements
    call shgi_gr_efield(efield,gr)

    ! transformation of gradient from W and D to A and B
    call shgi_gr_wd_to_ab(gr,dip_gr,0)

!   ASSERT(NAB == NEA * NEB)
!   call show ('dipder', dip_gr(1,:,:,6)/0.002)

    ! compute trace(PAB,dip_gr):
    do i=1,6
      VAB(i) = sum( dip_gr(:,:,:,i) * PAB(:,:,:) )
    enddo

    ! for symmetry adaptation
    call shgi_ve_abc_store(IU2,IE2,IU1,IE1,0,0,VAB,VECUAGR,+1)
    !
    ! Now I have confirm the sign of electronic contribution to the total gradients
    ! in presence of electric field.
    !

    deallocate(PAB,dip_gr,gr)


    call shgi_close_ab()
    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()
    DPRINT  'SHGI: shgi_gr_efield_drv: exit'

  end subroutine shgi_gr_efield_drv

  subroutine shgi_gr_efield(efield,gr)
    !
    ! it will return 6 gradient wrt W and D
    ! for gradient of all type dipole intergral, i.e. ll, ls and ss type
    !
    use shgi_common, only: W, S5, NORM
    use shgi_shr, only: SHR_D3C
    use shgi_common, only: WDAL, WDBL
    use shgi_utils, only: grFS
    implicit none
    !------------ Declaration of formal parameters ------------------
    real(RK), intent(in)  :: efield(3) ! Electric field strength
    real(RK), intent(inout) :: gr(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------
    integer(IK), parameter  :: LW = 1
    integer(IK)             :: lm
    integer(IK)             :: lma,ila
    integer(IK)             :: lmb,ilb
    integer(IK)             :: lmc
    integer(IK)             :: i,ma,mb
    real(RK)                :: D3C(NAB,(LW+1)**2,(LA+1)**2,(LB+1)**2,(LC+1)**2)
    real(RK)                :: scal(NAB,         (LA+1)**2,(LB+1)**2,(LC+1)**2)
    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    ! compute the derivatives D(lma)D(lmb)D(lmc) * C(lmw), LW == 1:
    call SHR_D3C(NAB,LW,LA,LB,LC,W,D3C)

    ! D(lma)D(lmb)D(lmc) -> DA(lma)DB(lmb)DC(lmc):
    do lmc=1,(LC+1)**2
      do lmb=1,(LB+1)**2
         ilb=lof(lmb)
        do lma=1,(LA+1)**2
           ila=lof(lma)
           if( ila + ilb > LW ) CYCLE ! high derivatives of vector are anyway zero

           do lm =1,(LW+1)**2
               D3C(:,lm,lma,lmb,lmc) = D3C(:,lm,lma,lmb,lmc) &
                                 * WDAL(:,1+ila) * WDBL(:,1+ilb)
           enddo
        enddo
      enddo
    enddo

    ! Coupling of Electric field with W
    scal(:,:,:,:) = efield(VX)*D3C(:,CX,:,:,:) &
                  + efield(VY)*D3C(:,CY,:,:,:) &
                  + efield(VZ)*D3C(:,CZ,:,:,:)

    ! Taking the gradient of the product of two function F S (including the product rule)
    gr = 0.0
    call grFS(scal,S5(:,:,:,:,C0,C0),gr)

    ! renorm the integrals (S5 was not normalized):
    do i=1,6
    do mb=1,2*LB+1
    do ma=1,2*LA+1
        gr(:,ma,mb,i) = gr(:,ma,mb,i) * NORM(:,2)
    enddo
    enddo
    enddo
  end subroutine shgi_gr_efield


  !--------------- End of module -------------------------------------
end module shgi_dip

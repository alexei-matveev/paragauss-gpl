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
module shgi_ep_ef
  !---------------------------------------------------------------
  !
  ! Maybe  used for  many point  charges  in PCM  and other  embedding
  ! models.
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
  use shgi_cntrl
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: shgi_pot_drv
  public :: shgi_field_drv

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

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !****************************************************************
  !****************** DRIVER FOR INTEGRALS ************************
  !****************************************************************

  subroutine shgi_pot_drv(IU1,IE1,IL1,IU2,IE2,IL2,uas,spp,POTEN)
    use unique_atom_module, only: uat=>unique_atom_type
    use potential_module, only: spt=>spacepoint_type, N_points
    use options_module, only: options_integral_expmax
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: upack2c, shgi_timing
    use shgi_pcm, only: shgi_pc
    implicit none
    type(uat)  , intent(in)  :: uas(:)            ! array of unique atoms, normally all of them
    type(spt), intent(in)    :: spp(:)            ! array of unique space points
    integer(IK), intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    real(RK), intent(out)    :: POTEN (:,:,:,:,:)
    ! *** end of interface ***

    integer(IK) :: up,ep,n_equal
    real(RK)    :: x(3)
    real(RK), parameter:: z=1.0_RK

    ! temp storage for electrostatic potential
    real(RK), allocatable :: PIEPOT(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    integer(IK) :: memstat

    DPRINT 'SHGI: shgi_pot_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)
    FPP_TIMER_START(pcs)
    FPP_TIMER_START(epot)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! NOT USED: ! allocate global storage of THIS module
    ! NOT USED: call shgi_glob_alloc(NAB,LA,LB)

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(0,0,0)

    ! compute overlap S5:
    call shgi_set_ovrl(LA,LB,0,0,0) !(2) shgi_pot_drv

    allocate(PIEPOT(NAB,2*LA+1,2*LB+1),stat=memstat)
    ASSERT(memstat==0)

    do up=1,N_points
       PIEPOT = 0.0_rk
       n_equal=spp(up)%n_equal_points
       do ep=1,n_equal
          x=spp(up)%position(:,ep)
          call shgi_pc(z,x,PIEPOT)
       end do
       call upack2c(PIEPOT,POTEN(:,:,up,:,:))
    end do

    deallocate(PIEPOT,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()
    ! NOT USED: ! deallocate global storage of THIS module:
    ! NOT USED: call shgi_glob_free()

    FPP_TIMER_STOP(epot)
    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_pot_drv

  subroutine shgi_field_drv(IU1,IE1,IL1,IU2,IE2,IL2,uas,spp,EFIELD)
    use unique_atom_module, only: uat=>unique_atom_type
    use elec_static_field_module, only: spt=>outward_normal, N_surface_points
    use options_module, only: options_integral_expmax
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_utils , only: shgi_fl_wd_store, shgi_timing
    use shgi_pcm, only: shgi_gr_pc
    implicit none
    type(uat)  , intent(in)  :: uas(:)            ! array of unique atoms, normally all of them
    type(spt), intent(in)    :: spp(:)            ! array of unique space points
    integer(IK), intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    real(RK), intent(out)    :: EFIELD (:,:,:,:,:)
    ! *** end of interface ***

    integer(IK) :: up,ep,n_equal
    real(RK)    :: x(3)
    real(RK), parameter:: z=1.0_RK

    ! temp storage for electrostatic potential
    real(RK), allocatable :: PIEFLD(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    integer(IK) :: memstat

    DPRINT 'SHGI: shgi_field_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2

    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)
    FPP_TIMER_START(pcs)
    FPP_TIMER_START(efld)

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
    call shgi_set_ovrl(LA,LB,1,0,0) !(2) shgi_pot_drv

    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    allocate(PIEFLD(NAB,2*LA+1,2*LB+1,6),stat=memstat)
    ASSERT(memstat==0)

    do up=1,N_surface_points
       n_equal=spp(up)%n_equal_points
       do ep=1,n_equal
          PIEFLD = 0.0_rk
          x=spp(up)%position(:,ep)
          call shgi_gr_pc(z,x,PIEFLD)
          call shgi_fl_wd_store(up,ep,PIEFLD,EFIELD,1)
       end do
    end do

    deallocate(PIEFLD,stat=memstat)
    ASSERT(memstat==0)

    call shgi_close_ab()

    FPP_TIMER_STOP(efld)
    FPP_TIMER_STOP(pcs)
    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    call shgi_timing()

  end subroutine shgi_field_drv
  !--------------- End of module ----------------------------------
end module shgi_ep_ef

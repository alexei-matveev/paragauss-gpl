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
!===============================================================
! Public interface of module
!===============================================================
module shgi_relfit
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  !
  !
  !  Author: ...
  !  Date: ...
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
! use CPU_TIME for timers:
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
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
  public :: shgi_nuc
  public :: shgi_nuc_so
  public :: shgi_nuc_sr
  public :: shgi_gr_nuc
  public :: shgi_gr_nuc_sr
  public :: shgi_ch
  public :: shgi_ch_so
  public :: shgi_ch_sr
  public :: shgi_r2
  public :: shgi_r2_so
  public :: shgi_r2_sr

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

  !------------ Subroutines ---------------------------------------
contains

  subroutine shgi_nuc(z,nuc,RAD)
    use shgi_common, only: W2, NORM, LAMBDA, YS
    use shgi_rad, only: doIL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: nuc(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(in)    :: RAD
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R

    if (Z == ZERO) then
       return
    endif

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_nuc: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      real(RK) :: IL(NAB,1+LA+LB)

      call doIL(LA+LB,W2,LAM, NRM, IL)
      call doRadAng(LA+LB,IL,YS(:,:,:,:,C00),nuc)
    end subroutine point

    subroutine finit()
      implicit none
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_nuc

  subroutine shgi_gr_nuc(z,nuc,RAD)
    ! Returns W,D gradients of dimension 1:6 !
    use shgi_common, only: NORM, LAMBDA
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: nuc(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: RAD
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R

    if (Z == ZERO) RETURN

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_gr_nuc: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: W2, GS
      use shgi_rad, only: doIL
      use shgi_utils, only: doRadAng
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      real(RK) :: IL(NAB,1+LA+LB+1)

      call doIL(LA+LB+1,W2,LAM, NRM, IL)

      ! GW:
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GWX),nuc(:,:,:,GWX))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GWY),nuc(:,:,:,GWY))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GWZ),nuc(:,:,:,GWZ))

      ! GD:
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GDX),nuc(:,:,:,GDX))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GDY),nuc(:,:,:,GDY))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,GDZ),nuc(:,:,:,GDZ))
    end subroutine point

    subroutine finit()
      implicit none
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_gr_nuc

  subroutine shgi_nuc_so(z,so,RAD)
    use shgi_common, only: W2, NORM, LAMBDA, YS
    use shgi_rad, only: doIL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: so(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(in)    :: RAD
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R

    if (Z == ZERO) then
       return
    endif

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_nuc_so: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      real(RK) :: IL(NAB,2+LA+LB) ! one more

      call doIL(LA+LB+1,W2,LAM, NRM, IL)
      call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVX),so(:,:,:,VX))
      call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVY),so(:,:,:,VY))
      call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVZ),so(:,:,:,VZ))
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_nuc_so

  subroutine shgi_nuc_sr(z,sr,RAD)
    use shgi_common, only: W2, NORM, LAMBDA, YS, WDA, WDB
    use shgi_rad, only: doIL, doP2IL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: sr(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(in)    :: RAD
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R

    if (Z == ZERO) then
       return
    endif

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_nuc_sr: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      real(RK) :: IL(NAB,3+LA+LB) ! two more

      call doIL(    LA+LB+2, W2, LAM, NRM, IL)
      call doRadAng(LA+LB  , IL, YS(:,:,:,:,CP2), sr(:,:,:))
      call doRadAng(LA+LB+1, IL, YS(:,:,:,:,CR2), sr(:,:,:), fact=WDB-WDA)
      call doP2IL(  LA+LB  , W2, IL)
      call doRadAng(LA+LB  , IL, YS(:,:,:,:,C00), sr(:,:,:), fact=WDB*WDA)
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_nuc_sr

  subroutine shgi_gr_nuc_sr(z,nuc,RAD)
    ! Returns W,D grads of dimension 1:6 !
    use shgi_common, only: NORM, LAMBDA
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: nuc(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: RAD
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R

    if (Z == ZERO) RETURN

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_gr_nuc_sr: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: W2, GS, WDA, WDB
      use shgi_rad, only: doIL, doP2IL
      use shgi_utils, only: doRadAng
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      real(RK) :: IL(NAB,1+LA+LB+3)

      call doIL(    LA+LB+3,W2,LAM, NRM, IL)

      ! F(0)*Sp2(0)
      ! GW:
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,12+GWX),nuc(:,:,:,GWX))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,12+GWY),nuc(:,:,:,GWY))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,12+GWZ),nuc(:,:,:,GWZ))
      ! GD:
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,12+GDX),nuc(:,:,:,GDX))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,12+GDY),nuc(:,:,:,GDY))
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:,12+GDZ),nuc(:,:,:,GDZ))

      ! F(x)*S(x) + F(y)*S(y) + F(y)*S(y):
      ! GW:
      call doRadAng(LA+LB+2,IL,GS(:,:,:,:, 6+GWX),nuc(:,:,:,GWX), fact=WDB-WDA)
      call doRadAng(LA+LB+2,IL,GS(:,:,:,:, 6+GWY),nuc(:,:,:,GWY), fact=WDB-WDA)
      call doRadAng(LA+LB+2,IL,GS(:,:,:,:, 6+GWZ),nuc(:,:,:,GWZ), fact=WDB-WDA)
      ! GD:
      call doRadAng(LA+LB+2,IL,GS(:,:,:,:, 6+GDX),nuc(:,:,:,GDX), fact=WDB-WDA)
      call doRadAng(LA+LB+2,IL,GS(:,:,:,:, 6+GDY),nuc(:,:,:,GDY), fact=WDB-WDA)
      call doRadAng(LA+LB+2,IL,GS(:,:,:,:, 6+GDZ),nuc(:,:,:,GDZ), fact=WDB-WDA)

      ! Fp2(0)*S(0):
      call doP2IL(  LA+LB+1,W2, IL)

      ! GW:
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:, 0+GWX),nuc(:,:,:,GWX), fact=WDB*WDA)
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:, 0+GWY),nuc(:,:,:,GWY), fact=WDB*WDA)
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:, 0+GWZ),nuc(:,:,:,GWZ), fact=WDB*WDA)
      ! GD:
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:, 0+GDX),nuc(:,:,:,GDX), fact=WDB*WDA)
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:, 0+GDY),nuc(:,:,:,GDY), fact=WDB*WDA)
      call doRadAng(LA+LB+1,IL,GS(:,:,:,:, 0+GDZ),nuc(:,:,:,GDZ), fact=WDB*WDA)
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_gr_nuc_sr

  subroutine shgi_ch(ch)
    use shgi_common, only: W2, NORMCH, LAMBDACH, YS
    use shgi_rad, only: doIL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(inout) :: ch(:,:,:,:) ! (NAB,NEC,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,1+LA+LB)
    integer(IK) :: ic

    do ic=1,size(LAMBDACH,2)
       call doIL(LA+LB,W2,LAMBDACH(:,ic), NORMCH(:,ic), IL)
       call doRadAng(LA+LB,IL,YS(:,:,:,:,C00),ch(:,ic,:,:))
    enddo
  end subroutine shgi_ch

  subroutine shgi_ch_so(so)
    use shgi_common, only: W2, NORMCH, LAMBDACH, YS
    use shgi_rad, only: doIL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(inout) :: so(:,:,:,:,:) ! (NAB,NEC,2*LA+1,2*LB+1,3)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,2+LA+LB) ! one more
    integer(IK) :: ic

    do ic=1,size(LAMBDACH,2)
       call doIL(LA+LB+1,W2,LAMBDACH(:,ic), NORMCH(:,ic), IL)
       call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVX),so(:,ic,:,:,VX))
       call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVY),so(:,ic,:,:,VY))
       call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVZ),so(:,ic,:,:,VZ))
    enddo
  end subroutine shgi_ch_so

  subroutine shgi_ch_sr(sr)
    use shgi_common, only: W2, NORMCH, LAMBDACH, YS, WDA, WDB
    use shgi_rad, only: doIL, doP2IL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(inout) :: sr(:,:,:,:) ! (NAB,NEC,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,3+LA+LB) ! two more
    integer(IK) :: ic

    do ic=1,size(LAMBDACH,2)
       call doIL(    LA+LB+2, W2, LAMBDACH(:,ic), NORMCH(:,ic), IL)
       call doRadAng(LA+LB  , IL, YS(:,:,:,:,CP2), sr(:,ic,:,:))
       call doRadAng(LA+LB+1, IL, YS(:,:,:,:,CR2), sr(:,ic,:,:), fact=WDB-WDA)
       call doP2IL(  LA+LB  , W2, IL)
       call doRadAng(LA+LB  , IL, YS(:,:,:,:,C00), sr(:,ic,:,:), fact=WDB*WDA)
    enddo
  end subroutine shgi_ch_sr

  subroutine shgi_r2(ch)
    use shgi_common, only: W2, NORMR2, LAMBDAR2, F132R2, YS
    use shgi_rad, only: doRL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(inout) :: ch(:,:,:,:) ! (NAB,NEC,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,1+LA+LB)
    integer(IK) :: ic

    do ic=1,size(LAMBDAR2,2)
       call doRL(LA+LB,W2,LAMBDAR2(:,ic), F132R2(:,ic), NORMR2(:,ic), IL)
       call doRadAng(LA+LB,IL,YS(:,:,:,:,C00),ch(:,ic,:,:))
    enddo
  end subroutine shgi_r2

  subroutine shgi_r2_so(so)
    use shgi_common, only: W2, NORMR2, LAMBDAR2, F132R2, YS
    use shgi_rad, only: doRL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(inout) :: so(:,:,:,:,:) ! (NAB,NEC,2*LA+1,2*LB+1,3)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,2+LA+LB) ! one more
    integer(IK) :: ic

    do ic=1,size(LAMBDAR2,2)
       call doRL(LA+LB+1,W2,LAMBDAR2(:,ic), F132R2(:,ic), NORMR2(:,ic), IL)
       call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVX),so(:,ic,:,:,VX))
       call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVY),so(:,ic,:,:,VY))
       call doRadAng(LA+LB+1,IL,YS(:,:,:,:,CVZ),so(:,ic,:,:,VZ))
    enddo
  end subroutine shgi_r2_so

  subroutine shgi_r2_sr(sr)
    use shgi_common, only: W2, NORMR2, LAMBDAR2, F132R2, YS, WDA, WDB
    use shgi_rad, only: doRL, doP2IL
    use shgi_utils , only: doRadAng
    implicit none
    real(RK)   , intent(inout) :: sr(:,:,:,:) ! (NAB,NEC,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK)    :: IL(NAB,3+LA+LB) ! two more
    integer(IK) :: ic

    do ic=1,size(LAMBDAR2,2)
       call doRL(    LA+LB+2, W2, LAMBDAR2(:,ic), F132R2(:,ic), NORMR2(:,ic), IL)
       call doRadAng(LA+LB  , IL, YS(:,:,:,:,CP2), sr(:,ic,:,:))
       call doRadAng(LA+LB+1, IL, YS(:,:,:,:,CR2), sr(:,ic,:,:), fact=WDB-WDA)
       call doP2IL(  LA+LB  , W2, IL)
       call doRadAng(LA+LB  , IL, YS(:,:,:,:,C00), sr(:,ic,:,:), fact=WDB*WDA)
    enddo
  end subroutine shgi_r2_sr

  !--------------- End of module ----------------------------------
end module shgi_relfit

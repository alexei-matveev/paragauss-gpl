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
module shgi_relnuc
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
  public :: shgi_rel_nuc   !(z,x,nu,sr,so,RAD) -- does not need pre-computed ang part
  public :: shgi_rel_gr_nuc!(z,x,nu,sr,   RAD)
  public :: shgi_rel_sd_nuc!(z,x,nu,sr,   RAD)

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

  subroutine shgi_rel_nuc(z,x,nu,sr,so,RAD)
    use shgi_common, only: W, NORM, LAMBDA
    implicit none
    real(RK)   , intent(in)    :: z, x(3)
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(inout) :: sr(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(inout) :: so(:,:,:,:)!(NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(in)    :: RAD
    optional :: nu
    optional :: sr
    optional :: so
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R
    logical  :: nucl
    logical  :: srel
    logical  :: spor
    real(RK) :: wc(NAB,3), w2(NAB)

    if (Z == ZERO) RETURN

    nucl = present(nu)
    srel = present(sr)
    spor = present(so)

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: S5, K4, WDA, WDB
      use shgi_rad, only: doIL, doP2IL
      use shgi_shr, only: SHR_D3Fv, SHR_PR2v
      use shgi_utils, only: shgi_dww_to_dab
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      integer(IK), parameter :: L1=1
      integer(IK) :: LX, ND
      real(RK)    :: IL(NAB,3+LA+LB) ! two more
      real(RK)    :: F3(NAB,(LA+1)**2,(LB+1)**2,4) ! (L1+1)**2 == 4

      LX = 0
      ND = 0
      if( srel .or. spor )then
        LX = 1 ! one more nabla for SR
        ND = 2 ! two more orders for SR
      endif

      ! radial ints, higher orders for SR/SO:
      call doIL(LA+LB+ND,w2,LAM,NRM,IL)

      ! derivatives D(lma)D(lmb)D(lmc) * F:
      call SHR_D3Fv(NAB,LA,LB,LX,wc,IL,F3)

      ! D(lma)D(lmb)D(lmc)F -> DA(lma)DB(lmb)D(lmc)F:
      call shgi_dww_to_dab(NAB,LA,LB,(LX+1)**2,F3)

      ! NUCLEAR ATTRACTION:
      if( nucl )then
        FPP_TIMER_START(tpr2)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,C0),S5(:,:,:,C0,C0,C0),nu)
        FPP_TIMER_STOP(tpr2)
      endif

      ! SCALAR RELATVISTIC MAT EL:
      if( srel )then
        FPP_TIMER_START(tpr2)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,C0),K4(:,:,:,C0,C0   ),sr)

        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CX),S5(:,:,:,CX,C0,C0),sr, fact=WDB-WDA)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CY),S5(:,:,:,CY,C0,C0),sr, fact=WDB-WDA)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CZ),S5(:,:,:,CZ,C0,C0),sr, fact=WDB-WDA)
        FPP_TIMER_STOP(tpr2)

        ! radial ints for laplace:
        call doP2IL(LA+LB,w2,IL)

        ! derivatives D(lma)D(lmb) * Fp2:
        call SHR_D3Fv(NAB,LA,LB,0,wc,IL,F3)

        ! D(lma)D(lmb)F -> DA(lma)DB(lmb)F:
        call shgi_dww_to_dab(NAB,LA,LB,1,F3)

        FPP_TIMER_START(tpr2)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,C0),S5(:,:,:,C0,C0,C0),sr, fact=WDB*WDA)
        FPP_TIMER_STOP(tpr2)
      endif

      ! SPIN-ORBIT MAT EL:
      if( spor )then
        FPP_TIMER_START(tpr2)
        ! x = yz-zy:
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CY),S5(:,:,:,CZ,C0,C0),so(:,:,:,VX), fact= WDB+WDA)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CZ),S5(:,:,:,CY,C0,C0),so(:,:,:,VX), fact=-WDB-WDA)
        ! y = zx-xz:
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CZ),S5(:,:,:,CX,C0,C0),so(:,:,:,VY), fact= WDB+WDA)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CX),S5(:,:,:,CZ,C0,C0),so(:,:,:,VY), fact=-WDB-WDA)
        ! z = xy-yx:
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CX),S5(:,:,:,CY,C0,C0),so(:,:,:,VZ), fact= WDB+WDA)
        call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CY),S5(:,:,:,CX,C0,C0),so(:,:,:,VZ), fact=-WDB-WDA)
        FPP_TIMER_STOP(tpr2)
      endif
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit
  end subroutine shgi_rel_nuc

  subroutine shgi_rel_gr_nuc(z,x,nu,sr,RAD)
    ! Returns W,D grads of dimension 1:6 !
    use shgi_common, only: W, NORM, LAMBDA
    implicit none
    real(RK)   , intent(in)    :: z, x(3)
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(inout) :: sr(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: RAD
    optional :: nu
    optional :: sr
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R
    logical  :: nucl
    logical  :: srel
#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
    real(RK) :: wc(NAB,3), w2(NAB)
#else
    real(RK), allocatable :: wc(:,:), w2(:)

    allocate(wc(NAB,3), w2(NAB))
#endif


    if (Z == ZERO) RETURN

    nucl = present(nu)
    srel = present(sr)

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
      deallocate(wc,w2)
#endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: S5, K4, WDA, WDB
      use shgi_rad, only: doIL, doP2IL
      use shgi_shr, only: SHR_D4Fv
      use shgi_utils, only: grFS
      use shgi_utils, only: shgi_dww_to_dab
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      integer(IK), parameter :: L1=1
      integer(IK) :: LX, ND
#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
      real(RK)    :: IL(NAB,1+LA+LB+3)
      real(RK)    :: F4(NAB,(LA+1)**2,(LB+1)**2,4,4) ! (L1+1)**2 == 4
#else
      real(RK), allocatable :: IL(:,:)
      real(RK), allocatable :: F4(:,:,:,:,:)

      allocate(IL(NAB,1+LA+LB+3),F4(NAB,(LA+1)**2,(LB+1)**2,4,4))
#endif

      LX = 0
      ND = 0
      if( srel )then
        LX = 1 ! one more nabla for SR
        ND = 2 ! two more orders for SR
      endif

      call doIL(LA+LB+L1+ND,w2,LAM, NRM, IL)

      ! derivatives D(lma)D(lmb)D(lm1)D(lmx) * F:
      call SHR_D4Fv(NAB,LA,LB,L1,LX,wc,IL,F4)

      ! D(lma)D(lmb)F -> DA(lma)DB(lmb)F:
      call shgi_dww_to_dab(NAB,LA,LB,4*(LX+1)**2,F4) ! 4 == (L1+1)**2

      ! NUCLEAR ATTRACTION:
      if( nucl )then
        ! GR(cg) : -(wa+wb)*F(g)*S(0)             , wa+wb=1, g=x,y,z
        ! GR(ag) :   wa    *F(g)*S(0) + F(0)*S(g)
        ! GR(bg) :      wb *F(g)*S(0) - F(0)*S(g)
        !                   \_ GW _/    \_ GD _/

        ! F(g)*S(0), F(0)*S(g):
        call grFS(F4(:,:,:,:,C0),S5(:,:,:,:,C0,C0),nu)
      endif

      ! SCALAR RELATIVISTIC MAT EL:
      if( srel )then
        ! F(g)*Sp2(0), F(0)*Sp2(g):
        call grFS(F4(:,:,:,:,C0),K4(:,:,:,:,C0),sr)

        ! F(x)*S(x) + F(y)*S(y) + F(y)*S(y):
        !   F(g)*S(0) -SR-> F(g,x)*S(0,x) + F(g,y)*S(0,y) + F(g,z)*S(0,z),
        !   F(0)*S(g) -SR-> F(0,x)*S(g,x) + F(0,y)*S(g,y) + F(0,z)*S(g,z):
        call grFS(F4(:,:,:,:,CX),S5(:,:,:,:,CX,C0),sr, fact=WDB-WDA)
        call grFS(F4(:,:,:,:,CY),S5(:,:,:,:,CY,C0),sr, fact=WDB-WDA)
        call grFS(F4(:,:,:,:,CZ),S5(:,:,:,:,CZ,C0),sr, fact=WDB-WDA)

        ! Fp2:
        call doP2IL(LA+LB+L1,w2,IL)

        ! derivatives D(lma)D(lmb)D(lm1)D(lmx) * Fp2:
        call SHR_D4Fv(NAB,LA,LB,L1,0,wc,IL,F4)

        ! D(lma)D(lmb)F -> DA(lma)DB(lmb)F:
        call shgi_dww_to_dab(NAB,LA,LB,4,F4) ! 4 == (L1+1)**2

        ! Fp2(g)*S(0), Fp2(0)*S(g):
        call grFS(F4(:,:,:,:,C0),S5(:,:,:,:,C0,C0),sr, fact=WDB*WDA)
      endif
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
      deallocate(IL,F4)
#endif
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit
  end subroutine shgi_rel_gr_nuc

  subroutine shgi_rel_sd_nuc(z,x,nu,sr,RAD)
    ! Returns W,D derivatives of dimension 6x6 !
    use shgi_common, only: W, NORM, LAMBDA
    implicit none
    real(RK)   , intent(in)    :: z, x(3)
    real(RK)   , intent(inout) :: nu(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK)   , intent(inout) :: sr(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK)   , intent(in)    :: RAD
    optional :: nu
    optional :: sr
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R
    logical  :: nucl
    logical  :: srel
    real(RK) :: wc(NAB,3), w2(NAB)

    if (Z == ZERO) RETURN

    nucl = present(nu)
    srel = present(sr)

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    wc(:,1) = W(:,1) - x(1)
    wc(:,2) = W(:,2) - x(2)
    wc(:,3) = W(:,3) - x(3)
    w2(:)   = wc(:,1)**2 + wc(:,2)**2 + wc(:,3)**2

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: S5, K4, WDA, WDB
      use shgi_rad, only: doIL, doP2IL
      use shgi_shr, only: SHR_D5Fv
      use shgi_utils, only: sdFS, sdSYM
      use shgi_utils, only: shgi_dww_to_dab
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      integer(IK), parameter :: L1=1
      integer(IK) :: LX, ND
#ifndef FPP_NO_BIG_AUTOMATIC_ARRAYS
      real(RK)    :: IL(NAB,1+LA+LB+4)
      real(RK)    :: F5(NAB,(LA+1)**2,(LB+1)**2,4,4,4) ! (L1+1)**2 == 4
#else
      real(RK), allocatable :: IL(:,:)
      real(RK), allocatable :: F5(:,:,:,:,:,:)

      allocate(IL(NAB,1+LA+LB+4),F5(NAB,(LA+1)**2,(LB+1)**2,4,4,4))
#endif

      LX = 0
      ND = 0
      if( srel )then
        LX = 1 ! one more nabla for SR
        ND = 2 ! two more orders for SR
      endif

      call doIL(LA+LB+L1+L1+ND,w2,LAM, NRM, IL)

      ! derivatives D(lma)D(lmb)D(lm1)D(lm1)D(lmx) * F:
      call SHR_D5Fv(NAB,LA,LB,L1,L1,LX,wc,IL,F5)

      ! D(lma)D(lmb)F -> DA(lma)DB(lmb)F:
      call shgi_dww_to_dab(NAB,LA,LB,4*4*(LX+1)**2,F5) ! 4 == (L1+1)**2

      ! NUCLEAR ATTRACTION:
      if( nucl )then
        ! F(i,j)*S(0), F(i,0)*S(0,j), F(0)*S(i,j):
        call sdFS(F5(:,:,:,:,:,C0),S5(:,:,:,:,:,C0),nu)

        call sdSYM(nu)
      endif

      ! SCALAR RELATIVISTIC MAT EL:
      if( srel )then
        ! F*Sp2:
        call sdFS(F5(:,:,:,:,:,C0),K4(:,:,:,:,:),sr)

        ! F(x)*S(x) + F(y)*S(y) + F(y)*S(y):
        call sdFS(F5(:,:,:,:,:,CX),S5(:,:,:,:,:,CX),sr, fact=WDB-WDA)
        call sdFS(F5(:,:,:,:,:,CY),S5(:,:,:,:,:,CY),sr, fact=WDB-WDA)
        call sdFS(F5(:,:,:,:,:,CZ),S5(:,:,:,:,:,CZ),sr, fact=WDB-WDA)

        ! Fp2:
        call doP2IL(LA+LB+L1+L1,W2,IL)

        ! D4 Fp2:
        ! derivatives D(lma)D(lmb)D(lm1)D(lm1) * Fp2:
        call SHR_D5Fv(NAB,LA,LB,L1,L1,0,wc,IL,F5)

        ! D(lma)D(lmb)F -> DA(lma)DB(lmb)F:
        call shgi_dww_to_dab(NAB,LA,LB,4*4,F5) ! 4 == (L1+1)**2

        ! Fp2*S:
        call sdFS(F5(:,:,:,:,:,C0),S5(:,:,:,:,:,C0),sr, fact=WDB*WDA)

        call sdSYM(sr)
      endif
#ifdef FPP_NO_BIG_AUTOMATIC_ARRAYS
      deallocate(IL,F5)
#endif
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit
  end subroutine shgi_rel_sd_nuc

#if 0
  subroutine shgi_rel_nuc(z,nu,sr,so,RAD)
    use shgi_common, only: W2, NORM, LAMBDA
    use shgi_rad, only: doIL, doP2IL
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: nu(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(inout) :: sr(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(inout) :: so(:,:,:,:)!(NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(in)    :: RAD
    optional :: nu
    optional :: sr
    optional :: so
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R
    logical  :: nucl
    logical  :: srel
    logical  :: spor

    if (Z == ZERO) RETURN

    nucl = present(nu)
    srel = present(sr)
    spor = present(so)

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_rel_nuc: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: YL, S5, K4, WDA, WDB
      use shgi_shr, only: SHR_PR2v
      use shgi_dnf, only: doD3F
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      integer(IK) :: LX=1, ND=2
      real(RK)    :: IL(NAB,3+LA+LB) ! two more
      real(RK)    :: F3(NAB,(LA+1)**2,(LB+1)**2,4) ! 4==(P+1)**2

      LX = 1 ! one more nabla for SR
      ND = 2 ! two more orders for SR
      if( .not. srel )then
         LX = 0
         ND = 0
      endif

      call doIL(LA+LB+ND,W2,LAM,NRM,IL)
      ASSERT(LC>=LX)
      call doD3F(LA,LB,LX,IL,YL(:,:,:,:,C0,C0,:),F3)

      ! *** NUCLEAR ATTRACTION ***
      if( .not. nucl ) goto 111
      FPP_TIMER_START(tpr2)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,C0),S5(:,:,:,C0,C0,C0),nu)
      FPP_TIMER_STOP(tpr2)

111   if( .not. srel ) goto 222

      ! *** SCALAR RELATVISTIC MAT EL ***
      FPP_TIMER_START(tpr2)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,C0),K4(:,:,:,C0,C0   ),sr)

      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CX),S5(:,:,:,CX,C0,C0),sr, fact=WDB-WDA)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CY),S5(:,:,:,CY,C0,C0),sr, fact=WDB-WDA)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CZ),S5(:,:,:,CZ,C0,C0),sr, fact=WDB-WDA)
      FPP_TIMER_STOP(tpr2)

      call doP2IL(LA+LB,W2,IL)
      call doD3F(LA,LB,0,IL,YL(:,:,:,:,C0,C0,:),F3)
      FPP_TIMER_START(tpr2)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,C0),S5(:,:,:,C0,C0,C0),sr, fact=WDB*WDA)
      FPP_TIMER_STOP(tpr2)

      ! *** SPIN-ORBIT MAT EL ***
222   if( .not. spor ) RETURN

      FPP_TIMER_START(tpr2)
      ! x = yz-zy:
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CY),S5(:,:,:,CZ,C0,C0),so(:,:,:,VX), fact= WDB+WDA)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CZ),S5(:,:,:,CY,C0,C0),so(:,:,:,VX), fact=-WDB-WDA)
      ! y = zx-xz:
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CZ),S5(:,:,:,CX,C0,C0),so(:,:,:,VY), fact= WDB+WDA)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CX),S5(:,:,:,CZ,C0,C0),so(:,:,:,VY), fact=-WDB-WDA)
      ! z = xy-yx:
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CX),S5(:,:,:,CY,C0,C0),so(:,:,:,VZ), fact= WDB+WDA)
      call SHR_PR2v(NAB,LA,LB,F3(:,:,:,CY),S5(:,:,:,CX,C0,C0),so(:,:,:,VZ), fact=-WDB-WDA)
      FPP_TIMER_STOP(tpr2)
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_rel_nuc

  subroutine shgi_rel_gr_nuc(z,nu,sr,RAD)
    ! Returns W,D grads of dimension 1:6 !
    use shgi_common, only: W2, NORM, LAMBDA
    use shgi_rad, only: doIL, doP2IL
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: nu(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(inout) :: sr(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(in)    :: RAD
    optional :: nu
    optional :: sr
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R
    logical  :: nucl
    logical  :: srel

    if (Z == ZERO) RETURN

    nucl = present(nu)
    srel = present(sr)

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    DPRINT 'SHGI: shgi_rel_gr_nuc: R=',R

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: YL, S5, K4, WDA, WDB
      use shgi_dnf, only: doD4F
      use shgi_utils, only: grFS
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      integer(IK) :: LX=1, ND=2
      real(RK)    :: IL(NAB,1+LA+LB+3)
      real(RK)    :: F4(NAB,(LA+1)**2,(LB+1)**2,4,4) ! 4==(P+1)**2

      LX = 1 ! one more nabla for SR
      ND = 2 ! two more orders for SR
      if( .not. srel )then
         LX = 0
         ND = 0
      endif

      call doIL(LA+LB+ND+1,W2,LAM, NRM, IL)
      ASSERT(LC>= 1)
      ASSERT(LD>=LX)
      call doD4F(LA,LB,1,LX,IL,YL(:,:,:,:,:,C0,:),F4)

      if( .not. nucl ) goto 111
      ! *** NUCLEAR ATTRACTION ***
      ! GR(cg) : -(wa+wb)*F(g)*S(0)             , wa+wb=1, g=x,y,z
      ! GR(ag) :   wa    *F(g)*S(0) + F(0)*S(g)
      ! GR(bg) :      wb *F(g)*S(0) - F(0)*S(g)
      !                   \_ GW _/    \_ GD _/

      ! F(g)*S(0), F(0)*S(g):
      call grFS(F4(:,:,:,:,C0),S5(:,:,:,:,C0,C0),nu)

111   if( .not. srel ) RETURN
      ! *** SCALAR RELATIVISTIC MAT EL ***

      ! F(g)*Sp2(0), F(0)*Sp2(g):
      call grFS(F4(:,:,:,:,C0),K4(:,:,:,:,C0),sr)

      ! F(x)*S(x) + F(y)*S(y) + F(y)*S(y):
      !   F(g)*S(0) -SR-> F(g,x)*S(0,x) + F(g,y)*S(0,y) + F(g,z)*S(0,z),
      !   F(0)*S(g) -SR-> F(0,x)*S(g,x) + F(0,y)*S(g,y) + F(0,z)*S(g,z):
      call grFS(F4(:,:,:,:,CX),S5(:,:,:,:,CX,C0),sr, fact=WDB-WDA)
      call grFS(F4(:,:,:,:,CY),S5(:,:,:,:,CY,C0),sr, fact=WDB-WDA)
      call grFS(F4(:,:,:,:,CZ),S5(:,:,:,:,CZ,C0),sr, fact=WDB-WDA)

      ! Fp2:
      call doP2IL(  LA+LB+1,W2, IL)
      ! D3 Fp2:
      call doD4F(LA,LB,1,0,IL,YL(:,:,:,:,:,C0,:),F4)

      ! Fp2(g)*S(0), Fp2(0)*S(g):
      call grFS(F4(:,:,:,:,C0),S5(:,:,:,:,C0,C0),sr, fact=WDB*WDA)
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_rel_gr_nuc

  subroutine shgi_rel_sd_nuc(z,nu,sr,RAD)
    ! Returns W,D derivatives of dimension 6x6 !
    use shgi_common, only: W2, NORM, LAMBDA
    use shgi_rad, only: doIL, doP2IL
    implicit none
    real(RK)   , intent(in)    :: z
    real(RK)   , intent(inout) :: nu(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK)   , intent(inout) :: sr(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK)   , intent(in)    :: RAD
    optional :: nu
    optional :: sr
    optional :: RAD
    ! *** end of interface ***

    real(RK) :: R
    logical  :: nucl
    logical  :: srel

    if (Z == ZERO) RETURN

    nucl = present(nu)
    srel = present(sr)

    R = ZERO
    if(present(RAD)) R = max(RAD,ZERO)

    if(R == ZERO)then
       call point(LAMBDA, z * NORM(:,3))
    else
       call finit()
    endif

  contains

    subroutine point(LAM,NRM)
      use shgi_common, only: YL, S5, K4, WDA, WDB
      use shgi_dnf, only: doD5F
      use shgi_utils, only: sdFS, sdSYM
      implicit none
      real(RK), intent(in) :: LAM(:), NRM(:) ! (NAB)
      ! *** end of interface **

      integer(IK) :: LX=1, ND=2
      real(RK)    :: IL(NAB,1+LA+LB+4)
      real(RK)    :: F5(NAB,(LA+1)**2,(LB+1)**2,4,4,4) ! 4 = (P+1)**2

      LX = 1 ! one more nabla for SR
      ND = 2 ! two more orders for SR
      if( .not. srel )then
         LX = 0
         ND = 0
      endif

      call doIL(LA+LB+ND+2,W2,LAM, NRM, IL)
      ASSERT(LC>=1)
      ASSERT(LD>=1)
      ASSERT(LE>=LX)
      call doD5F(LA,LB,1,1,LX,IL,YL,F5)

      if( .not. nucl ) goto 111
      ! *** NUCLEAR ATTRACTION ***

      ! F(i,j)*S(0), F(i,0)*S(0,j), F(0)*S(i,j):
      call sdFS(F5(:,:,:,:,:,C0),S5(:,:,:,:,:,C0),nu)

      call sdSYM(nu)

!     ! transform and add GW and GD to GA, GB, and GC
!     call shgi_sd_wd_to_abc(sd,nu)

111   if( .not. srel ) RETURN
      ! *** SCALAR RELATIVISTIC MAT EL ***

      ! F*Sp2:
      call sdFS(F5(:,:,:,:,:,C0),K4(:,:,:,:,:),sr)

      ! F(x)*S(x) + F(y)*S(y) + F(y)*S(y):
      call sdFS(F5(:,:,:,:,:,CX),S5(:,:,:,:,:,CX),sr, fact=WDB-WDA)
      call sdFS(F5(:,:,:,:,:,CY),S5(:,:,:,:,:,CY),sr, fact=WDB-WDA)
      call sdFS(F5(:,:,:,:,:,CZ),S5(:,:,:,:,:,CZ),sr, fact=WDB-WDA)

      ! Fp2:
      call doP2IL(LA+LB+2,W2, IL)
      ! D4 Fp2:
      call doD5F(LA,LB,1,1,0,IL,YL,F5)

      ! Fp2*S:
      call sdFS(F5(:,:,:,:,:,C0),S5(:,:,:,:,:,C0),sr, fact=WDB*WDA)

      call sdSYM(sr)

!      ! transform and add GW and GD to GA, GB, and GC
!      call shgi_gr_wd_to_abc(gr,sr)
    end subroutine point

    subroutine finit()
      ! *** end of interface **

      real(RK) :: GAM, LAM(NAB), NRM(NAB)

      GAM = (THREE/TWO) / R**2
      LAM = ( LAMBDA * GAM ) / ( LAMBDA + GAM )
      NRM = z * NORM(:,3) * SQRT( GAM / ( LAMBDA + GAM ))

      call point(LAM,NRM)
    end subroutine finit

  end subroutine shgi_rel_sd_nuc
#endif

  !--------------- End of module ----------------------------------
end module shgi_relnuc

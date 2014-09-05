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
module shgi_ab
  !-------------------------------------------------------------------
  !
  ! Computes integral preliminaries for a pair of shells.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
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
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind, & ! type specification parameters
       I8K=>i8_kind
  use shgi_cntrl
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: shgi_set_ab
  public :: shgi_close_ab
  public :: shgi_set_ovrl
  public :: shgi_kin
  public :: shgi_gr_kin
  public :: shgi_sd_kin
  public :: shgi_atomic_ints
  public :: shgi_atomic_coul

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

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi_set_ab(iLA,iLB,iXA,iXB,ea,eb)
    use shgi_common
!   use shgi_pseudo, only: AEXP, BEXP, XA, XB
    use constants
    implicit none
    integer(IK), intent(in)  :: iLA,iLB
    real(RK),    intent(in)  :: iXA(3), iXB(3)
    real(RK),    intent(in)  :: ea(:), eb(:)
    ! *** end of interface ***

!   integer(IK) :: warn_count = 0
!   integer(IK) :: warn_count1 = 0

    integer(IK) :: memstat
    integer(IK) :: ia,ib,iab
    integer(IK) :: L
    real(RK)    :: a,b

    ! set globals:
    LA   = iLA
    LB   = iLB
    LAB  = max(LA,LB)
    NEA  = size(ea)
    NEB  = size(eb)

    XA   = iXA ! only needed for PP
    XB   = iXB ! only needed for PP
    XD   = XA - XB
    XD2  = XD(1)**2 + XD(2)**2 + XD(3)**2

    ! store copies of ALPHA and BETA exponents, only needed for PP
     allocate(AEXP(NEA),BEXP(NEB))
     AEXP = ea
     BEXP = eb

    ! integral screening (CUTOFF) by overlap > exp(-MAXEXP):
    ! MAXEXP in shgi_cntrl is set by shgi_set_maxexp()
    ! prepare mask for screening integrals:
    allocate(CUTOFF(NEA,NEB),stat=memstat)
    ASSERT(memstat==0)

    NAB = 0
    do ib=1,NEB
       b = eb(ib)
       do ia=1,NEA
          a = ea(ia)
          if( ( a*b/(a+b) ) * XD2 < MAXEXP ) then
             CUTOFF(ia,ib) =  .true.
             NAB = NAB + 1
          else
             CUTOFF(ia,ib) =  .false.
          endif
       enddo
    enddo
    DPRINT ' NAB/(NEA*NEB)%=',100*real(NAB)/(NEA*NEB)

    ! collecting statistics (vars declared in shgi_cntrl):
    shgi_stat_num_ints = shgi_stat_num_ints + NEA * NEB * (2*LA+1) * (2*LB+1)
    shgi_stat_screened = shgi_stat_screened +    NAB    * (2*LA+1) * (2*LB+1)

!   if( NAB /= NEA*NEB )then
!     if( warn_count1 < 10 )then
!       WARN('***************** NAB /= NA*NB ********************')
!     else if( warn_count1 == 10 )then
!       WARN('***************** NAB /= NA*NB, will ignore further warings ...')
!     endif
!     warn_count1 = warn_count1 + 1
!   endif
!   if( NAB == 0 )then
!     if( warn_count < 10 )then
!       WARN('***************** NAB == 0 ********************')
!     else if( warn_count == 10 )then
!       WARN('***************** NAB == 0, will ignore further warings ...')
!     endif
!     warn_count = warn_count + 1
!   endif
    ! NAB was set

    ! FIXME: some are needed for the third center only:
    allocate(ZETA(NAB),LAMBDA(NAB),NORM(NAB,4),stat=memstat)
    ASSERT(memstat==0)

    allocate(W(NAB,3),WC(NAB,3),W2(NAB),WDA(NAB),WDB(NAB),stat=memstat)
    ASSERT(memstat==0)

    iab = 0
    do ib=1,NEB
       do ia=1,NEA
          if(.not.CUTOFF(ia,ib)) cycle
          iab = iab + 1

          b = eb(ib)
          a = ea(ia)

          ZETA(iab)   = ( a * b ) / ( a + b )
          LAMBDA(iab) =   a + b

          W(iab,:)  = ( a * xa + b * xb ) / ( a + b )
          WDA(iab)  = a / ( a + b )
          WDB(iab)  = b / ( a + b )

          ! common pre-factor for norms:
          NORM(iab,1)  = ((TWO*a)/PI)**(THREE/FOUR) * ((TWO*b)/PI)**(THREE/FOUR) &
                       / sqrt( a**LA * DFAC(LA)     *      b**LB * DFAC(LB)    )
          ! norm for overlap:
          NORM(iab,2)  = ( PI / ( a + b ) )**(THREE/TWO)      * NORM(iab,1)
          ! norm for Nuclear attraction
          NORM(iab,3)  = ( ( a + b ) / PI )**(  ONE/TWO) *TWO * NORM(iab,2)

          ! norm for 2c Coulomb norm integrals
          NORM(iab,4)  = TWO * PI**(FIVE/TWO) / ( a*b * SQRT(a+b) )
       enddo
    enddo

    ![[=== Factors for conversion ==================
    allocate(WDAL(NAB,1+LA),WDBL(NAB,1+LB),stat=memstat)
    ASSERT(memstat==0)

    ! D(lma)D(lmb) D(lmc)I -> DA(lma)DB(lmb) D(lmc)I
    WDAL(:,1) = ONE
    do L=2,1+LA
       WDAL(:,L) = WDAL(:,L-1) * WDA
    enddo

    WDBL(:,1) = ONE
    do L=2,1+LB
       WDBL(:,L) = WDBL(:,L-1) * WDB
    enddo
    !]]=============================================

    !MOVE OUT:   call shgi_set_lcde(LC,LD,LE)
    !MOVE OUT:   call shgi_2c()

    ! MOVE OUT: ! will be possibly moved out:
    ! MOVE OUT: call shgi_glob_alloc(NAB,LA,LB)
  end subroutine shgi_set_ab

  subroutine shgi_close_ab()
    use shgi_common
!   use shgi_pseudo, only: AEXP, BEXP
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat

    LA   = -1
    LB   = -1
    LAB  = -1
    LMAX = -1
    NEA  = -1
    NEB  = -1
    NAB  = -1

    if( allocated(W) )then
       deallocate(W,WC,W2,WDA,WDB,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(WDAL) )then
       deallocate(WDAL,WDBL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(CUTOFF) )then
       deallocate(CUTOFF,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(ZETA) )then
       deallocate(ZETA,LAMBDA,NORM,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(X5) )then
       deallocate(X5,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(S5) )then
       deallocate(S5,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(K4) )then
       deallocate(K4,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(YL) )then
       deallocate(YL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(YS) )then
       deallocate(YS,stat=memstat)
       ASSERT(memstat==0)
    endif

    ! only for GRADS: ------------------>
    if( allocated(GS) )then
       deallocate(GS,stat=memstat)
       ASSERT(memstat==0)
    endif
    !<-----------------------------------

    if( allocated(AEXP) )then
       deallocate(AEXP,BEXP,stat=memstat)
       ASSERT(memstat==0)
    endif

    ! MOVED OUT: ! will be possibly moved out:
    ! MOVED OUT: call shgi_glob_free()
  end subroutine shgi_close_ab

  subroutine shgi_set_ovrl(LA,LB,LC,LD,LE)
    ! Precomputes (unnormalized) overlap ints  S5(lma,lmb,lmc,lmd,lme)
    ! and corresponding laplace (kinetic) ints K4(lma,lmb,lmc,lmd)
    ! both are use-associated from shgi_common:
    use shgi_common, only: XD, XD2, ZETA, X5, S5, K4
    use shgi_shr, only: SHR_D5Av
    use shgi_rad, only: doSL, doP2SL
    use shgi_dnf, only: doD4S, doD5S
    implicit none
    integer(IK), intent(in) :: LA,LB,LC,LD,LE ! shadow global vars!
    ! *** end of interface ***

    real(RK) :: SL(NAB,1+LA+LB+LC+LD+LE+2) ! ..+2 needed for laplacian nabla^2

    integer(IK) :: memstat

    allocate(X5((LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2,1+LA+LB+LC+LD+LE),&
         stat=memstat)
    ASSERT(memstat==0)

    ! compute angluar part of D = ( A - B )
    call SHR_D5Av(1,LA,LB,LC,LD,LE,XD,X5)

    ! derivatives of overlap integral:
    call doSL(LA+LB+LC+LD+LE+2,XD2,ZETA,SL)

    allocate(S5(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2,(LE+1)**2),&
         stat=memstat)
    ASSERT(memstat==0)

    ! nablas in PVSP reduce to laplacian:
    allocate(K4(NAB,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2),&
         stat=memstat)
    ASSERT(memstat==0)

    call doD5S(LA,LB,LC,LD,LE,SL,X5,S5)

    !
    ! From this point on the (unnormalized) OVERLAP is in S5(...)
    !

    !============================================
    ! FXIME: the second part may be made optional:
    !
    ! now replace by laplacian of S:
    call doP2SL(LA+LB+LC+LD,XD2,SL)
    SL = -SL ! -nabla^2
    call doD4S(LA,LB,LC,LD,SL,X5(:,:,:,:,C0,:),K4)

    !
    ! From this point on the (unnormalized) -LAPLACE is in K4(...)
    !
  end subroutine shgi_set_ovrl

  subroutine shgi_kin(PIOVRL,PIKNTC)
    ! scales overlap ints S5 and Laplace ints K4 by
    ! normalization factors and stores final OVRL/KNTC ints
    use shgi_common, only: S5, K4, NORM
    implicit none
    real(RK), intent(out) :: PIOVRL(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK), intent(out) :: PIKNTC(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ma,mb

    FPP_TIMER_START(t2c)

    if(is_on(IAEQB).and.LA==LB)then
      ! A==B, LA==LB, fallthrough to atomic version:
      call shgi_atomic_kin(PIOVRL,PIKNTC)
      RETURN
    endif

    ASSERT(allocated(S5))
    do mb=1,2*LB+1
    do ma=1,2*LA+1
      PIOVRL(:,ma,mb) = S5(:,LA**2+ma,LB**2+mb,C0,C0,C0) * NORM(:,2)
    enddo
    enddo

    ASSERT(allocated(K4))
    do mb=1,2*LB+1
    do ma=1,2*LA+1
      PIKNTC(:,ma,mb) = K4(:,LA**2+ma,LB**2+mb,C0,C0) * NORM(:,2) / 2
    enddo
    enddo

    FPP_TIMER_STOP(t2c)
  end subroutine shgi_kin

  subroutine shgi_atomic_kin(PIOVRL,PIKNTC)
    ! employs diagonal structure in magentic quantum numbers
    ! (mainly to test atomic version)
    use shgi_common, only: CUTOFF,AEXP ! for atomic version
    implicit none
    real(RK), intent(out) :: PIOVRL(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK), intent(out) :: PIKNTC(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: m
    real(RK)    :: S(NEA,NEB),T(NEA,NEB) ! for atomic version

    ! A==B, LA==LB, make sure atomic version works:
    ASSERT(is_on(IAEQB))
    ASSERT(LA==LB)
    ASSERT(NEA==NEB)

    call shgi_atomic_ints(LA,AEXP,S,T)

    PIOVRL = 0.0
    do m=1,2*LA+1
      PIOVRL(:,m,m) = PACK(S,CUTOFF)
    enddo

    PIKNTC = 0.0
    do m=1,2*LA+1
      PIKNTC(:,m,m) = PACK(T,CUTOFF)
    enddo
  end subroutine shgi_atomic_kin

  subroutine shgi_gr_kin(GROVRL,GRKNTC)
    ! scales overlap ints S5 and Laplace ints K4 by
    ! normalization factors and stores final OVRL/KNTC grads
    use shgi_common, only: S5, K4, NORM
    use shgi_cntrl , only: is_on,IAEQB
    implicit none
    real(RK), intent(out) :: GROVRL(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK), intent(out) :: GRKNTC(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    ! *** end of interface ***

    integer(IK) :: ma,mb

    FPP_TIMER_START(t2c)

    if(is_on(IAEQB))then
      ! single-center (A==B) gradients sum up to zero
      ! anyway GA+GB==0 -- dont even compute them:
      GROVRL = 0.0
      GRKNTC = 0.0
      RETURN
    endif

    ASSERT(allocated(S5))
    ASSERT(LC>=1)
    do mb=1,2*LB+1
    do ma=1,2*LA+1
    ! FIXME: sign?
    GROVRL(:,ma,mb,GAX)   = - S5(:,LA**2+ma,LB**2+mb,CX,C0,C0) * NORM(:,2)
    GROVRL(:,ma,mb,GAY)   = - S5(:,LA**2+ma,LB**2+mb,CY,C0,C0) * NORM(:,2)
    GROVRL(:,ma,mb,GAZ)   = - S5(:,LA**2+ma,LB**2+mb,CZ,C0,C0) * NORM(:,2)
    enddo
    enddo

    ASSERT(allocated(K4))
    ASSERT(size(K4,4)>=4)
    do mb=1,2*LB+1
    do ma=1,2*LA+1
    ! FIXME: sign?
    GRKNTC(:,ma,mb,GAX)   = - K4(:,LA**2+ma,LB**2+mb,CX,C0) * NORM(:,2) / 2
    GRKNTC(:,ma,mb,GAY)   = - K4(:,LA**2+ma,LB**2+mb,CY,C0) * NORM(:,2) / 2
    GRKNTC(:,ma,mb,GAZ)   = - K4(:,LA**2+ma,LB**2+mb,CZ,C0) * NORM(:,2) / 2
    enddo
    enddo

    FPP_TIMER_STOP(t2c)
  end subroutine shgi_gr_kin

  subroutine shgi_sd_kin(SDOVRL,SDKNTC)
    ! scales overlap ints S5 and Laplace ints K4 by
    ! normalization factors and stores final OVRL/KNTC grads
    use shgi_common, only: S5, K4, NORM
    implicit none
    real(RK), intent(out) :: SDOVRL(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    real(RK), intent(out) :: SDKNTC(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    ! *** end of interface ***

    integer(IK) :: ma,mb

    FPP_TIMER_START(t2c)

    ASSERT(allocated(S5))
    ASSERT(LC>=1)
    ASSERT(LD>=1)
    do mb=1,2*LB+1
    do ma=1,2*LA+1
    SDOVRL(:,ma,mb,GAX,GAX) = S5(:,LA**2+ma,LB**2+mb,CX,CX,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAY,GAX) = S5(:,LA**2+ma,LB**2+mb,CY,CX,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAZ,GAX) = S5(:,LA**2+ma,LB**2+mb,CZ,CX,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAX,GAY) = S5(:,LA**2+ma,LB**2+mb,CX,CY,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAY,GAY) = S5(:,LA**2+ma,LB**2+mb,CY,CY,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAZ,GAY) = S5(:,LA**2+ma,LB**2+mb,CZ,CY,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAX,GAZ) = S5(:,LA**2+ma,LB**2+mb,CX,CZ,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAY,GAZ) = S5(:,LA**2+ma,LB**2+mb,CY,CZ,C0) * NORM(:,2)
    SDOVRL(:,ma,mb,GAZ,GAZ) = S5(:,LA**2+ma,LB**2+mb,CZ,CZ,C0) * NORM(:,2)
    enddo
    enddo

    ASSERT(allocated(K4))
    ASSERT(size(K4,4)>=4)
    ASSERT(size(K4,5)>=4)
    do mb=1,2*LB+1
    do ma=1,2*LA+1
    SDKNTC(:,ma,mb,GAX,GAX) = K4(:,LA**2+ma,LB**2+mb,CX,CX) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAY,GAX) = K4(:,LA**2+ma,LB**2+mb,CY,CX) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAZ,GAX) = K4(:,LA**2+ma,LB**2+mb,CZ,CX) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAX,GAY) = K4(:,LA**2+ma,LB**2+mb,CX,CY) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAY,GAY) = K4(:,LA**2+ma,LB**2+mb,CY,CY) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAZ,GAY) = K4(:,LA**2+ma,LB**2+mb,CZ,CY) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAX,GAZ) = K4(:,LA**2+ma,LB**2+mb,CX,CZ) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAY,GAZ) = K4(:,LA**2+ma,LB**2+mb,CY,CZ) * NORM(:,2) / 2
    SDKNTC(:,ma,mb,GAZ,GAZ) = K4(:,LA**2+ma,LB**2+mb,CZ,CZ) * NORM(:,2) / 2
    enddo
    enddo

    FPP_TIMER_STOP(t2c)
  end subroutine shgi_sd_kin

  subroutine shgi_atomic_ints(L,e,S,T,V,O,rad)
    !
    ! Computes atomic integrals
    !
    use constants ! DFAC(L)=(2L-1)!!, FAC(L)=L!
    implicit none
    integer(IK), intent(in)  :: L
    real(RK),    intent(in)  :: e(:)   ! (n_exp)
    real(RK),    intent(out) :: S(:,:) ! (n_exp,n_exp)
    real(RK),    intent(out) :: T(:,:) ! (n_exp,n_exp)
    real(RK),    intent(out) :: V(:,:) ! (n_exp,n_exp)
    real(RK),    intent(out) :: O(:,:) ! (n_exp,n_exp)
    real(RK),    intent(in)  :: rad    ! finite nucleus size
    optional :: V, O, rad
    ! *** end of interface ***

    integer(IK) :: ne,i,j
    real(RK)    :: a,b
    real(RK), dimension(size(e),size(e)) :: zet,lam
    real(RK), dimension(size(e))         :: norm
    real(RK)                             :: R
    real(RK)                             :: gam(1),cof(1)

    ne = size(e)
    ASSERT(ne==size(S,1))
    ASSERT(ne==size(S,2))
    ASSERT(ne==size(T,1))
    ASSERT(ne==size(T,2))
    if(present(V))then
    ASSERT(present(O))
    ASSERT(present(rad))
    ASSERT(ne==size(V,1))
    ASSERT(ne==size(V,2))
    ASSERT(ne==size(O,1))
    ASSERT(ne==size(O,2))
    endif

    do j=1,ne
      b = e(j)
    do i=1,ne
      a = e(i)

      zet(i,j) = ( a * b ) / ( a + b )
      lam(i,j) =   a + b
      ! zeta/lambda == a*b/(a+b)^2
    enddo
    enddo

    ! (unnormalized) ll-overlap:
    S = 1 / lam**(L+THREE/TWO)

    ! ll-kinetic:
    T = S * zet * (2*L+3)

    if(present(V))then
      ! potential and sr matrix elements:
      R = max(rad,0.0_rk)

      if( R == 0 )then ! point nucleus

        ! nuclear attraction:
        V = S * sqrt(lam)       &
              * gamma2(L+1,1)   &
              / gamma2(2*L+3,2) ! Gamma(L+1)/Gamma(L+3/2)

        ! scalar relativity:
        O = V * zet * 4 * ( 1 + L )

        ! higher momenta get an additional term:
        if( L > 0 )then
          O = O + V * lam
        endif

      else ! finite nucleus
        print *,'shgi_atomic: using finite radius R=',R

        ! exponent of the gaussian distribution, that corresponds to radius:
        gam = 3 / ( 2 * R**2 )      ! array of size 1

        ! use an s-Gaussian of unit charge:
        cof = (gam/PI)**(THREE/TWO) ! array of size 1

        ! the next call will INCREMENT them:
        V = 0
        O = 0
        call coul(L,lam,zet,gam,cof,S,V,O) ! computes Coulomb screening
      endif
    endif

    ! GTOs are _normalized_ (by square norm), record the norms:
    do i=1,ne
      norm(i) = 1/sqrt(S(i,i))
    enddo

    ! renorm the ovl/kin ints:
    do j=1,ne
    do i=1,ne
      S(i,j) = S(i,j) * (norm(i) * norm(j))
      T(i,j) = T(i,j) * (norm(i) * norm(j))
    enddo
    enddo
    if(.not.present(V)) RETURN

    ! renorm the nuc/sr-nuc ints:
    do j=1,ne
    do i=1,ne
      V(i,j) = V(i,j) * (norm(i) * norm(j))
      O(i,j) = O(i,j) * (norm(i) * norm(j))
    enddo
    enddo
  end subroutine shgi_atomic_ints

  subroutine shgi_atomic_coul(L,e,sexp,scof,V,O)
    ! Computes atomic integrals
    use constants ! DFAC(L)=(2L-1)!!, FAC(L)=L!
    implicit none
    integer(IK), intent(in)  :: L
    real(RK),    intent(in)  :: e(:)   ! (n_exp)
    real(RK),    intent(in)  :: sexp(:) ! (n_gam) fit exponents
    real(RK),    intent(in)  :: scof(:) ! (n_gam) fit coefficients
    real(RK),    intent(out) :: V(:,:) ! (n_exp,n_exp)
    real(RK),    intent(out) :: O(:,:) ! (n_exp,n_exp)
    ! *** end of interface ***

    integer(IK) :: ne,i,j
    real(RK)    :: a,b
    real(RK), dimension(size(e),size(e)) :: zet,lam,S
    real(RK), dimension(size(e))         :: norm

    ne = size(e)
    ASSERT(ne==size(V,1))
    ASSERT(ne==size(V,2))
    ASSERT(ne==size(O,1))
    ASSERT(ne==size(O,2))
    ASSERT(size(sexp)==size(scof))

    do j=1,ne
      b = E(j)
    do i=1,ne
      a = E(i)

      zet(i,j) = ( a * b ) / ( a + b )
      lam(i,j) =   a + b
      ! zeta/lambda == a*b/(a+b)^2
    enddo
    enddo

    ! (unnormalized) ll-overlap:
    S = 1 / lam**(L+THREE/TWO)

    ! the screening density parametrized by s-fit was provided --
    ! compute its Coulomb potential:
    call coul(L,lam,zet,sexp,scof,S,V,O)

    ! GTOs are _normalized_ (by square norm), record the norms:
    do i=1,ne
      norm(i) = 1/sqrt(S(i,i))
    enddo

    ! scale the ints:
    do j=1,ne
    do i=1,ne
      V(i,j) = V(i,j) * (norm(i) * norm(j))
      O(i,j) = O(i,j) * (norm(i) * norm(j))
    enddo
    enddo
  end subroutine shgi_atomic_coul

  subroutine coul(L,lam,zet,gam,cof,S,V,O)
    !
    ! provided with exponents and coefficients of the fit
    ! computes the (sum of the) Coulomb and SR-Coulomb ints
    !
    use constants
    implicit none
    integer(IK), intent(in)  :: L
    real(RK),    intent(in)  :: lam(:,:) ! (n_exp,n_exp): (a+b)
    real(RK),    intent(in)  :: zet(:,:) ! (n_exp,n_exp): a*b/(a+b)
    real(RK),    intent(in)  :: gam(:)   ! (n_fit)      : exponents of the fit
    real(RK),    intent(in)  :: cof(:)   ! (n_fit)      : coefficients of the fit
    real(RK),    intent(in)  :: S(:,:)   ! (n_exp,n_exp): overlap, a common factor INPUT!
    real(RK),    intent(out) :: V(:,:)   ! (n_exp,n_exp): V = S * ... (Coulomb)
    real(RK),    intent(out) :: O(:,:)   ! (n_exp,n_exp): O = S * ... (SR of Coulomb)
    ! *** end of interface ***

    real(RK), dimension(size(S,1),size(S,2)) :: U
    real(RK)    :: gmm,cff
    integer(IK) :: i

    ASSERT(size(gam)==size(cof))
    ASSERT(all(shape(lam)==shape(zet)))
    ASSERT(all(shape(lam)==shape(U)))
    ASSERT(all(shape(lam)==shape(V)))
    ASSERT(all(shape(lam)==shape(O)))

    ! first potential of the *point* nucleus:
    U = S * sqrt(lam)       &
          * gamma2(L+1,1)   &
          / gamma2(2*L+3,2) ! Gamma(L+1)/Gamma(L+3/2)

    V = 0
    O = 0
    ! add contributions of each fit-function:
    do i=1,size(gam) ! over fit-exponents
      gmm = gam(i)
      cff = cof(i)

      ! The expressions below are for NORMALIZED gaussians.
      ! The fit-functions are NOT normalized, therefore:
      cff = cff * (PI/gmm)**(THREE/TWO)

      ! scalar relativity:
      O = O + U * sqrt(gmm/(lam+gmm))    &
                * zet * 4 * ( 1 + L )    &
                * pol(L+1,lam/(lam+gmm)) & ! a polynomial
                * cff

      ! higher momenta get an additional term:
      if( L > 0 )then
        O = O + U * lam                          &
                  * sqrt(gmm/(lam+gmm))          &
                  * (                            &
                      pol(L-1, lam/(lam+gmm))    &
                    - 2 * L * (lam/(lam+gmm))**L &
                            * gamma2(2*L+1,2)    &
                            / gamma2(L+1,1)      &
                            / gamma2(1,2)        &
                    )                            &
                  * cff
      endif

      ! finally, the potential itself:
      V = V + U * sqrt(gmm/(lam+gmm))  &
                * pol(L,lam/(lam+gmm)) & ! a polynomial
                * cff
    enddo
  contains
    elemental function pol(l,z) result(f)
      !
      ! Computes the polynomial:
      !
      !            l                               l
      !           ___                             ___
      !           \       Gamma(k+1/2)        k   \   (2k-1)!!     k
      ! pol(l,z) = >  -------------------- * z  =  >  --------- * z
      !           /__ Gamma(k+1)Gamma(1/2)        /__  (2k)!!
      !
      !           k=0                             k=0
      !
      implicit none
      integer(IK), intent(in) :: l
      real(RK),    intent(in) :: z
      real(RK)                :: f ! result
      ! *** end of interface **

      integer(IK) :: k
      real(RK)    :: n,d,zk

      f = 1
      if( l == 0 ) RETURN

      n  = 1 ! numerator
      d  = 1 ! denominator
      zk = 1 ! z^k power
      do k=1,l
        n  = n * (2*k-1)
        d  = d * (2*k)
        zk = zk * z
        f  = f + n * zk / d
      enddo
    end function pol
  end subroutine coul

  function gamma2(a,b) result(gam)
    !
    ! gamma2(a,b) = Gamma(a/b)
    !
    ! currently only for integer and half-integer a/b
    !
    use constants, only: PI
    use debug,     only: NaN
    implicit none
    integer(IK), intent(in) :: a, b
    real(RK)                :: gam  ! result
    ! *** end of interface **

    integer(IK) :: x

    x   = a
    gam = 1

    ! recurrence untli (x/b) <= 1:
    do while ( x > b )
      ! Gamma(z) = (z-1) * Gamma(z-1):
      x = x - b
      gam = x * gam / b
    enddo
    ! multiply by Gamma(x/b) with argument (x/b) <= 1:
    if( x == b ) RETURN           ! Gamma(1) = 1

    if( x == 1 .and. b == 2 )then ! Gamma(1/2) = sqrt(pi)
      gam = gam * sqrt(PI)
      RETURN
    endif

    ! other cases not implemented, return error:
    gam = NaN()
  end function gamma2


  !--------------- End of module -------------------------------------
end module shgi_ab

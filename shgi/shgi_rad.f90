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
module shgi_rad
  !-------------------------------------------------------------------
  !
  ! Radial part for various integrals.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
  ! Copyright (c) 2007 Alexey Shor
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
  use shgi_cntrl ! only: whatis, IXEQC,IAEQC,IBEQC,MAXEXP and TIMERS
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: doSL
  public :: doEL
  public :: doEIL
  public :: doIL!(LL,w2,lambda,norm,IL)
  public :: doFL!(LL,w2,lambda,norm,IL,rad) -- wrapper for the prev for finite nuc size
  public :: doRL
#ifndef NO_PSEUDO /* uses unique_atom_module, needs much of PG sources */
  public :: doPQ1L
  public :: doPQ2L
#endif
  public :: doQSSL
  public :: doQSRL
  public :: doQRRL
  public :: doP2SL ! scalar w2
  public :: doP2IL ! vector w2(:)
  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine doSL(LL,d2,zeta,SL)
    ! computes the derivatives of the
    ! "overlap" integral
    !   S(0) := exp(-z*d2)=exp(-2z*x)
    !   S(L) := S^L(x) where x=d2/2
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: d2, zeta(:)
    real(RK)   , intent(out) :: SL(:,:)
    ! *** end of interface ***

    integer(IK) :: L

    ASSERT(1+LL<=size(SL,2))

    SL(:,1) = exp( - zeta * d2 )
    do L=2,1+LL ! 1+LA+LB
       SL(:,L) = ( -TWO * zeta ) * SL(:,L-1)
    enddo
  end subroutine doSL

  subroutine doEL(LL,w2,lambda,norm,EL)
    ! computes the derivatives of the
    ! "Exponent" integral
    !   E(0) := exp(-z*d2)=exp(-2z*x)
    !   E(L) := E^L(x) where x=d2/2
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2(:),lambda(:),norm(:)
    real(RK)   , intent(out) :: EL(:,:)
    ! *** end of interface ***

    real(RK)    :: a(size(w2))
    integer(IK) :: L

    ASSERT(1+LL<=size(EL,2))

    a=1.0_RK
    EL(:,1) = norm*exp( - lambda * w2 )
    do L=2,1+LL
       a=( -TWO * lambda )*a
       EL(:,L) = a*EL(:,1)
    enddo

  end subroutine doEL

  subroutine doEIL(LL,w2,lambda1,lambda2,norm,EIL)
    ! computes the derivatives of the integral
    ! F(x)=norm*exp(-lam1*d2)*I0(lam2*d2)=
    !     =norm*exp(-2*lam1*x)*I0(2*lam2*x)=
    ! x=d2/2
    use gamma_module, only: gamma
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2(:),lambda1(:),lambda2(:),norm(:)
    real(RK)   , intent(out) :: EIL(:,:)
    ! *** end of interface ***
    integer(IK) :: L,i,k
    real(RK) :: IL(size(EIL,1),size(EIL,2))
    real(RK) :: EL(size(EIL,1),size(EIL,2))
    real(RK) :: a1(size(w2)),a2(size(w2))

    FPP_TIMER_START(tgam)
    IL(:,:1+LL) = gamma(1+LL, lambda2 * w2 )
    FPP_TIMER_STOP(tgam)

    a1=1.0_RK
    EL(:,1) = exp( -lambda1 * w2 )
    a2=1.0_RK
    IL(:,1) = a2 * IL(:,1)
    do L=2,1+LL
       a1=( -TWO * lambda1 )*a1
       EL(:,L) = a1*EL(:,1)

       a2 = ( -TWO * lambda2 )*a2
       IL(:,L) = a2*IL(:,L)
    enddo

    EIL(:,1)=norm*EL(:,1)*IL(:,1)
    do L=2,1+LL
       EIL(:,L)=zero
       do i=1,L
          k=min(L-i,i-1)
          EIL(:,L)=EIL(:,L)+binom(L-1,k)*EL(:,L-i+1)*IL(:,i)*norm
       end do
    end do

  end subroutine doEIL

  !subroutine doKL(LL,d2,zeta,KL,norm)
    !
    ! FIXME: convert to doP2IL(LL<d2,inout(SL))
    !
    ! computes the derivatives of the
    ! "kinetic" integral:
    !
    ! call doSL  (LL+2,XD2,ZETA,SL)
    ! call doP2SL(LL  ,XD2,     SL)
    !
    ! SL = -SL for -nabla^2
  !end subroutine doKL

  subroutine doIL(LL,w2,lambda,norm,IL)
    !
    ! Compute derivatives
    !
    !         L
    !       d
    !       --   I  ( 2 * lambda * x )
    !       dx    0
    !
    ! with
    !            2
    !       x = w  /  2
    !
    use gamma_module, only: gamma
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2(:), lambda(:), norm(:)
    real(RK)   , intent(out) :: IL(:,:)
    ! *** end of interface ***

    real(RK)    :: a(size(w2))
    integer(IK) :: L

    ASSERT(size(w2)==size(lambda))
    ASSERT(size(w2)==size(norm))
    ASSERT(size(w2)==size(IL,1))
    ASSERT(1+LL<=size(IL,2))

    FPP_TIMER_START(tgam)
    IL(:,:1+LL) = gamma(1+LL, lambda * w2 )
    FPP_TIMER_STOP(tgam)

    a = norm
    IL(:,1) = a * IL(:,1)
    do L=2,1+LL
       a = ( - TWO * lambda ) * a
       IL(:,L) = a * IL(:,L)
    enddo
  end subroutine doIL

  subroutine doFL(LL,w2,lambda,norm,IL,rad)
    !
    ! wrapper for the previous, lambda and norm for
    ! finite nucleus size
    !
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2(:), lambda(:), norm(:)
    real(RK)   , intent(out) :: IL(:,:)
    real(RK)   , intent(in)  :: rad
    ! *** end of interface ***

    real(RK) :: gam,lam(size(norm)),nrm(size(norm))

    if( rad == 0.0 )then
      call doIL(LL,w2,lambda,norm,IL)
    else
      ! adjust lambda/norm and call the same sub:
      gam = (THREE/TWO) / rad**2
      lam = ( lambda * gam ) / ( lambda + gam )
      nrm = norm * sqrt( gam / ( lambda + gam ))
      call doIL(LL,w2,lam,nrm,IL)
    endif
  end subroutine doFL

  subroutine doRL(LL,w2,lambda,f132,norm,IL)
    use gamma_module, only: gamma
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2(:), lambda(:), f132(:), norm(:)
    real(RK)   , intent(out) :: IL(:,:)
    ! *** end of interface ***

    real(RK)    :: a(size(w2)), b(size(w2))
    real(RK)    :: gam(size(w2),2+LL) ! one more for r2
    integer(IK) :: L

    ASSERT(size(w2)==size(lambda))
    ASSERT(size(w2)==size(norm))
    ASSERT(size(w2)==size(f132))
    ASSERT(size(w2)==size(IL,1))
    ASSERT(1+LL==size(IL,2))

    FPP_TIMER_START(tgam)
    gam = gamma(2+LL, lambda * w2 )
    FPP_TIMER_STOP(tgam)

    a = lambda * norm
    b = f132   * norm
    IL(:,1) = a * w2 * gam(:,2) + b * gam(:,1)
    do L=2,1+LL
       b = ( - TWO * lambda ) * b + TWO * a
       a = ( - TWO * lambda ) * a
       IL(:,L) = a * w2 * gam(:,L+1) + b * gam(:,L)
    enddo
  end subroutine doRL

#ifndef NO_PSEUDO /* uses unique_atom_module, needs much of PG sources */
  !
  ! For pseudopotential integrals
  !  local part (bessel_1):
  !
  subroutine doPQ1L(LL,w2,lambda,pp,NRM,IL)
    !
    ! Computes the derivatives (and value) of
    !
    !    F(x) = 4 * PI * exp(-lambda*W2) * Q_0(k)
    !      x  = W2/2
    !      k  = 2*lambda*|W|
    !    Q_L  = radial integrals with modified bessel function
    ! (all derivatives w.r.t. x)
    !
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use bessel_module, only: bessel_1
    implicit none
    integer(IK), intent(in)    :: LL
    real(RK)   , intent(in)    :: w2(:), lambda(:)
    type(ppt)  , intent(in)    :: pp
    real(RK)   , intent(in)    :: NRM(:)
    real(RK)   , intent(inout) :: IL(:,:) ! INOUT!: adds PP contributions
    ! *** end of interface ***

    real(RK)    :: k(size(w2)), e(size(w2)), q(size(w2),1+LL)
    real(RK)    :: bess(size(w2),1+LL)
    real(RK)    :: gamma,coeff
    integer(IK) :: i,L,pow

    FPP_TIMER_START(tpq1)

    DPRINT 'shgi_pseu: LL=',LL

    ASSERT(size(w2)==size(lambda))
    ASSERT(size(w2)==size(IL,1))
    ASSERT(size(w2)==size(NRM))
    ASSERT(1+LL==size(IL,2))

    k = TWO * lambda * sqrt(w2)
    e = - lambda * w2

    q = zero
    do i=1,pp%n_exponents
       coeff = pp%coefficients(i)
       gamma = pp%exponents(i)
       pow   = pp%powers(i)
       DPRINT 'shgi::doPQ1L: gam=',gamma,' pow=',pow,' cf=',coeff

       FPP_TIMER_START(tbs1)
       ! now computing the series of functions:
       !  exp(-lambda*w2) * Q_L, L=0,LL
       !
       ! bessel_1() multiplies result by exp(e), which is otherwise
       ! not used:
       bess = bessel_1(pow,1+LL,lambda+gamma,k,e,1)
       FPP_TIMER_STOP(tbs1)

       ! FIXME: nuclear attraction is historically positive:
       ! therefore to sum up we change the sign.
       ! FIXME: common 4PI factor, maybe somwhere else?
       coeff = - FOUR * PI * coeff
       q = q + coeff * bess
    enddo

    ! we obtained the series   : exp(-lambda*w^2) Q_L
    ! now do the derivatives of: exp(-lambda*w^2) Q_0
    ! w.r.t. x=w2/2.
    ! We invoke binomial rule and note that:
    !  a) derivative of exp() gives  -2*lambda    * exp(...)
    !  b) derivative of Q_L   gives (-2*lambda)^2 * Q_(L+1)

    ! Polynomials of X=(-2*lambda) with coeffs given by A(n)=binom(L,n)*q(n):
    !
    ! I(L) = SUM(n=0:L) binom(L,n) * Q(n) * (-2*lambda)^(L+n)
    !
    do L=0,LL
       E = q(:,1+L) ! * binom(L,L) == 1
       do i=L-1,0,-1
          E = E * (-TWO*lambda) + q(:,1+i) * binom(L,i)
       enddo
       E = E * (-TWO*lambda)**L
       ! FIXME: normalize and add to what we get (inout!)
       IL(:,1+L) = IL(:,1+L) + E * NRM
    enddo
    FPP_TIMER_STOP(tpq1)
  end subroutine doPQ1L

  !
  !  non-local part (bessel_2):
  !

  subroutine doPQ2L(LP,LLA,LLB,ac2,bc2,ab2,a,b,pp,NRM,IL,WDA,WDB,LAMBDA,CUTOFF)
    !
    ! Computes the derivatives (and value) of
    !
    !  F(x,y) = (4 * PI)^2 (2ab)^L * exp(-a*a^2-b*b^2) * Q(ka,kb,eta)
    !      x  = a^2/2
    !      y  = b^2/2
    !      ka = 2*a*|a|
    !      kb = 2*b*|b|
    !       Q = radial integrals with two modified bessel functions
    ! (all derivatives w.r.t. x and y)
    !
    use type_module, only: i8_kind
    use unique_atom_module, only: ppt=>unique_atom_pseudopot_type
    use bessel_module, only: bessel2, bessel_1
    use constants, only: DFAC ! double factorial for LP<=16!
    implicit none
    integer(IK), intent(in)  :: LP,LLA,LLB
    real(RK)   , intent(in)  :: ac2, bc2, ab2 ! distances squared
    real(RK)   , intent(in)  :: a(:), b(:) ! exponents
    type(ppt)  , intent(in)  :: pp
    real(RK)   , intent(in)  :: NRM(:)
    real(RK)   , intent(out) :: IL(:,:,:) ! (NAB,1+LLA,1+LLB)
    real(RK)   , intent(in)  :: WDA(:),WDB(:),LAMBDA(:) ! (NAB)
    logical    , intent(in)  :: CUTOFF(:,:)
    ! *** end of interface ***

    ! FIXME: bessel2() assumes this funny b-a-a-b order:
    real(RK)    :: bess2(size(b),size(a),1+LLA,1+LLB)
    real(RK)    :: Q(size(IL,1),1+LLA,1+LLB)   ! packed version of the above
    real(RK)    :: ka(size(IL,1)), kb(size(IL,1)) ! bessel2 contains the factor ka^(lp+la) * kb^(lp+la)
    real(RK)    :: bess1(size(IL,1),1+LLA+LLB) ! 1+LLA+LP, or 1+LLB+LP formally
    ! BESS1 used in the limit cases when either:
    ! * B=C =>
    !    need only ints with LB=LP =>
    !      need only ders up to Q2[LP](LA+LG,0) =>
    !        that are same as ders Q1(LP:LP+LA+LG) { remember LP=LB }
    ! * A=C =>
    !    need only ints with LA=LP =>
    !      need only ders up to Q2[LP](0,LB+LG) =>
    !        that are same as ders Q1(LP:LP+LB+LG) { remember LP=LA }
    !
    ! Therefore, in both cases one needs 1+L[AB]+LG ders from the range
    !   Q1(LP:LA+LB+LG) { remember LP=L[BA] }
    ! {LG=0,1 for Energy and Gradients, respectively}
    !
    ! It seems to be safe to use BESS1 even if
    ! * A=B=C =>
    !    need only ints with LA=LB=LP =>
    !      need only value Q2[LP](0,0) =>
    !        that is the same as der Q1(LP:LP) { remember LP=LA=LB }
    !
    ! All forces are zero in the one-center case {note no LG}.
    real(RK)    :: K(size(IL,1)), e(size(IL,1))
    real(RK)    :: m2a(size(IL,1)), m2b(size(IL,1))
    real(RK)    :: gamma,coeff
    real(RK)    :: ac,bc,ab
    integer(IK) :: i,pow
    integer(IK) :: ia,ib,iab
    integer(IK) :: LA,LB ! shadow global vars
    integer (i8_kind) :: same
    real(RK)    :: RE2LP1 ! real  2*LP+1
    real(RK)    :: DFACLP ! real (2*LP+1)!!

    FPP_TIMER_START(tpq2)

    DPRINT 'shgi::doPQ2L: LP=',LP,' LLA=',LLA,' LLB=',LLB

    ASSERT(size(NRM)==size(IL,1))
    ASSERT(1+LLA==size(IL,2))
    ASSERT(1+LLB==size(IL,3))

    ac = sqrt(ac2)
    bc = sqrt(bc2)
    ab = sqrt(ab2) ! A-B distance

    same = whatis (IXEQC)
    DPRINT 'shgi::doPQ2L: same=',same,' ac=',ac,' bc=',bc

    select case(same)
    case (0) ! all centers different
       ka = TWO * WDA * LAMBDA * ac
       kb = TWO * WDB * LAMBDA * bc
       bess2 = zero
    case (IAEQC) ! AC is zero
       ! bessel_1 works with packed storage:
       K  = TWO * WDB * LAMBDA * bc
       e  =     - WDB * LAMBDA * bc**2
       bess1 = zero
    case (IBEQC) ! BC is zero
       ! bessel_1 works with packed storage:
       K  = TWO * WDA * LAMBDA * ac
       e  =     - WDA * LAMBDA * ac**2
       bess1 = zero
    case (IXEQC) ! AC and BC are zeros
       ! FIXME: bessel_1 seems to manage K=0:
       K = zero
       e = zero
       bess1 = zero
    end select

    do i=1,pp%n_exponents
       coeff = pp%coefficients(i)
       gamma = pp%exponents(i)
       pow   = pp%powers(i)
       DPRINT 'shgi::doPQ2L: gam=',gamma,' pow=',pow,' cf=',coeff

       select case(same)
       case (0) ! all centers different
          FPP_TIMER_START(tbs2)
          ! bessel2() adds a "factor=coeff" of current contribution to
          ! its first argument:
          call bessel2(bess2,coeff,pow,LP,LLA,LLB,a,b,gamma,ac,bc,ab,MAXEXP)
          ! FIXME: screening in bessel2() requires distance AB and MAXEXP
          FPP_TIMER_STOP(tbs2)
       case default
          FPP_TIMER_START(tbs1)
          bess1 = bess1 + bessel_1(pow,1+LLA+LLB,LAMBDA+gamma,K,e,1) * coeff
          FPP_TIMER_STOP(tbs1)
       end select
    enddo

    RE2LP1 = REAL(2*LP+1,RK) !  2*LP+1
    DFACLP = DFAC(1+LP)      ! (2*LP+1)!!

    select case(same)
    case (0)
       ! go to packed storage again:
       do LB=0,LLB
          do LA=0,LLA
             iab = 0
             do ib=1,size(b)
                do ia=1,size(a)
                   ! FIXME: cutoff in bessel2 and here the same?:
                   if(.not.CUTOFF(ia,ib)) cycle
                   iab = iab + 1
                   ! FIXME: b-a-a-b order:
                   Q(iab,1+LA,1+LB) = bess2(ib,ia,1+LA,1+LB)! &
!!$                        / ( ka(iab)**(LP+LA) * kb(iab)**(LP+LB) )
                   ! FIXME: bessel2 contains
                   ! the power ka^(lp+la) * kb^(lp+lb)
                enddo
             enddo
          enddo
       enddo
    case (IAEQC)
       ! AC == ZERO, LA==LP
       Q = zero
       Q(:,1,:1+LLB) = bess1(:,1+LP:1+LP+LLB) / DFACLP
       ! FIXME: lower derivatives discarded! Dont compute!
    case (IBEQC)
       ! BC == ZERO, LB==LP
       Q = zero
       Q(:,:1+LLA,1) = bess1(:,1+LP:1+LP+LLA) / DFACLP
       ! FIXME: lower derivatives discarded! Dont compute!
    case (IXEQC)
       ! AC == BC == ZERO, LA==LB==LP
       Q = zero
       Q(:,1,1)      = bess1(:,1+LP)          / DFACLP
       ! FIXME: lower and higher derivatives discarded! Dont compute!
       ! FIXME: would be (2*LP+1)!!^2 if integral_radial were used
    end select

    ! compute -2a and -2b from global variables:
    m2a = -TWO * WDA * LAMBDA ! - 2 alpha
    m2b = -TWO * WDB * LAMBDA ! - 2 beta

    ! FIXME: ranges, some are zeros
    do LB=1,1+LLB
       do LA=1,1+LLA
          ! FIXME: normalize, maybe somewhere else
          Q(:,LA,LB) = Q(:,LA,LB) * NRM
          ! FIXME: common (4PI)^2 factor, maybe somwhere else?
          Q(:,LA,LB) = Q(:,LA,LB) * (FOUR * PI) !**2
          ! FIXME: common (2*LP+1) factor, maybe somwhere else?
          Q(:,LA,LB) = Q(:,LA,LB) * RE2LP1
          ! FIXME: (4ab)^LP factor:
!!$          Q(:,LA,LB) = Q(:,LA,LB) * (FOUR * ZETA * LAMBDA)**LP
          Q(:,LA,LB) = Q(:,LA,LB) * ( m2a * m2b )**LP
          ! FIXME: nuclear attraction is historically positive:
          ! therefore to sum up we change the sign.
          Q(:,LA,LB) = - Q(:,LA,LB)
       enddo
    enddo

    ! we obtained the series   : exp(-a*a^2-b*b^2) Q_LaLb
    ! now do the derivatives of:  ...............  Q_00
    ! w.r.t. x=a^2/2 and y=b^2/2
    ! We invoke binomial rule and note that:
    !  a) derivative of exp() gives  -2*a    * exp(...)
    !  b) derivative of Q_L   gives (-2*a)^2 * Q_(L+1)

    ! FIXME: ranges, some are zeros:
    IL = zero
    do LB=0,LLB
       do LA=0,LLA
          do ib=0,LB
             do ia=0,LA
                ! FIXME: performance here:
!!$                DPRINT 'IL(',LA,LB,')<-',&
!!$                     ' (-2a)^',LA+ia,binom(LA,ia),&
!!$                     ' (-2b)^',LB+ib,binom(LB,ib),&
!!$                     ' Q(',ia,ib,')'
                IL(:,1+LA,1+LB) = IL(:,1+LA,1+LB)  &
                     + m2b**(LB+ib) * binom(LB,ib) &
                     * m2a**(LA+ia) * binom(LA,ia) &
                     * Q(:,1+ia,1+ib)
             enddo
          enddo
       enddo
    enddo
    FPP_TIMER_STOP(tpq2)
  end subroutine doPQ2L
#endif

  recursive function binom(n,r) result(res)
    implicit none
    integer(IK), intent(in) :: n,r
    real(RK)                :: res
    ! *** end of interface ***

    if(n==r .or. r==0) then
       res = ONE
    else if (r==1) then
       res = real(n,RK)
    else
       res = real(n,RK)/real(n-r,RK)*binom(n-1,r)
    end if
  end function binom

  !
  ! These are for the 2c coulomb norm integrals.
  ! Difer from the rest in that they want just one
  ! distance w2 = d^2 = ( a - b )^2.
  !
  ! doSSL  < s| g(r1,r2) | s> and dervs
  ! doSRL  < s| g(r1,r2) |r2> and  <r2| g(r1,r2) | s> and dervs
  ! doRRL  <r2| g(r1,r2) |r2> and no dervs
  !

  subroutine doQSSL(LL,w2,lambda,norm,IL)
    use gamma_module, only: gamma
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2, lambda(:), norm(:)
    real(RK)   , intent(out) :: IL(:,:)
    ! *** end of interface ***

    real(RK)    :: a(size(lambda))
    integer(IK) :: L

    ASSERT(size(norm)==size(lambda))
    ASSERT(size(norm)==size(IL,1))
    ASSERT(1+LL==size(IL,2))

    FPP_TIMER_START(tgam)
    IL(:,:1+LL) = gamma(1+LL, lambda * w2 )
    FPP_TIMER_STOP(tgam)

    a = norm
    IL(:,1) = a * IL(:,1)
    do L=2,1+LL
       a = ( - TWO * lambda ) * a
       IL(:,L) = a * IL(:,L)
    enddo
  end subroutine doQSSL

  subroutine doQSRL(LL,w2,lambda,norm,IL,WDA,WDB,typ)
    use gamma_module, only: gamma
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2, lambda(:), norm(:)
    real(RK)   , intent(out) :: IL(:,:)
    real(RK)   , intent(in)  :: WDA(:),WDB(:) ! (NAB)
    integer(IK), intent(in)  :: typ
    ! *** end of interface ***

    real(RK)    :: a(size(lambda)), b(size(lambda))
    real(RK)    :: gam(size(lambda),2+LL) ! one more for sr2
    integer(IK) :: L

    ASSERT(size(norm)==size(lambda))
    ASSERT(size(norm)==size(IL,1))
    ASSERT(1+LL==size(IL,2))

    FPP_TIMER_START(tgam)
    gam = gamma(2+LL, lambda * w2 )
    FPP_TIMER_STOP(tgam)

    select case(typ)
    case( 1)
       a = (WDA/WDB * lambda    ) * norm
       b = (WDA/WDB + THREE/TWO ) * norm
    case(-1)
       a = (WDB/WDA * lambda    ) * norm
       b = (WDB/WDA + THREE/TWO ) * norm
    case default
       ABORT('no such typ')
    end select

    IL(:,1) = a * w2 * gam(:,2) + b * gam(:,1)
    do L=2,1+LL
       b = ( - TWO * lambda ) * b + TWO * a
       a = ( - TWO * lambda ) * a
       IL(:,L) = a * w2 * gam(:,L+1) + b * gam(:,L)
    enddo
  end subroutine doQSRL

  subroutine doQRRL(LL,w2,lambda,norm,IL,WDA,WDB)
    use gamma_module, only: gamma
    implicit none
    integer(IK), intent(in)  :: LL
    real(RK)   , intent(in)  :: w2, lambda(:), norm(:)
    real(RK)   , intent(out) :: IL(:,:)
    real(RK)   , intent(in)  :: WDA(:),WDB(:) ! (NAB)
    ! *** end of interface ***

    real(RK)    :: gam(size(lambda),3+LL) ! two more for r2r2

    ASSERT(LL==0)

    ASSERT(size(norm)==size(lambda))
    ASSERT(size(norm)==size(IL,1))
    ASSERT(1+LL==size(IL,2))

    FPP_TIMER_START(tgam)
    gam = gamma(3+LL, lambda * w2 )
    FPP_TIMER_STOP(tgam)

    IL(:,1) = norm * &
         ( &
           (THREE/TWO)*( (WDA/WDB) + (WDB/WDA) + (FIVE/TWO) ) &
           * gam(:,1) &
         + (THREE/TWO)*( (WDA/WDB) + (WDB/WDA) ) &
           * gam(:,2) *  lambda * w2 &
         +   gam(:,3) * (lambda * w2)**2 &
         )
  end subroutine doQRRL

  subroutine doP2IL(LL,w2,IL)
    implicit none
    integer(IK), intent(in)    :: LL
    real(RK)   , intent(in)    :: w2(:)
    real(RK)   , intent(inout) :: IL(:,:)
    ! *** end of interface ***

    real(RK)    :: a, b
    integer(IK) :: L

    ASSERT(size(w2)==size(IL,1))
    ASSERT(1+LL+2==size(IL,2))

    a =   ONE
    b = THREE
    IL(:,1) = a * w2 * IL(:,3) + b * IL(:,2)
    do L=2,1+LL
       b = b + TWO * a
       a = a
       IL(:,L) = a * w2 * IL(:,L+2) + b * IL(:,L+1)
    enddo
  end subroutine doP2IL

  subroutine doP2SL(LL,w2,IL)
    implicit none
    integer(IK), intent(in)    :: LL
    real(RK)   , intent(in)    :: w2
    real(RK)   , intent(inout) :: IL(:,:)
    ! *** end of interface ***

    real(RK)    :: a, b
    integer(IK) :: L

    ASSERT(1+LL+2<=size(IL,2))

    a =   ONE
    b = THREE
    IL(:,1) = a * w2 * IL(:,3) + b * IL(:,2)
    do L=2,1+LL
       b = b + TWO * a
       a = a
       IL(:,L) = a * w2 * IL(:,L+2) + b * IL(:,L+1)
    enddo
  end subroutine doP2SL

  !--------------- End of module -------------------------------------
end module shgi_rad

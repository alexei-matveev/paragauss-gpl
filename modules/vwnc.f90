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
module vwnc
!---------------------------------------------------------------
!
!  Purpose: contains routines calculating exchange correlation
!           potentials
!
!
!  Module called by: xc_hamiltonian, post_scf_module, response_module
!
!  References: J.Chem.Phys. Vol.98 No.7, April 1993 pp. 5612
!              (Erratum: JCP Vol.1 No.10, Nov. 15, pp. 9202)
!
!  Author: MS
!  Date: 1/96
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: HH
! Date:   11/97
! Description: Added subroutines "vwn_resp_calc","xalpha_resp_calc"
!              needed by the response_module.
!              Both subroutines do calculate the second functional
!              derivative of the correlation/exchange functional
!              with respect to the density.
!              In the following denote spin up by "a", spin down
!              by "b" and the correlation/exchange potential
!              by "v".
!              Then the *_resp subroutines calculate (depending on
!              user parameters) the following arrays:
!              spin=1:
!                dvdrho(:,1)= 0.5*dv_a/drho_a
!                dvdrho(:,2)= 0.5*dv_a/drho_b
!              spin=2:
!                dvdrho(:,1)= dv_a/drho_a
!                dvdrho(:,2)= dv_a/drho_b
!                dvdrho(:,3)= dv_b/drho_b
!                dvdrho(:,4)= dv_b/drho_a
!              Note:
!              Depending on input parameters for the response module
!              for spin=2 maybe the following calculation is done:
!              spin=2:
!                dvdrho(:,1)= dv_a/drho_a
!                dvdrho(:,2)= dv_a/drho_b
!                dvdrho(:,3)= calculated later in response_module
!                dvdrho(:,4)= calculated later in response_module
!                dvdrho(:,5)= dv_b/drho_b
!                dvdrho(:,6)= dv_b/drho_a
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author:      Uwe Birkenheuer
! Date:        8/98
! Description: subroutine xalpha_calc extended such that
!              a) second functional derivatives are available
!                 (with a new storage mode for dvdrho !)
!              b) the original offset technique can be turned off
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------

  use type_module
  use cpksdervs_matrices, only: cpks_noxc
  implicit none
  private         ! by default, all names are private
  save

!== Interrupt end of public interface of module =================


!------------ public functions and subroutines ------------------
public vwn_calc, vwn_resp_calc, xalpha_calc, xalpha_resp_calc,vwn_calcMDA,xalpha_calcMDA

public :: vwn_ldac!(nv, ispin, n, Fc, dFdn, dFdndn)

!================================================================
! End of public interface of module
!================================================================

!------------ Declaration of constants and variables ----
  integer(kind=i4_kind),private :: vec_length

  real(kind=r8_kind) :: fact, ffact, dfact, alpha, two13, qa, qf, qp,&
         XXPX0P,XXFX0f,XXAx0a,ddfzero

    real(kind=r8_kind),parameter :: zero=0.0_r8_kind,&
         one=1.0_r8_kind,&
         two=2.0_r8_kind,&
         three=3.0_r8_kind,&
         four=4.0_r8_kind,&
         five=5.0_r8_kind,&
         six=6.0_r8_kind,&
         eight=8.0_r8_kind,&
         nine=9.0_r8_kind,&
         sixteen=16.0_r8_kind,&
         thirtysix=36.0_r8_kind,&
         half=0.5_r8_kind,&
         deps=1.0E-50_r8_kind,&
         depsMDA=1.0E-8_r8_kind,&
         third=one/three,&
         f2third=two/three,&
         f4thrd=four/three,&
         beta=0.0042_r8_kind,&
         pi=3.141592653589793240_r8_kind,&
        !     fit parameters for interpolation
    aa=-0.0337740_r8_kind/two,&
         ba=1.13107_r8_kind,&
         ca=13.0045_r8_kind,&
         x0a=-0.0047584_r8_kind,&
         af=0.0310907_r8_kind/two,&
         bf=7.06042_r8_kind,&
         cf=18.0578_r8_kind,&
         x0f=-0.32500_r8_kind,&
         ap=0.0621841_r8_kind/two,&
         bp=3.72744_r8_kind,&
         cp=12.9352_r8_kind,&
         x0p=-0.10498_r8_kind



!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
contains

  function c_lda(na, nb) result(fxc)
    use type_module, only: RK => r8_kind
    use ad2x2
    implicit none
    type(ad), intent(in) :: na, nb
    type(ad) :: fxc
    ! *** end of interface ***

    type(ad) :: sums, X, ZETA, ECPARA, ECFERR, ALPHAC
    type(ad) :: BETAC, DELTEC, ECOR

    ! These eclipse module vars (rather constants):
    real(RK), parameter :: TWO13  = TWO**THIRD
    real(RK), parameter :: FACT   = THREE/(FOUR*PI)
    real(RK), parameter :: FFACT  = HALF/(TWO13 - ONE)
    real(RK), parameter :: DFACT  = F4THRD*FFACT
    real(RK), parameter :: ALPHA  = (FOUR/(NINE*PI))**THIRD
    real(RK), parameter :: DDFZERO = THIRD * DFACT*TWO

    ! FIXME: a  better way to  avoid floating point  exceptions? Older
    ! code  used depsMDA =  1.0e-8 at  this place  for some  time, see
    ! vwn_calcMDA():
    sums = na + nb + deps

    X = sqrt((fact/sums)**third)

    ZETA = (na - nb) / sums

    ECPARA = AP * ec(X, BP, CP, X0P)
    ECFERR = AF * ec(X, BF, CF, X0F)
    ALPHAC = AA * ec(X, BA, CA, X0A)

    BETAC  = (DDFZERO/ALPHAC)*(ECFERR-ECPARA)-ONE
    DELTEC = ALPHAC * (F(ZETA) / DDFZERO) * (ONE + BETAC * ZETA**4)

    ECOR = ECPARA + DELTEC

    fxc = ECOR * sums
  contains

    function ec(X, B, C, X0)
      type(ad), intent(in) :: X
      real(RK), intent(in) :: B, C, X0
      type(ad) :: ec
      ! *** end of interface ***

      real(RK) :: Q, XX0
      type(ad) :: XXX, FXX

      Q = SQRT(4 * C - B**2)
      XX0 = (X0 + B) * X0 + C
      XXX = (X  + B) * X  + C
      FXX = ATAN(Q / (TWO * X + B))

      ec = LOG(X**2 / XXX) + (2 * B / Q) * FXX &
           - (B * X0 / XX0) * (LOG((X - X0)**2 / XXX) + ((2 * (B + 2 * X0)) / Q) * FXX)
    end function ec

    function f(y)
      type(ad), intent(in) :: y
      type(ad) :: f
      f = ffact * ((ONE+Y)**F4THRD + (ONE-Y)**F4THRD - TWO)
    end function f
  end function c_lda

  subroutine vwn_ldac(nv, ispin, n, Fc, dFdn, dFdndn)
    use type_module, only: IK=>i4_kind, RK=>r8_kind
    use ad2x2, only: ad, var, fix, val, fst, snd, operator(*)
    implicit none
    integer(IK), intent(in) :: nv ! vector length
    integer(IK), intent(in) :: ispin ! 1 or 2
    real(RK), dimension(:,:), intent(in) :: n ! electron density n(nv,1) or n(nv,2)
    real(RK), dimension(:), intent(inout) :: Fc ! c-functional Fc(nv)
    real(RK), dimension(:,:), intent(inout) :: dFdn !dF/dn dFdn(nv,1) or dFdn(nv,2)
    real(RK), dimension(:,:), intent(inout), optional :: dFdndn
    ! *** end of interface ***

    !
    ! VN  variant  (will  be  implemented  and  used  everywhere,  see
    ! e.g. exchange.f90).  NOTE for second derivatives packing:
    !
    !              1    2    3
    ! dfdg        AA   BB   AB
    ! dfdndn      AA   BB   AB
    !
    enum, bind(c)
       enumerator :: A = 1, B
       enumerator :: AA = 1, BB, AB
    end enum

    type(ad) :: f, na, nb, nn
    integer :: i

    select case(ispin)

    case (1)
       ! ispin == 1 rhoa = rhob = rho/2, rho as a single independent
       ! variable:
       do i = 1, nv
          nn = var(1, n(i, 1))

          f = c_lda(half * nn, half * nn)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, 1) = dFdn(i, 1) + fst(1, f)

          if (present(dFdndn)) then
             dFdndn(i, 1) = dFdndn(i, 1) + snd(1, 1, f)
          endif
       enddo

    case (2)
       ! ispin == 2 general case, rhoa /= rhob as two independent
       ! variables:
       do i = 1, nv
          na = var(1, n(i, A))
          nb = var(2, n(i, B))

          f = c_lda(na, nb)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, A) = dFdn(i, A) + fst(1, f)
          dFdn(i, B) = dFdn(i, B) + fst(2, f)

          if (present(dFdndn)) then
             dFdndn(i, AA) = dFdndn(i, AA) + snd(1, 1, f)
             dFdndn(i, BB) = dFdndn(i, BB) + snd(2, 2, f)
             dFdndn(i, AB) = dFdndn(i, AB) + snd(1, 2, f)
          endif
       enddo

    case default
       print *,'vwn_calc1: no such case: ispin=', ispin
       stop
    end select
  end subroutine vwn_ldac

  subroutine vwn_calc(rho,dfdrho,ispin,fxc,vec_length_act,eps)
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dfdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length_act
    real(kind=r8_kind),intent(inout),optional :: eps(:)
    !** End of interface *****************************************

    real(kind=r8_kind),dimension(vec_length_act)  :: ecpara,dep, alphac, ecor,sums,x
    real(kind=r8_kind),dimension(vec_length_act)  :: ecferr,def,dea,zeta,fe,ep,em,betac,zeta3,zeta4
    real(kind=r8_kind),dimension(vec_length_act)  :: deltec,u1,u2,u3,u4,u5,u6,u7,u8

    vec_length=vec_length_act

    TWO13  = TWO**THIRD

    QA = SQRT(FOUR*CA-BA*BA)
    QF = SQRT(FOUR*CF-BF*BF)
    QP = SQRT(FOUR*CP-BP*BP)

    FACT   = THREE/(FOUR*PI)
    FFACT  = HALF/(TWO13 - ONE)
    DFACT  = F4THRD*FFACT
    ALPHA  = (FOUR/(NINE*PI))**THIRD

    XXPX0P=x0p*x0p+bp*x0p+cp
    XXFX0f=x0f*x0f+bf*x0f+cf
    XXAx0a = x0a*x0a + BA*x0a + CA
    ddfzero= THIRD * DFACT*two
    !
    ! calculation starts here
    !
    if(ispin==1) then  ! non spin-polarized
       x=sqrt((fact/(rho(1:vec_length,1)+deps))**third)

       ECPARA = AP*(LOG(X*X/XXP(X))+(TWO*BP/QP)*FXP(X)-(BP*X0P/XXPX0P)*&
            (LOG((X-X0P)**TWO/XXP(X))+&
            ((TWO*(BP+TWO*X0P))/QP)*FXP(X)))

       DEP = AP*(TWO/X-(ONE-(BP*X0P)/XXPX0P)*(TWO*X+BP)/XXP(X)-&
            (TWO*BP*X0P)/((X-X0P)*XXPX0P)-FOUR*BP*(ONE-((BP+TWO*X0P)/&
            XXPX0P)*X0P)*(ONE/(QP*QP+(TWO*X+BP)**TWO)))

       dfdrho(1:vec_length,1)=dfdrho(1:vec_length,1)+ecpara-(x/six)*dep

       if(present(eps)) eps(1:vec_length)=eps(1:vec_length)+ecpara

       fxc(1:vec_length)=fxc(1:vec_length)+ecpara*rho(1:vec_length,1)

    else   ! spin-polarized
       sums=rho(1:vec_length,1)+rho(1:vec_length,2)+deps

       x=sqrt((fact/sums)**third)


       ZETA = ( RHO(1:vec_length,1)  - RHO(1:vec_length,2)  )/sums

       ECPARA = AP*(LOG(X*X/XXP(X))+(TWO*BP/QP)*FXP(X)-(BP*X0P/XXPX0P)*&
            (LOG((X-X0P)**TWO/XXP(X))+((TWO*(BP+TWO*X0P))/QP)*FXP(X)))

       ECFERR = AF*(LOG(X*X/XXF(X))+(TWO*BF/QF)*FXF(X)-(BF*X0F/XXFX0F)*&
            (LOG((X-X0F)**TWO/XXF(X))+((TWO*(BF+TWO*X0F))/QF)*FXF(X)))

       ALPHAC = AA*(LOG(X*X/XXA(X))+(TWO*BA/QA)*FXA(X)-(BA*X0A/XXAX0A)*&
            (LOG((X-X0A)**TWO/XXA(X))+((TWO*(BA+TWO*X0A))/QA)*FXA(X)))

       DEP = AP*(TWO/X-(ONE-(BP*X0P)/XXPX0P)*(TWO*X+BP)/XXP(X)-&
            (TWO*BP*X0P)/((X-X0P)*XXPX0P)-FOUR*BP*(ONE-((BP+TWO*X0P)/&
            XXPX0P)*X0P)*(ONE/(QP*QP+(TWO*X+BP)**TWO)))

       DEF = AF*(TWO/X-(ONE-(BF*X0F)/XXFX0F)*(TWO*X+BF)/XXF(X)-&
            (TWO*BF*X0F)/((X-X0F)*XXFX0F)-FOUR*BF*(ONE-((BF+TWO*X0F)/&
            XXFX0F)*X0F)*(ONE/(QF*QF+(TWO*X+BF)**TWO)))

       DEA = AA*(TWO/X-(ONE-(BA*X0A)/XXAX0A)*(TWO*X+BA)/XXA(X)-&
            (TWO*BA*X0A)/((X-X0A)*XXAX0A)-FOUR*BA*(ONE-((BA+TWO*X0A)/&
            XXAX0A)*X0A)*(ONE/(QA*QA+(TWO*X+BA)**TWO)))

       FE     = F(ZETA)
       EP     = ONE + ZETA
       EM     = ONE - ZETA
       BETAC  = (DDFZERO/ALPHAC)*(ECFERR-ECPARA)-ONE
       ZETA3  = ZETA*ZETA*ZETA
       ZETA4  = ZETA*ZETA3
       DELTEC = ALPHAC*(FE/DDFZERO)*(ONE+BETAC*(ZETA4))

       U1 = (ONE - ZETA4 * FE) * DEP
       U2 = ZETA4 * FE * DEF
       U3 = (ONE - ZETA4) * FE/DDFZERO * DEA
       U4 = + FOUR * (ONE - ZETA) * ALPHAC/DDFZERO
       U5 = ZETA3 * FE * BETAC
       U6 = ONE + BETAC * ZETA4
       U7 = DF(ZETA)/FOUR
       U8 = - FOUR * (ONE + ZETA) * ALPHAC/DDFZERO

       ECOR = ECPARA + DELTEC

       dfdrho(1:vec_length,1) =  dfdrho(1:vec_length,1)+ ECOR&
       - (X/SIX)*(U1 + U2 + U3) + U4*(U5 + U6*U7)

       dfdrho(1:vec_length,2) = dfdrho(1:vec_length,2) + &
         ECOR - (X/SIX)*(U1 + U2 + U3) + U8*(U5 + U6*U7)

       if(present(eps)) eps(1:vec_length)=ecor

       fxc(1:vec_length)= fxc(1:vec_length) + ECOR*sums

    end if
! no longer needed ?    u8=xxp(x)

  end subroutine vwn_calc

  subroutine vwn_calcMDA(rho,dfdrho,ispin,fxc,vec_length_act,eps,dvdrho)
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dfdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length_act
    real(kind=r8_kind),intent(inout),optional :: eps(:)
    real(kind=r8_kind),intent(inout),optional :: dvdrho(:,:)
    !** End of interface *****************************************

    real(kind=r8_kind),dimension(vec_length_act)  :: ecpara,ecferr,dep,def,dea,&
         zeta,alphac,fe,ep,em,betac,zeta3,zeta4,deltec,u1,u2,u3,u4,u5,u6,u7,u8,&
         ecor,sums,x,ddep
    real(kind=r8_kind),dimension(vec_length_act)  :: ddef,ddea,DBETAC, &
                      u1_z,u1_up,u1_dn,u2_z,u2_up,u2_dn,u3_z,u3_up,u3_dn,      &
                      u4_z,u4_up,u4_dn,u5_z,u5_up,u5_dn,u6_z,u6_up,u6_dn,      &
                      u7_z,u7_up,u7_dn,u8_z,u8_up,u8_dn,x_r,dfe


    vec_length=vec_length_act

    TWO13  = TWO**THIRD

    QA = SQRT(FOUR*CA-BA*BA)
    QF = SQRT(FOUR*CF-BF*BF)
    QP = SQRT(FOUR*CP-BP*BP)

    FACT   = THREE/(FOUR*PI)
    FFACT  = HALF/(TWO13 - ONE)
    DFACT  = F4THRD*FFACT
    ALPHA  = (FOUR/(NINE*PI))**THIRD

    XXPX0P=x0p*x0p+bp*x0p+cp
    XXFX0f=x0f*x0f+bf*x0f+cf
    XXAx0a = x0a*x0a + BA*x0a + CA
    ddfzero= THIRD * DFACT*two
    !
    ! calculation starts here
    !
       sums=sum(rho(:vec_length,:),2)+depsMDA
       x=sqrt((fact/sums)**third); x_r=-(x/six)/sums

       ECPARA = AP*(LOG(X*X/XXP(X))+(TWO*BP/QP)*FXP(X)-(BP*X0P/XXPX0P)*&
            (LOG((X-X0P)**TWO/XXP(X))+((TWO*(BP+TWO*X0P))/QP)*FXP(X)))

       DEP = AP*(TWO/X-(ONE-(BP*X0P)/XXPX0P)*(TWO*X+BP)/XXP(X)-&
            (TWO*BP*X0P)/((X-X0P)*XXPX0P)-FOUR*BP*(ONE-((BP+TWO*X0P)/&
            XXPX0P)*X0P)*(ONE/(QP*QP+(TWO*X+BP)**TWO)))


      if(ispin==1) then
       dfdrho(:vec_length,1)=dfdrho(:vec_length,1)+ecpara-(x/six)*dep

       if(present(eps)) eps(1:vec_length)=eps(1:vec_length)+ecpara

       fxc(:vec_length)=fxc(:vec_length)+ecpara*sums
      else


       ZETA = ( RHO(:vec_length,1)  - RHO(:vec_length,2)  )/sums

       ECFERR = AF*( LOG(X*X/XXF(X))&
                    +(TWO*BF/QF)*FXF(X) &
                    -(BF*X0F/XXFX0F)*&
                     ( LOG((X-X0F)**TWO/XXF(X)) + ((TWO*(BF+TWO*X0F))/QF)*FXF(X) ))
       DEF = AF*( TWO/X- &
                  (ONE-(BF*X0F)/XXFX0F)*(TWO*X+BF)/XXF(X)-&
                  TWO*BF*X0F/((X-X0F)*XXFX0F)-&
                  FOUR*BF*(ONE-((BF+TWO*X0F)/XXFX0F)*X0F)*(ONE/(QF*QF+(TWO*X+BF)**TWO)) )

       ALPHAC = AA*( LOG(X*X/XXA(X))    &
                    +(TWO*BA/QA)*FXA(X)-&
                     (BA*X0A/XXAX0A)*&
                     ( LOG((X-X0A)**TWO/XXA(X))+((TWO*(BA+TWO*X0A))/QA)*FXA(X) ))



       DEA = AA*( TWO/X- &
                  (ONE-(BA*X0A)/XXAX0A)*(TWO*X+BA)/XXA(X)-&
                  (TWO*BA*X0A)/((X-X0A)*XXAX0A)- &
                  FOUR*BA*(ONE-((BA+TWO*X0A)/XXAX0A)*X0A)*(ONE/(QA*QA+(TWO*X+BA)**TWO)))


       FE     = F(ZETA)
       EP     = ONE + ZETA
       EM     = ONE - ZETA

       BETAC  =  (DDFZERO/ALPHAC)*(ECFERR-ECPARA)-ONE


       ZETA3  = ZETA*ZETA*ZETA
       ZETA4  = ZETA*ZETA3
       DELTEC = ALPHAC*(FE/DDFZERO)*(ONE+BETAC*(ZETA4))

       U1 = (ONE - ZETA4 * FE) * DEP

       DFE=DF(ZETA)

       U2 = ZETA4 * FE * DEF
       U3 = (ONE - ZETA4) * FE/DDFZERO * DEA
       U4 = + FOUR * (ONE - ZETA) * ALPHAC/DDFZERO
       U5 = BETAC * ZETA3 * FE
       U6 = ONE + BETAC * ZETA4
       U7 = DF(ZETA)/FOUR

       U8 = - FOUR * (ONE + ZETA) * ALPHAC/DDFZERO

       ECOR = ECPARA + DELTEC

       dfdrho(:vec_length,1)=dfdrho(:vec_length,1) &
       +ECOR - (X/SIX)*(U1 + U2 + U3) + (U5 + U6*U7)*U4

       dfdrho(:vec_length,2)=dfdrho(:vec_length,2) &
       +ECOR - (X/SIX)*(U1 + U2 + U3) + (U5 + U6*U7)*U8

       if(present(eps)) eps(:vec_length)=ecor

       fxc(:vec_length)= fxc(:vec_length) + ECOR*sums
      endif


  if(present(dvdrho)) then
             ddep = ap*(-two/(x*x)-(one-(bp*x0p/xxpx0p))*(two/xxp(x))&
                  &+two*bp*x0p/(xxpx0p*(x-x0p)**two)&
                  &+(one-(bp*x0p/xxpx0p))*((two*x+bp)**two/(xxp(x))**two)&
                  &+sixteen*bp*(two*x+bp)*(one-(x0p*(two*x0p+bp)/xxpx0p))&
                  &/(((two*x+bp)**two+qp*qp)**two))
    if(ispin==1) then  ! non spin-polarized
             dvdrho(:vec_length,1)=dvdrho(:vec_length,1) &
               +(-(x*five/thirtysix)*dep+(x*x/thirtysix)*ddep)/sums
             ! For spin unpol. case we need dv/drho, but this equals
             !    dv/dvrho = 0.5*[dv_a/drho_a + dv_a/drho_b]
             ! evaluated at rho_a=rho_b=0.5*rho.

!             ! 0.5*dv_a/drho_a
!             dvdrho(i_vec,1)=dvdrho(i_vec,1)+help_vec+half*alphac/(rho(i_vec,1))
!             ! 0.5*dv_a/drho_b
!             dvdrho(i_vec,2)=dvdrho(i_vec,2)+help_vec-half*alphac/(rho(i_vec,1))


    else   ! spin-polarized
       DDEF = AF*(-TWO/X**2 &
                  -(ONE-(BF*X0F)/XXFX0F)*TWO/XXF(X)&
                  -(ONE-(BF*X0F)/XXFX0F)*(TWO*X+BF)/XXF(X)**2*DXXF(X)&
                  -TWO*BF*X0F/XXFX0F/(X-X0F)**2&
                  -FOUR*BF* (ONE-((BF+TWO*X0F)/XXFX0F) * X0F )*(ONE/(QF*QF+(X+X+BF)**2)**2)*(X+X+BF)*FOUR )
       DDEA = AA*(-TWO/X**2 &
                  -(ONE-(BA*X0A)/XXAX0A)*TWO/XXA(X)&
                  -(ONE-(BA*X0A)/XXAX0A)*(TWO*X+BA)/XXA(X)**2*DXXA(X)&
                  -TWO*BA*X0A/XXAX0A/(X-X0A)**2&
                  -FOUR*BA* (ONE-((BA+TWO*X0A)/XXAX0A) * X0A )*(ONE/(QA*QA+(X+X+BA)**2)**2)*(X+X+BA)*FOUR )

       DBETAC = -(DDFZERO/ALPHAC**2)*(ECFERR-ECPARA)*DEA &
                +(DDFZERO/ALPHAC)*(def-dep)

       U1_Z=-(four*ZETA3*FE+ZETA4*DFE)*DEP
       U1_UP= (ONE - ZETA4 * FE) * DDEP*x_r + U1_Z*EM/sums
       U1_DN= (ONE - ZETA4 * FE) * DDEP*x_r - U1_Z*EP/sums

       U2_Z = (FOUR*ZETA3*FE+ZETA4 *DFE) * DEF
       U2_UP = ZETA4 * FE * DDEF*x_r + U2_Z*EM/sums
       U2_DN = ZETA4 * FE * DDEF*x_r - U2_Z*EP/sums

       U3_Z = ( (ONE - ZETA4) *DFE-four*ZETA3*FE )/DDFZERO * DEA
       U3_UP = (ONE - ZETA4) * FE/DDFZERO * ddea*x_r + U3_Z*EM/sums
       U3_DN = (ONE - ZETA4) * FE/DDFZERO * ddea*x_r - U3_Z*EP/sums

       U4_Z = - FOUR * ALPHAC/DDFZERO
       U4_UP= FOUR * (ONE - ZETA) * DEA/DDFZERO*x_r + U4_Z*EM/sums
       U4_DN= FOUR * (ONE - ZETA) * DEA/DDFZERO*x_r - U4_Z*EP/sums

       U5_Z = (three*ZETA**2*FE+ZETA3 *DFE)*BETAC
       U5_UP=   U5_Z*EM/sums + ZETA3 * FE * DBETAC * x_r
       U5_DN= - U5_Z*EP/sums + ZETA3 * FE * DBETAC * x_r

       U6_Z = four*ZETA3*BETAC
       U6_UP= DBETAC * ZETA4*x_r + U6_Z*EM/sums
       U6_DN= DBETAC * ZETA4*x_r - U6_Z*EP/sums

       U7_Z = DDF(ZETA)/FOUR
       U7_UP= U7_Z*EM/sums
       U7_DN=-U7_Z*EP/sums

       U8_Z= - FOUR * ALPHAC/DDFZERO
       U8_UP=- FOUR * (ONE + ZETA) *DEA/DDFZERO*x_r+U8_Z*EM/sums
       U8_DN=- FOUR * (ONE + ZETA) *DEA/DDFZERO*x_r-U8_Z*EP/sums
!!!       dfdrho(:vec_length,1)=dfdrho(:vec_length,1)+ECOR-(X/SIX)*(U1+U2+U3)  &
!!!                           + U4*(U5 + U6*U7)
       dvdrho(:vec_length,1)=dvdrho(:vec_length,1) &
                            -((X/SIX)*(U1+U2+U3)-U4*(U5+U6*U7))/sums  &
                            -(x_r/SIX)*(U1+U2+U3)  &
                            -(X/SIX)*(U1_UP+U2_UP+U3_UP) &
                            +(U5_UP+U6_UP*U7+U6*U7_UP)*U4+U4_UP*(U5 + U6*U7)
!!!       dfdrho(:vec_length,2)=dfdrho(:vec_length,2)+ECOR &
!!!       - (X/SIX)*(U1 + U2 + U3) + U8*(U5 + U6*U7)
       dvdrho(:vec_length,2)=dvdrho(:vec_length,2) &
                            -((X/SIX)*(U1+U2+U3)-U8*(U5+U6*U7) )/sums   &
                            -(x_r/SIX)*(U1+U2+U3)-(X/SIX)*(U1_DN+U2_DN+U3_DN) &
                            + U8_DN*(U5 + U6*U7)+U8*(U5_DN+U6_DN*U7+U6*U7_DN)
       dvdrho(:vec_length,3)=dvdrho(:vec_length,3) &
                            -((X/SIX)*(U1+U2+U3)-U8*(U5+U6*U7) )/sums  &
                            -(x_r/SIX)*(U1+U2+U3) &
                            -(X/SIX)*(U1_DN+U2_DN+U3_DN) &
                            +(U5_DN+U6_DN*U7+U6*U7_DN)*U4+U4_DN*(U5 + U6*U7)

    endif
   endif
! no longer needed ?    u8=xxp(x)

  end subroutine vwn_calcMDA


  subroutine vwn_resp_calc(rho,ispin,vec_length_act,resp_deps,dvdrho,fxc)
    use error_module
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dvdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length_act
    real(kind=r8_kind),intent(in) :: resp_deps
    real(kind=r8_kind),dimension(:),optional,intent(inout) :: fxc
    !** End of interface *****************************************

    integer(kind=i4_kind) :: i_vec
    real(kind=r8_kind)  :: ecpara, dep, alphac, x, ddep, help_vec
!!    real(kind=r8_kind),dimension(vec_length_act)  :: ecferr,def,dea,zeta,fe,ep,em,betac,zeta3,zeta4
!!    real(kind=r8_kind),dimension(vec_length_act)  :: deltec,u1,u2,u3,u4,u5,u6,u7,u8,ecor,sums
!!    real(kind=r8_kind),dimension(vec_length_act)  :: ddea,ddef,dbetac,ddbetac,aux1,aux2,aux3,aux4
!!    real(kind=r8_kind),dimension(vec_length_act)  :: dfe,ddfe

    vec_length=vec_length_act

    TWO13  = TWO**THIRD

    QA = SQRT(FOUR*CA-BA*BA)
    QF = SQRT(FOUR*CF-BF*BF)
    QP = SQRT(FOUR*CP-BP*BP)

    FACT   = THREE/(FOUR*PI)
    FFACT  = HALF/(TWO13 - ONE)
    DFACT  = F4THRD*FFACT
    ALPHA  = (FOUR/(NINE*PI))**THIRD

    XXPX0P=x0p*x0p+bp*x0p+cp
    XXFX0f=x0f*x0f+bf*x0f+cf
    XXAx0a = x0a*x0a + BA*x0a + CA
    ddfzero= THIRD * DFACT*two
    !
    ! calculation starts here
    !

    if(ispin==1) then  ! non spin-polarized

       do i_vec=1, vec_length

          if( abs(rho(i_vec,1)) > resp_deps ) then

             x=sqrt((fact/(rho(i_vec,1)))**third)

             ECPARA = AP*(LOG(X*X/s_XXP(X))+(TWO*BP/QP)*s_FXP(X)-(BP*X0P/XXPX0P)*&
                  (LOG((X-X0P)**TWO/s_XXP(X))+&
                  ((TWO*(BP+TWO*X0P))/QP)*s_FXP(X)))

             DEP = AP*(TWO/X-(ONE-(BP*X0P)/XXPX0P)*(TWO*X+BP)/s_XXP(X)-&
                  (TWO*BP*X0P)/((X-X0P)*XXPX0P)-FOUR*BP*(ONE-((BP+TWO*X0P)/&
                  XXPX0P)*X0P)*(ONE/(QP*QP+(TWO*X+BP)**TWO)))

             ! calculation of derivatives dvdrho of vwn corr. potential starts here
             ALPHAC = AA*(LOG(X*X/s_XXA(X))+(TWO*BA/QA)*s_FXA(X)-(BA*X0A/XXAX0A)*&
                  (LOG((X-X0A)**TWO/s_XXA(X))+((TWO*(BA+TWO*X0A))/QA)*s_FXA(X)))

             ddep = ap*(-two/(x*x)-(one-(bp*x0p/xxpx0p))*(two/s_xxp(x))&
                  &+two*bp*x0p/(xxpx0p*(x-x0p)**two)&
                  &+(one-(bp*x0p/xxpx0p))*((two*x+bp)**two/(s_xxp(x))**two)&
                  &+sixteen*bp*(two*x+bp)*(one-(x0p*(two*x0p+bp)/xxpx0p))&
                  &/(((two*x+bp)**two+qp*qp)**two))

             help_vec=half*(-(x*five/thirtysix)*dep&
                  &+(x*x/thirtysix)*ddep)/(rho(i_vec,1))
             ! For spin unpol. case we need dv/drho, but this equals
             !    dv/dvrho = [dv_a/drho_a + dv_a/drho_b]
             ! evaluated at rho_a=rho_b=0.5*rho.
             ! dv_a/drho_a
             dvdrho(i_vec,1)=dvdrho(i_vec,1)+help_vec&
                  &+alphac/(rho(i_vec,1))
             ! dv_a/drho_b
             dvdrho(i_vec,2)=dvdrho(i_vec,2)+help_vec&
                  &-alphac/(rho(i_vec,1))
             ! Note: the elements
             !  dvdrho(:,j), j=3,4
             ! are calculated later when the contributions of all
             ! functionals to dvdrho have been collected.

             ! optional: evaluate the energy functional itself
             ! (may be used for debugging in grid module)
             if(present(fxc)) fxc(i_vec)=fxc(i_vec)+ecpara*rho(i_vec,1)

          end if

       end do

    else   ! spin-polarized
      call error_handler ("NO VWN for OPEN SHELL, PLEASE ADD VWN = FALSE IN RESPONSE FIELD")
      stop
    end if
! no longer needed ?     u8=xxp(x)

  end subroutine vwn_resp_calc


  subroutine xalpha_calc(rho,dfdrho,ispin,fxc,vec_length,eps,dvdrho,offset)

    real(kind=r8_kind),dimension(:,:),intent(in) :: rho
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dfdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    real(kind=r8_kind),intent(inout),dimension(:) :: fxc
    real(kind=r8_kind),intent(inout),dimension(:),optional :: eps
    logical           , intent(in   ), optional :: offset ! default = .TRUE.
    real(kind=r8_kind), intent(inout), optional :: dvdrho(:,:)
    ! spin-symmetric case
    !   d/drho(r) d/drho(r') E_xc[rho] = Z_xc(r) delta(r-r')
    !   dvdrho(:,1) = Z_xc(r_grid)
    ! spin-polarized case
    !   d/drho_s(r) d/drho_s'(r') E_xc[rho] = Z_xc,s,s'(r) delta(r-r')
    !   dvdrho(:,1) = Z_xc,up,up(r_grid)
    !   dvdrho(:,2) = Z_xc,dn,dn(r_grid)
    !   dvdrho(:,3) = Z_xc,up,dn(r_grid) = Z_xc,dn,up(r_grid)
    !** End of interface *****************************************
    real(kind=r8_kind),dimension(vec_length,ispin) :: tmp
    real(kind=r8_kind) :: factor
    logical :: use_deps

    if (present(offset)) then
       use_deps = offset
    else
       use_deps = .true.
    endif

    if(ispin==1) then  ! non spinpolarized case
!..............................................................................
!      for alpha = 2/3 :
!      f_xc   = - 3/4 * (3/pi)^1/3 * rho^4/3
!      eps_xc = - 3/4 * (3/pi)^1/3 * rho^1/3
!      V_xc   = -       (3/pi)^1/3 * rho^1/3
!      Z_xc   = - 1/3 * (3/pi)^1/3 / rho^2/3
!..............................................................................
       factor = - (3.0_r8_kind/pi)**third
       if (use_deps) then
       tmp(:,1) = factor * ( rho(1:vec_length,1) + deps )**third
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1) = dvdrho(1:vec_length,1) + &
                  third * tmp(:,1) / ( rho(1:vec_length,1) + deps )
          endif
       else
          tmp(:,1) = factor * exp( log( rho(1:vec_length,1) ) * third )
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1) = dvdrho(1:vec_length,1) + &
                                      third * tmp(:,1) / rho(1:vec_length,1)
          endif
       endif
       dfdrho(1:vec_length,1) = dfdrho(1:vec_length,1) + tmp(:,1)
       tmp(:,1) = 0.75_r8_kind * tmp(:,1)
       fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)*rho(1:vec_length,1)
!       if(cpks_noxc) then
!        dfdrho(1:vec_length,1) =0.0_r8_kind ! cpks checks
!        fxc(1:vec_length) = 0.0_r8_kind
!       endif
       if (present(eps)) then
          eps(1:vec_length) = eps(1:vec_length) + tmp(:,1)
       endif

    else              ! spinpolarized case
!..............................................................................
!      for alpha = 2/3 :
!      f_xc    = - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^4/3
!      eps_xc  = - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^4/3 / Sum(s) rho_s
!        and not - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^1/3 !
!      V_xc,s  = -       (6/pi)^1/3 * rho_s^1/3
!      Z_xc,st = - 1/3 * (6/pi)^1/3 / rho_s^2/3 delta_st
!..............................................................................
       factor = - (6.0_r8_kind/pi)**third
       if (use_deps) then
       tmp = factor * ( rho(1:vec_length,:) + deps )**third
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1:2) = dvdrho(1:vec_length,1:2) + &
                  third * tmp / ( rho(1:vec_length,:) + deps )
          endif
       else
          tmp = factor * exp( log( rho(1:vec_length,:) ) * third )
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1:2) = dvdrho(1:vec_length,1:2) + &
                                        third * tmp / rho(1:vec_length,:)
          endif
       endif
       dfdrho(1:vec_length,:) = dfdrho(1:vec_length,:) + tmp
       tmp = tmp * rho(1:vec_length,:)
       tmp(:,1) = 0.75_r8_kind * ( tmp(:,1) + tmp(:,2) )
       fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)
!:BUG  ! eps(rho) is  n o t  factor * Sum(s) rho_s^(1/3) !
       if (present(eps)) then
          if (use_deps) then
          tmp(:,2) = rho(1:vec_length,1) + rho(1:vec_length,2) + deps
          else
             tmp(:,2) = rho(1:vec_length,1) + rho(1:vec_length,2)
          endif
          eps(1:vec_length) = eps(1:vec_length) + tmp(:,1) / tmp(:,2)
       endif
    endif

  end subroutine xalpha_calc

  subroutine xalpha_calcMDA(rho,dfdrho,ispin,fxc,vec_length,eps,dvdrho,offset)

    real(kind=r8_kind),dimension(:,:),intent(in) :: rho
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dfdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    real(kind=r8_kind),intent(inout),dimension(:) :: fxc
    real(kind=r8_kind),intent(inout),dimension(:),optional :: eps
    logical           , intent(in   ), optional :: offset ! default = .TRUE.
    real(kind=r8_kind), intent(inout), optional :: dvdrho(:,:)
    ! spin-symmetric case
    !   d/drho(r) d/drho(r') E_xc[rho] = Z_xc(r) delta(r-r')
    !   dvdrho(:,1) = Z_xc(r_grid)
    ! spin-polarized case
    !   d/drho_s(r) d/drho_s'(r') E_xc[rho] = Z_xc,s,s'(r) delta(r-r')
    !   dvdrho(:,1) = Z_xc,up,up(r_grid)
    !   dvdrho(:,2) = Z_xc,dn,dn(r_grid)
    !   dvdrho(:,3) = Z_xc,up,dn(r_grid) = Z_xc,dn,up(r_grid)
    !** End of interface *****************************************
    real(kind=r8_kind),dimension(vec_length,ispin) :: tmp
    real(kind=r8_kind) :: factor
    logical :: use_deps

    if (present(offset)) then
       use_deps = offset
    else
       use_deps = .true.
    endif

    if(ispin==1) then  ! non spinpolarized case
!..............................................................................
!      for alpha = 2/3 :
!      f_xc   = - 3/4 * (3/pi)^1/3 * rho^4/3
!      eps_xc = - 3/4 * (3/pi)^1/3 * rho^1/3
!      V_xc   = -       (3/pi)^1/3 * rho^1/3
!      Z_xc   = - 1/3 * (3/pi)^1/3 / rho^2/3
!..............................................................................
       factor = - (3.0_r8_kind/pi)**third
       if (use_deps) then
       tmp(:,1) = factor * ( rho(1:vec_length,1) + depsMDA )**third
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1) = dvdrho(1:vec_length,1) + &
                  third * tmp(:,1) / ( rho(1:vec_length,1) + depsMDA )
          endif
       else
          tmp(:,1) = factor * exp( log( rho(1:vec_length,1) ) * third )
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1) = dvdrho(1:vec_length,1) + &
                   third * tmp(:,1) / (rho(1:vec_length,1) + depsMDA )
          endif
       endif
       dfdrho(1:vec_length,1) = dfdrho(1:vec_length,1) + tmp(:,1)
       tmp(:,1) = 0.75_r8_kind * tmp(:,1)
       fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)*rho(1:vec_length,1)

       if(cpks_noxc)  then
        fxc=0.0_r8_kind ! cpks checks
        dfdrho(1:vec_length,1) =0.0_r8_kind ! cpks checks
       if(present(dvdrho)) &
        dvdrho(1:vec_length,1) =0.0_r8_kind ! cpks checks
       endif

       if (present(eps)) then
          eps(1:vec_length) = eps(1:vec_length) + tmp(:,1)
       endif

    else              ! i.e. spinpolarized case
       if(cpks_noxc)  then
        fxc=0.0_r8_kind ! cpks checks
        dfdrho =0.0_r8_kind ! cpks checks
       if(present(dvdrho)) &
        dvdrho=0.0_r8_kind ! cpks checks
       else
!..............................................................................
!      for alpha = 2/3 :
!      f_xc    = - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^4/3
!      eps_xc  = - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^4/3 / Sum(s) rho_s
!        and not - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^1/3 !
!      V_xc,s  = -       (6/pi)^1/3 * rho_s^1/3
!      Z_xc,st = - 1/3 * (6/pi)^1/3 / rho_s^2/3 delta_st
!..............................................................................
       factor = - (6.0_r8_kind/pi)**third
       if (use_deps) then
       tmp = factor * ( rho(1:vec_length,:) + depsMDA )**third
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1:2) = dvdrho(1:vec_length,1:2) + &
                  third * tmp / ( rho(1:vec_length,:) + depsMDA )
          endif
       else
          tmp = factor * exp( log( rho(1:vec_length,:) ) * third )
          if (present(dvdrho)) then
             dvdrho(1:vec_length,1:2) = dvdrho(1:vec_length,1:2) + &
                                        third * tmp / rho(1:vec_length,:)
          endif
       endif
       dfdrho(1:vec_length,:) = dfdrho(1:vec_length,:) + tmp
       tmp = tmp * rho(1:vec_length,:)
       tmp(:,1) = 0.75_r8_kind * ( tmp(:,1) + tmp(:,2) )
       fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)
!:BUG  ! eps(rho) is  n o t  factor * Sum(s) rho_s^(1/3) !
       if (present(eps)) then
          if (use_deps) then
          tmp(:,2) = rho(1:vec_length,1) + rho(1:vec_length,2) + depsMDA
          else
             tmp(:,2) = rho(1:vec_length,1) + rho(1:vec_length,2)
          endif
          eps(1:vec_length) = eps(1:vec_length) + tmp(:,1) / tmp(:,2)
       endif

      endif
    endif

  end subroutine xalpha_calcMDA


  subroutine xalpha_resp_calc(rho,ispin,vec_length,resp_deps,dvdrho,fxc)

    real(kind=r8_kind),dimension(:,:),intent(in) :: rho
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dvdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    real(kind=r8_kind),intent(in) :: resp_deps
    real(kind=r8_kind),intent(inout),dimension(:),optional :: fxc
    !** End of interface *****************************************
    real(kind=r8_kind),dimension(vec_length,ispin) :: tmp
    real(kind=r8_kind) :: factor
!!    real(kind=r8_kind) :: temp
!!    integer(kind=i4_kind) :: i_vec

    factor = -((three/pi)**third)/three

    if(ispin==1) then  ! non spinpolarized case
      tmp(1:vec_length,1) = factor/(((rho(1:vec_length,1)/2.0_r8_kind) + deps)**f2third)
      dvdrho(1:vec_length,1) = dvdrho(1:vec_length,1) + tmp(1:vec_length,1)
      ! note: dv_a/drho_b = 0 for Slater functional
      ! => dvdrho(1:vec_length,2) = dvdrho(1:vec_length,2) + 0

      ! optional: evaluate the energy functional itself
      ! (may be used for debugging in grid module)
      if(present(fxc)) then
        factor = -(three/four)*((three/pi)**third)
        tmp(1:vec_length,1) = factor * ( rho(1:vec_length,1))**f4thrd
        fxc = fxc + tmp(1:vec_length,1)
      end if

    else              ! i.e. spinpolarized case

      tmp = factor / ( ( rho(1:vec_length,:) + deps )**f2third)
      ! dv_a/drho_a
      dvdrho(1:vec_length,1) = dvdrho(1:vec_length,1) + tmp(:,1)
      ! dv_b/drho_b
      if(size(dvdrho,2)==4) then
         ! no additional matrix elements calculated by response_module
         dvdrho(1:vec_length,3) = dvdrho(1:vec_length,3) + tmp(:,2)
         ! note: dv_a/drho_b=0
         ! => dvdrho(1:vec_length,4) = dvdrho(1:vec_length,4) + 0
      else
         ! dv_b/drho_b
         dvdrho(1:vec_length,5) = dvdrho(1:vec_length,5) + tmp(:,2)
         ! note:  dv_b/drho_a=0
         ! => dvdrho(1:vec_length,6) = dvdrho(1:vec_length,6) + 0
      endif

      ! optional: evaluate the energy functional itself
      ! (may be used for debugging in grid module)
      if(present(fxc)) then
         factor = - (six/pi)**third
         tmp = factor * ( rho(1:vec_length,:) + deps )**third
         tmp = tmp * rho(1:vec_length,:)
         tmp(:,1) = (three/four) * ( tmp(:,1) + tmp(:,2) )
         fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)
      end if

    endif

  end subroutine xalpha_resp_calc


    function pthird(y)
      real(kind=r8_kind)   :: pthird
      real(kind=r8_kind),intent(in) :: y
      pthird=sign(1.0_r8_kind,y)*exp(log(abs(y))*third)
    end function pthird

    function f(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: f
      f = ffact * ((ONE+Y)**F4THRD    + (ONE-Y)**F4THRD - TWO)
    end function f

    function df(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) ::df
      DF  = DFACT *       ((ONE+Y)**THIRD     - (ONE-Y)**THIRD)
    end function df

    function ddf(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: ddf
      DDF = THIRD * DFACT*((ONE+Y)**(-F2THIRD) + (ONE-Y)**(-F2THIRD))
    end function ddf

    function xxp(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: xxp
      XXP = y*y + BP*y + CP
    end function xxp

    function xxf(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: xxf
      xxf = y*y + bf*y +  cf
    end function xxf

    function dxxf(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: dxxf
      dxxf = y+y + bf
    end function dxxf

    function xxa(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: xxa
      XXA = y*y + BA*y + CA
    end function xxa

    function dxxa(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: dxxa
      DXXA = y+y + BA
    end function dxxa

    function fxp(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: fxp
      FXP = ATAN(QP/(TWO*y+BP))
    end function fxp

    function fxf(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: fxf
      FXF = ATAN(QF/(TWO*y+BF))
    end function fxf

    function dfxf(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: dfxf
      DFXF = -two*QF/(QF**2+(TWO*y+BF)**2)
    end function dfxf

    function fxa(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: fxa
      FXA = ATAN(QA/(TWO*y+BA))
    end function fxa

    function arsinh(y)
      real(kind=r8_kind),dimension(vec_length),intent(in) :: y
      real(kind=r8_kind),dimension(vec_length) :: arsinh
      ARSINH = LOG(y + SQRT(y**2 + ONE))
    end function arsinh

    function s_f(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_f
      s_f = ffact * ((ONE+Y)**F4THRD    + (ONE-Y)**F4THRD - TWO)
    end function s_f

    function s_df(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) ::s_df
      s_DF  = DFACT *       ((ONE+Y)**THIRD     - (ONE-Y)**THIRD)
    end function s_df

    function s_ddf(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_ddf
      s_DDF = THIRD * DFACT*((ONE+Y)**(-F2THIRD) + (ONE-Y)**(-F2THIRD))
    end function s_ddf

    function s_xxp(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_xxp
      s_XXP = y*y + BP*y + CP
    end function s_xxp

    function s_xxf(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_xxf
      s_xxf = y*y + bf*y +  cf
    end function s_xxf

    function s_xxa(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_xxa
      s_XXA = y*y + BA*y + CA
    end function s_xxa

    function s_fxp(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_fxp
      s_FXP = ATAN(QP/(TWO*y+BP))
    end function s_fxp

    function s_fxf(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_fxf
      s_FXF = ATAN(QF/(TWO*y+BF))
    end function s_fxf

    function s_fxa(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_fxa
      s_FXA = ATAN(QA/(TWO*y+BA))
    end function s_fxa

    function s_arsinh(y)
      real(kind=r8_kind),intent(in) :: y
      real(kind=r8_kind) :: s_arsinh
      s_ARSINH = LOG(y + SQRT(y**2 + ONE))
    end function s_arsinh

end module vwnc


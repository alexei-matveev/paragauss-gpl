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
module exchange
  !-------------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  ! Copyright (c) 2010-2012 Thomas Soini
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
# define AUTODIFF ad2x1
  use type_module, only: &
       ik=>i4_kind, &
       rk=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  integer(ik), parameter, public :: &
       X_PW91      =  1, &
       X_BECKE86   =  2, &
       X_PBE       =  3, &
       X_REVPBE    =  4, &
       X_PBESOL    =  5, &
       X_PBEN      =  6, &
       X_BECKE88   =  7, &
       X_ECMV92    =  8, &
       X_XALPHA    =  9, &
       X_VMT       = 10, &
       X_VT84      = 11, &
       X_VMTsol    = 12, &
       X_VT84sol   = 13, &
       X_LB94      = 14, &
       X_GGA_MASK  = 15    ! for IAND(X,X_GGA_MASK)

  integer(ik), parameter, public :: &
       X_LONG      = 16, &
       X_TRANS     = 32, &
       X_LONGTRANS = X_LONG + X_TRANS, & !=48
       X_RELAT     = X_LONG + X_TRANS    !=48

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: exchange_gga
  public :: exchange_lda
  public :: correlation_lda

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !! SB variant (no more in use)
  !! NOTE for second derivatives packing:
  !!              1    2    3    4    5    6
  !! dfdg        AA   BB   AB
  !! dfdndn      AA   AB   BB
  !! dfdndg     AAA  BBB  AAB  ABB  BAB  BAA
  !! dfdgdg    AAAA BBBB ABAB AABB AAAB BBAB

  !! VN variant (will be implemented and used everywhere)
  !! NOTE for second derivatives packing:
  !!              1    2    3    4    5    6
  !! dfdg        AA   BB   AB
  !! dfdndn      AA   BB   AB
  !! dfdndg     AAA  BBB  BAA  ABB  AAB  BAB
  !! dfdgdg    AAAA BBBB AABB AAAB BBAB ABAB


  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !********************************************************************
  subroutine correlation_lda(vl,ispin,n,F,dFdn)
    use AUTODIFF
    use constants, only: zero, one
    use vwnc, only: vwn_calc
    implicit none
    real(rk)   , intent(in)    :: n(:,:)    ! (vl,ispin)
    real(rk)   , intent(inout) :: F(:)      ! (vl)
    real(rk)   , intent(inout) :: dFdn(:,:) ! (vl,ispin)
    ! derivatives of rho with respect to rho
    integer(ik), intent(in)    :: ispin,vl
    !** End of interface *****************************************

    integer(IK) :: s
    type(ad)    :: rho(vl), fi0(vl)

    F(:vl)      = zero
    dFdn(:vl,:) = zero
    call vwn_calc(n, dFdn, ispin, F, vl)

    ! independent variable, rank-1, var(...) is elemental:
    rho = var(sum(n(:vl, 1:ispin), 2))

    ! want the correction only:
    fi0 = phic0(rbeta(rho)) - one

    !
    ! F = F * fi0
    !
    do s = 1, ispin
        ! Elemental val(...) returns the value, fst(...) is the derivative:
        dFdn(:vl, s) = dFdn(:vl, s) * val(fi0) + F(:vl) * fst(fi0)
    enddo
    F(:vl) = F(:vl) * val(fi0)
  end subroutine correlation_lda

  elemental function phic0(b) result(phi)
    use AUTODIFF
    use constants
    implicit none
    type(ad), intent(in) :: b
    type(ad)             :: phi
    ! *** end of interface ***

    type(ad) :: nomin, denom, b2, b3, b4, b7, logb

    ! parameters and formula from
    ! Schmid, Engel, Dreizler, Blaha, Schwarz
    ! "Full potential linearized-augmented-plane-wave calculaltions..."
    real(RK), parameter ::   &
         a1_ = -2.44968_rk,   &
         a2_ =  1.91853_rk,   &
         a3_ =  0.0718854_rk, &
         b1_ = -1.59583_rk,   &
         b2_ =  1.29176_rk,   &
         b3_ =  0.364044_rk
    real(RK), parameter :: &
         CLIGHT = 137.03604_rk

    real(RK) :: AA, BB
    ! VWN = - [ 0.062 * ln(beta*C) + 0.052], when rho -> Inf
    AA = 0.0621841_rk
    BB = AA*log(CLIGHT) + 0.05275877384_rk

    if ( b > zero ) then ! because of LOG(beta)
        b2 = b * b
        b3 = b2 * b
        b4 = b2 * b2
        b7 = b4 * b3
        logb = log(b)
        nomin = one + a1_ * b3 * logb               &
             &      + a2_ * b4                      &
             &      + a3_ * (one + b2)**2 * b4
        denom = one + b1_ * b3 * logb &
             &      + b2_ * b4 &
             &      + b3_ * (AA * logb + BB) * b7
        phi = nomin / denom
    else
       phi = fix(one)
    endif
  end function phic0

  !********************************************************************
  subroutine exchange_gga(itype, vl, ispin, n, gam &
                         , F, dFdn, dFdg, dFdndn, dFdndg, dFdgdg)
    !
    ! This code  uses the identity  between UKS functional  E(...) and
    ! RKS functional e(...) for exchange:
    !
    !                           1
    !    E(r , r , g  , g  ) = --- [e(2r , 4g ) + e(2r , 4g )]
    !       a   b   aa   bb     2       a    aa       b    bb
    !
    ! The function x_gga() implements the RKS functional e(...).
    !
    use constants, only: zero, half, two, four, eight
    implicit none
    integer(ik), intent(in)    :: itype
    real(rk)   , intent(in)    :: n(:,:)    ! (vl,ispin)
    real(rk)   , intent(in)    :: gam(:,:)  ! (vl,2*ispin-1)
    ! gam(:,1)=g_aa, gam(:,2)=g_bb, gam(:,3)=g_ab
    real(rk)   , intent(inout) :: F(:)      ! (vl)
    real(rk)   , intent(inout) :: dFdn(:,:) ! (vl,ispin)
    real(rk)   , intent(inout) :: dFdg(:,:) ! (vl,2*ispin-1)
    ! second derivatives
    real(rk)   , intent(inout), optional :: dFdndn(:,:)
    real(rk)   , intent(inout), optional :: dFdndg(:,:)
    real(rk)   , intent(inout), optional :: dFdgdg(:,:)
    ! derivatives of rho with respect to rho and gam
    integer(ik), intent(in)      :: ispin,vl
    !** End of interface *****************************************

    real(rk) :: GGA(vl, 6) ! [F,dn,dg,dndn,dndg,dgdg]

    select case(ispin)
    case (1)
       GGA = zero
       if ((present(dFdndn)).and.(present(dFdndg)).and.(present(dFdgdg))) then
          call x_gga(itype, n(:vl, 1), gam(:vl, 1) &
                    , GGA(:, 1), GGA(:, 2), GGA(:, 3), GGA(:, 4), GGA(:, 5), GGA(:, 6))
          dFdndn(:vl, 1) = dFdndn(:vl, 1) + GGA(:, 4)
          dFdndg(:vl, 1) = dFdndg(:vl, 1) + GGA(:, 5)
          dFdgdg(:vl, 1) = dFdgdg(:vl, 1) + GGA(:, 6)
       else
          call x_gga(itype, n(:vl, 1), gam(:vl, 1) &
                    , GGA(:, 1), GGA(:, 2), GGA(:, 3))
       end if
       F   (:vl)    = F   (:vl)    + GGA(:, 1)
       dFdn(:vl, 1) = dFdn(:vl, 1) + GGA(:, 2)
       dFdg(:vl, 1) = dFdg(:vl, 1) + GGA(:, 3)
    case (2)
       if ((present(dFdndn)).and.(present(dFdndg)).and.(present(dFdgdg))) then
          GGA = zero
          call x_gga(itype, TWO * n(:vl, 1), FOUR * gam(:vl, 1) &
                    , GGA(:, 1), GGA(:, 2), GGA(:, 3), GGA(:, 4), GGA(:, 5), GGA(:, 6))
          F     (:vl)    = F     (:vl)    + HALF  * GGA(:, 1)
          dFdn  (:vl, 1) = dFdn  (:vl, 1) +         GGA(:, 2)
          dFdg  (:vl, 1) = dFdg  (:vl, 1) + TWO   * GGA(:, 3)
          dFdndn(:vl, 1) = dFdndn(:vl, 1) + TWO   * GGA(:, 4)
          dFdndg(:vl, 1) = dFdndg(:vl, 1) + FOUR  * GGA(:, 5)
          dFdgdg(:vl, 1) = dFdgdg(:vl, 1) + EIGHT * GGA(:, 6)

          GGA = zero
          call x_gga(itype, TWO * n(:vl, 2), FOUR * gam(:vl, 2) &
                    , GGA(:, 1), GGA(:, 2), GGA(:, 3), GGA(:, 4), GGA(:, 5), GGA(:, 6))
          F     (:vl)    = F     (:vl)    + HALF  * GGA(:, 1)
          dFdn  (:vl, 2) = dFdn  (:vl, 2) +         GGA(:, 2)
          dFdg  (:vl, 2) = dFdg  (:vl, 2) + TWO   * GGA(:, 3)
          dFdndn(:vl, 2) = dFdndn(:vl, 2) + TWO   * GGA(:, 4)
          dFdndg(:vl, 2) = dFdndg(:vl, 2) + FOUR  * GGA(:, 5)
          dFdgdg(:vl, 2) = dFdgdg(:vl, 2) + EIGHT * GGA(:, 6)
       else
          GGA = zero
          call x_gga(itype, TWO * n(:vl, 1), FOUR * gam(:vl, 1) &
                    , GGA(:, 1), GGA(:, 2), GGA(:, 3))
          F     (:vl)    = F     (:vl)    + HALF  * GGA(:, 1)
          dFdn  (:vl, 1) = dFdn  (:vl, 1) +         GGA(:, 2)
          dFdg  (:vl, 1) = dFdg  (:vl, 1) + TWO   * GGA(:, 3)

          GGA = zero
          call x_gga(itype, TWO * n(:vl, 2), FOUR * gam(:vl, 2) &
                    , GGA(:, 1), GGA(:, 2), GGA(:, 3))
          F     (:vl)    = F     (:vl)    + HALF  * GGA(:, 1)
          dFdn  (:vl, 2) = dFdn  (:vl, 2) +         GGA(:, 2)
          dFdg  (:vl, 2) = dFdg  (:vl, 2) + TWO   * GGA(:, 3)
       end if
    case default
       ABORT('no such case')
    end select
  end subroutine exchange_gga
  !********************************************************************

#if 0
  subroutine exc_lb94(vl,ispin,n,gam,dFdn)
    use constants, only: zero,two,four,half
    implicit none
    real(rk)   , intent(in)    :: n(:,:)    ! (vl,ispin)
    real(rk)   , intent(in)    :: gam(:,:)  ! (vl,2*ispin-1)
    ! gam(:,1)=g_aa, gam(:,2)=g_bb, gam(:,3)=g_ab
    real(rk)   , intent(inout) :: dFdn(:,:) ! (vl,ispin)
    ! derivatives of rho with respect to rho and gam
    integer(ik), intent(in)    :: ispin,vl
    real(rk)   , dimension(vl) :: LB94 ![dn]
    !** End of interface ****************************************

    select case(ispin)
    case (1)
       LB94= zero
       call xc_lb94(vl,n(:vl,1),gam(:vl,1),LB94(:vl))
       dFdn(:vl,1) = dFdn(:vl,1) + LB94(:vl)
    case (2)
       LB94= zero
       call xc_lb94(vl,TWO*n(:vl,1),FOUR*gam(:vl,1),LB94(:vl))
       dFdn(:vl,1) = dFdn(:vl,1) + LB94(:vl)
       LB94 = zero
       call xc_lb94(vl,TWO*n(:vl,2),FOUR*gam(:vl,2),LB94(:vl))
       dFdn(:vl,2) = dFdn(:vl,2) + LB94(:vl)
    end select

  end subroutine exc_lb94
#endif

  elemental subroutine x_gga(itype, n, gam, F, dFdn, dFdg, dFdndn, dFdndg, dFdgdg)
    use AUTODIFF
    use constants
    implicit none
    integer(ik), intent(in)    :: itype
    real(rk)   , intent(in)    :: n
    real(rk)   , intent(in)    :: gam
    real(rk)   , intent(inout) :: F
    real(rk)   , intent(inout) :: dFdn
    real(rk)   , intent(inout) :: dFdg
    ! second derivatives
    real(rk)   , intent(inout), optional :: dFdndn
    real(rk)   , intent(inout), optional :: dFdndg
    real(rk)   , intent(inout), optional :: dFdgdg
    ! *** end of interface ***

    integer(ik), parameter :: &
         D0  = 0, & ! value
         DR  = 1, & ! deriv wrt RHO
         DRR = 2, & ! deriv wrt RHO RHO    (second derivatives)
         DG  = 3, & ! deriv wrt GAMMA
         DRG = 4    ! deriv wrt RHO GAMMA  (second derivatives)

    real(rk) :: xalpha, s2_fact
    real(rk) :: s2(D0:DRG) ! D0, DR, DRR, DG, DRG

    integer(ik) :: igga, irel
    type(ad) :: rho, ex, G, beta, fi0, fi2

    ! independent variable:
    rho = var(n)

    ! LDAX:
    ! -(3/4) * [(3/pi)^(1/3)] * rho^(4/3)
    xalpha  = - (three/four) * (three/pi)**(one/three)

    !
    ! The second derivative diverge as 0**(-1/3) = Inf:
    !
    if (rho > 0.0d0) then
       ex = xalpha * rho**(four/three)
    else
       !
       ! 0**(4/3) = 0**(1/3) = 0, higher derivatives are forced to 0:
       !
       ex = fix(0.0d0)
    endif

    s2_fact   = one / (four * (three*pi**2)**(two/three))
    if (n >= 1.0E-20_rk) then
       s2(DG)  = one / n**(8.D0/3.D0) * s2_fact
       s2(D0)  = gam * s2(DG) ! linear in gamma
       s2(DR)  = - ( 8.D0/3.D0) * s2(D0) / n
       s2(DRR) = - (11.D0/3.D0) * s2(DR) / n
       s2(DRG) = - ( 8.D0/3.D0) * s2(DG) / n
    else
       s2(:)  = zero
    endif

    !
    ! Compute GGA enhancement factor
    !

    ! gga bits:
    igga = IAND(itype, X_GGA_MASK)

    ! value, first and second derivatives wrt s2:
    G = g_factor(igga, var(s2(D0)))

    ! relativistic bits:
    irel = IAND(itype, X_RELAT)

    if( irel == 0 ) then
       ! no relativity

       F    = F    + val(ex) * val(G)    ! 1 of (1+G) is omitted

       dFdn = dFdn + fst(ex) * val(G) &  ! 1 of (1+G) is omitted
                   + val(ex) * fst(G) * s2(DR)

       dFdg = dFdg + val(ex) * fst(G) * s2(DG)

       ! second derivatives
       if ( present(dFdndn) ) then
!         ASSERT(present(dFdndg))
!         ASSERT(present(dFdgdg))

          dFdndn = dFdndn + snd(ex) * val(G) &
                          + fst(ex) * fst(G) * s2(DR) * two &
                          + val(ex) * snd(G) * s2(DR)**2    &
                          + val(ex) * fst(G) * s2(DRR)

          dFdndg = dFdndg + fst(ex) * fst(G) * s2(DG) &
                          + val(ex) * snd(G) * s2(DG) * s2(DR) &
                          + val(ex) * fst(G) * s2(DRG)

          dFdgdg = dFdgdg + val(ex) * snd(G) * s2(DG)**2
       end if

    else

       ! rho is the independent variable:
       beta = rbeta(rho)

       ! value, first and second derivatives wrt rho:
       fi0 = phi0(irel, beta)  ! pass X_LONG || X_TRANS bits
       fi2 = phi2(itype, beta) ! pass X_LONG || X_TRANS || X_GGA_MASK bits

       F    = F    + val(ex) * ( val(fi0) - one + val(fi2) * val(G) ) ! -1 is NO XALPHA
       dFdn = dFdn + fst(ex) * ( val(fi0) - one + val(fi2) * val(G) ) &
                   + val(ex) * ( fst(fi0)       + fst(fi2) * val(G)   &
                                                + val(fi2) * fst(G) * s2(DR) )
       dFdg = dFdg + val(ex) *                    val(fi2) * fst(G) * s2(DG)
    endif
  end subroutine x_gga

  elemental function g_factor(igga, s2) result(G)
    use AUTODIFF, only: ad, fix, nan
    implicit none
    integer(ik), intent(in) :: igga
    type(ad), intent(in)    :: s2
    type(ad)                :: G
    ! *** end of interface ***

    select case( igga )
    case (X_PW91)
        call g_pw91(s2, G)
    case (X_BECKE86)
        call g_becke86(s2, G)
    case (X_BECKE88)
        call g_becke88(s2, G)
    case (X_PBE)
        call g_pbe(s2, G)
    case (X_REVPBE)
        call g_revpbe(s2, G)
    case (X_PBESOL)
        call g_pbesol(s2, G)
    case (X_PBEN)
        call g_pben(s2, G)
    case (X_ECMV92)
        call g_ecmv92(s2, G)
    case (X_VMT)
        call g_vmt(s2, G)
    case (X_VT84)
        call g_vt84(s2, G)
    case (X_VMTsol)
        call g_vmtsol(s2, G)
    case (X_VT84sol)
        call g_vt84sol(s2, G)
    case (X_LB94)
        call g_lb94(s2, G)
    case (X_XALPHA)
        G = fix(0.0_rk) ! FIXME: lazy
    case default
        ! ABORT('g_factor: no such igga')
        G = nan(s2)
    end select
  end function g_factor

  subroutine xc_lb94(n, gam, dFdn)
    use AUTODIFF
    use constants
    implicit none
    real(rk)   , intent(in)    :: n
    real(rk)   , intent(in)    :: gam
    real(rk)   , intent(inout) :: dFdn
    ! *** end of interface ***

    real(rk) :: ex
    real(rk) :: s2
    real(rk) :: G(0:2)

    ex  = n**(one/three)
    if( abs(n) .GE. 1.0E-20_rk )then
       s2 = gam / n**(8.D0/3.D0)
    else
       s2 = zero
    endif

    G = taylor(g_factor(X_LB94, var(s2)))

    ! FIXME: this doesnt look right:
    dFdn   = dFdn + ex * G(0)
  end subroutine xc_lb94

  !********************************************************************
  subroutine exchange_lda(irel, vl, ispin, n, F, dFdn, dFdndn)
    use AUTODIFF
    use constants, only: zero, half, two, four
    implicit none
    integer(ik), intent(in)    :: irel
    real(rk)   , intent(in)    :: n(:,:)    ! (vl,ispin)
    real(rk)   , intent(inout) :: F(:)      ! (vl)
    real(rk)   , intent(inout) :: dFdn(:,:) ! (vl,ispin)
    ! second derivatives
    real(rk)   , intent(inout), optional :: dFdndn(:,:)
    ! derivatives of rho with respect to rho and gam
    integer(ik), intent(in)    :: ispin,vl
    !** End of interface *****************************************

    type(ad) :: LDA(vl)

    select case(ispin)
    case (1)
        ! total density:
        LDA = x_lda(irel, var(n(:vl, 1)))

        F(:vl) = F(:vl) + val(LDA)

        dFdn(:vl, 1) = dFdn(:vl, 1) + fst(LDA)

        if (present(dFdndn)) then
            dFdndn(:vl, 1) = dFdndn(:vl, 1) + snd(LDA)
        endif

    case (2)
        ! alpha density:
        LDA = x_lda(irel, var(2 * n(:vl, 1)))

        F(:vl) = F(:vl) + val(LDA) / 2

        dFdn(:vl, 1) = dFdn(:vl, 1) + fst(LDA)

        if (present(dFdndn)) then
            dFdndn(:vl, 1) = dFdndn(:vl, 1) + 2 * snd(LDA)
        endif

        ! beta density:
        LDA = x_lda(irel, var(2 * n(:vl, 2)))

        F(:vl) = F(:vl) + val(LDA) / 2

        dFdn(:vl, 2) = dFdn(:vl, 2) + fst(LDA)

        if (present(dFdndn)) then
            dFdndn(:vl, 2) = dFdndn(:vl, 2) + 2 * snd(LDA)
        endif
    case default
       ABORT('no such case')
    end select
  end subroutine exchange_lda

  elemental function x_lda(itype, rho) result(ex)
    use AUTODIFF
    use constants
    implicit none
    integer(ik), intent(in) :: itype
    type(ad), intent(in)    :: rho
    type(ad)                :: ex
    ! *** end of interface ***

    real(rk)    :: xalpha
    integer(ik) :: irel

    ! LDAX:
    ! -(3/4) * [(3/pi)^(1/3)] * rho^(4/3)
    xalpha  = - (three/four) * (three/pi)**(one/three)

    !
    ! The second derivative diverge as 0**(-1/3) = Inf:
    !
    if ( rho > 0.0d0 ) then
       ex = xalpha * rho**(four/three)
    else
       !
       ! 0**(4/3) = 0**(1/3) = 0, higher derivatives are forced to 0:
       !
       ex = fix(0.0d0)
    endif

    ! relativity
    irel = IAND(itype, X_RELAT)

    if( irel /= 0 ) then
       ex = ex * phi0(irel, rbeta(rho))
    endif
  end function x_lda

  elemental subroutine g_becke86(s2, G)
    use AUTODIFF
    implicit none
    type(ad)   , intent(in)  :: s2
    type(ad)   , intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter ::  &
         kappa = 0.967_rk, &
         mu    = 0.235_rk

    call g_becke86_form(s2, G, kappa, mu)
  end subroutine g_becke86

  elemental subroutine g_pbe(s2, G)
    use AUTODIFF
    use constants, only: pi, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter :: &
         kappa = 0.804_rk, & ! revised kappa=1.245
         beta  = 0.066725_rk, &
         mu    = beta * ( pi**2 / three)
    !
    !          = 0.21951645122089583
    !
    ! compare with
    !
    !       mu = 0.2195149727645171d0
    !
    ! e.g. in http://dft.uci.edu/pubs/PBEsol.html
    !
    call g_becke86_form(s2, G, kappa, mu)
  end subroutine g_pbe

  elemental subroutine g_revpbe(s2, G)
    use AUTODIFF
    use constants, only: pi, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter :: &
         kappa = 1.245_rk, & ! original kappa=0.804
         beta  = 0.066725_rk, & ! same as in PBE
         mu    = beta * ( pi**2 / three)

    call g_becke86_form(s2, G, kappa, mu)
  end subroutine g_revpbe

  elemental subroutine g_pbesol(s2, G)
    use AUTODIFF
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter :: kappa = 0.804_rk ! same as in PBE
    real(rk), parameter :: mu = 10.0_rk / 81.0_rk ! PBEsol
    !
    !                         = 0.12345679012345678
    !
    call g_becke86_form(s2, G, kappa, mu)
  end subroutine g_pbesol

  elemental subroutine g_becke88(s2, G)
    use AUTODIFF
    use constants, only: pi, three, six, two, one, four, nine
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk)            :: &
         B != two * (six*pi**2)**(one/three)

    ! two ways of parameter definitions differ in accuracy
    ! I (am) would persanally prefer BECKE88_2743 way,
    ! but PG was historically using BECKE88_0042
#ifdef BECKE88_2743
    ! BECKE88_2743, one way to define parameters:
    real(rk), parameter ::  &
         D = 0.2743_rk, &
         A = nine * D / (four * pi)
#else
    ! BECKE88_0042, another way:
    real(rk), parameter :: &
         beta = 0.0042_rk
    real(rk)            :: &
         A,D
    A = six * beta * two**(one/three) &
         &  * (three*pi**2)**(one/three) * two
    D = A   * (four * pi) / nine
#endif

    B = two * (six*pi**2)**(one/three)

    call g_becke88_form(s2, G, A, B, D)
  end subroutine g_becke88

  elemental subroutine g_lb94(s2, G)
    use AUTODIFF
    use constants, only: three, one
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter      :: &
         A = three * 0.05_rk, &
         B = one, &
         D = - 0.05_rk

    call g_becke88_form(s2, G, A, B, D)
  end subroutine g_lb94

  elemental subroutine g_vmt(s2, G)
    use AUTODIFF
    use constants, only: pi, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***
    ! Parameters beta / mu done identical to PBE for better comparison
    ! Parameter alpha: see J. Chem. Phys. 130, 244103
    real(rk), parameter :: beta  = 0.066725_RK                        &
                         , mu    = beta * ( pi * pi / three)          &
                         , alpha = 0.002762_RK

    call g_vmt_form(s2, G, alpha, mu)
  end subroutine g_vmt

  elemental subroutine g_vt84(s2, G)
    use AUTODIFF
    use constants, only: pi, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    type(ad) :: term1, term2

    ! Parameters beta / mu done identical to PBE for better comparison
    ! Parameter alpha: see Trickey presentation, september 2010
    real(rk), parameter :: beta  = 0.066725_RK                        &
                         , mu    = beta * ( pi * pi / three)          &
                         , alpha = 0.000069_RK

    call g_vmt_form(s2, term1, alpha, mu)

    call g_vt84_2nd_term(s2, term2, alpha)

    G = term1 + term2
  end subroutine g_vt84

  elemental subroutine g_vmtsol(s2, G)
    use AUTODIFF
    use constants, only: pi, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***
    ! Parameter done identical to PBEsol for better comparison
    ! Parameter alpha: see J. Chem. Phys. 130, 244103
    real(rk), parameter :: mu    = 10.0_RK / 81.0_RK                  &
                         , alpha = 0.0015532_RK

    call g_vmt_form(s2, G, alpha, mu)
  end subroutine g_vmtsol

  elemental subroutine g_vt84sol(s2, G)
    use AUTODIFF
    use constants, only: pi, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    type(ad) :: term1, term2

    ! Parameter done identical to PBEsol for better comparison
    ! Parameter alpha: see Trickey presentation, september 2010
    real(rk), parameter :: mu    = 10.0_RK / 81.0_RK                  &
                         , alpha = 0.000023_RK

    call g_vmt_form(s2, term1, alpha, mu)

    call g_vt84_2nd_term(s2, term2, alpha)

    G = term1 + term2
  end subroutine g_vt84sol

  elemental subroutine g_vmt_form(s2, G, alpha, mu)
      !
      !                               2
      !                   2   -alpha s
      !          2    mu s   e
      !       g(s ) = ------------------
      !                       2
      !                   mu s  + 1
      !
      use AUTODIFF
      use constants
      implicit none
      type(ad), intent(in)    :: s2
      type(ad), intent(out)   :: G
      real(rk), intent(in)    :: alpha, mu
      ! *** end of interface ***

      G = mu * s2 * exp( - alpha * s2) / (mu * s2 + one)
  end subroutine g_vmt_form

  elemental subroutine g_vt84_2nd_term(s2, G, alpha)
      !
      !     (1 - exp(-alpha s^(m/2))) (s^(-n/2) -1)
      !
      ! where m = 8, n = 4
      !
      use AUTODIFF
      use constants
      implicit none
      type(ad), intent(in)     :: s2
      type(ad), intent(out)    :: G
      real(rk), intent(in)     :: alpha
      ! *** end of interface ***

      if ( s2 > 1.0e-4_rk ) then
        G = (one - exp(- alpha * s2 * s2)) * (one / s2 - one)
      else
        !
        !                      2    1                  -4
        ! g(x) = (1 - exp(- a x )) (- - 1), with a < 10
        !                           x
        !
        ! taylor(g(x)/a, x, 0) =
        !
        !               3      4    2  5
        !        2   a x    a x    a  x
        ! = x - x  - ---- + ---- + ----- + . . .
        !             2      2       6
        !
        !                            2
        !     -4                    a    -20
        ! = 10   +    . . .      + --- 10
        !                           6
        !
        G = alpha * s2 &
          * (one + s2 * (-one + s2 * (-alpha/2 + s2 * (alpha/2 + s2 * (alpha**2/6)))))
      end if
  end subroutine g_vt84_2nd_term

  elemental subroutine g_becke86_form(s2, G, kappa, mu)
    use AUTODIFF
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    real(rk), intent(in)  :: kappa, mu
    ! *** end of interface ***

    G = kappa * mu * s2 / (kappa + mu * s2)
  end subroutine g_becke86_form

  elemental subroutine g_pben(s2, G)
    use AUTODIFF
    use constants, only: pi, one, three
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter :: &
         kappa = 0.804_rk, &
         beta  = 0.066725_rk, &
         mu    = beta * ( pi**2 / three)

    type(ad) :: x

    x = mu * s2 / kappa

    if ( x > 1.0e-4_rk ) then
        G = kappa * ( one - exp(-x) )
    else
        !
        ! f(x) = 1 - exp(-x) =
        !
        !              2    3    4    5
        !             x    x    x    x
        !      =  x - -- + -- - -- + --- - ...
        !             2    6    24   120
        !
        !      ~ 1e-4 +    ...     + 1e-22
        !
        G = kappa * x &
          * (one + x * (-one/2 + x * (one/6 + x * (-one/24 + x * (one/120)))))
    endif
  end subroutine g_pben

  elemental subroutine g_ecmv92(s2, G)
    use AUTODIFF
    use constants, only: one, zero
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter :: & ! ECMV92 set for GGA
         A(3) = (/ zero,  0.3402_rk, 5.9955_rk /), &
         B(3) = (/  one, 27.5026_rk, 5.7728_rk /)

    call pade22(s2, G, A(1), A(2), A(3), B(1), B(2), B(3))
  end subroutine g_ecmv92

  elemental subroutine g_pw91(s2, G)
    !
    !      2         2                2
    !     s  (D - C s  - E exp(- 100 s ))
    ! G = -------------------------------- =
    !                                  4
    !        1 + A s arcsinh(B s) + C s
    !
    !                2                                4
    !   = ( D - E ) s  + ((B A + 100) E - B A D - C) s  + . . .
    !
    use AUTODIFF
    use constants, only: one
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    ! *** end of interface ***

    real(rk), parameter ::  &
         A = 0.19645_rk, &
         B = 7.7956_rk, &
         C = 0.004_rk, &
         D = 0.2743_rk, &
         E = 0.1508_rk, &
         HNDRT = 100.0_rk

    real(rk), parameter :: T1 = D - E
    real(rk), parameter :: T2 = (B * A + HNDRT) * E - B * A * D - C

    type(ad) :: s
    type(ad) :: nomin
    type(ad) :: denom

    if ( s2 > 1.0e-4_rk ) then
        nomin = D - C * s2 - E * exp( -HNDRT * s2 )

        s = sqrt(s2)

        denom = one + A * s * arcsinh(B*s) + C * s2**2
        G = s2 * nomin / denom
    else
        G = s2 * (T1 + s2 * T2) ! + O(1.0e-24)
    endif
  end subroutine g_pw91

  elemental subroutine g_becke88_form(s2, G, A, B, D)
    !
    !                    2
    !                 D s
    !  G(s) = --------------------
    !         1 + A s arcsinh(B s)
    !
    !            2          4
    !       = D s  - B A D s  + . . .
    !
    use AUTODIFF
    use constants, only: one
    implicit none
    type(ad), intent(in)  :: s2
    type(ad), intent(out) :: G
    real(rk), intent(in)  :: A, B, D
    ! *** end of interface ***

    type(ad) :: s

    if ( s2 > 1.0e-4_rk ) then
        s = sqrt(s2)

        G = s2 * D / (one + A * s * arcsinh(B*s))
    else
        G = s2 * (D + s2 * (- B * A * D))
    endif
  end subroutine g_becke88_form

  elemental function rbeta(n) result(beta)
    ! safe way to compute relativistic beta
    use AUTODIFF
    use constants, C=>SPEED_OF_LIGHT
    implicit none
    type(ad), intent(in) :: n
    type(ad)             :: beta
    ! *** end of interface ***

    if ( n .ge. 1.0E-20_rk ) then
       beta = (three * pi**2 * n)**(one/three) / C
    else
       beta = fix(zero)
    endif
  end function rbeta

  elemental function phi0(irel, b) result(phi)
    use AUTODIFF
    use constants
    implicit none
    integer(ik), intent(in) :: irel ! used bits: X_LONG || X_TRANS
    type(ad), intent(in)  :: b
    type(ad)              :: phi
    ! *** end of interface ***

    type(ad) :: term1
    type(ad) :: term2
    type(ad) :: term3
    type(ad) :: term4
    type(ad) :: b2, b4, eta, ash

    b2 = b * b
    b4 = b2 * b2

    eta = sqrt(one + b2)
    ash = arcsinh(b)

    if( IAND(irel, X_LONGTRANS) .ne. X_LONGTRANS )then
    ! either X_LONG or X_TRANS
    ! only in this case compute term1, term2 and term3
    term1 =  (one/three) / b2

    term2 = (two/three) * eta * ash / b

    term3 = (two/three) * eta**4 * log(eta) / b4
    endif

    ! first not squared:
!   term4 = eta / b - ash / b2
    call t4(b,term4)

    ! now squared:
    term4 = term4 * term4

    select case( IAND(irel, X_LONGTRANS) )
    case (X_LONG)
       phi =   term1 + term2 - term3 - half * term4
       phi = phi + five/six
    case (X_TRANS)
       phi = - term1 - term2 + term3 -        term4
       phi = phi + one/six
    case (X_LONGTRANS)
       phi = - (three/two) * term4
       phi = phi + one
    case default
       ! ABORT('phi0: no such irel')
       phi = nan(b)
    end select

    ! FIXME: do it gracefully:
    if ( val(b) == 0.0_rk ) then
       phi = fix(0.0_rk)
    endif

  contains
    pure subroutine t4(b, term4)
      implicit none
      type(ad), intent(in)  :: b
      type(ad), intent(out) :: term4
      ! *** end of interface ***

      type(ad) :: b2, b4, res(3)

      if ( abs(val(b)) > 1.0E-4_rk ) then ! FIXME: double precision
         term4 = (eta - ash / b) / b
      else
         b2 = b * b
         b4 = b2 * b2

         res(1) =   ( 2.0_rk /   3.0_rk) * b    !  6.666666666666667E-005
         res(2) = - ( 1.0_rk /   5.0_rk) * b*b2 ! -2.000000000000000E-013
         res(3) =   ( 3.0_rk /  28.0_rk) * b*b4 !  1.071428571428572E-021
         !es(4) = - ( 5.0_rk /  72.0_rk) * b**7 ! -6.944444444444446E-030
         !es(5) = + (35.0_rk / 704.0_rk) * b**9

         term4 = res(1) + res(2) + res(3)
      end if
    end subroutine t4
  end function phi0

  elemental function phi2(irel, b) result(phi)
    use AUTODIFF
    use constants, only: one, zero
    implicit none
    integer(ik), intent(in)  :: irel ! used bits: X_LONG || X_TRANS || X_GGA_MASK
    type(ad), intent(in)  :: b
    type(ad)              :: phi
    ! *** end of interface ***

    real(rk), parameter :: & ! RECMV92 set
         AL92(3) = (/  one, 2.21259_rk, 0.669152_rk /), &
         BL92(3) = (/  one, 1.32998_rk, 0.794803_rk /), &
         AT92(3) = (/ zero, 3.48754_rk, 0.218599_rk /), &
         BT92(3) = (/  one, 1.15417_rk, 0.015802_rk /)

    real(rk), parameter :: & ! RB88 set
         AL88(3) = (/  one, 2.20848_rk, 0.668684_rk /), &
         BL88(3) = (/  one, 1.33075_rk, 0.795105_rk /), &
         AT88(3) = (/ zero, 3.48638_rk, 0.614753_rk /), &
         BT88(3) = (/  one, 1.32260_rk, 0.101805_rk /)

    type(ad) :: b2, phiX

    ! there are eventually two contributions:
    phi = fix(zero)

    b2 = b * b

    if( IAND(irel, X_LONG) /= 0 )then

       if( X_ECMV92 == IAND(irel, X_GGA_MASK) )then
          call pade22(b2, phiX, AL92(1), AL92(2), AL92(3), BL92(1), BL92(2), BL92(3))
       else
          call pade22(b2, phiX, AL88(1), AL88(2), AL88(3), BL88(1), BL88(2), BL88(3))
       endif

       phi = phiX
    endif

    if( IAND(irel, X_TRANS) /= 0 )then

       if( X_ECMV92 == IAND(irel, X_GGA_MASK) )then
          call pade22(b2, phiX, AT92(1), AT92(2), AT92(3), BT92(1), BT92(2), BT92(3))
       else
          call pade22(b2, phiX, AT88(1), AT88(2), AT88(3), BT88(1), BT88(2), BT88(3))
       endif

       phi = phi + phiX
    endif
  end function phi2

  elemental subroutine pade22(x, y, a0, a1, a2, b0, b1, b2)
    ! 2/2 Pade function:
    ! y = [ a0 + a1*x + a2*x^2]/[ b0 + b1*x + b2*x^2 ]
    use AUTODIFF
    implicit none
    type(ad), intent(in)  :: x
    type(ad), intent(out) :: y
    real(rk), intent(in)  :: a0, a1, a2, b0, b1, b2
    ! *** end of interface ***

    type(ad) :: nomin, denom

    nomin = a0 + x * (a1 + x * a2)
    denom = b0 + x * (b1 + x * b2)
    y = nomin / denom
  end subroutine pade22

  !--------------- End of module -------------------------------------
end module exchange

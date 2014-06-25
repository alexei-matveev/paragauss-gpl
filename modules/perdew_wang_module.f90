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
module perdew_wang_module
  !---------------------------------------------------------------
  !
  !  Purpose: contains routines calculating the exchande correlation
  !           potential from Perdew and Wang 91
  !
  !
  !
  !  Module called by: xc_hamiltonian, post_scf_module
  !
  !  References: J.P. Perdew et al., Phys Rev B, 46 (1992) 6671
  ! 
  !  Author: MS
  !  Date: 1/96
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: AM
  ! Date:   9/98
  ! Description: to avoid failure in an exceptional case when
  !              zeta (pw91c_calc) is 1.0 or -1.0
  ! Date:   10/98
  ! Description: revision of PW91-correlation functional, introduced
  !              a possibility to switch off "ungeneralized" term
  !
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

# include <def.h>
  use type_module
  implicit none
  private         ! by default, all names are private
  save

  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  public pw91x_calc, pw91c_calc


  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of constants and variables ----


  real(kind=r8_kind),parameter ::  ZERO  = 0.0_r8_kind, &  
       ONE   = 1.0_r8_kind, &  
       TWO   = 2.0_r8_kind, &  
       THREE = 3.0_r8_kind, &  
       FOUR  = 4.0_r8_kind, &  
       FIVE  = 5.0_r8_kind, &  
       SIX   = 6.0_r8_kind, &  
       seven = 7.0_r8_kind, &
       eight = 8.0_r8_kind, &
       NINE  = 9.0_r8_kind, &  
       ten  = 10.0_r8_kind, &  
       HALF  = 0.5_r8_kind, &  
       DEPS  = 1.0e-50_r8_kind, & 
       ZEPS  = 1.0e-06_r8_kind, & 
       YEPS  = 1.0e-06_r8_kind, & 
       HUNDT = 100.0_r8_kind, &
       THND = 1.0e-3_r8_kind


  real(kind=r8_kind),parameter ::  THIRD  = ONE/THREE , & 
       F2THRD = TWO/THREE, &
       F4THRD = FOUR/THREE, &  
       F7o6 = seven/SIX, &
       PI = 3.141592653589793240_r8_kind

contains

  subroutine pw91c_calc(rho,gamma,dfdrho,ispin,fxc,dfdgrarho,vl,&
       fxc_lda,dfdrho_lda, revised) 
    ! Purpose: calculate the correlation energy of PW91 functional
    ! Note: no LDA contributions are included, a call to a LDA 
    !       functional must preceed a call to this functional
    real(kind=r8_kind),dimension(:,:) :: rho  ! rho(vl,ispin)
    real(kind=r8_kind),dimension(:,:) :: gamma ! gamma(vl,(ispin-1)*2+1)
    ! gamma(:,1)=grarho_down*grarho_down
    ! gamma(:,2)=grarho_up*grarho_up
    ! gamma(:,3)=grarho_down*grarho_up
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc ! f according to JGP
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho,dfdgrarho
    ! dfdrho(vl,ispin), dfdgrarho(vl,(ispin-1)*2+1)
    ! derivatives of rho with respect to rho and gamma
    integer(kind=i4_kind),intent(in) :: ispin,vl
    real(kind=r8_kind),dimension(:),intent(in) :: fxc_lda
    ! fxc of a LDA functional; must have been called previously; eg. VWN
    real(kind=r8_kind),dimension(:,:),intent(in) :: dfdrho_lda
    ! dfxc_lda/drho; see above
    logical, intent(in), optional :: revised
    !** End of interface *****************************************  



    real(kind=r8_kind),parameter :: Cc0   =  0.004235_r8_kind , &     ! PW parameters
         Cx    = -0.001667_r8_kind      , &
         AL_PW =  0.09_r8_kind          , &
         C192  =  192.0_r8_kind         , &
         Cx37    = Cx*THREE/seven , &
         a_rg  = 23.266_r8_kind         , &
         b_rg  = 7.389e-3_r8_kind       , &
         c_rg  = 8.723_r8_kind          , &
         d_rg  = 0.472_r8_kind          , &
         e_rg  = 2.568_r8_kind          , &
         c_21  = -0.458_r8_kind         , &
         c_22  = -0.037_r8_kind         , &
         c_23  =  0.10_r8_kind

    real(kind=r8_kind), dimension(vl) ::  sums, fact, rs_fact, h1_fact, b_pw, &
         h0_fact, a_pw_f, fact_ks, t_fact, e_fact, eps_lda, rhothrd, rs, rs2, x, &
         gr, t_pw, t2, ex_p0, a_den, a_pw, a_pwt2, oa_pwt2, h0, t6, h0_a, f0_t, &
         t_g, f0_g, t_n, eps_n_lda, a_n, h0_n, f0_n, ex_p1, rg_nom, rg_dnom, cc_rs, &
         h1, h1_t, h1_g, f1_g, rs_n, rg_n_, rg_dn_, h_rs, h1_dn, h1_n, f1_n, sum2, &
         sum3, sum4, sum5, zeta, gr2, g3, g_, g2, g4, z_a, z_b, t_gi, &
         vexp, a_gi, eps_lda_na, eps_lda_nb, a_na, a_nb, h0_na, h0_nb, &
         f0_na, f0_nb, dg, h0_t, h0_gi, h1_gi, h1_na, h1_nb, f1_na, f1_nb, dz2, dz2_aa, &
         dz2_bb, dz2_ab, dz2_a_s, dz2_b_s, dz_dn, dz_dn_b, dz_dn_a, dz_dn_aa, dz_dn_bb, &
         dz_dn_ab, h2, f2, f2_nd, f2_h21a, f2_h21b, f2_dz_dn_a, f2_dz_dn_b, f2_h22a, &
         h2_g_ab, h2_g_bb, h2_g_aa, f2_nb, f2_na, f2_dz2_b, f2_dz2_a, f2_h22b, &
         h21, h22, dh21, dh22

    real(kind=r8_kind)                        :: zeta_eps
    real(kind=r8_kind), dimension(vl) :: r_2nadn, r_2nbdn

    logical :: revision = .false.

    ! FIX for (PW91 == PW91-IIA) bug. Always reset revision to PW91-plain:
    revision = .false.
    if(present(revised))then
       revision=revised
    endif

         
    FACT   = THREE/(FOUR*PI)
    rs_fact= (FACT)**THIRD
    H1_fact= (FOUR*FOUR/PI)*(THREE*PI*PI)**THIRD
    b_pw   = H1_fact*Cc0
    H0_fact = b_pw**2/(two*AL_PW)
    A_PW_F = two*AL_PW/b_pw
    fact_ks= sqrt((C192/PI)**THIRD)
    t_fact = one/(two*fact_ks)
    e_fact = hundt*four/(pi*(three*pi*pi)**THIRD)


    !-----------------------------------------------------------------------
    !     CALCULATIONAL ENTRY POINT
    !-----------------------------------------------------------------------
    if (ispin==1) then  ! ** SPIN RESTRICTED CASE  **

#ifdef FPP_DEBUG
# define DSHAPE(X) count( X < zero ),'<zero<',count( X > zero ),'of',size( X )
#endif

       sums = rho(1:vl,1)
       where(sums<deps) sums = deps
       ! local energy density
       eps_lda = fxc_lda(1:vl)/sums ! - deps
       DPRINT 'pw::pw91c: eps_lda = ', DSHAPE(eps_lda)
       rhothrd = sums**third
       rs  = rs_fact/rhothrd
       rs2  = rs*rs
       x   = sqrt(rs)

       gr(:)=sqrt(gamma(1:vl,1)) ! gradient of the density
       t_pw = t_fact*gr(:)/(sqrt(rhothrd)*sums)
       t2= t_pw*t_pw
       ex_p0 = -two*al_pw*eps_lda/b_pw**2
       DPRINT 'pw::pw91c: ex_p0 = ', DSHAPE(ex_p0)
       vexp = exp(ex_p0)
       A_den = (Vexp-ONE )
       where(A_den<deps) A_den = deps
       DPRINT 'pw::pw91c: A_den = ', DSHAPE(A_den)
       A_PW = A_PW_F/A_den
       A_PWt2 = A_PW*t2
       DPRINT 'pw::pw91c: A_PWt2 = ', DSHAPE(A_PWt2)
       OA_PWt2 = one+A_PWt2
#ifdef FPP_DEBUG
       sum2 = OA_PWt2+A_PWt2**2 ! used only in polarized case
       DPRINT 'pw::pw91c: OA_PWt2+A_PWt2**2 = ', DSHAPE(sum2)
#endif
       h0 = h0_fact*log(one+A_PW_F*(t2+A_PWt2*t2)/(OA_PWt2+A_PWt2**2))

       t6 = t2**3
       h0_a = -B_PW*t6*A_PW*(two+A_PWt2)&
            /((one+(A_PW_F+A_PW)*t2*OA_PWt2)*(OA_PWt2+A_PWt2**2))

       F0_t = sums*two*B_PW*t_pw*(OA_PWt2+A_PWt2)&
            /((one+(A_PW_F+A_PW)*t2*OA_PWt2)*(OA_PWt2+A_PWt2**2))

       t_g = half*t_fact/(sqrt(rhothrd)*sums*gr+deps)
       f0_g = f0_t*t_g

       t_n = -F7o6 * t_pw/sums
       ! calculate dervative of eps_lda with respect to n
       !eps_n = -(X/SIX)*DEP/sums
       eps_n_lda = (dfdrho_lda(1:vl,1)-eps_lda)/sums
       A_n = (A_PW_F**2/B_PW)*Vexp/A_den**2*eps_n_lda
       h0_n = H0_a*A_n+ F0_t/sums*t_n
       F0_n = H0+ sums*H0_a*A_n+F0_t*t_n

       ex_p1= -e_fact*t2/rhothrd
       rg_nom = e_rg+a_rg*RS+b_rg*RS2
       rg_dnom = one+c_rg*RS+d_rg*RS2+ten*b_rg*RS2*RS
       Cc_rs = thnd*rg_nom/rg_dnom-Cx


       h1 = h1_fact*(Cc_rs-Cc0-Cx37)*t2*exp(ex_p1)
       h1_t = h1_fact*(Cc_rs-Cc0-Cx37)*exp(ex_p1) &
            *two*t_pw*(one+ex_p1)

       H1_g = H1_t*t_g 
       F1_g = sums*H1_g

       RS_n  = -third*rs/sums
       rg_n_ =a_rg+two*b_rg*RS
       rg_dn_=c_rg+d_rg*two*RS+three*ten*b_rg*RS2

       h_rs = h1_fact*t2*exp(ex_p1)*thnd*(rg_n_*rg_dnom-rg_dn_*rg_nom)/rg_dnom**2

       h1_dn =  h1*third*t2*e_fact/(sums*rhothrd)

       h1_n = h_rs*rs_n + h1_t*t_n  + h1_dn
       f1_n = h1+ sums*h1_n


       ! contribution to energy and potentials
       fxc(1:vl)=fxc(1:vl)+(h0+h1)*sums
       dfdrho(1:vl,1)=dfdrho(1:vl,1) + (F1_n+F0_n)
       dfdgrarho(1:vl,1) = dfdgrarho(1:vl,1)+(f0_g+f1_g)  

    else  ! **  spin_polarized **

       sums = rho(1:vl,1) + rho(1:vl,2)
       where(sums<deps) sums = deps
       eps_lda = fxc_lda(1:vl)/sums
       sum2 = sums * sums 
       sum3 = sum2*sums
       sum4 = sum2*sum2
       sum5 = sum4*sums

       rhothrd = sums**THIRD
       RS  = Rs_fact/rhoTHRD
       X   = SQRT(RS)
       RS2 = RS **2

       ZETA = ( rho(1:vl,1)  - rho(1:vl,2)  )/(sums)

       r_2nadn=one+zeta
       r_2nbdn=one-zeta
       zeta_eps=sqrt(DEPS)
       where(r_2nadn<zeta_eps)r_2nadn=zeta_eps
       where(r_2nbdn<zeta_eps)r_2nbdn=zeta_eps
       ! everywhere 1.0-zeta, 1.0+zeta and 1.0-zeta**2
       ! to be replaced with their expressions through r_2n?dn

       gr2 = gamma(1:vl,1)+gamma(1:vl,2)+2.0_r8_kind*&
            gamma(1:vl,3)
       gr = sqrt(gr2)
       g_ = half*((r_2nadn)**f2thrd+(r_2nbdn)**f2thrd)
       g2 = g_**2
       g3 = g2*g_
       g4 = g2**2
       z_a =  (r_2nbdn)/sums
       z_b = -(r_2nadn)/sums 

       t_pw = t_fact*GR/(sqrt(RhoTHRD)*sums)/g_
       t2= t_pw*t_pw
       t_g = half*t_fact/(sqrt(RhoTHRD)*sums*GR+deps)/g_
       t_n = -F7o6 * t_pw/sums
       t_gi = -t_pw/g_
       ex_p0 = -two*AL_PW*eps_lda/B_PW**2/g3
       Vexp = exp(ex_p0)
       A_den = (Vexp-ONE)
       where(A_den<deps) A_den = deps
       A_PW = A_PW_F/A_den
       A_PWt2 = A_PW*t2
       OA_PWt2 = one+A_PWt2

       A_gi =  THREE* (A_PW_F/A_den**2)*Vexp*ex_p0*g3/g4
       eps_lda_na= ( dfdrho_lda(1:vl,1)-eps_lda(1:vl) )/sums
       eps_lda_nb= ( dfdrho_lda(1:vl,2)-eps_lda(1:vl) )/sums
       dg= third  * (one/(r_2nadn)**(third)-one/(r_2nbdn)**(THIRD) )
       A_na = (A_PW_F**2/(B_PW*g3))*Vexp/A_den**2*eps_lda_na &
            +  A_gi*dg*z_a

       A_nb = (A_PW_F**2/(B_PW*g3))*Vexp/A_den**2*eps_lda_nb &
            +  A_gi*dg*z_b


       H0 = g3*H0_fact*log(one+A_PW_F*(t2+A_PWt2*t2)/(OA_PWt2+A_PWt2**2))

       F0_t = g3*sums*two*B_PW*t_pw*(OA_PWt2+A_PWt2) &
            /((one+(A_PW_F+A_PW)*t2*OA_PWt2)*(OA_PWt2+A_PWt2**2))
       h0_t = F0_t/sums
       F0_g = F0_t*t_g


       t6 = t2**3
       H0_a = -g3*B_PW*t6*A_PW*(two+A_PWt2) &
            /((one+(A_PW_F+A_PW)*t2*OA_PWt2)*(OA_PWt2+A_PWt2**2))

       h0_gi = (THREE*g2/g3)*H0 

       h0_na = h0_gi*dg*z_a + H0_a*A_na &
            + h0_t*(t_n+t_gi*dg*z_a)

       h0_nb = h0_gi*dg*z_b + H0_a*A_nb &
            + h0_t*(t_n+t_gi*dg*z_b)

       f0_na = h0 + sums*h0_na
       f0_nb = h0 + sums*h0_nb
       ex_p1= - g4*e_fact*t2/RhoTHRD
       rg_nom = e_rg+a_rg*RS+b_rg*RS2
       rg_dnom = one+c_rg*RS+d_rg*RS2+ten*b_rg*RS2*RS
       Cc_rs = thnd*rg_nom/rg_dnom-Cx

       rg_n_ =a_rg+two*b_rg*RS
       rg_dn_=c_rg+d_rg*two*RS+three*ten*b_rg*RS2
       H1 =    g3*H1_fact*(Cc_rs-Cc0-Cx37)*t2*exp(ex_p1)
       H1_t = g3*H1_fact*(Cc_rs-Cc0-Cx37)   *exp(ex_p1) &
            *two*t_pw*(one+ex_p1)
       H1_gi = (THREE*g2/g3)*h1+h1*(four*g3/g4)*ex_p1

       H_rs = h1_fact*t2*exp(ex_p1) &
            *thnd*(rg_n_*rg_dnom-rg_dn_*rg_nom)/rg_dnom**2
       RS_n  = -THIRD*RS/sums


       H1_dn=  g4*H1*THIRD*t2*e_fact/(sums*RhoTHRD)

       H1_na = H1_t* t_gi* dg*z_a + H1_t * t_n & 
            + H1_gi* dg*z_a &
            + H_rs*RS_n &
            + H1_dn

       H1_nb = H1_t* t_gi* dg*z_b &
            + H1_t * t_n  & 
            + H1_gi* dg*z_b &
            + H_rs*RS_n &
            + H1_dn

       F1_na = H1+ sums*H1_na
       F1_nb = H1+ sums*H1_nb

       F1_g = sums*H1_t*t_g

       if( revision )then
          dz2 = four/(sum4+yeps)* &
               (gamma(1:vl,2)*rho(1:vl,1)**2 &
               +gamma(1:vl,1)*rho(1:vl,2)**2 &
               - two*rho(1:vl,2)*rho(1:vl,1)*gamma(1:vl,3)) 

          dz2_aa = (four/(sum4+yeps))*rho(1:vl,2)**2
          dz2_bb = (four/(sum4+yeps))*rho(1:vl,1)**2
          dz2_ab = - (eight/(sum4+yeps))*rho(1:vl,2)*rho(1:vl,1)

          dz2_a_s = four/(sum4+yeps) &
               *( ( two*gamma(1:vl,2)*rho(1:vl,1)-&
               two*rho(1:vl,2)*gamma(1:vl,3))*sums &
               - four*(gamma(1:vl,1)*rho(1:vl,2)**2 &
               +gamma(1:vl,2)*rho(1:vl,1)**2 &
               - two*rho(1:vl,1)*rho(1:vl,2)*gamma(1:vl,3)) )

          dz2_b_s = four/(sum4+yeps) &
               *( ( two*gamma(1:vl,1)*rho(1:vl,2)-&
               two*rho(1:vl,1)*gamma(1:vl,3))*sums &
               - four*(gamma(1:vl,1)*rho(1:vl,2)**2 &
               +gamma(1:vl,2)*rho(1:vl,1)**2 &
               - two*rho(1:vl,1)*rho(1:vl,2)*gamma(1:vl,3)) )



          dz_dn = two &
               *(rho(1:vl,2)*gamma(1:vl,1)-&
               rho(1:vl,1)*gamma(1:vl,2)-&
               gamma(1:vl,3)*sums*zeta)/(sum2+yeps)

          dz_dn_b =    (two / (sum3+yeps)) &
               *((gamma(1:vl,1) - gamma(1:vl,3)*zeta - &
               gamma(1:vl,3)*sums*z_b)*sums &
               - two*(rho(1:vl,2)*gamma(1:vl,1)&
               -rho(1:vl,1)*gamma(1:vl,2)-gamma(1:vl,3)*sums*zeta))

          dz_dn_a = - (two / (sum3+yeps)) &
               *((gamma(1:vl,2) + gamma(1:vl,3)*zeta + &
               gamma(1:vl,3)*sums*z_a)*sums &
               + two*(rho(1:vl,2)*gamma(1:vl,1)&
               -rho(1:vl,1)*gamma(1:vl,2)-gamma(1:vl,3)*sums*zeta)) 

          dz_dn_aa =  (two/(sum2+yeps))*rho(1:vl,2)
          dz_dn_bb = -(two/(sum2+yeps))*rho(1:vl,1)
          dz_dn_ab = -(two/(sums + yeps))*zeta

          h21= c_21*zeta/(r_2nadn*r_2nbdn)**THIRD
          h22= (c_22+c_23*zeta**2)/(r_2nadn*r_2nbdn)
          dh21= (c_21/(r_2nadn*r_2nbdn)**(F2THRD)) *((r_2nadn*r_2nbdn)**(THIRD) &
               -F2THRD*zeta**2/(r_2nadn*r_2nbdn)**(F2THRD) )
          dh22=(two*zeta/(r_2nadn*r_2nbdn)**2)*(c_23*(r_2nadn*r_2nbdn)+(c_22+c_23*zeta**2))
          H2 = (Cc0/RhoTHRD)*((h21/(sums + yeps))*dz_dn +h22*dz2)
          !                             ^                          ^
          !     two summed terms are of zero*infinity type when |zeta|==1.0

          F2 = H2*sums


          f2_nd = - (Cc0/RhoTHRD)*(h21/(sums + yeps))*dz_dn -THIRD*H2

          f2_h21a = (Cc0/RhoTHRD)*dz_dn*dh21*z_a
          f2_h21b = (Cc0/RhoTHRD)*dz_dn*dh21*z_b

          f2_dz_dn_a = (Cc0/RhoTHRD)*h21*dz_dn_a
          f2_dz_dn_b = (Cc0/RhoTHRD)*h21*dz_dn_b

          f2_h22a = (Cc0/RhoTHRD)*dz2*dh22*z_a*sums
          f2_h22b = (Cc0/RhoTHRD)*dz2*dh22*z_b*sums

          f2_dz2_a = (Cc0/RhoTHRD)*h22*dz2_a_s
          f2_dz2_b = (Cc0/RhoTHRD)*h22*dz2_b_s


          f2_na = h2 +  f2_nd  + f2_dz_dn_a + f2_h21a + f2_dz2_a  + f2_h22a
          f2_nb = h2  + f2_nd  + f2_dz_dn_b + f2_h21b + f2_dz2_b  + f2_h22b
          h2_g_aa = ( (h21/(sums + yeps))*dz_dn_aa+h22*dz2_aa ) &
               * Cc0   /RhoTHRD
          h2_g_bb = ( (h21/(sums + yeps))*dz_dn_bb+h22*dz2_bb ) &
               * Cc0   /RhoTHRD
          h2_g_ab = ( (h21/(sums + yeps))*dz_dn_ab+h22*dz2_ab ) &
               * Cc0   /RhoTHRD
       else ! i.e.: if revision == .false.
          h2=zero
          f2_na=zero
          f2_nb=zero
          h2_g_aa=zero
          h2_g_bb=zero
          h2_g_ab=zero
       endif! revision?
       fxc(1:vl) = fxc(1:vl) + (h0+h1+h2)*sums
       ! spin down
       dfdrho(1:vl,1) = dfdrho(1:vl,1) + f0_na + f1_na + f2_na
       ! spin up
       dfdrho(1:vl,2) = dfdrho(1:vl,2) + f0_nb + f1_nb + f2_nb
       ! now derivatives with respect to gamma
       dfdgrarho(1:vl,1) = dfdgrarho(1:vl,1)&
             + F0_g + F1_g  + h2_g_aa*sums
       dfdgrarho(1:vl,2) = dfdgrarho(1:vl,2)&
            + F0_g + F1_g  + h2_g_bb*sums
       dfdgrarho(1:vl,3) = dfdgrarho(1:vl,3)&
            + (F0_g + F1_g)*two + h2_g_ab*sums

    end if! case: 'spin unrestricted' done

  end subroutine pw91c_calc

  !*************************************************************

  !*************************************************************
  subroutine pw91x_calc(rho,gamma,dfdrho,ispin,fxc,dfdgrarho,vl) 
    ! Purpose: calculate the correlation energy of PW91 functional
    ! Note: no LDA contributions are included, a call to a LDA 
    !       functional must preceed a call to this functional
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho  ! rho(vl,ispin)
    real(kind=r8_kind),dimension(:,:),intent(in) :: gamma ! gamma(vl,(ispin-1)*2+1)
    ! gamma(:,1)=grarho_down*grarho_down
    ! gamma(:,2)=grarho_up*grarho_up
    ! gamma(:,3)=grarho_down*grarho_up
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc ! f according to JGP
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho,dfdgrarho
    ! dfdrho(vl,ispin), dfdgrarho(vl,(ispin-1)*2+1)
    ! derivatives of rho with respect to rho and gamma
    integer(kind=i4_kind),intent(in) :: ispin,vl
    !** End of interface ***************************************** 

    real(kind=r8_kind), parameter ::  alpha   = 0.19645_r8_kind, &  
         beta    = 7.7956_r8_kind, &
         gamma_pw   = 0.004_r8_kind, &
         delta   = 0.2743_r8_kind, &
         epsilon = 0.1508_r8_kind

    real(kind=r8_kind),dimension(vl)  :: e_fact, s_fact, sums, ros, rhothrd, &
         rhof4th, grhos, s_pw, s2, xgr, a_s_pw, s4, s2_h, ars_xgr, pw_nom, pw_denom, &
         egr, ds_n, ds_g, d_nom_s, d_denom_s, df_s, df_n, df_g, egrup, &
         egrdn


    e_fact= -THREE/(FOUR*PI)*(THREE*PI*PI)**THIRD
    s_fact= one/(TWO*(THREE*PI*PI)**THIRD)

    if (ispin==1)then ! ** SPIN RESTRICTED CASE **
       sums =abs(RHO(1:vl,1)) + DEPS
       ROs = sums
       RhoTHRD = ROs**THIRD
       RhoF4TH = RhoTHRD*ROs
       GRHOs = sqrt(gamma(1:vl,1))
       s_pw  = s_fact*GRHOs/rhoF4TH
       s2 = s_pw*s_pw
       XGR  = s_pw*BETA 
       a_s_pw = alpha*s_pw
       s4 = s2*s2
       s2_h = - s2*hundt
       ARS_XGR = ARSINH(XGR)

       pw_nom = s2*(delta-epsilon*exp(s2_h))-gamma_pw*s4
       pw_denom = one+a_s_pw*ARS_XGR+gamma_pw*s4

       EGR  = e_fact*RhoTHRD*pw_nom/pw_denom

       fxc(1:vl)=fxc(1:vl)+egr*sums
       ds_n = -F4THRD*s_pw/ROs
       ds_g = half*s_fact/(RhoF4TH*GRHOs+deps)

       d_nom_s = -four*gamma_pw*s2*s_pw+two*delta*s_pw &
            -two*epsilon*s_pw*exp(s2_h)*(one+s2_h)

       d_denom_s = four*gamma_pw*s2*s_pw &
            + BETA*a_s_pw/sqrt(one+XGR*XGR) &
            + alpha*ARS_XGR

       df_s = (d_nom_s*pw_denom-pw_nom*d_denom_s)/pw_denom**2

       df_n = EGR*F4THRD + e_fact*RhoF4TH*df_s*ds_n 

       df_g = e_fact*RhoF4TH*df_s*ds_g

       dfdrho(1:vl,1) = dfdrho(1:vl,1)+df_n
       dfdgrarho(1:vl,1)=dfdgrarho(1:vl,1)+df_g

    else  ! (ISPIN.NE.1)   ! *** SPIN POLARIZED CASE ***

       sums = RHO(1:vl,1) + RHO(1:vl,2)+ DEPS

       ROs = rho(1:vl,1)+rho(1:vl,1)
       rhothrd = ROs**THIRD
       RhoF4TH = RhoTHRD*ROs

       GRHOs = two*sqrt(gamma(1:vl,1))

       s_pw  = s_fact*GRHOs/(RhoF4TH+deps)
       s2 = s_pw*s_pw
       XGR  = s_pw*BETA
       a_s_pw = alpha*s_pw
       s4 = s2*s2
       s2_h = - s2*hundt
       ARS_XGR = ARSINH(XGR)

       pw_nom = s2*(delta-epsilon*exp(s2_h))-gamma_pw*s4
       pw_denom = one+a_s_pw*ARS_XGR+gamma_pw*s4
       EGRup  = e_fact*RhoF4TH*pw_nom/pw_denom

       ds_n = -F4THRD*s_pw/(ROs+deps)
       ds_g = half*s_fact/(RhoF4TH*GRHOs+deps)

       d_nom_s = -four*gamma_pw*s2*s_pw+two*delta*s_pw &
            -two*epsilon*s_pw*exp(s2_h)*(one+s2_h)

       d_denom_s = four*gamma_pw*s2*s_pw &
            + BETA*a_s_pw/sqrt(one+XGR*XGR) &
            + alpha*ARS_XGR

       df_s = (d_nom_s*pw_denom-pw_nom*d_denom_s)/pw_denom**2
       df_n = F4THRD*e_fact*RhoTHRD*pw_nom/pw_denom &
            + e_fact*RhoF4TH*df_s*ds_n

       df_g = e_fact*RhoF4TH*df_s*ds_g

       dfdrho(1:vl,1) = dfdrho(1:vl,1) + df_n
       dfdgrarho(1:vl,1) = dfdgrarho(1:vl,1) + df_g*two


       ROs = two*rho(1:vl,2)
       RhoTHRD = ROs**THIRD
       RhoF4TH = RhoTHRD*ROs

       GRHOs = two*sqrt(gamma(1:vl,2))

       s_pw  = s_fact*GRHOs/(RhoF4TH+deps)
       s2 = s_pw*s_pw
       XGR  = s_pw*BETA
       a_s_pw = alpha*s_pw
       s4 = s2*s2
       s2_h = - s2*hundt
       ARS_XGR = ARSINH(XGR)

       pw_nom = s2*(delta-epsilon*exp(s2_h))-gamma_pw*s4
       pw_denom = one+a_s_pw*ARS_XGR+gamma_pw*s4
       EGRdn  = e_fact*RhoF4TH*pw_nom/pw_denom

       ds_n = -F4THRD*s_pw/(ROs+deps)
       ds_g = half*s_fact/(RhoF4TH*GRHOs+deps)

       d_nom_s = -four*gamma_pw*s2*s_pw+two*delta*s_pw &
            -two*epsilon*s_pw*exp(s2_h)*(one+s2_h)

       d_denom_s = four*gamma_pw*s2*s_pw &
            + BETA*a_s_pw/sqrt(one+XGR*XGR) &
            + alpha*ARS_XGR

       df_s = (d_nom_s*pw_denom-pw_nom*d_denom_s)/pw_denom**2
       df_n = F4THRD*e_fact*RhoTHRD*pw_nom/pw_denom &
            + e_fact*RhoF4TH*df_s*ds_n

       df_g = e_fact*RhoF4TH*df_s*ds_g

       dfdrho(1:vl,2) = dfdrho(1:vl,2) + df_n
       dfdgrarho(1:vl,2) = dfdgrarho(1:vl,2) + df_g*two
       EGR    = HALF* (  EGRup+EGRdn )
       fxc(1:vl)=fxc(1:vl)+ egr


    end if! case: 'spin unrestricted' done


  end subroutine pw91x_calc
  !********************************************************************

  !********************************************************************
  function arsinh(y)
    ! Purpose: calculate arsinh for a vector y
    real(kind=r8_kind),intent(in),dimension(:) :: y
    real(kind=r8_kind),dimension(size(y)) :: arsinh
    arsinh=log(y+sqrt(y*y+1.0_r8_kind))
    return
  end function arsinh
  !********************************************************************
end module perdew_wang_module


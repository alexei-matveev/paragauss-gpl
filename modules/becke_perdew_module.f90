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
module becke_perdew_module
!---------------------------------------------------------------
!
!  Purpose: The module calculates the xc-functionals
!           from becke and perdew
!           The calculation follows the sheme of
!  Reference: Johnson, Gill, Pople, J. Chem. Phys. 98 (7) 1 april 1993
!
!           The used variables are as follows:
!           rho:  density
!           gamma: square of the norm of the density gradients
!           dfdrho: derivative of with respect to rho
!           fxc:    f
!           dfdgrarho: derivative of f with respect to grarho
!
!  Module called by: xc_hamiltonian
!
!  Author: MS
!  Date: 2/96
!
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!------------ Modules used --------------------------------------
  use type_module
  implicit none
  private
  save
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
  public :: perdew_calc,becke88_calc,becke88_calcMDA,perdew_calcMDA


  !===================================================================
  ! End of public interface of module
  !===================================================================

  real(kind=r8_kind), parameter :: ZERO  = 0.0_r8_kind ,&
       ONE   = 1.0_r8_kind ,&
       TWO   = 2.0_r8_kind ,&
       THREE = 3.0_r8_kind ,&
       FOUR  = 4.0_r8_kind ,&
       FIVE  = 5.0_r8_kind ,&
       SIX   = 6.0_r8_kind ,&
       SEVEN = 7.0_r8_kind ,&
       NINE  = 9.0_r8_kind ,&
       HALF  = 0.5_r8_kind ,&
       DEPS  = 1.0E-30_r8_kind ,&
       DEPSMDA=1.0E-08_r8_kind ,&
       THIRD  = ONE/THREE ,&
       F2THRD = TWO/THREE ,&
       F4THRD = FOUR/THREE ,&
       F5THRD = FIVE/THREE ,&
       F7THRD = SEVEN/THREE ,&
       F8THRD =  8.0_r8_kind/THREE ,&
       F11THRD = 11._r8_kind/THREE ,&
       F10THRD = 10._r8_kind/THREE ,&
       F5sixth = FIVE/SIX ,&
       F7sixth = SEVEN/SIX ,&
       F4NINETH = FOUR/NINE ,&
       F5twelf = F5sixth/two ,&
       PI = 3.14159265358979324  ,&
       alpha = 0.023266_r8_kind,&
       beta = 7.389e-6_r8_kind ,&
       gammas=8.723_r8_kind ,&
       delta = 0.472_r8_kind, &
       f_wave= 0.11_r8_kind, &
       CFAC = 0.001667_r8_kind, &
       OCFAC = 0.002568_r8_kind, &
       rs_fac = 0.62035049089940001667_r8_kind , &
       c_infi = CFAC + OCFAC, &
       fi_fac = 1.745_r8_kind*f_wave*c_infi , &
       e4_beta = beta*1.0e4_r8_kind , &
       two_delta = TWO * delta , &
       two_beta = TWO * beta , &
       f3_e4_beta = THREE * e4_beta, &
       ! // TWO_oth = TWO ** THIRD
       TWO_oth = 1.25992104989487316476_r8_kind , &
       TWO_tth = TWO_oth **2 , &
       s4_fac = F5sixth * TWO_tth

 integer(kind=i4_kind), parameter:: uup=1,ddn=2,udn=3,dup=4,udup=5,uddn=6
 integer(kind=i4_kind), parameter:: aa=1,bb=2,ab=3,ac=4,bc=5,cc=6



contains


  subroutine perdew_calc(rho,gamma,dfdrho,ispin,fxc,dfdgrarho,vec_length, &
                         df_drhodrho,df_drhodgamma,df_dgammadgamma)
    ! purpose : Calculation of the perdew86 functional. Formulas can be
    !           found in L. FAN,PHdthesis, Calgary or J.P. Perdew
    !           Phys. Rev., B34, 7046(1986) or Phys. Rev. B 33, 8822
    !           (1986)
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho,gamma
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho,dfdgrarho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    real(kind=r8_kind),optional,intent(inout):: df_drhodrho(:,:),df_drhodgamma(:,:)
    real(kind=r8_kind),optional,intent(inout):: df_dgammadgamma(:,:)
    !** End of interface *****************************************
    real(kind=r8_kind),dimension(vec_length) :: sums,rothrd,rof4th,rs,rs2,rs3,&
         gr2,cden,ovcf,c_rho,fi,exp_fi,dp_g,t1,t2,dc_drho,dfi_rho,&
         xi,dm,ddm_xi,dp_rhod,&
         gr,dp_rho,dp_rhou!,dfi_g
    real(kind=r8_kind),dimension(vec_length) :: RS_r,OVCF_r,CDEN_r, &
         dP_rho_r, dC_dRho_r, t1_r, t2_r, dfi_rho_r, &
         fi_g, dP_rho_g, dfi_rho_g, dP_g_g !, dfi_g_g, exp_fi_r
    real(kind=r8_kind),dimension(vec_length) :: dP_g_restr, &
                                                dP_g_restr_r,dP_g_rhoup,dP_g_rhodn, &
                                                dm_dP_rho,dP_rhou_u,dP_rhou_d,dP_rhod_d, &
                                                xi_u,xi_d,ddm_xi_xi,xi_u_u,xi_d_d,xi_u_d,P

    if(ispin==1) then ! spin restricted case
       sums =RHO(:vec_length,1)+deps
    else
       SUMS = sum(RHO(:vec_length,:),2) + DEPS

       xi = (RHO(:vec_length,1) - RHO(:vec_length,2)) / SUMS
       xi_u= (one-xi)/SUMS
       xi_d=-(one+xi)/SUMS
       dm =  (one/TWO_oth)/sqrt( ((one+xi)/two)**F5THRD + ((one-xi)/two)**F5THRD )

       ddm_xi = F5twelf*dm**3 *TWO_oth**2*( ((one-xi)/two)**F2THRD - ((one+xi)/two)**F2THRD )
    endif

       RoTHRD = sums ** THIRD
       RoF4TH = sums * RoTHRD
       RS = rs_fac / RoTHRD
       RS2 = RS*RS
       RS3 = RS2 * RS
    if(ispin==1) then ! spin restricted case
       GR2 = gamma(:VEC_LENGTH,1)
    else
       GR2 =gamma(:vec_length,1)+gamma(:vec_length,2) &
           +2.0_r8_kind*gamma(:vec_length,3)
    endif

       gr=sqrt(gr2)
       CDEN = (ONE+gammas*RS+ delta*RS2 + e4_beta*RS3)
       OVCF = OCFAC  + alpha*RS + beta*RS2
       C_Rho = CFAC +OVCF/CDEN

       fi_g =  fi_fac/(C_Rho*sums**F7sixth)
       fi = fi_g*gr
       exp_fi = exp(-fi)

!       dfi_g = fi / ( 2 * GR2 + deps ) !fi_g/(2*gr+deps)
!       dP_g_restr = (exp_fi * C_Rho/ RoF4TH)*(one - dfi_g * GR2)

       dP_g_restr = (exp_fi * C_Rho/ RoF4TH)*(one - fi_g*GR/2)

    if(ispin==2) then
       dP_g =   dP_g_restr * dm
       !___________________________________________________________
#if 1
       dfdgrarho(:vec_length,1) = dfdgrarho(:vec_length,1) + dP_g
       dfdgrarho(:vec_length,2) = dfdgrarho(:vec_length,2) + dP_g
       dfdgrarho(:vec_length,3) = dfdgrarho(:vec_length,3) + dP_g * 2
!print*,'vxc dfdgrarho',sum(dfdgrarho(:vec_length,1)),sum(dfdgrarho(:vec_length,2)),sum(dfdgrarho(:vec_length,3))
#endif
    endif

       !  Vc_NL Eq. 4.34 PhD Thesis by Fan **
       t1 = ( alpha + two_beta*RS)/CDEN
       t2 = (spread(gammas,1,vec_length) &
            + two_delta *RS + f3_e4_beta*RS2) * OVCF/(CDEN*CDEN)

       dC_dRho = - THIRD * RS/sums * ( t1 - t2 )
         RS_r=-third*RS/sums
         OVCF_r=alpha*RS_r+beta*two*RS*RS_r
         CDEN_r= gammas*RS_r+delta*two*RS*RS_r+e4_beta*three*RS2*RS_r
!        C_Rho_r=OVCF_r/CDEN-OVCF*CDEN_r/CDEN**2

       dfi_rho_g= - (fi_fac/(C_Rho*sums**F7sixth)**2)* &
                    sqrt(RoTHRD)*(dC_dRho*sums+F7sixth*C_Rho)
       dfi_rho = dfi_rho_g*gr
       P=C_Rho/RoF4TH*exp_fi*GR2
       dP_rho=(exp_fi/RoF4TH)*(dC_dRho-F4THRD*C_Rho/sums-C_Rho*dfi_rho)*GR2

    if(ispin==1) then
       dfdrho(:vec_length,1) = dfdrho(:vec_length,1) + dP_rho ! potential
       !__________________________________________________________________
       dfdgrarho(1:vec_length,1) = dfdgrarho(1:vec_length,1) + dP_g_restr
       fxc(1:vec_length)=fxc(1:vec_length)+P

    else
       dm_dP_rho = dP_rho * dm
       dP_rhou = dm_dP_rho + ddm_xi*xi_u * P
       dP_rhod = dm_dP_rho + ddm_xi*xi_d * P

       !_______________________________________________________
       dfdrho(:vec_length,1) = dfdrho(:vec_length,1) + dP_rhou
       dfdrho(:vec_length,2) = dfdrho(:vec_length,2) + dP_rhod
       fxc(1:vec_length) = fxc(1:vec_length)+P*dm
    endif

   if(present(df_drhodrho)) then

    if(ispin==2) then
     dP_g_restr_r=-dP_g_restr*dfi_rho
     dP_g_restr_r=dP_g_restr_r+(exp_fi*dC_dRho/RoF4TH)*(one - fi_g * GR/2)
     dP_g_restr_r=dP_g_restr_r-dP_g_restr*F4THrd/sums
     dP_g_restr_r=dP_g_restr_r-(exp_fi * C_Rho/ RoF4TH)*dfi_rho_g*GR/2

     dP_g_rhoup=dP_g_restr_r*dm+dP_g_restr*ddm_xi*xi_u
     dP_g_rhodn=dP_g_restr_r*dm+dP_g_restr*ddm_xi*xi_d

     !_______________________________________________________________________
     df_drhodgamma(:vec_length,uup)=df_drhodgamma(:vec_length,uup)+dP_g_rhoup
     df_drhodgamma(:vec_length,ddn)=df_drhodgamma(:vec_length,ddn)+dP_g_rhodn
     df_drhodgamma(:vec_length,udn)=df_drhodgamma(:vec_length,udn)+dP_g_rhodn
     df_drhodgamma(:vec_length,dup)=df_drhodgamma(:vec_length,dup)+dP_g_rhoup
     df_drhodgamma(:vec_length,udup)=df_drhodgamma(:vec_length,udup)+dP_g_rhoup*2
     df_drhodgamma(:vec_length,uddn)=df_drhodgamma(:vec_length,uddn)+dP_g_rhodn*2


    else ! i.e. restricted
!!       dP_rho = (GR2 * exp_fi / RoF4TH)*( dC_dRho - F4THRD*C_Rho/sums-C_Rho*dfi_rho )
       dP_rho_g=two*gr*exp_fi/RoF4TH*( dC_dRho - F4THRD*C_Rho/sums-C_Rho*dfi_rho ) !!(1) d->GR2
       dP_rho_g=dP_rho_g-dP_rho*fi_g                                               !!(2)d->exp_fi
       dP_rho_g=dP_rho_g-(GR2 * exp_fi / RoF4TH)*C_Rho*dfi_rho_g                   !!(3)d->dfi_rho
    !! ______________________________________________________________________________
       df_drhodgamma(:vec_length,1)=df_drhodgamma(:vec_length,1)+dP_rho_g*half/(gr+deps)
    endif

   endif

    if(present(df_drhodrho)) then
!!        exp_fi_r=-exp_fi*fi_r
!!        dP_rho=(GR2*exp_fi/RoF4TH)*( dC_dRho - F4THRD*C_Rho/sums-C_Rho*dfi_rho )
        dP_rho_r=-dfi_rho*dP_rho                                       !!(1) - exp_fi_r contrib
        dP_rho_r=dP_rho_r-F4THrd/sums*dP_rho                       !!(2) - RoF4TH_r contrib
!!       dC_dRho = - THIRD * RS/sums * ( t1 - t2 )
        dC_dRho_r=- THIRD * RS_r/sums * ( t1 - t2 )
        dC_dRho_r=dC_dRho_r-dC_dRho/sums
!!       t1 = ( alpha + two_beta*RS)/CDEN
         t1_r=(two_beta*RS_r-CDEN_r*t1)/CDEN
!!       t2 = ( spread(gammas,1,vec_length) &
!!            + two_delta *RS + f3_e4_beta*RS2 ) * OVCF/(CDEN*CDEN)
         t2_r= two_delta*RS_r*OVCF/(CDEN*CDEN)
         t2_r=t2_r+f3_e4_beta*two*RS*RS_r*OVCF/(CDEN*CDEN)
         t2_r=t2_r+( spread(gammas,1,vec_length) &
            + two_delta *RS + f3_e4_beta*RS2 ) * OVCF_r/(CDEN*CDEN)
         t2_r=t2_r-two*CDEN_r*t2/CDEN
         dC_dRho_r=dC_dRho_r-THIRD*RS/sums*( t1_r-t2_r)
        dP_rho_r=dP_rho_r+(GR2*exp_fi/RoF4TH)*dC_dRho_r                !!(3) - dC_dRho_r contrib
        dP_rho_r=dP_rho_r-(GR2 * exp_fi / RoF4TH)*F4THRD*dC_dRho/sums  !!(4) -C_Rho_r contrib
        dP_rho_r=dP_rho_r+(GR2 * exp_fi / RoF4TH)*F4THRD*C_Rho/sums**2 !!(5)
        dP_rho_r=dP_rho_r-(GR2 * exp_fi / RoF4TH)*dC_dRho*dfi_rho      !!(6)

!!       dfi_rho = - (fi_fac*gr/(C_Rho*sums**F7sixth)**2)&
!!            *sqrt(RoTHRD)*(dC_dRho*sums+F7sixth*C_Rho)
         dfi_rho_r = - two*dfi_rho/(C_Rho*sums**F7sixth)* &
                       (dC_dRho*sums**F7sixth+F7sixth*C_Rho*sqrt(RoTHRD))
         dfi_rho_r=dfi_rho_r+dfi_rho/sums/six
         dfi_rho_r=dfi_rho_r-(fi_fac*gr/(C_Rho*sums**F7sixth)**2)&
            *sqrt(RoTHRD)*(dC_dRho_r*sums+dC_dRho+F7sixth*dC_dRho)

         dP_rho_r=dP_rho_r-(GR2 * exp_fi / RoF4TH)*C_Rho*dfi_rho_r      !!(7)

!!       dP_rho = (GR2 * exp_fi / RoF4TH)*( dC_dRho - F4THRD*C_Rho/sums-C_Rho*dfi_rho )
    if(ispin==1) then ! spin restricted case
    !!_______________________________________________________________
       df_drhodrho(:vec_length,1)=df_drhodrho(:vec_length,1)+dP_rho_r
    else
!!!       dP_rhou = dm_dP_rho + ddm_xi*xi_u * P
!!!       dP_rhod = dm_dP_rho + ddm_xi*xi_d * P
!!!       ddm_xi = F5twelf*dm**3 *TWO_oth**2*( ((one-xi)/two)**F2THRD - ((one+xi)/two)**F2THRD )

       ddm_xi_xi = 3*F5twelf*dm**2*TWO_oth**2*( ((one-xi)/two)**F2THRD - ((one+xi)/two)**F2THRD )*ddm_xi &
      + F5twelf*dm**3 *TWO_oth**2*F2THRD*( -(two/(one-xi) )**third -(two/(one+xi))**third )/2
        xi_u_u=-2*xi_u/sums
        xi_d_d=-2*xi_d/sums
        xi_u_d=-(xi_d+xi_u)/sums
       dP_rhou_u=dP_rho_r*dm+dP_rho*ddm_xi*xi_u*2+(ddm_xi_xi*xi_u**2+ddm_xi*xi_u_u)*P
       dP_rhou_d=dP_rho_r*dm+dP_rho*ddm_xi*(xi_d+xi_u)+(ddm_xi_xi*xi_d*xi_u+ddm_xi*xi_u_d)*P
       dP_rhod_d=dP_rho_r*dm+dP_rho*ddm_xi*xi_d*2+(ddm_xi_xi*xi_d**2+ddm_xi*xi_d_d)*P
       !______________________________________________________________
       df_drhodrho(:vec_length,1)=df_drhodrho(:vec_length,1)+dP_rhou_u
       df_drhodrho(:vec_length,2)=df_drhodrho(:vec_length,2)+dP_rhod_d
       df_drhodrho(:vec_length,3)=df_drhodrho(:vec_length,3)+dP_rhou_d
    endif


!!     dP_g_restr = (exp_fi * C_Rho/ RoF4TH)*(one - fi_g*GR/2)
       dP_g_g=-dP_g_restr*fi_g
       dP_g_g=dP_g_g-(exp_fi * C_Rho/ RoF4TH)*(fi_g/2)
     if(ispin==1) then ! spin restricted case
!!      ______________________________________________________
       df_dgammadgamma(:vec_length,1)=df_dgammadgamma(:vec_length,1)+dP_g_g*half/(gr+deps)
     else
       ! dP_g=dm*dP_g_restr
       !______________________________________________________________________________________
       df_dgammadgamma(:vec_length,aa)=df_dgammadgamma(:vec_length,aa)+dP_g_g*half/(gr+deps)*dm
       df_dgammadgamma(:vec_length,bb)=df_dgammadgamma(:vec_length,bb)+dP_g_g*half/(gr+deps)*dm
       df_dgammadgamma(:vec_length,ab)=df_dgammadgamma(:vec_length,ab)+dP_g_g*half/(gr+deps)*dm
       df_dgammadgamma(:vec_length,ac)=df_dgammadgamma(:vec_length,ac)+dP_g_g/(gr+deps)*dm
       df_dgammadgamma(:vec_length,bc)=df_dgammadgamma(:vec_length,bc)+dP_g_g/(gr+deps)*dm
       df_dgammadgamma(:vec_length,cc)=df_dgammadgamma(:vec_length,cc)+dP_g_g*two/(gr+deps)*dm
     endif

    endif
    end subroutine perdew_calc

  subroutine perdew_calcMDA(rho,gamma,dfdrho,ispin,fxc,dfdgrarho,vec_length)
    ! purpose : Calculation of the perdew86 functional. Formulas can be
    !           found in L. FAN,PHdthesis, Calgary or J.P. Perdew
    !           Phys. Rev., B34, 7046(1986) or Phys. Rev. B 33, 8822
    !           (1986)
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho,gamma
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho,dfdgrarho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    !** End of interface *****************************************
    real(kind=r8_kind),dimension(vec_length) :: sums,rothrd,rof4th,rs,rs2,rs3,&
         gr2,cden,ovcf,c_rho,fi,exp_fi,dp_g,t1,t2,dc_drho,dfi_rho,&
         xi,dm,ddm_xi,dfi_gud,dfi_g,dp_gud,dp_rhod,&
         rof8th,gr,fi2,dp_rho,dp_rhou



    if(ispin==1) then ! spin restricted case

       sums =RHO(1:vec_length,1)+depsMDA   !(1)
       RoTHRD = sums ** THIRD
       RoF4TH = sums * RoTHRD
       RS = rs_fac / RoTHRD
       RS2 = RS*RS
       RS3 = RS2 * RS
       GR2 = gamma(1:VEC_LENGTH,1)
       gr=sqrt(gr2)
       CDEN = (ONE+gammas*RS+ delta*RS2 + e4_beta*RS3)
       OVCF = OCFAC  + alpha*RS + beta*RS2

       C_Rho = CFAC +OVCF/CDEN

       fi = (fi_fac/(C_Rho*sums**F7sixth))*gr
       exp_fi = exp(-fi)
       dfi_g = half*(fi_fac/(C_Rho*sums**F7sixth))/(gr+DEPSmda)   !(2)
       dP_g = (exp_fi * C_Rho/ RoF4TH)*(one - dfi_g * GR2)
       !  Vc_NL Eq. 4.34 PhD Thesis by Fan **
       t1 = ( alpha + two_beta*RS)/CDEN
       t2 = (spread(gammas,1,vec_length) &
            + two_delta *RS + f3_e4_beta*RS2) * OVCF/(CDEN*CDEN)

       dC_dRho = - THIRD * RS/sums * ( t1 - t2 )
       dfi_rho = - (fi_fac*gr/(C_Rho*sums**F7sixth)**2)&
            *sqrt(RoTHRD)*(dC_dRho*sums+F7sixth*C_Rho)

       dP_rho = (GR2 * exp_fi / RoF4TH)&
            *( dC_dRho - F4THRD*C_Rho/sums-C_Rho*dfi_rho )
       ! collecting the results
       dfdrho(1:vec_length,1) = dfdrho(1:vec_length,1) + dP_rho ! potential
       dfdgrarho(1:vec_length,1) = dfdgrarho(1:vec_length,1) + dP_g
       fxc(1:vec_length)=fxc(1:vec_length)+exp_fi * C_Rho * GR2 / RoF4TH

    end if

    if (ispin==2) then  ! spinpolarized case

       SUMS = sum(RHO(1:vec_length,:),2) + DEPSmda

       xi = (RHO(1:vec_length,1) - RHO(1:vec_length,2)) / SUMS
       dm =  (one/TWO_oth)/sqrt(((one+xi)/two)**F5THRD+((one-xi)/two)**F5THRD)

       ddm_xi = F5twelf*dm**3 *TWO_oth**2* &
                (((one-xi)/two)**F2THRD - ((one+xi)/two)**F2THRD )

       RoTHRD = sums ** THIRD
       RoF4TH = sums * RoTHRD
       RoF8TH = RoF4TH * RoF4TH
       RS = rs_fac / RoTHRD
       RS2 = RS*RS
       RS3 = RS2 * RS
       GR2 =gamma(1:vec_length,1)+gamma(1:vec_length,2)+2.0_r8_kind*&
            gamma(1:vec_length,3)
       GR=sqrt(GR2)
       CDEN = (ONE+gammas*RS+ delta*RS2 + e4_beta*RS3)
       OVCF = OCFAC  + alpha*RS + beta*RS2

       C_Rho = CFAC +OVCF/CDEN

       fi = (fi_fac/(C_Rho*sums**F7sixth))*GR
       dfi_gud = (fi_fac/(C_Rho*sums**F7sixth))/(GR+depsMDA)
       dfi_g = half*dfi_gud
       fi2 = fi**2
       exp_fi = exp(-fi)

       dP_gud = dm * (exp_fi * C_Rho/ RoF4TH)*(two - dfi_gud * GR2)
       dP_g =   dm * (exp_fi * C_Rho/ RoF4TH)*(one - dfi_g * GR2)
       dfdgrarho(1:vec_length,1) = dfdgrarho(1:vec_length,1) + dP_g
       dfdgrarho(1:vec_length,2) = dfdgrarho(1:vec_length,2) + dP_g
       dfdgrarho(1:vec_length,3) = dfdgrarho(1:vec_length,3) + dP_gud
       fxc(1:vec_length) = fxc(1:vec_length)+&
            dm * exp_fi * C_Rho * GR2 / RoF4TH

       !   **    Vc_NL Eq. 4.34 PhD Thesis by Fan **
       t1 = ( alpha + two_beta*RS)/CDEN
       t2 = (spread(gammas,1,vec_length)&
            + two_delta *RS + f3_e4_beta*RS2) * OVCF/(CDEN*CDEN)
       dC_dRho = - THIRD * RS * ( t1 - t2 ) / sums

       dfi_rho = - (fi_fac*GR/(C_Rho*sums**F7sixth)**2)&
            *sqrt(RoTHRD)*(dC_dRho*sums+F7sixth*C_Rho)

       dP_rho = dm * (GR2 * exp_fi / RoF4TH) &
            *( dC_dRho - F4THRD*C_Rho/sums-C_Rho*dfi_rho )
       dP_rhou = dP_rho + ddm_xi*(rho(1:vec_length,1)-xi)/sums &
            * exp_fi * C_Rho * GR2 / RoF4TH
       dP_rhod = dP_rho - ddm_xi*(rho(1:vec_length,2)+xi)/sums &
            * exp_fi * C_Rho * GR2 / RoF4TH

       dfdrho(1:vec_length,1) = dfdrho(1:vec_length,1) + dP_rhou
       dfdrho(1:vec_length,2) = dfdrho(1:vec_length,2) + dP_rhod


        ENDIF  ! case: 'spin unrestricted' done

      RETURN
    end subroutine perdew_calcMDA


    subroutine becke88_calc(rho,gamma,dfdrho,ispin,fxc,dfdgrarho,vec_length)
      ! purpose : Calculation of the Becke88 functional. Formulas can be
      !           found in A.D.Becke, Phys. Rev. A 38, 3098 (1988) or
      !           Johnson, Gill, Pople J. Chem. Phys. 98 (7) 1 april 1993
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho,gamma
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho,dfdgrarho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    !** End of interface *****************************************
    real(kind=r8_kind) :: two13
    real(kind=r8_kind),dimension(vec_length) :: sums,xgr,degr,g,egr,bx,dg,&
         dgm,dro,xgrup,xgrdn,degrup,degrdn,gup,gdn,bxup,bxdn,dgup,dgdn,&
         droup,drodn,dgmup,dgmdn
    real(kind=r8_kind),dimension(vec_length,ispin) :: gr,arsinhxgr

    two13=two**third

    if(ispin==1) then ! spin restricted

       sums=rho(1:vec_length,1)+deps
       gr(:,1)=sqrt(gamma(1:vec_length,1))

       XGR  = TWO13*gr(:,1)  /sums**F4THRD
       arsinhxgr(:,1)=arsinh(xgr)
       DEGR = ONE + SIX*0.0042_r8_kind*XGR*ARSINHXGR(:,1)
       g = 0.0042_r8_kind*XGR**2/DEGR
       EGR  = g/TWO13*sums**THIRD

       bX = 0.0042_r8_kind*XGR
       dG = ( SIX*bX**2* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) )- two*bX ) &
            / DEGR**2

       dGM =  HALF*dG/(gr(:,1)+deps)
       dRo =  F4THRD*rho(1:vec_length,1)**THIRD * (G+XGR*dG)/TWO13

       ! results
       dfdgrarho(1:vec_length,1) = dfdgrarho(1:vec_length,1)+dgm
       fxc(1:vec_length)=fxc(1:vec_length)-egr*sums
       dfdrho(1:vec_length,1)=dfdrho(1:vec_length,1)-dro

    else    ! spinpolarized case
       sums = rho(1:vec_length,1) + rho(1:vec_length,1)+ DEPS
       gr(:,1:2)=sqrt(gamma(1:vec_length,1:2))
       XGRUP  = GR(:,1)  /(RHO(1:vec_length,1) + DEPS)**F4THRD
       XGRDN  = GR(:,2)  /(RHO(1:vec_length,2) + DEPS)**F4THRD
       arsinhxgr(1:vec_length,1)=arsinh(xgrup)
       arsinhxgr(1:vec_length,2)=arsinh(xgrdn)
       DEGRUP = ONE + SIX*0.0042_r8_kind*XGRUP*ARSINHXGR(:,1)
       DEGRDN = ONE + SIX*0.0042_r8_kind*XGRDN*ARSINHXGR(:,2)
       Gup = 0.0042_r8_kind*XGRUP**2/DEGRUP
       Gdn = 0.0042_r8_kind*XGRDN**2/DEGRDN
       EGR    = rho(1:vec_length,1)**F4THRD*Gup+rho(1:vec_length,2)&
            **F4THRD*Gdn
       bXup = 0.0042_r8_kind*XGRUP
       bXdn = 0.0042_r8_kind*XGRDN
       dGup = ( six*bXup**2* ( XGRUP/sqrt(one+XGRUP**2)- ARSINHXGR(:,1) )- &
            two*bXup )/ DEGRUP**2
       dGdn = ( SIX*bXdn**2* ( XGRdn/sqrt(one+XGRdn**2) &
            - ARSINHXGR(:,2) ) - two*bXdn )/ DEGRdn**2
       dROup =  F4THRD*rho(1:vec_length,1)**THIRD* (Gup+XGRUP*dGup)
       dROdn =  F4THRD*rho(1:vec_length,2)**THIRD* (Gdn+XGRdn*dGdn)
       dGMup =  HALF*dGup/(GR(:,1)+deps)
       dGMdn =  HALF*dGdn/(GR(:,2)+deps)

       ! results
       fxc(1:vec_length) = fxc(1:vec_length) - egr
       dfdrho(1:vec_length,1) = -dROup+dfdrho(1:vec_length,1)
       dfdrho(1:vec_length,2) = -dROdn+dfdrho(1:vec_length,2)
       dfdgrarho(1:vec_length,1) = dGMup+dfdgrarho(1:vec_length,1)
       dfdgrarho(1:vec_length,2) = dGMdn+dfdgrarho(1:vec_length,2)
     end if
  end subroutine becke88_calc

    subroutine becke88_calcMDA(rho,gamma,dfdrho,ispin,fxc,dfdgrarho,vec_length,&
                      df_drhodrho,df_drhodgamma,df_dgammadgamma)
      ! purpose : Calculation of the Becke88 functional. Formulas can be
      !           found in A.D.Becke, Phys. Rev. A 38, 3098 (1988) or
      !           Johnson, Gill, Pople J. Chem. Phys. 98 (7) 1 april 1993
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho,gamma
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho,dfdgrarho
    real(kind=r8_kind),optional,intent(inout):: df_drhodrho(:,:), &
                  df_drhodgamma(:,:)
    real(kind=r8_kind),optional,intent(inout):: df_dgammadgamma(:,:)
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    !** End of interface *****************************************
    real(kind=r8_kind) :: two13
    real(kind=r8_kind),dimension(vec_length) :: sums,xgr,degr,g,egr,bx,dg,&
         dgm,dro,xgrup,xgrdn,degrup,degrdn,gup,gdn,bxup,bxdn,dgup,dgdn,&
         droup,drodn,dgmup,dgmdn
    real (r8_kind), dimension (vec_length) :: y, dydx, dydxdx, dxdr, dxdrdr, &
         ddegrdx, ddegrdxdx, dxdg, dgdx, dG_dx, dgdx_g, dg_g, dgmdx, dgmdx_g, &
         dxdgdr, dGM_g
    real(kind=r8_kind),dimension(vec_length,ispin) :: gr,arsinhxgr

    two13=two**third

       sums = sum(rho(:vec_length,:),2) + DEPSmda
    if(ispin==1) then ! spin restricted
       gr(:,1)=sqrt(gamma(:vec_length,1))
       XGR  = TWO13*gr(:,1)  /sums**F4THRD
       arsinhxgr(:,1)=arsinh(xgr)
       DEGR = ONE + SIX*0.0042_r8_kind*XGR*ARSINHXGR(:,1)
       g = 0.0042_r8_kind*XGR**2/DEGR
       EGR  = g/TWO13*sums**THIRD

       bX = 0.0042_r8_kind*XGR

       y=XGR+sqrt(one+XGR**2)
       dydx=(one+XGR/sqrt(one+XGR**2))

       ddegrdx=SIX*0.0042_r8_kind*ARSINHXGR(:,1)+SIX*bX/y*dydx

       dgdx=(two-XGR/degr*ddegrdx)*bx/degr
       dxdg=TWO13/sums**F4THRD

       !!!  df/dg
       dG = ( SIX*bX**2*( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) )- two*bX ) / DEGR**2

       dG_g= (SIX*bX*0.0042_r8_kind*dxdg* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) ) &
              - two*0.0042_r8_kind*dxdg)/DEGR**2

       dGM=dg_g*HALF
       dGM_g=half*(SIX*0.0042_r8_kind**2*dxdg**2*(XGR/sqrt(one+XGR**2)-ARSINHXGR(:,1) ) &
              - two*0.0042_r8_kind*dxdg/(gr(:vec_length,1)+depsMDA ))/DEGR**2

       dRo =  F4THRD*rho(:vec_length,1)**THIRD* (G+XGR*dG)/TWO13  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! results
       dfdgrarho(:vec_length,1) = dfdgrarho(:vec_length,1)+dgm
       fxc(:vec_length)=fxc(:vec_length)-egr*sums
       dfdrho(:vec_length,1)=dfdrho(:vec_length,1)-dro

       if(present(df_drhodrho)) then

        dydxdx=one/sqrt(one+XGR**2)-XGR**2/sqrt(one+XGR**2)/(one+XGR**2)
        dxdr=-TWO13*F4THRD*gr(:,1)/sums**f7THRD
        dxdrdr=TWO13*28_r8_kind/nine*gr(:,1)/sums**F10THRD
        ddegrdxdx=SIX*0.0042_r8_kind*(two/y*dydx-XGR*(dydx/y)**2+XGR/y*dydxdx )
        dxdgdr=-TWO13*F4THRD/sums**F7THRD

       dgdx_g=(two-XGR/degr*ddegrdx)*0.0042_r8_kind*dxdg/degr


       dG_dx=( SIX*TWO*bX*0.0042_r8_kind * ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) ) ) &
             / DEGR**2 &
            +( SIX*bX**2* (one/sqrt(one+XGR**2)-XGR**2/sqrt(one+XGR**2)/(one+XGR**2)  &
               - one/y*dydx))/DEGR**2 &
            - two*0.0042_r8_kind/DEGR**2 &
            - two*( SIX*bX**2* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) )- two*bX )/DEGR**3*ddegrdx


        dgmdx=half*SIX*bX*0.0042_r8_kind*dxdg*(two/sqrt(one+XGR**2) &
                                               -XGR**2/sqrt(one+XGR**2)/(one+XGR**2)) &   !!!!!
              /DEGR**2 &
             -half*SIX*0.0042_r8_kind**2*dxdg*ARSINHXGR(:,1)/DEGR**2 &
             -half*SIX*bX*0.0042_r8_kind*dxdg/y*dydx/DEGR**2 &
             -two*dGM/DEGR*ddegrdx

        dgmdx_g=half*SIX*0.0042_r8_kind**2*dxdg**2*(two/sqrt(one+XGR**2) &   !!!!!
                                               -XGR**2/sqrt(one+XGR**2)/(one+XGR**2)) &   !!!!!
                /DEGR**2 &
               -half*SIX*0.0042_r8_kind**2*dxdg*ARSINHXGR(:,1)/(gr(:vec_length,1)+depsMDA)/DEGR**2 &
               -half*SIX*0.0042_r8_kind**2*dxdg**2/y*dydx/DEGR**2 &
               -two*dGM_g/DEGR*ddegrdx

        df_drhodgamma(:vec_length,1)= df_drhodgamma(:vec_length,1) &
          -F4THRD*rho(:vec_length,1)**THIRD/TWO13*dxdg*(dgdx_g+dG_g+dxdg*dG_dx)/two

!        df_dgammadrho(:vec_length,1)= df_dgammadrho(:vec_length,1) &   !!!!!
!          -dgmdx*dxdr-half*(SIX*bX*0.0042_r8_kind*dxdgdr*(XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) ) &
!                           -two*0.0042_r8_kind*dxdgdr )/degr**2

        df_dgammadgamma(:vec_length,1)=df_dgammadgamma(:vec_length,1)+dgmdx_g*dxdg*half

        df_drhodrho(:vec_length,1)=df_drhodrho(:vec_length,1) &
          -F4THRD*THIRD/(rho(:vec_length,1)**F2THRD+depsMDA)*(G+XGR*dG)/TWO13 &
          -F4THRD*rho(:vec_length,1)**THIRD*dxdr*(dgdx+one*dG+XGR*dG_dx)/TWO13



                                                       !------
!  print*,'dgammadrho drhodgamma',sum(df_dgammadrho,1),sum(df_drhodgamma(:vec_length,1))



       endif

    else    ! spinpolarized case
       gr(:,1:2)=sqrt(gamma(:vec_length,1:2))

       XGRUP  = GR(:,1)  /(RHO(:vec_length,1) + DEPSmda)**F4THRD
       XGRDN  = GR(:,2)  /(RHO(:vec_length,2) + DEPSmda)**F4THRD
       arsinhxgr(:vec_length,1)=arsinh(xgrup)
       arsinhxgr(:vec_length,2)=arsinh(xgrdn)
       DEGRUP = ONE + SIX*0.0042_r8_kind*XGRUP*ARSINHXGR(:,1)
       DEGRDN = ONE + SIX*0.0042_r8_kind*XGRDN*ARSINHXGR(:,2)
       Gup = 0.0042_r8_kind*XGRUP**2/DEGRUP
       Gdn = 0.0042_r8_kind*XGRDN**2/DEGRDN
       EGR    = rho(1:vec_length,1)**F4THRD*Gup+rho(1:vec_length,2)**F4THRD*Gdn
       bXup = 0.0042_r8_kind*XGRUP
       bXdn = 0.0042_r8_kind*XGRDN

       dGup = ( six*bXup**2* ( XGRUP/sqrt(one+XGRUP**2)- ARSINHXGR(:,1) )- &
                two*bXup )/ DEGRUP**2
       dGdn = ( SIX*bXdn**2* ( XGRdn/sqrt(one+XGRdn**2)- ARSINHXGR(:,2) ) -&
                two*bXdn )/ DEGRdn**2
       dROup =  F4THRD*rho(:vec_length,1)**THIRD* (Gup+XGRUP*dGup)
       dROdn =  F4THRD*rho(:vec_length,2)**THIRD* (Gdn+XGRdn*dGdn)
       dGMup =  HALF*dGup/(GR(:,1)+depsMDA)
       dGMdn =  HALF*dGdn/(GR(:,2)+depsMDA)

       ! results
       fxc(:vec_length) = fxc(1:vec_length) - egr
       dfdrho(:vec_length,1) = -dROup+dfdrho(:vec_length,1)
       dfdrho(:vec_length,2) = -dROdn+dfdrho(:vec_length,2)
       dfdgrarho(:vec_length,1) = dGMup+dfdgrarho(:vec_length,1)
       dfdgrarho(:vec_length,2) = dGMdn+dfdgrarho(:vec_length,2)

       if(present(df_drhodrho)) then

       bX=bXup
       XGR=XGRup
       DEGR=DEGRUP
       dxdr=-F4THRD*gr(:,1)/(RHO(:vec_length,1) + DEPSmda)**f7THRD
       dxdg=one/(RHO(:vec_length,1)**F4THRD+DEPSmda)

       y=XGR+sqrt(one+XGR**2)
       dydx=(one+XGR/sqrt(one+XGR**2))
       ddegrdx=SIX*0.0042_r8_kind*ARSINHXGR(:,1)+SIX*bX/y*dydx
       dgdx=(two-XGR/degr*ddegrdx)*bx/degr
       dxdg=one/(RHO(:vec_length,1) + DEPSmda)**F4THRD

       dgdx_g=(two-XGR/degr*ddegrdx)*0.0042_r8_kind*dxdg/degr

       dG_dx=( SIX*TWO*bX*0.0042_r8_kind * ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) ) ) &
             / DEGR**2 &
            +( SIX*bX**2* (one/sqrt(one+XGR**2)-XGR**2/sqrt(one+XGR**2)/(one+XGR**2)  &
               - one/y*dydx))/DEGR**2 &
            - two*0.0042_r8_kind/DEGR**2 &
            - two*( SIX*bX**2* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) )- two*bX )/DEGR**3*ddegrdx

        dGM_g=half*(SIX*0.0042_r8_kind**2*dxdg**2*(XGR/sqrt(one+XGR**2)-ARSINHXGR(:,1) ) &
              - two*0.0042_r8_kind*dxdg/(gr(:vec_length,1)+depsMDA ))/DEGR**2

        dgmdx_g=half*SIX*0.0042_r8_kind**2*dxdg**2*(two/sqrt(one+XGR**2) &   !!!!!
                                               -XGR**2/sqrt(one+XGR**2)/(one+XGR**2)) &   !!!!!
                /DEGR**2 &
               -half*SIX*0.0042_r8_kind**2*dxdg*ARSINHXGR(:,1)/(gr(:vec_length,1)+depsMDA)/DEGR**2 &
               -half*SIX*0.0042_r8_kind**2*dxdg**2/y*dydx/DEGR**2 &
               -two*dGM_g/DEGR*ddegrdx

       dG_g= (SIX*bX*0.0042_r8_kind*dxdg* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,1) ) &
              - two*0.0042_r8_kind*dxdg)/DEGR**2


        !________________________________________________________
        df_drhodrho(:vec_length,1)=df_drhodrho(:vec_length,1) &
          -F4THRD*THIRD/(rho(:vec_length,1)**F2THRD+depsMDA)*(Gup+XGRup*dGup) &
          -F4THRD*rho(:vec_length,1)**THIRD*dxdr*(dgdx+dGup+XGR*dG_dx)
!!!       dROup =  F4THRD*rho(:vec_length,1)**THIRD* (Gup+XGRUP*dGup)
        df_drhodgamma(:vec_length,1)= df_drhodgamma(:vec_length,1) &
          -F4THRD*rho(:vec_length,1)**THIRD*dxdg*(dgdx_g+dG_g+dxdg*dG_dx)/2
        df_dgammadgamma(:vec_length,1)=df_dgammadgamma(:vec_length,1)+dgmdx_g*dxdg*half

       bX=bXdn
       XGR=XGRdn
       DEGR=DEGRdn
       dxdr=-F4THRD*gr(:,2)/(RHO(:vec_length,2) + DEPSmda)**f7THRD

       y=XGR+sqrt(one+XGR**2)
       dydx=(one+XGR/sqrt(one+XGR**2))
       ddegrdx=SIX*0.0042_r8_kind*ARSINHXGR(:,2)+SIX*bX/y*dydx
       dgdx=(two-XGR/degr*ddegrdx)*bx/degr
       dxdg=one/(RHO(:vec_length,2) + DEPSmda)**F4THRD

       dgdx_g=(two-XGR/degr*ddegrdx)*0.0042_r8_kind*dxdg/degr

       dG_dx=( SIX*TWO*bX*0.0042_r8_kind * ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,2) ) ) &
             / DEGR**2 &
            +( SIX*bX**2* (one/sqrt(one+XGR**2)-XGR**2/sqrt(one+XGR**2)/(one+XGR**2)  &
               - one/y*dydx))/DEGR**2 &
            - two*0.0042_r8_kind/DEGR**2 &
            - two*( SIX*bX**2* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,2) )- two*bX )/DEGR**3*ddegrdx

        dGM_g=half*(SIX*0.0042_r8_kind**2*dxdg**2*(XGR/sqrt(one+XGR**2)-ARSINHXGR(:,2) ) &
              - two*0.0042_r8_kind*dxdg/(gr(:vec_length,2)+depsMDA ))/DEGR**2

        dgmdx_g=half*SIX*0.0042_r8_kind**2*dxdg**2*(two/sqrt(one+XGR**2) &   !!!!!
                                               -XGR**2/sqrt(one+XGR**2)/(one+XGR**2)) &   !!!!!
                /DEGR**2 &
               -half*SIX*0.0042_r8_kind**2*dxdg*ARSINHXGR(:,2)/(gr(:vec_length,2)+depsMDA)/DEGR**2 &
               -half*SIX*0.0042_r8_kind**2*dxdg**2/y*dydx/DEGR**2 &
               -two*dGM_g/DEGR*ddegrdx


       dG_g= (SIX*bX*0.0042_r8_kind*dxdg* ( XGR/sqrt(one+XGR**2)- ARSINHXGR(:,2) ) &
              - two*0.0042_r8_kind*dxdg)/DEGR**2
        !____________________________________________________
        df_drhodrho(:vec_length,2)=df_drhodrho(:vec_length,2) &
          -F4THRD*THIRD/(rho(:vec_length,2)**F2THRD+depsMDA)*(Gdn+XGRdn*dGdn) &
          -F4THRD*rho(:vec_length,2)**THIRD*dxdr*(dgdx+dGdn+XGR*dG_dx)
        df_drhodgamma(:vec_length,2)= df_drhodgamma(:vec_length,2) &
          -F4THRD*rho(:vec_length,2)**THIRD*dxdg*(dgdx_g+dG_g+dxdg*dG_dx)/two
        df_dgammadgamma(:vec_length,2)=df_dgammadgamma(:vec_length,2)+dgmdx_g*dxdg*half

       endif
     end if
  end subroutine becke88_calcMDA

  !********************************************************************

  function arsinh(y)
    real(kind=r8_kind),intent(in),dimension(:) :: y
    real(kind=r8_kind),dimension(size(y)) :: arsinh
    arsinh=log(y+sqrt(y*y+1.0_r8_kind))
    return
  end function arsinh

end module becke_perdew_module

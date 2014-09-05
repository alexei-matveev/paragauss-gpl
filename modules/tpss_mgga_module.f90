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
module tpss_mgga_module
!----------------------------------------------------------------------------------------------------------------------------------+
!                                                                                                                                  |
!  Purpose: This module calculates the values of the TPSS                                                                          |
!           exchange-correlation energy functional based on                                                                        |
!                                                                                                                                  |
!  References: 1) J.Tao, J.P. Perdew, V.N. Staroverov and G.E. Scuseria                                                            |
!                 Phys. Rev. Lett. 91 (14), 2003, pages 146401-1 -- 146401-4                                                       |
!                                                                                                                                  |
!              2) J.Tao, J.P. Perdew, V.N. Staroverov and G.E. Scuseria                                                            |
!                 J. Chem .Phys. 120 (15), 2004, pages 6898 -- 6911                                                                |
!                                                                                                                                  |
!  Module called by: xc_hamiltonian                                                                                                |
!                                                                                                                                  |
!  Author: TS                                                                                                                      |
!  Date: August, 2009                                                                                                              |
!                                                                                                                                  |
!== Interrupt of public interface of module =======================================================================================+
!----------------------------------------------------------------------------------------------------------------------------------+
! Modifications                                                                                                                    |
!----------------------------------------------------------------------------------------------------------------------------------+
!                                                                                                                                  |
! Modification (Please copy before editing)                                                                                        |
! Author: ...                                                                                                                      |
! Date:   ...                                                                                                                      |
! Description: ...                                                                                                                 |
!                                                                                                                                  |
!------------ Modules used --------------------------------------------------------------------------------------------------------+
use type_module, RK=>r8_kind, IK=>i4_kind
implicit none
private
public :: tpss_mgga_X, tpss_mgga_C, pbe_c
!==================================================================================================================================+
! End of public interface of module                                                                                                |
!==================================================================================================================================+
!----------------------------------------------------------------------------------------------------------------------------------+
! local parameters                                                                                                                 |
!---------------------------------------------------------------+                                                                  |
real(RK),parameter  :: ZERO           = 0.00_RK               ,&!                                                                  |
                       quarter        = 0.25_RK               ,&!                                                                  |
                       third          = 1.0_RK/3.0_RK         ,&!                                                                  |
                       ninetwenty     = 9.0_RK/20.0_RK        ,&!                                                                  |
                       half           = 0.50_RK               ,&!                                                                  |
                       twothird       = 2.0_RK*third          ,&!                                                                  |
                       ONE            = 1.00_RK               ,&!                                                                  |
                       fourthird      = 4.0_RK*third          ,&!                                                                  |
                       fivethird      = 5.0_RK*third          ,&!                                                                  |
                       TWO            = 2.00_RK               ,&!                                                                  |
                       seventhird     = 7.0_RK*third          ,&!                                                                  |
                       eightthird     = 8.0_RK*third          ,&!                                                                  |
                       THREE          = 3.0_RK                ,&!                                                                  |
                       FOUR           = 4.0_RK                ,&!                                                                  |
                       ftthird        = 14.0_RK*third         ,&!                                                                  |
                       EIGHT          = 8.0_RK                ,&!                                                                  |
                       thld           = 1.0E-9_RK             ,&! Treshold value for evaluation cutoffs                            |
                       very_large     = 1.0E+30_RK            ,&! Treshold value for evaluation cutoffs                            |
                       c_LDA          = -0.738558766382022_RK   !                                                                  |
!---------------------------------------------------------------+                                                                  |
! subroutines                                                                                                                      |
!----------------------------------------------------------------------------------------------------------------------------------+
contains
  subroutine tpss_mgga_X(NGrid,ISpin,rho,gam,tau,F,dF_dr,dF_dg,dF_dt)
    implicit none
    !==============================================================================================================================+
    ! Public interface of subroutine                                                                                               |
    !====================================================+                                                                         |
    integer(IK),intent(in)                    :: NGrid ,&! Number of grid points                                                   |
                                                 ISpin   ! Number of spin directions                                               |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(in)        :: rho   ,&! Density. Array with dimensions (NGrid,ISpin)                            |
                                                 gam   ,&! Dot products of density gradients. (NGrid,2*ISpin-1)                    |
                                                 tau     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)         :: F       ! Kernel of xc-Functional                                                 |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(out)       :: dF_dr ,&! Functional derivatives                                                  |
                                                 dF_dg ,&!                                                                         |
                                                 dF_dt   !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of subroutine                                                                                        |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    real(RK),dimension(NGrid,2)               :: irho  ,&!                                                                         |
                                                 odFdr ,&!                                                                         |
                                                 itau  ,&!                                                                         |
                                                 odFdt   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,3)               :: igam  ,&!                                                                         |
                                                 odFdg   !                                                                         |
                                                         !-------------------------------------------------------------------------+

    if (ISpin .eq. 2) then
      call tpss(.TRUE., .FALSE., NGrid,2,rho,gam,tau,F,dF_dr,dF_dg,dF_dt)
    elseif(ISpin .eq. 1) then
      irho(:NGrid,1) = half           * rho(:NGrid,1)
      irho(:NGrid,2) = irho(:NGrid,1)
      igam(:NGrid,1) = quarter        * gam(:NGrid,1)
      igam(:NGrid,2) = igam(:NGrid,1)
      igam(:NGrid,3) = igam(:NGrid,1)
      itau(:NGrid,1) = half           * tau(:NGrid,1)
      itau(:NGrid,2) = itau(:NGrid,1)
      !
      call tpss(.TRUE., .FALSE.,NGrid,2,irho,igam,itau,F,odFdr,odFdg,odFdt)
      !
      dF_dr(:NGrid,1) = half    * sum(odFdr(:NGrid,:),dim=2)
      dF_dg(:NGrid,1) = quarter * sum(odFdg(:NGrid,:),dim=2)
      dF_dt(:NGrid,1) = half    * sum(odFdt(:NGrid,:),dim=2)
    end if
  end subroutine tpss_mgga_X

  subroutine tpss_mgga_C(NGrid,ISpin,rho,gam,tau,F,dF_dr,dF_dg,dF_dt)
    implicit none
    !==============================================================================================================================+
    ! Public interface of subroutine                                                                                               |
    !====================================================+                                                                         |
    integer(IK),intent(in)                    :: NGrid ,&! Number of grid points                                                   |
                                                 ISpin   ! Number of spin directions                                               |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(in)        :: rho   ,&! Density. Array with dimensions (NGrid,ISpin)                            |
                                                 gam   ,&! Dot products of density gradients. (NGrid,2*ISpin-1)                    |
                                                 tau     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)         :: F       ! Kernel of xc-Functional                                                 |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(out)       :: dF_dr ,&! Functional derivatives                                                  |
                                                 dF_dg ,&!                                                                         |
                                                 dF_dt   !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of subroutine                                                                                        |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    real(RK),dimension(NGrid,2)               :: irho  ,&!                                                                         |
                                                 odFdr ,&!                                                                         |
                                                 itau  ,&!                                                                         |
                                                 odFdt   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,3)               :: igam  ,&!                                                                         |
                                                 odFdg   !                                                                         |
                                                         !-------------------------------------------------------------------------+

    if (ISpin .eq. 2) then
      call tpss(.FALSE., .TRUE., NGrid,2,rho,gam,tau,F,dF_dr,dF_dg,dF_dt)
    elseif(ISpin .eq. 1) then
      irho(:NGrid,1) = half           * rho(:NGrid,1)
      irho(:NGrid,2) = irho(:NGrid,1)
      igam(:NGrid,1) = quarter        * gam(:NGrid,1)
      igam(:NGrid,2) = igam(:NGrid,1)
      igam(:NGrid,3) = igam(:NGrid,1)
      itau(:NGrid,1) = half           * tau(:NGrid,1)
      itau(:NGrid,2) = itau(:NGrid,1)
      !
      call tpss(.FALSE., .TRUE.,NGrid,2,irho,igam,itau,F,odFdr,odFdg,odFdt)
      !
      dF_dr(:NGrid,1) = half    * sum(odFdr(:NGrid,:),dim=2)
      dF_dg(:NGrid,1) = quarter * sum(odFdg(:NGrid,:),dim=2)
      dF_dt(:NGrid,1) = half    * sum(odFdt(:NGrid,:),dim=2)
    end if
  end subroutine tpss_mgga_C
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine tpss(exchon, corron, NGrid,ISpin,rho,gam,tau,F,dF_dr,dF_dg,dF_dt)
    implicit none
    !==============================================================================================================================+
    ! Public interface of subroutine                                                                                               |
    !====================================================+                                                                         |
    logical, intent(in)                       :: exchon,&!                                                                         |
                                                 corron  !                                                                         |
    integer(IK),intent(in)                    :: NGrid ,&! Number of grid points                                                   |
                                                 ISpin   ! Number of spin directions                                               |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(in)        :: rho   ,&! Density. Array with dimensions (NGrid,ISpin)                            |
                                                 gam   ,&! Dot products of density gradients. (NGrid,2*ISpin-1)                    |
                                                 tau     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)         :: F       ! Kernel of xc-Functional                                                 |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(out)       :: dF_dr ,&! Functional derivatives                                                  |
                                                 dF_dg ,&!                                                                         |
                                                 dF_dt   !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of subroutine                                                                                        |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !------------------------------------------------------+                                                                       |
    integer(IK)                               :: i1        !                                                                       |
                                                           !                                                                       |
    real(RK)                                  :: UKS_f     ! real value of spin (for scaling the density etc.)                     |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: p       ,&! reduced gradient                                                      |
                                                 z_x     ,&! tau^W / tau (for either alpha or beta)                                |
                                                 z_x_2   ,&! (z_x)^2                                                               |
                                                 q       ,&! substitute for reduced laplacian                                      |
                                                 x       ,&! inhomogeneity variable of exchange enhancement factor                 |
                                                 Fenx    ,&! MGGA enhancement factor                                               |
                                                 Xa      ,&! Xalpha factor of exchange integral kernel                             |
                                                 z_c     ,&! tau^W / tau (of total density)                                        |
                                                 z_c_2   ,&! (z_c)^2                                                               |
                                                 sumrho  ,&! total density                                                         |
                                                 sumgam  ,&! total squared gradient                                                |
                                                 sumtau  ,&! total kinetic energy density                                          |
                                                 ecPKZB  ,&! modified PKZB correlation energy density per particle                 |
                                                 ecTPSS  ,&! TPSS correlation energy density per particle                          |
                                                 av_5      ! Auxiliary variable                                                    |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: dzx_dr  ,&! d z_x / d rho                                                         |
                                                 dzx_dg  ,&! d z_x / d gam                                                         |
                                                 dzx_dt    ! d z_x / d tau                                                         |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: dp_dr   ,&! d p / d rho                                                           |
                                                 dp_dg     ! d p / d gam                                                           |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: dq_dp   ,&! d q / d p                                                             |
                                                 dq_dz     ! d q / d z_x                                                           |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: dx_dz2  ,&! d x / d z_x_2                                                         |
                                                 dx_dp   ,&! d x / d p                                                             |
                                                 dx_dq   ,&! d x / d q                                                             |
                                                 dx_dr   ,&! d x / d rho                                                           |
                                                 dx_dg   ,&! d x / d gam                                                           |
                                                 dx_dt     ! d x / d tau                                                           |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: dFenx_dx  ! d E_enhx / d x                                                        |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: dzc_dr  ,&! d z_c / d rho                                                         |
                                                 dzc_dg  ,&! d z_c / d gam                                                         |
                                                 dzc_dt    ! d z_c / d tau                                                         |
                                                           !                                                                       |
    real(RK),dimension(NGrid,ISpin)           :: de_dr     ! d ecPKZB / d rho                                                      |
                                                           !                                                                       |
    real(RK),dimension(NGrid,2*ISpin-1)       :: de_dg     ! d ecPKZB / d gam                                                      |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: de_dz2    ! d ecPKZB / d z_c_2                                                    |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: df_dz     ! d F_c / z                                                             |
                                                           !                                                                       |
    real(RK),dimension(NGrid)                 :: F_x     ,&! Exchange ...                                                          |
                                                 F_c       ! ... and Correlation part of the xc-Functional kernel                  |
    !------------------------------------------------------+-----------------------------------------------------------------------+
    ! local parameters                                                                                                             |
    !------------------------------------------------------+                                                                       |
    real(RK),parameter             :: kappa = 0.804_RK   ,&! Parameters                                                            |
                                      d_H   = 2.800_RK     !                                                                       |
    !------------------------------------------------------+                                                                       |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! initialize variables                                                                                                         |
    !------------------------------------------------------------------------------------------------------------------------------+
    F_x   = ZERO
    F_c   = ZERO
    F     = ZERO
    dF_dr = ZERO
    dF_dg = ZERO
    dF_dt = ZERO
    UKS_f = real(ISpin,RK)
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Exchange part: loop over polarizations                                                                                       |
    !------------------------------------------------------------------------------------------------------------------------------+
    do i1 = 1,ISpin
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call subroutine for the calculation of the variable z                                                                      |
      !----------------------------------------------------------------------------------------------------------------------------+
      call xc_z(rho(1:NGrid,i1),gam(1:NGrid,i1),tau(1:NGrid,i1),z_x,dzx_dr,dzx_dg,dzx_dt)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Square variable z                                                                                                          |
      !----------------------------------------------------------------------------------------------------------------------------+
      z_x_2 = z_x * z_x
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call subroutine for the calculation of the reduced gradient p                                                              |
      !----------------------------------------------------------------------------------------------------------------------------+
      call ex_p(ISpin,rho(1:NGrid,i1),gam(1:NGrid,i1),p,dp_dr,dp_dg)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call subroutine for the calculation of the reduced laplacian substitute q                                                  |
      !----------------------------------------------------------------------------------------------------------------------------+
      call ex_q(NGrid,p,z_x,q,dq_dp,dq_dz)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call subroutine for the calculation of the inhomogeneity variable x                                                        |
      !----------------------------------------------------------------------------------------------------------------------------+
      call ex_x(NGrid,z_x_2,p,q,x,dx_dz2,dx_dp,dx_dq)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate the exchange part of the integral kernel F_x                                                                     |
      !----------------------------------------------------------------------------------------------------------------------------+
      if (exchon) then
        where (rho(1:NGrid,i1) .ge. thld)
          !------------------------------------------------------------------------------------------------------------------------+
          ! Calculate the Xalpha factor of the exchange part                                                                       |
          !------------------------------------------------------------------------------------------------------------------------+
          Xa = c_LDA * rho(1:NGrid,i1) * (UKS_f * rho(1:NGrid,i1))**third
          !
          !------------------------------------------------------------------------------------------------------------------------+
          ! Calculate the exchange enhancement factor                                                                              |
          !------------------------------------------------------------------------------------------------------------------------+
          Fenx = ONE + kappa - kappa**2 / (kappa + x)
          !
          !------------------------------------------------------------------------------------------------------------------------+
          ! Sum up to total exchange part                                                                                          |
          !------------------------------------------------------------------------------------------------------------------------+
          F_x = F_x + Xa * Fenx
          !
          !------------------------------------------------------------------------------------------------------------------------+
          ! Calculate derivatives of exchange part                                                                                 |
          !------------------------------------------------------------------------------------------------------------------------+
          dFenx_dx = (kappa / (kappa + x))**2                                                                                      !
          !
          dx_dr = (dx_dq * dq_dz + TWO * dx_dz2 * z_x) * dzx_dr + (dx_dq * dq_dp + dx_dp) * dp_dr                                  !
          dx_dg = (dx_dq * dq_dz + TWO * dx_dz2 * z_x) * dzx_dg + (dx_dq * dq_dp + dx_dp) * dp_dg                                  !
          dx_dt = (dx_dq * dq_dz + TWO * dx_dz2 * z_x) * dzx_dt
          !
          dF_dr(:,i1) = (c_LDA * (fourthird * (UKS_f * rho(1:NGrid,i1))**third * Fenx) + Xa * dFenx_dx * dx_dr)
          dF_dg(:,i1) = Xa * dFenx_dx * dx_dg
          dF_dt(:,i1) = Xa * dFenx_dx * dx_dt
        end where
      end if ! exchon
          ! dx_dq * dq_dp + dx_dp                   ok (2x ref)
          ! dx_dq * dq_dz + TWO * dx_dz2 * z_x      ok (2x ref)
          ! dx_dr                                   ok (2x ref)
          ! dx_dg                                   ok (dp_dg and dzx_dg verfied numerically. REF: ~ (dx/dg) / |\nabla\rho|)
          ! dx_dt                                   ok (2x ref)
          ! dF_dr, dF_dg and dF_dt                  ok (== ref -- at the very end of the exchange part)
    end do
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Correlation part: total rho, gam and tau as well as spin resolved rho used                                                   |
    !------------------------------------------------------------------------------------------------------------------------------+
    if (corron) then
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate total rho, gam and tau                                                                                           |
      !----------------------------------------------------------------------------------------------------------------------------+
      sumrho = sum(rho(1:NGrid,1:ISpin),dim=2)
      sumgam = sum(gam(1:NGrid,1:(2*ISpin-1)),dim=2)
      sumtau = sum(tau(1:NGrid,1:ISpin),dim=2)
      !----------------------------------------------------------------------------------------------------------------------------+
      ! sumgam = (nabla sumrho)^2 = gam_aa + gam_bb + 2 * gam_ab                                                                   |
      !----------------------------------------------------------------------------------------------------------------------------+
      if (ISpin .eq. 2) then
        sumgam = sumgam + gam(1:NGrid,3)
      end if
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call subroutine for the calculation of the variable z_c                                                                    |
      !----------------------------------------------------------------------------------------------------------------------------+
      call xc_z(sumrho,sumgam,sumtau,z_c,dzc_dr,dzc_dg,dzc_dt)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Square variable z_c and calculate derivatives                                                                              |
      !----------------------------------------------------------------------------------------------------------------------------+
      where (sumrho .ge. thld .and. sumtau .ge. thld)
        z_c_2 = z_c * z_c
      elsewhere
        z_c_2 = ONE
      end where
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call subroutine for the calculation of the modified PBZB correlation energy density per particle                           |
      !----------------------------------------------------------------------------------------------------------------------------+
      call co_e(NGrid,ISpin,rho,sumrho,gam,z_c_2,ecPKZB,de_dr,de_dg,de_dz2)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate TPSS correlation energy density per particle                                                                     |
      !----------------------------------------------------------------------------------------------------------------------------+
      ecTPSS = ecPKZB * (ONE + d_H * ecPKZB * z_c_2 * z_c)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Sum up to total correlation part                                                                                           |
      !----------------------------------------------------------------------------------------------------------------------------+
      F_c = sumrho * ecTPSS
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate auxiliary variables for the derivatives of correlation part                                                      |
      !----------------------------------------------------------------------------------------------------------------------------+
      av_5  = sumrho * (ONE + d_H * TWO * ecPKZB * z_c_2 * z_c)                                                                    !
      df_dz = de_dz2 * TWO *  z_c * av_5 + sumrho * THREE * d_H * ecPKZB * ecPKZB * z_c_2                                          !
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate derivatives of correlation part                                                                                  |
      !----------------------------------------------------------------------------------------------------------------------------+
      dF_dr(:,1) = dF_dr(:,1) + df_dz * dzc_dr + de_dr(:,1) * av_5 + ecTPSS                                                        !
      !
      dF_dg(:,1) = dF_dg(:,1) + df_dz * dzc_dg + de_dg(:,1) * av_5                                                                 !
      !
      dF_dt(:,1) = dF_dt(:,1) + df_dz * dzc_dt                                                                                     !
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate additional derivatives in case of unrestricted calculation                                                       |
      !----------------------------------------------------------------------------------------------------------------------------+
      if (ISpin .eq. 2) then
        dF_dr(:,2) = dF_dr(:,2) + df_dz * dzc_dr       + de_dr(:,2) * av_5 + ecTPSS                                                !
        !
        dF_dg(:,2) = dF_dg(:,2) + df_dz * dzc_dg       + de_dg(:,2) * av_5                                                         !
        dF_dg(:,3) = dF_dg(:,3) + df_dz * dzc_dg * TWO + de_dg(:,3) * av_5                                                         !
        !
        dF_dt(:,2) = dF_dt(:,2) + df_dz * dzc_dt                                                                                   !
      end if
      !
    end if ! corron
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Sum up to total xc integral kernel F                                                                                         |
    !------------------------------------------------------------------------------------------------------------------------------+
    F = F_x + F_c
    !
  end subroutine tpss
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine xc_z(rho,gam,tau,z,dz_dr,dz_dg,dz_dt)
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    real(RK),dimension(:),intent(in)        :: rho     ,&! Density. Array with dimensions (NGrid,ISpin)                            |
                                               gam     ,&! Dot products of density gradients. (NGrid,2*ISpin-1)                    |
                                               tau       ! Kinetic energy density. Array with dimensions (NGrid,ISpin)             |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)       :: z       ,&! z - variable                                                            |
                                               dz_dr   ,&!                                                                         |
                                               dz_dg   ,&!                                                                         |
                                               dz_dt     !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    !----------------------------------------------------+                                                                         |
    ! local parameters                                                                                                             |
    !----------------------------------------------------+                                                                         |
    !----------------------------------------------------+                                                                         |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    where (rho .ge. thld .and. tau .ge. thld)
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate z_x  --  Note: z_x =< 1 for all real densities                                                                   |
      !----------------------------------------------------------------------------------------------------------------------------+
      dz_dg  = ONE / (EIGHT * rho * tau)
      z      = gam * dz_dg
      dz_dr  = - z / rho
      dz_dt  = - z / tau
    elsewhere
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Set limit values to 0                                                                                                      |
      !----------------------------------------------------------------------------------------------------------------------------+
      z      = ONE
      dz_dr  = ZERO
      dz_dg  = ZERO
      dz_dt  = ZERO
    end where
    !
  end subroutine xc_z
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine ex_p(ISpin,rho,gam,p,dp_dr,dp_dg)
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    integer(IK),intent(in)                  :: ISpin     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(in)        :: rho     ,&! Density. Array with dimensions (NGrid,ISpin)                            |
                                               gam       ! Dot products of density gradients. (NGrid,2*ISpin-1)                    |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)       :: p       ,&! Reduced density gradient                                                |
                                               dp_dr   ,&!                                                                         |
                                               dp_dg     !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    real(RK)                                :: c_dless   !                                                                         |
    !----------------------------------------------------+                                                                         |
    ! local parameters                                                                                                             |
    !----------------------------------------------------+                                                                         |
    real(RK),parameter :: c_rest = 38.2831200025092_RK ,&! 4 (3 pi^2)^(2/3) For the restricted case                                |
                          c_unrs = 60.7706649646079_RK   ! 4 (6 pi^2)^(2/3) For the unrestricted case                              |
    !----------------------------------------------------+                                                                         |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! set the parameter c_dless for either restricted or unrestricted case                                                         |
    !------------------------------------------------------------------------------------------------------------------------------+
    if (ISpin .eq. 1) then
      c_dless = c_rest
    else
      c_dless = c_unrs
    end if
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! calculate p and its derivatives  --  Note the use of dp_dg as auxiliary variable                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    where (rho .ge. thld)
      dp_dg = ONE / (c_dless * rho**eightthird)
      p     = gam * dp_dg
      dp_dr = - eightthird * p / rho
    elsewhere
      p     = ZERO
      dp_dr = ZERO
      dp_dg = ZERO
    end where
    !
  end subroutine ex_p
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine ex_q(NGrid,p,z,q,dq_dp,dq_dz)
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    integer(IK),intent(in)                    :: NGrid   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(in)          :: p     ,&!                                                                         |
                                                 z       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)         :: q     ,&!                                                                         |
                                                 dq_dp ,&!                                                                         |
                                                 dq_dz   !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    real(RK) ,dimension(NGrid)                :: av_1  ,&! auxiliary variable                                                      |
                                                 alpha ,&!                                                                         |
                                                 dq_da   !                                                                         |
    !----------------------------------------------------+                                                                         |
    ! local parameters                                                                                                             |
    !----------------------------------------------------+                                                                         |
    real(RK), parameter                   :: b = 0.4_RK  !                                                                         |
    !----------------------------------------------------+                                                                         |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Calculate q and its derivative                                                                                               |
    !------------------------------------------------------------------------------------------------------------------------------+
    where (z.ge. thld)
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate alpha, then the squareroot av_1 and finally q                                                                    |
      !----------------------------------------------------------------------------------------------------------------------------+
      alpha = abs(fivethird * p * (ONE / z - ONE))                                                                                 !
      av_1  = ONE / sqrt(ONE + b * alpha * (alpha - ONE))                                                                          !
      q     = ninetwenty * (alpha - ONE) * av_1 + twothird * p                                                                     !
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculate first derivative wrt alpha, then  d q / d p  and  d q / d z                                                      |
      !----------------------------------------------------------------------------------------------------------------------------+
      dq_da = ninetwenty * (ONE + half * b * (alpha - ONE)) * av_1**3                                                              !
      dq_dp = twothird + dq_da * fivethird * (ONE / z - ONE)                                                                       !
      dq_dz = - dq_da * fivethird * p / z**2                                                                                       !
    elsewhere
      !----------------------------------------------------------------------------------------------------------------------------+
      ! This should correspond to the case gam = 0                                                                                 |
      ! Use the limit from reference 2) and set dq_dz artificially to 0, because for  z = 0  also p may diverge                    |
      !----------------------------------------------------------------------------------------------------------------------------+
      q = ninetwenty / sqrt(b) + twothird * p
      dq_dz = ZERO
      dq_dp = twothird
      !
    end where
    !
  end subroutine ex_q
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine ex_x(NGrid,z2,p,q,x,dx_dz2,dx_dp,dx_dq)
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    integer(IK),intent(in)                    :: NGrid   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(in)          :: z2    ,&!                                                                         |
                                                 p     ,&!                                                                         |
                                                 q       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)         :: x     ,&!                                                                         |
                                                 dx_dz2,&!                                                                         |
                                                 dx_dp ,&!                                                                         |
                                                 dx_dq   !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !-----------------------------------------------------------+                                                                  |
    real(RK),dimension(NGrid)                :: av_2          ,&!                                                                  |
                                                av_3          ,&!                                                                  |
                                                term1         ,&!                                                                  |
                                                term2         ,&!                                                                  |
                                                term3         ,&!                                                                  |
                                                term4         ,&!                                                                  |
                                                term5         ,&!                                                                  |
                                                term6         ,&!                                                                  |
                                                denom           !                                                                  |
    !-----------------------------------------------------------+                                                                  |
    ! local parameters                                                                                                             |
    !-----------------------------------------------------------+                                                                  |
    real(RK), parameter  :: c10_81      = 10.0_RK/81.0_RK     ,&! 10/81                                                            |
                            c_c         = 1.59096_RK          ,&! c                                                                |
                            c146_2025   = 146.0_RK/2025.0_RK  ,&! 146 / 2025                                                       |
                            c73_405     = 73.0_RK/405.0_RK    ,&! 73 / 405                                                         |
                            threefifth2 = 9.0_RK/25.0_RK      ,&! (3 / 5)^2                                                        |
                            c_k         = ONE / 0.804_RK      ,&! 1/kappa                                                          |
                            c_e         = 1.23975804090960_RK ,&! sqrt(1.537)                                                      |
                            c_m         = 0.21951_RK            ! mu                                                               |
    !-----------------------------------------------------------+                                                                  |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Calculating different terms of x                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    av_2  = sqrt(half * (threefifth2 * z2 + p**2))
    av_3  = ONE / (ONE + z2)
    term1 = c10_81 + c_c * z2 * av_3**2
    term2 = c146_2025 * q
    term3 = - c73_405
    term4 = c_k * c10_81**2
    term5 = TWO * c_e * c10_81 * threefifth2
    term6 = c_e**2 * c_m * p
    denom = ONE / (ONE + c_e * p)
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Summarizing to x                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    x = (term1 * p + (term2 + term3 * av_2) * q + term5 * z2 + (term4 + term6) * p**2) * denom**2
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Calculation of derivatives                                                                                                   |
    !------------------------------------------------------------------------------------------------------------------------------+
    where(av_2 .ge. thld)
      dx_dz2 = (c_c * p * (ONE - TWO * z2 * av_3) * av_3**2 - quarter * c73_405 * q * threefifth2 / av_2 + term5) * denom**2       !
      dx_dp  = ((term1 + (half * term3 * q / av_2 + TWO * term4 + THREE * term6) * p) * denom - TWO * c_e * x) * denom             !
      dx_dq  = (TWO * term2 + term3 * av_2) * denom**2                                                                             !
    end where
    !
  end subroutine ex_x
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine co_c(NGrid,ISpin,rho,gam,c,dc_dr,dc_dg)
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    integer(IK),intent(in)                    :: NGrid ,&!                                                                         |
                                                 ISpin   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(in)        :: rho   ,&!                                                                         |
                                                 gam     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)         :: c       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(out)       :: dc_dr ,&!                                                                         |
                                                 dc_dg   !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    real(RK),dimension(NGrid)             :: sumrho    ,&!                                                                         |
                                             rezrho    ,&!                                                                         |
                                             av_4      ,&!                                                                         |
                                             zeta      ,&!                                                                         |
                                             zet2      ,&!                                                                         |
                                             xi2       ,&!                                                                         |
                                             cn        ,&!                                                                         |
                                             cd          !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2)           :: dzeta_dr  ,&!                                                                         |
                                             dxi2_dr     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,3)           :: dxi2_dg     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid)             :: dcn_dzet2 ,&!                                                                         |
                                             dcd_dzet  ,&!                                                                         |
                                             dcd_dxi2  ,&!                                                                         |
                                             dc_dzet   ,&!                                                                         |
                                             dc_dxi2     !                                                                         |
    !----------------------------------------------------+                                                                         |
    ! local parameters                                                                                                             |
    !----------------------------------------------------+                                                                         |
    real(RK), parameter   :: c0 = 0.53_RK              ,&!                                                                         |
                             c2 = 0.87_RK              ,&!                                                                         |
                             c4 = 0.50_RK              ,&!                                                                         |
                             c6 = 2.26_RK              ,&!                                                                         |
                             cx = 0.104484691940934_RK   ! cx = 1 / (3 * pi^2)^(2/3)                                               |
    !----------------------------------------------------+                                                                         |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    if (ISpin .eq. 1) then
      !----------------------------------------------------------------------------------------------------------------------------+
      ! In the restricted case, c reduces to the constant c0                                                                       |
      !----------------------------------------------------------------------------------------------------------------------------+
      c     = c0
      dc_dr = ZERO
      dc_dg = ZERO
    else
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Calculation of c for the unrestricted case                                                                                 |
      !----------------------------------------------------------------------------------------------------------------------------+
      sumrho = sum(rho(:,1:ISpin),dim=2)
      where (rho(:,1) .le. thld .or. rho(:,2) .le. thld .or. sumrho .le. thld)
        !--------------------------------------------------------------------------------------------------------------------------+
        ! Limits if at least one polarization of the density vanishes                                                              |
        !--------------------------------------------------------------------------------------------------------------------------+
        c          = ZERO
        dc_dr(:,1) = ZERO
        dc_dr(:,2) = ZERO
        dc_dg(:,1) = ZERO
        dc_dg(:,2) = ZERO
        dc_dg(:,3) = ZERO
      elsewhere
        rezrho = ONE / sumrho
        !--------------------------------------------------------------------------------------------------------------------------+
        ! Polarized case                                                                                                           |
        !--------------------------------------------------------------------------------------------------------------------------+
        av_4 = cx * rezrho**ftthird
        zeta = (rho(:,1) - rho(:,2)) * rezrho
        zet2 = zeta * zeta
        xi2  = (rho(:,2)**2 * gam(:,1) + rho(:,1)**2 * gam(:,2) - TWO * rho(:,1) * rho(:,2) * gam(:,3)) * av_4
        !
        !--------------------------------------------------------------------------------------------------------------------------+
        ! Calculation of Numerator and Denominator and c                                                                           |
        !--------------------------------------------------------------------------------------------------------------------------+
        cn   = c0 + (c2 + (c4 + c6 * zet2) * zet2) * zet2
        cd   = ONE / (ONE + half * xi2 * ((ONE + zeta)**(-fourthird) + (ONE - zeta)**(-fourthird)))
        c    = cn * cd**4
        !
        !--------------------------------------------------------------------------------------------------------------------------+
        ! Calculation of derivatives                                                                                               |
        !--------------------------------------------------------------------------------------------------------------------------+
        !
        dzeta_dr(:,1) =   TWO * rho(:,2) * rezrho**2
        dzeta_dr(:,2) = - TWO * rho(:,1) * rezrho**2
        dxi2_dr(:,1)  = (TWO * rho(:,1) * gam(:,2) - TWO * rho(:,2) * gam(:,3)) * av_4 - ftthird * xi2 * rezrho
        dxi2_dr(:,2)  = (TWO * rho(:,2) * gam(:,1) - TWO * rho(:,1) * gam(:,3)) * av_4 - ftthird * xi2 * rezrho
        dxi2_dg(:,1)  = rho(:,2)**2 * av_4
        dxi2_dg(:,2)  = rho(:,1)**2 * av_4
        dxi2_dg(:,3)  = - TWO * rho(:,1) * rho(:,2) * av_4
        !
        dcn_dzet2  = c2 + (TWO * c4 + THREE * c6 * zet2) * zet2                                                                    !
        dcd_dzet   = - twothird * xi2 * ((ONE + zeta)**(-seventhird) - (ONE - zeta)**(-seventhird))  ! wrong sign in reference !!!!!
        dcd_dxi2   = half * ((ONE + zeta)**(-fourthird) + (ONE - zeta)**(-fourthird))                                              !
        dc_dzet    = (TWO * zeta * dcn_dzet2 * cd**3 - FOUR * c * dcd_dzet) * cd                     ! refence corrected           !
        dc_dxi2    = - FOUR * c * dcd_dxi2 * cd                                                                                    !
        !
        dc_dr(:,1) = dc_dzet * dzeta_dr(:,1) + dc_dxi2 * dxi2_dr(:,1)                                                              !
        dc_dr(:,2) = dc_dzet * dzeta_dr(:,2) + dc_dxi2 * dxi2_dr(:,2)                                                              !
        dc_dg(:,1) = dc_dxi2 * dxi2_dg(:,1)                                                                                        !
        dc_dg(:,2) = dc_dxi2 * dxi2_dg(:,2)                                                                                        !
        dc_dg(:,3) = dc_dxi2 * dxi2_dg(:,3)                                                                                        !
        !
      end where
    end if
  end subroutine co_c
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine pbe_c(NGrid,ISpin,halfcalc,rho_in,sumrho,gam_in,e_cgga,degga_dr,degga_dg)
    use pw_ldac_module , only : pw_ldac
    use gga_response_module , only : gga_correlation
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    integer(IK),intent(in)               :: NGrid      ,&!                                                                         |
                                            ISpin        !                                                                         |
                                                         !                                                                         |
    logical,intent(in)                   :: halfcalc     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(in)     :: sumrho       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(in)   :: rho_in     ,&!                                                                         |
                                            gam_in       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)    :: e_cgga       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(out)  :: degga_dr   ,&!                                                                         |
                                            degga_dg     !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    real(RK),dimension(NGrid)            :: rezip_rho    !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,ISpin)      :: rho          !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2*ISpin-1)  :: gam          !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid)            :: F_clda     ,&!                                                                         |
                                            F_cgga       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,ISpin)      :: dF_clda_dr ,&!                                                                         |
                                            dF_cgga_dr   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2*ISpin-1)  :: dF_cgga_dg   !                                                                         |
    !----------------------------------------------------+                                                                         |
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Initialize values (pw_ldac and gga_correlation have  intent(inout)  for their outputs                                        |
    !------------------------------------------------------------------------------------------------------------------------------+
    F_clda     = ZERO
    F_cgga     = ZERO
    dF_clda_dr = ZERO
    dF_cgga_dr = ZERO
    dF_cgga_dg = ZERO
    !
    rho = rho_in ! + 1.0E-10_RK
    gam = gam_in ! + 1.0E-10_RK
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Call subroutines for the PBE-GGA correlation. First the LDA , then the GGA subroutine                                        |
    !------------------------------------------------------------------------------------------------------------------------------+
    call pw_ldac(NGrid,ISpin,rho,F_clda,dF_clda_dr)
    !
    if (halfcalc) then
      dF_clda_dr(:,2) = ZERO
    end if
    !
    call gga_correlation(3,ISpin,1,rho,gam,NGrid,F_cgga,F_clda,dF_cgga_dr,dF_cgga_dg,dF_clda_dr)
    !
    if (halfcalc) then
     dF_cgga_dr(:,2) = ZERO
     dF_cgga_dg(:,2:3) = ZERO
    end if
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Transform to GGA correlation energy density per particle  e_cgga  and calculate its derivatives                              |
    !------------------------------------------------------------------------------------------------------------------------------+
    where (sumrho .ge. thld)
      rezip_rho      = ONE / sumrho
      e_cgga         = (F_cgga + F_clda) * rezip_rho
      degga_dr(:,1)  = (dF_clda_dr(:,1) + dF_cgga_dr(:,1) - e_cgga) * rezip_rho
      degga_dr(:,2)  = (dF_clda_dr(:,2) + dF_cgga_dr(:,2) - e_cgga) * rezip_rho
      degga_dg(:,1)  = dF_cgga_dg(:,1) * rezip_rho
      degga_dg(:,2)  = dF_cgga_dg(:,2) * rezip_rho
      degga_dg(:,3)  = dF_cgga_dg(:,3) * rezip_rho
    elsewhere
      e_cgga         = ZERO
      degga_dr(:,1)  = ZERO
      degga_dr(:,2)  = ZERO
      degga_dg(:,1)  = ZERO
      degga_dg(:,2)  = ZERO
      degga_dg(:,3)  = ZERO
    end where
    !
  end subroutine pbe_c
  !================================================================================================================================+

  !================================================================================================================================+
  subroutine co_e(NGrid,ISpin,rho,sumrho,gam,zett2,epsPKZB,de_dr,de_dg,de_dz)
    implicit none
    !==============================================================================================================================+
    ! Public interface of module                                                                                                   |
    !====================================================+                                                                         |
    integer(IK),intent(in)                  :: NGrid   ,&!                                                                         |
                                               ISpin     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(in)        :: sumrho  ,&!                                                                         |
                                               zett2     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(in)      :: rho     ,&!                                                                         |
                                               gam       !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)       :: epsPKZB   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:,:),intent(out)     :: de_dr   ,&!                                                                         |
                                               de_dg     !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(:),intent(out)       :: de_dz     !                                                                         |
    !====================================================+                                                                         |
    ! End of public interface of module                                                                                            |
    !==============================================================================================================================+
    !------------------------------------------------------------------------------------------------------------------------------+
    ! local variables                                                                                                              |
    !----------------------------------------------------+                                                                         |
    integer(IK)                         :: is            ! index for spin                                                          |
                                                         !                                                                         |
    real(RK)                            :: RKS_f         !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid)           :: f_pbe_t     ,&! total pbe correlation energy density                                    |
                                           f_pbe_s     ,&!                                                                         |
                                           sumt        ,&!                                                                         |
                                           ct          ,&! C term                                                                  |
                                           rezip_rho     ! reziprocal density                                                      |
                                                         !                                                                         |
    real(RK),dimension(NGrid,Ispin)     :: df_pbe_t_dr   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2*Ispin-1) :: df_pbe_t_dg   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2)         :: h_rho_s     ,&!                                                                         |
                                           df_pbe_s_dr   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,3)         :: h_gam_s     ,&!                                                                         |
                                           df_pbe_s_dg   !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,Ispin)     :: epst          !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,ISpin,ISpin)     :: depst_dr!                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,ISpin,2*ISpin-1) :: depst_dg!                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,Ispin)     :: dc_dr         !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2*ISpin-1) :: dc_dg         !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,ISpin)     :: dsumt_dr      !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid,2*Ispin-1) :: dsumt_dg      !                                                                         |
                                                         !                                                                         |
    real(RK),dimension(NGrid)           :: term1       ,&!                                                                         |
                                           term2       ,&!                                                                         |
                                           term3         !                                                                         |
    !----------------------------------------------------+                                                                         |
    ! local parameters                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    ! executable lines                                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Initializing variables                                                                                                       |
    !------------------------------------------------------------------------------------------------------------------------------+
    depst_dr = ZERO
    depst_dg = ZERO
    !
    dsumt_dr = ZERO
    dsumt_dg = ZERO
    !
    epsPKZB  = ZERO
    de_dr    = ZERO
    de_dg    = ZERO
    de_dz    = ZERO
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Calculate reziprocal density                                                                                                 |
    !------------------------------------------------------------------------------------------------------------------------------+
    where (sumrho .ge. thld)
      rezip_rho = ONE / sumrho
    elsewhere
      rezip_rho = ZERO
    end where
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Factor: 0.5 for restricted calculation, 1 for unrestricted calculation                                                       |
    !------------------------------------------------------------------------------------------------------------------------------+
    RKS_f = ONE / real(3-ISpin,RK)
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Call PBE subroutine for PBE correlation energy of total density                                                              |
    !------------------------------------------------------------------------------------------------------------------------------+
    call pbe_c(NGrid,ISpin,.FALSE.,rho,sumrho,gam,f_pbe_t,df_pbe_t_dr,df_pbe_t_dg)
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Prepare variables for calculation PBE correlation energy of spin densities                                                   |
    !------------------------------------------------------------------------------------------------------------------------------+
    do is = 1,ISpin
      h_rho_s(1:NGrid,1) = rho(1:NGrid,is) * RKS_f
      h_rho_s(1:NGrid,2) = ZERO
      !
      h_gam_s(1:NGrid,1) = gam(1:NGrid,is) * RKS_f * RKS_f
      h_gam_s(1:NGrid,2) = ZERO
      h_gam_s(1:NGrid,3) = ZERO
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Call PBE subroutine for PBE correlation energy of spin densities                                                           |
      !----------------------------------------------------------------------------------------------------------------------------+
      call pbe_c(NGrid,2,.TRUE.,h_rho_s,h_rho_s(1:NGrid,1),h_gam_s,f_pbe_s,df_pbe_s_dr,df_pbe_s_dg)
      !
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Assign value to the "switch" term and calculate its derivatives                                                            |
      ! NB:         depst_dr(:,i,j) = d e^tilde_i / d rho_j                 depst_dg(:,i,j) = d e^tilde_i / d gam_j                |
      !----------------------------------------------------------------------------------------------------------------------------+
      where (f_pbe_t .ge. f_pbe_s)
        epst(1:NGrid,is)        = f_pbe_t(1:NGrid)
        depst_dr(1:NGrid,is,is) = df_pbe_t_dr(1:NGrid,1) ! is that correct?
        depst_dg(1:NGrid,is,is) = df_pbe_t_dg(1:NGrid,1)
      elsewhere
        epst(1:NGrid,is)        = f_pbe_s(1:NGrid)
        depst_dr(1:NGrid,is,is) = df_pbe_s_dr(1:NGrid,1)
        depst_dg(1:NGrid,is,is) = df_pbe_s_dg(1:NGrid,1)
      end where
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Add additional derivatives in case of unrestricted calculation                                                             |
      !----------------------------------------------------------------------------------------------------------------------------+
      if (ISpin == 2) then
        where (f_pbe_t .ge. f_pbe_s)
          depst_dr(1:NGrid,is,3-is) = df_pbe_t_dr(1:NGrid,3-is)
          depst_dg(1:NGrid,is,3-is) = df_pbe_t_dg(1:NGrid,3-is)
          depst_dg(1:NGrid,is,3)    = df_pbe_t_dg(1:NGrid,3)
        end where
      end if
      !
    end do
    !
    if (ISpin .eq. 1) then
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Assign value to the "switch" term of the beta polarization in case of a restricted calculation                             |
      !----------------------------------------------------------------------------------------------------------------------------+
      sumt = epst(1:NGrid,1)
    else
      !----------------------------------------------------------------------------------------------------------------------------+
      ! Explicitely sum up to sum-term of epsilon tilde                                                                            |
      !----------------------------------------------------------------------------------------------------------------------------+
      sumt = sum(rho(1:NGrid,1:ISpin) * epst(1:NGrid,1:ISpin),dim=2) * rezip_rho
    end if
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Call subroutine for the calculation of the C_term                                                                            |
    !------------------------------------------------------------------------------------------------------------------------------+
    call co_c(NGrid,ISpin,rho,gam,ct,dc_dr,dc_dg)                                                                                  !
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Some auxiliary variables                                                                                                     |
    !------------------------------------------------------------------------------------------------------------------------------+
    term1 = ONE + ct * zett2
    term2 = f_pbe_t * zett2
    term3 = ONE + ct
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Sum up to PKZB correlation energy density per particle                                                                       |
    !------------------------------------------------------------------------------------------------------------------------------+
    where (sumrho .ge. thld)
      epsPKZB   = f_pbe_t * term1 - term3 * zett2 * sumt
    elsewhere
      epsPKZB   = ZERO
    end where
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Calculation of first derivatives of sum-term                                                                                 |
    !------------------------------------------------------------------------------------------------------------------------------+
    if (ISpin .eq. 1) then
      where (sumrho .ge. thld)
        dsumt_dr(1:NGrid,1) = depst_dr(1:NGrid,1,1) * half
        dsumt_dg(1:NGrid,1) = depst_dg(1:NGrid,1,1) * quarter
      end where
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Derivatives in case of unrestricted calculation                                                                              |
    !------------------------------------------------------------------------------------------------------------------------------+
    else
      where (sumrho .ge. thld)
        dsumt_dr(1:NGrid,1) = (sum(rho(1:NGrid,1:ISpin) * depst_dr(1:NGrid,1:ISpin,1),dim=2) + epst(1:NGrid,1) - sumt) * rezip_rho
        dsumt_dr(1:NGrid,2) = (sum(rho(1:NGrid,1:ISpin) * depst_dr(1:NGrid,1:ISpin,2),dim=2) + epst(1:NGrid,2) - sumt) * rezip_rho
        dsumt_dg(1:NGrid,1) =  sum(rho(1:NGrid,1:ISpin) * depst_dg(1:NGrid,1:ISpin,1),dim=2)                     * rezip_rho
        dsumt_dg(1:NGrid,2) =  sum(rho(1:NGrid,1:ISpin) * depst_dg(1:NGrid,1:ISpin,2),dim=2)                     * rezip_rho
        dsumt_dg(1:NGrid,3) =  sum(rho(1:NGrid,1:ISpin) * depst_dg(1:NGrid,1:ISpin,3),dim=2)                     * rezip_rho
      end where
    end if
    !
    !------------------------------------------------------------------------------------------------------------------------------+
    ! Calculation of first derivatives                                                                                             |
    !------------------------------------------------------------------------------------------------------------------------------+
    where (sumrho .ge. thld)
      de_dr(:,1) = df_pbe_t_dr(:,1) * term1(:) + term2 * dc_dr(:,1) - zett2 * (dc_dr(:,1) * sumt + term3 * dsumt_dr(:,1))          !
      de_dg(:,1) = df_pbe_t_dg(:,1) * term1(:) + term2 * dc_dg(:,1) - zett2 * (dc_dg(:,1) * sumt + term3 * dsumt_dg(:,1))          !
    end where
    if (ISpin .eq. 2) then
      where (sumrho .ge. thld)
        de_dr(:,2) = df_pbe_t_dr(:,2) * term1(:) + term2 * dc_dr(:,2) - zett2 * (dc_dr(:,2) * sumt + term3 * dsumt_dr(:,2))        !
        de_dg(:,2) = df_pbe_t_dg(:,2) * term1(:) + term2 * dc_dg(:,2) - zett2 * (dc_dg(:,2) * sumt + term3 * dsumt_dg(:,2))        !
        de_dg(:,3) = df_pbe_t_dg(:,3) * term1(:) + term2 * dc_dg(:,3) - zett2 * (dc_dg(:,3) * sumt + term3 * dsumt_dg(:,3))        !
      end where
    end if
    where (sumrho .ge. thld)
      de_dz = f_pbe_t * ct - term3 * sumt                                                                                          !
    end where
    !
  end subroutine co_e
  !================================================================================================================================+
end module tpss_mgga_module

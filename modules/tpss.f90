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
module tpss
  !----------------------------------------------------------------------------+
  !                                                                            !
  !                                                                            !
  !                                                                            !
  !                                                                            !
  !                                                                            !
  !                                                                            !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  !  Purpose: ...                                                              !
  !                                                                            !
  !                                                                            !
  !  Module called by: ...                                                     !
  !                                                                            !
  !                                                                            !
  !  References: ...                                                           !
  !                                                                            !
  !                                                                            !
  !  Author: ...                                                               !
  !  Date: ...                                                                 !
  !                                                                            !
  !                                                                            !
  !----------------------------------------------------------------------------+
  ! Modifications                                                              !
  !----------------------------------------------------------------------------+
  !                                                                            !
  ! Modification (Please copy before editing)                                  !
  ! Author: ...                                                                !
  ! Date:   ...                                                                !
  ! Description: ...                                                           !
  !                                                                            !
  !----------------------------------------------------------------------------+
# include "def.h"
  use type_module,        only: IK => i4_kind, RK => r8_kind                   !
  !                                                                            !
  use ad1x7, only: ad &
                              , var                                            &
                              , fix                                            &
                              , val                                            &
                              , fst                                            &
                              , operator(+)                                    &
                              , operator(-)                                    &
                              , operator(*)                                    &
                              , operator(/)                                    &
                              , operator(**)                                   &
                              , exp                                            &
                              , log                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  public :: tpss_X                                                             &
          , tpss_C                                                             !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  integer(IK), parameter, public :: TPSS_key    = 1                            &
                                  , revTPSS_key = 2                            !
  !                                                                            !
  real(RK), parameter      :: DTol  = 1.0E-40_RK                               !
  real(RK), parameter      :: DTolr = 1.0E-30_RK                               !
  real(RK), parameter      :: pi    = 3.14159265358979323846_RK                !
  !                                                                            !
  real(RK), parameter      :: one   = 1.0_RK                                   !
  real(RK), parameter      :: n1d2  = 1.0_RK / 2.0_RK                          !
  real(RK), parameter      :: n2d3  = 2.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n4d3  = 4.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n14d3 = 14.0_RK / 3.0_RK                         !
  real(RK), parameter      :: n16d3 = 16.0_RK / 3.0_RK                         !
  real(RK), parameter      :: n56d3 = 56.0_RK / 3.0_RK                         !
  real(RK), parameter      :: cxi2  = 9.57078000062730_RK ! =(3 pi^2)^(2/3)    !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  contains                                                                     !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  !                           SUBROUTINES                                      !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine tpss_X(Flavor, NGrid, ISpin, rho, gr2, tau,                       &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! Flavor denotes type of functional:    1  -- TPSS                         !
    !                                       2  -- revTPSS                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: Flavor                                         !
    integer(IK), intent(in)  :: NGrid                                          !
    integer(IK), intent(in)  :: ISpin                                          !
    !                                                                          !
    real(RK),    intent(in)  :: rho(:,:)                                       !
    real(RK),    intent(in)  :: gr2(:,:)                                       !
    real(RK),    intent(in)  :: tau(:,:)                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),  intent(inout) :: F(:)                                           !
    real(RK),  intent(inout) :: dF_drho(:,:)                                   !
    real(RK),  intent(inout) :: dF_dgr2(:,:)                                   !
    real(RK),  intent(inout) :: dF_dtau(:,:)                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK)                 :: rho_tmp(NGrid, 2)                              !
    real(RK)                 :: gr2_tmp(NGrid, 3)                              !
    real(RK)                 :: tau_tmp(NGrid, 2)                              !
    real(RK)                 :: F_tmp(NGrid)                                   !
    real(RK)                 :: dF_drho_tmp(NGrid, 2)                          !
    real(RK)                 :: dF_dgr2_tmp(NGrid, 3)                          !
    real(RK)                 :: dF_dtau_tmp(NGrid, 2)                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    IF (ISpin .eq. 1) THEN                                                     !
      !------------------------------------------------------------------------+
      !                                                                        !
      rho_tmp(:, 1) = rho(:NGrid, 1) * 0.5_RK                                  !
      rho_tmp(:, 2) = rho(:NGrid, 1) * 0.5_RK                                  !
      gr2_tmp(:, 1) = gr2(:NGrid, 1) * 0.25_RK                                 !
      gr2_tmp(:, 2) = gr2(:NGrid, 1) * 0.25_RK                                 !
      gr2_tmp(:, 3) = gr2(:NGrid, 1) * 0.25_RK                                 !
      tau_tmp(:, 1) = tau(:NGrid, 1) * 0.5_RK                                  !
      tau_tmp(:, 2) = tau(:NGrid, 1) * 0.5_RK                                  !
      !                                                                        !
      call tpss_exchange(Flavor, NGrid, rho_tmp, gr2_tmp, tau_tmp,             &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp)  !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + dF_drho_tmp(:, 1)              !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + dF_dgr2_tmp(:, 1) * 0.5_RK     !
      dF_dtau(:NGrid, 1) = dF_dtau(:NGrid, 1) + dF_dtau_tmp(:, 1)              !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      call tpss_exchange( Flavor, NGrid                                        &
                        , rho(:NGrid,:2), gr2(:NGrid,:3), tau(:NGrid,:2),      &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp)  !
      !                                                                        !
      dF_drho(:NGrid,:2) = dF_drho(:NGrid,:2) + dF_drho_tmp                    !
      dF_dgr2(:NGrid,:3) = dF_dgr2(:NGrid,:3) + dF_dgr2_tmp                    !
      dF_dtau(:NGrid,:2) = dF_dtau(:NGrid,:2) + dF_dtau_tmp                    !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSE                                                                       !
      !------------------------------------------------------------------------+
      !                                                                        !
      WRITE(*,*) 'ISpin can only be either 1 or 2'                             !
      STOP                                                                     !
      !                                                                        !
      !------------------------------------------------------------------------+
    END IF                                                                     !
    !                                                                          !
    F(:NGrid)     = F(:NGrid) + F_tmp                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine tpss_X                                                        !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine tpss_C(PSet, NGrid, ISpin, rho, gr2, tau,                         &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! Flavor denotes type of functional:    1  -- TPSS                         !
    !                                       2  -- revTPSS                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: PSet                                           !
    integer(IK), intent(in)  :: NGrid                                          !
    integer(IK), intent(in)  :: ISpin                                          !
    !                                                                          !
    real(RK),    intent(in)  :: rho(:,:)                                       !
    real(RK),    intent(in)  :: gr2(:,:)                                       !
    real(RK),    intent(in)  :: tau(:,:)                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),  intent(inout) :: F(:)                                           !
    real(RK),  intent(inout) :: dF_drho(:,:)                                   !
    real(RK),  intent(inout) :: dF_dgr2(:,:)                                   !
    real(RK),  intent(inout) :: dF_dtau(:,:)                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+

    type(ad)                 :: r_a(NGrid) , r_b(NGrid)                        !
    type(ad)                 :: g_aa(NGrid), g_bb(NGrid), g_ab(NGrid)          !
    type(ad)                 :: t_a(NGrid) , t_b(NGrid)                        !
    type(ad)                 :: F_ad(NGrid)                                    !
    !                                                                          !
    integer(IK)              :: ind                                            !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    IF (ISpin .eq. 1) THEN                                                     !
      !------------------------------------------------------------------------+
      !                                                                        !
      r_a           = var( 1, rho(:NGrid, 1) * 0.50_RK )                       !
      r_b           = var( 2, rho(:NGrid, 1) * 0.50_RK )                       !
      g_aa          = var( 3, gr2(:NGrid, 1) * 0.25_RK )                       !
      g_bb          = var( 4, gr2(:NGrid, 1) * 0.25_RK )                       !
      g_ab          = var( 5, gr2(:NGrid, 1) * 0.25_RK )                       !
      t_a           = var( 6, tau(:NGrid, 1) * 0.50_RK )                       !
      t_b           = var( 7, tau(:NGrid, 1) * 0.50_RK )                       !
      !                                                                        !
      DO ind = 1, NGrid ! FIXME: make it elemental after verification          !
        F_ad(ind)   = tpss_correlation_ad( PSet                                &
                                         , r_a(ind) , r_b(ind)                 &
                                         , g_aa(ind), g_bb(ind), g_ab(ind)     &
                                         , t_a(ind) , t_b(ind)             )   !
      ENDDO                                                                    !
      !                                                                        !
      dF_drho(:, 1) = dF_drho(:, 1) + fst(1, F_ad)                             !
      dF_dgr2(:, 1) = dF_dgr2(:, 1) + fst(3, F_ad) * 0.5_RK                    !
      dF_dtau(:, 1) = dF_dtau(:, 1) + fst(6, F_ad)                             !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      r_a           = var( 1, rho(:NGrid, 1) )                                 !
      r_b           = var( 2, rho(:NGrid, 2) )                                 !
      g_aa          = var( 3, gr2(:NGrid, 1) )                                 !
      g_bb          = var( 4, gr2(:NGrid, 2) )                                 !
      g_ab          = var( 5, gr2(:NGrid, 3) )                                 !
      t_a           = var( 6, tau(:NGrid, 1) )                                 !
      t_b           = var( 7, tau(:NGrid, 2) )                                 !
      !                                                                        !
      DO ind = 1, NGrid                                                        !
        F_ad(ind)   = tpss_correlation_ad( PSet                                &
                                         , r_a(ind) , r_b(ind)                 &
                                         , g_aa(ind), g_bb(ind), g_ab(ind)     &
                                         , t_a(ind) , t_b(ind)             )   !
      ENDDO                                                                    !
      !                                                                        !
      dF_drho(:, 1) = dF_drho(:, 1) + fst(1, F_ad)                             !
      dF_drho(:, 2) = dF_drho(:, 2) + fst(2, F_ad)                             !
      dF_dgr2(:, 1) = dF_dgr2(:, 1) + fst(3, F_ad)                             !
      dF_dgr2(:, 2) = dF_dgr2(:, 2) + fst(4, F_ad)                             !
      dF_dgr2(:, 3) = dF_dgr2(:, 3) + fst(5, F_ad)                             !
      dF_dtau(:, 1) = dF_dtau(:, 1) + fst(6, F_ad)                             !
      dF_dtau(:, 2) = dF_dtau(:, 2) + fst(7, F_ad)                             !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSE                                                                       !
      !------------------------------------------------------------------------+
      !                                                                        !
      WRITE(*,*) 'ISpin can only be either 1 or 2'                             !
      STOP                                                                     !
      !                                                                        !
      !------------------------------------------------------------------------+
    END IF                                                                     !
    !                                                                          !
    F               = F + val(F_ad)                                            !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine tpss_C                                                        !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine tpss_exchange(Flavor, NGrid, rho, gr2, tau,                       &
                                            Fx, dFx_drho, dFx_dgr2, dFx_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    use pbe_gga_module,     only: eX_lsda                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: Flavor                                         !
    integer(IK), intent(in)  :: NGrid                                          !
    !                                                                          !
    real(RK),    intent(in)  :: rho(:, :) ! (NGrid, 2)
    real(RK),    intent(in)  :: gr2(:, :) ! (NGrid, 3)
    real(RK),    intent(in)  :: tau(:, :) ! (NGrid, 2)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: Fx(:) ! (NGrid)
    real(RK),    intent(out) :: dFx_drho(:, :) ! (NGrid, 2)
    real(RK),    intent(out) :: dFx_dgr2(:, :) ! (NGrid, 3)
    real(RK),    intent(out) :: dFx_dtau(:, :) ! (NGrid, 2)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)              :: s                                              !
    !                                                                          !
    real(RK)                 :: e(NGrid)                                       !
    real(RK)                 :: de_drho(NGrid)                                 !
    !                                                                          !
    real(RK)                 :: f(NGrid)                                       !
    real(RK)                 :: df_drho(NGrid)                                 !
    real(RK)                 :: df_dgr2(NGrid)                                 !
    real(RK)                 :: df_dtau(NGrid)                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    Fx       = 0.0_RK                                                          !
    dFx_drho = 0.0_RK                                                          !
    dFx_dgr2 = 0.0_RK                                                          !
    dFx_dtau = 0.0_RK                                                          !
    !                                                                          !
    spin: do s = 1, 2                                                          !
      !------------------------------------------------------------------------+
      !                                                                        !
      call eX_lsda(NGrid, rho(:NGrid, s), e, de_drho)                          !
      call exchef(Flavor, NGrid, rho(:NGrid, s), gr2(:NGrid, s), tau(:NGrid, s)&
                                               , f, df_drho, df_dgr2, df_dtau) !
      !                                                                        !
      Fx = Fx + e * f                                                          !
      !                                                                        !
      dFx_drho(1:NGrid, s) = de_drho * f                                       &
                           + e * df_drho                                       !
      dFx_dgr2(1:NGrid, s) = e * df_dgr2                                       !
      dFx_dtau(1:NGrid, s) = e * df_dtau                                       !
      !                                                                        !
      !------------------------------------------------------------------------+
    enddo spin                                                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine exchef(Flavor, NGrid, rho, gr2, tau,                            &
                                                 f, df_drho, df_dgr2, df_dtau) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: Flavor                                       !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(:) ! (NGrid)
      real(RK),    intent(in)  :: gr2(:) ! (NGrid)
      real(RK),    intent(in)  :: tau(:) ! (NGrid)
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: f(:) ! (NGrid)
      real(RK),    intent(out) :: df_drho(:) ! (NGrid)
      real(RK),    intent(out) :: df_dgr2(:) ! (NGrid)
      real(RK),    intent(out) :: df_dtau(:) ! (NGrid)
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                 :: x(NGrid)                                     !
      real(RK)                 :: dx_drho(NGrid)                               !
      real(RK)                 :: dx_dgr2(NGrid)                               !
      real(RK)                 :: dx_dtau(NGrid)                               !
      real(RK)                 :: df_dx(NGrid)                                 !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK), parameter      :: kappa = 0.804_RK                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      call mgga_inhomogeneity(Flavor, NGrid, rho, gr2, tau,                    &
                                                 x, dx_drho, dx_dgr2, dx_dtau) !
      !                                                                        !
      f     = 1.0_RK + kappa - kappa / (1.0_RK + x / kappa)                    !
      !                                                                        !
      df_dx = 1 / ((1.0_RK + x / kappa)**2)                                    !
      !                                                                        !
      df_drho = df_dx * dx_drho                                                !
      df_dgr2 = df_dx * dx_dgr2                                                !
      df_dtau = df_dx * dx_dtau                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine exchef                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine mgga_inhomogeneity(Flavor, NGrid, rho, gr2, tau,                &
                                                 x, dx_drho, dx_dgr2, dx_dtau) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: Flavor                                       !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      real(RK),    intent(in)  :: gr2(NGrid)                                   !
      real(RK),    intent(in)  :: tau(NGrid)                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: x(NGrid)                                     !
      real(RK),    intent(out) :: dx_drho(NGrid)                               !
      real(RK),    intent(out) :: dx_dgr2(NGrid)                               !
      real(RK),    intent(out) :: dx_dtau(NGrid)                               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                 :: z(NGrid)                                     !
      real(RK)                 :: dz_drho(NGrid)                               !
      real(RK)                 :: dz_dgr2(NGrid)                               !
      real(RK)                 :: dz_dtau(NGrid)                               !
      real(RK)                 :: p(NGrid)                                     !
      real(RK)                 :: dp_drho(NGrid)                               !
      real(RK)                 :: dp_dgr2(NGrid)                               !
      real(RK)                 :: a(NGrid)                                     !
      real(RK)                 :: da_drho(NGrid)                               !
      real(RK)                 :: da_dgr2(NGrid)                               !
      real(RK)                 :: da_dtau(NGrid)                               !
      real(RK)                 :: q(NGrid)                                     !
      real(RK)                 :: dq_da(NGrid)                                 !
      real(RK)                 :: dq_dp(NGrid)                                 !
      !                                                                        !
      real(RK)                 :: d1(NGrid)                                    !
      real(RK)                 :: t1(NGrid)                                    !
      real(RK)                 :: dt1_dz(NGrid)                                !
      real(RK)                 :: dt1_dp(NGrid)                                !
      real(RK)                 :: t2(NGrid)                                    !
      real(RK)                 :: dt2_dq(NGrid)                                !
      real(RK)                 :: s3(NGrid)                                    !
      real(RK)                 :: t3(NGrid)                                    !
      real(RK)                 :: dt3_dz(NGrid)                                !
      real(RK)                 :: dt3_dp(NGrid)                                !
      real(RK)                 :: dt3_dq(NGrid)                                !
      real(RK)                 :: t4(NGrid)                                    !
      real(RK)                 :: dt4_dp(NGrid)                                !
      real(RK)                 :: t5(NGrid)                                    !
      real(RK)                 :: dt5_dz(NGrid)                                !
      real(RK)                 :: t6(NGrid)                                    !
      real(RK)                 :: dt6_dp(NGrid)                                !
      real(RK)                 :: t7(NGrid)                                    !
      real(RK)                 :: dt7_dp(NGrid)                                !
      !                                                                        !
      real(RK)                 :: dx_dp(NGrid)                                 !
      real(RK)                 :: dx_dz(NGrid)                                 !
      real(RK)                 :: dx_dq(NGrid)                                 !
      !                                                                        !
      real(RK)                 :: c_par                                        !
      real(RK)                 :: e_par                                        !
      real(RK)                 :: mupar                                        !
      real(RK), parameter      :: kappa = 0.804_RK                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      call get_pars_tpss_exchange(Flavor, c_par, e_par, mupar)                 !
      !                                                                        !
      call get_z_var(NGrid, rho, gr2, tau, z, dz_drho, dz_dgr2, dz_dtau)       !
      call get_p_var(NGrid, rho, gr2, p, dp_drho, dp_dgr2)                     !
      call get_a_var(NGrid, rho, gr2, tau, a, da_drho, da_dgr2, da_dtau)       !
      call get_q_var(NGrid, a, p, q, dq_da, dq_dp)                             !
      !                                                                        !
      d1 = 1.0_RK / (1.0_RK + z**2)                                            !
      !                                                                        !
      if (Flavor == TPSS_key) then                                             !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        t1     = (10.0_RK / 81.0_RK + c_par * z**2 * d1**2) * p                !
        !                                                                      !
        dt1_dz = (2.0_RK * c_par * z    * d1**2                                &
                - 4.0_RK * c_par * z**3 * d1**3) * p                           !
        !                                                                      !
        dt1_dp = (10.0_RK / 81.0_RK + c_par * z**2 * d1**2)                    !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      elseif (Flavor == revTPSS_key) then                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        t1     = (10.0_RK / 81.0_RK + c_par * z**3 * d1**2) * p                !
        !                                                                      !
        dt1_dz = (3.0_RK * c_par * z**2 * d1**2                                &
                - 4.0_RK * c_par * z**4 * d1**3) * p                           !
        !                                                                      !
        dt1_dp = (10.0_RK / 81.0_RK + c_par * z**3 * d1**2)                    !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      endif                                                                    !
      !                                                                        !
      t2      = (146.0_RK / 2025.0_RK) * q**2                                  !
      !                                                                        !
      dt2_dq  = (146.0_RK / 2025.0_RK) * 2.0_RK * q                            !
      !                                                                        !
      s3      = sqrt(0.5_RK * (3.0_RK * z / 5.0_RK)**2 + 0.5_RK * p**2)        !
      !                                                                        !
      t3      = - (73.0_RK / 405.0_RK) * s3 * q                                !
      !                                                                        !
      where (s3 .gt. DTol)                                                     !
        !                                                                      !
        dt3_dz  = - q * (73.0_RK / (405.0_RK * s3)) * z * 4.5_RK / 25.0_RK     !
        !                                                                      !
        dt3_dp  = - q * (0.5_RK * 73.0_RK / (405.0_RK * s3)) * p               !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        dt3_dz = 0.0_RK                                                        !
        !                                                                      !
        dt3_dp = 0.0_RK                                                        !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      dt3_dq  = - (73.0_RK / 405.0_RK) * s3                                    !
      !                                                                        !
      t4      = (1.0_RK / kappa) * (10.0_RK / 81.0_RK)**2 * p**2               !
      !                                                                        !
      dt4_dp  = (1.0_RK / kappa) * (10.0_RK / 81.0_RK)**2 * 2.0_RK * p         !
      !                                                                        !
      t5      = 2.0_RK * sqrt(e_par) * (10.0_RK / 81.0_RK)                     &
              * (3.0_RK * z / 5.0_RK)**2                                       !
      !                                                                        !
      dt5_dz  = 2.0_RK * sqrt(e_par) * (10.0_RK / 81.0_RK)                     &
              * (3.0_RK / 5.0_RK)**2 * 2.0_RK * z                              !
      !                                                                        !
      t6      = e_par * mupar * p**3                                           !
      !                                                                        !
      dt6_dp  = e_par * mupar * 3.0_RK * p**2                                  !
      !                                                                        !
      t7      = 1.0_RK + sqrt(e_par) * p                                       !
      !                                                                        !
      dt7_dp  = sqrt(e_par)                                                    !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      x       = (t1 + t2 + t3 + t4 + t5 + t6) / (t7**2)                        !
      !                                                                        !
      dx_dz   = (dt1_dz + dt3_dz + dt5_dz) / (t7**2)                           !
      !                                                                        !
      dx_dp   = (dt1_dp + dt3_dp + dt4_dp + dt6_dp) / (t7**2)                  &
              - 2.0_RK * x * dt7_dp / t7                                       !
      !                                                                        !
      dx_dq   = (dt2_dq + dt3_dq) / (t7**2)                                    !
      !                                                                        !
      dx_drho = dx_dz * dz_drho                                                &
              + dx_dq * dq_da * da_drho                                        &
              + (dx_dp + dx_dq * dq_dp) * dp_drho                              !
      !                                                                        !
      dx_dgr2 = dx_dz * dz_dgr2                                                &
              + dx_dq * dq_da * da_dgr2                                        &
              + (dx_dp + dx_dq * dq_dp) * dp_dgr2                              !
      !                                                                        !
      dx_dtau = dx_dz * dz_dtau                                                &
              + dx_dq * dq_da * da_dtau                                        !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine mgga_inhomogeneity                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine get_pars_tpss_exchange(Flavor, c_par, e_par, mupar)             !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: Flavor                                       !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: c_par                                        !
      real(RK),    intent(out) :: e_par                                        !
      real(RK),    intent(out) :: mupar                                        !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (Flavor)                                                     !
        case(TPSS_key)  ! TPSS                                                 !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          c_par = 1.59096_RK                                                   !
          e_par = 1.537_RK                                                     !
          mupar = 0.21951_RK                                                   !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(revTPSS_key)  ! revTPSS                                           !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          c_par = 2.35204_RK                                                   !
          e_par = 2.1677_RK                                                    !
          mupar = 0.14_RK                                                      !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case default                                                           !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          c_par = 0.0_RK                                                       !
          e_par = 0.0_RK                                                       !
          mupar = 0.0_RK                                                       !
          write(*,*) 'unknown key in TPSS -> stop'; stop                       !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_tpss_exchange                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine get_z_var(NGrid, rho, gr2, tau, z, dz_drho, dz_dgr2, dz_dtau)   !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      real(RK),    intent(in)  :: gr2(NGrid)                                   !
      real(RK),    intent(in)  :: tau(NGrid)                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: z(NGrid)                                     !
      real(RK),    intent(out) :: dz_drho(NGrid)                               !
      real(RK),    intent(out) :: dz_dgr2(NGrid)                               !
      real(RK),    intent(out) :: dz_dtau(NGrid)                               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where(rho .ge. DTol .and. tau .ge. DTol)                                 !
        !                                                                      !
        z       = gr2 / (8.0_RK * rho * tau)                                   !
        !                                                                      !
        dz_drho = - gr2 / (8.0_RK * rho**2 * tau)                              !
        !                                                                      !
        dz_dgr2 = 1.0_RK / (8.0_RK * rho * tau)                                !
        !                                                                      !
        dz_dtau = - gr2 / (8.0_RK * rho * tau**2)                              !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        z       = 0.0_RK                                                       !
        !                                                                      !
        dz_drho = 0.0_RK                                                       !
        !                                                                      !
        dz_dgr2 = 0.0_RK                                                       !
        !                                                                      !
        dz_dtau = 0.0_RK                                                       !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_z_var                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine get_p_var(NGrid, rho, gr2, p, dp_drho, dp_dgr2)                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      real(RK),    intent(in)  :: gr2(NGrid)                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: p(NGrid)                                     !
      real(RK),    intent(out) :: dp_drho(NGrid)                               !
      real(RK),    intent(out) :: dp_dgr2(NGrid)                               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where(rho .ge. DTol)                                                     !
        !                                                                      !
        p       = gr2                                                          &
         / (4.0_RK * (6.0_RK * pi**2)**(2.0_RK/3.0_RK) * rho**(8.0_RK / 3.0_RK))
        !                                                                      !
        dp_drho = - (8.0_RK / 3.0_RK) * p / rho                                !
        !                                                                      !
        dp_dgr2 = 1.0_RK                                                       &
         / (4.0_RK * (6.0_RK * pi**2)**(2.0_RK/3.0_RK) * rho**(8.0_RK / 3.0_RK))
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        p       = 0.0_RK                                                       !
        !                                                                      !
        dp_drho = 0.0_RK                                                       !
        !                                                                      !
        dp_dgr2 = 0.0_RK                                                       !
        !                                                                      !
      endwhere                                                                 !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_p_var                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine get_q_var(NGrid, a, p, q, dq_da, dq_dp)                         !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: a(NGrid)                                     !
      real(RK),    intent(in)  :: p(NGrid)                                     !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: q(NGrid)                                     !
      real(RK),    intent(out) :: dq_da(NGrid)                                 !
      real(RK),    intent(out) :: dq_dp(NGrid)                                 !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK), parameter      :: b_par = 0.40_RK                              !
      !                                                                        !
      real(RK)                 :: denom(NGrid)                                 !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      denom = 1.0_RK / sqrt(1.0_RK + b_par * a * (a - 1.0_RK))                 !
      !                                                                        !
      q     = 0.45_RK * (a - 1.0_RK) * denom + (2.0_RK / 3.0_RK) * p           !
      !                                                                        !
      dq_da = 0.45_RK * (denom                                                 &
            - 0.50_RK * b_par * (a - 1.0_RK) * (2.0_RK * a - 1.0_RK) * denom**3)
      !                                                                        !
      dq_dp = 2.0_RK / 3.0_RK                                                  !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_q_var                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine get_a_var(NGrid, rho, gr2, tau, a, da_drho, da_dgr2, da_dtau)   !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      real(RK),    intent(in)  :: gr2(NGrid)                                   !
      real(RK),    intent(in)  :: tau(NGrid)                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: a(NGrid)                                     !
      real(RK),    intent(out) :: da_drho(NGrid)                               !
      real(RK),    intent(out) :: da_dgr2(NGrid)                               !
      real(RK),    intent(out) :: da_dtau(NGrid)                               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK), parameter      :: CT = 4.55779987234560_RK ! 3/10(6 pi^2)^(2/3)!
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (rho .gt. DTol)                                                    !
        !                                                                      !
        a       = (tau - gr2 / (8.0_RK * rho)) / (CT * rho**(5.0_RK/3.0_RK))   !
        !                                                                      !
        da_drho = (8.0_RK/3.0_RK) * gr2 / (8.0_RK * CT * rho**(11.0_RK/3.0_RK))&
                - (5.0_RK/3.0_RK) * tau / (CT * rho**(8.0_RK/3.0_RK))          !
        !                                                                      !
        da_dgr2 = - 1.0_RK / (8.0_RK * CT * rho **(8.0_RK/3.0_RK))             !
        !                                                                      !
        da_dtau = 1.0_RK / (CT * rho **(5.0_RK/3.0_RK))                        !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        a       = 0.0_RK                                                       !
        !                                                                      !
        da_drho = 0.0_RK                                                       !
        !                                                                      !
        da_dgr2 = 0.0_RK                                                       !
        !                                                                      !
        da_dtau = 0.0_RK                                                       !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_a_var                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine tpss_exchange                                                 !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  elemental function tpss_correlation_ad( PSet, r_a, r_b                       &
                                        , gaa, gbb, gab , t_a, t_b ) result(F) !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: PSet                                           !
    !                                                                          !
    type(ad),    intent(in)  :: r_a, r_b                                       !
    type(ad),    intent(in)  :: gaa, gbb, gab                                  !
    type(ad),    intent(in)  :: t_a, t_b                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: F                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: qe                                             !
    type(ad)                 :: qz                                             !
    type(ad)                 :: rt                                             !
    type(ad)                 :: tt                                             !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: d_par = 2.80_RK                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    rt = r_a + r_b                                                             !
    tt = t_a + t_b                                                             !
    !                                                                          !
    IF ( val(rt) < DTolr ) THEN                                                !
      !                                                                        !
      ! SCREENING: NO DENSITY, NO CORRELATION ENERGY                           !
      !                                                                        !
      F  = fix(0.0_RK)                                                         !
      !                                                                        !
    ELSE                                                                       !
      !                                                                        !
      qz = tauW_d_tau( rt, gaa, gbb, gab, tt )                                 !
      !                                                                        !
      qe = rev_pkzb_correlation( PSet, r_a, r_b, rt, gaa, gbb, gab, qz )       !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      ! ADDING EVERYTHING TO TPSS CORRELATION ENERGY DENSITY                   !
      !                                                                        !
      F  = rt * qe * ( one + d_par * qe * qz**3 )                              !
      !                                                                        !
    ENDIF                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function rev_pkzb_correlation( PSet, ra, rb, rt, ga, gb, gc, z ) &
                                                                     result(e) !
      !------------------------------------------------------------------------+
      !                                                                        !
      use pbe_gga_module,     only: pbe_correlation_ad                         &
                                  , PBE_key, PBE_revTPSS_key                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in) :: PSet                                          !
      !                                                                        !
      type(ad),    intent(in) :: ra, rb, rt                                    !
      type(ad),    intent(in) :: ga, gb, gc ! = g_aa, g_bb, g_ab               !
      type(ad),    intent(in) :: z                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: e                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK)             :: Key                                           !
      !                                                                        !
      type(ad)                :: qc                                            !
      type(ad)                :: qpo, qpa, qpb                                 !
      type(ad)                :: f0                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      f0    = fix( 0.0_RK )                                                    !
      !                                                                        !
      ! GET INTERMEDIATE QUANTITY C(zeta,xi)                                   !
      !                                                                        !
      qc    = c_term( PSet, ra, rb, rt, ga, gb, gc )                           !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      ! PBE TERMS: total (qpo), polarized a (qpa) and polarized b (qpb)        !
      !                                                                        !
      ! Set what to do in PBE                                                  !
      IF ( PSet == TPSS_key)    Key = PBE_key                                  !
      IF ( PSet == revTPSS_key) Key = PBE_revTPSS_key                          !
      !                                                                        !
      !                                                                        !
      ! PBE energy density per particle (total)                                !
      qpo   = pbe_correlation_ad( Key, ra, rb, ga, gb, gc )                    !
      !                                                                        !
      !                                                                        !
      ! PBE energy density per particle (alpha)                                !
      qpa   = pbe_correlation_ad( Key, ra, f0, ga, f0, f0 )                    !
      !                                                                        !
      !                                                                        !
      ! PBE energy density per particle (beta)                                 !
      qpb   = pbe_correlation_ad( Key, f0, rb, f0, gb, f0 )                    !
      !                                                                        !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      ! OPPOSITE SPIN PART                                                     !
      !                                                                        !
      e     = qpo * (one + qc * z**2)                                          !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      ! SAME SPIN PART ALPHA                                                   !
      !                                                                        !
      IF( val(qpa) > val(qpo) ) THEN                                           !
        e   = e - (one + qc) * z**2 * (ra / rt) * qpa                          !
      ELSE                                                                     !
        e   = e - (one + qc) * z**2 * (ra / rt) * qpo                          !
      END IF                                                                   !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      ! SAME SPIN PART BETA                                                    !
      !                                                                        !
      IF( val(qpb) > val(qpo) ) THEN                                           !
        e   = e - (one + qc) * z**2 * (rb / rt) * qpb                          !
      ELSE                                                                     !
        e   = e - (one + qc) * z**2 * (rb / rt) * qpo                          !
      END IF                                                                   !
      !                                                                        !
      e = e * 1.0_RK                                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function rev_pkzb_correlation                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function tauW_d_tau( rt, ga, gb, gc, tt ) result(z)              !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: rt                                            !
      type(ad),    intent(in) :: ga, gb, gc ! = g_aa, g_bb, g_ab               !
      type(ad),    intent(in) :: tt                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: z                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      IF( val(tt) > DTol ) THEN                                                !
        !                                                                      !
        z = ( ga + gb + 2.0_RK * gc ) / ( 8.0_RK * rt * tt )                   !
        !                                                                      !
      ELSE                                                                     !
        !                                                                      !
        z = fix(0.0_RK)                                                        !
        !                                                                      !
      END IF                                                                   !
      !                                                                        !
      z = z * 1.0_RK
      !                                                                        !
      !------------------------------------------------------------------------+
    end function tauW_d_tau                                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function c_term( PSet, ra, rb, rt, ga, gb, gc ) result(c)        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in) :: PSet                                          !
      !                                                                        !
      type(ad),    intent(in) :: ra, rb, rt                                    !
      type(ad),    intent(in) :: ga, gb, gc ! = g_aa, g_bb, g_ab               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: c                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                :: par_c1, par_c2, par_c3, par_c4                !
      type(ad)                :: qx2                                           !
      type(ad)                :: qy2                                           !
      type(ad)                :: cn, cd                                        !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      call get_C_term_parameters( PSet, par_c1, par_c2, par_c3, par_c4)        !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      ! Initially :  x2 = [ |grad(z)| / (2 (3 pi^2 rho)^(2/3)) ]^2             !
      ! Now:         x2 = ( a^2 * ga + b^2 * gb - 2 * a * b * gc ) / const     !
      ! with factor  rt^(14/3)  put outside to reduce divergencies             !
      !                                                                        !
      qx2 = spin_pol_gradient( ra, rb, ga, gb, gc )                            !
      !                                                                        !
      !                                                                        !
      ! y = (a - b) / (a + b)                                                  !
      qy2 = relative_spin_pol_squared( ra, rb, rt )                            !
      !                                                                        !
      ! cn = C(qy,0) * (a*b)^(16/3) * (a+b)^(56/3)                             !
      cn  = (par_c1 + (par_c2 + (par_c3 + par_c4 * qy2) * qy2) * qy2)          &
            * ( ra * rb )**n16d3 * rt**n56d3                                   !
      !                                                                        !
      cd  = (( ra * rb )**n4d3 * rt**n14d3                                     &
                 + qx2 * n1d2 * (rt * n1d2)**n4d3 * (ra**n4d3 + rb**n4d3 ))**4 &
                 + DTol                                                        !
      !                                                                        !
      !                                                                        !
      c = cn / cd ! SCREENING DONE ABOVE                                       !
      !                                                                        !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function c_term                                                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function relative_spin_pol_squared( ra, rb, rt ) result(y2)      !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: ra, rb, rt                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: y2                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      y2 = (( ra - rb ) / rt )**2 ! SCREENING DONE OUTSIDE                     !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function relative_spin_pol_squared                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function spin_pol_gradient( ra, rb, ga, gb, gc ) result(x2)      !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: ra, rb                                        !
      type(ad),    intent(in) :: ga, gb, gc ! = g_aa, g_bb, g_ab               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: x2                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      x2 = ( ra**2 * gb + rb**2 * ga - 2.0_RK * ra * rb * gc ) / cxi2          !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function spin_pol_gradient                                             !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    pure subroutine get_c_term_parameters( PSet, c1, c2, c3, c4 )              !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: c1, c2, c3, c4                               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case ( PSet )                                                     !
        case(TPSS_key)  ! TPSS                                                 !
          !                                                                    !
          c1 = 0.5300_RK                                                       !
          c2 = 0.8700_RK                                                       !
          c3 = 0.5000_RK                                                       !
          c4 = 2.2600_RK                                                       !
          !                                                                    !
        case(revTPSS_key)  ! revTPSS                                           !
          !                                                                    !
          c1 = 0.5900_RK                                                       !
          c2 = 0.9269_RK                                                       !
          c3 = 0.6225_RK                                                       !
          c4 = 2.1540_RK                                                       !
          !                                                                    !
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_c_term_parameters                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end function tpss_correlation_ad                                             !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine tpss_correlation(Flavor, NGrid, rho, gr2, tau,                    &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: Flavor                                         !
    integer(IK), intent(in)  :: NGrid                                          !
    !                                                                          !
    real(RK),    intent(in)  :: rho(NGrid, 2)                                  !
    real(RK),    intent(in)  :: gr2(NGrid, 3)                                  !
    real(RK),    intent(in)  :: tau(NGrid, 2)                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: F(NGrid)                                       !
    real(RK),    intent(out) :: dF_drho(NGrid, 2)                              !
    real(RK),    intent(out) :: dF_dgr2(NGrid, 3)                              !
    real(RK),    intent(out) :: dF_dtau(NGrid, 2)                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK)                 :: e(NGrid)                                       !
    real(RK)                 :: de_drho(NGrid, 2)                              !
    real(RK)                 :: de_dgr2(NGrid, 3)                              !
    real(RK)                 :: de_dtau(NGrid, 2)                              !
    !                                                                          !
    real(RK)                 :: z(NGrid)                                       !
    real(RK)                 :: dz_drho(NGrid, 2)                              !
    real(RK)                 :: dz_dgr2(NGrid, 3)                              !
    real(RK)                 :: dz_dtau(NGrid, 2)                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: d_par = 2.80_RK                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F       = 0.0_RK                                                           !
    dF_drho = 0.0_RK                                                           !
    dF_dgr2 = 0.0_RK                                                           !
    dF_dtau = 0.0_RK                                                           !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    call revPKZBcorrelation(Flavor, NGrid, rho, gr2, tau,                      &
                                                e, de_drho, de_dgr2, de_dtau)  !
    !                                                                          !
    call tauW_d_tau(NGrid, rho, gr2, tau, z, dz_drho, dz_dgr2, dz_dtau)        !
    !                                                                          !
    !                                                                          !
    F = sum(rho,2) * e * (1.0_RK + d_par * e * z**3)                           !
    !                                                                          !
    dF_drho(:,1) = e * (1.0_RK + d_par * e * z**3)                             &
                 + sum(rho,2) * de_drho(:,1) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_drho(:,1) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_drho(:,1) * z**2 * 3.0_RK    !
    dF_drho(:,2) = e * (1.0_RK + d_par * e * z**3)                             &
                 + sum(rho,2) * de_drho(:,2) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_drho(:,2) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_drho(:,2) * z**2 * 3.0_RK    !
    dF_dgr2(:,1) = sum(rho,2) * de_dgr2(:,1) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_dgr2(:,1) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_dgr2(:,1) * z**2 * 3.0_RK    !
    dF_dgr2(:,2) = sum(rho,2) * de_dgr2(:,2) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_dgr2(:,2) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_dgr2(:,2) * z**2 * 3.0_RK    !
    dF_dgr2(:,3) = sum(rho,2) * de_dgr2(:,3) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_dgr2(:,3) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_dgr2(:,3) * z**2 * 3.0_RK    !
    dF_dtau(:,1) = sum(rho,2) * de_dtau(:,1) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_dtau(:,1) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_dtau(:,1) * z**2 * 3.0_RK    !
    dF_dtau(:,2) = sum(rho,2) * de_dtau(:,2) * (1.0_RK + d_par * e * z**3)     &
                 + sum(rho,2) * e * d_par * de_dtau(:,2) * z**3                &
                 + sum(rho,2) * e**2 * d_par * dz_dtau(:,2) * z**2 * 3.0_RK    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine tauW_d_tau(NGrid, rho, gr2, tau, z, dz_drho, dz_dgr2, dz_dtau)  !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      real(RK),    intent(in)  :: gr2(NGrid, 3)                                !
      real(RK),    intent(in)  :: tau(NGrid, 2)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: z(NGrid)                                     !
      real(RK),    intent(out) :: dz_drho(NGrid, 2)                            !
      real(RK),    intent(out) :: dz_dgr2(NGrid, 3)                            !
      real(RK),    intent(out) :: dz_dtau(NGrid, 2)                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where(sum(tau,2) .ge. DTol .and. sum(rho,2) .ge. DTol)                   !
        !                                                                      !
        z = (gr2(:,1) + gr2(:,2) + 2.0_RK * gr2(:,3))                          &
                                         / (8.0_RK * sum(rho,2) * sum(tau,2))  !
        dz_drho(:,1) = - z / sum(rho,2)                                        !
        dz_drho(:,2) = dz_drho(:,1)                                            !
        dz_dgr2(:,1) = 1.0_RK / (8.0_RK * sum(rho,2) * sum(tau,2))             !
        dz_dgr2(:,2) = 1.0_RK / (8.0_RK * sum(rho,2) * sum(tau,2))             !
        dz_dgr2(:,3) = 2.0_RK / (8.0_RK * sum(rho,2) * sum(tau,2))             !
        dz_dtau(:,1) = - z / sum(tau,2)                                        !
        dz_dtau(:,2) = dz_dtau(:,1)                                            !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        z = 0.0_RK                                                             !
        dz_drho(:,1) = 0.0_RK                                                  !
        dz_drho(:,2) = 0.0_RK                                                  !
        dz_dgr2(:,1) = 0.0_RK                                                  !
        dz_dgr2(:,2) = 0.0_RK                                                  !
        dz_dgr2(:,3) = 0.0_RK                                                  !
        dz_dtau(:,1) = 0.0_RK                                                  !
        dz_dtau(:,2) = 0.0_RK                                                  !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine tauW_d_tau                                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine revPKZBcorrelation(Flavor, NGrid, rho, gr2, tau,                &
                                               e, de_drho, de_dgr2, de_dtau)   !
      !------------------------------------------------------------------------+
      !                                                                        !
      use pbe_gga_module,     only: pbe_correlation                            &
                                  , PBE_key, PBE_revTPSS_key                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: Flavor                                       !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      real(RK),    intent(in)  :: gr2(NGrid, 3)                                !
      real(RK),    intent(in)  :: tau(NGrid, 2)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: e(NGrid)                                     !
      real(RK),    intent(out) :: de_drho(NGrid, 2)                            !
      real(RK),    intent(out) :: de_dgr2(NGrid, 3)                            !
      real(RK),    intent(out) :: de_dtau(NGrid, 2)                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK)              :: s                                            !
      !                                                                        !
      real(RK)                 :: g(NGrid)                                     !
      real(RK)                 :: dg_drho(NGrid, 2)                            !
      real(RK)                 :: dg_dgr2(NGrid, 3)                            !
      !                                                                        !
      real(RK)                 :: z(NGrid)                                     !
      real(RK)                 :: dz_drho(NGrid, 2)                            !
      real(RK)                 :: dz_dgr2(NGrid, 3)                            !
      real(RK)                 :: dz_dtau(NGrid, 2)                            !
      !                                                                        !
      real(RK)                 :: c(NGrid)                                     !
      real(RK)                 :: dc_drho(NGrid, 2)                            !
      real(RK)                 :: dc_dgr2(NGrid, 3)                            !
      !                                                                        !
      real(RK)                 :: gs(NGrid)                                    !
      real(RK)                 :: dgs_drhos(NGrid, 2)                          !
      real(RK)                 :: dgs_dgr2s(NGrid, 3)                          !
      !                                                                        !
      real(RK)                 :: rhos(NGrid, 2)                               !
      real(RK)                 :: gr2s(NGrid, 3)                               !
      !                                                                        !
      integer(IK)              :: Key                                          !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      if (Flavor == TPSS_key) then                                             !
        Key = PBE_key                                                          !
      elseif (Flavor == revTPSS_key) then                                      !
        Key = PBE_revTPSS_key                                                  !
      else                                                                     !
        write(*,*) 'Unknown key for PBE -> stop'; stop                         !
      end if                                                                   !
      !                                                                        !
      call pbe_correlation(Key, .False., NGrid, rho, gr2,                      &
                                                    g, dg_drho, dg_dgr2)       !
      ! convert to energy density                                              !
      where(sum(rho,2) .ge. DTol)                                              !
        !                                                                      !
        g            = g / sum(rho,2)                                          !
        dg_drho(:,1) = (dg_drho(:,1) - g) / sum(rho,2)                         !
        dg_drho(:,2) = (dg_drho(:,2) - g) / sum(rho,2)                         !
        dg_dgr2(:,1) = dg_dgr2(:,1) / sum(rho,2)                               !
        dg_dgr2(:,2) = dg_dgr2(:,2) / sum(rho,2)                               !
        dg_dgr2(:,3) = dg_dgr2(:,3) / sum(rho,2)                               !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        g            = 0.0_RK                                                  !
        dg_drho(:,1) = 0.0_RK                                                  !
        dg_drho(:,2) = 0.0_RK                                                  !
        dg_dgr2(:,1) = 0.0_RK                                                  !
        dg_dgr2(:,2) = 0.0_RK                                                  !
        dg_dgr2(:,3) = 0.0_RK                                                  !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      call tauW_d_tau(NGrid, rho, gr2, tau, z, dz_drho, dz_dgr2, dz_dtau)      !
      !                                                                        !
      call C_term(Flavor, NGrid, rho, gr2, c, dc_drho, dc_dgr2)                !
      !                                                                        !
      ! opposite spin part                                                     !
      e = g * (1.0_RK + c * z**2)                                              !
      !                                                                        !
      de_drho(:,1) = dg_drho(:,1) * (1.0_RK + c * z**2)                        &
                   + g * dc_drho(:,1) * z**2                                   &
                   + g * c * z * 2.0_RK * dz_drho(:,1)                         !
      de_drho(:,2) = dg_drho(:,2) * (1.0_RK + c * z**2)                        &
                   + g * dc_drho(:,2) * z**2                                   &
                   + g * c * z * 2.0_RK * dz_drho(:,2)                         !
      de_dgr2(:,1) = dg_dgr2(:,1) * (1.0_RK + c * z**2)                        &
                   + g * dc_dgr2(:,1) * z**2                                   &
                   + g * c * z * 2.0_RK * dz_dgr2(:,1)                         !
      de_dgr2(:,2) = dg_dgr2(:,2) * (1.0_RK + c * z**2)                        &
                   + g * dc_dgr2(:,2) * z**2                                   &
                   + g * c * z * 2.0_RK * dz_dgr2(:,2)                         !
      de_dgr2(:,3) = dg_dgr2(:,3) * (1.0_RK + c * z**2)                        &
                   + g * dc_dgr2(:,3) * z**2                                   &
                   + g * c * z * 2.0_RK * dz_dgr2(:,3)                         !
      de_dtau(:,1) = g * c * z * 2.0_RK * dz_dtau(:,1)                         !
      de_dtau(:,2) = g * c * z * 2.0_RK * dz_dtau(:,2)                         !
      !                                                                        !
      do s = 1, 2                                                              !
        !                                                                      !
        rhos      = 0.0_RK                                                     !
        rhos(:,s) = rho(:,s)                                                   !
        gr2s      = 0.0_RK                                                     !
        gr2s(:,s) = gr2(:,s)                                                   !
        !                                                                      !
        gs        = 0.0_RK                                                     !
        dgs_drhos = 0.0_RK                                                     !
        dgs_dgr2s = 0.0_RK                                                     !
        !                                                                      !
        call pbe_correlation(Key, .True., NGrid, rhos, gr2s,                   &
                                                  gs, dgs_drhos, dgs_dgr2s)    !
        ! convert to energy density                                            !
        where(sum(rhos,2) .ge. DTol)                                           !
          !                                                                    !
          gs             = gs / sum(rhos,2)                                    !
          dgs_drhos(:,s) = (dgs_drhos(:,s) - gs) / sum(rhos,2)                 !
          dgs_dgr2s(:,s) = dgs_dgr2s(:,s) / sum(rhos,2)                        !
          !                                                                    !
        elsewhere                                                              !
          !                                                                    !
          gs             = 0.0_RK                                              !
          dgs_drhos(:,s) = 0.0_RK                                              !
          dgs_dgr2s(:,s) = 0.0_RK                                              !
          !                                                                    !
        end where                                                              !
        !                                                                      !
        where(g .gt. gs)                                                       !
          !                                                                    !
          gs             = g                                                   !
          dgs_drhos(:,1) = dg_drho(:,1)                                        !
          dgs_drhos(:,2) = dg_drho(:,2)                                        !
          dgs_dgr2s(:,1) = dg_dgr2(:,1)                                        !
          dgs_dgr2s(:,2) = dg_dgr2(:,2)                                        !
          dgs_dgr2s(:,3) = dg_dgr2(:,3)                                        !
          !                                                                    !
        end where                                                              !
        !                                                                      !
        where(sum(rho,2) .gt. DTol)                                            !
          !                                                                    !
          e = e - (1.0_RK + c) * z**2 * rho(:,s) / sum(rho,2) * gs             !
          !                                                                    !
          de_drho(:,1) = de_drho(:,1)                                          &
                       - rho(:,s) / sum(rho,2) * gs * dc_drho(:,1) * z**2      &
                       + (1.0_RK + c) * (z**2 * rho(:,s) / (sum(rho,2)**2) * gs&
                       - 2.0_RK * z * dz_drho(:,1) * rho(:,s) / sum(rho,2) * gs&
                       - z**2 * rho(:,s) / sum(rho,2) * dgs_drhos(:,1))        !
          de_drho(:,2) = de_drho(:,2)                                          &
                       - rho(:,s) / sum(rho,2) * gs * dc_drho(:,2) * z**2      &
                       + (1.0_RK + c) * (z**2 * rho(:,s) / (sum(rho,2)**2) * gs&
                       - 2.0_RK * z * dz_drho(:,2) * rho(:,s) / sum(rho,2) * gs&
                       - z**2 * rho(:,s) / sum(rho,2) * dgs_drhos(:,2))        !
          de_drho(:,s) = de_drho(:,s)                                          &
                       - (1.0_RK + c) * z**2 / sum(rho,2) * gs                 !
          !                                                                    !
          de_dgr2(:,1) = de_dgr2(:,1)                                          &
                       - dc_dgr2(:,1) * z**2 * rho(:,s) / sum(rho,2) * gs      &
                       - (1.0_RK + c) * 2.0_RK * z * dz_dgr2(:,1)              &
                                                  * rho(:,s) / sum(rho,2) * gs &
                       - (1.0_RK + c) * z**2 * rho(:,s) / sum(rho,2)           &
                                                              * dgs_dgr2s(:,1) !
          de_dgr2(:,2) = de_dgr2(:,2)                                          &
                       - dc_dgr2(:,2) * z**2 * rho(:,s) / sum(rho,2) * gs      &
                       - (1.0_RK + c) * 2.0_RK * z * dz_dgr2(:,2)              &
                                                  * rho(:,s) / sum(rho,2) * gs &
                       - (1.0_RK + c) * z**2 * rho(:,s) / sum(rho,2)           &
                                                              * dgs_dgr2s(:,2) !
          de_dgr2(:,3) = de_dgr2(:,3)                                          &
                       - dc_dgr2(:,3) * z**2 * rho(:,s) / sum(rho,2) * gs      &
                       - (1.0_RK + c) * 2.0_RK * z * dz_dgr2(:,3)              &
                                                  * rho(:,s) / sum(rho,2) * gs &
                       - (1.0_RK + c) * z**2 * rho(:,s) / sum(rho,2)           &
                                                              * dgs_dgr2s(:,3) !
          !                                                                    !
          de_dtau(:,1) = de_dtau(:,1)                                          &
                       - (1.0_RK + c) * 2.0_RK * z * dz_dtau(:,1)              &
                                                  * rho(:,s) / sum(rho,2) * gs !
          de_dtau(:,2) = de_dtau(:,2)                                          &
                       - (1.0_RK + c) * 2.0_RK * z * dz_dtau(:,2)              &
                                                  * rho(:,s) / sum(rho,2) * gs !
          !                                                                    !
        end where                                                              !
        !                                                                      !
      end do                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine revPKZBcorrelation                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine C_term(Flavor, NGrid, rho, gr2,            c, dc_drho, dc_dgr2) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: Flavor                                       !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      real(RK),    intent(in)  :: gr2(NGrid, 3)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: c(NGrid)                                     !
      real(RK),    intent(out) :: dc_drho(NGrid, 2)                            !
      real(RK),    intent(out) :: dc_dgr2(NGrid, 3)                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                 :: y(NGrid)                                     !
      real(RK)                 :: dy_drho(NGrid, 2)                            !
      !                                                                        !
      real(RK)                 :: x2(NGrid)                                    !
      real(RK)                 :: dx2_drho(NGrid, 2)                           !
      real(RK)                 :: dx2_dgr2(NGrid, 3)                           !
      !                                                                        !
      real(RK)                 :: c_pars(4)                                    !
      !                                                                        !
      real(RK)                 :: denom(NGrid)                                 !
      real(RK)                 :: ddenom_dy(NGrid)                             !
      real(RK)                 :: ddenom_dx2(NGrid)                            !
      !                                                                        !
      real(RK)                 :: c0(NGrid)                                    !
      real(RK)                 :: dc0_dy(NGrid)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      call relative_spin_polarization(NGrid, rho, y, dy_drho)                  !
      !                                                                        !
      call spin_polarized_gradient(NGrid, rho, gr2, x2, dx2_drho, dx2_dgr2)    !
      !                                                                        !
      call get_C_term_parameters(Flavor, c_pars)                               !
      !                                                                        !
      where(y .lt. 1.0_RK - DTol .and. y .gt. -1.0_RK + DTol)                  !
        !                                                                      !
        c0 = c_pars(1) + c_pars(2) * y**2 + c_pars(3) * y**4 + c_pars(4) * y**6!
        !                                                                      !
        dc0_dy = c_pars(2) * 2.0_RK * y + c_pars(3) * 4.0_RK * y**3            &
               + c_pars(4) * 6.0_RK * y**5                                     !
        !                                                                      !
        denom = 1.0 / (1.0_RK + x2 * ((1.0_RK + y + DTol)**(-4.0_RK/3.0_RK)    &
                                    + (1.0_RK - y + DTol)**(-4.0_RK/3.0_RK))   &
                                                         * 0.5_RK)             !
        ddenom_dy  = - denom**2 * 0.5_RK * x2                                  &
                   * ((-4.0_RK/3.0_RK) * (1.0_RK + y + DTol)**(-7.0_RK/3.0_RK) &
                     -(-4.0_RK/3.0_RK) * (1.0_RK - y + DTol)**(-7.0_RK/3.0_RK))!
        ddenom_dx2 = - denom**2 * 0.5_RK                                       &
                   * ((1.0_RK + y + DTol)**(-4.0_RK/3.0_RK)                    &
                     +(1.0_RK - y + DTol)**(-4.0_RK/3.0_RK))                   !
        !                                                                      !
        c = c0 * denom**4                                                      !
        !                                                                      !
        dc_drho(:,1) = denom**4 * dc0_dy * dy_drho(:,1)                        &
                     + 4.0_RK * c0 * denom**3 * ddenom_dy * dy_drho(:,1)       &
                     + 4.0_RK * c0 * denom**3 * ddenom_dx2 * dx2_drho(:,1)     !
        dc_drho(:,2) = denom**4 * dc0_dy * dy_drho(:,2)                        &
                     + 4.0_RK * c0 * denom**3 * ddenom_dy * dy_drho(:,2)       &
                     + 4.0_RK * c0 * denom**3 * ddenom_dx2 * dx2_drho(:,2)     !
        !                                                                      !
        dc_dgr2(:,1) = 4.0_RK * c0 * denom**3 * ddenom_dx2 * dx2_dgr2(:,1)     !
        dc_dgr2(:,2) = 4.0_RK * c0 * denom**3 * ddenom_dx2 * dx2_dgr2(:,2)     !
        dc_dgr2(:,3) = 4.0_RK * c0 * denom**3 * ddenom_dx2 * dx2_dgr2(:,3)     !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        c = 0.0_RK                                                             !
        dc_drho(:,1) = 0.0_RK                                                  !
        dc_drho(:,2) = 0.0_RK                                                  !
        dc_dgr2(:,1) = 0.0_RK                                                  !
        dc_dgr2(:,2) = 0.0_RK                                                  !
        dc_dgr2(:,3) = 0.0_RK                                                  !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine C_term                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine get_C_term_parameters(Flavor, c_pars)                           !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: Flavor                                       !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: c_pars(4)                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (Flavor)                                                     !
        case(TPSS_key)  ! TPSS                                                 !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          c_pars(1) = 0.53_RK                                                  !
          c_pars(2) = 0.87_RK                                                  !
          c_pars(3) = 0.50_RK                                                  !
          c_pars(4) = 2.26_RK                                                  !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(revTPSS_key)  ! revTPSS                                           !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          c_pars(1) = 0.59_RK                                                  !
          c_pars(2) = 0.9269_RK                                                !
          c_pars(3) = 0.6225_RK                                                !
          c_pars(4) = 2.1540_RK                                                !
       case default
          c_pars(:) = -1.0
          ABORT("no such case")
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_C_term_parameters                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine relative_spin_polarization(NGrid, rho,              y, dy_drho) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: y(NGrid)                                     !
      real(RK),    intent(out) :: dy_drho(NGrid, 2)                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (sum(rho, 2) .ge. DTol)                                            !
        !                                                                      !
        y = (rho(:,1) - rho(:,2)) / sum(rho, 2)                                !
        !                                                                      !
        dy_drho(:, 1) =   2.0_RK * rho(:,2) / ((sum(rho, 2))**2)               !
        dy_drho(:, 2) = - 2.0_RK * rho(:,1) / ((sum(rho, 2))**2)               !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        y = 0.0_RK                                                             !
        !                                                                      !
        dy_drho(:, 1) = 0.0_RK                                                 !
        dy_drho(:, 2) = 0.0_RK                                                 !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine relative_spin_polarization                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine spin_polarized_gradient(NGrid, rho, gr2, x2, dx2_drho, dx2_dgr2)!
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      real(RK),    intent(in)  :: gr2(NGrid, 3)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: x2(NGrid)                                    !
      real(RK),    intent(out) :: dx2_drho(NGrid, 2)                           !
      real(RK),    intent(out) :: dx2_dgr2(NGrid, 3)                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                 :: denom(NGrid)                                 !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (sum(rho, 2) .ge. DTol)                                            !
        !                                                                      !
        denom = 1.0_RK / ((3.0_RK * pi**2 * sum(rho, 2))**(2.0_RK / 3.0_RK)    &
                                                          * sum(rho, 2)**4)    !
        !                                                                      !
        x2 = (rho(:,1)**2 * gr2(:,2) + rho(:,2)**2 * gr2(:,1)                  &
                    - 2.0_RK * rho(:,1) * rho(:,2) * gr2(:,3)) * denom         !
        !                                                                      !
        dx2_drho(:, 1) = 2.0_RK * (rho(:,1) * gr2(:,2) - rho(:,2) * gr2(:,3))  &
                       * denom                                                 &
                       - (14.0_RK / 3.0_RK) * x2 / sum(rho, 2)                 !
        dx2_drho(:, 2) = 2.0_RK * (rho(:,2) * gr2(:,1) - rho(:,1) * gr2(:,3))  &
                       * denom                                                 &
                       - (14.0_RK / 3.0_RK) * x2 / sum(rho, 2)                 !
        dx2_dgr2(:, 1) = rho(:,2)**2 * denom                                   !
        dx2_dgr2(:, 2) = rho(:,1)**2 * denom                                   !
        dx2_dgr2(:, 3) = - 2.0_RK * rho(:,1) * rho(:,2) * denom                !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        x2 = 0.0_RK                                                            !
        !                                                                      !
        dx2_drho(:, 1) = 0.0_RK                                                !
        dx2_drho(:, 2) = 0.0_RK                                                !
        dx2_dgr2(:, 1) = 0.0_RK                                                !
        dx2_dgr2(:, 2) = 0.0_RK                                                !
        dx2_dgr2(:, 3) = 0.0_RK                                                !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine spin_polarized_gradient                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine tpss_correlation                                              !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module tpss
! Default options for vim:sw=2:expandtab:smarttab:autoindent:syntax

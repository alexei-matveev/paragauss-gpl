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
module m06
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
  use type_module,        only: IK => i4_kind, RK => r8_kind                   !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  public :: m06_X                                                              &
          , m06_C                                                              !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  real(RK), parameter      :: DTol = 1.0E-30_RK                                !
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
  subroutine m06_X(PSet, NGrid, ISpin, rho, gr2, tau,                          &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! PSet denotes parameter set:           1  -- M06-L                        !
    !                                       2  -- M06                          !
    !                                       3  -- M06-2X                       !
    !                                       4  -- M06-HF                       !
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    integer(IK), intent(in)  :: PSet                                           !
    integer(IK), intent(in)  :: NGrid                                          !
    integer(IK), intent(in)  :: ISpin                                          !
    real(RK),    intent(in)  :: rho(:,:)                                       !
    real(RK),    intent(in)  :: gr2(:,:)                                       !
    real(RK),    intent(in)  :: tau(:,:)                                       !
    !                                                                          !
    real(RK),  intent(inout) :: F(:)                                           !
    real(RK),  intent(inout) :: dF_drho(:,:)                                   !
    real(RK),  intent(inout) :: dF_dgr2(:,:)                                   !
    real(RK),  intent(inout) :: dF_dtau(:,:)                                   !
    !                                                                          !
    real(RK)                 :: rho_tmp(NGrid, 2)                              !
    real(RK)                 :: gr2_tmp(NGrid, 3)                              !
    real(RK)                 :: tau_tmp(NGrid, 2)                              !
    real(RK)                 :: F_tmp(NGrid)                                   !
    real(RK)                 :: dF_drho_tmp(NGrid, 2)                          !
    real(RK)                 :: dF_dgr2_tmp(NGrid, 3)                          !
    real(RK)                 :: dF_dtau_tmp(NGrid, 2)                          !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    IF (ISpin .eq. 1) THEN                                                     !
      !                                                                        !
      rho_tmp(:, 1) = rho(:NGrid, 1) * 0.5_RK                                  !
      rho_tmp(:, 2) = rho(:NGrid, 1) * 0.5_RK                                  !
      gr2_tmp(:, 1) = gr2(:NGrid, 1) * 0.25_RK                                 !
      gr2_tmp(:, 2) = gr2(:NGrid, 1) * 0.25_RK                                 !
      gr2_tmp(:, 3) = gr2(:NGrid, 1) * 0.25_RK                                 !
      tau_tmp(:, 1) = tau(:NGrid, 1) * 0.5_RK                                  !
      tau_tmp(:, 2) = tau(:NGrid, 1) * 0.5_RK                                  !
      !                                                                        !
      call m06_exchange( PSet, NGrid, 2, rho_tmp, gr2_tmp, tau_tmp,            &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      F(:NGrid)          = F(:NGrid) + F_tmp                                   !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + dF_drho_tmp(:, 1)              !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + dF_dgr2_tmp(:, 1) * 0.5_RK     !
      dF_dtau(:NGrid, 1) = dF_dtau(:NGrid, 1) + dF_dtau_tmp(:, 1)              !
      !                                                                        !
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !                                                                        !
      call m06_exchange( PSet, NGrid, 2                                        &
                       , rho(:NGrid, :2), gr2(:NGrid, :3), tau(:NGrid, :2),    &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      F(:NGrid)           = F(:NGrid) + F_tmp                                  !
      dF_drho(:NGrid, :2) = dF_drho(:NGrid, :2) + dF_drho_tmp                  !
      dF_dgr2(:NGrid, :3) = dF_dgr2(:NGrid, :3) + dF_dgr2_tmp                  !
      dF_dtau(:NGrid, :2) = dF_dtau(:NGrid, :2) + dF_dtau_tmp                  !
      !                                                                        !
    ELSE                                                                       !
      !                                                                        !
      WRITE(*,*) 'ISpin can only be either 1 or 2'                             !
      STOP                                                                     !
      !                                                                        !
    END IF                                                                     !
    !                                                                          !
  end subroutine m06_X                                                         !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine m06_C(PSet, NGrid, ISpin, rho, gr2, tau,                          &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! PSet denotes parameter set:           1  -- M06-L                        !
    !                                       2  -- M06                          !
    !                                       3  -- M06-2X                       !
    !                                       4  -- M06-HF                       !
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    integer(IK), intent(in)  :: PSet                                           !
    integer(IK), intent(in)  :: NGrid                                          !
    integer(IK), intent(in)  :: ISpin                                          !
    real(RK),    intent(in)  :: rho(:,:)                                       !
    real(RK),    intent(in)  :: gr2(:,:)                                       !
    real(RK),    intent(in)  :: tau(:,:)                                       !
    !                                                                          !
    real(RK),  intent(inout) :: F(:)                                           !
    real(RK),  intent(inout) :: dF_drho(:,:)                                   !
    real(RK),  intent(inout) :: dF_dgr2(:,:)                                   !
    real(RK),  intent(inout) :: dF_dtau(:,:)                                   !
    !                                                                          !
    real(RK)                 :: rho_tmp(NGrid, 2)                              !
    real(RK)                 :: gr2_tmp(NGrid, 3)                              !
    real(RK)                 :: tau_tmp(NGrid, 2)                              !
    real(RK)                 :: F_tmp(NGrid)                                   !
    real(RK)                 :: dF_drho_tmp(NGrid, 2)                          !
    real(RK)                 :: dF_dgr2_tmp(NGrid, 3)                          !
    real(RK)                 :: dF_dtau_tmp(NGrid, 2)                          !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    !                                                                          !
    IF (ISpin .eq. 1) THEN                                                     !
      !                                                                        !
      rho_tmp(:, 1) = rho(:NGrid, 1) * 0.5_RK                                  !
      rho_tmp(:, 2) = rho(:NGrid, 1) * 0.5_RK                                  !
      gr2_tmp(:, 1) = gr2(:NGrid, 1) * 0.25_RK                                 !
      gr2_tmp(:, 2) = gr2(:NGrid, 1) * 0.25_RK                                 !
      gr2_tmp(:, 3) = gr2(:NGrid, 1) * 0.25_RK                                 !
      tau_tmp(:, 1) = tau(:NGrid, 1) * 0.5_RK                                  !
      tau_tmp(:, 2) = tau(:NGrid, 1) * 0.5_RK                                  !
      !                                                                        !
      call m06_correlation(PSet, NGrid, 2, rho_tmp, gr2_tmp, tau_tmp,          &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp)  !
      !                                                                        !
      F(:NGrid)          = F(:NGrid) + F_tmp                                   !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + dF_drho_tmp(:, 1)              !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + dF_dgr2_tmp(:, 1) * 0.5_RK     !
      dF_dtau(:NGrid, 1) = dF_dtau(:NGrid, 1) + dF_dtau_tmp(:, 1)              !
      !                                                                        !
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !                                                                        !
      call m06_correlation( PSet, NGrid, 2                                     &
                          , rho(:NGrid, :2), gr2(:NGrid, :3), tau(:NGrid, :2), &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      F(:NGrid)           = F(:NGrid) + F_tmp                                  !
      dF_drho(:NGrid, :2) = dF_drho(:NGrid, :2) + dF_drho_tmp                  !
      dF_dgr2(:NGrid, :3) = dF_dgr2(:NGrid, :3) + dF_dgr2_tmp                  !
      dF_dtau(:NGrid, :2) = dF_dtau(:NGrid, :2) + dF_dtau_tmp                  !
      !                                                                        !
    ELSE                                                                       !
      !                                                                        !
      WRITE(*,*) 'ISpin can only be either 1 or 2'                             !
      STOP                                                                     !
      !                                                                        !
    END IF                                                                     !
    !                                                                          !
  end subroutine m06_C                                                         !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine m06_exchange(PSet, NGrid, ISpin, rho, gr2, tau,                   &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    use vsxc,               only: eX_lsda                                      &
                                , vsxc_exchange                                !
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
    real(RK),    intent(in)  :: rho(:, :) ! (NGrid, ISpin)
    real(RK),    intent(in)  :: gr2(:, :) ! (NGrid, 2*ISpin-1)
    real(RK),    intent(in)  :: tau(:, :) ! (NGrid, ISpin)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: F(:) ! (NGrid)
    real(RK),    intent(out) :: dF_drho(:, :) ! (NGrid, ISpin)
    real(RK),    intent(out) :: dF_dgr2(:, :) ! (NGrid, 2*ISpin-1)
    real(RK),    intent(out) :: dF_dtau(:, :) ! (NGrid, ISpin)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)              :: s                                              !
    !                                                                          !
    real(RK)                 :: a(0:11)                                        !
    !                                                                          !
    real(RK)                 :: e(NGrid)                                       !
    real(RK)                 :: de_drho(NGrid)                                 !
    !                                                                          !
    real(RK)                 :: p(NGrid)                                       !
    real(RK)                 :: dp_drho(NGrid)                                 !
    real(RK)                 :: dp_dgr2(NGrid)                                 !
    !                                                                          !
    real(RK)                 :: k(NGrid)                                       !
    real(RK)                 :: dk_drho(NGrid)                                 !
    real(RK)                 :: dk_dtau(NGrid)                                 !
    !                                                                          !
    real(RK)                 :: V(NGrid)                                       !
    real(RK)                 :: dV_drho(NGrid, ISpin)                          !
    real(RK)                 :: dV_dgr2(NGrid, 2*ISpin-1)                      !
    real(RK)                 :: dV_dtau(NGrid, ISpin)                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F       = 0.0_RK                                                           !
    dF_drho = 0.0_RK                                                           !
    dF_dgr2 = 0.0_RK                                                           !
    dF_dtau = 0.0_RK                                                           !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! M05 term                                                                 !
    !--------------------------------------------------------------------------+
    call get_pars_m06_exchange(PSET, a)                                        !
    !                                                                          !
    spin: DO s = 1, ISpin                                                      !
      !------------------------------------------------------------------------+
      !                                                                        !
      call eX_lsda(NGrid, rho(:NGrid, s), e, de_drho)                          !
      call pbeef(NGrid, rho(:NGrid, s), gr2(:NGrid, s), p, dp_drho, dp_dgr2)   !
      call kedef(NGrid, a, rho(:NGrid, s), tau(:NGrid, s), k, dk_drho, dk_dtau)!
      !                                                                        !
      F = F + e * p * k                                                        !
      !                                                                        !
      dF_drho(1:NGrid, s) = de_drho * p * k                                    &
                          + e * dp_drho * k                                    &
                          + e * p * dk_drho                                    !
      dF_dgr2(1:NGrid, s) = e * dp_dgr2 * k                                    !
      dF_dtau(1:NGrid, s) = e * p * dk_dtau                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO spin                                                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! VSXC term                                                                !
    !--------------------------------------------------------------------------+
    call vsxc_exchange(PSet, NGrid, ISpin, rho, gr2, tau,                      &
                                                V, dV_drho, dV_dgr2, dV_dtau)  !
    !                                                                          !
    F       = F       + V                                                      !
    dF_drho = dF_drho + dV_drho                                                !
    dF_dgr2 = dF_dgr2 + dV_dgr2                                                !
    dF_dtau = dF_dtau + dV_dtau                                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    subroutine get_pars_m06_exchange(PSET, a)                                  !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: a(0:11)                                      !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (PSet)                                                       !
        case(1)  ! M06-L                                                       !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          a = (/  3.987756E-01_RK                                              &
               ,  2.548219E-01_RK                                              &
               ,  3.923994E-01_RK                                              &
               , -2.103655E+00_RK                                              &
               , -6.302147E+00_RK                                              &
               ,  1.097615E+01_RK                                              &
               ,  3.097273E+01_RK                                              &
               , -2.318489E+01_RK                                              &
               , -5.673480E+01_RK                                              &
               ,  2.160364E+01_RK                                              &
               ,  3.421814E+01_RK                                              &
               , -9.049762E+00_RK /)                                           !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(2)  ! M06                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          a = (/  5.877943E-01_RK                                              &
               , -1.371776E-01_RK                                              &
               ,  2.682367E-01_RK                                              &
               , -2.515898E+00_RK                                              &
               , -2.978892E+00_RK                                              &
               ,  8.710679E+00_RK                                              &
               ,  1.688195E+01_RK                                              &
               , -4.489724E+00_RK                                              &
               , -3.299983E+01_RK                                              &
               , -1.449050E+01_RK                                              &
               ,  2.043747E+01_RK                                              &
               ,  1.256504E+01_RK /)                                           !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(3)  ! M06-2X                                                      !
          !                                                                    !
          a = (/  4.600000E-01_RK                                              &
               , -2.206052E-01_RK                                              &
               , -9.431788E-02_RK                                              &
               ,  2.164494E+00_RK                                              &
               , -2.556466E+00_RK                                              &
               , -1.422133E+01_RK                                              &
               ,  1.555044E+01_RK                                              &
               ,  3.598078E+01_RK                                              &
               , -2.722754E+01_RK                                              &
               , -3.924093E+01_RK                                              &
               ,  1.522808E+01_RK                                              &
               ,  1.522227E+01_RK /)                                           !
          !                                                                    !
        case(4)  ! M06-HF                                                      !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          a = (/  1.179732E-01_RK                                              &
               , -1.066708E+00_RK                                              &
               , -1.462405E-01_RK                                              &
               ,  7.481848E+00_RK                                              &
               ,  3.776679E+00_RK                                              &
               , -4.436118E+01_RK                                              &
               , -1.830962E+01_RK                                              &
               ,  1.003903E+02_RK                                              &
               ,  3.864360E+01_RK                                              &
               , -9.806018E+01_RK                                              &
               , -2.557716E+01_RK                                              &
               ,  3.590404E+01_RK /)                                           !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_m06_exchange                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine pbeef(NGrid, rho, gr2, p, dp_drho, dp_dgr2)                     !
      !------------------------------------------------------------------------+
      !                                                                        !
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
      real(RK)                 :: x2(NGrid)                                    !
      real(RK)                 :: dx2_drho(NGrid)                              !
      real(RK)                 :: dx2_dgr2(NGrid)                              !
      real(RK)                 :: dp_dx2(NGrid)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK), parameter      :: pbec2 = 4.49267E-3_RK                        !
      real(RK), parameter      :: pbec1 = 3.612108584108E-03_RK                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (rho .ge. Dtol)                                                    !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        x2       = gr2 / (rho**(8.0_RK / 3.0_RK))                              !
        dx2_drho = - (8.0_RK / 3.0_RK) * gr2 / (rho**(11.0_RK / 3.0_RK))       !
        dx2_dgr2 = 1.0_RK / (rho**(8.0_RK / 3.0_RK))                           !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      elsewhere
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        x2       = 0.0_RK                                                      !
        dx2_drho = 0.0_RK                                                      !
        dx2_dgr2 = 0.0_RK                                                      !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end where                                                                !
      !                                                                        !
      p       = 1.0_RK + pbec1 * x2 / (1.0_RK + pbec2 * x2)                    !
      !                                                                        !
      dp_dx2  = pbec1 / (1.0_RK + pbec2 * x2)                                  &
              - pbec1 * x2 * pbec2 / ((1.0_RK + pbec2 * x2)**2)                !
      dp_drho = dp_dx2 * dx2_drho                                              !
      dp_dgr2 = dp_dx2 * dx2_dgr2                                              !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine pbeef                                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine kedef(NGrid, a, rho, tau, k, dk_drho, dk_dtau)                  !
      !------------------------------------------------------------------------+
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: a(0:11)                                      !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      real(RK),    intent(in)  :: tau(NGrid)                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: k(NGrid)                                     !
      real(RK),    intent(out) :: dk_drho(NGrid)                               !
      real(RK),    intent(out) :: dk_dtau(NGrid)                               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK)              :: j                                            !
      !                                                                        !
      real(RK)                 :: t(NGrid)                                     !
      real(RK)                 :: dt_drho(NGrid)                               !
      real(RK)                 :: w(NGrid)                                     !
      real(RK)                 :: dw_dt(NGrid)                                 !
      real(RK)                 :: dw_dtau(NGrid)                               !
      real(RK)                 :: dk_dw(NGrid)                                 !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK), parameter      :: CT = 4.55779987234560_RK ! 3/10(6 pi^2)^(2/3)!
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      WHERE (rho .ge. Dtol .and. tau .ge. DTol)                                !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        t       = CT * rho**(5.0_RK/3.0_RK)                                    !
        dt_drho = CT * (5.0_RK/3.0_RK) * rho**(2.0_RK/3.0_RK)                  !
        !                                                                      !
        w       = (t - tau) / (t + tau)                                        !
        dw_dt   =   1.0_RK / (t + tau) - (t - tau) / ((t + tau)**2)            !
        dw_dtau = - 1.0_RK / (t + tau) - (t - tau) / ((t + tau)**2)            !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      ELSEWHERE                                                                !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        t       = 0.0_RK                                                       !
        dt_drho = 0.0_RK                                                       !
        !                                                                      !
        w       = 0.0_RK                                                       !
        dw_dt   = 0.0_RK                                                       !
        dw_dtau = 0.0_RK                                                       !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      END WHERE                                                                !
      !                                                                        !
      k       = a(11)                                                          !
      DO j = 10, 0, -1                                                         !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        k = k * w + a(j)                                                       !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      END DO                                                                   !
      !                                                                        !
      dk_dw   = 11.0_RK * a(11)                                                !
      DO j = 10, 1, -1                                                         !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        !                                                                      !
        dk_dw = dk_dw * w + real(j, RK) * a(j)                                 !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      END DO                                                                   !
      !                                                                        !
      dk_drho = dk_dw * dw_dt * dt_drho                                        !
      dk_dtau = dk_dw * dw_dtau                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine kedef                                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine m06_exchange                                                  !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine m06_correlation(PSet, NGrid, ISpin, rho, gr2, tau,                & 
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    use pw_ldac_module,     only: pw_ldac                                      !
    use vsxc,               only: Dsic_term                                    &
                                , vsxc_correlation                             !
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
    real(RK),    intent(in)  :: rho(:, :) ! (NGrid, ISpin)
    real(RK),    intent(in)  :: gr2(:, :) ! (NGrid, 2*ISpin-1)
    real(RK),    intent(in)  :: tau(:, :) ! (NGrid, ISpin)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: F(:) ! (NGrid)
    real(RK),    intent(out) :: dF_drho(:, :) ! (NGrid, ISpin)
    real(RK),    intent(out) :: dF_dgr2(:, :) ! (NGrid, 2*ISpin-1)
    real(RK),    intent(out) :: dF_dtau(:, :) ! (NGrid, ISpin)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)              :: s                                              !
    !                                                                          !
    real(RK)                 :: c0, c1, c2, c3, c4                             !
    !                                                                          !
    real(RK)                 :: eo(NGrid)                                      !
    real(RK)                 :: deo_drho(NGrid, ISpin)                         !
    real(RK)                 :: rho_d_2(NGrid, 2)                              !
    real(RK)                 :: es(NGrid)                                      !
    real(RK)                 :: des_drho(NGrid, 2)                             !
    real(RK)                 :: g(NGrid)                                       !
    real(RK)                 :: dg_drho(NGrid)                                 !
    real(RK)                 :: dg_dgr2(NGrid)                                 !
    real(RK)                 :: D(NGrid)                                       !
    real(RK)                 :: dD_drho(NGrid)                                 !
    real(RK)                 :: dD_dgr2(NGrid)                                 !
    real(RK)                 :: dD_dtau(NGrid)                                 !
    real(RK)                 :: dgos_drho(NGrid, ISpin)                        !
    real(RK)                 :: dgos_dgr2(NGrid, ISpin)                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: betss = 0.06_RK                                !
    real(RK), parameter      :: betos = 0.0031_RK                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F       = 0.0_RK                                                           !
    dF_drho = 0.0_RK                                                           !
    dF_dgr2 = 0.0_RK                                                           !
    dF_dtau = 0.0_RK                                                           !
    !                                                                          !
    eo       = 0.0_RK                                                          !
    deo_drho = 0.0_RK                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! VSXC part                                                                !
    !--------------------------------------------------------------------------+
    call vsxc_correlation(PSet, NGrid, ISpin, rho, gr2, tau,                   &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Opposite spin, LDA part                                                  !
    !--------------------------------------------------------------------------+
    call pw_ldac(NGrid, ISpin, rho, eo, deo_drho)                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Same spin part                                                           !
    !--------------------------------------------------------------------------+
    call get_pars_m06_correlation(PSET, 1, c0, c1, c2, c3, c4)                 !
    !                                                                          !
    spin: DO s = 1, ISpin                                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      ! LDA part                                                               !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      !                                                                        !
      rho_d_2(:, s)   = rho(1:NGrid, s)                                        !
      rho_d_2(:, 3-s) = 0.0_RK                                                 !
      !                                                                        !
      es       = 0.0_RK                                                        !
      des_drho = 0.0_RK                                                        !
      !                                                                        !
      call pw_ldac(NGrid, 2, rho_d_2, es, des_drho)                            !
      !                                                                        !
      ! Update opposite spin LDA part                                          !
      eo = eo - es                                                             !
      deo_drho(1:NGrid, s) = deo_drho(1:NGrid, s) - des_drho(1:NGrid, s)       !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      ! MGGA part                                                              !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      call bmk_term(NGrid, 1, betss, c0, c1, c2, c3, c4, rho(:, s), gr2(:, s), &
                                                          g, dg_drho, dg_dgr2) !
      !                                                                        !
      call Dsic_term(NGrid, rho(1:NGrid, s), gr2(1:NGrid, s), tau(1:NGrid, s), &
                                                 D, dD_drho, dD_dgr2, dD_dtau) !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      ! Summarize same spin contributions                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      F                   = F + es * g * D                                     !
      !                                                                        !
      dF_drho(1:NGrid, s) = dF_drho(1:NGrid, s)                                &
                          + des_drho(1:NGrid, s) * g * D                       &
                          + es * dg_drho * D                                   &
                          + es * g * dD_drho                                   !
      !                                                                        !
      dF_dgr2(1:NGrid, s) = dF_dgr2(1:NGrid, s)                                &
                          + es * dg_dgr2 * D                                   &
                          + es * g * dD_dgr2                                   !
      !                                                                        !
      dF_dtau(1:NGrid, s) = dF_dtau(1:NGrid, s)                                &
                          + es * g * dD_dtau                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO spin                                                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Opposite spin MGGA part                                                  !
    !--------------------------------------------------------------------------+
    call get_pars_m06_correlation(PSET, 2, c0, c1, c2, c3, c4)                 !
    !                                                                          !
    call bmk_term(NGrid, ISpin, betos, c0, c1, c2, c3, c4, rho, gr2,           &
                                                      g, dgos_drho, dgos_dgr2) !
    !                                                                          !
    F = F + eo * g                                                             !
    !                                                                          !
    DO s = 1, ISpin                                                            !
      !------------------------------------------------------------------------+
      !                                                                        !
      dF_drho(1:NGrid, s) = dF_drho(1:NGrid, s)                                &
                          + deo_drho(1:NGrid, s) * g                           &
                          + eo * dgos_drho(1:NGrid, s)                         !
      !                                                                        !
      dF_dgr2(1:NGrid, s) = dF_dgr2(1:NGrid, s)                                &
                          + eo * dgos_dgr2(1:NGrid, s)                         !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO                                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    subroutine get_pars_m06_correlation(PSET, term, c0, c1, c2, c3, c4)        !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      integer(IK), intent(in)  :: term  ! 1 = same spin,  2 = opposite spin    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: c0                                           !
      real(RK),    intent(out) :: c1                                           !
      real(RK),    intent(out) :: c2                                           !
      real(RK),    intent(out) :: c3                                           !
      real(RK),    intent(out) :: c4                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (PSet)                                                       !
        case(1)  ! M06-L                                                       !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            c0 =  5.349466E-01_RK                                              !
            c1 =  5.396620E-01_RK                                              !
            c2 = -3.161217E+01_RK                                              !
            c3 =  5.149592E+01_RK                                              !
            c4 = -2.919613E+01_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            c0 =  6.042374E-01_RK                                              !
            c1 =  1.776783E+02_RK                                              !
            c2 = -2.513252E+02_RK                                              !
            c3 =  7.635173E+01_RK                                              !
            c4 = -1.255699E+01_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(2)  ! M06                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            c0 =  5.094055E-01_RK                                              !
            c1 = -1.491085E+00_RK                                              !
            c2 =  1.723922E+01_RK                                              !
            c3 = -3.859018E+01_RK                                              !
            c4 =  2.845044E+01_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            c0 =  3.741539E+00_RK                                              !
            c1 =  2.187098E+02_RK                                              !
            c2 = -4.531252E+02_RK                                              !
            c3 =  2.936479E+02_RK                                              !
            c4 = -6.287470E+01_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(3)  ! M06-2X                                                      !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            c0 =  3.097855E-01_RK                                              !
            c1 = -5.528642E+00_RK                                              !
            c2 =  1.347420E+01_RK                                              !
            c3 = -3.213623E+01_RK                                              !
            c4 =  2.846742E+01_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            c0 =  8.833596E-01_RK                                              !
            c1 =  3.357972E+01_RK                                              !
            c2 = -7.043548E+01_RK                                              !
            c3 =  4.978271E+01_RK                                              !
            c4 = -1.852891E+01_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(4)  ! M06-HF                                                      !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            c0 =  1.023254E-01_RK                                              !
            c1 = -2.453783E+00_RK                                              !
            c2 =  2.913180E+01_RK                                              !
            c3 = -3.494358E+01_RK                                              !
            c4 =  2.315955E+01_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            c0 =  1.674634E+00_RK                                              !
            c1 =  5.732017E+01_RK                                              !
            c2 =  5.955416E+01_RK                                              !
            c3 = -2.311007E+02_RK                                              !
            c4 =  1.255199E+02_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_m06_correlation                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine m06_correlation                                               !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine bmk_term(NGrid, ISpin, beta, c0, c1, c2, c3, c4, rho, gr2,        &
                                                        g, dg_drho, dg_dgr2)   !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in) :: NGrid                                           !
    integer(IK), intent(in) :: ISpin                                           !
    !                                                                          !
    real(RK),    intent(in) :: beta                                            !
    real(RK),    intent(in) :: c0                                              !
    real(RK),    intent(in) :: c1                                              !
    real(RK),    intent(in) :: c2                                              !
    real(RK),    intent(in) :: c3                                              !
    real(RK),    intent(in) :: c4                                              !
    !                                                                          !
    real(RK),    intent(in)  :: rho(NGrid, ISpin)                              !
    real(RK),    intent(in)  :: gr2(NGrid, ISpin)                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: g(NGrid)                                       !
    real(RK),    intent(out) :: dg_drho(NGrid, ISpin)                          !
    real(RK),    intent(out) :: dg_dgr2(NGrid, ISpin)                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)              :: s                                              !
    !                                                                          !
    real(RK)                 :: x2s(NGrid, ISpin)                              !
    real(RK)                 :: x2(NGrid)                                      !
    real(RK)                 :: dx2_drho(NGrid, ISpin)                         !
    real(RK)                 :: dx2_dgr2(NGrid, ISpin)                         !
    real(RK)                 :: v(NGrid)                                       !
    real(RK)                 :: dv_dx2(NGrid)                                  !
    real(RK)                 :: dg_dv(NGrid)                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    where (rho .ge. Dtol)                                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      x2s      = gr2 / (rho**(8.0_RK / 3.0_RK))                                !
      dx2_drho = - (8.0_RK / 3.0_RK) * gr2 / (rho**(11.0_RK / 3.0_RK))         !
      dx2_dgr2 = 1.0_RK / (rho**(8.0_RK / 3.0_RK))                             !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    elsewhere                                                                  !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      x2s      = 0.0_RK                                                        !
      dx2_drho = 0.0_RK                                                        !
      dx2_dgr2 = 0.0_RK                                                        !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    end where                                                                  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    x2 = sum(x2s, dim=2)                                                       !
    !                                                                          !
    v       = beta * x2 / (1.0_RK + beta * x2)                                 !
    dv_dx2  = beta / ((1.0_RK + beta * x2)**2)                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    g       = c0 + (c1 + (c2 + (c3 + c4 * v) * v) * v) * v                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    dg_dv   = c1 + (c2 * 2.0_RK + (c3 * 3.0_RK + c4 * 4.0_RK * v) * v) * v     !
    !                                                                          !
    DO s = 1, ISpin                                                            !
      !------------------------------------------------------------------------+
      !                                                                        !
      dg_drho(:,s) = dg_dv * dv_dx2 * dx2_drho(:,s)                            !
      dg_dgr2(:,s) = dg_dv * dv_dx2 * dx2_dgr2(:,s)                            !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO                                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine bmk_term                                                      !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module m06
! Default options for vim:sw=2:expandtab:smarttab:autoindent:syntax

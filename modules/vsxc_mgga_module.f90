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
module vsxc
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
  !                                                                            !
  use type_module,        only: IK => i4_kind, RK => r8_kind                   !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  public :: vsxc_X                                                             &
          , vsxc_C                                                             &
          , vsxc_exchange                                                      &
          , vsxc_correlation                                                   &
          , eX_lsda                                                            &
          , Dsic_term                                                          !
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
  subroutine vsxc_X(NGrid, ISpin, rho, gr2, tau, F, dF_drho, dF_dgr2, dF_dtau) !
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
      call vsxc_exchange( 0, NGrid, 2, rho_tmp, gr2_tmp, tau_tmp,              &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + dF_drho_tmp(:, 1) * 0.5_RK     !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + dF_dgr2_tmp(:, 1) * 0.25_RK    !
      dF_dtau(:NGrid, 1) = dF_dtau(:NGrid, 1) + dF_dtau_tmp(:, 1) * 0.5_RK     !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      call vsxc_exchange( 0, NGrid, 2                                          &
                        , rho(:NGrid, :2), gr2(:NGrid, :3), tau(:NGrid, :2),   &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      dF_drho(:NGrid, :2) = dF_drho(:NGrid, :2) + dF_drho_tmp                  !
      dF_dgr2(:NGrid, :3) = dF_dgr2(:NGrid, :3) + dF_dgr2_tmp                  !
      dF_dtau(:NGrid, :2) = dF_dtau(:NGrid, :2) + dF_dtau_tmp                  !
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
    F(:NGrid) = F(:NGrid) + F_tmp                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine vsxc_X                                                        !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine vsxc_C(NGrid, ISpin, rho, gr2, tau, F, dF_drho, dF_dgr2, dF_dtau) !
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
      call vsxc_correlation( 0, NGrid, 2, rho_tmp, gr2_tmp, tau_tmp,           &
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + dF_drho_tmp(:, 1) * 0.5_RK     !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + dF_dgr2_tmp(:, 1) * 0.25_RK    !
      dF_dtau(:NGrid, 1) = dF_dtau(:NGrid, 1) + dF_dtau_tmp(:, 1) * 0.5_RK     !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      call vsxc_correlation( 0, NGrid, 2                                       &
                           , rho(:NGrid, :2), gr2(:NGrid, :3), tau(:NGrid, :2),&
                                F_tmp, dF_drho_tmp, dF_dgr2_tmp, dF_dtau_tmp ) !
      !                                                                        !
      dF_drho(:NGrid, :2) = dF_drho(:NGrid, :2) + dF_drho_tmp                  !
      dF_dgr2(:NGrid, :3) = dF_dgr2(:NGrid, :3) + dF_dgr2_tmp                  !
      dF_dtau(:NGrid, :2) = dF_dtau(:NGrid, :2) + dF_dtau_tmp                  !
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
    F(:NGrid) = F(:NGrid) + F_tmp                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine vsxc_C                                                        !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine vsxc_exchange(PSet, NGrid, ISpin, rho, gr2, tau,                  &
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
    real(RK),    intent(out) :: F(NGrid)                                       !
    real(RK),    intent(out) :: dF_drho(:, :) ! (NGrid, ISpin)
    real(RK),    intent(out) :: dF_dgr2(:, :) ! (NGrid, 2*ISpin-1)
    real(RK),    intent(out) :: dF_dtau(:, :) ! (NGrid, ISpin)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)              :: s                                              !
    !                                                                          !
    real(RK)                 :: d0, d1, d2, d3, d4, d5                         !
    !                                                                          !
    real(RK)                 :: e(NGrid)                                       !
    real(RK)                 :: de_drho(NGrid)                                 !
    !                                                                          !
    real(RK)                 :: h(NGrid)                                       !
    real(RK)                 :: dh_drho(NGrid)                                 !
    real(RK)                 :: dh_dgr2(NGrid)                                 !
    real(RK)                 :: dh_dtau(NGrid)                                 !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: alpha = 0.00186726_RK                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F       = 0.0_RK                                                           !
    dF_drho = 0.0_RK                                                           !
    dF_dgr2 = 0.0_RK                                                           !
    dF_dtau = 0.0_RK                                                           !
    !                                                                          !
    call get_pars_vsxc_exchange(PSET, d0, d1, d2, d3, d4, d5)                  !
    !                                                                          !
    spin: DO s = 1, ISpin                                                      !
      !------------------------------------------------------------------------+
      !                                                                        !
      call eX_lsda(NGrid, rho(1:NGrid, s), e(1:NGrid), de_drho(1:NGrid))       !
      !                                                                        !
      call gvt4_term(NGrid, 1, alpha, d0, d1, d2, d3, d4, d5,                  &
                     rho(1:NGrid, s), gr2(1:NGrid, s), tau(1:NGrid, s),        &
                                   h(:), dh_drho(:), dh_dgr2(:), dh_dtau(:))   !
      !                                                                        !
      F(1:NGrid) = F(1:NGrid) + e(1:NGrid) * h(1:NGrid)                        !
      !                                                                        !
      dF_drho(1:NGrid, s) =  de_drho(1:NGrid) * h(1:NGrid)                     &
                          +  e(1:NGrid) * dh_drho(1:NGrid)                     !
      !                                                                        !
      dF_dgr2(1:NGrid, s) =  e(1:NGrid) * dh_dgr2(1:NGrid)                     !
      !                                                                        !
      dF_dtau(1:NGrid, s) =  e(1:NGrid) * dh_dtau(1:NGrid)                     !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO spin                                                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    subroutine get_pars_vsxc_exchange(PSET, d0, d1, d2, d3, d4, d5)            !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: d0                                           !
      real(RK),    intent(out) :: d1                                           !
      real(RK),    intent(out) :: d2                                           !
      real(RK),    intent(out) :: d3                                           !
      real(RK),    intent(out) :: d4                                           !
      real(RK),    intent(out) :: d5                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (PSet)                                                       !
        case(0)  ! This should be original VS98XC                              !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          d0 =  1.05324147599106640E+00_RK  !-9.800683E-01_RK / LDAfactor      !
          d1 =  3.82234242542822070E-03_RK  !-3.556788E-03_RK / LDAfactor      !
          d2 = -6.71698348131996321E-03_RK  ! 6.250326E-03_RK / LDAfactor      !
          d3 =  2.53030938105796706E-05_RK  !-2.354518E-05_RK / LDAfactor      !
          d4 =  1.37850244210630278E-04_RK  !-1.282732E-04_RK / LDAfactor      !
          d5 = -3.84172286736070984E-04_RK  ! 3.574822E-04_RK / LDAfactor      !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(1)  ! M06-L                                                       !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          d0 =  6.012244E-01_RK                                                !
          d1 =  4.748822E-03_RK                                                !
          d2 = -8.635108E-03_RK                                                !
          d3 = -9.308062E-06_RK                                                !
          d4 =  4.482811E-05_RK                                                !
          d5 =  0.000000E+00_RK                                                !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(2)  ! M06                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          d0 =  1.422057E-01_RK                                                !
          d1 =  7.370319E-04_RK                                                !
          d2 = -1.601373E-02_RK                                                !
          d3 =  0.000000E+00_RK                                                !
          d4 =  0.000000E+00_RK                                                !
          d5 =  0.000000E+00_RK                                                !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(3)  ! M06-2X all parameters zero                                  !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          d0 =  0.000000E+00_RK                                                !
          d1 =  0.000000E+00_RK                                                !
          d2 =  0.000000E+00_RK                                                !
          d3 =  0.000000E+00_RK                                                !
          d4 =  0.000000E+00_RK                                                !
          d5 =  0.000000E+00_RK                                                !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(4)  ! M06-HF                                                      !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          d0 = -1.179732E-01_RK                                                !
          d1 = -2.500000E-03_RK                                                !
          d2 = -1.180065E-02_RK                                                !
          d3 =  0.000000E+00_RK                                                !
          d4 =  0.000000E+00_RK                                                !
          d5 =  0.000000E+00_RK                                                !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_vsxc_exchange                                      !
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine vsxc_exchange                                                 !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine vsxc_correlation(PSet, NGrid, ISpin, rho, gr2, tau,               &
                                                F, dF_drho, dF_dgr2, dF_dtau)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    use pw_ldac_module,     only: pw_ldac                                      !
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
    real(RK)                 :: d0, d1, d2, d3, d4, d5                         !
    !                                                                          !
    real(RK)                 :: eo(NGrid)                                      !
    real(RK)                 :: deo_drho(NGrid, ISpin)                         !
    real(RK)                 :: rho_d_2(NGrid, 2)                              !
    real(RK)                 :: es(NGrid)                                      !
    real(RK)                 :: des_drho(NGrid, 2)                             !
    real(RK)                 :: h(NGrid)                                       !
    real(RK)                 :: dh_drho(NGrid)                                 !
    real(RK)                 :: dh_dgr2(NGrid)                                 !
    real(RK)                 :: dh_dtau(NGrid)                                 !
    real(RK)                 :: D(NGrid)                                       !
    real(RK)                 :: dD_drho(NGrid)                                 !
    real(RK)                 :: dD_dgr2(NGrid)                                 !
    real(RK)                 :: dD_dtau(NGrid)                                 !
    real(RK)                 :: dhos_drho(NGrid, ISpin)                        !
    real(RK)                 :: dhos_dgr2(NGrid, ISpin)                        !
    real(RK)                 :: dhos_dtau(NGrid, ISpin)                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: alpss = 0.00515088E0_RK                        !
    real(RK), parameter      :: alpos = 0.00304966E0_RK                        !
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
    ! Opposite spin, LDA part                                                  !
    !--------------------------------------------------------------------------+
    call pw_ldac(NGrid, ISpin, rho, eo, deo_drho)                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Same spin part                                                           !
    !--------------------------------------------------------------------------+
    call get_pars_vsxc_correlation(PSET, 1, d0, d1, d2, d3, d4, d5)            !
    !                                                                          !
    spin: DO s = 1, ISpin                                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      ! LDA part                                                               !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
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
      call gvt4_term(NGrid, 1, alpss, d0, d1, d2, d3, d4, d5,                  &
                     rho(1:NGrid, s), gr2(1:NGrid, s), tau(1:NGrid, s),        &
                                                 h, dh_drho, dh_dgr2, dh_dtau) !
      !                                                                        !
      call Dsic_term(NGrid, rho(1:NGrid, s), gr2(1:NGrid, s), tau(1:NGrid, s), &
                                                 D, dD_drho, dD_dgr2, dD_dtau) !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      ! Summarize same spin contributions                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      where ( rho(:,s) .ge. DTol)                                              !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
        !                                                                      !
        F                   = F + es * h * D                                   !
        !                                                                      !
        dF_drho(1:NGrid, s) = des_drho(1:NGrid, s) * h * D                     &
                            + es * dh_drho * D                                 &
                            + es * h * dD_drho                                 !
        !                                                                      !
        dF_dgr2(1:NGrid, s) = es * dh_dgr2 * D                                 &
                            + es * h * dD_dgr2                                 !
        !                                                                      !
        dF_dtau(1:NGrid, s) = es * dh_dtau * D                                 &
                            + es * h * dD_dtau                                 !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      end where                                                                !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO spin
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Opposite spin MGGA part                                                  !
    !--------------------------------------------------------------------------+
    call get_pars_vsxc_correlation(PSET, 2, d0, d1, d2, d3, d4, d5)            !
    !                                                                          !
    call gvt4_term(NGrid, ISpin, alpos, d0, d1, d2, d3, d4, d5,                &
                      rho, gr2, tau,     h, dhos_drho, dhos_dgr2, dhos_dtau)   !
    !                                                                          !
    where(sum(rho,dim=2) .ge. DTol)                                            !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      !                                                                        !
      F = F + eo * h                                                           !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
    end where                                                                  !
    !                                                                          !
    DO s = 1, ISpin                                                            !
      where ( rho(:,s) .ge. DTol)                                              !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
        !                                                                      !
        dF_drho(1:NGrid, s) = dF_drho(1:NGrid, s)                              &
                            + deo_drho(1:NGrid, s) * h                         &
                            + eo * dhos_drho(1:NGrid, s)                       !
        !                                                                      !
        dF_dgr2(1:NGrid, s) = dF_dgr2(1:NGrid, s)                              &
                            + eo * dhos_dgr2(1:NGrid, s)                       !
        !                                                                      !
        dF_dtau(1:NGrid, s) = dF_dtau(1:NGrid, s)                              &
                            + eo * dhos_dtau(1:NGrid, s)                       !
        !                                                                      !
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -+
      end where                                                                !
      !                                                                        !
    END DO                                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    subroutine get_pars_vsxc_correlation(Pset, term, d0, d1, d2, d3, d4, d5)   !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      integer(IK), intent(in)  :: term  ! 1 = same spin,  2 = opposite spin    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: d0                                           !
      real(RK),    intent(out) :: d1                                           !
      real(RK),    intent(out) :: d2                                           !
      real(RK),    intent(out) :: d3                                           !
      real(RK),    intent(out) :: d4                                           !
      real(RK),    intent(out) :: d5                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (PSet)                                                       !
        case(0)  ! This should be original VS98XC                              !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            d0 =   3.270912E-01_RK                                             !
            d1 =  -3.228915E-02_RK                                             !
            d2 =  -2.942406E-02_RK                                             !
            d3 =   2.134222E-03_RK                                             !
            d4 =  -5.451559E-03_RK                                             !
            d5 =   1.577575E-02_RK                                             !
          ELSE IF (term .eq. 2) THEN                                           !
            d0 =   7.035010E-01_RK                                             !
            d1 =   7.694574E-03_RK                                             !
            d2 =   5.152765E-02_RK                                             !
            d3 =   3.394308E-05_RK                                             !
            d4 =  -1.269420E-03_RK                                             !
            d5 =   1.296118E-03_RK                                             !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(1)  ! M06-L                                                       !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            d0 =  4.650534E-01_RK                                              !
            d1 =  1.617589E-01_RK                                              !
            d2 =  1.833657E-01_RK                                              !
            d3 =  4.692100E-04_RK                                              !
            d4 = -4.990573E-03_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            d0 =  3.957626E-01_RK                                              !
            d1 = -5.614546E-01_RK                                              !
            d2 =  1.403963E-02_RK                                              !
            d3 =  9.831442E-04_RK                                              !
            d4 = -3.577176E-03_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(2)  ! M06                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            d0 =  4.905945E-01_RK                                              !
            d1 = -1.437348E-01_RK                                              !
            d2 =  2.357824E-01_RK                                              !
            d3 =  1.871015E-03_RK                                              !
            d4 = -3.788963E-03_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            d0 = -2.741539E+00_RK                                              !
            d1 = -6.720113E-01_RK                                              !
            d2 = -7.932688E-02_RK                                              !
            d3 =  1.918681E-03_RK                                              !
            d4 = -2.032902E-03_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(3)  ! M06-2X                                                      !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            d0 =  6.902145E-01_RK                                              !
            d1 =  9.847204E-02_RK                                              !
            d2 =  2.214797E-01_RK                                              !
            d3 = -1.968264E-03_RK                                              !
            d4 = -6.775479E-03_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            d0 =  1.166404E-01_RK                                              !
            d1 = -9.120847E-02_RK                                              !
            d2 = -6.726189E-02_RK                                              !
            d3 =  6.720580E-05_RK                                              !
            d4 =  8.448011E-04_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case(4)  ! M06-HF                                                      !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          IF (term .eq. 1) THEN                                                !
            d0 =  8.976746E-01_RK                                              !
            d1 = -2.345830E-01_RK                                              !
            d2 =  2.368173E-01_RK                                              !
            d3 = -9.913890E-04_RK                                              !
            d4 = -1.146165E-02_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          ELSE IF (term .eq. 2) THEN                                           !
            d0 = -6.746338E-01_RK                                              !
            d1 = -1.534002E-01_RK                                              !
            d2 = -9.021521E-02_RK                                              !
            d3 = -1.292037E-03_RK                                              !
            d4 = -2.352983E-04_RK                                              !
            d5 =  0.000000E+00_RK                                              !
          END IF                                                               !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_vsxc_correlation                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine vsxc_correlation                                              !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine eX_lsda(NGrid, rho, e, de_drho)                                   !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in) :: NGrid                                           !
    !                                                                          !
    real(RK),    intent(in)  :: rho(NGrid)                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: e(NGrid)                                       !
    real(RK),    intent(out) :: de_drho(NGrid)                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: CX = -0.930525736349100_RK                     !
    !--------------------------------------------------------------------------+
    !                                                                          !
    where (rho .ge. Dtol)                                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      e       = CX * rho**(4.0_RK / 3.0_RK)                                    !
      de_drho = CX * (4.0_RK / 3.0_RK) * rho**(1.0_RK / 3.0_RK)                !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    elsewhere                                                                  !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      e       = 0.0_RK                                                         !
      de_drho = 0.0_RK                                                         !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    end where                                                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine eX_lsda                                                       !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine Dsic_term(NGrid, rho, gr2, tau,     D, dD_drho, dD_dgr2, dD_dtau) !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in) :: NGrid                                           !
    !                                                                          !
    real(RK),    intent(in)  :: rho(NGrid)                                     !
    real(RK),    intent(in)  :: gr2(NGrid)                                     !
    real(RK),    intent(in)  :: tau(NGrid)                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: D(NGrid)                                       !
    real(RK),    intent(out) :: dD_drho(NGrid)                                 !
    real(RK),    intent(out) :: dD_dgr2(NGrid)                                 !
    real(RK),    intent(out) :: dD_dtau(NGrid)                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK)                 :: x2(NGrid)                                      !
    real(RK)                 :: dx2_drho(NGrid)                                !
    real(RK)                 :: dx2_dgr2(NGrid)                                !
    real(RK)                 :: z(NGrid)                                       !
    real(RK)                 :: dz_drho(NGrid)                                 !
    real(RK)                 :: dz_dtau(NGrid)                                 !
    real(RK)                 :: dD_dx2(NGrid)                                  !
    real(RK)                 :: dD_dz(NGrid)                                   !
    !                                                                          !
    real(RK), parameter      :: CF = 9.11559974469119_RK ! 3/5(6pi^2)^(2/3)    !
    !--------------------------------------------------------------------------+
    !                                                                          !
    where (rho .ge. DTol .and. tau .ge. DTol)                                  !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      x2      = gr2 / (rho**(8.0_RK / 3.0_RK))                                 !
      dx2_drho = - (8.0_RK / 3.0_RK) * x2 / rho                                !
      dx2_dgr2 = 1.0_RK / (rho**(8.0_RK / 3.0_RK))                             !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      z       = 2.0_RK * tau / (rho**(5.0_RK / 3.0_RK)) - CF                   !
      dz_drho = - (5.0_RK / 3.0_RK) * 2.0_RK * tau / (rho**(8.0_RK / 3.0_RK))  !
      dz_dtau = 2.0_RK / (rho**(5.0_RK / 3.0_RK))                              !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      D       = 1.0_RK - x2 / (4.0_RK * (z + CF))                              !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      dD_dx2  = - 1.0_RK  / (4.0_RK * (z + CF))                                !
      dD_dz   = x2 * 4.0_RK / (16.0_RK * (z + CF) * (z + CF))                  !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      dD_drho = dD_dx2 * dx2_drho + dD_dz * dz_drho                            !
      dD_dgr2 = dD_dx2 * dx2_dgr2                                              !
      dD_dtau = dD_dz  * dz_dtau                                               !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    elsewhere                                                                  !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      D       = 0.0_RK                                                         !
      dD_drho = 0.0_RK                                                         !
      dD_dgr2 = 0.0_RK                                                         !
      dD_dtau = 0.0_RK                                                         !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    end where                                                                  !
    !                                                                          !
    where (D .lt. 0.0_RK)                                                      !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      D       = 0.0_RK                                                         !
      dD_drho = 0.0_RK                                                         !
      dD_dgr2 = 0.0_RK                                                         !
      dD_dtau = 0.0_RK                                                         !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    end where                                                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine Dsic_term                                                     !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine gvt4_term(NGrid, ISpin,                                           &
                          alpha, d0, d1, d2, d3, d4, d5, rho, gr2, tau,        &
                                               h, dh_drho, dh_dgr2, dh_dtau)   !
    !--------------------------------------------------------------------------+
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in) :: NGrid                                           !
    integer(IK), intent(in) :: ISpin                                           !
    !                                                                          !
    real(RK),    intent(in) :: alpha                                           !
    real(RK),    intent(in) :: d0                                              !
    real(RK),    intent(in) :: d1                                              !
    real(RK),    intent(in) :: d2                                              !
    real(RK),    intent(in) :: d3                                              !
    real(RK),    intent(in) :: d4                                              !
    real(RK),    intent(in) :: d5                                              !
    !                                                                          !
    real(RK),    intent(in)  :: rho(NGrid, ISpin)                              !
    real(RK),    intent(in)  :: gr2(NGrid, ISpin)                              !
    real(RK),    intent(in)  :: tau(NGrid, ISpin)                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: h(NGrid)                                       !
    real(RK),    intent(out) :: dh_drho(NGrid, ISpin)                          !
    real(RK),    intent(out) :: dh_dgr2(NGrid, ISpin)                          !
    real(RK),    intent(out) :: dh_dtau(NGrid, ISpin)                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK)              :: s                                              !
    !                                                                          !
    real(RK)                 :: x2s(NGrid, ISpin)                              !
    real(RK)                 :: x2(NGrid)                                      !
    real(RK)                 :: dx2_drho(NGrid, ISpin)                         !
    real(RK)                 :: dx2_dgr2(NGrid, ISpin)                         !
    real(RK)                 :: zs(NGrid, ISpin)                               !
    real(RK)                 :: z(NGrid)                                       !
    real(RK)                 :: dz_drho(NGrid, ISpin)                          !
    real(RK)                 :: dz_dtau(NGrid, ISpin)                          !
    real(RK)                 :: g(NGrid)                                       !
    real(RK)                 :: dg_dx2(NGrid)                                  !
    real(RK)                 :: dg_dz(NGrid)                                   !
    real(RK)                 :: dh_dx2(NGrid)                                  !
    real(RK)                 :: dh_dz(NGrid)                                   !
    real(RK)                 :: dh_dg(NGrid)                                   !
    !                                                                          !
    real(RK), parameter      :: CF = 9.11559974469119_RK ! 3/5(6pi^2)^(2/3)    !
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
      !                                                                        !
      zs      = 2.0_RK * tau / (rho**(5.0_RK / 3.0_RK)) - CF                   !
      dz_drho = -2.0_RK * (5.0_RK / 3.0_RK) * tau / (rho**(8.0_RK / 3.0_RK))   !
      dz_dtau = 2.0_RK / (rho**(5.0_RK / 3.0_RK))                              !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    elsewhere                                                                  !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !                                                                        !
      x2s      = 0.0_RK                                                        !
      dx2_drho = 0.0_RK                                                        !
      dx2_dgr2 = 0.0_RK                                                        !
      !                                                                        !
      zs       = 0.0_RK                                                        !
      dz_drho  = 0.0_RK                                                        !
      dz_dtau  = 0.0_RK                                                        !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    end where                                                                  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    x2 = sum(x2s, dim=2)                                                       !
    z  = sum(zs,  dim=2)                                                       !
    !                                                                          !
    g       = 1.0_RK + alpha * (x2 + z)                                        !
    dg_dx2  = alpha                                                            !
    dg_dz   = alpha                                                            !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    h       = d0 / g                                                           &
            + (d1 * x2 + d2 * z) / (g * g)                                     &
            + (d3 * x2 * x2 + d4 * x2 * z + d5 * z * z) / (g * g * g)          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    dh_dx2  = d1 / (g * g)                                                     &
            + (d3 * 2.0_RK * x2 + d4 * z) / (g * g * g)                        !
    dh_dz   = d2 / (g * g)                                                     &
            + (d4 * x2 + d5 * 2.0_RK * z) / (g * g * g)                        !
    dh_dg   = - d0 / (g * g)                                                   &
            - 2.0_RK * (d1 * x2 + d2 * z) / (g * g * g)                        &
            - 3.0_RK * (d3 * x2 * x2 + d4 * x2 * z + d5 * z * z) / (g**4)      !
    !                                                                          !
    DO s = 1, ISpin                                                            !
      !------------------------------------------------------------------------+
      !                                                                        !
      dh_drho(:,s) = dh_dx2 * dx2_drho(:,s)                                    &
                   + dh_dz  * dz_drho(:,s)                                     &
                   + dh_dg  * dg_dx2  * dx2_drho(:,s)                          &
                   + dh_dg  * dg_dz   * dz_drho(:,s)                           !
      !                                                                        !
      dh_dgr2(:,s) = dh_dx2 * dx2_dgr2(:,s)                                    &
                   + dh_dg  * dg_dx2  * dx2_dgr2(:,s)                          !
      !                                                                        !
      dh_dtau(:,s) = dh_dz  * dz_dtau(:,s)                                     &
                   + dh_dg  * dg_dz   * dz_dtau(:,s)                           !
      !                                                                        !
      !------------------------------------------------------------------------+
    END DO                                                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine gvt4_term                                                     !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module vsxc
! Default options for vim:sw=2:expandtab:smarttab:autoindent:syntax

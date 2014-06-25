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
module lyp
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
  !  References: B. Miehlich, A. Savin, H. Stoll and H. Preuss;                !
  !              Chem. Phys. Lett. 157 (1989), p. 200                          !
  !                                                                            !
  !              NB: the original formulas in Phys. Rev. B 37 (1988), p. 785   !
  !                  are commonly referenced, but almost never used            !
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
  use ad1x5, only: ad &
                              , var                                            &
                              , fix                                            &
                              , val                                            &
                              , fst                                            &
                              , operator(+)                                    &
                              , operator(-)                                    &
                              , operator(*)                                    &
                              , operator(/)                                    &
                              , operator(**)                                   &
                              , operator(>) &
                              , exp                                            !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  public :: lyp_C                                                              !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  real(RK), parameter      :: DTol  = 1.0E-30_RK                               !
  real(RK), parameter      :: n1d3  = 1.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n2d3  = 2.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n4d3  = 4.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n8d3  = 8.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n11d3 = 11.0_RK / 3.0_RK                         !
  real(RK), parameter      :: n1d9  = 1.0_RK / 9.0_RK                          !
  real(RK), parameter      :: n7d9  = 7.0_RK / 9.0_RK                          !
  real(RK), parameter      :: n47d9 = 47.0_RK / 9.0_RK                         !
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
  subroutine lyp_C(PSet, NGrid, ISpin, rho, gr2,         F, dF_drho, dF_dgr2)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! PSet denotes parameter set:           1  -- LYP-98                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! Performs conversion to UKS case to autodiff variables                    !
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
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),  intent(inout) :: F(:)                                           !
    real(RK),  intent(inout) :: dF_drho(:,:)                                   !
    real(RK),  intent(inout) :: dF_dgr2(:,:)                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: r_a(NGrid) , r_b(NGrid)                        !
    type(ad)                 :: g_aa(NGrid), g_bb(NGrid), g_ab(NGrid)          !
    type(ad)                 :: F_ad(NGrid)                                    !
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
      !                                                                        !
      F_ad = lyp_corr( PSet, r_a, r_b, g_aa, g_bb, g_ab )                      !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + fst( 1, F_ad )                 !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + fst( 3, F_ad )                 !
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
      !                                                                        !
      F_ad = lyp_corr( PSet, r_a, r_b, g_aa, g_bb, g_ab )                      !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + fst( 1, F_ad )                 !
      dF_drho(:NGrid, 2) = dF_drho(:NGrid, 2) + fst( 2, F_ad )                 !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + fst( 3, F_ad )                 !
      dF_dgr2(:NGrid, 2) = dF_dgr2(:NGrid, 2) + fst( 4, F_ad )                 !
      dF_dgr2(:NGrid, 3) = dF_dgr2(:NGrid, 3) + fst( 5, F_ad )                 !
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
    F               = F + val( F_ad )                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine lyp_C                                                         !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  elemental function lyp_corr(PSet, r_a, r_b, g_aa, g_bb, g_ab) result(F)      !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: PSet                                           !
    !                                                                          !
    type(ad),    intent(in)  :: r_a , r_b                                      !
    type(ad),    intent(in)  :: g_aa, g_bb, g_ab                               !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: F                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK)                 :: a_par, b_par, c_par, d_par, e_par              !
    !                                                                          !
    type(ad)                 :: rsum, gsum                                     !
    type(ad)                 :: qZ, qw, qd                                     !
    type(ad)                 :: lda_t1, lda_t2                                 !
    type(ad)                 :: gga_aa, gga_ab, gga_bb                         !
    type(ad)                 :: f1, f2                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    rsum = r_a  + r_b                                                          !
    gsum = g_aa + g_bb + 2.0_RK * g_ab                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Get parameters of the Colle Salvetti formula                             !
    !--------------------------------------------------------------------------+
    call get_pars_lyp_correlation(PSet, a_par, b_par, c_par, d_par, e_par)     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Calculate some intermediates                                             !
    !--------------------------------------------------------------------------+

    IF ( rsum > DTol ) THEN                                               !
      qZ = rsum**n1d3 / ( d_par + rsum**n1d3 )                                 !
      qw = exp( -c_par / (rsum**n1d3) ) * qZ * rsum**(-n11d3)                  !
      qd = ( c_par + d_par * qZ ) / ( rsum**n1d3 )                             !
    ELSE                                                                       !
      qZ = fix(0.0_RK) ! not used in this case
      qw = fix(0.0_RK)                                                         !
      qd = fix(0.0_RK)                                                         !
    ENDIF                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! Some prefactors                                                          !
    !--------------------------------------------------------------------------+
    f1 = a_par * b_par * qw                                                    !
    f2 = f1 * r_a * r_b                                                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! LDA part                                                                 !
    !--------------------------------------------------------------------------+
    IF ( rsum > DTol ) THEN                                               !
      lda_t1 = 4.0_RK * a_par * r_a * r_b * qZ / rsum                          !
    ELSE                                                                       !
      lda_t1 = fix( 0.0_RK )                                                   !
    ENDIF                                                                      !
    !                                                                          !
    lda_t2 = f2 * 2.0_RK**n11d3 * e_par * (r_a**n8d3 + r_b**n8d3)              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    ! GGA part                                                                 !
    !--------------------------------------------------------------------------+
    IF ( rsum > DTol ) THEN                                               !
      gga_aa = ( f2 * ( n1d9 - n1d3*qd - n1d9*(qd - 11.0_RK)*r_a / rsum )      &
                 - f1 * r_b**2 ) * g_aa                                        !
    ELSE                                                                       !
      gga_aa = ( f2 * ( n1d9 - n1d3*qd ) - f1 * r_b**2 ) * g_aa                !
    ENDIF                                                                      !
    !                                                                          !
    IF ( rsum > DTol ) THEN                                               !
      gga_bb = ( f2 * ( n1d9 - n1d3*qd - n1d9*(qd - 11.0_RK)*r_b / rsum )      &
               - f1 * r_a**2 ) * g_bb                                          !
    ELSE                                                                       !
      gga_bb = ( f2 * ( n1d9 - n1d3*qd ) - f1 * r_a**2 ) * g_bb                !
    ENDIF                                                                      !
    !                                                                          !
    gga_ab   = (f2 * ( n47d9 - n7d9 * qd ) - f1 * n4d3 * rsum**2 ) * g_ab      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F  = - lda_t1 - lda_t2 - gga_aa - gga_bb - gga_ab                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    pure subroutine get_pars_lyp_correlation(PSET, a, b, c, d, e)              !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: a                                            !
      real(RK),    intent(out) :: b                                            !
      real(RK),    intent(out) :: c                                            !
      real(RK),    intent(out) :: d                                            !
      real(RK),    intent(out) :: e                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (PSet)                                                       !
        case(1)  ! LYP                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          !                                                                    !
          a = 0.04918_RK                                                       !
          b = 0.13200_RK                                                       !
          c = 0.25330_RK                                                       !
          d = 0.34900_RK                                                       !
          e = 2.87123400018819_RK ! = 3/10*(3*pi*pi)^(2/3)                     !
          !                                                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_lyp_correlation                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end function lyp_corr                                                        !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module lyp
!# Default options for vim:sw=2:expandtab:smarttab:autoindent:syntax

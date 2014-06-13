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
module pbe_gga_module
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
  !  Purpose: Autodiff implementation of the PBE functional                    !
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
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  public :: pbe_X                                                              !
  public :: pbe_C                                                              !
  public :: pbe_correlation                                                    !
  public :: pbe_correlation_ad                                                 !
  public :: eX_lsda                                                            !
  public :: eX_lsda_ad                                                         !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  integer(IK), parameter, public :: PBE_key         = 1                        &
                                  , PBE_revTPSS_key = 2                        !
  !                                                                            !
  real(RK), parameter      :: DTol = 1.0E-30_RK                                !
  !                                                                            !
  real(RK), parameter      :: pi   = 3.14159265358979323846_RK                 !
  real(RK), parameter      :: crs  = 0.620350490899400_RK                      !
  real(RK), parameter      :: ck2  = 3.938980087370790_RK                      !
  !                                                                            !
  real(RK), parameter      :: n1d3 = 1.0_RK / 3.0_RK                           !
  real(RK), parameter      :: n2d3 = 2.0_RK / 3.0_RK                           !
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
  subroutine pbe_X( PSet, NGrid, ISpin, rho, gr2,        F, dF_drho, dF_dgr2 ) !
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
    type(ad)                 :: g_aa(NGrid), g_bb(NGrid)                       !
    type(ad)                 :: F_ad(NGrid)                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    IF ( PSet .ne. 1 ) THEN                                                    !
      write(*,*) 'PSet /= 1'                                                   !
      stop                                                                     !
    ENDIF                                                                      !
    !                                                                          !
    IF (ISpin .eq. 1) THEN                                                     !
      !------------------------------------------------------------------------+
      !                                                                        !
      r_a           = var( 1, rho(:NGrid, 1) * 0.50_RK )                       !
      r_b           = var( 2, rho(:NGrid, 1) * 0.50_RK )                       !
      g_aa          = var( 3, gr2(:NGrid, 1) * 0.25_RK )                       !
      g_bb          = var( 4, gr2(:NGrid, 1) * 0.25_RK )                       !
      !                                                                        !
      F_ad = pbe_exchange( r_a, r_b, g_aa, g_bb )                              !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + fst( 1, F_ad )                 !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + fst( 3, F_ad ) * 0.50_RK       !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      r_a           = var( 1, rho(:NGrid, 1) )                                 !
      r_b           = var( 2, rho(:NGrid, 2) )                                 !
      g_aa          = var( 3, gr2(:NGrid, 1) )                                 !
      g_bb          = var( 4, gr2(:NGrid, 2) )                                 !
      !                                                                        !
      F_ad = pbe_exchange( r_a, r_b, g_aa, g_bb )                              !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + fst( 1, F_ad )                 !
      dF_drho(:NGrid, 2) = dF_drho(:NGrid, 2) + fst( 2, F_ad )                 !
      dF_dgr2(:NGrid, 1) = dF_dgr2(:NGrid, 1) + fst( 3, F_ad )                 !
      dF_dgr2(:NGrid, 2) = dF_dgr2(:NGrid, 2) + fst( 4, F_ad )                 !
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
    F(:NGrid)     = F(:NGrid) + val( F_ad )                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine pbe_X                                                         !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine pbe_C( PSet, NGrid, ISpin, rho, gr2,       F, dF_drho, dF_dgr2 )  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! PSet denotes parameter set:         1  -- PBE                            !
    !                                     2  -- modified PBE for revTPSS       !
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
      !                                                                        !
      DO ind = 1, NGrid                                                        !
        F_ad(ind) = pbe_correlation_ad( PSet                                   &
                                 , r_a(ind), r_b(ind)                          &
                                 , g_aa(ind), g_bb(ind), g_ab(ind) )           &
                  * ( r_a(ind) + r_b(ind) )                                    !
      ENDDO                                                                    !
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
      !                                                                        !
      DO ind = 1, NGrid                                                        !
        F_ad(ind) = pbe_correlation_ad( PSet                                   &
                                 , r_a(ind), r_b(ind)                          &
                                 , g_aa(ind), g_bb(ind), g_ab(ind) )           &
                  * ( r_a(ind) + r_b(ind) )                                    !
      ENDDO                                                                    !
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
    F(:NGrid) = F(:NGrid) + val( F_ad )                                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine pbe_C                                                         !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  elemental function pbe_exchange( r_a, r_b, g_aa, g_bb ) result(F)            !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    type(ad),    intent(in)  :: r_a , r_b                                      !
    type(ad),    intent(in)  :: g_aa, g_bb                                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: F                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: ea, eb                                         !
    type(ad)                 :: pa, pb                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ea = eX_lsda_ad( r_a )                                                     !
    eb = eX_lsda_ad( r_b )                                                     !
    !                                                                          !
    pa = pbeef( r_a, g_aa )                                                    !
    pb = pbeef( r_b, g_bb )                                                    !
    !                                                                          !
    F = ea * pa + eb * pb                                                      !
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function pbeef( rho, gr2 ) result(p)                             !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in)  :: rho                                          !
      type(ad),    intent(in)  :: gr2                                          !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                 :: p                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                 :: x2                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK), parameter      :: pbec2 = 4.49267E-3_RK                        !
      real(RK), parameter      :: pbec1 = 3.612108584108E-03_RK                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      IF ( val(rho) .ge. Dtol ) THEN                                           !
        !                                                                      !
        x2 = gr2 / (rho**(8.0_RK / 3.0_RK))                                    !
        !                                                                      !
      ELSE                                                                     !
        !                                                                      !
        x2 = fix(0.0_RK)                                                       !
        !                                                                      !
      ENDIF                                                                    !
      !                                                                        !
      p    = 1.0_RK + pbec1 * x2 / ( 1.0_RK + pbec2 * x2 )                     !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function pbeef                                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end function pbe_exchange                                                    !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  elemental function eX_lsda_ad( rho ) result(e)                               !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad),    intent(in)  :: rho                                            !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: e                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK), parameter      :: CX = -0.930525736349100_RK                     !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    IF (val(rho) .ge. Dtol ) THEN                                              !
      !                                                                        !
      e = CX * rho**(4.0_RK / 3.0_RK)                                          !
      !                                                                        !
    ELSE                                                                       !
      !                                                                        !
      e = fix( 0.0_RK )                                                        !
      !                                                                        !
    ENDIF                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end function eX_lsda_ad                                                      !
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
    real(RK),    intent(in)  :: rho(:) ! (NGrid)
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: e(:) ! (NGrid)
    real(RK),    intent(out) :: de_drho(:) ! (NGrid)
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
  elemental function pbe_correlation_ad( PSet                                  &
                                       , r_a, r_b, g_aa, g_bb, g_ab ) result(F)!
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! PSet denotes parameter set:         1  -- PBE                            !
    !                                     2  -- modified PBE for revTPSS       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    use type_module,        only: IK => i4_kind, RK => r8_kind                 !
    use pw_lda,             only: pw_lda_correlation                           &
                                , std_Key                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
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
    type(ad)                 :: rsum, gsum                                     !
    type(ad)                 :: qr                                             !
    type(ad)                 :: qp                                             !
    type(ad)                 :: qt2                                            !
    type(ad)                 :: qb                                             !
    !                                                                          !
    type(ad)                 :: L                                              !
    type(ad)                 :: H                                              !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! Some constants                                                           !
    real(RK), parameter      :: g_par = 0.031090690869654895_RK                !
    real(RK), parameter      :: b_par = 0.066725_RK                            !
    ! set beta to match with old one... reference:  beta = 0.06672455060314922 !
    real(RK), parameter      :: c_par = 0.1000_RK                              !
    real(RK), parameter      :: d_par = 0.1778_RK                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    rsum   = r_a  + r_b                                                        !
    gsum   = g_aa + g_bb + 2.0_RK * g_ab                                       !
    !                                                                          !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! INTERMEDIATE QUANTITIES                                                  !
    !                                                                          !
    ! r = (3 / (4 * pi * rho))^(1/3)                                           !
    qr     = seitz_radius( rsum )                                              !
    !                                                                          !
    ! p = ((1 + z)^(2/3) + (1 - z)^(2/3)) / 2                                  !
    qp     = spin_scaling_factor( r_a, r_b, rsum )                             !
    !                                                                          !
    ! t^2 = gr2 / (2 * p * k(rho) * rho)                                       !
    qt2    = dimless_gradient( rsum, gsum, qp )                                !
    !                                                                          !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! LDA CORRELATION ENERGY DENSITY PER PARTICLE                              !
    !                                                                          !
    L      = pw_lda_correlation( std_Key, r_a, r_b )                           !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! GGA CORRECTION TERM H                                                    !
    !                                                                          !
    ! beta is a parameter in PBE and TPSS but a function in revTPSS            !
    IF (PSet == PBE_key) THEN                                                  !
      qb   = fix(b_par)                                                        !
    ELSEIF (PSet == PBE_revTPSS_key) THEN                                      !
      qb   = d_par * (1.0_RK + c_par * qr) / (1.0_RK + d_par * qr)             !
    ENDIF                                                                      !
    !                                                                          !
    H      = gga_correction( L, qp, qt2, qb )                                  !
    !                                                                          !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! ADDING EVERYTHING TO PBE CORRELATION ENERGY DENSITY                      !
    !                                                                          !
    F      = L + H                                                             !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    elemental function seitz_radius( rt ) result(r)                            !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: rt                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: r                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      r = crs * ( rt + DTol )**(-n1d3)                                         !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function seitz_radius                                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function spin_scaling_factor( ra, rb, rt ) result(p)             !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: ra, rb, rt                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: p                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      IF( val(ra) > DTol .and. val(rb) > DTol ) THEN                           !
        p = 0.5_RK * ( (2.0_RK * ra / rt)**n2d3                                &
                     + (2.0_RK * rb / rt)**n2d3 )                              !
      ELSEIF ( val(ra) > DTol ) THEN                                           !
        p = 0.5_RK *   (2.0_RK * ra / rt)**n2d3                                !
      ELSEIF ( val(rb) > DTol ) THEN                                           !
        p = 0.5_RK *   (2.0_RK * rb / rt)**n2d3                                !
      ELSE                                                                     !
        p = fix(0.0_RK)                                                        !
      ENDIF                                                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function spin_scaling_factor                                           !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function screening_wave_number( rt ) result(k2)                  !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: rt                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: k2                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      k2 = ck2 * (rt + DTol )**n1d3                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function screening_wave_number                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function dimless_gradient( rt, gt, p ) result(t2)                !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: rt                                            !
      type(ad),    intent(in) :: gt                                            !
      type(ad),    intent(in) :: p                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: t2                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: qk2                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      qk2 = screening_wave_number( rt )                                        !
      !                                                                        !
      t2  = gt / ( qk2 * (2.0_RK * p * rt)**2 + DTol )                         !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function dimless_gradient                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function gga_correction( L, p, t2, b ) result(H)                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: L                                             !
      type(ad),    intent(in) :: p                                             !
      type(ad),    intent(in) :: t2                                            !
      type(ad),    intent(in) :: b                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: H                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: qA                                            !
      type(ad)                :: qI                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      qA = A_term( L, p, b )                                                   !
      !                                                                        !
      qI = (t2 + qA * t2**2) / (1.0_RK + qA * t2 + qA**2 * t2**2)              !
      !                                                                        !
      H  = g_par * p**3 * log(1.0_RK + b / g_par * qI )                        !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function gga_correction                                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function A_term( L, p, b ) result(A)                             !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: L                                             !
      type(ad),    intent(in) :: p                                             !
      type(ad),    intent(in) :: b                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: A                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: D                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      D = g_par * (exp(- L / (g_par * p**3 + DTol)) - 1.0_RK)                  !
      !                                                                        !
      A = b / ( D + sign(DTol, val(D)) )                                       !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function A_term                                                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end function pbe_correlation_ad                                              !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  subroutine pbe_correlation(Flavor, half_calc, NGrid, rho, gr2,               &
                                 F, dF_drho, dF_dgr2,                          &
                                 ddF_drho_drho, ddF_drho_dgr2, ddF_dgr2_dgr2)  !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! Flavor denotes parameter set:         1  -- PBE                          !
    !                                       2  -- modified PBE for revTPSS     !
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
    integer(IK), intent(in)  :: Flavor                                         !
    logical,     intent(in)  :: half_calc                                      !
    integer(IK), intent(in)  :: NGrid                                          !
    !                                                                          !
    real(RK),    intent(in)  :: rho(NGrid, 2)                                  !
    real(RK),    intent(in)  :: gr2(NGrid, 3)                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),    intent(out) :: F(NGrid)                                       !
    real(RK),    intent(out) :: dF_drho(NGrid, 2)                              !
    real(RK),    intent(out) :: dF_dgr2(NGrid, 3)                              !
    real(RK), optional, intent(out) :: ddF_drho_drho(NGrid, 3)                 !
    real(RK), optional, intent(out) :: ddF_drho_dgr2(NGrid, 6)                 !
    real(RK), optional, intent(out) :: ddF_dgr2_dgr2(NGrid, 6)                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK)                 :: rsum(NGrid) !                                  !
    !                                                                          !
    real(RK)                 :: r(NGrid)    ! seitz radius                     !
    real(RK)                 :: dr_drho(NGrid, 2)                              !
    !                                                                          !
    real(RK)                 :: z(NGrid)    ! relative spin polarization       !
    real(RK)                 :: dz_drho(NGrid, 2)                              !
    !                                                                          !
    real(RK)                 :: p(NGrid)    ! spin scaling factor              !
    real(RK)                 :: dp_dz(NGrid)                                   !
    !                                                                          !
    real(RK)                 :: t2(NGrid)    ! dimensionless density gradient  !
    real(RK)                 :: dt2_dgr2(NGrid,3)                              !
    real(RK)                 :: dt2_drho(NGrid)                                !
    real(RK)                 :: dt2_dp(NGrid)                                  !
    !                                                                          !
    real(RK)                 :: beta(NGrid)  ! constant in PBE                 !
    real(RK)                 :: dbeta_dr(NGrid)                                !
    !                                                                          !
    real(RK)                 :: L(NGrid)     ! LDA-correlation kernel          !
    real(RK)                 :: dL_drho(NGrid,2)                               !
    !                                                                          !
    real(RK)                 :: H(NGrid)     ! H-term                          !
    real(RK)                 :: dH_dL(NGrid)                                   !
    real(RK)                 :: dH_dp(NGrid)                                   !
    real(RK)                 :: dH_dt2(NGrid)                                  !
    real(RK)                 :: dH_dbeta(NGrid)                                !
    !                                                                          !
    real(RK), parameter      :: gam_par = 0.031090690869654895_RK              !
    real(RK), parameter      :: betapar = 0.066725_RK                          !
    ! set beta to match with old one... reference:  beta = 0.06672455060314922 !
    real(RK), parameter      :: a_par = 0.1000_RK                              !
    real(RK), parameter      :: b_par = 0.1778_RK                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F       = 0.0_RK                                                           !
    dF_drho = 0.0_RK                                                           !
    dF_dgr2 = 0.0_RK                                                           !
    !                                                                          !
    !                                                                          !
    rsum    = sum(rho,2)                                                       !
    !                                                                          !
    ! r = (3 / (4 * pi * rho))^(1/3)                                           !
    call seitz_radius(NGrid, rho, r, dr_drho)                                  !
    !                                                                          !
    ! z = (rho_a - rho_b) / (rho_a + rho_b)                                    !
    call relative_spin_polarization(NGrid, rho, z, dz_drho)                    !
    !                                                                          !
    ! p = ((1 + z)^(2/3) + (1 - z)^(2/3)) / 2                                  !
    call spin_scaling_factor(NGrid, z, p, dp_dz)                               !
    !                                                                          !
    ! t^2 = gr2 / (2 * p * k(rho) * rho)                                       !
    call dimless_density_gradient(NGrid, gr2, rsum, p,                         &
                                               t2, dt2_dgr2, dt2_drho, dt2_dp) !
    !                                                                          !
    !                                                                          !
    ! L = e_C(LDA) / rho                                                       !
    L = 0.0_RK                                                                 !
    dL_drho = 0.0_RK                                                           !
    call pw_ldac(NGrid, 2, rho, L, dL_drho)                                    !
    !                                                                          !
    if (half_calc) then                                                        !
      where(rho .lt. DTol)                                                     !
        dL_drho = 0.0_RK                                                       !
      endwhere                                                                 !
    endif                                                                      !
    !                                                                          !
    where (sum(rho,2) .ge. DTol)                                               !
      !                                                                        !
      L            = L / sum(rho,2)                                            !
      dL_drho(:,1) = (dL_drho(:,1) - L) / sum(rho,2)                           !
      dL_drho(:,2) = (dL_drho(:,2) - L) / sum(rho,2)                           !
      !                                                                        !
    elsewhere                                                                  !
      !                                                                        !
      L            = 0.0_RK                                                    !
      dL_drho(:,1) = 0.0_RK                                                    !
      dL_drho(:,2) = 0.0_RK                                                    !
      !                                                                        !
    end where                                                                  !
    !                                                                          !
    !                                                                          !
    if (Flavor ==  PBE_key) then                                               !
      !                                                                        !
      beta = betapar                                                           !
      dbeta_dr = 0.0_RK                                                        !
      !                                                                        !
    elseif (Flavor == PBE_revTPSS_key) then                                    !
      !                                                                        !
      beta = betapar * (1.0_RK + a_par * r) / (1.0_RK + b_par * r)             !
      dbeta_dr = betapar * (a_par/ (1.0_RK + b_par * r)                        &
               - beta * b_par / (1.0_RK + b_par * r))                          !
      !                                                                        !
    endif                                                                      !
    !                                                                          !
    call H_term(NGrid, L, p, t2, beta,       H, dH_dL, dH_dp, dH_dt2, dH_dbeta)!
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    F = sum(rho,2) * (L + H)                                                   !
    !                                                                          !
    dF_drho(:,1) = (L + H) + sum(rho,2) * (dL_drho(:,1)                        &
                 + dH_dL * dL_drho(:,1)                                        &
                 + dH_dp * dp_dz * dz_drho(:,1)                                &
                 + dH_dt2 * (dt2_drho + dt2_dp * dp_dz * dz_drho(:,1))         &
                 + dH_dbeta * dbeta_dr * dr_drho(:,1))                         !
    dF_drho(:,2) = (L + H) + sum(rho,2) * (dL_drho(:,2)                        &
                 + dH_dL * dL_drho(:,2)                                        &
                 + dH_dp * dp_dz * dz_drho(:,2)                                &
                 + dH_dt2 * (dt2_drho + dt2_dp * dp_dz * dz_drho(:,2))         &
                 + dH_dbeta * dbeta_dr * dr_drho(:,2))                         !
    dF_dgr2(:,1) = sum(rho,2) * dH_dt2 * dt2_dgr2(:,1)                         !
    dF_dgr2(:,2) = sum(rho,2) * dH_dt2 * dt2_dgr2(:,2)                         !
    dF_dgr2(:,3) = sum(rho,2) * dH_dt2 * dt2_dgr2(:,3)                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    if (present(ddF_drho_drho)) then                                           !
      ! for now                                                                !
      ddF_drho_drho = 0.0_RK                                                   !
      ddF_drho_dgr2 = 0.0_RK                                                   !
      ddF_dgr2_dgr2 = 0.0_RK                                                   !
      !                                                                        !
    end if                                                                     !
    !                                                                          !
    if (half_calc) then                                                        !
      where(rho .lt. DTol)                                                     !
        dF_drho = 0.0_RK                                                       !
      endwhere                                                                 !
      where(gr2 .lt. DTol)                                                     !
        dF_dgr2 = 0.0_RK                                                       !
      endwhere                                                                 !
    endif                                                                      !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    subroutine seitz_radius(NGrid, rho, r, dr_drho, ddr_drho_drho)             !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: r(NGrid)                                     !
      real(RK),    intent(out) :: dr_drho(NGrid, 2)                            !
      real(RK), optional, intent(out) :: ddr_drho_drho(NGrid)                  !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (sum(rho, 2) .ge. DTol)                                            !
        !                                                                      !
        r = (3.0_RK / (4.0_RK * pi * sum(rho, 2)))**(1.0_RK/3.0_RK)            !
        !                                                                      !
        dr_drho(:, 1) = - (1.0_RK / 3.0_RK) * r / sum(rho, 2)                  !
        dr_drho(:, 2) = dr_drho(:, 1)                                          !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        r = (3.0_RK / (4.0_RK * pi * DTol))**(1.0_RK/3.0_RK)                   !
        !                                                                      !
        dr_drho(:, 1) = - (1.0_RK / 3.0_RK) * r / DTol                         !
        dr_drho(:, 2) = dr_drho(:, 1)                                          !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !                                                                        !
      if (present(ddr_drho_drho)) then                                         !
        !                                                                      !
        ddr_drho_drho = 0.0_RK                                                 !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine seitz_radius                                                !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine relative_spin_polarization(NGrid, rho,              z, dz_drho, &
                                                                ddz_drho_drho) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid, 2)                                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: z(NGrid)                                     !
      real(RK),    intent(out) :: dz_drho(NGrid, 2)                            !
      real(RK), optional, intent(out) :: ddz_drho_drho(NGrid, 3)               !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (sum(rho, 2) .ge. DTol)                                            !
        !                                                                      !
        z = (rho(:,1) - rho(:,2)) / sum(rho, 2)                                !
        !                                                                      !
        dz_drho(:, 1) =   2.0_RK * rho(:,2) / ((sum(rho, 2))**2)               !
        dz_drho(:, 2) = - 2.0_RK * rho(:,1) / ((sum(rho, 2))**2)               !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        z = 0.0_RK                                                             !
        !                                                                      !
        dz_drho(:, 1) = 0.0_RK                                                 !
        dz_drho(:, 2) = 0.0_RK                                                 !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !                                                                        !
      if (present(ddz_drho_drho)) then                                         !
        !                                                                      !
        where (sum(rho, 2) .ge. DTol)                                          !
          ! (:, 1) = AA,  (:, 2) = BB,  (:, 3) = AB                            !
          ddz_drho_drho(:, 1) = - 2.0_RK * dz_drho(:, 1) / sum(rho, 2)         !
          ddz_drho_drho(:, 2) = - 2.0_RK * dz_drho(:, 2) / sum(rho, 2)         !
          ddz_drho_drho(:, 3) = 2.0_RK * z / ((sum(rho, 2))**2)                !
          !                                                                    !
        elsewhere                                                              !
          !                                                                    !
          ddz_drho_drho(:, 1) = 0.0_RK                                         !
          ddz_drho_drho(:, 2) = 0.0_RK                                         !
          ddz_drho_drho(:, 3) = 0.0_RK                                         !
          !                                                                    !
        end where                                                              !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine relative_spin_polarization                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine spin_scaling_factor(NGrid, z, p, dp_dz, ddp_dz_dz)              !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: z(NGrid)                                     !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),           intent(out) :: p(NGrid)                              !
      real(RK),           intent(out) :: dp_dz(NGrid)                          !
      real(RK), optional, intent(out) :: ddp_dz_dz(NGrid)                      !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      p = 0.5_RK * ((1.0_RK + z)**(2.0_RK/3.0_RK)                              &
                  + (1.0_RK - z)**(2.0_RK/3.0_RK))                             !
      !                                                                        !
      WHERE( 1.0_RK - abs(z) > DTol )                                          !
        dp_dz = (1.0_RK / 3.0_RK)                                              &
              * ((1.0_RK / (1.0_RK + z)**(1.0_RK/3.0_RK))                      &
               - (1.0_RK / (1.0_RK - z)**(1.0_RK/3.0_RK)))                     !
      ELSEWHERE                                                                !
        dp_dz = 0.0_RK                                                         !
      ENDWHERE                                                                 !
      !                                                                        !
      !                                                                        !
      if (present(ddp_dz_dz)) then                                             !
        !                                                                      !
        ddp_dz_dz = (1.0_RK / 9.0_RK)                                          &
                                    * ((1.0_RK + z + DTol)**(-4.0_RK/3.0_RK)   &
                                     + (1.0_RK - z + DTol)**(-4.0_RK/3.0_RK))  !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine spin_scaling_factor                                         !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine screening_wave_number(NGrid, rho, k2, dk2_drho, ddk2_drho_drho) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),           intent(out) :: k2(NGrid)                             !
      real(RK),           intent(out) :: dk2_drho(NGrid)                       !
      real(RK), optional, intent(out) :: ddk2_drho_drho(NGrid)                 !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      k2 = 4.0_RK * (rho * 3.0_RK / pi)**(1.0_RK / 3.0_RK)                     !
      !                                                                        !
      where (rho .ge. DTol)                                                    !
        !                                                                      !
        dk2_drho = (1.0_RK / 3.0_RK) * k2 / rho                                !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        dk2_drho = 0.0_RK                                                      !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      if (present(ddk2_drho_drho)) then                                        !
        !                                                                      !
        where (rho .ge. DTol)                                                  !
          !                                                                    !
          ddk2_drho_drho = - (2.0_RK / 3.0_RK) * dk2_drho / rho                !
          !                                                                    !
        elsewhere                                                              !
          !                                                                    !
          ddk2_drho_drho = 0.0_RK                                              !
          !                                                                    !
        end where                                                              !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine screening_wave_number                                       !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine dimless_density_gradient(NGrid, gr2, rho, p                     &
                                       ,t2, dt2_dgr2, dt2_drho, dt2_dp         &
                                       ,ddt2_dgr2_dgr2, ddt2_dgr2_drho         &
                                       ,ddt2_dgr2_dp, ddt2_drho_drho           &
                                       ,ddt2_drho_dp, ddt2_dp_dp)              !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: gr2(NGrid, 3)                                !
      real(RK),    intent(in)  :: rho(NGrid)                                   !
      real(RK),    intent(in)  :: p(NGrid)                                     !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),           intent(out) :: t2(NGrid)                             !
      real(RK),           intent(out) :: dt2_dgr2(NGrid, 3)                    !
      real(RK),           intent(out) :: dt2_drho(NGrid)                       !
      real(RK),           intent(out) :: dt2_dp(NGrid)                         !
      real(RK), optional, intent(out) :: ddt2_dgr2_dgr2(NGrid, 6)              !
      real(RK), optional, intent(out) :: ddt2_dgr2_drho(NGrid, 3)              !
      real(RK), optional, intent(out) :: ddt2_dgr2_dp(NGrid, 3)                !
      real(RK), optional, intent(out) :: ddt2_drho_drho(NGrid)                 !
      real(RK), optional, intent(out) :: ddt2_drho_dp(NGrid)                   !
      real(RK), optional, intent(out) :: ddt2_dp_dp(NGrid)                     !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                        :: k2(NGrid)                             !
      real(RK)                        :: dk2_drho(NGrid)                       !
      real(RK)                        :: ddk2_drho_drho(NGrid)                 !
      real(RK)                        :: dt2_dk2(NGrid)                        !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      if (.not. present(ddt2_dgr2_dgr2)) then                                  !
        !                                                                      !
        call screening_wave_number(NGrid, rho, k2, dk2_drho)                   !
        !                                                                      !
      else ! present(ddt2_dgr2_dgr2)                                           !
        !                                                                      !
        call screening_wave_number(NGrid, rho, k2, dk2_drho, ddk2_drho_drho)   !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !                                                                        !
      where (rho .ge. DTol)                                                    !
        !                                                                      !
        t2 = (gr2(:,1) + gr2(:,2) + 2.0_RK * gr2(:,3))                         &
                                                / (k2 * (2.0_RK * p * rho)**2) !
        !                                                                      !
        dt2_dk2 = - t2 / k2                                                    !
        !                                                                      !
        dt2_drho = - 2.0_RK * t2 / rho + dt2_dk2 * dk2_drho                    !
        !                                                                      !
        dt2_dgr2(:,1) = 1.0_RK / (k2 * (2.0_RK * p * rho)**2)                  !
        dt2_dgr2(:,2) = dt2_dgr2(:,1)                                          !
        dt2_dgr2(:,3) = dt2_dgr2(:,2) * 2.0_RK                                 !
        !                                                                      !
        dt2_dp   = - 2.0_RK * t2 / p                                           !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        t2 = 0.0_RK                                                            !
        !                                                                      !
        dt2_dgr2(:,1) = 0.0_RK                                                 !
        dt2_dgr2(:,2) = 0.0_RK                                                 !
        dt2_dgr2(:,3) = 0.0_RK                                                 !
        !                                                                      !
        dt2_drho = 0.0_RK                                                      !
        dt2_dp   = 0.0_RK                                                      !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      if (present(ddt2_drho_drho)) then                                        !
        ! for now                                                              !
        ddt2_dgr2_dgr2 = 0.0_RK                                                !
        ddt2_dgr2_drho = 0.0_RK                                                !
        ddt2_dgr2_dp   = 0.0_RK                                                !
        ddt2_drho_drho = 0.0_RK                                                !
        ddt2_drho_dp   = 0.0_RK                                                !
        ddt2_dp_dp     = 0.0_RK                                                !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine dimless_density_gradient                                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine H_term(NGrid, L, p, t2, beta, H, dH_dL, dH_dp, dH_dt2, dH_dbeta,&
                                              ddH_dL_dL, ddH_dL_dp, ddH_dL_dt2,&
                                           ddH_dL_dbeta, ddH_dp_dp, ddH_dp_dt2,&
                                      ddH_dp_dbeta, ddH_dt2_dt2, ddH_dt2_dbeta,&
                                                               ddH_dbeta_dbeta)!
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: L(NGrid)                                     !
      real(RK),    intent(in)  :: p(NGrid)                                     !
      real(RK),    intent(in)  :: t2(NGrid)                                    !
      real(RK),    intent(in)  :: beta(NGrid)                                  !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),           intent(out) :: H(NGrid)                              !
      real(RK),           intent(out) :: dH_dL(NGrid)                          !
      real(RK),           intent(out) :: dH_dp(NGrid)                          !
      real(RK),           intent(out) :: dH_dt2(NGrid)                         !
      real(RK),           intent(out) :: dH_dbeta(NGrid)                       !
      real(RK), optional, intent(out) :: ddH_dL_dL(NGrid)                      !
      real(RK), optional, intent(out) :: ddH_dL_dp(NGrid)                      !
      real(RK), optional, intent(out) :: ddH_dL_dt2(NGrid)                     !
      real(RK), optional, intent(out) :: ddH_dL_dbeta(NGrid)                   !
      real(RK), optional, intent(out) :: ddH_dp_dp(NGrid)                      !
      real(RK), optional, intent(out) :: ddH_dp_dt2(NGrid)                     !
      real(RK), optional, intent(out) :: ddH_dp_dbeta(NGrid)                   !
      real(RK), optional, intent(out) :: ddH_dt2_dt2(NGrid)                    !
      real(RK), optional, intent(out) :: ddH_dt2_dbeta(NGrid)                  !
      real(RK), optional, intent(out) :: ddH_dbeta_dbeta(NGrid)                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                        :: A(NGrid)                              !
      real(RK)                        :: dA_dL(NGrid)                          !
      real(RK)                        :: dA_dp(NGrid)                          !
      real(RK)                        :: dA_dbeta(NGrid)                       !
      real(RK)                        :: ddA_dL_dL(NGrid)                      !
      real(RK)                        :: ddA_dL_dp(NGrid)                      !
      real(RK)                        :: ddA_dL_dbeta(NGrid)                   !
      real(RK)                        :: ddA_dp_dp(NGrid)                      !
      real(RK)                        :: ddA_dp_dbeta(NGrid)                   !
      real(RK)                        :: ddA_dbeta_dbeta(NGrid)                !
      real(RK)                        :: one_d_At2A2t4(NGrid)                  !
      real(RK)                        :: arg(NGrid)                            !
      real(RK)                        :: darg_dA(NGrid)                        !
      real(RK)                        :: darg_dt2(NGrid)                       !
      real(RK)                        :: darg_dbeta(NGrid)                     !
      real(RK)                        :: ddarg_dA_dA(NGrid)                    !
      real(RK)                        :: ddarg_dA_dt2(NGrid)                   !
      real(RK)                        :: ddarg_dA_dbeta(NGrid)                 !
      real(RK)                        :: ddarg_dt2_dt2(NGrid)                  !
      real(RK)                        :: ddarg_dt2_dbeta(NGrid)                !
      real(RK)                        :: ddarg_dbeta_dbeta(NGrid)              !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      if (.not. present(ddH_dL_dL)) then                                       !
        !                                                                      !
        call A_term(NGrid, L, p, beta,              A, dA_dL, dA_dp, dA_dbeta) !
        !                                                                      !
      else ! present(ddH_dL_dL)                                                !
        !                                                                      !
        write(*,*) '2nd derivatives not checked -> stop'                       !
        stop                                                                   !
      ! call A_term(NGrid, L, p, beta,              A, dA_dL, dA_dp, dA_dbeta, &
      !                                    ddA_dL_dL, ddA_dL_dp, ddA_dL_dbeta, &
      !                              ddA_dp_dp, ddA_dp_dbeta, ddA_dbeta_dbeta) !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      one_d_At2A2t4 = 1.0_RK / (1.0_RK + A * t2 + A**2 * t2**2)                !
      !                                                                        !
      ! argument of logarithm                                                  !
      arg = 1.0_RK + (beta / gam_par) * (t2 + A * t2**2) * one_d_At2A2t4       !
      !                                                                        !
      darg_dA = beta / gam_par * (t2**2 * one_d_At2A2t4                        &
              - (t2 + A * t2**2) * (t2 + 2.0_RK * A * t2**2) * one_d_At2A2t4**2)
      darg_dt2 = beta / gam_par * ((1.0_RK + 2.0_RK * A * t2) * one_d_At2A2t4  &
              - (t2 + A * t2**2) * (A + 2.0_RK * A**2 * t2) * one_d_At2A2t4**2)!
      darg_dbeta = (1.0_RK / gam_par) * (t2 + A * t2**2) * one_d_At2A2t4       !
      !                                                                        !
      !                                                                        !
      if (.not. present(ddH_dL_dL)) then                                       !
        !                                                                      !
        ddarg_dA_dA = beta / gam_par * one_d_At2A2t4**2                        &
                    * (- t2 * (t2 + 2.0_RK * A * t2**2)                        &
                    - t2**2 * (t2 + 2.0_RK * A * t2**2)                        &
                    - (t2 + A * t2**2) * 2.0_RK * t2**2)                       !
        ddarg_dA_dt2 = beta / gam_par * (2.0_RK * t2 * one_d_At2A2t4           &
                     - t2**2 * one_d_At2A2t4**2 * (A + 2.0_RK * A**2 * t2)     &
                     + 2.0_RK * (t2 + A * t2**2) * (t2 + 2.0_RK * A * t2**2)   &
                     * one_d_At2A2t4**3 * (A + 2.0_RK * A**2 * t2)             &
                     - one_d_At2A2t4**2 * ((1.0_RK + 2.0_RK * A * t2)          &
                     * (t2 + 2.0_RK * A * t2**2)                               &
                     + (t2 + A * t2**2) * (1.0_RK + 4.0_RK * A * t2)))         !
        ddarg_dA_dbeta = 1.0_RK / gam_par * (t2**2 * one_d_At2A2t4             &
                       - (t2 + A * t2**2) * (t2 + 2.0_RK * A * t2**2)          &
                       * one_d_At2A2t4**2)                                     !
        ddarg_dt2_dt2 = beta / gam_par * 2.0_RK * (A * one_d_At2A2t4           &
                      - ((1.0_RK + 2.0_RK * A * t2) * (A + 2.0_RK * A**2 * t2) &
                      + A**2) * one_d_At2A2t4**2                               &
                      + (t2 + A * t2**2) * (A + 2.0_RK * A**2 * t2)**2         &
                      * one_d_At2A2t4**2)                                      !
        ddarg_dt2_dbeta = 1.0_RK / gam_par * ((1.0_RK + 2.0_RK * A * t2)       &
                        * one_d_At2A2t4                                        &
                        - (t2 + A * t2**2) * (A + 2.0_RK * A**2 * t2)          &
                        * one_d_At2A2t4**2)                                    !
        ddarg_dbeta_dbeta = 0.0_RK                                             !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !                                                                        !
      H = gam_par * p**3 * log(arg)                                            !
      !                                                                        !
      dH_dL = gam_par * p**3 / arg * darg_dA * dA_dL                           !
      dH_dp = 3.0_RK * gam_par * p**2 * log(arg)                               &
            + gam_par * p**3 / arg * darg_dA * dA_dp                           !
      dH_dt2 = gam_par * p**3 / arg * darg_dt2                                 !
      dH_dbeta = gam_par * p**3 / arg * (darg_dbeta + darg_dA * dA_dbeta)      !
      !                                                                        !
      if (present(ddH_dL_dL)) then                                             !
        !                                                                      !
        ddH_dL_dL = - gam_par * p**3 / arg**2 * (darg_dA * dA_dL)**2           &
                  + gam_par * p**3 / arg * ddarg_dA_dA * dA_dL **2             &
                  + gam_par * p**3 / arg * darg_dA * ddA_dL_dL                 !
        ddH_dL_dp = 3.0_RK * gam_par * p**2 / arg * darg_dA * dA_dL            &
                  + gam_par * p**3 / arg**2 * darg_dA**2 * dA_dL * dA_dp       &
                  + gam_par * p**3 / arg * ddarg_dA_dA * dA_dL * dA_dp         &
                  + gam_par * p**3 / arg * darg_dA * ddA_dL_dp                 !
        ddH_dL_dt2 = 3.0_RK * gam_par * p**2 / arg * darg_dA * dA_dL           &
                  - gam_par * p**3 / arg**2 * darg_dA * dA_dL * darg_dt2       !
        ddH_dL_dbeta = -gam_par * p**3 / arg**2 * darg_dA**2 * dA_dL * dA_dbeta&
                  + gam_par * p**3 / arg * ddarg_dA_dbeta * dA_dL              &
                  + gam_par * p**3 / arg * ddarg_dA_dA * dA_dL * dA_dbeta      &
                  + gam_par * p**3 / arg * darg_dA * ddA_dL_dbeta              !
        ddH_dp_dp = 6.0_RK * gam_par * p * log(arg)                            &
                  + gam_par * p**3 / arg**2 * (darg_dA * dA_dp)**2             &
                  + gam_par * p**3 / arg * ddarg_dA_dA * (dA_dp)**2            &
                  + gam_par * p**3 / arg * darg_dA * ddA_dp_dp                 !
        ddH_dp_dt2 = 3.0_RK * gam_par * p**2 / arg * darg_dA * dA_dp           &
                  - gam_par * p**3 / arg**2 * darg_dt2 * darg_dA * dA_dp       !
        ddH_dp_dbeta = 3.0_RK * gam_par * p**2 / arg * darg_dbeta              &
                  + 3.0_RK * gam_par * p**2 / arg * darg_dA * dA_dbeta         &
                  - gam_par * p**3 / arg**2 * darg_dA * dA_dp * darg_dbeta     &
                  - gam_par * p**3 / arg**2 * darg_dA**2 * dA_dp * dA_dbeta    &
                  + gam_par * p**3 / arg * ddarg_dA_dbeta * dA_dp              &
                  + gam_par * p**3 / arg * ddarg_dA_dA * dA_dp * dA_dbeta      &
                  + gam_par * p**3 / arg * darg_dA * ddA_dp_dbeta              !
        ddH_dt2_dt2 = - gam_par * p**3 / arg**2 * darg_dt2**2                  &
                  + gam_par * p**3 / arg * ddarg_dt2_dt2                       !
        ddH_dt2_dbeta = - gam_par * p**3 / arg**2 * darg_dt2 * darg_dbeta      &
                  - gam_par * p**3 / arg**2 * darg_dt2 * darg_dA * dA_dbeta    &
                  + gam_par * p**3 / arg * ddarg_dt2_dbeta                     &
                  + gam_par * p**3 / arg * ddarg_dA_dt2 * dA_dbeta             !
        ddH_dbeta_dbeta = - gam_par * p**3 / arg**2 * (darg_dbeta**2           &
                  + darg_dbeta * darg_dA * dA_dbeta)                           &
                  + gam_par * p**3 / arg * ddarg_dbeta_dbeta                   &
                  + gam_par * p**3 / arg * ddarg_dA_dbeta * dA_dbeta           &
                  - gam_par * p**3 / arg**2 * (darg_dA * dA_dbeta * darg_dbeta &
                  + darg_dA**2 * dA_dbeta**2)                                  &
                  + gam_par * p**3 / arg * (ddarg_dA_dbeta * dA_dbeta          &
                  + ddarg_dA_dA * dA_dbeta**2 + darg_dA * ddA_dbeta_dbeta)     !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine H_term                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    subroutine A_term(NGrid, L, p, beta,            A, dA_dL, dA_dp, dA_dbeta, &
                                           ddA_dL_dL, ddA_dL_dp, ddA_dL_dbeta, &
                                     ddA_dp_dp, ddA_dp_dbeta, ddA_dbeta_dbeta) !
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: NGrid                                        !
      !                                                                        !
      real(RK),    intent(in)  :: L(NGrid)                                     !
      real(RK),    intent(in)  :: p(NGrid)                                     !
      real(RK),    intent(in)  :: beta(NGrid)                                  !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),           intent(out) :: A(NGrid)                              !
      real(RK),           intent(out) :: dA_dL(NGrid)                          !
      real(RK),           intent(out) :: dA_dp(NGrid)                          !
      real(RK),           intent(out) :: dA_dbeta(NGrid)                       !
      real(RK), optional, intent(out) :: ddA_dL_dL(NGrid)                      !
      real(RK), optional, intent(out) :: ddA_dL_dp(NGrid)                      !
      real(RK), optional, intent(out) :: ddA_dL_dbeta(NGrid)                   !
      real(RK), optional, intent(out) :: ddA_dp_dp(NGrid)                      !
      real(RK), optional, intent(out) :: ddA_dp_dbeta(NGrid)                   !
      real(RK), optional, intent(out) :: ddA_dbeta_dbeta(NGrid)                !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK)                        :: denom(NGrid)                          !
      real(RK)                        :: ddenom_dL(NGrid)                      !
      real(RK)                        :: ddenom_dp(NGrid)                      !
      real(RK)                        :: dddenom_dL_dL(NGrid)                  !
      real(RK)                        :: dddenom_dL_dp(NGrid)                  !
      real(RK)                        :: dddenom_dp_dp(NGrid)                  !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      where (p .ge. DTol)                                                      !
        !                                                                      !
        denom = exp(- L / (gam_par * p**3)) - 1.0_RK                           !
        !                                                                      !
        ddenom_dL = - exp(- L / (gam_par * p**3)) / (gam_par * p**3)           !
        ddenom_dp = 3.0_RK * L / (gam_par * p**4) * exp(- L / (gam_par * p**3))!
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        denom = - 1.0_RK                                                       !
        ddenom_dL = 0.0_RK                                                     !
        ddenom_dp = 0.0_RK                                                     !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      if (present(ddA_dL_dL)) then                                             !
        !                                                                      !
        where (p .ge. DTol)                                                    !
          !                                                                    !
          dddenom_dL_dL = exp(- L / (gam_par * p**3)) / ((gam_par * p**3)**2)  !
          dddenom_dL_dp = - ddenom_dL * 3.0_RK * L / (gam_par * p**4)          &
                          - ddenom_dL * 3.0_RK / p                             !
          dddenom_dp_dp = - ddenom_dp * (4.0_RK / p                            &
                                              + 3.0_RK * L / (gam_par * p**3)) !
          !                                                                    !
        elsewhere                                                              !
          !                                                                    !
          dddenom_dL_dL = 0.0_RK                                               !
          dddenom_dL_dp = 0.0_RK                                               !
          dddenom_dp_dp = 0.0_RK                                               !
          !                                                                    !
        end where                                                              !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !                                                                        !
      where (denom .ge. DTol)                                                  !
        !                                                                      !
        A = beta / (gam_par * denom)                                           !
        !                                                                      !
        dA_dL = - A / denom * ddenom_dL                                        !
        dA_dp = - A / denom * ddenom_dp                                        !
        dA_dbeta = 1.0_RK / (gam_par * denom)                                  !
        !                                                                      !
      elsewhere                                                                !
        !                                                                      !
        A = 0.0_RK                                                             !
        dA_dL = 0.0_RK                                                         !
        dA_dp = 0.0_RK                                                         !
        dA_dbeta = 0.0_RK                                                      !
        !                                                                      !
      end where                                                                !
      !                                                                        !
      !                                                                        !
      if (present(ddA_dL_dL)) then                                             !
        !                                                                      !
        where (denom .ge. DTol)                                                !
          !                                                                    !
          ddA_dL_dL = 2.0_RK * A / (denom**2) * ddenom_dL**2                   &
                                                   - A / denom * dddenom_dL_dL !
          ddA_dL_dp = 2.0_RK * A / (denom**2) * ddenom_dL * ddenom_dp          &
                                                   - A / denom * dddenom_dL_dp !
          ddA_dL_dbeta = - 1.0_RK / (gam_par * denom) * ddenom_dL              !
          ddA_dp_dp = 2.0_RK * A / (denom**2) * ddenom_dp**2                   &
                                                   - A / denom * dddenom_dp_dp !
          ddA_dp_dbeta = - 1.0_RK / (gam_par * denom) * ddenom_dp              !
          ddA_dbeta_dbeta = 0.0_RK                                             !
          !                                                                    !
        elsewhere                                                              !
          !                                                                    !
          ddA_dL_dL = 0.0_RK                                                   !
          ddA_dL_dp = 0.0_RK                                                   !
          ddA_dL_dbeta = 0.0_RK                                                !
          ddA_dp_dp = 0.0_RK                                                   !
          ddA_dp_dbeta = 0.0_RK                                                !
          ddA_dbeta_dbeta = 0.0_RK                                             !
          !                                                                    !
        end where                                                              !
        !                                                                      !
      end if                                                                   !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine A_term                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end subroutine pbe_correlation                                               !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module pbe_gga_module
! Default options for vim:sw=2:expandtab:smarttab:autoindent:syntax

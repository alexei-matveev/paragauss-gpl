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
module pw_lda
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
  !  Purpose: Automatic derivative version of the Perdew Wang LDA functional   !
  !                                                                            !
  !                                                                            !
  !  Module called by: PBE automatic derivative version                        !
  !                                                                            !
  !                                                                            !
  !  References: J.P. Perdew, Y. Wang;                                         !
  !              Phys. Rev. B 45 (1992), p. 13244                              !
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
                              , log                                            !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  implicit none                                                                !
  private         ! by default, all names are private                          !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  public :: pw_lda_C                                                           !
  public :: pw_lda_correlation                                                 !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  integer(IK), parameter, public :: std_Key  = 0                               !
  integer(IK), parameter, public :: RPA_Key  = 3                               !
  !                                                                            !
  integer(IK), parameter         :: ECU_key  = 1                               &
                                  , ECP_key  = 2                               &
                                  , mALF_key = 3                               !
  !                                                                            !
  real(RK), parameter      :: DTol  = 1.0E-30_RK                               !
  !                                                                            !
  real(RK), parameter      :: n1d3  = 1.0_RK / 3.0_RK                          !
  real(RK), parameter      :: n4d3  = 4.0_RK / 3.0_RK                          !
  !                                                                            !
  real(RK), parameter      :: crs   = 0.620350490899400_RK                     !
  real(RK), parameter      :: cf    = 1.923661050931540_RK                     !
  real(RK), parameter      :: cfzz  = 1.709920934161370_RK                     !
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
  subroutine pw_lda_C(PSet, NGrid, ISpin, rho,       F, dF_drho, ddF_drhodrho) !
    !--------------------------------------------------------------------------+
    !                                                                          !
    ! PSet denotes parameter set:  std_key (0) -- standard PW91 LDA            !
    !                              RPA_key (3) -- RPA                          !
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
    integer(IK),            intent(in) :: PSet                                 !
    integer(IK),            intent(in) :: NGrid                                !
    integer(IK),            intent(in) :: ISpin                                !
    !                                                                          !
    real(RK),               intent(in) :: rho(:,:)                             !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK),            intent(inout) :: F(:)                                 !
    real(RK),            intent(inout) :: dF_drho(:,:)                         !
    real(RK), optional,  intent(inout) :: ddF_drhodrho(:,:)                    !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                           :: r_a(NGrid) , r_b(NGrid)              !
    type(ad)                           :: F_ad(NGrid)                          !
    !                                                                          !
    integer(IK)                        :: ind                                  !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    IF (ISpin .eq. 1) THEN                                                     !
      !------------------------------------------------------------------------+
      !                                                                        !
      r_a           = var( 1, rho(:NGrid, 1) * 0.50_RK )                       !
      r_b           = var( 2, rho(:NGrid, 1) * 0.50_RK )                       !
      !                                                                        !
      DO ind = 1, NGrid                                                        !
        F_ad(ind) = pw_lda_correlation( PSet, r_a(ind), r_b(ind) )             &
                  * ( r_a(ind) + r_b(ind) )                                    !
      ENDDO                                                                    !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + fst( 1, F_ad )                 !
      !                                                                        !
      IF ( present(ddF_drhodrho) ) THEN                                        !
        ! FIXME                                                                !
      ENDIF                                                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
    ELSEIF (ISpin .eq. 2) THEN                                                 !
      !------------------------------------------------------------------------+
      !                                                                        !
      r_a           = var( 1, rho(:NGrid, 1) )                                 !
      r_b           = var( 2, rho(:NGrid, 2) )                                 !
      !                                                                        !
      DO ind = 1, NGrid                                                        !
        F_ad(ind) = pw_lda_correlation( PSet, r_a(ind), r_b(ind) )             &
                  * ( r_a(ind) + r_b(ind) )                                    !
      ENDDO                                                                    !
      !                                                                        !
      dF_drho(:NGrid, 1) = dF_drho(:NGrid, 1) + fst( 1, F_ad )                 !
      dF_drho(:NGrid, 2) = dF_drho(:NGrid, 2) + fst( 2, F_ad )                 !
      !                                                                        !
      IF ( present(ddF_drhodrho) ) THEN                                        !
        ! FIXME                                                                !
      ENDIF                                                                    !
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
  end subroutine pw_lda_C                                                      !
  !                                                                            !
  !----------------------------------------------------------------------------+
  !                                                                            !
  elemental function pw_lda_correlation(PSet, r_a, r_b) result(F)              !
    !--------------------------------------------------------------------------+
    !                                                                          !
    implicit none                                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    integer(IK), intent(in)  :: PSet                                           !
    !                                                                          !
    type(ad),    intent(in)  :: r_a , r_b                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    type(ad)                 :: F                                              !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    real(RK)                 :: a_par, a1par                                   !
    real(RK)                 :: b1par, b2par, b3par, b4par                     !
    real(RK)                 :: p_par                                          !
    !                                                                          !
    type(ad)                 :: rsum                                           !
    type(ad)                 :: rs                                             !
    type(ad)                 :: qz4                                            !
    type(ad)                 :: qf                                             !
    type(ad)                 :: e_c_r_0                                        !
    type(ad)                 :: e_c_r_1                                        !
    type(ad)                 :: a_c_r                                          !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    rsum      = r_a + r_b                                                      !
    !                                                                          !
    IF ( val(rsum) < DTol ) THEN                                               !
      !                                                                        !
      F       = fix(0.0_RK)                                                    !
      !                                                                        !
    ELSE                                                                       !
      !                                                                        !
      rs      = seitz_radius( rsum )                                           !
      !                                                                        !
      qz4     = relative_spin_pol( r_a, r_b, rsum )                            !
      !                                                                        !
      qf      = spin_interpolation( r_a, r_b, rsum )                           !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      ! unpolarized limit:  e_C( r_s, 0 )                                      !
      !                                                                        !
      call get_pars_pw_lda_correlation( ECU_key+PSET, a_par, a1par             &
                                      , b1par, b2par, b3par, b4par, p_par )    !
      !                                                                        !
      e_c_r_0 = G_form( a_par, a1par, b1par, b2par, b3par, b4par, p_par, rs )  !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      ! fully polarized limit:  e_C( r_s, 1 )                                  !
      !                                                                        !
      call get_pars_pw_lda_correlation( ECP_key+PSET, a_par, a1par             &
                                      , b1par, b2par, b3par, b4par, p_par )    !
      !                                                                        !
      e_c_r_1 = G_form( a_par, a1par, b1par, b2par, b3par, b4par, p_par, rs )  !
      !                                                                        !
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      ! spin stiffness:  a_C( r_s )                                            !
      !                                                                        !
      call get_pars_pw_lda_correlation( mALF_key+PSET, a_par, a1par            &
                                      , b1par, b2par, b3par, b4par, p_par )    !
      !                                                                        !
      a_c_r   = G_form( a_par, a1par, b1par, b2par, b3par, b4par, p_par, rs )  !
      !                                                                        !
      !                                                                        !
      F       =    e_c_r_0                                                     &
                + (e_c_r_1 - e_c_r_0) * qf * qz4                               &
                -  a_c_r * qf * (1.0_RK - qz4) / cfzz                          !
      !                                                                        !
    ENDIF                                                                      !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    contains                                                                   !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    pure subroutine get_pars_pw_lda_correlation(PSET, a, a1, b1, b2, b3, b4, p)!
      !------------------------------------------------------------------------+
      !                                                                        !
      integer(IK), intent(in)  :: PSet                                         !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(out) :: a                                            !
      real(RK),    intent(out) :: a1                                           !
      real(RK),    intent(out) :: b1                                           !
      real(RK),    intent(out) :: b2                                           !
      real(RK),    intent(out) :: b3                                           !
      real(RK),    intent(out) :: b4                                           !
      real(RK),    intent(out) :: p                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      select case (PSet)                                                       !
        case( ECU_key+std_key )  ! ECU                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          a  = 0.031091_RK                                                     !
          a1 = 0.21370_RK                                                      !
          b1 = 7.5957_RK                                                       !
          b2 = 3.5876_RK                                                       !
          b3 = 1.6382_RK                                                       !
          b4 = 0.49294_RK                                                      !
          p  = 1.0_RK                                                          !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case( ECP_key+std_key )  ! ECP                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          a  = 0.015545_RK                                                     !
          a1 = 0.20548_RK                                                      !
          b1 = 14.1189_RK                                                      !
          b2 = 6.1977_RK                                                       !
          b3 = 3.3662_RK                                                       !
          b4 = 0.62517_RK                                                      !
          p  = 1.0_RK                                                          !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case( mALF_key+std_key )  ! mALF                                       !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          a  = 0.016887_RK                                                     !
          a1 = 0.11125_RK                                                      !
          b1 = 10.357_RK                                                       !
          b2 = 3.6231_RK                                                       !
          b3 = 0.88026_RK                                                      !
          b4 = 0.49671_RK                                                      !
          p  = 1.0_RK                                                          !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case( ECU_key+RPA_key )  ! ECU, RPA                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          a  = 0.031091_RK                                                     !
          a1 = 0.082477_RK                                                     !
          b1 = 5.1486_RK                                                       !
          b2 = 1.6483_RK                                                       !
          b3 = 0.23647_RK                                                      !
          b4 = 0.20614_RK                                                      !
          p  = 0.75_RK                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case( ECP_key+RPA_key )  ! ECP, RPA                                    !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          a  = 0.015545_RK                                                     !
          a1 = 0.035374_RK                                                     !
          b1 = 6.4869_RK                                                       !
          b2 = 1.3083_RK                                                       !
          b3 = 0.15180_RK                                                      !
          b4 = 0.082349_RK                                                     !
          p  = 0.75_RK                                                         !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
        case( mALF_key+RPA_key )  ! mALF, RPA                                  !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
          a  = 0.016887_RK                                                     !
          a1 = 0.028829_RK                                                     !
          b1 = 10.357_RK                                                       !
          b2 = 3.6231_RK                                                       !
          b3 = 0.47990_RK                                                      !
          b4 = 0.12279_RK                                                      !
          p  = 1.0_RK                                                          !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      end select                                                               !
      !                                                                        !
      !------------------------------------------------------------------------+
    end subroutine get_pars_pw_lda_correlation                                 !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
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
    elemental function relative_spin_pol( ra, rb, rt ) result(z4)              !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: ra, rb, rt                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: z4                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      z4 = (( ra - rb ) / rt )**4                                              !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function relative_spin_pol                                             !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function spin_interpolation( ra, rb, rt ) result(f)              !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad),    intent(in) :: ra, rb, rt                                    !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: f                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      f = ( ( 2.0_RK * ra / rt )**n4d3 + ( 2.0_RK * rb / rt )**n4d3 - 2.0_RK ) &
                                                      / ( 2.0**n4d3 - 2.0_RK ) !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function spin_interpolation                                            !
    !                                                                          !
    !--------------------------------------------------------------------------+
    !                                                                          !
    elemental function G_form(a, a1, b1, b2, b3, b4, p, rs) result(G)          !
      !------------------------------------------------------------------------+
      !                                                                        !
      real(RK),    intent(in) :: a                                             !
      real(RK),    intent(in) :: a1                                            !
      real(RK),    intent(in) :: b1                                            !
      real(RK),    intent(in) :: b2                                            !
      real(RK),    intent(in) :: b3                                            !
      real(RK),    intent(in) :: b4                                            !
      real(RK),    intent(in) :: p                                             !
      !                                                                        !
      type(ad),    intent(in) :: rs                                            !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: G                                             !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      type(ad)                :: den                                           !
      !                                                                        !
      !------------------------------------------------------------------------+
      !                                                                        !
      den = 2.0_RK * a * ( b1 * rs**0.5_RK + b2 * rs                           &
                                     + b3 * rs**1.5_RK + b4 * rs**(p+1.0_RK) ) &
          + DTol ! Screen !!                                                   !
      !                                                                        !
      G   = -2.0_RK * a * ( 1.0_RK + a1 * rs ) * log( 1.0_RK + 1.0_RK / den )  !
      !                                                                        !
      !------------------------------------------------------------------------+
    end function G_form                                                        !
    !                                                                          !
    !--------------------------------------------------------------------------+
  end function pw_lda_correlation                                              !
  !                                                                            !
  !----------------------------------------------------------------------------+
end module pw_lda

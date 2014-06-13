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
!===============================================================
! Public interface of module
!===============================================================
module libdftauto
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: ...
  !  Date: ...
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  ! real with 8 bytes, precision 15 decimal digits
  integer, parameter :: rk = selected_real_kind(15)

  ! integer with 4 bytes, range 9 decimal digits
  integer, parameter :: ik = selected_int_kind(9)

  integer(IK), parameter, public  :: &
    X_LDA      = 2**0 , &
    C_VWN5     = 2**1 , &
    C_VWN5RPA  = 2**2 , &
    X_B3       = 2**3 , &
    X_B88      = 2**4 , &
    X_FT97B    = 2**5 , &
    X_PBE      = 2**6 , &
    X_PW91     = 2**7 , &
    C_FT97     = 2**8 , &
    C_LYP      = 2**9 , &
    C_P86      = 2**10, &
    C_PBE      = 2**11, &
    C_PW91     = 2**12, &
    C_PW92     = 2**13, &
    C_PZ81     = 2**14, &
    XC_B3LYP   = 2**15, &
    XC_B97_1   = 2**16, &
    XC_B97_2   = 2**17, &
    XC_B97     = 2**18, &
    XC_EDF1    = 2**19, &
    XC_FT97    = 2**20, &
    XC_HCTH120 = 2**21, &
    XC_HCTH147 = 2**22, &
    XC_HCTH407 = 2**23, &
    XC_HCTH    = 2**24, &
    XC_PW91    = 2**25

  integer(IK), parameter, public  :: &
    XC_LDA_MASK = X_LDA + C_VWN5 + C_VWN5RPA, &
    XC_GGA_MASK = 2**26-1 - XC_LDA_MASK


  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: xc_lib

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----


  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine xc_lib(imode,vl,ispin,ider,ro,gm,f,fr,fg,frr,frg,fgg)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(ik), intent(in)    :: imode
    integer(ik), intent(in)    :: vl, ispin
    integer(ik), intent(in)    :: ider
    real(rk)   , intent(in)    :: ro(:,:)  ! (vl,ispin)
    real(rk)   , intent(in)    :: gm(:,:)  ! (vl,1/3)
    ! gm(:,1)=g_aa, gm(:,2)=g_bb, gm(:,3)=g_ab
    real(rk)   , intent(inout) :: f(:)     ! (vl)
    real(rk)   , intent(inout) :: fr(:,:)  ! (vl,1/2)
    real(rk)   , intent(inout) :: fg(:,:)  ! (vl,1/3)
    real(rk)   , intent(inout) :: frr(:,:) ! (vl,1/3)
    real(rk)   , intent(inout) :: frg(:,:) ! (vl,1/6)
    real(rk)   , intent(inout) :: fgg(:,:) ! (vl,1/6)
    optional :: frr, frg, fgg
    ! derivatives of rho with respect to rho and gam
    ! *** end of interface ***

    if(imode==0) RETURN

    select case (ispin)

    case (1)

      select case(ider)
      case (1)
        call rks1(imode,vl,ro(:,1),gm(:,1),f,fr(:,1),fg(:,1))
      case (2)
        call rks2(imode,vl,ro(:,1),gm(:,1),f,fr(:,1),fg(:,1) &
                 ,frr(:,1),frg(:,1),fgg(:,1)                 &
                 )
      case default
        stop "no such ider in rks"
      end select

    case (2)

      select case(ider)
      case (1)
        call uks1(imode,vl,ro(:,:),gm(:,:),f,fr(:,:),fg(:,:))
      case default
        stop "no such ider in uks"
      end select

    case default
      stop "no such ispin"
    end select
  end subroutine xc_lib

  subroutine rks1(imode,vl,ro,gm,f,fr,fg)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(ik), intent(in)    :: imode
    integer(ik), intent(in)    :: vl
    real(rk)   , intent(in)    :: ro(:)   ! (vl)
    real(rk)   , intent(in)    :: gm(:)   ! (vl)
    real(rk)   , intent(inout) :: f(:)    ! (vl)
    real(rk)   , intent(inout) :: fr(:) ! (vl)
    real(rk)   , intent(inout) :: fg(:) ! (vl)
    ! *** end of interface ***

    integer(ik), parameter  :: ider=1
    real(rk), dimension(vl) :: x0,xr,xg,xrr,xrg,xgg

    if ( is_on( X_LDA    ) ) then
      call  RKS_X_LDA (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr
    endif
    if ( is_on( C_VWN5   ) ) then
      call  RKS_C_VWN5(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr
    endif
    if ( is_on( X_B88    ) ) then
      call  RKS_X_B88 (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
    if ( is_on( C_P86    ) ) then
      call  RKS_C_P86 (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
    if ( is_on( X_PBE    ) ) then
      call  RKS_X_PBE (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
    if ( is_on( C_PBE    ) ) then
      call  RKS_C_PBE (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
!   if ( is_on( C_VWN5RPA) ) then
!   endif
!   if ( is_on( X_B3     ) ) then
!   endif
!   if ( is_on( X_FT97B  ) ) then
!   endif
    if ( is_on( X_PW91   ) ) then
      call  RKS_X_PW91(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
!   if ( is_on( C_FT97   ) ) then
!   endif
!   if ( is_on( C_LYP    ) ) then
!   endif
    if ( is_on( C_PW91   ) ) then
      call  RKS_C_PW91(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
!   if ( is_on( C_PW92   ) ) then
!   endif
!   if ( is_on( C_PZ81   ) ) then
!   endif
!   if ( is_on( XC_B3LYP ) ) then
!   endif
!   if ( is_on( XC_B97_1 ) ) then
!   endif
!   if ( is_on( XC_B97_2 ) ) then
!   endif
!   if ( is_on( XC_B97   ) ) then
!   endif
!   if ( is_on( XC_EDF1  ) ) then
!   endif
!   if ( is_on( XC_FT97  ) ) then
!   endif
    if ( is_on( XC_HCTH  ) ) then
      call  RKS_XC_HCTH(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
    if ( is_on( XC_HCTH120)) then
      call  RKS_XC_HCTH120(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
    if ( is_on( XC_HCTH147)) then
      call  RKS_XC_HCTH147(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
    if ( is_on( XC_HCTH407)) then
      call  RKS_XC_HCTH407(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
    endif
!   if ( is_on( XC_PW91  ) ) then
!   endif
    contains
      function is_on(xc) result(yes)
        implicit none
        integer(IK), intent(in) :: xc
        logical                 :: yes
        ! *** end of interface ***

        yes = IAND(imode,xc) /= 0
      end function is_on
  end subroutine rks1

  subroutine rks2(imode,vl,ro,gm,f,fr,fg,frr,frg,fgg)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(ik), intent(in)    :: imode
    integer(ik), intent(in)    :: vl
    real(rk)   , intent(in)    :: ro(:)  ! (vl)
    real(rk)   , intent(in)    :: gm(:)  ! (vl)
    real(rk)   , intent(inout) :: f(:)   ! (vl)
    real(rk)   , intent(inout) :: fr(:)  ! (vl)
    real(rk)   , intent(inout) :: fg(:)  ! (vl)
    real(rk)   , intent(inout) :: frr(:) ! (vl)
    real(rk)   , intent(inout) :: frg(:) ! (vl)
    real(rk)   , intent(inout) :: fgg(:) ! (vl)
    ! *** end of interface ***

    integer(ik), parameter  :: ider=2
    real(rk), dimension(vl) :: x0,xr,xg,xrr,xrg,xgg

    if ( is_on( X_LDA    ) ) then
      call  RKS_X_LDA (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr
      frr=frr+xrr/2
    endif
    if ( is_on( C_VWN5   ) ) then
      call  RKS_C_VWN5(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr
      frr=frr+xrr/2
    endif
    if ( is_on( X_B88    ) ) then
      call  RKS_X_B88 (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
    if ( is_on( C_P86    ) ) then
      call  RKS_C_P86 (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
    if ( is_on( X_PBE    ) ) then
      call  RKS_X_PBE (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
    if ( is_on( C_PBE    ) ) then
      call  RKS_C_PBE (ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
!   if ( is_on( C_VWN5RPA) ) then
!   endif
!   if ( is_on( X_B3     ) ) then
!   endif
!   if ( is_on( X_FT97B  ) ) then
!   endif
    if ( is_on( X_PW91   ) ) then
      call  RKS_X_PW91(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
!   if ( is_on( C_FT97   ) ) then
!   endif
!   if ( is_on( C_LYP    ) ) then
!   endif
    if ( is_on( C_PW91   ) ) then
      call  RKS_C_PW91(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
!   if ( is_on( C_PW92   ) ) then
!   endif
!   if ( is_on( C_PZ81   ) ) then
!   endif
!   if ( is_on( XC_B3LYP ) ) then
!   endif
!   if ( is_on( XC_B97_1 ) ) then
!   endif
!   if ( is_on( XC_B97_2 ) ) then
!   endif
!   if ( is_on( XC_B97   ) ) then
!   endif
!   if ( is_on( XC_EDF1  ) ) then
!   endif
!   if ( is_on( XC_FT97  ) ) then
!   endif
    if ( is_on( XC_HCTH  ) ) then
      call  RKS_XC_HCTH(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
    if ( is_on( XC_HCTH120)) then
      call  RKS_XC_HCTH120(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
    if ( is_on( XC_HCTH147)) then
      call  RKS_XC_HCTH147(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
    if ( is_on( XC_HCTH407)) then
      call  RKS_XC_HCTH407(ider,vl,ro,gm,x0,xr,xg,xrr,xrg,xgg)
      f=f+x0; fr=fr+xr; fg=fg+xg/4
      frr=frr+xrr/2; frg=frg+xrg/4; fgg=fgg+xgg/16
    endif
!   if ( is_on( XC_PW91  ) ) then
!   endif
    contains
      function is_on(xc) result(yes)
        implicit none
        integer(IK), intent(in) :: xc
        logical                 :: yes
        ! *** end of interface ***

        yes = IAND(imode,xc) /= 0
      end function is_on
  end subroutine rks2

  subroutine uks1(imode,vl,ro,gm,f,fr,fg)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(ik), intent(in)    :: imode
    integer(ik), intent(in)    :: vl
    real(rk)   , intent(in)    :: ro(:,:)   ! (vl,2)
    real(rk)   , intent(in)    :: gm(:,:)   ! (vl,3)
    real(rk)   , intent(inout) :: f(:)      ! (vl)
    real(rk)   , intent(inout) :: fr(:,:) ! (vl,2)
    real(rk)   , intent(inout) :: fg(:,:) ! (vl,3)
    ! *** end of interface ***

    integer(ik), parameter    :: ider=1
    real(rk), dimension(vl  ) :: x0
    real(rk), dimension(vl,2) :: xr
    real(rk), dimension(vl,3) :: xg
    real(rk), dimension(vl  ) :: a

    if ( is_on( X_LDA    ) ) then
      call  UKS_X_LDA (ider,vl,ro(:,1),ro(:,2),gm(:,1),gm(:,2),gm(:,3), &
                            x0,xr(:,1),xr(:,2), a,a,a,                  &
                            a,a,a, a,a,a,a,a,a, a,a,a,a,a,a )
      f=f+x0; fr(:,1:2)=fr(:,1:2)+xr
    endif
    if ( is_on( C_VWN5   ) ) then
      call  UKS_C_VWN5(ider,vl,ro(:,1),ro(:,2),gm(:,1),gm(:,2),gm(:,3), &
                            x0,xr(:,1),xr(:,2), a,a,a,                  &
                            a,a,a, a,a,a,a,a,a, a,a,a,a,a,a )
      f=f+x0; fr(:,1:2)=fr(:,1:2)+xr
    endif
    if ( is_on( X_B88    ) ) then
      call  UKS_X_B88(ider,vl,ro(:,1),ro(:,2),gm(:,1),gm(:,2),gm(:,3), &
                           x0,xr(:,1),xr(:,2),xg(:,1),xg(:,2),xg(:,3), &
                           a,a,a, a,a,a,a,a,a, a,a,a,a,a,a )
      f=f+x0; fr(:,1:2)=fr(:,1:2)+xr; fg(:,1:3)=fg(:,1:3)+xg
    endif
    if ( is_on( C_P86    ) ) then
      call  UKS_C_P86(ider,vl,ro(:,1),ro(:,2),gm(:,1),gm(:,2),gm(:,3), &
                           x0,xr(:,1),xr(:,2),xg(:,1),xg(:,2),xg(:,3), &
                           a,a,a, a,a,a,a,a,a, a,a,a,a,a,a )
      f=f+x0; fr(:,1:2)=fr(:,1:2)+xr; fg(:,1:3)=fg(:,1:3)+xg
    endif
    if ( is_on( X_PBE    ) ) then
      call  UKS_X_PBE(ider,vl,ro(:,1),ro(:,2),gm(:,1),gm(:,2),gm(:,3), &
                           x0,xr(:,1),xr(:,2),xg(:,1),xg(:,2),xg(:,3), &
                           a,a,a, a,a,a,a,a,a, a,a,a,a,a,a )
      f=f+x0; fr(:,1:2)=fr(:,1:2)+xr; fg(:,1:3)=fg(:,1:3)+xg
    endif
    if ( is_on( C_PBE    ) ) then
      call  UKS_C_PBE(ider,vl,ro(:,1),ro(:,2),gm(:,1),gm(:,2),gm(:,3), &
                           x0,xr(:,1),xr(:,2),xg(:,1),xg(:,2),xg(:,3), &
                           a,a,a, a,a,a,a,a,a, a,a,a,a,a,a )
      f=f+x0; fr(:,1:2)=fr(:,1:2)+xr; fg(:,1:3)=fg(:,1:3)+xg
    endif
!   if ( is_on( C_VWN5RPA) ) then
!   endif
!   if ( is_on( X_B3     ) ) then
!   endif
!   if ( is_on( X_FT97B  ) ) then
!   endif
!   if ( is_on( X_PW91   ) ) then
!   endif
!   if ( is_on( C_FT97   ) ) then
!   endif
!   if ( is_on( C_LYP    ) ) then
!   endif
!   if ( is_on( C_PW91   ) ) then
!   endif
!   if ( is_on( C_PW92   ) ) then
!   endif
!   if ( is_on( C_PZ81   ) ) then
!   endif
!   if ( is_on( XC_B3LYP ) ) then
!   endif
!   if ( is_on( XC_B97_1 ) ) then
!   endif
!   if ( is_on( XC_B97_2 ) ) then
!   endif
!   if ( is_on( XC_B97   ) ) then
!   endif
!   if ( is_on( XC_EDF1  ) ) then
!   endif
!   if ( is_on( XC_FT97  ) ) then
!   endif
!   if ( is_on( XC_HCTH120 ) then
!   endif
!   if ( is_on( XC_HCTH147 ) then
!   endif
!   if ( is_on( XC_HCTH407 ) then
!   endif
!   if ( is_on( XC_HCTH  ) ) then
!   endif
!   if ( is_on( XC_PW91  ) ) then
!   endif
    contains
      function is_on(xc) result(yes)
        implicit none
        integer(IK), intent(in) :: xc
        logical                 :: yes
        ! *** end of interface ***

        yes = IAND(imode,xc) /= 0
      end function is_on
  end subroutine uks1
  !--------------- End of module ----------------------------------
end module libdftauto

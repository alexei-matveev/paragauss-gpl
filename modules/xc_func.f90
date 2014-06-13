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
module xc_func
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  ! Copyright (c) Alexey Shor
  ! Copyright (c) Sergey Bosko
  ! Copyright (c) Thomas Soini

  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: SB
  ! Date:   03/2007
  ! Description: Second derivatives have been added, from now on
  ! EXCHANGE:
  !       Xalpha
  !       BECKE88
  !       PBE
  !       PW91
  !       PBEN
  !       revPBE
  ! CORRELATION:
  !       VWN
  !       PW_lda
  !       PERDEW
  !       PBE
  !       PW91
  !
  ! TO TURN THE OLD CODE ON UNCOMMENT THIS:
# define OLD_CODE
  !
  !----------------------------------------------------------------
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

# include "def.h"
  use type_module ! type specification parameters
  use xc_cntrl, only: xc_NXC, xc_NLDAXC
  USE_MEMLOG
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: xc_functionals
  public :: xc_func_reset
  public :: xc_func_reduce
  public :: xc_func_write
  public :: xc_func_get_exc

#ifdef WITH_GUILE
  public :: qm_xc!(xc, ra, rb, gaa, gab, gbb) result (f) bind(c)
#endif

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  type, private :: xc_pot
    real(r8_kind)     :: contrib(1:xc_NXC)
    character(len=32) :: name
  end type xc_pot

  !------------ Declaration of constants and variables ----
  logical :: &
       gga_mode, &
       tau_mode, &
       hyb_mode, &
       eny_mode, &
       scf_mode

  integer(i4_kind)   :: xc_flavors(1:xc_NXC) = 0
  ! to track e.g. VWN/PWLDA use for PW91c
  integer(i4_kind), parameter :: &
       VWN_BASED   = 1, &
       PWLDA_BASED = 2

  real(kind=r8_kind) :: energy(1:xc_NXC) ! xc_NXC from xc_cntrl
  ! array for the values of the functionals

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  real(r8_kind) function frac(op)
    ! wrapper for this one:
    use xc_cntrl, only: XC_CNTRL_ON=>IS_ON, xc_cntrl_frac                      &
                                          , xc_M06L_X   , xc_M06L_C            &
                                          , xc_M06_X    , xc_M06_C             &
                                          , xc_VS_X     , xc_VS_C              &
                                          , xc_TPSS_X   , xc_TPSS_C            &
                                          , xc_EXX      , xc_HF
    implicit none
    integer(i4_kind), intent(in) :: op
    ! *** end of interface ***

    !------------ Declaration of local variables -----------------
    logical                      :: is_on
    !------------ Executable code --------------------------------

    frac = 0.0_r8_kind
    if( .not.scf_mode )then
       is_on = .true.
       !
       ! dont care about total fxc and derivatives, call all of them,
       ! except that in SO calculations kinetic energy density required
       ! for MGGA is not yet available. FIXME:
       if( .not. tau_mode .and. (op == xc_M06L_X    .or. op == xc_M06L_C       &
                            .or. op == xc_M06_X     .or. op == xc_M06_C        &
                            .or. op == xc_VS_X      .or. op == xc_VS_C         &
                            .or. op == xc_TPSS_X    .or. op == xc_TPSS_C  )    &
         ) then
         is_on = .false.
       endif
         !
         ! The same holds for hybrid functionals
       if( .not. hyb_mode .and. (op == xc_EXX)     .or. op == xc_HF            &
         ) then
         is_on = .false.
       endif
       ! set full contribution FIXME: find a way to deal with the non-pure ones
       if ( is_on ) frac = 1.0_r8_kind
    else
       ! call only the requested ones:
       is_on = XC_CNTRL_ON(op)
       !
       ! set frac value according to xc_cntrl
       frac = xc_cntrl_frac( op )
       !
       ! check if options match
       ASSERT(frac/=0.0.eqv.is_on)
       !
    endif
    !
  end function frac

  subroutine xc_functionals (vla, ispin, rho, f, fr, gam, fg, grdwts, &
       frr, frg, fgg, tau, ft)
    !
    ! Purpose: Calculates the XC functional and derivatives
    !
    use xc_cntrl, only: XC_CNTRL_ON=>IS_ON, xc_EXX, xc_HF
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(in) :: vla, ispin
    real(r8_kind), intent(in) :: rho(:, :), gam(:, :), tau(:, :)
    real(r8_kind), intent(out) :: f(:), fr(:, :), fg(:, :), ft(:, :)
    real(r8_kind), intent(in) :: grdwts(:)
    real(r8_kind), intent(out) :: frr(:, :)
    real(r8_kind), intent(out) :: frg(:, :)
    real(r8_kind), intent(out) :: fgg(:, :)
    optional :: f, fr         ! for scf only (and PH/GRADS)
    optional :: gam, fg       ! for GGA only
    optional :: tau, ft       ! for MGGA only
    optional :: grdwts        ! for energy calc only
    optional :: frr, frg, fgg ! for sec ders only
    !** End of interface *****************************************

    !------------ Declaration of local variables -----------------
    real(r8_kind) :: x(vla), xr(vla,ispin), xg(vla,2*ispin-1),xt(vla,ispin)
    real(r8_kind) :: y(vla), yr(vla,ispin), yrr(vla,2*ispin-1) ! for e.g. perdewwang
    real(r8_kind) :: xc_frac
    real(r8_kind), parameter :: zero = 0.0_r8_kind
    integer(i4_kind) :: lda_flavor
    integer(i4_kind) :: xc, NXC
    integer(i4_kind) :: ider
    !------------ Executable code --------------------------------

    ASSERT(size(rho,1)>=vla)
    ASSERT(size(rho,2)>=ispin)

    ider = 1
    if(present(frr)) ider = 2

    gga_mode = present(gam)
    tau_mode = present(tau)
    hyb_mode = XC_CNTRL_ON(xc_EXX) .or. XC_CNTRL_ON(xc_HF)

    if(gga_mode)then
       ASSERT(size(gam,1)>=vla)
       ASSERT(size(gam,2)>=2*ispin-1)

       if(tau_mode)then
          ASSERT(size(tau,1)>=vla)
          ASSERT(size(tau,2)>=ispin)
       endif

       if( ider == 2 )then
         ASSERT(present(frg))
         ASSERT(present(fgg))
       endif
    endif

    eny_mode = present(grdwts)
    scf_mode = present(f)
    ASSERT(eny_mode.or.scf_mode)

    if(scf_mode)then
       ASSERT(present(fr))
       ASSERT(size(fr,1)>=vla)
       ASSERT(size(fr,2)>=ispin)
       if(gga_mode)then
          ASSERT(present(fg))
          ASSERT(size(fg,1)>=vla)
          ASSERT(size(fg,2)>=2*ispin-1)
          if(tau_mode)then
             ASSERT(present(ft))
             ASSERT(size(ft,1)>=vla)
             ASSERT(size(ft,2)>=ispin)
          endif
       endif
    endif

    ! "f" and its derivatives are initialized here, so this is intent(out)!
    if(present(f))   f(:vla)     = zero
    if(present(fr))  fr(:vla,:)  = zero
    if(present(fg))  fg(:vla,:)  = zero
    if(present(ft))  ft(:vla,:)  = zero
    if(present(frr)) frr(:vla,:) = zero
    if(present(frg)) frg(:vla,:) = zero
    if(present(fgg)) fgg(:vla,:) = zero

    if( gga_mode )then
       NXC = xc_NXC    ! from xc_cntrl
    else
       NXC = xc_NLDAXC ! from xc_cntrl
    endif

    DPRINT 'xcf: scf_mode=',scf_mode,' eny_mode=',eny_mode

    lda_flavor = 0

    do xc = 1, NXC ! not more than in the case list in xc_no()
       xc_frac = frac(xc)
       if( xc_frac /= 0.0_r8_kind ) then
          x  = zero
          xr = zero
          xg = zero
          xt = zero
          call xc_no(xc)
          if(eny_mode) call add_eny(energy(xc))
          if(scf_mode)then
             f(1:vla)          = f(1:vla)                   + x  * xc_frac
             fr(1:vla,1:ispin) = fr(1:vla,1:ispin)          + xr * xc_frac
             if(gga_mode)then
                fg(1:vla,:2*ispin-1) = fg(1:vla,:2*ispin-1) + xg * xc_frac
                if(tau_mode)then
                   ft(1:vla,1:ispin) = ft(1:vla,1:ispin)    + xt * xc_frac
                end if
             endif
          endif
       endif
    enddo

    DPRINT 'xcf: exit'
  contains

    subroutine add_eny(eny)
      ! adds energy contribution
      implicit none
      real(r8_kind), intent(inout) :: eny
      ! *** end of interface **

      eny = eny + sum(x(1:vla)*grdwts(1:vla))
    end subroutine add_eny

    subroutine xc_no(no)
      use xc_cntrl, xc_cntrl_is_on=>is_on
#ifdef WITH_LIBDFTAUTO
      use libdftauto
#else
      use vwnc, only: vwn_ldac
      use becke_perdew_module, only: perdew_calc
      use perdew_wang_module, only: pw91c_calc
      use pbe_ggcxc_module, only: pbe_ggcc
      use pw_ldac_module, only: pw_ldac
      use baerends_module, only: baerends94_calc
      use relxc, only: rel_vwn_calc, rpw91c_calc
      use hcth, only: hcth_x, hcth_c
      use m06 , only: m06_X, m06_C
      use lyp, only: lyp_C
      use vsxc, only: vsxc_X, vsxc_C
      use tpss, only: tpss_X, tpss_C, TPSS_key, revTPSS_key
      use tpss_mgga_module, only: tpss_mgga_X, tpss_mgga_C
      use exchange
#ifdef WITH_RESPONSE
      use gga_response_module
#endif
#endif
      implicit none
      integer(i4_kind), intent(in) :: no
      ! *** end of interface ***

      integer(i4_kind) :: IX

#ifndef WITH_LIBDFTAUTO

      select case (no)

      case (xc_vwn)
         DPRINT 'xcf: ',no,' call vwn_calc()'

         select case(ider)
         case (1)
           call vwn_ldac (vla, ispin, rho, x, xr)
         case (2)
           call vwn_ldac (vla, ispin, rho, x, xr, frr)
         end select

         ! in scf only one (if any) of two LDAc can be used:
         if(xc_cntrl_is_on(xc_vwn))then
            ! these maybe used later for GGAs
            y  = x
            yr = xr
            if( ider == 2 )then
              ! also save the second derivative wrt density,
              ! it has one component in unpolarized case (ispin == 1)
              ! and tree components aa, bb, ab in polarized case (ispin == 2)
              ! However I am not sure if there are even four redundant
              ! components, both ab and ba there:
              ASSERT(size(frr,2)==2*ispin-1)
              yrr(:,:2*ispin-1) = frr(:vla,:2*ispin-1)
            endif
            lda_flavor = VWN_BASED
         endif

      case (xc_pwldac)
         DPRINT 'xcf: ',no,' call pw_ldac()'

         select case(ider)
         case (1)
            call pw_ldac (vla, ispin, rho, x, xr)
         case (2)
            ! FIXME:  before it  was called  with ispin  + 2,  that is
            ! either with PWLDAC_RESP_R or with PWLDAC_RESP_R:
            call pw_ldac (vla, ispin, rho, x, xr, frr)
         end select

         ! in scf only one (if any) of two LDAc can be used:
         if(lda_flavor==0)then
            ! these maybe used later for GGAs
            y  = x
            yr = xr
            if( ider == 2 )then
              ! also save the second derivative wrt density,
              ! it has one component in unpolarized case (ispin == 1)
              ! and tree components aa, bb, ab in polarized case (ispin == 2)
              ! However I am not sure if there are even four redundant
              ! components, both ab and ba there:
              ASSERT(size(frr,2)==2*ispin-1)
              yrr(:,:2*ispin-1) = frr(:vla,:2*ispin-1)
            endif
            lda_flavor = PWLDA_BASED
         endif

      case (xc_perdewwang91c)
         ASSERT(lda_flavor/=0)
         DPRINT 'xcf: ',no,' call pw91c_calc()'

#ifdef OLD_CODE
         ASSERT(ider<2)
         call pw91c_calc (rho,gam,xr,ispin,x,xg,vla,&
                          y,yr)
#else
         !! 4 for PW91
         select case(ider)
         case (1)
            call gga_correlation(4,ispin,ider,rho,gam,vla,x,y,xr,xg,yr)
         case (2)
            call gga_correlation(4,ispin,ider,rho,gam,vla,x,y,xr,xg,yr,&
                 dfdndn=frr,dfdndg=frg,dfdgdg=fgg,dEcdndn=yrr)
         end select

#endif
         xc_flavors(xc_perdewwang91c) = lda_flavor

      case (xc_revPW91c)
         ASSERT(lda_flavor/=0)
         ASSERT(ider<2)
         DPRINT 'xcf: ',no,' call pw91c_calc(rev)'
         call pw91c_calc (rho,gam,xr,ispin,x,xg,vla,&
              & y,yr, revised=.true.)
         xc_flavors(xc_revPW91c) = lda_flavor

      case (xc_rperdewwang91c)
         ASSERT(lda_flavor/=0)
         ASSERT(ider<2)
         DPRINT 'xcf: ',no,' call rpw91c_calc()'
         ! NO: x = zero
         ! FIXME: do something about it:
         call pw91c_calc (rho,gam,xr,ispin,x,xg,vla,&
              & y,yr )
         ! this one scales the previous result:
         call rpw91c_calc(rho,gam,ispin,vla,x,xr,xg)
         xc_flavors(xc_rperdewwang91c) = lda_flavor

      case (xc_rvwn)
         select case(whatis(xc_rvwn))
         case(0,1)
            DPRINT 'xcf: ',no,' call rel_vwn_calc()'
            call rel_vwn_calc(rho,ispin,vla,x,xr)
         case(2)
            DPRINT 'xcf: ',no,' call correlation_lda()'
            call correlation_lda(vla,ispin,rho,x,xr)
         case default
            ABORT('no such xc_rvwn')
         end select

      case (xc_perdew)
         DPRINT 'xcf: ',no,' call perdew_calc()'

         select case(ider)
         case (1)
           call perdew_calc(rho,gam,xr,ispin,x,xg, vla)
         case (2)
           call perdew_calc(rho,gam,xr,ispin,x,xg, vla, &
                          df_drhodrho=frr, &
                          df_drhodgamma=frg, &
                          df_dgammadgamma=fgg)
         end select

      case (xc_pbec)  !pbe(1)
         ASSERT(lda_flavor/=0)
         DPRINT 'xcf: ',no,' call pbe_ggcc(PBE)'

#if 1
         select case(ider)
         case (1)
            call pbe_ggcc("pbe", vla, ispin, rho, gam, x, xr, xg)
         case (2)
            call pbe_ggcc("pbe", vla, ispin, rho, gam, x, xr, xg, frr, frg, fgg)
         end select
#else
         !! 3 for PBE
         select case(ider)
         case (1)
            call gga_correlation(3,ispin,ider,rho,gam,vla,x,y,xr,xg,yr)
         case (2)
            call gga_correlation(3,ispin,ider,rho,gam,vla,x,y,xr,xg,yr,&
                 dfdndn=frr,dfdndg=frg,dfdgdg=fgg,dEcdndn=yrr)
         end select
#endif
         xc_flavors(xc_pbec) = lda_flavor

      case (xc_pbesolc)
         ASSERT(lda_flavor/=0)
         DPRINT 'xcf: ', no, ' call pbe_ggcc(PBEsol)'

         select case(ider)
         case (1)
            call pbe_ggcc("pbesol", vla, ispin, rho, gam, x, xr, xg)
         case (2)
            call pbe_ggcc("pbesol", vla, ispin, rho, gam, x, xr, xg, frr, frg, fgg)
         end select

         xc_flavors(xc_pbesolc) = lda_flavor

      case (xc_baerends94)
         ASSERT(ider<2)
         ! has no energy functional:
         DPRINT 'xcf: ',no,' call baerends94_calc()'
         call baerends94_calc (rho,gam,xr,ispin,vla)
         !if(eny_mode) call add_eny(exc_baerends94)

      case (xc_hcth_x)
         ASSERT(ider<2)
         DPRINT 'xcf: ',no,' call hcth_x()'
         call hcth_x(rho,gam,ispin,vla,whatis(xc_hcth_version),x,xr,xg)

      case (xc_hcth_c)
         ASSERT(ider<2)
         DPRINT 'xcf: ',no,' call hcth_c()'
         call hcth_c(rho,gam,ispin,vla,whatis(xc_hcth_version),x,xr,xg)

      case (xc_xalpha)
         IX = X_XALPHA
         DPRINT 'xcf: ',no,' call exchange_lda(Xa',IX,')'

         select case(ider)
         case (1)
            call exchange_lda(IX,vla,ispin,rho,x,xr)
         case (2)
            call exchange_lda(IX,vla,ispin,rho,x,xr,frr)
         end select

      case (xc_rxalpha)
         ASSERT(ider<2)
         if(xc_cntrl_is_on(xc_rel))then
            IX = X_XALPHA + whatis(xc_rel)
         else
            IX = X_XALPHA + X_LONGTRANS
         endif
         DPRINT 'xcf: ',no,' call exchange_lda(RXa',IX,')'
         call exchange_lda(IX,vla,ispin,rho,x,xr)

      case (xc_becke88)
         IX = X_BECKE88
         DPRINT 'xcf: ',no,' call exchange_gga(BECKE88',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_rbecke88)
         IX = X_BECKE88 + X_LONGTRANS
         DPRINT 'xcf: ',no,' call exchange_gga(BECKE88+REL',IX,')'
         ASSERT(ider<2)
         call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)

      case (xc_perdewwang91x)
         IX = X_PW91
         DPRINT 'xcf: ',no,' call exchange_gga(PW91',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_rperdewwang91x)
         ASSERT(ider<2)
         IX = X_PW91 + X_LONGTRANS
         DPRINT 'xcf: ',no,' call exchange_gga(PW91+REL',IX,')'
         call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)

      case (xc_pbex) !pbe(2)
         IX = X_PBE + whatis(xc_rel)
         DPRINT 'xcf: ',no,' call exchange_gga(PBE',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_revPBEx)
         IX = X_REVPBE + whatis(xc_rel)
         DPRINT 'xcf: ',no,' call exchange_gga(revPBE',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_pbesolx)
         IX = X_PBESOL ! FIXME: + whatis(xc_rel), check that
         DPRINT 'xcf: ', no, ' call exchange_gga(PBESOLX', IX, ')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_PBENx)
         IX = X_PBEN + whatis(xc_rel)
         DPRINT 'xcf: ',no,' call exchange_gga(PBEN',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_ecmv92)
         IX = X_ECMV92 + X_LONGTRANS
         DPRINT 'xcf: ',no,' call exchange_gga(ECMV92+REL',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_nrecmv92)
         IX = X_ECMV92
         DPRINT 'xcf: ',no,' call exchange_gga(ECMV92',IX,')'

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_m06l_x)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_m06l_x),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call m06_X(1,vla,ispin,rho(:vla,:),gam(:vla,:),tau(:vla,:)&
                       ,x(:vla),xr(:vla,:),xg(:vla,:),xt(:vla,:))

      case (xc_m06l_c)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_m06l_c),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call m06_C(1,vla,ispin,rho(:vla,:),gam(:vla,:),tau(:vla,:)&
                       ,x(:vla),xr(:vla,:),xg(:vla,:),xt(:vla,:))

      case (xc_m06_x)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_m06_x),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call m06_X(2,vla,ispin,rho(:vla,:),gam(:vla,:),tau(:vla,:)&
                       ,x(:vla),xr(:vla,:),xg(:vla,:),xt(:vla,:))

      case (xc_m06_c)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_m06_c),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call m06_C(2,vla,ispin,rho(:vla,:),gam(:vla,:),tau(:vla,:)&
                       ,x(:vla),xr(:vla,:),xg(:vla,:),xt(:vla,:))

      case (xc_tpss_x)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_tpss_x),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call tpss_X(TPSS_key, vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)&
                       , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))

      case (xc_tpss_c)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_tpss_c),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
#if 0
         call tpss_C(TPSS_key, vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)&
                       , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))
#endif
         call tpss_mgga_C( vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)   &
                         , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))

!     case (xc_revtpss_x)
!        DPRINT   'xcf: call xc_meta_gga(',whatis(xc_revtpss_x),')'
!        ASSERT(gga_mode)
!        ASSERT(tau_mode)
!        ASSERT(ider.eq.0.or.ider.eq.1)
!        ! LEGACY VARIANT
!        WRITE(*,*) 'REVTPSS currently not supported'
#if 0
!        call tpss_mgga_X( vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)   &
!                        , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))
!        call tpss_X(revTPSS_key, vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)&
!                      , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))
#endif

!     case (xc_revtpss_c)
!        DPRINT   'xcf: call xc_meta_gga(',whatis(xc_revtpss_c),')'
!        ASSERT(gga_mode)
!        ASSERT(tau_mode)
!        ASSERT(ider.eq.0.or.ider.eq.1)
!        ! LEGACY VARIANT
!        WRITE(*,*) 'REVTPSS currently not supported'
#if 0
!        call tpss_mgga_C( vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)   &
!                        , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))
!        call tpss_C(revTPSS_key, vla, ispin, rho(:vla,:), gam(:vla,:), tau(:vla,:)&
!                      , x(:vla), xr(:vla,:), xg(:vla,:), xt(:vla,:))

#endif
      case (xc_vmt)
         IX = X_VMT
         DPRINT 'xcf: ',no,' call exchange_gga(VMT',IX,')'
         ASSERT(gga_mode)

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_vt84)
         IX = X_VT84
         DPRINT 'xcf: ',no,' call exchange_gga(VT84',IX,')'
         ASSERT(gga_mode)

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_vmtsol)
         IX = X_VMTsol
         DPRINT 'xcf: ',no,' call exchange_gga(VMTsol',IX,')'
         ASSERT(gga_mode)

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_vt84sol)
         IX = X_VT84sol
         DPRINT 'xcf: ',no,' call exchange_gga(VT84sol',IX,')'
         ASSERT(gga_mode)

         select case(ider)
         case (1)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg)
         case (2)
            call exchange_gga(IX,vla,ispin,rho,gam,x,xr,xg,frr,frg,fgg)
         end select

      case (xc_vs_x)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_vs_x),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call vsxc_X(vla,ispin,rho(:vla,:),gam(:vla,:),tau(:vla,:)&
                       ,x(:vla),xr(:vla,:),xg(:vla,:),xt(:vla,:))

      case (xc_vs_c)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_vs_c),')'
         ASSERT(gga_mode)
         ASSERT(tau_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call vsxc_C(vla,ispin,rho(:vla,:),gam(:vla,:),tau(:vla,:)&
                       ,x(:vla),xr(:vla,:),xg(:vla,:),xt(:vla,:))

      case (xc_lyp_c)
         DPRINT   'xcf: call xc_meta_gga(',whatis(xc_vs_c),')'
         ASSERT(gga_mode)
         ASSERT(ider.eq.0.or.ider.eq.1)
         call lyp_C( 1, vla, ispin, rho(:vla,:), gam(:vla,:)                   &
                                             , x(:vla), xr(:vla,:), xg(:vla,:) )

      case ( xc_EXX )
         DPRINT   'xcf: call nothing for EXX'
         ASSERT(hyb_mode)

      case ( xc_HF )
         DPRINT   'xcf: call nothing for HF'
         ASSERT(hyb_mode)

      case (xc_new)
        ! PLACE THE CALL TO YOUR TEST-FUNCTIONAL BELOW

      case default
         ABORT('no such case')
      end select

#else /* of ifndef WITH_LIBDFTAUTO */

      real(r8_kind)    :: NOg  (vla,2*ispin-1) ! 1 or 3
      real(r8_kind)    :: NOfg (vla,2*ispin-1) ! 1 or 3
      real(r8_kind)    :: NOfrg(vla,5*ispin-4) ! 1 or 6
      real(r8_kind)    :: NOfgg(vla,5*ispin-4) ! 1 or 6

      ! use external XC-implementations:


      IX = dftauto_no(no)

      if( IX == 0 ) RETURN

      if( .not. gga_mode )then
         select case(ider)
         case (1)
           call xc_lib(IX,vla,ispin,1,rho,NOg,x,fr,NOfg)
         case (2)
           call xc_lib(IX,vla,ispin,2,rho,NOg,x,fr,NOfg,frr,NOfrg,NOfgg)
         end select
      else
         select case(ider)
         case (1)
           call xc_lib(IX,vla,ispin,1,rho,gam,x,fr,fg)
         case (2)
           call xc_lib(IX,vla,ispin,2,rho,gam,x,fr,fg,frr,frg,fgg)
         end select
      endif
#endif /* of ifndef WITH_LIBDFTAUTO */

    end subroutine xc_no

#ifdef WITH_LIBDFTAUTO
    function dftauto_no(no) result(IX)
      use xc_cntrl
      use libdftauto
      implicit none
      integer(i4_kind), intent(in) :: no
      integer(i4_kind)             :: IX
      ! *** end of interface ***

      select case (no)

      ! VWN:
      case (xc_xalpha)
         IX = X_LDA
      case (xc_vwn)
         IX = C_VWN5
      ! BP:
      case (xc_becke88)
         IX = X_B88
      case (xc_perdew)
         IX = C_P86
      ! PBE:
      case (xc_pbex)
         IX = X_PBE
      case (xc_pbec)
         IX = C_PW91
      ! PW91:
      case (xc_perdewwang91x)
         IX = X_PW91
      case (xc_perdewwang91c)
         IX = C_PW91
      ! HCTH:
      case (xc_hcth_x)
         ! do nothing, wait until xc_hcth_c
         IX = 0
      case (xc_hcth_c) ! or xc_hcth_x
         select case(whatis(xc_hcth_version))
         case (1)
           IX = XC_HCTH
         case (2)
           IX = XC_HCTH120
         case (3)
           IX = XC_HCTH147
         case (4)
           IX = XC_HCTH407
         end select
      case default
         print *,'XC_NO(',no,'): not available'
         ABORT('recompile w/o -DWITH_LIBDFTAUTO')
      end select
    end function dftauto_no
#endif /* of ifndef WITH_LIBDFTAUTO */

  end subroutine xc_functionals

  subroutine xc_func_reset()
    implicit none
    ! *** end of interface ***

    energy = 0.0_r8_kind

    ! will be re-set in xc_functionals()
    scf_mode = .true.
    eny_mode = .false.
    gga_mode = .false.
    tau_mode = .false.
  end subroutine xc_func_reset

  subroutine xc_func_write(unit)
    use xc_cntrl, only: xc_cntrl_names
    implicit none
    integer(i4_kind), intent(in) :: unit
    ! *** end of interface ***

    integer(i4_kind) :: xc

    DPRINT 'xc_func_write: scf=',scf_mode,'eny=',eny_mode,'gga=',gga_mode

    write(unit,'(A6,A24)') 'PHXC: ','----------------'
    write(unit,*)
    write(unit,'(A6,A24)') 'PHXC: ','xc-energies:'
    do xc=1,xc_NXC
       write(unit,'(A6,A24,A1,F25.12)') 'PHXC: ',&
            trim(xc_cntrl_names(1,xc))//' '//trim(xc_cntrl_names(2,xc)), ':',&
            energy(xc)
    enddo
    write(unit,'(A6,A24)') 'PHXC: ','----------------'
    if( scf_mode )then
       write(unit,'(A6,A24,A1,F25.12)') 'PHXC: ','total',':',sum(energy)
    endif
  end subroutine xc_func_write

  function pot( name, xc1, xc2, xc3, xc4                                       &
                    , xc1_frac, xc2_frac, xc3_frac, xc4_frac ) result(p)
    ! xc_pot 'constructor'
    implicit none
    INTEGER(i4_kind), optional, intent(in) :: xc1
    INTEGER(i4_kind), optional, intent(in) :: xc2
    INTEGER(i4_kind), optional, intent(in) :: xc3
    INTEGER(i4_kind), optional, intent(in) :: xc4
    REAL(r8_kind)   , optional, intent(in) :: xc1_frac
    REAL(r8_kind)   , optional, intent(in) :: xc2_frac
    REAL(r8_kind)   , optional, intent(in) :: xc3_frac
    REAL(r8_kind)   , optional, intent(in) :: xc4_frac
    character(len=*), intent(in)  :: name

    type(xc_pot)                  :: p
    ! *** end of interface ***

    p%contrib      = 0.0_r8_kind
    p%name = repeat(' ',32)
    p%name = name

    if(present(xc1))then
       ASSERT(xc1>0)
       ASSERT(xc1<=xc_NXC)
       IF (present(xc1_frac)) THEN
         p%contrib(xc1) = xc1_frac
       ELSE
         p%contrib(xc1) = 1.0_r8_kind
       ENDIF
    endif
    if(present(xc2))then
       ASSERT(xc2>0)
       ASSERT(xc2<=xc_NXC)
       IF (present(xc2_frac)) THEN
         p%contrib(xc2) = xc2_frac
       ELSE
         p%contrib(xc2) = 1.0_r8_kind
       ENDIF
    endif
    if(present(xc3))then
       ASSERT(xc3>0)
       ASSERT(xc3<=xc_NXC)
       IF (present(xc3_frac)) THEN
         p%contrib(xc3) = xc3_frac
       ELSE
         p%contrib(xc3) = 1.0_r8_kind
       ENDIF
    endif
    if(present(xc4))then
       ASSERT(xc4>0)
       ASSERT(xc4<=xc_NXC)
       IF (present(xc4_frac)) THEN
         p%contrib(xc4) = xc4_frac
       ELSE
         p%contrib(xc4) = 1.0_r8_kind
       ENDIF
    endif

    DPRINT 'xcf::pot: constructed >'// p%name//'<',p%contrib
  end function pot

  function add_pot(p1,p2) result(psum)
    ! xc_pot 'constructor'
    use constants, only: zero, one
    implicit none
    type(xc_pot), intent(in) :: p1,p2
    type(xc_pot)             :: psum
    ! *** end of interface ***

    psum%contrib = p1%contrib    +    p2%contrib
    psum%name    = repeat(' ',32)
    psum%name    = trim(p1%name) //'+'// trim(p2%name)
    DPRINT 'xcf::add: constructed >'// psum%name//'<',psum%contrib
  end function add_pot

  subroutine xc_func_get_exc (a_xc_arr, name_xc_arr)
    ! moved from PH, was named post_hoc_get_exc
    ! use ph_cntrl
    use xc_cntrl, only: xc_cntrl_is_on=>is_on, xc_cntrl_frac                   &
                      , xc_xalpha        , xc_rxalpha                          &
                      , xc_pwldac        , xc_vwn           , xc_rvwn          &
                      , xc_becke88       , xc_rbecke88      , xc_perdew        &
                      , xc_perdewwang91x , xc_perdewwang91c , xc_revPW91c      &
                      , xc_rperdewwang91x, xc_rperdewwang91c                   &
                      , xc_pbex          , xc_pbec          , xc_pbenx         &
                      , xc_revpbex       , xc_pbesolx       , xc_pbesolc       &
                      , xc_vmt           , xc_vt84          , xc_lyp_c         &
                      , xc_ecmv92        , xc_nrecmv92                         &
                      , xc_HCTH_x        , xc_HCTH_c                           &
                      , xc_tpss_x        , xc_tpss_c                           &
                      , xc_vs_x          , xc_vs_c                             &
                      , xc_m06l_x        , xc_m06l_c
    !
    implicit none
    ! purpose : make the module private integrated xc_energyies
    !           for various functionals
    !           available to the calling unit
    real(r8_kind), allocatable :: a_xc_arr(:)
    character(len=32), allocatable  :: name_xc_arr(:)
    !** End of interface *****************************************

    integer(kind=i4_kind) :: counter,alloc_stat

    integer(i4_kind), parameter :: N_POT=36
    integer(i4_kind)            :: xc
    type(xc_pot)                :: pots(N_POT)
    type(xc_pot)                :: &
         VWNc, VWN, PWLDA, RVWN, &
         BP, BPW91, &
         PW91, PW91IIA, &
         PBE, revPBE, PBEsol, PBEN, &
         SCF

    pots(1)  = pot('XALPHA'     , xc_xalpha)
    pots(2)  = pot('RXALPHA'    , xc_rxalpha)

    VWNc     = pot('VWN'        , xc_vwn)
    VWN      = pot('VWN'        , xc_xalpha, xc_vwn)
    PWLDA    = pot('PWLDA'      , xc_xalpha, xc_pwldac)

    pots(3)  = VWN
    pots(4)  = PWLDA

    RVWN     = pot('RVWN'       , xc_rxalpha, xc_vwn, xc_rvwn)

    pots(5)  = RVWN

    BP       = pot('BP'         , xc_becke88,       xc_perdew)

    select case(xc_flavors(xc_perdewwang91c))
    case (PWLDA_BASED,0) ! zero == was not evaluated
       BPW91    = pot('BPW91'      , xc_becke88,       xc_perdewwang91c)
       PW91     = pot('PW91'       , xc_perdewwang91x, xc_perdewwang91c)
       PW91IIA  = pot('PW91IIA'    , xc_perdewwang91x, xc_revPW91c)
    case (VWN_BASED)
       BPW91    = pot('BPW91/VWN'  , xc_becke88,       xc_perdewwang91c)
       PW91     = pot('PW91/VWN'   , xc_perdewwang91x, xc_perdewwang91c)
       PW91IIA  = pot('PW91IIA/VWN', xc_perdewwang91x, xc_revPW91c)
    case default
       ABORT('no such case')
    end select

    pots(6)  = add_pot(VWN   , BP)
    pots(7)  = add_pot(VWN   , BPW91)
    pots(8)  = add_pot(PWLDA , BPW91)
    pots(9)  = add_pot(PWLDA , PW91)
    pots(10) = add_pot(PWLDA , PW91IIA)
    pots(11) = add_pot(VWN   , PW91)

    select case(xc_flavors(xc_pbec))
    case (PWLDA_BASED,0) ! zero == was not evaluated
       PBE     = pot('PBE'        , xc_pbex   , xc_pbec)
       revPBE  = pot('revPBE'     , xc_revPBEx, xc_pbec)
       PBEN    = pot('PBEN'       , xc_PBENx  , xc_pbec)
       PBEsol  = pot('PBEsol'     , xc_pbesolx, xc_pbesolc)
    case (VWN_BASED)
       PBE     = pot('PBE/VWN'    , xc_pbex   , xc_pbec)
       revPBE  = pot('revPBE/VWN' , xc_revPBEx, xc_pbec)
       PBEN    = pot('PBEN/VWN'   , xc_PBENx  , xc_pbec)
       PBEsol  = pot('PBEsol/VWN' , xc_pbesolx, xc_pbesolc)
    case default
       ABORT('no such case')
    end select

    pots(12) = add_pot(PWLDA , PBE)
    pots(13) = add_pot(PWLDA , revPBE)
    pots(14) = add_pot(PWLDA , PBEN)
    pots(15) = add_pot(PWLDA , PBEsol)

    pots(16) = add_pot(VWN   , PBE)
    pots(17) = add_pot(VWN   , revPBE)
    pots(18) = add_pot(VWN   , PBEN)
    pots(19) = add_pot(VWN   , PBEsol)

    pots(20) = pot('HCTH'       , xc_HCTH_x  , xc_HCTH_c)
    pots(21) = pot('ECMV92'     , xc_xalpha, xc_ecmv92)
    pots(22) = pot('nrECMV92'   , xc_xalpha, xc_nrecmv92)

    pots(23) = add_pot(RVWN , BPW91)

    select case(xc_flavors(xc_perdewwang91c))
    case (PWLDA_BASED,0) ! zero == was not evaluated
       pots(24) = add_pot(VWN , pot('BPW91[rx,c]'  , xc_rbecke88, xc_perdewwang91c))
       pots(25) = add_pot(VWN , pot('BPW91[rx,rc]' , xc_rbecke88, xc_rperdewwang91c))
       pots(26) = pot('PWLDA+PW91[rx,c]'   , xc_rperdewwang91x, xc_pwldac, xc_perdewwang91c)
       pots(27) = pot( 'RVWN+PW91[rx,c]'   , xc_rperdewwang91x, xc_vwn, xc_rvwn, xc_perdewwang91c)
    case (VWN_BASED)
       pots(24) =  add_pot(VWN , pot('BPW91[rx,c]/VWN' , xc_rbecke88      , xc_perdewwang91c))
       pots(25) =  add_pot(VWN , pot('BPW91[rx,rc]/VWN', xc_rbecke88      , xc_rperdewwang91c))
       pots(26) =  add_pot(VWN , pot('PW91[rx,c]/VWN'  , xc_rperdewwang91x, xc_perdewwang91c))
       pots(27) =  add_pot(pot('RVWN',xc_xalpha, xc_vwn, xc_rvwn) &
            &          , pot('PW91[rx,c]/VWN' , xc_rperdewwang91x, xc_perdewwang91c))
    case default
       ABORT('no such case')
    end select

    pots(28) = pot('M06L' , xc_m06l_x, xc_m06l_c)
    pots(29) = pot('TPSS' , xc_tpss_x, xc_tpss_c)

    pots(30) = pot('VMT'  , xc_vmt  , xc_pbec)
    pots(31) = pot('VT84' , xc_vt84 , xc_pbec)
    pots(32) = pot('VMTsol'  , xc_vmt  , xc_pbesolc)
    pots(33) = pot('VT84sol' , xc_vt84 , xc_pbesolc)

    pots(34) = pot('VSXC' , xc_vs_x, xc_vs_c)
    pots(35) = pot('BLYP' , xc_becke88, xc_lyp_c )

    SCF      = pot('scf')
    do xc=1,xc_NXC
       SCF%contrib(xc) = xc_cntrl_frac(xc)
    enddo
    pots(N_POT) = SCF

    ! Count the number of various functionals:
    call compute (counter)

    if (counter == 0) then
       allocate(a_xc_arr(1), name_xc_arr(1), stat=alloc_stat)
       ASSERT(alloc_stat==0)

       a_xc_arr(1) = 0.0
       name_xc_arr(1) = 'no functional'
    else
       allocate(a_xc_arr(counter), name_xc_arr(counter), stat=alloc_stat)
       ASSERT(alloc_stat==0)

       ! Fill the energies and names for various functionals:
       call compute (counter, a_xc_arr, name_xc_arr)
    endif

  contains

    subroutine compute (counter, a_xc_arr, name_xc_arr)
      implicit none
      integer(i4_kind), intent(out)   :: counter
      real(r8_kind), intent(out), optional :: a_xc_arr(:)
      character(len=32), intent(out),  optional :: name_xc_arr(:)
      ! *** end of interface ***

      real(r8_kind), parameter :: zero = 0.0_r8_kind
      integer(i4_kind) :: xc, xc_contrib
      logical :: have_all
      real(r8_kind) :: exc

      counter = 0

      do xc = 1, N_POT
         ! check if have all pre-requisites
         have_all = .true.
         do xc_contrib = 1, xc_NXC
            if (.not. abs(pots(xc)%contrib(xc_contrib)) > 0.0) cycle
            have_all = have_all .and. frac(xc_contrib) /= 0.0
         enddo

         if (have_all) then
            ! this one may be computed
            counter = counter + 1

            ! set the name of th XC functional:
            if (present(name_xc_arr)) then
               name_xc_arr(counter) = pots(xc)%name
            endif

            ! actually compute the net energy:
            if (present(a_xc_arr)) then
               exc = zero
               do xc_contrib = 1, xc_NXC
                  if (abs(pots(xc)%contrib(xc_contrib)) > 0.0) then
                     exc = exc + energy(xc_contrib) * pots(xc)%contrib(xc_contrib)
                  endif
               enddo
               a_xc_arr(counter) = exc
               DPRINT pots(xc)%name,' = ', pots(xc)%contrib
               DPRINT pots(xc)%name,' = ', exc
            endif
         endif
      enddo
    end subroutine compute

  end subroutine xc_func_get_exc

  subroutine xc_func_reduce()
    use comm, only: comm_reduce
    implicit none
    ! *** end of interface ***

    call comm_reduce(energy)
  end subroutine xc_func_reduce

#ifdef WITH_GUILE
  function qm_xc (xc, ra, rb, gaa, gbb, gab) result (f) bind(c)
    !
    ! Export XC functionals to the  Scheme world. Try to keep it close
    ! to what really used in PG.
    !
    use iso_c_binding!, only: bind(c)
    use scm, only: scm_t, scm_from, scm_to_double, &
         scm_to_stringbuf, scm_list
    use scheme, only: scheme_make_list
    use vwnc, only: xalpha_calcMDA, vwn_calcMDA, vwn_ldac
    use pw_ldac_module, only: pw_ldac
    use pbe_ggcxc_module, only: pbe_ggcc
    use becke_perdew_module, only: becke88_calcMDA
    use exchange, only: exchange_lda, exchange_gga, &
         X_XALPHA, X_BECKE88
    implicit none
    type(scm_t), intent(in), value :: xc, ra, rb, gaa, gbb, gab
    type(scm_t) :: f
    ! *** end of interface ***

    character(len=32) :: sxc
    integer :: slen
    type(scm_t) :: f1, f2

    integer, parameter :: vla = 1, ispin = 2

    !
    ! VN  variant  (will  be  implemented  and  used  everywhere,  see
    ! e.g. exchange.f90).  NOTE for second derivatives packing:
    !
    !              1    2    3    4    5    6
    ! dfdg        AA   BB   AB
    ! dfdndn      AA   BB   AB
    ! dfdndg     AAA  BBB  BAA  ABB  AAB  BAB
    ! dfdgdg    AAAA BBBB AABB AAAB BBAB ABAB
    !
    enum, bind(c)
       enumerator :: A = 1, B
       enumerator :: AA = 1, BB, AB
       enumerator :: AAA = 1, BBB,  BAA,  ABB,  AAB,  BAB
       enumerator :: AAAA = 1, BBBB, AABB, AAAB, BBAB, ABAB
    end enum

    real(r8_kind) :: rho(vla, ispin)
    real(r8_kind) :: gam(vla, 2 * ispin - 1)

    real(r8_kind) :: x(vla)

    ! 5 first derivatives, historically distributed over two arrays:
    real(r8_kind) :: xr(vla, ispin)
    real(r8_kind) :: xg(vla, 2 * ispin - 1)

    ! 15  = 5  *  (5 +  1)  / 2  of  second derivatives,  historically
    ! distributed over three arrays:
    real(r8_kind) :: xrr(vla, 2 * ispin - 1) ! 3
    real(r8_kind) :: xrg(vla, ispin * (2 * ispin - 1)) ! 6
    real(r8_kind) :: xgg(vla, 6) ! FIXME!
    real(r8_kind) :: x2(15)      ! all off them ...


    ! FIXME: pass args as a list (of variable length?):
    rho(1, A) = scm_to_double (ra)
    rho(1, B) = scm_to_double (rb)

    gam(1, AA) = scm_to_double (gaa)
    gam(1, BB) = scm_to_double (gbb) ! sic
    gam(1, AB) = scm_to_double (gab) ! sic

    ! Some  of the  subroutines may  increment the  result  instead of
    ! setting it:
    x = 0.0

    xr = 0.0
    xg = 0.0

    xrr = 0.0
    xrg = 0.0
    xgg = 0.0

    ! The  output string  is only  complete if  length <=  len(buf) on
    ! output:
    call scm_to_stringbuf (xc, sxc, slen)
    ASSERT(slen<=len(sxc))

    select case (sxc)

    case ("xalpha")
       call exchange_lda (X_XALPHA, vla, ispin, rho, x, xr, xrr)

    case ("xalpha-old")
       call xalpha_calcMDA (rho, xr, ispin, x, vla, dvdrho=xrr)

    case ("vwn")
       call vwn_ldac (vla, ispin, rho, x, xr, xrr)

    case ("vwn-old")
       call vwn_calcMDA (rho, xr, ispin, x, vla, dvdrho=xrr)

    case ("pwlda")
       ! FIXME:  before it  was called  with ispin  + 2,  that is
       ! either with PWLDAC_RESP_R or with PWLDAC_RESP_R:
       call pw_ldac (vla, ispin, rho, x, xr, xrr)

    case ("becke88")
       call exchange_gga (X_BECKE88, vla, ispin, rho, gam, x, xr, xg, xrr, xrg, xgg)

    case ("becke88-old")
       call becke88_calcMDA (rho, gam, xr, ispin, x, xg, vla, &
            df_drhodrho=xrr, &
            df_drhodgamma=xrg, &
            df_dgammadgamma=xgg)

    case ("pbec")
       call pbe_ggcc ("pbe", vla, ispin, rho, gam, x, xr, xg, xrr, xrg, xgg)

    case ("pbecsol")
       call pbe_ggcc ("pbesol", vla, ispin, rho, gam, x, xr, xg, xrr, xrg, xgg)

    case default
       ABORT("no such XC: "//sxc)
    end select

    ! Store the value here (temporarily):
    f = scm_from (x (1))

    ! First derivatives in the same order as the arguments:
    f1 = scm_list (scm_from (xr (1, A)),  & ! d / dra
                   scm_from (xr (1, B)),  & ! d / drb
                   scm_from (xg (1, AA)), & ! d / dgaa
                   scm_from (xg (1, BB)), & ! d / dgbb
                   scm_from (xg (1, AB)))   ! d / dgab

    ! It appears to be more  convenient to flatten the array of second
    ! derivatives first:
    x2 (1) = xrr (1, AA)        ! d2 / dra dra
    x2 (2) = xrr (1, AB)        ! d2 / dra drb
    x2 (3) = xrg (1, AAA)       ! d2 / dra dgaa
    x2 (4) = xrg (1, ABB)       ! d2 / dra dgbb
    x2 (5) = xrg (1, AAB)       ! d2 / dra dgab

    x2 (6) = xrr (1, BB)        ! d2 / drb drb
    x2 (7) = xrg (1, BAA)       ! d2 / drb dgaa
    x2 (8) = xrg (1, BBB)       ! d2 / drb dgbb
    x2 (9) = xrg (1, BAB)       ! d2 / drb dgab

    x2(10) = xgg (1, AAAA)      ! d2 / dgaa dgaa
    x2(11) = xgg (1, AABB)      ! d2 / dgaa dgbb
    x2(12) = xgg (1, AAAB)      ! d2 / dgaa dgab

    x2(13) = xgg (1, BBBB)      ! d2 / dgbb dgbb
    x2(14) = xgg (1, BBAB)      ! d2 / dgbb dgab

    x2(15) = xgg (1, ABAB)      ! d2 / dgab dgab

    ! This is  an auxilary func  to convert an  array of doubles  to a
    ! list:
    f2 = scheme_make_list (x2)

    ! Return  a  list  with  the   value  in  CAR  positions  and  the
    ! list-of-lists of derivatives in CDR:
    f = scm_list (f, f1, f2)
  end function qm_xc
#endif


  !--------------- End of module ----------------------------------
end module xc_func

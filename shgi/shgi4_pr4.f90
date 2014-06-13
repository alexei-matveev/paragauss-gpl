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
module shgi4_pr4
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

!# define FPP_TIMERS 2
# include "def.h"
  use type_module, only: &
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: SHR_PR4v ! for "poor man" 4-center ints

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine SHR_PR4v(NAB, NCD, LA, LB, LC, LD, F, SAB, SCD, SFS)
    !
    ! Couples 4-derivatives
    !
    !    F(lma,lmb,lmc,lmd)
    !
    ! with
    !
    !    SAB(lma,lmb) and SCD(lmc,lmd)
    !
    ! by nesting the product rules for solid harmonics:
    !
    !    FS(m) += PR(lm,lm1,lm2) * F(lm1) * S(lm2)
    !
    implicit none
    integer(IK), intent(in)  :: NAB,NCD
    integer(IK), intent(in)  :: LA,LB,LC,LD
    real(RK)   , intent(in)  :: F(NAB,NCD,(LA+1)**2,(LB+1)**2,(LC+1)**2,(LD+1)**2)
    real(RK)   , intent(in)  :: SAB(NAB  ,(LA+1)**2,(LB+1)**2)
    real(RK)   , intent(in)  :: SCD(NCD  ,(LC+1)**2,(LD+1)**2)
    real(RK)   , intent(out) :: SFS(NAB,NCD,2*LA+1,2*LB+1,2*LC+1,2*LD+1)
    ! NOTE: NABCD = NAB * NCD
    ! *** end of interface ***

    ! not sure if switching between two alternatives
    ! makes any improvement (or even makes it worse):
    if( (LC+1) * (LD+1) >= (LA+1) * (LB+1) )then
      ! PR over CD-indices first:
      SFS = AB(NAB, NCD, LA, LB, CD(NAB, NCD, LC, LD, F, SCD), SAB)
    else
      ! PR over AB-indices first:
      SFS = CD(NAB, NCD, LC, LD, AB(NAB, NCD, LA, LB, F, SAB), SCD)
    endif

  end subroutine SHR_PR4v

  function CD(NAB, NCD, LC, LD, F, S) result(FS)
    use solhrules_module, only: solhrules_product
    implicit none
    integer(IK), intent(in) :: NAB, NCD, LC, LD
    real(RK), intent(in) :: F(:, :, :, :, :, :)
    real(RK), intent(in) :: S(:, :, :)
    real(RK)             :: FS(NAB, NCD, size(F,3), size(F,4), 2*LC+1, 2*LD+1)
    ! *** end of interface ***

    integer(IK) :: mc,lmc,tc,lmc1,lmc2
    integer(IK) :: md,lmd,td,lmd1,lmd2
    real(RK)    :: cfc,cfd
    integer(IK) :: a,b
    integer(IK) :: i,j

    do md=1, 2*LD+1
       lmd = LD**2 + md

    do mc=1, 2*LC+1
       lmc = LC**2 + mc

    do b=1, size(F, 4)

    do a=1, size(F, 3)

    ! init (a, b, mc, md) integral:
    FS(:, :, a, b, mc, md) = 0.0

    do  td=1, solhrules_product(lmd)%n_summands
       cfd  = solhrules_product(lmd)%coef  (td)
       lmd1 = solhrules_product(lmd)%lm_sh1(td)
       lmd2 = solhrules_product(lmd)%lm_sh2(td)

    do  tc=1, solhrules_product(lmc)%n_summands
       cfc  = solhrules_product(lmc)%coef  (tc)
       lmc1 = solhrules_product(lmc)%lm_sh1(tc)
       lmc2 = solhrules_product(lmc)%lm_sh2(tc)

       do j=1, NCD
       do i=1, NAB
          FS(i, j, a, b, mc, md) = FS(i, j, a, b,  mc,   md ) &
                                  + F(i, j, a, b, lmc1, lmd1) &
                                  * S(   j,       lmc2, lmd2) &
                                  * cfc * cfd
       enddo ! i
       enddo ! j

    enddo ! tc
    enddo ! td

    enddo ! a
    enddo ! b
    enddo ! mc
    enddo ! md
  end function CD

  function AB(NAB, NCD, LA, LB, F, S) result(FS)
    use solhrules_module, only: solhrules_product
    implicit none
    integer(IK), intent(in) :: NAB, NCD, LA, LB
    real(RK), intent(in) :: F(:, :, :, :, :, :)
    real(RK), intent(in) :: S(:, :, :)
    real(RK)             :: FS(NAB, NCD, 2*LA+1, 2*LB+1, size(F,5), size(F,6))
    ! *** end of interface ***

    integer(IK) :: ma,lma,ta,lma1,lma2
    integer(IK) :: mb,lmb,tb,lmb1,lmb2
    real(RK)    :: cfa,cfb
    integer(IK) :: c,d
    integer(IK) :: i,j

    do d=1, size(F, 6)

    do c=1, size(F, 5)

    do mb=1, 2*LB+1
       lmb = LB**2 + mb

    do ma=1, 2*LA+1
       lma = LA**2 + ma

    ! init (ma, mb, c, d) integral:
    FS(:, :, ma, mb, c, d) = 0.0

    do  tb=1, solhrules_product(lmb)%n_summands
       cfb  = solhrules_product(lmb)%coef  (tb)
       lmb1 = solhrules_product(lmb)%lm_sh1(tb)
       lmb2 = solhrules_product(lmb)%lm_sh2(tb)

    do  ta=1, solhrules_product(lma)%n_summands
       cfa  = solhrules_product(lma)%coef  (ta)
       lma1 = solhrules_product(lma)%lm_sh1(ta)
       lma2 = solhrules_product(lma)%lm_sh2(ta)

       do j=1, NCD
       do i=1, NAB
          FS(i, j, ma, mb, c, d) = FS(i, j,  ma,   mb,  c, d) &
                                  + F(i, j, lma1, lmb1, c, d) &
                                  * S(i,    lma2, lmb2)       &
                                  * cfa * cfb
       enddo ! i
       enddo ! j

    enddo ! ta
    enddo ! tb

    enddo ! ma
    enddo ! mb
    enddo ! c
    enddo ! d
    end function AB


  !--------------- End of module ----------------------------------
end module shgi4_pr4

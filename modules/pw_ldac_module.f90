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
!=====================================================================
! Public interface of module
!=====================================================================
module pw_ldac_module
!---------------------------------------------------------------
!
!  Purpose: The module calculates local part c-functional
!           of Perdew and Wang
!           The calculation follows the scheme of
!  Reference: J.P. Perdew, Y. Wang, Phys. Rev. B, 45, 13244
!             (1992)
!
!           The used variables are as follows:
!           n:  density
!           G: squared densities` gradient
!           dFdn: derivative of Fx(c) with respect to n
!           Fc:    functionals
!
!  Module called by: xc_hamiltonian
!
!  Author: AM
!  Date: 10.98
!
!== Interrupt of public interface of module =====================
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: SB
! Date:   10.04
! Description: For TDDFT
!  Added dFdndn
!              spin=1:
!                dFdn(:,1) = dv_a/drho_a
!                dFdn(:,2) = dv_a/drho_b
!              spin=2:
!                dFdn(:,1) = dv_a/drho_a
!                dFdn(:,2) = dv_a/drho_b
!                dFdn(:,3) = dv_b/drho_b
!
!------------ Modules used --------------------------------------
!---------------------------------------------------------------------
! Modifications
!---------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!------------ Modules used --------------------------------------

# include "def.h"
  use type_module, RK=>r8_kind, IK=>i4_kind
  implicit none
  private
  save
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
  public :: pw_ldac

  public :: set_param
  ! what-to-do keys:
  integer(IK), parameter, public :: &
       do_ECU = 1, &
       do_ECP = 2, &
       do_mALF= 3

  integer(IK), parameter, PUBLIC :: &
       PWLDAC_RESTR     =  1, &   ! compat with ispin == 1
       PWLDAC_UNRESTR   =  2, &   ! compat with ispin == 2
       PWLDAC_ECU       = -do_ECU, & ! same meaning as ispin == 1
       PWLDAC_ECP       = -do_ECP, & ! fully polarized: rho == rhoa
       PWLDAC_RESP_R    =  3, & ! response, spin = 1
       PWLDAC_RESP_U    =  4 ! response, spin = 2

  !===================================================================
  ! End of public interface of module
  !===================================================================

  real(RK), parameter:: &
       ZERO  = 0.0_RK ,&
       HALF  = 0.5_RK ,&
       ONE   = 1.0_RK ,&
       TWO   = 2.0_RK ,&
       THREE = 3.0_RK ,&
       FOUR  = 4.0_RK ,&
       FIVE  = 5.0_RK ,&
       SIX   = 6.0_RK ,&
       SEVEN = 7.0_RK ,&
       EIGHT = 8.0_RK, &
       NINE  = 9.0_RK ,&
       TWELVE=12.0_RK, &
       DEPS  = 1.0e-30_RK ,&
       pi = 3.14159265358979324_RK

!+++++++++++++++++++++++++++++++

contains

  pure function c_lda(na, nb, RPA) result(Fc)
    !
    ! FIXME: code duplication, see pbe_ggcxc_module.f90
    !
    use ad2x2
    implicit none
    type(ad), intent(in) :: na, nb ! electron density
    logical, intent(in) :: RPA
    type(ad) :: Fc ! c-functional
    ! *** end of interface ***

    real(RK), parameter :: c1d3 = one/three, c4d3 = four/three
    real(RK), parameter :: crs  = (three / (four * pi))**c1d3
    real(RK), parameter :: cf   = one / (two**c4d3 - two)
    real(RK), parameter :: cfzz = two * c4d3 * c1d3 * cf

    type(ad) :: nn, rs, ECU, ECP, mALF, z, z2, z4, f, EC

    ! general case, na /= nb
    nn = na + nb

    ! historically:
    if (nn < DEPS) nn = fix(deps)

    rs = crs / nn**c1d3

    ECU = g(rs, do_ECU, RPA)
    ECP = g(rs, do_ECP, RPA)
    mALF = g(rs, do_mALF, RPA)

    z = (na - nb) / nn
    z2 = z * z
    z4 = z2 * z2

    f = cf * ((one + z + DEPS)**c4d3 + (one - z + DEPS)**c4d3 - two)

    ! mALF == - ALFA:
    EC = ECU + (ECP - ECU) * f * z4 - mALF * f * (one - z4) / cfzz

    Fc = nn * EC
  contains

    pure function g(rs, do_WHAT, RPA) result(E)
      !
      ! FIXME: code duplication, see pbe_ggcxc_module.f90
      !
      implicit none
      type(ad), intent(in) :: rs ! Seitz radius
      integer(IK), intent(in) :: do_WHAT
      logical, intent(in) :: RPA ! together: which of those three to do ?
      type(ad) :: E ! either ec(rs,0), ec(rs,1), -alpha(rs)
      ! *** end of interface ***

      real(RK) :: a, a1, b1, b2, b3, b4, p
      type(ad) :: rs1d2, rs3d2, rsp, Q0, Q1, Q2

      call set_param(a, a1, b1, b2, b3, b4, p, do_WHAT, RPA)

      rs1d2 = rs**half
      rs3d2 = rs * rs1d2
      rsp = rs**p

      Q0 = - two * a * (one + a1 * rs)
      Q1 = two * a * (b1 * rs1d2 + b2 * rs + b3 * rs3d2 + b4 * rs * rsp)
      Q2 = log(one + one / Q1)
      E = Q0 * Q2
    end function g
  end function c_lda

#if 1
  subroutine pw_ldac(nv, iop, n, Fc, dFdn, dFdndn, iRPA)
    use ad2x2, only: ad, var, fix, val, fst, snd, operator(*)
    implicit none
    integer(IK), intent(in) :: nv ! vector length
    integer(IK), intent(in) :: iop ! See PWDLAC_* keys
    real(RK), dimension(:,:), intent(in) :: n ! electron density n(nv,1) or n(nv,2)
    real(RK), dimension(:), intent(inout) :: Fc ! c-functional Fc(nv)
    real(RK), dimension(:,:), intent(inout) :: dFdn !dF/dn dFdn(nv,1) or dFdn(nv,2)
    real(RK), dimension(:,:), intent(inout), optional :: dFdndn
    !dF/dnn (spin = 1) dFdn(nv,1), dFdn(nv,2)  or
    !       (spin = 2) dFdn(nv,i), i = 1..4
    logical, optional, intent(in)            :: iRPA
    ! Random Phase Approximation switch: true == on
    ! *** end of interface ***

    !
    ! VN  variant  (will  be  implemented  and  used  everywhere,  see
    ! e.g. exchange.f90).  NOTE for second derivatives packing:
    !
    !              1    2    3
    ! dfdg        AA   BB   AB
    ! dfdndn      AA   BB   AB
    !
    enum, bind(c)
       enumerator :: A = 1, B
       enumerator :: AA = 1, BB, AB
    end enum

    type(ad) :: f, na, nb, nn
    logical :: RPA
    integer :: i

    RPA = .false.
    if(present(iRPA)) RPA = iRPA

    select case(iop)

    case (PWLDAC_UNRESTR)
       ! ispin == 2 general case, rhoa /= rhob as two independent
       ! variables:
       do i = 1, nv
          na = var(1, n(i, A))
          nb = var(2, n(i, B))

          f = c_lda(na, nb, RPA)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, A) = dFdn(i, A) + fst(1, f)
          dFdn(i, B) = dFdn(i, B) + fst(2, f)

          if (present(dFdndn)) then
             dFdndn(i, AA) = dFdndn(i, AA) + snd(1, 1, f)
             dFdndn(i, BB) = dFdndn(i, BB) + snd(2, 2, f)
             dFdndn(i, AB) = dFdndn(i, AB) + snd(1, 2, f)
          endif
       enddo

    case (PWLDAC_ECU, PWLDAC_RESTR)
       ! ispin == 1 rhoa = rhob = rho/2, rho as a single independent
       ! variable:
       do i = 1, nv
          nn = var(1, n(i, 1))

          f = c_lda(half * nn, half * nn, RPA)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, 1) = dFdn(i, 1) + fst(1, f)

          if (present(dFdndn)) then
             dFdndn(i, 1) = dFdndn(i, 1) + snd(1, 1, f)
          endif
       enddo

    case (PWLDAC_ECP)
       !
       ! No corresponding ispin. Fully polarized case, rho == rhoa, as
       ! the only independent variable.  FIXME: This is used from HCTH
       ! code exclusively.
       !
       do i = 1, nv
          na = var(1, n(i, 1))
          nb = fix(zero)

          f = c_lda(na, nb, RPA)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, 1) = dFdn(i, 1) + fst(1, f)

          if (present(dFdndn)) then
             dFdndn(i, 1) = dFdndn(i, 1) + snd(1, 1, f)
          endif
       enddo

    case (PWLDAC_RESP_R) ! response case dVdnadna and dVdnadnb
       ! ispin = 1 => rhoa = rhob = rho/2
       ! rhoa = rhob = rho/E2
       !
       ! FIXME: is one supposed to  treat rhoa and rhob as independent
       ! variables  here?  I  vaguely remember  something  about mixed
       ! second derivative  needed also for  (magnetic?)  solutions of
       ! TDDFT both in the restricted and unrestricted cases.
       !
       ABORT("needs work!")

    case (PWLDAC_RESP_U) ! ispin = 2
       ! FIXME: how is it different from PWLDAC_UNRESTR?
       ABORT("needs work!")

    case default
       print *,'pw_ldac: no such case: iop=',iop
       stop
    end select
  end subroutine pw_ldac
#else
  subroutine pw_ldac(nv,iop,n,Fc,dFdn,dFdndn,iRPA)
    implicit none
    integer(IK), intent(in)                  :: nv
    ! vector length
    integer(IK), intent(in)                  :: iop
    ! See PWDLAC_* keys
    real(RK), dimension(:,:), intent(in)     :: n
    !electron density n(nv,1) or n(nv,2)
    real(RK), dimension(:), intent(inout)    :: Fc
    !c-functional Fc(nv)
    real(RK), dimension(:,:), intent(inout)  :: dFdn
    !dF/dn dFdn(nv,1) or dFdn(nv,2)
    real(RK), dimension(:,:), intent(inout),optional  :: dFdndn
    !dF/dnn (spin = 1) dFdn(nv,1), dFdn(nv,2)  or
    !       (spin = 2) dFdn(nv,i), i = 1..4
    logical, optional, intent(in)            :: iRPA
    ! Random Phase Approximation switch: true == on
    ! *** end of interface ***

    real(RK), parameter :: &
         c1d3 = one/three, &
         c4d3 = four/three
    real(RK) :: crs,cf,cfzz
    real(RK), dimension(nv) :: nn,rs,z,z2,z3,z4,f,fz,EC,ECU,ECP,mALF,&
         & dECdrs,dECdrsdrs,dECUdrs, dECUdrsdrs, dECPdrs, dECPdrsdrs,&
         & dmALFdrs,dmALFdrsdrs,dECdz,dECdzdz,vcomm,drsdnn,drsdnndnn,&
         & dzda,dzdada,dzdb,dzdadb,fzz,dECdrsdz,dzdbdb
!!!    real(RK) :: NUM_DIFF
    logical :: RPA

    RPA = .false.
    if(present(iRPA)) RPA = iRPA

    crs  = (three/(four*pi))**c1d3
    cf   = one/(two**c4d3 - two)
    cfzz = two*c4d3*c1d3*cf

    select case(iop)

    case (PWLDAC_UNRESTR) ! ispin == 2
       ! general case, rhoa /= rhob
       nn = n(1:nv,1) + n(1:nv,2)
       where(nn<deps) nn = deps

       rs = crs/nn**c1d3

       call g(rs,ECU,dECUdrs,  dECUdrsdrs,nv, do_ECU, RPA)
       call g(rs,ECP,dECPdrs,  dECPdrsdrs,nv, do_ECP, RPA)
       call g(rs,mALF,dmALFdrs,dmALFdrsdrs,nv,do_mALF,RPA)

       z = (n(1:nv,1)-n(1:nv,2))/nn
       z2 = z*z
       z4 = z2*z2

       f = cf*((one+z+deps)**c4d3+(one-z+deps)**c4d3 - two)

       EC = ECU+(ECP-ECU)*f*z4 - mALF*f*(one-z4)/cfzz
       ! mALF == - ALFA !!!
       Fc(1:nv) = Fc(1:nv) + nn*EC
       ! end calculation of Fc

       ! now the potential:
       dECdrs = dECUdrs+(dECPdrs-dECUdrs)*f*z4 - dmALFdrs*f*(one-z4)/cfzz
       fz = cf*c4d3*((one+z+deps)**c1d3 - (one-z+deps)**c1d3)
       dECdz = (four*(z*z2)*f + z4*fz)*(ECP-ECU+mALF/cfzz) - fz*mALF/cfzz

       vcomm = EC - c1d3*rs*dECdrs - z*dECdz
       dFdn(1:nv,1) = dFdn(1:nv,1) + vcomm + dECdz
       dFdn(1:nv,2) = dFdn(1:nv,2) + vcomm - dECdz
       ! done with ppotential

    case (PWLDAC_ECU, PWLDAC_RESTR) ! ispin == 1
       ! rhoa = rhob = rho/2
       nn = n(1:nv,1)
       where(nn<deps) nn = deps

       rs = crs/nn**c1d3

       call g(rs,ECU,dECUdrs,  dECUdrsdrs,nv,do_ECU, RPA)

       Fc(1:nv) = Fc(1:nv) + nn*ECU
       ! done with Fc
       dFdn(1:nv,1) = dFdn(1:nv,1) + ECU - c1d3*rs*dECUdrs
       ! done with potential

    case (PWLDAC_ECP) ! no corresponding ispin
       ! fully polarized case, rho == rhoa
       nn = n(1:nv,1)
       where(nn<deps) nn = deps

       rs = crs/nn**c1d3

       call g(rs,ECP,dECPdrs,dECPdrsdrs,nv,do_ECP,RPA)

       Fc(1:nv) = Fc(1:nv) + nn*ECP
       ! done with Fc
       dFdn(1:nv,1) = dFdn(1:nv,1) + ECP - c1d3*rs*dECPdrs
       ! done with potential

    case (PWLDAC_RESP_R) ! response case dVdnadna and dVdnadnb
       ! ispin = 1 => rhoa = rhob = rho/2
       ! rhoa = rhob = rho/2
       nn = n(1:nv,1)
       where(nn<deps) nn = deps

       rs = crs/nn**c1d3
       drsdnn = -c1d3*rs/nn
       drsdnndnn  = c1d3*c4d3*rs/nn**two

       call g(rs,ECU,dECUdrs,  dECUdrsdrs,nv,do_ECU, RPA)
       call g(rs,ECP,dECPdrs,  dECPdrsdrs,nv,do_ECP, RPA)
       call g(rs,mALF,dmALFdrs,dmALFdrsdrs,nv,do_mALF,RPA)

       Fc(1:nv) = Fc(1:nv) + nn*ECU
       ! done with Fc
       dFdn(1:nv,1) = dFdn(1:nv,1) + ECU - c1d3*rs*dECUdrs
       ! done with potential

       !! SB: dVdrho_a
       z = zero
       dzda =  one/nn-z/nn
       dzdada = -two/(nn*nn)+two*z/(nn*nn)
       dzdb = -one/nn-z/nn
       dzdadb = two*z/(nn*nn)
       z2 = z*z
       z3 = z*z*z
       z4 = z2*z2
       f   = cf*((one+z+deps)**c4d3+(one-z+deps)**c4d3 - two)
       fz  = cf*c4d3*((one+z+deps)**c1d3 - (one-z+deps)**c1d3)
       fzz = cf*c4d3*c1d3*((one+z+deps)**(c1d3-one)+(one-z+deps)**(c1d3-one))

       dECdrs(1:nv)  = dECUdrs+(-dmALFdrs)*f*(one-z4)/cfzz+(dECPdrs-dECUdrs)*f*z4
       dECdz(1:nv)   = (-mALF)*fz*(one-z4)/cfzz-four*(-mALF)*f*z3/cfzz+(ECP-ECU)*fz*z4+four*(ECP-&
            & ECU)*f*z3

       dECdrsdrs(1:nv) = dECUdrsdrs-dmALFdrsdrs*f*(one - z4)/cfzz+ (dECPdrsdrs - dECUdrsdrs)*f*z4

       dECdrsdz(1:nv)  = -dmALFdrs*fz*(one-z4)/cfzz-four*(-dmALFdrs)*f*z3/cfzz+(dECPdrs-dECUdrs)*&
            & fz*z4+four*(dECPdrs-dECUdrs)*fz*z3

       dECdzdz(1:nv)   = -mALF*fzz*(one - z4)/cfzz- eight*(-mALF)*fz*z3/cfzz-twelve*(-mALF)*f*z2/&
            & cfzz+(ECP-ECU)*fzz*z4+eight*(ECP-ECU)*fz*z3+twelve*(ECP-ECU)*f*z2

       dFdndn(1:nv,1) = dFdndn(1:nv,1)&
            +(two*dECdrs*drsdnn+nn*(dECdrsdrs*drsdnn*drsdnn + dECdrs*drsdnndnn)&
            +two*nn*(dECdrsdz  *dzda*drsdnn)&
            +nn*dECdzdz*dzda*dzda) ! dVdnadna

       dFdndn(1:nv,2) = dFdndn(1:nv,2)&
            +(two*dECdrs*drsdnn+nn*(dECdrsdrs*drsdnn*drsdnn + dECdrs*drsdnndnn)&
            +nn*dECdzdz*dzda*dzdb&
            -two*z*dECdrsdz*drsdnn) ! dVdnadnb



    case (PWLDAC_RESP_U) ! ispin = 2
       nn = n(1:nv,1) + n(1:nv,2)
       where(nn<deps) nn = deps

       rs = crs/nn**c1d3
       drsdnn = -c1d3*rs/nn
       drsdnndnn  = c1d3*c4d3*rs/nn**two

       call g(rs,ECU,dECUdrs,  dECUdrsdrs,nv,do_ECU, RPA)
       call g(rs,ECP,dECPdrs,  dECPdrsdrs,nv,do_ECP, RPA)
       call g(rs,mALF,dmALFdrs,dmALFdrsdrs,nv,do_mALF,RPA)

       z = (n(1:nv,1)-n(1:nv,2))/nn
       z2 = z*z
       z4 = z2*z2

       f = cf*((one+z+deps)**c4d3+(one-z+deps)**c4d3 - two)

       EC = ECU+(ECP-ECU)*f*z4 - mALF*f*(one-z4)/cfzz
       ! mALF == - ALFA !!!
       Fc(1:nv) = Fc(1:nv) + nn*EC
       ! end calculation of Fc

       ! now the potential:
       dECdrs = dECUdrs+(dECPdrs-dECUdrs)*f*z4 - dmALFdrs*f*(one-z4)/cfzz
       fz = cf*c4d3*((one+z+deps)**c1d3 - (one-z+deps)**c1d3)
       dECdz = (four*(z*z2)*f + z4*fz)*(ECP-ECU+mALF/cfzz) - fz*mALF/cfzz

       vcomm = EC - c1d3*rs*dECdrs - z*dECdz
       dFdn(1:nv,1) = dFdn(1:nv,1) + vcomm + dECdz
       dFdn(1:nv,2) = dFdn(1:nv,2) + vcomm - dECdz
       ! done with ppotential

       !! SB: dVdrho_a
       z = (n(1:nv,1)-n(1:nv,2))/nn
       dzda =  one/nn-z/nn
       dzdada = -two/(nn*nn)+two*z/(nn*nn)
       dzdb = -one/nn-z/nn
       dzdadb = two*z/(nn*nn)
       dzdbdb = two/(nn*nn)+two*z/(nn*nn)
       z3 = z*z*z
       fzz = cf*c4d3*c1d3*((one+z+deps)**(c1d3-one)+(one-z+deps)**(c1d3-one))

       dECdrsdrs(1:nv) = dECUdrsdrs-dmALFdrsdrs*f*(one - z4)/cfzz + (dECPdrsdrs - dECUdrsdrs)*f*z4

       dECdrsdz(1:nv)  = -dmALFdrs*fz*(one-z4)/cfzz-four*(-dmALFdrs)*f*z3/cfzz+(dECPdrs-dECUdrs)*&
            & fz*z4+four*(dECPdrs-dECUdrs)*f*z3

       dECdzdz(1:nv)   = -mALF*fzz*(one - z4)/cfzz - eight*(-mALF)*fz*z3/cfzz-twelve*(-mALF)*f*z2/&
            & cfzz+(ECP-ECU)*fzz*z4+eight*(ECP-ECU)*fz*z3+twelve*(ECP-ECU)*f*z2

       dFdndn(1:nv,1) = dFdndn(1:nv,1)&
            +(two*dECdrs*drsdnn+nn*(dECdrsdrs*drsdnn*drsdnn + dECdrs*drsdnndnn)&
            +two*nn*(dECdrsdz  *dzda*drsdnn)&
            +nn*dECdzdz*dzda*dzda)  ! dVdnadna

       dFdndn(1:nv,2) = dFdndn(1:nv,2)&
            +(two*dECdrs*drsdnn+nn*(dECdrsdrs*drsdnn*drsdnn + dECdrs*drsdnndnn)&
            +two*nn*(dECdrsdz  *dzdb*drsdnn)&
            +nn*dECdzdz*dzdb*dzdb) ! dVdnbdnb

       dFdndn(1:nv,3) = dFdndn(1:nv,3)&
            +(two*dECdrs*drsdnn+nn*(dECdrsdrs*drsdnn*drsdnn + dECdrs*drsdnndnn)&
            +nn*dECdzdz*dzda*dzdb&
            -two*z*dECdrsdz*drsdnn) ! dVdnadnb

!!       ! same as dFdndn(1:nv,2)
!!       dFdndn(1:nv,4) = dFdndn(1:nv,4)&
!!            +(two*dECdrs*drsdnn+nn*(dECdrsdrs*drsdnn*drsdnn + dECdrs*drsdnndnn)&
!!            +nn*dECdzdz*dzda*dzdb&
!!            -two*z*dECdrsdz*drsdnn)! dVdnbdna
    case default
       print *,'pw_ldac: no such case: iop=',iop
       stop
    end select
  end subroutine pw_ldac

  subroutine g(rs,E,dEdrs,dEdrsdrs,nv,do_WHAT,RPA)
    ! dEdrsdrs - SB implementation
    integer(IK), intent(in)                :: nv
    ! dimension
    real(RK), dimension(nv), intent(in)    :: rs
    ! Seitz radius
    real(RK), dimension(nv), intent(out)   :: E
    ! one of the three: ec(rs,0), ec(rs,1), -alpha(rs)
    real(RK), dimension(nv), intent(out)   :: dEdrs
    ! dE/drs
    real(RK), dimension(nv), intent(out)   :: dEdrsdrs
    ! dE/drsdrs
    integer(IK), intent(in)                :: do_WHAT
    logical, intent(in)                    :: RPA
    ! together: which of those three to do ?
    ! *** end of interface ***

    real(RK) :: p1, a,a1,b1,b2,b3,b4,p
    real(RK), dimension(nv) :: rs1d2,rs3d2,rsp,Q0,Q1,Q2,Q3
    real(RK), dimension(nv) :: Z1,Z2,Z3,Z4

    call set_param(a,a1,b1,b2,b3,b4,p,do_WHAT,RPA)

    rs1d2 = rs**half
    rs3d2 = rs*rs1d2
    rsp = rs**p

    Q0 = -two*a*(one+a1*rs)
    Q1 = two*a*(b1*rs1d2+b2*rs+b3*rs3d2+b4*rs*rsp)
    Q2 = log(one+one/Q1)
    E = Q0*Q2
    ! end calculation of E
    p1 = p+one
    Q3 = a*(b1/rs1d2+two*b2+three*b3*rs1d2+two*b4*p1*rsp)
    dEdrs = -two*a*a1*Q2 - Q0*Q3/(Q1**two+Q1)
    ! end calculation of dE/drs
    Z1 = one+one/Q1
    Z2 = -b1/(four*rs3d2)+three*b3/(four*rs1d2)+b4*rs*rsp*p1*p1/(rs*rs)-b4*rs*rsp*p1/(rs*rs)
    Z3 = b1/(two*rs1d2)+b2+three*b3*rs1d2/two+b4*rs*rsp*p1/rs
    Z4 = b1*rs1d2+b2*rs+b3*rs3d2+b4*rs*rsp
    dEdrsdrs = two*a1*Z3/(Z4*Z4*Z1)-two*(one+a1*rs)*Z3*Z3/(Z4**three*Z1)+(one+a1*rs)*Z2/(Z4*Z4*Z1)+&
               (one+a1*rs)*Z3*Z3/(two*Z4**four*Z1*Z1*A)
  end subroutine g
#endif

  pure subroutine set_param(a, a1, b1, b2, b3, b4, p, do_WHAT, RPA)
    implicit none
    real(RK), intent(out)   :: a,a1,b1,b2,b3,b4,p
    integer(IK), intent(in) :: do_WHAT
    logical, intent(in)     :: RPA
    ! *** end of interface ***

    if( .not. RPA )then
       if     ( do_WHAT == do_ECU )then
          a  = 0.031091_RK
          a1 = 0.21370_RK
          b1 = 7.5957_RK
          b2 = 3.5876_RK
          b3 = 1.6382_RK
          b4 = 0.49294_RK
          p  = 1.0_RK
       else if( do_WHAT == do_ECP )then
          a  = 0.015545_RK
          a1 = 0.20548_RK
          b1 = 14.1189_RK
          b2 = 6.1977_RK
          b3 = 3.3662_RK
          b4 = 0.62517_RK
          p  = 1.0_RK
       else if( do_WHAT == do_mALF )then
          a  = 0.016887_RK
          a1 = 0.11125_RK
          b1 = 10.357_RK
          b2 = 3.6231_RK
          b3 = 0.88026_RK
          b4 = 0.49671_RK
          p  = 1.0_RK
       endif
    else ! i.e. if RPA
      if     ( do_WHAT == do_ECU )then
         a  = 0.031091_RK
         a1 = 0.082477_RK
         b1 = 5.1486_RK
         b2 = 1.6483_RK
         b3 = 0.23647_RK
         b4 = 0.20614_RK
         p  = 0.75_RK
      else if( do_WHAT == do_ECP )then
         a  = 0.015545_RK
         a1 = 0.035374_RK
         b1 = 6.4869_RK
         b2 = 1.3083_RK
         b3 = 0.15180_RK
         b4 = 0.082349_RK
         p  = 0.75_RK
      else if( do_WHAT == do_mALF )then
         a  = 0.016887_RK
         a1 = 0.028829_RK
         b1 = 10.357_RK
         b2 = 3.6231_RK
         b3 = 0.47990_RK
         b4 = 0.12279_RK
         p  = 1.0_RK
      endif
    endif
  end subroutine set_param

end module pw_ldac_module

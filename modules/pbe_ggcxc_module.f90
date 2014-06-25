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
!===============================================================
! Public interface of module
!===============================================================
module pbe_ggcxc_module
!---------------------------------------------------------------
!
!  Purpose: The module calculates the NON-LOCAL part of
!           xc-functionals of Perdew, Burke and Ernzerhof
!           The calculation follows the scheme of
!  Reference: 1) J.P.Perdew, K.Burke and M.Ernzerhof, Phys. Rev. Lett. 77,
!                3865, 1996.
!             2) B.Hammer, L.B.Hansen and J.K.N0rskov
!                to be published
!
!           The used variables are as follows:
!           n:  density
!           G: squared densities` gradient
!           dFdn: derivative of Fx(c) with respect to n
!           Fx,Fc:    functionals
!           dFdg: derivative of Fx(c) with respect to G
!
!           non-local part of xc (and derivatives) is added
!           to [previously set corresponding local term]
!
!  Module called by: xc_hamiltonian
!
!  Author: AM
!  Date: 9.98
!
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
!------------ Modules used --------------------------------------

# include "def.h"
  use type_module, only: IK=>i4_kind, RK=>r8_kind

  implicit none
  private
  save
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
  public :: pbe_ggcc

!================================================================
! End of public interface of module
!================================================================

  real(RK), parameter:: ZERO  = 0.0_RK ,&
       ONE   = 1.0_RK ,&
       TWO   = 2.0_RK ,&
       THREE = 3.0_RK ,&
       FOUR  = 4.0_RK ,&
       FIVE  = 5.0_RK ,&
       SIX   = 6.0_RK ,&
       SEVEN = 7.0_RK ,&
       NINE  = 9.0_RK ,&
       HALF  = 0.5_RK ,&
       DEPS  = 1.0e-30_RK ,&
       pi = 3.14159265358979324_RK

contains

  function c_gga(na, nb, g2, beta) result(Fc)
    use ad2x3
    implicit none
    type(ad), intent(in) :: na, nb  ! electron density
    type(ad), intent(in) :: g2  ! g2 = gaa + gbb + 2 * gab
    real(RK), intent(in) :: beta ! GGA parameter
    type(ad)  :: Fc ! d(c-functional)/dV
    ! *** end of interface ***

    real(RK), parameter :: gamma=0.03109069086965489503494086371_RK
    ! ==(ONE-log(TWO))/(pi**2)

    type(ad) :: n, ELDAC, za, zb, fi, fi2, fi3, A, t2, x, H, rat
    logical, parameter :: RPA = .false.

    n = na + nb

    if (n <= zero) then
       Fc = fix(zero)
    else
       ELDAC = c_lda(na, nb, RPA) / n

       za = na / n
       zb = nb / n

       if (za < deps) za = fix(deps)
       if (zb < deps) zb = fix(deps)

       fi = (za**(two/three) + zb**(two/three)) * (one / two**(one/three))

       fi2 = fi**2
       fi3 = fi2 * fi

       A = (beta / gamma) / (exp(- ELDAC / (gamma * fi3)) - one + deps)

       t2 = g2 / (fi2 * n**(seven/three) + deps) &
            * (pi / (four * four * (three * pi**2)**(one/three)))

       x = A * t2

       rat = (one + x) / (one + x + x**2)

       H = gamma * fi3 * log(one + (beta / gamma) * t2 * rat)
       Fc = n * H
    endif

  contains

    pure function c_lda(na, nb, RPA) result(Fc)
      !
      ! FIXME: code duplication, see pw_ldac_module.f90
      !
      use pw_ldac_module, only: do_ECU, do_ECP, do_mALF
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
    end function c_lda

    pure function g(rs, do_WHAT, RPA) result(E)
      !
      ! FIXME: code duplication, see pw_ldac_module.f90
      !
      use pw_ldac_module, only: set_param
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
  end function c_gga

  subroutine pbe_ggcc(op, vla, ispin, n, G, Fc, dFdn, dFdg, dFdndn, dFdndg, dFdgdg)
    use ad2x3, only: ad, var, val, fst, snd, operator(*)
    implicit none
    character(len=*), intent(in) :: op ! "pbe" or "pbesol"
    integer(IK), intent(in) :: vla, ispin ! actual vector length, ispin = 1 or 2
    real(RK), dimension(:,:), intent(in) :: n  ! electron density
    real(RK), dimension(:,:), intent(in) :: G  ! Gaa,Gbb,Gab
    real(RK), dimension(:), intent(inout) :: Fc ! d(c-functional)/dV
    real(RK), dimension(:,:), intent(inout):: dFdn ! dF / dn
    real(RK), dimension(:,:), intent(inout):: dFdg ! dF / dg
    real(RK), dimension(:,:), intent(inout), optional :: dFdndn ! d2F / dn dn
    real(RK), dimension(:,:), intent(inout), optional :: dFdndg ! d2F / dn dg
    real(RK), dimension(:,:), intent(inout), optional :: dFdgdg ! d2F / dg dg
    ! *** end of interface ***

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

    type(ad) :: na, nb, nn, g2, f
    integer :: i

    real(RK), parameter :: beta_pbe = 0.066725_RK
    !
    ! compare with             beta = 0.06672455060314922d0
    !
    ! in e.g. http://dft.uci.edu/pubs/PBEsol.html
    !
    real(RK), parameter :: beta_pbesol = 0.046_RK
    !
    ! compare with                beta = 0.046d0
    !
    ! in e.g. http://dft.uci.edu/pubs/PBEsol.html
    !
    real(RK) :: beta ! set to either beta_pbe or beta_pbesol at runtime

    select case(op)
    case ("pbe")
       beta = beta_pbe
    case ("pbesol")
       beta = beta_pbesol
    case default
       ABORT(op)
    end select

    select case (ispin)
    case (1)
       do i = 1, vla
          !
          ! Total density and total gamma as independent variables:
          !
          nn = var(1, n(i, 1))
          g2 = var(2, G(i, 1))

          f = c_gga(half * nn, half * nn, g2, beta)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, 1) = dFdn(i, 1) + fst(1, f)

          dFdg(i, 1) = dFdg(i, 1) + fst(2, f)

          if (present(dFdndn)) then
             dFdndn(i, 1) = dFdndn(i, 1) + snd(1, 1, f)
          endif

          if (present(dFdndg)) then
             dFdndg(i, 1) = dFdndg(i, 1) + snd(1, 2, f)
          endif

          if (present(dFdgdg)) then
             dFdgdg(i, 1) = dFdgdg(i, 1) + snd(2, 2, f)
          endif
       enddo

    case (2)
       do i = 1, vla
          !
          ! Two densities and total gamma as independent variables:
          !
          na = var(1, n(i, A))
          nb = var(2, n(i, B))
          g2 = var(3, G(i, AA) + G(i, BB) + 2 * G(i, AB))

          f = c_gga(na, nb, g2, beta)

          Fc(i) = Fc(i) + val(f)

          dFdn(i, A) = dFdn(i, A) + fst(1, f)
          dFdn(i, B) = dFdn(i, B) + fst(2, f)

          dFdg(i, AA) = dFdg(i, AA) + fst(3, f)
          dFdg(i, BB) = dFdg(i, BB) + fst(3, f)
          dFdg(i, AB) = dFdg(i, AB) + fst(3, f) * 2

          if (present(dFdndn)) then
             dFdndn(i, AA) = dFdndn(i, AA) + snd(1, 1, f)
             dFdndn(i, BB) = dFdndn(i, BB) + snd(2, 2, f)
             dFdndn(i, AB) = dFdndn(i, AB) + snd(1, 2, f)
          endif

          if (present(dFdndg)) then
             dFdndg(i, AAA) = dFdndg(i, AAA) + snd(1, 3, f)
             dFdndg(i, ABB) = dFdndg(i, ABB) + snd(1, 3, f)
             dFdndg(i, AAB) = dFdndg(i, AAB) + snd(1, 3, f) * 2
             dFdndg(i, BAA) = dFdndg(i, BAA) + snd(2, 3, f)
             dFdndg(i, BBB) = dFdndg(i, BBB) + snd(2, 3, f)
             dFdndg(i, BAB) = dFdndg(i, BAB) + snd(2, 3, f) * 2
          endif

          if (present(dFdgdg)) then
             dFdgdg(i, AAAA) = dFdgdg(i, AAAA) + snd(3, 3, f)
             dFdgdg(i, AABB) = dFdgdg(i, AABB) + snd(3, 3, f)
             dFdgdg(i, AAAB) = dFdgdg(i, AAAB) + snd(3, 3, f) * 2
             dFdgdg(i, BBBB) = dFdgdg(i, BBBB) + snd(3, 3, f)
             dFdgdg(i, BBAB) = dFdgdg(i, BBAB) + snd(3, 3, f) * 2
             dFdgdg(i, ABAB) = dFdgdg(i, ABAB) + snd(3, 3, f) * 4
          endif
       enddo
    case default
       ABORT("ispin")
    end select
  end subroutine pbe_ggcc

end module pbe_ggcxc_module

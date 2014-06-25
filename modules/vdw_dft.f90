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
module vdw_dft
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
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

  use type_module, only: rk => r8_kind ! type specification parameters
  use iso_c_binding
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

#ifdef WITH_GUILE
  public :: vdw_dft_phi!(SCM d1, SCM d2) -> SCM double
#endif

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  real(rk), parameter :: pi = 4 * atan(1.0d0)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains


  function phi(d1, d2)
    !
    ! Van  der Waals  Density  Functional for  General Geometries,  M.
    ! Dion, H.   Rydberg, E.  Schröder,  D.  C.  Langreth, and  B.  I.
    ! Lundqvist,    Phys.    Rev.     Lett.    92,    246401   (2004),
    ! http://dx.doi.org/10.1103/PhysRevLett.92.246401
    !
    implicit none
    real(rk), intent(in), value :: d1, d2
    real(rk) :: phi
    ! *** end of interface ***

    real(rk), parameter :: cut = 35.0 ! cutoff "radius" for a**2 + b**2
    integer, parameter :: nrad = 100
    integer, parameter :: nang = 50
    integer :: i, k
    real(rk) :: u, v, a, b
    real(rk) :: du, dv

    !
    ! FIXME: braindead integration  over the quater of the  plane a in
    ! [0, oo),  b in [0, oo)  using the symmetry of  the integrad with
    ! respect to exchange of a and b:
    !
    du = cut / nrad
    dv = pi / (4 * nang) ! (pi / 4) / nang

    phi = 0.0
    do i = 0, nrad - 1
       u = (2 * i + 1) * du / 2 ! (i + 1/2) du
       do k = 0, nang - 1
          v = (2 * k + 1) * dv / 2 ! (k + 1/2) dv
          a = u * cos(v)
          b = u * sin(v)
          phi = phi + f(a, b, d1, d2) * u
       enddo
    enddo
    phi = phi * du * dv * 2

    contains

      pure function T(w, x, y, z)
        implicit none
        real(rk), intent(in) :: w, x, y, z
        real(rk) :: T
        ! *** end of interface ***

        T = (1 / (w + x) + 1 / (y + z)) * (1 / ((w + y) * (x + z)) + 1 / ((w + z) * (y + x))) / 2
      end function T

      pure function w(a, b)
        !
        !            2 2
        ! w(a, b) = a b  W(a, b)
        !
        ! with W(a, b) introduced in the literature.
        !
        implicit none
        real(rk), intent(in) :: a, b
        real(rk) :: W
        ! *** end of interface ***

        w =  2 * ((3 - a**2) * cos(b) * sinc(a) &
                + (3 - b**2) * cos(a) * sinc(b) &
                + (a**2 + b**2 - 3) * sinc(a) * sinc(b) &
                - 3 * cos(a) * cos(b))
      end function w

      pure function nu(y, d)
        !
        !              2                                     2
        ! nu(y, d) = (y / 2) / (1 - exp(- (ETA / 2) * (y / d) ))
        !
        implicit none
        real(rk), intent(in) :: y, d
        real(rk) :: nu
        ! *** end of interface ***

        real(rk), parameter :: ETA = 8 * pi / 9
        real(rk) :: x

        x = (ETA / 2) * y**2 / d**2

        if ((ETA / 2) * y**2 > 0.01 * d**2) then
           !
           ! To avoid 0.0 * inf for d = 0.0 in the alternative branch:
           !
           nu = (y**2 / 2) / (1 - exp(-x))
        else
           !
           ! Use Taylor series for small y/d ratios:
           !
           nu = (d**2 / ETA) * nx(x)
        endif
      end function nu

      pure function nx(x)
        !
        ! nx(x) = x / (1 - exp(-x))
        !
        ! taylor(nx(x), x, 0, 8) =
        !
        !              2    4      6        8
        !         x   x    x      x        x
        ! =   1 + - + -- - --- + ----- - ------- + . . .
        !         2   12   720   30240   1209600
        !
        implicit none
        real(rk), intent(in) :: x
        real(rk) :: nx
        ! *** end of interface ***

        if (abs(x) > 0.01) then
           nx = x / (1 - exp(-x))
        else
           nx = 1 + x/2 + x**2/12 - x**4/720 + x**6/30240 - x**8/1209600
        endif
      end function nx

      pure function f(a, b, d1, d2)
        implicit none
        real(rk), intent(in) :: a, b, d1, d2
        real(rk) :: f
        ! *** end of interface ***

        !
        !            2 2
        ! w(a, b) = a b  W(a, b) here:
        !
        f = (2 / pi**2) * w(a, b) * T(nu(a, d1), nu(b, d1), nu(a, d2), nu(b, d2))
      end function f

      pure function sinc(x)
        !
        ! sinc(x) = sin(x)/x
        !
        implicit none
        real(rk), intent(in) :: x
        real(rk) :: sinc
        ! *** end of interface ***

        if (abs(x) > 0.01) then
           sinc = sin(x) / x
        else
           ! (%i1) taylor(sin(x)/x, x, 0, 8);
           ! >                           2    4      6       8
           ! >                          x    x      x       x
           ! > (%o1)/T/             1 - -- + --- - ---- + ------ + . . .
           ! >                          6    120   5040   362880
           !
           ! Below x  < 0.01, terms greater than  x**8 contribute less
           ! than 1e-16  and so  are unimportant for  double precision
           ! arithmetics.
           !
           sinc = 1 - x**2/6 + x**4/120 - x**6/5040 + x**8/362880
        endif
      end function sinc

    end function phi

#ifdef WITH_GUILE
  function vdw_dft_phi(d1, d2) bind(c)
    !
    ! Export phi(d1, d2) to Scheme.
    !
    use scm, only: scm_t, scm_from, scm_to_double
    implicit none
    type(scm_t), intent(in), value :: d1, d2
    type(scm_t) :: vdw_dft_phi
    ! *** end of interface ***

    vdw_dft_phi = scm_from (phi (scm_to_double (d1), scm_to_double (d2)))
  end function vdw_dft_phi
#endif

  !--------------- End of module ----------------------------------
end module vdw_dft

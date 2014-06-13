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
module becke_step_func
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

! define SLOW_REFERENCE_CODE
! define MANUALDIFF
# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------


  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: mu0
  public :: mu1
  public :: mu2
  public :: muab

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----



  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function mu0(m,a) result(m0)
    !
    ! computes the value of the
    ! odd polynomial of order 27
    !
    !       mu3(m) = p(p(p(m))), p(m) = (3m-m^3)/2
    !
    !       mu3(1) = p(1) = 1  =>  mu3(-1) = p(-1) = -1
    !
    ! If paramter ``a'' is not zero then
    !
    !       mu3(m) = p(p(p(n(m))
    !
    ! with           m - a
    !       n(m) = ----------,   n(1) = 1,   n(-1) = -1
    !               1 - a m
    !
    use ad0x1 ! for testing purposes only ...
    implicit none
    real(r8_kind), intent(in) :: m
    real(r8_kind), intent(in) :: a
    real(r8_kind)             :: m0
    ! *** end of interface ***

    integer(i4_kind) :: k
    type(ad)         :: p ! AD-variable with overloaded arithmetics

    p = var(m) ! independent variable

    ! account for different atomic sizes:
    ! increases the cell (mu<0) for big (art>0) atoms:
    ! by replacing m by n = (m-a)/(1-a*m) :
    p = (p - a) / (1.0D0 - a * p)

    ! compute p(p(p(m)), with p(m) = (3m-m^3)/2
    do k=1,3
      p = p * (3.0D0 - p*p) / 2.0D0
    enddo

    m0 = val(p) ! value
  end function mu0

  function mu1(m,a) result(m1)
    !
    ! computes the values and first derivative of the
    ! odd polynomial of order 27
    !
    !       mu3(m) = p(p(p(m))), p(m) = (3m-m^3)/2
    !
    !       mu3(1) = p(1) = 1  =>  mu3(-1) = p(-1) = -1
    !
    ! If paramter ``a'' is not zero then
    !
    !       mu3(m) = p(p(p(n(m))
    !
    ! with           m - a
    !       n(m) = ----------,   n(1) = 1,   n(-1) = -1
    !               1 - a m
    !
    use ad1x1
    implicit none
    real(r8_kind), intent(in) :: m
    real(r8_kind), intent(in) :: a
    real(r8_kind)             :: m1(0:1)
    ! *** end of interface ***

    integer(i4_kind) :: k
    type(ad)         :: p ! AD-variable with overloaded arithmetics

    p = var(m) ! independent variable

    ! account for different atomic sizes:
    ! increases the cell (mu<0) for big (art>0) atoms:
    ! by replacing m by n = (m-a)/(1-a*m) :
    p = (p - a) / (1.0D0 - a * p)

    ! compute p(p(p(m)), with p(m) = (3m-m^3)/2
    do k=1,3
      p = p * (3.0D0 - p*p) / 2.0D0
    enddo

    m1(0) = val(p) ! value
    m1(1) = fst(p) ! first derivative
  end function mu1

  function mu2(m,a) result(m2)
    !
    ! computes the values and derivatives of the
    ! odd polynomial of order 27
    !
    !       mu3(m) = p(p(p(m))), p(m) = (3m-m^3)/2
    !
    !       mu3(1) = p(1) = 1  =>  mu3(-1) = p(-1) = -1
    !
    ! If paramter ``a'' is not zero then
    !
    !       mu3(m) = p(p(p(n(m))
    !
    ! with           m - a
    !       n(m) = ----------,   n(1) = 1,   n(-1) = -1
    !               1 - a m
    !
    use ad2x1
    implicit none
    real(r8_kind), intent(in) :: m
    real(r8_kind), intent(in) :: a
    real(r8_kind)             :: m2(0:2)
    ! *** end of interface ***

    integer(i4_kind) :: k
    type(ad)         :: p ! AD-variable with overloaded arithmetics

    p = var(m) ! independent variable

    ! account for different atomic sizes:
    ! increases the cell (mu<0) for big (art>0) atoms:
    ! by replacing m by n = (m-a)/(1-a*m) :
    p = (p - a) / (1.0D0 - a * p)

    ! compute p(p(p(m)), with p(m) = (3m-m^3)/2
    do k=1,3
      p = p * (3.0D0 - p*p) / 2.0D0
    enddo

    m2(0) = val(p) ! value
    m2(1) = fst(p) ! first derivative
    m2(2) = snd(p) ! second derivative
    ! or shorter: m2 = taylor(p)
  end function mu2

#ifdef AUTODIFF

  subroutine muab(a,b,ma,mb,maa,mab,mbb)
    !
    ! derivatives of
    !
    !         |a| - |b|     |a|     |b|
    !  mu  = ---------- == ----- - -----
    !    ab   | a - b |    |a-b|   |b-a|
    !
    !                        |       |
    !                       (1)     (2)
    !
    use ad2x6 ! autogenerated with -DORD=2 -DNVAR=6 -DAUTODIFF=ad2x6, see Makefile
    implicit none
    real(r8_kind), intent(in)  :: a(3),b(3)
    real(r8_kind), intent(out) :: ma(3)    ! d/da
    real(r8_kind), intent(out) :: mb(3)    ! d/db
    real(r8_kind), intent(out) :: maa(3,3) ! d^2/da^2
    real(r8_kind), intent(out) :: mab(3,3) ! d^2/dadb
    real(r8_kind), intent(out) :: mbb(3,3) ! d^2/db^2
    ! *** end of interface ***

    type(ad) :: av(3),bv(3),cv(3),ra,rb,rc,m
    integer  :: i,j

    do i=1,3
      ! components of vector "a" are independent vars 1,2,3
      av(i) = var(i  ,a(i))
      bv(i) = var(i+3,b(i))
      ! components of vector "b" are independent vars 4,5,6
    enddo

    ra = sqrnorm(av) ! square norm, defined in autodiff-module
    rb = sqrnorm(bv) ! for vectors of any length, but may be computed
                     ! similarly to "rc" below
    cv = av - bv
    rc = sqrt( cv(1)**2 + cv(2)**2 + cv(3)**2 )
                     ! arithmetic operators and (some) functions are "overloaded"
    m  = (ra - rb) / rc

    ! now extract the results:
    ! mu = val(m) ! if needed
    do j=1,3
      ! first derivatives are computed each time as well:
      ma(j) = fst(j  ,m)          ! d/da
      mb(j) = fst(j+3,m)          ! d/db
      do i=1,3
        ! extract second derivatives:
        maa(i,j) = snd(i  ,j  ,m) ! d^2 /  da^2
        mab(i,j) = snd(i  ,j+3,m) ! d^2 / da db
        mbb(i,j) = snd(i+3,j+3,m) ! d^2 /  db^2
      enddo
    enddo
  end subroutine muab

#else

  subroutine muab(a,b,ma,mb,maa,mab,mbb)
    !
    ! derivatives of
    !
    !         |a| - |b|     |a|     |b|
    !  mu  = ---------- == ----- - -----
    !    ab   | a - b |    |a-b|   |b-a|
    !
    !                        |       |
    !                       (1)     (2)
    !
    implicit none
    real(r8_kind), intent(in)  :: a(3),b(3)
    real(r8_kind), intent(out) :: ma(3)    ! d/da
    real(r8_kind), intent(out) :: mb(3)    ! d/db
    real(r8_kind), intent(out) :: maa(3,3) ! d^2/da^2
    real(r8_kind), intent(out) :: mab(3,3) ! d^2/dadb
    real(r8_kind), intent(out) :: mbb(3,3) ! d^2/db^2
    ! *** end of interface ***

    real(r8_kind) :: m1a(3),m1b(3),m1aa(3,3),m1ab(3,3),m1bb(3,3)
    real(r8_kind) :: m2b(3),m2a(3),m2bb(3,3),m2ba(3,3),m2aa(3,3)

    call mua(a,b,m1a,m1b,m1aa,m1ab,m1bb)
    call mua(b,a,m2b,m2a,m2bb,m2ba,m2aa)

    ! mu = m1 - m2
    ma  = m1a - m2a
    mb  = m1b - m2b
    maa = m1aa - m2aa
    mbb = m1bb - m2bb
    mab = m1ab - transpose(m2ba)
  end subroutine muab

  subroutine mua(a,b,ma,mb,maa,mab,mbb)
    !
    ! derivatives of
    !
    !         |a|
    !  mu  = -----
    !    a   |a-b|
    !
    implicit none
    real(r8_kind), intent(in)  :: a(3),b(3)
    real(r8_kind), intent(out) :: ma(3)    ! d/da
    real(r8_kind), intent(out) :: mb(3)    ! d/db
    real(r8_kind), intent(out) :: maa(3,3) ! d^2/da^2
    real(r8_kind), intent(out) :: mab(3,3) ! d^2/dadb
    real(r8_kind), intent(out) :: mbb(3,3) ! d^2/db^2
    ! *** end of interface ***

    integer(i4_kind) :: i,j
    real(r8_kind)    :: ra,rb,rc,c(3)
    real(r8_kind)    :: u,v,ua(3),va(3),uij,vij

    c = a - b
    ra = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
    rb = sqrt(b(1)**2 + b(2)**2 + b(3)**2)
    rc = sqrt(c(1)**2 + c(2)**2 + c(3)**2)

    ! mu == u*v, u = |a|, v = 1/|a-b|
    u  = ra
    v  = 1/rc

    ua = a/ra         ! du/da
    va = - c / rc**3  ! dv/da

    ! first derivatives:
    ma = ua * v + u * va
    mb =        - u * va

    ! second derivatives:
    do j=1,3
    do i=1,3
      !                    2
      !   2        delta  a  - a a
      !  d u           ij        i j
      ! ------- == -----------------
      ! da da              3
      !   i  j            a
      !
      !                             2
      !   2        3 c c  - delta  c
      !  d v          i j        ij       
      ! ------- == ------------------
      ! da da              5
      !   i  j            c
      !

      uij =   - a(i) * a(j)
      vij = 3 * c(i) * c(j)

      if( i == j )then
        uij = uij + ra**2
        vij = vij - rc**2
      endif

      uij = uij / ra**3
      vij = vij / rc**5

      ! second derivatives of
      !
      !      m(a,b) = u(a) * v(a-b):
      !
      maa(i,j) = ua(i) * va(j) &
               + ua(j) * va(i) &
               + uij * v       &
               + u * vij
      ! note that: du/db == 0 and dv/db = -dv/da:
      mbb(i,j) =   u * vij
      mab(i,j) = - u * vij     &
                 - ua(i) * va(j)
    enddo
    enddo
  end subroutine mua
#endif

#ifdef MANUALDIFF
  subroutine mu3(m,a)
    !
    ! computes the values and derivatives of the
    ! odd polynomial of order 27
    !
    !       mu3(m) = p(p(p(m))), p(m) = (3m-m^3)/2
    !
    !       mu3(1) = p(1) = 1  =>  mu3(-1) = p(-1) = -1
    !
    ! If paramter ``a'' is not zero then
    !
    !       mu3(m) = p(p(p(n(m))
    !
    ! with           m - a
    !       n(m) = ----------,   n(1) = 1,   n(-1) = -1
    !               1 - a m
    !
    implicit none
    real(r8_kind), intent(inout) :: m(0:) ! m(0:0), m(0:1), or m(0:2)
    real(r8_kind), intent(in)    :: a     ! cell size ratio a = (Ri-Rj)/(Ri+Rj)
    ! *** end of interface ***

    integer(i4_kind) :: k
    real(r8_kind)    :: p(0:2)

    select case(size(m))
    case(1) ! value only

      ! account for different atomic sizes:
      if( a /= 0 )then
        ! increases the cell (mu<0) for big (art>0) atoms:
        ! by replacing m by n = (m-a)/(1-a*m) :
        m(0) = (m(0) - a) / (1 - a * m(0))
      endif

      ! compute p(p(p(m)), with p(m) = (3m-m^3)/2
      do k=1,3
        m(0) = m(0) * (3 - m(0)**2) / 2
      enddo

    case(2) ! value and derivative

      ! account for different atomic sizes:
      if( a /= 0 )then
        ! increases the cell (mu<0) for big (art>0) atoms:
        ! by replacing m by n = p(m) = (m-a)/(1-a*m) :
        p(0) = (    m(0) - a) / (1 - a * m(0))
        p(1) = (1 + a * p(0)) / (1 - a * m(0)) ! p(0) is not a typo!

        ! now account for the dependencies of mu ...
        m(1) = p(1) * m(1)
        m(0) = p(0)
      endif

      ! compute p(p(p(m)), with p(m) = (3m-m^3)/2
      do k=1,3
        ! derivatives of p(x) wrt independent x at x==m:
        p(0) = m(0) * (3 -   m(0)**2) / 2
        p(1) =        (3 - 3*m(0)**2) / 2

        ! now account for the dependencies of mu ...
        m(1) = p(1) * m(1)
        m(0) = p(0)
      enddo

    case(3) ! value, first- and second derivatives

      ! account for different atomic sizes:
      if( a /= 0 )then
        ! increases the cell (mu<0) for big (art>0) atoms:
        ! by replacing m by n = p(m) = (m-a)/(1-a*m) :
        p(0) = (    m(0) - a) / (1 - a * m(0))
        p(1) = (1 + a * p(0)) / (1 - a * m(0)) ! p(0) is not a typo!
        p(2) = (2 * a * p(1)) / (1 - a * m(0)) ! p(1) is not a typo!

        ! now account for the dependencies of mu ...
        ! Second derivative first! Because mu(1) on the RHS
        ! will be modified by the second statement:
        m(2) = p(2) * m(1) * m(1) &
             + p(1) * m(2)
        m(1) = p(1) * m(1)
        m(0) = p(0)
      endif

      ! compute p(p(p(m)), with p(m) = (3m-m^3)/2
      do k=1,3
        ! derivatives of p(x) wrt independent x at x==m:
        p(0) = m(0) * (3 -   m(0)**2) / 2
        p(1) =        (3 - 3*m(0)**2) / 2
        p(2) =           - 3*m(0)

        ! now account for the dependencies of mu ...
        ! Second derivative first! Because mu(1) on the RHS
        ! will be modified by the second statement:
        m(2) = p(2) * m(1) * m(1) &
             + p(1) * m(2)
        m(1) = p(1) * m(1)
        m(0) = p(0)
      enddo

    case default
      ABORT('not implemented')
    end select
  end subroutine mu3
#endif

  !--------------- End of module ----------------------------------
end module becke_step_func


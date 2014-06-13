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
!
! MODULE AUTODIFF, MASTER FILE
!
! Copyright (c) 2008-2013 Alexei Matveev
!
# include "def.h"

#ifndef NVAR
#  define NVAR 1
#endif

! by default compile second derivatives support:
#ifndef ORD
#  define ORD 2
#endif

#if ORD == 0
#  define autodiff autodiff0 ! no derivatives, for reference purposes only
#  define SND !
#  define FST !
#endif

#if ORD == 1
#  define autodiff autodiff1 ! first derivatives only
#  define SND !
#  define FST 
#endif

#if ORD == 2
#  define autodiff autodiff2 ! supports second derivative
#  define SND
#  define FST
#endif

! supply the name with -DAUTODIFF=modulename for other builds:
#ifdef AUTODIFF
#  undef autodiff
#  define autodiff AUTODIFF
#endif

! define EQUALITY, may be confusing, dont abuse
# define ORDERING
# define FUNCTIONS
! define ELEMENTAL elemental (but replace ABORTs by something appropriate first!)
# define ELEMENTAL 

module autodiff

  use type_module, only: &
    IK=>i4_kind, RK=>r8_kind ! type specification parameters
  implicit none
  ! save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  type ad
    private
    real(RK) :: d0
#if NVAR == 1
FST real(RK) :: d1
SND real(RK) :: d2
#else
#define TR(a,b) a+b*(b-1)/2
FST real(RK) :: d1(NVAR)
SND real(RK) :: d2(TR(NVAR,NVAR)) ! only (a,b) with a <= b is used
#endif
  end type 

! use only these macros in the code below:
#if NVAR == 1

#define D0(x) x%d0
#define D1(x) x%d1
#define D2(x) x%d2

#define DA(x) D1(x)
#define DB(x) D1(x)

#define DECL1 
#define DECL2 

#define CODE0(expr)     D0(z)=expr
#define CODE1(expr) FST D1(z)=expr
#define CODE2(expr) SND D2(z)=expr

#else

#define D0(x) x%d0
#define D1(x) x%d1(a)
#define D2(x) x%d2(TR(a,b))

#define DA(x) x%d1(a)
#define DB(x) x%d1(b)

#define DECL1 integer :: a
#define DECL2 integer :: a SND,b

#define CODE0(expr) D0(z)=expr
#define CODE1(expr) FST do a=1,NVAR;D1(z)=expr;enddo
#define CODE2(expr) SND do b=1,NVAR;do a=1,b;D2(z)=expr;enddo;enddo

#endif

  !------------ Declaration of constants and variables ------------


  !------------ Interface statements ------------------------------

  interface operator(+)
      module procedure add_cv
      module procedure add_vc
      module procedure add_vv
  end interface

  interface operator(-)
      module procedure sub_cv
      module procedure sub_vc
      module procedure sub_vv
      module procedure usub_v  ! Unary minus!
  end interface

  interface operator(*)
      module procedure mult_cv
      module procedure mult_vc
      module procedure mult_vv
  end interface

  interface recip
      module procedure recip_v
      module procedure recip_c
  end interface

  interface operator(/)
      module procedure div_cv
      module procedure div_vc
      module procedure div_vv
  end interface

#ifdef EQUALITY
  !
  ! These compare the values, derivatives may differ.
  ! You have been warned!
  !
  interface operator(==)
      module procedure eq_cv
      module procedure eq_vc
      module procedure eq_vv
  end interface

  interface operator(/=)
      module procedure ne_cv
      module procedure ne_vc
      module procedure ne_vv
  end interface
#endif

#ifdef ORDERING
  interface operator(>)
      module procedure gt_cv
      module procedure gt_vc
      module procedure gt_vv
  end interface

  interface operator(>=)
      module procedure ge_cv
      module procedure ge_vc
      module procedure ge_vv
  end interface

  interface operator(<)
      module procedure lt_cv
      module procedure lt_vc
      module procedure lt_vv
  end interface

  interface operator(<=)
      module procedure le_cv
      module procedure le_vc
      module procedure le_vv
  end interface
#endif

#ifdef FUNCTIONS
  interface operator(**)
      module procedure expon_cv
      module procedure expon_vc ! x ** n, real n
      module procedure expon_vv
      module procedure expon_vn ! x ** n, integer n
  end interface

  interface sin
      module procedure sin_v
  end interface

  interface cos
      module procedure cos_v
  end interface

  interface tan
      module procedure tan_v
  end interface

  interface atan
      module procedure atan_v
  end interface

  interface exp
      module procedure exp_v
  end interface

  interface log
      module procedure log_v
  end interface

  interface sqrt
      module procedure sqrt_v
  end interface

  interface sqrnorm
      module procedure sqrnorm_v
      module procedure sqrnorm_c
  end interface
#endif


  !------------ public functions and subroutines ------------------

    public :: ad

    public :: var ! value    -> type(ad) -- "independent variable" (val,1,0,...)
    public :: fix ! value    -> type(ad) -- "fixed variable"       (val,0,0,...)
    public :: val ! type(ad) -> value
FST public :: fst ! type(ad) -> first derivative
SND public :: snd ! type(ad) -> second derivative
    public :: nan ! ()       -> (NaN, NaN, ...)

  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)

#ifdef EQUALITY
  !
  ! These compare the values, derivatives may differ.
  ! You have been warned!
  !
  public :: operator(==)
  public :: operator(/=)
#endif

#ifdef ORDERING
  public :: operator(>)
  public :: operator(>=)
  public :: operator(<)
  public :: operator(<=)
#endif

#ifdef FUNCTIONS
  public :: operator(**)
  public :: sin
  public :: cos
  public :: tan
  public :: atan
  public :: exp
  public :: log
  public :: sqrt
  public :: arcsinh
  public :: sqrnorm
#endif

#if NVAR == 1
  public :: taylor ! ! type(ad) -> real array (0:ORD)
#endif

  public :: show!(ad)

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----



  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  ! derivvar
  !     Properly initiate an independent variable
  ! Arguments:
  !     val        Value for the variable
  ! Result:
  !     Independent variable with value "val" and dervative 1
  !

  elemental function fix(val) result(z)
    real(RK), intent(in) :: val
    type(ad)             :: z

    DECL2

    CODE0( val )
    CODE1(  0  )
    CODE2(  0  )
  end function fix

  elemental function nan(x) result(z)
    !
    ! It is not possible to STOP/ABORT in pure funcitons.
    ! Use NaNs instead of failing silently.
    !
    implicit none
    type(ad), intent(in) :: x ! is ignored
    type(ad)             :: z
    ! *** end of interface ***

    real(RK) :: negative, invalid

    DECL2

    negative = -1.0_rk
    invalid = sqrt(negative)

    CODE0( invalid )
    CODE1( invalid )
    CODE2( invalid )
  end function nan

#if NVAR == 1

  elemental function var(val) result(z)
    real(RK), intent(in) :: val
    type(ad)             :: z

    CODE0( val )
    CODE1(  1  )
    CODE2(  0  )
  end function var

#else

  elemental function var(a,val) result(z)
    integer,  intent(in) :: a
    real(RK), intent(in) :: val
    type(ad)             :: z

!   ASSERT(a>0)
!   ASSERT(a<=NVAR)
    z%d0    = val
FST z%d1    = 0
FST z%d1(a) = 1
SND z%d2    = 0
  end function var

#endif

  elemental function val(x) result(v)
    type(ad), intent(in) :: x
    real(RK)             :: v

    v = D0(x)
  end function val

#if NVAR == 1

  pure function taylor(x) result(series)
    type(ad), intent(in) :: x
    real(RK)             :: series(0:ORD)

    series(0) = val(x)
FST series(1) = fst(x)
SND series(2) = snd(x)
  end function taylor

FST elemental function fst(x) result(d)
FST   type(ad), intent(in) :: x
FST   real(RK)             :: d
FST
FST   d = D1(x)
FST end function fst

#else

FST elemental function fst(a,x) result(d)
FST   integer,  intent(in) :: a
FST   type(ad), intent(in) :: x
FST   real(RK)             :: d
FST
!ST   ASSERT(a>0)
!ST   ASSERT(a<=NVAR)
FST   d = D1(x)
FST end function fst

#endif

#if NVAR == 1

SND elemental function snd(x) result(d)
SND   type(ad), intent(in) :: x
SND   real(RK)             :: d
SND
SND   d = D2(x)
SND end function snd

#else

SND ELEMENTAL function snd(a,b,x) result(d)
SND   integer,  intent(in) :: a,b
SND   type(ad), intent(in) :: x
SND   real(RK)             :: d
SND
SND
SND   ASSERT(a>0)
SND   ASSERT(a<=NVAR)
SND   ASSERT(b>0)
SND   ASSERT(b<=NVAR)
SND   ! only the upper triangle is used:
SND   d = x%d2(TR(min(a,b),max(a,b)))
SND end function snd

#endif

  ! +, -, *, /, **
  !     Implementation of the arithmetic operations
  ! Arguments:
  !     x          First operand (constant or variable)
  !     y          Second operand (constant or variable)
  ! Result:
  !     The sum, etc. and its first (total) derivative
  ! Note:
  !     We need three versions to allow combinations
  !     with constants.
  !

  elemental function add_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( D0(x) + D0(y) )
    CODE1( D1(x) + D1(y) )
    CODE2( D2(x) + D2(y) )
  end function add_vv

  elemental function add_cv( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( x + D0(y) )
    CODE1(     D1(y) )
    CODE2(     D2(y) )
  end function add_cv

  elemental function add_vc( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( D0(x) + y )
    CODE1( D1(x)     )
    CODE2( D2(x)     )
  end function add_vc

  elemental function usub_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0( - D0(x) )
    CODE1( - D1(x) )
    CODE2( - D2(x) )
  end function usub_v

  elemental function sub_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( D0(x) - D0(y) )
    CODE1( D1(x) - D1(y) )
    CODE2( D2(x) - D2(y) )
  end function sub_vv

  elemental function sub_cv( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( x - D0(y) )
    CODE1(   - D1(y) )
    CODE2(   - D2(y) )
  end function sub_cv

  elemental function sub_vc( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( D0(x) - y )
    CODE1( D1(x)     )
    CODE2( D2(x)     )
  end function sub_vc

  elemental function mult_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0(D0(x)*D0(y))
    CODE1(D1(x)*D0(y)+D0(x)*D1(y))
    CODE2(D2(x)*D0(y)+D0(x)*D2(y)+DA(x)*DB(y)+DB(x)*DA(y))
  end function mult_vv

  elemental function mult_cv( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( x * D0(y) )
    CODE1( x * D1(y) )
    CODE2( x * D2(y) )
  end function mult_cv

  elemental function mult_vc( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( D0(x) * y )
    CODE1( D1(x) * y )
    CODE2( D2(x) * y )
  end function mult_vc

  elemental function recip_v( y ) result( z )
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0(       1 / D0(y) )
    CODE1( - D1(y) / D0(y)**2 )
    CODE2( - D2(y) / D0(y)**2 + 2 * DA(y) * DB(y) / D0(y)**3 )
  end function recip_v

  elemental function recip_c( y ) result( z )
    real(RK), intent(in) :: y
    real(RK)             :: z

    z = 1 / y
  end function recip_c

  elemental function div_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    z = x * recip(y)
  end function div_vv

  elemental function div_cv( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    z = x * recip(y)
  end function div_cv

  elemental function div_vc( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    type(ad)             :: z

    z = x * recip(y)
  end function div_vc

#ifdef FUNCTIONS
  ! sin, cos, tan, log, exp, sqrt
  !     Return the sine etc. of x
  ! Arguments:
  !     x          Variable in question
  ! Result:
  !     The function value and its derivatives at x
  ! Note:
  !     Variables (not constants) need to be initialised as
  !                x = var(x)
  !
  elemental function sin_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0( sin(D0(x)) )
    CODE1( cos(D0(x)) * D1(x) )
    CODE2( cos(D0(x)) * D2(x) - sin(D0(x)) * DA(x) * DB(x) )
  end function sin_v

  elemental function cos_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0(   cos(D0(x)) )
    CODE1( - sin(D0(x)) * D1(x) )
    CODE2( - sin(D0(x)) * D2(x) - cos(D0(x)) * DA(x) * DB(x) )
  end function cos_v

  ELEMENTAL function tan_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0( tan(D0(x)) )
    CODE1( D1(x) / (1 + D0(z)**2 ) )
    ABORT('fix')
    CODE2( 0 )
  end function tan_v

  elemental function atan_v( x ) result( z )
    !
    ! maxima> diff(atan(x), x);
    !
    !               1
    ! (%o1)       ------
    !              2
    !             x  + 1
    ! maxima> diff(atan(x), x, 2);
    !
    !                2 x
    ! (%o2)     - ---------
    !               2     2
    !             (x  + 1)
    !
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0(atan(D0(x)))
    CODE1(D1(x)/(1+D0(x)**2))
    CODE2((D2(x)-2*D0(x)*DA(x)*DB(x)/(1+D0(x)**2))/(1+D0(x)**2))
  end function atan_v

  elemental function exp_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0( exp(D0(x)) )
    CODE1( D0(z) *   D1(x) )
    CODE2( D0(z) * ( D2(x) + DA(x) * DB(x) ) )
  end function exp_v

  elemental function log_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0( log(D0(x)) )
    CODE1( D1(x) / D0(x) )
    CODE2( (D2(x) - DA(x) * DB(x) / D0(x)) / D0(x) )
  end function log_v

  elemental function sqrt_v( x ) result( z )
    type(ad), intent(in) :: x
    type(ad)             :: z

    DECL2

    CODE0( sqrt(D0(x)) )
    CODE1(  D1(x) / (2 * D0(z)) )
    CODE2(( D2(x) - DA(x) * DB(x) / (2 * D0(x)) ) / (2 * D0(z)) )
!   z = x**0.5D0
  end function sqrt_v

  ELEMENTAL function expon_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( D0(x)  ** D0(y) )
    CODE1( (D1(y) * log(D0(x)) + D0(y)/D0(x) * D1(x) ) * D0(z) )
    ABORT('fix')
    CODE2( 0 )
  end function expon_vv

  ELEMENTAL function expon_cv( x, y ) result( z )
    real(RK), intent(in)        :: x
    type(ad), intent(in) :: y
    type(ad)             :: z

    DECL2

    CODE0( x**D0(y) )
    CODE1( D1(y) * log(x) * D0(z) )
    ABORT('fix')
    CODE2( 0 )
  end function expon_cv

  elemental function expon_vc( x, n ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: n
    type(ad)             :: z

    DECL2

    CODE0(  D0(x)**n )
    CODE1(n*D0(x)**(n-1)*D1(x))
    CODE2(n*D0(x)**(n-1)*D2(x)+n*(n-1)*D0(x)**(n-2)*DA(x)*DB(x))
  end function expon_vc

  elemental function expon_vn( x, n ) result( z )
    type(ad), intent(in) :: x
    integer, intent(in)  :: n
    type(ad)             :: z

    DECL2

    CODE0(  D0(x)**n)
    CODE1(n*D0(x)**(n-1)*D1(x))
    CODE2(n*D0(x)**(n-1)*D2(x)+n*(n-1)*D0(x)**(n-2)*DA(x)*DB(x))
  end function expon_vn

  elemental function arcsinh(x) result(y)
    !
    !                                2
    ! arcsinh(x) = log(x + sqrt(1 + x )) =
    !
    !                   3      5      7
    !                  x    3 x    5 x
    !            = x - -- + ---- - ---- + . . .
    !                  6     40    112
    implicit none
    type(ad), intent(in) :: x
    type(ad)             :: y
    ! *** end of interface ***

    type(ad) :: x2, x4

    x2 = x * x
    if( abs(val(x)) >= 1.0e-4_rk )then
       y = log( x + sqrt( 1.0_rk + x2 ) )
    else
       x4 = x2 * x2

       ! taylor series with relative error of O(x^6):
       y = x                                 !  1.000000000000000E-004
       y = y - ( 1.0_rk /    6.0_rk) * x*x2  ! -1.666666666666667E-013
       y = y + ( 3.0_rk /   40.0_rk) * x*x4  !  7.500000000000001E-022
!      y = y - ( 5.0_rk /  112.0_rk) * x**7  ! -4.464285714285714E-030
!      y = y + (35.0_rk / 1152.0_rk) * x**9  !  3.038194444444445E-038
    endif
  end function arcsinh

  ! Square Norm of a Vector :: [type(ad)] -> type(ad)
  function sqrnorm_v( x ) result( z )
    type(ad), intent(in) :: x(:)
    type(ad)             :: z

    integer :: i

    z = fix(0D0)
    do i=1,size(x)
      z = z + x(i)**2
    enddo
    z = sqrt(z)
  end function sqrnorm_v

  ! Square Norm of a Vector :: [Double] -> Double
  function sqrnorm_c( x ) result( z )
    real(RK), intent(in) :: x(:)
    real(RK)             :: z

    integer :: i

    z = 0D0
    do i=1,size(x)
      z = z + x(i)**2
    enddo
    z = sqrt(z)
  end function sqrnorm_c
#endif

#ifdef EQUALITY
  !
  ! x == y, x /= y may be confusing, dont abuse
  !
  elemental function eq_cv( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    logical              :: z

    z = val(x) .eq. y
  end function eq_cv

  elemental function eq_vc( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = x .eq. val(y)
  end function eq_vc

  elemental function eq_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = val(x) .eq. val(y)
  end function eq_vv

  elemental function ne_cv( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    logical              :: z

    z = val(x) .ne. y
  end function ne_cv

  elemental function ne_vc( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = x .ne. val(y)
  end function ne_vc

  elemental function ne_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = val(x) .ne. val(y)
  end function ne_vv
#endif

#ifdef ORDERING
  !
  ! x > y
  !
  elemental function gt_cv( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    logical              :: z

    z = val(x) .gt. y
  end function gt_cv

  elemental function gt_vc( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = x .gt. val(y)
  end function gt_vc

  elemental function gt_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = val(x) .gt. val(y)
  end function gt_vv

  !
  ! x >= y
  !
  elemental function ge_cv( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    logical              :: z

    z = val(x) .ge. y
  end function ge_cv

  elemental function ge_vc( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = x .ge. val(y)
  end function ge_vc

  elemental function ge_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = val(x) .ge. val(y)
  end function ge_vv

  !
  ! x < y
  !
  elemental function lt_cv( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    logical              :: z

    z = val(x) .lt. y
  end function lt_cv

  elemental function lt_vc( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = x .lt. val(y)
  end function lt_vc

  elemental function lt_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = val(x) .lt. val(y)
  end function lt_vv

  !
  ! x <= y
  !
  elemental function le_cv( x, y ) result( z )
    type(ad), intent(in) :: x
    real(RK), intent(in) :: y
    logical              :: z

    z = val(x) .le. y
  end function le_cv

  elemental function le_vc( x, y ) result( z )
    real(RK), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = x .le. val(y)
  end function le_vc

  elemental function le_vv( x, y ) result( z )
    type(ad), intent(in) :: x
    type(ad), intent(in) :: y
    logical              :: z

    z = val(x) .le. val(y)
  end function le_vv
#endif

  subroutine show(v)
    implicit none
    type(ad), intent(in) :: v
    ! *** end of interface ***

#if NVAR != 1
    integer :: i
SND integer :: j
#endif

    print *, "val:", val(v)
#if NVAR == 1
FST print *, "fst:", fst(v)
SND print *, "snd:", snd(v)
#else
FST print *, "fst:", (fst(i, v), i=1,NVAR)
SND print *, "snd:", ((snd(i, j, v), j=1,i), i=1,NVAR)
#endif
  end subroutine show

  !--------------- End of module ----------------------------------
end module autodiff

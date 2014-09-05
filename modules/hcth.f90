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
module hcth
  !
  ! Copyright (c) Ajanta Deka
  ! Copyright (c) Alexei Matveev
  !
# include "def.h"
  use type_module, RK=>r8_kind, IK=>i4_kind
  implicit none

  save
  private

  public &
       hcth_x, &
       hcth_c

  integer(IK), parameter, public :: &
       IXCHNG      =  1, &
       ICORR       =  2, &
       IXCHNG_ORIG = -1, &
       ICORR_ORIG  = -2

  integer(IK), parameter :: &
       IVAL  = 0, &
       IRA   = 1, &
       IRB   = 2, &
       IZA   = 3, &
       IZB   = 4, &
       IZAB  = 5

  integer(IK), parameter :: max_pow_u = 4
  real (RK),   parameter :: &
       pi        = 3.14159265358979324_RK ,&
       DEPS      = 1.0e-15_RK ,&
       zero      = 0.0_RK ,&
       one       = 1.0_RK ,&
       two       = 2.0_RK ,&
       three     = 3.0_RK ,&
       four      = 4.0_RK ,&
       five      = 5.0_RK ,&
       six       = 6.0_RK ,&
       seven     = 7.0_RK ,&
       eght      = 8.0_RK ,&
       elvn      = 11.0_RK ,&
       sixtn     = 16.0_RK , &
       thrd      = one/three   ,&
       tothrd    = two/three   ,&
       frthrd    = four/three  ,&
       fvthrd    = five/three ,&
       svnthrd   = seven/three ,&
       eghthrd   = eght/three ,&
       elvnthrd  = elvn/three ,&
       sixtnthrd = sixtn/three

  real (RK),   parameter :: &
       HCTH_PARAMETERS(3*5*4) = (/ &
       ! HCTH_version==1
       ! please refer to these coeff's as THCH1/iterate-e750-g500-v1-m4-n4
       0.109320E+01_rk, & ! X   0
       0.222601E+00_rk, & ! Caa 0
       0.729974E+00_rk, & ! Cab 0
      -0.744056E+00_rk, & ! X   1
      -0.338622E-01_rk, & ! Caa 1
       0.335287E+01_rk, & ! Cab 1
       0.559920E+01_rk, & ! X   2
      -0.125170E-01_rk, & ! Caa 3
      -0.115430E+02_rk, & ! Cab 3
      -0.678549E+01_rk, & ! ...
      -0.802496E+00_rk, & !
       0.808564E+01_rk, & !
       0.449357E+01_rk, & !
       0.155396E+01_rk, & !
      -0.447857E+01_rk, & !
       ! HCTH_version==2
       0.109163E+01_rk, & !
       0.489508E+00_rk, & !
       0.514730E+00_rk, & !
      -0.747215E+00_rk, & !
      -0.260699E+00_rk, & !
       0.692982E+01_rk, & !
       0.507833E+01_rk, & !
       0.432917E+00_rk, & !
      -0.247073E+02_rk, & !
      -0.410746E+01_rk, & !
      -0.199247E+01_rk, & !
       0.231098E+02_rk, & !
       0.117173E+01_rk, & !
       0.248531E+01_rk, & !
      -0.113234E+02_rk, & !
       ! HCTH_version==3
       0.109025E+01_rk, & !
       0.562576E+00_rk, & !
       0.542352E+00_rk, & !
      -0.799194E+00_rk, & !
       0.171436E-01_rk, & !
       0.701464E+01_rk, & !
       0.557212E+01_rk, & !
      -0.130636E+01_rk, & !
      -0.283822E+02_rk, & !
      -0.586760E+01_rk, & !
       0.105747E+01_rk, & !
       0.350329E+02_rk, & !
       0.304544E+01_rk, & !
       0.885429E+00_rk, & !
      -0.204284E+02_rk, & !
       ! HCTH_version==4
       0.108184E+01_rk, & !
       0.118777E+01_rk, & !
       0.589076E+00_rk, & !
      -0.518339E+00_rk, & !
      -0.240292E+01_rk, & !
       0.442374E+01_rk, & !
       0.342562E+01_rk, & !
       0.561741E+01_rk, & !
      -0.192218E+02_rk, & !
      -0.262901E+01_rk, & !
      -0.917923E+01_rk, & !
       0.425721E+02_rk, & !
       0.228855E+01_rk, & !
       0.624798E+01_rk, & !
      -0.420052E+02_rk  & !
      /)

  ! USE THESE INSTEAD:
  real(RK), parameter :: &
       CHCTH(0:4,3,4) = reshape( HCTH_PARAMETERS, (/ 5, 3, 4 /), ORDER=(/2,1,3/))
       ! 0:max_pow_u == 0:4
       !    five terms corresponding
       !    to powers of u = s^2/(1+s^2)
       ! 1:3
       !    for eXchange (IPARX), AA=BB (IPARCAA,IPARCBB), and AB (IPARCAB) correlation
       ! 1:4
       !    four versions of HCTH parameters

  ! to address the second dimension of CHCTH:
  integer(IK), parameter :: &
       IPARX   = 1, &
       IPARCAA = 2, &
       IPARCBB = IPARCAA, &
       IPARCAB = 3

contains

  subroutine hcth_x(rho,gamma,ispin,vl,HCTH_version,F_x,dF_xdrho,dF_xdg)
    real(RK), intent(in)     :: rho(:,:)            ! (1:vl,1:ispin)
    real(RK), intent(in)     :: gamma(:,:)          ! (1:vl,1:1+2*(ispin-1))
    integer(IK), intent(in)  :: ispin,vl,HCTH_version
    real(RK), intent(inout)  :: F_x(:)              ! (1:vl)
    real(RK), intent(inout)  :: dF_xdrho(:,:)       ! (1:vl,1:ispin)
    real(RK), intent(inout)  :: dF_xdg(:,:)         ! (1:vl,1:1+2*(ispin-1))
    ! *** end of interface ***

    call pg_hcth(IXCHNG,rho,gamma,ispin,vl,HCTH_version,F_x,dF_xdrho,dF_xdg)
  end subroutine hcth_x

  subroutine hcth_c(rho,gamma,ispin,vl,HCTH_version,F_c,dF_cdrho,dF_cdg)
    real(RK), intent(in)     :: rho(:,:)            ! (1:vl,1:ispin)
    real(RK), intent(in)     :: gamma(:,:)          ! (1:vl,1:1+2*(ispin-1))
    integer(IK), intent(in)  :: ispin,vl,HCTH_version
    real(RK), intent(inout)  :: F_c(:)              ! (1:vl)
    real(RK), intent(inout)  :: dF_cdrho(:,:)       ! (1:vl,1:ispin)
    real(RK), intent(inout)  :: dF_cdg(:,:)         ! (1:vl,1:1+2*(ispin-1))
    ! *** end of interface ***

    call pg_hcth(ICORR,rho,gamma,ispin,vl,HCTH_version,F_c,dF_cdrho,dF_cdg)
  end subroutine hcth_c

  subroutine pg_hcth(what,rho,gamma,ispin,vl,HCTH_version,Fxc,dFxcdrho,dFxcdg)
    integer(IK), intent(in)  :: what
    real(RK), intent(in)     :: rho(:,:)            ! (1:vl,1:ispin)
    real(RK), intent(in)     :: gamma(:,:)          ! (1:vl,1:1+2*(ispin-1))
    integer(IK), intent(in)  :: ispin,vl,HCTH_version
    real(RK), intent(inout)  :: Fxc(:)              ! (1:vl)
    real(RK), intent(inout)  :: dFxcdrho(:,:)       ! (1:vl,1:ispin)
    real(RK), intent(inout)  :: dFxcdg(:,:)         ! (1:vl,1:1+2*(ispin-1))
    ! *** end of interface ***

    integer(IK) :: i
    real(RK)    :: r(5)   ! ra,rb,za,zb,zab
    real(RK)    :: f(0:5) ! val,dra,drb,dza,dzb,dzab

    select case(ispin)
    case(1)
       do i=1,vl
          r(IRA)  = rho(i,1)/two
          r(IRB)  = r(IRA)
          r(IZA)  = gamma(i,1)/four
          r(IZB)  = r(IZA)
          r(IZAB) = gamma(i,1)/four ! DOESNT MATTER, IS NOT USED
          call hcth_calc(what,r,HCTH_version,f)
          Fxc(i)        = Fxc(i)        + f(IVAL)
          dFxcdrho(i,1) = dFxcdrho(i,1) + (f(IRA) + f(IRB))/two
          dFxcdg(i,1)   = dFxcdg(i,1)   + (f(IZA) + f(IZB) + f(IZAB))/four
       enddo
    case(2)
       do i=1,vl
          r(IRA)  = rho(i,1)
          r(IRB)  = rho(i,2)
          r(IZA)  = gamma(i,1)
          r(IZB)  = gamma(i,2)
          r(IZAB) = gamma(i,3) ! DOESNT MATTER, IS NOT USED
          call hcth_calc(what,r,HCTH_version,f)
          Fxc(i)        = Fxc(i)        + f(IVAL)
          dFxcdrho(i,1) = dFxcdrho(i,1) + f(IRA)
          dFxcdrho(i,2) = dFxcdrho(i,2) + f(IRB)
          dFxcdg(i,1)   = dFxcdg(i,1)   + f(IZA)
          dFxcdg(i,2)   = dFxcdg(i,2)   + f(IZB)
          dFxcdg(i,3)   = dFxcdg(i,3)   + f(IZAB) ! DOESNT MATTER, IS ZERO
       enddo
    case default
       print *,'pg_hcth: no such case: ispin=',ispin
       stop
    end select
  end subroutine pg_hcth

  subroutine hcth_calc(what,rho,HCTH_version,fxc)
    implicit none
    integer(IK), intent(in)    :: what
    real(RK),intent(in)        :: rho(5) ! ra,rb,za,zb,zab
    integer(IK),intent(in)     :: HCTH_version
    real(RK),intent(out)       :: fxc(0:5) ! val,dra,drb,dza,dzb,dzab
    ! *** end of interface ***

    select case(what)
    case (IXCHNG)
       call hcth_x_calc(rho,HCTH_version,fxc)
!!$    case (IXCHNG_ORIG)
!!$       call ORIG_hcth_x_calc(rho,HCTH_version,fxc)
    case (ICORR)
       call hcth_c_calc(rho,HCTH_version,fxc)
!!$    case (ICORR_ORIG)
!!$       call ORIG_hcth_c_calc(rho,HCTH_version,fxc)
    case default
       print *,'hcth_calc: no such case:',what
       stop 'hcth_calc'
    end select
  end subroutine hcth_calc

  subroutine hcth_x_calc(rho,HCTH_version,fx)
    implicit none
    real(RK),intent(in)        :: rho(5) ! ra,rb,za,zb,zab
    integer(IK),intent(in)     :: HCTH_version
    real(RK),intent(out)       :: fx(0:5) ! val,dra,drb,dza,dzb,dzab
    ! *** end of interface ***

    real(RK) :: cx(0:max_pow_u)

    cx(0:4) = CHCTH(0:4,IPARX,HCTH_version)

    fx(IVAL) = zero
    fx(IZAB) = zero
    call calc(rho(IRA),rho(IZA),fx(IVAL),fx(IRA),fx(IZA))
    call calc(rho(IRB),rho(IZB),fx(IVAL),fx(IRB),fx(IZB))

  contains

    subroutine calc(rho,z2,f,dfdr,dfdz2)
      implicit none
      real(RK), intent(in)    :: rho,z2
      real(RK), intent(inout) :: f
      real(RK), intent(out)   :: dfdr,dfdz2
      ! *** end of interfacce ***

      integer(IK), parameter :: &
           DR=1,    &
           DZ=2
      real(RK) :: ex,dex,s2,u,du
      real(RK) :: ds2(DR:DZ),duu(DR:DZ)
      real(RK) :: xalpha

      if (rho < DEPS)then
         dfdr  = zero
         dfdz2 = zero
         return
      endif

      ! LDAX:
      ! -(3/2) * [(3/4pi)^(1/3)] * rho^(4/3)
      xalpha  = - (three/two) * (three/(four*pi))**(one/three)
      ex  = xalpha * rho**(four/three)
      dex = xalpha * rho**( one/three) * (four/three)

      s2      = 0.004D0 * z2  / rho**(8.D0/3.D0)
      ds2(DR) = - (8.D0/3.D0) * s2/rho
      ds2(DZ) = 0.004D0 * one / rho**(8.D0/3.D0)

      call g(cx,s2,u,du)

      duu(DR) = du * ds2(DR)
      duu(DZ) = du * ds2(DZ)

      !  construction of f_xc itself
      !  ----------------------------

      f = f + ex * u

      dfdr  = ex * duu(DR) + dex * u

      dfdz2 = ex * duu(DZ)
    end subroutine calc

    subroutine g(c,s2,g0,g1)
      ! ... same as in hcth_c_calc() ...
      ! computes g0 = g_c = SUM[n=0..4] c(n) * u^n
      ! and its derivative g1 = d/ds2 g0
      real(RK), intent(in)  :: c(0:max_pow_u), s2
      real(RK), intent(out) :: g0,g1
      ! *** end of interface ***

      integer(IK) :: n
      real(RK)    :: u,du

      ! u = s2/(1+s2)
      u  = s2/(one + s2)
      du = one/(one + s2)**2
      ! polinomial c1*u + c2*u^2 + ...
      g0 = zero
      g1 = zero
      do n=max_pow_u,1,-1
         g0 =     c(n) + u * g0
         g1 = n * c(n) + u * g1
      enddo
      g0 = c(0) + u * g0
      g1 = g1 * du
    end subroutine g

  end subroutine hcth_x_calc

  subroutine hcth_c_calc( rho, HCTH_version, fc )
    implicit none
    real(RK),intent(in)    :: rho(5)  ! ra,rb,ga,gb,gab
    integer(IK),intent(in) :: HCTH_version
    real(RK),intent(out)   :: fc(0:5) ! val,dra,drb,dga,dgb,dgab
    ! *** end of interface ***

    real(RK), dimension(0:max_pow_u) :: caa,cab ! (5) cbb==caa

    integer(IK), parameter :: &
         IAA = 1, &
         IBB = 2, &
         IAB = 3
    real(RK) :: s2(IAA:IAB), u(IAA:IAB), du(IAA:IAB)
    integer(IK), parameter :: &
         IAADRA = 1, &
         IBBDRB = 2, &
         IABDRA = 3, &
         IABDRB = 4, &
         IAADZA = 5, &
         IBBDZB = 6, &
         IABDZA = 7, &
         IABDZB = 8
    real(RK) :: ds2(8), duu(8)
    real(RK) :: e(IAA:IAB)        ! (1:3)
    real(RK) :: de(IAADRA:IABDRB) ! (1:4)

    fc = zero

    if (rho(IRA)+rho(IRB) < DEPS) return

    !  construction of the LDA functional and the derivatives
    !  ------------------------------------------------------

      ! PWLDAC F >
      call pwldac( &
           rho(IRA), rho(IRB), &
           e(IAA), e(IBB), e(IAB), &
           de(IAADRA), de(IBBDRB), &
           de(IABDRA), de(IABDRB)  &
           )
      ! < PWLDAC F

    !  construction of the u terms and their derivatives
    !  -------------------------------------------------------

      s2(IAA) = rho(IZA) / max(rho(IRA),DEPS)**(8.D0/3.D0)
      s2(IBB) = rho(IZB) / max(rho(IRB),DEPS)**(8.D0/3.D0)
      s2(IAB) = 0.5D0*(s2(IAA) + s2(IBB)) ! avg

      ds2(IAADRA) = - (8.D0/3.D0) * s2(IAA)/max(rho(IRA),DEPS)
      ds2(IBBDRB) = - (8.D0/3.D0) * s2(IBB)/max(rho(IRB),DEPS)
      ds2(IAADZA) = one / max(rho(IRA),DEPS)**(8.D0/3.D0)
      ds2(IBBDZB) = one / max(rho(IRB),DEPS)**(8.D0/3.D0)
      ds2(IABDRA) = 0.5D0 * ds2(IAADRA)
      ds2(IABDRB) = 0.5D0 * ds2(IBBDRB)
      ds2(IABDZA) = 0.5D0 * ds2(IAADZA)
      ds2(IABDZB) = 0.5D0 * ds2(IBBDZB)

      caa(0:4) = CHCTH(0:4,IPARCAA,HCTH_version)
      cab(0:4) = CHCTH(0:4,IPARCAB,HCTH_version)

      call g(caa,0.2D0   * s2(IAA),u(IAA),du(IAA))
      call g(caa,0.2D0   * s2(IBB),u(IBB),du(IBB)) ! caa==cbb
      call g(cab,0.006D0 * s2(IAB),u(IAB),du(IAB))


      duu(IAADRA) = du(IAA) * (0.2D0   * ds2(IAADRA)) ! 1
      duu(IBBDRB) = du(IBB) * (0.2D0   * ds2(IBBDRB)) ! 2
      duu(IABDRA) = du(IAB) * (0.006D0 * ds2(IABDRA)) ! 3
      duu(IABDRB) = du(IAB) * (0.006D0 * ds2(IABDRB)) ! 4
      duu(IAADZA) = du(IAA) * (0.2D0   * ds2(IAADZA)) ! 5
      duu(IBBDZB) = du(IBB) * (0.2D0   * ds2(IBBDZB)) ! 6
      duu(IABDZA) = du(IAB) * (0.006D0 * ds2(IABDZA)) ! 7
      duu(IABDZB) = du(IAB) * (0.006D0 * ds2(IABDZB)) ! 8

      !  construction of f_xc itself
      !  ----------------------------

      fc(IVAL) = &
           e(IAA) * u(IAA) + &
           e(IBB) * u(IBB) + &
           e(IAB) * u(IAB)

      ! now the derivatives:
       fc(IRA) = &
            de(IAADRA) * u(IAA) + &
            e(IAA) * duu(IAADRA) + &
            de(IABDRA) * u(IAB) + &
            e(IAB) * duu(IABDRA)
       fc(IRB) = &
            de(IBBDRB) * u(IBB) + &
            e(IBB) * duu(IBBDRB) + &
            de(IABDRB) * u(IAB) + &
            e(IAB) * duu(IABDRB)
       fc(IZA) = &
            e(IAA) * duu(IAADZA) + &
            e(IAB) * duu(IABDZA)
       fc(IZB) = &
            e(IBB) * duu(IBBDZB) + &
            e(IAB) * duu(IABDZB)
       fc(IZAB) = zero

  contains

    subroutine g(c,s2,g0,g1)
      ! ... same as in hcth_x_calc() ...
      ! computes g0 = g_c = SUM[n=0..4] c(n) * u^n
      ! and its derivative g1 = d/ds2 g0
      real(RK), intent(in)  :: c(0:max_pow_u), s2
      real(RK), intent(out) :: g0,g1
      ! *** end of interface ***

      integer(IK) :: n
      real(RK)    :: u,du

      ! u = s2/(1+s2)
      u  = s2/(one + s2)
      du = one/(one + s2)**2
      ! polinomial c1*u + c2*u^2 + ...
      g0 = zero
      g1 = zero
      do n=max_pow_u,1,-1
         g0 =     c(n) + u * g0
         g1 = n * c(n) + u * g1
      enddo
      g0 = c(0) + u * g0
      g1 = g1 * du
    end subroutine g

  end subroutine hcth_c_calc

      subroutine pwldac( &
           rhoa, rhob, &
           e_caa, e_cbb, e_cab, &
           de_caa_by_drhoa, de_cbb_by_drhob, &
           de_cab_by_drhoa, de_cab_by_drhob &
           )
      use pw_ldac_module
      implicit none
      REAL(RK), intent(in)  :: rhoa, rhob
      REAL(RK), intent(out) :: e_caa, e_cbb, e_cab
      REAL(RK), intent(out) :: de_caa_by_drhoa, de_cbb_by_drhob
      REAL(RK), intent(out) :: de_cab_by_drhoa, de_cab_by_drhob
      ! *** end of iterface ***

      real(RK) :: rho(1,2),e(1),de(1,2)

      rho(1,1) = rhoa
      rho(1,2) = rhob

      ! e(rho=rhoa,z=1)
      e  = zero
      de = zero
      call pw_ldac(1,PWLDAC_ECP,rho(:,1:1),e,de(:,1:1))
      e_caa           = e(1)
      de_caa_by_drhoa = de(1,1)

      ! e(rho=rhob,z=1)
      e  = zero
      de = zero
      call pw_ldac(1,PWLDAC_ECP,rho(:,2:2),e,de(:,2:2))
      e_cbb           = e(1)
      de_cbb_by_drhob = de(1,2)

      ! true correlation:
      e  = zero
      de = zero
      call pw_ldac(1,PWLDAC_UNRESTR,rho,e,de)
      e_cab           = e(1) ! will be modified later
      de_cab_by_drhoa = de(1,1)
      de_cab_by_drhob = de(1,2)

      ! subtract *ECP(a) *ECP(b)
      e_cab           = e_cab - e_caa - e_cbb
      de_cab_by_drhoa = de_cab_by_drhoa - de_caa_by_drhoa
      de_cab_by_drhob = de_cab_by_drhob - de_cbb_by_drhob
      end subroutine pwldac

      subroutine print_parameters()
        implicit none
        ! *** end of interface ***

        integer :: i,j,k
        do k=1,4
           write(*,*) '=== HCTH VERSION ',k,' === '
           do j=1,3
              select case(j)
              case (1)
                 write(*,*) 'COEFFS FOR EXCHANGE:'
              case (2)
                 write(*,*) 'COEFFS FOR aa AND bb CORRELATION:'
              case (3)
                 write(*,*) 'COEFFS FOR ab CORRELATION:'
              end select
              do i=0,4
                 write(*,'(F12.6)') CHCTH(i,j,k)
              enddo
           enddo
        enddo
      end subroutine print_parameters

end module hcth

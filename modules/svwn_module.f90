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
!=============================================================================
! Public interface of module
!=============================================================================
module svwn_module
!------------------------------------------------------------------------------
!
!  Purpose: Contains routines for evaluating Slater-VWN-LDA exchange-correlation 
!           functionals:
!
!  Subroutine svwn_xalpha_calc()  -> LDA exchange    of Slater
!  Subroutine svwn_vwn_calc()     -> LDA correlation of Vosko, Wilk Nusair
!
!  Definitions and notation:
!
!  spin-symmetric case
!                         d/drho(r ) E_xc[rho(r')] = V_xc(r) delta(r-r')
!               d/drho(r) d/drho(r') E_xc[rho(r')] = Z_xc(r) delta(r-r')
!    d/drho(r) d/drho(r') d/drho(r") E_xc[rho(r')] = Q_xc(r) delta(r-r')delta(r-r")
!  
!  spin-polarized case
!    d/drho_s(r) d/drho_s'(r') E_xc[rho] = Z_xc,s,s'(r) delta(r-r')
!    d/drho_s(r) d/drho_s'(r') E_xc[rho] = Z_xc,s,s'(r) delta(r-r')
!
!    rho(r)       electronic density
!    E_xc[rho]    exchange-correlation energy functional
!    f_xc[rho]
!  eps_xc[rho]
!    V_xc[rho]    exchange-correlation potential
!    Z_xc[rho]    exchange-correlation kernel
!
!
!  Module called by: response_module
!
!  References: 
!   1. Based on "vwnc.f90" module of MS.
!   2. B.G. Johnson, P.M.W. Gill and. J.A. Pople,
!      "The performance of a family of density functional methods",
!      J.Chem.Phys. Vol.98 No.7, April 1993 pp. 5612.
!      (Erratum: JCP Vol.1 No.10, Nov. 15, pp. 9202)
! 
!  Author: HH, UB
!  Date: 3/99
!
!------------------------------------------------------------------------------
!== Interrupt of public interface of module ===================================
!------------------------------------------------------------------------------
! Modifications
!------------------------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!------------------------------------------------------------------------------

  use type_module
  implicit none
  private         ! by default, all names are private
  save

!== Interrupt end of public interface of module ===============================


!------------ public functions and subroutines --------------------------------
  public svwn_vwn_calc, svwn_xalpha_calc

!==============================================================================
  ! End of public interface of module
!==============================================================================

!------------ Declaration of constants and variables --------------------------
  real   (kind=r8_kind),parameter ::&
       & zero      = 0.0_r8_kind,&
       & one       = 1.0_r8_kind,&
       & two       = 2.0_r8_kind,&
       & three     = 3.0_r8_kind,&
       & four      = 4.0_r8_kind,&
       & five      = 5.0_r8_kind,&
       & six       = 6.0_r8_kind,&
       & eight     = 8.0_r8_kind,&
       & sixteen   =16.0_r8_kind,&
       & thirtytwo =32.0_r8_kind,&
       & thirtyfive=35.0_r8_kind,&
       & thirtysix =36.0_r8_kind,&
       & half      = 0.5_r8_kind,&
       & third     = one/three,& 
       & deps      = 1.0E-50_r8_kind, &
       & pi        = 3.141592653589793240_r8_kind

!------------------------------------------------------------------------------
!------------ Subroutines -----------------------------------------------------

contains

  !****************************************************************************
  ! Note: Qxc is only evaluated for spin-unpolarized case
  subroutine svwn_xalpha_calc(rho,Vxc,ispin,fxc,vec_length,eps,Zxc,Qxc,offset)
    implicit none
    !---- declaration of dummy arguments --------------------------------------
    real   (kind=r8_kind),intent(in   )          :: rho(:,:)
    real   (kind=r8_kind),intent(inout)          :: Vxc(:,:)
    integer(kind=i4_kind),intent(in   )          :: ispin
    real   (kind=r8_kind),intent(inout)          :: fxc(:)
    integer(kind=i4_kind),intent(in   )          :: vec_length
    real   (kind=r8_kind),intent(inout),optional :: eps(:)
    real   (kind=r8_kind),intent(inout),optional :: Zxc(:,:)
    real   (kind=r8_kind),intent(inout),optional :: Qxc(:)
    logical              ,intent(in   ),optional :: offset  ! default: .true.
    ! spin-symmetric case
    !   d/drho(r) d/drho(r') E_xc[rho] = Z_xc(r) delta(r-r')
    !   Zxc(:,1) = Z_xc(r_grid)
    ! spin-polarized case
    !   d/drho_s(r) d/drho_s'(r') E_xc[rho] = Z_xc,s,s'(r) delta(r-r')
    !   Zxc(:,1) = Z_xc,up,up(r_grid)
    !   Zxc(:,2) = Z_xc,dn,dn(r_grid)
    !   Zxc(:,3) = Z_xc,up,dn(r_grid) = Z_xc,dn,up(r_grid)
    !** End of interface ******************************************************
    real(kind=r8_kind) :: tmp(vec_length,ispin)
    real(kind=r8_kind) :: factor
    logical            :: use_offset
!>>>>begin debug
!!$    real(kind=r8_kind) :: rho2(vec_length)
!!$    integer(kind=i4_kind) :: i
!<<<<end debug
    !---- executable part -----------------------------------------------------
    
    if(present(offset)) then
       use_offset = .true.   ! use rho+deps instead of rho like in vwnc.f90
    else
       use_offset = .false.  ! directly work with rho
    end if

    if(ispin==1) then  ! non spinpolarized case
!..............................................................................
!      for alpha = 2/3 :
!      f_xc   = - 3/4 * (3/pi)^1/3 * rho^4/3
!      eps_xc = - 3/4 * (3/pi)^1/3 * rho^1/3
!      V_xc   = -       (3/pi)^1/3 * rho^1/3
!      Z_xc   = - 1/3 * (3/pi)^1/3 / rho^2/3
!      Q_xc   = + 2/9 * (3/pi)^1/3 / rho^5/3
!..............................................................................
       factor = - (3.0_r8_kind/pi)**third
       if (use_offset) then
!>>>>>begin original
          tmp(:,1) = factor * ( rho(1:vec_length,1) + deps )**third
          if (present(Zxc)) then
             Zxc(1:vec_length,1) = Zxc(1:vec_length,1) + &
                  third * tmp(:,1) / ( rho(1:vec_length,1) + deps )
          endif
          if (present(Qxc)) then
             Qxc(1:vec_length) = Qxc(1:vec_length) - &
                  (2.0_r8_kind/9.0_r8_kind) * tmp(:,1) / ( (rho(1:vec_length,1) + deps)**2 )
          endif
!<<<<end original
!
! try to simulate EXACTLY behaviour in vwnc.f90
!
!!$!>>>>begin debug
!!$          rho2(1:vec_length)=rho(1:vec_length,1)
!!$          do i=1,vec_length
!!$             if(rho2(i)<1.0E-10_r8_kind) rho2(i)=1.0E-10_r8_kind
!!$          end do
!!$          tmp(:,1) = factor * rho2(1:vec_length)**third
!!$          if (present(Zxc)) then
!!$             Zxc(1:vec_length,1) = Zxc(1:vec_length,1) + &
!!$                  third * tmp(:,1) / rho2(1:vec_length)
!!$          endif
!!$          if (present(Qxc)) then
!!$             Qxc(1:vec_length) = Qxc(1:vec_length) - &
!!$                  (2.0_r8_kind/9.0_r8_kind) * tmp(:,1) / ( rho2(1:vec_length)**2 )
!!$          endif
!!$!<<<<end debug
       else
          tmp(:,1) = factor * exp( log( rho(1:vec_length,1) ) * third )
          if (present(Zxc)) then
             Zxc(1:vec_length,1) = Zxc(1:vec_length,1) + &
                                      third * tmp(:,1) / rho(1:vec_length,1)
          endif
          if (present(Qxc)) then
             Qxc(1:vec_length) = Qxc(1:vec_length) - &
                  (2.0_r8_kind/9.0_r8_kind)  * tmp(:,1) / ( rho(1:vec_length,1)**2 )
          endif
       endif
       Vxc(1:vec_length,1) = Vxc(1:vec_length,1) + tmp(:,1)
       tmp(:,1) = 0.75_r8_kind * tmp(:,1)
       fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)*rho(1:vec_length,1)
       if (present(eps)) then
          eps(1:vec_length) = eps(1:vec_length) + tmp(:,1)
       endif

    else              ! spinpolarized case
!..............................................................................
!      for alpha = 2/3 :
!      f_xc    = - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^4/3
!      eps_xc  = - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^4/3 / Sum(s) rho_s
!        and not - 3/4 * (6/pi)^1/3 * Sum(s) rho_s^1/3 !
!      V_xc,s  = -       (6/pi)^1/3 * rho_s^1/3
!      Z_xc,st = - 1/3 * (6/pi)^1/3 / rho_s^2/3 delta_st
!..............................................................................
       factor = - (6.0_r8_kind/pi)**third
       if (use_offset) then
          tmp = factor * ( rho(1:vec_length,:) + deps )**third
          if (present(Zxc)) then
             Zxc(1:vec_length,1:2) = Zxc(1:vec_length,1:2) + &
                  third * tmp / ( rho(1:vec_length,:) + deps )
          endif
       else
          tmp = factor * exp( log( rho(1:vec_length,:) ) * third )
          if (present(Zxc)) then
             Zxc(1:vec_length,1:2) = Zxc(1:vec_length,1:2) + &
                                        third * tmp / rho(1:vec_length,:)
          endif
       endif
       Vxc(1:vec_length,:) = Vxc(1:vec_length,:) + tmp
       tmp = tmp * rho(1:vec_length,:)
       tmp(:,1) = 0.75_r8_kind * ( tmp(:,1) + tmp(:,2) )
       fxc(1:vec_length) = fxc(1:vec_length) + tmp(:,1)
!:BUG  ! eps(rho) is  n o t  factor * Sum(s) rho_s^(1/3) !
       if (present(eps)) then
          if (use_offset) then
          tmp(:,2) = rho(1:vec_length,1) + rho(1:vec_length,2) + deps
          else
             tmp(:,2) = rho(1:vec_length,1) + rho(1:vec_length,2)
          endif
          eps(1:vec_length) = eps(1:vec_length) + tmp(:,1) / tmp(:,2)
       endif
    endif

  end subroutine svwn_xalpha_calc
  !****************************************************************************




  !****************************************************************************
  ! Note: Qxc is only evaluated for spin-unpolarized case
  subroutine svwn_vwn_calc(rho,Vxc,ispin,fxc,vec_length,eps,Zxc,Qxc,offset)
    implicit none
    !---- declaration of dummy arguments --------------------------------------
    real   (kind=r8_kind),intent(in   )          :: rho(:,:)
    real   (kind=r8_kind),intent(inout)          :: Vxc(:,:)
    integer(kind=i4_kind),intent(in   )          :: ispin
    real   (kind=r8_kind),intent(inout)          :: fxc(:)
    integer(kind=i4_kind),intent(in   )          :: vec_length
    real   (kind=r8_kind),intent(inout),optional :: eps(:)
    real   (kind=r8_kind),intent(inout),optional :: Zxc(:,:)
    real   (kind=r8_kind),intent(inout),optional :: Qxc(:)
    logical              ,intent(in   ),optional :: offset  ! default: .true.
    ! spin-symmetric case
    !   d/drho(r) d/drho(r') E_xc[rho] = Z_xc(r) delta(r-r')
    !   Zxc(:,1) = Zxc(:,1) + 0.5*Z_xc,up,up(r_grid)
    !   Zxc(:,2) = Zxc(:,2) + 0
    !   Zxc(:,3) = Zxc(:,3) + 0.5*Z_xc,up,dn(r_grid) = 0.5*Z_xc,dn,up(r_grid)
    !   Qxc(:)   = Q_xc(r_grid)
    ! spin-polarized case
    !   d/drho_s(r) d/drho_s'(r') E_xc[rho] = Z_xc,s,s'(r) delta(r-r')
    !   Zxc(:,1) = Z_xc,up,up(r_grid)
    !   Zxc(:,2) = Z_xc,dn,dn(r_grid)
    !   Zxc(:,3) = Z_xc,up,dn(r_grid) = Z_xc,dn,up(r_grid)
    !** End of interface ******************************************************
    logical            :: use_offset
    integer(kind=i4_kind) :: vlen
    real   (kind=r8_kind),dimension(vec_length) :: &
         & x                 , & ! (3/4*pi*rho)^(1/6)
         & zeta              , & ! (rho_alpha - rho_beta) / rho
         & a_eps_c           , &
         & p_eps_c           , &
         & f_eps_c           , &
         & a_eps_c_prime     , &
         & p_eps_c_prime     , &
         & f_eps_c_prime     , &
         & a_eps_c_dblprime  , &
         & p_eps_c_dblprime  , &
         & f_eps_c_dblprime  , &
         & p_eps_c_trplprime , &
         & helpvec
    real   (kind=r8_kind)            :: &
         & a_Q,&
         & p_Q,&
         & f_Q
    real   (kind=r8_kind),parameter  :: &
         & a_A  = -0.0337740_r8_kind/two,&
         & p_A  =  0.0621841_r8_kind/two,&
         & f_A  =  0.0310907_r8_kind/two,&
         & a_B  =  1.13107_r8_kind      ,&
         & p_B  =  3.72744_r8_kind      ,&
         & f_B  =  7.06042_r8_kind      ,&
         & a_C  = 13.0045_r8_kind       ,&
         & p_C  = 12.9352_r8_kind       ,&
         & f_C  = 18.0578_r8_kind       ,&
         & a_x0 = -0.0047584_r8_kind    ,&
         & p_x0 = -0.10498_r8_kind      ,&
         & f_x0 = -0.32500_r8_kind

    !---- external subroutines ------------------------------------------------
!!$    external error_handler
    !---- executable part -----------------------------------------------------
    
    vlen = vec_length    ! actual vector length for current rho(:) 

    if(present(offset)) then
       use_offset = .true.   ! use rho+deps instead of rho like in vwnc.f90
    else
       use_offset = .false.  ! directly work with rho
    end if

    ! preset some constants needed for VWN functional
    a_Q = sqrt(four*a_C - a_B*a_B)  
    p_Q = sqrt(four*p_C - p_B*p_B)  
    f_Q = sqrt(four*f_C - f_B*f_B)  

    if(ispin==1) then  ! non spinpolarized case

       ! x =  (3/4*pi*rho)^(1/6)
       x = three/four
       if(use_offset) then
          x = sqrt( (three/(four*pi*(rho(1:vlen,1)+deps)))**third )
       else
          x = sqrt( (three/(four*pi*rho(1:vlen,1)))**third )
       end if

       ! auxiliary function p_eps_c(x) -> needed for fxc, Vxc, Zxc, eps
       p_eps_c = p_A*(  log(x*x/bigX(x,"p"))               &
                    & + (two*p_B/p_Q)*atan(p_Q/(2*x+p_B))  &
                    & - (p_B*p_x0/bigX0(p_x0,"p"))          &
                    &     *(  log((x-p_x0)**2/bigX(x,"p")) &
                    &        + ( two*(two*p_x0+p_B)/ p_Q)  &
                    &           * atan( p_Q/(two*x+p_B) )  &
                    &      )  &
                    &)

       ! derivative of p_eps_c(x) wrto x -> needed for Vxc, Zxc
       p_eps_c_prime = p_A*(  two/x - (two*x+p_B)/bigX(x,"p")     &
                          &  - four*p_B/( (two*x+p_B)**2+p_Q**2 ) &
                          &  - (p_B*p_x0/bigX0(p_x0,"p"))          &
                          &    * (  two/(x-p_x0)                  &
                          &       - (two*x+p_B)/bigX(x,"p")       &
                          &       - four*(two*p_x0+p_B)/(         &
                          &               (two*x+p_B)**2+p_Q**2 ) &
                          &      )                                &
                          &)

       ! 2nd derivative of p_eps_c(x) wrto x -> needed for Zxc
       p_eps_c_dblprime =  p_A*( -two/x**2 -two/bigX(x,"p")            &
                              &  +(two*x+p_B)**2/bigX(x,"p")**2        &
                              &  +sixteen*p_B*(two*x+p_B)              &
                              &    /( ( (two*x+p_B)**2+p_Q**2 )**2 )   &
                              &  -(p_B*p_x0/bigX0(p_x0,"p"))           &
                              &    * (-two/(x-p_x0)**2                 &
                              &       -two/bigX(x,"p")                 &
                              &       +(two*x+p_B)**2/bigX(x,"p")**2   &
                              &       +sixteen*(two*p_x0+p_B)          &
                              &          *(two*x+p_B)                  &
                              &            /((two*x+p_B)**2+p_Q**2)**2 &
                              &      )                                 &
                              &)

       ! 3rd derivative of p_eps_c(x) wrto x -> needed for Qxc
       p_eps_c_trplprime =  p_A*( -four/x**3                                     &
                              &   -(p_B*p_x0/bigX0(p_x0,"p"))*(four/(x-p_x0)**3) &
                              &   +(one-(p_B*p_x0/bigX0(p_x0,"p")))*             &
                              &    (                                             &
                              &       six*(two*x+p_B)/bigX(x,"p")**2             &
                              &      -two*( (two*x+p_B)/bigX(x,"p") )**3         &
                              &    )                                             &
                              &   +(                                             &
                              &       thirtytwo*p_B/( (two*x+p_B)**2 + p_Q**2 )**2   &
                              &      -eight*sixteen*p_B*(two*x+p_B)**2           &
                              &             /( (two*x+p_B)**2 + p_Q**2 )**3      &
                              &    )                                             &
                              &   *(                                             &
                              &     one-(p_B*p_x0/bigX0(p_x0,"p"))*(two*p_x0+p_B)&
                              &    ))                                             

       ! auxiliary function a_eps_c(x) -> needed for Zxc
       a_eps_c = a_A*(  log(x*x/bigX(x,"a"))               &
                    & + (two*a_B/a_Q)*atan(a_Q/(2*x+a_B))  &
                    & - (a_B*a_x0/bigX0(a_x0,"a"))          &
                    &     *(  log((x-a_x0)**2/bigX(x,"a")) &
                    &        + ( two*(two*a_x0+a_B)/ a_Q)  &
                    &           * atan( a_Q/(two*x+a_B) )  &
                    &      )  &
                    &)


       ! calculate final output values
       fxc(1:vlen)   = fxc(1:vlen)   + p_eps_c*rho(1:vlen,1)
       Vxc(1:vlen,1) = Vxc(1:vlen,1) + p_eps_c - (x/six)*p_eps_c_prime

       if(present(eps)) then
          eps(1:vlen) = eps(1:vlen) + p_eps_c
       end if

       if(use_offset) then

          if(present(Zxc)) then

             helpvec = -(five*x/thirtysix)*p_eps_c_prime &
                  & +(x*x/thirtysix)*p_eps_c_dblprime

             !    Zxc(:,1) = Zxc(:,1) + Z_xc,up,up(r_grid)
             Zxc(1:vlen,1) = Zxc(1:vlen,1) &
                  & + half*( helpvec+a_eps_c )/(rho(1:vlen,1)+deps)
             
!!$             !    Zxc(:,2) = Zxc(:,2) + 0
!!$             Zxc(1:vlen,2) = Zxc(1:vlen,2)
!!$             
!!$             !    Zxc(:,3) = Zxc(:,3) + Z_xc,up,dn(r_grid) = Z_xc,dn,up(r_grid)
!!$             Zxc(1:vlen,3) = Zxc(1:vlen,3) &
             Zxc(1:vlen,2) = Zxc(1:vlen,2) &
                  & + half*( helpvec-a_eps_c )/(rho(1:vlen,1)+deps)
             
          end if

          if(present(Qxc)) then

             helpvec = thirtyfive*p_eps_c_prime -three*x*p_eps_c_dblprime &
                  &    -x*x*p_eps_c_trplprime

             Qxc(1:vlen) = Qxc(1:vlen) &
                  & + x*helpvec/( six*thirtysix*( (rho(1:vlen,1)+deps)**2) )
          end if

       else

!!$          call error_handler("svwn_vwn_calc: offset not yet implemented ! ")

       end if

!..............................................................................

    else              ! spinpolarized case

       ! zeta = 
!!$       call error_handler("svwn_vwn_calc: spin=2 not yet implemented ! ")

!..............................................................................
!..............................................................................

    endif

  contains

    !==========================================================================
    function bigX(y,ctype)
      real(kind=r8_kind) :: y(vlen)
      character(len=*)   :: ctype
      real(kind=r8_kind) :: bigX(vlen)
      !------------------------------------------------------------------------
      
      select case (ctype)
      case ("p")
         bigX = y*y + p_B*y + p_C 
      case ("a")
         bigX = y*y + a_B*y + a_C 
      case ("f")
         bigX = y*y + f_B*y + f_C 
      case default
!!$         call error_handler("svwn_vwn_cal:X: unrecognized dummy argument &
!!$              &ctype="//trim(ctype))
         print*,"svwn_vwn_cal:X: unrecognized dummy argument ctype="//trim(ctype)
         stop
      end select
    
    end function bigX
    !==========================================================================
  
    !==========================================================================
    function bigX0(y,ctype)
      real(kind=r8_kind) :: y
      character(len=*)   :: ctype
      real(kind=r8_kind) :: bigX0(vlen)
      !------------------------------------------------------------------------
      
      select case (ctype)
      case ("p")
         bigX0 = y*y + p_B*y + p_C 
      case ("a")
         bigX0 = y*y + a_B*y + a_C 
      case ("f")
         bigX0 = y*y + f_B*y + f_C 
      case default
!!$         call error_handler("svwn_vwn_cal:X0: unrecognized dummy argument &
!!$              &ctype="//trim(ctype))
         print*, "svwn_vwn_cal:X0: unrecognized dummy argument ctype="//trim(ctype)
         stop
      end select
    
    end function bigX0
    !==========================================================================
  
 
  end subroutine svwn_vwn_calc
  !****************************************************************************

end module svwn_module





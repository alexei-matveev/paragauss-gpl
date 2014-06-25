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
module relxc
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
  !      Reference:E.Engel,S.Keller and R.M. Dreizler, Phys. Rev A, 53, 1367,1996.
  !      The used variables are as follows:
  !      phizero   : Rel LDA correction term
  !      dphzdb    : Derivative of phizero with respect to beta
  !      phi2      : Rel GGA correction term
  !      dph2db    : Derivative of phi2 with respect to beta
  !      n,nn      : Density
  !      Ga        : Squared Density gradient
  !      g         : GGA correction factor
  !      deNRLDAdn : Derivative of e(x)NRLDA with respect to n
  !      deRLDAdn  : Derivative of e(x)RLDA with respect to n
  !      deRGGAdn  : Derivative of e(x)RGGA with respect to n
  !      deRGGAdGa : Derivative of e(x)RGGA with respect to Squared Density Gradient
  !      eNRLDA, eRLDA, eRGGA : exchange functionals
  !
  !----------------------------------------------------------------

# include "def.h"
  use type_module, &
       & only: r8_kind, i4_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: rel_vwn_calc, rpw91c_calc

  !================================================================
  ! End of public interface of module
  !================================================================


  !------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  real(r8_kind), parameter    :: pi           = 3.14159265358979324_r8_kind ,&
       c0           = 137.03604_r8_kind ,&
       one          = 1.0_r8_kind ,&
       two          = 2.0_r8_kind ,&
       three        = 3.0_r8_kind ,&
       four         = 4.0_r8_kind ,&
       five         = 5.0_r8_kind ,&
       six          = 6.0_r8_kind ,&
       seven        = 7.0_r8_kind ,&
       eight        = 8.0_r8_kind ,&
       nine         = 9.0_r8_kind ,&
       sixtn        = 16.0_r8_kind ,&
       tweneight    = 28.0_r8_kind ,&
       thirfive     = 35.0_r8_kind ,&
       forty        = 40.0_r8_kind ,&
       sevtwo       = 72.0_r8_kind ,&
       hundt        = 100.0_r8_kind ,&
       oneonetwo    = 112.0_r8_kind ,&
       onetwoeight  = 128.0_r8_kind ,&
       onesevfiv    = 175.0_r8_kind ,&
       threonefiv   = 315.0_r8_kind ,&
       sixforty     = 640.0_r8_kind ,&
       sevhunfor    = 704.0_r8_kind ,&
       onesixonesev = 1617.0_r8_kind ,&
       deps         = 1.0e-50_r8_kind ,&
       x            = 1.0e-45_r8_kind
  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

 subroutine rel_vwn_calc(nn,ispin,vl,ec,dec)
    real(r8_kind), intent(in)      :: nn(:,:)  ! (1:vl,1:ispin)
    integer(i4_kind), intent(in)   :: ispin,vl
    real(r8_kind), intent(inout)   :: ec(:)  ! (1:vl)
    real(r8_kind), intent(inout)   :: dec(:,:) ! (1:vl,1:ispin)

    ! *** end of interface ***


    integer(i4_kind) :: i !,s
    real(r8_kind)    :: n,ecc,decc !,ec_sum
     if(ispin==1)then
           do i=1,vl

              n=nn(i,1)
              call ec_RLDA(n,ecc,decc)
              ec(i) = ec(i) + ecc
              dec(i,1) = dec(i,1) + decc

            enddo

      else

            do i=1,vl
               ! correlation depends on the total density only:
               n = nn(i,1) + nn(i,2)
               call ec_RLDA(n,ecc,decc)
               ec(i)    = ec(i)    + ecc
               dec(i,1) = dec(i,1) + decc
               dec(i,2) = dec(i,2) + decc
           enddo
     endif
 end subroutine rel_vwn_calc

  subroutine rpw91c_calc(rho,gamma,ispin,vl,fxc,dfdrho,dfdgrarho)
     ! Purpose : Calculates the relativistic correction to the correlation functional pw91c
     ! Note : No LDA and NRGGA contributions are included, calls to LDA(VWN)
     !        and GGA(pw91c) functionals must preceed a call to this functional.

    use spin_orbit_module, only: c => speed_of_light
    implicit none
    real(kind=r8_kind),dimension(:,:) :: rho  ! rho(vl,ispin)
    real(kind=r8_kind),dimension(:,:) :: gamma ! gamma(vl,(ispin-1)*2+1)
    real(kind=r8_kind),dimension(:),intent(inout) :: fxc
    real(kind=r8_kind),dimension(:,:),intent(inout) :: dfdrho,dfdgrarho
    integer(kind=i4_kind),intent(in) :: ispin,vl

   !** End of interface *****************************************

    real(kind=r8_kind),parameter :: a1 = 1.9407_r8_kind, &
                                    a2 = 0.14435_r8_kind, &
                                    b1 = 0.28412_r8_kind, &
                                    b2 = 0.004723_r8_kind

   real(kind=r8_kind),dimension(vl) :: b,num,denom,drho_num,drho_denom,phi_c,drho_phi_c,sums,drho1_phi_c,drho2_phi_c

  if(ispin==1)then

       sums = rho(1:vl,1)
       where(sums<deps) sums = deps
       b = ((three*pi**2*sums)**(one/three))/c
       num = one + a1*b**2 + a2*b**4
       denom = one + b1*b**2 + b2*b**4
       drho_num = two*a1*b + four*a2*b**3
       drho_denom = two*b1*b + four*b2*b**3
       phi_c = num/denom
       drho_phi_c = (denom*drho_num - num*drho_denom)*b/(three*sums*denom**2)
       fxc(1:vl) = fxc(1:vl) * phi_c
       dfdrho(1:vl,1) = dfdrho(1:vl,1)*phi_c + fxc*drho_phi_c
       dfdgrarho(1:vl,1) = dfdgrarho(1:vl,1) * phi_c
   else

      sums = rho(1:vl,1) + rho (1:vl,2)
      where(sums<deps)sums = deps
      b = ((three*pi**2*sums)**(one/three))/c
      num = one + a1*b**2 + a2*b**4
      denom = one + b1*b**2 + b2*b**4
      drho_num = two*a1*b + four*a2*b**3
      drho_denom = two*b1*b + four*b2*b**3
      phi_c = num/denom
      drho1_phi_c = (denom*drho_num - num*drho_denom)*b/(three*sums*denom**2)
      drho2_phi_c = (denom*drho_num - num*drho_denom)*b/(three*sums*denom**2)
      fxc(1:vl) = fxc(1:vl) * phi_c
      dfdrho(1:vl,1) = dfdrho(1:vl,1)*phi_c + fxc*drho1_phi_c
      dfdrho(1:vl,2) = dfdrho(1:vl,2)*phi_c + fxc*drho2_phi_c
      dfdgrarho(1:vl,1) = dfdgrarho(1:vl,1) * phi_c
      dfdgrarho(1:vl,2) = dfdgrarho(1:vl,2) * phi_c
      dfdgrarho(1:vl,3) = dfdgrarho(1:vl,3) * phi_c
   endif

   end subroutine rpw91c_calc


 subroutine ec_RLDA(n,ecRLDA,decRLDA)
    ! This subroutine calculates the relativistic correction to LDA correlation.
    ! So a call to this subroutine must be preceeded by a call to the VWN functional.
    use spin_orbit_module, only: c => speed_of_light
    implicit none
    real(kind=r8_kind), intent(in) :: n
    real(kind=r8_kind), intent(out) :: ecRLDA
    real(kind=r8_kind), intent(out) :: decRLDA
    !** End of interface *****************************************

    real(kind=r8_kind) :: b,num,denom,d_num,d_denom

    real(kind=r8_kind),parameter :: a0 = -0.16797_r8_kind, &
                                    b0 = 0.149530_r8_kind, &
                                    c0 = -0.007568_r8_kind, &
                                    d0 = 0.0504295_r8_kind, &
                                    c1 = -0.007217_r8_kind, &
                                    c2 = 0.0064572_r8_kind, &
                                    c3 = -0.005352_r8_kind, &
                                    d1 = 0.0995769_r8_kind, &
                                    d2 = 0.110576_r8_kind, &
                                    d3 = 0.0136740_r8_kind, &
                                    d4 = 0.0311479_r8_kind

!!$    ! for 0.014 cutoff
!!$    real(kind=r8_kind),parameter ::&
!!$         AA = -992.9266410_r8_kind, &
!!$         BB = 52205.81276_r8_kind
    ! for 0.25 cutoff:
    real(kind=r8_kind),parameter ::&
         AA = -0.2670592864_r8_kind, &
         BB =  0.5321072256_r8_kind

      if( n<1.0E-7_r8_kind )then
         ecRLDa = 0.0
         decRLDA = 0.0
         return
      endif

      b = ((three*(pi**2)*n)**(one/three))/c
#ifdef MAYER_ADEKA
      if(b<0.014_r8_kind)then ! Mayer and ADeka seem to use this cutoff
         ecRLDA  = 0.0_r8_kind
         decRLDA = 0.0_r8_kind
#else
      if(b<0.250_r8_kind)then
         ! approximate the function as
         ! y = AA x^3 + BB x^4, set AA and BB to match
         ! value and derivatives at 0.25 \pm 0
         ecRLDA  =          AA * b**3 +        BB * b**4
         decRLDA =( three * AA * b**2 + four * BB * b**3 ) * b/(three*n)
#endif
      else
           num = c3*b**3 + c2*b**2 + c1*b + c0
           denom = d4*b**4 + d3*b**3 + d2*b**2 + d1*b + d0
           d_num = three*c3*b**2 + two*c2*b + c1
           d_denom = four*d4*b**3 + three*d3*b**2 + two*d2*b + d1
           ecRLDA = a0*b + b0 + num/denom
           decRLDA = (a0 + (denom*d_num - num*d_denom)/denom**2)*b/(three*n)

     endif

 end subroutine ec_RLDA

 !--------------- End of module ----------------------------------
 end module relxc




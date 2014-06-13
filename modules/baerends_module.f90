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
module baerends_module
!---------------------------------------------------------------
!
!  Purpose: The module calculates the xc-functionals
!           from van Leeuwen and Baerends
!           It is only a function for the potential
!           (dfdrho), but not for the energy (fxc)
!           
!  Reference: Phys. Rev. A 49 (2421) 1994
!
!           The used variables are as follows:
!           rho:  density
!           gamma: square of the norm of the density gradients
!           dfdrho: derivative of with respect to rho (potential)
!           fxc:    f 
!           dfdgrarho: derivative of f with respect to grarho 
!                      not used
!
!  Module called by: xc_hamiltonian
!                    
!  Author: MS
!  Date: 2/96
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
  use type_module
  use type_module
  implicit none
  private
  save
!== Interrupt end of public interface of module =================

!------------ public functions and subroutines ------------------
  public :: baerends94_calc


!================================================================
! End of public interface of module
!================================================================

  real(kind=r8_kind), parameter :: ZERO  = 0.0_r8_kind ,&
       one    = 1.0_r8_kind ,&
       two    = 2.0_r8_kind ,&
       three  = 3.0_r8_kind ,&
       four   = 4.0_r8_kind ,&
       deps   = 1.0E-50_r8_kind , &
       third  = 1.0_r8_kind/3.0_r8_kind, &
       f4thrd = four/three, &
       beta   = 0.05

contains 



  subroutine baerends94_calc(rho,gamma,dfdrho,ispin,vec_length)
    ! purpose : Calculation of the functional of van Leeuwen and Baerends
    real(kind=r8_kind),dimension(:,:),intent(in) :: rho,gamma
    real(kind=r8_kind),dimension(:,:),intent(out) :: dfdrho
    integer(kind=i4_kind),intent(in) :: ispin,vec_length
    !** End of interface *****************************************       
    real(kind=r8_kind) :: two13
    real(kind=r8_kind),dimension(vec_length) :: sums,xgr,degr,&
         dro,xgrup,xgrdn,degrup,degrdn,droup,drodn
    real(kind=r8_kind),dimension(vec_length,ispin) :: gr

    two13=two**third

    if(ispin==1) then ! spin restricted

       sums=rho(1:vec_length,1)+deps
       gr(:,1)=sqrt(gamma(1:vec_length,1))
       xgr  = two13*gr(:,1)  /sums**f4thrd 
       degr = one + three*beta*xgr*arsinh(xgr)
       dro =  beta*rho(1:vec_length,1)**third * xgr**2/degr
       dro=dro/two13
       ! results
       ! dfdgrarho(1:vec_length,1) not set
       ! fxc(1:vec_length) not set 
       dfdrho(1:vec_length,1)=dfdrho(1:vec_length,1)-dro

    else    ! spinpolarized case
       sums = rho(1:vec_length,1) + rho(1:vec_length,1)+ deps
       gr(:,1:2)=sqrt(gamma(1:vec_length,1:2))
       xgrup  = gr(:,1)  /(rho(1:vec_length,1) + deps)**f4thrd 
       xgrdn  = gr(:,2)  /(rho(1:vec_length,2) + deps)**f4thrd 
       degrup = one + three*beta*xgrup*arsinh(xgrup)
       degrdn = one + three*beta*xgrdn*arsinh(xgrdn)
       droup =  beta*rho(1:vec_length,1)**third*xgrup*xgrup/degrup
       drodn =  beta*rho(1:vec_length,2)**third*xgrdn*xgrdn/degrdn

       ! results
       ! fxc(1:vec_length) not set       
       dfdrho(1:vec_length,1) = -droup+dfdrho(1:vec_length,1)
       dfdrho(1:vec_length,2) = -drodn+dfdrho(1:vec_length,2)
       ! dfdgrarho(1:vec_length,1) not set
       ! dfdgrarho(1:vec_length,2) not set
    end if
  end subroutine baerends94_calc

  !********************************************************************

  function arsinh(y)
    real(kind=r8_kind),intent(in),dimension(:) :: y
    real(kind=r8_kind),dimension(size(y)) :: arsinh
    arsinh=log(y+sqrt(y*y+1.0_r8_kind))
    return
  end function arsinh

end module baerends_module

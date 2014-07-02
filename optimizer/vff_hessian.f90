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
module  vff_hessian 
  !-------------------------------------------------------------------
  !
  !  Purpose : returns the appropriate force parameters for
  !           calculation the force constants (diagonals of
  !           the hesse matrix in valence coordinates).
  !           The data of this routine come from the EARLIER
  !           paper of Schlegel.
  !  Author: FN
  !  Date  : -- 
  !  References:  H.B.Schlegel, Theor.Chim.Acta 66,333-340 (1984)
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  !  Modifications
  !-------------------------------------------------------------------
  !  Purpose: Contains new data for estimate initial hessian (at B3LYP 
  !  level of theory)
  !  Modification:   
  !  Author: VVP
  !  Date: 6/05
  !  References: Wittbrodt, J.M., Schlegel, B.H. Estimating stretching 
  !  force constants for geometry optimization // Theochem.- 1997. - V. 
  !  398-399. - P. 55-61.
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  public          ! by default, all names are private
  !== Interrupt end of public interface of module ================= 
  real(kind=r8_kind)   :: stretch_par(6,6),avff
  real(kind=r8_kind),parameter      :: torsion_a = 0.0023_r8_kind,&
                                       torsion_b = 0.07_r8_kind
  real(kind=r8_kind),parameter      :: oop_par = 0.045_r8_kind
  real(kind=r8_kind),dimension(2,2) :: angle_par
  data stretch_par  /-0.2573_r8_kind, 0.3401_r8_kind, 0.6937_r8_kind, &
                      0.7126_r8_kind, 0.8335_r8_kind, 0.9491_r8_kind, &
                      0.3401_r8_kind, 0.9652_r8_kind, 1.2843_r8_kind, &
                      1.4725_r8_kind, 1.6549_r8_kind, 1.7190_r8_kind, &
                      0.6937_r8_kind, 1.2843_r8_kind, 1.6925_r8_kind, &
                      1.8238_r8_kind, 2.1164_r8_kind, 2.3185_r8_kind, &
                      0.7126_r8_kind, 1.4725_r8_kind, 1.8238_r8_kind, &
                      2.0203_r8_kind, 2.2137_r8_kind, 2.5206_r8_kind, &
                      0.8335_r8_kind, 1.6549_r8_kind, 2.1164_r8_kind, &
                      2.2137_r8_kind, 2.3718_r8_kind, 2.5110_r8_kind, &
                      0.9491_r8_kind, 1.7190_r8_kind, 2.3185_r8_kind, &
                      2.5206_r8_kind, 2.5110_r8_kind, 1.0000_r8_kind /
  data angle_par    /0.160_r8_kind,0.160_r8_kind,&
                     0.160_r8_kind,0.250_r8_kind /
end module vff_hessian
 

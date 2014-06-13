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
module constants

  use type_module, only: rk => r8_kind, ik => i4_kind
  public

  !
  ! Mathematical constants:
  !
  real(rk), parameter ::         &
       ZERO         =    0.0_rk, &
       HALF         =    0.5_rk, &
       ONE          =    1.0_rk, &
       THREE_HALF   =    1.5_rk, &
       TWO          =    2.0_rk, &
       FIVE_HALF    =    2.5_rk, &
       THREE        =    3.0_rk, &
       PI           =    3.1415926535897932384626433832794_rk, &
       FOUR         =    4.0_rk, &
       FIVE         =    5.0_rk, &
       SIX          =    6.0_rk, &
       SEVEN        =    7.0_rk, &
       TWO_PI       =    6.2831853071795864769352867665588_rk, &
       EIGHT        =    8.0_rk, &
       NINE         =    9.0_rk, &
       TEN          =   10.0_rk, &
       TWELVE       =   12.0_rk, &
       TWENTYSEVEN  =   27.0_rk, &
       HUNDRED      =  100.0_rk, &
       THOUSAND     = 1000.0_rk, &
       TWO_THOUSAND = 2000.0_rk

  ! DOUBLE FACTORIAL (2*L-1)!! at addresses   L=0:17
  ! DOUBLE FACTORIAL (2*L+1)!! at addresses 1+L=1:16
  real(rk), parameter, dimension(0:17) :: &
       DFAC = (/                   &
       1.00000000000000_rk,        & ! DFAC(0) == -1!!
       1.00000000000000_rk,        & ! DFAC(1) ==  1!!
       3.00000000000000_rk,        & ! DFAC(2) ==  3!!
       15.0000000000000_rk,        & ! DFAC(3) ==  5!!
       105.000000000000_rk,        & ! ...
       945.000000000000_rk,        & ! DFAC(L) == (2L - 1)!!
       10395.0000000000_rk,        &
       135135.000000000_rk,        &
       2027025.00000000_rk,        &
       34459425.0000000_rk,        &
       654729075.000000_rk,        &
       13749310575.0000_rk,        &
       316234143225.000_rk,        &
       7905853580625.00_rk,        &
       213458046676875._rk,        &
       6.190283353629375E+015_rk,  &
       1.918987839625106E+017_rk,  &
       6.332659870762850E+018_rk  /)

  ! FACTORIAL L! at addresses L=0:17
  real(rk), parameter, dimension(0:17) :: &
       FAC = (/              &
       1.00000000000000_rk,  & ! 0 -> 0!
       1.00000000000000_rk,  & ! 1 -> 1!
       2.00000000000000_rk,  & ! 2 -> 2!
       6.00000000000000_rk,  &
       24.0000000000000_rk,  &
       120.000000000000_rk,  &
       720.000000000000_rk,  &
       5040.00000000000_rk,  &
       40320.0000000000_rk,  &
       362880.000000000_rk,  &
       3628800.00000000_rk,  &
       39916800.0000000_rk,  &
       479001600.000000_rk,  &
       6227020800.00000_rk,  &
       87178291200.0000_rk,  &
       1307674368000.00_rk,  &
       20922789888000.0_rk,  &
       355687428096000._rk  /)


  !
  ! Fundamental  constants and  conversion  factors. Unless  specified
  ! otherwise the constants are in atomic units. Here is the length of
  ! an angstrom in au:
  !
  real(rk), parameter, private :: bohr = 0.52917706_rk ! in Angstrom
  real(rk), parameter :: angstrom = 1 / bohr ! ~= 1.9 au
  !
  ! FIXME: there are other literal constants in use all over the code,
  ! notably    0.529177249,    0.529177    and   a    few    remaining
  ! 0.52917706. Those were not yet replaced by the above parameter. Do
  ! this when/if you  feel comfortable and/or if the  precision is not
  ! critical in that part of the code.
  !

  !
  ! 1 kcal/mol in  Hartree, used at several places  dealing with force
  ! fields.   FIXME: the  file  molmech/common_data_module.f90 uses  a
  ! different value, namely 1 / 627.50956:
  !
  real(rk), parameter :: kcal = 1 / 627.49_rk


  real(rk), parameter :: &
       SPEED_OF_LIGHT = 137.03604_rk           , &
       c_speedoflight = 299792458.0_rk         , &            ! speed of light in SI-units
       R_gas          = 8.314472_rk            , &                    ! Ideal gas constant
       k_B            = 1.3806504E-23_rk       , &                    ! Boltzmann constant
       h_Planck       = 6.62606896E-34_rk      , &                       ! Planck constant
       hbar_Planck    = 1.05457163E-34_rk      , &                   ! Planck constant/2pi
       u_Atom_Mass    = 1.660538782E-027_rk    , &      ! atomic mass constant in SI-units
       a_Bohr         = 5.291772108E-11_rk     , &               ! bohr-radius in SI-units
       E_h2Jmol       = 2.6254995E+06_rk                       ! Conversion Eh to SI-units

end module constants

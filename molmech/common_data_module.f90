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
module common_data_module
  !------------ Modules used --------------------------------------
  use type_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  real(kind=r8_kind), public, parameter :: zero=0.0_r8_kind
  real(kind=r8_kind), public, parameter :: quarter=0.25_r8_kind
  real(kind=r8_kind), public, parameter :: half=0.5_r8_kind
  real(kind=r8_kind), public, parameter :: quarter3=0.75_r8_kind
  real(kind=r8_kind), public, parameter :: one=1.0_r8_kind
  real(kind=r8_kind), public, parameter :: mone=-1.0_r8_kind
  real(kind=r8_kind), public, parameter :: two=2.0_r8_kind
  real(kind=r8_kind), public, parameter :: three=3.0_r8_kind
  real(kind=r8_kind), public, parameter :: four=4.0_r8_kind
  real(kind=r8_kind), public, parameter :: five=5.0_r8_kind
  real(kind=r8_kind), public, parameter :: six=6.0_r8_kind
  real(kind=r8_kind), public, parameter :: eight=8.0_r8_kind
  real(kind=r8_kind), public, parameter :: nine=9.0_r8_kind
  real(kind=r8_kind), public, parameter :: ten=10.0_r8_kind
  real(kind=r8_kind), public, parameter :: twelve=12.0_r8_kind
  real(kind=r8_kind), public, parameter :: fifteen=15.0_r8_kind
  real(kind=r8_kind), public, parameter :: hundred=100.0_r8_kind
  real(kind=r8_kind), public, parameter :: thousand=1000.0_r8_kind
  real(kind=r8_kind), public, parameter :: infinity=10.0e20_r8_kind
  real(kind=r8_kind), public, parameter :: default_value=1.0e-5_r8_kind
  real(kind=r8_kind), public, parameter :: default_value1=1.0e-4_r8_kind

  real(kind=r8_kind), public, parameter :: pi=3.1415926535897932368_r8_kind

  real(kind=r8_kind), public, parameter :: pi_degree=180.0_r8_kind
  real(kind=r8_kind), public, parameter :: deg2rad=pi/pi_degree
  real(kind=r8_kind), public, parameter :: rad2deg=pi_degree/pi

  real(kind=r8_kind), public, parameter :: b2a=0.529177249_r8_kind
  real(kind=r8_kind), public, parameter :: a2b=1.889725988579_r8_kind
  real(kind=r8_kind), public, parameter :: h2eV=27.2113961_r8_kind
  real(kind=r8_kind), public, parameter :: h2kcm=627.50956_r8_kind
  real(kind=r8_kind), public, parameter :: h2kJm=2625.499999904_r8_kind
  real(kind=r8_kind), public, parameter :: c2J=4.184_r8_kind
  real(kind=r8_kind), public, parameter :: J2c=0.23900574_r8_kind
  real(kind=r8_kind), public, parameter :: eV2kcm=23.060542638_r8_kind
  real(kind=r8_kind), public, parameter :: kcm2eV=0.043364114_r8_kind
  real(kind=r8_kind), public, parameter :: eV2kJm=96.485310432_r8_kind
  real(kind=r8_kind), public, parameter :: kJm2eV=0.010364272_r8_kind

  real(kind=r8_kind), public , parameter :: h2erg = 4.35926565522e-11_r8_kind
  real(kind=r8_kind), public , parameter :: au2dyn = 8.237e-3_r8_kind
  real(kind=r8_kind), public , parameter :: avogadro = 6.0221367e23_r8_kind  !mol^(-1)
  real(kind=r8_kind), public , parameter :: pressure = 1.01325e6_r8_kind      !dyn/cm^2
  real(kind=r8_kind), public , parameter :: boltzmann = 1.380658e-16_r8_kind   !erg/K

  real(kind=r8_kind), public, parameter :: coulomb_factor=1389.354845_r8_kind
  real(kind=r8_kind), public, parameter :: dipole_factor=14.39418_r8_kind
  real(kind=r8_kind), public, parameter :: dip_eps=1.5_r8_kind

  real(kind=r8_kind), public, parameter :: rcf=1.2_r8_kind

  integer(kind=i4_kind), public, parameter :: len_name=6
  integer(kind=i4_kind), public, parameter :: n_parameter=15
  integer(kind=i4_kind), public, parameter :: n_poten=16 !!!!!!!!!!!!!
  integer(kind=i4_kind), public, parameter :: n_elements=98
  real(kind=r8_kind), public, parameter :: small0=1.0e-6_r8_kind
  real(kind=r8_kind), public, parameter :: small01=1.0e-7_r8_kind
  real(kind=r8_kind), public, parameter :: small=1.0e-8_r8_kind
  real(kind=r8_kind), public, parameter :: small1=1.0e-12_r8_kind
  real(kind=r8_kind), public, parameter :: small2=1.0e-20_r8_kind

  integer(i4_kind) ,public, parameter :: max_gxatoms=99
  real(kind=r8_kind), public, parameter :: gx_dummy_atom=99.0_r8_kind

  integer(i4_kind) ,public, parameter :: ierf=0
  integer(i4_kind) ,public, parameter :: ierfc=1

  real(kind=r8_kind), public, parameter :: r_ext_pc=2.2_r8_kind
  real(kind=r8_kind), public, parameter :: eps_h2o=78.4_r8_kind
  !------------ public functions and subroutines ------------------
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines ---------------------------------------
end module common_data_module

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
!=====================================================================
! Public interface of module
!=====================================================================
module efp_data_module
  !-------------------------------------------------------------------
  !
  !  Purpose: stores EFP-DFT  parameters for effective fragments
  !           (H2O currently). All data taken from US-GAMESS source code
  !           and were optimized for B3LYP XC potential and Dunning-Hay
  !           basis set
  !
  !  Module called by: ...
  !
  !
  !  References: US_GAMESS - file "efinp.src"
  ! 
  !
  !  Author: AS
  !  Date: 4/10/2007
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

  use type_module ! type specification parameters
  use common_data_module
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------
  type, public :: charge
     real(r8_kind) :: q
     real(r8_kind) :: C
     real(r8_kind) :: A
  end type charge

  type, public :: multpol
     character(len=12)  :: name
     real(r8_kind) :: coor(3)
     real(r8_kind) :: mass
     type(charge)  :: z
     real(r8_kind) :: dip(3)
     real(r8_kind) :: quad(3,3)
     real(r8_kind) :: oct(3,3,3)
  end type multpol

  type, public :: polarization
     character(len=12)  :: name
     real(r8_kind) :: coor(3)
     real(r8_kind) :: alpha(3,3)
  end type polarization

  type, public :: qm_ef_rep
     character(len=12)  :: name
     real(r8_kind) :: coor(3)
     real(r8_kind) :: C(2)
     real(r8_kind) :: A(2)
  end type qm_ef_rep

  type, public :: ef_ef_rep
     character(len=12)  :: name
     real(r8_kind) :: coor(3)
     real(r8_kind) :: C(4)
     real(r8_kind) :: A(4)
  end type ef_ef_rep

  !------------ Declaration of constants and variables ---------------
  integer(i4_kind), public, parameter :: n_gx_points=3,n_nuc=3, n_emp=5, n_pol=5, n_rep_q=3, n_rep_f=4
  real(r8_kind), public, parameter :: gx_num(3)=(/8.01, 1.02, 1.03/)
  integer(i4_kind), public, parameter :: n_water_centers=11

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------
  real(r8_kind), parameter :: C_i        =zero
  real(r8_kind), parameter :: A_i        =zero
  real(r8_kind), parameter :: D_i(3)     =zero
  real(r8_kind), parameter :: Q_i(3,3)   =zero
  real(r8_kind), parameter :: O_i(3,3,3) =zero
  real(r8_kind), parameter :: Al_i(3,3)  =zero

  !names of multipole
  character(len=12), parameter :: mname1="O1          "
  character(len=12), parameter :: mname2="H2          "
  character(len=12), parameter :: mname3="H3          "
  character(len=12), parameter :: mname1e="O1e         "
  character(len=12), parameter :: mname2e="H2e         "
  character(len=12), parameter :: mname3e="H3e         "
  character(len=12), parameter :: mname4="B12         "
  character(len=12), parameter :: mname5="B13         "

  !coordinates of multipole centers
  real(r8_kind), parameter :: mcoor1(3)=(/ zero,             zero, -0.119151_r8_kind/)
  real(r8_kind), parameter :: mcoor2(3)=(/-1.431042_r8_kind, zero,  0.945510_r8_kind/)
  real(r8_kind), parameter :: mcoor3(3)=(/ 1.431042_r8_kind, zero,  0.945510_r8_kind/)
  real(r8_kind), parameter :: mcoor4(3)=(/-0.715521_r8_kind, zero,  0.413179_r8_kind/)
  real(r8_kind), parameter :: mcoor5(3)=(/ 0.715521_r8_kind, zero,  0.413179_r8_kind/)

  !mass
  real(r8_kind), parameter :: mass1=15.99491_r8_kind
  real(r8_kind), parameter :: mass2=1.007825_r8_kind
  real(r8_kind), parameter :: mass3=1.007825_r8_kind

  !charges
  real(r8_kind), parameter :: Zn1=eight
  real(r8_kind), parameter :: Zn2=one
  real(r8_kind), parameter :: Zn3=one
  real(r8_kind), parameter :: Ze1=-8.2245784_r8_kind
  real(r8_kind), parameter :: Ze2=-0.5790554_r8_kind
  real(r8_kind), parameter :: Ze3=-0.5790554_r8_kind
  real(r8_kind), parameter :: Ze4=-0.3086554_r8_kind
  real(r8_kind), parameter :: Ze5=-0.3086554_r8_kind

  !dipoles
  real(r8_kind), parameter :: dip1(3)=(/ zero,             zero,  0.439368_r8_kind/)
  real(r8_kind), parameter :: dip2(3)=(/-0.045030_r8_kind, zero,  0.019745_r8_kind/)
  real(r8_kind), parameter :: dip3(3)=(/ 0.045030_r8_kind, zero,  0.019745_r8_kind/)
  real(r8_kind), parameter :: dip4(3)=(/ 0.151206_r8_kind, zero, -0.116204_r8_kind/)
  real(r8_kind), parameter :: dip5(3)=(/-0.151206_r8_kind, zero, -0.116204_r8_kind/)

  !quadrupoles
  real(r8_kind), parameter :: quad1(3,3)=reshape(                             &
       (/-3.545741_r8_kind,  zero,              zero,                         &
         zero,              -4.529667_r8_kind,  zero,                         &
         zero,               zero,             -4.022426_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: quad2(3,3)=reshape(                             &
       (/-0.305693_r8_kind,  zero,             -0.005554_r8_kind,             &
          zero,             -0.313535_r8_kind,  zero,                         &
         -0.005554_r8_kind,  zero,             -0.304098_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: quad3(3,3)=reshape(                             &
       (/-0.305693_r8_kind,  zero,              0.005554_r8_kind,             &
          zero,             -0.313535_r8_kind,  zero,                         &
          0.005554_r8_kind,  zero,             -0.304098_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: quad4(3,3)=reshape(                             &
       (/-0.096367_r8_kind,  zero,             -0.011666_r8_kind,             &
          zero,             -0.130668_r8_kind,  zero,                         &
         -0.011666_r8_kind,  zero,             -0.106838_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: quad5(3,3)=reshape(                             &
       (/-0.096367_r8_kind,  zero,              0.011666_r8_kind,             &
          zero,             -0.130668_r8_kind,  zero,                         &
          0.011666_r8_kind,  zero,             -0.106838_r8_kind/), (/3,3/))

  !octopoles
  real(r8_kind), parameter :: oct1(3,3,3)=reshape(                &
       (/ zero,              zero,              0.313291_r8_kind, &
          zero,              zero,              zero,             &
          0.313291_r8_kind,  zero,              zero,             &
          zero,              zero,              zero,             &
          zero,              zero,              0.360352_r8_kind, &
          zero,              0.360352_r8_kind,  zero,             &
          0.313291_r8_kind,  zero,              zero,             &
          zero,              0.360352_r8_kind,  zero,             &
          zero,              zero,              1.073886_r8_kind/), (/3,3,3/))
  real(r8_kind), parameter :: oct2(3,3,3)=reshape(                &
       (/-0.044409_r8_kind,  zero,             -0.000817_r8_kind, &
          zero,             -0.018142_r8_kind,  zero,             &
         -0.000817_r8_kind,  zero,             -0.012450_r8_kind, &
          zero,             -0.018142_r8_kind,  zero,             &
         -0.018142_r8_kind,  zero,              0.006012_r8_kind, &
          zero,              0.006012_r8_kind,  zero,             &
         -0.000817_r8_kind,  zero,             -0.012450_r8_kind, &
          zero,              0.006012_r8_kind,  zero,             &
         -0.012450_r8_kind,  zero,              0.012008_r8_kind/), (/3,3,3/))
  real(r8_kind), parameter :: oct3(3,3,3)=reshape(                &
       (/ 0.044409_r8_kind,  zero,             -0.000817_r8_kind, &
          zero,              0.018142_r8_kind,  zero,             &
         -0.000817_r8_kind,  zero,              0.012450_r8_kind, &
          zero,              0.018142_r8_kind,  zero,             &
          0.018142_r8_kind,  zero,              0.006012_r8_kind, &
          zero,              0.006012_r8_kind,  zero,             &
         -0.000817_r8_kind,  zero,              0.012450_r8_kind, &
          zero,              0.006012_r8_kind,  zero,             &
          0.012450_r8_kind,  zero,              0.012008_r8_kind/), (/3,3,3/))
  real(r8_kind), parameter :: oct4(3,3,3)=reshape(                &
       (/ 0.467428_r8_kind,  zero,             -0.109886_r8_kind, &
          zero,              0.153403_r8_kind,  zero,             &
         -0.109886_r8_kind,  zero,              0.156548_r8_kind, &
          zero,              0.153403_r8_kind,  zero,             &
          0.153403_r8_kind,  zero,             -0.113583_r8_kind, &
          zero,             -0.113583_r8_kind,  zero,             &
         -0.109886_r8_kind,  zero,              0.156548_r8_kind, &
          zero,             -0.113583_r8_kind,  zero,             &
          0.156548_r8_kind,  zero,             -0.336104_r8_kind/), (/3,3,3/))
  real(r8_kind), parameter :: oct5(3,3,3)=reshape(                &
       (/-0.467428_r8_kind,  zero,             -0.109886_r8_kind, &
          zero,             -0.153403_r8_kind,  zero,             &
         -0.109886_r8_kind,  zero,             -0.156548_r8_kind, &
          zero,             -0.153403_r8_kind,  zero,             &
         -0.153403_r8_kind,  zero,             -0.113583_r8_kind, &
          zero,             -0.113583_r8_kind,  zero,             &
         -0.109886_r8_kind,  zero,             -0.156548_r8_kind, &
          zero,             -0.113583_r8_kind,  zero,             &
         -0.156548_r8_kind,  zero,             -0.336104_r8_kind/), (/3,3,3/))

  !Charge screening QM-EF
  real(r8_kind), parameter :: CsQmEf1= 0.186119_r8_kind
  real(r8_kind), parameter :: AsQmEf1= 0.549105_r8_kind
  real(r8_kind), parameter :: CsQmEf2= 0.112182_r8_kind
  real(r8_kind), parameter :: AsQmEf2= 0.389541_r8_kind
  real(r8_kind), parameter :: CsQmEf3= 0.112182_r8_kind
  real(r8_kind), parameter :: AsQmEf3= 0.389541_r8_kind
  real(r8_kind), parameter :: CsQmEf4=-0.717580_r8_kind
  real(r8_kind), parameter :: AsQmEf4= 0.962143_r8_kind
  real(r8_kind), parameter :: CsQmEf5=-0.717580_r8_kind
  real(r8_kind), parameter :: AsQmEf5= 0.962143_r8_kind

  !Charge screening EF-EF
  real(r8_kind), parameter :: CsEfEf1=one
  real(r8_kind), parameter :: AsEfEf1=1.960183_r8_kind
  real(r8_kind), parameter :: CsEfEf2=one
  real(r8_kind), parameter :: AsEfEf2=2.383508_r8_kind
  real(r8_kind), parameter :: CsEfEf3=one
  real(r8_kind), parameter :: AsEfEf3=2.383508_r8_kind
  real(r8_kind), parameter :: CsEfEf4=one
  real(r8_kind), parameter :: AsEfEf4=9.999913_r8_kind
  real(r8_kind), parameter :: CsEfEf5=one
  real(r8_kind), parameter :: AsEfEf5=9.999913_r8_kind

  !names of polarizable centers
  character(len=12), parameter :: pname1="LMO1        "
  character(len=12), parameter :: pname2="LMO2        "
  character(len=12), parameter :: pname3="LMO3        "
  character(len=12), parameter :: pname4="LMO4        "
  character(len=12), parameter :: pname5="LMO5        "

  !coordinates of polarizable centers
  real(r8_kind), parameter :: pcoor1(3)=(/ zero,              zero,             -0.118741_r8_kind/)
  real(r8_kind), parameter :: pcoor2(3)=(/ 0.767899_r8_kind,  zero,              0.494658_r8_kind/)
  real(r8_kind), parameter :: pcoor3(3)=(/-0.767899_r8_kind,  zero,              0.494658_r8_kind/)
  real(r8_kind), parameter :: pcoor4(3)=(/ zero,             -0.492163_r8_kind, -0.404375_r8_kind/)
  real(r8_kind), parameter :: pcoor5(3)=(/ zero,              0.492163_r8_kind, -0.404375_r8_kind/)

  !tensors of polarization
  real(r8_kind), parameter :: pol_ten1(3,3)=reshape(              &
       (/ 0.003101_r8_kind,  zero,              zero,             &
          zero,              0.004553_r8_kind,  zero,             &
          zero,              zero,              0.002822_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: pol_ten2(3,3)=reshape(              &
       (/ 2.109350_r8_kind,  zero,              0.934523_r8_kind, &
          zero,              0.835171_r8_kind,  zero,             &
          1.313981_r8_kind,  zero,              1.556153_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: pol_ten3(3,3)=reshape(              &
       (/ 2.109350_r8_kind,  zero,             -0.934523_r8_kind, &
          zero,              0.835171_r8_kind,  zero,             &
         -1.313981_r8_kind,  zero,              1.556153_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: pol_ten4(3,3)=reshape(              &
       (/ 1.519500_r8_kind,  zero            ,  zero,             &
          zero            ,  0.736181_r8_kind,  0.089060_r8_kind, &
          zero,              0.780917_r8_kind,  1.183064_r8_kind/), (/3,3/))
  real(r8_kind), parameter :: pol_ten5(3,3)=reshape(              &
       (/ 1.519500_r8_kind,  zero            ,  zero,             &
          zero,              0.736181_r8_kind, -0.089060_r8_kind, &
          zero,             -0.780917_r8_kind,  1.183064_r8_kind/), (/3,3/))

  !names of repulsive centers
  character(len=12), parameter :: rname1=mname1
  character(len=12), parameter :: rname2=mname2
  character(len=12), parameter :: rname3=mname3
  character(len=12), parameter :: rname4="CMS         "

  !coordinates of repulsive centers
  real(r8_kind), parameter :: rcoor1(3)=mcoor1
  real(r8_kind), parameter :: rcoor2(3)=mcoor2
  real(r8_kind), parameter :: rcoor3(3)=mcoor3
  real(r8_kind), parameter :: rcoor4(3)=(/zero, zero, zero/) ! center of masses

  !repulsion QM-EF
  real(r8_kind), parameter :: CrQmEf1(2)=(/-0.0012471382_r8_kind,17.1604999173_r8_kind/)
  real(r8_kind), parameter :: ArQmEf1(2)=(/ 0.0761496713_r8_kind, 0.9999999904_r8_kind/)
  real(r8_kind), parameter :: CrQmEf2(2)=(/-0.0002488372_r8_kind, 0.1063437802_r8_kind/)
  real(r8_kind), parameter :: ArQmEf2(2)=(/ 0.0949999961_r8_kind, 0.5999987638_r8_kind/)
  real(r8_kind), parameter :: CrQmEf3(2)=(/-0.0002488372_r8_kind, 0.1063437802_r8_kind/)
  real(r8_kind), parameter :: ArQmEf3(2)=(/ 0.0949999961_r8_kind, 0.5999987638_r8_kind/)
  
  !repulsion EF-EF
  real(r8_kind), parameter :: CrEfEf1(4)=(/-47.204255_r8_kind, 44.184769_r8_kind, 44.184769_r8_kind, 31.732616_r8_kind/)
  real(r8_kind), parameter :: ArEfEf1(4)=(/  1.502333_r8_kind,  2.731355_r8_kind,  2.731355_r8_kind,  1.518409_r8_kind/)
  real(r8_kind), parameter :: CrEfEf2(4)=(/ 44.184769_r8_kind,  0.120699_r8_kind,  0.120699_r8_kind, -0.061708_r8_kind/)
  real(r8_kind), parameter :: ArEfEf2(4)=(/  2.731355_r8_kind,  1.315321_r8_kind,  1.315321_r8_kind,  0.802301_r8_kind/)
  real(r8_kind), parameter :: CrEfEf3(4)=(/ 44.184769_r8_kind,  0.120699_r8_kind,  0.120699_r8_kind, -0.061708_r8_kind/)
  real(r8_kind), parameter :: ArEfEf3(4)=(/  2.731355_r8_kind,  1.315321_r8_kind,  1.315321_r8_kind,  0.802301_r8_kind/)
  real(r8_kind), parameter :: CrEfEf4(4)=(/ 31.732616_r8_kind, -0.061708_r8_kind, -0.061708_r8_kind,146.763464_r8_kind/)
  real(r8_kind), parameter :: ArEfEf4(4)=(/  1.518409_r8_kind,  0.802301_r8_kind,  0.802301_r8_kind,  2.083947_r8_kind/)

  integer(i4_kind) :: i

  !-------------------------------------------------------------------
  type(multpol), public :: dft_water_nuc(n_nuc)=(/                     &
                                multpol(                               &
                                        mname1,mcoor1,mass1,           &
                                        charge(Zn1,C_i,A_i),           &
                                        D_i,Q_i,O_i                    &
                                       ),                              &
                                multpol(                               &
                                        mname2,mcoor2,mass2,           &
                                        charge(Zn2,C_i,A_i),           &
                                        D_i,Q_i,O_i                    &
                                       ),                              &
                                multpol(                               &
                                        mname3,mcoor3,mass3,           &
                                        charge(Zn3,C_i,A_i),           &
                                        D_i,Q_i,O_i                    &
                                       )                               &
                                /)

  type(multpol), public :: dft_water_nuc_f(n_nuc)=(/                   &
                                multpol(                               &
                                        mname1,mcoor1,mass1,           &
                                        charge(Zn1,C_i,A_i),           &
                                        D_i,Q_i,O_i                    &
                                       ),                              &
                                multpol(                               &
                                        mname2,mcoor2,mass2,           &
                                        charge(Zn2,C_i,A_i),           &
                                        D_i,Q_i,O_i                    &
                                       ),                              &
                                multpol(                               &
                                        mname3,mcoor3,mass3,           &
                                        charge(Zn3,C_i,A_i),           &
                                        D_i,Q_i,O_i                    &
                                       )                               &
                                /)

  type(multpol), public :: dft_water_mp(n_emp)=(/                      &
                                multpol(                               &
                                        mname1e,mcoor1,zero,           &
                                        charge(Ze1,CsQmEf1,AsQmEf1),   &
                                        dip1,quad1,oct1                &
                                       ),                              &
                                multpol(                               &
                                        mname2e,mcoor2,zero,           &
                                        charge(Ze2,CsQmEf2,AsQmEf2),   &
                                        dip2,quad2,oct2                &
                                       ),                              &
                                multpol(                               &
                                        mname3e,mcoor3,zero,           &
                                        charge(Ze3,CsQmEf3,AsQmEf3),   &
                                        dip3,quad3,oct3                &
                                       ),                              &
                                multpol(                               &
                                        mname4,mcoor4,zero,            &
                                        charge(Ze4,CsQmEf4,AsQmEf4),   &
                                        dip4,quad4,oct4                &
                                       ),                              &
                                multpol(                               &
                                        mname5,mcoor5,zero,            &
                                        charge(Ze5,CsQmEf5,AsQmEf5),   &
                                        dip5,quad5,oct5                &
                                       )                               &
                                /)

  type(multpol), public :: dft_water_mp_f(n_emp)=(/                    &
                                multpol(                               &
                                        mname1e,mcoor1,zero,           &
                                        charge(Ze1,CsEfEf1,AsEfEf1),   &
                                        dip1,quad1,oct1                &
                                       ),                              &
                                multpol(                               &
                                        mname2e,mcoor2,zero,           &
                                        charge(Ze2,CsEfEf2,AsEfEf2),   &
                                        dip2,quad2,oct2                &
                                       ),                              &
                                multpol(                               &
                                        mname3e,mcoor3,zero,           &
                                        charge(Ze3,CsEfEf3,AsEfEf3),   &
                                        dip3,quad3,oct3                &
                                       ),                              &
                                multpol(                               &
                                        mname4,mcoor4,zero,            &
                                        charge(Ze4,CsEfEf4,AsEfEf4),   &
                                        dip4,quad4,oct4                &
                                       ),                              &
                                multpol(                               &
                                        mname5,mcoor5,zero,            &
                                        charge(Ze5,CsEfEf5,AsEfEf5),   &
                                        dip5,quad5,oct5                &
                                       )                               &
                                /)

  type(polarization), public :: dft_water_pol(n_pol)=(/               &
                                polarization(pname1,pcoor1,pol_ten1), &
                                polarization(pname2,pcoor2,pol_ten2), &
                                polarization(pname3,pcoor3,pol_ten3), &
                                polarization(pname4,pcoor4,pol_ten4), &
                                polarization(pname5,pcoor5,pol_ten5)  &
                                /)

  type(qm_ef_rep), public :: dft_water_qrep(n_rep_q)=(/                    &
                                 qm_ef_rep(rname1,rcoor1,CrQmEf1,ArQmEf1), &
                                 qm_ef_rep(rname2,rcoor2,CrQmEf2,ArQmEf2), &
                                 qm_ef_rep(rname3,rcoor3,CrQmEf3,ArQmEf3)  &
                                 /)

  type(ef_ef_rep), public :: dft_water_frep(n_rep_f)=(/                    &
                                 ef_ef_rep(rname1,rcoor1,CrEfEf1,ArEfEf1), &
                                 ef_ef_rep(rname2,rcoor2,CrEfEf2,ArEfEf2), &
                                 ef_ef_rep(rname3,rcoor3,CrEfEf3,ArEfEf3), &
                                 ef_ef_rep(rname4,rcoor4,CrEfEf4,ArEfEf4)  &
                                 /)
  !------------ Subroutines ------------------------------------------

  !--------------- End of module -------------------------------------
end module efp_data_module

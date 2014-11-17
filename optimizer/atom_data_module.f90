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
module  atom_data_module
  !-------------------------------------------------------------------
  !
  !  Purpose: contains all data for atoms: Symbols, Masses ...
  !
  !  Module called by: ...
  !
  !  References: Periodic table of Elements next to Room 63405 (masses)
  !              Program 'NABOR' by T.A.Halgren, City University of 
  !              NewYork (covalent radii)
  !
  !  Author: FN
  !  Date: ...
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !-------------------------------------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Declaration of types ---------------------------------
  
  !------------ Declaration of constants and variables ---------------
  integer(i4_kind)  ,parameter :: NELEM=103
  character(len=10) ,public    :: symbol(NELEM)
  real(kind=r8_kind)           :: mass(NELEM)
  real(kind=r8_kind),public    :: covrad(NELEM)

  real(kind=r8_kind),public,parameter   :: dummy_charge=99.0_r8_kind
  real(kind=r8_kind),public,parameter   :: dummy_mass=0.0_r8_kind
  
  data      mass/1.007825_r8_kind, 4.002603_r8_kind, 7.016005_r8_kind, 9.012182_r8_kind,&
       11.009305_r8_kind,12.000000_r8_kind,14.003074_r8_kind,15.994915_r8_kind, &
       18.998403_r8_kind,19.992439_r8_kind,22.989770_r8_kind,23.985045_r8_kind, &
       26.981541_r8_kind,27.976928_r8_kind,30.973763_r8_kind,31.973909_r8_kind, &
       34.968853_r8_kind,39.962383_r8_kind,38.963708_r8_kind,39.962591_r8_kind, &
       44.955914_r8_kind,47.947947_r8_kind,50.943963_r8_kind,51.940510_r8_kind, &
       54.938046_r8_kind,55.934939_r8_kind,58.933198_r8_kind,57.935347_r8_kind, &
       62.929599_r8_kind,63.929145_r8_kind,68.925581_r8_kind,73.921179_r8_kind, &
       74.921596_r8_kind,79.916521_r8_kind,78.918336_r8_kind,83.911506_r8_kind, &
       84.911800_r8_kind,87.956250_r8_kind,88.905856_r8_kind,89.904708_r8_kind, &
       92.906378_r8_kind,97.905405_r8_kind,97.907110_r8_kind,101.904348_r8_kind,&
       102.90550_r8_kind,106.903480_r8_kind,106.905095_r8_kind,113.903361_r8_kind, &
       114.90388_r8_kind,119.902199_r8_kind,120.903824_r8_kind,129.906230_r8_kind, &
       126.904477_r8_kind,131.90415_r8_kind,132.905770_r8_kind,137.905240_r8_kind, &
       138.906360_r8_kind,139.90544_r8_kind,140.907660_r8_kind,141.907730_r8_kind, &
       144.912691_r8_kind,151.91974_r8_kind,152.921240_r8_kind,157.924110_r8_kind, &
       158.92535_r8_kind,163.929180_r8_kind,164.930330_r8_kind,167.930310_r8_kind, &
       168.93423_r8_kind,173.938870_r8_kind,174.940790_r8_kind,177.943710_r8_kind, &
       180.94801_r8_kind,183.950950_r8_kind,186.955770_r8_kind,191.961490_r8_kind, &
       192.96294_r8_kind,194.964790_r8_kind,196.966560_r8_kind,201.970630_r8_kind, &
       204.97441_r8_kind,207.976640_r8_kind,208.980390_r8_kind,208.982420_r8_kind, &
       219.01130_r8_kind,222.017570_r8_kind,223.019734_r8_kind,226.025406_r8_kind, &
       227.02775_r8_kind,232.038050_r8_kind,231.035881_r8_kind,238.050786_r8_kind, &
       237.048169_r8_kind,244.06420_r8_kind,243.0614_r8_kind  ,247.0703_r8_kind  , &
       247.0703_r8_kind  ,251.0796_r8_kind ,252.0829_r8_kind  ,257.0951_r8_kind  , &
       258.0986_r8_kind  ,259.1009_r8_kind ,262.1100_r8_kind/
       ! older versions used mass 250 for elements > Np and < 99

  data symbol /'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na', &
       'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr', &
       'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In', &
       'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
       'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re', &
       'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
       'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','X ','Fm',      &
       'Md','No','Lr'/
       ! FIXME: X at 99 is in fact Es!

  data covrad /0.3200_r8_kind,0.9300_r8_kind,1.2300_r8_kind,0.9000_r8_kind,0.8200_r8_kind,&
       0.7700_r8_kind,0.7500_r8_kind,0.7300_r8_kind,0.7200_r8_kind,0.7100_r8_kind, &
       1.5400_r8_kind,1.3600_r8_kind,1.1800_r8_kind,1.1100_r8_kind,1.0600_r8_kind,&
       1.0200_r8_kind,0.9900_r8_kind,0.9800_r8_kind,2.0300_r8_kind,1.7400_r8_kind,&
       1.4400_r8_kind,1.3200_r8_kind,1.2200_r8_kind,1.1800_r8_kind,1.1700_r8_kind,&
       1.1700_r8_kind,1.1600_r8_kind,1.1500_r8_kind,1.1700_r8_kind,1.2500_r8_kind,&
       1.2600_r8_kind,1.2200_r8_kind,1.2000_r8_kind,1.1600_r8_kind,1.1400_r8_kind,&
       1.1200_r8_kind,2.1600_r8_kind,1.9100_r8_kind,1.6200_r8_kind,1.4500_r8_kind,&
       1.3400_r8_kind,1.3000_r8_kind,1.2700_r8_kind,1.2500_r8_kind,1.2500_r8_kind,&
       1.2800_r8_kind,1.3400_r8_kind,1.4800_r8_kind,1.4400_r8_kind,1.4100_r8_kind,&
       1.4000_r8_kind,1.3600_r8_kind,1.3300_r8_kind,1.3100_r8_kind,2.3500_r8_kind,&
       1.9800_r8_kind,1.6900_r8_kind,1.6500_r8_kind,1.6500_r8_kind,1.6400_r8_kind,&
       1.6300_r8_kind,1.6200_r8_kind,1.8500_r8_kind,1.6100_r8_kind,1.5900_r8_kind,&
       1.5900_r8_kind,1.5800_r8_kind,1.5700_r8_kind,1.5600_r8_kind,1.5600_r8_kind,&
       1.5600_r8_kind,1.4400_r8_kind,1.3400_r8_kind,1.3000_r8_kind,1.2800_r8_kind,&
       1.2600_r8_kind,1.2700_r8_kind,1.3000_r8_kind,1.3400_r8_kind,1.4900_r8_kind,&
       1.4800_r8_kind,1.4700_r8_kind,1.4600_r8_kind,1.4600_r8_kind,1.4500_r8_kind,&
       1.0000_r8_kind,1.0000_r8_kind,1.0000_r8_kind,1.0000_r8_kind,1.6500_r8_kind,&
       1.0000_r8_kind,1.4200_r8_kind,1.0000_r8_kind,1.0000_r8_kind,1.0000_r8_kind,&
       1.0000_r8_kind,1.0000_r8_kind,1.0000_r8_kind,1.0000_r8_kind,0.8000_r8_kind,&
       1.0000_r8_kind,1.0000_r8_kind,1.0000_r8_kind/
  !------------ public functions and subroutines ---------------------
  public :: get_row
  public :: nuc_mass
  public :: nuc_radius
  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables ---------------
  
  !-------------------------------------------------------------------
contains
  !------------ Subroutines ------------------------------------------
  function get_row(number)
    ! Purpose: return the row of the periodic table where the
    !          element with atomic number 'number' is located.
    !---------------------------------------------------------
    integer(kind=i4_kind) :: get_row
    integer(kind=i4_kind) :: number
    !---- Declaration of local variables ---------------------
    
    if (number==1.or.number==2) then
       get_row=1
    elseif (number>=3.and.number<=10) then
       get_row=2
    elseif (number>=11.and.number<=18) then
       get_row=3
    elseif (number>=19.and.number<=36) then
       get_row=4
    elseif (number>=37.and.number<=54) then
       get_row=5
    elseif (number>=55.and.number<=86) then   
       get_row=6
    elseif (number>=87.and.number<=103) then 
       get_row=7
    else
       stop ' get_row: rubbish entered'
    endif

  end function get_row

  function nuc_mass(Z) result(M)
    use comm, only: comm_rank
    implicit none
    integer(kind=i4_kind), intent(in) :: Z
    real(r8_kind)                     :: M ! result
    !---- Declaration of local variables ---------------------

    if(Z>0 .and. Z<=NELEM)then
       M = mass(Z)
    else
       if(comm_rank() == 0) then
          WARN('Z out of range')
       endif
       M = 999.0_r8_kind
    endif
    if(comm_rank() == 0) then
       DPRINT 'nuc_mass(',Z,')=',M
    endif
  end function nuc_mass

  function nuc_radius(Z) result(R)
    ! according to Mayers PHD
    implicit none
    integer(kind=i4_kind), intent(in) :: Z
    real(r8_kind)                     :: R ! result
    !---- Declaration of local variables ---------------------

    real(r8_kind), parameter :: third = 1.0_r8_kind/3.0_r8_kind
    real(r8_kind) :: M                
    
    M = nuc_mass(Z)
    R = 2.27E-05_r8_kind * M**third
    DPRINT 'nuc_radius(',Z,')=',R
  end function nuc_radius

  !******************************************************************
       
!--------------- End of module ----------------------------------
end module atom_data_module

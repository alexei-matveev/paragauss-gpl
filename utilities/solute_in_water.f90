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
program solute_in_water

  implicit none
!===============================================
  type :: water_atom
     character*2 :: name
     real*8 :: coor(3)
     real*8 :: r_vdw
  end type water_atom

  type :: efp_water
     type(water_atom) :: atom1,atom2,atom3
  end type efp_water

  integer*4, parameter :: n_efp=216

  type(efp_water) :: water_array0(n_efp)   !TIP4P
  type(efp_water) :: water_array1(n_efp)  !EFP

  type(efp_water), allocatable ::  water_array(:) !working array

  real*8,parameter :: box_size=18.6206d0

  character*2, parameter :: O1="O1", H2="H2", H3="H3"
  real*8, parameter :: r_O1=1.52d0, r_H2=1.20d0, r_H3=1.20d0
  real*8 :: EFP_O1(3)
  real*8 :: EFP_H2(3)
  real*8 :: EFP_H3(3)

  type :: vdw_data
     character*2 :: name
     real*8 :: radius
  end type vdw_data

  integer*4, parameter :: n_atoms=98
  type(vdw_data) :: r_vdW(n_atoms)

  integer*4, parameter :: max_solute_atoms=100

  type :: solute_atom
     character*2 :: name
     real*8 :: coor(3)
     real*8 :: r_vdw
  end type solute_atom

  type(solute_atom) :: solute(max_solute_atoms)
  integer*4 :: n_solute_atoms
  real*8 :: g_solu_center(3), center(3)

  type(efp_water), allocatable :: water_shell(:)
  real*8, allocatable :: water_center(:,:)

  real*8 :: R_sol, increment
  integer*4 :: n_waters
  character*3 :: solv !TIP, EFP
  character*1 :: w_atom ! A(Atom), M(Molecule)
  character*1 :: pos ! QM cluster position: C(Center), D(ranDom)

  integer*4,parameter :: max_param=6

  integer*4 :: n_waters_in_shell
  integer*4 :: n_efp_full

  integer*4 :: i,status

  character*100 :: buffer
  character*4 :: buf
  integer*4 :: n_param
  integer*4 :: pos1(max_param-1), pos2(1)
  integer*4 :: i_begin,i_end
!===============================================

  data (water_array0(i), i=1,20) &
       /                                                    &
!1
         efp_water(water_atom(O1, (/ 2.300,  6.280,  1.130/), r_O1) , &
                   water_atom(H2, (/ 1.370,  6.260,  1.500/), r_H2) , &
                   water_atom(H3, (/ 2.310,  5.890,  0.210/), r_H3)),  &
!2
         efp_water(water_atom(O1, (/ 2.250,  2.750, -8.660/), r_O1) , &
                   water_atom(H2, (/ 2.600,  2.580, -7.740/), r_H2) , &
                   water_atom(H3, (/ 1.370,  2.300, -8.780/), r_H3)), &
!3
         efp_water(water_atom(O1, (/ 0.190,  3.680,  6.470/), r_O1) , &
                   water_atom(H2, (/-0.630,  4.110,  6.860/), r_H2) , &
                   water_atom(H3, (/-0.090,  2.950,  5.840/), r_H3)), &
!4
         efp_water(water_atom(O1, (/ 5.690, -5.870, -6.970/), r_O1) , &
                   water_atom(H2, (/ 4.760, -5.940, -7.340/), r_H2) , &
                   water_atom(H3, (/ 5.800, -4.980, -6.530/), r_H3)), &
!5
         efp_water(water_atom(O1, (/-3.070, -3.510,  7.030/), r_O1) , &
                   water_atom(H2, (/-3.640, -3.670,  7.840/), r_H2) , &
                   water_atom(H3, (/-3.660, -3.410,  6.230/), r_H3)), &
!6
         efp_water(water_atom(O1, (/-1.190,  6.180,  8.560/), r_O1) , &
                   water_atom(H2, (/-0.860,  7.120,  8.560/), r_H2) , &
                   water_atom(H3, (/-0.680,  5.640,  9.220/), r_H3)), &
!7
         efp_water(water_atom(O1, (/-7.270,  7.030,  7.170/), r_O1) , &
                   water_atom(H2, (/-6.700,  7.810,  6.920/), r_H2) , &
                   water_atom(H3, (/-7.870,  7.290,  7.930/), r_H3)), &
!8
         efp_water(water_atom(O1, (/-1.070,  6.070,  2.310/), r_O1) , &
                   water_atom(H2, (/-1.190,  5.940,  1.320/), r_H2) , &
                   water_atom(H3, (/-1.370,  5.260,  2.800/), r_H3)), &
!9
         efp_water(water_atom(O1, (/ 7.680, -7.180, -8.390/), r_O1) , &
                   water_atom(H2, (/ 6.900, -7.010, -7.790/), r_H2) , &
                   water_atom(H3, (/ 8.020, -6.310, -8.750/), r_H3)), &
!10
         efp_water(water_atom(O1, (/ 8.500,  7.980, -0.390/), r_O1) , &
                   water_atom(H2, (/ 8.460,  8.740,  0.260/), r_H2) , &
                   water_atom(H3, (/ 8.720,  8.340, -1.300/), r_H3)), &
!11
         efp_water(water_atom(O1, (/ 6.850, -8.500,  6.650/), r_O1) , &
                   water_atom(H2, (/ 7.540, -8.660,  7.350/), r_H2) , &
                   water_atom(H3, (/ 6.120, -7.930,  7.030/), r_H3)), &
!12
         efp_water(water_atom(O1, (/ 6.860, -7.010, -0.590/), r_O1) , &
                   water_atom(H2, (/ 7.460, -6.220, -0.450/), r_H2) , &
                   water_atom(H3, (/ 6.000, -6.700, -1.000/), r_H3)), &
!13
         efp_water(water_atom(O1, (/ 3.350, -4.270, -8.010/), r_O1) , &
                   water_atom(H2, (/ 2.570, -4.580, -8.540/), r_H2) , &
                   water_atom(H3, (/ 3.930, -3.690, -8.580/), r_H3)), &
!14
         efp_water(water_atom(O1, (/-4.020, -3.570, -5.230/), r_O1) , &
                   water_atom(H2, (/-3.780, -2.630, -4.970/), r_H2) , &
                   water_atom(H3, (/-4.180, -4.110, -4.410/), r_H3)), &
!15
         efp_water(water_atom(O1, (/ 4.380,  3.920, -3.630/), r_O1) , &
                   water_atom(H2, (/ 5.200,  3.360, -3.540/), r_H2) , &
                   water_atom(H3, (/ 3.570,  3.340, -3.590/), r_H3)), &
!16
         efp_water(water_atom(O1, (/-2.590,  4.470,  7.370/), r_O1) , &
                   water_atom(H2, (/-3.330,  4.930,  6.870/), r_H2) , &
                   water_atom(H3, (/-2.080,  5.150,  7.900/), r_H3)), &
!17
         efp_water(water_atom(O1, (/ 2.310, -1.490,  4.830/), r_O1) , &
                   water_atom(H2, (/ 2.650, -0.720,  5.370/), r_H2) , &
                   water_atom(H3, (/ 2.750, -1.490,  3.930/), r_H3)), &
!18
         efp_water(water_atom(O1, (/-7.350, -5.210, -1.720/), r_O1) , &
                   water_atom(H2, (/-6.880, -5.210, -0.840/), r_H2) , &
                   water_atom(H3, (/-7.830, -6.080, -1.830/), r_H3)), &
!19
         efp_water(water_atom(O1, (/ 2.300, -4.280,  5.380/), r_O1) , &
                   water_atom(H2, (/ 2.040, -3.320,  5.380/), r_H2) , &
                   water_atom(H3, (/ 1.590, -4.820,  5.830/), r_H3)), &
!20
         efp_water(water_atom(O1, (/ 2.400, -7.710,  8.860/), r_O1) , &
                   water_atom(H2, (/ 2.540, -8.550,  9.380/), r_H2) , &
                   water_atom(H3, (/ 1.850, -7.070,  9.410/), r_H3))  &
       /

  data (water_array0(i), i=21,40) &
       /                                                     &
!21
         efp_water(water_atom(O1, (/ 6.200, -0.760, -4.230/), r_O1) , &
                   water_atom(H2, (/ 5.280, -0.930, -3.880/), r_H2) , &
                   water_atom(H3, (/ 6.480,  0.160, -3.970/), r_H3)), &
!22
         efp_water(water_atom(O1, (/ 6.060, -8.980,  1.230/), r_O1) , &
                   water_atom(H2, (/ 6.130, -8.140,  0.690/), r_H2) , &
                   water_atom(H3, (/ 6.520, -8.850,  2.110/), r_H3)), &
!23
         efp_water(water_atom(O1, (/-2.680,  1.140, -3.820/), r_O1) , &
                   water_atom(H2, (/-2.860,  1.810, -4.540/), r_H2) , &
                   water_atom(H3, (/-2.710,  1.600, -2.930/), r_H3)), &
!24
         efp_water(water_atom(O1, (/ 1.220,  6.430,  5.630/), r_O1) , &
                   water_atom(H2, (/ 0.770,  5.550,  5.800/), r_H2) , &
                   water_atom(H3, (/ 1.210,  6.970,  6.470/), r_H3)), &
!25
         efp_water(water_atom(O1, (/-0.200, -0.950,  3.590/), r_O1) , &
                   water_atom(H2, (/ 0.340, -1.240,  4.390/), r_H2) , &
                   water_atom(H3, (/ 0.100, -0.050,  3.300/), r_H3)), &
!26
         efp_water(water_atom(O1, (/ 0.270, -2.660,  1.170/), r_O1) , &
                   water_atom(H2, (/ 0.080, -3.620,  1.380/), r_H2) , &
                   water_atom(H3, (/-0.060, -2.080,  1.920/), r_H3)), &
!27
         efp_water(water_atom(O1, (/-1.730,  9.220,  6.120/), r_O1) , &
                   water_atom(H2, (/-0.780,  8.930,  6.200/), r_H2) , &
                   water_atom(H3, (/-1.810,  9.870,  5.370/), r_H3)), &
!28
         efp_water(water_atom(O1, (/-2.210, -7.540,  4.320/), r_O1) , &
                   water_atom(H2, (/-1.350, -7.520,  3.800/), r_H2) , &
                   water_atom(H3, (/-2.070, -7.070,  5.200/), r_H3)), &
!29
         efp_water(water_atom(O1, (/ 1.130,  7.370, -2.650/), r_O1) , &
                   water_atom(H2, (/ 2.010,  7.240, -2.200/), r_H2) , &
                   water_atom(H3, (/ 1.000,  8.340, -2.870/), r_H3)), &
!30
         efp_water(water_atom(O1, (/ 6.130, -4.970,  7.260/), r_O1) , &
                   water_atom(H2, (/ 5.640, -5.840,  7.350/), r_H2) , &
                   water_atom(H3, (/ 5.900, -4.540,  6.390/), r_H3)), &
!31
         efp_water(water_atom(O1, (/-5.690, -6.340, -4.390/), r_O1) , &
                   water_atom(H2, (/-5.320, -7.070, -4.970/), r_H2) , &
                   water_atom(H3, (/-5.170, -6.290, -3.540/), r_H3)), &
!32
         efp_water(water_atom(O1, (/ 8.090,  0.040,  5.020/), r_O1) , &
                   water_atom(H2, (/ 8.490,  0.950,  4.930/), r_H2) , &
                   water_atom(H3, (/ 7.090,  0.120,  5.080/), r_H3)), &
!33
         efp_water(water_atom(O1, (/ 1.970, -8.860, -5.980/), r_O1) , &
                   water_atom(H2, (/ 2.860, -9.310, -6.120/), r_H2) , &
                   water_atom(H3, (/ 1.240, -9.510, -6.170/), r_H3)), &
!34
         efp_water(water_atom(O1, (/-3.370, -8.630,  1.900/), r_O1) , &
                   water_atom(H2, (/-4.000, -9.390,  2.030/), r_H2) , &
                   water_atom(H3, (/-2.890, -8.450,  2.760/), r_H3)), &
!35
         efp_water(water_atom(O1, (/-6.750, -0.700, -2.460/), r_O1) , &
                   water_atom(H2, (/-6.510, -0.100, -3.220/), r_H2) , &
                   water_atom(H3, (/-6.680, -1.650, -2.760/), r_H3)), &
!36
         efp_water(water_atom(O1, (/ 3.170,  2.510, -0.610/), r_O1) , &
                   water_atom(H2, (/ 3.880,  3.220, -0.550/), r_H2) , &
                   water_atom(H3, (/ 2.290,  2.900, -0.330/), r_H3)), &
!37
         efp_water(water_atom(O1, (/-3.960, -4.450, -9.090/), r_O1) , &
                   water_atom(H2, (/-4.550, -4.390, -8.290/), r_H2) , &
                   water_atom(H3, (/-4.110, -5.330, -9.550/), r_H3)), &
!38
         efp_water(water_atom(O1, (/-1.950, -1.480,  5.720/), r_O1) , &
                   water_atom(H2, (/-2.360, -1.710,  4.840/), r_H2) , &
                   water_atom(H3, (/-2.130, -2.220,  6.370/), r_H3)), &
!39
         efp_water(water_atom(O1, (/ 5.980,  7.290,  2.700/), r_O1) , &
                   water_atom(H2, (/ 6.220,  7.980,  2.020/), r_H2) , &
                   water_atom(H3, (/ 5.200,  7.620,  3.240/), r_H3)), &
!40
         efp_water(water_atom(O1, (/-5.810,  3.450, -9.180/), r_O1) , &
                   water_atom(H2, (/-6.670,  2.950, -9.310/), r_H2) , &
                   water_atom(H3, (/-5.190,  2.910, -8.620/), r_H3))  &
       /

  data (water_array0(i), i=41,60) &
       /                                                     &
!41
         efp_water(water_atom(O1, (/-2.860, -2.000,  3.070/), r_O1) , &
                   water_atom(H2, (/-1.970, -1.540,  3.100/), r_H2) , &
                   water_atom(H3, (/-3.070, -2.240,  2.120/), r_H3)), &
!42
         efp_water(water_atom(O1, (/ 8.070,  6.050, -3.970/), r_O1) , &
                   water_atom(H2, (/ 7.600,  6.020, -3.080/), r_H2) , &
                   water_atom(H3, (/ 7.560,  5.500, -4.630/), r_H3)), &
!43
         efp_water(water_atom(O1, (/-4.680,  4.690, -1.880/), r_O1) , &
                   water_atom(H2, (/-4.880,  5.120, -1.000/), r_H2) , &
                   water_atom(H3, (/-3.900,  4.070, -1.790/), r_H3)), &
!44
         efp_water(water_atom(O1, (/-8.890,  8.900, -2.900/), r_O1) , &
                   water_atom(H2, (/-8.430,  8.060, -3.190/), r_H2) , &
                   water_atom(H3, (/-9.450,  9.240, -3.650/), r_H3)), &
!45
         efp_water(water_atom(O1, (/-8.710,  4.100, -6.200/), r_O1) , &
                   water_atom(H2, (/-9.480,  4.440, -5.660/), r_H2) , &
                   water_atom(H3, (/-9.050,  3.590, -6.990/), r_H3)), &
!46
         efp_water(water_atom(O1, (/-8.210,  7.010,  4.290/), r_O1) , &
                   water_atom(H2, (/-7.950,  6.970,  5.250/), r_H2) , &
                   water_atom(H3, (/-9.060,  6.500,  4.150/), r_H3)), &
!47
         efp_water(water_atom(O1, (/ 0.760,  8.110,  7.890/), r_O1) , &
                   water_atom(H2, (/ 1.750,  7.990,  7.980/), r_H2) , &
                   water_atom(H3, (/ 0.520,  9.060,  8.100/), r_H3)), &
!48
         efp_water(water_atom(O1, (/ 1.300, -0.410, -2.910/), r_O1) , &
                   water_atom(H2, (/ 1.200, -0.560, -1.920/), r_H2) , &
                   water_atom(H3, (/ 0.440, -0.050, -3.270/), r_H3)), &
!49
         efp_water(water_atom(O1, (/ 8.650,  3.480,  1.950/), r_O1) , &
                   water_atom(H2, (/ 9.240,  4.110,  1.460/), r_H2) , &
                   water_atom(H3, (/ 8.840,  2.540,  1.660/), r_H3)), &
!50
         efp_water(water_atom(O1, (/-1.430,  5.850, -0.310/), r_O1) , &
                   water_atom(H2, (/-1.690,  6.740, -0.670/), r_H2) , &
                   water_atom(H3, (/-1.450,  5.170, -1.040/), r_H3)), &
!51
         efp_water(water_atom(O1, (/-5.000, -7.180,  5.450/), r_O1) , &
                   water_atom(H2, (/-4.170, -7.470,  4.970/), r_H2) , &
                   water_atom(H3, (/-5.490, -6.510,  4.890/), r_H3)), &
!52
         efp_water(water_atom(O1, (/ 5.500,  1.960,  8.850/), r_O1) , &
                   water_atom(H2, (/ 5.450,  1.910,  9.850/), r_H2) , &
                   water_atom(H3, (/ 5.520,  2.920,  8.560/), r_H3)), &
!53
         efp_water(water_atom(O1, (/-8.540, -4.060,  4.770/), r_O1) , &
                   water_atom(H2, (/-9.000, -3.340,  4.250/), r_H2) , &
                   water_atom(H3, (/-8.580, -3.860,  5.750/), r_H3)), &
!54
         efp_water(water_atom(O1, (/ 3.510, -0.610,  8.530/), r_O1) , &
                   water_atom(H2, (/ 4.010, -1.470,  8.590/), r_H2) , &
                   water_atom(H3, (/ 4.160,  0.160,  8.500/), r_H3)), &
!55
         efp_water(water_atom(O1, (/-0.670, -7.960,  8.730/), r_O1) , &
                   water_atom(H2, (/-1.290, -8.110,  7.970/), r_H2) , &
                   water_atom(H3, (/-1.190, -7.850,  9.580/), r_H3)), &
!56
         efp_water(water_atom(O1, (/-6.350, -3.120, -3.560/), r_O1) , &
                   water_atom(H2, (/-6.290, -3.890, -2.920/), r_H2) , &
                   water_atom(H3, (/-6.870, -3.380, -4.360/), r_H3)), &
!57
         efp_water(water_atom(O1, (/ 3.210, -9.190,  2.420/), r_O1) , &
                   water_atom(H2, (/ 4.030, -8.800,  2.000/), r_H2) , &
                   water_atom(H3, (/ 2.940,-10.010,  1.930/), r_H3)), &
!58
         efp_water(water_atom(O1, (/-4.040,  7.350,  7.280/), r_O1) , &
                   water_atom(H2, (/-4.090,  6.700,  8.030/), r_H2) , &
                   water_atom(H3, (/-3.240,  7.940,  7.410/), r_H3)), &
!59
         efp_water(water_atom(O1, (/ 4.610, -5.960, -1.350/), r_O1) , &
                   water_atom(H2, (/ 4.110, -5.950, -2.210/), r_H2) , &
                   water_atom(H3, (/ 3.980, -6.140, -0.590/), r_H3)), &
!60
         efp_water(water_atom(O1, (/-7.510, -0.860,  2.370/), r_O1) , &
                   water_atom(H2, (/-8.110, -1.480,  2.870/), r_H2) , &
                   water_atom(H3, (/-7.200, -1.300,  1.520/), r_H3))  &
       /

  data (water_array0(i), i=61,80) &
       /                                                     &
!61
         efp_water(water_atom(O1, (/ 2.020,  2.850, -3.640/), r_O1) , &
                   water_atom(H2, (/ 1.220,  3.450, -3.770/), r_H2) , &
                   water_atom(H3, (/ 1.920,  2.360, -2.780/), r_H3)), &
!62
         efp_water(water_atom(O1, (/-2.300, -4.850,  0.810/), r_O1) , &
                   water_atom(H2, (/-2.620, -3.910,  0.710/), r_H2) , &
                   water_atom(H3, (/-3.060, -5.480,  0.690/), r_H3)), &
!63
         efp_water(water_atom(O1, (/ 4.640, -1.190,  3.230/), r_O1) , &
                   water_atom(H2, (/ 4.970, -0.800,  4.090/), r_H2) , &
                   water_atom(H3, (/ 5.400, -1.260,  2.580/), r_H3)), &
!64
         efp_water(water_atom(O1, (/-4.620,  1.070,  4.260/), r_O1) , &
                   water_atom(H2, (/-4.860,  0.700,  3.360/), r_H2) , &
                   water_atom(H3, (/-3.630,  1.230,  4.300/), r_H3)), &
!65
         efp_water(water_atom(O1, (/ 2.490, -0.770, -6.210/), r_O1) , &
                   water_atom(H2, (/ 3.060, -1.420, -5.710/), r_H2) , &
                   water_atom(H3, (/ 2.330, -1.100, -7.140/), r_H3)), &
!66
         efp_water(water_atom(O1, (/-9.220, -1.640,  9.040/), r_O1) , &
                   water_atom(H2, (/-8.420, -2.210,  9.250/), r_H2) , &
                   water_atom(H3, (/-9.710, -2.040,  8.270/), r_H3)), &
!67
         efp_water(water_atom(O1, (/ 3.820,  7.000,  4.800/), r_O1) , &
                   water_atom(H2, (/ 4.270,  6.100,  4.770/), r_H2) , &
                   water_atom(H3, (/ 2.880,  6.890,  5.130/), r_H3)), &
!68
         efp_water(water_atom(O1, (/-3.150,  2.220, -1.330/), r_O1) , &
                   water_atom(H2, (/-3.200,  2.590, -0.410/), r_H2) , &
                   water_atom(H3, (/-3.870,  1.530, -1.450/), r_H3)), &
!69
         efp_water(water_atom(O1, (/ 6.140,  1.220,  1.170/), r_O1) , &
                   water_atom(H2, (/ 7.120,  1.000,  1.240/), r_H2) , &
                   water_atom(H3, (/ 5.830,  1.050,  0.240/), r_H3)), &
!70
         efp_water(water_atom(O1, (/ 7.810,  2.640, -1.130/), r_O1) , &
                   water_atom(H2, (/ 8.480,  2.030, -0.700/), r_H2) , &
                   water_atom(H3, (/ 7.080,  2.830, -0.480/), r_H3)), &
!71
         efp_water(water_atom(O1, (/ 8.880, -3.480, -6.670/), r_O1) , &
                   water_atom(H2, (/ 8.650, -3.730, -7.610/), r_H2) , &
                   water_atom(H3, (/ 9.490, -4.170, -6.280/), r_H3)), &
!72
         efp_water(water_atom(O1, (/-5.110,  5.900, -4.290/), r_O1) , &
                   water_atom(H2, (/-4.830,  5.470, -3.440/), r_H2) , &
                   water_atom(H3, (/-4.860,  6.860, -4.280/), r_H3)), &
!73
         efp_water(water_atom(O1, (/ 8.030, -4.600,  9.240/), r_O1) , &
                   water_atom(H2, (/ 8.930, -4.460,  8.820/), r_H2) , &
                   water_atom(H3, (/ 7.320, -4.580,  8.530/), r_H3)), &
!74
         efp_water(water_atom(O1, (/ 9.220,  5.030,  8.990/), r_O1) , &
                   water_atom(H2, (/ 8.970,  4.940,  8.030/), r_H2) , &
                   water_atom(H3, (/ 9.700,  4.210,  9.300/), r_H3)), &
!75
         efp_water(water_atom(O1, (/ 5.390,  0.640,  5.120/), r_O1) , &
                   water_atom(H2, (/ 4.580,  0.650,  5.700/), r_H2) , &
                   water_atom(H3, (/ 5.420,  1.470,  4.570/), r_H3)), &
!76
         efp_water(water_atom(O1, (/-4.280, -6.740,  0.410/), r_O1) , &
                   water_atom(H2, (/-3.960, -7.500,  0.980/), r_H2) , &
                   water_atom(H3, (/-5.200, -6.470,  0.710/), r_H3)), &
!77
         efp_water(water_atom(O1, (/ 2.970,  0.350,  1.710/), r_O1) , &
                   water_atom(H2, (/ 3.460,  1.190,  1.500/), r_H2) , &
                   water_atom(H3, (/ 3.590, -0.300,  2.160/), r_H3)), &
!78
         efp_water(water_atom(O1, (/-9.270,  2.360,  4.800/), r_O1) , &
                   water_atom(H2, (/-9.750,  2.770,  4.020/), r_H2) , &
                   water_atom(H3, (/-8.280,  2.340,  4.610/), r_H3)), &
!79
         efp_water(water_atom(O1, (/-7.860,  6.830, -3.980/), r_O1) , &
                   water_atom(H2, (/-8.660,  6.220, -3.950/), r_H2) , &
                   water_atom(H3, (/-7.050,  6.300, -4.220/), r_H3)), &
!80
         efp_water(water_atom(O1, (/-6.350, -2.920,  7.930/), r_O1) , &
                   water_atom(H2, (/-6.140, -2.180,  7.280/), r_H2) , &
                   water_atom(H3, (/-5.670, -2.920,  8.660/), r_H3))  &
       /

  data (water_array0(i), i=81,100) &
       /                                                        &
!81
         efp_water(water_atom(O1, (/ 4.590, -7.100,  7.410/), r_O1) , &
                   water_atom(H2, (/ 3.880, -7.370,  8.060/), r_H2) , &
                   water_atom(H3, (/ 4.330, -7.380,  6.480/), r_H3)), &
!82
         efp_water(water_atom(O1, (/-5.910, -0.650,  5.910/), r_O1) , &
                   water_atom(H2, (/-5.470, -0.010,  5.270/), r_H2) , &
                   water_atom(H3, (/-6.410, -0.130,  6.610/), r_H3)), &
!83
         efp_water(water_atom(O1, (/-8.300,  5.490,  0.160/), r_O1) , &
                   water_atom(H2, (/-8.710,  6.310, -0.230/), r_H2) , &
                   water_atom(H3, (/-7.660,  5.750,  0.890/), r_H3)), &
!84
         efp_water(water_atom(O1, (/ 0.780,  5.560, -4.760/), r_O1) , &
                   water_atom(H2, (/ 1.700,  5.550, -5.170/), r_H2) , &
                   water_atom(H3, (/ 0.720,  6.300, -4.090/), r_H3)), &
!85
         efp_water(water_atom(O1, (/ 5.610,  2.220, -7.150/), r_O1) , &
                   water_atom(H2, (/ 5.990,  1.380, -6.780/), r_H2) , &
                   water_atom(H3, (/ 4.730,  2.410, -6.710/), r_H3)), &
!86
         efp_water(water_atom(O1, (/ 8.660,  4.540,  6.420/), r_O1) , &
                   water_atom(H2, (/ 8.340,  5.260,  5.800/), r_H2) , &
                   water_atom(H3, (/ 8.900,  3.730,  5.890/), r_H3)), &
!87
         efp_water(water_atom(O1, (/-8.450,  0.390,  7.530/), r_O1) , &
                   water_atom(H2, (/-9.170,  0.440,  6.840/), r_H2) , &
                   water_atom(H3, (/-8.690, -0.300,  8.220/), r_H3)), &
!88
         efp_water(water_atom(O1, (/-4.330, -6.890,  8.670/), r_O1) , &
                   water_atom(H2, (/-4.880, -7.730,  8.600/), r_H2) , &
                   water_atom(H3, (/-4.070, -6.600,  7.750/), r_H3)), &
!89
         efp_water(water_atom(O1, (/-3.960,  5.900, -8.700/), r_O1) , &
                   water_atom(H2, (/-4.260,  4.950, -8.630/), r_H2) , &
                   water_atom(H3, (/-3.230,  6.060, -8.040/), r_H3)), &
!90
         efp_water(water_atom(O1, (/-0.050,  8.330,  3.770/), r_O1) , &
                   water_atom(H2, (/ 0.370,  7.690,  4.410/), r_H2) , &
                   water_atom(H3, (/-0.430,  7.820,  2.990/), r_H3)), &
!91
         efp_water(water_atom(O1, (/ 4.880, -4.770,  1.740/), r_O1) , &
                   water_atom(H2, (/ 4.010, -4.920,  2.210/), r_H2) , &
                   water_atom(H3, (/ 4.710, -4.510,  0.790/), r_H3)), &
!92
         efp_water(water_atom(O1, (/-1.980, -5.820,  6.570/), r_O1) , &
                   water_atom(H2, (/-0.990, -5.740,  6.710/), r_H2) , &
                   water_atom(H3, (/-2.430, -4.980,  6.880/), r_H3)), &
!93
         efp_water(water_atom(O1, (/-4.720,  5.750,  0.780/), r_O1) , &
                   water_atom(H2, (/-5.260,  5.540,  1.590/), r_H2) , &
                   water_atom(H3, (/-3.810,  5.340,  0.870/), r_H3)), &
!94
         efp_water(water_atom(O1, (/ 5.270,  2.560,  3.280/), r_O1) , &
                   water_atom(H2, (/ 5.540,  1.970,  2.530/), r_H2) , &
                   water_atom(H3, (/ 5.270,  3.510,  2.970/), r_H3)), &
!95
         efp_water(water_atom(O1, (/-1.080, -6.390, -2.740/), r_O1) , &
                   water_atom(H2, (/-0.170, -6.780, -2.870/), r_H2) , &
                   water_atom(H3, (/-1.000, -5.430, -2.500/), r_H3)), &
!96
         efp_water(water_atom(O1, (/-7.980, -5.150, -5.220/), r_O1) , &
                   water_atom(H2, (/-8.780, -5.380, -4.670/), r_H2) , &
                   water_atom(H3, (/-7.150, -5.410, -4.730/), r_H3)), &
!97
         efp_water(water_atom(O1, (/-2.700, -2.330, -2.370/), r_O1) , &
                   water_atom(H2, (/-2.430, -1.990, -3.270/), r_H2) , &
                   water_atom(H3, (/-1.910, -2.710, -1.910/), r_H3)), &
!98
         efp_water(water_atom(O1, (/-7.510, -6.670, -7.620/), r_O1) , &
                   water_atom(H2, (/-7.910, -6.230, -6.810/), r_H2) , &
                   water_atom(H3, (/-7.920, -6.300, -8.450/), r_H3)), &
!99
         efp_water(water_atom(O1, (/-2.240, -7.630, -7.830/), r_O1) , &
                   water_atom(H2, (/-2.190, -6.820, -7.240/), r_H2) , &
                   water_atom(H3, (/-3.100, -7.610, -8.340/), r_H3)), &
!100
         efp_water(water_atom(O1, (/ 9.150,  0.890, -4.600/), r_O1) , &
                   water_atom(H2, (/ 9.400,  0.690, -5.550/), r_H2) , &
                   water_atom(H3, (/ 9.870,  1.450, -4.180/), r_H3))  &
       /

  data (water_array0(i), i=101,120) &
       /                                                     &
!101
         efp_water(water_atom(O1, (/-8.820, -7.460, -1.430/), r_O1) , &
                   water_atom(H2, (/-9.810, -7.400, -1.330/), r_H2) , &
                   water_atom(H3, (/-8.590, -8.260, -1.990/), r_H3)), &
!102
         efp_water(water_atom(O1, (/ 7.050, -8.120,  3.680/), r_O1) , &
                   water_atom(H2, (/ 6.910, -8.050,  4.670/), r_H2) , &
                   water_atom(H3, (/ 7.890, -8.630,  3.500/), r_H3)), &
!103
         efp_water(water_atom(O1, (/ 4.100,  8.130, -6.110/), r_O1) , &
                   water_atom(H2, (/ 4.960,  8.250, -5.610/), r_H2) , &
                   water_atom(H3, (/ 3.680,  7.260, -5.840/), r_H3)), &
!104
         efp_water(water_atom(O1, (/-5.880,  3.860, -6.000/), r_O1) , &
                   water_atom(H2, (/-5.670,  4.600, -5.360/), r_H2) , &
                   water_atom(H3, (/-6.770,  4.030, -6.430/), r_H3)), &
!105
         efp_water(water_atom(O1, (/ 0.640, -2.980, -5.310/), r_O1) , &
                   water_atom(H2, (/ 0.180, -2.160, -5.650/), r_H2) , &
                   water_atom(H3, (/ 1.620, -2.790, -5.220/), r_H3)), &
!106
         efp_water(water_atom(O1, (/ 3.670, -7.620,  5.010/), r_O1) , &
                   water_atom(H2, (/ 3.600, -6.790,  4.450/), r_H2) , &
                   water_atom(H3, (/ 3.710, -8.420,  4.410/), r_H3)), &
!107
         efp_water(water_atom(O1, (/ 5.660,  5.370,  8.650/), r_O1) , &
                   water_atom(H2, (/ 5.780,  6.030,  7.910/), r_H2) , &
                   water_atom(H3, (/ 6.120,  5.710,  9.480/), r_H3)), &
!108
         efp_water(water_atom(O1, (/-6.100, -5.140,  3.880/), r_O1) , &
                   water_atom(H2, (/-5.600, -4.370,  4.280/), r_H2) , &
                   water_atom(H3, (/-7.050, -5.120,  4.200/), r_H3)), &
!109
         efp_water(water_atom(O1, (/-5.900, -4.170, -7.200/), r_O1) , &
                   water_atom(H2, (/-5.430, -4.040, -6.330/), r_H2) , &
                   water_atom(H3, (/-6.560, -4.910, -7.110/), r_H3)), &
!110
         efp_water(water_atom(O1, (/-2.800,  6.390,  4.720/), r_O1) , &
                   water_atom(H2, (/-3.110,  7.000,  5.450/), r_H2) , &
                   water_atom(H3, (/-2.300,  6.910,  4.030/), r_H3)), &
!111
         efp_water(water_atom(O1, (/ 3.540, -3.520, -5.330/), r_O1) , &
                   water_atom(H2, (/ 3.330, -3.960, -6.200/), r_H2) , &
                   water_atom(H3, (/ 4.510, -3.260, -5.300/), r_H3)), &
!112
         efp_water(water_atom(O1, (/ 4.020,  7.510, -2.640/), r_O1) , &
                   water_atom(H2, (/ 4.700,  8.060, -3.110/), r_H2) , &
                   water_atom(H3, (/ 4.420,  6.630, -2.370/), r_H3)), &
!113
         efp_water(water_atom(O1, (/-2.750,  7.790, -1.920/), r_O1) , &
                   water_atom(H2, (/-3.670,  8.170, -1.970/), r_H2) , &
                   water_atom(H3, (/-2.150,  8.260, -2.570/), r_H3)), &
!114
         efp_water(water_atom(O1, (/-8.490,  1.050, -0.920/), r_O1) , &
                   water_atom(H2, (/-8.430,  1.900, -1.440/), r_H2) , &
                   water_atom(H3, (/-8.170,  0.290, -1.490/), r_H3)), &
!115
         efp_water(water_atom(O1, (/ 5.040,  0.500, -1.220/), r_O1) , &
                   water_atom(H2, (/ 4.620, -0.070, -1.920/), r_H2) , &
                   water_atom(H3, (/ 4.380,  1.190, -0.900/), r_H3)), &
!116
         efp_water(water_atom(O1, (/ 5.730,  8.700, -8.330/), r_O1) , &
                   water_atom(H2, (/ 6.170,  9.590, -8.420/), r_H2) , &
                   water_atom(H3, (/ 5.100,  8.700, -7.560/), r_H3)), &
!117
         efp_water(water_atom(O1, (/-5.020,  8.620, -8.170/), r_O1) , &
                   water_atom(H2, (/-5.770,  8.620, -8.830/), r_H2) , &
                   water_atom(H3, (/-4.650,  7.700, -8.080/), r_H3)), &
!118
         efp_water(water_atom(O1, (/-6.530,  5.250,  2.750/), r_O1) , &
                   water_atom(H2, (/-6.400,  4.410,  3.290/), r_H2) , &
                   water_atom(H3, (/-6.820,  5.990,  3.350/), r_H3)), &
!119
         efp_water(water_atom(O1, (/ 3.070,  2.130, -6.310/), r_O1) , &
                   water_atom(H2, (/ 2.840,  2.500, -5.410/), r_H2) , &
                   water_atom(H3, (/ 2.770,  1.180, -6.370/), r_H3)), &
!120
         efp_water(water_atom(O1, (/ 0.370, -5.520, -5.800/), r_O1) , &
                   water_atom(H2, (/ 0.900, -6.010, -5.120/), r_H2) , &
                   water_atom(H3, (/ 0.590, -4.540, -5.750/), r_H3))  &
       /

  data (water_array0(i), i=121,140) &
       /                                                     &
!121
         efp_water(water_atom(O1, (/ 7.320,  6.340, -7.980/), r_O1) , &
                   water_atom(H2, (/ 7.910,  6.080, -8.740/), r_H2) , &
                   water_atom(H3, (/ 7.040,  7.300, -8.090/), r_H3)), &
!122
         efp_water(water_atom(O1, (/-1.340, -9.270, -0.080/), r_O1) , &
                   water_atom(H2, (/-1.800, -9.340, -0.970/), r_H2) , &
                   water_atom(H3, (/-1.960, -8.830,  0.580/), r_H3)), &
!123
         efp_water(water_atom(O1, (/ 3.070,  0.630,  6.180/), r_O1) , &
                   water_atom(H2, (/ 2.960,  1.570,  6.510/), r_H2) , &
                   water_atom(H3, (/ 3.020, -0.000,  6.950/), r_H3)), &
!124
         efp_water(water_atom(O1, (/-2.400,  3.670,  3.740/), r_O1) , &
                   water_atom(H2, (/-2.380,  2.910,  4.380/), r_H2) , &
                   water_atom(H3, (/-2.880,  4.440,  4.140/), r_H3)), &
!125
         efp_water(water_atom(O1, (/-8.390,  7.660, -8.960/), r_O1) , &
                   water_atom(H2, (/-8.240,  7.870, -8.000/), r_H2) , &
                   water_atom(H3, (/-8.690,  6.710, -9.050/), r_H3)), &
!126
         efp_water(water_atom(O1, (/-8.820, -2.890, -1.620/), r_O1) , &
                   water_atom(H2, (/-9.020, -2.450, -2.500/), r_H2) , &
                   water_atom(H3, (/-8.430, -3.800, -1.780/), r_H3)), &
!127
         efp_water(water_atom(O1, (/-0.030, -3.440, -2.570/), r_O1) , &
                   water_atom(H2, (/ 0.110, -3.170, -3.520/), r_H2) , &
                   water_atom(H3, (/ 0.800, -3.220, -2.040/), r_H3)), &
!128
         efp_water(water_atom(O1, (/ 3.500,  8.980, -0.580/), r_O1) , &
                   water_atom(H2, (/ 4.260,  9.420, -0.100/), r_H2) , &
                   water_atom(H3, (/ 3.850,  8.510, -1.400/), r_H3)), &
!129
         efp_water(water_atom(O1, (/-3.220,  2.740,  1.250/), r_O1) , &
                   water_atom(H2, (/-3.830,  1.990,  1.480/), r_H2) , &
                   water_atom(H3, (/-3.000,  3.260,  2.080/), r_H3)), &
!130
         efp_water(water_atom(O1, (/-5.590,  8.380,  0.420/), r_O1) , &
                   water_atom(H2, (/-5.250,  7.450,  0.570/), r_H2) , &
                   water_atom(H3, (/-5.410,  8.650, -0.530/), r_H3)), &
!131
         efp_water(water_atom(O1, (/-7.940, -5.290,  8.490/), r_O1) , &
                   water_atom(H2, (/-7.870, -6.130,  7.940/), r_H2) , &
                   water_atom(H3, (/-7.320, -4.600,  8.130/), r_H3)), &
!132
         efp_water(water_atom(O1, (/ 3.190,  8.100, -9.130/), r_O1) , &
                   water_atom(H2, (/ 4.120,  8.460, -9.080/), r_H2) , &
                   water_atom(H3, (/ 3.130,  7.250, -8.610/), r_H3)), &
!133
         efp_water(water_atom(O1, (/ 3.390,  5.090, -8.560/), r_O1) , &
                   water_atom(H2, (/ 2.870,  4.260, -8.730/), r_H2) , &
                   water_atom(H3, (/ 4.160,  5.140, -9.200/), r_H3)), &
!134
         efp_water(water_atom(O1, (/ 5.110,  4.150, -0.540/), r_O1) , &
                   water_atom(H2, (/ 4.930,  4.600,  0.340/), r_H2) , &
                   water_atom(H3, (/ 5.530,  4.800, -1.170/), r_H3)), &
!135
         efp_water(water_atom(O1, (/-7.240,  3.800, -1.840/), r_O1) , &
                   water_atom(H2, (/-7.690,  4.430, -1.200/), r_H2) , &
                   water_atom(H3, (/-6.310,  4.110, -2.010/), r_H3)), &
!136
         efp_water(water_atom(O1, (/-7.020,  2.070, -3.850/), r_O1) , &
                   water_atom(H2, (/-7.020,  2.710, -3.080/), r_H2) , &
                   water_atom(H3, (/-6.740,  2.550, -4.680/), r_H3)), &
!137
         efp_water(water_atom(O1, (/ 0.080, -5.360,  2.000/), r_O1) , &
                   water_atom(H2, (/-0.850, -5.150,  1.690/), r_H2) , &
                   water_atom(H3, (/ 0.180, -6.350,  2.130/), r_H3)), &
!138
         efp_water(water_atom(O1, (/ 0.880, -0.610,  9.270/), r_O1) , &
                   water_atom(H2, (/ 0.460, -1.470,  9.000/), r_H2) , &
                   water_atom(H3, (/ 1.820, -0.580,  8.930/), r_H3)), &
!139
         efp_water(water_atom(O1, (/ 5.040, -2.940,  9.100/), r_O1) , &
                   water_atom(H2, (/ 5.700, -2.200,  9.190/), r_H2) , &
                   water_atom(H3, (/ 5.480, -3.730,  8.680/), r_H3)), &
!140
         efp_water(water_atom(O1, (/-8.600,  7.960, -6.240/), r_O1) , &
                   water_atom(H2, (/-8.190,  7.640, -5.380/), r_H2) , &
                   water_atom(H3, (/-9.560,  7.690, -6.270/), r_H3))  &
       /

  data (water_array0(i), i=141,160) &
       /                                                     &
!141
         efp_water(water_atom(O1, (/ 0.400,  5.440, -7.480/), r_O1) , &
                   water_atom(H2, (/ 1.250,  5.110, -7.890/), r_H2) , &
                   water_atom(H3, (/ 0.530,  5.590, -6.500/), r_H3)), &
!142
         efp_water(water_atom(O1, (/ 1.890,  5.200, -1.400/), r_O1) , &
                   water_atom(H2, (/ 2.480,  4.800, -2.100/), r_H2) , &
                   water_atom(H3, (/ 1.310,  5.910, -1.810/), r_H3)), &
!143
         efp_water(water_atom(O1, (/-4.930, -9.120, -2.020/), r_O1) , &
                   water_atom(H2, (/-4.540, -8.230, -1.820/), r_H2) , &
                   water_atom(H3, (/-4.830, -9.320, -2.990/), r_H3)), &
!144
         efp_water(water_atom(O1, (/ 8.150,  5.720,  3.250/), r_O1) , &
                   water_atom(H2, (/ 8.220,  4.830,  2.790/), r_H2) , &
                   water_atom(H3, (/ 7.210,  6.060,  3.170/), r_H3)), &
!145
         efp_water(water_atom(O1, (/-2.050,  6.040, -6.560/), r_O1) , &
                   water_atom(H2, (/-2.430,  5.350, -5.940/), r_H2) , &
                   water_atom(H3, (/-1.230,  5.680, -7.000/), r_H3)), &
!146
         efp_water(water_atom(O1, (/ 2.520, -2.980, -1.180/), r_O1) , &
                   water_atom(H2, (/ 2.220, -2.410, -0.420/), r_H2) , &
                   water_atom(H3, (/ 2.450, -3.950, -0.920/), r_H3)), &
!147
         efp_water(water_atom(O1, (/ 6.710,  4.640, -5.930/), r_O1) , &
                   water_atom(H2, (/ 6.370,  3.750, -6.230/), r_H2) , &
                   water_atom(H3, (/ 6.970,  5.180, -6.730/), r_H3)), &
!148
         efp_water(water_atom(O1, (/ 9.300, -1.840, -3.970/), r_O1) , &
                   water_atom(H2, (/ 9.060, -2.020, -4.920/), r_H2) , &
                   water_atom(H3, (/ 9.600, -0.900, -3.870/), r_H3)), &
!149
         efp_water(water_atom(O1, (/ 4.730,  5.000,  1.910/), r_O1) , &
                   water_atom(H2, (/ 5.340,  5.800,  1.950/), r_H2) , &
                   water_atom(H3, (/ 3.780,  5.310,  1.980/), r_H3)), &
!150
         efp_water(water_atom(O1, (/ 1.590, -7.250, -3.960/), r_O1) , &
                   water_atom(H2, (/ 1.810, -7.860, -3.200/), r_H2) , &
                   water_atom(H3, (/ 1.690, -7.740, -4.820/), r_H3)), &
!151
         efp_water(water_atom(O1, (/-5.150, -8.030, -6.280/), r_O1) , &
                   water_atom(H2, (/-4.910, -8.660, -7.020/), r_H2) , &
                   water_atom(H3, (/-6.050, -7.630, -6.460/), r_H3)), &
!152
         efp_water(water_atom(O1, (/-5.600,  8.550,  3.090/), r_O1) , &
                   water_atom(H2, (/-6.460,  8.240,  3.510/), r_H2) , &
                   water_atom(H3, (/-5.640,  8.410,  2.100/), r_H3)), &
!153
         efp_water(water_atom(O1, (/-1.030, -1.150, -7.080/), r_O1) , &
                   water_atom(H2, (/-0.420, -0.850, -7.810/), r_H2) , &
                   water_atom(H3, (/-1.410, -2.040, -7.300/), r_H3)), &
!154
         efp_water(water_atom(O1, (/-6.100, -1.310, -7.340/), r_O1) , &
                   water_atom(H2, (/-5.260, -1.260, -7.880/), r_H2) , &
                   water_atom(H3, (/-6.330, -2.270, -7.160/), r_H3)), &
!155
         efp_water(water_atom(O1, (/ 0.830, -6.040, -8.400/), r_O1) , &
                   water_atom(H2, (/ 0.780, -6.050, -7.400/), r_H2) , &
                   water_atom(H3, (/ 0.000, -6.450, -8.780/), r_H3)), &
!156
         efp_water(water_atom(O1, (/ 6.880, -2.000, -1.460/), r_O1) , &
                   water_atom(H2, (/ 6.320, -1.190, -1.370/), r_H2) , &
                   water_atom(H3, (/ 7.400, -1.960, -2.320/), r_H3)), &
!157
         efp_water(water_atom(O1, (/ 9.030,  0.860,  1.330/), r_O1) , &
                   water_atom(H2, (/ 9.540,  0.870,  0.470/), r_H2) , &
                   water_atom(H3, (/ 9.590,  0.440,  2.040/), r_H3)), &
!158
         efp_water(water_atom(O1, (/-1.360,  1.350,  5.230/), r_O1) , &
                   water_atom(H2, (/-0.630,  1.180,  4.560/), r_H2) , &
                   water_atom(H3, (/-1.670,  0.480,  5.610/), r_H3)), &
!159
         efp_water(water_atom(O1, (/-4.740, -2.890,  4.770/), r_O1) , &
                   water_atom(H2, (/-4.070, -2.770,  4.030/), r_H2) , &
                   water_atom(H3, (/-5.140, -2.000,  5.000/), r_H3)), &
!160
         efp_water(water_atom(O1, (/ 1.300, -0.680, -0.110/), r_O1) , &
                   water_atom(H2, (/ 0.890, -1.420,  0.420/), r_H2) , &
                   water_atom(H3, (/ 1.940, -0.170,  0.470/), r_H3))  &
       /

  data (water_array0(i), i=161,180) &
       /                                                     &
!161
         efp_water(water_atom(O1, (/-5.820,  9.270,  6.720/), r_O1) , &
                   water_atom(H2, (/-5.220,  8.460,  6.740/), r_H2) , &
                   water_atom(H3, (/-5.420,  9.960,  6.120/), r_H3)), &
!162
         efp_water(water_atom(O1, (/ 8.300, -5.890, -4.400/), r_O1) , &
                   water_atom(H2, (/ 8.250, -5.560, -3.450/), r_H2) , &
                   water_atom(H3, (/ 7.440, -5.700, -4.860/), r_H3)), &
!163
         efp_water(water_atom(O1, (/ 6.720, -2.460,  1.540/), r_O1) , &
                   water_atom(H2, (/ 6.810, -2.360,  0.550/), r_H2) , &
                   water_atom(H3, (/ 6.320, -3.350,  1.750/), r_H3)), &
!164
         efp_water(water_atom(O1, (/-2.120, -1.420, -4.680/), r_O1) , &
                   water_atom(H2, (/-1.590, -1.320, -5.520/), r_H2) , &
                   water_atom(H3, (/-2.390, -0.520, -4.340/), r_H3)), &
!165
         efp_water(water_atom(O1, (/-0.210,  1.750, -8.990/), r_O1) , &
                   water_atom(H2, (/ 0.180,  0.900, -9.350/), r_H2) , &
                   water_atom(H3, (/-1.190,  1.770, -9.180/), r_H3)), &
!166
         efp_water(water_atom(O1, (/ 2.630,  3.260,  7.200/), r_O1) , &
                   water_atom(H2, (/ 1.840,  3.770,  6.860/), r_H2) , &
                   water_atom(H3, (/ 2.540,  3.110,  8.180/), r_H3)), &
!167
         efp_water(water_atom(O1, (/-6.680, -2.500,  0.310/), r_O1) , &
                   water_atom(H2, (/-6.620, -3.430,  0.680/), r_H2) , &
                   water_atom(H3, (/-7.270, -2.500, -0.490/), r_H3)), &
!168
         efp_water(water_atom(O1, (/ 8.220, -8.600, -4.900/), r_O1) , &
                   water_atom(H2, (/ 8.620, -8.610, -5.820/), r_H2) , &
                   water_atom(H3, (/ 8.320, -7.680, -4.500/), r_H3)), &
!169
         efp_water(water_atom(O1, (/ 9.160,  9.100,  2.910/), r_O1) , &
                   water_atom(H2, (/ 9.790,  9.480,  2.230/), r_H2) , &
                   water_atom(H3, (/ 9.560,  8.270,  3.300/), r_H3)), &
!170
         efp_water(water_atom(O1, (/-3.580, -2.550,  0.440/), r_O1) , &
                   water_atom(H2, (/-4.500, -2.180,  0.510/), r_H2) , &
                   water_atom(H3, (/-3.200, -2.350, -0.460/), r_H3)), &
!171
         efp_water(water_atom(O1, (/ 3.720, -5.740, -3.720/), r_O1) , &
                   water_atom(H2, (/ 3.590, -4.810, -4.060/), r_H2) , &
                   water_atom(H3, (/ 2.880, -6.260, -3.850/), r_H3)), &
!172
         efp_water(water_atom(O1, (/-2.480, -5.700, -5.730/), r_O1) , &
                   water_atom(H2, (/-1.880, -5.670, -4.930/), r_H2) , &
                   water_atom(H3, (/-3.230, -5.060, -5.600/), r_H3)), &
!173
         efp_water(water_atom(O1, (/-8.230, -7.640,  6.960/), r_O1) , &
                   water_atom(H2, (/-8.930, -8.110,  7.500/), r_H2) , &
                   water_atom(H3, (/-7.640, -8.320,  6.530/), r_H3)), &
!174
         efp_water(water_atom(O1, (/-8.480,  2.360, -8.910/), r_O1) , &
                   water_atom(H2, (/-8.560,  2.000, -9.840/), r_H2) , &
                   water_atom(H3, (/-8.500,  1.600, -8.260/), r_H3)), &
!175
         efp_water(water_atom(O1, (/ 5.900, -3.750,  4.910/), r_O1) , &
                   water_atom(H2, (/ 6.320, -4.330,  4.210/), r_H2) , &
                   water_atom(H3, (/ 5.460, -2.960,  4.470/), r_H3)), &
!176
         efp_water(water_atom(O1, (/-1.530,  3.850, -4.810/), r_O1) , &
                   water_atom(H2, (/-0.800,  4.540, -4.770/), r_H2) , &
                   water_atom(H3, (/-1.250,  3.100, -5.400/), r_H3)), &
!177
         efp_water(water_atom(O1, (/ 2.550, -5.140,  2.900/), r_O1) , &
                   water_atom(H2, (/ 1.590, -5.130,  2.630/), r_H2) , &
                   water_atom(H3, (/ 2.670, -4.610,  3.740/), r_H3)), &
!178
         efp_water(water_atom(O1, (/ 1.050, -8.490, -1.360/), r_O1) , &
                   water_atom(H2, (/ 0.280, -8.820, -0.820/), r_H2) , &
                   water_atom(H3, (/ 1.900, -8.790, -0.940/), r_H3)), &
!179
         efp_water(water_atom(O1, (/ 6.720,  2.030, -3.730/), r_O1) , &
                   water_atom(H2, (/ 7.620,  1.870, -4.130/), r_H2) , &
                   water_atom(H3, (/ 6.800,  2.080, -2.740/), r_H3)), &
!180
         efp_water(water_atom(O1, (/ 0.750,  3.450,  0.330/), r_O1) , &
                   water_atom(H2, (/-0.170,  3.170,  0.040/), r_H2) , &
                   water_atom(H3, (/ 1.060,  4.220, -0.230/), r_H3))  &
       /

  data (water_array0(i), i=181,200) &
       /                                                     &
!181
         efp_water(water_atom(O1, (/-4.220,  8.560, -4.640/), r_O1) , &
                   water_atom(H2, (/-4.790,  9.080, -5.270/), r_H2) , &
                   water_atom(H3, (/-3.260,  8.680, -4.880/), r_H3)), &
!182
         efp_water(water_atom(O1, (/ 0.720,  1.660,  3.180/), r_O1) , &
                   water_atom(H2, (/ 0.550,  2.490,  2.640/), r_H2) , &
                   water_atom(H3, (/ 1.620,  1.290,  2.960/), r_H3)), &
!183
         efp_water(water_atom(O1, (/-6.790, -5.270,  1.190/), r_O1) , &
                   water_atom(H2, (/-7.780, -5.380,  1.210/), r_H2) , &
                   water_atom(H3, (/-6.450, -5.120,  2.120/), r_H3)), &
!184
         efp_water(water_atom(O1, (/ 6.130,  8.420, -4.310/), r_O1) , &
                   water_atom(H2, (/ 6.690,  9.230, -4.480/), r_H2) , &
                   water_atom(H3, (/ 6.720,  7.620, -4.280/), r_H3)), &
!185
         efp_water(water_atom(O1, (/-3.690, -0.950, -9.030/), r_O1) , &
                   water_atom(H2, (/-3.360, -0.310, -9.720/), r_H2) , &
                   water_atom(H3, (/-3.030, -1.010, -8.280/), r_H3)), &
!186
         efp_water(water_atom(O1, (/ 7.160,  5.650, -1.540/), r_O1) , &
                   water_atom(H2, (/ 7.350,  6.300, -0.800/), r_H2) , &
                   water_atom(H3, (/ 7.760,  4.850, -1.450/), r_H3)), &
!187
         efp_water(water_atom(O1, (/-4.120, -6.420, -2.290/), r_O1) , &
                   water_atom(H2, (/-4.210, -6.520, -1.300/), r_H2) , &
                   water_atom(H3, (/-3.160, -6.490, -2.550/), r_H3)), &
!188
         efp_water(water_atom(O1, (/ 3.900, -1.210, -3.020/), r_O1) , &
                   water_atom(H2, (/ 2.990, -0.800, -3.040/), r_H2) , &
                   water_atom(H3, (/ 3.830, -2.150, -2.700/), r_H3)), &
!189
         efp_water(water_atom(O1, (/-1.880,  8.830, -6.080/), r_O1) , &
                   water_atom(H2, (/-2.150,  7.940, -6.450/), r_H2) , &
                   water_atom(H3, (/-1.870,  9.510, -6.810/), r_H3)), &
!190
         efp_water(water_atom(O1, (/-6.370,  3.250,  4.490/), r_O1) , &
                   water_atom(H2, (/-5.720,  2.510,  4.380/), r_H2) , &
                   water_atom(H3, (/-6.170,  3.750,  5.330/), r_H3)), &
!191
         efp_water(water_atom(O1, (/ 5.940,  7.450,  6.520/), r_O1) , &
                   water_atom(H2, (/ 6.440,  8.300,  6.330/), r_H2) , &
                   water_atom(H3, (/ 5.060,  7.470,  6.040/), r_H3)), &
!192
         efp_water(water_atom(O1, (/-0.850,  3.420, -2.200/), r_O1) , &
                   water_atom(H2, (/-1.020,  3.730, -3.140/), r_H2) , &
                   water_atom(H3, (/-1.690,  3.050, -1.820/), r_H3)), &
!193
         efp_water(water_atom(O1, (/-1.320, -9.280, -3.450/), r_O1) , &
                   water_atom(H2, (/-0.940, -8.370, -3.300/), r_H2) , &
                   water_atom(H3, (/-1.400, -9.450, -4.440/), r_H3)), &
!194
         efp_water(water_atom(O1, (/ 8.590, -4.880,  0.160/), r_O1) , &
                   water_atom(H2, (/ 8.130, -4.730,  1.040/), r_H2) , &
                   water_atom(H3, (/ 9.030, -4.030, -0.140/), r_H3)), &
!195
         efp_water(water_atom(O1, (/ 6.610, -0.720, -9.090/), r_O1) , &
                   water_atom(H2, (/ 6.150,  0.160, -9.220/), r_H2) , &
                   water_atom(H3, (/ 7.600, -0.600, -9.160/), r_H3)), &
!196
         efp_water(water_atom(O1, (/-4.540, -0.110, -1.420/), r_O1) , &
                   water_atom(H2, (/-5.500, -0.220, -1.690/), r_H2) , &
                   water_atom(H3, (/-3.980, -0.780, -1.900/), r_H3)), &
!197
         efp_water(water_atom(O1, (/ 8.590, -9.060,  8.610/), r_O1) , &
                   water_atom(H2, (/ 9.130, -9.750,  9.090/), r_H2) , &
                   water_atom(H3, (/ 8.270, -8.370,  9.270/), r_H3)), &
!198
         efp_water(water_atom(O1, (/-7.790, -8.780,  0.870/), r_O1) , &
                   water_atom(H2, (/-8.020, -8.250,  0.050/), r_H2) , &
                   water_atom(H3, (/-6.980, -9.340,  0.680/), r_H3)), &
!199
         efp_water(water_atom(O1, (/-0.010, -2.930,  8.510/), r_O1) , &
                   water_atom(H2, (/-0.720, -3.050,  7.810/), r_H2) , &
                   water_atom(H3, (/ 0.000, -3.720,  9.110/), r_H3)), &
!200
         efp_water(water_atom(O1, (/ 2.210, -5.480, -0.180/), r_O1) , &
                   water_atom(H2, (/ 1.560, -6.210, -0.390/), r_H2) , &
                   water_atom(H3, (/ 2.250, -5.340,  0.800/), r_H3))  &
       /

  data (water_array0(i), i=201,216) &
       /                                                     &
!201
         efp_water(water_atom(O1, (/ 0.790, -6.220,  6.530/), r_O1) , &
                   water_atom(H2, (/ 0.780, -6.690,  7.410/), r_H2) , &
                   water_atom(H3, (/ 1.610, -6.500,  6.020/), r_H3)), &
!202
         efp_water(water_atom(O1, (/ 6.720, -4.710, -2.380/), r_O1) , &
                   water_atom(H2, (/ 5.940, -5.210, -2.000/), r_H2) , &
                   water_atom(H3, (/ 6.690, -3.760, -2.070/), r_H3)), &
!203
         efp_water(water_atom(O1, (/-0.380,  1.920, -6.350/), r_O1) , &
                   water_atom(H2, (/-0.420,  1.020, -5.910/), r_H2) , &
                   water_atom(H3, (/-0.350,  1.810, -7.340/), r_H3)), &
!204
         efp_water(water_atom(O1, (/ 4.280,  4.240,  5.200/), r_O1) , &
                   water_atom(H2, (/ 4.580,  3.520,  4.580/), r_H2) , &
                   water_atom(H3, (/ 3.890,  3.840,  6.030/), r_H3)), &
!205
         efp_water(water_atom(O1, (/-1.570, -3.750, -7.580/), r_O1) , &
                   water_atom(H2, (/-2.500, -4.000, -7.850/), r_H2) , &
                   water_atom(H3, (/-1.310, -4.250, -6.760/), r_H3)), &
!206
         efp_water(water_atom(O1, (/ 3.170,  5.470, -5.820/), r_O1) , &
                   water_atom(H2, (/ 3.550,  4.880, -5.100/), r_H2) , &
                   water_atom(H3, (/ 3.570,  5.210, -6.700/), r_H3)), &
!207
         efp_water(water_atom(O1, (/ 8.120, -2.760,  6.870/), r_O1) , &
                   water_atom(H2, (/ 8.440, -2.660,  5.930/), r_H2) , &
                   water_atom(H3, (/ 7.330, -3.380,  6.890/), r_H3)), &
!208
         efp_water(water_atom(O1, (/-4.380,  2.140, -7.500/), r_O1) , &
                   water_atom(H2, (/-3.860,  1.490, -6.950/), r_H2) , &
                   water_atom(H3, (/-4.870,  2.770, -6.890/), r_H3)), &
!209
         efp_water(water_atom(O1, (/-8.610,  0.340, -7.080/), r_O1) , &
                   water_atom(H2, (/-9.240, -0.380, -7.390/), r_H2) , &
                   water_atom(H3, (/-7.680, -0.020, -7.080/), r_H3)), &
!210
         efp_water(water_atom(O1, (/ 7.700, -5.320,  3.010/), r_O1) , &
                   water_atom(H2, (/ 7.240, -6.190,  3.180/), r_H2) , &
                   water_atom(H3, (/ 8.610, -5.350,  3.420/), r_H3)), &
!211
         efp_water(water_atom(O1, (/ 6.180, -2.950, -5.780/), r_O1) , &
                   water_atom(H2, (/ 6.130, -2.130, -5.210/), r_H2) , &
                   water_atom(H3, (/ 7.070, -2.980, -6.230/), r_H3)), &
!212
         efp_water(water_atom(O1, (/-5.100,  0.520,  1.680/), r_O1) , &
                   water_atom(H2, (/-4.750,  0.110,  0.840/), r_H2) , &
                   water_atom(H3, (/-6.000,  0.140,  1.880/), r_H3)), &
!213
         efp_water(water_atom(O1, (/-5.620,  4.530,  6.910/), r_O1) , &
                   water_atom(H2, (/-6.210,  5.330,  6.950/), r_H2) , &
                   water_atom(H3, (/-5.470,  4.180,  7.840/), r_H3)), &
!214
         efp_water(water_atom(O1, (/-2.690,  2.210,  8.820/), r_O1) , &
                   water_atom(H2, (/-3.530,  2.200,  9.360/), r_H2) , &
                   water_atom(H3, (/-2.670,  3.040,  8.260/), r_H3)), &
!215
         efp_water(water_atom(O1, (/ 0.390, -7.850,  3.000/), r_O1) , &
                   water_atom(H2, (/ 1.380, -7.960,  2.910/), r_H2) , &
                   water_atom(H3, (/-0.010, -8.710,  3.320/), r_H3)), &
!216
         efp_water(water_atom(O1, (/ 8.750, -2.160,  3.370/), r_O1) , &
                   water_atom(H2, (/ 7.980, -2.510,  2.830/), r_H2) , &
                   water_atom(H3, (/ 8.430, -1.450,  3.990/), r_H3))  &
       /

  data (water_array1(i), i=1,20) &
       /                                                    &
!1
         efp_water(water_atom(O1, (/ 1.47101773,     -0.69414993,     -0.56263210/), r_O1) , &
                   water_atom(H2, (/ 0.84679968,     -1.33025089,     -0.25182949/), r_H2) , &
                   water_atom(H3, (/ 1.89626675,     -0.33618544,      0.20020205/), r_H3)),  &
!2
         efp_water(water_atom(O1, (/ 0.74994304,      1.15648567,     -2.26062771/), r_O1) , &
                   water_atom(H2, (/ 0.84484740,      0.58546563,     -1.51509902/), r_H2) , &
                   water_atom(H3, (/-0.15769031,      1.11938601,     -2.51698932/), r_H3)), &
!3
         efp_water(water_atom(O1, (/-0.11572351,     -2.52078943,      0.52306203/), r_O1) , &
                   water_atom(H2, (/ 0.17520241,     -3.41691479,      0.57963406/), r_H2) , &
                   water_atom(H3, (/-0.08388190,     -2.16945887,      1.39852356/), r_H3)), &
!4
         efp_water(water_atom(O1, (/ 0.63821689,      1.40516265,      1.90421696/), r_O1) , &
                   water_atom(H2, (/ 0.47520344,      1.95978020,      1.15809011/), r_H2) , &
                   water_atom(H3, (/ 1.54318325,      1.14176074,      1.85388467/), r_H3)), &
!5
         efp_water(water_atom(O1, (/ 0.82672388,      3.06852313,     -0.09378427/), r_O1) , &
                   water_atom(H2, (/ 0.20133278,      3.24546395,     -0.77822113/), r_H2) , &
                   water_atom(H3, (/ 0.96568214,      3.88230999,      0.36373028/), r_H3)), &
!6
         efp_water(water_atom(O1, (/ 2.96274758,      0.47416172,      1.23742600/), r_O1) , &
                   water_atom(H2, (/ 3.43105695,      1.10984965,      0.72027581/), r_H2) , &
                   water_atom(H3, (/ 3.61235604,     -0.03199250,      1.69864372/), r_H3)), &
!7
         efp_water(water_atom(O1, (/-0.04887542,     -1.03828236,      2.63131815/), r_O1) , &
                   water_atom(H2, (/ 0.67865262,     -1.05991908,      3.23224784/), r_H2) , &
                   water_atom(H3, (/-0.07271689,     -0.16856077,      2.26540232/), r_H3)), &
!8
         efp_water(water_atom(O1, (/ 3.45393836,     -1.78109489,     -1.81542483/), r_O1) , &
                   water_atom(H2, (/ 2.69342220,     -1.64205879,     -1.27409581/), r_H2) , &
                   water_atom(H3, (/ 3.72393571,     -2.67524298,     -1.67944332/), r_H3)), &
!9
         efp_water(water_atom(O1, (/-2.88382767,      2.35063734,     -1.33623884/), r_O1) , &
                   water_atom(H2, (/-2.84744073,      2.60987326,     -0.42940381/), r_H2) , &
                   water_atom(H3, (/-3.05693103,      1.42283884,     -1.34641852/), r_H3)), &
!10
         efp_water(water_atom(O1, (/-1.85093038,     -2.03478432,     -1.73336629/), r_O1) , &
                   water_atom(H2, (/-1.80228164,     -2.07663120,     -2.67504566/), r_H2) , &
                   water_atom(H3, (/-1.01106766,     -2.31426554,     -1.40564271/), r_H3)), &
!11
         efp_water(water_atom(O1, (/ 3.37735630,      2.35960091,     -0.50986456/), r_O1) , &
                   water_atom(H2, (/ 3.83762462,      3.18313583,     -0.53854379/), r_H2) , &
                   water_atom(H3, (/ 2.45581814,      2.56366785,     -0.51147138/), r_H3)), &
!12
         efp_water(water_atom(O1, (/ 0.50536512,     -3.31803899,     -2.10176543/), r_O1) , &
                   water_atom(H2, (/ 0.99525557,     -3.02362530,     -2.85290294/), r_H2) , &
                   water_atom(H3, (/ 1.06916879,     -3.90435542,     -1.62297913/), r_H3)), &
!13
         efp_water(water_atom(O1, (/-2.55931110,     -1.86153057,      2.60685058/), r_O1) , &
                   water_atom(H2, (/-1.62939680,     -1.75214630,      2.72590052/), r_H2) , &
                   water_atom(H3, (/-2.68785778,     -2.13703959,      1.71329080/), r_H3)), &
!14
         efp_water(water_atom(O1, (/-0.98910493,      4.05379010,     -1.80275376/), r_O1) , &
                   water_atom(H2, (/-0.76667488,      4.32515256,     -2.67897572/), r_H2) , &
                   water_atom(H3, (/-1.67378656,      3.40920178,     -1.88396573/), r_H3)), &
!15
         efp_water(water_atom(O1, (/-3.33052687,     -2.79184834,      0.33736957/), r_O1) , &
                   water_atom(H2, (/-4.26718700,     -2.70292929,      0.41246623/), r_H2) , &
                   water_atom(H3, (/-3.09623618,     -2.47815123,     -0.52145550/), r_H3)), &
!16
         efp_water(water_atom(O1, (/-3.16293216,      2.74725331,      1.17684187/), r_O1) , &
                   water_atom(H2, (/-3.51387337,      1.92543538,      1.48070751/), r_H2) , &
                   water_atom(H3, (/-2.66095093,      3.10816010,      1.89003102/), r_H3)), &
!17
         efp_water(water_atom(O1, (/ 2.15124559,      3.41712877,     -3.08788931/), r_O1) , &
                   water_atom(H2, (/ 1.96861874,      4.17642920,     -2.55780493/), r_H2) , &
                   water_atom(H3, (/ 1.74514373,      2.68177963,     -2.65751418/), r_H3)), &
!18
         efp_water(water_atom(O1, (/ 2.46277593,     -0.21469374,     -3.73654774/), r_O1) , &
                   water_atom(H2, (/ 1.87883028,      0.34613547,     -3.25134302/), r_H2) , &
                   water_atom(H3, (/ 2.55354936,     -1.00973270,     -3.23600606/), r_H3)), &
!19
         efp_water(water_atom(O1, (/-2.70915253,      1.19069385,     -4.20054451/), r_O1) , &
                   water_atom(H2, (/-3.40945789,      1.42960132,     -4.78652584/), r_H2) , &
                   water_atom(H3, (/-2.87576846,      1.63361187,     -3.38387957/), r_H3)), &
!20
         efp_water(water_atom(O1, (/ 4.92929796,      0.25926647,     -1.08037432/), r_O1) , &
                   water_atom(H2, (/ 4.44326131,     -0.51348823,     -1.32015423/), r_H2) , &
                   water_atom(H3, (/ 4.35283180,      0.99461887,     -1.21386847/), r_H3))  &
       /

  data (water_array1(i), i=21,40) &
       /                                                     &
!21
         efp_water(water_atom(O1, (/-3.51239673,     -0.12037135,     -1.31785877/), r_O1) , &
                   water_atom(H2, (/-4.15226469,     -0.32368495,     -1.98126155/), r_H2) , &
                   water_atom(H3, (/-2.80389035,     -0.73434568,     -1.42710286/), r_H3)), &
!22
         efp_water(water_atom(O1, (/-2.24959503,     -1.46068111,     -4.24155719/), r_O1) , &
                   water_atom(H2, (/-1.72182705,     -1.53767826,     -5.02028210/), r_H2) , &
                   water_atom(H3, (/-2.50007454,     -0.55330511,     -4.17222385/), r_H3)), &
!23
         efp_water(water_atom(O1, (/-0.36458708,      2.17343082,      4.16515428/), r_O1) , &
                   water_atom(H2, (/ 0.10898907,      2.00231167,      3.36683004/), r_H2) , &
                   water_atom(H3, (/-0.55635605,      1.33412922,      4.55203379/), r_H3)), &
!24
         efp_water(water_atom(O1, (/-4.43108531,      0.38939119,      1.18791038/), r_O1) , &
                   water_atom(H2, (/-4.09220225,      0.07154898,      0.36634696/), r_H2) , &
                   water_atom(H3, (/-5.35576524,      0.53174018,      1.06301039/), r_H3)), &
!25
         efp_water(water_atom(O1, (/-2.00124531,     -5.02233355,      0.47294195/), r_O1) , &
                   water_atom(H2, (/-2.37030808,     -4.15908697,      0.57029067/), r_H2) , &
                   water_atom(H3, (/-2.72015620,     -5.63190560,      0.52261326/), r_H3)), &
!26
         efp_water(water_atom(O1, (/ 2.17498240,     -0.77753005,      3.99551492/), r_O1) , &
                   water_atom(H2, (/ 2.37397522,     -0.09742355,      4.61900405/), r_H2) , &
                   water_atom(H3, (/ 2.94536094,     -0.89028780,      3.46197058/), r_H3)), &
!27
         efp_water(water_atom(O1, (/ 2.32803577,      5.60203107,     -1.48703018/), r_O1) , &
                   water_atom(H2, (/ 3.17756435,      5.39229525,     -1.13322904/), r_H2) , &
                   water_atom(H3, (/ 2.44909957,      6.33849249,     -2.06482777/), r_H3)), &
!28
         efp_water(water_atom(O1, (/ 0.29409048,     -5.63787020,      1.35390380/), r_O1) , &
                   water_atom(H2, (/-0.61484722,     -5.46247366,      1.16965661/), r_H2) , &
                   water_atom(H3, (/ 0.35500896,     -6.55453054,      1.57047150/), r_H3)), &
!29
         efp_water(water_atom(O1, (/-3.13820696,      0.43930747,      3.68520638/), r_O1) , &
                   water_atom(H2, (/-3.17474480,     -0.15120527,      2.94979113/), r_H2) , &
                   water_atom(H3, (/-2.32205829,      0.27061032,      4.12824552/), r_H3)), &
!30
         efp_water(water_atom(O1, (/ 1.79688365,     -5.25124706,     -0.75847014/), r_O1) , &
                   water_atom(H2, (/ 1.67800276,     -6.05488666,     -1.23899828/), r_H2) , &
                   water_atom(H3, (/ 1.35537136,     -5.36174751,      0.06841934/), r_H3)), &
!31
         efp_water(water_atom(O1, (/-2.04325924,      3.99071634,      3.22018842/), r_O1) , &
                   water_atom(H2, (/-1.57679177,      3.38300019,      3.77152087/), r_H2) , &
                   water_atom(H3, (/-2.74724170,      4.34369942,      3.74046796/), r_H3)), &
!32
         efp_water(water_atom(O1, (/ 1.70735837,     -2.94873386,     -4.37735891/), r_O1) , &
                   water_atom(H2, (/ 1.47312003,     -2.31135262,     -5.03291784/), r_H2) , &
                   water_atom(H3, (/ 2.61330042,     -3.16891746,     -4.52455137/), r_H3)), &
!33
         efp_water(water_atom(O1, (/ 4.47995634,     -0.95643849,      2.72704787/), r_O1) , &
                   water_atom(H2, (/ 5.02538777,     -0.45296412,      3.31008172/), r_H2) , &
                   water_atom(H3, (/ 4.99136795,     -1.70053684,      2.45202081/), r_H3)), &
!34
         efp_water(water_atom(O1, (/-2.05218206,      5.99634932,     -0.27609371/), r_O1) , &
                   water_atom(H2, (/-2.21116849,      6.73897822,     -0.83654380/), r_H2) , &
                   water_atom(H3, (/-1.62918371,      5.34100941,     -0.80758285/), r_H3)), &
!35
         efp_water(water_atom(O1, (/-0.90605207,      6.05297197,      2.06044528/), r_O1) , &
                   water_atom(H2, (/-1.44942776,      6.21822631,      1.30656921/), r_H2) , &
                   water_atom(H3, (/-1.26984478,      5.29942015,      2.49712151/), r_H3)), &
!36
         efp_water(water_atom(O1, (/ 2.05098399,     -4.70320905,      2.94050003/), r_O1) , &
                   water_atom(H2, (/ 1.26320392,     -4.80783119,      2.43124334/), r_H2) , &
                   water_atom(H3, (/ 1.91199179,     -3.96013504,      3.50565418/), r_H3)), &
!37
         efp_water(water_atom(O1, (/-1.20227435,     -0.01774219,      5.34831591/), r_O1) , &
                   water_atom(H2, (/-0.81897359,     -0.87497305,      5.44378650/), r_H2) , &
                   water_atom(H3, (/-1.47869630,      0.25322614,      6.20915500/), r_H3)), &
!38
         efp_water(water_atom(O1, (/-0.47036549,      4.62778872,     -4.30435894/), r_O1) , &
                   water_atom(H2, (/ 0.39922028,      4.95919822,     -4.46206050/), r_H2) , &
                   water_atom(H3, (/-0.52994516,      3.79872194,     -4.75154615/), r_H3)), &
!39
         efp_water(water_atom(O1, (/ 6.15209884,      1.29493497,      0.98077213/), r_O1) , &
                   water_atom(H2, (/ 6.93981532,      0.85449213,      1.25717017/), r_H2) , &
                   water_atom(H3, (/ 5.68360022,      0.69093686,      0.42707609/), r_H3)), &
!40
         efp_water(water_atom(O1, (/ 1.31109970,      5.38740003,      0.94654924/), r_O1) , &
                   water_atom(H2, (/ 0.54758028,      5.77982225,      1.33887760/), r_H2) , &
                   water_atom(H3, (/ 1.46241727,      5.84221926,      0.13345597/), r_H3))  &
       /

  data (water_array1(i), i=41,60) &
       /                                                     &
!41
         efp_water(water_atom(O1, (/-4.93696531,      3.95792133,     -1.77469249/), r_O1) , &
                   water_atom(H2, (/-4.86695849,      4.44698118,     -0.97045524/), r_H2) , &
                   water_atom(H3, (/-4.22760276,      3.33528080,     -1.77619108/), r_H3)), &
!42
         efp_water(water_atom(O1, (/ 5.06979708,      2.86012896,      2.80269934/), r_O1) , &
                   water_atom(H2, (/ 5.38210460,      2.28441701,      2.12306926/), r_H2) , &
                   water_atom(H3, (/ 4.29156725,      3.27693118,      2.46876123/), r_H3)), &
!43
         efp_water(water_atom(O1, (/ 4.63010411,      3.88508614,     -3.89332475/), r_O1) , &
                   water_atom(H2, (/ 5.28900290,      3.23650139,     -3.70340704/), r_H2) , &
                   water_atom(H3, (/ 3.83170770,      3.58772976,     -3.48709723/), r_H3)), &
!44
         efp_water(water_atom(O1, (/-0.69850775,      2.33467612,     -5.56300298/), r_O1) , &
                   water_atom(H2, (/-1.38686901,      1.80113354,     -5.19915476/), r_H2) , &
                   water_atom(H3, (/-0.97680339,      2.58224323,     -6.43025482/), r_H3)), &
!45
         efp_water(water_atom(O1, (/ 4.58457747,     -4.94532858,      2.62963055/), r_O1) , &
                   water_atom(H2, (/ 3.65959954,     -4.84102639,      2.47338091/), r_H2) , &
                   water_atom(H3, (/ 4.77076068,     -5.86727879,      2.55075277/), r_H3)), &
!46
         efp_water(water_atom(O1, (/-1.42183469,     -5.26605689,     -2.12713028/), r_O1) , &
                   water_atom(H2, (/-1.48625514,     -5.20310526,     -1.18757516/), r_H2) , &
                   water_atom(H3, (/-0.84587268,     -4.57249787,     -2.40662437/), r_H3)), &
!47
         efp_water(water_atom(O1, (/ 4.61424811,      4.85923816,     -0.19648524/), r_O1) , &
                   water_atom(H2, (/ 4.37530241,      4.93099868,      0.71380635/), r_H2) , &
                   water_atom(H3, (/ 5.17142304,      5.59621585,     -0.38960914/), r_H3)), &
!48
         efp_water(water_atom(O1, (/ 1.50470511,     -0.88446803,     -6.09222161/), r_O1) , &
                   water_atom(H2, (/ 1.83438425,     -0.40409025,     -5.34963813/), r_H2) , &
                   water_atom(H3, (/ 1.99504578,     -0.59209847,     -6.84386026/), r_H3)), &
!49
         efp_water(water_atom(O1, (/ 1.90132486,      2.50689632,     -5.59571926/), r_O1) , &
                   water_atom(H2, (/ 2.05573665,      3.07381225,     -4.85704691/), r_H2) , &
                   water_atom(H3, (/ 1.04391938,      2.13124076,     -5.47479853/), r_H3)), &
!50
         efp_water(water_atom(O1, (/-4.19246082,     -2.86288013,      4.37510206/), r_O1) , &
                   water_atom(H2, (/-3.58654229,     -2.54507355,      3.72491779/), r_H2) , &
                   water_atom(H3, (/-4.64359324,     -2.10784182,      4.71755712/), r_H3)), &
!51
         efp_water(water_atom(O1, (/-3.87681951,     -3.48471901,     -4.38026627/), r_O1) , &
                   water_atom(H2, (/-3.39707927,     -2.67496845,     -4.30932835/), r_H2) , &
                   water_atom(H3, (/-4.55805190,     -3.45450182,     -3.72766492/), r_H3)), &
!52
         efp_water(water_atom(O1, (/ 2.30737613,     -3.09164249,      4.88935946/), r_O1) , &
                   water_atom(H2, (/ 2.08773378,     -2.17520572,      4.83664340/), r_H2) , &
                   water_atom(H3, (/ 1.91792814,     -3.41863872,      5.68452577/), r_H3)), &
!53
         efp_water(water_atom(O1, (/ 3.09104598,      0.62490818,      5.92533115/), r_O1) , &
                   water_atom(H2, (/ 2.71227564,      1.44478864,      6.19957145/), r_H2) , &
                   water_atom(H3, (/ 3.09683118,      0.06189120,      6.68286477/), r_H3)), &
!54
         efp_water(water_atom(O1, (/ 4.93153968,      0.29248547,     -4.28759253/), r_O1) , &
                   water_atom(H2, (/ 4.10311104,      0.11063181,     -3.87345442/), r_H2) , &
                   water_atom(H3, (/ 5.42629236,      0.82262388,     -3.68339464/), r_H3)), &
!55
         efp_water(water_atom(O1, (/-5.56550258,     -0.47719412,     -2.92730689/), r_O1) , &
                   water_atom(H2, (/-5.91198686,      0.10951330,     -3.58045378/), r_H2) , &
                   water_atom(H3, (/-5.86510924,     -1.34303274,     -3.15410510/), r_H3)), &
!56
         efp_water(water_atom(O1, (/ 6.01585958,     -2.87055882,      1.95974261/), r_O1) , &
                   water_atom(H2, (/ 6.17367620,     -2.91053316,      1.03001781/), r_H2) , &
                   water_atom(H3, (/ 5.69561826,     -3.72024991,      2.21733925/), r_H3)), &
!57
         efp_water(water_atom(O1, (/ 7.38902371,     -0.80122214,     -1.43455118/), r_O1) , &
                   water_atom(H2, (/ 6.63304244,     -0.23911478,     -1.49302090/), r_H2) , &
                   water_atom(H3, (/ 7.85571311,     -0.71374812,     -2.25020718/), r_H3)), &
!58
         efp_water(water_atom(O1, (/-4.53870313,      5.37601181,      0.33812165/), r_O1) , &
                   water_atom(H2, (/-4.71761970,      4.86570793,      1.11172208/), r_H2) , &
                   water_atom(H3, (/-3.61223542,      5.55637525,      0.34038064/), r_H3)), &
!59
         efp_water(water_atom(O1, (/ 5.33650731,      0.73151705,      4.52993612/), r_O1) , &
                   water_atom(H2, (/ 4.60412710,      0.72640570,      5.12530847/), r_H2) , &
                   water_atom(H3, (/ 5.33423380,      1.57396010,      4.10430331/), r_H3)), &
!60
         efp_water(water_atom(O1, (/ 0.05023806,      3.72206291,      6.18449264/), r_O1) , &
                   water_atom(H2, (/-0.79374203,      3.62081747,      6.59476983/), r_H2) , &
                   water_atom(H3, (/ 0.02064400,      3.23174797,      5.37851897/), r_H3))  &
       /

  data (water_array1(i), i=61,80) &
       /                                                     &
!61
         efp_water(water_atom(O1, (/ 3.51925838,      4.85367577,      2.16121899/), r_O1) , &
                   water_atom(H2, (/ 3.76107965,      5.57395669,      2.72121397/), r_H2) , &
                   water_atom(H3, (/ 2.63187984,      5.01322035,      1.88197054/), r_H3)), &
!62
         efp_water(water_atom(O1, (/-5.80562515,     -2.69819555,      1.19588443/), r_O1) , &
                   water_atom(H2, (/-5.97629604,     -3.54857752,      1.56817856/), r_H2) , &
                   water_atom(H3, (/-6.57013068,     -2.46303911,      0.69479205/), r_H3)), &
!63
         efp_water(water_atom(O1, (/ 4.18612251,     -4.31065864,     -1.34437501/), r_O1) , &
                   water_atom(H2, (/ 4.23779039,     -4.76251298,     -2.17143866/), r_H2) , &
                   water_atom(H3, (/ 3.41337660,     -4.63166150,     -0.90770434/), r_H3)), &
!64
         efp_water(water_atom(O1, (/ 1.69511790,      5.83904745,     -5.08926647/), r_O1) , &
                   water_atom(H2, (/ 2.48314563,      5.44533531,     -5.42820857/), r_H2) , &
                   water_atom(H3, (/ 1.96512105,      6.56207727,     -4.54593948/), r_H3)), &
!65
         efp_water(water_atom(O1, (/-0.98071314,     -1.52633823,     -6.50025356/), r_O1) , &
                   water_atom(H2, (/-0.08697974,     -1.22808230,     -6.44401697/), r_H2) , &
                   water_atom(H3, (/-0.96025386,     -2.36434625,     -6.93407680/), r_H3)), &
!66
         efp_water(water_atom(O1, (/ 4.39115841,     -5.48954730,     -3.59288539/), r_O1) , &
                   water_atom(H2, (/ 4.25195274,     -4.89102535,     -4.30931624/), r_H2) , &
                   water_atom(H3, (/ 3.77645243,     -6.19668125,     -3.70680443/), r_H3)), &
!67
         efp_water(water_atom(O1, (/ 4.03404895,     -3.57849443,     -5.25934589/), r_O1) , &
                   water_atom(H2, (/ 3.80587578,     -3.65984329,     -6.17159497/), r_H2) , &
                   water_atom(H3, (/ 4.73661419,     -2.95009342,     -5.21042742/), r_H3)), &
!68
         efp_water(water_atom(O1, (/ 0.46885488,     -5.18710493,     -4.69707740/), r_O1) , &
                   water_atom(H2, (/ 0.98546095,     -5.84333254,     -4.25734658/), r_H2) , &
                   water_atom(H3, (/ 0.81218518,     -4.34887474,     -4.43179587/), r_H3)), &
!69
         efp_water(water_atom(O1, (/ 4.39159292,      3.81336069,      5.59649846/), r_O1) , &
                   water_atom(H2, (/ 4.47772473,      3.40939405,      4.74781367/), r_H2) , &
                   water_atom(H3, (/ 3.74018812,      3.31871701,      6.06751332/), r_H3)), &
!70
         efp_water(water_atom(O1, (/-3.79032107,     -6.41514673,     -2.31579643/), r_O1) , &
                   water_atom(H2, (/-3.95127266,     -6.39707475,     -1.38593096/), r_H2) , &
                   water_atom(H3, (/-3.08626329,     -5.80863143,     -2.48109039/), r_H3)), &
!71
         efp_water(water_atom(O1, (/-6.64774882,     -0.97504130,      3.26879224/), r_O1) , &
                   water_atom(H2, (/-5.92967378,     -0.73002697,      3.83024595/), r_H2) , &
                   water_atom(H3, (/-6.27042560,     -1.40307297,      2.51693886/), r_H3)), &
!72
         efp_water(water_atom(O1, (/-2.06164423,     -5.19155807,     -5.26870789/), r_O1) , &
                   water_atom(H2, (/-1.26012270,     -5.31256359,     -4.78518193/), r_H2) , &
                   water_atom(H3, (/-2.65046682,     -4.71746678,     -4.70354819/), r_H3)), &
!73
         efp_water(water_atom(O1, (/ 6.39120755,      2.07918606,     -3.19958758/), r_O1) , &
                   water_atom(H2, (/ 7.22925001,      1.98967718,     -3.62449671/), r_H2) , &
                   water_atom(H3, (/ 6.56606865,      2.29702752,     -2.29800631/), r_H3)), &
!74
         efp_water(water_atom(O1, (/ 2.85959401,      7.49381037,     -3.32163598/), r_O1) , &
                   water_atom(H2, (/ 3.77798105,      7.31786475,     -3.45002521/), r_H2) , &
                   water_atom(H3, (/ 2.77464867,      8.42775440,     -3.21481527/), r_H3)), &
!75
         efp_water(water_atom(O1, (/ 6.61641330,     -3.17621791,     -0.55809594/), r_O1) , &
                   water_atom(H2, (/ 5.82759161,     -3.45798248,     -0.99402268/), r_H2) , &
                   water_atom(H3, (/ 6.89915212,     -2.38765183,     -0.99263138/), r_H3)), &
!76
         efp_water(water_atom(O1, (/-6.13417524,     -3.11151560,     -2.99801140/), r_O1) , &
                   water_atom(H2, (/-6.30212195,     -3.49186831,     -2.15065993/), r_H2) , &
                   water_atom(H3, (/-6.72114394,     -3.53255649,     -3.60552265/), r_H3)), &
!77
         efp_water(water_atom(O1, (/ 6.57689316,      3.11330060,     -0.83384032/), r_O1) , &
                   water_atom(H2, (/ 6.59977558,      2.51008049,     -0.10825064/), r_H2) , &
                   water_atom(H3, (/ 5.88666289,      3.73076355,     -0.65163303/), r_H3)), &
!78
         efp_water(water_atom(O1, (/-4.05052493,     -6.54272699,      0.30366097/), r_O1) , &
                   water_atom(H2, (/-3.88865789,     -7.39427826,      0.67720814/), r_H2) , &
                   water_atom(H3, (/-4.83089630,     -6.20916863,      0.71674513/), r_H3)), &
!79
         efp_water(water_atom(O1, (/-6.90299513,      1.53108610,     -4.17119661/), r_O1) , &
                   water_atom(H2, (/-7.38861342,      1.82415398,     -3.41676692/), r_H2) , &
                   water_atom(H3, (/-6.58412397,      2.30564817,     -4.60622902/), r_H3)), &
!80
         efp_water(water_atom(O1, (/ 4.83455611,     -3.77643124,      4.95698528/), r_O1) , &
                   water_atom(H2, (/ 4.96062844,     -4.28860034,      4.17425461/), r_H2) , &
                   water_atom(H3, (/ 4.01077315,     -3.32666570,      4.85713069/), r_H3))  &
       /

  data (water_array1(i), i=81,100) &
       /                                                        &
!81
         efp_water(water_atom(O1, (/-4.82206823,      2.35529815,      3.88152675/), r_O1) , &
                   water_atom(H2, (/-4.09383635,      1.76936022,      3.75025009/), r_H2) , &
                   water_atom(H3, (/-4.63163324,      2.85467557,      4.65949516/), r_H3)), &
!82
         efp_water(water_atom(O1, (/-7.37190050,      2.92527486,     -1.72949287/), r_O1) , &
                   water_atom(H2, (/-7.73832401,      3.36926554,     -0.98150874/), r_H2) , &
                   water_atom(H3, (/-6.48744724,      3.24051094,     -1.82565056/), r_H3)), &
!83
         efp_water(water_atom(O1, (/-4.89649690,     -0.41930782,      5.33104080/), r_O1) , &
                   water_atom(H2, (/-4.24578980,      0.05189670,      4.83563407/), r_H2) , &
                   water_atom(H3, (/-5.47339506,      0.22729662,      5.70516951/), r_H3)), &
!84
         efp_water(water_atom(O1, (/ 2.34292862,      2.83414766,      7.02005922/), r_O1) , &
                   water_atom(H2, (/ 1.47813903,      3.19446904,      6.90522927/), r_H2) , &
                   water_atom(H3, (/ 2.43908277,      2.64756475,      7.94028702/), r_H3)), &
!85
         efp_water(water_atom(O1, (/-6.82561922,      0.42076645,     -0.72593295/), r_O1) , &
                   water_atom(H2, (/-7.14483204,      1.27783321,     -0.95921186/), r_H2) , &
                   water_atom(H3, (/-6.27546012,      0.13338031,     -1.43699413/), r_H3)), &
!86
         efp_water(water_atom(O1, (/ 3.86237429,      4.87614285,     -6.12418594/), r_O1) , &
                   water_atom(H2, (/ 4.31348778,      4.32930492,     -5.50101727/), r_H2) , &
                   water_atom(H3, (/ 3.76713131,      4.36678774,     -6.91308636/), r_H3)), &
!87
         efp_water(water_atom(O1, (/-2.76132936,     -4.02642375,      6.15270595/), r_O1) , &
                   water_atom(H2, (/-3.23233010,     -4.44855875,      6.85330452/), r_H2) , &
                   water_atom(H3, (/-3.40657526,     -3.60751708,      5.60584733/), r_H3)), &
!88
         efp_water(water_atom(O1, (/-1.25806696,     -3.97634935,     -7.46070247/), r_O1) , &
                   water_atom(H2, (/-1.92021666,     -4.09919133,     -8.12202566/), r_H2) , &
                   water_atom(H3, (/-1.57926222,     -4.38947435,     -6.67518430/), r_H3)), &
!89
         efp_water(water_atom(O1, (/-4.28443637,      1.90238062,     -6.09740147/), r_O1) , &
                   water_atom(H2, (/-4.94422135,      1.26644562,     -6.32356053/), r_H2) , &
                   water_atom(H3, (/-4.74122479,      2.68442002,     -5.83162437/), r_H3)), &
!90
         efp_water(water_atom(O1, (/-4.01919143,      5.35434512,      4.28117473/), r_O1) , &
                   water_atom(H2, (/-3.77645067,      6.04965836,      4.87150853/), r_H2) , &
                   water_atom(H3, (/-4.70461958,      5.69452920,      3.72859829/), r_H3)), &
!91
         efp_water(water_atom(O1, (/ 0.58638046,     -8.24139809,      1.82395136/), r_O1) , &
                   water_atom(H2, (/ 1.42289587,     -8.56127252,      2.12194311/), r_H2) , &
                   water_atom(H3, (/ 0.32824991,     -8.78992411,      1.10051150/), r_H3)), &
!92
         efp_water(water_atom(O1, (/ 0.96358962,      6.16195238,      5.99528700/), r_O1) , &
                   water_atom(H2, (/ 0.63300486,      5.27938706,      6.04696901/), r_H2) , &
                   water_atom(H3, (/ 0.55450968,      6.64625401,      6.69456033/), r_H3)), &
!93
         efp_water(water_atom(O1, (/ 2.20289379,     -7.02584467,     -3.91324354/), r_O1) , &
                   water_atom(H2, (/ 1.85981913,     -7.59969921,     -3.24700544/), r_H2) , &
                   water_atom(H3, (/ 2.17720141,     -7.50673428,     -4.72500897/), r_H3)), &
!94
         efp_water(water_atom(O1, (/ 0.11805696,      7.64906280,      3.93675754/), r_O1) , &
                   water_atom(H2, (/ 0.52741276,      7.01208353,      4.50028769/), r_H2) , &
                   water_atom(H3, (/-0.27100235,      7.16950601,      3.22293656/), r_H3)), &
!95
         efp_water(water_atom(O1, (/-5.60334370,      5.44520965,     -3.72205534/), r_O1) , &
                   water_atom(H2, (/-5.15350701,      5.05989510,     -2.98716771/), r_H2) , &
                   water_atom(H3, (/-5.30011204,      6.33654601,     -3.78874803/), r_H3)), &
!96
         efp_water(water_atom(O1, (/-1.52080902,     -7.98494513,      3.43404607/), r_O1) , &
                   water_atom(H2, (/-0.65363194,     -8.05997788,      3.06900862/), r_H2) , &
                   water_atom(H3, (/-1.51553845,     -7.22326471,      3.99144507/), r_H3)), &
!97
         efp_water(water_atom(O1, (/ 1.30419086,     -7.81337796,     -1.55231865/), r_O1) , &
                   water_atom(H2, (/ 0.42714245,     -7.99451101,     -1.25424061/), r_H2) , &
                   water_atom(H3, (/ 1.87859019,     -8.38540714,     -1.06886044/), r_H3)), &
!98
         efp_water(water_atom(O1, (/-2.83376765,      8.09098689,     -1.65740357/), r_O1) , &
                   water_atom(H2, (/-3.75671280,      8.22100063,     -1.50858268/), r_H2) , &
                   water_atom(H3, (/-2.68293833,      8.24425145,     -2.57644616/), r_H3)), &
!99
         efp_water(water_atom(O1, (/-1.86201113,      6.44006142,     -5.68509503/), r_O1) , &
                   water_atom(H2, (/-1.52117730,      5.78524687,     -5.09693785/), r_H2) , &
                   water_atom(H3, (/-1.25257686,      6.50343899,     -6.40304441/), r_H3)), &
!100
         efp_water(water_atom(O1, (/-1.30422051,     -5.91486647,      5.10621506/), r_O1) , &
                   water_atom(H2, (/-0.44270786,     -5.56442815,      5.26704200/), r_H2) , &
                   water_atom(H3, (/-1.91814805,     -5.22377322,      5.29690642/), r_H3))  &
       /

  data (water_array1(i), i=101,120) &
       /                                                     &
!101
         efp_water(water_atom(O1, (/-6.25291784,     -5.66542378,      1.43676060/), r_O1) , &
                   water_atom(H2, (/-6.76256262,     -6.45206331,      1.32568294/), r_H2) , &
                   water_atom(H3, (/-6.24879323,     -5.47307612,      2.36080806/), r_H3)), &
!102
         efp_water(water_atom(O1, (/-5.98787098,     -4.72261081,      3.87166236/), r_O1) , &
                   water_atom(H2, (/-5.24087473,     -4.14792333,      3.92272237/), r_H2) , &
                   water_atom(H3, (/-6.73795829,     -4.21714918,      4.14141772/), r_H3)), &
!103
         efp_water(water_atom(O1, (/-5.75846694,      4.00465386,      2.14397542/), r_O1) , &
                   water_atom(H2, (/-5.27371818,      3.28523440,      2.51590581/), r_H2) , &
                   water_atom(H3, (/-6.23035841,      4.40968357,      2.85400891/), r_H3)), &
!104
         efp_water(water_atom(O1, (/ 4.26644009,      1.52885398,     -6.43921510/), r_O1) , &
                   water_atom(H2, (/ 4.55093161,      0.99148044,     -5.71729453/), r_H2) , &
                   water_atom(H3, (/ 3.43678563,      1.90271761,     -6.18865207/), r_H3)), &
!105
         efp_water(water_atom(O1, (/ 5.20511887,      6.45299038,     -3.54797167/), r_O1) , &
                   water_atom(H2, (/ 5.81615847,      6.72699298,     -4.21313433/), r_H2) , &
                   water_atom(H3, (/ 5.05265924,      5.53132328,     -3.68274438/), r_H3)), &
!106
         efp_water(water_atom(O1, (/ 1.48100917,      2.46752909,     -8.26279982/), r_O1) , &
                   water_atom(H2, (/ 1.59345005,      2.64406290,     -7.34243173/), r_H2) , &
                   water_atom(H3, (/ 0.60037021,      2.72609987,     -8.48304171/), r_H3)), &
!107
         efp_water(water_atom(O1, (/ 5.86534492,     -1.81254030,     -5.42156771/), r_O1) , &
                   water_atom(H2, (/ 5.66091664,     -1.00464968,     -4.97839761/), r_H2) , &
                   water_atom(H3, (/ 6.69359743,     -2.10799751,     -5.07866401/), r_H3)), &
!108
         efp_water(water_atom(O1, (/ 3.61492275,      6.28957250,      5.83298140/), r_O1) , &
                   water_atom(H2, (/ 3.84827402,      5.37876878,      5.75014354/), r_H2) , &
                   water_atom(H3, (/ 2.67789246,      6.32071474,      5.94198805/), r_H3)), &
!109
         efp_water(water_atom(O1, (/-0.56769507,     -2.56908692,      6.28259605/), r_O1) , &
                   water_atom(H2, (/-1.42451710,     -2.95509841,      6.37051082/), r_H2) , &
                   water_atom(H3, (/ 0.02836968,     -3.27379546,      6.08518712/), r_H3)), &
!110
         efp_water(water_atom(O1, (/ 0.20647562,      5.88112371,     -7.30104535/), r_O1) , &
                   water_atom(H2, (/ 0.77292772,      6.17755150,     -7.99540729/), r_H2) , &
                   water_atom(H3, (/ 0.72882051,      5.86001669,     -6.51517543/), r_H3)), &
!111
         efp_water(water_atom(O1, (/ 0.19784970,     -0.81081666,      8.24956252/), r_O1) , &
                   water_atom(H2, (/ 0.00692373,     -1.49930004,      7.63278389/), r_H2) , &
                   water_atom(H3, (/ 1.13720250,     -0.72224627,      8.27505909/), r_H3)), &
!112
         efp_water(water_atom(O1, (/-2.28321591,      3.53399075,      7.41003231/), r_O1) , &
                   water_atom(H2, (/-3.12443674,      3.40877736,      7.00070439/), r_H2) , &
                   water_atom(H3, (/-2.32218168,      4.36138482,      7.86255735/), r_H3)), &
!113
         efp_water(water_atom(O1, (/ 8.19137541,     -1.74114049,      2.91208222/), r_O1) , &
                   water_atom(H2, (/ 7.44606069,     -2.28406475,      2.71060924/), r_H2) , &
                   water_atom(H3, (/ 8.04018400,     -1.38091638,      3.77131471/), r_H3)), &
!114
         efp_water(water_atom(O1, (/-7.18545638,     -4.57301414,     -0.98641564/), r_O1) , &
                   water_atom(H2, (/-6.82396242,     -4.77622315,     -0.13853891/), r_H2) , &
                   water_atom(H3, (/-7.60302393,     -5.35969051,     -1.29890793/), r_H3)), &
!115
         efp_water(water_atom(O1, (/-1.01133734,      3.55402005,     -7.84732713/), r_O1) , &
                   water_atom(H2, (/-1.87723567,      3.63035105,     -8.21512181/), r_H2) , &
                   water_atom(H3, (/-0.70190035,      4.43299596,     -7.69723084/), r_H3)), &
!116
         efp_water(water_atom(O1, (/-1.02838453,     -8.75335927,     -0.52275062/), r_O1) , &
                   water_atom(H2, (/-1.53853341,     -8.96014887,     -1.28948423/), r_H2) , &
                   water_atom(H3, (/-1.64184923,     -8.53439115,      0.16031265/), r_H3)), &
!117
         efp_water(water_atom(O1, (/ 1.00764554,     -4.86928498,      6.23193497/), r_O1) , &
                   water_atom(H2, (/ 0.97359620,     -5.20572749,      7.11314225/), r_H2) , &
                   water_atom(H3, (/ 1.49969898,     -5.49297446,      5.72224982/), r_H3)), &
!118
         efp_water(water_atom(O1, (/ 7.58721335,     -0.44462886,      5.21070649/), r_O1) , &
                   water_atom(H2, (/ 8.24719616,      0.14233087,      5.54354672/), r_H2) , &
                   water_atom(H3, (/ 6.81658324,      0.07634909,      5.05073496/), r_H3)), &
!119
         efp_water(water_atom(O1, (/ 8.31982924,      4.38616073,     -2.42935060/), r_O1) , &
                   water_atom(H2, (/ 8.65816776,      5.12311067,     -1.94632299/), r_H2) , &
                   water_atom(H3, (/ 7.80496019,      3.87782202,     -1.82324054/), r_H3)), &
!120
         efp_water(water_atom(O1, (/ 2.89648032,     -0.48622468,      8.28159145/), r_O1) , &
                   water_atom(H2, (/ 3.43563223,     -1.20608077,      8.56794718/), r_H2) , &
                   water_atom(H3, (/ 3.12524672,      0.25481518,      8.81955386/), r_H3))  &
       /

  data (water_array1(i), i=121,140) &
       /                                                     &
!121
         efp_water(water_atom(O1, (/-5.96968445,      3.65823114,     -5.46989371/), r_O1) , &
                   water_atom(H2, (/-5.75512355,      4.42306268,     -4.96011469/), r_H2) , &
                   water_atom(H3, (/-6.48756886,      3.95584364,     -6.20071552/), r_H3)), &
!122
         efp_water(water_atom(O1, (/-7.89318987,     -1.88497556,     -0.07672653/), r_O1) , &
                   water_atom(H2, (/-7.75466197,     -1.00330841,     -0.38389169/), r_H2) , &
                   water_atom(H3, (/-8.11748572,     -2.40356107,     -0.83279725/), r_H3)), &
!123
         efp_water(water_atom(O1, (/ 8.44633406,      0.17333166,      1.24077027/), r_O1) , &
                   water_atom(H2, (/ 8.58127129,     -0.12711584,      0.35624519/), r_H2) , &
                   water_atom(H3, (/ 8.48021257,     -0.59011786,      1.79473605/), r_H3)), &
!124
         efp_water(water_atom(O1, (/ 5.78306824,     -6.67817847,      0.05854016/), r_O1) , &
                   water_atom(H2, (/ 6.58070790,     -6.23176692,      0.29369801/), r_H2) , &
                   water_atom(H3, (/ 5.20114919,     -6.02509440,     -0.29593935/), r_H3)), &
!125
         efp_water(water_atom(O1, (/-6.06193781,     -6.87318167,     -3.55105440/), r_O1) , &
                   water_atom(H2, (/-5.84839532,     -7.15917758,     -4.42483102/), r_H2) , &
                   water_atom(H3, (/-5.25903350,     -6.56390314,     -3.16302325/), r_H3)), &
!126
         efp_water(water_atom(O1, (/-1.84248661,      0.77407999,      7.76095006/), r_O1) , &
                   water_atom(H2, (/-1.06783408,      0.47095680,      8.20693770/), r_H2) , &
                   water_atom(H3, (/-1.88203804,      1.70818677,      7.89040615/), r_H3)), &
!127
         efp_water(water_atom(O1, (/ 2.25533355,     -6.69677277,      4.60979800/), r_O1) , &
                   water_atom(H2, (/ 2.38505080,     -6.02847665,      3.95601421/), r_H2) , &
                   water_atom(H3, (/ 2.29196727,     -7.52622082,      4.16085228/), r_H3)), &
!128
         efp_water(water_atom(O1, (/-5.25037349,     -1.39988552,     -7.77926781/), r_O1) , &
                   water_atom(H2, (/-4.39045743,     -1.01157870,     -7.80457866/), r_H2) , &
                   water_atom(H3, (/-5.14223587,     -2.27719200,     -7.44833268/), r_H3)), &
!129
         efp_water(water_atom(O1, (/-2.92928920,     -0.08272120,     -7.58245236/), r_O1) , &
                   water_atom(H2, (/-3.02837267,      0.72882846,     -7.11080919/), r_H2) , &
                   water_atom(H3, (/-2.16948789,     -0.51343978,     -7.22459197/), r_H3)), &
!130
         efp_water(water_atom(O1, (/ 8.96665801,      2.46130988,      2.68960316/), r_O1) , &
                   water_atom(H2, (/ 8.97817240,      3.18319992,      2.08163885/), r_H2) , &
                   water_atom(H3, (/ 8.90077631,      1.67345310,      2.17402607/), r_H3)), &
!131
         efp_water(water_atom(O1, (/-2.95615921,     -8.61904025,      1.40513742/), r_O1) , &
                   water_atom(H2, (/-3.30645749,     -9.48747905,      1.52339012/), r_H2) , &
                   water_atom(H3, (/-2.56802536,     -8.36867151,      2.22827328/), r_H3)), &
!132
         efp_water(water_atom(O1, (/ 7.82359234,     -5.25917077,      0.91475688/), r_O1) , &
                   water_atom(H2, (/ 7.79863222,     -5.10749287,      1.84602165/), r_H2) , &
                   water_atom(H3, (/ 7.73469145,     -4.41594997,      0.50009161/), r_H3)), &
!133
         efp_water(water_atom(O1, (/ 6.08164287,      6.77713072,     -1.07456692/), r_O1) , &
                   water_atom(H2, (/ 7.01707616,      6.69669757,     -0.97775121/), r_H2) , &
                   water_atom(H3, (/ 5.89961304,      6.73020970,     -1.99952239/), r_H3)), &
!134
         efp_water(water_atom(O1, (/-6.68986848,      6.85045923,      0.38946424/), r_O1) , &
                   water_atom(H2, (/-5.83957179,      6.46938939,      0.23893302/), r_H2) , &
                   water_atom(H3, (/-6.90395348,      7.34962822,     -0.38246677/), r_H3)), &
!135
         efp_water(water_atom(O1, (/-2.83226198,     -8.87665545,     -2.54788679/), r_O1) , &
                   water_atom(H2, (/-3.07439337,     -7.96675406,     -2.61368952/), r_H2) , &
                   water_atom(H3, (/-2.63709275,     -9.16946784,     -3.42370301/), r_H3)), &
!136
         efp_water(water_atom(O1, (/ 3.14306460,     -3.47375016,     -7.75109450/), r_O1) , &
                   water_atom(H2, (/ 2.43737146,     -4.10047241,     -7.74137420/), r_H2) , &
                   water_atom(H3, (/ 2.75881084,     -2.63565900,     -7.95316904/), r_H3)), &
!137
         efp_water(water_atom(O1, (/ 4.60239409,      6.83593655,      3.54584452/), r_O1) , &
                   water_atom(H2, (/ 4.60864434,      7.72074409,      3.21727929/), r_H2) , &
                   water_atom(H3, (/ 4.25815050,      6.87501057,      4.42382271/), r_H3)), &
!138
         efp_water(water_atom(O1, (/-5.71712345,     -6.54205982,      5.72275347/), r_O1) , &
                   water_atom(H2, (/-5.68301986,     -7.41107697,      5.35597669/), r_H2) , &
                   water_atom(H3, (/-5.78434554,     -5.94338446,      4.99615468/), r_H3)), &
!139
         efp_water(water_atom(O1, (/-7.80732100,      4.31837898,      0.42890202/), r_O1) , &
                   water_atom(H2, (/-7.66392814,      5.24793042,      0.34984143/), r_H2) , &
                   water_atom(H3, (/-7.14809872,      3.99150815,      1.02005356/), r_H3)), &
!140
         efp_water(water_atom(O1, (/-4.94928059,     -3.82105907,     -6.78155510/), r_O1) , &
                   water_atom(H2, (/-4.52964392,     -3.78649850,     -5.93681375/), r_H2) , &
                   water_atom(H3, (/-5.74228846,     -4.32069588,     -6.67028857/), r_H3))  &
       /

  data (water_array1(i), i=141,160) &
       /                                                     &
!141
         efp_water(water_atom(O1, (/ 6.57257072,      5.03348004,      3.16737269/), r_O1) , &
                   water_atom(H2, (/ 6.04715624,      4.27431889,      2.97117066/), r_H2) , &
                   water_atom(H3, (/ 5.97270839,      5.74325439,      3.33248168/), r_H3)), &
!142
         efp_water(water_atom(O1, (/ 2.82632665,     -8.81885440,      3.03562896/), r_O1) , &
                   water_atom(H2, (/ 3.65762807,     -8.48196981,      2.74181784/), r_H2) , &
                   water_atom(H3, (/ 2.99542601,     -9.66756493,      3.41242481/), r_H3)), &
!143
         efp_water(water_atom(O1, (/ 7.30960443,     -5.04891818,      3.58618347/), r_O1) , &
                   water_atom(H2, (/ 7.18292953,     -5.88832171,      3.99878034/), r_H2) , &
                   water_atom(H3, (/ 7.61690216,     -4.46046173,      4.25711998/), r_H3)), &
!144
         efp_water(water_atom(O1, (/ 6.54637872,      2.91502567,     -6.58817590/), r_O1) , &
                   water_atom(H2, (/ 5.70523057,      2.49357650,     -6.51249113/), r_H2) , &
                   water_atom(H3, (/ 6.62237589,      3.21049183,     -7.48137346/), r_H3)), &
!145
         efp_water(water_atom(O1, (/ 1.13432668,     -5.16365729,     -7.26593697/), r_O1) , &
                   water_atom(H2, (/ 0.96533273,     -5.13537195,     -6.33775619/), r_H2) , &
                   water_atom(H3, (/ 0.33871352,     -4.89193122,     -7.69493917/), r_H3)), &
!146
         efp_water(water_atom(O1, (/-5.52242810,     -8.99535818,     -1.58035948/), r_O1) , &
                   water_atom(H2, (/-5.76096552,     -8.44595758,     -2.30982678/), r_H2) , &
                   water_atom(H3, (/-4.61623598,     -9.22960951,     -1.70213687/), r_H3)), &
!147
         efp_water(water_atom(O1, (/-4.42564101,     -1.44913291,      7.69773348/), r_O1) , &
                   water_atom(H2, (/-4.34824996,     -1.05364461,      6.84422513/), r_H2) , &
                   water_atom(H3, (/-3.60934017,     -1.29174298,      8.14468805/), r_H3)), &
!148
         efp_water(water_atom(O1, (/-4.73795761,      3.26430377,      6.39944364/), r_O1) , &
                   water_atom(H2, (/-5.24794758,      4.05769770,      6.43570820/), r_H2) , &
                   water_atom(H3, (/-5.30864690,      2.56242487,      6.66880789/), r_H3)), &
!149
         efp_water(water_atom(O1, (/ 8.66884643,     -0.67794792,     -3.73888206/), r_O1) , &
                   water_atom(H2, (/ 8.63242794,     -1.50068358,     -4.20003077/), r_H2) , &
                   water_atom(H3, (/ 8.52862487,     -0.00057820,     -4.38105557/), r_H3)), &
!150
         efp_water(water_atom(O1, (/ 5.90024740,     -2.80322142,     -7.89491178/), r_O1) , &
                   water_atom(H2, (/ 5.05147581,     -3.17561787,     -8.07318945/), r_H2) , &
                   water_atom(H3, (/ 5.79963632,     -2.24887437,     -7.13764211/), r_H3)), &
!151
         efp_water(water_atom(O1, (/-7.38378733,      1.69880370,      4.04206224/), r_O1) , &
                   water_atom(H2, (/-7.63891958,      0.93196379,      3.55447324/), r_H2) , &
                   water_atom(H3, (/-6.47486196,      1.85700412,      3.84280070/), r_H3)), &
!152
         efp_water(water_atom(O1, (/ 7.83586119,      4.78305376,     -5.04066694/), r_O1) , &
                   water_atom(H2, (/ 7.91074709,      4.55400601,     -4.12808428/), r_H2) , &
                   water_atom(H3, (/ 7.27123139,      4.13864132,     -5.43664804/), r_H3)), &
!153
         efp_water(water_atom(O1, (/ 4.58019801,     -4.69918689,      7.35337282/), r_O1) , &
                   water_atom(H2, (/ 4.11008220,     -5.51694796,      7.38708782/), r_H2) , &
                   water_atom(H3, (/ 4.66931136,     -4.47492637,      6.44087999/), r_H3)), &
!154
         efp_water(water_atom(O1, (/ 3.60367634,      3.99383642,     -8.54811192/), r_O1) , &
                   water_atom(H2, (/ 2.94591506,      3.34411315,     -8.73806496/), r_H2) , &
                   water_atom(H3, (/ 4.40805805,      3.68663962,     -8.93475626/), r_H3)), &
!155
         efp_water(water_atom(O1, (/ 8.27425718,      1.69421998,     -4.99790182/), r_O1) , &
                   water_atom(H2, (/ 7.76476380,      2.00019158,     -5.73114803/), r_H2) , &
                   water_atom(H3, (/ 9.07536640,      2.19332033,     -4.99710798/), r_H3)), &
!156
         efp_water(water_atom(O1, (/-6.43053642,      6.72119188,      3.15150129/), r_O1) , &
                   water_atom(H2, (/-7.12649966,      6.11001955,      3.33310570/), r_H2) , &
                   water_atom(H3, (/-6.55238017,      7.01481224,      2.26278287/), r_H3)), &
!157
         efp_water(water_atom(O1, (/-7.70937303,     -3.01656760,      4.76267599/), r_O1) , &
                   water_atom(H2, (/-7.57221326,     -2.22351885,      4.26958717/), r_H2) , &
                   water_atom(H3, (/-7.42962029,     -2.83838682,      5.64634803/), r_H3)), &
!158
         efp_water(water_atom(O1, (/ 5.00863044,      7.03588040,     -7.11748588/), r_O1) , &
                   water_atom(H2, (/ 5.71154441,      7.29779078,     -6.54459878/), r_H2) , &
                   water_atom(H3, (/ 4.57524222,      6.30861343,     -6.70019287/), r_H3)), &
!159
         efp_water(water_atom(O1, (/ 7.05820365,     -5.07628150,     -3.23015753/), r_O1) , &
                   water_atom(H2, (/ 7.26984962,     -5.07722673,     -2.31023227/), r_H2) , &
                   water_atom(H3, (/ 6.15713739,     -5.34808125,     -3.30191250/), r_H3)), &
!160
         efp_water(water_atom(O1, (/ 7.24472852,     -3.10675966,      5.74812892/), r_O1) , &
                   water_atom(H2, (/ 7.47360317,     -2.19233246,      5.69998258/), r_H2) , &
                   water_atom(H3, (/ 6.34007108,     -3.17406704,      5.48748575/), r_H3))  &
       /

  data (water_array1(i), i=161,180) &
       /                                                     &
!161
         efp_water(water_atom(O1, (/-7.74388540,     -4.93515256,     -4.23589992/), r_O1) , &
                   water_atom(H2, (/-8.54658582,     -4.99831514,     -3.74339399/), r_H2) , &
                   water_atom(H3, (/-7.17641484,     -5.61941923,     -3.91867068/), r_H3)), &
!162
         efp_water(water_atom(O1, (/ 5.00731378,     -7.60163762,      2.39247312/), r_O1) , &
                   water_atom(H2, (/ 5.35622953,     -7.52265033,      1.51903412/), r_H2) , &
                   water_atom(H3, (/ 5.74516100,     -7.71439896,      2.97018152/), r_H3)), &
!163
         efp_water(water_atom(O1, (/ 3.18392411,      1.93184661,      9.39857968/), r_O1) , &
                   water_atom(H2, (/ 2.81972439,      2.03277059,     10.26347878/), r_H2) , &
                   water_atom(H3, (/ 3.98975420,      2.42308900,      9.38453492/), r_H3)), &
!164
         efp_water(water_atom(O1, (/-2.28829156,      6.12777510,      8.35654093/), r_O1) , &
                   water_atom(H2, (/-1.50313321,      6.64952372,      8.30981622/), r_H2) , &
                   water_atom(H3, (/-2.68197544,      6.30462861,      9.19595355/), r_H3)), &
!165
         efp_water(water_atom(O1, (/-5.27689655,      7.98677251,     -4.04453660/), r_O1) , &
                   water_atom(H2, (/-5.63186674,      8.07235524,     -4.91491003/), r_H2) , &
                   water_atom(H3, (/-4.39706801,      8.32707821,     -4.07571492/), r_H3)), &
!166
         efp_water(water_atom(O1, (/-4.10704968,      6.01505464,     -7.12052193/), r_O1) , &
                   water_atom(H2, (/-4.05201015,      5.15447091,     -7.50424730/), r_H2) , &
                   water_atom(H3, (/-3.37012877,      6.09945343,     -6.53682647/), r_H3)), &
!167
         efp_water(water_atom(O1, (/ 4.36599996,     -2.58891739,      8.92150382/), r_O1) , &
                   water_atom(H2, (/ 5.27232127,     -2.48834845,      9.16511966/), r_H2) , &
                   water_atom(H3, (/ 4.31761239,     -3.34924557,      8.36433343/), r_H3)), &
!168
         efp_water(water_atom(O1, (/-8.17719035,      5.14699635,     -3.41126509/), r_O1) , &
                   water_atom(H2, (/-8.24384621,      4.40625792,     -2.83010960/), r_H2) , &
                   water_atom(H3, (/-7.25810864,      5.27397455,     -3.58459177/), r_H3)), &
!169
         efp_water(water_atom(O1, (/-3.47335224,     -4.68300437,     -8.74947034/), r_O1) , &
                   water_atom(H2, (/-4.09986127,     -4.32085233,     -8.14348847/), r_H2) , &
                   water_atom(H3, (/-3.78752145,     -4.48252306,     -9.61663840/), r_H3)), &
!170
         efp_water(water_atom(O1, (/-3.26148800,     -6.86423854,     -6.95753151/), r_O1) , &
                   water_atom(H2, (/-2.76021305,     -6.39031782,     -6.31332748/), r_H2) , &
                   water_atom(H3, (/-3.23148706,     -6.35594928,     -7.75227693/), r_H3)), &
!171
         efp_water(water_atom(O1, (/-2.66063552,      8.49812911,     -4.30259212/), r_O1) , &
                   water_atom(H2, (/-2.33800311,      7.80822081,     -4.86009810/), r_H2) , &
                   water_atom(H3, (/-2.35121745,      9.31071083,     -4.66981098/), r_H3)), &
!172
         efp_water(water_atom(O1, (/-6.40968830,      0.53992915,     -6.70597580/), r_O1) , &
                   water_atom(H2, (/-6.98119211,      0.50924544,     -5.95542722/), r_H2) , &
                   water_atom(H3, (/-6.23443457,     -0.35477189,     -6.95026083/), r_H3)), &
!173
         efp_water(water_atom(O1, (/-1.34546799,      8.91150039,      5.58310673/), r_O1) , &
                   water_atom(H2, (/-0.83807533,      8.56234846,      4.86789817/), r_H2) , &
                   water_atom(H3, (/-1.54279163,      9.80799541,      5.36347598/), r_H3)), &
!174
         efp_water(water_atom(O1, (/ 3.25929272,     -6.98515523,      7.05304852/), r_O1) , &
                   water_atom(H2, (/ 2.48439290,     -7.08605144,      7.58241401/), r_H2) , &
                   water_atom(H3, (/ 2.97068770,     -6.93975728,      6.15553812/), r_H3)), &
!175
         efp_water(water_atom(O1, (/-3.37426019,      7.24091629,      6.07132000/), r_O1) , &
                   water_atom(H2, (/-3.13306521,      6.81915954,      6.88053410/), r_H2) , &
                   water_atom(H3, (/-2.68361256,      7.85094253,      5.86699063/), r_H3)), &
!176
         efp_water(water_atom(O1, (/-3.68397021,      3.47780563,     -8.23105654/), r_O1) , &
                   water_atom(H2, (/-4.27345036,      3.33306861,     -8.95383655/), r_H2) , &
                   water_atom(H3, (/-3.87890957,      2.81733381,     -7.58556209/), r_H3)), &
!177
         efp_water(water_atom(O1, (/ 2.07330649,     -7.64615140,     -6.47979907/), r_O1) , &
                   water_atom(H2, (/ 2.80431578,     -7.97650263,     -6.97718459/), r_H2) , &
                   water_atom(H3, (/ 1.74147445,     -6.89785766,     -6.94970620/), r_H3)), &
!178
         efp_water(water_atom(O1, (/ 2.79613242,      0.03617807,     -8.17114197/), r_O1) , &
                   water_atom(H2, (/ 3.57084828,      0.42934319,     -7.80221533/), r_H2) , &
                   water_atom(H3, (/ 2.25664617,      0.74365999,     -8.48625843/), r_H3)), &
!179
         efp_water(water_atom(O1, (/ 8.09307114,      3.77202477,      4.97861916/), r_O1) , &
                   water_atom(H2, (/ 7.52854832,      4.31550139,      4.45246888/), r_H2) , &
                   water_atom(H3, (/ 8.43707703,      3.10457422,      4.40674290/), r_H3)), &
!180
         efp_water(water_atom(O1, (/ 7.04023967,      7.18379150,     -5.36414012/), r_O1) , &
                   water_atom(H2, (/ 7.68187042,      7.85786752,     -5.20661923/), r_H2) , &
                   water_atom(H3, (/ 7.50459653,      6.36205469,     -5.36467052/), r_H3))  &
       /

  data (water_array1(i), i=181,200) &
       /                                                     &
!181
         efp_water(water_atom(O1, (/ 5.30377842,      5.98263662,      8.04719679/), r_O1) , &
                   water_atom(H2, (/ 6.00435226,      5.50418050,      7.63348739/), r_H2) , &
                   water_atom(H3, (/ 4.72554997,      6.26736129,      7.35765942/), r_H3)), &
!182
         efp_water(water_atom(O1, (/-6.51971083,      1.28009814,      6.48711365/), r_O1) , &
                   water_atom(H2, (/-7.17639718,      1.44907121,      5.83053480/), r_H2) , &
                   water_atom(H3, (/-6.94526013,      0.79940264,      7.17900528/), r_H3)), &
!183
         efp_water(water_atom(O1, (/-0.32333493,      7.81806477,      7.72144355/), r_O1) , &
                   water_atom(H2, (/ 0.22618398,      8.33598400,      8.28771755/), r_H2) , &
                   water_atom(H3, (/-0.65953057,      8.40253842,      7.06095759/), r_H3)), &
!184
         efp_water(water_atom(O1, (/-7.14851507,     -7.83212467,      0.32329059/), r_O1) , &
                   water_atom(H2, (/-7.73450674,     -7.58458507,     -0.37399364/), r_H2) , &
                   water_atom(H3, (/-6.45483643,     -8.33605566,     -0.07136074/), r_H3)), &
!185
         efp_water(water_atom(O1, (/ 6.77151080,     -7.59602173,      4.35895184/), r_O1) , &
                   water_atom(H2, (/ 6.35484155,     -7.73888217,      5.19372409/), r_H2) , &
                   water_atom(H3, (/ 7.56444297,     -8.10801320,      4.36039720/), r_H3)), &
!186
         efp_water(water_atom(O1, (/-5.69042932,     -7.32331157,     -6.05345140/), r_O1) , &
                   water_atom(H2, (/-4.82271824,     -7.24863915,     -6.41729104/), r_H2) , &
                   water_atom(H3, (/-6.23372552,     -6.70774716,     -6.51905679/), r_H3)), &
!187
         efp_water(water_atom(O1, (/-7.50373658,      4.48673649,      4.15736229/), r_O1) , &
                   water_atom(H2, (/-7.26478028,      4.58545559,      5.06512542/), r_H2) , &
                   water_atom(H3, (/-7.83479434,      3.60844415,      4.05795281/), r_H3)), &
!188
         efp_water(water_atom(O1, (/-7.61540734,      4.21639420,     -7.44031359/), r_O1) , &
                   water_atom(H2, (/-8.47474607,      4.12197660,     -7.06149304/), r_H2) , &
                   water_atom(H3, (/-7.52014619,      3.52255761,     -8.07308196/), r_H3)), &
!189
         efp_water(water_atom(O1, (/-4.01058012,     -5.80445736,      7.67853763/), r_O1) , &
                   water_atom(H2, (/-4.58605109,     -5.63847832,      8.40803419/), r_H2) , &
                   water_atom(H3, (/-4.54393244,     -6.16677494,      6.98923543/), r_H3)), &
!190
         efp_water(water_atom(O1, (/-1.53963635,     -7.09216392,      7.56427542/), r_O1) , &
                   water_atom(H2, (/-1.46420819,     -6.90664939,      6.64190129/), r_H2) , &
                   water_atom(H3, (/-2.36002684,     -6.72245304,      7.84914745/), r_H3)), &
!191
         efp_water(water_atom(O1, (/ 5.20824223,      3.50683585,      9.20239204/), r_O1) , &
                   water_atom(H2, (/ 4.95357174,      4.40211010,      9.04584868/), r_H2) , &
                   water_atom(H3, (/ 5.92643485,      3.32298141,      8.61820453/), r_H3)), &
!192
         efp_water(water_atom(O1, (/ 8.03417296,     -3.06738764,     -4.74904858/), r_O1) , &
                   water_atom(H2, (/ 8.26400849,     -3.48839522,     -5.56195631/), r_H2) , &
                   water_atom(H3, (/ 7.78276714,     -3.75240803,     -4.15035782/), r_H3)), &
!193
         efp_water(water_atom(O1, (/ 8.75084100,      6.36017937,     -0.61188882/), r_O1) , &
                   water_atom(H2, (/ 8.82770921,      5.96263078,      0.24070841/), r_H2) , &
                   water_atom(H3, (/ 9.38559269,      7.05798513,     -0.6440778/), r_H3)), &
!194
         efp_water(water_atom(O1, (/-8.21566584,     -7.00818047,     -2.07609942/), r_O1) , &
                   water_atom(H2, (/-8.96575164,     -7.50429221,     -2.36268162/), r_H2) , &
                   water_atom(H3, (/-7.50812008,     -7.22296998,     -2.66270918/), r_H3)), &
!195
         efp_water(water_atom(O1, (/ 0.86593290,     -6.42518952,      8.37528360/), r_O1) , &
                   water_atom(H2, (/-0.00807377,     -6.76415130,      8.26530901/), r_H2) , &
                   water_atom(H3, (/ 0.94936330,     -6.17935012,      9.28274141/), r_H3)), &
!196
         efp_water(water_atom(O1, (/ 2.40812675,      6.41924459,     -8.74002406/), r_O1) , &
                   water_atom(H2, (/ 3.00581424,      7.06985682,     -8.40783283/), r_H2) , &
                   water_atom(H3, (/ 2.89383336,      5.61253563,     -8.80474100/), r_H3)), &
!197
         efp_water(water_atom(O1, (/-6.67207609,     -2.98703297,      7.24960356/), r_O1) , &
                   water_atom(H2, (/-6.74275764,     -3.87614755,      7.55840569/), r_H2) , &
                   water_atom(H3, (/-5.78383028,     -2.71287291,      7.41309158/), r_H3)), &
!198
         efp_water(water_atom(O1, (/-6.54568734,      2.20599115,     -8.78314648/), r_O1) , &
                   water_atom(H2, (/-6.73406078,      1.74291198,     -9.58374246/), r_H2) , &
                   water_atom(H3, (/-6.51336621,      1.55899193,     -8.09668756/), r_H3)), &
!199
         efp_water(water_atom(O1, (/-5.99006324,      7.57284207,      5.66080915/), r_O1) , &
                   water_atom(H2, (/-5.06059052,      7.68770230,      5.77813869/), r_H2) , &
                   water_atom(H3, (/-6.14495382,      7.54772904,      4.73007982/), r_H3)), &
!200
         efp_water(water_atom(O1, (/-7.24314589,     -5.18598559,     -6.76489669/), r_O1) , &
                   water_atom(H2, (/-7.70615378,     -5.03738209,     -5.95593340/), r_H2) , &
                   water_atom(H3, (/-7.85209367,     -5.00550305,     -7.46309811/), r_H3))  &
       /

  data (water_array1(i), i=201,216) &
       /                                                     &
!201
         efp_water(water_atom(O1, (/-6.58676789,      5.23155764,      6.57083519/), r_O1) , &
                   water_atom(H2, (/-6.39677852,      6.12370651,      6.32823926/), r_H2) , &
                   water_atom(H3, (/-7.12822712,      5.26982824,      7.34299796/), r_H3)), &
!202
         efp_water(water_atom(O1, (/ 7.13525304,     -4.50493866,      8.05326789/), r_O1) , &
                   water_atom(H2, (/ 7.43753088,     -3.99949167,      7.31568057/), r_H2) , &
                   water_atom(H3, (/ 6.20614115,     -4.62445542,      7.93775714/), r_H3)), &
!203
         efp_water(water_atom(O1, (/ 7.98065247,     -6.45060534,     -5.60084360/), r_O1) , &
                   water_atom(H2, (/ 7.93184149,     -5.71271876,     -6.18737577/), r_H2) , &
                   water_atom(H3, (/ 7.78923086,     -6.12323617,     -4.73651578/), r_H3)), &
!204
         efp_water(water_atom(O1, (/-6.93129308,     -0.46385411,      8.36907928/), r_O1) , &
                   water_atom(H2, (/-6.00892652,     -0.65588507,      8.42601634/), r_H2) , &
                   water_atom(H3, (/-7.36117743,     -1.26778601,      8.12459894/), r_H3)), &
!205
         efp_water(water_atom(O1, (/-5.92921252,      7.91609248,     -6.52962846/), r_O1) , &
                   water_atom(H2, (/-6.80368204,      7.70530494,     -6.81556281/), r_H2) , &
                   water_atom(H3, (/-5.36156079,      7.24659779,     -6.87664365/), r_H3)), &
!206
         efp_water(water_atom(O1, (/ 5.76463947,     -8.06721778,      6.71008527/), r_O1) , &
                   water_atom(H2, (/ 6.36448565,     -7.76782811,      7.37448589/), r_H2) , &
                   water_atom(H3, (/ 4.90618491,     -7.76945262,      6.96556624/), r_H3)), &
!207
         efp_water(water_atom(O1, (/ 6.22378871,      3.69205692,     -9.16649134/), r_O1) , &
                   water_atom(H2, (/ 6.62348211,      3.21989723,     -9.87936179/), r_H2) , &
                   water_atom(H3, (/ 6.33582585,      4.61053581,     -9.35282898/), r_H3)), &
!208
         efp_water(water_atom(O1, (/-7.21717191,      7.75556178,     -2.10018810/), r_O1) , &
                   water_atom(H2, (/-7.75145029,      7.07267687,     -2.47314195/), r_H2) , &
                   water_atom(H3, (/-6.55766298,      7.96509271,     -2.74207698/), r_H3)), &
!209
         efp_water(water_atom(O1, (/ 8.03307128,     -4.15908092,     -7.18612708/), r_O1) , &
                   water_atom(H2, (/ 7.25971397,     -3.76522091,     -7.55716018/), r_H2) , &
                   water_atom(H3, (/ 8.68910864,     -4.16031684,     -7.86472841/), r_H3)), &
!210
         efp_water(water_atom(O1, (/-9.50045797,      5.03745829,     -5.68988327/), r_O1) , &
                   water_atom(H2, (/-9.06821735,      5.05259061,     -4.85094403/), r_H2) , &
                   water_atom(H3, (/-10.4271722,      4.99140731,     -5.51679959/), r_H3)), &
!211
         efp_water(water_atom(O1, (/ 8.78232489,      5.21077765,      1.77794113/), r_O1) , &
                   water_atom(H2, (/ 9.43968951,      5.40877047,      2.42567161/), r_H2) , &
                   water_atom(H3, (/ 7.95004013,      5.22523927,      2.22287540/), r_H3)), &
!212
         efp_water(water_atom(O1, (/ 6.06631778,      6.27118966,     -9.40155719/), r_O1) , &
                   water_atom(H2, (/ 6.66625884,      6.93317648,     -9.70606161/), r_H2) , &
                   water_atom(H3, (/ 5.67257283,      6.60476114,     -8.61125633/), r_H3)), &
!213
         efp_water(water_atom(O1, (/-6.91585090,     -5.59071855,      7.87340507/), r_O1) , &
                   water_atom(H2, (/-7.73963059,     -5.90122594,      8.21376124/), r_H2) , &
                   water_atom(H3, (/-6.75442292,     -6.06875883,      7.07572504/), r_H3)), &
!214
         efp_water(water_atom(O1, (/ 7.24868154,      4.06992745,      7.48256396/), r_O1) , &
                   water_atom(H2, (/ 7.44932384,      3.81491379,      6.59622830/), r_H2) , &
                   water_atom(H3, (/ 8.07173148,      4.26355393,      7.90205528/), r_H3)), &
!215
         efp_water(water_atom(O1, (/-8.25070710,      6.87215115,     -7.35824029/), r_O1) , &
                   water_atom(H2, (/-8.84913820,      6.57991861,     -6.68939207/), r_H2) , &
                   water_atom(H3, (/-7.88427097,      6.09658442,     -7.75206749/), r_H3)), &
!216
         efp_water(water_atom(O1, (/ 7.66974287,     -7.06432521,      8.27562028/), r_O1) , &
                   water_atom(H2, (/ 7.89974140,     -7.28231428,      9.16469777/), r_H2) , &
                   water_atom(H3, (/ 7.64515721,     -6.12214155,      8.22499251/), r_H3))  &
       /

  data(r_vdw(i), i=1,20)      &
       /                      &
         vdw_data("H ",1.20), &
         vdw_data("HE",1.40), &
         vdw_data("LI",1.82), &
         vdw_data("BE",0.52), &
         vdw_data("B ",1.70), &
         vdw_data("C ",1.70), &
         vdw_data("N ",1.55), &
         vdw_data("O ",1.52), &
         vdw_data("F ",1.47), &
         vdw_data("NE",1.54), &
         vdw_data("NA",2.27), &
         vdw_data("MG",1.73), &
         vdw_data("AL",2.01), &
         vdw_data("SI",2.10), &
         vdw_data("P ",1.80), &
         vdw_data("S ",1.80), &
         vdw_data("CL",1.75), &
         vdw_data("AR",1.88), &
         vdw_data("K ",2.75), &
         vdw_data("CA",1.48)  &
       /

  data(r_vdw(i), i=21,40)     &
       /                      & 
         vdw_data("SC",2.15), &
         vdw_data("TI",2.19), &
         vdw_data("V ",1.99), &
         vdw_data("CR",2.01), &
         vdw_data("MN",2.01), &
         vdw_data("FE",2.00), &
         vdw_data("CO",1.99), &
         vdw_data("NI",1.63), &
         vdw_data("CU",1.49), &
         vdw_data("ZN",1.39), &
         vdw_data("GA",1.87), &
         vdw_data("GE",1.75), &
         vdw_data("AS",1.85), &
         vdw_data("SE",1.90), &
         vdw_data("BR",1.85), &
         vdw_data("KR",2.06), &
         vdw_data("RB",2.19), &
         vdw_data("SR",1.67), &
         vdw_data("Y ",2.66), &
         vdw_data("ZR",2.33)  &
       /

  data(r_vdw(i), i=41,60)     &
       /                      &
         vdw_data("NB",2.21), &
         vdw_data("MO",2.19), &
         vdw_data("TC",2.01), &
         vdw_data("RU",2.09), &
         vdw_data("RH",2.16), &
         vdw_data("PD",1.63), &
         vdw_data("AG",1.72), &
         vdw_data("CD",1.58), &
         vdw_data("IN",1.93), &
         vdw_data("SN",2.17), &
         vdw_data("SB",2.17), &
         vdw_data("TE",2.06), &
         vdw_data("I ",1.98), &
         vdw_data("XE",2.16), &
         vdw_data("CS",2.49), &
         vdw_data("BA",2.00), &
         vdw_data("LA",2.79), &
         vdw_data("CE",2.73), &
         vdw_data("PR",2.72), &
         vdw_data("ND",2.70)  &
       /

  data(r_vdw(i), i=61,80)     &
       /                      &
         vdw_data("PM",2.69), &
         vdw_data("SM",2.69), &
         vdw_data("EU",2.97), &
         vdw_data("GD",2.67), &
         vdw_data("TB",2.63), &
         vdw_data("DY",2.61), &
         vdw_data("HO",2.60), &
         vdw_data("ER",2.58), &
         vdw_data("TM",2.57), &
         vdw_data("YB",2.90), &
         vdw_data("LU",2.57), &
         vdw_data("HF",2.34), &
         vdw_data("TA",2.13), &
         vdw_data("W ",2.04), &
         vdw_data("RE",2.01), &
         vdw_data("OS",2.04), &
         vdw_data("IR",1.97), &
         vdw_data("PT",1.75), &
         vdw_data("AU",1.66), &
         vdw_data("HG",1.55)  &
       /

  data(r_vdw(i), i=81,98)     &
       /                      &
         vdw_data("TL",1.96), &
         vdw_data("PB",2.02), &
         vdw_data("BI",2.03), &
         vdw_data("PO",2.51), &
         vdw_data("AT",1.45), &
         vdw_data("RN",1.90), &
         vdw_data("FR",1.90), &
         vdw_data("RA",2.84), &
         vdw_data("AC",2.81), &
         vdw_data("TH",2.67), &
         vdw_data("PA",2.40), &
         vdw_data("U ",1.86), &
         vdw_data("NP",2.31), &
         vdw_data("PU",2.28), &
         vdw_data("AM",2.25), &
         vdw_data("CM",0.00), &
         vdw_data("BK",0.00), &
         vdw_data("CF",0.00)  &
       /
!===============================================

  R_sol=0.0d0; n_waters=1000; increment=1.0d0; solv='TIP'; w_atom='A'
  pos='C'
  n_param=0

  ! Read in input file (modified XYZ file)
  read(*,'(a100)',iostat=status) buffer
  if(status /=0) stop 'Wrong input line 1'

  call upcase(buffer)
  if(check_string(buffer,'R=')) then
     n_param=n_param+1
  else
     stop 'Radius of solvent shell has to be defined'
  end if
  if(check_string(buffer,'N=')) n_param=n_param+1
  if(check_string(buffer,'I=')) n_param=n_param+1
  if(check_string(buffer,'S=')) n_param=n_param+1
  if(check_string(buffer,'W=')) n_param=n_param+1
  if(check_string(buffer,'L=')) n_param=n_param+1
  if(n_param > max_param) stop 'Wrong input line 1'
  
  if(n_param > 1) then 
     pos1=find_character(buffer,';',n_param-1)

     if(check_string(buffer,'R=')) then
        pos2=find_character(buffer,'R',1)
        i_begin=pos2(1)
        do i=1,n_param-1
           if(i==1 .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
           if(i==n_param-1 .and. i_begin > pos1(i)) then
              i_end=len(trim(buffer))
              exit
           end if
           if(i_begin > pos1(i-1) .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
        end do
        read(buffer(i_begin:i_end),'(2x,f5.2)') R_sol
     end if

     if(check_string(buffer,'N=')) then
        pos2=find_character(buffer,'N',1)
        i_begin=pos2(1)
        do i=1,n_param-1
           if(i==1 .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
           if(i==n_param-1 .and. i_begin > pos1(i)) then
              i_end=len(trim(buffer))
              exit
           end if
           if(i_begin > pos1(i-1) .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
        end do
        read(buffer(i_begin:i_end),'(2x,i4)') n_waters
     end if

     if(check_string(buffer,'I=')) then
        pos2=find_character(buffer,'I',1)
        i_begin=pos2(1)
        do i=1,n_param-1
           if(i==1 .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
           if(i==n_param-1 .and. i_begin > pos1(i)) then
              i_end=len(trim(buffer))
              exit
           end if
           if(i_begin > pos1(i-1) .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
        end do
        read(buffer(i_begin:i_end),'(2x,f4.2)') increment
     end if

     if(check_string(buffer,'S=')) then
        pos2=find_character(buffer,'S',1)
        i_begin=pos2(1)
        do i=1,n_param-1
           if(i==1 .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
           if(i==n_param-1 .and. i_begin > pos1(i)) then
              i_end=len(trim(buffer))
              exit
           end if
           if(i_begin > pos1(i-1) .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
        end do
        read(buffer(i_begin:i_end),'(2x,a3)') solv
        if(trim(solv) /= "TIP" .and. trim(solv) /= "EFP") &
             stop 'Wrong solvent name (neither TIP nor EFP)'
     end if

     if(check_string(buffer,'W=')) then
        pos2=find_character(buffer,'W',1)
        i_begin=pos2(1)
        do i=1,n_param-1
           if(i==1 .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
           if(i==n_param-1 .and. i_begin > pos1(i)) then
              i_end=len(trim(buffer))
              exit
           end if
           if(i_begin > pos1(i-1) .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
        end do
        read(buffer(i_begin:i_end),'(2x,a1)') w_atom
        if(w_atom /= 'A' .and. w_atom /= 'M') &
             stop 'Wrong excluded water parameter (neither A nor M)'
     end if

     if(check_string(buffer,'L=')) then
        pos2=find_character(buffer,'L',1)
        i_begin=pos2(1)
        do i=1,n_param-1
           if(i==1 .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
           if(i==n_param-1 .and. i_begin > pos1(i)) then
              i_end=len(trim(buffer))
              exit
           end if
           if(i_begin > pos1(i-1) .and. i_begin < pos1(i)) then
              i_end=pos1(i)-1
              exit
           end if
        end do
        read(buffer(i_begin:i_end),'(2x,a1)') pos
        if(pos /= 'C' .and. pos /= 'D') &
             stop 'Wrong position parameter of QM cluster (neither C nor D)'
     end if
  else
     if(check_string(buffer,'R=')) then
        pos2=find_character(buffer,'R',1)
        i_begin=pos2(1)
        i_end=len(trim(buffer))
        read(buffer(i_begin:i_end),'(2x,f5.2)') R_sol
     end if
  end if

  print*,'Shell raduis: ',R_sol
  print*,'Maximal number of EFP water molecules in shell: ',n_waters
!!$  if(increment >= 1.0d0) &
!!$       print*,'Van der Waals radii of the solute atoms will be increased by ',increment
!!$  if(increment < 1.0d0) &
!!$       print*,'Van der Waals radii of the solute atoms will be reduced by ',increment
  print*,'Solvent type: ',solv
  if(w_atom == 'A') print*,'Remove water if even one atom overlaps  solute'
  if(w_atom == 'M') print*,'Remove water if even all atoms overlaps  solute'
  if(pos == 'C') then 
     center=0.0d0
     print*, 'Solute position - center of box: ',center
  else if(pos == 'D') then 
     call calc_solute_position(center)
     print*, 'Solute position - random: ',center
  end if

  read(*,*,iostat=status) n_solute_atoms
  if(status /=0) stop 'Wrong input line 2'
  print*,'Number of atoms within QM cluster: ',n_solute_atoms
  if(n_solute_atoms > max_solute_atoms) & 
       stop 'Number of atoms in input file exceeds maximal value (100)'

  g_solu_center=0.0d0  
  do i=1,n_solute_atoms
     read(*,*,iostat=status) solute(i)%name, solute(i)%coor
     if(status /=0) then
        write(buf,'(i4)') i+2
        stop 'Wrong input line'//buf
     end if
     call upcase(solute(i)%name)
     g_solu_center=g_solu_center+solute(i)%coor/n_solute_atoms
  end do

  !Move the solute cluster to the center of coordinates or random position
  do i=1,n_solute_atoms
     solute(i)%coor=solute(i)%coor-g_solu_center+center
     solute(i)%r_vdW=def_vdW_radius(solute(i)%name)!*increment
  end do

  if(gen_super_box()) then 
     n_efp_full=27*n_efp
  else
     n_efp_full=n_efp
  end if
  allocate(water_array(n_efp_full),water_shell(n_efp_full),water_center(3,n_efp_full))

  if(solv=='TIP') then
     water_array(1:n_efp) = water_array0
  else if(solv=='EFP') then
     water_array(1:n_efp) = water_array1
  end if

  if(gen_super_box()) call fill_super_box()

  call calc_water_shell(n_waters_in_shell)

  call sort_water_shell(n_waters_in_shell)

  call save_result(n_waters_in_shell)

!===============================================
  
  
contains
!===============================================

  function gen_super_box() result(yes)

    logical :: yes
    real*8 :: drxyz(2,3),half_bs

    yes=.false.

    half_bs=box_size/2.0d0

    drxyz(1,:)=half_bs-center(:); drxyz(2,:)=center(:)-(-half_bs)

    if(R_sol > drxyz(1,1) .or. &
       R_sol > drxyz(2,1) .or. &
       R_sol > drxyz(1,2) .or. &
       R_sol > drxyz(2,2) .or. &
       R_sol > drxyz(1,3) .or. &
       R_sol > drxyz(2,3)) yes = .true. 
    
  end function gen_super_box
!===============================================

  subroutine fill_super_box()

    integer*4 :: i,j,k,l,m

    m=n_efp
    do i=-1,1
       do j=-1,1
          do k=-1,1
             if(i==0 .and. j==0 .and. k==0) cycle

             do l=1,n_efp
                m=m+1
                water_array(m)%atom1%name   =water_array(l)%atom1%name
                water_array(m)%atom1%r_vdw  =water_array(l)%atom1%r_vdw
                water_array(m)%atom1%coor(1)=water_array(l)%atom1%coor(1)+real(i*box_size)
                water_array(m)%atom1%coor(2)=water_array(l)%atom1%coor(2)+real(j*box_size)
                water_array(m)%atom1%coor(3)=water_array(l)%atom1%coor(3)+real(k*box_size)

                water_array(m)%atom2%name   =water_array(l)%atom2%name
                water_array(m)%atom2%r_vdw  =water_array(l)%atom2%r_vdw
                water_array(m)%atom2%coor(1)=water_array(l)%atom2%coor(1)+real(i*box_size)
                water_array(m)%atom2%coor(2)=water_array(l)%atom2%coor(2)+real(j*box_size)
                water_array(m)%atom2%coor(3)=water_array(l)%atom2%coor(3)+real(k*box_size)

                water_array(m)%atom3%name   =water_array(l)%atom3%name
                water_array(m)%atom3%r_vdw  =water_array(l)%atom3%r_vdw
                water_array(m)%atom3%coor(1)=water_array(l)%atom3%coor(1)+real(i*box_size)
                water_array(m)%atom3%coor(2)=water_array(l)%atom3%coor(2)+real(j*box_size)
                water_array(m)%atom3%coor(3)=water_array(l)%atom3%coor(3)+real(k*box_size)
             end do

          end do
       end do
    end do

  end subroutine fill_super_box
!===============================================

  subroutine calc_solute_position(position)

    real*8 ::  position(3)
    real*8 :: harvest(3)
    integer*4 :: values(8)
    integer*4 :: seed(1)

    call date_and_time(values=values)

    seed(1)=values(6)*60.0d0*1000.0d0+values(7)*1000.0d0+values(8)

    call random_seed(put=seed)
    call random_number(harvest)

    position=(harvest-0.5d0)*box_size

  end subroutine calc_solute_position
!===============================================

  subroutine upcase(string)

    character(len=*) :: string
    integer*4 :: i,ln,ich

    ln = len(trim(string))
    do i=1,ln
       ich=iachar(string(i:i))
       if(ich>=97 .and. ich<=122) then
          ich=ich-32
          string(i:i)=achar(ich)
       end if
    end do

  end subroutine upcase
!===============================================

  function check_string(string,word)

    character(len=*), intent(in) :: string
    character(len=*), intent(in) :: word

    logical :: check_string

    if(index(string,word) /= 0) then
       check_string = .true.
    else
       check_string = .false.
    end if

  end function check_string
!===============================================

  function find_character(string,char,n_char) result(chars)

    integer*4 :: chars(n_char)
    integer*4 , intent(in) :: n_char
    character(len=*), intent(in) :: string
    character(len=1), intent(in) :: char

    integer*4 :: i,j

    j=1
    do i=1,n_char
       chars(i)=index(string(j:),char)+j-1
       j=chars(i)+1
    end do

  end function find_character
!===============================================

  function def_vdW_radius(name) result(radius)

    character*2,intent(in) :: name
    real*8 :: radius

    integer*4 :: i

    radius=0.0d0

    do i=1,n_atoms
       if(trim(name)==trim(r_vdW(i)%name)) radius=r_vdW(i)%radius
    end do

  end function def_vdW_radius
!===============================================

  subroutine calc_water_shell(n_h2o)

    integer*4 :: n_h2o

    integer*4 :: i,k,l
    real*8 :: r1(3),r2(3),r3(3),dr1,dr2,dr3,r_water(3,3)
    real*8 :: r_vdw_water(3),vdw_sum1,vdw_sum2,vdw_sum3
    real*8 :: w_center(3)
    integer*4 :: nn

    call def_efp_water()

    l=0
    first_cycle:do i=1,n_efp_full
       call produce_EFP_water(i,w_center)

       r_water(:,1)=water_array(i)%atom1%coor
       r_water(:,2)=water_array(i)%atom2%coor
       r_water(:,3)=water_array(i)%atom3%coor

       r_vdw_water(1)=water_array(i)%atom1%r_vdw
       r_vdw_water(2)=water_array(i)%atom2%r_vdw
       r_vdw_water(3)=water_array(i)%atom3%r_vdw

       nn=0
       third_cycle:do k=1,n_solute_atoms

          r1=r_water(:,1)-solute(k)%coor; dr1=sqrt(dot_product(r1,r1))
          r2=r_water(:,2)-solute(k)%coor; dr2=sqrt(dot_product(r2,r2))
          r3=r_water(:,3)-solute(k)%coor; dr3=sqrt(dot_product(r3,r3))
          if(dr1 > R_sol .and. dr2 > R_sol .and. dr3 > R_sol) then 
             nn=nn+1
             cycle third_cycle
          end if

          vdw_sum1=r_vdw_water(1)+solute(k)%r_vdw
          vdw_sum2=r_vdw_water(2)+solute(k)%r_vdw
          vdw_sum3=r_vdw_water(3)+solute(k)%r_vdw

          if(w_atom == 'A') then
             if(dr1 < vdw_sum1 .or. dr2 < vdw_sum2 .or. dr3 < vdw_sum3) cycle first_cycle
          else if(w_atom == 'M') then
             if(dr1 < vdw_sum1 .and. dr2 < vdw_sum2 .and. dr3 < vdw_sum3) cycle first_cycle
          end if
       end do third_cycle
       if(nn==n_solute_atoms) cycle first_cycle
       l=l+1
       water_shell(l)=water_array(i)
       water_center(:,l)=w_center
    end do first_cycle

    n_h2o=l
    print*,'Calculated number of EFP water molecules within shell: ',n_h2o

  end subroutine calc_water_shell
!===============================================

  subroutine def_efp_water()

    EFP_O1=(/ 0.0d0     , 0.0d0, -0.119151d0/)
    EFP_H2=(/-1.431042d0, 0.0d0,  0.945510d0/)
    EFP_H3=(/ 1.431042d0, 0.0d0,  0.945510d0/)

  end subroutine def_efp_water
!===============================================

  function get_efp_OH() result(length)

    real*8 :: length
    real*8 :: r1(3),r2(3),r3(3)

    r1=EFP_O1; r2=EFP_H2; r3=EFP_H3
    length=sqrt(dot_product(r2-r1,r2-r1))
    
  end function get_efp_OH
!===============================================

  function get_efp_HOH() result(angle)

    real*8 :: angle
    real*8 :: r1(3),r2(3),r3(3),d12,d13,cos_a

    r1=EFP_O1; r2=EFP_H2; r3=EFP_H3

    d12=sqrt(dot_product(r2-r1,r2-r1))
    d13=sqrt(dot_product(r3-r1,r3-r1))

    cos_a=dot_product(r2-r1,r3-r1)/(d12*d13)
    angle=acos(cos_a)
    
  end function get_efp_HOH
!===============================================

  subroutine add_1_4_atom(ra,rb,internals)

    real*8 :: ra(3),rb(3,3)
    real*8 :: internals(3)

    real*8 :: bond,angle,dihedral
    real*8 :: cos_a,sin_a,cos_d,sin_d
    real*8 :: r21(3),d21,r32(3),d32,r1(3),d1,r2(3),d2,r3(3),r4(3),d4

    bond=internals(1)
    angle=internals(2)
    dihedral=internals(3)

    cos_a=cos(angle); sin_a=sin(angle)
    cos_d=cos(dihedral); sin_d=sin(dihedral)

    r21=rb(:,2)-rb(:,1); r32=rb(:,3)-rb(:,2)
    d21=sqrt(dot_product(r21,r21))
    r21=r21/d21
    d32=sqrt(dot_product(r32,r32))
    r32=r32/d32

    r1=vector_product(r21,r32)
    d1=sqrt(dot_product(r1,r1))
    r1=r1/d1

    r2=vector_product(r1,r32)
    d2=sqrt(dot_product(r2,r2))
    r2=r2/d2

    r3=r2*cos_d+r1*sin_d

    r4=-r32*cos_a+r3*sin_a
    d4=sqrt(dot_product(r4,r4))
    r4=r4*bond/d4

    ra=rb(:,3)+r4

  end subroutine add_1_4_atom
!===============================================

  function vector_product(v1,v2)

    real*8 :: vector_product(3)
    real*8 :: v1(3),v2(3)

    vector_product(1)=v1(2)*v2(3)-v1(3)*v2(2)
    vector_product(2)=v1(3)*v2(1)-v1(1)*v2(3)
    vector_product(3)=v1(1)*v2(2)-v1(2)*v2(1)

  end function vector_product
!===============================================

  subroutine produce_efp_water(num,ccenter)

    integer*4 :: num
    real*8 :: ccenter(3)

    real*8 :: OH,HOH
    real*8 :: r(3,3),r12(3),r13(3),d12,d13,cos_a,angle
    real*8 :: rm(3,3),zmat(3),rq(3),r1q(3),d1q
    real*8, parameter :: pi=3.1415926535897932368d0

    OH=get_efp_OH()*0.529177249d0 ; HOH=get_efp_HOH()

    r(:,1)=water_array(num)%atom1%coor
    r(:,2)=water_array(num)%atom2%coor
    r(:,3)=water_array(num)%atom3%coor

    r12=r(:,2)-r(:,1); d12=sqrt(dot_product(r12,r12))
    r13=r(:,3)-r(:,1); d13=sqrt(dot_product(r13,r13))

    !bond correction
    if(d12 /= OH) then
       r12=r12*OH/d12; r(:,2)=r(:,1)+r12
    end if
    if(d13 /= OH) then
       r13=r13*OH/d13; r(:,3)=r(:,1)+r13
    end if

    !angle correction
    cos_a=dot_product(r12,r13)/(d12*d13); angle=acos(cos_a)
    if(angle /= HOH) then
       rq=(r(:,2)+r(:,3))/2.0d0
       r1q=rq-r(:,1); d1q=sqrt(dot_product(r1q,r1q))
       
       rm(:,1)=r(:,3); rm(:,2)=rq; rm(:,3)=r(:,1)
       zmat(1)=OH; zmat(2)=HOH/2.0d0; zmat(3)=pi
       call add_1_4_atom(r(:,2),rm,zmat)

       rm(:,1)=r(:,2); rm(:,2)=rq; rm(:,3)=r(:,1)
       zmat(1)=OH; zmat(2)=HOH/2.0d0; zmat(3)=pi
       call add_1_4_atom(r(:,3),rm,zmat)
    end if

    water_array(num)%atom1%coor=r(:,1)
    water_array(num)%atom2%coor=r(:,2)
    water_array(num)%atom3%coor=r(:,3)

    ccenter=(r(:,1)+r(:,2)+r(:,3))/3.0d0

    r12=r(:,2)-r(:,1); d12=sqrt(dot_product(r12,r12))
    r13=r(:,3)-r(:,1); d13=sqrt(dot_product(r13,r13))
    cos_a=dot_product(r12,r13)/(d12*d13); angle=acos(cos_a)

  end subroutine produce_efp_water
!===============================================

  subroutine sort_water_shell(n_h2o)

    integer*4 :: n_h2o

    integer*4 :: i,j,k
    real*8 :: dr1,dr2,ccenter(3)
    type(efp_water) :: water_buf

    first_cycle:do i=2,n_h2o
       dr1=sqrt(dot_product(water_center(:,i)-center,water_center(:,i)-center))

       second_cycle:do j=1,i-1
          dr2=sqrt(dot_product(water_center(:,j)-center,water_center(:,j)-center))

          if(dr1 >= dr2) cycle second_cycle

          water_buf=water_shell(i)
          ccenter=water_center(:,i)

          third_cycle:do k=i-1,j,-1
             water_shell(k+1)=water_shell(k)
             water_center(:,k+1)=water_center(:,k)
          end do third_cycle

          water_shell(j)=water_buf
          water_center(:,j)=ccenter
          
          cycle first_cycle

       end do second_cycle
    end do first_cycle

  end subroutine sort_water_shell
!===============================================

  subroutine save_result(n_h2o)

    integer*4 :: n_h2o

    integer*4, parameter :: id_gx=10, id_inp=11, id_xyz=12
    integer*4 :: n_xyz
    integer*4 :: i,j
    real*8 ::at_num

    if(n_h2o > n_waters) n_h2o=n_waters
    print*,'Saved number of EFP water molecules within shell: ',n_h2o

    n_xyz=n_solute_atoms+n_h2o*3

    !save XYZ
    open(id_xyz,file="QM_shell.xyz")

    write(id_xyz,'(i6)') n_xyz
    write(id_xyz, '(a12)') '!!!!!!!!!!!!'
    do i=1,n_solute_atoms
       write(id_xyz,'(a2,4x,3f15.8)') solute(i)%name, solute(i)%coor
    end do
    do i=1,n_h2o
       write(id_xyz,'(a1,5x,3f15.8)') & 
            water_shell(i)%atom1%name, water_shell(i)%atom1%coor
       write(id_xyz,'(a1,5x,3f15.8)') & 
            water_shell(i)%atom2%name, water_shell(i)%atom2%coor
       write(id_xyz,'(a1,5x,3f15.8)') & 
            water_shell(i)%atom3%name, water_shell(i)%atom3%coor       
    end do

    close(id_xyz)

    !save GX
    open(id_gx,file="gx.QM_shell")

    do i=1,n_solute_atoms
       at_num=get_atom_number(solute(i)%name)
       write(id_gx,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4)') &
            at_num,solute(i)%coor/0.529177249d0,i,i,0,0,0,0,0,0
    end do
    j=i-1
    do i=1,n_h2o
       j=j+1
       write(id_gx,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4)') &
            8.01, water_shell(i)%atom1%coor/0.529177249d0,j,j,0,0,0,0,0,0
       write(id_gx,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4)') &
            1.02, water_shell(i)%atom2%coor/0.529177249d0,j,j,0,0,0,0,0,0
       write(id_gx,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4)') &
            1.03, water_shell(i)%atom3%coor/0.529177249d0,j,j,0,0,0,0,0,0
    end do
    write(id_gx,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4,2x,i4)') &
         -1.0d0,0,0d0,0.0d0,0.0d0,0,0,0,0,0,0,0,0,0

    close(id_gx)

    !save input
    open(id_inp,file="QM_shell.inp")

    write(id_inp,'(1x,a1)')                 ' '
    write(id_inp,'(1x,a15)')                '&SYMMETRY_GROUP'
    write(id_inp,'(1x,a23)')                '   POINT_GROUP = "C1  "'
    write(id_inp,'(1x,a15)')                '/SYMMETRY_GROUP'
    write(id_inp,'(1x,a1)')                 ' '
    write(id_inp,'(1x,a19)')                '&UNIQUE_ATOM_NUMBER'
    write(id_inp,'(1x,a19,i4)')             '   N_UNIQUE_ATOMS =',n_solute_atoms  
    write(id_inp,'(1x,a19)')                '/UNIQUE_ATOM_NUMBER'
    write(id_inp,'(1x,a1)')                 ' '
    do i=1,n_solute_atoms
       at_num=get_atom_number(solute(i)%name)
       write(id_inp,'(1x,a12)')             '&UNIQUE_ATOM'
       write(id_inp,'(1x,a20,a2,a1)')       '   NAME          = "',solute(i)%name,'"'
       write(id_inp,'(1x,a20,f6.3)')        '   Z             =  ',at_num
       write(id_inp,'(1x,a20,i4)')          '   N_EQUAL_ATOMS =  ',1
       write(id_inp,'(1x,a12)')             '/UNIQUE_ATOM'
       write(id_inp,'(3f21.15)')             solute(i)%coor/0.529177249d0
       write(id_inp,'(1x,a1)')               ' '
    end do

    write(id_inp,'(1x,a43)')                 '#Do not forget to add basis sets !!!!!!!!!!'
    write(id_inp,'(1x,a1)')                  ' '

    write(id_inp,'(1x,a9)')                  '&EFP_DATA'
    write(id_inp,'(1x,a17,i4)')              '   N_EFP        =',n_h2o
    write(id_inp,'(1x,a22)')                 '   FIXED_REGION = "NO"'
    write(id_inp,'(1x,a9)')                  '/EFP_DATA'
    write(id_inp,'(1x,a1)')                  ' '
    do i=1,n_h2o
       write(id_inp,'(1x,a10)')              '&EFP_WATER'
       write(id_inp,'(1x,a12,a2,a1)')        '   NAME1 = "',water_shell(i)%atom1%name,'"'
       write(id_inp,'(1x,a12,3(f15.10,a1))') '   COOR1 =  ',water_shell(i)%atom1%coor(1)/0.529177249d0,',', &
                                                            water_shell(i)%atom1%coor(2)/0.529177249d0,',', &
                                                            water_shell(i)%atom1%coor(3)/0.529177249d0,','
       write(id_inp,'(1x,a12,a2,a1)')        '   NAME2 = "',water_shell(i)%atom2%name,'"'
       write(id_inp,'(1x,a12,3(f15.10,a1))') '   COOR2 =  ',water_shell(i)%atom2%coor(1)/0.529177249d0,',', &
                                                            water_shell(i)%atom2%coor(2)/0.529177249d0,',', &
                                                            water_shell(i)%atom2%coor(3)/0.529177249d0,','
       write(id_inp,'(1x,a12,a2,a1)')        '   NAME3 = "',water_shell(i)%atom3%name,'"'
       write(id_inp,'(1x,a12,3(f15.10,a1))') '   COOR3 =  ',water_shell(i)%atom3%coor(1)/0.529177249d0,',', &
                                                            water_shell(i)%atom3%coor(2)/0.529177249d0,',', &
                                                            water_shell(i)%atom3%coor(3)/0.529177249d0,','
       write(id_inp,'(1x,a10)')              '/EFP_WATER'
       write(id_inp,'(1x,a1)')               ' '
    end do

    close(id_inp)

  end subroutine save_result
!===============================================

  function get_atom_number(name) result(atom_number)

    real*8 :: atom_number
    character*2 :: name

    integer*4 :: i

    do i=1,n_atoms
       if(trim(name)==trim(r_vdW(i)%name)) then
          atom_number=i
          exit
       end if
    end do

  end function get_atom_number
!===============================================
end program solute_in_water

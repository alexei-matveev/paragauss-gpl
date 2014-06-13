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
module element_data_module
  !------------ Modules used --------------------------------------
  use type_module

  implicit none
  private       
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of public constants and variables -----
  type, public :: eledata
     character(len=2) :: name
     real(kind=r8_kind) :: r_coval
     real(kind=r8_kind) :: r_vdw_bondi
     real(kind=r8_kind) :: r0_dr
     real(kind=r8_kind) :: k_dr
     real(kind=r8_kind) :: solv_scale_fac
     real(kind=r8_kind) :: solv_scale_fac_smooth
  end type eledata
  type(eledata), public :: atom_data(98)=(/                                                   &
                                           eledata("H_",0.37, 1.20, 1.20, 1.00, 1.00, 0.9375), &
                                           eledata("HE",1.60, 1.40, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("LI",0.68, 1.82, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("BE",0.35, 0.52, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("B_",0.83, 1.70, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("C_",0.77, 1.70, 1.70, 1.00, 1.20, 1.125), &
                                           eledata("N_",0.75, 1.55, 1.60, 1.18, 1.20, 1.125), &
                                           eledata("O_",0.73, 1.52, 1.50, 1.36, 1.20, 1.125), &
                                           eledata("F_",0.71, 1.47, 1.45, 1.50, 1.20, 1.125), &
                                           eledata("NE",1.12, 1.54, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("NA",0.97, 2.27, 1.20, 1.40, 1.20, 1.125), &
                                           eledata("MG",1.10, 1.73, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("AL",1.35, 2.01, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("SI",1.20, 2.10, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("P_",1.05, 1.80, 1.85, 2.10, 1.20, 1.125), &
                                           eledata("S_",1.02, 1.80, 1.80, 2.40, 1.20, 1.125), &
                                           eledata("CL",0.99, 1.75, 1.76, 2.10, 1.20, 1.125), &
                                           eledata("AR",1.54, 1.88, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("K_",1.33, 2.75, 1.46, 2.90, 1.20, 1.125), &
                                           !                   *
                                           eledata("CA",0.99, 1.48, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("SC",1.44, 2.15, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("TI",1.47, 2.19, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("V_",1.33, 1.99, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("CR",1.35, 2.01, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("MN",1.35, 2.01, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("FE",1.34, 2.00, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("CO",1.33, 1.99, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("NI",1.50, 1.63, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("CU",1.52, 1.40, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("ZN",1.45, 1.39, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("GA",1.22, 1.87, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("GE",1.17, 1.75, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("AS",1.21, 1.85, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("SE",1.22, 1.90, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("BR",1.21, 1.85, 1.85, 2.40, 1.20, 1.125), &
                                           eledata("KR",1.89, 2.06, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("RB",1.47, 2.19, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("SR",1.12, 1.67, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("Y_",1.78, 2.66, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("ZR",1.56, 2.33, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("NB",1.48, 2.21, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("MO",1.47, 2.19, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("TC",1.35, 2.01, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("RU",1.40, 2.09, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("RH",1.45, 2.16, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("PD",1.54, 1.63, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("AG",1.59, 1.72, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("CD",1.69, 1.58, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("IN",1.63, 1.93, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("SN",1.46, 2.17, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("SB",1.46, 2.17, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("TE",1.47, 2.06, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("I_",1.40, 1.98, 1.96, 3.20, 1.20, 1.125), &
                                           eledata("XE",1.33, 2.16, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("CS",1.67, 2.49, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("BA",1.34, 2.00, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("LA",1.87, 2.79, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("CE",1.83, 2.73, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("PR",1.82, 2.72, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("ND",1.81, 2.70, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("PM",1.80, 2.69, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("SM",1.80, 2.69, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("EU",1.99, 2.97, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("GD",1.79, 2.67, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("TB",1.76, 2.63, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("DY",1.75, 2.61, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("HO",1.74, 2.60, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("ER",1.73, 2.58, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("TM",1.72, 2.57, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("YB",1.94, 2.90, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("LU",1.72, 2.57, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("HF",1.57, 2.34, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("TA",1.43, 2.13, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("W_",1.37, 2.04, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("RE",1.35, 2.01, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("OS",1.37, 2.04, 0.00, 0.00, 1.20, 1.125), & 
                                           !                   *
                                           eledata("IR",1.32, 1.97, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("PT",1.50, 1.75, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("AU",1.50, 1.66, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("HG",1.70, 1.55, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("TL",1.55, 1.96, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("PB",1.54, 2.02, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("BI",1.54, 2.03, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("PO",1.68, 2.51, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("AT",1.45, 0.00, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("RN",1.90, 0.00, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("FR",0.00, 0.00, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("RA",1.90, 2.84, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("AC",1.88, 2.81, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("TH",1.79, 2.67, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("PA",1.61, 2.40, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("U_",1.58, 1.86, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("NP",1.55, 2.31, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("PU",1.53, 2.28, 0.00, 0.00, 1.20, 1.125), &
                                           !                   *
                                           eledata("AM",1.51, 2.25, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("CM",0.00, 0.00, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("BK",0.00, 0.00, 0.00, 0.00, 1.20, 1.125), &
                                           eledata("CF",0.00, 0.00, 0.00, 0.00, 1.20, 1.125)  /)

  !------------ public functions and subroutines ------------------
  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of private constants and variables ----
  !------------ Subroutines -----------------------------------------
end module element_data_module

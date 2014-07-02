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
module atoms_data_module
!
! Author: AS
! Date: 03 2002
!
! This module collect set of atomic data which are used 
! at the moment in solvation_module, disp_rep_module and
! potential_calc_module 
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!=====================================================================
# include "def.h"
  use type_module ! type specification parameters
  implicit none

  save            ! save all variables defined in this module
  private         ! by default, all names are private
!== Interrupt end of public interface of module =================
!
!-- Public variables --------------------------------------------
!!!  character(len=2),parameter,public :: &
#ifdef _LINUX
  character(len=2),parameter,public :: &
# else
  character(len=2),public :: &
#endif   
       atom_name(98)= (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne", &
                        "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca", &
                        "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
                        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr", &
                        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn", &
                        "Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                        "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tu","Yb", &
                        "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
                        "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", &
                        "Pa","U ","Np","Pu","Am","Cm","Bk","Cf"/)

  ! default van der Waals radii used in solvation_module to build a solvent 
  ! excluding molecular surfaces (A.Bondi, J.Phys.Chem., V.68, N.3, P.441, 1964)
!!!  real(kind=r8_kind),parameter,public :: &
#ifdef _LINUX
  real(kind=r8_kind),parameter,public :: &
# else
  real(kind=r8_kind),public :: &
#endif
       vdW_radius(98)=(/1.20,1.40,1.82,0.00,0.00,1.70,1.55,1.52,1.47,1.54, &
                        ! H    He   Li   Be   B    C    N    O    F    Ne 
                        2.27,1.73,0.00,2.10,1.80,1.80,1.75,1.88,2.75,0.00, &
                        ! Na   Mg   Al   Si   P    S    Cl   Ar   K    Ca
                        0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.63,1.40,1.39, &
                        ! Sc   Ti   V    Cr   Mn   Fe   Co   Ni   Cu   Zn
                        1.87,0.00,1.85,1.90,1.85,2.06,0.00,0.00,0.00,0.00, &
                        ! Ga   Ge   As   Se   Br   Kr   Rb   Sr   Y    Zr
                        0.00,0.00,0.00,0.00,0.00,1.63,1.72,1.58,1.93,2.17, &
                        ! Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn
                        0.00,2.06,1.98,2.16,0.00,0.00,0.00,0.00,0.00,0.00, &
                        ! Sb   Te   I    Xe   Cs   Ba   La   Ce   Pr   Nd
                        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                        ! Pm   Sm   Eu   Gd   Tb   Dy   Ho   Er   Tu   Yb
                        0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.75,1.66,1.55, &
                        ! Lu   Hf   Ta   W    Re   Os   Ir   Pt   Au   Hg
                        1.96,2.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                        ! Tl   Pb   Bi   Po   At   Rn   Fr   Ra   Ac   Th
                        0.00,1.86,0.00,0.00,0.00,0.00,0.00,0.00/)
                        ! Pa   U    Np   Pu   Am   Cm   Bk   Cf

  ! default van der Waals radii and K parameters used in disp_rep_module 
  ! to calculate the short range contribution to the free solvation energy
  ! (J.Caillet, P.Claverie, B.Pullman, Acta Cryst. B, V.32, P.2740, 1976
  !  J.Caillet, P.Claverie, Acta Cryst. A, V.31, P.448, 1975
  !  F.Vign�-Maeder, P.Claverie, J.Am.Chem.Soc., V.109, P.24, 1987
  !  J.Caillet, P.Claverie, B.Pullman, Acta Cryst. B, V.36, P.3266, 1978)
!!!  real(kind=r8_kind),parameter,public :: &
#ifdef _LINUX
  real(kind=r8_kind),parameter,public :: &
# else
  real(kind=r8_kind),public :: &
#endif
       R0_def(98)=(/1.20,0.00,0.00,0.00,0.00,1.70,1.60,1.50,1.45,0.00, &
                    ! H    He   Li   Be   B    C    N    O    F    Ne 
                    1.20,0.00,0.00,0.00,1.85,1.80,1.76,0.00,1.46,0.00, &
                    ! Na   Mg   Al   Si   P    S    Cl   Ar   K    Ca
                    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                    ! Sc   Ti   V    Cr   Mn   Fe   Co   Ni   Cu   Zn
                    0.00,0.00,0.00,0.00,1.85,0.00,0.00,0.00,0.00,0.00, &
                    ! Ga   Ge   As   Se   Br   Kr   Rb   Sr   Y    Zr
                    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                    ! Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn
                    0.00,0.00,1.96,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                    ! Sb   Te   I    Xe   Cs   Ba   La   Ce   Pr   Nd
                    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                    ! Pm   Sm   Eu   Gd   Tb   Dy   Ho   Er   Tu   Yb
                    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                    ! Lu   Hf   Ta   W    Re   Os   Ir   Pt   Au   Hg
                    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                    ! Tl   Pb   Bi   Po   At   Rn   Fr   Ra   Ac   Th
                    0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
                    ! Pa   U    Np   Pu   Am   Cm   Bk   Cf

!!!  real(kind=r8_kind),parameter,public :: & 
#ifdef _LINUX
  real(kind=r8_kind),parameter,public :: &
# else
  real(kind=r8_kind),public :: &
#endif
       K_def(98)=(/1.00,0.00,0.00,0.00,0.00,1.00,1.18,1.36,1.50,0.00, &
                   ! H    He   Li   Be   B    C    N    O    F    Ne 
                   1.40,0.00,0.00,0.00,2.10,2.40,2.10,0.00,2.90,0.00, &
                   ! Na   Mg   Al   Si   P    S    Cl   Ar   K    Ca
                   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                   ! Sc   Ti   V    Cr   Mn   Fe   Co   Ni   Cu   Zn
                   0.00,0.00,0.00,0.00,2.40,0.00,0.00,0.00,0.00,0.00, &
                   ! Ga   Ge   As   Se   Br   Kr   Rb   Sr   Y    Zr
                   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                   ! Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn
                   0.00,0.00,3.20,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                   ! Sb   Te   I    Xe   Cs   Ba   La   Ce   Pr   Nd
                   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                   ! Pm   Sm   Eu   Gd   Tb   Dy   Ho   Er   Tu   Yb
                   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                   ! Lu   Hf   Ta   W    Re   Os   Ir   Pt   Au   Hg
                   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, &
                   ! Tl   Pb   Bi   Po   At   Rn   Fr   Ra   Ac   Th
                   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00/)
                   ! Pa   U    Np   Pu   Am   Cm   Bk   Cf

  ! (A.K.Rappe, C.C. Casewit, K.S. Colwell, W.A.Goddard, W.M.Skiff, 
  ! J.Am.Chem.Soc., V.114, P.10024, 1992)
!!!  real(kind=r8_kind),parameter,public :: &
#ifdef _LINUX
  real(kind=r8_kind),parameter,public :: &
# else
  real(kind=r8_kind),public :: &
#endif
       R_def_rap(103)=(/2.886,2.362,2.451,2.745,4.083,3.851,3.660,3.500, &
                        !  H     He    Li    Be    B     C     N     O
                        3.364,3.243,2.983,3.021,4.499,4.295,4.147,4.035, &
                        !  F     Ne    Na    Mg    Al    Si    P     S
                        3.947,3.868,3.812,3.399,3.295,3.175,3.144,3.023, &
                        !  Cl    Ar    K     Ca    Sc    Ti    V     Cr
                        2.961,2.912,2.872,2.834,3.495,2.763,4.383,4.280, &
                        !  Mn    Fe    Co    Ni    Cu    Zn    Ga    Ge
                        4.230,4.205,4.189,4.141,4.114,3.641,3.345,3.124, &
                        !  As    Se    Br    Kr    Rb    Sr    Y     Zr
                        3.165,3.052,2.998,2.963,2.929,2.899,3.148,2.848, &
                        !  Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
                        4.463,4.392,4.420,4.470,4.500,4.404,4.517,3.703, &
                        !  In    Sn    Sb    Te    J     Xe    Cs    Ba
                        3.522,3.556,3.606,3.575,3.547,3.520,3.493,3.368, &
                        !  La    Ce    Pr    Nd    Pm    Sm    Eu    Gd
                        3.451,3.428,3.409,3.391,3.374,3.355,3.640,3.141, &
                        !  Tb    Dy    Ho    Er    Tm    Yb    Lu    Hf
                        3.170,3.069,2.954,3.120,2.840,2.754,3.293,2.705, &
                        !  Ta    W     Re    Os    Ir    Pt    Au    Hg
                        4.347,4.297,4.370,4.709,4.750,4.765,4.900,3.677, &
                        !  Tl    Pb    Bi    Po    At    Rn    Fr    Ra
                        3.478,3.396,3.424,3.395,3.424,3.424,3.381,3.326, &
                        !  Ac    Th    Pa    U     Np    Pu    Am    Cm
                        3.339,3.313,3.299,3.286,3.274,3.248,3.236/)
                        !  Bk    Cf    Es    Fm    Md    No    Lw
  !!!real(kind=r8_kind),parameter,public :: &
#ifdef _LINUX
  real(kind=r8_kind),parameter,public :: &
# else
  real(kind=r8_kind),public :: &
#endif
       D_def_rap(103)=(/0.044,0.056,0.025,0.085,0.180,0.105,0.069,0.060, &
                        !  H     He    Li    Be    B     C     N     O
                        0.050,0.042,0.030,0.111,0.505,0.402,0.305,0.274, &
                        !  F     Ne    Na    Mg    Al    Si    P     S
                        0.227,0.185,0.035,0.238,0.019,0.017,0.016,0.015, &
                        !  Cl    Ar    K     Ca    Sc    Ti    V     Cr
                        0.013,0.013,0.014,0.015,0.005,0.124,0.415,0.379, &
                        !  Mn    Fe    Co    Ni    Cu    Zn    Ga    Ge
                        0.309,0.291,0.251,0.220,0.040,0.235,0.072,0.069, &
                        !  As    Se    Br    Kr    Rb    Sr    Y     Zr
                        0.059,0.056,0.048,0.056,0.053,0.048,0.036,0.228, &
                        !  Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
                        0.599,0.567,0.449,0.398,0.339,0.332,0.045,0.364, &
                        !  In    Sn    Sb    Te    J     Xe    Cs    Ba
                        0.017,0.013,0.010,0.010,0.009,0.008,0.008,0.009, &
                        !  La    Ce    Pr    Nd    Pm    Sm    Eu    Gd
                        0.007,0.007,0.007,0.007,0.006,0.228,0.041,0.072, &
                        !  Tb    Dy    Ho    Er    Tm    Yb    Lu    Hf
                        0.081,0.067,0.066,0.037,0.073,0.080,0.039,0.385, &
                        !  Ta    W     Re    Os    Ir    Pt    Au    Hg
                        0.680,0.663,0.518,0.325,0.284,0.248,0.050,0.404, &
                        !  Tl    Pb    Bi    Po    At    Rn    Fr    Ra
                        0.033,0.026,0.022,0.022,0.019,0.016,0.014,0.013, &
                        !  Ac    Th    Pa    U     Np    Pu    Am    Cm
                        0.013,0.013,0.012,0.012,0.011,0.011,0.011/)
                        !  Bk    Cf    Es    Fm    Md    No    Lw

!!!  real(kind=r8_kind),parameter,public :: &
#ifdef _LINUX
  real(kind=r8_kind),parameter,public :: &
# else
  real(kind=r8_kind),public :: &
#endif
       zeta_def_rap(103)=(/12.00 ,15.24 ,12.00 ,12.00 ,12.052,12.73 ,13.407,14.085, &
                           !  H      He     Li     Be     B      C      N      O
                           14.762,15.440,12.0  ,12.0  ,11.278,12.175,13.072,13.969, &
                           !  F      Ne     Na     Mg     Al     Si     P      S
                           14.866,15.763,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  Cl     Ar     K      Ca     Sc     Ti     V      Cr
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,11.0  ,12.0  , &
                           !  Mn     Fe     Co     Ni     Cu     Zn     Ga     Ge
                           13.0  ,14.0  ,15.0  ,16.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  As     Se     Br     Kr     Rb     Sr     Y      Zr
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  Nb     Mo     Tc     Ru     Rh     Pd     Ag     Cd
                           11.0  ,12.0  ,13.0  ,14.0  ,15.0  ,12.0  ,12.0  ,12.0  , &
                           !  In     Sn     Sb     Te     J      Xe     Cs     Ba
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  La     Ce     Pr     Nd     Pm     Sm     Eu     Gd
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  Tb     Dy     Ho     Er     Tm     Yb     Lu     Hf
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  Ta     W      Re     Os     Ir     Pt     Au     Hg
                           11.0  ,12.0  ,13.0  ,14.0  ,15.0  ,16.0  ,12.0  ,12.0  , &
                           !  Tl     Pb     Bi     Po     At     Rn     Fr     Ra
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  , &
                           !  Ac     Th     Pa     U      Np     Pu     Am     Cm
                           12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0  ,12.0/)
                           !  Bk     Cf     Es     Fm     Md     No     Lw

end module atoms_data_module



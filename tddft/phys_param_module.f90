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
MODULE  phys_param_module
  !---------------------------------------------------------------
  !
  !  Purpose: 
  !  Database to hold values of physical constants.
  !  
  !  
  !  Module called by: nearly every other module/subroutine
  !
  !  Author: HH
  !  Date:   12/98
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
  !
  !----------------------------------------------------------------

  USE type_module ! type specification parameters

  IMPLICIT NONE

  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of constants and variables ------------

  ! NIST reference values used throughout. 
  ! Note: (xx) = uncertainty in last 2 digits is xx
  ! 1 hartree (i.e. a.u.) = 27.2113961(81) eV
  ! 1 Angstroem^-1        = 455.63352672(54) hartree
!!$REAL(KIND=r8_kind),PARAMETER,PUBLIC  :: hartree2ev = 27.2113961_r8_kind,&
!!$     &invang2hartree = 455.63352672_r8_kind
  !
  ! Note: the values below were used in tdfrt_V* and are needed for comparison 
  REAL(KIND=r8_kind),PARAMETER,PUBLIC  :: & 
       & hartree2ev     = 27.211652000_r8_kind,&
       & ev2hartree     = 0.0367493090_r8_kind,&
       & invang2hartree = 455.63352672_r8_kind
       
END MODULE phys_param_module

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
	module	latt_module
	real*8, dimension(3,3),target	::	A
	integer ::NATOMS,NPOINT,NQT,KSD,KSR,np_ions=0
	integer, dimension(2,3):: npd=0	! limits of LCGTO array
	integer, parameter:: ngxat=200	! limit of No of gx-centers
	real*8, parameter:: angSAU=0.529177D0
	logical pcmerge		! true to make or use merge of the PC arr
	logical make_pcmerge	! true to make or use merge of the PC arr
        integer ::      numion=0
	real*8 :: c1=1.0,c2=1.0
	integer*4 nopc_arr	! number of atoms in the merged PC arr
	logical other_pc_arr    ! file with merged PC array exist
	real*8, allocatable, dimension(:,:):: rrpm	! merged PC arr
	logical:: xmol_form=.false.
	end	module latt_module

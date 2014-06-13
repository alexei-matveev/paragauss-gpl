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
	module smlcom
	integer, parameter:: r8_kind=selected_real_kind(15)
	integer, parameter :: i4_kind = selected_int_kind(9)
	integer, parameter:: kamax = 150
	integer, parameter:: nlocv = 400
	integer, parameter:: iu152 = 92
	integer, parameter:: nlocm = ((nlocv+1)*nlocv)/2
	integer, dimension(nlocv)::	mu_type
	real*8,  dimension(nlocv)::	sxo
	integer, dimension(kamax)::	iepe
	real(kind=r8_kind),parameter:: angsau=0.52917706_r8_kind
	real(kind=r8_kind) pmste
	integer(kind=i4_kind) :: igener(nlocv)
	end module smlcom

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
	module sar_module
	integer, parameter:: ndqt=20, ndvt=40 ,ndpt=400
	character*4, dimension(400)::    symvl
	integer, dimension(400)	::	IQT,MT,MQ,NA
	real*8,	dimension(3,400):: 	RATOM,RP
	real*8,	dimension(ndqt,ndpt)::	PSUM,PSUMnp
	real*8, dimension(400,ndqt)::	VM
	real(8), dimension(ndpt)::	VMS,PMADL,PSUMA,vel
     &         ,vpcm	! potenytials due to mergrd arr's
	real(8), dimension(ndvt)::	SUMA,SUMB,sump
	real	ZATOM(20)
	end	module sar_module

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
 real*8,dimension(4)::re
 integer*4,dimension(9)::in
1001 format((f5.2,3(2x,f21.12),2i3,2x,3I3,2X,3I3,i5))
	open(1,file='gxfile',form='formatted')
	open(2,file='gxfile_opt',form='formatted')
	re(1)=1.0
	n_atoms=0
	do while(re(1).gt.0.0)
	read(1,*) re,in
	write(2,1001) re,in
	if(in(1).ne.0) n_atoms=n_atoms+1
	enddo
	read(1,*) re(1:2),in(1:3)
	write(2,'(2f21.12,i5,2i5)') re(1:2),in(1:3)
	do n=1,n_atoms
	read(1,*) in(1), re(1:3)
	write(2,'(i5,5x,3f17.12)') in(1), re(1:3)
	enddo
	end

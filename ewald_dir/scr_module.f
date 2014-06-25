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
	module scr_module
	real*8, dimension(3,3) :: rr 
	real*8, dimension(3) :: ds != (/0.0, 0.0, 0.0/) 
	logical	:: lhds
	data	rr /1.0, 0.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/
	data ds /0.0, 0.0, 0.0/, lhds/.false./
	contains 
	subroutine  scrot_ev(NN,ABC,is)
	real*8, intent(inout), dimension(3,*) :: 	ABC
	integer, intent(in) :: nn, is
	real*8, dimension(3) :: rw
	      DO N=1,NN

	if(is.lt.0) then
	abc(1:3,n)=abc(1:3,n)-ds(1:3)
	rw(1:3)=0.0 
	else
	rw(1:3)=ds(1:3)
	endif

      DO  J=1,3
	if(is.lt.0) then
	RW(J)=RW(J)+dot_product(ABC(1:3,n),RR(1:3,j))
	else
	RW(J)=RW(J)+dot_product(ABC(1:3,n),RR(J,1:3))
	endif
	enddo	! 	J=1,3			

        ABC(1:3,N)=RW(1:3)
	enddo	! N=1,NN
      RETURN
      E N D	subroutine  scrot_ev
      end module scr_module

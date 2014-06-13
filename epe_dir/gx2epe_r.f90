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
	implicit real*8 (a-h,o-z)
        real*8,	dimension(3,300):: r_gx
        real*8, allocatable, dimension(:,:):: r_ewa
        integer*4 ::ind(9)
       open(7,file='ewald.pcr')
       read(7,*) n_ewpc
       allocate(r_ewa(3,n_ewpc))
	do i=1,n_ewpc
	 read(7,*) r_ewa(:,i)
        enddo
	
	ico_r=0
       do ico=1,300
	read(5,*) an,r_gx(:,ico),ind(1:9)
	
	if(an.lt.0.5) then
	 exit
	else if(ind(1).ne.0.and.ind(9).eq.-1) then
	do i=1,n_ewpc
	 if(dot_product(r_gx(:,ico)-r_ewa(:,i),r_gx(:,ico)-r_ewa(:,i)).lt.0.9) then
          ico_r=ico_r+1
	  print*,r_ewa(:,i),ico_r
	 endif
	enddo	
	  	
	endif
       enddo
    end

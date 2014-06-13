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
      integer*4 ind(9,200)
      real*8 a(200),xyz(3,200)
      do i=1,200
      read(5,*) a(i), xyz(:,i),ind(:,i)
      if(a(i).lt.0.5) exit
      enddo

      j=0
      do i=1,200
      if(ind(6,i).ne.0) then
      j=j+1
      ind(6,i)=j
      endif
      if(a(i).lt.0.5) exit
      enddo

      do i=1,200
      if(ind(7,i).ne.0) then
      j=j+1
      ind(7,i)=j
      endif
      if(a(i).lt.0.5) exit
      enddo
      
      do i=1,200
      if(ind(8,i).ne.0) then
      j=j+1
      ind(8,i)=j
      endif
      if(a(i).lt.0.5) exit
      enddo

      do i=1,200
      write(6,'(f7.2,3f20.12,2i4, i5,2i4,i5,2i4,i5)') a(i), xyz(:,i),ind(:,i)
      if(a(i).lt.0.5) exit
      enddo
      
      end

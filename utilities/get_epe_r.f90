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
program get_epe_r

  !Small utilite to prepare epe.r from gx file
  ! get_epe_r < gxfile

  character*256 :: buffer
  integer*4 :: i,n1(8),n2
  real*8 :: an,x,y,z

  open(10,file="epe.r")

  i=0
  do
     i=i+1
     read(5,'(a256)',err=100) buffer
     read(buffer,'(f5.2)',err=100) an
     if(an <= 0.0d0) exit
     read(buffer,*,err=100) an,x,y,z,n1,n2
     if(n2 == 0) cycle
     if(n2 == -1) write(10,'(3f25.12)') x,y,z
  end do

  close(10)
  stop

100 print*,'Error of reading in GXFILE. Line ',i

end program get_epe_r

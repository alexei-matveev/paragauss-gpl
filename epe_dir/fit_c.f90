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
program main
real*8 x(10),f(10),x_ref(10),f_ref(10),f_o(10)
call system('rm dfp*')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u+0.4')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u+0.3')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u+0.2')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u+0.1')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u+0.0')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u-0.1')
call system('/home/bin/sim_V2.0vn_mpi alt-o-alt-u-0.2')
call system('grep " EPE correcred total energy" dfp* > grep.res')
open(5,file='./grep.res',form='formatted')
open(7,file='al-o-alt-u.res.dat',form='formatted')
do i=1,7
read(5,'(15x,f4.2,1x,f20.5)'),x(i),f(i)
print*,x(i),f(i)
enddo
fm=f(2)
do i=1,7
f(i)=f(i)-fm
print*,x(i),f(i)
enddo
do i=1,7
read(7,*) x_ref(i),f_ref(i)
print*,x_ref(i),f_ref(i)
enddo
f_o(1)=f(7)
f_o(2)=f(6)
f_o(3:7)=f(1:5)
sqerr=sqerr+dot_product(f_o(1:7)-f_ref(1:7),f_o(1:7)-f_ref(1:7))
print*,sqerr
do i=1,7
print*,x_ref(i),f_o(i)
enddo

end program main


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
! driver for numerical estimation of FF hessian
character(len=3) symv,ssymv
logical:: converged
call system('calc_ffh; rm conv;cat gxfile')
print *,'input number of internals'
read(5,*) n_intern
do i=1,n_intern*2+1
call system('optimizer_FFH FFH FFH FFH')
write(symv,'(i3)') i
ssymv=adjustl(symv)
print *,'cp flepo.stat flepo'//ssymv
call system('cp ./flepo.stat flepo'//ssymv)
inquire(file='conv',exist=converged)
if(converged) exit
enddo
call system('make_optimization;echo "please privide gx and  recover files"')
end

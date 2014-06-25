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
program epe2xyz
!Convert epe environment (region IIa) in XYZ format
!To be used files epe.pcs and epe.pcc has to be in the current directory
  character*2 :: name(98)= &
       (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne", &
       "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca", &
       "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
       "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr", &
       "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn", &
       "Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd", &
       "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tu","Yb", &
       "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
       "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", &
       "Pa","U ","Np","Pu","Am","Cm","Bk","Cf"/)
  
  open(10,file='epe.pcs')
  open(11,file='epe.pcc')
      
  read(10,'(i5)') NN
  read(11,'(i5)') NN
  write(6,'(i5)') NN
  write(6,*) '!!!!!!!!!!!!'
  
  do i=1,nn
     read(10,'(3f17.8,29x,i2)') xx,yy,zz,na
     read(11,'(3f17.8)') xx,yy,zz
     write(6,'(a2,3f17.8)') name(na),xx*0.529177249d0,yy*0.529177249d0,zz*0.529177249d0
  enddo

  close(10)
  close(11)

end program epe2xyz

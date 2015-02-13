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
program gx2xyz
!Usage: gx2xyz_a [-99] [-efp] < gxfile > aaa.xyz
  real*8 :: an,xx,yy,zz
  dimension :: an(10000),xx(10000),yy(10000),zz(10000)
  logical :: dummy,efp
  integer*4 :: i,j,k,na,n_args
  character(len=4) :: args(2)
  character(len=2) :: name(99)= &
       (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne", &
         "Na","Mg","Al","Si","P ","S ","Cl","Ar","K ","Ca", &
         "Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
         "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr", &
         "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn", &
         "Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd", &
         "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tu","Yb", &
         "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg", &
         "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", &
         "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","XX"/)

  args='   '
  dummy=.false.; efp=.true.
  n_args=COMMAND_ARGUMENT_COUNT()
  if(n_args == 1) call get_command_argument(n_args,args(1))
  if(n_args == 2) call get_command_argument(n_args,args(2))
  if(trim(args(1)) == '-99' .or. trim(args(2)) == '-99') dummy=.true.
  if(trim(args(1)) == '-no_efp' .or. trim(args(2)) == '-no_efp') efp=.false.

  i=1
  do 
     read(5,*,err=2) an(i),xx(i),yy(i),zz(i)
     if(.not.efp .and. int(an(i)*100.0d0)==801) cycle
     if(.not.efp .and. int(an(i)*100.0d0)==102) cycle
     if(.not.efp .and. int(an(i)*100.0d0)==103) cycle
     na=int(an(i))
     if(na == 99 .and. .not. dummy) cycle
     if(na <= 0) exit
     i=i+1
  enddo

2 k=i-1
  write(6,'(i4)') k
  write(6,*) '!!!!!!!!!!!!!!!!'
  do j=1,k
     na=int(an(j))
     write(6,'(a2,3f21.12)') name(na),xx(j)*0.529177249d0,yy(j)*0.529177249d0,zz(j)*0.529177249d0
  enddo

end program gx2xyz

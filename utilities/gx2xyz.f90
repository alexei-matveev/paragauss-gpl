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
!Usage: gx2xyz [-99] < gxfile > aaa.xyz
        dimension :: an(10000),xx(10000),yy(10000),zz(10000)
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
        character(len=255) :: cmd
        character(len=3) :: key99

        key99='  '
        n_args=COMMAND_ARGUMENT_COUNT()
        if(n_args > 0) call get_command_argument(n_args,key99)
 
        i=1
        do 
                read(5,*) an(i),xx(i),yy(i),zz(i)
                na=int(an(i))
                if(na == 99.and.key99 /= '-99') cycle
                if(na <= 0) exit
                i=i+1
        enddo

        write(6,'(i4)') i-1
        write(6,*) '!!!!!!!!!!!!!!!!'
        do j=1,i-1
                na=int(an(j))
                write(6,'(a2,3f17.8)') name(na),xx(j)*0.529177249d0,yy(j)*0.529177249d0,zz(j)*0.529177249d0
        enddo
end

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
program xyz2gx
  
  implicit none

  integer*4 :: Ngx,nat_gx(3000)
  integer*4 :: Nxyz
  real*8 :: coor_xyz(3000,3)
  character*300 :: xyz_name,gx_name
  character*2 :: atm_name(3000)
  character*1 :: yes_no
  logical :: high
  
  character*80 :: comment
  integer*4 :: i,j,k

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

  write(6,*) 'Number of atoms in GX file'
  read(5,*) ngx
  write(6,*) 'Specify numbers of atoms in XYZ file'
  read(5,*,err=1) (nat_gx(i),i=1,ngx)

  write(6,*) 'Name of XYZ file'
  read(5,*) xyz_name

  open(10,file=trim(xyz_name))
  read(10,*) Nxyz
  read(10,*) comment

  do i=1,Nxyz
     read(10,*) atm_name(i),coor_xyz(i,:)
  end do
  close (10)

  goto 2
1 ngx=Nxyz

2 write(6,*) 'Name of GX file'
  read(5,*) gx_name
3 write(6,*) 'High precision?'
  read(5,*,err=3) yes_no
  if (trim(yes_no)=='Y' .or. trim(yes_no)=='y') then
     high=.true.
  else if (trim(yes_no)=='N' .or. trim(yes_no)=='n') then
     high=.false.
  else
     goto 3
  endif

  open(10,file=trim(gx_name))

  do i=1,ngx
     if(ngx == Nxyz) then
        j=i
     else
        j=nat_gx(i)
     end if
     do k=1,98
        if(atm_name(j) == name(k)) exit
     end do
     if (high) then
        write(10,'(f5.2,3f23.12,2i4,2x,3i3,2x,3i3,i5)') dble(k),coor_xyz(j,:)/0.529177d0,i,i,0,0,0,0,0,0,0
     else
        write(10,'(f5.2,3f15.7,2i4,1x,3i4,1x,3i4,i5)') dble(k),coor_xyz(j,:)/0.529177d0,i,i,0,0,0,0,0,0,-1
     end if
  end do
  if (high) then
     write(10,'(f5.2,3f23.12,2i4,2x,3i3,2x,3i3,i5)') -1.0,0.0,0.0,0.0,0,0,0,0,0,0,0,0,0
  else
     write(10,'(f5.2,3f15.7,2i4,2x,3i3,2x,3i3,i5)') -0.0,0.0,0.0,0.0,0,0,0,0,0,0,0,0,1
  end if

  close(10)

end program xyz2gx




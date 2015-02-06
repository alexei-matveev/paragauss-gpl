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
program xyz2gx
  
  implicit none

  integer*4 :: Ngx,nat_gx(3000),nskip,n_sk(3000)
  integer*4 :: Nxyz
  real*8 :: coor_xyz(3000,3)
  character*300 :: xyz_name,gx_name
  character*2 :: atm_name(3000)
  character*1 :: yes_no
  logical :: high,stat
  
  character*80 :: comment
  integer*4 :: i,j,k,l,what_to_do

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
  character*2 :: name1(98)= &
       (/"H ","HE","LI","BE","B ","C ","N ","O ","F ","NE", &
         "NA","MG","AL","SI","P ","S ","CL","AR","K ","CA", &
         "SC","TI","V ","CR","MN","FE","CO","NI","CU","ZN", &
         "GA","GE","AS","SE","BR","KR","RB","SR","Y ","ZR", &
         "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN", &
         "SB","TE","I ","XE","CS","BA","LA","CE","PR","ND", &
         "PM","SM","EU","GD","TB","DY","HO","ER","TU","YB", &
         "LU","HF","TA","W ","RE","OS","IR","PT","AU","HG", &
         "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH", &
         "PA","U ","NP","PU","AM","CM","BK","CF"/)

  stat=.true.

  write(6,*) 'Name of XYZ file'
  read(5,*) xyz_name

  open(10,file=trim(xyz_name))
  read(10,*) Nxyz
  read(10,*) !comment

  do j=1,Nxyz
     read(10,*) atm_name(j),coor_xyz(j,:)
  end do
  close (10)

1 write(6,*) 'All atoms to GX(1), some atoms to GX(2), some atoms skip(3)'
  read(5,*) what_to_do

  if(what_to_do==1) then
     stat=.false.
     ngx=Nxyz
  elseif(what_to_do==2) then
     write(6,*) 'Number of atoms in GX file'
     read(5,*) ngx
     write(6,*) 'Specify atomic numbers in XYZ file'
     read(5,*) (nat_gx(i),i=1,ngx)
  elseif(what_to_do==3) then
     write(6,*) 'Number of skiped atoms'
     read(5,*) nskip
     write(6,*) 'Specify skiped atomic numbers in XYZ file'
     read(5,*) (n_sk(i),i=1,nskip)
     ngx=Nxyz-nskip
     k=1;l=1
     do 
        mm:do i=l,Nxyz
           do j=1,nskip
              if(i==n_sk(j)) cycle mm
           enddo
           nat_gx(k)=i
           l=i+1
           exit mm
        enddo mm
        if(k==ngx) exit
        k=k+1
     enddo
  else
     goto 1
  endif

  write(6,*) 'Name of GX file'
  read(5,*) gx_name

  open(10,file=trim(gx_name))

  do i=1,ngx
     if(.not.stat) then
        j=i
     else
        j=nat_gx(i)
     end if
     do k=1,98
        if(atm_name(j) == name(k).or.atm_name(j) == name1(k)) exit
     end do
     write(10,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4,i5)') dble(k),coor_xyz(j,:)/0.529177249d0,i,i,0,0,0,0,0,0,0
  end do
  write(10,'(f5.2,3f23.12,2i4,2x,3i4,2x,3i4,i5)') -1.0,0.0,0.0,0.0,0,0,0,0,0,0,0,0,0

  close(10)

end program xyz2gx




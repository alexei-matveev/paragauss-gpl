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
program check_epe
  !The program checks coinsidence of atomic centers (regular
  !positions) collected in epe.r, epe.pcr and ewald.pcr
  !files

  real*8 :: r_epe_r(300,3),r_epe_pcr(1000,3),r_ewald(10000,3)
  real*8 :: r1(3),r2(3),dr,dmin
  integer*4 :: n_r,n_pcr,n_ew,i_r,i_pcr
  integer*4 :: i,j,imin
  real*8, parameter :: small=1.0d-2,large=1.0d4
  logical :: check

  r_epe_r=0.0d0; r_epe_pcr=0.0d0; r_ewald=0.0d0
  
  open(10,file='epe.r')
  open(11,file='epe.pcr')
  open(12,file='ewald.pcr')

  n_r=1
  do
     read(10,*,end=100) r_epe_r(n_r,:) 
     n_r=n_r+1
  enddo
100 n_r=n_r-1
  
  read(11,*) n_pcr
  do i=1,n_pcr
     read(11,*) r_epe_pcr(i,:)
  enddo
  
  read(12,*) n_ew
  do i=1,n_ew
     read(12,*)r_ewald(i,:)
  enddo
  
  close(10)
  close(11)
  close(12)

  print*,'Number centers in epe.r :',n_r
  print*,'Number centers in epe.pcr :',n_pcr
  print*,'Number centers in ewald.pcr :',n_ew

  i_r=0
  m1: do i=1,n_r
     check=.false.
     dmin=large
     r1=r_epe_r(i,:)
     m2: do j=1,n_ew
        r2=r_ewald(j,:)
        dr=dot_product(r1-r2,r1-r2)
        if(dr < dmin) then
           dmin=dr
           imin=j
        endif
        if(dr < small) then
           i_r=i_r+1
           check=.true.
           cycle m1
        endif
     enddo m2
     if(.not. check) then
        print*,'epe.r center ',i,' does not coinside with any ewald.pcr center'
        print*,'closest ewald.pcr center is ',imin,', on the distance: ', sqrt(dmin)
     endif
  enddo m1
  print*,i_r,' epe.r centers coinside with ewald.pcr centers'

  i_pcr=0
  l1: do i=1,n_pcr
     check=.false.
     dmin=large
     r1=r_epe_pcr(i,:)
     l2: do j=1,n_ew
        r2=r_ewald(j,:)
        dr=dot_product(r1-r2,r1-r2)
        if(dr < dmin) then
           dmin=dr
           imin=j
        endif
        if(dr < small) then
           i_pcr=i_pcr+1
           check=.true.
           cycle l1
        endif
     enddo l2
     if(.not. check) then
        print*,'epe.pcr center ',i,' does not coinside with any ewald.pcr center'
        print*,'closest ewald.pcr center is ',imin,', on the distance: ', sqrt(dmin)
     endif
  enddo l1
  print*,i_pcr,' epe.pcr centers coinside with ewald.pcr centers'
  print*,'Total number of coinsided centers: ',i_pcr+i_r
  
end program check_epe

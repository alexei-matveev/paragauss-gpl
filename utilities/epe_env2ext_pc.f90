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
program epe_env2ext_pc
  !this routine combines epe environment files (epe.pcc, epe.pcs, ewald.pcr)
  !and converts them to MolMech format as external point charges

  integer*4 :: i,j
  integer*4 :: nn,nn1,nn2,ii(5),ntypes
  real*8, allocatable :: r_ref(:,:)
  real*8 :: r(3),q,rr(3)
  character*5 :: name
  character*5 :: buf
  real*8 :: start_name,types(999),curr_type
  real*8, parameter :: b2a=0.529177249d0
  real*8, parameter :: small=1.0d-4

  open(10,file="ext_pc")
  write(10,'(a8)') "&pc_coor"

  open(11,file="epe.pcs")
  read(11,*) nn
  do i=1,nn
     read(11,'(4f17.8,1x,5i2,1x,a5)') r,q,ii(1:5),name
     buf=adjustl(name)
     write(10,'(a5,4f17.8)') buf,r*b2a,q
  end do
  close(11)
  print*,"epe.pcs",nn

  start_name=0d0
  ntypes=0

  open(12,file="epe.pcc")
  read(12,*) nn
  print*,"epe.pcc",nn
  nn2=0
  do i=1,nn
     read(12,'(4f17.8,1x,5i2)') r,q,ii(1:5)
     if(q == 0d0) cycle
     curr_type=get_type(q)
     write(buf,'(f5.3)') curr_type
     write(10,'(a5,4f17.8)')adjustl(buf),r*b2a,q
     nn2=nn2+1
  end do
  close(12)
  print*,"epe.pcc",nn2,"saved,",nn-nn2,"zero charged droped"

  open(13,file="epe.pcr")
  read(13,*) nn
  allocate(r_ref(3,nn+200))
  do i=1,nn
     read(13,'(3f17.8)') r_ref(1:3,i)
  end do
  close(13)
  print*,"epe.pcr",nn

  open(14,file="epe.r")
  i=nn
  do
     i=i+1
     read(14,*,end=10) r_ref(1:3,i)
  end do
10 print*,"epe.r",i-1-nn
  nn1=i-1
  print*,"epe.r+epe.pcr",nn1
  close(14)

  open(15,file="ewald.pcr")
  read(15,*) nn
  print*,"ewald.pcr ",nn
  nn2=0
  l1: do i=1,nn
     read(15,'(4f15.8,1x,5i2)') r,q,ii(1:5)
     do j=1,nn1
        rr=r-r_ref(:,j)
        if(dot_product(rr,rr) <= small) then
           nn2=nn2+1
           cycle l1
        end if
     end do
     if(q == 0d0) cycle
     curr_type=get_type(q)
     write(buf,'(f5.3)') curr_type
     write(10,'(a5,4f17.8)')adjustl(buf),r*b2a,q
  end do l1
  deallocate(r_ref)
  close(15)
  print*,"ewald.pcr",nn-nn2," saved,",nn2,"droped"

  write(10,'(a8)') "/pc_coor"
  close(10)

contains

  function get_type(qi)
    real*8 :: get_type
    real*8 :: qi
    integer*4 :: ig

    if(ntypes==0) then
       ntypes=1
       get_type=ntypes*1.0d-3
       types(1)=qi
    else
       do ig=1,ntypes
          if (qi == types(ig)) then
             get_type=ig*1.0d-3
             return
          end if
       end do
       get_type=(ntypes+1)*1.0d-3
       types(ntypes+1)=qi
       ntypes=ntypes+1
    end if

  end function get_type

end program epe_env2ext_pc

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
program build_mcm41
!Generate MCM-41 unit cell of different size of hexagonal channel
  implicit none

  real*8, parameter :: pi=3.14159265355897932368
  real*8, parameter :: tet_angle_r=1.9106332362490185960979261
  real*8, parameter :: tet_angle_d=109.4712206355644070754351560
  real*8, parameter :: angle_120=120.0
  real*8, parameter :: angle_60=60.0
  real*8, parameter :: angle_30=30.0
  real*8, parameter :: Si_O_bond=1.626
  real*8, parameter :: zero=0.0
  real*8, parameter :: small=1.0e-7

  real*8 :: ri(5,3),angle_120_r,angle_60_r,angle_30_r
  real*8 :: trans,length,X,Y,t_x,t_y,plane
  integer*4 :: wall
  character*2 :: atom_nmi(5)
  real*8 :: trans_vec(3,3)
  integer*4 :: save_format

  real*8, allocatable :: xyz(:,:),coor(:,:)
  character*2, allocatable :: atom_name(:),at_nm_c(:)

  integer*4 :: ind,ind1,i,j,k,l

  angle_120_r=pi*angle_120/180.0_8
  angle_60_r=pi*angle_60/180.0_8
  angle_30_r=pi*angle_30/180.0_8

  !generating the coordinates of the initial tetragedron
  ri(1,1)=zero
  ri(1,2)=zero
  ri(1,3)=zero

  ri(2,1)=Si_O_bond
  ri(2,2)=zero
  ri(2,3)=zero

  ri(3,1)=ri(2,1)*cos(tet_angle_r)+ri(2,3)*sin(tet_angle_r)
  ri(3,2)=zero
  ri(3,3)=-ri(2,1)*sin(tet_angle_r)+ri(2,3)*cos(tet_angle_r)

  ri(4,1)=ri(3,1)
  ri(4,2)=ri(3,2)*cos(angle_120_r)+ri(3,3)*sin(angle_120_r)
  ri(4,3)=-ri(3,2)*sin(angle_120_r)+ri(3,3)*cos(angle_120_r)

  ri(5,1)=ri(3,1)
  ri(5,2)=ri(3,2)*cos(2*angle_120_r)+ri(3,3)*sin(2*angle_120_r)
  ri(5,3)=-ri(3,2)*sin(2*angle_120_r)+ri(3,3)*cos(2*angle_120_r)

  ri(:,3)=ri(:,3)-ri(3,3)

  atom_nmi(1)='Si'
  atom_nmi(2:5)='O '

  trans=abs(ri(4,2)-ri(5,2))
  plane=ri(4,3)

  do
     write(*,*) 'The number of tetragedrons into the MCM-41 walls (>= 3)'
     read(5,*) wall
     if(wall >= 3) exit
  enddo

  length=wall*trans
  X=length*cos(angle_30_r)
  Y=length*sin(angle_30_r)
  t_x=X-ri(4,1)
  t_y=Y+ri(4,2)

  ri(:,1)=ri(:,1)+t_x
  ri(:,2)=ri(:,2)-t_y

  !definition of the vectors of a translation
  ! a - vector
  trans_vec(1,1)=ri(2,1)*2
  trans_vec(1,2)=zero
  trans_vec(1,3)=zero
  ! b - vector
  trans_vec(2,1)=trans_vec(1,1)*cos(-angle_120_r)+trans_vec(1,2)*sin(-angle_120_r)
  trans_vec(2,2)=-trans_vec(1,1)*sin(-angle_120_r)+trans_vec(1,2)*cos(-angle_120_r)
  trans_vec(2,3)=zero
  ! c - vector
  trans_vec(3,1)=zero
  trans_vec(3,2)=zero
  trans_vec(3,3)=(ri(3,3)+2*(plane-ri(3,3)))*2

  !generating MCM-41 double ring
  allocate(xyz(6*wall*5*2,3),coor(60*wall,3))
  allocate(atom_name(6*wall*5*2),at_nm_c(60*wall))

  lab1: do l=1,6
     lab2: do i=1,wall
        lab3: do k=1,2
           lab4: do j=1,5
              ind=10*wall*(l-1)+10*(i-1)+5*(k-1)+j
              atom_name(ind)=atom_nmi(j)
              if(mod(i,2) /= 0 .and. k==2 .and. j==2) then
                 xyz(ind,:)=zero
                 cycle lab4
              endif
              if(mod(i,2) == 0 .and. k==1 .and. j==2) then
                 xyz(ind,:)=zero
                 cycle lab4
              endif
              if(mod(i,2) == 0 .and. mod(i,4) /= 0 .and.k==1 .and. j==3) then
                 xyz(ind,:)=zero
                 cycle lab4
              endif
              if(mod(i,4) == 0 .and.k==2 .and. j==3) then
                 xyz(ind,:)=zero
                 cycle lab4
              endif
              xyz(ind,1)=ri(j,1)*cos((l-1)*angle_60_r)+ &
                   (ri(j,2)+trans*(i-1))*sin((l-1)*angle_60_r)
              xyz(ind,2)=-ri(j,1)*sin((l-1)*angle_60_r)+ & 
                   (ri(j,2)+trans*(i-1))*cos((l-1)*angle_60_r)
              if(mod(wall,2) == 0) then
                 if(mod(i,2) /= 0 .and. (j /= 4 .and. j /= 5)) then
                    xyz(ind,3)=ri(j,3)
                 else
                    xyz(ind,3)=ri(j,3)+2*(plane-ri(j,3))
                 endif
              else
                 if(j /= 2) then
                    if(mod(i,2) /= 0 .and. (j /= 4 .and. j /= 5)) then
                       xyz(ind,3)=ri(j,3)
                    else
                       xyz(ind,3)=ri(j,3)+2*(plane-ri(j,3))
                    endif
                 else
                    if(mod(l,2) /= 0) then
                       if(mod(i,2) /= 0) then
                          xyz(ind,3)=ri(j,3)
                       else
                          xyz(ind,3)=ri(j,3)+2*(plane-ri(j,3))
                       endif
                    else
                       if(mod(i,2) /= 0) then
                          xyz(ind,3)=ri(j,3)-4*(plane-ri(j,3))
                       else
                          xyz(ind,3)=ri(j,3)-6*(plane-ri(j,3))
                       endif
                    endif
                 endif
              endif
              if(k==2) xyz(ind,3)=-xyz(ind,3)
           enddo lab4
        enddo lab3
     enddo lab2
  enddo lab1

  !excluding coincided atoms 
  ind1=1
  coor(ind1,:)=xyz(1,:)
  at_nm_c(ind1)=atom_name(1)
  met1: do j=2,ind
     if(xyz(j,1) == zero .and. &
        xyz(j,2) == zero .and. &
        xyz(j,3) == zero) cycle met1
     do k=1,ind1
        if(abs(xyz(j,1)-coor(k,1)) <= small .and. &
           abs(xyz(j,2)-coor(k,2)) <= small .and. &
           abs(xyz(j,3)-coor(k,3)) <= small) cycle met1
     enddo
     ind1=ind1+1
     coor(ind1,:)=xyz(j,:)
     at_nm_c(ind1)=atom_name(j)
  enddo met1

  do i=1,ind1
     coor(i,:)=coor(i,:)+trans_vec(1,:)+trans_vec(3,:)+trans_vec(2,:)
  enddo

  !saving result
  do 
     write(*,*) 'save rezult in'
     write(*,*) 'XYZ[1], DL_POLY COOR[2], DL_POLY CONFIG[3]:'
     read(5,*) save_format
     if(save_format==1 .or. &
        save_format==2 .or. &
        save_format==3) exit
  enddo

  select case(save_format)
     case(1)
        !saving result in XYZ format
        open(10,file='mcm_41.xyz')
        write(10,'(i4)') ind1
        write(10,'(a6)') 'MCM-41'
        do i=1,ind1
           write(10,'(a2,3(3x,f13.9))') at_nm_c(i),coor(i,:)
        enddo
        close(10)
     case(2)
        !saving result in DL_POLY COOR format
        open(10,file='COOR')
        write(10,'(a4)') 'CART'
        write(10,'(a6)') 'MCM-41'
        write(10,'(i4)') ind1
        do i=1,3
           write(10,'(3f20.8)') trans_vec(i,:)
        enddo
        write(10,'(3i8)') 1,1,1
        do i=1,ind1
           write(10,'(a2,6x,3f20.8)') at_nm_c(i),coor(i,:)
        enddo
        close(10)
     case(3)
        !saving result in DL_POLY CONFIG format
        open(10,file='COOR_CONF')
        write(10,'(a6)') 'MCM-41'
        write(10,'(2i10)')0,3
        do i=1,3
           write(10,'(3f20.8)') trans_vec(i,:)
        enddo
        write(10,'(3i8)') 1,1,1
        do i=1,ind1
           write(10,'(a2,6x,i10/,3f20.8)') at_nm_c(i),1,coor(i,:)
        enddo
        close(10)
     endselect

end program build_mcm41


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
program punch2cube

  integer:: i_p,i_t,N_at,N_con,N,N_grid(3),N_lines
  character(len=128) :: buffer1
  character(len=128) :: buffer2
  character(len=5) :: Atom_name(1000),Atoms(5)
  integer :: Atom_num(5),At_num,i1,i2,i3,count,count1
  integer, parameter :: N_tot_at=5
  real*8 :: core_charge(5),Core_chg
  character(len=20) :: name
  real*8 :: xyz(3),x0(3),x1(3),x2(3),x3(3),d1(3),d2(3),d3(3),p_d,c_d
  real*8, allocatable :: punch_data(:), cube_data(:,:,:)
  real*8, parameter :: b2a=0.529177249d0

  read (5,*) buffer1
  read (5,*) buffer1
  write(6,'(a37)') 'ParaGauss Cube File - standard format'
  write(6,'(a37)') 'ParaGauss Cube File - standard format'

  read (5,'(i5,3F12.6)'),N_at,x0
  write(6,'(i5,3F12.6)'),N_at,x0

  read (5,'(i5,3F12.6)') N_grid(1),d1
  read (5,'(i5,3F12.6)') N_grid(2),d2
  read (5,'(i5,3F12.6)') N_grid(3),d3
  write(6,'(i5,3F12.6)') N_grid(1),d1
  write(6,'(i5,3F12.6)') N_grid(2),d2
  write(6,'(i5,3F12.6)') N_grid(3),d3

  do i_p=1,N_at
     read (5,'(i5,4F12.6)') At_num,Core_chg,xyz(1:3)
     write(6,'(i5,4F12.6)') At_num,Core_chg,xyz(1:3)
  end do

  allocate(punch_data(N_grid(1)*N_grid(2)*N_grid(3)))
  N_lines = N_grid(1)*N_grid(2)*N_grid(3)
  do i_p=1,N_lines
     read(5,'(E13.5)',err=20) punch_data(i_p)
     cycle
20   punch_data(i_p)=9999.d0
  end do

  allocate(cube_data(N_grid(1),N_grid(2),N_grid(3)))
  count=0
  do i1=1,N_grid(1)
     do i2=1,N_grid(2)
        do i3=1,N_grid(3)
           count=count+1
           cube_data(i1,i2,i3)=punch_data(count)
        end do
     end do
  end do

  deallocate(punch_data)
  do i1=1,N_grid(1)
     do i2=1,N_grid(2)
        write(6,'(6E13.5)') (cube_data(i1,i2,i3), i3=1,N_grid(3))
     end do
  end do

  deallocate(cube_data)

end program punch2cube

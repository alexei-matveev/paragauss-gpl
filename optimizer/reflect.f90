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
implicit none
real*8 point_to_reflect(3), origin(3), reflect_plane_vector(3)
real*8 housholder_matrix(3,3), reflected_point(3),r2
write(*,*) 'point_to_reflect'
read(*,*) point_to_reflect
write(*,*) 'origin'
read(*,*) origin
write(*,*) 'reflect_plane_vector'
read(*,*) reflect_plane_vector

reflect_plane_vector=reflect_plane_vector-origin
point_to_reflect=point_to_reflect-origin
r2=sum( reflect_plane_vector**2)
housholder_matrix(1,1)=1-2*reflect_plane_vector(1)**2/r2
housholder_matrix(2,2)=1-2*reflect_plane_vector(2)**2/r2
housholder_matrix(3,3)=1-2*reflect_plane_vector(3)**2/r2
housholder_matrix(1,2)=-2*reflect_plane_vector(1)*reflect_plane_vector(2)/r2
housholder_matrix(2,1)=-2*reflect_plane_vector(1)*reflect_plane_vector(2)/r2
housholder_matrix(1,3)=-2*reflect_plane_vector(1)*reflect_plane_vector(3)/r2
housholder_matrix(3,1)=-2*reflect_plane_vector(1)*reflect_plane_vector(3)/r2
housholder_matrix(2,3)=-2*reflect_plane_vector(2)*reflect_plane_vector(3)/r2
housholder_matrix(3,2)=-2*reflect_plane_vector(2)*reflect_plane_vector(3)/r2
reflected_point=matmul(point_to_reflect,housholder_matrix)
reflected_point=reflected_point+origin
write(*,*) 'reflected_point'
write(*,'(5x,3f23.12)') reflected_point
end

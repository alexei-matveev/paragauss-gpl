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
real*8 point_to_move(3), origin(3), direction_vector(3)
real*8 diplacement_length, moved_point(3),r2
write(*,*) 'point_to_move'
read(*,*) point_to_move
write(*,*) 'origin'
read(*,*) origin
write(*,*) 'direction_vector'
read(*,*) direction_vector
write(*,*) 'diplacement_length'
read(*,*) diplacement_length

direction_vector=direction_vector-origin
direction_vector=direction_vector/sqrt( sum(direction_vector**2) )
moved_point=point_to_move+diplacement_length*direction_vector
write(*,*) 'moved_point'
write(*,'(5x,3f23.12)') moved_point
end

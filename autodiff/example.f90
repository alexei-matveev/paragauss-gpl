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
program example
    use ad2x3
    implicit none

    ! AD-variables:
    type(ad) :: x, y, z
    type(ad) :: a, b, c, f

    ! independent variables, with IDs and values:
    x = var(1, 2.0D0)
    y = var(2, 3.0D0)
    z = var(3, 5.0D0)

    a = sqrt(x*x + y*y + z*z)

    b = 1.0D0 / ( 1.0D0 + exp(a) )

    c = sin(cos(b))

    f = a + b * c

    f = f**2

    print *, "f=", val(f)

    print *, "df/dx=", fst(1, f)
    print *, "df/dy=", fst(2, f)
    print *, "df/dz=", fst(3, f)

    print *, "d2f/dxx=", snd(1, 1, f)
    print *, "d2f/dyx=", snd(2, 1, f)
    print *, "d2f/dzx=", snd(3, 1, f)
    print *, "d2f/dxy=", snd(1, 2, f)
    print *, "d2f/dyy=", snd(2, 2, f)
    print *, "d2f/dzy=", snd(3, 2, f)
    print *, "d2f/dxz=", snd(1, 3, f)
    print *, "d2f/dyz=", snd(2, 3, f)
    print *, "d2f/dzz=", snd(3, 3, f)
end program example

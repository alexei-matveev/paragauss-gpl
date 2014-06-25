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
module time_module

	use options_module


	type timeObj
		real(kind=r8_kind) :: t
	end type timeObj


contains

  type(timeObj) function timeconvert(out_date, out_time)
    !# Function arguments
    character(len=8), intent(out), optional  :: out_date
    character(len=10), intent(out), optional :: out_time

    !# Function variables
    character(len=8)  :: str_date
    character(len=10) :: str_time

    !# Subroutine variables
    integer :: dummy


    call date_and_time( str_date, str_time )

    !# Read miliseconds
    read( str_time(8:10), fmt='( I3 )' ) dummy
    timeconvert%t = .001_r8_kind * dble(dummy)

    !# Read seconds
    read( str_time(5:6), fmt='( I2 )' ) dummy
    timeconvert%t = timeconvert%t + dble(dummy)

    !# Read minutes
    read( str_time(3:4), fmt='( I2 )' ) dummy
    timeconvert%t = timeconvert%t + dble(dummy)*60._r8_kind

    !# Read hours
    read( str_time(1:2), fmt='( I2 )' ) dummy
    timeconvert%t = timeconvert%t + dble(dummy)*3600._r8_kind

    if( present(out_date) ) out_date = str_date
    if( present(out_time) ) out_time = str_time

  end function timeconvert


end module time_module

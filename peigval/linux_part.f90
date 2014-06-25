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


  module linux_part

  use type_module, only: IK=>i4_kind, RK=>r8_kind, RK_=>r4_kind

  implicit none

  private

  public :: time_fit, nrow, ncol, k_opt, permitted, k0, f

  contains

!============================================================================
!============================================================================

     real(kind=RK_) function time_fit( procs, dim )

     implicit none

!    Number of processors an dimension of the matrix:
     integer(kind=IK), intent(in) :: procs, dim

     if( procs == 1 ) then
        if ( dim < 200 ) then
           time_fit = 0.0023 * dim
        else if( dim >= 200 .and. dim < 300 ) then
           time_fit = -2.04 + 0.0125 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -8.715 + 0.03475 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -34.565 + 0.08645 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -102.78 + 0.1839 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -340.28 + 0.4214 * dim
        else if( dim >= 1500  ) then
           time_fit = -953.21 + 0.83002 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 2 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = -1.38 + 0.0102 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -3.315 + 0.01665 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -13.265 + 0.03655 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -34.79 + 0.0673 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -152.31 + 0.18482 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -455.66 + 0.38672 * dim
        else if( dim >= 2000  ) then
           time_fit = -1306.24 + 0.81226 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 3 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = -1.9 + 0.0171 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -4.945 + 0.02725 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -13.22 + 0.0438 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -18.47 + 0.0513 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -104.03 + 0.13686 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -120.36 + 0.14774 * dim
        else if( dim >= 2000 ) then
           time_fit = -596.82 + 0.38597 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 4 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = -1.8 + 0.0155 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -5.865 + 0.02905 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -11.565 + 0.04045 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -36.42667 + 0.07597 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -84.08 + 0.12362 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -239.12 + 0.22698 * dim
        else if( dim >= 2000 ) then
           time_fit = -875.94 + 0.54539 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 5 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = -1.3 + 0.0155 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -8.11 + 0.0382 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -4.135 + 0.03025 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -32.38 + 0.0706 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -63.92 + 0.10214 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -185.99 + 0.18352 * dim
        else if( dim >= 2000 ) then
           time_fit = -645.61 + 0.41333 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 6 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = -0.58 + 0.0129 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -8.185 + 0.03825 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = 2.165 + 0.01755 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -34.29333 + 0.06963 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -96.48 + 0.13182 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -120.36 + 0.14774 * dim
        else if( dim >= 2000 ) then
           time_fit = -596.82 + 0.38597 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 8 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = -0.25 + 0.0124 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -1.645 + 0.01705 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -5.92 + 0.0256 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -16.863 + 0.04123 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -106.47 + 0.13084 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -54.63 + 0.09628 * dim
        else if( dim >= 2000 ) then
           time_fit = -699.75 + 0.41884 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 9 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = 0.14 + 0.0115 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -0.88 + 0.0149 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -6.58 + 0.0263 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -28.373 + 0.05743 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -74.6 + 0.10366 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -253.52 + 0.22294 * dim
        else if( dim >= 2000 ) then
           time_fit = -198.36 + 0.19536 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 10 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = 0.45 + 0.01 * dim
        else if( dim >= 300 .and. dim < 500 ) then
           time_fit = -4.41 + 0.0262 * dim
        else if( dim >= 500 .and. dim < 700 ) then
           time_fit = -13.06 + 0.0435 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -16.8633 + 0.048 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -30.75 + 0.06286 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -78.33 + 0.094 * dim
        else if( dim >= 2000 ) then
           time_fit = -280.23 + 0.19549 * dim
        else
           time_fit = 0.
        end if
     else
       time_fit = 0.
     end if

     end function time_fit
!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
     integer(kind=IK) function k_opt( dim )

!    Dimension of the matrix:
     integer(kind=IK), intent(in) :: dim

!    Corresponding optimal number of processors;
!    here it is considered to be 8 for all matrices having dimension
!    above 500 and 1 for the rest.

     if( dim >= 300 .and. dim < 1500  ) then
        k_opt = 6
     else if( dim >= 1500 ) then
        k_opt = 10
     else
        k_opt = 1
     end if

     end function k_opt
!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
     logical function permitted( procs )

     implicit none

     integer(kind=IK), intent(in) :: procs

     if( procs == 1 .or. procs == 2 .or. procs == 3 .or. procs == 4 .or. procs == 5 &
       .or.  procs == 6 .or. procs == 8 .or. procs == 9 .or. procs == 10 ) then
        permitted = .true.
     else
        permitted = .false.
     end if

     end function permitted

!===========================================================================
!===========================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===========================================================================
!===========================================================================

      integer(kind=IK) function k0( dim )

      implicit none

      integer(kind=IK), intent(in) :: dim

      if( dim >= 300 .and. dim < 1500 ) then
         k0 = 4
      else if( dim >= 1500 ) then
         k0 = 3
      end if

      end function k0


!===========================================================================
!===========================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===========================================================================
!===========================================================================
     integer(kind=IK) function f( count, steps )

     implicit none

     integer(kind=IK), intent(in) :: count
     logical, intent(in) :: steps
     if( steps ) then
        f = 1
     else
        f = count
     end if

     end function f
!============================================================================
!============================================================================

!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
      integer(kind=IK) function NROW( procs )  !theo5 case

      implicit none

      integer(kind=IK), intent(in) :: procs
      if ( procs == 2 ) then
           NROW = 2
      else if ( procs == 3 .or. procs == 6 .or. procs == 9 ) then
           NROW = 3
      else if( procs == 4 ) then
           NROW = 4
      else if( procs == 5 .or. procs == 10 ) then
           NROW = 5
      else if( procs == 8 ) then
           NROW = 4
      else
           NROW = 1
      end if

      end function NROW
!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
      integer(kind=IK) function NCOL( procs ) !theo5 case

      implicit none

      integer(kind=IK), intent(in) :: procs
      if ( procs == 2 ) then
           NCOL = 1
      else if ( procs == 3 .or. procs == 4 .or. procs == 5 ) then
           NCOL = 1
      else if( procs == 6 .or. procs == 8 ) then
           NCOL = 2
      else if( procs == 9 ) then
           NCOL = 3
      else if( procs == 10 ) then
           NCOL = 2
      else
           NCOL = procs
      end if

      end function NCOL
!============================================================================
!============================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end module linux_part


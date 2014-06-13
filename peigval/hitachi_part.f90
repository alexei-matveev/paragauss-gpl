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


  module hitachi_part

  use type_module, only: IK=>i4_kind, RK=>r8_kind, RK_=>r4_kind

  implicit none

  private

  public :: time_fit, nrow, ncol, k_opt, permitted, k0, f

  contains

!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
     real(kind=RK_) function time_fit( procs, dim )

     implicit none

!    Number of processors an dimension of the matrix:
     integer(kind=IK), intent(in) :: procs, dim

     if( procs == 1 ) then
        if( dim >= 200 .and. dim < 300 ) then
           time_fit = 0.03 + 0.0019 * dim
        else if( dim >= 100 .and. dim < 200 ) then
           time_fit = -0.11 + 0.0026 * dim
        else if( dim >= 50 .and. dim < 100 ) then
           time_fit = -0.09 + 0.0024 * dim
        else if( dim >= 300 .and. dim < 480 ) then
           time_fit = -1.9333 + 0.0084 * dim
        else if( dim >= 480 .and. dim < 500 ) then
           time_fit = -3.64 + 0.012 * dim
        else if( dim >= 500 .and. dim < 600 ) then
           time_fit = -1.09 + 0.0159 * dim
        else if( dim >= 600 .and. dim < 700 ) then
           time_fit = -9.01 + 0.0216 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -19.0666 + 0.0359 * dim
        else if( dim >= 1000 .and. dim < 1500 ) then
           time_fit = -66.34 + 0.0832 * dim
        else if( dim >= 1500 .and. dim < 2000 ) then
           time_fit = -136.24 + 0.1298 * dim
        else if( dim >= 2000 ) then
           time_fit = -444.96 + 0.2842 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 3 ) then
        if( dim >= 500 .and. dim < 700 ) then
           time_fit = -5.295 + 0.0161 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -14.01 + 0.0286 * dim
        else if( dim >= 1000 .and. dim < 2000 ) then
           time_fit = -65.38 + 0.0799 * dim
        else if( dim >= 2000 ) then
           time_fit = -109.58 + 0.2041 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 4 ) then
        if( dim >= 500 .and. dim < 700 ) then
           time_fit = -4.19 + 0.0129 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -11.61 + 0.0235 * dim
        else if( dim >= 1000 .and. dim < 2000 ) then
           time_fit = -55.08 + 0.0669 * dim
        else if( dim >= 2000 ) then
           time_fit = -239.8 + 0.1593 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 6 ) then
        if( dim >= 500 .and. dim < 700 ) then
           time_fit = -3.43 + 0.0109 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -10.48 + 0.0178 * dim
        else if( dim >= 1000 .and. dim < 2000 ) then
           time_fit = -39.53 + 0.0490 * dim
        else if( dim >= 2000 ) then
           time_fit = -176.15 + 0.1174 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 8 ) then
        if( dim >= 500 .and. dim < 700 ) then
           time_fit = -2.76 + 0.0091 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -6.86 + 0.149 * dim
        else if( dim >= 1000 .and. dim < 2000 ) then
           time_fit = -32.1 + 0.0402 * dim
        else if( dim >= 2000 ) then
           time_fit = -139.8 + 0.0940 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 9 ) then
        if( dim >= 500 .and. dim < 700 ) then
           time_fit = -2.61 + 0.0087 * dim
        else if( dim >= 500 .and. dim < 1000 ) then
           time_fit = -8.46 + 0.0171 * dim
        else if( dim >= 1000 .and. dim < 2000 ) then
           time_fit = -27.13 + 0.0357 * dim
        else if( dim >= 2000 .and. dim < 3000 ) then
           time_fit = -127.85 + 0.0861 * dim
        else
           time_fit = 0.
        end if
     else if( procs == 16 ) then
        if( dim >= 500 .and. dim < 700 ) then
           time_fit = -2.07 + 0.0073 * dim
        else if( dim >= 700 .and. dim < 1000 ) then
           time_fit = -3.3033 + 0.0159 * dim
        else if( dim >= 1000 .and. dim < 2000 ) then
           time_fit = -17.8 + 0.0256 * dim
        else if( dim >= 2000 .and. dim < 3000 ) then
           time_fit = -79.92 + 0.0567 * dim
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

     if( dim >= 500 .and. dim < 750  ) then
        k_opt = 9
     else if( dim >= 750 ) then
        k_opt = 16
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

     if( procs == 1 .or. procs == 3 .or. procs == 4 .or. procs == 6 &
         .or. procs == 8 .or. procs == 9 .or. procs == 16 ) then
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

      if( dim >= 500 .and. dim < 700 ) then
         k0 = 4
      else if( dim >= 700 ) then
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
     integer(kind=IK) function NROW( procs )

     implicit none

     integer(kind=IK), intent(in) :: procs
     if( procs == 3 .or. procs == 6 .or. procs == 9 ) then
          NROW = 3
     else if( procs == 4 ) then
          NROW = 2
     else if( procs == 8 .or. procs == 16 ) then
          NROW = 4
     else if( procs == 1 ) then
     else if( procs == 1 ) then
          NROW = 1
     end if

     end function NROW
!============================================================================
!============================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!============================================================================
!============================================================================
     integer(kind=IK) function NCOL( procs )

     implicit none

     integer(kind=IK), intent(in) :: procs
     if( procs == 3 ) then
          NCOL = 1
     else if( procs == 4 .or. procs == 6 .or. procs == 8 ) then
          NCOL = 2
     else if( procs == 9 ) then
          NCOL = 3
     else if( procs == 16 ) then
          NCOL = 4
     else if( procs == 1 ) then
          NCOL = 1
     end if

     end function NCOL
!============================================================================
!============================================================================

  end module hitachi_part


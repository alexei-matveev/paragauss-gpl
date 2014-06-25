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
!===============================================================
! Public interface of module
!===============================================================
module gamma_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: calculation of the uncomplete gamma function
  !
  !  Remark of MF:
  !  More precise: gamma(x,j) = int_0^1 t^2(j-1) exp(-x t^2)
  !  i.e. directly the formula for I_(j-1)[x] used in the
  !  theses of Goerling, Belling etc., it is in fact NOT the 
  !  gamma function (e.g. from Abramowich), but of course related
  !  to it (see also remark in Goerlings thesis).
  !
  !  Called by: All routines calculating primitive integrals        
  !
  !  Reference: old lcgto,
  !             McMurchie and Davidson,J. Comp. Phys. Vol. 26 (1978)
  !
  !    MS 3.5.1996
  !
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !----------------------------------------------------------------
  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module
  implicit none
  private
  save
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  public :: gamma, gamma_setup, gamma_close, gamma_is_closed


!================================================================
! End of public interface of module
!================================================================


  real(kind=r8_kind),allocatable :: fgam(:,:), i2_real(:), inv_jx_real(:)
  integer(kind=i4_kind), parameter  :: lup=6
  integer(kind=i4_kind) :: numj = -1 ! was 17, negative means not initialized
                                     ! set up to jmax + 1 in gamma_setup(jmax)
  real(kind=r8_kind) :: inv_i(lup)

contains

  !******************************************************************

  subroutine gamma_close()
    implicit none
    !** End of interface ***************************************

    deallocate(fgam, i2_real, inv_jx_real)

    numj = -1
  end subroutine gamma_close

  !******************************************************************

  subroutine gamma_setup(jmax)
    !
    !     THIS SUBROUTINE CALCULATES A TABLE OF VALUES
    !     OF THE INCOMPLETE GAMMA FUNKTION F (T)
    !                                       J
    !
    !     THESE VALUES ARE TABULATED FOR THE VALUES OF
    !     J = JMIN, JMIN+1, ..., JMIN+NUMJ-1; AND FOR VALUES OF
    !     F= XSTRT, XSTRT+XINT, ..., XSTRT+(NPNT-1)*XINT
    !
    !     NMAX EQUALS THE MAXIMUM NUMBER OF TERMS ALLOWED IN THE
    !     ASYMPTOTIC EXPANSION FOR F    (T)
    !                               JMAX
    !
    !     NUMJ MUST BE 7 MORE THAN MAX Q:
    !
    implicit none
    integer(i4_kind), intent(in) :: jmax
    !** End of interface ***************************************

    integer(kind=i4_kind),parameter :: npnt=121,nmax=100,jmin=0
    real(kind=r8_kind),parameter :: xstrt=0.0_r8_kind,xint=0.1_r8_kind,&
         precis=1.0e-16_r8_kind
    real(kind=r8_kind)  :: term(nmax), nless, jmax2, r, r2, ex, sums
    integer(kind=i4_kind) :: jk,i,j,i2,jx,j2

    ! maybe we should gamma_close() right here if numj >= 0:
    ASSERT(numj==-1)

    ! set module global variable, also indicates
    ! that the module is initialized:
    numj = jmax - jmin + 1 ! == jmax + 1, since jmin == 0

    nless = nmax - 1

    allocate(fgam(numj,npnt), i2_real(2:numj), inv_jx_real(numj-1))

    ! LOOP OVER INDIVIDUAL POINTS

    do i=1,npnt
       jmax2 = 2*jmax + 1
       r = xint * real(I-1,r8_kind) + xstrt
       r2 = 2.0_r8_kind * r
       ex = exp(-r)
       term(1) = 1.0_r8_kind / real(jmax2,r8_kind)

       ! calculate terms in asymptotoc expansion
       jk = 1
       do
          jmax2 = jmax2 + 2
          term(jk+1) = term(jk) * r2/ real(jmax2,r8_kind)
          jk = jk + 1
          if (term(jk)<precis.or.jk>nless) exit
       end do

       ! sum terms in ascending (reverse) order

       sums = 0.0_r8_kind
       do j=jk,1,-1
          sums = sums + term(j)
       enddo
       fgam(numj,i) = sums * ex
       !use downward recursion to claculate other values of fgam(*,i)
       jmax2 = 2*jmax - 1
       jk = numj - 1
       do j=jk,1,-1
          fgam(j,i) = (r2 * fgam(j+1,i) + ex)/real(jmax2,r8_kind)
          jmax2 = jmax2 - 2
       end do
    end do

    i2 = -1
    do  i=2,numj
       i2 = i2 + 2
       i2_real(i) = real(i2,r8_kind)
    end do

    jx = 2*numj - 3
    do j2=numj-1,1,-1
       inv_jx_real(j2) = 1.0_r8_kind / real(JX,r8_kind)
       jx = jx - 2
    end do

    do  i=lup,1,-1
       inv_i(i) = 1.0_r8_kind / real(i,r8_kind)
    end do
  end subroutine gamma_setup


  function gamma(j,t) result(f)
    !     THIS ROUTINE EVALUATES THE INCOMPLETE GAMMA FUNCTION
    !     INTEGRALS F (T) THROUGH F (T).
    !                0             J
    integer(kind=i4_kind),intent(in) :: j ! maximum order of the uncomplete
    ! gamma function plus 1
    real(kind=r8_kind),intent(in)    :: t(:) ! arguments
    real(kind=r8_kind)               :: f(size(t),j) ! results
    !** End of interface ***************************************

    ASSERT(numj>=0)

    call igamma(j,t,f)
  end function gamma

  subroutine igamma(j,t,f)
    !     THIS ROUTINE EVALUATES THE INCOMPLETE GAMMA FUNCTION
    !     INTEGRALS F (T) THROUGH F (T).
    !                0             J
    integer(kind=i4_kind),intent(in) :: j ! maximum order of the uncomplete
    ! gamma function plus 1
    real(kind=r8_kind),intent(in)    :: t(:) ! arguments
    real(kind=r8_kind),intent(out)   :: f(:,:) ! (size(t),j) ! results
    !** End of interface ***************************************
    !     LUP LESS THAN OR EQUAL TO 6:
    real(kind=r8_kind),parameter ::&
         pihlf= 0.1772453850905516e+01_r8_kind, &
         pihlf_half = pihlf * 0.5_r8_kind
    real(kind=r8_kind) :: tdif, fj, g, x, expt(size(t))
    integer(kind=i4_kind) :: ll, i, it, j2, sizet, idintt, igo
#ifdef _VECTOR
! has been !FPP:VPP! commented[[
    real(kind=r8_kind), allocatable :: x_array(:), fj_array(:),  &
         tdif_array(:), exptp(:), sqrt_array(:)
    logical, allocatable :: mask(:)
    integer(kind=i4_kind), allocatable :: it_array(:)
    integer(kind=i4_kind) :: n_x, i_x
! ]]
#endif

    ASSERT(size(t)==size(f,1))
    ASSERT(     j <=size(f,2))

    if ( j+lup .gt. numj ) then
       print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
       print *, "XXXXXXXXX G A M M A  M O D U L E XXXXXXXXXXX"
       print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
       print *, "XXXXXXXXXX W A R N I N G ! ! ! XXXXXXXXXXXXX"
       print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
       print *, "XXXXXXX Derivative order exceeded! XXXXXXXXX"
       print *, "XXXX Provide a better JMAX estimete to XXXXX"
       print *, "XXXXXXXXXXX gamma_setup(JMAX) XXXXXXXXXXXXXX"
       print *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
       print *, "numj=", numj
       print *, "j, lup, j + lup =", j, lup, j + lup

       call gamma_close()
       ! WAS: numj = j + lup
       call gamma_setup(j + lup - 1)
    endif

    expt = exp(-t)

    sizet = size(t)

#ifdef _VECTOR
! has been !FPP:VPP! commented[[
    if (sizet .ge. 32) then
       ! vectorised code
       allocate(x_array(sizet), mask(sizet), sqrt_array(sizet))

       ! t > 12
       mask = t >= 12.0_r8_kind
       where (mask)
          x_array = 1.0_r8_kind / t
       elsewhere
          x_array = 1.0_r8_kind
       endwhere
       sqrt_array = sqrt(x_array)
       f(:,1) = pihlf_half * sqrt_array - expt * x_array * &
            ( 0.4999489092_r8_kind + x_array * (-0.2473631686_r8_kind + &
            x_array * (0.321180909_r8_kind - x_array * 0.3811559346_r8_kind)) )
       ! USE UPWARD RECURSION TO CALCULATE THE OTHER VALUES FOR F(*)
       x_array = x_array * 0.5_r8_kind
       do  i=2,j
             f(:,i) = x_array * (i2_real(i) * f(:,i-1) - expt)
       end do


       ! T < 12
       ! CALCULATE F(T) USING TAYLOR SERIES EXPANSION AND TABULATED VALUES
       mask = t < 12.0_r8_kind
       n_x = count(mask)
       if ( n_x .ge. 128 ) then
          x_array = pack(t,mask)
          allocate(tdif_array(n_x), fj_array(n_x), it_array(n_x), exptp(n_x))
          exptp = exp(-x_array(1:n_x))
          it_array = int(x_array(1:n_x) * 10.0_r8_kind, i4_kind) + 1
          tdif_array = real(it_array - 1, r8_kind) / 10.0_r8_kind - x_array(1:n_x)
          do i_x = 1, n_x
             fj_array(i_x) = fgam(j+lup,it_array(i_x))
          enddo
          do  i=lup,1,-1
             do i_x = 1, n_x
                fj_array(i_x) = fj_array(i_x) * tdif_array(i_x) * inv_i(i) + &
                     fgam(j-1+i,it_array(i_x))
             enddo
          enddo
          f(:,j) = unpack(fj_array, mask, f(:,j))
          !  USE DOWNWARD RECURSION TO CALCULATE THE NEEDED VALUES
          !  OF F(I), I= 1, J
          x_array(1:n_x) = 2.0_r8_kind * x_array(1:n_x)
          do j2=j-1,1,-1
             fj_array = (x_array(1:n_x) * fj_array + exptp) * inv_jx_real(j2)
             f(:,j2) = unpack(fj_array, mask, f(:,j2))
          end do
          deallocate(tdif_array, fj_array, it_array, exptp)
       elseif ( n_x .gt. 0 ) then
          do ll = 1, sizet
             if ( mask(ll) ) then
                it = int( t(ll) * 10.0_r8_kind ,i4_kind) + 1
                tdif = real(it - 1, r8_kind) / 10.0_r8_kind - t(ll)
                fj = fgam(j+lup,it)
                do  i=lup,1,-1
                   fj = fj * tdif * inv_i(i) + fgam(j-1+i,it)
                end do
                f(ll,j) = fj
                !  USE DOWNWARD RECURSION TO CALCULATE THE NEEDED VALUES
                !  OF F(I), I= 1, J
                x = 2.0_r8_kind * t(ll)
                do j2=j-1,1,-1
                   f(ll,j2) = (x * f(ll,j2+1) + expt(ll)) * inv_jx_real(j2)
                end do
             endif
          enddo
       endif

       deallocate(x_array, mask, sqrt_array)

    else
! ]]
#endif

       ! optimized scalar code

       elements: do ll=1,sizet
          if ( t(ll) <= 1.0e-14_r8_kind) then
             !  if t = 0, just define f(*) from fgam(*,1)
             do i=1,J
                f(ll,i) = fgam(i,1)
             end do
          else if (t(ll) >= 12.0_r8_kind) then
             ! T  >= 12
             x = 1.0_r8_kind / t(ll)
             if (t(ll) < 30.0_r8_kind) then
                idintt = int (t(ll))
                igo = (idintt - 9) / 3
                select case(igo)
                case (1)
                   g = 0.4999489092_r8_kind + x * (-0.2473631686_r8_kind + &
                        x * (0.321180909_r8_kind - x * 0.3811559346_r8_kind))
                case (2)
                   g = 0.4998436875_r8_kind + x * (-0.24249438_r8_kind + &
                        x * 0.24642845_r8_kind)
                case (3)
                   g = 0.499093162_r8_kind - x * 0.2152832_r8_kind
                case (4)
                   g = 0.499093162_r8_kind - x * 0.2152832_r8_kind
                case (5)
                   g = 0.490_r8_kind
                case (6)
                   g = 0.490_r8_kind
                case default
                   g = 0.0 ! to make the compiler happy
                   ABORT('no such case')
                end select
                f(ll,1) = pihlf_half * sqrt(x) - expt(ll) * g * x
             else
                f(ll,1) = pihlf_half * sqrt(x)
             end if
             ! USE UPWARD RECURSION TO CALCULATE THE OTHER VALUES FOR F(*)
             x = x * 0.5_r8_kind
             do  i=2,j
                f(ll,i) = x * (i2_real(i) * f(ll,i-1) - expt(ll))
             end do
          else
             ! T < 12
             ! CALCULATE F(T) USING TAYLOR SERIES EXPANSION AND TABULATED VALUES
             it = int( t(ll) * 10.0_r8_kind ,i4_kind) + 1
             tdif = real(it - 1, r8_kind) / 10.0_r8_kind - t(ll)
             fj = fgam(j+lup,it)
             do  i=lup,1,-1
                fj = fj * tdif * inv_i(i) + fgam(j-1+i,it)
             end do
             f(ll,j) = fj
             !  USE DOWNWARD RECURSION TO CALCULATE THE NEEDED VALUES
             !  OF F(I), I= 1, J
             x = 2.0_r8_kind * t(ll)
             do j2=j-1,1,-1
                f(ll,j2) = (x * f(ll,j2+1) + expt(ll)) * inv_jx_real(j2)
             end do
          endif
       end do elements

#ifdef _VECTOR
! has been !FPP:VPP! commented[[
   endif
! ]]
#endif

  end subroutine igamma

  function gamma_is_closed()
    implicit none
    logical :: gamma_is_closed
    ! *** end of interface ***

    gamma_is_closed = (numj >= 0)
  end function gamma_is_closed

end module gamma_module

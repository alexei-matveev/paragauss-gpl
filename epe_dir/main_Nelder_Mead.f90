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
PROGRAM tminim
use NELDER_MEAD 

! Use minim to maximize the objective function:

! 1/{1 + (x-y)^2} + sin(pi.y.z/2) + exp[-{(x+z)/y - 2}^2]

! with respect to x, y and z.
! We will actually minimize its negative.
! The maximum occurs at x = y = z = +/-sqrt(4n+1) for any integer n.

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp)          :: object, p(3), simp, step(3), stopcr, var(3)
INTEGER            :: ier, iprint, iquad, maxf, nloop, nop
LOGICAL            :: first
external objfun,pgobjfun

! Set up starting values & step sizes.
! Parameters p(1), p(2), p(3) are x, y & z.

open(11,file='mim_res')
p(1) = +229.8
p(2) = 0.3378
p(3) = 20.18

step = 1.2D0
step(2)=0.005

nop = 3

! Set max. no. of function evaluations = 50, print every 10.

maxf = 70
iprint = 1

! Set value for stopping criterion.   Stopping occurs when the
! standard deviation of the values of the objective function at
! the points of the current simplex < stopcr.

stopcr = 1.d-05
nloop = 6

! Fit a quadratic surface to be sure a minimum has been found.

iquad = 1

! As function value is being evaluated in REAL (dp), it
! should be accurate to about 15 decimals.   If we set simp = 1.d-6,
! we should get about 9 dec. digits accuracy in fitting the surface.

simp = 1.d-7

! Now call MINIM to do the work.

first = .true.
DO
  CALL minim(p, step, nop, object, maxf, iprint, stopcr, nloop,   &
             iquad, simp, var, pgobjfun, ier)

! If ier > 0, try a few more function evaluations.

  IF (ier == 0) EXIT
  IF (.NOT. first) STOP
  first = .false.
  maxf = 100
END DO

! Successful termination.

WRITE(*, 900) object, p
900 FORMAT(' Success !'/' Objective function = ', f12.6/ ' at: ', 6f12.6)
WRITE(*, 910) var
910 FORMAT(' Elements of var = ', 4f12.6)


contains
 subroutine make_simcv_script
	print*,'make_simcv_script'
 end subroutine make_simcv_script
END PROGRAM tminim



SUBROUTINE objfun(p, func)

! This is the subroutine which the user must write.
! Remember that we are minimizing the negative of the function we
! really want to maximize.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp), INTENT(IN)  :: p(:)
REAL (dp), INTENT(OUT) :: func

!     Local variables

REAL (dp), PARAMETER :: half = 0.5D0, one = 1.D0, pi = 3.14159265357989D0,  &
                        two = 2.D0
REAL (dp)            :: x, y, z

x = p(1)
y = p(2)
z = p(3)
func = -one/(one + (x-y)**2) - SIN(half*pi*y*z) - EXP(-((x+z)/y - two)**2)
RETURN
END SUBROUTINE objfun

SUBROUTINE pgobjfun(p, func)

! This is the subroutine which the user must write.
! Remember that we are minimizing the negative of the function we
! really want to maximize.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp), INTENT(IN)  :: p(:)
REAL (dp), INTENT(OUT) :: func

!     Local variables
integer*4::i,npoints=7
REAL (dp):: x(7),curv(7),curv_c(7),e
character(len=14) psymv(3)
write(psymv(1),'(1PE14.6)') p(1)
write(psymv(2),'(1PE14.6)') p(2)
write(psymv(3),'(1PE14.6)') p(3)
print*,psymv 
call system('sed "s/-130.000/'//psymv(1)//'/g" epe_simol_parameters_base'//&
     '|sed "s/0.3012000/'//psymv(2)//'/g"'//&
     '|sed "s/20.11000/'//psymv(3)//'/g">epe_simol_parameters')
!call system('./bcurvmSiH4_F')
call system('./bsim')
open (10,file='prelim_curve')
open (12,file='curv_comp')
print *,'contents of curv_comp'
do i=1,npoints
 read(10,*) x(i)
 read(10,*) e     
 write(12,*) x(i),e
 print*,x(i),e
enddo
close (10)
close (12)
open (10,file='curv_o10mg26')
print*,'reference curve'
do i=1,npoints
read(10,*) x(i),curv(i)
print *,x(i),curv(i)
end do
curv=curv-curv(4)
close (10)
open (10,file='curv_comp')
print*,'variable corve'
do i=1,npoints
read(10,*) x(i),curv_c(i)
end do
curv_c=curv_c-curv_c(4)
do i=1,npoints
print*, curv_c(i)/4.0,curv(i)
end do
print*, 'var curve'
do i=1,npoints
print*,x(i), curv_c(i)/4.0
end do
print*, 'ref curve'
do i=1,npoints
print*,x(i),curv(i)
end do
func=   sum((curv_c(1:npoints)/4.d0-curv(1:npoints))**2)
!func=   sum((curv_c-curv)**2)

print*,'func1',func
close (10)

RETURN
END SUBROUTINE pgobjfun

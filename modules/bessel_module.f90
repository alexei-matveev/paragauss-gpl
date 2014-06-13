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
!===============================================================
! Public interface of module
!===============================================================
module bessel_module
  !-------------- Module specification ---------------------------
  !
  !  Purpose: calculation of the modified spherical bessel integrals
  !
  !  Called by:  routines calculating primitive integrals of 
  !              pseudopotentials       
  !
  !  Reference: 
  !             McMurchie and Davidson,J. Comp. Phys. Vol. 26 (1978)
  !
  !    A.Hu   on 25th, Aug, 1998
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
  public :: bessel_1,bessel_setup,bessel_close,integral_radial
  public :: bessel2

!================================================================
! End of public interface of module
!================================================================


!  real(kind=r8_kind)             :: alpha_1,rk_1
!  real(kind=r8_kind)             :: t_1
!  real(kind=r8_kind)             :: argab,expab
  real(kind=r8_kind)             :: fac(1:17),fprod(1:9,1:9),ddfac(1:17)

  integer(i4_kind), parameter    :: low=5,medium=10,high=20
  !nteger(i4_kind), parameter    :: low=3,medium= 6,high=12
  !nteger(i4_kind), parameter    :: low=2,medium= 4,high= 8
  !nteger(i4_kind), parameter    :: low=1,medium= 2,high= 4
  ! Orders of numerical Gauss-Hermite quadrature.

  ! The orders we use, seem to be too high: for Pd19(Stutt)
  ! I didnt observe ANY change for all orders listed above.
  ! For the problematic F2(Stutt) case I get same
  !      10^-12 relative accuracy in total energy and
  !    5*10^ -8 relative accuracy in gradients
  ! for the highest and lowest quadrature orders from the list above.
  ! TODO: dynamic damping of the order would be great

  real(kind=r8_kind)             :: ptpow(1:high,1:13)   ! (high,1+la+lb)
  real(kind=r8_kind)             :: abess(1:high,1:7), & ! (high,1+la   ) for bess(lm)..bess(lm+la)
                                    bbess(1:high,1:7)    ! (high,1   +lb) for bess(lm)..bess(lm+lb)
  integer(kind=i4_kind)          :: ilo,ihi

  integer(i4_kind), parameter    :: NDFAC = 64
  real(kind=r8_kind)             :: dfac(NDFAC)

! real(kind=r8_kind),parameter,dimension(1:9) :: tmin=(/  &
!      31.0_r8_kind, 28.0_r8_kind, 25.0_r8_kind, 23.0_r8_kind, 22.0_r8_kind, &
!      20.0_r8_kind, 19.0_r8_kind, 18.0_r8_kind, 15.0_r8_kind /)
  real(kind=r8_kind),dimension(1:9) :: tmin=(/  &
       31.0_r8_kind, 28.0_r8_kind, 25.0_r8_kind, 23.0_r8_kind, 22.0_r8_kind, &
       20.0_r8_kind, 19.0_r8_kind, 18.0_r8_kind, 15.0_r8_kind /)
  integer(kind=i4_kind), parameter  ::  nzero = 0_i4_kind
  real(kind=r8_kind),parameter   :: zero    = 0.0_r8_kind
  real(kind=r8_kind),parameter   :: half    = 0.5_r8_kind
  real(kind=r8_kind),parameter   :: one     = 1.0_r8_kind
  real(kind=r8_kind),parameter   :: two     = 2.0_r8_kind
  real(kind=r8_kind),parameter   :: three   = 3.0_r8_kind
  real(kind=r8_kind),parameter   :: pi      = 3.14159265358979324_r8_kind
  real(kind=r8_kind),parameter   :: four    = 4.0_r8_kind
  real(kind=r8_kind),parameter   :: five    = 5.0_r8_kind
  real(kind=r8_kind),parameter   :: sqpi    = 1.77245385090551603_r8_kind
  real(kind=r8_kind),parameter   :: fpi     = 12.56637061435917295_r8_kind
  real(kind=r8_kind),parameter   :: sqpi2   = 1.25331413731550025_r8_kind
  real(kind=r8_kind),parameter   :: hundred = 100.0_r8_kind

  integer(i4_kind), private :: REFCOUNT = 0

contains

  !******************************************************************

  subroutine bessel_close()
!    !** End of interface ***************************************
!    integer(kind=i4_kind)          :: alloc_stat
!    deallocate(dfac,stat=alloc_stat)
!     if(alloc_stat.ne.0) call error_handler &
!       ("bessel_close: dfac deallocation failed ")

    ASSERT(REFCOUNT>0)
    REFCOUNT = REFCOUNT - 1
  end subroutine bessel_close

  !******************************************************************

  subroutine ptprep(nj,lm,la,lb,rka,rkb,alpha_1)
    use gaussq_data, only: LOC,ZW
    implicit none
    integer(i4_kind)  , intent(in) :: nj,lm,la,lb
    real(kind=r8_kind), intent(in) :: rka,rkb,alpha_1
    ! *** end of interface ***
    !  set up the points and weights for the calculating redial integrals
    real(kind=r8_kind), parameter  :: const1 = 1.0e+5_r8_kind
    real(kind=r8_kind), parameter  :: const2 = 1.0e+3_r8_kind
    real(kind=r8_kind)             :: c, sqalp
    integer(kind=i4_kind)          :: i, L, P
    real(kind=r8_kind)             :: argsum
    integer(i4_kind)               :: N
    real(r8_kind)                  :: R(high) ! only (:N) is used
    !-- End of interface ---------------------------------------------------

    sqalp=sqrt(alpha_1)
    c=(rka+rkb)/(two*sqalp)
    argsum = ((rka+rkb)**2)/(two*alpha_1)

    ! select the order of the quadrature:
    if(argsum.le.const1) then    ! <= 1.0e+5
       if(argsum.gt.const2) then ! >  1.0e+3
         N = medium
       else
          N = high
       end if
    else
       N = low
    end if
    ! location of the zeros/weights inside of an array:
    ilo = LOC(N)
    ihi = ilo + N - 1

    ASSERT(la<=6)
    ASSERT(lb<=6)

    ! abscissas for the quadrature:
    R(:N) = ( c + ZW(ilo:ihi,1) ) / sqalp

    ! FIXME: so far bess() retruns zero for negative args
    !        but if that changes ...

    ! tabulate the values of bess(ka*r)
    ! WARNING: Bessel of order lm+L is stored at abess(:,1+L)
    FORALL( i=1:N, L=0:la )
       abess(i,1+L) = bess(rka*R(i),lm+L)
    END FORALL

    ! tabulate the values of bess(kb*r)
    ! WARNING: Bessel of order lm+L is stored at abess(:,1+L)
    FORALL( i=1:N, L=0:lb )
       bbess(i,1+L) = bess(rkb*R(i),lm+L)
    END FORALL

    ! compute necessary powers R^(nj+2*lm+2*la+2*lb), la=0,LA, lb=0,LB
    ptpow(:N,1) = R(:N)**(nj+2*lm)/sqalp ! FIXME: common factor ``scalp'' here!
    R(:N) = R(:N)**2
    do P=1,la+lb
       ptpow(:N,1+P) = ptpow(:N,P) * R(:N)
    end do

  end subroutine ptprep
 !********************************************************************************

  function ptwt(n,lm,la,lb) result(f)
  ! computing the redial integrals Q(n,la,lb)
    use gaussq_data, only: ZW
    implicit none
  !-------------------------------------------------------------------------------
  integer(kind=i4_kind), intent(in) :: n
  integer(kind=i4_kind), intent(in) :: lm
  integer(kind=i4_kind), intent(in) :: la
  integer(kind=i4_kind), intent(in) :: lb
  real(kind=r8_kind)                :: f
  !--  end of interface ----------------------------------------------------------
  integer(kind=i4_kind)             :: i, i_1 
  real(kind=r8_kind)                :: func
  !--  end of local variables ----------------------------------------------------

  f = zero
  do i= ilo, ihi
     i_1=i-ilo+1
     ! ptpow at position (1+la+lb) contains R^(n+2*lm + 2*la + 2*lb)
     ! i.e. all powers together including those of bess(lm+la) and bess(lm+lb):
     ! WARNING: Bessel of order lm+L is stored at abess(:,1+L)
     func=ptpow(i_1,1+la+lb)*abess(i_1,1+la)*bbess(i_1,1+lb)
     f = f + ZW(i,2)*func
  end do
  f = f / ( DFAC(2*(lm+la)+3) * DFAC(2*(lm+lb)+3))
  ! factors due to use of ``normalized'' bessel functions

 end function ptwt
 !********************************************************************************
  function qpasy(n,la1,lb1,alp,xka1,xkb1,iflag) result(f)
  ! partially asymptotic form for q(n,la1,lb1)
  ! resutl includes a factor, exp(-xka*xkb/(4.0*alpha))
  ! to help prevent overflow
  !-------------------------------------------------------------------------------
  integer(kind=i4_kind), intent(in)   :: n ! the order of function q(n,la1,lb1)
  integer(kind=i4_kind), intent(in)   :: la1 ! the order of bessel function
  integer(kind=i4_kind), intent(in)   :: lb1 ! the order of bessel function
  real(kind=r8_kind)   , intent(in)   :: alp ! alpha+beta+gamma
  real(kind=r8_kind)   , intent(in)   :: xka1 ! 2*alpha|c-a|
  real(kind=r8_kind)   , intent(in)   :: xkb1 ! 2*beat |c-b|
  integer(kind=i4_kind), intent(in)   :: iflag ! flag for the size of xka1 and xkb1
  real(kind=r8_kind)                  :: f
  !--- end of interface ----------------------------------------------------------
  real(kind=r8_kind)                  :: xka,xkb,tk,qold1,sum,qnew,term,coe,qold2&
                                         , perfac
  integer(kind=i4_kind)               :: la,lb,j,xxxx,yyyy,nprime
  real(kind=r8_kind), parameter       :: const1=1.0e-13_r8_kind
  real(kind=r8_kind)                  :: t_1
  real(kind=r8_kind)                  :: rk_1,alpha_1
  !--- end of local variables -----------------------------------------------------
  !  first set up xkb as largest

  if (iflag.ne.3) then
     xka=xka1
     xkb=xkb1
     la=la1
     lb=lb1
  else
     xka=xkb1
     xkb=xka1
     la=lb1
     lb=la1
  end if
  ! to set up parameters for qcomp using xkb
  alpha_1 = one
  rk_1=xkb/sqrt(alp)
  t_1=rk_1*rk_1/four
  ! now run power serise using xka, obtaining initial q(n,l)
  ! from qcomp then recurring upwards.
  tk=xka*xka/(two*alp)
  ! j= 0 term in sum
  qold1=qcomp(n+la,lb,rk_1,alpha_1,t_1) / sqrt(alp)**lb ! FIXME!
  sum=qold1/dfac(la+la+3)

  if (tk.ne.zero) then
    ! j=1 term in sum
    nprime=n+la+2
    qnew=qcomp(nprime,lb,rk_1,alpha_1,t_1) / sqrt(alp)**lb ! FIXME!
    term=qnew*tk/dfac(la+la+5)
    sum=sum+term
    ! j=2 term in sum
    j=2
    nprime=n+la+j+j
    coe=one/dfac(la+la+3)
    qold2=qold1*coe
    qold1=qnew*coe
    xxxx = nprime+nprime-5
    yyyy = (lb-nprime+4)*(lb+nprime-3)              
    qnew = (t_1+real(xxxx,r8_kind)/two)*qold1 +&
           real(yyyy,r8_kind)*qold2/four
    xxxx = (j-1)*(la+la+j+j-1)
    yyyy = j*(la+la+j+j+1)             
    term = qnew*tk*tk/(real(xxxx,r8_kind)*&
           real(yyyy,r8_kind))
    sum=sum+term
    if (abs(term/sum).gt.const1) then
       ! increment j for next term
       do_block_loop:   do
          j=j+1
          nprime=n+la+j+j
          ! compute j-2 factor to be absorbed into qold1 and qold2
          xxxx = (j-2)*(la+la+j+j-3)
          coe  = tk/real(xxxx,r8_kind)
          qold2=qold1*coe
          qold1=qnew*coe
          xxxx = nprime+nprime-5
          yyyy = (lb-nprime+4)*(lb+nprime-3)        
          qnew = (t_1+real(xxxx,r8_kind)/two)*qold1 +&
                  real(yyyy,r8_kind)*qold2/four
          ! now correct qnew for j-1 and j factors
          xxxx = (j-1)*(la+la+j+j-1)
          yyyy = j*(la+la+j+j+1)       
          term = qnew*tk*tk/(real(xxxx,r8_kind)*&
                 real(yyyy,r8_kind)) 
          sum=sum+term
          if (abs(term/sum).lt.const1) exit do_block_loop
          if(j.gt.100) then
             write(1,*)' in qpasy j is greater 100'
             stop
          end if
          cycle do_block_loop
        end do do_block_loop
    end if ! abs(term/sum)<=const1
  end if ! tk/=0
!!$  if (la.eq.0) perfac=one/sqrt(alp**(n+la+1))
!!$  if (la.ne.0) perfac=(xka**la)/sqrt(alp**(n+la+1))
  perfac=one/sqrt(alp**(n+la+1))
  f=perfac*sum
 end function qpasy

 !*********************************************************************
  PURE ELEMENTAL function bess(z,l) result(f)
  ! calculates the normalized form
  !      M (z)
  !       l
  ! of the modified shperical Bessel function of the first kind
  !
  !               z^l * exp(z)
  !      I (z) = ------------- * M (z)
  !       l        ( 2l+1 )!!     l
  !
  !-------------------------------------------------------------------------------
  integer(kind=i4_kind), intent(in) :: l ! the order of Bessel function
  real(kind=r8_kind)   , intent(in) :: z ! the variable of Bessel function
  !--- end of interface ----------------------------------------------------------
  real(kind=r8_kind), parameter     :: epsx = 5.0e-14_r8_kind
  real(kind=r8_kind)                :: term, rp, rm, zp
  integer(kind=i4_kind)             :: j,l1,k
  real(kind=r8_kind)                :: f
  !--- end of local variables ----------------------------------------------------

      IF (Z.LT.ZERO) THEN
         f = ZERO
      ELSE IF(Z.EQ.ZERO) THEN
         IF(L.NE.0) THEN
            f = ZERO
         ELSE
            f = ONE
         END IF
      ELSE IF((Z.GT.ZERO).AND.(Z.LE.FIVE)) THEN
         ZP   = Z*Z/TWO
         TERM = ONE ! BEFORE: (Z**L)/DFAC(L+L+3)
         f = TERM
         J    = 0
    do k=1,16 ! it takes 15 iter for Z=5.0
         J = J+1
         TERM   = TERM*ZP/(J*(L+L+J+J+1))
         f   = f +TERM
         ! all terms are positive:
         IF ( TERM/f .LE. EPSX) exit
    enddo
         f=f*EXP(-Z)
!!$         if(J.ge.16) then
!!$         write(1,*)' cycle in bess J is greater 16'
!!$         stop
!!$         end if
      ELSE IF((Z.GT.FIVE).AND.(Z.LE.16.1_r8_kind)) THEN
         ZP = HALF/Z
         L1  = L+1
         RP = FPROD(L1,L1)
         RM = FPROD(L1,L1)
         DO  K=L,1,-1
            RP=FPROD(K,L1) + ZP*RP
            RM=FPROD(K,L1) - ZP*RM
         END DO
         f=( RM - ((-1)**L)*RP*EXP(-TWO*Z) )*ZP * DFAC(L+L+3)/Z**L
         ! BEFORE: no DFAC(L+L+3)/Z**L
      ELSE
         ZP = HALF/Z
         L1  = L+1
         RM  = FPROD(L1,L1)
         DO K=L,1,-1
            RM=FPROD(K,L1) - ZP*RM
         END DO
         f=RM*ZP * DFAC(L+L+3)/Z**L
         ! BEFORE: no DFAC(L+L+3)/Z**L
      END IF
 end function bess
 !**********************************************************************

  subroutine bessel_setup( )
    !
    !     This routine controls the computation of Q(n,l) ...
    !     arguments are alpha, xk, and tx = xk**2 / (4.0*alpha).
    !     there are no restrictions on the magnitude of tx.
    !     In order to prevent overflow, the result return in bessel_1
    !     is actually exp( - tx ) * Q(n,l).
    !     qalt and qpow are restricted to (n+l).le.26, (2l).le.24,
    !     but this may be increased by simply increasing the dfac array.
    !     qasy is valid for arbitrary ( n,l )
    !
    
    integer(kind=i4_kind) :: i,j,k,xxxi
    real(kind=r8_kind)    :: xi2

    if( REFCOUNT > 0 ) return
    REFCOUNT = REFCOUNT + 1

    dfac(1) = one
    dfac(2) = one
    dfac(3) = one

    do i = 4, NDFAC
         xi2 = real(i,r8_kind) - two
         dfac ( i ) = dfac ( i - 2 ) * xi2
    enddo

    fac(1) = one
    do i = 1, 16
       xxxi = i
       fac(i+1) = fac(i) * real(xxxi,r8_kind)
    enddo

    do i = 1, 9
       do j = 1, i
          k = j - 1
          fprod(j,i) = fac(i+k) / (fac(j)*fac(i-j+1))
       enddo
    enddo

    ddfac(1)=one
    do i=2,17
       ddfac(i)=ddfac(i-1)*real(2*i-1,r8_kind)
    enddo

  end subroutine bessel_setup

 !**************************************************************************
  function bessel_1(nj,lj_in,eta,rk,abc_fact,check_abc) result(f)
    ! This subroutine is used to evluate the integrals Q^(nj)_(lj)[eta,k] of
    ! Bessel function combined with algebric and exponential 
    ! functions. The integrals deal with the pseudopotential
    ! with the maximum angular momentum plus one
    integer(kind=i4_kind),intent(in) :: nj     ! the symmetry of basis set
    integer(kind=i4_kind),intent(in) :: lj_in  ! l1 + l2 + 1
    integer(kind=i4_kind),intent(in) :: check_abc !
    real(kind=r8_kind),intent(in)    :: eta(:) ! alpha(i) + beta(j) + gamma
    real(kind=r8_kind),intent(in)    :: rk(:)  ! scalar k =2|alpha(i)(c-a)+beta(j)*(c-b)|
    real(kind=r8_kind),intent(in)    :: abc_fact(:)
    real(kind=r8_kind)               :: f(size(eta),1:lj_in)

    !xxxxxxxxx end of interface xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
    !xxxxxxxxx local variables  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  
    
    real(kind=r8_kind), allocatable :: q(:,:,:) ! Q(:,N,L)
    integer(kind=i4_kind)           :: lj1,lhi,lr,i,l,xxxx,yyyy,nbeg,&
                                       nmax,lt,alloc_stat, temp, lj
    real(kind=r8_kind), dimension(size(eta))  :: talph, t
    real(kind=r8_kind)               :: talph_1
    real(kind=r8_kind)               :: t_1
    real(kind=r8_kind)               :: rk_1,alpha_1
   !xxxxxxxxx end of local variables xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
   !xxxxxxxx executable codes xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   !call bessel_setup()

   lt = size(eta)
   lj = lj_in - 1

   allocate(q(lt,nj+lj+3,lj+1),stat=alloc_stat)
     if ( alloc_stat .ne. 0) call error_handler &
        ( " Bessel_1: allocation(1) failed" )

    q=0.0_r8_kind
   t = rk * rk / ( four * eta )
  
    if ( check_abc == nzero )  then     ! the centers a, b and c are all the same
      nmax = nj + lj
      q(:,1,1) = sqrt ( pi / eta(1:lt) ) * half  
      q(:,2,1) = one / ( two * eta(1:lt) )    
      do i = 2, nmax   
         temp = i - 1
         q(:,i+1,1) = temp * q(:,i-1,1) / ( two * eta(:) ) 
      enddo

      do i = 1, lj_in
           f(:,i) = exp(abc_fact(:)+t(:))*q(:,i+nj,1)  
      enddo
    elseif(lj==0)then
        do i = 1, lt
          alpha_1   = eta   ( i )
          rk_1      = rk    ( i )
          t_1       = t     ( i )
          q(i,nj+1,1)=qcomp(nj,0,rk_1,alpha_1,t_1)
        end do
        f(:,1)=exp(abc_fact(:)+t(:))*q(:,nj+1,1)     
    elseif(lj==1)then 
        do i = 1, lt
           alpha_1   = eta   ( i )
           rk_1      = rk    ( i )
           t_1       = t     ( i )
           if (rk(i)/=0.0_r8_kind) then
               q(i,nj+1,1)=qcomp(nj,0,rk_1,alpha_1,t_1)
               q(i,nj+2,2)=qcomp(nj+1,1,rk_1,alpha_1,t_1) !/rk(i)
           else
               q(i:i,nj+1,1)=integral_radial(nj,eta(i:i))
               q(i:i,nj+2,2)=integral_radial(nj+2,eta(i:i))/three
           endif
        enddo
        f(:,1)=exp(abc_fact(:)+t(:))*q(:,nj+1,1) 
        f(:,2)=exp(abc_fact(:)+t(:))*q(:,nj+2,2) 
    elseif(lj.ge.2)then
     talph = two * eta
    select case ( nj )                       ! decide nj = 0, 1, 2
      case ( 0 )                       ! in the case of nj = 0
      do i = 1, lt
         alpha_1   = eta   ( i )
         rk_1      = rk    ( i )
         t_1       = t     ( i )
         talph_1   = talph ( i ) 
      if (rk(i)/=0.0_r8_kind) then
         q(i,3,1) = qcomp ( 2, 0,rk_1,alpha_1, t_1 )     
         lj1 = lj -1 
         do l = 1, lj1   
           q(i,l+3,l+1) = q(i,l+2,l) / talph_1
         enddo
         nbeg=4
         if (t_1.le.three ) then               ! DOWNWARDS RECURRENCE FOR T .LE. 3.
            q(i,lj+1,lj+1) = qcomp  ( lj,lj,rk_1,alpha_1,t_1 )  
            lhi            = lj - 1
            do lr = nzero, lhi 
               l    = lhi - lr  
               xxxx = l + l + 1  
               q(i,l+1,l+1) = ( talph_1 * q(i,l+3,l+1) &
                           - rk_1**2 * q(i,l+2,l+2))/real(xxxx,r8_kind)
            enddo
         else                                  !  UPWARDS RECURRENCE FOR T .GT. 3.
            q(i,1,1) = qcomp(0,0,rk_1,alpha_1,t_1)   
            do l = 1, lj     
                xxxx = l + l - 1       
                q(i,l+1,l+1) = ( talph_1 * q(i,l+2,l) & 
                         -real(xxxx,r8_kind) *q(i,l,l))/rk_1**2
            enddo 
         endif
!!$         do l=0,lj
!!$            q(i,l+1,l+1)=q(i,l+1,l+1)/rk(i)**l
!!$         enddo
      else
         do l=0,lj
            q(i:i,l+1,l+1) = integral_radial(nj+2*l,eta(i:i))/ddfac(l+1)
         enddo
      endif
      enddo                            ! end of case nj = 0
      case ( 1 )                       ! beginning of the case nj = 1
      do i = 1, lt
         alpha_1   = eta   ( i )
         rk_1      = rk    ( i )
         t_1       = t     ( i )
         talph_1   = talph ( i )
         if (rk(i)/=0.0_r8_kind) then
             nbeg=3
             if (t_1.le.three ) then        !  DOWNWARDS RECURRENCE FOR T .LE. 3. 
                q(i, lj+2, lj+1) = qcomp( lj+1, lj,rk_1,alpha_1, t_1)   
                q(i, lj+1, lj  ) = qcomp( lj  , lj - 1,rk_1,alpha_1, t_1)
                lhi=lj-2 
                do lr = nzero, lhi    
                    L          = LHI-LR
                    XXXX       = L+L+3
                    YYYY       = L+L+2
                    q(i,l+2,l+1) = (talph_1*rk_1**2 * q(i,l+4,l+3)-&
                                   (rk_1**2-real(xxxx,r8_kind)*talph_1) &
                                   *q(i,l+3,l+2))/real(yyyy,r8_kind)
                enddo
             else                        !  UPWARDS RECURRENCE FOR T .GT. 3.
                q(i,2,1)=qcomp(1,0,rk_1,alpha_1,t_1)
                q(i,3,2)=qcomp(2,1,rk_1,alpha_1,t_1)
                do l = 2, lj                      
                   XXXX       = L+L-2
                   YYYY       = L+L-1
                   q(i,l+2,l+1) = (real(xxxx,r8_kind)/rk_1**2 * q(i,l,l-1)&
                        + (one - real(yyyy,r8_kind)*talph_1/rk_1**2) *  q(i,l+1,l))/talph_1
                enddo
             endif
!!$             do l=0,lj
!!$                q(i,l+2,l+1)=q(i,l+2,l+1)/rk(i)**l
!!$             enddo
         else
            do l=0,lj
               q(i:i,l+2,l+1) = integral_radial(nj+2*l,eta(i:i))/ddfac(l+1)
            enddo
         endif
      enddo                           ! end of the case nj = 1
     case ( 2 )                      ! beginning of the case nj = 2
     do i = 1, lt
          alpha_1   = eta   ( i )
          rk_1      = rk    ( i )
          t_1       = t     ( i )
          talph_1   = talph ( i )
          if (rk(i)/=0.0_r8_kind) then
             q(i,3,1) = qcomp(2,0,rk_1,alpha_1,t_1)
             do l = 1, lj
                q(i,l+3,l+1) = q(i,l+2,l)/talph_1
             enddo
!!$             do l=0,lj
!!$                q(i,l+3,l+1)=q(i,l+3,l+1)/rk(i)**l
!!$             enddo
          else
             do l=0,lj
                q(i:i,l+3,l+1) = integral_radial(nj+2*l,eta(i:i))/ddfac(l+1)
             enddo
          endif
     end do
   case default
     ABORT('bessel_1: n out of range')
   end select 
      do i = 1, lj_in
           f(:,i) = exp(abc_fact(:)+t(:))*q(:,i+nj,i) 
      enddo
  end if                         ! end of lj>2

    deallocate(q, stat = alloc_stat )
    if (alloc_stat /= 0 ) call error_handler &
       ("Bessel_1: deallocation (1) failed ")


  end function bessel_1

  subroutine bessel2(f,factor,nj,lm,la,lb,alpha,beta,gamma,ca,cb,ab,expmax)
   ! This functional subroutine is used to evluate the integrals
   ! Q^(nj)_(lm+la)_)lm+lb)(eta, k) of two Bessel functions
   ! combined with algebric and exponentail
   ! functions. The integrals deal with the pseudopotentials of
   ! semi-local forms
   ! ------------------------------------------------------------------------------
   integer(kind=i4_kind), intent(in)    :: nj  ! the power of pseudopotentials
   integer(kind=i4_kind), intent(in)    :: lm  ! angular quantum number 
   integer(kind=i4_kind), intent(in)    :: la  ! angular quantum number of basis a
   integer(kind=i4_kind), intent(in)    :: lb  ! angular quantum number of basis b
   real(kind=r8_kind),    intent(in)    :: alpha(:) ! exponents alpha
   real(kind=r8_kind),    intent(in)    :: beta(:)  ! exponents beta
   real(kind=r8_kind),    intent(in)    :: gamma    ! exponents of pseudopotentials
   real(kind=r8_kind),    intent(in)    :: ca       ! distance |a-c|
   real(kind=r8_kind),    intent(in)    :: cb       ! distance |b-c|
   real(kind=r8_kind),    intent(in)    :: ab       ! distance |a-b|
   real(kind=r8_kind),    intent(in)    :: expmax   ! cutoff parameter
   real(kind=r8_kind),    intent(in)    :: factor
   real(kind=r8_kind),    intent(inout) :: f(size(beta),size(alpha),&
                                           1:la+1,1:lb+1) 
   !--- end of interface ---------------------------------------------------------
   integer(kind=i4_kind)                :: na,nb,j,k,&
                                           iflag,lama,lbmb
   real(kind=r8_kind)                   :: arg,qq,alpha_2
   real(kind=r8_kind)                   :: rka,rkb,argsum
   real(kind=r8_kind)                   :: alpha_1
   !--- end of local variables ---------------------------------------------------

   ASSERT(   la<=6) ! abess
   ASSERT(   lb<=6) ! bbess

   na = size(alpha)
   nb = size(beta)

   do k=1,nb
     do j=1,na
       if (((alpha(j)*beta(k))/(alpha(j)+beta(k)))* ab**2 .lt. &
           expmax) then
         rka = two*alpha(j)*ca
         rkb = two*beta(k)*cb
         alpha_1 = alpha(j)+beta(k)+gamma
         alpha_2 = alpha_1
         argsum = ((rka+rkb)**2)/(two*alpha_1)

         if(argsum.gt.hundred) then
           iflag = 1
           arg=(-alpha(j)*beta(k)*(ca-cb)**2 &
                -gamma*(alpha(j)*ca*ca+beta(k)*cb*cb))/alpha_1
           call ptprep(nj,lm,la,lb,rka,rkb,alpha_1)
         elseif(rka.le.rkb) then
           iflag = 2
           arg=(-alpha(j)*ca*ca*alpha_1-beta(k)*cb*cb*&
                (alpha(j)+gamma))/alpha_1
         else
           iflag = 3
           arg=(-alpha(j)*ca*ca*(beta(k)+gamma)-beta(k)*cb*cb*&
                alpha_1)/alpha_1
         end if
           
           do lama = 0,la
              do lbmb = 0,lb
                 if (iflag.ne.1) then
                    qq=qpasy(nj+lama+lbmb,lm+lama,lm+lbmb,alpha_2,&
                             rka,rkb,iflag)
                 else
                    qq=ptwt(nj,lm,lama,lbmb)
                 end if
                    f(k,j,1+lama,1+lbmb)=f(k,j,1+lama,1+lbmb)+factor*&
                    qq*exp(arg) ! exparg

              end do
           end do
       end if                   
      end do
   end do
 end subroutine bessel2
  !******************************************************************************

   function qcomp(n,l,rk_1,alpha_1,t_1) result(f)
     !  This function subroutine is used to evaluate
     !  the integrals of Bessel function combined with
     !  algebric and exponential functions. The integrals
     !  result in the hypergeometric series that can be
     !  separated into qalt, qasy and qpow , which are
     !  asymptotic, analytic forms and power series.
     integer(kind=i4_kind),  intent(in)  :: n    !  the index of variable r
     integer(kind=i4_kind),  intent(in)  :: l    !  the orderr of Bessel function
     real(kind=r8_kind)   ,  intent(in)  :: rk_1
     real(kind=r8_kind)   ,  intent(in)  :: alpha_1
     real(kind=r8_kind)   ,  intent(in)  :: t_1
     real(kind=r8_kind)                  :: f    !  results
     !xxxxx end of interface xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!    if (iand(n+l,1).eq.0) then
!       if (n.gt.l) then
!         f = qalt ( n,l)
!       end if
!    else 
!     if ( n.ge.9 ) then
!        if ( t_1 > 15_r8_kind ) then
!           f = qasy ( n,l )
!        else
!           f = qpow ( n,l )
!        endif
!     else
!        if ( t_1 > tmin(n+1) ) then
!           f = qasy ( n,l )
!        else
!           f = qpow ( n,l )
!        endif
!      end if
!    endif
      IF ( IAND(N+L,1).NE.0) GO TO 100
      IF(N.LE.L) GO TO 100
      f=QALT(N,L,rk_1,alpha_1,t_1)
      RETURN
100   CONTINUE
      IF(N.LT.9) GO TO 110
      IF(T_1.LT.15.0_r8_kind) GO TO 130
      GO TO 120
110   CONTINUE
      IF(T_1.LT.TMIN(N+1)) GO TO 130
120   CONTINUE
      f=QASY(N,L,rk_1,alpha_1,t_1)
      RETURN
130   CONTINUE
      f=QPOW(N,L,rk_1,alpha_1,t_1)
      RETURN
  end function qcomp

  function qpow(n,l,rk_1,alpha_1,t_1) result(f)
    !  The power series is used to evaluate integrals
    integer(kind=i4_kind),  intent(in)  :: n    !  the index of variable r
    integer(kind=i4_kind),  intent(in)  :: l    !  the orderr of Bessel function
    real(kind=r8_kind)   ,  intent(in)  :: rk_1
    real(kind=r8_kind)   ,  intent(in)  :: alpha_1
    real(kind=r8_kind)   ,  intent(in)  :: t_1
    real(kind=r8_kind)                  :: f    !  results

    !xxxxx end of interface xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    !xxxxx local variables  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    real(kind=r8_kind)                  :: term, sum, perfac,&
                                           xnum, xden, xj
    integer(kind=i4_kind)               :: isum
    real(kind=r8_kind), parameter       :: con1 = 1.0e-14_r8_kind
    real(kind=r8_kind), parameter       :: con2 = 0.4e-14_r8_kind
    !xxxxx end of local variables xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!!$    if ( l == 0 ) xkp = one
!!$    if ( l /= 0 ) xkp = rk_1**l

    perfac = exp(-t_1) / sqrt ((two*alpha_1)**(n+l+1))
    xnum = real(l+n-1,r8_kind)
    xden = real(l+l+1,r8_kind)
    term = dfac(l+n+1)/dfac(l+l+3)
    sum  = term
    xj   = zero
    isum = 1
    !xxxxx being program loop xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    do_block_loop:  do
                    xnum = xnum + two
                    xden = xden + two
                    xj   = xj   + one
                    term = term * t_1 * xnum / ( xj  *  xden ) 
                    sum  = sum + term
                    if (( term/sum).le.con1 ) exit do_block_loop
                     isum = isum + 1
                     if(isum.gt.100)then
                     write(1,*)" in qpow, sum is ",sum
                     stop
                     endif
                    cycle do_block_loop
    end do do_block_loop

    !xxxxx end of program loop xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    f = perfac * sum
    if ( mod((l+n),2) == 0) f = f * sqpi2
    f = f * ( one + xj * con2 )
    !write(1,fmt='("in qpow, l = ",i4,"  n = , ",i4," f = , ",e10.4)') l,n,f
  end function qpow
    
  function qalt(n,l,rk_1,alpha_1,t_1) result(f)
    !  The power series is used to evaluate integrals
    integer(kind=i4_kind),  intent(in)  :: n    !  the index of variable r
    integer(kind=i4_kind),  intent(in)  :: l    !  the orderr of Bessel function
    real(kind=r8_kind)   ,  intent(in)  :: rk_1
    real(kind=r8_kind)   ,  intent(in)  :: alpha_1
    real(kind=r8_kind)   ,  intent(in)  :: t_1
    real(kind=r8_kind)                  :: f    !  results

    !xxxxx end of interface xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    !xxxxx local variables  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    real(kind=r8_kind)                  :: term, sum, perfac, xc
    integer(kind=i4_kind)               :: xnum, num, xden,isum

    !xxxxx end of local variables xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!!$    if ( l == 0 ) then
!!$       xkp = one
!!$    else
!!$       xkp = rk_1 ** L
!!$    endif
    perfac = sqpi2 * dfac( n+l+1 )/ &
             (sqrt((two*alpha_1)**(n+l+1)) * dfac ( l + l + 3 )) 
    num    = l - n + 2
    xden   = l + l + 3
    term   = one
    sum    = term
    xc     = - one
    isum = 1
   !xxxxx end of local variables xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   do_block_loop:   do
                    if ( num == 0 ) exit do_block_loop
                    xnum  = num
                    term  = term * real( xnum ) * t_1 / &
                            (real(xden,r8_kind) * real(xc,r8_kind))
                    xc    = xc - one
                    sum   = sum + term
                    num   = num + 2
                    xden  = xden + 2
                     isum = isum + 1
                     if(isum.gt.100)then
                     write(1,*)" in qalt, sum is ",sum
                     write(1,*)" num          is ",num
                     stop
                     endif
                    cycle do_block_loop
    end do do_block_loop

    !xxxxx end of program loop xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    f = perfac * sum

    !write(1,fmt='("in qalt, l = ",i4,"  n = , ",i4," f = , ",e10.4)') l,n,f
    !write(1,*)' in qalt result =, ',f
  end function qalt

  function qasy(n,l,rk_1,alpha_1,t_1) result(f)
    !  The power series is used to evaluate integrals
    integer(kind=i4_kind),  intent(in)  :: n    !  the index of variable r
    integer(kind=i4_kind),  intent(in)  :: l    !  the orderr of Bessel function
    real(kind=r8_kind)   ,  intent(in)  :: rk_1
    real(kind=r8_kind)   ,  intent(in)  :: alpha_1
    real(kind=r8_kind)   ,  intent(in)  :: t_1
    real(kind=r8_kind)                  :: f    !  results

    !xxxxx end of interface xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    !xxxxx local variables  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    real(kind=r8_kind)                  :: xkp, sum, perfac, xc, fac1, &
                                           fac2, tn, to
    real(kind=r8_kind), parameter       :: eps  = 1.0e-12_r8_kind
    real(kind=r8_kind), parameter       :: ep2  = 0.40e-14_r8_kind
    integer(kind=i4_kind)               :: isum
    !xxxxx end of local variables xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!!$    print *,'qasy: n=',n,' l=',l,' rk=',rk_1
    xkp  = rk_1 ** ( n - 2 ) / rk_1**l
    perfac = xkp * sqpi2 / sqrt ((two*alpha_1)**(2*n-1))
    sum   =  one
    to    =  one
    fac1  =  real(l-n+2,r8_kind)
    fac2  =  real(1-l-n,r8_kind)
    xc    =  one
    isum = 1
    do_block_loop:   do
                     tn = to * fac1 * fac2 / ( four * xc * t_1 )
                     if ( tn == zero ) exit do_block_loop
                     sum  = sum + tn
                     if ( abs(tn/sum) < eps ) exit do_block_loop
                     fac1 = fac1 + two
                     fac2 = fac2 + two
                     xc   = xc   + one
                     to   = tn
                     isum = isum + 1
                     if(isum.gt.100)then
                     write(1,*)" in qasy, sum is ",sum
                     stop
                     endif
                     cycle do_block_loop
    end do do_block_loop

    !xxxxx end of program loop xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    f = perfac * sum

    f = f * ( one + xc * ep2 )
    !write(1,fmt='("in qasy, l = ",i4,"  n = , ",i4," f = , ",e10.4)') l,n,f
  end function qasy

  !xxxxx function integral_radial xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  function integral_radial(n_j,alpha) result(f)
    !
    ! This functional subroutine is used to evluate the integrals
    ! which are simple exponential integrals.
    !
    implicit none
    integer(kind=i4_kind), intent(in) :: n_j                ! index of power
    real(kind=r8_kind)   , intent(in) :: alpha(:)           ! exponential
    real(kind=r8_kind)                :: f(size(alpha))
    ! *** end of interface ***

    integer(kind=i4_kind) :: n, k

    ! Target series member is f(n_j + 1),
    ! Set "n" to 1 or 2 for odd and even (n_j + 1):
    n = 2 - mod( n_j + 1, 2)

    ! first two members in the series:
    if( n == 1 )then
      ! n_j + 1 is odd, set: f(1):
      f = sqrt ( pi / alpha ) * half
    else ! i.e. n == 2
      ! n_j + 1 is even, set f(2):
      f = one / ( two * alpha )
    endif

    !
    ! Now use recurrence
    !
    !    f(k) = ( k - 2 ) * f( k - 2 ) / 2a
    !
    ! from k = 3 or 4 up to k = n_j + 1
    !
    ! start from k=3 for odd  (n_j + 1)
    ! or    from k=4 for even (n_j + 1)
    do k = n + 2, n_j + 1, 2
       ! f(k) = ( k - 2 ) * f( k - 2 ) / 2a:
       f = ( k - 2 ) * f / ( two * alpha )
    end do
    ! FIXME: one may accumulate the product of factors (k-2)/2
    ! and together with the correspnding power "p" of the array alpha**p
    ! compute result in closed form. But I dont want to check
    ! accuracy implications now, also this funciton is hardly performance
    ! critical, I guess.
  end function integral_radial

end module bessel_module

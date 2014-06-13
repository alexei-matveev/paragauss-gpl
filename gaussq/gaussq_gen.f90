module gaussq_gen

  implicit none
  ! number precision:
!!$  integer, parameter :: DP  = selected_real_kind(15)
  integer, parameter :: DP  = selected_real_kind(23)
!!$  integer, parameter :: SP  = selected_real_kind(6)
  integer, parameter :: SP  = selected_real_kind(15)
  integer, parameter :: I4B = selected_int_kind(9)

contains

  subroutine f90mod(n)
    integer(I4B), intent(in) :: n
    ! *** end of interface ***

    integer  :: i

    do i=1,n
       call tabulate(i)
    enddo
  end subroutine f90mod

  subroutine tabulate(n)
    integer(I4B), intent(in) :: n
    ! *** end of interface ***

    real(DP) :: x(n), w(n)
    integer  :: i,loc

    call GAUHER(x,w)
100 FORMAT(X,'! GAUSS-HERMITE ABSCISSAS AND WEIGHTS')
101 FORMAT(X,'! N=',I4)
102 FORMAT(X,'DATA ((ZW(i,j),j=1,2),i=',I4,',',I4,'), LOC(',I2,') / &')
103 FORMAT(F30.22,'_RP,',1PE30.20,'_RP, &')
104 FORMAT(I8,' /')
    ! write(*,100)
    write(*,101) n
    loc = 1+n*(n-1)/2
    write(*,102) loc, loc+n-1, n

    do i=1,(n+1)/2
       if(i<=n/2)then
          ! dont print zero twice!
          write(*,103) -x(i),w(i)
          write(*,103)  x(i),w(i)
       else
          write(*,103)  x(i),w(i)
       endif
    enddo

!!$    do i=n,1,-1
!!$       write(*,103)  x(i),w(i)
!!$    enddo    

    write(*,104) loc
  end subroutine tabulate

  SUBROUTINE gauher(x,w)
!!$  USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
!!$    REAL(DP), PARAMETER :: EPS=3.0e-13_dp
    REAL(DP), PARAMETER :: EPS=3.0e-23_dp
    REAL(DP), PARAMETER :: PI=3.1415926535897932384626433832794_dp
    REAL(DP)            :: PIM4 !=PI**(0.25_dp) !=0.7511255444649425_dp
    !This routine returns arrays x and w of length N containing the abscissas and weights of
    !the N-point Gauss-Hermite quadrature formula. The abscissas are returned in descending
    !order. Note that internal computations are done in double precision.
    !Parameters: EPS is the relative precision, PIM4 = 1/(PI^1/4).
    
    INTEGER(I4B) :: its,j,m,n
    INTEGER(I4B), PARAMETER :: MAXIT=10
    REAL(SP) :: anu
    REAL(SP), PARAMETER :: C1=9.084064e-01_sp,C2=5.214976e-02_sp,&
         C3=2.579930e-03_sp,C4=3.986126e-03_sp
    REAL(SP), DIMENSION((size(x)+1)/2) :: rhs,r2,r3,theta
    REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    LOGICAL , DIMENSION((size(x)+1)/2) :: unfinished

!!$  n=assert_eq(size(x),size(w),’gauher’)
    n = size(x)

    m=(n+1)/2 ! The roots are symmetric about the origin, so we have to
    ! find only half of them
    anu=2.0_sp*n+1.0_sp
!!$  rhs=arth(3,4,m)*PI/anu ! arth(first,increment,n)
    ! returns an arithmetic progression as an array.
    do j=1,m
       rhs(j) = ( 3 + 4*(j-1) )*PI/anu
    enddo

    r3=rhs**(1.0_sp/3.0_sp)
    r2=r3**2
    theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
    z=sqrt(anu)*cos(theta) ! Initial approximations to the roots.

!!$    print *,'DIFF PI=',PI-4.0_DP*ATAN(1.0_DP)
    PIM4 = PI**(-0.25_dp)
!!$    print *,'PIM4=',PIM4

    unfinished=.true.
    do its=1,MAXIT         ! Newton’s method carried out simultaneously on the roots.
       where (unfinished)
          p1=PIM4
          p2=0.0
       end where
       do j=1,n           ! Loop up the recurrence relation to get the Hermite poly-
          ! nomials evaluated at z.
          where (unfinished) 
             p3=p2
             p2=p1
             p1=z*sqrt(2.0_dp/j)*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
          end where
       end do
       ! p1 now contains the desired Hermite polynomials. We next compute pp, the derivatives,
       ! by the relation (4.5.21) using p2, the polynomials of one lower order.
       where (unfinished)
          pp=sqrt(2.0_dp*n)*p2
          z1=z
          z=z1-p1/pp ! Newton’s formula.
          unfinished=(abs(z-z1) > EPS)
       end where
       if (.not. any(unfinished)) exit
    end do

    if (its == MAXIT+1) then
       print *,"too many iterations in gauher, >",MAXIT
       stop
    endif

    x(1:m)=z             ! Store the root
    x(n:n-m+1:-1)=-z     ! and its symmetric counterpart.
    w(1:m)=2.0_dp/pp**2  ! Compute the weight
    w(n:n-m+1:-1)=w(1:m) ! and its symmetric counterpart.
  END SUBROUTINE gauher

end module gaussq_gen

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
module math_module
  !---------------------------------------------------------------
  !
  !  Purpose: module contains subroutines for various
  !           mathematical tasks
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
#include <def.h>
  use type_module ! type specification parameters
  implicit none
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ public functions and subroutines ------------------
  public cross_product,abs_value,rotmat,invert_matrix,gen_invert_matrix, &
       dot_prod,binomial_coeff,schmidt,round,eigvec_sign_check,&
       print_matrix,equal_vector,ortho2,gauss_solve

  public :: invert_matrix_sd

  !================================================================
  ! End of public interface of module
  !================================================================


  real(kind=r8_kind),parameter,public :: one=1.0_r8_kind,zero=0.0_r8_kind, &
       two=2.0_r8_kind,three=3.0_r8_kind,&
       four=4.0_r8_kind,five=5.0_r8_kind,&
       half=0.5_r8_kind
  real(kind=r8_kind),parameter,public :: small=1.0e-10_r8_kind
  real(kind=r8_kind),parameter,public :: pi=3.14159265358979324_r8_kind, &
       convert1 = 180.0_r8_kind / pi

contains

  function cross_product(vec1,vec2)
    real(kind=r8_kind),intent(in)   :: vec1(3),vec2(3)
    real(kind=r8_kind)              :: cross_product(3)
    !** End of interface *****************************************

    cross_product(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_product(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_product(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  end function cross_product

  function dot_prod(vec1,vec2)
    real(kind=r8_kind),intent(in)   :: vec1(3),vec2(3)
    real(kind=r8_kind)              :: dot_prod
    !** End of interface *****************************************

    dot_prod = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
  end function dot_prod


  function abs_value(vec1)
    real(kind=r8_kind)             :: abs_value
    real(kind=r8_kind),intent(in)  :: vec1(:)
    !** End of interface *****************************************

    abs_value = sqrt(sum(vec1**2))
  end function abs_value

  function binomial_coeff(n,k)
    real(kind=r8_kind)     :: binomial_coeff
    integer(kind=i4_kind)  :: n,k
    integer(kind=i4_kind)  :: ik,k_fac,nom
    !** End of interface *****************************************

    k_fac=1
    do ik=1,k
       k_fac=k_fac*ik
    enddo
    nom = 1
    do ik=0,k-1
       nom = nom*(n-ik)
    enddo
    binomial_coeff = nom/k_fac
  end function binomial_coeff

  function rotmat(axis,angle)
    integer(kind=i4_kind)  :: axis ! 1=x,2=y,3=z
    real(kind=r8_kind)     :: angle
    real(kind=r8_kind)     :: rotmat(3,3)
    !** End of interface *****************************************

    real(kind=r8_kind)     :: c_a,s_a

    c_a = cos(angle)
    s_a = sin(angle)

    if (axis == 1 ) then
       rotmat(1,1) = one
       rotmat(1,2) = zero
       rotmat(1,3) = zero
       rotmat(2,1) = zero
       rotmat(2,2) = c_a
       rotmat(2,3) = -s_a
       rotmat(3,1) = zero
       rotmat(3,2) = s_a
       rotmat(3,3) = c_a
       elseif (axis == 2 ) then
       rotmat(1,1) = c_a
       rotmat(1,2) = zero
       rotmat(1,3) = -s_a
       rotmat(2,1) = zero
       rotmat(2,2) = one
       rotmat(2,3) = zero
       rotmat(3,1) = c_a
       rotmat(3,2) = zero
       rotmat(3,3) = s_a
       elseif ( axis == 3) then
       rotmat(1,1) = c_a
       rotmat(1,2) = -s_a
       rotmat(1,3) = zero
       rotmat(2,1) = s_a
       rotmat(2,2) = c_a
       rotmat(2,3) = zero
       rotmat(3,1) = zero
       rotmat(3,2) = zero
       rotmat(3,3) = one
    else
       call error_handler(' rotmat : stupid axis specified')
    endif

  end function rotmat

  function round(number,digit)
    real(kind=r8_kind)     :: number,round
    integer(kind=i4_kind)  :: digit
    !** End of interface *****************************************
    !-----------------------
    real(kind=r8_kind) :: hilf
    integer(kind=i4_kind)  :: hilf_int
    if (number<=zero) then
       hilf = number*10.0_r8_kind**(real(digit))-half
    else
       hilf = number*10.0_r8_kind**(real(digit))+half
    endif
    hilf_int = int(hilf,kind=i4_kind)
    round=real(hilf_int,kind=r8_kind)/10.0_r8_kind**(real(digit))
  end function round

  function equal_vector(vec1,vec2)
    logical            :: equal_vector
    real(kind=r8_kind) :: vec1(:),vec2(:)
    !** End of interface *****************************************
    ! ----------------------------------
    integer(kind=i4_kind)  :: dim1,dim2,summ,i

    equal_vector=.false.
    dim1=ubound(vec1,1)
    dim2=ubound(vec2,1)
    if (dim1/=dim2) return
    summ=0
    do i=1,dim1
       if ( abs(vec1(i)-vec2(i))<small ) summ=summ+1
    enddo
    if (summ==dim1) then
       equal_vector=.true.
    endif
    return
  end function equal_vector

  subroutine invert_matrix(dimen,matrix)
    ! Purpose : wrapper for the dge*-LAPACK Routines for
    !          the special case of quadratic matrices (i.e. dim1=dim2).
    !          Up to now only the double-precision general-matrix
    !          routines are implemented, although  for the Hessian
    !          which is supposed to be symmetric (and positive definite
    !          if a minimum is searched) other routines would
    !          be appropriate. This has to be chagned lateron.
    !-------------------------------------------------------
    use f77_lapack, only: dgetrf,dgetri
    integer(kind=i4_kind),intent(in)  :: dimen
    real(kind=r8_kind),intent(inout)  :: matrix(:,:)
    !** End of interface *****************************************
    ! ------- Declaration of local variables ---------------
    integer(kind=i4_kind)    :: alloc_stat,info
    integer(kind=i4_kind),allocatable :: ipiv(:)
    real(kind=r8_kind),allocatable    :: work(:)
    ! ------ Executable code -------------------------------
    allocate(ipiv(dimen),STAT=alloc_stat)
    if (alloc_stat/=0) &
         call error_handler(' invert_matrix : allocation (1) failed')

    allocate(work(dimen),STAT=alloc_stat)
    if (alloc_stat/=0) &
         call error_handler(' invert_matrix : allocation (2) failed')
    info=0
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(dimen,dimen,matrix,dimen,ipiv,info)
    if (info < 0 ) then
       call error_handler&
            ("invert_matrix: dgetrf exited with an illegal value as argument")
    elseif (info > 0) then
       write(*,*)" dgetrf:  element ",info," is exactly zero in the LU-factorization"
       write(*,*)"          input matrix is singular"
       call error_handler("invert_matrix: dgetrf exited")
    endif
    info=0

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call dgetri(dimen,matrix,dimen,ipiv,work,dimen,info)
    if (info < 0 ) then
       write(*,*)' invert_bmat : warning - dgetri exited with info < 0'
       write(*,*)" argument ",-info," had an illegal value"
       call error_handler(" invert_matrix: exit")
    elseif (info > 0 ) then
       write(*,*)" invert_matrix: the matrix is singular and cannot be inverted"
       call error_handler(" invert_matrix: exit")
    endif
    deallocate(work,ipiv,STAT=alloc_stat)
    if (alloc_stat/=0) &
         call error_handler('invert_matrix : deallocation (1) failed')

  end subroutine invert_matrix

  subroutine gen_invert_matrix(dimen,matrix)
    ! Purpose : generalized procedure of matrix inversion
    ! Mainly used to invert Cartesian Hessian
    !-------------------------------------------------------
    integer(i4_kind),intent(in)  :: dimen
    real(r8_kind),intent(inout)  :: matrix(:,:)
    !** End of interface *****************************************
    ! ------- Declaration of local variables ---------------
    integer(i4_kind)    :: status,info,i,j,k
    real(r8_kind), allocatable :: work(:),eval(:),mat(:,:)
    real(r8_kind), parameter :: small=1.0e-6_r8_kind
    ! ------ Executable code -------------------------------

    allocate(work(3*dimen-1),eval(dimen),mat(dimen,dimen),stat=status)
    ASSERT(status == 0)
    ASSERT(dimen==size(matrix,1))
    mat=matrix

    call dsyev('V','L',dimen,mat,dimen,eval,work,3*dimen-1,info)
    ASSERT(info == 0)

    deallocate(work,stat=status)
    ASSERT(status == 0)

    do i=1,dimen
       if(abs(eval(i)) <= small) then
          eval(i)=zero
       else
          eval(i)=one/eval(i)
       end if
    end do

    do i=1,dimen
       do j=1,dimen
          matrix(i,j)=zero
          do k=1,dimen
             matrix(i,j)=matrix(i,j)+mat(i,k)*mat(k,j)*eval(k)
          end do
       end do
    end do

    deallocate(eval,stat=status)
    ASSERT(status == 0)
    deallocate(mat,stat=status)
    ASSERT(status == 0)

  end subroutine gen_invert_matrix

  !*************************************************************
  subroutine invert_matrix_sd(dimen,matrix)
    !VVP:
    !Purpose: It's new procedure for  inversion of hessian (positive definite symmetric matrix).
    !Uses Cholesk(y|i) decomposition (H=matmul(transpose(L),L), where L is lower triangular matrix).
    !"Second order necessary condition". Please, see Fletcher, R. Practical methods
    !of optimization: A Wiley-Interscience publication,1987, P. 14

     use type_module
     implicit none
     integer(kind=i4_kind),intent(in)  :: dimen
     real(kind=r8_kind),intent(inout)  :: matrix(:,:)
     real(kind=r8_kind),allocatable    :: work(:),ainv(:)
     integer(kind=i4_kind)              :: alloc_stat,IERR,n
     n=dimen*(dimen+1)/2
     allocate(work(n),STAT=alloc_stat)
     if (alloc_stat/=0) &
     call error_handler(' invert_matrix_sd : allocation (2) failed')
     call COMPACT(matrix,dimen,work)
     allocate(ainv(n),STAT=alloc_stat)
     if (alloc_stat/=0) &
     call error_handler(' invert_matrix_sd : allocation (2) failed')
     call AIH1D(work,ainv,dimen,IERR)
     call COMPLETE(ainv,dimen,matrix)
     deallocate(work,STAT=alloc_stat)
     if (alloc_stat/=0) &
     call error_handler('invert_matrix_sd : deallocation (1) failed')
     deallocate(ainv,STAT=alloc_stat)
     if (alloc_stat/=0) &
     call error_handler('invert_matrix_sd : deallocation (1) failed')
     return
   contains

     SUBROUTINE AIH1D(A,AINV,N,IERR)
       use type_module
       implicit none
       integer(kind=i4_kind), intent(in)    :: N
       real(kind=r8_kind),intent(inout)     :: A(N*(N+1)/2)
       real(kind=r8_kind),intent(out)      :: AINV(N*(N+1)/2)
       integer(kind=i4_kind)               :: IERR,L,N1,LN,I,J,I1
       real(kind=r8_kind)                  :: ZERO,ONE
       DATA  ZERO,ONE/0.0_r8_kind,1.0_r8_kind/
       IERR=0
       L=1
       N1=N-1
       LN=N
       CALL AFH1D(A,N,IERR)
       IF(IERR.NE.0) GO TO 3
       DO 2 I=1,N
          DO 1 J=L,LN
             AINV(J)=ZERO
1         CONTINUE
          I1=L+I-1
          AINV(I1)=ONE
          CALL ASH1DT(A,AINV(L),AINV(L),N)
          L=L+I
          LN=L+N1
2      CONTINUE
       GO TO 4
3      call error_handler&
            ("FATAL ERROR AIH1D:The set matrix is not positive definite,therefore triangular decomposition of this matrix &
            & it is impossible")
4      RETURN
     END SUBROUTINE AIH1D

     SUBROUTINE AFH1D(A,N,IERR)
       use type_module
       integer(kind=i4_kind), intent(in)    :: N
       real(kind=r8_kind),intent(inout)     :: A(N*(N+1)/2)
       INTEGER(kind=i4_kind)               :: IERR,IP,IQ,IR,IP1,J,I,K
       real(kind=r8_kind)                  :: SIXTN,X,RN,DSQRT,ONE
       DATA  SIXTN/16.0_r8_kind/
       DATA  ONE/1.0_r8_kind/
       RN=ONE/(DFLOAT(N)*SIXTN)
       IP=1
       IERR=0
       DO 6 I=1,N
          IQ=IP
          IR=1
          DO 5 J=1,I
             X=A(IP)
             IF(J.EQ.1) GO TO 2
             DO 1 K=IQ,IP1
                X=X-A(K)*A(IR)
                IR=IR+1
1            CONTINUE
2            IF(I.NE.J) GO TO 3
             IF(A(IP)+X*RN.LE.A(IP)) GO TO 7
             A(IP)=DSQRT(X)
             GO TO 4
3            A(IP)=X/A(IR)
4            IP1=IP
             IP=IP+1
             IR=IR+1
5         CONTINUE
6      CONTINUE
       GO TO 8
7      IERR=65
       IF (IERR.NE.1) THEN
          call error_handler&
               ("FATAL ERROR AFH1D:The set matrix is not positive definite,therefore triangular decomposition of this matrix &
               &it is impossible")
       END IF
8      RETURN
     END SUBROUTINE AFH1D

     SUBROUTINE ASH1DT(A,B,X,N)
       integer(kind=i4_kind),intent(in)     :: N
       real(kind=r8_kind)                  :: A(N*(N+1)/2),B(N*(N+1)/2),X(N*(N+1)/2)
       integer(kind=i4_kind)               :: IP,IW,I,IM1,K,N1,II,IS,IQ,KK
       real(kind=r8_kind)                  :: T,ZERO
       DATA  ZERO/0.0_r8_kind/
       IP=1
       IW=0
       DO 4 I=1,N
          T=B(I)
          IM1=I-1
          IF(IW.EQ.0) GO TO 2
          IP=IP+IW-1
          DO 1 K=IW,IM1
             T=T-A(IP)*X(K)
             IP=IP+1
1         CONTINUE
          GO TO 3
2         IF(T.NE.ZERO) IW=I
          IP=IP+IM1
3         X(I)=T/A(IP)
          IP=IP+1
4      CONTINUE
       N1=N+1
       DO 7 I=1,N
          II=N1-I
          IP=IP-1
          IS=IP
          IQ=II+1
          T=X(II)
          IF(N.LT.IQ) GO TO 6
          KK=N
          DO 5 K=IQ,N
             T=T-A(IS)*X(KK)
             KK=KK-1
             IS=IS-KK
5         CONTINUE
6         X(II)=T/A(IS)
7      CONTINUE
       RETURN
     END SUBROUTINE ASH1DT

     SUBROUTINE COMPACT(A,N,B)
       use type_module
       implicit none
       integer(kind=i4_kind), intent(in) :: N
       real(kind=r8_kind),intent(in)     :: A(N,N)
       real(kind=r8_kind),intent(out)    :: B(N*(N+1)/2)
       integer(kind=i4_kind)             :: I,K,J
       K=1
       DO  I=1,N
          DO  J=1,I
             B(K)=A(J,I)
             K=K+1
          end do
       end do
       RETURN
     END SUBROUTINE COMPACT

     SUBROUTINE COMPLETE(A,N,B)
       use type_module
       implicit  none
       integer(kind=i4_kind), intent(in) :: N
       real(kind=r8_kind),intent(in)     :: A(N*(N+1)/2)
       real(kind=r8_kind),intent(out)    :: B(N,N)
       INTEGER(kind=i4_kind)            :: I1,J,JP1,K,I
       I1=(N*(N+1))/2
       J=N
1      JP1=J+1
       DO 2 K=1,J
          I=JP1-K
          B(I,J)=A(I1)
          I1=I1-1
2      CONTINUE
       J=J-1
       IF(J.GE.1) GO TO 1
       IF(N.LT.2) GO TO 5
       DO 4 I=2,N
          I1=I-1
          DO 3 J=1,I1
             B(I,J)=B(J,I)
3         CONTINUE
4      CONTINUE
5      RETURN
     END SUBROUTINE COMPLETE
  end subroutine invert_matrix_sd


  !*************************************************************
  subroutine eigvec_sign_check(dimen,eigvec,eigvec_prev)
    ! Purpose: compare the signs of the previous eigenvectors
    !          to the actual ones. If  ALL signs on an eigenvector
    !          have changed, reset them to those of the previous
    !          eigenvector.
    !          This is necessary for the backtransformation from
    !          delocalized internal coordinates to cartesians
    !          as well as for the comparison of internals and
    !          gradients between two geometry steps.
    !
    ! Routine called by: generate_bmat_reduced
    ! --------------------------------------------------------
    integer(kind=i4_kind),intent(in)    :: dimen
    real(kind=r8_kind),intent(inout)    :: eigvec(:,:)
    real(kind=r8_kind),intent(in)       :: eigvec_prev(:,:)
    !** End of interface *****************************************
    ! --- Declaration of local variables ---------------------
    real(kind=r8_kind)    :: dir_1(dimen),dir_2(dimen),c_phi
    integer(kind=i4_kind) :: i

    dir_1=zero
    dir_2=zero
    do i=1,dimen
       where (eigvec(:,i)>zero) dir_1 = one
       where ( eigvec(:,i)<zero) dir_1 = -one
       where ( abs(eigvec(:,i))<small) dir_1 = zero
       where (eigvec_prev(:,i)>zero) dir_2 = one
       where ( eigvec_prev(:,i)<zero) dir_2 = -one
       where ( abs(eigvec_prev(:,i))<small) dir_2 = zero

       dir_1 = dir_1/abs_value(dir_1)
       dir_2 = dir_2/abs_value(dir_2)
       c_phi = dot_product(dir_1,dir_2)
       if (abs(abs(c_phi)-one)<small .and. c_phi<zero ) then
          eigvec(:,i) = -eigvec(:,i)
       endif
    enddo
  end subroutine eigvec_sign_check

  !*************************************************************

  subroutine schmidt(n_cons,c_cons,vec_in,vec_out)
    ! Purpose: Perform a Schmidt-Orthogonalization using
    !          <n_cons> -already orthogonal- starting vectors
    !          <c_cons>.
    !          Care has to be taken for the case that one of the
    !          projected constraint vectors is equal to one
    !          of the vectors in 'vec_in'. In that case the
    !          vector has to be skipped in order to prevent its
    !          orthogonalization to zero.
    ! --------------------------------------------------
    integer(kind=i4_kind),intent(in) :: n_cons
    real(kind=r8_kind),intent(in)    :: c_cons(:,:)
    real(kind=r8_kind),intent(in)    :: vec_in(:,:)
    real(kind=r8_kind),intent(out)   :: vec_out(:,:)
    !** End of interface *****************************************
    ! ---------------------------------------------------------
    real(kind=r8_kind),parameter   :: small_ortho=1.0e-6_r8_kind
    integer(kind=i4_kind)          :: n_prim,n_int,alloc_stat,&
         i_vec,i_help,i_cons,counter,ii,n_reject
    real(kind=r8_kind),allocatable :: help_vec(:)
    real(kind=r8_kind)             :: test_prod
    logical,allocatable            :: list(:)

    n_prim = ubound(vec_in,1)
    n_int = ubound(vec_in,2)
    if (ubound(c_cons,1)/=n_prim) call error_handler &
         (" Schmidt: error in dimensions")

    vec_out=zero
    allocate(help_vec(n_prim),list(n_int),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler &
         ("Schmidt: allocation (1) failed")
    ! first generate the projected constraint vectors.
    do i_cons=1,n_cons
       help_vec=zero
       do i_help=1,n_int
          help_vec = help_vec + &
               dot_product(c_cons(:,i_cons),vec_in(:,i_help))*&
               vec_in(:,i_help)
       enddo
       vec_out(:,i_cons) = help_vec/abs_value(help_vec)
    enddo


    write(*,*)" The projected constraint vectors are:"
    do i_vec=1,n_prim
       write(*,'(6(3x,f9. 5))')(vec_out(i_vec,i_cons),i_cons=1,n_cons)
    enddo
    write(*,*)" "
    do i_vec=1,n_cons
       do ii=1,n_cons
          if (i_vec==ii) cycle
          test_prod = dot_product(vec_out(:,i_vec),vec_out(:,ii))
          if (abs(test_prod)>=small_ortho) then
             write(*,*)" schmidt: after projection onto the active subspace"
             write(*,*)"          the constraint vectors are no longer orthogonal."
             write(*,*)"          This means that your constraints imply further "
             write(*,*)"          constraints, that are not specified.            "
             write(*,*)"          Constraint vectors ",i_vec," and ",ii," are no longer orthogonal"
             !call error_handler(" Please revise your constraints")
          endif
       enddo
    enddo
    ! set a logical mask to true for those vectors of vec_in which are
    ! NOT equal to any of the (projected) constraint vectors.
    list = .true.
    do i_cons =1,n_cons
       do i_vec=1,n_int
          help_vec = vec_out(:,i_cons) - vec_in(:,i_vec)
          if (abs_value(help_vec)<small) list(i_vec) = .false.
       enddo
    enddo
    if (any(.not.list)) then
       write(*,*)" schmidt: one of the original vectors in UMAT equals "
       write(*,*)"          one of the constraint vectors. Skip this   "
       write(*,*)"          one in the Schmidt-Orthogonalization"
       write(*,*)" "
    endif

    counter=1
!    counter=n_cons
    n_reject=0
    vectors: do i_vec = 1,n_int
       if (list(i_vec)) then
          counter=counter+1
       else
          cycle vectors
       endif
       if (counter>n_int) exit vectors
       help_vec = zero
       do i_help = 1,counter-1
          help_vec = help_vec + &
               dot_product(vec_in(:,i_vec),vec_out(:,i_help))*&
               vec_out(:,i_help)
       enddo
       vec_out(:,counter) = (vec_in(:,i_vec) - help_vec)
       if ( abs_value(vec_out(:,counter))<=small ) then
          write(*,*)" schmidt: a vector dropped out of the "
          write(*,*)"          orthogonalization scheme"
          write(*,*)" "
          counter=counter-1
          cycle vectors
       endif
       vec_out(:,counter) = vec_out(:,counter)/abs_value(vec_out(:,counter))
       do ii=1,counter-1
          ! test if the latest vector is orthogonal to the
          ! previous ones
          test_prod = dot_product(vec_out(:,counter),vec_out(:,ii))
          if (test_prod>=small_ortho) then
             write(*,*)" schmidt: a vector was found not to be "
             write(*,*)"          orthogonal to the previous ones."
             write(*,*)"          Reject vector ",counter
             counter=counter-1
             n_reject=n_reject+1
             cycle vectors
          endif
       enddo
    enddo vectors
    print*,' Counter is now                   :',counter
    print*,' Total number of vectors required :',n_int
    print*,' Number of constraint vectors     :',n_cons
    print*,' Number of rejected vectors       :',n_reject

    deallocate(help_vec,list,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler&
         (" Schmidt: deallocation (1) failed")
  end subroutine schmidt

  !*************************************************************

  !*************************************************************
  subroutine ortho2(c_cons,vec_in,vec_out)
    use f77_lapack, only: dgels
    implicit none
    real(kind=r8_kind)   :: c_cons(:,:),vec_in(:,:),vec_out(:,:)
    ! ----------------------------------------------------------
    integer(kind=i4_kind)               :: n_constraint,n_internal,n_prim
    real(kind=r8_kind),allocatable      :: matrix(:,:),mat(:,:),&
                                           coeff(:,:),rhs(:),&
                                           work(:),help_vec(:),&
                                           eigval(:),eigvec(:,:),help_mat(:,:)
    integer(kind=i4_kind)               :: start,dimi,alloc_stat,&
                                           counter_j,counter_k,&
                                           i,j,k,counter_test,kk,&
                                           n_int,n_cons,counter_equal,&
                                           counter_cons
    integer(kind=i4_kind),allocatable   :: ind(:),ind_cons(:),ind_sort(:)
    integer                             :: m,n,lda,ldb,info,lwork,&
                                           nrhs
    real(kind=r8_kind)                  :: abs_val,minni
    real(kind=r8_kind),parameter        :: mat_test=1.0e5_r8_kind
!!$    external dgels

    n_constraint = ubound(c_cons,2) ! Anzahl constraints
    n_internal = ubound(vec_in,2)  ! Anzahl Vektoren im active set
    n_prim = ubound(vec_in,1) ! Anzahl primitiver Koordinate

    allocate(ind_sort(n_constraint),ind_cons(n_constraint),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler("ortho2: allocation (1) failed")
    ind_cons=0_i4_kind
    ind_sort=0_i4_kind

    ! Stelle zuerst fest, ob es im active set Vektoren gibt, die mit einem
    ! der constraint vektoren identisch sind.
    ! 'ind_sort' : enthaelt Indizes der Eigenvektoren, fuer die es identische
    !              unter den Constraint Vektoren gibt.
    ! 'ind_cons' : enthaelt Indizes der Constraint Vektoren, zu denen es
    !              KEINEN identischen unter den Eigenvektoren gibt.
    counter_equal=0
    counter_cons=0
    do k=1,n_constraint
       do j=1,n_internal
          if (equal_vector(c_cons(:,k),vec_in(:,j))) then
             counter_equal = counter_equal + 1
             ind_sort(counter_equal)=j
          else
             counter_cons=counter_cons+1
             ind_cons(counter_cons)=k
          endif
       enddo
    enddo
    ! Dies ist die Anzahl von Basisvektoren, die zur Konstruktion
    ! des neuen active sets benoetigt werden. Gehoert einer der Vektoren des alten
    ! active set gleichzeitig zu den Constraints, so kann er vom Raum, aus dem der Rest
    ! des neuen active sets konstruiert werden soll, ausgeschlossen werden.
    n_int = n_internal-counter_equal
    ! Gleichzeitig muss er auch von den restlichen constraint vektoren ausgeschlossen
    ! werden.
    n_cons=n_constraint-counter_equal
    if (n_cons/=n_constraint) then
       write(*,*)" ortho2: Among the Constraint Vectors ",counter_equal," were found "
       write(*,*)"         to be identical to one of the vectors of the active set.   "
       write(*,*)" ortho2: "
       write(*,*)" This is a list of indices of those eigenvectors of the (original) "
       write(*,*)" active set for which an identical constraint vector has been found:"
       write(*,*)" ----  ",ind_sort(1:counter_equal)," ---- "
       write(*,*)" "
    endif

    allocate(ind(n_int),coeff(n_int,n_int-n_cons),help_vec(n_prim),STAT=alloc_stat)
    if(alloc_stat/=0)call error_handler("ortho2: allocation (2) failed")
    coeff=zero
    ind=0_i4_kind

    if (n_cons==0) then ! das heisst es gab nur einen constraint vektor
       !                 und der war auch noch identisch mit einem der eigenvektoren
       vec_out=vec_in
       return
    endif

    i=0 ! Zaehlindex der Vektoren v
    vector_v: do
       i=i+1
       if (i > (n_int-n_cons) )exit vector_v  ! Wir brauchen n_int-n_cons
       !                                     ! Vektoren v
       ! Wie gross wird das zu loesende GLS sein?
       dimi = n_cons+i-1
       ! Die Koeffiziente von 1 ... start sollen frei gewaehlt werden,
       ! die von start+1 ... n_int sollen bestimmt werden:
       start = n_int-dimi

       allocate(matrix(dimi,n_int),mat(dimi,dimi),rhs(dimi),STAT=alloc_stat)
       if(alloc_stat/=0)call error_handler(" ortho2: allocation (3) failed")
       matrix=zero

       ! Hier wird 'ind_sort' lediglich dazu benutzt, um aus
       ! dem vollen active set die Eigenvektoren herauszusuchen,
       ! die mit keinem der Constraint vektoren gleich sind.
       k=n_internal+1
       counter_j=0
       full_active: do j=1,n_internal
          k=k-1
          do kk=1,counter_equal
             if (k==ind_sort(kk)) then
                cycle full_active
             endif
          enddo
          counter_j=counter_j+1
          ind(counter_j)=k
       enddo full_active

       ! Um groesstmoegliche Aehnlichkeit mit den input Vektoren
       ! zu erreichen, setze nur einen Koeffizienten auf eins und
       ! alle anderen frei waehlbaren auf null:
       counter_test=0
       matrix_test: do
          counter_test=counter_test+1
          if (start/=0 .and. counter_test>start+1) then
             write(*,*)" ortho2: No satisfactory solution to the LES could be found"
             write(*,*)"         Try to re-formulate your constraints OR change the"
             write(*,*)"         variable MAT_TEST in routine ORTHO2 to a higher   "
             write(*,*)"         value"
             call error_handler(" ")
          endif
          coeff(:,i) = zero
          ! This is how it worked allus fine
          if (start-counter_test+1>0) then
             coeff(ind(start-counter_test+1),i) = one
          endif
          if (start==0) coeff(ind(n_int),i)=one
          ! Zunaechst (wird spaeter weggelassen) wird die gesamt Matrix
          ! allociert und besetzt (matrix). Daraus wird die Matrix 'mat'
          ! isoliert und die rechte Seite ausgerechnet:'rhs'
          matrix=zero
          ! erst die Zeilen ,die der Bedingung v*c_cons=0 entsprechen
          ! Was kann hier Schlimmes passieren? Es koennte constraint vektoren
          ! geben, die vollstaendig ausserhalb des active set liegen!
          ! Dann waere eine Zeile von 'matrix' = 0
          do j=1,n_int
             do k=1,n_cons
                matrix(k,j) = dot_product(vec_in(:,ind(j)),c_cons(:,ind_cons(k)))
             enddo
          enddo
          if (i==1) then
             do k=1,n_cons
                if (abs(maxval(matrix(k,:)))<small .and. &
                     abs(minval(matrix(k,:)))<small ) then
                   call error_handler(" ortho_constraint: One of the constraint vectors&
                        & lies entirely out of the active set")
                endif
             enddo
          endif
          ! jetzt die Zeilen ,die der Bedingung entsprechen, das der aktuelle
          ! Vektor auf allen zuvor berechneten orthogonal ist
          do j=1,n_int
             counter_k=0
             do k=n_cons+1,dimi
                counter_k=counter_k+1
                matrix(k,j) = coeff(ind(j),counter_k)
             enddo
          enddo

          ! Jetzt berechne die rechte Seite:
          rhs=zero
          do k=1,dimi
             do j=1,start
                rhs(k) = rhs(k) + &
                     matrix(k,j)*coeff(ind(j),i)
             enddo
             rhs(k)=-rhs(k)
          enddo
          ! und jetzt setze 'mat' auf:
          mat=zero
          do k=1,dimi
             counter_j=0
             do j=start+1,n_int
                counter_j=counter_j+1
                if (counter_j>dimi) call error_handler&
                     ("ortho2: sth. fishy (1). Seek technical assistance")
                mat(k,counter_j) = matrix(k,j)
             enddo
          enddo

          m=int(dimi)
          n=m
          lda=m
          ldb=m
          info=0
          lwork=2*m
          allocate(work(lwork),STAT=alloc_stat)
          if(alloc_stat/=0) call error_handler(" ortho2: allocation (4) failed")
          work=zero
          nrhs=1
          call dgels('N',m,n,nrhs,mat,lda,rhs,ldb,work,lwork,info)
          if (info<0) then
             write(*,*)" ortho_constraint: on entry to the routine DGELS"
             write(*,*)"                   the ",-info," th argument    "
             write(*,*)"                   had an illegal value         "
             call error_handler(" ")
          endif
          deallocate(work,STAT=alloc_stat)
          if (alloc_stat/=0) call error_handler(" ortho2: deallocation (1) failed")

          !Wenn die Loesung zu grosse Komponenten hat, d.h. , wenn der
          !anfaenglich gewaehlte Vektor (durch Wahl von Coeff(1:start,i)
          !bestimmt) zu sehr *verdreht* werden musste, waehle
          !die STartkoeffizienten anders und mache den Mist noch mal
          if (abs(maxval(rhs))<mat_test) then
             exit matrix_test
          endif
          write(*,*)"ortho2: The solution to the LES just found is dropped due to "
          write(*,*)"        too high values in the solution vector. Try again with"
          write(*,*)"        a different initial choice of the arbitrary coefficients"

       enddo matrix_test

       ! Die Loesung des GLS steckt jetzt in rhs - umschaufeln
       ! nach coeff(start+1:n_int,i)
       counter_j=start
       do j=1,dimi
          counter_j=counter_j+1
          coeff(ind(counter_j),i)=rhs(j)
       enddo

       deallocate(matrix,mat,rhs,STAT=alloc_stat)
       if (alloc_stat/=0) call error_handler(" ortho2: deallocation (2) failed")

       ! normiere den so erhaltenen Vektor
       help_vec = zero
       do j=1,n_int
          help_vec=help_vec+coeff(ind(j),i)*vec_in(:,ind(j))
       enddo
       abs_val = abs_value(help_vec)
       if (abs_val<=small) then
          write(*,*)"ortho2: Sth. fishy (2). The computed vector has length zero"
          call error_handler(" ")
       endif
       coeff(:,i)=coeff(:,i)/abs_val

    enddo vector_v
    vec_out = zero
    do k=1,n_cons
       vec_out(:,k) = c_cons(:,k)
    enddo
    counter_k=0
    do k=n_cons+1,n_int
       counter_k=counter_k+1
       do j=1,n_int
          vec_out(:,k) = vec_out(:,k) + coeff(j,counter_k)*vec_in(:,j)
       enddo
       abs_val = abs_value(vec_out(:,k))
       if (abs(abs_val-one)>small) call error_handler &
            (" ortho2: Sth. fishy (3). The computed vector does not have length one")
       do kk=1,k-1
          abs_val=dot_product(vec_out(:,k),vec_out(:,kk))
          if (abs(abs_val)>1.0e-8) then
             write(*,*)" Ortho2: Two vectors are found to be non-orthogonal"
             write(*,*)"         The calculated dotproduct is :",abs_val
             call error_handler(" ")
          endif
       enddo
    enddo

    ! In order to ensure the the projection of the new active set onto
    ! the original active set does not introduce new linear dependencies
    ! calculate the eigenvalues of the new set
    allocate(eigval(n_internal),eigvec(n_internal,n_internal),&
         help_mat(n_internal,n_internal),STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler(" ortho2: allocation (5) failed")
    eigval=zero
    eigvec=zero
    help_mat=zero
    help_mat=matmul(transpose(vec_out),vec_in)
    call eigensolver(help_mat,n_internal,eigval,eigvec)
    minni=minval(abs(eigval))
    write(*,*)"ortho2: the smallest eigenvalue the projection of the new active"
    write(*,*)"        set onto the old active set ",minni
    write(*,*)"        If this values is too close to zero, this means that your"
    write(*,*)"        constraints contain linear dependencies "
    if (minni < 1.0e-8_r8_kind ) then
       call error_handler("ortho2: Linear dependencies detected in the constraints")
    endif

    deallocate (eigval,eigvec,help_mat,STAT=alloc_stat)
    if (alloc_stat/=0) call error_handler(" ortho2: deallocation (3) failed")


  end subroutine ortho2
    !** End of interface *****************************************
  !*************************************************************

  subroutine print_matrix(matrix,n,m,n_col)
    ! Purpose: print out a matrix if dimension n x m in a pretty
    !          format with n_col columns.
    !          n : numbers rows
    !          m : number of columns
    ! ----------------------------------------------------------
    real(kind=r8_kind),intent(in)      :: matrix(:,:)
    integer(kind=i4_kind),intent(in)   :: n,m,n_col
    !** End of interface *****************************************
    ! --- declaration of local variables -----------------------
    integer(kind=i4_kind) :: n_blocks,n_rest,counter,i_start,i_end,&
         j,i,k

    n_blocks=ceiling(real(m,kind=r8_kind)/real(n_col,kind=r8_kind))
    n_rest = mod(m,n_col)
    counter=1_i4_kind
    do j=1,n_blocks
       write(*,1202)(i,i=counter,counter+n_col-1)
       counter=counter+n_col
       i_start=(j-1)*n_col+1
!       i_end=min(i_start+n_col-1,(n_blocks-1)*n_col+n_rest)
       i_end=min(i_start+n_col-1,m-(j-1)*n_col)
       do i=1,n
          write(*,1201)i,(matrix(i,k),k=i_start,i_end)
       enddo
       write(*,*)
    enddo

1201 format(i4,20(2x,f10.7))
1202 format(4x,20(8x,i3))

  end subroutine print_matrix

  subroutine gauss_solve(a,b,x)
# include "def.h"
    implicit none
    real(kind=r8_kind), intent (in)    :: a(:,:)
    real(kind=r8_kind), intent (in)    :: b(:)
    real(kind=r8_kind), intent (out)   :: x(:)

    integer(kind=i4_kind)              :: mloc(1)
    integer(kind=i4_kind)              :: k,n,alloc_stat
    real(kind=r8_kind),allocatable     :: d(:,:),drow(:)
    n=size(a,dim=1)
    allocate (d(n,n+1),STAT=alloc_stat)
    ASSERT(alloc_stat==0)
    allocate(drow(n+1),STAT=alloc_stat)
    ASSERT(alloc_stat==0)

    d(:,1:n)=a
    d(:,n+1)=b

!triangularization phase
    do k=1,n
!pivoting
       mloc=maxloc(abs(d(k:n,k)))
       drow(k:n+1)=d(k,k:n+1)
       d(k,k:n+1)=d(k-1+mloc(1),k:n+1)
       d(k-1+mloc(1),k:n+1)=drow(k:n+1)

       if(abs(d(k,k))<epsilon(d(1,1))) then
         return
       endif
       d(k,k:n+1)=d(k,k:n+1)/d(k,k)
       d(k+1:n,k+1:n+1)=d(k+1:n,k+1:n+1)-spread(d(k,k+1:n+1),1,n-k)&
       *spread(d(k+1:n,k),2,n-k+1)
    enddo

!back substitution phase
    do k=n,1,-1
       x(k)=d(k,n+1)-sum(d(k,k+1:n)*x(k+1:n))
    enddo
    deallocate (d,drow)
  end subroutine gauss_solve

!!$  !*************************************************************
!!$
!!$  !subroutine eigensolver(matrix,dimen,eigval,eigvec,destroy)
!!$    ! Purpose: wrapper for 'evvrsp'-Routine from EISPACK.
!!$    !          This routine is used for compatibilty with
!!$    !          A.Voityuk. Do not change it without having a good
!!$    !          reason for doing so.
!!$    !
!!$    ! ----------------------------------------------------------
!!$    !integer(kind=i4_kind),intent(in)   :: dimen
!!$    !real(kind=r8_kind),intent(inout)   :: matrix(dimen,dimen)
!!$    !real(kind=r8_kind),intent(inout)   :: eigval(dimen)
!!$    !real(kind=r8_kind),intent(inout)   :: eigvec(dimen,dimen)
!!$    !logical,optional                   :: destroy
!!$    !** End of interface *****************************************
!!$    ! -----------------------------------------------------------
!!$    !real(kind=r8_kind)             :: b(dimen,9)
!!$    !real(kind=r8_kind),allocatable :: help_mat(:,:)
!!$    integer(kind=i4_kind)          :: iwork(dimen)
!!$    integer(kind=i4_kind)          :: ierr,alloc_stat
!!$    logical                        :: local
!!$
!!$    if (present(destroy)) then
!!$       if (destroy) then
!!$          local = .true.
!!$       else
!!$          local= .false.
!!$       endif
!!$    else
!!$       local = .false.
!!$    endif
!!$
!!$    if (.not.local) then
!!$       allocate(help_mat(dimen,dimen),STAT=alloc_stat)
!!$       if (alloc_stat/=0) then
!!$          write(*,*)"eigensolver: allocation (1) failed"
!!$          stop 1
!!$       endif
!!$       help_mat = matrix
!!$    endif
!!$
!!$    eigval = zero
!!$    eigvec=zero
!!$    iwork = 0_i4_kind
!!$    b = zero
!!$    call evvrsp(dimen,dimen,dimen,matrix,b,iwork,eigval,eigvec,ierr)
!!$
!!$    if (.not.local) then
!!$       matrix=help_mat
!!$       deallocate(help_mat,STAT=alloc_stat)
!!$       if (alloc_stat/=0) then
!!$          write(*,*)" eigensolver: deallocation (1) failed"
!!$          stop 1
!!$       endif
!!$    endif
!!$    if(ierr/=0) then
!!$       if (ierr<0) then
!!$          write(*,*)"eigensolver: Iteration for eigevector ",ierr," failed"
!!$          stop
!!$       else
!!$          write(*,*)"eigensolver: Iteration for eigenvalue ",ierr," failed"
!!$       endif
!!$    endif
!!$
!!$  end subroutine eigensolver

  !*************************************************************

  !--------------- End of module ----------------------------------
end module math_module

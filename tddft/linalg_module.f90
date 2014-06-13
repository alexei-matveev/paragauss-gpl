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
MODULE  linalg_module
!---------------------------------------------------------------
!
!  Purpose: 
!  Contains wrap around subroutines for some LAPACK F77/F66 
!  subroutines:
!  (1) invert_sym_matrix: compute inverse of a real symmetric 
!      matrix via two LAPACK subroutines:
!      - DSYTRF - compute the factorization of a real symmetric
!        matrix A using the Bunch-Kaufman diagonal pivoting method
!      - DSYTRI - compute the inverse of a real symmetric indefinite
!        matrix A using the factorization A = U*D*U**T or 
!        A = L*D*L**T computed by DSYTRF
!  (2) DGEMM: BLAS Level 3 subroutine for matrix-matrix operation
!        C := alpha*op(A)*op(B) + beta*C
!      where alpha, beta = scalars
!            A,B,C are matrices, with op(A) a m by k matrix,
!            op(B) a k by n matrix and C a m by n matrix.
!  (3) DGEMM: BLAS Level 3 subroutine for matrix-vector operation
!        y := alpha*A*x + beta*y or
!        y := alpha*transpose(A)*x + beta*y or
!        y := alpha*conjugate(transpose(A))*x + beta*y
!      where alpha, beta = scalars, x,y = vectors of size n and
!            A is an m by n matrix
!
!  Module called by: nearly every module in which matrix operations
!                    are performed
!
!  References: LAPACK 2.0 reference guide, Blas 3.0 ref. guide
! 
!
!  Author: HH
!
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification 
! Author: HH
! Date:   3/98
! Description: Added wrapper for DGEMM/DGEMV BLAS Level 3 subroutines
!
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
# include "def.h"
  USE type_module ! type specification parameters
  IMPLICIT NONE
  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================


!------------ public functions and subroutines ------------------
  PUBLIC invert_sym_matrix, matmatmul, matvecmul

!================================================================
! End of public interface of module
!================================================================

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
CONTAINS


  !*************************************************************
  SUBROUTINE invert_sym_matrix(A_matrix)
    !  Purpose: 
    !  Wrap around for the LAPACK routine used to invert the 
    !  square real symmetric matrix A_matrix.
    !  Inversion of the matrix is done in two steps by the LAPACK
    !  subroutines 
    ! (1) DSYTRF - compute the factorization of a real symmetric
    !     matrix A using the Bunch-Kaufman diagonal pivoting method
    ! (2) DSYTRI - compute the inverse of a real symmetric indefinite
    !     matrix A using the factorization A = U*D*U**T or 
    !     A = L*D*L**T computed by DSYTRF
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    REAL(KIND=r8_kind),INTENT(inout) :: A_matrix(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    ! -- Variables used as arguments for LAPACK subroutine
    CHARACTER                         :: uplo
    INTEGER(KIND=i4_kind)             :: info,ldA,lwork,N
    INTEGER(KIND=i4_kind),ALLOCATABLE :: ipiv(:)
    REAL(KIND=r8_kind),ALLOCATABLE    :: work(:)
  
    ! -- Variables used in this subroutine only 
    INTEGER(KIND=i4_kind)             :: alloc_stat,i,j
    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler,dsytri,dsytrf
    !------------ Executable code -----------------------------------
    
    !-- First initialize argument variables for the LAPACK subroutines
  
    uplo  = 'L' ! lower triangular part of matrix is stored
    N   = SIZE(A_matrix,1)  ! number of rowws of square matrix

    if (N.eq.0) return

    ldA = N                 ! assume:leading dimensions of A = number of rows
    lwork=64*N              ! is recommended by NAG for these LAPACK subr

    ! allocate memory for the various work arrays
    ALLOCATE(work(lwork),ipiv(N),STAT=alloc_stat)
    IF (alloc_stat/=0) &
         CALL error_handler("invert_sym_matrix: allocation of work &
         &arrays for LAPACK failed")

    !-- LAPACK - compute the factorization of the real symmetric matrix 
    !   A  using the Bunch-Kaufman diagonal pivoting method
    CALL dsytrf( uplo, N, A_matrix, ldA, ipiv, work, lwork, info )
    IF (info<0) THEN
       CALL error_handler &
            ("invert_sym_matrix: dsytrf detected illegal argument")
    ELSEIF (info>0) THEN
       CALL error_handler &
            ("invert_sym_matrix: unknown error situation after dsytrf")
    ENDIF

    ! LAPACK - compute the inverse of the real symmetric indefinite
    ! matrix A_matrix.
    ! On exit the lower triangular part of A_matrix contains the
    ! inverse of A_matrix.
    CALL dsytri( uplo, N, A_matrix, ldA, ipiv, work, info )
    IF (info<0) THEN
       CALL error_handler &
            ("invert_sym_matrix: dsytri detected illegal argument")
    ELSEIF (info>0) THEN
       CALL error_handler &
            ("invert_sym_matrix: unknown error situation after dsytri")
    ENDIF
    ! release the memory
    DEALLOCATE(ipiv,work, stat=alloc_stat)
    IF (alloc_stat/=0) &
         CALL error_handler("invert_sym_matrix: deallocation of work &
         &arrays for LAPACK failed")

    ! Now copy lower triangular part to upper triangular part to
    ! get the full inverse matrix
    DO i=2, N
       DO j=1, i-1
          A_matrix(j,i)=A_matrix(i,j)
       END DO
    END DO

  END SUBROUTINE invert_sym_matrix
  !*************************************************************


  !*************************************************************
  SUBROUTINE matmatmul(A, B, C, trans_A, trans_B, alpha, beta)
    !  Purpose: 
    !  Wrap around for the BLAS Level 3 subroutine DGEMM:
    !  Calculate the matrix product
    !
    !   C := alpha*op(A)*op(B) + beta*C
    !
    !  where alpha, beta are scalars, and A, B and C are matrices, 
    !  with op(A) an m by k matrix, op(B) a k by n matrix and C an
    !  m by n matrix.
    !  Parameter trans_A (ditto for trans_B):
    !     `N' or `n'   Use m-by-k matrix A
    !     `T' or `t'   Use A-transpose, where A is a k-by-m matrix
    !     `C' or `c'   Use A-conjugate-transpose,
    !
    !  OPTIONAL PARAMETERS:
    !  trans_A,trans_B   (default: trans_A = "N" = trans_B)
    !  alpha,beta        (default: alpha = 1.0_r8_kind = beta)
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    REAL(KIND=r8_kind),INTENT(in   ) :: A(:,:),B(:,:)
    REAL(KIND=r8_kind),INTENT(inout) :: C(:,:)
    CHARACTER(LEN=1),  INTENT(in   ),OPTIONAL :: trans_A, trans_B
    REAL(KIND=r8_kind),INTENT(in   ),OPTIONAL :: alpha,beta
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    ! -- Variables used as arguments for BLAS subroutine
    CHARACTER(LEN=1)                  :: transa, transb
    REAL(KIND=r8_kind)                :: alpha_in, beta_in
    INTEGER(KIND=i4_kind)             :: m, n, k, lda, ldb, ldc
  
    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler,dgemm
    !------------ Executable code -----------------------------------
    
    !-- First initialize argument variables for the BLAS subroutine

    ! check optinal arguments
    IF (PRESENT(trans_A)) THEN
       transa = trans_A
    ELSE
       transa = "N"      ! default is op() = id()
    END IF
    IF (PRESENT(trans_B)) THEN
       transb = trans_B
    ELSE
       transb = "N"      ! default is op() = id()
    END IF
    IF(PRESENT(alpha)) THEN
       alpha_in = alpha
    ELSE
       alpha_in = 1.0_r8_kind
    END IF
    IF(PRESENT(beta)) THEN
       beta_in = beta
    ELSE
       ! ** old and wrong default; I meant beta=0 but did set beta=1
       ! beta_in = 1.0_r8_kind
       ! new and correct default: C is **NOT** added to the product
       beta_in=0.0_r8_kind
    END IF
  
    ! set dimensions of assumed shape dummy arrays
    lda = SIZE(A,1)    ! number of rows of matrix A
    if (lda.eq.0) return
    ldb = SIZE(B,1)    ! number of rows of matrix B
    if (ldb.eq.0) return
    ldc = SIZE(C,1)    ! number of rows of matrix C
    if (ldc.eq.0) return

    ! Set the dimensions m,k and n which are defined
    ! as follows:
    ! op(A) -> m times k
    ! op(B) -> k times n
    !    C  -> k times k
    IF(transa=='N' .OR. transa=='n') THEN 
       m = lda
       k = SIZE(A,2)     ! number of columns of matrix A
    ELSE IF(transa=='T' .OR. transa=='t' .OR. &
          & transa=='C' .OR. transa=='c') THEN
       m = SIZE(A,2)     ! number of rows of transpose(A)
       k = lda
    ELSE
       CALL error_handler("matmatmul: invalid character &
            &in parameter trans_A !")
    END IF
    IF(transb=='N' .OR. transb=='n') THEN 
       n = SIZE(B,2)     ! number of columns of matrix B
    ELSE IF(transb=='T' .OR. transb=='t' .OR. &
          & transb=='C' .OR. transb=='c') THEN
       n = SIZE(B,1)     ! number of columns of transpose(B)
    ELSE
       CALL error_handler("matmatmul: invalid character &
            &in parameter trans_A !")
    END IF

    !-- Call Fortran 77 BLAS Level 3 subroutine
    CALL dgemm( transa, transb, m, n, k, alpha_in, A, lda,&
              & B, ldb, beta_in, C, ldc )

  END SUBROUTINE matmatmul
  !*************************************************************

  !*************************************************************
  SUBROUTINE matvecmul(A,x,y,trans,alpha,beta)
    !  Purpose: 
    !  Wrap around for the BLAS Level 3 subroutine DGEMV:
    !  Perform one of the matrix-vector operations:
    !
    !    y := alpha*A*x + beta*y
    !    y := alpha*transpose(A)*x + beta*y
    !    y := alpha*conjugate(transpose(A))*x + beta*y
    !
    !  where alpha and beta are scalars, x and y are vectors and
    !  A is an m by n matrix.
    !  OPTIONAL PARAMETERS:
    !  trans      (default: trans ="N")
    !  alpha,beta (default: alpha = 1.0_r8_kind = beta)
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    REAL(KIND=r8_kind),INTENT(in   ) :: x(:),A(:,:)
    REAL(KIND=r8_kind),INTENT(inout) :: y(:)
    CHARACTER(LEN=1),  INTENT(in   ),OPTIONAL :: trans
    REAL(KIND=r8_kind),INTENT(in   ),OPTIONAL :: alpha,beta
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    ! -- Variables used as arguments for BLAS subroutine
    CHARACTER(LEN=1)                 :: trans_in
    REAL(KIND=r8_kind)               :: alpha_in, beta_in
    INTEGER(KIND=i4_kind)            :: n, m, lda, incx, incy
  
    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler,dgemv
    !------------ Executable code -----------------------------------
    
    !-- First initialize argument variables for the BLAS subroutine

    ! check optinal arguments
    IF (PRESENT(trans)) THEN
       trans_in = trans
    ELSE
       trans_in = "N"      ! default is op() = id()
    END IF
    IF(PRESENT(alpha)) THEN
       alpha_in = alpha
    ELSE
       alpha_in = 1.0_r8_kind
    END IF
    IF(PRESENT(beta)) THEN
       beta_in = beta
    ELSE
       beta_in = 1.0_r8_kind
    END IF
  
    ! set dimensions of assumed shape dummy arrays
    lda = SIZE(A,1)    ! number of rows    of matrix A
    if (lda.eq.0) return
    m   = lda
    n   = SIZE(A,2)    ! number of columns of matrix A
    if (n.eq.0) return

    ! set increment for the elements of vectors x and y
    incx = 1_i4_kind
    incy = 1_i4_kind


    !-- Call Fortran 77 BLAS Level 3 subroutine
    CALL dgemv( trans_in, m, n, alpha_in, A, lda, x, incx, beta_in, y, incy )


  END SUBROUTINE matvecmul
  !*************************************************************

!--------------- End of module ----------------------------------
END MODULE linalg_module

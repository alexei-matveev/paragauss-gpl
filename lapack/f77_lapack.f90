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
MODULE F77_LAPACK
  !
  ! Interface declarations of (selected) LAPACK
  ! subroutines
  !
  ! Everybody who needs other than listed here
  ! is encouraged to extend this interface
  !
  ! Usage scheme:
  !  subroutine sub()
  !   use f77_lapack, only: dgexx
  !   implicit none
  !
  !   call dgexx(...)
  !  end subroutine sub
  !
  ! AM,05/2000

  IMPLICIT NONE

  INTERFACE

     FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       INTEGER(IP)                  :: ILAENV
       CHARACTER(LEN=*), INTENT(IN) :: NAME, OPTS
       INTEGER(IP), INTENT(IN)      :: ISPEC, N1, N2, N3, N4
     END FUNCTION ILAENV

     ! wrapper for DGETRF, DGETRS:
     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       INTEGER(IP), INTENT(IN)      :: N, NRHS, LDA, LDB
       REAL(DP),    INTENT(INOUT)   :: A(LDA,*)
       INTEGER(IP), INTENT(OUT)     :: IPIV(*)
       REAL(DP),    INTENT(INOUT)   :: B(LDB,*)
       INTEGER(IP), INTENT(OUT)     :: INFO
     END SUBROUTINE DGESV

     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       INTEGER(IP), INTENT(IN)    :: LDA, M, N
       INTEGER(IP), INTENT(OUT)   :: INFO
       INTEGER(IP), INTENT(OUT)   :: IPIV( * )
       REAL(DP),    INTENT(INOUT) :: A( LDA, * )
     END SUBROUTINE DGETRF

     SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: TRANS
       INTEGER(IP), INTENT(IN)      :: N, NRHS, LDA, LDB
       INTEGER(IP), INTENT(OUT)     :: INFO
       INTEGER(IP), INTENT(IN)      :: IPIV( * )
       REAL(DP),    INTENT(IN)      :: A( LDA, * )
       REAL(DP),    INTENT(INOUT)   :: B( LDB, * )
     END SUBROUTINE DGETRS

     SUBROUTINE ZGETRF( M, N, A, LDA, PIV, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       INTEGER(IP), INTENT(IN)    :: LDA, M, N
       INTEGER(IP), INTENT(OUT)   :: INFO
       INTEGER(IP), INTENT(OUT)   :: PIV( * )
       COMPLEX(DP), INTENT(INOUT) :: A( LDA, * )
     END SUBROUTINE ZGETRF

     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       INTEGER(IP), INTENT(IN)    :: LDA, LWORK, N
       INTEGER(IP), INTENT(OUT)   :: INFO
       INTEGER(IP), INTENT(IN)    :: IPIV(*)
       REAL(DP),    INTENT(OUT)   :: WORK(LWORK)
       REAL(DP),    INTENT(INOUT) :: A(LDA,*)
     END SUBROUTINE DGETRI

     SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       INTEGER(IP), INTENT(IN)    :: LDA, LWORK, N
       INTEGER(IP), INTENT(OUT)   :: INFO
       INTEGER(IP), INTENT(IN)    :: IPIV(*)
       COMPLEX(DP), INTENT(OUT)   :: WORK(LWORK)
       COMPLEX(DP), INTENT(INOUT) :: A(LDA,*)
     END SUBROUTINE ZGETRI

     SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,&
          &                   INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: TRANS
       INTEGER(IP), INTENT(IN)      :: NRHS, M, N, LDA, LDB, LWORK
       INTEGER(IP), INTENT(OUT)     :: INFO
       REAL(DP),    INTENT(INOUT)   :: A(LDA,*), B(LDB,*)
       REAL(DP),    INTENT(OUT)     :: WORK(*)
     END SUBROUTINE DGELS

     SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,&
          &                   INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: TRANS
       INTEGER(IP), INTENT(IN)    :: NRHS, M, N, LDA, LDB, LWORK
       INTEGER(IP), INTENT(OUT)   :: INFO
       COMPLEX(DP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
       COMPLEX(DP), INTENT(OUT)   :: WORK(*)
     END SUBROUTINE ZGELS

     SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: UPLO
       INTEGER(IP), INTENT(IN)      :: ITYPE, LDA, LDB, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       REAL(DP), INTENT(IN)         :: B(LDB,*)
       REAL(DP), INTENT(INOUT)      :: A(LDA,*)
     END SUBROUTINE DSYGST

     SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: UPLO
       INTEGER(IP), INTENT(IN)      :: ITYPE, LDA, LDB, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       COMPLEX(DP), INTENT(IN)      :: B(LDB,*)
       COMPLEX(DP), INTENT(INOUT)   :: A(LDA,*)
     END SUBROUTINE ZHEGST

     SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
       INTEGER(IP), INTENT(IN)      :: LDA, LWORK, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       REAL(DP), INTENT(INOUT)      :: A(LDA,*)
       REAL(DP), INTENT(OUT)        :: W(*)
       REAL(DP), INTENT(OUT)        :: WORK(*)
     END SUBROUTINE DSYEV

     SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  &
          &                   INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
       INTEGER(IP), INTENT(IN)      :: LDA, LWORK, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       COMPLEX(DP), INTENT(INOUT)   :: A(LDA,*)
       REAL(DP), INTENT(OUT)        :: W(*), RWORK(*)
       COMPLEX(DP), INTENT(OUT)     :: WORK(*)
     END SUBROUTINE ZHEEV

     SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
          &                   LWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
       INTEGER(IP), INTENT(IN)      :: ITYPE, LDA, LDB, LWORK, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       REAL(DP), INTENT(OUT)        :: W(*)
       REAL(DP), INTENT(INOUT)      :: A(LDA,*), B(LDB,*)
       REAL(DP), INTENT(OUT)        :: WORK(*)
     END SUBROUTINE DSYGV

     SUBROUTINE ZHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
          &                   LWORK, RWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
       INTEGER(IP), INTENT(IN)      :: ITYPE, LDA, LDB, LWORK, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       REAL(DP), INTENT(OUT)        :: W(*), RWORK(*)
       COMPLEX(DP), INTENT(INOUT)   :: A(LDA,*), B(LDB,*)
       COMPLEX(DP), INTENT(OUT)     :: WORK(*)
     END SUBROUTINE ZHEGV

     SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&
          &             WORK, LWORK, RWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBVT
       INTEGER(IP), INTENT(IN)      :: LDA, LDU, LDVT, LWORK, M, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       COMPLEX(DP), INTENT(IN)      :: A(LDA,*)
       COMPLEX(DP), INTENT(OUT)     :: U(LDU,*), VT(LDVT,*)
       REAL(DP), INTENT(OUT)        :: S(*)
       COMPLEX(DP), INTENT(OUT)     :: WORK(*)
       REAL(DP), INTENT(OUT)        :: RWORK(*)
     END SUBROUTINE ZGESVD

     SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&
          &             WORK, LWORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBVT
       INTEGER(IP), INTENT(IN)      :: LDA, LDU, LDVT, LWORK, M, N
       INTEGER(IP), INTENT(OUT)     :: INFO
       COMPLEX(DP), INTENT(IN)      :: A(LDA,*)
       COMPLEX(DP), INTENT(OUT)     :: U(LDU,*), VT(LDVT,*)
       REAL(DP), INTENT(OUT)        :: S(*)
       COMPLEX(DP), INTENT(OUT)     :: WORK(*)
     END SUBROUTINE DGESVD


     SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: UPLO
       INTEGER(IP), INTENT(IN)      :: N
       REAL(DP), INTENT(INOUT)      :: AP(*)
       INTEGER(IP), INTENT(OUT)     :: IPIV(*)
       INTEGER(IP), INTENT(OUT)     :: INFO
     END SUBROUTINE DSPTRF

     SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: UPLO
       INTEGER(IP), INTENT(IN)      :: N
       REAL(DP), INTENT(INOUT)      :: AP(*)
       INTEGER(IP), INTENT(IN)     :: IPIV(*)
       REAL(DP), INTENT(OUT)      :: WORK(*)
       INTEGER(IP), INTENT(OUT)     :: INFO
     END SUBROUTINE DSPTRI

  END INTERFACE

END MODULE F77_LAPACK

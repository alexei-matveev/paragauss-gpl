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
MODULE F77_BLAS

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
          &             BETA, C, LDC )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), intent(in) :: TRANSA, TRANSB
       INTEGER(IP), intent(in)      :: M, N, K, LDA, LDB, LDC
       REAL(DP), intent(in)         :: ALPHA, BETA
       REAL(DP), intent(in)         :: A( LDA, * ), B( LDB, * )
       REAL(DP), intent(inout)      :: C( LDC, * )
     END SUBROUTINE DGEMM

     SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
          &             BETA, C, LDC )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1), intent(in) :: TRANSA, TRANSB
       INTEGER(IP), intent(in)      :: M, N, K, LDA, LDB, LDC
       COMPLEX(DP), intent(in)      :: ALPHA, BETA
       COMPLEX(DP), intent(in)      :: A( LDA, * ), B( LDB, * )
       COMPLEX(DP), intent(inout)   :: C( LDC, * )
     END SUBROUTINE ZGEMM

     SUBROUTINE DSYRK( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       CHARACTER(LEN=1) ::   UPLO, TRANS
       INTEGER(IP)      ::   N, K, LDA, LDC
       REAL(DP)         ::   ALPHA, BETA
       REAL(DP)         ::   A( LDA,*), C(LDC,*)
     END SUBROUTINE DSYRK

     SUBROUTINE ZSYRK( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       CHARACTER(LEN=1) ::   UPLO, TRANS
       INTEGER(IP)      ::   N, K, LDA, LDC
       COMPLEX(DP)      ::   ALPHA, BETA
       COMPLEX(DP)      ::   A( LDA,*), C(LDC,*)
     END SUBROUTINE ZSYRK

     SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&
          &             B, LDB )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1) :: SIDE, UPLO, TRANSA, DIAG
       INTEGER(IP)      :: M, N, LDA, LDB
       REAL(DP)         :: ALPHA
       REAL(DP)         :: A( LDA, * ), B( LDB, * )
     END SUBROUTINE DTRSM

     SUBROUTINE ZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,&
          &             B, LDB )
       USE TYPE_MODULE, ONLY: DP=>double_precision_kind,IP=>integer_kind
       IMPLICIT NONE
       CHARACTER(LEN=1) :: SIDE, UPLO, TRANSA, DIAG
       INTEGER(IP)      :: M, N, LDA, LDB
       COMPLEX(DP)      :: ALPHA
       COMPLEX(DP)      :: A( LDA, * ), B( LDB, * )
     END SUBROUTINE ZTRSM

  END INTERFACE

END MODULE F77_BLAS

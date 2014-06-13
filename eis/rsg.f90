!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE RSG(NM,N,A,B,W,MATZ,Z,IERR)
      use type_module
      implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) N,NM,IERR,MATZ
      REAL(kind=r8_kind) A(:,:),B(:,:),W(:),Z(:,:)
      REAL(kind=r8_kind),allocatable :: W_loc(:),Z_loc(:,:),FV1(:),FV2(:)
!
!    THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
!    SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
!    TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
!    FOR THE REAL SYMMETRI!GENERALIZED EIGENPROBLEM  AX = (LAMBDA)BX.
!
!    ON INPUT
!
!       NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
!       ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!       DIMENSION STATEMENT.
!
!       N  IS THE ORDER OF THE MATRICES  A  AND  B.
!
!       A  CONTAINS A REAL SYMMETRI!MATRIX.
!
!       B  CONTAINS A POSITIVE DEFINITE REAL SYMMETRIC MATRIX.
!
!       MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
!       ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
!       ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
!
!    ON OUTPUT
!
!       W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
!
!       Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
!
!       IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
!          COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
!          AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
!
!    QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!    MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!   XS THIS VERSION DATED AUGUST 1983.
!   
!   adapted to f90 and parameter list modified TB 5/97
!   requires interface block in calling code !!!
!
!     ------------------------------------------------------------------
        allocate (W_loc(N),Z_loc(NM,N),FV1(N),FV2(N))
!
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
!
   10 CALL  REDUC(NM,N,A,B,FV2,IERR)
      IF (IERR .NE. 0) GO TO 50
      IF (MATZ .NE. 0) GO TO 20
!     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W_loc,FV1,FV2)
      CALL  TQLRAT(N,W_loc,FV2,IERR)
      GO TO 50

   20 CALL  TRED2(NM,N,A,W_loc,FV1,Z_loc)
      CALL  TQL2(NM,N,W_loc,FV1,Z_loc,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  REBAK(NM,N,B,FV2,N,Z_loc)
   50 continue
      Z = Z_loc(1:NM,1:N)
      W = W_loc(1:N)
        deallocate(W_loc,Z_loc,FV1,FV2)
      RETURN

      END
!

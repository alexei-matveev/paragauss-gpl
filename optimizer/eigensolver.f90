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
subroutine eigensolver(matrix,dimen,eigval,eigvec)
  ! Purpose: wrapper for 'evvrsp'-Routine from EISPACK.
  !          This routine is used for compatibilty with
  !          A.Voityuk. Do not change it without having a good 
  !          reason for doing so.
  ! 
  ! ----------------------------------------------------------
  use type_module
  use math_module
  integer(kind=i4_kind),intent(in)   :: dimen      
  real(kind=r8_kind),intent(inout)   :: matrix(dimen,dimen)
  real(kind=r8_kind),intent(inout)   :: eigval(dimen)
  real(kind=r8_kind),intent(inout)   :: eigvec(dimen,dimen)
  !================================================================
  ! End of public interface of module
  !================================================================
  ! -----------------------------------------------------------
  real(kind=r8_kind)             :: b(dimen,9)
  real(kind=r8_kind),allocatable :: help_mat(:,:)
  integer(kind=i4_kind)          :: iwork(dimen)
  integer(kind=i4_kind)          :: ierr,alloc_stat

  allocate(help_mat(dimen,dimen),STAT=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       ("eigensolver: allocation (1) failed")
  help_mat = matrix

  eigval = zero
  eigvec=zero
  iwork = 0_i4_kind
  b = zero
  call evvrsp(dimen,dimen,dimen,matrix,b,iwork,eigval,eigvec,ierr)

  
  matrix=help_mat
  deallocate(help_mat,STAT=alloc_stat)
  if (alloc_stat/=0) call error_handler &
       (" eigensolver: deallocation (1) failed")
  if(ierr/=0) then
     if (ierr<0) then
        write(*,*)"eigensolver: Iteration for eigevector ",ierr," failed"
     else
        write(*,*)"eigensolver: Iteration for eigenvalue ",ierr," failed"
     endif
  endif

  return
end subroutine eigensolver

!  contains 
          SUBROUTINE EVVRSP (N,NVECT,NV,A,B,IWORK,ROOT,VECT,IERR)
        use type_module
!     *
!     DIAGONALIZATION PROCEDURE USING EISPACK-BASED ROUTINES.
!     *
!     ORIGINAL AUTHOR:  S.T.ELBERT, AMES LABORATORY-USDOE, JUNE 1985.
!     REVISED: MAY 1987.
!     THE ORIGINAL ROUTINE REFERS TO A SYMMETRIC PACKED MATRIX.
!
!     PURPOSE -
!        FINDS (ALL) EIGENVALUES  AND  (SOME OR ALL) EIGENVECTORS
!        OF A SYMMETRIC MATRIX.
!
!     METHOD -
!        THE METHOD AS PRESENTED IN THIS ROUTINE CONSISTS OF FOUR STEPS:
!        FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
!        HOUSEHOLDER TECHNIQUE (ORTHOGONAL SIMILARITY TRANSFORMATIONS).
!        SECOND, THE ROOTS ARE LOCATED USING THE RATIONAL QL METHOD.
!        THIRD, THE VECTORS OF THE TRIDIAGONAL FORM ARE EVALUATED BY THE
!        INVERSE ITERATION TECHNIQUE.  VECTORS FOR DEGENERATE OR NEAR-
!        DEGENERATE ROOTS ARE FORCED TO BE ORTHOGONAL.
!        FOURTH, THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
!        ORIGINAL ARRAY.
!
!        THESE ROUTINES ARE MODIFICATIONS OF THE EISPACK 3
!        ROUTINES TRED1, TQLRAT, TINVIT AND TRBAK1
!
!     REFERENCE:
!        STEPHEN T. ELBERT,
!        THEOR. CHIM. ACTA (1987) 71:169-186.
!
!        FOR FURTHER DETAILS, SEE EISPACK USERS GUIDE, B. T. SMITH
!        ET AL, SPRINGER-VERLAG, LECTURE NOTES IN COMPUTER SCIENCE,
!        VOL. 6, 2-ND EDITION, 1976.  ANOTHER GOOD REFERENCE IS
!        THE SYMMETRIC EIGENVALUE PROBLEM BY B. N. PARLETT
!        PUBLISHED BY PRENTICE-HALL, INC., ENGLEWOOD CLIFFS, N.J. (1980)
!
!     ON ENTRY -
!        N     - INTEGER
!                ORDER OF MATRIX A.
!        NVECT - INTEGER
!                NUMBER OF VECTORS DESIRED.  0.LE.NVECT AND NVECT.LE.N
!        NV    - INTEGER
!                ROW DIMENSION OF A AND VECT IN CALLING ROUTINE. N.LE.NV
!        A     - WORKING PRECISION (NV,NV)
!                INPUT MATRIX, 2-D SYMMETRIC MATRIX
!        B     - WORKING PRECISION (NV,9)
!                SCRATCH ARRAY, 9*NV ELEMENTS
!        IWORK - INTEGER
!                SCRATCH ARRAY, NV ELEMENTS
!
!     ON EXIT  -
!        A     - DESTROYED.  NOW HOLDS REFLECTION OPERATORS.
!        ROOT  - WORKING PRECISION (N)
!                ALL EIGENVALUES IN ASCENDING ORDER.
!        VECT  - WORKING PRECISION (NV,NVECT)
!                EIGENVECTORS FOR ROOT(1), ..., ROOT(NVECT).
!        IERR  - INTEGER
!                = 0 IF NO ERROR DETECTED,
!                = K IF ITERATION FOR K-TH EIGENVALUE FAILED,
!                =-K IF ITERATION FOR K-TH EIGENVECTOR FAILED.
!
!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit real(kind=r8_kind) (a-h,o-z)
      implicit integer(kind=i4_kind) (i-n)
!RAY  DIMENSION A(NV,NV),B(NV,9),ROOT(NVECT),VECT(NV,NVECT)
      DIMENSION A(NV,NV),B(NV,8),ROOT(NVECT),VECT(NV,NVECT)
      DIMENSION IWORK(NV)
!
! *** CHECK FOR SIMPLE ERRORS
      IERR = N-1
      IF(N.LE.0) THEN
         WRITE(6,910)
         WRITE(6,900) N,NVECT,NV,IERR
         STOP 'EVVRSP'
      ENDIF
      IF(NV.LT.N) THEN
         WRITE(6,920)
         WRITE(6,900) N,NVECT,NV,IERR
         STOP 'EVVRSP'
      ENDIF
      IERR = N+1
! *** REDUCE SYMMETRIC MATRIX A TO TRIDIAGONAL FORM
!     CALL ETRED3(N,LENA,A,B(1,1),B(1,2),B(1,3))
      CALL tred1_loc(NV,N,A,B(1,1),B(1,2),B(1,3))
!RAY  DO 10 I=1,N
!  10 B(I,9) = B(I,2)
! *** FIND ALL EIGENVALUES OF TRIDIAGONAL MATRIX
      CALL EQLRAT(N,B(1,1),B(1,2),B(1,3),ROOT,IWORK,IERR,B(1,4))
      IF(IERR.NE.0) THEN
         WRITE(6,930) IERR
         WRITE(6,900) N,NVECT,NV,IERR
         STOP 'EVVRSP'
      ENDIF
      IF(NVECT.LE.0) RETURN
! *** FIND EIGENVECTORS OF TRI-DIAGONAL MATRIX VIA INVERSE ITERATION
      IERR   = 6
      B(1,3) = 0.0_r8_kind
!RAY  IF(NVECT.EQ.N) THEN
         CALL EINVIT(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,IWORK,&
              VECT,IERR,B(1,4),B(1,5),B(1,6),B(1,7),B(1,8))
!RAY  ELSE
!        CALL TINVIT(NV,N,B(1,1),B(1,2),B(1,3),NVECT,ROOT,IWORK,
!    1               VECT,IERR,B(1,4),B(1,5),B(1,6),B(1,7),B(1,8))
!     ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,940) -IERR
         WRITE(6,900) N,NVECT,NV,IERR
         ! ATTENTION: the stop has been commented out to get the eigenvectors
         !            even if EINVIT did not converge.
         !         STOP 'EVVRSP'
      ENDIF
! *** FIND EIGENVECTORS OF SYMMETRIC MATRIX VIA BACK TRANSFORMATION
!     CALL ETRBK3(NV,N,LENA,A,NVECT,VECT)
!RAY  CALL TRBAK1(NV,N,A,B(1,9),NVECT,VECT)
      CALL TRBAK1(NV,N,A,B(1,2),NVECT,VECT)
      RETURN
  900 FORMAT(/1X,'*** EVVRSP PARAMETERS ***',&
           /1X,'***      N = ',I8,' ***',&
           /1X,'***  NVECT = ',I8,' ***',&
           /1X,'***     NV = ',I8,' ***',&
           /1X,'***   IERR = ',I8,' ***')


  910 FORMAT(1X,'VALUE OF N IS LESS THAN OR EQUAL ZERO')
  920 FORMAT(1X,'NV IS LESS THAN N')
  930 FORMAT(1X,'EQLRAT HAS FAILED TO CONVERGE FOR ROOT',I5)
  940 FORMAT(1X,'EINVIT HAS FAILED TO CONVERGE FOR VECTOR',I5)
    END SUBROUTINE EVVRSP
    !     ******************************************************************
    SUBROUTINE tred1_loc(NM,N,A,D,E,E2)
      use type_module
!
      integer(kind=i4_kind) I,J,K,L,N,II,NM,JP1
      real(kind=r8_kind) A(NM,N),D(N),E(N),E2(N)
      real(kind=r8_kind) F,G,H,SCALE
      real(kind=r8_kind) ZERO
      PARAMETER (ZERO=0.0_r8_kind)
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE tred1_loc,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
!     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
!          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
!
!     THIS VERSION DATED AUGUST 1983.
!
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(D(K))
         IF (SCALE .NE. ZERO) GO TO 140
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = ZERO
  125    CONTINUE
  130    E(I) = ZERO
         E2(I) = ZERO
         GO TO 300
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = ZERO
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
  220       E(J) = G
  240    CONTINUE
!     .......... FORM P ..........
         F = ZERO
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
         H = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
  280    CONTINUE
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
  300 CONTINUE
      RETURN
   END 
!     ******************************************************************
      SUBROUTINE EQLRAT(N,DIAG,E,E2IN,D,IND,IERR,E2)
      use type_module
!
!     AUTHORS -
!        THIS IS A MODIFICATION OF ROUTINE TQLRAT FROM EISPACK EDITION 3
!        DATED AUGUST 1983.
!        TQLRAT IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
!        ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!        THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE).
!        SOME MINOR MODIFICATIONS BY W. THIEL, JANUARY 1991.
!
!     PURPOSE -
!        THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!        TRIDIAGONAL MATRIX
!
!     METHOD -
!        RATIONAL QL
!
!     ON ENTRY -
!        N      - INTEGER
!                 THE ORDER OF THE MATRIX.
!        D      - W.P. REAL (N)
!                 CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!        E2     - W.P. REAL (N)
!                 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF
!                 THE INPUT MATRIX IN ITS LAST N-1 POSITIONS.
!                 E2(1) IS ARBITRARY.
!
!      ON EXIT -
!        D      - W.P. REAL (N)
!                 CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!                 ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!                 ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!                 THE SMALLEST EIGENVALUES.
!        E2     - W.P. REAL (N)
!                 DESTROYED.
!        IERR   - INTEGER
!                 0          FOR NORMAL RETURN,
!                 J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                            DETERMINED AFTER 30 ITERATIONS.
!
!     DIFFERENCES FROM EISPACK 3 - DUE TO S. T. ELBERT
!        G=G+B INSTEAD OF IF(G.EQ.0) G=B ; B=B/64
!        F77 BACKWARD LOOPS INSTEAD OF F66 CONSTRUCT
!        GENERIC INTRINSIC FUNCTIONS
!        ARRARY  IND  ADDED FOR USE BY EINVIT
!
!     EXTERNAL ROUTINES -
!        epsilon
!        INTRINSIC--ABS, SIGN, SQRT
!
!     MODIFICATIONS BY WALTER THIEL, JANUARY 1991.
!        REPLACE SOME GO-TO-CONSTRUCTS BY BLOCK-IF-CONSTRUCTS
!        COSMETIC CHANGES IN COMMENTS ETC.
!
      integer(kind=i4_kind) I,J,L,M,N,II,L1,IERR
      integer(kind=i4_kind) IND(N)
      real(kind=r8_kind) D(N),E(N),E2(N),DIAG(N),E2IN(N)
      real(kind=r8_kind) B,C,F,G,H,P,R,S,T,epsilon
      real(kind=r8_kind) SCALE,ZERO,ONE,TWO
      PARAMETER (ZERO = 0.0_r8_kind, ONE = 1.0_r8_kind)
      PARAMETER (TWO = 2.0_r8_kind)
      PARAMETER (SCALE= 1.0_r8_kind/64.0_r8_kind)
!
      IERR = 0
      D(1)=DIAG(1)
      IND(1) = 1
      K = 0
      ITAG = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
         D(I)=DIAG(I)
  100 E2(I-1) = E2IN(I)
      F = ZERO
      T = ZERO
      B = epsilon(ONE)
      C = B *B
      B = B * SCALE
      E2(N) = ZERO
!     .......... MAIN LOOP ...................
      DO 290 L = 1, N
         H = ABS(D(L)) + ABS(E(L))
         IF (T .LT. H) THEN
            T = H
            B = epsilon(T)
            C = B * B
            B = B * SCALE
         ENDIF
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
         M = L - 1
  110    M = M + 1
         IF (E2(M) .GT. C) GO TO 110
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS AN EXIT
!                FROM THE LOOP ..........
         IF (M .GT. K) THEN
            IF (M .NE. N) E2IN(M+1) = ZERO
            K = M
            ITAG = ITAG + 1
         ENDIF
         IF (M .EQ. L) GO TO 210
!        ....... ITERATE .......
         DO 205 J = 1, 30
!              .......... FORM SHIFT ..........
            L1 = L + 1
            S = SQRT(E2(L))
            G = D(L)
            P = (D(L1) - G) / (TWO * S)
            R = SQRT(P*P+ONE)
            D(L) = S / (P + SIGN(R,P))
            H = G - D(L)
            DO 140 I = L1, N
  140       D(I) = D(I) - H
            F = F + H
!              .......... RATIONAL QL TRANSFORMATION ..........
            G = D(M) + B
            H = G
            S = ZERO
            DO 200 I = M-1,L,-1
               P = G * H
               R = P + E2(I)
               E2(I+1) = S * R
               S = E2(I) / R
               D(I+1) = H + S * (H + D(I))
               G = D(I) - E2(I) / G   + B
               H = G * P / R
  200       CONTINUE
            E2(L) = S * G
            D(L) = H
!              .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST
            IF (H .EQ. ZERO) GO TO 210
            IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
            E2(L) = H * E2(L)
            IF (E2(L) .EQ. ZERO) GO TO 210
  205    CONTINUE
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
      IERR = L
      GO TO 1001
!     .......... CONVERGED ............
  210    P = D(L) + F
!           .......... ORDER EIGENVALUES ..........
         I = 1
         IF (L .EQ. 1) GO TO 250
            IF (P .LT. D(1)) GO TO 230
               I = L
!           .......... LOOP TO FIND ORDERED POSITION
  220          I = I - 1
               IF (P .LT. D(I)) GO TO 220
               I = I + 1
               IF (I .EQ. L) GO TO 250
  230       CONTINUE
            DO 240 II = L, I+1, -1
               D(II) = D(II-1)
               IND(II) = IND(II-1)
  240       CONTINUE
  250    CONTINUE
         D(I) = P
         IND(I) = ITAG
  290 CONTINUE
 1001 RETURN
      END
!     ******************************************************************
      SUBROUTINE EINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
      use type_module
!
!     AUTHORS-
!        THIS IS A MODIFICATION OF ROUTINE TINVIT FROM EISPACK EDITION 3
!        DATED AUGUST 1983.
!        TINVIT IS A TRANSLATION OF THE INVERSE ITERATION TECHNIQUE
!        IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
!        HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!        THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE).
!        SOME MINOR MODIFICATIONS BY W. THIEL, JANUARY 1991.
!
!     PURPOSE -
!        THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
!        SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES.
!
!     METHOD -
!        INVERSE ITERATION.
!
!     ON ENTRY -
!        NM     - INTEGER
!                 MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!                 ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE
!                 DIMENSION STATEMENT.
!        N      - INTEGER
!        D      - W.P. REAL (N)
!                 CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!        E      - W.P. REAL (N)
!                 CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!                 IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!        E2     - W.P. REAL (N)
!                 CONTAINS THE SQUARES OF CORRESPONDING ELEMENTS OF E,
!                 WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
!                 E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
!                 THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE
!                 SUM OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST
!                 CONTAIN 0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER,
!                 OR 2.0 IF THE EIGENVALUES ARE IN DESCENDING ORDER.
!                 IF TQLRAT, BISECT, TRIDIB, OR IMTQLV
!                 HAS BEEN USED TO FIND THE EIGENVALUES, THEIR
!                 OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
!        M      - INTEGER
!                 THE NUMBER OF SPECIFIED EIGENVECTORS.
!        W      - W.P. REAL (M)
!                 CONTAINS THE M EIGENVALUES IN ASCENDING
!                 OR DESCENDING ORDER.
!        IND    - INTEGER (M)
!                 CONTAINS IN FIRST M POSITIONS THE SUBMATRIX INDICES
!                 ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
!                 1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX
!                 FROM THE TOP, 2 FOR THOSE BELONGING TO THE SECOND
!                 SUBMATRIX, ETC.
!        IERR   - INTEGER (LOGICAL UNIT NUMBER)
!                 LOGICAL UNIT FOR ERROR MESSAGES
!
!     ON EXIT -
!        ALL INPUT ARRAYS ARE UNALTERED.
!        Z      - W.P. REAL (NM,M)
!                 CONTAINS THE ASSOCIATED SET OF ORTHONORMAL
!                 EIGENVECTORS. ANY VECTOR WHICH WHICH FAILS TO CONVERGE
!                 IS LEFT AS IS (BUT NORMALIZED) WHEN ITERATING STOPPED.
!        IERR   - INTEGER
!                  0      FOR NORMAL RETURN,
!                 -R      IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
!                         EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
!                         (ONLY LAST FAILURE TO CONVERGE IS REPORTED)
!
!        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
!
!        RV1    - W.P. REAL (N)
!                 DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
!        RV2    - W.P. REAL (N)
!                 SUPER(1)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
!        RV3    - W.P. REAL (N)
!                 SUPER(2)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
!        RV4    - W.P. REAL (N)
!                 ELEMENTS DEFINING L IN LU DECOMPOSITION
!        RV6    - W.P. REAL (N)
!                 APPROXIMATE EIGENVECTOR
!
!     DIFFERENCES FROM EISPACK 3 -  DUE TO S. T. ELBERT
!        EPS3 IS SCALED BY  EPSCAL  (ENHANCES CONVERGENCE, BUT
!           LOWERS ACCURACY)!
!        ONE MORE ITERATION (MINIMUM 2) IS PERFORMED AFTER CONVERGENCE
!           (ENHANCES ACCURACY)!
!        REPLACE LOOP WITH PYTHAG WITH SINGLE CALL TO DNRM2!
!        IF NOT CONVERGED, USE PERFORMANCE INDEX TO DECIDE ON ERROR
!           VALUE SETTING, BUT DO NOT STOP!
!        L.U. FOR ERROR MESSAGES PASSED THROUGH IERR
!        USE PARAMETER STATEMENTS AND GENERIC INTRINSIC FUNCTIONS
!        USE LEVEL 1 BLAS
!        USE IF-THEN-ELSE TO CLARIFY LOGIC
!        LOOP OVER SUBSPACES MADE INTO DO LOOP.
!        LOOP OVER INVERSE ITERATIONS MADE INTO DO LOOP
!        ZERO ONLY REQUIRED PORTIONS OF OUTPUT VECTOR
!
!        EXTERNAL ROUTINES -
!        epsilon
!        BLAS(1)--DASUM, DAXPY, DDOT, DNRM2, DSCAL
!        INTRINSIC FUNCTIONS - ABS, MAX, SQRT
!
!     MODIFICATIONS BY WALTER THIEL, JANUARY 1991.
!        REMOVE BLAS(1) REFERENCES--DASUM,DAXPY,DDOT,DSCAL
!           AND REPLACE THEM BY FORTRAN CODE
!           LOOK FOR 'CBLAS' TO SEE THE REPLACEMENTS
!        COSMETIC CHANGES IN COMMENT CARDS ETC.
!
!
      LOGICAL CONVGD
      integer(kind=i4_kind) GROUP,I,IERR,ITS,J,JJ,M,N,NM
      integer(kind=i4_kind) P,Q,R,S,SUBMAT,TAG
      integer(kind=i4_kind) IND(M)
      real(kind=r8_kind) D(N),E(N),E2(N),W(M),Z(NM,M)
      real(kind=r8_kind) RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      real(kind=r8_kind) ANORM,EPS2,EPS3,EPS4,NORM,ORDER,RHO,U,UK,V
      real(kind=r8_kind) X0,X1,XU
      real(kind=r8_kind) EPSCAL,GRPTOL,HUNDRD,ONE,TEN,ZERO
      real(kind=r8_kind) epsilon, ESTPI1, dnrm2_loc
!BLAS real(kind=r8_kind) epsilon, ESTPI1, DASUM, DDOT, dnrm2_loc
!
      PARAMETER (ZERO = 0.0_r8_kind, ONE = 1.0_r8_kind)
      PARAMETER( GRPTOL = 0.001_r8_kind)
      PARAMETER (EPSCAL = 0.5_r8_kind, HUNDRD = 100.0_r8_kind)
      PARAMETER( TEN = 10.0_r8_kind)
!
  001 FORMAT(' EIGENVECTOR ROUTINE EINVIT DID NOT CONVERGE FOR VECTOR'&
           ,I5,'.  NORM =',1PE10.2,' PERFORMANCE INDEX =',E10.2/&
           ' (AN ERROR HALT WILL OCCUR IF THE PI IS GREATER THAN 100)')
!
!-----------------------------------------------------------------------
!
      LUEMSG = IERR
      IERR = 0
      X0 = ZERO
      UK = ZERO
      NORM = ZERO
      EPS2 = ZERO
      EPS3 = ZERO
      EPS4 = ZERO
      GROUP = 0
      TAG = 0
      ORDER = ONE - E2(1)
      Q = 0
      DO 930 SUBMAT = 1, N
         P = Q + 1
!        .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
         DO 120 Q = P, N-1
            IF (E2(Q+1) .EQ. ZERO) GO TO 140
  120    CONTINUE
         Q = N
!        .......... FIND VECTORS BY INVERSE ITERATION ..........
  140    CONTINUE
         TAG = TAG + 1
         ANORM = ZERO
         S = 0
         DO 920 R = 1, M
            IF (IND(R) .NE. TAG) GO TO 920
            ITS = 1
            X1 = W(R)
            IF (S .NE. 0) GO TO 510
!        .......... CHECK FOR ISOLATED ROOT ..........
            XU = ONE
            IF (P .EQ. Q) THEN
               RV6(P) = ONE
               CONVGD = .TRUE.
               GO TO 860
            END IF
            NORM = ABS(D(P))
            DO 500 I = P+1, Q
               NORM = MAX( NORM, ABS(D(I)) + ABS(E(I)) )
  500       CONTINUE
!        .......... EPS2 IS THE CRITERION FOR GROUPING,
!                   EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                   ROOTS ARE MODIFIED BY EPS3,
!                   EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .........
            EPS2 = GRPTOL * NORM
            EPS3 = EPSCAL * epsilon(NORM)
            UK = Q - P + 1
            EPS4 = UK * EPS3
            UK = EPS4 / SQRT(UK)
            S = P
            GROUP = 0
            GO TO 520
!        .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510       IF (ABS(X1-X0) .GE. EPS2) THEN
!                 ROOTS ARE SEPARATE
               GROUP = 0
            ELSE
!                 ROOTS ARE CLOSE
               GROUP = GROUP + 1
               IF (ORDER * (X1 - X0) .LE. EPS3) X1 = X0 + ORDER * EPS3
            END IF
!        .......... ELIMINATION WITH INTERCHANGES AND
!                   INITIALIZATION OF VECTOR ..........
  520       CONTINUE
            U = D(P) - X1
            V = E(P+1)
            RV6(P) = UK
            DO 550 I = P+1, Q
               RV6(I) = UK
               IF (ABS(E(I)) .GT. ABS(U)) THEN
!                 EXCHANGE ROWS BEFORE ELIMINATION
!                  *** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
!                      E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY .......
                  XU = U / E(I)
                  RV4(I) = XU
                  RV1(I-1) = E(I)
                  RV2(I-1) = D(I) - X1
                  IF(I.NE.Q) RV3(I-1) = E(I+1)
                  U = V - XU * RV2(I-1)
                  V = -XU * RV3(I-1)
               ELSE
!                 STRAIGHT ELIMINATION
                  XU = E(I) / U
                  RV4(I) = XU
                  RV1(I-1) = U
                  RV2(I-1) = V
                  RV3(I-1) = ZERO
                  U = D(I) - X1 - XU * V
                  IF(I.NE.Q) V = E(I+1)
               END IF
  550       CONTINUE
            IF (ABS(U) .LE. EPS3) U = EPS3
            RV1(Q) = U
            RV2(Q) = ZERO
            RV3(Q) = ZERO
!     ........ DO INVERSE ITERATIONS
            CONVGD = .FALSE.
!            DO 800 ITS = 1, 5
            DO 800 ITS = 1, 10
               IF (ITS .EQ. 1) GO TO 600
!                    .......... FORWARD SUBSTITUTION ..........
                  IF (NORM .EQ. ZERO) THEN
                     RV6(S) = EPS4
                     S = S + 1
                     IF (S .GT. Q) S = P
                  ELSE
                     XU = EPS4 / NORM
!BLAS                CALL DSCAL (Q-P+1, XU, RV6(P), 1)
                     DO 560 I = P, Q
  560                RV6(I) = RV6(I) * XU
                  END IF
!                     ... ELIMINATION OPERATIONS ON NEXT VECTOR
                  DO 590 I = P+1, Q
                     U = RV6(I)
!                         IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
!                         WAS PERFORMED EARLIER IN THE
!                         TRIANGULARIZATION PROCESS ..........
                     IF (RV1(I-1) .EQ. E(I)) THEN
                        U = RV6(I-1)
                        RV6(I-1) = RV6(I)
                     ELSE
                        U = RV6(I)
                     END IF
                     RV6(I) = U - RV4(I) * RV6(I-1)
  590             CONTINUE
  600          CONTINUE
!           .......... BACK SUBSTITUTION
               RV6(Q) = RV6(Q) / RV1(Q)
               V = U
               U = RV6(Q)
               NORM = ABS(U)
               DO 620 I = Q-1, P, -1
                  RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
                  V = U
                  U = RV6(I)
                  NORM = NORM + ABS(U)
  620          CONTINUE
               IF (GROUP .EQ. 0) GO TO 700
!                 ....... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
!                         MEMBERS OF GROUP ..........
                  J = R
                  DO 680 JJ = 1, GROUP
  630                J = J - 1
                     IF (IND(J) .NE. TAG) GO TO 630
!BLAS                CALL DAXPY(Q-P+1, -DDOT(Q-P+1,RV6(P),1,Z(P,J),1),
!BLAS*                          Z(P,J),1,RV6(P),1)
                     XU = ZERO
                     DO 640 I = P, Q
  640                XU = XU + RV6(I) * Z(I,J)
                     DO 660 I = P, Q
  660                RV6(I) = RV6(I) - XU * Z(I,J)
  680             CONTINUE
!BLAS             NORM = DASUM(Q-P+1, RV6(P), 1)
                  NORM = ZERO
                  DO 720 I = P, Q
  720             NORM = NORM + ABS(RV6(I))
  700          CONTINUE
               IF (CONVGD) GO TO 840
               IF (NORM .GE. ONE) CONVGD = .TRUE.
  800       CONTINUE
!        .......... NORMALIZE SO THAT SUM OF SQUARES IS
!                   1 AND EXPAND TO FULL ORDER ..........
  840       CONTINUE
            XU = ONE / dnrm2_loc(Q-P+1,RV6(P),1)
  860       CONTINUE
            DO 870 I = 1, P-1
               Z(I,R) = ZERO
  870       CONTINUE
            DO 890 I = P,Q
               Z(I,R) = RV6(I) * XU
  890       CONTINUE
            DO 900 I = Q+1, N
               Z(I,R) = ZERO
  900       CONTINUE
            IF (.NOT.CONVGD) THEN
               RHO = ESTPI1(Q-P+1,X1,D(P),E(P),Z(P,R),ANORM)
               IF (RHO .GE. TEN  .AND.  LUEMSG .GT. 0)&
                    WRITE(LUEMSG,001) R,NORM,RHO
!               *** SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
               IF (RHO .GT. HUNDRD) IERR = -R
            END IF
            X0 = X1
  920    CONTINUE
         IF (Q .EQ. N) GO TO 940
  930 CONTINUE
  940 CONTINUE
      RETURN
      END
!     ******************************************************************
      SUBROUTINE TRBAK1(NM,N,A,E,M,Z)
      use type_module
!
      integer(kind=i4_kind) I,J,K,L,M,N,NM
      real(kind=r8_kind) A(NM,N),E(N),Z(NM,M)
      real(kind=r8_kind) S
      real(kind=r8_kind) ZERO
      PARAMETER (ZERO=0.0_r8_kind)
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
!     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  tred1_loc.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION BY  tred1_loc
!          IN ITS STRICT LOWER TRIANGLE.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
!
!        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!          IN ITS FIRST M COLUMNS.
!
!     ON OUTPUT
!
!        Z CONTAINS THE TRANSFORMED EIGENVECTORS
!          IN ITS FIRST M COLUMNS.
!
!     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS.
!
!     THIS VERSION DATED AUGUST 1983.
!
!     MODIFICATION BY WALTER THIEL, JANUARY 1991.
!        USE UPPER TRIANGLE OF MATRIX A IN ORDER TO LOOP OVER COLUMNS.
!        LOWER TRIANGLE IS COPIED TO UPPER TRIANGLE IN DO 105 LOOP.
!        IN THE REMAINDER OF THE CODE, A(I,J) IS CHANGED TO A(J,I) ETC.
!
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
      DO 105 I = 2, N
! FIXME: if you want directives back rm one "!", may conflict with OpenMP
!!$DIR IVDEP
!!$DIR NO_RECURRENCE
!*VDIR NODEP
!*VOCL LOOP,NOVREC
      DO 100 J = 1, I-1
  100 A(J,I) = A(I,J)
  105 CONTINUE
!     .......... MAIN LOOP ...................
      DO 140 I = 2, N
         L = I - 1
         IF (E(I) .EQ. ZERO) GO TO 140
         DO 130 J = 1, M
            S = ZERO
            DO 110 K = 1, L
  110       S = S + A(K,I) * Z(K,J)
!     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN tred1_loc.
!                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / A(L,I)) / E(I)
            DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(K,I)
  130    CONTINUE
  140 CONTINUE
  200 RETURN
      END
!     ******************************************************************
      FUNCTION epsilon (X)
      use type_module
!
!     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
!
      real(kind=r8_kind) X,epsilon
      real(kind=r8_kind) A,B,C,EPS
      real(kind=r8_kind) ZERO,ONE,THREE,FOUR
      PARAMETER (ZERO=0.0_r8_kind,ONE=1.0_r8_kind)
      PARAMETER( THREE=3.0_r8_kind,FOUR=4.0_r8_kind)
!
!     THIS PROGRAM SHOULD FUNCTION PROPERLY ON ALL SYSTEMS
!     SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
!        1.  THE BASE USED IN REPRESENTING FLOATING POINT
!            NUMBERS IS NOT A POWER OF THREE.
!        2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
!            THE ACCURACY USED IN FLOATING POINT VARIABLES
!            THAT ARE STORED IN MEMORY.
!     THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
!     FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
!     ASSUMPTION 2.
!     UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
!            A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
!            B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
!            C  IS NOT EXACTLY EQUAL TO ONE,
!            EPS  MEASURES THE SEPARATION OF 1.0 FROM
!                 THE NEXT LARGER FLOATING POINT NUMBER.
!     THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
!     ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
!
!     THIS VERSION DATED 4/6/83.
!
      A = FOUR/THREE
   10 B = A - ONE
      C = B + B + B
      EPS = ABS(C-ONE)
      IF (EPS .EQ. ZERO) GO TO 10
      epsilon = EPS*ABS(X)
      RETURN
      END
!     ******************************************************************
      FUNCTION dnrm2_loc (N,DX,INCX)
      use type_module
      integer(kind=i4_kind) NEXT
      real(kind=r8_kind) dnrm2_loc
      real(kind=r8_kind) DX(N),CUTLO,CUTHI,HITEST,SUM,XMAX
      real(kind=r8_kind) ZERO,ONE
      PARAMETER (ZERO=0.0_r8_kind,ONE=1.0_r8_kind)
!
!     EUCLIDEAN NORM OF THE VECTOR DX(N) WITH STORAGE INCREMENT INCX.
!     ***** INCX=1 IN THIS VERSION *****
!     dnrm2_loc IS ZERO FOR N.LE.0.
!
!     AUTHOR          C.L.LAWSON, 08 JAN 1978
!     MODIFICATIONS   W.THIEL, 11 JAN 1991
!
!     FOUR PHASE METHOD USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!     BRIEF OUTLINE OF ALGORITHM..
!
!     PHASE 1 SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
! *** CHECK FOR VECTOR LENGTH.
      IF(N.LE.0) THEN
         dnrm2_loc  = ZERO
         GO TO 300
      ENDIF
! *** INITIALIZATION.
      HITEST = CUTHI/N
      NEXT   = 1
      SUM    = ZERO
! *** BEGIN MAIN LOOP.
      I      = 1
   20 GO TO (30,50,70,80), NEXT
   30 IF( ABS(DX(I)) .GT. CUTLO) GO TO 85
      NEXT   = 2
      XMAX   = ZERO
!     PHASE 1. SCAN ZERO COMPONENTS.
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( ABS(DX(I)) .GT. CUTLO) GO TO 85
!     PREPARE FOR PHASE 2.
      NEXT   = 3
      XMAX   = ABS(DX(I))
      SUM    = SUM + (DX(I)/XMAX)**2
      GO TO 200
!     PHASE 2. SUM IS SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
   70 IF( ABS(DX(I)) .GT. CUTLO ) THEN
         SUM = (SUM * XMAX) * XMAX
         GO TO 85
      ENDIF
!     PHASE 4. SUM IS LARGE. SCALE TO AVOID OVERFLOW.
   80 IF( ABS(DX(I)) .GT. XMAX ) THEN
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX= ABS(DX(I))
      ELSE
         SUM = SUM + (DX(I)/XMAX)**2
      ENDIF
      GO TO 200
!     PHASE 3. COMPONENT IS .GT. CUTLO.
   85 CONTINUE
      DO 95 J=I,N
      IF(ABS(DX(J)) .LT. HITEST) THEN
!        PHASE 3. SUM IS MID-RANGE. NO SCALING.
         SUM    = SUM + DX(J)**2
      ELSE
         GO TO 100
      ENDIF
   95 CONTINUE
      dnrm2_loc  = SQRT(SUM)
      GO TO 300
!     PREPARE FOR PHASE 4.
  100 I      = J
      NEXT   = 4
      SUM    = (SUM / DX(I)) / DX(I)
      XMAX   = ABS(DX(I))
      SUM    = SUM + (DX(I)/XMAX)**2
  200 CONTINUE
      I      = I + 1
      IF ( I .LE. N ) GO TO 20
! *** END OF MAIN LOOP.
!     COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
      dnrm2_loc  = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
!     ******************************************************************
      FUNCTION ESTPI1 (N,EVAL,D,E,X,ANORM)
      use type_module
!
!     AUTHOR -
!        STEPHEN T. ELBERT (AMES LABORATORY-USDOE), 5 DEC 1986.
!        MINOR COSMETIC CHANGES BY W. THIEL, JANUARY 1991.
!
!     PURPOSE -
!        EVALUATE SYMMETRIC TRIDIAGONAL MATRIX PERFORMANCE INDEX
!        FOR 1 EIGENVECTOR
!
!     METHOD -
!        THIS ROUTINE FORMS THE 1-NORM OF THE RESIDUAL MATRIX A*X-X*EVAL
!        WHERE  A  IS A SYMMETRIC TRIDIAGONAL MATRIX STORED
!        IN THE DIAGONAL (D) AND SUB-DIAGONAL (E) VECTORS, EVAL IS THE
!        EIGENVALUE OF AN EIGENVECTOR OF  A,  NAMELY  X.
!        THIS NORM IS SCALED BY MACHINE ACCURACY FOR THE PROBLEM SIZE.
!        ALL NORMS APPEARING IN THE COMMENTS BELOW ARE 1-NORMS.
!
!     ON ENTRY -
!        N      - INTEGER
!                 THE ORDER OF THE MATRIX  A.
!        EVAL   - W.P. REAL
!                 THE EIGENVALUE CORRESPONDING TO VECTOR  X.
!        D      - W.P. REAL (N)
!                 THE DIAGONAL VECTOR OF  A.
!        E      - W.P. REAL (N)
!                 THE SUB-DIAGONAL VECTOR OF  A.
!        X      - W.P. REAL (N)
!                 AN EIGENVECTOR OF  A.
!        ANORM  - W.P. REAL
!                 THE NORM OF  A  IF IT HAS BEEN PREVIOUSLY COMPUTED.
!
!     ON EXIT -
!        ANORM  - W.P. REAL
!                 THE NORM OF  A, COMPUTED IF INITIALLY ZERO.
!        ESTPI1 - W.P. REAL
!           !!A*X-X*EVAL!! / (epsilon(10*N)*!!A!!*!!X!!);
!           WHERE epsilon(X) IS THE SMALLEST NUMBER SUCH THAT
!              X + epsilon(X) .NE. X
!
!           ESTPI1 .LT. 1 == SATISFACTORY PERFORMANCE
!                  .GE. 1 AND .LE. 100 == MARGINAL PERFORMANCE
!                  .GT. 100 == POOR PERFORMANCE
!           (SEE LECT. NOTES IN COMP. SCI. VOL.6 PP 124-125)
!
      real(kind=r8_kind) ESTPI1
      real(kind=r8_kind) ANORM,EVAL,RNORM,SIZE,XNORM
      real(kind=r8_kind) D(N), E(N), X(N)
      real(kind=r8_kind) epsilon, ONE, ZERO
      PARAMETER (ZERO = 0.0_r8_kind, ONE = 1.0_r8_kind)
!
      ESTPI1 = ZERO
      IF( N .LE. 1 ) RETURN
      SIZE = 10 * N
      IF (ANORM .EQ. ZERO) THEN
!        COMPUTE NORM OF  A
         ANORM = MAX( ABS(D(1)) + ABS(E(2))&
              ,ABS(D(N)) + ABS(E(N)))
         DO 110 I = 2, N-1
            ANORM = MAX( ANORM, ABS(E(I))+ABS(D(I))+ABS(E(I+1)))
  110    CONTINUE
         IF(ANORM .EQ. ZERO) ANORM = ONE
      END IF
!     COMPUTE NORMS OF RESIDUAL AND EIGENVECTOR
      XNORM = ABS(X(1)) + ABS(X(N))
      RNORM = ABS( (D(1)-EVAL)*X(1) + E(2)*X(2))&
           +ABS( (D(N)-EVAL)*X(N) + E(N)*X(N-1))
      DO 120 I = 2, N-1
         XNORM = XNORM + ABS(X(I))
         RNORM = RNORM + ABS(E(I)*X(I-1) + (D(I)-EVAL)*X(I)&
              + E(I+1)*X(I+1))
  120 CONTINUE
      ESTPI1 = RNORM / (epsilon(SIZE)*ANORM*XNORM)
      RETURN
      END



  !*************************************************************

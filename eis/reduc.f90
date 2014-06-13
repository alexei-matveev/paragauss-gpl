! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE REDUC(NM,N,A,B,DL,IERR)
      use type_module
        implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) I,J,K,N,I1,J1,NM,NN,IERR
      real(kind=r8_kind) A(NM,N),B(NM,N),DL(N)
      real(kind=r8_kind) X,Y
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE REDUC1,
!     NUM. MATH. 11, 99-110(1968) BY MARTIN AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
!
!     THIS SUBROUTINE REDUCES THE GENERALIZED SYMMETRIC EIGENPROBLEM
!     AX=(LAMBDA)BX, WHERE B IS POSITIVE DEFINITE, TO THE STANDARD
!     SYMMETRIC EIGENPROBLEM USING THE CHOLESKY FACTORIZATION OF B.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRICES A AND B.  IF THE CHOLESKY
!          FACTOR L OF B IS ALREADY AVAILABLE, N SHOULD BE PREFIXED
!          WITH A MINUS SIGN.
!
!        A AND B !ONTAIN THE REAL SYMMETRI! INPUT MATRICES.  ONLY THE
!          FULL UPPER TRIANGLES OF THE MATRICES NEED BE SUPPLIED.  IF
!          N IS NEGATIVE, THE STRICT LOWER TRIANGLE OF B CONTAINS,
!          INSTEAD, THE STRICT LOWER TRIANGLE OF ITS CHOLESKY FACTOR L.
!
!        DL CONTAINS, IF N IS NEGATIVE, THE DIAGONAL ELEMENTS OF L.
!
!     ON OUTPUT
!
!        A CONTAINS IN ITS FULL LOWER TRIANGLE THE FULL LOWER TRIANGLE
!          OF THE SYMMETRIC MATRIX DERIVED FROM THE REDUCTION TO THE
!          STANDARD FORM.  THE STRICT UPPER TRIANGLE OF A IS UNALTERED.
!
!        B CONTAINS IN ITS STRICT LOWER TRIANGLE THE STRICT LOWER
!          TRIANGLE OF ITS CHOLESKY FACTOR L.  THE FULL UPPER
!          TRIANGLE OF B IS UNALTERED.
!
!        DL CONTAINS THE DIAGONAL ELEMENTS OF L.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          7*N+1      IF B IS NOT POSITIVE DEFINITE.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      NN = IABS(N)
      IF (N .LT. 0) GO TO 100
!     .......... FORM L IN THE ARRAYS B AND DL ..........
      DO 80 I = 1, N
         I1 = I - 1
!
         DO 80 J = I, N
            X = B(I,J)
            IF (I .EQ. 1) GO TO 40
!
            DO 20 K = 1, I1
   20       X = X - B(I,K) * B(J,K)
!
   40       IF (J .NE. I) GO TO 60
            IF (X .LE. 0.0d0) GO TO 1000
            Y = dSQRT(X)
            DL(I) = Y
            GO TO 80
   60       B(J,I) = X / Y
   80 CONTINUE
!     .......... FORM THE TRANSPOSE OF THE UPPER TRIANGLE OF INV(L)*A
!                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  100 DO 200 I = 1, NN
         I1 = I - 1
         Y = DL(I)
!
         DO 200 J = I, NN
            X = A(I,J)
            IF (I .EQ. 1) GO TO 180
!
            DO 160 K = 1, I1
  160       X = X - B(I,K) * A(J,K)
!
  180       A(J,I) = X / Y
  200 CONTINUE
!     .......... PRE-MULTIPLY BY INV(L) AND OVERWRITE ..........
      DO 300 J = 1, NN
         J1 = J - 1
!
         DO 300 I = J, NN
            X = A(I,J)
            IF (I .EQ. J) GO TO 240
            I1 = I - 1
!
            DO 220 K = J, I1
  220       X = X - A(K,J) * B(I,K)
!
  240       IF (J .EQ. 1) GO TO 280
!
            DO 260 K = 1, J1
  260       X = X - A(J,K) * B(I,K)
!
  280       A(I,J) = X / DL(I)
  300 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
 1000 IERR = 7 * N + 1
 1001 RETURN
      END
!

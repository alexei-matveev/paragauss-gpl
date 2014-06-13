! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE TRED1(NM,N,A,D,E,E2)
      use type_module
        implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) I,J,K,L,N,II,NM,JP1
      REAL(kind=r8_kind) A(NM,N),D(N),E(N),E2(N)
      REAL(kind=r8_kind) F,G,H,SCALE
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
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
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0d0
         SCALE = 0.0d0
         IF (L .LT. 1) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + dABS(D(K))
!
         IF (SCALE .NE. 0.0d0) GO TO 140
!
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0d0
  125    CONTINUE
!
  130    E(I) = 0.0d0
         E2(I) = 0.0d0
         GO TO 300
!
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
!
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -dSIGN(dsqrt(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0d0
!
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
!
  220       E(J) = G
  240    CONTINUE
!     .......... FORM P ..........
         F = 0.0d0
!
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
!
         H = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
!
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
!
  280    CONTINUE
!
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
!
  300 CONTINUE
!
      RETURN
      END
!

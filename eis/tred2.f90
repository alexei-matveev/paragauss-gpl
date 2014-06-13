! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      use type_module
        implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) I,J,K,L,N,II,NM,JP1
      real(kind=r8_kind) A(NM,N),D(N),E(N),Z(NM,N)
      real(kind=r8_kind) F,G,H,HH,SCALE
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
!     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
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
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
!
!        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
!          PRODUCED IN THE REDUCTION.
!
!        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      DO 100 I = 1, N
!
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
!
         D(I) = A(N,I)
  100 CONTINUE
!
      IF (N .EQ. 1) GO TO 510
!     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0d0
         SCALE = 0.0d0
         IF (L .LT. 2) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE +dABS(D(K))
!
         IF (SCALE .NE. 0.0d0) GO TO 140
  130    E(I) = D(L)
!
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0d0
            Z(J,I) = 0.0d0
  135    CONTINUE
!
         GO TO 290
!
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
!
         F = D(L)
         G = -dSIGN(dSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
!     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0d0
!
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
!
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
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
         HH = F / (H + H)
!     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
!     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
!
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
!
            D(J) = Z(L,J)
            Z(I,J) = 0.0d0
  280    CONTINUE
!
  290    D(I) = H
  300 CONTINUE
!     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0d0
         H = D(I)
         IF (H .EQ. 0.0d0) GO TO 380
!
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
!
         DO 360 J = 1, L
            G = 0.0d0
!
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
!
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
!
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0d0
!
  500 CONTINUE
!
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0d0
  520 CONTINUE
!
      Z(N,N) = 1.0d0
      E(1) = 0.0d0
      RETURN
      END
!

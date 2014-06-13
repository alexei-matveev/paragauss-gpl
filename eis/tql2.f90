!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
      use type_module
        implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      real(kind=r8_kind) D(N),E(N),Z(NM,N)
      real(kind=r8_kind) C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1.
!
!        E HAS BEEN DESTROYED.
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES.
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
!
      DO 100 I = 2, N
  100 E(I-1) = E(I)
!
      F = 0.0d0
      TST1 = 0.0d0
      E(N) = 0.0d0
!
      DO 240 L = 1, N
         J = 0
         H =dABS(D(L)) + dABS(E(L))
         IF (TST1 .LT. H) TST1 = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 + dabs(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0d0 * E(L))
         R = PYTHAG(P,1.0d0)
         D(L) = E(L) / (P + dSIGN(R,P))
         D(L1) = E(L) * (P + dSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
!
         DO 140 I = L2, N
  140    D(I) = D(I) - H
!
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0d0
         C2 = C
         EL1 = E(L1)
         S = 0.0d0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + dABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
!
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
!
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
!

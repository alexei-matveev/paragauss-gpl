!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE TQLRAT (N,D,E2,IERR)
      use type_module
        implicit real(kind=r8_kind) (a-h,o-z)
!
      INTEGER(kind=i4_kind) I,J,L,M,N,II,L1,MML,IERR
      real(kind=r8_kind) D(N),E2(N)
      real(kind=r8_kind) B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
!     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
!
!     ON INPUT
!
!        N IS THE ORDER OF THE MATRIX.
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!
!        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
!          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
!
!      ON OUTPUT
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!          THE SMALLEST EIGENVALUES.
!
!        E2 HAS BEEN DESTROYED.
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
  100 E2(I-1) = E2(I)
!
      F = 0.0d0
      T = 0.0d0
      E2(N) = 0.0d0
!
      DO 290 L = 1, N
         J = 0
         H = dABS(D(L)) + dSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     .......... FORM SHIFT ..........
         L1 = L + 1
         S = dSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0d0 * S)
         R = PYTHAG(P,1.0d0)
         D(L) = S / (P + dsign(R,P))
         H = G - D(L)
!
         DO 140 I = L1, N
  140    D(I) = D(I) - H
!
         F = F + H
!     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0d0) G = B
         H = G
         S = 0.0d0
         MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0d0) G = B
            H = G * P / R
  200    CONTINUE
!
         E2(L) = S * G
         D(L) = H
!     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0d0) GO TO 210
         IF (dABS(E2(L)) .LE. dABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0d0) GO TO 130
  210    P = D(L) + F
!     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
!     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
!
  250    I = 1
  270    D(I) = P
  290 CONTINUE
!
      GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
!

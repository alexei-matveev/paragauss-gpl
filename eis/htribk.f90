! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
      use type_module
        implicit real(kind=r8_kind) (a-h,o-z)
!
      integer(kind=i4_kind) I,J,K,L,M,N,NM
      real(kind=r8_kind) AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,M),ZI(NM,M)
      real(kind=r8_kind) H,S,SI
!
!     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
!     THE ALGOL PROCEDURE TRBAK1, NUM. MATH. 11, 181-195(1968)
!     BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX HERMITIAN
!     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     REAL SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  HTRIDI.
!
!     ON INPUT
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT.
!
!        N IS THE ORDER OF THE MATRIX.
!
!        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
!          FORMATIONS USED IN THE REDU!TION BY  HTRIDI  IN THEIR
!          FULL LOWER TRIANGLES EXCEPT FOR THE DIAGONAL OF AR.
!
!        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
!
!        M IS THE NUMBER OF EIGENVECTORS TO BE BA!K TRANSFORMED.
!
!        ZR CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!          IN ITS FIRST M COLUMNS.
!
!     ON OUTPUT
!
!        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
!          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
!          IN THEIR FIRST M COLUMNS.
!
!     NOTE THAT THE LAST COMPONENT OF EACH RETURNED VECTOR
!     IS REAL AND THAT VECTOR EUCLIDEAN NORMS ARE PRESERVED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
!     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
!
!     THIS VERSION DATED AUGUST 1983.
!
!     ------------------------------------------------------------------
!
      IF (M .EQ. 0) GO TO 200
!     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
!                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
!                TRIDIAGONAL MATRIX. ..........
      DO 50 K = 1, N
!
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
!
      IF (N .EQ. 1) GO TO 200
!     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. 0.0d0) GO TO 140
!
         DO 130 J = 1, M
            S = 0.0d0
            SI = 0.0d0
!
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
!     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H
!
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE
!
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
      END
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      FUNCTION PYTHAG(A,B)
        use type_module
        real(kind=r8_kind) :: PYTHAG
!       implicit real(kind=r8_kind) (a-h,o-z)
        real(kind=r8_kind) A,B
!
!     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
!
      real(kind=r8_kind) P,R,S,T,U
      P =dMAX1(dABS(A),dABS(B))
      IF (P .EQ. 0.0d0) GO TO 20
      R = (dMIN1(dABS(A),dABS(B))/P)**2
   10 CONTINUE
         T = 4.0d0 + R
         IF (T .EQ. 4.0d0) GO TO 20
         S = R/T
         U = 1.0d0 + 2.0d0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      FUNCTION EPSLON (X)
      use type_module
      real(kind=r8_kind) :: EPSLON 
!       implicit real(kind=r8_kind) (a-h,o-z)
      real(kind=r8_kind) X
!
!     ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
!
      real(kind=r8_kind) A,B,C,EPS
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
      A = 4.0d0/3.0d0
   10 B = A - 1.0d0
      C = B + B + B
      EPS = dabs(C-1.0d0)
      IF (EPS .EQ. 0.0d0) GO TO 10
      EPSLON = EPS*dabs(X)
      RETURN
      END

!:DEC>>
!       subroutine nalloc(msize,iaddr,ier)
!       integer*8 msize,ier,iaddr
!       external malloc
!       iaddr=malloc(msize)
!       ier=iaddr
!       return
!       end

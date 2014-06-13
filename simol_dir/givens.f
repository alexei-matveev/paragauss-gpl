        SUBROUTINE GIVENS( NDIM, NROOTS, NM,AL,E,V)
C
C       Title:
C       -----      Module HQRII1 to diagonalize dense real-symmetric
C                  matrices.
C
C       Abstract:
C       --------   An input real-symmetric matrix A is tridiagonalized
C                  by the method of Householder.  The eigenvalues of
C                  the tridiaqonal matrix T are found by the QR method
C                  with origin shift.  Optionally, some or all
C                  eigenvectors of matrix T are determined by inverse
C                  iteration.  The eigenvectors of matrix A are found
C                  by left-multiplying the eigenvectors of matrix T
C                  times the orthogonal matrix which brings A to
C                  tridiagonal form.
C
C     Environment: Standard Fortran 77.
C     -----------
C
C     Copyright by Annik Vivier Bunge, Carlos F. Bunge, Yoshitaka Beppu,
C     ---------    Ichizo Ninomiya And Zdenko A. Tomasic, 1986.
C
C     Reference: A.V.Bunge and C.F.Bunge,
C     ---------  Comput. Chem. 1986, v.10, N4, p.259-268
C
C
C     Machine dependent parameter:
C     ---------------------------
C
C     MACHEP       is the smallest number which makes 1. + MACHEP
C                  greater than 1. For A VAX computer MACHEP=1.4D-17.
C
C     Parameter values:
C     ----------------
C
C     ZERO, HALF, ONE  are trivial constants.
C
C     BOUNDER      is a small number involved in deciding whether two
C                  eigenvalues are close enough to be considered
C                  degenerate while computing eigenvectors.
C
C     GOLDRT       GOLDRT=0.5D0*SQRT(5.D0)-1.D0) is a random number
C                  between 0 and 1.  Other choices may affect the
C                  accuracy of nearly degenerate eigenvectors, for
C                  better or for worse, depending on the particular
C                  case considered.
C
C     Argument values:
C     ----------------
C
C     NDIM      -  order of the matrix to be diagonalized in a given
C                  run.
C
C     IV        -  index of first wanted eigenvector.  Default is
C                  IV=1.
C
C     LV        -  index of last wanted eigenvector.  If LV .LT. IV,
C                  no eigenvectors are calculated.
C
C     IORD      -  controls the order of eigenvalues.  IORD .GE. 0
C                  will cause eigenvalues to be given in
C                  non-increasing order.  IORD .LT. 0 specifies
C                  non-decreasing order.
C
C     AL        -  real symmetric matrix of dimension at lest
C                  N1=(N*N+N)/2 given in rowwise lower triangular
C                  form, viz., AL(IJ)= AA(I,J) where IJ=J+(I*I-I)/2
C                  and AA is the usual square symmetric form.
C
C     E         -  array of dimension at least N holding the
C                  eigenvalues in the order determined by variable
C                  IORD.
C
C     NVX       -  column dimension of array V in calling program. If
C                  V is a one dimensional array in the calling
C                  program, most likely NVX=N.
C
C     V         -  array of dimension NVX*N which contains the J-th
C                  eigenvector ln positions V(IJ) starting at
C                  IJ=(J-1)*NVX+1.
C
C     SCHMID    -  =.TRUE.  for normal use.  If SCHMID=.FALSE.
C                  eigenvectors of degenerate or nearly degenerate
C                  eigenvalues are not orthogonalized among
C                  themselves.
C
C     IER       -  is a condition indicator.
C
C    Working arrays (for efficiency. working arrays must be local ones):
C    --------------
C
C     IX        -  is a working array of dimension NX, first defined
C                  by IX(I) = (I*I-I)/2 so that the lower triangle
C                  matrix AL is indexed by IJ=IX(I)+J.  Further on, IX
C                  is defined by IX(1)=0 and IX(I)=IX(I-l)+N-I+1 for
C                  I.GT.1 so that the upper triangle matrix AU is
C                  indexed by IJ=IX(I)+J.
C
C     Wl-W6     -  are six working vectors each of dimension NX.
C
C     Accuracy:
C     --------     for the largest (in absolute value) eigenvalues
C                  E(I) which are not highly degenerate, machine
C                  accuracy for MR(I),
C                       MR(I) = MODR(I)/E(I)
C                  where MODR is the modulus of the residual vector
C                  R(I),
C                       R(I) = (A*V)(I) - E(I)*V(I),
C                  and V(I) denotes the I-th eigenvector.  Up to M
C                  figures may be lost for the smallest eigenvalues
C                  (and eigenvectors), where M=
C                  ALOG10(LARGEST/SMALLEST) with LARGEST and SMALLEST
C                  being the largest and smallest (ln absolute value)
C                  eigenvalues, respectively.  In presence of highly
C                  degenerate eigenvalues, users may enter
C                  SCHMID=.FALSE., which will produce very accurate
C                  eigenvectors, but then lt becomes their
C                  responsibility to produce an orthogonal set.
C
C
C     Caveat:
C     ------       matrices formed exclusively by very small or very
C                  large matrix elements will give an arithmetic
C                  overflow.
C
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8           MACHEP
        LOGICAL   SCHMID
C ----
C
C
        PARAMETER( ZERO=0.D0, HALF=0.5D0, ONE=1.D0, BOUNDR=1.D-6,
     -             MACHEP=2.22D-16, GOLDRT=0.618033988749894D0 )
C
        DIMENSION  AL(*), E(*), V(*),IX(400)
        DIMENSION W1(400),W2(400),W3(400),W4(400),W5(400),W6(400)
C
C ** PS **
        NL = NDIM
        IVL = 1
        LVL = NROOTS
        NVXL = NM
        IORDL = -1
        SCHMID = .TRUE.
        IF( IVL.LE.0 ) IVL = 1
        IF( NL.LT.1 .OR. ( NL.GT.NVXL .AND. LVL.GE.IVL ) .OR.
     -      LVL .GT. NL )  THEN
                WRITE(*,*)' N OR LV OR NVX OUTSIDE PERMISSIBLE BOUNDS'
                IER = 129
                STOP
        ENDIF
C
C*      call vecprt( al, nroots )
C
        IER = 0
        IF( NL.EQ.1 ) THEN
                E(1) = AL(1)
                V(1) = ONE
                GOTO 700
        ENDIF
C
        IX(1) = 0
        DO J=2,NL
                IX(J) = IX(J-1) + J - 1
        ENDDO

        NM1 = NL - 1
        NVF = (LVL -1) * NVXL
        IF( NVF.LT.0 ) NVF = 0
C

        NM2 = NL - 2

        IF( NL.GT.2 ) THEN
C
C       Householder transformation.
C
                DO K=1,NM2
                        KP1 = K + 1
                        W2(K) = AL(K+IX(K))
                        SCALE = ZERO
                        DO J = KP1 ,NL
                                SCALE = ABS( AL( IX(J)+K ) ) + SCALE
                        ENDDO
                        W1(K) = AL( IX(KP1)+K )
                        IF( SCALE.GT.ZERO ) THEN
                                SCALEI = ONE/SCALE
                                SUM = ZERO
                                DO J=KP1, NL
                                        W2(J) = AL( IX(J)+K ) * SCALEI
                                        SUM = W2(J)**2 + SUM
                                ENDDO
                                S = SIGN( SQRT(SUM),W2(KP1) )
                                W1(K) = -S*SCALE
                                W2(KP1) = W2(KP1) + S
                                AL( IX(KP1)+K ) = W2(KP1)*SCALE
                                H = W2(KP1)*S
                                HUNS = (H*SCALE)*SCALE
                                HI = ONE/H
                                SUMM = ZERO
                                DO II=KP1,NL
                                        W5(II) = ZERO
                                ENDDO
                                DO I=KP1,NL
                                        IM1 = I - 1
                                        SUM = ZERO
                                        I0 = IX(I)
                                        W2I = W2(I)
                                        NREST = MOD( IM1-KP1+1,6 )
                                        DO J=KP1, KP1+NREST-1
                                                SUM = W2(J)*AL(I0+J)+SUM
                                                W5(J)=W5(J)+W2I*AL(I0+J)
                                        ENDDO
                                        DO J=KP1+NREST,IM1,6
                                                 SUM = SUM
     -                                              + W2(J)*AL(I0+J)
     -                                              + W2(J+1)*AL(I0+J+1)
     -                                              + W2(J+2)*AL(I0+J+2)
     -                                              + W2(J+3)*AL(I0+J+3)
     -                                              + W2(J+4)*AL(I0+J+4)
     -                                              + W2(J+5)*AL(I0+J+5)
                                        ENDDO

                                        DO J=KP1+NREST,IM1,6
                                          W5(J)=W5(J)+W2I*AL(I0+J)
                                          W5(J+1)=W5(J+1)+W2I*AL(I0+J+1)
                                          W5(J+2)=W5(J+2)+W2I*AL(I0+J+2)
                                          W5(J+3)=W5(J+3)+W2I*AL(I0+J+3)
                                          W5(J+4)=W5(J+4)+W2I*AL(I0+J+4)
                                          W5(J+5)=W5(J+5)+W2I*AL(I0+J+5)
                                        ENDDO
                                        W6(I) = W2I*AL(I0+I) + SUM
                                ENDDO
                                DO I=KP1,NL
                                        W1(I) = ( W5(I)+W6(I) ) * HI
                                        SUMM = W1(I)*W2(I) + SUMM
                                ENDDO
                                U = HALF * SUMM * HI
                                DO I=KP1,NL
                                        I0 = IX(I)
                                        W1(I) = W2(I)*U - W1(I)
                                        W1I = W1(I)
                                        W2I = W2(I)
                                        NREST = MOD(I-KP1+1,6)
                                        DO J=KP1,KP1+NREST-1
                                                AL(I0+J) = W2I * W1(J)
     -                                             + W1I * W2(J)
     -                                             + AL(I0+J)
                                        ENDDO
                                        DO J=KP1+NREST,I,6
                                          AL(I0+J) = W2I * W1(J)
     -                                             + W1I * W2(J)
     -                                             + AL(I0+J)
                                          AL(I0+J+1) = W2I * W1(J+1)
     -                                               + W1I * W2(J+1)
     -                                               + AL(I0+J+1)
                                          AL(I0+J+2) = W2I * W1(J+2)
     -                                               + W1I * W2(J+2)
     -                                               + AL(I0+J+2)
                                          AL(I0+J+3) = W2I * W1(J+3)
     -                                               + W1I * W2(J+3)
     -                                               + AL(I0+J+3)
                                          AL(I0+J+4) = W2I * W1(J+4)
     -                                               + W1I * W2(J+4)
     -                                               + AL(I0+J+4)
                                          AL(I0+J+5) = W2I * W1(J+5)
     -                                               + W1I * W2(J+5)
     -                                               + AL(I0+J+5)
                                        ENDDO
                                ENDDO
                        ELSE
                                HUNS = ZERO
                        ENDIF
                        AL( IX(K)+K ) = HUNS
                ENDDO
        ENDIF
C
        NM1NM1 = IX(NM1) + NM1
        NM1N = IX(NL) + NM1
        NN = NM1N + 1
        W2(NM1) = AL(NM1NM1)
        W2(NL ) = AL(NN)
        W1(NM1) = AL(NM1N)
        W1(NL ) = ZERO
        GERSCH = ABS( W2(1) ) + ABS( W1(1) )
        DO I=1,NM1
                GERSCH = MAX( ABS(W2(I+1)) + ABS(W1(I)) + ABS(W1(I+1)),
     -                   GERSCH )
        ENDDO
C
C     Trap null natrix before it is too late.
C
        IF( GERSCH.EQ.ZERO ) THEN
                WRITE(*,*) ' - NULL MATRIX IN DIAGONALIZATION ROUTINE'
                IER = 33
                GOTO 700
        ENDIF
C
        SUMD  = ZERO
        SUMCOD= ZERO
        DO I=1,NL
                SUMCOD= ABS( W1(I) ) + SUMCOD
                SUMD =  ABS( W2(I) ) + SUMD
        ENDDO
        SCALE = SUMD + SUMCOD
        SCALEI= ONE/SCALE
        DO I=1,NL
                W1(I) = W1(I)*SCALEI
                W2(I) = W2(I)*SCALEI
                W3(I) = W1(I)
                E(I) = W2(I)
                V(I+NVF) = W2(I)
        ENDDO
        EPS = SQRT( MACHEP )
        GERSCH= GERSCH*SCALEI
        DEL = GERSCH*EPS
        DELW5 = GERSCH*MACHEP
        IF( (SUMCOD .EQ. ZERO ) .OR. (SUMD/SUMCOD .GT. DEL) ) THEN
C
C       QR method with original shift.
C
                DO K=NL,2,-1
110                     L = K
120                     IF( ABS(W3(L-1)).GT.DEL ) THEN
                                L = L - 1
                                IF( L.GT.1 ) GOTO 120
                        ENDIF
                        IF( L.NE.K ) THEN
                                WW = ( E(K-1)+E(K) ) * HALF
                                R = E(K)-WW
                                Z = WW-SIGN( SQRT(W3(K-1)**2 + R*R),WW )
                                EE = E(L) - Z
                                E(L)= EE
                                FF = W3(L)
                                R = SQRT (EE*EE + FF*FF)
                                RI = ONE/R
                                C = E(L)*RI
                                S = W3(L)*RI
                                WW = E(L+1) - Z
                                E(L)=( FF*C + WW*S ) * S + EE + Z
                                E(L+1)= C*WW - S*FF
                                DO J=L+1,K-1
                                        R = SQRT( E(J)**2 + W3(J)**2 )
                                        RI = ONE/R
                                        W3(J-1) = S*R
                                        EE = E(J)*C
                                        FF = W3(J)*C
                                        C = E(J)*RI
                                        S = W3(J)*RI
                                        WW = E(J+1) - Z
                                        E(J)=( FF*C + WW*S )*S + EE + Z
                                        E(J+1)= C*WW - S*FF
                                ENDDO
                                W3(K-1) = E(K)*S
                                E(K) = E(K)*C + Z
                                GOTO 110
                        ENDIF
                ENDDO
C
C       Straight selection sort of eigenvalues.
C

                SORTER = ONE
C*
C               DO I=1,NL
C                       WRITE(*,*) 'I', I, 'E:', E(I)
C               ENDDO
C*
                IF( IORDL.LT.0 ) SORTER = -ONE
                J = NL
170             L = 1
                II= 1
                LL= 1
                DO I=2,J
                        IF( (E(I)-E(L))*SORTER .LE. ZERO ) THEN
                                L = I
                        ELSE
                                II= I
                                LL= L
                        ENDIF
                ENDDO
C
C               write(*,*) 'ii,ll', ii, ll
C
                IF( II.NE.LL ) THEN
                        WW = E(LL)
                        E(LL)= E(II)
                        E(II)= WW
                ENDIF
                J = II - 1
C               write(*,*) 'j=', j
                IF( J.GT.1 ) GOTO 170
C               write(*,*) 'ky-ky'
        ENDIF
C*      write(*,*) 'ky-ky 1'

C

        IF( SUMCOD.EQ.ZERO ) THEN
C
C       Special case of diagonal matrix
C
                IF (LVL .GE. IVL) THEN

                        IF = (IVL-1)*NVXL
                        DO I=IVL,LVL
                                DO J=1,I-1
                                        V( IF+J ) = ZERO
                                ENDDO
                                V( IF + I ) = ONE
                                DO J=I+1,NL
                                        V( IF+J ) = ZERO
                                ENDDO
                                IF  = IF + NVXL
                        ENDDO
                ENDIF

        ELSEIF (LVL .GE. IVL) THEN
C
C       Inverse iteration for eigenvectors.
C
                FN  = FLOAT(NL)
                EPS1= SQRT(FN)*EPS
                SEPS= SQRT( EPS )
                EPS2=( GERSCH * BOUNDR ) / ( FN * SEPS )
                RN  = ZERO
                RA  = EPS*GOLDRT
                IF = (IVL-2) * NVXL
                DO I=IVL,LVL
                        IF = IF + NVXL
                        DO J=1,NL
                                W3(J) = ZERO
                                W4(J) = W1(J)
                                W5(J) = V(NVF+J) - E(I)
                                RN = RN + RA
                                IF( RN.GE.EPS ) RN = RN - EPS
                                W6(J) = RN
                        ENDDO
                        DO J=1,NM1
                                IF( ABS(W5(J)).LE.ABS(W1(J)) ) THEN
                                        IF( W1(J).EQ.ZERO ) W1(J) = DEL
                                        W2(J) = -W5(J)/W1(J)
                                        W5(J) =  W1(J)
                                        T = W5(J+1)
                                        W5(J+1)=  W4(J)
                                        W4(J) =  T
                                        W3(J) =  W4(J+1)
                                        IF( W3(J).EQ.ZERO ) W3(J) = DEL
                                        W4(J+1) = ZERO
                                ELSE
                                        W2(J) = -W1(J) / W5(J)
                                ENDIF
                                W4(J+1) = W3(J) * W2(J) + W4(J+1)
                                W5(J+1) = W4(J) * W2(J) + W5(J+1)
                        ENDDO
                        IF( W5(NL).EQ.ZERO ) W5(NL) = DELW5
                        WNM15I = ONE / W5(NM1)
                        WN5I = ONE / W5(NL)
                        DO ITERE = 1,2
                                IF( ITERE.NE.1 ) THEN
                                        DO J=1,NM1
                                            IF( W5(J).EQ.W1(J) ) THEN
                                                T   = W6(J)
                                                W6(J) = W6(J+1)
                                                W6(J+1) = T
                                            ENDIF
                                            W6(J+1)=W6(J)*W2(J)+W6(J+1)
                                        ENDDO
                                ENDIF
                                W6(NL) = W6(NL )*WN5I
                                W6(NM1) = ( W6(NM1) - W6(NL)*W4(NM1) )
     -                                  * WNM15I
                                VN = MAX( ABS(W6(NL)), ABS(W6(NM1)) )
                                DO K = NM2,1,-1
                                        W6(K) = (
     -                                        W6(K) - W6(K+1)*W4(K)
     -                                      - W6(K+2)*W3(K)
     -                                          ) / W5(K)
                                        VN = MAX( ABS(W6(K)), VN )
                                ENDDO
                                S = EPS1 / VN
                                NREST = MOD(NL,6)
                                DO J=1,NREST
                                        W6(J) = S*W6(J)
                                ENDDO
                                DO J=NREST+1,NL,6
                                        W6(J) = S*W6(J)
                                        W6(J+1) = S*W6(J+1)
                                        W6(J+2) = S*W6(J+2)
                                        W6(J+3) = S*W6(J+3)
                                        W6(J+4) = S*W6(J+4)
                                        W6(J+5) = S*W6(J+5)
                                ENDDO
                        ENDDO
C
                        DO J=1,NL
                                V(IF+J) = W6(J)
                        ENDDO
                ENDDO
C
C           Back transformation of eigenvectors.
C
                IG = 1
                IF = (IVL-2)*NVXL
                DO I=IVL,LVL
                        IF  = IF + NVXL
                        DO J=1,NL
                                W6(J) = V(IF+J)
                        ENDDO
                        IM1 = I - 1
                        IF( NL.GT.2 ) THEN
                                DO J=1,NM2
                                        K = NL - J - 1
C
C                                       PUT TRIANGLE LINE K FROM AL TO W1
C
                                        K0 = (K*(K+1))/2
                                        W1( K ) = AL( K0 )
                                        DO L = K+1, NL
                                                K0 = K0 + L - 1
                                                W1( L ) = AL( K0 )
                                        ENDDO

                                        IF( W1(K).NE.ZERO ) THEN
                                                KP1 = K + 1
                                                SUM = ZERO
                                                NREST = MOD(NL-KP1+1,6)
                                                DO KK=KP1,KP1+NREST-1
                                                        SUM = SUM +
     -                                                  W1(KK)*W6(KK)
                                                ENDDO
                                                DO KK=KP1+NREST,NL,6
                                                  SUM = SUM
     -                                            + W1(KK)*W6(KK)
     -                                            + W1(KK+1)*W6(KK+1)
     -                                            + W1(KK+2)*W6(KK+2)
     -                                            + W1(KK+3)*W6(KK+3)
     -                                            + W1(KK+4)*W6(KK+4)
     -                                            + W1(KK+5)*W6(KK+5)
                                                ENDDO
                                                S = -SUM/W1(K)
                                                DO KK=KP1,KP1+NREST-1
                                                   W6(KK) = W6(KK)
     -                                               + S*W1(KK)
                                                ENDDO
                                                DO KK=KP1+NREST,NL,6
                                                   W6(KK) = W6(KK)
     -                                               + S*W1(KK)
                                                   W6(KK+1) = W6(KK+1)
     -                                               + S*W1(KK+1)
                                                   W6(KK+2) = W6(KK+2)
     -                                               + S*W1(KK+2)
                                                   W6(KK+3) = W6(KK+3)
     -                                               + S*W1(KK+3)
                                                   W6(KK+4) = W6(KK+4)
     -                                               + S*W1(KK+4)
                                                   W6(KK+5) = W6(KK+5)
     -                                               + S*W1(KK+5)
                                                ENDDO
                                        ENDIF
                                ENDDO
                        ENDIF
C
                        DO J=IG,I
                                JJ = J
                                IF( ABS(E(J)-E(I)).LT.EPS2 ) GOTO 380
                        ENDDO
380                     IG = JJ
                        NREST = MOD(NL,6)
C
                        IF( IG.NE.I .AND. SCHMID ) THEN
C
C       Degenerate eigenvalues. First, orthogonalize.
C
                                KF = (IG-2)*NVXL
                                DO K=IG,IM1
                                        KF = KF + NVXL
                                        SUM = ZERO
                                        DO J=1,NREST
                                                SUM = V(KF+J)*W6(J) +SUM
                                        ENDDO
                                        DO J=1+NREST,NL,6
                                                SUM = SUM
     -                                            + V(KF+J) * W6(J)
     -                                            + V(KF+J+1) * W6(J+1)
     -                                            + V(KF+J+2) * W6(J+2)
     -                                            + V(KF+J+3) * W6(J+3)
     -                                            + V(KF+J+4) * W6(J+4)
     -                                            + V(KF+J+5) * W6(J+5)
                                        ENDDO
                                        S = -SUM
                                        DO J=1,NREST
                                                W6(J) = S*V(KF+J)+W6(J)
                                        ENDDO
                                         DO J=NREST+1,NL,6
                                            W6(J) = S*V(KF+J) + W6(J)
                                            W6(J+1)= S*V(KF+J+1)+W6(J+1)
                                            W6(J+2)= S*V(KF+J+2)+W6(J+2)
                                            W6(J+3)= S*V(KF+J+3)+W6(J+3)
                                            W6(J+4)= S*V(KF+J+4)+W6(J+4)
                                            W6(J+5)= S*V(KF+J+5)+W6(J+5)
                                        ENDDO
                                ENDDO
                        ENDIF
C
C       Normalization.
C
                        SUM = ZERO
                        DO J=1,NREST
                                SUM = W6(J)**2 + SUM
                        ENDDO
                        DO J=NREST+1,NL,6
                                SUM = W6(J)**2 + W6(J+1)**2
     -                              + W6(J+2)**2 + W6(J+3)**2
     -                              + W6(J+4)**2 + W6(J+5)**2
     -                              + SUM
                        ENDDO
                        S = ONE / SQRT( SUM )
                        DO J=1,NREST
                                V( IF+J ) = S * W6( J )
                        ENDDO
                        DO J=NREST+1,NL,6
                                V(IF+J ) = S * W6(J)
                                V(IF+J+1) = S*W6(J+1)
                                V(IF+J+2) = S*W6(J+2)
                                V(IF+J+3) = S*W6(J+3)
                                V(IF+J+4) = S*W6(J+4)
                                V(IF+J+5) = S*W6(J+5)
                        ENDDO
                ENDDO
        ENDIF
C
        DO I=1,NL
                E(I) = SCALE * E(I)
        ENDDO
C
700     RETURN
        END


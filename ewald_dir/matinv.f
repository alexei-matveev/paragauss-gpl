      SUBROUTINE matinv (A,N,D,nn,ipv)
	implicit real*8 (A-H,O-Z)
	dimension A(nn,nn),IPV(nn,3)

      D=1.D0
      TOL=1.0D-30
      DO 10 J=1,N
   10 IPV(J,3)=0

      DO 100 I=1,N
      AMAX=0.D0
      DO 25 J=1,N
      IF(IPV(J,3).EQ.1) GO TO 25

      DO 20 K=1,N
      IF(IPV(K,3).EQ.1) GO TO 20
      IF(AMAX.GE.ABS(A(J,K))) GO TO 20
      IROW=J
      ICOLUM=K
      AMAX=ABS(A(J,K))
   20 CONTINUE

   25 CONTINUE
      IF(AMAX.LE.TOL) GO TO 300
      IPV(ICOLUM,3)=1
      IPV(I,1)=IROW
      IPV(I,2)=ICOLUM
      IF(IROW.EQ.ICOLUM) GO TO 35
      DO 30 L=1,N
      SWAP=A(IROW,L)
      A(IROW,L)=A(ICOLUM,L)
   30 A(ICOLUM,L)=SWAP
   35 PIVOT=A(ICOLUM,ICOLUM)
      D=D*PIVOT
      A(ICOLUM,ICOLUM)=1.d0
      DO 40 L=1,N
   40 A(ICOLUM,L)=A(ICOLUM,L)/PIVOT
      DO 55 L1=1,N
      IF(L1.EQ.ICOLUM) GO TO 55
      T=A(L1,ICOLUM)
      A(L1,ICOLUM)=0.D0
      DO 50 L=1,N
   50 A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
   55 CONTINUE
  100 CONTINUE
      NSWAP=0
      DO 250 I=1,N
      L=N-I+1
      IF(IPV(L,1).EQ.IPV(L,2)) GO TO 250
      JROW=IPV(L,1)
      JCOLUM=IPV(L,2)
      NSWAP=NSWAP+1
      DO 200 K=1,N
      SWAP=A(K,JROW)
      A(K,JROW)=A(K,JCOLUM)
      A(K,JCOLUM)=SWAP
  200 CONTINUE
  250 CONTINUE
c	if(mod(NSWAP,2).ne.0) D=-D
      D=D*((-1.d0)**NSWAP)
      RETURN
  300 D=0.0D0
      RETURN
      E N D

      SUBROUTINE QNSTEP(ITV,MODE,U,EIGVAL,F,S,SNORM,G,N,SMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.D0)
      PARAMETER (FACTOR=4.D0)
      PARAMETER (SMOD= 1.54d0)
c      PARAMETER (SMAX=0.25)
      PARAMETER (VLIMIT=0.5d0)
      DIMENSION U(N,N),EIGVAL(N),F(N),S(n),G(n)
C     MODIFIED QUASI-NEWTON STEP

	if(mode.eq.0) itv=0   ! geometry optimization

	if(ITV.ne.1) then
		if(EIGVAL(1).gt.0) then
			aln0=0.0d0
			daln=3.0d0
			alnl=aln0-daln
1			aln=0.0d0
			do i=1,n
				aln=aln+F(I)**2/(alnl-EIGVAL(I))
			enddo !=1,n
			if(aln.lt.alnl) then
				daln=daln+2
				alnl=aln0-daln
				goto 1
			endif
		else
			daln=0.0001
			aln0=EIGVAL(1)-daln
2			aln=0.0d0
                        do i=1,n
                                aln=aln+F(I)**2/(aln0-EIGVAL(I))
                        enddo !=1,n
                        if(aln.gt.aln0) then
                                daln=daln/10
				aln0=EIGVAL(1)-daln
                                goto 2
                        endif
                        daln=3.0d0
                        alnl=aln0-daln
3			aln=0.0d0
                        do i=1,n
                                aln=aln+F(I)**2/(alnl-EIGVAL(I))
                        enddo !=1,n
                        if(aln.lt.alnl) then
                                daln=daln+2
                                alnl=aln0-daln
                                goto 3
                        endif
		endif !EIGVAL(1).gt.0
	else
		daln=0.0001
		aln0=EIGVAL(2)-daln
4			aln=0.0d0
                        do i=1,n
                                aln=aln+F(I)**2/(aln0-EIGVAL(I))
                        enddo !=1,n
                        if(aln.gt.aln0) then
                                daln=daln/10
                                aln0=EIGVAL(2)-daln
                                goto 4
                        endif
                daln=0.0001
                alnl=EIGVAL(1)+daln
5			aln=0.0d0
                        do i=1,n
                                aln=aln+F(I)**2/(alnl-EIGVAL(I))
                        enddo !=1,n
                        if(aln.lt.alnl) then
                                daln=daln/10
                                alnl=EIGVAL(1)+daln
                                goto 5
                        endif
	endif !ITV.ne.1


	do k=1,50
		aln1=(aln0+alnl)/2.d0
		aln=0.0d0
                do i=1,n
	                aln=aln+F(I)**2/(aln1-EIGVAL(I))
                enddo !=1,n
		if(aln.lt.aln1) aln0=aln1
		if(aln.gt.aln1) alnl=aln1
	enddo !k=1,50


	write(*,*) 'QNSTEP: ALn,ALn0,ALn '
        write(*,*) aln, aln0,ALn
	
      DO I=1,N
        S(I)=ZERO
      ENDDO

      DO I=1,N
C      QUASI-NEWTON STEP P=-F(I)/EIGVAL(I) is SUBSTITUTED BY
C      MODIFIED ONE
C      
	AEV=EIGVAL(I)
          IF(MODE.gt.0)THEN
       AEV=ABS(EIGVAL(I))
       AL=(AEV-SQRT(AEV**2+ FACTOR*F(I)**2))/2.D0
	elseif(I.eq.ITV) then 
	AL=(AEV+SQRT(AEV**2+ FACTOR*F(I)**2))/2.D0
	write(*,*) 'ALp ',AL
          ELSE
        AL=aln
          ENDIF
C      TO MINIMIZE ENERGY     
       IF(I.NE.ITV.or.mode.le.0)THEN
          P=-F(I)/(AEV-AL)
	  
       ELSE
C      FOR TRANSITION VECTOR
          P= F(I)/(AEV-AL)
       ENDIF
       WRITE(6,'(5x,''I,F,P '',I5,4F9.4)')I,F(I),P
       DO J=1,N
         S(J)=S(J)+P*U(J,I)
       ENDDO
c	write(*,*) S
      ENDDO   ! I=1,N

      SNORM=ZERO
      DO I =1,N
         AG=ABS(G(I))
         SF=1.0D0-0.90 D0*AG/(SMOD+AG)
         S(I) = S(I)*SF
C
c 	IF(MODE.gt.0)THEN
c           AS = ABS(S(I))
c           IF(AS.GT.VLIMIT) THEN
c           S(I)=S(I)*VLIMIT/AS 
c          ENDIF
c 	ENDIF ! MODE.gt.0
c       WRITE(6,'(5x,''I,SF,GG  '',I5,4F12.3)')I,SF
      SNORM=SNORM+S(I)*S(I)
      ENDDO
      SNORM=SQRT(SNORM)
      IF(SNORM.GT.SMAX)THEN
         FF=SMAX/SNORM
        DO I=1,N
          S(I)=S(I)*FF        
        ENDDO
      SNORM=SMAX
         ENDIF
      RETURN
      E N D



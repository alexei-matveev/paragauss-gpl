C***********************************************************************
C
      SUBROUTINE DGELS(B,A,M,N)
C
C***********************************************************************
C
C  LINEAR EQUATION SOLVER FOR SYMMETRIC MATRICES
C
C  THIS SUBROUTINE REPLACES THE OLD DGELS SUBROUTINE BY APPROPRIATE
C  CALLS TO LINPACK SUBROUTINES. ON CONVEX COMPUTERS THESE ROUTINES
C  ARE LOCATED IN THE VECLIB LIBRARY, ON CRAY COMPUTERS IN THE
C  SCILIB LIBRARY.
C
C  HISTORY:
C
C    OH 24.06.92  DATE OF CREATION
C
C***********************************************************************
C
C  ON INPUT:
C
C    M               DIMENSION OF LINEAR EQUATION SYSTEM
C    N               NUMBER OF RIGHT HAND SIDES TO SOLVE
C    B(M,N)          RIGHT HAND SIDES
C    A(M*(M+1)/2)    SYMMETRIC MATRIX IN PACKED STORAGE
C
C  ON OUTPUT:
C
C    B(M,N)          SOLUTIONS TO THE RIGHT HAND SIDES
C
C***********************************************************************
C
      use type_module
      IMPLICIT   none
      integer(kind=i4_kind)    M,N,K,IAUX(M)

      real(kind=r8_kind)     B(M,N),A(M*(M+1)/2),FAUX(M)

      integer(kind=i4_kind)    INFO
      real(kind=r8_kind)     ONE,RCOND
C       EXTERNAL   SPPCO,SPPSL,SSPFA,SSPSL
C
C***********************************************************************
C
      ONE = 1.0_r8_kind
C      ESTIMATE CONDITION AND FACTOR MATRIX
       CALL DPPCO(A,M,RCOND,FAUX,INFO)

C       WRITE (6,*) ' DGELS/SPPCO: CONDITION NUMBER ESTIMATE: ',RCOND

       IF (INFO.NE.0) GOTO 1000    ! MATRIX NOT POSITIVE DEFINITE

       IF (ONE+RCOND.EQ.ONE) THEN
         WRITE (6,*) 
     > ' ***WARNING*** DGELS/SPPCO: MATRIX NEARLY SINGULAR!'
       ENDIF

C      SOLVE LINEAR EQUATION SYSTEM FOR EACH RIGHT HAND SIDE

       DO K=1,N
         CALL DPPSL(A,M,B(1,K))
       ENDDO

       RETURN

 1000  CONTINUE

C      MATRIX IS NOT POSITIVE DEFINITE

       WRITE (6,*) '***WARNING*** DGELS/SPPCO: MATRIX NOT',
     &             ' POSITIVE DEFINITE'

C      FACTOR MATRIX WITHOUT CONDITION ESTIMATE (NORMALLY THIS PART
C      OF THE SUBROUTINE SHOULD NEVER BE ENTERED, AS THE MATRICES
C      OF THE FIT EQUATIONS ARE POSITIVE DEFINITE)
       CALL SSPFA(A,M,IAUX,INFO)

       IF (INFO.NE.0) THEN
         WRITE (6,*)
         WRITE (6,*) ' ***ERROR*** DGELS/SSPFA: MATRIX SINGULAR'
         WRITE (6,*)
         STOP
       ENDIF

C      SOLVE LINEAR EQUATION SYSTEM FOR EACH RIGHT HAND SIDE

       DO K=1,N
         CALL DSPSL(A,M,IAUX,B(1,K))
       ENDDO

       RETURN

       E N D


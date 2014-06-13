
      SUBROUTINE VECNRM(N,VEC,VNORM)
      IMPLICIT REAL*8 (A-H,O-Z)
C     *
      DIMENSION VEC(*)
      DATA ZERO /0.0 D0/
C 
         VNORM=ZERO    
         DO I=1,N    
            VNORM=VNORM+VEC(I)**2   
         ENDDO
        VNORM=SQRT(VNORM)
      IF(VNORM.EQ.ZERO)THEN
         WRITE(6,'(//20x,''***   VNORM=0 !  ***''//)')
       RETURN
      ENDIF
      DO I=1,N    
         VEC(I)=VEC(I)/VNORM
      ENDDO  
      RETURN
      END

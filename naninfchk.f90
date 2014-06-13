!***********************************************************************
!
!                        naninfchk.f
!
!  	*****************************************************************
! 	* 								*
!	* 	Absoft Corporation 					* 
! 	*	2781 Bond Street					*
!	*	Rochester Hills, MI  48309				*
!	*								*
!	*	This file contains example code for demonstration	*
!	*	purposes only.  Absoft makes no warranty of the		* 
!	*	suitability of this code for any purpose.		*
!	*								*
!	*	In no event shall Absoft be liable for any incidental,	*
!	*	indirect, special, or consequential damages arising	*
!	*	out of the use of this code.				*
!	*								*
!	***************************************************************** 
!
! Routines to test real and double values against NaN and INF
!
!            NANCHK(X) - tests REAL*4 value X against NaN
!            DNANCHK(X) - tests REAL*8 value X against NaN
!            INFCHK(X) - tests REAL*4 value X against INF
!            DINFCHK(X) - test REAL*8 value X against INF
!
! For little endian machines (Intel x86), compile with
!
!      f77 -c -DBYTE_SWAPPED=1 naninfchk.f
!	or
!      f90 -c -DBYTE_SWAPPED=1 naninfchk.f -YBOZTYPE=INT
!
! For big endian machines (PowerPC), compile with
!
!      f77 -c naninfchk.f
!	or
!      f90 -c naninfchk.f -YBOZTYPE=INT
!
!***********************************************************************
#define BYTE_SWAPPED

        LOGICAL FUNCTION NANCHK(X)
        IMPLICIT NONE
        REAL X,Y
        INTEGER I
        EQUIVALENCE(Y,I)
        Y = X
        NANCHK = ((I .AND. z'7f80 0000') .EQ. z'7f80 0000') .AND.     &
     &           ((I .AND. z'007f ffff') .NE. z'0000 0000')
        RETURN
        END

        LOGICAL FUNCTION DNANCHK(X)
        IMPLICIT NONE
        REAL*8 X,Y
        INTEGER I(2)
        EQUIVALENCE(Y,I)
        Y = X
#ifdef BYTE_SWAPPED
        DNANCHK = ((I(2) .AND. z'7ff0 0000') .EQ. z'7ff0 0000') .AND. &
     &           (((I(2) .AND. z'000f ffff') .NE. z'0000 0000') .OR.  &
     &             (I(1) .NE. 0))
#else
        DNANCHK = ((I(1) .AND. z'7ff0 0000') .EQ. z'7ff0 0000') .AND. &
     &           (((I(1) .AND. z'000f ffff') .NE. z'0000 0000') .OR.  &
     &             (I(2) .NE. 0))
#endif
        RETURN
        END

        LOGICAL FUNCTION INFCHK(X)
        IMPLICIT NONE
        REAL X,Y
        INTEGER I
        EQUIVALENCE(Y,I)
        Y = X
        INFCHK = ((I .AND. z'7f80 0000') .EQ. z'7f80 0000') .AND.    &
     &           ((I .AND. z'007f ffff') .EQ. z'0000 0000')
        RETURN
        END

        LOGICAL FUNCTION DINFCHK(X)
        IMPLICIT NONE
        REAL*8 X,Y
        INTEGER I(2)
        EQUIVALENCE(Y,I)
        Y = X
#ifdef BYTE_SWAPPED
        DINFCHK = ((I(2) .AND. z'7ff0 0000') .EQ. z'7ff0 0000') .AND. &
     &           (((I(2) .AND. z'000f ffff') .EQ. z'0000 0000') .AND. &
     &             (I(1) .EQ. 0))
#else
        DINFCHK = ((I(1) .AND. z'7ff0 0000') .EQ. z'7ff0 0000') .AND. &
     &           (((I(1) .AND. z'000f ffff') .EQ. z'0000 0000') .AND. &
     &             (I(2) .EQ. 0))
#endif
        RETURN
        END

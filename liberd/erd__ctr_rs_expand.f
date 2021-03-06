C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
         SUBROUTINE  ERD__CTR_RS_EXPAND
     +
     +                    ( NXYZT,NRS,NTU,
     +                      NR,NS,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__CTR_RS_EXPAND
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine expands the rs contraction indices:
C
C                      x (nxyzt,r>=s,tu) --> y (nxyzt,rs,tu)
C
C
C                  Input:
C
C                      NXYZT      =  total # of cartesian monomial
C                                    quadruplets (invariant indices)
C                      NRS(TU)    =  total # of contraction index pairs
C                                    for contraction shells rs(tu)
C                      NR(S)      =  # of contractions for contraction
C                                    shells r(s)
C                      X          =  original integral batch with
C                                    contraction index restriction
C                                    r>=s
C
C                  Output:
C
C                      Y          =  new integral batch with complete
C                                    set of contraction indices rs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     INITRS
         INTEGER     N,R,S
         INTEGER     NR,NS
         INTEGER     NRS,NTU
         INTEGER     NXYZT
         INTEGER     RS,TU

         DOUBLE PRECISION  X (1:NXYZT,1:NRS,1:NTU)
         DOUBLE PRECISION  Y (1:NXYZT,1:NR,1:NS,1:NTU)
C
C
C------------------------------------------------------------------------
C
C
C             ...do the expansion.
C
C
         INITRS = NRS + 1

         DO 100 TU = 1,NTU

            RS = INITRS

            DO 200 S = NS,1,-1

               DO 300 R = NR,S+1,-1

                  RS = RS - 1

                  DO N = 1,NXYZT
                     Y (N,R,S,TU) = X (N,RS,TU)
                  END DO

                  DO N = 1,NXYZT
                     Y (N,S,R,TU) = X (N,RS,TU)
                  END DO

  300          CONTINUE

               RS = RS - 1
               DO N = 1,NXYZT
                  Y (N,S,S,TU) = X (N,RS,TU)
               END DO

  200       CONTINUE

  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END

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
         SUBROUTINE  ERD__MEMORY_ERI_DERV_BATCH
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      DER3X,DER3Y,DER3Z,
     +                      DER4X,DER4Y,DER4Z,
     +                      ALPHA,CC,
     +                      SPHERIC,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__MEMORY_ERI_DERV_BATCH
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__MEMORY_DERV_CSGTO
C  DESCRIPTION : Main operation that calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted derivative electron repulsion integrals
C                on up to four different centers between cartesian or
C                spherical gaussian type shells.
C
C
C                  Input (x = 1,2,3 and 4):
C
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2,3 and 4
C                    DERyp        =  the order of differentiation on
C                                    centers y = 1,2,3,4 with respect
C                                    to the p = x,y,z coordinates
C                    ALPHA        =  primitive exponents for csh
C                                    1,2,3,4 in that order
C                    CC           =  contraction coefficient for csh
C                                    1,2,3,4 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C
C                  Output:
C
C                    IMIN,IOPT    =  minimum/optimum integer memory
C                    ZMIN,ZOPT    =  minimum/optimum flp memory
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

         INCLUDE     'erd__tuning.inc'

         LOGICAL     SPHERIC

         INTEGER     DER1X,DER2X,DER3X,DER4X
         INTEGER     DER1Y,DER2Y,DER3Y,DER4Y
         INTEGER     DER1Z,DER2Z,DER3Z,DER4Z
         INTEGER     IMIN,IOPT
         INTEGER     NALPHA,NCOEFF
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4
         INTEGER     ZMIN,ZOPT

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
C
C
C------------------------------------------------------------------------
C
C
C             ...call memory routine.
C
C
         CALL  ERD__MEMORY_DERV_CSGTO
     +
     +              ( NALPHA,NCOEFF,
     +                NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                SHELL1,SHELL2,SHELL3,SHELL4,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
     +                ALPHA,CC,
     +                L1CACHE,NCTROW,
     +                SPHERIC,
     +
     +                          IMIN,IOPT,
     +                          ZMIN,ZOPT )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END

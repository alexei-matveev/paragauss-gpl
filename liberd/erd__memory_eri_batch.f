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
         SUBROUTINE  ERD__MEMORY_ERI_BATCH
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      ALPHA,CC,
     +                      SPHERIC,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__MEMORY_ERI_BATCH
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__MEMORY_1111_CSGTO
C                ERD__MEMORY_CSGTO
C  DESCRIPTION : Main operation that calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted electron repulsion integrals on up to
C                four different centers between cartesian or spherical
C                gaussian type shells.
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
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INCLUDE     'erd__tuning.inc'

         LOGICAL     SPHERIC

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
C             ...call special memory routine for only s- and p-type
C                integrals.
C
C
         IF (MAX0(SHELL1,SHELL2,SHELL3,SHELL4).LT.2) THEN

             CALL  ERD__MEMORY_1111_CSGTO
     +
     +                  ( NALPHA,NCOEFF,
     +                    NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                    NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                    SHELL1,SHELL2,SHELL3,SHELL4,
     +                    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                    ALPHA,CC,
     +                    L1CACHE,NCTROW,
     +
     +                              IMIN,IOPT,
     +                              ZMIN,ZOPT )
     +
     +
         ELSE

             CALL  ERD__MEMORY_CSGTO
     +
     +                  ( NALPHA,NCOEFF,
     +                    NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                    NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                    SHELL1,SHELL2,SHELL3,SHELL4,
     +                    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                    ALPHA,CC,
     +                    L1CACHE,NCTROW,
     +                    SPHERIC,
     +
     +                              IMIN,IOPT,
     +                              ZMIN,ZOPT )
     +
     +
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END

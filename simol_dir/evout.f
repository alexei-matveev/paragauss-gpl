!
! ParaGauss, a program package for high-performance computations
! of molecular systems
! Copyright (C) 2014
! T. Belling, T. Grauschopf, S. Krüger, F. Nörtemann, M. Staufer,
! M. Mayer, V. A. Nasluzov, U. Birkenheuer, A. Hu, A. V. Matveev,
! A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman, D. I. Ganyushin,
! T. Kerdcharoen, A. Woiterski, A. B. Gordienko, S. Majumder,
! M. H. i Rotllant, R. Ramakrishnan, G. Dixit, A. Nikodem, T. Soini,
! M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License version 2 as published
! by the Free Software Foundation [1].
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
      SUBROUTINE EVOUT( U,E,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C     *
C     *
      DIMENSION E(N),U(N,N)
C
C     WRITE(6,140)
      NSTEP = 8 
      IB=-(NSTEP-1)
   10 IB=IB+NSTEP
      IF(IB.GT.N) GO TO 50
      NMO=NSTEP
      IF(IB+NMO.GT.N) NMO=N-IB+1
      IQ=IB+NMO-1
      WRITE(6,150) (E(K),K=IB,IQ)
      WRITE(6,160) 
      DO 20 M=1,N
   20 WRITE(6,170) M,(U(M,L),L=IB,IQ)
      WRITE(6,180)
      GO TO 10
   50 RETURN
  140 FORMAT(//18X,'EIGENVALUES AND EIGENVECTORS'//)
  150 FORMAT(//5X,10F10.3)
  160 FORMAT(/)
  170 FORMAT(I5,10F10.3)
  180 FORMAT(///)
      END


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
#include <def.h>
!===============================================================
! Public interface of module
!===============================================================
MODULE  lan_solve_module
  !---------------------------------------------------------------
  !
  !  Purpose:
  !  Contains wrap around subroutines for the F77
  !  Lancsoz eigensolver in "laso/dnlaso.f":
  !  (1) eigenvalues and -vectors of a real, symmetrical matrix M
  !      lan_solve_main    -> wrapper for laso/dnlaso.f
  !
  !  Conceptual ideas:
  !  The overall aim is to calculate the lowest eigenvalues
  !  (and, maybe, eigenvectors) of a NxN square matrix
  !      M = eps**2 + L*T*L^T
  !  by using an iterative, Lancsoz-type eigensolver.
  !
  !  This eigensolver takes a start vector of size N (typically
  !  the diagonal elements of M) and produces a set of test
  !  eigenvectors from this start vector.
  !  The eigensolver iteratively calls an auxiliary block
  !  multiplication subroutine which multiplies the test vectors
  !  with M and use the results to produce new test eigenvectors
  !  for the next cycle until the test eigenvectors are
  !  converged towards the ``true'' eigenvectors within a certain
  !  limit.
  !
  !  This wrapper does the following:
  !  1. Allocate work memory for the davidson driver module
  !     "laso/dnlaso.f"
  !  2. Call F77 dnlaso() subroutine within "dnlaso".
  !
  !  This construction must be used to get rid of the F77 COMMON blocks.
  !
  !
  !  Input:
  !  ------
  !  Diag(M) -> calculated in "df_data_module" or "open_shell_module" and
  !   stored "global_data_module" before this module is called
  !
  !  Output:
  !  -------
  !  Eigenvectors/eigenvalues of M:
  !  The eigenvalues and eigenvectors calculated by the Davidson
  !  subroutine are stored in the "global_data_module" in the
  !  arrays - global_eigvec(length of one vector,number of eigenvectors)
  !         - global_eigval(number of eigenvalues/vectors)
  !  They are subsequently processed by other subroutines called
  !  from "main_master" after THIS module is finished.
  !
  !
  !  Module called by: main_master()
  !
  !  References:
  !  (1) dnlaso.f
  !      http://www.cfm.brown.edu/people/gk/AM117/Reading/Lanczos/
  !  (2) dnlaso.f
  !      http://www.netlib.org/laso/
  !
  !  Author: SB
  !  Date:   11/04
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  USE type_module ! type specification parameters
  USE iounitadmin_module 

  IMPLICIT NONE
  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  PUBLIC lan_solve_main

  !================================================================
  ! End of public interface of module
  !================================================================

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
  REAL   (KIND=r8_kind),POINTER :: QS(:,:)
CONTAINS



  !*************************************************************
  SUBROUTINE lan_solve_main(DIAG)
    !  Purpose: 
    !  Wrap around for the f77 DNLASO.F LANCZOS-eigensolver
    !  for symmetric matrices for iterative calculation of
    !  eigenvalues of
    !     [M] = [eps**2 - L*T*L^T]
    !  L,T          defined in global_data_module
    !  eps2         squares of KS eigenvalue differences eps
    !  debug        =T -> print info like optimum workspace etc.
    !------------ Modules used ----------------------------------
    USE phys_param_module,  ONLY: hartree2ev
    USE error_module,       ONLY: MyID
    USE debug
    USE comm, only: comm_barrier
    USE global_module
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    REAL   (KIND=r8_kind) :: DIAG(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------

    ! -- Variables used as arguments for DVDSON subroutine
    INTEGER(KIND=i4_kind) :: N,NFIG,NMVAL,NMVEC,NBLOCK
    INTEGER(KIND=i4_kind) :: MAXOP,MAXJ,NVAL,NPERM,IERR
    INTEGER(KIND=i4_kind) :: SPACE, alloc_stat, i, j
    INTEGER(KIND=i4_kind) :: IND(20)
    REAL   (KIND=r8_kind), ALLOCATABLE :: VAL(:,:),VEC(:,:),WORK(:)
    !------------ Executable code -----------------------------------

    print *,MyID," LANCSZOS:  IN "

    !** set parameters for dnlaso subroutine
    N      = size(DIAG) ! order of matrix M

    ! THE NUMBER OF DECIMAL DIGITS OF ACCURACY REQUESTED IN THE EIGENVALUES 
    NFIG   = ABS(NINT(log(gl_eig_crite)/log(10.0_r8_kind)))

    ! THE ROW DIMENSION OF VAL (10 IN ALL CASES)
    NMVAL  = min(gl_NLow,N)

    NBLOCK = 2 !!min(NMVAL,16)
    NMVEC  = N !! + NBLOCK

    ! THE MAXIMUM NUMBER OF MATRIX MULTIPLIES ALLOWED BEFORE TERMINATION
    MAXOP  = max(N,gl_MaxIt)
    ! THE MAXIMUM SIZE OF A KRYLOV SUBSPACE USED
    MAXJ   = (6*NBLOCK)**2

    ! (FIND THE global_num_eigen SMALLEST) USED IN DNLASO ONL
    NVAL  = -min(N,gl_NLow)
    ! THE NUMBER OF KNOWN EIGENPAIRS ON ENTRY (ZERO IN ALL CASES)
    NPERM = 0

    allocate(QS(NMVEC,MAXJ),stat=alloc_stat)
    ASSERT(alloc_stat == 0)
    QS = 0.0_r8_kind

    allocate(VEC(NMVEC,NMVAL),stat=alloc_stat) !!NMVEC
    ASSERT(alloc_stat==0)
    VEC = 0.0_r8_kind

    allocate(VAL(NMVAL,4),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    VAL = 0.0_r8_kind

    ! SET STARTING VECTORS
    SPACE = NBLOCK*(3*N + 2*NBLOCK) &
         + MAXJ*(3*NBLOCK + ABS(NVAL) + 6) &
         + 3*ABS(NVAL)

    allocate(WORK(SPACE),stat=alloc_stat)
    ASSERT(alloc_stat==0)
    WORK = 0.0_r8_kind
    do I=1,N
       WORK(I)=1.0_r8_kind
    end do

    ! ** call eigensolver
    CALL write_to_trace_unit("    lan_solve_main:&
         & call Lanczos iterative eigensolver...")
    CALL write_to_output_units("  lan_solve_main:&
         & call Lanczos iterative eigensolver...")

    print *,MyID,"LANCSZOS PARAMETERS:"
    print *,"       * N      = ",N
    print *,"       * NFIG   = ",NFIG
    print *,"       * NMVAL  = ",NMVAL
    print *,"       * NBLOCK = ",NBLOCK
    print *,"       * NMVEC  = ",NMVEC
    print *,"       * MAXOP  = ",MAXOP
    print *,"       * MAXJ   = ",MAXJ
    print *,"       * NVAL   = ",NVAL
    print *,"       * NPERM  = ",NPERM

    call comm_barrier()

    CALL DNLASO(OP,IOVECT,N,NVAL,NFIG,NPERM,NMVAL,VAL,NMVEC,VEC,&
         NBLOCK,MAXOP,MAXJ,WORK,IND,IERR)

    CALL write_to_trace_unit  ("    lan_solve_main: after call of DNLASO")
    CALL write_to_output_units("    lan_solve_main: after call of DNLASO")

    !!    WRITE OUTPUT
    WRITE(6,100)IERR,IND(1)
100 FORMAT(//'1 IERR =',I5,'   NUMBER OF MATRIX ACCESSES =',I6)

    WRITE(6,200)((VAL(I,J),J=1,4),I=1,NPERM)
200 FORMAT(//7X,'  EIGENVALUE ',10X,'RESIDUAL NORM',2X,&
         'VALUE ERROR',4X,'VECTOR ERROR'//(D25.15,3D15.5))

    WRITE(6,300)(J,(VEC(I,J),I=1,6),J=1,NPERM)
300 FORMAT(//' FIRST SIX COMPONENTS OF THE EIGENVECTORS'//(I5,6D12.3))

    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) '              Omega                  Omega'
    WRITE(*,*) '             [a.u.]                   [eV]'
    WRITE(*,*)
    DO I=1,NPERM
       IF(VAL(I,1) .GE. 0.0_r8_kind) THEN
          WRITE(*,'(2X,I2,2(2X,D20.10))') &
               &     I,DSQRT(VAL(I,1)),DSQRT(VAL(I,1))*hartree2ev
       ELSE
          WRITE(*,*) I,'  omega^2=',val(i,1)
       END IF
    END DO
    WRITE(*,*)

    print*,"SHAPE(VEC)    = ",SHAPE(VEC)
    print*,"SHAPE(gl_VEC) = ",SHAPE(gl_VEC)

    gl_vec = VEC

    print*,"SHAPE(VAL)    = ",SHAPE(VAL)
    print*,"SHAPE(gl_VAL) = ",SHAPE(gl_VAL)
    gl_val = VAL(:,1)

    do i = 1,10
       print *,"VAL",i," = ",gl_val(i)
    end do

    !** Finally deallocate data
    DEALLOCATE(QS,VEC,VAL,WORK,stat=alloc_stat)
    ASSERT(alloc_stat==0)

    print *,MyID," LANCSZOS: OUT "

  END SUBROUTINE lan_solve_main
  !*************************************************************

  !*************************************************************
  SUBROUTINE OP(N,M,P,Q)
    !
    ! THIS IS THE MATRIX MULTIPLY FOR THE MATRIX M.
    use eigenblock_module, only: eigenblock_mult
    use error_module,      only: MyID
    use comm, only: comm_barrier
    INTEGER(KIND=i4_kind), intent(IN   ) :: N,M
    REAL   (KIND=r8_kind), intent(IN   ) :: P(N,M)
    REAL   (KIND=r8_kind), intent(INOUT) :: Q(N,M)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    INTEGER(KIND=i4_kind) :: counter, j
    REAL   (KIND=r8_kind) :: A(N*M), D(N*M)

    call comm_barrier()

    print *,MyID," LANCSZOS: OP subroutine "

    A = RESHAPE(P,(/N*M/))

    CALL eigenblock_mult(N, M, A, D)

    counter = 1_i4_kind
    DO j=1,M
       Q(:,j) = D(counter:counter+N-1)
       counter = counter + N
    END DO

  END SUBROUTINE OP
  !*************************************************************

  !*************************************************************
  SUBROUTINE IOVECT(N,M,Q,J,K)
    !
    ! THIS SUBROUTINE STORES AND RECALLS VECTORS.  IT STORES THE VECTORS
    ! IN CORE ALTHOUGH MOST REAL PROBLEMS WOULD USE A DISK.
    REAL   (KIND=r8_kind) :: Q(N,M)
    INTEGER(KIND=i4_kind),intent(in) :: K,J,N,M
    INTEGER(KIND=i4_kind) :: L,L1,I

    IF(K.EQ.1)GO TO 30
    DO L=1,M
       L1=J-M+L
       DO I=1,N
          QS(I,L1)=Q(I,L)
       END DO
    END DO
    RETURN
30  DO L=1,M
       L1=J-M+L
       DO I=1,N
          Q(I,L)=QS(I,L1)
       END DO
    END DO
    RETURN
  END SUBROUTINE IOVECT
  !*************************************************************

  !--------------- End of module ----------------------------------
END MODULE lan_solve_module

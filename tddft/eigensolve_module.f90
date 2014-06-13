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
MODULE  eigensolve_module
!---------------------------------------------------------------
!
!  Purpose:
!  Contains wrap around subroutines for the F77
!  davidson eigensolver in "dvdson_module":
!  (1) eigenvalues and -vectors of a real, symmetrical matrix M
!      dav_solve_main    -> wrapper for dvdson.f
!
!  Conceptual ideas:
!  The overall aim is to calculate the lowest eigenvalues
!  (and, maybe, eigenvectors) of a NxN square matrix
!      M = eps**2 + L*T*L^T
!  by using an iterative, Davidson-type eigensolver.
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
!     "dvdson_module"
!  2. Call modified F77 dvdson() subroutine within "dvdson_module".
!     Main modification is that this routine no longer requires
!     an external subroutine for block matrix operations as
!     ACTUAL PARAMETER, but calls the subroutine dav_block_mult()
!     from the "dav_block_module".
!     That routine does a parallel multiplication of the testvectors
!     with M on all processors.
!
!  This construction must be used to get rid of the F77 COMMON blocks.
!
!
!  Input:
!  ------
!  Diag(M) -> calculated in "df_data_module" and stored
!             "global_data_module" before this module is called
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
!  (1) dvdson.f
!      A DAVIDSON PROGRAM FOR FINDING A FEW SELECTED EXTREME
!      EIGENPAIRS OF A LARGE, SPARSE, REAL, SYMMETRIC MATRIX.
!      Andreas Stathopoulos, Charlotte F. Fischer
!      REF. IN COMP. PHYS. COMMUN. 79 (1994) 268
!
!
!  Author: HH
!  Date:   12/98
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
  USE dvdson_module,ONLY: dvdson       ! contains davidson solver
  USE phys_param_module,  ONLY: hartree2ev
  USE constants,          ONLY: ZERO, ONE
  USE global_module
  USE comm_module

  IMPLICIT NONE
  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================


!------------ public functions and subroutines ------------------
PUBLIC eigensolve_main

!================================================================
! End of public interface of module
!================================================================

!----------------------------------------------------------------
!------------ Subroutines ---------------------------------------
CONTAINS



   !*************************************************************
  SUBROUTINE eigensolve_main(debug,DIAG)
    !  Purpose:
    !  Wrap around for the f77 DVDSON.F Davidson-eigensolver
    !  for symmetric matrices for iterative calculation of
    !  eigenvalues of
    !     [M] = [eps**2 - L*T*L^T]
    !  L,T          defined in global_data_module
    !  eps2         squares of KS eigenvalue differences eps
    !  debug        =T -> print info like optimum workspace etc.
    !------------ Modules used ----------------------------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    ! blocksize for the iterative davidson eigensolver
!!$    INTEGER(KIND=i4_kind) ::&
!!$         & global_LIM           = 60_i4_kind    ! upper limit on expand. basis
!!$
!!$    & global_mblock        = 1_i4_kind,&    ! blocksize for Davidson
!!$         & global_maxiter       = 200_i4_kind     ! max. number of iterations
!!$    ! convergence parameters for the iterative davidson eigensolver
!!$    REAL(KIND=r8_kind) ::&
!!$         & global_crite = 1.0E-10_r8_kind, & ! convergence criteria (see Doc)
!!$         & global_critc = 1.0E-10_r8_kind, & ! convergence criteria (see Doc)
!!$         & global_critr = 1.0E-10_r8_kind, & ! convergence criteria (see Doc)
!!$         & global_ortho = 1.0E-10_r8_kind    ! convergence criteria (see Doc)

    LOGICAL,INTENT(in   ),OPTIONAL :: debug
    REAL   (KIND=r8_kind) :: DIAG(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------

    ! -- Variables used as arguments for DVDSON subroutine
    LOGICAL               :: HIEND

    INTEGER(KIND=i4_kind) :: N,NUME,LIM,ILOW,IHIGH,NIV,MBLOCK,MAXITER,&
         IRWSZ,IIWSZ,NMV,NLOOPS,IERR,N1

    REAL   (KIND=r8_kind) :: CRITE,CRITC,CRITR,ORTHO

    INTEGER(KIND=i4_kind),ALLOCATABLE :: ISELEC(:),IWORK(:)
    REAL   (KIND=r8_kind),ALLOCATABLE ::            WORK(:)

    ! -- Variables used in this subroutine only
    LOGICAL               :: output_info
    INTEGER(KIND=i4_kind) :: alloc_stat, i, icount, j
    !------------ Declaration of external procedures ----------------
    EXTERNAL error_handler
    !------------ Executable code -----------------------------------

    gl_Iter = 0

    ! do we have to print output in this subroutine ?
    output_info = .false.
    IF(PRESENT(debug)) THEN
       IF(debug) output_info=.true.
    END IF

    N1 = size(DIAG)

    NUME = min(gl_NLow,N1)

    IF(gl_calcall .or. (N1 .le. gl_Nlow)) THEN
       MBLOCK = N1
       LIM    = N1
    else
       MBLOCK = min(NUME,16)
       LIM    = min(NUME+2*MBLOCK,N1) !! FIXME: should be variable!
    END IF

    !** set parameters for DVDSON subroutine
    N       = N1              ! order of matrix M
    ILOW    = 1               ! index of lowest  eigenpair to be computed
    IHIGH   = NUME            ! index of highest eigenpair to be computed
    NIV     = 0               ! Let subroutine generate inital estimates
    MAXITER = gl_MaxIt        ! max. number of iterations

    ! convergence thresholds
    CRITE   = gl_eig_crite   ! for eigenvalues
    CRITC   = gl_eig_crite   ! for coeff. of last added basis vector
    CRITR   = gl_eig_crite   ! for residual vector norms
    ORTHO   = gl_eig_crite   ! over which loss of orthogonality is assumed
    ! usually 10*CRITR; set to 1.0D3 to disable

    ! sizes of the scratch workspace
    IRWSZ   = 2*N*LIM + LIM*LIM + (NUME+10)*LIM + NUME !!+ 100_i4_kind
    IIWSZ   = 6*LIM + NUME !! + 100_i4_kind

    ! allocate work arrays needed by the DVDSON-subroutine
    ALLOCATE(ISELEC(LIM),WORK(IRWSZ),IWORK(IIWSZ),&
         & stat=alloc_stat)
    ASSERT(alloc_stat==0)
    iselec = 0
    iwork  = 0
    work   = ZERO

    ! ** call eigensolver
    if (comm_i_am_master()) then
       CALL write_to_trace_unit("    dav_solve_main:&
            & call Davidson iterative eigensolver...")
       CALL write_to_output_units("    dav_solve_main:&
            & call Davidson iterative eigensolver...")
    end if

#if 0
    print *,"N       = ",N      ,"(vector dimension)"
    print *,"LIM     = ",LIM    ,"(upper limit on dimension of expanding basis)"
    print *,"MBLOCK  = ",MBLOCK ,"(number of vectors to be targeted in each iteration)"
    print *,"NUME    = ",NUME   ,"(lowest transitions)"
    print *,"MAXITER = ",MAXITER,"(number of iterations)"
    print *,"CRITE,CRITC,CRITR,ORTHO = ",CRITE,"(accuracy)"
#endif

    !
    ! This is not the original DVDSON, but modified with matrix
    ! vector multiplicaiton procedure hardcoded:
    !
    CALL DVDSON(N,LIM,DIAG,             &
         &             ILOW,IHIGH,ISELEC,NIV,MBLOCK,   &
         &             CRITE,CRITC,CRITR,ORTHO,MAXITER,&
         &             WORK,IRWSZ,IWORK,IIWSZ,HIEND,NLOOPS,NMV,IERR)

    if (comm_i_am_master()) then
       CALL write_to_trace_unit  ("    dav_solve_main: after call of DVDSON")
       CALL write_to_output_units("    dav_solve_main: after call of DVDSON")
    end if

    ! First check if we had a fatal error and now eigenvalues
    ! were computed. If we had only to many iterations then
    ! we still want to print the result.
    if (comm_i_am_master()) then
       IF((IERR/=0) .AND. (IERR/=2048)) call eigensolve_errorhandler(IERR)
    end if

    if (comm_i_am_master()) then
       ! print some info to output units
       CALL write_to_trace_unit(&
            & "       Davidson eigensolver converged:")
       CALL write_to_trace_unit(&
            & "       Matrix accesses       :",INTE=NLOOPS)
       CALL write_to_trace_unit(&
            & "       Matrix-Vector products:",INTE=NMV)
       CALL write_to_output_units(&
            & "       Davidson eigensolver converged:")
       CALL write_to_output_units(&
            & "       Matrix accesses       :",INTE=NLOOPS)
       CALL write_to_output_units(&
            & "       Matrix-Vector products:",INTE=NMV)
    end if
    ! if requested then print debug output


    if (comm_i_am_master()) then
       WRITE(output_unit,2000) ((WORK(i), i=NUME*N+j,(N+3)*NUME,NUME),j=1,MIN(30,NUME))
2000   FORMAT(//9X,'Eigenvalues',8X,'Eigval Differences',6X,&
            &    'Residuals',//(D25.15,2D20.10))

    end if

    ! here we check again for errors:
    ! in case of a not-so-severe error the first
    ! call to the errorhandler above is skipped
    ! and first the debug output is printed before this
    ! call catches...
    IF(IERR==2048) then
       if (comm_i_am_master()) then
          CALL write_to_trace_unit("")
          CALL write_to_trace_unit("***********************************************")
          CALL write_to_trace_unit("*                !!! WARNING !!!              *")
          CALL write_to_trace_unit("* Maximum number of iterations reached before *")
          CALL write_to_trace_unit("* convergence criteria were satisfied !       *")
          CALL write_to_trace_unit("*                                             *")
          CALL write_to_trace_unit("***********************************************")
          CALL write_to_trace_unit("")
          CALL write_to_output_units("")
          CALL write_to_output_units("***********************************************")
          CALL write_to_output_units("*                !!! WARNING !!!              *")
          CALL write_to_output_units("* Maximum number of iterations reached before *")
          CALL write_to_output_units("* convergence criteria were satisfied !       *")
          CALL write_to_output_units("*                                             *")
          CALL write_to_output_units("***********************************************")
          CALL write_to_output_units("")
       end if
    end if

    icount = 0
    DO I=NUME*N+1,NUME*N+NUME
       icount = icount + 1
       gl_val(icount) = WORK(I)
    END DO

    icount=0_i4_kind
    DO i=1, NUME
       gl_vec(:,i) = WORK(1+icount:N+icount)
       icount = icount + N
    END DO

    !** Finally deallocate data
    DEALLOCATE(ISELEC,WORK,IWORK,stat=alloc_stat)
    IF(alloc_stat/=0) CALL error_handler("dav_solve_main: deallocation&
         & of work arrays failed !")

  END SUBROUTINE eigensolve_main
  !************************************************************


   !*************************************************************
   SUBROUTINE eigensolve_errorhandler(ierr)
     !  Purpose:
     !  Analyze error flag from dvdson() subroutine
     IMPLICIT NONE
     !------------ Declaration of formal parameters ---------------
     INTEGER(KIND=i4_kind),INTENT(in) :: ierr
     !** End of interface *****************************************
     !------------ Declaration of local variables -----------------
     INTEGER(KIND=i4_kind) :: my_ierr, number, j
     !------------ Declaration of external procedures ----------------
     EXTERNAL error_handler
     !------------ Executable code -----------------------------------

     CALL write_to_output_units("    dav_solve_errorhandler:&
          &ERROR after DVDSON():IERR=", INTE=IERR)
     CALL write_to_trace_unit("    dav_solve_errorhandler:&
          &ERROR after DVDSON():IERR=", INTE=IERR)

     IF(IERR<0) THEN
        CALL write_to_output_units("    ABS(IERR) eigenvalues not &
             &converged in DSPEVX")
        CALL write_to_trace_unit("    ABS(IERR) eigenvalues not &
             &converged in DSPEVX")
        CALL error_handler("    dav_solve_errorhandler: unknown error situation after&
             & DVDSON() Davidson eigensolver !")
     ELSE IF (ierr<=2048) THEN

        my_ierr = ierr
        number=1_i4_kind
        DO j=11,0,-1
           number=2**j
           IF( (my_ierr/number)>0 ) THEN
              CALL print_msg(number)
              my_ierr = my_ierr - number
           END IF
        END DO

        CALL error_handler("    dav_solve_errorhandler: error after&
             & DVDSON() Davidson eigensolver !")

     ELSE
        CALL error_handler("    dav_solve_errorhandler: unknown error situation after&
             & DVDSON() Davidson eigensolver !")
     END IF

   CONTAINS

     SUBROUTINE print_msg(msg_num)
       IMPLICIT NONE
       INTEGER,INTENT(in) :: msg_num

       SELECT CASE (msg_num)
       CASE (1)
          CALL myprint("    error    1: N < LIM")
       CASE (2)
          CALL myprint("    error    2: LIM < 1")
       CASE (4)
          CALL myprint("    error    4: ISELEC(1)<1, and no range specified")
       CASE (8)
          CALL myprint("    error    8: IHIGH>N (in range of ISELEC)")
       CASE (16)
          CALL myprint("    error   16: IHIGH<ILOW (Invalid range)")
       CASE (32)
          CALL myprint("    error   32: K>LIM (Too many wanted)")
       CASE (64)
          CALL myprint("    error   64: duplication in ISELEC")
       CASE (128)
          CALL myprint("    error  128: NUME>LIM")
       CASE (256)
          CALL myprint("    error  256: MBLOCK is out of bounds")
       CASE (512)
          CALL myprint("    error  512: IWRSZ or IWISZ is not enough")
       CASE (1024)
          CALL myprint("    error 1024: orthogonalization failed")
       CASE (2048)
          CALL myprint("    error 2048: NLOOPS>MAXITER")
       CASE default
          CALL error_handler("    dav_solve_errorhandler: unknown error message")
       END SELECT
     END SUBROUTINE print_msg

     SUBROUTINE myprint(str)
       IMPLICIT NONE
       CHARACTER(LEN=*) :: str
       CALL write_to_output_units(TRIM(str))
       CALL write_to_trace_unit(TRIM(str))
     END SUBROUTINE myprint

   END SUBROUTINE eigensolve_errorhandler
  !*************************************************************


!--------------- End of module ----------------------------------
 END MODULE eigensolve_module

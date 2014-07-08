!
! ParaGauss,  a program package  for high-performance  computations of
! molecular systems
!
! Copyright (C) 2014     T. Belling,     T. Grauschopf,     S. Krüger,
! F. Nörtemann, M. Staufer,  M. Mayer, V. A. Nasluzov, U. Birkenheuer,
! A. Hu, A. V. Matveev, A. V. Shor, M. S. K. Fuchs-Rohr, K. M. Neyman,
! D. I. Ganyushin,   T. Kerdcharoen,   A. Woiterski,  A. B. Gordienko,
! S. Majumder,     M. H. i Rotllant,     R. Ramakrishnan,    G. Dixit,
! A. Nikodem, T. Soini, M. Roderus, N. Rösch
!
! This program is free software; you can redistribute it and/or modify
! it under  the terms of the  GNU General Public License  version 2 as
! published by the Free Software Foundation [1].
!
! This program is distributed in the  hope that it will be useful, but
! WITHOUT  ANY   WARRANTY;  without  even  the   implied  warranty  of
! MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE. See  the GNU
! General Public License for more details.
!
! [1] http://www.gnu.org/licenses/gpl-2.0.html
!
! Please see the accompanying LICENSE file for further information.
!
!=====================================================================
! Public interface of module
!=====================================================================
MODULE  tddft_diag
  !-------------------------------------------------------------------
  !
  !  Purpose:
  !
  !  Detalization:
  !  --------------
  !
  !  Module called by: main_master
  !
  !
  !  References: ...
  !  Remarks   : As usual NOTHING is stored within this module.
  !              Any important information of needed by other parts
  !              parts of the program will be stored in the
  !              "global_data_module.f90"
  !
  !  Author: SB
  !  Date:   01/06
  !
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------
! define FPP_TIMERS 2
# include "def.h"
  USE type_module          ! contains standard data types
  USE iounitadmin_module   ! to open output units
  USE output_module        ! defines amount of output
  USE msgtag_module
  USE comm_module
  USE datatype, only: arrmat2
  USE debug
  USE xpack

  IMPLICIT NONE
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ public functions and subroutines ---------------------
  public :: diag_init

  INTEGER(KIND=i4_kind), PUBLIC, PARAMETER :: UP = 1, DN = 2
  !===================================================================
  ! End of public interface of module
  !===================================================================

  FPP_TIMER_DECL(diag_timer)
  FPP_TIMER_DECL(dvdson_timer)
  FPP_TIMER_DECL(output_timer)
  FPP_TIMER_DECL(dvdDiag_all)

  !------------ Declaration of constants and variables ---------------

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
CONTAINS
  !*************************************************************
  SUBROUTINE diag_init()
    !  Purpose:
    ! MASTER AND SLAVE SHOULD COME HERE
    !------------ Modules used ----------------------------------
  USE filename_module,   ONLY: outfile
    USE eigensolve_module, ONLY: eigensolve_main
    USE lan_solve_module,  ONLY: lan_solve_main
    USE result_module,     ONLY: result_main
    USE resp_util_module,  ONLY: resp_util_buildi, resp_util_buildr
    USE error_module,      ONLY: MyID
    USE nto_module
    USE global_module
    USE read_module
    USE comm, only: comm_bcast, comm_barrier
    IMPLICIT NONE
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------

    REAL   (KIND=r8_kind),ALLOCATABLE :: QM1(:,:)
    REAL   (KIND=r8_kind),ALLOCATABLE :: auxdiag(:)
    INTEGER(KIND=i4_kind)             :: status
    INTEGER(KIND=i4_kind)             :: NK, i, N1, N2
    INTEGER(KIND=i4_kind)             :: n_sp, n_ir, i_ir
    INTEGER(KIND=i4_kind)             :: DIMUP, DIMDN, DIMAL
    INTEGER(KIND=i4_kind)             :: NS,NF
    REAL   (KIND=r8_kind),ALLOCATABLE :: KERN(:,:), DSTART(:)
    REAL   (KIND=r8_kind),ALLOCATABLE :: eps(:), eta(:), diag(:)
    INTEGER(KIND=i4_kind),ALLOCATABLE :: MO_ALL(:,:),IRR(:,:)
    INTEGER(KIND=i4_kind)             :: i_sp
    ! LOGICAL                           :: it_exists
    REAL (r8_kind) :: CC(4)
    INTEGER(KIND=i4_kind)             :: DIMWR(2) ,&
         DIMMS(gl_N_spin), DIMSL(gl_N_spin)
    LOGICAL                           :: tSS
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code ------------------------------------

    print *, MyID, "diag_init: entering"

    if (comm_i_am_master()) CALL system("rm "//TRIM(outfile("nto.tmp")))

    call write_to_output_units &
         ('  diag_init: start calculation of diagonal')

   !if(gl_NTO.and.comm_i_am_master()) then
   !    inquire(file=TRIM(outfile("nto.tmp")), exist=it_exists)
   !    ! Skip TDDFT calculation because the results are on file nto.tmp
   !    if(it_exists) then
   !        print *, MyID, "diag_init: master found nto.tmp file, skipping TDDFT calculation"
   !        call nto_module_main()
   !        print *, MyID, "diag_init: nto done"
   !        ! *** DO NO MORE! ***
   !        RETURN
   !    end if
   !end if

    FPP_TIMER_START(dvdDiag_all)

    n_sp = gl_N_spin
    n_ir = gl_N_irr

    gl_S_APP_XC = .false.

    i_ir_: do i_ir = 1, n_ir

       FPP_TIMER_START(diag_timer)
       gl_ir = i_ir
       NK    = gl_N_k(i_ir)
       DIMUP = gl_N_as(i_ir, UP)
       DIMDN = gl_N_as(i_ir, DN)
       DIMAL = DIMUP + DIMDN

       DIMWR = 0
       DO i = 1, n_sp
          DIMWR(i) = (size(gl_eps(i_ir,i)%m,1))
          DIMMS(i) = gl_N_as_mastr(i_ir,i)
          DIMSL(i) = gl_N_as_slave(i_ir,i)
       END DO

!DBG::       print *,"MyID: ",comm_myindex()," For concrete ir c : ", i_ir


       call write_to_output_units &
            ('  diag_init: read of two-center Coulomb for irrep ',inte=i_ir)

       !! Coul part
       call global_alloc_M('Qlm',NK)
       !! XC part
       call global_alloc_M('glT',NK,NK)

       ALLOCATE (QM1(NK,NK),STAT=status)
       ASSERT(status==0)

       call Q_calc(i_ir,gl_Q,QM1)
       if( gl_XC ) then
          call write_to_output_units &
               ('  diag_init: read of two-center XC for irrep', inte=i_ir)
          call XC_calc     (i_ir,n_sp,gl_T,QM1)
       end if

       DEALLOCATE(QM1,STAT=status)
       ASSERT(status==0)

       call global_alloc_M('glL',DIMWR(UP),DIMWR(DN),NK)
       NS = 1
       NF = DIMUP
       ALLOCATE(DIAG(DIMAL),STAT=status)
       ASSERT(status==0)
       DIAG = 0.0
       ALLOCATE(KERN(NK,NK), STAT=status)
       ASSERT(status==0)

       call write_to_output_units &
            ('  diag_init: read of three-center Coulomb and calculation of diagonal for irrep ',inte=i_ir)
       i_sp_: do i_sp = 1, n_sp
          call write_to_output_units &
               ('  diag_init: read of three-center Coulomb for spin',inte=i_sp)
          call read_C( i_ir,i_sp,&
               gl_eta (i_ir,i_sp)%m, gl_eps(i_ir,i_sp)%m,&
               gl_N_as(i_ir,i_sp),   NK,&
               gl_what_N_as(i_ir,i_sp)%m,gl_L(i_sp)%m)

          !! HERE WOULD BE LIKE THAT: THE TWO CENTER INTEGRALS WILL **NOT**
          !! BE INCLUDED INTO CALCULATIONS FOR noRI CASE
          !! IT MEANS, THAT THE DIAGONAL FOR NORI WILL BE BASED
          !! ONLY ON COULOMB INTERCHANGE INTEGRALS
          !! THIS CASE IS HERE BECAUSE OF THE TEST OF THE
          !! NUMERICAL INSTABILITY DUE TO THE XC INTEGRALS

          if (gl_noRI .or. gl_S_App) then
             CC = (/1,0,0,0/) !! Q UPUP UPDN DNDN
          else
             if (n_sp == 1) then       !! CS
                CC = (/1,1,1,0/)
             elseif (i_sp == 1) then   !! OS
                CC = (/1,1,0,0/)
             elseif (i_sp == 2) then
                CC = (/1,0,0,1/)
             end if
          end if

          KERN = CC(1) * gl_Q      &
               + CC(2) * gl_T(1)%m &
               + CC(3) * gl_T(2)%m
          if(i_sp==2) KERN = KERN + CC(4) * gl_T(3)%m

          gl_SS_: if (gl_SS) then

             ALLOCATE(DSTART(DIMWR(i_sp)),STAT=status)
             ASSERT(status==0)
             DSTART = gl_eps(i_ir,i_sp)%m ** 2

             ALLOCATE(auxdiag(gl_N_as(i_ir,i_sp)),STAT=status)
             ASSERT(status==0)
             auxdiag = 0.0

             call diag_assembl(DSTART,gl_L(i_sp)%m,KERN,&
                  &            DIMAL,DIMMS(i_sp),DIMSL(i_sp),NK,&
                  &            auxdiag)

             call comm_bcast(auxdiag)

             DIAG(NS:NF) = auxdiag
             NS = NF + 1
             NF = DIMAL
             DEALLOCATE(auxdiag,DSTART,STAT=status)
             ASSERT(status==0)

          end if gl_SS_

       end do i_sp_

       DEALLOCATE(KERN,STAT=status)
       ASSERT(status==0)

       if (gl_LOWESTeV .ne. 0) then
          gl_NLow = NLOW(diag,gl_LOWESTeV)
          gl_calcall = .false.
       end if

       NS = MIN(gl_NLow,DIMAL)
       call global_alloc_M('val',NS)
       call global_alloc_M('vec',DIMAL,NS)

       if ((size(DIAG) .eq. 0) .or. (NS .eq. 0)) then
          deallocate(diag,STAT=status)
          ASSERT(status==0)
          call global_dealloc_M
          cycle
       end if

       if (gl_SS) then

           call write_to_output_units&
                ('  diag_init: calculation of davidson for irrep ',inte=i_ir)

          FPP_TIMER_STOP (diag_timer)
          FPP_TIMER_START(dvdson_timer)
          if (gl_lanczos) then
             call lan_solve_main(DIAG)
          else
             call eigensolve_main(.true.,DIAG)
          end if
          FPP_TIMER_STOP (dvdson_timer)
          FPP_TIMER_START(output_timer)

       end if

       ALLOCATE(   eps(DIMAL),  eta(DIMAL),&
            &   MO_ALL(DIMAL,2),IRR(DIMAL,2), &
            STAT = status)
       ASSERT(status==0)

       NS = 1
       NF = DIMUP
       do i = 1, n_sp
          call resp_util_buildr(DIMMS(i),DIMSL(i),&
               gl_eps(i_ir,i)%m,eps(NS:NF))
          call resp_util_buildr(DIMMS(i),DIMSL(i),&
               gl_eta(i_ir,i)%m,eta(NS:NF))

          call resp_util_buildi(DIMMS(i),DIMSL(i),&
               gl_MO(i_ir,i)%m(:,1),MO_ALL(NS:NF,1))
          call resp_util_buildi(DIMMS(i),DIMSL(i),&
               gl_MO(i_ir,i)%m(:,2),MO_ALL(NS:NF,2))

          call resp_util_buildi(DIMMS(i),DIMSL(i),&
               gl_IRR(i_ir,i)%m(:,1),IRR(NS:NF,1))
          call resp_util_buildi(DIMMS(i),DIMSL(i),&
               gl_IRR(i_ir,i)%m(:,2),IRR(NS:NF,2))
          NS = NF + 1
          NF = DIMAL
       end do

       if (comm_i_am_master()) then
           if(gl_SS) call result_main(i_ir,DIMAL,eps,eta,MO_ALL,IRR)
       end if

       FPP_TIMER_STOP (output_timer)

       gl_ST_: if (gl_ST) then
          tSS = gl_SS
          FPP_TIMER_START(diag_timer)
          if (tSS) gl_SS = .false.
          !! DIAG ASSEMBL
          ALLOCATE(auxdiag(DIMUP),STAT=status)
          ASSERT(status==0)
          auxdiag = 0.0

          ALLOCATE(DSTART(DIMWR(UP)),STAT=status)
          ASSERT(status==0)
          DSTART = gl_eps(i_ir,UP)%m ** 2
          !! PREVIOUSLY WASN'T USE IF noRI IS ON
          call diag_assembl(DSTART,gl_L(UP)%m,gl_T(1)%m - gl_T(2)%m, &
               &            DIMAL,DIMMS(UP),DIMSL(UP),NK,&
               &            auxdiag)


          call comm_bcast(auxdiag)
          DIAG = auxdiag
          DEALLOCATE(auxdiag,DSTART,STAT=status)
          ASSERT(status==0)

          if (size(DIAG)==0) then
             deallocate(diag,STAT=status)
             ASSERT(status==0)
             call global_dealloc_M
             cycle
          end if

          FPP_TIMER_STOP(diag_timer)
          FPP_TIMER_START(dvdson_timer)

          if (gl_lanczos) then
             call lan_solve_main(DIAG)
          else
             call eigensolve_main(.true.,DIAG)
          end if
          FPP_TIMER_STOP(dvdson_timer)

          FPP_TIMER_START(output_timer)
          if (comm_i_am_master()) then
               call result_main(i_ir,DIMAL,eps,eta,MO_ALL,IRR)
          end if
          FPP_TIMER_STOP(output_timer)
          if (tSS) gl_SS = .true.
       end if gl_ST_

       !! Successive Approximations
       if (gl_S_App) then

          if(gl_ST) call write_to_output_units("WARNING!!! SA for S->T")

!DBG:          PRINT *,"Successive Approximation... XC"
          gl_S_APP_XC = .true.

          N1 = size(gl_vec,1)
          N2 = size(gl_vec,2)
          ASSERT(N1==N2)

          !! DIAG ASSEMBL
          ALLOCATE(auxdiag(DIMUP),STAT=status)
          ASSERT(status==0)
          auxdiag = 0.0

          ALLOCATE(DSTART(DIMWR(UP)),STAT=status)
          ASSERT(status==0)
          DSTART  = 0.0
          if (gl_SS) then
             call diag_assembl(DSTART,gl_L(UP)%m,gl_T(1)%m + gl_T(2)%m, &
                  &            DIMAL,DIMMS(UP),DIMSL(UP),NK,     &
                  &            auxdiag)
          else
             call diag_assembl(DSTART,gl_L(UP)%m,gl_T(1)%m - gl_T(2)%m, &
                  &            DIMAL,DIMMS(UP),DIMSL(UP), NK,    &
                  &            auxdiag)
          end if
          DEALLOCATE(DSTART,STAT=status)
          ASSERT(status==0)

          call comm_bcast(auxdiag)
          DIAG = auxdiag
          DEALLOCATE(auxdiag,STAT=status)
          ASSERT(status==0)

          if (size(DIAG)==0) then
             deallocate(diag,STAT=status)
             ASSERT(status==0)
             call global_dealloc_M
             cycle
          end if

          if (gl_lanczos) then
             call lan_solve_main(DIAG)
          else
             call eigensolve_main(.true.,DIAG)
          end if

          if (comm_i_am_master()) then
               call result_main(i_ir,DIMAL,eps,eta,MO_ALL,IRR)
          end if
       end if

       FPP_TIMER_START(output_timer)
       deallocate(diag,eps,eta,MO_ALL,IRR,STAT=status)
       ASSERT(status==0)

       call global_dealloc_M
       FPP_TIMER_STOP (output_timer)
    end do i_ir_

    FPP_TIMER_START(output_timer)
    call global_dealloc

    if (comm_i_am_master()) call missing_irreps
    FPP_TIMER_STOP (output_timer)

    FPP_TIMER_STOP(dvdDiag_all)

#ifdef FPP_TIMERS
    block
       real (r8_kind) ::  tt
       tt = FPP_TIMER_VALUE(diag_timer)
       WRITE (*,*) MyID, "TDDFT DIAG TIMER "
       WRITE (*,*) "   * DIAGONAL        ", tt
       tt = FPP_TIMER_VALUE(dvdson_timer)
       WRITE (*,*) "   * DAVIDSON/FULLXC ", tt
       tt = FPP_TIMER_VALUE(output_timer)
       WRITE (*,*) "   * OUTPUT          ", tt
       tt = FPP_TIMER_VALUE(dvdDiag_all)
       WRITE (*,*) "   * SUMMARY         ", tt
    end block
#endif

    print *, MyID, "diag_init: entering the barrier BEFORE"
    call comm_barrier()
    print *, MyID, "diag_init: passing the barrier BEFORE"

    print *, MyID, "diag_init: about to start NTO"
    ! MH: entry point for nto_module
    if(gl_NTO.and.comm_i_am_master()) then
        print *, MyID, "diag_init: master calling NTO"
        call nto_module_main()
        print *, MyID, "diag_init: master finished NTO"
    end if
    print *, MyID, "diag_init: continue after NTO"

    print *, MyID, "diag_init: entering the barrier AFTER"
    call comm_barrier()
    print *, MyID, "diag_init: passing the barrier AFTER"

  END SUBROUTINE diag_init
  !*************************************************************

  function NLOW(x,a) result(i)
    real(kind = r8_kind)    :: x(:)
    real(kind = r8_kind)    :: a
    integer(kind = i4_kind) :: i

    integer(kind = i4_kind) :: j

    i = 0
    do j = 1, size(x)
       if (x(j) < a) i = i + 1
    end do

  end function NLOW

  !*************************************************************
  subroutine diag_assembl(DIAG_AUX, AX, KERN, NAS, NM, NS, NK, DIAG)
    !  Purpose: ..
    !------------ Modules used ----------------------------------
    use global_module
    use comm_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    REAL   (KIND=r8_kind),INTENT(INOUT) :: DIAG_AUX(:)
    REAL   (KIND=r8_kind),INTENT(IN)    :: KERN(:,:), AX(:,:)
    INTEGER(KIND=i4_kind),INTENT(IN)    :: NM, NS, NK
    REAL   (KIND=r8_kind),INTENT(INOUT) :: DIAG(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    REAL(KIND=r8_kind),ALLOCATABLE      :: DA(:)
    INTEGER(KIND=i4_kind) :: status, i, ias, j, i_proc, &
         NS_UP, NF_UP, Nas, NN
    !------------ Executable code ------------------------------------

    if (comm_i_am_master()) then
       NN = NM
    else
       NN = NS
    end if

    do ias = 1, NN
       do i = 1, NK
          do j = 1, NK
             DIAG_AUX(ias) = AX  (ias,i) &
                  &        * KERN(i,j) &
                  &        * AX  (ias,j) &
                  &        + DIAG_AUX(ias)
          end do
       end do
    end do

    if (comm_i_am_master()) then

       DIAG(1:NM) = DIAG(1:NM) + DIAG_AUX(:)

       if (NS /= 0) then
          ALLOCATE(DA(NS),STAT=status)
          ASSERT(status==0)
          do i_proc = 2, comm_get_n_processors()
             call comm_save_recv(i_proc,msgtag_tddft_diagonl)
             call upck(DA)
             NS_UP = NM + (i_proc-2)*NS + 1
             NF_UP = NM + (i_proc-1)*NS
             DIAG(NS_UP:NF_UP) = DIAG(NS_UP:NF_UP) + DA(:)
          end do
          DEALLOCATE(DA,STAT=status)
          ASSERT(status==0)
       end if
    else
       if (NS /= 0) then
          call comm_init_send(comm_master_host,msgtag_tddft_diagonl)
          call pck(DIAG_AUX)
          call comm_send()
       end if
    end if

  end subroutine diag_assembl
  !*************************************************************

  !*************************************************************
  subroutine Q_calc(ir,Qm,Qdm1)
    !  Purpose:
    !
    !  Output: Qm   =  Nsp * gl_Q^(-1)xgl_Qx(gl_Q^(-1))^T = Nsp * gl_Q^(-1)
    !          Qdm1 =  gl_Q^(-1)
    !------------ Modules used ----------------------------------
    USE global_module, ONLY: gl_N_spin
    use read_module,   ONLY: read_Q
    USE linalg_module, ONLY: invert_sym_matrix!!, matmatmul
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(IN   ) :: ir
    real   (r8_kind), intent(OUT)   :: Qm(:,:), Qdm1(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
!!    integer(i4_kind) :: NSP
!!    real   (r8_kind) :: AP(size(Qm,1),size(Qm,1))
    !------------ Executable code ------------------------------------

    CALL read_Q(ir,Qm) ! read from file Qm:=gl_Q

    CALL invert_sym_matrix(Qm) ! Qdm1:=Qm^-1
    Qdm1 = Qm
    Qm = gl_N_spin * Qm

#if 0
    ! AP = Qdm1xQm

#if 0
    call matmatmul(Qdm1,Qm,AP)
    call matmatmul(AP,Qdm1,Qm,'N','T')
#endif

    !!FIXME: Workaround baceause of strange bug in DGEMM
    !!FIXME: QM = NSP * Qdm1 should be enough
    AP = matmul(Qdm1,Qm)
    Qm = matmul(AP,transpose(Qdm1))
#endif

  end subroutine Q_calc
  !***********************************************************

  !*************************************************************
  subroutine XC_calc(ir,N_SPIN,XC,Qdm1)
    !  Purpose: ..
    !------------ Modules used ----------------------------------
    use read_module,   ONLY: read_R
!!    use linalg_module, ONLY: matmatmul
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: ir,N_SPIN
    TYPE(arrmat2),         intent(INOUT) :: XC(:)
    real(kind=r8_kind),    intent(IN   ) :: Qdm1(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind) :: i
    real(kind=r8_kind)    :: AP(size(Qdm1,1),size(Qdm1,1))
    !------------ Executable code ------------------------------------

    call read_R(ir,XC)

    do i = 1,2*N_SPIN
       XC(i)%m = N_SPIN * XC(i)%m

       !! Workaround because of strange bug
       AP      = matmul(Qdm1,XC(i)%m)
       XC(i)%m = matmul(AP,transpose(Qdm1))

#if 0
       ! Q^-1xXC
       call matmatmul(Qdm1,XC(i)%m,AP)
       ! (Q^-1xXC)xQ^(-1)
       call matmatmul(AP,Qdm1,XC(i)%m,'N','T')
#endif
    end do


  end subroutine XC_calc
  !***********************************************************

  !*************************************************************
  subroutine missing_irreps
    !  Purpose: ..
    !------------ Modules used ----------------------------------
    use clebsch_gordan,      only: cgr => cg_reordered
    use global_module,       only: nirr => gl_N_irr
    use iounitadmin_module,  only: output_unit
    use symmetry_data_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables ---------------------
    integer(kind=i4_kind)                :: irc, ira, irb, mlt
    !------------ Executable code ------------------------------------

    do irc = nirr+1, size(cgr,1)
       do ira = 1, nirr
          do irb = 1, nirr
             mlt = cgr(irc,ira,irb)%mult
             if (mlt .ge. 1) then
                write(output_unit,*) "WARNING: ********************************** "
                write(output_unit,*) "WARNING: transitions from ",&
                     trim(ssym%name(ira))," -> ", trim(ssym%name(irb))," within ",irc,&
                     " has been omitted."
                write(output_unit,*) &
                     "WARNING: those transitions are equal to the difference "
                write(output_unit,*) &
                     "WARNING: between target unoccupied orbitals and source occupied orbitals"
             end if
          end do
       end do
    end do
  end subroutine missing_irreps
  !***********************************************************


!!$  !*************************************************************
!!$  subroutine tddft_
!!$  !  Purpose: ..
!!$  !------------ Modules used ----------------------------------
!!$    use
!!$    implicit none
!!$    !------------ Declaration of formal parameters ---------------
!!$    integer(kind=i4_kind), intent(     ) ::
!!$    real(kind=r8_kind),    intent(     ) ::
!!$    logical,               intent(     ) ::
!!$    character,             intent(     ) ::
!!$    !** End of interface *****************************************
!!$    !------------ Declaration of local variables -----------------
!!$    integer(kind=i4_kind)                ::
!!$    real(kind=r8_kind)                   ::
!!$    logical                              ::
!!$    character                            ::
!!$    !------------ Executable code --------------------------------
!!$
!!$
!!$  end subroutine tddft_
!!$  !***********************************************************



  !--------------- End of module -------------------------------------
END MODULE tddft_diag

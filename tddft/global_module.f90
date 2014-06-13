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
!===============================================================
! Public interface of module
!===============================================================
MODULE  global_module
  !---------------------------------------------------------------
  !
  !  Purpose:
  !  Database to hold GENERAL information needed in various parts
  !  of the program:
  !   + input parameters describing what to calculate and how to
  !     calculate it
  !   + all LARGE arrays needed in the program
  !     --> these will be handled as pointers
  !   +
  !   +
  !
  !
  !  Module called by: nearly every other module/subroutine
  !
  !  Author: HH
  !  Date:   11/98
  !
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification
  ! Author: SB
  ! Date:   11/04
  ! Description: totally rebuilded
  !
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
#include <def.h>
  USE type_module ! type specification parameters
  USE iounitadmin_module
  USE output_module
  USE comm_module
  USE msgtag_module
  USE datatype

  IMPLICIT NONE
  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Declaration of constants and variables ------------

  LOGICAL, PUBLIC               :: &
       gl_calcall, gl_oscstr, &
       gl_XC, gl_SS, gl_ST, &
       gl_S_App, gl_S_APP_XC, gl_noRI, &
       gl_lanczos, &
       gl_NTO                     !!!!!!MH NTO CALCULATIONS

  INTEGER(KIND=i4_kind), PUBLIC :: gl_N_spin, gl_N_irr  !! spin & irrep

  INTEGER(KIND=i4_kind), PUBLIC :: gl_ir, gl_X, gl_C

  INTEGER(KIND=i4_kind), PUBLIC, ALLOCATABLE :: gl_N_k(:), &
       gl_N_as(:,:), &
       gl_N_as_slave(:,:), &
       gl_N_as_mastr(:,:)

  type(arrmat1), ALLOCATABLE, PUBLIC :: gl_eps(:,:), gl_eta(:,:)
  type(intmat2), ALLOCATABLE, PUBLIC :: gl_MO (:,:), gl_IRR(:,:)

  type(intmat1), ALLOCATABLE, PUBLIC :: gl_what_processor(:,:), &
       &                                gl_what_N_as     (:,:)

  REAL(KIND=r8_kind), PUBLIC, ALLOCATABLE :: gl_vec(:,:), gl_val(:)
  INTEGER(KIND=i4_kind), PUBLIC :: gl_Nlow, gl_max_sp_trans
  INTEGER(KIND=i4_kind), PUBLIC :: gl_Iter
  INTEGER(KIND=i4_kind), PUBLIC :: gl_MaxIt

  REAL(KIND=r8_kind), PUBLIC              :: gl_eig_crite, gl_LOWESTeV

  TYPE(ARRMAT2),      PUBLIC, ALLOCATABLE :: gl_L(:)
  TYPE(ARRMAT2),      PUBLIC, ALLOCATABLE :: gl_T(:)
  REAL(KIND=r8_kind), PUBLIC, ALLOCATABLE :: gl_Q(:,:)


  !------------ public functions and subroutines ------------------
  PUBLIC global_alloc,global_alloc_M,global_dealloc_M,global_dealloc


  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of constants and variables ----

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
CONTAINS

  !*************************************************************
  SUBROUTINE global_alloc(what,n_irrep,dim1,dim2)
    ! Purpose: allocate 'what'
    !------------ Modules used ------------------- ---------------

    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------

    INTEGER(KIND=i4_kind), PARAMETER :: ZEROi = 0_i4_kind
    REAL(KIND=r8_kind), PARAMETER    :: ZEROr = 0.0_r8_kind

    CHARACTER(3), intent(in) :: what
    INTEGER(KIND=i4_kind), intent(in), optional  :: dim1(:,:)
    INTEGER(KIND=i4_kind), intent(in), optional  :: dim2(:)
    INTEGER(KIND=i4_kind), intent(in)  :: n_irrep
    !------------ Declaration of local variables --------------------
    INTEGER(KIND=i4_kind)  :: status, n_spin, n_procs
    INTEGER(KIND=i4_kind)  :: i_sp, i_ir, dm
    INTEGER(KIND=i4_kind)  :: NS, NF, i, nas, ias, my_proc
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code -----------------------------------

    n_spin  = gl_N_spin
    n_procs = comm_get_n_processors()

    SELECT CASE (what)

    case ('ALL')
       ALLOCATE(gl_eps(n_irrep, n_spin), &
            &   gl_eta(n_irrep, n_spin), &
            &   gl_MO (n_irrep, n_spin), &
            &   gl_IRR(n_irrep, n_spin), STAT = status)
       ASSERT(status==0)
       do i_sp = 1, n_spin
          do i_ir = 1, n_irrep
             dm = dim1(i_ir,i_sp)
             ALLOCATE(gl_eps(i_ir,i_sp)%m (dm),  STAT=status)
             ASSERT(status==0)
             ALLOCATE(gl_eta(i_ir,i_sp)%m (dm),  STAT=status)
             ASSERT(status==0)
             ALLOCATE(gl_MO (i_ir,i_sp)%m(dm,2),STAT=status)
             ASSERT(status==0)
             ALLOCATE(gl_IRR(i_ir,i_sp)%m(dm,2),STAT=status)
             ASSERT(status==0)
          end do
       end do

    CASE ('chr')
       ALLOCATE(gl_N_k(n_irrep),STAT=status)
       ASSERT(status==0)
       gl_N_k = 0_i4_kind

    CASE ('nas')
       ALLOCATE(gl_N_as      (n_irrep,2),STAT=status)
       ASSERT(status==0)
       gl_N_as = 0_i4_kind
       ALLOCATE(gl_N_as_slave(n_irrep,n_spin),STAT=status)
       ASSERT(status==0)
       ALLOCATE(gl_N_as_mastr(n_irrep,n_spin),STAT=status)
       ASSERT(status==0)

    CASE ('prc')
       ALLOCATE(gl_what_processor(n_irrep,n_spin),&
            &   gl_what_N_as     (n_irrep,n_spin),STAT=status)
       ASSERT(status==0)
       do i_sp = 1, n_spin
          do i_ir = 1, n_irrep
             dm = dim1(i_ir,i_sp)
             ALLOCATE(gl_what_processor(i_ir,i_sp)%m(dm),&
                  &   gl_what_N_as     (i_ir,i_sp)%m(dm), STAT=status)
             ASSERT(status==0)
             gl_what_processor(i_ir,i_sp)%m = -1_i4_kind
             gl_what_N_as     (i_ir,i_sp)%m = -1_i4_kind

             NF = 0
             do i = 1, n_procs
                NS = NF + 1
                NF = gl_N_as_mastr(i_ir,i_sp) &
                     + (i-1)*gl_N_as_slave(i_ir,i_sp)
                gl_what_processor(i_ir,i_sp)%m(NS:NF) = i
             end do

             my_proc = comm_myindex()
             nas   = 0
             do ias = 1, dm
                if (my_proc /= gl_what_processor(i_ir,i_sp)%m(ias)) then
                   cycle
                else
                   nas = nas + 1
                end if
                gl_what_N_as(i_ir,i_sp)%m(ias) = nas
             end do

#if 0
             print *,"MyID = ",my_proc," i_ir = ",i_ir," i_sp = ",i_sp
             print *,"MyID = ",my_proc," SPLITTING gl_what_processor : ", &
                  gl_what_processor(i_ir,i_sp)%m(:)
             print *,"MyID = ",my_proc," SPLITTING gl_what_N_as      : ", &
                  gl_what_N_as(i_ir,i_sp)%m(:)
#endif
          end do
       end do

    END SELECT

  END SUBROUTINE global_alloc
  !*************************************************************

  !*************************************************************
  SUBROUTINE global_alloc_M(what,Ndim,Ndim2,NK)
    ! Purpose: allocate 'what'
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind), PARAMETER :: ZEROi = 0_i4_kind
    REAL(KIND=r8_kind),    PARAMETER :: ZEROr = 0.0_r8_kind

    CHARACTER(3),          intent(in)             :: what
    INTEGER(KIND=i4_kind), intent(in)             :: Ndim
    INTEGER(KIND=i4_kind), intent(inout),optional :: Ndim2,NK
    !------------ Declaration of local vaiables --------------------
    INTEGER(KIND=i4_kind)  :: status, n_spin, i
    INTEGER(KIND=i4_kind)  :: DM1, DM2
    !------------ Declaration of subroutines used ----------------
    EXTERNAL error_handler
    !------------ Executable code -----------------------------------

    n_spin = gl_N_spin

    SELECT CASE (what)

    CASE ('glT')
       if (Ndim2 == 0) Ndim2 = Ndim
       ALLOCATE(gl_T(n_spin * 2), STAT=status)
       ASSERT(status==0)
       do i = 1, n_spin * 2
          select case(i)
          case (1)
             DM1 = Ndim
             DM2 = Ndim
          case (2)
             DM1 = Ndim
             DM2 = Ndim2
          case (3)
             DM1 = Ndim2
             DM2 = Ndim2
          case (4)
             DM1 = Ndim2
             DM2 = Ndim
          end select
          ALLOCATE(gl_T(i)%m(DM1,DM2),STAT =status)
          ASSERT(status==0)
          gl_T(i)%m = ZEROr
       end do

     CASE ('glL')
        ALLOCATE(gl_L(n_spin), STAT=status)
        ASSERT(status==0)
        do i = 1, n_spin
           select case (i)
           case (1)
              DM1 = NDIM
           case (2)
              DM1 = NDIM2
           end select
           ALLOCATE(gl_L(i)%m(DM1,NK),STAT =status)
           ASSERT(status==0)
           gl_L(i)%m = ZEROr
        end do

    CASE ('Qlm')
       ALLOCATE(&
            gl_Q(Ndim,Ndim),&
            STAT = status)

    CASE ('vec')
       ALLOCATE(gl_vec(Ndim,Ndim2),STAT = status)
       ASSERT(status==0)
       gl_vec = ZEROr

    CASE ('val')
       ALLOCATE(gl_val(Ndim),STAT = status)
       ASSERT(status==0)
       gl_val = ZEROr

    END SELECT
  END SUBROUTINE global_alloc_M
  !*************************************************************

  !*************************************************************
  SUBROUTINE global_dealloc_M
    USE debug, only:show
    ! Purpose: allocate 'what'
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of local vaiables --------------------
    INTEGER(KIND=i4_kind)  :: status, n_spin, i
    !------------ Declaration of subroutines used ----------------

    n_spin = gl_N_spin

#if 0
    if(comm_i_am_master() .and. comm_parallel()) then
       call comm_init_send(comm_master_host,msgtag_tddft_dealloM)
       call comm_send()
    end if
#endif

    if (allocated(gl_L)) then
       do i = 1,n_spin
          DEALLOCATE(gl_L(i)%m,STAT=status)
          ASSERT(status==0)
       end do
       DEALLOCATE(gl_L,STAT=status)
       ASSERT(status==0)
    end if

    if (allocated(gl_vec)) then
       deallocate(gl_vec,STAT=status)
       ASSERT(status==0)
    end if

    if (allocated(gl_val)) then
       deallocate(gl_val,STAT=status)
       ASSERT(status==0)
    end if

    if (allocated(gl_Q)) then
       deallocate(gl_Q,STAT=status)
       ASSERT(status==0)
    end if

    if (allocated(gl_T)) then
       do i = 1, n_spin * 2
          deallocate(gl_T(i)%m,STAT=status)
          ASSERT(status==0)
       end do
       deallocate(gl_T,STAT=status)
       ASSERT(status==0)
    end if

  END SUBROUTINE global_dealloc_M
  !*************************************************************

  SUBROUTINE global_dealloc()
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of local vaiables --------------------
    INTEGER(KIND=i4_kind)  :: status, i_ir, i_sp
    !------------ Declaration of subroutines used ----------------


    do i_ir = 1, gl_N_irr
       do i_sp = 1, gl_N_spin
          DEALLOCATE(gl_what_processor(i_ir,i_sp)%m,STAT=status)
          ASSERT(status==0)
       end do
    end do
    DEALLOCATE(gl_what_processor,STAT=status)
    ASSERT(status==0)

    do i_ir = 1, gl_N_irr
       do i_sp = 1, gl_N_spin
          DEALLOCATE(gl_what_N_as(i_ir,i_sp)%m,STAT=status)
          ASSERT(status==0)
       end do
    end do
    DEALLOCATE(gl_what_N_as,STAT=status)
    ASSERT(status==0)

    do i_ir = 1, gl_N_irr
       do i_sp = 1, gl_N_spin
          DEALLOCATE(gl_eps(i_ir,i_sp)%m,STAT=status)
          ASSERT(status==0)
       end do
    end do
    DEALLOCATE(gl_eps, STAT = status)
    ASSERT(status==0)

    do i_ir = 1, gl_N_irr
       do i_sp = 1, gl_N_spin
          DEALLOCATE(gl_eta(i_ir,i_sp)%m,STAT=status)
          ASSERT(status==0)
       end do
    end do
    DEALLOCATE(gl_eta, STAT = status)
    ASSERT(status==0)

    do i_ir = 1, gl_N_irr
       do i_sp = 1, gl_N_spin
          DEALLOCATE(gl_IRR(i_ir,i_sp)%m,STAT=status)
          ASSERT(status==0)
       end do
    end do
    DEALLOCATE(gl_IRR, STAT = status)
    ASSERT(status==0)

    if (allocated(gl_N_k)) then
       DEALLOCATE(gl_N_k,STAT=status)
       ASSERT(status==0)
    end if

    if (allocated(gl_N_as)) then
       DEALLOCATE(gl_N_as,STAT=status)
       ASSERT(status==0)
    end if

    do i_ir = 1, gl_N_irr
       do i_sp = 1, gl_N_spin
          DEALLOCATE(gl_MO(i_ir,i_sp)%m,STAT=status)
          ASSERT(status==0)
       end do
    end do
    DEALLOCATE(gl_MO,STAT=status)
    ASSERT(status==0)

  END SUBROUTINE global_dealloc



!!$  !*************************************************************
!!$  INTEGER(KIND=i4_kind) FUNCTION global_data_get_S_index()
!!$    !  Purpose:
!!$    !  Gives back the array value
!!$    !     global_data_S_index(combined_norm_index)
!!$    !------------ Modules used ------------------- ---------------
!!$    IMPLICIT NONE
!!$    !------------ Declaration of formal parameters ---------------
!!$    INTEGER(KIND=i4_kind),INTENT(IN) :: combined_norm_index
!!$    !------------ Declaration of local variables -----------------
!!$    !------------ Declaration of subroutines used ----------------
!!$    EXTERNAL error_handler
!!$    !------------ Executable code --------------------------------
!!$
!!$
!!$  END FUNCTION global_data_get_S_index
!!$  !*************************************************************

  !--------------- End of module ----------------------------------
END MODULE global_module

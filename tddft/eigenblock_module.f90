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
!===============================================================
! Public interface of module
!===============================================================
MODULE  eigenblock_module
  !---------------------------------------------------------------
  !
  ! Author: SB
  ! Date: 07/2005
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
# include "def.h"
  USE type_module ! type specification parameters
  USE iounitadmin_module 
  USE output_module
  USE comm_module
  USE debug

  IMPLICIT NONE
  SAVE            ! save all variables defined in this module
  PRIVATE         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  PUBLIC eigenblock_mult


  !================================================================
  ! End of public interface of module
  !================================================================

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
CONTAINS

  SUBROUTINE eigenblock_mult(N, K, A, D)
    !  Purpose: 
    !  Computes the product of matrix M=[eps^2+L*T*L^T] with a 
    !  block of vectors  A(N,K), i.e.
    !------------ Modules used ------------------- ---------------
    USE global_module
    USE noRI_module, only: noRI_4C
    USE constants, only: ZERO, ONE
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind),INTENT(in   )  :: N,K
    REAL   (KIND=r8_kind),INTENT(in   )  :: A(:)
    REAL   (KIND=r8_kind),INTENT(inout)  :: D(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    REAL   (KIND=r8_kind), ALLOCATABLE :: A_2d_up(:,:), A_2d_dn(:,:)
    REAL   (KIND=r8_kind), ALLOCATABLE :: aD(:,:)
    REAL   (KIND=r8_kind), ALLOCATABLE :: D_2d(:,:), KERN(:,:)
    INTEGER(KIND=i4_kind)              :: status, ias, ias_dn, counter, j, NK
    INTEGER(KIND=i4_kind)              :: CQ,CX, n_spin, DIM

    INTEGER(KIND=i4_kind)              :: i_ir, NSP(gl_N_spin)

    INTEGER(KIND=i4_kind),PARAMETER    :: &
         UP = 1, &
         DN = 2

    INTEGER(KIND=i4_kind),PARAMETER    :: &
         UPUP = 1, &
         UPDN = 2, &
         DNDN = 3, &
         DNUP = 4
    
    !------------ Executable code -----------------------------------

    i_ir   = gl_ir
    n_spin = gl_N_spin

    do j = 1, n_spin
       NSP(j) = gl_N_as(i_ir,j)
    end do

    if(comm_i_am_master()) then
       ! first increment running counter for current number of iterations
       ! and display in output/traceoutput for showing progress
       CALL eigenblock_showprogress(gl_Iter)
    end if

    if (comm_i_am_master()) then
       ALLOCATE(A_2d_up(gl_N_as_mastr(i_ir,UP),K),stat=status)
       ASSERT(status==0)
    else
       ALLOCATE(A_2d_up(gl_N_as_slave(i_ir,UP),K),stat=status)
       ASSERT(status==0)
    end if
    A_2d_up    = ZERO

    if (gl_N_spin == 2) then
       if (comm_i_am_master()) then
          ALLOCATE(A_2d_dn(gl_N_as_mastr(i_ir,DN),K),stat=status)
          ASSERT(status==0)
       else
          ALLOCATE(A_2d_dn(gl_N_as_slave(i_ir,DN),K),stat=status)
          ASSERT(status==0)
       end if
       A_2d_dn = ZERO
    end if

    counter = 1_i4_kind   ! map 1d array C on 2d matrix A_2d

    DO j=1,K
       do ias = 1, N
          if (ias <= NSP(UP)) then
             if (gl_what_N_as(i_ir,UP)%m(ias) == -1 ) cycle
             A_2d_up(gl_what_N_as(i_ir,UP)%m(ias),j) = A(counter+ias-1)
          else             
             ias_dn =ias - NSP(UP)
             if (gl_what_N_as(i_ir,DN)%m(ias_dn) == -1 ) cycle
             A_2d_dn(gl_what_N_as(i_ir,DN)%m(ias_dn),j) = A(counter+ias-1)
          end if
       end do
       counter = counter + N
    END DO

    if (gl_N_spin == 1) then
       DIM = NSP(UP)
    else
       DIM = NSP(UP)+NSP(DN)
    end if

    ALLOCATE(D_2d(DIM,K),STAT=status)
    ASSERT(status==0)
    D_2d = ZERO

    NK = size(gl_Q,1)
    ALLOCATE(KERN(NK,NK),STAT = status)
    ASSERT(status==0)

    if (gl_S_App_XC) then
       CQ = 0 
       CX = 1
    elseif(gl_S_App) then
       CQ = 1
       CX = 0
    else
       CQ = 1
       CX = 1 
    endif

    if (gl_noRI) then
       CQ = 1
       CX = 0
    end if

    if (gl_N_spin == 1) then !! CS

       if (gl_SS) then
          KERN = CQ * gl_Q + CX * (gl_T(UPUP)%m + gl_T(UPDN)%m) !!SS
       else
          KERN =             CX * (gl_T(UPUP)%m - gl_T(UPDN)%m) !!ST
       end if

       call eigenblock_paramult(gl_L(UP)%m,&
            &              KERN,gl_L(UP)%m,&
            &              A_2d_up,D_2d,NSP(UP),&
            &              gl_N_as_mastr(i_ir,UP),gl_N_as_slave(i_ir,UP),&
            &              gl_eps(i_ir,UP)%m)
    else !! OS

       ALLOCATE(aD(NSP(UP),K),STAT = status)
       ASSERT(status==0)
       aD = ZERO

       !! *** upup ***
       KERN = CQ * gl_Q + CX * gl_T(UPUP)%m
        call eigenblock_paramult(gl_L(UP)%m,&
             &              KERN,gl_L(UP)%m,&
             &              A_2d_up,aD,NSP(UP),&
             &              gl_N_as_mastr(i_ir,UP),gl_N_as_slave(i_ir,UP),&
             &              gl_eps(i_ir,UP)%m)
       D_2d(1:NSP(UP),:) = aD
       !! *** updn ***
       KERN = CQ * gl_Q + CX * gl_T(UPDN)%m
       call eigenblock_paramult(gl_L(UP)%m,&
            &              KERN,gl_L(DN)%m,&
            &              A_2d_dn,aD,NSP(UP),&
            &              gl_N_as_mastr(i_ir,UP),gl_N_as_slave(i_ir,UP))
       D_2d(1:NSP(UP),:) = D_2d(1:NSP(UP),:) + aD
       DEALLOCATE(aD,STAT=status)
       ASSERT(status==0)

       ALLOCATE(aD(NSP(DN),K),STAT = status)
       ASSERT(status==0)
       aD = ZERO
       !! *** dnup ***
       KERN = CQ * gl_Q + CX * gl_T(DNUP)%m
       call eigenblock_paramult(gl_L(DN)%m,&
            &              KERN,gl_L(UP)%m,&
            &              A_2d_up,aD,NSP(DN),&
            &              gl_N_as_mastr(i_ir,DN),gl_N_as_slave(i_ir,DN))       
       D_2d(NSP(UP)+1:NSP(UP)+NSP(DN),:) = aD

       !! *** dndn ***
       KERN = CQ * gl_Q + CX * gl_T(DNDN)%m
       call eigenblock_paramult(gl_L(DN)%m,&
            &              KERN,gl_L(DN)%m,&
            &              A_2d_dn,aD,NSP(DN),&
            &              gl_N_as_mastr(i_ir,DN),gl_N_as_slave(i_ir,DN),&
            &              gl_eps(i_ir,DN)%m)
       D_2d(NSP(UP)+1:NSP(UP)+NSP(DN),:) = &
            D_2d(NSP(UP)+1:NSP(UP)+NSP(DN),:) &
            + aD
       DEALLOCATE(aD,STAT=status)
       ASSERT(status==0)

       DEALLOCATE(A_2d_dn,STAT=status)
       ASSERT(status==0)

    end if
    DEALLOCATE(A_2d_up,STAT=status)
    ASSERT(status==0)

!!$    call show('D_IN',D_2d)

    counter = 1
    DO j=1,K
       D(counter:counter+N-1) = D_2d(:,j)
       counter = counter + N
    END DO

!!    call show('EIGENBLOCK: D',D)
    
    if (gl_noRI .and. (NK .ne. 0)) call noRI_4C(gl_X,gl_C,i_ir,N,NSP,K,gl_eps,gl_eta,A,D)

    DEALLOCATE(KERN,STAT=status)
    ASSERT(status==0)

    DEALLOCATE(D_2d,STAT=status)
    ASSERT(status==0)

  contains

    subroutine eigenblock_paramult(leftMatr,cntrMatr,rghtMatr,&
         multMatr,rsltMatr, dim, dim_ma_spl, dim_sl_spl, eps)
      !*************************************************************
      !PURPOSE: parallel multiplication eps^2xA + [lMxcMxrM]xA results is on master
      ! step1: 
      !     a. auxM1 = rM x A = [K,Nas_spl]x[Nas_spl,K_prime]
      !     b. collect all [K,K_prime] on the master and summed up
      ! step2 (on the master):
      !     a. auxM2 = cM x auxM1 = [K,K_prime]
      !     b. send the result to all slaves
      ! step3:
      !     a. auxM3 = lM x auxM2 = [Nas_spl,K] x [K,K_prime] 
      ! step4:
      !     a. if eps is present then
      !     b. auxM3 = auxM3 + eps^2 x A = [Nas_spl,K_prime]
      !     c. collect auxM3 on the master and combine rsltMatr
      !     d. send rsltMatr to all slaves
      !*************************************************************
      USE comm, only: comm_reduce, comm_bcast, comm_barrier
      USE xpack
      USE msgtag_module
      USE linalg_module, only: matmatmul
      IMPLICIT NONE
      !------------ Declaration of formal parameters ---------------
      REAL   (KIND=r8_kind),INTENT(in   )            :: leftMatr(:,:)
      REAL   (KIND=r8_kind),INTENT(in   )            :: cntrMatr(:,:)
      REAL   (KIND=r8_kind),INTENT(in   )            :: rghtMatr(:,:)
      REAL   (KIND=r8_kind),INTENT(in   )            :: multMatr(:,:)
      REAL   (KIND=r8_kind),INTENT(inout)            :: rsltMatr(:,:)
      INTEGER(KIND=i4_kind),INTENT(IN   )            :: dim, &
           dim_ma_spl, dim_sl_spl
      REAL   (KIND=r8_kind),INTENT(in   ), OPTIONAL  :: eps(:)
      !** End of interface *****************************************
      !------------ Declaration of local variables -----------------
      REAL   (KIND=r8_kind), PARAMETER   :: ZERO = 0.0_r8_kind
      REAL   (KIND=r8_kind), ALLOCATABLE :: auxM1(:,:), auxM2(:,:), auxM3(:,:)
      INTEGER(KIND=i4_kind)              :: K1, K2, status, i_proc, N1, ias, &
           Ns, Nf
      !------------ Executable code --------------------------------

     ! print *,"eigenblock_paramult: barrier in  1"
      call comm_barrier()
     ! print *,"eigenblock_paramult: barrier out 1"

      K1 = size(leftMatr,2)
      K2 = size(rsltMatr,2)

      allocate(auxM1(K1,K2),STAT=status)
      ASSERT(status==0)
      auxM1 = ZERO
      allocate(auxM2(K1,K2),STAT=status)
      ASSERT(status==0)
      auxM2 = ZERO

      call matmatmul(rghtMatr, multMatr, auxM1, "T", "N")

      if (dim_sl_spl /=0 ) then
         call comm_reduce(auxM1)
         call comm_bcast(auxM1)
      end if

      call matmatmul(cntrMatr, auxM1, auxM2)
      DEALLOCATE(auxM1,STAT=status)
      ASSERT(status==0)

      N1 = size(leftMatr,1)
      allocate(auxM3(N1,K2),STAT=status)
      auxM3 = ZERO
      call matmatmul(leftMatr, auxM2, auxM3)
      DEALLOCATE(auxM2,STAT=status)
      ASSERT(status==0)


      if (present(eps)) then
         do ias = 1, N1
            auxM3(ias,:) = auxM3(ias,:) + multMatr(ias,:)*&
                 (eps(ias)*eps(ias))
         end do
      end if

      if (comm_i_am_master()) rsltMatr(1:dim_ma_spl,:) = auxM3 

      if (dim_sl_spl /=0 ) then
         if (comm_i_am_master()) then
            ALLOCATE(auxM1(dim_sl_spl,K2),STAT=status)
            ASSERT(status==0)
            do i_proc = 2, comm_get_n_processors()
               call comm_save_recv(i_proc,msgtag_tddft_parmlt3)
               call upck(auxM1)
               Ns = dim_ma_spl + (i_proc - 2) * dim_sl_spl + 1
               Nf = dim_ma_spl + (i_proc - 1) * dim_sl_spl
               rsltMatr(Ns:Nf,:) = auxM1
            end do
            DEALLOCATE(auxM1,STAT=status)
            ASSERT(status==0)
         else
            call comm_init_send(comm_master_host,msgtag_tddft_parmlt3)
            call pck(auxM3)
            call comm_send()  
         end if

         call comm_bcast(rsltMatr)
      end if

      DEALLOCATE(auxM3,STAT=status)
      ASSERT(status==0)

    end subroutine eigenblock_paramult

  END SUBROUTINE eigenblock_mult

  !*************************************************************
  SUBROUTINE eigenblock_showprogress(iteration_cycle)
    !  Purpose: 
    !   - increment number of iterations
    !   - write to trace file for runtime information
    !   - write to output file
    !------------ Modules used ------------------- ---------------
    IMPLICIT NONE
    !------------ Declaration of formal parameters ---------------
    INTEGER(KIND=i4_kind),INTENT(INOUT) :: iteration_cycle
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !------------ Executable code -----------------------------------

    ! increment number of iterations
    iteration_cycle = iteration_cycle + 1_i4_kind

    ! write to trace file for runtime information
    CALL write_to_trace_unit("Iteration ",iteration_cycle)

    ! write to output file
    CALL write_to_output_units("Iteration ",iteration_cycle)

  END SUBROUTINE eigenblock_showprogress
  !*************************************************************

  !--------------- End of module ----------------------------------
END MODULE eigenblock_module

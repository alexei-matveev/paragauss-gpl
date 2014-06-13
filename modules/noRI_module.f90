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
module noRI_module
  !---------------------------------------------------------------
  !
  !  Purpose: ...
  !
  !
  !  Module called by: ...
  !
  !
  !  References: ...
  ! 
  !
  !  Author: ...
  !  Date: ...
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
! define FPP_TIMERS 2
#include "def.h"
  use iounitadmin_module, only: write_to_output_units,output_unit
  use type_module ! type specification parameters
  use datatype
  use orbitalstore_module, only: orbital_type, orbital_gradient_type
  use machineparameters_module, only: vec_length => machineparameters_veclen
  use symmetry_data_module
  use linalg_module
  use clebsch_gordan,      only: cg=>cg_eliminated, prod_bas
  use grid_module, only: more_grid_atom, atomicweight
  use grid_module, only: grid_loop_setup_atom
  use debug
  use msgtag_module
  use comm_module
  use ch_response_module
  use orbital_module
  use resp_util_module
  use error_module, only: MyID

  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ public functions and subroutines ------------------
  public noRI_2c, noRI_4c

  !================================================================
  ! End of public interface of module
  !================================================================
  integer(i4_kind), parameter :: UP = 1, DN = 2
  integer(i4_kind), parameter :: C2 = 1, C4 = 2 !!TWO CENTER / FOUR CENTER
  integer(i4_kind), parameter :: & !! VERY IMPORTANT
       & UPUP = 1, &
       & UPDN = 2, &
       & DNDN = 4, &
       & DNUP = 3

  integer(i4_kind), parameter :: & !! VERY IMPORTANT
       & X_NONE   = 0, &
       & C_NONE   = 0, &
       & C_VWN    = 1, &
       & C_PWlda  = 2, &
       & C_PBE    = 3, &
       & C_PW91   = 4, &
       & C_PERDEW = 5

  integer(i4_kind), parameter ::&
       & A = 1, &
       & B = 2

  integer(i4_kind), parameter ::&
       & AA = 1, &
       & AB = 2, &
       & BB = 4, &
       & BA = 3

  integer(i4_kind), parameter ::&
       & XCAA = 1, &
       & XCBB = 2, &
       & XCAB = 3
  !!         & BA = 4

  integer(i4_kind), parameter ::&
       & AAA = 1, &
       & BBB = 2, &
       & BAA = 3, &
       & ABB = 4, &
       & AAB = 5, &
       & BAB = 6

  integer(i4_kind), parameter ::&
       & AAAA = 1, &
       & BBBB = 2, &
       & AABB = 3, &
       & AAAB = 4, &
       & BBAB = 5, &
       & ABAB = 6

  !----------------------------------------------------------------
  type(orbital_type),pointer             :: fcts_orbs_ch(:)
  type(orbital_gradient_type),pointer    :: grads_ch(:)

  type(orbital_type),pointer             :: orbs_ob(:) 
  type(orbital_gradient_type),pointer    :: orbs_grads(:)

  type(arrmat4),     allocatable, target :: phi (:)
  type(arrmat5),     allocatable, target :: gphi(:)
  !------------ Subroutines ---------------------------------------

  FPP_TIMER_DECL(XC_all)
  FPP_TIMER_DECL(XC_input)
  FPP_TIMER_DECL(XC_pb_sym)
  FPP_TIMER_DECL(XC_pb_sym_OC)
  FPP_TIMER_DECL(XC_pb_sym_PC)
  FPP_TIMER_DECL(XC_pb_sym_dRC)
  FPP_TIMER_DECL(XC_functional)
  FPP_TIMER_DECL(XC_output)
  FPP_TIMER_DECL(XC_int)
  FPP_TIMER_DECL(XC_int_dV)
  FPP_TIMER_DECL(XC_int_NB)
  FPP_TIMER_DECL(XC_int_B)

contains

  !*************************************************************
  subroutine noRI_2c(X,C)
    use msgtag_module
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use constants, only: zero
    use comm, only: comm_bcast
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind), intent(IN),optional :: X,C
    !** End of interface *****************************************
    integer(i4_kind)             :: n_spin, n_ir, i_ir, n_dim, i, X_KIND, C_KIND
    type(arrmat2), allocatable   :: RM(:)
    integer(i4_kind)             :: status
    real   (r8_kind)             :: tt
!!$
!!$    FPP_TIMER_ZERO(input)
!!$    FPP_TIMER_ZERO(OB_calc)
!!$    FPP_TIMER_ZERO(XC_calc)
!!$    FPP_TIMER_ZERO(XC_assm)
!!$    FPP_TIMER_ZERO(XC_assm_NOBLAS)
!!$    FPP_TIMER_ZERO(XC_assm_BLAS)
!!$    FPP_TIMER_ZERO(output)
!!$    FPP_TIMER_ZERO(all)

    if (comm_i_am_master()) then
       X_KIND = X
       C_KIND = C
       if(comm_parallel()) then  ! is this a parallel run ?
          call comm_init_send(comm_all_other_hosts,msgtag_nori_2c_send)
          call comm_send()
       end if
    end if

    call comm_bcast(X_KIND)
    call comm_bcast(C_KIND)

    FPP_TIMER_START(XC_all)
    FPP_TIMER_START(XC_input)

    n_ir   = symmetry_data_n_irreps()
    n_spin = symmetry_data_n_spin()

    call fit_fct_allocate_response(fcts_ch = fcts_orbs_ch, grads = grads_ch)

#if 0
    call orbital_setup_response(vec_length,.true.,.true.)
    call orbital_allocate (orbs_ob,orbs_grads)
#endif

    ALLOCATE(RM(2*n_spin),STAT = status)
    ASSERT(status==0)

    FPP_TIMER_START(XC_input)

    do i_ir = 1, n_ir

       n_dim = dimension_of_fit_ch(i_ir)
       if (n_dim==0) cycle

       do i = 1, 2*n_spin 
          ALLOCATE(RM(i)%m(n_dim,n_dim),STAT=status)
          ASSERT(status==0)
          RM(i)%m = ZERO
       end do
       FPP_TIMER_STOP(XC_input)

       call noRI_XC(X_KIND,C_KIND,C2,i_ir,RM)

       FPP_TIMER_START(XC_output)
       if (comm_i_am_master()) call noRI_totape(C2,i_ir,n_spin,RM)
       do i = 1, 2*n_spin
          DEALLOCATE(RM(i)%m,STAT = status)
          ASSERT(status==0)
       end do
    end do

    DEALLOCATE(RM,STAT=status)
    ASSERT(status==0)

#if 0
    call orbital_free(orbs_ob)
#endif

    call fit_fct_free_response(fcts = fcts_orbs_ch, grads = grads_ch)
    FPP_TIMER_STOP(XC_output)
    FPP_TIMER_STOP(XC_all)

    WRITE (*,*) MyID, "noRI TIMER      "
    tt = FPP_TIMER_VALUE(XC_all)
    WRITE (*,*) "   SUMMARY            ", tt
    tt = FPP_TIMER_VALUE(XC_input)
    WRITE (*,*) "   |- PREPARATION     ", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym)
    WRITE (*,*) "   |- XC  PROD. BAS.SYM. ", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym_OC)
    WRITE (*,*) "      |- XC  PROD. BAS.SYM. ORB.  CALC.", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym_PC)
    WRITE (*,*) "      |- XC  PROD. BAS.SYM. PHI.  CALC.", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym_dRC)
    WRITE (*,*) "      \- XC  PROD. BAS.SYM. dRHO. CALC.", tt
    tt = FPP_TIMER_VALUE(XC_functional)
    WRITE (*,*) "   |- XC  FUNCTIONAL  ", tt
    tt = FPP_TIMER_VALUE(XC_output)
    WRITE (*,*) "   |- OUTPUT          ", tt
    tt = FPP_TIMER_VALUE(XC_int)
    WRITE (*,*) "   \- XC  INTEGRATION ", tt
    tt = FPP_TIMER_VALUE(XC_int_dV)
    WRITE (*,*) "      |-  XC deltaV PART ", tt
    tt = FPP_TIMER_VALUE(XC_int_NB)
    WRITE (*,*) "      |-  XC NOBLAS PART ", tt
    tt = FPP_TIMER_VALUE(XC_int_B)
    WRITE (*,*) "      \-  XC   BLAS PART ", tt

  end subroutine noRI_2c
  !*************************************************************\

  !*************************************************************
  subroutine noRI_4c(X_KIND,C_KIND,i_ir,N,NSP,K,eps,eta,A,RES)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use constants,        only: zero, one, two
    use resp_util_module, only: resp_util_buildr
    use global_module,    only: gl_N_as_mastr, gl_N_as_slave, gl_SS, gl_ST, gl_N_spin
    use orbital_module
    use comm, only: comm_bcast

    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),intent( IN) :: i_ir, N, NSP(:), K, X_KIND,C_KIND
    real   (r8_kind),intent( IN) :: A(:)
    type   (arrmat1),intent( IN) :: eps(:,:), eta(:,:)
    real   (r8_kind),intent(OUT) :: RES(:)
    !** End of interface *****************************************
    integer(i4_kind)             :: n_spin
    integer(i4_kind)             :: status, i, j, ii
    integer(i4_kind)             :: counter
    real   (r8_kind)             :: CC, tt
    real   (r8_kind),ALLOCATABLE :: D(:,:)
    real   (r8_kind),ALLOCATABLE :: eps_full(:), eta_full(:)
    type   (arrmat2),ALLOCATABLE :: RM(:), DSP(:)
    type   (arrmat1),ALLOCATABLE :: SV(:)

    IF ((X_KIND == X_NONE) .AND. (C_KIND==C_NONE)) RETURN

!!$    FPP_TIMER_ZERO(input)
!!$    FPP_TIMER_ZERO(OB_calc)
!!$    FPP_TIMER_ZERO(XC_calc)
!!$    FPP_TIMER_ZERO(XC_assm)
!!$    FPP_TIMER_ZERO(XC_assm_NOBLAS)
!!$    FPP_TIMER_ZERO(XC_assm_BLAS)
!!$    FPP_TIMER_ZERO(output)
!!$    FPP_TIMER_ZERO(all)

    FPP_TIMER_START(XC_all)
    FPP_TIMER_START(XC_input)

    n_spin = symmetry_data_n_spin()

    call orbital_setup(vec_length,.true.,.true.)
    call orbital_allocate (orbs_ob,orbs_grads)!!,phi_ob=phi_ob,orb_ch=orbs_ch)

    if (n_spin == 1) then
       CC = TWO
    else
       CC = ONE
    end if

    ALLOCATE(SV(n_spin),DSP(n_spin),STAT=status)
    ASSERT(status==0)

    do i = 1, n_spin

       ALLOCATE(eps_full(NSP(i)),eta_full(NSP(i)),STAT = status)
       ASSERT(status == 0)

       !! ASSEMBLE eps and ets
       call resp_util_buildr(gl_N_as_mastr(i_ir,i),&
            gl_N_as_slave(i_ir,i),&
            eps(i_ir,i)%m,eps_full)
       call comm_bcast(eps_full)
       call resp_util_buildr(gl_N_as_mastr(i_ir,i),&
            gl_N_as_slave(i_ir,i),&
            eta(i_ir,i)%m,eta_full)
       call comm_bcast(eta_full)

       ALLOCATE(SV(i)%m(NSP(i)),DSP(i)%m(NSP(i),K),STAT=status)
       ASSERT(status==0)

       SV(i)%m  = SQRT(REAL(CC * eps_full * eta_full,r8_kind))

       DEALLOCATE(eps_full,eta_full, STAT = status)
       ASSERT(status==0)
    end do

    call noRI_DDIV(N, NSP, K, A, SV, DSP)

    ALLOCATE(RM(n_spin*2),STAT=status)
    ASSERT(status==0)

    ALLOCATE(RM(UPUP)%m(NSP(UP),K),RM(UPDN)%m(NSP(UP),K),STAT=status)
    ASSERT(status==0)
    RM(UPUP)%m = ZERO
    RM(UPDN)%m = ZERO
    if (n_spin==2) then
       ALLOCATE(RM(DNUP)%m(NSP(DN),K),RM(DNDN)%m(NSP(DN),K),STAT=status)
       ASSERT(status==0)
       RM(DNDN)%m = ZERO
       RM(DNUP)%m = ZERO
    end if

    FPP_TIMER_STOP(XC_input)

    call noRI_XC(X_KIND,C_KIND,C4,i_ir,RM,NSP,DSP,K)

    FPP_TIMER_START(XC_output)

    do ii = 1, size(RM(UPUP)%m,1)
       RM(UPUP)%m(ii,:) = SV(UP)%m(ii) * RM(UPUP)%m(ii,:)
       RM(UPDN)%m(ii,:) = SV(UP)%m(ii) * RM(UPDN)%m(ii,:)
    end do
    if (n_spin==2) then
       do ii = 1, size(RM(DNDN)%m,1)
          RM(DNUP)%m(ii,:) = SV(DN)%m(ii) * RM(DNUP)%m(ii,:)
          RM(DNDN)%m(ii,:) = SV(DN)%m(ii) * RM(DNDN)%m(ii,:)
       end do
    end if

#if 0
    call show("RM(UPUP)%m",RM(UPUP)%m)
    call show("RM(UPDN)%m",RM(UPDN)%m)
    if (n_spin == 2) then
       call show("RM(DNDN)%m",RM(DNDN)%m)
       call show("RM(DNUP)%m",RM(DNUP)%m)
    end if
#endif

    ALLOCATE(D(N,K),STAT=status)
    ASSERT(status==0)

    !! IT SEEMS I FORGOT TO IMPLEMENT THE S-T CASE (FIXED)
    if ((gl_N_spin==2) .or. (gl_SS)) then
       D(1:NSP(UP),:) = RM(UPUP)%m + RM(UPDN)%m
       if (n_spin == 2) D(NSP(UP)+1:N,:) = RM(DNDN)%m + RM(DNUP)%m
    elseif (gl_ST) then
       D(1:NSP(UP),:) = RM(UPUP)%m - RM(UPDN)%m
    else
       ABORT('DEAD BRANCH')
    end if

    counter = 1  
    DO j=1,K
       RES(counter:counter+N-1) = RES(counter:counter+N-1) + D(:,j)
       counter = counter + N
    END DO

    do i = 1, n_spin
       DEALLOCATE(SV(i)%m,DSP(i)%m, STAT=status)
       ASSERT(status==0)
    end do
    DEALLOCATE(SV,DSP,D, STAT=status)
    ASSERT(status==0)

    call orbital_free(orbs_ob)

    FPP_TIMER_STOP(XC_output)
    FPP_TIMER_STOP(XC_all)

    WRITE (*,*) MyID, "  RI TIMER      "
    tt = FPP_TIMER_VALUE(XC_all)
    WRITE (*,*) "   SUMMARY            ", tt
    tt = FPP_TIMER_VALUE(XC_input)
    WRITE (*,*) "   |- PREPARATION     ", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym)
    WRITE (*,*) "   |- XC  PROD. BAS.SYM. ", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym_OC)
    WRITE (*,*) "      |- XC  PROD. BAS.SYM. ORB.  CALC.", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym_PC)
    WRITE (*,*) "      |- XC  PROD. BAS.SYM. PHI.  CALC.", tt
    tt = FPP_TIMER_VALUE(XC_pb_sym_dRC)
    WRITE (*,*) "      \- XC  PROD. BAS.SYM. dRHO. CALC.", tt
    tt = FPP_TIMER_VALUE(XC_functional)
    WRITE (*,*) "   |- XC  FUNCTIONAL  ", tt
    tt = FPP_TIMER_VALUE(XC_output)
    WRITE (*,*) "   |- OUTPUT          ", tt
    tt = FPP_TIMER_VALUE(XC_int)
    WRITE (*,*) "   \- XC  INTEGRATION ", tt
    tt = FPP_TIMER_VALUE(XC_int_dV)
    WRITE (*,*) "      |-  XC deltaV PART ", tt
    tt = FPP_TIMER_VALUE(XC_int_NB)
    WRITE (*,*) "      |-  XC NOBLAS PART ", tt
    tt = FPP_TIMER_VALUE(XC_int_B)
    WRITE (*,*) "      \-  XC   BLAS PART ", tt

  end subroutine noRI_4C
  !*************************************************************

  !*************************************************************
  subroutine noRI_XC(X_KIND,C_KIND,X,i_ir,M4i,NSP,DSP,K)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use iounitadmin_module,       only: write_to_output_units,write_to_trace_unit
    use unique_atom_module,       only: unique_atoms
    use constants
    use comm_module
    use comm, only: comm_bcast
    use density_calc_module!!,      only: density_calc_nl
    use xpack
    use exchange
    use ch_response_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind)                      :: X, i_ir, X_KIND, C_KIND
    type   (arrmat2), intent(inout)       :: M4i(:) ! assumed to be allocated
    integer(i4_kind)            ,OPTIONAL :: NSP(:)
    type   (arrmat2),intent( IN),OPTIONAL :: DSP(:)
    integer(i4_kind)            ,OPTIONAL :: K
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------

    integer(i4_kind)                    :: vla, i_atom, n_spin
    integer(i4_kind)                    :: i_sp, i_ira
    integer(i4_kind)                    :: status

    real(kind=r8_kind), dimension(vec_length,2) :: rho
    real(kind=r8_kind), dimension(vec_length,2) :: nn

    real(kind=r8_kind) :: gamma(vec_length,3), gw(vec_length),grarho(vec_length,3,2)

    !!    real(kind=r8_kind), dimension(vec_length,2) :: &
    !!         dFdn

    real(kind=r8_kind), dimension(vec_length,4) :: &
         dFdndn, dFdg

    real(kind=r8_kind), dimension(vec_length,6) :: &
         dFdndg, dFdgdg

    ! gridpoints and weigths
    real(kind=r8_kind),pointer      :: grdpts(:,:),grdwts(:)

    integer(i4_kind) :: i_dm, idx, npc, n_irr, nk, npa, ndim, i_nk, i,j,s

    real(r8_kind), ALLOCATABLE :: deltaV(:,:), deltaRHO(:,:)

    real(r8_kind), ALLOCATABLE :: deltaVg(:,:,:), gdeltaRHO(:,:,:), gdeltaRHO2(:,:,:)

    real(r8_kind)    :: CG2(2),CG3(4),temp(vec_length)
    integer(i4_kind) :: GG2(2),DG2(2),DG3(4),GG3(4),GG4(4)

    logical :: XLDA,CLDA,LDAON

    !------------ Executable code --------------------------------

    vla    = vec_length                      !! GRID
    n_spin = symmetry_data_n_spin()          !! SPIN
    npc    = symmetry_data_n_partners(i_ir)  !! PARTNER
    n_irr  = symmetry_data_n_irreps()        !! IRREP
    nk     = size(DSP(1)%m,2)                !! 

    select case (X_KIND)
    case (X_NONE)
       XLDA = .true.
    case (X_XALPHA)
       XLDA = .true.
    case default
       XLDA = .false.
    end select

    select case (C_KIND)
    case (C_NONE)
       CLDA = .true.
    case (C_VWN )
       CLDA = .true.
    case (C_PWlda)
       CLDA = .true.
    case default
       CLDA = .false.
    end select

    if (XLDA .and. CLDA) then 
       LDAON = .true.
    else
       LDAON = .false.
    end if

    !! **** ALLOCATION BLOCK ****

    !! **** ALLOCATE OF deltaV ****
    allocate(deltaV(vla,npc),deltaRHO(vla,npc),STAT=status)
    ASSERT(status==0)

    if(.not. LDAON) then
       allocate(gdeltaRHO(vla,npc,3),gdeltaRHO2(vla,npc,n_spin),STAT=status)
       ASSERT(status==0)
       allocate(deltaVg(vla,npc,3),STAT=status)
       ASSERT(status==0)
    end if

    !! **** ALLOCATION OF PHI ****
    allocate(phi(n_irr),STAT=status)
    ASSERT(status==0)

    if (.not. LDAON) then
       allocate(gphi(n_irr),STAT=status)
       ASSERT(status==0)
    end if

    do i_ira = 1, n_irr
       npa  = symmetry_data_n_partners(i_ira)
       ndim = ssym%dim(i_ira)

       allocate(phi(i_ira)%m(vla, ndim, npa, n_spin), STAT=status)
       ASSERT(status==0)
       phi(i_ira)%m = zero

       if (.not. LDAON) then
          allocate(gphi(i_ira)%m(vla,ndim,npa,3,n_spin),STAT=status)
          ASSERT(status==0)
          gphi(i_ira)%m = zero
       end if
    end do

    !
    ! Loop over grid points:
    !
    call grid_loop_setup_atom()

    grid_points_: do while( more_grid_atom(vec_length, i_atom, grdpts, grdwts)) ! loop over gridpoints
       ASSERT(i_atom/=0)

       ! get vector_length portion of the grid on this processor
       vla=size(grdpts,1)

       FPP_TIMER_START(XC_pb_sym)

       ! calculate gridweights
       gw(1:vla) = grdwts(1:vla) &
            &    * atomicweight(i_atom, grdpts(1:vla, :)) &
            &    * unique_atoms(i_atom)%n_equal_atoms


       FPP_TIMER_START(XC_pb_sym_OC)
       call orbital_calculate(grdpts(1:vla,1:3),vla,orbs_ob,orbs_grads)
       FPP_TIMER_STOP(XC_pb_sym_OC)

       FPP_TIMER_STOP(XC_pb_sym)

       FPP_TIMER_START(XC_functional)

       !! **** <<<< RHO CALCULATION ****
       nn     = zero
       rho    = zero
       gamma  = zero
       grarho = zero

       call density_calc_nl(vla,nn,gamma,grarho,orbs_ob,orbs_grads)

       if (n_spin==1) then
          rho(:,1) = nn(:,1)/TWO
          rho(:,2) = nn(:,1)/TWO

          gamma(:,2) = gamma(:,1)/FOUR
          gamma(:,3) = gamma(:,1)/FOUR
          gamma(:,1) = gamma(:,1)/FOUR
       else
          rho(:,1) = nn(:,1)
          rho(:,2) = nn(:,2)
       end if

       !! CUT OFF
       where ( abs(rho) < 1.0E-16_r8_kind ) !! rho_cutoff
          rho = 1.0E-16_r8_kind !!rho_cutoff
       endwhere
       where ( abs(gamma) < 1.0E-32_r8_kind ) !! rho_cutoff
          gamma = 1.0E-32_r8_kind !!rho_cutoff
       endwhere
       !! **** >>>> RHO CALCULATION ****

       !! **** <<<< XC KERNEL CALCULATION ****
       dFdg   = zero
       dFdndn = zero
       dFdndg = zero
       dFdgdg = zero

       call noRI_XC_calc(X_KIND,C_KIND,X,n_spin,vla,gw,rho,gamma,dFdg,dFdndn,dFdndg,dFdgdg)
       !! **** >>>> XC KERNEL CALCULATION ****

       FPP_TIMER_STOP(XC_functional)

       !! **** <<<< calculation of PHI
       FPP_TIMER_START(XC_pb_sym)
       FPP_TIMER_START(XC_pb_sym_PC)

       select case(X)
       case (C2)
          call fit_fct_calculate_response(grdpts(1:vla,1:3),vla,fcts_orbs_ch, grads_ch)
       case (C4)
          if (LDAON) then
             call noRI_phicalc(vla, n_irr, n_spin, orbs_ob, phi)
          else
             call noRI_phicalc(vla, n_irr, n_spin, orbs_ob, phi, orbs_grads, gphi)
          end if
       end select

       FPP_TIMER_STOP(XC_pb_sym_PC)
       FPP_TIMER_STOP(XC_pb_sym)

       !! **** <<<< MAIN CALCULATION ****

       FPP_TIMER_START(XC_int)
       i_nk_: do i_nk = 1, nk

          !! sigma'
          i_dm_:  do i_dm = 1, 2
             !!                i_spp = i_dm
             !!                if (n_spin==1) i_spp = 1

             if (((n_spin==1) .and. (i_dm == 1)) .or. (n_spin==2)) then

                !! **** CALCULATE RESPONSE OF THE POTENTIAL - DELTAV ****
                FPP_TIMER_STOP(XC_int)
                FPP_TIMER_START(XC_pb_sym)
                FPP_TIMER_START(XC_pb_sym_dRC)

                if (LDAON) then
                   call noRI_deltaRHOcalc(LDAON, i_ir,  i_dm, vla,& !! i_spp
                        DSP(i_dm)%m(:,i_nk), phi, deltaRHO)
                else
                   call noRI_deltaRHOcalc(LDAON, i_ir,  i_dm, vla,& !! i_spp
                        DSP(i_dm)%m(:,i_nk), phi, deltaRHO, gphi, gdeltaRHO)
                end if

                FPP_TIMER_STOP(XC_pb_sym_dRC)
                FPP_TIMER_STOP(XC_pb_sym)
                FPP_TIMER_START(XC_int)

             end if

             !! sigma
             i_sp_: do i_sp = 1,n_spin

                idx = (i_sp - 1) * 2 + i_dm

                deltaV = multderF(dFdndn(:,idx),deltaRHO)

                if (.not. LDAON) then
                   !! multiplication of dFdg x gra(phia x phis) x Xas

                   !! closed shell
                   CG2 = (/ONE, HALF/)
                   CG3 = (/ONE, HALF, HALF, (ONE/FOUR) /)
                   GG2 = (/  A,    A/)
                   GG3 = (/  A,    A,    A,    A/)
                   GG4 = (/  A,    A,    A,    A/)
                   select case (idx)
                   case (AA)
                      DG2 = (/AAA,  AAB/)
                      DG3 = (/AAAA, AAAB, AAAB, ABAB/)
                   case (AB)
                      DG2 = (/ABB,  AAB/)
                      DG3 = (/AABB, AAAB, AAAB, ABAB/)
                   end select

                   if(n_spin== 2) then !! open shell
                      CG2 = (/TWO, ONE/) !! dndg
                      CG3 = (/FOUR,  TWO,  TWO,  ONE/)
                      select case (idx)
                      case (AA)
                         DG2 = (/AAA, AAB/)
                         GG2 = (/  A,   B/)

                         DG3 = (/AAAA, AAAB, AAAB, ABAB/)
                         GG3 = (/   A,    A,    B,    B/)
                         GG4 = (/   A,    B,    A,    B/)
                      case (AB)
                         DG2 = (/ABB, AAB/)
                         GG2 = (/  B,   A/)

                         DG3 = (/AABB, AAAB, BBAB, ABAB/)
                         GG3 = (/   A,    A,    B,    A/)
                         GG4 = (/   B,    A,    B,    B/)
                      case (BB)
                         DG2 = (/BBB, BAB/)
                         GG2 = (/  B,   A/)

                         DG3 = (/BBBB, BBAB, BBAB, ABAB/)
                         GG3 = (/   B,    B,    A,    A/)
                         GG4 = (/   B,    A,    B,    A/)
                      case (BA)
                         DG2 = (/ABB, AAB/)
                         GG2 = (/  B,   A/)

                         DG3 = (/AABB, AAAB, BBAB, ABAB/)
                         GG3 = (/   A,    A,    B,    A/)
                         GG4 = (/   B,    A,    B,    B/)
                      end select
                   end if

                   !! scalar multiplication of (gRHO,gdeltaRHO)

                   gdeltaRHO2 = 0.0_r8_kind

                   do i = 1,3
                      do s = 1, n_spin
                         gdeltaRHO2(1:vla,1:npc,s) = &
                              gdeltaRHO2(1:vla,1:npc,s) + multderF(grarho(:,i,s),gdeltaRHO(:,:,i))
                      end do
                   end do

                   !! calculation Rlda + c

                   do i = 1,2
                         deltaV(:,:) = deltaV(:,:) &
                              + multderF(CG2(i) * dFdndg(:,DG2(i)),gdeltaRHO2(:,:,GG2(i)))
                      end do

                      !! calculation vec(a)+vec(b)+vec(d)

                      do i = 1,3
                         deltaVg(:,:,i) = multderF(TWO * dFdg(:,idx),gdeltaRHO(:,:,i))

                         do j = 1,2

                            temp(:) = CG2(j) * grarho(:,i,GG2(j)) * dFdndg(:,DG2(j))

                            deltaVg(:,:,i) = deltaVg(:,:,i) + multderF(temp,deltaRHO)
                         end do

                         do j = 1,4

                            temp(:) = CG3(j) * grarho(:,i,GG4(j)) * dFdgdg(:,DG3(j))

                            deltaVg(:,:,i) = deltaVg(:,:,i) + multderF(temp,gdeltaRHO2(:,:,GG3(j)))
                         end do

                      end do

                   end if

                   !! **** INTEGRATION ****
                   if (LDAON) then
                      call noRI_XCCOMBINE_4C(LDAON, i_ir, i_sp, vla, &
                           deltaV, phi, M4i(idx)%m(:,i_nk))
                   else
                      call noRI_XCCOMBINE_4C(LDAON, i_ir, i_sp, vla, &
                           deltaV, phi, M4i(idx)%m(:,i_nk), &
                           gphi, deltaVg)
                   end if

                end do i_sp_
             end do i_dm_
          end do i_nk_

          FPP_TIMER_STOP(XC_int)

       end do grid_points_ ! loop over gridpoints

    FPP_TIMER_START(XC_output)
    ! collect all results on the master
    if(comm_i_am_master()) then
       ! receiving results from the slaves

       if (comm_parallel()) then
!!$          call write_to_output_units('noRI_XC: receiving results from the slaves')
          call noRI_receive(n_spin,M4i)
       end if

    else
!!$       call write_to_output_units('noRI_XC: sending the result to the master')
       ! sending the result to the master
       call noRI_tomaster(n_spin,M4i)
    end if

#if 0
    if (comm_i_am_master()) then 
       call show("M4i(AA)",M4i(AA)%m)
       call show("M4i(AB)",M4i(AB)%m)
       if (n_spin==2) then
          call show("M4i(BB)",M4i(BB)%m)
          call show("M4i(BA)",M4i(BA)%m)
       end if
    end if
#endif

    do idx = 1,size(M4i)
       call comm_bcast(M4i(idx)%m)
    end do

    FPP_TIMER_STOP(XC_output)

    ! **** DEALLOCATION BLOCK ****

    !! **** DEALLOCATE deltaV, deltaRHO ****
    DEALLOCATE(deltaV,deltaRHO,STAT=status)
    ASSERT(status==0)

    if (.not. LDAON) then
       deallocate(gdeltaRHO,gdeltaRHO2,STAT=status)
       ASSERT(status==0)
       deallocate(deltaVg,STAT=status)
       ASSERT(status==0)
    end if

    !! **** DEALLOCATION OF PHI ****
    do i_ira = 1, n_irr
       deallocate(phi(i_ira)%m, STAT=status)
       ASSERT(status==0)
    end do

    deallocate(phi, STAT=status)
    ASSERT(status==0)

    if (.not. LDAON) then
       do i_ira = 1, n_irr
          deallocate(gphi(i_ira)%m, STAT=status)
          ASSERT(status==0)
       end do

       deallocate(gphi, STAT=status)
       ASSERT(status==0) 
    end if

  end subroutine noRI_XC
  !*************************************************************

  function multderF(a,b) result(c)
    real(kind=r8_kind) :: a(:), b(:,:)
    real(kind=r8_kind) :: c(size(b,1),size(b,2))

    integer(i4_kind)   :: i

    FPP_TIMER_START(XC_int_dV)
    do i = 1, size(b,2)
       c(:,i) = a(:) * b(:,i)
    end do
    FPP_TIMER_STOP(XC_int_dV)
  end function multderF

  !*************************************************************
  subroutine noRI_deltaRHOcalc(LDAON,i_ir_c,i_sp,vla,DSP,&
       phi,deltaRHO,gphi,gdeltaRHO)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use constants
    implicit none
    !------------ Declaration of formal parameters ---------------
    logical,                intent(IN )   :: LDAON
    integer(i4_kind),       intent(IN )   :: i_ir_c, vla, i_sp
    type(arrmat4), target,  intent(IN )   :: phi (:)
    real(r8_kind),          intent(IN )   :: DSP(:)
    real(r8_kind),          intent(INOUT) :: deltaRHO(:,:)
    type(arrmat5),optional,target, intent(IN )   :: gphi(:)
    real(r8_kind),optional, intent(INOUT) :: gdeltaRHO(:,:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)      :: n_irr, as, as_tmp
    integer(kind=i4_kind)      :: i_ir_a, i_ir_s
    integer(kind=i4_kind)      :: oo_dim

    real(kind=r8_kind),pointer :: phi_a(:,:,:), phi_s(:,:,:) !! (vl,exp,npc,sp)
    real(kind=r8_kind),pointer :: gphi_a(:,:,:,:), gphi_s(:,:,:,:) !!(vl,exp,npc,1:3,sp)

    integer(kind=i4_kind)      :: occs,occe,unoccs,unocce
    integer(kind=i4_kind)      :: i_mlt, na, ns, i

    real(r8_kind),pointer      :: pcg(:,:,:)
    !------------ Executable code --------------------------------

    n_irr  = symmetry_data_n_irreps() ! number of irreps

    as = 0
    deltaRHO = 0.0
    if (.not. LDAON) gdeltaRHO = 0.0

    i_ir_a_: do i_ir_a = 1, n_irr
       phi_a => phi(i_ir_a)%m(:,:,:,i_sp)
       if (.not. LDAON) gphi_a => gphi(i_ir_a)%m(:,:,:,:,i_sp)

       i_ir_s_: do i_ir_s = 1, n_irr
          phi_s => phi(i_ir_s)%m(:,:,:,i_sp)
          if (.not. LDAON) gphi_s => gphi(i_ir_s)%m(:,:,:,:,i_sp)

          call resp_util_calc_transitions_v2(i_ir_a, i_ir_s, i_ir_c, &
               i_sp, oo_dim)
          if (oo_dim == 0) cycle

          call resp_util_borders(i_ir_a,i_ir_s,i_sp,&
               &                 occs,occe,unoccs,unocce)

          na = occe - occs + 1
          ns = unocce - unoccs + 1

          i_mlt_: do i_mlt = 1, cg(i_ir_c,i_ir_a,i_ir_s)%mult
             pcg => cg(i_ir_c,i_ir_a,i_ir_s)%sub(i_mlt)%c

             as_tmp = as
             call symRHO2(vla, na, ns, pcg, &
                  phi_a(:, occs:, :), phi_s, unoccs, DSP, as, deltaRHO)

             if (.not. LDAON) then
                do i = 1,3
                   as = as_tmp
                   call symRHO2(vla,na,ns,pcg,&
                        gphi_a(:,occs:,:,i),phi_s,unoccs,DSP,as,gdeltaRHO(:,:,i))
                   as = as_tmp
                   call symRHO2(vla,na,ns,pcg,&
                        phi_a(:,occs:,:),gphi_s(:,:,:,i),unoccs,DSP,as,gdeltaRHO(:,:,i))
                end do
             end if

          end do i_mlt_
       end do i_ir_s_
    end do i_ir_a_

  end subroutine noRI_deltaRHOcalc
  !*************************************************************

  !*************************************************************
  subroutine symRHO2(vla, na, ns, pcg, phia, phis, us, Uin, as, Uout)
    use cpks_grid_utils, only: rho2
    use constants
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: vla,na,ns
    real(kind=r8_kind),    intent(IN   ) :: pcg(:,:,:)
    real(kind=r8_kind),    intent(IN   ) :: phia(:,:,:), phis(:,:,:)
    real(kind=r8_kind),    intent(IN   ) :: Uin(:)
    integer(kind=i4_kind), intent(IN   ) :: us
    integer(kind=i4_kind), intent(INOUT) :: as
    real(kind=r8_kind),    intent(INOUT) :: Uout(:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i_pa_c,i_pa_a,i_pa_s
    integer(kind=i4_kind) :: npa,nps,npc
    integer(kind=i4_kind) :: a, s
    real(kind=r8_kind)    :: Xas(na,ns)
    real(kind=r8_kind)    :: r1(vla)
    !------------ Executable code --------------------------------
    npc = size(pcg,1)
    npa = size(pcg,2)
    nps = size(pcg,3)

    i_a_: do a = 1, na
       i_s_: do s = 1, ns
          as = as + 1
          Xas(a,s) = Uin(as)
       end do i_s_
    end do i_a_

    i_pa_a_: do i_pa_a = 1, npa
       i_pa_s_: do i_pa_s = 1, nps
          r1(:) =  rho2(vla, us-1, phia(:, :, i_pa_a), phis(:, :, i_pa_s), Xas)
          i_pa_c_: do i_pa_c = 1, npc

             if ( abs(pcg(i_pa_c, i_pa_a, i_pa_s)) < 1e-7 ) CYCLE

             Uout(:vla, i_pa_c) = Uout(:vla, i_pa_c) + &
                  pcg(i_pa_c, i_pa_a, i_pa_s) / 2.0_r8_kind * r1(:vla)

          end do i_pa_c_
       end do i_pa_s_
    end do i_pa_a_

  end subroutine symRHO2
  !*************************************************************

  !*************************************************************
  subroutine noRI_XCCOMBINE_4C(LDAON, i_ir_c, i_sp, vla, &
       deltaV, phi, res_mat, gphi, deltaVg )
    !  Purpose: ..
    !------------ Modules used ----------------------------------
    use constants
    implicit none
    !------------ Declaration of formal parameters ---------------
    logical,               intent(IN   )           :: LDAON
    integer(kind=i4_kind), intent(IN   )           :: i_ir_c, i_sp, vla
    real(kind=r8_kind),    intent(IN   )           :: deltaV(:,:)
    type(arrmat4), target, intent(IN   )           :: phi(:)
    real(kind=r8_kind),    intent(INOUT)           :: res_mat(:)
    !! optional because GGA
    type(arrmat5), target, intent(IN   ), optional :: gphi(:)
    real(kind=r8_kind),    intent(IN   ), optional :: deltaVg  (:,:,:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)          :: n_irr, as_tmp, n_spin
    integer(kind=i4_kind)          :: as, i_ir_a
    integer(kind=i4_kind)          ::     i_ir_s, oo_dim
    real(kind=r8_kind),pointer     :: phi_a(:,:,:), phi_s(:,:,:) !! (vl,na,npa)
    real(kind=r8_kind),pointer     :: gphi_a(:,:,:,:), gphi_s(:,:,:,:) !! (vl,na,npa,1:3)
    integer(kind=i4_kind)          :: occs,occe,unoccs,unocce
    integer(kind=i4_kind)          :: na, ns
    integer(kind=i4_kind)          :: i_mlt
    real(kind=r8_kind),pointer     :: pcg(:,:,:)
    !------------ Executable code --------------------------------
    n_irr  = symmetry_data_n_irreps() ! number of irrps
    n_spin = symmetry_data_n_spin()   ! number of spins

    as = 0
    i_ir_a_: do i_ir_a = 1, n_irr
       phi_a => phi(i_ir_a)%m(:,:,:,i_sp)
       if (.not. LDAON) gphi_a => gphi(i_ir_a)%m(:,:,:,:,i_sp)

       i_ir_s_: do i_ir_s = 1, n_irr

          call resp_util_calc_transitions_v2(i_ir_a, i_ir_s, i_ir_c, &
               i_sp, oo_dim)
          if (oo_dim == 0) cycle

          phi_s => phi(i_ir_s)%m(:,:,:,i_sp)
          if (.not. LDAON) gphi_s => gphi(i_ir_s)%m(:,:,:,:,i_sp)

          call resp_util_borders(i_ir_a,i_ir_s,i_sp,&
               &                 occs,occe,unoccs,unocce)

          na = occe   - occs   + 1
          ns = unocce - unoccs + 1

          i_mlt_: do i_mlt = 1, cg(i_ir_c,i_ir_a,i_ir_s)%mult
             pcg  => cg(i_ir_c,i_ir_a,i_ir_s)%sub(i_mlt)%c

             as_tmp = as
             call symGridInt2(vla, na, ns, pcg, &
                  phi_a, occs - 1, phi_s, unoccs - 1, deltaV, as, res_mat)

             if (.not. LDAON) then
#if 0
                do i = 1,3
                   as = as_tmp
                   call symGridInt2(vla,na,ns,pcg,&
                        gphi_a(:,:,:,i),occs - 1,phi_s,unoccs - 1,deltaVg(:,:,i),as,res_mat)
                   as = as_tmp
                   call symGridInt2(vla,na,ns,pcg,&
                        phi_a,occs - 1,gphi_s(:,:,:,i),unoccs - 1,deltaVg(:,:,i),as,res_mat)
                end do
#endif

                as = as_tmp
                call symGridInt2a(vla,na,ns,pcg,&
                     gphi_a(:,:,:,:),occs - 1,phi_s,unoccs - 1,deltaVg(:,:,:),as,res_mat)
                
                as = as_tmp
                call symGridInt2s(vla,na,ns,pcg,&
                     phi_a,occs - 1,gphi_s(:,:,:,:),unoccs - 1,deltaVg(:,:,:),as,res_mat)
                
             end if

          end do i_mlt_
       end do i_ir_s_
    end do i_ir_a_

  end subroutine noRI_XCCOMBINE_4C
  !*************************************************************

  !*************************************************************
  subroutine symGridInt2(vla, na, ns, pcg, phia, oa, phis, os, Uin, as, Uout)
    use cpks_grid_utils, only: gridInt2
    use constants
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: vla,na,ns
    real(kind=r8_kind),    intent(IN   ) :: pcg(:,:,:)
    real(kind=r8_kind),    intent(IN   ) :: phia(:,:,:), phis(:,:,:) 
    real(kind=r8_kind),    intent(IN   ) :: Uin(:,:)
    integer(kind=i4_kind), intent(IN   ) :: oa, os
    integer(kind=i4_kind), intent(INOUT) :: as
    real(kind=r8_kind),    intent(INOUT) :: Uout(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i_pa_c,i_pa_a,i_pa_s
    integer(kind=i4_kind) :: npa,nps,npc
    integer(kind=i4_kind) :: a, s, i, status
    real(kind=r8_kind)    :: LDA0(na,ns)
    real(kind=r8_kind)    :: coeff

    real(kind=r8_kind),allocatable    :: phiv(:,:)
    !------------ Executable code --------------------------------
    npc = size(pcg,1)
    npa = size(pcg,2)
    nps = size(pcg,3)

    LDA0 =  0.0

    if (nps .le. npa) then
       allocate(phiv(vla,na),STAT=status)
       ASSERT(status==0)
       i_pa_s_: do i_pa_s = 1, nps
          phiv = 0.0
          i_pa_a_: do i_pa_a = 1, npa
             i_pa_c_: do i_pa_c = 1, npc
                FPP_TIMER_START(XC_int_NB)
                coeff = pcg(i_pa_c,i_pa_a,i_pa_s)
                if( abs(coeff) < 1e-7 ) CYCLE 

                do i = 1, na
                   phiv(:vla,i) = phiv(:vla,i) &
                        + phia(:vla,oa+i,i_pa_a) &
                        * Uin(:vla,i_pa_c)       &
                        * coeff
                end do
                FPP_TIMER_STOP (XC_int_NB)
             end do i_pa_c_
          end do i_pa_a_

          FPP_TIMER_START(XC_int_B)
          call dgemm( 't', 'n', size(LDA0,1), size(LDA0,2) &
               , vla, 1.0_r8_kind  &
               , phiv(1,1),           size(phiv,1) &
               , phis(1,os+1,i_pa_s), size(phis,1), 1.0_r8_kind &
               , LDA0(1,1) ,          size(LDA0,1)             &
               )
          FPP_TIMER_STOP (XC_int_B)
       end do i_pa_s_
       deallocate(phiv,STAT=status)
       ASSERT(status==0)
    else
       allocate(phiv(vla,ns),STAT=status)
       ASSERT(status==0)
       do i_pa_a = 1, npa
          phiv = 0.0
          do i_pa_s = 1, nps
             do i_pa_c = 1, npc
                FPP_TIMER_START(XC_int_NB)
                coeff = pcg(i_pa_c,i_pa_a,i_pa_s)
                if( abs(coeff) < 1e-7 ) CYCLE 

                do i = 1, ns
                   phiv(:vla,i) = phiv(:vla,i) &
                        + phis(:vla,os+i,i_pa_s) &
                        * Uin(:vla,i_pa_c)       &
                        * coeff
                end do
                FPP_TIMER_STOP (XC_int_NB)
             end do
          end do

          FPP_TIMER_START(XC_int_B)
          call dgemm( 't', 'n', size(LDA0,1), size(LDA0,2) &
               , vla,  1.0_r8_kind&
               , phia(1,oa+1,i_pa_a), size(phia,1) &
               , phiv(1,1),           size(phiv,1), 1.0_r8_kind &
               , LDA0(1,1) ,          size(LDA0,1)             &
               )
          FPP_TIMER_STOP (XC_int_B)
       end do
       deallocate(phiv,STAT=status)
       ASSERT(status==0)
    end if

    FPP_TIMER_START(XC_int_NB)
    i_a_: do a = 1, na
       i_s_: do s = 1, ns
          as = as + 1
          Uout(as) = Uout(as) + LDA0(a,s) / npc
       end do i_s_
    end do i_a_
    FPP_TIMER_STOP (XC_int_NB)

  end subroutine symGridInt2
  !*************************************************************

  !*************************************************************
  subroutine symGridInt2a(vla, na, ns, pcg, gphia, oa, phis, os, Uin, as, Uout)
    use cpks_grid_utils, only: gridInt2
    use constants
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: vla,na,ns
    real(kind=r8_kind),    intent(IN   ) :: pcg(:,:,:)
    real(kind=r8_kind),    intent(IN   ) :: gphia(:,:,:,:), phis(:,:,:) 
    real(kind=r8_kind),    intent(IN   ) :: Uin(:,:,:)
    integer(kind=i4_kind), intent(IN   ) :: oa, os
    integer(kind=i4_kind), intent(INOUT) :: as
    real(kind=r8_kind),    intent(INOUT) :: Uout(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i_pa_c,i_pa_a,i_pa_s
    integer(kind=i4_kind) :: npa,nps,npc
    integer(kind=i4_kind) :: a, s, i, status,r
    real(kind=r8_kind)    :: LDA0(na,ns)
    real(kind=r8_kind)    :: coeff

    real(kind=r8_kind),allocatable    :: phiv(:,:)
    !------------ Executable code --------------------------------
    npc = size(pcg,1)
    npa = size(pcg,2)
    nps = size(pcg,3)

    LDA0 =  0.0

    allocate(phiv(vla,na),STAT=status)
    ASSERT(status==0)
    i_pa_s_: do i_pa_s = 1, nps
       phiv = 0.0
       i_pa_a_: do i_pa_a = 1, npa
          i_pa_c_: do i_pa_c = 1, npc
             FPP_TIMER_START(XC_int_NB)
             coeff = pcg(i_pa_c,i_pa_a,i_pa_s)
             if( abs(coeff) < 1e-7 ) CYCLE 

             do i = 1, na
                do r = 1,3
                   phiv(:vla,i) = phiv(:vla,i) &
                        + gphia(:vla,oa+i,i_pa_a,r) &
                        * Uin(:vla,i_pa_c,r)       &
                        * coeff
                end do
             end do
             FPP_TIMER_STOP (XC_int_NB)
          end do i_pa_c_
       end do i_pa_a_

       FPP_TIMER_START(XC_int_B)
       call dgemm( 't', 'n', size(LDA0,1), size(LDA0,2) &
            , vla, 1.0_r8_kind  &
            , phiv(1,1),           size(phiv,1) &
            , phis(1,os+1,i_pa_s), size(phis,1), 1.0_r8_kind &
            , LDA0(1,1) ,          size(LDA0,1)             &
            )
       FPP_TIMER_STOP (XC_int_B)
    end do i_pa_s_
    deallocate(phiv,STAT=status)
    ASSERT(status==0)

    FPP_TIMER_START(XC_int_NB)
    i_a_: do a = 1, na
       i_s_: do s = 1, ns
          as = as + 1
          Uout(as) = Uout(as) + LDA0(a,s) / npc
       end do i_s_
    end do i_a_
    FPP_TIMER_STOP (XC_int_NB)

  end subroutine symGridInt2a
  !*************************************************************

  !*************************************************************
  subroutine symGridInt2s(vla,na,ns,pcg,phia,oa,gphis,os,Uin,as,Uout)
    use cpks_grid_utils, only: gridInt2
    use constants
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind), intent(IN   ) :: vla,na,ns
    real(kind=r8_kind),    intent(IN   ) :: pcg(:,:,:)
    real(kind=r8_kind),    intent(IN   ) :: phia(:,:,:), gphis(:,:,:,:)
    real(kind=r8_kind),    intent(IN   ) :: Uin(:,:,:)
    integer(kind=i4_kind), intent(IN   ) :: oa, os
    integer(kind=i4_kind), intent(INOUT) :: as
    real(kind=r8_kind),    intent(INOUT) :: Uout(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i_pa_c,i_pa_a,i_pa_s
    integer(kind=i4_kind) :: npa,nps,npc,r
    integer(kind=i4_kind) :: a, s, i, status
    real(kind=r8_kind)    :: LDA0(na,ns)
    real(kind=r8_kind)    :: coeff

    real(kind=r8_kind),allocatable    :: phiv(:,:)
    !------------ Executable code --------------------------------
    npc = size(pcg,1)
    npa = size(pcg,2)
    nps = size(pcg,3)

    LDA0 =  0.0

    allocate(phiv(vla,ns),STAT=status)
    ASSERT(status==0)
    do i_pa_a = 1, npa
       phiv = 0.0
       do i_pa_s = 1, nps
          do i_pa_c = 1, npc
             FPP_TIMER_START(XC_int_NB)
             coeff = pcg(i_pa_c,i_pa_a,i_pa_s)
             if( abs(coeff) < 1e-7 ) CYCLE 

             do i = 1, ns
                do r = 1,3
                   phiv(:vla,i) = phiv(:vla,i) &
                        + gphis(:vla,os+i,i_pa_s,r) &
                        * Uin(:vla,i_pa_c,r)       &
                        * coeff
                end do
             end do
             FPP_TIMER_STOP (XC_int_NB)
          end do
       end do

       FPP_TIMER_START(XC_int_B)
       call dgemm( 't', 'n', size(LDA0,1), size(LDA0,2) &
            , vla,  1.0_r8_kind&
            , phia(1,oa+1,i_pa_a), size(phia,1) &
            , phiv(1,1),           size(phiv,1), 1.0_r8_kind &
            , LDA0(1,1) ,          size(LDA0,1)             &
            )
       FPP_TIMER_STOP (XC_int_B)
    end do
    deallocate(phiv,STAT=status)
    ASSERT(status==0)

    FPP_TIMER_START(XC_int_NB)
    i_a_: do a = 1, na
       i_s_: do s = 1, ns
          as = as + 1
          Uout(as) = Uout(as) + LDA0(a,s) / npc
       end do i_s_
    end do i_a_
    FPP_TIMER_STOP (XC_int_NB)

  end subroutine symGridInt2s
  !*************************************************************

  !*************************************************************
  subroutine noRI_receive(n_spin,M4i)
    !  Purpose: receive the results from the slaves
    !  Values of 2-index-integrals on the grid portion of
    !  accessable to each slave are collected from the slaves    
    !  and added up to yield the final matrix
    !------------ Modules used ------------------- ---------------
    use xpack,       only: upck
    use constants,   only: zero
    use comm_module
    implicit none
    !** End of interface *****************************************
    integer(i4_kind), INTENT(IN)     :: n_spin
    type(arrmat2),    INTENT(INOUT)  :: M4i(:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: status,n_procs,i_proc
    integer(kind=i4_kind) :: i_spin,i_dim,dim_ou_a,dim_ou_b
    integer(kind=i4_kind) :: idx
    real(kind=r8_kind),allocatable :: help_mat(:,:)
    !------------ Executable code --------------------------------


    ! loop over all slave processors
    n_procs=comm_get_n_processors()
    do i_proc=2,n_procs

       do i_spin = 1,n_spin
          do i_dim=1, 2 !! dim_factor
             idx = (i_spin-1)*2+i_dim
             dim_ou_a = size(M4i(idx)%m,1)
             dim_ou_b = size(M4i(idx)%m,2)
             allocate(help_mat(dim_ou_a,dim_ou_b),STAT=status)
             ASSERT(status==0)
             help_mat = zero

             call comm_save_recv(i_proc,msgtag_response_2in_send)
             if(comm_msgtag()/=msgtag_response_2in_send) &
                  call error_handler('Wrong msgtag in response_2index_receive')
             call upck(help_mat)

             ! add it to the final matrix
             M4i(idx)%m = M4i(idx)%m + help_mat

             deallocate(help_mat,STAT=status)
             ASSERT(status==0)  
          end do
       end do
    end do

  end subroutine noRI_receive
  !*************************************************************

  !*************************************************************
  subroutine noRI_tomaster(n_spin,M4i)
    !  Purpose: send results from the slave to the master
    !  Values of 2-index-integrals on the grid portion of
    !  accessable to each slave are send to the master        
    !------------ Modules used ------------------- ---------------
    use xpack, only: pck
    use comm_module
    implicit none
    !** End of interface *****************************************
    integer(i4_kind), INTENT(IN)     :: n_spin
    type(arrmat2),    INTENT(INOUT)  :: M4i(:)
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind) :: i_dim,i_spin
    integer(kind=i4_kind) :: idx
    !------------ Executable code --------------------------------

    do i_spin = 1,n_spin
       do i_dim=1, 2 !! dim_factor
          idx = (i_spin-1)*2+i_dim
          call comm_init_send(comm_master_host,msgtag_response_2in_send)
          call pck(M4i(idx)%m)
          call comm_send()
       end do
    end do

  end subroutine noRI_tomaster
  !*************************************************************



  !*************************************************************
  subroutine noRI_totape(X,i_ir,n_spin,M4i)
    !  Purpose: 
    !  Write 2index matrix row by row to a linear tape.
    !------------ Modules used ------------------- ---------------
    !!    use readwriteblocked_module
    use io,                 only: write_buffer
    use filename_module ,   only: resp_dir
    use resp_util_module,   only: resp_util_fname
    implicit none
    integer(i4_kind), INTENT(IN)  :: X, i_ir, n_spin
    type(arrmat2),    INTENT(IN)  :: M4i(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    !!    type(readwriteblocked_tapehandle) :: th_2index
    integer(kind=i4_kind)             :: i_spin,i_dim,idx
    !------------ Executable code --------------------------------

    do i_spin=1, n_spin
       do i_dim=1, 2
          idx = (i_spin-1)*2+i_dim
          select case(X)
          case(C2)
             call write_buffer( trim(resp_dir)//'/'//resp_util_fname('xc_2c',i_ir,idx) &
                  , triconvert(M4i(idx)%m)                                 &
                  )
          case(C4)
             call write_buffer( trim(resp_dir)//'/'//resp_util_fname('XC_4i',i_ir,idx) &
                  , M4i(idx)%m                                             &
                  )
          case default
             ABORT('noRI: NO SUCH CASE')
          end select
       end do
    end do

  contains
    function triconvert(QM) result(TM)
      ! convert from quadratic matrix (QM) to triangle (TM)
      implicit none
      real(r8_kind), intent(in)  :: QM(:,:)
      real(r8_kind)              :: TM( (size(QM,1)*(size(QM,1)+1))/2 )
      integer(i4_kind)           :: idx, i, j
      ! *** end of interface ***
      ASSERT(size(QM,1)==size(QM,2))
      idx = 0
      do i = 1, size(QM,1)
         do j = 1, i
            idx = idx + 1
            TM(idx) = QM(i,j) !!(QM(i,j) + QM(j,i)) / 2
         end do
      end do
    end function triconvert
  end subroutine noRI_totape
  !*************************************************************


  !*************************************************************
  subroutine noRI_XC_calc(X_KIND,C_KIND,X,n_spin,vla,gw,rho,gamma,dFdg,dFdndn,dFdndg,dFdgdg)
    !  Purpose: calculate XC functional
    !------------ Modules used ------------------- ---------------
    use exchange
    use gga_response_module
    use constants
    use pw_ldac_module
    use vwnc
    use becke_perdew_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(kind=i4_kind),             intent(IN ) :: X_KIND,C_KIND,X,n_spin,vla
    real(kind=r8_kind),                intent(IN ) :: gw(:), rho(:,:), gamma(:,:)
    real(kind=r8_kind),dimension(:,:), intent(OUT) :: dFdg,dFdndn,dFdndg,dFdgdg
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    real(kind=r8_kind), dimension(vla  ) :: &
         Fx, Fc, Flda, eps
    real(kind=r8_kind), dimension(vla,2) :: &
         dFxdn, dFcdn, dFldadn
    real(kind=r8_kind), dimension(vla,3) :: &
         dFxdg, dFcdg, dFxdndn, dFcdndn, dFldadndn
    real(kind=r8_kind), dimension(vla,6) :: &
         dFxdndg, dFcdndg
    real(kind=r8_kind), dimension(vla,6) :: &
         dFxdgdg, dFcdgdg
    integer(kind=i4_kind) :: i
    real(kind=r8_kind)    :: CC
    !------------ Executable code --------

    !! EXCHANGE PART
    Flda      = zero
    dFldadn   = zero
    dFldadndn = zero

    Fx        = zero
    dFxdn     = zero
    dFxdg     = zero
    dFxdndn   = zero
    dFxdndg   = zero
    dFxdgdg   = zero

    select case(X_KIND)
    case ( X_XALPHA )  
       call exchange_lda(X_XALPHA,vla,2,rho,Fx,dFxdn,dFxdndn)
    case ( X_BECKE88 )
       call exchange_lda(X_XALPHA ,vla,2,rho,Fx,dFxdn,dFxdndn)
       call exchange_gga(X_BECKE88,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
    case ( X_PBE    )
       call exchange_lda(X_XALPHA,vla,2,rho,Fx,dFxdn,dFxdndn)
       call exchange_gga(X_PBE   ,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
    case ( X_PW91    )
       call exchange_lda(X_XALPHA,vla,2,rho,Fx,dFxdn,dFxdndn)
       call exchange_gga(X_PW91  ,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
    case ( X_PBEN    )
       call exchange_lda(X_XALPHA,vla,2,rho,Fx,dFxdn,dFxdndn)
       call exchange_gga(X_PBEN  ,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
    case ( X_REVPBE  )
       call exchange_lda(X_XALPHA,vla,2,rho,Fx,dFxdn,dFxdndn)
       call exchange_gga(X_REVPBE,vla,2,rho,gamma,Fx,dFxdn,dFxdg,dFxdndn,dFxdndg,dFxdgdg)
    end select

    !! CORRELATION PART        
    Flda      = zero
    dFldadn   = zero
    dFldadndn = zero

    Fc        = zero
    dFcdn     = zero
    dFcdg     = zero
    dFcdndn   = zero
    dFcdndg   = zero
    dFcdgdg   = zero

    select case (C_KIND)
    case ( C_VWN )  
       eps(1:vla) = 1.0E-16_r8_kind
       call vwn_calcMDA(  rho,dFcdn,2,Fc,vla,eps,dFcdndn)
    case ( C_PWlda )
       call pw_ldac   (vla,4,rho,Fc,dFcdn,dFcdndn) !! 4 = unrestr case
    case ( C_PERDEW )
       eps(1:vla) = 1.0E-16_r8_kind
       call vwn_calcMDA(rho,dFldadn,2,Flda,vla,eps,dFldadndn)
       call perdew_calc(rho,gamma,dFcdn,2,Fc,dFcdg, vla, &
            dFcdndn, &
            dFcdndg, &
            dFcdgdg)
       Fc      = Fc      + Flda
       dFcdn   = dFcdn   + dFldadn
       dFcdndn = dFcdndn + dFldadndn
    case ( C_PBE   )
       call pw_ldac   (vla,4,rho,Flda,dFldadn,dFldadndn)
       call gga_correlation(C_PBE,2,2,rho,gamma,vla,&
            Fc,Flda,&
            dFcdn,dFcdg,dFldadn,&
            dFcdndn,dFcdndg,dFcdgdg,dFldadndn)
       Fc      = Fc      + Flda
       dFcdn   = dFcdn   + dFldadn
       dFcdndn = dFcdndn + dFldadndn
    case ( C_PW91   )
       call pw_ldac   (vla,4,rho,Flda,dFldadn,dFldadndn)
       call gga_correlation(C_PW91,2,2,rho,gamma,vla,&
            Fc,Flda,&
            dFcdn,dFcdg,dFldadn,&
            dFcdndn,dFcdndg,dFcdgdg,dFldadndn)
       Fc      = Fc      + Flda
       dFcdn   = dFcdn   + dFldadn
       dFcdndn = dFcdndn + dFldadndn  
    end select

    if (n_spin == 1) then
       CC = TWO
    else
       select case(X)
       case(C2)
          CC = ONE
       case(C4)
          CC = ONE / TWO
       end select
    end if

    dFxdg   = (dFxdg  +dFcdg  ) / CC
    dFxdndn = (dFxdndn+dFcdndn) / CC
    dFxdndg = (dFxdndg+dFcdndg) / CC
    dFxdgdg = (dFxdgdg+dFcdgdg) / CC

    dFdg   = zero
    dFdndn = zero

    select case (n_spin)
    case (1)
       dFdg  (1:vla,AA) = gw(1:vla) * dFxdg(1:vla,1)
       dFdg  (1:vla,AB) = gw(1:vla) * dFxdg(1:vla,3) / TWO !! will be x2 later

       dFdndn(1:vla,AA) = gw(1:vla) * dFxdndn(1:vla,XCAA)
       dFdndn(1:vla,AB) = gw(1:vla) * dFxdndn(1:vla,XCAB)
    case (2) 
       dFdg  (1:vla,AA) = gw(1:vla) * dFxdg(1:vla,XCAA)   
       dFdg  (1:vla,AB) = gw(1:vla) * dFxdg(1:vla,XCAB) / TWO !! will be x2 later
       dFdg  (1:vla,BB) = gw(1:vla) * dFxdg(1:vla,XCBB)   
       dFdg  (1:vla,BA) = gw(1:vla) * dFxdg(1:vla,XCAB) / TWO !! will be x2 later

       dFdndn(1:vla,AA) = gw(1:vla) * dFxdndn(1:vla,XCAA)
       dFdndn(1:vla,AB) = gw(1:vla) * dFxdndn(1:vla,XCAB)
       dFdndn(1:vla,BB) = gw(1:vla) * dFxdndn(1:vla,XCBB)
       dFdndn(1:vla,BA) = gw(1:vla) * dFxdndn(1:vla,XCAB) 
    end select

    dFdndg = zero
    dFdgdg = zero

    do i = 1,size(dFdndg,2)
       dFdndg(1:vla,i)  = gw(1:vla) * dFxdndg(1:vla,i) 
       dFdgdg(1:vla,i)  = gw(1:vla) * dFxdgdg(1:vla,i) 
    end do

  end subroutine noRI_XC_calc
  !*************************************************************

  !*************************************************************
  subroutine noRI_phicalc(vla, n_irr, n_spin, orbs, phi, gorbs, gphi)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use eigen_data_module,      only: eigvec 
    use constants,              only: ZERO,ONE
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),                    intent(IN ) :: vla,n_irr,n_spin
    type(orbital_type),                  intent(IN ) :: orbs(:)
    type(arrmat4),  target,              intent(inout) :: phi(:) ! comes allocated
    type(orbital_gradient_type),optional,intent(IN ) :: gorbs(:)
    type(arrmat5),optional, target,      intent(inout) :: gphi(:) ! comes allocated
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)      :: i_ir_a, i_sp, npa, i_pa,i
    real(kind=r8_kind),pointer :: eigv(:,:,:)
    real(kind=r8_kind),pointer :: phi_p(:,:,:,:), orbs_p(:,:,:)
    real(kind=r8_kind),pointer :: gphi_p(:,:,:,:,:), gorbs_p(:,:,:,:)
    !------------ Executable code --------------------------------
    do i_ir_a = 1, n_irr

       npa     =  symmetry_data_n_partners(i_ir_a)
       orbs_p  => orbs(i_ir_a)%o
       if (present(gorbs)) gorbs_p => gorbs(i_ir_a)%o
       eigv    => eigvec(i_ir_a)%m

       phi(i_ir_a)%m = zero
       phi_p => phi(i_ir_a)%m

       if (present(gorbs)) then
          gphi(i_ir_a)%m = zero
          gphi_p => gphi(i_ir_a)%m
       end if

       i_pa_: do i_pa = 1, npa
          i_sp_: do i_sp = 1, n_spin
!            phi_p(1:vla, :, i_pa, i_sp) = matmul(orbs_p(1:vla, :, i_pa), eigv(:, :, i_sp))
             call matmatmul(orbs_p(1:vla, :, i_pa), eigv(:, :, i_sp), &
                  phi_p(1:vla, :, i_pa, i_sp), "N", "N", ONE)

             if (present(gorbs)) then
                do i =1,3
                   call matmatmul(gorbs_p(1:vla, i, :, i_pa), eigv(:, :, i_sp), &
                        gphi_p(1:vla, :, i_pa, i, i_sp), "N", "N", ONE)
                end do
             end if

          end do i_sp_
       end do i_pa_
    end do
  end subroutine noRI_phicalc
  !*************************************************************

  !*************************************************************
  subroutine noRI_DDIV(N, NSP, K, D, EPSETA, DSP)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(i4_kind),  intent(   IN) :: N, NSP(:), K
    real   (r8_kind),  intent(   IN) :: D(:)
    type   (arrmat1),  intent(   IN) :: EPSETA(:)
    type   (arrmat2),  intent(INOUT) :: DSP(:)
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)      :: j, counter, ias, ias_dn
    !------------ Executable code -------------------------------- 

    counter = 1

    !! DIVIDE DSP(UP) = sqrt(eps*eta)[UP] * D(UP,N)
    !!        DSP(DN) = sqrt(eps*eta)[DN] * D(DN,N)

    DO j=1,K
       do ias = 1, N
          if (ias <= NSP(UP)) then
             DSP(UP)%m(ias   ,j) = D(counter+ias-1) * EPSETA(UP)%m(ias)
          else              
             ias_dn =ias - NSP(UP)
             DSP(DN)%m(ias_dn,j) = D(counter+ias-1) * EPSETA(DN)%m(ias_dn)
          end if
       end do
       counter = counter + N
    END DO

  end subroutine noRI_DDIV
  !*************************************************************

!!$  !*************************************************************
!!$  subroutine noRI_
!!$    !  Purpose: ..
!!$    !------------ Modules used ------------------- ---------------
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
!!$  end subroutine noRI_
!!$  !*************************************************************


  !--------------- End of module ----------------------------------
end module noRI_module

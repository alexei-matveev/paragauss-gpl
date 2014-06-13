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
module reltrafo
  !---------------------------------------------------------------
  !
  !  Purpose: This module implement relativistic SO-DKH trafo, for the
  !  SR counterpart see relgrads.f90
  !
  !  Reads from disk:
  !
  !     START_DIR/BASENAME_real.datN
  !     START_DIR/BASENAME_imag.datN
  !
  !  with N  being the decimal  integer in the range  (1:n_irreps) and
  !  BASENAME  being  one  of   the  string  parameters  in  following
  !  invokations:
  !
  !     call get_hmatrix('overlap_rel', N, ...)
  !     call get_hmatrix('kin', N, ...)
  !     call get_hmatrix('nuc', N, ...)
  !     call get_hmatrix('pvsp', N, ...)
  !     call get_hmatrix('pvxp', N, ...)
  !
  !  See back_trafo_tapes.f90.  The input  data seems to  be available
  !  only  on  the host  running  the  master  process. The  procedure
  !  get_hmatrix()  takes  care  of  the  broadcasing  the  result  to
  !  everyone else.
  !
  !  Writes to disk (all irreps in one file):
  !
  !     START_DIR/overlap_real.dat
  !     START_DIR/overlap_imag.dat
  !     START_DIR/ham_kin_nuc_real.dat
  !     START_DIR/ham_kin_nuc_imag.dat
  !
  !  These  are  relativistically  transformed,  contracted  matrices.
  !  Moreover for legacy reasons kin- and nuc- matrices go to a single
  !  file.  Additional  output such as  intermediate matrices required
  !  by other  parts of the code  are computed and written  to disk by
  !  all workers.
  !
  !  FIXME: prefere integrals_on_file = .false. branch
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
  ! Author: Sonjoy Majumder
  ! Date:   23/03/03
  ! Description: BSOA
  !
  !----------------------------------------------------------------
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind  ! type specification parameters
  use symmetry_data_module, only: ssym
  use matrix_module, only: cmatrix
  use dimensions, only: dim_irr => IrrUBasDimSpor, &
                      & dim_irr_c => IrrBasDimSpor
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Interface statements ------------------------------

  interface p2_diag
     module procedure p2_diag_plain
     module procedure p2_diag_typed
  end interface

  !------------ public functions and subroutines ------------------

  public :: rel_trafo ! only for SPOR
  public :: p2_diag ! called from hfc_module

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of constants and variables ----
  logical,parameter  :: debug=.true.
!!$  real(RK)           :: speed_of_light = 137.03604_RK
  real(RK),parameter :: au2ev          = 27.211658_RK
  real(RK),parameter ::&
       & zero = 0.0_RK,&
       & half = 0.5_RK,&
       & one  = 1.0_RK,&
       & two  = 2.0_RK,&
       & eight= 8.0_RK,&
       & sixteen = 16.0_RK
  real(RK),parameter ::&
       & EPS_vc = 1.0E-04_RK
  integer(IK),parameter ::&
       & on  = 1,&
       & off = 0

  integer(IK),parameter            :: NO_IRREP = -1

  integer(IK) ::&
       & io_overlap_real,&
       & io_overlap_imag,&
       & io_ham_kin_nuc_real,&
       & io_ham_kin_nuc_imag

  type(cmatrix), allocatable ::    V_coul(:) ! (n_irr)
  type(cmatrix), allocatable :: PVYP_coul(:) ! (n_irr)

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  !*************************************************************
  subroutine rel_trafo()
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    use comm_module
    implicit none
    !------------ Declaration of formal parameters ---------------
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------

    logical :: i_am_master

    DPRINT 'rel_trafo: entered'

    if(comm_parallel())then
       i_am_master = comm_i_am_master()
    else
       i_am_master = .true.
    endif

    call init(i_am_master)
    call main()
    call done(i_am_master)

  end subroutine rel_trafo
  !*************************************************************

  !*************************************************************
  subroutine init(i_am_master)
    !  Purpose: ..
    use error_module
    use matrix_module, only: mm_init=>init
    implicit none
    logical,intent(in) :: i_am_master
    !** End of interface ***
#ifdef FPP_DEBUG
    integer(IK) :: i
#endif

    DPRINT MyID,'rt/init: entered'

    if(i_am_master)then

       call mm_init(verb=.false.)

!!$       if(debug)then
#ifdef FPP_DEBUG
          do i=1,size(dim_irr)
             write(*,'(A12,"(",I2,") u_dim ",I4," c_dim ",I4)')&
                  & ssym%name_proj(i), i, dim_irr(i),dim_irr_c(i)
          enddo
#endif
!!$       endif
    endif! i_am_master

    call open_output_tapes(i_am_master=i_am_master)
  end subroutine init

  subroutine done(i_am_master)
    use error_module
    implicit none
    logical,intent(in) :: i_am_master
    ! *** end of interface **

    call close_output_tapes(i_am_master=i_am_master)
  end subroutine done

  !*************************************************************

  subroutine main()
    !
    ! Driver  for   spin-orbit  relativistic  transformation   of  the
    ! hamiltonian. Executed in parallel context.
    !
    use error_module
    use spin_orbit_module, only: &
         is_on, whatis, &
         op_FitTrafo, op_BackTrafo
    implicit none
    ! *** end of interface ***

    integer :: irr, n_irr
    integer :: memstat

    n_irr = size(dim_irr)

    ! RelTrafo in Vnuc+Vcoul:
    if(whatis(op_BackTrafo).eq.2)then
       allocate(V_coul(n_irr), PVYP_coul(n_irr), stat=memstat)
       ASSERT(memstat==0)
       call build_coul(V_coul, PVYP_coul) ! allocates storage
    endif

    do irr = 1, size(dim_irr) ! n_irr

       DPRINT MyID,'rt/main: ...'
       DPRINT MyID,'rt/main: processing irrep ', irr, dim_irr(irr), dim_irr_c(irr)

       ! Pass  irrep, uncontracted  and  contracted matrix  dimensions
       ! e.g. for automatic arrays:
       call do_one_block(irr, dim_irr(irr), dim_irr_c(irr))
    enddo

    ! RelTrafo in Vnuc+Vcoul:
    if(whatis(op_BackTrafo).eq.2)then
       ! deallocates recursively:
       deallocate(V_coul, PVYP_coul, stat=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine main

  subroutine do_one_block(irr, n, n_c)
    !
    ! Executed in parallel context.
    !
    use comm, only: comm_rank
    use matrix_module, only: tr, sim, alloc, cmatrix, rdmatrix, &
        operator(*),assignment(=),operator(-),operator(+)
    use error_module
    use spin_orbit_module, only: is_on, whatis, &
         op_SpinOrbit, op_FitTrafo, op_NoPVxP, op_BSOA, &
         op_BackTrafo
    use back_trafo_tapes, only: get_hmatrix, put_matrix
    use options_module, only:options_kinematic_factors
    use operations_module,only :operations_gtensor, operations_hfc
    implicit none
    integer(IK), intent(in) :: irr,n,n_c
    !*** end of interface ***

    integer :: rank
    type(cmatrix)  :: UF,UB,Nuc,PVYP,V_rel
    type(cmatrix)  :: VP ! Nuc in momentum space
    type(cmatrix)  :: S,T,PVSP,PVXP
    type(cmatrix)  :: S_c,T_c,V_rel_c !<<< contracted
    type(rdmatrix) :: t_diag,Tp

    FPP_TIMER_DECL(ge)

    rank = comm_rank()

    DPRINT MyID,'rt/do_one_block: overlap contraction:'

    ! get overlap:
    call alloc(n, S)
    call get_hmatrix('overlap_rel', irr, S)

    ! --- contract overlap S:
    call alloc(n_c, S_c)
    call contract_cm(irr, S, S_c)

    ! get kin-energy:
    call alloc(n, T)
    call get_hmatrix('kin', irr, T)

    print *,'call momBas(n,T,S,t_diag,UF,UB)'
    FPP_TIMER_START(ge)
    !
    ! Compute t_diag, UF, UB
    ! (forward and backwards trafo matrices)
    ! using gneralized eigensolver:
    !
    call momBas(n, T, S, t_diag, UF, UB)

    ! p**2 = 2 * T_kin:
    t_diag%d = TWO * t_diag%d
    ! FIXME: is t_diag = 2 * t_diag%d supposed to work?

    FPP_TIMER_STOP(ge)
    print *,'done momBas(...)'

    DPRINT 'generalized: cpu=',FPP_TIMER_VALUE(ge)

    ! <<< UF, UB - forward and backward to p-space
    !     t_diag - p2
    !---------------------

    if (is_on(op_FitTrafo)) then
       DPRINT MyID,'rt/do_one_block: store transformation on files:'
       call put_matrix('p2diag', irr, t_diag)
       call put_matrix('forward', irr, UF)
       call put_matrix('backward', irr, UB)
    endif

    !DG
    if (rank == 0) then ! FIXME: historical assymetry
       if( operations_gtensor .or. operations_hfc )then

          ! write to file overlap matrix
          call save_matrix(S_c%re,irr,"gt_overlap_real")
          call save_matrix(S_c%im,irr,"gt_overlap_imag")

          if( options_kinematic_factors )then
             call save_vector(t_diag%d,irr,"gt_t_diag" )

             !write to file forward matrix
             call save_matrix(UF%re,irr,"forward_real" )
             call save_matrix(UF%im,irr,"forward_imag" )

             !write to file backward matrix
             call save_matrix(UB%re,irr,"U_backward_real" )
             call save_matrix(UB%im,irr,"U_backward_imag" )
          end if
       end if
    endif
    !DG

    DPRINT MyID,'rt/do_one_block: to momentum space:'
    !------------------------
    ! transformation of matricies to the momentum space >>>
    !
    call alloc(n, Nuc)
    call get_hmatrix('nuc', irr, Nuc)

#if FPP_AIX_XLF
    Nuc%re = - Nuc%re
    Nuc%im = - Nuc%im
#else
    Nuc = - Nuc  !<<< V := -V
#endif

    ! RelTrafo in Vnuc+Vcoul:
    if(whatis(op_BackTrafo).eq.2)then
       print *,'ADD V_COUL to V_IN(',irr,')'
       Nuc = Nuc + V_coul(irr)
    endif

    VP = sim(Nuc, UF) ! tr(UF) * Nuc * UF

    if (is_on(op_FitTrafo)) then
       DPRINT MyID,'rt/do_one_block: store Nuc_p on file:'
       call put_matrix('pnuc', irr, VP)
    endif

    call alloc(n, PVSP, PVXP)
    call get_hmatrix('pvsp', irr, PVSP)
    call get_hmatrix('pvxp', irr, PVXP)

    if(is_on(op_BSOA))then
       call bsoa(irr, PVXP)
    endif

    if(.not.is_on(op_NoPVxP))then
#if FPP_AIX_XLF
       call alloc(n, PVYP)
       PVYP%re = - PVSP%re - PVXP%re
       PVYP%im = - PVSP%im - PVXP%im
#else
       PVYP = - PVSP - PVXP  !<<< V := -V
#endif
    else
#if FPP_AIX_XLF
       call alloc(n, PVYP)
       PVYP%re = - PVSP%re
       PVYP%im = - PVSP%im
#else
       PVYP = - PVSP
#endif
    endif

    ! RelTrafo in Vnuc+Vcoul:
    if(whatis(op_BackTrafo).eq.2)then
       print *,'ADD V_PVYP to V_IN(',irr,')'
       PVYP = PVYP + PVYP_coul(irr)
    endif

    PVYP = sim(PVYP, UF)

    if (rank == 0) then ! FIXME: historical assymetry
       if( operations_gtensor .or. operations_hfc )then
          if (options_kinematic_factors) then
             call save_matrix(VP%re  ,irr,"gt_Nuc_real")
             call save_matrix(VP%im  ,irr,"gt_Nuc_imag")
             call save_matrix(PVYP%re,irr,"gt_PVYP_real")
             call save_matrix(PVYP%im,irr,"gt_PVYP_imag")
          end if
       end if
    endif

    !
    ! <<< stop transformations here
    !-------------------------

    !
    ! Compute Tp, V_rel:
    !
    select case(whatis(op_SpinOrbit))
    case (1)
       print *, MyID,'rt/do_one_block: call dkh1(t_diag,VP,PVYP,V_rel)'
       call dkh1(t_diag, VP, PVYP, Tp, V_rel)
    case (2)
       print *, MyID,'rt/do_one_block: call dkh2(t_diag,VP,PVYP,V_rel)'
       call dkh2(t_diag, VP, PVYP, Tp, V_rel)
    case (3)
       print *, MyID,'rt/do_one_block: call dkh3(t_diag,VP,PVYP,V_rel)'
       call dkh3(t_diag, VP, PVYP, Tp, V_rel)
    case default
       print *,'SpinOrbit=', whatis(op_SpinOrbit)
       ABORT('no such SpinOrbit')
    end select

    ! this is the differential effect of relativity:
    V_rel = V_rel - VP
    Tp%d  = Tp%d  - half * t_diag%d

    DPRINT MyID,'rt/do_one_block: back to real space:'
    ! back to real space:
    ! Nuc still holds the untransformed potential in real space,
    ! T   still holds the untransformed kin. eny. in real space:
    V_rel = tr(UB) * V_rel * UB + Nuc
    T     = tr(UB) * Tp    * UB + T

!!$    ! FIXME: coulomb
!!$    print *,'ADD V_COUL to V_REL(',irr,')'
!!$    V_rel = V_rel + V_coul(irr)

    DPRINT MyID,'rt/do_one_block: contract matrices:'
    ! contract matrices:
    call alloc(n_c, T_c, V_rel_c)
    call contract_cm(irr, T, T_c)
    call contract_cm(irr, V_rel, V_rel_c)

    !
    ! FIXME:  master/slave assymetry here.  Slaves only  save overlap,
    ! while master also saves kin and nuc:
    !
    if (rank == 0) then
       DPRINT MyID,'rt/do_one_block: store results:'
       ! store hamiltonian
       call write_rel_ham(S=S_c, T=T_c, V=V_rel_c)
    else
       call write_rel_ham(S=S_c)
    endif

    DPRINT MyID,'rt/do_one_block: exit'
  contains

    subroutine save_matrix(A, i_ir, filename)
      !  Purpose: writes matrix to a file
      !----------------------------------------------------------------
      use quadrupel_fname, only: qfilename
      use io, only: write_buffer
      implicit none
      integer(IK), intent(in) :: i_ir   ! irrep
      real(RK)   , intent(in) :: A(:,:) ! Matrix
      character(len=*) :: filename
      ! *** end of interface ***

      call write_buffer(qfilename(filename, i_ir, "dat"), A)
    end subroutine save_matrix

    subroutine save_vector(A, i_ir, filename)
      !  Purpose: writes matrix to a file
      !----------------------------------------------------------------
      use quadrupel_fname, only: qfilename
      use io, only: write_buffer
      implicit none
      integer(IK), intent(in) :: i_ir   ! irrep
      real(RK)   , intent(in) :: A(:)   ! Vector
      character(len=*) :: filename
      ! *** end of interface ***

      call write_buffer(qfilename(filename, i_ir, "dat"), A)
    end subroutine save_vector

  end subroutine do_one_block

#if 1
  subroutine momBas(n, T, S, t_diag, UF, UB)
    use matrix_module
    implicit none
    integer(IK), intent(in) :: n
    type(cmatrix), intent(in) :: T, S
    type(rdmatrix), intent(out) :: t_diag
    type(cmatrix), intent(out) :: UF, UB
    ! *** end of interface ***

    DPRINT 'momBas: n=',n

    call geigs(T, S, t_diag, UF)
    UB = tr(UF) * S
    !
    ! <<< UF, UB - forward and backward to p-space
    !     t_diag - p2
    DPRINT 'momBas: exit'
  end subroutine momBas
#else
  subroutine momBas(n, T, S, t_diag, UF, UB)
    use matrix_module
    implicit none
    integer(IK), intent(in) :: n
    type(cmatrix), intent(in) :: T, S
    type(rdmatrix), intent(out) :: t_diag
    type(cmatrix) , intent(out) :: UF, UB
    ! *** end of interface ***

    type(rdmatrix) :: s_diag
    type(cmatrix)  :: U, V, QF, QB, TF

    DPRINT 'momBas: n=',n

    !
    ! first canonical orthogonalization >>>
    !
    !    diagonalize overlap matrix >>>
    !
    call eigs(S, s_diag, U)

    s_diag%d = sqrt(s_diag%d)

    QB = s_diag * tr(U)
    QF = U * s_diag**(-1)
    !
    ! <<< done orthogonalization, QF,QB - forward and backward
    !     orthogonality transformations
    !

    !
    ! now prepare transformations to momentum space >>>
    !
    ! first diagonalize p**2
    !

    TF = sim(T, QF) ! <<< similarity transformation Q(+) * T * Q

    call eigs(TF, t_diag, V) ! <<< eigenvalue problem

    !
    ! now determine forward and backward transformation
    ! matrices UF and UB to the momentum space >>>
    !

    UF = QF * V
    UB = tr(V) * QB

    !
    ! <<< UF, UB - forward and backward to p-space
    !     t_diag - p2
    !
    DPRINT 'momBas: exit'
  end subroutine momBas
#endif

  subroutine dkh1(p2, V, PVYP, Tp, V_rel)
    use matrix_module
    implicit none
    type(rdmatrix),intent(in) :: p2
    type(cmatrix), intent(in) :: V, PVYP
    type(rdmatrix), intent(out) :: Tp
    type(cmatrix), intent(out) :: V_rel
    ! *** end of interface ***

    type(rdmatrix) :: Ep, Ap, Kp, ApKp, K2p2

    call p2_diag(p2, Tp, Ep, Ap, Kp, ApKp, K2p2)

    V_rel = mult(Ap, V, Ap) + mult(ApKp, PVYP, ApKp)
  end subroutine dkh1

  subroutine dkh2(p2, V, PVYP, Tp, V_rel)
    use matrix_module
    implicit none
    type(rdmatrix),intent(in) :: p2
    type(cmatrix), intent(in) :: V, PVYP
    type(rdmatrix), intent(out) :: Tp
    type(cmatrix), intent(out) :: V_rel
    ! *** end of interface ***

    integer(IK)   :: n
    type(rdmatrix) :: Ep, Ap, Kp, ApKp, K2p2

    type(cmatrix)  :: AVA, ARVRA
    real(RK)       :: DeltaE
    integer(IK)    :: i, j

    type(cmatrix)  :: RW, W2
    type(rdmatrix) :: ERm2

    n = size(p2%d)

    !
    ! now do relativistic corrections ... >>>
    !
    call p2_diag(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)

    !
    ! now construct ARVRA and AVA
    !
    AVA   = mult(Ap,V,Ap)

    ARVRA = mult(ApKp,PVYP,ApKp)

    !
    ! now construct first step of relativistic potential
    !

    V_rel = AVA + ARVRA

    !
    ! second step for constructing
    ! relativistic potential >>>
    !
    ! first  V -> V_tilda == Vpp`/(Ep + Ep`)
    !
    do j=1,n
       do i=1,n
          DeltaE = Ep%d(i)+Ep%d(j)
          AVA%re(i,j)   = AVA%re(i,j)/DeltaE
          AVA%im(i,j)   = AVA%im(i,j)/DeltaE

          ARVRA%re(i,j) = ARVRA%re(i,j)/DeltaE
          ARVRA%im(i,j) = ARVRA%im(i,j)/DeltaE
       enddo
    enddo

    ! now the sum:
    !
    ! - 1/2*{( W * 2Ep * W) - (Ep * W^2) - (W^2 * Ep)}
    !
    ! insert between every two W factors identity: 1 = (R * R) / K2p2
    !

    RW = K2p2 * AVA - ARVRA

    ERm2 = Ep * K2p2**(-1)

    V_rel = V_rel + tr(RW) * ERm2 * RW  !<<< sign

    W2 = - tr(RW) * K2p2**(-1) * RW  !<<< sign: tr(R*W) = -W*R

    V_rel = V_rel - half * ( Ep * W2 + W2 * Ep )
  end subroutine dkh2

  subroutine dkh3(p2, V, PVYP, Tp, V_rel)
    use matrix_module
    implicit none
    type(rdmatrix),intent(in) :: p2
    type(cmatrix), intent(in) :: V, PVYP
    type(rdmatrix), intent(out) :: Tp
    type(cmatrix), intent(out) :: V_rel
    ! *** end of interface ***

    integer(IK)   :: n
    type(rdmatrix) :: Ep,Ap,Kp,ApKp,K2p2

    type(cmatrix)  :: AVA,ARVRA
    type(cmatrix)  :: E1,E3
    real(RK)       :: DeltaE
    integer(IK)    :: i,j

    type(cmatrix)  :: RW,W2 !,wrk,H
    type(rdmatrix) :: ERm2

    n = size(p2%d)

    !
    ! now do relativistic corrections ... >>>
    !
    call p2_diag(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)

    !
    ! now construct ARVRA and AVA
    !
    AVA   = mult(Ap,V,Ap)

    ARVRA = mult(ApKp,PVYP,ApKp)

    !
    ! now construct first step of relativistic potential
    !

    E1 = AVA + ARVRA
    E3 = AVA + mult(K2p2**(-1), ARVRA, K2p2**(-1))

    !
    ! second step for constructing
    ! relativistic potential >>>
    !
    ! first  V -> V_tilda == Vpp`/(Ep + Ep`)
    !
    do j=1,n
       do i=1,n
          DeltaE = Ep%d(i)+Ep%d(j)
          AVA%re(i,j)   = AVA%re(i,j)/DeltaE
          AVA%im(i,j)   = AVA%im(i,j)/DeltaE

          ARVRA%re(i,j) = ARVRA%re(i,j)/DeltaE
          ARVRA%im(i,j) = ARVRA%im(i,j)/DeltaE
       enddo
    enddo

    ! now the sum:
    !
    ! - 1/2*{( W * 2Ep * W) - (Ep * W^2) - (W^2 * Ep)}
    !
    ! insert between every two W factors identity: 1 = (R * R) / K2p2
    !
    RW = K2p2 * AVA - ARVRA

    W2 = - tr(RW) * K2p2**(-1) * RW  !<<< sign: tr(R*W) = -W*R

    ERm2 = Ep * K2p2**(-1)

    V_rel = E1    + tr(RW) * ERm2 * RW  !<<< -W*E*W
    V_rel = V_rel - half * ( Ep * W2 + W2 * Ep )

    V_rel = V_rel + tr(RW) *  E3  * RW  !<<< -W*E1*W
    V_rel = V_rel + half * ( E1 * W2 + W2 * E1 )
  end subroutine dkh3

  subroutine p2_diag_plain(p2, Tp, Ep, Ap, Kp, ApKp, K2p2)
    use spin_orbit_module, only: speed_of_light
    implicit none
    real(RK), dimension(:) :: p2,Tp,Ep,Ap,Kp,ApKp,K2p2
    intent(in)                p2
    intent(out)                  Tp,Ep,Ap,Kp,ApKp,K2p2
    ! *** end of intarface ***

    real(RK) :: c,c2,c4

    c=speed_of_light; c2=c**2; c4=c2**2
    where ((p2/c2)<EPS_vc)
       Tp = c2*(p2/(two*c2)&
            &   - (p2/c2)**2/eight&
            &   + (p2/c2)**3/sixteen&
            &  )
       Ep = Tp + c2
    elsewhere
       Ep = sqrt(c4+c2*p2)
       Tp = Ep - c2
    end where
    Ap = sqrt((Ep+c2)/(2*Ep))
    Kp = c/(Ep+c2)
    ApKp = Ap*Kp
    K2P2 = p2*(Kp**2)
  end subroutine p2_diag_plain

  subroutine p2_diag_typed(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)
    use matrix_module, only: rdmatrix, alloc
    implicit none
    type(rdmatrix) :: p2,Tp,Ep,Ap,Kp,ApKp,K2p2
    intent(in)        p2
    intent(out)          Tp,Ep,Ap,Kp,ApKp,K2p2
    ! *** end of intarface ***

    call alloc(size(p2%d), Tp, Ep, Ap, Kp, ApKp, K2p2)

    call p2_diag(p2%d, Tp%d, Ep%d, Ap%d, Kp%d, ApKp%d, K2p2%d)
  end subroutine p2_diag_typed

  subroutine contract_cm(irr,um,cm)
    use matrix_module
    use contraction_module
    implicit none
    integer(IK),intent(in)      :: irr
    type(cmatrix),intent(in)    :: um
    type(cmatrix),intent(inout) :: cm
    ! *** end of interface ***

    if(irr.eq.NO_IRREP)return

    ASSERT(square(um))
    ASSERT(square(cm))

    call contract(irr, um%re, cm%re)
    call contract(irr, um%im, cm%im)
  end subroutine contract_cm

  subroutine write_rel_ham(S,T,V)
    use matrix_module
    use options_module, only      : options_integrals_on_file
    use integralstore_module, only: integralstore_2cob_ol                      &
                                  , integralstore_2cob_kin                     &
                                  , integralstore_2cob_nuc
    implicit none
    type(cmatrix),intent(in),optional :: S,T,V
    ! *** end of interface ***

    integer(IK) :: sz,counter,n,m

    sz = -1
    if ( present(S) ) then
        ASSERT(square(S))
        sz = size(S%re, 1)
    endif
    if ( present(T) ) then
        ASSERT(square(T))
        sz = size(T%re, 1)
    endif
    if(present(T).or.present(V))then
       ASSERT(present(T).and.present(V))
    endif
    ASSERT(sz.ne.-1)

    counter = 1

    if(options_integrals_on_file()) then
       ! overlap:
       do m=1,sz
          do n=1,m
             if(present(S))then
                ! overlap:
                write(io_overlap_real) S%re(m,n)
                write(io_overlap_imag) S%im(m,n)
             endif
             if(present(T))then
                ! kin-nuc, nuc with "-" historically:
                write(io_ham_kin_nuc_real) T%re(m,n), -V%re(m,n)
                write(io_ham_kin_nuc_imag) T%im(m,n), -V%im(m,n)
             endif
          end do
       end do
    else
       do m = 1,sz
          do n = 1,m
             ! overlap:
             if(present(S)) integralstore_2cob_ol(counter)= S%re(m,n)
             if(present(T))then
                ! kin-nuc:
                integralstore_2cob_kin(counter)=T%re(m,n)
                integralstore_2cob_nuc(counter)=V%re(m,n)
             endif
             counter=counter+1
             ! overlap:
             if(present(S)) integralstore_2cob_ol(counter)= S%im(m,n)
             if(present(T))then
                ! kin-nuc:
                integralstore_2cob_kin(counter)=T%im(m,n)
                integralstore_2cob_nuc(counter)=V%im(m,n)
             endif
             counter=counter+1
          enddo
       enddo
       if(counter>size(integralstore_2cob_kin)) counter=1
    endif
  end subroutine write_rel_ham

  subroutine open_output_tapes(i_am_master)
    use error_module
    use filename_module, only: tmpfile
    use options_module, only: options_integrals_on_file
    use iounitadmin_module
    implicit none
    logical,intent(in) :: i_am_master
    ! *** end of interface ***
    if (options_integrals_on_file()) then
       ! overlap:
       io_overlap_real=openget_iounit&
            (trim(tmpfile("overlap_real.dat")), form='unformatted',&
            status='unknown')
       io_overlap_imag=openget_iounit&
            (trim(tmpfile("overlap_imag.dat")), form='unformatted',&
            status='unknown')
       if(i_am_master)then
          ! kin-nuc:
          io_ham_kin_nuc_real = openget_iounit&
               (trim(tmpfile("ham_kin_nuc_real.dat")), form='unformatted',status='unknown')
          io_ham_kin_nuc_imag = openget_iounit&
               (trim(tmpfile("ham_kin_nuc_imag.dat")), form='unformatted',status='unknown')
       endif! i_am_master
    else
       call error("rt/open_output_tapes: not yet supported")
    end if
  end subroutine open_output_tapes

  subroutine close_output_tapes(i_am_master)
    use error_module
    use iounitadmin_module
    use options_module, only: options_integrals_on_file
    implicit none
    logical,intent(in) :: i_am_master
    ! *** end of interface ***
    if (options_integrals_on_file()) then
       ! overlap:
       call returnclose_iounit(io_overlap_real)
       call returnclose_iounit(io_overlap_imag)
       if(i_am_master)then
          ! kin-nuc:
          call returnclose_iounit(io_ham_kin_nuc_real)
          call returnclose_iounit(io_ham_kin_nuc_imag)
       endif! i_am_master
    else
       call error("rt/close_output_tapes: not yet supported")
    end if
  end subroutine close_output_tapes

!***********************************************
! Subroutine written by Sonjoy Majumder on 24/03/03

  subroutine QLi(irr,QL)
    use unique_atom_module
    use dimensions
    implicit none
    integer(IK), intent(in)             :: irr
    real(RK), dimension(:), intent(out) :: QL
    ! *** end of interface ***

    integer :: uaL,n,i_ua,L
    type(unique_atom_type), pointer :: ua
    real(RK) :: Zi, QLj
    integer(IK) :: n_rad,l_ua,i_rad

    integer(IK), parameter :: L_max = 10
    real(RK)               :: Z_scr(0:L_max)

    do L = 0, L_max
       n = ((L - 1) * L * (L + 1)) / 3 + (L * (L + 1)) / 2
       Z_scr(L) = 2 * n
    enddo
    DPRINT 'BSOA: magic numbers:', Z_scr

    uaL = 0
    do i_ua = 1, N_unique_atoms
       ua => unique_atoms(i_ua)
       Zi = ua%Z

       do L = 0, ua%lmax_ob

!!$          QLj=sqrt(ua%l_ob(L)%bsoa_par/Zi)
          ASSERT(L.LE.L_max)
          QLj = sqrt (Z_scr(L) / Zi)

          n_rad = ua%l_ob(L)%n_exponents
          l_ua = SymDimSpor(i_ua)%LM(L, irr)
          do i_rad = 1, n_rad * l_ua
             uaL = uaL + 1
             QL(uaL) = QLj
          enddo
       enddo
    enddo
    ASSERT(uaL.eq.size(QL))
  end subroutine QLi

  subroutine bsoa(irr,V)
    use matrix_module
    implicit none
    integer(IK), intent(in)      :: irr
    type(cmatrix), intent(inout) :: V
    ! *** end of interface ***

    type(rdmatrix) :: QL

    ASSERT(square(V))

    call alloc(size(V%re, 1), QL)

    call QLi(irr, QL%d)

    V = V - mult(QL, V, QL)
  end subroutine bsoa

  subroutine build_coul(V,VY)
    use filename_module, only: recover_dir
    use readwriteblocked_module, only: &
         readwriteblocked_startread, &
         readwriteblocked_read, &
         readwriteblocked_stopread, &
         readwriteblocked_tapehandle
    use fit_coeff_module, only: fit, get_fit ! dimenstions
    implicit none
    type(cmatrix), intent(out) :: V(:), VY(:) ! (n_irr)
    optional :: VY
    ! *** end of interface ***

    real(RK), allocatable :: ak(:) ! (nff)
    type(fit)             :: nfit

    character(len=14)     :: scf_file = 'saved_scfstate'
    type(readwriteblocked_tapehandle) :: th

    integer(IK) :: memstat

    print *,'rt::build_coul: entered'

    call get_fit(nfit)

    allocate(ak(nfit%n_ch),stat=memstat)
    ASSERT(memstat==0)

    call readwriteblocked_startread (&
         trim (recover_dir) // "/" // trim (scf_file) // ".dat", &
         th, variable_length=.true.)

    call readwriteblocked_read(ak,th)
    print *,'rt::build_coul: r COE=',sum(ak)

    call readwriteblocked_stopread(th)

    if(present(VY))then
       call do_build_coul(ak,V,VY)
    else
       call do_build_coul(ak,V)
    endif
  end subroutine build_coul

  subroutine do_build_coul(ak,V,VY)
    use matrix_module, only: cmatrix, alloc
    use fit_trafo_tapes, only: &
         init_fit_trafo_tapes, done_fit_trafo_tapes, &
         OpenTapes,CloseTapes, &
         TH_KIJ,TH_COUL,TH_SS,TH_SX,TH_RS,TH_RX
    use fit_coeff_module, only: &
         ff_map=>fit_coeff_ff_map
    implicit none
    real(RK), intent(in) :: ak(:) ! (nff) ! all fit func
    type(cmatrix), intent(out) :: V(:), VY(:) ! (n_irr)
    optional :: VY
    ! *** end of interface ***

    integer(IK) :: irr,n_ob
    integer(IK) :: nff,sff,rff
    real(RK), dimension(size(ak)) :: as ! (:sff) ! s-fitfunc
    real(RK), dimension(size(ak)) :: ar ! (:rff) ! r2-fitfunc

    print *,'rt::do_build_coul: entered'
    print *,'rt::do_build_coul: dims=', dim_irr

    nff = size(ak)
    sff = count(ff_map(:)%L.eq.-1) ! yes!
    rff = count(ff_map(:)%L.eq. 0) ! yes!
    print *,'rt::do_build_coul: nff=',nff,' sff=',sff,' rff=',rff

    as(:sff) = pack(ak,mask=(ff_map(:)%L.eq.-1))
    ar(:rff) = pack(ak,mask=(ff_map(:)%L.eq. 0))

    call init_fit_trafo_tapes()
    call OpenTapes(TH_KIJ)

    do irr = 1, size(dim_irr)
       n_ob = dim_irr(irr)
       print *,'rt::do_build_coul: irr=',irr,' dim=',n_ob

       ! Coulomb itself:
       call alloc(n_ob, V(irr))
       V(irr)%re = 0.0_rk
       V(irr)%im = 0.0_rk
       call do_build_F(n_ob, ak, TH_COUL, V(irr))

       if(present(VY))then
          ! Rel Ints of Coulomb:
          call alloc(n_ob, VY(irr))
          VY(irr)%re = 0.0_rk
          VY(irr)%im = 0.0_rk
          call do_build_F(n_ob, as(:sff), TH_SS, VY(irr))
          call do_build_F(n_ob, ar(:rff), TH_RS, VY(irr))
          call do_build_F(n_ob, as(:sff), TH_SX, VY(irr))
          call do_build_F(n_ob, ar(:rff), TH_RX, VY(irr))
       endif
    enddo

    print *,'call CloseTapes(TH_KIJ)'
    call CloseTapes(TH_KIJ)
     print *,'call done_fit_trafo_tapes()'
    call done_fit_trafo_tapes()
  end subroutine do_build_coul

  subroutine do_build_F(n_ob,ak,TH,V)
    use fit_trafo_tapes, only: ReadTape
    USE DEBUG, only: disp
    implicit none
    integer(IK), intent(in) :: n_ob
    real(RK)   , intent(in) :: ak(:) ! (nff)
    integer(IK), intent(in)      :: TH
    type(cmatrix), intent(inout) :: V
    ! *** end of interface ***

    real(RK), dimension(n_ob*(n_ob+1)/2,2) :: ch !_re, ch_im
    real(RK), dimension(size(ak),2) :: coul !_re, coul_im
    integer(IK) :: i_mn,m,n,k

    print *,'rt::do_build_F: (',TH, n_ob,')'

!!$    i_mn = 0
!!$    do m=1,n_ob
!!$       do n=1,m
!!$          i_mn = i_mn + 1
    do i_mn=1,n_ob*(n_ob+1)/2
          ! read in 3-Z integrals
          call ReadTape(TH, 1, coul(:, 1), coul(:, 2))
          ch(i_mn,1) = 0.0_rk
          ch(i_mn,2) = 0.0_rk
          do k=1,size(ak)
             ch(i_mn,1) = ch(i_mn,1) + ak(k) * coul(k,1)
             ch(i_mn,2) = ch(i_mn,2) + ak(k) * coul(k,2)
          enddo! k-loop
    enddo
!!$       enddo! n-loop
!!$    enddo!m-loop

    ! same code as in ham_calc_ch()
    i_mn = 1
    do m=1,n_ob
       do n=1,m-1
          V%re(n,m) = V%re(n,m) + ch(i_mn,1)
          V%re(m,n) = V%re(m,n) + ch(i_mn,1)
          ! this choice of sign is strange and must be considered furtheron !
          ! FIXME: sign changed
          V%im(n,m) = V%im(n,m) + ch(i_mn,2)
          V%im(m,n) = V%im(m,n) - ch(i_mn,2)
          i_mn = i_mn + 1
       enddo
       V%re(m,m)    = V%re(m,m) + ch(i_mn,1)
       V%im(m,m)    = V%im(m,m) + ch(i_mn,2)
       ASSERT(abs(V%im(m,m))<1.0E-7)
       i_mn = i_mn + 1
    enddo
    print *,'V%re=',maxval(abs(V%re))
    print *,'V%im=',maxval(abs(V%im))
!!$    call disp("IMAGE",V%im)
  end subroutine do_build_F

  !--------------- End of module ----------------------------------
end module reltrafo

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
module relgrads
  !---------------------------------------------------------------
  !
  !  Purpose: This module implement relativistic SR-DKH trafo, for the
  !  SO counterpart see reltrafo.f90
  !
  !  Reads from disk:
  !
  !     ... uncontracted kin- nuc- ovl- integrals and
  !     their first- and second derivatives.
  !     See realgrads_store.f90 for details.
  !
  !  Writes to disk:
  !
  !     START_DIR/overlap.dat
  !     START_DIR/ham_kin_nuc.dat
  !
  !  These are relativistically transformed and contracted matrices
  !  For legacy reasons kin- and nuc- go to a single file.
  !
  !  FIXME: prefer integrals_on_file = false branch
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
  ! Author:
  ! Date:
  ! Description:
  !
  !----------------------------------------------------------------
! define FPP_TIMERS 2
# include "def.h"
#ifdef WITH_MATRIX_PARALLEL
#define MATRIX_IMPL matrix_parallel
#else
#define MATRIX_IMPL matrix_module
#endif
  use type_module, only: IK => i4_kind, RK => r8_kind  ! type specification parameters
  use MATRIX_IMPL, only: rmatrix, rdmatrix
  use symmetry_data_module, only: ssym
  use error_module, only: MyID
  use dimensions, only: dim_irr => IrrUBasDim, &
                        dim_irr_c => IrrBasDim
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public :: rel_trafo_sr

  !================================================================
  ! End of public interface of module
  !================================================================

  interface bk
    module procedure bk_r_r
    module procedure bk_rd_r
    module procedure bk_r
  end interface

  !------------ Declaration of constants and variables ----

  integer(IK), parameter, public :: &
       IDKHN =     3, & ! (bitmask) DKH level 1, 2, or 3
       IENRG =     8, & ! need rel. Hamiltonian?
       IGRAD =    16, & ! need rel. Gradients?
       ISDER =    32    ! need rel. Second Derivatives?

  integer(IK) :: MODE = 2+IENRG ! second order DKH

!!$  real(RK)           :: speed_of_light = 137.03604_RK
  real(RK),parameter :: au2ev          = 27.211658_RK
  real(RK),parameter ::&
       & zero = 0.0_RK,&
       & half = 0.5_RK,&
       & one  = 1.0_RK,&
       & two  = 2.0_RK,&
       & eight= 8.0_RK,&
       & sixteen = 16.0_RK

  integer(IK),parameter            :: NO_IRREP = -1

  integer(IK) ::&
       & io_overlap,&
       & io_ham_kin_nuc,&
       & io_gt_overlap

  !
  ! We will store  O(N^2) matrices and O(N) diagonals  in memory. They
  ! are required again in e.g. gradient runs. IO is costly:
  !
  type(rmatrix), allocatable :: rmatrices(:, :) ! (4, n_irr)

  integer, parameter :: RM_FORWARD = 1
  integer, parameter :: RM_BACKWARD = 2
  integer, parameter :: RM_PNUC = 3
  integer, parameter :: RM_PPVSP = 4

  type(rdmatrix), allocatable :: dmatrices(:, :)  ! (1, n_irr)

  integer, parameter :: RD_P2DIAG = 1

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function is_on(OP) result(yes)
    implicit none
    integer(IK), intent(in) :: OP
    logical                 :: yes
    ! *** end of interface ***

    yes = ( IAND(MODE,OP) /= 0 )
  end function is_on

  function whatis(OP) result(and)
    implicit none
    integer(IK), intent(in) :: OP
    integer(IK)             :: and
    ! *** end of interface ***

    and = IAND(MODE,OP)
  end function whatis

  !*************************************************************
  subroutine rel_trafo_sr(IMODE, GRSTO, SDSTO)
    !  Purpose: ..
    !------------ Modules used ------------------- ---------------
    implicit none
    !------------ Declaration of formal parameters ---------------
    integer(IK), intent(in) :: IMODE
    real(RK), intent(inout) :: GRSTO(:) ! == gradient_totalsym in gradient_data_module
    real(RK), intent(inout) :: SDSTO(:,:) !  == dervs_totalsym in gradient_data_module
    optional :: GRSTO, SDSTO
    !** End of interface *****************************************
    !------------ Declaration of local variables -----------------

    integer(IK) :: jobs

    TRACE("rel_trafo_sr/entered")
    DPRINT MyID,'rt/rel_trafo_sr: entered, MODE=',IMODE
    MODE = IMODE

    if(is_on(IENRG)) call open_output_tapes()

    jobs = 0

    if ( is_on(IENRG) ) then
       ! if is_on(IGRAD) it will save some files for future use
       DPRINT MyID,'rt/rel_trafo_sr: 0) call ham_trafo()'
       call ham_trafo()
       jobs = jobs + 1
    endif

    if ( is_on(IGRAD) .and. .not. is_on(IENRG) ) then
       ! if is_on(ISDER) it will save grad mat into CPKS
       DPRINT MyID,'rt/rel_trafo_sr: 1) call gra_trafo()'
       ASSERT(present(GRSTO))
       call gra_trafo(GRSTO)
       jobs = jobs + 1
    endif

#ifdef WITH_SECDER
    if ( is_on(ISDER) ) then
       ASSERT(is_on(IGRAD))
       DPRINT MyID,'rt/rel_trafo_sr: 2) call sdr_trafo()'
       ASSERT(present(SDSTO))
       call sdr_trafo(SDSTO)
       jobs = jobs + 1
    endif
#endif

    if ( jobs == 0 ) then
       DPRINT 'rel_trafo_sr: MODE=',MODE
       ABORT('no such mode')
    endif

    if(is_on(IENRG)) call close_output_tapes()

    if (.not. is_on(IENRG)) then
       !
       ! This  is  a second  (post-scf)  call.  Gradients (and  second
       ! derivatives)  have  been  transformed.  It  is  now  safe  to
       ! dellocate global vars (see ham_trafo() for allocataion):
       !
       if (allocated(rmatrices)) then
          deallocate(rmatrices, dmatrices)
       endif
    endif

    DPRINT MyID,'rt/rel_trafo_sr: exit'
    TRACE("rel_trafo_sr/return")
  end subroutine rel_trafo_sr

  subroutine ham_trafo()
    !
    ! Driver   for   scalar   relativistic   transformation   of   the
    ! hamiltonian. Executed in parallel context.
    !
    implicit none
    ! *** end of interface ***

    integer(IK) :: irr
    intrinsic :: size

    if (is_on(IGRAD)) then
       !
       ! Allocate global vars for saving/restoring data across SCF:
       !
       allocate(rmatrices(4, size(dim_irr)), dmatrices(1, size(dim_irr)))
    endif

    do irr = 1, size(dim_irr) ! n_irr

       DPRINT MyID,'rt/main: ...'
       DPRINT MyID,'rt/main: processing irrep ', irr, dim_irr(irr), dim_irr_c(irr)

       ! Pass  irrep, uncontracted  and  contracted matrix  dimensions
       ! e.g. for automatic arrays:
       call do_one_block(irr, dim_irr(irr), dim_irr_c(irr))
    enddo
  end subroutine ham_trafo

  subroutine do_one_block(irr, n, n_c)
    !
    ! Executed in parallel context.
    !
    use comm, only: comm_rank
    use MATRIX_IMPL, only: tr, sim, rmatrix, rdmatrix, &
        operator(*), operator(**), operator(-), operator(+)
    use relgrads_store, only: rg_flag, RGPSEU
    implicit none
    integer(IK), intent(in) :: irr,n,n_c
    !*** end of interface ***

    integer :: rank
    type(rmatrix)  :: UF,UB,Nuc,V_rel
    type(rmatrix)  :: VP ! Nuc in momentum space
    type(rmatrix)  :: S,T,PVSP
    type(rdmatrix) :: t_diag,Tp

    FPP_TIMER_DECL(tot)
    FPP_TIMER_DECL(gev)
    FPP_TIMER_DECL(dk2)
    FPP_TIMER_DECL(dio)

    FPP_TIMER_START(tot)

    rank = comm_rank()

    FPP_TIMER_START(dio)
    !
    ! Allocate and get the data for all four matrices:
    !
    call read_relham(irr, n, S, T, Nuc, PVSP)
    FPP_TIMER_STOP(dio)

    TRACE("rel_trafo_sr/block_started")

    FPP_TIMER_START(gev)

    !
    ! Compute t_diag, UF, UB,
    ! use generalized eigensolver, when possible:
    !
    call momBas(T, S, t_diag, UF, UB)

    FPP_TIMER_STOP(gev)

    ! p**2 = 2 * T_kin:
    t_diag = TWO * t_diag

    ! <<< UF, UB - forward and backward to p-space
    !     t_diag - p2
    !---------------------

    DPRINT MyID,'rt/do_one_block: to momentum space:'
    !------------------------
    ! transformation of matricies to the momentum space >>>
    !

    ! FIXME: Nuc = - Nuc  !<<< V := -V
!!$    Nuc%m = - Nuc%m ! Regatta f90 workaround
    Nuc = (-1.0d0) * Nuc ! FIXME: unary minus?

    VP = sim(Nuc, UF) ! tr(UF) * Nuc * UF

    ! FIXME: PVSP = - PVSP
!!$    PVSP%m = - PVSP%m ! Regatta f90 workaround
    PVSP = (-1.0d0) * PVSP ! FIXME: unary minus?

    PVSP = sim(PVSP, UF)

    !
    ! <<< stop transformations here
    !-------------------------

    if (is_on(IGRAD)) then
       !
       ! Save  to  memory  for  the  future  use  in  gradients/second
       ! derivatives:
       !
       call put_1(RD_P2DIAG, irr, t_diag)
       call put_2(RM_FORWARD, irr, UF)
       call put_2(RM_BACKWARD, irr, UB)
       call put_2(RM_PNUC, irr, VP)
       call put_2(RM_PPVSP, irr, PVSP)
    endif

    FPP_TIMER_START(dk2)
    !
    ! Compute Tp, V_rel:
    !
    select case(whatis(IDKHN))
    case (1)
       DPRINT MyID,'rt/do_one_block: call dkh1(t_diag,VP,PVYP,V_rel)'
       call dkh1(t_diag, VP, PVSP, Tp, V_rel)
    case (2)
       DPRINT MyID,'rt/do_one_block: call dkh2(t_diag,VP,PVYP,V_rel)'
       call dkh2(t_diag, VP, PVSP, Tp, V_rel)
    case (3)
       DPRINT MyID,'rt/do_one_block: call dkh3(t_diag,VP,PVYP,V_rel)'
       call dkh3(t_diag, VP, PVSP, Tp, V_rel)
    case default
       print *,'DKH=', whatis(IDKHN)
       ABORT('no such DKH')
    end select
    FPP_TIMER_STOP(dk2)

    ! this is the differential effect of relativity:
    V_rel = V_rel - VP
    Tp = Tp  - half * t_diag

    DPRINT MyID,'rt/do_one_block: back to real space:'
    ! back to real space:
    ! Nuc still holds the untransformed potential in real space,
    ! T   still holds the untransformed kin. eny. in real space:
    V_rel = tr(UB) * V_rel * UB + Nuc
    T     = tr(UB) * Tp    * UB + T

    !
    ! Add the untransformed part of hamiltonian (pseudo):
    !
    if (rg_flag(RGPSEU)) then
       V_rel = V_rel - read_one(irr, n, RGPSEU) ! FIXME: sign
    endif

    TRACE("rel_trafo_sr/block_finished")

    FPP_TIMER_START(dio)
    DPRINT MyID,'rt/do_one_block: contract matrices and store results:'
    !
    ! Contract matrices and store contracted hamiltonian:
    !
    call write_rel_ham(irr, &
         contract(irr, n_c, S), &
         contract(irr, n_c, T), &
         contract(irr, n_c, V_rel))
    FPP_TIMER_STOP(dio)

    FPP_TIMER_STOP(tot)
#ifdef FPP_TIMERS
    print*,'RELTR: timing(',irr,')=',FPP_TIMER_VALUE(tot),'[',FPP_TIMER_SLICE(tot),']'
    print*,'RELTR: |-Mom. Bas.    :',FPP_TIMER_VALUE(gev),'[',FPP_TIMER_SLICE(gev),']'
    print*,'RELTR: |-DKH2         :',FPP_TIMER_VALUE(dk2),'[',FPP_TIMER_SLICE(dk2),']'
    print*,'RELTR: +-IO           :',FPP_TIMER_VALUE(dio),'[',FPP_TIMER_SLICE(dio),']'
#endif
    DPRINT MyID,'rt/do_one_block: exit'
  end subroutine do_one_block

  subroutine momBas(T, S, t_diag, UF, UB)
    use MATRIX_IMPL, only: rmatrix, rdmatrix, geigs, tr, &
         operator(*)
    implicit none
    type(rmatrix) , intent(in) :: T, S
    type(rdmatrix), intent(out) :: t_diag
    type(rmatrix) , intent(out) :: UF, UB
    ! *** end of interface ***

    DPRINT 'momBas2: entered'

    !
    ! Compute t_diag, UF:
    !
    call geigs(T, S, t_diag, UF)

    UB = tr(UF) * S
    !
    ! UF, UB - forward and backward to p-space
    ! t_diag - p2
    !
    DPRINT 'momBas2: exit'
  end subroutine momBas

  subroutine dkh1(p2, V, PVYP, Tp, V_rel)
    use MATRIX_IMPL
    implicit none
    type(rdmatrix),intent(in) :: p2
    type(rmatrix), intent(in) :: V, PVYP
    type(rdmatrix),intent(out) :: Tp
    type(rmatrix), intent(out) :: V_rel
    ! *** end of interface ***

    type(rdmatrix) :: Ep,Ap,Kp,ApKp,K2p2

    call p2_diag(p2, Tp, Ep, Ap, Kp, ApKp, K2p2)

    V_rel = mult(Ap, V, Ap) + mult(ApKp, PVYP, ApKp)
  end subroutine dkh1

  subroutine dkh2(p2, V, PVYP, Tp, V_rel)
    use MATRIX_IMPL
    implicit none
    type(rdmatrix),intent(in) :: p2
    type(rmatrix), intent(in) :: V, PVYP
    type(rdmatrix),intent(out) :: Tp
    type(rmatrix), intent(out) :: V_rel
    ! *** end of interface ***

    type(rdmatrix) :: Ep, Ap, Kp, ApKp, K2p2

    type(rmatrix) :: AVA, ARVRA

    type(rmatrix) :: RW, E2

    !
    ! now do relativistic corrections:
    !
    call p2_diag(p2, Tp, Ep, Ap, Kp, ApKp, K2p2)

    !
    ! now construct ARVRA and AVA
    !
    AVA   = mult(Ap, V, Ap)

    ARVRA = mult(ApKp, PVYP, ApKp)

    !
    ! now construct first step of relativistic potential
    !

    V_rel = AVA + ARVRA

    !
    ! second step for constructing
    ! relativistic potential:
    !

    ! now the sum:
    !
    ! - 1/2*{( W * 2Ep * W) - (Ep * W^2) - (W^2 * Ep)}
    !
    ! insert between every two W factors identity: 1 = (R * R) / K2p2
    !

    RW = K2p2 * AVA - ARVRA

    ! now  V -> V_tilda == Vpp`/(Ep + Ep`)
    RW = rpt(RW, Ep)

    ! FIXME: one may use RW as temp storage instead!
    E2 = tr( K2p2**(-1) * ( Ep * RW + RW * Ep ) ) * RW ! .. + h.c.)/2, see next:

    V_rel = V_rel + HALF * ( E2 + tr(E2) )
  end subroutine dkh2

  subroutine dkh3(p2, V, PVYP, Tp, V_rel)
    use MATRIX_IMPL
    implicit none
    type(rdmatrix),intent(in) :: p2
    type(rmatrix), intent(in) :: V, PVYP
    type(rdmatrix),intent(out) :: Tp
    type(rmatrix), intent(out) :: V_rel
    ! *** end of interface ***

    type(rdmatrix) :: Ep, Ap, Kp, ApKp, K2p2

    type(rmatrix)  :: AVA, ARVRA
    type(rmatrix)  :: E1, E3

    type(rmatrix)  :: RW, W2
    type(rdmatrix) :: ERm2

    !
    ! now do relativistic corrections:
    !
    call p2_diag(p2, Tp, Ep, Ap, Kp, ApKp, K2p2)

    !---------------------------------------
    ! now construct ARVRA and AVA
    !
    AVA   = mult(Ap, V, Ap)

    ARVRA = mult(ApKp, PVYP, ApKp)

    !
    ! now construct first step of relativistic potential
    !

    E1 = AVA + ARVRA
    E3 = AVA + mult(K2p2**(-1), ARVRA, K2p2**(-1))

    !
    ! second step for constructing
    ! relativistic potential >>>
    !

    ! now the sum:
    !
    ! - 1/2*{( W * 2Ep * W) - (Ep * W^2) - (W^2 * Ep)}
    !
    ! insert between every two W factors identity: 1 = (R * R) / K2p2
    !
    RW = K2p2 * AVA - ARVRA

    ! now  V -> V_tilda == Vpp`/(Ep + Ep`)
    RW = rpt(RW, Ep)

    W2 = - tr(RW) * K2p2**(-1) * RW  !<<< sign: tr(R*W) = -W*R

    ERm2 = Ep * K2p2**(-1)

    V_rel = E1    + tr(RW) * ERm2 * RW  !<<< -W*E*W
    V_rel = V_rel - half * ( Ep * W2 + W2 * Ep )

    V_rel = V_rel + tr(RW) *  E3  * RW  !<<< -W*E1*W
    V_rel = V_rel + half * ( E1 * W2 + W2 * E1 )
  end subroutine dkh3

  function RPT(V, e) result(VT)
    !
    ! Relativistic ``Perturbation Theory'' Weighting:
    !
    !    V(p,q) -> V(p,q)/ ( e(p) + e(q) )
    !
    use matrix_module, only: rmatrix, rdmatrix, size
    implicit none
    type(rmatrix), intent(in) :: V
    type(rdmatrix), intent(in) :: e
    type(rmatrix) :: VT
    ! *** end of interface ***

    integer(IK) :: p, q

    ! copy, instead of allocating:
    VT = V

    do q = 1, size(e)
       do p = 1, size(e)
          VT%m(p, q) = VT%m(p, q) / ( e%d(p) + e%d(q) )
       enddo
    enddo
  end function RPT

  function RPT_gr_weight(VT, e, gre) result(grVT)
    !
    ! Relativistic ``Perturbation Theory'' Weighting Gradients,
    !
    !    grVT(p,q) = - VT(p,q) * ( e'(p) + e'(q) ) / ( e(p) + e(q) )
    !
    ! WARNING 1 : gradients DUE TO DENOMINATOR ONLY,
    !             as if V were is a constant and we wanted to know
    !             gradietns of V-TILDE.
    ! WARNING 2 : expected input is a V-TILDE
    !             i.e. output VT of the previous sub:
    !                   VT = RPT(V).
    !
    ! This allows re-use of aVa/arVra storage:
    !      aVa(0) =  ...
    !      aVa(1) =  ...
    !         and now re-use the storage for tilde-counterparts:
    !      aVa(0) = RPT(aVa(0))
    !      aVa(1) = RPT(aVa(1)) + RPT_gr_weight(aVa(0))
    !                                                  ^^^^^^
    !
    use matrix_module, only: rmatrix, rdmatrix, size
    implicit none
    type(rmatrix), intent(in) :: VT
    type(rdmatrix), intent(in) :: e, gre
    type(rmatrix) :: grVT
    ! *** end of interface ***

    integer(IK) :: p, q

    ! copy, instead of allocating:
    grVT = VT

    do q = 1, size(e)
       do p = 1, size(e)
          grVT%m(p, q) = - grVT%m(p, q) * ( gre%d(p) + gre%d(q) ) / ( e%d(p) + e%d(q) ) ! NO **2
       enddo
    enddo
  end function RPT_gr_weight

  subroutine t_ders_s(nd,p2,T,E,A,K,AK,R)
    !
    ! kinematic factors f(p2) and their derivatives
    ! wrt p2:
    !
    ! T := E - c2
    ! R := K2 * p2 = T / ( E + c2 ) = T / ( 2*c2 + T )
    !
    use spin_orbit_module, only: c=>speed_of_light
    implicit none
    integer(IK), intent(in)  :: nd ! derivative order
    real(RK)                 :: p2
    real(RK),dimension(0:)   ::    T,E,A,K,AK,R ! (0:der)
    intent(in)                  p2
    intent(out)                    T,E,A,K,AK,R
    ! *** end of intarface ***

    integer(IK) :: d
    real(RK)    :: x
    real(RK), parameter :: EPS = 1.0E-04_rk

    ASSERT(nd>=0)

    x  = p2/c**2

    ! VALUES OF RELATIVISTIC KINEMATIC FACTORS:
    if( x > EPS ) then
       E(0) = sqrt( 1 + x )
       T(0) = E(0) - 1
       R(0) = T(0) / ( E(0) + 1 )
    else
       T(0) =       x    /   2 &
              -     x**2 /   8 &
              +     x**3 /  16   ! 0.6E-13 for x=1E-4
       E(0) = T(0) + 1
       R(0) =       x    /   4 &
              -     x**2 /   8 &
              + 5 * x**3 /  64   ! 0.8E-13 for x=1E-4
    endif
    A (0) = sqrt(( E(0) + 1 ) / ( 2 * E(0) ))
    K (0) =  1 / ( E(0) + 1 )
    AK(0) = A(0) * K(0)

    if(nd==0) goto 999

    ! FIRST DERIVATIVES:
    if( x > EPS ) then
       E(1) =  1 / ( 2 * E(0) )
       T(1) = E(1)
       R(1) = K(0)**2 / E(0)
    else
       ! FIXME: do I need expansion for O(1) ???
       T(1) =        ONE  /   2 &
              -      x    /   4 & ! FIXME: truncate as above?
              +  3 * x**2 /  16 & ! 0.2E-8  for x=1E-4
              -  5 * x**3 /  32   ! 0.2E-12
       E(1) = T(1)
       R(1) =        ONE  /   4 &
              -      x    /   4 &
              + 15 * x**2 /  64 &
              -  7 * x**3 /  32
    endif
    A (1) = - 1       / ( 8 * E(0)**3 * A(0) )
    K (1) = - K(0)**2 / ( 2 * E(0) )
    AK(1) =   A(1) * K(0) + A(0) * K(1)

    if(nd==1) goto 999

    ! SECOND DERIVATIVES:
    if( x > EPS ) then
       E(2) = - 1 / ( 4 * E(0)**3 )
       T(2) =   E(2)
       R(2) = - ( 3 - 2 * K(0) ) * K(0)**2 / ( 2 * E(0)**3 )
    else
       ! FIXME: do I need expansion for O(1) ???
       T(2) = -      ONE  /   4 &
              +  3 * x    /   8 &
              - 15 * x**2 /  32 &
              + 35 * x**3 /  64
       E(2) = T(2)
       R(2) = -      ONE  /   4 &
              + 15 * x    /  32 &
              - 21 * x**2 /  32 &
              +105 * x**3 / 128
    endif
    A (2) =  ( 6 -     K(0) )           / ( 32 * E(0)**5 * A(0) )
    K (2) =  ( 3 - 2 * K(0) ) * K(0)**2 / (  4 * E(0)**3 )
    AK(2) =  A(2) * K(0) + 2 * A(1) * K(1) + A(0) * K(2)

    ASSERT(nd<=2)

999 CONTINUE ! SCALE DIMENSIONLESS FACTORS
    ! x := p2/c2   =>    d/dp2 = (1/c2) d/dx
    do d=0,nd
       E (d) =  E(d) * c**2 / c**(2*d)
       T (d) =  T(d) * c**2 / c**(2*d)
       K (d) =  K(d) / c    / c**(2*d)
       AK(d) = AK(d) / c    / c**(2*d)
       A (d) =  A(d)        / c**(2*d)
       R (d) =  R(d)        / c**(2*d)
    enddo
  end subroutine t_ders_s

  subroutine t_ders(nd,p2,Tp,Ep,Ap,Kp,ApKp,K2p2)
    ! wrapper for the prev
    implicit none
    integer(IK), intent(in) :: nd ! derivative order
    real(RK),dimension(:)    :: p2
    real(RK),dimension(:,0:) ::    Tp,Ep,Ap,Kp,ApKp,K2p2 ! (:,0:nd+)
    intent(in)                 p2
    intent(out)                   Tp,Ep,Ap,Kp,ApKp,K2p2
    ! *** end of intarface ***

    integer(IK) :: i, sp2(1)
    intrinsic :: size
    sp2 = shape( p2 )
    do i=1,sp2(1) ! over p2 values
       call t_ders_s(nd,p2(i),Tp(i,:),Ep(i,:),Ap(i,:),Kp(i,:),ApKp(i,:),K2p2(i,:))
    enddo
  end subroutine t_ders

  subroutine p2_diag_ders(nd, p2, Tp, Ep, Ap, Kp, ApKp, K2p2)
    use matrix_module, only: rdmatrix, assignment(=)
    implicit none
    integer(IK), intent(in)       :: nd ! derivative order
    type(rdmatrix), dimension(0:) :: p2,Tp,Ep,Ap,Kp,ApKp,K2p2
    ! 0: value, 1: X-gradient, 2: Y-gradient, 3: XY-derivative.
    intent(in)                       p2
    intent(out)                         Tp,Ep,Ap,Kp,ApKp,K2p2
    ! *** end of intarface ***

    integer(IK), parameter :: X=1, Y=2, XY=3
    integer(IK)            :: tmpsize(1)
    real(RK), dimension(size(p2(0)%d), 0:nd) :: T,E,A,K,AK,R2
    intrinsic :: size

    ASSERT(nd>=0)
    ASSERT(nd<=2)

    ! compute (T, dT/dp2,..), (E, dE/dp2,..) at p2:
    call t_ders(nd, p2(0)%d, T, E, A, K, AK, R2)

    ! copy out-args (use of custom assignment here):
    Tp  (0) = T (:,0)
    Ep  (0) = E (:,0)
    Ap  (0) = A (:,0)
    Kp  (0) = K (:,0)
    ApKp(0) = AK(:,0)
    K2p2(0) = R2(:,0)
    if( nd == 0 ) RETURN

    tmpsize=shape(p2)
    ASSERT(tmpsize(1)>=2)
    tmpsize=shape(Tp)
    ASSERT(tmpsize(1)>=2)
    ASSERT(X==1)

    Tp  (X) = T (:,1) * p2(X)%d
    Ep  (X) = E (:,1) * p2(X)%d
    Ap  (X) = A (:,1) * p2(X)%d
    Kp  (X) = K (:,1) * p2(X)%d
    ApKp(X) = AK(:,1) * p2(X)%d
    K2p2(X) = R2(:,1) * p2(X)%d
    if( nd == 1 ) RETURN

    ASSERT(nd==2)
    tmpsize=shape(p2)
    ASSERT(tmpsize(1)==4)
    tmpsize=shape(Tp)
    ASSERT(tmpsize(1)==4)

    Tp  (Y) = T (:,1) * p2(Y)%d
    Ep  (Y) = E (:,1) * p2(Y)%d
    Ap  (Y) = A (:,1) * p2(Y)%d
    Kp  (Y) = K (:,1) * p2(Y)%d
    ApKp(Y) = AK(:,1) * p2(Y)%d
    K2p2(Y) = R2(:,1) * p2(Y)%d

    ! E(p)^xy = E` * p^xy + E`` * p^x * p^y:
    Tp  (XY) = T (:,1) * p2(XY)%d + T (:,2) * p2(X)%d * p2(Y)%d
    Ep  (XY) = E (:,1) * p2(XY)%d + E (:,2) * p2(X)%d * p2(Y)%d
    Ap  (XY) = A (:,1) * p2(XY)%d + A (:,2) * p2(X)%d * p2(Y)%d
    Kp  (XY) = K (:,1) * p2(XY)%d + K (:,2) * p2(X)%d * p2(Y)%d
    ApKp(XY) = AK(:,1) * p2(XY)%d + AK(:,2) * p2(X)%d * p2(Y)%d
    K2p2(XY) = R2(:,1) * p2(XY)%d + R2(:,2) * p2(X)%d * p2(Y)%d
  end subroutine p2_diag_ders

  subroutine p2_diag(p2, Tp, Ep, Ap, Kp, ApKp, K2p2)
    !
    ! Special case of the above, for values only.
    !
    use MATRIX_IMPL, only: rdmatrix, matrix, array, size
    implicit none
    type(rdmatrix) :: p2, Tp, Ep, Ap, Kp, ApKp, K2p2
    intent(in)     :: p2
    intent(out)    ::     Tp, Ep, Ap, Kp, ApKp, K2p2
    ! *** end of intarface ***

    real(RK), dimension(size(p2), 0:0) :: T, E, A, K, AK, R2

    ! Compute T, E, ... at p2:
    call t_ders(0, array(p2), T, E, A, K, AK, R2)

    ! Copy out-args:
    Tp    = matrix(T (:, 0))
    Ep    = matrix(E (:, 0))
    Ap    = matrix(A (:, 0))
    Kp    = matrix(K (:, 0))
    ApKp  = matrix(AK(:, 0))
    K2p2  = matrix(R2(:, 0))
  end subroutine p2_diag

  subroutine gra_trafo(grsto)
    use comm, only: comm_rank
    use gradient_data_module, only: gradient_rel_index!(gradient_data_n_gradients,n_irreps)
    implicit none
    real(RK), intent(inout) :: grsto(:) ! == gradient_totalsym in gradient_data_module
    ! *** end of interface ***

    integer(IK) :: irr,n_irr, n_u, n_c

    ! parallel context:

    n_irr = ssym%n_irrep
    DPRINT MyID,'rt/gra: n_irr=', n_irr, 'rank=', comm_rank()

    do irr = 1,n_irr
       ! number of gradients to transform:
       n_u  = dim_irr(irr)   ! uncontracted dim
       n_c  = dim_irr_c(irr) ! contracted dim

       DPRINT MyID,'rt/gra: processing irrep ',irr,dim_irr(irr),dim_irr_c(irr)
       DPRINT MyID,'rt/gra:    mask(',irr,')=',gradient_rel_index(:,irr)

       call do_gra_trafo(irr,n_u,n_c,grsto &
                        , mask = (gradient_rel_index(:, irr) == 1 + comm_rank()))
    enddo
  end subroutine gra_trafo

#ifdef WITH_SECDER
  subroutine sdr_trafo(sdsto)
    use comm, only: comm_rank
    use gradient_data_module, only: sdind => gradient_der_index!(n_gr, n_gr, n_irr)
    implicit none
    real(RK), intent(inout) :: sdsto(:,:) ! == dervs_totalsym in gradient_data_module
    ! *** end of interface ***

    integer(IK) :: irr,n_irr, n_u, n_c
    real(RK)    :: sd(size(sdsto,1),size(sdsto,2))

    integer(IK) :: n_gr, tmpsize(2)
    intrinsic :: size

    ! parallel context:

    tmpsize = shape( sdsto )
    n_gr    = tmpsize(1)
    ASSERT(n_gr==tmpsize(2))

    n_irr = ssym%n_irrep
    DPRINT MyID,'rt/sdr: n_irr=', n_irr, 'rank=', comm_rank()

    sd = zero

    do irr = 1,n_irr
       n_u  = dim_irr(irr)   ! uncontracted dim
       n_c  = dim_irr_c(irr) ! contracted dim

       DPRINT MyID,'rt/sdr: processing irrep ',irr,dim_irr(irr),dim_irr_c(irr)

       call do_sdr_trafo( irr, n_u, n_c, n_gr, sd          &
                        , mask=(sdind(:, :, irr) == 1 + comm_rank()))
       ! mask derivatives assigned to this processor!
    enddo

    sdsto = sdsto + sd
  end subroutine sdr_trafo
#endif

  subroutine do_gra_trafo(irr,n,n_c,grsto,mask)
    use matrix_module, only: tr, sim, alloc, rmatrix, rdmatrix, &
        operator(*), operator(**), assignment(=), operator(-), &
        operator(+), diag
    use relgrads_store, only: rg_flag, RGPSEU
#ifdef WITH_SECDER
    use cpks_utils, only: cpks_h1_store
#endif
    USE STRINGS, only: itoa
    implicit none
    integer(IK), intent(in)    :: irr,n,n_c
    real(RK)   , intent(inout) :: grsto(:) ! gradient storage
    logical    , intent(in)    :: mask(:)  ! (n_gr) mask for this proc
    !*** end of interface ***

    type(rmatrix)                  :: UF, UB
    type(rmatrix)                  :: P

    ! 0: value, 1: first derivative:
    type(rdmatrix), dimension(0:1) :: p2
    type(rmatrix) , dimension(0:1) :: S, T, V, O, VR ! O=PVSP
    type(rdmatrix), dimension(0:1) :: tp,e,a,k,ak,r2,m2
    type(rmatrix) , dimension(0:1) :: aVa, akOka
    type(rmatrix) , dimension(0:1) :: RW, XX, E2
    type(rmatrix) , dimension(0:1) :: VRC ! contracted
    type(rmatrix)                  :: DMAT
    type(rmatrix) , dimension(0:1) :: NR, NRP ! NR in real and mom space

    integer(IK) :: tmpsize(1), tmpsiz2(1)
    integer(IK) :: nd
    integer(IK) :: j, igr
    real(RK)    :: gradient

    FPP_TIMER_DECL(tot)
    FPP_TIMER_DECL(loo)
    FPP_TIMER_DECL(nd0)
    FPP_TIMER_DECL(gen)
    FPP_TIMER_DECL(frw)
    FPP_TIMER_DECL(dk1)
    FPP_TIMER_DECL(rwm)
    FPP_TIMER_DECL(dk2)
    FPP_TIMER_DECL(bck)
    FPP_TIMER_DECL(dio)

    FPP_TIMER_START(tot)

    DPRINT MyID,'do_gra_trafo:   irr=',irr,' dim=',n,' cdim=',n_c
    DPRINT MyID,'do_gra_trafo:  mask=',mask

    ! count gradients to be processed by me:
    igr = count(mask)
    if( igr == 0 )then
      WARN('ZERO grads case!')
      ! DO  NOT  RETURN  early,  get_2(), requires  communication  for
      ! distributed matrix impl.  FIXME: all the rest is  done in vain
      ! on this worker, of course.
    endif

    ! Prepare (contracted) density matrix (in orbital space).
    ! (density_data seems to be shut down)
    call alloc(n_c, DMAT)
    call densmat(irr, DMAT)

    !
    ! Recover matrices saved in do_one_block() called from ham_trafo()
    ! in pre-scf.  These are already in mom-space, get_2() may require
    ! communication between all processes in MPI world:
    !
    p2(0) = rdmatrix(get_1(n, RD_P2DIAG, irr))
    UF = rmatrix(get_2(n, RM_FORWARD, irr))
    UB = rmatrix(get_2(n, RM_BACKWARD, irr))
    V(0) = rmatrix(get_2(n, RM_PNUC, irr))
    O(0) = rmatrix(get_2(n, RM_PPVSP, irr))

    V(0) = - V(0) ! FIXME: sign
    O(0) = - O(0) ! FIXME: sign

    ! FIXME: maybe redundant allocations here?
    call alloc(n, V(1), O(1))

    ![[=== Rel-Effect-Only: ===
    ! Compute NR-base in mom space: NRP(0) = V(0) - half * p2(0)
    NRP(0) = V(0)
    tmpsize = shape( p2(0)%d )
    do j = 1, tmpsize(1)
      ! FIXME: sign:
      NRP(0)%m(j,j) = NRP(0)%m(j,j) - half * p2(0)%d(j)
    enddo
    !]]========================

    ! storage for gradient of overlap:
    call alloc(n, S(1)) ! S(0) is implicitly a unity

    ! storage for gradient of kin. energy:
    call alloc(n, T(1)) ! T(0) is implicitly a diagonal p2

    ! contracted:
    call alloc(n_c, VRC(1))

    ! The iteration 0 (zero) of the following loop is
    ! is used to PRECOMPUTE the zero-order matrices
    ! that will be used in further gradient calculations.
    ! The idea was to keep the energy/gradient code close
    ! enough for cross-cheking. Any better solutions?
    tmpsize = shape( mask )
    do igr=0, tmpsize(1)
                                       FPP_TIMER_START(loo)
       nd = 1              ! gradients
       if( igr == 0 ) nd = 0 ! energy
       DPRINT '>>>>>>>>>>>> i=',i,' ND=',nd,'>>>>>>>>>>>>>'

!!!GR: >>>>>>>>>>>>>>>>>> gradients >>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 1 )then
                                       FPP_TIMER_START(gen)
       ! process only gradients that were assigned to this processor
       if( .not. mask(igr) ) CYCLE ! to next gradient

                                       FPP_TIMER_START(dio)
       ! fetch gradients of overlap, kinetic, of V and PVSP :
       call read_grads(irr, igr, S(1), T(1), V(1), O(1))
                                       FPP_TIMER_STOP(dio)

       ![[=== Rel-Effect-Only: ===
       ! Compute NR-base in real space:
       NR(1) = V(1) + T(1) ! FIXME: sign!
       !]]========================

       ! to the momentum basis:
       ! FIXME: S(1) and T(1) are very sparse matrices.
       !        This can be used for speed in similarity trafo.
       ! FIXME: historically there is "a minus" between
       !        S(1), T(1) and gradients of S(0), T(0)
       S(1)  = - sim(S(1),UF)
       T(1)  = - sim(T(1),UF) * 2.0_rk ! since p^2 = 2 * T

       ! build the response generator, and response of eigenvalues:
       call PT(p2(0), T(1), S(1), p2(1), P)
                                       FPP_TIMER_STOP(gen)

                                       FPP_TIMER_START(frw)
       ! 1) DIRECT: to the momentum basis,
       ! 2) PULLAY: due to changes in momentum basis
       V(1) = sim(V(1), UF) + bk(V(0), P)
       O(1) = sim(O(1), UF) + bk(O(0), P)
                                       FPP_TIMER_STOP(frw)

       ![[=== Rel-Effect-Only: ===
       ! Compute NR-base in mom space: NRP(1) = V(1) - half * p2(1)
       NRP(1) = V(1)
       tmpsiz2 = shape( p2(1)%d )
       do j = 1, tmpsiz2(1)
         ! FIXME: sign:
         NRP(1)%m(j,j) = NRP(1)%m(j,j) - half * p2(1)%d(j)
       enddo
       !]]========================

       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!GR: >>>>>>>>>>>>>> energy and gradients >>>>>>>>>>>>>>>>

       ! compute values and derivative of kinematic factors:
       call p2_diag_ders(nd, p2, tp, e, a, k, ak, r2)

       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!EN: >>>>>>>>>>>>>>>>>>>> energy >>>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 0 )then
                                       FPP_TIMER_START(nd0)
       ! DKH1 level: V := aVa + akOka
       aVa(0)   =  a(0) * V(0) *  a(0)

       akOka(0) = ak(0) * O(0) * ak(0)

       VR(0)    = aVa(0) + akOka(0)
                                       FPP_TIMER_STOP(nd0)
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!GR: >>>>>>>>>>>>>>>>>> gradients >>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 1 )then
                                       FPP_TIMER_START(dk1)
       ! DKH1 level: V := aVa + akOka
       aVa(1)   =  a(0) * V(1) *  a(0) &
                +  a(1) * V(0) *  a(0) &
                +  a(0) * V(0) *  a(1)

       akOka(1) = ak(0) * O(1) * ak(0) &
                + ak(1) * O(0) * ak(0) &
                + ak(0) * O(0) * ak(1)

       VR(1)    = aVa(1) + akOka(1)
                                       FPP_TIMER_STOP(dk1)
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       ! DKH2 level:
       ! we are going to compute
       !   E2 = - WEW - {W^2,E}/2
       ! using the identity
       !   W^2 = WR * (R^-2) * RW = - (RW)^T * 1/(k2p2) * RW
       ! and finally:
       !  E2 = (1/2)[ ( R^-2 * ( E * RW + RW * E ) )^T * RW + h.c. ]
!!!EN: >>>>>>>>>>>>>>>>>>>> energy >>>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 0 )then
                                       FPP_TIMER_START(nd0)
       ! RW = r2 * aVa - arVra:
       RW(0) = r2(0) * aVa(0) - akOka(0)
       RW(0) = RPT(RW(0), e(0))
                                       FPP_TIMER_STOP(nd0)
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!GR: >>>>>>>>>>>>>>>>>> gradients >>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 1 )then
                                       FPP_TIMER_START(rwm)
       ! RW = r2 * aVa - arVra:
       RW(1) = r2(0) * aVa(1) - akOka(1) &
             + r2(1) * aVa(0)

       ! WARNING: RW(0) already redefined to RW-TILDE!
       !          RPT_gr_weight() expects that
       RW(1) =           RPT( RW(1)  ,e(0)      ) &
             + RPT_gr_weight( RW(0)  ,e(0), e(1))
                                       FPP_TIMER_STOP(rwm)
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!EN: >>>>>>>>>>>>>>>>>>>> energy >>>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 0 )then
                                       FPP_TIMER_START(nd0)
       ! r2: R^ 2 =   k2p2
       ! m2: R^-2 = 1/k2p2, read as ``R minus 2''
       m2(0) =            r2(0)**(-1)
!!!GR: m2(1) =  - r2(1) * r2(0)**(-2)

       ! anticommutator {e,RW}:
       XX(0) = RW(0) * e(0) + e(0) * RW(0)

       E2(0) = tr( m2(0) * XX(0) ) * RW(0) ! +h.c.)/2, see below:

       VR(0) = VR(0) - HALF * ( E2(0) + tr(E2(0)) )

       ! add diagonal kinetic energy:
       do j=1,n
          VR(0)%m(j,j) =  - tp(0)%d(j) + VR(0)%m(j,j)
       enddo

       ![[=== Rel-Effect-Only: ===
       ! replace rel-ham by rel-eff:
       VR(0) = VR(0) - NRP(0)
       !]]========================
                                       FPP_TIMER_STOP(nd0)
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!GR: >>>>>>>>>>>>>>>>>> gradients >>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 1 )then
                                       FPP_TIMER_START(dk2)
       ! r2: R^ 2 =   k2p2
       ! m2: R^-2 = 1/k2p2, read as ``R minus 2''
!!!EN: m2(0) =            r2(0)**(-1)
!      m2(1) =  - r2(1) * r2(0)**(-2)
       m2(1) =  - ( r2(1) * m2(0) ) * m2(0) ! FIXME: R^-4 may be too big a number!


       XX(1) = RW(1) * e(0) + e(0) * RW(1) + RW(0) * e(1) + e(1) * RW(0)

       E2(1) = tr( m2(1) * XX(0) ) * RW(0) &
             + tr( m2(0) * XX(1) ) * RW(0) &
             + tr( m2(0) * XX(0) ) * RW(1) ! +h.c.)/2, see below:

       VR(1) = VR(1) - HALF * ( E2(1) + tr(E2(1)) ) ! FIXME: sign?

       ! add diagonal kinetic energy:
       do j=1,n
          VR(1)%m(j,j) =  - tp(1)%d(j) + VR(1)%m(j,j) ! FIXME: sign?
       enddo
                                       FPP_TIMER_STOP(dk2)
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!!GR: >>>>>>>>>>>>>>>>>> gradients >>>>>>>>>>>>>>>>>>>>>>>
       if( nd == 1 )then
                                       FPP_TIMER_START(bck)

       ![[=== Rel-Effect-Only: ===
       ! compute the differential effect in mom space:
       ! eff = SR    - NR
       VR(1) = VR(1) - NRP(1)
       !]]========================

       ! PULLAY CORRECTIONS ON BACK-TARFO:
       ! Rel-Effect-Only: only that due to rel-chage in VR(0),
       !                  NR-base NRP(1) already removed above
       VR(1) = VR(1) - bk(VR(0), P)

       ! BACK TO ORBITAL SPACE:
       VR(1) = sim(VR(1),UB)
!!!EN: VR(0) = sim(VR(0),UB) ! unneded!

       ![[=== Rel-Effect-Only: ===
       ! compute the ham matrix in real space:
       ! SR  = eff   + NR
       VR(1) = VR(1) + NR(1)
       !]]========================

       ![[=== add the untransformed part (pseudo): ===
       if( rg_flag(RGPSEU) )then
         call read_grads(irr, igr, NR(1))
         VR(1) = VR(1) + NR(1) ! FIXME: sign
       endif
       ![[============================================

       ! contract the basis:
       call contract_rm(irr, VR(1), VRC(1))

       gradient = trace2(VRC(1), DMAT)
       DPRINT 'GRADIENT(',irr,igr,')=',gradient

       ! add up this gradient:
       tmpsiz2 = shape( grsto )
       ASSERT(igr<=tmpsiz2(1))
       grsto(igr) = grsto(igr) + gradient
                                       FPP_TIMER_STOP(bck)
#ifdef WITH_SECDER
       if ( is_on(ISDER) ) then
         DPRINT  'RGR: call cpks_h1_store(',irr,igr,',VRC(1)%m,+1)'
         call cpks_h1_store(irr, igr, VRC(1)%m, +1)
       endif
#endif
       endif ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                       FPP_TIMER_STOP(loo)
    enddo

    DPRINT MyID,'do_gra_trafo: exit'
    FPP_TIMER_STOP(tot)
#ifdef FPP_TIMERS
    print*,'RELGR: timing(',irr,')=',FPP_TIMER_VALUE(tot),'[',FPP_TIMER_SLICE(tot),']'
    print*,'RELGR: loop           :',FPP_TIMER_VALUE(loo),'[',FPP_TIMER_SLICE(loo),']'
    print*,'RELGR: |-generator    :',FPP_TIMER_VALUE(gen),'[',FPP_TIMER_SLICE(gen),']'
    print*,'RELGR: |-forward      :',FPP_TIMER_VALUE(frw),'[',FPP_TIMER_SLICE(frw),']'
    print*,'RELGR: |-DKH1         :',FPP_TIMER_VALUE(dk1),'[',FPP_TIMER_SLICE(dk1),']'
    print*,'RELGR: |-RW matrix    :',FPP_TIMER_VALUE(rwm),'[',FPP_TIMER_SLICE(rwm),']'
    print*,'RELGR: |-DKH2         :',FPP_TIMER_VALUE(dk2),'[',FPP_TIMER_SLICE(dk2),']'
    print*,'RELGR: |-backward     :',FPP_TIMER_VALUE(bck),'[',FPP_TIMER_SLICE(bck),']'
    print*,'RELGR: +-zero-order   :',FPP_TIMER_VALUE(nd0),'[',FPP_TIMER_SLICE(nd0),']'
    print*,'RELGR: +-disk io      :',FPP_TIMER_VALUE(dio),'[',FPP_TIMER_SLICE(dio),']'
#endif
  end subroutine do_gra_trafo

#ifdef WITH_SECDER
  subroutine do_sdr_trafo(irr,n,n_c,n_gr,sdsto,mask)
    use matrix_module, only: tr, sim, alloc, rmatrix, rdmatrix, &
        operator(*), operator(**), assignment(=), operator(-), &
        operator(+), diag, mm_init=>init
    USE STRINGS, only: itoa
    implicit none
    integer(IK), intent(in)    :: irr,n,n_c,n_gr
    real(RK)   , intent(inout) :: sdsto(:,:) ! sec der storage
    logical    , intent(in)    :: mask(:,:)  ! (n_gr,n_gr) proc ids for each sec der
    !*** end of interface ***

    integer(IK)   , parameter       :: X=1,Y=2,XY=3
    ! 0: value; X,Y: first derivatives; XY: (mixed) second derivative.

    type(rmatrix)                   :: UF, UB
    type(rmatrix) , dimension(0:XY) :: S, T, V, O, VR ! O=PVSP
    type(rmatrix) , dimension(0:XY) :: P
    type(rdmatrix), dimension(0:XY) :: p2
    type(rdmatrix), dimension(0:XY) :: tp,e,a,k,ak,r2,m2
    type(rmatrix) , dimension(0:XY) :: aVa, akOka
    type(rmatrix) , dimension(0:XY) :: RW
    type(rmatrix)                   :: VRC ! contracted
    type(rmatrix)                   :: DMAT

    integer(IK) :: i
    integer(IK) :: gx,gy
    real(RK)    :: drv(0:XY) ! forces (derivatives)

    FPP_TIMER_DECL(tot)
    FPP_TIMER_DECL(cod2)
    FPP_TIMER_DECL(cod1)
    FPP_TIMER_DECL(cod0)
    FPP_TIMER_DECL(dio)

    FPP_TIMER_START(tot)
    DPRINT MyID,'do_sdr_trafo:   irr=',irr,' dim=',n,' cdim=',n_c

    ! Prepare (contracted) density matrix (in orbital space).
    ! (density_data seems to be shut down)
    call alloc(n_c, DMAT)
    call densmat(irr, DMAT)

    do i = 0, XY
       ! storage  for  ovrl,  kin,  nuc  and  pvsp.  FIXME:  redundant
       ! allocation of V(0), O(0) here:
       call alloc( n, S(i), T(i), V(i), O(i) )
    enddo

    ! contracted:
    call alloc( n_c, VRC )

    !
    ! Recover matrices saved in do_one_block() called from ham_trafo()
    ! in pre-scf. These are  already in mom-space, get_2() may require
    ! communication between all processes in MPI world:
    !
    p2(0) = rdmatrix(get_1(n, RD_P2DIAG, irr))
    UF = rmatrix(get_2(n, RM_FORWARD, irr))
    UB = rmatrix(get_2(n, RM_BACKWARD, irr))
    V(0) = rmatrix(get_2(n, RM_PNUC, irr))
    O(0) = rmatrix(get_2(n, RM_PPVSP, irr))

    V(0) = - V(0) ! FIXME: sign
    O(0) = - O(0) ! FIXME: sign

    ! trafo of DKH hamiltonian in internal sub:
    call code0(irr)

    do gy=1,n_gr

       ! trafo of Y gradient in internal sub:
       call code1(irr,gy,Y) ! see contains

    do gx=1,gy !n_gr
       !if ( gx < gy ) CYCLE
       if( .not. mask(gx,gy) ) CYCLE

       DPRINT   MyID,'XXX: SECDER(',gx,gy,')'

       ! trafo of X gradient in internal sub:
       call code1(irr,gx,X) ! see contains

       ! trafo of XY sec der in internal sub:
       call code2(irr,gx,X,gy,Y) ! see contains

       ![[==== Contract and sum up XY-derivative: ==================
       ! contract the basis:
       call contract_rm(irr, VR(XY), VRC)
       drv(XY) = trace2(VRC, DMAT)

#if 0 /* DEBUG ONLY */
       ! for debug only:
       VR(X) = VR(X) - bk(VR(0), P(X))
       VR(X) = sim( VR(X), UB )
       call contract_rm(irr, VR(X), VRC)
       drv(X) = trace2(VRC,DMAT)

       print *, MyID,'DDD: b) GRADIENT(',irr,gx,   ')=',drv(X)
       print *, MyID,'DDD: b)   SECDER(',irr,gy,gx,')=',drv(XY)
#endif /* DEBUG ONLY */

       sdsto(gy,gx) = sdsto(gy,gx) - drv(XY) ! FIXME: again V=-V?
       if( gx /= gy )then
       sdsto(gx,gy) = sdsto(gx,gy) - drv(XY) ! FIXME: again V=-V?
       endif
       !]]==========================================================
    enddo ! X == gx == i
    enddo ! Y == gy == j

    DPRINT MyID,'do_sdr_trafo: exit'
    FPP_TIMER_STOP(tot)
#ifdef FPP_TIMERS
    print *,MyID,'RELGR: timing(',irr,')=',FPP_TIMER_VALUE(tot)
    print *,MyID,'RELGR: |- code 2      :',FPP_TIMER_VALUE(cod2)
    print *,MyID,'RELGR: |- code 1      :',FPP_TIMER_VALUE(cod1)
    print *,MyID,'RELGR: |- code 0      :',FPP_TIMER_VALUE(cod0)
    print *,MyID,'RELGR: disk io        :',FPP_TIMER_VALUE(dio)
#endif
    contains ! code0, code1 and code2

      subroutine code0(irr)
       ! internal sub that processes (in host storage)
       ! the DKH hamiltonian
       implicit none
       integer(IK), intent(in) :: irr
       ! *** end of interface ***

       integer(IK) :: ip

       FPP_TIMER_START(cod0)
       DPRINT MyID,'RELSD: code0(',irr,')'

       ![[==== Kinematic factors (prep) ============================
       ! compute values kinematic factors:
       call p2_diag_ders(0, p2, tp, e, a, k, ak, r2)
       !]]==========================================================

       ![[==== DKH1 hamiltonian (prep) =============================
       ! DKH1 level: V := aVa + akOka
       aVa(0)   =  a(0) * V(0) *  a(0)

       akOka(0) = ak(0) * O(0) * ak(0)

       VR(0)    = aVa(0) + akOka(0)
       !]]==========================================================

       ![[==== DKH2 hamiltonian (prep) =============================
       ! DKH2 level:
       ! we are going to compute
       !   E2 = - WEW - {W^2,E}/2
       ! using the identity
       !   W^2 = WR * (R^-2) * RW = - (RW)^T * 1/(k2p2) * RW
       ! and finally:
       !  E2 = (1/2)[ ( R^-2 * ( E * RW + RW * E ) )^T * RW + h.c. ]

       ! RW = r2 * aVa - arVra:
       RW(0) = r2(0) * aVa(0) - akOka(0)

       RW(0) = RPT(RW(0), e(0))

       ! r2: R^ 2 =   k2p2
       ! m2: R^-2 = 1/k2p2, read as ``R minus 2''
       m2(0) =            r2(0)**(-1)
      !m2(1) =  - r2(1) * r2(0)**(-2)

       ! FIXME: S(0), T(0) used as a temp storage:
       ! anticommutator {e,RW}:
       S(0) = RW(0) * e(0) + e(0) * RW(0)

       T(0) = tr( m2(0) * S(0) )  * RW(0) ! +h.c.)/2, see below:

       VR(0) = VR(0) - HALF * ( T(0) + tr(T(0)) )

       ! add diagonal kinetic energy:
       do ip=1,n
          VR(0)%m(ip,ip) =  - tp(0)%d(ip) + VR(0)%m(ip,ip)
       enddo
       !]]==========================================================
       FPP_TIMER_STOP(cod0)
      end subroutine code0

      subroutine code1(irr,gx,X)
       ! internal sub that processes (in host storage)
       ! the X-gradient
       implicit none
       integer(IK), intent(in) :: irr,gx,X
       ! *** end of interface ***

       integer(IK) :: ip

       FPP_TIMER_START(cod1)
       DPRINT MyID,'RELSD: code1(',irr,gx,X,')'

       ![[==== PT1 for X ===========================================
                                       FPP_TIMER_START(dio)
       ! fetch gradients of overlap, kinetic, of V and PVSP :
       call read_grads(irr, gx, S(X), T(X), V(X), O(X))
                                       FPP_TIMER_STOP(dio)

       ! to the momentum basis:
       S(X)  = - sim(S(X),UF)
       T(X)  = - sim(T(X),UF) * 2.0_rk ! since p^2 = 2 * T
       V(X)  =   sim(V(X),UF)
       O(X)  =   sim(O(X),UF)

       ! build the response generator, and response of eigenvalues:
       call PT( p2(0), T(X), S(X), p2(X), P(X) )
       ! NOTE: S(X), T(X) now may be used as temp storage invalidated
       !       at the next X-loop
       !]]==========================================================

       ![[==== X-grads of V in mom. bas.: ==========================
       ! A0' = A0
       ! A1' = A1 + [A0,W1]
       V(X)   = V(X)  +       bk( V(0), P(X) )
       O(X)   = O(X)  +       bk( O(0), P(X) )
       !]]==========================================================

       ![[==== Kinematic factors (prep) ============================
       ! compute values and X-grads of kinematic factors:
       call p2_diag_ders(1,p2(0:X:X),tp(0:X:X),e(0:X:X),a(0:X:X),k(0:X:X),ak(0:X:X),r2(0:X:X))
       !]]==========================================================

       ![[==== DKH1 X-gradient =====================================
       ! DKH1 level: V := aVa + akOka
       aVa(X)   =  a(0) * V(X) *  a(0) &
                +  a(X) * V(0) *  a(0) &
                +  a(0) * V(0) *  a(X)

       akOka(X) = ak(0) * O(X) * ak(0) &
                + ak(X) * O(0) * ak(0) &
                + ak(0) * O(0) * ak(X)

       VR(X)    = aVa(X) + akOka(X)
       !]]==========================================================

       ![[==== DKH2 X-gradient =====================================
       ! RW = r2 * aVa - arVra:
       RW(X) = r2(0) * aVa(X) - akOka(X) &
             + r2(X) * aVa(0)

       ! WARNING: RW(0) already redefined to RW-TILDE!
       !          RPT_gr_weight() expects that
       RW(X) =           RPT( RW(X)  ,e(0)      ) &
             + RPT_gr_weight( RW(0)  ,e(0), e(X))

       ! R2 : R^ 2 =   k2p2
       ! Rm2: R^-2 = 1/k2p2, read as ``R minus 2''
      !m2(0) =            r2(0)**(-1)
      !m2(X) =  - r2(X) * r2(0)**(-2)
       m2(X) =  - ( r2(X) * m2(0) ) * m2(0) ! FIXME: R^-4 may be too big a number!

       ! FIXME: S(X), T(X) are used as temp storage:
       S(X) = RW(X) * e(0) + e(0) * RW(X) + RW(0) * e(X) + e(X) * RW(0)

       ! FIXME: S(0) still holds RW * e + e * RW:
       T(X)  = tr( m2(X) * S(0) ) * RW(0) &
             + tr( m2(0) * S(X) ) * RW(0) &
             + tr( m2(0) * S(0) ) * RW(X) ! +h.c.)/2, see below:

       VR(X) = VR(X) - HALF * ( T(X) + tr(T(X)) )

       ! add diagonal kinetic energy:
       do ip=1,n
          VR(X)%m(ip,ip) =  - tp(X)%d(ip) + VR(X)%m(ip,ip)
       enddo
       !]]==========================================================
       FPP_TIMER_STOP(cod1)
      end subroutine code1

      subroutine code2(irr,gx,X,gy,Y)
       ! internal sub that processes (in host storage)
       ! the XY-derivative
       implicit none
       integer(IK), intent(in) :: irr,gx,X,gy,Y
       ! *** end of interface ***

       integer(IK) :: ip

       FPP_TIMER_START(cod2)
       DPRINT MyID,'RELSD: code2(',irr,gx,X,gy,Y,')'

       ![[==== PT2 for XY ==========================================
       ! fetch derivatives of overlap, kinetic, of V and PVSP :
       DPRINT  MyID,'RELSD: read dervs',gx,gy
       ASSERT(gx<=gy)
                                       FPP_TIMER_START(dio)
       call read_dervs(irr, gx, gy, S(XY), T(XY), V(XY), O(XY))
                                       FPP_TIMER_STOP(dio)

       ! to the momentum basis:
       S(XY)  =   sim(S(XY),UF)
       T(XY)  =   sim(T(XY),UF) * 2.0_rk ! since p^2 = 2 * T
       V(XY)  = - sim(V(XY),UF)
       O(XY)  = - sim(O(XY),UF)

       ! H2' = T2 + [e1,W1] - [[E,W1],W1]/2 + [E,W2]: off=0, diag=e2
       !       ^^^^^^^^^^^^^\/^^^^^^^^^^^^^      ??
       T(XY)  = T(XY) +     bk(p2(X), P(Y) )             &
                      +     bk(p2(Y), P(X) )             &
                      - bk( bk(p2(0), P(X) ), P(Y) ) * 0.5_rk &
                      - bk( bk(p2(0), P(Y) ), P(X) ) * 0.5_rk

       ! S2' = S2           - [[1,W1],W1]/2 + [1,W2]: off=0, diag=0
       !       ^^^^^^^^^^^^^\/^^^^^^^^^^^^^      ??
       S(XY)  = S(XY) - bk( bk( P(X)), P(Y) ) * 0.5_rk &
                      - bk( bk( P(Y)), P(X) ) * 0.5_rk

       ! solve for W2 = P(XY) and e2 = p2(XY):
       call PT( p2(0), T(XY), S(XY), p2(XY), P(XY) )
       ! NOTE: S(XY), T(XY) now may be used as temp storage
       !]]==========================================================

       ![[==== XY-dervs of V in mom. bas.: =========================
       !
       ! A2' = A2 + [A1 ,W1] + [[A0,W1],W1]/2 + [A0,W2] =
       !     = A2 + [A1',W1] - [[A0,W1],W1]/2 + [A0,W2].
       !               ^     ^
       ! Note that the first derivatives have been already overwritten
       ! by their corrected counterparts in MOVING mom basis!
       ! => minus by double "brakets"!

       V(XY)  = V(XY) +     bk( V(X), P(Y) )             &
                      +     bk( V(Y), P(X) )             &
                      - bk( bk( V(0), P(X) ), P(Y) ) * 0.5_rk &
                      - bk( bk( V(0), P(Y) ), P(X) ) * 0.5_rk &
                      +     bk( V(0), P(XY))

       O(XY)  = O(XY) +     bk( O(X), P(Y) )             &
                      +     bk( O(Y), P(X) )             &
                      - bk( bk( O(0), P(X) ), P(Y) ) * 0.5_rk &
                      - bk( bk( O(0), P(Y) ), P(X) ) * 0.5_rk &
                      +     bk( O(0), P(XY))
       !]]==========================================================

       ![[==== X-, Y-, and XY-derivatives of kin. fac.: ============
       ! compute values and derivative of kinematic factors:
       call p2_diag_ders(2, p2, tp, e, a, k, ak, r2)
       !]]==========================================================

       ![[==== DKH1 XY-derivative ==================================
       ! DKH1 level: V := aVa + akOka
       aVa(XY)  =  a(0)   * V(XY) *  a(0) &
                +  a(X)   * V(Y)  *  a(0) &
                +  a(0)   * V(Y)  *  a(X) &
                +  a(Y)   * V(X)  *  a(0) &
                +  a(XY)  * V(0)  *  a(0) &
                +  a(Y)   * V(0)  *  a(X) &
                +  a(0)   * V(X)  *  a(Y) &
                +  a(X)   * V(0)  *  a(Y) &
                +  a(0)   * V(0)  *  a(XY)

       akOka(XY)=  ak(0)  * O(XY) *  ak(0) &
                +  ak(X)  * O(Y)  *  ak(0) &
                +  ak(0)  * O(Y)  *  ak(X) &
                +  ak(Y)  * O(X)  *  ak(0) &
                +  ak(XY) * O(0)  *  ak(0) &
                +  ak(Y)  * O(0)  *  ak(X) &
                +  ak(0)  * O(X)  *  ak(Y) &
                +  ak(X)  * O(0)  *  ak(Y) &
                +  ak(0)  * O(0)  *  ak(XY)

       VR(XY)    = aVa(XY) + akOka(XY)
       !]]==========================================================

       ![[==== DKH2 XY-derivative ==================================
       ! RW = r2 * aVa - arVra:
       RW(XY)= r2(0)  * aVa(XY) - akOka(XY) &
             + r2(X)  * aVa(Y)              &
             + r2(Y)  * aVa(X)              &
             + r2(XY) * aVa(0)

       ! WARNING: RW(0) and RW(1) already redefined to RW-TILDE!
       !          RPT_gr_weight() expects that
       RW(XY)=           RPT( RW(XY), e(0)       ) &
             + RPT_gr_weight( RW(X) , e(0), e(Y) ) &
             + RPT_gr_weight( RW(Y) , e(0), e(X) ) &
             + RPT_gr_weight( RW(0) , e(0), e(XY))

       ! r2: R^ 2 =   k2p2
       ! m2: R^-2 = 1/k2p2, read as ``R minus 2''
       m2(0) =             r2(0)**(-1)
      !m2(X) =  - r2(X)  * r2(0)**(-2)
      !m2(Y) =  - r2(Y)  * r2(0)**(-2)
      !m2(XY)=  - r2(XY) * r2(0)**(-2) + 2 * (r2(X) * r2(Y) * r2(0)**(-3))
       m2(X) =  - ( r2(X) * m2(0) ) * m2(0) ! FIXME: R^-4 may be too big a number
       m2(Y) =  - ( r2(Y) * m2(0) ) * m2(0) ! FIXME: R^-4 may be too big a number
       m2(XY)= (-  r2(XY) * m2(0) + 2 * ( r2(X) * m2(0) ) * ( r2(Y) * m2(0) ) ) * m2(0)

       ! temp := RW * e + e * RW:
       ! FIXME: S(XY), T(XY) are used as temp storage:
       S(XY)= RW(XY) * e(0)    +    e(0)  * RW(XY) &
            + RW(0)  * e(XY)   +    e(XY) * RW(0)  &
            + RW(X)  * e(Y)    +    e(Y)  * RW(X)  &
            + RW(Y)  * e(X)    +    e(X)  * RW(Y)

       ! FIXME: S(0), S(1) still hold RW * e + e * RW ???
       T(XY) = tr( m2(XY) * S(0)  ) * RW(0) &
             + tr( m2(Y)  * S(X)  ) * RW(0) &
             + tr( m2(Y)  * S(0)  ) * RW(X) &
             + tr( m2(X)  * S(Y)  ) * RW(0) &
             + tr( m2(0)  * S(XY) ) * RW(0) &
             + tr( m2(0)  * S(Y)  ) * RW(X) &
             + tr( m2(X)  * S(0)  ) * RW(Y) &
             + tr( m2(0)  * S(X)  ) * RW(Y) &
             + tr( m2(0)  * S(0)  ) * RW(XY) ! +h.c.)/2, see below:

       VR(XY) = VR(XY) - HALF * ( T(XY) + tr(T(XY)) )

       ! add diagonal kinetic energy:
       do ip=1,n
          VR(XY)%m(ip,ip) =  - tp(XY)%d(ip) + VR(XY)%m(ip,ip)
       enddo
       !]]==========================================================

       ![[==== XY-dervs back to real space: ========================
       !
       ! First account for the changes in mom basis:
       !
       ! A2  = A2' - [A1',W1] + [[A0',W1],W1]/2 - [A0',W2]
       !       ( W => -W of the above counterpart)
       ! Note that the lower derivatives on the rhs
       ! are still in MOVING mom basis!

       VR(XY) = VR(XY) -     bk( VR(X), P(Y) )             &
                       -     bk( VR(Y), P(X) )             &
                       + bk( bk( VR(0), P(X) ), P(Y) ) * 0.5_rk &
                       + bk( bk( VR(0), P(Y) ), P(X) ) * 0.5_rk &
                       -     bk( VR(0), P(XY))

       VR(XY) = sim( VR(XY), UB )
       !]]==========================================================
       FPP_TIMER_STOP(cod2)
      end subroutine code2
  end subroutine do_sdr_trafo
#endif

  subroutine PT_gen(n,e,V,S,d,W)
    !
    ! Perturbation Theory Response Generator:
    !   p /= q:
    !          W(p,q) = [ V(p,q) - S(p,q) * e(q) ] / ( e(q) - e(p) )
    !   p == q:
    !          W(q,q) = - S(q,q) / 2
    !
    ! And response of Eigenvalues:
    !          d(q)   = V(q,q) - S(q,q) * e(q)
    !
    implicit none
    integer(IK), intent(in)  :: n
    real(RK)   , intent(in)  :: V(n,n), S(n,n), e(n)
    real(RK)   , intent(out) :: d(n)   ! result
    real(RK)   , intent(out) :: W(n,n) ! result
    ! *** end of interface ***

    integer(IK) :: p,q

    do q=1,n
       do p=1,q-1
          if( .not. equal(e(q),e(p)) )then
             W(p,q) = ( V(p,q) - S(p,q) * e(q) ) / ( e(q) - e(p) )
             W(q,p) = - W(p,q) - S(p,q)
          else
             W(p,q) = zero
             W(q,p) = zero
          endif
       enddo
       W(q,q) =        - S(q,q) / 2
       d(q)   = V(q,q) - S(q,q) * e(q)
    enddo
  contains
    function equal(r1,r2) result(yes)
      ! same as in relativisitc_trafo_module:
      implicit none
      real(RK),intent(in) :: r1,r2
      logical             :: yes !<<<result
      ! *** end of interface ***

      real(RK),parameter :: eps = 1.0E-7_rk
      real(RK)           :: r1r2

      r1r2 = r1*r2
      if(r1r2/=zero)then
         yes = ( abs((r1-r2)**2/r1r2) < eps**2 )
      else
         yes = ( r1.eq.r2 )
      endif
    end function equal
  end subroutine PT_gen

  subroutine PT(e, V, S, d, W)
    !
    ! Perturbation Theory Response Generator,
    ! and response of Eigenvalues.
    !
    use matrix_module, only: rmatrix, rdmatrix, alloc
    implicit none
    type(rdmatrix), intent(in) :: e ! 0-in
    type(rmatrix), intent(in) :: V, S
    type(rdmatrix), intent(out) :: d ! 1-out
    type(rmatrix), intent(out) :: W
    ! *** end of interface ***

    integer(IK) :: n(1)

    n = shape(e%d)

    call alloc(n(1), d)
    call alloc(n(1), W)

    call PT_gen(n(1), e%d, V%m, S%m, d%d, W%m)
  end subroutine PT

  function trace2(V,P) result(res)
    !
    ! SUM(i) SUM(j) { V(i,j) * P(j,i) }
    !
    use matrix_module, only: rmatrix
    implicit none
    type(rmatrix),intent(in) :: V,P
    real(RK)                 :: res
    ! *** end of interface ***

    ! FIXME: only for symmetric?
    res = SUM( V%m * P%m )
  end function trace2

  subroutine densmat(irr,DMAT)
    ! fetch gradients of overlap and kinetic energy,
    ! fetch gradients of V and PVSP.
    use matrix_module , only: rmatrix
    use density_data_module, only: gendensmat_tot
!   USE STRINGS, only: itoa
!   USE IO
    implicit none
    integer(IK)     , intent(in)    :: irr
    type(rmatrix)   , intent(inout) :: DMAT
    ! *** end of interface ***

!      call io_set_error_handler(IO_CONTINUE)
!      call read_buffer('/home/matveev/dmat.dat'//trim(itoa(irr)),DMAT%m)
!      if ( io_status() /= 0 )then
    call gendensmat_tot(irr,DMAT%m)
!      call write_buffer('/home/matveev/dmat.dat'//trim(itoa(irr)),DMAT%m)
!      print *,'DDD: densmat written to dmat.dat'//trim(itoa(irr))
!      else
!      print *,'DDD: densmat taken from dmat.dat'//trim(itoa(irr))
!      endif
!      call io_set_error_handler()
  end subroutine densmat

  subroutine read_relham(irr, n, S, T, V, O)
    !
    ! Fetch matrices of overlap, kinetic energy, V and PVSP.
    !
    use MATRIX_IMPL, only: rmatrix, matrix
    use relgrads_store, only: rg_get
    implicit none
    integer(IK), intent(in) :: irr, n
    type(rmatrix), intent(out) :: S, T, V, O
    ! *** end of interface ***

    real(RK), dimension(n, n) :: S_, T_, V_, O_

    call rg_get(irr, S_, T_, V_, O_)

    ! return rmatrices:
    S = matrix(S_)
    T = matrix(T_)
    V = matrix(V_)
    O = matrix(O_)
  end  subroutine read_relham

  function read_one(irr, n, typ) result(V)
    !
    ! Fetch matrix of (pseudpotential) matrix elements.
    !
    use MATRIX_IMPL, only: rmatrix, matrix
    use relgrads_store, only: rg_flag, rg_xx_get
    implicit none
    integer(IK), intent(in) :: irr, n, typ
    type(rmatrix) :: V
    ! *** end of interface ***

    real(RK) :: V_(n, n)

    ASSERT(rg_flag(typ))

    call rg_xx_get(0, irr, 0, 0, V_, typ=typ)

    ! return an rmatrix:
    V = matrix(V_)
  end function read_one

  subroutine read_grads(irr,igr,S,T,V,O)
    ! fetch gradients of overlap and kinetic energy,
    ! fetch gradients of V and PVSP.
    use matrix_module , only: rmatrix
    use relgrads_store, only: rg_gr_get, rg_flag, rg_xx_get, RGPSEU
    implicit none
    integer(IK)     , intent(in)    :: irr, igr
    type(rmatrix)   , intent(inout) :: S,T,V,O
    optional T,V,O
    ! *** end of interface ***

    if( present(T) )then
      ASSERT(present(V))
      ASSERT(present(O))
      call rg_gr_get(irr,igr,S%m,T%m,V%m,O%m)
    else
      if( rg_flag(RGPSEU) )then
        call rg_xx_get(1,irr,igr,0,S%m,typ=RGPSEU)
      endif
    endif
  end subroutine read_grads

  subroutine read_dervs(irr,x,y,S,T,V,O)
    ! fetch gradients of overlap and kinetic energy,
    ! fetch gradients of V and PVSP.
    use matrix_module , only: rmatrix
    use relgrads_store, only: rg_sd_get, rg_flag, rg_xx_get, RGPSEU
    implicit none
    integer(IK)     , intent(in)    :: irr, x, y
    type(rmatrix)   , intent(inout) :: S,T,V,O
    optional T,V,O
    ! *** end of interface ***

    if( present(T) )then
      ASSERT(present(V))
      ASSERT(present(O))
      call rg_sd_get(irr,x,y,S%m,T%m,V%m,O%m)
    else
      if( rg_flag(RGPSEU) )then
        call rg_xx_get(2,irr,x,y,S%m,typ=RGPSEU)
      endif
    endif
  end subroutine read_dervs

  function contract(irr, nc, um) result(cm)
    !
    ! Functional  version, uses  corresponding subroutine.   Returns a
    ! plain array  as a contraction  procedure is not  yet implemented
    ! for  distributed  matrices and  conversion  matrix  -> array  ->
    ! matrix is communication expensinve.
    !
    use MATRIX_IMPL, only: rmatrix, array
    use contraction_module, only: do_contract => contract
    implicit none
    integer(IK), intent(in) :: irr, nc
    type(rmatrix), intent(in) :: um
    real(RK) :: cm(nc, nc)
    ! *** end of interface ***

    !
    ! This one works with plain arrays:
    !
    call do_contract(irr, array(um), cm)
  end function contract

  subroutine contract_rm(irr, um, cm)
    !
    ! Subroutine version, uses array version.
    !
    use matrix_module, only: rmatrix, square
    use contraction_module, only: contract
    implicit none
    integer(IK), intent(in) :: irr
    type(rmatrix), intent(in) :: um
    type(rmatrix), intent(inout) :: cm ! needs to be allocated
    ! *** end of interface ***

    if(irr.eq.NO_IRREP)return

    ASSERT(square(um))
    ASSERT(square(cm))
    call contract(irr, um%m, cm%m)
  end subroutine contract_rm

  subroutine write_rel_ham(irr, S, T, V)
    !
    ! Store relativistic matrices at the location the rest of the code
    ! is  goiing to  look  for them.   NOTE:  cumulative counter  runs
    ! through the joint storage array of triagonal irrep matices.
    !
    ! FIXME:  master/slave assymetry here.  Slaves only  save overlap,
    ! while master also saves kin and nuc:
    !
    use comm, only: comm_rank
    use options_module, only: options_integrals_on_file
    use integralstore_module, only: integralstore_2cob_ol , &
         integralstore_2cob_kin , integralstore_2cob_nuc
    implicit none
    integer(IK), intent(in) :: irr
    real(RK), intent(in) :: S(:, :), T(:, :), V(:, :)
    ! *** end of interface ***

    integer(IK) :: sz(2), n, m, rank
    integer(IK), save :: last_irr = 0
    integer(IK), save :: counter = 0

    rank = comm_rank()

    sz = shape( S )
    ASSERT(sz(1)==sz(2))

    ! reset the initial address:
    if (irr == 1) then
      last_irr = 0
      counter  = 0
    endif

    ! make sure irreps come in order:
    ASSERT(irr==last_irr+1)
    last_irr = irr

    if(options_integrals_on_file()) then
       !
       ! Integrals on disk:
       !
       do m = 1, sz(1)
          do n = 1, m
             ! overlap:
             write(io_overlap) S(m, n)
             if (rank == 0) then
                ! kin-nuc, nuc with "-" historically:
                write(io_ham_kin_nuc) T(m, n), -V(m, n)
             endif
          end do
       end do
    else
       !
       ! Integrals in memory:
       !
ASSERT(allocated(integralstore_2cob_ol))
if (rank == 0) then
ASSERT(allocated(integralstore_2cob_kin))
ASSERT(allocated(integralstore_2cob_nuc))
endif

       do m = 1, sz(1)
          do n = 1, m
             counter = counter + 1
             ! overlap:
             integralstore_2cob_ol(counter) = S(m, n)
             if (rank == 0) then
                ! kin-nuc:
                integralstore_2cob_kin(counter) =  T(m, n)
                integralstore_2cob_nuc(counter) = -V(m, n)
             endif
          enddo
       enddo
    endif
  end subroutine write_rel_ham

  subroutine open_output_tapes()
    use comm, only: comm_rank
    use filename_module, only: tmpfile
    use options_module, only: options_integrals_on_file
    use iounitadmin_module, only: openget_iounit
    implicit none
    ! *** end of interface ***

    if (options_integrals_on_file()) then
       ! overlap:
       io_overlap=openget_iounit&
            (trim(tmpfile("overlap.dat")), form='unformatted',&
            status='unknown')
       if (comm_rank() == 0) then
          ! kin-nuc:
          io_ham_kin_nuc = openget_iounit&
               (trim(tmpfile("ham_kin_nuc.dat")), form='unformatted',status='unknown')
       endif
    end if
  end subroutine open_output_tapes

  subroutine close_output_tapes()
    use comm, only: comm_rank
    use iounitadmin_module, only: returnclose_iounit
    use options_module, only: options_integrals_on_file
    implicit none
    ! *** end of interface ***

    if (options_integrals_on_file()) then
       ! overlap:
       call returnclose_iounit(io_overlap)
       if (comm_rank() == 0) then
          ! kin-nuc:
          call returnclose_iounit(io_ham_kin_nuc)
       endif
    end if
  end subroutine close_output_tapes

  !**************************************************************
  !
  ! interface bk
  !   module procedure bk_r_r  ( rmatrix , rmatrix )
  !   module procedure bk_rd_r ( rdmatrix, rmatrix )
  !   module procedure bk_r    (           rmatrix ) = bk(1, rmatrix)
  ! end interface
  !
  ! Perturbation Theory Response Brackets ( ``Commutator'' ):
  !
  !   exp(W)^T * V * exp(W) = V + ( W^T * V + V * W ) + ...
  !                               \_______. ._______/
  !                                        V
  !                               linear response dV
  ! For symmetric V:
  !
  !   VW = V * W + transpose thereof
  !
  ! We also used to have anti-symetric W (in prev. impl.)
  !
  !    W = - W^T                      (not the case!)
  !
  ! and computed
  !
  !   VW = V * W - W * V = [V,W]      (not the case!)
  !

  function bk_r_r(V, W) result(VW)
    ! Perturbation Theory Response Brackets ( ``Commutator'' ):
    use matrix_module, only: rmatrix, operator(*), tr, square
    use matrix_module, only: operator(+)
    implicit none
    type(rmatrix), intent(in) :: V, W
    type(rmatrix) :: VW
    ! *** end of interface ***

    ASSERT(square(V))
    ASSERT(square(W))

    VW = V * W
    VW = VW + tr(VW)
  end function bk_r_r

  function bk_rd_r(V, W) result(VW)
    ! Perturbation Theory Response Brackets ( ``Commutator'' ):
    use matrix_module, only: rdmatrix, rmatrix, operator(*), tr, square
    use matrix_module, only: operator(+)
    implicit none
    type(rdmatrix), intent(in) :: V
    type(rmatrix), intent(in) :: W
    type(rmatrix) :: VW
    ! *** end of interface ***

    ASSERT(square(W))

    VW = V * W
    VW = VW + tr(VW)
  end function bk_rd_r

  function bk_r(W) result(VW)
    ! Perturbation Theory Response Brackets ( ``Commutator'' ):
    use matrix_module, only: rmatrix, operator(+), tr, square
    implicit none
    type(rmatrix), intent(in) :: W
    type(rmatrix) :: VW
    ! *** end of interface ***

    ASSERT(square(W))

    VW = W + tr(W)
  end function bk_r

  subroutine put_1(id, irr, d)
    use MATRIX_IMPL, only: rdmatrix
    implicit none
    integer, intent(in) :: id, irr
    type(rdmatrix), intent(in) :: d
    ! *** end of interface ***

    ASSERT(id==RD_P2DIAG)
    dmatrices(id, irr) = d
  end subroutine put_1

  subroutine put_2(id, irr, m)
    use MATRIX_IMPL, only: rmatrix
    implicit none
    integer, intent(in) :: id, irr
    type(rmatrix), intent(in) :: m
    ! *** end of interface ***

    rmatrices(id, irr) = m
  end subroutine put_2

  function get_1(n, id, irr) result(d)
    !
    ! FIXME: depending on MATRIX_IMPL requires converting a
    ! distributed object.
    !
    use MATRIX_IMPL, only: array
    implicit none
    integer, intent(in) :: n, id, irr
    real(RK) :: d(n)
    ! *** end of interface ***

    ASSERT(id==RD_P2DIAG)
    d = array(dmatrices(id, irr))
  end function get_1

  function get_2(n, id, irr) result(m)
    !
    ! FIXME: depending on MATRIX_IMPL requires converting a
    ! distributed object.
    !
    use MATRIX_IMPL, only: array
    implicit none
    integer, intent(in) :: n, id, irr
    real(RK) :: m(n, n)
    ! *** end of interface ***

    m = array(rmatrices(id, irr))
  end function get_2

  !--------------- End of module ----------------------------------
end module relgrads

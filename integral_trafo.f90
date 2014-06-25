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
subroutine integral_trafo(mode)
  !----------------------------------------------------------------
  !
  !  Purpose: an entree to IntegralTransformations
  !           for all processors
  !
  !
  !  Subroutine called by:
  !     main_integral()
  !     (on slaves via main_slave() deamon)
  !
  ! Copyright (c) Alexei Matveev
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

  !------------ Modules used --------------------------------------
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  use error_module
  use integralpar_module, only: integralpar_spor
  use spin_orbit_module,  only: is_on,op_FitTrafo
  use fit_trafo_module,   only: fit_trafo
  use comm, only: comm_barrier
  use options_module   , only: options_debug_key &
                             , options_relativistic
  use relgrads, only: rel_trafo_sr, IENRG, IGRAD, ISDER
  use reltrafo, only: rel_trafo ! spin-orbit version
  use gradient_data_module, only: gradient_totalsym &
                                , dervs_totalsym
  use interfaces, only: RELENRG, RELGRAD, RELSDER ! accept as args
  use relgrads_store, only: rg_close, rg_bcast, rg_distr
  implicit none
  integer(IK), intent(in) :: mode ! RELENRG, RELGRAD, RELSDER

  !== Interrupt end of public interface of module =================
  !------------ Declaration of formal parameters ------------------

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of subroutines ------------------------

  !------------ Declaration of local constants --------------------

  !------------ Declaration of local variables --------------------

  integer(IK) :: iop
  logical     :: cleanup

  FPP_TIMER_DECL(rt)
  !----------------------------------------------------------------
  !------------ Executable code -----------------------------------


  DPRINT MyID,"integral_trafo: entered"

  ASSERT(options_relativistic)

  DPRINT MyID,"integral_trafo: entree sync ..."
  call comm_barrier()
  DPRINT MyID,"integral_trafo:             ... synced"
  !---------------------------------
  !
  ! RELATIVISTIC TRANSFORMATION:
  ! (call moved from main_integral)
  !
  FPP_TIMER_START(rt)

  ![[=== Integral transformations: ===

  ! deallocate integrals after trafos by default:
  cleanup = .true.

  if( integralpar_spor )then
     ! not yet parallelized:
     call rel_trafo()
  else

     iop = 2 ! second order DKH

     if( IAND( mode, RELENRG ) /= 0 )then
        iop = iop + IENRG

        ! do not deallocate integrals till after SCF
        ! in a gradients run:
        if( IAND( mode, RELGRAD ) /= 0 ) cleanup = .false.
     endif

     if( IAND( mode, RELGRAD ) /= 0 )then
        iop = iop + IGRAD
        if( IAND( mode, RELSDER ) /= 0 )then
           iop = iop + ISDER
           ! gradient+secder trafo in a secder run :

           ! DONT call rg_bcast(0) ! has alredy been distributed before SCF!

           ! broadcast first derivatives to all processors:
           call rg_bcast(1)
           ! distribute second derivatives as required later in rel_trafo_sr:
           call rg_distr(2) ! DONT call rg_bcast(2)

           call rel_trafo_sr( iop, GRSTO=gradient_totalsym &
                                 , SDSTO=dervs_totalsym    )
        else
           ! a) energy trafo in a gradient run (before SCF):
           ! b) gradient trafo in a gradient run (after SCF):

           if( IAND( mode, RELENRG ) /= 0 )then
             ! Relativistic Derivative Storage: re-distribute quads
             call rg_bcast(0)
           else
             call rg_bcast(1)
           endif

           call rel_trafo_sr( iop, GRSTO=gradient_totalsym )
        endif
     else
           ! normal reltrafo in a normal run (before SCF):

           ! Relativistic Derivative Storage: re-distribute quads
            call rg_bcast(0)

           call rel_trafo_sr( iop )
     endif
  endif

  ! Relativistic Derivative Storage: clean up
  if( cleanup )then
    ! clean up the storage for relativistic integrals
    call rg_close('d') ! delete
  endif
  !]]================================

  FPP_TIMER_STOP(rt)
  FPP_TIMER_PRINT(rt)

  !
  ! RELATIVISTIC FIT TRANSFORMATION:
  !
  if(is_on(op_FitTrafo)) then
     DPRINT MyID,"integral_trafo: call fit_trafo()"
     call fit_trafo()
     DPRINT MyID,"integral_trafo: ."
  endif

  DPRINT MyID,"integral_trafo: exit"
end subroutine integral_trafo

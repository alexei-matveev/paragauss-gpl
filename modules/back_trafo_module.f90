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
module back_trafo_module
  !-------------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  !
  !-------------------------------------------------------------------
  !== Interrupt of public interface of module ========================
  !-------------------------------------------------------------------
  ! Modifications
  !-------------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: ...
  ! Date:   ...
  ! Description: ...
  !
  !-------------------------------------------------------------------

# include "def.h"
  use type_module, only: IK=>i4_kind, RK=>r8_kind ! type specification parameters
  use symmetry_data_module, only: sym
  use dimensions, only: dim_irr_L => IrrUBasDimSporL, &
                        dim_irr_S => IrrUBasDimSporS
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: back_trafo!()

  !===================================================================
  ! End of public interface of module
  !===================================================================

  !------------ Declaration of constants and variables ---------------

  integer(IK), parameter :: NO_IRREP = -1

  !-------------------------------------------------------------------
  integer(IK), parameter :: LL=1, SS=2, LS=1, SL=2

  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************
  subroutine back_trafo()
    !
    ! Preparations  for A  legacy  way  to perform  DKH  trafo of  the
    ! Coulomb  potential.  Relies on  basis representation  of momentm
    ! vector p. Therefore fragile. Runs in parallel context.
    !
    use error_module, only: error
#if FPP_DEBUG
    use error_module, only: MyID
#endif
    use symmetry_data_module, only: ssym, symmetry_data_get_pcoupling
    implicit none
    ! *** end of interface ***

    integer (IK), allocatable :: pcoupl(:)
    logical, allocatable :: processed(:)
    integer (IK) :: ix, n_irr, irrs(2), dims(2)
    integer (IK) :: memstat

    DPRINT MyID,'back_trafo: entered'

    n_irr = ssym % n_proj_irrep

    allocate (pcoupl(n_irr), processed(n_irr), STAT=memstat)
    call error(memstat,"bt/main: alloc failed")

    pcoupl = symmetry_data_get_pcoupling()
    processed = .false.

    DPRINT '...'
    DPRINT 'bt/main: pcoupling (reduced) is:'
    DWRITE(*,'(10I3)') (ix,ix=1,n_irr)
    DWRITE(*,'(10I3)') pcoupl

    where( pcoupl>n_irr ) pcoupl = NO_IRREP

    DPRINT '...'
    DPRINT 'bt/main: pcoupling (modified) is:'
    DWRITE(*,'(10I3)') (ix,ix=1,n_irr)
    DWRITE(*,'(10I3)') pcoupl

    ix_: do ix = 1,n_irr

       if(processed(ix))cycle

       DPRINT 'bt/main: ...'

       if(ix == pcoupl(ix))then
          !
          ! trivial coupling
          !
          DPRINT 'bt/main: processing irrep ',ix,' no coupling'

          call do_one_block(ix,dim_irr_L(ix))

          processed(ix) = .true.

       else
          !
          !i.e.: coupling is not trivial
          !
          if(pcoupl(ix) /= NO_IRREP)then
             ! basis for small components is not empty
             !
             DPRINT 'bt/main: processing irrep ',ix,', it couples with',pcoupl(ix)

             irrs(1) = ix
             irrs(2) = pcoupl(ix)

             if(pcoupl(irrs(2))/=irrs(1))&
                  & call error("bt/main: strange coupling")

             dims(1) = dim_irr_L(irrs(1))
             dims(2) = dim_irr_S(irrs(1))
             !
             ! dims = dim_irr_L(irrs) should be the same:
             if(any(dims/=dim_irr_L(irrs)))&
                  & call error("bt/main: Leha, your theory is completely wrong")

             call do_two_blocks(irrs,dims)

             processed(irrs(1)) = .true.
             processed(irrs(2)) = .true.

          else
             !
             !i.e.: pcoupl(ix) == NO_IRREP
             !
             DPRINT 'bt/main: processing irrep ',ix,', it couples with an empty irrep'

             irrs(1) = ix
             irrs(2) = NO_IRREP !<<< lousy solution

             dims(1) = dim_irr_L(irrs(1))
             dims(2) = 0        !<<< lousy solution

             call do_two_blocks(irrs,dims)

             processed(irrs(1)) = .true.
          endif
       endif
    enddo ix_

    deallocate(pcoupl,processed,STAT=memstat)
    call error(memstat,"bt/main: dealloc failed")

    DPRINT MyID,'back_trafo: exit'
  end subroutine back_trafo

  subroutine do_two_blocks(irrs,dims)
    use spin_orbit_module, only: is_on,op_FitTrafo
    use back_trafo_tapes, only: put_matrix
    use matrix_module
    implicit none
    integer(IK), intent(in) :: irrs(2),dims(2)
    !*** end of interface ***

    integer(IK)           :: b,irrep1,irrep2,nL,nS

    type(cmatrix), dimension(2) :: UF, SP

    !
    ! WE ARE ALWAYS THINKING OF _ONE_ IRREP
    ! BUT DOING _TWO_ IN PARALLEL
    ! because LL block of one irrep is
    ! the same as SS block of another
    !

    !---------------------------------
    ! store matrices for FitTrafo part:

       irrep1 = irrs(1)
       irrep2 = irrs(2)

       nL = dims(LL)
       nS = dims(SS)

       do b=LL,SS
          call alloc(dims(b), UF(b))
          call get('forward', irrs(b), UF(b))
       enddo

       call alloc(nL, nS, SP(LS))
       call alloc(nS, nL, SP(SL))

       call get('alphap', irrep1, SP(LS))
       call get('alphap', irrep2, SP(SL))
       !
       ! in a simple case following must be done:
       ! SP = tr(UF) * SP * UF
       !
       ! now:
       !
       SP(LS) = tr(UF(LL)) * SP(LS) * UF(SS)
       SP(SL) = tr(UF(SS)) * SP(SL) * UF(LL)

       do b=LS,SL
          if(irrs(b).ne.NO_IRREP)then
             call put_matrix('palphap', irrs(b), SP(b))
          endif
       enddo
  end subroutine do_two_blocks

  subroutine do_one_block(irr,n)
    use spin_orbit_module, only: is_on,op_FitTrafo
    use back_trafo_tapes, only: put_matrix
    use matrix_module
    implicit none
    integer(IK), intent(in) :: irr,n
    !*** end of interface ***

    type(cmatrix) :: UF, SP

    !---------------------------------
    ! store matrices for FitTrafo
       call alloc(n, UF)
       call get('forward', irr, UF)

       call alloc(n, SP)
       call get('alphap', irr, SP) ! LS-block of irr

       SP = tr(UF) * SP * UF ! == sim(SP,UF)

       if(irr.ne.NO_IRREP)then
          call put_matrix('palphap', irr, SP)
       endif
  end subroutine do_one_block

  subroutine get(fn, irr, MAT)
    !
    ! Master  executes  get_matrix, and  bcasts  result.  Executed  in
    ! parallel context. FIXME: make input available on all hosts.
    !
    use matrix_module, only: cmatrix
    use back_trafo_tapes, only: get_matrix
    use comm, only: comm_rank, comm_bcast
    implicit none
    character(len=*), intent(in) :: fn
    integer(IK), intent(in) :: irr
    type(cmatrix), intent(inout) :: MAT ! needs to be allocated
    ! *** end of interface ***

    if (comm_rank() == 0) then
       call get_matrix(fn, irr, MAT)
    endif
    call comm_bcast(MAT%re)
    call comm_bcast(MAT%im)
  end subroutine get

  !--------------- End of module -------------------------------------
end module back_trafo_module

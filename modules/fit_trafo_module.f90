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
module fit_trafo_module
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
  ! Author: SM
  ! Date:   21/04/2004
  ! Description: inclusion of analytic integral
  !
  !----------------------------------------------------------------
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
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  implicit none
  !#save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  public fit_trafo

  !================================================================
  ! End of public interface of module
  !================================================================

  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  real(RK),parameter    ::&
       & half = 0.5_rk

  integer(IK),parameter ::&
       & DKH1  = 1,&
       & oDKH2 = 2,&
       & hDKH2 = 3
  integer(IK)           ::&
       & DKH
  integer(IK),parameter ::&
       & NO_IRREP = -1

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  subroutine init_fit_trafo_module()
    use fit_trafo_tapes
    implicit none
    ! *** end of interface ***

    call init_fit_trafo_tapes()

    call Split3cInt()
  end subroutine init_fit_trafo_module

  subroutine done_fit_trafo_module()
    use fit_trafo_tapes
    implicit none
    ! *** end of interface ***

    call Merge3cInt()

    call done_fit_trafo_tapes()
  end subroutine done_fit_trafo_module

  subroutine fit_trafo()
    use comm, only: comm_barrier
    use back_trafo_module, only: back_trafo
    implicit none
    ! *** end of interface ***

    !
    ! Not parallelized:
    !
    DPRINT MyID, "fit_trafo:: call back_trafo()"
    call back_trafo()
    DPRINT MyID, "fit_trafo: ."

    call init_fit_trafo_module()

    call comm_barrier()

    call main()
    call done_fit_trafo_module()
  end subroutine fit_trafo

  subroutine main()
    use error_module
    use dimensions,&
         & IrrDimL=>IrrUBasDimSporL,&
         & IrrDimS=>IrrUBasDimSporS
    use spin_orbit_module, only:&
         & WhatIs,op_FitTrafo
    use symmetry_data_module
    use spin_orbit_module, only:&
         & is_on,&
         & op_RelFit
    implicit none

    integer (IK), allocatable :: pcoupl(:)
    logical, allocatable :: processed(:)
    integer (IK) :: ix, n_irr, irrs(2)
    integer :: memstat

    n_irr = size (IrrDimL)

    allocate (pcoupl(n_irr), processed(n_irr), STAT=memstat)
    call error(memstat,"ftm/main: alloc failed")

    pcoupl = symmetry_data_get_pcoupling()
    processed = .false.

!!$    print *,'...'
!!$    print *,'ftm/main: pcoupling (reduced) is:'
!!$    write(*,'(10I3)') (ix,ix=1,n_irr)
!!$    write(*,'(10I3)') pcoupl

    where(pcoupl>n_irr)pcoupl = NO_IRREP

!!$    print *,'...'
!!$    print *,'ftm/main: pcoupling (modified) is:'
!!$    write(*,'(10I3)') (ix,ix=1,n_irr)
!!$    write(*,'(10I3)') pcoupl

    !---------------------------------
    ! which transformation? :
    !
    select case(WhatIs(op_FitTrafo))
    case(1)
       DKH = DKH1
       print *,MyID,'ftm/main: DKH is set to DKH(1)'
    case(2)
       DKH = oDKH2
       print *,MyID,'ftm/main: DKH is set to oDKH(2)'
    case(3)
       DKH = hDKH2
       print *,MyID,'ftm/main: DKH is set to hDKH(2)'
    case default
       print *,MyID,'ftm/main: op_FitTrafo ~ ',WhatIs(op_FitTrafo)
       call error("ftm/main: what`s that?")
    end select

    ix_: do ix = 1,n_irr

    Rel_fit : if (is_on(op_RelFit)) then
          DPRINT MyID,'ftm/main: processing irrep ',ix
          call ft_1c_ana(ix)

    else !Rel_fit

       if(processed(ix))cycle

       if(ix == pcoupl(ix))then
          !
          ! trivial coupling
          !
          DPRINT MyID,'ftm/main: processing irrep ',ix, '; no coupling'

          call ft_1b(ix)

          processed(ix) = .true.
       else
          !
          ! i.e.: coupling is not trivial
          !
          if(pcoupl(ix) /= NO_IRREP)then
             !
             ! basis for small components is not empty
             !
             DPRINT MyID,'ftm/main: processing irrep ',ix, ', it couples with',pcoupl(ix)

             irrs(1) = ix
             irrs(2) = pcoupl(ix)

             if(pcoupl(irrs(2))/=irrs(1))&
                  & call error("ftm/main: strange coupling")

             DPRINT MyID,'ftm/main:  call ft_2b(irrs)'
             call ft_2b(irrs)

             processed(irrs(1)) = .true.
             processed(irrs(2)) = .true.

          else
             !
             ! i.e.: pcoupl(ix) == NO_IRREP:
             ! basis for small components IS empty
             !
             DPRINT MyID,'ftm/main: processing irrep ',ix, ', it couples with an empty irrep'

             irrs(1) = ix
             irrs(2) = NO_IRREP !<<< lousy solution

             call ft_2b(irrs)

             processed(irrs(1)) = .true.
          endif
       endif
    endif Rel_fit
    enddo ix_

    deallocate(pcoupl,processed,STAT=memstat)
    call error(memstat,"ftm/main: dealloc failed")
  end subroutine main

  subroutine ft_1c_ana(irr)
    use error_module, only: warn, error
#ifdef FPP_DEBUG
    use error_module, only: MyID
#endif
    use unique_atom_module
    use matrix_module, mm_init => init
    use dimensions, only:&
         & IrrUDim=>IrrUBasDimSpor,&
         & IrrCDim=>IrrBasDimSpor
    use int_send_2cob3c_spor, only:&
         IX_CHFIT, &
         IX_CHFIT_S, &
         IX_CHFIT_R2, &
         MaxFitIndex, &
         MyFitIndex
    use fit_trafo_tapes
    use back_trafo_tapes
    use reltrafo, only: p2_diag
    use spin_orbit_module, only:&
         is_on,whatis,&
         op_NoPVxP,op_NoPFyP,op_RelFit, &
         relfit_minexp, relfit_maxexp, &
         relfit_minZ
    use symmetry_data_module, only: get_totalsymmetric_irrep
    use fit_coeff_module, only: ff_map => fit_coeff_ff_map
    USE_DEBUG
    implicit none
    integer(IK),intent(in) :: irr
    ! *** end of interface ***

    integer(IK)           :: n_ff,ff,n_u,n_c
    integer(IK)           :: n_s,n_r2,r2_ff,s_ff,L,ua

    integer(IK)           :: my_ff

    type(cmatrix)         ::&
         & UF,UB,Ve,&
         & AFA,ARFRA,We_R,R_Wc,wrk,&
         & F_i,F_o,&
         & FAB,FKAB,PVSP,PVXP,PVYP,PFYP

    type(rdmatrix)        ::&
         & p2,Tp,Ep,Ap,Kp,ApKp,K2p2,&
         & invK2p2,wEp

    type(chmatrix)        ::&
         & uFpckd,cFpckd, &
         & PFSSP,PFSXP

    integer(IK) :: doPFyP
    logical     :: do_trafo
    real(RK)    :: aexp, minexp, maxexp
    real(RK)    :: Z, minZ
    character   :: uplo='U'
    !----------------------------------

    doPFyP = whatis(op_NoPFyP)
    minexp = relfit_minexp
    maxexp = relfit_maxexp
    minZ   = relfit_minZ

    print *,'ftm::ana: doPFyP=',doPFyP,' minexp=',minexp,' maxexp=',maxexp,' minZ=',minZ
    select case(doPFyP)
    case (0)
       ! scalar and vector part:
       print *,'ftm: both vector and scalar parta'
    case (1)
       ! no scalar part:
       print *,'ftm: only vector part'
    case (2)
       ! no vector part:
       print *,'ftm: only scalar part'
    case (3)
       ! ignore both
       print *,'ftm: ignore both'
    case default
       ABORT('no such case')
    end select

    n_ff = MaxFitIndex()
    if ( is_on(op_RelFit) ) then
       n_s = MaxFitIndex(IX_CHFIT_S)
       n_r2 = MaxFitIndex(IX_CHFIT_R2)
    endif

    n_u = IrrUDim(irr)
    n_c = IrrCDim(irr)

    call alloc(n_u,UF,UB)
    call alloc(n_u,p2,Tp,Ep,Ap,Kp,ApKp,K2p2)

    call get_matrix('forward' ,irr,UF)
    call get_matrix('backward',irr,UB)

    DPRINT 'ftm: DKH=',DKH
    if(DKH > DKH1)then
       DPRINT 'ftm: get Ve'
       call alloc(n_u,Ve)
       call get_matrix('pnuc',irr,Ve)
    endif

    call get_matrix('p2diag',irr,p2)

    !---------------------------------
    ! prepare some matrices:

    call p2_diag(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)

    if(DKH > DKH1)then
       DPRINT 'ftm: do PVYP'
       call alloc(n_u,PVSP,PVXP,PVYP)

       ! FIXME: re-use matrices:
       call get_hmatrix('pvsp',irr,PVSP)
       call get_hmatrix('pvxp',irr,PVXP)

       if(.not.is_on(op_NoPVxP))then
          ! FIXME: FPP_AIX_XLF as in reltrafo?
#ifdef FPP_AIX_XLF
          PVYP%re = - PVSP%re - PVXP%re
          PVYP%im = - PVSP%im - PVXP%im
#else
          PVYP = - PVSP - PVXP  !<<< V := -V
#endif
       else
#ifdef FPP_AIX_XLF
          PVYP%re = - PVSP%re
          PVYP%im = - PVSP%im
#else
          PVYP = - PVSP
#endif
       endif

       call pgfree(PVSP,PVXP)
       PVYP = sim(PVYP,UF)
    endif

    if(DKH<=DKH1)then
       call alloc(n_u,FAB,FKAB)
       FAB  = UF * Ap * UB    ! F_orward * A  * B_ackward
       FKAB = UF * ApKp * UB  ! F_orward * KA * B_ackward
       ! FIXME: DKH2?
       call pgfree(UF,UB)
    else
       DPRINT 'ftm: do We_R'
       call alloc(n_u,AFA,ARFRA,We_R,R_Wc,wrk)
       call alloc(n_u, wEp)

       invK2p2 = K2P2**(-1)
       wEp%d     = Ep%d * invK2p2%d

       AFA   = Ap * Ve * Ap

     ! ARFRA = AR * (Ve * tr(AR))
       ARFRA = mult(ApKp,PVYP,ApKp)

       call perturb(Ep,AFA,ARFRA)

       ! that is what we need:
       We_R = ARFRA - AFA * K2p2

       call pgfree(Ve,PVYP)
    endif

    ! <<< stop preparations
    !---------------------------------

    call alloc(n_u,F_i)
    call alloc(n_c,F_o)

    call alloc(n_u,uFpckd)
    call alloc(n_c,cFpckd)
    call alloc(n_u,PFSSP)
    call alloc(n_u,PFSXP)
    call alloc(n_u,PFYP)

    call OpenTapes(TH_IJK,irr)

    r2_ff = 0
    s_ff  = 0
    my_ff = 0

    ff_: do ff=1,size(ff_map) ! total N_FF
       if(.not.MyFitIndex(IX_CHFIT,ff) )then
          DPRINT MyID// 'rf: ff=',ff, ' skipping'
          cycle ff_ ! do ff=1,..
       else
          DPRINT MyID// 'rf: ff=',ff,' continue'
       endif
       my_ff = my_ff + 1
       ! orbital momentum of this ff:
       L    = ff_map(ff)%L
       ! "average" exponent of this ff:
       aexp = ff_map(ff)%aexp
       ! location of the ff:
       ua   = ff_map(ff)%ua
       ! charge of the atom ff is located at:
       Z    = ff_map(ff)%Z

       call ReadTape(TH_COUL, irr, uFpckd%re, uFpckd%im)

       ! unpack chmatrix:
       F_i = convert(uFpckd,uplo)

       if (L.eq.-1) then
          s_ff = s_ff+1
          call ReadTape(TH_SS, irr, PFSSP%re, PFSSP%im)
          call ReadTape(TH_SX, irr, PFSXP%re, PFSXP%im)
       elseif (L.eq.0) then
          r2_ff = r2_ff+1
          call ReadTape(TH_RS, irr, PFSSP%re, PFSSP%im)
          call ReadTape(TH_RX, irr, PFSXP%re, PFSXP%im)
       endif

       do_trafo = (L.eq.-1.or.L.eq.0)
       do_trafo = ( aexp > minexp ) .and. do_trafo
       do_trafo = ( aexp < maxexp ) .and. do_trafo
       do_trafo = (   Z >= minZ   ) .and. do_trafo

       if ( do_trafo ) then
          print '(A15,3I4,F20.3," (",F5.1,")")','ftm: TRAFO ',ff,ua,L,aexp,Z

          select case(doPFyP)
          case (0)
             ! scalar and vector part:
             PFYP = convert(PFSSP,uplo)
             PFYP = PFYP + convert(PFSXP,uplo)
          case (1)
             ! no scalar part:
             PFYP = convert(PFSXP,uplo)
          case (2)
             ! no vector part:
             PFYP = convert(PFSSP,uplo)
          case (3)
             ! ignore both
             PFYP%re = 0.0
             PFYP%im = 0.0
          case default
             ABORT('no such case')
          end select

          ! transform it:
          call transform(F_i,PFYP)

       else
          print '(A15,3I4,F20.3," (",F5.1,")")','ftm: skip ',ff,ua,L,aexp,Z
       endif

       call contract_cm(irr,F_i,F_o)

       ! pack chmatrix:
       cFpckd = convert(F_o,uplo)

       call WriteTape(irr, cFpckd%re, cFpckd%im)

       ! report how much is done:
       !call trace(my_ff,n_ff)
    enddo ff_

    call CloseTapes(TH_IJK,irr)

    !---------------------------------
    ! clean up >>>

    if(DKH<=DKH1)then
       call pgfree(FAB,FKAB)
    else
       DPRINT 'ftm: clean DKH2'
       call pgfree(AFA,ARFRA,We_R,R_Wc,wrk)
       call pgfree(invK2p2,wEp)
       call pgfree(UF,UB)
    endif

    call pgfree(F_i,F_o)

    call pgfree(uFpckd)
    call pgfree(cFpckd)
    call pgfree(PFSSP)
    call pgfree(PFSXP)
    call pgfree(PFYP)

    call pgfree(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)
  contains

    subroutine transform(F,PFYP)
      implicit none
      type(cmatrix),intent(inout) :: F
      type(cmatrix),intent(inout) :: PFYP
      ! *** end of interface ***

      if(DKH.eq.DKH1)then
         !    F = tr(FAB) * F * FAB + tr(FRAB) * F * FRAB
!!$         ! transform in p-space:
!!$         F = Ap * F * Ap + AR * (F * tr(AR))
!!$         F = Ap * F * Ap + mult(ApKp,PFYP,ApKp)

         F = tr(FAB) * F * FAB + tr(FKAB) * PFYP * FKAB

      else
         ! to momentum space:
         F    = tr(UF) *  F   * UF
         PFYP = tr(UF) * PFYP * UF

         ! i.e. if DKH(2)
         AFA   = Ap * F * Ap
         ! ARFRA = AR * (F * tr(AR))
         ARFRA = mult(ApKp,PFYP,ApKp)

         F = AFA + ARFRA

         ! F -> F_tilda:
         call perturb(Ep,AFA,ARFRA)

         ! R * W_c:
         R_Wc = K2p2* AFA - ARFRA

         select case(DKH)
         case ( hDKH2 )
            !---------------------------------
            ! linear in F expansion of hamiltonian:
            ! W_e * W_c:
            ! we insert 1 = R*R/K2p2 identity
            wrk = We_R * invK2p2 * R_Wc

            ! {W_e,W_c} = W_e * W_c + W_c * W_e:
            wrk = wrk + tr(wrk)

            ! 1/2 {W_e,W_c}*Ep:
            wrk = half * wrk * Ep
            !---------------------------------
         case ( oDKH2 )
            !---------------------------------
            ! direct transformation of F:
            ! W_e * W_c * Ep:
            ! we insert 1 = R*R/K2p2 identity
            wrk = We_R * invK2p2 * R_Wc * Ep
            !---------------------------------
         case default
            ABORT('no such case')
         end select

         ! + W_e * Ep * W_c:
         wrk = wrk + We_R * wEp * R_Wc !<<< sign?

         ! F = AFA + ARFRA - 1/2 {{W_e,W_c},Ep} - W_e*Ep*W_c - W_c*Ep*W_e:
         F = F - (wrk + tr(wrk)) !<<< sign?

         ! back to real space:
         F    = tr(UB) *  F   * UB
      endif
    end subroutine transform
  end subroutine ft_1c_ana

  subroutine ft_1b(irr)
    use error_module, only: warn, error
    use matrix_module
    use dimensions, only:&
         & IrrUDim=>IrrUBasDimSpor,&
         & IrrCDim=>IrrBasDimSpor
    use int_send_2cob3c_spor, only:&
         & MaxFitIndex
    use fit_trafo_tapes
    use back_trafo_tapes
    use reltrafo, only: p2_diag
    implicit none
    integer(IK),intent(in) :: irr
    ! *** end of interface ***

    integer(IK)           :: n_ff,ff,n_u,n_c

    type(cmatrix)         ::&
         & UF,UB,SP,Ve,&
         & AR,AFA,ARFRA,We_R,R_Wc,wrk,&
         & F_i,F_o,&
         & FAB,FRAB

    type(rdmatrix)        ::&
         & p2,Tp,Ep,Ap,Kp,ApKp,K2p2,&
         & invK2p2,wEp

    type(chmatrix)        ::&
         & uFpckd,cFpckd

    n_ff = MaxFitIndex()

    n_u = IrrUDim(irr)
    n_c = IrrCDim(irr)

    call alloc(n_u,UF,UB)
    call alloc(n_u,p2,Tp,Ep,Ap,Kp,ApKp,K2p2)

    ! this one is also square in this case:
    call alloc(n_u,SP)

    call get_matrix('forward' ,irr,UF)
    call get_matrix('backward',irr,UB)
    call get_matrix('palphap' ,irr,SP)

    if(DKH > DKH1)then
       call alloc(n_u,Ve)
       call get_matrix('pnuc',irr,Ve)
    endif

    call get_matrix('p2diag',irr,p2)

    !---------------------------------
    ! prepare some matices:

    call p2_diag(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)

    call alloc(n_u,AR)

    AR  = ApKp * SP

    call pgfree(SP)

    if(DKH<=DKH1)then
       call alloc(n_u,FAB,FRAB)
       FRAB = UF * tr(AR) * UB
       call pgfree(AR)
       FAB  = UF * Ap * UB
       call pgfree(UF,UB)
    else
       call alloc(n_u,AFA,ARFRA,We_R,R_Wc,wrk)
       call alloc(n_u, wEp)

       invK2p2 = K2P2**(-1)
       wEp%d     = Ep%d * invK2p2%d

       AFA   = Ap * Ve * Ap

       ARFRA = AR * (Ve * tr(AR))

       call perturb(Ep,AFA,ARFRA)

       ! that is what we need:
       We_R = ARFRA - AFA * K2p2

       call pgfree(Ve)
    endif

    ! <<< stop preparations
    !---------------------------------

    call alloc(n_u,F_i)
    call alloc(n_c,F_o)

    call alloc(n_u,uFpckd)
    call alloc(n_c,cFpckd)

    call OpenTapes(TH_IJK,irr)

    ff_:do ff=1,n_ff
       call ReadTape(TH_COUL, irr, uFpckd%re, uFpckd%im)

       ! unpack chmatrix:
       F_i = convert(uFpckd)

       ! transform it:
       call transform(F_i)

       call contract_cm(irr,F_i,F_o)

       ! pack chmatrix:
       cFpckd = convert(F_o)

       call WriteTape(irr, cFpckd%re, cFpckd%im)

       ! report how much is done:
       call trace(ff,n_ff)
    enddo ff_

    call CloseTapes(TH_IJK,irr)

    !---------------------------------
    ! clean up >>>

    if(DKH<=DKH1)then
       call pgfree(FAB,FRAB)
    else
       call pgfree(AR)
       call pgfree(UF,UB)
       call pgfree(AFA,ARFRA,We_R,R_Wc,wrk)
       call pgfree(invK2p2,wEp)
    endif

    call pgfree(F_i,F_o)

    call pgfree(uFpckd)
    call pgfree(cFpckd)

    call pgfree(p2,Tp,Ep,Ap,Kp,ApKp,K2p2)
  contains

    subroutine transform(F)
      implicit none
      type(cmatrix),intent(inout) :: F
      ! *** end of interface ***

      if(DKH.eq.DKH1)then
         F = tr(FAB) * F * FAB + tr(FRAB) * F * FRAB
!!$         ! transform in p-space:
!!$         F = Ap * F * Ap + AR * (F * tr(AR))
      else

         ! to momentum space:
         F = tr(UF) * F * UF

         ! i.e. if DKH(2)
         AFA   = Ap * F * Ap
         ARFRA = AR * (F * tr(AR))

         F = AFA + ARFRA

         ! F -> F_tilda:
         call perturb(Ep,AFA,ARFRA)

         ! R * W_c:
         R_Wc = K2p2* AFA - ARFRA

         select case(DKH)
         case ( hDKH2 )
            !---------------------------------
            ! linear in F expansion of hamiltonian:
            ! W_e * W_c:
            ! we insert 1 = R*R/K2p2 identity
            wrk = We_R * invK2p2 * R_Wc

            ! {W_e,W_c} = W_e * W_c + W_c * W_e:
            wrk = wrk + tr(wrk)

            ! 1/2 {W_e,W_c}*Ep:
            wrk = half * wrk * Ep
            !---------------------------------
         case ( oDKH2 )
            !---------------------------------
            ! direct transformation of F:
            ! W_e * W_c * Ep:
            ! we insert 1 = R*R/K2p2 identity
            wrk = We_R * invK2p2 * R_Wc * Ep
            !---------------------------------
         case default
            ABORT('no such case')
         end select

         ! + W_e * Ep * W_c:
         wrk = wrk + We_R * wEp * R_Wc !<<< sign?

         ! F = AFA + ARFRA - 1/2 {{W_e,W_c},Ep} - W_e*Ep*W_c - W_c*Ep*W_e:
         F = F - (wrk + tr(wrk)) !<<< sign?

         ! back to real space:
         F = tr(UB) * F * UB
      endif
    end subroutine transform
  end subroutine ft_1b

  subroutine ft_2b(irrs)
    use error_module, only: warn, error
    use matrix_module
    use dimensions, only:&
         & IrrUDimL=>IrrUBasDimSporL,&
         & IrrCDimL=>IrrBasDimSporL
    use int_send_2cob3c_spor, only:&
         & MaxFitIndex
    use fit_trafo_tapes
    use back_trafo_tapes
    use reltrafo, only: p2_diag
    implicit none
    integer(IK),intent(in)      :: irrs(2)
    ! *** end of interface ***

    integer(IK),parameter       :: LL=1,SS=2,LS=1,SL=2
    integer(IK)                 :: b,z
    integer(IK)                 :: udims(2),cdims(2),nL,nS
    integer(IK)                 :: irr,n_ff,ff,n_u,n_c

    type(cmatrix),dimension(1:2) ::&
         & UF,UB,SP,&
         & AR,AFA,ARFRA,Ve,We_R,R_Wc,wrk,&
         & F_i,F_o,&
         & BAF,BARF

    type(chmatrix),dimension(1:2) ::&
         & uFpckd,cFpckd

    type(rdmatrix),dimension(1:2) ::&
         & p2,Tp,Ep,Ap,Kp,ApKp,K2p2,&
         & invK2p2,wEp

    n_ff = MaxFitIndex()

    do b=1,2
       if(irrs(b).ne.NO_IRREP)then
          udims(b) = IrrUDimL(irrs(b))
          cdims(b) = IrrCDimL(irrs(b))
       else
          udims(b) = 0
          cdims(b) = 0
       endif
    enddo

    do b=LL,SS
       n_u = udims(b)

       call alloc(n_u,UF(b),UB(b))
       call alloc(n_u,p2(b),Tp(b),Ep(b),Ap(b),Kp(b),ApKp(b),K2p2(b))

       if(DKH > DKH1)then
          call alloc(n_u,Ve(b))
       endif
    enddo

    nL = udims(LL)
    nS = udims(SS)

    call alloc(nL,nS,SP(LS))
    call alloc(nS,nL,SP(SL))

    call alloc(nL,nS,AR(LS))
    call alloc(nS,nL,AR(SL))

    do b=1,2
       irr = irrs(b)
       call get_matrix('forward' ,irr,UF(b))
       call get_matrix('backward',irr,UB(b))
       call get_matrix('palphap' ,irr,SP(b))
       call get_matrix('p2diag'  ,irr,p2(b))

       if(DKH > DKH1)then
          call get_matrix('pnuc',irr,Ve(b))
       endif
    enddo

    do b=LL,SS
       call p2_diag(p2(b),Tp(b),Ep(b),Ap(b),Kp(b),ApKp(b),K2p2(b))
    enddo

    do b=1,2
       ! LS = LL * LS:
       AR(b) = ApKp(b) * SP(b)
       call pgfree(SP(b))
    enddo

    if(DKH<=DKH1)then
       do b=LL,SS
          call alloc(udims(b),BAF(b))
          BAF(b) = tr(UF(b) * Ap(b) * UB(b))
       enddo

       call alloc(nL,nS,BARF(LS))
       call alloc(nS,nL,BARF(SL))

       ! LS = LL * LS * SS:
       BARF(LS) = tr( UF(SS) * tr(AR(LS)) * UB(LL) )
       BARF(SL) = tr( UF(LL) * tr(AR(SL)) * UB(SS) )

       do b=1,2
          call pgfree(UF(b),UB(b),AR(b))
       enddo
    else
       do b=1,2
          n_u = udims(b)

          call alloc(n_u,AFA(b),ARFRA(b))
          call alloc(n_u,We_R(b),R_Wc(b))
          call alloc(n_u, wEp(b))

          invK2p2(b) = K2P2(b)**(-1)
          wEp(b)%d     = Ep(b)%d * invK2p2(b)%d

          if(b.eq.LL) z = SS
          if(b.eq.SS) z = LL

          AFA(b)   = Ap(b) * Ve(b) * Ap(b)

          ARFRA(b) = AR(b) * (Ve(z) * tr(AR(b)))

          call perturb(Ep(b),AFA(b),ARFRA(b))

          ! that is what we need:
          We_R(b) = ARFRA(b) - AFA(b) * K2p2(b)
       enddo

       do b=1,2
          call pgfree(Ve(b))
       enddo
    endif
    do b=LL,SS
       call alloc(udims(b),wrk(b))
    enddo
    ! <<< stop preparations
    !---------------------------------

    do b=LL,SS
       n_u = udims(b)
       n_c = cdims(b)

       call alloc(n_u,F_i(b))
       call alloc(n_c,F_o(b))

       call alloc(n_u,uFpckd(b))
       call alloc(n_c,cFpckd(b))
    enddo

    do b=1,2
       call OpenTapes(TH_IJK,irrs(b))
    enddo

    ff_:do ff=1,n_ff

       do b=LL,SS
          call ReadTape(TH_COUL, irrs(b), uFpckd(b)%re, uFpckd(b)%im)

          ! unpack chmatrix:
          F_i(b) = convert(uFpckd(b))
       enddo

       ! transform it:
       call transform(F_i)

       do b=LL,SS
          call contract_cm(irrs(b),F_i(b),F_o(b))

          ! pack chmatrix:
          cFpckd(b) = convert(F_o(b))

          call WriteTape(irrs(b), cFpckd(b)%re, cFpckd(b)%im)
       enddo

       ! report how much is done:
       call trace(ff,n_ff)
    enddo ff_

    do b=1,2
       call CloseTapes(TH_IJK,irrs(b))
    enddo

    !---------------------------------
    ! clean up >>>

    do b=1,2
       call pgfree(F_i(b),F_o(b))

       call pgfree(p2(b),Tp(b),Ep(b),Ap(b),Kp(b),ApKp(b),K2p2(b))

       call pgfree(uFpckd(b))
       call pgfree(cFpckd(b))

       if(DKH<=DKH1)then
          call pgfree(BAF(b),BARF(b))
       else
          call pgfree(UF(b),UB(b),AR(b))
          call pgfree(AFA(b),ARFRA(b))
          call pgfree(We_R(b),R_Wc(b))
          call pgfree(invK2p2(b),wEp(b))
       endif
       call pgfree(wrk(b))
    enddo
  contains

    subroutine transform(F)
      implicit none
      type(cmatrix),intent(inout) :: F(2)
      ! *** end of interface ***

      integer(IK) :: x,z

      if(DKH.eq.DKH1)then
         !
         ! DKH(1) transformation:
         !
         wrk(LL) = BAF(LL) * (F(LL) * tr(BAF(LL)))&
              &  + BARF(LS) * (F(SS) * tr(BARF(LS)))

         wrk(SS) = BAF(SS) * (F(SS) * tr(BAF(SS)))&
              &  + BARF(SL) * (F(LL) * tr(BARF(SL)))

         do x=LL,SS
            F(x) = wrk(x)
         enddo
      else
         ! to momentum space:
         do x=LL,SS
            F(x) = tr(UF(x)) * F(x) * UF(x)
         enddo
         !
         ! xDKH(2) transformation:
         !
         do x=LL,SS
            if(x.eq.LL) z = SS
            if(x.eq.SS) z = LL

            AFA(x)   = Ap(x) * F(x) * Ap(x)

            ARFRA(x) = AR(x) * (F(z) * tr(AR(x)))
         enddo

         do x=LL,SS

            F(x) = AFA(x) + ARFRA(x)

            ! F -> F_tilda:
            call perturb(Ep(x),AFA(x),ARFRA(x))

            ! R * W_c:
            R_Wc(x) = K2p2(x)* AFA(x) - ARFRA(x)

            select case(DKH)
            case ( hDKH2 )
               !---------------------------------
               ! linear in F hamiltonian expansion:
               ! W_e * W_c:
               ! we insert 1 = R*R/K2p2 identity
               wrk(x) = We_R(x) * invK2p2(x) * R_Wc(x)

               ! {W_e,W_c} = W_e * W_c + W_c * W_e:
               wrk(x) = wrk(x) + tr(wrk(x))

               ! 1/2 {W_e,W_c}*Ep:
               wrk(x) = half * wrk(x) * Ep(x)
               !---------------------------------
            case ( oDKH2 )
               !---------------------------------
               ! direct transformation of F:
               ! W_e * W_c * Ep:
               ! we insert 1 = R*R/K2p2 identity
               wrk(x) = We_R(x) * invK2p2(x) * R_Wc(x) * Ep(x)
               !---------------------------------
            case default
               ABORT('no such case')
            end select

            ! + W_e * Ep * W_c:
            wrk(x) = wrk(x) + We_R(x) * wEp(x) * R_Wc(x) !<<< sign?

            ! F = AFA + ARFRA - 1/2 {(W_e,W_c),Ep} - W_e*Ep*W_c - W_c*Ep*W_e:
            F(x) = F(x) - (wrk(x) + tr(wrk(x))) !<<< sign?
         enddo

         ! back to real space:
         do x=LL,SS
            F(x) = tr(UB(x)) * F(x) * UB(x)
         enddo
      endif
    end subroutine transform
  end subroutine ft_2b

  subroutine perturb(Ep,AVA,ARVRA)
    use matrix_module
    implicit none
    type(rdmatrix),intent(in)    :: Ep
    type(cmatrix), intent(inout) :: AVA,ARVRA
    ! *** end of interface ***

    integer(IK) :: i,j
    real(RK)    :: DeltaE

    do i = 1, size(Ep%d)
       do j = 1, size(Ep%d)
          DeltaE = Ep%d(j)+Ep%d(i)
          AVA%re(j,i)   = AVA%re(j,i)/DeltaE
          AVA%im(j,i)   = AVA%im(j,i)/DeltaE

          ARVRA%re(j,i) = ARVRA%re(j,i)/DeltaE
          ARVRA%im(j,i) = ARVRA%im(j,i)/DeltaE
       enddo
    enddo
  end subroutine perturb

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

  subroutine trace(ff,n_ff)
    use error_module, only: MyID
    implicit none
    integer(IK),intent(in) :: ff,n_ff
    ! *** end of interface ***

    integer(IK),parameter :: n_msg=10,many=1

    if(n_ff/n_msg >= many)then
       if(mod(ff,n_ff/n_msg).eq.0.or.ff.eq.n_ff)then
          write(*,'(A,A,I4,A)') MyID,'ftm/trace: done ',(100*ff)/n_ff,' %'
       endif
    endif
  end subroutine trace

  !--------------- End of module ----------------------------------
end module fit_trafo_module

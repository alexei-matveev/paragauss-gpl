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
module s2_expect
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
  use type_module, only: &
       & IK=>i4_kind, &
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public s2_calc

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  !*************************************************************
  subroutine s2_calc()
    !
    ! Runs on all workers.
    !
    use overlap_module, only: overlap
    use eigen_data_module, only: eigvec
    use occupation_module, only: occ_num
    use comm, only: comm_rank
    implicit none
    !** End of interface *****************************************

    integer (IK) :: n

    ! MOVED TO BEFORE SCF LOOP: call read_overlap()
    if (comm_rank() == 0) then
       n = size (occ_num)

       ! Does not seem to be allocated on slaves:
       ASSERT(allocated(overlap))

       ASSERT(n==size(overlap))
       ASSERT(n==size(eigvec))

       ! FIXME:  O(N^3) serial  step of  computing  alpha/beta overlap
       ! here:
       call s2_expectation (occ_num, overlap, eigvec)
    endif
    ! DONE AFTER hamiltonian_shutdown(): call dealloc_overlap()
  end subroutine s2_calc
  !*************************************************************


  !*************************************************************
  subroutine s2_expectation (occ, ovl, ev)
    !
    ! Computes  the  expectation value  of  the  S2  according to  the
    ! formula:
    !
    !     S2 = S(S+1) + ( Nb - SUM(ja,jb) |Sab(ja,jb)|^2 )
    !
    ! with Na  >= Nb, S =  (Na - Nb)/2,  and Sab(ja, jb) which  is the
    ! overlap of two occupied orbitals with different spins
    !
    ! Does no communication but writes  to output. Runs on master only
    ! so far. FIXME: O(N^3) here!
    !
    use datatype, only: arrmat2, arrmat3
    use iounitadmin_module, only: output_unit
    implicit none
    type(arrmat2), intent(in) :: occ(:) ! (n_irr)
    type(arrmat2), intent(in) :: ovl(:) ! (n_irr)
    type(arrmat3), intent(in) :: ev(:) ! (n_irr)
    !** End of interface *****************************************

    real(RK), parameter :: one=1.0_RK, two=2.0_RK
    integer(IK) :: irr, n_irr
    integer(IK) :: iorb,ilev,n_orb,ia,ib
    integer(IK), dimension(size(occ)) :: oa,ob ! (n_irr)
    integer(IK) :: occa, occb
    real(RK),    dimension(size(occ)) :: na,nb ! (n_irr)
    real(RK)    :: nna, nnb
    integer(IK) :: memstat
    real(RK), allocatable :: S(:,:), EVA(:,:), EVB(:,:)
    real(RK), allocatable :: nela(:), nelb(:)
    real(RK)              :: nel
    real(RK) :: s2_exact, contrib, s2, sz
    logical  :: flip
    real(RK), parameter :: thresh=0.5_rk
    integer(IK) :: ALF=1,BET=2
    !------------ Executable code ------------------------------------


    DPRINT 's2::s2_expectation: entered'
    n_irr = size(occ)

    WRITE(output_unit,'(A)') " =================================================="
    WRITE(output_unit,'(A)') " =           EXPECTATION VALUE OF S^2             ="
    WRITE(output_unit,'(A)') " = S(S+1) + N_{minor} - \Sum_{ij} |S^{ab}_{ij}|^2 ="
    WRITE(output_unit,'(A)') " =================================================="
    WRITE(output_unit,'(A)')
    WRITE(output_unit,'(" level is considered occupied if nel >= ",F6.1)') thresh

    ! determine minor major spins:
    nna  = 0.0_rk
    nnb  = 0.0_rk
    ALF = 1
    BET = 2
    do irr=1, n_irr
       ASSERT(size(occ(irr)%m,2)==2)
       nna  = nna  + sum(pack(occ(irr)%m(:,ALF), mask=(occ(irr)%m(:,ALF) .ge. thresh) ))
       nnb  = nnb  + sum(pack(occ(irr)%m(:,BET), mask=(occ(irr)%m(:,BET) .ge. thresh) ))
    enddo
    flip = ( nna < nnb )
    if( flip )then
       ALF = 2
       BET = 1
    endif

    ! compute occupations (nab) and number of occupied levels (oab):
    oa = 0
    ob = 0
    na = 0.0_rk
    nb = 0.0_rk
    do irr=1, n_irr
       oa(irr) = oa(irr)  + count( occ(irr)%m(:,ALF) .ge. thresh )
       ob(irr) = ob(irr)  + count( occ(irr)%m(:,BET) .ge. thresh )
       na(irr) = na(irr)  + sum(pack(occ(irr)%m(:,ALF), mask=(occ(irr)%m(:,ALF) .ge. thresh) ))
       nb(irr) = nb(irr)  + sum(pack(occ(irr)%m(:,BET), mask=(occ(irr)%m(:,BET) .ge. thresh) ))
    enddo
    occa = sum(oa)
    occb = sum(ob)
    nna  = sum(na)
    nnb  = sum(nb)

    DPRINT 's2::s2_expectation: UP=',nna,' on ',occa
    DPRINT 's2::s2_expectation: DN=',nnb,' on ',occb

    WRITE(output_unit,'(" A ",F6.1," ELECTRONS ON ",I4," LEVELS")') nna, occa
    WRITE(output_unit,'(" B ",F6.1," ELECTRONS ON ",I4," LEVELS")') nnb, occb
    WRITE(output_unit,'(" A ELECTRONS(IRR)  ",20F6.1)') na
    WRITE(output_unit,'(" B ELECTRONS(IRR)  ",20F6.1)') nb
    WRITE(output_unit,'(" A    LEVELS(IRR)",  20I6)'  ) oa
    WRITE(output_unit,'(" B    LEVELS(IRR)",  20I6)'  ) ob

    ! trivial contribution from the A-B difference:
    sz =  (nna-nnb)/two
    s2_exact = sz * (sz + one)

    DPRINT 's2::s2_expectation: s2(exact)=',s2_exact

    WRITE(output_unit,'(A)') " INITIAL GUESS AND IRREP CONTRIBUTIONS ARE:"
    WRITE(output_unit,'(" S(S+1)   ",     F10.6)') s2_exact

    s2 = s2_exact

    ! for each irrep compute contributions due to AB overlap:
    do irr=1,n_irr
       occa = oa(irr)
       occb = ob(irr)
       n_orb = size(ev(irr)%m,1)
       DPRINT 's2::s2_expectation: irr=',irr,' a=',occa,' b=',occb
       ASSERT(occa<=size(ev(irr)%m,2))
       ASSERT(occb<=size(ev(irr)%m,2))

       allocate( S(occa,occb), EVA(n_orb,occa), EVB(n_orb,occb), STAT=memstat )
       ASSERT(memstat==0)
       allocate( nela(occa), nelb(occb), STAT=memstat )
       ASSERT(memstat==0)

       ! copy ALFA egenvectors and occupation nums:
       ilev = 0
       do iorb=1,size(ev(irr)%m,1)
          nel = occ(irr)%m(iorb,ALF)
          if( nel .lt. thresh ) cycle
          ilev = ilev + 1
          EVA(:,ilev) = ev(irr)%m(:,iorb,ALF)
          nela(ilev)  = nel
       enddo
       ASSERT(ilev==occa)

       ! copy BETA eigenvectors and occupation nums:
       ilev = 0
       do iorb=1,size(ev(irr)%m,1)
          nel = occ(irr)%m(iorb,BET)
          if( nel .lt. thresh ) cycle
          ilev = ilev + 1
          EVB(:,ilev) = ev(irr)%m(:,iorb,BET)
          nelb(ilev)  = nel
       enddo
       ASSERT(ilev==occb)

       DPRINT 'NELA=',nela
       DPRINT 'NELB=',nelb

       if( ANY(nela.NE.REAL(NINT(nela),RK)) )then
          WRITE(output_unit,'(" WARNING: fractional ALPHA occ: ",8F6.1,/:,(8F6.1))') nela
       endif
       if( ANY(nelb.NE.REAL(NINT(nelb),RK)) )then
          WRITE(output_unit,'(" WARNING: fractional BETA  occ: ",8F6.1,/:,(8F6.1))') nelb
       endif

       ! Overlap of eigenfunctions with different spin
       ! orientations
       ! S = Ca^T * OVL * Cb :
       S = matmul( transpose(EVA), matmul( ovl(irr)%m, EVB ) )

       ! debug>>>
       DCALL show(S)

       S = S**2
       ! take "degeneracies" into account:
       ! SUM(ij) S2(ij) = SUM(kl) n(k) S2(k,l) n(l)
       ! partners considered orthogonal: S2(p1,p2) ~ delta(p1,p2)
       ! i assume if there are 2 electrons on the 3-fold p level
       ! then e.g. px and py are occupied and pz is unocc.
       ! (better suggestions?)
       do ia=1,occa
          do ib=1,occb
             S(ia,ib) = min(nela(ia),nelb(ib)) * S(ia,ib)
          enddo
       enddo
       contrib = nb(irr) - sum(S)

       DPRINT 's2::s2_expectation: contrib(',irr,')=',contrib

       WRITE(output_unit,'(" IRREP(",I2,")",F10.6)') irr,contrib

       s2 = s2 + contrib

       deallocate(S,EVA,EVB,nela,nelb,STAT=memstat)
       ASSERT(memstat==0)
    enddo

    DPRINT 's2::s2_expectation: s2(real)=',s2

    WRITE(output_unit,'(A)') " FINAL:"
    WRITE(output_unit,'(" < S^2 >  ",     F10.6)') s2
    WRITE(output_unit,'(A)') " =================================================="
  end subroutine s2_expectation
  !*************************************************************

#ifdef FPP_DEBUG
  subroutine show(MM)
    real(RK),dimension(:,:),intent(in) :: MM
    ! *** end of interface ***

    integer(IK),parameter       :: many=7
    integer(IK)           :: i,j,n,m,d_n,d_m

    n = size(MM,1)
    m = size(MM,2)

    d_n = 1
    d_m = 1
    if(n>many) d_n = max(n/many,2)
    if(m>many) d_m = max(m/many,2)

    DPRINT's2::show: shape=',n,'x',m

    ! print a block>>>
    write(*,'(A)') '\matrix>>>'
    write(*,'(8I16)') (i,i=1,m,d_m)
    do i=1,n,d_n
       write(*,'(I3," ",8(E15.8," "))')&
            & i, (MM(i,j), j=1,m,d_m)
    enddo
  end subroutine show
#endif


  !--------------- End of module -------------------------------------
end module s2_expect

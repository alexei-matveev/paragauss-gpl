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
module shgi4_driver
  !-------------------------------------------------------------------
  !
  ! Driver for two-electron integrals. Unfinished, unused.
  !
  ! Copyright (c) 2010 Alexei Matveev
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
! define FPP_TIMERS 1
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  type :: array1 ! isomorph to arrmat in datatype.f90
    real(RK), allocatable :: m(:)
  end type

  type :: array2 ! isomorph to arrmat2 in datatype.f90
    real(RK), allocatable :: m(:, :)
  end type

  public :: array1, array2 ! for callers to prepare input

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: shgi4_direct

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  FPP_TIMER_DECL(tot)
  FPP_TIMER_DECL(erd)

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi4_direct(ibasis, iposition, lvalues, exponents, contractions, positions)
    !
    ! Shell structure is represented by indexing into arrays parametrizing
    ! bases and positions of atoms. Specifically, for a shell number "ishl"
    !
    !   ibas = ibasis(ishl)
    !
    ! is an index into arrays (lvalues, exponents, contractions)
    ! to access L-value, exponents and contractions for the shell "ishl"
    !
    !   L = lvalues(ibas)
    !   E = exponents(ibas)%m(:)
    !   C = contractions(ibas)%m(:, :)
    !
    ! Similarly,
    !
    !   ipos = iposition(ishl)
    !
    ! is an index into array of positions:
    !
    !   X = positions(:, ipos)
    !
    use shgi4_eri, only: shgi4_eri_batch
    use erd, only: erd_batch
    use comm, only: comm_rank, comm_size, comm_reduce
    implicit none
    integer(IK),  intent(in) :: ibasis(:), iposition(:) ! (nshl)
    integer(IK),  intent(in) :: lvalues(:)              ! (nbas)
    type(array1), intent(in) :: exponents(:)            ! (nbas)%m(nexp)
    type(array2), intent(in) :: contractions(:)         ! (nbas)%m(nexp, ncnt)
    real(RK),     intent(in) :: positions(:, :)         ! (3, natm)
    ! *** end of interface ***

    !$ integer, external  :: omp_get_num_threads
    !$ integer(IK)        :: chunk

    integer(IK)           :: rank, np

    integer(IK)           :: ia, ib, ic, id, ip, iq, ipq
    integer(IK)           :: ibasa, ibasb, ibasc, ibasd
    integer(IK)           :: iposa, iposb, iposc, iposd
    integer(IK)           :: la, lb, lc, ld
    integer(IK)           :: nea, neb, nec, ned
    integer(IK)           :: nints
    integer(IK)           :: nshells, npairs, nquarts
    real(RK), allocatable :: batch(:)
!   real(RK)              :: diff ! for debugging
    real(RK)              :: chksum

    ! wrappers around MPI or just placeholders:
    rank = comm_rank()
    np = comm_size()

    ! number of shells:
    nshells = size(ibasis)
    ASSERT(nshells==size(iposition))

    ! four-byte long number of shell quartets may overflow
    ! with (not so) large number of shells:
    ASSERT(nshells<=300)
    !
    !         361 * (361 + 1) / 2 == 65 341
    !
    !     65341 * (65341 + 1) / 2 == 2 134 755 811
    !
    !                        2^16 == 65 536
    !
    !                        2^31 == 2 147 483 648
    !

    ! number of shell pairs:
    npairs = nshells * ( nshells + 1 ) / 2

    ! number of shell quartets:
    nquarts = npairs * ( npairs + 1 ) / 2

    ! for debugging:
    chksum = 0.0

    !
    ! An intuitive way to loop over unique brackets (ia, ib | ic, id):
    !
    !     ip = 0
    !     do ia = 1, nshells ! A-shell counter
    !       do ib = 1, ia    ! B-shell counter
    !         ! only unique bras: (ia, ib |
    !
    !         ip = ip + 1 ! AB-shell pair counter
    !         ASSERT(ip==ia*(ia-1)/2+ib)
    !
    !         iq = 0
    !         do ic = 1, nshells ! C-shell counter
    !           do id = 1, ic    ! D-shell counter
    !             ! only unique kets: | ic, id)
    !
    !             iq = iq + 1 ! CD-shell pair counter
    !             ASSERT(iq==ic*(ic-1)/2+id)
    !
    !             if( iq > ip ) cycle ! only unique bra-kets (ia, ib | ic, id)
    !
    !             ...
    !
    !           enddo
    !         enddo
    !       enddo
    !     enddo
    !

    !
    ! However we will do it differently. To allow for "omp parallel do"
    ! we will manually collapse the four loops into one.
    !

    !
    ! Single loop over (unique) brackets (ia, ib | ic, id):
    !

    ! with OMP you need firstprivate attribute if you want initialization to be
    ! preserved:
    ia = 1
    ib = 1
    ic = 1
    id = 1
    ip = 1
    iq = 1

    !$
    !$ chunk = max(nquarts / omp_get_num_threads() / 1000, 1)
    !$ print *, "OMP chunk size =", chunk, "of", nquarts
    !$
    !$omp parallel do reduction(+: chksum),        &
    !$omp firstprivate(ip, iq, ia, ib, ic, id),    &
    !$omp private( nints, batch                    &
    !$omp        , ibasa, la, nea, iposa           &
    !$omp        , ibasb, lb, neb, iposb           &
    !$omp        , ibasc, lc, nec, iposc           &
    !$omp        , ibasd, ld, ned, iposd ),        &
    !$omp schedule(dynamic, chunk)
    do ipq = 1, nquarts
      ! in MPI-version, decide if this quartet is scheduled for the current
      ! process:
      if ( MOD(ipq - 1, np) /= rank ) cycle ! round robin

      !
      ! Update shell pair and shell indices:
      !
      !   ipq -> (ip, iq), iq <= ip
      !   iq  -> (ic, id), id <= ic
      !   ip  -> (ia, ib), ib <= ia
      !
      call moveto(ipq, ip, iq)
      call moveto(iq, ic, id)
      call moveto(ip, ia, ib)

      ASSERT(ip==ia*(ia-1)/2+ib)
      ASSERT(iq==ic*(ic-1)/2+id)
      ASSERT(ipq==ip*(ip-1)/2+iq)

      ! get basis set indices for four shells:
      ibasa = ibasis(ia)
      ibasb = ibasis(ib)
      ibasc = ibasis(ic)
      ibasd = ibasis(id)

      ! get angular momenta for four shells:
      la = lvalues(ibasa)
      lb = lvalues(ibasb)
      lc = lvalues(ibasc)
      ld = lvalues(ibasd)

      ! basis set size (number of exponents) for four shells:
      nea = size(exponents(ibasa)%m)
      neb = size(exponents(ibasb)%m)
      nec = size(exponents(ibasc)%m)
      ned = size(exponents(ibasd)%m)

      ! get center indices for four shells:
      iposa = iposition(ia)
      iposb = iposition(ib)
      iposc = iposition(ic)
      iposd = iposition(id)

      ! number of integrals for this quartet:
      nints = (2 * la + 1) * nea &
            * (2 * lb + 1) * neb &
            * (2 * lc + 1) * nec &
            * (2 * ld + 1) * ned

      !
      ! MEMORY: (dd|dd) integral with 16 exponents takes
      !         41M doubles or 312MB:
      !
      allocate(batch(nints))

!100   format("Quad = (", 2I3, " |", 2I3, "), LLLL=", 4I3, I6, " of ", I6)
!      write(*, 100) ia, ib, ic, id, la, lb, lc, ld, ipq, nquarts

#if 1
      FPP_TIMER_START(tot)
      call shgi4_eri_batch( la, lb, lc, ld      &
                          , positions(:, iposa) &
                          , positions(:, iposb) &
                          , positions(:, iposc) &
                          , positions(:, iposd) &
                          , exponents(ibasa)%m  &
                          , exponents(ibasb)%m  &
                          , exponents(ibasc)%m  &
                          , exponents(ibasd)%m  &
                          , batch, nints        &
                          )
      FPP_TIMER_STOP(tot)

      ! for debugging:
      chksum = chksum + sum(batch)

!     write(*,'("(", 2I3, " |", 2I3, ") = ", 2G10.4)') &
!              la, lb, lc, ld        &
!            , FPP_TIMER_SLICE(tot)  &
!            , FPP_TIMER_VALUE(tot)

#else
      FPP_TIMER_START(erd)
      ! FIXME: LIBERD may add ints to batch with fake minus sign for comparison,
      !        consult your sources. In such case the content of the batch prior
      !        to entry does matter:
      ! batch = 0.0
      call erd_batch(   la, lb, lc, ld        &
                      , positions(:, iposa)   &
                      , positions(:, iposb)   &
                      , positions(:, iposc)   &
                      , positions(:, iposd)   &
                      , exponents(ibasa)%m    &
                      , exponents(ibasb)%m    &
                      , exponents(ibasc)%m    &
                      , exponents(ibasd)%m    &
                      , contractions(ibasa)%m &
                      , contractions(ibasb)%m &
                      , contractions(ibasc)%m &
                      , contractions(ibasd)%m &
                      , 0                     & ! do not contract
                      , batch, nints          &
                      )
      ! LIBERD returns nints == 0 and possibly junk instead of integrals if
      ! the selection rules tell that all integrals vanish.
      FPP_TIMER_STOP(erd)

      ! for debugging:
      if (nints > 0) then
        chksum = chksum + sum(batch)
      endif

!     write(*,'("(", 2I3, " |", 2I3, ") = ", 2G10.4)') &
!              la, lb, lc, ld       &
!            , FPP_TIMER_SLICE(erd) &
!            , FPP_TIMER_VALUE(erd)
#endif

!     diff = maxval(abs(batch))
!     print *,'int diff maxval(abs(batch))=',diff
!     ASSERT(diff<1.0E-10)

      deallocate(batch)

    enddo ! ipq
    !$omp end parallel do

    ! reduce partial check sums on master:
    call comm_reduce(chksum)

    if ( rank == 0 ) then
      ! print checksums, for debugging only:
      print *, "chksum =", chksum
    endif

  contains

    subroutine moveto(N, I, J)
      !
      ! Updates (I, J), J <= I, so that
      !
      !      I * (I - 1)
      !  N = ----------- + J
      !           2
      !
      ! holds.
      !
      ! FIXME: a better, more direct, implementation?
      !
      implicit none
      integer(IK), intent(in)    :: N
      integer(IK), intent(inout) :: I, J
      ! *** end of interface ***

      integer(IK) :: K

      !   K = 0
      !   IJ: do I = 1, N ! N is very conservatve upper limit
      !     do J = 1, I
      !       K = K + 1
      !       if ( K == N ) exit IJ
      !     enddo
      !   enddo IJ

      !
      ! Find X so that X * (X - 1) / 2 + X < N,
      ! start with X = I - 1, store X at I:
      !
      I = I - 1
      K = I * (I + 1) / 2
      do while ( K >= N )
        K = K - I
        I = I - 1
      enddo
!     ASSERT(I*(I-1)/2+I<N)

      !
      ! It may be that the final I = X + 1,
      ! unless X is too low.
      !

      !
      ! Find Y so that Y * (Y - 1) / 2 + 1 > N,
      ! start with Y = X + 2, store Y at I:
      !
      I = I + 2
      K = I * (I - 1) / 2 + 1
      do while ( K <= N )
        K = K + I
        I = I + 1
      enddo
!     ASSERT(I*(I-1)/2+1>N)

      I = I - 1
      J = N - K + I + 1

!     ASSERT(N==I*(I-1)/2+J)
!     ASSERT(J<=I)
    end subroutine moveto

  end subroutine shgi4_direct

#if 0
  subroutine display( ls, ns, ints)
    implicit none
    integer(IK), intent(in)  :: ls(4), ns(4)
    real(RK),    intent(in)  :: &
    ints( ( 2 * ls(1) + 1) * ns(1), ( 2 * ls(2) + 1) * ns(2) &
        , ( 2 * ls(3) + 1) * ns(3), ( 2 * ls(4) + 1) * ns(4) )
    ! *** end of interface ***

    integer :: i, j, k, l

    do l=1, ( 2 * ls(4) + 1) * ns(4)
    do k=1, ( 2 * ls(3) + 1) * ns(3)
    do j=1, ( 2 * ls(2) + 1) * ns(2)
    do i=1, ( 2 * ls(1) + 1) * ns(1)

      write(*, 1000) ' [', i, j, '|', k, l, '] = ', ints(i, j, k, l)

    enddo
    enddo
    enddo
    enddo
    1000 format (A2, 2I3, A1, 2I3, A4, F20.10)
  end subroutine display

  subroutine DIVMOD3(N, D, I, J, K, L)
    !
    ! Compute (I, J, K, L) such that;
    !
    !         N = ((I * D + J) * D + K) * D + L
    !
    ! e.g.
    !
    !         DIMOD3( 7503, 10, i, j, k ,l )
    !
    ! returns
    !
    !         i, j, k, l == 7, 5, 0, 3
    !
    implicit none
    integer(IK), intent(in)  :: N, D
    integer(IK), intent(out) :: I, J, K, L
    ! *** end of interface ***

    integer(IK) :: a, b

    call DIVMOD(N, D, a, L)
    call DIVMOD(a, D, b, K)
    call DIVMOD(b, D, I, J)
    ! ASSERT(N==((I*D+J)*D+K)*D+L)
  end subroutine DIVMOD3

  subroutine DIVMOD(N, D, A, B)
    !
    ! A = DIV(N, D), B = MOD(N, D)
    !
    implicit none
    integer(IK), intent(in)  :: N, D
    integer(IK), intent(out) :: A, B
    ! *** end of interface ***

    A = N / D
    B = N - A * D
  end subroutine DIVMOD
#endif

  !--------------- End of module -------------------------------------
end module shgi4_driver

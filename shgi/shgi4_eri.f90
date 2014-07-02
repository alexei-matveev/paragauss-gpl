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
module shgi4_eri
  !-------------------------------------------------------------------
  !
  ! Processing a batch of two-electron integrals. Unfinished, unused.
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

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: shgi4_eri_batch

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  real(RK), parameter :: MAXEXP = 50.0

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  ! a copy from shgi_shr.f90:
  integer(IK)   :: L_,M_
  integer(IK), parameter :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)
  integer(IK), parameter :: mof( (MAXL+1)**2 ) = (/((M_,M_=1,2*L_+1),L_=0,MAXL)/)

  FPP_TIMER_DECL(ovl)
  FPP_TIMER_DECL(gam)
  FPP_TIMER_DECL(f4)
  FPP_TIMER_DECL(f4ab)
  FPP_TIMER_DECL(pr4)

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi4_eri_batch( LA, LB, LC, LD &
                            , xa, xb, xc, xd &
                            , ea, eb, ec, ed &
                            , chksum         )
    !
    ! View the integral storage as an array with four axes:
    !
    implicit none
    integer(IK), intent(in)    :: LA, LB, LC, LD
    real(RK),    intent(in)    :: xa(3), xb(3), xc(3), xd(3)
    real(RK),    intent(in)    :: ea(:), eb(:), ec(:), ed(:)
    real(RK),    intent(inout) :: chksum
!   real(RK),    intent(out) :: batch( (2 * LA + 1) * size(ea) &
!                                    , (2 * LB + 1) * size(eb) &
!                                    , (2 * LC + 1) * size(ec) &
!                                    , (2 * LD + 1) * size(ed) )
!   integer(IK), intent(out) :: nints    ! number of ints
    ! *** end of interface ***

    ! FIXME: we always compute the integrals, even if they are zeroes by summetry:
!   nints = size(batch)

    call shgi4_eri_batch0( LA, LB, LC, LD &
                         , xa, xb, xc, xd &
                         , ea, eb, ec, ed &
                         , chksum         )
  end subroutine shgi4_eri_batch

  subroutine shgi4_eri_batch1( LA, LB, LC, LD &
                             , xa, xb, xc, xd &
                             , ea, eb, ec, ed &
                             , chksum         )
    use constants, only: five, two, pi
    use shgi_rad, only: doIL
    use shgi_shr, only: SHR_D4Fv
    use shgi4_pr4, only: SHR_PR4v
!   use debug, only: countNaN
    implicit none
    integer(IK), intent(in)  :: LA, LB, LC, LD
    real(RK),    intent(in)  :: xa(3), xb(3), xc(3), xd(3)
    real(RK),    intent(in)  :: ea(:), eb(:), ec(:), ed(:)
    real(RK),    intent(out) :: chksum !(:, :, :, :) ! see copy()
    ! batch( (2 * LA + 1) * size(ea), (2 * LB + 1) * size(eb) &
    !      , (2 * LC + 1) * size(ec), (2 * LD + 1) * size(ed) )
    ! *** end of interface ***

    logical               :: cutoff12(size(ea),size(eb))
    logical               :: cutoff34(size(ec),size(ed))
    real(RK)              :: ab(3), ab2, cd(3), cd2
    integer               :: nab, ncd, nabcd
    real(RK), allocatable :: zeta1(:),lambda1(:),norm1(:)!ovrl1(:) ! (nab)
    real(RK), allocatable :: zeta2(:),lambda2(:),norm2(:)!ovrl2(:) ! (ncd)
    real(RK), allocatable :: zeta(:),lambda(:),norm(:)             ! (nabcd)
    real(RK), allocatable :: wa(:),wb(:)                           ! (nab)
    real(RK), allocatable :: wc(:),wd(:)                           ! (ncd)
    real(RK), allocatable :: SAB(:,:,:)                            ! (nab, (LA+1)**2, (LB+1)**2)
    real(RK), allocatable :: SCD(:,:,:)                            ! (ncd, (LC+1)**2, (LD+1)**2)
    real(RK), allocatable :: p(:,:), q(:,:)                        ! (nab,3), (ncd,3)
    real(RK), allocatable :: pq2(:), pq(:,:)                       ! (nabcd), (nabcd,3)
    real(RK), allocatable :: IL(:,:)                               ! (nabcd, 1+LA+LB+LC+LD)
    real(RK), allocatable :: F4(:,:,:,:,:)                         ! (nabcd, (LA+1)**2, ..., (LD+1)**2)
    real(RK), allocatable :: I4(:,:,:,:,:)                         ! (nabcd, 2*LA+1, ..., 2*LD+1)
    integer(IK)           :: ls(4), nm(4), ne(4)
    integer(IK)           :: i, j, ij
    integer(IK)           :: mem
    integer  :: memstat

    ls(1) = LA
    ls(2) = LB
    ls(3) = LC
    ls(4) = LD

    ne(1) = size(ea)
    ne(2) = size(eb)
    ne(3) = size(ec)
    ne(4) = size(ed)

    ! angular multiplicity:
    nm(:) = 2 * ls(:) + 1

    mem = 0

    FPP_TIMER_START(ovl)
    ab   = xa - xb
    ab2  = ab(1)**2 + ab(2)**2 + ab(3)**2

    cd   = xc - xd
    cd2  = cd(1)**2 + cd(2)**2 + cd(3)**2

    call shgi4_cutoff(ea, eb, ab2, MAXEXP, cutoff12, nab)

    call shgi4_cutoff(ec, ed, cd2, MAXEXP, cutoff34, ncd)

    !
    ! MEMORY: for (dd|dd) integral with 16 exponents one needs
    !         16^2 * (2 + 1)^4 = 20736 doubles or 162KB
    !
    ! ============= AB-factors ===============================
    ! FIXME: once we have allocatable arguments ...
    allocate( zeta1(nab), lambda1(nab), norm1(nab) & !, ovrl1(nab) &
            , wa(nab), wb(nab), SAB(nab,(LA+1)**2,(LB+1)**2)   &
            , stat=memstat)
    ASSERT(memstat==0)
    mem = mem + nab * 5 + size(SAB)

    call shgi4_factors(LA, LB, ea, eb, cutoff12 &
                      , lambda1, zeta1 &
                      , norm1 & !, ovrl1   &
                      , wa, wb         &
                      )

    ! Compute AB-overlap ints:
    call shgi4_ovrl(LA, LB, ab, zeta1, norm1, SAB)

!   ASSERT(countNaN(SAB)==0)

    ! ============= CD-factors ===============================
    ! FIXME: once we have allocatable arguments ...
    allocate( zeta2(ncd), lambda2(ncd), norm2(ncd) & !, ovrl2(ncd) &
            , wc(ncd), wd(ncd), SCD(ncd,(LC+1)**2,(LD+1)**2)   &
            , stat=memstat)
    ASSERT(memstat==0)
    mem = mem + ncd * 5 + size(SCD)

    call shgi4_factors(LC, LD, ec, ed, cutoff34 &
                      , lambda2, zeta2 &
                      , norm2 & !, ovrl2   &
                      , wc, wd         &
                      )

    ! Compute CD-overlap ints:
    call shgi4_ovrl(LC, LD, cd, zeta2, norm2, SCD)

!   ASSERT(countNaN(SCD)==0)

    FPP_TIMER_STOP(ovl)

    !
    ! locations of mass centers of gaussian products:
    !

    allocate(p(nab,3),q(ncd,3),stat=memstat)
    ASSERT(memstat==0)

    mem = mem + size(p) + size(q)

    do i=1,3
      p(:,i) = wa(:) * xa(i) + wb(:) * xb(i) ! (nab)
      q(:,i) = wc(:) * xc(i) + wd(:) * xd(i) ! (ncd)
    enddo
!   ASSERT(countNaN(p)==0)
!   ASSERT(countNaN(q)==0)

    !
    ! Now to [ijkl] quartets:
    !

    ! no PQ screening so far:
    nabcd = nab * ncd

    !
    ! MEMORY: for (dd|dd) integral with 16 exponents one needs
    !         storage for 16^4 * (2 + 1)^8 = 430M doubles or
    !         3.3GB:
    !
    ! F4 is huge:
    allocate( zeta(nabcd), lambda(nabcd), norm(nabcd)               &
            , pq(nabcd, 3), pq2(nabcd)                              &
            , IL(nabcd, 1+LA+LB+LC+LD)                              &
            , F4(nabcd, (LA+1)**2, (LB+1)**2, (LC+1)**2, (LD+1)**2) &
            , I4(nabcd, 2*LA+1, 2*LB+1, 2*LC+1, 2*LD+1)             &
            , stat=memstat)
    ASSERT(memstat==0)
    mem = mem + nabcd * ( 3 + 3 + 1 ) + size(F4) + size(IL) + size(I4)

    ij = 0
    do j = 1, ncd
      do i = 1, nab
        ij = ij + 1

        ! a + b + g + d:
        lambda(ij) = lambda1(i) + lambda2(j)

        ! ( a + b ) * ( g + d ) / ( a + b + g + d ):
        zeta(ij)   = lambda1(i) * lambda2(j) / lambda(ij)

        norm(ij)   = 2 * pi**(five/two) / ( lambda1(i) * lambda2(j) * sqrt(lambda(ij)) )

        ! PQ vector, and its square:
        pq(ij,:) = p(i,:) - q(j,:)
        pq2(ij)  = pq(ij,1)**2 + pq(ij,2)**2 + pq(ij,3)**2
      enddo
    enddo

    FPP_TIMER_START(gam)

    ! scalar part:
    call doIL(LA + LB + LC + LD, pq2, zeta, norm, IL)

    FPP_TIMER_STOP(gam)

!   ASSERT(countNaN(IL)==0)

    FPP_TIMER_START(f4)

    ! 4-derivatives wrt PQ-harmonics D(lma)D(lmb)D(lmc)D(lmd) * F:
    call SHR_D4Fv(nabcd, LA, LB, LC, LD, pq, IL, F4) ! compute intensive! 1 of 2

    FPP_TIMER_STOP(f4)

!   ASSERT(countNaN(F4)==0)

    FPP_TIMER_START(f4ab)
    !
    ! One needs to convert derivatives wrt PQ to derivatives wrt A, B, C and D:
    ! D(lma)D(lmb)... F => DA(lma)DB(lmb)... F
    ! But instead of scaling a long FABCD, apply the corresponding factors
    ! to SAB and SCD:
    !
    call shgi4_ovrl_wawb(nab, LA, LB, +wa, +wb, SAB)
    call shgi4_ovrl_wawb(ncd, LC, LD, -wc, -wd, SCD)

    FPP_TIMER_STOP(f4ab)

!   ASSERT(countNaN(F4)==0)

    FPP_TIMER_START(pr4)

    ! couple overlap and coulomb derivatives:
    call SHR_PR4v(nab, ncd, LA, LB, LC, LD, F4, SAB, SCD, I4) ! compute intensive! 2 of 2

    FPP_TIMER_STOP(pr4)

!   ASSERT(countNaN(I4)==0)

    ! copy/unpack result:
!   call copy(nm, ne, cutoff12, cutoff34, I4, batch)
    chksum = chksum + sum(I4)

!   ASSERT(countNaN(batch)==0)

!   print *,'shgi4: nabcd =',nabcd,'=',nab,'*',ncd,' size(F4)=',size(F4),'mem=',mem
  end subroutine shgi4_eri_batch1

  subroutine copy(nm, ne, cutoff12, cutoff34, I4, batch)
    !
    ! Piece of code for eathier addressing into batch(*)
    !
    implicit none
    integer(IK), intent(in) :: nm(4), ne(4)
    logical,     intent(in) :: cutoff12(:, :) ! (ne(1), ne(2))
    logical,     intent(in) :: cutoff34(:, :) ! (ne(3), ne(4))

    real(RK),    intent(in) :: I4(:, :, :, :, :)
    ! I4( n1234, nm(1), nm(2), nm(3), nm(4) )

    real(RK), intent(out) :: batch(:, :, :, :)
    ! batch( nm(1) * ne(1),  nm(2) * ne(2),  nm(3) * ne(3),  nm(4) * ne(4) )

    ! *** end of interface ***

    !
    ! This declaration would exceed the maximum number of array axes (8 > 7):
    !
    ! batch( nm(1), NR(1) , nm(2), NR(2) , nm(3), NR(3) ,  nm(4), NR(4) )
    !

    integer :: ma, mb, mc, md
    integer :: i, j, k, l, ijkl

    ! unpack integrals in (1, 2, 3, 4) (ang, exp) order:
    do md = 1, nm(4)
    do mc = 1, nm(3)
    do mb = 1, nm(2)
    do ma = 1, nm(1)

      ijkl = 0
      do l = 1, ne(4)
      do k = 1, ne(3)
      do j = 1, ne(2)
      do i = 1, ne(1)

        if ( cutoff12(i, j) .and. cutoff34(k, l) ) then
          ijkl = ijkl + 1

          batch( ma + (i - 1) * nm(1) &
               , mb + (j - 1) * nm(2) &
               , mc + (k - 1) * nm(3) &
               , md + (l - 1) * nm(4) &
               ) = I4(ijkl, ma, mb, mc, md)
        else
          batch( ma + (i - 1) * nm(1) &
               , mb + (j - 1) * nm(2) &
               , mc + (k - 1) * nm(3) &
               , md + (l - 1) * nm(4) &
               ) = 0.0
        endif

      enddo
      enddo
      enddo
      enddo

    enddo
    enddo
    enddo
    enddo
  end subroutine copy

  subroutine shgi4_ovrl_wawb(nab, LA, LB, wa, wb, SAB)
    !
    ! Computes overlap ints  SAB(:nab,lma,lmb)
    !
    use constants, only: ONE
    implicit none
    integer(IK), intent(in) :: nab, LA, LB
    real(RK), intent(in)    :: wa(:), wb(:)   ! (nab)
    real(RK), intent(out)   :: SAB(:,:,:)     ! (nab,(LA+1)**2,(LB+1)**2)
    ! *** end of interface ***

    integer(IK) :: lma, lmb
    integer(IK) :: ila, ilb
    integer(IK) :: i, L
    real(RK)    :: waL(nab, 1+LA), wbL(nab, 1+LB)

    ! powers of "wa" and "wb":
    waL(:, 1) = ONE
    do L=2, 1+LA
       waL(:, L) = waL(:, L-1) * wa(:)
    enddo

    wbL(:, 1) = ONE
    do L=2, 1+LB
       wbL(:, L) = wbL(:, L-1) * wb(:)
    enddo

    !
    ! convert derivatives wrt PQ to derivatives wrt A, B, C and D:
    !
    !   D(lma)D(lmb)... F => DA(lma)DB(lmb)... F
    !

    do lmb=1,(LB+1)**2
      ilb=lof(lmb)

    do lma=1,(LA+1)**2
      ila=lof(lma)

      do i=1,nab
        !
        !        a         b         g         d
        ! PQ = ----- A + ----- B - ----- C - ----- D
        !      a + b     a + b     g + d     g + d
        !
        ! so that
        !
        !        d          d
        !       --- = wa * ---,  etc ...
        !       dA         dPQ
        !
        ! and similarly for spherical harmonics:
        !
        !                  l
        !       DA(lm) = wa  * D(lm)
        !
        SAB(i, lma, lmb) = SAB(i, lma, lmb)     &
                         * waL(i, 1 + LA - ila) &
                         * wbL(i, 1 + LB - ilb)
      enddo

    enddo
    enddo
  end subroutine shgi4_ovrl_wawb

  subroutine shgi4_ovrl(LA, LB, xd, zeta, norm, SAB)
    !
    ! Computes overlap ints  SAB(:nab,lma,lmb)
    !
    use shgi_rad, only: doSL
    use shgi_shr, only: SHR_D2Av
    use shgi_dnf, only: doD2S
    implicit none
    integer(IK), intent(in) :: LA,LB
    real(RK), intent(in)    :: xd(3)
    real(RK), intent(in)    :: zeta(:), norm(:)   ! (nab)
    real(RK), intent(out)   :: SAB(:,:,:)          ! (nab,(LA+1)**2,(LB+1)**2)
    ! *** end of interface ***

    real(RK)    :: SL(size(zeta),1+LA+LB)
    real(RK)    :: D2A((LA+1)**2,(LB+1)**2,1+LA+LB)
    integer(IK) :: L
    real(RK)    :: d2

    d2 = xd(1)**2 + xd(2)**2 + xd(3)**2

    ! derivatives of overlap integral:
    call doSL(LA+LB, d2, zeta, SL)

    ! FIXME: add norm arg to doSL:
    do L=1,1+LA+LB
      SL(:,L) = SL(:,L) * norm(:) ! FIXME: why not? * (-1)**(LA + LB)
    enddo
    !
    ! Note also the overall sign change by (-1)^(LA+LB)
    ! which is, strictly speaking, not related to change of
    ! differentiation argument but is rather to account for
    ! the relation for the parameter differentiation:
    !
    !   CLM( r - A ) f( ( r - A )^2 / 2 ) = (-1)^l CLM( d/dA ) f^(l)( ( r - A )^2 / 2 )
    !


    ! compute angular part double derivatives D(lma)D(lmb) S:
    call SHR_D2Av(1, LA, LB, xd, D2A)

    ! couple angular and radial parts, also accounts for
    ! the sign change:
    !
    !   DB(lmb) = (-1)^lb D(lmb)
    !
    call doD2S(LA, LB, SL, D2A, SAB)
  end subroutine shgi4_ovrl

  subroutine shgi4_cutoff(ea, eb, d2, maxexp, cutoff, nab)
    !
    ! Screen integrals by pair overlap
    !
    implicit none
    real(RK),    intent(in)   :: ea(:), eb(:)
    real(RK),    intent(in)   :: d2
    real(RK),    intent(in)   :: maxexp
    logical,     intent(out)  :: cutoff(:,:) ! (NEA,NEB)
    integer(IK), intent(out)  :: nab         ! count(CUTOFF)
    ! *** end of interface ***

    integer(IK) :: i,j
    real(RK)    :: a,b

    ! integral screening (CUTOFF) by overlap > exp(-MAXEXP):
    ! MAXEXP in shgi_cntrl is set by shgi_set_maxexp()
    ! prepare mask for screening integrals:

    nab = 0
    do j=1,size(eb)
       b = eb(j)
       do i=1,size(ea)
          a = ea(i)
          if( ( a*b/(a+b) ) * d2 < maxexp ) then
             cutoff(i,j) =  .true.
             nab = nab + 1
          else
             cutoff(i,j) =  .false.
          endif
       enddo
    enddo
  end subroutine shgi4_cutoff

  subroutine shgi4_factors(LA, LB, ea, eb, cutoff &
                          , lambda, zeta &
                          , norm & !, ovrl   &
                          , wa, wb       &
                          )
    !
    ! FIXME:
    !
    use constants, only: TWO, THREE, FOUR, PI, DFAC
    implicit none
    integer(IK), intent(in) :: LA, LB
    real(RK), intent(in)    :: ea(:), eb(:) ! (nea), (neb)
    logical, intent(in)     :: cutoff(:,:)  ! (nea,neb)
    real(RK), intent(out)   :: lambda(:), zeta(:), norm(:) !, ovrl(:) ! (nab)
    real(RK), intent(out)   :: wa(:), wb(:)       ! wx(nab)
    ! *** end of interface ***

    integer(IK) :: i, j, ij
    real(RK)    :: a, b

    ij = 0
    do j=1,size(eb)
       do i=1,size(ea)
          if(.not.cutoff(i,j)) cycle
          ij = ij + 1

          b = eb(j)
          a = ea(i)

          lambda(ij) = a + b
          zeta(ij)   = a * b / lambda(ij)

          ! common pre-factor for norms:
          norm(ij)  = ((2 * a) / PI)**(THREE/FOUR) * ((2 * b ) / PI)**(THREE/FOUR) &
                    / sqrt( a**LA * DFAC(LA) * b**LB * DFAC(LB)    )
!         ! norm for overlap:
!         ovrl(ij)  = ( PI / lambda(ij) )**(THREE/TWO) * norm(ij)

          ! Factors for conversion D(lma)D(lmb) D(lmc) -> DA(lma)DB(lmb) D(lmc)

          wa(ij) = a / ( a + b )
          wb(ij) = b / ( a + b )
       enddo
    enddo
  end subroutine shgi4_factors

  recursive subroutine shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb, ec, ed, chksum)
    !
    ! Divide (and maybe conquer) the exponent range until it is
    ! sensible to delegate the real processing to
    !
    !         shgi4_eri_batch1()
    !
    implicit none
    integer(IK), intent(in) :: LA, LB, LC, LD
    real(RK), intent(in)    :: xa(3), xb(3), xc(3), xd(3)
    real(RK), intent(in)    :: ea(:), eb(:), ec(:), ed(:)
    real(RK), intent(inout) :: chksum !(:, :, :, :)
    ! batch( (2 * LA + 1) * size(ea), (2 * LB + 1) * size(eb) &
    !      , (2 * LC + 1) * size(ec), (2 * LD + 1) * size(ed) )
    ! *** end of interface ***

    integer(IK), parameter :: AX = 4, BX = 3, CX = 2, DX = 1
    integer(IK)            :: ne(4), ls(4), ns(4), loc(1)
    integer(IK)            :: axis, num
    integer(IK)            :: nx, lx
    integer(IK)            :: n, k
!   integer(IK), save      :: level = 0

    !
    ! maxloc(ne(:)) returns the leftmost position, to divide on the rightmost index
    ! of batch(:, :, :, :) more often, reverse the order,
    ! see values of AX, BX, CX, DX:
    !
    ne(AX) = size(ea)
    ne(BX) = size(eb)
    ne(CX) = size(ec)
    ne(DX) = size(ed)

    ls(AX) = LA
    ls(BX) = LB
    ls(CX) = LC
    ls(DX) = LD

    ! an estimate for the size of temporaries is a product of four numbers:
    ns(:) = ne(:) * (ls(:) + 1)**2
    num = product(ns)

    if ( num <= 1 * 2**17 ) then ! <= 1 MB for double precision

      ! actually compute them:
      call shgi4_eri_batch1(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb, ec, ed, chksum)

    else

      ! the axis to divide:
      loc = maxloc(ne)
      axis = loc(1)

      ! number of exponents and angular momentum for this axis:
      nx = ne(axis)
      lx = ls(axis)

      !
      ! Divide and conquer: (1:nx) -> (1:n) ++ (1+n:nx)
      !
      n = nx / 2
      k = (2 * lx + 1) * n

      ! if memory threshold is too low, this can occur:
      ASSERT(n>0)

!     level = level + 1
!     write(*, '("shgi4: RECLEV", I2, " LLLL=", 4I2, " NNNN=", 4I3, " AXIS", I2, " NEXP", I3, " =", I3, " +", I3)') &
!           level, ls, ne, axis, nx, n, nx - n

      select case(axis)
        case (AX)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea(1 : n ), eb, ec, ed, chksum)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea(1 + n:), eb, ec, ed, chksum)
        case (BX)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb(1 : n ), ec, ed, chksum)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb(1 + n:), ec, ed, chksum)
        case (CX)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb, ec(1 : n ), ed, chksum)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb, ec(1 + n:), ed, chksum)
        case (DX)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb, ec, ed(1 : n ), chksum)
          call shgi4_eri_batch0(LA, LB, LC, LD, xa, xb, xc, xd, ea, eb, ec, ed(1 + n:), chksum)
      end select
!     level = level - 1
    endif
  end subroutine shgi4_eri_batch0

  !--------------- End of module -------------------------------------
end module shgi4_eri

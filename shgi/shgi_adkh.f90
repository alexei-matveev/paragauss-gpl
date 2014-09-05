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
module shgi_adkh
  !-------------------------------------------------------------------
  !
  ! Performs  Atomic  DKH  in  each subspace  with  different  angular
  ! momenta.
  !
  ! Performs  DKH transformations  of primitive  integrals in  so that
  ! electrons   inside  of  heavy   atoms  become   relativistic,  but
  ! inter-atomic forces are not affected.
  !
  ! See also modules/relgrads.f90
  !
  ! Copyright (c) 2006-2013 Alexei Matveev
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

! define FPP_TIMERS 2
# include "def.h"
  use type_module, only: &
      IK => i4_kind,  &
      RK => r8_kind  ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================

! if you want ``compatibility'' mode of ADKH then:
!define ADKH_COMPAT


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------

  public :: shgi_adkh_atom
  public :: shgi_adkh_kin
  public :: shgi_adkh_nuc
  public :: shgi_adkh_gr_kin
  public :: shgi_adkh_gr_nuc
  public :: shgi_adkh_sd_kin
  public :: shgi_adkh_sd_nuc
  public :: shgi_adkh_close!() -- clean up module

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  type shell
    ! to hold relativistic ``contraction'' coefficients
    real(RK), allocatable :: t(:)    ! kinetic energy eigenvalues
    real(RK), allocatable :: UF(:,:) ! momentum eigenvectors
    real(RK), allocatable :: UB(:,:) ! inverse P^-1
    real(RK), allocatable :: A1(:,:) ! A1 = A - 1
    real(RK), allocatable :: B1(:,:) ! B1 = B - 1
  end type
  type(shell), allocatable :: shells(:) ! (n_uaL)
  ! uaL is is a linear index for (ua,L) shells:
  !
  ! uaL = 0
  ! do ua=1,n_ua
  !  do L=1,LMAX(ua)
  !    uaL = uaL + 1

  !------------ Declaration of constants and variables ---------------

  !-------------------------------------------------------------------
  !------------ Subroutines ------------------------------------------
contains

  subroutine shgi_adkh_atom(L, e, Z, rad, sexp, scof, uaL)
    !
    ! compute atomic ints and generate rel contractions
    !
    !------------ Modules used ------------------- ---------------
    use shgi_ab, only: shgi_atomic_ints, shgi_atomic_coul
    implicit none
    integer(IK), intent(in)    :: L
    real(RK)   , intent(in)    :: e(:) ! orbital exponents
    real(RK)   , intent(in)    :: Z    ! nucl charge
    real(RK)   , intent(in)    :: rad  ! nucl radius
    real(RK)   , allocatable   :: sexp(:) ! exponents of s-density fit
    real(RK)   , allocatable   :: scof(:) ! coefficients of s-density fit
    integer(IK), intent(in)    :: uaL ! shell index
    ! *** end of interface ***

    real(RK), dimension(size(e),size(e)) :: S,T
    real(RK), dimension(size(e),size(e)) :: V,O
    real(RK), dimension(size(e),size(e)) :: U,W ! coulomb

    if(ready(uaL)) RETURN

    ! compute non-rel (S,T):
    call shgi_atomic_ints(L, e, S, T, V, O, rad)

    ! scale by nucl charge:
    V = V * Z
    O = O * Z

    if(allocated(sexp).and.allocated(scof))then
      WARN('shgi_adkh_atomic: using sreened potential')
      ! compute the coulomb field of fitted density:
      call shgi_atomic_coul(L, e, sexp, scof, U, W)

      ! screened nuclear field:
      V = V - U
      O = O - W
    endif

    call atom(size(E), S, T, V, O, uaL)
  end subroutine shgi_adkh_atom

  subroutine atom(n,S,T,V,O,uaL)
    !
    ! Generate and store rel-contractions
    !
    !------------ Modules used ------------------- ---------------
    use matrix_module, only: sim, mult
    implicit none
    integer(IK), intent(in)    :: n
    real(RK)   , intent(in)    :: S(n,n)
    real(RK)   , intent(inout) :: T(n,n)
    real(RK)   , intent(inout) :: V(n,n)
    real(RK)   , intent(inout) :: O(n,n)
    integer(IK), intent(in)    :: uaL ! shell index
    ! *** end of interface ***

    real(RK), dimension(n)   :: td,trel
    real(RK), dimension(n,n) :: UF,UB
    real(RK), dimension(n,n) :: VP,OP
    real(RK), dimension(n,n) :: Vrel
    real(RK), dimension(n,n) :: A1, B1 ! proj coeffs, full matrices in general

    call momBas2(n,T,S,td,UF,UB)

    ! to momentum basis (CHANGE SIGN!):
    VP = - sim(V,UF) ! tr(UF) * V * UF
    OP = - sim(O,UF) ! tr(UF) * O * UF

    ! needs td = p2/2:
    call dkh2(n,td,VP,OP,trel,Vrel,A1,B1)
    ! ignore output trel, Vrel ...

    ! obtain the transformation matrix for a sequence of ...
    !   1) trafo to mom-space
    !   2) applying projectors A or B (kinematic factors in fpFW case)
    !   3) trafo back to real space
    ! by similarity transform into the ``rest'' frame:
    A1 = mult( mult( UF, A1 ), UB ) ! A = P * a * P^-1
    B1 = mult( mult( UF, B1 ), UB ) ! B = P * b * P^-1

    ! NOTE: we apply the forward/backward loop only to
    !       the rel-corrections A1 = A - 1 and B1 = B - 1

    ! save for use with off-diagonal elements:
    call load('w', uaL, td, UF=UF, UB=UB, A1=A1, B1=B1)
  end subroutine atom

  subroutine shgi_adkh_kin(S,T,uaL,ubL)
    !
    ! Computes and transforms (S,T) -> (Srel,Trel)
    !
    ! So far Srel /= S only for off-diag
    ! (in atom-shell indices) elements.
    !
    !
    !------------ Modules used ------------------- ---------------
#ifdef FPP_TIMERS
    use shgi_cntrl,  only: FPP_TIMER_VARS(adkh)
#endif
    use shgi_cntrl,  only: NEA,NEB,NAB,LA,LB
    use shgi_cntrl,  only: is_on,IAEQB ! needs to know if A==B
    use shgi_common, only: CUTOFF
    use shgi_common, only: AEXP!,BEXP
    use shgi_ab,     only: shgi_atomic_ints, shgi_kin
    use options_module, only: option!(adkh.inf)
    implicit none
    real(RK)   , intent(out) :: S(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK)   , intent(out) :: T(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    integer(IK), intent(in)  :: uaL, ubL ! shell indices
    ! *** end of interface ***

    integer(IK) :: m,n
    real(RK)    :: SX(NEA,NEB),TX(NEA,NEB) ! unpacked storage

    FPP_TIMER_START(adkh)
    ASSERT(NAB==size(S,1))
    ASSERT(NAB==size(T,1))
    ASSERT(2*LA+1==size(S,2))
    ASSERT(2*LB+1==size(S,3))
    ASSERT(2*LA+1==size(T,2))
    ASSERT(2*LB+1==size(T,3))

    if( is_on(IAEQB) .and. uaL == ubL )then
      ! diagonal in atomic centers (A,B)
      ! and shell indices (uaL,ubL)
      ASSERT(LA==LB)
      ASSERT(NEA==NEB)
      ASSERT(NAB==NEA*NEB)

      ! FIXME: we (possibly) repeat this computation
      !        if it was already ``prepare''d for some <A|x|B>

      ! compute non-rel (S,T):
      call shgi_atomic_ints(LA,AEXP,SX,TX) ! LA==LB, AEXP==BEXP

      if( .not.option('adkh.inf') )then
        ! Legacy case of truncated DKH2 transformation
        ! applied to intra-atomic integrals (compat-mode).

        ! returns TX with rel-corrections...
        call kin1(NEA,SX,TX,uaL) ! NEA==NEB
      else
        ! General (unitary, not truncated) transformation
        ! here of intra-atomic integrals.

        ! now apply rel-contraction coeffs:
        call trafo2(NEA,NEB,NEA*NEB,1,TX,SX,uaL,ubL,'kin')
        ! NOTE: shape of TX and SX is always square, no packing ever!
        !       Also the order of T and S is different!
      endif

      ! ... and spread diagonal:
      do m=1,2*LA+1 ! == 2*LB+1
        ! off-diagonal (zeroes) are left unmodified:
        T(:,:,m) = 0.0
        S(:,:,m) = 0.0
        T(:,m,m) = PACK(TX,CUTOFF)
        S(:,m,m) = PACK(SX,CUTOFF) ! FIXME: remove it, if S(diag) is preserved!
      enddo

    else ! off-diagonal:

      ! compute the non-relativistic (S,T)
      ! by invoking general-purpose sub:
      call shgi_kin(S,T)

      n = (2*LA+1)*(2*LB+1) ! number of integrals

      ! now apply rel-contraction coeffs:
      call trafo2(NEA,NEB,NAB,n,T,S,uaL,ubL,'kin')
      ! NOTE: the order of T and S is different!
    endif
    FPP_TIMER_STOP(adkh)
  end subroutine shgi_adkh_kin

!ifdef ADKH_COMPAT
  subroutine kin1(n,S,T,uaL)
    !
    ! Add rel-corrections to T (and none to S so far)
    !
    !------------ Modules used ------------------- ---------------
    use matrix_module, only: mult
    implicit none
    integer(IK), intent(in)    :: n
    real(RK)   , intent(in)    :: S(n,n)
    real(RK)   , intent(inout) :: T(n,n)
    integer(IK), intent(in)    :: uaL ! shell index
    ! *** end of interface ***

    real(RK), dimension(n)   :: td
    real(RK), dimension(n)   :: e,t1,a1,b1,x1
    real(RK), dimension(n,n) :: UF,UB

    ASSERT(ready(uaL))
    ! kinetic energy eigenvalues were saved
    ! when preparing projection coeffs:

    call load('r', uaL, td, UF=UF, UB=UB)

    ! compute kinematic factors:
    call t_diag(td,e,t1,a1,b1,x1)

    ! Trel = (td + t1), t1 -- rel correction:
    T = mult( mult(transpose(UB),t1) , UB ) + T

    ! S is not changed as DKH is (almost) unitary
  end subroutine kin1
!endif

  subroutine shgi_adkh_nuc(z,x,V,rad,uaL,ubL)
    !
    ! First evaluates, then transforms < A | V(of C) | B >,
    ! Result is added to the total (input/output) V.
    !
    !------------ Modules used ------------------- ---------------
#ifdef FPP_TIMERS
    use shgi_cntrl,  only: FPP_TIMER_VARS(adkh)
#endif
    use shgi_cntrl,  only: NEA,NEB,NAB,LA,LB
    use shgi_cntrl,  only: is_on,IAEQC,IBEQC ! needs to know if A==B==C
    use shgi_common, only: CUTOFF, AEXP!==BEXP
    use shgi_relnuc, only: shgi_rel_nuc
    use shgi_ab,     only: shgi_atomic_ints
    use options_module, only: option!(adkh.inf)
    implicit none
    real(RK), intent(in)    :: z, x(3)
    real(RK), intent(inout) :: V(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK), intent(in)    :: rad
    integer(IK), intent(in) :: uaL,ubL ! shell index for storage only
    ! *** end of interface ***

    real(RK), dimension(NAB,2*LA+1,2*LB+1) :: W,O
    real(RK), dimension(NEA,NEB)           :: SX,TX
    real(RK), dimension(NEA,NEB)           :: WX,OX
    integer(IK) :: m,n

    FPP_TIMER_START(adkh)
    ASSERT(NAB==size(V,1))
    ASSERT(2*LA+1==size(V,2))
    ASSERT(2*LB+1==size(V,3))

    ! needs to know if A==B==C:
    if( is_on(IAEQC,IBEQC) .and. uaL == ubL )then
      ! diagonal in atomic centers (A,B,C)
      ! and shell indices (uaL,ubL): A==B==C, LA==LB
      ASSERT(LA==LB)
      ASSERT(NEA==NEB)
      ASSERT(NAB==NEA*NEB)
      ! then no pack/unpack is necessary!

      ! compute atomic ints:
      call shgi_atomic_ints(LA,AEXP,SX,TX,WX,OX,rad)

      ! scale by actual nuc-charge:
      WX = WX * z
      OX = OX * z

      if( .not.option('adkh.inf') )then
        ! Legacy case of truncated DKH2 transformation
        ! applied to intra-atomic integrals (compat-mode).

        ! apply relativistic corrections to V(A==B==C):
        call nuc1(NEA,WX,OX,uaL)
      else
        ! General (unitary, not truncated) transformation
        ! here of intra-atomic integrals.

        ! apply relativistic corrections to off-diag V:
        call trafo2(NEA,NEB,NEA*NEB,1,WX,OX,uaL,ubL,'nuc')
        ! NOTE: shape of WX and OX is always square, no packing ever!
      endif

      ! add the corrected (diagonal) matrices to the cumulative potential:
      do m=1,2*LA+1
         ! off-diag elements (ma/=mb) are not affected by V(A==B==C) contribution:
         V(:,m,m) = V(:,m,m) + PACK(WX,CUTOFF)
      enddo

    else ! off-diagonal, or C/=A or C/=B

      ! instead of just calling
      !           call shgi_rel_nuc(z,V,rad=rad)
      ! which would mean that in ADKH NUCL(of C/=A or C/=B)
      ! is non-relativistic, we evaluate the contribution,
      ! scale it and add to the total ...

      ! FIXME: transforming each contribution of different
      !        centers C is more expensive than transforming
      !        the final SUM!

      ! in ADKH NUCL(of C/=A or C/=B) is FW-transformed!
      W = 0.0
      O = 0.0
      call shgi_rel_nuc(z,x,W,O,rad=rad)

      n = (2*LA+1)*(2*LB+1) ! number of integrals

      ! apply relativistic corrections to off-diag V:
      call trafo2(NEA,NEB,NAB,n,W,O,uaL,ubL,'nuc')
      V = V + W
    endif
    FPP_TIMER_STOP(adkh)
  end subroutine shgi_adkh_nuc

!ifdef ADKH_COMPAT
  subroutine nuc1(n,V,O,uaL)
    !
    ! Adds relativistic corrections to one-center nuclear
    ! attraction integrals.
    !
    !------------ Modules used ------------------- ---------------
    use matrix_module, only: sim, mult
    implicit none
    integer(IK), intent(in)    :: n
    real(RK)   , intent(inout) :: V(n,n)
    real(RK)   , intent(in)    :: O(n,n)
    integer(IK), intent(in)    :: uaL ! shell index for storage only
    ! *** end of interface ***

    real(RK), dimension(n)   :: td,trel
    real(RK), dimension(n,n) :: UF,UB
    real(RK), dimension(n,n) :: VP,OP
    real(RK), dimension(n,n) :: Vrel

    call load('r', uaL, td, UF=UF, UB=UB)

    ! to momentum basis (CHANGE SIGN!):
    VP = - sim(V,UF) ! tr(UF) * V * UF
    OP = - sim(O,UF) ! tr(UF) * O * UF

    ! needs td = p2/2:
    call dkh2(n,td,VP,OP,trel,Vrel)
    ! returns differential effect of relativity!

    ! back to real space (CHANGE SIGN BACK!):
    VP = - mult( mult( transpose(UB), Vrel), UB )

    ! add rel corrections:
    V = V + VP
  end subroutine nuc1
!endif

  subroutine dkh2(n,t,V,O,Trel,Vrel,A1M,B1M)
    ! returns rel corrections Trel := Trel - Tnrel
    !                         Vrel := Vrel - Vnrel
    use matrix_module, only: mult
    use matrix_functions, only: funm
    use spin_orbit_module, only: c=>speed_of_light
    use options_module, only: option!(adkh.inf)
    USE_DEBUG, only: octave
    implicit none
    integer(IK), intent(in)  :: n
    real(RK)   , intent(in)  :: t(:)
    real(RK)   , intent(in)  :: V(:,:),O(:,:)
    real(RK)   , intent(out) :: Trel(:)
    real(RK)   , intent(out) :: Vrel(:,:)
    real(RK)   , intent(out) :: A1M(:,:), B1M(:,:) ! proj coeffs
    optional :: A1M, B1M
    ! *** end of interface ***

    real(RK), dimension(n)   :: e,t1,a1,b1,x1
    real(RK), dimension(n)   :: p,co,si
    real(RK), dimension(n,n) :: X,Y,Z,W
    real(RK), dimension(n,n) :: ELL,ESS,OSL

    integer(IK) :: i

    ! kinematic factors (some offset by nr-limit 1):
    call t_diag(t,e,t1,a1,b1,x1)

    DCALL octave('e',e)

    ! momentum eigenvalue:
    p = sqrt(2*t) / c ! in rel units!

    ! scaled and renormalized O:
    X  = mult( 1/(p*c), O, 1/(p*c) ) ! == mult( 2/p, O/4/c**2, 2/p )

    ! kinematic factors:
    co = (1+a1)         ! cos = a
    si = (1+b1) * p / 2 ! sin = b * p / 2

    ! even/odd parts of the FW-transformed potential
    Y = ( mult(a1,V,a1) + mult(si,X,si) ) + ( mult(a1,V) + mult(V,a1) ) ! + V omitted!
    Z =   mult(co,X,si) - mult(si,V,co)

    ELL = mult(co,V,co) + mult(si,X,si)
    ESS = mult(co,X,co) + mult(si,V,si)
    OSL = mult(co,X,si) - mult(si,V,co)

    DCALL octave( 'ELL', ELL )
    DCALL octave( 'ESS', ESS )
    DCALL octave( 'OSL', OSL )

    ! NOTE: as a difference of two big values:
    !       Y = mult(co,V,co) + mult(si,X,si) - V

    ! PT-weight O1-coupling by energy gap:
    W = RPT(e,Z) ! W1(p,q) = O1(p,q)/(e(p)+e(q))

    ! re-use O1 for subexpressions:
    Z = mult(Z,W,'tn')
    ! upper-left block of four-commutator [O1,W1]:
    Z = ( Z + transpose(Z) ) / 2

    ! we return only rel-corrections:
    Trel = t1
    Vrel = Y + Z ! + V omitted!

    if(.not.present(A1M)) RETURN
    ! otherwise proceed to projection coeffs ...

    ASSERT(present(B1M))
    ASSERT(n==size(A1M,1))
    ASSERT(n==size(A1M,2))
    ASSERT(n==size(B1M,1))
    ASSERT(n==size(B1M,2))

    ! set the diagonal, zero off-diagonal:
    do i=1,n
      A1M(:,i) = 0
      B1M(:,i) = 0
      A1M(i,i) = a1(i)
      B1M(i,i) = b1(i)
    enddo

    ! if fpFW then RETURN from here, otherwise proceed beyond fpFW ...
    if(option("adkh.fpfw"))then
      ! the inter-atomic integrals will be only fpFW transformed!
      WARN('adkh.fpfw')
      RETURN
    endif
    ! ... proceed beyond fpFW:

    !
    ! A product of
    !
    !      | c  -s |                | z  -w'|
    ! fw = |       |    and    dk = |       |
    !      | s   c |                | w   y |
    !
    ! where e.g.
    !          z = 1 - w'w/2  or  z = sqrt( 1 - w'w )
    ! and
    !          y = 1 - ww'/2  or  y = sqrt( 1 - ww' )
    !
    ! is given by
    !
    !       | cz - sw    - sy - cw'|    | COS ... |
    ! dk2 = |                      | =: |         |
    !       | sz + cw      cy - sw'|    | SIN ... !
    !
    ! In the non-orthogonal basis of the modified Dirac equation
    ! the final unitary matrix is sim-transformed by G^1/2:
    !
    !                               | A  ... |
    ! DK2 = G^-1/2 * dk2 * G^1/2 =: |        |
    !                               | B  ... |
    !
    ! here
    !
    !        A = COS    and    B = (2/p) * SIN
    !
    ! because the matrix G^1/2 is as simple as
    !
    !         | 1   0 |               | 1   0 |
    ! G^1/2 = |       | and  G^-1/2 = |       |
    !         | 0 p/2 |               | 0 2/p |
    !

    ! dump SL-block of generator for analysis:
    DCALL octave("W",W)

    ![[=== variant 1: truncated DKH-2 ===
    ! z-1: dk-renormalization offset by 1:
    Z = - mult( W, W, 'tn' ) / 2 ! + 1 on the diagonal is omitted!

    DCALL octave("Z1",Z)
    !=================================]]]

    ![[=== variant 2: unitary DKH-2 ===
    ! z-1: dk-renormalization offset by 1 as
    !
    !     zz - 1 = SQRT( 1 + 2 z ) - 1
    !
    !            = SQRT( 1 - w'w ) - 1
    !
    Z = funm( 2*Z, sqrt1 ) ! refines (Z-1) of 2nd order DKH

    ! NOTE: funm(M,fun) computes function of a square symmetric matrix
    !       via diagonal representation

    DCALL octave("Z2",Z)
    !===============================]]]

    ![[=== variant 3: exact decoupling ===
    if( option('adkh.inf') )then
      ! needs all blocks of the 4-hamiltonian ( HLL, HSS, HSL==HLS' )
      ! to compute SL-coupling W and LL-renormalizaiton Z (offset by 1):
      call exact(e,OSL,ELL,ESS,Z,W)

      DCALL octave("Z3",Z)
    endif
    !==================================]]]

    ! A1M already contains a - 1 == c - 1:
    A1M = A1M           &       ! a - 1, i.e. fpFW approx
        - mult( si, W ) &       ! 1st order in V
        + mult( co, Z )         ! 2nd and higher orders

    ! B1M already contains b - 1 == (2/p)*s - 1:
    B1M = B1M                 & ! b - 1, i.e. fpFW approx
        + mult( (2/p)*co, W ) & ! 1st order in V
        + mult( (2/p)*si, Z )   ! 2nd and higher orders
  end subroutine dkh2

  subroutine exact(e,O,E1,E2,A1,B)
    !
    ! Solves for X in the four-matrix
    !
    !                  | 1  -X'|
    !         1 + X4 = |       |
    !                  | X   1 |
    !
    ! that block-diagonalizes four-hamiltonian
    !
    !               | H11  H12 |
    !          H4 = |          | ( H4' == H4 )
    !               | H21  H22 |
    ! with
    !          H11 =   e + E1
    !          H22 = - e + E2
    !          H21 =   O
    !
    ! This requires solving Ricotti-like equation
    !
    !    H21 - X * H11 + H22 * X - X * H12 * X = 0
    !
    ! and sets
    !
    !         A = 1 / sqrt( 1 + X'X ) =~ 1 - X'X / 2
    !
    !         B = X * A
    !
    ! ``A'' is the LL-block of the renormalization matrix
    !
    !        A4 = 1 / sqrt( 1 - X4^2 )
    !
    ! The procedure returns the difference
    !
    !        A1 = A - 1
    !
    use matrix_module, only: mult
    use matrix_functions, only: funm
    implicit none
    real(RK), intent(in)  :: e(:)
    real(RK), intent(in)  :: O(:,:)
    real(RK), intent(in)  :: E1(:,:), E2(:,:)
    real(RK), intent(out) :: A1(:,:), B(:,:)
    ! *** end of interface ***

    integer(IK) :: n

    n = size(e)

    ASSERT(n==size(O,1))
    ASSERT(n==size(O,2))
    ASSERT(n==size(E1,1))
    ASSERT(n==size(E1,2))
    ASSERT(n==size(E2,1))
    ASSERT(n==size(E2,2))
    ASSERT(n==size(A1,1))
    ASSERT(n==size(A1,2))
    ASSERT(n==size(B,1))
    ASSERT(n==size(B,2))

    ! find such X to eliminate the coupling O:
    B  = rico(e,O,E1,E2)

    ! estimate offset ``renormalization'': A1 = A - 1  =~ ( 1 - X'X / 2 ) - 1:
    A1 = - mult( B, B, 'tn' ) / 2 ! + 1 omitted

    ! refine by A1 = 1 / sqrt( 1 + B'B ) - 1:
    A1 = funm( 2*A1, invsqrt1 ) ! invsqrt1(x) == 1 / sqrt( 1 - x ) - 1

    ! renormalize B = X * A == X * ( 1 + A1 ):
    B = B + mult( B, A1 )
  end subroutine exact

  function rico(e,O,E1,E2) result(X)
    !
    ! Solves iteratively Ricotty-like equations:
    !
    !      O - X * E11 + E22 * X - X * O'* X = 0
    !
    ! using the diagonal domination of
    !
    !          E11 = + e + E1
    ! and
    !          E22 = - e + E2
    !
    ! by resolving for the l.h.s:
    !
    !      X * e + e * X = O
    !                    - X * E1 + E2 * X
    !                    - X * O' * X
    !
    ! or in subscripted form
    !
    !      X(p,q) * ( e(q) + e(p) ) = <p| O + ... |q>
    !
    use matrix_module, only: mult
    implicit none
    real(RK), intent(in) :: e(:)
    real(RK), intent(in) :: O(:,:)
    real(RK), intent(in) :: E1(:,:), E2(:,:)
    real(RK)             :: X(size(e),size(e)) ! result
    ! *** end of interface ***

    integer(IK) :: n,iter
    real(RK)    :: O1(size(e),size(e))
    real(RK)    :: X1(size(e),size(e))
    real(RK)    :: cond,tol

    n = size(e)

    ASSERT(n==size(O,1))
    ASSERT(n==size(O,2))
    ASSERT(n==size(E1,1))
    ASSERT(n==size(E1,2))
    ASSERT(n==size(E2,1))
    ASSERT(n==size(E2,2))

    ! 1) first-order approximation:
    X = RPT(e,O)

    ! ``condition'' number, should be better small:
    cond = maxval(abs(X))
      print *,'rico( dim=',n,'): cond=',cond

    ASSERT(cond<1)
    if( cond == 0 )then
      ! |X| = 0, no need to iterate, return immediately
      ! happens e.g. for Z=0 in BSSE calculations
      RETURN
    endif
    ASSERT(cond>0)

    tol  = 999.0
    iter = 0

    ! iterate ``tol'' to zero, do at least one full iteration:
    do while( tol/cond > 10*epsilon(tol) .or. iter < 3 )
      iter = iter + 1
      ASSERT(iter<50)

      ! 2) second-order in first iteration:
      O1 = O - mult( X, E1 ) + mult( E2, X )

      if( iter > 1 )then
        ! 3) third-order in third approximation:
        O1 = O1 - mult( X, mult( O, X, 'tn') )
      endif

      X1 = X ! save prev value
      X  = RPT(e,O1)

      tol = maxval(abs(X-X1))
      print *,'rico( dim=',n,'): iter=',iter,'tol=',tol,'tol/cond=',tol/cond
    enddo
  end function rico

  subroutine momBas2(n,T,S,t_diag,UF,UB)
    use matrix_module, only: tr, geigs, mult
    implicit none
    integer(IK), intent(in)  :: n
    real(RK)   , intent(in)  :: T(:,:),S(:,:)
    real(RK)   , intent(out) :: t_diag(:)
    real(RK)   , intent(out) :: UF(:,:),UB(:,:)
    ! *** end of interface ***

    call geigs(T,S,t_diag,UF)
    UB = mult(tr(UF),S)
    !
    ! <<< UF, UB - forward and backward to p-space
    !     t_diag - p2
    !---------------------
  end subroutine momBas2

  subroutine t_diag(t,E,T1,A1,B1,X1)
    ! wrapper for the next
    implicit none
    real(RK), dimension(:), intent(in)  :: t
    real(RK), dimension(:), intent(out) ::   E,T1,A1,B1,X1
    ! *** end of interface ***

    integer(IK) :: i

    do i=1,size(t) ! over t eigenvalues
       call t_diag_s(t(i),E(i),T1(i),A1(i),B1(i),X1(i))
    enddo
  end subroutine t_diag

  subroutine t_diag_s(Tnr,E,T1,A1,B1,X1)
    !
    ! kinematic factors f(t), t=p2/2
    !
    ! in c2-units:
    ! T1 := E - 1 - t
    ! A1 := A - 1
    ! B1 := B - 1
    ! X1 := X - 1
    ! where
    ! E  = sqrt( 1 + 2*t )
    ! A  = sqrt( (E+1)/(2*E) )
    ! B  = 2*A/(E+1)
    ! X  = B/A = (E-1)/t
    !
    use spin_orbit_module, only: c=>speed_of_light
    implicit none
    real(RK), intent(in)     :: Tnr
    real(RK), intent(out)    ::     E,T1,A1,B1,X1
    ! *** end of interface ***

    real(RK) :: t
    real(RK), parameter :: EPS = 1.0E-04_rk

    ! compute in rel units with c == 1:
    t  = Tnr / c**2

    ! VALUES OF RELATIVISTIC KINEMATIC FACTORS:
    if( t > EPS ) then
       E  = sqrt( 1 + 2*t )
       X1 = 2 / ( E + 1 )               ! ratio of B/A and (E-1)/t == Trel/Tnrel
       ! this is *correction* to the energy:
       T1 = t * ( X1 - 1 )              ! E = 1 + x*t
       A1 = sqrt(( E + 1 ) / ( 2 * E ))
       B1 = A1 * X1
       ! now offset them by nr-limit 1:
       A1 = A1 - 1
       B1 = B1 - 1
       X1 = X1 - 1
    else
       X1 =  -       t    /   2 & ! 5E-5
             +(      t**2 /   2 &
             -(   5* t**3 /   8 & ! 6E-13
             +(   7* t**4 /   8 & ! 9E-17
             -(  21* t**5 /  16 ))))

       ! this is *correction* to the energy:
       T1 =          t * X1

       A1 = -        t    /   4 & ! 3E-5
            +(   11* t**2 /  32 &
            -(   69* t**3 / 128 &
            +( 1843* t**4 /2048 & ! 9E-17
            -(12767* t**5 /8192 ))))

       B1 = -     3* t    /   4 & ! 7E-5
            +(   31* t**2 /  32 &
            -(  187* t**3 / 128 &
            +( 4859* t**4 /2048 & ! 2E-16
            -(32965* t**5 /8192 ))))

       E  = (T1 + t) + 1
    endif
    ! back to atomic units with c /= 1:
    E  =  E  * c**2
    T1 =  T1 * c**2
  end subroutine t_diag_s

  function sqrt1(x) result(f)
    !
    ! sqrt1(x) =  sqrt( 1 + x ) - 1 =
    !
    !            2    3      4      5
    !      x   x    x    5 x    7 x
    !   =  - - -- + -- - ---- + ---- + . . .
    !      2   8    16   128    256
    !
    ! ( here called only with x < 0 )
    real(RK), intent(in) :: x
    real(RK)             :: f ! result
    ! *** end of interface ***

    ASSERT(1+x>=0)

    if( abs(x) > 1.0D-4 )then
      f = sqrt( 1 + x ) - 1
    else
      f =     x    / 2   & ! ~ -0.5D-4
        -     x**2 / 8   &
        +     x**3 / 16  &
        - 5 * x**4 / 128 & ! ~  0.04D-16
        + 7 * x**5 / 256   ! ~ -0.03D-20
    endif
  end function sqrt1

  function invsqrt1(x) result(f)
    !
    ! invsqrt1(x) =  1 / sqrt( 1 - x ) - 1 =
    !
    !             2      3       4       5
    !      x   3 x    5 x    35 x    63 x
    !   =  - + ---- + ---- + ----- + ----- + . . .
    !      2    8      16     128     256
    !
    ! ( here called only with x < 0 )
    real(RK), intent(in) :: x
    real(RK)             :: f ! result
    ! *** end of interface ***

    ASSERT(1-x>=0)

    if( abs(x) > 1.0D-4 )then
      f = 1 / sqrt( 1 - x ) - 1
    else
      f =     x    / 2   & ! ~ -0.5D-4
        + 3 * x**2 / 8   &
        + 5 * x**3 / 16  &
        + 35* x**4 / 128 & ! ~  0.3D-16
        + 63* x**5 / 256   ! ~ -0.2D-20
    endif
  end function invsqrt1

  function RPT(e,V) result(VT)
    !
    ! Relativistic ``Perturbation Theory'' Weighting:
    !
    !    VT(p,q) = V(p,q) / ( e(p) + e(q) )
    !
    implicit none
    real(RK)   , intent(in) :: e(:), V(:,:)
    real(RK)                :: VT(size(e),size(e))
    ! *** end of interface ***

    integer(IK) :: p,q,n

    n = size(e)
    ASSERT(n==size(V,1))
    ASSERT(n==size(V,2))
    do q=1,n
       do p=1,n
          VT(p,q) = V(p,q) / ( e(p) + e(q) )
       enddo
    enddo
  end function RPT

  subroutine shgi_adkh_gr_kin(S,T,uaL,ubL)
    !
    ! Transforms gradients of (S,T) -> (Srel,Trel)
    ! Expects NR gradients on input ...
    !
    ! So far Srel /= S only for off-diag
    ! (in atom-shell indices) elements.
    !
    !
    !------------ Modules used ------------------- ---------------
#ifdef FPP_TIMERS
    use shgi_cntrl,  only: FPP_TIMER_VARS(adkh)
#endif
    use shgi_cntrl,  only: NEA,NEB,NAB,LA,LB
    use shgi_cntrl,  only: is_on,IAEQB ! for checking only
    implicit none
    real(RK)   , intent(inout) :: S(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    real(RK)   , intent(inout) :: T(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
                                             ! grads wrt A (and -B)
    integer(IK), intent(in)    :: uaL, ubL   ! shell indices
    ! *** end of interface ***

    integer(IK) :: n

    FPP_TIMER_START(adkh)
    ASSERT(NAB==size(S,1))
    ASSERT(NAB==size(T,1))
    ASSERT(2*LA+1==size(S,2))
    ASSERT(2*LB+1==size(S,3))
    ASSERT(2*LA+1==size(T,2))
    ASSERT(2*LB+1==size(T,3))
    ASSERT(3==size(S,4))
    ASSERT(3==size(T,4))

    if(is_on(IAEQB))then
      ! single-center (A==B) gradients sum up to zero
      ! anyway GA+GB==0 -- dont even compute them:
      ABORT('why calling?')
    endif

    n = (2*LA+1)*(2*LB+1)*3 ! number of integrals

    ! now apply rel-contraction coeffs:
    call trafo2(NEA,NEB,NAB,n,T,S,uaL,ubL,'kin')
    ! NOTE: the order of T and S is different!
    FPP_TIMER_STOP(adkh)
  end subroutine shgi_adkh_gr_kin

  subroutine shgi_adkh_gr_nuc(V,Y,uaL,ubL)
    !
    ! Transforms
    !         V(a,b) = < A |  V  | B >
    ! using SR ints
    !         Y(a,b) = < A | pVp | B >
    !
    ! Result is in (input/output) V.
    !
    !------------ Modules used ------------------- ---------------
#ifdef FPP_TIMERS
    use shgi_cntrl,  only: FPP_TIMER_VARS(adkh)
#endif
    use shgi_cntrl,  only: NEA,NEB,NAB,LA,LB
    use shgi_cntrl,  only: is_on,IAEQC,IBEQC ! needs to know if A==B==C
    implicit none
    real(RK), intent(inout) :: V(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    real(RK), intent(inout) :: Y(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9) not modified
                                          ! grads wrt A, B, and C
    integer(IK), intent(in) :: uaL,ubL ! shell index for storage only
    ! *** end of interface ***

    integer(IK) :: n

    FPP_TIMER_START(adkh)
    ASSERT(NAB==size(V,1))
    ASSERT(2*LA+1==size(V,2))
    ASSERT(2*LB+1==size(V,3))
    ASSERT(9==size(V,4))
    ASSERT(NAB==size(Y,1))
    ASSERT(2*LA+1==size(Y,2))
    ASSERT(2*LB+1==size(Y,3))
    ASSERT(9==size(Y,4))

    if( is_on(IAEQC,IBEQC) )then
      ! gradient sums to zero if A==B==C:
      ABORT('why calling?')
    endif

    n = (2*LA+1)*(2*LB+1)*9 ! number of integrals

    ! now apply rel-contraction coeffs:
    call trafo2(NEA,NEB,NAB,n,V,Y,uaL,ubL,'nuc')
    FPP_TIMER_STOP(adkh)
  end subroutine shgi_adkh_gr_nuc

  subroutine shgi_adkh_sd_kin(S,T,uaL,ubL)
    !
    ! Transforms gradients of (S,T) -> (Srel,Trel)
    ! Expects NR gradients on input ...
    !
    ! So far Srel /= S only for off-diag
    ! (in atom-shell indices) elements.
    !
    !
    !------------ Modules used ------------------- ---------------
#ifdef FPP_TIMERS
    use shgi_cntrl,  only: FPP_TIMER_VARS(adkh)
#endif
    use shgi_cntrl,  only: NEA,NEB,NAB,LA,LB
    use shgi_cntrl,  only: is_on,IAEQB ! for checking only
    implicit none
    real(RK)   , intent(inout) :: S(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    real(RK)   , intent(inout) :: T(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
                                               ! dervs wrt A (and -B)
    integer(IK), intent(in)    :: uaL, ubL     ! shell indices
    ! *** end of interface ***

    integer(IK) :: n

    FPP_TIMER_START(adkh)
    ASSERT(NAB==size(S,1))
    ASSERT(NAB==size(T,1))
    ASSERT(2*LA+1==size(S,2))
    ASSERT(2*LB+1==size(S,3))
    ASSERT(2*LA+1==size(T,2))
    ASSERT(2*LB+1==size(T,3))
    ASSERT(3==size(S,4))
    ASSERT(3==size(T,4))
    ASSERT(3==size(S,5))
    ASSERT(3==size(T,5))

    if(is_on(IAEQB))then
      ! single-center (A==B) gradients sum up to zero
      ! anyway GA+GB==0 -- dont even compute them:
      ABORT('why calling?')
    endif

    n = (2*LA+1)*(2*LB+1)*3*3 ! number of integrals

    ! now apply rel-contraction coeffs:
    call trafo2(NEA,NEB,NAB,n,T,S,uaL,ubL,'kin')
    ! NOTE: the order of T and S is different!
    FPP_TIMER_STOP(adkh)
  end subroutine shgi_adkh_sd_kin

  subroutine shgi_adkh_sd_nuc(V,Y,uaL,ubL)
    !
    ! Transforms
    !         V(a,b) = < A |  V  | B >
    ! using SR ints
    !         Y(a,b) = < A | pVp | B >
    !
    ! Result is in (input/output) V.
    !
    !------------ Modules used ------------------- ---------------
#ifdef FPP_TIMERS
    use shgi_cntrl,  only: FPP_TIMER_VARS(adkh)
#endif
    use shgi_cntrl,  only: NEA,NEB,NAB,LA,LB
    use shgi_cntrl,  only: is_on,IAEQC,IBEQC ! needs to know if A==B==C
    implicit none
    real(RK), intent(inout) :: V(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9,9)
    real(RK), intent(inout) :: Y(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9,9) not modified
                                            ! dervs wrt A, B, and C
    integer(IK), intent(in) :: uaL,ubL ! shell index for storage only
    ! *** end of interface ***

    integer(IK) :: n

    FPP_TIMER_START(adkh)
    ASSERT(NAB==size(V,1))
    ASSERT(2*LA+1==size(V,2))
    ASSERT(2*LB+1==size(V,3))
    ASSERT(9==size(V,4))
    ASSERT(9==size(V,5))
    ASSERT(NAB==size(Y,1))
    ASSERT(2*LA+1==size(Y,2))
    ASSERT(2*LB+1==size(Y,3))
    ASSERT(9==size(Y,4))
    ASSERT(9==size(Y,5))

    if( is_on(IAEQC,IBEQC) )then
      ! gradient sums to zero if A==B==C:
      ABORT('why calling?')
    endif

    n = (2*LA+1)*(2*LB+1)*9*9 ! number of integrals

    ! now apply rel-contraction coeffs:
    call trafo2(NEA,NEB,NAB,n,V,Y,uaL,ubL,'nuc')
    FPP_TIMER_STOP(adkh)
  end subroutine shgi_adkh_sd_nuc

  subroutine trafo2(na,nb,nab,n,TV,SW,uaL,ubL,sub)
    !
    !  applies projection transformation to ints between two shells
    !
    !  case 'kin'
    !      TV==T, SW==S
    !
    !  case 'nuc'
    !      TV==V, SW==pVp (=:W)
    !
    !------------ Modules used ------------------- ---------------
    use matrix_module, only: mult
    use spin_orbit_module, only: c=>speed_of_light
    implicit none
    integer(IK),      intent(in)    :: na,nb
    integer(IK),      intent(in)    :: nab,n
    real(RK)   ,      intent(inout) :: TV(nab,n) ! packed ints
    real(RK)   ,      intent(inout) :: SW(nab,n) ! packed ints
    integer(IK),      intent(in)    :: uaL,ubL ! shell indices
    character(len=3), intent(in)    :: sub ! kin or nuc
    ! *** end of interface ***

    real(RK), dimension(na,na) :: A1a,B1a
    real(RK), dimension(nb,nb) :: A1b,B1b
    integer(IK) :: i
    integer(IK), save :: warncount = 0

    ! check availability of contractions:
    ASSERT(ready(uaL))
    ASSERT(ready(ubL))

    ! get contractions for A and B:
    call load('r', uaL, A1=A1a, B1=B1a)
    call load('r', ubL, A1=A1b, B1=B1b)
    ! NOTE: A1 = A - 1 and B1 = B - 1

!   ASSERT(nab==na*nb)
    if(nab==na*nb)then
      ! no unpack is necessary:
      do i=1,n
        ! direct call to kin() or nuc():
        call dcall(sub,TV(:,i),SW(:,i))
      enddo
    else
      if( warncount <= 10 )then
        WARN('nab/=na*nb in ADKH!')
        if( warncount == 10 )then
          WARN('will ignore further warnings')
        endif
        warncount = warncount + 1
      endif
      do i=1,n
        ! call through pack/unpack:
        call pcall(sub,TV(:,i),SW(:,i))
      enddo
    endif
  contains

    subroutine kin(T,S)
      !
      ! apply fpFW (or beyond) to T and S
      !
      ! Arguments declared with explicit shape!
      ! so that array (nab) becomes (na,nb)
      !
      implicit none
      real(RK), intent(inout) :: T(na,nb)
      real(RK), intent(inout) :: S(na,nb)
      ! *** end of interface ***

      real(RK), dimension(na,nb) :: X,Y

      !
      ! S := a'Sa + b'Tb/2        ( here ' == transpose )
      !
      !    = S + a1'*S + S*a1 + a1'*S*a1
      !     (T + b1'*T + S*b1 + b1'*T*b1) / 2
      !
      ! computing corrections is numerically safer ...

      ! SS-metric block in rel-units:
      X = T / 2 / c**2

      ! correction to <a|S|b> approximately grouped by magnitude:
      Y = ( mult( A1a, S, 'tn' )              &
          + mult( S, A1b, 'nn' )              &
          + X                                 &
          )                                   &
        +(( mult( A1a, mult( S, A1b ), 'tn' ) &
          + mult( B1a, X, 'tn' )              &
          + mult( X, B1b, 'nn' )              &
          )                                   &
         +  mult( B1a, mult( X, B1b ), 'tn' ) &
         )

      !
      ! T := a'Tb + b'Ta - b'Tb == a'Ta - (a-b)'T(a-b)
      !
      !    = T + ( a1'*T + T*a1 ) + ( a1'*T*b1 + b1'*T*a1 ) - b1'*T*b1
      !

      ! correction to <a|T|b> approximately grouped by magnitude:
      X = ( mult( A1a, T, 'tn' )              &
          + mult( T, A1b, 'nn' )              &
          )                                   &
        + ( mult( A1a, mult( T, B1b ), 'tn' ) &
          + mult( B1a, mult( T, A1b ), 'tn' ) &
          - mult( B1a, mult( T, B1b ), 'tn' ) &
          )

      T = T + X
      S = S + Y
    end subroutine kin

    subroutine nuc(V,O)
      ! apply fpFW to V using O
      !
      ! Arguments declared with explicit shape!
      ! so that array (nab) becomes (na,nb)
      !
      implicit none
      real(RK), intent(inout) :: V(na,nb)
      real(RK), intent(inout) :: O(na,nb)
      ! *** end of interface ***

      real(RK) :: X(na,nb)

      ! V := a'Va + b'Ob/4
      !
      !    =   V + a1'*V + V*a1 + a1'*V*a1
      !
      !    + ( O + b1'*O + O*b1 + b1'*O*b1 ) / 4
      !

      ! scale O and transform to ``rel-units'':
      O = O / 4 / c**2

      ! correction to <a|V|b> approximately grouped by magnitude:
      X = ( O                                   &
          + mult( A1a, V, 'tn' )                &
          + mult( V, A1b, 'nn' )                &
          )                                     &
        + ( ( mult( A1a, mult( V, A1b ), 'tn' ) &
            + mult( B1a, O, 'tn' )              &
            + mult( O, B1b, 'nn' )              &
            )                                   &
          + mult( B1a, mult( O, B1b ), 'tn' )   &
          )

      V = V + X
    end subroutine nuc

    subroutine dcall(sub,T,S)
      implicit none
      character(len=3), intent(in)    :: sub ! kin or nuc
      real(RK),         intent(inout) :: T(na,nb)
      real(RK),         intent(inout) :: S(na,nb)
      ! *** end of interface ***

      select case(sub)
      case('kin')
        call kin(T,S)
      case('nuc')
        call nuc(T,S)
      case default
        ABORT('no such sub')
      end select
    end subroutine dcall

    subroutine pcall(sub,T,S)
      ! unpack matrices to rectangular form
      ! and pass through the transformation
      use shgi_common, only: CUTOFF
      implicit none
      character(len=3), intent(in)    :: sub ! kin or nuc
      real(RK),         intent(inout) :: T(nab)
      real(RK),         intent(inout) :: S(nab)
      ! *** end of interface ***

      real(RK) :: X(na,nb),Y(na,nb)

      X = unpack(T,CUTOFF,0.0_rk)
      Y = unpack(S,CUTOFF,0.0_rk)

      call dcall(sub,X,Y)

      T = pack(X,CUTOFF)
!     if( sub == 'nuc' ) RETURN
      S = pack(Y,CUTOFF) ! overlap may be also modified!
    end subroutine pcall
  end subroutine trafo2

  subroutine load(rw, uaL, t, UF, UB, A1, B1)
    !
    ! save atomic data for future use
    !
    use shgi_utils, only: shellNum!(ua,L)
    implicit none
    character  , intent(in)    :: rw ! r or w
    integer(IK), intent(in)    :: uaL
    real(RK)   , intent(inout) :: t(:)
    real(RK)   , intent(inout) :: UF(:,:) ! forward
    real(RK)   , intent(inout) :: UB(:,:) ! backward
    real(RK)   , intent(inout) :: A1(:,:) ! == A - 1
    real(RK)   , intent(inout) :: B1(:,:) ! == B - 1
    optional :: t,UF,UB,A1,B1
    ! *** end of interface ***

    integer(IK), save :: uaL_max
    integer(IK)       :: n

    DPRINT 'shgi_adkh: load(',rw,uaL,')'

    if(.not.allocated(shells))then
      ! on first entry...

      ! count number of shells:
      uaL_max = shellNum(999,999) ! counts all of them
      allocate(shells(uaL_max))
    endif

    ASSERT(uaL>0)
    ASSERT(uaL<=uaL_max)

    if(present(t))then
      n = size(t)
    else
      n = size(shells(uaL)%t)
    endif

    if(present(t))then
      ASSERT(present(UF))
      ASSERT(present(UB))
      ASSERT(n==size(UF,1))
      ASSERT(n==size(UF,2))
      ASSERT(n==size(UB,1))
      ASSERT(n==size(UB,2))
    endif

    if(present(A1))then
      ASSERT(present(B1))
      ASSERT(n==size(A1,1))
      ASSERT(n==size(A1,2))
      ASSERT(n==size(B1,1))
      ASSERT(n==size(B1,2))
    endif

    select case( rw )
    case( 'w' )
      ASSERT(present(t))
      ASSERT(present(UF))
      ASSERT(present(UB))
      ASSERT(present(A1))
      ASSERT(present(B1))
      ASSERT(.not.allocated(shells(uaL)%t))
      ASSERT(.not.allocated(shells(uaL)%UF))
      ASSERT(.not.allocated(shells(uaL)%UB))
      ASSERT(.not.allocated(shells(uaL)%A1))
      ASSERT(.not.allocated(shells(uaL)%B1))

      shells(uaL) = shell(t, UF, UB, A1, B1)

    case( 'r' )
      ASSERT(allocated(shells(uaL)%t))
      ASSERT(allocated(shells(uaL)%UF))
      ASSERT(allocated(shells(uaL)%UB))
      ASSERT(allocated(shells(uaL)%A1))
      ASSERT(allocated(shells(uaL)%B1))

      ASSERT(n==size(shells(uaL)%t))

      if(present(t))  t  = shells(uaL)%t
      if(present(UF)) UF = shells(uaL)%UF
      if(present(UB)) UB = shells(uaL)%UB
      if(present(A1)) A1 = shells(uaL)%A1
      if(present(B1)) B1 = shells(uaL)%B1
    case default
      ABORT('no such case')
    end select
  end subroutine load

  function ready(uaL) result(yes)
    !
    ! returns true if rel contractions for (ua,L) are ready
    !
    implicit none
    integer(IK), intent(in)    :: uaL
    logical                    :: yes ! result
    ! *** end of interface ***

    yes = .false.
    if(.not.allocated(shells)) return

    ASSERT(uaL>0)
    ASSERT(uaL<=size(shells))

    yes = allocated(shells(uaL)%t)
  end function ready

  subroutine shgi_adkh_close()
    !
    ! Idempotent, does not do anything if the module was not used.
    !
    implicit none
    ! *** end of interface ***

    !
    ! Nested deallocation:
    !
    if ( allocated(shells) ) then
       deallocate(shells)
    endif
  end subroutine shgi_adkh_close

  !--------------- End of module -------------------------------------
end module shgi_adkh

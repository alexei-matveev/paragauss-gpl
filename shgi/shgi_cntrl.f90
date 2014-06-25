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
module shgi_cntrl
  !---------------------------------------------------------------
  !
  ! Handling integral settings as a bitmap.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
  ! Copyright (c) 2006 Vladimir Nasluzov
  ! Copyright (c) 2006-2008 Alexey Shor
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
! use CPU_TIME for timers:
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind, & ! type specification parameters
       I8K=>i8_kind
! use constants -- use directly!
  implicit none
  save            ! save all variables defined in this module
  public          ! by default, all names are public
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  integer(I8K), parameter, private :: b=2_i8k ! Binary Base
  integer(I8K), parameter, public ::  &
       INUCL = b**0                 , & ! Nuclear Attraction
       INUSR = b**1                 , & ! Nuclear Attraction: Scalar Relativistic
       INUSO = b**2                 , & ! Nuclear Attraction: Spin-Orbit
       INUXX = INUCL + INUSR + INUSO, & ! Nuclear Attraction: Any
       ICHRG = b**3                 , & ! Charge Fit
       ICHSR = b**4                 , & ! Charge Fit: Scalar Relativistic
       ICHSO = b**5                 , & ! Charge Fit: Spin-Orbit
       ICHXX = ICHRG + ICHSR + ICHSO, & ! Charge Fit: Any
       IPSEU = b**6                 , & ! Pseudopotential
       IFINA = b**7                 , & ! Finite Nucleus
       IXXSR = INUSR + ICHSR        , & ! Any Scalar Relativistic
       IXXSO = INUSO + ICHSO        , & ! Any Spin-Orbit
       INUGR = b**8                 , & ! Nuclear Attraction Gradient
       IPSGR = b**9                 , & ! Pseudoputential Gradient
       ISRGR = b**10                , & ! Scalar Relativistic Gradient
       IXXGR = INUGR + IPSGR + ISRGR, & ! Any Gradient
       INUSD = b**11                , & ! Nuclear Attraction Second Derivatives
       ISRSD = b**12                , & ! Scalar Relativistic Second Derivatives
       IPSSD = b**13                , & ! Pseudoputential Second Derivatives
       IYSNU = b**14                , & ! Angular Part: Nuclear Attraction
       IYSSR = b**15                , & ! Angular Part: Scalar Relativistic
       IYSSO = b**16                , & ! Angular Part: Spin-Orbit
       IGSNU = b**17                , & ! Angular Part: Nuclear Attraction Gradient
       IGSSR = b**18                , & ! Angular Part: Scalar Relativistic Gradient
       ! options for checking if A =? B =? C:
       IAEQB = b**19                , & ! A == B
       IAEQC = b**20                , & ! A == C
       IBEQC = b**21                , & ! B == C
       IXEQC = IAEQC + IBEQC        , & ! A == C OR B == C
       IOVRL = b**22                , & ! overlap ints
       IKNTC = b**23                , & ! kinetic ints
       IOVGR = b**24                , & ! overlap grads
       IKNGR = b**25                , & ! kinetic grads
       IOVSD = b**26                , & ! overlap sec ders
       IKNSD = b**27                , & ! kinetic sec ders
       IPCSD = b**28                , & ! Point Charges: Second Derivatives
       ! AS: to calculate gardients of solvation
       ISLGR = b**29                , & ! Solvation Gradients
       ISLGRT= b**30                , & ! FIXME: Solvation Gradients
       IADKH = b**32                    ! Atomic DKH

  integer(I8K), private :: MODE ! bitmap of those above

  real(RK)               :: MAXEXP = HUGE(1.0_rk) ! use options_integral_expmax() when possible!
  integer(IK), parameter :: R2TYPE=-1, STYPE=0

  integer(IK), parameter :: &
       ! indices for SH:
       ! - C*  for global enumeration
       ! - M*  for local  enumeration within same L
       C1Y =  4, CY = 4, MY = 3, & ! PY
       C1X =  3, CX = 3, MX = 2, & ! PX
       C1Z =  2, CZ = 2, MZ = 1, & ! PZ
       C00 =  1, C0 = 1, M0 = 1, & ! C0 and also C0*C0
       ! indices for products of two SH, first acting
       ! on I, second on S in D(1)D(2) I*S
       ! [ and the above C00=C0*C0 ]        ! for NUCL, etc.
       CR2 =  0, & ! CX*CX + CY*CY + CZ*CZ  ! for SR
       CP2 = -1, & ! C0 * SC(ij)DA(i)DB(j)  ! for SR
       CVZ = -2, & ! CX*CY - CY*CX          ! for SO
       CVY = -3, & ! CZ*CX - CX*CZ          ! for SO
       CVX = -4    ! CY*CZ - CZ*CY          ! for SO

  ! indices for vectors:
  integer(IK), parameter :: &
       VX  = 1, & ! /= CVX
       VY  = 2, & ! /= CVY
       VZ  = 3    ! /= CVZ

  ! use to denote gradients dA_i, dB_i, i=1,2,3
  integer(IK), parameter :: &
       GAX =  1, & ! dAx: WA*CX*CO + C0*CX  ! for GRADS
       GAY =  2, & ! dAy: WA*CY*CO + C0*CY  ! for GRADS
       GAZ =  3, & ! dAz: WA*CZ*CO + C0*CZ  ! for GRADS
       GBX =  4, & ! dBx: WB*CX*CO - C0*CX  ! for GRADS
       GBY =  5, & ! dBy: WB*CY*CO - C0*CY  ! for GRADS
       GBZ =  6, & ! dBz: WB*CZ*CO - C0*CZ  ! for GRADS
       GCX =  7, & ! dCx:    CX*CO          ! for GRADS
       GCY =  8, & ! dCy:    CY*CO          ! for GRADS
       GCZ =  9    ! dCz:    CZ*CO          ! for GRADS

  ! There are two independent gradients of translationally invariant
  ! matrix element <A|B|C>, GD and GW w.r.t. vectors
  !
  !      D = A - B
  !      W = wa * A + wb * B - C
  !
  ! grads w.r.t. A,B, and C are given by
  !
  !     GC =    - GW
  !     GA = wa * GW + GD
  !     GB = wb * GW - GD
  !
  integer(IK), parameter :: &
       GWX =  1, &
       GWY =  2, &
       GWZ =  3, &
       GDX =  4, &
       GDY =  5, &
       GDZ =  6
  ! However, for PP matrix elements that depend rather
  ! on
  !
  !     AC = A - C
  !     BC = B - C
  !
  ! two independent grads, GA and GB, w.r.t. A and B maybe more
  ! natural

  integer(IK)            :: LA,LB
  integer(IK)            :: LC=0 ! FIXME: raise LC=1 for SR/SO and GRADS
  integer(IK)            :: LD=0 ! FIXME: raise LD=1 for SR gradients
  integer(IK)            :: LE=0 ! FIXME: raise LE=1 for SR sec. ders
  integer(IK)            :: LF=0
  integer(IK)            :: LAB  ! max(LA,LB)
  integer(IK)            :: LMAX ! max(LA,LB,LC,LD)
  integer(IK)            :: NEA,NEB,NAB

  integer(IK), parameter :: PC=1
  integer(IK), parameter :: PD=2
  integer(IK), parameter :: PQ=3
  integer(IK), parameter :: PO=4
  integer(IK), parameter :: IPD=5
  integer(IK), parameter :: RC=6
#ifdef WITH_EFP
  logical                :: do_efp_grad=.false.
#endif

  !----------------------------------------------------------------

  ! for collecting statistics:
  integer(IK) :: shgi_stat_num_ints = 0 ! save!
  integer(IK) :: shgi_stat_screened = 0 ! save!

  ! timers:
  FPP_TIMER_DECL(tot)
  FPP_TIMER_DECL(totI)
  FPP_TIMER_DECL(totG)
  FPP_TIMER_DECL(totD)
  FPP_TIMER_DECL(t2c)
  FPP_TIMER_DECL(adkh)
  FPP_TIMER_DECL(tsc)
  FPP_TIMER_DECL(tra)
  FPP_TIMER_DECL(td3a)
  FPP_TIMER_DECL(td3f)
  FPP_TIMER_DECL(tpr2)
  FPP_TIMER_DECL(tyls)
  FPP_TIMER_DECL(tgam)
  FPP_TIMER_DECL(tpq1)
  FPP_TIMER_DECL(tbs1)
  FPP_TIMER_DECL(tpq2)
  FPP_TIMER_DECL(tbs2)
  FPP_TIMER_DECL(tpsc)
  FPP_TIMER_DECL(tra2)
  FPP_TIMER_DECL(pcs)
  FPP_TIMER_DECL(tpcr)
  FPP_TIMER_DECL(pd)
  FPP_TIMER_DECL(tpd)
  FPP_TIMER_DECL(pq)
  FPP_TIMER_DECL(tpq)
  FPP_TIMER_DECL(po)
  FPP_TIMER_DECL(tpo)
  FPP_TIMER_DECL(epot)
  FPP_TIMER_DECL(efld)
  FPP_TIMER_DECL(slgr)
  FPP_TIMER_DECL(slsd)
  !------------ Subroutines ---------------------------------------
contains

  function is_on(OP,OP1,OP2) result(yes)
    ! returns true if ...
    ! - ... any of the bits in OP  are set
    ! - AND any of the bits in OP1 are set
    ! - AND any of the bits in OP2 are set
    ! (of course only if the latter are present)
    implicit none
    integer(I8K), intent(in) :: OP
    integer(I8K), intent(in) :: OP1
    integer(I8K), intent(in) :: OP2
    logical                 :: yes
    optional :: OP1,OP2
    ! *** end of interface ***

    yes = ( IAND(MODE,OP) /= 0 )
    if( present(OP1) )then
       yes = yes .and. ( IAND(MODE,OP1) /= 0 )
    endif
    if( present(OP2) )then
       yes = yes .and. ( IAND(MODE,OP2) /= 0 )
    endif
  end function is_on

  function whatis(OP) result(and)
    implicit none
    integer(I8K), intent(in) :: OP
    integer(I8K)             :: and
    ! *** end of interface ***

    and = IAND(MODE,OP)
  end function whatis

  subroutine setif(OP,val)
    ! set bits to true/false, leave others as is
    implicit none
    integer(I8K), intent(in) :: OP
    logical    , intent(in) :: val
    optional :: val
    ! *** end of interface ***

    logical :: set

    set = .true.
    if( present(val) ) set = val

    ! set bits to true/false, leave others as is
    if( set )then
       ! set OP bits to true:
       MODE = IOR(MODE,OP)
    else
       ! set OP bits to false:
       MODE = IAND(MODE,NOT(OP))
       if(OP==-1)then ! NOT(-1)=0 holds everywhere?
          ASSERT(MODE==0)
       endif
    endif
  end subroutine setif

  subroutine shgi_set_xeqy(UA,EA,UB,EB,UC,EC)
    implicit none
    integer(IK), intent(in) :: UA,EA,UB,EB,UC,EC
    ! *** end of interface ***

    call setif(IAEQB, (UA==UB).and.(EA==EB) )
    call setif(IAEQC, (UA==UC).and.(EA==EC) )
    call setif(IBEQC, (UB==UC).and.(EB==EC) )
    DPRINT UA,EA,UB,EB,UC,EC,'XEQY',is_on(IAEQB),is_on(IAEQC),is_on(IBEQC)
  end subroutine shgi_set_xeqy

  subroutine shgi_set_maxexp(expmax)
    implicit none
    real(RK), intent(in) :: expmax
    ! *** end of interface ***

    MAXEXP = expmax
  end subroutine shgi_set_maxexp

  !--------------- End of module ----------------------------------
end module shgi_cntrl

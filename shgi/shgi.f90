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
module shgi
  !---------------------------------------------------------------
  !
  ! Drivers for one-electron and relativistic density fitting integral
  ! evaluation.  Module   called  by  integral_calc_quad_2cob3c()  and
  ! friends.
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
  use constants
  use shgi_cntrl
  ! FIXME: only for re-export:
  use shgi_slv, only: shgi_gr_solv_drv, shgi_gr_Q_solv_drv, shgi_gr_solv_drv_vtn, &
       shgi_sd_solv_drv
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: shgi_drv_coul2
  public :: shgi_drv
  public :: shgi_gr_drv
  public :: shgi_sd_drv
  ! re-export from shgi_slv:
  public :: shgi_gr_solv_drv
  public :: shgi_gr_solv_drv_vtn
  public :: shgi_gr_Q_solv_drv
  public :: shgi_sd_solv_drv

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  ! ALL INTEGER CONSTANTS, LIKE GAX,GAY,GAZ ARE NOW IN
  !                     shgi_cntrl.f90
  ! THIS HAS BEEN DONE TO ENABLE ITS USE IN OTHER MODULES

  ! MOST GLOBAL VARIABLES HOLDING ``ANGULAR'' FACTORS WERE MOVED TO
  !                     shgi_common.f90
  ! THIS HAS BEEN DONE TO SPLIT THIS FILE INTO PARTS LATER

  ! *** KEEP GLOBALS TO MINUMUM!
  ! *** USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!

  ! temp storage for nuclear attraction:
  real(RK), allocatable :: PINUCL(:,:,:)               ! (NAB,      2*LA+1,2*LB+1  )
  real(RK), allocatable :: PINUSR(:,:,:)               ! (NAB,      2*LA+1,2*LB+1  )
  real(RK), allocatable :: PINUSO(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,3)
  real(RK), allocatable :: PIPSEU(:,:,:)               ! (NAB,      2*LA+1,2*LB+1  )
  real(RK), allocatable :: PIOVRL(:,:,:)               ! (NAB,      2*LA+1,2*LB+1  )
  real(RK), allocatable :: PIKNTC(:,:,:)               ! (NAB,      2*LA+1,2*LB+1  )

  ! for gradients:
  real(RK), allocatable :: GRNUCL(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,6), W,D
  real(RK), allocatable :: GRNUSR(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,6), W,D
  real(RK), allocatable :: GRPSEU(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,9), A,B,C
  real(RK), allocatable :: GROVRL(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,3), A and B=-A
  real(RK), allocatable :: GRKNTC(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,3), A and B=-A

  ! for second derivatives:
  real(RK), allocatable :: SDNUCL(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,6,6)
  real(RK), allocatable :: SDPSEU(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,9,9)
  real(RK), allocatable :: SDNUSR(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,6,6)
  real(RK), allocatable :: SDOVRL(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,3,3)
  real(RK), allocatable :: SDKNTC(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,3,3)

  ! temp storage for charge fit (s-, and r2-):
  real(RK), allocatable :: PICHRG(:,:,:,:)             ! (NAB,NECS ,2*LA+1,2*LB+1  )
  real(RK), allocatable :: PIR2RG(:,:,:,:)             ! (NAB,NECR2,2*LA+1,2*LB+1  )
  real(RK), allocatable :: PICHSR(:,:,:,:)             ! (NAB,NECS ,2*LA+1,2*LB+1  )
  real(RK), allocatable :: PIR2SR(:,:,:,:)             ! (NAB,NECR2,2*LA+1,2*LB+1  )
  real(RK), allocatable :: PICHSO(:,:,:,:,:)           ! (NAB,NECS ,2*LA+1,2*LB+1,3)
  real(RK), allocatable :: PIR2SO(:,:,:,:,:)           ! (NAB,NECR2,2*LA+1,2*LB+1,3)

  ! temp storage for grads wrt (A,B,C), used in ADKH:
  real(RK), allocatable :: GRABC1(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,9), A,B,C
  real(RK), allocatable :: GRABC2(:,:,:,:)             ! (NAB,      2*LA+1,2*LB+1,9), A,B,C
  ! temp storage for dervs wrt (A,B,C), used in ADKH:
  real(RK), allocatable :: SDABC1(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,9,9), A,B,C
  real(RK), allocatable :: SDABC2(:,:,:,:,:)           ! (NAB,      2*LA+1,2*LB+1,9,9), A,B,C

  !----------------------------------------------------------------
  ! a copy from shgi_shr.f90:
  integer(IK)   :: L_,M_
  integer(IK), parameter :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)
  integer(IK), parameter :: mof( (MAXL+1)**2 ) = (/((M_,M_=1,2*L_+1),L_=0,MAXL)/)
  !------------ Subroutines ---------------------------------------
contains

  !****************************************************************
  !****************** DRIVER FOR INTEGRALS ************************
  !****************************************************************

  subroutine shgi_drv(IU1,IE1,IL1,IU2,IE2,IL2, uas, pcs, &
       NUCL, &
       NUSR, &
       NUSO, &
       CHSR, &
       R2SR, &
       CHSO, &
       R2SO, &
       PSEU, &
       OVRL, &
       KNTC, &
       IMOD  &
       )
    use unique_atom_module, only: uat=>unique_atom_type
    use datatype          , only: pct=>pointcharge_type
    use shgi_ang, only: shgi_set_c, shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_ab,  only: shgi_kin
    use shgi_utils , only: upack2c, upack3c, shgi_timing
    use shgi_utils , only: shellNum!(ua,L)
    use shgi_relfit, only: &
         shgi_nuc, shgi_nuc_so, shgi_nuc_sr, &
         shgi_ch , shgi_ch_so , shgi_ch_sr , &
         shgi_r2 , shgi_r2_so , shgi_r2_sr
    use shgi_relnuc, only: shgi_rel_nuc
    use shgi_pseudo, only: &
         shgi_pseu_set_ab, shgi_pseu_set_abc, shgi_pseu_close_abc, shgi_pseu
    use shgi_pcm, only: shgi_pcs
    use shgi_adkh, only: shgi_adkh_kin, shgi_adkh_nuc, shgi_adkh_atom
    use shgi_ext_c, only: shgi_ext
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)  :: uas(:)            ! array of unique atoms, normally all of them
    type(pct)  , intent(in)  :: pcs(:)            ! array of point charges
    real(RK)   , intent(out) :: NUCL(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
    real(RK)   , intent(out) :: NUSR(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
    real(RK)   , intent(out) :: NUSO(:,:,:,:,:)   ! (N2E,N1E,    2*L2+1,2*L1+1,3)
    real(RK)   , intent(out) :: CHSR(:,:,:,:,:)   ! (N2E,N1E,NFF,2*L2+1,2*L1+1  )
    real(RK)   , intent(out) :: R2SR(:,:,:,:,:)   ! (N2E,N1E,NFF,2*L2+1,2*L1+1  )
    real(RK)   , intent(out) :: CHSO(:,:,:,:,:,:) ! (N2E,N1E,NFF,2*L2+1,2*L1+1,3)
    real(RK)   , intent(out) :: R2SO(:,:,:,:,:,:) ! (N2E,N1E,NFF,2*L2+1,2*L1+1,3)
    real(RK)   , intent(out) :: PSEU(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
    real(RK)   , intent(out) :: OVRL(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
    real(RK)   , intent(out) :: KNTC(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
    integer(I8K), intent(in)  :: IMOD
    optional :: pcs
    optional :: NUCL, NUSR, NUSO
    optional :: CHSR, CHSO
    optional :: R2SR, R2SO
    optional :: PSEU
    optional :: OVRL, KNTC
    optional :: IMOD
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------

    integer(IK) :: ua,eq   ! Indices of the third center C.
!   integer(IK) :: xua,xeq ! Just loop variables, DONT USE! Use ``ua'' and ``eq'' instead!
!                          ! --- The iteration order over C may be not obvious ...
    integer(IK) :: ich,ir2,cnt ! ff counters
    real(RK)    :: Z,R     ! charge and radius of the nucleus
    real(RK)    :: xc(3)

    integer(IK) :: LPNUM
    real(RK)    :: ZCORE
    integer(IK) :: npc  ! size(pcs)

    logical     :: diagonal  ! flags if UA==UB, EA==EB, and LA=LB
    integer(IK) :: sameUA    ! ==IU2==IU1 if diagonal, -1 otherwise
    integer(IK) :: sameEA    ! ==IE2==IE1 if diagonal, -1 otherwise
    integer(IK) :: sameL     ! ==IL2==IL1 if diagonal, -1 otherwise
    integer(IK) :: uaL,ubL   ! linear index of (UA,LA) and (OB,LB) shells

    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    DPRINT 'SHGI: shgi_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2
    FPP_TIMER_START(tot)
    FPP_TIMER_START(totI)

    call setif(-1_i8k,.false.) ! zero all bits
    if( present(IMOD) )then
      call setif(IMOD)
    else
      call setif(INUCL,present(NUCL))
      call setif(INUSR,present(NUSR))
      call setif(INUSO,present(NUSO))
      call setif(ICHSR,present(CHSR))
      call setif(ICHSO,present(CHSO))
      call setif(IPSEU,present(PSEU))
      call setif(IOVRL,present(OVRL))
      call setif(IKNTC,present(KNTC))
    endif

    if(is_on(INUCL).and..not.present(NUCL))then
      ABORT('NUCL?')
    endif
    if(is_on(INUSR).and..not.present(NUSR))then
      ABORT('NUSR?')
    endif
    if(is_on(IOVRL).and..not.present(OVRL))then
      ABORT('OVRL?')
    endif
    if(is_on(IKNTC).and..not.(present(KNTC).or.present(NUCL)))then
      ABORT('KNGR?')
    endif
    if(is_on(IPSEU).and..not.(present(PSEU).or.present(NUCL)))then
      ABORT('PSGR?')
    endif
    if(is_on(ICHSR).and..not.present(CHSR))then
      ABORT('CHSR?')
    endif
    if(is_on(ICHSO).and..not.present(CHSO))then
      ABORT('CHSO?')
    endif

    if(is_on(IADKH))then
      ASSERT(is_on(IOVRL))
      ASSERT(is_on(IKNTC))
      ASSERT(is_on(INUCL))
      ! both do not make sense:
      ASSERT(.not.is_on(INUSR))
      ! nevertheless set as ADKH may need for SR-integrals:
      call setif(INUSR)
      ASSERT(is_on(INUSR))
    endif

!    DPRINT 'SHGI: shgi_drv: MODE=',MODE

    if( is_on(IXXSR+IXXSO) )then
       ! if ANY (bitmask!) of them is set then raise LC
       LC = 1
       ! FIXME: when doing p-,d-,.. fit functions
    else
       LC = 0
    endif

    ! RelFit [and possibly PPs] require precomputed angular factor:
    if( is_on(ICHSR) ) call setif(IYSNU)
    if( is_on(ICHSR) ) call setif(IYSSR)
    if( is_on(ICHSO) ) call setif(IYSSO)
#if SHGI_USE_ANGULAR
    if( is_on(IPSEU) ) call setif(IYSNU)
    ! for PPs it may be faster not to precompute
    ! the angular factor of LOCAL part!
#endif

    ! set bits to flag coinciding centers: A==B?:
    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    ! indices of (ua,L) shells:
    uaL = shellNum(IU2,IL2)
    ubL = shellNum(IU1,IL1)

    ! set sameXX to the index of coinciding centers, -1 otherwise:
    sameUA = -1
    sameEA = -1
    sameL  = -1
    if( IU2 == IU1 ) sameUA = IU1 ! == IU2
    if( IE2 == IE1 ) sameEA = IE1 ! == IE2
    if( IL2 == IL1 ) sameL  = IL1 ! == IL2
    diagonal = (sameUA>0) .and. (sameEA>0) .and. (sameL>=0)

    ! For some reason even conservative screening of
    ! integrals (MAXEXP==50) causes the MO overlap matrix
    ! to be non-positive definite -- dpptrf() returns INFO>0.
    ! So in normal integral run dont use any screening:
    call shgi_set_maxexp(HUGE(1.0_rk))
    ! OTOH, in "property" runs, eg. gradients, second derivatives
    ! it may be safer to apply integral screening.
    ! See shgi_gr_drv(), shgi_sd_drv(), etc ...

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! allocate global storage of THIS module
    call shgi_glob_alloc(NAB,LA,LB)

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(LC,0,0)

    ! 2c overlap and kinetic integrals:
    call shgi_set_ovrl(LA,LB,LC,0,0)  ! (1) shgi_drv

    if(is_on(IADKH))then
      ! prepare and store relativistic projection coefficients:
      call shgi_adkh_atom(IL2, uas(IU2)%l_ob(IL2)%exponents &
                             , uas(IU2)%Z                   &
                             , uas(IU2)%nuclear_radius      &
                             , uas(IU2)%s_core%exponents    &
                             , uas(IU2)%s_core%contractions &
                             , uaL                          &
                             )
      call shgi_adkh_atom(IL1, uas(IU1)%l_ob(IL1)%exponents &
                             , uas(IU1)%Z                   &
                             , uas(IU1)%nuclear_radius      &
                             , uas(IU1)%s_core%exponents    &
                             , uas(IU1)%s_core%contractions &
                             , ubL                          &
                             )
      ! note that a==2, b==1 ...
    endif

    if(is_on(IKNTC))then
      ASSERT(is_on(IOVRL))
      ! kinetic energy goes always with overlap:
      if(is_on(IADKH))then
        ! compute kinetic energy with relativistic corrections:
!       write(*,"(' shgi_adkh: INTS <',3I3,'|     KIN     |',3I3,'>')") IU2,IE2,IL2,IU1,IE1,IL1
        call shgi_adkh_kin(PIOVRL,PIKNTC,uaL,ubL)
      else
        call shgi_kin(PIOVRL,PIKNTC)
      endif
    endif

    ! unpack overlap (maybe already with rel corrections):
    if(is_on(IOVRL))then
       ASSERT(present(OVRL))
       call upack2c(PIOVRL,OVRL)
    endif

    ! unpack (maybe relativistic) kinetic energy:
    if(is_on(IKNTC))then
       ASSERT(present(KNTC))
       call upack2c(PIKNTC,KNTC)
    endif

    if(is_on(INUCL))then
       PINUCL = ZERO
    endif

    if(is_on(INUSR))then
       PINUSR = ZERO
    endif

    if(is_on(INUSO))then
       PINUSO = ZERO
    endif

    if(is_on(IPSEU))then
       PIPSEU = ZERO
       call shgi_pseu_set_ab() ! just calls bessel_setup()
    endif

    ! reset counters for the s- and r2-fitfunctions
    ich = 0
    ir2 = 0

    ! == LOOP: UNIQUE ATOMS (C)
    do ua=1,size(uas)

       if(is_on(INUXX))then

          Z = uas(ua)%Z
          ! ZCORE is handled by shgi_pseu0()
          R = radius(ua)
       endif

       if(is_on(IPSEU))then
          LPNUM  = uas(ua)%lmax_pseudo
          ZCORE  = uas(ua)%Z - uas(ua)%ZC
          if( LPNUM > 0 )then
             Z = zero ! disable NUCL
             ! prepare staff common for all EAs:
             call shgi_pseu_set_abc(LPNUM-1,0) ! case=0 i.e. Energy
          endif
       else
          LPNUM = -1 ! PP is not computed, even if there are some?
       endif

       if(is_on(ICHXX))then
          call shgi_set_ec(STYPE , uas(ua)%l_ch(0)%exponents)
          call shgi_set_ec(R2TYPE, uas(ua)%r2_ch%exponents  )

          if(is_on(ICHSR))then
             PICHSR = ZERO
             PIR2SR = ZERO
          endif

          if(is_on(ICHSO))then
             PICHSO = ZERO
             PIR2SO = ZERO
          endif
       endif

       ! == LOOP: EQUIVALENT ATOMS (C)
       do eq=1,uas(ua)%N_equal_atoms

          ! set bits to flag coinciding centers: A==B? A==C? B==C?:
          call shgi_set_xeqy(IU2,IE2,IU1,IE1,ua,eq)
          xc = uas(ua)%position(:,eq)

          if( is_on(ICHXX) )then
            ! in case of the fit integrals it makes sense to
            ! precompute angular part which is done by:
            call shgi_set_c(xc)
          endif

          ![[==========================================
          ! NUCLEAR ATTRACTION AND ITS SR/SO MAT. ELs.:
          if( is_on(IYSNU) )then ! angular part available:

            if(is_on(INUCL))then
               ! NR,SR, and SO runs:
               call shgi_nuc(   Z,PINUCL,RAD=R)
            endif

            if(is_on(INUSR,IYSSR))then
               ! SR, and SO runs:
               call shgi_nuc_sr(Z,PINUSR,RAD=R)
            else
               ABORT('IYSSR?')
            endif

            if(is_on(INUSO,IYSSO))then
               ! SO runs:
               call shgi_nuc_so(Z,PINUSO,RAD=R)
            else
               ABORT('IYSSO?')
            endif

          else ! angular part not available:

            if(     is_on(INUCL,INUSR,INUSO) )then
              ! SO runs:
              call shgi_rel_nuc(Z,xc,PINUCL,PINUSR,PINUSO,RAD=R)

            elseif( is_on(IADKH) )then

              ! increments by the non-rel NUCL(of C) if not (A==B==C):
!             write(*,"(' shgi_adkh: INTS <',3I3,'|',2I3,F7.1,'|',3I3,'>')") IU2,IE2,IL2,ua,eq,Z,IU1,IE1,IL1
              call shgi_adkh_nuc(Z,xc,PINUCL,R,uaL,ubL)

            elseif( is_on(INUCL,INUSR) )then
              ! SR runs:
              call shgi_rel_nuc(Z,xc,PINUCL,PINUSR,RAD=R)

            elseif( is_on(INUCL) )then
              ! NR runs:
              call shgi_rel_nuc(Z,xc,PINUCL,RAD=R)
            endif
          endif ! is_on(IYSNU)
          !]]==========================================

          ![[==========================================
          ! SR/SO MAT. ELs. OF REL-FIT:
          if(is_on(ICHSR))then
             call shgi_ch_sr(PICHSR)
             call shgi_r2_sr(PIR2SR)
          endif

          if(is_on(ICHSO))then
             call shgi_ch_so(PICHSO)
             call shgi_r2_so(PIR2SO)
          endif
          !]]==========================================

          if( LPNUM > 0 )then
             ! do pseudopotentials if needed:
             call shgi_pseu(ZCORE,LPNUM-1,uas(ua)%l_pseudopot,xc,PIPSEU)
          endif

       enddo ! loop over equivalent atoms eq

       ! unpack the s-ints starting at address 'ich'
       cnt = ich
       if(is_on(ICHSR))then
          cnt = ich
          call upack3c(cnt,PICHSR,CHSR)
       endif

       if(is_on(ICHSO))then
          cnt = ich
          call upack3c(cnt,PICHSO(:,:,:,:,1),CHSO(:,:,:,:,:,1))
          cnt = ich
          call upack3c(cnt,PICHSO(:,:,:,:,2),CHSO(:,:,:,:,:,2))
          cnt = ich
          call upack3c(cnt,PICHSO(:,:,:,:,3),CHSO(:,:,:,:,:,3))
       endif
       ich = cnt

       ! unpack the r2-ints starting at address 'ir2'
       cnt = ir2
       if(is_on(ICHSR))then
          cnt = ir2
          call upack3c(cnt,PIR2SR,R2SR)
       endif

       if(is_on(ICHSO))then
          cnt = ir2
          call upack3c(cnt,PIR2SO(:,:,:,:,1),R2SO(:,:,:,:,:,1))
          cnt = ir2
          call upack3c(cnt,PIR2SO(:,:,:,:,2),R2SO(:,:,:,:,:,2))
          cnt = ir2
          call upack3c(cnt,PIR2SO(:,:,:,:,3),R2SO(:,:,:,:,:,3))
       endif
       ir2 = cnt

       call shgi_close_abc()
       if( LPNUM > 0 )then
         call shgi_pseu_close_abc()
        endif
    enddo ! loop over unique atoms ua

    if(is_on(ICHSR))then
       ASSERT(ich==size(CHSR,3))
       ASSERT(ir2==size(R2SR,3))
    endif

    if(is_on(ICHSO))then
       ASSERT(ich==size(CHSO,3))
       ASSERT(ir2==size(R2SO,3))
    endif

    ![[=== add field of PCs to nuclear attraction: ===
    npc = 0
    if( present(pcs) ) npc = size(pcs)
    if( npc > 0 )then

       ! WARNING: since no SR-ints for PCs re-setting LC:
       call shgi_set_lcde(0,0,0)
       ! FIXME: do I need to reshape S5?

       if( is_on(IPSEU) )then
         call shgi_pcs(pcs,PIPSEU)
       else
         ASSERT(is_on(INUCL))
         call shgi_pcs(pcs,PINUCL)
       endif
    endif
    ![[===============================================


    ![[=== add field of external centers to nuclear attraction: ===
    if( is_on(IPSEU) )then
       call shgi_ext(PIPSEU)
    else
       ASSERT(is_on(INUCL))
       call shgi_ext(PINUCL)
    endif
    ![[============================================================

    ! unpack (maybe relativistic) kinetic energy:
    if(is_on(IKNTC))then
       ASSERT(present(KNTC))
       call upack2c(PIKNTC,KNTC)
    endif

    ! unpack nuclear attraction:
    if(is_on(INUCL))then
       ! at least with EPE the NUCL storage already
       ! contains something, therefore add, not overwrite:
       call upack2c(PINUCL,NUCL,+1)
!!!    call upack2c(PINUCL,NUCL)
    endif

    ! unpack pseudopotential:
    if(is_on(IPSEU))then
       if( present(PSEU) )then
          ! in a realtivistic computation
          ! pseudopotential needs to be separated
          ! from NUCL:
          call upack2c(PIPSEU,PSEU)
       else
          ! just add it to the NUCL:
          ASSERT(present(NUCL))
          call upack2c(PIPSEU,NUCL,+1)
       endif
    endif

    ! unpack SR of nuclear attraction:
    if(is_on(INUSR) .and..not.is_on(IADKH))then
       call upack2c(PINUSR,NUSR)
    endif

    ! unpack SO of nuclear attraction:
    if(is_on(INUSO))then
       call upack2c(PINUSO(:,:,:,1),NUSO(:,:,:,:,1))
       call upack2c(PINUSO(:,:,:,2),NUSO(:,:,:,:,2))
       call upack2c(PINUSO(:,:,:,3),NUSO(:,:,:,:,3))
    endif

    call shgi_close_ab()
    ! deallocate global storage of THIS module:
    call shgi_glob_free()
    FPP_TIMER_STOP(totI)
    FPP_TIMER_STOP(tot)

    DCALL shgi_timing()
    DPRINT  'SHGI: shgi_drv: exit'

  contains
    function radius(ua) result(r)
      !
      ! returns uas(ua)%nuclear_radius
      ! or zero if the latter is negative
      !
      implicit none
      integer(IK), intent(in) :: ua
      real(RK)                :: r
      ! *** end of interface ***

      ASSERT(ua>0)
      ASSERT(ua<=size(uas))

      ! uas(:) is host-associated from enclosing sub:
      r = uas(ua)%nuclear_radius
      r = max(r,0.0_rk) ! avoid negative
    end function radius
  end subroutine shgi_drv

  !****************************************************************
  !****************** DRIVER FOR GRADIENTS ************************
  !****************************************************************

  subroutine shgi_gr_drv(IU1,IE1,IL1,IU2,IE2,IL2, uas, pcs, &
       NUGR, &
       SRGR, &
       PSGR, &
       OVGR, &
       KNGR, &
       PSSD, &
       IMOD  &
       )
    use unique_atom_module, only: uat=>unique_atom_type
    use datatype          , only: pct=>pointcharge_type
    use datatype, only: arrmat4 ! PG storage of grs
    use options_module, only: options_integral_expmax
!   use shgi_ang, only: shgi_set_c
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_ab,  only: shgi_gr_kin
    use shgi_utils , only: shgi_gr_wd_store, shgi_timing
    use shgi_utils , only: shgi_gr_wd_to_abc, shgi_gr_abc_store
    use shgi_utils , only: shellNum!(ua,L)
    use shgi_relfit, only: shgi_gr_nuc, shgi_gr_nuc_sr
    use shgi_relnuc, only: shgi_rel_gr_nuc
    use shgi_adkh  , only: shgi_adkh_gr_kin, shgi_adkh_gr_nuc, shgi_adkh_atom
    use shgi_pseudo, only: &
         shgi_pseu_set_ab, shgi_pseu_set_abc, shgi_pseu_close_abc, shgi_gr_pseu
    use shgi_pcm, only: shgi_gr_pcs
    use shgi_pcm, only: shgi_sd_pcs
    use shgi_ext_c, only: shgi_ext_gr
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)       :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)       :: uas(:)  ! array of unique atoms, normally all of them
    type(pct)  , intent(in)       :: pcs(:)  ! array of point charges
    type(arrmat4) , intent(inout) :: NUGR(:) ! (num of tot-sym grads)
    type(arrmat4) , intent(inout) :: SRGR(:) ! (num of tot-sym grads)
    type(arrmat4) , intent(inout) :: PSGR(:) ! (num of tot-sym grads)
    type(arrmat4) , intent(inout) :: OVGR(:) ! (num of tot-sym grads)
    type(arrmat4) , intent(inout) :: KNGR(:) ! (num of tot-sym grads)
    type(arrmat4) , intent(inout) :: PSSD(:,:) ! (num of tot-sym grads)
    ! each %m(:,:,:,:) is of the shape         (N2E,N1E,2*L2+1,2*L1+1)
    integer(I8K), intent(in)       :: IMOD
    optional :: pcs
    optional :: NUGR
    optional :: SRGR
    optional :: PSGR
    optional :: OVGR
    optional :: KNGR
    optional :: PSSD
    optional :: IMOD
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------

    integer(IK) :: ua,eq
    real(RK)    :: Z,R
    real(RK)    :: xc(3)

    integer(IK) :: LPNUM
    real(RK)    :: ZCORE

    integer(IK) :: npc  ! size(pcs)

    logical     :: diagonal  ! flags if UA==UB, EA==EB, and LA=LB
    integer(IK) :: sameUA    ! ==IU2==IU1 if diagonal, -1 otherwise
    integer(IK) :: sameEA    ! ==IE2==IE1 if diagonal, -1 otherwise
    integer(IK) :: sameL     ! ==IL2==IL1 if diagonal, -1 otherwise
    integer(IK) :: uaL,ubL   ! linear index of (UA,LA) and (OB,LB) shells

!   integer(IK):: i,j

    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    DPRINT  'SHGI: shgi_gr_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2
    FPP_TIMER_START(tot)
    FPP_TIMER_START(totG)

    call setif(-1_i8k,.false.) ! zero all bits
    if( present(IMOD) )then
      call setif(IMOD)
    else
      call setif(INUGR,present(NUGR))
      call setif(ISRGR,present(SRGR))
      call setif(IPSGR,present(PSGR))
      call setif(IOVGR,present(OVGR))
      call setif(IKNGR,present(KNGR))
      call setif(IPSSD,present(PSSD))
    endif

    if(is_on(INUGR).and..not.present(NUGR))then
      ABORT('NUGR?')
    endif
    if(is_on(ISRGR).and..not.present(SRGR))then
      ABORT('SRGR?')
    endif
    if(is_on(IOVGR).and..not.present(OVGR))then
      ABORT('OVGR?')
    endif
    if(is_on(IKNGR).and..not.(present(KNGR).or.present(NUGR)))then
      ABORT('KNGR?')
    endif
    if(is_on(IPSGR).and..not.(present(PSGR).or.present(NUGR)))then
      ABORT('PSGR?')
    endif
    if(is_on(IPSSD).and..not.present(PSSD))then
      ABORT('PSSD?')
    endif

    if(is_on(IPSSD))then
      WARN('shgi_gr_drv: ignore PSSD!')
      call setif(IPSSD,.false.)
    endif
    if(is_on(IPCSD))then
      WARN('shgi_gr_drv: ignore PCSD!')
      call setif(IPCSD,.false.)
    endif

    if(is_on(IADKH))then
      ASSERT(is_on(IOVGR))
      ASSERT(is_on(IKNGR))
      ASSERT(is_on(INUGR))
      ! both do not make sense:
      ASSERT(.not.is_on(ISRGR))
      ! nevertheless set as ADKH may need for SR-integrals:
      call setif(ISRGR)
    endif

!    DPRINT 'SHGI: shgi_gr_drv: MODE=',MODE

    ! FIXME: raise more when doing p-,d-,.. fit functions
    LC = 1
    LD = 0

    if( is_on(ISRGR+IPSSD+IPCSD) )then
       ! for gradients of SR mat. el. one more derivative:
       ! or if PC second derivatives or pseudo second derivative
       LD = 1
    endif

#if SHGI_USE_ANGULAR
    ! PPs possibly require precomputed angular factor:
    if( is_on(IPSGR) ) call setif(IGSNU)
    ! for PPs it may be faster not to precompute
    ! the angular factor of LOCAL part!
#endif

    ! set bits to flag coinciding centers: A==B?:
    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    ! indices of (ua,L) shells:
    uaL = shellNum(IU2,IL2)
    ubL = shellNum(IU1,IL1)

    ! set sameXX to the index of coinciding centers, -1 otherwise:
    sameUA = -1
    sameEA = -1
    sameL  = -1
    if( IU2 == IU1 ) sameUA = IU1 ! == IU2
    if( IE2 == IE1 ) sameEA = IE1 ! == IE2
    if( IL2 == IL1 ) sameL  = IL1 ! == IL2
    diagonal = (sameUA>0) .and. (sameEA>0) .and. (sameL>=0)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! allocate global storage of THIS module
    call shgi_glob_alloc(NAB,LA,LB)

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(LC,LD,0)

    ! 2c overlap and kinetic gradients:
    call shgi_set_ovrl(LA,LB,LC,LD,0)  ! (2) in shgi_gr_drv

    if(is_on(IADKH))then
      ! prepare and store relativistic projection coefficients:
      ! ( in case any of A or B wasnt touched on this processor before scf ... )
      call shgi_adkh_atom(IL2, uas(IU2)%l_ob(IL2)%exponents &
                             , uas(IU2)%Z                   &
                             , uas(IU2)%nuclear_radius      &
                             , uas(IU2)%s_core%exponents    &
                             , uas(IU2)%s_core%contractions &
                             , uaL                          &
                             )
      call shgi_adkh_atom(IL1, uas(IU1)%l_ob(IL1)%exponents &
                             , uas(IU1)%Z                   &
                             , uas(IU1)%nuclear_radius      &
                             , uas(IU1)%s_core%exponents    &
                             , uas(IU1)%s_core%contractions &
                             , ubL                          &
                             )
      ! note that a==2, b==1 ...
    endif

    if(is_on(IKNGR))then
      ASSERT(is_on(IOVGR))
      ! kinetic energy goes always with overlap:
      call shgi_gr_kin(GROVRL,GRKNTC)
      ! NOTE: if A==B then gradients sum to zero G(A)+G(B)==0
      !       and need not to be evaluated
      !       (they are set to zero inside of the sub)

      if(is_on(IADKH).and..not.is_on(IAEQB))then
!       write(*,"(' shgi_adkh: GRAD <',3I3,'|     KIN     |',3I3,'>')") IU2,IE2,IL2,IU1,IE1,IL1
        ! apply relativistic corrections kinetic energy (and overlap):
        call shgi_adkh_gr_kin(GROVRL,GRKNTC,uaL,ubL)
        ! NOTE: if A==B then gradients sum to zero G(A)+G(B)==0
        !       and no trafo is necessary
      endif
    endif

    if(is_on(IOVGR))then
       ASSERT(present(OVGR))
       ! FIXME: note the order of UAs:
       call shgi_gr_abc_store(IU1,IE1,IU2,IE2,0,2,GROVRL,OVGR, 0) ! 0==overwrite
       ! (0,2) is a special case of 2c-grad into ITS OWN storage.
    endif

    if(is_on(IKNGR))then
       if(present(KNGR))then
         ! FIXME: note the order of UAs:
         call shgi_gr_abc_store(IU1,IE1,IU2,IE2,0,2,GRKNTC,KNGR, 0) ! 0==overwrite
         ! (0,2) is a special case of 2c-grad into ITS OWN storage.
       else
         ! - ( KIN + NUC ) grads both go to NUGR:
         ASSERT(present(NUGR))
         call shgi_gr_abc_store(IU2,IE2,IU1,IE1,0,1,GRKNTC,NUGR,-1)
         ! (0,1) is a special case of 2c-grad into COMMON storage.
       endif
    endif

    if(is_on(IPSGR))then
       call shgi_pseu_set_ab() ! just calls bessel_setup()
    endif

    ! == LOOP: UNIQUE ATOMS (C)
    do ua=1,size(uas)
       if(is_on(INUGR+ISRGR))then ! if any of them

          Z = uas(ua)%Z
          ! ZCORE is handled by shgi_pseu0()
          Z = Z - uas(ua)%ZC
          R = radius(ua)
       endif

       if(is_on(IPSGR))then
          LPNUM  = uas(ua)%lmax_pseudo
          ZCORE  = uas(ua)%Z - uas(ua)%ZC
          if( LPNUM > 0 )then
             Z = zero ! disable NUCL
             ASSERT(.not.is_on(IPSSD))
             !(use shgi_sd_drv(...) for that)
             ! prepare staff common for all EAs:
             call shgi_pseu_set_abc(LPNUM-1,1) ! case=1, i.e. Gradients
          endif
       else
          LPNUM = -1 ! PP is not computed, even if there are some?
       endif

       ! == LOOP: EQUIVALENT ATOMS (C)
       do eq = 1, uas(ua) % N_equal_atoms
          ! set bits to flag coinciding centers: A==B? A==C? B==C?:
          call shgi_set_xeqy(IU2,IE2,IU1,IE1,ua,eq)

          if (whatis (IXEQC) == IXEQC) cycle ! eq atom loop, as gradients are zero

          xc = uas(ua)%position(:,eq)

          ! FXIME: only when precomputing angular part
          !        brings something:
!         call shgi_set_c(xc)

          if(is_on(INUGR,ISRGR))then
             ASSERT(present(NUGR))
             ! both required, do them together:
             GRNUCL = ZERO
             GRNUSR = ZERO
             call shgi_rel_gr_nuc(Z,xc,GRNUCL,GRNUSR,RAD=R)

             if(is_on(IADKH))then
               ! transform and store the grads of V and SR ints...
!              write(*,"(' shgi_adkh: GRAD <',3I3,'|',2I3,F7.1,'|',3I3,'>')") IU2,IE2,IL2,ua,eq,Z,IU1,IE1,IL1
               ASSERT(allocated(GRABC1))
               ASSERT(allocated(GRABC2))

               ! transform grads wrt (W,D) to grads wrt (A,B,C):
               call shgi_gr_wd_to_abc(GRNUCL,GRABC1,0)
               call shgi_gr_wd_to_abc(GRNUSR,GRABC2,0)
               ! NOTE: need to do it to apply rel contractions over
               !       alpha- and beta-exponents, but vector
               !                W = wa*A + wb*B - C
               !       depends on them. Vectors A, B, and C -- dont.

               ! applies rel contraction to grads of V (in GRABC1)
               ! using grads of SR ints (in GRABC2), outputs grads of Vrel in-place:
               call shgi_adkh_gr_nuc(GRABC1,GRABC2,uaL,ubL)

               ! rotates and stores ABC-gradients in PG structures:
               call shgi_gr_abc_store(IU2,IE2,IU1,IE1,ua,eq,GRABC1,NUGR,-1)
               ! ... but ignore no more necessary SR ints (in GRABC2)
             else
               ! true DKH: just store the grads of V and SR ints...
               ASSERT(present(SRGR))
               ! rotates and stores ABC-gradients in PG structures:
               call shgi_gr_wd_store(IU2,IE2,IU1,IE1,ua,eq,GRNUCL,NUGR,-1)
               call shgi_gr_wd_store(IU2,IE2,IU1,IE1,ua,eq,GRNUSR,SRGR,-1)
             endif
          else
             ! if any required, do it:
             if(is_on(INUGR))then
                GRNUCL = ZERO

                if(is_on(IGSNU))then
                   call shgi_gr_nuc(Z,GRNUCL,RAD=R)
                else
                   call shgi_rel_gr_nuc(Z,xc,GRNUCL,RAD=R)
                endif

                ! rotates and stores ABC-gradients in PG structures:
                call shgi_gr_wd_store(IU2,IE2,IU1,IE1,ua,eq,GRNUCL,NUGR,-1)
             endif

             if(is_on(ISRGR))then
                ASSERT(present(SRGR))
                GRNUSR = ZERO

                if(is_on(IGSSR,IGSNU))then
                   call shgi_gr_nuc_sr(Z,GRNUSR,RAD=R)
                else
                   call shgi_rel_gr_nuc(Z,xc,SR=GRNUSR,RAD=R)
                endif

                ! rotates and stores ABC-gradients in PG structures:
                call shgi_gr_wd_store(IU2,IE2,IU1,IE1,ua,eq,GRNUSR,SRGR,-1)
             endif
          endif ! is_on(INUCL,INUSR)

          if( LPNUM > 0 )then
             ! do pseudopotential gradients if needed:
             ! w.r.t A, B, and C:


             GRPSEU = ZERO
             call shgi_gr_pseu(ZCORE,LPNUM-1,uas(ua)%l_pseudopot,xc,GRPSEU)!,FIXME: GRPSEU)

             if(present(PSGR))then
                ! in a realtivistic computation
                ! pseudopotential needs to be separated from NUCL:
                call shgi_gr_abc_store(IU2,IE2,IU1,IE1,ua,eq,GRPSEU,PSGR,-1)
             else
                ! just add it to the NUCL:
                ASSERT(present(NUGR))
                call shgi_gr_abc_store(IU2,IE2,IU1,IE1,ua,eq,GRPSEU,NUGR,-1)
             endif

             ! FIXME: first accumulate A and B and do this once!
             if(is_on(IPSSD)) then
              ABORT('use shgi_sd_drv')
!             SDPSEU=ZERO
!             call shgi_sd_pseu(ZCORE,LPNUM-1,uas(ua)%l_pseudopot,xc,SDPSEU)
!             call shgi_sd_abc_store(IU2,IE2,IU1,IE1,ua,eq,SDPSEU,PSSD,-1)
             endif

          endif

       enddo ! loop over equivalent atoms eq

       call shgi_close_abc()
       if( LPNUM > 0 )then
          call shgi_pseu_close_abc()
       endif
    enddo ! loop over unique atoms ua

    ![[=== add field of PCs to nuclear attraction: ===
    npc = 0
    if( present(pcs) ) npc = size(pcs)
    if( npc > 0 )then

       ! set bits to flag coinciding centers: A==B? A==C? B==C?:
       ! cannot use U2==U1 this way:
       call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

       ! WARNING: since no SR-ints for PCs re-setting LD:
       LD=0
       if(is_on(IPCSD)) LD=1
       call shgi_set_lcde(LC,LD,0)
       ! FIXME: do I need to reshape S5?

       if( present(PSGR) )then
         ASSERT(is_on(IPSGR))
         ASSERT(allocated(GRPSEU))
         GRPSEU = zero
         call shgi_gr_pcs(pcs,GRPSEU)
         ! rotates and stores ABC-gradients in PG structures:
         DPRINT  'PSGR> shgi_gr_wd_store(',IU2,IE2,IU1,IE1,0,0,')'
         call shgi_gr_wd_store(IU2,IE2,IU1,IE1,0,0,GRPSEU,PSGR,-1)
       else
         ASSERT(is_on(INUGR))
         ASSERT(allocated(GRNUCL))
         GRNUCL = zero
         call shgi_gr_pcs(pcs,GRNUCL)
         ! rotates and stores ABC-gradients in PG structures:
         ASSERT(present(NUGR))
         DPRINT  'NUGR> shgi_gr_wd_store(',IU2,IE2,IU1,IE1,0,0,')'
         call shgi_gr_wd_store(IU2,IE2,IU1,IE1,0,0,GRNUCL,NUGR,-1)
       endif
    endif
    ![[===============================================

    ![[=== add field of external centers to nuclear attraction: ===
    if(present(PSGR))then
       call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)
       GRPSEU = zero
       call shgi_ext_gr(GRPSEU)
       call shgi_gr_wd_store(IU2,IE2,IU1,IE1,0,0,GRPSEU,PSGR,-1)
    else
       call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)
       GRNUCL = zero
       call shgi_ext_gr(GRNUCL)
       call shgi_gr_wd_store(IU2,IE2,IU1,IE1,0,0,GRNUCL,NUGR,-1)
    endif
    ![[============================================================


    ! FIXME: first accumulate [AB] and store them once here:
    ! a) rotate the [AB]-gradients along the directions
    !    of (totally) symmetric modes,
    ! b) then add them to PG datastructures:

    call shgi_close_ab()
    ! deallocate global storage of THIS module:
    call shgi_glob_free()
    FPP_TIMER_STOP(totG)
    FPP_TIMER_STOP(tot)
    DCALL shgi_timing()
    DPRINT  'SHGI: shgi_gr_drv: exit'

  contains
    function radius(ua) result(r)
      !
      ! returns uas(ua)%nuclear_radius
      ! or zero if the latter is negative
      !
      implicit none
      integer(IK), intent(in) :: ua
      real(RK)                :: r
      ! *** end of interface ***

      ASSERT(ua>0)
      ASSERT(ua<=size(uas))

      ! uas(:) is host-associated from enclosing sub:
      r = uas(ua)%nuclear_radius
      r = max(r,0.0_rk) ! avoid negative
    end function radius
  end subroutine shgi_gr_drv

  !****************************************************************
  !************* DRIVER FOR SECOND DERIVATIVES ********************
  !****************************************************************

  subroutine shgi_sd_drv(IU1,IE1,IL1,IU2,IE2,IL2, uas, pcs, &
       NUSD, &
       SRSD, &
       OVSD, &
       KNSD, &
       PSSD, &
       IMOD  &
       )
    use unique_atom_module, only: uat=>unique_atom_type
    use datatype          , only: pct=>pointcharge_type
    use datatype, only: arrmat4 ! PG storage of grs
    use options_module, only: options_integral_expmax
!   use shgi_ang, only: shgi_set_c
    use shgi_ang, only: shgi_set_lcde
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab, shgi_set_ovrl
    use shgi_ab,  only: shgi_sd_kin
    use shgi_adkh,   only: shgi_adkh_sd_kin, shgi_adkh_sd_nuc, shgi_adkh_atom
    use shgi_utils , only: shgi_sd_wd_store, shgi_timing
    use shgi_utils , only: shgi_sd_wd_to_abc, shgi_sd_abc_store
    use shgi_utils , only: shellNum!(ua,L)
    use shgi_relnuc, only: shgi_rel_sd_nuc
    use shgi_pseudo, only: &
         shgi_pseu_set_ab, shgi_pseu_set_abc, shgi_pseu_close_abc, shgi_sd_pseu
    use shgi_pcm,    only: shgi_sd_pcs
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)       :: IU1,IE1,IL1,IU2,IE2,IL2
    type(uat)  , intent(in)       :: uas(:)    ! array of unique atoms, normally all of them
    type(pct)  , intent(in)       :: pcs(:)    ! array of point charges
    type(arrmat4) , intent(inout) :: NUSD(:,:) ! (num of tot-sym grads^2)
    type(arrmat4) , intent(inout) :: SRSD(:,:) ! (num of tot-sym grads^2)
    type(arrmat4) , intent(inout) :: OVSD(:,:) ! (num of tot-sym grads^2)
    type(arrmat4) , intent(inout) :: KNSD(:,:) ! (num of tot-sym grads^2)
    type(arrmat4) , intent(inout) :: PSSD(:,:) ! (num of tot-sym grads^2)
    ! each %m(:,:,:,:) is of the shape         (N2E,N1E,2*L2+1,2*L1+1)
    integer(I8K), intent(in)       :: IMOD
    optional :: pcs
    optional :: NUSD
    optional :: SRSD
    optional :: OVSD
    optional :: KNSD
    optional :: PSSD
    optional :: IMOD
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------

    integer(IK) :: ua,eq
    real(RK)    :: Z,R
    real(RK)    :: xc(3)

    integer(IK) :: LPNUM
    real(RK)    :: ZCORE

    integer(IK) :: npc  ! size(pcs)

    logical     :: diagonal  ! flags if UA==UB, EA==EB, and LA=LB
    integer(IK) :: sameUA    ! ==IU2==IU1 if diagonal, -1 otherwise
    integer(IK) :: sameEA    ! ==IE2==IE1 if diagonal, -1 otherwise
    integer(IK) :: sameL     ! ==IL2==IL1 if diagonal, -1 otherwise
    integer(IK) :: uaL,ubL   ! linear index of (UA,LA) and (OB,LB) shells
    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    DPRINT  'SHGI: shgi_sd_drv: UAs=',IU1,IU2,' EAs=',IE1,IE2,' Ls=',IL1,IL2
    FPP_TIMER_START(tot)
    FPP_TIMER_START(totD)

    call setif(-1_i8k,.false.) ! zero all bits
    if( present(IMOD) )then
      call setif(IMOD)
    else
      call setif(INUSD,present(NUSD))
      call setif(ISRSD,present(SRSD))
      call setif(IOVSD,present(OVSD))
      call setif(IKNSD,present(KNSD))
      call setif(IPSSD,present(PSSD))
    endif

    if(is_on(INUSD).and..not.present(NUSD))then
      ABORT('NUSD?')
    endif
    if(is_on(ISRSD).and..not.present(SRSD))then
      ABORT('SRSD?')
    endif
    if(is_on(IOVSD).and..not.present(OVSD))then
      ABORT('OVSD?')
    endif
    if(is_on(IKNSD).and..not.(present(KNSD).or.present(NUSD)))then
      ABORT('KNSD?')
    endif
    if(is_on(IPSSD).and..not.(present(PSSD).or.present(NUSD)))then
      ABORT('PSSD?')
    endif

    if(is_on(IADKH))then
      ASSERT(is_on(IOVSD))
      ASSERT(is_on(IKNSD))
      ASSERT(is_on(INUSD))
      ! both do not make sense:
      ASSERT(.not.is_on(ISRSD))
      ! nevertheless set as ADKH may need for SR-integrals:
      call setif(ISRSD)
    endif

    DPRINT 'SHGI: shgi_sd_drv: MODE=',is_on(INUSD),is_on(ISRSD),is_on(IPSSD),whatis(-1_8)
    DPRINT 'SHGI: shgi_sd_drv: present(pcs)=',present(pcs)

    ! FIXME: raise more when doing p-,d-,.. fit functions
    LC = 1
    LD = 1

    if( is_on(ISRSD) )then
       ! for gradients of SR mat. el. one more derivative:
       LE = 1
    endif

    ! PPs require precomputed angular factor:
!   if( is_on(IPSGR) ) call setif(IGSNU)
    ! FIXME: for PPs it may be faster not to precompute
    !        the angular factor of LOCAL part!
    !        Fix shgi_gr_pseu0() .

    ! set bits to flag coinciding centers: A==B?:
    call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

    ! indices of (ua,L) shells:
    uaL = shellNum(IU2,IL2)
    ubL = shellNum(IU1,IL1)

    ! set sameXX to the index of coinciding centers, -1 otherwise:
    sameUA = -1
    sameEA = -1
    sameL  = -1
    if( IU2 == IU1 ) sameUA = IU1 ! == IU2
    if( IE2 == IE1 ) sameEA = IE1 ! == IE2
    if( IL2 == IL1 ) sameL  = IL1 ! == IL2
    diagonal = (sameUA>0) .and. (sameEA>0) .and. (sameL>=0)

    ! Use default integral screening:
    call shgi_set_maxexp(options_integral_expmax())

    call shgi_set_ab( IL2, IL1, &
         uas(IU2)%position(:,IE2)    , &
         uas(IU1)%position(:,IE1)    , &
         uas(IU2)%l_ob(IL2)%exponents, &
         uas(IU1)%l_ob(IL1)%exponents  &
         )
    ! <- also sets LA, LB, and NAB

    ! allocate global storage of THIS module
    call shgi_glob_alloc(NAB,LA,LB)

    ! set global LC, LD, LE and allocate angular vars:
    call shgi_set_lcde(LC,LD,LE)

    ! 2c overlap and kinetic dervs:
    call shgi_set_ovrl(LA,LB,LC,LD,LE) ! (4) in shgi_sd_drv

    if(is_on(IADKH))then
      ! prepare and store relativistic projection coefficients:
      ! ( in case any of A or B wasnt touched on this processor before scf ... )
      call shgi_adkh_atom(IL2, uas(IU2)%l_ob(IL2)%exponents &
                             , uas(IU2)%Z                   &
                             , uas(IU2)%nuclear_radius      &
                             , uas(IU2)%s_core%exponents    &
                             , uas(IU2)%s_core%contractions &
                             , uaL                          &
                             )
      call shgi_adkh_atom(IL1, uas(IU1)%l_ob(IL1)%exponents &
                             , uas(IU1)%Z                   &
                             , uas(IU1)%nuclear_radius      &
                             , uas(IU1)%s_core%exponents    &
                             , uas(IU1)%s_core%contractions &
                             , ubL                          &
                             )
      ! note that a==2, b==1 ...
    endif

    if(is_on(IKNSD))then
      ASSERT(is_on(IOVSD))
      ! kinetic energy goes always with overlap:
      call shgi_sd_kin(SDOVRL,SDKNTC)
      ! NOTE: if A==B then gradients sum to zero G(A)+G(B)==0
      !       and need not to be evaluated
      !       (they are set to zero inside of the sub)

      if(is_on(IADKH).and..not.is_on(IAEQB))then
!       write(*,"(' shgi_adkh: DERV <',3I3,'|     KIN     |',3I3,'>')") IU2,IE2,IL2,IU1,IE1,IL1
        ! apply relativistic corrections kinetic energy (and overlap):
        call shgi_adkh_sd_kin(SDOVRL,SDKNTC,uaL,ubL)
        ! NOTE: if A==B then gradients sum to zero G(A)+G(B)==0
        !       and no trafo is necessary
      endif
    endif

    if(is_on(IOVSD))then
       ASSERT(present(OVSD))
       ! note the order of UAs:
       call shgi_sd_abc_store(IU1,IE1,IU2,IE2,0,2,SDOVRL,OVSD, 0) ! 0==overwrite
       ! (0,2) is a special case of 2c-grad into ITS OWN storage.
    endif

    if(is_on(IKNSD))then
       if(present(KNSD))then
         ! note the order of UAs:
         call shgi_sd_abc_store(IU1,IE1,IU2,IE2,0,2,SDKNTC,KNSD, 0) ! 0==overwrite
         ! (0,2) is a special case of 2c-grad into ITS OWN storage.
       else
         ! - ( KIN + NUC ) dervs both go to NUSD:
         ASSERT(present(NUSD))
         call shgi_sd_abc_store(IU2,IE2,IU1,IE1,0,1,SDKNTC,NUSD,+1)
         ! (0,1) is a special case of 2c-grad into COMMON storage.
       endif
    endif

    if(is_on(IPSSD))then
       call shgi_pseu_set_ab() ! just calls bessel_setup()
    endif

    ! == LOOP: UNIQUE ATOMS (C)
    do ua=1,size(uas)
       if(is_on(INUSD+ISRSD))then ! if any of them

          Z = uas(ua)%Z
          ! ZCORE is handled by shgi_pseu0()
          Z = Z - uas(ua)%ZC
          R = radius(ua)
       endif

       if(is_on(IPSSD))then
          LPNUM  = uas(ua)%lmax_pseudo
          ZCORE  = uas(ua)%Z - uas(ua)%ZC
          if( LPNUM > 0 )then
             Z = zero ! disable NUCL
             ! prepare staff common for all EAs:
             call shgi_pseu_set_abc(LPNUM-1,2) ! case=2, i.e. SecDers
          endif
       else
          LPNUM = -1 ! PP is not computed, even if there are some?
       endif

       ! == LOOP: EQUIVALENT ATOMS (C)
       do eq = 1, uas(ua) % N_equal_atoms
          ! set bits to flag coinciding centers: A==B? A==C? B==C?:
          call shgi_set_xeqy(IU2,IE2,IU1,IE1,ua,eq)

          if (whatis (IXEQC) == IXEQC) cycle ! eq atom loop, as gradients are zero

          xc = uas(ua)%position(:,eq)

          ! FXIME: only when precomputing angular part
          !        brings something:
!         call shgi_set_c(xc)

          if(is_on(INUSD,ISRSD))then
             ASSERT(present(NUSD))
             ! both required, do them together:
             SDNUCL = ZERO
             SDNUSR = ZERO
             call shgi_rel_sd_nuc(Z,xc,SDNUCL,SDNUSR,RAD=R)

             if(is_on(IADKH))then
               ! transform and store the dervs of V and SR ints...
!              write(*,"(' shgi_adkh: DERV <',3I3,'|',2I3,F7.1,'|',3I3,'>')") IU2,IE2,IL2,ua,eq,Z,IU1,IE1,IL1
               ASSERT(allocated(SDABC1))
               ASSERT(allocated(SDABC2))

               ! transform grads wrt (W,D) to grads wrt (A,B,C):
               call shgi_sd_wd_to_abc(SDNUCL,SDABC1,0)
               call shgi_sd_wd_to_abc(SDNUSR,SDABC2,0)
               ! NOTE: need to do it to apply rel contractions over
               !       alpha- and beta-exponents, but vector
               !                W = wa*A + wb*B - C
               !       depends on them. Vectors A, B, and C -- dont.

               ! applies rel contraction to dervs of V (in SDABC1)
               ! using dervs of SR ints (in SDABC2), outputs grads of Vrel in-place:
               call shgi_adkh_sd_nuc(SDABC1,SDABC2,uaL,ubL)

               ! rotates and stores ABC-gradients in PG structures:
               call shgi_sd_abc_store(IU2,IE2,IU1,IE1,ua,eq,SDABC1,NUSD,-1)
               ! ... but ignore no more necessary SR ints (in SDABC2)
             else
               ! true DKH: just store the grads of V and SR ints...
               ASSERT(present(SRSD))
               ! rotates and stores derivatives in PG structures:
               call shgi_sd_wd_store(IU2,IE2,IU1,IE1,ua,eq,SDNUCL,NUSD,-1)
               call shgi_sd_wd_store(IU2,IE2,IU1,IE1,ua,eq,SDNUSR,SRSD,-1)
             endif
          else
             ! if any required, do it:
             if(is_on(INUSD))then
                SDNUCL = ZERO
                call shgi_rel_sd_nuc(Z,xc,SDNUCL,RAD=R)
                ! rotates and stores derivatives in PG structures:
                call shgi_sd_wd_store(IU2,IE2,IU1,IE1,ua,eq,SDNUCL,NUSD,-1)
             endif

             if(is_on(ISRSD))then
                SDNUSR = ZERO
                call shgi_rel_sd_nuc(Z,xc,SR=SDNUSR,RAD=R)
                ! rotates and stores derivatives in PG structures:
                call shgi_sd_wd_store(IU2,IE2,IU1,IE1,ua,eq,SDNUSR,SRSD,-1)
             endif
          endif ! is_on(INUSD,INUSD)

          if( LPNUM > 0 )then
             ! do pseudopotential gradients if needed:
             ! w.r.t A, B, and C:
             SDPSEU=ZERO
             call shgi_sd_pseu(ZCORE,LPNUM-1,uas(ua)%l_pseudopot,xc,SDPSEU)

             if(present(PSSD))then
                ! in a realtivistic computation
                ! pseudopotential needs to be separated from NUCL:
                call shgi_sd_abc_store(IU2,IE2,IU1,IE1,ua,eq,SDPSEU,PSSD,-1)
             else
                ! just add it to the NUCL:
                ASSERT(present(NUSD))
                call shgi_sd_abc_store(IU2,IE2,IU1,IE1,ua,eq,SDPSEU,NUSD,-1)
             endif
          endif

       enddo ! loop over equivalent atoms eq

       call shgi_close_abc()
       if( LPNUM > 0 )then
          call shgi_pseu_close_abc()
       endif
    enddo ! loop over unique atoms ua

    ![[=== Add field of PCs to nuclear attraction: ===
    npc = 0
    if( present(pcs) ) npc = size(pcs)
    if( npc > 0 )then

       ! set bits to flag coinciding centers: A==B? A==C? B==C?:
       ! can not use U2==U1 this way:
       call shgi_set_xeqy(IU2,IE2,IU1,IE1,0,0)

       call shgi_set_lcde(1,1,0)
       ! FIXME: do I need to reshape S5?

       ASSERT(is_on(INUSD))
       ASSERT(allocated(SDNUCL))
       ASSERT(present(NUSD))

       SDNUCL = zero
       call shgi_sd_pcs(pcs,SDNUCL)

       if( present(PSSD) )then
         ! put PCSD together with pseudopotential

         ! FIXME: SDPSEU is of different shape 9x9, so use SDNUCL!

         ! rotates and stores ABC-gradients in PG structures:
         call shgi_sd_wd_store(IU2,IE2,IU1,IE1,0,0,SDNUCL,PSSD,-1)
       else
         ! put PCSD together with nuclear attraction

         ! rotates and stores ABC-gradients in PG structures:
         call shgi_sd_wd_store(IU2,IE2,IU1,IE1,0,0,SDNUCL,NUSD,-1)
       endif
    endif
    !]]===============================================

    call shgi_close_ab()
    ! deallocate global storage of THIS module:
    call shgi_glob_free()
    FPP_TIMER_STOP(totD)
    FPP_TIMER_STOP(tot)
    DCALL shgi_timing()
    DPRINT  'SHGI: shgi_sd_drv: exit'

  contains
    function radius(ua) result(r)
      !
      ! returns uas(ua)%nuclear_radius
      ! or zero if the latter is negative
      !
      implicit none
      integer(IK), intent(in) :: ua
      real(RK)                :: r
      ! *** end of interface ***

      ASSERT(ua>0)
      ASSERT(ua<=size(uas))

      ! uas(:) is host-associated from enclosing sub:
      r = uas(ua)%nuclear_radius
      r = max(r,0.0_rk) ! avoid negative
    end function radius
  end subroutine shgi_sd_drv

  subroutine shgi_drv_coul2(IU1,IE1,IL1,IU2,IE2,IL2, CO)
    use unique_atom_module, only: uas=>unique_atoms
    implicit none
    !------------ Declaration of formal parameters ------------------
    integer(IK), intent(in)  :: IU1,IE1,IL1,IU2,IE2,IL2
    real(RK)   , intent(out) :: CO(:,:,:,:)     ! (N2E,N1E,    2*L2+1,2*L1+1  )
    ! *** end of interface ***

    !------------ Declaration of local variables --------------------

    integer(IK) :: typ

    !----------------------------------------------------------------
    !------------ Executable code -----------------------------------

    DPRINT 'SHGI: shgi_drv_coul2: ',IL1,IL2,shape(CO)
    FPP_TIMER_START(tot)
    FPP_TIMER_START(totD)

    if( IL1 >= 0 .and. IL2 >= 0)then
       typ = 0
       call shgi_coul2( typ, IL1, IL2      &
            , uas(IU1)%position(:,IE1)     &
            , uas(IU2)%position(:,IE2)     &
            , uas(IU1)%l_ch(IL1)%exponents &
            , uas(IU2)%l_ch(IL2)%exponents &
            , CO                           &
            )
    endif
    if( IL1 >= 0 .and. IL2 < 0)then
       typ =  1
       call shgi_coul2( typ, IL1, 0        &
            , uas(IU1)%position(:,IE1)     &
            , uas(IU2)%position(:,IE2)     &
            , uas(IU1)%l_ch(IL1)%exponents &
            , uas(IU2)%r2_ch%exponents     &
            , CO                           &
            )
    endif
    if( IL1 < 0 .and. IL2 >= 0)then
       typ = - 1
       call shgi_coul2( typ, 0 , IL2       &
            , uas(IU1)%position(:,IE1)     &
            , uas(IU2)%position(:,IE2)     &
            , uas(IU1)%r2_ch%exponents     &
            , uas(IU2)%l_ch(IL2)%exponents &
            , CO                           &
            )
    endif
    if( IL1 < 0 .and. IL2 < 0)then
       typ = 2
       call shgi_coul2( typ, 0, 0          &
            , uas(IU1)%position(:,IE1)     &
            , uas(IU2)%position(:,IE2)     &
            , uas(IU1)%r2_ch%exponents     &
            , uas(IU2)%r2_ch%exponents     &
            , CO                           &
            )
    endif

    FPP_TIMER_STOP(totD)
    FPP_TIMER_STOP(tot)
  end subroutine shgi_drv_coul2

  subroutine shgi_glob_alloc(NAB,LA,LB)
    ! allocate global storage of THIS module
    ! (check that no mods are USEd!)
    !
    ! KEEP GLOBALS TO MINUMUM!
    ! USE PRIVATE SUBROUTINE VARIABLES WHERE POSSIBLE!
    !
    implicit none
    integer(IK), intent(in) :: NAB,LA,LB
    ! *** end of interface ***

    integer(IK) :: memstat

    if(is_on(IOVRL))then
       allocate(PIOVRL(NAB,2*LA+1,2*LB+1),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IKNTC))then
       allocate(PIKNTC(NAB,2*LA+1,2*LB+1),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(INUCL))then
       allocate(PINUCL(NAB,2*LA+1,2*LB+1),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(INUSR))then
       allocate(PINUSR(NAB,2*LA+1,2*LB+1),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(INUSO))then
       allocate(PINUSO(NAB,2*LA+1,2*LB+1,3),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IPSEU))then
       allocate(PIPSEU(NAB,2*LA+1,2*LB+1),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IOVGR))then
       allocate(GROVRL(NAB,2*LA+1,2*LB+1,3),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IKNGR))then
       allocate(GRKNTC(NAB,2*LA+1,2*LB+1,3),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(INUGR))then
       allocate(GRNUCL(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(ISRGR))then
       allocate(GRNUSR(NAB,2*LA+1,2*LB+1,6),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IPSGR))then
       allocate(GRPSEU(NAB,2*LA+1,2*LB+1,9),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IADKH,INUGR,ISRGR))then
       allocate(GRABC1(NAB,2*LA+1,2*LB+1,9),stat=memstat)
       ASSERT(memstat==0)
       allocate(GRABC2(NAB,2*LA+1,2*LB+1,9),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IADKH,INUSD,ISRSD))then
       allocate(SDABC1(NAB,2*LA+1,2*LB+1,9,9),stat=memstat)
       ASSERT(memstat==0)
       allocate(SDABC2(NAB,2*LA+1,2*LB+1,9,9),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IOVSD))then
       allocate(SDOVRL(NAB,2*LA+1,2*LB+1,3,3),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IKNSD))then
       allocate(SDKNTC(NAB,2*LA+1,2*LB+1,3,3),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(INUSD))then
       allocate(SDNUCL(NAB,2*LA+1,2*LB+1,6,6),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(ISRSD))then
       allocate(SDNUSR(NAB,2*LA+1,2*LB+1,6,6),stat=memstat)
       ASSERT(memstat==0)
    endif

    if(is_on(IPSSD))then
       allocate(SDPSEU(NAB,2*LA+1,2*LB+1,9,9),stat=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine shgi_glob_alloc

  subroutine shgi_glob_free()
    ! no modules used -- only deallocates
    ! global vars of THIS module!
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat

    if( allocated(PIOVRL) )then
       deallocate(PIOVRL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PIKNTC) )then
       deallocate(PIKNTC,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PINUCL) )then
       deallocate(PINUCL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PINUSR) )then
       deallocate(PINUSR,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PINUSO) )then
       deallocate(PINUSO,stat=memstat)
       ASSERT(memstat==0)
    endif

    ! only for GRADS: ------------------>
    if( allocated(GROVRL) )then
       deallocate(GROVRL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(GRKNTC) )then
       deallocate(GRKNTC,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(GRNUCL) )then
       deallocate(GRNUCL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(GRNUSR) )then
       deallocate(GRNUSR,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(GRABC1) )then
       deallocate(GRABC1,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(GRABC2) )then
       deallocate(GRABC2,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(SDABC1) )then
       deallocate(SDABC1,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(SDABC2) )then
       deallocate(SDABC2,stat=memstat)
       ASSERT(memstat==0)
    endif
    !<-----------------------------------

    ! only for SECDER: ----------------->
    if( allocated(SDOVRL) )then
       deallocate(SDOVRL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(SDKNTC) )then
       deallocate(SDKNTC,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(SDNUCL) )then
       deallocate(SDNUCL,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(SDNUSR) )then
       deallocate(SDNUSR,stat=memstat)
       ASSERT(memstat==0)
    endif
    !<-----------------------------------

    ! only for PP: -------------------->
    if( allocated(PIPSEU) )then
       deallocate(PIPSEU,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(GRPSEU) )then
       deallocate(GRPSEU,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(SDPSEU) )then
       deallocate(SDPSEU,stat=memstat)
       ASSERT(memstat==0)
    endif
    !<-----------------------------------
  end subroutine shgi_glob_free

  subroutine shgi_close_abc()
    use shgi_common
    implicit none
    ! *** end of interface ***

    integer(IK) :: memstat

    if( allocated(LAMBDACH) )then
       deallocate(LAMBDACH,NORMCH,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(LAMBDAR2) )then
       deallocate(LAMBDAR2,NORMR2,F132R2,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PICHRG) )then
       deallocate(PICHRG,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PICHSO) )then
       deallocate(PICHSO,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PICHSR) )then
       deallocate(PICHSR,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PIR2RG) )then
       deallocate(PIR2RG,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PIR2SR) )then
       deallocate(PIR2SR,stat=memstat)
       ASSERT(memstat==0)
    endif

    if( allocated(PIR2SO) )then
       deallocate(PIR2SO,stat=memstat)
       ASSERT(memstat==0)
    endif
  end subroutine shgi_close_abc

  subroutine shgi_set_ec(LTYP,ec)
    use shgi_common
    implicit none
    integer(IK), intent(in) :: LTYP
    real(RK)   , intent(in) :: ec(:)
    ! *** end of interface ***

    integer(IK) :: memstat
    integer(IK) :: NEC,iab,ic
    real(RK)    :: eta

    DPRINT 'shgi_set_ec: LTYP=',LTYP,' NEC=',size(ec)

    ! local var:
    NEC = size(ec)

    select case(LTYP)
    case (STYPE)
       allocate(LAMBDACH(NAB,NEC),NORMCH(NAB,NEC),stat=memstat)
       ASSERT(memstat==0)

       do iab=1,NAB
          do ic=1,NEC
             ! eta = ( a + b )/( a + b + c)
             eta = LAMBDA(iab)/( LAMBDA(iab) + ec(ic) )
             ! lambda for fit ints = eta * c
             LAMBDACH(iab,ic) = eta * ec(ic)
             NORMCH(iab,ic)   = TWO * PI * sqrt(eta) / ec(ic) * NORM(iab,2)
          enddo
       enddo

       if(is_on(ICHRG))then
          allocate(PICHRG(NAB,NEC,2*LA+1,2*LB+1),stat=memstat)
          ASSERT(memstat==0)
       endif
       if(is_on(ICHSR))then
          allocate(PICHSR(NAB,NEC,2*LA+1,2*LB+1),stat=memstat)
          ASSERT(memstat==0)
       endif
       if(is_on(ICHSO))then
          allocate(PICHSO(NAB,NEC,2*LA+1,2*LB+1,3),stat=memstat)
          ASSERT(memstat==0)
       endif

    case (R2TYPE)
       allocate(LAMBDAR2(NAB,NEC),NORMR2(NAB,NEC),stat=memstat)
       ASSERT(memstat==0)
       ! only for R2 type:
       allocate(F132R2(NAB,NEC),stat=memstat)
       ASSERT(memstat==0)

       do iab=1,NAB
          do ic=1,NEC
             ! eta = ( a + b )/( a + b + c)
             eta = LAMBDA(iab)/( LAMBDA(iab) + ec(ic) )
             ! lambda for fit ints = eta * c
             LAMBDAR2(iab,ic) = eta * ec(ic)
             NORMR2(iab,ic)   = TWO * PI * eta**(THREE/TWO) / ec(ic)**2 * NORM(iab,2)
             F132R2(iab,ic)   = ONE + (THREE/TWO) * ec(ic) / LAMBDA(iab)
          enddo
       enddo

       if(is_on(ICHRG))then
          allocate(PIR2RG(NAB,NEC,2*LA+1,2*LB+1),stat=memstat)
          ASSERT(memstat==0)
       endif
       if(is_on(ICHSR))then
          allocate(PIR2SR(NAB,NEC,2*LA+1,2*LB+1),stat=memstat)
          ASSERT(memstat==0)
       endif
       if(is_on(ICHSO))then
          allocate(PIR2SO(NAB,NEC,2*LA+1,2*LB+1,3),stat=memstat)
          ASSERT(memstat==0)
       endif

    case default
       ABORT('LTYP')
    end select
  end subroutine shgi_set_ec

  subroutine shgi_coul2(typ,LA,LB,xa,xb,ea,eb,co)
    use shgi_common, only: XD,XD2,ZETA,NORM,WDA,WDB,LAMBDA
    use shgi_rad, only: doQSSL, doQSRL, doQRRL
    use shgi_ab,  only: shgi_set_ab, shgi_close_ab
    use shgi_shr, only: SHR_D2Av
    use shgi_utils , only: radang, upack2c
    implicit none
    integer(IK), intent(in)    :: typ
    integer(IK), intent(in)    :: LA,LB ! shadow globals
    real(RK)   , intent(in)    :: xa(:),xb(:) ! xa(3), xb(3)
    real(RK)   , intent(in)    :: ea(:),eb(:) ! ea(NEXPA), eb(NEXPB)
    real(RK)   , intent(inout) :: co(:,:,:,:) ! (NA,NB,2*LA+1,2*LB+1)
    ! *** end of interface ***

    real(RK), allocatable  :: IL(:,:)     ! (NAB,1+LA+LB)
    real(RK), allocatable  :: X2(:,:,:)   ! (    2*LA+1,2*LB+1,1+LA+LB)
    real(RK), allocatable  :: pckd(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    integer(IK), parameter :: &
         SS =  0, &
         SR =  1, &
         RS = -1, &
         RR =  2

    call setif(-1_i8k,.false.) ! zero all bits

    ! For long-range coulomb "overlap" integrals dont use screening:
    call shgi_set_maxexp(HUGE(1.0_rk))

    ! copy exponents, precompute 2c norms...:
    call shgi_set_ab(LA,LB,xa,xb,ea,eb)
    ! NOT USED: call shgi_glob_alloc(NAB,LA,LB)

    ! NAB was set by masking NAEXP x NBEXP in shgi_set_ab:
    allocate( IL(NAB,1+LA+LB)                 &
            , X2((LA+1)**2,(LB+1)**2,1+LA+LB) &
            , pckd(NAB,2*LA+1,2*LB+1)         &
            )

    ! compute angluar part of D = ( A - B ):
    call SHR_D2Av(1,LA,LB,XD,X2)

    ! compute radial part of D^2 = ( A - B )^2:
    select case(typ)
    case (SS)
       call doQSSL(LA+LB,XD2,ZETA,NORM(:,4)          , IL)
    case (SR)
       call doQSRL(LA+LB,XD2,ZETA,NORM(:,4)/LAMBDA   , IL, WDA, WDB, SR)
    case (RS)
       call doQSRL(LA+LB,XD2,ZETA,NORM(:,4)/LAMBDA   , IL, WDA, WDB, RS)
    case (RR)
       call doQRRL(LA+LB,XD2,ZETA,NORM(:,4)/LAMBDA**2, IL, WDA, WDB)
    case default
       ABORT('no such type')
    end select

    ! couple radial and angular parts:
    call radang(LA,LB,IL,X2,pckd)

    call upack2c(pckd,co)

    deallocate(IL,pckd,X2)

    call shgi_close_ab()
    ! NOT USED: call shgi_glob_free()
  end subroutine shgi_coul2

  !--------------- End of module ----------------------------------
end module shgi

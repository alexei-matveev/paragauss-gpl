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
module shgi_utils
  !-------------------------------------------------------------------
  !
  ! Coupling harmonic derivarives of F and  S to those of FS, also for
  ! gradients,  also  for  second  derivatives.  Coupling  radial  and
  ! angular parts. Conversion from harmonic derivatives wrt W and D to
  ! those wrt  A, B, and C. Packing/unpacking  screened integrals. And
  ! more ...
  !
  ! Kitchen sink for everithing without a better place to go.
  !
  ! Copyright (c) 2005-2013 Alexei Matveev
  ! Copyright (c) 2006 Vladimir Nasluzov
  ! Copyright (c) 2006-2008 Alexey Shor
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
! use CPU_TIME for timers:
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       IK=>i4_kind, RK=>r8_kind ! type specification parameters
  use constants
  use shgi_cntrl
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module ====================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  !------------ Interface statements ---------------------------------

  !------------ public functions and subroutines ---------------------
  public :: shellNum!(ua,L) returns linear index of (ua,L) shell
  public :: shgi_dww_to_dab
  public :: grFS
  public :: grFS1
  public :: sdFS, sdSYM
  public :: doRadAng
  public :: radang
  public :: upack2c
  public :: upack3c
  public :: shgi_gr_wd_store
  public :: shgi_ve_abc_store
  public :: shgi_gr_abc_store
  public :: shgi_gr_wd_to_ab
  public :: shgi_gr_wd_to_abc
  public :: shgi_sd_wd_to_abc
  public :: shgi_sd_wd_store

  public :: shgi_fl_wd_store
  public :: shgi_X_wd_store

  public :: shgi_slvd_abc_store
  public :: shgi_slvd_abc_store1
  public :: shgi_slvd_abc_store2

  public :: shgi_vd_abc_store
  public :: shgi_sd_abc_store

  public :: shgi_timing

  !===================================================================
  ! End of public interface of module
  !===================================================================


  !------------ Declaration of types ---------------------------------

  !------------ Declaration of constants and variables ---------------

  ! ALL INTEGER CONSTANTS, LIKE GAX,GAY,GAZ ARE NOW IN
  !                     shgi_cntrl.f90
  ! THIS HAS BEEN DONE TO ENABLE ITS USE IN OTHER MODULES

  ! MOST GLOBAL VARIABLES HOLDING ``ANGULAR'' FACTORS WERE MOVED TO
  !                     shgi_common.f90
  ! THIS HAS BEEN DONE TO SPLIT THIS FILE INTO PARTS LATER

  ! a copy from shgi_shr.f90:
  integer(IK)   :: L_,M_
  integer(IK), parameter :: MAXL = 6 ! s,p,d,f,g,h,i
  integer(IK), parameter :: lof( (MAXL+1)**2 ) = (/((L_,M_=1,2*L_+1),L_=0,MAXL)/)
  integer(IK), parameter :: mof( (MAXL+1)**2 ) = (/((M_,M_=1,2*L_+1),L_=0,MAXL)/)

  !------------ Subroutines ------------------------------------------
contains

  function shellNum(ua,L) result(uaL)
    ! returns cumulative index of (ua,L) shell,
    ! if called with wheird arguments, returns
    ! the total number of shells
    use unique_atom_module, only: uas => unique_atoms!(:)%lmax_ob
    implicit none
    integer(IK), intent(in) :: ua,L
    integer(IK)             :: uaL ! result
    ! *** end of interface ***

    integer(IK) :: i,j

    uaL = 0
    do i=1,size(uas)
      do j=0,uas(i)%lmax_ob
        uaL = uaL + 1
        if( i==ua .and. j==L ) RETURN
      enddo
    enddo
  end function shellNum

  subroutine grFS(F,S,gr,fact)
    ! Computes GW and GD gradients of
    !         F(W^2/2)*S(D^2/2)
    ! w.r.t. to W and D
    implicit none
    real(RK), intent(in)    :: F(:,:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2,4)
    real(RK), intent(in)    :: S(:,:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2,4)
    real(RK), intent(inout) :: gr(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), intent(in)    :: fact(:)     ! (NAB)
    optional :: fact
    ! *** end of interface **

    logical :: prefact

    prefact = present(fact)

    ! F(g)*S(0):
    call PR2(F(:,:,:,CX),S(:,:,:,C0),gr(:,:,:,GWX))
    call PR2(F(:,:,:,CY),S(:,:,:,C0),gr(:,:,:,GWY))
    call PR2(F(:,:,:,CZ),S(:,:,:,C0),gr(:,:,:,GWZ))

    ! F(0)*S(g):
    call PR2(F(:,:,:,C0),S(:,:,:,CX),gr(:,:,:,GDX))
    call PR2(F(:,:,:,C0),S(:,:,:,CY),gr(:,:,:,GDY))
    call PR2(F(:,:,:,C0),S(:,:,:,CZ),gr(:,:,:,GDZ))
  contains
    subroutine PR2(F,S,fs)
      use shgi_shr, only: SHR_PR2v
      implicit none
      real(RK), intent(in)    :: F(:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2)
      real(RK), intent(in)    :: S(:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2)
      real(RK), intent(inout) :: fs(:,:,:) ! (NAB,2*LA+1,2*LB+1)
      ! *** end of interface ***

      if(.not.prefact)then
        call SHR_PR2v(NAB,LA,LB,F(:,:,:),S(:,:,:),fs(:,:,:))
      else
        call SHR_PR2v(NAB,LA,LB,F(:,:,:),S(:,:,:),fs(:,:,:), fact=fact)
      endif
    end subroutine PR2
  end subroutine grFS

  subroutine grFS1(F,S,gr)
    ! Computes GW and GD gradients of
    !         F(W^2/2)*S(D^2/2)
    ! w.r.t. to W and D
    implicit none
    real(RK), intent(in)    :: F(:,:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2,4)
    real(RK), intent(in)    :: S(:,:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2,4)
    real(RK), intent(inout) :: gr(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    ! *** end of interface **


    ! F(g)*S(0):
    call PR2(F(:,:,:,CX),S(:,:,:,C0),gr(:,:,:,GWX))
    call PR2(F(:,:,:,CY),S(:,:,:,C0),gr(:,:,:,GWY))
    call PR2(F(:,:,:,CZ),S(:,:,:,C0),gr(:,:,:,GWZ))

  contains
    subroutine PR2(F,S,fs)
      use shgi_shr, only: SHR_PR2v
      implicit none
      real(RK), intent(in)    :: F(:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2)
      real(RK), intent(in)    :: S(:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2)
      real(RK), intent(inout) :: fs(:,:,:) ! (NAB,2*LA+1,2*LB+1)
      ! *** end of interface ***

      call SHR_PR2v(NAB,LA,LB,F(:,:,:),S(:,:,:),fs(:,:,:))

    end subroutine PR2
  end subroutine grFS1

  subroutine sdFS(F,S,sd,fact)
    ! Computes GW and GD derivatives of
    !         F(W^2/2)*S(D^2/2)
    ! w.r.t. to W and D
    ! DONT FORGET TO CALL sdSYM(sd)!
    use shgi_shr, only: SHR_PR2v
    implicit none
    real(RK), intent(in)    :: F(:,:,:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2,4,4)
    real(RK), intent(in)    :: S(:,:,:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2,4,4)
    real(RK), intent(inout) :: sd(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK), intent(in)    :: fact(:)       ! (NAB)
    optional :: fact
    ! *** end of interface **

    logical :: prefact

    prefact = present(fact)

    ! computing triangle sd(i,j), i >= j
    !      w w w   d d d
    !      x y z   x y z
    !
    ! w x  1
    ! w y  2 4
    ! w z  3 5 6

    ! d x  7 0 3   6
    ! d y  8 1 4   7 9
    ! d z  9 2 5   8 0 1

    ! F(i,j)*S(0):                                   W   W
    call PR2(F(:,:,:,CX,CX),S(:,:,:,C0,C0),sd(:,:,:,GWX,GWX)) ! 1
    call PR2(F(:,:,:,CY,CX),S(:,:,:,C0,C0),sd(:,:,:,GWY,GWX)) ! 2
    call PR2(F(:,:,:,CZ,CX),S(:,:,:,C0,C0),sd(:,:,:,GWZ,GWX)) ! 3

    call PR2(F(:,:,:,CY,CY),S(:,:,:,C0,C0),sd(:,:,:,GWY,GWY)) ! 4
    call PR2(F(:,:,:,CZ,CY),S(:,:,:,C0,C0),sd(:,:,:,GWZ,GWY)) ! 5

    call PR2(F(:,:,:,CZ,CZ),S(:,:,:,C0,C0),sd(:,:,:,GWZ,GWZ)) ! 6

    ! F(i,0)*S(0,j):                                 D   W
    call PR2(F(:,:,:,CX,C0),S(:,:,:,CX,C0),sd(:,:,:,GDX,GWX)) ! 7
    call PR2(F(:,:,:,CX,C0),S(:,:,:,CY,C0),sd(:,:,:,GDY,GWX)) ! 8
    call PR2(F(:,:,:,CX,C0),S(:,:,:,CZ,C0),sd(:,:,:,GDZ,GWX)) ! 9

    call PR2(F(:,:,:,CY,C0),S(:,:,:,CX,C0),sd(:,:,:,GDX,GWY)) !10
    call PR2(F(:,:,:,CY,C0),S(:,:,:,CY,C0),sd(:,:,:,GDY,GWY)) !11
    call PR2(F(:,:,:,CY,C0),S(:,:,:,CZ,C0),sd(:,:,:,GDZ,GWY)) !12

    call PR2(F(:,:,:,CZ,C0),S(:,:,:,CX,C0),sd(:,:,:,GDX,GWZ)) !13
    call PR2(F(:,:,:,CZ,C0),S(:,:,:,CY,C0),sd(:,:,:,GDY,GWZ)) !14
    call PR2(F(:,:,:,CZ,C0),S(:,:,:,CZ,C0),sd(:,:,:,GDZ,GWZ)) !15

    ! F(0)*S(i,j):                                   D   D
    call PR2(F(:,:,:,C0,C0),S(:,:,:,CX,CX),sd(:,:,:,GDX,GDX)) !16
    call PR2(F(:,:,:,C0,C0),S(:,:,:,CY,CX),sd(:,:,:,GDY,GDX)) !17
    call PR2(F(:,:,:,C0,C0),S(:,:,:,CZ,CX),sd(:,:,:,GDZ,GDX)) !18

    call PR2(F(:,:,:,C0,C0),S(:,:,:,CY,CY),sd(:,:,:,GDY,GDY)) !19
    call PR2(F(:,:,:,C0,C0),S(:,:,:,CZ,CY),sd(:,:,:,GDZ,GDY)) !20

    call PR2(F(:,:,:,C0,C0),S(:,:,:,CZ,CZ),sd(:,:,:,GDZ,GDZ)) !21
    ! 6 * ( 6 + 1 ) / 2 = 21

    ! DONT FORGET TO CALL sdSYM(sd)!
  contains
    subroutine PR2(F,S,fs)
      use shgi_shr, only: SHR_PR2v
      implicit none
      real(RK), intent(in)    :: F(:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2)
      real(RK), intent(in)    :: S(:,:,:)  ! (NAB,(LA+1)**2,(LB+1)**2)
      real(RK), intent(inout) :: fs(:,:,:) ! (NAB,2*LA+1,2*LB+1)
      ! *** end of interface ***

      if(.not.prefact)then
        call SHR_PR2v(NAB,LA,LB,F(:,:,:),S(:,:,:),fs(:,:,:))
      else
        call SHR_PR2v(NAB,LA,LB,F(:,:,:),S(:,:,:),fs(:,:,:), fact=fact)
      endif
    end subroutine PR2
  end subroutine sdFS

  subroutine sdSYM(sd)
    ! suncs lower and upper triangles
    implicit none
    real(RK), intent(inout) :: sd(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    ! *** end of interface **

    ! copy computed triangle sd(i,j), i >= j
    !                     to sd(j,i)
    !      w w w   d d d
    !      x y z   x y z
    !
    ! w x  1
    ! w y  2 4
    ! w z  3 5 6

    ! d x  7 0 3   6
    ! d y  8 1 4   7 9
    ! d z  9 2 5   8 0 1

    sd(:,:,:,GWX,GWY) = sd(:,:,:,GWY,GWX) !  2
    sd(:,:,:,GWX,GWZ) = sd(:,:,:,GWZ,GWX) !  3
    sd(:,:,:,GWY,GWZ) = sd(:,:,:,GWZ,GWY) !  5

    sd(:,:,:,GWX,GDX) = sd(:,:,:,GDX,GWX) !  7
    sd(:,:,:,GWX,GDY) = sd(:,:,:,GDY,GWX) !  8
    sd(:,:,:,GWX,GDZ) = sd(:,:,:,GDZ,GWX) !  9
    sd(:,:,:,GWY,GDX) = sd(:,:,:,GDX,GWY) ! 10
    sd(:,:,:,GWY,GDY) = sd(:,:,:,GDY,GWY) ! 11
    sd(:,:,:,GWY,GDZ) = sd(:,:,:,GDZ,GWY) ! 12
    sd(:,:,:,GWZ,GDX) = sd(:,:,:,GDX,GWZ) ! 13
    sd(:,:,:,GWZ,GDY) = sd(:,:,:,GDY,GWZ) ! 14
    sd(:,:,:,GWZ,GDZ) = sd(:,:,:,GDZ,GWZ) ! 15

    sd(:,:,:,GDX,GDY) = sd(:,:,:,GDY,GDX) ! 17
    sd(:,:,:,GDX,GDZ) = sd(:,:,:,GDZ,GDX) ! 18
    sd(:,:,:,GDY,GDZ) = sd(:,:,:,GDZ,GDY) ! 20
    ! 6 * 6 = 36 = 21 + 15
  end subroutine sdSYM

  subroutine shgi_dww_to_dab(NAB,LA,LB,NLM,DF)
    !
    !   D(lma)D(lmb)... F => DA(lma)DB(lmb)... F
    !
    use shgi_common, only: WDAL, WDBL
    implicit none
    integer(IK), intent(in)    :: NAB,LA,LB,NLM
    real(RK)   , intent(inout) :: DF(NAB,(LA+1)**2,(LB+1)**2,NLM)
    ! *** end of interface ***

    integer(IK) :: lm
    integer(IK) :: lma,ila
    integer(IK) :: lmb,ilb

    do lm =1,NLM
    do lmb=1,(LB+1)**2
       ilb=lof(lmb)
    do lma=1,(LA+1)**2
       ila=lof(lma)

       DF(:,lma,lmb,lm) = DF(:,lma,lmb,lm) &
                        * WDAL(:,1+ila) * WDBL(:,1+ilb)
    enddo
    enddo
    enddo
  end subroutine shgi_dww_to_dab

  subroutine doRadAng(LL,IL,Y,NUC,fact)
    ! couples radial and "angular"
    ! adds integral to NUC
    ! DONT FORGET TO CLEAR NUC!
    implicit none
    integer(IK), intent(in) :: LL
    real(RK), intent(in)    :: IL(:,:)    ! (NAB,              1+LL)
    real(RK), intent(in)    :: Y(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,1+LL)
    real(RK), intent(inout) :: NUC(:,:,:) ! (NAB,2*LA+1,2*LB+1)
    real(RK), intent(in)    :: fact(:)    ! (NAB)
    optional :: fact
    ! *** end of interface ***

    integer(IK) :: ma,mb,L

    FPP_TIMER_START(tra)

    ASSERT(NAB==size(Y,1))
    ASSERT(size(IL,1)==size(Y,1))
    ASSERT(size( IL,1)==size(NUC,1))
    ASSERT(size( IL,2)>=1+LL)
    ASSERT(size(  Y,4)>=1+LL)
    ASSERT(size(  Y,2)==2*LA+1)
    ASSERT(size(  Y,3)==2*LB+1)
    ASSERT(size(NUC,2)==2*LA+1)
    ASSERT(size(NUC,3)==2*LB+1)

    if( present(fact) )then
       ASSERT(NAB==size(fact))
       do mb=1,2*LB+1
          do ma=1,2*LA+1
             do L=1,1+LL
                NUC(:,ma,mb) =  NUC(:,ma,mb) &
                     + fact(:) * IL(:,L) * Y(:,ma,mb,L)
             enddo
          enddo
       enddo
    else
       do mb=1,2*LB+1
          do ma=1,2*LA+1
             do L=1,1+LL
                NUC(:,ma,mb) =  NUC(:,ma,mb) &
                     + IL(:,L) * Y(:,ma,mb,L)
             enddo
          enddo
       enddo
    endif

    FPP_TIMER_STOP(tra)
  end subroutine doRadAng

  subroutine radang(LA,LB,SL,X,pckd)
    implicit none
    integer(IK),intent(in) :: LA,LB
    real(RK), intent(in)   :: SL(:,:)
    real(RK), intent(in)   :: X(:,:,:)
    real(RK), intent(out)  :: pckd(:,:,:)
    ! *** end of interface ***

    integer(IK) :: ma,mb,lma,lmb,L

    ASSERT(1+LA+LB<=size(SL,2))

    do mb=1,2*LB+1
       lmb = LB**2 + mb
       do ma=1,2*LA+1
          lma = LA**2 + ma
          ! FIXME: nullified:
          pckd(:,ma,mb) = ZERO
          ! FIXME: range:
          !do L=1,1+LA+LB
          do L=1+max(LA,LB),1+LA+LB
             pckd(:,ma,mb) =  pckd(:,ma,mb) &
                  + (-1)**LB * SL(:,L) * X(lma,lmb,L)
          enddo
       enddo
    enddo
  end subroutine radang

  subroutine shgi_gr_wd_to_ab(gwd,gab,iop)
    !
    ! converts gradients ``gwd'' w.r.t. to
    !
    !         W = wa * A + wb * B
    !         D = A - B
    !
    ! to gradients w.r.t. A, B adding them to ``gab''
    !
    ! GR(ag) :   wa    *F(g)*S(0) + F(0)*S(g)
    ! GR(bg) :      wb *F(g)*S(0) - F(0)*S(g)
    !                   \_ gw _/    \_ gd _/
    !
    use shgi_cntrl , only: LA, LB, NAB
    use shgi_common, only: WDA, WDB
    implicit none
    real(RK), intent(in)    :: gwd(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), intent(inout) :: gab(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    integer(IK), intent(in) :: iop ! +1,0,-1
    ! *** end of interface **

    integer(IK) :: ab,ma,mb

    select case (iop)
    case(0)

  !   FORALL( ab=1:NAB, ma=1:2*LA+1, mb=1:2*LB+1 )
      do mb=1,2*LB+1
      do ma=1,2*LA+1
      do ab=1,NAB
         gab(ab,ma,mb,GAX) = + gwd(ab,ma,mb,GWX) * WDA(ab) + gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GAY) = + gwd(ab,ma,mb,GWY) * WDA(ab) + gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GAZ) = + gwd(ab,ma,mb,GWZ) * WDA(ab) + gwd(ab,ma,mb,GDZ)

         gab(ab,ma,mb,GBX) = + gwd(ab,ma,mb,GWX) * WDB(ab) - gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GBY) = + gwd(ab,ma,mb,GWY) * WDB(ab) - gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GBZ) = + gwd(ab,ma,mb,GWZ) * WDB(ab) - gwd(ab,ma,mb,GDZ)
      enddo
      enddo
      enddo
  !   END FORALL

    case(+1)

  !   FORALL( ab=1:NAB, ma=1:2*LA+1, mb=1:2*LB+1 )
      do mb=1,2*LB+1
      do ma=1,2*LA+1
      do ab=1,NAB
         gab(ab,ma,mb,GAX) = gab(ab,ma,mb,GAX) + gwd(ab,ma,mb,GWX) * WDA(ab) + gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GAY) = gab(ab,ma,mb,GAY) + gwd(ab,ma,mb,GWY) * WDA(ab) + gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GAZ) = gab(ab,ma,mb,GAZ) + gwd(ab,ma,mb,GWZ) * WDA(ab) + gwd(ab,ma,mb,GDZ)

         gab(ab,ma,mb,GBX) = gab(ab,ma,mb,GBX) + gwd(ab,ma,mb,GWX) * WDB(ab) - gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GBY) = gab(ab,ma,mb,GBY) + gwd(ab,ma,mb,GWY) * WDB(ab) - gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GBZ) = gab(ab,ma,mb,GBZ) + gwd(ab,ma,mb,GWZ) * WDB(ab) - gwd(ab,ma,mb,GDZ)
      enddo
      enddo
      enddo
  !   END FORALL
    case default
      ABORT('FIXME: implement -1')
    end select
  end subroutine shgi_gr_wd_to_ab

  subroutine shgi_gr_wd_to_abc(gwd,gab,iop)
    !
    ! converts gradients ``gwd'' w.r.t. to
    !
    !         W = wa * A + wb * B - C
    !         D = A - B
    !
    ! to gradients w.r.t. A, B, and C adding them to ``gab''
    !
    ! GR(cg) : -(wa+wb)*F(g)*S(0)             , wa+wb=1, g=x,y,z
    ! GR(ag) :   wa    *F(g)*S(0) + F(0)*S(g)
    ! GR(bg) :      wb *F(g)*S(0) - F(0)*S(g)
    !                   \_ gw _/    \_ gd _/
    !
    use shgi_cntrl , only: LA, LB, NAB
    use shgi_common, only: WDA, WDB
    implicit none
    real(RK), intent(in)    :: gwd(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), intent(inout) :: gab(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    integer(IK), intent(in) :: iop ! +1,0,-1
    ! *** end of interface **

    integer(IK) :: ab,ma,mb

    select case (iop)
    case(0)
      ! overwrite ...
      gab(:,:,:,GCX)       = - gwd(:,:,:,GWX)  ! *(WDA+WDB)
      gab(:,:,:,GCY)       = - gwd(:,:,:,GWY)  ! *(WDA+WDB)
      gab(:,:,:,GCZ)       = - gwd(:,:,:,GWZ)  ! *(WDA+WDB)

  !   FORALL( ab=1:NAB, ma=1:2*LA+1, mb=1:2*LB+1 )
      do mb=1,2*LB+1
      do ma=1,2*LA+1
      do ab=1,NAB
         gab(ab,ma,mb,GAX) = + gwd(ab,ma,mb,GWX) * WDA(ab) + gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GAY) = + gwd(ab,ma,mb,GWY) * WDA(ab) + gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GAZ) = + gwd(ab,ma,mb,GWZ) * WDA(ab) + gwd(ab,ma,mb,GDZ)

         gab(ab,ma,mb,GBX) = + gwd(ab,ma,mb,GWX) * WDB(ab) - gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GBY) = + gwd(ab,ma,mb,GWY) * WDB(ab) - gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GBZ) = + gwd(ab,ma,mb,GWZ) * WDB(ab) - gwd(ab,ma,mb,GDZ)
      enddo
      enddo
      enddo
  !   END FORALL

    case(+1)
      ! add ...
      gab(:,:,:,GCX)       = gab(:,:,:,GCX)    - gwd(:,:,:,GWX)  ! *(WDA+WDB)
      gab(:,:,:,GCY)       = gab(:,:,:,GCY)    - gwd(:,:,:,GWY)  ! *(WDA+WDB)
      gab(:,:,:,GCZ)       = gab(:,:,:,GCZ)    - gwd(:,:,:,GWZ)  ! *(WDA+WDB)

  !   FORALL( ab=1:NAB, ma=1:2*LA+1, mb=1:2*LB+1 )
      do mb=1,2*LB+1
      do ma=1,2*LA+1
      do ab=1,NAB
         gab(ab,ma,mb,GAX) = gab(ab,ma,mb,GAX) + gwd(ab,ma,mb,GWX) * WDA(ab) + gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GAY) = gab(ab,ma,mb,GAY) + gwd(ab,ma,mb,GWY) * WDA(ab) + gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GAZ) = gab(ab,ma,mb,GAZ) + gwd(ab,ma,mb,GWZ) * WDA(ab) + gwd(ab,ma,mb,GDZ)

         gab(ab,ma,mb,GBX) = gab(ab,ma,mb,GBX) + gwd(ab,ma,mb,GWX) * WDB(ab) - gwd(ab,ma,mb,GDX)
         gab(ab,ma,mb,GBY) = gab(ab,ma,mb,GBY) + gwd(ab,ma,mb,GWY) * WDB(ab) - gwd(ab,ma,mb,GDY)
         gab(ab,ma,mb,GBZ) = gab(ab,ma,mb,GBZ) + gwd(ab,ma,mb,GWZ) * WDB(ab) - gwd(ab,ma,mb,GDZ)
      enddo
      enddo
      enddo
  !   END FORALL
    case default
      ABORT('FIXME: implement -1')
    end select
  end subroutine shgi_gr_wd_to_abc

  subroutine shgi_sd_wd_to_abc(gwd,gab,iop)
    !
    ! converts second derivatives ``gwd'' to ``gab''
    !
    use shgi_cntrl , only: LA, LB, NAB
    implicit none
    real(RK), intent(in)    :: gwd(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    real(RK), intent(inout) :: gab(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9,9)
    integer(IK), intent(in) :: iop ! +1,0,-1
    ! *** end of interface **

    real(RK)    :: gxx(NAB,2*LA+1,2*LB+1,9,6) ! temp storage
    integer(IK) :: i

    do i=1,6
    call shgi_gr_wd_to_abc(gwd(:,:,:,:,i),gxx(:,:,:,:,i),  0)
    enddo
    do i=1,9
    call shgi_gr_wd_to_abc(gxx(:,:,:,i,:),gab(:,:,:,i,:),iop)
    enddo
  end subroutine shgi_sd_wd_to_abc

  subroutine shgi_gr_wd_store(UA,EA,UB,EB,UC,EC,GWD,STO,iop)
    ! Convert W,D derivatives to A,B,C derivatives and store
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in) :: UA,EA
    integer(IK), intent(in) :: UB,EB
    integer(IK), intent(in) :: UC,EC
    real(RK)   , intent(in) :: GWD(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    type(arrmat4) , intent(inout) :: STO(:)   ! (num of tot-sym grads)
    integer(IK), intent(in) :: iop
    ! *** end of interface ***

    real(RK) :: GAB(NAB,2*LA+1,2*LB+1,9)

    call shgi_gr_wd_to_abc(GWD,GAB,0)
    call shgi_gr_abc_store(UA,EA,UB,EB,UC,EC,GAB,STO,iop)
  end subroutine shgi_gr_wd_store

  subroutine shgi_sd_wd_store(UA,EA,UB,EB,UC,EC,GWD,STO,iop)
    ! Convert W,D derivatives to A,B,C derivatives and store
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in) :: UA,EA
    integer(IK), intent(in) :: UB,EB
    integer(IK), intent(in) :: UC,EC
    real(RK)   , intent(in) :: GWD(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6,6)
    type(arrmat4) , intent(inout) :: STO(:,:)   ! (num of tot-sym grads^2)
    integer(IK), intent(in) :: iop
    ! *** end of interface ***

    real(RK) :: GAB(NAB,2*LA+1,2*LB+1,9,9)

    call shgi_sd_wd_to_abc(GWD,GAB,0)
    call shgi_sd_abc_store(UA,EA,UB,EB,UC,EC,GAB,STO,iop)
  end subroutine shgi_sd_wd_store

  subroutine shgi_X_wd_store(Xc,uc,ec,GWD,STO,iop)
    ! X centers (PC, PD, PQ, PO, RC) grads - Convert W derivatives to C derivatives and store
    use datatype, only: arrmat4 ! PG storage of grs
    use shgi_sym, only: shgi_sym_X_rot
    implicit none

    real(RK)   , intent(in) :: GWD(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    type(arrmat4) , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in) :: Xc,uc,ec,iop
    ! *** end of interface ***

    real(RK) :: GAB(NAB,2*LA+1,2*LB+1,9)
    integer(IK) :: ngr ! number of tot sym grads
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3) ! only ..,:ngr) is defined

    GAB = 0.0_rk
    call shgi_fl_wd_to_abc(GWD,GAB)
    call shgi_sym_X_rot(Xc,uc,ec,GAB(:,:,:,GCX:GCZ),SYGR,ngr) ! returns ngr
    call shgi_sym_X_add(Xc,uc,SYGR,ngr,STO,iop) ! uses ngr
  end subroutine shgi_X_wd_store

  subroutine shgi_fl_wd_store(uc,ec,GWD,STO,iop)
    ! Electrostatic field - Convert W derivatives to C derivatives and store
    use shgi_sym, only: shgi_sym_fl_rot
    implicit none

    real(RK)   , intent(in) :: GWD(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK)   , intent(inout) :: STO(:,:,:,:,:)
    !                                     ^ - num of tot-sym grads
    integer(IK), intent(in) :: uc,ec,iop
    ! *** end of interface ***

    real(RK) :: GAB(NAB,2*LA+1,2*LB+1,9)
    integer(IK) :: ngr ! number of tot sym grads
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3) ! only ..,:ngr) is defined

    GAB = 0.0_rk
    call shgi_fl_wd_to_abc(GWD,GAB)
    call shgi_sym_fl_rot(uc,ec,GAB(:,:,:,GCX:GCZ),SYGR,ngr) ! returns ngr
    call shgi_sym_fl_add(uc,SYGR,ngr,STO,iop) ! uses ngr
  end subroutine shgi_fl_wd_store

  subroutine shgi_fl_wd_to_abc(gwd,gab)
    implicit none
    real(RK), intent(in)    :: gwd(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,6)
    real(RK), intent(inout) :: gab(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    ! *** end of interface **

    gab(:,:,:,GCX)       = gab(:,:,:,GCX)    - gwd(:,:,:,GWX)  ! *(WDA+WDB)
    gab(:,:,:,GCY)       = gab(:,:,:,GCY)    - gwd(:,:,:,GWY)  ! *(WDA+WDB)
    gab(:,:,:,GCZ)       = gab(:,:,:,GCZ)    - gwd(:,:,:,GWZ)  ! *(WDA+WDB)
  end subroutine shgi_fl_wd_to_abc

  subroutine shgi_ve_abc_store(UA,EA,UB,EB,UC,EC,GR,STO,iop)
    use solv_cavity_module, only: VTN
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    integer(IK), intent(in)    :: UC,EC
    real(RK)   , intent(in)    :: GR(:)  ! (9) or (6)
    real(RK)   , intent(inout) :: STO(:) ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    if( UC == 0 .and. EC /= 0 )then
      ! keep special cases to minimum:
      ABORT('must die')
    endif

    ! The code below also handles the case UC=0 as if there is no C
    ! (e.g. point charges)
    ASSERT(size(GR)>=6)

    if(UA>0) call shgi_sym_ve_sto(UA,EA,GR(GAX:GAZ),STO,iop)
    if(UB>0) call shgi_sym_ve_sto(UB,EB,GR(GBX:GBZ),STO,iop)
!   if(UC>0) call shgi_sym_ve_sto(UC,EC,GR(GCX:GCZ),STO,iop)
    if(UC>0)then
      ASSERT(size(GR)==9)
      if(.not.is_on(ISLGR) .and. .not.is_on(ISLGRT))then
             call shgi_sym_ve_sto(UC,EC,GR(GCX:GCZ),STO,iop)
      else
         if(VTN) then
            call shgi_sym_ve_sto(UC,EC,GR(GCX:GCZ),STO,iop)
         else
            ! NOTE: symmetry type is passed as EC=ISM:
            call shgi_sym_slve_sto(UC,EC,GR(GCX:GCZ),STO,iop)
         end if
      endif
    endif
  end subroutine shgi_ve_abc_store

  subroutine shgi_gr_abc_store(UA,EA,UB,EB,UC,EC,GR,STO,iop)
    use datatype, only: arrmat4 ! PG storage of grs
    use shgi_cntrl, only: is_on,ISLGR,ISLGRT
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    integer(IK), intent(in)    :: UC,EC
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9)
    type(arrmat4) , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

!    integer(IK) :: i
!    do i=1,6
!     STO(i)%m=0.0_rk
!    enddo



    ! Rotate the C-gradient along the directions
    ! of (totally) symmetric modes:

    if( UC == 0 .and. EC == 1 )then
      ! special case: Ga(3) and Gb(3) = - Ga(3) to be added
      ! to the total grad:
      call shgi_gr_ab1_store(UA,EA,UB,EB,GR,STO,iop)
      RETURN
    endif

    if( UC == 0 .and. EC == 2 )then
      ! special case: Ga(3) and Gb(3) = - Ga(3) to be stored
      ! in its own storage of max. 6 grads (MUSTDIE):
      call shgi_gr_ab2_store(UA,EA,UB,EB,GR,STO,iop)
      RETURN
    endif

    ! The code below also handles the case UC=0 as if there is no C
    ! (e.g. point charges)
    ASSERT(size(GR,4)==9)

    if(UA>0) call shgi_sym_gr_sto(UA,EA,GR(:,:,:,GAX:GAZ),STO,iop)
    if(UB>0) call shgi_sym_gr_sto(UB,EB,GR(:,:,:,GBX:GBZ),STO,iop)
    if(UC>0)then
      if(.not.is_on(ISLGR) .and. .not.is_on(ISLGRT))then
             call shgi_sym_gr_sto(UC,EC,GR(:,:,:,GCX:GCZ),STO,iop)
      else
           ! NOTE: symmetry type is passed as EC=ISM:
           call shgi_sym_slgr_sto(UC,EC,GR(:,:,:,GCX:GCZ),STO,iop)
      endif
    endif
  end subroutine shgi_gr_abc_store

  subroutine shgi_gr_ab1_store(UA,EA,UB,EB,GR,STO,iop)
    ! properly add kinetic Ga(1:3) = - Gb(1:3) to the total grad
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    type(arrmat4) , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    ! Rotate the C-gradient along the directions
    ! of (totally) symmetric modes:

    ASSERT(size(GR,4)==3)
    ASSERT(iop/=0)
    ! for 2c kinetic energy only G(a) is computed
    ! G(b) = - G(a), G(c) == 0
    call shgi_sym_gr_sto(UA,EA,GR,STO,+iop)
    call shgi_sym_gr_sto(UB,EB,GR,STO,-iop)
  end subroutine shgi_gr_ab1_store

  subroutine shgi_gr_ab2_store(UA,EA,UB,EB,GR,STO,iop)
    ! properly add kinetic Ga(1:3) = - Gb(1:3) to THEIR OWN STORAGE (MUSTDIE)
    use datatype, only: arrmat4 ! PG storage of grs
    use shgi_cntrl, only: NAB, LA, LB
    use shgi_sym, only: shgi_sym_gr_rot
    implicit none
    integer(IK), intent(in) :: UA,EA
    integer(IK), intent(in) :: UB,EB
    real(RK)   , intent(in) :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    type(arrmat4) , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in) :: iop
    ! *** end of interface ***

    integer(IK) :: adr,ngr
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3) ! only ..,:ngr) is defined

    ! so far for all 2c-ints GR(b) = - GR(a):
    ASSERT(size(GR,4)==3)

    ! overlap/kinetic energy have only two grads:
    ! w.r.t. A and B. Storage is allocated only for
    ! totally symmetric modes: NGR(UA) + NGR(UB)
    ! ( even if UA=UB ! )

    ! start at STO(A)%m:
    adr = 1
    call shgi_sym_gr_rot(UA,EA,  GR,SYGR,ngr) ! returns ngr

    call gr_add(adr,ngr,SYGR,STO,iop) ! overwrite STO

    ! point at STO(B)%m:
    adr = adr + ngr

    call shgi_sym_gr_rot(UB,EB, -GR,SYGR,ngr) ! returns ngr

    ! start at STO(B)%m:
    call gr_add(adr,ngr,SYGR,STO,iop) ! overwrite STO
    adr = adr + ngr
    ASSERT(adr==size(STO)+1)
  end subroutine shgi_gr_ab2_store

  subroutine shgi_vd_abc_store(UA,EA,UB,EB,UC,EC,GR,STO,iop)
    !
    ! Rotate the A,B,C derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    ! (vector version)
    !
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    integer(IK), intent(in)    :: UC,EC
    real(RK)   , intent(in)    :: GR(:,:)  ! (9,9)
    real(RK)   , intent(inout) :: STO(:,:) ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    if( UC == 0 .and. EC /= 0 )then
      ! keep special cases to minimum:
      ABORT('must die')
    endif

    ! FIXME: The code below also must handle the case UC=0, EC=0 as if there is no C
    ! (e.g. point charges)
    ASSERT(size(GR,1)==9)
    ASSERT(size(GR,2)==9)

    call shgi_sym_vd_sto( UA,EA, UA,EA, GR( GAX:GAZ, GAX:GAZ ),STO,iop)
    call shgi_sym_vd_sto( UB,EB, UA,EA, GR( GBX:GBZ, GAX:GAZ ),STO,iop)
    call shgi_sym_vd_sto( UA,EA, UB,EB, GR( GAX:GAZ, GBX:GBZ ),STO,iop)
    call shgi_sym_vd_sto( UB,EB, UB,EB, GR( GBX:GBZ, GBX:GBZ ),STO,iop)

    if(UC.gt.0) then
       call shgi_sym_vd_sto( UA,EA, UC,EC, GR( GAX:GAZ, GCX:GCZ ),STO,iop) !3x3
       call shgi_sym_vd_sto( UC,EC, UA,EA, GR( GCX:GCZ, GAX:GAZ ),STO,iop)

       call shgi_sym_vd_sto( UB,EB, UC,EC, GR( GBX:GBZ, GCX:GCZ ),STO,iop) !3x3
       call shgi_sym_vd_sto( UC,EC, UB,EB, GR( GCX:GCZ, GBX:GBZ ),STO,iop)

       call shgi_sym_vd_sto( UC,EC, UC,EC, GR( GCX:GCZ, GCX:GCZ ),STO,iop) !3x3
    endif
  end subroutine shgi_vd_abc_store

  subroutine shgi_sd_abc_store(UA,EA,UB,EB,UC,EC,GR,STO,iop)
    ! Rotate the A,B,C derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    integer(IK), intent(in)    :: UC,EC
    real(RK)   , intent(inout) :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,9,9)
    type(arrmat4) , intent(inout) :: STO(:,:)   ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    if( UC == 0 .and. EC == 1 )then
      ! special case: Ga(3) and Gb(3) = - Ga(3) to be added
      ! to the total grad:
      call shgi_sd_ab1_store(UA,EA,UB,EB,GR,STO,iop)
      RETURN
    endif

    if( UC == 0 .and. EC == 2 )then
      ! special case: Ga(3) and Gb(3) = - Ga(3) to be stored
      ! in its own storage of max. 6 grads (MUSTDIE):
      call shgi_sd_ab2_store(UA,EA,UB,EB,GR,STO,iop)
      RETURN
    endif

    ! FIXME: The code below also must handle the case UC=0, EC=0 as if there is no C
    ! (e.g. point charges)
    ASSERT(size(GR,4)==9)
    ASSERT(size(GR,5)==9)

    call shgi_sym_sd_sto( UA,EA, UA,EA, GR(:,:,:, GAX:GAZ, GAX:GAZ ),STO,iop)
    call shgi_sym_sd_sto( UB,EB, UA,EA, GR(:,:,:, GBX:GBZ, GAX:GAZ ),STO,iop)
    call shgi_sym_sd_sto( UA,EA, UB,EB, GR(:,:,:, GAX:GAZ, GBX:GBZ ),STO,iop)
    call shgi_sym_sd_sto( UB,EB, UB,EB, GR(:,:,:, GBX:GBZ, GBX:GBZ ),STO,iop)

    if(UC.gt.0) then
       call shgi_sym_sd_sto( UA,EA, UC,EC, GR(:,:,:, GAX:GAZ, GCX:GCZ ),STO,iop) !3x3
       call shgi_sym_sd_sto( UC,EC, UA,EA, GR(:,:,:, GCX:GCZ, GAX:GAZ ),STO,iop)

       call shgi_sym_sd_sto( UB,EB, UC,EC, GR(:,:,:, GBX:GBZ, GCX:GCZ ),STO,iop) !3x3
       call shgi_sym_sd_sto( UC,EC, UB,EB, GR(:,:,:, GCX:GCZ, GBX:GBZ ),STO,iop)

       call shgi_sym_sd_sto( UC,EC, UC,EC, GR(:,:,:, GCX:GCZ, GCX:GCZ ),STO,iop) !3x3
    endif
  end subroutine shgi_sd_abc_store

  subroutine shgi_slvd_abc_store(UA,EA,UB,EB,UC,EC,GR,STO,iop)
    !
    ! Rotate the A,B,C solvation derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    ! (vector version)
    !
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    integer(IK), intent(in)    :: UC,EC
    real(RK)   , intent(in)    :: GR(:,:)  ! (9,9)
    real(RK)   , intent(inout) :: STO(:,:) ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    call shgi_sym_slvd_sto( UA,EA, UC,EC, 2, GR( GAX:GAZ, GCX:GCZ ),STO,iop) !3x3
    call shgi_sym_slvd_sto( UC,EC, UA,EA, 1, GR( GCX:GCZ, GAX:GAZ ),STO,iop)

    call shgi_sym_slvd_sto( UB,EB, UC,EC, 2, GR( GBX:GBZ, GCX:GCZ ),STO,iop) !3x3
    call shgi_sym_slvd_sto( UC,EC, UB,EB, 1, GR( GCX:GCZ, GBX:GBZ ),STO,iop)
  end subroutine shgi_slvd_abc_store

  subroutine shgi_slvd_abc_store1(UC1,UC2,EC,GR,STO,iop)
    ! Rotate the A,B,C solvation derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UC1,UC2,EC
    real(RK)   , intent(in)    :: GR(:,:) ! (9,9)
    real(RK)   , intent(inout) :: STO(:,:)   ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    call shgi_sym_slvd_sto( UC1,EC, UC2,EC, 3, GR( GCX:GCZ, GCX:GCZ ),STO,iop) !3x3
  end subroutine shgi_slvd_abc_store1

  subroutine shgi_slvd_abc_store2(UC1,UC2,EC,GR,STO,iop)
    ! Rotate the A,B,C solvation derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    implicit none
    integer(IK), intent(in)    :: UC1,UC2,EC
    real(RK)   , intent(in)    :: GR(:,:)  ! (9,9)
    real(RK)   , intent(inout) :: STO(:,:) ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    ! passes only (GCX:GCZ,GCX:GCZ) of which only vector (GCX,GCX:GCZ) is actually used:
    call shgi_sym_slvd_sto( UC1,EC, UC2,EC, 4, GR( GCX:GCZ, GCX:GCZ ),STO,iop) !3x3
  end subroutine shgi_slvd_abc_store2

  subroutine shgi_sd_ab1_store(UA,EA,UB,EB,GR,STO,iop)
    ! Rotate the A,B,C derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UA,EA
    integer(IK), intent(in)    :: UB,EB
    real(RK)   , intent(inout) :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    type(arrmat4) , intent(inout) :: STO(:,:)   ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    ASSERT(size(GR,4)==3)
    ASSERT(size(GR,5)==3)
    ASSERT(iop/=0)
    ! for 2c kinetic energy only G(a) is computed
    ! G(b) = - G(a), G(c) == 0
    call shgi_sym_sd_sto( UA,EA, UA,EA, GR,STO,+iop)
    call shgi_sym_sd_sto( UB,EB, UA,EA, GR,STO,-iop)
    call shgi_sym_sd_sto( UA,EA, UB,EB, GR,STO,-iop)
    call shgi_sym_sd_sto( UB,EB, UB,EB, GR,STO,+iop)
  end subroutine shgi_sd_ab1_store

  subroutine shgi_sd_ab2_store(UA,EA,UB,EB,GR,STO,iop)
    ! Rotate the A,B derivatives along the directions
    ! of (totally) symmetric modes, and store in PG
    use datatype, only: arrmat4 ! PG storage of grs
    use shgi_sym, only: shgi_sym_sd_rot
    implicit none
    integer(IK), intent(in) :: UA,EA
    integer(IK), intent(in) :: UB,EB
    real(RK)   , intent(in) :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    type(arrmat4) , intent(inout) :: STO(:,:)   ! (num of tot-sym grads^2)
    integer(IK), intent(in) :: iop
    ! *** end of interface ***

    integer(IK) :: ng1,ng2 ! number of tot sym grads
    integer(IK) :: nga,ngb ! number of tot sym grads
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3,3) ! only ..,:ngr) is defined

    ! so far for all 2c-ints GR(b) = - GR(a):
    ASSERT(size(GR,4)==3)
    ASSERT(size(GR,5)==3)

    call shgi_sym_sd_rot(UA,EA,UA,EA,  GR,SYGR,ng1,ng2) ! returns ngr
    ASSERT(ng1==ng2)
    nga = ng1

    call sd_add(    1,    1,ng1,ng2,SYGR,STO,iop)

    call shgi_sym_sd_rot(UB,EB,UB,EB,  GR,SYGR,ng1,ng2) ! returns ngr
    ASSERT(ng1==ng2)
    ngb = ng1

    call sd_add(nga+1,nga+1,ng1,ng2,SYGR,STO,iop)

    call shgi_sym_sd_rot(UB,EB,UA,EA, -GR,SYGR,ng1,ng2) ! returns ngr

    call sd_add(nga+1,    1,ng1,ng2,SYGR,STO,iop)

    call shgi_sym_sd_rot(UA,EA,UB,EB, -GR,SYGR,ng1,ng2) ! returns ngr

    call sd_add(    1,nga+1,ng1,ng2,SYGR,STO,iop)
  end subroutine shgi_sd_ab2_store

  subroutine shgi_sym_ve_sto(UA,EA,GR,STO,iop)
    !
    !  Purpose: Rotate and Store a 3-vector
    !
    use shgi_sym, only: shgi_sym_ve_rot
    implicit none
    integer(IK), intent(in)    :: UA,EA
    real(RK)   , intent(in)    :: GR(:) ! (3)
    real(RK)   , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ngr ! number of tot sym grads
    real(RK)    :: SYGR(3) ! only ..,:ngr) is defined

    ASSERT(UA>0)
    ASSERT(EA>0)

    call shgi_sym_ve_rot(UA,EA,GR,SYGR,ngr) ! returns ngr
    call shgi_sym_ve_add(UA,      SYGR,ngr,STO,iop) ! uses ngr
  end subroutine shgi_sym_ve_sto

  subroutine shgi_sym_gr_sto(UA,EA,GR,STO,iop)
    !
    !  Purpose: Rotate and Store a matrix gradient
    !
    use shgi_cntrl, only: NAB,LA,LB
    use shgi_sym, only: shgi_sym_gr_rot
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UA,EA
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    type(arrmat4) , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ngr ! number of tot sym grads
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3) ! only ..,:ngr) is defined

    ASSERT(UA>0)
    ASSERT(EA>0)

    call shgi_sym_gr_rot(UA,EA,GR,SYGR,ngr) ! returns ngr
    call shgi_sym_gr_add(UA,      SYGR,ngr,STO,iop) ! uses ngr
  end subroutine shgi_sym_gr_sto

  subroutine shgi_sym_slve_sto(UA,ISM,GR,STO,iop)
    !
    !  Purpose: Rotate and Store solv gards (C center)
    !           Vector version.
    !
    use shgi_sym, only: shgi_sym_slve_rot
    implicit none
    integer(IK), intent(in)    :: UA,ISM
    real(RK)   , intent(in)    :: GR(:)  ! (3)
    real(RK)   , intent(inout) :: STO(:) ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ngr ! number of tot sym grads
    real(RK)    :: SYGR(3) ! only ..,:ngr) is defined

    ASSERT(UA>0)

    call shgi_sym_slve_rot(UA,ISM,GR,SYGR,ngr) ! returns ngr
    call shgi_sym_slve_add(UA,       SYGR,ngr,STO,iop) ! uses ngr
  end subroutine shgi_sym_slve_sto

  subroutine shgi_sym_slgr_sto(UA,ISM,GR,STO,iop)
    !
    !  Purpose: Rotate and Store solv gards (C center)
    !           Matrix version.
    !
    use shgi_cntrl, only: NAB,LA,LB
    use shgi_sym, only: shgi_sym_slgr_rot
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UA,ISM
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3)
    type(arrmat4) , intent(inout) :: STO(:)  ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ngr ! number of tot sym grads
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3) ! only ..,:ngr) is defined

    ASSERT(UA>0)

    call shgi_sym_slgr_rot(UA,ISM,GR,SYGR,ngr) ! returns ngr
    call shgi_sym_slgr_add(UA,       SYGR,ngr,STO,iop) ! uses ngr
  end subroutine shgi_sym_slgr_sto

  subroutine shgi_sym_vd_sto(U1,E1,U2,E2,GR,STO,iop)
    !  Purpose: Rotate and Store
    !
    use shgi_sym, only: shgi_sym_vd_rot
    implicit none
    integer(IK), intent(in)    :: U1,E1
    integer(IK), intent(in)    :: U2,E2
    real(RK)   , intent(in)    :: GR(:,:)  ! (3,3)
    real(RK)   , intent(inout) :: STO(:,:) ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ng1,ng2 ! number of tot sym grads
    real(RK)    :: SYGR(3,3) ! only ..,:ngr) is defined

    ASSERT(U1>0)
    ASSERT(E1>0)
    ASSERT(U2>0)
    ASSERT(E2>0)

    call shgi_sym_vd_rot(U1,E1,U2,E2,GR,SYGR,ng1,ng2)     ! returns ngr
    call shgi_sym_vd_add(U1,   U2,      SYGR,ng1,ng2,STO,iop) ! uses ngr
  end subroutine shgi_sym_vd_sto

  subroutine shgi_sym_sd_sto(U1,E1,U2,E2,GR,STO,iop)
    !  Purpose: Rotate and Store
    !
    use shgi_cntrl, only: NAB,LA,LB
    use shgi_sym, only: shgi_sym_sd_rot
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: U1,E1
    integer(IK), intent(in)    :: U2,E2
    real(RK)   , intent(in)    :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    type(arrmat4) , intent(inout) :: STO(:,:) ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ng1,ng2 ! number of tot sym grads
    real(RK)    :: SYGR(NAB,2*LA+1,2*LB+1,3,3) ! only ..,:ngr) is defined

    ASSERT(U1>0)
    ASSERT(E1>0)
    ASSERT(U2>0)
    ASSERT(E2>0)

    call shgi_sym_sd_rot(U1,E1,U2,E2,GR,SYGR,ng1,ng2)     ! returns ngr
    call shgi_sym_sd_add(U1,   U2,      SYGR,ng1,ng2,STO,iop) ! uses ngr
  end subroutine shgi_sym_sd_sto

  subroutine shgi_sym_slvd_sto(U1,E1,U2,E2,IC,GR,STO,iop)
    !
    !  Purpose: Rotate and Store (vector version)
    !
    use shgi_sym, only: shgi_sym_slvd_rot
    implicit none
    integer(IK), intent(in)    :: U1,E1
    integer(IK), intent(in)    :: U2,E2
    integer(IK), intent(in)    :: IC
    real(RK)   , intent(in)    :: GR(:,:) ! (NAB,2*LA+1,2*LB+1,3,3)
    real(RK)   , intent(inout) :: STO(:,:) ! (num of tot-sym grads^2)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ng1,ng2 ! number of tot sym grads
    real(RK)    :: SYGR(3,3) ! only ..,:ngr is defined

    ASSERT(U1>0)
    ASSERT(E1>0)
    ASSERT(U2>0)
    ASSERT(E2>0)

    call shgi_sym_slvd_rot(U1,E1,U2,E2,IC,GR,SYGR,ng1,ng2)     ! returns ngr
    call shgi_sym_vd_add(U1,   U2,      SYGR,ng1,ng2,STO,iop) ! uses ngr
  end subroutine shgi_sym_slvd_sto

  subroutine shgi_sym_X_add(Xc,UC,GR,ngr,STO,iop)
    !  Purpose: store the X grads in PG structure
    use datatype, only: arrmat4 ! PG storage of grs
    use pointcharge_module, only: unique_pc_index
    use point_dqo_module, only: dipoles_grad_index
    use point_dqo_module, only: qpoles_grad_index
    use point_dqo_module, only: opoles_grad_index
    use point_dqo_module, only: rep_grad_index
    use induced_dipoles_module, only: ind_dip_grad_index
    implicit none
    integer(IK), intent(in)    :: Xc,UC, NGR
    real(RK)   , intent(inout) :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3)
    type(arrmat4) , intent(inout) :: STO(:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***
    integer(IK) :: ADR

    if(Xc==PC)  ADR = unique_pc_index(UC)
    if(Xc==PD)  ADR = dipoles_grad_index(UC)
    if(Xc==PQ)  ADR = qpoles_grad_index(UC)
    if(Xc==PO)  ADR = opoles_grad_index(UC)
    if(Xc==IPD) ADR = ind_dip_grad_index(UC)
    if(Xc==RC)  ADR = rep_grad_index(UC)

    call X_add(ADR,NGR,GR,STO,iop)
  end subroutine shgi_sym_X_add

  subroutine shgi_sym_fl_add(UC,EF,ngr,STO,iop)
    !  Purpose: store the electrostatic field in PG structure
    use elec_static_field_module, only: surf_points_grad_index
    implicit none
    integer(IK), intent(in)    :: UC, NGR
    real(RK)   , intent(inout) :: EF(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3)
    real(RK)   , intent(inout) :: STO(:,:,:,:,:)
    !                                     ^ - num of tot-sym grads
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***
    integer(IK) :: ADR

    ADR = surf_points_grad_index(UC)

    call fl_add(ADR,NGR,EF,STO,iop)
  end subroutine shgi_sym_fl_add

  subroutine shgi_sym_ve_add(UA,GR,NGR,STO,iop)
    !
    !  Purpose: store the gradient vector in PG array
    !
    use gradient_data_module, only: gradient_index
#ifdef WITH_EFP
    use pointcharge_module, only: unique_pc_index
    use shgi_cntrl, only: do_efp_grad
#endif
    implicit none
    integer(IK), intent(in)    :: UA, NGR
    real(RK)   , intent(in)    :: GR(:) ! (:NGR) used, NGR<=3
    real(RK)   , intent(inout) :: STO(:)! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ADR,a,igr

    ! FIXME: split gradients/moving atoms need more work:
    ASSERT(UA>0)
#ifdef WITH_EFP
    if(do_efp_grad) then
       ASSERT(UA<=size(unique_pc_index))

       ADR = unique_pc_index(UA)
    else
       ASSERT(UA<=size(gradient_index))

       ADR = gradient_index(UA)
    end if
#else
    ASSERT(UA<=size(gradient_index))

    ADR = gradient_index(UA)
#endif

    do igr=1,NGR
       ! offset into array of totally-symmetric gradients:
       a = ADR + igr - 1
       ASSERT(a<=size(STO))

       ! overwrite, add, or subtract depending on iop:
       select case( iop )
       case (  0 )
         ! overwrite:
         STO(a) =   GR(igr)
       case ( +1 )
         ! add:
         STO(a) = + GR(igr) + STO(a)
       case ( -1 )
         ! subtract:
         STO(a) = - GR(igr) + STO(a)
       case default
         ABORT('no such iop')
       end select
    enddo
  end subroutine shgi_sym_ve_add

  subroutine shgi_sym_gr_add(UA,GR,NGR,STO,iop)
    !
    !  Purpose: store the gradient matrix in PG structure
    !
    !  Note: consider first contracting the matrix with
    !        density matrix in property runs
    !
    use gradient_data_module, only: gradient_index
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: UA, NGR
    real(RK)   , intent(in)       :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<=3)
    type(arrmat4) , intent(inout) :: STO(:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ADR

    ! FIXME: split gradients/moving atoms need more work:
    ASSERT(UA<=size(gradient_index))

    ADR = gradient_index(UA)
    ! FIXME: V=-V different for energy and grads:
    call gr_add(ADR,NGR,GR,STO,iop)
  end subroutine shgi_sym_gr_add

  subroutine shgi_sym_slve_add(UA,GR,NGR,STO,iop)
    !  Purpose: store the gradient in PG structure
    !
    use gradient_data_module, only: gradient_index
    use unique_atom_module, only: n_unique_atoms
    use elec_static_field_module, only: surf_points_grad_index
    implicit none
    integer(IK), intent(in)    :: UA, NGR
    real(RK)   , intent(in)    :: GR(:)  ! (3) only (:NGR) used
    real(RK)   , intent(inout) :: STO(:) ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ADR,a,igr

    if(UA <= n_unique_atoms) then
       ADR = gradient_index(UA)
    else
       ADR = surf_points_grad_index(UA-n_unique_atoms)
    end if

    do igr=1,NGR
       ! offset into array of totally-symmetric gradients:
       a = ADR + igr - 1
       ASSERT(a<=size(STO))

       ! overwrite, add, or subtract depending on iop:
       select case( iop )
       case (  0 )
         ! overwrite:
         STO(a) =   GR(igr)
       case ( +1 )
         ! add:
         STO(a) = + GR(igr) + STO(a)
       case ( -1 )
         ! subtract:
         STO(a) = - GR(igr) + STO(a)
       case default
         ABORT('no such iop')
       end select
    enddo
  end subroutine shgi_sym_slve_add

  subroutine shgi_sym_slgr_add(UA,GR,NGR,STO,iop)
    !  Purpose: store the gradient in PG structure
    !
    use gradient_data_module, only: gradient_index
    use datatype, only: arrmat4 ! PG storage of grs
    use unique_atom_module, only: n_unique_atoms
    use elec_static_field_module, only: surf_points_grad_index
    implicit none
    integer(IK), intent(in)    :: UA, NGR
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3)
    type(arrmat4) , intent(inout) :: STO(:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ADR

    if(UA <= n_unique_atoms) then
       ADR = gradient_index(UA)
    else
       ADR = surf_points_grad_index(UA-n_unique_atoms)
    end if
    ! FIXME: V=-V different for energy and grads:
    call gr_add(ADR,NGR,GR,STO,iop)
  end subroutine shgi_sym_slgr_add

  subroutine X_add(ADR,NGR,GR,STO,iop)
    !  Purpose: store the X grads in PG structure
    !
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: ADR      ! start at ADR
    integer(IK), intent(in)    :: NGR      ! num of grads
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3)
    type(arrmat4) , intent(inout) :: STO(:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: igr,a

    do igr=1,NGR
       a = ADR + igr - 1
       ASSERT(a<=size(STO))
       ! overwrite, add, or subtract depending on iop:
       call upack2c(GR(:,:,:,igr),STO(a)%m,iop)
    end do
  end subroutine X_add

  subroutine fl_add(ADR,NGR,GR,STO,iop)
    !  Purpose: store the electrostatic field in PG structure
    !
    implicit none
    integer(IK), intent(in)    :: ADR      ! start at ADR
    integer(IK), intent(in)    :: NGR      ! num of grads
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3)
    real(RK)   , intent(inout) :: STO(:,:,:,:,:)
    !                                     ^ - num of tot-sym grads
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: igr,a

    do igr=1,NGR
       a = ADR + igr - 1
       ASSERT(a<=size(STO,3))
       ! overwrite, add, or subtract depending on iop:
       call upack2c(GR(:,:,:,igr),STO(:,:,a,:,:),iop)
    end do
  end subroutine fl_add

  subroutine gr_add(ADR,NGR,GR,STO,iop)
    !  Purpose: store the gradient in PG structure
    !
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: ADR      ! start at ADR
    integer(IK), intent(in)    :: NGR      ! num of grads
    real(RK)   , intent(in)    :: GR(:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3)
    type(arrmat4) , intent(inout) :: STO(:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: igr,a

    do igr=1,NGR
       a = ADR + igr - 1
       ASSERT(a<=size(STO))
       ! overwrite, add, or subtract depending on iop:
       call upack2c(GR(:,:,:,igr),STO(a)%m,iop)
    end do
  end subroutine gr_add

  subroutine shgi_sym_vd_add(U1,U2,GR,NG1,NG2,STO,iop)
    !
    !  Purpose: store the gradient in PG structure
    !           (vector version)
    !
    use gradient_data_module, only: gradient_index
    implicit none
    integer(IK), intent(in)    :: U1, NG1
    integer(IK), intent(in)    :: U2, NG2
    real(RK)   , intent(in)    :: GR(:,:)  ! (<=3,<=3)
    real(RK)   , intent(inout) :: STO(:,:) ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: ADR1
    integer(IK) :: ADR2
    integer(IK) :: igr1,a1
    integer(IK) :: igr2,a2

    ! FIXME: split gradients/moving atoms need more work:
    ASSERT(U1<=size(gradient_index))
    ASSERT(U2<=size(gradient_index))

    ADR1 = gradient_index(U1)
    ADR2 = gradient_index(U2)

!   print *,'shape(STO)=',shape(STO)
!   print *,'ADR1,NG1=',ADR1,NG1
!   print *,'ADR2,NG2=',ADR2,NG2

    do igr2=1,NG2
       a2 = ADR2 + igr2 - 1
    do igr1=1,NG1
       a1 = ADR1 + igr1 - 1
       ASSERT(a1<=size(STO,1))
       ASSERT(a2<=size(STO,2))
       ! overwrite, add, or subtract depending on iop::
       select case( iop )
       case (  0 )
         ! overwrite:
         STO(a1,a2) =   GR(igr1,igr2)
       case ( +1 )
         ! add:
         STO(a1,a2) = + GR(igr1,igr2) + STO(a1,a2)
       case ( -1 )
         ! subtract:
         STO(a1,a2) = - GR(igr1,igr2) + STO(a1,a2)
       case default
         ABORT('no such iop')
       end select
    enddo
    enddo
  end subroutine shgi_sym_vd_add

  subroutine shgi_sym_sd_add(U1,U2,GR,NG1,NG2,STO,iop)
    !
    !  Purpose: store the gradient in PG structure
    !           (matrix version)
    !
    use gradient_data_module, only: gradient_index
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: U1, NG1
    integer(IK), intent(in)    :: U2, NG2
    real(RK)   , intent(in)    :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3,<3)
    type(arrmat4) , intent(inout) :: STO(:,:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: adr1
    integer(IK) :: adr2

    ! FIXME: split gradients/moving atoms need more work:
    ASSERT(U1<=size(gradient_index))
    ASSERT(U2<=size(gradient_index))

    adr1 = gradient_index(U1)
    adr2 = gradient_index(U2)

    call sd_add(ADR1,ADR2,NG1,NG2,GR,STO,iop)
  end subroutine shgi_sym_sd_add

  subroutine sd_add(ADR1,ADR2,NG1,NG2,GR,STO,iop)
    !  Purpose: store the gradient in PG structure
    !
    use datatype, only: arrmat4 ! PG storage of grs
    implicit none
    integer(IK), intent(in)    :: ADR1, NG1
    integer(IK), intent(in)    :: ADR2, NG2
    real(RK)   , intent(in)    :: GR(:,:,:,:,:) ! (NAB,2*LA+1,2*LB+1,<3,<3)
    type(arrmat4) , intent(inout) :: STO(:,:)   ! (num of tot-sym grads)
    integer(IK), intent(in)    :: iop
    ! *** end of interface ***

    integer(IK) :: igr1,a1
    integer(IK) :: igr2,a2

!   print *,'shape(STO)=',shape(STO)
!   print *,'ADR1,NG1=',ADR1,NG1
!   print *,'ADR2,NG2=',ADR2,NG2

    do igr2=1,NG2
       a2 = ADR2 + igr2 - 1
    do igr1=1,NG1
       a1 = ADR1 + igr1 - 1
       ASSERT(a1<=size(STO,1))
       ASSERT(a2<=size(STO,2))
       ! overwrite, add, or subtract depending on iop::
       call upack2c(GR(:,:,:,igr1,igr2),STO(a1,a2)%m,iop)
    enddo
    enddo
  end subroutine sd_add

  subroutine upack2c(pckd,upckd,iop)
    implicit none
    real(RK), intent(in)  :: pckd(:,:,:)
    real(RK), intent(out) :: upckd(:,:,:,:)
    integer(IK), intent(in), optional :: iop
    ! *** end of interface ***

    integer(IK) :: sel
    integer(IK) :: ma,mb

    sel = 0
    if(present(iop)) sel = iop

    do mb=1,2*LB+1
       do ma=1,2*LA+1
          call upackSS(pckd(:,ma,mb),upckd(:,:,ma,mb),sel)
       enddo
    enddo
  end subroutine upack2c

  subroutine upack3c(icc,pckd,upckd)
    implicit none
    integer(IK), intent(inout) :: icc
    real(RK)   , intent(in)    :: pckd(:,:,:,:)    ! (  NAB  ,NEC,2*LA+1,2*LB+1)
    real(RK)   , intent(out)   :: upckd(:,:,:,:,:) ! (NEA,NEB,NEC,2*LA+1,2*LB+1)
    ! *** end of interface ***

    integer(IK) :: ma,mb,ic,ff

    ASSERT(size(pckd,2)<=size(upckd,3))
    ASSERT(size(pckd,3)==size(upckd,4))
    ASSERT(size(pckd,4)==size(upckd,5))

    do mb=1,2*LB+1
       do ma=1,2*LA+1
          ff = icc
          do ic=1,size(pckd,2)
             ff = ff + 1
             call upackSS(pckd(:,ic,ma,mb),upckd(:,:,ff,ma,mb),0)
             !upckd(:,:,ff,ma,mb) = UNPACK(pckd(:,ic,ma,mb),CUTOFF,ZERO)
          enddo
       enddo
    enddo
    icc = ff
  end subroutine upack3c

  subroutine upackSS(pckd,upckd,iop)
    use shgi_common, only: CUTOFF
    implicit none
    real(RK), intent(in)    :: pckd(:)    ! (NAB)
    real(RK), intent(inout) :: upckd(:,:) ! (NEA,NEB)
    integer(IK), intent(in) :: iop
    ! *** end of interface ***

    select case( iop )

    case (  0 )
      ! overwrite:
      upckd(:,:) =   UNPACK(pckd(:),CUTOFF,ZERO)

    case ( +1 )
      ! add:
      upckd(:,:) = + UNPACK(pckd(:),CUTOFF,ZERO) + upckd(:,:)

    case ( -1 )
      ! subtract:
      upckd(:,:) = - UNPACK(pckd(:),CUTOFF,ZERO) + upckd(:,:)

    case default
    ABORT('no such iop')
    end select
  end subroutine upackSS

  subroutine shgi_timing(res)
    use error_module, only: MyID
    real(RK), intent(in), optional :: res!olution
    ! *** end of interface ***

    real(RK), parameter :: interval = 60.0_rk
    real(RK), save :: told = zero
    real(RK)       :: t,cmp

    t = FPP_TIMER_VALUE(tot)

    cmp = interval
    if(present(res)) cmp = res

    if( t-told <= cmp ) return ! not too frequently
    told = t

    print *,MyID//'SHGI: =============='
    print *,MyID//'SHGI: total    time:', FPP_TIMER_VALUE(tot)
    if( FPP_TIMER_VALUE(totG) > 0 )then
    print *,MyID//'SHGI: |-0 ints time:', FPP_TIMER_VALUE(totI)
    print *,MyID//'SHGI: |-1 grad time:', FPP_TIMER_VALUE(totG)
    if( FPP_TIMER_VALUE(totD) > 0 ) &
    print *,MyID//'SHGI: `-2 derv time:', FPP_TIMER_VALUE(totD)
    endif
    print *,MyID//'SHGI: |-2cent  time:', FPP_TIMER_VALUE(t2c)
    print *,MyID//'SHGI: |-gamma  time:', FPP_TIMER_VALUE(tgam)
    print *,MyID//'SHGI: |-radang time:', FPP_TIMER_VALUE(tra)
    print *,MyID//'SHGI: |-set_c  time:', FPP_TIMER_VALUE(tsc)
    print *,MyID//'SHGI: | |-D3A  time:', FPP_TIMER_VALUE(td3a)
    print *,MyID//'SHGI: | `-YLS  time:', FPP_TIMER_VALUE(tyls)
    print *,MyID//'SHGI: |---PR2  time:', FPP_TIMER_VALUE(tpr2)
    if( FPP_TIMER_VALUE(adkh) > 0 )then![[
    print *,MyID//'SHGI: `-adkh   time:', FPP_TIMER_VALUE(adkh)
    endif!]]
    if(is_on(IPSEU+IPSGR+IPSSD))then
    print *,MyID//'SHGI: P-Q1     time:', FPP_TIMER_VALUE(tpq1)
    print *,MyID//'SHGI: | `-bs1  time:', FPP_TIMER_VALUE(tbs1)
    print *,MyID//'SHGI: |-Q2     time:', FPP_TIMER_VALUE(tpq2)
    print *,MyID//'SHGI: | `-bs2  time:', FPP_TIMER_VALUE(tbs2)
    print *,MyID//'SHGI: |-set_c  time:', FPP_TIMER_VALUE(tpsc)
    print *,MyID//'SHGI: `-radang time:', FPP_TIMER_VALUE(tra2)
    endif
    if( FPP_TIMER_VALUE(pcs) > 0 )then
    print *,MyID//'SHGI: PC-tot   time:', FPP_TIMER_VALUE(pcs)
    print *,MyID//'SHGI:   |-rad  time:', FPP_TIMER_VALUE(tpcr)
    print *,MyID//'SHGI:   |-D3F  time:', FPP_TIMER_VALUE(td3f)
    if( FPP_TIMER_VALUE(epot) > 0 ) &
    print *,MyID//'SHGI:   |-slv0 time:', FPP_TIMER_VALUE(epot)
    print *,MyID//'SHGI:   |-slv1 time:', FPP_TIMER_VALUE(slgr)
    print *,MyID//'SHGI:   `-slv2 time:', FPP_TIMER_VALUE(slsd)
    endif
    print *,MyID//'SHGI: int screening:',shgi_stat_screened,'/',shgi_stat_num_ints &
                               ,'=',real(shgi_stat_screened) /  shgi_stat_num_ints
  end subroutine shgi_timing

  !--------------- End of module -------------------------------------
end module shgi_utils

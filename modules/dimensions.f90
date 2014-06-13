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
module dimensions
  !---------------------------------------------------------------
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

# include "def.h"
  use type_module,only:&
       & IK=>i4_kind,&
       & RK=>r8_kind ! type specification parameters
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !=======================================================
  !
  ! A STRUCT TO HOLD DIMENSIONS >>>
  !

  type,public :: SubSpaceDim
     !
     integer(IK)         :: lmax, n_irr
     !
     integer(IK), allocatable :: LM(:,:)
     ! LM(0:lmax,n_irr) ~ sum(JM,2) (symbolic relation)
     ! i.e. sum of dimensions of j = L + 1/2, and j = L - 1/2
     ! (if the last exists, of course)
     integer(IK), allocatable :: n_rad(:) ! n_rad(lmax)
     ! number of radial functions (exponents,contractions)
     ! historically LM(:,L) _IS_ already scaled by n_rad(L)
     ! (by default n_rad == 1)
  end type SubSpaceDim

  !------------ Declaration of constants and variables ------------
  integer(IK), parameter, public :: UNCONTRACTED = 0
  integer(IK), parameter, public :: CONTRACTED   = 1
  integer(IK), parameter, public :: ANGULAR      = -1

  !
  ! FIXME: these are not protected, because they are set from symm_adapt_module:
  !
  ! Y_LM splitting:
  type(SubSpaceDim), public, allocatable :: SymDim(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable :: SymDimSpor(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable :: SymDimSporL(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable :: SymDimSporS(:) ! (n_ua)

  ! Y_LM (*) RADIAL CONTRACTED splitting:
  type(SubSpaceDim), public, allocatable, protected :: BasDim(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable, protected :: BasDimSpor(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable, protected :: BasDimSporL(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable, protected :: BasDimSporS(:) ! (n_ua)

  ! Y_LM (*) RADIAL UNCONTRACTED splitting:
  type(SubSpaceDim), public, allocatable, protected :: UBasDim(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable, protected :: UBasDimSpor(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable, protected :: UBasDimSporL(:) ! (n_ua)
  type(SubSpaceDim), public, allocatable, protected :: UBasDimSporS(:) ! (n_ua)

  ! dimensions of UA subspaces:
  ! UA_space = SUM[L=0,L(ua)]{ UA_L_space }
  integer(IK), allocatable, public, protected :: UABasDimSpor(:, :) ! ...(n_irr, n_ua)

  ! dimensions of reduced basis:
  integer(IK), allocatable, public, protected :: IrrBasDim(:)
  integer(IK), allocatable, public, protected :: IrrBasDimSpor(:)
  integer(IK), allocatable, public, protected :: IrrBasDimSporL(:)
  integer(IK), allocatable, public, protected :: IrrBasDimSporS(:)
  ! IrrBasDimXXX = SUM[ua=1,N_ua;L=0,L(ua)]{ BasDimXXX(ua,L) }
  ! i.e. full dimensions of irreps

  ! dimensions of reduced uncontracted basis:
  integer(IK), allocatable, public, protected :: IrrUBasDim(:)
  integer(IK), allocatable, public, protected :: IrrUBasDimSpor(:)
  integer(IK), allocatable, public, protected :: IrrUBasDimSporL(:)
  integer(IK), allocatable, public, protected :: IrrUBasDimSporS(:)
  ! IrrUBasDimXXX = SUM[ua=1,N_ua;L=0,L(ua)]{ UBasDimXXX(ua,L) }
  ! i.e. full dimensions of uncontruncted irreps

  integer(IK), public, protected :: number_of_irreps = -1

  !
  ! FIXME: "shell" is probably a better prefix for the public names
  !         than the cryptic "uaL". Rename them some day ...
  !
  ! uaL is a metaindex (ua, L) == (unique atom ua, angular momentum L)
  ! which is given by a metacode:
  !
  ! uaL = 0
  ! do ua=1,n_ua
  !    do L=1,Lmax(ua)
  !       uaL = uaL + 1
  !       Whatever(uaL) = AnythingElse(ua, L)
  !    enddo
  ! enddo
  !
  integer(IK),              public, protected :: uaL_max = -1
  integer(IK), allocatable, public, protected :: uaL_nrad(:) ! (uaL_max)
  integer(IK), allocatable, public, protected :: uaL_vec_mult(:, :) ! (uaL_max, n_vec_irrep)
  integer(IK), allocatable, public, protected :: uaL_vec_dims(:, :) ! (uaL_max, n_vec_irrep)
  integer(IK), allocatable, public, protected :: uaL_proj_dims(:, :) ! (uaL_max, n_proj_irrep)

  !------------ Interface statements ------------------------------

  interface dimens
     module procedure dimens_sym     ! (SYM)
     module procedure dimens_sym_u_l ! (SYM,U,L)
  end interface

  interface dimoff
     module procedure dimoff_sym_u_l ! (SYM,U,L)
  end interface

  interface alloc
     module procedure alloc_SubSpaceDim
  end interface

  interface dealloc
     module procedure dealloc_SubSpaceDim
     module procedure dealloc_arr_SubSpaceDim
  end interface

  !------------ public functions and subroutines ------------------

  public :: dimens
  public :: dimoff
  public :: dimensions_init
  public :: dimensions_free
  public :: alloc !,dealloc

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

! /* integer parameter bug on SGI */
#ifdef _SGI
# define ENUM_T character
# define VInitLevelNone     '0'
# define VInitLevelAlloc    '1'
# define VInitLevelComplete '2'
# define VInitLevelDone     '3'
#else
# define ENUM_T integer(IK)
# define VInitLevelNone     -1
# define VInitLevelAlloc    10
# define VInitLevelComplete 20
# define VInitLevelDone    100
#endif

  ENUM_T, parameter, PUBLIC ::&
       & InitLevelNone     =  VInitLevelNone,&
       & InitLevelAlloc    =  VInitLevelAlloc,&
       & InitLevelComplete =  VInitLevelComplete,&
       & InitLevelDone     =  VInitLevelDone

  ENUM_T, private :: initialized = InitLevelNone

  logical :: SpinOrbit
  logical :: Small

  integer(IK), private :: REFCOUNT = 0

  !----------------------------------------------------------------
contains
  !------------ Subroutines ---------------------------------------

  !*************************************************************
  function dimens_sym(SYM, c) result(N)
    ! returns dimension of SYM irrep
    ! FIXME: contraction, SO, ...
    use unique_atom_module, only: uas=>unique_atoms
    implicit none
    integer(IK), intent(in) :: SYM
    integer(IK), intent(in) :: c ! 1 for contracted, 0 for uncontracted
    integer(IK)             :: N ! result
    ! *** end of interface ***

    integer(IK) :: U, L

    N = 0
    do U = 1, size(uas)
       do L = 0, uas(U)%lmax_ob
          N = N + dimens(SYM, U, L, c)
       enddo
    enddo
  end function dimens_sym
  !*************************************************************

  !*************************************************************
  function dimens_sym_u_l(SYM, U, L, c) result(N)
    use unique_atom_module, only: uas=>unique_atoms
    implicit none
    integer(IK), intent(in) :: SYM,U,L
    integer(IK), intent(in) :: c ! 1 for contracted, 0 for uncontracted
    integer(IK)             :: N ! result
    ! *** end of interface ***

    integer(IK) :: NIF, NEXP, NCNT, NUNC

    ASSERT(U>0)
    ASSERT(U<=size(uas))
    ASSERT(SYM>0)
    ! NIF is used as temp for NIRR and LMAX:
    NIF = size(uas(U)%symadapt_partner,1)
    ASSERT(SYM<=NIF)
    ASSERT(L>=0)
    NIF = ubound(uas(U)%symadapt_partner,2)
    ASSERT(L<=NIF)
    NIF = ubound(uas(U)%l_ob,1)
    ASSERT(L<=NIF)

    ! FIXME: options for contracted, SO, ...
    NIF  = uas(U)%symadapt_partner(SYM,L)%n_independent_fcts
    NEXP = uas(U)%l_ob(L)%n_exponents
    NCNT = uas(U)%l_ob(L)%n_contracted_fcts
    NUNC = uas(U)%l_ob(L)%n_uncontracted_fcts

    select case ( c )
    case (UNCONTRACTED)
      ! uncontracted:
      N = NIF * NEXP
    case (CONTRACTED)
      ! contracted:
      N = NIF * ( NCNT + NUNC )
    case (ANGULAR)
      ! only one radial function
      N = NIF
    case default
      N = -1
      print*,'no such c=',c
      ABORT('no such case, see tty')
    end select
  end function dimens_sym_u_l
  !*************************************************************

  !*************************************************************
  function dimoff_sym_u_l(SYM,U,L,c) result(N)
    ! returns the offset of (SYM,U,L) --
    ! a partial sum of all subspace dimensions preceeding
    ! the triple (SYM,U,L)
    use unique_atom_module, only: uas=>unique_atoms
    implicit none
    integer(IK), intent(in) :: SYM,U,L
    integer(IK), intent(in) :: c ! 1 for contracted, 0 for uncontracted
    integer(IK)             :: N ! result
    ! *** end of interface ***

    integer(IK) :: uu,ll

    N = 0
    uaL: do uu=1,size(uas)
       do ll=0,uas(uu)%lmax_ob
          if( uu==U .and. ll==L ) exit uaL
          N = N + dimens(SYM,uu,ll,c)
       enddo
    enddo uaL
  end function dimoff_sym_u_l
  !*************************************************************

  subroutine dimensions_init(InitLevel)
    use error_module, only: MyID
    use unique_atom_module, only:&
         & n_ua=>n_unique_atoms
    use spin_orbit_module, only: &
         & is_on,&
         & op_FitTrafo,&
         & op_SpinOrbit
    implicit none
    ENUM_T, intent(in)          :: InitLevel
    ! *** end of interface ***

    integer :: memstat,a,L

    SpinOrbit = is_on(op_SpinOrbit)
    Small = is_on(op_FitTrafo)

    select case(InitLevel)

    case(InitLevelAlloc)

       ASSERT(initialized==InitLevelNone)

       allocate(SymDim(n_ua), STAT=memstat)
       ASSERT(memstat==0)

       if(SpinOrbit)then
          allocate(&
               & SymDimSpor(n_ua),&
               & SymDimSporL(n_ua),&
               & SymDimSporS(n_ua),&
               & STAT=memstat)
          ASSERT(memstat==0)
       endif

       allocate(BasDim(n_ua), STAT=memstat)
       ASSERT(memstat==0)

       if(SpinOrbit)then
          allocate(&
               & BasDimSpor(n_ua),&
               & BasDimSporL(n_ua),&
               & BasDimSporS(n_ua),&
               & STAT=memstat)
          ASSERT(memstat==0)
       endif

       allocate(&
            & UBasDim(n_ua),&
            & STAT=memstat)
       ASSERT(memstat==0)

       if(SpinOrbit)then
          allocate(&
               & UBasDimSpor(n_ua),&
               & UBasDimSporL(n_ua),&
               & UBasDimSporS(n_ua),&
               & STAT=memstat)
          ASSERT(memstat==0)
       endif

       initialized = InitLevelAlloc

    case (InitLevelComplete)
       ASSERT(initialized==InitLevelAlloc)

       DPRINT 'dimesions: call CalcBasDim(SymDim,UBasDim,IrrUBasDim ,CNT=.false.)'
       call CalcBasDim(SymDim, UBasDim, IrrUBasDim, CNT=.false.)
       call CalcBasDim(SymDim, BasDim, IrrBasDim,  CNT=.true.)


       if(SpinOrbit)then
          call CalcBasDim(SymDimSpor ,UBasDimSpor ,IrrUBasDimSpor,CNT=.false.)
          call CalcBasDim(SymDimSpor , BasDimSpor , IrrBasDimSpor,CNT=.true.)

          call CalcBasDim(SymDimSporL,UBasDimSporL,IrrUBasDimSporL,CNT=.false.) ! dublicate
          call CalcBasDim(SymDimSporL, BasDimSporL, IrrBasDimSporL,CNT=.true.)  ! with above

          if(is_on(op_FitTrafo))then
             call CalcBasDim(SymDimSporS,UBasDimSporS,IrrUBasDimSporS,CNT=.false.)
             call CalcBasDim(SymDimSporS, BasDimSporS, IrrBasDimSporS,CNT=.true.)
          endif
       endif

       if(SpinOrbit)then
          number_of_irreps = size(IrrBasDimSpor)
       else
          number_of_irreps = size(IrrBasDim)
       endif

       DPRINT MyID,'dim/Init: call CalcUALDims()'
       call CalcUALDims()
       DPRINT MyID,'dim/Init: done CalcUALDims()'


       if(SpinOrbit)then
          DPRINT MyID,'dim/Init: calculating UABasDimSpor'
          allocate(UABasDimSpor(number_of_irreps, size(BasDimSpor)))

          do a = 1, size(BasDimSpor)
             UABasDimSpor(:, a) = 0
             do L = 0, BasDimSpor(a)%lmax
                UABasDimSpor(:, a) = UABasDimSpor(:, a) + BasDimSpor(a)%LM(L, :)
             enddo
             DPRINT MyID,'dim/Init: dims of subspace UA=',a,' : ',UABasDimSpor(:, a)
          enddo
          ! ASSERT SymDimSpor == SymDimSporL, identically
       endif

#ifdef FPP_DEBUG
       call dimensions_print()
#endif

       REFCOUNT = REFCOUNT + 1
       ASSERT(REFCOUNT==1)

    case default
       print *,MyID,'InitLevel=',InitLevel
       ABORT("no action for this key")
    end select
  end subroutine dimensions_init
 
  subroutine dimensions_free()
    use spin_orbit_module, only: is_on,op_FitTrafo
    implicit none
    ! *** end of interface ***

    integer :: memstat

    if(REFCOUNT==0)then
       WARN('dimensions_free: not used, ignoring')
       return
    endif

    call dealloc(SymDim)

    if(SpinOrbit)then
       call dealloc(SymDimSpor)
       call dealloc(SymDimSporL)
       if(Small) call dealloc(SymDimSporS)
    endif

    call dealloc(UBasDim)

    if(SpinOrbit)then
       call dealloc(UBasDimSpor)
       call dealloc(UBasDimSporL)
       if(Small) call dealloc(UBasDimSporS)
    endif

    call dealloc(BasDim)

    if(SpinOrbit)then
       call dealloc(BasDimSpor)
       call dealloc(BasDimSporL)
       if(Small) call dealloc(BasDimSporS)
    endif

    deallocate(IrrUBasDim)

    if(SpinOrbit)then
       deallocate(IrrUBasDimSpor)
       deallocate(IrrUBasDimSporL)
       if(Small) deallocate(IrrUBasDimSporS)
    endif

    deallocate(IrrBasDim)

    if(SpinOrbit)then
       deallocate(IrrBasDimSpor)
       deallocate(IrrBasDimSporL)
       if(Small) deallocate(IrrBasDimSporS)
    endif

    deallocate(SymDim, BasDim, UBasDim, STAT=memstat)
    ASSERT(memstat==0)

    if(SpinOrbit)then
       deallocate(&
            & SymDimSpor, SymDimSporL, SymDimSporS,&
            & BasDimSpor, BasDimSporL, BasDimSporS,&
            & UBasDimSpor,UBasDimSporL,UBasDimSporS,&
            & STAT=memstat)
       ASSERT(memstat==0)
    endif

    if ( allocated(UABasDimSpor) ) then
       deallocate(UABasDimSpor, STAT=memstat)
       ASSERT(memstat==0)
    endif

    if ( allocated(uaL_nrad) ) then
       deallocate(uaL_nrad, STAT=memstat)
       ASSERT(memstat==0)
    endif

    if ( allocated(uaL_vec_mult) ) then
       deallocate(uaL_vec_mult, STAT=memstat)
       ASSERT(memstat==0)
    endif

    if ( allocated(uaL_vec_dims) ) then
       deallocate(uaL_vec_dims, STAT=memstat)
       ASSERT(memstat==0)
    endif

    if ( allocated(uaL_proj_dims) ) then
        deallocate(uaL_proj_dims, STAT=memstat)
        ASSERT(memstat==0)
    endif

    uaL_max = -1

    initialized = InitLevelNone
    REFCOUNT = REFCOUNT - 1
  end subroutine dimensions_free

#ifdef FPP_DEBUG
  subroutine dimensions_print()
    use error_module, only: MyID
    use spin_orbit_module, only: is_on,op_FitTrafo
    implicit none
    ! *** end of interface ***

    integer(IK) :: a

    print *,MyID,'dim/Init: SymDim >>>'
    call print_arr_SubSpaceDim(SymDim)
    if(SpinOrbit)then
       print *,MyID,'dim/Init: SymDimSpor >>>'
       call print_arr_SubSpaceDim(SymDimSpor)
       print *,MyID,'dim/Init: BasDimSpor >>>'
       call print_arr_SubSpaceDim(BasDimSpor)
       print *,MyID,'dim/Init: SymDimSporL >>>'
       call print_arr_SubSpaceDim(SymDimSporL)
       if(is_on(op_FitTrafo))then
          print *,MyID,'dim/Init: SymDimSporS >>>'
          call print_arr_SubSpaceDim(SymDimSporS)
       endif
    endif

    print *,MyID,'dim/Init: Vector:'
    do a = 1, size(IrrBasDim)
       print *,MyID,'dim/Init: irr=',a,'IrrDim =',IrrBasDim(a)
    enddo

    if(SpinOrbit)then
       print *,MyID,'dim/Init: Projective: irr,IrrDim,IrrDimL,IrrDimS'
       do a = 1, size(IrrBasDimSpor)
          if(is_on(op_FitTrafo))then
             print *,MyID,'dim/Init: ',a, &
                  IrrBasDimSpor(a), &
                  IrrBasDimSporL(a), &
                  IrrBasDimSporS(a)
          else
             print *,MyID,'dim/Init: ',a, &
                  IrrBasDimSpor(a), &
                  IrrBasDimSporL(a)
          endif
       enddo
    endif
  end subroutine dimensions_print
#endif

  subroutine CalcUALDims()
    implicit none
    ! *** end of interface ***

    integer(IK) :: n_ua, n_v_irr, n_p_irr, ua, L, uaL, memstat

    n_ua = size(BasDim)
    ASSERT(n_ua>0)
    n_v_irr = BasDim(1)%n_irr ! FIXME: ugly

    !
    ! Count shells:
    !
    uaL = 0
    do ua = 1, n_ua
       do L = 0, BasDim(ua)%lmax
          uaL = uaL + 1
       enddo
    enddo
    uaL_max = uaL

    !
    ! Shell dimensions (regular, no spin-orbit case):
    !
    allocate(uaL_nrad(uaL_max), STAT=memstat)
    ASSERT(memstat==0)

    allocate(uaL_vec_mult(uaL_max, n_v_irr), STAT=memstat)
    ASSERT(memstat==0)

    allocate(uaL_vec_dims(uaL_max, n_v_irr), STAT=memstat)
    ASSERT(memstat==0)

    uaL = 0
    do ua = 1, n_ua
       do L = 0, BasDim(ua)%lmax ! SymDim(ua)%lmax may differ
          uaL = uaL + 1

          !
          ! uaL_vec_dims(uaL, irrep) =
          !     = uaL_vec_mult(uaL, irrep) * uaL_nrad(uaL)
          !

          !
          ! Radial (e.g. contractions):
          !
          uaL_nrad(uaL) = BasDim(ua)%n_rad(L)

          !
          ! This is the multiplicity of the shell, that is
          ! the number of times an irrep occurs:
          !
          uaL_vec_mult(uaL, :) = SymDim(ua)%LM(L, :)

          !
          ! Total shell dimension:
          !
          uaL_vec_dims(uaL, :) = uaL_vec_mult(uaL, :) * uaL_nrad(uaL)
       enddo
    enddo

    if ( SpinOrbit ) then
       n_p_irr = BasDimSpor(1)%n_irr

       allocate(uaL_proj_dims(uaL_max, n_p_irr), STAT=memstat)
       ASSERT(memstat==0)

       !
       ! Shell dimensions, spin-orbit:
       !
       uaL = 0
       do ua = 1, n_ua
          do L = 0, BasDimSpor(ua)%lmax
             uaL = uaL + 1
             uaL_proj_dims(uaL, :) = BasDimSpor(ua)%LM(L, :)
          enddo
       enddo
       ASSERT(uaL==size(uaL_proj_dims,1))
    endif
  end subroutine CalcUALDims

  subroutine CalcBasDim(sd, bd, irrs, CNT)
    use unique_atom_module, only: unique_atoms
    implicit none
    type(SubSpaceDim),intent(in)    :: sd(:)
    type(SubSpaceDim),intent(inout) :: bd(:)
    integer(IK), allocatable, intent(out) :: irrs(:)
    logical,optional,intent(in)     :: CNT
    ! *** end of interface ***

    integer(IK) :: n_irr,lmax,l,a,n_rad
    logical     :: contracted

    ASSERT(size(sd)==size(bd))

    if(present(CNT))then
       contracted = CNT
    else
       contracted = .false.
    endif

    n_irr = sd(1)%n_irr

    allocate(irrs(n_irr))

    irrs(:) = 0

    do a=1,size(sd)
       !mdf>>
       lmax  = sd(a)%lmax
       if(lmax/=unique_atoms(a)%lmax_ob)then
          ! FIXME: not a real lmax:
          lmax = unique_atoms(a)%lmax_ob
       endif


       call alloc(lmax,n_irr,bd(a))

       bd(a)%LM = sd(a)%LM(:lmax,:)

      ! sd(A)%LM(L,Irr)
      ! is a independent function number for irrep Irr
      ! (dim(Irr)) of the LM shell of
      ! atom A after symmetry reduction:

       do l=0,lmax

          if(contracted)then
             n_rad = unique_atoms(a)%l_ob(l)%n_contracted_fcts&
                  & +unique_atoms(a)%l_ob(l)%n_uncontracted_fcts
          else
             n_rad = unique_atoms(a)%l_ob(l)%n_exponents
          endif

          !added>>>
          bd(a)%n_rad(l) = n_rad
          !<<<added
          bd(a)%LM(l,:)  = bd(a)%LM(l,:) * n_rad

          irrs(:) = irrs(:) + bd(a)%LM(l,:)
       enddo
    enddo
  end subroutine CalcBasDim

  subroutine print_arr_SubSpaceDim(d)
    implicit none
    type(SubSpaceDim),intent(in) :: d(:)
    ! *** end of interface ***

    integer(IK) :: i,n_ua,irr,l !!$,s,twoj

    print *, 'dim/print_arr_SubSpaceDim: entered'

    n_ua = size(d)

    do i=1,n_ua
             write(*,*)
             write(*,'("UNIQUE ATOM",I5)') i
       
       do l=0,d(i)%lmax
             write(*,'("---")')
             write(*,'("SPLITTING OF L =",I5)') l
             write(*,'("Irrep ",(10I5))') (irr,irr=1,d(i)%n_irr)
             write(*,'("Dim   ",(10I5))') d(i)%lm(l,:)
       enddo
    enddo
  end subroutine print_arr_SubSpaceDim

  !-------------------------------------------------------
  !
  ! memory managment for SubSpaceDim type >>>
  !

  subroutine alloc_SubSpaceDim(lmax,n_irr,bd)
    implicit none
    integer(IK),intent(in)          :: lmax,n_irr
    type(SubSpaceDim),intent(inout) :: bd
    ! *** end of interface ***

    integer :: memstat

    bd%lmax  = lmax
    bd%n_irr = n_irr

    allocate(bd%LM(0:lmax,n_irr), STAT=memstat)
    ASSERT(memstat==0)

    allocate(bd%n_rad(0:lmax),STAT=memstat)
    ASSERT(memstat==0)
    bd%n_rad = 1

    bd%LM = 0
  end subroutine alloc_SubSpaceDim

  subroutine dealloc_SubSpaceDim(bd)
    implicit none
    type(SubSpaceDim),intent(inout) :: bd
    ! *** end of interface ***

    integer :: memstat

    bd%lmax  = -1
    bd%n_irr = -1

    deallocate(bd%LM, STAT=memstat)
    ASSERT(memstat==0)

    deallocate(bd%n_rad,STAT=memstat)
    ASSERT(memstat==0)
  end subroutine dealloc_SubSpaceDim

  subroutine dealloc_arr_SubSpaceDim(bd)
    implicit none
    type(SubSpaceDim),intent(inout) :: bd(:)
    ! *** end of interface ***

    integer(IK) :: a

    do a=1,size(bd)
       call dealloc(bd(a))
    enddo
  end subroutine dealloc_arr_SubSpaceDim

  !--------------- End of module ----------------------------------
end module dimensions

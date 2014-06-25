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
module xc_cntrl
  !---------------------------------------------------------------
  !
  ! Copyright (c) Alexei Matveev
  ! Copyright (c) Thomas Soini
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
  use type_module, only: IK=>i4_kind, RK=>r8_kind ! type specification parameters
  use strings,     only: strlen=>stringlength_string, rtoa
  implicit none
  save            ! save all variables defined in this module
  private         ! by default, all names are private
  !== Interrupt end of public interface of module =================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ------------
  integer(IK)                    :: k
  ! FIXME: this was a bad idea ...
  integer(IK), parameter, public ::&
       xc_XAlpha              =  1, &
       xc_VWN                 =  2, &
       xc_PWLDAc              =  3, &
       xc_RXalpha             =  4, &
       xc_RVWN                =  5, &
       ! this is a PARAMETER, not an option:
       xc_NLDAXC              =  5, & ! Number of LDA XC functionals
       ! options continue:
       xc_Becke88             =  6, &
       xc_RBecke88            =  7, & ! IS A TOTAL, NOT A CORRECTION
       xc_PerdewWang91x       =  8, &
       xc_RPerdewwang91x      =  9, &
       xc_PBEx                = 10, &
       xc_revPBEx             = 11, &
       xc_pbesolx             = 12, &
       xc_PBENx               = 13, &
       xc_Perdew              = 14, &
       xc_PerdewWang91c       = 15, & ! may be the pucking bug i=29 on theo1
       xc_revPW91c            = 16, &
       xc_RPerdewwang91c      = 17, &
       xc_PBEc                = 18, &
       xc_pbesolc             = 19, &
       xc_Baerends94          = 20, &
       xc_ECMV92              = 21, & ! IS A TOTAL, NOT A CORRECTION
       xc_NRECMV92            = 22, & ! IS A TOTAL, NOT A CORRECTION
       xc_HCTH_X              = 23, &
       xc_HCTH_C              = 24, &
       xc_M06L_X              = 25, & ! is a total, not a correction
       xc_M06L_C              = 26, & ! is a total, not a correction
       xc_M06_X               = 27, & ! is a total, not a correction
       xc_M06_C               = 28, & ! is a total, not a correction
       xc_TPSS_X              = 29, & ! IS A TOTAL, NOT A CORRECTION
       xc_TPSS_C              = 30, & ! IS A TOTAL, NOT A CORRECTION
!      xc_revTPSS_X           = 31, & ! IS A TOTAL, NOT A CORRECTION
!      xc_revTPSS_C           = 32, & ! IS A TOTAL, NOT A CORRECTION
       xc_VMT                 = 31, & !
       xc_VT84                = 32, & !
       xc_VMTsol              = 33, & !
       xc_VT84sol             = 34, & !
       xc_VS_X                = 35, & ! is a total, not a correction
       xc_VS_C                = 36, & ! is a total, not a correction
       xc_LYP_C               = 37, & ! is a total, not a correction
       xc_EXX                 = 38, &
       xc_HF                  = 39, & ! mainly to set exact coulomb for HF
       xc_new                 = 40, & ! THE LAST OF XC, copy to default(xc_NXC)
       ! this is a PARAMETER, not an option:
       xc_NXC                 = 40, & ! Number of XC functionals
       ! options continue:
       xc_GGA                 = 41, &
       xc_GGA_version         = 42, &
       xc_MGGA                = 43, &
       xc_HYBRID              = 44, & ! may be obsolete or can substitute xc_EXX
       xc_SO_Spatial          = 45, &
       xc_SO_orb_to_sporb     = 46, &
       xc_orbcalc_version     = 47, &
       xc_denscalc_version    = 48, & ! bug? xc_denscalc_version /= xc_NRECMV92
       xc_sop_version         = 49, &
       xc_HCTH_version        = 50, &
       xc_nuc_sec_der_version = 51, &
       xc_rel                 = 52, &
       xc_ANY                 = 53, &
       xc_NOp                 = 54  ! Number of Options, incuding this one

  integer(IK), parameter, private :: df_Options(xc_NOp) = (/&
       1,                    & !   1) XALPHA (set!)
       1,                    & !   2) VWN    (set!)
       ( 0, k = 3, xc_NXC ), & ! All others disabled !!!!
       0,                    & !  41) GGA (set if one of the GGA is set)
       2,                    & !  42) GGA_version
       0,                    & !  43) MGGA (set if one of the MGGA is set)
       0,                    & !  44) HYBRID (set if one of the HYBRIDS is set)
       0,                    & !  45) Spatial
       0,                    & !  46) orb_to_sporb
       3,                    & !  47) orbcalc_version (changed to 3 in r1.4)
       2,                    & !  48) denscalc_version
       2,                    & !  49) sop_version
       4,                    & !  50) hcth_version
       2,                    & !  51) nuc_sec_der_version
       0,                    & !  52) rel
       1,                    & !  53) is XC on at all?
       xc_NOp               /) !  54) number of options, including this one

  integer(IK), private :: Options(xc_NOp) = df_Options

  real(RK),    private :: frac(xc_NOp)    = (/ 1.0_RK               & ! XAlpha
                                             , 1.0_RK               & ! VWN
                                             , (0.0_RK, k=3,xc_NOp) /)! others

  real(RK),    private :: full         = 1.0_RK ! For "pure" functionals
  real(RK),    private :: off          = 0.0_RK ! Check if really needed
  real(RK),    private :: frac_tmp

  integer(IK), private :: iyes = 1

  logical, public   :: xc_nl_calc ! == is_on(xc_GGA), DEPRECATED

  real(RK), public  :: xc_sdens_cutoff = 1.0E-20_rk
  real(RK), private :: df_sdens_cutoff = 1.0E-20_rk

  character(len=16), public ::&
       xc_cntrl_names(2, xc_NXC) = reshape( (/ &
       'xalpha          ' , '[X L N ]        ', &
       'vwn             ' , '[ CL N ]        ', &
       'pwldac          ' , '[ CL N ]        ', &
       'rxalpha         ' , '[X    R]        ', &
       'rvwn            ' , '[ C   R]        ', &
       'becke88         ' , '[X  GN ]        ', &
       'rbecke88        ' , '[X  G R]        ', &
       'perdewwang91x   ' , '[X  GN ]        ', &
       'rperdewwang91x  ' , '[X LGNR]        ', &
       'pbex            ' , '[X  GN ]        ', &
       'revpbex         ' , '[X  GN ]        ', &
       'pbesolx         ' , '[X  GN ]        ', &
       'pbenx           ' , '[X  GN ]        ', &
       'perdew          ' , '[ C GN ]        ', &
       'perdewang91c    ' , '[ C GN ]        ', &
       'revpw91c        ' , '[ C GN ]        ', &
       'rperdewwang91c  ' , '[ C GNR]        ', &
       'pbec            ' , '[ C GN ]        ', &
       'pbesolc         ' , '[ C GN ]        ', &
       'baerends94      ' , '[XCLGN ]        ', &
       'ecmv92          ' , '[X LGNR]        ', &
       'nrecmv92        ' , '[X LGN ]        ', &
       'hcth_x          ' , '[X LGN ]        ', &
       'hcth_c          ' , '[ CLGN ]        ', &
       'm06l_x          ' , '[X LGN ]        ', &
       'm06l_c          ' , '[ CLGN ]        ', &
       'm06_x           ' , '[X LGN ]        ', &
       'm06_c           ' , '[ CLGN ]        ', &
       'tpss_x          ' , '[X LGN ]        ', &
       'tpss_c          ' , '[ CLGN ]        ', &
       'vmt             ' , '[X LGN ]        ', &
       'vt84            ' , '[X LGN ]        ', &
       'vmtsol          ' , '[X LGN ]        ', &
       'vt84sol         ' , '[X LGN ]        ', &
       'vs_x            ' , '[X LGN ]        ', &
       'vs_c            ' , '[ CLGN ]        ', &
       'lyp_c           ' , '[ CLGN ]        ', &
       'exx             ' , '[X LGN ]        ', &
       'hf              ' , '[X LGN ]        ', &
       'new             ' , '[      ]        '  &
       /), (/ 2, xc_NXC /) )

  !------------ Interface statements ------------------------------

  !------------ public functions and subroutines ------------------
  public :: xc_read_input
  public :: xc_write_input
  public :: xc_input_bcast
  public :: is_on
  public :: whatis
  public :: xc_cntrl_give
  public :: xc_cntrl_frac

  !================================================================
  ! End of public interface of module
  !================================================================


  !------------ Declaration of types ------------------------------

  !------------ Declaration of constants and variables ----

  ! ----------- XC control parameters ---------------------------
  ! USE OF THESE VARIABLES DEPRECATED:
  logical :: becke88, &
             rbecke88, &
             perdewwang91x , &
             perdewwang91c ,&
             rperdewwang91c ,&
             baerends94, &
             perdew , &
             xalpha , &
             vwn,rvwn , &
             pbex, pbec, revPBEx, PBENx, &
             revPW91c, pwldac, &
             rxalpha, &
             ecmv92, &
             nrecmv92, &
             rperdewwang91x, &
             hcth_x, &
             hcth_c

  ! is going to replace the above:
  character(len=strlen) :: XC = "undef"
  character(len=strlen) :: xc_setting = "undef"

  real(RK), parameter, private :: zero = 0.0_rk

  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function get_contrib(s, t) result(ROp)
    use strings, only : sindex, takeword
    !
    ! Decodes the <XC>:<frac> for custom functionals
    !
    implicit none
    character(len=*),   intent(in) :: s, t
    real(RK)                       :: ROp
    ! *** end of interface ***
    integer(IK)                    :: pos, pose
    !
    pos = sindex( s, t ) + len(t)
    pose= index(s(pos:),' ')
    read(s(pos:pos+pose-1),*) ROp
    !
  end function get_contrib

  subroutine set_contrib(Op, IOp, ROp)
    !
    ! Associates an integer IOp and a real ROp with
    ! an option Op preserving constrains:
    !
    !   is_on(Op) <=> (whatis(Op) /= 0) <=> (xc_cntrl_frac(Op) /= 0.0)
    !
    implicit none
    integer(IK),        intent(in) :: Op
    integer(IK),        intent(in) :: IOp
    real(RK), optional, intent(in) :: ROp
    ! *** end of interface ***
    !
    ASSERT(Op>0.and.Op<=xc_NXC)
    ASSERT(IOp>=0)
    !
    Options( Op ) = IOp
    if     ( IOp > 0 .and. .not. present( ROp ) ) then
      ! default case: frac = 1.0
      frac( Op ) = full
      !
    elseif ( IOp > 0 .and.       present( ROp ) ) then
      ASSERT(ROp/=0.0_RK)
      ! ROp present
      frac( Op ) = frac( Op ) + ROp
      !
    else
      ASSERT(IOp==0)
      ! no ROp
      frac( Op ) = off
    endif
    !
  end subroutine set_contrib

  function is_on(Op) result(Yes)
    implicit none
    integer(IK), intent(in) :: Op
    logical :: Yes
    ! *** end of interface ***

    DPRINT 'xccntl::is_on: Op(',Op,')=',Options(Op)

    !
    ! FIXME: should it rather be an
    !
    !           ASSERT(Op>=1.and.Op<=xc_NXC)
    !
    !        instead?
    !
    ASSERT(Op>=1.and.Op<=xc_NOp)

    Yes = (whatis(Op).ne.0)
  end function is_on

  function xc_cntrl_frac(Op) result(ROp)
    implicit none
    integer(IK), intent(in) :: Op
    real(RK) :: ROp
    ! *** end of interface ***

    ASSERT(Op>=1.and.Op<=xc_NOp)

    ROp = frac(Op)

    ASSERT(ROp/=0.0.eqv.is_on(Op))
  end function xc_cntrl_frac

  function whatis(Op) result(IOp)
    implicit none
    integer(IK), intent(in) :: Op
    integer(IK)             :: IOp ! <<< result
    ! *** end of interface ***

    ASSERT(Op>=1.and.Op<=xc_NOp)

    IOp = Options(Op)
  end function whatis

  subroutine set_functional_type()
    ! sets Options(xc_GGA) , Options(xc_MGGA), etc.
    implicit none
    ! *** end of interface ***

    if(      is_on(xc_Perdew)              .or. is_on(xc_Becke88)              &
        .or. is_on(xc_RBecke88)            .or. is_on(xc_PerdewWang91x)        &
        .or. is_on(xc_RPerdewWang91x)      .or. is_on(xc_PerdewWang91c)        &
        .or. is_on(xc_RPerdewWang91c)      .or. is_on(xc_Baerends94)           &
        .or. is_on(xc_PBEx)                .or. is_on(xc_PBEc)                 &
        .or. is_on(xc_pbesolx)             .or. is_on(xc_pbesolc)              &
        .or. is_on(xc_revPBEx)             .or. is_on(xc_PBENx)                &
        .or. is_on(xc_revPW91c)            .or. is_on(xc_HCTH_x)               &
        .or. is_on(xc_HCTH_c)              .or. is_on(xc_ECMV92)               &
        .or. is_on(xc_NRECMV92)            .or. is_on(xc_M06L_x)               &
        .or. is_on(xc_M06L_c)              .or. is_on(xc_M06_x)                &
        .or. is_on(xc_M06_c)               .or. is_on(xc_TPSS_x)               &
        .or. is_on(xc_TPSS_c)              .or. is_on(xc_VMT)                  &
        .or. is_on(xc_VT84)                .or. is_on(xc_VMTsol)               &
        .or. is_on(xc_VT84sol)             .or. is_on(xc_VS_x)                 &
        .or. is_on(xc_VS_c)                .or. is_on(xc_LYP_c)                &
        .or. is_on(xc_new)                                                    )&
             Options(xc_GGA) = iyes

    if(      is_on(xc_M06L_x)              .or. is_on(xc_M06L_c)               &
        .or. is_on(xc_M06_x)               .or. is_on(xc_M06_c)                &
        .or. is_on(xc_TPSS_x)              .or. is_on(xc_TPSS_c)               &
        .or. is_on(xc_VS_x)                .or. is_on(xc_VS_c)                 &
        .or. is_on(xc_new)                                                    )&
             Options(xc_MGGA) = iyes

    if(      is_on(xc_EXX)                 .or. is_on(xc_HF)                  )&
             Options(xc_HYBRID) = iyes

  end subroutine set_functional_type

  subroutine set_global()
    ! compatibility
    implicit none
    ! *** end of interface ***

    DPRINT 'xccntl::set_global: entered'
    xalpha        = is_on(xc_xalpha)
    vwn           = is_on(xc_vwn)
    rvwn          = is_on(xc_rvwn)
    pwldac        = is_on(xc_pwldac)
    rxalpha       = is_on(xc_rxalpha)

    perdew        = is_on(xc_perdew)
    becke88       = is_on(xc_becke88)
    rbecke88      = is_on(xc_rbecke88)
    perdewwang91x = is_on(xc_perdewwang91x)
    perdewwang91c = is_on(xc_perdewwang91c)
    rperdewwang91c = is_on(xc_rperdewwang91c)
    baerends94    = is_on(xc_baerends94)
    pbex          = is_on(xc_pbex)
    pbec          = is_on(xc_pbec)
    ! no legacy flag for PBESOLC, use XC="PBESOL" in input
    revPBEx       = is_on(xc_revPBEx)
    ! no legacy flag for PBESOLX, use XC="PBESOL" in input
    PBENx         = is_on(xc_PBENx)
    revPW91c      = is_on(xc_revPW91c)
    ecmv92        = is_on(xc_ecmv92)
    nrecmv92      = is_on(xc_nrecmv92)
    rperdewwang91x= is_on(xc_rperdewwang91x)
    hcth_x        = is_on(xc_hcth_x)
    hcth_c        = is_on(xc_hcth_c) 

    call set_functional_type()

    xc_nl_calc = is_on(xc_GGA)
    DPRINT 'xccntl::set_global: Options=',Options
    DPRINT 'xccntl::set_global: xc_nl_calc=',xc_nl_calc
    DPRINT 'xccntl::set_global: exit'
  end subroutine set_global

  subroutine set_options()
    ! compatibility
    implicit none
    ! *** end of interface ***

    !MF bugfix:
    Options(:xc_GGA) = 0
    if (xalpha)         call set_contrib( xc_xalpha        , iyes )
    if (vwn)            call set_contrib( xc_vwn           , iyes )
    if (rvwn)           call set_contrib( xc_rvwn          , iyes )
    if (perdew)         call set_contrib( xc_perdew        , iyes )
    if (becke88)        call set_contrib( xc_becke88       , iyes )
    if (rbecke88)       call set_contrib( xc_rbecke88      , iyes )
    if (perdewwang91x)  call set_contrib( xc_perdewwang91x , iyes )
    if (perdewwang91c)  call set_contrib( xc_perdewwang91c , iyes )
    if (baerends94)     call set_contrib( xc_baerends94    , iyes )
    if (pbex)           call set_contrib( xc_pbex          , iyes )
    if (pbec)           call set_contrib( xc_pbec          , iyes )
    ! no legacy flag for PBESOLC, use XC=" PBESOL" in input
    if (revPBEx)        call set_contrib( xc_revPBEx       , iyes )
    ! no legacy flag for PBESOLX, use XC=" PBESOL" in input
    if (PBENx)          call set_contrib( xc_PBENx         , iyes )
    if (revPW91c)       call set_contrib( xc_revPW91c      , iyes )
    if (pwldac)         call set_contrib( xc_pwldac        , iyes )
    if (rxalpha)        call set_contrib( xc_rxalpha       , iyes )
    if (ecmv92)         call set_contrib( xc_ecmv92        , iyes )
    if (nrecmv92)       call set_contrib( xc_nrecmv92      , iyes )
    if (rperdewwang91x) call set_contrib( xc_rperdewwang91x, iyes )
    if (rperdewwang91c) call set_contrib( xc_rperdewwang91c, iyes )
    if (hcth_x)         call set_contrib( xc_hcth_x        , iyes )
    if (hcth_c)         call set_contrib( xc_hcth_c        , iyes )

    call set_functional_type()

    xc_nl_calc = is_on(xc_GGA)
  end subroutine set_options

  subroutine xc_read_input()
    ! Purpose : Reading the input concerning the xc_control
    !           The input consits of logicals which determine if
    !           a special functional has to be calculated or not.
    !           if a nonlocal functional is selected the logical
    !           nl_clac is set to true.
    use iounitadmin_module, only: write_to_output_units
    use input_module
    use strings,            only: wpresent, spresent
    use exchange, only: &
       X_PW91, X_BECKE86, X_PBE, X_REVPBE, X_PBESOL, X_PBEN, X_BECKE88, X_ECMV92, X_XALPHA, &
       X_LONG, X_TRANS, X_LONGTRANS
    implicit none
    !** End of interface **************************************

    integer(kind=IK) :: unit,status
    logical,parameter     :: yes=.true., no=.false.
    character(len=strlen) :: cbuf
    real(rk)              :: sdens_cutoff
    integer(IK)           :: found
    integer(IK)           :: xc_found, xc_fracs_found
    integer(IK)           :: longtrans
    integer(IK)           :: ixc

    namelist /xc_control/ &
         & XC,&
         & xalpha,&
         & rxalpha, &
         & vwn,&
         & rvwn,&
         & pwldac,&
         & becke88,&
         & rbecke88,&
         & perdew,&
         & perdewwang91x,&
         & rperdewwang91x,&
         & perdewwang91c,&
         & rperdewwang91c,&
         & revPW91c,&
         & baerends94,&
         & pbex,&
         & pbec,&
         & revPBEx,&
         & PBENx,&
         & ecmv92, &
         & nrecmv92, &
         & hcth_x, &
         & hcth_c,&
         & sdens_cutoff

    DPRINT 'xccntl::xc_read: entered'

    found = 0

    ! default values
    call set_defaults( yes )
    DPRINT 'xccntl::xc_read: defaults=',Options

    if ( input_line_is_namelist("xc_control") ) then
       call input_read_to_intermediate
       unit= input_intermediate_unit()
       read(unit, nml=xc_control, iostat=status)
       if (status .gt. 0) call input_error( &
            "xc_read: namelist xc_control")
    endif

    if (  xalpha.neqv.(df_Options(xc_xalpha) > 0) .or.                         &
          rxalpha.neqv.(df_Options(xc_rxalpha) > 0) .or.                       &
          vwn.neqv.(df_Options(xc_vwn) > 0) .or.                               &
          rvwn.neqv.(df_Options(xc_rvwn) > 0) .or.                             &
          pwldac.neqv.(df_Options(xc_pwldac) > 0) .or.                         &
          becke88.neqv.(df_Options(xc_becke88) > 0) .or.                       &
          rbecke88.neqv.(df_Options(xc_rbecke88) > 0) .or.                     &
          perdew.neqv.(df_Options(xc_perdew) > 0) .or.                         &
          perdewwang91x.neqv.(df_Options(xc_perdewwang91x) > 0) .or.           &
          rperdewwang91x.neqv.(df_Options(xc_rperdewwang91x) > 0) .or.         &
          perdewwang91c.neqv.(df_Options(xc_perdewwang91c) > 0) .or.           &
          rperdewwang91c.neqv.(df_Options(xc_rperdewwang91c) > 0) .or.         &
          revPW91c.neqv.(df_Options(xc_revPW91c) > 0) .or.                     &
          baerends94.neqv.(df_Options(xc_baerends94) > 0) .or.                 &
          pbex.neqv.(df_Options(xc_pbex) > 0) .or.                             &
          pbec.neqv.(df_Options(xc_pbec) > 0) .or.                             &
          revPBEx.neqv.(df_Options(xc_revPBEx) > 0) .or.                       &
          PBENx.neqv.(df_Options(xc_PBENx) > 0) .or.                           &
          ecmv92.neqv.(df_Options(xc_ecmv92) > 0) .or.                         &
          nrecmv92.neqv.(df_Options(xc_nrecmv92) > 0) .or.                     &
          hcth_x.neqv.(df_Options(xc_hcth_x) > 0) .or.                         &
          hcth_c.neqv.(df_Options(xc_hcth_c) > 0) ) then
      write(*,*) 'WARNING: Explicit functional keyword obsolete!'
      write(*,*) '         Please use XC = combination.'
    endif

    xc_setting = " "
    if( XC .ne. "undef" )then
       ! ignore old fashioned input:
       call set_defaults( no )

       xc_found = 0
       xc_fracs_found = 0

       if(wpresent(XC,"XALPHA"))then
          call set_contrib( xc_xalpha       , iyes )
          xc_found = xc_found + 1
          xc_setting = "XALPHA"
       endif

       if(wpresent(XC,"VWN"))then
          call set_contrib( xc_xalpha       , iyes )
          call set_contrib( xc_vwn          , iyes )
          xc_found = xc_found + 1
          xc_setting = "VWN"
       endif

       if(    wpresent(XC,"PWLDA" ) .or.&
            & wpresent(XC,"PW-LDA") .or.&
            & wpresent(XC,"PWLSD" ) .or.&
            & wpresent(XC,"PW-LSD") )then
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_pwldac        , iyes )
          xc_found = xc_found + 1
          xc_setting = "PWLDA"
       endif

       if(wpresent(XC,"RVWN"))then
          call set_contrib( xc_xalpha        , 0    )
          call set_contrib( xc_rxalpha       , iyes )
          call set_contrib( xc_vwn           , iyes )
          if(wpresent(XC,"RVWN:Engel"))then
            call set_contrib( xc_rvwn        , 2    )
            xc_setting = "RVWN_ENGEL"
          else
            call set_contrib( xc_rvwn        , 1    )
            xc_setting = "RVWN"
          endif
          xc_found = xc_found + 1
       endif

       if(wpresent(XC,"RVWN_ENGEL"))then
          call set_contrib( xc_xalpha        , 0    )
          call set_contrib( xc_rxalpha       , iyes )
          call set_contrib( xc_vwn           , iyes )
          call set_contrib( xc_rvwn          , 2    )
          xc_found = xc_found + 1
          xc_setting = "RVWN_ENGEL"
       endif

       if(wpresent(XC,"LB94"))then
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_vwn           , iyes )
          call set_contrib( xc_baerends94    , iyes )
          xc_found = xc_found + 1
          xc_setting = "LB94"
       endif

       if(wpresent(XC,"BP"))then
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_xalpha        , 0    )
          call set_contrib( xc_vwn           , 0    )
#else
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_vwn           , iyes )
#endif
          call set_contrib( xc_becke88       , iyes )
          call set_contrib( xc_perdew        , iyes )
          xc_found = xc_found + 1
          xc_setting = "BP"
       endif

       if(wpresent(XC,"BLYP").or.spresent(XC,"B-LYP"))then
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_xalpha        , 0    )
#else
          call set_contrib( xc_xalpha        , iyes )
#endif
          call set_contrib( xc_becke88       , iyes )
          call set_contrib( xc_lyp_c         , iyes )
          xc_found = xc_found + 1
          xc_setting = "BLYP"
       endif

       if(wpresent(XC,"B88"))then
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_xalpha        , 0    )
#else
          call set_contrib( xc_xalpha        , iyes )
#endif
          call set_contrib( xc_becke88       , iyes )
          xc_found = xc_found + 1
          xc_setting = "B88"
       endif

       if(wpresent(XC,"RBP"))then
          ! "RBP" means relativistic Becke for the exchange part
          ! and non-relativistic Perdew for the correlation part.
#ifdef WITH_LIBDFTAUTO
          ABORT('recompile w/o -DWITH_LIBDFTAUTO')
#else
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_vwn           , iyes )
          call set_contrib( xc_becke88       , 0    )
          call set_contrib( xc_rbecke88      , iyes )
          call set_contrib( xc_perdew        , iyes )
#endif
          xc_found = xc_found + 1
          xc_setting = "RBP"
       endif

       if(wpresent(XC,"BPW91"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_becke88       , iyes )
          call set_contrib( xc_perdewwang91c , iyes )
          xc_found = xc_found + 1
          xc_setting = "BPW91"
       endif

       if(wpresent(XC,"PW91"))then
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_xalpha        , 0    )
          call set_contrib( xc_pwldac        , 0    )
#else
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
#endif
          call set_contrib( xc_perdewwang91c , iyes )
          call set_contrib( xc_perdewwang91x , iyes )
          xc_found = xc_found + 1
          xc_setting = "PW91"
       endif

       if(wpresent(XC,"PW91c"))then
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_perdewwang91c , iyes )
          xc_found = xc_found + 1
          xc_setting = "PW91c"
       endif

       if(wpresent(XC,"PW91x"))then
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_perdewwang91x , iyes )
          xc_found = xc_found + 1
          xc_setting = "PW91x"
       endif

       if(wpresent(XC,"RPW91x"))then
          ! "RPW91x" is relativistic PW91x.
          call set_contrib( xc_rperdewwang91x, iyes ) ! IS A TOTAL
          xc_found = xc_found + 1
          xc_setting = "RPW91x"
       endif

       if(wpresent(XC,"RPW91c"))then
          ! "RPW91" is relativistic PW91c
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_perdewwang91c , iyes )
          call set_contrib( xc_rperdewwang91c, iyes )
          xc_found = xc_found + 1
          xc_setting = "RPW91c"
       endif

       if(wpresent(XC,"RPW91"))then
          ! "RPW91" is relativistic PW91c and relativistic PW91x.
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_rperdewwang91x, iyes ) ! IS A TOTAL
          call set_contrib( xc_perdewwang91c , iyes )
          call set_contrib( xc_rperdewwang91c, iyes ) ! IS A CORRECTION
          xc_found = xc_found + 1
          xc_setting = "RPW91"
       endif

       if(wpresent(XC,"NRECMV92"))then
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_nrecmv92      , iyes )
          xc_found = xc_found + 1
          xc_setting = "NRECMV92"
       endif

       if(wpresent(XC,"ECMV92"))then
          call set_contrib( xc_xalpha        , iyes )
          call set_contrib( xc_ecmv92        , iyes )
          xc_found = xc_found + 1
          xc_setting = "ECMV92"
       endif

       if(wpresent(XC,"PBE"))then
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_xalpha        , 0    )
          call set_contrib( xc_pwldac        , 0    )
          call set_contrib( xc_pbex          , iyes )
          call set_contrib( xc_pbec          , iyes )
#else
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_pbex          , X_PBE)
          call set_contrib( xc_pbec          , iyes )
#endif
          xc_found = xc_found + 1
          xc_setting = "PBE"
       endif

       if(wpresent(XC,"PBEX"))then
         if(spresent(XC,"PBEX:"))then
           frac_tmp = get_contrib( XC, 'PBEX:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_xalpha        , 0    )
          call set_contrib( xc_pbex          , iyes , frac_tmp )
#else
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes , frac_tmp )
          endif
          call set_contrib( xc_pbex          , X_PBE, frac_tmp )
#endif
          xc_found = xc_found + 1
          xc_setting = "PBEx"//trim(xc_setting)
       endif

       if(wpresent(XC,"PBEC"))then
         if(spresent(XC,"PBEC:"))then
           frac_tmp = get_contrib( XC, 'PBEC:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
#ifdef WITH_LIBDFTAUTO
          call set_contrib( xc_pwldac        , 0    )
          call set_contrib( xc_pbec          , iyes , frac_tmp )
#else
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes , frac_tmp )
          endif
          call set_contrib( xc_pbec          , iyes , frac_tmp )
#endif
          xc_found = xc_found + 1
          xc_setting = "PBEC"//trim(xc_setting)
       endif

       if(wpresent(XC,"PBEN"))then
         if(spresent(XC,"PBEN:"))then
           frac_tmp = get_contrib( XC, 'PBEN:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
         if(.not.is_on(xc_rxalpha))then
           call set_contrib( xc_xalpha       , iyes  , frac_tmp )
         endif
         if(.not.is_on(xc_vwn))then
           call set_contrib( xc_pwldac       , iyes  , frac_tmp )
         endif
         call set_contrib( xc_pbenx          , X_PBEN, frac_tmp )
         call set_contrib( xc_pbec           , iyes  , frac_tmp )
         xc_found = xc_found + 1
         xc_setting = "PBEN"//trim(xc_setting)
       endif

       if(wpresent(XC,"PBENX"))then
         if(spresent(XC,"PBENx:"))then
           frac_tmp = get_contrib( XC, 'PBENx:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes  , frac_tmp )
          endif
          call set_contrib( xc_pbenx         , X_PBEN, frac_tmp )
          xc_found = xc_found + 1
          xc_setting = "PBENx"//trim(xc_setting)
       endif

       if(wpresent(XC,"PBENC"))then
         if(spresent(XC,"PBENC:"))then
           frac_tmp = get_contrib( XC, 'PBENC:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
         if(.not.is_on(xc_vwn))then
           call set_contrib( xc_pwldac       , iyes  , frac_tmp )
         endif
         call set_contrib( xc_pbec           , iyes  , frac_tmp )
         xc_found = xc_found + 1
         xc_setting = "PBENc"//trim(xc_setting)
       endif

       if(wpresent(XC,"revPBE"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_revpbex       , X_REVPBE)
          call set_contrib( xc_pbec          , iyes )
          xc_found = xc_found + 1
          xc_setting = "revPBE"
       endif

       if(wpresent(XC,"revPBEx"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          call set_contrib( xc_revpbex       , X_REVPBE)
          xc_found = xc_found + 1
          xc_setting = "revPBEx"
       endif

       if(wpresent(XC, "PBESOL"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_pbesolx       , X_PBESOL)
          call set_contrib( xc_pbesolc       , iyes )
          xc_found = xc_found + 1
          xc_setting = "PBEsol"
       endif

       if(wpresent(XC,"HCTHv4").or.wpresent(XC,"HCTH407"))then
          call set_contrib( xc_hcth_x        , iyes )
          call set_contrib( xc_hcth_c        , iyes )
          Options(xc_HCTH_version)  = 4
          xc_found = xc_found + 1
          xc_setting = "HCTHv4"
       endif

       if(wpresent(XC,"HCTHv3").or.wpresent(XC,"HCTH147"))then
          call set_contrib( xc_hcth_x        , iyes )
          call set_contrib( xc_hcth_c        , iyes )
          Options(xc_HCTH_version)  = 3
          xc_found = xc_found + 1
          xc_setting = "HCTHv3"
       endif

       if(wpresent(XC,"HCTHv2").or.wpresent(XC,"HCTH120"))then
          call set_contrib( xc_hcth_x        , iyes )
          call set_contrib( xc_hcth_c        , iyes )
          Options(xc_HCTH_version)  = 2
          xc_found = xc_found + 1
          xc_setting = "HCTHv2"
       endif

       if(wpresent(XC,"HCTHv1").or.wpresent(XC,"HCTH"))then
          call set_contrib( xc_hcth_x        , iyes )
          call set_contrib( xc_hcth_c        , iyes )
          Options(xc_HCTH_version)  = 1
          xc_found = xc_found + 1
          xc_setting = "HCTHv1"
       endif

       if(wpresent(XC,"HCTHx"))then
          call set_contrib( xc_hcth_x        , iyes )
          xc_found = xc_found + 1
          xc_setting = "HCTHx"
       endif

       if(wpresent(XC,"HCTHc"))then
          call set_contrib( xc_hcth_c        , iyes )
          xc_found = xc_found + 1
          xc_setting = "HCTHc"
       endif

       if(wpresent(XC,"TPSS"))then
          call set_contrib( xc_tpss_x        , iyes )
          call set_contrib( xc_tpss_c        , iyes )
          xc_found = xc_found + 1
          xc_setting = "TPSS"
       endif

       if(wpresent(XC,"TPSSX"))then
         if(spresent(XC,"TPSSx:"))then
           frac_tmp = get_contrib( XC, 'TPSSx:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
         xc_found = xc_found + 1
         xc_setting = "TPSSx"//trim(xc_setting)
         call set_contrib( xc_tpss_x         , iyes, frac_tmp )
       endif

       if(wpresent(XC,"TPSSC"))then
         if(spresent(XC,"TPSSc:"))then
           frac_tmp = get_contrib( XC, 'TPSSc:' )
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
           xc_fracs_found = xc_fracs_found + 1
         else
           frac_tmp = 1.0_RK
         endif
         xc_found = xc_found + 1
         xc_setting = "TPSSc"//trim(xc_setting)
         call set_contrib( xc_tpss_c         , iyes, frac_tmp )
       endif

!      if(wpresent(XC,"revTPSS"))then
!         call set_contrib( xc_revtpss_x     , iyes )
!         call set_contrib( xc_revtpss_c     , iyes )
!         xc_found = xc_found + 1
!         xc_setting = "revTPSS"
!      endif

       if(wpresent(XC,"M06L"))then
          call set_contrib( xc_m06l_x        , iyes )
          call set_contrib( xc_m06l_c        , iyes )
          xc_found = xc_found + 1
          xc_setting = "M06L"
       endif

       if(wpresent(XC,"M06"))then
         call set_contrib( xc_m06_x          , iyes )
         call set_contrib( xc_m06_c          , iyes )
         call set_contrib( xc_EXX            , iyes, 0.27_RK  )
         xc_found = xc_found + 1
         xc_setting = "M06"
       endif

       if(wpresent(XC,"VMT").or.wpresent(XC,"VMT1"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_vmt           , iyes )
          call set_contrib( xc_pbec          , iyes )
          xc_found = xc_found + 1
          xc_setting = "VMT"
       endif

       if(wpresent(XC,"VT84").or.spresent(XC,"VT{84}"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_vt84          , iyes )
          call set_contrib( xc_pbec          , iyes )
          xc_found = xc_found + 1
          xc_setting = "VT84"
       endif

       if(wpresent(XC,"VMTsol").or.wpresent(XC,"VMT1sol"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_vmtsol        , iyes )
          call set_contrib( xc_pbec          , iyes )
          xc_found = xc_found + 1
          xc_setting = "VMTsol"
       endif

       if(wpresent(XC,"VT84sol").or.spresent(XC,"VT{84}sol"))then
          if(.not.is_on(xc_rxalpha))then
            call set_contrib( xc_xalpha      , iyes )
          endif
          if(.not.is_on(xc_vwn))then
            call set_contrib( xc_pwldac      , iyes )
          endif
          call set_contrib( xc_vt84sol       , iyes )
          call set_contrib( xc_pbec          , iyes )
          xc_found = xc_found + 1
          xc_setting = "VT84sol"
       endif

       if(wpresent(XC,"VSXC").or.wpresent(XC,"VS98"))then
          call set_contrib( xc_vs_x          , iyes )
          call set_contrib( xc_vs_c          , iyes )
          xc_found = xc_found + 1
          xc_setting = "VS98"
       endif

       if(wpresent(XC,"HF").or.wpresent(XC,"HartreeFock"))then
          call set_contrib( xc_HF            , iyes )
          xc_found = xc_found + 1
          xc_setting = "HartreeFock"
       endif

       if(wpresent(XC,"B3LYP").or.wpresent(XC,"B3-LYP"))then
         ! Exc[B3LYP] = Exc[LDA] + a * ( Ex[EXX] - Ex[LDA] )
         !                       + b *   Ex[B88]              <- Correction
         !                       + c * ( Ec[LYP] - Ec[VWN] )
         !            = (1-a-b)*Ex[LDA] + a*Ex[EXX] + b*Ex[B88]
         !              + (1-c)*Ex[VWN] + c*Ec[LYP]
         !
         ! with :     a = 0.20 ,  b = 0.72 ,  c = 0.81
         !
#ifdef WITH_LIBDFTAUTO
         ! There is 0.72*Ex[LDA] already present in B88
         call set_contrib( xc_XAlpha  , iyes,  0.08_RK ) ! = 1 - a - b
#else
         call set_contrib( xc_XAlpha  , iyes,  0.80_RK ) ! = 1 - a
#endif
         call set_contrib( xc_Becke88 , iyes,  0.72_RK ) ! = b
         call set_contrib( xc_VWN     , iyes,  0.19_RK ) ! = 1 - c
         call set_contrib( xc_LYP_c   , iyes,  0.81_RK ) ! = c
         call set_contrib( xc_EXX     , iyes,  0.20_RK ) ! = a
         xc_found = xc_found + 1
         xc_setting = "B3LYP"
       endif

       if(wpresent(XC,"TPSSh").or.wpresent(XC,"TPSS1TPSS"))then
         call set_contrib( xc_tpss_x        , iyes, 0.9_RK )
         call set_contrib( xc_EXX           , iyes, 0.1_RK )
         call set_contrib( xc_tpss_c        , iyes )
         xc_found = xc_found + 1
         xc_setting = "TPSSh"
       endif

       if(wpresent(XC,"PBE0").or.wpresent(XC,"PBE1PBE"))then
         ! This is PBE0 method of Adamo and Barone (J. Chem. Phys. 110)
         ! with 25% exact exchange contribution.
#ifdef WITH_LIBDFTAUTO
         call set_contrib( xc_XAlpha  , 0              )
         call set_contrib( xc_pwldac  , 0              )
         call set_contrib( xc_pbex    , iyes, 0.75_RK  )
         call set_contrib( xc_pbec    , iyes           )
#else
         if(.not.is_on(xc_rxalpha))then
           call set_contrib( xc_XAlpha, iyes, 0.75_RK  )
         endif
         if(.not.is_on(xc_vwn))then
           call set_contrib( xc_pwldac, iyes           )
         endif
         call set_contrib( xc_pbex,    X_PBE, 0.75_RK  )
         call set_contrib( xc_pbec,     iyes           )
#endif
         call set_contrib( xc_EXX     , iyes, 0.25_RK  )
         xc_found = xc_found + 1
         xc_setting = "PBE0"
       endif

       if(wpresent(XC,"EXX"))then
         ! This is pure exact exchange intended for custom settings of methods
         if(spresent(XC,"EXX:"))then
           frac_tmp = get_contrib( XC, 'EXX:' )
           xc_fracs_found = xc_fracs_found + 1
           xc_setting = ':'//trim(rtoa(frac_tmp,8))//' '//trim(xc_setting)
         else
           frac_tmp = 1.0_RK
         endif
         call set_contrib( xc_EXX     , iyes, frac_tmp )
         xc_setting = "EXX"//trim(xc_setting)
         xc_found = xc_found + 1
       endif

       longtrans = X_LONGTRANS
       if(wpresent(XC,"LONG"))then
          longtrans = X_LONG
       endif

       if(wpresent(XC,"REL").or.wpresent(XC,"RLDA"))then
          ! promote to relativistic, if possible

          Options(xc_rel) = longtrans

          found = 0
          if(is_on(xc_xalpha))then
            call set_contrib( xc_rxalpha     , iyes )
            call set_contrib( xc_xalpha      , 0    )
            found = found + 1
            print *,'xc_cntrl: Xalpha promoted to RXalpha'
            print *,'xc_cntrl: Xalpha turned off'
          endif

          if(is_on(xc_vwn))then
             if(wpresent(XC,"Engel"))then
               call set_contrib( xc_rvwn     , 2    )
                print *,'xc_cntrl: VWNc promoted to RVWNc(Engel)'
             else
               call set_contrib( xc_rvwn     , 1    )
               print *,'xc_cntrl: VWNc promoted to RVWNc'
             endif
             found = found + 1
          endif

          if(is_on(xc_pwldac))then
             print *,'xc_cntrl: no RPWLDAC available, maybe RVWN?'
          endif
       endif

       if(wpresent(XC,"REL").or.wpresent(XC,"RGGA").or.wpresent(XC,"RGGAX"))then
          ! promote to relativistic, if possible

          if(is_on(xc_becke88))then
            call set_contrib( xc_xalpha      , iyes )
            call set_contrib( xc_rxalpha     , 0    ) ! RBECKE88 includes RXALPHA (but not XALPHA itself)
            call set_contrib( xc_rbecke88    , iyes )
            call set_contrib( xc_becke88     , 0    )
            found = found + 1
            print *,'xc_cntrl: BECKE88 promoted to RBECKE88'
            print *,'xc_cntrl: XALPHA switched on'
            print *,'xc_cntrl: RXALPHA switched off'
          endif

          if(is_on(xc_perdewwang91x))then
            call set_contrib( xc_xalpha      , iyes )
            call set_contrib( xc_rxalpha     , 0    )
            call set_contrib( xc_rperdewwang91x, iyes )
            call set_contrib( xc_perdewwang91x , 0    ) ! rperdewwang91x is a TOTAL
            found = found + 1
            print *,'xc_cntrl: PW91x promoted to RPW91x'
            print *,'xc_cntrl: XALPHA switched on'
            print *,'xc_cntrl: RXALPHA switched off'
          endif

          if(is_on(xc_pbex))then
            call set_contrib( xc_xalpha      , iyes )
            call set_contrib( xc_rxalpha     , 0    )
            call set_contrib( xc_pbex        , X_PBE + longtrans )
            found = found + 1
            print *,'xc_cntrl: PBEx promoted to REL:PBEx'
            print *,'xc_cntrl: XALPHA switched on'
            print *,'xc_cntrl: RXALPHA switched off'
          endif

          if(is_on(xc_pbenx))then
            call set_contrib( xc_xalpha      , iyes )
            call set_contrib( xc_rxalpha     , 0    )
            call set_contrib( xc_pbenx       , X_PBEN + longtrans )
            found = found + 1
            print *,'xc_cntrl: PBENx promoted to REL:PBENx'
            print *,'xc_cntrl: XALPHA switched on'
            print *,'xc_cntrl: RXALPHA switched off'
          endif

          if(is_on(xc_revpbex))then
            call set_contrib( xc_xalpha      , iyes )
            call set_contrib( xc_rxalpha     , 0    )
            call set_contrib( xc_revpbex     , X_REVPBE + longtrans )
            found = found + 1
            print *,'xc_cntrl: revPBEx promoted to REL:revPBEx'
            print *,'xc_cntrl: XALPHA switched on'
            print *,'xc_cntrl: RXALPHA switched off'
          endif
       endif

       if(wpresent(XC,"REL").or.wpresent(XC,"RGGA").or.wpresent(XC,"RGGAC"))then
          ! promote to relativistic, if possible

          if(is_on(xc_perdewwang91c))then
            call set_contrib( xc_rvwn          , 0    ) ! FIXME: do I need it?
            call set_contrib( xc_rperdewwang91c, iyes )
            call set_contrib( xc_perdewwang91c , 0    ) ! rperdewwang91x is NOW a TOTAL
            found = found + 1
            print *,'xc_cntrl: PW91c promoted to RPW91c'
            print *,'xc_cntrl: RVWNc switched off'
          endif
       endif

       if(wpresent(XC,"off").or.wpresent(XC,"null").or.wpresent(XC,"none"))then
          ! turn everything off:
          Options(:xc_NXC) = 0
          frac             = off
          WARN('xc was turned OFF!')
          xc_found = xc_found + 1
          xc_setting = 'none" #   "WARNING: XC disabled manually'
       endif
       !
       ! more cases wanted ...
       !
       if(wpresent(XC,"spatial"))then
          ! evaluate xc for spin-orbit 
          ! in orbital (as opposit to spinor) space
          Options(xc_so_spatial) = iyes ! <<< global
       endif
       if(wpresent(XC,"orb->sporb"))then
          ! use (spatial0 orbital to compute spinors,
          ! as opposed to independent evaluation of both
          Options(xc_so_orb_to_sporb) = iyes ! <<< global
       endif
       if(wpresent(XC,"GGAv1"))then
          ! older version of GGA subroutine
          Options(xc_GGA_version) = 1
       endif
       if(wpresent(XC,"ORBv1"))then
          ! deprecated
          ! older version of orbitals_calculate:
          Options(xc_orbcalc_version) = 1
       endif
       if(wpresent(XC,"ORBv2"))then
          ! deprecated:
          Options(xc_orbcalc_version) = 2
       endif
       if(wpresent(XC,"ORBv3"))then
          ! latest:
          Options(xc_orbcalc_version) = 3
       endif
       if(wpresent(XC,"DENSv1"))then
          ! older version of density_calc_nl:
          Options(xc_denscalc_version) = 1
       endif
       if(wpresent(XC,"SOPv1"))then
          ! older version of calc_xcks_polarized
          Options(xc_sop_version) = 1
       endif
       if(wpresent(XC,"SOPv2"))then
          ! newer version of calc_xcks_polarized
          Options(xc_sop_version) = 2
       endif
       if(wpresent(XC,"NSDv2"))then
          ! newer version of calculate_nuc_sec_der
          Options(xc_nuc_sec_der_version) = 2
       endif
       if(wpresent(XC,"NSDv1"))then
          ! older version of calculate_nuc_sec_der
          Options(xc_nuc_sec_der_version) = 1
       endif

       if(wpresent(XC,"new").or.wpresent(XC,"test"))then
         call set_contrib( xc_new, iyes )
         xc_found = xc_found + 1
       endif

       if(sdens_cutoff.ne.df_sdens_cutoff)then
          ! spin-density cutoff changed
          xc_sdens_cutoff = sdens_cutoff
       endif
       DPRINT 'xccntl::xc_read: options=',Options

       if ( xc_found == 0 ) then
         call input_error("xc_read: unknown method set in keyword XC")
       elseif ( xc_found > 1 .and. xc_found .ne. xc_fracs_found ) then
         call input_error("xc_read: ambiguous settings in keyword XC")
       endif

       call set_global()
    else
       call set_options()
    endif

    ! Now determine the type of functional i.e. the step in jakobs ladder
    call set_functional_type()

    ! check if any XC was set:
    Options(xc_ANY)  = 0
    do ixc=1,xc_NXC
      if( is_on(ixc) ) Options(xc_ANY) = iyes
    enddo
    DPRINT  'xc_cntrl: xc_ANY=',Options(xc_ANY)

    ! Check consistency of user input (not complete):
    call check()

    call say ("--> XC settings -->")
    if (XC /= "undef") &
         call say ('XC = "' // trim (xc_setting) // '"  (original: "' // trim (XC) // '")')
    if (is_on (xc_so_spatial)) call say ("using orbital space for XC evaluating")
    if (is_on (xc_so_orb_to_sporb)) call say ("using orbitals for evaluation of spinors")

    if (is_on (xc_GGA))then
       if (whatis (xc_orbcalc_version) == 1) call say ("ORBv1")
       if (whatis (xc_orbcalc_version) == 2) call say ("ORBv2")
       if (whatis (xc_orbcalc_version) == 3) call say ("ORBv3")
       if (whatis (xc_denscalc_version) == 1) call say ("DENSv1")
       if (whatis (xc_denscalc_version) == 2) call say ("DENSv2")
       if (whatis (xc_GGA_version) == 1) call say ("GGAv1")
       if (whatis (xc_GGA_version) == 2) call say ("GGAv2")
    endif
    if (whatis (xc_sop_version) == 1) call say ("SOPv1")
    if (whatis (xc_sop_version) == 2) call say ("SOPv2")
    if (whatis (xc_nuc_sec_der_version) == 1) call say ("NSDv1")
    if (whatis (xc_nuc_sec_der_version) == 2) call say ("NSDv2")

    if (xc_sdens_cutoff /= df_sdens_cutoff) then
       write (cbuf, '(E12.5)') xc_sdens_cutoff
       call say ("spin-density cutoff changed --")
       call say ("new value: " // trim (cbuf))
    endif

    call say ("<--" )

    DPRINT "xccntl::xc_read: options   : ",Options
    DPRINT "xccntl::xc_read: xc_nl_calc: ",xc_nl_calc
    
  contains

    subroutine say (phrase)
      use comm, only: comm_rank
      use iounitadmin_module, only: write_to_output_units
      implicit none
      character (len=*), intent (in) :: phrase
      ! *** end of interface ***

      if (comm_rank() == 0) then
         call write_to_output_units ("XC: " // phrase)
      endif
    end subroutine say

    subroutine set_defaults( set_df )
      implicit none
      logical, intent(in) :: set_df
      ! *** end of interface ***

      if( set_df )then
        ! set default values
        Options            = df_Options
        frac(:xc_NXC)      = real( df_Options(:xc_NXC), RK )
        frac(xc_NXC+1:)    = 0.0_RK ! those have no meaning at all
      else
        ! set all XC settings to zero
        Options(:xc_NXC)   = 0
        Options(xc_NXC+1:) = df_Options(xc_NXC+1:)
        frac               = 0.0_RK
      endif

      sdens_cutoff = df_sdens_cutoff

      call set_global()
    end subroutine set_defaults

    subroutine check()
      implicit none
      ! *** end of interface ***

      ! Check consistency (not complete):
#ifdef WITH_LIBDFTAUTO
          ! dont check
#else
      DPRINT 'xcntl::check: Options=',Options
      if ( is_on(xc_vwn).and.is_on(xc_pwldac))&
           & call input_error("Two different LDA-c specified") 

      if (.not.(vwn.or.pwldac).and.(perdewwang91c.or.pbec)) &
           & call input_error&
           & ("PW91c, PBEC require PWLDAc or VWN for the LDA part of Ec")
      ! do not other functio
#endif

      if(xc_sdens_cutoff.lt.zero)&
           & call input_error&
           & ("Spin density cutoff is negative (&xc_control). Aborting...")
    end subroutine check
  end subroutine xc_read_input

  !*************************************************************

  subroutine xc_write_input(iounit)
    !
    ! Purpose: Writing the input concerning xc_control to the output
    !
    use echo_input_module, only: start, real, flag, intg, strng, stop, &
         echo_level_full
    use operations_module, only: operations_echo_input_level
    implicit none
    integer(kind=IK), intent(in) :: iounit
    !** End of interface *****************************************

    ! ----------- Default values for input parameters -------------
    logical ::&
         df_becke88       = df_Options(xc_becke88)       .ne.0 , &
         df_rbecke88      = df_Options(xc_rbecke88)      .ne.0 , &
         df_perdew        = df_Options(xc_perdew)        .ne.0 , &
         df_perdewwang91x = df_Options(xc_perdewwang91x) .ne.0 , &
         df_rperdewwang91x= df_Options(xc_rperdewwang91x).ne.0 , &
         df_perdewwang91c = df_Options(xc_perdewwang91c) .ne.0 , &
         df_rperdewwang91c= df_Options(xc_rperdewwang91c).ne.0 , &
         df_baerends94    = df_Options(xc_baerends94)    .ne.0 , &
         df_xalpha        = df_Options(xc_xalpha)        .ne.0 , &
         df_vwn           = df_Options(xc_vwn)           .ne.0 , &
         df_rvwn          = df_Options(xc_rvwn)          .ne.0 , & 
         df_pbex          = df_Options(xc_pbex)          .ne.0 , & 
         df_pbec          = df_Options(xc_pbec)          .ne.0 , &
         df_revPBEx       = df_Options(xc_revpbex)       .ne.0 , &
         df_PBENx         = df_Options(xc_pbenx)         .ne.0 , &
         df_revPW91c      = df_Options(xc_revpw91c)      .ne.0 , &
         df_pwldac        = df_Options(xc_pwldac)        .ne.0 , &
         df_rxalpha       = df_Options(xc_rxalpha)       .ne.0 , &
         df_hcth_x        = df_Options(xc_hcth_x)        .ne.0 , &
         df_hcth_c        = df_Options(xc_hcth_c)        .ne.0 , &
         df_ecmv92        = df_Options(xc_ecmv92)        .ne.0 , &
         df_nrecmv92      = df_options(xc_nrecmv92)      .ne.0

    call start("XC_CONTROL","XC_WRITE_INPUT", &
         iounit,operations_echo_input_level)
    call strng('XC            ',trim(xc_setting)//'"  # original: "'//trim(XC)&
              ,'undef"  # original: "undef')
    call real("SDENS_CUTOFF  ",xc_sdens_cutoff,df_sdens_cutoff) 
    call flag("XALPHA        ",xalpha ,df_xalpha )
    call flag("VWN           ",vwn    ,df_vwn    )
    call flag("RVWN          ",rvwn   ,df_rvwn    )
    call flag("PWLDAc        ",pwldac ,df_pwldac    )
    call flag("PERDEW        ",perdew ,df_perdew )
    call flag("BECKE88       ",becke88,df_becke88)
    call flag("RBECKE88       ",rbecke88,df_rbecke88)
    call flag("PERDEWWANG91C ",perdewwang91c,df_perdewwang91c)
    call flag("RPERDEWWANG91c",rperdewwang91c,df_rperdewwang91c) 
    call flag("PERDEWWANG91X ",perdewwang91x,df_perdewwang91x)
    call flag("RPERDEWWANG91x ",rperdewwang91x,df_rperdewwang91x)  
    call flag("BAERENDS94    ",baerends94,df_baerends94)
    call flag("PBEx          ",pbex,df_pbex)
    call flag("PBEc          ",pbec,df_pbec)
    call flag("revPBEx       ",revPBEx,df_revPBEx)
    call flag("PBENx         ",PBENx,df_PBENx)
    call flag("revPW91c      ",revPW91c,df_revPW91c)
    call flag("RXALPHA       ",rxalpha ,df_rxalpha)
    call flag("HCTH_X        ",hcth_x ,df_hcth_x)
    call flag("HCTH_C        ",hcth_c ,df_hcth_c)
    call flag("ECMV92        ",ecmv92,df_ecmv92)
    call flag("NRECMV92      ",nrecmv92,df_nrecmv92)
    call stop()

  end subroutine xc_write_input

  !*****************************************************************************

  subroutine xc_input_bcast()
    ! Purpose : broadcasting the xc_input
    use comm, only: comm_bcast                                                 &
                  , comm_rank
    implicit none
    !** End of interface *******************************************************

    call comm_bcast( Options         )
    call comm_bcast( frac            )
    call comm_bcast( xc_sdens_cutoff )

    if ( comm_rank() /= 0 ) call set_global()

  end subroutine xc_input_bcast

  !*****************************************************************************

  subroutine xc_cntrl_give( &
       lbecke88, &
       lrbecke88, &
       lperdew, &
       lxalpha, &
       lvwn,&
       lrvwn,&
       lperdewwang91x, &
       lrperdewwang91x, &
       lperdewwang91c, &
       lrperdewwang91c, &
       lbaerends94,&
       lrxalpha, &
       lhcth_x, &
       lhcth_c, &
       lecmv92, &
       lnrecmv92 )
    implicit none
    logical, intent(out) :: &
         lbecke88, &
         lrbecke88, &
         lperdew, &
         lxalpha, &
         lvwn, &
         lrvwn, &
         lperdewwang91x, &
         lrperdewwang91x, &
         lperdewwang91c, &
         lrperdewwang91c, &
         lbaerends94, &
         lrxalpha, &
         lhcth_x, &
         lhcth_c, &
         lecmv92, &
         lnrecmv92
    ! *** end of interface ***

    lbecke88       = is_on(xc_becke88)
    lrbecke88      = is_on(xc_rbecke88)
    lperdew        = is_on(xc_perdew)
    lxalpha        = is_on(xc_xalpha)
    lvwn           = is_on(xc_vwn)
    lrvwn          = is_on(xc_rvwn) 
    lperdewwang91x = is_on(xc_perdewwang91x)
    lrperdewwang91x= is_on(xc_rperdewwang91x)
    lperdewwang91c = is_on(xc_perdewwang91c)
    lrperdewwang91c= is_on(xc_rperdewwang91c)
    lbaerends94    = is_on(xc_baerends94)
    lrxalpha       = is_on(xc_rxalpha)  
    lhcth_x        = is_on(xc_hcth_x)
    lhcth_c        = is_on(xc_hcth_c)  
    lecmv92        = is_on(xc_ecmv92)
    lnrecmv92      = is_on(xc_nrecmv92)
  end subroutine xc_cntrl_give

  !--------------- End of module ----------------------------------
end module xc_cntrl

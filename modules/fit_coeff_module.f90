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
module  fit_coeff_module
  !---------------------------------------------------------------
  !-------------- Module specification ---------------------------
  !---------------------------------------------------------------
  !
  !  Purpose: contains the following data and setup routines
  !           connected with them:
  !  Types:
  !           fit             -> contains fitbasis dimensions
  !                              for exchange and charge fit
  !                              basis. In case of numerical
  !                              exch.corr. potential, the 
  !                              dimension for exch. is 0
  !  Variables:
  !           coeff_charge    -> fit coefficients for total and spin density
  !           coeff_spin         this should later be included in
  !                              a chargefit_module 
!:UB[ added
  ! 
  !           coeff_proj      -> charge density projection coefficients
  !                              as required for the potential extended MDA
!:UB]
  !           
  !        coeff_charge_old   -> stores old coeff_charge and coeff_spin
  !        coeff_spin_old        coefficients for mixing
  ! 
  !             charge_norm   -> contains the normalization constants
  !                              of the charge fit functions.
  !
  !           coeff_xc
  !           coeff_en_xc     -> fit coefficients for exchange potential
  !                              they also might go in a some kind
  !                              of exchangefit_module
  !
  !         coeff_xcmda      -> Coulomb-type fit coefficients for the exchange
  !                             potential in case the MDA approach is used
  !                             they subsitute coeff_xc (no mixing is required)
  !
  !        coeff_xc_old      -> stores old exchange fir coefficients
  !                             coeff_xc for the mixing procedure
  !
!:UB[ modified
  !         coeff_deltarho   -> Charge density difference fit coefficients
  !                             as required for the potential extended MDA 
  !
!:UB]
  !  Module called by: prescf, initialize_master, chargefit, main_integral, ...
  !
  !  Author: FN
  !  Date: 10/95
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   12/95
  ! Description: packing and unpacking routines added, restructured
  !
  ! Modification (Please copy before editing)
  ! Author: FN
  ! Date:   8/96
  ! Description: added a routine fit_dimensions_calc that
  !              calculates the size of the fitbases using
  !              data from the unique_atoms_module.
  !
  ! Modification (Please copy before editing)
  ! Author: TB
  ! Date:   9/96
  ! Description: Added fit_coeff_calc_chargefitnorm() and
  !              for debugging fit_coeff_read_chargefitnorm()
  !              that was derived from read_chargefitnorm().
  !              Removed obsolent subroutines read_chargefitnorm(),
  !              fit_read_dimensions() and fit_write_input()
  !
  ! Modification (Please copy before editing)
  ! Author: UB
  ! Date:   6/97
  ! Description: a) replace options_numeric_exch by the extended function
  !                 options_xcmode()
  !              b) In case of the model density approach charge fit
  !                 coefficients for the spin density are required as well
  !              c) Status variables CH_INITIALIZED and XC_INITIALIZED for
  !                 the "DISCARD_INITialized fit coefficients after
  !                 the first scf cycle" option introduced.
  !              d) store and recovering of fit coefficients now done with
  !                 readwriteblocked file "saved_scfstate.dat"
  !
  ! Modification (Please copy before editing) 
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: treatment of exchange coefficients coeff_xc and coeff_en_xc
  !              adapted to the extended MDA option
!:UB[ added
  !              coeff_proj and coeff_deltarho for the extended MDA introduced
!:UB]
  !
  ! Modification (Please copy before editing)
  ! Author: AS
  ! Date:   7/97
  ! Description: ...
  !
  ! Modification (Please copy before editing)
  ! Author: 
  ! Date:   ...
  ! Description: ...
  !
  !----------------------------------------------------------------

  !------------ Modules used --------------------------------------
# include "def.h"
  use type_module ! type specification parameters
  use constants
#ifdef WITH_CORE_DENS
  use operations_module, only: operations_core_density
#endif
#ifdef FPP_AIX_XLF
  use matrix_module, only: matmult
# define MATMUL(a,b) matmult(a,b)
#endif
#define FPP_FIT_COEFF_INIT /* use %functions(:, :) of s- and r2-bases */

  implicit none
  private
  save
  !== Interrupt end of public interface of module =================

  ! ----------- Declaration of Types ----------------------------
  type, public :: fit
     integer(kind=i4_kind)    :: n_ch = -1 ! fitbasis dim. CHARGE
     integer(kind=i4_kind)    :: n_s  = -1 ! S fitbasis dim. CHARGE
     integer(kind=i4_kind)    :: n_r2 = -1 ! R2 fitbasis dim. CHARGE
     integer(kind=i4_kind)    :: n_xc = -1 ! fitbasis dim. EXCH.
     integer(kind=i4_kind)    :: n_cd = -1 ! fitbasis dim. CORE DENSITY
  end type fit

  !------------ Declaration of public constants and variables ----

  integer(i4_kind), parameter, public :: &
       ICHFIT = 1, & ! charge (density) fit
       IXCFIT = 2    ! XC fit, used with MDA

  real(kind=r8_kind),allocatable,target,public, protected :: coeff_charge(:)
  ! coeff_charge(n_fit%n_ch)

  real(kind=r8_kind),allocatable,       public :: coeff_charge_veff(:)
  !
  ! coeff_charge_veff(size(coeff_charge))
  !
  ! used to "offset" the charge coeff, e.g. to redefine the nuclear
  ! attraction V_nuc to the "screened by Coulomb" nuclear attraction
  ! to be used in relativistic transformations of the Coulomb potential

  real(kind=r8_kind),allocatable,target,public :: coeff_charge_eperef(:)
  ! coeff_charge_eperef(n_fit%n_ch)
  real(kind=r8_kind),allocatable,target,public :: coeff_core(:)
  ! coeff_core(n_fit%n_ch)
  real(kind=r8_kind),allocatable,target,public :: coeff_spin(:)
  ! coeff_spin(n_fit%n_ch)
  real(kind=r8_kind),allocatable,target,public :: coeff_charge_old(:)
  ! coeff_charge_old(n_fit%n_xc)
  real(kind=r8_kind),allocatable,target,public :: coeff_spin_old(:)
  ! coeff_spin_old(n_fit%n_xc)
  real(kind=r8_kind), allocatable, public, protected :: charge_norm(:)
  ! charge_norm(n_fit%n_ch)
  real(kind=r8_kind),allocatable,public     :: coeff_xc(:,:)
  ! coeff_xc(n_fit%n_xc,options_n_spin())
  real(kind=r8_kind),allocatable,public     :: coeff_xcmda(:,:)
  ! coeff_xcmda(n_fit%n_ch,options_n_spin())
!!! MF for testing
  real(kind=r8_kind),allocatable,public     :: coeff_xcmda_old(:,:)
  !total/spin contribution
  real(kind=r8_kind),allocatable,public     :: coeff_xcmda_ts(:,:)
  ! coeff_xcmda(n_fit%n_ch,options_n_spin())
  real(kind=r8_kind),allocatable,public     :: coeff_xc_old(:,:)
  ! coeff_xc_old(n_fit%n_xc,options_n_spin())
  real(kind=r8_kind),allocatable,public     :: coeff_en_xc(:)
  ! coeff_en_xc(n_fit%n_xc)
!:UB[ modified 
  real(kind=r8_kind),allocatable,public     :: coeff_proj(:,:)
  ! coeff_proj(n_fit%n_xc,options_n_spin())
  real(kind=r8_kind),allocatable,public     :: coeff_deltarho(:,:)
  ! coeff_deltarho(n_fit%n_xc,options_n_spin())
  real(kind=r8_kind),allocatable,public     :: coeff_deltarho_old(:,:)

  logical, public :: pv_initialized
!:UB]

  logical,allocatable,public                :: fitfct_map(:)
  ! fitfct_map(n_fit%n_ch)
  logical, public :: dr_initialized, drold_initialized
  logical, public :: ch_initialized, chold_initialized
  logical, public :: sp_initialized, spold_initialized
  logical, public :: xc_initialized, xcold_initialized

  type, public :: ff_info
     integer(i4_kind) :: UA = -1
     integer(i4_kind) ::  L = -2 ! -1, 0 , l, ... for S, R2, L, ...
     real(r8_kind)    :: COUL_NORM = 0.0
     ! COUL_NORM(i) = [fi||fi], integrals used to be scaled
     ! with  1/sqrt([fi||fi])
     ! if scaled, then COUL_NORM==1.0, identically
     integer(i4_kind) :: next_ual = -1
     ! shows where the next UAxL block in fit_coeff_ff_map(:)
     ! starts, not valid if > size(fit_coeff_ff_map)

     real(r8_kind)    :: aexp = -1.0 ! approximate exponent for use
     ! in cutoffs/screening
     real(r8_kind)    :: Z = -1.0 ! nuclear charge of UA for use
     ! in cutoffs/screening
  end type ff_info

  type(ff_info), allocatable :: fit_coeff_ff_map(:) ! (N_FF)
  ! mapping of ff -> (UA,L)
  ! UA(ff)        == fit_coeff_ff_map(ff)%UA
  !  L(ff)        == fit_coeff_ff_map(ff)%L
  ! Coul_Norm(ff) == fit_coeff_ff_map(ff)%coul_norm
  ! for quick search:
  ! ff_next       == fit_coeff_ff_map(ff)%next_ual
  ! starting index of the next UAxL block
  public :: fit_coeff_ff_map

  !------------ public functions and subroutines ------------------
  public :: get_fit
  public :: fit_coeff_n_ch
  public :: fit_coeff_n_s
  public :: fit_coeff_n_r2
  public :: fit_coeff_n_xc
  public :: copy_coeff_charge_old
  public :: copy_coeff_xc_old
  public :: print_coeff_charge
  public :: free_coeff_charge_old
  public :: free_coeff_xc_old
  public :: fit_coeff_calc_chargefitnorm
  public :: fit_coeff_normalize
  public :: fit_coeff_initialize
  public :: fit_coeff_dimensions_calc
  public :: fit_coeff_charge_norm
  public :: fit_coeff_shutdown
  public :: fit_coeff_store
  public :: fit_coeff_recover
  public :: fit_coeff_n_cd
  public :: fit_coeff_set_ch

#ifdef WITH_CORE_DENS
  public :: write_coeff_core
#endif

  public :: fit_coeff_sndrcv                  ! from parallel context
  public :: fit_coeff_send, fit_coeff_receive ! MUSTDIE: from different contexts

  public :: fit_coeff_store_norm
  public :: fit_coeff_bcast!(), do not confuse with fit_coeff_module_bcast
  public :: fit_coeff_setup!()

!================================================================
! End of public interface of module
!================================================================

  !------------ Declaration of private constants and variables ----
  type(fit)                      :: n_fit 
  !----------------------------------------------------------------
  !------------ Subroutines ---------------------------------------
contains

  function ff_charge(L,alpha) result(charge)
    !
    ! Integral charge of a fitfunciton, e.g. for s-type:
    !
    !                /
    !               |           2   3    / pi \ 3/2
    ! charge(a,S) = | exp( - a r ) d r =| ---- |
    !               |                    \  a /
    !              /
    !
    ! L = s or r (for r^2)
    !
    ! This is NOT the Coulomb self-interaction norm!
    !
    implicit none
    character(len=1), intent(in) :: L
    real(r8_kind), intent(in)    :: alpha
    real(r8_kind)                :: charge ! result
    ! *** end of interface ***

    select case(L)
    case('s') ! s-type
       charge = (PI/alpha)**(three/two)
    case('r') ! for r^2
       charge = (three/two) * PI**(three/two) * (one/alpha)**(five/two)
    case default
       charge = 0
       ABORT('why asking?')
    end select
  end function ff_charge

  subroutine fit_coeff_store_norm(triang)
    ! triang: [fi||fj] ints in tanglular storage
    ! picks up diagonal elements and stores them
    ! in fit_coeff_ff_map
    implicit none
    real(r8_kind), intent(in) :: triang(:)
    ! *** end of interface ***

    integer(I4_kind) :: i,j,ij,ii,n_ch,n_tri

    DPRINT 'fit_coeff_store_norm: entered'

    n_ch = n_fit%n_ch
    n_tri = (n_ch * (n_ch+1)) / 2

    ASSERT(size(triang)==n_tri)
    ASSERT(allocated(fit_coeff_ff_map))
    ASSERT(n_ch==size(fit_coeff_ff_map))

    ij = 0
    do i=1,n_ch
       do j=1,i
          ij = ij + 1
          if(i/=j) cycle
          ii = i ! == j
          fit_coeff_ff_map(ii)%coul_norm  = triang(ij)
       enddo
    enddo
  end subroutine fit_coeff_store_norm

  !*************************************************************
  subroutine get_fit(fit_dummy)
    !  Purpose: make the private variable n_fit accessible
    !           to the calling unit
    !------------ Declaration of formal parameters ---------------
    type(fit),intent(out)      :: fit_dummy
    !** End of interface *****************************************
    !------------ Executable code --------------------------------

    fit_dummy = n_fit ! copy all fields, FIXME: if ever pointers
  end subroutine get_fit
  !*************************************************************

  !*************************************************************
  integer function fit_coeff_n_ch()
    ! Purpose: returns dimension of charge fitfunctions
    !** End of interface ***************************************
    fit_coeff_n_ch = n_fit%n_ch
  end function fit_coeff_n_ch
  !*************************************************************

  !*************************************************************
  integer function fit_coeff_n_s()
    ! Purpose: returns dimension of S charge fitfunctions
    !** End of interface ***************************************
    fit_coeff_n_s = n_fit%n_s
  end function fit_coeff_n_s
  !*************************************************************

  !*************************************************************
  integer function fit_coeff_n_r2()
    ! Purpose: returns dimension of R2 charge fitfunctions
    !** End of interface ***************************************
    fit_coeff_n_r2 = n_fit%n_r2
  end function fit_coeff_n_r2
  !*************************************************************

  !*************************************************************
  integer function fit_coeff_n_xc()
    ! Purpose: returns dimension of exchange fitfunctions
    !** End of interface ***************************************
    fit_coeff_n_xc = n_fit%n_xc
  end function fit_coeff_n_xc
  !*************************************************************

  !*************************************************************
  integer function fit_coeff_n_cd()
    ! Purpose: returns dimension of core density fitfunctions
    !** End of interface ***************************************
    fit_coeff_n_cd = n_fit%n_cd
  end function fit_coeff_n_cd
  !*************************************************************

  subroutine fit_coeff_setup()
    use comm, only: comm_rank
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        !
        ! This fills the structure that is later broadcasted ...
        !
        call fit_coeff_dimensions_calc()
    endif

    call fit_coeff_module_bcast()
  end subroutine fit_coeff_setup

  subroutine fit_coeff_module_bcast()
    !
    ! Do not confuse with fit_coeff_bcast()
    !
    use comm, only: comm_rank
    use comm_module, only: comm_init_send, comm_send, &
        comm_all_other_hosts, comm_save_recv, &
        comm_master_host
    use msgtag_module, only: msgtag_packed_message
    implicit none
    ! *** end of interface ***

    if ( comm_rank() == 0 ) then
        call comm_init_send(comm_all_other_hosts, msgtag_packed_message)
        call fit_coeff_pack_n_fit()
        call comm_send()
    else
        call comm_save_recv(comm_master_host, msgtag_packed_message)
        call fit_coeff_unpack_n_fit()
    endif
  end subroutine fit_coeff_module_bcast

  !*************************************************************
  subroutine fit_coeff_pack_n_fit()
    !  Purpose: packs dimensions of fitfunctions
    use xpack
    !** End of interface ***************************************

    call pck(n_fit%n_xc)
    call pck(n_fit%n_ch)
    call pck(n_fit%n_cd)
    call pck(n_fit%n_s)
    call pck(n_fit%n_r2)

    ASSERT(allocated(fit_coeff_ff_map))
    call pck(fit_coeff_ff_map(:)%UA)
    call pck(fit_coeff_ff_map(:)%L)
    call pck(fit_coeff_ff_map(:)%next_ual)
    call pck(fit_coeff_ff_map(:)%aexp)
    call pck(fit_coeff_ff_map(:)%Z)
    ! FIXME: norms?

  end subroutine fit_coeff_pack_n_fit
  !*************************************************************


  !*************************************************************
  subroutine fit_coeff_unpack_n_fit()
    !  Purpose: packs dimensions of fitfunctions
    use xpack
    !** End of interface ***************************************

    integer(i4_kind) :: memstat

    call upck(n_fit%n_xc)
    call upck(n_fit%n_ch)
    call upck(n_fit%n_cd)
    call upck(n_fit%n_s)
    call upck(n_fit%n_r2)

    ! i may be a slave only?
    allocate(fit_coeff_ff_map(n_fit%n_ch),stat=memstat)
    ASSERT(memstat==0)

    call upck(fit_coeff_ff_map(:)%UA)
    call upck(fit_coeff_ff_map(:)%L)
    call upck(fit_coeff_ff_map(:)%next_ual)
    call upck(fit_coeff_ff_map(:)%aexp)
    call upck(fit_coeff_ff_map(:)%Z)
    ! FIXME: norms?

  end subroutine fit_coeff_unpack_n_fit
  !*************************************************************

  !*************************************************************
  subroutine fit_coeff_bcast()
    ! broadcase data to slaves
    ! called from "integral_shutdown_2cff"
    ! at the time where coulomb norms should be ready
    use comm, only: comm_bcast
    implicit none

    DPRINT 'fit_coeff_bcast: entered'
    ASSERT(allocated(fit_coeff_ff_map))

    !
    ! FIXME: array temporary is created here:
    !
    call comm_bcast(fit_coeff_ff_map(:)%COUL_NORM)
  end subroutine fit_coeff_bcast
  !*************************************************************

  !*************************************************************
  subroutine fit_coeff_calc_chargefitnorm()
    ! Purpose : calculates the one-center integral
    !          <1|fk>, the norm of the fit functions.
    !          This norm is used in the charge fit
    !          The result is stored in the variable
    !          charge_norm of this module. This variable
    !          is also allocated here
    ! Subroutine called by: main_integral
    ! Author: TB, 8/96
    use unique_atom_module
    use symmetry_data_module, only: get_totalsymmetric_irrep
    implicit none
    !** End of interface *****************************************

    !------------ Declaration of local variables --------------
    integer(kind=i4_kind) :: i_ir,i_ua,i_l,alloc_stat,n_exp, &
         i_cn,i_cn_new,n_uc,n_c,i_gc,i_cf,i_exp
    type(unique_atom_type), pointer :: ua
    type(unique_atom_basis_type), pointer :: uab
    type(unique_atom_partner_type), pointer :: uap
    type(unique_atom_glob_con_type), pointer :: uag
    real(kind=r8_kind), pointer :: exps(:)
    real(kind=r8_kind) :: exp, coef, sum, const, const2
    real(kind=r8_kind), allocatable :: intermediate(:)
    !----------- Executable code ------------------------------

    const=three/two*pi*sqrt(pi)
    const2=pi*sqrt(pi)

    allocate(charge_norm(n_fit%n_ch),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("fit_coeff_calc_chargefitnorm: allocation (1) failed")

    i_ir = get_totalsymmetric_irrep()
    i_cn = 1
    unique_atom: do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       ! --- first the s-types ---------------------------
       uab => ua%l_ch(0)
       exps => uab%exponents
       n_exp = uab%N_exponents
       n_uc = uab%N_uncontracted_fcts
       n_c = uab%N_contracted_fcts
       allocate(intermediate(n_exp),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_calc_chargefitnorm: allocation (2) failed")
       intermediate = const2 / ( exps * sqrt(exps) )
       intermediate=intermediate*ua%N_equal_atoms
       i_cn_new = i_cn + n_uc
       charge_norm(i_cn:i_cn_new-1) = &
            intermediate(1:n_uc)
       i_cn = i_cn_new

       if ( n_c .gt. 0 ) then
          i_cn_new = i_cn + n_c
          charge_norm(i_cn:i_cn_new-1) = &
               & MATMUL(intermediate(:),uab%contractions(:,:))
          i_cn = i_cn_new
       endif
       deallocate(intermediate,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_calc_chargefitnorm: deallocation (1) failed")
       ! --- now the r2-types ----------------------------
       uab => ua%r2_ch
       exps => uab%exponents
       n_exp = uab%N_exponents
       n_uc = uab%N_uncontracted_fcts
       n_c = uab%N_contracted_fcts
       allocate(intermediate(n_exp),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_calc_chargefitnorm: allocation (3) failed")
       intermediate = const / ( exps * exps * sqrt(exps) )
       ! symmetry adaption
       intermediate=intermediate*ua%N_equal_atoms
       i_cn_new = i_cn + n_uc
       charge_norm(i_cn:i_cn_new-1) = &
            intermediate(1:n_uc)
     !  print*,'c_exp in fit_coeff_calc_chargefitnorm=',n_uc,uar%c_exp(1:n_uc)
       i_cn = i_cn_new
       if ( n_c .gt. 0 ) then
          i_cn_new = i_cn + n_c
          charge_norm(i_cn:i_cn_new-1) = &
               MATMUL(intermediate,uab%contractions)
          i_cn = i_cn_new
       endif
       deallocate(intermediate,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_calc_chargefitnorm: deallocation (2) failed")
       ! --- fill higher l with zero  ----------------------
       i_cn = i_cn_new
       do i_l = 1, ua%lmax_ch
          uap => ua%symadapt_partner(i_ir,i_l)
          uab => ua%l_ch(i_l)
          i_cn_new = i_cn_new + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
               * uap%N_independent_fcts
       enddo
       charge_norm(i_cn:i_cn_new-1) = zero
       i_cn = i_cn_new
    enddo unique_atom
    ! global contractions
    ua_gc: do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       gc: do i_gc = 1, ua%N_glob_cons_ch
          uag => ua%glob_con_ch(i_gc)
          sum = zero
          cf: do i_cf = 1, uag%N_contributing_fcts
             select case (uag%l(i_cf))
             case (-1)
                i_exp = uag%index_exp(i_cf)
                exp = ua%l_ch(0)%exponents(i_exp)
                coef = uag%coefs(i_cf) * ua%N_equal_atoms
                sum = sum + coef * const2 / ( exp * sqrt(exp) )
             case (0)
                i_exp = uag%index_exp(i_cf)
                exp = ua%r2_ch%exponents(i_exp)
                coef = uag%coefs(i_cf) * ua%N_equal_atoms
                sum = sum + coef * const / ( exp * exp * sqrt(exp) )
             end select
          enddo cf
          charge_norm(i_cn) = sum
          i_cn = i_cn + 1
       enddo gc
    enddo ua_gc
  end subroutine fit_coeff_calc_chargefitnorm
  !*************************************************************

!:UB[ added
  !*************************************************************
  subroutine fit_coeff_allocate()
    !  Purpose: allocates all permanent fit coefficients
    !           called by: fit_coeff_read and fit_coeff_initialized
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use options_module
    use spin_orbit_module, only: whatis, op_BackTrafo
    !------------ Declaration of local variables  ----------------
    integer(kind=i4_kind) :: alloc_stat, n_spin
    logical               :: exchange_fit, model_density, extended_mda
    external error_handler
    !------------ Executable code --------------------------------
    n_spin        = options_n_spin()
    exchange_fit  = options_xcmode() == xcmode_exchange_fit
    model_density = options_xcmode() == xcmode_model_density
    extended_mda  = options_xcmode() == xcmode_extended_mda

    allocate(coeff_charge(n_fit%n_ch),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("FIT_COEFF_ALLOCATE: allocation of coeff_charge failed!")

    if(whatis(op_BackTrafo).eq.2)then
       allocate(coeff_charge_veff(n_fit%n_ch),STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    endif

#ifdef WITH_CORE_DENS
    if (operations_core_density) then
       allocate(coeff_core(n_fit%n_ch),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_READ: allocation (1c) failed!")
       allocate(fitfct_map(n_fit%n_ch),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_READ: allocation (1f) failed!")
       call calc_fitfct_map(fitfct_map)
    endif
#endif
    if (model_density .or. extended_mda) then
       if (n_spin > 1) then
          allocate(coeff_spin(n_fit%n_ch),STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("FIT_COEFF_ALLOCATE: allocation of coeff_spin failed!")
       end if
       allocate(coeff_xcmda(n_fit%n_ch,n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_xcmda failed!")
       allocate(coeff_xcmda_ts(n_fit%n_ch,n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_xcmda failed!")
       coeff_xcmda_ts=zero
!!! MF for testing
     if(.not. allocated(coeff_xcmda_old)) then
       allocate(coeff_xcmda_old(n_fit%n_ch,n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_xcmda_old failed!")
       coeff_xcmda_old=zero
     endif
    endif
    if (extended_mda) then
       allocate(coeff_proj(n_fit%n_xc,n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_proj failed!")
       allocate(coeff_deltarho(n_fit%n_xc,n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_deltarho failed!")
    endif
    if (exchange_fit .or. extended_mda) then
       allocate(coeff_xc(n_fit%n_xc,n_spin),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_xc failed!")
    endif
    if (exchange_fit) then
       allocate(coeff_en_xc(n_fit%n_xc),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("FIT_COEFF_ALLOCATE: allocation of coeff_en_xc failed!")
    endif
  end subroutine fit_coeff_allocate
  !*************************************************************

  !*************************************************************
  subroutine fit_coeff_store(th,mode)
    ! Purpose: stores the fit coefficients on a readwriteblocked file
    !          in case the input switch 'save_scfstate' is set.
    !
    ! << th >>     << mode >>   action
    ! PRESENT      NOT PRESENT  store present state immediatelly
    ! NOT PRESENT  PRESENT      keep present state for later storage
    ! PRESENT      PRESENT      store the previously saved state
    !
    ! subroutine called by: 'main_scf'
    !
    ! UB, 9/97
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use options_module
    use readwriteblocked_module
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module     , only: output_main_scf, output_data_saved
    use spin_orbit_module , only: whatis, op_BackTrafo 
    !------------ Declaration of formal parameters ---------------
    type(readwriteblocked_tapehandle), optional, intent(inout) :: th
    integer(kind=i4_kind), optional, intent(in) :: mode
    !------------ Declaration of local variables  ----------------
    integer(kind=i4_kind) :: s,n_spin
    real(kind=r8_kind) :: yes(1) = (/1.0_r8_kind/), no(1) = (/0.0_r8_kind/)
    real(r8_kind), save, pointer :: coeff_charge_kept(:), coeff_spin_kept(:), &
                                    coeff_xc_kept(:,:), coeff_en_xc_kept(:)
    logical, save :: charge_and_spin_kept = .false., &
                     xc_and_en_xc_kept    = .false.
    logical  :: reset, spin_required, exchange_fit
    integer  :: status
    external error_handler
    !------------ Executable code --------------------------------

    n_spin = options_n_spin()
    exchange_fit = options_xcmode() == xcmode_exchange_fit
    spin_required = ( options_xcmode() ==  xcmode_model_density .or. &
                      options_xcmode() ==  xcmode_extended_mda ) .and. &
                    n_spin > 1

    if (.not.present(th)) then ! RETURN INSIDE!
       ! keep modus
       if (mode == recover_eigenvec) then
          ! save coeff_charge and coeff_spin
          if (.not.charge_and_spin_kept) then
             allocate(coeff_charge_kept(n_fit%n_ch),stat=status)
             if (status /= 0) call error_handler &
                  ("FIT_COEFF_STORE: allocation of coeff_charge_kept failed")
          endif
          coeff_charge_kept = coeff_charge
          if (spin_required) then
             if (.not.charge_and_spin_kept) then
                allocate(coeff_spin_kept(n_fit%n_ch),stat=status)
                if (status /= 0) call error_handler &
                     ("FIT_COEFF_STORE: allocation of coeff_spin_kept failed")
             endif
             coeff_spin_kept = coeff_spin
          endif
          charge_and_spin_kept = .true.
          if (exchange_fit) then
             ! save coeff_xc and coeff_en_xc
             if (.not.xc_and_en_xc_kept) then
                allocate(coeff_xc_kept(n_fit%n_xc,n_spin),stat=status)
                if (status /= 0) call error_handler &
                     ("FIT_COEFF_STORE: allocation of coeff_xc_kept failed")
                allocate(coeff_en_xc_kept(n_fit%n_xc),stat=status)
                if (status /= 0) call error_handler &
                     ("FIT_COEFF_STORE: allocation of coeff_en_xc_kept failed")
             endif
             coeff_xc_kept    = coeff_xc
             coeff_en_xc_kept = coeff_en_xc
             xc_and_en_xc_kept = .true.
          endif
       endif
       ! ==== RETURN POINT ====
       return
    endif

    if (output_main_scf) call write_to_output_units &
         ("FIT_COEFF_STORE: saving fit coefficients")

    if (present(mode)) then
       reset = mode == recover_eigenvec
    else
       reset = .false.
    endif
    if (output_data_saved) then
       write(output_unit,'(/ a     )')'Stored fit coefficients :'
    endif

    if (reset) then
       charge_and_spin_kept = .false.
       call readwriteblocked_write(coeff_charge_kept,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'coeff_charge'
          write(output_unit,'(4es20.13)')coeff_charge_kept
       endif
       deallocate(coeff_charge_kept,stat=status)
       if (status /= 0) call error_handler &
            ("FIT_COEFF_STORE: dallocation of coeff_charge_kept failed")
    else

       if(whatis(op_BackTrafo).eq.2)then
          print *,'w2 COE=',sum(coeff_charge),maxval(abs(coeff_charge))

          WARN('OFFSET coeff_charge by coeff_charge_veff')
          coeff_charge = coeff_charge + coeff_charge_veff

          print *,'w2 COE=',sum(coeff_charge),maxval(abs(coeff_charge))
       endif

       call readwriteblocked_write(coeff_charge,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'coeff_charge'
          write(output_unit,'(4es20.13)')coeff_charge
       endif
    endif

    if (.not.spin_required) then
       call readwriteblocked_write(no,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'coeff_spin exist ?'
          write(output_unit,'(4es20.13)')no
       endif
    else
       call readwriteblocked_write(yes,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'coeff_spin exist ?'
          write(output_unit,'(4es20.13)')yes
       endif
       if (reset) then
          call readwriteblocked_write(coeff_spin_kept,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'coeff_spin'
             write(output_unit,'(4es20.13)')coeff_spin_kept
          endif
          deallocate(coeff_spin_kept,stat=status)
          if (status /= 0) call error_handler &
               ("FIT_COEFF_STORE: dallocation of coeff_spin_kept failed")
       else
          call readwriteblocked_write(coeff_spin,th)
          if (output_data_saved) then
             write(output_unit,'(  a     )')'coeff_spin'
             write(output_unit,'(4es20.13)')coeff_spin
          endif
       endif
    endif

    if (.not.exchange_fit) then
       call readwriteblocked_write(no,th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'exchange fit used ?'
          write(output_unit,'(4es20.13)')no
       endif
    else
       call readwriteblocked_write(yes,th)
       call readwriteblocked_write((/real(n_spin,r8_kind), &
                                     real(n_fit%n_xc,r8_kind)/),th)
       if (output_data_saved) then
          write(output_unit,'(  a     )')'exchange fit used ?'
          write(output_unit,'(4es20.13)')yes
          write(output_unit,'(  a     )')'n_spin, n_xc'
          write(output_unit,'(4es20.13)')(/real(n_spin,r8_kind), &
                                           real(n_fit%n_xc,r8_kind)/)
       endif
       if (reset) then
          do s=1,n_spin
             call readwriteblocked_write(coeff_xc_kept(:,s),th)
          enddo
          call readwriteblocked_write(coeff_en_xc_kept,th)
          if (output_data_saved) then
             do s=1,n_spin
                write(output_unit,'( a,i1,a )')'coeff_xc(:,',s,')'
                write(output_unit,'(4es20.13)')coeff_xc_kept(:,s)
             enddo
             write(output_unit,'(  a     )')'coeff_en_xc'
             write(output_unit,'(4es20.13)')coeff_en_xc_kept
          endif
          xc_and_en_xc_kept = .false.
          deallocate(coeff_xc_kept,stat=status)
          if (status /= 0) call error_handler &
               ("FIT_COEFF_STORE: deallocation of coeff_xc_kept failed")
          deallocate(coeff_en_xc_kept,stat=status)
          if (status /= 0) call error_handler &
               ("FIT_COEFF_STORE: deallocation of coeff_en_xc_kept failed")
       else
          do s=1,n_spin
             call readwriteblocked_write(coeff_xc(:,s),th)
          enddo
          call readwriteblocked_write(coeff_en_xc,th)
          if (output_data_saved) then
             do s=1,n_spin
                write(output_unit,'( a,i1,a )')'coeff_xc(:,',s,')'
                write(output_unit,'(4es20.13)')coeff_xc(:,s)
             enddo
             write(output_unit,'(  a     )')'coeff_en_xc'
             write(output_unit,'(4es20.13)')coeff_en_xc
          endif
       endif
    endif
  end subroutine fit_coeff_store

!*************************************************************
! record  1: coeff_charge
! record  2: spin_required
! record  3: coeff_spin        [if spin_required]
! record  4: exchange_fit
! record  5: n_spin, n_xc      [if exchange_fit]
! record  6: coeff_xc(spin=1)  [if exchange_fit]
! record  7: coeff_xc(spin=2)  [if exchange_fit & n_spin > 1]
! record  8: coeff_en_xc       [if exchange_fit]
!*************************************************************

  subroutine fit_coeff_recover(th)
    ! Purpose: recovers the fit coefficients from a readwriteblocked
    !          file in case the input switch 'read_scfstate' is set.
    !
    ! routine called by: main_scf
    !** End of interface *****************************************
    !------------ Modules ----------------------------------------
    use options_module
    use readwriteblocked_module
    use iounitadmin_module, only: write_to_output_units, output_unit
    use output_module     , only: output_main_scf, output_data_read
    use spin_orbit_module , only: whatis, op_BackTrafo
    !------------ Declaration of formal paramaters ---------------
    type(readwriteblocked_tapehandle), intent(inout) :: th
    !------------ Declaration of local variables  ----------------
    allocatable           :: buffer
    integer(kind=i4_kind) :: n_spin,spin_stored,n_xc_stored
    real(kind=r8_kind)    :: spin_nxc(2), spin_req(1), exch_fit(1), &
                             buffer(:), zero = 0.0_r8_kind, half = 0.5_r8_kind
    logical               :: spin_required, exchange_fit
    integer               :: alloc_stat

    ! must correspond to arrays "yes" and "no" from fit_coeff_store:
    real(r8_kind), parameter :: yes = 1.0_r8_kind, no = 0.0_r8_kind
    logical               :: valid_saved_scfstate
    !------------ Executable code --------------------------------

    n_spin = options_n_spin()
    exchange_fit = options_xcmode() == xcmode_exchange_fit
    spin_required = ( options_xcmode() == xcmode_model_density .or. &
                      options_xcmode() == xcmode_extended_mda ) .and. &
                    n_spin > 1

    if (output_main_scf) call write_to_output_units &
         ("FIT_COEFF_RECOVER: reading fit coefficients")

    if ( .not. allocated(coeff_charge) ) then
        !
        ! FIXME: When called from properties_module, coeff_charge
        ! is not allocated:
        !
        WARN('not allocated(coeff_charge)')
        allocate(coeff_charge(fit_coeff_n_ch()))
    endif
    ASSERT(allocated(coeff_charge))
    call readwriteblocked_read(coeff_charge,th)

    if(whatis(op_BackTrafo).eq.2)then
       WARN('STORE coeff_charge from saved_scfstate')
       coeff_charge_veff = coeff_charge
       WARN('ZERO  coeff_charge')
       coeff_charge      = 0.0_r8_kind
    endif

    ch_initialized = .false.
    call readwriteblocked_read(spin_req,th)

    ! check sanity, maybe this is a different saved_scfstate?
    valid_saved_scfstate = spin_req(1) == yes .or. spin_req(1) == no
    ASSERT(valid_saved_scfstate)

    if (output_data_read) then
       write(output_unit,'(/ a     )')'Recovered fit coefficients :'
       write(output_unit,'(  a     )')'coeff_charge'
       write(output_unit,'(4es20.13)')coeff_charge
       write(output_unit,'(  a     )')'coeff_spin exist ?'
       write(output_unit,'(4es20.13)')spin_req(1)
    endif

    if (spin_req(1) == yes ) then
       if (spin_required) then
          call readwriteblocked_read(coeff_spin,th)
          sp_initialized = .false.
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_spin'
             write(output_unit,'(4es20.13)')coeff_spin
          endif
       else
          if (output_main_scf) call write_to_output_units &
               ("FIT_COEFF_RECOVER: coeff_spin ignored")
          call readwriteblocked_skipread(n_fit%n_ch,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_spin skipped'
          endif
       endif
    elseif (spin_required) then
       if (output_main_scf) call write_to_output_units &
            ("FIT_COEFF_RECOVER: coeff_spin used as initialized")
       sp_initialized = .true.
    endif
   
    call readwriteblocked_read(exch_fit,th)

    ! check sanity, maybe this is a different saved_scfstate?
    valid_saved_scfstate = exch_fit(1) == yes .or. exch_fit(1) == no
    ASSERT(valid_saved_scfstate)

    if (output_data_read) then
       write(output_unit,'(  a     )')'exchange fit used ?'
       write(output_unit,'(4es20.13)')exch_fit(1)
    endif
    if (exch_fit(1) == yes) then
       call readwriteblocked_read(spin_nxc,th)

       spin_stored = int(spin_nxc(1),i4_kind)
       n_xc_stored = int(spin_nxc(2),i4_kind)

       ! check sanity, maybe this is a different saved_scfstate?
       valid_saved_scfstate = &
               spin_nxc(1) - spin_stored == 0 .and. & ! integer
               spin_nxc(2) - n_xc_stored == 0 .and. & ! integer
               (spin_stored == 1 .or. spin_stored == 2) .and. & ! n_spin
               n_xc_stored > 0 ! n_xc
       ASSERT(valid_saved_scfstate)

       if (output_data_read) then
          write(output_unit,'(  a     )')'n_spin, n_xc'
          write(output_unit,'(4es20.13)')spin_nxc(1:2)
       endif
       if (exchange_fit) then
          call readwriteblocked_read(coeff_xc(:,1),th)
          xc_initialized = .false.
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_xc(:,1)'
             write(output_unit,'(4es20.13)')coeff_xc(:,1)
          endif
       else
          if (output_main_scf) call write_to_output_units &
               ("FIT_COEFF_RECOVER: coeff_xc ignored")
          call readwriteblocked_skipread(n_xc_stored,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_xc(:,1) skipped'
          endif
       endif
    else
       spin_stored = 0
       n_xc_stored = 0
       if (exchange_fit) then
          if (output_main_scf) call write_to_output_units &
               ("FIT_COEFF_RECOVER: coeff_xc used as initialized")
          xc_initialized = .true.
       endif
    endif
    if (exch_fit(1) == yes .and. spin_stored > 1) then
       if (exchange_fit .and. n_spin > 1) then
          call readwriteblocked_read(coeff_xc(:,2),th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_xc(:,2)'
             write(output_unit,'(4es20.13)')coeff_xc(:,2)
          endif
       elseif (exchange_fit) then
          if (output_main_scf) then
             call write_to_output_units("FIT_COEFF_RECOVER: "// &
                  "trying to convert spin-polarized data from")
             call write_to_output_units("                   "// &
                  "tape into the spin-restricted data required")
          endif
          allocate(buffer(n_fit%n_xc),STAT=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("FIT_COEFF_RECOVER: buffer allocation failed!")
          call readwriteblocked_read(buffer,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_xc(:,2)'
             write(output_unit,'(4es20.13)')buffer
          endif
          coeff_xc(:,1) = ( coeff_xc(:,1) + buffer ) * half
          deallocate(buffer,STAT=alloc_stat)
          if (alloc_stat /= 0) call error_handler &
               ("FIT_COEFF_RECOVER: buffer deallocation failed!")
       else
          call readwriteblocked_skipread(n_xc_stored,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_xc(:,2) skipped'
          endif
       endif
    elseif (exchange_fit .and. n_spin > 1) then
       if (exch_fit(1) /= zero) then
          if (output_main_scf) call write_to_output_units &
               ("FIT_COEFF_RECOVER: coeff_xc(spin) set to zero")
          coeff_xc(:,2) = coeff_xc(:,1)
          xc_initialized = .true. ! because the spin polarization is missing
       else
         ! used as initialized
       endif
    endif

    if (exch_fit(1) == yes) then
       if (exchange_fit) then
          call readwriteblocked_read(coeff_en_xc,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_en_xc'
             write(output_unit,'(4es20.13)')coeff_en_xc
          endif
       else
          call readwriteblocked_skipread(n_xc_stored,th)
          if (output_data_read) then
             write(output_unit,'(  a     )')'coeff_en_xc skipped'
          endif
       endif
    elseif (exchange_fit) then
       ! used as initialized
    endif
    
  end subroutine fit_coeff_recover
  !*************************************************************

#ifdef FPP_FIT_COEFF_INIT /* alternative initial guess */

  subroutine fit_coeff_initialize()
    !
    ! Purpose: if s- and r2-bases for charge fit contain non-empty
    !          array %funcitons(:, :) then %funcitons(:, 1) is used
    !          to initialize the charge fitting coefficients.
    !          The rest is set oto zero
    !
    ! Subroutine called by: prescf
    !
    use unique_atom_module, only: unique_atoms
    use orbitalprojection_module, only: orbitalprojection_ch
    use options_module
    implicit none
    ! *** end of interface ***

    integer(i4_kind)  :: ua, start, n_bas
    logical           :: spin_coeff_required

    call fit_coeff_allocate()

    ! most elements are zero ...
    coeff_charge = 0.0

    !
    ! ... but entries for s- and r2-functions may be found in the fit basis:
    !

    ! loop over uniques
    do ua = 1, size(unique_atoms)

       ! This is the r2-type.
       if ( size(unique_atoms(ua)%r2_ch%functions) > 0 ) then
          start = orbitalprojection_ch(0, ua)
          n_bas = size(unique_atoms(ua)%r2_ch%functions, 1)
          coeff_charge(start:start+n_bas-1) = unique_atoms(ua)%r2_ch%functions(:, 1)
       endif

       ! the same procedure for s-type...
       if ( size(unique_atoms(ua)%l_ch(0)%functions) > 0 ) then
          start = orbitalprojection_ch(-1, ua)
          n_bas = size(unique_atoms(ua)%l_ch(0)%functions, 1)
          coeff_charge(start:start+n_bas-1) = unique_atoms(ua)%l_ch(0)%functions(:, 1)
       endif
    enddo
    ch_initialized = .true.

    spin_coeff_required = ( options_xcmode() == xcmode_model_density .or. &
                            options_xcmode() == xcmode_extended_mda  ) &
                          .and. options_n_spin() > 1

    if (spin_coeff_required) then
       ! FIXME: why?
       coeff_spin = coeff_charge
       sp_initialized = .true.
    endif

    if (options_xcmode() == xcmode_exchange_fit) then
       coeff_xc = 0.0
       coeff_en_xc = 0.0
       xc_initialized = .true.
    endif

    if (options_xcmode() == xcmode_extended_mda) then
       coeff_proj = 0.0
       pv_initialized = .true.
    endif
  end subroutine fit_coeff_initialize

#else /* legacy LCGTO-inspired version */

  subroutine fit_coeff_initialize()
    ! Purpose: initializes the fit coefficients as it is done 
    !          in the old lcgto:
    !          Only those fit coefficients where the 
    !          corresponding function
    !          is either of s- or of r2-type
    !          are set to one. Of course we have to look for
    !          s- or r2-type functions among the globally
    !          contracted functions, too.
    !          All others are set to zero.
    ! Subroutine called by: prescf
    ! Author: FN
    ! Date: 8/96
    !** End of interface *****************************************
    !------------ Modules used ----------------------------------
    use unique_atom_module, only : unique_atoms,&
         N_unique_atoms
    use orbitalprojection_module, only: &
         orbitalprojection_ch, &
         orbitalprojection_xc, &
         orbitalprojection_globcontr_ch, &
         orbitalprojection_globcontr_xc
    use options_module
    !------------ Declaration of local variables -----------------
    integer(kind=i4_kind)  :: na,i_l,l_help,start,alloc_stat,i,&
         n_spin,i_spin,counter,i_glob,n_contributing,i_cont
    logical                :: spin_coeff_required
    logical,allocatable    :: mask_ch(:),mask_xc(:)
    !------------ Executable code --------------------------------
    n_spin = options_n_spin()
    spin_coeff_required = ( options_xcmode() == xcmode_model_density .or. &
                            options_xcmode() == xcmode_extended_mda  ) &
                          .and. n_spin > 1

    call fit_coeff_allocate()

    coeff_charge=zero
    allocate(mask_ch(n_fit%n_ch),STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("INITIALIZE_FITCOEFFS: allocation (2) failed")
    mask_ch=.false.

    if (options_xcmode() == xcmode_exchange_fit) then
       coeff_xc=zero
       coeff_en_xc=zero
       allocate(mask_xc(n_fit%n_ch),STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("INITIALIZE_FITCOEFFS: allocation (4) failed")
       mask_xc=.false.
    endif

    ! loop over uniques
    uniques: do na=1,N_unique_atoms

       ! This is the r2-type.
       i_l=-1
       l_help=0
       start = orbitalprojection_ch(l_help,na)       
       do i=start,start+unique_atoms(na)%r2_ch%N_uncontracted_fcts+&
            unique_atoms(na)%r2_ch%N_contracted_fcts-1
          mask_ch(i) = .true.
       enddo
       if (options_xcmode() == xcmode_exchange_fit) then
          start=orbitalprojection_xc(l_help,na)
          do i=start,start+unique_atoms(na)%r2_xc%N_uncontracted_fcts+&
               unique_atoms(na)%r2_xc%N_contracted_fcts-1  
             mask_xc(i)=.true.
          enddo
       endif

       ! the same procedure for s-type...
       i_l=0
       l_help=-1
       start=orbitalprojection_ch(l_help,na)
       do i=start,start+unique_atoms(na)%l_ch(i_l)%N_uncontracted_fcts+&
            unique_atoms(na)%l_ch(i_l)%N_contracted_fcts-1
          mask_ch(i)=.true.
       enddo
       if (options_xcmode() == xcmode_exchange_fit) then
          start=orbitalprojection_xc(l_help,na)
          do i=start,start+unique_atoms(na)%l_xc(i_l)%N_uncontracted_fcts+&
               unique_atoms(na)%l_xc(i_l)%N_contracted_fcts-1
             mask_xc(i)=.true.
          enddo
       endif

       ! ... and now the global contractions
        if (unique_atoms(na)%N_glob_cons_ch.gt.0) then
          counter = orbitalprojection_globcontr_ch(na)
          do i_glob=1,unique_atoms(na)%N_glob_cons_ch
             n_contributing = &
                  unique_atoms(na)%glob_con_ch(i_glob)%N_contributing_fcts 
             contributing_ch: do i_cont=1,n_contributing
                i_l = unique_atoms(na)%glob_con_ch(i_glob)%l(i_cont)
                if (i_l.eq.0.or.i_l.eq.-1) then
                   mask_ch(counter)=.true.
                   exit contributing_ch
                endif
             enddo contributing_ch
             counter=counter+1
          enddo
       endif
       if (options_xcmode() == xcmode_exchange_fit .and. &
           unique_atoms(na)%N_glob_cons_xc.gt.0) then
          counter = orbitalprojection_globcontr_xc(na) 
          do i_glob=1,unique_atoms(na)%N_glob_cons_xc
             n_contributing = &
                  unique_atoms(na)%glob_con_xc(i_glob)%N_contributing_fcts 
             contributing_xc: do i_cont=1,n_contributing
                i_l = unique_atoms(na)%glob_con_xc(i_glob)%l(i_cont)
                if (i_l.eq.0.or.i_l.eq.-1) then
                   mask_xc(counter)=.true.
                   exit contributing_xc
                endif
             enddo contributing_xc
             counter=counter+1
          enddo
       endif
                
    enddo uniques

    ! now use the mask to set the fitcoefficients.
    where (mask_ch) coeff_charge = one
    ch_initialized = .true.
    if (spin_coeff_required) then
       coeff_spin = coeff_charge
       sp_initialized = .true.
    endif
    if (options_xcmode() == xcmode_exchange_fit) then
       do i_spin=1,n_spin
          where (mask_xc) coeff_xc(:,i_spin) = one
       enddo
       xc_initialized = .true.
       coeff_en_xc = coeff_xc(:,1)
    endif
    if (options_xcmode() == xcmode_extended_mda) then
       pv_initialized = .true.
       coeff_proj = zero
    endif

    ! deallocation of masks
    deallocate(mask_ch,STAT=alloc_stat)
    if (alloc_stat.ne.0) call error_handler &
         ("INITIALIZE_FITCOEFFS: deallocation (1) failed")
    if (options_xcmode() == xcmode_exchange_fit) then
       deallocate(mask_xc,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("INITIALIZE_FITCOEFFS: deallocation (1) failed")
    endif

  end subroutine fit_coeff_initialize

#endif /* ifdef FPP_FIT_COEFF_INIT */

  !*************************************************************
  subroutine fit_coeff_normalize(spin_coeff)
    !
    !  Purpose: normalize the chargefit coefficients to the
    !           total number of electrons
    !
    ! NOTE: normalization of the density fitting coefficients is
    !       disabled. In a situation when an initial charge density
    !       is provided for C atom in CH4 but not for H-atoms
    !       it is equally wrong to scale the density on H atom up
    !       by ~66% to match the total number of electron as not
    !       to scale at all.
    !
    !       In a regular case, charge fit procedure is supposed to output
    !       normalized fitting coefficients.
    !
    !  Subroutine called by: prescf
    !
    !------------ Modules used -----------------------------------
    use occupation_module, only: get_n_elec, get_spin_diff, get_min_spin_diff
    use options_module
    implicit none
    logical,optional,intent(in)      :: spin_coeff
    !  if present and true the spin density coefficients are normalized too
    ! *** end of interface ***

    real(r8_kind) :: summ, n_elec, factor, spin_diff

    summ = dot_product(coeff_charge, charge_norm)

    call get_n_elec(n_elec)

    if ( summ == 0.0 ) then
       WARN('fitted charge is zero')
    else
       ! summ might be zero if Veff is used
       factor = n_elec / summ

! DONT WARN('skip renormalization of charge_coeff')
! DONT coeff_charge = coeff_charge * factor
    endif

    DPRINT   'fit_coeff_normalize: number of electrons', n_elec
    DPRINT   'fit_coeff_normalize:       fitted charge', summ
    DPRINT   "fit_coeff_normalize:                diff", summ - n_elec

    if (present(spin_coeff) )then
       if ( spin_coeff .and. options_n_spin() > 1 .and. &
            ( options_xcmode() == xcmode_model_density .or. &
              options_xcmode() == xcmode_extended_mda ) ) then

          call get_spin_diff(spin_diff)

          ! to avoid devision by zero in case spin_diff /= 0 but summ == 0
          if (spin_diff == zero) then
             coeff_spin = zero
          else
             summ = dot_product(coeff_spin, charge_norm)

             if (summ >= get_min_spin_diff()) then
                factor = spin_diff / summ

                coeff_spin = coeff_spin * factor
             endif
          endif
       endif
    endif
  end subroutine fit_coeff_normalize


  subroutine fit_coeff_charge_norm()
    !  Purpose: check the charge-fit norm 
    !  Subroutine called by: chargefit
    !** End of interface *****************************************
    use occupation_module, only: get_n_elec, get_spin_diff
    use options_module
    implicit none
    real(kind=r8_kind)           :: summ, n_elec, spin_diff
    !------------ Executable code -----------------------------------

    summ = dot_product(coeff_charge, charge_norm)

    call get_n_elec(n_elec)

    if ( abs(summ - n_elec) > 1.0E-7 ) then
        WARN('fit quality or large system?')
    endif

    if ( ( options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda ) .and. &
         options_n_spin() > 1 ) then

       summ = dot_product(coeff_spin, charge_norm)

       call get_spin_diff(spin_diff)

       if( size(coeff_charge) > 0 ) then
           ASSERT(abs(spin_diff-summ)<1.0E-7)
       endif
    endif
  end subroutine fit_coeff_charge_norm
  !*************************************************************


  !*************************************************************
  subroutine dimensions_calc(n_fit,map)
    ! Purpose: calculate the dimension of the sym.adapted,
    !          contracted fitbases. This was formerly done in the
    !          unique_atom_module. 
    !
    ! Routine called by: initialize_master
    !** End of interface **************************************
    !------------ Modules used --------------------------------
    use unique_atom_module, only: unique_atoms, &
         unique_atom_type,unique_atom_basis_type,&
         unique_atom_partner_type,N_unique_atoms
#ifdef WITH_CORE_DENS
     use unique_atom_module, only: &
         unique_atom_atomic_dens_type,pseudopot_present, &
         core_density_setup
#endif
    use symmetry_data_module, only: get_totalsymmetric_irrep
    use options_module, only: options_xcmode, &
                              xcmode_exchange_fit, xcmode_extended_mda
    type(fit), intent(inout) :: n_fit
    type(ff_info), intent(out), optional :: map(:) ! (N_FF)
    !------------ Declaration of local variables --------------
    integer(kind=i4_kind)  :: dim,i_ir,i_ua,i_l
    type(unique_atom_type),pointer :: ua
    type(unique_atom_basis_type),pointer :: uab
#ifdef WITH_CORE_DENS
    type(unique_atom_atomic_dens_type), pointer :: uac
#endif
    type(unique_atom_partner_type),pointer :: uap
    integer(i4_kind) :: dims(-1:2) ! dimenstons: s, r2, L>0, glob_contr
    integer(i4_kind) :: counter
    integer(i4_kind) :: nff,ff,nffu,nffc,nffi,ie,ii
    !----------- Executable code ------------------------------

    ! charge fitfunctions
    i_ir = get_totalsymmetric_irrep()
    dims    = 0
    counter = 0
    do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       do i_l = -1, ua%lmax_ch ! s, r2, p, ...
          select case(i_l)
          case (-1) ! s-
             uab => ua%l_ch(0)
             uap => ua%symadapt_partner(i_ir,0)
          case (0)  ! r2-
             uab => ua%r2_ch
             uap => ua%symadapt_partner(i_ir,0)
          case (1:) ! p-, d-, ...
             uab => ua%l_ch(i_l)
             uap => ua%symadapt_partner(i_ir,i_l)
          case default
             ABORT('no such case')
          end select
          nffu = uab%N_uncontracted_fcts
          nffc = uab%N_contracted_fcts
          nffi = uap%N_independent_fcts
          nff  = (nffu + nffc) * nffi
          ! -1=s, 0=r2, 1=L>0:
          dims(sgn(i_l)) = dims(sgn(i_l)) + nff

          ! store the mapping ff -> (ua,L) if requested
          if( present(map) )then
             ASSERT(counter+nff<=size(map))
             map(counter+1:counter+nff)%UA   = i_ua
             map(counter+1:counter+nff)%Z    = ua%Z
             map(counter+1:counter+nff)%L    = i_l
             map(counter+1:counter+nff)%next_ual = counter+nff+1
             ! ASSERT(nffc==0)
             ! ASSERT(nffu==size(uab%exponents))
             !
             ! This code saves an exponent for each symmetry adapted
             ! fit funciton that only used in the relativistic fit
             ! trafo code for (optional) screening:
             !
             ff = 0
             do ie = 1, nffu ! uncontracted exponents first
                do ii = 1, nffi
                   ff = ff + 1
                   map(counter + ff)%aexp = uab%exponents(ie)
                enddo
             enddo
             do ie = 1, nffc ! contracted functions second
                do ii = 1, nffi
                   ff = ff + 1
                   ! FIXME: use a different measure of extent
                   map(counter + ff)%aexp = minval(uab%exponents)
                enddo
             enddo
             ASSERT(ff==nff)
          endif
          counter = counter + nff

       enddo
       if ( ua%N_glob_cons_ch .gt. 0 )then
          dims(2) = dims(2)+ ua%N_glob_cons_ch
          ASSERT(.not.present(map))
       endif
    enddo
    n_fit%n_ch  = sum(dims)
    n_fit%n_s   = dims(-1)
    n_fit%n_r2  = dims( 0)

    dim = sum(dims)
    ASSERT(counter==dim)

#ifdef WITH_CORE_DENS
    n_fit%n_cd = 0
    if (operations_core_density .or. (pseudopot_present.and.core_density_setup)) then
       dim = 0
       do i_ua = 1,N_unique_atoms
          ua => unique_atoms(i_ua)
          do i_l = -1, 0
             select case (i_l)
             case (-1) ! s-type
                uac => ua%s_core
             case (0)  ! r^2-type
                uac => ua%r2_core
             end select
             ! uap%N_independent_fcts is always 1 here, thus
             ASSERT(uac%N_exponents>=0)
             dim = dim + uac%N_exponents     
          enddo
       enddo
       n_fit%n_cd = dim
       !...  duplicate code removed ...
    end if
#endif
 
    if ( options_xcmode() == xcmode_exchange_fit .or. &
         options_xcmode() == xcmode_extended_mda ) then
    ! exchange fitfunctions
    dim = 0
    do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       do i_l = 0, ua%lmax_xc
          uap => ua%symadapt_partner(i_ir,i_l)
          uab => ua%l_xc(i_l)
          dim = dim + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
               * uap%N_independent_fcts
       enddo
       uap => ua%symadapt_partner(i_ir,0)
       uab => ua%r2_xc
       dim = dim + (uab%N_uncontracted_fcts + uab%N_contracted_fcts) &
            * uap%N_independent_fcts
       if ( ua%N_glob_cons_xc .gt. 0 ) dim = dim + ua%N_glob_cons_xc
    enddo
    n_fit%n_xc = dim
    else
       n_fit%n_xc = 0
    endif

  contains

    integer(i4_kind) function sgn(x)
      implicit none
      integer(i4_kind), intent(in) :: x
      ! *** end of interface ***

      select case(x)
      case( 0)
         sgn =  0
      case(-1)
         sgn = -1
      case(1:)
         sgn = +1
      case default
         sgn = -999
         ABORT('shlould not be')
      end select
    end function sgn

  end subroutine dimensions_calc
  !*************************************************************

  subroutine fit_coeff_dimensions_calc()
    ! Wraper for the above sub.
    ! Purpose: calculate the dimension of the sym.adapted,
    !          contracted fitbases. This was formerly done in the
    !          unique_atom_module. 
    !
    ! Routine called by: initialize_master
    implicit none
    !** End of interface **************************************

    integer(i4_kind) :: memstat
    integer(i4_kind) :: n_s, n_r2 ! for debug only

    ! compute dimensions and strore them in global
    ! struct:
    call dimensions_calc(n_fit) ! "dry" run, compute dimensions

    DPRINT 'fit_coeff_dimensions_calc: allocate(fit_coeff_ff_map)'
    allocate(fit_coeff_ff_map(n_fit%n_ch),stat=memstat)
    ASSERT(memstat==0)

    ! fill in the map
    call dimensions_calc(n_fit,fit_coeff_ff_map)

    ! for debug only:
    n_s  = count(fit_coeff_ff_map(:)%L==-1)
    n_r2 = count(fit_coeff_ff_map(:)%L==0 )
    ASSERT(n_s==n_fit%n_s)
    ASSERT(n_r2==n_fit%n_r2)
  end subroutine fit_coeff_dimensions_calc

#ifdef WITH_CORE_DENS
    subroutine calc_fitfct_map(fitfct_map)
    ! purpose: establishes the mapping array between the fitting
    !          functions and the core density fitting functions
    ! only called if operations_core_density is true
    !** End of interface **************************************
    !------------ Modules used --------------------------------
    use unique_atom_module, only: unique_atoms, N_unique_atoms, &
         unique_atom_type,unique_atom_atomic_dens_type
    use orbitalprojection_module, only: orbitalprojection_ch, &
                                        orbitalprojection_globcontr_ch ! for testing
    implicit none
    logical, intent(out) :: fitfct_map(:) ! (n_ch)
    !------------ Declaration of local variables --------------
    integer(kind=i4_kind)  :: i_ua,i_bas,i_cexp,i_fexp
    type(unique_atom_type),pointer :: ua
    type(unique_atom_atomic_dens_type), pointer :: uac
    !----------- Executable code ------------------------------

    fitfct_map = .false.
    do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       if (ua%zc /= 0.0_r8_kind) then
          do i_bas = -1, 0 ! just s- and r^2-type fitting functions
             select case(i_bas)
             case (-1) ! s-type
                uac => ua%s_core
             case (0) ! r2-type
                uac => ua%r2_core
             end select
             do i_cexp = 1, uac%N_exponents
                i_fexp = uac%linked_fitfcts(i_cexp) &
                       + orbitalprojection_ch(i_bas,i_ua) - 1
                fitfct_map(i_fexp) = .true.
             enddo
          enddo
!:TST[
    write(1,*)'orbitalprojection_ch : [s,r2,p,...]'
    write(1,*) orbitalprojection_ch(:,i_ua)
    write(1,*)'orbitalprojection_globcontr_ch :'
    write(1,*) orbitalprojection_globcontr_ch(i_ua)
!:TST]
       endif
    enddo
!:TST[
    write(1,*)'fitfct_map :'
    write(1,'(40L2)')fitfct_map
!:TST]
  end subroutine calc_fitfct_map

    subroutine write_coeff_core(loop)
    ! purpose: write the final core density charge fit functions 
    !          to file (excluding atomic re-normalization)
    ! only called if operations_core_density is true
    !** End of interface **************************************
    !------------ Modules used --------------------------------
    use unique_atom_module, only: unique_atoms, N_unique_atoms, &
         unique_atom_type,unique_atom_atomic_dens_type
    use iounitadmin_module
    use options_module
    use filename_module, only: outfile
    !------------ Declaration of local variables --------------
    integer(kind=i4_kind),intent(in),optional  :: loop
    integer(kind=i4_kind)  :: i_ua,i_bas,i_cexp,i_sa,n_exps, &
                              io_u, count=0
    type(unique_atom_type),pointer :: ua
    type(unique_atom_atomic_dens_type), pointer :: uac
    !----------- Executable code ------------------------------
    100 format('  # unique atom ... ',a,'-type core density')
    101 format('  &UNIQUE_ATOM_CORE_DENSITY'/ &
               '    N_EXPONENTS = ',i5/ &
               '  /UNIQUE_ATOM_CORE_DENSITY')
    102 format('  # ',a)
    103 format(3ES25.15:"  %")

    io_u = get_iounit()
    open(io_u,status='unknown',form='formatted', &
         position='append', &
         file=trim(outfile('coeff_core_charge')))

    if(.not.present(loop)) then
       count=count+1
    endif
    write(io_u,*)' '
    if(present(loop)) then
       write(io_u,*)' +++++++++++++ Loop ',loop,' +++++++++++++++'
    else
       write(io_u,*)' +++++++++++++ Loop ',count,' +++++++++++++++'
    endif

    i_sa = 0
    do i_ua = 1,N_unique_atoms
       ua => unique_atoms(i_ua)
       if (ua%zc /= 0.0_r8_kind) then
          do i_bas = -1, 0 ! just s- and r^2-type fitting functions
             select case(i_bas)
             case (-1) ! s-type
                write(io_u,100) 's'
                uac => ua%s_core
             case (0) ! r2-type
                write(io_u,100) 'r2'
                uac => ua%r2_core
             end select

             n_exps = uac%N_exponents
             write(io_u,101) n_exps
             write(io_u,102) 'exponents'
             write(io_u,103) uac%exponents
             write(io_u,102) 'contraction (non-renormalized)'
             write(io_u,103) (coeff_core(i_sa+uac%linked_fitfcts(i_cexp)) &
                           ,i_cexp=1,n_exps)
             i_sa = i_sa + n_exps
          enddo
       endif
    enddo

    close(io_u)
    call return_iounit(io_u)
  end subroutine write_coeff_core
 !****************************************************************************
#endif

  !***************************************************************************
  subroutine copy_coeff_charge_old()
    ! purpose: copies coeff_charge into coeff_charge_old
    !          after checking if coeff_charge_old is
    !          allocated. If not, coeff_charge_old will
    !          be allocated by this routine.
    !** End of interface *************************************** 
    !------------ Modules used --------------------------------
    use options_module
    !------------ Declaration of local variables --------------
    integer(kind=i4_kind)    :: alloc_stat, n_spin
    external error_handler
    ! ----- executable code -----------------------------------
!   the old charge coefficients are now deallocated directly 
!   after the call of mixing_ch in chargefit, thus:
    n_spin = options_n_spin()
    allocate(coeff_charge_old(n_fit%n_ch),stat=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("COPY_COEFF_CHARGE_OLD: allocation of coeff_charge_old failed")
    coeff_charge_old = coeff_charge
    chold_initialized = ch_initialized
    if ( ( options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda ) .and. n_spin > 1 ) then
       allocate(coeff_spin_old(n_fit%n_ch),stat=alloc_stat)
       if(alloc_stat.ne.0) call error_handler &
            ("COPY_COEFF_CHARGE_OLD: allocation of coeff_spin_old failed")
       coeff_spin_old = coeff_spin
       spold_initialized = sp_initialized
    endif
  end subroutine copy_coeff_charge_old
  !*************************************************************
  

  !*************************************************************
  subroutine copy_coeff_xc_old()
    ! purpose: as in 'copy_coeff_charge_old', simply
    !          for the exchange fit coefficients
    !** End of interface *************************************** 
    use options_module, only: options_n_spin
    integer(kind=i4_kind)     :: alloc_stat,n_spin
    external error_handler
    ! ----- executable code -----------------------------------
!   the old exchange coefficients are now deallocated directly 
!   after the call of mixing_xc in build_xcfit, thus:
    n_spin = options_n_spin()
    allocate(coeff_xc_old(n_fit%n_xc,n_spin),stat=alloc_stat)
    if(alloc_stat.ne.0) call error_handler &
         ("COPY_COEFF_XC_OLD: allocation failed")
    coeff_xc_old = coeff_xc
    xcold_initialized = xc_initialized
  end subroutine copy_coeff_xc_old
  !*************************************************************


  !*************************************************************
  subroutine free_coeff_charge_old()
    ! purpose: deallocates coeff_charge_old 
    !** End of interface *************************************** 
    !------------ Modules used --------------------------------
    use options_module
    !------------ Declaration of local variables --------------
    integer(kind=i4_kind)    :: dealloc_stat
    external error_handler
    ! ----- executable code -----------------------------------
    deallocate(coeff_charge_old,stat=dealloc_stat)
    if(dealloc_stat.ne.0) call error_handler &
         ("FREE_COEFF_CHARGE_OLD: deallocation of coeff_charge_old failed")
    if ( ( options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda ) .and. &
         options_n_spin() > 1 ) then
       deallocate(coeff_spin_old,stat=dealloc_stat)
       if(dealloc_stat.ne.0) call error_handler &
            ("FREE_COEFF_CHARGE_OLD: deallocation of coeff_spin_old failed")
    endif
  end subroutine free_coeff_charge_old

  !*************************************************************
  subroutine free_coeff_xc_old()
    ! purpose: deallocates coeff_xc_old()
    !** End of interface *************************************** 
    use options_module, only: options_n_spin
    integer(kind=i4_kind)     :: dealloc_stat
    external error_handler
    ! ----- executable code -----------------------------------
    deallocate(coeff_xc_old,stat=dealloc_stat)
    if(dealloc_stat.ne.0) call error_handler &
         ("free_COEFF_XC_OLD: deallocation failed")
  end subroutine free_coeff_xc_old
  !*************************************************************


  !*************************************************************
  subroutine print_coeff_charge(loop)
    ! purpose print out the charge fitcoefficients at any
    !         point of the program to file 'coeff_charge'
    use iounitadmin_module
    use options_module
    use filename_module, only: outfile
    !------------ Declaration of formal parameters -------------
    integer(kind=i4_kind),intent(in),optional  :: loop
    !** End of interface *************************************** 
    integer(kind=i4_kind)   :: io_u, count=0, i

    io_u=get_iounit()
    open(io_u,status='unknown',form='formatted', &
         position='append', &
         file=trim(outfile('coeff_charge')))

    if(.not.present(loop)) then
       count=count+1
    endif
    write(io_u,*)' '
    if(present(loop)) then
       write(io_u,*)' +++++++++++++ Loop ',loop,' +++++++++++++++'
    else
       write(io_u,*)' +++++++++++++ Loop ',count,' +++++++++++++++'
    endif
    write(io_u,*)' '
    write(io_u,*)(coeff_charge(i), i=1,n_fit%n_ch)
    if ( ( options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda ) .and. &
         options_n_spin() > 1 ) then
       write(io_u,*)' '
       write(io_u,*)(coeff_spin(i), i=1,n_fit%n_ch)
    endif
    close(io_u)
    call return_iounit(io_u)
  end subroutine print_coeff_charge
  !*************************************************************

  subroutine fit_coeff_sndrcv (IMOD)
    !
    ! To be used in parallel context. Handles serial case ok.
    !
    use comm_module
    use msgtag_module, only: msgtag_fit_coeff_send
    implicit none
    integer (i4_kind), intent(in) :: IMOD
    ! *** end of interface ***

    if (.not. comm_parallel ()) return

    if (comm_i_am_master()) then
       call fit_coeff_send (IMOD)
    else
       call comm_save_recv (comm_master_host, msgtag_fit_coeff_send)
       ! unpack, actualy:
       call fit_coeff_receive()
    endif
  end subroutine fit_coeff_sndrcv

  !*************************************************************
  subroutine fit_coeff_send (IMOD)
    !
    ! Send fit_coeffs to  the slaves; the fitcoeffs are  needed on the
    ! slaves during the gradient run  and during the SCF SCF procedure
    ! as well if the model density approach is used.
    !
    use comm_module
    use msgtag_module
    use xpack, only: pck
    use options_module
    implicit none
    integer (i4_kind), intent (in), optional :: IMOD
    !** End of interface ***************************************

    integer(kind=i4_kind) :: n_spin
    logical               :: send_rho_coeff, send_pot_coeff, extended_mda
    integer(i4_kind) :: IMODE

    IMODE = ICHFIT ! send density fit, by default
    if (present (IMOD)) IMODE = IMOD

    send_rho_coeff = IAND (IMODE, ICHFIT) /= 0
    send_pot_coeff = IAND (IMODE, IXCFIT) /= 0

    extended_mda = options_xcmode() == xcmode_extended_mda

    call comm_init_send(comm_all_other_hosts,msgtag_fit_coeff_send) 
    call pck(n_fit%n_ch)
    call pck(send_rho_coeff)
    if (send_rho_coeff) then
       call pck(coeff_charge)
    endif
    if (options_xcmode() == xcmode_model_density .or. extended_mda) then
       n_spin = options_n_spin() 
       if (n_spin > 1 .and. send_rho_coeff) then
          call pck(coeff_spin)
       endif
       call pck(send_pot_coeff)
       if (send_pot_coeff) then
          call pck(coeff_xcmda)
          call pck(coeff_xcmda_ts)
          if (extended_mda) then
             call pck(n_fit%n_xc)
             call pck(coeff_xc)
          endif
       endif
    endif
    call comm_send()     
    
  end subroutine fit_coeff_send
  !*************************************************************

  !*************************************************************
  subroutine fit_coeff_receive()
    ! purpose: Receive fit_coeffs from  the master; the fitcoeffs are needed
    !          on the slaves during the gradient run and during the SCF
    !          SCF procedure as well if the model density approach is used.
    use comm_module
    use xpack, only: upck
    use options_module
    !------------ Declaration of formal parameters -------------
    !** End of interface ***************************************
    integer(kind=i4_kind) :: alloc_stat, n_spin
    logical               :: receive_rho_coeff, receive_pot_coeff, &
                             extended_mda

    extended_mda = options_xcmode() == xcmode_extended_mda

    call upck(n_fit%n_ch)
    call upck(receive_rho_coeff)
    if (receive_rho_coeff) then
       allocate(coeff_charge(n_fit%n_ch),stat=alloc_stat)
       if(alloc_stat/=0) call error_handler(&
            'fit_coeff_receive: allocation of coeff_charge failed')
       call upck(coeff_charge)
    endif
    if (options_xcmode() == xcmode_model_density .or. extended_mda) then
       n_spin = options_n_spin()
       if (n_spin > 1 .and. receive_rho_coeff) then
          allocate(coeff_spin(n_fit%n_ch),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'fit_coeff_receive: allocation of coeff_spin failed')
          call upck(coeff_spin)
       endif
       call upck(receive_pot_coeff)
       if (receive_pot_coeff) then
          allocate(coeff_xcmda(n_fit%n_ch,n_spin),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'fit_coeff_receive: allocation of coeff_xcmda failed')
          call upck(coeff_xcmda)
          allocate(coeff_xcmda_ts(n_fit%n_ch,n_spin),stat=alloc_stat)
          if(alloc_stat/=0) call error_handler(&
               'fit_coeff_receive: allocation of coeff_xcmda failed')
          call upck(coeff_xcmda_ts)
          if (extended_mda) then
             call upck(n_fit%n_xc)
             allocate(coeff_xc(n_fit%n_xc,n_spin),stat=alloc_stat)
             if(alloc_stat/=0) call error_handler(&
                  'fit_coeff_receive: allocation of coeff_xc failed')
             call upck(coeff_xc)
          endif
       endif
    endif
  end subroutine fit_coeff_receive

  subroutine fit_coeff_shutdown()
    !
    ! Purpose: deallocate all variables of this module
    !
    ! Executed on all workers in parallel context.
    !
    !------------ Declaration of formal parameters -------------
    use options_module
    !** End of interface ***************************************
!:UB[ modified
    integer(kind=i4_kind) :: alloc_stat, n_spin
    logical               :: exchange_fit, model_density, extended_mda
    external error_handler
    !------------ Executable code --------------------------------
    n_spin        = options_n_spin()
    exchange_fit  = options_xcmode() == xcmode_exchange_fit
    model_density = options_xcmode() == xcmode_model_density
    extended_mda  = options_xcmode() == xcmode_extended_mda
!:UB]

    if (allocated(coeff_charge)) then
       deallocate(coeff_charge,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_shutdown: deallocation of coeff_charge failed")

    if (allocated(coeff_charge_veff)) then
       deallocate(coeff_charge_veff,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    endif

#ifdef WITH_CORE_DENS
    if (operations_core_density) then
       deallocate(coeff_core,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_shutdown: deallocation (1c) failed ")
       deallocate(fitfct_map,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_shutdown: deallocation (1f) failed ")
    endif
#endif
    endif
    if (allocated(charge_norm)) then
       deallocate(charge_norm,STAT=alloc_stat)
       if (alloc_stat.ne.0) call error_handler &
            ("fit_coeff_shutdown: deallocation of charge_norm failed ")
    endif
    if (model_density .or. extended_mda) then
       if (n_spin > 1) then
          if (allocated(coeff_spin)) then
             deallocate(coeff_spin,STAT=alloc_stat)
             if (alloc_stat.ne.0) call error_handler &
                  ("fit_coeff_shutdown: deallocation of coeff_spin failed")
          endif
       endif
       if (allocated(coeff_xcmda)) then
          deallocate(coeff_xcmda,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("fit_coeff_shutdown: deallocation of coeff_xcmda failed")
       endif
       if (allocated(coeff_xcmda_ts)) then
          deallocate(coeff_xcmda_ts,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("fit_coeff_shutdown: deallocation of coeff_xcmda failed")
       endif
    endif
    if (extended_mda) then
       if (allocated(coeff_proj)) then
          deallocate(coeff_proj,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("fit_coeff_shutdown: deallocation of coeff_proj failed")
       endif
       if (allocated(coeff_deltarho)) then
          deallocate(coeff_deltarho,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("fit_coeff_shutdown: deallocation of coeff_deltarho failed")
       endif
    endif
    if (exchange_fit .or. extended_mda) then
       if (allocated(coeff_xc)) then
          deallocate(coeff_xc,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("fit_coeff_shutdown: deallocation of coeff_xc failed ")
       endif
    endif
    if (exchange_fit) then
       if (allocated(coeff_en_xc)) then
          deallocate(coeff_en_xc,STAT=alloc_stat)
          if (alloc_stat.ne.0) call error_handler &
               ("fit_coeff_shutdown: deallocation of coeff_en_xc failed ")
       endif
    endif
!:UB]

    if( allocated(fit_coeff_ff_map) )then
       DPRINT 'fit_coeff_shutdown: deallocate(fit_coeff_ff_map)'
       deallocate(fit_coeff_ff_map,STAT=alloc_stat)
       ASSERT(alloc_stat==0)
    else
       ABORT('map not allocated')
    endif

  end subroutine fit_coeff_shutdown

  subroutine fit_coeff_set_ch( coeff_charge_new )
    ! resets coeff_charge array. needed to make it private
    implicit none
    real(r8_kind), intent(in) :: coeff_charge_new(:)
    ! *** end of interface ***
    integer(i4_kind)          :: nnch, n_ch
    ! ************************
    ASSERT(allocated(coeff_charge))
    nnch = size(coeff_charge_new)
    n_ch = size(coeff_charge)
    ASSERT(nnch==n_ch)

    coeff_charge = coeff_charge_new

  end subroutine fit_coeff_set_ch
    
  !--------------- End of module ----------------------------------
end module fit_coeff_module

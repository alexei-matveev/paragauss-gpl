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
module  int_data_2cff_module
  !---------------------------------------------------------------
  !
  !  Purpose: This is the module containing the data used when
  !           calculating the 2 center fitfuncyion
  !           integrals to one specific quadrupel
  !           ( unique atom 1, l 1, unique atom 2, l 2 ).
  !           Which qudrupel is actually processed and the
  !           main loop variables and pointers to relevant data
  !           structures are to be contained here as well.
  !           Variables to hold primitive, contracted
  !           and symmetry adapted integarals are contained.
  !           Variables of primitive and contracted integrals
  !           only hold the values for one specific non-symmetry-
  !           equivalent connection vector.
  !
  !           Methods for allocating and freeing all necessary
  !           integals (int_data_2cff_setup and 
  !           int_data_2cff_shutdown) are included.
  !           int_data_2cff_setup should be called after
  !           quadrupel is set to present value.
  !
  !           To shorten list of use statements of calling
  !           routines and modules, this module is public.
  !
  !  Module called by: integral_calc_2cff, 
  !                    fitcontract_2c_module
  !                    int_send_2cff_module
  !
  !  References: Publisher Document: Concepts of Integral Part
  !
  !  Author: TB
  !  Date: 5/96
  !
  !----------------------------------------------------------------
  !== Interrupt of public interface of module =====================
  !----------------------------------------------------------------
  ! Modifications
  !----------------------------------------------------------------
  ! Modification (Please copy before editing)
  ! Author: FN
  ! Date:   8/96
  ! Description: added Variables containing the locally contracted
  !              and parts of the globally contracted 2-center
  !              fitintegrals
  ! Modification (Please copy before editing)
  ! Author:      Uwe Birkenheuer
  ! Date:        8/98
  ! Description: intermediate arrays for mixed overlap integrals introduced
  !              global control variable need_ch and need_xc replaced by
  !              need_ch_ch, need_xc_xc and new control variables need_xc_ch
  !              and need_ch_xc introduced.    
  !
  ! Modification (Please copy before editing)
  ! Author:      SB
  ! Date:        03/05
  ! Description: response case with all irreps
  !
  ! Modification (Please copy before editing)
  ! Author:      ...
  ! Date:        ...
  ! Description: ...
  !----------------------------------------------------------------
  !== Interrupt end of public interface of module =================
  !------------ all modules used are public !! ------------------
# include <def.h>
  use type_module
  use quadrupel_module
  use unique_atom_module
  use symmetry_data_module, only: get_totalsymmetric_irrep, symmetry_data_n_irreps
  use integralpar_module
  use fit_coeff_module
  use datatype, only: arrmat4
  !== Interrupt of public interface of module =====================
  implicit none
  save
  !== Interrupt end of public interface of module =================
  !------------ Declaration of constants and variables ----------
  type(quadrupel_type), public :: quadrupel
  ! quadrupel ( unique atom 1, l 1, unique atom 2, l 2 )
  ! that is processed right now
  type(unique_atom_symequivatoms_type), public, pointer :: symequivatoms
  ! different sym.equiv. atom pairs to be processed
  integer(kind=i4_kind), public  :: i_symequivatoms
  ! index of not-symmetry-equivalent connection vector
  ! that is processed right now
  type(unique_atom_type), pointer, public  :: ua1, ua2
  ! unique atoms processed right now
  type(unique_atom_basis_type), pointer, public :: &
       ua1_basis_ch, ua2_basis_ch, ua1_basis_xc, ua2_basis_xc
  ! basis of unique atoms processed right now

  integer(kind=i4_kind), public  :: &
       n_ch_exp1, n_ch_exp2, n_xc_exp1, n_xc_exp2
  ! number of exponents in ua1_basis and ua2_basis
  integer(kind=i4_kind), public  :: &
       n_ch_contr1, n_ch_contr2, n_xc_contr1, n_xc_contr2
  ! number of contractions in ua1_basis and ua2_basis
  integer(kind=i4_kind), public  :: &
       n_ch_uncontr1, n_ch_uncontr2, n_xc_uncontr1, n_xc_uncontr2
  ! number of uncontracted functions in ua1_basis and ua2_basis
  integer(kind=i4_kind), public  :: &
       n_ch_gc1, n_ch_gc2, n_xc_gc1, n_xc_gc2
  ! number of globally contracted functions in ua1_basis and ua2_basis
  integer(kind=i4_kind), public  :: &
       n_ch_c1, n_ch_c2, n_xc_c1, n_xc_c2
  ! sum of n_xx_uncontry and n_xx_contry:
  ! Total number of functions after local contraction for given ua, l.
  logical, public  :: diagonal
  ! indicates that a diagonal element is to be processed
  logical, public  :: need_ch_ch
  ! if true, [f_k|f_l] integrals are to be calculated
  logical, public  :: need_ch_xc
  ! if true, <f_k|g_l> integrals are to be calculated
  logical, public  :: need_xc_ch
  ! if true, <g_k|f_l> integrals are to be calculated
  logical, public  :: need_xc_xc
  ! if true, <g_k|g_l> integrals are to be calculated
  integer(kind=i4_kind), public  :: n_indep1,n_indep2
  ! number of independent functions of each quadrupel


  ! symmetry-adapted integrals:
  !---- locally contracted and locally contracted
  ! dimensions: (i_exp2,i_exp1,i_if2,i_if1) 
  real(kind=r8_kind), allocatable, public, target  :: &
       loc_loc_int_ch   (:,:,:,:), loc_loc_int_xc   (:,:,:,:), &
       loc_loc_int_ch_xc(:,:,:,:), loc_loc_int_xc_ch(:,:,:,:), &
       loc_loc_int_pre(:,:,:,:)
  !---- globally contracted and locally contracted
  ! dimensions: (i_exp2,i_if2,i_gc1) 
  real(kind=r8_kind), allocatable, public, target  :: &
       glob_loc_int_ch   (:,:,:), glob_loc_int_xc   (:,:,:), &
       glob_loc_int_ch_xc(:,:,:), glob_loc_int_xc_ch(:,:,:), &
       glob_loc_int_pre(:,:,:)
  !---- locally contracted and globally contracted
  ! dimensions: (i_gc2,i_exp1,i_if1) 
  real(kind=r8_kind), allocatable, public, target  :: &
       loc_glob_int_ch   (:,:,:), loc_glob_int_xc   (:,:,:), &
       loc_glob_int_ch_xc(:,:,:), loc_glob_int_xc_ch(:,:,:), &
       loc_glob_int_pre(:,:,:)
  !---- globally contracted and globally contracted
  ! dimensions: (i_gc2,i_gc1) 
  real(kind=r8_kind), allocatable, public, target  :: &
       glob_glob_int_ch   (:,:), glob_glob_int_xc   (:,:), &
       glob_glob_int_ch_xc(:,:), glob_glob_int_xc_ch(:,:), &
       glob_glob_int_pre(:,:)


  !---- counter of number of contributing summands  --
  !     to global contractions 1 in this quadrupel
  ! dimensions: glob_loc_contrib_xx(n_gcx)
  integer(kind=i4_kind), allocatable, public, target :: &
       glob_loc_contrib_ch   (:), glob_loc_contrib_xc   (:), &
       glob_loc_contrib_ch_xc(:), glob_loc_contrib_xc_ch(:), &
       loc_glob_contrib_ch   (:), loc_glob_contrib_xc   (:), &
       loc_glob_contrib_ch_xc(:), loc_glob_contrib_xc_ch(:), &
       glob_loc_contrib_pre(:), loc_glob_contrib_pre(:)

  !------------ public functions and subroutines ----------------
  public int_data_2cff_setup, int_data_2cff_shutdown

  !================================================================
  ! End of public interface of module
  !================================================================
  !------------ Subroutines -------------------------------------

contains

  !**************************************************************
  subroutine int_data_2cff_setup()
    !  Purpose: allocating of storage required all over the module
    !           and setting of public variables
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind)   :: status,i_ir
    !------------ Executable code -------------------------------

    ! init public variables and error checks
    if ( quadrupel%ua1 .lt. 1 .or. quadrupel%ua1 .gt. N_unique_atoms) &
         call error_handler &
         ("int_data_2cff_setup: wrong number of unique atom 1")
    if ( quadrupel%ua2 .lt. 1 .or. quadrupel%ua2 .gt. N_unique_atoms) &
         call error_handler &
         ("int_data_2cff_setup: wrong number of unique atom 2")
    ua1 => unique_atoms(quadrupel%ua1)
    ua2 => unique_atoms(quadrupel%ua2)
    symequivatoms => unique_atom_symequiv(quadrupel%ua1,quadrupel%ua2)

    ! set a flag to tell that a diagonal quadrupel is processed
    ! => evaluation of renormalization coefficient and special storeage modes
    diagonal = (quadrupel%ua1 .eq. quadrupel%ua2) .and. &
               (quadrupel%l1  .eq. quadrupel%l2)

    ! check if one of the angular momenta is above
    ! the maximal angular momenta of the fitbases
    need_ch_ch = integralpar_2cch_no .and. &
         quadrupel%l1 <= ua1%lmax_ch .and. quadrupel%l2 <= ua2%lmax_ch
    need_ch_xc = integralpar_2c_mixed .and. &
         quadrupel%l1 <= ua1%lmax_ch .and. quadrupel%l2 <= ua2%lmax_xc
    need_xc_ch = integralpar_2c_mixed .and. &
         quadrupel%l1 <= ua1%lmax_xc .and. quadrupel%l2 <= ua2%lmax_ch
    need_xc_xc = integralpar_2cxc_no .and. &
         quadrupel%l1 <= ua1%lmax_xc .and. quadrupel%l2 <= ua2%lmax_xc

    ! ---------------------------------------------------------

    if ( need_ch_ch .or. need_ch_xc ) then
       ! associate first set of arguments with the charge fit basis
       select case(quadrupel%l1)
       case(-1)
          ua1_basis_ch => ua1%l_ch(0)
       case(0)
          ua1_basis_ch => ua1%r2_ch
       case default
          ua1_basis_ch => ua1%l_ch(quadrupel%l1)
       end select
       n_ch_exp1 = ua1_basis_ch%N_exponents
       n_ch_contr1 = ua1_basis_ch%N_contracted_fcts
       n_ch_uncontr1 = ua1_basis_ch%N_uncontracted_fcts
       n_ch_c1 = n_ch_uncontr1 + n_ch_contr1
       n_ch_gc1 = ua1%N_glob_cons_ch
    endif
    if ( need_ch_ch .or. need_xc_ch ) then
       ! associate second set of arguments with the charge fit basis
       select case(quadrupel%l2)
       case(-1)
          ua2_basis_ch => ua2%l_ch(0)
       case(0)
          ua2_basis_ch => ua2%r2_ch
       case default
          ua2_basis_ch => ua2%l_ch(quadrupel%l2)
       end select
       n_ch_exp2 = ua2_basis_ch%N_exponents
       n_ch_contr2 = ua2_basis_ch%N_contracted_fcts
       n_ch_uncontr2 = ua2_basis_ch%N_uncontracted_fcts
       n_ch_c2 = n_ch_uncontr2 + n_ch_contr2
       n_ch_gc2 = ua2%N_glob_cons_ch
    endif
    if ( need_xc_xc .or. need_xc_ch ) then
       ! associate first set of arguments with the exchange fit basis
       select case(quadrupel%l1)
       case(-1)
          ua1_basis_xc => ua1%l_xc(0)
       case(0)
          ua1_basis_xc => ua1%r2_xc
       case default
          ua1_basis_xc => ua1%l_xc(quadrupel%l1)
       end select
       n_xc_exp1 = ua1_basis_xc%N_exponents
       n_xc_contr1 = ua1_basis_xc%N_contracted_fcts
       n_xc_uncontr1 = ua1_basis_xc%N_uncontracted_fcts
       n_xc_c1 = n_xc_uncontr1 + n_xc_contr1
       n_xc_gc1 = ua1%N_glob_cons_xc
    endif
    if ( need_xc_xc .or. need_ch_xc ) then
       ! associate second set of arguments with the exchange fit basis
       select case(quadrupel%l2)
       case(-1)
          ua2_basis_xc => ua2%l_xc(0)
       case(0)
          ua2_basis_xc => ua2%r2_xc
       case default
          ua2_basis_xc => ua2%l_xc(quadrupel%l2)
       end select
       n_xc_exp2 = ua2_basis_xc%N_exponents
       n_xc_contr2 = ua2_basis_xc%N_contracted_fcts
       n_xc_uncontr2 = ua2_basis_xc%N_uncontracted_fcts
       n_xc_c2 = n_xc_uncontr2 + n_xc_contr2
       n_xc_gc2 = ua2%N_glob_cons_xc
    endif

    i_ir = get_totalsymmetric_irrep()
    if (quadrupel%l1.eq.0 .or. quadrupel%l1.eq.-1) then
       n_indep1 = 1
    else
       n_indep1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%n_independent_fcts
    endif
    if (quadrupel%l2.eq.0 .or. quadrupel%l2.eq.-1) then
       n_indep2 = 1
    else
       n_indep2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%n_independent_fcts
    endif

    ! allocate the contracted and symmetry adapted integrals
    if (need_ch_ch .and. integralpar_2cch_no) then
       allocate( &
            loc_loc_int_ch(n_ch_c2,n_ch_c1,n_indep2,n_indep1), &
            glob_loc_int_ch(n_ch_c2,n_indep2,n_ch_gc1), &
            loc_glob_int_ch(n_ch_gc2,n_ch_c1,n_indep1), &
            glob_glob_int_ch(n_ch_gc2,n_ch_gc1), &
            glob_loc_contrib_ch(n_ch_gc1), &
            loc_glob_contrib_ch(n_ch_gc2), &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SETUP: allocation (1) failed")
    endif
    if (integralpar_2cch_pre .and. need_ch_ch ) then
       allocate( &
            loc_loc_int_pre(n_ch_c2,n_ch_c1,n_indep2,n_indep1), &
            glob_loc_int_pre(n_ch_c2,n_indep2,n_ch_gc1), &
            loc_glob_int_pre(n_ch_gc2,n_ch_c1,n_indep1), &
            glob_glob_int_pre(n_ch_gc2,n_ch_gc1), &
            glob_loc_contrib_pre(n_ch_gc1), &
            loc_glob_contrib_pre(n_ch_gc2), &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SETUP: allocation (1) failed")
    endif

    if (need_xc_xc .and. integralpar_2cxc_no) then
       allocate( &
            loc_loc_int_xc(n_xc_c2,n_xc_c1,n_indep2,n_indep1), &
            glob_loc_int_xc(n_xc_c2,n_indep2,n_xc_gc1), &
            loc_glob_int_xc(n_xc_gc2,n_xc_c1,n_indep1), &
            glob_glob_int_xc(n_xc_gc2,n_xc_gc1), &
            glob_loc_contrib_xc(n_xc_gc1), &
            loc_glob_contrib_xc(n_xc_gc2), &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SETUP: allocation (2) failed")
    endif
    if (need_ch_xc) then
       allocate( &
            loc_loc_int_ch_xc(n_xc_c2,n_ch_c1,n_indep2,n_indep1), &
            glob_loc_int_ch_xc(n_xc_c2,n_indep2,n_ch_gc1), &
            loc_glob_int_ch_xc(n_xc_gc2,n_ch_c1,n_indep1), &
            glob_glob_int_ch_xc(n_xc_gc2,n_ch_gc1), &
            glob_loc_contrib_ch_xc(n_ch_gc1), &
            loc_glob_contrib_ch_xc(n_xc_gc2), &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SETUP: allocation (3) failed")
    endif
    if (need_xc_ch) then
       allocate( &
            loc_loc_int_xc_ch(n_ch_c2,n_xc_c1,n_indep2,n_indep1), &
            glob_loc_int_xc_ch(n_ch_c2,n_indep2,n_xc_gc1), &
            loc_glob_int_xc_ch(n_ch_gc2,n_xc_c1 ,n_indep1), &
            glob_glob_int_xc_ch(n_ch_gc2,n_xc_gc1), &
            glob_loc_contrib_xc_ch(n_xc_gc1), &
            loc_glob_contrib_xc_ch(n_ch_gc2), &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SETUP: allocation (4) failed")
    endif

  end subroutine int_data_2cff_setup
  !**************************************************************


  !**************************************************************
  subroutine int_data_2cff_shutdown()
    !  Purpose: deallocatimg of storare required all over the module
    !** End of interface ****************************************
    implicit none
    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind)                :: status
    !------------ Executable code -------------------------------

    ! deallocate the contracted and symmetry adapted integrals
    if (need_ch_ch .and. integralpar_2cch_no) then
       deallocate( &
            loc_loc_int_ch, &
            glob_loc_int_ch, &
            loc_glob_int_ch, &
            glob_glob_int_ch, &
            glob_loc_contrib_ch, &
            loc_glob_contrib_ch, &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SHUTDOWN: deallocation (1) failed")
    endif
    if (integralpar_2cch_pre .and. need_ch_ch ) then
       deallocate( &
            loc_loc_int_pre, &
            glob_loc_int_pre, &
            loc_glob_int_pre, &
            glob_glob_int_pre, &
            glob_loc_contrib_pre, &
            loc_glob_contrib_pre, &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SHUTDOWN: deallocation (1) failed")
    endif
    if (need_xc_xc .and. integralpar_2cxc_no) then
       deallocate( &
            loc_loc_int_xc, &
            glob_loc_int_xc, &
            loc_glob_int_xc, &
            glob_glob_int_xc, &
            glob_loc_contrib_xc, &
            loc_glob_contrib_xc, &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SHUTDOWN: deallocation (2) failed")
    endif
    if (need_ch_xc) then
       deallocate( &
            loc_loc_int_ch_xc, &
            glob_loc_int_ch_xc, &
            loc_glob_int_ch_xc, &
            glob_glob_int_ch_xc, &
            glob_loc_contrib_ch_xc, &
            loc_glob_contrib_ch_xc, &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SHUTDOWN: deallocation (3) failed")
    endif
    if (need_xc_ch) then
       deallocate( &
            loc_loc_int_xc_ch, &
            glob_loc_int_xc_ch, &
            loc_glob_int_xc_ch, &
            glob_glob_int_xc_ch, &
            glob_loc_contrib_xc_ch, &
            loc_glob_contrib_xc_ch, &
            STAT=status)
       if(status.ne.0) call error_handler &
            ("INT_DATA_2CFF_SHUTDOWN: deallocation (4) failed")
    endif
  end subroutine int_data_2cff_shutdown
  !**************************************************************


end module int_data_2cff_module

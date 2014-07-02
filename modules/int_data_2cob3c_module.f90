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
!=========================================================
! Public interface of module
!=====================================================================
module  int_data_2cob3c_module
!---------------------------------------------------------------
!
!  Purpose: This is the module containing the data used when
!           calculating the 2 center orbital and three center
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
!           integals (int_data_2cob3c_setup and
!           int_data_2cob3c_shutdown) are included.
!           int_data_2cob3c_setup should be called after
!           quadrupel is set to present value.
!           There is a reset method (int_data_2cob3c_reset)
!           to initialize all required primitive and
!           contracted integrals with zero.
!
!
!  Module called by: ...
!
!
!  References: Publisher Document: Concepts of Integral Part
!
!
!  Author: TB
!  Date: 5/96
!
!
!----------------------------------------------------------------
!== Interrupt of public interface of module =====================
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification
! Author: MM
! Date:   12/96
! Description: Relativistic matrixelements added
!
!----------------------------------------------------------------
!----------------------------------------------------------------
! Modifications
!----------------------------------------------------------------
!
! Modification (Please copy before editing)
! Author: MS
! Date:   3/97
! Description: New datastructures for keeping of the gradients of
!              primitive matrixelements and their contractions
!              have been added
!
! Modification
! Author: MM
! Date:   8/97
! Description: extension of matrix elements to complex spin orbit
!              integrals
!
! Modification (Please copy before editing)
! Author: Uwe Birkenheuer
! Date:   June 1998
! Description: Moving_Unique_Atom concept introduced
!              Split_Gradients concept introduced
!              Gradients for Model_Density_Approach introduced
! Modification (Please copy before editing)
! Author: AS
! Date:   11-12/99
! Description: integrals of electrostatic potential have been added
!
! Modification (Please copy before editing)
! Author: AS
! Date:   3/00
! Description: integrals of electrostatic field have been added
!
!
! Modification (Please copy before editing)
! Author: AS
! Date:   5/00
! Description: integrals to calculate derivatives due to solvent effect
!              has been added
!
! Modification (Please copy before editing)
! Modification
! Author: AM
! Date:   first: Tue Mar  2 1999
! Description: ...
!
!
! Modification Finite Nuclei
! Author: AM
! Date:   Aug-Sep 2004
! Description: to introduce a finite nuclei  (keeping some of the nuclei point)
!  into a *relativistic*  calculation one needs:
!
!     1) V_{fin}(ua)     a FinNuc counterpart of V(ua) to replace some of the latter
!     2) V(ua)           potential of point nucs (for point nucs, eventualy present in a calc)
!     3) PV_{fin}SP(ua)  scalar relativistic FinNuc counterpart of PVSP
!     4) PVSP(ua)        for (still) point nucs
!     5) PV_{fin}XP(ua)  (vector-type) SO-couterpart of PVSP
!     6) PVXP(ua)        (vector-type) for (still) point nucs
!
!  all these ints are treated as 3c, where 3c == ua. After ints are computed
!  one *merges* two models ( SUM[ua] ) depending on which model was chosen for
!  particular ua:
!
!     1) V_{fin}    = MERGE[ua] V_{fin}(ua)    or V(ua)
!     2) PV_{fin}SP = MERGE[ua] PV_{fin}SP(ua) or PVSP(ua)
!     3) PV_{fin}XP = MERGE[ua] PV_{fin}XP(ua) or PVXP(ua)
!
!  these (now 2c) ints are meant to replace their *standard* counterparts
!  V, PVSP, PVXP
!  The three *scalar* int types computed by good old ??_calculate()
!     V(ua), PVSP(ua), and V_{fin}(ua)
!  are treated as triple-sized
!     INT[3*(ua-1)+1]
!     INT[3*(ua-1)+2]
!     INT[3*(ua-1)+3]
!
! Modification (Please copy before editing)
! Author: ...
! Date:   ...
! Description: ...
!
!----------------------------------------------------------------
!== Interrupt end of public interface of module =================

  !------------ all modules used are public !! ------------------
! define FPP_TIMERS 2
# include "def.h"
  use type_module, only:&
       & IK=>i4_kind,&
       & RK=>r8_kind
  use datatype
  use quadrupel_module
  use unique_atom_module
  use symmetry_data_module
  use integralpar_module
  use fit_coeff_module
  use gradient_data_module, only: gradient_data_n_gradients       &
                                , gradient_data_n_spin_gradients  &
                                , gradient_number                 &
                                , gradient_index                  &
                                , calc_cluster_epe_energy
  use options_module, only: options_split_gradients, options_n_spin, &
                            xcmode_model_density, xcmode_extended_mda, &
                            options_xcmode
#ifdef WITH_EPE
  use ewaldpc_module,only: ewpc_n,ewa_allocstat
#endif
  use potential_module, only: N_points
  use elec_static_field_module, only: N_surface_points,totsym_field_length
  use solv_cavity_module, only: with_pc,fixed_pc, to_calc_grads
  use pointcharge_module, only: totsym_grad_pc_length
  use point_dqo_module, only: totsym_grad_dip_length, &
       totsym_grad_quad_length,totsym_grad_oct_length,totsym_grad_rep_length
  use induced_dipoles_module, only: totsym_grad_idip_length
  use symm_adapt_int,&
       & symadapt_totsym_2c_int_type => sa_int_block,&
       & symadapt_totsym_3c_int_type => sa_3c_int_block
  USE_MEMLOG
  use prim_int_store, only: &
       prim_int, &
       prim_3c_int, &
       prim_vec_int, &
       prim_3cv_int

 use cpksdervs_matrices
 use calc3c_switches

!== Interrupt of public interface of module =====================
  implicit none
  save ! save all variables defined in this module
  private
!== Interrupt end of public interface of module =================

  !------------ Declaration of types ----------------------------
  !
  !  DECLARATION OF TYPES MOVED TO SYMM_ADAPT_INT MODULE,
  !  BUT THEY ARE MADE VISIBLE FROM THIS MODULE TOO ...
  !
  !
  !  type, public :: symadapt_totsym_2c_int_type
  !       ! to hold symmetry-adapted total symmetric
  !       ! two center integrals of one IRREP
  !
  !     real(RK), pointer :: int(:,:,:,:)
  !       ! int(n_c2, n_c1, n_if2, n_if1)
  !       ! for relativistic integrals, n_c = n_exp
  !
  !     real(RK), pointer :: int_imag(:,:,:,:)
  !       ! int(n_c2, n_c1, n_if2, n_if1)
  !       ! for relativistic integrals, n_c = n_exp
  !       ! imaginary part of integrals in case of spin orbit
  !  end type symadapt_totsym_2c_int_type
  public symadapt_totsym_2c_int_type

  !  type, public :: symadapt_nottotsym_2c_int_type
  !       ! to hold symmetry-adapted not total symmetric
  !       ! two center integrals of one IRREP
  !
  !     real(RK), pointer :: int(:,:,:,:,:)
  !       ! int(n_c2, n_c1, n_if2, n_if1, n_partner)
  !       ! for relativistic integrals, n_c = n_exp
  !  end type symadapt_nottotsym_2c_int_type
  public symadapt_nottotsym_2c_int_type

  !  type, public :: symadapt_totsym_3c_int_type
  !       ! to hold symmetry-adapted total symmetric
  !       ! three center integrals of one IRREP
  !
  !     real(RK), pointer :: int(:,:,:,:,:)
  !       ! int(n_c2, n_c1, n_ff, n_if2, n_if1)
  !       ! for relativistic integrals, n_c = n_exp
  !
  !     real(RK), pointer :: int_imag(:,:,:,:,:)
  !       ! int(n_c2, n_c1, n_ff, n_if2, n_if1)
  !       ! for relativistic integrals, n_c = n_exp
  !       ! imaginary part of integrals in case of spin orbit
  !  end type symadapt_totsym_3c_int_type
  public symadapt_totsym_3c_int_type

  !  type, public :: symadapt_nottotsym_3cint_type
  !       ! to hold symmetry-adapted not total symmetric
  !       ! three center integrals of one IRREP
  !
  !     real(RK), pointer :: int(:,:,:,:,:,:)
  !       ! int(n_c2, n_c1, n_ff, n_if2, n_if1, n_partner)
  !       ! for relativistic integrals, n_c = n_exp
  !  end type symadapt_nottotsym_3cint_type
  public symadapt_nottotsym_3cint_type


  !------------ Declaration of constants and variables ----------
  type(quadrupel_type), public  :: quadrupel
    ! quadrupel ( unique atom 1, l 1, unique atom 2, l 2 )
    ! that is processed right now
  type(unique_atom_type), pointer, public  :: ua1, ua2
    ! unique atoms processed right now
  type(unique_atom_basis_type), pointer, public :: ua1_basis, ua2_basis
    ! basis of unique atoms processed right now
  real(RK), pointer, public :: contractions1(:,:), contractions2(:,:)
    ! contractions(N_exponents:N_contracted_fcts)
    ! contraction matrices of unique atoms basis 1 and 2
  real(RK), dimension(3), public :: center1, center2
    ! coordinates of centers
  integer(IK),             public  :: n_m1, n_m2
    ! number of mangnetical quantum numbers for l in quadrupel
  integer(IK),             public  :: n_exp1, n_exp2
    ! number of primitive exponents in ua1_basis and ua2_basis
  integer(IK),             public  :: n_contr1, n_contr2
    ! number of contractions in ua1_basis and ua2_basis
  integer(IK),             public  :: n_uncontr1, n_uncontr2
    ! number of uncontracted exponents used as trivial contractions
    ! in ua1_basis and ua2_basis
  integer(IK),             public  :: n_c1, n_c2
    ! number of contracted basis functions in ua1_basis and ua2_basis
    ! n_cx = n_contrx + n_uncontrx
    ! attention: this is not a dimension for symmetry_adapted
    ! relativistic integrals !!
  logical,                           public :: diagonal
    ! true if (ua1 == ua2) .and. (l1 == l2) for quadrupel

  integer(IK),             public  :: grad_dim_index(6)
  ! overlap <a|1|b> or kinetic <a|p^2|b> gradients and sec derivatives
  ! have at most 6=3+3 totally symmetric components.
  ! So the storage is indexed by the shorter index (cf. GRAD_DIM vs N_SPIN_GRADIENTS)
  ! This variable connects both:
  ! grad_dim_index: i -> g = grad_dim_index(i)
  ! where "g" is a true global (tot. symm.) mode index


  !---- Primitive integrals
  ! dimensions:
  ! 2 center integrals : ( n_exp2, n_exp1, n_m2, n_m1 )
  ! 3 center integrals : ( n_exp2, n_exp1, fitfunction, n_m2, n_m1 )
  ! or some:
  ! 3 center integrals : ( n_exp2, n_exp1, n_ua, n_m2, n_m1, n_integrals)
  real(RK), allocatable, public   :: prim_int_dens_mat(:,:,:,:)
  real(RK), pointer, public   :: prim_int_2cob_kin(:,:,:,:)
    ! 2-center orbital integral of kinetic energy
  real(RK), pointer, public   :: prim_int_2cob_nuc(:,:,:,:)
    ! 2-center orbital integral of nuclear attraction
  real(RK), pointer, public   :: prim_int_2cob_nuc_pseudo(:,:,:,:)
    ! 2-center orbital integral of the pseudopot operator
  real(RK), pointer, public   :: prim_int_2cob_ol(:,:,:,:)
    ! 2-center orbital integral of overlap
  real(RK), pointer, public   :: prim_int_2cob_poten(:,:,:,:,:)
  ! 2-center orbital integral of potential
  real(RK), pointer, public   :: prim_int_2cob_field(:,:,:,:,:)
  ! 2-center orbital integral of electrostatic field
  real(RK), pointer, public   :: prim_int_2cob_pvsp(:,:,:,:)
    ! 2-center orbital integral of pv scalar p
  real(RK), pointer, public   :: prim_int_3c_xc(:,:,:,:,:)
    ! 3-center exchange integral
  real(RK), pointer, public   :: prim_int_3c_co(:,:,:,:,:)
    ! 3-center coulomb integral

  !---- Gradients of primitive integrals
  ! root: gradient_data_n_gradients
  ! branch: (n_exp2,nexp1,n_m2,n_m1)
  type(arrmat4),public,pointer :: prim_int_2cob_epe
    ! 2-center orbital integral of attraction to epe

  type(arrmat4),public,pointer :: prim_int_2cob_ol_dervs(:,:)
  type(arrmat4),public,pointer :: prim_int_2cob_ol_grad(:)
  ! gradients of 2-center orbital integral of overlap
  ! (n_exp2,nexp1,n_m2,n_m1)

  type(arrmat4),public,pointer :: prim_int_3cob_grad(:)
  ! sum gradients of 3-center orbitals and nuc and kin
  ! (n_exp2,nexp1,n_m2,n_m1)

  type(arrmat4),public,pointer :: prim_int_coul_dervs(:,:)
  ! sum dervs of 3-center orbitals and nuc and kin
  ! (NGR,NGR)%m(n_exp2,nexp1,n_m2,n_m1)

  type(arrmat4),public,allocatable :: prim_int_2cob_kin_dervs(:,:)
  type(arrmat4),public,allocatable :: prim_int_nucl_dervs(:,:)
  type(arrmat4),public,allocatable :: prim_int_pvsp_dervs(:,:)
  ! 2nd dervs of KIN, NUC, and SR mat. el. PVSP

  type(arrmat4),public,pointer :: prim_int_2cob_ks_grad(:)
  ! orbital-center gradients of the integral of the kohn sham operator
  ! (n_exp2,nexp1,n_m2,n_m1)
  type(arrmat4),public,pointer :: prim_int_3cob_nuc_grad(:)
  ! nuclear-center gradients of the integral of the nuclear potential
  ! (n_exp2,nexp1,n_m2,n_m1)
   type(arrmat4),public,pointer :: prim_int_3cob_epe
  !  integral of the epe center atraction

  ! solvation gradients
  type(arrmat4),public,pointer :: prim_int_3cob_solv_grad(:)      !!!!!!!!!!!!!!!!
  type(arrmat4),public,pointer :: prim_int_3cob_solv_grad_pc(:)      !!!!!!!!!!!!!!!!

  ! gradients and torques on eXternal centers: PC, PD, PQ, PO, IPD
  type(arrmat4),public,pointer :: prim_int_2cob_pc_grad(:)

  type(arrmat4),public,pointer :: prim_int_2cob_pd_grad(:)
  type(arrmat4),public,pointer :: prim_int_2cob_pd_torq(:)

  type(arrmat4),public,pointer :: prim_int_2cob_pq_grad(:)
  type(arrmat4),public,pointer :: prim_int_2cob_pq_torq(:)

  type(arrmat4),public,pointer :: prim_int_2cob_po_grad(:)
  type(arrmat4),public,pointer :: prim_int_2cob_po_torq(:)

  type(arrmat4),public,pointer :: prim_int_2cob_ipd_grad(:)
  type(arrmat4),public,pointer :: prim_int_2cob_ipd_torq(:)

  type(arrmat4),public,pointer :: prim_int_2cob_rc_grad(:)

#ifndef no_cpks_coul_grads
  type(arrmat5),public,pointer :: prim_int_cpks_coul_grad(:)
#endif
  type(arrmat4),public,pointer :: prim_int_3cob_coul_grad(:)
  ! nuclear-center gradients of the integral of the Coulomb potential
  ! (n_exp2,nexp1,n_m2,n_m1)
  type(arrmat4),public,pointer :: prim_int_2cob_nuc_grad(:)
  ! gradients of nuclear attraction ; only used for relativistic gradients
  ! (n_exp2,nexp1,n_m2,n_m1)
  type(arrmat4),public,pointer :: prim_int_2cob_pseudo_grad(:)
  ! gradients of PP operators matrix elements; only used for relativistic gradients
  ! (n_exp2,nexp1,n_m2,n_m1)
  type(arrmat4),public,pointer :: prim_int_2cob_pvsp_grad(:)
  ! gradients of nuclear attraction ; only used for relativistic gradients
  ! (n_exp2,nexp1,n_m2,n_m1)
  type(arrmat4),public,pointer :: prim_int_2cob_kin_grad(:)
  ! gradients of kinetic ; only used for relativistic gradients
  ! (n_exp2,nexp1,n_m2,n_m1)

  !---- Symmetry-adapted integrals
  ! dimension: Irreps
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_kin(:)
    ! 2-center orbital integral of kinetic energy
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_nuc(:)
    ! 2-center orbital integral of nuclear attraction
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_nuc_pseudo(:)
    ! 2-center orbital integral of pseudopot operators
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pvsp(:)
  ! 2-center orbital integral of pvsp
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ol(:)
    ! 2-center orbital integral of overlap
  type(symadapt_totsym_3c_int_type), pointer, public :: symadapt_int_2cob_poten(:)
    ! 2-center orbital integral of potential
  type(symadapt_totsym_3c_int_type), pointer, public :: symadapt_int_2cob_field(:)
    ! 2-center orbital integral of of el.static field
  type(symadapt_totsym_3c_int_type), pointer, public :: symadapt_int_3c_xc(:)
    ! 3-center exchange integral
  type(symadapt_totsym_3c_int_type), allocatable, public :: symadapt_int_3c_co(:) ! (n_irr)
    ! 3-center coulomb integral

  !
  ! SB: New type of 3-center Coulomb integral need for response(ira,irb,irc)%m(ia,ib,ic,ifa,ifb)
  !
  ! For excitaitons one needs also
  !
  !     [ ij | k ]
  !
  ! integrals with all of i, j and k of different symmetries. We store only
  ! the symmtry reduced combinatios where summations over symmetry equivalent
  ! functions (irrep partners) of i, j, and k and corresponding Clebsch-Gordan
  ! coefficients has been performed. In fact, the integrals for totally
  ! symmetric fit functions k are just a special case of these integrals.
  !
  ! FIXME: is this redundacy exploited anywhere?
  !
  type, public::  multiplicity
     !
     ! The size of this array depends on the three irrep indices (ira, irb, irc)
     ! and corresponds to the number of ways the direct product of the three irreps
     ! can be reduced to the totally symmetric combination. This is the same
     ! multiplicity as that of the Clebsch-Gordan coefficients cg(ira, irb, irc)%mult.
     !
     type(symadapt_totsym_3c_int_type), allocatable :: mult(:)
  end type multiplicity
  type(multiplicity), allocatable, public :: symadapt_int_3c_co_resp(:,:,:) ! (n_irr, n_irr, n_irr)

    ! 3-center coulomb integral

  !---- Gradients of symmetry-adapted integrals
  ! dimension: Irreps, number of gradients
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ol_grad(:,:)
  ! gradient of  overlap
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ol_dervs(:,:,:)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_3cob_grad(:,:)
  ! gradient of  three center orbitals and kin and nuc

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_coul_dervs(:,:,:)
  ! (NGR,NGR,N_IRR)%m(n_c2,n_c1,n_if2,n_if1)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_olu_dervs(:,:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_kin_dervs(:,:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_nucl_dervs(:,:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_pvsp_dervs(:,:,:)
  ! ONLY IN RELATIVISTIC CALCULATION:
  !   (uncontracted) overlap, kinetic, nuclear and pvsp sec. der.
  ! (NGR,NGR,N_IRR)%m(n_e2,n_e1,n_if2,n_if1)

  ! orbital gradient of kohn sham hamiltonian (only used for split_gradients)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ks_grad(:,:)
  ! fitfct gradient of nuclear potential (only used for split_gradients)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_3cob_nuc_grad(:,:)
  ! epe atraction potential (only used for ewpc_n.ne.0)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_3cob_epe(:,:)
  ! symmetry adapted solvation gradients

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_3cob_solv_grad(:,:)  !!!!!!!!!!!!
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_3cob_solv_grad_pc(:,:)  !!!!!!!!!!!!

  ! gradients and torques on eXternal centers: PC, PD, PQ, PO, IPD
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pc_grad(:,:)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pd_grad(:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pd_torq(:,:)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pq_grad(:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pq_torq(:,:)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_po_grad(:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_po_torq(:,:)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ipd_grad(:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ipd_torq(:,:)

  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_rc_grad(:,:)

  ! fitfct gradient of Coulomb potential (only used for split_gradients)
#ifndef no_cpks_coul_grads
  type(symadapt_totsym_3c_int_type), pointer, public :: symadapt_int_cpks_coul_grad(:,:)
#endif
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_3cob_coul_grad(:,:)
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_kin_grad(:,:)
  ! gradient of kin attraction; only used for relativistic gradients
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_nuc_grad(:,:)
  ! gradient of  nuclear attraction; only used for relativistic gradients
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pseudo_grad(:,:)
  ! gradient of  PP operator matrix elements; only used for relativistic gradients
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_epe(:)
  ! energy of interection with epe; only used for relativistic gradients
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_pvsp_grad(:,:)
  ! gradient of pvsp integral; only used for relativistic gradients
  type(symadapt_totsym_2c_int_type), pointer, public :: symadapt_int_2cob_ol_rel_grad(:,:)
  ! gradient of kin attraction; only used for relativistic gradients

  ! FIXME: the rest wants here too:
  type(prim_int),     allocatable, public :: prim_2c_ints(:)
  type(prim_vec_int), allocatable, public :: prim_2cv_ints(:)
  type(prim_3c_int),  allocatable, public :: prim_3c_ints(:)
  type(prim_3cv_int), allocatable, public :: prim_3cv_ints(:)

  type(symadapt_totsym_2c_int_type), allocatable, public  :: sa_2c_ints(:, :) ! (n_irr, n_2c_ints)
  type(symadapt_totsym_3c_int_type), allocatable, public  :: sa_3c_ints(:, :) ! (n_irr, n_3c_ints)
  type(symadapt_totsym_2c_int_type), allocatable, public  :: sa_alphap(:, :) ! (n_irr, 2)

  ! integer (uninitialized) handles for 2c-SA ints:
  integer(IK),              public  ::&
       & int_sa_kin  = -1,&
       & int_sa_nuc  = -1,&
       & int_sa_ol   = -1,&
       & int_sa_pvsp = -1,&
       & int_sa_pvxp = -1,&
       & int_sa_sigp = -1

  ! integer (uninitialized) handles for 3c-SA ints:
  integer(IK),              public  ::&
       & int_sa_xc   = -1,&
       & int_sa_co   = -1,&
       & int_sa_rcoul_pvsp = -1, &
       & int_sa_rcoul_pvxp = -1, &
       & int_sa_r2_pvsp = -1, &
       & int_sa_r2_pvxp = -1

  ! integer (uninitialized) handles for 2c-PRIM scalar-type ints:
  integer(IK),              public  ::&
       & int_prim_kin  = -1,&
       & int_prim_nuc  = -1,&
       & int_prim_ol   = -1,&
       & int_prim_pvsp = -1,&
       & int_prim_nucfin = -1,&
       & int_prim_pfinsp = -1

  ! integer (uninitialized) handles for 2c-PRIM vector-type ints:
  integer(IK),              public  ::&
       & int_prim_pvxp = -1,&
       & int_prim_pvec = -1,&
       & int_prim_pfinxp = -1

  ! integer (uninitialized) handles for 3c-PRIM scalar-type ints:
  integer(IK),              public  ::&
       & int_prim_xc   = -1,&
       & int_prim_co   = -1,&
       & int_prim_rcoul_pvsp = -1, &
       & int_prim_r2_pvsp = -1, &
       & int_prim_many_3c = -1, &
       & int_prim_split_pfinsp = -1


  ! integer (uninitialized) handles for 3c-PRIM vector-type ints:
  integer(IK),              public  ::&
       & int_prim_rcoul_pvxp = -1, &
       & int_prim_r2_pvxp = -1, &
       & int_prim_split_pvxp = -1, &
       & int_prim_split_pfinxp = -1

  integer(i4_kind), parameter, public :: &
       OFF_PVSP   = 1, &
       OFF_V      = 2, &
       OFF_VFIN   = 3, &
       OFF_STRIDE = 3 ! multiplicity of the many_3c ints

  !------------ public functions and subroutines ----------------
  public &
       int_data_2cob3c_setup, &
       int_data_2cob3c_shutdown, &
       allocate_primitives


  !===================================================================
  ! End of public interface of module
  !===================================================================


  real(r8_kind) :: zero = 0.0_r8_kind
!!$  integer(IK) :: N_pc !!!!!!!!!!!!1

  !--------------------------------------------------------------
  !------------ Subroutines -------------------------------------
contains


  subroutine int_data_2cob3c_setup()
    !  Purpose: setting of public variables and
    !           allocation of storage required all over the module
    !           (in particular symmetry-adapted integrals)
    !** End of interface ****************************************
    implicit none

    !------------ Declaration of local variables ----------------
    integer(IK) :: n_ir, n_ff_ch, n_ff_xc, dim_c1, dim_c2,&
         & n_ff_s,n_ff_r2
    !! SB: variables
    integer(IK) :: i_ir_a, i_ir_b, n_if1_resp, n_if2_resp
    integer(IK) :: k2dr, imult
    !------------ Executable code -------------------------------

    ! init public variables and error checks
    if ( quadrupel%ua1 .lt. 1 .or. quadrupel%ua1 .gt. N_unique_atoms) &
         call error_handler("int_data_2cob3c_setup: wrong number of unique atom 1")
    if ( quadrupel%ua2 .lt. 1 .or. quadrupel%ua2 .gt. N_unique_atoms) &
         call error_handler("int_data_2cob3c_setup: wrong number of unique atom 2")
    ua1 => unique_atoms(quadrupel%ua1)
    ua2 => unique_atoms(quadrupel%ua2)

    if ( quadrupel%l1 .lt. 0 .or. quadrupel%l1 .gt. ua1%lmax_ob) &
         call error_handler("int_data_2cob3c_setup: wrong l for unique atom 1")
    if ( quadrupel%l2 .lt. 0 .or. quadrupel%l2 .gt. ua2%lmax_ob) &
         call error_handler("int_data_2cob3c_setup: wrong l for unique atom 2")
    ua1_basis => ua1%l_ob(quadrupel%l1)
    ua2_basis => ua2%l_ob(quadrupel%l2)
    contractions1 => ua1_basis%contractions
    contractions2 => ua2_basis%contractions

    n_m1 = ( 2 * quadrupel%l1 ) + 1
    n_m2 = ( 2 * quadrupel%l2 ) + 1
    n_exp1 = ua1_basis%N_exponents
    n_exp2 = ua2_basis%N_exponents
    n_contr1 = ua1_basis%N_contracted_fcts
    n_contr2 = ua2_basis%N_contracted_fcts
    n_uncontr1 = ua1_basis%N_uncontracted_fcts
    n_uncontr2 = ua2_basis%N_uncontracted_fcts
    n_c1 = n_contr1 + n_uncontr1
    n_c2 = n_contr2 + n_uncontr2
    diagonal = (quadrupel%ua1 == quadrupel%ua2) .and. &
               (quadrupel%l1 == quadrupel%l2)

    ! local variables required for allocation
    n_ff_ch  = fit_coeff_n_ch()
    n_ff_s = fit_coeff_n_s()
    n_ff_r2 = fit_coeff_n_r2()
    n_ff_xc  = fit_coeff_n_xc()
    if ( integralpar_relativistic) then
       dim_c1 = n_exp1
       dim_c2 = n_exp2
    else
       dim_c1 = n_c1
       dim_c2 = n_c2
    endif
    n_ir = symmetry_data_n_irreps()

    if(integralpar_spor)then
       call alloc_sa_int_spor()
    else
       call alloc_sa_int()
    endif

!   print*,'call allocate_primitives',INT_PRIM_ALLOCATE
    call allocate_primitives('allocate')

  contains

    subroutine alloc_sa_int()
    use fit_coeff_module, only: get_fit
#ifdef WITH_RESPONSE
      use ch_response_module, only: dimension_of_fit_ch
      use clebsch_gordan, only: cg=>cg_eliminated
#endif
      implicit none
    type(fit) :: n_fit

      ! internal sub
      ! *** end of interface ***
      integer(IK) :: status, ma1, ma2, i_ir, n_if2, n_if1,&
           & n_pa, i_grad, grad_dim, spin_grad_dim
      integer(IK) :: mult_irrep
      logical :: model_density, split_gradients

      !---------------------------------
      !
      ! STANDARD SCF (NO SPIN ORBIT)
      !

      !=================================
      !
      ! ALLOCATE TWO CENTER INTEGRALS >>>
      !

!     print*,'call alloc sa int'
      if ( integralpar_2cob_kin ) then
         allocate( symadapt_int_2cob_kin(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_2cob_kin failed" )
      endif
      if ( integralpar_2cob_nuc ) then
         allocate( symadapt_int_2cob_nuc(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_2cob_nuc failed" )
         if((pseudopot_present.and.integralpar_relativistic) .or. &
              integralpar_pseudo) & !AS
              allocate( symadapt_int_2cob_nuc_pseudo(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_2cob_nuc_pseudo failed" )
      endif
      if ( integralpar_2cob_pvsp ) then
         allocate( symadapt_int_2cob_pvsp(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_2cob_pvsp failed" )
      end if
      if ( integralpar_2cob_ol ) then
         allocate( symadapt_int_2cob_ol(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_2cob_ol failed" )
      endif

      split_gradients = options_split_gradients()
      model_density = options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda
      if (integralpar_gradients .or. integralpar_rel_gradients) then
         ma1 = ua1%moving_atom
         ma2 = ua2%moving_atom
         ASSERT(ma1>0)
         ASSERT(ma2>0)

         ! setup grad_dim_index:
         grad_dim_index(:) = 0
         grad_dim = 0

         ! over a-grads
         do i_grad=1,gradient_number(ma1)
           grad_dim = grad_dim + 1
           grad_dim_index(grad_dim) = gradient_index(ma1) + i_grad - 1
         enddo

         ! over b-grads
         do i_grad=1,gradient_number(ma2)
           grad_dim = grad_dim + 1
           grad_dim_index(grad_dim) = gradient_index(ma2) + i_grad - 1
         enddo
         ASSERT(grad_dim<=6)

         if (model_density) then
            spin_grad_dim = grad_dim * options_n_spin()
         else
            spin_grad_dim = grad_dim
         endif
      endif

      if ( integralpar_2cob_ol_grad ) then

         allocate(prim_int_2cob_ol_grad (grad_dim), &
                  symadapt_int_2cob_ol_grad(grad_dim,n_ir), &
                                        stat=status )
            ASSERT(status.eq.0)

         if(integralpar_dervs) then
          allocate(prim_int_2cob_ol_dervs(grad_dim,grad_dim), &          !!! alloc_sa
                   symadapt_int_2cob_ol_dervs(grad_dim,grad_dim,n_ir), &
                                            stat=cpksalloc(45))
          ASSERT(cpksalloc(45).eq.0)
          ol_dervs_size=0
         endif

      endif

      if ( integralpar_2cob_ol_rel_grad ) then
         allocate( symadapt_int_2cob_ol_rel_grad(grad_dim,&
              n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of symadapt_int_2cob_ol_rel_grad failed" )
      end if



      !=======================================================
      !
      ! ALLOCATE THREE CENTER INTEGRALS
      !
      !---------------------------------

      if ( integralpar_3c_xc) then
         allocate( symadapt_int_3c_xc(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_3c_xc failed" )
      endif
      if ( integralpar_3c_co) then
         allocate( symadapt_int_3c_co(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_3c_co failed" )
      endif
      !! SB: allocate for response
      if ( integralpar_3c_co_resp) then
         allocate( symadapt_int_3c_co_resp(n_ir, n_ir, n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_3c_co_resp failed" )
      endif
      if ( integralpar_2cob_potential ) then
         allocate( symadapt_int_2cob_poten(n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_setup: allocate of symadapt_int_2cob_poten failed" )
      endif
      if ( integralpar_2cob_field ) then
         allocate( symadapt_int_2cob_field(n_ir), stat=status )
         ASSERT(status==0)
      endif

      if (integralpar_2cob_nuc_grad) then
         allocate( symadapt_int_2cob_nuc_grad(gradient_data_n_gradients,&
              n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of symadapt_int_2cob_nuc_grad failed" )
         allocate( prim_int_2cob_nuc_grad(gradient_data_n_gradients)&
              , stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of prim_int_2cob_nuc_grad failed" )
         if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
            allocate( symadapt_int_2cob_pseudo_grad(gradient_data_n_gradients,&
                 n_ir), stat=status )
            if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
                 &allocate of symadapt_int_2cob_nuc_pseudo_grad failed" )
            allocate( prim_int_2cob_pseudo_grad(gradient_data_n_gradients)&
                 , stat=status )
            if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
                 &allocate of prim_int_2cob_pseudo_grad failed" )
         endif! pseudopot_presen
      end if
      if (integralpar_2cob_kin_grad) then
         allocate( symadapt_int_2cob_kin_grad(grad_dim,&
              n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of symadapt_int_2cob_kin_grad failed" )
         allocate( prim_int_2cob_kin_grad(grad_dim)&
              , stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of prim_int_2cob_kin_grad failed" )
      end if
      if (integralpar_2cob_pvsp_grad) then
         allocate( symadapt_int_2cob_pvsp_grad(gradient_data_n_gradients,&
              n_ir), stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of symadapt_int_2cob_pvsp_grad failed" )
         allocate( prim_int_2cob_pvsp_grad(gradient_data_n_gradients)&
              , stat=status )
         if ( status .ne. 0 ) call error_handler("int_data_2cob3c_setup: &
              &allocate of prim_int_2cob_pvsp_grad failed" )
      end if
      if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
         DPRINT  '2cob3c: allocate prim_int_3cob_solv_grad ...'
         allocate( prim_int_3cob_solv_grad(gradient_data_n_spin_gradients), &
              stat=status )
         if (status /= 0) call error_handler("int_data_2cob3c_setup: &
              &allocate of prim_int_3cob_solv_grad failed" )
         allocate( symadapt_int_3cob_solv_grad(gradient_data_n_spin_gradients,n_ir),&
              stat=status )
         if (status /= 0) call error_handler("int_data_2cob3c_setup: &
              &allocate of symadapt_int_3cob_nuc_grad failed" )

         if(with_pc .and. .not.fixed_pc) then
            allocate( prim_int_3cob_solv_grad_pc(totsym_field_length), &
                 stat=status )
            if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                 &allocate of prim_int_3cob_solv_grad_pc failed" )
            allocate( symadapt_int_3cob_solv_grad_pc(totsym_field_length,n_ir),&
                 stat=status )
            if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                 &allocate of symadapt_int_3cob_nuc_grad_pc failed" )
         end if
      endif

      if(integralpar_2cob_pc_grad) then
         allocate( prim_int_2cob_pc_grad(totsym_grad_pc_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_pc_grad(totsym_grad_pc_length,&
              n_ir), stat=status )
         ASSERT(status==0)
      end if
      if(integralpar_2cob_X_grad) then
         allocate( prim_int_2cob_pd_grad(totsym_grad_dip_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_pd_grad(totsym_grad_dip_length,&
              n_ir), stat=status )
         ASSERT(status==0)
         allocate( prim_int_2cob_pd_torq(totsym_grad_dip_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_pd_torq(totsym_grad_dip_length,&
              n_ir), stat=status )
         ASSERT(status==0)

         allocate( prim_int_2cob_pq_grad(totsym_grad_quad_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_pq_grad(totsym_grad_quad_length,&
              n_ir), stat=status )
         ASSERT(status==0)
         allocate( prim_int_2cob_pq_torq(totsym_grad_quad_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_pq_torq(totsym_grad_quad_length,&
              n_ir), stat=status )
         ASSERT(status==0)

         allocate( prim_int_2cob_po_grad(totsym_grad_oct_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_po_grad(totsym_grad_oct_length,&
              n_ir), stat=status )
         ASSERT(status==0)
         allocate( prim_int_2cob_po_torq(totsym_grad_oct_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_po_torq(totsym_grad_oct_length,&
              n_ir), stat=status )
         ASSERT(status==0)

         allocate( prim_int_2cob_rc_grad(totsym_grad_rep_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_rc_grad(totsym_grad_rep_length,&
              n_ir), stat=status )
         ASSERT(status==0)
      end if
      if(integralpar_2cob_ipd_grad) then
         allocate( prim_int_2cob_ipd_grad(totsym_grad_idip_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_ipd_grad(totsym_grad_idip_length,&
              n_ir), stat=status )
         ASSERT(status==0)
         allocate( prim_int_2cob_ipd_torq(totsym_grad_idip_length)&
              , stat=status )
         ASSERT(status==0)
         allocate( symadapt_int_2cob_ipd_torq(totsym_grad_idip_length,&
              n_ir), stat=status )
         ASSERT(status==0)
      end if

      grad_3cob: if ( integralpar_3cob_grad ) then
#ifdef WITH_EPE
         if(ewpc_n.ne.0.and.calc_cluster_epe_energy) then
            allocate(prim_int_3cob_epe &
                    ,symadapt_int_3cob_epe(1,n_ir), &
              stat=ewa_allocstat(5))
              ASSERT(ewa_allocstat(5).eq.0)
              ewa_allocstat(5)=1
         endif
#endif
          allocate( symadapt_int_3cob_grad(gradient_data_n_spin_gradients,n_ir)   &
                   ,prim_int_3cob_grad(gradient_data_n_spin_gradients)  &
                   ,stat=status)
          ASSERT(status.eq.0)

          if(integralpar_2dervs) then
             allocate( prim_int_coul_dervs &
                       (gradient_data_n_spin_gradients,gradient_data_n_spin_gradients)&
                      ,symadapt_int_coul_dervs &
                       (gradient_data_n_spin_gradients,gradient_data_n_spin_gradients,n_ir)&
                      ,stat=cpksalloc(36))
                  ASSERT(cpksalloc(36).eq.0)
          size_coul_dervs=0
          endif

          if(integralpar_2dervs.and.integralpar_relativistic) then
             allocate( prim_int_nucl_dervs &
                       (gradient_data_n_spin_gradients,gradient_data_n_spin_gradients)&
                      ,symadapt_int_nucl_dervs &
                       (gradient_data_n_spin_gradients,gradient_data_n_spin_gradients,n_ir)&
                      ,stat=status)
             ASSERT(status.eq.0)
             allocate( prim_int_pvsp_dervs &
                       (gradient_data_n_spin_gradients,gradient_data_n_spin_gradients)&
                      ,symadapt_int_pvsp_dervs &
                       (gradient_data_n_spin_gradients,gradient_data_n_spin_gradients,n_ir)&
                      ,stat=status)
             ASSERT(status.eq.0)
             ! no need for separate (uncontracted) primitive overlap:
             allocate( prim_int_2cob_kin_dervs(grad_dim,grad_dim)&
                      ,symadapt_int_2cob_kin_dervs(grad_dim,grad_dim,n_ir)&
                      ,symadapt_int_2cob_olu_dervs(grad_dim,grad_dim,n_ir)&
                      ,stat=status)
             ASSERT(status.eq.0)
          endif

         splitg: if (split_gradients) then
            allocate( prim_int_2cob_ks_grad(spin_grad_dim), &
                 stat=status )
            if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                 &allocate of prim_int_2cob_ks_grad failed" )
            allocate( symadapt_int_2cob_ks_grad(spin_grad_dim,n_ir), &
                 stat=status )
            if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                 &allocate of symadapt_int_2cob_ks_grad failed" )
            allocate( prim_int_3cob_nuc_grad(gradient_data_n_gradients), &
                 stat=status )
            if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                 &allocate of prim_int_3cob_nuc_grad failed" )
            allocate( symadapt_int_3cob_nuc_grad(gradient_data_n_gradients,n_ir),&
                 stat=status )
            ASSERT(status.eq.0)

            if (model_density) then
               allocate( prim_int_3cob_coul_grad(gradient_data_n_gradients),stat=status )
                    ASSERT(status.eq.0)
               allocate( symadapt_int_3cob_coul_grad(gradient_data_n_gradients, n_ir),stat=status )
                    ASSERT(status.eq.0)
            endif
         elseif(integralpar_cpksdervs) then
          call get_fit(n_fit)
#ifndef no_cpks_coul_grads
          allocate( prim_int_cpks_coul_grad(gradient_data_n_gradients), &
                    symadapt_int_cpks_coul_grad(gradient_data_n_gradients, n_ir), &
                    stat=cpksalloc(17))
          ASSERT(cpksalloc(17).eq.0)
#endif
          cpks_coul_grad_size=0

         endif splitg
      endif grad_3cob

      do i_ir = 1, n_ir
         n_if1 = ua1%symadapt_partner(i_ir,quadrupel%l1)%N_independent_fcts
         n_if2 = ua2%symadapt_partner(i_ir,quadrupel%l2)%N_independent_fcts
         n_pa = symmetry_data_n_partners(i_ir)
         ! in case of scalarrelativistic calculations
         if ( integralpar_2cob_kin ) then
            allocate( symadapt_int_2cob_kin(i_ir)%int(dim_c2,dim_c1,n_if2,n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_2cob_kin%int failed" )
            symadapt_int_2cob_kin(i_ir)%int = 0.0_r8_kind
         endif
         if ( integralpar_2cob_nuc ) then
            allocate( symadapt_int_2cob_nuc(i_ir)%int(dim_c2,dim_c1,n_if2,n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_2cob_nuc%int failed" )
            symadapt_int_2cob_nuc(i_ir)%int = 0.0_r8_kind
            if ((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
               allocate( symadapt_int_2cob_nuc_pseudo(i_ir)%int(dim_c2,dim_c1,n_if2,n_if1), &
                    stat=status )
               symadapt_int_2cob_nuc_pseudo(i_ir)%int = 0.0_r8_kind
            endif! pseudopot_present.and.integralpar_relativistic
         endif
         if ( integralpar_2cob_pvsp ) then
            allocate( symadapt_int_2cob_pvsp(i_ir)%int(dim_c2,dim_c1,n_if2,n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_2cob_pvsp%int failed" )
            symadapt_int_2cob_pvsp(i_ir)%int = 0.0_r8_kind
         endif
         if ( integralpar_2cob_ol ) then
            allocate( symadapt_int_2cob_ol(i_ir)%int(dim_c2,dim_c1,n_if2,n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_2cob_ol%int failed" )
            symadapt_int_2cob_ol(i_ir)%int = 0.0_r8_kind
         endif
         if ( integralpar_2cob_ol_grad ) then
            do i_grad=1,grad_dim
               allocate( symadapt_int_2cob_ol_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status.eq.0)
               symadapt_int_2cob_ol_grad(i_grad,i_ir)%int = 0.0_r8_kind

               if(integralpar_dervs) then
                do k2dr=1,grad_dim
                    allocate( symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%&    !!!alloc_sa
                              int(n_c2,n_c1,n_if2,n_if1), stat=cpksalloc(55) )
                ASSERT(cpksalloc(55).eq.0)
                MEMLOG(size(symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int))
                ol_dervs_size=ol_dervs_size+size(symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int)
                 symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int = 0.0_r8_kind
                enddo
               endif
            end do
           if(integralpar_dervs) then
            if(ol_dervs_size.gt.max_ol_dervs_size) then
             max_ol_dervs_size=ol_dervs_size
             DPRINT 'max_ol_dervs_size ',max_ol_dervs_size
            endif
           endif
         endif
         if ( integralpar_2cob_ol_rel_grad ) then
            do i_grad=1,grad_dim
               allocate( symadapt_int_2cob_ol_rel_grad(i_grad,i_ir)%&
                    int(dim_c2,dim_c1,n_if2,n_if1), stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: allocate of &
                    &symadapt_int_2cob_ol_grad%int failed" )
               symadapt_int_2cob_ol_rel_grad(i_grad,i_ir)%int = 0.0_r8_kind
            end do
         endif
         if ( integralpar_2cob_kin_grad ) then
            do i_grad=1,grad_dim
               allocate( symadapt_int_2cob_kin_grad(i_grad,i_ir)%&
                    int(dim_c2,dim_c1,n_if2,n_if1), stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: allocate of &
                    &symadapt_int_2cob_kin_grad%int failed" )
               symadapt_int_2cob_kin_grad(i_grad,i_ir)%int = 0.0_r8_kind
            end do
         endif
         if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
            do i_grad=1,gradient_data_n_spin_gradients
               allocate( symadapt_int_3cob_solv_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                    &allocate of symadapt_int_3cob_solv_grad%int failed" )
               symadapt_int_3cob_solv_grad(i_grad,i_ir)%int = 0.0_r8_kind
            end do

            if(with_pc .and. .not.fixed_pc) then
               do i_grad=1,totsym_field_length
                  allocate( symadapt_int_3cob_solv_grad_pc(i_grad,i_ir)%&
                       int(n_c2,n_c1,n_if2,n_if1), stat=status )
                  if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                       &allocate of symadapt_int_3cob_solv_grad_pc%int failed" )
                  symadapt_int_3cob_solv_grad_pc(i_grad,i_ir)%int = 0.0_r8_kind
               end do
            end if
         endif

         grad_3cob1:if ( integralpar_3cob_grad ) then

#ifdef WITH_EPE
        if(ewpc_n.ne.0.and.calc_cluster_epe_energy) then
             allocate(symadapt_int_3cob_epe(1,i_ir)%int(n_c2,n_c1,n_if2,n_if1), &
                           stat=ewa_allocstat(8) )
          ASSERT(ewa_allocstat(8).eq.0)
                 ewa_allocstat(8)=1
                  symadapt_int_3cob_epe(1,i_ir)%int=0.0_r8_kind
         endif
#endif

       do i_grad=1,gradient_data_n_spin_gradients

         allocate( symadapt_int_3cob_grad(i_grad,i_ir)%int(n_c2,n_c1,n_if2,n_if1), &
                   stat=status )
         ASSERT(status.eq.0)
         symadapt_int_3cob_grad(i_grad,i_ir)%int = 0.0_r8_kind
!         print*,'setup: symadapt_int_3cob_grad', &
!                 shape(symadapt_int_3cob_grad(i_grad,i_ir)%int),i_grad,i_ir

            if(integralpar_2dervs) then
             do k2dr=1,gradient_data_n_spin_gradients
               allocate( symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int(n_c2,n_c1,n_if2,n_if1), &
                         stat=cpksalloc(38))
               ASSERT(cpksalloc(38).eq.0)
               MEMLOG(size(symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int))
               symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int = 0.0_r8_kind
               size_coul_dervs=size_coul_dervs+size(symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int)
             enddo
            endif
       enddo
       if(integralpar_2dervs) then
        if(size_coul_dervs.gt.max_size_coul_dervs) then
         max_size_coul_dervs=size_coul_dervs
         DPRINT 'int_data_2cob3c_setup: max_size_coul_dervs ',max_size_coul_dervs
        endif
       endif

       ! FIXME: mv out of if(integralpar_3cob_grad) ???
       if(integralpar_2dervs.and.integralpar_relativistic)then
       do i_grad=1,gradient_data_n_spin_gradients
         do k2dr=1,gradient_data_n_spin_gradients

           allocate( symadapt_int_nucl_dervs(i_grad,k2dr,i_ir)&
                    %int(n_exp2,n_exp1,n_if2,n_if1)           &
                    ,stat=status)
           ASSERT(status.eq.0)
           symadapt_int_nucl_dervs(i_grad,k2dr,i_ir)%int = 0.0_r8_kind

           allocate( symadapt_int_pvsp_dervs(i_grad,k2dr,i_ir)&
                    %int(n_exp2,n_exp1,n_if2,n_if1)           &
                    ,stat=status)
           ASSERT(status.eq.0)
           symadapt_int_pvsp_dervs(i_grad,k2dr,i_ir)%int = 0.0_r8_kind
         enddo
       enddo
       do i_grad=1,size(symadapt_int_2cob_kin_dervs,1)
         do k2dr=1,size(symadapt_int_2cob_kin_dervs,2)

           allocate( symadapt_int_2cob_olu_dervs(i_grad,k2dr,i_ir)&
                    %int(n_exp2,n_exp1,n_if2,n_if1)           &
                    ,stat=status)
           ASSERT(status.eq.0)
           symadapt_int_2cob_olu_dervs(i_grad,k2dr,i_ir)%int = 0.0_r8_kind

           allocate( symadapt_int_2cob_kin_dervs(i_grad,k2dr,i_ir)&
                    %int(n_exp2,n_exp1,n_if2,n_if1)           &
                    ,stat=status)
           ASSERT(status.eq.0)
           symadapt_int_2cob_kin_dervs(i_grad,k2dr,i_ir)%int = 0.0_r8_kind
         enddo
       enddo
       endif

            splitg1:if (split_gradients) then
               do i_grad=1,spin_grad_dim
                  allocate( symadapt_int_2cob_ks_grad(i_grad,i_ir)%&
                       int(n_c2,n_c1,n_if2,n_if1), stat=status )
                       ASSERT(status.eq.0)
                  symadapt_int_2cob_ks_grad(i_grad,i_ir)%int = 0.0_r8_kind
               end do
!!!! MF: bug-fix: moved after merge (AG)
!!!!           do i_grad=1,gradient_data_n_spin_gradients
               do i_grad=1,gradient_data_n_gradients
                  allocate( symadapt_int_3cob_nuc_grad(i_grad,i_ir)%&
                       int(n_c2,n_c1,n_if2,n_if1), stat=status )
                  if (status /= 0) call error_handler("int_data_2cob3c_setup: &
                       &allocate of symadapt_int_3cob_nuc_grad%int failed" )
                  symadapt_int_3cob_nuc_grad(i_grad,i_ir)%int = 0.0_r8_kind
               end do
               if (model_density) then
!!!! MF: bug-fix: moved after merge (AG)
!!!!            do i_grad=1,gradient_data_n_spin_gradients
                  do i_grad=1,gradient_data_n_gradients
                     allocate( symadapt_int_3cob_coul_grad(i_grad,i_ir)%&
                          int(n_c2,n_c1,n_if2,n_if1), stat=status )
                          ASSERT(status.eq.0)
                     symadapt_int_3cob_coul_grad(i_grad,i_ir)%int = 0.0_r8_kind
                  end do
               endif
#ifndef no_cpks_coul_grads
            elseif(integralpar_cpksdervs) then
                  do i_grad=1,gradient_data_n_gradients
                   allocate( symadapt_int_cpks_coul_grad(i_grad,i_ir)%int &
                                     (n_c2,n_c1,n_fit%n_ch,n_if2,n_if1), stat=cpksalloc(23) )
             ASSERT(cpksalloc(23).eq.0)
             MEMLOG(size(symadapt_int_cpks_coul_grad(i_grad,i_ir)%int))
             cpks_coul_grad_size=cpks_coul_grad_size+ &
                                              size(symadapt_int_cpks_coul_grad(i_grad,i_ir)%int)
                     symadapt_int_cpks_coul_grad(i_grad,i_ir)%int = 0.0_r8_kind
                  enddo
      if(cpks_coul_grad_size.gt.max_cpks_coul_grad_size) then
        max_cpks_coul_grad_size=cpks_coul_grad_size
        DPRINT 'max_cpks_coul_grad_size :::', max_cpks_coul_grad_size
      endif
#endif
            endif splitg1
         endif grad_3cob1

         if (integralpar_2cob_nuc_grad) then
            do i_grad=1,gradient_data_n_gradients
               allocate( symadapt_int_2cob_nuc_grad(i_grad,i_ir)%&
                    int(dim_c2,dim_c1,n_if2,n_if1), stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: allocate of &
                    &symadapt_int_2cob_nuc_gr%int failed" )
               symadapt_int_2cob_nuc_grad(i_grad,i_ir)%int = 0.0_r8_kind
               if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
                  allocate( symadapt_int_2cob_pseudo_grad(i_grad,i_ir)%&
                       int(dim_c2,dim_c1,n_if2,n_if1), stat=status )
                  if ( status .ne. 0 ) call error_handler( &
                       "int_data_2cob3c_setup: allocate of &
                       &symadapt_int_2cob_pseudo_gr%int failed" )
                  symadapt_int_2cob_pseudo_grad(i_grad,i_ir)%int = 0.0_r8_kind
               endif
            end do
         end if

         if(integralpar_2cob_pc_grad) then
            do i_grad=1,totsym_grad_pc_length
               allocate( symadapt_int_2cob_pc_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_pc_grad(i_grad,i_ir)%int = 0.0_r8_kind
            end do
         end if

         if(integralpar_2cob_X_grad) then
            do i_grad=1,totsym_grad_dip_length
               allocate( symadapt_int_2cob_pd_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_pd_grad(i_grad,i_ir)%int = 0.0_r8_kind
               allocate( symadapt_int_2cob_pd_torq(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_pd_torq(i_grad,i_ir)%int = 0.0_r8_kind
            end do

            do i_grad=1,totsym_grad_quad_length
               allocate( symadapt_int_2cob_pq_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_pq_grad(i_grad,i_ir)%int = 0.0_r8_kind
               allocate( symadapt_int_2cob_pq_torq(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_pq_torq(i_grad,i_ir)%int = 0.0_r8_kind
            end do

            do i_grad=1,totsym_grad_oct_length
               allocate( symadapt_int_2cob_po_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_po_grad(i_grad,i_ir)%int = 0.0_r8_kind
               allocate( symadapt_int_2cob_po_torq(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_po_torq(i_grad,i_ir)%int = 0.0_r8_kind
            end do

            do i_grad=1,totsym_grad_rep_length
               allocate( symadapt_int_2cob_rc_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_rc_grad(i_grad,i_ir)%int = 0.0_r8_kind
            end do
         end if

         if(integralpar_2cob_ipd_grad) then
            do i_grad=1,totsym_grad_idip_length
               allocate( symadapt_int_2cob_ipd_grad(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_ipd_grad(i_grad,i_ir)%int = 0.0_r8_kind
               allocate( symadapt_int_2cob_ipd_torq(i_grad,i_ir)%&
                    int(n_c2,n_c1,n_if2,n_if1), stat=status )
               ASSERT(status==0)
               symadapt_int_2cob_ipd_torq(i_grad,i_ir)%int = 0.0_r8_kind
            end do
         end if

         if (integralpar_2cob_pvsp_grad) then
            do i_grad=1,gradient_data_n_gradients
               allocate( symadapt_int_2cob_pvsp_grad(i_grad,i_ir)%&
                    int(dim_c2,dim_c1,n_if2,n_if1), stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: allocate of &
                    &symadapt_int_2cob_pvsp_grd%int failed" )
               symadapt_int_2cob_pvsp_grad(i_grad,i_ir)%int = 0.0_r8_kind
            end do
         end if
         if ( integralpar_3c_xc ) then
            allocate( symadapt_int_3c_xc(i_ir)%int(n_c2,n_c1,n_ff_xc,n_if2,n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_3c_xc%int failed" )
            symadapt_int_3c_xc(i_ir)%int = 0.0_r8_kind
         endif
         if ( integralpar_3c_co ) then
            allocate( symadapt_int_3c_co(i_ir)%int(n_c2, n_c1, n_ff_ch, n_if2, n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_3c_co%int failed" )
            symadapt_int_3c_co(i_ir)%int = 0.0_r8_kind
         endif
#ifdef WITH_RESPONSE
         !! SB: allocation
         if ( integralpar_3c_co_resp ) then

            ! FIXME: dont special case on zeros:
            if (dimension_of_fit_ch(i_ir) == 0) cycle

            do i_ir_a = 1,n_ir
               n_if1_resp = ua1%symadapt_partner(i_ir_a, quadrupel%l1)%N_independent_fcts
               ! FIXME: dont special case on zeros:
               if (n_if1_resp == 0) cycle

               do i_ir_b = 1,n_ir
                  n_if2_resp = ua2%symadapt_partner(i_ir_b, quadrupel%l2)%N_independent_fcts
                  ! FIXME: dont special case on zeros:
                  if (n_if2_resp == 0) cycle

!!$                  print *,"data_2cob3c: ", i_ir, i_ir_a, i_ir_b
!!$                  print *,"data_2cob3c: checking... ", n_if1_resp, n_if2_resp, dimension_of_fit_ch(i_ir)

                  mult_irrep = cg(i_ir, i_ir_a, i_ir_b)%mult
!!$                  print *,"data_2cob3c: checking... mult = ", mult_irrep

                  ! FIXME: dont special case on zeros:
                  if (mult_irrep .eq. 0) cycle

!!$                  print *,"data_2cob3c: ", i_ir, i_ir_a, i_ir_b," was allocated"
                  allocate( symadapt_int_3c_co_resp(i_ir, i_ir_a, i_ir_b)%mult(mult_irrep), stat=status )

                  do imult = 1, mult_irrep
                     allocate( symadapt_int_3c_co_resp(i_ir, i_ir_a, i_ir_b)%mult(imult)%int(&
                          n_c2, n_c1, dimension_of_fit_ch(i_ir), n_if2_resp, n_if1_resp), stat=status )
                     if ( status .ne. 0 ) call error_handler( &
                          "int_data_2cob3c_setup: allocate of symadapt_int_3c_co_resp%int failed" )
!!$                     print *,"data_2cob3c: shape = ",shape(symadapt_int_3c_co_resp(i_ir,i_ir_a,i_ir_b)%mult(imult)%int)
                     symadapt_int_3c_co_resp(i_ir, i_ir_a, i_ir_b)%mult(imult)%int = 0.0_r8_kind
                  end do

               end do
            end do
         endif
#endif
         if ( integralpar_2cob_potential ) then
            allocate( symadapt_int_2cob_poten(i_ir)%int(n_c2,n_c1,N_points,n_if2,n_if1), &
                 stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_setup: allocate of symadapt_int_2cob_poten%int failed" )
            symadapt_int_2cob_poten(i_ir)%int = 0.0_r8_kind
         endif
         if ( integralpar_2cob_field) then
            allocate( symadapt_int_2cob_field(i_ir)%int(n_c2,n_c1,totsym_field_length,n_if2,n_if1), &
                 stat=status )
            ASSERT(status==0)
            symadapt_int_2cob_field(i_ir)%int = 0.0_r8_kind
         endif
      enddo
    end subroutine alloc_sa_int

    subroutine alloc_sa_int_spor()
      use spin_orbit_module, only: &
           is_on, &
           op_FitTrafo
      use error_module
      use symm_adapt_int, &
           sa_alloc => typealloc
      use prim_int_store, &
           prim_alloc => alloc
      implicit none
      ! *** end of interface ***

      integer(IK) :: memstat
      integer(IK) :: n_proj_ir,a2,a1,L2,L1
      integer(IK) :: ne2,ne1,nc2,nc1,nch2,nch1,nff_ch,nff_xc,nff_s,nff_r2
      integer(IK) :: i_sa
      integer(IK) :: &
           n_2c_ints, &   ! SA integrals, sum of two following (?)
           n_2c_pints, &  ! PRIM 2c integrals
           n_2cv_pints, & ! PRIM 2c vector integrals (pvxp)
           n_3c_ints, &   ! ... same for 3c
           n_3c_pints, &
           n_3cv_pints


      a2  = quadrupel%ua2
      a1  = quadrupel%ua1

      L2  = quadrupel%L2
      L1  = quadrupel%L1

      ne2 = n_exp2
      ne1 = n_exp1

      nc2 = n_c2
      nc1 = n_c1

      nff_ch  = n_ff_ch
      nff_s = n_ff_s
      nff_r2 = n_ff_r2
      nff_xc  = n_ff_xc

      n_proj_ir = symmetry_data_n_proj_irreps()

      if(is_on(op_FitTrafo))then
         nch2 = ne2
         nch1 = ne1
      else
         nch2 = nc2
         nch1 = nc1
      endif

      !=======================================================
      !
      ! ALLOCATE TWO CENTER INTEGRALS >>>
      !
      !-------------------------------------------------------

      n_2c_pints  = count((/ &
           integralpar_2cob_kin,  &
           integralpar_2cob_nuc, &
           integralpar_2cob_pvsp, &
           integralpar_2cob_ol &
           /))

      n_2cv_pints  = count((/ &
           integralpar_2cob_pvxp,&
           integralpar_2cob_pvec &
           /))

      n_2c_ints = n_2c_pints + n_2cv_pints
      DPRINT 'int: 2cp=',n_2c_pints
      DPRINT 'int: 2vcp=',n_2cv_pints
      DPRINT 'int: 2c=',n_2c_ints

      ! ==============================================
      ! Symmetry adapted ints:
      !
      ! (for SO) are ALL allocated here:
      allocate(sa_2c_ints(n_proj_ir, n_2c_ints), STAT=memstat)
      ASSERT(memstat==0)

      do i_sa=1,n_2c_ints
         call sa_alloc(a2, a1, L2, L1, ne2, ne1, sa_2c_ints(:, i_sa))
      enddo

      if(is_on(op_FitTrafo))then
         allocate(sa_alphap(n_proj_ir, 2), STAT=memstat)
         ASSERT(memstat==0)

         call sa_alloc(a2, a1, L2, L1, ne2, ne1, sa_alphap(:, 1), WHAT=LS_INT_BLOCK)
         call sa_alloc(a2, a1, L2, L1, ne2, ne1, sa_alphap(:, 2), WHAT=SL_INT_BLOCK)
      endif

      ! initialize index-references to particular ints,
      i_sa = 0
      if ( integralpar_2cob_kin ) then
         i_sa  = i_sa + 1
         int_sa_kin  = i_sa
         DPRINT 'int: sa_kin=',int_sa_kin
      endif
      if ( integralpar_2cob_nuc ) then
         i_sa  = i_sa + 1
         int_sa_nuc  = i_sa
         DPRINT 'int: sa_nuc=',int_sa_nuc
      endif
      if ( integralpar_2cob_pvsp ) then
         i_sa  = i_sa + 1
         int_sa_pvsp = i_sa
         DPRINT 'int: sa_pvsp=',int_sa_pvsp
      end if
      if ( integralpar_2cob_pvxp ) then
         i_sa  = i_sa + 1
         int_sa_pvxp = i_sa
         DPRINT 'int: sa_pvxp=',int_sa_pvxp
      endif
      if ( integralpar_2cob_pvec ) then
         i_sa = i_sa + 1
         int_sa_sigp = i_sa
         DPRINT 'int: sa_sigp=',int_sa_sigp
      endif
      if ( integralpar_2cob_ol ) then
         i_sa = i_sa + 1
         int_sa_ol    = i_sa
      endif
      ASSERT(i_sa==n_2c_ints)

      !=======================================================
      !
      ! ALLOCATE THREE CENTER INTEGRALS >>>
      !
      !-------------------------------------------------------

      n_3c_pints = count((/ &
           integralpar_3c_xc, &
           integralpar_3c_co, &
           integralpar_3c_rcoul_pvsp, &
           integralpar_3c_r2_pvsp &
           /))

      n_3cv_pints = count((/ &
           integralpar_3c_rcoul_pvxp, &
           integralpar_3c_r2_pvxp &
           /))

      n_3c_ints = n_3c_pints + n_3cv_pints

      DPRINT 'int: 3cp=',n_3c_pints
      DPRINT 'int: 3vcp=',n_3cv_pints
      DPRINT 'int: 3c=',n_3c_ints

      allocate(sa_3c_ints(n_proj_ir, n_3c_ints),STAT=memstat)
      call error(memstat,"idm/alloc_sa_int_spor: alloc 3c failed")

      i_sa = 0
      if ( integralpar_3c_xc ) then
         i_sa = i_sa + 1
         call sa_alloc(a2, a1, L2, L1, nc2, nc1, nff_xc, sa_3c_ints(:, i_sa))
         int_sa_xc    = i_sa
      endif
      if ( integralpar_3c_co ) then
         i_sa = i_sa + 1
         call sa_alloc(a2, a1, L2, L1, nch2, nch1, nff_ch, sa_3c_ints(:, i_sa))
         int_sa_co    = i_sa
      endif
      if ( integralpar_3c_rcoul_pvsp ) then
         i_sa = i_sa + 1
         call sa_alloc(a2, a1, L2, L1, nch2, nch1, nff_s, sa_3c_ints(:, i_sa))
         int_sa_rcoul_pvsp    = i_sa
      endif
      if ( integralpar_3c_rcoul_pvxp ) then
         i_sa = i_sa + 1
         call sa_alloc(a2, a1, L2, L1, nch2, nch1, nff_s, sa_3c_ints(:, i_sa))
         int_sa_rcoul_pvxp    = i_sa
      endif
      if ( integralpar_3c_r2_pvsp ) then
         i_sa = i_sa + 1
         call sa_alloc(a2, a1, L2, L1, nch2, nch1, nff_r2, sa_3c_ints(:, i_sa))
         int_sa_r2_pvsp    = i_sa
      endif
      if ( integralpar_3c_r2_pvxp ) then
         i_sa = i_sa + 1
         call sa_alloc(a2, a1, L2, L1, nch2, nch1, nff_r2, sa_3c_ints(:, i_sa))
         int_sa_r2_pvxp    = i_sa
      endif

    end subroutine alloc_sa_int_spor

  end subroutine int_data_2cob3c_setup
  !**************************************************************


  !**************************************************************
  subroutine int_data_2cob3c_shutdown()
    !  Purpose: shutdown work + deallocating
    !           (in particular symmetry-adapted integrals)
    !** End of interface ****************************************
    implicit none


!   print*,'call allocate_primitives shutdown'
    call allocate_primitives('deallocate')

    ! deallocate symmetry adapted integrals
    if (integralpar_spor) then
       !
       ! SPIN ORBIT
       !
       call free_sa_int_spor()
    else
       !
       ! STANDARD SCF (NO SPIN ORBIT)
       !
       call free_sa_int()
    endif

  contains

    subroutine free_sa_int()
#ifdef WITH_RESPONSE
      use clebsch_gordan, only: cg=>cg_eliminated
      use ch_response_module, only: dimension_of_fit_ch
#endif
      implicit none
      ! internal sub
      ! *** end of interface ***
      !------------ Declaration of local variables ----------------
      integer(IK) :: status, i_ir,i_grad
      integer(IK) :: k2dr
      logical     :: model_density, split_gradients

      !!SB: variables
      integer(IK) :: i_ir_a, i_ir_b, mult_irrep, imult
      integer(IK) :: n_if1_resp, n_if2_resp

      !------------ Executable code -------------------------------
      !
      ! STANDARD SCF (NO SPIN ORBIT)
      !

      split_gradients = options_split_gradients()
      model_density = options_xcmode() == xcmode_model_density .or. &
           options_xcmode() == xcmode_extended_mda


      do i_ir = 1, symmetry_data_n_irreps()
         if ( integralpar_2cob_kin ) then
            deallocate( symadapt_int_2cob_kin(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_kin%int failed" )
         endif
         if ( integralpar_2cob_nuc ) then
            deallocate( symadapt_int_2cob_nuc(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_nuc%int failed" )
            if((pseudopot_present.and.integralpar_relativistic) .or. &
                 integralpar_pseudo) &
                 deallocate( symadapt_int_2cob_nuc_pseudo(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_nuc_pseudo%int failed" )
         endif
         if ( integralpar_2cob_pvsp ) then
            deallocate( symadapt_int_2cob_pvsp(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_pvsp%int failed" )
         endif
         if ( integralpar_2cob_ol ) then
            deallocate( symadapt_int_2cob_ol(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_ol%int failed" )
         endif

         if ( integralpar_2cob_ol_grad ) then
            do i_grad=1,size(symadapt_int_2cob_ol_grad,1)
             deallocate( symadapt_int_2cob_ol_grad(i_grad,i_ir)%int,stat=status )
            ASSERT(status.eq.0)
          if(integralpar_dervs) then

            do k2dr=1,size(symadapt_int_2cob_ol_grad,1)
            ol_dervs_size=ol_dervs_size-size(symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int)
             MEMLOG(-size(symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int))
             deallocate(symadapt_int_2cob_ol_dervs(i_grad,k2dr,i_ir)%int, &   !!!free_sa
                                                        stat=cpksalloc(55))
            ASSERT(cpksalloc(55).eq.0)
            cpksalloc(55)=1
            enddo
          endif
            end do
         endif

         if ( integralpar_2cob_ol_rel_grad ) then
            do i_grad=1,size(symadapt_int_2cob_ol_rel_grad,1)
               deallocate( symadapt_int_2cob_ol_rel_grad(i_grad,i_ir)%int,&
                    & stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: deallocate of symadapt_int_2cob_ol_rel_grad%int failed" )
            end do
         endif
         if ( integralpar_2cob_kin_grad ) then
            do i_grad=1,size(symadapt_int_2cob_kin_grad,1)
               deallocate( symadapt_int_2cob_kin_grad(i_grad,i_ir)%int,&
                    & stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: deallocate of symadapt_int_2cob_kin_grad%int failed" )
            end do
         endif
         if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
            do i_grad=1,gradient_data_n_spin_gradients
               deallocate( symadapt_int_3cob_solv_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status.eq.0)
            end do
            if(with_pc .and. .not.fixed_pc) then
               do i_grad=1,totsym_field_length
                  deallocate( symadapt_int_3cob_solv_grad_pc(i_grad,i_ir)%int,&
                       & stat=status )
                  if ( status .ne. 0 ) call error_handler( &
                       "int_data_2cob3c_setup: deallocate of symadapt_int_3cob_solv_grad_pc%int failed" )
               end do
            end if
         endif

         if ( integralpar_3cob_grad ) then
#ifdef WITH_EPE
             if(ewpc_n.ne.0.and.calc_cluster_epe_energy)  then
              deallocate( symadapt_int_3cob_epe(1,i_ir)%int,stat=ewa_allocstat(8))
              ASSERT(ewa_allocstat(8).eq.0)
             endif
#endif

            do i_grad=1,gradient_data_n_spin_gradients
             deallocate( symadapt_int_3cob_grad(i_grad,i_ir)%int,stat=status)
               ASSERT(status.eq.0)



            if(integralpar_2dervs) then
             do k2dr=1,gradient_data_n_spin_gradients
               MEMLOG(-size(symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int))
               size_coul_dervs=size_coul_dervs-size(symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int)
               deallocate( symadapt_int_coul_dervs(i_grad,k2dr,i_ir)%int,stat=cpksalloc(38))
               ASSERT(cpksalloc(38).eq.0)
               cpksalloc(38)=1
             enddo
            endif
            enddo

#ifndef no_cpks_coul_grads
                if(integralpar_cpksdervs) then
                  do i_grad=1,size(symadapt_int_cpks_coul_grad,1)
                   MEMLOG(-size(symadapt_int_cpks_coul_grad(i_grad,i_ir)%int))
                   if(cpks_coul_grad_size.gt.max_cpks_coul_grad_size) max_cpks_coul_grad_size=cpks_coul_grad_size
                      cpks_coul_grad_size=cpks_coul_grad_size- &
                                                    size(symadapt_int_cpks_coul_grad(i_grad,i_ir)%int)
                   deallocate( symadapt_int_cpks_coul_grad(i_grad,i_ir)%int,stat=cpksalloc(23) )
                   ASSERT(cpksalloc(23).eq.0)
                   cpksalloc(23)=1
                  enddo
                 endif
#endif

            ! FIXME: mv out of if(integralpar_3cob_grad) ???
            if(integralpar_2dervs.and.integralpar_relativistic)then
            do i_grad=1,gradient_data_n_spin_gradients
              do k2dr=1,gradient_data_n_spin_gradients

                deallocate( symadapt_int_nucl_dervs(i_grad,k2dr,i_ir)%int,stat=status)
                ASSERT(status.eq.0)

                deallocate( symadapt_int_pvsp_dervs(i_grad,k2dr,i_ir)%int,stat=status)
                ASSERT(status.eq.0)
              enddo
            enddo
            do i_grad=1,size(symadapt_int_2cob_kin_dervs,1)
              do k2dr=1,size(symadapt_int_2cob_kin_dervs,2)

                deallocate( symadapt_int_2cob_olu_dervs(i_grad,k2dr,i_ir)%int,stat=status)
                ASSERT(status.eq.0)

                deallocate( symadapt_int_2cob_kin_dervs(i_grad,k2dr,i_ir)%int,stat=status)
                ASSERT(status.eq.0)
              enddo
            enddo
            endif

            if (split_gradients) then
               do i_grad=1,size(symadapt_int_2cob_ks_grad,1)
                  deallocate( symadapt_int_2cob_ks_grad(i_grad,i_ir)%int,&
                       & stat=status )
                  if (status /= 0) call error_handler(&
                       & "int_data_2cob3c_setup: deallocate of symadapt_int_2cob_ks_grad%int failed" )
               end do
               do i_grad=1,gradient_data_n_gradients
                  deallocate( symadapt_int_3cob_nuc_grad(i_grad,i_ir)%int,&
                       & stat=status )
                  if (status /= 0) call error_handler(&
                       & "int_data_2cob3c_setup: deallocate of symadapt_int_3cob_nuc_grad%int failed" )
               end do
               if (model_density) then
                  do i_grad=1,gradient_data_n_gradients
                     deallocate( symadapt_int_3cob_coul_grad(i_grad,i_ir)%int,&
                          & stat=status )
                     if (status/=0)call error_handler(&
                          & "int_data_2cob3c_setup: deallocate of symadapt_int_3cob_coul_grad%int failed")
                  end do
               endif
            endif
         endif
         if (integralpar_2cob_nuc_grad) then
            do i_grad=1,gradient_data_n_gradients
               deallocate( symadapt_int_2cob_nuc_grad(i_grad,i_ir)%int,&
                    & stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: deallocate of symadapt_int_2cob_nuc_grad%int failed" )
               if((pseudopot_present.and.integralpar_relativistic) .or.integralpar_pseudo) then
                  deallocate( symadapt_int_2cob_pseudo_grad(i_grad,i_ir)%int,&
                       & stat=status )
                  if ( status .ne. 0 ) call error_handler( &
                       "int_data_2cob3c_setup: deallocate of symadapt_int_2cob_pseudo_grad%int failed" )
               endif!
            end do

         end if

         if (integralpar_2cob_pc_grad) then
            do i_grad=1,totsym_grad_pc_length
               deallocate( symadapt_int_2cob_pc_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
            end do
         end if

         if (integralpar_2cob_X_grad) then
            do i_grad=1,totsym_grad_dip_length
               deallocate( symadapt_int_2cob_pd_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
               deallocate( symadapt_int_2cob_pd_torq(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
            end do

            do i_grad=1,totsym_grad_quad_length
               deallocate( symadapt_int_2cob_pq_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
               deallocate( symadapt_int_2cob_pq_torq(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
            end do

            do i_grad=1,totsym_grad_oct_length
               deallocate( symadapt_int_2cob_po_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
               deallocate( symadapt_int_2cob_po_torq(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
            end do

            do i_grad=1,totsym_grad_rep_length
               deallocate( symadapt_int_2cob_rc_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
            end do
         end if

         if (integralpar_2cob_ipd_grad) then
            do i_grad=1,totsym_grad_idip_length
               deallocate( symadapt_int_2cob_ipd_grad(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
               deallocate( symadapt_int_2cob_ipd_torq(i_grad,i_ir)%int,&
                    & stat=status )
               ASSERT(status==0)
            end do
         end if

         if (integralpar_2cob_pvsp_grad) then
            do i_grad=1,gradient_data_n_gradients
               deallocate( symadapt_int_2cob_pvsp_grad(i_grad,i_ir)%int,&
                    & stat=status )
               if ( status .ne. 0 ) call error_handler( &
                    "int_data_2cob3c_setup: deallocate of symadapt_int_2cob_pvsp_grad%int failed" )
            end do
         end if

         if ( integralpar_3c_xc ) then
            deallocate( symadapt_int_3c_xc(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_3c_xc%int failed" )
         endif
         if ( integralpar_3c_co ) then
            deallocate( symadapt_int_3c_co(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_3c_co%int failed" )
         endif
#ifdef WITH_RESPONSE
         !! SB: dealloc
         if ( integralpar_3c_co_resp ) then

            ! FIXME: dont special case on zeros:
            if (dimension_of_fit_ch(i_ir) == 0) cycle

            do i_ir_a = 1, symmetry_data_n_irreps()
               n_if1_resp = ua1%symadapt_partner(i_ir_a,quadrupel%l1)%N_independent_fcts
               if (n_if1_resp == 0) cycle
               do i_ir_b = 1, symmetry_data_n_irreps()
                  n_if2_resp = ua2%symadapt_partner(i_ir_b,quadrupel%l2)%N_independent_fcts
                  if (n_if2_resp == 0) cycle

                  mult_irrep = cg(i_ir,i_ir_a,i_ir_b)%mult

                  ! FIXME: dont special case on zeros:
                  if (mult_irrep .eq. 0) cycle

!!$                  print *,"int_data: sa_3c_co_resp dealloc:", i_ir, i_ir_a, i_ir_b
                  do imult = 1, mult_irrep

                     deallocate( symadapt_int_3c_co_resp(i_ir, i_ir_a, i_ir_b)%mult(imult)%int, stat=status )
                     if ( status .ne. 0 ) call error_handler( &
                          "int_data_2cob3c_shutdown: deallocate of symadapt_int_3c_co_resp%int failed" )
                  end do

                  deallocate( symadapt_int_3c_co_resp(i_ir, i_ir_a, i_ir_b)%mult, stat=status )
                  ASSERT(status == 0)

               end do
            end do
         endif
#endif
         if ( integralpar_2cob_potential ) then
            deallocate( symadapt_int_2cob_poten(i_ir)%int, stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_poten%int failed" )
         endif
         if ( integralpar_2cob_field ) then
            deallocate( symadapt_int_2cob_field(i_ir)%int, stat=status )
            ASSERT(status==0)
         endif
      enddo

      if ( integralpar_2cob_kin ) then
         deallocate( symadapt_int_2cob_kin, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_kin failed" )
      endif
      if ( integralpar_2cob_nuc ) then
         deallocate( symadapt_int_2cob_nuc, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_nuc failed" )
         if((pseudopot_present.and.integralpar_relativistic) .or. &
              integralpar_pseudo) &
              deallocate( symadapt_int_2cob_nuc_pseudo,stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_nuc_pseudo failed" )
      endif
      if ( integralpar_2cob_pvsp ) then
         deallocate( symadapt_int_2cob_pvsp, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_pvsp failed" )
      endif
      if ( integralpar_2cob_ol ) then
         deallocate( symadapt_int_2cob_ol, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_ol failed" )
      endif
      if ( integralpar_2cob_potential ) then
         deallocate( symadapt_int_2cob_poten, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_poten failed" )
      endif
      if ( integralpar_2cob_field ) then
         deallocate( symadapt_int_2cob_field, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_field failed" )
      endif
      if ( integralpar_3c_xc) then
         deallocate( symadapt_int_3c_xc, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_3c_xc failed" )
      endif
      if ( integralpar_3c_co) then
         deallocate( symadapt_int_3c_co, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_3c_co failed" )
      endif
      !! SB : dealloc
      if ( integralpar_3c_co_resp) then
         deallocate( symadapt_int_3c_co_resp, stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_3c_co_resp failed" )
      endif

      if ( integralpar_2cob_ol_grad ) then
         deallocate( symadapt_int_2cob_ol_grad,&
              prim_int_2cob_ol_grad &
             ,stat=status )
              ASSERT(status.eq.0)

       if(integralpar_dervs) then
        deallocate(prim_int_2cob_ol_dervs, symadapt_int_2cob_ol_dervs, &     !!! free_sa
                         stat=cpksalloc(45))
        ASSERT(cpksalloc(45).eq.0)
        cpksalloc(45)=1
       endif

      endif
      if ( integralpar_2cob_ol_rel_grad ) then
         deallocate( symadapt_int_2cob_ol_rel_grad &
              ,stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_ol_rel_grad failed" )
      endif
      if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
         DPRINT  '2cob3c: deallocate prim_int_3cob_solv_grad ...'
         deallocate( prim_int_3cob_solv_grad, &                            !!!!!!!!!!!!
              symadapt_int_3cob_solv_grad, stat=status )                   !!!!!!!!!!!!
         if (status /= 0) call error_handler(&
              & "int_data_2cob3c_shutdown: deallocate of int_3cob_solv_grad failed" )
         if(with_pc .and. .not.fixed_pc) then
            deallocate( prim_int_3cob_solv_grad_pc, &                            !!!!!!!!!!!!
                 symadapt_int_3cob_solv_grad_pc, stat=status )                   !!!!!!!!!!!!
            if (status /= 0) call error_handler(&
                 & "int_data_2cob3c_shutdown: deallocate of int_3cob_solv_grad_pc failed" )
         end if
      endif

      grad_3cob: if ( integralpar_3cob_grad ) then
#ifdef WITH_EPE
          if(ewpc_n.ne.0.and.calc_cluster_epe_energy) then
             deallocate(prim_int_3cob_epe &
                       ,symadapt_int_3cob_epe, &
                         stat=ewa_allocstat(5))
            if (ewa_allocstat(5).ne.0) call error_handler(&
              & "int_data_2cob3c_shutdown: deallocate of prim_int_3cob_epe failed" )
          endif
#endif
         deallocate( symadapt_int_3cob_grad, prim_int_3cob_grad,stat=status )
          ASSERT(status.eq.0)

          if(integralpar_2dervs) then
            deallocate(prim_int_coul_dervs,symadapt_int_coul_dervs,stat=cpksalloc(36) )
             ASSERT(cpksalloc(36).eq.0)
             cpksalloc(36)=1
          endif

          ! FIXME: mv out of if(integralpar_3cob_grad) ???
          if(integralpar_2dervs.and.integralpar_relativistic) then
            deallocate(    prim_int_nucl_dervs &
                      ,symadapt_int_nucl_dervs &
                      ,stat=status )
            ASSERT(status.eq.0)
            deallocate(    prim_int_pvsp_dervs &
                      ,symadapt_int_pvsp_dervs &
                      ,stat=status )
            ASSERT(status.eq.0)
            deallocate(    prim_int_2cob_kin_dervs &
                      ,symadapt_int_2cob_kin_dervs &
                      ,symadapt_int_2cob_olu_dervs &
                      ,stat=status )
            ASSERT(status.eq.0)
          endif

         splitg: if (split_gradients) then
            deallocate( prim_int_2cob_ks_grad, &
                 symadapt_int_2cob_ks_grad, stat=status )
            if (status /= 0) call error_handler(&
                 & "int_data_2cob3c_shutdown: deallocate of int_2cob_ks_grad failed" )
            deallocate( prim_int_3cob_nuc_grad &
                       ,symadapt_int_3cob_nuc_grad, stat=status )
            ASSERT(status.eq.0)
            if (model_density) then
               deallocate( prim_int_3cob_coul_grad &
                          ,symadapt_int_3cob_coul_grad, stat=status )
                    ASSERT(status.eq.0)
            endif

#ifndef no_cpks_coul_grads
         elseif(integralpar_cpksdervs) then
          deallocate(prim_int_cpks_coul_grad,symadapt_int_cpks_coul_grad,stat=cpksalloc(17))
          ASSERT(cpksalloc(17).eq.0)
          cpksalloc(17)=1
#endif
         endif splitg
      endif grad_3cob

      if ( integralpar_2cob_pvsp_grad ) then
         deallocate( symadapt_int_2cob_pvsp_grad,&
              prim_int_2cob_pvsp_grad, stat=status )
              ASSERT(status.eq.0)
      endif
      if ( integralpar_2cob_nuc_grad ) then
         deallocate( symadapt_int_2cob_nuc_grad,&
              prim_int_2cob_nuc_grad,stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_nuc_grad failed" )
         if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
            deallocate( symadapt_int_2cob_pseudo_grad, prim_int_2cob_pseudo_grad,stat=status )
            if ( status .ne. 0 ) call error_handler( &
                 "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_pseudo_grad failed" )
         endif! pseudopot_present
      endif

      if ( integralpar_2cob_pc_grad ) then
         deallocate( symadapt_int_2cob_pc_grad,&
              prim_int_2cob_pc_grad,stat=status )
         ASSERT(status==0)
      end if

      if ( integralpar_2cob_X_grad ) then
         deallocate( symadapt_int_2cob_pd_grad,&
              prim_int_2cob_pd_grad,stat=status )
         ASSERT(status==0)
         deallocate( symadapt_int_2cob_pd_torq,&
              prim_int_2cob_pd_torq,stat=status )
         ASSERT(status==0)

         deallocate( symadapt_int_2cob_pq_grad,&
              prim_int_2cob_pq_grad,stat=status )
         ASSERT(status==0)
         deallocate( symadapt_int_2cob_pq_torq,&
              prim_int_2cob_pq_torq,stat=status )
         ASSERT(status==0)

         deallocate( symadapt_int_2cob_po_grad,&
              prim_int_2cob_po_grad,stat=status )
         ASSERT(status==0)
         deallocate( symadapt_int_2cob_po_torq,&
              prim_int_2cob_po_torq,stat=status )
         ASSERT(status==0)

         deallocate( symadapt_int_2cob_rc_grad,&
              prim_int_2cob_rc_grad,stat=status )
         ASSERT(status==0)
      end if

      if ( integralpar_2cob_ipd_grad ) then
         deallocate( symadapt_int_2cob_ipd_grad,&
              prim_int_2cob_ipd_grad,stat=status )
         ASSERT(status==0)
         deallocate( symadapt_int_2cob_ipd_torq,&
              prim_int_2cob_ipd_torq,stat=status )
         ASSERT(status==0)
      end if

      if ( integralpar_2cob_kin_grad ) then
         deallocate( symadapt_int_2cob_kin_grad,&
              prim_int_2cob_kin_grad,stat=status )
         if ( status .ne. 0 ) call error_handler( &
              "int_data_2cob3c_shutdown: deallocate of symadapt_int_2cob_kin_grad failed" )
      endif
    end subroutine free_sa_int

    subroutine free_sa_int_spor()
      use error_module
      use spin_orbit_module, only: &
           is_on, &
           op_FitTrafo
      use symm_adapt_int, sa_dealloc => typedealloc
      use prim_int_store, prim_free => dealloc
      implicit none
      ! *** end of interface ***

      integer(IK) :: i, j, memstat

      !=======================================================
      !
      ! DEALLOCATE TWO-CENTER INTEGRALS >>>
      !

      !------------------------------------
      ! Symmetry adapted:
      do j = 1, size(sa_2c_ints, 2)
        do i = 1, size(sa_2c_ints, 1)
          call sa_dealloc(sa_2c_ints(i, j))
        enddo
      enddo

      deallocate(sa_2c_ints, STAT=memstat)
      ASSERT(memstat==0)

      if(is_on(op_FitTrafo))then
         do j = 1, size(sa_alphap, 2)
           do i = 1, size(sa_alphap, 1)
             call sa_dealloc(sa_alphap(i, j))
           enddo
         enddo

         deallocate(sa_alphap, STAT=memstat)
         ASSERT(memstat==0)
      endif

      int_sa_kin  = -1
      int_sa_nuc  = -1
      int_sa_pvsp = -1
      int_sa_pvxp = -1
      int_sa_sigp = -1
      int_sa_ol   = -1

      !=======================================================
      !
      ! DEALLOCATE THREE-CENTER INTEGRALS >>>
      !

      do j = 1, size(sa_3c_ints, 2)
        do i = 1, size(sa_3c_ints, 1)
          call sa_dealloc(sa_3c_ints(i, j))
        enddo
      enddo

      deallocate(sa_3c_ints, STAT=memstat)
      ASSERT(memstat==0)

      int_sa_xc = -1
      int_sa_co = -1
      int_sa_rcoul_pvsp = -1
      int_sa_rcoul_pvxp = -1
      int_sa_r2_pvsp = -1
      int_sa_r2_pvxp = -1

    end subroutine free_sa_int_spor

  end subroutine int_data_2cob3c_shutdown
  !**************************************************************

  !**************************************************************
  subroutine allocate_primitives(do_what)
    ! Purpose: allocates/deallocates storage for primitive integrals
    use symmetry_data_module, only: symmetry_data_n_irreps,  &
         symmetry_data_n_partners,&
         get_totalsymmetric_irrep
#ifdef WITH_RESPONSE
    use ch_response_module, only: dimension_of_fit_ch
#endif
    use prim_int_store, &
         prim_alloc => alloc, &
         prim_free => dealloc
    use spin_orbit_module, only: &
         is_on, &
         op_FinNuc
    use unique_atom_module, only: &
         n_ua => n_unique_atoms
    use cpksdervs_matrices,only: cpksalloc
    implicit none
    character(len=*), intent(in) :: do_what
    ! *** end of interface ***

    !------------ Declaration of local variables ----------------
    integer(kind=i4_kind) :: status,i_grad,k2dr
    logical :: alloc_, dealloc_, init_
    integer(i4_kind) :: inc
    integer(i4_kind) :: memstat
    integer(i4_kind) :: &
         i_pi, i_p2c, i_p2cv, i_p3c, i_p3cv, &
         n_2c_pints, &  ! PRIM 2c integrals
         n_2cv_pints, & ! PRIM 2c vector integrals (pvxp)
         n_3c_pints, &
         n_3cv_pints
    integer(i4_kind) :: n_dim_ch, i_ir
    !------------ Executable code ------------------------------

    select case ( do_what )
    case ( 'allocate' )
       alloc_   = .true.
       init_    = .true.
       dealloc_ = .false.
       inc      = +1
    case ( 'deallocate' )
       alloc_   = .false.
       init_    = .false.
       dealloc_ = .true.
       inc      = -1
    case ( 'initialize' )
       alloc_   = .false.
       init_    = .true.
       dealloc_ = .false.
       inc      =  0
    case default
       print *,'allocate_primitives: do_what = >',do_what,'<'
       ABORT('no such case')
    end select


    DPRINT '###(',do_what,')',inc,quadrupel%ua1,quadrupel%l1,quadrupel%ua2,quadrupel%l2

    n_2c_pints  = count((/ &
         integralpar_2cob_kin,  &
         integralpar_2cob_nuc, &
         integralpar_2cob_pvsp, &
         integralpar_2cob_ol &
         /))

    n_2cv_pints  = count((/ &
         integralpar_2cob_pvxp,&
         integralpar_2cob_pvec &
         /))

    n_3c_pints = count((/ &
         integralpar_3c_xc, &
         integralpar_3c_co, &
         integralpar_3c_rcoul_pvsp, &
         integralpar_3c_r2_pvsp &
         /))

    n_3cv_pints = count((/ &
         integralpar_3c_rcoul_pvxp, &
         integralpar_3c_r2_pvxp &
         /))

    if( is_on(op_FinNuc) )then
       ! In total:
       !      1(3) + 1 = 2(4) scalar 3c(ua) ints
       ! and  1    + 1 = 2    vector 3c(ua) ints

       ! + split PVSP(ua), V(ua), V_{fin}(ua) ints
       ! in one three-fold aray (check the dims!):
       n_3c_pints  = n_3c_pints  + 1 ! (3)

       ! + split PV_{fin}SP(ua) ints:
       n_3c_pints  = n_3c_pints  + 1

       ! + split PVXP(ua) ints:
       n_3cv_pints = n_3cv_pints + 1

       ! + split PV_{fin}XP(ua) ints:
       n_3cv_pints = n_3cv_pints + 1

       ! + merged 2c V_{???} ints:
       n_2c_pints  = n_2c_pints  + 1

       ! + merged 2c PV_{???}SP ints:
       n_2c_pints  = n_2c_pints  + 1

       ! + merged 2c PV_{???}XP ints:
       n_2cv_pints = n_2cv_pints + 1
    endif

    DPRINT 'ap: 2cp =',n_2c_pints
    DPRINT 'ap: 2cvp=',n_2cv_pints
    DPRINT 'ap: 3cp =',n_3c_pints
    DPRINT 'ap: 3cvp=',n_3cv_pints


    !---------------------------------------
    ! 2c scalar-type integrals:
    if( alloc_ )then
       allocate(prim_2c_ints(n_2c_pints),STAT=memstat)
       ASSERT(memstat==0)
       do i_pi=1,n_2c_pints
          ! all of same dimensions:
          call prim_alloc(n_exp2,n_exp1,n_m2,n_m1,prim_2c_ints(i_pi))
       end do
    endif

    if( dealloc_ )then
       do i_pi=1,size(prim_2c_ints)
          call prim_free(prim_2c_ints(i_pi))
       end do
       deallocate(prim_2c_ints,STAT=memstat)
       ASSERT(memstat==0)
    endif


    !---------------------------------------
    ! 2c vector-type integrals:
    if( alloc_ )then
       ! prim_vec_ints are allocated here:
       allocate(prim_2cv_ints(n_2cv_pints),STAT=memstat)
       ASSERT(memstat==0)
       do i_pi=1,n_2cv_pints
          ! all of same dimensions:
          call prim_alloc(n_exp2,n_exp1,n_m2,n_m1,prim_2cv_ints(i_pi))
       end do
    endif

    if( dealloc_ )then
       ! prim_vec_ints are deallocated here:
       do i_pi=1,size(prim_2cv_ints)
          call prim_free(prim_2cv_ints(i_pi))
       end do
       deallocate(prim_2cv_ints,STAT=memstat)
       ASSERT(memstat==0)
    endif


    !---------------------------------------
    ! 3c scalar-type integrals:
    if( alloc_ )then
       ! prim_3c_ints are allocated here:
       allocate(prim_3c_ints(n_3c_pints),STAT=memstat)
       ASSERT(memstat==0)
       ! 3c are of different dimensions, allocate storage later
    endif

    if( dealloc_ )then
       ! prim_3c_ints are deallocated here:
       do i_pi=1,size(prim_3c_ints)
          call prim_free(prim_3c_ints(i_pi))
       end do
       deallocate(prim_3c_ints,STAT=memstat)
       ASSERT(memstat==0)
    endif

    !---------------------------------------
    ! 3C vector-type integrals:
    if( alloc_ )then
       ! prim_3cv_ints are allocated here:
       allocate(prim_3cv_ints(n_3cv_pints),STAT=memstat)
       ASSERT(memstat==0)
       ! 3c are of different dimensions, allocate storage later
    endif

    if( dealloc_ )then
       ! prim_3c_ints are deallocated here:
       do i_pi=1,size(prim_3cv_ints)
          call prim_free(prim_3cv_ints(i_pi))
       end do
       deallocate(prim_3cv_ints,STAT=memstat)
       ASSERT(memstat==0)
    endif

    ! for compatibility set the pointers ...
    i_p2c  = 0
    i_p2cv = 0
    i_p3c  = 0
    i_p3cv = 0
    ! ... of (1) 2c-scalar-type ints:
    if ( integralpar_2cob_kin ) then
       i_p2c = i_p2c + 1
       int_prim_kin = i_p2c
       if( alloc_ ) prim_int_2cob_kin => prim_2c_ints(int_prim_kin)%prim
       if( dealloc_ ) int_prim_kin = -1
    endif

    if ( integralpar_2cob_nuc ) then
       i_p2c = i_p2c + 1
       int_prim_nuc = i_p2c
       if( alloc_ ) prim_int_2cob_nuc => prim_2c_ints(int_prim_nuc)%prim
       if( dealloc_ ) int_prim_nuc = -1

       ! or do the "old" way ...
       if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
          if( alloc_ )then
             allocate( prim_int_2cob_nuc_pseudo(n_exp2,n_exp1,n_m2,n_m1), &
                  stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_nuc_pseudo))

          if( dealloc_ )then
             deallocate( prim_int_2cob_nuc_pseudo, &
                  stat=status )
             ASSERT(status==0)
          endif
       endif
    endif

    if ( integralpar_2cob_ol ) then
       i_p2c = i_p2c + 1
       int_prim_ol = i_p2c
       if( alloc_ ) prim_int_2cob_ol => prim_2c_ints(i_p2c)%prim
       if( dealloc_ ) int_prim_ol = -1
    endif

    if ( integralpar_2cob_pvsp ) then
       i_p2c = i_p2c + 1
       int_prim_pvsp = i_p2c
       if( alloc_ ) prim_int_2cob_pvsp => prim_2c_ints(i_p2c)%prim
       if( dealloc_ ) int_prim_pvsp = -1
    endif

    if( is_on(op_FinNuc) )then
       i_p2c = i_p2c + 1
       int_prim_nucfin = i_p2c
       if( dealloc_ ) int_prim_nucfin = -1

       i_p2c = i_p2c + 1
       int_prim_pfinsp = i_p2c
       if( dealloc_ ) int_prim_pfinsp = -1
    endif

    ! ... of (2) 2c-vector-type ints:
    if ( integralpar_2cob_pvxp ) then
       i_p2cv = i_p2cv + 1
       int_prim_pvxp = i_p2cv
       if( dealloc_ ) int_prim_pvxp = -1
    endif

    if ( integralpar_2cob_pvec ) then
       i_p2cv = i_p2cv + 1
       int_prim_pvec = i_p2cv
       if( dealloc_ ) int_prim_pvec = -1
    endif

    if ( is_on(op_FinNuc) ) then
       i_p2cv = i_p2cv + 1
       int_prim_pfinxp = i_p2cv
       if( dealloc_ ) int_prim_pfinxp = -1
    endif

    ! ... of (3) 3c-scalar-type ints:
    if ( integralpar_3c_xc ) then
       i_p3c = i_p3c + 1
       int_prim_xc = i_p3c
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,fit_coeff_n_xc(),n_m2,n_m1, &
               prim_3c_ints(i_p3c) &
               )
          prim_int_3c_xc => prim_3c_ints(i_p3c)%prim
       endif
       if( dealloc_ ) int_prim_xc = -1
       ! storage already deallocated in a loop above
    endif

#ifdef WITH_RESPONSE
    if ( integralpar_3c_co_resp ) then
       i_p3c = i_p3c + 1
       int_prim_co = i_p3c

       n_dim_ch = 0
       do i_ir = 1,symmetry_data_n_irreps()
          n_dim_ch = n_dim_ch+dimension_of_fit_ch(i_ir)*&
               symmetry_data_n_partners(i_ir)
       end do

       if( alloc_ )then
          !
          ! Here n_dim_ch is the total number of fit functions including
          ! those that are not totally symmetric.
          !
          call prim_alloc(n_exp2, n_exp1, n_dim_ch, n_m2, n_m1, &
               prim_3c_ints(i_p3c) &
               )
          prim_int_3c_co => prim_3c_ints(i_p3c)%prim
       endif
       if( dealloc_ ) int_prim_co = -1

       ! storage already deallocated in a loop above
       !!SB added
    endif
#endif

    !
    ! FIXME: why? Possibly because the storage for the primitive integrals
    !        including those with the totally symmetric fit functions was
    !        already reserved above.
    !
    if ( integralpar_3c_co .and. (.not. integralpar_3c_co_resp) ) then
       i_p3c = i_p3c + 1
       int_prim_co = i_p3c
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,fit_coeff_n_ch(),n_m2,n_m1, &
               prim_3c_ints(i_p3c) &
               )
          prim_int_3c_co => prim_3c_ints(i_p3c)%prim
       endif
       if( dealloc_ ) int_prim_co = -1
       ! storage already deallocated in a loop above
    endif

    if ( integralpar_3c_rcoul_pvsp ) then
       i_p3c = i_p3c + 1
       int_prim_rcoul_pvsp = i_p3c
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,fit_coeff_n_s(),n_m2,n_m1, &
               prim_3c_ints(i_p3c) &
               )
       endif
       if( dealloc_ ) int_prim_rcoul_pvsp = -1
       ! storage already deallocated in a loop above
    endif

    if ( integralpar_3c_r2_pvsp ) then
       i_p3c = i_p3c + 1
       int_prim_r2_pvsp = i_p3c
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,fit_coeff_n_r2(),n_m2,n_m1, &
               prim_3c_ints(i_p3c) &
               )
       endif
       if( dealloc_ ) int_prim_r2_pvsp = -1
       ! storage already deallocated in a loop above
    endif

    if( is_on(op_FinNuc) )then
       i_p3c = i_p3c + 1
       int_prim_many_3c = i_p3c
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,OFF_STRIDE*n_ua,n_m2,n_m1, &
               prim_3c_ints(i_p3c) &
               )
       endif
       if( dealloc_ ) int_prim_many_3c = -1
       ! storage already deallocated in a loop above

       i_p3c = i_p3c + 1
       int_prim_split_pfinsp = i_p3c
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,n_ua,n_m2,n_m1, &
               prim_3c_ints(i_p3c) &
               )
       endif
       if( dealloc_ ) int_prim_split_pfinsp = -1
       ! storage already deallocated in a loop above
    endif

    ! ... of (4) 3c-vector-type ints:
    if ( integralpar_3c_rcoul_pvxp ) then
       i_p3cv = i_p3cv + 1
       int_prim_rcoul_pvxp = i_p3cv
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,fit_coeff_n_s(),n_m2,n_m1, &
               prim_3cv_ints(i_p3cv) &
               )
       endif
       if( dealloc_ ) int_prim_rcoul_pvxp = -1
       ! storage already deallocated in a loop above
    endif

    if ( integralpar_3c_r2_pvxp ) then
       i_p3cv = i_p3cv + 1
       int_prim_r2_pvxp = i_p3cv
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,fit_coeff_n_r2(),n_m2,n_m1, &
               prim_3cv_ints(i_p3cv) &
               )
       endif
       if( dealloc_ ) int_prim_r2_pvxp = -1
       ! storage already deallocated in a loop above
    endif

    if( is_on(op_FinNuc) )then
       i_p3cv = i_p3cv + 1
       int_prim_split_pvxp = i_p3cv
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,n_ua,n_m2,n_m1, &
               prim_3cv_ints(i_p3cv) &
               )
       endif
       if( dealloc_ ) int_prim_split_pvxp = -1
       ! storage already deallocated in a loop above

       i_p3cv = i_p3cv + 1
       int_prim_split_pfinxp = i_p3cv
       if( alloc_ )then
          call prim_alloc(n_exp2,n_exp1,n_ua,n_m2,n_m1, &
               prim_3cv_ints(i_p3cv) &
               )
       endif
       if( dealloc_ ) int_prim_split_pfinxp = -1
       ! storage already deallocated in a loop above
    endif

    ! ... and do the rest traditionally:
    if ( integralpar_2cob_potential ) then
       if( alloc_ )then
          allocate( prim_int_2cob_poten(n_exp2,n_exp1,N_points,n_m2,n_m1), &
               stat=status )
          ASSERT(status==0)
       endif

       MEMLOG(inc*size(prim_int_2cob_poten))

       if( dealloc_ )then
          deallocate( prim_int_2cob_poten, &
               stat=status )
          ASSERT(status==0)
       endif
    endif

    if ( integralpar_2cob_field )then
       if( alloc_ )then
          allocate( prim_int_2cob_field(n_exp2,n_exp1,totsym_field_length,n_m2,n_m1), &
                  stat=status )
          ASSERT(status==0)
       endif

       MEMLOG(inc*size(prim_int_2cob_field))

       if( dealloc_ )then
          deallocate( prim_int_2cob_field, &
               stat=status )
          ASSERT(status==0)
       endif
    endif

    if ( integralpar_2cob_ol_grad ) then
       do i_grad=1,size(prim_int_2cob_ol_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_ol_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_ol_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_ol_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       if(integralpar_dervs) then
          do i_grad=1,size(prim_int_2cob_ol_dervs,1)
             do k2dr=1,size(prim_int_2cob_ol_dervs,2)
                if( alloc_ )then
                   allocate( prim_int_2cob_ol_dervs(i_grad,k2dr)%m(n_exp2,n_exp1,n_m2,n_m1)&  !alloc prim
                           ,stat=cpksalloc(158) )
                   ASSERT(cpksalloc(158).eq.0)
                endif

                MEMLOG(inc*size(prim_int_2cob_ol_dervs(i_grad,k2dr)%m))
                ol_dervs_size=ol_dervs_size+inc*size(prim_int_2cob_ol_dervs(i_grad,k2dr)%m)

                if( dealloc_ )then
                   deallocate( prim_int_2cob_ol_dervs(i_grad,k2dr)%m,stat=cpksalloc(158) )
                   ASSERT(cpksalloc(158).eq.0)
                   cpksalloc(158)=1
                endif
             enddo
          enddo
       if(ol_dervs_size.gt.max_ol_dervs_size) then
        max_ol_dervs_size=ol_dervs_size
        DPRINT 'max_ol_dervs_size ',max_ol_dervs_size
       endif
       endif
    endif

    if ( integralpar_2cob_kin_grad ) then
       do i_grad=1,size(prim_int_2cob_kin_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_kin_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_kin_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_kin_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo
    endif

    if ( integralpar_2cob_nuc_grad ) then
       do i_grad=1,size(prim_int_2cob_nuc_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_nuc_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_nuc_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_nuc_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
          do i_grad=1,size(prim_int_2cob_pseudo_grad)
             if( alloc_ )then
                allocate( prim_int_2cob_pseudo_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                     ,stat=status )
                ASSERT(status==0)
             endif

             MEMLOG(inc*size(prim_int_2cob_pseudo_grad(i_grad)%m))

             if( dealloc_ )then
                deallocate( prim_int_2cob_pseudo_grad(i_grad)%m&
                     ,stat=status )
                ASSERT(status==0)
             endif
          enddo
       endif !pseudopot_present
    endif ! integralpar_2cob_nuc_grad

    if ( integralpar_2cob_pc_grad ) then
       do i_grad=1,size(prim_int_2cob_pc_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_pc_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_pc_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_pc_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo
    end if

    if ( integralpar_2cob_X_grad ) then
       do i_grad=1,size(prim_int_2cob_pd_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_pd_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
             allocate( prim_int_2cob_pd_torq(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_pd_grad(i_grad)%m))
          MEMLOG(inc*size(prim_int_2cob_pd_torq(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_pd_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
             deallocate( prim_int_2cob_pd_torq(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       do i_grad=1,size(prim_int_2cob_pq_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_pq_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
             allocate( prim_int_2cob_pq_torq(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_pq_grad(i_grad)%m))
          MEMLOG(inc*size(prim_int_2cob_pq_torq(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_pq_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
             deallocate( prim_int_2cob_pq_torq(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       do i_grad=1,size(prim_int_2cob_po_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_po_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
             allocate( prim_int_2cob_po_torq(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_po_grad(i_grad)%m))
          MEMLOG(inc*size(prim_int_2cob_po_torq(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_po_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
             deallocate( prim_int_2cob_po_torq(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       do i_grad=1,size(prim_int_2cob_rc_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_rc_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_rc_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_rc_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo
    end if

    if ( integralpar_2cob_ipd_grad ) then
       do i_grad=1,size(prim_int_2cob_ipd_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_ipd_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
             allocate( prim_int_2cob_ipd_torq(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_ipd_grad(i_grad)%m))
          MEMLOG(inc*size(prim_int_2cob_ipd_torq(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_ipd_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
             deallocate( prim_int_2cob_ipd_torq(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo
    end if

    if ( integralpar_2cob_pvsp_grad ) then
       do i_grad=1,size(prim_int_2cob_pvsp_grad)
          if( alloc_ )then
             allocate( prim_int_2cob_pvsp_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_2cob_pvsp_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_2cob_pvsp_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo
    endif

    if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
       do i_grad=1,gradient_data_n_spin_gradients
          if( alloc_ )then
             allocate( prim_int_3cob_solv_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_3cob_solv_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_3cob_solv_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       if (with_pc .and. .not.fixed_pc) then
          do i_grad=1,totsym_field_length
             if( alloc_ )then
                allocate( prim_int_3cob_solv_grad_pc(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                     ,stat=status )
                ASSERT(status==0)
             endif

             ! FIXME: move init to init_storage():
             if ( init_ ) prim_int_3cob_solv_grad_pc(i_grad)%m = 0.0_r8_kind

             MEMLOG(inc*size(prim_int_3cob_solv_grad_pc(i_grad)%m))

             if( dealloc_ )then
                deallocate( prim_int_3cob_solv_grad_pc(i_grad)%m&
                     ,stat=status )
                ASSERT(status==0)
             endif
          enddo
       end if
    endif

    if ( integralpar_3cob_grad ) then

#ifdef WITH_EPE
       if(ewpc_n.ne.0.and.calc_cluster_epe_energy) then
          if( alloc_ )then
             allocate(prim_int_3cob_epe%m(n_exp2,n_exp1,n_m2,n_m1),stat=ewa_allocstat(6))
             ASSERT(ewa_allocstat(6)==0)
             ewa_allocstat(6)=1
          endif

          MEMLOG(inc*size(prim_int_3cob_epe%m))

          if( dealloc_ )then
             deallocate(prim_int_3cob_epe%m,stat=ewa_allocstat(6))
             ASSERT(ewa_allocstat(6)==0)
          endif
       endif
#endif

       do i_grad=1,gradient_data_n_spin_gradients
          if( alloc_ )then
             allocate( prim_int_3cob_grad(i_grad)%m(n_exp2,n_exp1,n_m2,n_m1)&
                  ,stat=status )
             ASSERT(status==0)
          endif

          MEMLOG(inc*size(prim_int_3cob_grad(i_grad)%m))

          if( dealloc_ )then
             deallocate( prim_int_3cob_grad(i_grad)%m&
                  ,stat=status )
             ASSERT(status==0)
          endif
       enddo

       if(integralpar_2dervs) then
          do i_grad=1,gradient_data_n_spin_gradients
            do k2dr  =1,gradient_data_n_spin_gradients
               if( alloc_ )then
                  allocate( prim_int_coul_dervs(i_grad,k2dr)%m(n_exp2,n_exp1,n_m2,n_m1)&
                     ,stat=cpksalloc(35) )
                    ASSERT(cpksalloc(35).eq.0)
                endif

                MEMLOG(inc*size(prim_int_coul_dervs(i_grad,k2dr)%m))
                size_coul_dervs=size_coul_dervs+inc*size(prim_int_coul_dervs(i_grad,k2dr)%m)

                if( dealloc_ )then
                   deallocate(prim_int_coul_dervs(i_grad,k2dr)%m,stat=cpksalloc(35))
                   ASSERT(cpksalloc(35).eq.0)
                   cpksalloc(35)=1
                endif
            enddo
          enddo
         if(size_coul_dervs.gt.max_size_coul_dervs) then
          max_size_coul_dervs=size_coul_dervs
          DPRINT 'int_data_2cob3c_setup: max_size_coul_dervs ',max_size_coul_dervs
         endif
       endif

       if(integralpar_2dervs.and.integralpar_relativistic) then
          do i_grad=1,gradient_data_n_spin_gradients
            do k2dr  =1,gradient_data_n_spin_gradients
               if( alloc_ )then
                  allocate( prim_int_nucl_dervs(i_grad,k2dr)%m(n_exp2,n_exp1,n_m2,n_m1)&
                     ,stat=status )
                    ASSERT(status.eq.0)
                  allocate( prim_int_pvsp_dervs(i_grad,k2dr)%m(n_exp2,n_exp1,n_m2,n_m1)&
                     ,stat=status )
                    ASSERT(status.eq.0)
                endif

                MEMLOG(2*inc*size(prim_int_nucl_dervs(i_grad,k2dr)%m))

                if( dealloc_ )then
                   deallocate(prim_int_nucl_dervs(i_grad,k2dr)%m,stat=status)
                   ASSERT(status.eq.0)
                   deallocate(prim_int_pvsp_dervs(i_grad,k2dr)%m,stat=status)
                   ASSERT(status.eq.0)
                endif
            enddo
          enddo
          do i_grad=1,size(prim_int_2cob_kin_dervs,1)
            do k2dr  =1,size(prim_int_2cob_kin_dervs,2)
               if( alloc_ )then
                  allocate( prim_int_2cob_kin_dervs(i_grad,k2dr)%m(n_exp2,n_exp1,n_m2,n_m1)&
                     ,stat=status )
                    ASSERT(status.eq.0)
                endif

                MEMLOG(1*inc*size(prim_int_2cob_kin_dervs(i_grad,k2dr)%m))

                if( dealloc_ )then
                   deallocate(prim_int_2cob_kin_dervs(i_grad,k2dr)%m,stat=status)
                   ASSERT(status.eq.0)
                endif
            enddo
          enddo
       endif

       splitg: if (options_split_gradients()) then
          do i_grad=1,size(prim_int_2cob_ks_grad)
             if( alloc_ )then
                allocate( prim_int_2cob_ks_grad(i_grad)%&
                     &m(n_exp2,n_exp1,n_m2,n_m1), stat=status )
                ASSERT(status==0)
             endif

             MEMLOG(inc*size(prim_int_2cob_ks_grad(i_grad)%m))

             if( dealloc_ )then
                deallocate( prim_int_2cob_ks_grad(i_grad)%&
                     &m, stat=status )
                ASSERT(status==0)
             endif
          enddo

          do i_grad=1,gradient_data_n_gradients
             if( alloc_ )then
                allocate( prim_int_3cob_nuc_grad(i_grad)%&
                     &m(n_exp2,n_exp1,n_m2,n_m1) , stat=status )
                ASSERT(status==0)
             endif

             MEMLOG(inc*size(prim_int_3cob_nuc_grad(i_grad)%m))

             if( dealloc_ )then
                deallocate( prim_int_3cob_nuc_grad(i_grad)%&
                     &m , stat=status )
                ASSERT(status==0)
             endif
          enddo

          if ( options_xcmode() == xcmode_model_density .or. &
               options_xcmode() == xcmode_extended_mda ) then
             do i_grad=1,gradient_data_n_gradients
                if( alloc_ )then
                   allocate( prim_int_3cob_coul_grad(i_grad)%&
                        &m(n_exp2,n_exp1,n_m2,n_m1) , stat=status )
                   ASSERT(status==0)
                endif

                MEMLOG(inc*size(prim_int_3cob_coul_grad(i_grad)%m))

                if( dealloc_ )then
                   deallocate( prim_int_3cob_coul_grad(i_grad)%&
                        &m , stat=status )
                   ASSERT(status==0)
                endif
             enddo
          endif !options_xcmode

       elseif(allocated(cpks)) then
#ifndef no_cpks_coul_grads
         do i_grad=1,gradient_data_n_gradients
          if( alloc_ )then
           allocate(prim_int_cpks_coul_grad(i_grad)%m &
                           (n_exp2,n_exp1,fit_coeff_n_ch(),n_m2,n_m1) &
                                                  ,stat=cpksalloc(15) )
           ASSERT(cpksalloc(15).eq.0)
          endif

          MEMLOG(inc*size(prim_int_cpks_coul_grad(i_grad)%m))
          cpks_coul_grad_size=cpks_coul_grad_size+inc*size(prim_int_cpks_coul_grad(i_grad)%m)

           if( dealloc_ )then
            deallocate(prim_int_cpks_coul_grad(i_grad)%m,stat=cpksalloc(15))
            ASSERT(cpksalloc(15).eq.0)
            cpksalloc(15)=1
           endif
         enddo
        if(cpks_coul_grad_size.gt.max_cpks_coul_grad_size) then
         max_cpks_coul_grad_size=cpks_coul_grad_size
         DPRINT 'max_cpks_coul_grad_size ::',max_cpks_coul_grad_size
        endif
#endif

       endif splitg
    endif ! integralpar_3cob_grad

    if( integralpar_gradients )then ! not always used, but costs nothing
       ! (de)allocate density matrix section in uncontracted AO basis
       ! needed in property runs:
       if( alloc_ )then
         DPRINT  '2cob3c: allocate(prim_int_dens_mat(',n_exp2,n_exp1,n_m2,n_m1,'))'
         allocate(prim_int_dens_mat(n_exp2,n_exp1,n_m2,n_m1), &
                  stat=cpksalloc(170) )
         ASSERT(cpksalloc(170).eq.0)
       endif
       MEMLOG(inc*size(prim_int_dens_mat))
       if( dealloc_ )then
         DPRINT  '2cob3c: deallocate(prim_int_dens_mat)'
         deallocate(prim_int_dens_mat,stat=cpksalloc(170) )
         ASSERT(cpksalloc(170).eq.0)
         cpksalloc(170)=1
       endif
     endif


    ASSERT(i_p2c==n_2c_pints)
    ASSERT(i_p2cv==n_2cv_pints)
    ASSERT(i_p3c==n_3c_pints)
    ASSERT(i_p3cv==n_3cv_pints)

    if( init_ )then
      call init_storage()
    endif
  contains

    subroutine init_storage()
      implicit none
      ! *** end of interface ***

      ! to intentially break things:
      real(r8_kind), parameter :: INITVAL = HUGE(1.0_r8_kind)
      ! those who do not survive the test:
      real(r8_kind), parameter :: FIXME  = 0.0_r8_kind
      FPP_TIMER_DECL(ini)

      ! this initialization consumes lots
      ! of time on Altix, skip it.
      ! The most ll_calculate subs do it anyway.
      ! If not, please fix the subs!
#define FAILSAFE
#ifndef FAILSAFE
      integer(i4_kind),save :: count = 0
      if(count < 10 )then
      WARN('allocate_primitives: skip initialization')
      elseif( count==10 )then
      WARN('allocate_primitives: SKIP INITIALIZATION')
      WARN('allocate_primitives: AND SKIP WARNINGS ABOUT IT')
      endif
      count = count + 1
      return
#endif

      FPP_TIMER_START(ini)

      ! FIXME: intialization appears to be required as soon
      !        as one disables ss/ls/ll_calculate and lets
      !        SHGI compute (nuclear?) integrals:

      ! Must be fixed, test it with INITVAL=HUGE
      ! in case of doubts
      ! (?pseudo, potential, field?)

      ! initialize the values:
      do i_pi=1,size(prim_2c_ints)
         prim_2c_ints(i_pi)%prim = FIXME
      enddo

      !
      ! Those between goto/continue seem not to require initialization:
      !
      goto 999

      do i_pi=1,size(prim_2cv_ints)
         prim_2cv_ints(i_pi)%prim = INITVAL
      enddo
      do i_pi=1,size(prim_3c_ints)
         prim_3c_ints(i_pi)%prim = INITVAL
      enddo
      do i_pi=1,size(prim_3cv_ints)
         prim_3cv_ints(i_pi)%prim = INITVAL
      enddo

      if ( integralpar_2cob_nuc ) then
         ! or do the "old" way ...
         if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
            prim_int_2cob_nuc_pseudo = INITVAL
         endif
      endif

      ! ... and do the rest traditionally:
      if ( integralpar_2cob_potential ) then
         prim_int_2cob_poten = INITVAL
      endif

999   continue ! Gradients depend on input:
      ! FIX THEM!

      if ( integralpar_2cob_field )then
            prim_int_2cob_field = 0.0_r8_kind
      endif

      if ( integralpar_2cob_ol_grad ) then
         do i_grad=1,size(prim_int_2cob_ol_grad)
            prim_int_2cob_ol_grad(i_grad)%m = FIXME
         enddo

         if(integralpar_dervs) then
            do i_grad=1,size(prim_int_2cob_ol_dervs,1)
               do k2dr=1,size(prim_int_2cob_ol_dervs,2)
                  prim_int_2cob_ol_dervs(i_grad,k2dr)%m = 0.0_r8_kind
               enddo
            enddo
         endif
      endif

      if ( integralpar_2cob_kin_grad ) then
         do i_grad=1,size(prim_int_2cob_kin_grad)
            prim_int_2cob_kin_grad(i_grad)%m = FIXME
         enddo
      endif

      if ( integralpar_2cob_nuc_grad ) then
         do i_grad=1,size(prim_int_2cob_nuc_grad)
            prim_int_2cob_nuc_grad(i_grad)%m = FIXME
         enddo

         if((pseudopot_present.and.integralpar_relativistic) .or. integralpar_pseudo) then
            do i_grad=1,size(prim_int_2cob_pseudo_grad)
               prim_int_2cob_pseudo_grad(i_grad)%m = FIXME
            enddo
         endif !pseudopot_present
      endif ! integralpar_2cob_nuc_grad

      if ( integralpar_2cob_pc_grad ) then
         do i_grad=1,size(prim_int_2cob_pc_grad)
            prim_int_2cob_pc_grad(i_grad)%m = FIXME
         enddo
      endif

      if ( integralpar_2cob_X_grad ) then
         do i_grad=1,size(prim_int_2cob_pd_grad)
            prim_int_2cob_pd_grad(i_grad)%m = FIXME
            prim_int_2cob_pd_torq(i_grad)%m = FIXME
         enddo

         do i_grad=1,size(prim_int_2cob_pq_grad)
            prim_int_2cob_pq_grad(i_grad)%m = FIXME
            prim_int_2cob_pq_torq(i_grad)%m = FIXME
         enddo

         do i_grad=1,size(prim_int_2cob_po_grad)
            prim_int_2cob_po_grad(i_grad)%m = FIXME
            prim_int_2cob_po_torq(i_grad)%m = FIXME
         enddo

         do i_grad=1,size(prim_int_2cob_rc_grad)
            prim_int_2cob_rc_grad(i_grad)%m = FIXME
         enddo
      endif

      if ( integralpar_2cob_ipd_grad ) then
         do i_grad=1,size(prim_int_2cob_ipd_grad)
            prim_int_2cob_ipd_grad(i_grad)%m = FIXME
            prim_int_2cob_ipd_torq(i_grad)%m = FIXME
         enddo
      endif

      if ( integralpar_2cob_pvsp_grad ) then
         do i_grad=1,size(prim_int_2cob_pvsp_grad)
            prim_int_2cob_pvsp_grad(i_grad)%m = FIXME
         enddo
      endif

      if ( integralpar_solv_grad .and. integralpar_cpksdervs ) then
         do i_grad=1,gradient_data_n_spin_gradients
            prim_int_3cob_solv_grad(i_grad)%m = FIXME
         enddo
      endif

      if ( integralpar_3cob_grad ) then
#ifdef WITH_EPE
         if(ewpc_n.ne.0.and.calc_cluster_epe_energy) then
            prim_int_3cob_epe%m=FIXME
         endif
#endif

         do i_grad=1,gradient_data_n_spin_gradients
            prim_int_3cob_grad(i_grad)%m = FIXME
         enddo

         if(integralpar_2dervs) then
            do i_grad=1,size(prim_int_coul_dervs,1)
               do k2dr=1,size(prim_int_coul_dervs,2)
                  prim_int_coul_dervs(i_grad,k2dr)%m = 0.0_r8_kind
               enddo
            enddo
         endif

         if(integralpar_2dervs.and.integralpar_relativistic) then
            do i_grad=1,size(prim_int_nucl_dervs,1)
               do k2dr=1,size(prim_int_nucl_dervs,2)
                  prim_int_nucl_dervs(i_grad,k2dr)%m = 0.0_r8_kind
                  prim_int_pvsp_dervs(i_grad,k2dr)%m = 0.0_r8_kind
               enddo
            enddo
            do i_grad=1,size(prim_int_2cob_kin_dervs,1)
               do k2dr=1,size(prim_int_2cob_kin_dervs,2)
                  prim_int_2cob_kin_dervs(i_grad,k2dr)%m = 0.0_r8_kind
               enddo
            enddo
         endif

         splitg: if (options_split_gradients()) then
            do i_grad=1,size(prim_int_2cob_ks_grad)
               prim_int_2cob_ks_grad(i_grad)%m = FIXME
            enddo

            do i_grad=1,gradient_data_n_gradients
               prim_int_3cob_nuc_grad(i_grad)%m = FIXME
            enddo

            if ( options_xcmode() == xcmode_model_density .or. &
                 options_xcmode() == xcmode_extended_mda ) then
               do i_grad=1,gradient_data_n_gradients
                  prim_int_3cob_coul_grad(i_grad)%m = FIXME
               enddo
            endif !options_xcmode
#ifndef no_cpks_coul_grads
         elseif(allocated(cpks)) then
            do i_grad=1,gradient_data_n_gradients
               prim_int_cpks_coul_grad(i_grad)%m=FIXME
            enddo
#endif
         endif splitg
      endif ! integralpar_3cob_grad

      FPP_TIMER_STOP(ini)
      DPRINT  'TIMER: init ints =',FPP_TIMER_VALUE(ini)
    end subroutine init_storage

  end subroutine allocate_primitives
  !**************************************************************

end module int_data_2cob3c_module
